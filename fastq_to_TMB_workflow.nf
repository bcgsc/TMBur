#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {copy_reference; bwa_index_reference; samtools_index_reference; gatk_index_reference; create_RTG_reference} from './fastq_to_TMB_processes.nf'
include {count_fasta_bases; count_CDS_bases} from './fastq_to_TMB_processes.nf'
include {trim_pair; align_reads; sort_bam; merge_bams; mark_duplicates} from './fastq_to_TMB_processes.nf'
include {MSIsensor2; manta; strelka; create_pass_vcfs_strelka; rtg_intersect_calls} from './fastq_to_TMB_processes.nf'
include {annotate_small_variants; create_report; create_signatures; create_panel_report} from './fastq_to_TMB_processes.nf'
include {mutect2_wf; create_split_coords_wf} from './mutect2_workflow.nf'


/*
* pipeline input parameters
*/
	params.samples_file = "/projects/rcorbettprj2/mutationalBurden/PROFYLE_container/2p0/test_data/samples.csv"
	params.out_dir = "/projects/rcorbettprj2/mutationalBurden/PROFYLE_container/2p0/test_data/output"
    params.release = "0_5_0"
		
	log.info """\
    TMB estimation pipeline
    ===================================
    samples_file    : ${params.samples_file}
    out_dir         : ${params.out_dir}
    reference       : ${params.reference_file}
    annotation      : ${params.annotation}
    """
	.stripIndent()


//main workflow
workflow {

    //Load in the samples file
    samples = Channel
        .fromPath(params.samples_file)
        .splitCsv(header:true)
        .map{ row-> tuple(file(row.read1).baseName, row.patient, row.tissue, file(row.read1), file(row.read2)) }
    
    //print out the CSV data
    //samples.view()
    
    //Set up the reference files 
    reference_ch = copy_reference(params.reference_file)
    bwa_index = bwa_index_reference(reference_ch)
    fai_index = samtools_index_reference(reference_ch)
    dict_index = gatk_index_reference(reference_ch)
    base_count_file = count_fasta_bases(reference_ch)
    count_CDS_bases(params.annotation)
    RTG_reference = create_RTG_reference(reference_ch)

    //preprocess fastq files
    trim_pair(samples)
    //trim_pair.out.fastqs.view()

    //align files
    align_reads(trim_pair.out.fastqs,reference_ch,bwa_index)
    sort_bam(align_reads.out)

    //merge + dupmark
    // split into 2 channels, one for merging bams, the other uses bams from one
    // fastq pair
    bam_lists = sort_bam.out.groupTuple(by:[0, 1]).branch {
        single: it[3].size() <= 1
        multiple: it[3].size() > 1
    } 
    //merge where multiple bams from same patient-tissue combination
    merge_bams(bam_lists.multiple)
    //reformat non-merged entries to match output of merge_bam
    single_bams = bam_lists.single.map {
        patient, tissue, an_id, bam -> [patient, tissue, bam[0]]
    }

    //mark_duplicates on all single and merged bams 
    mark_duplicates(single_bams.mix(merge_bams.out))

    //to set up for somatic analysis get the cross product of the tumours and normals per patient
    split_tissues = mark_duplicates.out.branch{
	    tumour : it[1] =~ /^T.*/
	    normal : it[1] =~ /^N.*/
	}
    //split_tissues.tumour.view { "[T/N] $it is a tumour" }
    //split_tissues.normal.view { "[T/N] $it is a normal" }

    bams_for_somatic_ch = split_tissues.tumour.combine(split_tissues.normal, by : 0)
	//bams_for_somatic_ch.view { "[crossed] $it" }

    //call all the somatic analyses
    MSIsensor2(bams_for_somatic_ch)
    manta(bams_for_somatic_ch, fai_index, reference_ch) // only works for deep genomes
    strelka(manta.out, fai_index, reference_ch)
    create_pass_vcfs_strelka(strelka.out)
    mutect2_wf(bams_for_somatic_ch,
        reference_ch,
        dict_index,
		fai_index)
        
    //merge Strelka and GATK SNVs and INDELs
    rtg_intersect_calls(create_pass_vcfs_strelka.out.join(mutect2_wf.out, by: [1,2,3]), RTG_reference)

    //annotate
    annotate_small_variants(rtg_intersect_calls.out.joined_calls)

    //merge all the tools into one data structure
    all_results = MSIsensor2.out.join(annotate_small_variants.out.annnotations, by:[1,2,3])
    
    //produce AF report
    create_report(base_count_file, count_CDS_bases.out.CDS_size_file, count_CDS_bases.out.CDS_bed, all_results, params.release)

    //produce signature results
    create_signatures(rtg_intersect_calls.out.joined_calls)

    //panel_foundation_one (laura)
    create_panel_report(base_count_file, count_CDS_bases.out.CDS_size_file, count_CDS_bases.out.CDS_bed, params.cosmic_vcf, all_results)
    
    //QC?
     

}

//will run at the end of the analysis.
workflow.onComplete {
	//final_results_ch.view{ println "[complete] $it"}
	println "Pipeline $workflow.scriptName completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}