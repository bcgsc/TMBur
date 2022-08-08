#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    copyReference;
    bwaIndexReference;
    samtoolsIndexReference;
    gatkIndexReference;
    createRtgReference
} from './fastq_to_TMB_processes.nf'

include {
    countFastaBases;
    countCdsBases
} from './fastq_to_TMB_processes.nf'

include {
    trimPair;
    alignReads;
    sortBam;
    mergeBams;
    markDuplicates
} from './fastq_to_TMB_processes.nf'

include {
    msiSensor2;
    manta;
    strelka;
    createPassVcfsStrelka;
    rtgIntersectCalls
} from './fastq_to_TMB_processes.nf'

include {
    annotateSmallVariants;
    createReport;
    createSignatures;
    createPanelReport
} from './fastq_to_TMB_processes.nf'

include {
    mutect2Wf;
    createSplitCoordsWf
} from './mutect2_workflow.nf'

// pipeline input parameters
params.samples_file = 'samples.csv'
params.out_dir = 'output'

log.info """\
TMB estimation pipeline ${params.release}
===================================
samples file  : ${params.samples_file}
output dir    : ${params.out_dir}
reference     : ${params.reference_file}
annotation    : ${params.annotation}
GATK bin size : ${params.gatk_bin_size}
COSMIC VCF    : ${params.cosmic_vcf}
"""

// main workflow
workflow {
    // load in the samples file
    samples = Channel
        .fromPath(params.samples_file)
        .splitCsv(header: true)
        .map{ row-> tuple(file(row.read1).baseName, row.patient, row.tissue, file(row.read1), file(row.read2)) }

    // set up the reference files
    reference_ch = copyReference(params.reference_file)
    bwa_index = bwaIndexReference(reference_ch)
    fai_index = samtoolsIndexReference(reference_ch)
    dict_index = gatkIndexReference(reference_ch)
    base_count_file = countFastaBases(reference_ch)
    countCdsBases(params.annotation)
    rtg_reference = createRtgReference(reference_ch)

    // preprocess fastq files
    trimPair(samples)

    // align files
    alignReads(trimPair.out.fastqs, reference_ch, bwa_index)
    sortBam(alignReads.out)

    // split into 2 channels, one for merging bams, the other uses bams from one fastq pair
    bam_lists = sortBam.out.groupTuple(by: [0, 1]).branch {
        single: it[3].size() <= 1
        multiple: it[3].size() > 1
    }
    // merge where multiple bams from same patient-tissue combination
    mergeBams(bam_lists.multiple)

    // reformat non-merged entries to match output of merge_bam
    single_bams = bam_lists.single.map { patient, tissue, an_id, bam -> [patient, tissue, bam[0]] }

    // markDuplicates on all single and merged bams
    markDuplicates(single_bams.mix(mergeBams.out))

    // to set up for somatic analysis get the cross product of the tumours and normals per patient
    split_tissues = markDuplicates.out.branch{
	    tumour : it[1] =~ /^T.*/
	    normal : it[1] =~ /^N.*/
	}

    bams_for_somatic_ch = split_tissues.tumour.combine(split_tissues.normal, by: 0)

    // call all the somatic analyses
    msiSensor2(bams_for_somatic_ch)
    manta(bams_for_somatic_ch, fai_index, reference_ch) // only works for deep genomes
    strelka(manta.out, fai_index, reference_ch)
    createPassVcfsStrelka(strelka.out)
    mutect2Wf(bams_for_somatic_ch, reference_ch, dict_index, fai_index)

    // merge Strelka and GATK SNVs and INDELs
    rtgIntersectCalls(createPassVcfsStrelka.out.join(mutect2Wf.out, by: [1, 2, 3]), rtg_reference)

    // annotate
    annotateSmallVariants(rtgIntersectCalls.out.joined_calls)

    // merge all the tools into one data structure
    all_results = msiSensor2.out.join(annotateSmallVariants.out.annnotations, by: [1, 2, 3])

    // produce AF report
    createReport(base_count_file, countCdsBases.out.cds_size_file, countCdsBases.out.cds_bed, all_results, params.release)

    // produce signature results
    createSignatures(rtgIntersectCalls.out.joined_calls)

    // panel_foundation_one (laura)
    createPanelReport(base_count_file, countCdsBases.out.cds_size_file, countCdsBases.out.cds_bed, params.cosmic_vcf, all_results)
}

// will run at the end of the analysis.
workflow.onComplete {
	println "Pipeline $workflow.scriptName completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
