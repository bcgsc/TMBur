#!/usr/bin/env nextflow
nextflow.preview.dsl=2
include {split_reference; mutect2; merge_vcf_and_stats_files; mark_pass_vcfs} from './mutect2_processes.nf'
include {merge_vcf_files; create_pass_vcfs; split_vcfs_into_SNVs_and_INDELs; gatk_index_feature} from './mutect2_processes.nf'

/*
* pipeline input parameters
*/

//defaults
params.bin_size = 50000000
params.out_dir = "output"
	
log.info """\
Call Somatic Variants Using Mutect2
===================================
samples_file			: ${params.samples_file}
output_folder			: ${params.output_folder}
reference_file			: ${params.reference_file}
reference_index			: ${params.reference_index}
reference_dict			: ${params.reference_dict}
panel_of_normals		: ${params.PON}
panel_of_normals_index	: ${params.PON_index}
germline_resource		: ${params.germline_resource}
germline_resource_index	: ${params.germline_resource_index}
bin_size				: ${params.bin_size}
"""
.stripIndent()


//main workflow
workflow {

	samples = Channel
        .fromPath(params.samples_file)
        .splitCsv(header:true)
		.map{ row-> tuple(
			row.ID, 
			row.tumour_bam,
			row.normal_bam,
			row.tumour_index,
			row.normal_index,
			row.tumour_sample,
			row.normal_sample ) }
    //samples.view()

	//ostensibly constants
	reference_infos = Channel.of([params.reference_file,
        params.reference_dict,
		params.reference_index,
        params.PON,
		params.PON_index,
        params.germline_resource,
		params.germline_resource_index])

	//set up the coordinate ranges
	split_bed = split_reference(params.reference_file, params.reference_dict, params.reference_index, params.bin_size)
	bed_entries = split_bed
		.splitCsv(sep:"\t", header:['chr', 'start', 'end', 'NA1', 'NA2', 'NA3'])
		.map{row  -> tuple(row.chr, row.start, row.end)}

	//make a big matrix combining all of the needed files and a bed range to pass into each Mutect2 job
	//put all the information into a tuple to cross with the bed entries
	for_mutect = samples.combine(reference_infos).combine(bed_entries)

	//run the variant calling jobs
	split_vcfs = mutect2(for_mutect)

	//collect VCFs from the same samples
	collected_vcfs = split_vcfs.groupTuple()
	//collected_vcfs.view()

	//merge and filter the results
	merged_vcfs = merge_vcf_files(collected_vcfs)
	marked_vcfs = mark_pass_vcfs(merged_vcfs,params.reference_file,params.reference_index,params.reference_dict)
	pass_vcfs = create_pass_vcfs(marked_vcfs)
	index_vcfs = gatk_index_feature(pass_vcfs)
	split_vcfs = split_vcfs_into_SNVs_and_INDELs(index_vcfs)

}
