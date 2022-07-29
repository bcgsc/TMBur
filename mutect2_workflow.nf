#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include {
	split_reference;
	mutect2;
	merge_vcf_and_stats_files;
	mark_pass_vcfs
} from './mutect2_processes.nf'

include {
	merge_vcf_files;
	create_pass_vcfs;
	split_vcfs_into_SNVs_and_INDELs;
	gatk_index_feature
} from './mutect2_processes.nf'

//main workflow
workflow mutect2_wf {
	take:
		samples
        reference_ch
        dict_index
		fai_index

	main:
		//get the coords for the gatk split work
    	create_split_coords_wf(reference_ch,
			dict_index,
			fai_index,
			params.gatk_bin_size)

		//make a big matrix combining all of the needed files and a bed range to pass into each Mutect2 job
		//put all the information into a tuple to cross with the bed entries
		for_mutect = samples.combine(create_split_coords_wf.out)

		//run the variant calling jobs
		split_vcfs = mutect2(for_mutect, reference_ch, dict_index, fai_index)

		//collect VCFs from the same samples
		collected_vcfs = split_vcfs.groupTuple(by: [0, 1, 2])
		//collected_vcfs.view()

		//merge and filter the results
		merged_vcfs = merge_vcf_files(collected_vcfs)
		marked_vcfs = mark_pass_vcfs(merged_vcfs, reference_ch, fai_index, dict_index)
		pass_vcfs = create_pass_vcfs(marked_vcfs)
		index_vcfs = gatk_index_feature(pass_vcfs)
		split_vcfs = split_vcfs_into_SNVs_and_INDELs(index_vcfs)

	emit:
		split_vcfs
}

//creates a channel where each entry is the bed coords of a region to use for splitting
workflow create_split_coords_wf {
	take:
		reference_ch
        dict_index
		fai_index
        gatk_bin_size
	main:
		//set up the coordinate ranges
		split_bed = split_reference(reference_ch, dict_index, fai_index, gatk_bin_size)
		bed_entries = split_bed
			.splitCsv(sep:"\t", header:['chr', 'start', 'end', 'NA1', 'NA2', 'NA3'])
			.map{row  -> tuple(row.chr, row.start, row.end)}
	emit:
		bed_entries
}
