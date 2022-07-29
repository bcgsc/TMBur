include {
	splitReference;
	mutect2;
	mergeVcfAndStatsFiles;
	markPassVcfs
} from './mutect2_processes.nf'

include {
	mergeVcfFiles;
	createPassVcfs;
	splitVcfsIntoSnvsAndIndels;
	gatkIndexFeature
} from './mutect2_processes.nf'

// main workflow
workflow mutect2Wf {
	take:
		samples
        reference_ch
        dict_index
		fai_index

	main:
		// get the coords for the gatk split work
    	createSplitCoordsWf(reference_ch, dict_index, fai_index, params.gatk_bin_size)

		/* Make a big matrix combining all of the needed files and a bed range to pass into each
		   Mutect2 job. Put all the information into a tuple to cross with the bed entries. */
		for_mutect = samples.combine(createSplitCoordsWf.out)

		// run the variant calling jobs
		split_vcfs = mutect2(for_mutect, reference_ch, dict_index, fai_index)

		// collect VCFs from the same samples
		collected_vcfs = split_vcfs.groupTuple(by: [0, 1, 2])

		// merge and filter the results
		merged_vcfs = mergeVcfFiles(collected_vcfs)
		marked_vcfs = markPassVcfs(merged_vcfs, reference_ch, fai_index, dict_index)
		pass_vcfs = createPassVcfs(marked_vcfs)
		index_vcfs = gatkIndexFeature(pass_vcfs)
		split_vcfs = splitVcfsIntoSnvsAndIndels(index_vcfs)

	emit:
		split_vcfs
}

// creates a channel where each entry is the bed coords of a region to use for splitting
workflow createSplitCoordsWf {
	take:
		reference_ch
        dict_index
		fai_index
        gatk_bin_size
	main:
		// set up the coordinate ranges
		split_bed = splitReference(reference_ch, dict_index, fai_index, gatk_bin_size)
		bed_entries = split_bed
			.splitCsv(sep: "\t", header: ['chr', 'start', 'end', 'NA1', 'NA2', 'NA3'])
			.map{ row  -> tuple(row.chr, row.start, row.end) }
	emit:
		bed_entries
}
