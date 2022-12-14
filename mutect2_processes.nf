/*
 * This is the main job that runs Mutect2 calling. It is set up to run on a specific region defined
 * by chromosome, start, and end.
 */
process mutect2 {
    tag "${patient}_${T}_${N}_${chromosome}_${start}_${end}"

    input:
        tuple val(patient),
            val(T),
            path(T_bam),
            path(T_bai),
            val(N),
            path(N_bam),
            path(N_bai),
            val(chromosome),
            val(start),
            val(end)
        path reference_file
        path reference_dict
        path reference_index

    output:
        tuple val(patient),
            val(T),
            val(N),
            path("mutect2_${patient}_${T}_${N}_${chromosome}_${start}_${end}.vcf")

    /* Unlike the coordinate system used by other standards such as GFF, the system used by
       the BED format is zero-based for the coordinate start and one-based for the coordinate
       end. Thus, the nucleotide with the coordinate 1 in a genome will have a value of 0 in
       column 2 and a value of 1 in column 3. */
    script:
        """
        start_plus_1=\$((${start}+1))

        /usr/TMB/gatk-4.0.10.0/gatk Mutect2 \
            -R ${reference_file} \
            -I ${T_bam} \
            -tumor ${patient}_${T} \
            -I ${N_bam} \
            -normal ${patient}_${N} \
            -O mutect2_${patient}_${T}_${N}_${chromosome}_${start}_${end}.vcf \
            --intervals ${chromosome}:\${start_plus_1}-${end}
        """
}

// use the fasta file and split it into X sized bins
process splitReference {
    input:
        path reference_fa
        path reference_dict
        path reference_index
        val bin_size

    output:
        path "${reference_fa.baseName}_${bin_size}.bed"

    script:
        """
        /usr/TMB/gatk-4.0.10.0/gatk PreprocessIntervals \
            -R ${reference_fa} \
            -O ${reference_fa.baseName}_${bin_size}.ilist \
            --bin-length ${bin_size} \
            --padding 0

        /usr/TMB/gatk-4.0.10.0/gatk IntervalListToBed \
            -I ${reference_fa.baseName}_${bin_size}.ilist \
            -O ${reference_fa.baseName}_${bin_size}.bed

        grep -v hs37d5 ${reference_fa.baseName}_${bin_size}.bed > tmp.tmp
        mv tmp.tmp ${reference_fa.baseName}_${bin_size}.bed
        """
}

/*
 * Just takes a list of vcfs (in a file) and makes another VCF containing a merge from the list
 * i.e. the "gather" step of a scatter-gather workflow.  Does the same for the stats files which
 * are needed for GATK 4.1+
 */
process mergeVcfAndStatsFiles {
    tag "${patient}_${T}_${N}"

    input:
        tuple val(id), path(unmerged_vcfs_files), path(unmerged_stats_files)

    output:
        tuple val(id),
            path("${patient}_${T}_${N}_merged_unfiltered.vcf"),
            path("${patient}_${T}_${N}_merged_unfiltered.vcf.stats")

    script:
        vcf_file_list = unmerged_vcfs_files.collect{ "-I ${it} " }.join(' ')
        stats_file_list = unmerged_stats_files.collect{ "-stats ${it} " }.join(' ')
        """
        /usr/TMB/gatk-4.0.10.0/gatk MergeVcfs\
            ${vcf_file_list} \
            -O ${patient}_${T}_${N}_merged_unfiltered.vcf

        /usr/TMB/gatk-4.0.10.0/gatk MergeMutectStats\
            ${stats_file_list} \
            -O ${patient}_${T}_${N}_merged_unfiltered.vcf.stats
        """
}

/*
 * Just takes a list of vcfs (in a file) and makes another VCF containing a merge from the list
 * i.e. the "gather" step of a scatter-gather workflow.
 */
process mergeVcfFiles {
    tag "${params.out_dir}/${patient}_${T}_${N}/mutect_raw/"

    input:
        tuple val(patient), val(T), val(N), path(unmerged_vcfs_files)

    output:
        tuple val(patient), val(T), val(N), path("${patient}_${T}_${N}_merged_unfiltered.vcf")

    script:
        vcf_file_list = unmerged_vcfs_files.collect{ "-I ${it} " }.join(' ')
        """
        /usr/TMB/gatk-4.0.10.0/gatk MergeVcfs\
            ${vcf_file_list} \
            -O ${patient}_${T}_${N}_merged_unfiltered.vcf
        """
}

// filter the final file - this just marks variants as PASS or not.
process markPassVcfs {
    tag "${patient}_${T}_${N}"
    publishDir "${params.out_dir}/${patient}_${T}_${N}", mode: 'copy'

    input:
        tuple val(patient),
            val(T),
            val(N),
            path(merged_vcf_file)
        path reference_fa
        path reference_index
        path reference_dict  // needs to be included to be linked in the folder

    output:
        tuple val(patient), val(T), val(N), path("${merged_vcf_file.simpleName}.passmarked.vcf")

    script:
        """
        /usr/TMB/gatk-4.0.10.0/gatk FilterMutectCalls \
            -V ${merged_vcf_file} \
            -O ${merged_vcf_file.simpleName}.passmarked.vcf \
            -R ${reference_fa}
        """
}

// create a file with just the passed variants.
process createPassVcfs {
    tag "${patient}_${T}_${N}"
    publishDir "${params.out_dir}/${patient}_${T}_${N}/mutect_passed/", mode: 'copy'

    input:
        tuple val(patient), val(T), val(N), path(vcf_file)

    output:
        tuple val(patient), val(T), val(N), path("${patient}_${T}_${N}.PASS.vcf")

    script:
        """
        java -Xmx${task.memory.toGiga()}G -jar /usr/TMB/snpEff/SnpSift.jar filter \
            "(FILTER = 'PASS')" \
            ${vcf_file} \
            > ${patient}_${T}_${N}.PASS.vcf
        """
}

process gatkIndexFeature {
    tag "${patient}_${T}_${N}"
    publishDir "${params.out_dir}/${patient}_${T}_${N}/mutect_indexed/", mode: 'copy'

    input:
        tuple val(patient), val(T), val(N), path(vcf_file)

    output:
        tuple val(patient), val(T), val(N), path(vcf_file), path("${vcf_file}.idx")

    script:
        """
        /usr/TMB/gatk-4.0.10.0/gatk IndexFeatureFile -F ${vcf_file}
        """
}

// Split the vcf file into SNVs and INDELs. May omit complex variants that combine SNVs and INDELs.
process splitVcfsIntoSnvsAndIndels {
    tag "${patient}_${T}_${N}"
    publishDir "${params.out_dir}/${patient}_${T}_${N}/mutect_snv_indel_split/", mode: 'copy'

    input:
        tuple val(patient), val(T), val(N), path(vcf_file), path(vcf_index)

    output:
        tuple val('Mutect2'),
            val(patient),
            val(T),
            val(N),
            path("mutect_${patient}_${T}_${N}.PASS.indel.vcf.gz"),
            path("mutect_${patient}_${T}_${N}.PASS.snv.vcf.gz"),
            path("mutect_${patient}_${T}_${N}.PASS.indel.vcf.gz.tbi"),
            path("mutect_${patient}_${T}_${N}.PASS.snv.vcf.gz.tbi")

    script:
        """
        /usr/TMB/gatk-4.0.10.0/gatk SplitVcfs \
            -I ${vcf_file} \
            --INDEL_OUTPUT mutect_${patient}_${T}_${N}.PASS.indel.vcf.gz \
            --SNP_OUTPUT mutect_${patient}_${T}_${N}.PASS.snv.vcf.gz \
            --STRICT false
        """
}
