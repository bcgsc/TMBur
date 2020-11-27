#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
binaries
*/
gatk="/gsc/software/linux-x86_64-centos6/gatk-4.0.10.0/gatk" //production
//gatk="/home/rcorbett/bin/gatk-4.1.9.0/gatk" //latest
java_path="/gsc/software/linux-x86_64-centos6/jdk1.8.0_162/bin"
sambamba="/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5"
snpSift="/gsc/software/linux-x86_64-centos7/snpEff-4.3/SnpSift.jar"


//This is the main job that runs Mutect2 calling.
//It is set up to run on a specific region defined by 
//chromosome, start, and end.
process mutect2 {
    tag "${id}_${chromosome}_${start}_${end}"
    cpus 4 //seems to use about 400% CPU per job

    input:
        tuple val(id), 
            path(tumour_bam),
            path(normal_bam),
            path(tumour_bam_index),
            path(normal_bam_index),
            val(tumour_sample_name),
            val(normal_sample_name),
            path(reference_file),
            path(reference_dict),
            path(reference_index),
            path(panel_of_normals),
            path(panel_of_normals_index),
            path(germline_resource),
            path(germline_resource_index),
            val(chromosome),
            val(start), 
            val(end)

    output:
        tuple val(id), 
            path("mutect2_${id}_${chromosome}_${start}_${end}.vcf")
            
    /* Unlike the coordinate system used by other standards such as GFF, the system used by 
        the BED format is zero-based for the coordinate start and one-based for the coordinate 
        end. Thus, the nucleotide with the coordinate 1 in a genome will have a value of 0 in 
        column 2 and a value of 1 in column 3. */
    script:
    """
        start_plus_1=\$((${start}+1))

        export PATH=${java_path}:\$PATH;
        ${gatk} Mutect2 \
        -R ${reference_file} \
        -I ${tumour_bam} \
        -tumor ${tumour_sample_name} \
        -I ${normal_bam} \
        -normal ${normal_sample_name} \
        -O mutect2_${id}_${chromosome}_${start}_${end}.vcf \
        --germline-resource  ${germline_resource}  \
        --af-of-alleles-not-in-resource 0.0000322 \
        --panel-of-normals ${panel_of_normals} \
        --intervals ${chromosome}:\${start_plus_1}-${end}
    """
}

//use the fasta file and split it into X sized bins
process split_reference {   
    input:
        path reference_fa
        path reference_dict
        path reference_index
        val bin_size

    output:
        path("${reference_fa.baseName}_${bin_size}.bed")

    script:
    """
        export PATH=${java_path}:\$PATH;
        ${gatk} PreprocessIntervals -R ${reference_fa} -O ${reference_fa.baseName}_${bin_size}.ilist --bin-length ${bin_size} --padding 0
        ${gatk} IntervalListToBed -I ${reference_fa.baseName}_${bin_size}.ilist -O ${reference_fa.baseName}_${bin_size}.bed
    """
}


//just takes a list of vcfs (in a file) and makes another VCF containing a merge from the list
//ie. the "gather" step of a scatter-gather workflow.  Does the same for the stats files
//which are needed for GATK 4.1+
process merge_vcf_and_stats_files {
    tag "${id}"

    input:
        tuple val(id), path(unmerged_vcfs_files), path(unmerged_stats_files)

    output:
        tuple val(id), path("${id}_merged_unfiltered.vcf"), path("${id}_merged_unfiltered.vcf.stats")

    script:
    vcf_file_list = unmerged_vcfs_files.collect{ "-I ${it} " }.join(' ')
    stats_file_list = unmerged_stats_files.collect{ "-stats ${it} " }.join(' ')
    """
        export PATH=${java_path}:\$PATH;
        ${gatk} MergeVcfs\
        ${vcf_file_list} \
        -O ${id}_merged_unfiltered.vcf

        ${gatk} MergeMutectStats\
        ${stats_file_list} \
        -O ${id}_merged_unfiltered.vcf.stats
    """
}

//just takes a list of vcfs (in a file) and makes another VCF containing a merge from the list
//ie. the "gather" step of a scatter-gather workflow. 
process merge_vcf_files {
    tag "${id}"

    input:
        tuple val(id), path(unmerged_vcfs_files)

    output:
        tuple val(id), path("${id}_merged_unfiltered.vcf")

    script:
    vcf_file_list = unmerged_vcfs_files.collect{ "-I ${it} " }.join(' ')
    """
        export PATH=${java_path}:\$PATH;
        ${gatk} MergeVcfs\
        ${vcf_file_list} \
        -O ${id}_merged_unfiltered.vcf

    """
}

//filter the final file - this just marks variants as PASS or not.
process mark_pass_vcfs {
    tag "${id}"
    publishDir "${params.output_folder}/${id}", mode: 'copy'

    input:
        tuple val(id), path(merged_vcf_file)
        path(reference_fa)
        path(reference_index)
        path(reference_dict) // needs to be included to be linked in the folder

    output:
        tuple val(id), path("${merged_vcf_file.simpleName}.passmarked.vcf")

    script:
    """
        export PATH=${java_path}:\$PATH;
        ${gatk} FilterMutectCalls\
        -V ${merged_vcf_file} \
        -O ${merged_vcf_file.simpleName}.passmarked.vcf \
        -R ${reference_fa}
    """  
}

//create a file with just the passed variants.
process create_pass_vcfs {
    tag "${id}"
    publishDir "${params.output_folder}/${id}", mode: 'copy'

    input:
        tuple val(id), path(vcf_file)
        
    output:
        tuple val(id), path("${id}.PASS.vcf")

    script:
    """
        export PATH=${java_path}:\$PATH;
        java -Xmx2g -jar /gsc/software/linux-x86_64-centos7/snpEff-4.3/SnpSift.jar filter "(FILTER = 'PASS')"  ${vcf_file} > ${id}.PASS.vcf
    """  
}

process gatk_index_feature {
    tag "${id}"
    publishDir "${params.output_folder}/${id}", mode: 'copy'

    input:
        tuple val(id), path(vcf_file)
        
    output:
        tuple val(id), path(vcf_file), path("${vcf_file}.idx")

    script:
    """
        export PATH=${java_path}:\$PATH;
        ${gatk} IndexFeatureFile -F ${vcf_file}
    """  
}

//split the vcf file into SNVs and INDELs.
process split_vcfs_into_SNVs_and_INDELs {
    tag "${id}"
    publishDir "${params.output_folder}/${id}", mode: 'copy'

    input:
        tuple val(id), path(vcf_file), path(vcf_index)
        
    output:
        tuple path("${id}.PASS.indel.vcf.gz"), path("${id}.PASS.snv.vcf.gz")

    script:
    """
        export PATH=${java_path}:\$PATH;
        ${gatk} SplitVcfs \
        -I ${vcf_file} \
        --INDEL_OUTPUT ${id}.PASS.indel.vcf.gz \
        --SNP_OUTPUT ${id}.PASS.snv.vcf.gz \
        --STRICT false
    """  
}