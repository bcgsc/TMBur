#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//an ugly way to access the file that is distributed in the container
process copy_reference {
	tag "$fasta_name"
	input:
		val(fasta_name) 

	output:
		path("*.fa") 

	script:
	"""
		echo ${fasta_name} | xargs -i cp {} .
	"""
}

process bwa_index_reference {
	tag "$fasta"
	input:
		path(fasta) 

	output:
		path("${fasta}.*") 
		
	script:
	"""
		/usr/TMB/bwa index ${fasta}
	"""
}

process samtools_index_reference {
	tag "$fasta"
	input:
		path(fasta) 

	output:
		path("${fasta}.fai") 
        
	script:
	"""
		/usr/TMB/samtools faidx ${fasta}
    """
}

process gatk_index_reference {
	tag "$fasta"
	input:
		path(fasta) 

	output:
        path("${fasta.baseName}.dict")

	script:
	"""
        /usr/TMB/gatk-4.0.10.0/gatk CreateSequenceDictionary -R ${fasta}
	"""
}

process create_RTG_reference {
	tag "$fasta"
	input:
		path(fasta) 

	output:
        path("SDF")

	script:
	"""
        /usr/TMB/rtg-tools-3.11/rtg format -o SDF ${fasta}
	"""
}

process count_fasta_bases {
	tag "$fasta"
	input:
		path(fasta) 

	output:
		path("base_count.txt") 

	script:
	"""
		seqtk comp ${fasta} | grep -E '^[123456789XY]{1,2}\\s' | awk '{ print \$3+\$4+\$5+\$6 }' | paste -sd+ | bc > base_count.txt
	"""
}

process count_CDS_bases {
	tag "$anno"
	input:
		val(anno) 

	output:
		path("CDS_size.txt")

	script:
	"""
		java -jar /usr/TMB/snpEff/snpEff.jar  dump -v -bed $anno > ${anno}.bed
		grep -E '^[1234567890XY]{1,2}\\s' ${anno}.bed | grep CDS | bedtools sort | bedtools merge | awk '{ print \$3-\$2 }' | paste -sd+ | bc > CDS_size.txt
	"""
}

process trim_pair {
    tag "$an_id"
	cpus 10

	input:
		tuple val(an_id), val(patient), val(tissue), path(reads1), path(reads2)

	output:
        tuple val(an_id), val(patient), val(tissue), path("${reads1}.fastp.trimmed.fq.gz"),
			path("${reads2}.fastp.trimmed.fq.gz"),  emit: fastqs
        path("${an_id}.fastp.json"),                emit: jsons
		
	script:
 	"""
		fastp \
		-i ${reads1} \
		-I ${reads2} \
		-o ${reads1}.fastp.trimmed.fq.gz \
		-O ${reads2}.fastp.trimmed.fq.gz \
		--json ${an_id}.fastp.json
 	"""
}

process align_reads {
	tag "$an_id"
	cpus 48
	memory "48 GB"

	input:
		tuple val(an_id), val(patient), val(tissue), path(trim1), path(trim2) 
		path(reference_fasta) 
		path(bwa_index) 

	output:
		tuple val(an_id), val(patient), val(tissue), path("${an_id}.bam") 

	script:
	"""
		bwa mem \
		${reference_fasta} \
		-R '@RG\\tID:${an_id}_${patient}_${tissue}\\tSM:${patient}_${tissue}\\tLB:${an_id}_${patient}_${tissue}\\tPL:ILLUMINA' \
		-t 48 ${trim1} ${trim2}  | /usr/TMB/samtools view -S -b - > ${an_id}.bam
	"""
}

process sort_bam {
	tag "$an_id"
	cpus 10

	input:
		tuple val(an_id), val(patient), val(tissue), path(bam_file) 

	output: 
		tuple val(patient), val(tissue), val(an_id), path("${bam_file.baseName}.sorted.bam") 

	script:
	"""
		sambamba sort --tmpdir . ${bam_file}
	"""
}

process merge_bams {
	tag "${patient}_${tissue}"
	memory "64 GB"
	cpus 16

	input:
		tuple val(patient), val(tissue), val(an_id), path(bams_to_merge)

	output:
		tuple val(patient), val(tissue), path("${patient}-${tissue}_merged.bam") 

	script:
	"""
		sambamba merge ${patient}-${tissue}_merged.bam $bams_to_merge
	"""
}

process mark_duplicates {
	tag "${patient}_${tissue}"
	publishDir "${params.out_dir}/bams/${patient}_bams"
	memory "128 GB"
	cpus 8

	input:
		tuple val(patient), val(tissue), path(bam_file)

	output:
		tuple val(patient), val(tissue), path("${bam_file.baseName}.dup.bam"),
		path("${bam_file.baseName}.dup.bam.bai") 

	script:
	"""
		sambamba markdup \
		--nthreads 8 \
		--tmpdir . \
		--hash-table-size 5000000 \
		--overflow-list-size 5000000 \
		--io-buffer-size 1024 ${bam_file} ${bam_file.baseName}.dup.bam
	"""
}

process MSIsensor2 {
	tag "${patient}_${T}_${N}"
	publishDir "${params.out_dir}/${patient}_${T}_${N}/MSIsensor2"
	cpus 32

	input:
		tuple val(patient), 
			val(T), 
			path(T_bam), 
			path(T_bai), 
			val(N), 
			path(N_bam),
			path(N_bai) 
		
	output:
		tuple val("msisensor2"), val(patient), val(T), val(N), 
			path("msisensor2_${patient}_${T}_${N}.txt") 

	script: //Had to take the "2" off the binary name due to folder naming clashes in the container
	"""
		msisensor msi \
		-t ${T_bam} \
		-n ${N_bam} \
		-d /usr/TMB/msisensor2/models_b37_HumanG1Kv37/1030c0aa35ca5c263daeae866ad18632 \
		-b 32 \
		-o msisensor2_${patient}_${T}_${N}.txt 2> msisensor2_out_${patient}_${T}_${N}.log
	"""
}

//Needs to be run on Deep bams.  If testing on shallow bams (single-digit X coverage)
//this may fail due to no candidate events being detected
process manta {
	tag "${patient}_${T}_${N}"
	publishDir "${params.out_dir}/${patient}_${T}_${N}/Manta"
	cpus 48

	input:
		tuple val(patient), 
			val(T), 
			path(T_bam), 
			path(T_bai), 
			val(N), 
			path(N_bam),
			path(N_bai) 
		path(fai_file)
		path(reference) 

	output:
		tuple val("Manta"), 
			val(patient), 
			val(T), 
			path(T_bam), 
			path(T_bai), 
			val(N), 
			path(N_bam),
			path(N_bai), 
			path("Manta_${patient}_${T}_vs_${N}.candidateSmallIndels.vcf.gz"), 
			path("Manta_${patient}_${T}_vs_${N}.candidateSmallIndels.vcf.gz.tbi"),
			path("Manta_${patient}_${T}_vs_${N}.candidateSV.vcf.gz"), 
			path("Manta_${patient}_${T}_vs_${N}.candidateSV.vcf.gz.tbi"),
			path("Manta_${patient}_${T}_vs_${N}.somaticSV.vcf.gz"), 
			path("Manta_${patient}_${T}_vs_${N}.somaticSV.vcf.gz.tbi")

 	script:
 	"""
		/usr/TMB/manta-1.6.0.centos6_x86_64/bin/configManta.py \
			--normalBam=${N_bam} \
			--tumorBam=${T_bam} \
			--referenceFasta=${reference}  \
			--runDir Manta
		python Manta/runWorkflow.py -m local -j 48

		mv Manta/results/variants/candidateSmallIndels.vcf.gz \
			Manta_${patient}_${T}_vs_${N}.candidateSmallIndels.vcf.gz
		mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
			Manta_${patient}_${T}_vs_${N}.candidateSmallIndels.vcf.gz.tbi
		mv Manta/results/variants/candidateSV.vcf.gz \
			Manta_${patient}_${T}_vs_${N}.candidateSV.vcf.gz
		mv Manta/results/variants/candidateSV.vcf.gz.tbi \
			Manta_${patient}_${T}_vs_${N}.candidateSV.vcf.gz.tbi
		mv Manta/results/variants/somaticSV.vcf.gz \
			Manta_${patient}_${T}_vs_${N}.somaticSV.vcf.gz
		mv Manta/results/variants/somaticSV.vcf.gz.tbi \
			Manta_${patient}_${T}_vs_${N}.somaticSV.vcf.gz.tbi
 """
}


process strelka {
	tag "${patient}_${T}_${N}"
	publishDir "${params.out_dir}/${patient}_${T}_${N}/strelka"
	cpus 48

	input:
		tuple val(tool), 
			val(patient), 
			val(T), 
			path(T_bam), 
			path(T_bai), 
			val(N), 
			path(N_bam),
			path(N_bai), 
			path(SmallIndels), 
			path(SmallIndels_index),
			path(candidateSV), 
			path(candidateSV_index),
			path(somaticSV), 
			path(somaticSV_index)
		path(fai_file)
		path(reference) 

	output:
		tuple val("Strelka2"), val(patient), val(T), val(N), 
			path("Strelka_${patient}_${T}_vs_${N}_somatic_indels.vcf.gz"),
			path("Strelka_${patient}_${T}_vs_${N}_somatic_snvs.vcf.gz"),
			path("Strelka_${patient}_${T}_vs_${N}_somatic_indels.vcf.gz.tbi"),
			path("Strelka_${patient}_${T}_vs_${N}_somatic_snvs.vcf.gz.tbi")


	script:
	"""
		/usr/TMB/strelka-2.9.2.centos6_x86_64/bin//configureStrelkaSomaticWorkflow.py \
		    --tumor ${T_bam} \
			--normal ${N_bam} \
			--referenceFasta ${reference} \
			--runDir Strelka \
			--indelCandidates ${SmallIndels}
		python Strelka/runWorkflow.py -m local -j 48

		rm -f *gz *tbi

		mv Strelka/results/variants/somatic.indels.vcf.gz \
			Strelka_${patient}_${T}_vs_${N}_somatic_indels.vcf.gz
		mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
			Strelka_${patient}_${T}_vs_${N}_somatic_indels.vcf.gz.tbi
		mv Strelka/results/variants/somatic.snvs.vcf.gz \
			Strelka_${patient}_${T}_vs_${N}_somatic_snvs.vcf.gz
		mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
			Strelka_${patient}_${T}_vs_${N}_somatic_snvs.vcf.gz.tbi
	"""
}

//create VCF files with just the passed variants.
//Will be zipped/tabix
process create_pass_vcfs_strelka {
    tag "${patient}_${T}_${N}"
    publishDir "${params.out_dir}/${patient}_${T}_${N}/strelka_passed", mode: 'copy'

    input:
        tuple val(caller), val(patient), val(T), val(N), 
			path(strelka_indels),
			path(strelka_snvs),
			path(strelka_indels_index),
			path(strelka_snvs_index)
        
    output:
        tuple val("${caller}"), val(patient), val(T), val(N), 
			path("strelka_${patient}_${T}_${N}.PASS.indel.vcf.gz"),
			path("strelka_${patient}_${T}_${N}.PASS.snv.vcf.gz"),
			path("strelka_${patient}_${T}_${N}.PASS.indel.vcf.gz.tbi"),
			path("strelka_${patient}_${T}_${N}.PASS.snv.vcf.gz.tbi")

    script:
    """
        java -Xmx2g -jar /usr/TMB/snpEff/SnpSift.jar filter "(FILTER = 'PASS')"  ${strelka_indels} > strelka_${patient}_${T}_${N}.PASS.indel.vcf
		java -Xmx2g -jar /usr/TMB/snpEff/SnpSift.jar filter "(FILTER = 'PASS')"  ${strelka_snvs} > strelka_${patient}_${T}_${N}.PASS.snv.vcf 
		/usr/TMB/rtg-tools-3.11/rtg bgzip strelka_${patient}_${T}_${N}.PASS.indel.vcf 
		/usr/TMB/rtg-tools-3.11/rtg bgzip strelka_${patient}_${T}_${N}.PASS.snv.vcf 
        /usr/TMB/rtg-tools-3.11/rtg index strelka_${patient}_${T}_${N}.PASS.indel.vcf.gz
		/usr/TMB/rtg-tools-3.11/rtg index strelka_${patient}_${T}_${N}.PASS.snv.vcf.gz
    """  
}


//Use RTGtools to intersect the variants.   SNVs and INDELs are intersected separately
process rtg_intersect_calls {
   	tag "${patient}_${T}_${N}"
    publishDir "${params.out_dir}/${patient}_${T}_${N}/variant_intersect", mode: 'copy'
	cpus 10
	memory "16 GB"

    input:
		tuple val(patient), val(T), val(N), 
			val(strelka_name),
			path(strelka_indels),
			path(strelka_snvs),
			path(strelka_indels_index),
			path(strelka_snvs_index),
			val(mutect_name),
			path(mutect_indels),
			path(mutect_snvs),
			path(mutect_indels_index),
			path(mutect_snvs_index)
		path(rtg_reference)
    
    output:
        tuple val("intersected_variants"), val(patient), val(T), val(N), 
			path("${patient}_${T}_${N}_snv_strelka_only.vcf.gz"),
			path("${patient}_${T}_${N}_snv_mutect_only.vcf.gz"),
			path("${patient}_${T}_${N}_snv_both.vcf.gz"),
			path("${patient}_${T}_${N}_snv_strelka_only.vcf.gz.tbi"),
			path("${patient}_${T}_${N}_snv_mutect_only.vcf.gz.tbi"),
			path("${patient}_${T}_${N}_snv_both.vcf.gz.tbi"),
			path("${patient}_${T}_${N}_indel_strelka_only.vcf.gz"),
			path("${patient}_${T}_${N}_indel_mutect_only.vcf.gz"),
			path("${patient}_${T}_${N}_indel_both.vcf.gz"),
			path("${patient}_${T}_${N}_indel_strelka_only.vcf.gz.tbi"),
			path("${patient}_${T}_${N}_indel_mutect_only.vcf.gz.tbi"),
			path("${patient}_${T}_${N}_indel_both.vcf.gz.tbi"), emit: all_calls
		tuple val("intersected_variants"), val(patient), val(T), val(N), 
			path("${patient}_${T}_${N}_snv_both.vcf.gz"),
			path("${patient}_${T}_${N}_snv_both.vcf.gz.tbi"),
			path("${patient}_${T}_${N}_indel_both.vcf.gz"),
			path("${patient}_${T}_${N}_indel_both.vcf.gz.tbi"), emit: joined_calls

			
    script:
    """
		/usr/TMB/rtg-tools-3.11/rtg RTG_MEM=16G vcfeval \
		-b ${strelka_snvs} \
		-c ${mutect_snvs} \
		-t ${rtg_reference} \
		-o SNV_intersect \
		--sample ALT,ALT \
		--squash-ploidy \
		--vcf-score-field INFO.QSS 

		mv SNV_intersect/fp.vcf.gz ${patient}_${T}_${N}_snv_mutect_only.vcf.gz
		mv SNV_intersect/fn.vcf.gz ${patient}_${T}_${N}_snv_strelka_only.vcf.gz
		mv SNV_intersect/tp.vcf.gz ${patient}_${T}_${N}_snv_both.vcf.gz
		mv SNV_intersect/fp.vcf.gz.tbi ${patient}_${T}_${N}_snv_mutect_only.vcf.gz.tbi
		mv SNV_intersect/fn.vcf.gz.tbi ${patient}_${T}_${N}_snv_strelka_only.vcf.gz.tbi
		mv SNV_intersect/tp.vcf.gz.tbi ${patient}_${T}_${N}_snv_both.vcf.gz.tbi

		/usr/TMB/rtg-tools-3.11/rtg RTG_MEM=16G vcfeval \
		-b ${strelka_indels} \
		-c ${mutect_indels} \
		-t ${rtg_reference} \
		-o INDEL_intersect \
		--squash-ploidy \
		--sample ALT,ALT \
		--vcf-score-field INFO.QSS

		mv INDEL_intersect/fp.vcf.gz ${patient}_${T}_${N}_indel_mutect_only.vcf.gz
		mv INDEL_intersect/fn.vcf.gz ${patient}_${T}_${N}_indel_strelka_only.vcf.gz
		mv INDEL_intersect/tp.vcf.gz ${patient}_${T}_${N}_indel_both.vcf.gz
		mv INDEL_intersect/fp.vcf.gz.tbi ${patient}_${T}_${N}_indel_mutect_only.vcf.gz.tbi
		mv INDEL_intersect/fn.vcf.gz.tbi ${patient}_${T}_${N}_indel_strelka_only.vcf.gz.tbi
		mv INDEL_intersect/tp.vcf.gz.tbi ${patient}_${T}_${N}_indel_both.vcf.gz.tbi
    """  

}


//annotate the final calls with SNPEff
process annotate_small_variants {
	tag "${patient}_${T}_${N}"
	publishDir "${params.out_dir}/${patient}_${T}_${N}/annotated_variants", mode: 'copy'
	memory '48 GB'

	input:
		tuple val(names), val(patient), val(T), val(N), 
			path(snv_calls),
			path(snv_index),
			path(indel_calls),
			path(indel_index)

	output:
		tuple val("snpEff"), val(patient), val(T), val(N), 
			path("${snv_calls.baseName}.snpEff.vcf"),
			path("${indel_calls.baseName}.snpEff.vcf"), emit: annnotations
		tuple val("snpEff"), val(patient), val(T), val(N), 
			path("${patient}_${T}_${N}_somatic.snv.html"),
			path("${patient}_${T}_${N}_somatic.indel.html"), emit: stats_files

	script:
	"""
		java -Xmx48g -jar /usr/TMB/snpEff/snpEff.jar GRCh37.75 -s ${patient}_${T}_${N}_somatic.snv.html ${snv_calls} > ${snv_calls.baseName}.snpEff.vcf
		java -Xmx48g -jar /usr/TMB/snpEff/snpEff.jar GRCh37.75 -s ${patient}_${T}_${N}_somatic.indel.html ${indel_calls} > ${indel_calls.baseName}.snpEff.vcf
	"""
}

process create_signatures {
	tag "${patient}_${T}_${N}"
	publishDir "${params.out_dir}/${patient}_${T}_${N}/signatures", mode: 'copy'
	
	input:
		tuple val(names), val(patient), val(T), val(N), 
			path(snv_calls),
			path(snv_index),
			path(indel_calls),
			path(indel_index)

	output:
		tuple val("SigProfiler"), val(patient), val(T), val(N), 
			path("mutation_signature_output"), emit: mutsig_output
		
	script:
	"""
		/usr/TMB/gatk-4.0.10.0/gatk MergeVcfs -I ${snv_calls} -I ${indel_calls} -O all_passed_merged_variants.vcf
		mkdir vcf_input
		mv all_passed_merged_variants.vcf vcf_input
		python3 /usr/TMB/plot_mutation_spectrum.py -v vcf_input -n ${patient}_${T}_${N}
		mv vcf_input/output ./mutation_signature_output
	"""
}

process create_report {
	tag "${patient}"
	publishDir "${params.out_dir}/${patient}_${T}_${N}/report", mode: 'copy', overwrite: true
	
	input:
		file(base_count) 
		file(cds_count)
		tuple val(patient), val(T), val(N), val(dont_use1),
			path(msi_out),
			val(dont_use2),
			path(snv_vcf),
			path(indel_vcf)
		
	output:
		tuple val(patient), val(T), val(N), path("TMB_counts.txt"), path("passed_SNV_AF_counts.txt"),
			path("passed_SNV_coding_AF_counts.txt") 

	script:
	"""
		total_bases=`cat ${base_count}`
		total_CDS_bases=`cat ${cds_count}`

		#python /usr/TMB/allele_counts_from_strelka.py -v ${snv_vcf}
		java -jar /usr/TMB/snpEff//SnpSift.jar extractFields ${snv_vcf} GEN[1].AF > AF.csv
		for i in \$(seq 0 0.1 0.9); do \
			count=`awk -v var=\$i '{if (\$1 > var && \$1 <= var+0.1 ) print \$1}' AF.csv | wc -l`; 
			printf "%s-%s %d\\n"  \$i \$(echo \$i + 0.1 | bc) \$count; 
		done > passed_SNV_AF_counts.txt

		java -Xmx16g -jar /usr/TMB/snpEff/SnpSift.jar filter "(EFF[*].IMPACT = 'MODERATE') | (EFF[*].IMPACT = 'HIGH')" \
			${snv_vcf} | grep -E '^[1234567890XY]{1,2}\\s' | \
			java -jar /usr/TMB/snpEff/SnpSift.jar extractFields - GEN[1].AF > AF_coding.csv

		for i in \$(seq 0 0.1 0.9); do \
			count=`awk -v var=\$i '{if (\$1 > var && \$1 <= var+0.1 ) print \$1}' AF_coding.csv | wc -l`; 
			printf "%s-%s %d\\n"  \$i \$(echo \$i + 0.1 | bc) \$count; 
		done > passed_SNV_coding_AF_counts.txt

		total_SNVs=`cat ${snv_vcf} | grep -v ^# | wc -l`
		total_Indels=`cat ${indel_vcf} | grep -v ^# | wc -l`
		coding_SNVs=`java -Xmx16g -jar /usr/TMB/snpEff/SnpSift.jar filter \
			"(EFF[*].IMPACT = 'MODERATE') | (EFF[*].IMPACT = 'HIGH')" \
			${snv_vcf} | grep -E '^[1234567890XY]{1,2}\\s' | wc -l`
  		coding_Indels=`java -Xmx16g -jar /usr/TMB/snpEff/SnpSift.jar filter \
		  	"(EFF[*].IMPACT = 'MODERATE') | (EFF[*].IMPACT = 'HIGH')" \
			${indel_vcf} | grep -E '^[1234567890XY]{1,2}\\s' | wc -l`
		msi_score=`awk 'NR==2 { print \$NF }' ${msi_out}`
    	
		echo "TMB Pipeline V0.2.1" > TMB_counts.txt
		echo "Config file ${params.samples_file}"
		echo "Patient: ${patient}"
		echo "Tumour: ${T}" >> TMB_counts.txt
		echo "Normal: ${N}" >> TMB_counts.txt
		printf "Non-N bases in 1-22,X,Y: %d\\n" \${total_bases} >> TMB_counts.txt  
		printf "CDS bases in 1-22,X,Y: %d\\n" \${total_CDS_bases} >> TMB_counts.txt  
		printf "Total genome SNVs: %d\\n" \${total_SNVs} >> TMB_counts.txt  
		printf "Total genome Indels: %d\\n" \${total_Indels} >> TMB_counts.txt  
		printf "Coding SNVs: %d\\n" \${coding_SNVs} >> TMB_counts.txt  
		printf "Coding Indels: %d\\n" \${coding_Indels} >> TMB_counts.txt  
		printf "==============================\n" >> TMB_counts.txt  
		printf "Genome SNV TMB: %3.2f\\n" `echo "scale=8; \${total_SNVs} * 1000000 / \${total_bases}" | bc` >> TMB_counts.txt  
		printf "Genome Indel TMB: %3.2f\\n" `echo "scale=8; \${total_Indels} * 1000000 / \${total_bases}" | bc` >> TMB_counts.txt  
		printf "Coding SNV TMB: %3.2f\\n" `echo "scale=8; \${coding_SNVs} * 1000000 / \${total_CDS_bases}" | bc` >> TMB_counts.txt  
		printf "Coding Indel TMB: %3.2f\\n" `echo "scale=8; \${coding_Indels} * 1000000 / \${total_CDS_bases}" | bc` >> TMB_counts.txt  
		printf "MSI score: %3.2f\\n" \${msi_score} >> TMB_counts.txt  
		printf "Report Complete!" >> TMB_counts.txt  
	"""


}