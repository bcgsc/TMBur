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
	publishDir "${params.out_dir}/${patient}_bams"
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
	publishDir "${params.out_dir}/${patient}_${T}_${N}"
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
	publishDir "${params.out_dir}/${patient}_${T}_${N}"
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
		tuple val("Manta"), val(patient), val(T), path(T_bam), path(T_bai), val(N), path(N_bam), path(N_bai),
			path("*.vcf.gz"), path("*.vcf.gz.tbi") 
		tuple val(patient), val(T), path(T_bam), path(T_bai), val(N), path(N_bam), path(N_bai),
			path("Manta_${patient}_${T}_vs_${N}.candidateSmallIndels.vcf.gz"),
			path("Manta_${patient}_${T}_vs_${N}.candidateSmallIndels.vcf.gz.tbi") 

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
		mv Manta/results/variants/diploidSV.vcf.gz \
			Manta_${patient}_${T}_vs_${N}.diploidSV.vcf.gz
		mv Manta/results/variants/diploidSV.vcf.gz.tbi \
			Manta_${patient}_${T}_vs_${N}.diploidSV.vcf.gz.tbi
		mv Manta/results/variants/somaticSV.vcf.gz \
			Manta_${patient}_${T}_vs_${N}.somaticSV.vcf.gz
		mv Manta/results/variants/somaticSV.vcf.gz.tbi \
			Manta_${patient}_${T}_vs_${N}.somaticSV.vcf.gz.tbi
 """
}


process strelka {
	tag "${patient}_${T}_${N}"
	publishDir "${params.out_dir}/${patient}_${T}_${N}"
	cpus 48

	input:
		tuple val(patient), val(T), path(T_bam), path(T_bai), val(N), path(N_bam), path(N_bai),
			path("candidateSmallIndels.vcf.gz"),
			path("candidateSmallIndels.vcf.gz.tbi") 
		path(fai_file) 
		path(reference) 

	output:
		tuple val("Strelka"), val(patient), val(T), path(T_bam), path(T_bai), val(N), path(N_bam), path(N_bai),
			path("*.vcf.gz"), path("*.vcf.gz.tbi") 

	script:
	"""
		/usr/TMB/strelka-2.9.2.centos6_x86_64/bin//configureStrelkaSomaticWorkflow.py \
		    --tumor ${T_bam} \
			--normal ${N_bam} \
			--referenceFasta ${reference} \
			--runDir Strelka
		python Strelka/runWorkflow.py -m local -j 48

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