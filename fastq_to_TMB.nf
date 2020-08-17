#!/usr/bin/env nextflow

/*
* pipeline input parameters
*/
	params.samples_file = "/projects/rcorbettprj2/mutationalBurden/PROFYLE_container/2p0/test_data/samples.csv"
	params.out_dir = "/projects/rcorbettprj2/mutationalBurden/PROFYLE_container/2p0/test_data/output"
	params.annotation = "GRCh37.75"
	params.reference = "/projects/alignment_references/Homo_sapiens/hg19a/genome/bwa_64/hg19a.fa"
	
	log.info """\
	TMB estimation pipeline
	===================================
	samples_file : ${params.samples_file}
	out_dir      : ${params.out_dir}
	annotations  : ${params.annotation}
	reference    : ${params.reference}
	"""
	.stripIndent()

/*
binaries
*/

/*
bbmerge="~/bin/bbmap/bbmerge.sh"
bwa="/gsc/software/linux-x86_64-centos7/bwa-0.7.17/bwa"
skewer="/gsc/software/linux-x86_64-centos6/skewer-0.1.127/skewer"
sambamba="/gsc/software/linux-x86_64/sambamba-0.5.5/sambamba_v0.5.5"
//sambamba="/home/rcorbett/bin/sambamba-0.7.1-linux-static"  //better error reporting than 0.5.5, but MUCH slower
manta="/projects/vleblanc_prj/tools/manta-1.6.0.centos6_x86_64/bin/configManta.py"
strelka="/projects/da_workspace/DA-228/strelka-2.9.2.centos6_x86_64/bin//configureStrelkaSomaticWorkflow.py"
python="~/python-2.7-venv/bin/python"
snpSift="/home/rcorbett/bin/snpEff/snpEff/SnpSift.jar"
snpEff="/home/rcorbett/bin/snpEff/snpEff.jar"
seqtk="/home/rcorbett/bin/seqtk/seqtk-1.3/seqtk"
bedtools="/home/rcorbett/bin/bedtools2//bin/bedtools"
msisensor2="/gsc/software/linux-x86_64-centos7/msisensor2-0.1/msisensor2"
msisensor_resources="/gsc/resources/pipeline/msisensor/models_b37_HumanG1Kv37/1030c0aa35ca5c263daeae866ad18632"
fastp="/home/rcorbett/bin/fastp"
*/

//Load in the samples file
Channel
	.fromPath(params.samples_file)
	.splitCsv(header:true)
	.map{ row-> tuple(file(row.read1).baseName, row.patient, row.tissue, file(row.read1), file(row.read2)) }
	.into {samples_ch; samples_ch2; samples_4print }

samples_4print.subscribe onNext: { println "[samples] $it"}

//value channels can be re-used as many times as we like
ref_ch = Channel.value(file(params.reference))
annotation_ch = Channel.value(params.annotation)
bwa_index_ch = Channel.value(file("${params.reference}.*"))
fasta_index_ch = Channel.value(file("${params.reference}.fai"))

process count_fasta_bases {
	tag "$fasta"
	input:
		file(fasta) from ref_ch

	output:
		file("base_count.txt") into base_count_ch

	script:
	"""
		seqtk comp ${fasta} | grep -E '^[123456789XY]{1,2}\\s' | awk '{ print \$3+\$4+\$5+\$6 }' | paste -sd+ | bc > "base_count.txt"
	"""
}

process count_CDS_bases {
	tag "$anno"
	input:
		val(anno) from annotation_ch

	output:
		file("CDS_size.txt") into CDS_count_ch

	script:
	"""
		java -jar /usr/TMB/snpEff/snpEff.jar  dump -v -bed $anno > ${anno}.bed
		grep -E '^[1234567890XY]{1,2}\\s' ${anno}.bed | grep CDS | bedtools sort | bedtools merge | awk '{ print \$3-\$2 }' | paste -sd+ | bc > CDS_size.txt
	"""
}

process trim_pair {
    tag "$an_id"
	cpus 10
	publishDir "${params.out_dir}", mode: 'copy', overwrite: true

	input:
		tuple an_id, patient, tissue, file(reads1), file(reads2) from samples_ch

	output:
        tuple an_id, patient, tissue, file("${reads1}.fastp.trimmed.fq.gz"),
			file("${reads2}.fastp.trimmed.fq.gz") into trimmed_reads_ch
        file("${an_id}.fastp.json") into paired_trimmed_report_ch
		
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

	input:
		tuple an_id, patient, tissue, file(trim1), file(trim2) from trimmed_reads_ch
		file(reference_fasta) from ref_ch
		file(bwa_index) from bwa_index_ch

	output:
		tuple an_id, patient, tissue, file("${an_id}.sam") into sam_file_ch

	script:
	"""
		bwa mem ${reference_fasta} -t 48 ${trim1} ${trim2} > ${an_id}.sam
	"""
}

(sam_file_ch1, sam_file_ch2) = sam_file_ch.into(2)
//sam_file_ch2.subscribe onNext: { println "[sam_file] ${it}" }

process sam_to_bam {
	tag "$an_id"
	cpus 10

	input:
		tuple an_id, patient, tissue, file(sam_file) from sam_file_ch1

	output:
		tuple patient, tissue, an_id, file("${an_id}.bam") into bam_file_ch

	script:
 	"""
		sambamba view -S -f bam ${sam_file} > ${an_id}.bam
 	"""
}

process sort_bam {
	tag "$an_id"
	cpus 10

	input:
		tuple patient, tissue, an_id, file(unsorted_bam_file) from bam_file_ch

	output:
		tuple patient, tissue, an_id, file("${unsorted_bam_file.baseName}.sorted.bam") into sorted_bam_file_ch

	script:
	"""
		sambamba sort --tmpdir . ${unsorted_bam_file}
	"""
}

sorted_bam_file_ch.into { sorted_bam_file_ch1; sorted_bam_file_ch2 }
sorted_bam_file_ch1.subscribe onNext: { println "[sorted bam file] : $it" }

// split into 2 channels, one for merging bams, the other uses bams from one
// fastq pair
single_bam_ch = Channel.create()
multiple_bam_ch = Channel.create()

sorted_bam_file_ch2.groupTuple(by:[0, 1])
	.choice(single_bam_ch, multiple_bam_ch) {it[3].size() > 1 ? 1 : 0}

single_bam_ch = single_bam_ch.map {
	patient, tissue, an_id, bam -> [patient, tissue, bam]
}

process merge_bam {
	tag "${patient}_${tissue}"
	memory "64 GB"
	cpus 16

	input:
		tuple patient, tissue, an_id, file(bams_to_merge) from multiple_bam_ch

	output:
		tuple patient, tissue, file("${patient}-${tissue}_merged.bam") into merged_bams_ch

	script:
	"""
		sambamba merge ${patient}-${tissue}_merged.bam $bams_to_merge
	"""
}

//concatenate together the bams that required merging with the bams that were
//single pairs of fastqs only
merged_bams_ch.mix(single_bam_ch).into {bams_for_analysis_ch; bams_for_analysis_ch2 }
bams_for_analysis_ch2.subscribe onNext: { println "[bams_to_analyze] $it" }

process mark_duplicates {
	tag "${patient}_${tissue}"
	publishDir "${params.out_dir}/${patient}_bams"
	memory "128 GB"
	cpus 16

	input:
		tuple patient, tissue, file(bam_file) from bams_for_analysis_ch

	output:
		tuple patient, tissue, file("${bam_file.baseName}.dup.bam"),
		file("${bam_file.baseName}.dup.bam.bai") into duped_bams_ch

	script:
	"""
		sambamba markdup \
		--nthreads 8 \
		--tmpdir /var/tmp \
		--hash-table-size 5000000 \
		--overflow-list-size 5000000 \
		--io-buffer-size 1024 ${bam_file} ${bam_file.baseName}.dup.bam
	"""
}


//get the cross product of the tumours and normals per patient
duped_bams_ch.branch{
	tumour : it[1] =~ /^T.*/
	normal : it[1] =~ /^N.*/
	}.set { split_tissues }

split_tissues.tumour.view { "[T/N] $it is a tumour" }.set {tumours_ch}
split_tissues.normal.view { "[T/N] $it is a normal" }.set {normals_ch}

tumours_ch.combine(normals_ch, by : 0).view { println "[crossed] $it" }
	.into { bams_for_somatic_ch; bams_for_somatic_ch2; bams_for_somatic_ch3 }

// much of this is pulled from
// https://github.com/nf-core/sarek/blob/master/main.nf
process manta {
	tag "${patient}_${T}_${N}"
	publishDir "${params.out_dir}/${patient}_${T}_${N}"
	cpus 48

	input:
		tuple patient, T, file(T_bam), file(T_bai), N, file(N_bam),
			file(N_bai) from bams_for_somatic_ch
		file(fai_file) from fasta_index_ch
		file(reference) from ref_ch

	output:
		tuple val("Manta"), patient, T, file(T_bam), file(T_bai), N, file(N_bam), file(N_bai),
			file("*.vcf.gz"), file("*.vcf.gz.tbi") into Manta_results_ch
		tuple patient, T, file(T_bam), file(T_bai), N, file(N_bam), file(N_bai),
			file("Manta_${patient}_${T}_vs_${N}.candidateSmallIndels.vcf.gz"),
			file("Manta_${patient}_${T}_vs_${N}.candidateSmallIndels.vcf.gz.tbi") into manta_for_strelka_ch

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
		tuple patient, T, file(T_bam), file(T_bai), N, file(N_bam), file(N_bai),
			file(small_indels_from_manta),
			file(small_indels_from_manta_index) from manta_for_strelka_ch
		file(fai_file) from fasta_index_ch
		file(reference) from ref_ch

	output:
		tuple val("Strelka"), patient, T, file(T_bam), file(T_bai), N, file(N_bam), file(N_bai),
			file("*.vcf.gz"), file("*.vcf.gz.tbi") into strelka_results_ch

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

//combine channels so all the VCF files for annotation are in 1 place
strelka_results_ch.into{ strelka_for_indels_ch; strelka_for_snvs_ch}

//MSIsensor2
process MSIsensor2 {
	tag "${patient}_${T}_${N}"
	publishDir "${params.out_dir}/${patient}_${T}_${N}"
	cpus 32

	input:
		tuple patient, T, file(T_bam), file(T_bai), N, file(N_bam),
			file(N_bai) from bams_for_somatic_ch2
		
	output:
		tuple val("msisensor2"), patient, T, N, 
			file("msisensor2_${patient}_${T}_${N}.txt") into msi_results_ch

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

vcf_to_annotate = Channel.empty().mix(
	Manta_results_ch.map {  //somatic SV
		tool, patient, T, T_bam, T_bai, N, N_bam, N_bai, vcf_files, vcf_index_files ->
			["manta_SV", patient, T, N, vcf_files[3]]
	},
	strelka_for_indels_ch.map { //strela indels
		tool, patient, T, T_bam, T_bai, N, N_bam, N_bai, vcf_files, vcf_index_files ->
			["strelka_Indel", patient, T, N, vcf_files[0]]
	},
	strelka_for_snvs_ch.map { // strelka SNVs
		tool, patient, T, T_bam, T_bai, N, N_bam, N_bai, vcf_files, vcf_index_files ->
			["strelka_SNV", patient, T, N, vcf_files[1]]
	}
)

//vcf_to_annotate.view { "[annotation_vcfs] $it" }

//make passed variant VCF files
process get_passed_variants {
	tag "${patient}_${tool}_${T}_${N}"

	input:
		tuple tool, patient, T, N, file(vcf) from vcf_to_annotate

	output:
		tuple tool, patient, T, N, file("${vcf.baseName}.pass.vcf") into vcf_pass_ch

	script:
	"""
		java -Xmx2g -jar /usr/TMB/snpEff/SnpSift.jar filter "(FILTER = 'PASS')" ${vcf} > ${vcf.baseName}.pass.vcf
	"""
}

// annotate with SNPEff
process annotate {
	tag "${patient}_${tool}_${T}_${N}"
	publishDir "${params.out_dir}/${patient}_${T}_${N}"

	memory '48 GB'

	input:
		tuple tool, patient, T, N, file(vcf) from vcf_pass_ch

	output:
		tuple tool, patient, T, N, file("${vcf.baseName}.snpEff.vcf") into vcf_annotated_ch

	script:
	"""
		java -Xmx48g -jar /usr/TMB/snpEff/snpEff.jar GRCh37.75 -s somatic.snvs.PASS.vcf ${vcf} > ${vcf.baseName}.snpEff.vcf
	"""
}

vcf_annotated_ch.into { annotated_vcf_ch; for_print_ch}
for_print_ch.view{ println "[vcfs] $it"}

//join together the annotated VCFs for each of the unique Patient, T, N combinations
//using the sort option after join to keep the VCFs in a consistent order.
//This order is used to create the report in print_report.
//Dumped the tool name so it doesn't clash with the filename sort order
annotated_vcf_ch.mix(msi_results_ch)
	.groupTuple(by:[1,2,3])
	.map { tool, patient, T, N, vcf -> [patient, T, N, vcf.sort{it.simpleName} ] }
	.into{vcf_for_report_ch; vcf_for_print} 
vcf_for_print.view{ println "[collapsed] $it"}


process print_report {
	tag "${patient}"
	publishDir "${params.out_dir}/${patient}_${T}_${N}", mode: 'copy', overwrite: true
	
	input:
		file(base_count) from base_count_ch
		file(cds_count) from CDS_count_ch
		tuple patient, T, N, file(vcf_files) from vcf_for_report_ch
		/* vcf files are pre-sorted to have these indices
		0 - Manta SV
		1 - Strelka Indel
		2 - Strelka SNV
		3 - MSIsensor report
		*/

	output:
		tuple patient, T, N, file("TMB_counts.txt") into final_results_ch

	script:
	"""
		total_bases=`cat ${base_count}`
		total_CDS_bases=`cat ${cds_count}`

		total_SNVs=`cat ${vcf_files[2]} | grep -v ^# | wc -l`
		total_Indels=`cat ${vcf_files[1]} | grep -v ^# | wc -l`
		coding_SNVs=`java -Xmx16g -jar /usr/TMB/snpEff/SnpSift.jar filter \
			"(EFF[*].IMPACT = 'MODERATE') | (EFF[*].IMPACT = 'HIGH')" \
			${vcf_files[2]} | grep -E '^[1234567890XY]{1,2}\\s' | wc -l`
  		coding_Indels=`java -Xmx16g -jar /usr/TMB/snpEff/SnpSift.jar filter \
		  	"(EFF[*].IMPACT = 'MODERATE') | (EFF[*].IMPACT = 'HIGH')" \
			${vcf_files[1]} | grep -E '^[1234567890XY]{1,2}\\s' | wc -l`
		msi_score=`awk 'NR==2 { print \$NF }' ${vcf_files[3]}`
    	
		echo "TMB Pipeline V0.2" > TMB_counts.txt
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

//will run at the end of the analysis.
workflow.onComplete {
	final_results_ch.view{ println "[complete] $it"}
	println "Pipeline $workflow.scriptName completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
