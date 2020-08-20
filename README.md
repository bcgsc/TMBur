# PROFYLE_TMB

This project is part PROFYLE 
https://www.tfri.ca/profyle

In order to harmonize Tumour Mutation Burden (TMB) estimates across PROFYLE sites,  
a unified analysis pipepline was desired to go from fastq files to TMB estimates using  
software that can be deployed to all of the PROFYLE locations.  

This project uses the following technologies to make the code reprodicible and portable:
+ [nextflow] (https://www.nextflow.io/)
+ [singularity] https://sylabs.io/singularity/

The container and workflow leverage the following bioinformatics tools for analysis:
+ [fastp] (https://github.com/OpenGene/fastp)
+ [BWA] (https://github.com/lh3/bwa)
+ [sambamba] (https://github.com/biod/sambamba)
+ [Manta] (https://github.com/Illumina/manta)
+ [Strelka2] (https://github.com/Illumina/strelka)
+ [MSIsensor2] (https://github.com/niu-lab/msisensor2)
+ [SnpEff] (https://pcingola.github.io/SnpEff/)

  
&nbsp;  
## Getting started

### Install Nextflow
First, you need to be able to run nextflow.  Fortunately, it is easy to install:  
[with details here] (https://www.nextflow.io/),  
but if you already have java 1.8 installed you can just cd to the folder where you want  
to install it and run  
`curl -s https://get.nextflow.io | bash `.  
It is recommended that you install nextflow in a central locaton and not in the same  
folder as the repository files so that you can run this simple test:  
  
To run a quick test of your installation, run this in another folder   
(away from the `nextflow.config` file in this repository):    
```
/path/to/nextflow run hello
```


### Samples CSV file
Next you need to set up a csv file that describes that fastqs to use.  
**Make sure to list the residing path for the fastq files**.  Providing paths  
that are sym-links to files elsewhere on your filesystem may cause problems.  

An example csv file:
```
do_not_use,patient,tissue,read1,read2
1,patient1,T1,291150_S1_L001_R1_001_125bp_158388_chastity_passed.fastq.gz,291150_S1_L001_R2_001_125bp_158388_chastity_passed.fastq.gz
2,patient1,T1,291150_S1_L002_R1_001_125bp_158389_chastity_passed.fastq.gz,291150_S1_L002_R2_001_125bp_158389_chastity_passed.fastq.gz
3,patient1,N1,291647_S2_L003_R1_001_125bp_158387_chastity_passed.fastq.gz,291647_S2_L003_R2_001_125bp_158387_chastity_passed.fastq.gz
```
Note that the patient column is used to match the tumour and normal samples.    
Make sure to use the same tissue id for fastqs that came from the same  
source where tumour samples start with 'T' and normal samples start with 'N'.  
Providing multiple Tumour or Normal samples from the same patient should result  
in a set of analsys results for each possible pair.  

You can set up a csv file that contains many samples if you want to run them  
all at once.

&nbsp;  
## Running the pipeline
In this repository there is a file called `nextflow.config` that contains some   
important parameters to control how the pipeline is run.

Perhaps the only line you need to look at initially is the `executor`, and possibly the  
`queue` line.   The executor instructs nextflow about how to run the commands.   Unless you  
are running analysis in the cloud, you will want to set `executor` to either `local` (to   
run all the commands on the current machine) or `slurm` to run all the jobs using a  
slurm scheduler.   If you use `slurm` you want to run the nextflow command on your cluster   
head node and if you want to use a specific partition on your cluster you need to set   
`queue` appropriately.

To run the pipeline you can start with this command:  
```
nextflow run fastq_to_TMB.nf   
    --samples_file /path/to/samples.csv   
    --reference /projects/rcorbettprj2/mutationalBurden/PROFYLE_container/PROFYLE_tests/hs37d5.fa  
    --out_dir ./TMB_out  
```
**The reference needs to be indexed with BWA.**
`--out_dir` will be where the final results are copied.   During the run a folder called `work`  
will be created in the current directory that contains all of the intermediate files.

Some parameters that you might want to use:  
1. These parameters will create some figures and reports about the resources used during the run  
`-with-report -with-timeline -with-trace `  
2. If you are on slurm and you notice that only 20ish jobs run at a time, try using `-qs 50 ` to allow   
up to 50 jobs to run at once.  
3. `-w /path/to/folder` can be used to have all of the intermidiate analysis files stored somewhere  
other than the `work` folder in the current directory.  
4. `-resume` can be used to pick up analysis where it left off.  This is useful if for some reason something   
crashes at run time, using `-resume` will make sure that the analysis that was done upstream doesn't   
get re-run.   If you want to add another sample or pair of fastq files to the CSV file after you have   
processed other samples `-resume` will ensure that only the new samples get analyzed. 
  
&nbsp;  
## Storage Space
The `work` folder can get quite large. It would be good practise to ensure you have ~1.5Tb per  
sample that you want to analyze.  Once the analysis is complete you can delete the `work`   
folder as the final results files have been copied over to the `--out_dir`.  *The `--out_dir`  
folder contains symlinks to the VCF and BAM files in case you want to copy them from the `work`  
folder before cleaning up.*
  
&nbsp;  
## Output
The output folder, with name format `sample_tumour_normal` will contain a file named `TMB_counts.txt` where   
the following will be reported:

Field | Comment
----- | -------
 Non-N bases in 1-22,X,Y |   Count of the bases used as the whole genome calculation denominator
 CDS bases in 1-22,X,Y |      Count of the unique CDS bases used in the coding calculation denominator
 Total genome SNVs |          Number of passed SNVs called by Strelka 2
 xTotal genome Indels |        Number of passed Indels called by Strelka 2
 Coding SNVs |                 Number of passed Coding SNVs called by Strelka 2
 Coding Indels |               Number of passed Coding Indels called by Strelka 2
 Genome SNV TMB |              total_SNVs * 1000000 / total_bases
 Genome Indel TMB |            total_Indels * 1000000 / total_bases
 Coding SNV TMB |              coding_SNVs * 1000000 / CDS_bases
 Coding Indel TMB |           coding_Indels * 1000000 / CDS_bases
 MSI score |                   Fraction of sites reported as MSI by MSIsensor2
  
&nbsp;  
*BETA*  Variant Allele Fractions
The allele fraction the somatic variants are split into bins ranging in increments of 0.05 and are  
listed in these two files:
1. `passed_SNV_AF_counts.txt`
1. `passed_SNV_coding_AF_counts.txt`
The sum of the counts in each file amount to the counts reported in `TMB_counts.txt`    
These values may be of use when considering clonality for TMB estimates.