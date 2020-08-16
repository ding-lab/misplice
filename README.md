# MiSplice pipeline v1.2 #

Pipeline for discovering splice-site creating mutations (SCMs) from RNA-Seq. It works on LSF job scheduler and can run multiple jobs in parallel. 

HG38 reference. For old HG19 reference, please use v1.1 branch


## Citation 

Systematic Analysis of Splice-Site-Creating Mutations in Cancer

Jayasinghe RG*, Cao S*, Gao Q, Wendl MC, Vo NS, Reynolds SM, Zhao Y, Climente-González H, Chai S, Wang F, Varghese R, Huang M, Liang W, Wyczalkowski MA, Sengupta S, Li Z, Payne SH, Fenyö D, Miner JH, Walter MJ, The Cancer Genome Atlas Research Network, Vincent B, Eyras E, Chen K, Shmulevich I, Chen F, Ding L,
Cell Rep 2018 Apr 3;23(1):270-281 

## Install the third-party software

Basic Linux:

sudo apt-get install -y samtools.

With Conda environment:

conda install -c bioconda samtools

## Usage

perl misplice.pl --rdir run_folder --bed bed_file --ref ref_file --q q_name --maf f_maf --bamlist f_bam --step step_number

<run_folder> = full path of the folder holding mutation annotation file (maf) file for all mutations (misplice.input.maf) and file for the sample list (currently named as Samples).

<bed_file> bed file for HG38 refernece: please follow ./resource/GRCh38/work_log_gtf to generate it or download from https://drive.google.com/file/d/1M2yLLjxBwOjUOu-uIxwtUrD_DHTTpioJ/view?usp=sharing

<ref_file> HG38 refernece file
 
<q_name> which bsub quenue for submitting job

<f_maf> Somatic mutation file (in maf format), standard output from vcf2maf script

Example can be found from misplice.example.maf 

<f_bam> file for input bam list

Example can be found from misplice.rnabam.tsv, which contains four columns (sample_name, cancer_type, rna_bam_path, and chr_status). 

chr_status indicates whether the chr in the bam contains chr or not. If it has chr, type chr1 in the four column. Otherwise, type 1.
 

** Make sure the sample_name is consistent with name in the maf file **


<step_number> run this pipeline step by step

[1] Split mutation annotation file (MAF) into multiple files for processing

[2] Run discovery step of mutation induced splice creating events

[3] Generate control MAF

[4] Calculate number of supporting reads for control samples

[5] Combine the read count for case and control samples
        
[6] Calculate the splice score (MaxEntScan)

[7] Split novel sites for bam readcounts 

[8] Do bam readcounts

[9] Run rc, hla and expression filtering 

 


## Developers  

Song Cao, Reyka Jayasinghe, Qingsong Gao, Mike Wendl

Pipeline leading contact: scao@wustl.edu

