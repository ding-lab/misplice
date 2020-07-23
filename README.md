# MiSplice pipeline v1.2 #

HG38 reference, for old HG19 version, please go v1.1 branch

Song Cao, Reyka Jayasinghe, Qingsong Gao, Mike Wendl

Pipeline for discovering splice-site creating mutations (SCMs) from RNA-Seq. It works on LSF job scheduler and can run multiple jobs in parallel.
 
Please cite the following paper if you use MiSplice pipeline:

Systematic Analysis of Splice-Site-Creating Mutations in Cancer

Jayasinghe RG*, Cao S*, Gao Q, Wendl MC, Vo NS, Reynolds SM, Zhao Y, Climente-González H, Chai S, Wang F, Varghese R, Huang M, Liang W, Wyczalkowski MA, Sengupta S, Li Z, Payne SH, Fenyö D, Miner JH, Walter MJ, The Cancer Genome Atlas Research Network, Vincent B, Eyras E, Chen K, Shmulevich I, Chen F, Ding L,
Cell Rep 2018 Apr 3;23(1):270-281


Usage: perl misplice.pl --rdir run_folder --bed bed_file --ref ref_file --q q_name --step <step_numuber>

<run_folder> = full path of the folder holding mutation annotation file (maf) file for all mutations (misplice.input.maf) and file for the sample list (currently named as Samples).

<bed_file> bed file for HG38 refernece: please follow ./resource/GRCh38/work_log_gtf to generate it or download from https://drive.google.com/file/d/1M2yLLjxBwOjUOu-uIxwtUrD_DHTTpioJ/view?usp=sharing

<ref_file> HG38 refernece file
 
<q_name> which bsub quenue for submitting job

<step_number> run this pipeline step by step (Either run step 1 to step 9 separately) OR run step 11, then steps 7,8,9 separately.

[1] Split mutation annotation file (MAF) into multiple files for processing

[2] Run discovery step of mutation induced splice creating events

[3] Generate control MAF

[4] Calculate number of supporting reads for control samples

[5] Combine the read count for case and control samples
        
[6] Calculate the splice score (MaxEntScan)

[7] Split novel sites for bam readcounts 

[8] Do bam readcounts

[9] Run rc, hla and expression filtering 

File format for input MAF (shortened):
 
     1	GENE	VPS13D
    
	 2	CHR	1

     3	START	12418559

     4	STOP	12418559

     5	STRAND	+

     6	MUTATION TYPE	Missense_Mutation

     7	TYPE	SNP

     8	REFALLELE	C

     9	TUMORALLELE1	T

    10	TUMORALLELE2	T

    11	TUMOR SAMPLE	TCGA-02-0047-01A-01D-1490-08

    12	NORMAL SAMPLE	TCGA-02-0047-10A-01D-1490-08

    13	CANCER TYPE gbm

File format for samples as following (table-delimited with sample id and cancer type):
	TCGA-02-0047    gbm

	TCGA-02-0055    gbm

	TCGA-02-2483    gbm

	TCGA-02-2485    gbm

	TCGA-02-2486    gbm

