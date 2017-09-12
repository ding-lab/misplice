# MiSplice pipeline #

Pipeline for discovering mutation-induced splice creating events from RNA-Seq

Before running change lines 58,75,78 in MiSplice.pl to reflect personal directories.

Usage: perl MiSplice.pl <run_folder> <step_number>

<run_folder> = full path of the folder holding mutation annotation file(maf) file for all mutations (misplice.input.maf) and file for the sample list (currently named as Samples).

<step_number> run this pipeline step by step (Either run step 1 to step 9 separately) OR run step 11, then steps 7,8,9 separately.

[1] Split mutation annotation file (MAF) into multiple files for processing

[2] Run discovery of mutation induced splice creating events

[3] Generate control MAF

[4] Calculate number of supporting reads for control samples

[5] Combine the read count for case and control samples
        
[6] Calculate the splice score (MaxEntScan)

[7] Split novel sites for bam readcounts 

[8] Do bam readcounts

[9] Run rc, hla and expression filtering 

File format for input MAF (shortened)
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

