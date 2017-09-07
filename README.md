# MiSplice pipeline #

Pipeline for discovering mutation-induced splice creating events from RNA-Seq

Before running change lines 58,75,78 in MiSplice.pl to reflect personal directories.

Usage: perl MiSplice.pl <run_folder> <step_number>

<run_folder> = full path of the folder holding mutation annotation file(maf) file for all mutations (misplice.input.maf) and file for the sample list (currently named as Samples).

<step_number> run this pipeline step by step. (running the whole pipeline if step number is 11, once 11 completes run filtering steps (step 12))

[1] Split mutation annotation file (MAF) into multiple files for processing

[2] Run discovery of mutation induced splice creating events

[3] Generate control MAF

[4] Calculate number of supporting reads for control samples

[5] Combine the read count for case and control samples
        
[6] Calculate the splice score (MaxEntScan)

[7] Split novel sites for bam readcounts 

[8] Do bam readcounts

[9] Run rc, hla and expression filtering 

File format for samples as following (table-delimited with sample id and cancer type):

TCGA-02-0047    gbm

TCGA-02-0055    gbm

TCGA-02-2483    gbm

TCGA-02-2485    gbm

TCGA-02-2486    gbm

File format for Mutation Annotation File (MAF)
     1	VPS13D
     2	0
     3	.
     4	GRCh37
     5	1
     6	12418559
     7	12418559
     8	+
     9	Missense_Mutation
    10	SNP
    11	C
    12	T
    13	T
    14	novel
    15
    16	TCGA-02-0047-01A-01D-1490-08
    17	TCGA-02-0047-10A-01D-1490-08
