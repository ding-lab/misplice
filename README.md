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

[11] Run steps 1, 2 and 3

[12] Run steps 4, 5, 6, 7, 8 and 9

File format for samples as following (table-delimited with sample id and cancer type):

TCGA-02-0047    gbm

TCGA-02-0055    gbm

TCGA-02-2483    gbm

TCGA-02-2485    gbm

TCGA-02-2486    gbm

File format for Mutation Annotation File


