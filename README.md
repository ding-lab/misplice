# misplice pipeline #

The pipeline will detect mutation-induced splice events from cancer patients.

Usage: perl misplice.pl <run_folder> <step_number>

<run_folder> = full path of the folder holding maf file for all mutations (currently named as spliceinator.splice.score) and file for the sample list (currently named as Samples).

<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

[1] Split maf

[2] Run discovery for splice creating events

[3] Generate control maf

[4] Calculate # of supporting reads for control samples

[5] Run expression filtering step
        
[6] Calculate the splice score 

[7] Split novel sites for bam readcounts 

[8] Do bam readcounts

[9] Run rc, hla and expression filtering 

File format for samples as following (table-delimited with sample id and cancer type):

TCGA-02-0047    gbm

TCGA-02-0055    gbm

TCGA-02-2483    gbm

TCGA-02-2485    gbm

TCGA-02-2486    gbm

