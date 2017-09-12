######### Song Cao###########
## Pipeline for discovering mutation-induced splice creating events from RNA-Seq ##
#   MiSplice.pl #
### updated date: 04/24/2017 ###
### 04/26/2017 ####

#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;

my $version = 1.0;
#color code
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m";
#usage information

(my $usage = <<OUT) =~ s/\t+//g;

Pipeline for discovering mutation-induced splice creating events from RNA-Seq

Refer to github page for more info: https://github.com/ding-lab/misplice

Pipeline version: $version
$yellow     Usage: perl $0 <run_folder> <step_number> $normal

<run_folder> = full path of the folder holding files for this sequence run

<step_number> run this pipeline step by step. (Either run step 1 to step 9 separately) OR run step 11, then steps 7,8,9 separately.

$red	[1] Split maf
	[2] Run discovery for splice creating events

$gray	[3] Generate control maf
	[4] Calculate # of supporting reads for control samples

$green	[5] Generate table for # of supporting reads in case and control 
		
$purple	[6] Calculate the splice score 

$cyan		[7] Split novel sites for bam readcounts 
	[8] Do bam readcounts
	[9] Run rc, hla and expression filtering
$normal
OUT

die $usage unless @ARGV == 2;
my ( $run_dir, $step_number ) = @ARGV;
if ($run_dir =~/(.+)\/$/) {
    $run_dir = $1;
}
die $usage unless (($step_number >=0)&&($step_number <= 9) || ($step_number==11) || ($step_number==12));
my $email = "reyka\@wustl\.edu";
# everything else below should be automated
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-2];
my $HOME1="/gscmnt/gc2509/dinglab/reyka/MC3_splice";
#store job files here
if (! -d $HOME1."/tmp_misplice") {
    `mkdir $HOME1"/tmp_misplice"`;
}
my $job_files_dir = $HOME1."/tmp_misplice";
#store SGE output and error files here
if (! -d $HOME1."/LSF_DIR_MISPLICE") {
    `mkdir $HOME1"/LSF_DIR_MISPLICE"`;
}

my $file_number_of_spliting = 200;

my $lsf_file_dir = $HOME1."/LSF_DIR_MISPLICE" or die;

### USER can define them ###  

my $script_dir = "/gscmnt/gc2706/dinglab/medseq/LabCode/Reyka/misplice";
my $fmaf = $run_dir."/misplice.input.maf";
my $rcbam = $script_dir."/resource/bampath.txt";

#############

#my $run_script_path = `dirname $0`;

my $run_script_path=$script_dir;
chomp $run_script_path;

$run_script_path = "/usr/bin/perl ".$run_script_path."/";
#print $run_script_path,"\n";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";
$current_job_file="";

if($step_number==11)
 	{
	&bsub_maf_split(); 
	&bsub_job_array_ns();
	&bsub_control();
	} elsif($step_number==12)
	{
	&bsub_array_control();
	&bsub_support_reads_table();
	&bsub_splice_score();
	&bsub_novels_split();
	&bsub_job_array_rc();
	&bsub_rc_hla_filter();	
	} elsif ($step_number == 1) 
	{
     &bsub_maf_split(1);
    }elsif ($step_number == 2)
    {
     &bsub_job_array_ns(1);
    }elsif ($step_number == 3)
    {
     &bsub_control(1);
    }elsif ($step_number == 4)
    {
     &bsub_array_control(1);
    }elsif ($step_number == 5)
    {
     &bsub_support_reads_table(1);
    }elsif ($step_number == 6)
    {
     &bsub_splice_score(1);
    } elsif ($step_number == 7)
    {
     &bsub_novels_split(1);
    }elsif ($step_number == 8)
    {
     &bsub_job_array_rc(1);
    }elsif ($step_number == 9)
    {
     &bsub_rc_hla_filter(1);
	}

sub bsub_maf_split{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";

    $current_job_file = "j1_maf_split.sh";
	
	open(NS, ">$job_files_dir/$current_job_file") or die $!;

    print NS "#!/bin/bash\n";
    print NS "#BSUB -n 1\n";
    print NS "#BSUB -R \"rusage[mem=30000]\"","\n";
    print NS "#BSUB -M 30000000\n";
    print NS "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print NS "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print NS "#BSUB -J $current_job_file\n";
	print NS "MAF=".$run_dir."/misplice.input.maf\n";
	print NS "DIR_OUT=".$run_dir."/NS_CASE\n";
	print NS "if [ ! -d \${DIR_OUT} ]\n";
    print NS "then\n";
	print NS "mkdir \${DIR_OUT}\n";
	print NS "fi\n";
	print NS "     ".$run_script_path."split_maf.pl \${MAF} \${DIR_OUT} $file_number_of_spliting"."\n";  
	close NS;
	$bsub_com = "bsub < $job_files_dir/$current_job_file";	
	#$bsub_com = "qsub -V -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";

	system ($bsub_com);

	}	

sub bsub_job_array_ns {
	my ($step_by_step) = @_;
	if ($step_by_step) {
		$hold_job_file = "";
	}else{
		$hold_job_file = $current_job_file;
	}
	$current_job_file = "j2"."_array_ns".".sh";
	open (ANS, ">$job_files_dir/$current_job_file") or die $!;
	print ANS "#!/bin/bash\n";
	print ANS "#BSUB -n 1\n";
	#print RM "#BSUB -q ding-lab\n";
	print ANS "#BSUB -R \"span[hosts=1] rusage[mem=10000]\"","\n";
    #print RM "#BSUB -R \"rusage[mem=40000]\"","\n";
    print ANS "#BSUB -M 10000000\n";
    print ANS "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print ANS "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print ANS "#BSUB -J $current_job_file\[1-$file_number_of_spliting\]\n";
	print ANS "#BSUB -w \"$hold_job_file\"","\n";	
	#print ANS "#BSUB -q ding-lab\n";
	#####################
	#print RM "#\$ -t 1-$file_number_of_RepeatMasker:1","\n";
	print ANS "ANS_OUT=".$run_dir."/NS_CASE/misplice.input.maf.".'${LSB_JOBINDEX}'.".v2.filtered.5\n";
	print ANS "ANS_IN=".$run_dir."/NS_CASE/misplice.input.maf.".'${LSB_JOBINDEX}'."\n";
	#print ANS 'if [ ! -f $ANS_OUT]',"\n";
    print ANS "if [ ! -f \${ANS_OUT} ]\n";
	print ANS "then\n";
	print ANS "     ".$run_script_path."in_silico_ns.v4.pl \${ANS_IN} \${ANS_OUT} $script_dir"."\n";
	print ANS "fi\n";
	close ANS;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
	system ($bsub_com);
}

sub bsub_control {

  	my ($step_by_step) = @_;

    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j3_control.sh";

    open(CONTR, ">$job_files_dir/$current_job_file") or die $!;
    print CONTR "#!/bin/bash\n";
    print CONTR "#BSUB -n 1\n";
    print CONTR "#BSUB -R \"rusage[mem=30000]\"","\n";
    print CONTR "#BSUB -M 30000000\n";
    print CONTR "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print CONTR "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print CONTR "#BSUB -J $current_job_file\n";
	print CONTR "#BSUB -w \"$hold_job_file\"","\n";
    #print CONTR "#BSUB -q ding-lab\n";
	print CONTR "DCONTR=".$run_dir."/Controls\n"; 
    print CONTR "FS=".$run_dir."/Samples\n";	
	print CONTR "FJ=".$run_dir."/novel.junctions.filtered\n";
	print CONTR "FOUT=".$run_dir."/Controls/novel.junctions.filtered.controls\n";
	print CONTR "if [ ! -d \${DCONTR} ]\n";
	print CONTR "then\n";
	print CONTR "mkdir \${DCONTR}\n"; 
	print CONTR "fi\n";
	print CONTR "cat $run_dir/NS_CASE/*.filtered.5 > \${FJ}\n";
	print CONTR "     ".$run_script_path."controls.v2.pl $fmaf \${FJ} \${FOUT} \${FS}\n";
	close CONTR; 
	$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);
}


sub bsub_array_control {

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }
    $current_job_file = "j4_array_control.sh";

	my $f_contr=$run_dir."/Controls/novel.junctions.filtered.controls";
	my $cout=`wc -l $f_contr`;  
	my @temp=split(" ",$cout);
 
	my $NUM_CONTR=$temp[0];
	
	open(IN,"<$f_contr"); 
	my $cc=1; 
	while(<IN>)
	{
		my $line=$_; 
		chomp($line); 
		my $f_out=$f_contr.".$cc"; 
		open(OUT,">$f_out");
		print OUT $line,"\n";
		close OUT;  
		$cc++;
	}
	close IN; 

    open(CONTRS, ">$job_files_dir/$current_job_file") or die $!;

    print CONTRS "#!/bin/bash\n";
    print CONTRS "#BSUB -n 1\n";
    print CONTRS "#BSUB -R \"rusage[mem=30000]\"","\n";
    print CONTRS "#BSUB -M 30000000\n";
    print CONTRS "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print CONTRS "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print CONTRS "#BSUB -J $current_job_file\[1-$NUM_CONTR\]\n";
    print CONTRS "#BSUB -w \"$hold_job_file\"","\n";
    #print CONTRS "#BSUB -q ding-lab\n"; 
	print CONTRS "CONT_OUT=".$run_dir."/Controls/novel.junctions.filtered.controls.".'${LSB_JOBINDEX}'.".v2.filtered.5\n";
    print CONTRS "CONT_IN=".$run_dir."/Controls/novel.junctions.filtered.controls.".'${LSB_JOBINDEX}'."\n";    
    print CONTRS "     ".$run_script_path."in_silico_ns.control.v4.pl \${CONT_IN} \${CONT_OUT} $script_dir\n";
    close CONTRS;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);

}


sub bsub_support_reads_table{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }
    $current_job_file = "j5_support_reads_table.sh";
    my $f_contr=$run_dir."/Controls/novel.junctions.filtered.controls";
    my $cout=`wc -l $f_contr`;
    my @temp=split(" ",$cout);
    my $NUM_CONTR=$temp[0];

    open(EXPF, ">$job_files_dir/$current_job_file") or die $!;
    print EXPF "#!/bin/bash\n";
    print EXPF "#BSUB -n 1\n";
    print EXPF "#BSUB -R \"rusage[mem=30000]\"","\n";
    print EXPF "#BSUB -M 30000000\n";
    print EXPF "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print EXPF "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
	print EXPF "#BSUB -J $current_job_file\n";
    print EXPF "#BSUB -w \"$hold_job_file\"","\n";
   # print EXPF "#BSUB -q ding-lab\n";
	print EXPF "FS=".$run_dir."/Samples\n";
	print EXPF "OUTF=".$run_dir."/case.control.distributionmethod\n";
	print EXPF "OUTF_K=".$run_dir."/case.control.distributionmethod.withkey\n";
	print EXPF "SUFFIX="."v2.filtered.5\n";
	#print EXPF "     ".$run_script_path."case_control_rc_5.pl \${FS} $f_contr \${SUFFIX} $NUM_CONTR > \${OUTF}\n";
	print EXPF "     ".$run_script_path."case_control_rc_5.pl $f_contr \${SUFFIX} $NUM_CONTR > \${OUTF}\n";
	print EXPF "cat \${OUTF} | awk \'{print substr(\$3,1,12)\"_\"\$5\"\t\"\$0}\' > \${OUTF_K}\n";
    close EXPF;

    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);

 }

	


sub bsub_splice_score {

    my ($step_by_step) = @_;

    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j6_splice_score.sh";

    open(SSCORE, ">$job_files_dir/$current_job_file") or die $!;
    print SSCORE "#!/bin/bash\n";
    print SSCORE "#BSUB -n 1\n";
    print SSCORE "#BSUB -R \"rusage[mem=30000]\"","\n";
    print SSCORE "#BSUB -M 30000000\n";
    print SSCORE "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print SSCORE "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print SSCORE "#BSUB -J $current_job_file\n";
    print SSCORE "#BSUB -w \"$hold_job_file\"","\n";
    #print SSCORE "#BSUB -q ding-lab\n";
	print SSCORE "DIRC=".$run_dir."/NS_CASE\n";
	print SSCORE "INFO1=*.v2.filtered.5.detailed.alignment.5\n";
	print SSCORE "INFO2=misplice.input.maf.*.v2.filtered.5\n";
	print SSCORE "     ".$run_script_path."prefilter.pl \${DIRC} \${INFO1}\n";
	print SSCORE "     ".$run_script_path."splice_score_novel.pl $script_dir \${DIRC} \${INFO2}\n";
	print SSCORE "mv \${DIRC}/novel.splice.scores $run_dir/novel.splice.scores\n";
	close SSCORE; 	

    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);

	}


sub bsub_novels_split{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";

    my ($step_by_step) = @_;

    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j7_novels_split.sh";

    open(NS2, ">$job_files_dir/$current_job_file") or die $!;

    print NS2 "#!/bin/bash\n";
    print NS2 "#BSUB -n 1\n";
    print NS2 "#BSUB -R \"rusage[mem=30000]\"","\n";
    print NS2 "#BSUB -M 30000000\n";
    print NS2 "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print NS2 "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print NS2 "#BSUB -J $current_job_file\n";
    #print NS2 "#BSUB -q ding-lab\n";
    print NS2 "MAF=".$run_dir."/novel.splice.scores\n";
    print NS2 "DIR_OUT=".$run_dir."/RC\n";
    print NS2 "if [ ! -d \${DIR_OUT} ]\n";
    print NS2 "then\n";
    print NS2 "mkdir \${DIR_OUT}\n";
    print NS2 "fi\n";
    print NS2 "     ".$run_script_path."split_novels.pl \${MAF} \${DIR_OUT} $file_number_of_spliting"."\n";
    close NS2;
    $bsub_com = "bsub < $job_files_dir/$current_job_file";
    #$bsub_com = "qsub -V -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";

    system ($bsub_com);

    }

sub bsub_job_array_rc {
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }
    $current_job_file = "j8"."_array_rc".".sh";
    open (RC, ">$job_files_dir/$current_job_file") or die $!;
    print RC "#!/bin/bash\n";
    print RC "#BSUB -n 1\n";
    #print RM "#BSUB -q ding-lab\n";
    print RC "#BSUB -R \"span[hosts=1] rusage[mem=10000]\"","\n";
    #print RM "#BSUB -R \"rusage[mem=40000]\"","\n";
    print RC "#BSUB -M 10000000\n";
    print RC "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print RC "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print RC "#BSUB -J $current_job_file\[1-$file_number_of_spliting\]\n";
    print RC "#BSUB -w \"$hold_job_file\"","\n";
    #print RC "#BSUB -q ding-lab\n";
    #####################
    print RC "RC_OUT=".$run_dir."/RC/novel.splice.scores.".'${LSB_JOBINDEX}'.".rc\n";
    print RC "RC_IN=".$run_dir."/RC/novel.splice.scores.".'${LSB_JOBINDEX}'."\n";
    print RC "     ".$run_script_path."bamreadcount_RNA_expfile_samtools.pl \${RC_IN} $rcbam > \${RC_OUT}"."\n";
    #print ANS "fi\n";
    close RC;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);
}

sub bsub_rc_hla_filter{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";

    my ($step_by_step) = @_;

    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j9_rc_hla_filter.sh";

    open(RCHLA, ">$job_files_dir/$current_job_file") or die $!;

    print RCHLA "#!/bin/bash\n";
    print RCHLA "#BSUB -n 1\n";
    print RCHLA "#BSUB -R \"rusage[mem=30000]\"","\n";
    print RCHLA "#BSUB -M 30000000\n";
    print RCHLA "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print RCHLA "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print RCHLA "#BSUB -J $current_job_file\n";
   # print RCHLA "#BSUB -q ding-lab\n";
    print RCHLA "RCT=".$run_dir."/novel.splice.scores.rc\n";
	print RCHLA "RCTK=".$run_dir."/novel.splice.scores.rc.key\n";
	print RCHLA "RCTC=".$run_dir."/novel.splice.scores.rc.key.combined\n";
	print RCHLA "HLA=".$run_dir."/novel.splice.scores.rc.key.combined.noHLA\n";
	print RCHLA "HLAVAF=".$run_dir."/novel.splice.scores.rc.key.combined.noHLA.vaf\n";
	print RCHLA "HLAVAFH=".$run_dir."/novel.splice.scores.rc.key.combined.noHLA.vaf.highexp\n";
	print RCHLA "EXPK=".$run_dir."/case.control.distributionmethod.withkey\n";
	print RCHLA "MAFK=".$run_dir."/maf.key\n";
	print RCHLA "FINAL=".$run_dir."/novel.splice.scores.rc.key.combined.vaf.noHLA.highexp.maf\n";
	print RCHLA "cat $run_dir/RC/novel.splice.scores.*.rc > \${RCT}\n";
	print RCHLA "awk \'{print substr(\$1,1,12)\"_\"\$2\"_\"\$3\"_\"\$4\"_\"\$5\"\\t\"\$0}\' \${RCT} > \${RCTK}\n";
	print RCHLA "     ".$run_script_path."combine.pl \${RCTK} \${EXPK} 0 0 > \${RCTC}\n";
	print RCHLA "grep -v \'HLA-\' \${RCTC} > \${HLA}\n";
	print RCHLA "awk -F\' \' \'{ tmp=(\$21/(\$14))*100; print \$1\"\\t\"tmp\"\\t\"\$0}\' \${HLA}|awk \'\$2>5\' > \${HLAVAF}\n"; 
	print RCHLA "grep high_expression \${HLAVAF} > \${HLAVAFH}\n";
#	print RCHLA "cp novel.splice.scores.rc.key.combined.noHLA.vaf.highexp novel_splice_creating_variants.final"
	print RCHLA "cat ".$script_dir."/resource/header novel.splice.scores.rc.key.combined.noHLA.vaf.highexp > novel_splice_creating_variants.final";
#	print RCHLA "cat $run_dir/NS_CASE/misplice.input.maf.*.filtered.5 |awk -F\'\\t\' \'{print substr(\$16,1,12)\"_\"\$5\"_\"\$6\"_\"\$11\"_\"\$13\"\\t\"\$0}' > \${MAFK}\n";	
#	print RCHLA "     ".$run_script_path."combine.pl \${MAFK} \${HLAVAFH} 0 0 > \${FINAL}\n";
	
	close RCHLA;
	$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);
}



#mv case.control.distributionmethod* tmp_misplice/.
#mv novel.splice.scores.rc.key tmp_misplice/.
#mv novel.splice.scores.rc.key.combined tmp_misplice/.
#mv novel.splice.scores.rc.key.combined.noHLA tmp_misplice/.
#mv novel.splice.scores.rc.key.combined.noHLA.vaf tmp_misplice/.
#mv novel.splice.scores.rc.key.combined.noHLA.vaf.highexp tmp_misplice/.


