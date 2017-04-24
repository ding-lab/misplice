#!/gsc/bin/perl
use strict;
#rjayasin@genome.wustl.edu
#Edited May 1st, 2016
#Changes: Updated usage, added comments, changed format of extracting controls, updated restricted sample list usage
my $usage =<<USAGE;
 Usage: $0 <All Mutations File> <Mutation List of Interest> <Output File> <Restricted Sample List>
	<All Mutations File> - A file containing all mutations in all genes, this list will be u
sed to determine samples that shouldn't be used as controls.
	<Mutation List of Interest> - A file containing all mutations you are interested in for
further analysis, control samples will be added to the last column of this file.
	<Output File> - Name of the output file generated with mutation info and control information.
	<Restricted Sample List> - A list of samples that you want to restrict your control list to.

#EXAMPLE:perl controls_3.pl all_the_mutations.maf tcga_maf_of_interest.maf tcga_control_stringent.maf samples_list

USAGE
die $usage unless @ARGV==4;

my (@ss,@ss1,%rnasamples,$id,$cancer,$sample,$gene,$file,$linest,$line,$line1,%mutation);
open(my $OGMAF,'<',$ARGV[0]) or die "INPUT MAF not found!";
######################
#DEFINING SAMPLE LIST#
######################
####Read <Restricted Sample List> then read file
if($ARGV[3]=~/\w+/){
print "Reading <Restricted Sample List>\t$ARGV[3]\n";
	open(my $SPLICESAMPLES,'<',$ARGV[3]) or die "Restricted Sample List not found!";
		while(<$SPLICESAMPLES>){
   			chomp ($line1=$_);
			my @temp=split(/\t/,$line1);
			my $sn = $temp[0];
			my $cancertype=lc($temp[1]);
			my $shortname=$1 if $sn=~/(^TCGA-\w{2}-\w{4})/;
    		$rnasamples{$cancertype}{$shortname}=1;
		}
close SPLICESAMPLES;
}
###############################
#DETERMINING GENES OF INTEREST#
###############################
my $counter=0;
my %genehash;
#Extract all the genes names from the MAF of interest
print "\tFinding all genes of interest from <Mutation List of Interest>...\n";
open(my $MAF,'<',$ARGV[1]) or die "NEW MAF not found!";
while (<$MAF>){
#print header and add Controls to header line
    if ($_=~/^Hugo_Symbol/){
        next;
    }   
    chomp ($line=$_);
    my @line=split(/\t/,$line);
    my $gene=$line[0];#gene (Hugo Symbol)
    my $cancertype=$line[90];#cancertype
	my $id=$cancertype.$gene;
	$genehash{$cancertype}{$gene}=();	
}
close $MAF;
############################################
#DETERMING SAMPLES and GENES with MUTATIONS#
############################################
#Go through file that contains mutations for all genes and store sample name, cancer type, and gene name into hash if the sample contains a mutation. 
print "Defining samples and genes with mutations\n";
while(<$OGMAF>){
    next if /Hugo_Symbol/;
    $counter++;
    chomp ($linest=$_);
    @ss=split(/\t/,$linest);
    $sample=$ss[15];
	my $shortsnmaf=$1 if $sample=~/(^TCGA-\w{2}-\w{4})/;
    $gene=$ss[0];
    $cancer=$ss[90];
#Only store mutation info about gene of interest if gene is found in our <Mutation List of Interest>MAF
	if (exists $genehash{$cancer}{$gene}){
#Currently list of mutation types does not include Silent Mutations
	    if ($ss[8]=~/Translation_Start_Site|Targeted_Region|IGR|Splice_Site|Nonsense_Mutation|Frame_Shift_Del|Frame_Shift_Ins|3'Flank|3'UTR|5'Flank|5'UTR|In_Frame_Del|In_Frame_Ins|Intron|Missense_Mutation|Nonstop_Mutation|RNA/){
   	    	$mutation{$cancer}{$gene}{$shortsnmaf}=1;
   		 }   
#Exceptions would include any other mutation type than those specified in line 49 
### Song Marked out
	    #else {
    	 #   print "Not considering gene,mutation type:$gene,$ss[8]\n";
       		#$mutation{$cancer}{$gene}{$sample}=0;
    	#}	
	}
	else{next;}   
}
close $OGMAF;
#####################################################
#DEFINING CONTROL LIST FOR EACH GENE AND CANCER TYPE#
#####################################################
print "Defining control list for the following cancer types and genes:\n";
my (@controls,$short,$counter);
my %controls;
####For each gene in the genehash derived from the <Mutation List of Interest>, find out which samples from the <Restricted Sample List> have mutations in the gene of interest, and create a control list for each gene.
foreach my $c (sort keys %genehash){
	foreach my $g (sort keys %{$genehash{$c}}){		
		my $id=$c.$g;
	#	print "$id\n";
#Go through all samples for cancer type of interest
		foreach my $s (sort keys %{$rnasamples{$c}}){
#If the sample has a mutation found in  <All Mutations File> then this sample can't be used as a control
			if (exists $mutation{$c}{$g}{$s}){
				next;
			}
#Sample does not have a mutation in the gene of interest, add to control array for cancer type and gene
			else{
				if (exists $controls{$id}){
					push (@{$controls{$id}},$s);
				}
				else{
					my @controls;
					push @controls,$s;
					$controls{$id}=[@controls];
				}	
			}
		}
	}
}
#########################
##ADDING CONTROLS TO MAF#
#########################
#Finally go through mutation file containing mutations of interest and add all samples that do not contain mutations in gene of interest, for that same cancer type.
print "Going through MAF to add controls to....\n";
open(my $MAF,'<',$ARGV[1]) or die "NEW MAF not found!";
open(my $OUT,'>',$ARGV[2]) or die "OUTPUT not created!";
while (<$MAF>){
#print header and add Controls to header line
	if ($_=~/^Hugo_Symbol/){
		chomp $_;
		print $OUT "$_\tControls\n";
		next;
	}
	$counter++;
	chomp ($line=$_);
	@ss1=split(/\t/,$line);
#Extract sample name from MAF file
	my $samplemaf=$ss1[15];
#shorten sample name to only include unique identifiers TCGA-XX-XXXX
	my $shortsamplemaf=$1 if $samplemaf=~/(^TCGA-\w{2}-\w{4})/;
	chomp $shortsamplemaf;
	my $cancertype=$ss1[90];#cancertype
	my $genemaf=$ss1[0];#gene (Hugo Symbol)	
#Extract controls for cancer type and gene of interest
	my $idmaf=$cancertype.$genemaf;
	if (exists $rnasamples{$cancertype}{$shortsamplemaf}){
		if (exists $controls{$idmaf}){
			my @controlarray=@{$controls{$idmaf}};
			#Remove current sample from control array #just a double check
			my @controlarrayfinal = grep !/$shortsamplemaf/, @controlarray;
			print "$shortsamplemaf\t@controlarrayfinal\n";
			my $controlstring=join(',',@controlarrayfinal);
			print $OUT "$line\t$controlstring\n";	
		}
		else{
			print $OUT "$line\tNo Controls\n";
		}
	}
}
close $MAF;
print "Finished, $OUT created\n";
