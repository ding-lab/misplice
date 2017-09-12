#!/gsc/bin/perl
use strict;
#rjayasin@genome.wustl.edu
#Edited May 1st, 2016
#Changes: Updated usage, added comments, changed format of extracting controls, updated restricted sample list usage
### song changes the control samples to any samples without the mutation found in the case #
#### 5/4/2017 ####
my $usage =<<USAGE;
 Usage: $0 <All Mutations File> <Mutation List of Interest> <Output File> <Restricted Sample List>
	<All Mutations File> - A file containing all mutations in all genes, this list will be u
sed to determine samples that shouldn't be used as controls.
	<Mutation List of Interest> - A file containing all mutations you are interested in for
further analysis, control samples will be added to the last column of this file.
	<Output File> - Name of the output file generated with mutation info and control information.
	<Restricted Sample List> - A list of samples that you want to restrict your control list to.

#EXAMPLE:perl controls.v2..pl all_the_mutations.maf tcga_maf_of_interest.maf tcga_control_stringent.maf samples_list

USAGE
die $usage unless @ARGV==4;

my (@ss,@ss1,%rnasamples,$id,$cancer,$sample,$gene,$file,$linest,$line,$line1,%mutation);
my %sn_2_ct=();

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
			$sn_2_ct{$shortname}=$cancertype; 		
		}
close SPLICESAMPLES;
}
###############################
#DETERMINING GENES OF INTEREST#
###############################
my $counter=0;
my %casehash;
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
    #my $gene=$line[0];#gene (Hugo Symbol)
    #if(length($line[15])<15) { $line[15].="-01"; }
    if(length($line[10])<15) { $line[10].="-01"; }
    my $sn=substr($line[10],0,12);
	#my $sn=
    my $cancertype=$sn_2_ct{$sn};#cancertype
	#my $id=$cancertype.$gene;
	my $chr=$line[1]; 
	my $pos=$line[2];
	my $id=$chr."_".$pos;
	$casehash{$cancertype}{$id}=1;	
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
    $sample=$ss[10];
	my $shortsnmaf=$1 if $sample=~/(^TCGA-\w{2}-\w{4})/;
    $gene=$ss[0];
	if(defined $sn_2_ct{$shortsnmaf}) 
	{
		$cancer=$sn_2_ct{$shortsnmaf}; 
#Only store mutation info about gene of interest if gene is found in our <Mutation List of Interest>MAF
	 my $chr=$ss[1];
    my $pos=$ss[2];
    my $id=$chr."_".$pos;
	if (exists $casehash{$cancer}{$id}){
 	    	$mutation{$cancer}{$id}{$shortsnmaf}=1;
   		 }   
	}
}
close $OGMAF;

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
	my $samplemaf=$ss1[11];
#shorten sample name to only include unique identifiers TCGA-XX-XXXX
	my $shortsamplemaf=$1 if $samplemaf=~/(^TCGA-\w{2}-\w{4})/;
	chomp $shortsamplemaf;			
	my $cancertype=$sn_2_ct{$shortsamplemaf};#cancertype
    my $chr=$ss1[1];
    my $pos=$ss1[2];
    my $id=$chr."_".$pos;
	my $controlstring="";
	foreach my $s (sort keys %{$rnasamples{$cancertype}})
		{
		if (exists $mutation{$cancertype}{$id}{$s}){
		 next; 
			}
		else{
			if($controlstring eq "")
			{
			 $controlstring=$s;
			}	
			else { $controlstring.=",".$s; }
			}
		}
		if($controlstring ne "") 
		{
		  print $OUT $line,"\t",$controlstring,"\n";
		
		}
		else { print $OUT $line,"\t","No Controls\n";}
	}
close $MAF;
print "Finished, $OUT created\n";
