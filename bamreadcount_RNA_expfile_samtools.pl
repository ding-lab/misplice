use strict;

my $usage =<<USAGE;
 Usage: $0 <splice.score.format>
#Sample Chr Pos Ref_snp Alt_snp Strand  Junc    Ref_seq Alt_seq Ref_score   Alt_score
#TCGA-05-4244-01    17  42169750    C   A   5ss -   42169763    CAGGTGCCG   CAGGTGCCG6.46   6.46
#TCGA-05-4249-01    2   175111467   C   A   5ss -   175111467   TTGGTATTG   TTGTTATTG1.45   -7.05

USAGE
    die $usage unless @ARGV==2;
my (@array1,$sample,$file,$linest);

#my ($MAF,$RNAbamlist)=@ARGV; 

open(my $MAF,'<',$ARGV[0]) or die "INPUT MAF not found!";

my %sample_hash;
my $counter=0;
#Some cancer types use chr prefix in the bam file and this changes whether or not "chr" is included in the sites file and it changes the rereference fasta used for determining readcount calculations
#This initial section will load in information for each sample derived from samtools view -H to determine whether or not "chr" needs to be used as a prefix.

open(my $RNAbamlist,'<',$ARGV[1]);
 
#TCGA-06-0649-01B-01R-1849-01 gbm /gscmnt/gc2741/ding/TCGA_RNA/GBM/TCGA-06-0649-01B-01R-1849-01.bam 1
#TCGA-41-3915-01A-01R-1850-01 gbm /gscmnt/gc2741/ding/TCGA_RNA/GBM/TCGA-41-3915-01A-01R-1850-01.bam chr1
#TCGA-06-0190-01A-01R-1849-01 gbm /gscmnt/gc2741/ding/TCGA_RNA/GBM/TCGA-06-0190-01A-01R-1849-01.bam 1

#open(my $NOBAM,'>',"NObamfile.log"); 
my %chrprefixhash;
my %bamhash;

while(my $bamline=<$RNAbamlist>){
	chomp $bamline;
	my @baminfo=split(/ /,$bamline);
	my $bsample=$baminfo[0];
	my $sn=substr($bsample,0,15);
#	print "$sn\n";
#Skip all normal tissue samples
	next if $sn=~/TCGA-\w{2}-\w{4}-10/;
	next if $sn=~/TCGA-\w{2}-\w{4}-11/;
	my $snshort=substr($bsample,0,12);
	my $c=$baminfo[1];
	my $bam=$baminfo[2];
	my $chrprefix=$baminfo[3];
	$chrprefixhash{$snshort}=$chrprefix;
	$bamhash{$sn}=$bam;
	$bamhash{$snshort}=$bam;
}
close $RNAbamlist;

#GOING THROUGH EACH LINE OF SCORING METHOD TO EVALUATE TOTAL READS AT NOVEL JUNCTION SITE
#Sample	Chr	Pos	Ref_snp	Alt_snp	Strand	Junc	Ref_seq	Alt_seq	Ref_score	Alt_score
#TCGA-05-4244-01	17	42169750	C	A	5ss	-	42169763	CAGGTGCCG	CAGGTGCCG6.46	6.46
#TCGA-05-4249-01	2	175111467	C	A	5ss	-	175111467	TTGGTATTG	TTGTTATTG1.45	-7.05

while(<$MAF>){
		chomp ($linest=$_);
		next if $linest=~/^Sample/;
		my @array1=split(/\t/,$linest);
		my $sample=$array1[0];
		my $chr=$array1[1];
		my $pos=$array1[7];	

		my $shortsn=substr($sample,0,12);
		my $position=$chr.":".$pos."-".$pos;
        #my $bam=$bamhash{$shortsn};
        my $bam=$bamhash{$sample};
		my $chrprefixvalue=$chrprefixhash{$shortsn};
		my $reads;
		#print "$bam\t$position\n";
        if ($chrprefixvalue=~/chr/){
			my $posfinal="chr".$position;
			$reads = `samtools view $bam $posfinal|cut -f1,2|sort|uniq|wc -l`;				
        }
        else {
        	$reads = `samtools view $bam $position|cut -f1,2|sort|uniq|wc -l`;
        }
		chomp $reads;
		print "$array1[0]\t$array1[1]\t$array1[2]\t$array1[3]\t$array1[4]\t$array1[5]\t$array1[6]\t$array1[7]\t$array1[8]\t$array1[9]\t$array1[10]\t$array1[11]\t$reads\n";

}
close $MAF;
