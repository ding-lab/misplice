##! /gsc/bin/perl
################################################################
## parse_gtf.pl - Takes in a gtf and creates a bed file       ##
##      AUTHOR:vyellapa                                       ##
##      Last updated:Dec 10 2014                              ##
################################################################

use strict;
use warnings;
use FileHandle;
#use List::MoreUtils qw/ uniq /;

#my $fasta="/gscmnt/gc2534/dinglab/resources/hs37d5.fa";
if(!(-e $ARGV[0]))
{
        die "USAGE: perl $0 GTF\n";
}

open(GTF, "<$ARGV[0]") or die "Can't open sample file: $!\n";
my %gtf = ();
my %out = ();
my (@temp,@enst, @exon, @trans, @ensg, @tname, @gname);
my ($line, $site, $val);


warn "Reading $ARGV[0]...\n";
#Load lookup table into a hash
while(<GTF>) {
 chomp;
 $line=$_;
#Skip header and super contigs, patches
 next if(($line=~/^#/ || $line =~/^HG/ || $line=~/^HS/ || $line=~/^GL/ || $line=~/^MT/));
 
 @temp=split('\t', $line);
#Use exons only
 next if($temp[2] ne "exon");

#Parse required elements
 @trans=split(';',$temp[8]);
 @ensg=split(/"/,$trans[0]);
 @enst=split(/"/,$trans[2]);
 @exon=split(/"/,$trans[4]);
 @gname=split(/"/,$trans[5]);
 @tname=split(/"/,$trans[6]);

#key is gene_name transcript_name strand
my $key= join("\t",$gname[1],$enst[1],$temp[6]);

#value is exon1; chr:start-stop
#print $trans[2],"\t",$exon[1],"\n"; <STDIN>;
my $val="e" . $exon[1] . ";" . $temp[0] . ":" . $temp[3] . "-" . $temp[4];

#populate hash
unless ($out{$key}) {$out{$key}=$val;next;}

#If value is present, push new exon value into hash
if($out{$key}) {
my $vold=$out{$key};
my @t =();
push (@t, ($vold, $val));
$val=join("\t", @t); 
$out{$key}=$val;next;
}

#print "$gname[1]\t$tname[1]\t$temp[6]\te$exon[1];$temp[0]:$temp[3]:$temp[4]\t\n";
#print "$gname[1]\t$tname[1]\t$temp[6]\t$pos\n";
}

foreach my $s1 (keys(%out)) {
print "$s1\t$out{$s1}\n";
}

close(GTF);

