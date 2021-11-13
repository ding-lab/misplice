#!/usr/bin/perl

use strict;
use warnings;

(my $usage = <<OUT) =~ s/\t+//g;
perl jaf_filter.pl f_in f_out
OUT

die $usage unless @ARGV == 3;

my ($f_in, $f_out1,$f_out2) = @ARGV;
open(OUT,">$f_out1");
open(OUT2,">$f_out2");
foreach my $l (`cat $f_in`) 
{
my $ltr=$l; 
chomp($ltr); 
my @t=split("\t",$ltr); 
my $jaf=$t[20]/$t[13];
if($jaf>0.05) 
{
print OUT $t[0],"\t",$jaf,"\t",$ltr,"\n";
}
else { print OUT2 $t[0],"\t",$jaf,"\t",$ltr,"\n"; }

}
close OUT;
