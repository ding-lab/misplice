##Split bed file to a ###

### Song Cao ####

#!/usr/bin/perl

use strict;
use warnings;
(my $usage = <<OUT) =~ s/\t+//g;
perl split_maf.pl file_in dir_out max_num
OUT

die $usage unless @ARGV == 3;
#max_num is the number of splited files ##
my ($file_in, $dir_out, $max_num) = @ARGV;

#my $f_bed=$dir_in."/output/all_splice.bed.sort";

#my $f_out=$dir_out."/output/all_splice.bed.sort1";

my $file_name= (split(/\//,$file_in))[-1];

my $f_out=$dir_out."/".$file_name.".1";


my $n=`wc -l $file_in`;
my $n1;
 
if($n=~/^(\d+)/) { $n1=$1; }

if($n1<$max_num) { `cp $file_in $f_out`; exit; }

my $max_var = int($n1/$max_num); 
my $cc=0;
 
foreach my $l (`cat $file_in`) 
	{

		if($l=~/^Hugo/) { next; }

		my $ind=int($cc/$max_var);
 
		if($cc-$ind*$max_var==0) 
		{

			if($cc!=0 && $ind<$max_num) 
			{ close OUT; }

			if($ind<$max_num)
			{ 
			#print $ind,"\n";
			my $ind2=$ind+1; 
			my	$f_out=$dir_out."/".$file_name.".".$ind2; 
			open(OUT,">$f_out");
			}
		}
		print OUT $l; 
		$cc++
	}

