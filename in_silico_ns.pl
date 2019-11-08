
### Song Cao ####

### get supporting reads for novel junctions due to somatic mutation ### 

### Oct 25, 2016###

## last updated: Jan 16, 2017 ###

#!/usr/bin/perl

use strict;
use warnings;
(my $usage = <<OUT) =~ s/\t+//g;
perl filter_fp_ns.pl f_in f_bed f_out
OUT

die $usage unless @ARGV == 3;
my ($f_in, $f_bed, $f_out) = @ARGV;

#my $f_bam_list="/gscmnt/gc2521/dinglab/scao/mmy/scm/worklog/rna_bam_path.v1.tsv"; 
my @tempdir=split(/\//,$f_in);

my $f_bam_list=$tempdir[0];

for(my $i=1;$i<scalar @tempdir-2;$i++)
{
    $f_bam_list=$f_bam_list."/".$tempdir[$i];
}

$f_bam_list=$f_bam_list."/rna_bam_path.v1.tsv";
print $f_bam_list,"\n";

#my $f_bed="/gscuser/scao/gc2524/dinglab/bed_maker/E75_bed_v3.sort.tsv";

##0-based##

my %bampath=();
my %bampathchr=();
my %known_junc=();
my %known_junc_s=(); 
open(OUT,">$f_out");

my $f_out2=$f_out.".detailed.alignment.5"; 
open(OUT2,">$f_out2"); 

foreach my $l (`cat $f_bed`) 
	{
		my $ltr=$l; 
		chomp($ltr); 
		my @temp=split("\t",$ltr); 
		my @temp2=split(":",$temp[3]); 		

		if($temp2[3]=~/^i/) { 
		#print $temp[0],"\t",$temp[1],"\t",$temp[2],"\n"; 
		#<STDIN>;
		$known_junc{$temp[0]}{$temp[1]}{$temp[2]}=1;
		$known_junc_s{$temp[0]}{$temp[1]}=1;
		$known_junc_s{$temp[0]}{$temp[2]}=1; 

			}

	}

foreach my $l (`cat $f_bam_list`) 
	{
		my $ltr=$l; 
		chomp($ltr); 
	 	my @temp=split(" ",$ltr); 
		#print $temp[2],"\n";
		#print $temp[4],"\n";	
		my $sn=$temp[0]; 
		#print $sn,"\n";
		#<STDIN>;
		if($temp[4]=~/chr/) { $bampathchr{$sn}=$temp[2]; }
		else { $bampath{$sn}=$temp[2]; }
	}

foreach my $l (`cat $f_in`)
	{
		my $ltr=$l; chomp($ltr); 
		my @temp=split("\t",$ltr);
	#	if(length($temp[15])<15) { $temp[12].="-01"; }
		my $sn=$temp[15];
		$sn=~s/\_T//g;  
		#print $sn,"\n";
		my $chr=$temp[4]; 
		my $pos=$temp[5];
		my $ref=$temp[10];
		my $var=$temp[12];
		my $dellen=0; 
		if($var eq "-") { $dellen=length($ref); } 
		#print $chr,"\t",$pos,"\n";
		my $chr_pos; 
		my $bam="NULL";
		my $chr_bk=$chr; 
		my $pos_bk=$pos; 
		if(defined $bampathchr{$sn} && (-e $bampathchr{$sn}))
		{
			my $left_pos=$pos-20; 
			my $right_pos=$pos+20; 
			$chr_pos="chr".$chr.":".$left_pos."-".$right_pos;
			$bam=$bampathchr{$sn};
		}
 
		if(defined $bampath{$sn} && (-e $bampath{$sn}))
        {
            my $left_pos=$pos-20;
            my $right_pos=$pos+20;
            $chr_pos=$chr.":".$left_pos."-".$right_pos;
			$bam=$bampath{$sn};
        } 
			
		if($bam ne "NULL" && (-e $bam)) 
		{
		# print $bam,"\n";
		#	my $com=`samtools view $bam \"$chr_pos\"`;
		#	print $bam,"\n";
		#	print $chr_pos,"\n"; 
			#<STDIN>;
		#	my @temp=split("\n",$com); 
			my %count_read=();
			foreach my $t (`samtools view $bam \"$chr_pos\"`)
			{
				chomp($t);
				my @temp2=split("\t",$t); 
				if($temp2[5]=~/^(\d+)M(\d+)N(\d+)M$/)
				{

					#print $temp2[5],"\n";
					#<STDIN>;
					my $chr=$temp2[2]; 
					my $start_pos=$temp2[3];
					my $id=$temp2[0]; 
					my $flag=$temp2[1];
					my $r1;
					my $r2; 
					my $rid=$id; 
					my $x=$start_pos+$1-1;
                    my $y=$start_pos+$1+$2-2;
					my $mapq=$temp2[4];
					if($id=~/\/2$/ || ($flag & 0x80)) { $rid=$id; $rid=~s/\/2$//g; $rid.="\/2";  }
					if($id=~/\/1$/ || ($flag & 0x40)) { $rid=$id; $rid=~s/\/1$//g; $rid.="\/1"; } 
					#my $x=$start_pos+$1-1; 
					#my $y=$start_pos+$1+$2-2; 
					#print $t,"\n";	
					#print $rid,"\n";
					#<STDIN>;
					#print $start_pos,"\n";
					#print $1,"\n";
					#print $2,"\n";
					#print $t,"\n";
					#print $x,"\t",$y,"\t",$mapq,"\n";
					#<STDIN>;			
					if($chr=~/^chr/) { $chr=$chr; }
                    else { $chr="chr".$chr; }
					#print $chr,"\t",$x,"\t",$y,"\t",$temp2[6],"\t",$mapq,"\n";
					#<STDIN>;
					if(!defined $known_junc{$chr}{$x}{$y} && (!($y-$x+1==$dellen)) && ($temp2[6] eq "=") && (($x>=$pos-20 && $x<=$pos+20 && (!defined $known_junc_s{$chr}{$x})) || ($y>=$pos-20 && $y<=$pos+20 && (!defined $known_junc_s{$chr}{$y}))))			 						
					{ 
				      	#print $t,"\n";  
                    	#print $rid,"\n";
                    	#print $start_pos,"\n";
                    	#print $1,"\n";
                    	#print $2,"\n"; 
						#print $chr,"\n";
						#print $x,"\t",$y,"\n";
						#<STDIN>;	
						$count_read{$rid}=$t; 
					}	
				}		
			}
		
			my $n_key=keys %count_read;

			if($n_key>=5) 
			{
			 my $ltr=$l;
			 chomp($ltr); 
			 print OUT $ltr,"\t",$n_key,"\n";  
			 print OUT2 $sn,"\t",$chr_bk,"\t",$pos_bk,"\t",$ref,"\t",$var,"\n";
			 foreach my $rr (sort keys %count_read)
				{
					print OUT2 $count_read{$rr},"\n";
				}
			}
		}	
		
	}
close OUT;
close OUT2; 
