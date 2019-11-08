#!/usr/bin/perl

use strict;
use warnings;

(my $usage = <<OUT) =~ s/\t+//g;
perl splice_score_novel.pl dir_s dir_in info f_ref f_bed
OUT

die $usage unless @ARGV == 5;

my ($dir_s,$dir_in, $info,$f_ref,$f_bed) = @ARGV;

print "Reading genome.......";
my $genome = $f_ref;

#print $genome,"\n";
#my ($dir_in) = @ARGV;

open(REF,"$genome");

my %gen;
my $id;
while(<REF>)
{
        chomp;
        my $first = substr($_,0,1);
        if($first eq ">")
        {
                my @line = split(/\s+/, substr($_,1));
                $id=$line[0];
		## remove chr if it has ##
		$id=~s/chr//g; 
        }else
        {
                $gen{$id}.=uc($_);
        }
}
print "Done!\n";

print "Reading annotation...";
my $bed = $f_bed;

my %strand;
open(IN, "$bed");
while(<IN>)
{
	chomp;
	next if $_ =~m/^#/;
	my @l = split(/\t/,);
#	next if $l[2] ne "gene";
#	my @m = split(/\;/, $l[8]);
#	for(my $i=0;$i<scalar @m;$i++) 
#	{
#	 if($m[$i]=~/gene_name/)
#	{
#	my @m2=split(/\"/,$m[$i]);	
	$strand{uc($l[0])} = $l[2];
#	}
#	}
}
print "Done!\n";

my $DIR = "$dir_s/CB_splice/";

##### calculate 5’ splice strength #####

my $tmpfile1=$DIR."me2x5";
my %me2x5 = &makescorematrix($tmpfile1);
my $tmpfile2=$DIR."splicemodels/splice5sequences";
my %seq = &makesequencematrix($tmpfile2);

my %bgd;
$bgd{'A'} = 0.27;
$bgd{'C'} = 0.23;
$bgd{'G'} = 0.23;
$bgd{'T'} = 0.27; 
  
sub makesequencematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) { 
	chomp;
	$_=~ s/\s//;
	$matrix{$_} = $n;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}

sub makescorematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) { 
	chomp;
	$_=~ s/\s//;
	$matrix{$n} = $_;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}

sub getrest5{
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    return $seqa[0].$seqa[1].$seqa[2].$seqa[5].$seqa[6].$seqa[7].$seqa[8];
}

sub scoreconsensus5{
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    my %bgd; 
    $bgd{'A'} = 0.27; $bgd{'C'} = 0.23; $bgd{'G'} = 0.23; $bgd{'T'} = 0.27;  
    my %cons1;
    $cons1{'A'} = 0.004; $cons1{'C'} = 0.0032; $cons1{'G'} = 0.9896; $cons1{'T'} = 0.0032;
    my %cons2;
    $cons2{'A'} = 0.0034; $cons2{'C'} = 0.0039; $cons2{'G'} = 0.0042; $cons2{'T'} = 0.9884;
    my $addscore = $cons1{$seqa[3]}*$cons2{$seqa[4]}/($bgd{$seqa[3]}*$bgd{$seqa[4]}); 
    return $addscore;
}

sub log2{
    my ($val) = @_;
    return log($val)/log(2);
}
##### calculate 5’ splice strength #####

##### calculate 3’ splice strength #####
sub hashseq{
    #returns hash of sequence in base 4
    # &hashseq('CAGAAGT') returns 4619
    my $seq = shift;
    $seq = uc($seq);
    $seq =~ tr/ACGT/0123/;
    my @seqa = split(//,$seq);
    my $sum = 0;
    my $len = length($seq);
    my @four = (1,4,16,64,256,1024,4096,16384);
    my $i=0;
    while ($i<$len) {
        $sum+= $seqa[$i] * $four[$len - $i -1] ;
	$i++;
    }
    return $sum;
}

my @metables = &makemaxentscores;
sub makemaxentscores{
    my $dir = $DIR."splicemodels/";
    my @list = ('me2x3acc1','me2x3acc2','me2x3acc3','me2x3acc4','me2x3acc5','me2x3acc6','me2x3acc7','me2x3acc8','me2x3acc9');
    my @metables;
    my $num = 0 ;
    foreach my $file (@list) {
	my $n = 0;
	open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
	while(<SCOREF>) {
	    chomp;
	    $_=~ s/\s//;
	    $metables[$num]{$n} = $_;
	    $n++;
	}
	close(SCOREF);
	#print STDERR $file."\t".$num."\t".$n."\n";
	$num++;
    }
    return @metables;
}
sub makewmmscores{
    my $dir = $DIR."splicemodels/";
    my @list = ('me1s0acc1','me1s0acc2','me1s0acc3','me1s0acc4','me1s0acc5','me1s0acc6','me1s0acc7','me1s0acc8','me1s0acc9');
    my @metables;
    my $num = 0 ;
    foreach my $file (@list) {
	my $n = 0;
	open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
	while(<SCOREF>) {
	    chomp;
	    $_=~ s/\s//;
	    $metables[$num]{$n} = $_;
	    $n++;
	}
	close(SCOREF);
	#print STDERR $file."\t".$num."\t".$n."\n";
	$num++;
    }
    return @metables;
}
sub makemmscores{
    my $dir = $DIR."splicemodels/";
    my @list = ('me2s0acc1','me2s0acc2','me2s0acc3','me2s0acc4','me2s0acc5','me2s0acc6','me2s0acc7','me2s0acc8','me2s0acc9');
    my @metables;
    my $num = 0 ;
    foreach my $file (@list) {
	my $n = 0;
	open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
	while(<SCOREF>) {
	    chomp;
	    $_=~ s/\s//;
	    $metables[$num]{$n} = $_;
	    $n++;
	}
	close(SCOREF);
	#print STDERR $file."\t".$num."\t".$n."\n";
	$num++;
    }
    return @metables;
}
sub maxentscore{
    my $seq = shift;
    my $table_ref = shift;
    my @metables = @$table_ref;
    my @sc;
    $sc[0] = $metables[0]{&hashseq(substr($seq,0,7))};
    $sc[1] = $metables[1]{&hashseq(substr($seq,7,7))};
    $sc[2] = $metables[2]{&hashseq(substr($seq,14,7))};
    $sc[3] = $metables[3]{&hashseq(substr($seq,4,7))};
    $sc[4] = $metables[4]{&hashseq(substr($seq,11,7))};
    $sc[5] = $metables[5]{&hashseq(substr($seq,4,3))};
    $sc[6] = $metables[6]{&hashseq(substr($seq,7,4))};
    $sc[7] = $metables[7]{&hashseq(substr($seq,11,3))};
    $sc[8] = $metables[8]{&hashseq(substr($seq,14,4))};
    my $finalscore = $sc[0] * $sc[1] * $sc[2] * $sc[3] * $sc[4] / ($sc[5] * $sc[6] * $sc[7] * $sc[8]);
    return $finalscore;
}    
    
sub getrest3{
    my $seq = shift;
    my $seq_noconsensus = substr($seq,0,18).substr($seq,20,3);
    return $seq_noconsensus;
}

sub scoreconsensus3{
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    my %bgd; 
    $bgd{'A'} = 0.27; $bgd{'C'} = 0.23; $bgd{'G'} = 0.23; $bgd{'T'} = 0.27;  
    my %cons1;
    $cons1{'A'} = 0.9903; $cons1{'C'} = 0.0032; $cons1{'G'} = 0.0034; $cons1{'T'} = 0.0030;
    my %cons2;
    $cons2{'A'} = 0.0027; $cons2{'C'} = 0.0037; $cons2{'G'} = 0.9905; $cons2{'T'} = 0.0030;
    my $addscore = $cons1{$seqa[18]} * $cons2{$seqa[19]}/ ($bgd{$seqa[18]} * $bgd{$seqa[19]}); 
    return $addscore;
}
##### calculate 3’ splice strength #####

sub revcom
{
        my $in = shift;
        $in=~tr/[AGCT]/[TCGA]/;
        return reverse($in);
}

open(OUT, ">$dir_in/novel.splice.scores");
print OUT "Sample\tChr\tPos\tRef_snp\tAlt_snp\tStrand\tJunc\tRef_seq\tAlt_seq\tRef_score\tAlt_score\n";

my $indel_skipped = 0;
print "Processing intputs.......";
my $dir = $dir_in;
#foreach my $input (glob("$dir/spliceinator.splice.score.*.v2.filtered.5"))
foreach my $input (glob("$dir/$info"))
{
	my %hash;
	open(IN, "$input");
	while(<IN>)
	{
		chomp;
		my @line = split(/\t/,);
		my $id = join("\t", $line[4], $line[5]);
		next if $line[0] eq "Unknown";
		$hash{$id} = $strand{uc($line[0])};
	}

	open(DATA, "$input.detailed.alignment.5.max");
	#print "while","\n";
	#print $input,"\n";
	while(<DATA>)
	{
		#my $line=$_; 
		#chomp($line);
		chomp;
		my @templ=split(/\t/,);
		#print $line,"\n";
		my $nl=scalar @templ; 
		#print $nl,"\n";
		#<STDIN>;
		if($nl==5)
		{
			my @l = split(/\t/,);
			#### using length 20 as a cut-off ####
			#print $l[3],"\t",$l[4],"\n";
	 		if(length($l[3])>=20 || length($l[4])>=20) { next; }
			my $nxt = <DATA>;
			my @m = split(/\t/, $nxt);
			my @number=($m[5]=~m/(\d+)\w/g);
			my $d1 = abs($m[3] + $number[0]- $l[2]);
			my $d2 = abs($m[3] + $number[0] + $number[1] - $l[2]);
			my $jd = join("\t", $l[1], $l[2]);
			next if !exists $hash{$jd};
			my ($type, $start, $ref, $mut, $junc);
			$l[3]=~s/-//;
			$l[4]=~s/-//;
			my $dif=length($l[3])-length($l[4]);
			if($hash{$jd} eq "+")
			{
				if($d1 > $d2)
				{
					$type = "3ss";
					$start = $m[3] + $number[0] + $number[1] -21;
					$junc = $m[3] + $number[0] + $number[1] -1;
					$ref = substr($gen{$l[1]}, $start, 23);
					if($l[2]>$start && $l[2]<=$start+23)
					{
						if($l[2]<$m[3] + $number[0] + $number[1]) #intron
						{
							if($dif==0)
							{#tested
								$mut = substr($gen{$l[1]}, $start-$dif, $l[2]-$start-1+$dif).$l[4].substr($gen{$l[1]}, $l[2]+length($l[4])-1, $start+24-$l[2]-length($l[4]));
							}elsif($dif>0)
							{
								if($start+24-$l[2]-$dif<=0)
								{
									$mut="AAAAAAAAAAAAAAAAAAAAAAA";
								}else
								{
									$mut = substr($gen{$l[1]}, $start-$dif, $l[2]-$start-1+$dif).$l[4].substr($gen{$l[1]}, $l[2]+$dif-1, $start+24-$l[2]-$dif);
								}
							}else
							{#tested
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start).$l[4].substr($gen{$l[1]}, $l[2], $start+23-$l[2]);
								$mut = substr($mut, length($mut)-23);
							}
						}else
						{#rare
							if($dif==0)
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start-1).$l[4].substr($gen{$l[1]}, $l[2]+length($l[4])-1, $start+24-$l[2]-length($l[4]));
							}elsif($dif>0)
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start-1).$l[4].substr($gen{$l[1]}, $l[2]+$dif-1, $start+24-$l[2]);
							}else
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start).$l[4].substr($gen{$l[1]}, $l[2], $start+23-$l[2]);
								$mut = substr($mut, 0, 23);
							}
						}
						
					}else
					{
						$mut = $ref;
					}
				}else	
				{
					$type = "5ss";
					$start = $m[3] + $number[0] - 4;
					$junc = $m[3] + $number[0];
					$ref = substr($gen{$l[1]}, $start, 9);
					if($l[2]>$start && $l[2]<=$start+9)		
					{
						if($l[2]<$m[3]+$number[0])
						{#rare
							if($dif==0)
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start-1).$l[4].substr($gen{$l[1]}, $l[2]+length($l[4])-1, $start+10-$l[2]-length($l[4]));
							}elsif($dif>0)
							{
								$mut = substr($gen{$l[1]}, $start-$dif, $l[2]-$start-1+$dif).$l[4].substr($gen{$l[1]}, $l[2]+$dif-1, $start+10-$l[2]);
							}else
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start).$l[4].substr($gen{$l[1]}, $l[2], $start+9-$l[2]);
								$mut = substr($mut, length($mut)-9);
							}
						}else
						{
							if($dif==0)
							{#tested
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start-1).$l[4].substr($gen{$l[1]}, $l[2]+length($l[4])-1, $start+10-$l[2]-length($l[4]));
							}elsif($dif>0)
							{#tested
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start-1).$l[4].substr($gen{$l[1]}, $l[2]+$dif-1, $start+10-$l[2]);
							}else
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start).$l[4].substr($gen{$l[1]}, $l[2], $start+9-$l[2]);
								$mut = substr($mut, 0, 9);
							}
						}
					}else
					{
						$mut = $ref;
					}
				}
			}else
			{
				if($d1 > $d2)
				{
					$type = "5ss";
					$junc = $m[3] + $number[0] + $number[1] -1;
					$start = $m[3] + $number[0] + $number[1] - 7;
					$ref = substr($gen{$l[1]}, $start, 9);	
					if($l[2]+length($l[3])-1>$start && $l[2]<=$start+9)
					{
						if($l[2]<$m[3] + $number[0] + $number[1])
						{
							if($dif==0)
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start-1).$l[4].substr($gen{$l[1]}, $l[2]+length($l[4])-1, $start+10-$l[2]-length($l[4]));
							}elsif($dif>0)
							{#tested
								$mut = substr($gen{$l[1]}, $start-$dif, $l[2]-$start-1+$dif).$l[4].substr($gen{$l[1]}, $l[2]+$dif-1, $start+10-$l[2]-$dif);
							}else
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start).$l[4].substr($gen{$l[1]}, $l[2], $start+9-$l[2]);
								$mut = substr($mut, length($mut)-9);
							}
						}else
						{
							if($dif==0)
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start-1).$l[4].substr($gen{$l[1]}, $l[2]+length($l[4])-1, $start+10-$l[2]-length($l[4]));
							}elsif($dif>0)
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start-1).$l[4].substr($gen{$l[1]}, $l[2]+$dif-1, $start+10-$l[2]);
							}else
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start).$l[4].substr($gen{$l[1]}, $l[2], $start+9-$l[2]);
								$mut = substr($mut, 0, 9);
							}
						}
					}else
					{
						$mut = $ref;
					}
				}else
				{
					$type = "3ss";
					$junc = $m[3] + $number[0];
					$start = $m[3] + $number[0] - 4;
					$ref = substr($gen{$l[1]}, $start, 23);
					if($l[2]>$start && $l[2]<=$start+23)
                                        {
						if($l[2]<$m[3]+$number[0])
						{
							if($dif==0)
							{
								$mut = substr($gen{$l[1]}, $start-$dif, $l[2]-$start-1+$dif).$l[4].substr($gen{$l[1]}, $l[2]+length($l[4])-1, $start+24-$l[2]-length($l[4]));
							}elsif($dif>0)
							{
								$mut = substr($gen{$l[1]}, $start-$dif, $l[2]-$start-1+$dif).$l[4].substr($gen{$l[1]}, $l[2]+$dif-1, $start+24-$l[2]-$dif);
							}else
							{#tested
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start).$l[4].substr($gen{$l[1]}, $l[2], $start+23-$l[2]);
								$mut = substr($mut, length($mut)-23);
							}
						}else
						{
							if($dif==0)
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start-1).$l[4].substr($gen{$l[1]}, $l[2]+length($l[4])-1, $start+24-$l[2]+$dif-length($l[4]));
							}elsif($dif>0)
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start-1).$l[4].substr($gen{$l[1]}, $l[2]+$dif-1, $start+24-$l[2]);
							}else
							{
								$mut = substr($gen{$l[1]}, $start, $l[2]-$start).$l[4].substr($gen{$l[1]}, $l[2], $start+23-$l[2]);
								$mut = substr($mut, 0, 23);
							}
						}
                                        }else
                                        {
                                                $mut = $ref;
                                        }
				}
				$mut = revcom($mut);
				$ref = revcom($ref);
			}
			my ($mut_score, $ref_score);
			if($type eq "5ss")
			{
				$ref_score = sprintf("%.2f",&log2(&scoreconsensus5($ref)*$me2x5{$seq{&getrest5($ref)}}));
				$mut_score = sprintf("%.2f",&log2(&scoreconsensus5($mut)*$me2x5{$seq{&getrest5($mut)}}));
			}else
			{
				$ref_score = sprintf("%.2f", &log2(&scoreconsensus3($ref)*&maxentscore(&getrest3($ref),\@metables)));
				$mut_score = sprintf("%.2f", &log2(&scoreconsensus3($mut)*&maxentscore(&getrest3($mut),\@metables)));
			}
			print OUT "$_\t$type\t$hash{$jd}\t$junc\t$ref\t$mut\t$ref_score\t$mut_score\n";
		}
	}
}
print "Done!\n";


