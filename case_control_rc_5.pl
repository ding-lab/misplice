use strict;
use warnings;

my $usage =<<USAGE;
 Usage: $0 <prefix> <largest_file_#>
  Example: perl case_control_rc_5.pl prefix suffix tot_number

USAGE
    die $usage unless @ARGV==3;

my $total_number=$ARGV[2];
my $prefix=$ARGV[0];
my $suffix=$ARGV[1];
for(1...$total_number){

#novel.junctions.allgenes.filtered.5.controls.1.num2.run2
	my $control_file=$prefix.".".$_.".$suffix";
	my $control_maf=$prefix.".".$_;
#	print "Processing:$_\t$control_file\n";

####CONTROL FILE FORMAT (ENDS WITH .num)####
#TCGA-05-4384	1
#TCGA-05-4396	1
#TCGA-05-4397	1
#TCGA-05-4422	1

	open(my $CON,'<',$control_file) or die "Can't open $control_file!";

	#Counter to keep track of the number of controls and number of reads supporting this alternative junction for this particular site and gene
	my $sum=0;
	my $controlnumber=0;
	my @control_reads;
	my @control_samples;

	#Go through each line of the control file concatenate all junction supporting reads and control info into an array

	while(my $cline=<$CON>){
		chomp $cline;
		my ($sample,$reads)=split(/\t/,$cline);	
		$sum=$reads+$sum;
		$controlnumber++;
		
		push @control_samples,$sample;
		push @control_reads,$reads;
		#print "$sum(+$reads)\n";	
	} 

	close $CON;
	
	my $sizecontrol=scalar @control_reads;
	
	if ($sizecontrol==0){
		print STDERR "$control_maf\n";
	}

	my $readarray = join ',', @control_reads;
	my $samplearray = join ',', @control_samples;
	

#	print "$readarray\n$samplearray\n";

	open(my $MAF,'<', $control_maf) or die "Can't open $control_maf!";
	my ($casereads,$cancer,$sampleid,$genename,$Variant);
	while(my $mline=<$MAF>){
		chomp $mline;
		my @mafline=split(/\t/,$mline);
		$casereads=$mafline[103];
		$cancer=$mafline[90];
		$sampleid=$mafline[15];
		$genename=$mafline[0];
		$Variant=$mafline[4]."_".$mafline[5]."_".$mafline[10]."_".$mafline[12];
	
	}
	close $MAF;
	my $RSEM=$cancer."\t".$sampleid."\t".$genename."\t".$Variant."\t".$casereads."\t".$readarray;

	my ($rsem_data)=parse_gene_rsem($RSEM);

	my ($pvalue)=gene_expression($rsem_data);

	print "$pvalue\t$RSEM\n";
	
	#print "$casereads\n$readarray\nSUM:$sum,#CON:$controlnumber\n";

#Format for running gene expression test
#Cancer_Type	Sample_ID	Gene_Name	Variant	RSEM (case)	RSEM (control)
#blca	TCGA-E5-A2PC-01A-11D-A202-08	STAG2	X_123182853_A_C	2392.9020	2284.3111,2521.0872
}

######SUBROUTINES########
#written by Mike Wendl, Edited by Reyka Jayasinghe
#Adapted from Module1.pm from SpliceInator

sub parse_gene_rsem {
   #my ($file) = @_;
   my $rsem_data = {};
   my $data=$_[0];
#__OPEN FILE AND READ
   #open (F, $file) || die "cant open file $file";
   #while (<F>) {

   #__FILTERING
      next if /^#/ || /^Cancer\_Type/;
      next if /^\s+/;
      chomp;

   #__PARSE
      my ($cancer, $sample, $gene, $variant, $rsem_case, $rsem_control) = split /\t/,$data;
      my @control_values = split /,/, $rsem_control;

	
   #__STORE
      $rsem_data->{$cancer}->{$gene}->{$sample} = {
         'case' => $rsem_case,
         'controls' => [@control_values],
      };
   #}
   #close (F);

#__RETURN DATA
   return $rsem_data;
}



#  =======================
#  GENE EXPRESSION (RSEM)
#  =======================

sub gene_expression {
   my ($rsem_data) = @_;

#__INITIALIZE
   my $diagnostics = {};
   my ($group_i_count, $group_i_hits, $mode) = (0, {}, "gene expression");
#   my ($value) = (0, {}, "gene expression");

#__DATA STRUCTURES
   my ($control_check, $data) = ({}, {});
	my $value;
#__ITERATE OVER ALL WHOLE-GENE CASE-CONTROLS TO PROPERLY GROUP THE DATA
   foreach my $cancer (keys %{$rsem_data}){
      foreach my $gene (keys %{$rsem_data->{$cancer}}){
         foreach my $sample (keys %{$rsem_data->{$cancer}->{$gene}}){

         #__EXTRACT WHOLE GENE CONTROL RSEM AND CASE RSEM VALUES
            my $local_hash = $rsem_data->{$cancer}->{$gene}->{$sample};
            my $control_vals = $local_hash->{'controls'}; # list reference
            my $case_val = $local_hash->{'case'};         # scalar value

         #__DIAGNOSTICS
            $diagnostics->{$cancer}->{$gene}->{$sample}->{$mode} = $case_val;

         #__NUMBER OF ELEMENTS IN DISTRIBUTION
            my $num_vals = scalar @{$control_vals};
         #__NEED AT LEAST 3 ELEMENTS
            $group_i_hits->{$cancer}->{$gene}->{$sample} =
                  "distribution_too_small" unless $num_vals >= 3;

         #__ORDER THE LIST FROM LOWEST TO HIGHEST
            my $list;
            @{$list} = sort _numerical_ @{$control_vals};

			#Create hash to store #reads supporing controls
			my %control_reads;
			#Go through list of controls and sort them into bins
			#based on number of reads supporting novel junction
			foreach (@{$list}){
				my $supporting_reads=$_;
			
			#Create hash that as a key stores the value of the reads supporting the novel junction as the value stores the total number of samples with that number of supporting reads
			#EXAMPLE:
			#0 ==> 19
			#1 ==> 22
			#2 ==> 13
			#3 ==> 16
				if (exists $control_reads{$supporting_reads}){
					my $number_samples=$control_reads{$supporting_reads};
					my $new_number_samples=$number_samples+1;
					$control_reads{$supporting_reads}=$new_number_samples;
				}
				else {
					$control_reads{$supporting_reads}=1;
				}
			} 


#Number of elements in control distribution.
			my $total = scalar @{$control_vals};	
			my $tally=0;
			my $cutoff;
#FIND CUTOFF FOR TOP 5%
			foreach my $key (sort _numerical_ keys %control_reads){
				$tally += $control_reads{$key};
				my $fractionofsamples=$tally/$total;
				#print "$control_reads:$tally\t";
				#if the number of samples up to this point is greater than 95% of samples then save this value as cutoff and breakout of loop
				if ( $fractionofsamples > 0.95 ){
					$cutoff = $key;
					#print "reads:$key\tsamples:$control_reads{$key}\t$fractionofsamples\tCUTOFF\n";
					last;
				}	
			}

#FIND CUTOFF FOR BOTTOM 5%
			my $lowercutoff;
			my $tallylower;
		   foreach my $key (sort _numerical_ keys %control_reads){
                $tallylower += $control_reads{$key};
                my $fractionofsamples=$tallylower/$total;
                #print "$control_reads:$tally\t";
                #if the number of samples up to this point is greater than 95% of samples then save this value as cutoff and breakout of loop
                if ( $fractionofsamples < 0.05 ){
                    $lowercutoff = $key;
                    #print "reads:$key\tsamples:$control_reads{$key}\t$fractionofsamples\tCUTOFF\n";
                    last;
                }
            }

#If lower cutoff is not set this could indicate that many reads at some particular value in the control disribution such that the fraction of samples for a particular read value is never less than 0.05	
			#	print STDERR "$lowercutoff\t$cutoff\t$case_val\n";
			#if (($lowercutoff=~/^$/g) && ($cutoff=~/\d*/g)){
			if ( ( !defined $lowercutoff) && (defined $cutoff) ){
			#This means that 95% of samples have a 0 readcount value - This means it is likely not expressed in this particular tissue type	
				if( ($case_val == 0) && ($cutoff == 0) ){
					$value="not_expressed_in_tissue";
				}
			#If the case value is 0 and less than the high cutoff it is considered not expressed
				elsif ( ($case_val < $cutoff) && ($case_val == 0) ){
					$value="not_expressed";
				}
			#If the case value is greater than the set high cutoff - then it is still more expressed than 95% of samples in the contorl distribution
				elsif (( $case_val >= $cutoff)){
					$value="high_expression";
				}
			#If the case val is less than the high cutoff - then place it in the average expression bin since no lower cutoff can be justified
				elsif (( $case_val < $cutoff)){
					$value="avg_expression";
				}
			#This should never occur
				else{
					die;
					$value="PEWPEWPEWERROR:".$case_val."_".$lowercutoff."_".$cutoff;
				}
			}
##############
#Cutoff Check#
##############
        #__INFERENCE / CLASSIFICATION
		#If both the lower and higher cutoff are set
			elsif ( (defined $cutoff) && (defined $lowercutoff)){
				if ( ($case_val >= $cutoff) && ($case_val > $lowercutoff)){
    	                $value="high_expression";
        	    }
				elsif(($case_val <= $lowercutoff) && ($case_val == 0) && ($case_val < $cutoff) && ($lowercutoff != $cutoff)){
					$value="not_expressed";

				}
				elsif( ($case_val <= $lowercutoff) && ($case_val != 0) && ($case_val < $cutoff) && ($lowercutoff != $cutoff)){
						$value="low_expression";
				} 
				elsif (($case_val < $cutoff) && ($case_val > $lowercutoff) && ($lowercutoff != $cutoff)) {
						$value="avg_expression";
				}
				else{
					die;
					$value="MEOWMEOWMEOWERROR:".$case_val."_".$lowercutoff."_".$cutoff;
				}
			}

			else{
				die;
				$value="PEWMEOWERROR:".$case_val."_".$lowercutoff."_".$cutoff;
			}
########
#die "END OF TEST";
########
         }
      }
   }

         ####################
         #  CLASSIFICATION  #
         ####################

   #return ($group_i_hits, $group_i_count, $diagnostics);
   return ($value);
}

####UTILITY####
sub _numerical_ {$a <=> $b}

