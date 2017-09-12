use strict;

my $usage =<<USAGE;
 Usage: $0 <file 1> <file 2> <col to match file 1> <col to match file 2>

	col to match file 1 and 2 should be 0 based.

USAGE
    die $usage unless @ARGV==4;

#INPUT 
my $file1=$ARGV[0];
my $colfile1=$ARGV[2];
my $file2=$ARGV[1];
my $colfile2=$ARGV[3];

#OPEN FILE 1 STORE DATA
my %file1data;
open(my $F1,'<',$file1) or die "Can't open $file1!";
	while(my $line=<$F1>){
		chomp $line;
		my @linearray=split(/\t/,$line);
		$file1data{$linearray[$colfile1]}=$line;
	} 
close $F1;

#OPEN FILE 2 STORE DATA
my %file2data;
open(my $F2,'<',$file2) or die "Can't open $file2!";
	while (my $line2=<$F2>){
        chomp $line2;
		my @linearray2=split(/\t/,$line2);
        $file2data{$linearray2[$colfile2]}=$line2;
    }
close $F2;

#PRINT OUT RESULTS
foreach my $sites (keys %file1data){
	if (exists $file1data{$sites}){
		if (exists $file2data{$sites}){
			print "$file1data{$sites}\t$file2data{$sites}\n";
		}
	}
}

	
