#!/usr/bin/perl
use strict;use warnings;

my $CHROMOSOME_NR   = $ARGV[0];
my $CHR_SIZE = $ARGV[1];
my $FN_IN = $ARGV[2];
my $SEQ_PATH = $ARGV[3];

open(IN,"<",$FN_IN);
my $LINE = <IN>;
my @SEQ;
my $linenr=0;
while($LINE=<IN>){
	my ($s)=$LINE=~m/^([A-Za-z]*)\n$/;
	$s =~ s/A|a/1/g;
	$s =~ s/C|c/2/g;
	$s =~ s/T|t/3/g;
	$s =~ s/G|g/4/g;
	$s =~ s/N|n/0/g;
	push(@SEQ,split(//,$s));
	#$linenr++;
	#if($linenr==1000){last;}
	}

print(scalar(@SEQ),"\n");

my $count=0;
my @BIN_SEQ;
foreach my $s (@SEQ){	
	if ($count==0){
		#$bin_nr++;
		my $FN_OUT=join('',$SEQ_PATH,"/Chr",$CHROMOSOME_NR,"_BINALL",".perl.ascii");
		print("$FN_OUT\n");
		open(OUT,">",$FN_OUT);
		$count++;
		}
	elsif ($count<=$CHR_SIZE){	
			print OUT "$s\n";
			$count++;
			}
	else{close(OUT);$count=0;}
	}

