#!/usr/bin/perl -w
## perl get_pair_value.pl thematrix_clustered.chip.txt 
##
open MATRIX,"$ARGV[0]"||die "error open file";
open PAIR,">pair_value.txt"||die "error write file\n";
my %sample_pair;
while (<MATRIX>){
	if($_=~/^Sample/){
		chomp;
		@header=split(/\t/,$_);
		$count=0;
	}
	else{
		chomp;
		#my %sample_pair;
		@row=split(/\t/,$_);
		$rowname=$row[0];
		shift(@row);
		foreach $number(@row){
			$count++;
			#$group1="$rowname:$header[$count]";
			#$group2="$header[$count]:$rowname";
			#$sample_pair{$group1}=$number;
			#$sample_pair{$group2}=$number;	
	print PAIR "$rowname:$header[$count]\t$number\n";}
		}

	$count='';
}
#}

	#if((exists $sample_pair{$group1}) || (exists $sample_pair{$group2})){
	#       print "$group1\t$sample_pair{$group1}\n";}



close MATRIX;



