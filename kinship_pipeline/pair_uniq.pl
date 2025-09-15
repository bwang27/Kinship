#!/usr/bin/perl -w

open A,"pair_value.txt"||die $!;
my %pairs;
while(<A>){
	chomp;
	my ($pair,$number)=split(/\t/,$_);
	my ($node1,$node2)=split(/\:/,$pair);
	my @nodes=sort($node1, $node2);
	my $canonical_pair="$nodes[0]:$nodes[1]";

	unless(exists $pairs{$canonical_pair}){
		$pairs{$canonical_pair}=$number;
	}
}

close A;

open OUT,">final_pairs_4plot.txt"||die "errror write file";
foreach my $pair(sort keys%pairs){
	print OUT "$pair\t$pairs{$pair}\n";

}

close OUT;

