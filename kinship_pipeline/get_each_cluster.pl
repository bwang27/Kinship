#!/usr/bin/perl -w

open MAT,"cluster.txt"||die $!;
$count_pos=0;
$count_kinship=0;
$number_of_cluster=0;
my $ele;
my $id;
open OUT,">cluster.sample.val.txt"||die $!;
while(<MAT>){
	chomp;
	if($_=~/^Sample/){
		$head=$_;
		@header=split(/\t/,$head);
		foreach $ele(@header){
			$count_pos++;
		}
	}
	else{
		$count_kinship=0;
		$cluster_kin=0;
		my @cluster;
		my @sample;
		@each=split(/\t/,$_);
		$id=shift(@each);
		foreach $value(@each){
			$count_kinship++;
			if($value<10){
			$cluster_kin++;
			push(@cluster,$value);
			push(@sample,$header[$count_kinship]);
		}
		}
		if(@cluster>1){
			$number_of_cluster++;
			$sample_of_this_cluster=join(":",@sample);
			$value_of_kinship_in_this_cluster=join(":",@cluster);
		print OUT "Cluster:$number_of_cluster\t$sample_of_this_cluster\t$value_of_kinship_in_this_cluster\n";
		}
		#}
	undef @sample;
	undef @cluster;
}

}
#print "\n";

#}
close MAT;
close OUT;

system "awk -F '\t' '{print \$2}' cluster.sample.val.txt >tmp";
system "awk '!seen[\$0]++' tmp >samples.in.each.cluster";

open AA,"samples.in.each.cluster"||die $!;
$count=0;
$all=join'',<AA>;
@aa=split(/\n/,$all);
foreach $pair(@aa){
	$count++;
	print "cluster$count\t$pair\n";}
close AA;


