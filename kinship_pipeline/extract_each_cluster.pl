#!/usr/bin/perl -w

open MARKER,"horizontal_somatic_mutations.txt"||die $!;
$all=join'',<MARKER>;
@marker=split(/\n/,$all);
$header=shift(@marker);
close MARKER;
#print "$header\n";

open NAMEMAP,"namemap.txt"||die $!;
$NAME=join'',<NAMEMAP>;
@namemap=split(/\n/,$NAME);
close NAMEMAP;

@number=(S1,S2,S3,S4,S5);

open CLU,"93cluster"||die $!;
my @arr;

while(<CLU>){
	$HEAD=$header;
	chomp;
	push(@arr,$_);
	#print "$header\n";
	foreach $line(@namemap){
		@aa=split(/\t/,$line);
		@bb=split(/\//,$aa[1]);
		$sample=pop(@bb);
		$sample=~s/.bam//g;
		if($_ eq $aa[2]){
		push(@arr,"$aa[3]:$sample");
		#print "$aa[3]\:$sample\t";
	}
}

	#push(@arr,"\n");
	print "@arr\n";
	shift@arr;
	foreach $ss(@arr){
		#print "$ss\t";
		@dd=split(/\:/,$ss);
		$search_s=$dd[0];
		$replace_s=$dd[1];
		$HEAD=~s/$search_s/$replace_s/g;
	}
	@ww=split(/\t/,$HEAD);
	shift@ww;
	$rest=join("\t",@ww);
	print "cluster_id\tmarker_NO.\t$rest\n";
	$count=0;
	foreach $ele(@marker){
		@cc=split(/\t/,$ele);
		if($_ eq $cc[0]){
			$count++;
			shift@cc;
			$new=join("\t",@cc);
			print "$_\tmarker$count\t$new\n";}
	}
	@arr=();
	#$header;
}
#}
