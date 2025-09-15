#!/usr/bin/perl -w
##usage: perl creat_job.pl CHIPSeq.lst;

$prefix="bsub -q standard -eo eo.txt -oo oo.txt -M8000 Rscript ";

$r="/research/groups/xmagrp/projects/SensitiveDetection/common/bwang/kinship/20250825/pipeline/tmp/identify_genetic_dup_1_v_all_low.r";

open CHIP,"$ARGV[0]"||die $!;
#open WGS,"$ARGV[1]"||die $!;
#$pdep=join'',<WGS>;
#@each_file=split(/\n/,$pdep);
#open OUT,">CHIP.jobs.sh"||die $!;
while(<CHIP>){
	chomp;
	@sample=split(/\//,$_);
	$sample_id=pop(@sample);
	($name=$sample_id)=~s/.txt//g;;
	#	foreach $file(@each_file){
	print  "$prefix\t$r\t$_\tCHIP_WGS.lst\tout.$name.txt\n";
}

close CHIP;
