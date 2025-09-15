#!/usr/bin/perl -w

#usage: perl faster_add_more_infotokinship.pl merge_v_WGS.txt final_pairs_4plot.txt PedDepAcc_deliver_v11.txt runtable_03152025.txt;
#usage: perl faster_add_more_infotokinship_chip.pl merge_v_CHIPWGS.txt final_pairs_4plot.txt wgs_chip.bam.size PedDepAcc_deliver_v11.txt;
# perl faster_add_more_infotokinship_chip.pl merge_v_CHIPWGS.txt final_pairs_4plot.txt wgs_chip.bam.sizewgs.chip.cellline.txt

$header="SampleA\tSampleB\tAp_Bp\tAp_Bn\tAn_Bp\tAn_Bn\tkinship\tCategory\n";

open APBPANBN,"$ARGV[0]"||die $!;
$all=join'',<APBPANBN>;
@carrier=split(/\n/,$all);
my %APBP;
foreach $element(@carrier){
	@new=split(/\t/,$element);
	$name=join("\t",$new[0],$new[1]);
	$value=join("\t",$new[2],$new[3],$new[4],$new[5]);
	$APBP{$name}=$value;
}

close APBPANBN;

open OUT,">step1_faster.out"||die "error write output\n";
print OUT $header;

open PAIR,"$ARGV[1]"||die $!;
while(<PAIR>){
	chomp;
	@aa=split(/\t/,$_);
	$aa[0]=~s/\./\-/g;
	$aa[1]=~s/\./\-/g;
	$pair=join("\t",$aa[0],$aa[1]);
	if (exists $APBP{$pair}){
		print OUT "$pair\t$APBP{$pair}\t$aa[2]";
	}
	if($aa[2]==0 && $aa[0] ne $aa[1]){
		print OUT "\tduplicates";}
	elsif($aa[2]==0 && $aa[0] eq $aa[1]){
		print OUT "\tsample.itself";}
	elsif($aa[2]>0 && $aa[2]<10){
		print OUT "\tdiagnosis\/relapse";}
	elsif($aa[2]>=10){
		print OUT "\tunique";}
	print OUT "\n";
}

close PAIR;
close OUT;

open BAM,"$ARGV[2]"||die $!;
my %bam_size;
while(<BAM>){
	chomp;
	@comp1=split(/\t/,$_);
	$bam_size{$comp1[1]}=$comp1[0];
}

$header_line="SampleA\tSampleB\tAp_Bp\tAp_Bn\tAn_Bp\tAn_Bn\tkinship\tCategory\tbam.size.SampleA\tbam.size.SampleB\n";
open FINAL,">step2.out"||die $!;
print FINAL $header_line;

open STEP1,"step1_faster.out"||die $!;
while(<STEP1>){
	chomp;
	if($_!~/^Sample/){
	print FINAL "$_";
	@sample_id=split(/\t/,$_);
	#	@sampleid1=split(/\_/,$sample_id[0]);
	#@sampleid2=split(/\_/,$sample_id[1]);
	if (exists $bam_size{$sample_id[0]}){
		print FINAL "\t$bam_size{$sample_id[0]}";}
	else{print FINAL "\tNA";}
	if (exists $bam_size{$sample_id[1]}){
                print FINAL "\t$bam_size{$sample_id[1]}";}
	else{print FINAL "\tNA";}

#print FINAL "\n";

}
	print FINAL "\n";
}

close FINAL;
close STEP1;

open PEDDEP,"$ARGV[3]"||die $!;
my %sample_celline;
while(<PEDDEP>){
	chomp;
	@comp=split(/\t/,$_);
	$batch_cellline=join("\t",$comp[1],$comp[2]);
	$sample=$comp[0];
	$sample_celline{$sample}=$batch_cellline;
}

$header_line="SampleA\tSampleB\tAp_Bp\tAp_Bn\tAn_Bp\tAn_Bn\tkinship\tCategory\tbam.size.SampleA\tbam.size.SampleB\tbatch.SampleA\tCelllinename.SampleA\tbatch.SampleB\tCelllinename.SampleB\n";

open FINAL3,">step3.out"||die $!;
print FINAL3 $header_line;

open STEP3,"step2.out"||die $!;
while(<STEP3>){
	chomp;
	if(($_!~/^Sample/)||($_!~/^epi/)){
	print FINAL3 "$_";
	@sample_id=split(/\t/,$_);
	#@sampleid1=split(/\_/,$sample_id[0]);
	#@sampleid2=split(/\_/,$sample_id[1]);
	

	if (exists $sample_celline{$sample_id[0]}){
		print FINAL3 "\t$sample_celline{$sample_id[0]}";}
	else{print FINAL3 "\tNA\tNA";}
	if (exists $sample_celline{$sample_id[1]}){
                print FINAL3 "\t$sample_celline{$sample_id[1]}";}
	else{print FINAL3 "\tNA\tNA";}

#print FINAL "\n";

}
	elsif($_=~/^epi/){
		print FINAL3 "$_";
		@sample_id=split(/\t/,$_);
		#@sampleid1=split(/\_/,$sample_id[0]);
		#@sampleid2=split(/\_/,$sample_id[1]);
		if (exists $sample_celline{$sample_id[0]}){
                	print FINAL3 "\t$sample_celline{$sample_id[0]}";}
		else{print FINAL3 "\tNA\tNA";}
		if (exists $sample_celline{$sample_id[1]}){
                print FINAL3 "\t$sample_celline{$sample_id[1]}";}
		else {print FINAL3 "\tNA\tNA";}
}
	print FINAL3 "\n";
}

close FINAL3;
close STEP1;

###select CHIP WGS pairs only, excluding CHIP:CHIP and wgs:wgs PAIRS;
#
open CHIP_WGS,"step3.out"||die $!;
open CHIP_WGS_OUT,">CHIP_WGS_pair.txt"||die $!;
while(<CHIP_WGS>){
	if($_=~/^Sample/){
		print CHIP_WGS_OUT "$_";}
	else{
		chomp;
		@chip_wgs_pair=split(/\t/,$_);
		$chip_wgs_pair[0]=~s/epig_//g;
		$chip_wgs_pair[1]=~s/epig_//g;
		@ele1=split(/\_/,$chip_wgs_pair[0]);
		@ele2=split(/\_/,$chip_wgs_pair[1]);
		if(@ele1>2 && @ele2==2){
			print CHIP_WGS_OUT "$_\n";}
		elsif(@ele1==2 && @ele2>2){
			print CHIP_WGS_OUT "$_\n";}
	}
	#	print "\n";
}

close CHIP_WGS;
close CHIP_WGS_OUT;

############################################################
##this part maynot be used;
#
open PEDDEPID,"PedDep_Data_IDs.txt"||die $!;
$ped=join'',<PEDDEPID>;
@pedchip=split(/\n/,$ped);
my %pedchip_cell;
foreach $exp(@pedchip){
	@exp_aa=split(/\t/,$exp);
	$pedchip_cell{$exp_aa[4]}=join("\t",$exp_aa[1],$exp_aa[2]);
}

open STEP3,"step3.out"||die $!;
while(<STEP3>){
	chomp;
	if($_=~/^Sample/){
		print "$_\t";}
	else{
		@step3_s=split(/\t/,$_);
		$prefix=join("\t",$step3_s[0],$step3_s[1],$step3_s[2],$step3_s[3],$step3_s[4],$step3_s[5],$step3_s[6],$step3_s[7],$step3_s[8],$step3_s[9]);
		print "$prefix\t";
		$s1=join("\t",$step3_s[10],$step3_s[11]);
		$s2=join("\t",$step3_s[12],$step3_s[13]);
		foreach $dataid (keys%pedchip_cell){
			if($step3_s[0]=~/$dataid/){
				print "$pedchip_cell{$dataid}\t";
			next;}
			else{print "$s1\t";}
			last;
		#	print "$s1\t";}
		}
		foreach $data(keys%pedchip_cell){
			if($step3_s[1]=~/$data/){
				print "$pedchip_cell{$data}\t";
			next;}
			else{print "$s2\t";}
			last;
			#	print "$s2\n";}
		}
	}
	print "\ n";
}
#}


close STEP3;
close PEDDEPID;

