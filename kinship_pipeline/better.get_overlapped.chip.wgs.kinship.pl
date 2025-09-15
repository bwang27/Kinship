#!/usr/bin/perl -w
##perl get_overlapped.chip.wgs.kinship.pl new.CHIP_WGS_pair_afterstep4.txt overlapped.wgs.ids;
#$ARGV[0] here is 'CHIP_WGS_pair.txt' from faster_add_more_infotokinship_chip.pl step;
#
open KINSHIP,"$ARGV[0]"||die $!;
$kinship=join'',<KINSHIP>;
@pair=split(/\n/,$kinship);
my %chip_wgs;
my %wgs_chip;

$header="SampleA	SampleB	Ap_Bp	Ap_Bn	An_Bp	An_Bn	kinship	Category	bam.size.SampleA	bam.size.SampleB	batch.SampleA	Celllinename.SampleA	batch.SampleB	Celllinename.SampleB	same.cell.line.or.not	sampleA.No.	sampleB.No.	same.sampleid.or.not";

open WGS,"$ARGV[1]"||die $!;
open OUT_WGS,">v2_ovlap_WGS_kinship_pair.txt"||die $!;

print OUT_WGS "$header\n";

while(<WGS>){
	chomp;
	foreach $line(@pair){
		@aa=split(/\t/,$line);
		if(($_ eq $aa[0])||($_ eq $aa[1])){
			if($aa[10] eq $aa[13]){
				$aa[14]=1;
			}
			elsif($aa[10] ne $aa[13]){
			 	$aa[14]=0;}
			@nus1=split(/\_/,$aa[0]);
			@nus2=split(/\_/,$aa[1]);
			$nus1[0]=~s/[^0-9]//g;
			$nus2[0]=~s/[^0-9]//g;
			if($nus1[0]==$nus2[0]){
				$aa[15]=1;
			}
			elsif($nus1[0]!=$nus2[0]){
				$aa[15]=0;}
		
			print OUT_WGS "$line\t$aa[14]\t$nus1[0]\t$nus2[0]\t$aa[15]\n";}
	}
}

close WGS;
close KINSHIP;

