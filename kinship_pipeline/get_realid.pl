#!/usr/bin/perl -w
#
#usage: perl get_realid.pl namekey/namekey.txt tmpv2_ovlap_WGS_kinship_pair.txt 
#
open HE,"$ARGV[0]"||die $!;
$all=join'',<HE>;
@aa=split(/\n/,$all);
my %namekey;
foreach $ele(@aa){
	@name=split(/\t/,$ele);
	$val=join("\t",$name[1],$name[3]);
	$namekey{$name[0]}=$val;
}

open OO,">newfinal_pair.txt"||die $!;

$header="SampleA	SampleB	Ap_Bp	Ap_Bn	An_Bp	An_Bn	kinship	Category	bam.size.SampleA	bam.size.SampleB	batch.SampleA	Celllinename.SampleA	batch.SampleB	Celllinename.SampleB	same.cell.line.or.not	sampleA.No.	sampleB.No.	same.sampleid.or.not";

print OO  "$header\n";

open KINSHIP,"$ARGV[1]"||die $!;
while(<KINSHIP>){
	chomp;
	next if /^Sample/;
	@bb=split(/\t/,$_);
	if (exists $namekey{$bb[0]}){
		@new=split(/\t/,$namekey{$bb[0]});
		$bb[0]=$new[0];
		$bb[10]=$new[1];
	}

		if($bb[10] eq $bb[12]){
			$cell=1;
		}
		elsif($bb[10] ne $bb[12]){
		 	$cell=0;}
		@nus1=split(/\_/,$bb[0]);
		@nus2=split(/\_/,$bb[1]);
		$nus1[0]=~s/[^0-9]//g;
		$nus2[0]=~s/[^0-9]//g;
		if($nus1[0]==$nus2[0]){
			$sampleid=1;
		}
		elsif($nus1[0]!=$nus2[0]){
			$sampleid=0;
		}

#			print OUT_WGS "$line\t$aa[14]\t$nus1[0]\t$nus2[0]\t$aa[15]\n";}
	$replace=join("\t",@bb);
	print OO "$replace\t$cell\t$nus1[0]\t$nus2[0]\t$sampleid\n";
}

close HE;
close OO;
close KINSHIP;

	
