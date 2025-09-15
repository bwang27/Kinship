#!/usr/bin/perl -w

## usage: perl get_pair_value.pl thematrix_clustered.chip.txt 

## usage: perl merged.pl thematrix_clustered.txt merge_v_CHIPWGS.txt final_pairs_4plot.txt /research/groups/xmagrp/projects/SensitiveDetection/common/bwang/kinship/20250825/feature/wgs_chip.bam.size /research/groups/xmagrp/projects/SensitiveDetection/common/bwang/kinship/20250825/feature/wgs.chip.cellline.txt /research/groups/xmagrp/projects/SensitiveDetection/common/bwang/kinship/20250825/feature/namekey.txt
#

open MATRIX,"$ARGV[0]"||die "error open file\n";
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
		print PAIR "$rowname:$header[$count]\t$number\n";}
		}

	$count='';
}


close MATRIX;

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

open OUT,">final_pairs_4plot.txt"||die "error write file";
foreach my $pairr(sort keys%pairs){
	print OUT "$pairr\t$pairs{$pairr}\n";

}

close OUT;

#perl get_CHIP_WGS_pair.pl final_pairs_4plot.txt >tmp;
#open PAIR,"$ARGV[0]"||die $!;
open PAIR,"final_pairs_4plot.txt"||die "error open file";
open TMP,">tmp"||die $!;
while(<PAIR>){
	chomp;
	@st=split(/\t/,$_);
	@nd=split(/\:/,$st[0]);
	@nd1=split(/\_/,$nd[0]);
	@nd2=split(/\_/,$nd[1]);
	if(@nd1>2 && @nd2==2){
		$new=join("\t",$nd[0],$nd[1]);
		print TMP "$new\t$st[1]\n";
	}elsif(@nd1==2 && @nd2>2){
		$new=join("\t",$nd[1],$nd[0]);
  		print TMP  "$new\t$st[1]\n";

}
}

close PAIR;

system "mv tmp final_pairs_4plot.txt";

########### this part is to add bam size, cell line information to each kinship chip and wgs pairs#########
#
$header="SampleA\tSampleB\tAp_Bp\tAp_Bn\tAn_Bp\tAn_Bn\tkinship\tCategory\n";

open APBPANBN,"$ARGV[1]"||die $!;
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

open PAIR,"$ARGV[2]"||die $!;
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

open BAM,"$ARGV[3]"||die $!;
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

}
	print FINAL "\n";
}

close FINAL;
close STEP1;

open PEDDEP,"$ARGV[4]"||die $!;
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
	if(($_!~/^Sample/)||($_!~/^epi/)||($_!~/^\s*$/)){
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

############### select CHIP WGS pairs only, excluding CHIP:CHIP and wgs:wgs PAIRS ############
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

open KINSHIP,"CHIP_WGS_pair.txt"||die $!;
$kinship=join'',<KINSHIP>;
@pair=split(/\n/,$kinship);
my %chip_wgs;
my %wgs_chip;

$header="SampleA	SampleB	Ap_Bp	Ap_Bn	An_Bp	An_Bn	kinship	Category	bam.size.SampleA	bam.size.SampleB	batch.SampleA	Celllinename.SampleA	batch.SampleB	Celllinename.SampleB	same.cell.line.or.not	sampleA.No.	sampleB.No.	same.sampleid.or.not";

################ here needs to get the overlapped wgs samples ids first ###################
#
open WGS,"overlapped.wgs.ids"||die $!;
open OUT_WGS,">v2_ovlap_WGS_kinship_pair.txt"||die $!;

print OUT_WGS "$header\n";

while(<WGS>){
	chomp;
	foreach $line(@pair){
		@aa=split(/\t/,$line);
		if(($_ eq $aa[0])||($_ eq $aa[1])){
			if($aa[10] eq $aa[12]){
				$aa[14]=1;
			}
			elsif($aa[10] ne $aa[12]){
			 	$aa[14]=0;}
			@nus1=split(/\_/,$aa[0]);
			@nus2=split(/\_/,$aa[1]);
			$nus1[0]=~s/[^0-9]//g;
			$nus2[0]=~s/[^0-9]//g;
			if(($nus1[0]=~/\d+/)&&($nus2[0]=~/\d+/)&&($nus1[0]==$nus2[0])){
				$aa[15]=1;
			}
			#elsif($nus1[0]!=$nus2[0]){
			#	$aa[15]=0;}
			else {$aa[15]=0;}
			print OUT_WGS "$line\t$aa[14]\t$nus1[0]\t$nus2[0]\t$aa[15]\n";}
	}
}

close WGS;
close KINSHIP;

open HE,"$ARGV[5]"||die $!;
$all=join'',<HE>;
@aa=split(/\n/,$all);
my %namekey;
foreach $ele(@aa){
	@name=split(/\t/,$ele);
	$val=join("\t",$name[1],$name[3]);
	$namekey{$name[0]}=$val;
}

open OO,">newfinal_pair.txt"||die $!;

$header="SampleA	SampleB	Ap_Bp	Ap_Bn	An_Bp	An_Bn	kinship	Category	bam.size.SampleA	bam.size.SampleB	batch.SampleA	Celllinename.SampleA	batch.SampleB	Celllinename.SampleB	same.cell.line.or.not	sampleA.No.	sampleB.No.	same.sampleid.or.not	final.same.cell.line.or.not	sampleA.No.	sampleB.No.	final.same.sampleid.or.not";

print OO  "$header\n";

open KINSHIP,"v2_ovlap_WGS_kinship_pair.txt"||die $!;
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

#########this code is used to find the 2nd rank of kinship, to compare with the first rank #############

system "awk '\$19==0 {print}' newfinal_pair.txt >diffcellline.lookfor2nd.txt";
system "awk '\$19==1 && \$22==1 && \$7<20 {print \$1}' newfinal_pair.txt | awk '!seen[\$0]++' >corrected_ids";
system "awk '\$19==1 && \$22==1 && \$7<20 {print}' newfinal_pair.txt >summary.txt";

open KINSHIP,"diffcellline.lookfor2nd.txt"||die $!;
#my %chipp_wgs;
$all=join'',<KINSHIP>;
@kin=split(/\n/,$all);

open CHIP,"corrected_ids"||die $!;
open OUT,">corrected_samples_with_closest_kinship_in_diffcellline.txt"||die $!;
my @kinship_value;
my @kinship_sample;
while(<CHIP>){
	chomp;
	print OUT "$_\t";
	foreach $line(@kin){
		@aa=split(/\t/,$line);
		if($aa[0] eq $_){
			push @kinship_value,$aa[6];
			push @kinship_sample,$aa[1];
		}
	}

	$min=$kinship_value[0];
	foreach $num(@kinship_value){
		if($num<$min){
			$min=$num;}
	}

	for (my $i = 0; $i < @kinship_value; ++$i) {
 		   if ($kinship_value[$i] == $min) {
			print OUT "$kinship_sample[$i]\t";
    			last;
		}
    #	last;
}

		print OUT "$min\n";
	@kinship_value=();
	@kinship_sample=();

}

close CHIP;
close OUT;
close KINSHIP;

open CORRECT,"corrected_samples_with_closest_kinship_in_diffcellline.txt"||die $!;
$all=join'',<CORRECT>;
@nd=split(/\n/,$all);

open SUMMARY,"summary.txt"||die $!;
open FINAL,">final.corrected_samples_with_closest_kinship_in_diffcellline.txt"||die $!;

$header="SampleA	SampleB	Ap_Bp	Ap_Bn	An_Bp	An_Bn	kinship	Category	bam.size.SampleA	bam.size.SampleB	batch.SampleA	Celllinename.SampleA	batch.SampleB	Celllinename.SampleB	same.cell.line.or.not	sampleA.No.	sampleB.No.	same.sampleid.or.not	final.same.cell.line.or.not	sampleA.No.	sampleB.No.	final.same.sampleid.or.not	Sample.of.K2	K2";

print FINAL "$header\n";

while(<SUMMARY>){
	chomp;
	print FINAL "$_\t";
	@qry=split(/\t/,$_);
	foreach $ele(@nd){
		@qq=split(/\t/,$ele);
		if($qq[0] eq $qry[0]){
			print FINAL "$qq[1]\t$qq[2]";
		}
	}
	print FINAL "\n";
}

close SUMMARY;
close CORRECT;
close FINAL;


