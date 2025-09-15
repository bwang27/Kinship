#!/usr/bin/perl -w

###this code is used to find the 2nd rank of kinship, to compare with the first rank.

system "awk '\$15==0 {print}' newfinal_pair.txt >diffcellline.lookfor2nd.txt";
system "awk '\$15==1 && \$18==1 {print \$1}' newfinal_pair.txt | awk '!seen[\$0]++' >corrected_ids";
system "awk '\$15==1 && \$18==1 {print}' newfinal_pair.txt >summary.txt";

open KINSHIP,"diffcellline.lookfor2nd.txt"||die $!;
my %chip_wgs;
$all=join'',<KINSHIP>;
@kin=split(/\n/,$all);

#open CHIP,"432corrected.chip.ids"||die $!;
open CHIP,"corrected_ids"||die $!;
open OUT,">corrected_samples_with_closest_kinship_in_diffcellline.txt"||die $!;
#open OUT,">523_with_closest_kinship_in_diffcellline.txt"||die $!;
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
#open FINAL,">final.523_with_closest_kinship_in_diffcellline.txt"||die $!;
open FINAL,">final.corrected_samples_with_closest_kinship_in_diffcellline.txt"||die $!;

$header="SampleA	SampleB	Ap_Bp	Ap_Bn	An_Bp	An_Bn	kinship	Category	bam.size.SampleA	bam.size.SampleB	batch.SampleA	Celllinename.SampleA	batch.SampleB	Celllinename.SampleB	same.cell.line.or.not	sampleA.No.	sampleB.No.	same.sampleid.or.not	Sample.of.K2	K2";

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

