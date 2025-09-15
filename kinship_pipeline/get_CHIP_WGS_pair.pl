#!/usr/bin/perl -w

open PAIR,"$ARGV[0]"||die $!;
while(<PAIR>){
	chomp;
	@st=split(/\t/,$_);
	@nd=split(/\:/,$st[0]);
	@nd1=split(/\_/,$nd[0]);
	@nd2=split(/\_/,$nd[1]);
	if(@nd1>2 && @nd2==2){
		$new=join("\t",$nd[0],$nd[1]);
		print "$new\t$st[1]\n";
	}elsif(@nd1==2 && @nd2>2){
		$new=join("\t",$nd[1],$nd[0]);
  		print "$new\t$st[1]\n";

}

}


close PAIR;
