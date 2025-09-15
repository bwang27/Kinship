
use Data::Dumper;

sub convert_to_matrix {
  my $ifn = shift;
  my $ofn = shift;
  print $ifn, "\t", $ofn, "\n";
  
  my $samples= {};
  open(my $ifh, "<$ifn");
  while (my $line = <$ifh>) {
    chomp $line;
    next if($line =~ /Ap_Bp/);
    #$line =~ s/.txt//g;
    my @eles = split /\t/, $line;
    $samples->{$eles[0]} = 1;
    $samples->{$eles[1]} = 1;
  }
  close($ifh);
  my $mat = {};
  open(my $ifh, "<$ifn");  
  while (my $line = <$ifh>) {
    chomp $line;
    my @eles = split /\t/, $line;
    my $diff = $eles[3] + $eles[4];
    my $NN = $eles[3] + $eles[4] + $eles[2] + $eles[5];
    my $ratio = 1000*$diff / (0.1+$NN);
    $ratio = int($ratio)/10;
    $mat->{$eles[0]}->{$eles[1]} = $ratio;
    $mat->{$eles[1]}->{$eles[0]} = $ratio;
  }
  close($ifh);    

  open(my $ofh, ">$ofn");
  print $ofh "Sample";
  foreach my $smp (keys %{$samples}) {
    print $ofh "\t",$smp;
  }
  print $ofh "\n";
  foreach my $smp (keys %{$samples}) {
    print $ofh $smp;
    foreach my $smp2 (keys %{$samples}) {
      print $ofh "\t", $mat->{$smp}->{$smp2};
    }
    print $ofh "\n";
  }  
  close($ofh);
}

my $ifn=$ARGV[0];
my $ofn=$ARGV[1];
convert_to_matrix($ifn,$ofn)


