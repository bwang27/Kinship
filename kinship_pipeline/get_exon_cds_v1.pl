#!/usr/bin/perl
use strict;
use warnings;

# Input GTF file
my $gtf_file = 'mm10.gtf'; 

# Data structures
my %genes;  # gene_id -> hash of gene info
my %transcripts; # transcript_id -> gene_id
my %exons;  # transcript_id -> list of exon regions
my %cds;    # transcript_id -> list of CDS regions
my %gene_regions; # gene_id -> list of gene regions (exons + gene span)
my %transcript_info; # transcript_id -> hash of info

# Read GTF file
open(my $fh, '<', $gtf_file) or die "Cannot open $gtf_file: $!";

while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^#/;  # Skip comments
    my @fields = split("\t", $line);
    next unless @fields >= 9;

    my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attributes) = @fields;

    # Parse attributes
    my %attr;
    foreach my $attr (split /;\s*/, $attributes) {
        next unless $attr =~ /(\S+)\s+"([^"]+)"/;
        $attr{$1} = $2;
    }

    # Only process features of interest
    if ($feature eq 'exon') {
        # Get transcript_id
        my $transcript_id = $attr{'transcript_id'} or next;
        # Get gene_id
        my $gene_id = $attr{'gene_id'} or next;

        # Store transcript to gene mapping
        $transcripts{$transcript_id} = $gene_id;

        # Initialize data structures
        $exons{$transcript_id} ||= [];
        push @{$exons{$transcript_id}}, [$start, $end];

        # Store gene info
        unless (exists $genes{$gene_id}) {
            $genes{$gene_id} = {
                gene_name => $attr{'gene_name'} // 'NA',
                chr => $chr,
                strand => $strand,
            };
        }
    }
    elsif ($feature eq 'CDS') {
        my $transcript_id = $attr{'transcript_id'} or next;
        # Store CDS regions
        $cds{$transcript_id} ||= [];
        push @{$cds{$transcript_id}}, [$start, $end];

        # Ensure transcript to gene mapping exists
        unless (exists $transcripts{$transcript_id}) {
            # Might be missing transcript info; skip
            next;
        }
    }
}

close($fh);

# Now, process each gene and its transcripts
print join("\t",
    "gene_name",
    "transcript_id",
    "RefSeq_ID",
    "gene_region",
    "exon_regions",
    "coding_regions",
    "num_exons",
    "num_coding_exons",
    "gene_length",
    "coding_length"
), "\n";

# For each gene
for my $gene_id (keys %genes) {
    my $gene_name = $genes{$gene_id}->{gene_name} // 'NA';
    my $chr = $genes{$gene_id}->{chr};
    my $strand = $genes{$gene_id}->{strand};

    # Find all transcripts for this gene
    my @transcripts_for_gene = grep { $transcripts{$_} eq $gene_id } keys %transcripts;

    # For each transcript
    for my $transcript_id (@transcripts_for_gene) {
        # Get transcript-specific info
        my $refseq_id = $transcript_id; # Assuming transcript_id as RefSeq ID; modify if needed
        my @exon_regions = exists $exons{$transcript_id} ? sort_regions($exons{$transcript_id}) : ();
        my @cds_regions = exists $cds{$transcript_id} ? sort_regions($cds{$transcript_id}) : ();

        # Calculate gene region (min start to max end of exons)
        my ($gene_start, $gene_end) = get_region_bounds(@exon_regions);
        my $gene_region_str = defined($gene_start) ? "$chr:$gene_start-$gene_end" : 'NA';

        # Exon regions string
        my $exon_regions_str = @exon_regions ? join(";", map { "$$_[0]-$$_[1]" } @exon_regions) : 'NA';

        # Coding regions string
        my $coding_regions_str = @cds_regions ? join(";", map { "$$_[0]-$$_[1]" } @cds_regions) : 'NA';

        # Number of exons
        my $num_exons = scalar @exon_regions;

        # Number of coding exons (CDS exons)
        my $num_coding_exons = scalar @cds_regions;

        # Calculate gene length
        my $gene_length = defined($gene_start) ? ($gene_end - $gene_start + 1) : 'NA';

        # Calculate coding length
        my $coding_length = 0;
        if (@cds_regions) {
            ($coding_length) = calculate_total_length(@cds_regions);
        } else {
            $coding_length = 'NA';
        }

        print join("\t",
            $gene_name,
            $transcript_id,
            $refseq_id,
            $gene_region_str,
            $exon_regions_str,
            $coding_regions_str,
            $num_exons,
            $num_coding_exons,
            $gene_length,
            $coding_length
        ), "\n";
    }
}

# Helper functions

# Sort regions by start
sub sort_regions {
    my ($regions_ref) = @_;
    my @sorted = sort { $a->[0] <=> $b->[0] } @$regions_ref;
    return @sorted;
}

# Get overall region bounds
sub get_region_bounds {
    my @regions = @_;
    return (undef, undef) unless @regions;
    my $min = $regions[0][0];
    my $max = $regions[0][1];
    for my $region (@regions) {
        $min = $region->[0] if $region->[0] < $min;
        $max = $region->[1] if $region->[1] > $max;
    }
    return ($min, $max);
}

# Calculate total length of overlapping regions
sub calculate_total_length {
    my @regions = @_;
    # Merge overlapping regions
    my @sorted = sort { $a->[0] <=> $b->[0] } @regions;
    my @merged;
    my ($start, $end) = @{$sorted[0]};
    for my $region (@sorted[1..$#sorted]) {
        if ($region->[0] <= $end + 1) {
            $end = $region->[1] if $region->[1] > $end;
        } else {
            push @merged, [$start, $end];
            ($start, $end) = @$region;
        }
    }
    push @merged, [$start, $end];

    # Sum lengths
    my $total = 0;
    for my $mr (@merged) {
        $total += ($mr->[1] - $mr->[0] + 1);
    }
    return ($total);
}
