#!/usr/bin/env perl
# cpg_islands.pl - detect CpG islands in a FASTA file (Takai & Jones 2002 criteria)
#   length >= 500 bp, GC content >= 55%, observed/expected CpG ratio >= 0.65
# Usage: perl cpg_islands.pl [--window 500] [--step 1] [--gc 0.55] [--oe 0.65] genome.fa > islands.bed
# Output: BED6-style lines (chrom, start, end, name, length, GC, obs/exp CpG)
use strict;
use warnings;
use Getopt::Long;

my ($window, $step, $min_gc, $min_oe) = (500, 1, 0.55, 0.65);
GetOptions('window=i' => \$window, 'step=i' => \$step,
           'gc=f' => \$min_gc, 'oe=f' => \$min_oe)
    or die "Usage: $0 [--window N] [--step N] [--gc F] [--oe F] in.fa\n";
my $fasta = shift @ARGV or die "Usage: $0 [options] in.fa\n";

# Read FASTA into a hash: id => sequence (uppercased)
my (%seqs, @order, $id);
open my $fh, '<', $fasta or die "Cannot open $fasta: $!\n";
while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    if ($line =~ /^>(\S+)/) { $id = $1; push @order, $id; $seqs{$id} = ''; next; }
    die "Sequence data before first header\n" unless defined $id;
    $seqs{$id} .= uc $line;
}
close $fh;

# Composition stats for one window
sub stats {
    my ($s) = @_;
    my $len = length $s;
    my $c   = ($s =~ tr/C//);
    my $g   = ($s =~ tr/G//);
    my $cg  = () = $s =~ /CG/g;
    my $gc  = ($c + $g) / $len;
    my $oe  = ($c && $g) ? ($cg * $len) / ($c * $g) : 0;
    return ($gc, $oe);
}

my $n = 0;
for my $chrom (@order) {
    my $seq = $seqs{$chrom};
    my $L   = length $seq;
    next if $L < $window;
    my ($start, $end);   # current merged island (0-based, half-open)
    for (my $i = 0; $i + $window <= $L; $i += $step) {
        my ($gc, $oe) = stats(substr($seq, $i, $window));
        if ($gc >= $min_gc && $oe >= $min_oe) {
            if (defined $end && $i <= $end) { $end = $i + $window; }   # extend
            else { report($chrom, $start, $end, $seq) if defined $start;
                   ($start, $end) = ($i, $i + $window); }
        }
    }
    report($chrom, $start, $end, $seq) if defined $start;
}

sub report {
    my ($chrom, $s, $e, $seq) = @_;
    my ($gc, $oe) = stats(substr($seq, $s, $e - $s));
    # merged regions are re-checked as a whole, as in Takai & Jones
    return unless $gc >= $min_gc && $oe >= $min_oe;
    $n++;
    printf "%s\t%d\t%d\tCGI_%d\t%d\t%.3f\t%.3f\n", $chrom, $s, $e, $n, $e - $s, $gc, $oe;
}
