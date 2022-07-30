#!/usr/bin/perl

use strict;
use warnings;

MAIN: {
    my $indir = $ARGV[0] || usage();
    my $outdir = $ARGV[1] || usage();
    my $min  = $ARGV[2] || usage();
    my $count = 0;

    check_outdir($outdir);

    opendir DIR, $indir or die "cannot open $indir:$!";
    my @files = readdir DIR;
    foreach my $f (@files) {
        open IN, "$indir/$f" or die "cannot open $indir/$f:$!";
        my $count = 0;
        my $seqs = '';
        my %species = ();
        while (my $line = <IN>) {
            $seqs .= $line;
            next unless ($line =~ m/^>([^|]+)/);
            my $sp = $1;
            $count++ unless ($species{$sp});
            $species{$sp}++;
        }
        write_seqs($outdir,$f,$seqs) if ($count >= $min);
    }
}

sub check_outdir {
    my $outdir = shift;
    if (-d $outdir) {
        opendir OUTDIR, $outdir or die "cannot read $outdir:$!";
        my @existing = grep {!/^\.\.?$/} readdir OUTDIR;
        foreach my $e (@existing) {
            warn "warning: $outdir exists and includes $e\n";
        }
    } else {
        mkdir $outdir or die "cannot open $outdir";
    }
}

sub write_seqs {
    my $dir = shift;
    my $file = shift;
    my $seqs = shift;
    open OUT, ">$dir/$file" or die "cannot open >$dir/$file:$!";
    print OUT $seqs;
    close OUT; 
}

sub usage {
    die "usage: $0 INDIR OUTDIR MINIMUM_SEQS\n";    
}


