#!/usr/bin/perl

use strict;
use warnings;
use File::Copy;

MAIN: {
    my $fa_dir  = $ARGV[0] || usage();
    my $treedir = $ARGV[1] || usage();
    my $outdir  = $ARGV[2] || usage();
    my $min     = $ARGV[3] || usage();
    my $count   = 0;

    check_outdir($outdir);

    opendir DIR, $fa_dir or die "cannot open $fa_dir:$!";
    my @files = readdir DIR;
    foreach my $f (@files) {
        next if ($f =~ m/SpeciesTreeAlignment.fa/);
        open IN, "$fa_dir/$f" or die "cannot open $fa_dir/$f:$!";
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
        if ($count >= $min) {
            write_seqs($outdir,$f,$seqs);
            write_trees($treedir,$outdir,$f);
        }
    }
}

sub write_trees {
    my $tdir   = shift;
    my $outdir = shift;
    my $fasta  = shift;

    $fasta =~ m/^([^\/]+).fa$/ or die "unexpected format of fasta: $fasta";
    my $id = $1;
    File::Copy::copy("$tdir/${id}_tree.txt","$outdir/$id.tree");
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
    die "usage: $0 FASTA_DIR TREE_DIR OUTDIR MINIMUM_SEQS\n";    
}


