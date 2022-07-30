#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Copy;

our $VERSION = 0.01;

MAIN: {
    my $rh_o    = get_options();

    my $fa_dir  = $rh_o->{'fa_dir'};
    my $treedir = $rh_o->{'tree_dir'};
    my $outdir  = $rh_o->{'out_dir'};
    my $min     = $rh_o->{'min_taxa'};

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
    die "usage: $0 --fa_dir=FASTA_DIRECTORY --tree_dir=TREE_DIRECTORY --out_dir=OUT_DIRECTORY --min_taxa=MINIMUM_TAXA\n";
}

sub get_options {
    my $rh_opts = {};
    my $opt_results = Getopt::Long::GetOptions(
                                      'version' => \$rh_opts->{'version'},
                                     'fa_dir=s' => \$rh_opts->{'fa_dir'},
                                   'tree_dir=s' => \$rh_opts->{'tree_dir'},
                                    'out_dir=s' => \$rh_opts->{'out_dir'},
                                   'min_taxa=s' => \$rh_opts->{'min_taxa'});

    die "$0 version $VERSION\n" if ($rh_opts->{'version'});
    print "missing --fa_dir\n" unless ($rh_opts->{'fa_dir'});
    print "missing --tree_dir\n" unless ($rh_opts->{'tree_dir'});
    print "missing --out_dir\n" unless ($rh_opts->{'out_dir'});
    print "missing --min_taxa\n" unless ($rh_opts->{'min_taxa'});
    usage() unless ($rh_opts->{'fa_dir'});
    usage() unless ($rh_opts->{'tree_dir'});
    usage() unless ($rh_opts->{'out_dir'});
    usage() unless ($rh_opts->{'min_taxa'});
    return $rh_opts;
}
