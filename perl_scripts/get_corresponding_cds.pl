#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use JFR::Fasta;

our $VERSION = 0.02;

MAIN: {
    my $cds_dir = $ARGV[0] or die "usage: $0 CDS_DIR PRUNED_AA_DIR OUTDIR\n";
    my $aa_dir  = $ARGV[1] or die "usage: $0 CDS_DIR PRUNED_AA_DIR OUTDIR\n";
    my $outdir = $ARGV[2] or die "usage: $0 CDS_DIR PRUNED_AA_DIR OUTDIR\n";

    unless (-d $outdir) {
        mkdir $outdir or die "cannot make $outdir:$!";
    }
    my $ra_cds = get_files($cds_dir,'cds');
    my $ra_aa  = get_files($aa_dir,'fa');
    my $rh_cds = get_cds_hash($cds_dir,$ra_cds);
    make_pruned_cds_files($rh_cds,$ra_aa,$aa_dir,$outdir);
}

sub make_pruned_cds_files {
    my $rh_cds = shift;
    my $ra_aa  = shift;
    my $aa_dir = shift;
    my $dir    = shift;
    foreach my $file (@{$ra_aa}) {
        $file =~ m/(.*)\.fa$/;
        my $outcds = "$dir/$1.cds.fa";
        open OUT, ">$outcds" or die "cannot open $outcds:$!";
        my $fp = JFR::Fasta->new("$aa_dir/$file");
        while (my $rec = $fp->get_record()) {
            if ($rh_cds->{$rec->{'def'}}) {
                print OUT "$rec->{'def'}\n$rh_cds->{$rec->{'def'}}\n";
            } else {
                die "cannot find $rec->{'def'}\nWARNING: $outcds is corrupt";
            }
        }
    }
}

sub get_cds_hash {
    my $dir    = shift;
    my $ra_cds = shift;
    my %cds    = ();
    foreach my $cds (@{$ra_cds}) {
        my $fp = JFR::Fasta->new("$dir/$cds");
        while (my $rec = $fp->get_record()) {
            $rec->{'def'} =~ s/:/_/g; 
            $rec->{'def'} =~ s/ .*//;
            $cds{$rec->{'def'}} = $rec->{'seq'};
        }
    }
    return \%cds;
}

sub get_files {
    my $dir = shift;
    my $suf = shift;
    opendir DIR, $dir or die "cannot opendir $dir:$!";
    my @files = grep { /\.$suf/ } readdir DIR;
    return \@files;
}
