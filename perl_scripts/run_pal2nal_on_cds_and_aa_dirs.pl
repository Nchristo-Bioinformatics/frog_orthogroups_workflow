#!/usr/bin/perl

use strict;
use warnings; 
use Getopt::Long;
use Data::Dumper;

our $VERSION = '0.02';

MAIN: {
    my $rh_o = get_options();
    my $aa   = $rh_o->{'aa_dir'};
    my $cds  = $rh_o->{'cds_dir'};
    my $ra_c = get_files($cds);
    my $ra_a = get_files($aa);
    my $out  = $rh_o->{'outdir'};

    unless (-d $out) {
        mkdir $out or die "cannot make $out:$!";
    }

    for (my $i = 0; $i < @{$ra_c}; $i++) {
        print "pal2nal.pl $aa/$ra_a->[$i] $cds/$ra_c->[$i] -output paml -nomismatch -nogap -codontable 1 > $out/$ra_c->[$i]_align.phy 2> $out/$ra_c->[$i].paml.stderr\n";
        system "pal2nal.pl $aa/$ra_a->[$i] $cds/$ra_c->[$i] -output paml -nomismatch -nogap -codontable 1 > $out/$ra_c->[$i]_align.phy 2> $out/$ra_c->[$i].paml.stderr\n";
        print "pal2nal.pl $aa/$ra_a->[$i] $cds/$ra_c->[$i] -output fasta -nomismatch -nogap -codontable 1 > $out/$ra_c->[$i]_align.fa 2> $out/$ra_c->[$i].fasta.stderr\n";
        system "pal2nal.pl $aa/$ra_a->[$i] $cds/$ra_c->[$i] -output fasta -nomismatch -nogap -codontable 1 > $out/$ra_c->[$i]_align.fa 2> $out/$ra_c->[$i].fasta.stderr\n";
    }
}

sub usage {
    print "usage: $0 --cds_dir=CDS_DIR --aa_dir=AA_DIR --outdir=OUTDIR [--version]\n";
    exit;
}

sub get_files {
    my $dir = shift;
    opendir DIR, $dir or die "cannot opendir $dir:$!";
    my @files = grep { !/^\./ } readdir DIR;
    return \@files;
}


sub get_options {
    my $rh_opts = {};
    my $opt_results = Getopt::Long::GetOptions(
                                 'version'   => \$rh_opts->{'version'},
                                 'cds_dir=s' => \$rh_opts->{'cds_dir'},
                                 'aa_dir=s'  => \$rh_opts->{'aa_dir'},
                                 'outdir=s'  => \$rh_opts->{'outdir'},
                                 "help"      => \$rh_opts->{'help'});

    die "$0 version $VERSION\n" if ($rh_opts->{'version'});
    print "missing --cds_dir\n" unless ($rh_opts->{'cds_dir'});
    print "missing --aa_dir\n" unless ($rh_opts->{'aa_dir'});
    print "missing --outdir\n" unless ($rh_opts->{'outdir'});
    usage() unless ($rh_opts->{'cds_dir'});
    usage() unless ($rh_opts->{'aa_dir'});
    usage() unless ($rh_opts->{'outdir'});
    return $rh_opts;
}

