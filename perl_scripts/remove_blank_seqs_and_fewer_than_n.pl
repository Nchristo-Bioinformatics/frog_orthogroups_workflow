#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Copy;

our $VERSION = 0.01;

MAIN: {
    my $rh_o     = get_options();
    my $alndir  = $rh_o->{'aln_dir'};
    my $outdir  = $rh_o->{'out_dir'};
    my $minseq  = $rh_o->{'min_seq'};
    my $ra_alns  = get_files($alndir);

    mkdir $outdir or die "cannot open $outdir:$!";
    foreach my $aln (@{$ra_alns}) {
        if (check_aln("$alndir/$aln",$minseq)) {
            File::Copy::copy("$alndir/$aln","$outdir/$aln")
                or die "copy failed: $alndir/$aln $outdir/$aln:$!";
        }
    }
}

sub check_aln {
    my $file = shift;
    my $min  = shift;
    open IN, $file or die "cannot open $file:$!";
    my $count = 0;
    my $seq = '';
    while (my $line = <IN>) {
        if ($line =~ m/^>/) {
            return 0 if ($count && !$seq);
            $count++;
            $seq = '';
        } else {
            chomp $line;
            $line =~ s/\s+//g;
            $seq .= $line;
        }
    }
    return 0 unless ($count >= $min && $seq);
    return 1;
}

sub get_files {
    my $dir = shift;
    opendir DIR, $dir or die "cannot opendir $dir:$!";
    my @files = grep { !/^\./ } readdir DIR;
    return \@files;
}

sub usage {
    print "usage: $0 --aln_dir=ALIGNMENT_DIR --out_dir=OUT_DIR --min_seq=MINIMUM_NUMBER_OF_SEQS [--version]\n";
    exit;
}

sub get_options {
    my $rh_opts = {};
    my $opt_results = Getopt::Long::GetOptions(
                                 'version'   => \$rh_opts->{'version'},
                                 'aln_dir=s' => \$rh_opts->{'aln_dir'},
                                 'out_dir=s' => \$rh_opts->{'out_dir'},
                                 'min_seq=s' => \$rh_opts->{'min_seq'});

    die "$0 version $VERSION\n" if ($rh_opts->{'version'});
    print "missing --aln_dir\n" unless ($rh_opts->{'aln_dir'});
    print "missing --out_dir\n" unless ($rh_opts->{'out_dir'});
    print "missing --min_seq\n" unless ($rh_opts->{'min_seq'});
    usage() unless ($rh_opts->{'aln_dir'});
    usage() unless ($rh_opts->{'out_dir'});
    usage() unless ($rh_opts->{'min_seq'});
    return $rh_opts;
}
