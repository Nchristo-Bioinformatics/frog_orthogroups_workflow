#!/usr/bin/perl

$|++;

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

our $VERSION = 0.02;

MAIN: {
    my $rh_o = get_options();
    my $tree = $rh_o->{'tree'};
    my $aln_suf = $rh_o->{'aln_suf'};

    opendir DIR, "." or die "cannot opendir .:$!";
    my @files = grep { /\.$aln_suf$/ } readdir DIR;

    foreach my $file (@files) {
        if ($rh_o->{'alt'}) {
            write_ctl_file(0,'alt',$file,$tree);
            print "running: codeml ${file}.alt.ctl > ${file}.alt.out 2> ${file}.alt.err...";
            system "codeml ${file}.alt.ctl > ${file}.alt.out 2> ${file}.alt.err";
            print "\n";
        } 
        if ($rh_o->{'null'}) {
            write_ctl_file(1,'null',$file,$tree);
            system "codeml ${file}.null.ctl > ${file}.null.out 2> ${file}.null.err";
        } 
    }
}

sub write_ctl_file {
    my $fix_omega = shift;
    my $type      = shift;
    my $file      = shift;
    my $tree      = shift;

    my $ctl_file  = "${file}.$type.ctl";
    open OUT, ">$ctl_file" or die "cannot open $ctl_file:$!";
    print OUT qq~seqfile = $file
treefile = $tree
outfile = ${file}.${type}.codeml
noisy = 4
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
estFreq = 0
ndata = 1
clock = 0
aaDist = 0
model = 2
NSsites = 2
icode = 0
Mgene = 0
fix_kappa = 0
kappa = 2
fix_omega = $fix_omega
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 5
getSE = 0
RateAncestor = 0
Small_Diff = 5e-7
cleandata = 0
fix_blength = 1
method = 0
~;

}

sub usage {
    print "usage: run_codeml.pl --tree=TREEFILE --aln_suf=SUFFIX_OF_ALIGNMENT_FILES --null --alt [--version]\n";
    print "use --null to run null models and/or --alt to run alternative models\n";
    print "alignment files should be in the directory from which you are running the script.\n";
    exit;
}

sub get_options {
    my $rh_opts = {};
    my $opt_results = Getopt::Long::GetOptions(
                                 'version'   => \$rh_opts->{'version'},
                                 'tree=s'    => \$rh_opts->{'tree'},
                                 'aln_suf=s' => \$rh_opts->{'aln_suf'},
                                 'alt'       => \$rh_opts->{'alt'},
                                 'null'      => \$rh_opts->{'null'},
                                 "help"      => \$rh_opts->{'help'});
    die "$0 version $VERSION\n" if ($rh_opts->{'version'});
    print "missing --tree\n" unless ($rh_opts->{'tree'});
    print "missing --aln_suf\n" unless ($rh_opts->{'aln_suf'});
    print "--alt and/or --null is required\n" unless ($rh_opts->{'alt'});
    usage() unless ($rh_opts->{'tree'});
    usage() unless ($rh_opts->{'aln_suf'});
    usage() unless ($rh_opts->{'alt'} || $rh_opts->{'null'});
    return $rh_opts;
}
