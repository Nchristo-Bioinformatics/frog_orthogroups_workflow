#!/usr/bin/perl

# note: if running on your own system, need to update $HYPHY_LIB
# when using conda, this directory can be found here:
# $CONDA_ENV_PATH/lib/hyphy/TemplateBatchFiles/SelectionAnalyses
#
# Add, delete, or change hyphy analyses by editing %HYPHY_BF

$|++;

use strict;
use warnings; 
use File::Spec;
use Getopt::Long;
use Data::Dumper;

our $VERSION  = '0.02';

our $HYPHY_LIB = '/data1/vitorpavinato/00-CONDA/polar/share/hyphy/TemplateBatchFiles/SelectionAnalyses';

our $HYPHY    = 'hyphy';
our %HYPHY_BF = ('busted' => "$HYPHY_LIB/BUSTED.bf",
                 'absrel' => "$HYPHY_LIB/aBSREL.bf",
                 'meme'   => "$HYPHY_LIB/MEME.bf");

# UNTESTED: To run other tests, add from below to above
# our %HYPHY_BF = (
# 'aBSREL'        => "$HYPHY_LIB/aBSREL.bf",
# 'BranchSiteREL' => "$HYPHY_LIB/BranchSiteREL.bf",
# 'BUSTED'        => "$HYPHY_LIB/BUSTED.bf",
# 'contrast-fel'  => "$HYPHY_LIB/contrast-fel.bf",
# 'FADE'          => "$HYPHY_LIB/FADE.bf",
# 'FEL'           => "$HYPHY_LIB/FEL.bf",
# 'FitMultiModel' => "$HYPHY_LIB/FitMultiModel.bf",
# 'FUBAR'         => "$HYPHY_LIB/FUBAR.bf",
# 'MEME'          => "$HYPHY_LIB/MEME.bf",
# 'PRIME'         => "$HYPHY_LIB/PRIME.bf",
# 'RELAX'         => "$HYPHY_LIB/RELAX.bf",
# 'RELAX-Groups'  => "$HYPHY_LIB/RELAX-Groups.bf",
# 'SingleOmega'   => "$HYPHY_LIB/SingleOmega.bf",
# 'SLAC'          => "$HYPHY_LIB/SLAC.bf");


MAIN: {
    my $rh_o     = get_options();
    my $aln_dir  = $rh_o->{'aln_dir'};
    my $out_dir  = $rh_o->{'out_dir'};
    my $ra_alns  = get_files($aln_dir);
    my $pre      = $rh_o->{'pre'};
    my $tree     = $rh_o->{'tree'};

    mkdir $out_dir or die "cannot open $out_dir:$!";
    foreach my $prog (keys %HYPHY_BF) {
        mkdir "$pre.$prog" or die "cannot make $pre.$prog:$!";
        run_hyphy($ra_alns,$aln_dir,$out_dir,$prog,$pre,$tree);
    }
}

sub create_control {
    my $bf   = shift;
    my $aln  = shift;
    my $tree = shift;
    my $fg   = shift;
    my $out  = shift;

    # Hyphy requires full paths. We use File::Spec to make sure full
    my $abs_aln = File::Spec->rel2abs($aln);
    my $abs_tree = File::Spec->rel2abs($tree);
    
    open OUT, ">$out" or die "cannot open >$out";
    print OUT qq~fileToExe = "$bf";
inputRedirect = {};
inputRedirect["01"]="$abs_aln"; // codon data
inputRedirect["02"]="$abs_tree"; // tree
inputRedirect["03"]="$fg"; // Test for selection on a branch
inputRedirect["04"]=""; // complete selection

ExecuteAFile( fileToExe, inputRedirect);
~;
}

sub run_hyphy {
    my $ra_a = shift;
    my $adir = shift;
    my $odir = shift;
    my $prog = shift;
    my $pre  = shift;
    my $tree = shift;

    foreach my $aln (@{$ra_a}) {
        my $ctl = "$pre.$prog/$aln.$prog";
        create_control($HYPHY_BF{$prog},"$adir/$aln",$tree,'Foreground',$ctl);
        print "$HYPHY $ctl > $odir/$aln.$prog.out 2> $odir/$aln.$prog.err";
        system "$HYPHY $ctl > $odir/$aln.$prog.out 2> $odir/$aln.$prog.err";
        print "\n";
    }
}

sub usage {
    print "usage: $0 --aln_dir=ALIGNMENT_DIR --out_dir=OUT_DIR --tree=TREEFILE --pre=PREFIX_FOR_OUTPUT [--version]\n";
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
                                 'aln_dir=s' => \$rh_opts->{'aln_dir'},
                                 'out_dir=s' => \$rh_opts->{'out_dir'},
                                    'tree=s' => \$rh_opts->{'tree'},
                                    'pre=s'  => \$rh_opts->{'pre'});

    die "$0 version $VERSION\n" if ($rh_opts->{'version'});
    print "missing --aln_dir\n" unless ($rh_opts->{'aln_dir'});
    print "missing --out_dir\n" unless ($rh_opts->{'out_dir'});
    print "missing --tree\n" unless ($rh_opts->{'tree'});
    print "missing --pre\n" unless ($rh_opts->{'pre'});
    usage() unless ($rh_opts->{'aln_dir'});
    usage() unless ($rh_opts->{'out_dir'});
    usage() unless ($rh_opts->{'tree'});
    usage() unless ($rh_opts->{'pre'});
    return $rh_opts;
}

