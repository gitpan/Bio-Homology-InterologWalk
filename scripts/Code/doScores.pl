#!/usr/local/bin/perl
#
# Copyright (c) 2010 Giuseppe Gallone
#
# See COPYRIGHT section in walk.pm for usage and distribution rights.
#
# Example file to document typical usage of Bio::Homology::InterologWalk::Scores.
# This file uses Getopt::Long for simple management of command line arguments,
# and Term::AnsiColor for clearer console output.

#USAGE perl doScores.pl -tsvfile='yourfile.06out' -intactfile='yourfile.direct.02'

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor;
use Bio::Homology::InterologWalk;
use Carp qw(croak);

my $RC;
my $in_path;
my $out_path;
my $score_path;
my $m_mtaxa;

my $work_dir = '../Data/';
my $infilename; #actual dataset to score
my $intactfile; #direct interactions obtained with getDirectInteractions.pl
my $psimi_ont; #psi mi obo ontology from HUPO. See http://www.psidev.info/index.php?q=node/277
GetOptions( 'tsvfile=s'     => \$infilename,
            'intactfile=s'  => \$intactfile,
            'ontology=s'    => \$psimi_ont);

#filenames and files===============================================
if ( !$infilename ) {
    $infilename = 'R58_mmus_test.06out';
    print "doScores.pl: no filename specified..Trying default: $infilename\n";
}
$in_path = $work_dir . $infilename;
$infilename =~ s/(.*)\..*/$1/;

if ( !$intactfile ) {
     $intactfile = 'R58_mmus_test.direct.02';
     print "doScores.pl: no Intact direct interaction filename specified..Trying default: $intactfile\n";
}
my $intact_path = $work_dir . $intactfile;

if ( !$psimi_ont ) {
     $psimi_ont = 'psi-mi.obo';
     print "doScores.pl: no psi-mi obo ontology specified..Trying default: $psimi_ont\n";
}
my $ont_path   = $work_dir . $psimi_ont;

#output file (full datafile plus scores column)
my $out_filename = $infilename . $Bio::Homology::InterologWalk::OUTEX7;
$out_path = $work_dir . $out_filename;
#raw scores file (just for convenience)
my $score_filename = $infilename . $Bio::Homology::InterologWalk::OUTEX_SCORES;
$score_path = $work_dir . $score_filename;
#=================================================================



#==================================================================
#Computing Mean Multiple taxa score
#==================================================================
#WARNING: this might take a long time
$m_mtaxa = Bio::Homology::InterologWalk::Scores::compute_multiple_taxa_mean(
                                                            ds_size   => 10,          #number of ids per dataset, eg 500
                                                            ds_number => 2,           #max is 7, equal to the number of taxa
                                                            datadir   => $work_dir    #for the path
                                                            );
if ( !$m_mtaxa ) {
    print "There were errors. Stopping..\n";
    exit;
}

#create a Go:Parser graph to explore the ontology.
my $onto_graph = Bio::Homology::InterologWalk::Scores::parse_ontology($ont_path);
if ( !$onto_graph ) {
    print "There were errors. Stopping..\n";
    exit;
}

#3)#Process the direct interactions data file to retrieve the mean scores for
# interaction type, interaction detection method, experimental method and multiple detection method
#all of these can be obtained from a dataset of direct interactions (ie no orthology projections)
my ( $m_em, $m_it, $m_dm, $m_mdm ) =
  Bio::Homology::InterologWalk::Scores::get_mean_scores( $intact_path, $onto_graph );


#4) compute actual scores
print colored ( "Computing putative interaction scores...", 'green' ), "\n";
my $start_run = time;
$RC = Bio::Homology::InterologWalk::Scores::compute_scores(
                                             input_path        => $in_path,
                                             score_path        => $score_path,
                                             output_path       => $out_path,
                                             term_graph        => $onto_graph,
                                             meanscore_em      => $m_em,
                                             meanscore_it      => $m_it,
                                             meanscore_dm      => $m_dm,
                                             meanscore_me_dm   => $m_mdm,
                                             meanscore_me_taxa => $m_mtaxa
);
if ( !$RC ) {
    print "There were errors. Stopping..\n" ;
    exit;
}
my $end_run  = time;
my $run_time = $end_run - $start_run;
print "\n\n*FINISHED* Job took $run_time seconds\n";
