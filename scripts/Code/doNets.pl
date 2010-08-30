#!/usr/local/bin/perl
#
# Copyright (c) 2010 Giuseppe Gallone
#
# See COPYRIGHT section in walk.pm for usage and distribution rights.
#
# Example file to document typical usage of  Bio::Homology::InterologWalk::Networks. 
# This file uses Getopt::Long for simple management of command line arguments,
# and Term::AnsiColor for clearer console output.

#USAGE perl doNets.pl -filename='yourfile.07out' -sourceorg='Mus musculus'
#or
#USAGE perl doNets.pl -filename='yourfile.06out' -sourceorg='Mus musculus'

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor;
use Bio::Homology::InterologWalk;
use Carp qw(croak);

my $work_dir = "../Data/";

my $infilename;
my $outfilename;
my $sourceorg;
my $orthtype;
my $direct;
GetOptions(
    "filename=s"  => \$infilename,
    "sourceorg=s" => \$sourceorg,
    "orthtype=s"  => \$orthtype
);

#filenames and CLAs===============================================
if ( !$infilename ) {
    $infilename = "R58_mmus_test.06out";
    print "doNets.pl: no filename specified..Trying default: $infilename\n";
}
else {
    print "doNets.pl: Using input file: $infilename\n";
}

#What input datafile are we dealing with? If it's direct interactions, we won't need to deal with orthologs, etc.
if ( ( $infilename =~ /direct/i ) or ( $infilename =~ /intact/i ) ) {
    print
      "doNets.pl: input file seems to be a direct-interactions data file..\n";
    $direct   = 1;
    $orthtype = "direct";
}
else {
    print
      "doNets.pl: input file seems to be a putative-interactions data file..\n";
    if ( !$orthtype ) {
        $orthtype = "allortho";
        print
"doNets.pl: no orthology class specified..Using default: $orthtype..\n";
    }
    else {
        print "doNets.pl: Querying orthology class: $orthtype.\n";
    }
}

if ( !$sourceorg ) {    #then use some default
    #$sourceorg = "Drosophila melanogaster";
    $sourceorg = "Mus musculus";
    #$sourceorg = 'Caenorhabditis elegans';
    print
      "doNets.pl: no source organism specified..Using default: $sourceorg\n";
}
else {
    print "doNets.pl: Using source organism: $sourceorg.\n";
}

#==================================================================

print("\n=========================\n");
my $in_path = $work_dir . $infilename;

$infilename =~ s/(.*)\.(.*)/$1\_$2/;
if ($direct) {
    $outfilename =
      $Bio::Homology::InterologWalk::NETPREF . $infilename . $Bio::Homology::InterologWalk::OUTEX_SIF;
}
else {
    $outfilename =
      $Bio::Homology::InterologWalk::NETPREF . $infilename . $Bio::Homology::InterologWalk::OUTEX_SIF;
}
my $out_path = $work_dir . $outfilename;

my $ensembl_db = 'all';
my $registry = Bio::Homology::InterologWalk::setup_ensembl_adaptor(
                                                  connect_to_db => $ensembl_db,
                                                  source_org    => $sourceorg
);
if ( !$registry ) {
    print(
"\nThere were problems setting up the connection to Ensembl. Aborting..\n"
    );
    exit;
}

#get actual network
print colored ( "Building network file from PPI/Putative PPI data...",
    'green' ), "\n";
my $start_run = time();
my $RC1        = Bio::Homology::InterologWalk::Networks::do_network(
                                                       registry    => $registry,
                                                       input_path  => $in_path,
                                                       output_path => $out_path,
                                                       source_org  => $sourceorg,
                                                       #orthology_type => $orthtype,
                                                       #expand_taxa    => 1,
                                                       #ensembl_db     => $ensembl_db,
);
if ( !$RC1 ) {
    print("There were errors. Stopping..\n");
    exit;
}
my $end_run  = time();
my $run_time = $end_run - $start_run;
print "*FINISHED* Job took $run_time seconds\n";

if ($direct) {
    $outfilename =
      $Bio::Homology::InterologWalk::NETPREF . $infilename . $Bio::Homology::InterologWalk::OUTEX_NOA;
}
else {
    $outfilename =
      $Bio::Homology::InterologWalk::NETPREF . $infilename . $Bio::Homology::InterologWalk::OUTEX_NOA;
}
$out_path = $work_dir . $outfilename;

#create cytoscape attribute files for the .sif file you just obtained
print colored ( "Creating attribute file for sif network...", 'green' ), "\n";
$start_run = time();

my $RC2       = Bio::Homology::InterologWalk::Networks::do_attributes(
                                                       registry    => $registry,
                                                       input_path  => $in_path,
                                                       output_path => $out_path,
                                                       source_org  => $sourceorg,
                                                       #label_type    => 'descr' #options are 'extname' / 'description'
);
if ( !$RC2 ) {
    print("There were errors. Stopping..\n");
    exit;
}
$registry->clear();                                                                                                                 
$end_run  = time();
$run_time = $end_run - $start_run;
print "*FINISHED* Job took $run_time seconds\n";
