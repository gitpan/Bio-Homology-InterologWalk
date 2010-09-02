package Bio::Homology::InterologWalk;

use 5.008006;
use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use REST::Client;
use GO::Parser;
use DBI;
use Carp qw(croak);

require Exporter;
use base qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [ qw(
        
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.03';



#################### main pod documentation begins ###################

=head1 NAME

Bio::Homology::InterologWalk - Retrieve, score and visualize putative Protein-Protein Interactions through the orthology-walk method

=head1 VERSION

This document describes version 0.02 of Bio::Homology::InterologWalk released August 31st, 2010

=head1 SYNOPSIS

  use Bio::Homology::InterologWalk;

First, obtain Intact Interactions for the dataset (see example in C<getDirectInteractions.pl>):


  #get a registry from Ensembl
  my $registry = Bio::Homology::InterologWalk::setup_ensembl_adaptor(
                                                     connect_to_db  => $ensembl_db,
                                                     source_org     => $sourceorg,
                                                     verbose        => 1
                                                     );
  
  
  #query direct interactions
  $RC = Bio::Homology::InterologWalk::Direct::get_direct_interactions(
                                                      registry         => $registry,
                                                      source_org       => $sourceorg,
                                                      input_path       => $in_path,
                                                      output_path      => $out_path,
                                                      url              => $url,
                                                      );




do some postprocessing (see L</do_counts> and L</extract_unseen_ids> ) and then run the actual interolog walk on the dataset
with the following sequence of three methods.


get orthologues of starting set:

  $RC = Bio::Homology::InterologWalk::get_forward_orthologies(
                                              registry        => $registry,
                                              ensembl_db      => $ensembl_db,
                                              input_path      => $in_path,
                                              output_path     => $out_path,
                                              source_org      => $sourceorg,
                                              dest_org        => $destorg,
                                              );



add interactors of  orthologues found by C<get_forward_orthologies()>:



  $RC = Bio::Homology::InterologWalk::get_interactions(
                                       input_path    => $in_path,
                                       output_path   => $out_path,
                                       url           => $url,
                                       );


add orthologues of interactors found by C<get_interactions()>:

  $RC = Bio::Homology::InterologWalk::get_backward_orthologies(
                                               registry    => $registry,
                                               ensembl_db  => $ensembl_db,
                                               input_path  => $in_path,
                                               output_path => $out_path,
                                               error_path  => $err_path,
                                               source_org  => $sourceorg,  
                                               );

do some postprocessing (see L</remove_duplicate_rows>, L</do_counts>, L</extract_unseen_ids>)
and then optionally compute a composite score
for the putative interactions obtained:

  $RC = Bio::Homology::InterologWalk::Scores::compute_scores(
                                             input_path        => $in_path,
                                             score_path        => $score_path,
                                             output_path       => $out_path,
                                             term_graph        => $onto_graph,
                                             meanscore_it      => $m_it,
                                             meanscore_dm      => $m_dm,
                                             meanscore_me_dm   => $m_mdm,
                                             meanscore_me_taxa => $m_mtaxa
                                             );


get some networks and network attributes which you can then visualise with cytoscape

   $RC = Bio::Homology::InterologWalk::Networks::do_network(
                                            registry        => $registry,
                                            input_path      => $in_path,
                                            output_path     => $out_path,
                                            source_org      => $sourceorg
                                            );
                                               
   $RC = Bio::Homology::InterologWalk::Networks::do_attributes(
                                               registry      => $registry,
                                               input_path    => $in_path,
                                               output_path   => $out_path,
                                               source_org    => $sourceorg,
                                               label_type    => 'external name'
                                               );

I<The synopsis above only lists the major methods and parameters.>

=head1 DESCRIPTION



A common activity in computational biology is to mine protein-protein interactions 
from publicly available databases to build I<Protein-Protein Interaction> (PPI) datasets.  
In many instances, however, the number of experimentally obtained annotated PPIs is very scarce
and it would be helpful to enrich the experimental dataset with high-quality, computationally-inferred PPIs.
Such computationally-obtained dataset can extend, support or enrich experimental PPI datasets, and are
of crucial importance in high-throughput gene prioritization studies, i.e. to drive hypotheses and restrict the
dimensionality of functional discovery problems.
This Perl Module, C<Bio::Homology::InterologWalk>, is aimed at building 
putative PPI datasets on the basis of  a number of comparative biology paradigms: the module implements 
a collection of computational biology algorithms based on the concept of "orthology projection". 
If interacting proteins A and B in organism X have orthologues A' and B' in organism Y, under certain
conditions one can assume that the interaction will be conserved in organism Y, i.e. the A-B 
interaction can be "projected through the orthologies" to obtain a putative A'-B' interaction.
The pair of interactions (A-B) and (A'-B') are named "Interologs".



C<Bio::Homology::InterologWalk> collects, analyses and collates gene orthology data
provided by the Ensembl Consortium as well as PPI data provided by
EBI Intact. It provides the user with the possibility of rating the quality and reliability of 
the putative interactions collected by means of  confidence scores,
and optionally outputs network representations of the datasets, 
compatible with the biological network representation standard, Cytoscape.

=head1 USAGE

=head2 Rationale


                              \EBI Intact API/
         .--------------.            |             .-------------.
     (2) | A(e.g. mouse)|<------------------------>|   B(mouse)  |  (3)
         `--------------'          <PPI>           `-------------'
                ^                                         |
   /Ensembl\    | <Orthology>                 <Orthology> | \ Ensembl /
  / Compara \   |                                         |  \Compara/
 /    Api    \  |                                         |   \ Api /
                |                                         | 
         .--------------.                           .-------------.
     (1) | A'(e.g. fly) |. . . . . . . . . . . . .  |   B'(fly)   | (4)
         `--------------'     [SCORED]PUTATIVE PPI  `-------------'
                    (Output of Bio::Homology::InterologWalk)


In order to carry out an interolog walk we start with a set of gene identifiers in one organism of interest (1).
We query those ids against a number of comparative biology  databases to retrieve a list of orthologues for the
gene ids of interest, in one or more species (2). In the next step we rely instead on PPI databases 
to retrieve the list of available interactors for the protein ids obtained in (2). The output at this stage consists of
a list of interactors of the orthologues of the initial gene set, plus several fields of ancillary data (whose importance will
be explained later) (3). 
In the last step of the process we will need to project the interactions in (3) - again using orthology data - back
to the original species of interest. 
The final output is a list of B<putative interactors> for the initial gene set, plus several
fields of supporting data. 

C<Bio::Homology::InterologWalk> provides three main functions to carry out the basic walk, C<get_forward_orthologies()> , C<get_interactions()> and  
C<get_backward_orthologies()>. These functions must be called strictly in sequential order in the user's script, as they process, analyse and attach data to the output 
in a pipeline-like fashion, i.e. working on the output of the preceding function.


=over 4

=item  get_forward_orthologies

This methods queries the initial gene list against one or more Ensembl DBs (using the Ensembl Perl API) and retrieves their orthologues, 
plus a number of ancillary data fields ( conservation data, distance from ancestor, orthology type, etc) 

=item get_interactions

This queries the orthology list built in the previous stage against PSICQUIC-enabled PPI DBs using Rest. 
This step will enrich the dataset built through C<get_forward_orthologies> with the interactors of those orthologues, if any, plus ancillary data (including several
parameters describing the quality, nature and origin of the annotated interaction).

=item get_backward_orthologies

This queries the interactor list built in the previous stage against one or more Ensembl DBs (again  using the Ensembl Perl API)
to find orthologues back in the original species of interest. It will also adds a number of supplementary information fields, specularly to what done in C<get_forward_orthologies>.  

=back

The output of this sequence of subroutines will be a TSV file containing zero or more entries, closely resembling  the MITAB tab delimited data exchange format from the  HUPO PSI (Proteomics Standards Initiative).
Each row in the data file represents a binary putative interaction, plus currently 37 supplementary data fields.

This basic output can then be further processed with the help of other methods in the module: one can scan the results to compute counts, 
to check for duplicates, to verify the presence of new gene ids that were not present in the original dataset and save them in another datafile, and so on.

Most importantly, the user could need to further process the putative PPIs dataset to do one or more of the following:


=over

=item 1.

Compute a global confidence score to obtain a metric for the reliability of the each binary putative interaction

=item 2.

Extract the binary putative PPIs from the dataset and save them in a format compatible with Cytoscape. This helps providing a visual quality to the result as 
one could then apply network analysis tools to discover motifs, clusters, well-connected subnetworks, look for GO functional enrichment, and more. 
The format chosen for the network representation of the dataset is currently C<.sif>. (see http://cytoscape.wodaklab.org/wiki/Cytoscape_User_Manual#Supported_Network_File_Formats) 
The generation of node attributes is also possible, to allow for visualisation of node tags in terms of simpler human readable labels instead of database IDs.

=item 3.

Obtain a dataset of experimental/direct PPIs (i.e. just plain interactors, 
no orthology mapping across other taxa involved) from the gene list used as the input to the orthology walk. The reasons why this might be useful are several. 
The user might want to compare this dataset with the putative PPI dataset also generated by the module to see if/where the two overlap, 
what is the intersection/difference set, and more. See L</get_direct_interactions> for documentation relative to this function. Please also notice a dataset of 
direct interactions will also be pre-requisite if the user intends to compute confidence values for the putative PPI dataset: the direct PPI dataset is required to
compute score normalisation means.

=back


=head1 EXAMPLES

In order to demonstrate one way of using the module, four example perl scripts are provided in the C<scripts/Code> directory. 
Each sample script utilises the module and uses/reuses subroutines in a pipeline fashion. The workflow suggested with the scripts is as follows:

B<User Input>: a textfile containing one gene ID per row. All gene IDs must belong to the same species. All gene IDs must be current Ensembl gene IDs.


=over

=item 1. B<Mine Direct Interactions.>

Generate a dataset of direct PPIs based on the input ID list. See example in C<getDirectInteractions.pl>

=item 2. B<Run the basic Interolog-Walk Pipeline.>

Generate dataset of projected putative PPIs following the paradigm explained earlier. Do some postprocessing on the dataset. See example in C<doInterologWalk.pl>

=item 3. B<Compute confidence scores for putative PPIs.>

Score the dataset obtained in (2.) using the dataset obtained in (1.) to normalise the score components values. See example in C<doScores.pl>

=item 4. B<Extract network and attributes for the two PPI datasets.>

For each of the two datasets obtained from (1) and (2) (putative PPIs) or from (1) and (3) (scored putative PPIs) extract a text file containing a network representation and 
a text file of node attributes. 

See example in C<doNets.pl>

=back


=head1 DEPENDENCIES 

C<Bio::Homology::InterologWalk> relies on the following prerequisite packages.

=head2 Ensembl API

The Ensembl project is currently branched in two sub-projects:

=over

=item The Ensembl Vertebrates project 

This is of interest to you if you work with vertebrate genomes (although it also includes data from a few non-vertebrate common model organisms).
See http://www.ensembl.org/index.html for further details.

=item The Ensembl Genomes project 

This utilises the Ensembl software infrastructure (originally developed in the Ensembl Core project) to provide 
access to genome-scale data from non-vertebrate species. This is of interest to you if your species  is a non-vertebrate, or if
your species  is a vertebrate but you I<also want to obtain results mapped from non-vertebrates>. C<Bio::Homology::InterologWalk> currently only supports
the B<metazoa> sub-site from the Ensembl Genomes Project. See http://metazoa.ensembl.org/index.html for further details.

=back 



B<IMPORTANT> You will need to decide which Ensembl-DB set you will need B<prior> to installing C<Bio::Homology::InterologWalk>.  
The module requests that 

Ensembl API Version == Ensembl-DB set version. 

This means that if you install e.g. API V.58, 
you will only be able to get data from Ensembl Vertebrates / Metazoa databases V. 58. As the EnsemblGenomes DB releases are 
B<one version behind> the Ensembl Vertebrate DB release, if you install the bleeding-edge Ensembl Vertebrate API, I<a matching EnsemblGenomes DB release might
not be available yet>: you will still be able to use C<Bio::Homology::InterologWalk> to run an orthology walk using exclusively Ensembl Vertebrate DBs, but you
will get an error if you try to choose metazoan databases. See L</setup_ensembl_adaptor> for further information.


Therefore, before installing C<Bio::Homology::InterologWalk>, you are faced with the following choice: 

=over

=item a)

If you are exclusively interested in B<vertebrates> (plus the few
non-vertebrate model organisms still present in Ensembl Vertebrates) then obtain the APIs and
set up the environment by following the steps described on the Ensembl Vertebrates API installation pages:

http://www.ensembl.org/info/docs/api/api_installation.html

or alternatively

http://www.ensembl.org/info/docs/api/api_cvs.html

This option allows you to get the B<most recent> datasets provided by Ensembl Core. However, you might not be able to query EnsemblCompara data.

=item b) 

If you are interested in querying/getting back data from vertebrate + metazoan genomes, then obtain the APIs
and set up the environment by following the steps described on the Ensembl Metazoa API installation pages: 
(this allows you to query across a wider selection of taxa)

http://metazoa.ensembl.org/info/docs/api/api_installation.html

or alternatively

http://metazoa.ensembl.org/info/docs/api/api_cvs.html

This option will probably not use the most recent API+DBs, but will guarantee functionality across both Vertebrate and Metazoan genomes.

=back

Option (b) is the B<recommended> one.

NOTE 1: All the API components  (C<ensembl>, C<ensembl-compara>, C<ensembl-variation>, C<ensembl-functgenomics>) are required.

NOTE 2: The module has been tested on Ensembl Vertebrates API & DB v. 58 and v. 59 and EnsemblGenomes API & DB  v. 5 (58). 

=head2 Bioperl

Ensembl provides a customised Bioperl installation tailored to its API, v. 1.2.3. 
Should version 1.2.3 be no more available through Ensembl, please obtain release 1.6.x from CPAN. (while not officially supported by the
Ensembl Project it will work fine when using the API within the scope of the present module)

=head2 Additional Perl Modules

The following modules (including all dependencies) from CPAN are also required:

=over

=item 1. C<REST::Client>

=item 2. C<GO::Parser>

=item 3. C<DBD::CSV> (requires Perl DBI)

=item 4. C<String::Approx>

=back

I<See the README file for further information.>


=head1 INTERFACE

=cut


#################### main pod documentation ends ###################






#globals#############################
my $ENSEMBLIDFAILED;

our $ERREX          = '.04err';
our $OUTEX1         = '.01out';
our $OUTEX2         = '.02out';
our $OUTEX3         = '.03out';
our $OUTEX4         = '.04out';
our $OUTEX5         = '.05out';
our $OUTEX6         = '.06out';
our $OUTEX7         = '.07out';
our $OUTEX_NEW      = '.newID';
our $OUTEX_NEW_DIR  = '.direct.newID';
our $OUTEX_SIF      = '.sif';
our $OUTEX_SCORES   = '.scores';
our $OUTEX_NOA      = '-name.noa';
our $NETPREF        = 'NET-';
our $INTACTEX0      = '.direct.00';
our $INTACTEX1      = '.direct.01';
our $INTACTEX2      = '.direct.02';
my $ENSEMBLVERSION  =  Bio::EnsEMBL::Registry->software_version();
our $VERSIONEX      = "R" . "$ENSEMBLVERSION" . "_";

#headers#######################################
my $FN_initid                 =    'INIT_ID';
my $FN_orthologue_id          =    'ORTHOLOGUE_ID';
my $FN_oname_1                =    'ORTHOLOGUE_NAME_1';     
my $FN_odesc_1                =    'ODESCRIPTION_1';   
my $FN_opi_1                  =    'OPI_1';  
my $FN_dnds_1                 =    'DN_DS_1';     
my $FN_nndist_1               =    'NODE_NODE_DIST_1'; 
my $FN_fsa_orig_species_1     =    'FSA_INITIAL_1';
my $FN_fsa_dest_species_1     =    'FSA_ORTHOLOG_1';

my $HEADER_FWD_ORTH = join("\t", $FN_initid,$FN_orthologue_id, 
                                 $FN_oname_1, $FN_odesc_1, 
                                 $FN_opi_1, $FN_dnds_1, 
                                 $FN_nndist_1, $FN_fsa_orig_species_1, 
                                 $FN_fsa_dest_species_1);

my $FN_interaction_id         =    'INTERACTION_ID';
my $FN_acc_numb_a             =    'ACCESSION_NUMBER_A';    
my $FN_acc_numb_b             =    'ACCESSION_NUMBER_B';
my $FN_alt_id_a               =    'ALT_ID_A';
my $FN_alt_id_b               =    'ALT_ID_B';
my $FN_name_a                 =    'NAME_A'; 
my $FN_name_b                 =    'NAME_B';
my $FN_taxon_a                =    'TAXON_A';
my $FN_taxon_b                =    'TAXON_B';
my $FN_pub                    =    'PUBLICATION';
my $FN_int_type               =    'INTERACTION_TYPE'; 
my $FN_det_method             =    'DET_METHOD';
my $FN_exp_method             =    'EXP_METHOD';

my $HEADER_INT = join("\t", $FN_interaction_id, $FN_acc_numb_a, $FN_acc_numb_b,
                            $FN_alt_id_a, $FN_alt_id_b, $FN_name_a, $FN_name_b,
                            $FN_taxon_a, $FN_taxon_b, $FN_pub, $FN_int_type,
                            $FN_det_method, $FN_exp_method);

my $FN_interactor_id          =    'INTERACTOR';  
my $FN_oname_2                =    'ORTHOLOGUE_NAME_2';
my $FN_odesc_2                =    'ODESCRIPTION_2';
my $FN_opi_2                  =    'OPI_2';
my $FN_dnds_2                 =    'DN_DS_2';
my $FN_nndist_2               =    'NODE_NODE_DIST_2'; 
my $FN_fsa_orig_species_2     =    'FSA_INITIAL_2';    
my $FN_fsa_dest_species_2     =    'FSA_ORTHOLOG_2';

my $HEADER_BWD_ORTH = join("\t", $FN_interactor_id, $FN_oname_2,
                                 $FN_odesc_2, $FN_opi_2, $FN_dnds_2,
                                 $FN_nndist_2, $FN_fsa_orig_species_2, $FN_fsa_dest_species_2);

my $FN_compound_score         =    'COMPOUND_SCORE';

my $HEADER_SCORE              =    $FN_compound_score; 

my $FN_multiple_dm            =    'MULTIPLE_DM'; 
my $FN_multiple_taxa          =    'MULTIPLE_TAXA';    
my $FN_times_seen             =    'TIMES_SEEN';
my $FN_same_init_final_id     =    'SAME_INIT_FINAL_ID';    
my $FN_same_walk_ints         =    'SAME_WALK_INTERACTORS';

my $HEADER_MISC = join("\t", $FN_multiple_dm, $FN_multiple_taxa, $FN_times_seen,
                             $FN_same_init_final_id, $FN_same_walk_ints);
                                        
#my $HEADER_SIF_SCORES       =    join("\t", $FN_initid, 'INT_TYPE', $FN_compound_score, $FN_interactor_id);
# when no score is provided and we want one edge per taxon
my $HEADER_SIF               =    join("\t", $FN_initid, 'INT_TYPE', $FN_interactor_id);
#REMOVED IT FOR NOW
#my $FN_alt_ids               =    "OTHER_IDS";

#header for the direct interaction data file
my $HEADER_DIRECT             =    join("\t", $FN_initid, $FN_interaction_id,
                                              $FN_acc_numb_a, $FN_acc_numb_b,
                                              $FN_alt_id_a, $FN_alt_id_b,
                                              $FN_name_a, $FN_name_b,
                                              $FN_taxon_a, $FN_taxon_b,
                                              $FN_pub, $FN_int_type,
                                              $FN_det_method, $FN_exp_method,
                                              $FN_interactor_id
                                           #, $FN_alt_ids
                                              );
#######################################################

=head2 setup_ensembl_adaptor

 Usage     : $registry = Bio::Homology::InterologWalk::setup_ensembl_adaptor(
                                                             connect_to_db   => $ensembl_db,
                                                             source_org      => $sourceorg,
                                                             dest_org        => $destorg,
                                                             verbose         => 1
                                                             );
 Purpose   : This subroutine sets up the registry for connection to the Ensembl API and also gets 
             a species-dependent adaptor out of it
 Returns   : An Ensembl Registry object if successful, undefined in all other cases
 Argument  : -connect_to_db: ensembl db to connect to. Choices currently are:
                 a. 'multi' :  vertebrate compara (see http://www.ensembl.org/)
                 b. 'pan_homology' : pan taxonomic compara db, a selection of species from both 
                    Ensembl Compara and EnsemblGenomes Compara
                    (see http://nar.oxfordjournals.org/cgi/content/full/38/suppl_1/D563 )
                 c. 'metazoa' : EnsemblGenomes compara, metazoa db 
                    (see http://metazoa.ensembl.org/index.html). 
                 d. 'all'  : multi + metazoa.
                 Default is 'multi'.
             -source_org: the initial species for the interolog walk. This MUST match with your 
              choice of db. Exception is raised if not
             -(OPTIONAL) dest_org: the destination species to use for the interolog walk. This 
              MUST exist in your choice of db. "all" 
               chooses all the taxa offered by Ensemlb in that DB. Default is 'all'
             -(OPTIONAL) verbose: boolean, shows/hides connection info provided by Ensembl. 
              Default is '0' 
 Throws    : -
 Comment   : Currently the FULL SCIENTIFIC NAME of both the source organism and the destination 
             organism, as specified in Ensembl, is required.
             E.g.: 'Homo sapiens', 'Mus musculus', 'Drosophila melanogaster', etc. 
             Soon to be expanded to support short mnemonic names 
             (e.g.: 'Mmus' instead of 'Mus musculus')

See Also   : 

=cut


sub setup_ensembl_adaptor{
     my %args = @_;
     my $db         = $args{connect_to_db};
     my $sourceorg  = $args{source_org};
     my $destorg    = $args{dest_org};
     my $verbose    = $args{verbose};
     
     #defaults----
     $db = 'multi' if(!$db);
     $destorg = 'all' if(!$destorg); 
     $verbose = 0 if(!$verbose);
     #------------
     
     if(!$sourceorg){
          print("setup_ensembl_adaptor(): no source organism specified. Aborting..\n");
          return;
     }
     my $registry = 'Bio::EnsEMBL::Registry';
     $db = lc($db);

     if($db eq 'multi'){ #connect to main db default
          $registry->load_registry_from_db(
                    -host          => 'ensembldb.ensembl.org',
                    -user          => 'anonymous',
                    -verbose       => $verbose,
                    );  
          my @dba = @{ $registry ->get_all_DBAdaptors() };
          if(scalar(@dba) eq 0){ #it should never get in here
               print "\n\nsetup_ensembl_adaptor() - ERROR: no Ensembl Vertebrates dbs available for connection.\n";
               print "API VERSION: V. $ENSEMBLVERSION\n";
               print "*PLEASE make sure suitable databases for this API version have been released by Ensemb.\n";
               print "Aborting..\n";
               return;
          }
     }elsif(($db eq 'pan_homology') || ($db eq 'metazoa') ){
          $registry->load_registry_from_db(
                    -host          =>   'mysql.ebi.ac.uk',
                    -port          =>   4157,
                    -user          =>   'anonymous',
                    -verbose       =>   $verbose,
                    );
          my @dba = @{ $registry ->get_all_DBAdaptors() };
          if(scalar(@dba) eq 0){# it should get in here if you have eg api v.59 and ensemblgenomes project is at v.58
               print "\n\nsetup_ensembl_adaptor() - ERROR: no EnsemblGenomes Metazoa dbs available for connection.\n";
               print "Your installed API is version: V. $ENSEMBLVERSION.\n";
               print "This API might be too new for the current EnsemblGenomes Metazoa db release.\n";
               print "Please install an earlier API version from metazoa.ensembl.org/info/docs/api/api_cvs.html\n";
               print "Aborting..\n";
               return;
          }
     }elsif($db eq 'all'){ #vertebrates + metazoa
          $registry->load_registry_from_multiple_dbs(
          {    #VERTEBRATES
                  -host       => 'ensembldb.ensembl.org',
                  -user       => 'anonymous',
                  -verbose    => $verbose,
                  -VERSION    => $ENSEMBLVERSION
              },
          {     #METAZOA
                     -host    => 'mysql.ebi.ac.uk',
                     -user    => 'anonymous',
                     -port    => 4157,
                     -verbose => $verbose,
                     -VERSION => $ENSEMBLVERSION
              }
          );
          #If only vertebrate ensembl is correctly registered while genomes is not, this will not fail,
          #while I want it to fail now and tell the use to install an older api (compatible with genomes)
          #i can retrieve one significant dbadaptor per connection (multi/metazoa) using a species which I know
          #must be in there (human / worm)
          my $dba_multi = $registry->get_DBAdaptor("Homo sapiens", "core");
          my $dba_metazoa = $registry->get_DBAdaptor("Anopheles gambiae", "core");
          unless($dba_multi && $dba_metazoa){
               print "\n\nsetup_ensembl_adaptor() - ERROR: some databases are not available for connection:\n";
               print "Your installed API is version: V. $ENSEMBLVERSION.\n";
               if(!$dba_multi){
                    print "Vertebrate databases not found.\n";
               }
               if(!$dba_metazoa){
                    print "Metazoa databases not found.\n";
                    print "This API might be too new for the current EnsemblGenomes Metazoa db release.\n";
                    print "Please install an earlier API version from metazoa.ensembl.org/info/docs/api/api_cvs.html\n";
               }
               print "Aborting..\n";
               return;
          }

     }else{
          print("setup_ensembl_adaptor(): db identifier: $db  not recognised. Please double check. Exiting..\n");
          return;
     }
     
     #We have to check the species chosen by the user against his choice of dbs..
     if(!($registry->get_adaptor($sourceorg, "core", "Gene"))){
          print("setup_ensembl_adaptor(): The selected source organism ($sourceorg) is not available in the selected Ensembl Compara DB ( $db ).\n");
          print("*Please make sure you choose a meaningful combination of species and db.* Aborting..\n");
          return;
     }
     unless($destorg eq 'all'){ #all will mean "all the species included in the selected db"
          if(!($registry->get_adaptor($destorg, 'core', 'Gene'))){
               print("setup_ensembl_adaptor(): The selected destination organism ($destorg) is not available in the selected Ensembl Compara DB ( $db ).\n");
               print("*Please make sure you choose a meaningful combination of species and db. Aborting..\n");
               return;
          }
     }
    return $registry;
}



=head2 remove_duplicate_rows

 Usage     : $RC = Bio::Homology::InterologWalk::remove_duplicate_rows(
                                                       input_path   => $in_path,
                                                       output_path  => $out_path,
                                                       header       => 'standard',
                                                       );
 Purpose   : This is used to clean up a TSV data file of duplicate entries. 
             This routine will make sure no such duplicates are kept. A new datafile 
             is built. The number of unique data rows is updated. 
 Returns   : success/error
 Argument  : -input_path : path to input file. Input file for this subroutine is 
              a TSV file of PPIs. It can be one of the following two:
                 1. the output of get_backward_orthologies(). In this case please 
                    specify 'standard' header below.
                 2. the output of get_direct_interactions(). In this case please 
                    specify 'direct' header below.             
             -output_path : where you want the routine to write the data. Data is in 
              TSV format. 
             -(OPTIONAL)header : Header type is one of the following:
                 1. 'standard': when the routine is used to clean up an interolog-walk 
                    file (the header will be longer)
                 2. 'direct':   when the routine is used to clean up a file of real db 
                    interactions (the header is shorter)
               No field provided: default is 'standard'
 Throws    : -
 Comment   : -

See Also   : L</get_backward_orthologies>, L</get_direct_interactions>

=cut


sub remove_duplicate_rows{
     my %args = @_;
     my $in_path          = $args{input_path};
     my $out_path         = $args{output_path};
     my $header_type      = $args{header};
     
     $header_type = 'standard' if(!$header_type);
     #MANAGE FILES
     my $infile_handle = _setup_dbi_connection($in_path);
     open (my $out_data,  q{>}, $out_path) or croak("Unable to open $out_path : $!");
     #============
     
     #huge hash that contains all the entries. Should be improved upon
     my %seen = ();
     my $duplicates = 0;
     
     #header====
     if($header_type eq 'standard'){
          print $out_data $HEADER_FWD_ORTH, "\t", $HEADER_INT, "\t", $HEADER_BWD_ORTH, "\n";
     }elsif( ($header_type eq 'direct') ){
          print $out_data $HEADER_DIRECT, "\n";
     }else{
          print("remove_duplicate_rows(): header: $header_type not recognised, aborting..");
          return;
     }
     #==========
     my ($query) = "SELECT * FROM int";

     my $rowCount = $infile_handle->selectrow_array(
                                  qq{
                                        SELECT count(*)
                                        FROM int
                                    });
    if($rowCount == 0){
       print("remove_duplicate_rows(): input file appears to be empty. Aborting..\n");
       return;
    }          
     my $sth = $infile_handle->prepare($query);
     $sth->execute() or die "Cannot execute: " . $sth->errstr();
     
     while (my @data = $sth->fetchrow_array()) {
          my $datarow = join("\t", @data);
     
          if(!$seen{$datarow}){
               $seen{$datarow} = 1;
               print $out_data $datarow, "\n";
          }else{
               $duplicates += 1;
          }
     }
     
     my $number = keys %seen;
     print("\n======================\n");
     print("remove_duplicate_rows(): number of rows in input file: $rowCount\n");
     print("remove_duplicate_rows(): number of duplicate rows: $duplicates\n");
     print("remove_duplicate_rows(): new data files contains $number UNIQUE interactions.\n");
     
     $sth->finish();
     $infile_handle->disconnect();
     close($out_data);
     return 1;
}


=head2 get_forward_orthologies

 Usage     : $RC = Bio::Homology::InterologWalk::get_forward_orthologies(
                                                         registry      => $registry,
                                                         ensembl_db    => $ensembl_db,
                                                         input_path    => $in_path,
                                                         output_path   => $out_path,
                                                         source_org    => $sourceorg,
                                                         dest_org      => $destorg,
                                                         hq_only       => 1
                                                         no_output     => 0
                                                         );
 Purpose   : This is the core function to perform the orthology retrieval step of the 
             Interolog mapping algorithm. It will set up some important Ensembl components 
             and then proceed with the composition/computation of the values
 Returns   : success/error code
 Argument  : -registry object to connect to ensembl
             -ensembl db to connect to. Choices currently are:
                 a. 'multi' :  vertebrate compara (see http://www.ensembl.org/)
                 b. 'pan_homology' : pan taxonomic compara db, a selection of species 
                    from both Ensembl Compara and Ensembl Genomes 
                    (see http://nar.oxfordjournals.org/cgi/content/full/38/suppl_1/D563 )
                 c. 'metazoa' : ensembl compara genomes, metazoa db 
                    (see http://metazoa.ensembl.org/index.html). 
                 d. 'all'  : multi + metazoa.
                 Default is 'multi'.
             -input_path : path to input file. Input file MUST be a text file with one entry
              per row, each entry containing an up-to-date gene ID recognised by the Ensembl 
              consortium (http://www.ensembl.org/) followed by a new line char.
             -output_path : where you want the routine to write the data. Data is in TSV 
              format.
             -source organism name (eg: 'Mus musculus')
             -(OPTIONAL)destination organism name (eg 'Drosophila melanogaster'). Set this is 
              if you want to carry out the mapping through one 
              specific species, rather than all those available in Ensembl. Default : 'all'
             -(OPTIONAL)hq_only: discards one-to-many, many-to-one, many-to-many orthologues. 
              Only keeps one-to-one orthologues, i.e. where no duplication event has happened 
              after the speciation. One-to-one orthologues are ideally associated with higher 
              functional conservation (while paralogues often cause neo/sub-functionalisation). 
              For further information see
              http://www.ensembl.org/info/docs/compara/homology_method.html 
             -(OPTIONAL) no_output :  suppresses screen output. Used for clearer output during 
              test. Default is 0.
 Throws    : -
 Comment   : 1)Currently the FULL SCIENTIFIC NAME of both the source species and the destination 
               species, as specified in Ensembl, is required.
               E.g.: 'Homo sapiens', 'Mus musculus', 'Drosophila melanogaster', etc. 
               Soon to be expanded to support short mnemonic names (e.g.: 'Mmus' instead of 
               'Mus musculus')
             2)EXPERIMENTAL: early support for human readable gene names in the input file has been 
               added. Such gene names will be checked against Ensembl so they must be recognisable 
               by it.


See Also   : 

=cut


sub get_forward_orthologies{
     my %args = @_;

     my $registry       = $args{registry};
     my $db             = $args{ensembl_db};
     my $sourceorg      = $args{source_org};
     my $destorg        = $args{dest_org};
     my $in_path        = $args{input_path};
     my $out_path       = $args{output_path};
     my $onetoone_only  = $args{hq_only};
     my $no_output      = $args{no_output};
     
     #checks
     if(!$registry){
       print("get_forward_orthologies(): no registry object found. Aborting..\n");
       return;
     }
     if(!$sourceorg){
      print("get_forward_orthologies(): no source organism specified. Aborting..\n");
       return;
     }
     $destorg = 'all' if(!$destorg);
     
     $sourceorg = ucfirst(lc $sourceorg);
     $destorg = ucfirst(lc $destorg);
     
     #MANAGE FILES
     open (my $in_data,  q{<}, $in_path) or croak("Unable to open $in_path : $!");
     open (my $out_data,  q{>}, $out_path) or croak("Unable to open $out_path : $!");
     #============
     
     #defaults
     if(!$db){
          printf("get_forward_orthologies(): no db class specified. Using default: multi (vertebrates)\n") unless($no_output);
          $db = "multi" ;
     }
     
     my @ensembl_dbs;
     my $counter; 
     my $counter_db = 0; 
     my $global_count = 0;
     
     if($db eq "all"){
          push(@ensembl_dbs, "multi");
          push(@ensembl_dbs, "metazoa");
     }else{
          push(@ensembl_dbs,$db);
     }
     #set up header
     print $out_data $HEADER_FWD_ORTH, "\n";
     
     my $source_species_gene_adaptor = $registry->get_adaptor($sourceorg, 'core', 'Gene');
     #some control on the source organism name
     if(!$source_species_gene_adaptor){
          print("get_forward_orthologies(): source organism string: $sourceorg not recognised by Ensembl. Please double check. Aborting..\n");
          return;
     }
     
     foreach my $ensembl_db (@ensembl_dbs){
          print("\n----Querying Ensembl Compara ($ensembl_db) for orthologues----\n") unless($no_output);
          my $orthologues_mlss;
          my %genome_taxon_ids;
          $counter_db = 0;
          
          my $member_adaptor = $registry->get_adaptor($ensembl_db, 'compara', 'Member');
          my $homology_adaptor = $registry->get_adaptor($ensembl_db, 'compara', 'Homology');
          my $proteintree_adaptor = $registry->get_adaptor($ensembl_db, "compara", "ProteinTree");  
          my $genome_db_adaptor = $registry->get_adaptor($ensembl_db, 'compara', 'GenomeDB');
          my $NCBI_taxon_adaptor = $registry->get_adaptor($ensembl_db, "compara", "NCBITaxon");
          
          my $source_taxon = $NCBI_taxon_adaptor->fetch_node_by_name($sourceorg);
          my $source_NCBI_taxon_ID = $source_taxon->ncbi_taxid;
          
          #We need to verify that the genome exists in the db
          my $all_genome_dbs = $genome_db_adaptor->fetch_all();
          foreach my $genome (@{$all_genome_dbs}){
               $genome_taxon_ids{$genome->taxon_id} = 1;
          }

          unless($genome_taxon_ids{$source_NCBI_taxon_ID}){
               print "Genome name: $sourceorg ($source_NCBI_taxon_ID)  not recognised in ensembl db: $ensembl_db. Skipping this db..\n\n";
               next;
          }
          
          if($destorg eq "All"){#all species available
               my $num_of_genomes = @{$all_genome_dbs};
               print("\n$num_of_genomes genomes considered in database: $ensembl_db. \n") unless($no_output);
          }else{    #one specific species selected
               my $dest_taxon = $NCBI_taxon_adaptor->fetch_node_by_name($destorg);
               my $dest_NCBI_taxon_ID = $dest_taxon->ncbi_taxid;
               
               unless($genome_taxon_ids{$dest_NCBI_taxon_ID}){
                    print "Genome name: $destorg ($dest_NCBI_taxon_ID) is not recognised in ensembl db: $ensembl_db. Skipping this db..\n";
                    next;
               }
               my $method_link_species_set_adaptor = $registry->get_adaptor($ensembl_db, "compara", "MethodLinkSpeciesSet");
               my $gdb1; my $gdb2;

               eval { local $SIG{'__DIE__'}; $gdb1 = $genome_db_adaptor->fetch_by_taxon_id($source_NCBI_taxon_ID); };    warn $@ if $@;
               eval { local $SIG{'__DIE__'}; $gdb2 = $genome_db_adaptor->fetch_by_taxon_id($dest_NCBI_taxon_ID);   };    warn $@ if $@;

               $orthologues_mlss = $method_link_species_set_adaptor->fetch_by_method_link_type_GenomeDBs("ENSEMBL_ORTHOLOGUES",[$gdb1,$gdb2]);
               if(!$orthologues_mlss){
                    print "fetch_by_method_link_type_GenomeDBs for $gdb1, $gdb2 returns undefined. Aborting..\n";
                    return;
               }
          }
          
          while (<$in_data>){
               my ($ID) = $_;
               chomp $ID;
               next if ($ID eq '');
               
               $counter = 0;
               my $gene;my @genes;
               
               #stable id and display label
               $gene = $source_species_gene_adaptor->fetch_by_stable_id($ID);
               if(!$gene){
                    $gene = $source_species_gene_adaptor->fetch_by_display_label($ID);
               }
               if($gene){
                    push(@genes,$gene);
               }else{
                    @genes = @{$source_species_gene_adaptor->fetch_all_by_external_name($ID)};
               }

               foreach my $gene (@genes){
                    my $all_homologies;
                    my $gid = $gene->stable_id;
                    my $member = $member_adaptor->fetch_by_source_stable_id("ENSEMBLGENE", $gid);
               
                    if (defined $member){
                         if($destorg eq "All"){ #all destination genomes
                              $all_homologies = $homology_adaptor->fetch_all_by_Member($member);
                         }else{ #one destination genome
                              $all_homologies = $homology_adaptor->fetch_all_by_Member_MethodLinkSpeciesSet($member, $orthologues_mlss);
                         }
                    }else{
                         print "$gid ..no member object defined in Ensembl, skipping..\n";
                         next;
                    }
                    next if (scalar(@$all_homologies) == 0);
                    
                    print $gid, " ", ($gene->external_name || '-'), "\t" unless($no_output);
                    $counter = _process_homologies(homology_query_id    =>    $gid, 
                                                   homology_vector       =>   $all_homologies, 
                                                   protein_adaptor       =>   $proteintree_adaptor,
                                                   outfile               =>   $out_data,
                                                   hq_only               =>   $onetoone_only
                                                   );
                    print "..$counter orthologue(s).\n" unless($no_output);
                    $counter_db += $counter;
               }
          }
          print "Found $counter_db orthologues in database: $ensembl_db\n" unless($no_output);
          $global_count += $counter_db;
          seek($in_data,0,0);
    }
    close($in_data);
    close($out_data);   
    print "\n**Found $global_count orthologues in all databases**\n" unless($no_output);
    if($global_count == 0){
         unlink($out_path);
         print("No orthologues found for $sourceorg. Exiting..\n");
         return;
    }
    return 1;
}



=head2 get_interactions

 Usage     : $RC = Bio::Homology::InterologWalk::get_interactions(
                                                  input_path     => $in_path,
                                                  output_path    => $out_path,
                                                  url            => $url,
                                                  no_spoke       => 1, 
                                                  exp_only       => 1, 
                                                  physical_only  => 1,
                                                  no_output      => 0 
                                                  );
 Purpose   : this methods allows  to query the Intact database using the REST interface. 
             IntAct is the Molecular Interaction database at the European Bioinformatics 
             Institute (UK). The Intact project offers programmatic access to their data 
             through the PSICQUIC specification 
             (see http://code.google.com/p/psicquic/wiki/PsicquicSpecification).
             This subroutine interrogates via Rest the Intact PPI db with a list of ensembl
             gene ids (obtained usually from get_forward_orthologies()), obtains data in 
             the PSI-MI TAB format (see http://code.google.com/p/psimi/wiki/PsimiTabFormat), 
             processes it and appends it to the input data. 
 Returns   : success/failure code
 Argument  : -input_path : path to input file. Input file for this subroutine is the 
              output of get_forward_orthologies()
             -output_path : where you want the routine to write the data. Data is in TSV 
              format.
             -url : url for the REST service to query (currently only EBI Intact PSICQUIC 
              Rest)
             -(OPTIONAL) no_spoke: if set, interactions obtained from the expansion of 
              complexes through the SPOKE method 
              (see http://nar.oxfordjournals.org/cgi/content/full/38/suppl_1/D525)
              will be ignored
             -(OPTIONAL) exp_only: if set, only interactions whose MITAB25 field "Interaction 
              Detection Method" (MI:0001 in the PSI-MI controlled vocabulary) is at 
              least "experimental interaction detection" 
              (MI:0045 in the PSI-MI controlled vocabulary) will be retained. I.e. if set, 
              this flag only allows experimentally detected interactions to be retained and 
              stored in the data file
             -(OPTIONAL) physical_only: if set, only interactions whose MITAB25 field 
              "Interaction Type" (MI:0190 in the PSI-MI controlled vocabulary) is at least 
              "physical association" 
              (MI:0915 in the PSI-MI controlled vocabulary) will be retained. I.e. if set, 
              this flag only allows physically associated PPIs to be retained and stored 
              in the data file: colocalizations and genetic interactions will be discarded
             -(OPTIONAL) no_output :  suppresses screen output. Used for clearer output 
              during test. Default is 0.
 Throws    : -
 Comment   : -will soon be extended to work with other PSICQUIC-enabled protein interaction 
              dbs (for a list, see 
              http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=STATUS)
             -need to merge with get_direct_interactions. Maybe create core sub, then share.


See Also   : L</get_forward_orthologies>

=cut


sub get_interactions{
     my %args = @_;
     
     my $in_path         = $args{input_path};
     my $out_path        = $args{output_path};
     my $url             = $args{url};
     my $no_spokes       = $args{no_spoke};
     my $exp_only        = $args{exp_only};
     my $physical_only   = $args{physical_only};
     my $no_output       = $args{no_output};
     
     if(!$in_path){
          print("get_direct_interactions(): no PSICQUIC url specified. Aborting..\n");
          return;
     }
     if(!$url){
          print("get_direct_interactions(): no PSICQUIC url specified. Aborting..\n");
          return;
     }
     
     #MANAGE FILES
     my $dbh = _setup_dbi_connection($in_path);
     open (my $out_data,  q{>}, $out_path) or croak("Unable to open $out_path : $!");
     #============
     
     #interactor search string
     my $int_search_string = "search/interactor/";
     #GLOBAL search string
     my $glob_search_string = "search/query/";

     my $atleast_one_entry;
     my $totalSpokes = 0;
     my $client = REST::Client->new();
     
     my $DF_interaction_id;
     my $DF_acc_numb_a;  
     my $DF_acc_numb_b;
     my $DF_alt_id_a;
     my $DF_alt_id_b;
     my $DF_name_a; 
     my $DF_name_b;
     my $DF_taxon_a;
     my $DF_taxon_b;
     my $DF_pub;
     my $DF_int_type;    
     my $DF_det_method;
     my $DF_exp_method;
     
     ##################set up query to input file
     my ($query) = "SELECT * FROM int";
     my $sthHash = $dbh->prepare($query);
     my $sthVec = $dbh->prepare($query);
     $sthHash->execute() or die "Cannot execute: " . $sthHash->errstr();
     $sthVec->execute() or die "Cannot execute: " . $sthVec->errstr();
     ###################
     
     #Header
     print $out_data $HEADER_FWD_ORTH, "\t", $HEADER_INT, "\n";
     my $options = _build_query_options($no_spokes, $exp_only, $physical_only);
     
     while (my $rowHash = $sthHash->fetchrow_hashref){
          my @oldDataVec = $sthVec->fetchrow_array(); 
     
          my $ID = $rowHash->{$FN_initid};
          my $orthologueID = $rowHash->{$FN_orthologue_id};
     
          # do line-by-line processing.
          print "$ID: Querying IntAct WS for $orthologueID.." unless($no_output);
     
          my $queryTerm = $orthologueID;
          my $request = $url . $int_search_string . $queryTerm . $options ;
     
          $client->GET($request);
          print "(", $client->responseCode(), ")" unless($no_output);
  
          my $responseContent = $client->responseContent();
          if(!$responseContent){
               #Let's try
               #a global search, to search for the id in non-standard data fields:
               $request = $url . $glob_search_string . $queryTerm . $options;
               $client->GET($request);
               
               #print "(", $client->responseCode(), ")";
               $responseContent = $client->responseContent();
               if(!$responseContent){
                    print("..nothing..\n") unless($no_output);
                    next;
               }
          }
          $atleast_one_entry = 1;
          my @responsetoparse = split(/\n/,$responseContent);
          my $interactionsRetrieved = scalar @responsetoparse;
          print "..Interactions found: ", $interactionsRetrieved, "\n" unless($no_output);
     
          foreach my $intactInteraction (@responsetoparse){
               my @MITABDataRow = split("\t",$intactInteraction);
               
               #Field values rely on MITAB specification: if it changes, this won't work anymore
               $DF_interaction_id  = _get_intact_id($MITABDataRow[13]);
               next if(!$DF_interaction_id);
               $DF_acc_numb_a      = _get_interactor_uniprot_id($MITABDataRow[0]);
               $DF_acc_numb_b      = _get_interactor_uniprot_id($MITABDataRow[1]);
               $DF_alt_id_a        = _get_interactor_aliases($MITABDataRow[2]);
               $DF_alt_id_b        = _get_interactor_aliases($MITABDataRow[3]);
               $DF_name_a          = _get_interactor_name($MITABDataRow[4]);
               $DF_name_b          = _get_interactor_name($MITABDataRow[5]);
               $DF_taxon_a         = _get_interactor_taxon($MITABDataRow[9]);
               $DF_taxon_b         = _get_interactor_taxon($MITABDataRow[10]);
               $DF_pub             = $MITABDataRow[8];
               $DF_int_type        = $MITABDataRow[11];
               $DF_det_method      = $MITABDataRow[6];
               $DF_exp_method      = $MITABDataRow[24];
          
               print("Interaction ($DF_interaction_id): $DF_acc_numb_a <-> $DF_acc_numb_b\n") unless($no_output);
          
               my $fullDataRow = join("\t",@oldDataVec,$DF_interaction_id,$DF_acc_numb_a,$DF_acc_numb_b,$DF_alt_id_a,$DF_alt_id_b,
                                             $DF_name_a,$DF_name_b,$DF_taxon_a,$DF_taxon_b,$DF_pub,$DF_int_type,$DF_det_method,$DF_exp_method);
               print $out_data $fullDataRow, "\n";
          }
     }
     $sthVec->finish();
     $sthHash->finish();
     $dbh->disconnect();
     close($out_data);
     if(!$atleast_one_entry){
        unlink($out_path);
        print("\n**No interactions found. Exiting..**\n");
        return;
     }
     return 1;
}


=head2 get_backward_orthologies

 Usage     : $RC = Bio::Homology::InterologWalk::get_backward_orthologies(
                                                          registry      => $registry,
                                                          ensembl_db    => $ensembl_db,
                                                          input_path    => $in_path,
                                                          output_path   => $out_path,
                                                          error_path    => $err_path,
                                                          source_org    => $sourceorg,    
                                                          hq_only       => $onetoone,
                                                          no_output     => 0
                                                          );
 Purpose   : this routine mines orthologues back into the organism of interest. It accepts 
             as an input a data file containing interactions in the destination organism(s) 
             and maps those back to the source organism through orthology. 
             Such orthologues represent the putative interactors of the original genes as 
             requested.
 Returns   : success/error
 Argument  : -registry: registry object for ensembl connection
             -ensembl db to connect to. Choices currently are:
                 a. 'multi' :  vertebrate compara (see http://www.ensembl.org/)
                 b. 'pan_homology' : pan taxonomic compara db, a selection of species from 
                    both Ensembl Compara and Ensembl Genomes 
                    (see http://nar.oxfordjournals.org/cgi/content/full/38/suppl_1/D563 )
                 c. 'metazoa' : ensembl compara genomes, metazoa db 
                    (see http://metazoa.ensembl.org/index.html). 
                 d. 'all'  : multi + metazoa.
                 Default is 'multi'.
             -input_path : path to input file. Input file for this subroutine is the output 
              of get_interactions().
             -output_path : where you want the routine to write the data. Data is in TSV 
              format.
             -(OPTIONAL)error_path: each query to intact through psicquic returns a data 
               entry including a binary protein interaction. The two ids returned are, most 
               of the times, uniprotkb ids. Sometimes, however, Intact annotates its binary 
               interactions using an internal, proprietary ID (e.g.: EBI-1080281 ). While the 
               Ensembl API recognises UniprotKB IDs,it won't recognise these Intact IDs. Entries 
               annotated in such a way cannot therefore be completed. If error_path is present, 
               it indicates a file where the routine will dump all such failed entries for later 
               manual inspection.
             -source organism name (eg: "Mus musculus")
             -(OPTIONAL)hq_only: discards one-to-many, many-to-one, many-to-many orthologues. 
              Only keeps one-to-one orthologues, i.e. where no duplication event has happened 
              after the speciation. One-to-one orthologues are ideally associated with higher 
              functional conservation (while paralogues often cause neo/sub-functionalisation). 
              For further information see
              http://www.ensembl.org/info/docs/compara/homology_method.html 
             -(OPTIONAL) no_output :  suppresses screen output. Used for clearer output during 
              test. Default is 0.
 Throws    : -
 Comment   : Destination species is automatically dealt with on a case-to-case basis.
           : 'ensembl_db' must be the same  for all the other subroutines in the pipeline

See Also   : L</get_interactions>

=cut


sub get_backward_orthologies{
     my %args = @_;
     my $registry        = $args{registry};
     my $db              = $args{ensembl_db};
     my $sourceorg       = $args{source_org};
     my $in_path         = $args{input_path};
     my $out_path        = $args{output_path};
     my $err_path        = $args{error_path};
     my $onetoone_only   = $args{hq_only};
     my $no_output       = $args{no_output};
     
     if(!$registry){
       print("get_backward_orthologies(): no registry object found. Aborting..\n");
       return;
     }
     if(!$sourceorg){
      print("get_backward_orthologies(): no source organism specified. Aborting..\n");
       return;
     }
     $db = "multi" if (!$db);

     $sourceorg = ucfirst(lc $sourceorg);
     
     my $out_data; 
     my $err_data;
     #MANAGE FILES
     my $dbh = _setup_dbi_connection($in_path);
     open ($out_data,  q{>}, $out_path) or croak("Unable to open $out_path : $!");
     if($err_path){
        open ($err_data,  q{>}, $err_path) or croak("Unable to open $err_path : $!");
     }
     #============
     
     my $entries_err_found;
     
     $ENSEMBLIDFAILED = undef; 
     my $failedcounter = 0;
     my $global_count = 0;
     my $orthology_count = 0;
     my @ensembl_dbs;
     
     my @NCBItaxa = ();  
     my %NCBItaxa_seen = ();
     my $unique_taxa = 0;
     my %NCBI_taxa_global = ();
     
     if($db eq "all"){
      #the order of elements in this vector MUST BE in sync with the order of calling the dbs in the registry routine
          push(@ensembl_dbs, "multi");
          push(@ensembl_dbs, "metazoa");
     }else{
          push(@ensembl_dbs,$db);
     }

     print $err_data $HEADER_FWD_ORTH, "\t", $HEADER_INT, "\n" if($err_path);
     print $out_data $HEADER_FWD_ORTH, "\t", $HEADER_INT, "\t", $HEADER_BWD_ORTH, "\n";
     #collect unique taxa and basic statistics from the input file===============
     my ($taxonquery) = "SELECT $FN_taxon_a, $FN_taxon_b
                              FROM int";
                         
     my $sth = $dbh->prepare($taxonquery);
     $sth->execute() or die "Cannot execute: " . $sth->errstr();
     
     my $rowCount = $dbh->selectrow_array(
                                 qq{
                                      SELECT count(*)
                                         FROM int
                                });
     print "Total number of interactions: $rowCount.\n" unless($no_output);
     
     my $DF_taxon_a;
     my $DF_taxon_b;
     while (my $row = $sth->fetchrow_hashref) {
          $DF_taxon_a = $row->{$FN_taxon_a};
          $DF_taxon_b = $row->{$FN_taxon_b};
          #I don't want empty fields or hyphens or other crap
          $NCBItaxa_seen{$DF_taxon_a} = 1 if($DF_taxon_a =~ /\d+/);
          $NCBItaxa_seen{$DF_taxon_b} = 1 if($DF_taxon_b =~ /\d+/);
     }
     
     @NCBItaxa = sort(keys(%NCBItaxa_seen));
     $unique_taxa += 1 foreach (@NCBItaxa);
     print "Total Number of Unique NCBI taxa IDs: $unique_taxa.\n" unless($no_output);
     $sth->finish();
     #=============================================================================
     foreach my $ensembl_db (@ensembl_dbs){ # max 2 at the mom
          my %genome_taxon_ids;
          my $all_genome_dbs;
          
          print("\n----Querying Ensembl Compara ($ensembl_db) for orthologues back in $sourceorg----\n") unless($no_output);
          
          #CHECK THE GENOME ID/DB COMBO====
          #first we get the id from the source organism name
          my $NCBI_taxon_adaptor = $registry->get_adaptor($ensembl_db, "compara", "NCBITaxon");
          my $source_taxon = $NCBI_taxon_adaptor->fetch_node_by_name($sourceorg);
          my $source_NCBI_taxon_ID = $source_taxon->ncbi_taxid;
          
          my $genomedb_adaptor = $registry->get_adaptor($ensembl_db, "compara", "GenomeDB");  
          $all_genome_dbs = $genomedb_adaptor->fetch_all();
          foreach my $genome (@{$all_genome_dbs}){
               $genome_taxon_ids{$genome->taxon_id} = 1;
          }
          #If the given species does not exist in that db, the db is  not worth continuing
          unless($genome_taxon_ids{$source_NCBI_taxon_ID}){ 
               print "Genome name: $sourceorg ($source_NCBI_taxon_ID) not recognised in ensembl db: $ensembl_db. Skipping this db..\n";
               next;
          }
          #destination species must be the same as start species in the forward orthology script
          my $gdb2;
          eval { local $SIG{'__DIE__'}; $gdb2 = $genomedb_adaptor->fetch_by_taxon_id($source_NCBI_taxon_ID); };     warn $@ if $@;
          #=============================
          
          #SET UP FIXED ENSEMBLE VARIABLES====================================================
          my $member_adaptor =
          $registry->get_adaptor($ensembl_db, "compara", "Member");
          my $homology_adaptor =
          $registry->get_adaptor($ensembl_db, "compara", "Homology");
          my $method_link_species_set_adaptor =
          $registry->get_adaptor($ensembl_db, "compara", "MethodLinkSpeciesSet");
          my $proteintree_adaptor = 
               $registry->get_adaptor($ensembl_db, "compara", "ProteinTree");
          #===================================================================================
                    
          foreach my $NCBItaxonID (@NCBItaxa){
               next if($NCBI_taxa_global{$NCBItaxonID});
               $NCBI_taxa_global{$NCBItaxonID} = 1;
               #If I obtained the data for this taxon from the previous db, I skip it
               
               #query===========================
               my ($query) = "SELECT * 
                                   FROM int
                                   WHERE ($FN_taxon_a='$NCBItaxonID' AND $FN_taxon_b = $FN_taxon_a)
                                   OR ($FN_taxon_a='$NCBItaxonID' AND $FN_taxon_b = '-1')
                                   OR ($FN_taxon_a='$NCBItaxonID' AND $FN_taxon_b = '-2')
                                   OR ($FN_taxon_b='$NCBItaxonID' AND $FN_taxon_a = '-1')
                                   OR ($FN_taxon_b='$NCBItaxonID' AND $FN_taxon_a = '-2')";
               #intact uses numerical codes instead of NCBI taxon ids sometimes. To my knowledge, these are -1 and -2
               #('in vitro' and 'chemical synthesis')
               #I added the previous 3 lines cause I want such cases to be kept
               my $sthHash = $dbh->prepare($query);
               my $sthVec = $dbh->prepare($query);
               $sthHash->execute() or die "Cannot execute: " . $sthHash->errstr();
               $sthVec->execute() or die "Cannot execute: " . $sthVec->errstr();
               #===========================================================
     
               next if (($sthHash->rows) == 0);
               #SET UP SPECIES-DEPENDENT ENSEMBLE VARIABLES==============================
               my $NCBItaxon = $NCBI_taxon_adaptor->fetch_node_by_taxon_id($NCBItaxonID);
               next if (!$NCBItaxon);
               my $NCBItaxon_scientific_name = $NCBItaxon->binomial;
               print "\n============\n" unless($no_output);
               print "$NCBItaxon_scientific_name ($NCBItaxonID)\n" unless($no_output);
               print "============\n" unless($no_output);
          
               my $gene_adaptor = $registry->get_adaptor("$NCBItaxon_scientific_name", "core", "Gene");
               if(!$gene_adaptor){
                    print "gene_adaptor for $NCBItaxon_scientific_name not defined in Bio::EnsEMBL::Registry\n";
                    next;
               }
               
               if(!$genome_taxon_ids{$NCBItaxonID}){
                    print "Genome name: $NCBItaxon_scientific_name ($NCBItaxonID) not recognised in ensembl db: $ensembl_db. Skipping this taxon..\n";
                    next;
               }

               my $gdb1;
               eval { local $SIG{'__DIE__'}; $gdb1 = $genomedb_adaptor->fetch_by_taxon_id($NCBItaxonID); };  warn $@ if $@;
               
               my $orthologues_mlss = 
               $method_link_species_set_adaptor->fetch_by_method_link_type_GenomeDBs("ENSEMBL_ORTHOLOGUES",[$gdb1,$gdb2]);
               if(!$method_link_species_set_adaptor){
               print "fetch_by_method_link_type_GenomeDBs for $gdb1, $gdb2 returns undefined\n";
               next;
               }
               #==========================================================================
               while (my $row = $sthHash->fetchrow_hashref) {
                    my @oldDataVec = $sthVec->fetchrow_array(); #take the next row in a vector
                    my $entry = join("\t",@oldDataVec);
                    
                    my $candidate; my $member;
                    my $interactoridA; my $interactoridB;
                    
                    my $DF_orthologue_id   = $row->{$FN_orthologue_id};
                    my $DF_interaction_id  = $row->{$FN_interaction_id};
                    my $DF_name_a          = $row->{$FN_name_a}; 
                    my $DF_name_b          = $row->{$FN_name_b};
                    my $DF_alt_id_a        = $row->{$FN_alt_id_a};
                    my $DF_alt_id_b        = $row->{$FN_alt_id_b}; 
                    my $DF_acc_numb_a      = $row->{$FN_acc_numb_a};
                    my $DF_acc_numb_b      = $row->{$FN_acc_numb_b};
                    
                    #uniprotkb->ensemble stable id translation
                    if($DF_acc_numb_a =~ /^ENS/){
                         $interactoridA = $DF_acc_numb_a;
                    }else{
                         $interactoridA = _get_ensembl_id(adaptor   => $gene_adaptor, 
                                                     ebi_id         => $DF_interaction_id, 
                                                     ortho_id       => $DF_orthologue_id, 
                                                     acc_numb       => $DF_acc_numb_a, 
                                                     protein_name   => $DF_name_a, 
                                                     aliases        => $DF_alt_id_a,
                                                     no_output      => $no_output);
                    }
                    
                    if($DF_acc_numb_b =~ /^ENS/){
                         $interactoridB = $DF_acc_numb_b;
                    }else{
                         $interactoridB = _get_ensembl_id(adaptor    => $gene_adaptor, 
                                                        ebi_id       => $DF_interaction_id, 
                                                        ortho_id     => $DF_orthologue_id, 
                                                        acc_numb     => $DF_acc_numb_b, 
                                                        protein_name => $DF_name_b, 
                                                        aliases      => $DF_alt_id_b,
                                                        no_output    => $no_output);
                    }

               
                    if($ENSEMBLIDFAILED){#this is only to test if, for either or both of the two runs of getensembl id, i had to resort to
                    #using name, or aliases to try and get the ensembl id. this is only for benchmarking my back-up technique
                         $failedcounter += 1;
                         $ENSEMBLIDFAILED = undef;
                    }
                    
                    if($interactoridA && $interactoridB ){
                         if($DF_orthologue_id eq $interactoridA){
                              $candidate = $interactoridB;
                              $member = $member_adaptor->fetch_by_source_stable_id("ENSEMBLGENE", $candidate);
                         }elsif($DF_orthologue_id eq $interactoridB){
                              $candidate = $interactoridA;
                              $member = $member_adaptor->fetch_by_source_stable_id("ENSEMBLGENE", $candidate);
                         }else{
                              print("\nWARNING id mismatch between $DF_orthologue_id, $interactoridA, $interactoridB.\n");
                         }
                    }else{
                         #you might want to pack this into a sub
                         print("Converting all IDs in Uniprot KB IDs..\n") unless($no_output);
                         
                         $candidate = _compare_uniprotkbids($gene_adaptor, $DF_orthologue_id, $DF_acc_numb_a, $DF_acc_numb_b);
                         if(!$candidate){
                             if($err_path){
                                print $err_data $entry, "\n";
                                 $entries_err_found = 1;
                             }
                             next;
                         }
                              
                         $member = $member_adaptor->fetch_by_source_stable_id("Uniprot/SWISSPROT", $candidate);
                         if(!$member){
                              $member = $member_adaptor->fetch_by_source_stable_id("Uniprot/SPTREMBL", $candidate);
                         }
                    }
                    
                    if(!$member){
                         print "$DF_interaction_id-($candidate): member_adaptor->fetch_by_source_stable_id returns no member object.\n" unless($no_output);
                         if($err_path){
                            print $err_data $entry, "\n";
                            $entries_err_found = 1;
                         }
                         next;     
                    }
                    
                    my $all_homologies = $homology_adaptor->fetch_all_by_Member_MethodLinkSpeciesSet($member, $orthologues_mlss);
                    next if (scalar(@$all_homologies) == 0);
                    
                    $orthology_count = _process_homologies(homology_query_id   =>   $candidate, 
                                                           homology_vector     =>   $all_homologies, 
                                                           protein_adaptor     =>   $proteintree_adaptor,
                                                           outfile             =>   $out_data,
                                                           datavector          =>   $entry, 
                                                           hq_only             =>   $onetoone_only,
                                                           no_output           =>   $no_output
                                                           );
                    $global_count += $orthology_count;
               }
               $sthHash->finish();
               $sthVec->finish();
          }
     }
     $dbh->disconnect();
     close $out_data;
     if($err_path){
          close $err_data;
          unlink($err_path) if(!$entries_err_found);
     }
     if($global_count == 0){
          unlink($out_path);
          print("No orthologues found back in $sourceorg. Exiting..\n");
          return;
     }
     return 1;
}



=head2 do_counts

 Usage     : $RC = Bio::Homology::InterologWalk::do_counts(
                                           input_path  => $in_path,
                                           output_path => $out_path,
                                           header      => 'standard',
                                           no_output   => 0
                                           );
 Purpose   : The purpose of this routine is to scan the data produced by get_backward_orthologies() 
             or get_direct_interactions() (optionally cleaned up of duplicates by 
             remove_duplicate_rows() ) and compute counts/statistics useful for scoring purposes.
             In short, the subroutine:
             1)evaluates if an interaction has been obtained through more than one detection method
             2)evaluates if an interaction has been obtained through more than one taxon
             3)COUNTS the number of *unique* putative interactions found: remember that the same 
               interaction can be retrieved through several different interacting destination-species 
               orthologues. This script also adds the retrieved "number seen" number and appends it 
               to the TSV file. 
             4)flags the entry with Y if the putative interaction is an autointeraction
             5)flags the entry if the real interaction in the destination species (the one we are 
               mapping from) is an autointeraction
             The routine rewrites the input file in a new file, adding 1 or more data fields 
             (depending on the 'header' argument) containing the results of the count.
 Returns    : success/fail
 Argument   : -input_path : path to input file. Input file for this subroutine is a TSV file of PPIs. 
              It can be one of the following two:
                 1. the output of get_backward_orthologies(). In this case please specify 'standard' 
                    header below.
                 2. the output of get_direct_interactions(). In this case please specify 'direct' 
                    header below.   
              It is advisable to pre-process the input by using remove_duplicate_rows() prior to 
              this routine.          
             -output_path : where you want the routine to write the data. Data is in TSV format. 
             -(OPTIONAL) header : Header type is one of the following:
                 1. 'standard': when the routine is used to compute counts on an interolog-walk file 
                    (the header will be longer)
                 2. 'direct':   when the routine is used to compute counts on a real db interactions 
                    file (the header is shorter)
               No field provided: default is 'standard'
             -(OPTIONAL) no_output :  suppresses screen output. Used for clearer output during test. 
              Default is 0.
 Throws    : -
 Comment   : -

See Also   : L</get_backward_orthologies>, L</get_direct_interactions>, L</remove_duplicate_rows>

=cut


sub do_counts{
     my %args = @_;
     my $in_path       = $args{input_path};
     my $out_path      = $args{output_path};
     my $header_type   = $args{header};
     my $no_output     = $args{no_output};
     
     $header_type = 'standard' if(!$header_type);
     #MANAGE FILES
     my $dbh = _setup_dbi_connection($in_path);
     open (my $out_data,  q{>}, $out_path) or croak("Unable to open $out_path : $!");
     #============
     
     my %idSeen = ();
     my %intactInteractionSeen = ();#I want to know: has this Intact interaction been confirmed through other detection methods?
     my %seenThroughTaxa = ();#I want to know: has this putative interaction been obtained across different taxa?
     
     #header====
     if($header_type eq 'standard'){
          print $out_data $HEADER_FWD_ORTH, "\t", $HEADER_INT, "\t", $HEADER_BWD_ORTH, "\t", $HEADER_MISC,  "\n";
     }elsif($header_type eq 'direct'){
          print $out_data $HEADER_DIRECT, "\t", $FN_multiple_dm,  "\n";
     }else{
          print("do_counts(): header: $header_type string not recognised, aborting..");
          return;
     }
     #==================================
     
     my ($query) = "SELECT * FROM int";
     my $rowCount = $dbh->selectrow_array(
                                  qq{
                                        SELECT count(*)
                                        FROM int
                                    });
     my $sthHash = $dbh->prepare($query);
     my $sthVec = $dbh->prepare($query);
     $sthHash->execute() or die "Cannot execute: " . $sthHash->errstr();
     #====================================
     
     #======================
     #Gathering Statistics
     #======================
     while (my $row = $sthHash->fetchrow_hashref) {
          my $id_in       = $row->{$FN_initid};
          my $id_out      = $row->{$FN_interactor_id};
          my $accA        = $row->{$FN_acc_numb_a};
          my $accB        = $row->{$FN_acc_numb_b};
          my $detMethod   = $row->{$FN_det_method};
          my $taxon       = $row->{$FN_taxon_a};
          #my $taxonB = $row->{$FN_taxon_b};
          #I actually filter them to be equal in get_backward_orthologies so taxona = taxonb for my dataset
     
          ############################
          # Checking for reinforcement
          # of interaction through several det methods
          ############################
          my $intactPair = join("-", sort($accA, $accB));
          #hash of hash
          $intactInteractionSeen{$intactPair}{$detMethod} = 1;
          ############################
          #Checking for reinforcement 
          #of interaction in multiple taxa
          ############################
          my $idpair = join('-', sort($id_in, $id_out));
          #hash of hash
          $seenThroughTaxa{$idpair}{$taxon} = 1;
          
          if(exists($idSeen{$idpair})){
               $idSeen{$idpair} += 1;
          }else{
               $idSeen{$idpair} =  1;
          }
     }
     my $idNumber = keys %idSeen;
     #foreach my $interaction (sort keys %idSeen){
     #    print $interaction, " : ", $idSeen{$interaction}, "\n";
     #}
     
     print("\nTotal number of interactions (including duplicates): $rowCount\n") unless($no_output);
     print("Total number of interactions involving UNIQUE (initial_id, final_id) pairs: $idNumber\n") unless($no_output);
          
     #======================
     #Printing new data file
     #======================
     $sthHash->execute() or die "Cannot execute: " . $sthHash->errstr();
     $sthVec->execute() or die "Cannot execute: " . $sthVec->errstr();
     while (my $rowHash = $sthHash->fetchrow_hashref){
          my @rowVec = $sthVec->fetchrow_array();
          
          my $idpair = ""; my $accpair = "";
          my $numberseen = 0;
          my $sameid = "-"; my $sameaccnumb = "-";
          my $reconfirmed_DM;
          my $reconfirmed_TAXA;
          
          my $id_in = $rowHash->{$FN_initid};
          my $id_out = $rowHash->{$FN_interactor_id};
          my $acc_A = $rowHash->{$FN_acc_numb_a};
          my $acc_B = $rowHash->{$FN_acc_numb_b};
          
          $sameid = "Y" if($id_in eq $id_out);
          $sameaccnumb = "Y" if($acc_A eq $acc_B);
          
          $idpair = join('-', sort($id_in, $id_out));
          $numberseen = $idSeen{$idpair};
#         if($numberseen eq 0){
#              print "\nERROR: ", $idpair, "\n";
#              exit 1;
#         }
          $accpair = join("-", sort($acc_A, $acc_B));
          for my $detMethod ( keys %{ $intactInteractionSeen{$accpair} } ) {
               $reconfirmed_DM += 1;
          }
          for my $taxon ( keys %{ $seenThroughTaxa{$idpair} } ){
              $reconfirmed_TAXA +=1;
          }
         
          my $newRow;
          #to be improved upon. If I'm in the second case I'm computing 4 values I'm not using. I should really be making a separate routine
          if($header_type eq 'standard'){
               $newRow = join("\t", @rowVec, $reconfirmed_DM,$reconfirmed_TAXA, $numberseen, $sameid, $sameaccnumb);
          }elsif($header_type eq 'direct'){
               $newRow = join("\t", @rowVec, $reconfirmed_DM);
          }
          print $out_data $newRow, "\n";
     }
     $sthHash->finish();
     $sthVec->finish();
     $dbh->disconnect();
     close($out_data);
     return 1;
}


=head2 extract_unseen_ids

 Usage     : $RC = Bio::Homology::InterologWalk::extract_unseen_ids(
                                                    start_path    => $start_data_path,
                                                    input_path    => $in_path,
                                                    output_path   => $out_path,
                                                    hq_only       => $onetoone,
                                                    );
 Purpose   : it is often desirable to know if the interolog procedure found new ids at all 
             (i.e. not present in the starting dataset). Such new ids can then be analysed 
             further, ie. sent through GO term enrichment analysis, etc, to provide some 
             validation, see if they have been know before to belong to some specific process, 
             check if no function is associated to them at all.
             This script will create a simple textfile containing all the new ids discovered.
             This script is meant to be employed as a last step in the pipeline. It also 
             computes some simple statistics as follows:
             1. The list of NEW ids, ie those not present in the initial data file
             2. The frequencies, new vs total, old vs total
             3. the frequencies of new when the Expansion Method is not spoke and when orthology 
             is one_to_one (i.e.: new ids with high reliability)
 Returns   : Success/Fail
 Argument  : -start_path: path to the original text file with the ids of interest (the same file 
              given to get_forward_orthologies() as input)
             -input_path : path to input file. Input file for this subroutine is the output of 
              do_counts().
             -output_path : where you want the routine to write the data. Data is in TSV format.
             -(OPTIONAL)hq_only : if this is set, only entries mapped exclusively through 
              one-to-one orthologies will be taken into account.
 Throws    : -
 Comment   : -

See Also   : L</get_forward_orthologies>, L</do_counts>

=cut


sub extract_unseen_ids{
     my %args = @_;
     my $start_data_path   = $args{start_path};
     my $in_path           = $args{input_path};
     my $out_path          = $args{output_path};
     my $onetoone_only     = $args{hq_only};
     
     #MANAGE FILES
     my $dbh = _setup_dbi_connection($in_path);
     open (my $start_data, q{<}, $start_data_path) or croak("Unable to open $start_data_path : $!");
     open (my $out_data, q{>}, $out_path) or croak("Unable to open $out_path : $!");
     #============

     my $query;
     my %old_id_present = ();
     my %new_id_set = ();
     my %all_seen = (); #total number of nodes
     
     #we slurp all the ids in the initial file into a hash========
     my %start_data_set;
     while(<$start_data>) {
          chomp;
          next if($_ eq '');
          $start_data_set{$_} = 1;
     }
     my $number_of_elements_start_ds = keys %start_data_set;
     print("Number of unique ids in original dataset: $number_of_elements_start_ds\n");
     #=========================================================
     
     #now we get both interactors from the putative data set and we look at those that dont appear in the %start_data_set hash
     if($onetoone_only){
          $query = "SELECT    $FN_initid, $FN_interactor_id
                    FROM int
                    WHERE ($FN_odesc_1 like '%one2one%') 
                    AND ($FN_odesc_2 like '%one2one%')";
     }else{ # all orthology types
          $query = "SELECT    $FN_initid, $FN_interactor_id
                    FROM int";
     }
     
     my $sth = $dbh->prepare($query);
     $sth->execute() or die "Cannot execute: " . $sth->errstr();
     
     my $rownumb = $sth->rows;             
     return if ($rownumb == 0); 
     #you need to do this with all subroutines
     
     while (my $row = $sth->fetchrow_hashref) {
          my $idIN = $row->{$FN_initid};
          my $idOUT = $row->{$FN_interactor_id};
     
          #first I count all of the unique ids present (ie network nodes)
          $all_seen{$idIN} = 1;    $all_seen{$idOUT} = 1;
          
          #then I count what fraction of the original ids actually made it to feature any putative interactions
          $old_id_present{$idIN} = 1 if(exists($start_data_set{$idIN}));
          $old_id_present{$idOUT} = 1 if(exists($start_data_set{$idOUT}));
     
          #lastly I store those new ids never seen in the starting dataset
          $new_id_set{$idOUT} = 1 unless(exists($start_data_set{$idOUT})); 
     }
     my $number_of_old_IDs = keys %old_id_present;
     my $percentage = ($number_of_old_IDs / $number_of_elements_start_ds) * 100;
     print("Number of IDs from the original dataset that appear in the network: $number_of_old_IDs\n");
     print("Percentage of IDs from the original dataset that appear in the final dataset: $percentage\n");
     
     my $number_of_new_IDs = keys %new_id_set;
     my $number_of_network_nodes = keys %all_seen;
     $percentage= ($number_of_new_IDs / $number_of_network_nodes) * 100;
     if($onetoone_only){
          print("Number of total UNIQUE IDs in interaction dataset (considering ONE-TO-ONE ortologies only): $number_of_network_nodes\n");
     }else{
          print("Number of total UNIQUE IDs in interaction dataset: $number_of_network_nodes\n");
     }
     print("Number of NEW ids (e.g. not seen in starting data set): $number_of_new_IDs\n");
     print("Percentage of new ids over the total: $percentage\n");
     
     
     #I save all the new ids in a flat file. This might be useful to do some analysis of their functional annotation
     foreach my $id (sort keys %new_id_set){
          #print ("$id\t$new_id_set{$id}\n");
          print $out_data $id, "\n";
     }
     $sth->finish();
     $dbh->disconnect();
     close($start_data);
     close($out_data);
     return 1;
}



####################################
#HELPER routines for internal usage
####################################

#
#_compare_uniprotkbids
#
# Usage     : How to use this function/method
# Purpose   : Suppose we query Intact with the Ensembl ID of a gene of interest, to look for interactors. Intact will not return a simple
#             list of interactors as desidered. It will instead return 0 or more lines of tsv data, each one describing a binary interaction, in MITAB format. 
#             Each lines will contain 2 ids, in uniprot format: one of them being our query id. As the returned format is Uniprot, while the query id was Ensembl,
#             we cannot possibly know which of the two is the query id and which is the interactor. Therefore we must carry out a conversion.
#             This routine looks for a uniprot id of the query id in Ensembl format. It then compares it to both retrieved id to decide which one is the output interactor
#             This is used by get_backward_orthologies().
# Returns   : a candidate id in uniprot format or UNDEF
# Argument  : gene adaptor, Ensembl Id of the query gene, the two uniprotKB accession number returned by Intact
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 
sub _compare_uniprotkbids{
     my ($adaptor, $ensembl_id, $uniprot_id_a, $uniprot_id_b) = @_;
     my %uniprotANseen = ();
     my $output_id;
     
     my $query_gene = $adaptor->fetch_by_stable_id($ensembl_id);
     my $uniprot_links = $query_gene->get_all_DBLinks("Uniprot%");
     
     foreach my $link (@$uniprot_links) {
          my $item = $link->primary_id;
     $uniprotANseen{$item} = 1;
     }
     
     if(exists $uniprotANseen{$uniprot_id_a}){
          $output_id = $uniprot_id_b;
     }elsif(exists $uniprotANseen{$uniprot_id_b}){
          $output_id = $uniprot_id_a;
     }else{
          print("WARNING uniprot-ids of $ensembl_id contain neither $uniprot_id_a nor $uniprot_id_b. Skipping..\n");
     }
     return $output_id;
}

#
#_build_query_options
#
# Usage     : How to use this function/method
# Purpose   : this is used to build an options string to append to the url when doing the REST query
# Returns   : an options string (may be empty)
# Argument  : three boolean values (1/undef): whether we want interactions coming from spoke-expanded complexes, whether we want physical interactions only 
#             or not, whether we want experimentally obtained interactions only or not.
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 
sub _build_query_options{
     my ($no_spoke, $exponly, $physonly) = @_;
     my $outputstring;
     
     #eg " AND type:\"physical*\" AND detmethod:\"experimental*\" AND NOT expansion:\"spoke\"";
     if($no_spoke || $exponly || $physonly){
          my @outputitems;
          
          push(@outputitems, "type:\"physical*\"") if($physonly);
          push(@outputitems, "detmethod:\"experimental*\"") if($exponly);
          push(@outputitems, "NOT expansion:\"spoke\"") if($no_spoke);
          
          $outputstring = join(" AND ", @outputitems);
          $outputstring = " AND " . $outputstring;
          
          return $outputstring;
     }else{
          $outputstring = '';
          return $outputstring;
     }
}

#
#_get_intact_id
#
# Usage     : How to use this function/method
# Purpose   : This is used by get_interactions() to obtain the intact code for the interaction
# Returns   : Id in the form EBI-545553
# Argument  : full id, e.g. intact:EBI-434443
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 
sub _get_intact_id{
     my ($id) = @_;
     
     return if(!$id);
     
     if($id =~ /^(intact):(EBI-\d+)/){
          return $2;     
     }else{
          print "_get_intact_id(): unrecognised entry format: '$id'\n";
          return;
     }
}

#
#_get_interactor_name
#
# Usage     : How to use this function/method
# Purpose   : This is used by get_interactions() to clean up the name field retrieved from the db
#                eg uniprotkb:alp-1 ----> alp-1
# Returns   : The name of the protein
# Argument  : string in DB:name format
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 
sub _get_interactor_name{
     my ($nameString) = @_;
     my @namesVec;
     #no name: exit
     if($nameString eq "-"){
          return $nameString;
     }
     #multiple names:only get the first
     if($nameString =~ /\|+/){
          @namesVec = split(/\|/, $nameString);
          $nameString = $namesVec[0]; #I IGNORE all but the first name when there are several
          return $nameString;
     }
     
     if($nameString =~ /^(uniprotkb):(.*)/){ 
          return $2;
     }else{
          print "_get_interactor_name(): unrecognised entry format: '$nameString. Leaving as it is\n";
          return $nameString;
     }
}

#
#_get_interactor_uniprot_id
#
# Usage     : How to use this function/method
# Purpose   : This is used by get_interactions() to clean up the accession number field retrieved from the db
#                eg uniprotkb:Q9BXW9-2|intact:EBI-596878 --> Q9BXW9-2
# Returns   : The name of the protein
# Argument  : string in DB:name format
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 
sub _get_interactor_uniprot_id{
     my ($id) = @_;
     my $value;
     #ids in order of preference
     #uniprot accession number preferred to intact cause it's recognised by ensembl
     if($id =~ /(ensembl):(\w+)\|(.*)/){ #es ensembl:ENSG00000122592
          $value = $2;
     }elsif ($id =~ /(uniprotkb):(\w+)\|(.*)/){ #es uniprotkb:Q9BXW9
          $value = $2;
     }elsif($id =~ /(uniprotkb):(\w+-\d+)\|(.*)/){#es uniprotkb:Q9BXW9-2|intact:EBI-596878
          $value = $2;
     }elsif($id =~ /(intact):(EBI-\d+)$/){ #if there's no uniprot, get intact-ebi id
          $value = $2;
     }else{
          print "_get_interactor_uniprot_id(): unrecognised entry format: '$id'.\n";
          $value = "-";
     }
     return $value;
}

#
#_get_interactor_aliases
#
# Usage     : How to use this function/method
# Purpose   : This is used by get_interactions() to clean up the accession number field retrieved from the db
#                eg uniprotkb:q0pby3_camje(shortlabel) -> uniprotkb:q0pby3_camje
# Returns   : string in DB:name format
# Argument  : string in DB:name format(other)
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 
sub _get_interactor_aliases{
     #I'm fine with this string, save for the recently introduced explanations in brackets. I
     #want to remove everything in brackets
     my ($aliasString) = @_;
     my $cleanAlias;
     my @cleanAliasVec;
     my $result;
     
     my @aliasvec = split(/\|/, $aliasString);
     
     foreach my $alias (@aliasvec){ #db:string
          if($alias =~ /(\w+):(.+)\((.+)\)/){ #eg.: uniprotkb:q0pby3_camje(shortlabel)
               $cleanAlias = $1 . ":" . $2;
               push(@cleanAliasVec, $cleanAlias);
          }else{#just give up and save the original
               push(@cleanAliasVec, $alias);
          }
     }
     $result = join('|', @cleanAliasVec);
     return $result;
}

#
#_get_interactor_taxon
#
# Usage     : How to use this function/method
# Purpose   : This is used by get_interactions() to clean up the taxon field retrieved from the db geting only the NCBI code
#                eg taxid:6239(Caenorhabditis elegans) ----> 6239
# Returns   : NCBI code for taxon
# Argument  : string in taxid:NCBICODE(Scientific name)
# Throws    : -
# Comment   : you might want to change this so that it also shows the scientific name (you'd have to mod orthofly.pl though)
#
#See Also   : 
sub _get_interactor_taxon{ 
     my ($id) = @_;
     
     return $id if($id eq "-");
     return $id if( ($id eq "-2") or ($id eq "-1") ); # "chemical synthesis" or "in vitro"
     
     #there are unclean entries like "taxid:-1"
     
     if(($id =~ /(taxid):(-\d+)/) || ($id =~ /^(taxid):(\d+)(\(.*\))/)){
          return $2;
     }else{
          print "_get_interactor_taxon(): unrecognised entry format: '$id'\n";
          return "-";  
     }

}

#
#_get_vector_from_string
#
# Usage     : How to use this function/method
# Purpose   : This is used by get_backward_orthologies() and by get_real_interaction 
#             to clean up  the string containing the aliases. It will remove the "db:" bit
#             and save each alias as a variable in a vector
# Returns   : a vector containing the clean aliases
# Argument  : string in "dbid:alias|dbid:alias..." format
# Throws    : -
# Comment   : needs refining
#
#See Also   : 
sub _get_vector_from_string{
     my ($datastring) = @_;
     my @result = ();
     
     my @datavec = split(/\|/, $datastring);
     
     foreach my $item (@datavec){ #db:string
          if($item =~ /(\w+):(.+)/){
               push(@result, $2);
          }
     }
     return @result;
}

#
#_get_ensembl_id
#
# Usage     : How to use this function/method
# Purpose   : This is used by get_backward_orthologies() to obtain the Ensembl ID of a gene out of a UniprotKB id of a protein. A loss of 
#             information results from the conversion as isoforms are eliminated. It's acompromise and should be improved
# Returns   : the Ensembl ID of the gene 
# Argument  : a gene_adaptor from ensembl, the intact interaction id, the id of the orthologue of the original query gene, and
#             id, name, aliases of the gene whose ensembl id we're looking for
# Throws    : -
# Comment   : 
#
#See Also   : 
sub _get_ensembl_id{
     my %args = @_;
     my $gene_adaptor   = $args{adaptor};
     my $ebi_id         = $args{ebi_id};
     my $orthologueid   = $args{ortho_id};
     my $interactorid   = $args{acc_numb};
     my $name           = $args{protein_name};
     my $aliasvector    = $args{aliases};
     my $no_output      = $args{no_output};

     my $ensemblid;
     
     if($interactorid =~ /(\w+)(-\d{1})$/){
          $interactorid = $1;
     }
     
     #accessionnumber
     $ensemblid = _query_generic_id($gene_adaptor, $orthologueid, $interactorid);
     if(!$ensemblid){  #name
          #I want to know how many I couldn't get immediately by querying simply the id.
          $ENSEMBLIDFAILED = 1;    
          print "$ebi_id-($interactorid): query by accession number $interactorid ambiguous/empty, trying name..\n" unless($no_output);
          $ensemblid = _query_generic_id($gene_adaptor, $orthologueid, $name);
          if(!$ensemblid){ #aliases
               print "$ebi_id-($interactorid): query by name $name ambiguous/empty, trying aliases..\n" unless($no_output);
               my @aliases = _get_vector_from_string($aliasvector);
               foreach my $alias (@aliases){
                    $ensemblid = _query_generic_id($gene_adaptor, $orthologueid, $alias);
                    if($ensemblid){
                         print "Found by alias: $alias\n" unless($no_output);
                         return $ensemblid;
                    }
               }
               
               foreach my $alias (@aliases){
                    my $stableid = $gene_adaptor->fetch_by_stable_id($alias);
                    if (defined $stableid){
                         print("Found by alias: $alias\n") unless($no_output);
                         $ensemblid = $stableid->stable_id;
                         return $ensemblid;
                    }
               }
               
               if(!$ensemblid){
                    print "$ebi_id-($interactorid): name/aliases did not help, trying external_db string..\n" unless($no_output);
                    #last chance. I'll choose  'curated'/'automatic' over 'clone based'
                    my $gArray = $gene_adaptor->fetch_all_by_external_name($interactorid);
                    foreach my $gene (@$gArray){
                         my $externaldb = $gene->external_db;
                         if ($externaldb =~ /curated/){
                              return ($gene->stable_id);
                         }elsif ($externaldb =~ /automatic/){
                              return ($gene->stable_id);
                         }
                    }
                    print("$ebi_id-($interactorid): no way to distinguish multiple genes in array. Switching to UniprotKB ids..\n") unless($no_output);
                    return;
               }
          }
     }
     #print "Found by name: $name\n";
     return $ensemblid;
}

#
#_query_generic_id
#
# Usage     : How to use this function/method
# Purpose   : This is used by _get_ensembl_id() to query an id, or a name, or an alias, against Ensembl
# Returns   : the Ensembl ID of the gene or undefined
# Argument  : a gene adaptor, the ebi id of the interaction, the id of the initial interactor, the id we want to translate
# Throws    : -
# Comment   : 
#
#See Also   : 
sub _query_generic_id{
     my ($gene_adaptor, $orthologueid, $identifier) = @_;
     my $candidateID;
     
     return if($identifier eq '-');
     
     my $genesArray = $gene_adaptor->fetch_all_by_external_name($identifier);
          
     if(scalar(@$genesArray) == 1){
          $candidateID = @$genesArray[0]->stable_id;
          return  ($candidateID) if (defined $candidateID);
     }elsif(scalar(@$genesArray) == 0){
          return;
     }else{ #several gene objects are retrieved.
          foreach my $gene (@$genesArray){
               $candidateID = $gene->stable_id;
               if ($candidateID eq $orthologueid){
                   #check if one of the two is the initial interactor:
                   #if that is the case, we pick it and we throw away all the alternatives: what we need is the other one in the pair 
                   #anyway
                    return $candidateID;
               }
               #if we're here, there's ambiguity. Skip
          }
          return;
     }
     return;
}

#
#_process_homologies
#
# Usage     : How to use this function/method
# Purpose   : this is called by both get_forward_orthologies and by get_backward_orthologies and it is used to process the homology object
#             returned by ensembl, once for each gene id. This means in here the data fields are build and returned if necessary
# Returns   : a counter holding the number of orthologues found for the given id
# Argument  : the source id, the vector of homology members, a protein adaptor, the output file handle and (OPTIONAL) the data vector with the previous mined data
#             (used only if the routine is called by get_backward_orthologies) and an indicator saying we want only one-to-one orthologues
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 

#you have to return a counter
sub _process_homologies{
     my %args = @_;
     my $init_id        = $args{homology_query_id};
     #init_id is the original id provided by the user (if this routine is called in forward_orthologies)
     #or the interacting gene in the destination species (if the routine is called in backward orthologies)
     my $homologies     = $args{homology_vector};
     my $pt_adaptor     = $args{protein_adaptor};
     my $output_file    = $args{outfile};
     my $old_data_row   = $args{datavector};
     my $oo_only        = $args{hq_only};
     my $no_output      = $args{no_output};
     
     my $counter = 0;
     foreach my $homology (@{$homologies}){  
          my $data_line;
          my $ancestor;
          my $DF_orthologue_id;
          my $DF_oname; 
          my $DF_odesc;
          my $DF_dnds; 
          my $DF_opi;
          my $DF_fsa_x; 
          my $DF_fsa_y;
          my $DF_nndist;

          #the following filters the orthologies on the basis of their description             
          $DF_odesc = $homology->description;
          next if ($DF_odesc =~ /paralog/); #if you want paralogues instead work here
          next if (($oo_only) && (!($DF_odesc =~ /one2one/)));
          
          $counter += 1;
          #$osubtype = $homology->subtype; #not sure i need it for now
          #DNDS ratio
          $DF_dnds = $homology->dnds_ratio;#do not apply threshold for now
               
          my $genelist = $homology->gene_list;# genelist is a member object
          my $node_x; my $node_y;

          #genelist will contain ALL the genes in the orthology group: this means it will
          #also contain the original gene I made my query with.
          foreach my $homology_member (@{$genelist}){
               $DF_orthologue_id = $homology_member->stable_id;
               #my $mid = $homology_member->get_longest_peptide_Member->member_id; OBSOLETE
               my $mid = $homology_member->get_canonical_peptide_Member->member_id; #USE WITH V. > 57
               
               if($DF_orthologue_id eq $init_id){ #node_x will contain the one I started with
                    $node_x = $pt_adaptor->fetch_AlignedMember_by_member_id_root_id($mid,1);
               }else{
                    $node_y = $pt_adaptor->fetch_AlignedMember_by_member_id_root_id($mid,1);
               }
          }
                    
          if((defined($node_x)) && (defined($node_y) )){
               my $root = $node_x->subroot; #method from nested set
                    if (defined ($root)){
               eval{
                    $root->merge_node_via_shared_ancestor($node_y);
                    $ancestor = $node_x->find_first_shared_ancestor($node_y);   
                    };
                    if($@){
                         print "\nAn error occurred ($@), continuing..\n";
                    }
               if (defined $ancestor){
                    $DF_fsa_x = $node_x->distance_to_ancestor($ancestor);
                         $DF_fsa_y = $node_y->distance_to_ancestor($ancestor);
                         $DF_nndist = $node_x->distance_to_node($node_y);
                    }
               }
          }

          foreach my $homology_member (@{$genelist}){
               $DF_orthologue_id = $homology_member->stable_id;
               next if ($DF_orthologue_id eq $init_id);#I dont want to print again the gene name
          
               $DF_oname = $homology_member->display_label;

               #OPI
               my $pairwise_alignment_from_multiple = $homology->get_SimpleAlign;
               $DF_opi = $pairwise_alignment_from_multiple->overall_percentage_identity;
               #$opi = sprintf("%.3f", $overall_pid); #rounded
          
               $DF_orthologue_id = '-' if(!$DF_orthologue_id);
               $DF_oname         = '-' if(!$DF_oname); 
               $DF_odesc         = '-' if(!$DF_odesc);
               $DF_dnds          = '-' if(!$DF_dnds); 
               $DF_opi           = '-' if(!$DF_opi);
               $DF_fsa_x         = '-' if(!$DF_fsa_x); 
               $DF_fsa_y         = '-' if(!$DF_fsa_y);
               $DF_nndist        = '-' if(!$DF_nndist);     
          
               if($old_data_row){ #get_backward_orthologies - species_y is our species of interest
                    $data_line = join("\t",$old_data_row, 
                                   $DF_orthologue_id,$DF_oname, $DF_odesc, 
                                   $DF_opi, $DF_dnds, $DF_nndist,
                                   $DF_fsa_y, $DF_fsa_x);
                    print "Putative interactor found: $DF_oname ($DF_orthologue_id)\n" unless($no_output);
               }else{ #get_forward_orthologies - species_x is our species of interest
                    $data_line = join("\t",$init_id, 
                                   $DF_orthologue_id,$DF_oname, $DF_odesc, 
                                   $DF_opi, $DF_dnds, $DF_nndist,
                                   $DF_fsa_x, $DF_fsa_y);
               }
               print $output_file  $data_line, "\n";
               #$atLeastOneRowPrinted = 1;
          }
     }
     return $counter;
}

#
#_setup_dbi_connection
#
# Usage     : $dbh = _setup_dbi_connection($in_path);
# Purpose   : This subroutine sets up a DBD::CSV connection with the input file, which must be in TSV 
#             (Tab Separated Values) format. The file content is then accessed through DBI just like a 
#             Relational Database
# Returns   : A database handle
# Argument  : A string indicating the path to the file
# Throws    : -
# Comment   : -
#
#See Also   : 
sub _setup_dbi_connection{
     my ($path) = @_;
     
     my $dbh = DBI->connect("DBI:CSV:f_dir=.;csv_eol=\n;csv_sep_char=\t;csv_quote_char=;csv_escape_char=") 
     or die("Cannot connect: " . $DBI::errstr);
     
     $dbh->{'RaiseError'} = 1;
     local $@ = '';
     $dbh->{'csv_tables'}->{'int'} = { 'file' => $path};
     
     
     $dbh->{FetchHashKeyName} = 'NAME_uc';
     #weird behaviour..fetchrow_hashref might convert hash keys to lower case
     #got this from http://search.cpan.org/~rehsack/SQL-Statement-1.27/lib/SQL/Statement.pm
     
     return $dbh;
}

#############################################
#############################################
package Bio::Homology::InterologWalk::Scores;
#############################################
#############################################
use List::Util qw[min max sum];
use File::Glob;

#MITAB ontology core nodes
my $interaction_detection_method_acc = "MI:0001";
my $interaction_type_acc = "MI:0190";

########################################################
#look up table for the interaction scores.
my %first_level_hash = (
     #
     #INTERACTION TYPES
     #
     "MI:0403" =>   "0.1",    #colocalization
     "MI:0208" =>   "0",      #genetic interaction
     "MI:0914" =>   "0.5",    #association
     "MI:0915" =>   "0.8",    #physical association
     "MI:0407" =>   "1.1",    #direct interaction
     #
     #INTERACTION DET METHODS
     #
     "MI:0045" =>   "0.5",    #experimental interaction detection
     "MI:0362" =>   "0.2",    #inference
     "MI:0063" =>   "0.2",    #interaction prediction
     "MI:0686" =>   "0",      #unspecified method
);
#experimental method
my $SPOKE = 0;
my $NONSPOKE = 5;
my $INTERACTIONWEIGHT = 1;
my $ONTOLOGYWEIGHT = 1;
my $NNDIST_THRESHOLD = 100;
########################################################

=head2 parse_ontology

 Usage     : $onto_graph = 
                 Bio::Homology::InterologWalk::Scores::parse_ontology(
                                                                    $ont_path
                                                                    );
 Purpose   : This subroutine accepts one input, a path to a PSI-MI ontology file. 
             It uses GO::Parser to parse the file and returns a graph object of 
             the ontology: a structured graph-representation of it, that we can 
             walk and explore. This is useful when we need to look at the detection 
             method and at the interaction type for each entry. E.g. we might be in-
             terested in all interactions tagged generically with "experimental 
             detection method" but also in all the interactions tagged with 
             a *specific* detection method (ie a specialised subclass of the concept 
             "experimental detection method").
             Analysing the structure of the ontology through the graph returned by 
             this method helps in doing that.
 Returns   : A graph object containing the structured ontology
 Argument  : The path to the psi-mi ontology file
 Throws    : -
 Comment   : -

See Also   : L</get_mean_scores>

=cut


sub parse_ontology{
     my ($path) = @_;
     
     my $parser = GO::Parser->new({handler=>'obj'});
     $parser->parse($path); # parse file -> objects
     return ($parser->handler->graph);  # get L<GO::Model::Graph> object
}


=head2 get_mean_scores

 Usage     : ($m_em, $m_it, $m_dm, $m_mdm) = 
             Bio::Homology::InterologWalk::Scores::get_mean_scores(
                                                             $intact_path,
                                                             $onto_graph
                                                             );
 Purpose   : This is used to compute suitable mean values to normalise the components of 
             the score for
             - interaction type, 
             - interaction detection method
             - experimental method
             - multiple detection method.
             Each value is the MEAN value for the corresponding score computed on the set of 
             direct experimental interactions for the initial dataset. These are used to 
             normalise the scores obtained for the corresponding putative interactions.
 Returns   : a list of four numbers: 
             1. mean experimental method score
             2. mean interaction type score, 
             3. mean detection method score
             4. mean multiple dm score
 Argument  : 1) path to a tsv file of REAL intact interactions 
                (generated by get_direct_interactions())
             2) a graph representation of the obo PSI MI ontology 
                (generated by parse_ontology() )
 Throws    : -
 Comment   : -

See Also   : L</parse_ontology>, L</get_direct_interactions>

=cut


sub get_mean_scores{
     my ($file_path, $graph) = @_;
     
     my $mean_det_method_score = 0;
     my $mean_int_type_score = 0;
     my $mean_exp_method_score = 0;
     my $mean_ME_DM_score = 0; #multiple evidence, detection method
     my $interaction_count = 0;
     
     my %DMhash     = (); 
     my %IThash     = (); 
     my %EMhash     = (); 
     my %ME_DM_hash = ();
     
     my $dbh = Bio::Homology::InterologWalk::_setup_dbi_connection($file_path);
     #query--------
     my $query = "SELECT $FN_int_type, $FN_det_method, $FN_exp_method, $FN_multiple_dm
                              FROM int";
     my $sth = $dbh->prepare($query);
     $sth->execute() or die "Cannot execute: " . $sth->errstr();
     #-------------
     
     while (my $row = $sth->fetchrow_hashref){
          my $int_type;my $det_method;my $exp_method;my $multiple_det_meth;
          my $int_type_score;my $det_method_score;my $exp_method_score;
          $interaction_count += 1;
          
          $int_type                = $row->{$FN_int_type};
          $det_method              = $row->{$FN_det_method};
          $exp_method              = $row->{$FN_exp_method};
          $multiple_det_meth  = $row->{$FN_multiple_dm};
          
          #Multiple evidence through detection method - Mean ================
          $mean_ME_DM_score += $multiple_det_meth;
          $ME_DM_hash{$multiple_det_meth} += 1;
          
          #Expansion method============================
          #reward non-spoke-interactions
          if($exp_method =~ /spoke/i){
               $exp_method_score = $SPOKE;
          }else{
               $exp_method_score = $NONSPOKE;
          }
          $mean_exp_method_score += $exp_method_score;
          $EMhash{$exp_method_score} += 1;
          
          #Interaction type============================
          #we first want to look at the interaction type, in order to score the interaction detection method
          #within  its specific interaction class.
          $int_type_score = _score_interaction($graph, $int_type);
          $mean_int_type_score += $int_type_score;
          $IThash{$int_type_score} += 1;
          
          #Detection method============================
          #this should be combined to the interaction method routine to build the first 
          #component of the score.
          $det_method_score = _score_interaction($graph, $det_method);
          $mean_det_method_score += $det_method_score;
          $DMhash{$det_method_score} += 1;
     }
     
     $mean_exp_method_score = $mean_exp_method_score / $interaction_count; 
     $mean_int_type_score = $mean_int_type_score / $interaction_count;
     $mean_det_method_score = $mean_det_method_score / $interaction_count;
     $mean_ME_DM_score = $mean_ME_DM_score / $interaction_count;
     
     print "Mean Interaction Type Score: ", $mean_int_type_score, "\n";
     print "Mean Detection Method Score: ", $mean_det_method_score, "\n";
     print "Mean Experimental Method Score: ", $mean_exp_method_score, "\n";
     print "Mean Multiple DM Score: ", $mean_ME_DM_score, "\n";
     
     print "Interaction type Hist\n";
     foreach my $score (sort keys %IThash){
          $IThash{$score} = $IThash{$score}/$interaction_count;
          print "$score\t$IThash{$score}\n";
     }
     print "Detection Method Hist\n";
     foreach my $score (sort keys %DMhash){
          $DMhash{$score} = $DMhash{$score}/$interaction_count;
          print "$score\t$DMhash{$score}\n";
     }
     print "Experimental Method Hist\n";
     foreach my $score (sort keys %EMhash){
          $EMhash{$score} = $EMhash{$score}/$interaction_count;
          print "$score\t$EMhash{$score}\n";
     }
     
     print "Multiple DM Evidence Hist\n";
     foreach my $score (sort keys %ME_DM_hash){
          $ME_DM_hash{$score} = $ME_DM_hash{$score}/$interaction_count;
          print "$score\t$ME_DM_hash{$score}\n";
     }
     $sth->finish;
     $dbh->disconnect;
     return ($mean_exp_method_score, $mean_int_type_score, $mean_det_method_score, $mean_ME_DM_score);
}

#
#_get_multiple_taxa_mean_score
#
# Usage     : How to use this function/method
# Purpose   : This is used to compute a mean value for the "multiple taxa" score. RATIONALE: if the same putative interaction is obtained
#             by projection of a real interaction in MORE than one species, its MT score should be higher
# Returns   : a mean multiple taxa score, used to weight the actual value for the putative interaction
# Argument  : a db handle to a tsv file of intact interactions
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 
#TODO try to fuse this with the former, maybe passing two db handles, maybe making this optional
sub _get_multiple_taxa_mean_score{
     my ($dbh) = @_;
     my $mean_ME_TAXA_score = 0;
     my $interaction_count = 0;
     my %ME_TAXA_hash = (); #multiple evidence, taxa
     
     my $query = "SELECT $FN_multiple_taxa
                    FROM int";
     my $sth = $dbh->prepare($query);
     $sth->execute() or die "Cannot execute: " . $sth->errstr();
     
     if(!$sth->rows){
          print("_get_multiple_taxa_mean_score(): no data in input file. Skipping..");
          return;
     } # TODO temp solution do a SELECT COUNT etc etc
     
     while (my $row = $sth->fetchrow_hashref){
          $interaction_count += 1;
          my $multiple_taxa_score = $row->{$FN_multiple_taxa};
          $mean_ME_TAXA_score += $multiple_taxa_score;
          $ME_TAXA_hash{$multiple_taxa_score} += 1;
     }
     
     $mean_ME_TAXA_score = $mean_ME_TAXA_score / $interaction_count;
     print "Mean Multiple TAXA Score: ", $mean_ME_TAXA_score, "\n";

     print "Multiple TAXA Evidence Hist\n";
     foreach my $score (sort keys %ME_TAXA_hash){
          $ME_TAXA_hash{$score} = $ME_TAXA_hash{$score}/$interaction_count;
          print "$score\t$ME_TAXA_hash{$score}\n";
     }
     print "\n";
     $sth->finish;
     return $mean_ME_TAXA_score;
}

=head2 compute_scores

 Usage     : $RC = Bio::Homology::InterologWalk::Scores::compute_scores(
                                                        input_path        => $in_path,
                                                        score_path        => $score_path,
                                                        output_path       => $out_path,
                                                        term_graph        => $onto_graph,
                                                        meanscore_em      => $m_em,
                                                        meanscore_it      => $m_it,
                                                        meanscore_dm      => $m_dm,
                                                        meanscore_me_dm   => $m_mdm,
                                                        meanscore_me_taxa => $m_mtaxa,
                                                        no_output         => 0
                                                        );
 Purpose   : This is used to analyse several ancillary data fields obtained alongside the actual 
             putative PPI IDs and collate them into a global confidence score, which should provide 
             a measure of the reliability of each putative PPI. The score will take into account 
             a number of variables related to each of the steps involved in the orthology walk. 
             We can divide the meta-score related components in two broad classes:
             - parameters related to the interaction. These include: Interaction Type, Interaction 
               Detection Method, Interaction coming from a SPOKE-expanded complex, interaction recon-
               firmed through multiple taxa, interaction reconfirmed through multiple detection methods
             - parameters related to the two orthology mappings. These include: orthology type 
               (one-to-one, one-to-many, many-to-one, many-to-many), OPI (percentage identity of the 
               conserved columns - see Bio::SimpleAlign), node to node distance, distance from the 
               first shared ancestor, (under development) dN/dS ratio
             The score computation will also involve a normalisation stage. The subroutine requires 
             five arguments (meanscore_x) representing mean values to be used for normalisation.
             The actually means are computed in get_mean_scores(), which is pre-requisite to 
             compute_scores().
 Returns   : success/failure
 Argument  : -input_path : path to the input tsv file. A suitable input for this subroutine is the 
             final output of the orthology walk pipeline (see doInterologWalk.pl for usage guidelines).
              input file should have .06out extension
             -(OPTIONAL) score_path : path to  text file where scores will be saved one per row 
              (useful for looking at score distributions eg through matlab)
              output textfile has a .scores extension 
             -output_path : where you want the routine to write the data. Data is in TSV format. 
              File extension is .07out
             -term_graph :  a Go::Parser graph object obtained from parse_ontology() containing a 
              network representation of the PSI-MI controlled vocabulary of terms.
             -meanscore_em : mean experimental method score for normalisation
             -meanscore_it : mean interaction type score for normalisation  
             -meanscore_dm : mean detection method score for normalisation   
             -meanscore_me_dm : mean 'multiple detection methods' score for normalisation
             -meanscore_me_taxa : mean 'multiple taxa' score for normalisation
             -(OPTIONAL) no_output :  suppresses screen output. Used for clearer output during test. 
              Default is 0.
 Throws    : -
 Comment   : -

See Also   : L<http://search.cpan.org/~cjfields/BioPerl-1.6.1/Bio/SimpleAlign.pm#overall_percentage_identity>, L</get_mean_scores>, C<doScores.pl> for sample usage

=cut


sub compute_scores{
     my %args = @_;
     my $in_path                   = $args{input_path};
     my $out_path                  = $args{output_path};
     my $score_path                = $args{score_path};
     my $graph                     = $args{term_graph};
     my $mean_exp_method_score     = $args{meanscore_em};
     my $mean_int_type_score       = $args{meanscore_it};
     my $mean_det_method_score     = $args{meanscore_dm};
     my $mean_multiple_dm_score    = $args{meanscore_me_dm}; #multiple evidence, detection method 
     my $mean_multiple_taxa_score  = $args{meanscore_me_taxa}; #multiple evidence, taxa
     my $no_output                 = $args{no_output};
     
     if(!$graph){
          print("compute_scores(): PSI-MI graph representation not found. Aborting..\n");
          return;
     }
     unless($mean_exp_method_score  && 
            $mean_int_type_score    && 
            $mean_det_method_score  && 
            $mean_multiple_dm_score &&
            $mean_multiple_taxa_score){
                 
                 print("compute_scores(): missing mean scores. Aborting..\n");
                 return;
     }
     
     
     my @nodeNodeDist = ();
     my $MAXNNDISTANCE = 0;
     my $totalDIST = 0;
     my $score_data;
     
     #MANAGE FILES
     my $dbh = Bio::Homology::InterologWalk::_setup_dbi_connection($in_path);
     open (my $out_data, q{>}, $out_path) or croak("Unable to open $out_path : $!");
     if($score_path){
         open ($score_data,    q{>}, $score_path) or croak("Unable to open $score_path : $!"); 
     }
     #============
     
     #header----
     print $out_data $HEADER_FWD_ORTH,  "\t", 
                     $HEADER_INT,        "\t", 
                     $HEADER_BWD_ORTH,   "\t", 
                     $HEADER_MISC,       "\t",
                     $HEADER_SCORE,      "\n";
     
     ##########
     # 1st pass
     ##########
     # 1)need to retrieve the highest node to node distance value for weighting purposes
     # 2)need to retrieve some max/mean value for the multiple taxa score.
     my $query = "SELECT $FN_nndist_1, $FN_nndist_2, $FN_multiple_taxa
                         FROM int";
                    
     my $sth = $dbh->prepare($query);
     $sth->execute() or die "Cannot execute: " . $sth->errstr();
     
     #Ensembl might contain spurious outliers in the node to node dist (100.000).
     while (my $row = $sth->fetchrow_hashref){
          my $nnDist_1 = $row->{$FN_nndist_1};
          push(@nodeNodeDist, $nnDist_1) if ($nnDist_1 < $NNDIST_THRESHOLD); 
          my $nnDist_2 = $row->{$FN_nndist_2};
          push(@nodeNodeDist, $nnDist_2) if ($nnDist_2 < $NNDIST_THRESHOLD);
          $totalDIST += 2;
     }
     #print scalar @nodeNodeDist, " ", $totalDIST, "\n";
     $MAXNNDISTANCE = max(@nodeNodeDist);
     $sth->finish;
     
     ##########
     # actual score computations
     ##########
     $query = "SELECT * FROM int";
                         
     my $sthHash = $dbh->prepare($query);
     my $sthVec = $dbh->prepare($query);
     $sthHash->execute() or die "Cannot execute: " . $sthHash->errstr();
     $sthVec->execute() or die "Cannot execute: " . $sthVec->errstr();
     
     while (my $rowHash = $sthHash->fetchrow_hashref){ 
          my @rowVec = $sthVec->fetchrow_array();
          my $compoundScore = 0;
     
          #I'll need to gather all score-related fields
          my $opi_1        = $rowHash->{$FN_opi_1};
          my $nnDist_1     = $rowHash->{$FN_nndist_1};
          my $intType      = $rowHash->{$FN_int_type};
          my $detMethod    = $rowHash->{$FN_det_method};
          my $expMethod    = $rowHash->{$FN_exp_method};
          my $opi_2        = $rowHash->{$FN_opi_2};
          my $nnDist_2     = $rowHash->{$FN_nndist_2};
          my $multipleDM   = $rowHash->{$FN_multiple_dm};
          my $multipleTaxa = $rowHash->{$FN_multiple_taxa};
          #my $timesSeen   = $rowHash->{TIMES_SEEN}; 
          
          #####################
          # INTERACTIONS
          #####################
          
          #Interaction type============================
          my $intTypeScore    = _score_interaction($graph, $intType, $no_output);
          $intTypeScore       = $intTypeScore / $mean_int_type_score;
          
          #Detection method============================
          my $detScore   = _score_interaction($graph, $detMethod, $no_output);
          $detScore      = $detScore / $mean_det_method_score;
          
          my $multipleDMScore = $multipleDM / $mean_multiple_dm_score;
          my $multipleTaxaScore = $multipleTaxa / $mean_multiple_taxa_score;
          #print SCORE_DATA $multipleTaxaScore, "\n";
          
          #Expansion method============================
          my $expMethodScore;
          #penalise spoke-interactions
          if($expMethod =~ /spoke/i){
               $expMethodScore = $SPOKE;
          }else{
               $expMethodScore = $NONSPOKE;
          }
          # I leave it unweighted for now, to generate a "step" in the distribution
          $expMethodScore = $expMethodScore / $mean_exp_method_score;
          
          my $S_I = $expMethodScore + $intTypeScore + $detScore + $multipleDMScore + $multipleTaxaScore;
          
          
          #################################
          # ORTHOLOGIES
          #################################
          #Overall Percentage Identity================
          $opi_1 = $opi_1 / 100;
          $opi_2 = $opi_2 / 100;
          my $jointSeqIdentity = sqrt ($opi_1 * $opi_2); 
          
          #Node to node distance
          my $highest_nn_dist = max($nnDist_1, $nnDist_2);
          $highest_nn_dist = $highest_nn_dist / $MAXNNDISTANCE;
          #after normalisation I should have all values less than one, apart from outliers, that I discard:
          $highest_nn_dist = 1 if ($highest_nn_dist > 1);
          
          my $jointNodeDifference = 1 - $highest_nn_dist;
          
          my $S_O   =    $jointNodeDifference + $jointSeqIdentity;
          
          # SUMMING UP
          $compoundScore = $INTERACTIONWEIGHT * $S_I + $ONTOLOGYWEIGHT * $S_O;
          
          #prints scores in data file for printing histograms
          my $dataline = join("\t",$intTypeScore,$detScore,$multipleDMScore,$multipleTaxaScore, $jointSeqIdentity, $jointNodeDifference);
          #print SCORE_DATA $dataline, "\n";
          print $score_data $compoundScore, "\n" if($score_data);
          
          my $newRow = join("\t", @rowVec, $compoundScore);
          print $out_data $newRow, "\n";
     }
     $sthHash->finish;
     $sthVec->finish;
     $dbh->disconnect;
     close $out_data;
     close $score_data if($score_data);
     return 1;
}


=head2 compute_multiple_taxa_mean

 Usage     : $m_mtaxa = Bio::Homology::InterologWalk::Scores::compute_multiple_taxa_mean(
                                                  ds_size    => 500,   
                                                  ds_number  => 3,    
                                                  datadir    => $path     
                                                  );
 Purpose   : Suppose you want to run an interolog walk starting from initial interactor X. 
             You get a final putative interactor Y.
             Now suppose the putative interactor is the output of more than an interolog walk,
             each one based on a interaction annotated in a different organism.
             When building the score, we would like to account for the fact that a putative 
             interaction obtained through interactions in multiple species is more reliable than 
             one obtained through only one species. In order to weight the Multiple_Taxa Score, 
             however, a long procedure is required. It is not possible to use a mean taken from 
             the direct Intact interactions data file, for obvious reasons. 
             My solution can be currently summarised as follows:
             1. choose n<7 random taxa from a vector containing 7 well-supported NCBI taxa id;
             2. choose m (~500) random genes for each of the n taxa
             3. run the full orthology walk using the methods from Bio::Homology::InterologWalk 
                on each of the n datasets
             4. compute a mean_multiple taxa score for each of them
             5. Final Mean_Multiple_Taxa_Score = mean(M_1,M_2,M_3) 
             The procedure is LONG and SLOW and might lead to Ensembl refusing connections 
             in some instances.
 Returns   : the mean multiple taxa score, i.e. the global mean of the multiple taxa scores 
             obtained for each random dataset
 Argument  : ds_size : number of ids per dataset, eg 500
             ds_number : a number between 1 and 7, equal to the number of taxa to randomly 
                         pick from
             datadir :  work directory
 Throws    : -
 Comment   : #TODO should randomised data be saved and reused?
             Lots of hard coded stuff in here. Intact url is hard coded, etc. 
             Need to review.


See Also   : 

=cut


sub compute_multiple_taxa_mean{
     my %args = @_;
     my $dataset_size         = $args{ds_size};
     my $dataset_number       = $args{ds_number};
     my $dataset_dir          = $args{datadir};
     
     if(!$dataset_size){
          print("compute_multiple_taxa_mean(): no dataset size specified. Aborting\n");
          return;
     }
     
     my $in_path; my $out_path;
     my $partial_score;
     my $SCORE;
     my @scores;
    
     my $rand = "Rand_";
     my $registry = 'Bio::EnsEMBL::Registry';
     $registry->load_registry_from_db(
               -host          => 'ensembldb.ensembl.org',
               -user          => 'anonymous',
               );
     
     #first generate datasets
     my @selected_taxa;
     my %seen = (); #to avoid replacement

     my @taxa = (
          #"6239", #Cele #this has been removed from Ensembl Multi Compara. If you need this data you'll want to get a registry
          #through the multiple call.
          "7227", #Dele
          "9606", #Hsap
          "10090", #Mmus
          "10116", #Rnor
          "4932", # Scer
          "7955", #Drer
     );
     my $NCBI_taxon_adaptor = $registry->get_adaptor("Multi", "compara", "NCBITaxon");
    
     my $total_taxa = scalar(@taxa);
     if($dataset_number  > $total_taxa){
         print("compute_MTAXA_mean(): ds_number is higher than the number of available taxa ($total_taxa). Aborting..");
         exit 0;
     }
     my $i = 0;
     do{
          my $number = int(rand($total_taxa));
          unless($seen{$number}){
               push(@selected_taxa, $taxa[$number]);
               $seen{$number} = 1;
               $i++;
          }
     }while($i < $dataset_number);
     
     my $tempfiles = $dataset_dir . "Rand*";
     unlink glob(($tempfiles));
     
     foreach my $NCBItaxonID (@selected_taxa){
          my $dataset; #this is the input dataset to use for the orthology walk
          my $NCBItaxon = $NCBI_taxon_adaptor->fetch_node_by_taxon_id($NCBItaxonID);
          my $NCBItaxon_name = $NCBItaxon->binomial;
          print "======================\n";
          print "RANDOM DATASET: ", $NCBItaxon_name, "\n";
          print "======================\n";
          my $gene_adaptor = $registry->get_adaptor($NCBItaxon_name, "core", "Gene");
          if(!$gene_adaptor){
               print("No gene adaptor found for $NCBItaxon_name. Skipping..");
               next;
          }
          my @genes = @{$gene_adaptor->list_stable_ids()};# these are ALL the gene ids for this taxon
          my $total_genes = (scalar(@genes) - 1);
          for(my $j = 0; $j < $dataset_size; $j++){
               #print $genes[int(rand($total_genes))], "\n"; 
               push(@$dataset, $genes[int(rand($total_genes))]);
          }
          
          #let's put the dataset in a file============
          my $gene_data_path = $dataset_dir . $rand . $NCBItaxon->short_name . ".txt";
          open (my $gene_data,  q{>}, $gene_data_path) or croak("Unable to open $gene_data_path : $!");
          foreach my $id (@{$dataset}){
               print $gene_data $id, "\n";
          }
          close $gene_data;
          #===========================================
          $out_path = $dataset_dir . $rand . $NCBItaxon->short_name . ".out01";
          #I have the dataset, I'll do the full interolog walk on it and also compute its mean taxa score
          my $RC1 = Bio::Homology::InterologWalk::get_forward_orthologies(registry          => $registry,
                                                ensembl_db        => 'multi',
                                                input_path        => $gene_data_path,
                                                output_path       => $out_path,
                                                source_org        => $NCBItaxon_name
                                                );
          if(!$RC1){
               print "compute_multiple_taxa_mean(): get_forward_orthologies() returned errors. Stopping..\n";
               return;
          }
          $in_path = $out_path;
          $out_path = $dataset_dir . $rand . $NCBItaxon->short_name . ".out02";
          my $dbh = Bio::Homology::InterologWalk::_setup_dbi_connection($in_path);
          my $url = "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/";
          my $RC2 = Bio::Homology::InterologWalk::get_interactions(input_path    => $in_path,
                                         output_path   => $out_path,
                                         url           => $url,
                                         );
          if(!$RC2){
               print "compute_multiple_taxa_mean(): get_interactions() returned errors. Stopping..\n";
               return;
          }
          unlink $in_path;
          $in_path = $out_path;
          $out_path = $dataset_dir . $rand . $NCBItaxon->short_name . ".out03";
          my $RC3 = Bio::Homology::InterologWalk::get_backward_orthologies(registry     => $registry,
                                                 ensembl_db   => 'multi',
                                                 input_path   => $in_path,
                                                 output_path  => $out_path,
                                                 source_org   => $NCBItaxon_name,
                                                 );
          if(!$RC3){
               print "compute_multiple_taxa_mean(): get_backward_orthologies() returned errors. Stopping..\n";
               return;
          }
          unlink $in_path;
          $in_path = $out_path;
          $out_path = $dataset_dir . $rand . $NCBItaxon->short_name . ".out04";
          my $RC4 = Bio::Homology::InterologWalk::remove_duplicate_rows(input_path    => $in_path,
                                              output_path   => $out_path,
                                              header        => 'standard',
                                              );
          if(!$RC4){
               print "compute_multiple_taxa_mean(): remove_duplicate_rows() returned errors. Stopping..\n";
               return;
          }
          unlink $in_path;
          $in_path = $out_path;
          $out_path = $dataset_dir . $rand . $NCBItaxon->short_name . ".out";
          my $RC5 = Bio::Homology::InterologWalk::do_counts(input_path    => $in_path,
                                  output_path   => $out_path,
                                  header        => 'standard',
                                  );
          if(!$RC5){
               print "compute_multiple_taxa_mean(): do_counts() returned errors. Stopping..\n";
               return;
          }
          unlink $in_path;
          $in_path = $out_path;
          $dbh = Bio::Homology::InterologWalk::_setup_dbi_connection($in_path);                                                         
          $partial_score = Bio::Homology::InterologWalk::Scores::_get_multiple_taxa_mean_score($dbh);
          push(@scores, $partial_score) if($partial_score);
          $dbh->disconnect;
     }
     $SCORE = _average(\@scores);
     return $SCORE;
}



####################################
#HELPER routines for internal usage
####################################

sub _average { 
     #@_ == 1 or die ('Sub usage: $average = average(\@array);'); 
     my ($array_ref) = @_; 
     my $sum; 

     my $count = scalar @$array_ref; 
     foreach (@$array_ref) { $sum += $_; } 
     
     return $sum / $count; 
} 

#
#_get_ont_id
#
# Usage     : How to use this function/method
# Purpose   : ####write
# Returns   : A gene adaptor for the selected species
# Argument  : A string with the full scientific name of the organism as recognised by Ensembl
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 
sub _get_ont_id{
     my ($string) = @_;

     if($string =~ /(.+)(MI:\d+)(.+)/){
          return $2;     
     }
     return $string;
}

#
#_score_interaction
#
# Usage     : How to use this function/method
# Purpose   : ####write
# Returns   : A gene adaptor for the selected species
# Argument  : A string with the full scientific name of the organism as recognised by Ensembl
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 
sub _score_interaction{
     my ($graph, $termToScore, $no_output) = @_;
     my $globalScore = 0;
     my @pathScores = (); #in case there are more paths leading up to different metaconcepts. I'll extract the max from here
     
     #get the concept from the interaction and build and store each of the paths that go up
     #to one of the main four subclasses
     my $termID = _get_ont_id($termToScore);
     #if no detection method has been annotated at Intact, then I cannot take further action
     if($termID eq "-"){
          return 0;
     }
     #if the term is already one of the scored ones, there's no need to traverse the ontology. Just get the score
     if(exists($first_level_hash{$termID})){
          $globalScore += $first_level_hash{$termID};
          unless($no_output){
               printf "\n\n%s -> ", $termID;
               print "FOUND\n\n";
          }
          return $globalScore; 
     }

     my $term = $graph->get_term($termID);
     if(!$term){
          print "WARNING: the ontology does not return a term object for term: $termID. Skipping..\n";
          return 0;#it could be that the ontologyID string is not totally clean
     }
     my $termName = $term->name;
     my $paths = $graph->paths_to_top($termID);
     
     # a concept could have 1 or more parents
     #I'm going to look for one of the root concepts in all of them
     #the following can happen
     #1)no path to top contains any of the four concepts: this should never happen, but if it happens, trap
     #2)more than one path has the root  concepts:
     #    2.a)all are the same (e.g. all "experimental"). It means the branching is below. Then use anyone of them
     #    2.b)they're different (e.g. "experimental" and "inferred"). Use the one with the highest confidence (experimental)
     #only one path has the root concept: straightforward.
     print "\nClimbing ontology in search of known meta-concepts..\n" unless($no_output);
     
     if(scalar @$paths == 1){
          print "Only one path to top exists. Traversing it..\n\n" unless($no_output);
          my $tempAcc;
          $graph->iterate(sub {
                                   $term=shift->term;
                                   $tempAcc = $term->acc;
                                   unless($no_output){
                                       printf "%s %s -> \n", $tempAcc,$term->name; 
                                   }
                                   if($first_level_hash{$tempAcc}){
                                        $globalScore += $first_level_hash{$tempAcc};
                                        print "FOUND\n" unless($no_output);
                                        return;
                                   }
                              },
                         {direction=>'up',acc=>$termID});
                         
     }elsif(scalar @$paths > 1){
          print "More than one path to top. Traversing all of them..\n\n" unless($no_output);   
          for my $pathtoTop (@{$paths}){
               my $termlist = $pathtoTop->term_list;
               print "$termName", " ->> " unless($no_output);
               for my $terminPath (@{$termlist}){
                    unless($no_output){
                         printf "%s %s -> \n", $terminPath->acc,$terminPath->name;   
                    }
                    my $tempAcc = $terminPath->acc;
                    if($first_level_hash{$tempAcc}){
                         push(@pathScores, $first_level_hash{$tempAcc});
                         #$global_DetMethod_Score += $first_level_hash{$tempAcc};
                         print "FOUND. Stopping search on this path.\n" unless($no_output);
                         $globalScore = max(@pathScores);
                         return $globalScore;
                    }
               }
          }
          $globalScore = max(@pathScores);
          
     }else{#no paths
          print "No paths to top..Skipping\n" unless($no_output);
          return 0;
     }
     return $globalScore;
}



###############################################
###############################################
package Bio::Homology::InterologWalk::Networks;
###############################################
###############################################
use Carp qw(croak);

=head2 do_network

 Usage     : $RC = Bio::Homology::InterologWalk::Networks::do_network(
                                                      registry       => $registry,
                                                      input_path     => $in_path,
                                                      output_path    => $out_path,
                                                      source_org     => $sourceorg,
                                                      orthology_type => $orthtype,
                                                      expand_taxa    => 1,
                                                      ensembl_db     => $ensembl_db,
                                                      no_output      => 0
                                                      );
 Purpose   : This function  writes a .SIF file according to this cytoscape specification in:
             http://cytoscape.org/cgi-bin/moin.cgi/Cytoscape_User_Manual/Network_Formats.
             For each input data row, the subroutine will extract the initial id, the 
             (putative) PPI found, taxon information(optional), score (if present). 
             It will look up the ID pair on the Ensemble API and obtain gene names. 
             It will output a TSV file with .sif extension.
             The routine can expand taxa information: it might be useful in some cases to 
             know the taxon from which the putative PPI has been mapped
             Eg.: instead of  A--B (default behaviour) one can decide to get A-mouse-B, 
             A-human-B, A-fly-B, etc.
             The routine should work both with a putative PPIs data file and with a direct
             interactions data file 
             (it will look at the input file header to decide what it is dealing with).
 Returns   : success/ failure
 Argument  : -registry: ensembl registry object to connect to. Needed to retrieve up-to-date
              human readable gene names from Ensembl for the IDs in the input data file
             -input_path : path to input file. Input file for this subroutine is a tsv file 
              containing at least the fields INIT_ID and INTERACTOR.
              (output of get_backward_orthologies() or get_direct_interactions will work, 
              although output of remove_duplicate_rows() is recommended). 
              Optionally, if interaction scores are desired, input fill will have to be the 
              output of the scoring pipeline (see example file doScores.pl)
             -output_path : where you want the routine to write the data. Data is a TSV .sif 
              cytoscape file.
             -source organism name (eg: "Mus musculus")
             -(OPTIONAL) orthology_type: can be set to 'onetoone': if so, only entries 
              obtained through "one to one" orthology projections will be retained in the output. 
              Default: all orthologies retained.
             -(OPTIONAL) expand_taxa: if true, information related to the species from which 
              the putative PPI has been projected will be retained.
              if this is set to true, ensemb_db MUST also be set
             -(OPTIONAL) ensembl_db: only required if expand_taxa is set. For allowed values, 
              see get_forward_orthologies()
             -(OPTIONAL) no_output :  suppresses screen output. Used for clearer output during
              test. Default is 0.
 Throws    : -
 Comment   : in order to account for the fact that the edges of the network are undirected 
             (eg A-B = B-A), I can either
             1) cache both couples (eg (A,B) and (B,A)) so I'll find them both when I look up
             2) do a lexicographic sorting before caching and before looking-up the cache
             I use the second option at the moment.

See Also   : L</remove_duplicate_rows>, L</compute_scores>, L</get_forward_orthologies>

=cut

sub do_network{
     my %args = @_;
     my $registry           = $args{registry};
     my $in_path            = $args{input_path};
     my $out_path           = $args{output_path};
     my $sourceorg          = $args{source_org}; 
     my $ortology_class     = $args{orthology_type};
     my $expand_taxa        = $args{expand_taxa}; #default no     
     my $ensembl_db         = $args{ensembl_db};
     my $no_output          = $args{no_output};
     
     #checks
     if(!$registry){
          print("do_network(): no registry defined. Aborting..\n");
          return;
     }
     if(!$sourceorg){
          print("do_network(): no source organism specified. Aborting..\n");
          return;
     }
     if($expand_taxa){
          if(!$ensembl_db){
               print("do_network(): No ensembl db specified. An ensembl db is needed for if expand_taxa is set. Aborting..\n");
               return;
          }
     }
     $ortology_class = 'allortho' if(!$ortology_class);
     
     my $query;
     my $scored;
     my $it = "pp"; #following cytoscape standard
     my $dbh; my $sth;
     
     #MANAGE FILES
     open (my $in_data,   q{<}, $in_path) or croak("Unable to open $in_path : $!");
     open (my $out_data,  q{>}, $out_path) or croak("Unable to open $out_path : $!");
     $dbh = Bio::Homology::InterologWalk::_setup_dbi_connection($in_path);
     #============
     
     #read the header, check whether the compound score is present, and create the following query accordingly
     #I'll read the header and see if there's a field containing SCORE. To be improved upon
     print("Checking for the presence of a compound score field in the input TSV data..\n") unless($no_output);
     my $old_header = <$in_data>;
     close $in_data;
     if($old_header =~ /$FN_compound_score/){ #need to check if the score had been computed
          $scored = 1;
          print("Compound score present, including it in the cytoscape network file..\n") unless($no_output);
          print $out_data $HEADER_SIF, "\t", $FN_compound_score, "\n";
          
          if($expand_taxa){ #composite score, taxon info
              print("Expanding binary interactions through taxon information..\n\n") unless($no_output);
              if($ortology_class eq "onetoone"){
                   $query = "SELECT $FN_initid,$FN_interactor_id,$FN_compound_score, $FN_taxon_a
                             FROM int
                             WHERE ($FN_odesc_1 like '%one2one%') 
                             AND ($FN_odesc_2 like '%one2one%')
                             AND ($FN_taxon_a = $FN_taxon_b)";
              }elsif( ($ortology_class eq "allortho") or ($ortology_class eq "direct") ){
                   $query = "SELECT $FN_initid, $FN_interactor_id, $FN_compound_score, $FN_taxon_a
                             FROM int
                             WHERE ($FN_taxon_a = $FN_taxon_b)";     
              }else{
                  print "do_network(): impossible to select a query, bad command line argument: $ortology_class. Aborting..\n";
                  return;
              }
          }else{ #composite score, no taxon
              print("Ignoring taxon information..\n\n") unless($no_output);
              if($ortology_class eq "onetoone"){
                   $query = "SELECT $FN_initid,$FN_interactor_id,$FN_compound_score
                             FROM int
                             WHERE ($FN_odesc_1 like '%one2one%') 
                             AND ($FN_odesc_2 like '%one2one%')";
              }elsif( ($ortology_class eq "allortho") or ($ortology_class eq "direct") ){
                   $query = "SELECT $FN_initid, $FN_interactor_id, $FN_compound_score
                             FROM int";     
              }else{
                  print "do_network(): impossible to select a query, bad command line argument: $ortology_class. Aborting..\n";
                  return;
              }
          }
     }else{ 
          print("Compound score not found, proceeding without..\n\n") unless($no_output);
          print $out_data $HEADER_SIF, "\n";
          
          if($expand_taxa){ #no composite score, taxon info
               if($ortology_class eq "onetoone"){
                    $query = "SELECT $FN_initid,$FN_taxon_a, $FN_interactor_id
                         FROM int
                         WHERE ($FN_odesc_1 like '%one2one%') 
                         AND ($FN_odesc_2 like '%one2one%')
                         AND ($FN_taxon_a = $FN_taxon_b)";
               }elsif(($ortology_class eq "allortho") or ($ortology_class eq "direct") ){
                    $query = "SELECT $FN_initid,$FN_taxon_a, $FN_interactor_id
                              FROM int
                              WHERE ($FN_taxon_a = $FN_taxon_b)";     
               }else{
                    print "do_network(): impossible to select a query, bad command line argument: $ortology_class. Aborting..\n";
                    return;
               }
          }else{ #no composite score, no taxon info
               if($ortology_class eq "onetoone"){
                    $query = "SELECT $FN_initid, $FN_interactor_id
                         FROM int
                         WHERE ($FN_odesc_1 like '%one2one%') 
                         AND ($FN_odesc_2 like '%one2one%')";
               }elsif(($ortology_class eq "allortho") or ($ortology_class eq "direct") ){
                    $query = "SELECT $FN_initid, $FN_interactor_id
                              FROM int";     
               }else{
                    print "do_network(): impossible to select a query, bad command line argument: $ortology_class. Aborting..\n";
                    return;
               }
          }
     }
     
     #if we have found a score column, we'll need to cache the entries to only keep, for each A-B (B-A) the entry with the 
     #highest score. This is based on the insight that a putative PPI is as reliable as the most highly scoring walk it has been
     #found in.
     #
     #PRE-CACHING
     #
     my %pair_max_score = ();
     if($scored){
          $sth = $dbh->prepare($query);
          $sth->execute() or die "Cannot execute: " . $sth->errstr();
          
          while (my $row = $sth->fetchrow_hashref) {
               my $id_in = $row->{$FN_initid};
               my $id_out = $row->{$FN_interactor_id};
               my $score = $row->{$FN_compound_score};
               my $ncbi_taxon_id; 
               $ncbi_taxon_id = $row->{$FN_taxon_a} if($expand_taxa);               

               my $ab_entry = join("",sort($id_in, $id_out));
               $ab_entry = join("-",$ab_entry, $ncbi_taxon_id) if($expand_taxa);
               
               if(!$pair_max_score{$ab_entry}){
                    $pair_max_score{$ab_entry} = $score;
               }else{
                    if($score > $pair_max_score{$ab_entry}){
                         $pair_max_score{$ab_entry} = $score;
                    }
               }
          }
          $sth->finish();
     }
     
     #
     # 2nd PASS
     #
     $sth = $dbh->prepare($query);
     $sth->execute() or die "Cannot execute: " . $sth->errstr();
     #TCHECK FOR the case of NO entries (eg no one to one orthologies)
     
     #ensembl adaptors---------- 
     my $source_gene_adaptor = $registry->get_adaptor($sourceorg, 'core', 'Gene');
     my $ncbi_taxon_adaptor;
     $ncbi_taxon_adaptor = $registry->get_adaptor($ensembl_db, "compara", "NCBITaxon") if($expand_taxa);
     
     my %coupleexists = ();
     while (my $row = $sth->fetchrow_hashref) {
          my $ncbi_shortname;
          my $ncbi_taxon_id;
          my $new_row;
          my $ab_entry;
          
          my $id_in = $row->{$FN_initid};
          my $id_out = $row->{$FN_interactor_id};
          $ncbi_taxon_id = $row->{$FN_taxon_a} if($expand_taxa);
          
          $ab_entry = join("",sort($id_in, $id_out));
          $ab_entry = join("-",$ab_entry, $ncbi_taxon_id) if($expand_taxa);
          next if($coupleexists{$ab_entry});
          $coupleexists{$ab_entry} = 1;

          if($scored){
               if($expand_taxa){
                    my $ncbi_taxon = $ncbi_taxon_adaptor->fetch_node_by_taxon_id($ncbi_taxon_id);
                    if (!$ncbi_taxon){
                         print "do_network(): can't find a taxon object for this: $ncbi_taxon_id\n";
                         $ncbi_shortname = "-";
                    }else{
                         $ncbi_shortname = $ncbi_taxon->short_name;
                    }
                    $new_row = join("\t", $id_in,  $ncbi_shortname,  $id_out, $pair_max_score{$ab_entry});
               }else{
                    $new_row = join("\t", $id_in,$it,$id_out, $pair_max_score{$ab_entry});
               }
          }else{ #not scored
               if($expand_taxa){
                    my $ncbi_taxon = $ncbi_taxon_adaptor->fetch_node_by_taxon_id($ncbi_taxon_id);
                    if (!$ncbi_taxon){
                         print "do_network(): can't find a taxon object for this: $ncbi_taxon_id\n";
                         $ncbi_shortname = "-";
                    }else{
                         $ncbi_shortname = $ncbi_taxon->short_name;
                    }
                    $new_row = join("\t", $id_in,  $ncbi_shortname,  $id_out);
               }else{
                    $new_row = join("\t", $id_in,$it,$id_out);
               }
          }
          print $out_data $new_row, "\n";
          print $new_row, "\n" unless($no_output);
     }
     $sth->finish();
     $dbh->disconnect();
     close $out_data;
     return 1;
}


=head2 do_attributes

 Usage     : $RC = Bio::Homology::InterologWalk::Networks::do_attributes(
                                                         registry      => $registry,
                                                         input_path    => $in_path,
                                                         output_path   => $out_path,
                                                         source_org    => $sourceorg,
                                                         label_type    => 'extname',
                                                         no_output      => 0
                                                         );
 Purpose   : This is needed to create a node attribute file to go with the .sif 
             network created by do_networks(). 
             For a definition of node attribute file, see
             http://cytoscape.wodaklab.org/wiki/Cytoscape_User_Manual#Node_and_Edge_Attributes 
             The routine associates, for each stable id in the sif file, a human-readable 
             gene name/description obtained from Ensembl
 Returns   : a boolean value for success/failure
 Argument  : -registry: ensembl registry object to connect to. Needed to retrieve up-to-date
              human readable gene names from Ensembl for the IDs in the input data file
             -input_path : path to input file. Input file for this subroutine is a tsv file 
              containing at least the fields INIT_ID and INTERACTOR.
              (output of get_backward_orthologies() or get_direct_interactions will work, 
              although output of remove_duplicate_rows() is recommended). 
             -output_path : where you want the routine to write the data. Data is a .noa 
              cytoscape file.
             -source_org : source organism name (eg: "Mus musculus")
             -(OPTIONAL) label_type: what kind of human readable string to employ. Options 
              are 'extname' (external name) and 'description'. Default is 'extname'
             -(OPTIONAL) no_output :  suppresses screen output. Used for clearer output 
              during test. Default is 0.
 Throws    : -
 Comment   : -

See Also   : L</do_network>

=cut

sub do_attributes{
     my %args = @_;
     my $registry         = $args{registry};
     my $in_path          = $args{input_path};
     my $out_path         = $args{output_path};
     my $sourceorg        = $args{source_org}; 
     my $label_type       = $args{label_type};
     my $no_output        = $args{no_output};
     
     if(!$registry){
          print("do_attributes(): no registry defined. Aborting..\n");
          return;
     }
     if(!$sourceorg){
          print("do_attributes(): no source organism specified. Aborting..\n");
          return;
     }
     $label_type = 'ext' if(!$label_type);
     
     my %seen = ();
     #MANAGE FILES
     my $dbh = Bio::Homology::InterologWalk::_setup_dbi_connection($in_path);
     open (my $out_data,  q{>}, $out_path) or croak("Unable to open $out_path : $!");
     #============
     
     my $gene_adaptor = $registry->get_adaptor($sourceorg, 'core', 'Gene');
     print $out_data "Ensembl_Name\n"; #attribute name
     
     my $query = "SELECT $FN_initid, $FN_interactor_id FROM int";
     my $sth = $dbh->prepare($query);
     $sth->execute() or die "Cannot execute: " . $sth->errstr();
     
     while (my $row = $sth->fetchrow_hashref) {
          my $id_in = $row->{$FN_initid};
          my $id_out = $row->{$FN_interactor_id};
          
          $seen{$id_in}  = 1 if(!$seen{$id_in});
          $seen{$id_out} = 1 if(!$seen{$id_out});
     }

     foreach my $stableid (sort keys %seen){
          my $label;
          my $genesArray = ();
          my $gene = $gene_adaptor->fetch_by_stable_id($stableid);
          
          if(!$gene){
               $gene = $gene_adaptor->fetch_by_display_label($stableid);
               if(!$gene){
                    $genesArray = $gene_adaptor->fetch_all_by_external_name($stableid);
                    if(scalar(@$genesArray) == 1){ 
                         $gene = @$genesArray[0];
                    }elsif(scalar(@$genesArray) > 1){
                         print("do_attributes(): warning - $stableid has ambiguous gene object in Ensemlb. Choosing the first.\n") unless($no_output);
                         $gene = @$genesArray[0];
                    }else{ #no gene object, leaving ID
                         print("do_attributes(): no gene object found for id: $stableid. Leaving as is.\n") unless($no_output);
                         $label = $stableid;
                         print "$stableid\t=\t$label\n" unless($no_output);
                         print $out_data $stableid, "\t", "=", "\t", $label, "\n";
                         next;
                    }
               }
          }
          
          if($label_type =~ /desc/){
               $label = $gene->description;
          }elsif($label_type =~ /ext/){
               $label = $gene->external_name || $gene->display_xref->display_id;
          }else{
               print("do_attributes(): label type unrecognised: $label_type. Aborting..\n");
               return;
          }
          $label = $stableid if(!$label);
          print "$stableid\t=\t$label\n" unless($no_output);
          print $out_data $stableid, "\t", "=", "\t", $label, "\n";
     }
     close $out_data;
     $dbh->disconnect();
     return 1;
}


##############################################
##############################################
package Bio::Homology::InterologWalk::Direct;
##############################################
##############################################
use String::Approx 'amatch';
use Carp qw(croak);

#
#_clean_string
#
# Usage     : How to use this function/method
# Purpose   : this is used to simply clean up the intact alternative id list or the intact properties list
# Returns   : a string of alternative ids or properties lacking the explanations in brackets
# Argument  : a string with the list of alternative ids or properties
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 

sub _clean_string{
     #I'm fine with this string, save for the recently introduced explanations in brackets. I
     #want to remove everything in brackets
     my ($data_string) = @_;
     my $clean_item;
     my @clean_item_vec;
     my $result;
     
     my @data_vec = split(/\|/, $data_string);
     
     foreach my $item (@data_vec){ #db:string
          if($item =~ /(\w+):(.+)\((.+)\)/){
               $clean_item = $1 . ":" . $2;
               push(@clean_item_vec, $clean_item);
          }else{
               push(@clean_item_vec, $item); #leave it unchanged
          }
     }
     $result = join('|', @clean_item_vec);
     return $result;
}

#
#_get_ens_id
#
# Usage     : How to use this function/method
# Purpose   : This is used to get an ensembl id from the properties field of a binary intact interaction.
#             It's necessary because the binary interaction returned by intact will be in uniprotkb format, while I want to pull out
#             a gene-level representation of the interactor in ensembl-recognisable format
# Returns   : either undefined or an ensembl id
# Argument  : aproperties string
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 

sub _get_ens_id{
     my ($prop_string) = @_;
     
     my @prop_vector = split(/\|/, $prop_string);
     
     foreach my $item (@prop_vector){
          if($item =~ /^(ensembl):(.+)/){
               return $2;
          }
     }
     return;
}

#
#_query_generic_identifier
#
# Usage     : How to use this function/method
# Purpose   : used to query some form of id towards ensembl and see if it recognises it
# Returns   : either undefined or the ensembl id for the input id
# Argument  : gene adaptor for the taxon considered, the id to query
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 
sub _query_generic_identifier{
     my ($adaptor, $ID) = @_;
     my $result;
     my $genesArray = ();
     my $gene;
     
     if($ID eq '-'){
          #print "_query_generic_identifier(): empty identifier, backing up..\n";
          return;
     }
     
     $gene = $adaptor->fetch_by_stable_id($ID);
     return $gene->stable_id if(defined $gene);
     
     $gene = $adaptor->fetch_by_display_label($ID);
     return $gene->stable_id if(defined $gene);
     
     $genesArray = $adaptor->fetch_all_by_external_name($ID);
     
     if(scalar(@$genesArray) == 1){ 
          $gene = @$genesArray[0];
          return $gene->stable_id if(defined $gene);
     }elsif(scalar(@$genesArray) == 0){
          return;
     }else{ #more than one gene object in genesArray
          foreach my $gn (@$genesArray){
               if (defined $gn){ #TODO I return the first defined one, improve this
                    $result = $gn->stable_id;
                    return $result;
               }
          }
     }
     return;
}


#_get_ensid_from_uniprotkb
#
# Usage     : How to use this function/method
# Purpose   : this is used when we have to resort to ensembl to recognise an external id: eg, when the intact data line does not contain
#             clues about the gene related to the uniprokb id they provide. It's based on a "backing up" procedure whereby I try to obtain
#             the id from progressively less reliable information provided by intact: (the accession number itself,
#             the protein name, the protein properties, the protein alternative ids)
# Returns   : either undefined or the ensembl id for the input id
# Argument  : gene adaptor for the taxon considered, id to query, name, properties string, aliases string (from intact)
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 

sub _get_ensid_from_uniprotkb{
     my ($adaptor, $input_ID, $input_Name, $input_Aliases, $input_Props) = @_;
     my $outputID;
     
     #I'm getting read of the ISOFORM : improve on this
     if($input_ID =~ /(\w+)(-\d{1})$/){
           $input_ID = $1;
     }
          
     #accession number?
     $outputID = _query_generic_identifier($adaptor, $input_ID);
     return $outputID if($outputID);
     #name?
     print "\t$input_ID: query by accession number ambiguous/empty, trying name..\n";
     $outputID = _query_generic_identifier($adaptor, $input_Name);
     return $outputID if($outputID);
     #properties?
     print "\t\t$input_ID: query by name ambiguous/empty, trying properties field..\n";
     my @properties = Bio::Homology::InterologWalk::_get_vector_from_string($input_Props);
     foreach my $property (@properties){
          $outputID = _query_generic_identifier($adaptor, $property);
          if($outputID){
               print "\t\t$outputID - Found by property: $property\n";
               return $outputID;
          }
     }
     #aliases?
     print "\t\t\t$input_ID: query by properties ambiguous/empty, trying aliases..\n";
     my @aliases = Bio::Homology::InterologWalk::_get_vector_from_string($input_Aliases); 
     foreach my $alias (@aliases){
          $outputID = _query_generic_identifier($adaptor, $alias);
          if($outputID){
               print "\t\t\t$outputID - Found by alias: $alias\n";
               return $outputID;
          }
     }

     return;
}


#
#_fuzzy_match
#
# Usage     : How to use this function/method
# Purpose   : sometimes, even when an ensembl-compatible output id is found in the intact datafile, it is in a different format that the one
#             we are dealing with. EG: FBgnxxxx and CGxxx are both recognised by ensembl. If I get a CGxxx I will have duplicates, and ideally
#             I want not only an output id in a format recognised by ensembl, but also one in the same format as the input one. This function
#             does a fuzzy match between a "signature" I've extracted from the input id and and the output id. Please adjust the parameters
#             to best match your needs according to the string::approx module guide on cpan:
#             http://search.cpan.org/~jhi/String-Approx-3.26/Approx.pm
# Returns   : return code for error/success
# Argument  : signature extracted from the input id (eg: ENSMUSG ) and output id
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 

sub _fuzzy_match {
     my ($sig, $id_to_test) = @_;
 
     return amatch($sig, [ # this array sets match options:
                              "i",    # match case-insensitively
                         ], $id_to_test);
}


#_is_primary_id
#
# Usage     : How to use this function/method
# Purpose   : ids obtained by Intact CAN be obsolete/secondary. This routine checks an output id found in an intact data entry against ensembl,
#             to see if it recognises it. If it recognises it, but the id is secondary, the routine will return the primary id.
# Returns   : undefined or the primary id for the id given as input
# Argument  : id to check against ensembl, gene adaptor
# Throws    : -
# Comment   : This is a sample subroutine header.
#           : -
#
#See Also   : 

sub _is_primary_id {
     my ($id_to_test, $gene_adaptor) = @_;
     my $candidateID;
     
     my $gene = $gene_adaptor->fetch_by_stable_id($id_to_test);
     return 1 if(defined $gene);
     
     my $genesArray = $gene_adaptor->fetch_all_by_external_name($id_to_test);
          
     if(scalar(@$genesArray) == 1){
          return 1;
     }elsif(scalar(@$genesArray) == 0){
          print("_is_primary_id(): $id_to_test found in Intact not recognised by ensembl. Backing up..\n");
          return;
     }else{
          foreach my $gene (@$genesArray){
               $candidateID = $gene->stable_id;
               return 1 if($candidateID eq $id_to_test);
          }
          return; #multiple gene objects and none matches our id...
     }
     return;
}


=head2 get_direct_interactions

 Usage     : $RC = Bio::Homology::InterologWalk::Direct::get_direct_interactions(
                                                                 registry        => $registry,
                                                                 source_org      => $sourceorg,
                                                                 input_path      => $in_path,
                                                                 output_path     => $out_path,
                                                                 url             => $url,
                                                                 check_ids       => 1,   
                                                                 no_spoke        => 1, 
                                                                 exp_only        => 1, 
                                                                 physical_only   => 1, 
                                                                 no_output       => 0 
                                                                 );
 Purpose   : this methods allows  to query the Intact database using the REST interface. 
             IntAct is the Molecular Interaction database at the European Bioinformatics 
             Institute (UK). The Intact project offers programmatic access to their data 
             through the PSICQUIC specification (see 
             http://code.google.com/p/psicquic/wiki/PsicquicSpecification).
             This routine is different and more complex than get_interactions() from the 
             main module. This one is meant to query intact directly with the ids provided 
             by the user: no intermediate orthologues from ensembl are collected.
             The bulk of the script is used for the following reason: each query to intact 
             through psicquic returns a data entry including a binary protein interaction, 
             and the the two ids returned are uniprotkb or other protein ids. 
             We need to
                a- convert both to a format recognised by ensembl
                b- identify which of the two corresponds to our initial id
                c- convert the other one to ensembl and store it in the file
             This conversion is not trivial as the possibility of ambiguities/errors/wrong 
             matches between ensembl gene representations and uniprot protein representations 
             is high.
 Returns   : return code for error/success 
 Argument  : -registry: registry object to connect to Ensembl
             -source_org : source organism name (eg: "Mus musculus")
             -input_path : path to input file. Input file MUST be a text file with one entry 
              per row, each entry containing an up-to-date
              gene ID recognised by the Ensembl consortium (http://www.ensembl.org/) followed 
              by a new line char.
             -output_path : where you want the routine to write the data. Data is in TSV format.
             -url : url for the REST service to query (currently only EBI Intact PSICQUIC Rest)
             -(OPTIONAL) check_ids : if true, every interactor id found in intact data will 
              be double checked against ensembl.
              this is useful because intact dbs sometimes contain obsolete versions of some 
              ids. However chosing true will significantly slow down the processing
             -(OPTIONAL) no_spoke: if set, interactions obtained from the expansion of 
               complexes through the SPOKE method 
              (see http://nar.oxfordjournals.org/cgi/content/full/38/suppl_1/D525)
              will be ignored
             -(OPTIONAL) exp_only: if set, only interactions whose MITAB25 field 
              "Interaction Detection Method" 
              (MI:0001 in the PSI-MI controlled vocabulary) is at least "experimental 
              interaction detection" 
              (MI:0045 in the PSI-MI controlled vocabulary) will be retained. I.e. if set, 
              this flag only allows 
              experimentally detected interactions to be retained and stored in the data file
             -(OPTIONAL) physical_only: if set, only interactions whose MITAB25 field 
              "Interaction Type" 
              (MI:0190 in the PSI-MI controlled vocabulary) is at least "physical association" 
              (MI:0915 in the PSI-MI controlled vocabulary) will be retained. I.e. 
              if set, this flag only allows 
              physically associated PPIs to be retained and stored in the data file: 
              colocalizations and genetic interactions will be discarded
             -(OPTIONAL) no_output :  suppresses screen output. Used for clearer output 
              during test. Default is 0.
 Throws    : -
 Comment   : -

See Also   : L</get_interactions>

=cut

sub get_direct_interactions{
     my %args = @_;
     
     my $registry        = $args{registry};
     my $sourceorg       = $args{source_org};
     my $in_path         = $args{input_path};
     my $out_path        = $args{output_path};
     my $url             = $args{url};
     my $check_ids       = $args{check_ids};
     my $no_spokes       = $args{no_spoke};
     my $exp_only        = $args{exp_only};
     my $physical_only   = $args{physical_only};
     my $no_output       = $args{no_output};
     
     if(!$registry){
          print("get_direct_interactions(): no registry defined. Aborting..\n");
          return;
     }
     if(!$sourceorg){
          print("get_direct_interactions(): no source organism specified. Aborting..\n");
          return;
     }
     if(!$url){
          print("get_direct_interactions(): no PSICQUIC url specified. Aborting..\n");
          return;
     }

     #MANAGE FILES
     open (my $in_data,  q{<}, $in_path) or croak("Unable to open $in_path : $!");
     open (my $out_data,  q{>}, $out_path) or croak("Unable to open $out_path : $!");
     #============
     
     my $client = REST::Client->new();
     my $gene_adaptor = $registry->get_adaptor($sourceorg, 'core', 'Gene'); 
     
     my $atleast_one_entry;
     
     my $DF_interaction_id;
     my $DF_acc_numb_a;  
     my $DF_acc_numb_b;
     my $DF_alt_id_a;
     my $DF_alt_id_b;
     my $DF_name_a; 
     my $DF_name_b;
     my $DF_taxon_a;
     my $DF_taxon_b;
     my $DF_props_a;
     my $DF_props_b;
     my $DF_pub;
     my $DF_int_type;    
     my $DF_det_method;
     my $DF_exp_method;
     
     #PISCQUIC GLOBAL search string
     my $glob_search_string = "search/query/";
     #Header
     print $out_data $HEADER_DIRECT, "\n";
     
     my $options = Bio::Homology::InterologWalk::_build_query_options($no_spokes, $exp_only, $physical_only);
     #eg " AND type:\"physical*\" AND detmethod:\"experimental*\" AND NOT expansion:\"spoke\"";

     my $missed = 0;
     while (<$in_data>){
          my ($ID) = $_;
          chomp $ID;
          next if ($ID eq '');
          
          my $idsignature;
          #get a "signature" to spot the kind of id we are dealing with.
          #current solution involves getting all the letters starting from the beginning, if there's at least two.
          #otherwise get the initial three characters whatever they are, and then do a fuzzy regex matching using string::approx
          #this will be needed in order to be sure to get the same kind of id back.
          #eg "IPR006259" ----> "IPR"
          if($ID =~ /^([a-z]{2,})(.+)/i){
               $idsignature = $1;
          }else{
               $idsignature = substr($ID, 0, 1) . substr($ID, 1, 1) . substr($ID, 1, 1);
          }
          
          print "$ID: Querying IntAct web service for binary interactions.." unless($no_output);
          my $request = $url . $glob_search_string .  $ID . $options ;
          $client->GET($request);
          print "(", $client->responseCode(), ")" unless($no_output);
          my $responseContent = $client->responseContent();
          if(!$responseContent){
               print("..nothing..\n") unless($no_output);
               next;
          }
          $atleast_one_entry = 1;
          my @responsetoparse = split(/\n/,$responseContent);
          my $interactionsRetrieved = scalar @responsetoparse;
          print "..Interactions found: ", $interactionsRetrieved, "\n" unless($no_output);
          
          foreach my $intactInteraction (@responsetoparse){
               my @MITABDataRow = split("\t",$intactInteraction);
               
               #qui devi sfruttare le regole del MITAB
               $DF_interaction_id = Bio::Homology::InterologWalk::_get_intact_id($MITABDataRow[13]);
               next if(!$DF_interaction_id);
               
               $DF_acc_numb_a = Bio::Homology::InterologWalk::_get_interactor_uniprot_id($MITABDataRow[0]);
               $DF_acc_numb_b = Bio::Homology::InterologWalk::_get_interactor_uniprot_id($MITABDataRow[1]);
               $DF_alt_id_a   = _clean_string($MITABDataRow[2]);
               $DF_alt_id_b   = _clean_string($MITABDataRow[3]);
               $DF_name_a     = Bio::Homology::InterologWalk::_get_interactor_name($MITABDataRow[4]);
               $DF_name_b     = Bio::Homology::InterologWalk::_get_interactor_name($MITABDataRow[5]);
               $DF_taxon_a    = Bio::Homology::InterologWalk::_get_interactor_taxon($MITABDataRow[9]);
               $DF_taxon_b    = Bio::Homology::InterologWalk::_get_interactor_taxon($MITABDataRow[10]);
               
               next if($DF_taxon_a ne $DF_taxon_b);
               
               $DF_pub        = $MITABDataRow[8];
               $DF_int_type   = $MITABDataRow[11];
               $DF_det_method = $MITABDataRow[6];
               $DF_exp_method = $MITABDataRow[24];
               $DF_props_a    = _clean_string($MITABDataRow[19]);
               $DF_props_b    = _clean_string($MITABDataRow[20]);
               
               my $candidate;
               my $ID_OUT; #the object of our search
               
               #first of all let's consider the (rare) case in which the two ids are in the same format as the initial ids:
               if( ($ID eq $DF_acc_numb_a) or ($ID  eq $DF_acc_numb_b)  ){
                    $candidate = $DF_acc_numb_b if($ID eq $DF_acc_numb_a);
                    $candidate = $DF_acc_numb_a if($ID eq $DF_acc_numb_b); 
                    
                    $ID_OUT = _query_generic_identifier($gene_adaptor,$candidate);
                    if($ID_OUT){
                         print("Interaction ($DF_interaction_id): $ID <--> $ID_OUT\n") unless($no_output);
               
                         my $fullDataRow = join("\t",$ID,$DF_interaction_id,
                                             $DF_acc_numb_a, $DF_acc_numb_b,
                                             $DF_alt_id_a, $DF_alt_id_b,
                                             $DF_name_a, $DF_name_b,
                                             $DF_taxon_a, $DF_taxon_b,
                                             $DF_pub, $DF_int_type,$DF_det_method,
                                             $DF_exp_method,$ID_OUT);
                         print $out_data $fullDataRow, "\n";
                         next;
                    }
               }

               #I need ensembl ids for each interactor pair. They must be of the same kind of those present in the initial set
               #However, I only get UniprotKB accession numbers from Intact. There are several ways to obtain
               #ensembl ids from those. 
               #We could use the CG alternative id for each UNIPROTKB stored in "alt_id_A/alt_id_B" and then obtain the fbid through ensembl.
               #Intact MIGHT provide the fb ids of the interactors in the properties field (19 and 20). However this will not always happen
               #we can also query ensemble with the UNIPROTKB and see if we get back a fbid stable_id. The second way is probably better, however, ensembl
               #won't recognise isoforms (eg queries in the form "Q24312-2") which we have to clean out (Q24312) therefore losing some information.
               #CURRENT SOLUTION:
               #1)look for ids in the "properties" field, if none
               #2)look for ids in the "aliases" field, if none
               #3)query ensembl with the uniprotkb (slower and possibility of redundancy and ambiguities)
               
               #1
               my $target_PropsRow;my $target_AccNumb;my $target_Aliases;my $target_Name;
               #TODO autointeractions
               
               #if at least one of them is recognisable, I know I'll have to analyse the other
               if( ($DF_props_a =~ /($ID)/) or ($DF_props_b =~ /($ID)/)  ){
                    if( $DF_props_a =~ /($ID)/){
                         $target_PropsRow = $DF_props_b;
                         $target_AccNumb = $DF_acc_numb_b;
                         $target_Aliases = $DF_alt_id_b;
                         $target_Name = $DF_name_b;
                    }
                    if( $DF_props_b =~ /($ID)/){
                         $target_PropsRow = $DF_props_a;
                         $target_AccNumb = $DF_acc_numb_a;
                         $target_Aliases = $DF_alt_id_a;
                         $target_Name = $DF_name_a;
                    }
                    
                    $ID_OUT = _get_ens_id($target_PropsRow);
                    
                    #fuzzy string matching: this ID_OUT should be of the same kind as the original id. Does it feature the same initial
                    #id signature or something very close?
                    if(!$ID_OUT){
                         $ID_OUT = _get_ensid_from_uniprotkb($gene_adaptor, $target_AccNumb, $target_Name, $target_Aliases, $target_PropsRow);
                    }
                    
                    if($check_ids and (!_is_primary_id($ID_OUT, $gene_adaptor)) ){
                         $ID_OUT = _get_ensid_from_uniprotkb($gene_adaptor, $target_AccNumb, $target_Name, $target_Aliases, $target_PropsRow);
                    }
                         
                    if(!_fuzzy_match($idsignature, $ID_OUT)){
                         $ID_OUT = _get_ensid_from_uniprotkb($gene_adaptor, $target_AccNumb, $target_Name, $target_Aliases, $target_PropsRow);
                    }
                    

                    if($ID_OUT){
                         print("Interaction ($DF_interaction_id): $ID <--> $ID_OUT\n") unless($no_output);
               
                         my $fullDataRow = join("\t",$ID,$DF_interaction_id,
                                             $DF_acc_numb_a, $DF_acc_numb_b,
                                             $DF_alt_id_a, $DF_alt_id_b,
                                             $DF_name_a, $DF_name_b,
                                             $DF_taxon_a, $DF_taxon_b,
                                             $DF_pub, $DF_int_type,$DF_det_method,
                                             $DF_exp_method,$ID_OUT);
                         print $out_data $fullDataRow, "\n";
                         next;
                    }
               }
               
               #If none of the two was recognisable, I'll have to process them both
               print("Nothing found in Intact Data...Backing up \n");

               my $interactoridA = _get_ensid_from_uniprotkb($gene_adaptor, $DF_acc_numb_a, $DF_name_a, $DF_alt_id_a, $DF_props_a);
               my $interactoridB = _get_ensid_from_uniprotkb($gene_adaptor, $DF_acc_numb_b, $DF_name_b, $DF_alt_id_b, $DF_props_b);
               
               unless($interactoridA and $interactoridB){
                    print("One of the two identifier could not be retrieved. Skipping..\n");
                    $missed += 1;
                    next;
               }

               if($interactoridA eq $ID){
                    $ID_OUT = $interactoridB;
               }elsif($ID eq $interactoridB){
                    $ID_OUT = $interactoridA;
               }else{
                    print("\nWARNING id mismatch between $ID, $interactoridA, $interactoridB. Skipping..\n");
                    $missed += 1;
                    next;
               }

               print("Interaction ($DF_interaction_id): $ID <--> $ID_OUT\n") unless($no_output);
               my $fullDataRow = join("\t",$ID,$DF_interaction_id,
                                             $DF_acc_numb_a, $DF_acc_numb_b,
                                             $DF_alt_id_a, $DF_alt_id_b,
                                             $DF_name_a, $DF_name_b,
                                             $DF_taxon_a, $DF_taxon_b,
                                             $DF_pub, $DF_int_type,$DF_det_method,
                                             $DF_exp_method,$ID_OUT);
               print $out_data $fullDataRow, "\n";
          }
     }
     print("Missed: $missed\n") unless($no_output);
     close $in_data;
     if(!$atleast_one_entry){
        unlink($out_path);
        print("\n**No interactions found. Exiting..**\n");
        return;
     }
     close $out_data;
     return 1;
}


##################secondary pod documentation begins##############


=head1 BUGS AND LIMITATIONS

This is B<ALPHA> software. There will be bugs. The interface may change. Please be careful. Do not rely on it for anything mission-critical.


Please report any bugs you find, bug reports and any other feedback are most welcome. 

-Currently only the  EBI Intact DB is available for PPI retrieval. This will be expanded to account for  all available PSICQUIC-compliant PPI dbs transparently. 
This includes MINT, STRING, BioGrid and many more. For a full list of compliant DBs and for the status of the PSICQUIC service, check

http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=STATUS 



=head1 AUTHOR

Giuseppe Gallone  E<lt>ggallone@cpan.orgE<gt>

CPAN ID: GGALLONE

University of Edinburgh    

=head1 LICENSE AND COPYRIGHT

Bio::Homology::InterologWalk is Copyright (c) 2010 Giuseppe Gallone
All rights reserved.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.1 or,
at your option, any later version of Perl 5 you may have available.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

=head1 ACKNOWLEDGEMENTS

The author would like to thank the following individuals and organisations for their invaluable support and priceless suggestions.

Andrew Jarman, Ian Simpson, Douglas Armstrong and all the Jarman Lab,  University of Edinburgh

Javier Herrero, Albert Vilella, Andy Yates, Glenn Proctor, Michael Han, Gautier Koscielny and all the Ensembl Team

Samuel Kerrien & Bruno Aranda and all the EBI-Intact Team

Dave Messina, BioPerl list

=cut

#################### secondary pod documentation end ###################

1;
# The preceding line will help the module return a true value

