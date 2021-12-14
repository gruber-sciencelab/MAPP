###############################################################################
#
#   Select only poly(A) sites which may be categorized as proximal/distal.
#
#   AUTHOR: Ralf_Schmidt
#   MODIFIED BY: Mihaela Zavolan
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: mihaela.zavolan@unibas.ch
#   CREATED: XX-XX-XXXX
#   MODIFIED: 01-02-2021
#   LICENSE: Apache_2.0
#
###############################################################################

use strict;
use warnings;
use Getopt::Long;

my $annofile;
my $help;
my $minLevel;
my $nonredundant;
my $offset = 0;
my $locusExtension = 100;
my @types;
my $type_id;
GetOptions(
  "annotation=s" => \$annofile,
  "minLevel=f"   => \$minLevel,
  "type_id=s"      => \$type_id, #ENSEMBL (biotype)/Gencode (type) field in gtf
  "type=s"       => \@types, #types (e.g. protein-coding, lncRNA, etc.)
  "help"         => \$help,
  "h"            => \$help,
  "nonredundant" => \$nonredundant,
  "offset=i"     => \$offset, #3' extension (in nts) of annotated transcript
  "locusExtension=i"  => \$locusExtension #up/down extension of transcript locus to catch reads that start/end beyond transcript boundaries
);

# process options
my $showhelp = 0;

$showhelp = 1 if ( not defined $annofile );
$showhelp = 1 if ( not defined $ARGV[0] );

if ($showhelp) {
  print STDERR "Usage: $0 --annotation=gencode.gtf.gz in.bed > out.bed\n";
  exit 2;
}

if ( not -e $annofile ) {
  print STDERR "[ERROR] File '$annofile' does not exist.\n";
  exit 2;
}
if ( -d $annofile ) {
  print STDERR "[ERROR] '$annofile' is not a file.\n";
  exit 2;
}

if ( not -e $ARGV[0] ) {
  print STDERR "[ERROR] File '$ARGV[0]' does not exist.\n";
  exit 2;
}
if ( -d $ARGV[0] ) {
  print STDERR "[ERROR] '$ARGV[0]' is not a file.\n";
  exit 2;
}

if ( $#types > -1 && not defined $type_id) {
  print STDERR "[ERROR] type filter is given but the name of the type_id (e.g. 'transcript_biotype') is missing.\n";
  exit 2;
}

if ( not defined $type_id) {
  print STDERR "[WARNING] No type id given (--type_id). All exons are considered\n";
}

# check if gunzip is installed
my $gunzip    = `which gunzip`;
my $no_gunzip = 0;
if ( $gunzip =~ m/no\sgunzip\sin/ ) {
  $no_gunzip = 1;
  print STDERR "[INFO] Unzipping of .gz files not possible.\n";
}

# make types hash
my $checktypes = 0;
my %TYPES      = ();
if ( $#types > -1 ) {
  foreach my $type (@types) {
    $TYPES{$type} = 1;
  }
  $checktypes = 1;
}

my $fh;
if ( $annofile =~ m/\.gz$/ ) {
  if ( $no_gunzip == 1 ) {
    print STDERR "[INFO] Unzipping is not possible.\n";
    exit 1;
  }
  open( $fh, "gunzip -c $annofile |" );
} else {
  open( $fh, $annofile );
}

# process the annotation file

my %transcripts     = ();
my %genes           = ();
my %T2G             = ();
my %geneCoordinates = ();
while ( my $line = <$fh> ) {
  next if ( $line =~ m/^#/ );
  chomp $line;
  my @F = split( /\t/, $line );

  #####
  # only consider exon entries
  next if ( $F[2] ne 'exon' );
  my $gene_id         = '';
  my $transcript_id   = '';
  my $transcript_type = '';

  ####
  # get IDS for gene and transcript, gene type of the current exon
  if ( $F[8] =~ m/gene_id\s+"([^"]+)";/ ) {
    $gene_id = $1;
    $gene_id =~ s/\.\d+$//;
  }
  if( defined $type_id){
    if ( $F[8] =~ m/$type_id\s+"([^"]+)";/ ) {
      $transcript_type = $1;
    }
  }
  
  if ( $F[8] =~ m/transcript_id\s+"([^"]+)";/ ) {
    $transcript_id = $1;
    $transcript_id =~ s/\.\d+$//;
  }

  # only process this entry if type matches what was given
  # as argument
  if ( $checktypes == 1 ) {
    next if ($transcript_type eq '');
    next if ( not defined $TYPES{$transcript_type} );
  }

  #####
  # error handling
  #####
  if ( $transcript_id eq '' ) {
    print STDERR "[ERROR] Did not succeed in inferring transcript_id.\n";
    print STDERR "[ERROR] Offending line: '$line'\n";
    exit;
  }

  if ( $gene_id eq '' ) {
    print STDERR "[ERROR] Did not succeed in inferring gene_id.\n";
    print STDERR "[ERROR] Offending line: '$line'\n";
    exit;
  }

  if ( $gene_id ne '' ) {
    $genes{$gene_id} = [] if ( not defined $genes{$gene_id} );

    # save each transcript_id for the corresponding gene_id
    push @{ $genes{$gene_id} }, $transcript_id;

    # save the gene_id for every transcript_id
    $T2G{$transcript_id} = $gene_id;
  }

  # save each exon of each transcript
  # transcripts hash contains lists of exons for each transcript
  $transcripts{$transcript_id} = [] if ( not defined $transcripts{$transcript_id} );

  # save chr, start, stop, and strand of the current exon for each transcript
  # ensure that the exons are included 5' -> 3' in the direction of transcription
  if( $#{ $transcripts{$transcript_id} } == -1) {
    # first exon that we see for this transcript
    push @{ $transcripts{$transcript_id} }, [ $F[0], $F[3], $F[4], $F[6] ];
  } elsif ( $F[6] eq "-" && $F[3] > $transcripts{$transcript_id}->[$#{$transcripts{$transcript_id}}]->[2] ) {
    # negative strand transcript
    # pre-pend the current exon (obviously assumes that the exons come in sorted order
    unshift @{ $transcripts{$transcript_id} }, [ $F[0], $F[3], $F[4], $F[6] ];
  } else {
    # positive strand transcript
    push @{ $transcripts{$transcript_id} }, [ $F[0], $F[3], $F[4], $F[6] ];
  }
}
close($fh);

#ensure that the exons are in sorted order 5' to 3' for both plus and minus strand exons
foreach my $transcript_id (keys %transcripts) {
    if($transcripts{$transcript_id}->[0]->[3] eq '+') {
	@{ $transcripts{$transcript_id} } = sort {sort_by_start($a, $b)} @{ $transcripts{$transcript_id} }
    }
    else {
	@{ $transcripts{$transcript_id} } = sort {sort_by_start($b, $a)} @{ $transcripts{$transcript_id} } 
    }
}

my $n = scalar keys %transcripts;
print STDERR "[INFO] Collected $n transcript entries.\n";

$n = scalar keys %genes;
print STDERR "[INFO] Collected $n gene entries.\n";

########
# NOW, READ THE POLY(A)CLUSTERS FILE
########

my $fhCL;
if ( $ARGV[0] =~ m/\.gz$/ ) {
  if ( $no_gunzip == 1 ) {
    print STDERR "[INFO] Unzipping is not possible.\n";
    exit 1;
  }
  open( $fhCL, "gunzip -c $ARGV[0] |" );
} else {
  open( $fhCL, $ARGV[0] );
}

my %clusters = ();
my %sortIdx  = ();
my $sortN    = 0;
$n = 0;
while ( my $line = <$fhCL> ) {
  next if ( $line =~ m/^#/ );
  chomp($line);
  my @F = split( /\t/, $line );

  # check validity of the cluster
  _checkBED( \@F, $ARGV[0] );
  # Allow int and float for the score
  if ( $F[4] !~ m/^\d+$/ ) {
    print STDERR "[ERROR] Wrong support level annotation.\n";
    print STDERR "[ERROR] Offending line: '$line'.\n";
    exit;
  }

  # only consider clusters above the provided min support level (score)
  next if ( $F[4] < $minLevel );

  $n++;
  my $key = "$F[0]:$F[5]"; #key is chr:str
  if ( not defined $sortIdx{$key} ) {
    $sortN++;
    $sortIdx{$key} = $sortN;
  }

  # for each cluster save: [ start, stop, ID, score ]
  push @{ $clusters{$key} }, [ $F[1], $F[2], $F[3], $F[4] ];
}
close($fhCL);
print STDERR "[INFO] Collected $n PAS clusters.\n";

#######
# assign the clusters to the transcript
#######

# count the assigned clusters
my %assignedCL = ();

my %transcriptsCL = ();
$n = 0;
print STDERR "[INFO] ";
foreach my $t ( keys %transcripts ) {

  my @PAS   = ();
  my $n_PAS = 0;

  # for each exon the following array is available: [  chr, start, stop, strand ]

  #extend 3' end of most 3' exon by $offset
  # expects the gtf file to be in increasing genomic position order
  my $chr    = $transcripts{$t}->[0]->[0];
  my $strand = $transcripts{$t}->[0]->[3];
  my $key    = "$chr:$strand";

  if ( $strand eq '+' ) {
    $transcripts{$t}->[ $#{ $transcripts{$t} } ]->[2] += $offset;
  } else {
    $transcripts{$t}->[ $#{ $transcripts{$t} } ]->[1] -= $offset;
  }

  my @exons = @{ $transcripts{$t} };

  # determine min and max (==genomic region of the transcript)
  my $min = 10e12;
  my $max = 0;

  #helper saves the label of each genomic position
  #which is the number of the exon to which the position belongs
  my %helper = ();
  foreach my $idx ( 0 .. $#exons ) {
    my $start = $exons[$idx]->[1];
    my $end   = $exons[$idx]->[2];
    $min = $start if ( $start < $min );
    $min = $end   if ( $end < $min );
    $max = $start if ( $start > $max );
    $max = $end   if ( $end > $max );
    foreach my $pos ( $start .. $end ) {
      $helper{$pos} = $idx;
    }

    # initialize an array of PASes for each exon
    push @PAS, [];
  }

  # start and stop of the transcript is extended by X nts
  # to ensure that reads which start outside of the first exon
  # can be considered as well as long as the end position of the
  # read falls into a desired feature
  $min = $min - $locusExtension;
  $max = $max + $locusExtension;

  # now collect all polyA clusters that are located within these exons
  foreach my $cluster ( @{ $clusters{$key} } ) {

    # for each cluster saved: [ start, stop, name, score ]
    my $cl_start = $cluster->[0];
    my $cl_end   = $cluster->[1];

    # check if the current cluster is in the region of the current transcript
    next if ( $cl_start < $min and $cl_end < $min );
    next if ( $cl_start > $max and $cl_end > $max );

    my %tmp = ();
    foreach my $pos ( $cl_start .. $cl_end ) {

      # increment counter for an exon if a position of the cluster is in this exon
	if ( defined $helper{$pos} ) {
	    if( defined $tmp{ $helper{$pos} }) {
		$tmp{ $helper{$pos} }++;
	    }
	    else {
		$tmp{ $helper{$pos} } = 1;
	    }
	}
    }

    my @K = keys %tmp;

    #cluster is only considered if it belongs to exactly one exon
    if ( $#K == 0 ) {
      $n_PAS++;

      # push the current cluster at the apropriate position
      #of the @PAS array that saves the PAS per exon for each transcript
      # ($PAS[0] contains the pAs of the first exon)
      push @{ $PAS[ $K[0] ] },
        [ $chr, $cluster->[0], $cluster->[1], $cluster->[2], $cluster->[3], $strand ];
    }
    else {
	    print STDERR "Cluster " . $cluster->[2] . " is part of multiple exons\n";
    }
  }

  # write to the log file each 5000 transcripts
  $n++;
  print STDERR "." if ( $n % 5000 == 0 );

# if at least one exon contains one cluster, save the list of exons with their corresponding clusters
  if ( $n_PAS > 0 ) {
    $transcriptsCL{$t} = \@PAS;

    #foreach my $idx ( 0 .. $#PAS ) {
    #  foreach my $entry ( @{$PAS[$idx]} ) {
    #	print "$idx\t",join("\t",@{$entry}),"\n";
    # }
    #}
    #print "\n";
  }
}

print STDERR "\n[INFO] Clusters to transcript assignment done.\n";

# what we have so far: %transcriptsCL, that saves an array for each transcript
# each array position corresponds to one exon (ordered from 5' to 3')
# the exon-arrays contain references to pAs located on that exon

my %out = ();

foreach my $t ( keys %transcriptsCL ) {
  #data structure with PASes corresponding to individual exons in the transcript $t
  #it's an array of lists, index in the array is the exon index in the transcript
  my @exons = @{ $transcriptsCL{$t} };

  # each entry contains: [ $chr, $cluster->[0], $cluster->[1], $cluster->[2], $cluster->[3], $strand ]
  # the cluster->[] entries correspond to: [0]:start, [1]:stop, [2]:name, [3]:score ]
  foreach my $idx ( 0 .. $#exons ) {
    my @entries = @{ $exons[$idx] };

    #current exon contains more than one polyAsite
    if ( $#entries > 0 ) {

      # sort entries from 5' to 3'
      my $chr    = $entries[0]->[0];
      my $strand = $entries[0]->[5];
      if ( $strand eq '+' ) {
        @entries = sort { $a->[1] <=> $b->[1] } @entries;
      } else {
        @entries = sort { $b->[1] <=> $a->[1] } @entries;
      }

      my $ln = $#entries + 1; #nr PASes in this exon
      foreach my $entryIDX ( 0 .. $#entries ) {
        my $l = $entryIDX + 1;

        # save number of exon (considering the strand)
        # i.e. most 3' exon is first for plus strand and last for minus strand
        my $kn = $#exons + 1;
        my $k;

        # new gtf convention: exons are already sorted 5' -> 3' in the gtf
        # not according to coordinates anymore!
        $k = $idx + 1;
        my $G = ( defined $T2G{$t} ) ? $T2G{$t} : '.';
        $out{"$chr:$strand"} = [] if ( not defined $out{"$chr:$strand"} );

        # as last entry, save how many exons follow for this transcript after the one
        # processed currently
        my $distance_to_3p = $kn - $k;
	      # create output 
	      # PASinfo PASIndexInExon #nrPASesInExon "transcriptID:exonInTranscript:totalExons" exonStart exonEnd geneID exonsLeftDownstream
        push @{ $out{"$chr:$strand"} }, [
          @{ $entries[$entryIDX] },
          $l, $ln,
          "$t:$k:$kn:" . $transcripts{$t}->[$idx]->[1] . ":" . $transcripts{$t}->[$idx]->[2],
          $G, $distance_to_3p, "nan"
          ];
      }
    }
  }
}

#if one exon is included several times only the entry with the most polyAsites is kept
if ($nonredundant) {

  foreach my $key ( keys %out ) {
    my %G = ();

    #iterate over all pAs of this key
    foreach my $entry ( @{ $out{$key} } ) {

      #use exonID as key for %G
      $G{ $entry->[8] } = [] if ( not defined $G{ $entry->[8] } );
      push @{ $G{ $entry->[8] } }, $entry;
    }

    my @new  = ();
    my %used = ();

    #iterate over all exonIDs of current chr:strand (sorted after decreasing number of tandem pAs
    # and after position of the exon within th transcript (prefer more distal exons)
    # and after the transcript support level (prioritize transcripts with lower support level) )
    foreach my $groupid (
      sort { $G{$b}->[0]->[7] <=> $G{$a}->[0]->[7] or $G{$a}->[0]->[10] <=> $G{$b}->[0]->[10] }
      keys %G
      ) {

      my $keep = 1;

      #skip current exon if one of the pAs was already used for another exon
      foreach my $entry ( @{ $G{$groupid} } ) {
        $keep = 0 if ( defined $used{ $entry->[3] } );
      }
      if ( $keep == 1 ) {
        foreach my $entry ( @{ $G{$groupid} } ) {

          # delete last two element of the entry since it was only necessary for sorting purposes
          pop @{$entry};
          pop @{$entry};
          push @new, $entry;
          $used{ $entry->[3] } = 1;
        }
      }
    }
    $out{$key} = \@new;
  }
}

#print out all salected PASes
foreach my $key ( sort { $sortIdx{$a} <=> $sortIdx{$b} } keys %sortIdx ) {
    if ( not defined $out{$key} ) {
	print STDERR "$key discarded in the non-redundant exon selection step\n";
    }
    else {
	my @O = @{ $out{$key} };
	foreach my $entry ( sort { $a->[9] cmp $b->[9] } @O ) {
	    print join( "\t", @{$entry} ), "\n";
	}
    }
}

sub _checkBED {
  my @F    = @{ $_[0] };
  my $file = $_[1];
  if ( not defined $F[5] ) {
    print STDERR "[ERROR] Content of $file does not seem to be in BED format.\n";
    print STDERR "[ERROR] Offending line: '", join( "\t", @F ), "'.\n";
    exit;
  } else {
    if ( $F[5] !~ m/^(\+|-)$/ ) {
      print STDERR "[ERROR] Content of $file does not seem to be in BED format.\n";
      print STDERR "[ERROR] Offending line: '", join( "\t", @F ), "'.\n";
      exit;
    }
  }
}

sub sort_by_start {
    my ($listref1, $listref2) = @_;
    return $listref1->[1] <=> $listref2->[1];
}
