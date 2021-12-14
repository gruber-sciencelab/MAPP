###############################################################################
#
#   Intersect tandem poly(A) sites of terminal exons with the annotation.
#   Only retain exons that can be associated to a gene unambiguously.
#
#   AUTHOR: Ralf_Schmidt
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: mihaela.zavolan@unibas.ch
#   CREATED: XX-XX-XXXX
#   LICENSE: Apache_2.0
#
###############################################################################

use strict;
use warnings;
use Getopt::Long;

my $tandemPAS;
my $downstream_extend = 0;

GetOptions( "tandemPAS=s" => \$tandemPAS, "downstream_extend=i" => \$downstream_extend );

#gtf file is $ARGV[0] file with PASes is passed in as argument
if ( not defined $tandemPAS or not defined $ARGV[0] ) {
  print STDERR "[ERROR] Usage: perl $0 --tandemPAS=tandem.pas.bed --downstream_extend=100 annotation.ensembl.gtf.gz\n";
  exit;
}

if ( !-e $ARGV[0] ) {
  print STDERR "[ERROR] Input file not found\n";
  exit(2);
}

my %pos           = ();
my %exons         = ();
my $exon_number   = 1;
my %exon_order    = ();
my %E2G           = ();
my %pos_per_exons = ();
my %processed     = ();

# stats
my $deleted_cnt = 0;

# first:
# delete overlapping exons
# store location of each not overlapping exon position wise
open( F, "< $tandemPAS" ) or die "Can't open $tandemPAS\n";
while (<F>) {
  chomp;
  my @F = split("\t");

  #$F[8] is the exonID (transcript exonNr exonsInTranscript start end)
  $exons{ $F[8] } = [] if (not defined $exons{ $F[8] } );
  push @{ $exons{ $F[8] } }, \@F;

  next if ( defined $processed{ $F[8] } ); #if exon processed, go to next line in input file
  $exon_order{ $F[8] } = $exon_number; #keep track of all unique exons?
  $exon_number++;

  #$key is chromosome
  my $key = "$F[0]";
  my $strand = $F[5];

  ( undef, undef, undef, my $ex_start, my $ex_end ) = split( ":", $F[8] );
  if( $strand eq "+" ) {
    $ex_end += $downstream_extend;
  } else {
    $ex_start -= $downstream_extend;
  }

  my $to_delete = 0;
  $pos_per_exons{ $F[8] } = [];

  # save the exon-gene-association
  $E2G{ $F[8] } = $F[9];
  foreach my $i ( $ex_start .. $ex_end ) {
    $pos{$key} = {} if ( not defined $pos{$key} );
    push @{ $pos_per_exons{ $F[8] } }, $i;
    if ( defined $pos{$key}->{$i} ) {
      if ( defined $pos_per_exons{ $pos{$key}->{$i} } ) {
        delete $pos_per_exons{ $pos{$key}->{$i} };
        delete $E2G{ $pos{$key}->{$i} };
	$deleted_cnt++;
        print STDERR "[DEL] Deleted ", $pos{$key}->{$i}, " because of overlap with $F[8]\n";
      }
      $to_delete = 1;
    } else {
      $pos{$key}->{$i} = $F[8];
    }
  }
  if ($to_delete) {
    delete $pos_per_exons{ $F[8] };
    delete $E2G{ $F[8] };
    $deleted_cnt++;
    print STDERR "[DEL] Also deleted $F[8]\n";
  }
  $processed{ $F[8] } = 1;
}
close(F);

print STDERR "[INFO] Deleted exons due to overlap with any other terminal exon in the list: $deleted_cnt\n";
my $remaining_exons = scalar keys %pos_per_exons;
print STDERR "[INFO] Remaining exons: $remaining_exons\n";

my %new_exons = ();
%pos = ();
# re-initialize the pos-hash with the coordinates of
# all remaining not-overlapping exons
foreach my $ex (keys %pos_per_exons) {
  #$ex is exonID
  #$exons{$ex} is the reference to the list of PASes, for each PAS there is a list of fields (read from input file)
  $new_exons{ $ex } = $exons{ $ex };
  my $chrom = $exons{ $ex }->[0]->[0];
  $pos{ $chrom } = {} if( not defined $pos{ $chrom } );
  foreach my $p ( @{ $pos_per_exons{ $ex } } ) {
    $pos{ $chrom }->{$p} = $ex;
  }
}
%exons = %new_exons;
%new_exons = ();


# now go through the annotation file and check if any annotated gene overlaps
# with a region saved above
print STDERR "[INFO] poly(A) clusters read and filtered. Intersect with gtf now\n";

my @curr_tr = ();

my $fh;
if( $ARGV[0] =~ /gz$/ ){
  open( $fh, "gzip -dc $ARGV[0] | " );
} else {
  open( $fh, "< $ARGV[0]");
}

ENTRY:
while (<$fh>) {
  next if ( $_ =~ /^#/ );
  chomp;
  #print STDERR "." if ( $. % 100000 == 0 );
  my @F   = split("\t");
  my $key = "$F[0]";

  # only consider transcripts and their corresponding exons
  next if ( $F[2] ne "transcript" and $F[2] ne "exon" );

  if( $F[2] eq "transcript" ) {
    my $geneID;
    if ( $F[8] =~ /gene_id\s"([^"]+)"/ ) {
	$geneID = $1;
	$geneID =~ s/\.\d+$//;
    } else {
      print STDERR "[ERROR] Could not infer gene ID for $key:$F[2]:$F[3]:$F[4]\n";
      exit;
    }
    $geneID =~ s/\.\d+$//;

    my $transcriptID;
    if ( $F[8] =~ /transcript_id\s"([^"]+)"/ ) {
	$transcriptID = $1;
	$transcriptID =~ s/\.\d+$//;
    } else {
      print STDERR "[ERROR] Could not infer transcript ID for $key:$F[2]:$F[3]:$F[4]\n";
      exit;
    }
    my $curr_strand = $F[6];

    if( $#curr_tr == -1 ) {
      # initialize first transcript
      @curr_tr = ($transcriptID, $geneID, , $key, $F[3], $F[4], $curr_strand, -1,-1);
    } else {
      # process previous transcript
      my $tr_start = $curr_tr[3];
      my $tr_end = $curr_tr[4];
      my $chrom = $curr_tr[2];
      foreach my $i ( $tr_start .. $tr_end ) {
	if ( defined $pos{$chrom}->{$i} ) {
	  my $curr_exon = $pos{$chrom}->{$i};
	  if( $E2G{ $pos{$chrom}->{$i} } ne $curr_tr[1] ) {

	    # poly(A) site not unambiguously mappable to a gene
	    # delete exon with the poly(A) site (together with all other poly(A) sites
	    # of this exon

	    if ( defined $exons{$curr_exon} ) {
	      delete $exons{$curr_exon};
	      delete $exon_order{$curr_exon};
	      print STDERR "[DEL] Deleted $curr_exon due to overlap with gene $curr_tr[1]\n";
	      $deleted_cnt++;
	    }
	  }
	}
      }
      # update entry for current trancript
      @curr_tr = ($transcriptID, $geneID, $key, $F[3], $F[4], $curr_strand, -1, -1);
    }
  } elsif( $F[2] eq "exon" ) {
    # since gtf is sorted 5' to 3' the exon entry for the transcript 
    # is simply for every exon overwritten. In this way, finally the terminal exon is stored

    # consistency check: exon must belong to current transcript
    my $transcriptID;
    if ( $F[8] =~ /transcript_id\s"([^"]+)"/ ) {
	$transcriptID = $1;
	$transcriptID =~ s/\.\d+$//;
    } else {
      print STDERR "[ERROR] Could not infer transcript ID for exon $key:$F[2]:$F[3]:$F[4]\n";
      exit;
    }

    if ($transcriptID ne $curr_tr[0] ){
      print STDERR "######\n";
      print STDERR "[ERROR] Current exon does not belong to current transcript\n";
      print STDERR "[ERROR] Current exon's transcript: $transcriptID\n";
      print STDERR "[ERROR] Current transcript: $curr_tr[0]\n";
    }

    $curr_tr[6] = $F[3];
    $curr_tr[7] = $F[4];
  }
}
close($fh);

# process last transcript
if( $#curr_tr > -1 ) {
  my $tr_start = $curr_tr[3];
  my $tr_end = $curr_tr[4];
  my $chrom = $curr_tr[2];
  foreach my $i ( $tr_start .. $tr_end ) {
    if ( defined $pos{$chrom}->{$i} ) {
      my $curr_exon = $pos{$chrom}->{$i};
      if( $E2G{ $pos{$chrom}->{$i} } ne $curr_tr[1] ) {

	# poly(A) site not unambiguously mappable to a gene
	# delete exon with the poly(A) site (together with all other poly(A) sites
	# of this exon

	#test
	#if ( $E2G{ $pos{$chrom}->{$i} } eq "ENSG00000031698" ) {
	#  print STDERR "overlapper fuer SARS: $geneID\n";
	#  exit;
	#}

	if ( defined $exons{$curr_exon} ) {
	  delete $exons{$curr_exon};
	  delete $exon_order{$curr_exon};
	  print STDERR "[DEL] Deleted $curr_exon due to overlap with gene $curr_tr[1]\n";
	  $deleted_cnt++;
	}

      } else {
	# current transcript and the terminal exon stored for this pos belong to the same gene
	# however:
	# also check if the current transcript's terminal exon
	# starts downstream of the current exons end (in strand orientation)
	# in this case: don't use the terminal exon stored because it is internally some other transcript
	# of the same gene

	(undef, undef, undef, my $ex_start, my $ex_end) = split(":", $curr_exon);
	my $term_ex_start = $curr_tr[6];
	my $term_ex_end = $curr_tr[7];
	my $strand = $curr_tr[5];
	if( ($strand eq "+" and $term_ex_start > $ex_end) or
	    ($strand eq "-" and $term_ex_end < $ex_start) ) {
	  # current transcript has a terminal exon further downstream than the one with tandem PAS
	  if ( defined $exons{$curr_exon} ) {
	    delete $exons{$curr_exon};
	    delete $exon_order{$curr_exon};
	    print STDERR "[DEL] Deleted $curr_exon due to the terminal exon ($chrom:$term_ex_start-$term_ex_end) of a transcript of $curr_tr[0] ($chrom:$tr_start-$tr_end) that proceeds further downstream than the terminal exon\n";
	    $deleted_cnt++;
	  }
	}
      }
    }
  }
}

print STDERR "\n";

# print all exons with unambiguous poly(A) sites
foreach my $e ( sort { $exon_order{$a} <=> $exon_order{$b} } keys %exons ) {
  foreach my $site ( @{ $exons{$e} } ) {
    print join( "\t", @{$site} ), "\n";
  }
}

# print stats
print STDERR "[INFO] Total number of deleted exons: $deleted_cnt\n";
