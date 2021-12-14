###############################################################################
#
#   Filter tandem polyA sites by length
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

my $upstream   = 150;
my $downstream = 150;

GetOptions(
    "upstream=i"   => \$upstream,
    "downstream=i" => \$downstream
);

if ( not defined $ARGV[0] ) {
    print STDERR
"[ERROR] Usage: perl $0 --upstream=100 --downstream=100 human.tandem.unfiltered.tsv\n\n";
    exit;
}

my %exons  = ();
my @header = ();

open( F, "< $ARGV[0]" ) or die "Can't open $ARGV[0]\n";
while (<F>) {
    chomp;
    my @F = split("\t");
    if ( $. == 1 ) {

        # check the header line
        if (   $F[0] ne "chrom"
            or $F[1] ne "start"
            or $F[2] ne "end"
            or $F[3] ne "pas"
            or $F[4] ne "score"
            or $F[5] ne "strand"
            or $F[6] ne "polyAsite_exon_idx"
            or $F[7] ne "nr_polyAsites_on_exon"
            or $F[8] ne "exon"
            or $F[9] ne "gene" )
        {
            print STDERR
              "[ERROR] Header in file $ARGV[0] has unexpected entries\n";
            my @example = (
                "chrom", "start", "end", "pas", "score", "strand",
                "polyAsite_exon_idx", "nr_polyAsites_on_exon", "exon", "gene"
            );
            print STDERR "[ERROR] expected fields: ", join( "\t", @example ),
              "\n";
            exit(2);
        }
        if ( $#F < 10 ) {
            print STDERR "[ERROR] Sample names in header are missing\n";
            exit(2);
        }
        @header = @F;
        next;
    }
    $exons{ $F[8] } = [] if ( not defined $exons{ $F[8] } );

    # check consistent number of entries in the current line and in the header
    if ( $#F < $#header ) {
        print STDERR "[ERROR] Inconsistent number of entries in line $.\n";
        print STDERR "[ERROR] Expected: ", $#header, ", given: ", $#F, "\n";
        exit(2);
    }

    #save the line of each pA as array-ref in the exons hash
    push @{ $exons{ $F[8] } }, \@F;
}
close(F);

#store all sites that are too close to each other
my %discard = ();

# print header
print join( "\t", @header ), "\tPAS_overlap\n";

#iterate over all exons
foreach my $ex ( keys %exons ) {

    #delete all %discard keys
    %discard = ();

    #store the tandem polyAsites of this exon
    my @pAs = @{ $exons{$ex} };
    foreach my $i ( 0 .. $#pAs - 1 ) {
        ( my $site1 = $pAs[$i]->[3] ) =~ s/.+:(\d+):[\+-]$/$1/;
        ( my $site2 = $pAs[ $i + 1 ]->[3] ) =~ s/.+:(\d+):[\+-]$/$1/;
        if ( abs( $site1 - $site2 ) < $downstream + $upstream ) {
            $discard{$i} = 1;
            $discard{ $i + 1 } = 1;
        }
    }
    foreach my $i ( 0 .. $#pAs ) {
        if ( not defined $discard{$i} ) {
            print join( "\t", @{ $pAs[$i] } ), "\tOK\n";
        }
        else {
            print join( "\t", @{ $pAs[$i] } ), "\tOVERLAP\n";
        }
    }
}
