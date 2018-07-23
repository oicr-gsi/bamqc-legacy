package GSI::bamqc;
use strict;
use warnings;

BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

    $VERSION = 1.00;
    @ISA     = qw(Exporter);

    #@EXPORT	=	qw();
    @EXPORT_OK = qw(read_bed assess_start_point assess_flag cigar_stats md_stats
      onTarget addRunningBaseCoverage runningBaseCoverage HistStats insertMapping
      load_json toPhred generate_jsonHash
      generate_mismatch_rate generate_indel_rate generate_softclip_rate
      generate_hardclip_rate generate_error_rate get_barcode get_group
      get_raw_reads get_raw_yield get_map_percent get_ontarget_percent
      get_est_yield get_est_coverage);

    #%EXPORT_TAGS = ();
}

sub new {
    my ( $class, %parameters ) = @_;
    my $self = bless( {}, ref($class) || $class );
    return $self;
}

#################### main pod documentation begin ###################

=begin html

<label for="show-menu" class="show-menu">Show Menu</label>
<input type="checkbox" id="show-menu" role="button">

<h1>Genome Sequence Informatics</h1>

=end html


=head1 BamQC

=head2 NAME

bamqc - Generate quality control statistics from BAM files.

=head2 SYNOPSIS

  use GSI::bamqc;

=head2 DESCRIPTION

This library's whole function
is to get enough information to feed to L</"generate_jsonHash($stats,$p)"> and
produce a JSON file with lots and lots of information about the BAM file.

=head2 AUTHOR

L<Genome Sequence Informatics|https://gsi.oicr.on.ca>,
L<Ontario Institute for Cancer Research|https://oicr.on.ca>.
On Github at L<https://github.com/oicr-gsi/bamqc>.

=head2 COPYRIGHT

Copyright (C) 2018 The Ontario Institute for Cancer Research

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at
    your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.


=head1 Subroutines

=cut

#################### main pod documentation end ###################

use Exporter;
use JSON::PP;
use File::Basename;
### for debug
#use Data::Dumper;
#(open my $TTY,"/dev/tty") || die "unable to open keyboard input";

my @month = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );

=for html <hr>

=head2 assess_flag($flag,$stats,$qual,$qcut)

Assesses the sam flag and stores information in the $stats reference hash

B<Arguments>

=over

=item $flag : the sam record flag;

=item $stats : a reference to a hash, which will be modified;

=item $qual : the alignment quality;

=item $param : the parameter hash;

=back

B<Returns>

1 if the read is mapped, 0 if not (modifies $stats with statuses).
The following labels will be incremented in $stats if found.

=over

=item "non primary reads" = 256

=item "unmapped reads = 4;

=item "qual fail reads" = $qual < $qcut ;

=item "mapped reads" = if none of the above are true

=item "paired reads" = 1

=item "mate unmaped reads" [sic] = 8;

=item "properly paired reads" = 2

=back

=cut

sub assess_flag {
    my ( $flag, $stats, $qual, $qcut ) = @_;
    my $mapped = 0;
    ### processing ALL reads, not just sampled
    $$stats{"total reads"}++;

    if ( $flag & 256 ) {
        $$stats{"non primary reads"}++;
    }
    elsif ( $flag & 4 ) {
        $$stats{"unmapped reads"}++;
    }
    elsif ( $qual < $qcut ) {
        $$stats{"qual fail reads"}++;
    }
    else {
        $mapped = 1;
        $$stats{"mapped reads"}++;
        if ( $flag & 1 ) {
            $$stats{"paired reads"}++;
            if ( $flag & 8 ) {
                $$stats{"mate unmaped reads"}++;
            }
        }
        if ( $flag & 2 ) {
            $$stats{"properly paired reads"}++;
        }
    }
    return $mapped;
}

=for html <hr>

=head2 assess_start_point($chrom,$s1,$s2,$sphash)

Keeps a running calculation of ReadsPerStartPoint, updating
value with each read. It is dependent on the stream alignment data
being properly sorted.

Reads per start point is recalculated when the leftstart changes from the previous record


B<Arguments>

=over

=item $chrom = the chromosome of the current record;

=item $s1	 = the position of the current record = the leftmost start point;

=item $s2    = the position of the paired record = the rightmost start point;

=item $sphash	 = the current StartPoint hash, holding the current position and the running ReadsPerStartPoint value;

=back


B<Returns>

A revised StartPoint $sphash hash, with updated stats:

=over

=item "count" : count each distinct pairedstart as a distinct start point

=item "current" : record the current left start point

=item "RPSP" : reads per start point


=cut

sub assess_start_point {
    my ( $chrom, $s1, $s2, $sphash ) =
      @_;   ### receives the left/right start points and a reference to the hash

    my $leftstart = "$chrom\t$s1";
    my $pairStart = "$chrom\t$s1\t$s2";

    #### Recalculate ReadsPerStartPoint when the leftstart changes from the previous record
    if ( $leftstart ne $$sphash{current} ) {
        ### process all stored pairStarts
        for my $sp ( keys %{ $$sphash{pairStarts} } ) {
            $$sphash{count}++
              ;   #### count each distinct pairedstart as a distinct start point

            ### this is a running total adjusting positively or negativeley depending on the number of time the pair is seen
            $$sphash{RPSP} +=
              ( $$sphash{pairStarts}{$sp} - $$sphash{RPSP} ) / $$sphash{count};

            ### delete this pairStart after processing
            delete $$sphash{pairStarts}{$sp};
        }
        ## reset the current position
        $$sphash{current} = $leftstart;
    }

    ### store the pairStart
    $$sphash{pairStarts}{$pairStart}++;    #### collect the paired starts
    return $sphash;                        ### return the modified hash
}


=for html <hr>

=back

=head2 read_bed($file)

Reads in a bed file, storing in C<intervals> in an array of hashes
indicating, Start/Stop and size


B<Arguments>

=over

=item $file = the name of the bed file, containing the intervals

=back

B<Errors>

=over

=item if the bed file can't be opened

=item if the bed file contains an interval of size 0

=back

B<Returns>

Hash structure containing the bed intervals

 %bed {
  @intervals { Start, Stop, Size }
  targetSize
  numberOfTargets
 }

=cut



sub read_bed {
    my ($file) = @_;
    my %bed;
    my $targetCount = 0;
    open( my $BEDFILE, "<", $file )
      or return ("ERROR : Couldn't open target file: $file.\n");
    while (<$BEDFILE>) {
        chomp;
        next if (/^#/);
        my @f = split /\t/;

        $targetCount++;

        my $interval_size = $f[2] - $f[1];

        ### interval_size needs to be larger than 0
        if ( $interval_size < 1 ) {
            return (
"ERROR : the bedfile $file contains an interval with an invalid size $_ \n"
            );
        }

        ### also want to check that the sort order is correct

        push(
            @{ $bed{intervals}{ $f[0] } },
            { Start => $f[1], Stop => $f[2], Size => $interval_size }
        );
        $bed{targetSize} += $interval_size;

        #$bed{Hist}{"$f[0]\t$f[1]\t$f[2]"} = 0;

    }
    close $BEDFILE;

    $bed{numberOfTargets} = $targetCount;
    warn "Loaded " . $bed{numberOfTargets} . " targets.\n\n";

    return \%bed;
}

=for html <hr>

=head2 cigar_stats($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)

Processes the cigarString by breaking into pieces and assessing how
the read is mapped to the reference


B<Arguments>

=over

=item $chrom = the chromosome to which the current read maps;

=item $start1 = the leftmost startpoint;

=item	$start2 = the rightmost startpoint for the paired end;

=item	$R = the read (R1,R2,R?);

=item $strand = the strand to which the read maps, from the sam flag;

=item $cigar = the cigar string;

=item $stats = a reference to the stats hash, which will be modified;

=item	$p = a reference to the parameter hash;

=back

B<Returns>

A list with two values, $readLength and $mappedBases. Increments $stats with values:

=over

=item C<alignedCount>, C<hardClipCount>, C<softClipCount>, C<insertCount>, C<deletionCount>

=item C<byCycle> > C<aligned> > read > cycle

=item C<byCycle> > C<hardClip> > read > cycle

=item C<byCycle> > C<softClip> > read > cycle

=item C<byCycle> > C<insertion> > read > cycle

=item C<byCycle> > C<deletion> > read > cycle

=back

B<Errors>

=over

=item if the cigar string contains an unknown operation

=back

=cut

sub cigar_stats {
    my ( $chrom, $start1, $start2, $R, $strand, $cigar, $stats, $p ) = @_;

    ## the fragment positioning on the chromosome
    my $pairStart = "$chrom\t$start1\t$start2";

    ## split the cigar string, and process each piece
    my @Op = procCigar( $cigar, $strand );

    ## initialize these values
    my $readLength  = 0;
    my $cycle       = 1;    ### this will keep track of the current cycle
    my $posOffset   = 0;
    my $mappedBases = 0;

    foreach my $c (@Op) {
        ### MATCH, capture the number of bases in $1
        if ( $c =~ /(.*)M/ ) {
            $$stats{alignedCount} += $1;

            #for (my $i = 0; $i < $1; $i++){
            for my $i ( 1 .. $1 ) {
                $$stats{byCycle}{aligned}{$R}{$cycle}++;
                $cycle++;
            }
            $readLength  += $1;
            $mappedBases += $1;

            ### HARDCLIPPED BASES
        }
        elsif ( $c =~ /(.*)H/ ) {
            $$stats{hardClipCount} += $1;
            for my $i ( 1 .. $1 ) {
                $$stats{byCycle}{hardClip}{$R}{$cycle}++;
                $cycle++;
            }
            $readLength += $1;
            ### SOFTCLIPPED BASES
        }
        elsif ( $c =~ /(.*)S/ ) {
            $$stats{softClipCount} += $1;
            for my $i ( 1 .. $1 ) {
                $$stats{byCycle}{softClip}{$R}{$cycle}++;
                $cycle++;
            }
            $readLength += $1;
            ### INSERTED BASES
        }
        elsif ( $c =~ /(.*)I/ ) {
            $$stats{insertCount} += $1;
            for my $i ( 1 .. $1 ) {
                $$stats{byCycle}{insertion}{$R}{$cycle}++;
                $cycle++;
            }
            $readLength += $1;
            ### DELETED BASES
        }
        elsif ( $c =~ /(.*)D/ ) {
            $$stats{deletionCount} += $1;
            for my $i ( 1 .. $1 ) {
                $$stats{byCycle}{deletion}{$R}{$cycle}++;
            }
            $mappedBases += $1;
        }
        else {
            die "Can't handle CIGAR operation: $c\n";
        }
    }
    return ( $readLength, $mappedBases );

}

=for html <hr>

=head2 md_stats($R,$mdstring,$strand,$stats)

Processes the MD:Z tag in the SAM.


B<Arguments>

=over

=item $R = the read; one of R1,R2,R?;

=item $mdstring = the mismatch string, extracted from the MD:Z tag in the sam record;

=item $strand = to which strand does the read map;

=item $stats = reference to the stats hash, values will be modified;

=back

B<Returns>

The number of mismatches defined by the MD string.
Increments the following in the $stats hash, if encountered.

=over

=item C<mismatchCount>

=item C<byCycle> > C<mismatch> > Read > cycle

=back


=cut

sub md_stats {
    my ( $R, $mdstring, $strand, $stats ) = @_;
    my @Op    = procMD( $mdstring, $strand );
    my $cycle = 1;
    my $mm    = 0;

    for my $md (@Op) {
        if ( $md =~ /^([0-9]+)/ ) {

            # ignore matching bases
            $cycle += $1;
        }
        elsif ( $md =~ /^(\^[A-Z]+)/ ) {

            # ignore deletions
        }
        elsif ( $md =~ /^([A-Z]+)/ ) {    # mismatch!
            foreach my $i ( split( //, $1 ) ) {
                $$stats{mismatchCount}++;
                $$stats{byCycle}{mismatch}{$R}{$cycle}++;
                $cycle++;
                $mm++;
            }
        }
        else {
            die "Couldn't handle MD operation: $md\n";
        }
    }

    return $mm;
}

=for html <hr>

=head2 onTarget($chrom,$start,$mapped,$stats)

Determines if the read overlaps to any degree intervals contained
						within the bedfile record.  The bedfile record is modified to
						indicated mapping to this interval. One read may map to multiple intervals.

B<What is an on-target read?>

A read that overlaps the target region by any number of bases (0=false; 1=true).

 ---------------XXXXXXXXXXXX------------
   0000000000
           1111111111
              1111111111111111
                  11111111
              1111111111111
                               0000000000


B<Arguments>

=over

=item $chrom = the chromosome to which the read maps;

=item $start = the leftmost start position;

=item $mapped = the number of bases mapped;

=item $bed = a reference to the bed hash;

=item $stats = a reference to the stats hash containing the bed intervals, which will be modified;

=back

B<Returns>

1 if mapped to Target, 0 if not.
Increments the following in the $stats hash, if encountered.

=over

=item C<bed> > C<intervals> > chromosome > index > C<hist> : adds number of mapped

=back

=cut

sub onTarget {
    my ( $chrom, $start, $mapped, $stats ) = @_;
    my $onTarget = 0;

    ## bed intervals are hash:array:hash
    ## get array of intervals from the current chromosome
    if ( $$stats{bed}{intervals}{$chrom} ) {
        my @intervals = @{ $$stats{bed}{intervals}{$chrom} };
        ### extract the indices of intervals which meet the necessary conditions
        ### map start position is less than the interval stop
        ### map start + mapped bases is greater than the interval start
        my @idxs = grep {
                 ( $start <= $intervals[$_]{Stop} )
              && ( ( $start + $mapped ) >= $intervals[$_]{Start} )
        } ( 0 .. $#intervals );

        my $idxcount =
          scalar @idxs;   ## should only be ONE interval that contains the start

        ### the read may map to multiple intervals
        ### can either build hist on the 1st interval, or on all intervals
        ### I think the original script only did the first
        for my $idx ( sort { $a <=> $b } @idxs ) {
            next
              if ($onTarget)
              ;    ### will turn this off if it should be mapping to all
            ## there is overlap between the read and the interval, but not all bases are necessarily mapped to the interval
            ## the logic here indicates that all mapp
            $$stats{bed}{intervals}{$chrom}[$idx]{hist} += $mapped;
            $onTarget = 1;

        }
    }
    return $onTarget;
}

=for html <hr>

=head2 addRunningBaseCoverage($chrom,$start1,$start2,$cigar,$strand,$stats)

Adds fragments to the runningBaseCoverage collection.  Each fragment
is added to all positions to which the read maps.
This collection will be continuously process and cleared of all
positions that precede the current read start

B<Arguments>

=over

=item $chrom = the chromosome to which the current read maps;

=item $start1 = the left-most start of the paired end reads;

=item $start2 = the right-most start of the paired end reads;

=item $cigar = the cigar string;

=item $strand = the strand to which the read maps;

=item $stats = a reference to the stats hash;

=back

B<Returns>

The number of positions to which the current fragment is stored/mapped.
Increments the following in the $stats hash, if encountered.

=over

=item C<runningBaseCoverage> > chromosome > start + mismatch/deletion offset > "$chrom\t$start1\t$start2"

=back

=cut

sub addRunningBaseCoverage {
    my ( $chrom, $start1, $start2, $cigar, $strand, $stats ) = @_;
    my @Op = procCigar( $cigar, $strand );

    my $posOffset = 0;
    my $fragment  = "$chrom\t$start1\t$start2";

    ### first review the cigar string, for matches or deletions
    foreach my $c (@Op) {
        if ( $c =~ /(.*)[M|D]/ ) {    ## MISMATCH or DELETIONS
            for my $i ( 1 .. $1 ) {
                $$stats{runningBaseCoverage}{$chrom}{ $start1 + $posOffset }
                  {$fragment}++;
                $posOffset++;
            }
        }
    }
    return $posOffset;
}

=for html <hr>

=head2 runningBaseCoverage($stats,$chrom,$startpos)

Calculates base coverage as a running total. This requires that the
						streaming records are from a sorted bam file.
						This will produce a hash with chromosomal keys, each value being a
						hash of fragment positions and counts, by parsing the cigar string
						for the current record.
						When a new mapping position is found, stats are generated on all
						stored positions, and then the hash is cleared

B<Arguments>

=over

=item $stats = a reference to the stats hash, holding the runningBaseCoverage hash.  Other keys in the stats hash will be modified;

=item $chrom = the current chromosome, all recorded positions not on this chromosome will be processed and cleared.  If no value, then the whole hash is cleared;

=item $pos = the current position, all recorded positions on the current chromosome, before this position, will be cleared;

=back

B<Returns>

The number of cleared positions from the runningBaseCoverage hash.
Increments the following in the $stats hash, if encountered.

=over

=item C<collapsedCoverageHist> > index

=item C<nonCollapsedCoverageHist> > index

=back

=cut

sub runningBaseCoverage {
    my ( $stats, $chrom, $startpos ) = @_;

    my $cleared_positions = 0;
    ### will clear out all startpoints stored that precede the current chromosome and position
    for my $chr ( keys %{ $$stats{runningBaseCoverage} } ) {
        for my $pos ( sort { $a <=> $b }
            keys %{ $$stats{runningBaseCoverage}{$chr} } )
        {    ### check each start position

            ##process if not the current chromosome or current chromosome + before current position, or if the $chrom OR $pos variables are not supplied
            if (   !$chr
                || !$startpos
                || ( $chr ne $chrom )
                || ( $pos < $startpos ) )
            {
                my %startPoints = %{ $$stats{runningBaseCoverage}{$chr}{$pos} };

                ### how many start points/unique fragments are stored on this chromosome positionn
                my $numStartPoints = scalar keys %startPoints;

                ### collapsedCoverageHist : histogram of number of Unique fragments at a given position
                for my $i ( 1 .. $numStartPoints ) {
                    $$stats{collapsedCoverageHist}{$i}++;
                }

                ### nonCollapsedCoverageHist : histogram of number of fragments at a given position
                my $totalDepth = 0;
                map { $totalDepth += $_ } values %startPoints;

                for my $i ( 1 .. $totalDepth ) {
                    $$stats{nonCollapsedCoverageHist}{$i}++;
                }

                delete $$stats{runningBaseCoverage}{$chr}{$pos}
                  ;   ### delete this position from the running BaseCoverageHash
                $cleared_positions++;
            }

        }
        delete $$stats{runningBaseCoverage}{$chr}
          if ( !$chr || ( ( defined $chrom ) and ( $chr ne $chrom ) ) )
          ;           ### delete previous chromosome keys
    }

    return $cleared_positions;
}

=for html <hr>

=head2 HistStats(%val)

Calculate the mean and standard deviation of insert sizes in a hash

B<Arguments>

=over

=item %val = a hash, with keys = insert size, values = count for each insert size

=back

B<Returns>

Mean and standard deviation of the insert sizes.

=cut

sub HistStats {

    #my %val = %{ $_[0] };

    my ($val) = @_;
    my ( $mean, $stdv ) = ( 0, 0 );

    if ($val) {
        my ( $sum, $count ) = ( 0, 0 );
        map {
            $sum   += ( $_ * $$val{$_} );
            $count += $$val{$_};
        } keys %$val;
        $mean = $count > 0 ? $sum / $count : 0;

        my ($squareDiff);
        map {
            $squareDiff += ( ( ( $_ - $mean ) * ( $_ - $mean ) ) * $$val{$_} );
        } keys %$val;

        $stdv = $count > 0 ? sqrt( $squareDiff / $count ) : 0;
    }

    return ( $mean, $stdv );

}

=for html <hr>

=head2 insertMapping($tlen,$rnext,$hash,$p)

Identifies the paired-end insert size as being within the normal
						range, abnormally far, or on different chromosomes

B<Arguments>

=over

=item $tlen = template length;

=item $rnext = the chromosome of the other in the pair, = indicates the same;

=item $hash = reference to a hash collecting statistics;

=item $p = parameter hash;

=back

B<Returns>

A description of the insert mapping.
  Modifies the reference hash at the following key and returns a string, one of:

=over

=item "normalInsertSize" : if the template length is less than  "normalInsertMax" from $p

Also stores the template length in the hash.

=item "pairsMappedAbnormallyFar" : if the template length is longer than "normalInsertMax" from $p

=item "pairsMappedToDifferentChr" : if the pair is mapped to a different chromosome

=back

=cut

sub insertMapping {
    my ( $tlen, $rnext, $hash, $p ) = @_;
    my $class = "";
    if ( ( $tlen > 0 ) && ( $rnext eq "=" ) ) {  ### maps to the same chromosome
        if ( $tlen < $$p{normalInsertMax} ) {
            ### within expected size , will store the template length for histogram construction
            $$hash{normalInsertSizes}{$tlen}++;
            $class = "normalInsertSize";
        }
        else {
            $$hash{pairsMappedAbnormallyFar}++;
            $class = "pairsMappedAbnormallyFar";
        }
    }
    elsif ( $rnext ne "*" ) {    ### mapped to a different chromosome
        $$hash{pairsMappedToDifferentChr}++;
        $class = "pairsMappedToDifferentChr";
    }
    return $class;
}

=for html <hr>

=head2 load_json(@files)

Open, decode, and store each JSON file in a hash with filename keys

B<Arguments>

=over

=item @files = a list of json file paths

=back

B<Returns>

A hash with basename(filename) keys > decoded JSON hash

=cut

sub load_json {
    my @files = @_;
    my %json_hash;
    for my $file (@files) {
        print STDERR "reading from $file\n";
        open( my $FILE, "<", $file ) or die "Couldn't open $file.\n";
        if ( my $line = <$FILE> ) {
            $json_hash{ basename($file) }           = decode_json($line);
            $json_hash{ basename($file) }{basename} = basename($file);
            $json_hash{ basename($file) }{dirname}  = dirname($file);
        }
        else {
            warn "No data found in $file!\n";
        }
    }
    return %json_hash;
}

=for html <hr>

=head2 toPhred($char)

Convert a character to a phred score (ascii value - 33)


B<Arguments>

=over

=item $char = the single character to convert

=back

B<Returns>

Phred score

=cut

sub toPhred {
    my $char   = $_[0];
    my $ascii  = ord($char);
    my $offset = 33;
    return $ascii - $offset;
}


=for html <hr>

=head2 generate_jsonHash($stats,$p)

Transforms information in the stats hash to the appropriate keys and
values for the jsonHash.

=over

=item "aligned bases"

The (estimated) number of bases that are aligned. alignedCount + insertCount
from L</"cigar_stats($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)">
multiplied by the sampling rate (bamqc.pl).

=item "average read length"

Float. the average read length is the sum of observed bases
(aligned + soft clipped + hard clipped + inserted) / # of mapped reads,
multiplied by the sampling rate (bamqc.pl).

=item "collapsed bases covered"

Integer. See L</"runningBaseCoverage($stats,$chrom,$startpos)">

=item "deleted bases"

The (estimated) number of bases that are hard-clipped. deletionCount
from L</"cigar_stats($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)">
multiplied by the sampling rate (bamqc.pl).

=item "hard clip bases"

The (estimated) number of bases that are hard-clipped. hardClipCount
from L</"cigar_stats($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)">
multiplied by the sampling rate (bamqc.pl).

=item "insert histogram"

See L</"insertMapping($tlen,$rnext,$hash,$p)">. Or empty set {}

=item "insert mean"

Float. See L</"HistStats(%val)">. Float if exists or String "0" (bamqc.pl).

=item "insert stdev"

Float. See L</"HistStats(%val)">. Float if exists or String "0" (bamqc.pl).

=item "inserted bases"

The (estimated) number of bases that are inserted. insertionCount
from L</"cigar_stats($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)">
multiplied by the sampling rate (bamqc.pl).

=item "mapped reads"

Integer. A read is considered mapped if it does not have flags: non-primary
(256) or unmapped (4), or fell below the quality cutoff.
See L</"assess_flag($flag,$stats,$qual,$qcut)">

=item "mate unmaped reads" [sic]

Integer. A mate is considered unmapped if it has the flags: paired (1) and mate unmapped (8).
See L</"assess_flag($flag,$stats,$qual,$qcut)">

=item "mismatch bases"

=item "non collapsed bases covered"

nonCollapsedCoverageHist

=item "non primary reads"

Integer. See L</"assess_flag($flag,$stats,$qual,$qcut)">.

=item "number of ends"

String. "paired end" if there are more than 0 properly paired reads, "single end" otherwise.

=item "number of targets"

Integer. Set in L</"read_bed($file)">. Corresponds to number of lines in the BED
file.

=item "paired reads"

Integer. A read is considered paired if it is a mapped read and has sam flag:
paired (1). See L</"assess_flag($flag,$stats,$qual,$qcut)">.

=item "properly paired reads"

Integer. A read is considered properly paired if it is a mapped read and has sam
flag: properly paired (2). See L</"assess_flag($flag,$stats,$qual,$qcut)">

=item "qual cut"

Integer. Quality cutoff for reads. Supplied to bamqc.pl by -q or default 30.

=item "qual fail reads"

Integer. A read is considered qual fail if it falls below the quality cutoff,
supplied by -q. Default 30. See L</"assess_flag($flag,$stats,$qual,$qcut)">

=item "$read aligned by cycle"

List. See L</"cigar_stats($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)">.

=item "$read average length"

=item "$read deletion by cycle"

List. See L</"cigar_stats($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)">.

=item "$read hard clip by cycle"

List. See L</"cigar_stats($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)">.

=item "$read insertion by cycle"

List. See L</"cigar_stats($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)">.

=item "$read length histogram"

Hash. Contains buckets for every read length (bamqc.pl).

=item "$read mismatch by cycle"

Hash. Searches for "MD:Z:*" string and passes it to
L</"md_stats($R,$mdstring,$strand,$stats)">.

=item "$read quality by cycle"

Hash. (bamqc.pl).

=item "$read quality histogram"

Hash. (bamqc.pl).

=item "$read soft clip by cycle"

List. See L</"cigar_stats($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)">.

=item "reads on target"

Integer. Incremented if L</"onTarget($chrom,$start,$mapped,$stats)"> in bamqc.pl.

=item "reads per start point"

Float. See L</"assess_start_point($chrom,$s1,$s2,$sphash)">.

=item "soft clip bases"

The (estimated) number of bases that are soft-clipped. softClipCount
from L</"cigar_stats($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)">
multiplied by the sampling rate.

=item "target file"

Path to the BED target file.

=item "target size"

The total sum of the target sizes from the target file. See L</"read_bed($file)">.
Targets may be overlapping and so the target space could be higher than the
actual footprint.

=item "total reads"

Integer. Total number of reads * 1.

=item "unmapped reads"

Integer. See L</"assess_flag($flag,$stats,$qual,$qcut)">

=back

B<Arguments>

=over

=item $stats = reference to the stats hash, where data is stored;

=item $p = reference to the parameters hash;

=back

B<Returns>

The jsonHash, ready to print.


=cut

sub generate_jsonHash {
    my ( $stats, $p ) = @_;

    my %jsonHash = map { ( $_, $$stats{$_} ) } (
        "number of ends",
        "total reads",
        "mapped reads",
        "unmapped reads",
        "non primary reads",
        "paired reads",
        "properly paired reads",
        "mate unmaped reads",
        "qual fail reads",
        "qual cut",
        "aligned bases",
        "mismatch bases",
        "inserted bases",
        "deleted bases",
        "soft clip bases",
        "hard clip bases",
        "reads per start point",
        "reads on target",
        "target file",
        "target size",
        "number of targets",
    );

    $jsonHash{"total reads"} = $$stats{"total reads"} * 1;

    $jsonHash{"average read length"} = $$stats{averageReadLength}{overall};
    $jsonHash{"insert mean"}         = $$stats{meanInsert}
      || "0";    ## this is quoted to be consistent with samStats.pl
    $jsonHash{"insert stdev"} = $$stats{stdevInsert}
      || "0";    ## this is quoted to be consistent with samStats.pl
    ### capture any jsonHash elements storeed in the parameters

    map { $jsonHash{$_} = $$p{jsonHash}{$_}; } keys %{ $$p{jsonHash} };

    if ( $$p{reportBasesCovered} ) {
        $jsonHash{"non collapsed bases covered"} =
          $$stats{nonCollapsedCoverageHist};
        $jsonHash{"collapsed bases covered"} = $$stats{collapsedCoverageHist};
    }
    $jsonHash{"insert histogram"} = $$stats{normalInsertSizes} || {};

    for my $R (qw/R1 R2 R?/) {
        ( my $read = $R ) =~ s/R/read /;
        $jsonHash{"$read quality histogram"} = $$stats{qualHist}{$R} || {};
        $jsonHash{"$read quality by cycle"} = $$stats{qualLine}{$R}
          ; ### this will get an undef value, not empty hash, to match samStats.pl
        $jsonHash{"$read length histogram"} = $$stats{readLengthHist}{$R} || {};
        $jsonHash{"$read average length"} = $$stats{averageReadLength}{$R} || 0;
        $jsonHash{"$read aligned by cycle"} =
          $$stats{byCycle}{aligned}{$R} || {};
        $jsonHash{"$read mismatch by cycle"} =
          $$stats{byCycle}{mismatch}{$R} || {};
        $jsonHash{"$read insertion by cycle"} =
          $$stats{byCycle}{insertion}{$R} || {};
        $jsonHash{"$read deletion by cycle"} =
          $$stats{byCycle}{deletion}{$R} || {};
        $jsonHash{"$read soft clip by cycle"} =
          $$stats{byCycle}{softClip}{$R} || {};
        $jsonHash{"$read hard clip by cycle"} =
          $$stats{byCycle}{hardClip}{$R} || {};
    }

    $jsonHash{"number of ends"} =
      $$stats{"properly paired reads"} > 0 ? "paired end" : "single end";

    $jsonHash{"number of targets"}     = $$p{bed}{numberOfTargets};
    $jsonHash{"qual cut"}              = $$p{qualCut};
    $jsonHash{"reads per start point"} = $$stats{startPoint}{RPSP};
    $jsonHash{"target size"}           = $$p{bed}{targetSize};
    $jsonHash{"target file"}           = $$p{bedFile};

    return %jsonHash;

}

=for html <hr>

=head2 generate_error_rate( $hash, $prefix )

Compute the error rate with (mismatch+insertion+deletion)/aligned

B<Arguments>

=over

=item $hash = The hash containing the JSON file contents to analyse;

=item $prefix = A prefix to use when accessing the keys of hash;

=back

B<Returns>

A formatted percentage

=cut

sub generate_error_rate {
    my ( $hash, $prefix ) = @_;
    return generateRatePercent( $hash, $prefix,
        [ "mismatch", "insertion", "deletion" ],
        ["aligned"] );
}

=for html <hr>

=head2 generate_mismatch_rate( $hash, $prefix )

Compute the mismatch rate with mismatch/aligned

B<Arguments>

=over

=item $hash = The hash containing the JSON file contents to analyse;

=item $prefix = A prefix to use when accessing the keys of hash;

=back

B<Returns>

A formatted percentage

=cut

sub generate_mismatch_rate {
    my ( $hash, $prefix ) = @_;
    return generateRatePercent( $hash, $prefix, ["mismatch"], ["aligned"] );
}

=for html <hr>

=head2 generate_indel_rate( $hash, $prefix )

Compute the indel rate with (insertion+deletion)/aligned

B<Arguments>

=over

=item $hash = The hash containing the JSON file contents to analyse;

=item  $prefix = A prefix to use when accessing the keys of hash;

=back

B<Returns>

A formatted percentage

=cut

sub generate_indel_rate {
    my ( $hash, $prefix ) = @_;
    return generateRatePercent( $hash, $prefix, [ "insertion", "deletion" ],
        ["aligned"] );
}

=for html <hr>

=head2 generate_softclip_rate( $hash, $prefix )

Compute the softclip rate with soft clip/(soft clip + aligned)

B<Arguments>

=over

=item $hash = The hash containing the JSON file contents to analyse;

=item $prefix = A prefix to use when accessing the keys of hash;

=back

B<Returns>

A formatted percentage

=cut

sub generate_softclip_rate {
    my ( $hash, $prefix ) = @_;
    return generateRatePercent( $hash, $prefix, ["soft clip"],
        [ "soft clip", "aligned" ] );
}

=for html <hr>

=head2 generate_hardclip_rate( $hash, $prefix )

Compute the hardclip rate with hard clip/(hard clip + soft clip + aligned)

B<Arguments>

=over

=item $hash = The hash containing the JSON file contents to analyse;

=item $prefix = A prefix to use when accessing the keys of hash;

=back

B<Returns>

A formatted percentage

=cut

sub generate_hardclip_rate {
    my ( $hash, $prefix ) = @_;
    return generateRatePercent( $hash, $prefix, ["hard clip"],
        [ "soft clip", "aligned", "hard clip" ] );
}


=for html <hr>

=head2 get_barcode( $jsonHash)

Get the sequencing index / barcode

B<Arguments>

=over

=item $jsonHash = The hash containing the JSON file contents to analyse.

=back

B<Returns>

The barcode if it exists; otherwise 'NoIndex'

=cut

sub get_barcode {
    my ($jsonHash) = @_;
    return exists $jsonHash->{barcode} ? $jsonHash->{barcode} : 'NoIndex';
}

=for html <hr>

=head2 get_group( $jsonHash)

Get the group id and group id description

B<Arguments>

=over

=item $jsonHash = The hash containing the JSON file contents to analyse.

=back

B<Returns>

The group id and group id description separated by a space, if they exist.
If they don't exist, then 'na'.


=cut

sub get_group {
    my ($jsonHash) = @_;
    return
        exists $jsonHash->{"group id"}
      ? exists $jsonHash->{"group id description"}
          ? $jsonHash->{"group id description"} . " " . $jsonHash->{'group id'}
          : $jsonHash->{'group id'}
      : "na";
}

=for html <hr>

=head2 get_raw_reads( $jsonHash)

Get the total raw reads, counted by summing mapped, unmapped and qual fail reads

B<Arguments>

=over

=item $jsonHash = The hash containing the JSON file contents to analyse.

=back

B<Returns>

The total number of reads (int).

=cut

sub get_raw_reads {
    my ($jsonHash) = @_;
    return
      int( $jsonHash->{"mapped reads"} +
          $jsonHash->{"unmapped reads"} +
          $jsonHash->{"qual fail reads"} );
}

=for html <hr>

=head2 get_raw_yield( $jsonHash)

Get the total raw yield, multipying get_raw_reads by the average read length

B<Arguments>

=over

=item $jsonHash = The hash containing the JSON file contents to analyse.

=back

B<Returns>

The total yield (int).

=cut

sub get_raw_yield {
    my ($jsonHash) = @_;
    return int( get_raw_reads($jsonHash) * $jsonHash->{"average read length"} );
}


=for html <hr>

=head2 get_map_percent( $jsonHash)

Get the total map percentage, calculated by dividing mapped reads by
						total number of reads and multiplying by 100. If total reads is 0,
						treat as 1.

B<Arguments>

=over

=item $jsonHash = The hash containing the JSON file contents to analyse.

=back

B<Returns>

The map percentage as an integer.

=cut

sub get_map_percent {
    my ($jsonHash) = @_;
    return 100 *
      ( $jsonHash->{"mapped reads"} / ( get_raw_reads($jsonHash) || 1 ) );
}


=for html <hr>

=head2 get_ontarget_percent( $jsonHash)

Get the total on target percentage by dividing reads on target by
						the number of mapped reads, and multiplying by 100. If mapped reads
						is 0, treat as 1.

B<Arguments>

=over

=item $jsonHash = The hash containing the JSON file contents to analyse.

=back

B<Returns>

The on target percentage as an integer

=cut

sub get_ontarget_percent {
    my ($jsonHash) = @_;
    return 100 *
      ( $jsonHash->{"reads on target"} / ( $jsonHash->{"mapped reads"} || 1 ) )
      ;    # $rawReads) * 100;   # could argue using this either way
}

=for html <hr>

=head2 get_est_yield( $jsonHash)

Get the estimated total yield by multiplying total aligned based by
						the on target percentage, and dividing that by reads per start point.
						If reads per start point is is 0, treat as 1.

B<Arguments>

=over

=item $collapse = whether to use reads per start point to collapse down the coverage;

=item $jsonHash = The hash containing the JSON file contents to analyse.

=back

B<Returns>

The estimated yield as an integer

=cut

sub get_est_yield {
    my ( $jsonHash, $collapse ) = @_;
    my $denom = 1;
    if ($collapse) {
        $denom = ( $jsonHash->{"reads per start point"} || 1 );
    }
    return
      int( $jsonHash->{"aligned bases"} *
          ( get_ontarget_percent($jsonHash) / 100 ) /
          $denom );
}

=for html <hr>

=head2 get_est_coverage( $jsonHash)

Get the estimated total coverage by dividing estimated yield by
						the target size.
						If the target is is 0, treat as 1.

B<Arguments>

=over

=item $collapse = whether to use reads per start point to collapse down the coverage;

=item $jsonHash = The hash containing the JSON file contents to analyse.

=back

B<Returns>

The estimated coverage as an integer.

=cut

sub get_est_coverage {
    my ( $jsonHash, $collapse ) = @_;
    return get_est_yield( $jsonHash, $collapse ) /
      ( $jsonHash->{"target size"} || 1 );
}

=for html <hr>

=head2 findStart($cigarOp, $start)

Find the start index of the read using the cigar string

B<Arguments>

=over

=item $cigarOp = the cigar string

=back

B<Returns>

The new start position, adjusted to take soft clipping into account

=cut

sub findStart {
    my @cigarOp = @{ $_[0] };
    my $start   = $_[1];

    if ( $cigarOp[0] =~
        /(.*)S/ )  # if first cigar operation is a soft clip, adjust start point
    {
        $start -= $1;
    }

    return $start;
}

=for html <hr>

=head2 findEnd($cigarOp, $end)

Find the start index of the read using the cigar string

B<Arguments>

=over

=item $cigarOp = the cigar string;

=item $end = the assumed end of the string;

=back

B<Returns>

The new start position, adjusted to take soft clipping into account

=cut

sub findEnd {
    my @cigarOp = @{ $_[0] };
    my $end     = $_[1];
    if ( $cigarOp[0] =~ /(.*)S/ )
    {    # if first cigar operation is a soft clip, adjust start point
        $end -= $1;
    }
    foreach my $cig (@cigarOp) {
        if ( $cig =~ /(.*)S/ ) {
            $end += $1;
        }
        elsif ( $cig =~ /(.*)M/ ) {
            $end += $1;
        }
        elsif ( $cig =~ /(.*)D/ ) {
            $end += $1;
        }
    }
    return $end;
}


##################### INTERNAL Subroutines ##############################
##get_est_yield
## INTERNAL procCigar()
## Arguments:
## 		$cigar : the cigar string
## 		$strand : the strand to which the read maps
## Returns : an array of cigar pieces
## Description : breaks the cigar string into pieces denoted by length/type, accounting for strand
sub procCigar {
    my ( $cigar, $strand ) = @_;
    my @cigarOp;

    while ( $cigar =~ /^([0-9]+[MIDNSHPX=]).*$/ ) {
        push( @cigarOp, $1 );
        $cigar =~ s/$1//;
    }
    @cigarOp = reverse(@cigarOp) if ( $strand eq "-" );

    return @cigarOp;
}

## INTERNAL procMD()
## Arguments>
##	$md : the mdstring
##  $strand : the strand to which the read maps
## Returns> : an array of mdstring pieces
## Descriptions : breaks the mdstring into pieces, accounting for strand
sub procMD {
    my ( $md, $strand ) = @_;
    my @mdOp;
    while ( $md ne "" ) {
        if ( $md =~ /^([0-9]+)/ ) {
            push( @mdOp, $1 );
            $md =~ s/^$1//;
        }
        if ( $md =~ /^([A-Z]+)/ ) {
            push( @mdOp, $1 );
            $md =~ s/^$1//;
        }
        if ( $md =~ /^(\^)([A-Z]+)/ ) {
            push( @mdOp, "^$2" );
            $md =~ s/^\^$2//;
        }
    }

    @mdOp = reverse(@mdOp) if ( $strand eq "-" );

    return @mdOp;
}

# Compute a rate from the properties of a run.
#
# 1[HashRef]: The hash containing the JSON file contents to analyse.
# 2[Str]: A prefix to use when accessing the keys of hash.
# 3[ArrayRef[Str]]: The names of the "by cycle" keys that belong in the numerator.
# 4[ArrayRef[Str]]: The names of the "by cycle" keys that belong in the denominator.
# Returns[Str]: A formatted percentage of all the sum of all numerator items divided by the sum of all the denominator items.
sub generateRatePercent {
    my ( $hash, $prefix, $numerator_names, $denominator_names ) = @_;
    my $numerator = 0;
    for my $name ( @{$numerator_names} ) {
        $numerator +=
          byCycleToCount( $hash->{ $prefix .' '. $name . " by cycle" } );
    }
    my $denominator = 0;
    for my $name ( @{$denominator_names} ) {
        $denominator +=
          byCycleToCount( $hash->{ $prefix .' '. $name . " by cycle" } );
    }
    return $numerator * 100.0 / ( $denominator || 1 );
}



sub byCycleToCount {
    my $histRef = $_[0];

    my $sum = 0;
    for my $i ( keys %{$histRef} )
    {    # TODO(apmasell): How is this not just values?
        $sum += $histRef->{$i};
    }

    return $sum;
}

######################## END INTERNAL Subroutines########################

1;
