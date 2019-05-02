#!/usr/bin/env perl
# Copyright (C) 2017 The Ontario Institute for Cancer Research
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
#      your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.

=head1 NAME

bamqc.pl - generation of bam file statistics

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

Processes a stream of data generated by "samtools view", producing a variety of QC statistics in json format

=over 2

=item * Usage: bamqc.pl 

=head1 AUTHOR

Original Script : Rob Denroche (samStats.pl)
Modifications : Lawrence Heisler (bamqc.pl,bamqc.pm) 
<< <lheisler.oicr.on.ca> >>
Last modified : 2014-12-08
The Ontario Institute for Cancer Research

=cut

use strict;
use warnings;

use Getopt::Std;
use vars qw/ %opt /;
use File::Basename;

use GSI::bamqc;

use JSON::PP;    # imports encode_json, decode_json, to_json and from_json

### ASSESS options and generate the parameter hash
my $opt_string = "s:i:l:m:r:b:j:o:q:ch:H:D";
getopts( $opt_string, \%opt ) or usage("Incorrect arguments.");

### DEBUG, comment out when not in use
#use Data::Dumper;  ## will allow the code to be modified inline to show data dumps
#(open my $TTY,"/dev/tty") || die "unable to open keyboard input"; ## will allow code to be modified with <$TTY>, to hold for keyboard input.  <STDIN> will not work as this script reads from a stream

my %p = validate_opts(%opt);

### Check for streaming input
usage(
    "Missing Input stream.  This script must receive input from samtools view")
  if ( -t STDIN and not @ARGV );

### load the bed file, store in the parameter hash
$p{bed} = GSI::bamqc::read_bed( $p{bedFile} ) if ( $p{bedFile} );
my $targets_message = "Loaded " . $p{bed}{numberOfTargets} . " targets.\n\n";
if ( $p{bed} =~ /ERROR/ ) {
    print usage( $p{bed} );
}

### read histogram and markDuplicates input (if any) into parameter hash
if ( $p{histFile} ) {
    %{ $p{jsonHash}{"non collapsed bases covered"} } =
      GSI::bamqc::get_hist_cvg( $p{histFile}, "noncollapsed", 2000 )
      ;    ### add to the jsonhash
    %{ $p{jsonHash}{"collapsed bases covered"} } =
      GSI::bamqc::get_hist_cvg( $p{histFile}, "collapsed", 2000 );
}

if ( $p{markDupMetrics} ) {
    $p{markDuplicatesHash} = GSI::bamqc::read_markdup_metrics($p{markDupMetrics});
}

$p{sortedChars} = [
    split " ",
q(! " # $ % & ' \( \) * + , - . / 0 1 2 3 4 5 6 7 8 9 : ; < = > ? @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z [ \\ ] ^ _ ` a b c d e f g h i j k l m n o p q r s t u v w x y z { | } ~)
];

#<$TTY>;

### initalize entries to 0
### form the json hash at the end
### keep all information in stats hash
my %stats = map { ( $_, 0 ) } ( GSI::bamqc::bam_stats_keys() );

$stats{bed} = $p{bed};
### this will track the reads per startpoint
$stats{startPoint} = { count => 0, current => "null", RPSP => 1 };

my $timeToSampleCount
  ; ### this is a flag to identify when it is time to do analysis for the sampled read
while (<STDIN>)
{    ### need to check that STDIN exists, and run usage if it doesn't

    chomp;
    ### split the bam line and extract necessary fields
    my @f = split /\t/;

    my (
        $flag,   $chrom,           $start1,
        $qual,   $cigarstring,     $chrom2,
        $start2, $template_length, $qualString
    ) = @f[ 1, 2, 3, 4, 5, 6, 7, 8, 10 ];

    #### assess the flags
    my $mapped = GSI::bamqc::assess_flag( $flag, \%stats, $qual, $p{qualCut} );
    ### do not proceed with this record, unless mapped
    next unless ($mapped);

    ###  assess the start point of the mapped read,using the %startPoint hash which tracks this information
    $stats{startPoint} =
      GSI::bamqc::assess_start_point( $chrom, $start1, $start2,
        $stats{startPoint} );

    ### the remaining statistics are only collected on the sample data
    $timeToSampleCount++;
    my $sampleflag = 0;
    if ( ( $timeToSampleCount >= $p{sampleRate} ) ) {
        $timeToSampleCount = 0;    ### reset this value
        $stats{"sampled reads"}++;
        $sampleflag = 1;
    }
    next unless ($sampleflag);   ### don't proceed unless the sample flag is set

    my $mappedstrand = $flag & 16 ? "-" : "+";
    my $R = ( $flag & 64 ) ? "R1" : ( $flag & 128 ) ? "R2" : "R?";
    $qualString = reverse($qualString) if ( $mappedstrand eq "-" );
    my ( $readLength, $mapped_bases ) =
      GSI::bamqc::cigar_stats( $chrom, $start1, $start2, $R, $mappedstrand,
        $cigarstring, \%stats, \%p );    ### get cigar stats for this string

    #### mismatch statistics
    if ( $_ =~ /MD:Z:(.*?)\t/ ) {
        my $mmcount = GSI::bamqc::md_stats( $R, $1, $mappedstrand, \%stats );
    }
    else {
        $stats{readsMissingMDtags}++;
    }

    ### is the read onTarget
    my $onTarget =
      GSI::bamqc::onTarget( $chrom, $start1, $mapped_bases, \%stats );
    $stats{"reads on target"}++ if ($onTarget);

    #### RUNNING BASE COVERAGE
    if ( $onTarget && $p{reportBasesCovered} ) {
        my $added =
          GSI::bamqc::addRunningBaseCoverage( $chrom, $start1, $start2,
            $cigarstring, $mappedstrand, \%stats );
        my $cleared =
          GSI::bamqc::runningBaseCoverage( \%stats, $chrom, $start1 );
    }

    ### READ LENGTH HISTOGRAM
    $stats{readLengthHist}{$R}{$readLength}++;

    ### INSERT MAPPING DETAILS
    my $class =
      GSI::bamqc::insertMapping( $template_length, $chrom2, \%stats, \%p )
      if ( $template_length >= 0 );   ## ignore negative templates = 2nd in pair

    ## QUALITY,BY CYCLE hash
    my @qual = split //, $qualString;
    map {
        my $cyc = $_ + 1;
        $stats{qualByCycle}{$R}{ $qual[$_] }{$cyc}++

          #$stats{qualByCycle}{$R}{$cyc}{$qual[$_]}++

    } ( 0 .. $#qual );

}

$stats{ends} =
  $stats{"properly paired reads"} > 0 ? "paired end" : "single end";

######## Generation of summary stats from collected data ########

### the average read length is the sum of observed bases, characterized in different ways / # of mapped reads
$stats{averageReadLength}{overall} =
  $stats{"mapped reads"} > 0
  ? ( $stats{alignedCount} +
      $stats{softClipCount} +
      $stats{hardClipCount} +
      $stats{insertCount} ) / $stats{"mapped reads"} * $p{sampleRate}
  : 0;

$stats{"aligned bases"} =
  ( $stats{alignedCount} + $stats{insertCount} ) * $p{sampleRate};
$stats{"mismatch bases"}  = $stats{mismatchCount} * $p{sampleRate};
$stats{"inserted bases"}  = $stats{insertCount} * $p{sampleRate};
$stats{"deleted bases"}   = $stats{deletionCount} * $p{sampleRate};
$stats{"soft clip bases"} = $stats{softClipCount} * $p{sampleRate};
$stats{"hard clip bases"} = $stats{hardClipCount} * $p{sampleRate};

for my $R (qw/R1 R2 R?/) {
    my ( $sum, $count ) = ( 0, 0 );
    for my $l ( keys %{ $stats{readLengthHist}{$R} } ) {
        $sum   += $stats{readLengthHist}{$R}{$l} * $l;
        $count += $stats{readLengthHist}{$R}{$l};

        ### correct for sample rate
        $stats{readLengthHist}{$R}{$l} *= $p{sampleRate};
    }
    $stats{averageReadLength}{$R} = $count ? $sum / $count : 0;
}

if ( $stats{ends} eq "paired end" ) {
    ### stats on the Insert Sizes
#($stats{meanInsert},$stats{stdevInsert})=$stats{normalInsertSizes} ? HistStats($stats{normalInsertSizes}) : (0,0);
    ( $stats{meanInsert}, $stats{stdevInsert} ) =
      GSI::bamqc::HistStats( $stats{normalInsertSizes} );

    ### min and max insert sizes
    my @insertSizes = sort { $a <=> $b } keys %{ $stats{normalInsertSizes} };
    $stats{minInsert} = shift @insertSizes;
    $stats{maxInsert} = pop @insertSizes;
    for my $iSize ( $stats{minInsert} .. $stats{maxInsert} ) {
        $stats{normalInsertSizes}{$iSize} =
          $stats{normalInsertSizes}{$iSize} || 0;
        $stats{normalInsertSizes}{$iSize} *= $p{sampleRate};
    }
}

### determine the maximum read length for each read

for my $R (qw/R1 R2 R?/) {
    my @read_lengths = sort { $b <=> $a } keys %{ $stats{readLengthHist}{$R} };
    $stats{maxReadLength}{$R} = @read_lengths ? shift @read_lengths : 0;

    ### adjust the byCycle stats
    next
      unless ( $stats{maxReadLength}{$R} )
      ;    ### don't process unless this read has length

    #for my $key(keys %{$stats{byCycle}}){
    for my $key (qw/aligned mismatch insertion deletion softClip hardClip/)
    {      ### will process all if the read exists, setting to 0
        for my $cyc ( 1 .. $stats{maxReadLength}{$R} ) {
            $stats{byCycle}{$key}{$R}{$cyc} =
              $stats{byCycle}{$key}{$R}{$cyc} || 0;
            $stats{byCycle}{$key}{$R}{$cyc} *= $p{sampleRate};
        }
    }

    for my $cyc ( 1 .. $stats{maxReadLength}{$R} ) {
        my ( $qSum, $qCount ) = ( 0, 0 );

        #for my $q(keys %{$stats{qualByCycle}{$R}}){
        for my $q ( @{ $p{sortedChars} } ) {
            my $qscore = GSI::bamqc::toPhred($q);

            if ( my $count = $stats{qualByCycle}{$R}{$q}{$cyc} ) {
                $stats{qualHist}{$R}{$qscore} += $count;    # * $p{sampleRate};
                $qCount += $count;
                $qSum += $count * $qscore;
            }
        }
        $stats{qualLine}{$R}{$cyc} = $qCount ? int( $qSum / $qCount ) : 0;
    }

    ## correct qualHist for samplerate
    $stats{qualHist}{$R}{$_} *= $p{sampleRate}
      for keys %{ $stats{qualHist}{$R} };

}

### cut the number of significant figures down to 6 after the decimal point
$stats{startPoint}{RPSP} = ( sprintf "%.6f", $stats{startPoint}{RPSP} );
$stats{stdevInsert}      = ( sprintf "%.6f", $stats{stdevInsert} );
$stats{meanInsert}       = ( sprintf "%.6f", $stats{meanInsert} );
$stats{averageReadLength}{overall} =
  ( sprintf "%.6f", $stats{averageReadLength}{overall} );

### adjust sample rate on normalInsertSizes
#map{$stats{normalInsertSizes}{$_}*=$p{sampleRate}} keys %{$stats{normalInsertSizes}};

### clear out remaining runningBaseCoverage
my $cleared = GSI::bamqc::runningBaseCoverage( \%stats );

my $out;
if ($p{outputText}) {
    open $out, ">", $p{outputText} || die "Cannot open text output ".$p{outputText}.": $!";
} else {
    $out = *STDERR;
}
### write a human-readable text summary indicating various statistics

print $out $p{reportBasesMessage};
print $out $targets_message;
print $out "\n\n";

print $out "Total reads: " . $stats{"total reads"} . "\n";
print $out "Mapped reads: " . $stats{"mapped reads"} . "\n";
print $out "Non primary reads: " . $stats{"non primary reads"} . "\n";
print $out "MAPQ < $p{qualCut} reads: " . $stats{"qual fail reads"} . "\n";
print $out "Unmapped reads: " . $stats{"unmapped reads"} . "\n";

print $out "Sampled reads: " . $stats{"sampled reads"} . "\n";
print $out "Reads on target: " . $stats{"reads on target"} . "\n\n";

print $out "Reads missing MD tags!: " . $stats{"readsMissingMDtags"} . "\n\n";

print $out "Aligned bases: " . $stats{"alignedCount"} . "\n";
print $out "Soft clipped bases: " . $stats{softClipCount} . "\n";
print $out "Hard clipped bases: " . $stats{hardClipCount} . "\n";
print $out "Mismatched bases: " . $stats{mismatchCount} . "\n";
print $out "Deleted base count: " . $stats{deletionCount} . "\n";
print $out "Inserted base count: " . $stats{insertCount} . "\n";
print $out "Average read length: " . $stats{averageReadLength}{overall} . "\n\n";

print $out "Mean insert: " . $stats{meanInsert} . "\n";
print $out "Stdev insert: " . $stats{stdevInsert} . "\n";
print $out "Pairs with insert longer than $p{normalInsertMax}: "
  . $stats{pairsMappedAbnormallyFar} . "\n";
print $out "Pairs mapped to different chromosomes: "
  . $stats{pairsMappedToDifferentChr} . "\n\n";

print $out "Reads per start point: " . $stats{startPoint}{RPSP} . "\n\n";

if ($p{outputText}) {
    close $out || die "Cannot close text output ".$p{outputText}.": $!";
}

### this is adjusted after warnings, to only show the pre-correction on target numbers
$stats{"reads on target"} *= $p{sampleRate};

############  generate the jsonHash from sspecivicadd stats fields to the jsonHash, that will be encoded and printed out

my %jsonHash = GSI::bamqc::generate_jsonHash( \%stats, \%p );

### for debug
#for my $key(sort keys %jsonHash){   next unless($key=~/target/);
#	print STDERR "$key\n";<$TTY>;
#	print STDERR Dumper($jsonHash{$key});<$TTY>;
#}

print encode_json( \%jsonHash );

### SUBROUTINES specific to this script

sub validate_opts {
    my (%opt) = @_;
    usage("Help requested.") if ( $opt{h} );

    my %param;

    #default sampleRate is 1001 to catch both R1s and R2s more evenly
    $param{sampleRate}      = $opt{'s'} || 1001;
    $param{normalInsertMax} = $opt{'i'} || 1500;
    $param{outputText}      = $opt{'o'}; # if not defined, will default to STDERR
    $param{bedFile}         = $opt{'r'}
      || "/oicr/data/genomes/homo_sapiens/UCSC/Genomic/UCSC_hg19_random/hg19_random.genome.sizes.bed";
    $param{qualCut} = $opt{'q'} || 30;

    if ( $opt{'j'} ) {
        if ( -e $opt{'j'} ) {    #### is it a file
            my $jsonString;
            ( open my $JSON, "<", $opt{'j'} )
              || die "unable to open json formatted string in $opt{j}";
            while (<$JSON>) {
                chomp;
                $jsonString .= $_;
            }
            close $JSON;
            $param{jsonHash} = decode_json($jsonString);
        }
        elsif ( $opt{'j'} =~ /^\{.*\}$/ ) {

            $param{jsonHash} = decode_json( $opt{'j'} )
              ;    ### appears to be the json hash passed directly

        }
        else {
            usage("file containing json formatted string $opt{j} not found");
        }
    }

    #%{$param{jsonHash}}=%{ decode_json($opt{'j'}) } if($opt{'j'});

    if ( $opt{'m'} ) {
	my $path = $opt{'m'};
	if ( -e $path) {
	    $param{markDupMetrics} = $path;
	} else {
	    die "Duplicate metrics file $path does not exist.";
	}
    }

    if ( $opt{'b'} ) {
        $param{bamPath}                   = $opt{'b'};
        $param{jsonHash}{"bam path"}      = $opt{'b'};
        $param{jsonHash}{"last modified"} = ( stat( $opt{'b'} ) )[9];
    }
    if ( $opt{'c'} ) {
        $param{reportBasesCovered} = 1;
        $param{sampleRate}         = 1;
	$param{reportBasesMessage} = "Sampling every read and calculating bases covered.\n";
    }
    else {
        $param{reportBasesCovered} = 0;
	$param{reportBasesMessage} = "Only sampling every $param{sampleRate} reads.\n";
    }

    if ( $opt{'D'} ) {
        use Data::Dumper;
        ( open( $param{TTY} ), "/dev/tty" )
          || die "unable to open keyboard input"
          ; ## will allow code to be modified with <$TTY>, to hold for keyboard input.  <STDIN> will not work as this script reads from a stream
    }

    if ( $opt{'H'} )
    { ### extract additional information from a Bedtools coverage histogram file
        usage("bedtools histogram file $opt{H} not found") if ( !-e $opt{'H'} );
        $param{histFile} = $opt{'H'};

    }

    return %param;
}

sub usage {
    print "\nUsage is samtools view aligned.sorted.bam | bamqc.pl [options].\n";
    print "Bam file must be sorted by coordinates.\n";
    print "Options are as follows:\n";
    print
"\t-s sample rate.  Defines how often to sample the reads (default $p{sampleRate}).\n";
    print
"\t-i normal insert max.  Defines upper limit to what is considered a normal insert (default $p{normalInsertMax}).\n";
    print
"\t-m mark duplicates metric file. Text file output by Picard MarkDuplicates. Optional.";
    print
"\t-o path for human-readable summary text output. Optional, defaults to STDERR.\n";
    print
"\t-q mapping quality cut.  Reads that map with a quality worse than this will not be considered \"uniquely mapped\" (default $p{qualCut}).\n";
    print
"\t-r target.bed.  Bed file containing targets to calculate coverage against (default $p{bedFile}).\n";
    print
"\t\tNOTE: target.bed file MUST be sorted in the same order as the bam file.\n";
    print
"\t-c enables % base covered reporting. Sets sample rate to 1 (overrides -s), and runs long.\n";
    print "\t-b bam path (for recording path and timestamp).\n";
    print
"\t-j file containing additional JSON formatted data string. e.g. '{\"sample\":\"TLCR2C10n\",\"library\":\"TLCR_2C10_nn_n_PE_501_nn\",\"barcode\":\"TAGCTT\",\"instrument\":\"h802\",\"run name\":\"110616_SN802_0060_AC00FWACXX\",\"lane\":\"4\"}'\n";
    print "\t-H bedtools histogram file, for coverage analysis\n";
    print "\t-D debug mode\n";
    print "\t-h displays this usage message.\n";

    die "\n@_\n\n";
}

