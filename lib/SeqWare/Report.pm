package SeqWare::Report;
use strict;
use warnings;

BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
		@ISA		=	qw(Exporter);
    @EXPORT_OK   = qw();
}

sub new{
    my ($class, %parameters) = @_;
    my $self = bless ({}, ref ($class) || $class);
    return $self;
}


use SeqWare::Html;

sub commaFormat {
    my ($str) = @_;
    $str =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
    return $str;
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

# Compute a fraction from two values in a hash.
# 1[HashRef]: The hash to analyse.
# 2[Str]: The key containing the numerator.
# 3[Str]: The key containing the denominator.
# 4[Numeric or Undef]: A multiplication factor.
sub fractionFromHash {
    my ( $hash, $numerator, $denominator, $factor ) = @_;
    return ( $hash->{$numerator} // 0 ) *
      ( $factor // 1 ) /
      ( $hash->{$denominator} || 1 );

# denominator is different do squash strings that are "0", which // considers non-empty
}

# Generate the link name from a JSON file.
# 1[HashRef]: The JSON data.
sub nameForFile {
    my ($file) = @_;
    return join( "_",
        grep { defined $_ && $_ ne "NoIndex" }
          ( $file->{"run name"}, $file->{lane}, $file->{barcode} ) );
}

# Format a fraction as a percentage.
sub percent {
    my ($num) = @_;
    return sprintf( "%.2f%%", $num * 100 );
}

# Create a run summary table.
# 1[Bool]: Whether to show the “Run” column.
# 2[Bool]: Whether to make a paired-end report or a single-end only report
# 3[Bool]: Whether to show hard clipping results
# 3[Bool]: Whether to keep duplicate reads.
# ...[HashRef]: The runs to show.
# Returns[Str<HTML>]: The table
sub runSummary {
    my $showRun      = shift(@_);
    my $reportHasPE  = shift(@_);
    my $showHardClip = shift(@_);
    my $noCollapse   = shift(@_);

    my $totalReads    = 0;
    my $totalRawYield = 0;
    my $totalEstYield = 0;

    my @cols = (
        $showRun ? "Run" : undef,
        "Lane",
        "Barcode",
        "Library",
        $reportHasPE ? "Insert Mean (SD)" : undef,
        "Read Length",
        "PF Reads",
        "PF Yield",
        "Map %",
        $reportHasPE ? "R1 Mismatch %"  : "Mismatch %",
        $reportHasPE ? "R2 Mismatch %"  : undef,
        $reportHasPE ? "R1 Indel %"     : "Indel %",
        $reportHasPE ? "R2 Indel %"     : undef,
        $reportHasPE ? "R1 Soft Clip %" : "Soft Clip %",
        $reportHasPE ? "R2 Soft Clip %" : undef,
        $showHardClip
        ? ( $reportHasPE ? "R1 Hard Clip %" : "Hard Clip %" )
        : undef,
        $reportHasPE && $showHardClip ? "R2 Hard Clip %" : undef,
        "Reads/SP",
        "% mapped on Target",
        "Estimated Yield*",
        "Coverage*"
    );

    my @rows;
    for my $file (@_) {
        my $run_is_paired_end = $file->{"number of ends"} eq "paired end";

        my $rawReads =
          $file->{"mapped reads"} +
          $file->{"unmapped reads"} +
          $file->{"qual fail reads"};

        my $rawYield =
          int( $rawReads * $file->{"average read length"} )
          ; # cast to int because average length is only from mapped reads (and may result in ugly decimal)

        my $readLength;
        if ($run_is_paired_end) {
            $readLength =
                $file->{"read 1 average length"} . ", "
              . $file->{"read 2 average length"};
        }
        else {
            $readLength =
              SeqWare::Report::twoDP( $file->{"read ? average length"} );
        }

        my @rates;
        if ( $file->{"aligned bases"} > 0 ) {

            if ($run_is_paired_end) {
                @rates = (
                    GSI::bamqc::generate_mismatch_rate($file, 'read 1'),
                    GSI::bamqc::generate_mismatch_rate($file, 'read 2'),
                    GSI::bamqc::generate_indel_rate($file, 'read 1'),
                    GSI::bamqc::generate_indel_rate($file, 'read 2'),
                    GSI::bamqc::generate_softclip_rate($file, 'read 1'),
                    GSI::bamqc::generate_softclip_rate($file, 'read 2'),
                    $showHardClip ? GSI::bamqc::generate_hardclip_rate($file, 'read 1') : undef,
                    $showHardClip ? GSI::bamqc::generate_hardclip_rate($file, 'read 2') : undef,
                );
            }
            else    # single end
            {
                @rates = (
                    GSI::bamqc::generate_mismatch_rate( $file,""),
                    $reportHasPE ? "n/a" : undef,
                    GSI::bamqc::generate_indel_rate( $file,""),
                    $reportHasPE ? "n/a" : undef,
                    GSI::bamqc::generate_softclip_rate( $file,""),
                    $reportHasPE ? "n/a" : undef,
                    $showHardClip ? GSI::bamqc::generate_hardclip_rate( $file,"") : undef,
                    $reportHasPE && $showHardClip ? "n/a" : undef
                );
            }
        }
        else {
            @rates = (
                "n/a",
                $reportHasPE ? "n/a" : undef,
                "n/a",
                $reportHasPE ? "n/a" : undef,
                "n/a",
                $reportHasPE  ? "n/a" : undef,
                $showHardClip ? "n/a" : undef,
                $reportHasPE && $showHardClip ? "n/a" : undef
            );
        }

        my $onTargetRate =
          SeqWare::Report::fractionFromHash( $file, "reads on target",
            "mapped reads" );

        my $estimatedYield;
        if ($noCollapse) {
            $estimatedYield = int( $file->{"aligned bases"} * $onTargetRate );
        }
        else {
            $estimatedYield = int(
                SeqWare::Report::fractionFromHash(
                    $file,                   "aligned bases",
                    "reads per start point", $onTargetRate
                )
            );
        }

        $totalRawYield += $rawYield;
        $totalReads    += $rawReads;
        $totalEstYield += $estimatedYield;
        my $linkName = nameForFile($file);

        my @row = (
            $showRun ? $file->{'run name'} : undef,
            "<a href=\"#$linkName\">$file->{lane}</a>",
            exists $file->{barcode} && $file->{barcode} ne "NoIndex"
            ? "<a href=\"#$linkName\">$file->{barcode}</a>"
            : "none",
            "<a href=\"#$linkName\">$file->{library}</a>",
            $reportHasPE
            ? (
                $run_is_paired_end
                ?

                  SeqWare::Report::twoDP( $file->{"insert mean"} ) . " ("
                  . SeqWare::Report::twoDP( $file->{"insert stdev"} ) . ")"
                : "n/a"
              )
            : undef,
            $readLength,
            SeqWare::Report::commaFormat($rawReads),
            SeqWare::Report::commaFormat($rawYield),
            SeqWare::Report::percent(
                $file->{"mapped reads"} / ( $rawReads // 1 )
            ),
            @rates,
            SeqWare::Report::twoDP( $file->{"reads per start point"} ),
            SeqWare::Report::percent($onTargetRate),
            SeqWare::Report::commaFormat($estimatedYield),
            SeqWare::Report::twoDP(

                $estimatedYield / $file->{"target size"}
              )
              . "×"
        );
        push( @rows, [ grep { defined $_ } @row ] );
    }

    my @footer = (
        "Total",
        $showRun     ? "" : undef,
        $reportHasPE ? "" : undef,
        "",
        "",
        "",
        SeqWare::Report::commaFormat($totalReads),
        SeqWare::Report::commaFormat($totalRawYield),
        "",
        $reportHasPE  ? "" : undef,
        $reportHasPE  ? "" : undef,
        $reportHasPE  ? "" : undef,
        $showHardClip ? "" : undef,
        $reportHasPE && $showHardClip ? "" : undef,
        "",
        "",
        "",
        "",
        SeqWare::Report::commaFormat($totalEstYield),
        ""
    );

    return SeqWare::Html::table(
        SeqWare::Html::tableHeader( grep { defined $_ } @cols ),
        [ grep                           { defined $_ } @footer ], @rows );

}

# Format a number to have two decimal places.
sub twoDP {
    my ($num) = @_;
    return sprintf( "%.2f", $num );
}

1;
