package GSI::RunReport;
use strict;
use warnings;
use SeqWare::Html;

BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

    @ISA       = qw(Exporter);
    @EXPORT_OK = qw(get_possible_plot_names get_possible_headers plot_data
      data_table coverage_table graph_table lane_info write_tsv
      assemble_run_report);
}

sub new {
    my ( $class, %parameters ) = @_;
    my $self = bless( {}, ref($class) || $class );
    return $self;
}

#################### main pod documentation begin ###################

=head1 NAME

RunReport - Use the JSON produced by GSI::bamqc to generate an HTML page

=head1 SYNOPSIS

  use GSI::RunReport;


=head1 DESCRIPTION


=head1 USAGE

use SeqWare::Html;

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

Genome Sequence Informatics
Ontario Institute for Cancer Research
https://github.com/oicr-gsi

=head1 COPYRIGHT

Copyright (C) 2017 The Ontario Institute for Cancer Research

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

=head1 SEE ALSO

perl(1).

=cut

#################### main pod documentation end ###################

my $NA = "<a class='na'>na</a>";

# Mapping the short names to the header titles. See get_possible_headers
my %table_headers = (
    data => {
        lane               => "Lane",
        barcode            => "Barcode",
        groupid            => "Group ID",
        ext_name           => "External Name",
        library            => "Library",
        insert_mean        => "Insert Mean",
        ins_stddev         => "Insert Stddev",
        read_length        => "Read Length",
        raw_reads          => "PF Reads",
        raw_yield          => "PF Yield",
        mapped             => "Map %",
        mismatch1          => "R1 Mismatch %",
        mismatch2          => "R2 Mismatch %",
        indel1             => "R1 Indel %",
        indel2             => "R2 Indel %",
        error1             => "R1 Error %",
        error2             => "R2 Error %",
        softclip1          => "R1 Soft Clip %",
        softclip2          => "R2 Soft Clip %",
        hardclip1          => "R1 Hard Clip %",
        hardclip2          => "R2 Hard Clip %",
        rpsp               => "Reads/SP",
        ontarget           => "% Mapped on Target",
        yield              => "Estimated Yield*",
        collapsed_yield    => "Estimated Yield*",
        coverage           => "Coverage*",
        collapsed_coverage => "Coverage*"

    },
    tsv => {
        lane               => "Lane",
        barcode            => "Barcode",
        groupid            => "Group ID",
        ext_name           => "External Name",
        library            => "Library",
        insert_mean        => "Insert Mean",
        ins_stddev         => "Insert Stddev",
        read_length        => "Read Length",
        raw_reads          => "PF Reads",
        raw_yield          => "PF Yield",
        mapped             => "Map Percent",
        mismatch1          => "R1 Mismatch Percent",
        mismatch2          => "R2 Mismatch Percent",
        indel1             => "R1 Indel Percent",
        indel2             => "R2 Indel Percent",
        error1             => "R1 Error Percent",
        error2             => "R2 Error Percent",
        softclip1          => "R1 Soft Clip Percent",
        softclip2          => "R2 Soft Clip Percent",
        hardclip1          => "R1 Hard Clip Percent",
        hardclip2          => "R2 Hard Clip Percent",
        rpsp               => "Reads/SP",
        ontarget           => "Percent Mapped on Target",
        yield              => "Estimated Yield (uncollapsed)",
        collapsed_yield    => "Estimated Yield (collapsed)",
        coverage           => "Coverage (uncollapsed)",
        collapsed_coverage => "Coverage (collapsed)",
        target_size        => "Target Size (bp)",
        num_targets        => "Number of Targets",
        target_file        => "Target File"
    },
    graph => {
        lane              => "Lane",
        barcode           => "Barcode",
        library           => "Library",
        read_breakdown    => "Read Breakdown",
        insert_distr      => "Insert Distribution",
        qual_hist         => "Quality Histogram",
        qual_by_cycle     => "Quality by Cycle",
        mismatch_by_cycle => "Mismatch by Cycle",
        indel_by_cycle    => "Indels by Cycle",
        softclip_by_cycle => "Soft Clip by Cycle"
    }
);
my %plot_names = (
    read_breakdown    => "readPie.png",
    insert_distr      => "insert.png",
    read_length_dist  => "readLength.png",
    qual_hist         => "qualHist.png",
    qual_by_cycle     => "qualCycle.png",
    mismatch_by_cycle => "misCycle.png",
    indel_by_cycle    => "indelCycle.png",
    softclip_by_cycle => "softCycle.png",
    hardclip_by_cycle => "hardCycle.png"
);

###########################################

=head2 get_possible_headers()

  Returns the headers that the RunReport supports, mapped to their human-friendly
  name. The keys in this hash are used by the various reporting methods to create
  their tables. In order to get tables to print different sets of headers
  change the parameter passed in: p{table_columns}{data|graph|tsv}.

  Returns: a hash mapping shortnames to human-friendly names for the table header

=cut

sub get_possible_headers {
    return \%table_headers;
}

=head2 get_possible_plot_names()

  Returns the image short names that the RunReport supports, mapped to their
  actual PNG filename. The keys in this hash are used by the various reporting
  methods to create their tables. In order to get tables to print different sets of headers
  change the parameter passed in: p{plotnames}.

  Returns: a hash mapping shortnames to plot filenames for the table header

=cut

sub get_possible_plot_names {
    return \%plot_names;
}

=head2 plot_data($jsonHash,$scriptPath)

  Creates the graphs for the report using jsonToGraphs.pl in the same directory
  as each JSON file.

  Arguments: $jsonHash = The hash containing the JSON file contents to analyse;
             $scriptPath = the directory where jsonToGraphs.pl is located.

=cut

sub plot_data {
    my ( $jsonHash, $scriptPath ) = @_;
    for my $report ( keys %$jsonHash ) {
        warn "graphing $report\n";
        my $dest_dir = $jsonHash->{$report}->{dirname} . "/" . $report;
        system("$scriptPath/jsonToGraphs.pl $dest_dir");
        if ( $? != 0 ) {
            warn "graphing failed for $dest_dir: $?\n";
        }
    }
}

=head2 data_table($p,$jsonHash,$sorted_lane_list,$run)

  Returns the HTML for the main metrics table.

  Arguments: $p=parameters used for the data_table. Required keys:
                $p->{table_columns}{data}-the ordered array of shortname for the header
                $p->{noCollapse}-whether the estimated yield and coverage should be collapsed by RPSP;
             $jsonHash=A hash containing the JSON files contents to analyse;
             $sorted_lane_list=the list of json files in sorted lane order;
             $run=the name of the sequencer run.

  Returns: the HTML for the main metrics table.

  See: get_possible_headers

=cut

sub data_table {
    my ( $p, $jsonHash, $sorted_lane_list, $run ) = @_;
    my @cols = @{ $p->{table_columns}{data} };

    #header
    my @headings = @{ $p->{table_headers}{data} }{@cols};
    my $header   = SeqWare::Html::tableHeader(@headings);

    #rows
    my %stats;
    my @rows;
    for my $report (@$sorted_lane_list) {
        my ( $row_stats, @data_row ) = data_row( $p, $jsonHash->{$report} );
        push( @rows, \@data_row );
        map { $stats{$_} += $$row_stats{$_} } qw/RawYield Reads EstYield/;
        map { $stats{qualCuts}{$_}++ } keys %{ $$row_stats{qualCut} };
    }

    #footer
    map {
        $stats{$_} =~
          s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
    } qw/RawYield Reads EstYield/;
    my %hash;
    map { $hash{$_} = "" } @cols;
    @hash{
        (
            $cols[0], "raw_reads",
            "raw_yield", ( $p->{noCollapse} ? 'yield' : 'collapsed_yield' )
        )
    } = ( "Total", $stats{Reads}, $stats{RawYield}, $stats{EstYield} );
    my @footer = @hash{@cols};

    #shove it together
    my $html = SeqWare::Html::table( $header, \@footer, @rows );

    # asterisk below the table
    my $qualCut = join( ", ", sort { $b <=> $a } keys %{ $stats{qualCuts} } );
    $html .= SeqWare::Html::p(
        "* Estimates exclude unmapped, off target, non-primary or ",
        "MAPQ < $qualCut reads ",
        $p->{noCollapse}
        ? "(uncollapsed)"
        : "and use reads per start point to approximate loss due to collapse"
    );

    return $html;
}

sub data_row {
    my ( $p, $jsonHashReport ) = @_;
    my @cols = @{ $p->{table_columns}{data} };
    ### the %row hash will hold the formatted contents of each possible cell in the row
    ### the %row hash key = column identifier
    ### the final list of cells to include in the row is @cols, extracted from %p

    my $rowHash = get_all_fields($jsonHashReport);

    ### these are cummulative totals, need to be added and returne....
    # as part of parameters?  or caputred
    my %stats = (
        RawYield => $rowHash->{raw_yield},
        Reads    => $rowHash->{raw_reads},
        EstYield => (
            $p->{noCollapse} ? $rowHash->{yield} : $rowHash->{collapsed_yield}
        )
    );

    format_na($rowHash);
    format_big_numbers($rowHash);
    format_2decimals($rowHash);
    add_hrefs($rowHash);

    $stats{qualCut}{ $jsonHashReport->{"qual cut"} } = 1;

    my @row = @{$rowHash}{@cols};
    return ( \%stats, @row );

}

=head2 coverage_table($p,$jsonHash,$sorted_lane_list,$run)

  Returns the HTML for the coverage table.

  Arguments: $p=parameters used for the coverage_table. Required keys:
                $p->{coverageXs}-an array with the N X values to report for coverage
             $jsonHash=A hash containing the JSON files contents to analyse;
             $sorted_lane_list=the list of json files in sorted lane order;
             $run=the name of the sequencer run.

  Returns: the HTML for the coverage table.


=cut

sub coverage_table {

    my ( $p, $jsonHash, $sorted_lane_list, $run ) = @_;

    my @covX   = @{ $p->{coverageXs} };
    my $levels = scalar @covX;

    #header has all kinds of funny bits, so just make it here
    my $header = "<tr><th \"col-md-1\" colspan=\"5\"></th>";
    $header .=
        "<th \"col-md-1\" colspan=\""
      . $levels
      . "\">Non-collapsed % bases covered</th>";
    $header .=
        "<th \"col-md-1\" colspan=\""
      . $levels
      . "\">Collapsed % bases covered</th>";
    $header .= "</tr>\n";
    $header .= "<tr>";
    $header .= "<th \"col-md-1\">Lane</th>";
    $header .= "<th \"col-md-1\">Barcode</th>";
    $header .= "<th \"col-md-1\">Library</th>";
    $header .= "<th \"col-md-1\">Target Size (bp)</th>";
    $header .= "<th \"col-md-1\"># Targets</th>";
    my $covXheadings = join( "", map { "<th>${_}x</th>" } @covX );
    $header .= $covXheadings . $covXheadings;
    $header .= "</tr>";

    #rows
    my @rows;
    my @cols = ( 'lane', 'barcode', 'library', 'target_size', 'num_targets' );
    for my $report (@$sorted_lane_list) {
        my $rowHash    = get_all_fields( $jsonHash->{$report} );
        my $targetSize = $rowHash->{"target_size"};
        format_big_numbers($rowHash);
        add_hrefs($rowHash);
        my @row = @{$rowHash}{@cols};
        for my $collapsed ( "non collapsed", "collapsed" ) {
            for my $i (@covX) {
                my $basesCovered =
                  exists $jsonHash->{$report}{"$collapsed bases covered"}{$i}
                  ? $jsonHash->{$report}{"$collapsed bases covered"}{$i}
                  : 0;
                $basesCovered = $basesCovered / $targetSize;
                $basesCovered = sprintf "%.2f%%", $basesCovered * 100;
                push( @row, $basesCovered );
            }
        }
        push( @rows, \@row );
    }

    return SeqWare::Html::table( $header, undef, @rows );

}

=head2 graph_table($p,$jsonHash,$sorted_lane_list,$run)

  Returns the HTML for the graph table. Note that this does not do the graphing.
  plot_data must be called to generate graphs.

  Arguments: $p=parameters used for the graph_table. Required keys:
                  $p->{table_columns}{graph}-the ordered array of shortname for the header
             $jsonHash=A hash containing the JSON files contents to analyse;
             $sorted_lane_list=the list of json files in sorted lane order;
             $run=the name of the sequencer run.

  Returns: the HTML for the graph table.

  See: plot_data, get_possible_headers

=cut

sub graph_table {
    my ( $p, $jsonHash, $sorted_lane_list, $run ) = @_;

    my @cols     = @{ $p->{table_columns}{graph} };
    my @headings = @{ $p->{table_headers}{graph} }{@cols};
    my $header   = SeqWare::Html::tableHeader(@headings);

    my @rows;
    for my $report (@$sorted_lane_list) {
        my $rowHash = get_all_fields( $jsonHash->{$report} );
        format_na($rowHash);
        format_big_numbers($rowHash);
        format_2decimals($rowHash);
        add_hrefs($rowHash);
        make_thumbnail($rowHash);

        my @row = @{$rowHash}{@cols};
        push( @rows, \@row );
    }

    return SeqWare::Html::table( $header, undef, @rows );
}

=head2 lane_info($p,$jsonHash,$sorted_lane_list,$run)

  Returns the HTML for the target block.

  Arguments: $p=parameters used for the graph_table. Required keys:
                  $p->{table_columns}{graph}-the ordered array of shortname for the graph table header
                  $p->{printAllImages}-whether or not to print all of the images in large
             $jsonHash=A hash containing the JSON files contents to analyse;
             $sorted_lane_list=the list of json files in sorted lane order;
             $run=the name of the sequencer run.

  Returns: the HTML for the lane block.

  See: get_possible_headers

=cut

sub lane_info {
    my ( $p, $jsonHash, $sorted_lane_list, $run ) = @_;
    my @html_lines;
    for my $report (@$sorted_lane_list) {
        my $lane = $jsonHash->{$report}{"lane"};
        my $barcode_line =
          exists $jsonHash->{$report}{"barcode"}
          ? "<h2><a name=\"${run}_${lane}_$jsonHash->{$report}{\"barcode\"}\">$run Lane: $lane Barcode: $jsonHash->{$report}{\"barcode\"}</a></h2>"
          : "<h2><a name=\"${run}_${lane}\">$run Lane: $lane</a></h2>";
        push( @html_lines, $barcode_line );
        push( @html_lines, "<p><a href=\"#$run\">Back to $run.</a></p>" );
        push( @html_lines, "<ul>" );
        push( @html_lines,
            "<li>Library: $jsonHash->{$report}{\"library\"}</li>" );
        push( @html_lines,
            "<li>Target file: $jsonHash->{$report}{\"target file\"}</li>" );
        my $targetSize = format_number( $jsonHash->{$report}{"target size"} );
        my $numberOfTargets =
          format_number( $jsonHash->{$report}{"number of targets"} );
        push( @html_lines, "<li>Target size: $targetSize bp</li>" );
        push( @html_lines, "<li>Number of targets: $numberOfTargets</li>" );

        push( @html_lines,
            "<li>Workflow Name: $jsonHash->{$report}{'workflow name'}</li>" )
          if ( $jsonHash->{$report}{"workflow name"} );
        push( @html_lines,
"<li>Workflow Version: $jsonHash->{$report}{'workflow version'}</li>"
        ) if ( $jsonHash->{$report}{"workflow version"} );
        push( @html_lines,
"<li>Workflow Run Creation Timestamp: $jsonHash->{$report}{'workflow run creation timestamp'}</li>"
        ) if ( $jsonHash->{$report}{"workflow run creation timestamp"} );
        push( @html_lines,
"<li>Bam File Creation Timestamp: $jsonHash->{$report}{'bam file creation timestamp'}</li>"
        ) if ( $jsonHash->{$report}{"bam file creation timestamp"} );
        push( @html_lines, "</ul>" );

        if ( $p->{printAllImages} ) {
            for my $col ( @{ $p->{table_columns}{graph} } )
            { ### the columns of the graph table, show these plots if they exists
                if ( my $plot = $p->{plotnames}{$col} ) {
                    push( @html_lines,
                        "<img src=\"${report}.graphs/$plot\"/>" );
                }
            }
        }
    }
    my $html = join( "\n", @html_lines );
    return $html;
}

=head2 write_tsv($p,$jsonHash,$sorted_lane_list,$run)

  Writes a TSV with all "data_table" metrics plus the information from
  "lane_info". Arbitrarily picks the last JSON file's directory to write the
  TSV into.

  Arguments: $p=parameters used for the data_table. Required keys:
                $p->{table_columns}{tsv}-the ordered array of shortname for the header
             $jsonHash=A hash containing the JSON files contents to analyse;
             $sorted_lane_list=the list of json files in sorted lane order;
             $run=the name of the sequencer run.

  Returns: the HTML pointing to the TSV file

  See: get_possible_headers

=cut

sub write_tsv {
    my ( $p, $jsonHash, $sorted_lane_list, $run ) = @_;

    my @cols = @{ $p->{table_columns}{tsv} };

    my @headings = @{ $p->{table_headers}{tsv} }{@cols};
    my $header = join( "\t", @headings );

    my $tsv = $header . "\n";

    my $output_dir;
    for my $report (@$sorted_lane_list) {
        my $rowHash = get_all_fields( $jsonHash->{$report} );

        format_big_numbers($rowHash);
        format_2decimals($rowHash);

        my @row = @{$rowHash}{@cols};
        $tsv .= join( "\t", @row ) . "\n";

       #arbitrarily picking the directory of the last JSON to write the TSV into
        $output_dir = $jsonHash->{$report}->{dirname};
    }

    #handle the tsv file
    my $tsv_file = $output_dir . '/' . $run . '_report.tsv';
    open TSV, '>', $tsv_file or die "Couldn't open output tsv file: $!";
    print TSV "$tsv";
    close TSV;
    my $html = "<a href=\"$tsv_file\" download>${run}_report.tsv</a>\n";
    return $html;
}

sub assemble_run_report {
    my ($html) = @_;
    my $date = `date`;    #TODO(apmasell): do this without fork
    chomp $date;
    return SeqWare::Html::document( "$date Generic Run Report",
        SeqWare::Html::NAV_PROJECTS, "../../../web/seqwareBrowser", $html,
        get_custom_head() );

}

sub get_custom_head {
    use Cwd qw(abs_path);
    use File::Basename qw(dirname);
    my $scriptPath = dirname( abs_path(__FILE__) );
    use File::Slurp;
    my $htmlHead = "<style type=\"text/css\">\n";
    $htmlHead .= read_file( $scriptPath . "/css/runreport.css" );
    $htmlHead .= "</style>\n";
    return $htmlHead;
}

#################################Useful methods

# Finds the string 'na' and formats it with an HTML class
#if anything is na, grey it out
sub format_na {
    my $rowHash = shift;
    map { $rowHash->{$_} = $NA if ( $rowHash->{$_} eq 'na' ) } keys %$rowHash;
    return $rowHash;
}

#format floats from specific headers to have two decimal places
sub format_2decimals {
    my $rowHash     = shift;
    my @printf_cols = (
        'mapped',             'softclip1',
        'softclip2',          'hardclip1',
        'hardclip2',          'mismatch1',
        'mismatch2',          'error1',
        'error2',             'rpsp',
        'ontarget',           'coverage',
        'collapsed_coverage', 'insert_mean',
        'ins_stddev',         'indel1',
        'indel2'
    );
    foreach my $col (@printf_cols) {
        $rowHash->{$col} = sprintf( "%.2f", $rowHash->{$col} )
          if ( ( $rowHash->{$col} ne $NA ) || ( $rowHash->{$col} ne 'na' ) );
    }
    return $rowHash;
}

#add commas to big numbers in specific headers
sub format_big_numbers {
    my $rowHash            = shift;
    my @format_number_cols = (
        'raw_reads',   'raw_yield', 'yield', 'collapsed_yield',
        'target_size', 'num_targets'
    );
    foreach (@format_number_cols) {
        $rowHash->{$_} = format_number( $rowHash->{$_} );
    }
    return $rowHash;
}

#formats other big numbers
sub format_number {
    my ($n) = @_;
    $n =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
    return $n;
}

#add hyperlinks for run name, barcode, lane to go to the lane info table
sub add_hrefs {
    my $rowHash = shift;
    ### add hyperlinks to the required cells
    my $linkName =
      exists $rowHash->{barcode}
      ? $rowHash->{"run_name"} . "_$rowHash->{lane}_$rowHash->{barcode}"
      : $rowHash->{"run_name"} . "_$rowHash->{lane}";

    my @linked_columns = ( 'lane', 'barcode', 'library' );
    foreach (@linked_columns) {
        $rowHash->{$_} = "<a href=\"#$linkName\">$rowHash->{$_}</a>";
    }
    return $rowHash;
}

#add a <td></td> tag to everything in a list
sub add_td_tags {
    my $rowHash = shift;

    #add table tags
    map { $rowHash->{$_} = "<td>$rowHash->{$_}</td>" } keys %$rowHash;
    return $rowHash;
}

#for all of the plots in %plot_names, iterate through a hash and format any
#matching string to be a thumbnail image.
sub make_thumbnail {
    my $rowHash = shift;
    map {
        $rowHash->{$_} =
"<a href=\"$rowHash->{$_}\"><img src=\"$rowHash->{$_}\" width=\"100\" height=\"100\"/></a> "
    } keys %plot_names;

}

# map all of the possible shortnames to their actual data. We regenerate this any
# time we need this information, and then filter based on the specific header of
# interest.
sub get_all_fields {
    my ($jsonHash) = @_;

    my %row =
      map { $_ => "na" } ( keys %{ $table_headers{'data'} }, keys %plot_names );
    $row{run_name} = $jsonHash->{"run name"};
    $row{lane}     = $jsonHash->{lane};
    $row{barcode}  = GSI::bamqc::get_barcode($jsonHash);
    $row{groupid}  = GSI::bamqc::get_group($jsonHash);
    $row{ext_name} = $jsonHash->{'external name'}
      if exists $jsonHash->{"external name"};
    $row{library}   = $jsonHash->{library};
    $row{raw_reads} = GSI::bamqc::get_raw_reads($jsonHash);
    $row{raw_yield} = GSI::bamqc::get_raw_yield($jsonHash);
    $row{read_length} =
        $jsonHash->{"number of ends"} eq "paired end"
      ? $jsonHash->{"read 1 average length"} . ", "
      . $jsonHash->{"read 2 average length"}
      : sprintf "%.2f", $jsonHash->{"read ? average length"};
    $row{mapped}    = GSI::bamqc::get_map_percent($jsonHash);
    $row{softclip1} = GSI::bamqc::generate_softclip_rate( $jsonHash, "read 1" );
    $row{softclip2} = GSI::bamqc::generate_softclip_rate( $jsonHash, "read 2" );

    $row{indel1} = GSI::bamqc::generate_indel_rate( $jsonHash, "read 1" );
    $row{indel2} = GSI::bamqc::generate_indel_rate( $jsonHash, "read 2" );

    $row{hardclip1} = GSI::bamqc::generate_hardclip_rate( $jsonHash, "read 1" );
    $row{hardclip2} = GSI::bamqc::generate_hardclip_rate( $jsonHash, "read 2" );

    $row{mismatch1} = GSI::bamqc::generate_mismatch_rate( $jsonHash, "read 1" );
    $row{mismatch2} = GSI::bamqc::generate_mismatch_rate( $jsonHash, "read 2" );
    $row{error1} = GSI::bamqc::generate_error_rate( $jsonHash, "read 1" );
    $row{error2} = GSI::bamqc::generate_error_rate( $jsonHash, "read 2" );

    $row{rpsp}               = $jsonHash->{"reads per start point"};
    $row{ontarget}           = GSI::bamqc::get_ontarget_percent($jsonHash);
    $row{yield}              = GSI::bamqc::get_est_yield( $jsonHash, 0 );
    $row{coverage}           = GSI::bamqc::get_est_coverage( $jsonHash, 0 );
    $row{collapsed_yield}    = GSI::bamqc::get_est_yield( $jsonHash, 1 );
    $row{collapsed_coverage} = GSI::bamqc::get_est_coverage( $jsonHash, 1 );
    $row{insert_mean}        = $jsonHash->{"insert mean"}
      if $jsonHash->{"number of ends"} eq "paired end";
    $row{ins_stddev} = $jsonHash->{"insert stdev"}
      if $jsonHash->{"number of ends"} eq "paired end";

    $row{target_size} = $jsonHash->{"target size"};
    $row{num_targets} = $jsonHash->{"number of targets"};
    $row{target_file} = $jsonHash->{"target file"};

    ### the remaining columns are all plots
    foreach ( keys %plot_names ) {
        $row{$_} = "$jsonHash->{basename}.graphs/$plot_names{$_}";
    }

    return \%row;
}

1;
