package GSI::jsonToGraphs;

#Original Script : Rob Denroche
#Modifications : Lawrence Heisler << <lheisler.oicr.on.ca> >>, Iain Bancarz << <ibancarz.oicr.on.ca> >>
#Last modified : 2019-10-09
#Updated to process inputs from BamQC Niassa workflow 3.0 as well as 2.x
#
# Copyright (C) 2015-2019 The Ontario Institute for Cancer Research
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

use strict;
use warnings;

BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

    $VERSION = 1.00;
    @ISA     = qw(Exporter);
    @EXPORT_OK =
      qw(readmap_piechart quality_histogram collapsed_base_coverage noncollapsed_base_coverage
      readlength_histogram insert_graph quality_by_cycle mismatch_by_cycle
      indel_by_cycle softclip_by_cycle hardclip_by_cycle coverage_by_depth);
}

sub new {
    my ( $class, %parameters ) = @_;
    my $self = bless( {}, ref($class) || $class );
    return $self;
}

use Data::Dumper;

#### coverage by depth chart, requires that the json has contains coverage information
sub coverage_by_depth {
    my ( $jHash, $p, $path, $title ) = @_;

    #print STDERR Dumper(keys %$jHash);
    #<STDIN>;
    #print STDERR Dumper($$jHash{"collapsed bases covered"});
    #<STDIN>;

}

##### pieChart, read mapping
sub readmap_piechart {
    my ( $jHash, $p, $path, $title ) = @_;
    my %pieChart = ( values => [], labels => [], colours => [] );
    $pieChart{values} = [
        $$jHash{"reads on target"},
        $$jHash{"mapped reads"} - $$jHash{"reads on target"},
        $$jHash{"qual fail reads"},
        $$jHash{"unmapped reads"}
    ];
    $pieChart{labels} =
      [ "on target", "off target", "repeat/low quality", "unmapped" ];
    $pieChart{colours} = [
        $$p{colours}{good},     $$p{colours}{mid},
        $$p{colours}{other}[0], $$p{colours}{bad}
    ];
    pieGraph( $path, "readPie.png", $pieChart{values}, $pieChart{labels},
        $pieChart{colours}, "$title Read Breakdown" );
    return 1;

}

##### barChart, quality histogram
sub quality_histogram {
    my ( $jHash, $p, $path, $title ) = @_;
    my %barChart = ( values => [], labels => [], colours => [] );
    for my $r ( @{ $$p{reads} } ) {
        my %histHash = %{ $$jHash{"$r quality histogram"} };
        my @tempArray = sort { $a <=> $b } keys %histHash;
        if ( ( scalar @tempArray ) > 0 )    ### do not proceed if no keys
        {
            my $histMax = $tempArray[$#tempArray];    ### the maximum value
            if ( $histMax > 0 ) {
                for ( my $i = 0 ; $i <= $histMax ; $i++ ) {
                    my $val = $histHash{$i} || 0;
                    push( @{ $barChart{values} }, $val );
                    push( @{ $barChart{labels} }, $i );
                    my $colour_code =
                        $i < $$p{qualcut}{low}  ? "bad"
                      : $i < $$p{qualcut}{high} ? "mid"
                      :                           "good";
                    push( @{ $barChart{colours} }, $$p{colours}{$colour_code} );
                }
            }
        }
    }
    barGraph(
        $path,              "qualHist.png",
        $barChart{values},  $barChart{labels},
        $barChart{colours}, "$title Quality Histogram",
        "Base Quality",     "# Bases",
        ""
    );
}

##### barChart, non-collapsed bases covered
sub noncollapsed_base_coverage {
    my ( $jHash, $p, $path, $title ) = @_;
    my %barChart = ( values => [], labels => [], colours => [] );
    my $targetSize = $$jHash{"target size"};
    if ( defined $$jHash{"non collapsed bases covered"} ) {
        my %histHash = %{ $$jHash{"non collapsed bases covered"} };

        if ( ( scalar keys %histHash ) > 0 ) {
            for ( my $i = 1 ; $i <= $$p{basecoverage}{displaymax} ; $i++ ) {

#my $percentCovered=exists $histHash{$i} ? ($histHash{$i} / $targetSize) * 100 : 0;
                my $percentCovered = $histHash{$i} ? $histHash{$i} : 0;

                push( @{ $barChart{values} }, $percentCovered );
                push( @{ $barChart{labels} }, "\"${i}x\"" );
                my $colourCode =
                    $percentCovered < $$p{basecoverage}{low}  ? "bad"
                  : $percentCovered < $$p{basecoverage}{high} ? "mid"
                  :                                             "good";
                push( @{ $barChart{colours} }, $$p{colours}{$colourCode} );
            }
        }
        barGraph(
            $path,              "nonCollapsedCover.png",
            $barChart{values},  $barChart{labels},
            $barChart{colours}, "$title Non Collapsed Coverage",
            "Coverage Depth",   "% Target Covered",
            ", ylim=c(0, 100)"
        );
        return 1;
    }
    else {
        return 0;
    }
}

##### collapsed base coverage
sub collapsed_base_coverage {
    my ( $jHash, $p, $path, $title ) = @_;
    my $targetSize = $$jHash{"target size"};
    my ( @values, @labels, @colours );
    if ( defined $$jHash{"collapsed bases covered"} ) {
        my %histHash = %{ $$jHash{"collapsed bases covered"} };

        if ( ( scalar keys %histHash ) > 0 ) {
            for ( my $i = 1 ; $i <= $$p{basecoverage}{displaymax} ; $i++ ) {

#my $percentCovered=exists $histHash{$i} ? ($histHash{$i} / $targetSize) * 100 : 0;
                my $percentCovered = $histHash{$i} ? $histHash{$i} : 0;

                push( @values, $percentCovered );
                push( @labels, "\"${i}x\"" );

                my $colour_code =
                    $percentCovered < $$p{basecoverage}{low}  ? "bad"
                  : $percentCovered < $$p{basecoverage}{high} ? "mid"
                  :                                             "good";
                push( @colours, $$p{colours}{$colour_code} );
            }
        }
        barGraph(
            $path,            "collapsedCover.png",
            \@values,         \@labels,
            \@colours,        "$title Collapsed Coverage",
            "Coverage Depth", "% Target Covered",
            ", ylim=c(0, 100)"
        );
        return 1;
    }
    else {
        return 0;
    }
}

##### read length histogram
sub readlength_histogram {
    my ( $jHash, $p, $path, $title ) = @_;
    my %lineHash = ();
    my ( @lineX, @lineY, @colours );
    for my $r ( @{ $$p{reads} } ) {
        if ( defined $$jHash{"$r length histogram"} ) {
            for my $l ( keys %{ $$jHash{"$r length histogram"} } ) {
                $lineHash{$l} += $$jHash{"$r length histogram"}{$l};
            }
        }
        for my $l ( sort { $a <=> $b } keys %lineHash ) {
            push( @lineX,   $l );
            push( @lineY,   $lineHash{$l} );
            push( @colours, $$p{colours}{good} );
        }
    }
    lineGraph(
        $path, "readLength.png", \@lineX, \@lineY, \@colours,
        "$title Read Length Histogram",
        "Read Length (bp)",
        "Number of reads", ""
    );
    return 1;
}

##### insert size distribution
sub insert_graph {
    my ( $jHash, $p, $path, $title ) = @_;
    my $is_single = 0;
    if (defined($$jHash{"paired end"}) && $$jHash{"paired end"} == 0) {
	$is_single = 1; # 3.0+
    } elsif (defined($$jHash{"number of ends"}) && $$jHash{"number of ends"} eq "single end") {
	$is_single = 1; # 2.x
    }
    return 0 if ($is_single);
    # insert graph
    my @lineX      = ();
    my @lineY      = ();
    my @colours    = ();
    my $insertMean = $$jHash{"insert mean"};
    my %lineHash;
    if (defined $$jHash{"insert histogram"}) {
	%lineHash   = %{ $$jHash{"insert histogram"} }; # 2.x
    } else {
	%lineHash   = %{ $$jHash{"insert size histogram"} }; # 3.0+
    }

    for my $i ( sort { $a <=> $b } keys %lineHash ) {
        if ( $i < $$p{insertMax} ) {
            push( @lineX, $i );
            push( @lineY, $lineHash{$i} );
            my $insertStep = $$p{insertStep};
            my $colour_key =
              (
                     ( $i < ( $insertMean - ( 2 * $insertStep ) ) )
                  or ( $i > ( $insertMean + ( 2 * $insertStep ) ) )
              ) ? "bad"
              : (    ( $i < ( $insertMean - ( 1 * $insertStep ) ) )
                  or ( $i > ( $insertMean + ( 1 * $insertStep ) ) ) ) ? "mid"
              : "good";
            push( @colours, $$p{colours}{$colour_key} );
        }
    }
    lineGraph(
        $path, "insert.png", \@lineX, \@lineY, \@colours,
        "$title Insert Distribution",
        "Insert Size (bp)",
        "Pairs", ""
    );
    return 1;
}

### note, these "by cycle" graphs can likely be a single function, with a specific argument provided to indicate what should be plotted.
### perhaps create a general by_cycle plot, that each can call with specific information

##### quality by cycle graph
sub quality_by_cycle {
    my ( $jHash, $p, $path, $title ) = @_;
    my ( @lineX, @lineY, @colours );
    my $read1max = 0;
    for my $r ( @{ $$p{reads} } ) {
        if ( defined $$jHash{"$r quality by cycle"} ) {
            my %lineHash = %{ $$jHash{"$r quality by cycle"} };
            for my $cyc ( sort { $a <=> $b } keys %lineHash ) {
                my $qual = $lineHash{$cyc};
                if ( $$jHash{"number of ends"} eq "single end" ) {
                    push( @lineX, $cyc );
                    push( @lineY, $qual );
                }
                else {
                    if ( $r eq "read 1" ) {
                        push( @lineX, $cyc );
                        push( @lineY, $qual );
                        $read1max++
                          ; ### increment the read max, this will determine cycle for read 2
                    }
                    elsif ( $r eq "read 2" ) {
                        push( @lineX, $cyc + $read1max );
                        push( @lineY, $qual );
                    }
                }
                my $colour_key =
                    $qual < $$p{qualcut}{low}  ? "bad"
                  : $qual < $$p{qualcut}{high} ? "mid"
                  :                              "good";
                push( @colours, $$p{colours}{$colour_key} );

            }
        }
        if ( $r eq "read 1"
          ) ### do this after goingthrough the quality by cycle for read 1.  introduces a gap
        {
            push( @lineX,   $read1max );
            push( @lineY,   0 );
            push( @colours, "white" );
            $read1max++;
        }
    }
    lineGraph( $path, "qualCycle.png", \@lineX, \@lineY, \@colours,
        "$title Quality by Cycle",
        "Cycle", "Mean Quality", ", ylim=c(0, $$p{qualLineMax})" );
    return 1;
}

#####
sub mismatch_by_cycle {
    my ( $jHash, $p, $path, $title ) = @_;
    my ( @lineX, @lineY, @colours ) = ();
    my $read1max = 0;
    my ( %alignedHash, %insertHash, %softClipHash );
    for my $r ( @{ $$p{reads} } ) {
        $alignedHash{$r}  = \%{ $$jHash{"$r aligned by cycle"} };
        $insertHash{$r}   = \%{ $$jHash{"$r insertion by cycle"} };
        $softClipHash{$r} = \%{ $$jHash{"$r soft clip by cycle"} };

        if ( ( scalar keys %{ $$jHash{"$r mismatch by cycle"} } ) > 0 ) {
            my %errorHash = %{ $$jHash{"$r mismatch by cycle"} };
            for my $cyc ( sort { $a <=> $b } keys %errorHash ) {
                if ( $$jHash{"number of ends"} eq "single end" ) {
                    push( @lineX, $cyc );
                }
                else {
                    if ( $r eq "read 1" ) {
                        push( @lineX, $cyc );
                        $read1max++;
                    }
                    elsif ( $r eq "read 2" ) {
                        push( @lineX, $cyc + $read1max );
                    }
                }
                if ( ( $alignedHash{$r}{$cyc} + $insertHash{$r}{$cyc} ) > 0 ) {
                    push(
                        @lineY,
                        (
                            (
                                $errorHash{$cyc} / (
                                    $alignedHash{$r}{$cyc} +
                                      $insertHash{$r}{$cyc}
                                )
                            ) * 100
                        )
                    );
                }
                else {
                    push( @lineY, 0 );
                }
                push( @colours, $$p{colours}{bad} );
            }
            if ( $r eq "read 1" ) {
                push( @lineX,   $read1max );
                push( @lineY,   0 );
                push( @colours, "white" );
                $read1max++;
            }
        }
    }
    lineGraph(
        $path, "misCycle.png", \@lineX, \@lineY, \@colours,
        "$title Mismatches by Cycle",
        "Cycle",
        "% Bases Mismatched",
        ", ylim=c(0, 10)"
    );
    return 1;
}

##### indel by cycle
sub indel_by_cycle {
    my ( $jHash, $p, $path, $title ) = @_;
    my ( @lineX, @lineY, @colours );
    my $read1max = 0;
    my ( %alignedHash, %insertHash, %softClipHash );
    for my $r ( @{ $$p{reads} } ) {
        $alignedHash{$r}  = \%{ $$jHash{"$r aligned by cycle"} };
        $insertHash{$r}   = \%{ $$jHash{"$r insertion by cycle"} };
        $softClipHash{$r} = \%{ $$jHash{"$r soft clip by cycle"} };
        if ( ( scalar keys %{ $$jHash{"$r deletion by cycle"} } ) > 0 ) {
            my %errorHash = %{ $$jHash{"$r deletion by cycle"} };
            for my $i ( sort { $a <=> $b } keys %errorHash ) {
                if ( $$jHash{"number of ends"} eq "single end" ) {
                    push( @lineX, $i );
                }
                else {
                    if ( $r eq "read 1" ) {
                        push( @lineX, $i );
                        $read1max++;
                    }
                    elsif ( $r eq "read 2" ) {
                        push( @lineX, $i + $read1max );
                    }
                }
                if ( ( $alignedHash{$r}{$i} + $insertHash{$r}{$i} ) > 0 ) {
                    push(
                        @lineY,
                        (
                            (
                                ( $errorHash{$i} + $insertHash{$r}{$i} ) / (
                                    $alignedHash{$r}{$i} + $insertHash{$r}{$i}
                                )
                            ) * 100
                        )
                    );
                }
                else {
                    push( @lineY, 0 );
                }
                push( @colours, $$p{colours}{bad} );
            }
            if ( $r eq "read 1" ) {
                push( @lineX,   $read1max );
                push( @lineY,   0 );
                push( @colours, "white" );
                $read1max++;
            }
        }
    }
    lineGraph(
        $path,     "indelCycle.png",
        \@lineX,   \@lineY,
        \@colours, "$title Indels by Cycle",
        "Cycle",   "% Bases Inserted/Deleted",
        ", ylim=c(0, 10)"
    );
    return 1;
}

##### soft clip by cycle
sub softclip_by_cycle {
    my ( $jHash, $p, $path, $title ) = @_;
    my ( @lineX, @lineY, @colours );
    my $read1max = 0;

    my ( %alignedHash, %insertHash, %softClipHash );
    for my $r ( @{ $$p{reads} } ) {
        $alignedHash{$r}  = \%{ $$jHash{"$r aligned by cycle"} };
        $insertHash{$r}   = \%{ $$jHash{"$r insertion by cycle"} };
        $softClipHash{$r} = \%{ $$jHash{"$r soft clip by cycle"} };

        if ( ( scalar keys %{ $$jHash{"$r soft clip by cycle"} } ) > 0 ) {
            my %errorHash = %{ $$jHash{"$r soft clip by cycle"} };
            for my $i ( sort { $a <=> $b } keys %errorHash ) {
                if ( $$jHash{"number of ends"} eq "single end" ) {
                    push( @lineX, $i );
                }
                else {
                    if ( $r eq "read 1" ) {
                        push( @lineX, $i );
                        $read1max++;
                    }
                    elsif ( $r eq "read 2" ) {
                        push( @lineX, $i + $read1max );
                    }
                }
                if (
                    (
                        $alignedHash{$r}{$i} +
                        $insertHash{$r}{$i} +
                        $errorHash{$i}
                    ) > 0
                  )
                {
                    push(
                        @lineY,
                        (
                            (
                                $errorHash{$i} / (
                                    $alignedHash{$r}{$i} +
                                      $errorHash{$i} +
                                      $insertHash{$r}{$i}
                                )
                            ) * 100
                        )
                    );
                }
                else {
                    push( @lineY, 0 );
                }
                push( @colours, $$p{colours}{bad} );
            }
            if ( $r eq "read 1" ) {
                push( @lineX,   $read1max );
                push( @lineY,   0 );
                push( @colours, "white" );
                $read1max++;
            }
        }
    }
    lineGraph(
        $path, "softCycle.png", \@lineX, \@lineY, \@colours,
        "$title Soft Clips by Cycle",
        "Cycle",
        "% Bases Soft Clipped",
        ", ylim=c(0, 100)"
    );
    return 1;
}

##### hard clip by cycle
sub hardclip_by_cycle {
    my ( $jHash, $p, $path, $title ) = @_;
    my ( @lineX, @lineY, @colours );
    my $read1max = 0;

    my ( %alignedHash, %insertHash, %softClipHash );
    for my $r ( @{ $$p{reads} } ) {
        $alignedHash{$r}  = \%{ $$jHash{"$r aligned by cycle"} };
        $insertHash{$r}   = \%{ $$jHash{"$r insertion by cycle"} };
        $softClipHash{$r} = \%{ $$jHash{"$r soft clip by cycle"} };
        if ( ( scalar keys %{ $$jHash{"$r hard clip by cycle"} } ) > 0 ) {
            my %errorHash = %{ $$jHash{"$r hard clip by cycle"} };
            for my $i ( sort { $a <=> $b } keys %errorHash ) {
                if ( $$jHash{"number of ends"} eq "single end" ) {
                    push( @lineX, $i );
                }
                else {
                    if ( $r eq "read 1" ) {
                        push( @lineX, $i );
                        $read1max++;
                    }
                    elsif ( $r eq "read 2" ) {
                        push( @lineX, $i + $read1max );
                    }
                }
                if (
                    (
                        $alignedHash{$r}{$i} +
                        $insertHash{$r}{$i} +
                        $errorHash{$i}
                    ) > 0
                  )
                {
                    push(
                        @lineY,
                        (
                            (
                                $errorHash{$i} / (
                                    $alignedHash{$r}{$i} +
                                      $errorHash{$i} +
                                      $insertHash{$r}{$i} +
                                      $softClipHash{$r}{$i}
                                )
                            ) * 100
                        )
                    );
                }
                else {
                    push( @lineY, 0 );
                }
                push( @colours, $$p{colours}{bad} );
            }
            if ( $r eq "read 1" ) {
                push( @lineX,   $read1max );
                push( @lineY,   0 );
                push( @colours, "white" );
                $read1max++;
            }
        }
    }
    lineGraph(
        $path, "hardCycle.png", \@lineX, \@lineY, \@colours,
        "$title Hard Clips by Cycle",
        "Cycle",
        "% Bases Hard Clipped",
        ", ylim=c(0, 100)"
    );
    return 1;
}

######
###  GENERIC plotting scripts, called by the SPECIFIC plots

sub pieGraph {
    my $path    = $_[0];
    my $name    = $_[1];
    my @values  = @{ $_[2] };
    my @labels  = @{ $_[3] };
    my @colours = @{ $_[4] };
    my $title   = $_[5];

    open( RFILE, ">${path}/${name}.Rcode" )
      or die "Couldn't open ${path}/${name}.Rcode.\n";

    print RFILE "slices <- c($values[0]";
    for ( my $i = 1 ; $i < scalar @values ; $i++ ) {
        print RFILE ", $values[$i]";
    }
    print RFILE ")\n";

    print RFILE "lbls <- c(\"$labels[0]\"";
    for ( my $i = 1 ; $i < scalar @labels ; $i++ ) {
        print RFILE ", \"$labels[$i]\"";
    }
    print RFILE ")\n";

    print RFILE "cols <- c(\"$colours[0]\"";
    for ( my $i = 1 ; $i < scalar @colours ; $i++ ) {
        print RFILE ", \"$colours[$i]\"";
    }
    print RFILE ")\n";

    print RFILE "pct <- round(slices/sum(slices)*100)\n";
    print RFILE "lbls <- paste(lbls, pct)\n";
    print RFILE "lbls <- paste(lbls,\"%\",sep=\"\")\n";

    print RFILE
      "png(filename = \"${path}/${name}\", width = 640, height = 640)\n";
    print RFILE
      "pie(slices, labels = lbls, col=cols, main=\"$title\", border=NA)\n";
    print RFILE "dev.off()\n";

    close RFILE;

    `Rscript ${path}/${name}.Rcode`;
}

sub lineGraph {
    my $path             = $_[0];
    my $name             = $_[1];
    my @xVal             = @{ $_[2] };
    my @yVal             = @{ $_[3] };
    my @colours          = @{ $_[4] };
    my $title            = $_[5];
    my $xlab             = $_[6];
    my $ylab             = $_[7];
    my $additionalParams = $_[8];

    open( RFILE, ">${path}/${name}.Rcode" )
      or die "Couldn't open ${path}/${name}.Rcode.\n";

    print RFILE "xvals <- c($xVal[0]";
    for ( my $i = 1 ; $i < scalar @xVal ; $i++ ) {
        print RFILE ", $xVal[$i]";
    }
    print RFILE ")\n";

    print RFILE "yvals <- c($yVal[0]";
    for ( my $i = 1 ; $i < scalar @yVal ; $i++ ) {
        print RFILE ", $yVal[$i]";
    }
    print RFILE ")\n";

    print RFILE "cols <- c(\"$colours[0]\"";
    for ( my $i = 1 ; $i < scalar @colours ; $i++ ) {
        print RFILE ", \"$colours[$i]\"";
    }
    print RFILE ")\n";

    print RFILE
      "png(filename = \"${path}/${name}\", width = 640, height = 640)\n";
    print RFILE
"plot(xvals, yvals, main=\"$title\", type=\"n\", col=\"black\", xlab=\"$xlab\", ylab=\"$ylab\"$additionalParams)\n";
    print RFILE "for(i in 1:(length(yvals)-1))\n{\n";

#	print RFILE "polygon(c(xvals[i], xvals[i], mean(c(xvals[i], xvals[i+1])), mean(c(xvals[i], xvals[i+1]))), c(0, yvals[i], yvals[i], 0), col=cols[i], border=NA)\n";   # c(0, yvals[i], mean(c(yvals[i], yvals[i+1])), 0), col=cols[i], border=NA)\n";
#	print RFILE "polygon(c(mean(c(xvals[i], xvals[i+1])), mean(c(xvals[i], xvals[i+1])), xvals[i+1], xvals[i+1]), c(0, yvals[i+1], yvals[i+1], 0), col=cols[i+1], border=NA)\n";  # c(0, mean(c(yvals[i], yvals[i+1])), yvals[i+1], 0), col=cols[i+1], border=NA)\n";
    print RFILE
"polygon(c(xvals[i] - 0.5, xvals[i] - 0.5, xvals[i] + 0.5, xvals[i] + 0.5), c(0, yvals[i], yvals[i], 0), col=cols[i], border=NA)\n"
      ; # c(0, yvals[i], mean(c(yvals[i], yvals[i+1])), 0), col=cols[i], border=NA)\n";
    print RFILE "}\n";
    print RFILE "dev.off()\n";

    close RFILE;
    `Rscript ${path}/${name}.Rcode`;
}

sub barGraph {
    my $path             = $_[0];
    my $name             = $_[1];
    my @values           = @{ $_[2] };
    my @labels           = @{ $_[3] };
    my @colours          = @{ $_[4] };
    my $title            = $_[5];
    my $xlab             = $_[6];
    my $ylab             = $_[7];
    my $additionalParams = $_[8];

    open( RFILE, ">${path}/${name}.Rcode" )
      or die "Couldn't open ${path}/${name}.Rcode.\n";

    print RFILE "values <- c($values[0]";
    for ( my $i = 1 ; $i < scalar @values ; $i++ ) {
        print RFILE ", $values[$i]";
    }
    print RFILE ")\n";

    print RFILE "labels <- c($labels[0]";
    for ( my $i = 1 ; $i < scalar @labels ; $i++ ) {
        print RFILE ", $labels[$i]";
    }
    print RFILE ")\n";

    print RFILE "cols <- c(\"$colours[0]\"";
    for ( my $i = 1 ; $i < scalar @colours ; $i++ ) {
        print RFILE ", \"$colours[$i]\"";
    }
    print RFILE ")\n";

    print RFILE
      "png(filename = \"${path}/${name}\", width = 640, height = 640)\n";
    print RFILE
"barplot(values, col=cols, names.arg=labels, main=\"$title\", xlab=\"$xlab\", ylab=\"$ylab\", border=cols$additionalParams)\n";
    print RFILE "dev.off()\n";

    close RFILE;
    `Rscript ${path}/${name}.Rcode`;
}

1;
