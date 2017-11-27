package GSI::RunReport;
use strict;
use warnings;

BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

		@ISA		=	qw(Exporter);
		@EXPORT_OK		=	qw(get_possible_plot_names get_possible_headers plot_data
                        data_table coverage_table graph_table lane_info write_tsv);
}

my $NA="<a class='na'>na</a>";

my %table_headers=(
	data=>{
		lane		=>	"Lane",
		barcode		=>	"Barcode",
		groupid		=>	"Group ID",
		ext_name	=>	"External Name",
		library		=>	"Library",
		insert_mean	=>	"Insert Mean",
    ins_stddev  =>  "Insert Stddev",
		read_length	=>	"Read Length",
		raw_reads	=>	"PF Reads",
		raw_yield	=>	"PF Yield",
		mapped		=>	"Map %",
    mismatch1   =>  "R1 Mismatch %",
    mismatch2   =>  "R2 Mismatch %",
    indel1      =>  "R1 Indel %",
    indel2      =>  "R2 Indel %",
		error1		=>	"R1 Error %",
    error2      =>  "R2 Error %",
		softclip1	=>	"R1 Soft Clip %",
    softclip2   =>  "R2 Soft Clip %",
		hardclip1	=>  "R1 Hard Clip %",
    hardclip2    => "R2 Hard Clip %",
		rpsp		=>	"Reads/SP",
		ontarget	=>	"% Mapped on Target",
		yield		=>	"Estimated Yield*",
		coverage	=>	"Coverage*"
	},
  tsv=>{
    lane		=>	"Lane",
    barcode		=>	"Barcode",
    groupid		=>	"Group ID",
    ext_name	=>	"External Name",
    library		=>	"Library",
    insert_mean	=>	"Insert Mean",
    ins_stddev  =>  "Insert Stddev",
    read_length	=>	"Read Length",
    raw_reads	=>	"PF Reads",
    raw_yield	=>	"PF Yield",
    mapped		=>	"Map Percent",
    mismatch1   =>  "R1 Mismatch Percent",
    mismatch2   =>  "R2 Mismatch Percent",
    indel1      =>  "R1 Indel Percent",
    indel2      =>  "R2 Indel Percent",
    error1		=>	"R1 Error Percent",
    error2      =>  "R2 Error Percent",
    softclip1	=>	"R1 Soft Clip Percent",
    softclip2   =>  "R2 Soft Clip Percent",
    hardclip1	=>  "R1 Hard Clip Percent",
    hardclip2    => "R2 Hard Clip Percent",
    rpsp		=>	"Reads/SP",
    ontarget	=>	"Percent Mapped on Target",
    yield		=>	"Estimated Yield",
    coverage	=>	"Coverage",
    target_size => "Target Size (bp)",
    num_targets => "Number of Targets",
    target_file => "Target File"
  },
	graph=>{
		lane		=>	"Lane",
		barcode		=>	"Barcode",
		library		=>	"Library",
		read_breakdown	=> "Read Breakdown",
		insert_distr	=> "Insert Distribution",
		qual_hist	=> "Quality Histogram",
		qual_by_cycle	=> "Quality by Cycle",
		mismatch_by_cycle	=> "Mismatch by Cycle",
		indel_by_cycle	=> "Indels by Cycle",
		softclip_by_cycle	=> "Soft Clip by Cycle"
	},
	merge=>{
		sample		=>	"Sample",
		lanes		=>	"# Lanes",
		last_mod	=>	"Last Modified",
		avg_read_length	=>	"Avg Read Length",
		error_pc		=>	"Error %",
		soft_clip_pc	=>	"Soft Clip %",
		reads_on_target_pc	=>	"% Reads on Target",
		reads_on_target		=>	"# Reads on Target",
		bases_on_target		=>	"# Bases on Target",
		coverage	=>	"Coverage",
		coverage8x	=>	"% at 8x",
		coverage90	=>	"90% covered at",
		coverageGraph	=>	"Graph",
	}
);
my %plot_names=(
	read_breakdown		=> "readPie.png",
	insert_distr		=> "insert.png",
	read_length_dist	=> "readLength.png",
	qual_hist	 		=> "qualHist.png",
	qual_by_cycle 		=> "qualCycle.png",
	mismatch_by_cycle 	=> "misCycle.png",
	indel_by_cycle		=> "indelCycle.png",
	softclip_by_cycle	=> "softCycle.png",
	hardclip_by_cycle	=> "hardCycle.png"
);

###########################################


sub get_possible_headers {
    return \%table_headers;
}

sub get_possible_plot_names {
    return \%plot_names;
}

sub plot_data{
	my ($jsonHash,$scriptPath)=@_;
	for my $report (keys %$jsonHash){
	   warn "graphing $report\n";
     my $dest_dir=$jsonHash->{$report}->{dirname}."/".$report;
		 system("$scriptPath/jsonToGraphs.pl $dest_dir");
     if ($? != 0) {
       warn "graphing failed for $dest_dir: $?\n";
     }
	}
}

use SeqWare::Html;

sub data_table{
	my($p,$jsonHash,$sorted_lane_list,$run)=@_;
	# my $html="<table border=\"1\" class=\"sortable\">\n";
  my @cols=@{$p->{table_columns}{data}};

  #header
  my @headings=@{$p->{table_headers}{data}}{@cols};
  my $header=SeqWare::Html::tableHeader(@headings);

  #rows
  my %stats;
  my @rows;
	for my $report(@$sorted_lane_list){
      my ($row_stats,@data_row)=data_row($p,$jsonHash->{$report});
      push( @rows, \@data_row);
		  map{$stats{$_}+=$$row_stats{$_}} qw/RawYield Reads EstYield/;
		  map{$stats{qualCuts}{$_}++} keys %{$$row_stats{qualCut}};
	}

  #footer
  map{ $stats{$_}=~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g; } qw/RawYield Reads EstYield/;
  my %hash;
	map{$hash{$_}=""} @cols;
  @hash{($cols[0],"raw_reads","raw_yield","yield")}=("Total",$stats{Reads},$stats{RawYield},$stats{EstYield});
  my @footer=@hash{@cols};


  #shove it together
  my $html=SeqWare::Html::table($header, \@footer, @rows);

  # asterisk below the table
	my $qualCut = join(", ",sort{$b<=>$a} keys %{$stats{qualCuts}});
	$html .= SeqWare::Html::p("* Estimates exclude unmapped, off target, non-primary or ",
           "MAPQ < $qualCut reads ",
           $p->{noCollapse} ? "" : "and use reads per start point to approximate loss due to collapse"
           );

	return $html;
}

sub data_row{
  my($p,$jsonHashReport)=@_;
  my @cols=@{$p->{table_columns}{data}};
  ### the %row hash will hold the formatted contents of each possible cell in the row
  ### the %row hash key = column identifier
  ### the final list of cells to include in the row is @cols, extracted from %p

  my $rowHash = get_all_fields($jsonHashReport);

  ### these are cummulative totals, need to be added and returne....
  # as part of parameters?  or caputred
  my %stats=(RawYield=>$rowHash->{raw_yield},Reads=>$rowHash->{raw_reads},EstYield=>$rowHash->{yield});

  format_na($rowHash);
  format_big_numbers($rowHash);
  format_2decimals($rowHash);
  add_hrefs($rowHash);

  $stats{qualCut}{$jsonHashReport->{"qual cut"}} = 1;

  my @row=@{$rowHash}{@cols};
  return (\%stats,@row);

}

sub coverage_table{

	my($p,$jsonHash,$sorted_lane_list,$run)=@_;

	my @covX=@{$p->{coverageXs}};
	my $levels=scalar @covX;

  #header has all kinds of funny bits, so just make it here
	my $header="<tr><th \"col-md-1\" colspan=\"5\"></th>";
	$header.="<th \"col-md-1\" colspan=\"" . $levels . "\">Non-collapsed % bases covered</th>";
	$header.="<th \"col-md-1\" colspan=\"" . $levels . "\">Collapsed % bases covered</th>";
	$header.="</tr>\n";
	$header.="<tr>";
	$header.="<th \"col-md-1\">Lane</th>";
	$header.="<th \"col-md-1\">Barcode</th>";
	$header.="<th \"col-md-1\">Library</th>";
	$header.="<th \"col-md-1\">Target Size (bp)</th>";
	$header.="<th \"col-md-1\"># Targets</th>";
	my $covXheadings=join("",map{"<th>${_}x</th>"} @covX);
	$header.=$covXheadings.$covXheadings;
  $header.="</tr>";

  #rows
  my @rows;
  my @cols=( 'lane','barcode', 'library', 'target_size', 'num_targets');
	for my $report (@$sorted_lane_list){
        my $rowHash = get_all_fields($jsonHash->{$report});
        my $targetSize = $rowHash->{"target_size"};
        format_big_numbers($rowHash);
        add_hrefs($rowHash);
        my @row=@{$rowHash}{@cols};
        # $html.="<tr>\n";
        # foreach my $col (@cols){
        #     $html.="<td>$rowHash->{$col}</td>";
        # }
		    for my $collapsed ("non collapsed", "collapsed"){
			       for my $i (@covX){
				           my $basesCovered = exists $jsonHash->{$report}{"$collapsed bases covered"}{$i}
                        ? $jsonHash->{$report}{"$collapsed bases covered"}{$i} : 0;
				           $basesCovered = $basesCovered / $targetSize;
				           $basesCovered = sprintf "%.2f%%", $basesCovered * 100;
				           push (@row, $basesCovered);
			       }
		    }
        push (@rows, \@row);
	}

  return SeqWare::Html::table($header, undef, @rows);

}

sub graph_table{
	my($p,$jsonHash,$sorted_lane_list,$run)=@_;

	my @cols=@{$p->{table_columns}{graph}};
	my @headings=@{$p->{table_headers}{graph}}{@cols};
  my $header=SeqWare::Html::tableHeader(@headings);

	my @rows;
	for my $report (@$sorted_lane_list){
    my $rowHash = get_all_fields($jsonHash->{$report});
    format_na($rowHash);
    format_big_numbers($rowHash);
    format_2decimals($rowHash);
    add_hrefs($rowHash);
    make_thumbnail($rowHash);

    my @row=@{$rowHash}{@cols};
    push(@rows, \@row);
	}

  return SeqWare::Html::table($header, undef, @rows);
}

sub lane_info{
	my($p,$jsonHash,$sorted_lane_list,$run)=@_;
	my @html_lines;
	for my $report (@$sorted_lane_list){
		my $lane=$jsonHash->{$report}{"lane"};
		my $barcode_line=exists $jsonHash->{$report}{"barcode"} ? "<h2><a name=\"${run}_${lane}_$jsonHash->{$report}{\"barcode\"}\">$run Lane: $lane Barcode: $jsonHash->{$report}{\"barcode\"}</a></h2>"
														 : "<h2><a name=\"${run}_${lane}\">$run Lane: $lane</a></h2>";
		push(@html_lines,$barcode_line);
		push(@html_lines,"<p><a href=\"#$run\">Back to $run.</a></p>");
		push(@html_lines,"<ul>");
		push(@html_lines,"<li>Library: $jsonHash->{$report}{\"library\"}</li>");
		push(@html_lines,"<li>Target file: $jsonHash->{$report}{\"target file\"}</li>");
		my $targetSize = format_number($jsonHash->{$report}{"target size"});
		my $numberOfTargets = format_number($jsonHash->{$report}{"number of targets"});
		push(@html_lines,"<li>Target size: $targetSize bp</li>");
		push(@html_lines,"<li>Number of targets: $numberOfTargets</li>");

		push(@html_lines,"<li>Workflow Name: $jsonHash->{$report}{'workflow name'}</li>")
			if($jsonHash->{$report}{"workflow name"});
		push(@html_lines,"<li>Workflow Version: $jsonHash->{$report}{'workflow version'}</li>")
			if($jsonHash->{$report}{"workflow version"}) ;
		push(@html_lines,"<li>Workflow Run Creation Timestamp: $jsonHash->{$report}{'workflow run creation timestamp'}</li>")
			if($jsonHash->{$report}{"workflow run creation timestamp"});
        push(@html_lines,"<li>Bam File Creation Timestamp: $jsonHash->{$report}{'bam file creation timestamp'}</li>")
          	if($jsonHash->{$report}{"bam file creation timestamp"});
		push(@html_lines,"</ul>");

		if($p->{printAllImages}){
			for my $col(@{$p->{table_columns}{graph}}){  ### the columns of the graph table, show these plots if they exists
				if(my $plot=$p->{plotnames}{$col}){
					push(@html_lines,"<img src=\"${report}.graphs/$plot\"/>");
				}
			}
		}
    }
	my $html=join("\n",@html_lines);
	return $html;
}

=head2 write_tsv

  Writes a TSV with all "data_table" metrics plus the informationg from
  "lane_info". Arbitrarily picks the last JSON file's directory to write the
  TSV into.

=cut
sub write_tsv{
  	my($p,$jsonHash,$sorted_lane_list,$run)=@_;

  	my @cols=@{$p->{table_columns}{tsv}};

  	my @headings=@{$p->{table_headers}{tsv}}{@cols};
    my $header=join("\t",@headings);

    my $tsv=$header."\n";

    my $output_dir;
  	for my $report (@$sorted_lane_list){
      my $rowHash = get_all_fields($jsonHash->{$report});

      format_big_numbers($rowHash);
      format_2decimals($rowHash);

      my @row=@{$rowHash}{@cols};
      $tsv.=join("\t",@row)."\n";
      #arbitrarily picking the directory of the last JSON to write the TSV into
      $output_dir=$jsonHash->{$report}->{dirname};
  	}

    #handle the tsv file
  	my $tsv_file = $output_dir.'/'.$run.'_report.tsv';
  	open TSV, '>', $tsv_file or die "Couldn't open output tsv file: $!";
    print TSV "$tsv";
    close TSV;
    my $html = "<a href=\"$tsv_file\" download>${run}_report.tsv</a>\n";
    return $html;
}

#################################Useful junk

sub format_na {
    my $rowHash=shift;
    #if anything is na, grey it out
    map { $rowHash->{$_}=$NA if ($rowHash->{$_} eq 'na') } keys %$rowHash;
    return $rowHash;
}
sub format_2decimals {
    my $rowHash=shift;
    #format floats to have two decimal places
    my @printf_cols = ('mapped', 'softclip1', 'softclip2', 'hardclip1', 'hardclip2', 'mismatch1', 'mismatch2',
                       'error1', 'error2', 'rpsp', 'ontarget', 'coverage', 'insert_mean', 'ins_stddev', 'indel1', 'indel2');
    foreach my $col (@printf_cols) {
       $rowHash->{$col}=sprintf ("%.2f", $rowHash->{$col})
          if (($rowHash->{$col} ne $NA) || ($rowHash->{$col} ne 'na')) ;
    }
    return $rowHash;
}

sub format_big_numbers {
    my $rowHash=shift;
    #add commas to big numbers
    my @format_number_cols = ('raw_reads', 'raw_yield', 'yield','target_size','num_targets');
    foreach (@format_number_cols) {
       $rowHash->{$_}=format_number($rowHash->{$_});
    }
    return $rowHash;
}

sub format_number{
	my ($n)=@_;
	$n=~s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
	return $n;
}

sub add_hrefs {
    my $rowHash=shift;
    ### add hyperlinks to the required cells
	my $linkName=exists $rowHash->{barcode}
           ? $rowHash->{"run_name"} ."_$rowHash->{lane}_$rowHash->{barcode}"
           : $rowHash->{"run_name"}."_$rowHash->{lane}";

    my @linked_columns = ( 'lane', 'barcode', 'library' ) ;
    foreach (@linked_columns) {
        $rowHash->{$_}="<a href=\"#$linkName\">$rowHash->{$_}</a>";
    }
    return $rowHash;
}

sub add_td_tags {
    my $rowHash=shift;
    #add table tags
    map { $rowHash->{$_}="<td>$rowHash->{$_}</td>" } keys %$rowHash;
    return $rowHash;
}

sub make_thumbnail {
    my $rowHash=shift;
    map {
        $rowHash->{$_}="<a href=\"$rowHash->{$_}\"><img src=\"$rowHash->{$_}\" width=\"100\" height=\"100\"/></a> "
    } keys %plot_names;

}


sub get_all_fields {
    my ($jsonHash) =@_;

    my %row=map{ $_ => "na"} (keys %{$table_headers{'data'}}, keys %plot_names);
    $row{run_name}=$jsonHash->{"run name"};
    $row{lane}=$jsonHash->{lane};
    $row{barcode}=GSI::bamqc::get_barcode($jsonHash);
    $row{groupid}=GSI::bamqc::get_group($jsonHash);
    $row{ext_name}=$jsonHash->{'external name'} if exists $jsonHash->{"external name"};
    $row{library}=$jsonHash->{library};
    $row{raw_reads}=GSI::bamqc::get_raw_reads($jsonHash);
    $row{raw_yield}=GSI::bamqc::get_raw_yield($jsonHash);
    $row{read_length}=$jsonHash->{"number of ends"} eq "paired end"
                       ? $jsonHash->{"read 1 average length"} . ", " . $jsonHash->{"read 2 average length"}
                       : sprintf "%.2f", $jsonHash->{"read ? average length"};
	  $row{mapped}=GSI::bamqc::get_map_percent($jsonHash);
    $row{softclip1}=GSI::bamqc::generate_softclip_rate($jsonHash, "read 1");
    $row{softclip2}=GSI::bamqc::generate_softclip_rate($jsonHash, "read 2");

    $row{indel1}=GSI::bamqc::generate_indel_rate($jsonHash, "read 1");
    $row{indel2}=GSI::bamqc::generate_indel_rate($jsonHash, "read 2");

    $row{hardclip1}=GSI::bamqc::generate_hardclip_rate($jsonHash, "read 1");
    $row{hardclip2}=GSI::bamqc::generate_hardclip_rate($jsonHash, "read 2");

    $row{mismatch1}=GSI::bamqc::generate_mismatch_rate($jsonHash, "read 1");
    $row{mismatch2}=GSI::bamqc::generate_mismatch_rate($jsonHash, "read 2");
    $row{error1}=GSI::bamqc::generate_error_rate($jsonHash, "read 1");
    $row{error2}=GSI::bamqc::generate_error_rate($jsonHash, "read 2");

    $row{rpsp}=$jsonHash->{"reads per start point"};
    $row{ontarget}=GSI::bamqc::get_ontarget_percent($jsonHash);
    $row{yield}=GSI::bamqc::get_est_yield($jsonHash);
    $row{coverage}=GSI::bamqc::get_est_coverage($jsonHash);
    $row{insert_mean}=$jsonHash->{"insert mean"} if $jsonHash->{"number of ends"} eq "paired end";
    $row{ins_stddev}=$jsonHash->{"insert stdev"} if $jsonHash->{"number of ends"} eq "paired end";

    $row{target_size}=$jsonHash->{"target size"};
    $row{num_targets} = $jsonHash->{"number of targets"};
    $row{target_file} = $jsonHash->{"target file"};

    ### the remaining columns are all plots
    foreach (keys %plot_names) {
        $row{$_}="$jsonHash->{basename}.graphs/$plot_names{$_}";
    }

    return \%row;
}


1;
