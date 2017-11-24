package GSI::RunReport;
use strict;
use warnings;

BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

		@ISA		=	qw(Exporter);
		@EXPORT_OK		=	qw(get_possible_plot_names get_possible_headers plot_data
                        data_table coverage_table graph_table lane_info);
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
	my ($jsonHash,$scriptPath, $jsonDirs)=@_;
	for my $report (keys %$jsonHash){
	warn "graphing $report\n";
		my $rv=`$scriptPath/jsonToGraphs.pl $jsonDirs->{$report}/$report`;
		return $report;
	}
}

use SeqWare::Html;

sub data_table{
	my($p,$jsonHash,$sorted_lane_list,$run)=@_;
	# my $html="<table border=\"1\" class=\"sortable\">\n";

  my $header=get_header($p);


  my %stats;
  my @rows;
	for my $report(@$sorted_lane_list){
      my ($row_stats,@data_row)=data_row($p,$jsonHash->{$report});
      push( @rows, \@data_row);
		  map{$stats{$_}+=$$row_stats{$_}} qw/RawYield Reads EstYield/;
		  map{$stats{qualCuts}{$_}++} keys %{$$row_stats{qualCut}};
	}

  my @footer=get_footer($p,\%stats);



  my $html=SeqWare::Html::table($header, \@footer, @rows);

	my $qualCut = join(", ",sort{$b<=>$a} keys %{$stats{qualCuts}});

	$html .= SeqWare::Html::p("* Estimates exclude unmapped, off target, non-primary or ",
           "MAPQ < $qualCut reads ",
           $p->{noCollapse} ? "" : "and use reads per start point to approximate loss due to collapse"
           );

	return $html;
}

sub get_header{
    my($p)=@_;
    my @cols=@{$p->{table_columns}{data}};
    my @headings=@{$p->{table_headers}{data}}{@cols};
    return SeqWare::Html::tableHeader(@headings);
}

sub get_footer{
  my($p,$stats)=@_;
  my @cols=@{$p->{table_columns}{data}};
  map{ $stats->{$_}=~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g; } qw/RawYield Reads EstYield/;
  my %hash;
	map{$hash{$_}=""} @cols;
  @hash{($cols[0],"raw_reads","raw_yield","yield")}=("Total",$stats->{Reads},$stats->{RawYield},$stats->{EstYield});
	my @footings=@hash{@cols};
  return @footings;
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
  # add_td_tags($rowHash);

  $stats{qualCut}{$jsonHashReport->{"qual cut"}} = 1;

  # my $html="<tr>\n".join("\n",@{$rowHash}{@cols})."</tr>";
  my @row=@{$rowHash}{@cols};
  return (\%stats,@row);

}

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
       $rowHash->{$col}=sprintf ("%.2f", $rowHash->{$col}) if ($rowHash->{$col} ne $NA);
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

sub get_all_fields {
    my ($jsonHash) =@_;

    my %row=map{ $_ => $NA} (keys %{$table_headers{'data'}}, keys %plot_names);
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

    ### the remaining columns are all plots
    foreach (keys %plot_names) {
        $row{$_}="$jsonHash->{basename}.graphs/$plot_names{$_}";
    }

    return \%row;
}


1;
