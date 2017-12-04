#!/usr/bin/perl
#Original Script : Rob Denroche
#Modifications : Genome Sequence Informatics https://github.com/oicr-gsi
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

use strict;
use warnings;

use Cwd qw(abs_path);
use Getopt::Std;
use File::Basename;
use vars qw/ %opt /;

use GSI::bamqc 'load_json';
use GSI::RunReport;

use Data::Dumper;


my $scriptPath = dirname(abs_path($0));

my $plot_names=GSI::RunReport::get_possible_plot_names();
my $table_headers = GSI::RunReport::get_possible_headers();

my %param=(

	#### ordered data_table_columns, modify to limit to specific columns
	table_headers=>$table_headers,
	table_columns=>{   #### will show all possible columns, removing those that need to be explicitly indicated
		data	=> [qw/lane barcode groupid ext_name library insert_mean ins_stddev read_length raw_reads raw_yield mapped mismatch1 mismatch2 indel1 indel2 softclip1 softclip2 rpsp ontarget collapsed_yield collapsed_coverage/],
		graph	=> [qw/lane barcode library read_breakdown insert_distr qual_hist qual_by_cycle mismatch_by_cycle indel_by_cycle softclip_by_cycle hardclip_by_cycle/],
		tsv	=> [qw/lane barcode groupid ext_name library insert_mean ins_stddev read_length raw_reads raw_yield mapped mismatch1 mismatch2 indel1 indel2 softclip1 softclip2 rpsp ontarget collapsed_yield collapsed_coverage target_size num_targets target_file/]
	},
	plotnames=>$plot_names,
	coverageXs => [qw(1 4 8 15 30 50 100 200)],
	showDataTable 		=> 1,
	showGraphTable 		=> 1,
	showCoverageTable 	=> 0,
	showLaneInfo 		=> 1,
	printAllImages 		=> 0,
	plotData 			=> 0,
	noCollapse		=>0
);






#### process options
my $opt_string = "crpgnHh";
getopts ($opt_string, \%opt) or usage("Incorrect arguments.");
usage("Help requested.") if (exists $opt{h});

$param{showCoverageTable} = 1 if (exists $opt{c});
$param{printAllImages} = 1 if (exists $opt{p});
$param{plotData} = 1 if(exists $opt{g});

if (exists $opt{n}) {
	$param{noCollapse} = 1;
	#replace the collapsed value with the uncollapsed value
	foreach my $table (keys %{$param{table_columns}}) {
			my @cols=@{$param{table_columns}{$table}};
			s/collapsed_// for @cols;
			$param{table_columns}{$table}=\@cols;
	}
}

if(! exists $opt{r}){
	### remove the read length histogram
	@{$param{table_columns}{graph}}	=	grep(!/^insert_distr$/, 		@{$param{table_columns}{graph}}	);

}
if(! exists $opt{H}){
	### remove hard clipping info
	@{$param{table_columns}{data}}	=	grep(!/^hardclip$/,				@{$param{table_columns}{data}}	);
	@{$param{table_columns}{graph}}	=	grep(!/^hardclip_by_cycle$/,	@{$param{table_columns}{graph}}	);
}


#### get list of files
my @jsonFiles = @ARGV;
my %jsonHash=GSI::bamqc::load_json(@jsonFiles); #### contains decoded json data. Hash keys are filename, values are json hashes



my %runList;  ### contains a list of runs and lanes within the run
map{
	my $run =$jsonHash{$_}{"run name"};
	my $lane=$jsonHash{$_}{"lane"};
	$runList{$run}{$lane}++;
}keys %jsonHash;


GSI::RunReport::plot_data(\%jsonHash,$scriptPath) if($param{plotData});

#used for tsv file

my $html;
for my $run (sort keys %runList)
{
	#header with name of run
	$html.="<h1><a name=\"$run\">$run</a></h1>\n";

	### for this run, get a list of reports
	### do NOT want to send the run list to each of the functions
	my @sorted_lanes=sort{$jsonHash{$a}{lane}<=>$jsonHash{$b}{lane}} grep{$jsonHash{$_}{"run name"} eq $run} keys %jsonHash;

	$html.=GSI::RunReport::write_tsv(\%param,\%jsonHash,\@sorted_lanes,$run);
	$html.=GSI::RunReport::data_table(\%param,\%jsonHash,\@sorted_lanes,$run);  ### parameters, jsonHash, jsonfiles - in order, id
	$html.=GSI::RunReport::coverage_table(\%param,\%jsonHash,\@sorted_lanes,$run) 	if($param{showCoverageTable}	);
	$html.=GSI::RunReport::graph_table(\%param,\%jsonHash,\@sorted_lanes,$run) 		if($param{showGraphTable}	&& $param{plotData} );   ### set this as a parameter than can be turned off, on by default
	$html.=GSI::RunReport::lane_info(\%param,\%jsonHash,\@sorted_lanes,$run) 		if($param{showLaneInfo}		);

}

my $page = GSI::RunReport::assemble_run_report($html);
print $page;


exit;

sub usage{
        print "\nUsage is jsonToGenericRunReport.pl [options] path/to/*.json\n";
        print "Options are as follows:\n";
				print "\t -c show bases covered data.  Default is no bases covered.\n";
        print "\t-r show read length histogram.  Default is no histogram.\n";
        print "\t-p print all images.  Default is to only show thumbnails (with links).\n";
        print "\t-g plot data.  Default it to not plot the data\n";
				print "\t-H show hard clip stats and graph.  Default is off.\n";
				print "\t-n no collapse estimate!\n";
        print "\t-h displays this usage message.\n";

        die "\n@_\n\n";
}
