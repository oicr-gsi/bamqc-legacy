#!/usr/bin/perl
#Original Script : Rob Denroche
#Modifications : Lawrence Heisler << <lheisler.oicr.on.ca> >>
#Last modified : 2015-01-07
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


use Cwd 'abs_path';
use Getopt::Std;
use File::Basename;
use vars qw/ %opt /;

use lib dirname (__FILE__);
use bamqc;


use Data::Dumper;

my $scriptPath = dirname(abs_path($0));

my %table_headers=(
	data=>{
		lane		=>	"Lane",
		barcode		=>	"Barcode",
		groupid		=>	"Group ID",
		ext_name	=>	"External Name",
		library		=>	"Library",
		insert_mean	=>	"Insert Mean (SD)",
		read_length	=>	"Read Length",
		raw_reads	=>	"Raw Reads",
		raw_yield	=>	"Raw Yield",
		mapped		=>	"Map %",
		error		=>	"Error %",
		softclip	=>	"Soft Clip %",
		hardclip	=>  "Hard Clip %",
		rpsp		=>	"Reads/SP",
		ontarget	=>	"% on Target",
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
	}
);

my %plot_names=(
	read_breakdown		=>"readPie.png",
	insert_distr		=>"insert.png",
	read_length			=>"readLength.png",
	qual_hist	 		=>"qualHist.png",
	qual_by_cycle 		=>"qualCycle.png",
	mismatch_by_cycle 	=> "misCycle.png",
	indel_by_cycle		=>	"indelCycle.png", 
	softclip_by_cycle	=>	"softCycle.png",
	hardclip_by_cycle	=>	"hardCycle.png"
);



my %param=(
	
	#### ordered data_table_columns, modify to limit to specific columns
	table_headers=>\%table_headers,
	table_columns=>{   #### will show all possible columns, removing those that need to be explicitly indicated
		data	=> [qw/lane barcode groupid ext_name library insert_mean read_length raw_reads raw_yield mapped error softclip hardclip rpsp ontarget yield coverage/],
		graph	=> [qw/lane barcode library read_breakdown insert_distr qual_hist qual_by_cycle mismatch_by_cycle indel_by_cycle softclip_by_cycle hardclip_by_cycle/],
	},
	plotnames=>\%plot_names,
	coverageXs => [qw(1 4 8 15 30 50 100 200)],
	showDataTable 		=> 1,
	showGraphTable 		=> 1,
	showCoverageTable 	=> 0,
	showLaneInfo 		=> 1,
	printAllImages 		=> 0,
	plotData 			=> 0,
);






#### process options
my $opt_string = "crpgHh";
getopts ($opt_string, \%opt) or usage("Incorrect arguments.");
usage("Help requested.") if (exists $opt{h});

$param{showCoverageTable} = 1 if (exists $opt{c});
$param{printAllImages} = 1 if (exists $opt{p});
$param{plotData} = 1 if(exists $opt{g});

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
my %jsonHash = load_json(@jsonFiles); #### contains decoded json data. Hash keys are filename, values are json hashes



### map to get runList
#print STDERR Dumper(keys %{$jsonHash{"runlevel.jsonReport/SWID_1252183_CPCG_0329_Pr_P_PE_656_WG_141007_SN203_0253_AC5NM1ACXX_NoIndex_L008_R1_001.annotated.bam.BamQC.json"}});exit;

my %runList;  ### contains a list of runs and lanes within the run
map{
	my $run =$jsonHash{$_}{"run name"};
	my $lane=$jsonHash{$_}{"lane"};
	$runList{$run}{$lane}++;
}keys %jsonHash;
  
my $date = `date`;
chomp $date;
unless (-e "sorttable.js")
{
	`ln -s /u/lheisler/git/spb-analysis-tools/bamqc/sorttable.js`;
}


my $html;
$html.="<html>\n<head>\n";
$html.="<script src=\"./sorttable.js\"></script>\n";
$html.="<style type=\"text/css\">\n.na { color: #ccc; }\nth, td {\n  padding: 3px !important;\n}\ntable\n{\nborder-collapse:collapse;\n}\n/* Sortable tables */\ntable.sortable thead {\n\tbackground-color:#eee;\n\tcolor:#000000;\n\tfont-weight: bold;\n\tcursor: default;\n}\n</style>\n";
$html.="</head>\n<body>\n";
$html.="<p>Generic run report generated on $date.</p>\n";
for my $run (sort keys %runList)
{
	### for this run, get a list of reports
	### do NOT want to send the run list to each of the functions
	
	
	my $plotted=plot_data(\%jsonHash,$scriptPath) if($param{plotData});
	$html.="<h1><a name=\"$run\">$run</a></h1>\n";
	
	my @json=sort{$jsonHash{$a}{lane}<=>$jsonHash{$b}{lane}} grep{$jsonHash{$_}{"run name"} eq $run} keys %jsonHash;
	
	$html.=data_table(\%param,\%jsonHash,\@json,$run);  ### parameters, jsonHash, jsonfiles - in order, id
	$html.=coverage_table(\%param,\%jsonHash,\@json,$run) 	if($param{showCoverageTable}	);
	$html.=graph_table(\%param,\%jsonHash,\@json,$run) 		if($param{showGraphTable}	);   ### set this as a parameter than can be turned off, on by default
	$html.=lane_info(\%param,\%jsonHash,\@json,$run) 		if($param{showLaneInfo}		);
}
$html.="</body>\n</html>\n";	
print $html;


exit;


sub usage{
        print "\nUsage is jsonToGenericRunReport.pl [options] path/to/*.json\n";
        print "Options are as follows:\n";
		print "\t -c show bases covered data.  Default is no bases covered.\n";
        print "\t-r show read length histogram.  Default is no histogram.\n";
        print "\t-p print all images.  Default is to only show thumbnails (with links).\n";
        print "\t-g plot data.  Default it to not plot the data\n";
		print "\t-H show hard clip stats and graph.  Default is off.\n";
        print "\t-h displays this usage message.\n";

        die "\n@_\n\n";
}




