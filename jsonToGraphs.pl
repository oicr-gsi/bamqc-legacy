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
use Data::Dumper;


use JSON::PP;

#### identify location of this script, which should also contain the library jsonToGraphs.pm
#### this can be changes should the library already exist in @INC
use File::Basename;
use lib dirname (__FILE__);
use jsonToGraphs;

### currently on argument, the path to the json file
# add usage and ability to select subset of graphs
# add ability to change location of output directory

my $jsonFile = $ARGV[0];		
my $plot_dir="${jsonFile}.graphs";

my %jsonHash;
open (my $JSON,"<",$jsonFile) or die "Couldn't open $jsonFile.\n";
if (my $line = <$JSON>){
	%jsonHash = %{ decode_json($line) };
	mkdir $plot_dir unless(-d $plot_dir);
}else{
	warn "No data found in $jsonFile!\n";
}
close $JSON;

### param will contain information to be sent to all functions
my %param=(
	# draw graphs
	colours=>{
		good 	=> "forestgreen",
		bad  	=> "firebrick",
		mid  	=> "goldenrod",
		other 	=>  [qw(darkslategray3 darkorchid4 green2 red mediumblue orange grey yellow pink forestgreen goldenrod firebrick)]
	},
	qualcut=>{
		low     => 20,
		high	=> 30,
	},
	basecoverage=>{
		low		=> 80,
		high	=> 90,
		displaymax	=> 200
	},
	reads=>["read 1","read 2","read ?"],
	insertMax => 650, # should really pass this in
	qualLineMax =>50, # seems reasonable for phred scale
	insertStep => 50,
);

### plot title is take from the jsonHash, will be modified if there is a barcode indicated
my $title;
if($jsonHash{"library"} eq "merged"){
	$title=$jsonHash{"sample"} ."\\nmerged";
}else{
	$title = exists $jsonHash{"barcode"} ?
				$jsonHash{"run name"} . " Lane: " . $jsonHash{"lane"} . " Barcode: " . $jsonHash{"barcode"} . "\\n" . $jsonHash{"library"}
			: 	$jsonHash{"run name"} . " Lane: " . $jsonHash{"lane"} . "\\n" . $jsonHash{"library"};
}


## dispatch table, with named references to each subroutine
my %subs=(
	readmap_piechart 			=> \&readmap_piechart,
	quality_histogram 			=> \&quality_histogram,
	collapsed_base_coverage 	=> \&collapsed_base_coverage,
	noncollapsed_base_coverage 	=> \&noncollapsed_base_coverage,
	readlength_histogram 		=> \&readlength_histogram,
	insert_graph 				=> \&insert_graph,
	quality_by_cycle 			=> \&quality_by_cycle,
	mismatch_by_cycle 			=> \&mismatch_by_cycle,
	indel_by_cycle 				=> \&indel_by_cycle,
	softclip_by_cycle 			=> \&softclip_by_cycle,
	hardclip_by_cycle 			=> \&hardclip_by_cycle,
	coverage_by_depth			=> \&coverage_by_depth,
);

### list of subroutines to call
### this should be modifiable by including list on the command line
### default is to run all of these
my @plots=qw/readmap_piechart quality_histogram collapsed_base_coverage noncollapsed_base_coverage
				readlength_histogram insert_graph quality_by_cycle mismatch_by_cycle
				indel_by_cycle softclip_by_cycle hardclip_by_cycle/;
				
#my @plots=qw/coverage_by_depth/;



for my $plot(@plots){
	### each sub is called with 4 arguments, references to the jsonHash, param Hash, the directory to place plots and the plot title
	my $rv=&{$subs{$plot}}(\%jsonHash,\%param,$plot_dir,$title);
	
	### indicate failed plots
	print STDERR "$jsonFile : no plot generated for $plot\n" unless($rv);
}









