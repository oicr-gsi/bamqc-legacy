package GSI::report;

#Copyright (C) 2017 The Ontario Institute for Cancer Research
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

bamqc.pm - Perl Library for OICR bamqc and reporting

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

Subroutines in support of bamqc analysis, and construction of html reports

=over 2

=item * Usage: use bamqc 

=back

=head1 AUTHOR

Lawrence Heisler << <lheisler.oicr.on.ca> >>
Last modified : 2014-12-08
The Ontario Institute for Cancer Research
Incorporating code from the original samStats.pl script by Rob Denroche


=cut

use strict;
use warnings;
use Exporter;
use JSON::PP;
use File::Basename;


### for debug

#use Data::Dumper;
#(open my $TTY,"/dev/tty") || die "unable to open keyboard input";

my @month = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );



use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION	=	1.00;
@ISA		=	qw(Exporter);
@EXPORT		=	qw( load_param plot_data data_table coverage_table graph_table merge_table 
					lane_info );

###### reporting scripts, exported
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
	read_length			=> "readLength.png",
	qual_hist	 		=> "qualHist.png",
	qual_by_cycle 		=> "qualCycle.png",
	mismatch_by_cycle 	=> "misCycle.png",
	indel_by_cycle		=> "indelCycle.png", 
	softclip_by_cycle	=> "softCycle.png",
	hardclip_by_cycle	=> "hardCycle.png"
);


sub load_param{
	
	
	my %param=(
		#### ordered data_table_columns, modify to limit to specific columns
		table_headers=>\%table_headers,
		table_columns=>{   #### will show all possible columns, removing those that need to be explicitly indicated
			data	=> [qw/lane barcode groupid ext_name library insert_mean read_length raw_reads raw_yield mapped error softclip hardclip rpsp ontarget yield coverage/],
			graph	=> [qw/lane barcode library read_breakdown insert_distr qual_hist qual_by_cycle mismatch_by_cycle indel_by_cycle softclip_by_cycle hardclip_by_cycle/],
			merge	=> [qw/sample lanes last_mod avg_read_length error_pc soft_clip_pc reads_on_target_pc reads_on_target bases_on_target coverage coverage8x coverage90 coverageGraph/],
		},
		plotnames=>\%plot_names,
	);
	return %param;  ### return/capture as a hash
};


sub plot_data{
	my ($j,$scriptPath, $jsonDirs)=@_;
	for my $rpt (keys %$j){
		#my $title=exists $$j{$rpt}{"barcode"}? 
		#		$$j{$rpt}{"run name"} . " Lane: " . $$j{$rpt}{"lane"} . " Barcode: " . $$j{$rpt}{"barcode"} . "\\n" . $$j{$rpt}{"library"}
		#		$$j{$rpt}{"run name"} . " Lane: " . $$j{$rpt}{"lane"} . "\\n" . $$j{$rpt}{"library"};
		warn "graphing $rpt\n";
		my $rv=`$scriptPath/jsonToGraphs.pl $jsonDirs->{$rpt}/$rpt`;
		return $rpt;
	}
}
sub data_table{
	my($p,$j,$json,$run)=@_;
	my $html="<table border=\"1\" class=\"sortable\">\n";
	my %stats;

	####   header
	$html.="<thead>\n";
	$html.="<tr>\n";
	$html.=data_row($p,"header")."\n";
	$html.="</thead>\n";
	
	#### body, 1 row per lane
	$html.="<tbody>\n";
	for my $rpt(@$json){
		my ($data_row,$row_stats)=data_row($p,"body",$$j{$rpt});
		$html.="$data_row\n";
		map{$stats{$_}+=$$row_stats{$_}} qw/RawYield Reads EstYield/;
		map{$stats{qualCuts}{$_}++} keys %{$$row_stats{qualCut}};
	}

	$html.="</tbody>\n";

	####  footer	
	$html.="<tfoot>\n<tr>\n";
	$html.=data_row($p,"footer",\%stats);
	$html.="</tr></tfoot>\n";
	$html.="</table>\n";
	
	my $qualCut = join(", ",sort{$b<=>$a} keys %{$stats{qualCuts}});
	$html.="<p>* Estimates exclude unmapped, off target, non-primary or MAPQ < $qualCut reads and use reads per start point to approximate loss due to collapse.</p>\n";
	
	
	return $html;							
}
sub coverage_table{
		
	my($p,$j,$json,$run)=@_;
	
	my @covX=@{$$p{coverageXs}};
	my $levels=scalar @covX;

	my $html="<table border=\"1\" class=\"sortable\">\n<thead>\n<tr>\n";
	$html.="<th colspan=\"5\"></th>";
	$html.="<th colspan=\"" . $levels . "\">Non-collapsed % bases covered</th>";
	$html.="<th colspan=\"" . $levels . "\">Collapsed % bases covered</th>";

	$html.="</tr>\n";
	$html.="<tr>";
	$html.="<th>Lane</th>";
	$html.="<th>Barcode</th>";
	$html.="<th>Library</th>";
	$html.="<th>Target Size (bp)</th>";
	$html.="<th># Targets</th>";

	my $covXheadings=join("",map{"<th>${_}x</th>"} @covX);
	$html.=$covXheadings.$covXheadings;
	$html.="\n</tr>\n</thead>\n<tbody>\n";


	
	for my $rpt (@$json){
		my $lane=$$j{$rpt}{"lane"};	
		
		my $targetSize = $$j{$rpt}{"target size"};
		$targetSize =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
		my $numberOfTargets = $$j{$rpt}{"number of targets"};
		$numberOfTargets =~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;

		my $linkName=exists $$j{$rpt}{barcode} ? "${run}_${lane}_$$j{$rpt}{barcode}" : "${run}_${lane}";
		$html.="<tr>\n";
		$html.="<td><a href=\"#$linkName\">$$j{$rpt}{lane}</a></td>";
			
		$html.=exists $$j{$rpt}{barcode} ? "<td><a href=\"#$linkName\">$$j{$rpt}{barcode}</a></td>" : "<td>none</td>";
		$html.="<td><a href=\"#$linkName\">$$j{$rpt}{library}</a></td>";
		$html.="<td>$targetSize</td>";
		$html.="<td>$numberOfTargets</td>";
		$targetSize = $$j{$rpt}{"target size"};		# so it isn't formatted with commas
		for my $collapsed ("non collapsed", "collapsed"){
			for my $i (@covX){
				my $basesCovered = exists $$j{$rpt}{"$collapsed bases covered"}{$i} ? $$j{$rpt}{"$collapsed bases covered"}{$i} : 0;
				$basesCovered = $basesCovered / $targetSize;
				$basesCovered = sprintf "%.2f%%", $basesCovered * 100;
				$html.="<td>$basesCovered</td>";
			}
		}
	}
	$html.="</tr>\n";
	$html.="</table>\n";
	$html.="<p></p>\n";
	
	return $html;
}
sub graph_table{
	my($p,$j,$json,$run)=@_;
	
	my @cols=@{$$p{table_columns}{graph}};

	# print image thumbnail table
	my $html="<table border=\"1\" class=\"sortable\">\n<thead>\n<tr>\n";
	
	my @headings=@{$$p{table_headers}{graph}}{@cols};
	
	my @html_headings=map{"<th>$_</th>"} @headings;
	$html.="<thead><tr>".join("",@html_headings)."<\tr></thead>\n";
	
	$html.="<tbody>\n";
	for my $rpt (@$json){
		$html.="<tr>\n";
		my $lane=$$j{$rpt}{"lane"};	
		my $linkName=exists $$j{$rpt}{barcode} ? "${run}_${lane}_$$j{$rpt}{barcode}" : "${run}_${lane}";
		# read_breakdown insert_distr qual_hist qual_by_cycle mismatch_by_cycle indel_by_cycle softclip_by_cycle
		my %td;
			
		$td{run_name}="<td>".$$j{$rpt}{"run name"}."</td>";
		$td{lane}="<td>$lane</td>";
		$td{barcode}=exists $$j{$rpt}{barcode} ? "<td>$$j{$rpt}{barcode}</td>" : "<td></td>";
		$td{library}="<td>$$j{$rpt}{library}</td>";

					
		if($$p{linked_columns}){
			my $tag_open="<a href=\"#$linkName\">";
			my $tag_close="</a>";
			for my $col(@{$$p{linked_columns}}){
				if($td{$col}){
					$td{$col}=~s/<td>/<td>$tag_open/ ;
					$td{$col}=~s/<\/td>/$tag_close<\/td>/;
				}
			}
		}				
			
		### the remaining columns are all plots
		for my $col(@cols){
			next if($td{$col});  ## skip if column is already set (ie. a non-plot columns)
			my $plot=$$p{plotnames}{$col} || "";
			$td{$col} ="<td><a href=\"${rpt}.graphs/$plot\">";
			$td{$col}.="<img src=\"${rpt}.graphs/$plot\" width=\"100\" height=\"100\"/>";
			$td{$col}.="</a></td>";
		}

		$html.=join("",@td{@cols});
		$html.="</tr>\n";
	}
$html.="</tbody>\n</table>\n";
}
sub merge_table{
	my($p,$j,$json,$run,$lanes)=@_;
	
	my @cols=@{$$p{table_columns}{merge}};
	
	if($$p{seq_type} ne "exome"){
		## remove the exome columns
		@cols = splice(@cols,0,-3);  ### remove last three columns
	}

	
	
	# print image thumbnail table
	my $html="<table border=\"1\" class=\"sortable\">\n<thead>\n<tr>\n";
	my @headings=@{$$p{table_headers}{merge}}{@cols};
	my @html_headings=map{"<th>$_</th>"} @headings;
	$html.="<thead><tr>".join("",@html_headings)."<\tr></thead>\n";
	$html.="<tbody>\n";
 	for my $rpt (@$json){
 		$html.="<tr>\n";
 		my %td = map{ ($_,"<td></td>") } @cols;  ### initialize as empty cells
 
 		$td{sample}="<td>".$$j{$rpt}{sample}."</td>";
 		$td{lanes}="<td>".$lanes."</td>";
 		
 		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($$j{$rpt}{"last modified"});
		$year += 1900;
		$td{last_mod}="<td>$month[$mon] $mday, $year</td>";

		### exome will also include coverage8x coverage90 coverageGraph
		my $paired= $$j{$rpt}{"number of ends"} eq "paired end" ? 1 : 0;

		my $reads_on_target=$$j{$rpt}{"reads on target"};
		$td{reads_on_target}="<td>".format_number($reads_on_target)."</td>";
		
		if($paired){
			my ($R1,$R2)=($$j{$rpt}{"read 1 average length"},$$j{$rpt}{"read 2 average length"});
			
			$td{bases_on_target}="<td>". format_number($reads_on_target * (($R1+$R2)/2))."</td>";
			($R1,$R2)= map{ m/^\d+\.\d+$/ ? sprintf "%.2f",$_ : $_ } ($R1,$R2);
			$td{avg_read_length}="<td>$R1,$R2</td>";

		}else{
			my $R1=$$j{$rpt}{"read 1 average length"};
			$td{bases_on_target}="<td>". format_number($reads_on_target * $R1)."</td>";
			($R1)= map{ m/^\d+\.\d+$/ ? sprintf "%.2f",$_ : $_ } ($R1);
			$td{avg_read_length}="<td>$R1</td>";
		}
		my $percentOnTarget=$reads_on_target/$$j{$rpt}{"mapped reads"}*100;
		$td{reads_on_target_pc}="<td>". (sprintf "%.2f",$percentOnTarget )   ."</td>";

 	   my ($errorRate,$softClipRate,$hardClipRate) = ("0%","0%","0%");
       if ($$j{$rpt}{"aligned bases"} > 0){
			$errorRate    = sprintf "%.2f%%", (($$j{$rpt}{"mismatch bases"} + $$j{$rpt}{"inserted bases"} + $$j{$rpt}{"deleted bases"}) / $$j{$rpt}{"aligned bases"}) * 100;
            $softClipRate = sprintf "%.2f%%", $$j{$rpt}{"soft clip bases"} / ($$j{$rpt}{"aligned bases"} + $$j{$rpt}{"soft clip bases"} + $$j{$rpt}{"hard clip bases"}) * 100;
            $hardClipRate = sprintf "%.2f%%", $$j{$rpt}{"hard clip bases"} / ($$j{$rpt}{"aligned bases"} + $$j{$rpt}{"soft clip bases"} + $$j{$rpt}{"hard clip bases"}) * 100;
		}
		$td{error_pc}="<td>$errorRate</td>";
		$td{soft_clip_pc}="<td>$softClipRate</td>";
		$td{hard_clip_pc}="<td>$hardClipRate</td>";
		
		$td{RPSP}     = "<td>".(sprintf "%.2f", $$j{$rpt}{"reads per start point"})."</td>";
		my $estimatedYield = int($$j{$rpt}{"aligned bases"} * ($percentOnTarget / 100));
		my $estimatedCoverage = $estimatedYield / $$j{$rpt}{"target size"};
		$td{coverage} ="<td>".(sprintf "%.2f", $estimatedCoverage)."</td>";

		my $groupid_title=$$j{$rpt}{"group id description"} ? "title='" . $$j{$rpt}{"group id description"} : "";
		
		$td{groupid}=exists $$j{$rpt}{"group id"} ? "<td $groupid_title >".$$j{$rpt}{"group id"}."</td>"
												 : "<td class='na'>na</td>";
		$td{ext_name}=exists $$j{$rpt}{"external name"} ? "<td>".$$j{$rpt}{"external name"}."</td>"
												 : "<td class='na'>na</td>";

		### check if the column coverage8x is requested.  ### exome data!
		if($$p{seq_type} eq "exome"){
			if(exists $$j{$rpt}{"collapsed bases covered"}{8}){
				$td{coverage8x}= "<td>".$$j{$rpt}{"collapsed bases covered"}{8}."%</td>";
	 			### get al levels >=90, sort hi to low				     
				my @cov90=grep{$$j{$rpt}{"collapsed bases covered"}{$_}>90}  sort {$b<=>$a}  keys %{$${$rpt}{"collapsed bases covered"}};
				my $cov90=shift @cov90;  ### get the hightest
				$cov90=$cov90>50 ? $cov90 : ">50";
				$td{coverage90}="<td>${cov90}x</td>";
				my $img=$$j{$rpt}{"exome histogram image"};
				$td{coverageGraph}="<td><a href=\"../$img\"><img src=\"$img\" width=\"27\" height=\"20\">";
			}else{
				$td{coverage8x}="<td colspan=\"3\">Bed coverage not found</td>";
			}
		}	
 		$html.=join("",@td{@cols});
 		$html.="</tr>\n";
 	}
	$html.="</tbody>\n</table>\n";
	
	return $html;
}
sub lane_info{
	my($p,$j,$json,$run)=@_;
	my @html_lines;
	for my $rpt (@$json){
		my $lane=$$j{$rpt}{"lane"};	
		my $barcode_line=exists $$j{$rpt}{"barcode"} ? "<h2><a name=\"${run}_${lane}_$$j{$rpt}{\"barcode\"}\">$run Lane: $lane Barcode: $$j{$rpt}{\"barcode\"}</a></h2>"
														 : "<h2><a name=\"${run}_${lane}\">$run Lane: $lane</a></h2>";
		push(@html_lines,$barcode_line);
		push(@html_lines,"<p><a href=\"#$run\">Back to $run.</a></p>");
		push(@html_lines,"<ul>");
		push(@html_lines,"<li>Library: $$j{$rpt}{\"library\"}</li>");
		push(@html_lines,"<li>Target file: $$j{$rpt}{\"target file\"}</li>");
		my $targetSize = format_number($$j{$rpt}{"target size"});
		my $numberOfTargets = format_number($$j{$rpt}{"number of targets"});
		push(@html_lines,"<li>Target size: $targetSize bp</li>");
		push(@html_lines,"<li>Number of targets: $numberOfTargets</li>");
		
		push(@html_lines,"<li>Workflow Name: $$j{$rpt}{'workflow name'}</li>") 
			if($$j{$rpt}{"workflow name"}); 
		push(@html_lines,"<li>Workflow Version: $$j{$rpt}{'workflow version'}</li>") 
			if($$j{$rpt}{"workflow version"}) ;
		push(@html_lines,"<li>Workflow Run Creation Timestamp: $$j{$rpt}{'workflow run creation timestamp'}</li>")
			if($$j{$rpt}{"workflow run creation timestamp"}); 
        push(@html_lines,"<li>Bam File Creation Timestamp: $$j{$rpt}{'bam file creation timestamp'}</li>")
          	if($$j{$rpt}{"bam file creation timestamp"});
		push(@html_lines,"</ul>");
			
		if($$p{printAllImages}){
			for my $col(@{$$p{table_columns}{graph}}){  ### the columns of the graph table, show these plots if they exists
				if(my $plot=$$p{plotnames}{$col}){
					push(@html_lines,"<img src=\"${rpt}.graphs/$plot\"/>");
				}
			}
		}
	}
	
	my $html=join("\n",@html_lines);
	return $html;
	
}

#####  reporting scripts, internal
sub format_number{
	my ($n)=@_;
	$n=~s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g;
	return $n;
}
sub data_row{
	my($p,$section,$j)=@_;  ### parameter hash, jsonhash or header indication
	
	my $html;my %stats;
	my @cols=@{$$p{table_columns}{data}};
	
	if($section eq "header"){
		my @headings=@{$$p{table_headers}{data}}{@cols};
		my @html_headings=map{"<th>$_</th>"} @headings;
		$html=join("",@html_headings);
		return $html;

	}elsif($section eq "footer"){
		map{ $$j{$_}=~ s/(^[-+]?\d+?(?=(?>(?:\d{3})+)(?!\d))|\G\d{3}(?=\d))/$1,/g; } qw/RawYield Reads EstYield/;
		
		my %hash;
		map{$hash{$_}=""} @cols;
		
		@hash{($cols[0],"raw_reads","raw_yield","yield")}=("Total",$$j{Reads},$$j{RawYield},$$j{EstYield});
		my @footings=@hash{@cols};
		my @html_footings=map{"<th>$_</th>"} @footings;
		$html=join("",@html_footings);
		return $html;
	
	}elsif($section eq "body"){
	
	
		### the %td hash will hold the formatted contents of each possible cell in the row
		### the %td hash key = column identifier
		### the final list of cells to include in the row is @cols, extracted from %p
		
		
		my %td;   ### hash to hold table cells <td>..</td>
		### RUN NAME
		$td{run_name}="<td>".$$j{"run name"}."</td>";
		
		## LANE
		#$td{lane}="<td><a href=\"#$linkName\">$$j{lane}</a></td>";
		$td{lane}="<td>$$j{lane}</td>";
		## BARCODE OR NONE
		#$td{barcode}= exists $$j{barcode} ? "<td><a href=\"#$linkName\">$$j{barcode}</a></td>"
		#							: "<td>none</td>";
		$td{barcode}= exists $$j{barcode} ? "<td>$$j{barcode}</td>"
									       : "<td>none</td>";
		
		
		## GROUPID if it exists
		$td{groupid}= exists $$j{"group id"} ? 
               						exists $$j{"group id description"} ?  "<td title='".$$j{"group id description"}."'>$$j{'group id'}</td>"
             															: "<td>$$j{'group id'}</td>"
 									: "<td class='na'>na</td>";
 		## EXTERNAL NAME if it exists
 		$td{ext_name}=exists $$j{"external name"} ? "<td>$$j{'external name'}</td>"
											: "<td class='na'>na</td>";
		## LIBRARY
		#$td{library}="<td><a href=\"#$linkName\">$$j{library}</a></td>";
		$td{library}="<td>$$j{library}</td>";
		
		## rawReads = mapped reads + unmapped reads + reads with failed quality
		my $rawReads = $$j{"mapped reads"} + $$j{"unmapped reads"} + $$j{"qual fail reads"};
		$td{raw_reads}="<td>".format_number($rawReads)."</td>";
		
		## rawYield = rawReads X average read length = bases
		my $rawYield = int($rawReads * $$j{"average read length"});		# cast to int because average length is only from mapped reads (and may result in ugly decimal)
		$td{raw_yield}="<td>".format_number($rawYield)."</td>";
		
		## readlength, X,Y if paired X if not
		my $readLength=$$j{"number of ends"} eq "paired end" ? $$j{"read 1 average length"} . ", " . $$j{"read 2 average length"}
							   								: sprintf "%.2f", $$j{"read ? average length"};
		$td{read_length}="<td>$readLength</td>";

		### map rate = # mapped/ # rawReads
		my $mapRate = $rawReads>0 ? sprintf "%.2f%%",($$j{"mapped reads"} / $rawReads) * 100
								: "0%";
		$td{mapped}="<td>$mapRate</td>";

		### errorRate = mismatches+inserted+deleted / aligned
		my $errorRate= $$j{"aligned bases"}>0 ? sprintf "%.2f%%", (($$j{"mismatch bases"} + $$j{"inserted bases"} + $$j{"deleted bases"}) / $$j{"aligned bases"}) * 100
											  : "0%";
		$td{error}="<td>$errorRate</td>";
		
		### softClipRate = soft clip bases / aligned+softclip+hardclip		
		my $softClipRate= $$j{"aligned bases"}>0 ? sprintf "%.2f%%", $$j{"soft clip bases"} / ($$j{"aligned bases"} + $$j{"soft clip bases"} + $$j{"hard clip bases"}) * 100
											  : "0%";
		$td{softclip}="<td>$softClipRate</td>";
		
		### hardClipRate = hardclip / aligned+softclip+hardclip
		my $hardClipRate= $$j{"aligned bases"}>0 ? sprintf "%.2f%%", $$j{"hard clip bases"} / ($$j{"aligned bases"} + $$j{"soft clip bases"} + $$j{"hard clip bases"}) * 100
											  : "0%";
		$td{hardclip}="<td>$hardClipRate</td>";
											  
		### reads per startpoint									  
		my $readsPerStartPoint = sprintf "%.2f", $$j{"reads per start point"};
		$td{rpsp}="<td>$readsPerStartPoint</td>";

		### onTargetRate = onTarget/mapped (could use onTarget/raw)
		my $onTargetRate = $rawReads > 0 ? ($$j{"reads on target"} / $$j{"mapped reads"}) * 100  # $rawReads) * 100;   # could argue using this either way
										: 0;
		### estimatedYield = alignedbases * onTarget rate / reads per start point
		my $estimatedYield = int(($$j{"aligned bases"} * ($onTargetRate / 100)) / $$j{"reads per start point"});
		$td{yield}="<td>".format_number($estimatedYield)."</td>";
		
		$td{ontarget}="<td>".(sprintf "%.2f%%",$onTargetRate)."</td>";
		
		my $estimatedCoverage = $estimatedYield / $$j{"target size"};
		$td{coverage}="<td>".(sprintf "%2f",${estimatedCoverage})."x</td>\n";

		$td{insert_mean} = $$j{"number of ends"} eq "paired end" ?
								"<td>".sprintf ("%.2f", $$j{"insert mean"}) . " (" .sprintf("%.2f", $$j{"insert stdev"}) . ")</td>"
								: "<td></td>";

		
		#### add hyperlinks to the required cells
		my $linkName= exists $$j{barcode} ? $$j{"run name"} ."_$$j{lane}_$$j{barcode}"
											: $$j{"run name"}."_$$j{lane}";
		
		if($$p{linked_columns}){
			my $tag_open="<a href=\"#$linkName\">";
			my $tag_close="</a>";
			for my $col(@{$$p{linked_columns}}){
				$td{$col}=~s/<td>/<td>$tag_open/;
				$td{$col}=~s/<\/td>/$tag_close<\/td>/;
			}
		}
		
		### these are cummulative totals, need to be added and returne....as part of parameters?  or caputred
		%stats=(RawYield=>$rawYield,Reads=>$rawReads,EstYield=>$estimatedYield);
		
		$stats{qualCut}{$$j{"qual cut"}} = 1;

		$html="<tr>\n".join("\n",@td{@cols})."</tr>";
		
		return ($html,\%stats);
	}
	
}

1;
