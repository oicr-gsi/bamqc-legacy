package bamqc;

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
@EXPORT		=	qw(read_bed assess_start_point assess_flag cigar_stats  md_stats
					onTarget addRunningBaseCoverage runningBaseCoverage HistStats insertMapping
					load_param plot_data data_table coverage_table graph_table merge_table 
					lane_info load_json toPhred generate_jsonHash);

=pod

=over 2

=item assess_flag()

=back

=over 4

=item Arguments: 

$flag : the sam record flag

$stats : a reference to a hash, which will be modified

$qual : the alignment quality

$param : the parameter hash

=item Returns: 

1 if the read is mapped, 0 if not

=back

Assess the sam flag and stores information in the reference hash
   
=cut
sub assess_flag{
	my($flag,$stats,$qual,$qcut)=@_;
	my $mapped=0;
	### processing ALL reads, not just sampled
	$$stats{"total reads"}++;
	
	
	if ($flag & 256){								# not primary alignment
		$$stats{"non primary reads"}++;
	}elsif ($flag & 4){								# read unmapped
		$$stats{"unmapped reads"}++;
	}elsif ($qual < $qcut	){					# below qual cut
		$$stats{"qual fail reads"}++;
	}else{
		$mapped=1;								# read is mapped (I hope)
		$$stats{"mapped reads"}++;				# increment mapped reads
		if ($flag & 1){						# read is paired
			$$stats{"paired reads"}++;
			if ($flag & 8){
				$$stats{"mate unmaped reads"}++	# mate is unmapped
			}
		}
		if ($flag & 2){
			$$stats{"properly paired reads"}++	# read is properly paired
		}
	}
	return $mapped;
}

=pod

=head1 Subroutines

=over 2

=item assess_start_point()

=back

=over 4

=item Arguments: 

$chrom : the chromosome of the current record

$s1	 : the position of the current record = the leftmost start point

$s2    : the position of the paired record = the rightmost start point

$sphash	 : the current StartPoint hash, holding the current position and the running ReadsPerStartPoint value

=item Returns: 

A revised StartPoint hash, with updated stats

=back

This subroutine keeps a running calculation of ReadsPerStartPoint, updating this value with each read.  It is dependent on the stream alignment data being properly sorted
   
=cut
sub assess_start_point{
	my ($chrom,$s1,$s2,$sphash)=@_;   ### receives the left/right start points and a reference to the hash
	   
	my $leftstart="$chrom\t$s1";
	my $pairStart="$chrom\t$s1\t$s2";
	
	#### Recalculate ReadsPerStartPoint when the leftstart changes from the previous record
	if ($leftstart ne $$sphash{current}){         
		### process all stored pairStarts
		for my $sp (keys %{$$sphash{pairStarts}}){
			$$sphash{count}++;   #### count each distinct pairedstart as a distinct start point
			
			### this is a running total adjusting positively or negativeley depending on the number of time the pair is seen
			$$sphash{RPSP} += ($$sphash{pairStarts}{$sp} - $$sphash{RPSP} ) / $$sphash{count};  
			
			### delete this pairStart after processing
			delete $$sphash{pairStarts}{$sp};
		}
		## reset the current position
		$$sphash{current} = $leftstart;
	}
	
	### store the pairStart 
	$$sphash{pairStarts}{$pairStart}++;     #### collect the paired starts
	return $sphash;  ### return the modified hash
}




=pod

=over 2

=item initializeCycleHash()

=back

=over 4

=item Arguments: 

=item Returns: 

=back

Description
   
=cut

sub initializeCycleHash{
	my $hash = $_[0];
	my $max = $_[1];

	for (my $i = 1; $i < $max; $i++)
	{
		$hash->{$i} = 0;
	}

}

=pod

=over 2

=item printHist()

=back

=over 4

=item Arguments: 

=item Returns: 

=back

Description
   
=cut

sub printHist{
	my %hist = %{ $_[0] };

	for my $i (sort {$a <=> $b} keys %hist)
	{
#		warn " $hist{$i}";
	}
}

=pod

=over 2

=item toPhred()

=back

=over 4

=item Arguments: 

=item Returns: 

=back

Description
   
=cut

sub toPhred{
	my $char = $_[0];
	my $ascii = ord($char);
	my $offset = 33;
	return $ascii - $offset;
}
sub mean{
	my $val = $_[0];
	my $sum = 0;
	my $count = 0;
	for my $v (@{ $val })
	{
		$sum += $v;
		$count++;
	}
	if ($count > 0)
	{
		return $sum / $count;
	}
	else
	{
		return 0;
	}
}
sub stdev{
	my $val = $_[0];
	my $mean = mean($val);
	my $squareDiff = 0;
	my $count = 0;
	for my $v (@{ $val })
	{
		$squareDiff += (($v - $mean)*($v - $mean));
		$count++;
	}
	if ($count > 0)
	{
		return sqrt($squareDiff / $count);
	}
	else
	{
		return 0;
	}
}

=pod

=over 2

=item HistStats()

=back

=over 4

=item Arguments: 

%val : a hash, with keys = insert size, values = count for each insert size

=item Returns: 

mean and standard deviation of the insert sizes

=back

Description
   
=cut


sub HistStats{
	#my %val = %{ $_[0] };

	my ($val)=@_;
	my ($mean,$stdv)=(0,0);
	
	if($val){
		my($sum,$count)=(0,0);
		map{
			$sum += ($_ * $$val{$_});
			$count += $$val{$_};
		} keys %$val;
		$mean = $count>0 ? $sum/$count : 0;
	
		my ($squareDiff);
		map{
			$squareDiff += ((($_ - $mean)*($_ - $mean)) * $$val{$_});
		} keys %$val;
	
		$stdv = $count>0 ? sqrt($squareDiff/$count) : 0;
	}
	
	return ($mean,$stdv);
		
}




=pod

=over 2

=item findStart()

=back

=over 4

=item Arguments: 

=item Returns: 

=back

Description
   
=cut

sub findStart{
        my @cigarOp = @{ $_[0] };
        my $start = $_[1];

        if ($cigarOp[0] =~ /(.*)S/)     # if first cigar operation is a soft clip, adjust start point
        {
                $start -= $1;
        }

        return $start;
}

=pod

=over 2

=item findEnd()

=back

=over 4

=item Arguments: 

=item Returns: 

=back

Description
   
=cut

sub findEnd{
        my @cigarOp = @{ $_[0] };

        my $end = $_[1];

        if ($cigarOp[0] =~ /(.*)S/)     # if first cigar operation is a soft clip, adjust start point
        {
                $end -= $1;
        }

        foreach my $cig (@cigarOp)
        {
                if ($cig =~ /(.*)S/)
                {
                        $end += $1;
                }
                elsif ($cig =~ /(.*)M/)
                {
                        $end += $1;
                }
                elsif ($cig =~ /(.*)D/)
                {
                        $end += $1;
                }
        }

        return $end;
}

=pod

=over 2

=item read_bed()

=back

=over 4

=item Arguments: 

$file : the name of the bed file, containing the intervals

=item Returns: 

Hash structure containing the bed intervals

=back

Reads in a bed file, storing in the intervals in an array of hashes indicating, Start/Stop and size
   
=cut

sub read_bed{
	my ($file)=@_;
	my %bed;
	my $targetCount=0;
	open (my $BEDFILE,"<",$file) or return("ERROR : Couldn't open target file: $file.\n");
	while (<$BEDFILE>){
		chomp;
		next if(/^#/);
		my @f = split /\t/;
	
		$targetCount++;
		
		my $interval_size=$f[2]-$f[1];
		
		### interval_size needs to be larger than 0
		if($interval_size<1){
			return("ERROR : the bedfile $file contains an interval with an invalid size $_ \n");
		} 
		
		### also want to check that the sort order is correct
		
		
		push(@{$bed{intervals}{$f[0]}},{Start=>$f[1],Stop=>$f[2],Size=>$interval_size});
		$bed{targetSize} += $interval_size;
	    #$bed{Hist}{"$f[0]\t$f[1]\t$f[2]"} = 0;
	    
	    
	    
	    
	    
	}
	close $BEDFILE;
	
	$bed{numberOfTargets} = $targetCount;
	warn "Loaded " . $bed{numberOfTargets} . " targets.\n\n";
	
	
	
	return \%bed;
}

=pod

=over 2

=item cigar_stats()

=back

=over 4

=item Arguments: 

$chrom : the chromosome to which the current read maps

$start1 : the leftmost startpoint

$start2 : the rightmost startpoint for the paired end

$R : the read (R1,R2,R?)

$strand : the strand to which the read maps, from the sam flag

$cigar : the cigar string

$stats : a reference to the stats hash, which will be modified

$p : a reference to the parameter hash

=item Returns: 

A list with two values, $readLength and $mappedBases

=back

Processes the cigarString, by breaking into pieces and assessing how the read is mapped to the reference
   
=cut

sub cigar_stats{
	my($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)=@_;
	
	## the fragment positioning on the chromosome
	my $pairStart="$chrom\t$start1\t$start2";
	
	## split the cigar string, and process each piece 
	my @Op=procCigar($cigar,$strand);

    ## initialize these values
	my $readLength = 0;
	my $cycle = 1;   ### this will keep track of the current cycle
	my $posOffset = 0;
	my $mappedBases = 0;
	
	foreach my $c (@Op){
		### MATCH, capture the number of bases in $1
		if ($c =~ /(.*)M/){       
			$$stats{alignedCount} += $1;
			#for (my $i = 0; $i < $1; $i++){
			for my $i(1..$1){
				$$stats{byCycle}{aligned}{$R}{$cycle}++;
				$cycle++;
			}
			$readLength += $1;
			$mappedBases += $1;

		### HARDCLIPPED BASES
		}elsif ($c =~ /(.*)H/){
			$$stats{hardClipCount} += $1;
			for my $i(1..$1){
				$$stats{byCycle}{hardClip}{$R}{$cycle}++;
				$cycle++;
			}
			$readLength += $1;
		### SOFTCLIPPED BASES
		}elsif ($c =~ /(.*)S/){
			$$stats{softClipCount} += $1;
			for my $i(1..$1){
				$$stats{byCycle}{softClip}{$R}{$cycle}++;
				$cycle++;
			}
			$readLength += $1;
		### INSERTED BASES	
		}elsif ($c =~ /(.*)I/){	   
			$$stats{insertCount} += $1;
			for my $i(1..$1){
				$$stats{byCycle}{insertion}{$R}{$cycle}++;
				$cycle++;
			}
			$readLength += $1;
		### DELETED BASES	
		}elsif ($c =~ /(.*)D/){ 
			$$stats{deletionCount} += $1;
			for my $i(1..$1){
				$$stats{byCycle}{deletion}{$R}{$cycle}++;
			}
			$mappedBases += $1;
		}else{
			die "Can't handle CIGAR operation: $c\n";
		}
	}
	return ($readLength,$mappedBases);
		
}

=pod

=over 2

=item md_stats()

=back

=over 4

=item Arguments: 

$R : R1,R2,R?

$mdstring : the mismatch string, extracted from the MD:Z tag in the sam record

$strand : to which strand does the read map

$stats : reference to the stats hash, values will be modified

=item Returns: 

the number of mismatches defined by the MD string

=back

Description
   
=cut

sub md_stats{
	my ($R,$mdstring,$strand,$stats)=@_;
	my @Op=procMD($mdstring,$strand);
	my $cycle=1;
	my $mm=0;
	
	for my $md(@Op){
		if ($md =~ /^([0-9]+)/){
			# ignore matching bases
			$cycle += $1;
		}
		elsif ($md =~ /^(\^[A-Z]+)/){
			# ignore deletions
		}
		elsif ($md =~ /^([A-Z]+)/){		# mismatch!
			foreach my $i (split(//, $1)){
				$$stats{mismatchCount}++;
				$$stats{byCycle}{mismatch}{$R}{$cycle}++;
				$cycle++;
			   	$mm++;
			}
		}
		else{
				die "Couldn't handle MD operation: $md\n";
		}
	}
	
	return $mm;
}

=pod

=over 2

=item onTarget()

=back

=over 4

=item Arguments: 

$chrom : the chromosome to which the read maps

$start : the leftmost start position

$mapped : the number of bases mapped

$bed : a reference to the bed hash

$stats : a reference to the stats hash containing the bed intervals, which will be modified

=item Returns: 

1 if mapped to Target, 0 if not

=back

Determines if the read overlaps to any degree intervals contained within the bedfile record.  The bedfile record is modified to indicated mapping to this interval
   
=cut

sub onTarget{
	my($chrom,$start,$mapped,$stats)=@_;
	my $onTarget=0;
	
######### WHAT IS AN ONTARGET READ
#####---------------XXXXXXXXXXXX------------
#####  0000000000
#####          1111111111
#####             1111111111111111
#####                 11111111
##### 					  1111111111111
##### 								0000000000
	
	## bed intervals are hash:array:hash 
	## get array of intervals from the current chromosome
	if($$stats{bed}{intervals}{$chrom}){
		my @intervals=@{$$stats{bed}{intervals}{$chrom}};
		### extract the indices of intervals which meet the necessary conditions
		### map start position is less than the interval stop
		### map start + mapped bases is greater than the interval start
		my @idxs=grep{ 
					($start <=$intervals[$_]{Stop}) && ( ($start+$mapped) >= $intervals[$_]{Start}) 
				} (0..$#intervals);

		my $idxcount=scalar @idxs;  ## should only be ONE interval that contains the start
		
		
		### the read may map to multiple intervals
		### can either build hist on the 1st interval, or on all intervals
		### I think the original script only did the first
		for my $idx(sort{$a<=>$b} @idxs){
			next if($onTarget);  ### will turn this off if it should be mapping to all
			## there is overlap between the read and the interval, but not all bases are necessarily mapped to the interval
			## the logic here indicates that all mapp
			$$stats{bed}{intervals}{$chrom}[$idx]{hist}+=$mapped;
			$onTarget=1;
	
		}
	}
	return $onTarget;
}

=pod

=over 2

=item addRunningBaseCoverage()

=back

=over 4

=item Arguments: 

$chrom : the chromosome to which the current read maps

$start1 : the left-most start of the paired end reads

$start2 : the right-most start of the paired end reads

$cigar : the cigar string

$strand : the strand to which the read maps

$stats : a reference to the stats hash

=item Returns: the number of positions to which the current fragment is stored/mapped

=back

Adds fragments to the runningBaseCoverage collection.  Each fragment is added to all positions to which the read maps.
This collection will be continuously process and cleared of all positions that precede the current read start
   
=cut

sub addRunningBaseCoverage{  
	my ($chrom,$start1,$start2,$cigar,$strand,$stats)=@_;
	my @Op=procCigar($cigar,$strand);
	
	my $posOffset=0;
	my $fragment="$chrom\t$start1\t$start2";

	### first review the cigar string, for matches or deletions
	foreach my $c (@Op){
		if ($c =~ /(.*)[M|D]/){   ## MISMATCH or DELETIONS
			for my $i(1..$1){
				$$stats{runningBaseCoverage}{$chrom}{$start1 + $posOffset}{$fragment}++;
				$posOffset++;
			}
		}
	}
	return $posOffset;
}
	


=pod

=over 2

=item RunningBaseCoverage()

=back

=over 4

=item Arguments: 

$stats : a reference to the stats hash, holding the runningBaseCoverage hash.  Other keys in the stats hash will be modified

$chrom : the current chromosome, all recorded positions not on this chromosome will be processed and cleared.  If no value, then the whole hash is cleared

$pos : the current position, all recorded positions on the current chromosome, before this position, will be cleared

=item Returns: The number of cleared positions from the runningBaseCoverage hash

=back

Calculates base coverage as a running total. This requires that the streaming records are from a sorted bam file.
This will produce a hash with chromosomal keys, each value being a hash of fragment positions and counts, by parsing the cigar string for the current record.
When a new mapping position is found, stats are generated on all stored positions, and then the hash is cleared

   
=cut

sub runningBaseCoverage{  
	my ($stats,$chrom,$startpos)=@_;

	my $cleared_positions=0;
	### will clear out all startpoints stored that precede the current chromosome and position
	for my $chr (keys %{$$stats{runningBaseCoverage}}){
		for my $pos (sort {$a<=>$b} keys %{ $$stats{runningBaseCoverage}{$chr} }){  ### check each start position
			
			##process if not the current chromosome or current chromosome + before current position, or if the $chrom OR $pos variables are not supplied
			if( !$chr || !$startpos || ($chr ne $chrom) || ($pos<$startpos) ) {
				my %startPoints=%{$$stats{runningBaseCoverage}{$chr}{$pos}};

			    ### how many start points/unique fragments are stored on this chromosome positionn
				my $numStartPoints = scalar keys %startPoints;
				
				### collapsedCoverageHist : histogram of number of Unique fragments at a given position
				for my $i(1..$numStartPoints){
					$$stats{collapsedCoverageHist}{$i}++;
				}
				
				### nonCollapsedCoverageHist : histogram of number of fragments at a given position
				my $totalDepth=0;
				map{ $totalDepth+=$_ } values %startPoints;
				
				for my $i(1..$totalDepth){
					$$stats{nonCollapsedCoverageHist}{$i}++;
				}
				
				delete $$stats{runningBaseCoverage}{$chr}{$pos};  ### delete this position from the running BaseCoverageHash
				$cleared_positions++;
			}
			
			
		}
		delete $$stats{runningBaseCoverage}{$chr} if(!$chr || ($chr ne $chrom));  ### delete previous chromosome keys
	}
	
	return $cleared_positions;
}

=pod

=over 2

=item insertMapping()

=back

=over 4

=item Arguments: 

$tlen = template length

$rnext = the chromosome of the other in the pair, = indicates the same

$hash = reference to a hash collecting statistics

$p = parameter hash

=item Returns: 

A description of the insert mapping (normlInsertSize/pairsMappedAbnormallyFar/pairsMappedToDifferentChr).  Modifies the reference hash

=back

Identifies the paired-end insert size as being within the normal range, abnormally far, or on different chromosomes
   
=cut
	
sub insertMapping{
	my($tlen,$rnext,$hash,$p)=@_;    
	my $class="";
	if(($tlen>0) && ($rnext eq "=")){  ### maps to the same chromosome
		if($tlen<$$p{normalInsertMax}){   
			### within expected size , will store the template length for histogram construction
			$$hash{normalInsertSizes}{$tlen}++;
			$class="normalInsertSize";
		}else{
			$$hash{pairsMappedAbnormallyFar}++;
			$class="pairsMappedAbnormallyFar";
		}	
	}elsif($rnext ne "*"){              ### mapped to a different chromosome
		$$hash{pairsMappedToDifferentChr}++;  
		$class="pairsMappedToDifferentChr"; 
	}
	return $class;
}

=pod

=over 2

=item generate_jsonHash()

=back

=over 4

=item Arguments: 

$stats = reference to the stats hash, where data is stored

$p = reference to the parameters hash

=item Returns: 

The jsonHash, ready to print

=back

Transforms information in the stats hash to the appropriate keys and values for the jsonHash
   
=cut

sub generate_jsonHash{
	my($stats,$p)=@_;
	
	
	my %jsonHash=map{ ($_,$$stats{$_} ) } (	"number of ends",
  								"total reads","mapped reads","unmapped reads","non primary reads","paired reads",
						    	"properly paired reads","mate unmaped reads","qual fail reads","qual cut",
						    	"aligned bases","mismatch bases","inserted bases","deleted bases","soft clip bases","hard clip bases",
								"reads per start point","reads on target","target file","target size","number of targets",
							);
							
	$jsonHash{"total reads"} = $$stats{"total reads"} * 1;						
							
	$jsonHash{"average read length"} = $$stats{averageReadLength}{overall};	
	$jsonHash{"insert mean"} = $$stats{meanInsert} || "0";      ## this is quoted to be consistent with samStats.pl
	$jsonHash{"insert stdev"} = $$stats{stdevInsert} || "0";	## this is quoted to be consistent with samStats.pl				
	### capture any jsonHash elements storeed in the parameters
	
	map{
		$jsonHash{$_}=$$p{jsonHash}{$_};
	} keys %{$$p{jsonHash}};

	if($$p{reportBasesCovered}){
		$jsonHash{"non collapsed bases covered"} = $$stats{nonCollapsedCoverageHist};
		$jsonHash{"collapsed bases covered"}     = $$stats{collapsedCoverageHist};
	}
	$jsonHash{"insert histogram"} = $$stats{normalInsertSizes} || {};
	
	for my $R(qw/R1 R2 R?/){
		(my $read=$R)=~s/R/read /;
		$jsonHash{"$read quality histogram"}  	= $$stats{qualHist}{$R} || {};
		$jsonHash{"$read quality by cycle"}   	= $$stats{qualLine}{$R};   ### this will get an undef value, not empty hash, to match samStats.pl
		$jsonHash{"$read length histogram"} 	= $$stats{readLengthHist}{$R} || {};
		$jsonHash{"$read average length"} 		= $$stats{averageReadLength}{$R} || 0;
		$jsonHash{"$read aligned by cycle"} 	= $$stats{byCycle}{aligned}{$R} || {};
		$jsonHash{"$read mismatch by cycle"} 	= $$stats{byCycle}{mismatch}{$R} || {};
		$jsonHash{"$read insertion by cycle"} 	= $$stats{byCycle}{insertion}{$R} || {};
		$jsonHash{"$read deletion by cycle"} 	= $$stats{byCycle}{deletion}{$R} || {};
		$jsonHash{"$read soft clip by cycle"} 	= $$stats{byCycle}{softClip}{$R} || {};
		$jsonHash{"$read hard clip by cycle"} 	= $$stats{byCycle}{hardClip}{$R} || {};
	}
	
	$jsonHash{"number of ends"} = $$stats{"properly paired reads"}>0 ? "paired end" : "single end";
	
	$jsonHash{"number of targets"} = $$p{bed}{numberOfTargets};
	$jsonHash{"qual cut"} = $$p{qualCut};
	$jsonHash{"reads per start point"} = $$stats{startPoint}{RPSP};
	$jsonHash{"target size"} = $$p{bed}{targetSize};
	$jsonHash{"target file"} = $$p{bedFile};
	
	return %jsonHash;	

}





## INTERNAL procCigar()
## Arguments: 
## 		$cigar : the cigar string
## 		$strand : the strand to which the read maps
## Returns : an array of cigar pieces
## Description : breaks the cigar string into pieces denoted by length/type, accounting for strand
sub procCigar{
    my ($cigar,$strand) = @_;
    my @cigarOp;

    while ($cigar =~ /^([0-9]+[MIDNSHPX=]).*$/){
        push (@cigarOp, $1);
        $cigar =~ s/$1//;
    }
    @cigarOp=reverse(@cigarOp) if($strand eq "-");

    return @cigarOp;
}

## INTERNAL procMD()
## Arguments
##	$md : the mdstring
##  $strand : the strand to which the read maps
## Returns : an array of mdstring pieces
## Descriptions : breaks the mdstring into pieces, accounting for strand
sub procMD{
    my ($md,$strand) = @_;
    my @mdOp;
    while ($md ne ""){
        if ($md =~ /^([0-9]+)/)
        {
            push(@mdOp, $1);
            $md =~ s/^$1//;
        }
        if ($md =~ /^([A-Z]+)/)
        {
            push(@mdOp, $1);
            $md =~ s/^$1//;
        }
        if ($md =~ /^(\^)([A-Z]+)/)
        {
            push(@mdOp, "^$2");
            $md =~ s/^\^$2//;
        }
    }
    
    @mdOp=reverse(@mdOp) if($strand eq "-");

    return @mdOp;
}






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
sub load_json{
	
	my @files=@_;
	my %hash;
	for my $file (@files){
		print STDERR "reading from $file\n";
		open (my $FILE,"<",$file) or die "Couldn't open $file.\n";
		if (my $line = <$FILE>){
			$hash{basename($file)} = decode_json($line);
		}else{
			warn "No data found in $file!\n";
		}
	}
	return %hash;
}
sub plot_data{
	my ($j,$scriptPath)=@_;
	for my $rpt (keys %$j){
		#my $title=exists $$j{$rpt}{"barcode"}? 
		#		$$j{$rpt}{"run name"} . " Lane: " . $$j{$rpt}{"lane"} . " Barcode: " . $$j{$rpt}{"barcode"} . "\\n" . $$j{$rpt}{"library"}
		#		$$j{$rpt}{"run name"} . " Lane: " . $$j{$rpt}{"lane"} . "\\n" . $$j{$rpt}{"library"};
		warn "graphing $rpt\n";
		my $rv=`$scriptPath/jsonToGraphs.pl $rpt`;
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
