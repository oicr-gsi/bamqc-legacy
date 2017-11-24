package GSI::bamqc;
use strict;
use warnings;

BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

		$VERSION	=	1.00;
		@ISA		=	qw(Exporter);
		#@EXPORT	=	qw();
    @EXPORT_OK   = qw(read_bed assess_start_point assess_flag cigar_stats md_stats
							onTarget addRunningBaseCoverage runningBaseCoverage HistStats insertMapping
							load_json load_json_and_dirs toPhred generate_jsonHash
              generate_mismatch_rate generate_indel_rate generate_softclip_rate
              generate_hardclip_rate generate_error_rate get_barcode get_group
              get_raw_reads get_raw_yield get_map_percent get_ontarget_percent
              get_est_yield get_est_coverage);
    #%EXPORT_TAGS = ();
}

sub new{
    my ($class, %parameters) = @_;
    my $self = bless ({}, ref ($class) || $class);
    return $self;
}


#################### main pod documentation begin ###################


=head1 NAME

bamqc - Generate quality control statistics from BAM files

=head1 SYNOPSIS

  use GSI::bamqc;


=head1 DESCRIPTION


=head1 USAGE

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


use Exporter;
use JSON::PP;
use File::Basename;
### for debug
#use Data::Dumper;
#(open my $TTY,"/dev/tty") || die "unable to open keyboard input";

my @month = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );



=head2 assess_flag($flag,$stats,$qual,$qcut)

Assesses the sam flag and stores information in the reference hash

Returns   : 1 if the read is mapped, 0 if not

Argument  :	$flag = the sam record flag;
						$stats = a reference to a hash, which will be modified;
						$qual = the alignment quality;
						$param = the parameter hash;

Comment   : conditions assessed are as follows:
						256="non primary reads";
						4="unmapped read";
						$qual < $qcut="qual fail reads";
						if none of the above="mapped reads";
						1="paired reads";
						8="mate unmaped reads" [sic];
						2="properly paired reads"

=cut

sub assess_flag{
	my($flag,$stats,$qual,$qcut)=@_;
	my $mapped=0;
	### processing ALL reads, not just sampled
	$$stats{"total reads"}++;


	if ($flag & 256) {
		$$stats{"non primary reads"}++;
	}elsif ($flag & 4){
		$$stats{"unmapped reads"}++;
	}elsif ($qual < $qcut	){
		$$stats{"qual fail reads"}++;
	}else{
		$mapped=1;
		$$stats{"mapped reads"}++;
		if ($flag & 1){
			$$stats{"paired reads"}++;
			if ($flag & 8){
				$$stats{"mate unmaped reads"}++
			}
		}
		if ($flag & 2){
			$$stats{"properly paired reads"}++
		}
	}
	return $mapped;
}


=head2 assess_start_point($chrom,$s1,$s2,$sphash)

Keeps a running calculation of ReadsPerStartPoint, updating this
						value with each read.  It is dependent on the stream alignment data
						being properly sorted

Returns   : A revised StartPoint hash, with updated stats

Argument  :	$chrom = the chromosome of the current record;
						$s1	 = the position of the current record = the leftmost start point;
						$s2    = the position of the paired record = the rightmost start point;
						$sphash	 = the current StartPoint hash, holding the current position;
										and the running ReadsPerStartPoint value;

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





=head2 read_bed($file)

Reads in a bed file, storing in the intervals in an array of hashes
						indicating, Start/Stop and size

Returns   : Hash structure containing the bed intervals

Argument  :	$file = the name of the bed file, containing the intervals

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

=head2 cigar_stats($chrom,$start1,$start2,$R,$strand,$cigar,$stats,$p)

Processes the cigarString, by breaking into pieces and assessing how
						the read is mapped to the reference

Returns   : A list with two values, $readLength and $mappedBases

Argument  :	$chrom = the chromosome to which the current read maps;
						$start1 = the leftmost startpoint;
						$start2 = the rightmost startpoint for the paired end;
						$R = the read (R1,R2,R?);
						$strand = the strand to which the read maps, from the sam flag;
						$cigar = the cigar string;
						$stats = a reference to the stats hash, which will be modified;
						$p = a reference to the parameter hash;

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

=head2 md_stats($R,$mdstring,$strand,$stats)

Processes the cigarString, by breaking into pieces and assessing how
						the read is mapped to the reference

Returns   : the number of mismatches defined by the MD string

Argument  :	$R = R1,R2,R?;
						$mdstring = the mismatch string, extracted from the MD:Z tag in the sam record;
						$strand = to which strand does the read map;
						$stats = reference to the stats hash, values will be modified;

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

=head2 onTarget($chrom,$start,$mapped,$stats)

Determines if the read overlaps to any degree intervals contained
						within the bedfile record.  The bedfile record is modified to
						indicated mapping to this interval

Returns   : 1 if mapped to Target, 0 if not

Argument  :	$chrom = the chromosome to which the read maps;
						$start = the leftmost start position;
						$mapped = the number of bases mapped;
						$bed = a reference to the bed hash;
						$stats = a reference to the stats hash containing the bed intervals, which will be modified;

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


=head2 addRunningBaseCoverage($chrom,$start1,$start2,$cigar,$strand,$stats)

Adds fragments to the runningBaseCoverage collection.  Each fragment
						is added to all positions to which the read maps.
						This collection will be continuously process and cleared of all
						positions that precede the current read start

Returns   : the number of positions to which the current fragment is stored/mapped

Argument  :	$chrom = the chromosome to which the current read maps;
						$start1 = the left-most start of the paired end reads;
						$start2 = the right-most start of the paired end reads;
						$cigar = the cigar string;
						$strand = the strand to which the read maps;
						$stats = a reference to the stats hash;

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


=head2 runningBaseCoverage($stats,$chrom,$startpos)

Calculates base coverage as a running total. This requires that the
						streaming records are from a sorted bam file.
						This will produce a hash with chromosomal keys, each value being a
						hash of fragment positions and counts, by parsing the cigar string
						for the current record.
						When a new mapping position is found, stats are generated on all
						stored positions, and then the hash is cleared

Returns   : The number of cleared positions from the runningBaseCoverage hash

Argument  :	$stats = a reference to the stats hash, holding the runningBaseCoverage hash.  Other keys in the stats hash will be modified;
						$chrom = the current chromosome, all recorded positions not on this chromosome will be processed and cleared.  If no value, then the whole hash is cleared;
						$pos = the current position, all recorded positions on the current chromosome, before this position, will be cleared;


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
		delete $$stats{runningBaseCoverage}{$chr} if(!$chr || ((defined $chrom) and ($chr ne $chrom)));  ### delete previous chromosome keys
	}

	return $cleared_positions;
}

=head2 HistStats(%val)

Calculate the mean and standard deviation of insert sizes in a hash

Returns   : mean and standard deviation of the insert sizes

Argument  :	%val = a hash, with keys = insert size, values = count for each insert size

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


=head2 insertMapping($tlen,$rnext,$hash,$p)

Identifies the paired-end insert size as being within the normal
						range, abnormally far, or on different chromosomes

Returns   : A description of the insert mapping (normlInsertSize/
						pairsMappedAbnormallyFar/ pairsMappedToDifferentChr).
						Modifies the reference hash

Argument  :	$tlen = template length;
						$rnext = the chromosome of the other in the pair, = indicates the same;
						$hash = reference to a hash collecting statistics;
						$p = parameter hash;

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








=head2 load_json(@files)

Open, decode, and store each JSON file in a hash with filename keys

Returns   : a hash with filename keys > decoded JSON hash

Argument  :	@files = a list of json file paths

=cut
sub load_json{

	my @files=@_;
	my %json_hash;
	for my $file (@files){
		print STDERR "reading from $file\n";
		open (my $FILE,"<",$file) or die "Couldn't open $file.\n";
		if (my $line = <$FILE>){
			$json_hash{basename($file)} = decode_json($line);
            $json_hash{basename($file)}{basename}=basename($file);
            $json_hash{basename($file)}{dirname}=dirname($file);
		}else{
			warn "No data found in $file!\n";
		}
	}
	return %json_hash;
}

=head2 load_json_and_dirs(@files)

Open, decode, and store each JSON file in a hash with filename keys

Returns   : two hashes: a hash with filename keys > decoded JSON hash; and
												a hash with filename keys > directory

Argument  :	@files = a list of file paths to open

=cut
sub load_json_and_dirs{

	my @files=@_;
	my %json_hash;
  my %dir_hash;
	for my $file (@files){
		print STDERR "reading from $file\n";
		open (my $FILE,"<",$file) or die "Couldn't open $file.\n";
		if (my $line = <$FILE>){
			$json_hash{basename($file)} = decode_json($line);
            $json_hash{basename($file)}{basename}=basename($file);
            $json_hash{basename($file)}{dirname}=dirname($file);
            $dir_hash{basename($file)} = dirname($file);
		}else{
			warn "No data found in $file!\n";
		}
	}
	return (\%json_hash, \%dir_hash);
}

=head2 toPhred($char)

Convert a character to a phred score (ascii value - 33)

Returns   : Phred score

Argument  :	$char = the single character to convert

=cut
sub toPhred{
	my $char = $_[0];
	my $ascii = ord($char);
	my $offset = 33;
	return $ascii - $offset;
}

=head2 generate_jsonHash($stats,$p)

Transforms information in the stats hash to the appropriate keys and
						values for the jsonHash

Returns   : The jsonHash, ready to print

Argument  :	$stats = reference to the stats hash, where data is stored;
						$p = reference to the parameters hash;

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

=head2 generate_error_rate( $hash, $prefix )

Compute the error rate with (mismatch+insertion+deletion)/aligned

Returns   : A formatted percentage

Argument  :	$hash = The hash containing the JSON file contents to analyse;
						$prefix = A prefix to use when accessing the keys of hash;

=cut
sub generate_error_rate {
    my ( $hash, $prefix ) = @_;
    return generateRatePercent( $hash, $prefix,
             ["mismatch", "insertion", "deletion" ],
            ["aligned"]
        );
}

=head2 generate_mismatch_rate( $hash, $prefix )

Compute the mismatch rate with mismatch/aligned

Returns   : A formatted percentage

Argument  :	$hash = The hash containing the JSON file contents to analyse;
						$prefix = A prefix to use when accessing the keys of hash;

=cut
sub generate_mismatch_rate {
    my ( $hash, $prefix ) = @_;
    return generateRatePercent( $hash, $prefix,
            ["mismatch"],
            ["aligned"] );
}

=head2 generate_indel_rate( $hash, $prefix )

Compute the indel rate with (insertion+deletion)/aligned

Returns   : A formatted percentage

Argument  :	$hash = The hash containing the JSON file contents to analyse;
						$prefix = A prefix to use when accessing the keys of hash;

=cut
sub generate_indel_rate {
    my ( $hash, $prefix ) = @_;
   return generateRatePercent( $hash, $prefix,
            [ "insertion", "deletion" ],
            ["aligned"]);
}

=head2 generate_softclip_rate( $hash, $prefix )

Compute the softclip rate with soft clip/(soft clip + aligned)

Returns   : A formatted percentage

Argument  :	$hash = The hash containing the JSON file contents to analyse;
						$prefix = A prefix to use when accessing the keys of hash;

=cut
sub generate_softclip_rate {
    my ( $hash, $prefix ) = @_;
    return generateRatePercent( $hash, $prefix,
            ["soft clip"],
            [ "soft clip", "aligned" ]);
}

=head2 generate_hardclip_rate( $hash, $prefix )

Compute the hardclip rate with hard clip/(hard clip + soft clip + aligned)

Returns   : A formatted percentage

Argument  :	$hash = The hash containing the JSON file contents to analyse;
						$prefix = A prefix to use when accessing the keys of hash;

=cut
sub generate_hardclip_rate {
    my ( $hash, $prefix ) = @_;
    return generateRatePercent( $hash, $prefix,
            ["hard clip"],
            [ "soft clip", "aligned", "hard clip" ]);
}

=head2 get_barcode( $jsonHash)

Get the sequencing index / barcode

Returns   : the barcode if it exists; otherwise 'NoIndex'

Argument  :	$jsonHash = The hash containing the JSON file contents to analyse.

=cut
sub get_barcode {
    my ( $jsonHash ) =@_;
    return exists $jsonHash->{barcode} ? $jsonHash->{barcode} : 'NoIndex';
}

=head2 get_group( $jsonHash)

Get the group id and group id description

Returns   : the group id and group id description separated by a space, if they
						exist. If they don't exist, then 'na'

Argument  :	$jsonHash = The hash containing the JSON file contents to analyse.

=cut
sub get_group {
    my ($jsonHash)=@_;
    return exists $jsonHash->{"group id"}
            ? exists $jsonHash->{"group id description"}
                ? $jsonHash->{"group id description"}." ".$jsonHash->{'group id'}
                : $jsonHash->{'group id'}
            : "na";
}

=head2 get_raw_reads( $jsonHash)

Get the total raw reads, counted by summing mapped, unmapped and qual fail reads

Returns   : the total number of reads (int)

Argument  :	$jsonHash = The hash containing the JSON file contents to analyse.

=cut
sub get_raw_reads {
    my ($jsonHash)=@_;
    return int($jsonHash->{"mapped reads"} + $jsonHash->{"unmapped reads"} + $jsonHash->{"qual fail reads"});
}

=head2 get_raw_yield( $jsonHash)

Get the total raw yield, multipying get_raw_reads by the average read length

Returns   : the total yield (int)

Argument  :	$jsonHash = The hash containing the JSON file contents to analyse.

=cut
sub get_raw_yield {
    my ($jsonHash)=@_;
    return int(GSI::bamqc::get_raw_reads($jsonHash) * $jsonHash->{"average read length"});
}

=head2 get_map_percent( $jsonHash)

Get the total map percentage, calculated by dividing mapped reads by
						total number of reads and multiplying by 100. If total reads is 0,
						treat as 1.

Returns   : the map percentage as an integer

Argument  :	$jsonHash = The hash containing the JSON file contents to analyse.

=cut
sub get_map_percent {
    my ($jsonHash)=@_;
    return 100 * ($jsonHash->{"mapped reads"} / ( GSI::bamqc::get_raw_reads($jsonHash) || 1 ));
}

=head2 get_ontarget_percent( $jsonHash)

Get the total on target percentage by dividing reads on target by
						the number of mapped reads, and multiplying by 100. If mapped reads
						is 0, treat as 1.

Returns   : the on target percentage as an integer

Argument  :	$jsonHash = The hash containing the JSON file contents to analyse.

=cut
sub get_ontarget_percent {
    my ($jsonHash)=@_;
    return 100 * ($jsonHash->{"reads on target"} / ($jsonHash->{"mapped reads"} || 1) ); # $rawReads) * 100;   # could argue using this either way
}

=head2 get_est_yield( $jsonHash)

Get the estimated total yield by multiplying total aligned based by
						the on target percentage, and dividing that by reads per start point.
						If reads per start point is is 0, treat as 1.

Returns   : the estimated yield as an integer

Argument  :	$jsonHash = The hash containing the JSON file contents to analyse.

=cut
sub get_est_yield {
    my ($jsonHash)=@_;
    return int($jsonHash->{"aligned bases"} *
            (GSI::bamqc::get_ontarget_percent($jsonHash) / 100) /
            ($jsonHash->{"reads per start point"} || 1) );
}

=head2 get_est_coverage( $jsonHash)

Get the estimated total coverage by dividing estimated yield by
						the target size.
						If the target is is 0, treat as 1.

Returns   : the estimated coverage as an integer

Argument  :	$jsonHash = The hash containing the JSON file contents to analyse.

=cut
sub get_est_coverage {
    my ($jsonHash)=@_;
    return GSI::bamqc::get_est_yield($jsonHash) / ($jsonHash->{"target size"} || 1);
}

##################### INTERNAL Subroutines ##############################

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

# Compute a rate from the properties of a run.
#
# 1[HashRef]: The hash containing the JSON file contents to analyse.
# 2[Str]: A prefix to use when accessing the keys of hash.
# 3[ArrayRef[Str]]: The names of the "by cycle" keys that belong in the numerator.
# 4[ArrayRef[Str]]: The names of the "by cycle" keys that belong in the denominator.
# Returns[Str]: A formatted percentage of all the sum of all numerator items divided by the sum of all the denominator items.
sub generateRatePercent {
    my ( $hash, $prefix, $numerator_names, $denominator_names ) = @_;

    my $numerator = 0;
    for my $name ( @{$numerator_names} ) {
        $numerator +=
          byCycleToCount( $hash->{ $prefix . $name . " by cycle" } );
    }
    my $denominator = 0;
    for my $name ( @{$denominator_names} ) {
        $denominator +=
          byCycleToCount( $hash->{ $prefix . $name . " by cycle" } );
    }
    return $numerator * 100.0 / ( $denominator || 1 ) ;
}

=head2 findStart($cigarOp, $start)

Find the start index of the read using the cigar string

Returns   : the new start position, adjusted to take soft clipping into account

Argument  :	$cigarOp = the cigar string

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

=head2 findEnd($cigarOp, $end)

Find the start index of the read using the cigar string

Returns   : the new start position, adjusted to take soft clipping into account

Argument  :	$cigarOp = the cigar string;
						$end = the assumed end of the string;

=cut
sub findEnd{
        my @cigarOp = @{ $_[0] };
        my $end = $_[1];
        if ($cigarOp[0] =~ /(.*)S/) {     # if first cigar operation is a soft clip, adjust start point
                $end -= $1;
        }
        foreach my $cig (@cigarOp) {
                if ($cig =~ /(.*)S/) {
                        $end += $1;
                }
                elsif ($cig =~ /(.*)M/){
                        $end += $1;
                }
                elsif ($cig =~ /(.*)D/){
                        $end += $1;
                }
        }
        return $end;
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

######################## END INTERNAL Subroutines########################

1;
