#!/usr/bin/perl -w
use strict;
use warnings;
use Test::More;
use Text::Diff;
use GSI::bamqc 'load_json';
use File::Path 'remove_tree';
use File::Slurp 'read_file';
use File::Basename 'dirname';

my $script_dir=dirname($0);

sub run_bamqc {
    my($opts,$output_dir)=@_;
    my $TESTS_FAIL=0;
    
 
    my $actual_json = "$output_dir/actual_output.json";
    my $actual_txt = "$output_dir/actual_output.txt";
    my $expected_json = "$output_dir/expected_output.json";
    my $expected_txt = "$output_dir/expected_output.txt";

    my $bamqctest="samtools view t/test/neat_5x_EX_hg19_chr21.bam | perl bamqc.pl -r t/test/SureSelect_All_Exon_V4_Covered_Sorted_chr21.bed $opts > $actual_json 2> $actual_txt";
    print "################# Running bamqc in $output_dir with options $opts\n#\t$bamqctest\n";   

    is ( system($bamqctest), 0, "bamqc.pl test with opts: '$opts' returns 0 exit status");
    ok ( -e $actual_json , "the json output exists at $actual_json");
    ok ( -e $actual_txt, "the text output exists at $actual_txt");
    
    #expected error has full path instead of relative path, so get rid of that
    system("perl -p -i -e 's#$script_dir/##g' $actual_txt");
    #check that the expected and actual text files are the same
    my $diff = diff $actual_txt, $expected_txt, { CONTEXT=>1, STYLE=>'OldStyle' };
    $TESTS_FAIL=1 if not is ( $diff, '', "text files are the same: $actual_txt and $expected_txt");
    
    #check that the expected and actual JSON files are the same
    my %jsons=load_json($actual_json, $expected_json);
    $TESTS_FAIL=1 if not is_deeply ( $jsons{$actual_json}, $jsons{$expected_json}, "json files are the same: $actual_json and $expected_json"); 
    
    #if the file comparison tests fail, save the results. Otherwise, remove them
    if ($TESTS_FAIL) {
        print "Differing file is saved at $actual_txt and $actual_json \n";
        BAIL_OUT("Fix the above problems before proceeding with further tests.");
    } else {
        unlink $actual_txt;
        unlink $actual_json;
    }

    return $TESTS_FAIL;
}

#### test start

my $metadata_file='t/test/metadata.json';
my $metadata_str=read_file($metadata_file);

run_bamqc("-s 2000 -i 1500 -q 30 -j $metadata_file ","t/test/bamqc_vanilla");

run_bamqc("-s 2000 -i 1500 -q 30 -j \"$metadata_str\" ","t/test/bamqc_vanilla_jsonstring"); 

run_bamqc("-c -i 1500 -q 30 -j $metadata_file","t/test/bamqc_coverage");

run_bamqc("-H t/test/bamqc_histo_coverage/test.hist -i 1500 -q 30 -s 2000 -j $metadata_file ", "t/test/bamqc_histo_coverage");

done_testing();
