#!/usr/bin/perl -w
use strict;
use warnings;
use Test::More;
use Text::Diff;
use bamqc 'load_json';
use File::Path 'remove_tree';
use File::Slurp 'read_file';

sub run_bamqc {
    my($opts,$output_dir)=@_;
    my $TESTS_FAIL=0;
    
    print "################# Running bamqc in $output_dir with options $opts\n";   
 
    my $actual_json = "$output_dir/actual_output.json";
    my $actual_txt = "$output_dir/actual_output.txt";
    my $expected_json = "$output_dir/expected_output.json";
    my $expected_txt = "$output_dir/expected_output.txt";

    my $bamqctest="samtools view test/neat_5x_EX_hg19_chr21.bam | perl bamqc.pl -r test/SureSelect_All_Exon_V4_Covered_Sorted_chr21.bed $opts > $actual_json 2> $actual_txt";

    is ( system($bamqctest), 0, "bamqc.pl test with opts: '$opts' returns 0 exit status");
    ok ( -e $actual_json , 'the json output exists');
    ok ( -e $actual_txt, 'the text output exists');

    #check that the expected and actual text files are the same
    my $diff = diff $actual_txt, $expected_txt, { CONTEXT=>1, STYLE=>'OldStyle' };
    $TESTS_FAIL=1 if not is ( $diff, '', "text files are the same");
    
    #check that the expected and actual JSON files are the same
    my %jsons=load_json($actual_json, $expected_json);
    $TESTS_FAIL=1 if not is_deeply ( $jsons{$actual_json}, $jsons{$expected_json}, "json files are the same"); 
    
    #if the file comparison tests fail, save the results. Otherwise, remove them
    if ($TESTS_FAIL) {
        print "Differing file is saved at $actual_txt and $actual_json \n";
    } else {
        unlink $actual_txt;
        unlink $actual_json;
    }

    return $TESTS_FAIL;
}

sub run_genericrunreport {
    my($run_opts,$output_dir)=@_;
    my $TESTS_FAIL=0;

    print "########## Running report in $output_dir with options $run_opts\n";
    
    my $test_json = "$output_dir/test.json";
    my $actual_html="$output_dir/actual_output.html";
    my $expected_html="$output_dir/expected_output.html";

    my $reporttest="perl jsonToGenericRunReport.pl $run_opts $test_json | grep -v 'Generic run report generated on' > $actual_html";
    
    is ( system($reporttest), 0, "jsonToGenericRunReport.pl test with opts: '$run_opts' returns 0 exit status");
    ok ( -e $actual_html , 'the html output exists');
    
    my $diff = diff $actual_html, $expected_html, { CONTEXT=>1, STYLE=>'OldStyle' };
    $TESTS_FAIL=1 if not is ( $diff, '', "html files are the same");
    
    if ($TESTS_FAIL) {
        print "Differing files are saved at $actual_html\n";
    } else {
        unlink $actual_html;
    }
}

sub test_genericrunreport_graphs {
    my($output_dir,$has_coverage)=@_;
    my $TESTS_FAIL=0;
    print "####### testing report graphs in $output_dir\n";

    my $actual_dir="$output_dir/test.json.graphs";
    my $expected_dir="$output_dir/expected.graphs";

    ok ((-e $actual_dir && -d $actual_dir), "Graphs directory exists: $actual_dir");

    my @graphs=("hardCycle.png","misCycle.png","qualHist.png", "softCycle.png", "insert.png", "qualCycle.png","readPie.png", "indelCycle.png", "readLength.png");

    my @optgraphs=();
    
    my $diff;
    foreach my $graph (@graphs) {
        $TESTS_FAIL=test_image_rcode("$actual_dir/$graph", "$expected_dir/$graph");
    }
   
    if ($has_coverage) {
        my @optgraphs=("collapsedCover.png", "nonCollapsedCover.png");
        foreach my $graph (@optgraphs) {
            $TESTS_FAIL=test_image_rcode("$actual_dir/$graph", "$expected_dir/$graph");
        }
    }
 
    if ($TESTS_FAIL) {
        print "Differing files are saved at $actual_dir\n";
    } else {
        remove_tree($actual_dir, safe=>0);
    }

}

sub test_image_rcode {
    my($actual, $expected)=@_;
    my $TESTS_FAIL=0;
    if (ok (-e $actual && -e "$actual.Rcode" , "the image and Rcode $actual exists")) {
        my $diff= diff "$actual.Rcode", "$expected.Rcode", { CONTEXT=>1, STYLE=>'OldStyle' };
        $TESTS_FAIL=1 if not is ( $diff, '', "$actual.Rcode files are the same as expected");
    }
    return $TESTS_FAIL;
}

#### test start

my $metadata_file='test/metadata.json';
my $metadata_str=read_file($metadata_file);

run_bamqc("-s 2000 -i 1500 -q 30 -j $metadata_file ","test/bamqc_vanilla");

run_bamqc("-s 2000 -i 1500 -q 30 -j \"$metadata_str\" ","test/bamqc_vanilla_jsonstring"); 

run_genericrunreport("","test/report_vanilla");

run_genericrunreport("-r -g -p -H","test/report_most_opts");
test_genericrunreport_graphs("test/report_most_opts",0);

#run_bamqc("-c -i 1500 -q 30 -j $metadata_file","test/bamqc_coverage");

run_genericrunreport("-c","test/report_coverage");

run_genericrunreport("-g","test/report_coverage_graphs");
test_genericrunreport_graphs("test/report_coverage_graphs",1);

run_bamqc("-H test/bamqc_histo_coverage/test.hist -i 1500 -q 30 -s 2000 -j $metadata_file ", "test/bamqc_histo_coverage");
run_genericrunreport("-g","test/report_histo_coverage");
test_genericrunreport_graphs("test/report_histo_coverage",1);
done_testing();
