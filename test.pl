#!/usr/bin/perl -w
use strict;
use warnings;
use Test::More;
use Text::Diff;
use bamqc 'load_json';
use Cwd 'abs_path'; 


my $actual_json = abs_path('test/actual_output.json');
my $actual_txt = abs_path('test/actual_output.txt');
my $expected_json = abs_path('test/expected_output.json');
my $expected_txt = abs_path('test/expected_output.txt');


#########################################################

### BAM QC TEST

#########################################################

my $TESTS_FAIL=0;

my $vanilla_test="samtools view test/neat_5x_EX_hg19_chr21.bam | perl bamqc.pl -r test/SureSelect_All_Exon_V4_Covered_Sorted_chr21.bed -j test/metadata.json -s 2000 -i 1500 -q 30 > $actual_json 2> $actual_txt";

is ( system($vanilla_test), 0, 'vanilla bamqc.pl test returns 0 exit status');
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
    print "Differing file is saved at $actual_txt\n";
} else {
    unlink $actual_txt;

}


###########################################################

### Run Report test

###########################################################
$TESTS_FAIL=0;

my $actual_html=abs_path('test/actual_output.html');
my $expected_html=abs_path('test/expected_output.html');

$vanilla_test="perl jsonToGenericRunReport.pl $actual_json | grep -v 'Generic run report generated on' > $actual_html";

is ( system($vanilla_test), 0, 'vanilla jsonToGenericRunReport.pl test returns 0 exit status');
ok ( -e $actual_html , 'the html output exists');

$diff = diff $actual_html, $expected_html, { CONTEXT=>1, STYLE=>'OldStyle' };
$TESTS_FAIL=1 if not is ( $diff, '', "text files are the same");

if ($TESTS_FAIL) {
    print "Differing files are saved at $actual_json and $actual_html\n";
} else {
    unlink $actual_html;
    unlink $actual_json;

}



done_testing();
