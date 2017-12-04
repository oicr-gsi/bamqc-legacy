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

sub run_genericrunreport {
    my($run_opts,$output_dir)=@_;
    my $TESTS_FAIL=0;


    my $test_json = "$output_dir/test.json";
    my $actual_html="$output_dir/actual_output.html";
    my $actual_tsv="$output_dir/Simulated_Run_report.tsv";
    my $expected_html="$output_dir/expected_output.html";
    my $expected_tsv="$output_dir/expected_report.tsv";
    my $json_dir=dirname($actual_html);

    my $reporttest="perl jsonToGenericRunReport.pl $run_opts $test_json | grep -v 'Generic Run Report' | grep -v 'SeqWare Browser generated ' | grep -v '<!-- Executed as:' | perl -p -e 's#$json_dir/##g' > $actual_html";
    print "########## Running report in $output_dir with options $run_opts\n#\t$reporttest\n";

    is ( system($reporttest), 0, "jsonToGenericRunReport.pl test with opts: '$run_opts' returns 0 exit status");
    ok ( -e $actual_html , "the html output exists at $actual_html");

    my $diff = diff $actual_html, $expected_html, { CONTEXT=>1, STYLE=>'OldStyle' };
    $TESTS_FAIL=1 if not is ( $diff, '', "html files are the same: $actual_html and $expected_html");

    ok ( -e $actual_tsv , "the html output exists at $actual_tsv");
    $diff = diff $actual_tsv, $expected_tsv, { CONTEXT=>1, STYLE=>'OldStyle' };
    $TESTS_FAIL=1 if not is ( $diff, '', "tsv files are the same: $actual_tsv and $expected_tsv");

    if ($TESTS_FAIL) {
        print "Differing files are saved at $actual_html and $actual_tsv\n";
    } else {
        unlink $actual_html, $actual_tsv;
    }

    if ($run_opts =~ m/g/) {
        my $has_cov=($run_opts =~ m/c/);
        test_genericrunreport_graphs($output_dir,$has_cov);
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

run_genericrunreport("","t/test/report_vanilla");

run_genericrunreport("-r -g -p -H -n","t/test/report_most_opts");

run_genericrunreport("-c","t/test/report_coverage");

run_genericrunreport("-c -g","t/test/report_coverage_graphs");

# runs on a json file that has coverage stats
run_genericrunreport("-g","t/test/report_histo_coverage");

run_genericrunreport("-cgn","t/test/report_nocollapse");

done_testing();
