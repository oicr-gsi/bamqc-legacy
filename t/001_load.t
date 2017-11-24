# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 8;

BEGIN {
    use_ok( 'GSI::bamqc' );
    use_ok( 'GSI::report' );
    use_ok( 'GSI::jsonToGraphs' );
    use_ok( 'SeqWare::Report');
}

my $object = GSI::bamqc->new ();
isa_ok ($object, 'GSI::bamqc');
$object = GSI::jsonToGraphs->new();
isa_ok ($object, 'GSI::jsonToGraphs');
$object = GSI::report->new();
isa_ok ($object, 'GSI::report');
$object = SeqWare::Report->new();
isa_ok ($object, 'SeqWare::Report');

