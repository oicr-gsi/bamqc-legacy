# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 6;

BEGIN {
    use_ok( 'GSI::bamqc' );
    use_ok( 'GSI::jsonToGraphs' );
    use_ok( 'GSI::RunReport');
}

my $object = GSI::bamqc->new ();
isa_ok ($object, 'GSI::bamqc');
$object = GSI::jsonToGraphs->new();
isa_ok ($object, 'GSI::jsonToGraphs');
$object = GSI::RunReport->new();
isa_ok ($object, 'GSI::RunReport');
