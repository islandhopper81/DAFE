use Test::More tests => 6;

BEGIN {
use_ok( 'Aggregate' );
use_ok( 'DaFilter' );
use_ok( 'DaTable' );
use_ok( 'DecoupleDa' );
use_ok( 'Justify' );
use_ok( 'Param_handler' );
}

diag( "Testing DAFE $DAFE::VERSION" );
