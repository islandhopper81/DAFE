use strict;
use warnings;

use Test::More tests => 15;
use Test::Exception;
use Test::Warn;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
my $logger = get_logger();

# others to include
#use File::Temp qw/ tempfile tempdir /;
#use Cwd 'abs_path';
#use File::Basename;

# helper subroutines



### Begin Testing ###
BEGIN { use_ok( 'DAFE::AbEstimator' ); }

# variables for testing

# test constructor
my $obj = undef;
{
    # test for errors with bad parameters
    # -- there are none because there are no parameters
    
    # test with correct parameters
    lives_ok( sub { $obj = DAFE::AbEstimator->new() },
             "expected to live" );
}

# test private method _mean
{
    my @arr = (2,4);
        
    is( $obj->_mean(\@arr), 3, "_mean()" );
}

# test private method _median
{
    my @arr = (1,2,3,4);
    
    is( $obj->_median(\@arr), 2.5, "_median()" );
    
    my @arr2 = (2,1,3,4);
    is( $obj->_median(\@arr2), 2.5, "_median(unsorted)" );
}

# create some tables for testing
my $tbl1 = Table->new();
$tbl1->_set_col_count(4);
my @col_names = ("A", "B", "C", "D");
$tbl1->_set_col_names(\@col_names);
my @row_names = ("Z", "Y");
foreach my $r ( @row_names ) {
    my @row_vals = (1,2,3,4);
    $tbl1->add_row($r, \@row_vals);
}

my $tbl2 = $tbl1->copy();
my @new_row = (4,3,2,10);
$tbl2->add_row("X", \@new_row);

# tbl2
# Z	1	2	3	4
# Y	1	2	3	4
# X	4	3	2	10

# test calc_abund_est
{
    
    
    # test for mean
    is( $obj->calc_abund_est($tbl1, "mean"), 2.5, "calc_abund_est(mean)" );
    is( $obj->calc_abund_est($tbl2, "mean"), 3.25, "calc_abund_est(mean)" );
    
    # test for median
    is( $obj->calc_abund_est($tbl1, "median"), 2.5, "calc_abund_est(median)" );
    is( $obj->calc_abund_est($tbl2, "median"), 2.5, "calc_abund_est(median)" );
}

# test calc_abund_est_v2
{
    # test for mean
    is( $obj->calc_abund_est_v2($tbl1, "mean"), 2.5, "calc_abund_est(mean)" );
    is( $obj->calc_abund_est_v2($tbl2, "mean"), 3.25, "calc_abund_est(mean)" );
    
    # test for median
    is( $obj->calc_abund_est_v2($tbl1, "median"), 2.5, "calc_abund_est(median)" );
    is( $obj->calc_abund_est_v2($tbl2, "median"), 3, "calc_abund_est(median)" );
}

# test calc_marker_abund
{
	my @ans1 = (1,2,3,4);
	my @ans2 = (2,7/3,8/3,6);
	is_deeply( $obj->calc_marker_abund($tbl1, "mean"), \@ans1, "calc_marker_abund(mean)" );
	is_deeply( $obj->calc_marker_abund($tbl2, "mean"), \@ans2, "calc_marker_abund(mean)" );
}
