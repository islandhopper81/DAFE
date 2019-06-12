use strict;
use warnings;

use Test::More tests => 18;
use Test::Exception;
use Test::Warn;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
my $logger = get_logger();

# others to include
use File::Temp qw/ tempfile tempdir /;
use Cwd 'abs_path';
use File::Basename;

# helper subroutines



### Begin Testing ###
BEGIN { use_ok( 'DAFE::Count::OutputRestruct' ); }

# variables for testing
my $ref_aref = ["2503538034", "99000"];
my $sample_aref = ["1119441", "1119442"];
my $perc_ids_aref = ["60", "70", "80", "90", "95"];
my $out_dir = dirname(abs_path($0)) . "/../../../t/test_dir/test_OutputRestruct/ref_out/";
my $new_out_dir = tempdir();
print "new_out_dir: $new_out_dir\n";
my $htseq_file_prefix = "htseq_count";

# test constructor
my $obj = undef;
{
    # test for errors with bad parameters
    throws_ok( sub { DAFE::Count::OutputRestruct->new({}) },
              'MyX::Generic::Undef::Param', "new() - caught" );
    
    # test with correct parameters
    lives_ok( sub { $obj = DAFE::Count::OutputRestruct->new({
                                out_dir => $out_dir,
                                new_out_dir => $new_out_dir,
                                ref_aref => $ref_aref,
                                sample_aref => $sample_aref,
                                perc_ids_aref => $perc_ids_aref,
                                htseq_file_prefix => $htseq_file_prefix}) },
             "expected to live" );
}

# Test the simple getter methods
{
    lives_ok( sub { $obj->get_param("out_dir") }, "expected to live" );
    is ($obj->get_param("out_dir"), $out_dir, "get_param(out_dir)" );
    
    lives_ok( sub { $obj->get_param("new_out_dir") }, "expected to live" );
    is ($obj->get_param("new_out_dir"), $new_out_dir, "get_param(new_out_dir)" );
    
    lives_ok( sub { $obj->get_param("ref_aref") }, "expected to live" );
    is ($obj->get_param("ref_aref"), $ref_aref, "get_param(ref_aref)" );
    
    lives_ok( sub { $obj->get_param("sample_aref") }, "expected to live" );
    is ($obj->get_param("sample_aref"), $sample_aref, "get_param(sample_aref)" );
    
    lives_ok( sub { $obj->get_param("perc_ids_aref") }, "expected to live" );
    is ($obj->get_param("perc_ids_aref"), $perc_ids_aref, "get_param(perc_ids_aref)" );
    
    lives_ok( sub { $obj->get_param("htseq_file_prefix") }, "expected to live" );
    is ($obj->get_param("htseq_file_prefix"), $htseq_file_prefix,
        "get_param(htseq_file_prefix)" );
}

# Test the simple setter methods
{
    throws_ok( sub { $obj->set_param() },
              'MyX::Generic::Undef::Param', "set_param() - caught" );
    lives_ok( sub { $obj->set_param("out_dir", "test") },
            "expected to live" );
    is( $obj->get_param("out_dir"), "test", "set_param(out_dir, test)" );
    
    # now reset the out_dir parma
    $obj->set_param("out_dir", $out_dir);
}

# try reformating
{
    $obj->restructure();
}

