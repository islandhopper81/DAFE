use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'DAFE',
    license             => 'artistic2',
    dist_author         => 'Scott Yourstone and Nick Colaianni <scott.yourstone81@gmail.com>',
    dist_version   => '0.0.1',
    requires => {
        'Test::More' => 0,
        'version'    => 0,
    },
    add_to_cleanup      => [ 'DAFE-*' ],
    recursive_test_files => 1,
);

$builder->create_build_script();
