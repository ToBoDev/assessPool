{
package Test;
use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT     =qw(is not_exists ok is_smaller print_teststat);

my $testcounter;
my $failcounter;

sub print_teststat
{
    my $test=$testcounter-1;
    my $failed=$failcounter;
    my $passed=$test-$failed;
    if($failed)
    {
        print "Unit tests end status: FAILED\n";
    }
    else
    {
        print "Unit tests end status: OK\n"
    }
    
    print "Run $test Unit-tests; Passed Unit-tests $passed\n"
    
}

sub is
    {

        my $a=shift;
        my $b=shift;
        my $msg=shift;
        _initialize();
        
        if($a eq $b)
        {
            print "$testcounter - OK: $a = $b; $msg\n";
        }
        else
        {
            print "$testcounter - FAILED: $a = $b; $msg\n";
            $failcounter++;
        }
        $testcounter++;
    }
    
sub _initialize
{
        $testcounter=1 unless $testcounter;
        $failcounter=0 unless $failcounter;
}
    
sub is_smaller
    {

        my $a=shift;
        my $b=shift;
        my $msg=shift;
        _initialize();
        
        if($a < $b)
        {
            print "$testcounter - OK: $a < $b; $msg\n";
        }
        else
        {
            print "$testcounter - FAILED: $a < $b; $msg\n";
            $failcounter++;
        }
        $testcounter++;
    }
    
    sub not_exists
    {
        my $a=shift;
        my $msg=shift;
        _initialize();
        
        
        if($a)
        {
            print "$testcounter - FAILED: $msg\n";
            $failcounter++;
        }
        else
        {
            print "$testcounter - OK: $msg\n";

        }
        $testcounter++;
    }

    sub ok
    {
        my $a=shift;
        my $msg=shift;
        _initialize();
        
        if($a)
        {
            print "$testcounter - OK: $msg\n";
        }
        else
        {
            print "$testcounter - FAILED: $msg\n";
            $failcounter++;
        }
        $testcounter++;
    }
    



}

1;