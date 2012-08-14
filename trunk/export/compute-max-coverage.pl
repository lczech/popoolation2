#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path;
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use MaxCoverage;
use Test;


# Author: Ram Vinay Pandey
# Author: Robert Kofler


# Define the variables
my $input;
my $help=0;
my $test=0;
my $usermaxcoverage;

GetOptions(
    "input=s"	    =>\$input,
    "max-coverage=s"=>\$usermaxcoverage,
    "test"          =>\$test,
    "help"	    =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help;
MaxCoverageTest::runTests() && exit if $test;
pod2usage(-msg=>"A input file has to be provided\n",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Maximum coverage has to be provided",-verbose=>1) unless $usermaxcoverage;
pod2usage(-msg=>"Maximum coverage has to be provided as perentage") unless $usermaxcoverage=~/%$/;

my $maxcoverage=get_max_coverage($input,$usermaxcoverage);

exit(0);

{
    package MaxCoverageTest;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use Test;
    use Pileup;
    use Test::TMaxCoverage;
    use MajorAlleles;
    
    sub runTests
    {
        run_MaxCoverageTests();
        exit;
    }
}



=head1 NAME

compute-max-coverage.pl - compute the maximum coverage thresholds in basepairs given a threshold in percent

=head1 SYNOPSIS

 perl compute-max-coverage.pl --input input.sync  --max-coverage 2%

=head1 OPTIONS

=over 4

=item B<--input>

The input file has to be synchronized pileup file. Mandatory parameter

=item B<--max-coverage>

the maximum coverage threshold in percent that should be converted into a base pair threshold

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

Input is a single tab delimited file which contains a lightwight representation of every pileup file.
Every pileup file represents a population and will be parsed into a list of A-count:T-count:C-count:G-count:N-count:*-count

 2L	5002	G	0:0:0:17:0:0	0:0:0:28:0:0	0:0:0:31:0:0	0:0:0:35:0:0	0:1:0:33:0:0	0:3:0:31:0:0
 2L	5009	A	16:0:0:0:0:0	26:0:0:0:0:0	29:0:1:0:0:0	36:0:0:0:0:0	34:0:0:0:0:0	32:0:1:0:0:0
 2L	5233	G	0:0:5:46:0:0	0:0:0:43:0:0	0:0:0:60:0:0	0:0:3:61:0:0	0:0:0:56:0:0	0:0:0:48:0:0
 
 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: reference genome base
 col 4: population 1
 col 5: population 2
 col n: population n-3
 
 population data are in the form
 A:T:C:G:N:*
 A: count of character A
 T: count of character T
 C: count of character C
 G: count of character G
 N: count of character N
 *: deletion, count of deletion
 
=head2 Output

the maximum coverage threshold in basepairs for every generation:

 use: '--max-coverage 121,76,77'

The output may be provided to any PoPoolation2 scripts requiring a maximum coverage. This has the advantage of reduced computation time (parsing the sync file to determine the max-coverage is time consuming)

=head1 AUTHORS

 Robert Kofler
 Viola Nolte
 Christian Schloetterer

=cut
