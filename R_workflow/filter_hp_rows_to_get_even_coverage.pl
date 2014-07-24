#! /usr/bin/perl -w

use strict;
use IO::File;

die "Require input directory [dir] as input" if (scalar(@ARGV) != 1);
my $parent_dir = $ARGV[0];
die "Could not open parent directory $parent_dir: does not exist" if(! -e $parent_dir);

my $unfiltered = $parent_dir."/rawData/flowvaluestolengths_for_overlapping_rle.txt";
my $un_fh = IO::File->new($unfiltered, "r");

my %genomic_locations;
my $total_rows = 0;

while(my $line = <$un_fh>)
{
	chomp($line);
	my ($ReadID, $FlowVal, $FlowPos, $RLEPos, $ReadPos, $RefRLEPos, $RefRealPos, $RefLen, $RefBase, $Strand, $ValidFV, $OOP, $CallLen) = split(/\t/, $line);

	$genomic_locations{$RefRealPos} = (exists $genomic_locations{$RefRealPos})? $genomic_locations{$RefRealPos} + 1 : 1;
	$total_rows++;	
}

my $max_rows_per_coordinate = $total_rows / int((scalar(keys %genomic_locations)));

$max_rows_per_coordinate = 30 if($max_rows_per_coordinate > 30);

print "Going to try and get  $max_rows_per_coordinate\n";

my $filtered = IO::File->new($parent_dir."/rawData/flowvaluestolengths_for_overlapping_rle_even_coverage.txt", "w"); ## coverage bias should not look at this even coverage data

## re-open unfiltered.

$un_fh = IO::File->new($unfiltered, "r");

my %current_cover;

while(my $line = <$un_fh>)
{
	chomp($line);

	my ($ReadID, $FlowVal, $FlowPos, $RLEPos, $ReadPos, $RefRLEPos, $RefRealPos, $RefLen, $RefBase, $Strand, $ValidFV, $OOP, $CallLen) = split(/\t/, $line);

	if(! exists $current_cover{$RefRealPos} or $current_cover{$RefRealPos} <= $max_rows_per_coordinate)
	{
		$filtered->write($line."\n");
		$current_cover{$RefRealPos} = (exists $current_cover{$RefRealPos})? $current_cover{$RefRealPos} + 1 : 1;
	}
}

