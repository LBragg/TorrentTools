#! /usr/bin/perl -w

use strict;
use IO::File;

#Wonder why this data is required again?
die "Require [file] [outfile] as arguments" if (scalar(@ARGV) != 2);

my $input = IO::File->new($ARGV[0], "r");
my $output = IO::File->new($ARGV[1], "w");

while(my $line = <$input>)
{
	chomp($line);
	my (undef,$ReadID, $callLen, $readpos, $errors, $validFV, $oop) = split(",", $line);

	if($line =~ /ReadID/gi)
	{
		$output->print(join(",",("ReadID", "BasePosition", "Type", "ValidFV", "OOP"))."\n");
		next;
	}

	if($callLen =~ /[A-Za-z]/gi)
	{
		$callLen = 0;
	}

	my @basepositions = $readpos..($readpos + $callLen - 1);
	my @type = ("Correct") x scalar(@basepositions);

	if($errors < 0)
	{
		push @type, ("Deletion") x (-1 * $errors);
		
		my $lastBase = $basepositions[-1];

		if(! defined $lastBase)
		{
			$lastBase = $readpos;
		}

		push @basepositions, ($lastBase) x (-1 * $errors);
	}
	elsif($errors > 0)
	{
		my @indicies = (scalar(@basepositions) - $errors)..scalar(@basepositions - 1);
		@type[@indicies] = ("Insertion")x scalar(@indicies);
	}

	my @id = ($ReadID) x scalar(@basepositions);

	for(my $i = 0; $i < scalar(@id); $i++)
	{
		if(! defined $id[$i] or ! defined $basepositions[$i] or ! defined $type[$i])
		{
			warn  "Not defined id or base position or type\n";
		}
		$output->print(join(",", $id[$i], $basepositions[$i],$type[$i], $validFV, $oop)."\n")
	}
}
