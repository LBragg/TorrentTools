#! /usr/bin/perl -w 

use strict;
use IO::File;

sub processSubDir($$$);

#This is really a redundant script, the read length stuff should be incorporated into the join script at some point.
die "Require [main output directory]" if scalar(@ARGV != 1);

my $dirToProcess = $ARGV[0];
$dirToProcess = $dirToProcess."/" if ($dirToProcess !~ /\/$/gi);

#original tables that are normally loaded
my $output_dir = $dirToProcess."/rawData/";

if (! -e $output_dir)
{
        system("mkdir $output_dir");
}

opendir DIR, $dirToProcess;
my @subDirToProcess;

my $output = IO::File->new($output_dir."/read_lengths.tsv", "w");

while(my $entry = readdir DIR)
{
#       print "Entry: $entry\n";
        my $fullDirPath = $dirToProcess.$entry;
        if($entry =~ /read_dir_([0-9]+)_for_/gi)
        {
                print "Processing sub dir $fullDirPath\n";
		processSubDir($output, $fullDirPath, $1);
        }
}


sub processSubDir($$$)
{
	my ($output, $path_to_read_dir, $partition) = @_;
	my ($max_read_pos, $max_rle_pos) = (0,0);

	my $input = IO::File->new($path_to_read_dir."/reads.flow", "r");
	my $prev_id = undef;

	while(my $line = <$input>)
	{
		chomp($line);
		next if($line =~ /FLOW_VAL/gi);
                my ($read_id, $flow_pos, $rle_pos, $read_pos, $nucleotide, $flow_val, $usable, $oop, $called_len) = split(/\t/, $line);
		if(defined $prev_id and $prev_id ne $read_id)
		{
			$output->print(join("\t", $prev_id, $max_read_pos, $max_rle_pos)."\n");
			$max_read_pos = 0;
			$max_rle_pos = 0;
			$prev_id = $read_id;
		}
		
		if(! defined $prev_id)
		{		
			$prev_id = $read_id;
		}	

		if($max_read_pos < $read_pos)
		{
			$max_read_pos = $read_pos;
		}
	
		if($max_rle_pos < $rle_pos)
		{
			$max_rle_pos = $rle_pos;
		}
	}

	if(defined $prev_id)
	{
		$output->print(join("\t", $prev_id, $max_read_pos, $max_rle_pos)."\n");
	}
	
}



