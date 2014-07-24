#! /usr/bin/perl -w

use strict;
use IO::File;

use lib '/home/bra427/Projects/IonTorrentBenchmarking/profiling_pipeline/IonTorrentPipelineGit/R_workflow';
require IT_config_utility;

die "Require [output dir] <config file> | [ <config file> ] ... as arguments\n" if (scalar(@ARGV) == 0);

my $output_dir = shift @ARGV;
$output_dir = $output_dir."/" if ($output_dir !~ /\/$/gi);

my @configFiles = @ARGV;

my %allConfigs;
my %allAttributes; #collect the attributes across the files.

#input files

#actually need to check the config for the number of input files??

my $accuracy_prefix = "error_rate_by_flowpos_reflen_first"; #.csv
my $accuracy_zero_prefix = "error_rate_for_zero_first"; #.csv


my $rscript = IO::File->new($output_dir."aggregate_accuracy_by_hp_length.R", "w");

$rscript->print("res <- c()\n");

my @phase_str = ("inphase", "oop");

foreach my $config (@configFiles)
{
 	my %configuration;
        IT_config_utility::loadConfig($config, \%configuration);
	
	#input path
        my $fullPath = $configuration{"PARAM"}->{"DIR"};

	$fullPath = $fullPath."/" if ($fullPath !~ /\/$/gi);

        $allConfigs{$config} = \%configuration; #check this works

        foreach my $attr (keys %{$configuration{"ATTRIBUTE"}})
        {
                $allAttributes{$attr} = 1;
        }

	my $min = $configuration{"PARAM"}->{"MIN_WINDOW_SEGMENT"};
	my $max = $configuration{"PARAM"}->{"MAX_WINDOW_SEGMENT"};
	my $interval = $configuration{"PARAM"}->{"WINDOW_SEGMENT_INCREMENT"};

	my %files;
	my %filesZero;
	

	for(my $i = $min; $i <= $max; $i+= $interval)
	{
		foreach my $phase (@phase_str)
		{	
			$files{$i}->{$phase} = $fullPath.$accuracy_prefix.$i.$phase.".csv";
			$filesZero{$i}->{$phase} = $fullPath.$accuracy_zero_prefix.$i.$phase.".csv";
		}
	}

	my $counter = 0;

	foreach my $endPoint (keys %files)
	{
		foreach my $phase (keys %{$files{$endPoint}})
		{
			my $file = $files{$endPoint}->{$phase};
			my $tmpName = "my".$endPoint;
			
			$rscript->print($tmpName." <- read.table(\"$file\", header=TRUE)\n");
			$rscript->print($tmpName."\$bases <- \"1 - $endPoint\"\n");
			$rscript->print($tmpName."\$phase <- \"$phase\"\n");

			if($counter == 0)
			{
				$rscript->print("my.df <- data.frame(\"RefLen\" = $tmpName\$RefLen, \"FlowPos\" = $tmpName\$FlowPos, \"count.above\" = $tmpName\$count.above, \"count.at\" = $tmpName\$count.at, \"count.below\" = $tmpName\$count.below, \"bases\" = $tmpName\$bases, \"phase\" = $tmpName\$phase)\n");
			}
			else
			{
				$rscript->print("my.df <- rbind(my.df, $tmpName)\n");
			}
        		$rscript->print("rm($tmpName)\n");
			$counter++;
		}
	}
	
	my $counterZero = 0;

	foreach my $endPoint (keys %filesZero)
	{
		foreach my $phase (keys %{$files{$endPoint}})
		{
			my $file = $filesZero{$endPoint}->{$phase};	
			my $tmpName = "my".$endPoint."Zero";
	        	$rscript->print("$tmpName <- read.table(\"$file\", header=TRUE)\n");
			$rscript->print("$tmpName = subset($tmpName, select=-c(N, pos_calls)) \n");
			$rscript->print("$tmpName\$RefLen = 0\n");
	        	$rscript->print("$tmpName\$count.below <- 0\n");
	        	$rscript->print("$tmpName\$bases <- \"1 - $endPoint\"\n");
			$rscript->print("$tmpName\$phase = \"$phase\"\n");
	        	$rscript->print("my.df <- rbind(my.df, $tmpName)\n");
			$rscript->print("rm($tmpName)\n");
		}
	}

	foreach my $attr (keys %{$configuration{"ATTRIBUTE"}})
	{
		my $value = $configuration{"ATTRIBUTE"}->{$attr};
		$rscript->print("my.df\$$attr <- \"$value\"\n");
	}
	$rscript->print("res <- rbind(res, my.df)\n");		
	$rscript->print("rm(my.df)\n");
}

$rscript->print("save(res, file=\"$output_dir"."combined_accuracy_of_hp_firstXbases.RData\")\n");	
