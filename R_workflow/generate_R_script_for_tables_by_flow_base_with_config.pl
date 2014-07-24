#! /usr/bin/perl -w

use strict;
use IO::File;

use lib '/home/bra427/Projects/IonTorrentBenchmarking/profiling_pipeline/IonTorrentPipelineGit/R_workflow/';
require IT_config_utility;

die "Require [output dir] <config file> | [ <config file> ] ... as arguments\n" if (scalar(@ARGV) == 0);

my $output_dir = shift @ARGV;
$output_dir = $output_dir."/" if ($output_dir !~ /\/$/gi);

my @configFiles = @ARGV;

my %allConfigs;
my %allAttributes; #collect the attributes across the files.

## you will probably need a phasing thing 
## you know that the num errors per read and qual scores per read now have phasing information
my $flowFile_str = "error_rate_by_flowpos_phased.csv";
my $baseFile_str = "error_rate_by_base_position.csv";
my $errorsPerRead_str = "num_errors_per_read_with_phase_fv.csv";
my $qualScores_str  = "quality_scores_by_phase_fv.csv";

## a couple of these are by flow value or phase. 
## these need to be new types of analysis.
my $rScriptFlowFile = IO::File->new($output_dir."plot_flow_errors_for_all.R", "w");
my $rScriptBaseFile = IO::File->new($output_dir."plot_base_errors_for_all.R", "w");
my $rScriptErrorsFile = IO::File->new($output_dir."plot_errors_per_read.R", "w");
my $rScriptQualFile =  IO::File->new($output_dir."plot_qualities_for_read.R", "w");

my @inputOutput = ([$flowFile_str, $rScriptFlowFile], [$baseFile_str, $rScriptBaseFile], [$errorsPerRead_str, $rScriptErrorsFile], [$qualScores_str, $rScriptQualFile]);

foreach my $rscript ($rScriptFlowFile,$rScriptBaseFile,$rScriptErrorsFile, $rScriptQualFile)
{
	$rscript->print("res <- c()\n"); #set up that res object.
}

foreach my $config (@configFiles)
{
	my %configuration;
	IT_config_utility::loadConfig($config, \%configuration);
	$allConfigs{$config} = \%configuration; #check this works
        foreach my $attr (keys %{$configuration{"ATTRIBUTE"}})
        {
                $allAttributes{$attr} = 1;
        }
}

foreach my $config (@configFiles)
{
        my %configuration = %{$allConfigs{$config}};
        my $fullPath = $configuration{"PARAM"}->{"DIR"};
	foreach my $arry(@inputOutput)
	{
		my ($input_file, $rscript) = @{$arry};
		#might make this conditional on Rscriptflowfile, but see how we go
		my $full_input_path = $fullPath.$input_file;
		if($rscript == $rScriptFlowFile)
		{			
			$rscript->print("tmp = read.table(\"$full_input_path\", sep=\" \", header=TRUE)\n"); 
		}
		elsif($rscript == $rScriptErrorsFile) #errors per read.
		{
			$rscript->print("tmp = read.table(\"$full_input_path\", sep=\" \", header=TRUE)\n");		
		}
		elsif($rscript == $rScriptBaseFile)
		{
			$rscript->print("tmp <- read.table(\"$full_input_path\", sep=\" \", header=TRUE)\n"); #works.
		}
		else ###must be qual file? 
		{
			$rscript->print("tmp <- read.table(\"$full_input_path\", sep=\",\", header=TRUE)\n");
		}

		foreach my $attr (sort {$a cmp $b} keys %allAttributes)
		{
			my $value = (exists $configuration{"ATTRIBUTE"}->{$attr})? $configuration{"ATTRIBUTE"}->{$attr} : "N/A";
		
			$rscript->print("tmp\$$attr <- \"$value\"\n");
		}

		$rscript->print("res <- rbind(res, tmp)\n");		
		$rscript->print("rm(tmp)\n");
	}
}

foreach my $pair ([$rScriptFlowFile, "error_rate_by_flow_values.csv"], [$rScriptBaseFile, "error_rate_by_base_value.csv"], [$rScriptErrorsFile, "error_rate_by_error_distribution.csv"], [$rScriptQualFile, "qual_accuracy_across_datasets.csv"] )
{
	my ($rscript, $output_file) = @{$pair};
 	$rscript->print("write.csv(res, file=\"$output_file\")");	
}
