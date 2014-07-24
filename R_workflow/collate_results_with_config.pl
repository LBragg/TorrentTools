#! /usr/bin/perl -w

use strict;
use IO::File;

use lib '/home/bra427/Projects/IonTorrentBenchmarking/profiling_pipeline/IonTorrentPipelineGit/R_workflow';

require IT_config_utility;

die "Require [output dir] <config file> | [ <config file> ] ... as arguments\n" if (scalar(@ARGV) == 0);

my $output_dir = shift @ARGV;
$output_dir = $output_dir."/" if ($output_dir !~ /\/$/gi);

if(!(-e $output_dir and -d $output_dir))
{
	die "Specified output directory does not exist or is not a directory\n";
}


my @configFiles = @ARGV;

my $output_jobs =  0;
my $script_count = 0;
my %errors;
my %subMat;
my %subRateByRLE;
my %polyLocs;
my %allConfigs;
my %allAttributes; #collect the attributes across the files.

foreach my $config (@configFiles)
{
	my %configuration;

	IT_config_utility::loadConfig($config, \%configuration);

        my $fullPath = $configuration{"PARAM"}->{"DIR"};

	$allConfigs{$config} = \%configuration; #check this works

	foreach my $attr (keys %{$configuration{"ATTRIBUTE"}})
	{
		$allAttributes{$attr} = 1;
	}

	#just realised that the error rate will be slightly different, need to include the substitutions as 'compared' homopolymers.
	#what about the called snps? can I be bothered???

	
	my $insDelFP = $fullPath."number_insertion_deletion_errors.csv";
	my $insertionDeletionErrors = IO::File->new($insDelFP, "r") or die "Could not open file for $insDelFP\n";

	 <$insertionDeletionErrors>; #remove title line
	chomp(my $resLine = <$insertionDeletionErrors>);

	my (undef, $insertErrors, $deleteErrors, $correctHPs, $total ,$insertRate, $deleteRate) = split(",", $resLine);

	#ignore those rates for now

	$errors{$config} = {"insertionErr" => $insertErrors, "deletionErr" => $deleteErrors, "correctHomopolymers" => $correctHPs};


	my $subErrFP = $fullPath."substitution_rates.csv";
	my $substitutionErrors = IO::File->new($subErrFP, "r") or die "Could not open file for $subErrFP\n";

	<$substitutionErrors>; #remove header

	chomp($resLine = <$substitutionErrors>);
	my (undef, $numSubs, $covHPnoSNP, $subRate, $snpLocs, $snpinstances, $point51, $chisq, $no_flow_oop, $no_flow_ip, $usable_flow_oop, $usable_flow_ip) = split(",",  $resLine);

	$errors{$config}->{"substitutionErr"} = $numSubs;
	$errors{$config}->{"correctHomopolymerMinusPoly"} = $covHPnoSNP; #not including columns in reference which were considered polymorphisms as opposed to errors
	$errors{$config}->{"substitutionErrRate"} = $subRate;
	$errors{$config}->{"subSNPLocs"} = $snpLocs;
	$errors{$config}->{"subSNPCount"} = $snpinstances;
	$errors{$config}->{"num51SNPerr"} = $point51;
	$errors{$config}->{"chisq"} = $chisq;
	$errors{$config}->{"no_flow_oop"} = $no_flow_oop;
	$errors{$config}->{"no_flow_ip"} = $no_flow_ip;
	$errors{$config}->{"usable_flow_oop"} = $usable_flow_oop;
	$errors{$config}->{"usable_flow_ip"} = $usable_flow_ip;

	my $rawVersusPhased = $fullPath."errors_after_phase_corr.csv";
	my $rvp_fh = IO::File->new($rawVersusPhased, "r");
	<$rvp_fh>; #remove header
	chomp($resLine = <$rvp_fh>);

	my (undef, $lessErrInRaw, $equalRaw, $greaterErrInRaw, $pvalue) = split(",", $resLine);

	$errors{$config}->{'reads_with_fewer_errors_before_phasing'} = $lessErrInRaw;
	$errors{$config}->{'reads_with_no_improvement_after_phasing'} = $equalRaw;
	$errors{$config}->{'reads_with_fewer_errors_after_phasing'} = $greaterErrInRaw;
	$errors{$config}->{'pvalue_raw_versus_phased'} = $pvalue;

	#other stuff to consolidate, for fun, of course?

	my $propSubFP = $fullPath."proportion_of_substitutions.csv";

	my $propSub = IO::File->new($propSubFP, "r") or die "Could not open file for $propSubFP\n";
	chomp(my $header = <$propSub>);
	
	$header =~ s/\"//gi;	

	my @colOrder = split(",", $header);
	
	#ref is row,
	#col is read.
	while(my $line = <$propSub>)
	{
		chomp($line);
	
		my @res = split(/,/, $line);
		my $base = shift @res;
		$base =~ s/\"//gi;
	
		#submat {from} -> {to}
		for(my $i = 0; $i < scalar(@res); $i++)
		{
			#ref -> read
			$subMat{$config}->{$base}->{$colOrder[$i]} = $res[$i]; #double check it is the correct 'from' -> 'to'...
		}
	}


	#this needs up updated too. How are we going to do this?

	my $subRateFP = $fullPath."SNP_error_rate_by_read_pos.csv";
	my $subRateByHP = IO::File->new($subRateFP, "r") or die "Could not open $subRateFP";
	<$subRateByHP>; #remove header;

	while(my $line = <$subRateByHP>)
	{
		chomp($line);	
		my ($undef, $rlepos, $hpcov, $snpfreq, $snprate) =  split(",", $line);
		$rlepos =~ s/\"//gi;

		$subRateByRLE{$config}->{$rlepos} = [$hpcov, $snpfreq, $snprate];
	}
	
	my $plocFP = $fullPath."polymorphism_locations.csv";
	my $plocfh = IO::File->new($plocFP, "r") or die "Could not open $plocFP\n";
	while(my $line = <$plocfh>)
	{
		chomp($line);
		if($line eq "\"x\"")	
		{		
			next;
		}
		
		my (undef, $location) = split(",", $line);
		$polyLocs{$location}->{$config} = 1;
	}
}



#time to output all our data (could have done this while we were going, but i like to print at end.


#this is prepared with all attributes provided for the header.
my @header_base = sort {$a cmp $b} keys %allAttributes;
my $consolidatedErrs = IO::File->new($output_dir."consolidated_error_counts.csv", "w");
my @sortedKeys = ();
my @consolidatedHeader;

push @consolidatedHeader, @header_base;

#header
#this consolidated stuff won't work, as it is tied up in the config. Also need to include all the fields, not just selected fields.
foreach my $config_file (keys %errors)
{
	my @outputLine;
	if(scalar(@sortedKeys) == 0)
	{
		@sortedKeys = sort {$a cmp $b} keys %{$errors{$config_file}};
		push @consolidatedHeader, @sortedKeys;
		$consolidatedErrs->print(join(",",@consolidatedHeader)."\n");
	}

	foreach my $attr (@header_base)
	{
		my $value = (exists $allConfigs{$config_file}->{"ATTRIBUTE"}->{$attr})? $allConfigs{$config_file}->{"ATTRIBUTE"}->{$attr}:  "N/A";
		push @outputLine, $value;
	}

	foreach my $key (@sortedKeys)
	{
		push @outputLine, $errors{$config_file}->{$key};
	}
	$consolidatedErrs->print(join(",", @outputLine)."\n");	
}

my $subPref = IO::File->new($output_dir."consolidated_substitution_rates_between_bases.csv", "w");

my @subPrefHeader;
push @subPrefHeader, @header_base;
push @subPrefHeader, ("RefBase", "ReadBase", "Freq");

$subPref->print(join(",", @subPrefHeader)."\n");
$subPref->print(); #why? 

foreach my $config_file(keys %subMat)
{
	my @attributes;

	foreach my $attr (@header_base)
	{
		my $value = (exists $allConfigs{$config_file}->{"ATTRIBUTE"}->{$attr})? $allConfigs{$config_file}->{"ATTRIBUTE"}->{$attr}:  "N/A";
                push @attributes, $value;
	}

	foreach my $baseOld ("A", "T", "G", "C")
	{
		foreach my $baseNew ("A", "T", "G", "C")
		{
			my @outputLine;
			push @outputLine, @attributes;

			if($baseOld ne $baseNew)
			{
				push @outputLine, ($baseOld, $baseNew, $subMat{$config_file}->{$baseOld}->{$baseNew});
				$subPref->print(join(",", @outputLine)."\n");
			}
		}
	}
}

my $subRate = IO::File->new($output_dir."consolidated_snp_rate_by_base_pos.csv", "w");

my @subRateHeader;

push @subRateHeader,@header_base;
push @subRateHeader, ("BasePos", "HPCov", "SNPFreq", "SNPRate");


foreach my $config_file (keys %subRateByRLE)
{
	my @attributes;
	
	foreach my $attr (@header_base)
        {
                my $value = (exists $allConfigs{$config_file}->{"ATTRIBUTE"}->{$attr})? $allConfigs{$config_file}->{"ATTRIBUTE"}->{$attr}:  "N/A";
                push @attributes, $value;
        }
	
	foreach my $rlepos (sort {int($a) <=> int($b)} keys  %{$subRateByRLE{$config_file}})
	{
		my ($hpCov, $snpFreq, $snpRate) = @{$subRateByRLE{$config_file}->{$rlepos}};
		
		my @outputLine;

		push @outputLine, @attributes;
		push @outputLine, ($rlepos, $hpCov, $snpFreq, $snpRate);
		$subRate->print(join(",", @outputLine)."\n");
	}
}

my $subLoc = IO::File->new($output_dir."consolidated_sub_locations.csv", "w");
my @subLocHeader;
push @subLocHeader, @header_base;
push @subLocHeader, ("referenceLocation");

$subLoc->print(join("\t", @subLocHeader)."\n");

foreach my $location (keys %polyLocs)
{
	foreach my $config_file (keys %{$polyLocs{$location}})
	{
		my @outputLine;
		
		my @attributes;

	        foreach my $attr (@header_base)
	        {
	                my $value = (exists $allConfigs{$config_file}->{"ATTRIBUTE"}->{$attr})? $allConfigs{$config_file}->{"ATTRIBUTE"}->{$attr}:  "N/A";
	                push @attributes, $value;
	        }

		push @outputLine, @attributes;
		push @outputLine, $location;

		$subLoc->print(join(",", @outputLine)."\n");
	}
}
