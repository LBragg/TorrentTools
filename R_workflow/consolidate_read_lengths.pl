#! /usr/bin/perl -w 
use strict;
use IO::File;

use lib '/home/bra427/Projects/IonTorrentBenchmarking/profiling_pipeline/IonTorrentPipelineGit/R_workflow/';

require IT_config_utility;
my %paramsToScripts;
sub processFile($$$$);

die "Require [output dir] <config file> | [ <config file> ] ... as arguments\n" if (scalar(@ARGV) == 0);

my $output_dir = shift @ARGV;
$output_dir = $output_dir."/" if ($output_dir !~ /\/$/gi);

my @configFiles = @ARGV;

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
}

my $outfilesizes = IO::File->new($output_dir."consolidated_read_lengths.csv", "w");
my @sorted_attr = sort {$a cmp $b} keys %allAttributes;
my @header;

push @header, @sorted_attr;
push @header, ("READ_LENGTH", "FREQ", "TOTAL");

$outfilesizes->print(join(",", @header)."\n");

foreach my $config (@configFiles)
{
	my %configuration = %{$allConfigs{$config}};
	my $dirname = $configuration{"PARAM"}->{"DIR"};
	opendir DIR, $dirname or die "Could not open dir $!";

	my $seen = 0;
	my $sizeGC = undef;
	my $bam = undef;

	my $label = $configuration{"PARAM"}->{"LABEL"};

	print "Processing dir : ".$dirname."\n";

	my $rawlibbam = "rawlib.basecaller.bam";

	print "Looking for file $rawlibbam\n";

	foreach my $file (readdir DIR)
	{
		print "File: $file\n";
		if($file =~ /\.read_lengths/gi)
		{
			my $full_path = $dirname."/".$file;
			$sizeGC = $full_path;
		}
		elsif($file eq  $rawlibbam)
		{
			$bam = $dirname."/".$file;
		}
	}

	my $dataName = $configuration{"PARAM"}->{"LABEL"};

	if(! defined $sizeGC or -s $sizeGC < 100)
	{
		my $gcFile = $dirname.$dataName.".read_lengths";

		system("python /home/bra427/Scripts/calculate_sequence_size_from_BAM.py -b $bam -o $gcFile"); 
	
		if(! -e $gcFile or -s $gcFile <= 68)
		{
			die "Could not source gc file, nor create it\n";
		}
	}

	if(! defined $sizeGC)
	{
		print "Not seen $dirname\n";
	}
	else
	{
		processFile($sizeGC, \%configuration, \@sorted_attr, $outfilesizes);
	}
}


sub processFile($$$$)
{
	my ($gc, $config, $sorted_attr_arry, $outfile) = @_;
	my @sorted_attr = @{$sorted_attr_arry};

	my %numReadsForLength;
	my $inputFile = IO::File->new($gc, "r");
	my $total = 0;

	while(my $line = <$inputFile>)
	{
		chomp($line);

		my ($readid, $length) = split(/\t/, $line);
		my $shortID = ($readid =~ /^([^\s]+)/gi)[0];

		$numReadsForLength{$length} = 0 if (! exists $numReadsForLength{$length});
		$numReadsForLength{$length} = $numReadsForLength{$length} + 1;
		$total++;
	}

	my @outputPrefix;
	foreach my $attr (@sorted_attr)
	{
		my $value = (exists $config->{"ATTRIBUTE"}->{$attr})?  $config->{"ATTRIBUTE"}->{$attr} : "N/A";
		push @outputPrefix, $value;
	}
	
	foreach my $readLength (sort {$a <=> $b} keys %numReadsForLength)
	{
		my @outputLine;
		push @outputLine, @outputPrefix;
		push @outputLine, ($readLength, $numReadsForLength{$readLength}, $total);
	
		$outfile->print(join(",", @outputLine)."\n");
	}
}
