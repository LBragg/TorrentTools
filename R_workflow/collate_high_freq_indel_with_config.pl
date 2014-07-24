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
my $gagErrOut = IO::File->new($output_dir."high_occ_indels.csv", "w");

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

my @sorted_attr = sort {$a cmp $b} keys %allAttributes;

# are these enough headers.
my @other_header = ("READ_ID", "FLOW_VAL", "FLOW_POS", "RLE_POS", "READ_POS", "REF_RLE_POS", "REF_REAL_POS", "REF_LEN", "REF_BASE", "STRAND", "VALID_FV", "OOP", "CALL_LEN", "BOPS");


my @full_header;

push @full_header, @sorted_attr;
push @full_header, @other_header;

$gagErrOut->print(join(",", @full_header)."\n");

foreach my $config (@configFiles)
{
	my %configuration = %{$allConfigs{$config}};

	#just realised that the error rate will be slightly different, need to include the substitutions as 'compared' homopolymers.
	#what about the called snps? can I be bothered???
	my $fullPath = $configuration{"PARAM"}->{"DIR"};

	my @attribute_line;

	foreach my $attr (@sorted_attr)
	{
		my $value = (exists $configuration{"ATTRIBUTE"}->{$attr})? $configuration{"ATTRIBUTE"}->{$attr} : "N/A";
	
		push @attribute_line, $value;
	}

	my $errsInReference = IO::File->new($fullPath."/possible_indel_errors_in_ref.txt", "r");
	 <$errsInReference>; #remove title line


	#"ReadID","FlowVal","FlowPos","RLEPos","ReadPos","RefRLEPos","RefRealPos","RefLen","RefBase","Strand","ValidFV","OOP","CallLen","Error","Errors","baseOnPosStrand"
	while(my $resLine = <$errsInReference>)
	{
		chomp($resLine);
		my ($undef,$read_id, $flow_val,$flow_pos, $rle_pos, $read_pos, $ref_rle_pos, $ref_real_pos, $reflen, $refbase, $strand, $validfv, $oop, $callLen, $error, $errors, $bops) = split(",", $resLine);

		my @valuesOfInt = ($read_id, $flow_val,$flow_pos, $rle_pos, $read_pos, $ref_rle_pos, $ref_real_pos, $reflen, $refbase, $strand, $validfv, $oop, $callLen,$bops);
		my @line;

		push @line, @attribute_line;
		push @line, @valuesOfInt;

	$gagErrOut->print(join(",", @line)."\n");
	}
}
