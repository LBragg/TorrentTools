#! /usr/bin/perl -w

use strict;
use IO::File;
use File::Basename;
use Getopt::Std;
use vars qw/ %opt /;

sub generateEmpiricalTab($$);
sub generateErrRatesByFlow($$$$$);
sub generateGC($$$$$);
sub generateErrRates($$);
sub generateQualScores($$$$$);
sub generatePoly($$$$$$$$);
sub generateSub($$$$$);
sub generateIndelErrorRate($$);
sub generateConfig($);
sub loadConfig($$);

use constant {
    INPUT_DIR => "DIR",
    PARAM => "PARAM",
    ATTRIBUTE => "ATTRIBUTE",
    LABEL => "LABEL",
    START => "START",
    _END => "END", #clashes with perl reserved word
    REPEAT_FILE => "REPEAT_FILE",
    REF_BASES => "REF_BASES",
    WINDOW_SIZE => "WINDOW_SIZE",
    MIN_WINDOW_SEGMENT => "MIN_WINDOW_SEGMENT",
    MAX_WINDOW_SEGMENT => "MAX_WINDOW_SEGMENT",
    WINDOW_SEGMENT_INCREMENT => "WINDOW_SEGMENT_INCREMENT",
    MIN_READ_LENGTH => "MIN_READ_LENGTH",
    MAX_READ_LENGTH => "MAX_READ_LENGTH",
    MIN_LENGTH_ALIGNED => "MIN_LENGTH_ALIGNED",
    KIT => "KIT",
    CHIP => "CHIP",
    DESC => "DESC",
    SPECIES => "SPECIES"
    };


sub usage()
{
       print "Usage: $0 -c [config file] -g [config file]\n";
       print "-c config file: the config to process\n";
       print "-g config file: output a dummy config to this location\n";
       print "-h: this help message\n";
       exit;
}

sub init()
{
        my $opt_string = 'c:g:h';
        (@ARGV and getopts( "$opt_string", \%opt )) or usage();
        usage() if ($opt{'h'} or !($opt{'c'} and ! $opt{'g'}  or ($opt{'g'} and ! $opt{'c'})));
}
my $path_to_templates = dirname(__FILE__); #assume they are in the same location

MAIN:
{
	init();
	
	if(exists $opt{'g'})
	{
		generateConfig($opt{'g'});
		exit 0;
	}

	my %config;
	loadConfig($opt{'c'}, \%config);

	my $parent_dir = $config{PARAM()}->{INPUT_DIR()}; #this is null?
	$parent_dir = $parent_dir."/" if ($parent_dir !~ /\/$/gi);

	my $baseStart = $config{PARAM()}->{START()};
	my $baseEnd = $config{PARAM()}->{_END()};
	my $repeatMask = $config{PARAM()}->{REPEAT_FILE()};
	my $refBases = $config{PARAM()}->{REF_BASES()};
	my $windowSize = $config{PARAM()}->{WINDOW_SIZE()};
	my $minBase = $config{PARAM()}->{MIN_WINDOW_SEGMENT()};
	my $maxBase = $config{PARAM()}->{MAX_WINDOW_SEGMENT()};
	my $increment = $config{PARAM()}->{WINDOW_SEGMENT_INCREMENT()};
	## addition of parameters for filtering too short or too long reads. TODO -- make it part of the configuration file and used in the poly prep file.
	my $minReadLength = $config{PARAM()}->{MIN_READ_LENGTH()};
	my $maxReadLength = $config{PARAM()}->{MAX_READ_LENGTH()};
	my $minLengthAligned = $config{PARAM()}->{MIN_LENGTH_ALIGNED()};

	my $raw_dir = $parent_dir."rawData/";
	my $screened_dir = $parent_dir."screenedData";

	if(! -e $raw_dir)
	{
		die "Raw directory: $raw_dir, does not exist. Raw input files should go there\n";
	}

	if(! -e $screened_dir)
	{
		system("mkdir $screened_dir");
	}

	#first script to run is
	#prepare_screened_datasets_filter_poly.R (renamed poly script)
	#need to generate the poly script...
	my $polyR_fn = $parent_dir."custom_filter_genuine_poly_before_analysis.R";
	my $custom_poly = IO::File->new($polyR_fn, "w");

	generatePoly($custom_poly, $baseStart, $baseEnd, $repeatMask, $parent_dir, $minReadLength, $maxReadLength, $minLengthAligned);
	$custom_poly->close();

	#system("R CMD BATCH --no-restore --no-save $polyR_fn ".$polyR_fn."out");

	#okay, the screened data should now be available in screenedData.

#	if(! -e $screened_dir."/all_base_obs_filtered.RData")#whatever that file is called
#	{
#		die "Did not successfully produce the screened output\n";
#	}

	#successful, so move to script 2.
	#might take care of the expansion of the flow data for now
	my $flowDataPath = $parent_dir."flow_data.csv";
	my $flowDataExpandedPath = $parent_dir."flow_data_expanded.csv";

	#system("/srv/whitlam/projects1/data_modelling/IonTorrent/CleanedWorkflow/expand_flow_data.pl $flowDataPath $flowDataExpandedPath");

#	if(! -e $flowDataExpandedPath)
#	{
#		die "Failed to expand the flow data file\n";
#	}

	#think it is safe to analyse substitutions now...
	my $subsR_fn = $parent_dir."custom_analyse_substitutions.R";
	my $custom_subs = IO::File->new($subsR_fn, "w");
	generateSub($custom_subs, $baseStart, $baseEnd, $repeatMask, $parent_dir);
	$custom_subs->close();
	#system("R CMD BATCH --no-restore --no-save $subsR_fn $subsR_fn"."out");

	#output the indel rate...
	my $indelErrorRateR_fn = $parent_dir."custom_output_indel_rate.R";
	my $custom_indelErrorRate = IO::File->new($indelErrorRateR_fn, "w");
	generateIndelErrorRate($custom_indelErrorRate, $parent_dir);
	$custom_indelErrorRate->close();

	#system("R CMD BATCH --no-restore --no-save $indelErrorRateR_fn $indelErrorRateR_fn"."out");


	#quality scores, requires many of those variables again...
	my $qualScores_fn = $parent_dir."custom_analyse_quality_scores.R";
	my $custom_qualscores = IO::File->new($qualScores_fn, "w");
	generateQualScores($custom_qualscores, $baseStart, $baseEnd, $repeatMask, $parent_dir);
	$custom_qualscores->close();
	
	#system("R CMD BATCH  --no-restore --no-save $qualScores_fn $qualScores_fn"."out");				
	#finished, should probably check that the output files are being created... wonder if R will return a fail value?
	my $error_rates_fn = $parent_dir."custom_output_error_rates.R";
	my $custom_error_rates = IO::File->new($error_rates_fn, "w");
		
	generateErrRates($custom_error_rates, $parent_dir);
	$custom_error_rates->close();

	#system("R CMD BATCH  --no-restore --no-save $error_rates_fn $error_rates_fn"."out");

	#next script...
	my $gc_fn = $parent_dir."custom_output_gc.R";
	my $custom_gc = IO::File->new($gc_fn, "w");
	
	generateGC($custom_gc, $parent_dir, $repeatMask, $refBases, $windowSize);
	$custom_gc->close();

	#system("R CMD BATCH --no-restore --no-save $gc_fn $gc_fn"."out");

	#next!
	my $errRateByFlow_fn = $parent_dir."custom_errRateByFlow.R";
	my $custom_errRateByFlow = IO::File->new($errRateByFlow_fn, "w");

	generateErrRatesByFlow($custom_errRateByFlow, $parent_dir, $minBase, $maxBase, $increment);
	$custom_errRateByFlow->close();
	#system("R CMD BATCH --no-restore --no-save $errRateByFlow_fn $errRateByFlow_fn"."out");


	#assume it all went well for now...
	my $empirical_tab_fn = $parent_dir."custom_tabling_emp.R";
	my $custom_empirical = IO::File->new($empirical_tab_fn, "w");
	generateEmpiricalTab($custom_empirical, $parent_dir);
	$custom_empirical->close();

	#system("R CMD BATCH --no-restore --no-save $empirical_tab_fn $empirical_tab_fn"."out");


	#that completes the R scripts for now, excluding the modelling of flow values.. everything that remains is perl
	#generally, the perl scripts relate to collation. Looking for perl scripts that are run on individual datasets
	#below are used for trimming comparison. Probably not necessary.
#       system("perl /srv/whitlam/home/projects/DATA_MODELLING/IonTorrent/FixedInsertionsOnly/execute_laurens_trimming_mechanism.pl flow_data.csv flow_data_expanded.csv $qualityClipFile $maxFlow trimmedFlowData_corr_bounds.csv
      #system("perl /srv/whitlam/home/projects/DATA_MODELLING/IonTorrent/FixedInsertionsOnly/generate_simplified_data_for_R.pl trimmedFlowData.csv &");

	
	#this can be run, but is later amalgamated.

	#system("perl /srv/whitlam/projects1/data_modelling/IonTorrent/CleanedWorkflow/process_error_rates_by_flow.pl ".$parent_dir."error_rate_by_flowpos_reflen.csv ".$parent_dir."error_rate_for_zero.csv $parent_dir");


	#all the collate scripts should tap into the configuration files...

	


}

sub generateEmpiricalTab($$)
{
	my ($fh, $parent_dir) = @_;

	$fh->print("parentDir <- \"$parent_dir\"\n");

	my $template_emp =  IO::File->new($path_to_templates."/tabling_empirical.R", "r");

        while(my $line = <$template_emp>)
        {
                $fh->print($line);
        }

}

sub generateErrRatesByFlow($$$$$)
{
	my ($fh, $parent_dir, $minBase, $maxBase, $increment) = @_;
 
	my @lims;

	for(my $i = $minBase; $i <= $maxBase; $i+= $increment)
	{
		push @lims, $i;
	}
	
	$fh->print("parentDir <- \"$parent_dir\"\n");
	$fh->print("maxBounds <- c(".(join(",", @lims)).")\n");

	my $template_RBF = IO::File->new($path_to_templates."/compute_error_rates_by_flowpos.R", "r");
	
        while(my $line = <$template_RBF>)
        {
                $fh->print($line);
        }
}

sub generateGC($$$$$)
{
	my ($fh, $parent_dir, $repeatLocations, $refBases, $windowSize) = @_;

	#Reminder
	#genomeFileSul <- "~/IT/SulfobusTokodaii_bases.gc"
	#genomeFileBac <- "~/IT/BacillusAmyloliquefaciens_bases.gc"
	#might need to look for the script that generated these.
        $fh->print("parentDir <- \"$parent_dir\"\n");
       	$fh->print("repeatLocations <- read.table(\"$repeatLocations\", header=TRUE, sep=\"\\t\")\n");
	$fh->print("genomeFile <- \"$refBases\"\n");
	$fh->print("windowSize <- $windowSize\n");

	my $template_gc = IO::File->new($path_to_templates."/gcScript.R", "r");

	while(my $line = <$template_gc>)
        {
                $fh->print($line);
        }
}


sub generateErrRates($$)
{
        my ($fh, $parent_dir) = @_;
	
	$fh->print("parentDir <- \"$parent_dir\"\n");

	my $template_errRate = IO::File->new($path_to_templates."/output_error_rates.R", "r");

        while(my $line = <$template_errRate>)
        {
                $fh->print($line);
        }
}


sub generateQualScores($$$$$)
{
	my ($fh, $start, $truncate, $repeat_mask, $parent_dir) = @_;

        $fh->print("startPos <- $start\n");
        $fh->print("truncatePos <- $truncate\n"); #change this for the different kits
        $fh->print("repeatLocations <- read.table(\"$repeat_mask\", header=TRUE, sep=\"\\t\")\n");
        $fh->print("parentDir <- \"$parent_dir\"\n");

	my $template_qual = IO::File->new($path_to_templates."/analyse_quality_scores.R", "r");

        while(my $line = <$template_qual>)
        {
                $fh->print($line);
        }
}


#when should I append information about the dataset?
sub generatePoly($$$$$$$$)
{
	my ($fh, $start, $truncate, $repeat_mask, $parent_dir, $minReadLength, $maxReadLength, $minLengthAligned) = @_;

        $fh->print("startPos <- $start\n");
	$fh->print("truncatePos <- $truncate\n"); #change this for the different kits
	$fh->print("repeatLocations <- read.table(\"$repeat_mask\", header=TRUE, sep=\"\\t\")\n");
	$fh->print("parentDir <- \"$parent_dir\"\n");
        $fh->print("minrl <- \"$minReadLength\"\n");
	$fh->print("maxrl <- \"$maxReadLength\"\n");
	$fh->print("minLAlign <- \"$minLengthAligned\"\n");
	
	#as a reminder, the previous repeat masks are stored here:
	#repeatLocationsSul <- "/srv/whitlam/home/projects/DATA_MODELLING/IonTorrent/RAnalysisIT/Sulfolobus_repeat_locations.txt"
	#repeatLocationsBac <- "/srv/whitlam/home/projects/DATA_MODELLING/IonTorrent/RAnalysisIT/Bacillus_repeat_locations.txt"

	my $template_poly = IO::File->new($path_to_templates."/filter_genuine_poly_before_analysis.R", "r");

	while(my $line = <$template_poly>)
	{
		$fh->print($line);
	}
}

sub generateSub($$$$$)
{
	my ($fh, $start, $truncate, $repeat_mask, $parent_dir) = @_;

        $fh->print("startPos <- $start\n");
        $fh->print("truncatePos <- $truncate\n"); #change this for the different kits
        $fh->print("repeatLocations <- read.table(\"$repeat_mask\", header=TRUE, sep=\"\\t\")\n");
        $fh->print("parentDir <- \"$parent_dir\"\n");

	my $template_sub =  IO::File->new($path_to_templates."/analyse_substitutions.R", "r");

	while(my $line = <$template_sub>)
        {
                $fh->print($line);
        }
}

sub generateIndelErrorRate($$)
{
	my ($fh, $parent_dir) = @_;
	$fh->print("parentDir <- \"$parent_dir\"\n");
	
	my $template_indel_err = IO::File->new($path_to_templates."/output_indel_rate.R", "r");
	
	while(my $line = <$template_indel_err>)
        {
                $fh->print($line);
        }
}

sub generateConfig($)
{
	my $file = $_[0];
	my $newConfig = IO::File->new($file, "w");

	$newConfig->print(PARAM()."\t".INPUT_DIR().":VALUE\n");
	$newConfig->print(PARAM()."\t".LABEL().":VALUE\n");
	$newConfig->print(PARAM()."\t".START().":10\n");
	$newConfig->print(PARAM()."\t"._END().":VALUE\n");
	$newConfig->print(PARAM()."\t".REPEAT_FILE().":VALUE\n");
	$newConfig->print(PARAM()."\t".REF_BASES().":VALUE\n");
	$newConfig->print(PARAM()."\t".WINDOW_SIZE().":VALUE\n");
       	$newConfig->print(PARAM()."\t".MIN_WINDOW_SEGMENT().":200\n");
       	$newConfig->print(PARAM()."\t".MAX_WINDOW_SEGMENT().":300\n");
       	$newConfig->print(PARAM()."\t".WINDOW_SEGMENT_INCREMENT().":100\n");
	$newConfig->print(PARAM()."\t".MIN_READ_LENGTH().":30\n");
 	$newConfig->print(PARAM()."\t".MAX_READ_LENGTH().":500\n");
	$newConfig->print(PARAM()."\t".MIN_LENGTH_ALIGNED().":30\n");
	$newConfig->print(ATTRIBUTE()."\t".KIT().":VALUE\n");
	$newConfig->print(ATTRIBUTE()."\t".CHIP().":VALUE\n");
	$newConfig->print(ATTRIBUTE()."\t".DESC().":VALUE\n");
	$newConfig->print(ATTRIBUTE()."\t".SPECIES().":VALUE\n");
	
	$newConfig->close();
}

sub loadConfig($$)
{
	my ($file_name, $config) = @_;

	my $file = IO::File->new($file_name, "r");

	while(my $line = <$file>)
	{
		chomp($line);
		
		my ($name, $value) = split(/\t/, $line);
		my ($att, $att_val) = split(":", $value);

		print "Storing '$name' '$att' '$att_val'\n";

		$config->{$name}->{$att} = $att_val;
	}
}

