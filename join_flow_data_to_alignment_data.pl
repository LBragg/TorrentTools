#! /usr/bin/perl -w

use strict;
use IO::File;

# This script joins the alignment information from parsing segemehl to the flow and qual data for the reads.
# Uses a fair bit of RAM.

die "Require [main output directory]" if scalar(@ARGV != 1);

sub processAll($$$$$$$$$);
sub loadFlows($$$$);
sub processReference($$);
sub loadSubBases($$);
sub loadInsBases($$);
sub loadCorrRLE($$$);

my $dirToProcess = $ARGV[0];

$dirToProcess = $dirToProcess."/" if ($dirToProcess !~ /\/$/gi);

my %refbase2rle;
my %corrBases;
my %substitutions;
my %insertions; 
my %deletions;
my $base2qualfile;
my $flow2basefile;

my $output_dir = $dirToProcess."/rawData/";

if (! -e $output_dir)
{
	system("mkdir $output_dir");
}

opendir DIR, $dirToProcess;
my @subDirToProcess;

while(my $entry = readdir DIR)
{
	my $fullDirPath = $dirToProcess.$entry;
	if($entry =~ /read_dir_([0-9]+)_for_/gi)
	{
		print "Processing sub dir $fullDirPath\n";
		push @subDirToProcess, [$fullDirPath, $1];
	}
	elsif($entry =~ /^reference_rle_for_/gi)
	{
		print "Processing reference dir $fullDirPath\n";
		processReference($fullDirPath, \%refbase2rle);
	}
	elsif($entry =~ /reads_qual\.sql/gi)
	{
		print "Found base2qual $fullDirPath\n";
		$base2qualfile = $fullDirPath;
	}
}

#output files might need to be set up here..
my $output_base2qual = $output_dir."base2qual.txt";

#substitution differences
my $outfile_sic = $output_dir."substitution_information_complete.txt";

#aligned bases (in RLE format)
my $outfile_fvtl = $output_dir."flowvaluestolengths_for_overlapping_rle.txt";

# Read starting positions in the genome 
my $outfile_rsp = $output_dir."read_starting_positions_genome.txt";

#coverage RLE
my $outfile_cr = $output_dir."coverageRLE.tsv";

#insertion differences
my $outfile_fvi = $output_dir."flowvalues_insertions_full.txt";

## not used.
my $outfile_c = $output_dir."reads_incomplete_extension.txt";

print "Creating output files...\n";
print "$output_base2qual\n$outfile_sic\n$outfile_fvtl\n$outfile_rsp\n$outfile_cr\n$outfile_fvi\n$outfile_c\n";

my $outfile_sub_in_complete = IO::File->new($outfile_sic, "w");
my $outfile_flowvalues_to_len = IO::File->new($outfile_fvtl, "w");
my $outfile_read_start_pos = IO::File->new($outfile_rsp, "w");
my $outfile_coverage_rle =  IO::File->new($outfile_cr, "w");
my $outfile_flow_val_ins = IO::File->new($outfile_fvi, "w");
my $outfile_cf = IO::File->new($outfile_c, "w");

#only requires a rename.
if(! -e $output_base2qual)
{
	print "Trying cp $base2qualfile $output_base2qual\n";
	system("cp $base2qualfile $output_base2qual");
}

if(! -e $output_base2qual)
{
	warn  "It does not look like base2qual copied successfully\n";
}

foreach my $ref (@subDirToProcess)
{
	processAll($ref->[0], $ref->[1],  \%refbase2rle, $outfile_sub_in_complete, $outfile_flowvalues_to_len, 
		$outfile_read_start_pos, $outfile_coverage_rle, $outfile_flow_val_ins, $outfile_cf);
}

print "Done";

exit;

sub processAll($$$$$$$$$)
{
	#okay what outputs do you need...

	my ($dirPath, $subset, $refbase2rle, $outfile_sub_in_complete, $outfile_flowvalues_to_len,
                $outfile_read_start_pos, $outfile_coverage_rle, $outfile_flow_val_ins, $outfile_cf ) = @_;

	my %corrBases;
	my %substitutions;
	my %insertions;
	my %deletions;
	my %flow2base;
	my %cfie_pos;
        my %bounds_corr_bases;

        $dirPath  = $dirPath."/" if ($dirPath !~ /\/$/gi);

	opendir DIR, $dirPath or die "Could not open the subdirectory $dirPath\n";

	my ($corrRLE, $insBases, $subBases, $delBases);

	foreach my $element (readdir DIR)
	{
		if($element =~ /_corrRLEBases\.sql$/gi)
		{
			$corrRLE = $dirPath.$element;
		}
		elsif ($element =~ /_insBases\.sql$/gi)
		{
			$insBases = $dirPath.$element;
		}
		elsif ($element =~ /_subBases\.sql$/gi)
		{
			$subBases = $dirPath.$element;
		}
		elsif ($element =~ /_delBases\.sql$/gi)
		{
			$delBases = $dirPath.$element;
		}
	}

	my $flow2basefile = $dirPath."reads.flow";
	print "Flow to base: $flow2basefile\n";

	#The RLE bases have already been allocated a flow. Only flows which yield a base are recorded.
	
        loadCorrRLE($corrRLE, \%corrBases, \%bounds_corr_bases);

	#load only stuff relevant the current directory to process.
	#now, load flows, only based on the sequence IDs contained in corr bases.
	loadFlows($flow2basefile, \%flow2base, \%corrBases, \%bounds_corr_bases); #assume that every read has at least a correct base, duh!
        loadSubBases($subBases, \%substitutions);
	foreach my $s_read_id (keys %substitutions)
	{
		foreach my $s_rle_pos (keys %{$substitutions{$s_read_id}})
		{
			my ($s_strand, $s_ref_rle_pos, $s_nucleotide_ref, $s_nucleotide_read) = @{$substitutions{$s_read_id}->{$s_rle_pos}};
			my ($r_real_pos, $r_rle_len, $r_nucleotide) = @{$refbase2rle->{$s_ref_rle_pos}};
			my ($f_flow_pos, $f_read_pos, $f_nucleotide, $f_flow_val, $f_usable, $f_oop, $f_called_len);

			if(exists $flow2base{$s_read_id} and exists $flow2base{$s_read_id}->{$s_rle_pos})
			{
				($f_flow_pos, $f_read_pos, $f_nucleotide, $f_flow_val, $f_usable, $f_oop, $f_called_len) = @{$flow2base{$s_read_id}->{$s_rle_pos}};
			}
			else
			{
				($f_flow_pos, $f_read_pos, $f_nucleotide, $f_flow_val, $f_usable, $f_oop, $f_called_len) = ("NA")x 7;
			}
			my @res = ($s_read_id, $s_rle_pos, $f_read_pos, $f_flow_pos, $s_nucleotide_ref, $s_nucleotide_read, $f_flow_val, $r_real_pos, $s_strand, $r_rle_len, $f_usable, $f_oop, $f_called_len);
		
			$outfile_sub_in_complete->print(join("\t", @res)."\n");
		}
	}
	
	undef %substitutions;


	my %alignmentStartPos;
	foreach my $c_read_id (keys %corrBases)
	{
		foreach my $c_rle_pos (keys %{$corrBases{$c_read_id}})
		{
			my ($c_strand, $c_ref_rle_pos, $c_nucleotide) = @{$corrBases{$c_read_id}->{$c_rle_pos}};
			my ($r_real_pos, $r_rle_len, $r_nucleotide) = @{$refbase2rle->{$c_ref_rle_pos}};  #this did not stuff up
			my ($f_flow_pos, $f_read_pos, $f_nucleotide, $f_flow_val, $f_usable, $f_oop, $f_called_len);

                       #Have tried to correct for this by filtering out reads that have duplicated flow cycle bs.
                       if(exists $flow2base{$c_read_id} and exists $flow2base{$c_read_id}->{$c_rle_pos})
			{
				($f_flow_pos, $f_read_pos, $f_nucleotide, $f_flow_val, $f_usable, $f_oop, $f_called_len) =  @{$flow2base{$c_read_id}->{$c_rle_pos}};
			}		
			else
			{
				 ($f_flow_pos, $f_read_pos, $f_nucleotide, $f_flow_val, $f_usable, $f_oop, $f_called_len) =  ("NA") x 7; #ERROR: can't find read id or rle pos
			}

			my @res = ($c_read_id, $f_flow_val, $f_flow_pos, $c_rle_pos, $f_read_pos, $c_ref_rle_pos, $r_real_pos, $r_rle_len, $r_nucleotide, $c_strand, $f_usable, $f_oop, $f_called_len);

			$outfile_flowvalues_to_len->print(join("\t", @res)."\n");		
		
			if($c_rle_pos == 4)
			{
				$alignmentStartPos{$r_real_pos} = (exists $alignmentStartPos{$r_real_pos})?  $alignmentStartPos{$r_real_pos} + 1 : 1;
			}
		}
	}

	foreach my $r_real_pos (keys %alignmentStartPos)
	{
		my $count = $alignmentStartPos{$r_real_pos};
		$outfile_read_start_pos->print($count, "\t", $r_real_pos,"\n");
	}

	undef %alignmentStartPos;

        loadInsBases($insBases, \%insertions);
	foreach my $i_read_id (keys %insertions)
	{
		foreach my $i_rle_pos (keys %{$insertions{$i_read_id}})
		{
			my ($i_strand, $i_ref_pos_before, $i_ref_pos_after, $i_nucleotide_inserted, $i_nuc_before, $i_nuc_after, $i_ins_number) = @{$insertions{$i_read_id}->{$i_rle_pos}};
			my ($r1_real_pos, $r1_rle_len, $r1_nucleotide);

			if(! exists $refbase2rle->{$i_ref_pos_before})
			{
				($r1_real_pos, $r1_rle_len, $r1_nucleotide) = (-1,-1, "X");
			}
			else
			{
				($r1_real_pos, $r1_rle_len, $r1_nucleotide) = @{$refbase2rle->{$i_ref_pos_before}};
			}

			my ($r2_real_pos, $r2_rle_len, $r2_nucleotide);
			if(! exists $refbase2rle->{$i_ref_pos_after})
			{
				($r2_real_pos, $r2_rle_len, $r2_nucleotide) = (-1,-1,"X");
			}
			else
			{
				($r2_real_pos, $r2_rle_len, $r2_nucleotide) = @{$refbase2rle->{$i_ref_pos_after}};
			}

   		 	my ($f_flow_pos, $f_read_pos, $f_nucleotide, $f_flow_val, $f_usable, $f_oop, $f_called_len);

			if(exists $flow2base{$i_read_id} and exists $flow2base{$i_read_id}->{$i_rle_pos})
			{
				($f_flow_pos, $f_read_pos, $f_nucleotide, $f_flow_val, $f_usable, $f_oop, $f_called_len) =  @{$flow2base{$i_read_id}->{$i_rle_pos}};
			}
			else
			{
				($f_flow_pos, $f_read_pos, $f_nucleotide, $f_flow_val, $f_usable, $f_oop, $f_called_len) = ("NA") x 7; 
			}

			my @res = ($i_read_id, $f_flow_val, $f_flow_pos, $i_rle_pos, $f_read_pos, $r1_real_pos, $r2_real_pos, $i_ins_number, $i_nucleotide_inserted, $i_strand, $f_usable, $f_oop, $f_called_len);	

			$outfile_flow_val_ins->print(join("\t", @res), "\n");
		}
	}

	undef %insertions;
	undef %flow2base;

       my %temp;
       foreach my $c_read_id (keys %corrBases)
       {	
                foreach my $c_rle_pos (keys %{$corrBases{$c_read_id}})
                {
                        if($c_rle_pos < 4)
                        {
                                next;
                        }

                        my ($c_strand, $c_ref_rle_pos, $c_nucleotide) = @{$corrBases{$c_read_id}->{$c_rle_pos}};
                        my ($r_real_pos, $r_rle_len, $r_nucleotide) = @{$refbase2rle->{$c_ref_rle_pos}};

                        $temp{$r_real_pos}->{$r_nucleotide} = 0 if (! exists $temp{$r_real_pos}->{$r_nucleotide});
                        $temp{$r_real_pos}->{$r_nucleotide} = $temp{$r_real_pos}->{$r_nucleotide}  + 1;
                }
        }

        foreach my $ref_real_pos (sort {$a <=> $b} keys %temp)
        {
                foreach my $nucleotide (keys %{$temp{$ref_real_pos}})
                {
                        my @res = ($ref_real_pos, $nucleotide, $temp{$ref_real_pos}->{$nucleotide});
                        $outfile_coverage_rle->print(join("\t", @res)."\n");
                }
        }
        undef %temp;
}

sub loadFlows($$$$)
{
	my ($fullDirPath, $flow2base, $corrBases, $corrBasesBounds) = @_;

	my $fh = IO::File->new($fullDirPath, "r");
	my %earliest_corr_base;

	while(my $line =<$fh>)
	{
		chomp($line);

		next if($line =~ /FLOW_VAL/gi);
		
		my ($read_id, $flow_pos, $rle_pos, $read_pos, $nucleotide, $flow_val, $usable, $oop, $called_len) = split(/\t/, $line);

		my $start_for_align = $corrBasesBounds->{$read_id}->{'start'};
		my $end_for_align = $corrBasesBounds->{$read_id}->{'end'};

		if($flow_val >= 0.5 and exists $corrBases->{$read_id} and $rle_pos >= $start_for_align and $rle_pos <= $end_for_align) #store, its unique, otherwise requires another layer in the hash
		{
			$flow2base->{$read_id}->{$rle_pos} = [$flow_pos, $read_pos, $nucleotide, $flow_val, $usable, $oop, $called_len];
		}
	}
}

sub processReference($$)
{
	my ($dirPath, $refbase2rle) = @_;


	$dirPath = $dirPath."/" if ($dirPath !~ /\/$/gi);
	$dirPath = $dirPath."reference_rle.sql";

	my $fh = IO::File->new($dirPath, "r");

	while(my $line = <$fh>)
	{
		chomp($line);

		my ($reference_id, $rle_pos, $real_pos, $rle_len, $nucleotide) = split(/\t/, $line);
		$refbase2rle->{$rle_pos} = [$real_pos, $rle_len, $nucleotide];
	}
}


sub loadSubBases($$)
{
	my ($file, $substitutions) = @_;

	my $fh = IO::File->new($file, "r") or die "Could not open file $file\n";

	while(my $line = <$fh>)
	{
		chomp($line);
		my ($read_id, $reference, $alignment, $strand, $rle_pos, $ref_rle_pos, $nucleotide_ref, $nucleotide_read) = split(/\t/, $line);

		if($alignment == 0 and $rle_pos >= 4)
		{
			$substitutions->{$read_id}->{$rle_pos} = [$strand, $ref_rle_pos, $nucleotide_ref, $nucleotide_read];
		}

	}
}


sub loadInsBases($$)
{
	my ($file, $insertions) =  @_;
	my $fh = IO::File->new($file, "r") or die "Could not open file $file\n";

	while(my $line = <$fh>)
	{
		chomp($line);
		my ($read_id, $reference, $alignment, $strand, $rle_pos, $ref_pos_before, $ref_pos_after, $nucleotide_inserted, $nuc_before, $nuc_after, $ins_number) = split(/\t/, $line);
		
		if($alignment == 0 and $rle_pos >= 4)
		{
			$insertions->{$read_id}->{$rle_pos} = [$strand, $ref_pos_before,$ref_pos_after, $nucleotide_inserted, $nuc_before, $nuc_after, $ins_number];
		}
	}
}



sub loadCorrRLE($$$)
{
	my ($file, $corrBases, $corrBasesBounds) = @_;
	my $fh = IO::File->new($file, "r") or die "Could not open file $file \n";

	while(my $line = <$fh>)
	{
		chomp($line);
		my ($readID, $reference, $align_no, $strand, $rle_pos, $ref_rle_pos, $nucleotide) = split ("\t", $line);

		if($align_no == 0 and $rle_pos >= 4)
		{
			$corrBases->{$readID}->{$rle_pos} = [$strand, $ref_rle_pos, $nucleotide];
			if(! exists $corrBasesBounds->{$readID}->{'start'} or $corrBasesBounds->{$readID}->{'start'} > $rle_pos)
			{
				$corrBasesBounds->{$readID}->{'start'} = $rle_pos;
			}
		
			if(! exists $corrBasesBounds->{$readID}->{'end'} or $corrBasesBounds->{$readID}->{'end'} < $rle_pos)
			{
				$corrBasesBounds->{$readID}->{'end'} = $rle_pos;
			}
		}
	}
}	
