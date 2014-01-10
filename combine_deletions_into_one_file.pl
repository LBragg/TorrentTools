#! /usr/bin/perl -w

use strict;
use IO::File;
sub loadSubDirDeletions($$);
sub processDeletions($$$);

# This script is require dfor combining the files together. Actually think alot of the code is no longer necessary...
die "Require [top directory] [outfile] to start" if (scalar(@ARGV) != 2);

my $dirname = $ARGV[0];
my $outfile = $ARGV[1];

$dirname = $dirname."/" if ($dirname !~ /\/$/gi);

opendir (my $dir, $dirname) or die "Could not open directory\n";
my $out_fh = IO::File->new($outfile, "w");


my $refFH;
my @toProcess;

while(my $entry = readdir $dir)
{
	if(-d $entry)
	{
		$entry = $entry."/";
	}

	if($entry =~ /read\_dir\_([0-9]+)/gi)
	{
		my $captured = $1;
		my $fullDirPath = $dirname."/".$entry;

		opendir (my $sub_dir, $fullDirPath) or die "Could not open subdirectory $fullDirPath";

		my $deletionsFile = undef;
		my $flowFile = undef;

		while(my $inner = readdir $sub_dir)
		{
			if ($inner =~ /delBases\.sql$/gi)
			{
				$deletionsFile = $inner;
			}
			elsif($inner =~ /\.flow$/gi)
			{
				$flowFile = $inner;
			}
		}
		my $fullPathDel = $dirname.$entry."/".$deletionsFile;
		my $fullPathFlow = $dirname.$entry."/".$flowFile;
		push @toProcess, [$fullPathDel, $fullPathFlow];
		
	}

	if($entry =~ /reference\_rle\_for/gi)
	{
		$refFH = $dirname.$entry."/"."reference_rle.sql";
	}
}

my @refToBase;
my $rFH = IO::File->new($refFH, "r");

#store every base in the reference.

while(my $line = <$rFH>)
{
	chomp($line);
	my ($id, $rlePos, $basePos, $rleLength, $base) = (split(/\t/, $line));
	push @refToBase, [$basePos, $rleLength, $base]; #they start from zero also
}

#for each deletions file to process...

foreach my $ref (@toProcess)
{
	my ($deletionsFile, $flowFile) = @{$ref};

	my %deletions; #seem to process each sub directory separately.

	#want to load the deletions file into the deletions hash???
	loadSubDirDeletions($deletionsFile, \%deletions);
	processDeletions($flowFile, \%deletions, $out_fh);
}

sub processDeletions($$$)
{
	my ($flowFilePath, $deletionsHash, $out_fh) = @_;
	my %deletions = %{$deletionsHash};
	my $flowFile = IO::File->new($flowFilePath, "r");
	my $currRead = undef;
	my %startPosToFlow;
	my %endPosToFlow;

	#first pass.
	my %rlePosToBasePos;
	my %readToOOP;
	
while(my $line = <$flowFile>)
{	
	chomp($line);
        my ($read_id, $flow_position, $rle_position, $read_position, $nucleotide, $flow_value, $usable, $oop) = split(/\t/, $line);
	$readToOOP{$read_id} = $oop;

	if($read_id eq "READ_ID")
	{
		next;
	}
	if(exists $deletions{$read_id})
	{
		$rlePosToBasePos{$read_id}->{$rle_position} = $read_position;

		foreach my $startPos (keys %{$deletions{$read_id}})
		{
			if($rle_position == $startPos and $flow_value > 0.5)
			{
				$startPosToFlow{$read_id}->{$rle_position} = $flow_position;
			}
		
                        foreach my $endPos(keys %{$deletions{$read_id}->{$startPos}})
                        {
  				if($rle_position == $endPos and $flow_value >= 0.5)
				{
					$endPosToFlow{$read_id}->{$rle_position} = $flow_position;
				}
			}
		}
	}
}

#second pass.
$flowFile = IO::File->new($flowFilePath, "r");
my %processed;
my %assigned;

#output, and check whether things were not processed.
#READ_ID FLOW_VALUE      FLOW_POSITION   RLE_POSITION    READ_POSITION   REF_RLE_POSITION        REF_REAL_POSITION       RLE_LENGTH      NUCLEOTIDE_DELETED
#printing to standard out
foreach my $read_id (keys %deletions)
{
	foreach my $start (keys %{$deletions{$read_id}})
	{
		foreach my $end (keys %{$deletions{$read_id}->{$start}})
		{
			foreach my $delNumber (keys %{$deletions{$read_id}->{$start}->{$end}})
			{
					my ($baseDeleted, $refpos, $strand) = @{$deletions{$read_id}->{$start}->{$end}->{$delNumber}};
					my ($refRealPosition, $rleLength, $base) = @{$refToBase[$refpos]};
					my $readPos = $rlePosToBasePos{$read_id}->{$start} || "N/A";					
					my $oop = $readToOOP{$read_id};

					$out_fh->print(join("\t", $read_id, "N/A", "N/A", $start, $readPos, $refpos, $refRealPosition, $rleLength, $base, $strand, "False", $oop)."\n");
			}
		}
	}
}
}


#does not look like this code actually does anything.
sub loadSubDirDeletions($$)
{
	my ($file, $delhash) = @_;
	my $fh = IO::File->new($file, "r");

	while(my $line = <$fh>)
	{
		chomp($line);

		my ($read, $ref, $alignmentNumber, $strand, $refpos, $basePosBefore, $basePosAfter, $baseDeleted, $baseBefore, $baseAfter, $delNumber) = split(/\t/, $line);

		if($alignmentNumber != 0)
		{
			next;
		}
	

		#lets perform our corrections
		if($strand == "-1")	
		{
			$basePosBefore = $basePosBefore - 2; #changed to 3... somereason its quite off
			$basePosAfter = $basePosAfter - 2;

			my $tmp = $basePosBefore;
			$basePosBefore = $basePosAfter;
			$basePosAfter = $tmp;	

		
			#they've been swapped.
		}
		else
		{
			$basePosBefore = $basePosBefore;
			$basePosAfter = $basePosAfter;
		}

		$delhash->{$read}->{$basePosBefore}->{$basePosAfter}->{$delNumber} = [$baseDeleted, $refpos, $strand];
	}
}

