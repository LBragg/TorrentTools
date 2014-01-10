#! /usr/bin/perl -w

use strict;
use IO::File;
use Bio::SeqIO;

## Parses segemehl alignments of run-length-encoded reads against a run-length-encoded reference

die "Require [segeheml output] [read sequences] [reference sequence] [output_prefix] [verbose] \n" if (scalar(@ARGV) != 5);

my $input = IO::File->new($ARGV[0], "r");
my $sizeGC = IO::File->new($ARGV[1], "r");

my $output_prefix = $ARGV[3];

#load these in.
my $reference_fh = Bio::SeqIO->new(-file => $ARGV[2], -format => "FASTA");
my $ref_seq = $reference_fh->next_seq();
my $old_id =  $ref_seq->id();

my $new_id;
if($old_id =~ /(NC_[0-9]+(?:\.[0-9]+))/gi)
{
	$new_id = $1;
}
else
{
	$new_id = $old_id;
	$new_id =~ s/\s+//gi;
}

my $verbose = $ARGV[4] =~ /(True|T)/gi;


my %reads;
my $read_fh = Bio::SeqIO->new(-file => $ARGV[1], -format => "FASTA");

while(my $seq = $read_fh->next_seq())
{
	$reads{$seq->id()} = $seq;
}

my $sql_corr = IO::File->new($output_prefix."_corrRLEBases.sql", "w");
my $sql_ins = IO::File->new($output_prefix."_insBases.sql", "w");
my $sql_subs = IO::File->new($output_prefix."_subBases.sql", "w");
my $sql_dels = IO::File->new($output_prefix."_delBases.sql", "w");


my %numAlignments;

LINE:while(my $line = <$input>)
{
	chomp($line);
	if($line =~ /^@/gi)
	{
		next LINE;
	}

	my @fields = split(/\t/, $line);
	my ($desc, $flag, $rname, $pos, $mappingQuality, $cigar, $rnext, $pnext, $tlen, $seq, $qual, $identity_field) =  split(/\t/, $line);
	my $editDist = (split(":", $identity_field))[2];
	my $alignmentNumber = 0;

	if(exists $numAlignments{$desc})
	{
		$alignmentNumber = $numAlignments{$desc} + 1;
	}

	$numAlignments{$desc} = $alignmentNumber;


	#pos is start of alignment in read.
	my @aligns;
	my $strand = 1;

	if($flag == 0x10)
	{
		$strand = -1;
	}

	#need separate because of I and D
	my $left_start_read = undef;
	my $right_end_read = undef;

	#need separate because of I and D
	my $left_start_ref = undef;
	my $right_end_ref = undef;

	my $total_ref = 0;
	my $total_read = 0;
	my $special_case = 0;


	my @operations;

	my @tmpOp;

	my $opSum = 0;

	while($cigar =~ /([0-9]+)([MIDNSHP=X])/gi)
	{
		push @tmpOp, [$1, $2];
		
		$opSum += $1;
	}

	


#	if($strand != 1)
#	{
#		@operations = reverse @tmpOp;
#	}
#	else
#	{
		@operations = @tmpOp;
#	}



	foreach my $operation (@operations)
	{
		my ($length ,$type) = @{$operation};

		#double check this is true. 
		if($type !~ /[HSD]/gi and ! defined $left_start_read)
		{
			$left_start_read = $total_read;
		}
		
		if($type =~ /[HS]/gi and defined $left_start_read)
		{
			$right_end_read = $total_read;
		}
		
		if($type !~ /[HSI]/gi and ! defined $left_start_ref)
		{
			$left_start_ref = $total_ref;
		}
	
		if($type =~ /[HS]/gi and defined $left_start_ref)
		{
			$right_end_ref = $total_ref;
		}
	
		if($verbose)
		{
			print "Out: $length $type\n";
		}

		push @aligns, [$length, $type] ;

		if($type ne "I")
		{
			$total_ref += $length;
		}
		
		if($type ne "D")
		{
			$total_read += $length;
		}
	}

	if(! defined $right_end_read)
	{
		$right_end_read = $total_read;
	}
	if(! defined $right_end_ref)
	{
		$right_end_ref = $total_ref;
	}

	#one of these probably needs + 1; #remember using 0 coordinate system with other stuff.

	my $corrReadStart = ($strand eq -1)? $total_read - $right_end_read + 1  : $left_start_read + 1;
	my $corrReadEnd = ($strand eq -1)? $total_read - $left_start_read : $right_end_read;

	if($verbose)
 	{
		print "Corr read start: $corrReadStart, Corr Read End: $corrReadEnd\n";
	}


	my $refStart = $pos + $left_start_ref;
	my $refEnd = $pos + $right_end_ref - 1;	# also maybe - 1. or + 1.

	if($verbose)
	{
		print "Ref start: $refStart, Ref end: $refEnd\n";
	}

	my $id = ($desc =~ /^([^\s+]+)/gi)[0];

	my $total_match = abs($corrReadEnd - $corrReadStart);

	if($verbose)
	{
		print "The total match for the read is $corrReadEnd - $corrReadStart\n";
	}


	if($alignmentNumber > 0) #ignore those with multiple matches
	{
		#we don't consider follow matches.
	}
	else #uniquely located
	{
		#lets take a look at the sequences

		my $read = $reads{$id};

		if(! defined $read)
		{
			warn "Read is not defined.\n"

		}
 
		my $subRead;
	
		$subRead  = $read->subseq($corrReadStart, $corrReadEnd);
	
		if($verbose)
		{
			print "Aligned subsequence: <$subRead>\n";
		}

		my $subRef = $ref_seq->subseq($refStart, $refEnd);

		#so we are doing the reverse compliment, remember that!!!	
		if($strand eq -1)
		{
			$subRead =~ tr/[ATGCatgc]/[TACGtacg]/;
			$subRead = reverse $subRead;
		}

		#handle stranded ness
		my $offset = ($strand == 1)? 1:  -1;
	
		#changed this from + to 1;

		my $currPosReadArray = 0;
		my $currPosRefArray = 0;			

		my @subReadBases = split(//, $subRead);
		my @subRefBases = split(//, $subRef);

		#finish off here.		
		
		my $mismatchCount = 0;

		foreach my $operation (@operations)
		{
			my ($length, $type) = @{$operation};
			
			if($type eq "M") #doesn't necessarily mean it matches, could be substitution
			{
				for(my $i =0; $i < $length; $i++)
				{

					my $actualPosRead = ($strand eq "-1")? $total_read - $currPosReadArray - 1 : $corrReadStart + $currPosReadArray - 1 ;

					if($verbose)
					{
						print "M: $actualPosRead ".$subReadBases[$currPosReadArray]."\n";
					}

					#this is weird, would the strand effect the reference position?	
					my $posInRef = 	$refStart + $currPosRefArray - 1;				

					#RLE coordinate system starts from zero in database, therefore this must start from zero also.
					#reporting base on reverse compliemented strand.

					if(uc($subReadBases[$currPosReadArray]) ne uc($subRefBases[$currPosRefArray]))
					{
						$sql_subs->print(join("\t", $id, $new_id, $alignmentNumber, $strand, $actualPosRead, $posInRef, $subRefBases[$currPosRefArray], $subReadBases[$currPosReadArray])."\n");

						$mismatchCount++;
					}
					else
					{
						#warn "Matches: $actualPosRead ".($subReadBases[$currPosReadArray]).", pos ref ".($subRefBases[$currPosRefArray])."\n";
						$sql_corr->print(join("\t", $id, $new_id, $alignmentNumber, $strand, $actualPosRead, $posInRef, $subRefBases[$currPosRefArray])."\n");

					}
					$currPosReadArray++;
					$currPosRefArray++;
				}
			}
			elsif ($type eq "I")
			{
				$mismatchCount += $length;
		
				for(my $i = 0; $i < $length; $i++)
				{
					my $insInRead = ($strand eq "-1")? $total_read - $currPosReadArray - 1 : $corrReadStart + $currPosReadArray - 1;
					my $base = $subReadBases[$currPosReadArray];
			
					my $posInRefBefore =  $currPosRefArray - 1; #base before
					my $posInRefAfter =  $currPosRefArray;
	
					my $actualPosRefBefore = $refStart + $posInRefBefore - 1;
					my $actualPosRefAfter = $refStart + $posInRefBefore;

					my $baseBefore = ($posInRefBefore >= 0)? $subRefBases[$posInRefBefore]: "X";
					my $baseAfter = ($posInRefAfter >= 0)? $subRefBases[$posInRefAfter] : "X";

					if($posInRefAfter >= scalar(@subRefBases))
					{
						$baseAfter = "X";
					}
		
					if($verbose)
					{
						print "I: $base $i $insInRead $baseBefore $baseAfter\n";
					}
					
	
					$sql_ins->print(join("\t", $id, $new_id,$alignmentNumber, $strand, $insInRead, $actualPosRefBefore, $actualPosRefAfter, $base, $baseBefore, $baseAfter, $i)."\n");
			
					$currPosReadArray++;
				}
			}
			elsif($type eq "D")
			{
				$mismatchCount += $length;
				for(my $i = 0; $i < $length; $i++)
				{
					my $delInRef = $refStart + $currPosRefArray - 1;

					my $base = $subRefBases[$currPosRefArray];
					
					my $posInReadBefore = $currPosReadArray - 1;
					my $posInReadAfter = $currPosReadArray;

					
					my $actualPosInReadBefore = ($strand eq "-1")? $total_read - $posInReadBefore + 1 : $corrReadStart + $posInReadBefore - 1;
					my $actualPosInReadAfter = ($strand eq "-1") ? $total_read - $posInReadAfter + 1 : $corrReadStart + $posInReadAfter - 1; #changed the second bit to make sense too??

					my $baseBefore = ($posInReadBefore >= 0)? $subReadBases[$posInReadBefore] : "X";
					my $baseAfter = ($posInReadAfter >= 0) ? $subReadBases[$posInReadAfter] : "X";
		
			
					if($verbose)
					{
						print "D: $base $actualPosInReadBefore $actualPosInReadAfter\n"; 
					}

					$sql_dels->print(join("\t",$id, $new_id, $alignmentNumber, $strand, $delInRef, $actualPosInReadBefore, $actualPosInReadAfter,$base, $baseBefore, $baseAfter, $i )."\n");
					$currPosRefArray++;			

				}
				
			}
			elsif($type eq "H" or $type eq "S")
			{
				#just ignore.
			}
		}

		if($mismatchCount != $editDist)
		{
		}
	}
}
