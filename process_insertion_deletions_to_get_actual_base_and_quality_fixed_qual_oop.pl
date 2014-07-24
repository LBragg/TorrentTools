#! /usr/bin/perl -w

use strict;
use IO::File;
use Bio::SeqIO;

die "Require [base2qual.txt] [flowvals_insertions] [deletions] [corrbases] [substitution information] [output path]" if(scalar(@ARGV) != 6);

my $base2qualfh = IO::File->new($ARGV[0], "r");
my $ins = IO::File->new($ARGV[1], "r");
my $del = IO::File->new($ARGV[2], "r");
my $corrBases = IO::File->new($ARGV[3], "r");
my $subs = IO::File->new($ARGV[4], "r"); #fix this later

my $output_path = $ARGV[5];
$output_path = $output_path."/" if($output_path !~ /\/$/gi);

my %read2insertions;
my %read2del;
my %corrBasesHash;

my %base2qual; #for other analyses
my $out_error2base = IO::File->new($output_path."errors2basepos.txt", "w"); # we omit the fixed qual, because it is always fixed now.

while(my $line = <$base2qualfh>)
{
	chomp($line);
	if($line =~ /READ/gi)
	{
		next;
	}	

	my ($read, $base, $qual) = split (/\t/, $line);
	$base2qual{$read}->{$base} = $qual; #load base 2 qual.
}

## Add OOP, valid FV.

$out_error2base->print(join("\t", "Read", "BasePosition", "FlowPosition", "RLEPosition", "RefPosition", "Quality", "Strand", "Type", "ValidFV", "OOP"), "\n");


{
while(my $line = <$ins>)
{
	chomp($line);

	if($line =~ /^READ_ID/)
	{
		next;
	}
	else
	{
		# Add OOP, valid FV, calledLen.
		my ($read_id, $flow_value, $flow_pos, $rle_pos, $read_pos, $ref_pos_before, $ref_pos_after, $ins_number, $nucleotide_inserted, $strand, $usable, $oop, $callLen) = split(/\t/, $line);
			
		#need to know the start and end position of the flow stuff.
		my $numBasesInserted = $callLen;
		my $startRead = $read_pos;

		if($numBasesInserted eq "NA")
		{
			next;
		}
	
		for(my $i = 0; $i < $numBasesInserted; $i++)
		{
			my $pos = $startRead + $i;
			my $qual = $base2qual{$read_id}->{$pos};
			if(! defined $qual)
			{
				warn "$read_id $pos undefined for INSERTION \n";
		
				$qual = "N/A";
			}
			$out_error2base->print(join("\t", $read_id, $pos, $flow_pos, $rle_pos, $ref_pos_before, $qual, $strand, $usable, $oop, "INSERTION"), "\n");
		}
	}
}


}

#load process deletions pos

{
while(my $line = <$del>)
{
	chomp($line);
	
	if($line =~ /^READ_ID/gi)
	{
		next;
	}	
	else
	{
		my ($read_id, $flow_val, $flow_pos, $rle_pos, $read_pos, $ref_rle_pos, $ref_pos, $rle_length, $nucleotide_deleted, $strand, $usable, $oop) = split(/\t/, $line);
		#no call len, but it is zero.
	
		#nucleotide deleted is not very useful. But I can take the read position thing...
	#	$read2delPos{$read_id}->{$read_pos} = [$flow_pos, $rle_pos, $ref_pos]; #how will I usee this information.

		#why would the quality for a base be undefined??
		my $qual = $base2qual{$read_id}->{$read_pos};

		if(! defined $qual)
                {
                        warn "$read_id $read_pos undefined for DELETION\n";

                	$qual = "N/A";
                }

		#put 1 in for positive strand
		$out_error2base->print(join("\t", $read_id, $read_pos, $flow_pos, $rle_pos, $ref_pos, $qual, "1", $usable, $oop, "DELETION"), "\n");
	}
}

}




#load process deletions neg 

#I want an output file that says an error occurred at this base position and this flow
{
while(my $line = <$corrBases>)
{
	chomp($line);
	
	my ($read_id, $flow_value, $flow_position, $read_rle_position, $read_position, $ref_rle_position, $ref_real_position,  $rle_length, $nucleotide,$strand, $usable, $oop, $callLen)  = split(/\t/, $line);

	if($line =~ /READ_ID/gi)
	{
		next;
	}

	if($flow_value eq "NA") #this is now occurring as a consequence of the stupid changes in how flow-values are reported by Ion Torrent
	{
		next; 
	}

	
	my $start = $read_position;
	my $end = $callLen + $read_position - 1;
	
	#problem can consider errors at flow or base level
	#one error at flow level can appear as multiple errors at base level
	my $numErrors = ($rle_length - $callLen); #deletions will only be recorded once.
	
	#what do you want to use this for?
	if($numErrors < 0)
	{
		$numErrors = abs($numErrors); #these are all insertion errors

		for(my $pos = $end; $pos >= $start ; $pos--)
		{
			if($numErrors > 0)
			{
			       	my $qual = $base2qual{$read_id}->{$pos};
		               	if(! defined $qual)
        			{
 	                		warn "$read_id $pos was undefined for INSERTION 2\n";
                       			$qual = "N/A";
               			}
 		               $out_error2base->print(join("\t", $read_id, $pos, $flow_position, $read_rle_position, $ref_real_position, $qual, $strand, $usable, $oop, "INSERTION"), "\n");
			}
			else
			{	

			       my $qual = $base2qual{$read_id}->{$pos};

			       if(! defined $qual)
                               {
                                       warn "$read_id $pos was undefined for CORRECT 1\n";
                                       $qual = "N/A";
                               }


                               $out_error2base->print(join("\t", $read_id, $pos, $flow_position, $read_rle_position, $ref_real_position, $qual, $strand, $usable, $oop, "CORRECT"), "\n");

			}
			$numErrors--;
		}
	}
	elsif($numErrors > 0)
	{
		#deletion occurs at position
		my $pos = $start;
	
		for($pos = $start; $pos <= $end; $pos++)
		{
			  my $qual = $base2qual{$read_id}->{$pos};
			  if(! defined $qual)
                          {
                          	warn "$read_id $pos was undefined for CORRECT 2\n";
                                $qual = "N/A";
                          }
		
                         $out_error2base->print(join("\t", $read_id, $pos, $flow_position, $read_rle_position, $ref_real_position, $qual, $strand, $usable, $oop, "CORRECT"), "\n");
		}

		for(my $i = 0; $i < $numErrors; $i++)
		{
                	$out_error2base->print(join("\t", $read_id, $pos, $flow_position, $read_rle_position, $ref_real_position + $i + 1, "N/A", $strand, $usable, $oop, "DELETION"), "\n");
		}
	}	
	else
	{	
		for(my $pos = $start; $pos <= $end; $pos++)
		{
			 my $qual = $base2qual{$read_id}->{$pos};

			  if(! defined $qual)
                          {
                                warn "$read_id $pos was undefined for CORRECT 3\n";
                                $qual = "N/A";
                          }


                         $out_error2base->print(join("\t", $read_id, $pos, $flow_position, $read_rle_position, $ref_real_position, $qual, $strand, $usable, $oop, "CORRECT"), "\n");
		}
	}
}
}




my %subs;
while(my $line = <$subs>)
{
	chomp($line);
	if($line =~ /^READ_ID/gi)
	{
		next;
	}
	else
	{
		#READ_ID RLE_POSITION    NUCLEOTIDE_REF  NUCLEOTIDE_READ FLOW_VALUE      FLOW_POSITION
		my ($read_id, $rle_position, $base_position, $flow_position, $nucleotide_ref, $nucleotide_read, $flow_value, $ref_base_position, $strand, $ref_rle_len, $usable, $oop, $callLen) = split(/\t/, $line);
		my $qual = $base2qual{$read_id}->{$base_position};
	
		if(! defined $qual)
                {
                	warn "$read_id $base_position was undefined for substitution\n";
                        $qual = "N/A";
                }

		$out_error2base->print(join("\t", $read_id, $base_position, $flow_position, $rle_position, $ref_base_position, $qual, $strand, $usable, $oop, "SUSTITUTION"), "\n");		
	}

}

