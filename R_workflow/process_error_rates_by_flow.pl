#! /usr/bin/perl -w

use strict;
use IO::File;


die "Require [error rate] [error rate for zero] [output dir] [output_prefix] as input " if (scalar(@ARGV) != 4);

#This appears to be aggregated for Acacia models, but stuffed if I know.

my $inputGeneral = IO::File->new($ARGV[0], "r");
my $inputZero = IO::File->new($ARGV[1], "r");

my %ref2flow;
my $maxFlowPos = 0;

my $output_dir = $ARGV[2];

$output_dir = $output_dir."/" if ($output_dir !~ /\/$/gi);

my $output_prefix = $ARGV[3];

while(my $line = <$inputGeneral>)
{
	chomp($line);
	
	if($line =~ /FlowPos/gi)
	{
		next;
	}

	my ($undef, $ref_len, $flow_pos, $count_below, $count_at, $count_above) = split(/\s+/, $line);

	if($flow_pos !~ /[0-9]+/gi)
	{
		print "Flow Pos $flow_pos is not numeric\n";
		next;
	}	

	$ref2flow{$ref_len}->{$flow_pos} = [$count_below, $count_at, $count_above];

	if($flow_pos > $maxFlowPos)
	{
		$maxFlowPos = $flow_pos;
	}
}

while(my $line = <$inputZero>)
{
	if($line =~ /RefLen/)
	{
		next;
	}
	
	my ($undef, $ref_len, $flow_pos,$count_above, $count_at) = split(/\s+/, $line);

	if($flow_pos !~ /[0-9]+/gi)
	{
		print "Flow pos $flow_pos is NA\n";
		next;
	}

	
	$ref2flow{$ref_len}->{$flow_pos} = [0, $count_at, $count_above];
}



my $flowCycle = "TACGTACGTCTGAGCATCGATCGATGTACAGC";

#with flow cycle, models that take into account

# REF_LEN, FLOW SEGMENT, NUCLEOTIDE
# REF_LEN, FLOW_SEGMENT, FLOW_POSITION

my %amalgamatedByNuc;
my %amalgamatedByFlowPos;

#do this tomorrow

if(!exists $ref2flow{0})
{
	$ref2flow{0} = undef; #exists but stores nothing.
}


foreach my $reflen (keys %ref2flow)
{
	for(my $flow_pos = 0; $flow_pos < $maxFlowPos; $flow_pos++)
	{
		my $flowSegment = int($flow_pos / length($flowCycle));
		my $posInFlow = $flow_pos % length($flowCycle);
		my $nucleotide = (split(//, $flowCycle))[$posInFlow];

		my $base = "AT";
		
		if($nucleotide =~ /[GCgc]/)
		{
			$base = "GC";
		}

		if(! exists $ref2flow{$reflen}->{$flow_pos} and $reflen == 0)
		{
			$amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'below'}=  0;
			$amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'at'}=  1;
			$amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'above'}=  0;

			$amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'below'} = 0;
			$amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'at'} = 1;
			$amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'above'} = 0;
		}
		elsif(exists $ref2flow{$reflen}->{$flow_pos})
		{
			if(exists $amalgamatedByNuc{$reflen}->{$flowSegment}->{$base})
			{
				 $amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'below'}=  $amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'below'} + $ref2flow{$reflen}->{$flow_pos}[0];
				 $amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'at'}= $amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'at'} + $ref2flow{$reflen}->{$flow_pos}[1];
				 $amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'above'}= $amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'above'} + $ref2flow{$reflen}->{$flow_pos}[2];
			}
			else
			{	
				$amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'below'} =  $ref2flow{$reflen}->{$flow_pos}[0];
				$amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'at'} = $ref2flow{$reflen}->{$flow_pos}[1];
				$amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'above'} =   $ref2flow{$reflen}->{$flow_pos}[2];
			}

			if(exists $amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow})
			{
				$amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'below'} =  $amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'below'} +  $ref2flow{$reflen}->{$flow_pos}[0];
				$amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'at'} =  $amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'at'} +  $ref2flow{$reflen}->{$flow_pos}[1];
				$amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'above'} =  $amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'above'} +  $ref2flow{$reflen}->{$flow_pos}[2];
			}
			else
			{
				$amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'below'} =  $ref2flow{$reflen}->{$flow_pos}[0];
				$amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'at'} =  $ref2flow{$reflen}->{$flow_pos}[1];
				$amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'above'} = $ref2flow{$reflen}->{$flow_pos}[2];	
			}
		}
	}
}


#firstly, prepare the output file

my $output_by_nuc = IO::File->new($output_dir.$output_prefix."_error_rate_by_flow_segment_and_nuc.csv", "w");
my $output_by_flow_segment_and_flow_pos = IO::File->new($output_dir.$output_prefix."error_rate_by_flow_segment_and_pos_in_cycle.csv", "w");


foreach my $reflen (sort {$a <=> $b} keys %amalgamatedByFlowPos)
{	
	foreach my $flowSegment (sort {$a <=> $b} keys %{$amalgamatedByFlowPos{$reflen}})
	{
		foreach my $posInFlow (sort {$a <=> $b} keys %{$amalgamatedByFlowPos{$reflen}->{$flowSegment}})
		{
			my $below = $amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'below'};
			my $at =  $amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'at'};
			my $above = $amalgamatedByFlowPos{$reflen}->{$flowSegment}->{$posInFlow}->{'above'};
			my $total = $below + $at + $above;
		
			my $propBelow = $below / $total;
			my $propAt = $at / $total;
			my $propAbove = $above / $total;

			$output_by_flow_segment_and_flow_pos->print(join(",", ($reflen, $flowSegment, $posInFlow, $propBelow, $propAt, $propAbove, $total))."\n");
		}
	}
}

foreach my $reflen (sort {$a <=> $b} keys %amalgamatedByNuc)
{
        foreach my $flowSegment (sort {$a <=> $b} keys %{$amalgamatedByNuc{$reflen}})
        {
                foreach my $base (sort {$a cmp $b} keys %{$amalgamatedByNuc{$reflen}->{$flowSegment}})
                {
                        my $below = $amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'below'};
                        my $at =  $amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'at'};
                        my $above = $amalgamatedByNuc{$reflen}->{$flowSegment}->{$base}->{'above'};
                        my $total = $below + $at + $above;
                        my $propBelow = $below / $total;
                        my $propAt = $at / $total;
                        my $propAbove = $above / $total;
                        $output_by_nuc->print(join(",", ($reflen, $flowSegment, $base, $propBelow, $propAt, $propAbove, $total))."\n");
                }
        }
}








