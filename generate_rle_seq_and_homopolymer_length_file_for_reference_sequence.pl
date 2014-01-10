#! /usr/bin/perl -w

use strict;
use IO::File;
use Bio::SeqIO;
use Bio::Seq;

sub convertToRLE($);
sub writeHomopolymers($$$$);

## Given a FASTA file generates the RLE sequence, and corresponding HP length file.

die "Require [input_fasta] [output prefix] as arguments \n" if (scalar(@ARGV) != 2);

my $input_fasta = Bio::SeqIO->new(-file => $ARGV[0], -format => "FASTA");
my $output_fasta_rle = Bio::SeqIO->new(-file => ">".$ARGV[1]."_rle.fasta", -format => "FASTA");
my $output_hp_sql = IO::File->new($ARGV[1]."_rle.sql", "w");

while(my $seq = $input_fasta->next_seq())
{
	my $seq_id = $seq->id();
	if($seq_id =~ /(NC_[0-9]+(?:\.[0-9]+))/gi)
	{
		$seq_id = $1;
	}
	my ($rle_seq, $hp_counts)  = @{convertToRLE($seq)};
	$output_fasta_rle->write_seq($rle_seq);

	writeHomopolymers($output_hp_sql, $seq_id, $hp_counts, $rle_seq);
}

sub convertToRLE($)
{
	my $seq = shift;
	my $sequence = $seq->seq();
	
	my @RLESeq;
	my @RLEhp;

	while($sequence =~ /([ATGNC])(\1*)/gi)
	{
		my $captured = $1.$2;
		my $length_captured = length ($captured);
		push @RLESeq, (split(//, $captured))[0];
		push @RLEhp, $length_captured;
	}
			
	my $newSeq = Bio::Seq->new(-id => $seq->id(), -seq => join("", @RLESeq));

	return [$newSeq, \@RLEhp];
} 

#this is for the reference sequence.
sub writeHomopolymers($$$$)
{
	my ($fh, $seq_id, $counts, $rle_seq) = @_;
	my $readPosition = 0;

	my @bases = split(//, $rle_seq->seq());

	for(my $index = 0; $index < scalar(@{$counts}); $index++)
	{
		my $length = $counts->[$index];
		my $base = $bases[$index];
		
		$fh->print( join("\t", $seq_id, $index, $readPosition, $length, $base), "\n");
		$readPosition += $length;
	}
}
