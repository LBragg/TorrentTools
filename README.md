TorrentTools
============

Scripts for dealing with Ion Torrent BAM/SFF files.

For most readers, the only tools of interest are the BAM parser and associated libraries. This suite of tools is designed to 
extract PGM reads from a BAM file, align the reads in RLE format (using a 3rd party aligner) against a reference, and 
characterise the error occurrence with respect to factors such as base and flow position.

With recent changes in Ion Torrent PGM output, the BAM only supplies normalised but not phase-corrected flow-values. These 
do not correspond to the called sequence in the BAM file, thus making it difficult to assess error-rates at the base and flow
level.

This suite of tools attempts to identify as many valid flow-values as possible, while marking called bases with no corresponding
flow as an 'InvalidFlow'. In many cases, when there is no obvious reason for the modified base-call (due to phase-correction),
the flowgram is considered out of phase (OOP) with the bases as called by Ion Torrent phase-correction suite.

In this instance, no FlowValues can be provided for these bases. Once the flows are OOP with the read, the remaining read bases
are assigned to 'fake' flow positions with negative indexing.

The remainder of the workflow consists of R analysis scripts for characterising PGM data.

