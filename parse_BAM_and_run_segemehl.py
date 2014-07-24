'''
Created on 10/12/2013

@author: bra427
'''

#! /usr/bin/python
import argparse
import os, sys
from lib.ITBAMIterator import *
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from subprocess import call
from subprocess import check_output
from subprocess import list2cmdline
from subprocess import Popen
from subprocess import Popen, PIPE
import shlex

#CONSTANTS
segemehl_exe = "/opt/segemehl_0_1_6/segemehl_0_1_7/segemehl/segemehl.x" ## segemehl should be in the path too. #TODO remove hard-coded location.
path2scripts = os.path.split(os.path.abspath(sys.argv[0]))[0] + "/" #was hardcoded path before, need to check this works

print path2scripts

def loadOptions():
	parser = argparse.ArgumentParser(description='Parse a PGM BAM, generating alignments of a random subset of the reads.') 
	parser.add_argument("-b", "--bam", dest="bam", required=True, help="The input BAM file")
	parser.add_argument("-r", "--reference", dest="ref", required=True, help="The reference FASTA file")
	parser.add_argument("-n", "--name", dest="name", required=True, help="The dataset prefix name")
	parser.add_argument("-o", "--output-dir", dest="dir", required=True, help="The top output directory")
	parser.add_argument("-m", "--num-reads-max", dest="read_max", required=False, help="The maximum number of reads to process")
	parser.add_argument("-s", "--subset-size", dest="sub_size", required=False, help="The number of reads per subset (max)")
	parser.add_argument("-p", "--pass-sff", dest="skip_bam", action="store_true", required=False, default=False, help="Skip the BAM parsing bit (files already generated)")
	parser.add_argument("-t", "--terminate", dest="term_before_align", action="store_true", required=False, default=False, help="Terminate before performing sequence indexing and alignment")
	return parser

# p is an argument parser.
def processOptions(p): 
	userOpts = vars(p.parse_args())
	try:
		with open(userOpts['bam']): pass  # open file
	except IOError as e:
		print "Unable to open file " + str(e)  # Does not exist OR no read permissions
		p.print_help()
	return userOpts

def generate_formatted_files_from_bam(read_obj, fasta, qual, flow, log):
	#the FASTA is in RLE format.	 
	rle_seq = "" #string build it from iterating over read_obj
	base_pos = 0 #including the key.
	rle_pos = 0
	
	for used_flow_pos in read_obj.used_fp:
		base_call = read_obj.used_fp_to_bc[used_flow_pos]
		rle_seq = rle_seq + base_call[0]
		valid_fv = False if used_flow_pos in read_obj.used_fp_to_valid_fv else True 
		called_len = len(base_call)
		#same length as the called sequence...
		for index in xrange(len(read_obj.qual[used_flow_pos])):
			qual_val = read_obj.qual[used_flow_pos][index]				
			qual.write("\t".join([read_obj.id, str(base_pos), str(qual_val)]) + "\n")
			if index == 0: #only output for the first position.
				flow.write("\t".join([read_obj.id, str(used_flow_pos), str(rle_pos), str(base_pos), base_call[0], 
										str(read_obj.used_fp_to_fv[used_flow_pos]), str(valid_fv), str(read_obj.out_of_phase), str(called_len)]) + "\n")
			base_pos = base_pos + 1
		rle_pos = rle_pos + 1
			
	#flow file has flow calls for zero flows, but we are no longer interested in pursuing those (they are probably made up anyhow).	
	record = SeqRecord(Seq(rle_seq, IUPAC.Alphabet.DNAAlphabet),id=read_obj.id)
	SeqIO.write(record, fasta, "fasta")
	
	#adapter information should come from the BAM Parser.
def create_segemehl_database(output_dir, reference_rle,log):	
	log.write("Creating segemehl database") 
	segemehl_ref_db = output_dir + "reference.idx"
	
	#print "ref db %s reference rle %s" % (segemehl_ref_db, reference_rle)
	
	if(not os.path.isfile(segemehl_ref_db)):
		call([segemehl_exe, "-x", segemehl_ref_db, "-d", reference_rle])
		
	return segemehl_ref_db
	

def run_segemehl(read_rle_file, ref_rle_file, ref_index_file, output_file_mapped, output_file_unmapped, log):
	#print "Run command %s %s %s %s %s %s " % (segemehl_exe, ref_index_file, ref_rle_file, read_rle_file, output_file_mapped, output_file_unmapped) 
	infile = open(read_rle_file)
	
	read_rle_cleaned = read_rle_file + ".bk"
	
	outfile = open(read_rle_cleaned, 'w') 
	
	for line in infile:
		line = line.replace(" <unknown description>", "")
		outfile.write(line)
	infile.close()
	outfile.close()

	if(not os.path.isfile(output_file_mapped)):
		call([segemehl_exe, "-i", ref_index_file, "-d", ref_rle_file, "-q", read_rle_cleaned, "-o", output_file_mapped, "-u", output_file_unmapped, "-A", "85"])
	
	if(not os.path.isfile(output_file_mapped)):
		print "The output file does not seem to exist %s" % (output_file_mapped)
	
	
	
def parse_segemehl(align_out, read_rle, reference_rle, prefix, log):
	log.write("Parsing segemehl output")
	#print "perl %sparse_segemehl_output_final.pl %s %s %s %s false" % (path2scripts, align_out, read_rle, reference_rle, prefix)	
	perl_comm = 'perl %s/parse_segemehl_output_final.pl %s %s %s %s %s' % (path2scripts, align_out, read_rle, reference_rle, prefix, 'false')
	cmd = shlex.split(perl_comm)
	proc = Popen(cmd)
	

def prepare_reference_rle(output_dir, dataset_name, reference_file, log):
	log.write("Preparing reference RLE sequence")
	reference_dir = "%sreference_rle_for_%s/" % (output_dir,dataset_name)
	
	if not os.path.isdir(reference_dir):
		call(["mkdir", reference_dir])
	reference_prefix = reference_dir + "reference";
	
	#expected output files
	output_ref_rle = reference_prefix + "_rle.fasta";
	output_ref_hp_lengths = reference_prefix + "_rle.sql";

	retVal = 0
	#overall reference
	if(not os.path.isfile(output_ref_rle) or os.stat(output_ref_rle).st_size < 10): #something went wrong
		path_perl = os.path.abspath("/usr/bin/perl")
		path_script =  os.path.abspath(path2scripts + "generate_rle_seq_and_homopolymer_length_file_for_reference_sequence.pl")
		ref_file = os.path.abspath(reference_file)
		ref_prefix = os.path.abspath(reference_prefix)
		
		retVal = call([path_perl, path_script, ref_file, ref_prefix]);

	if(retVal != 0):
		log.write("Error: Failed to finish generate rle seq and hp length for reference\n");
		log.close()
		sys.exit(1)
	return [output_ref_rle, output_ref_hp_lengths] #return the file names for the reference fasta and sql


### MAIN ###
p = loadOptions()
opts = processOptions(p)

output_dir = opts['dir']
if output_dir[-1] != "/":
	output_dir = output_dir + "/"


if not os.path.exists(output_dir):
	os.makedirs(output_dir)

ref_file = opts['ref']

if not os.path.isfile(ref_file):
	raise "Reference file does not exist"

bam_file = opts['bam']

if not os.path.isfile(bam_file):
	raise "BAM file does not exist"

data_name = opts['name']

log = open(output_dir + "/" + data_name + ".log", "w")


### We are through the parameter checking.
# Need to produce FAST out, QUAL out, FLOW out, adapter locations (if any)

num_reads = int(check_output([os.path.abspath("/usr/bin/samtools"), "view", "-c", os.path.abspath(bam_file)])) # this does not make sense.

print "Num reads %d" % num_reads

read_max = -1
sub_size = -1

if opts['read_max']:
	read_max = int(opts['read_max'])
else:
	read_max = num_reads

if opts['sub_size']:
	sub_size = int(opts['sub_size'])
else:
	sub_size = 100000
	
#check that the outputs are available.
it = ITBAMIterator(bam_file, False) #False is for debug

num_subsets = 1 if read_max < sub_size else (read_max / sub_size) 

## this has an error when there is only 1 sequence.

prob_for_not_used = (num_reads - read_max) / float(num_reads)
prob_for_each_subset = ((1 - prob_for_not_used) / float(num_subsets))

print "Prob for not used portion %0.4f, prob for each subset: %0.4f" % (prob_for_not_used, prob_for_each_subset)

lower_bounds = [prob_for_not_used]

read_subset_dir_names = ["%sread_dir_%d_for_%s" % (output_dir,i,data_name)  for i in xrange(num_subsets)]
print read_subset_dir_names


if not opts['skip_bam']:
	#create the output directories
	[os.mkdir(d) if not os.path.isdir(d) else '' for d in read_subset_dir_names]
	
	#create the output files
	
	fasta_files = [open("%sread_dir_%d_for_%s/reads.fasta" % (output_dir,i,data_name), 'w') for i in xrange(num_subsets)]
	qual_files = [open("%sread_dir_%d_for_%s/reads.qual" % (output_dir,i,data_name), 'w') for i in xrange(num_subsets)]
	flow_files = [open("%sread_dir_%d_for_%s/reads.flow" % (output_dir,i,data_name), 'w') for i in xrange(num_subsets)]
	
	for f in qual_files:
		f.write("\t".join(["SEQUENCE","BASE_POSITION","QUALITY"]) + "\n")
		
	for f in flow_files:
		f.write("\t".join(["READ_ID", "FLOW_POSITION", "RLE_POSITION", "READ_POSITION", "NUCLEOTIDE", "FLOW_VALUE", "VALID_FV", "OUT_OF_PHASE", "CALLED_LEN"]) + "\n")
	
	for i in xrange(num_subsets):
		lower_bounds.append(prob_for_not_used + ((i + 1) * prob_for_each_subset))
	
	print lower_bounds
	reads_processed = 0
	
	## We can consider processing the entire file.
	## Wonder  whether it's worth it though?
	
	while True:
		try:
			read_object = it.get_next_read_obj()
			groupr = random.random()
			groupa = None
			
			if reads_processed >= read_max:
				break
			for i in xrange(len(lower_bounds)):
				if groupr < lower_bounds[i]:
						groupa = i
						break
			if(groupa == None):
				groupa = len(lower_bounds) - 1
				
			if(groupa > 0):
				groupa = groupa - 1 #to map it to the files indexes.
				generate_formatted_files_from_bam(read_object, fasta_files[groupa], qual_files[groupa], flow_files[groupa], log) # Error here
				reads_processed = reads_processed + 1
			else:
				pass #not used	 
		except StopIteration:
			break
	
	#all the files are complete.
	for i in xrange(num_subsets):
		fasta_files[i].close()
		qual_files[i].close()
		flow_files[i].close()	 


if opts['term_before_align']:
	sys.exit(0)


(output_ref_rle, output_ref_hp_len) = prepare_reference_rle(output_dir, data_name, ref_file, log)
segemehl_ref_index = create_segemehl_database(output_dir, output_ref_rle, log)

for i in xrange(len(read_subset_dir_names)):
	output_prefix = read_subset_dir_names[i]
	output_file_mapped = output_prefix + "/reads_versus_ref.map"
	output_file_unmapped = output_prefix + "/reads_not_mapping.txt"
	fasta_file = output_prefix + "/reads.fasta"
	run_segemehl(fasta_file, output_ref_rle, segemehl_ref_index, output_file_mapped, output_file_unmapped, log)
	parsed_prefix = output_prefix + "/" + data_name
	parse_segemehl(output_file_mapped, fasta_file, output_ref_rle, parsed_prefix, log)

		
if __name__ == '__main__':
	pass
