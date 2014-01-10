'''
Created on 14/11/2013
@author: bra427
'''
#! /usr/bin/python

import pysam, re
import argparse
import sys
import warnings
from FlowIterator import *
from BaseIterator import *

USED_FLOW_POSITIONS = "used_fp"
USED_FP_TO_BC = "used_fp2bc"
USED_FP_TO_FV = "used_fp2fv"
USED_FP_TO_VALID_FV = "used_fp2validfv"
LAST_PHASE_ERROR = "lpe"
CONSUMED_BASES = "consumed_bases"
OUT_OF_PHASE = "oop"


class AlignedFlowSeq(object):

    def __init__(self):
            self.id,self.qual, self.used_fp, self.used_fp_to_bc, self.used_fp_to_fv, self.used_fp_to_valid_fv, self.last_phase_error, self.consumed_bases, self.out_of_phase, self.five_prime_clip, self.three_prime_adapter_clip, self.insert_size = [None] * 12

class ITBAMIterator(object):
    bamfile = None
    debug = None
    tech, run_ID, key_seq,flow_order, flow_order_arr, trim_adapter = [None] * 6
    

    
    def __init__(self, bamfile, debug):
        self.debug = debug
        self.bamfile = pysam.Samfile(bamfile, "rb", check_sq=False)
        header = self.bamfile.header
        hd = header['HD']
        sq = header['SQ']
        rg_list = header['RG']
        pg = header['PG']  
           
        for my_dict in rg_list:
        #   print my_dict
            if("PL" in my_dict):
                self.tech = my_dict["PL"]
            if("KS" in my_dict):
                self.key_seq = my_dict["KS"]
            if("ID" in my_dict):
                self.run_ID = my_dict["ID"]
            if("FO" in my_dict):
                self.flow_order = my_dict["FO"]
        
        if self.tech != "IONTORRENT":
            raise Exception("Not IonTorrent Data: " + self.tech)

        self.flow_order_arr = list(self.flow_order)
        
        for my_dict in pg:
            if("CL" in my_dict):
                funCall = my_dict["CL"]
                if re.search("calibration-file", funCall):
                    raise Exception("Looks like the file was calibrated against a reference genome!: " + funCall); # don't want that for amplicon benchmarking.
                adapter_re = re.search("\-\-trim-adapter\s([ATGCN]+)", funCall)
                self.trim_adapter = adapter_re.group(1)
        
        
    def _getHeaderElement(self, key):
        if key in self.bamfile.header:
            return self.bamfile.header[key]

    def _processReadHelper(self, flow_it, base_it, last_phase_error, out_of_phase, fs, consumed_bases=0):  # will need to return a special result object, as integers are not mutable.

        # both flow_it and base_it are initialised to point to the position of the current unprocessed base.
        used_fp2fv = {}
        used_fp2bc = {}
        used_fp2validfv = {} 
        used_fp = []  # this is so I can change the validity of the last one.
        no_flow_ctr = -1  # this gets decremented every time we see a value that cannot be assigned a flow. The order of the flows are deterimed by 'used_fp'.
        outer_loop_condition = not (base_it.is_finished() or out_of_phase)   # for clarity

        # iterate through the flow-values
        while outer_loop_condition:  # need to check the iterators somewhere.
            
            fv_i = flow_it.get_curr_corrected_value()
            bases = int(fv_i + 0.5)
            lconsumed_bases = consumed_bases  # i made a local copy, I am going to modify it, just to make code clearer
            increment_flow = True
    
            if self.debug: print "Flow value %0.2f, curr flow position: %d, flow base: %s called_bases %d" % (fv_i, flow_it.get_curr_flow(), flow_it.get_curr_base(), bases)
            if(fv_i > 0): 
                for j in range(0, bases - consumed_bases):  # process the bases in the flow-call
                    if(flow_it.get_curr_base() == base_it.get_curr_base()):  # easy, perfect match
                        if self.debug: print "Perfect match"
                        if(lconsumed_bases > 0 and j == 0):  # we want to record that this flow_index will not match the flow-value.        
                            used_fp2validfv[flow_it.get_curr_flow()] = False  # this is not necessarily a phase correction either... if the carried over bases can only be from key sequence?
                            
                        
                        if(j == 0):    
                            used_fp.append(flow_it.get_curr_flow())
                            used_fp2fv[flow_it.get_curr_flow()] = fv_i
                    
                            if self.debug: print "Appending to flow: %d " % (flow_it.get_curr_flow())
                            
                        # actual bases are recorded to this flow.
                        if flow_it.get_curr_flow() in used_fp2bc :           
                            used_fp2bc[flow_it.get_curr_flow()] = used_fp2bc[flow_it.get_curr_flow()] + flow_it.get_curr_base()  # append
                        else:
                            if self.debug: print "Storing flow %d = %s" % (flow_it.get_curr_flow(), flow_it.get_curr_base())
                            used_fp2bc[flow_it.get_curr_flow()] = flow_it.get_curr_base()  # add this one
                        lconsumed_bases += 1  # number of bases used from this flow-call
                     
                        if self.debug: print "Base it position: %d nuc: %s " % (base_it.get_curr_base_pos(),base_it.get_curr_base())
                        base_it.increment()  # this base was processed.
                        
                        if base_it.is_finished():  # only break if it is finished.
                            if self.debug: print "Base iterator finished"
                            outer_loop_condition = False  # there are no more bases, and we've assigned *this* one      
                            # we are at the end of the template already... check whether there are more bases from this flow that were not included in seq.
                            if j + 1 < (bases - lconsumed_bases):
                                # this flow was truncated too, the flow-value is useless, but not-exactly a phase correction, so not an error in that sense.
                                used_fp2validfv[flow_it.get_curr_flow()] = False
                            break
                    else:  # # they do not match.
                        if self.debug: print("Base called but not matching template - flow %d, flow base called %s,  curr base in template %s, flow value %0.2f" % (flow_it.get_curr_flow(), flow_it.get_curr_base(), base_it.get_curr_base(), fv_i))  # this should not happen for the key
                        # all our special cases here.
                        # 1. This base equals the previous base, assume no corresponding flow-call.
                        next_pos_flow = flow_it.get_next_positive_flow_pos()  # get the flow position of the next 'positive' flow call
                        if(base_it.has_bases_before_curr_pos() and base_it.prev_base() == base_it.get_curr_base()):
                            if self.debug: print "Previous base in called sequence equals current base in called sequence"
                            # 1a. This next sequence base matches this flow, so it's just an extension of the previous call.
                            if(base_it.has_more_bases() and base_it.next_base() == flow_it.get_curr_base()):
                                # prev_base equals this base, so likely the flow-call needed to be extended
                                # the next base equals the current flow, so we shouldn't increment the flow position.
                                if self.debug: print "Next base in called sequence equals this flow"
                                increment_flow = False
                                # get the last positive call position, it will no longer match up with its FV.
                                
                                #TODO: causing error.
                                last_used_fp = used_fp[-1]
                                used_fp2bc[last_used_fp] = used_fp2bc[last_used_fp] + base_it.get_curr_base()
                                used_fp2validfv[last_used_fp] = False
                                # now since the next base equals this 'flow', we can move forward a base
                                base_it.increment()
                                # this was a phase error, so check whether there were 'too' many phase errors in this location
                                if last_phase_error != None and flow_it.get_curr_flow() - last_phase_error < 2:  # we've two errors in close proximity... we're out of phase. 
                                    out_of_phase = True
                                last_phase_error = flow_it.get_curr_flow()
                                break  # we didn't use this flow, so don't keep iterating through the flow-call bases.
                            # 1b. This base equals the last base, but the next base either doesn't exist, or doesn't line up with this flow.
                            else: 
                                if self.debug: print "Either no bases left after this OR next base in sequence does not match this flow"
                                last_used_fp = used_fp[-1]  # whats the last flow position we used.
                                
                                used_fp2bc[last_used_fp] = used_fp2bc[last_used_fp] + base_it.get_curr_base()  # add nuc
                                used_fp2fv[last_used_fp] = used_fp2fv[last_used_fp] + 1.0  # add one base                            
                                
                                # handle phase error business.
                                if last_phase_error != None and flow_it.get_curr_flow != last_phase_error and (flow_it.get_curr_flow() - last_phase_error) < 2:
                                    out_of_phase = True  
                                
                                last_phase_error = flow_it.get_curr_flow()
                                base_it.increment()
                                # special case.
                                if base_it.is_finished():
                                    if self.debug: print "No bases left"
                                    outer_loop_condition = False
                                    break
                        # 2. This base does not match the previous base, but it matches a future flow, assume this flow is an 'error'       
                        elif(next_pos_flow != None and flow_it.get_base_for_flow(next_pos_flow) == base_it.get_curr_base()):
                            if self.debug: print "This base matches a future flow, increment the flow not the base"
                            if (not last_phase_error == None) and flow_it.get_curr_flow() - last_phase_error < 2:
                                out_of_phase = True
                            last_phase_error = flow_it.get_curr_flow()
                            break 
                            # this flow-value is ignored (but recorded as phase error) and increment flow is still true.  
                        # 3. This base doesn't match this flow, doesn't match the previous base and there are no flows left
                        elif(next_pos_flow == None):       
                            # we've run out of positive flow values, and this base call does not match a previous base-call
                            if self.debug: print "No more +ve flow values, and this base does not match a previous base-call"
                            out_of_phase = True
                            last_phase_error = flow_it.get_curr_flow()
                            increment_flow = False
                        # 4.
                        else:
                            if self.debug: print "This base has no corresponding 'flow' must have been invented by IT"
                            # This base doesn't match the previous base
                            # This base doesn't match the next pos flow, and there is next pos flow.
                            # this base must be an insertion 'introduced' by the TS. It 'could' be multiple bases, though unlikely.
                            if last_phase_error != None and flow_it.get_curr_flow() - last_phase_error < 2:
                                out_of_phase = True
                                last_phase_error = flow_it.get_curr_flow()
                                break
                            else:  # my special case, an invented flow-value
                                last_phase_error = flow_it.get_curr_flow()
                                used_fp2fv[no_flow_ctr] = 1.0
                                used_fp2validfv[no_flow_ctr] = False  # this is not necessarily a phase-correction either... if the carried over bases are assumedly from the key?
                                used_fp.append(no_flow_ctr)
                                used_fp2bc[no_flow_ctr] = base_it.get_curr_base()
                                
                                if(base_it.has_more_bases):
                                    base_it.increment()
                                else:
                                    outer_loop_condition = False
    
                                no_flow_ctr -= 1  # decrement the counter
                                break    
            else:
                if self.debug: print "The flow-value is not +ve"
                pass  # the flow value is not positive.
           
            if lconsumed_bases != 0 and lconsumed_bases < bases:
                if self.debug: print "The number of consumed bases was less than the flow-call"
                if(last_phase_error != None and flow_it.get_curr_flow() - last_phase_error < 2):
                    out_of_phase = True
                last_phase_error = flow_it.get_curr_flow()
                # we did not consume all of the bases for this flow.
                # this is also a phase error.
                lconsumed_bases = 0  # we reset, cause we know the next base doesn't match.
         
            if out_of_phase:
                if self.debug: print "We are out of phase"
                # we could iterate through the rest of the read, and just assign each of the bases a truncated flow value.
                break
            # at the end, if the flow position should be incremented, it will be.     
     
            if flow_it.is_finished():  # this condition neesds to be fixed.
                break
      
            if increment_flow:
                if flow_it.has_more_flows():
                    flow_it.increment()
                    lconsumed_bases = 0  # new flow.
                else:
                    break
                
            if base_it.is_finished(): #need to make sure we increment the flow if it was used.
                break
      
        # outside the loop
        if not base_it.is_finished():
            while True:    
                if(base_it.get_curr_base() == base_it.prev_base()):
                    last_used_fp = used_fp[-1]
                    used_fp2bc[last_used_fp] = used_fp2bc[last_used_fp] + base_it.get_curr_base()  # add nuc
                    used_fp2fv[last_used_fp] = used_fp2fv[last_used_fp] + 1.0  # add one base
                else:
                    used_fp2fv[no_flow_ctr] = 1.0
                    used_fp2validfv[no_flow_ctr] = False  # this is not necessarily a phase correction either... if the carried over bases are assumedly from the key?
                    used_fp.append(no_flow_ctr)
                    used_fp2bc[no_flow_ctr] = base_it.get_curr_base()
                    no_flow_ctr -= 1  # decrement the counter
                # do stuff!
                if not base_it.has_more_bases():
                    break
                else:
                    base_it.increment()     
        fs.used_fp = used_fp
        fs.used_fp_to_bc = used_fp2bc
        fs.used_fp_to_fv = used_fp2fv
        fs.used_fp_to_valid_fv = used_fp2validfv
        fs.last_phase_error = last_phase_error
        fs.consumed_bases = consumed_bases
        fs.out_of_phase = out_of_phase
        return fs        

    # # this function is getting super-unwieldy.
    def _processRead(self, record):
        seq = record.seq

        flow_vals =  None
        aligned_fs = AlignedFlowSeq()
        aligned_fs.id = record.qname
        for key, value in record.tags:
            if(key == "ZA"):
                aligned_fs.three_prime_adapter_clip = value
            if(key == "ZG"):
                aligned_fs.insert_size = value
            if(key == "ZM"):
                flow_vals = value
            if(key == "ZF"):
                aligned_fs.five_prime_clip = value #after key
        flow_it = FlowIterator(self.flow_order, flow_vals)
       
        last_phase_error = None
        out_of_phase = False
        consumed_bases = 0
        comb_seq = self.key_seq + seq #rather do this instead.
        base_it_comb = BaseIterator(comb_seq)
       
        self._processReadHelper(flow_it, base_it_comb, last_phase_error, out_of_phase, aligned_fs, consumed_bases)
        #qual values should be a nested list
        aligned_fs.qual = self._qual_values_by_flow(record.qual, aligned_fs)
        return(aligned_fs)


    def _qual_values_by_flow(self, qual, aligned_flow_read_obj):
        pos_in_qual = 0
        qual_obj = {} #dictionary      
        pos_in_seq = 0
        
        verbose = False
        key_length = 4
        
        for fp in aligned_flow_read_obj.used_fp: #this is an array.  
            bc = aligned_flow_read_obj.used_fp_to_bc[fp] #these are the called bases.
            
            if(verbose):
                fv = aligned_flow_read_obj.used_fp_to_fv[fp]
                usable = not fp in aligned_flow_read_obj.used_fp_to_valid_fv
                print "1. Flow pos: %d base call: %s Flow value %s  Usable %s #Used flows %d" % (fp, bc, fv, usable, len(aligned_flow_read_obj.used_fp))
                
            qual_obj[fp] = []
            if(verbose):    
                print "2. Length qual: %d " % (len(qual))
            
            
            for nuc in xrange(len(bc)):
                
                if(pos_in_seq >= key_length and  fp >= aligned_flow_read_obj.five_prime_clip and pos_in_qual < len(qual)):
                    if(fp == 7 and nuc == 0): #special case, part of the key
                        qual_obj[fp].append("NA")
                        continue
                        
                    my_qual = ord(qual[pos_in_qual]) - 33
                    
                    if verbose:
                        print "3a. Qual pos: %d, QUAL: %d" % (pos_in_qual, my_qual)
                    
                    qual_obj[fp].append(my_qual) #add the quality for this flow-call to the array.
                    pos_in_qual = pos_in_qual + 1
                elif fp < 0 and pos_in_seq >= key_length:
                    my_qual = ord(qual[pos_in_qual]) - 33
                    if verbose:
                        print "3b. Qual pos: %d, QUAL: %d" % (pos_in_qual, my_qual)
                    
                    qual_obj[fp].append(my_qual)
                    pos_in_qual = pos_in_qual + 1
                    
                else:
                    qual_obj[fp].append("NA")
                pos_in_seq = pos_in_seq + len(bc)       
        if verbose:
            sys.exit(1)        
                
        return qual_obj
                    
    def get_next_read_obj(self):
        record = self.bamfile.next()
        #print record
        return self._processRead(record)
    
if __name__ == '__main__':
    pass
