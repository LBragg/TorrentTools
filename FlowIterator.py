'''
Created on 14/11/2013

@author: bra427
'''

class FlowIterator(object):
    '''
    classdocs
    '''
    fo, fv = [None] * 2
    flow_pos = 0

    def __init__(self, fo, fv):
        self.fo = fo
        self.fv = fv
        
    def increment(self):
        self.flow_pos+= 1
    
    def move_to_pos(self, fp):
        self.flow_pos = fp
    
    def get_prev_positive_flow_pos(self):
        if(self.flow_pos - 1 >= 0):
           for tmp_fp in range(self.flow_pos - 1, 0):
               if self.get_corrected_value_for_fp(tmp_fp) > 0:
                   return tmp_fp
           else:
               return None 
        return None
       
            
    def get_next_positive_flow_pos(self):    
       if(self.flow_pos + 1 < len(self.fv)):
           for tmp_fp in range(self.flow_pos + 1, len(self.fv)):
               if self.get_corrected_value_for_fp(tmp_fp) > 0:
                   return tmp_fp
           else:
              return None
       return None
               
                   
                   
    def get_corrected_value_for_fp(self, flow_pos):
        cv = self.fv[flow_pos] / float(256)
        if (cv < 0): cv = 0
        return cv
    
    def get_curr_flow(self):
       return self.flow_pos
    
    def get_curr_corrected_value(self):
        return self.get_corrected_value_for_fp(self.flow_pos)
    
    def is_finished(self):
        return self.flow_pos == len(self.fv)
   
    def has_more_flows(self):
        return self.flow_pos + 1 < len(self.fv) 
    def get_base_for_flow(self, flow_pos):
        assert(flow_pos < len(self.fo))
        return self.fo[flow_pos]
    
    def get_curr_base(self):
        return self.fo[self.flow_pos]
                       
