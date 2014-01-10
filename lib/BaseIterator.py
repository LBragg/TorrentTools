'''
Created on 14/11/2013

@author: bra427
'''

class BaseIterator(object):
    '''
    classdocs
    '''
    seq, base_pos = [None] * 2

    def __init__(self, seq):
        '''
        Constructor
        '''
        self.seq = seq
        self.base_pos = 0
        
    def get_base_at_position(self, position):
        assert(position < len(self.seq))
        return self.seq[position]
    
    def increment(self):
        assert(self.base_pos < len(self.seq))
        self.base_pos += 1

    def prev_base(self):
       assert(self.base_pos - 1 >= 0)
       return self.seq[self.base_pos - 1]

    def next_base(self):
      assert(self.base_pos + 1 < len(self.seq))
      return self.seq[self.base_pos + 1]
   
    def get_curr_base(self):
       return self.seq[self.base_pos]
   
    def get_curr_base_pos(self):
        return self.base_pos   
   
    def is_finished(self):
        return  self.base_pos == len(self.seq)
    
    def has_more_bases(self):
        return self.base_pos + 1 < len(self.seq)
    
    def has_bases_before_curr_pos(self):
        return self.base_pos > 0
