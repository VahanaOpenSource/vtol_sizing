#====================================================================
# Motor class that remembers parameters for sizing
#====================================================================

class motors:

   def __init__(self, data={}, nseg=0):
      self.ngroups      = 0
      self.groups       = {}

      ngrp              = 0
#loop over motor groups, find design parameters
      
      if(bool(data)):
         for key in sorted(data):
            self.groups[ngrp]  = motor_group(data[key], key, nseg)
            ngrp               = ngrp + 1

      self.ngroups      = ngrp
      self.nmotors      = 0       # total number of motors
      return None

#====================================================================
# Motor group details
#====================================================================

class motor_group:

#====================================================================
# function to initialize all segments
#====================================================================

   def __init__(self, data, key, nseg):

#=======================================================================
# geometric parameters
#=======================================================================

      self.cruise_efficiency  = data['cruise_efficiency']
      self.hover_efficiency   = data['hover_efficiency']
      self.key                = key 
      self.rotor_group_ids    = []      # which rotor groups are used to size this motor?
      self.p_ins              = numpy.zeros(nseg)
      self.nmotors            = 0
