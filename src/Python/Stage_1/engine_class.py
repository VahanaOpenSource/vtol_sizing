#====================================================================
# Fuel-burning engine class that remembers parameters for sizing
#====================================================================

class all_engines:

   def __init__(self, data, nseg):
      self.ngroups      = 0
      self.groups       = {}

      ngrp              = 0
#loop over engine groups, find design parameters
      for key in sorted(data):
         self.groups[ngrp]  = engine_group(data[key], key, nseg)
         ngrp               = ngrp + 1

      self.ngroups      = ngrp
      self.nengines     = 0       # total number of motors
      return None

#====================================================================
# Motor group details
#====================================================================

class engine_group:

#====================================================================
# function to initialize all segments
#====================================================================

   def __init__(self, data, key, nseg):

#=======================================================================
#geometric parameters
#=======================================================================

      self.key                = key 
      self.rotor_group_ids    = []      # which rotor groups are driven by the engine
      self.motor_group_ids    = []
      self.p_ins              = numpy.zeros(nseg)
      self.nengines           = 0
