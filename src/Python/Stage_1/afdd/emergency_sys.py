#====================================================================
# Weight model for emergency systems
#====================================================================

from conversions import *
def emergency_sys(vehicle_parameters):

   Wt       = vehicle_parameters['gtow']        # in lbs

#====================================================================
# parametric model: fraction of GTOW
#====================================================================

   wt       = 0.0233*Wt*lb2kg
   return wt 