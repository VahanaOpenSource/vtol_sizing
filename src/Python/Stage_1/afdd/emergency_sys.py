#====================================================================
# Weight model for emergency systems
#====================================================================

from conversions import *
def emergency_sys(mass, tech_factor):

#====================================================================
# parametric model: fraction of GTOW
#====================================================================

   mass_brs       = 0.0233*mass*tech_factor
   return mass_brs