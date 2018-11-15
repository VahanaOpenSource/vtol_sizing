#====================================================================
#
# ALIGHTING GEAR GROUP
#
# Basic structure, retraction, crashworthiness structure
#
#====================================================================

#import sys
#f_lgret = 0.00       # retraction weight (0.08)
#f_lgcw  = 0.14       # crashworthiness weight
f_lg    = 0.0497      # landing gear weight fraction, alpha
f_fa 	= 0.015 	  # fairing weight / vehicle weight, fraction
W_S     = 1.0

from conversions import *
def alighting_weight(vehicle_parameters, bfile=None):

   gtow       = vehicle_parameters['gtow']      # in lbs
   fac        = vehicle_parameters['tech_factors'].landing_gear

#====================================================================
# fractional weight based on max GTOW:
#====================================================================

   mass_lg    = gtow*f_lg*fac*lb2kg
   mass_fair  = gtow*f_fa*lb2kg
   total 	  = mass_lg + mass_fair
   wght_lg 	  = {'gear_structure':mass_lg, 'fairing': mass_fair, 'total': total}
   return wght_lg

#====================================================================
# qbp
#====================================================================

#   if (aircraftID ==5):
#      num_lg = 4

#   wght_lg = 0.4013*gtow**0.6662*num_lg**0.536*W_S**0.1525

   #wght_lg    = f_lg*WMTO
#   wght_lgret = f_lgret * wght_lg
#   wght_lgcw  = f_lgcw  * (wght_lg + wght_lgret)
