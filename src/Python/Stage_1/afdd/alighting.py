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
f_fa 	  = 0.015 	  # fairing weight / vehicle weight, fraction

def alighting_weight(mass, tech_factor):

   """
   This function calculates the weight of the alighting gear group

   Input: 
   1. mass : vehicle maximum take-off mass, kgs 
   2. tech_factor: a scaling factor that enables modeling heavier/lighter system
   relative to the trendline
   """
#====================================================================
# fractional weight based on max GTOW:
#====================================================================

   mass_lg    = mass*f_lg*tech_factor
   mass_fair  = mass*f_fa
   total 	  = mass_lg + mass_fair
   wght_lg 	  = {'structure':mass_lg, 'fairing': mass_fair, 'total': total}
   return wght_lg

#====================================================================
# qbp
#====================================================================

# W_S     = 1.0
#   if (aircraftID ==5):
#      num_lg = 4

#   wght_lg = 0.4013*gtow**0.6662*num_lg**0.536*W_S**0.1525

   #wght_lg    = f_lg*WMTO
#   wght_lgret = f_lgret * wght_lg
#   wght_lgcw  = f_lgcw  * (wght_lg + wght_lgret)
