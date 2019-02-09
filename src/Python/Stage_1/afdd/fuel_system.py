#====================================================================
# MODIFIED FUEL SYSTEM MODEL BY AS (tweak of NDARC)
#
# Tanks and plumbing for liquid fuel based
#====================================================================

import sys
from conversions import *
f_cw          = 1.3131  # 1.3131 for crashworthiness
f_bt          = 1.8     # for ballistic survivability
#k0_plumb      = 120    # constant weight; default
#k1_plumb      = 3      # proportional weight; default
nint          = 1       # number of internal tanks
fuel_scaling  = 1.42    # allow extra fuel capacity for (a) unused (b) mission flexibility
density       = 6.5     # fuel density in lb/gallons; value 6.7 is for jet fuel

#====================================================================
# begin routine
#====================================================================

def fuel_system(vehicle_parameters):

#====================================================================
# unpack input dictionary
#====================================================================

   GTOW             = vehicle_parameters['gtow'] 
   flow_rate        = vehicle_parameters['flow_rate']
   w_fuel           = vehicle_parameters['fuel']
   Nengine          = vehicle_parameters['nengine']
   fac              = vehicle_parameters['tech_factors'].fuel_system
   tank_vol         = fuel_scaling * w_fuel / density           # volume, gallons

#====================================================================
# From tank volumne, find tank weight from statistical fit
#====================================================================

#   print 'fuel mass is ',m_fuel,' lbs'
#   print 'VOLUME IS ',tank_vol, 'gallons'
   wght_tank  = ( 0.4341 * tank_vol**0.7717 * 
                  nint**0.5897 * f_cw * f_bt**1.9491 )

#====================================================================
# find plumbing system weight from flow rate equation
#====================================================================

   k0_plumb   = 0.011*GTOW               # constant part: % of weight (was 120 lb)
   if k0_plumb > 120:
      k0_plumb  = 120.0
   k1_plumb   = k0_plumb/40.0

   gps_engine = flow_rate/float(Nengine)#*1.3    # 30% margin so pumps dont fail at max capacity
   wght_plumb = k0_plumb + k1_plumb *                               \
                (0.01*nint + 0.06*Nengine) * (gps_engine**0.866 )

#   print k0_plumb,k1_plumb, Nengine, gps_engine

#====================================================================
# apply tech factor, calculate total 
#====================================================================

   wght_plumb = wght_plumb*fac 
   wght_tank  = wght_tank *fac 

   total      = wght_plumb + wght_tank

#====================================================================
# create dictionary with outputs and return it
#====================================================================
  
   fuel_sys  = {'plumbing': wght_plumb*lb2kg, 'tank': wght_tank*lb2kg, 
                   'total': total*lb2kg}

#====================================================================
# if > 0.5 fuel weight, saturate
#====================================================================

#   print total, w_fuel
   if total > 0.5*w_fuel:
      fuel_sys    = {'total': 0.5*w_fuel*fac*lb2kg} 

   return fuel_sys