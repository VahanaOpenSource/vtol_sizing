#====================================================================
# Ananth Sridharan, Jan 24th 2017
#
# AIRCRAFT WING GROUP: 29-1.2 SECTION IN NDARC THEORY MANUAL V1_11.PDF, pg 249
#
# AFDD 93 Model implemented by AS
# 
#====================================================================

f_LGloc     = 1.0                # 1.7247 if landing gear on wing, 1.0 otherwise
k_elect     = 0.5 # proportionality constants
k_rotor     = 0.5
from conversions import *
def fixed_wing_wt(vehicle_parameters):

   Wt       = vehicle_parameters['gtow']        # in lbs
   nwings   = vehicle_parameters['nwing']       # number of fixed wings
   fac      = vehicle_parameters['tech_factors'].wing
   wing     = vehicle_parameters['wing']
   nr       = vehicle_parameters['nrotor']
   P        = vehicle_parameters['pwr_installed']*hp2kw
   l_fus    = vehicle_parameters['l_fus']*f2m
   b_fus    = vehicle_parameters['b_fus']*f2m
   red      = vehicle_parameters['wt_redund']

#====================================================================
# empirical parameters
#====================================================================

   tau_w    = 0.158                             # wing t/c (NOT TILTROTOR WING)
   Lamda    = 0.0                               # no sweep
   lamda    = 1.0                               # taper ratio
   nz       = 3.8                               # load factor 
   b_fold   = 0.0                               # fraction of span that folds (0 = not folding)

#====================================================================
# calculate weight of power wires and signal wires
#====================================================================
   
   mass_sw     = 3.04/17.5*(l_fus)*red['wires']*10

#====================================================================
#battery to HVDB
#====================================================================

   mass_pw     = 1.75*wing.nwings*l_fus/5.8*P/304*red['wires']

#   print(mass_pw);quit('OK?')
#====================================================================
# expression for weight per wing, lbs
# then multiply by # wings, sum over groups
#====================================================================

   wing_wt     = 0.0
   actuator_wt = 0.0 
   area        = 0.0
   act_cost    = 0.0
   str_cost    = 0.0
   # print(P);quit('ok?')
   # print(wing.ngroups)
   if nwings > 0:
      for i in range(wing.ngroups):

         group       = wing.groups[i]                 
         nwing       = group.nwings                   # number of wings
         Aw          = group.aspectratio              # aspect ratio 
         Sw          = group.area * m2f * m2f         # in sq.ft 
         fL          = group.lift_frac                # wing lift fraction for this group (0 to 1)
         W           = fL*Wt/nwing                    # thrust carried by each wing, lbs; this is both wings together!
         
         wt          = 5.66411*f_LGloc* (W*0.001/cos(Lamda))**0.847 * (nz**0.3958) * (Sw**0.21754)        \
                       *sqrt(Aw)* ((1.0+lamda)/tau_w)**0.09359 #* (1.0 - b_fold)**(-0.14356)
         wt          = wt * nwing
         wing_wt     = wing_wt + wt
         area        = area + Sw*nwing             # in sq.m
         group.wt    = wt*lb2kg/nwing

#====================================================================
# Add signal wire and power wire weights here; ideal place because we
# loop over wings in this function, so piggyback. 
#====================================================================

         mass_sw     = mass_sw + group.span*nwing*3.04/17.5
         mass_pw     = mass_pw + nwing*8.5*P/304*group.span/5.8

         group.wire  = group.span*(3.04/17.5 + 8.5*P/(304*5.8)) # wires in each wing, kg

#====================================================================
# for control surfaces, take total area and total weight; scaling law adapted by AS
# Note: right now redundancy is there in cost but not weight!
#====================================================================

      actuator_wt    = actuator_wt + 0.01735*(Wt**0.6435)*(area**0.40952)

#====================================================================
# weight effect of redundany actuators
#====================================================================

   actuator_wt       = actuator_wt * red['wing_flap']

#====================================================================
# tilt actuators
#====================================================================
      
   tilt_wt        = 0.1*wing_wt*red['tilt_actuator']

#====================================================================
# apply tech factor, take total weight of wing+actuators on wing
#====================================================================

   actuator_wt       = actuator_wt * fac
   wing_wt           = wing_wt     * fac
   tilt_wt           = tilt_wt     * fac 
   total             = actuator_wt + wing_wt + tilt_wt 

#====================================================================
#convert to kg, return dict for wing structure, actuators and wires
#====================================================================

   wt       = {'structure': wing_wt*lb2kg,'actuators':actuator_wt*lb2kg, 
               'total': total*lb2kg, 'tilters': tilt_wt*lb2kg}

   wires    = {'signal_wire': mass_sw, 'power_wire': mass_pw, 'total': mass_sw + mass_pw} 

   return wt, wires