#====================================================================
# Cruise model of the aircraft
# Performs basic sizing values and obtains the power required
#====================================================================
import hydra

import sys, os
from numpy import arctan2,cos,pi,sin,sqrt
#from conversions         import *

#====================================================================

class _cruisemodel:
  def cruisemodel(self, segID, update, use_bemt):
 
   icontinue = 1

#====================================================================
# assign pointer shortcuts
#====================================================================

   mission        = self.mission 
   wing           = self.wing 
   rotor          = self.rotor 
   motor          = self.motor 
   prop           = self.prop 
   emp            = self.emp_data
   masstakeoff    = self.massTakeoff
   aircraft       = self.all_dict['aircraft']
   std            = self.constants

   bfus           = emp.Geometry.fuselage_width
   
#====================================================================
# Aircraft properties
#====================================================================

   tipMachMax     = 0.9
   nPropeller     = aircraft['npropeller'] 
   if(nPropeller > 0):
      quit('BANG: add cruise propeller and associated motor group')
   effPropeller   = prop.eta 

#====================================================================
# Universal constants
#====================================================================

   grav           = std.grav                 # m/s^2
   kg2lb          = std.kg2lb
   f2m            = std.f2m   
   kts2mps        = std.kts2mps

#====================================================================
# Segment information
#====================================================================

   iseg           = segID - 1
   segment        = mission.segment[iseg]
   ainf           = segment.ainf 
   VtipMax        = tipMachMax * ainf           # max adv. tip speed wrt still air

#====================================================================
# define hover mass based on value at end of previous segment
#====================================================================

   if (segID > 1):
      massSegmentPrev      = mission.segment[iseg-1].mass
   else:
      massSegmentPrev      = masstakeoff 

#====================================================================
# subtract half the fuel of the segment (from previous iteration)
# to calculate thrust and weight levels (avg. over segment)
#====================================================================

   massSegmentPrev         = massSegmentPrev - 0.5*segment.fuel_mass 

#====================================================================
# Error trap for numbers that dont make sense (e..g NaN or negative mass)
#====================================================================

   if massSegmentPrev < 0:
      print ('WARNING: initial guess was OFF! DUH..',iseg,massSegmentPrev, segment.fuel_mass)
      massSegmentPrev = 400.0

#   if isinf(massSegmentPrev):
#      massSegmentPrev = 1e5
#      quit ('ERROR: infinity triggered for mass')

#====================================================================
# set local atmospheric conditions
#====================================================================

   rho                  = segment.rho 
   Vcruise              = segment.cruisespeed*kts2mps   # in m/s
   fsdyn                = 0.5 * segment.rho * Vcruise*Vcruise
   W                    = massSegmentPrev*grav

#====================================================================
# Size wings (if required) + calculate performance
#====================================================================

   Wings_lift, Wings_f  = self.wing.cruise_sizing(W, segment, False, bfus)
   Wing_drag            = fsdyn*Wings_f*1.05             # add 20% interference

#====================================================================
# Find drag of all edgewise rotors (non-prop-rotors)
#====================================================================

   Blade_drag           = self.rotor.blade_drag(Vcruise, VtipMax, rho)

#====================================================================
# find parasitic drag of frame, hubs and spinners
#====================================================================
   
   self.flat_plate_wrapper(segment, update)
   Fus_drag             = fsdyn*self.f_plate

#====================================================================
# Compute total drag from all wings, edgewise rotors and fuselage
#====================================================================

   Drag                 = Fus_drag + Wing_drag + Blade_drag
   FxRotor              = Drag                  # drag force

#====================================================================
# Remember cruise lift-to-drag ratio for output
#====================================================================

   if(update):
      self.LbyD         = W/Drag

#====================================================================
# calculate rotor vertical thrust required in cruise
# and catch error for negative thrust
#====================================================================

   FzRotor              = W - Wings_lift
   if FzRotor < -0.1*massSegmentPrev*grav:
      print('Rotor vertical force gone!',FzRotor, Wings_lift,massSegmentPrev*grav)
      quit('BANG: LOGIC FAIL IN CRUISE')

#====================================================================
# include effect of climb on shaft tilt
#====================================================================

#   alpha             = alpha + gamma_FP

#====================================================================
# Compute total rotor thrust and effective rotor plane tilt angle
# Assume equal thrust distribution on all rotors
#====================================================================

   NR_total          = self.rotor.nrotors
   Total_thrust      = sqrt(FxRotor*FxRotor + FzRotor*FzRotor)
   Total_thrust      = Total_thrust/NR_total

#====================================================================
# calculate rotor shaft power in cruise from all groups
#====================================================================
   
   ngroups              = rotor.ngroups
   P_shafts             = 0.0
   for i in range(ngroups):
      group             = rotor.groups[i]
      NR                = group.nrotors 
      rotor_fx          = FxRotor*group.Arat_thrust/NR
      rotor_fz          = FzRotor*group.Arat_lift/NR
      alpha             = arctan2(rotor_fz, rotor_fx) 

#====================================================================
# cruise propellers, or tilting rotors operating within +/- 10 deg 
# use propulsive efficiency parameter
#====================================================================

      if(group.type == 'tilting' or group.type == 'cruise') and alpha < 0.17:
         if use_bemt:
            eta         = group.rotor_aero_eta[iseg]
            quit('bang: need multi-rotor BEMT')
         else:
            eta         = effPropeller
         rotorPower     = rotor_fx*Vcruise/eta       # in watts
         rotorPower     = rotorPower*group.kint      # add rotor-rotor interference 

#====================================================================
# edgewise rotors (including tilting type rotors in not-quite-axial flow)
#====================================================================

      else:
         print(group.type,rotor_fx,rotor_fz,Wings_lift)
         quit('should not be here! cruise.py, edgewise rotors detected!')
         rotorPower     = hydra.rotorpower(Vcruise,  \
                             alpha, Total_thrust, cruise_tipspeed,rho,rotor.radius,rotor.solidity, rotor.cd0)

#====================================================================
# Remember rotor thrust, torque and power 
#====================================================================
      
      wing_ids          = group.wing_group_ids 
      omega_cruise      = group.tipspeed*group.RPM_ratio 
      mgid              = group.motor_group_id 
      eta_motor         = motor.groups[mgid].cruise_efficiency
      for wid in wing_ids:
         wing.groups[wid].rotor_thrust[iseg] = Total_thrust
         wing.groups[wid].rotor_power[iseg]  = rotorPower
         wing.groups[wid].motor_power[iseg]  = rotorPower/eta_motor     
         wing.groups[wid].rotor_torque[iseg] = rotorPower/omega_cruise

#====================================================================
# increment total shaft power
#====================================================================

      P_shafts          = P_shafts + rotorPower*NR 

#====================================================================
# Time spent in cruise, SI (min)
# Vcruise is in m/s and distance is in kms
#====================================================================

   if (segment.distance == 0.0):
      timeCruise             = segment.time
      segment.distance       = Vcruise * 1.e-3 * timeCruise * 60.0
   else:
      distance               = segment.distance 
      timeCruise             = (1000 * distance)/(Vcruise*60.0)
      segment.time           = timeCruise

   return icontinue