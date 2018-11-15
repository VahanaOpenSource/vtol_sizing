#====================================================================
# Hover model of the aircraft
# Performs basic sizing values and obtains the power required
#====================================================================

import math
import os,sys
from numpy import pi, sqrt 

#====================================================================
# Function body
#====================================================================

class _hovermodel:
  def hovermodel(self, segID, update, use_bemt):

   icontinue = 1

#====================================================================
# assign pointer shortcuts
#====================================================================

   mission        = self.mission 
   aircraft       = self.all_dict['aircraft']
   engine         = self.emp_data.Engines
   motor          = self.motor 
   grav           = self.constants.grav
   wing           = self.wing
   rotor          = self.rotor 
   ngroups        = wing.ngroups
   aircraftID     = aircraft['aircraftID'] 
   do_sizing      = update
   effHoverPower  = engine.eff_hover_power 
   powerHoverAccs = engine.power_accs 
   self.engine.fracPowerMCP  = 1.0

   iseg           = segID - 1
   rho            = mission.segment[iseg].rho

#====================================================================
# get mass of system from end of previous segment
#====================================================================

   if (segID > 1 and not mission.span_driven):
      massSegmentPrev      = mission.segment[iseg-1].mass
   else:
      massSegmentPrev      = self.massTakeoff

#====================================================================
# subtract half the fuel of the segment from previous iteration
#====================================================================

   massSegmentPrev         = massSegmentPrev - 0.5*mission.segment[iseg].fuel_mass 

#====================================================================
# error trap: very large or very small value
#====================================================================

   if massSegmentPrev < 0 or massSegmentPrev > 1e9:
      massSegmentPrev = 1e9

   W           = massSegmentPrev*grav 
   Rmax        = self.constraints.max_rotor_radius
   geom        = self.emp_data.Geometry
   clearance   = geom.clearance           # rotor clearance as fraction of radius
   bfus        = geom.fuselage_width      # in meters 

#====================================================================
# loop over wing groups; associated with each wing group,
# we have a rotor group ID; find that rotor group 
#====================================================================

   for i in range(ngroups):
      wing_group     = wing.groups[i]
      rid            = wing_group.rotor_group_id
      rotor_group    = rotor.groups[rid]
      nrotors        = wing_group.nrotors       # left+right wings!
      lift_frac      = wing_group.lift_frac

      hvrDwldFactor  = rotor_group.hvr_dwld 

      if rotor_group.thrust_share == 'equal':
         thrust      = W*(1+hvrDwldFactor)/rotor.nrotors
         print('hacking hover thrust to be equal for all rotors')
      else:
         thrust      = W*lift_frac*(1.0+hvrDwldFactor)/nrotors       # thrust per rotor in the group

#====================================================================
# get motor group information (efficiencies)
#====================================================================

      mgid                    = rotor_group.motor_group_id 
      motor_group             = motor.groups[mgid]
      eta_motor               = motor_group.hover_efficiency 

#====================================================================
# perform sizing calculations if update is required
#====================================================================

#      print('segment # ',iseg,thrust)
      if do_sizing:                                   # perform sizing calculations
         rotor_group.hover_sizing(thrust, rho, Rmax, self.wing.groups, clearance, bfus)
         rotor_group.ainf        = mission.segment[iseg].ainf 

#====================================================================
# calculate thrust/torque/RPM margins
#====================================================================

         rotor_group.one_rotor_out(self.wing.groups)
         do_sizing               = False           # ensures that we wont run sizing several times unnecessarily

#====================================================================
# Compute the power required in hover (watts)
#====================================================================

      FM                   = wing_group.rotor_aero_eta[iseg]     # default 0.75; replaced by BEMT calc
      Phover               = rotor_group.hover_power(thrust, rho, FM, use_bemt, update)
      omega                = rotor_group.tipspeed/rotor_group.radius 

#====================================================================
# Total power in hover due to MR and accs.
#====================================================================

      rotor_power    = (Phover)*(1.0/effHoverPower + powerHoverAccs)
      engine_power   = rotor_power/eta_motor

#====================================================================
# power in kilowatts and torque in N-m; this is for a given rotor set
# belonging to a wing group, added over all rotors and all wings in 
# that group
#====================================================================

      wing_group.motor_power[iseg]    = engine_power         # power  per motor, watts 
      wing_group.rotor_power[iseg]    = rotor_power          # power  per rotor, watts
      wing_group.rotor_torque[iseg]   = Phover/omega         # torque per rotor, Nm
      wing_group.rotor_thrust[iseg]   = thrust

#      print(rotor_group.Q_overload,engine_power/1000)
#====================================================================
#adding tail rotor power for single MR type aircraft: 10% additional 
#====================================================================

   if (aircraftID == 1):
      mission.segment[iseg].p_req *= 1.1

   return icontinue