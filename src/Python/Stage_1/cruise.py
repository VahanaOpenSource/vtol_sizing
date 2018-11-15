#====================================================================
# Cruise model of the aircraft
# Performs basic sizing values and obtains the power required
#====================================================================
import hydra

from math                import atan2,isinf
from run_CSD             import run_CSD
from set_PRASADUM_inputs import set_linear_blade
from csd_rotor           import set_CSD_inputs
from flat_plate_buildup  import find_f
import sys, os
from numpy import arctan2,cos,pi,sin,sqrt
from conversions         import *
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

#====================================================================
# Aircraft properties
#====================================================================

   aircraftID     = aircraft['aircraftID'] 
   tipMachMax     = 0.9
   nPropeller     = aircraft['npropeller'] 
   if(nPropeller > 0):
      quit('BANG: add cruise propeller and associated motor group')
   effPropeller   = prop.eta 
   flatPlateFactor= aircraft['fdrag'] 

   ngroups        = wing.ngroups

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

   if (segID > 1 and not mission.span_driven):
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

   if isinf(massSegmentPrev):
      massSegmentPrev = 1e5
      quit ('ERROR: infinity triggered for mass')

#====================================================================
# set local atmospheric conditions
#====================================================================

   alt                     = segment.startalt
   temp                    = segment.deltatempisa
   density_alt             = segment.density_alt
   rho                     = segment.rho 
   
#====================================================================
# cruise speed in m/s, and mach number
#====================================================================
 
   Vcruise                 = segment.cruisespeed*kts2mps

#====================================================================
# climb angle (flight path angle gamma)
#====================================================================

   Vclimb                  = segment.rateofclimb/60.0 # convert to m/s
   if abs(Vclimb) > 1.e-3:
      gamma_FP             = atan2(Vclimb,Vcruise)
   else:
      gamma_FP             = 0.0

#====================================================================
# Wing-related parameters: size the wing based on the lift req.
#====================================================================
   
   fsdyn                = 0.5 * segment.rho * Vcruise*Vcruise
   Wings_lift           = 0.0
   Wings_f              = 0.0
   Blade_drag           = 0.0
   wing.nwings          = 0
   
   for i in range(ngroups):
      group             = wing.groups[i]
      loadingFrac       = group.lift_frac
      nwings            = group.nwings 
      wing_lift         = massSegmentPrev*grav*loadingFrac/nwings
      wing.nwings       = wing.nwings + nwings

#====================================================================
# obtain span and mean chord from wing cl
# also check for stall speed and set upper limit
#====================================================================

      V_ratio        = Vcruise/group.stall_speed
      CLmax          = group.CLmax*(V_ratio*V_ratio)

#====================================================================
# perform wing sizing if required
#====================================================================

      wing_cl        = min(group.cl, CLmax)
      if update:
         group.size_wing(wing_lift, fsdyn, wing_cl)

#====================================================================
# Non-sizing cruise flight condition
# Saturate wing cl at CLmax, calculate everything else
#====================================================================

      else:
         wing_lift   = group.calculate_lift(wing_lift, fsdyn, CLmax)

#====================================================================
# Compute drag coefficient of wing * area = equiv. flat plate area, sqm
#====================================================================

      fWing             = group.wing_f(wing_cl)
      Wings_f           = Wings_f + fWing*nwings 
      Wings_lift        = Wings_lift + wing_lift*nwings
#   print(massSegmentPrev,'lift of wings is ',Wings_lift)

#====================================================================
# now loop over all rotor groups, and calculate rotor drag
#====================================================================

   for i in range(rotor.ngroups):
      group             = rotor.groups[i]
      NR                = group.nrotors

#====================================================================
# edgewise rotor: get average profile drag coefficient of sections
#====================================================================

      cruise_tipspeed   = group.tipspeed*group.RPM_ratio

      if group.type == 'edgewise':
         if group.tipspeed + Vcruise > VtipMax:
            cruise_tipspeed   = VtipMax - Vcruise
            group.RPM_ratio   = cruise_tipspeed/group.tipspeed
         else:
            cruise_tipspeed   = group.tipspeed

         muCruise             = Vcruise/cruise_tipspeed
#         if muCruise > 1.0:
#            cd0               = 1.5*cd0
#         elif muCruise > 0.3:
#            cd0               = cd0*(1.0 + (muCruise-0.3)/0.7e0)
         cd0 = 0.012

#====================================================================
# continued for edgewise rotor: find blade drag in wind direction
#====================================================================

         D_blade        = group.solidity*cd0/8.0*(3.1*muCruise)
         D_blade        = D_blade * segment.rho * group.area * cruise_tipspeed**2
         D_blade        = D_blade * NR

         Blade_drag     = Blade_drag + D_blade 

#====================================================================
# tilting rotor: check helical tip mach number
#====================================================================

      elif group.type == 'tilting':

         Vmax               = sqrt(cruise_tipspeed*cruise_tipspeed + Vcruise*Vcruise)
         if(Vmax > VtipMax):
            cruise_tipspeed = sqrt(Vmax*Vmax - Vcruise*Vcruise)

         if cruise_tipspeed < 0.2*group.tipspeed:
            print ('warning: SLOWED ROTOR BELOW 20% RPM')
            print ('hitting min limit: 20% HOVER TIP SPEED')
            cruise_tipspeed = 0.2*group.tipspeed

#====================================================================
# continued: tilting rotor, get spinner + motor fairing drag
#====================================================================

         rspin           = 0.15                # spinner radius, meters
         Sf              = pi*rspin*rspin*NR   # frontal area of all spinners
         f_nac           = Sf*0.05             # drag coefficient with spinner head
         D_blade         = f_nac * fsdyn       # drag of hub
         Blade_drag      = Blade_drag + D_blade*0  # *0 is temp expression to match old code

#====================================================================
# unknown rotor type
#====================================================================

      else:
         quit('I dont know this rotor type')

#====================================================================
# AS: modified eqn. from JGL textbook: K.W^0.5 
# fparasitic in [m2]
#====================================================================

   if not update:
      fParasitic      = self.f_plate
   else:
      fParasitic      = flatPlateFactor*(masstakeoff*kg2lb*1.e-3)**0.5e0*f2m*f2m
      
#====================================================================
# estimate flat-plate area from component build-up 
#====================================================================

      if flatPlateFactor <= 1.e-3:
         inputs = {'W': masstakeoff*kg2lb,         # take-off weight in lbs 
                   'V': segment.cruisespeed,       # cruise speed, knots
                   'P': self.p_ins*kw2hp,          # installed power, hp
              'config': aircraftID,                # 1 = smr, 2 = tilt rotor, 3 = coaxial
               'Nprop': self.prop.num,             # # aux propellers
               'wings': self.wing,
           'gear_type': 'skids'}

         f_all, Vol     = find_f(inputs)
         self.Volume    = Vol
         f_all['ext']   = segment.add_f             # extra flat-plate area
         f_all['total'] = f_all['total'] + f_all['ext']
         # print f_all['total']
#         for key in f_all.keys():
#            print (key, f_all[key])#/f_all['total']*100
         # print 'updating flat plate area ',iseg, f_all['total']
         aircraft['f_breakdown']       = f_all

         fParasitic                    = f_all['total']*f2m*f2m

#====================================================================
# extra flat-plate area 
#====================================================================

      else:
         fParasitic     = fParasitic + segment.add_f*f2m*f2m

#====================================================================
# update flat plate factor and flat plate area for output
#====================================================================

      aircraft['flatplatefactor']   = flatPlateFactor
      self.f_plate                  = fParasitic

#====================================================================
# Compute total fuselage drag
#====================================================================

   Fus_drag       = 0.5*segment.rho*Vcruise**2*fParasitic

#====================================================================
# Get rotor drag due to blades (rotating + forward motion through space)
# and account for multiple edgewise moving rotors
#====================================================================

#====================================================================
# calculate rotor thrust in cruise
#====================================================================

   FzRotor           = massSegmentPrev*grav - Wings_lift
   if FzRotor < -0.1*massSegmentPrev*grav:
      print('Rotor vertical force gone!',FzRotor, Wings_lift,massSegmentPrev*grav)
      quit('BANG: LOGIC FAIL IN CRUISE')

#====================================================================
# Compute total drag from all wings 
# add interference factor of 1.05
#====================================================================

   Wing_drag         = fsdyn*Wings_f*1.05            # drag of all wings
   Drag              = Fus_drag + Wing_drag + Blade_drag

   if(update):
      self.LbyD      = massSegmentPrev*9.8/Drag

#====================================================================
# Propeller thrust calculation: use it to overcome all H-force
# Don't tilt rotor forwards when exclusive thrust propeller is present
#====================================================================

   if nPropeller > 0:
      prop.thrust    = Drag                  # Use propeller to counter drag
      FxRotor        = 0.0                   # when present; don't use rotor
      alpha          =-0.011                 # shaft angle relative to cruise flight direction
      prop.thrust    = prop.thrust - FzRotor*sin(alpha) 

   else:
      prop.thrust    = 0.0                   # Tilt rotor into wind to overcome
      FxRotor        = Drag                  # drag force
      alpha          = atan2(FxRotor, FzRotor)
   prop_power        = prop.thrust * Vcruise / effPropeller

#====================================================================
# < 15 deg tilt wrt vertical: add blade drag to horizontal force
#====================================================================

   if (abs(alpha) < 15.0*pi/180.0):
      print('rotor force demands are ',FxRotor,FzRotor,massSegmentPrev*grav)
      quit('BANG: cruise.py line 350')
#      muCruise    = Vcruise*cos(alpha)/cruise_tipspeed
#      Blade_drag  = Blade_drag + rotor.solidity*cd0/8.0*(3.1*muCruise)
#      Blade_drag  = Blade_drag * segment.rho * rotor.area * cruise_tipspeed**2
#      FxRotor     = FxRotor + Blade_drag*rotor.num

#====================================================================
# include effect of climb on shaft tilt
#====================================================================

   alpha             = alpha + gamma_FP

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

#====================================================================
# Tilting rotors and lifting wings: ideal!
#====================================================================

      if(group.type == 'tilting'and FzRotor <= 0.01*massSegmentPrev*grav):
         if use_bemt:
            eta         = group.rotor_aero_eta[iseg]
            quit('bang: need multi-rotor BEMT')
         else:
            eta         = effPropeller
         rotorPower     = Total_thrust*Vcruise/eta       # in watts
         rotorPower     = rotorPower*group.kint          # add rotor-rotor interference 

#====================================================================
# edgewise rotors (including tilting type rotors in not-quite-axial flow)
#====================================================================

      else:
         print(group.type,FzRotor,Wings_lift)
         quit('should not be here!')
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

      P_shafts          = P_shafts + rotorPower*NR_total 
      # print(Total_thrust,Vcruise,eta);quit()
#====================================================================
# add propeller shaft power too in case of dedicated cruise prop
#====================================================================

   P_shafts             =  P_shafts + prop_power 

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