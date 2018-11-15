#====================================================================
# Engine model (courtesy of Dylan Jude and Brent Mills from
#               2017 AHS Elysium design team)
#====================================================================

import os,sys
import numpy as np

from turboshaft         import turboshaft
from electric_motor     import electric_motor
from turbo_electric     import turbo_electric
from piston_electric    import piston_electric
from diesel             import diesel 
from stuttgart_diesel   import stuttgart_diesel
from piston             import piston 

#====================================================================
# part of hydra class (i.e., self refers to hydra)
#====================================================================

class _engines:
   
   # initialization
   def initEngine(self):
      
      etype = self.etype
   
      if etype == 'turboshaft':
         self.engine = turboshaft(self.emp_data)
      elif etype == 'hybrid':
         self.engine = hybrid(self.emp_data)
      elif etype == 'electric_motor':
         self.engine = electric_motor(self.emp_data)
      elif etype == 'turbo_electric':
         self.engine = turbo_electric(self.emp_data)
      elif etype == 'piston_electric':
         self.engine = piston_electric(self.emp_data)
      elif etype == 'piston':
         self.engine = piston(self.emp_data)
      elif etype == 'diesel':
         self.engine = diesel(self.emp_data)
      elif etype == 'stuttgart':
         self.engine = stuttgart_diesel(self.emp_data)
      else:
         print ('unknown engine type: what to do? SAVE ME ! CRITICAL ERROR')
         sys.exit(1)
 
#====================================================================
# wrapper function
#====================================================================

   def updateFuelAndSegmentWeight(self,segID):

      iseg = segID - 1

      mission      = self.mission                        # mission segment data structure
      aircraft     = self.all_dict['aircraft']
      fuel_tech    = self.emp_data.Tech_factors.Weight_scaling.fuel
      rotor        = self.rotor
      wing         = self.wing
      motor        = self.motor 
      eta_xmsn     = self.engine.eta_xmsn

#====================================================================
# segment properties
#====================================================================

      segment      = mission.segment[iseg]
      time         = segment.time          # time of segment, minutes
      flightMode   = segment.flightmode    # flight mode type (character string)
      theta        = segment.theta         # temperature ratio wrt MSL
      delta        = segment.delta         # pressure    ratio wrt MSL

#====================================================================
# Total energy used in this segment: use to size the battery 
#====================================================================

      P_rotors     = 0.0
      P_engine     = 0.0
      ngroups      = wing.ngroups
      for i in range(ngroups):
         group             = wing.groups[i]
         rotor_group_ids   = group.rotor_group_id 
         nrotors           = group.nrotors 

#====================================================================
# assume all rotors in this group operate at the same power setting
#====================================================================

         P_rotors          = P_rotors + group.rotor_power[iseg]*nrotors      
         P_engine          = P_engine + group.motor_power[iseg]*nrotors
      segment.p_req        = P_rotors*0.001
      segment.p_eng        = P_engine
      segment.energy       = P_engine*time*0.001/60.0        # kWh

#====================================================================
# rotor power required estimates: look at each rotor group
#====================================================================

      ngroups      = rotor.ngroups

      for i in range(ngroups):
         group     = rotor.groups[i]
         wgrp_id   = group.wing_group_ids          # all wing groups carrying this rotor group
         p_req     = 0.0

#====================================================================
# see max power required for a rotor from this group, based on all the 
# wings that carry a rotor from this group
# Note that this value is segment-specific so will keep changing
#====================================================================

         for iwg in wgrp_id:
            p_req  = max(p_req,wing.groups[iwg].rotor_power[iseg])

         p_req             = 0.001*p_req                       # kW, per rotor
         group.p_req       = p_req                             # power reqd per unit, for each rotor group

#====================================================================
# calculate power output required from the engine
# find motor efficiency
#====================================================================

         imotor            = group.motor_group_id
         if(flightMode == 'hover'):
            xmsn_eta       = motor.groups[imotor].hover_efficiency
         else:
            xmsn_eta       = motor.groups[imotor].cruise_efficiency 

#         group.powerEng    = p_req/self.engine.eta_xmsn       # power output required from source, kW
         group.powerEng    = p_req/xmsn_eta                    # power output required from source, kW
         group.powerTOP    = group.powerEng*group.Q_overload   # take-off power to install, kW (includes margins for safety/OMI)

#====================================================================
# Now loop over motor groups and see what max power output is reqd.
# to drive all the rotors associated with the group
#====================================================================

      ngroups         = motor.ngroups
      total_pins      = 0.0
      for i in range(ngroups):
         group        = motor.groups[i]
         rgrp_ids     = group.rotor_group_ids
         nmotor       = group.nmotors

#====================================================================
# Find maximum installed power for this motor, when its required to 
# power all of its associated rotor groups
#====================================================================

         p_ins        = 0.0
         for irg in rgrp_ids:
            p_ins           = max(p_ins,rotor.groups[irg].powerTOP)     # still from one rotor

         group.p_ins[iseg]  = p_ins 
         total_pins         = total_pins + p_ins*nmotor
      self.p_ins            = max(self.p_ins,total_pins)

#      print('installed power is ',self.p_ins)
         # print(p_ins);x=input('?')
#====================================================================
# find total installed power, SFC and fuel mass and flow rate of all
# engines
#====================================================================

      massFuel     = 0.0
      fuelFlowRate = 0.0

#====================================================================
# Apply tech factor for fuel weight: interpret as better SFC for 
# engine
#====================================================================

      massFuel     = massFuel * fuel_tech

#====================================================================
# update segment values for power, energy, sfc and fuel consumed
#====================================================================

      segment.sfc             = 1000.0/self.engine.sp_energy      # kg/kW-hr; not burned but used for calculations
      segment.fuel_mass       = massFuel               # fuel burn

#====================================================================
# update mission parameters
#====================================================================

      mission.totalenergyreqd += segment.energy
      mission.maxfuelflowrate = max(self.maxfuelflowrate,fuelFlowRate)
      
#====================================================================
# update mass segment
#====================================================================

      if (segID > 1):
         massSegmentPrev      = mission.segment[iseg-1].mass
      else:
         massSegmentPrev      = self.massTakeoff

#====================================================================
# added payload = negative of payload dropped; use a fractional value
# instead of absolute - relevant for sizing to  a given GTOW; note
# the negative sign between addPayload and frac_payload_drop 
#====================================================================

      addPayload =-mission.payload * mission.segment[iseg].frac_payload_drop

#====================================================================
#get mass at the end of the segment (for next segment power calcs) 
#====================================================================

      mission.segment[iseg].mass = massSegmentPrev - massFuel + addPayload

#====================================================================
# if weight changes over mission profile, remember this fact
# it will decide how many times we iterate for convergence 
#====================================================================

      if abs(addPayload - massFuel)>0.001:
         mission.same_weight     = False 

      return None

#====================================================================
# calculate engine weight
#====================================================================

   def getEngineWeight(self):

#====================================================================
# Setup quantities for each motor group
#====================================================================
      
      motor          = self.motor 
      ngroups        = motor.ngroups
      mission        = self.mission
      tech_factors   = self.emp_data.Tech_factors.Weight_scaling
      totalEnergy    = mission.totalenergyreqd     # in kWh
      fac            = tech_factors.powerplant     # motor technology factor
      battery        = tech_factors.battery        # battery technology factor
      mwts           = {}
      total          = 0.0
      for i in range(ngroups):
         group       = motor.groups[i]
         nmotor      = group.nmotors
         P_req       = np.amax(group.p_ins)

         inputs      = {'powerReq': P_req   , 'tech_fac': fac}
         motor_mass  = self.engine.getWeight(inputs)
         # print('rated power is ',P_req,'motor mass is ',motor_mass)
         motor_mass  = motor_mass*nmotor
         key         = 'motor_group' + str(i) 
         k2          = key + '_mount'
         mwts[key]   = motor_mass 
         mwts[k2]    = 0.15*motor_mass       # fraction based on Alpha

         total       = total+motor_mass+mwts[k2]

      mwts['total']  = total

#====================================================================
#setup input dictionary and calculate battery pack weight
#====================================================================

      inputs            = {'totalEnergy': totalEnergy,'tech_bat': battery , 'mission': mission}
      self.engine.battery_weight(inputs)
      # print(self.engine.m_batt);quit()
      return mwts