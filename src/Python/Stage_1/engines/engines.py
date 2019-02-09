#====================================================================
# Function to calculate fuel weight per segment and installed power req
#====================================================================

import os,sys
import numpy as np

#====================================================================
# part of hydra class (i.e., self refers to hydra)
#====================================================================

class _engines:
   
#====================================================================
# wrapper function
#====================================================================

   def updateFuelAndSegmentWeight(self,segID,i_priority):

      iseg = segID - 1

      mission        = self.mission                        # mission segment data structure
      aircraft       = self.all_dict['aircraft']
      fuel_tech      = self.emp_data.Tech_factors.Weight_scaling.fuel
      rotor          = self.rotor
      transmission   = self.transmission
      powerplant     = self.powerplant
      wing           = self.wing 
      fuselage       = self.fuselage 
      emp_data       = self.emp_data

#====================================================================
# segment properties
#====================================================================

      segment      = mission.segment[iseg]
      time         = segment.time          # time of segment, minutes
      flightMode   = segment.flightmode    # flight mode type (character string)
      theta        = segment.theta         # temperature ratio wrt MSL
      delta        = segment.delta         # pressure    ratio wrt MSL

#====================================================================
# Look at all the transmissions in the system, and identify the most 
# power-hungry/torque-hungry rotor group; use that value to size the 
# electric motor for an electric transmission
#====================================================================

      ngroups      = transmission.ngroups

      for i in range(ngroups):
         tgroup     = transmission.groups[i]
         
         if(tgroup.type == 'electric'):

#====================================================================
# loop over all the rotor groups powered by this electric transmission
# (defined as motor + wires, and optionally a generator)
#====================================================================

            for rgid in tgroup.rotor_group_ids:
               p_req     = 0.0
               q_req     = 0.0

               rgroup    = rotor.groups[rgid]

#====================================================================
# Examine the expected thrust/torque for the rotor, looking at all 
# support structures (wings, fuselage) that a rotor from this group 
# is mounted on
#====================================================================

               for iwg in rotor.groups[rgid].wing_group_ids:
                  p_req  = max(p_req, wing.groups[iwg].rotor_power[iseg] )
                  q_req  = max(q_req, wing.groups[iwg].rotor_torque[iseg])

#====================================================================
# check if this rotor is mounted on the fuselage too..
# if so, look at power and torque of those rotors to calculate max P,Q
#====================================================================

               if(rgroup.fuse_group_id != -1):
                  p_req  = max(p_req,fuselage.rotor_power[iseg])
                  q_req  = max(q_req,fuselage.rotor_torque[iseg])
   
            p_req             = 0.001*p_req           # kW, per rotor

            etainv      = 1.0/tgroup.eta
            if(segment.flightmode == 'hover'):
               etainv   = etainv/tgroup.motor_eta_hover
            else:
               etainv   = etainv/tgroup.motor_eta_cruise

#====================================================================
# set motor installed power, torque for this transmission
#====================================================================

            tgroup.motor_p_inp[iseg] = p_req * rgroup.Q_overload*etainv  # power  input to motor
            tgroup.motor_p_out[iseg] = p_req * rgroup.Q_overload         # power  output 
            tgroup.motor_q_out[iseg] = q_req * rgroup.Q_overload         # torque output

#====================================================================
# assume only one type of fuel used
#====================================================================

      Fuel           = 0.0
      TotalFlowRate  = 0.0 

#====================================================================
# loop over powerplants in the system
# identify linked transmissions, loop over those transmissions
# calculate total power output to all rotors, divide by xmsn efficiency
# also find the total fuel, flow rate and energy needed from this powerplant group
#====================================================================
      
      for pgid in range(powerplant.ngroups):
         pgroup    = powerplant.groups[pgid]

         P_rotors       = 0.0
         Q_rotors       = 0.0
         Energy         = 0.0 

         for xgid in pgroup.xmsn_group_ids:

            tgroup    = transmission.groups[xgid]
            etainv    = 1.0/tgroup.eta

            if(tgroup.type == 'electric'):
               if(segment.flightmode == 'hover'):
                  etainv   = etainv/tgroup.motor_eta_hover
               else:
                  etainv   = etainv/tgroup.motor_eta_cruise

            
            for irg in tgroup.rotor_group_ids:

#====================================================================
# adding tail rotor power for single MR type aircraft: 10% additional 
# Note that the torque is not modified, because the tail transmission
# is accounted for separately.
#====================================================================
               
               if (aircraft['aircraftID'] == 1 and segment.flightmode == 'hover'):
                  factor            = 1.1         
               else:
                  factor            = 1.0

               P_group        = factor*tgroup.power_req[ iseg,irg]*etainv
               Q_group        = factor*tgroup.torque_req[iseg,irg]*etainv

               P_rotors       = P_rotors + P_group
               Q_rotors       = Q_rotors + Q_group

#====================================================================
# then find rated power/torque for the powerplant and total energy 
# in this segment; this part specifically is the powerplant output power
# Note: loss_filter is a fraction representing how much power is lost 
# at engine output due to filter-related losses for the air intake
# powerHoverAccs is a fraction of the main rotor power diverted to running
# instrumentation and avionics
#====================================================================

         if(pgroup.type == 'turboshaft'):
            lossFilter        = emp_data.Engines.loss_filter
            powerHoverAccs    = emp_data.Engines.powerHoverAccs
         else:
            lossFilter        = 0.0 
            powerHoverAccs    = 0.0 

         P_plant              = P_rotors/(1.0-lossFilter)*(1.0+powerHoverAccs)
         Energy               = P_plant*time*0.001/60.0       # in kWh
         Q_plant              = Q_rotors

         # print(tgroup.power_req[ iseg,irg],P_plant)
         # print(Energy); input('ok energy for first segment?')
#====================================================================
# for turboshaft engines, account for losses here due to filters
# and additional power draw due to accessories
#====================================================================

         pgroup.P_rated[iseg] = P_plant      # power  output required from powerplant in that segment
         pgroup.Q_rated[iseg] = Q_plant      # torque output required from powerplant in that segment
         pgroup.E_rated[iseg] = Energy       # energy reqd from powerplant in that segment

#====================================================================
# if this powerplant is an engine, find effective SFC, total fuel mass
# and total flow rate 
#====================================================================

         if(pgroup.type == 'turboshaft' or pgroup.type == 'piston_engine'):
            # print(iseg, time, P_plant/1e3, flightMode, theta, delta)
            massFuel, fuelFlowRate  = pgroup.engine_wrapper(iseg, time, P_plant, flightMode, theta, delta)
            Fuel                    = Fuel + massFuel*fuel_tech
            TotalFlowRate           = TotalFlowRate + fuelFlowRate

      powerplant.seg_fuel[iseg]     = Fuel

# remember maximum flow rate for fuel system sizing
      powerplant.flow_rate          = max(powerplant.flow_rate, TotalFlowRate)

#====================================================================
# remember total fuel consumed in this segment for mass estimation in 
# next segment
#====================================================================

      # segment.sfc             = 1000.0/self.engine.sp_energy      # kg/kW-hr; not burned but used for calculations
      segment.fuel_mass       = Fuel               # fuel burned in the segment, kg

#====================================================================
# update mission parameters
#====================================================================

#      mission.totalenergyreqd += segment.energy
      
#====================================================================
# update mass segment
#====================================================================

      if (segID > 1 and i_priority >= 0):
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

      mission.segment[iseg].mass = massSegmentPrev - Fuel + addPayload
      # print(segID,massSegmentPrev,mission.segment[iseg].mass)
#====================================================================
# if weight changes over mission profile, remember this fact
# it will decide how many times we iterate for convergence 
#====================================================================

      if abs(addPayload - Fuel)>0.001:
         mission.same_weight     = False 

      return None

#====================================================================
# calculate engine weight
#====================================================================

      return mwts