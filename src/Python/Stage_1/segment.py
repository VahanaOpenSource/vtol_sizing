
#====================================================================
# Routine to iterate multiple designs 
# Accommodates changes in AR, nblade, tip speed etc.
#====================================================================

#====================================================================
# Access fortran routines
#====================================================================

import numpy 
import hydra
#====================================================================
# Begin function
#====================================================================

class all_segments:

# ===================================================================
# function to initialize all segments
# ===================================================================

   def __init__(self, mission_data, payload):

      self.nseg            = int(mission_data['nsegments'])
      self.segment         = {}
      self.totalenergyreqd = 0.0     
      self.max_Peng        = 0.0
      self.max_c_rating    = 0.0
      self.duration        = 0.0
      self.same_weight     = True  
      n                    = self.nseg
      lst                  = mission_data['sizing_order']
      sizing_seg_ids       = []
      non_sizing_seg_ids   = []
      self.sizing_order    = []
      for i in range(n):                       # loop over sorted vals
         size_id           = lst[i]
         if size_id > 0:                        # actual sizing 
            sizing_seg_ids.append(i)
         else:
            non_sizing_seg_ids.append(i)
      lst                  = numpy.asarray(lst)
      sizing_ord_lst       = numpy.asarray(lst[sizing_seg_ids])      # segment ids used for sizing
      sizing_ord_ids       = numpy.argsort(sizing_ord_lst)           # list indices of segments used for sizing
      nsize                = len(sizing_ord_ids)
      for idval in sizing_ord_ids:
         sizing_seg_id     = sizing_seg_ids[idval] 
         self.sizing_order.append(sizing_seg_id)
      for idval in non_sizing_seg_ids:
         self.sizing_order.append(idval)
      for i in range(n):
         self.segment[i]   = mission_segment(mission_data, i, payload)

         self.duration    += self.segment[i].time 

#====================================================================
# see which segment comes first: cruise or hover
# If cruise comes first, then do span-driven sizing of rotors
# do not honor the rotor disk loading or radius inputs in that case
#====================================================================

      self.hover_first     = False
      self.cruise_first    = False 
      self.span_driven     = False 

      for i_priority in range(n):
         iseg       = self.sizing_order[i_priority]
         flt_mode   = self.segment[iseg].flightmode 
         if(flt_mode == 'hover'):
            print('\n CHOOSING CLASSICAL ROTORCRAFT SIZING: HOVER BEFORE CRUISE \n')
            self.hover_first  = True 
            break
         else:
            print('\n CHOOSING E-VTOL STYLE SIZING: CRUISE BEFORE HOVER, SPAN-DRIVEN ROTOR RADIUS \n')
            self.cruise_first = True 
            self.span_driven  = True 
            break

#====================================================================
# fixed gtow mode
#====================================================================

      if 'fixed_GTOW' in mission_data:
         self.sizing_mode        = 2
         self.gtow_target        = mission_data['fixed_GTOW']
         self.payload_tar        = 0.0
      else:
         self.sizing_mode        = 1
         self.payload            = payload
         self.payload_tar        = payload         # target payload

#====================================================================
# function to reset segment information
# not yet activated! 
#====================================================================

   def reset_segments(self):
      n                          = self.nseg 

      for i in range(n):
         segment              = self.segment[i]
         segment.fuel_mass    = 0.0e0
         segment.mass         = 0.0e0
         segment.sfc          = 0.0e0
         segment.energy       = 0.0 

#=======================================================================
# Takeoff mass: generate initial guess
#=======================================================================

   def mass_takeoff_guess(self):

#=======================================================================
# when sizing to fixed payload, set take off mass = 6 x payload mass
#=======================================================================

      if self.sizing_mode == 1:
         if self.payload >= 0.0:
            massTakeoff = self.payload*7.0 + 0.05
         else:
            massTakeoff = 1000.e0
      
#=======================================================================
# when sizing to fixed take off mass, set payload mass = 1/6 x GTOW
#=======================================================================

      elif (self.sizing_mode == 2):
         massTakeoff                = self.gtow_target
         self.payload               = self.gtow_target/6.0e0 + 0.05
         # print('setting payload weight based on take off weight guess',self.payload)
      return massTakeoff

#=======================================================================
# This class holds information for each mission segment
#=======================================================================

class mission_segment:
   def __init__(self, mission_data, j, payload):

# ===================================================================
# segment information
# ===================================================================

      self.flightmode            =     mission_data['flight_mode'][j]
      self.time                  =     mission_data['time_seg'][j]
      self.startalt              =     mission_data['start_altitude'][j]
      self.endalt                =     mission_data['end_altitude'][j]
      self.deltatempisa          =     mission_data['delta_temp_isa'][j]
      self.rateofclimb           =     mission_data['rate_of_climb'][j]
      self.cruisespeed           =     mission_data['cruise_speed'][j]
      self.addpayload            =     mission_data['add_payload'][j]
      try:
         self.distance           =     mission_data['distance'][j]
      except:
         self.distance           =     0.0
      try:
         self.add_f              =     mission_data['external_store_f'][j]
      except:
#         print ('external stores not present for segment # ',j+1)
         self.add_f              =     0.0

#====================================================================
# set default efficiencies for rotor in axial flight
#====================================================================

      # if self.flightmode == 1:            # hover: rotor FM 
      #    self.aero_eta           = 0.75e0

      # if self.flightmode == 3:            # cruise: prop eta
      #    self.aero_eta           = 0.85e0
         
#====================================================================
# calculate percent of payload dropped per segment 
#====================================================================

      if payload > 1.e-3:
         self.frac_payload_drop  =-self.addpayload/payload      # desingularize
      else:
         self.frac_payload_drop  = 0.e0

#====================================================================
# also update segment temperature and pressure ratios (one time calc)
#====================================================================

      # current_vals            = hydra.currentvalues
      alt                     = self.startalt
      temp                    = self.deltatempisa
      localatm                = hydra.isa.calculatedensityaltitude( temp, alt)
      self.rho                = localatm.rho          # kg/cu.m
      self.theta              = localatm.theta        # temperature ratio  
      self.delta              = localatm.delta        # pressure ratio
      self.density_alt        = localatm.altdens      # density altitude
      self.ainf               = localatm.speedsound   # speed of sound
      self.temperature        = localatm.temp         # temperature, deg C

#====================================================================
# these numbers are calculated
#====================================================================

      self.energy             = 0.0 
      self.fuel_mass          = 0.0e0
      self.mass               = 0.0e0
      self.P_req              = 0.0e0           # power draw, total, watts
      self.p_eng              = 0.0e0           # power draw at engine, watts
      self.c_rating           = 0.0e0           # C-rating for each segment for the PACK
      self.cell_temp          = 0.0e0           # cell temperature in deg C
      self.sfc                = 0.0e0
      self.cell_Crating       = 0.0e0           # C-rating for the cell in each segment 
      self.type               = mission_data['segment_type'][j]

#======================================================================
# reserve segments do not feature in the operating costs but still used
# in sizing of vehicle, battery/fuel 
# profit is calculated with regular mission profile
#======================================================================

      if(self.type == 'reserve'):
         self.op_cost_factor  = 0.0e0
      else:
         self.op_cost_factor  = 1.0e0

#======================================================================
# reset segment information
#======================================================================