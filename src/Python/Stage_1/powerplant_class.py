#====================================================================
# Powerplant class that "holds" all powerplants in dictionaries
# Note: need not be present but highly unlikely you leave without it
#====================================================================
import numpy
from turboshaft    import turboshaft
from piston_engine import piston_engine
from battery_pack  import battery_pack
from cost_class    import dict_accumulation
from conversions   import *

#====================================================================
# constants for turboshaft engine powerplants
#====================================================================

# engine accessories
k_air_engine   = 0.06
f_lub          = 1.0 
f_airind       = 0.3
f_pylon        = 0.00     

# fuel system 
f_cw          = 1.3131  # 1.3131 for crashworthiness
f_bt          = 1.8     # for ballistic survivability
nint          = 2       # number of internal tanks
fuel_scaling  = 1.42    # allow extra fuel capacity for (a) unused (b) mission flexibility
density       = 6.5     # fuel density in lb/gallons; value 6.7 is for jet fuel
#k0_plumb      = 120    # constant weight; default
#k1_plumb      = 3      # proportional weight; default
f_pump        = 1.3

#====================================================================
# begin class definition
#====================================================================

class powerplant:

   def __init__(self, all_dict, nseg):
      self.ngroups         = 0
      self.nengines        = 0
      self.groups          = {}
      self.flow_rate       = 0.0
      self.mass_fuel       = 0.0
      self.mass_battery    = 0.0
      self.seg_fuel        = numpy.zeros(nseg)
      ngrp                 = 0

      data              = all_dict['sizing']['Powerplant']

# loop over powerplant groups      
      for key in sorted(data):

# set powerplant-specific inputs through 'extras'
         if(data[key]['type'] == 'turboshaft' or data[key]['type'] == 'piston_engine'):
            extras          = all_dict['empirical']['Engines']
         elif(data[key]['type'] == 'battery'):
            extras          = all_dict['empirical']['Battery']
         else:
            quit('WHAAAAAT: powerplant_class.py unknown option')

# initialize powerplant group class
         self.groups[ngrp]  = powerplant_group(data[key], key, nseg, extras)
         ngrp               = ngrp + 1

      self.ngroups      = ngrp
      self.nmotors      = 0       # total number of motors
      return None

#====================================================================
# function to calculate total fuel used up over all segments
#====================================================================

   def fuel_rollup(self):

      self.mass_fuel       = numpy.sum(self.seg_fuel)
      self.energy_stored   = 0.0                      # rated capacity
      for groupid,group in self.groups.items():
         if(group.type == 'battery'):
            group.energy_rollup()
            self.energy_stored = self.energy_stored + group.energy_stored
      return None

#====================================================================
# this function sums the energy used up during operations from all batteries
# and calculates total battery rated capacity, as well as volume
#====================================================================

   def summarize_batteries(self):

      self.rated_energy    = 0.0      
      self.ops_energy      = 0.0                      # used in regular operations (no contingency segments)
      self.battery_vol     = 0.0
      self.mass_battery    = 0.0
      for groupid,group in self.groups.items():
         if(group.type == 'battery'):
            self.rated_energy = self.rated_energy + group.methods.Eins           # kWh
            self.ops_energy   = self.ops_energy   + group.methods.E_operations   # kWh
            self.battery_vol  = self.battery_vol  + group.methods.Vbatt          # cu.m 
            self.mass_battery = self.mass_battery + group.methods.m_batt         # kg
            
#====================================================================
# function to reset installed power for all powerplants
#====================================================================

   def reset(self):

      for i in range(self.ngroups):
         self.groups[i].reset()
      
      self.mass_fuel    = 0.0 
      self.flow_rate    = 0.0

#====================================================================
# function to calculate fuel system weights
#====================================================================

   def fuel_system(self, GTOW, fac=1.0):
      """
      this function is a modified version of the afdd model that calculates
      weights of the fuel tank, pumps and plumbing required to store and 
      transport fuel from the tank to the combustion engine

      Inputs: 
      1. GTOW: take-off weight in lbs 
      2. fac : tech_factor for fuel system 
      """

# check: if no fuel is present, return zero for mass of fuel handling
      if(self.mass_fuel == 0):
         return {'total': 0.0}

# proceed with fuel system weight estimation
      flow_rate        = self.flow_rate*kg2lb      # rated mass flow rate in lb/s
      w_fuel           = self.mass_fuel*kg2lb
      Nengine          = self.nengines
      tank_vol         = fuel_scaling * w_fuel / density           # volume, gallons

# From tank volume, find tank weight from statistical fit
      wght_tank  = ( 0.4341 * tank_vol**0.7717 * 
                     nint**0.5897 * f_cw * f_bt**1.9491 )

# find plumbing system weight from flow rate equation
      k0_plumb   = 0.011*GTOW               # constant part: % of weight (was 120 lb)
      if k0_plumb > 120:
         k0_plumb  = 120.0
      k1_plumb   = k0_plumb/40.0

      gps_engine = flow_rate/float(Nengine)*f_pump    # add margin so pumps dont fail at max capacity
      wght_plumb = k0_plumb + k1_plumb *                               \
                  (0.01*nint + 0.06*Nengine) * (gps_engine**0.866 )

# apply tech factor, calculate total 
      wght_plumb = wght_plumb*fac 
      wght_tank  = wght_tank *fac 

      total      = wght_plumb + wght_tank

# create dictionary with outputs and return it
      fuel_sys  = {'plumbing': wght_plumb*lb2kg, 
                       'tank':  wght_tank*lb2kg, 
                      'total':      total*lb2kg}

# if > 0.5 fuel weight, saturate
      if total > 0.5*w_fuel:
         fuel_sys    = {'total': 0.5*w_fuel*fac*lb2kg} 

      return fuel_sys      

#====================================================================
# function to perform weight rollup of powerplant
#====================================================================
   
   def weight_rollup(self, tech_factors, mission):

      """
      This function performs a weight roll-up for all powerplants in the system
      by sizing each powerplant group.

      Inputs:
      1. tech_factors: object containing weight multipliers (technology factors in NDARC parlance)
                       for various vehicle components
      2. mission     : class containing details of the mission (segment durations, mostly) for battery sizing

      Output:
      1. pp_wts      : dictionary containing breakdown of weights and total weight of all powerplant-related
                       components in the aircraft. Includes batteries, engines+accessories and anti-icing
      """

      pp_wts            = {}
      ngroups           = self.ngroups
      for i in range(ngroups):
         group           = self.groups[i]
         if(group.type == 'turboshaft'):
            tech_fac      = tech_factors.powerplant 
         elif(group.type == 'battery'):
            tech_fac      = tech_factors.battery 
         else:
            quit('unknown powerplant group type')

# calculate engine weight
         engine          = group.weight(tech_fac, mission)
         key             = 'group'+str(i) + 'engine'
         pp_wts[key]     = engine['total']

#add air intake anti-icing and accessories for fuel-burning engines
         if(group.type in ['turboshaft','piston']):
            anti_icing    = group.icing_weight(tech_factors.anti_icing)
            key2          = 'group'+str(i) + 'intake_heater'
            pp_wts[key2]  = anti_icing

            accessories   = group.engine_accessories()
            for k,v in accessories.items():
              key3          = 'group'+str(i) + k
              pp_wts[key3]  = v

#accumulate total masses for this powerplant group and return dictionary
      pp_wts['total'] = dict_accumulation(pp_wts)

# call summarization function for batteries 
      self.summarize_batteries()
      
      return pp_wts

#====================================================================
# Powerplant group details
# Supported powerplant types are 
# 1. turboshaft engines
# 2. piston engines
# 3. batteries 
#====================================================================

class powerplant_group:

   """
   this class contains details of a powerplant group
   can be a fuel burning engine
   """

   def __init__(self, data, key, nseg, additional_data):

      """
      Inputs
      1. data           : dictionary containing some details 
      2. key            : powerplant group name 
      3. nseg           : number of mission segments 
      4. additional_data: powerplant-specific information
      """

      self.key                = key     # what its called
      self.xmsn_group_ids     = []      # which transmission groups draw power from here
      self.mass               = {}      # mass breakdown dictionary
      self.P_rated            = numpy.zeros(nseg)
      self.powerTOP           = numpy.zeros(nseg)
      self.Q_rated            = numpy.zeros(nseg)
      self.E_rated            = numpy.zeros(nseg)
      self.P_ins              = 0.0                   # installed sea level standard condition
      self.npowerplant        = data['num']
      self.type               = data['type']          
      self.redundancy         = data['redundancy']    # how many can fail and unit still work

# for turboshaft engine
      if(self.type == 'turboshaft'):
         self.methods         = turboshaft(additional_data)

# for piston engine
      elif(self.type == 'piston_engine'):
         self.methods         = piston_engine(additional_data)

# for battery
      elif(self.type == 'battery'):
         self.methods         = battery_pack(additional_data, nseg)

      else:
         quit('unknown powerplant type: must be turboshaft, piston_engine or battery')

#====================================================================
# function to reset installed power for a powerplant 
#====================================================================

   def reset(self):
      self.P_ins              = 0.0

   def energy_rollup(self):
      self.energy_stored      = numpy.sum(self.E_rated)

#====================================================================
# function calculate weight of the powerplant
#====================================================================

   def weight(self, tech_factor, mission):
      """
      function to calculate weight of a powerplant
      Inputs are 
      1. tech factor : A scaling factor that makes the powerplant heavier or lighter
                       to account for manufacturing technology improvements over the 
                       data used to build the presently used model.
      2. mission     : class containing mission details 

      Outputs are 
      1. weight_dict : A dictionary with weight breakdown for the powerplant
      """

# turboshaft engine: call weights calculation function
      if(self.type == 'turboshaft'):
         self.P_ins           = numpy.amax(self.powerTOP)*0.001         # from W to kW
         powerplant_mass      = self.methods.getWeight(self.npowerplant, self.P_ins, tech_factor)

# Li-ion battery: call weights calculation function with relevant inputs
      elif(self.type == 'battery'):
         inputs               = {'tech_bat': tech_factor, 'mission': mission, 'P_rated': self.P_rated, \
                                 'E_stored': self.energy_stored, 'segment_energy': self.E_rated}
                                 
         powerplant_mass      = self.methods.getWeight(inputs)

# error trap
      else:
         print('engine type is ',self.type)
         quit('DONT KNOW HOW TO GET WEIGHT OF THIS POWERPLANT')

      self.mass['engine']  = powerplant_mass['total']

      return powerplant_mass

#====================================================================
# function to calculate weight of anti-icing group
#====================================================================

   def icing_weight(self, tech_factor):
      """
      function to calculate weight of engine anti-icing
      Inputs are 
      1. tech_factor: multiplier for rotor system weights to scale results from 
         parametric models up or down

      Output is 
      1. deicing equipment mass for engine air intake, kgs
      """

      mass_deicing            = self.mass['engine']*k_air_engine
      self.mass['anti_icing'] = mass_deicing
      return mass_deicing

#====================================================================
# wrapper function to evaluate take-off power, fuel mass, fuel flow rate
# and SFC for a segment
#====================================================================

   def engine_wrapper(self, iseg, time, P_req, flightMode, theta, delta):

      """
      This function takes the flight time, flight mode, temperature ratio and 
      pressure ratio of ambient air wrt mean sea level, to compute take-off power, 
      sfc, fuel mass and fuel flow rate 

      Note: when iseg = 1, the take-off power is reset to zero to avoid
            errorneously over-sizing the powerplant
      """
      
      P_engine_output      = P_req/self.npowerplant
      Pmax                 = max(self.P_ins, P_engine_output)
      # print('inputs to get fuel weight',P_engine_output/1e3, time,flightMode,theta,delta,Pmax)
      powerTOP,sfc,massFuel,fuelFlowRate = self.methods.getFuelWeight( \
         P_engine_output, time,flightMode,theta,delta,Pmax)

      self.powerTOP[iseg]  = powerTOP
      massFuel             = massFuel*self.npowerplant
      fuelFlowRate         = fuelFlowRate*self.npowerplant

      return massFuel, fuelFlowRate

#====================================================================
# wrapper function for calculating weights of engine accessories 
# used for turboshaft engines only!
#====================================================================

   def engine_accessories(self):
      """
      this function is an adaptation of the afdd engine accessories weights model
      """
      nengine         = self.npowerplant
      powerplant_tot  = self.mass['engine']*kg2lb/nengine
      s_nac           = 50.0/1600.0*self.P_ins*kw2hp*self.npowerplant
      wght_eng        = 0.667*powerplant_tot     # take 2/3 of engine weight

# engine lubrication: included in engine weight models (wet weight, not dry weight)
#      wght_lub        = 2.088 * f_lub * (wght_eng)**0.5919 * nengine**0.7858

#====================================================================
# weight of engine support structure, engine cowling, pylon 
# support structure and air induction components
#====================================================================

      wght_supt       = ( 0.0412    * (1-f_airind) * 
                        wght_eng**1.1433 * nengine**1.3762 )
   
      wght_cowl       = 0.2315*s_nac**1.3476  
   
      wght_pylon      = 0.0          # done separately!
   
      wght_airind     = 0.0412 * f_airind * wght_eng**1.1433 * nengine**1.3762

# pack breakdown into a dictionary and return
      wght_supt       = wght_supt  
      wght_cowl       = wght_cowl  
      wght_pylon      = wght_pylon 
      wght_airind     = wght_airind
   
      total           = wght_supt + wght_cowl + wght_pylon + wght_airind #+ wght_acc

      engine_acc  = {'engine_support': wght_supt *lb2kg, 'cowling':wght_cowl  *lb2kg, 
                     'pylon_support' : wght_pylon*lb2kg, 'air_ind':wght_airind*lb2kg}#,
                     # 'total': total*lb2kg}
    
      return engine_acc
