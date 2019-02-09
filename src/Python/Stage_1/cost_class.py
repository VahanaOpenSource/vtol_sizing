#====================================================================
# physical constants
#====================================================================

class costs:

   def __init__(self, data, acq_data, beta_fac, mission, nrotors, npax):

#====================================================================
# memory initialization for group totals and breakdown
#====================================================================

      self.acquisition              = 0.0                         # frame acquisition cost
      self.fixed_costs              = 0.0                         # fixed operating costs
      self.variable_costs           = 0.0                         # variable operating costs
      self.fix_breakdown            = {}
      self.acq_breakdown            = {}
      self.var_breakdown            = {}

#====================================================================
# input parameters
#====================================================================

      self.tmission                 = mission.cost_duration/60.0  # mission duration in hours
      self.dmission                 = mission.cost_range          # range in km 
      self.Annual                   = data.Annual                 #
      self.Hourly                   = data.Hourly
      self.Battery                  = data.Battery 
      self.Taxi                     = data.Taxi 
      self.Vertiport                = data.Vertiport
      self.Acquisition              = acq_data 
      self.beta_fac                 = beta_fac 
      self.nrotors                  = nrotors 
      self.npax                     = npax 
      return None

#====================================================================
# function to calculate acquisition cost elements that scale with 
# max rated power (given in kW)
#====================================================================

   def power_scaling_elements(self, Pins):

      cost                 = self.acq_breakdown
      scaling              = self.Acquisition['Scaling_cost']
      cost['power_dist']   = scaling['power_dist']*Pins 

#====================================================================
# function to calculate acquisition cost elements that scale with 
# take-off mass (given in kgs)
#====================================================================

   def mass_scaling_elements(self, massTakeoff, massEmpty, total_empty_mass):
        
      cost              = self.acq_breakdown
      scaling           = self.Acquisition['Scaling_cost']

#====================================================================
# final assembly, testing and BRS costs scale with weight
#====================================================================

#      print(total_empty_mass)
      cost['final_assem_line'] = scaling['final_assem_line']*total_empty_mass
      cost['BRS']              = scaling['BRS']             *massTakeoff   # Ballistic recovery system 

#====================================================================
# motor: cost scales with weight of motor+ESC
#====================================================================

      cost['motors']           = scaling['motors']*massEmpty['drive_system']['total']

#====================================================================
# fuselage: cost scales with weight of fuselage structure
#====================================================================

      cost['fuselage']         = scaling['fuselage']*massEmpty['fuselage']['total']

#====================================================================
# landing gear cost scales with weight of its structure
#====================================================================

      cost['landing_gear']     = scaling['landing_gear']*massEmpty['alighting_gear']['total'] 

#====================================================================
# landing gear cost scales with weight of its structure
#====================================================================

      wing_structure           = 0.0 
      wing_actuators           = 0.0
      wing_tilters             = 0.0
      for k,v in massEmpty['wing'].items():
         if k.endswith('structure'):
            wing_structure     = wing_structure + v 
         elif k.endswith('tilters'):
            wing_tilters       = wing_tilters   + v 
         elif k.endswith('actuators'):
            wing_actuators     = wing_actuators + v 
      cost['wing_structure']   = scaling['wing_structure']*                 \
                               ( wing_structure + massEmpty['empennage']['total'])

#====================================================================
# weight of wires: $11.5/kg of cable based on AWG-1 rating; 
# use $23.9/kg based on Alpha vehicle estimate
#====================================================================

      cost['wires']           = scaling['wires']*massEmpty['wires']['total']

#====================================================================
#actuators: flap/ailerons and tilt mechanisms (with redundancy)
# estimate from Herve's spreadsheet puts it at $10920 for 2 wing actuators
# and $20583 for one tilt actuator -> this is per wing
# estimate includes redundancy and manufacturing in bulk
# AS modification: scale cost with actuator weights too..
# reference is Alpha vehicle actuator weight (5.7 kg for wing CS, 9.6kg for tilt)
#====================================================================

      cost['wing_flaps']       = scaling['wing_flap']    *wing_actuators
      cost['tilt_actuators']   = scaling['tilt_actuator']*wing_tilters

#====================================================================
# Loop over entries in fixed cost and store in dictionary
#====================================================================

      for key,value in self.Acquisition['Fixed_cost'].items():
         cost[key]    = value 
      return None

#====================================================================
# Acquisition cost for rotors
#====================================================================

   def rotor_acq_cost(self, rotor_groups):

      cost                    = self.acq_breakdown
      scaling                 = self.Acquisition['Scaling_cost']
      blade_costs             = 0.0 
      hub_costs               = 0.0
      for gid,group in rotor_groups.items():
         RNR                  = group.nrotors* group.radius
         AbNR                 = group.chord  * group.nblade * RNR
         blade_costs          = blade_costs  + AbNR*scaling['rotor_blade']
         hub_costs            = hub_costs    +  RNR*scaling['rotor_hub']

      cost['rotor_blade']     = blade_costs
      cost['rotor_hub']       = hub_costs
      return None 

#=============================================================================
# calculate costs for ground-only travel
#=============================================================================

   def taxi_analysis(self):
      speed                = self.taxi_speed(self.dmission)          # km/hr output, km input
      self.Taxi_time       = self.dmission/speed*60.0                # time in minutes
      self.taxi_cost       = self.Taxi.Distance_rate*self.dmission   # taxi cost in USD

      self.taxi_cost       = self.taxi_cost + self.Taxi_time*self.Taxi.Time_rate # add extra due to time rate

      self.Taxi_time       = self.Taxi_time + self.Taxi.Padding_time

#====================================================================
# for UAM: last leg distance: 2 km
#====================================================================

      distance                   = self.Vertiport.Ground_distance    # km for last leg
      speed                      = self.taxi_speed(distance)         # km/hr for last leg
      time                       = distance/speed                    # hours spent for last leg      
      taxi_cost                  = distance*self.Taxi.Distance_rate + \
                                       time*self.Taxi.Time_rate*60

      self.UAM_cost              = self.variable_costs*self.tmission + \
                                   taxi_cost 

      self.UAM_time              = (time + self.tmission)*60 + \
                                   self.Vertiport.Padding_time       # total time

      self.dollar_per_min  = (self.UAM_cost - self.taxi_cost)/(self.Taxi_time - self.UAM_time)
      if(self.npax > 0):
         self.dollar_per_min  = self.dollar_per_min/self.npax 

      return None 

#=============================================================================
# function to calculate average speed of a taxi in downtown areas
# input: distance in km, 
# output: speed in km/hr
#=============================================================================

   def taxi_speed(self, distance):
      speed             = 16.6*distance/(15.0+distance)*3.6   # ground speed, kmph
      return speed 

#====================================================================
# Fixed costs: insurance, depreciation, liability, yearly inspections
#====================================================================

   def calc_fixed_costs(self):

      Annual              = self.Annual 
      fcb                 = self.fix_breakdown
      fcb['liability']    = Annual.Liability
      fcb['inspection']   = Annual.Inspection
      fcb['insurance']    = Annual.Insurance_percent   *self.acquisition*0.01
      fcb['depreciation'] = Annual.Depreciation_percent*self.acquisition*0.01

      return None 

#====================================================================
# Variable costs 
# input is Emission - total battery charging energy to complete mission, kWh
#====================================================================

   def calc_variable_costs(self, Emission, Etotal):

#====================================================================
# initialize memory
#====================================================================

      Battery                    = self.Battery
      var                        = self.var_breakdown
      Pavg                       = Emission/self.tmission 

#====================================================================
# find electricity costs per flight hour
#====================================================================

      rate                       = Battery.Electricity     # electricity cost, $/kWh
      var['electricity']         = Pavg*rate               # cost/flight hour, ($/hr) 

#====================================================================
# Frame maintenance cost/flight hour
#====================================================================

      ops                        = self.Hourly
      var['frame_overhaul']      = ops.Frame_maintenance              # maintenance, inspection and repair cost,$/flight hr

#====================================================================
# rotor and motor inspection: averaged per flight hour
#====================================================================

      var['vpf_overhaul']    = (ops.Rotor_inspection + ops.Motor_inspection)*self.nrotors  # inspection cost for rotors, $/flight hr

#====================================================================
# Battery operating cost per flight hour
#====================================================================

      Ncycles                    = Battery.Cycles                 # number of charge/discharge cycles
      rate_b                     = Battery.Cost_per_kwh           # $/kWh of energy

#====================================================================
# assume batteries go through one charge/discharge cycle per flight
#====================================================================

      flights_per_hour           = 1.0/self.tmission              # tmission is in hours
      cycles_per_hour            = flights_per_hour               # one cycle of charge/discharge per flight
      Flight_hours               = int(Ncycles/cycles_per_hour)   # flight hours that battery lasts for
      var['battery_use']         = rate_b*Etotal/(Flight_hours)   # battery operating cost, $/flight hour 

#====================================================================
# landing fees
#====================================================================

      var['landing_fees']        = self.Vertiport.Landing_fees*flights_per_hour

      return None

#=============================================================================
# add up costs from different groups
#=============================================================================

   def cost_accumulation(self, powerplant, massEmpty):

#====================================================================
# Apply acquisition cost scaling factors going from Alpha -> Beta 
#====================================================================

      cost                       = self.acq_breakdown
#      print(self.beta_fac.keys())
#      print(cost.keys())
      for key in cost.keys():
         cost[key] = cost[key]*self.beta_fac[key]

#=============================================================================
# get acquisition cost total
#=============================================================================

      self.acquisition        = dict_accumulation(self.acq_breakdown)   # total acquisition costs
#      self.acquisition        = 550*massEmpty

#=============================================================================
# get fixed costs (yearly rate) => depends on acquisition costs
#=============================================================================

      self.calc_fixed_costs()                                           # calculate fixed costs
      self.fixed_costs        = dict_accumulation(self.fix_breakdown)   # depends on total acquisition cost

#====================================================================
# also calculate per flight hour cost from fixed cost total
#====================================================================

      self.var_breakdown['fixed_annual']  = self.fixed_costs/self.Annual.Flight_hours

#=============================================================================
# get variable cost components 
#=============================================================================

      self.calc_variable_costs(powerplant.ops_energy, powerplant.rated_energy)                        # component build up
      self.variable_costs     = dict_accumulation(self.var_breakdown)   # depends on mission
      self.taxi_analysis()

#=============================================================================
# add up key value entries and return total 
#=============================================================================

def dict_accumulation(group):

      total     = 0.0
      for key in group:
         total += group[key]

      return total