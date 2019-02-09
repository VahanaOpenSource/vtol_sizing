#=========================================================================   
# Extract data to write to file summary.dat with important numbers for 
# design sizing 
# Note: uses memory from a mix of fortran modules for output
#=========================================================================   

class _postprocessor:

  def postprocessor(self):

   data              = {}                    # initialize output dict

#=========================================================================   
# set pointers
#=========================================================================   

   aircraft          = self.all_dict['aircraft'] 
   rotor             = self.rotor.groups[0]
   wing              = self.wing 
   mission           = self.mission 
   nrotors           = self.rotor.nrotors 
   
#=========================================================================   
# conversion constants
#=========================================================================   

   m2ft              = 1.0/0.3048
   kg2lb             = 2.2046
   kW2hp             = 1.0/0.746

#=========================================================================   
# set actual data here  
#=========================================================================   

   data['Wt']        = self.massTakeoff                     # in kg
   data['Power']     = self.p_ins                           # in hp
   data['Radius']    = rotor.radius                         # in meters

   max_span          = 0.0
   min_span          = 5e9
   max_crd           = 0.0 
   min_crd           = 5e9
   for i in range(wing.ngroups):
      max_span          = max(wing.groups[i].span,  max_span)
      min_span          = min(wing.groups[i].span,  min_span)
      max_crd           = max(wing.groups[i].chord, max_crd)
      min_crd           = min(wing.groups[i].chord, min_crd)
   data['span']      = max_span                             # in meters    
   data['w_crd']     = max_crd                              # in meters
      
#=========================================================================   
# for non-electric systems
#=========================================================================   

   data['Fuel']      = self.powerplant.mass_fuel                          # Fuel  weight, kg

#=========================================================================   
# for electric systems: add battery weight to "fuel"
#=========================================================================   

   data['Fuel']      = data['Fuel'] + self.mass_battery

#=========================================================================   
# 
#=========================================================================   

   data['Empty']        = self.massempty                 # Empty weight, kg

#=========================================================================   
# for invalid designs, set CTsigma to 10.0
# unrealistically high to guide postprocessor away from this design
#=========================================================================   

   data['Ctsigma']   = rotor.ctsigma                  # Blade loading 
   if self.valid:
      data['valid']     = 1
   else:
      data['valid']     = 0
   data['payload']      = self.mission.payload           # payload in kg

#=========================================================================   
# acquisition, fixed and operating cost calculations: 
# $/vehicle, $/year, $/flight hour
#=========================================================================   

   if(bool(self.costs)):
      self.costs.mass_scaling_elements(self.massTakeoff, self.massEmptyGroup, self.massempty)
      self.costs.power_scaling_elements(self.powerplant.p_ins)
      self.costs.rotor_acq_cost(self.rotor.groups)
      self.costs.cost_accumulation(self.powerplant, self.massempty)
   
#=========================================================================   
# calculate operating cost of vehicle in USD/flight hour
#=========================================================================   

      data['op_cost']      = self.costs.variable_costs
   else:
      data['op_cost']      = 0.0
      
#=========================================================================   
# store data in class
#=========================================================================   
   
   self.summary         = data 

#=========================================================================   
# return information
#=========================================================================   

   return None