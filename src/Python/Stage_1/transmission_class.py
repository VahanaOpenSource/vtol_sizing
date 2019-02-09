#====================================================================
# Transmission class that transfers energy from powerplants to rotors
#====================================================================
import numpy
from drive                 import drivesys_weight
from electric_transmission import electric_transmission
from DC_motor              import DC_motor 
from cost_class            import dict_accumulation
from conversions           import *
class transmission:

   def __init__(self, all_dict, nseg=0, nrotor_groups=0):
      self.ngroups      = 0
      self.groups       = {}
      ngrp              = 0
      data              = all_dict['sizing']['Transmission']

# loop over powerplant groups, initialize group class
      for key in sorted(data):
         self.groups[ngrp]  = transmission_group(data[key], key, nseg, nrotor_groups)
         ngrp               = ngrp + 1

      self.ngroups      = ngrp
      self.nmotors      = 0       # total number of motors
      return None

#====================================================================
# function to perform weight roll-up for all transmission
#====================================================================

   def weight_rollup(self, tech_factors, rotor, powerplant):

      """
      this function performs weight estimation for mechanical and electric transmissions
      and returns a dictionary with the breakdown of component weights as well as the total
      weight of all transmissions in the vehicle

      Inputs:
      1. tech_factors: object containing weight multipliers (technology factors in NDARC parlance)
                       for various vehicle components
      2. rotor       : class containing information for all rotor groups 
      3. powerplant  : class containing information for all powerplants

      Outputs:
      1. xmsn        : dictionary containing breakdown and total of transmission component weights 
      2. p_DSlimit   : driveshaft rated power (kW) - used to size tail rotor, if present
      3. omega       : main rotor speed (rad/s)    - used to size tail rotor, if present
      """

      xmsn          = {}
      ngroups       = self.ngroups
      p_DSlimit     = 0.0
      omega         = 0.0
      for i in range(ngroups):
         group      = self.groups[i]

# setup inputs for mechanical transmission
         for irg in group.rotor_group_ids:
            if(group.type == 'mechanical'):
               rg         = rotor.groups[irg]      # rotor group
               wblade     = rg.mass_blades         # blade mass in kg for all rotors
               nrotor     = rg.nrotors             # number of rotors
               len_ds     = rg.radius* 1.25        # driveshaft length, meters
               omega      = rg.tipspeed/rg.radius  # rotor speed (rad/s)
               omega      = omega*30/numpy.pi      # rotor speed (RPM) - because NDARC uses it for calculating 
                                                   # an equivalent rated torque for mechanical transmissions

# assume only one powerplant group connected to a mechanical transmission
               ipg        = group.powerplant_group_ids[0]
               p_DSlimit  = powerplant.groups[ipg].P_ins*powerplant.groups[ipg].npowerplant
               q_DSlimit  = numpy.amax(powerplant.groups[ipg].P_rated)
               q_DSlimit  = q_DSlimit/omega*0.001*group.eta*kw2hp
               omega      = omega*numpy.pi/30.0
               vparams    = {'wblade': wblade*kg2lb, 'nrotor'   : nrotor, 
                             'len_ds': len_ds * m2f, 'nprop'    : 0, 
                             'omega' : omega,        'fac'      : tech_factors.drive_system,
                          'q_DSlimit': q_DSlimit,    'P_DSlimit': p_DSlimit*kw2hp }

# electric transmission: calculate motor weights here (wire weights done with fixed wing weight estimation)
            elif(group.type == 'electric'):
               vparams       = {'fac': tech_factors.drive_system}

# error trap: unknown transmission type
            else:
               quit('unknown transmission type: must be mechanical or electric')

# transmission components: call method from class
            xmsn_grp_mass    = group.weight(vparams)

# append group id to transmission elements and remember breakdown,total in a dictionary
            for key,value in xmsn_grp_mass.items():
               k2             = 'group'+str(i)+'xmsn_'+key 
               xmsn[k2]       = value 

# add motor weight to rotor assembly mass for electric transmissions
            if(group.type == 'electric'):
               for irg in group.rotor_group_ids:
                  rotor.groups[irg].mass_assembly += xmsn_grp_mass['motors']/rotor.groups[irg].nrotors

# get total weight and return dictionary
      xmsn['total']     = dict_accumulation(xmsn)

      return xmsn, p_DSlimit, omega

#====================================================================
# Class that remembers transmission group details
#====================================================================

class transmission_group:

   def __init__(self, data, key, nseg, nrotor_groups):

      """
      this function initializes the class for a transmission group
      Inputs are 
      1. data:          details for the transmission group (dictionary)
      2.  key:          transmission group name 
      3. nseg:          number of mission segments 
      4. nrotor_groups: number of rotor groups this transmission transmits power to
      """

# initialize transmission details
      self.key                   = key     # what its called
      self.rotor_group_ids       = []      # which rotor groups are powered by this transmission
      self.powerplant_group_ids  = []      # which powerplants supply energy to this transmission
      self.powerplant_power_rat  = []      # ratio of rated power provided by a particular powerplant
      self.nmotors               = 0
      self.power_req             = numpy.zeros((nseg,nrotor_groups))
      self.torque_req            = numpy.zeros((nseg,nrotor_groups))

# motor rated power and torque
      self.motor_p_inp           = numpy.zeros(nseg)
      self.motor_p_out           = numpy.zeros(nseg)
      self.motor_q_out           = numpy.zeros(nseg)

      self.type                  = data['type']
      self.eta                   = data['eta']

# redundancy: multiplier for weight, cost groups
      try:
         self.redundancy         = data['redundancy']
      except:
         self.redundancy         = 1.0

# for electric motors, remember hover and cruise efficiencies
      if(self.type == 'electric'):
         self.motor_eta_hover    = data['Motors']['hover_efficiency']
         self.motor_eta_cruise   = data['Motors']['cruise_efficiency']

#=======================================================================
# mass estimation wrapper function for a transmission group
#=======================================================================

   def weight(self, vparams):
      """
      function to calculate mass of mechanical transmission components
      Input:
      1. vparams : dictionary, array or float64 passed to weight calculation function

      Output:
      1. A dictionary containing breakdown of transmission elements 
      """

      if(self.type == 'mechanical'):         
         return drivesys_weight(vparams)
      elif(self.type == 'electric'):
         motor_mass  = {'motors':self.motor_sizing(vparams['fac'])}
         return motor_mass
      else:
         quit('unknown transmission type')

#====================================================================
# function to size motors using motor power input 
#====================================================================
   
   def motor_sizing(self, tech_factor):
      self.line_rated_power      = numpy.amax(self.motor_p_out)   # should be motor_p_out, kW
      P_req                      = self.line_rated_power
      motor_mass                 = DC_motor(P_req)*tech_factor
      print(P_req/motor_mass)
      return motor_mass*self.nmotors
      