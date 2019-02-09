#======================================================================
# Main driver file for the design code
# Original version in fortran by Bharath Govindarajan
# new version in python by Ananth Sridharan
#======================================================================

import os, shutil, numpy
import fea
from blade_properties      import blade_properties
from test_bemt             import run_bemt_once

class _run_hydra:

#====================================================================
# Initialize data (once only)
#====================================================================

   def initialize_bemt(self):
   
      use_bemt       = self.all_dict['sizing']['use_bemt']
      self.options   = {'bemt': use_bemt}
      if use_bemt:
         run_bemt_once(False)
      return None 
         
#====================================================================
# run hydra (this routine is done within the parameter loop)
#====================================================================

   def run_hydra(self):
   
      ifea           = self.all_dict['sizing']['ifea']
      use_bemt       = self.all_dict['sizing']['use_bemt']

#====================================================================
# reset variables before new run
#====================================================================

      aircraft       = self.all_dict['aircraft'] 
      aid            = aircraft['aircraftID']
      wing           = self.wing
      prop           = self.prop
      if(ifea==1):
         self.resetMaterialProperties()

      self.p_ins        = 1.e-3

# reset segment information
      mission           = self.mission
      self.massTakeoff  = mission.mass_takeoff_guess()

# reset powerplant properties
      self.powerplant.reset() 
      
# interpret inputs and populate class structures: done once per case
      self.interpret_inputs()

# perform iterations till convergence
      this_blade     = {}
      self.iterate_design(this_blade, False) 

      mission.reset_segments(self.massTakeoff)

#====================================================================
# for BEMT calculations, calculate power for a given thrust
# works only for quad bi plane right now.. also only use FM, eta
#====================================================================

      if use_bemt and (aid==5 or aid==7) and self.valid:
         #print 'calling bemt model'
         self.bemt_model()
         #print 'ran BEMT model: now what? redo sizing? ',valid2

#======================================================================
# rerun sizing with new FM, eta (based on "best" rotor)
#======================================================================

         if self.valid:
            self.iterate_design(this_blade, False) 

#======================================================================
#get natural frequencies here after sizing
#======================================================================

      if self.valid:
         if self.all_dict['sizing']['ifea'] ==1:
            fea.fea_types.iodata.get_modes = 1        # enable mode calculation
            temp                           = int(1)
            fea.optimizestructure(temp)
            fea.fea_types.iodata.get_modes = 0        # disable mode calculation

#======================================================================
# calculate additional details
#======================================================================

         Wmax     = self.constraints.max_gtow                  # in kgs
         Rmax     = self.constraints.max_rotor_radius/0.3048   # max radius, ft
         weight   = self.massTakeoff                     # in kgs

         if (weight > Wmax):
            print ('something is not right in error checking.. run_hydra.py')
            print (weight, Wmax, Radius, Rmax)
            sys.exit(1)

#=======================================================================
# get data directly from memory to write summary file
#=======================================================================

      self.postprocessor()

#======================================================================
# End of operations
#======================================================================
   
      return None 

#====================================================================
# function to convert input values to intermediate quantities used 
# for rotorcraft sizing; values interpreted are rotor, wing, body
# and propeller information
#====================================================================

   def interpret_inputs(self):

#====================================================================
# Assign pointer shortcuts
#====================================================================

      wing              = self.wing
      mission           = self.mission
      aircraft          = self.all_dict['aircraft']
      empirical         = self.emp_data
      prop              = self.prop 
      rotor             = self.rotor
      aid               = aircraft['aircraftID']
      nseg              = mission.nseg

#====================================================================
# conversion factors 
#====================================================================
   
      lb2kg             = self.constants.lb2kg
      grav              = self.constants.grav
      f2m               = self.constants.f2m 
         
#====================================================================
# propeller properties
#====================================================================

      try:
         Propellers        = empirical.Aerodynamics.Propellers 
         prop.eta          = Propellers.eta 
      except:
         prop.eta          = 0.82         # assign default value

#====================================================================
# rotor performance: empirical models
#====================================================================

      emp_rotors        = empirical.Aerodynamics.Rotors

#====================================================================
# for rotors mounted on fuselage: set lift fraction
#====================================================================

      if(self.fuselage.nrotors >0):
         self.fuselage.lift_frac    = aircraft['fuselage']['liftfraction']

         for iseg in range(nseg):
            segment     = mission.segment[iseg]
            flightmode  = segment.flightmode
            if(flightmode == 'hover'):
               self.fuselage.rotor_aero_eta[iseg]  = emp_rotors.FM
            else:
               self.fuselage.rotor_aero_eta[iseg]  = prop.eta
         
#====================================================================
# Wing operating condition and geometry
#====================================================================

#empirical parameters for wing
      if(self.wing.ngroups >0):
         emp_wings         = empirical.Aerodynamics.Wings 
   
#from sizing inputs   
         ngroups           = wing.ngroups
         for i in range(ngroups):
            k2             = 'group'+str(i)
            Wing           = aircraft['wing'][k2]
            w              = wing.groups[i]

            w.aspectratio  = Wing['aspectratio']

            w.oswald       = emp_wings.oswald

            w.cd0          = emp_wings.cd0
            w.nwings       = Wing['nwing']
            w.cl           = Wing['cl']

#====================================================================
# Loop over mission segments and set rotor aerodynamic efficiency defaults
#====================================================================

            for iseg in range(nseg):
               segment     = mission.segment[iseg]
               flightmode  = segment.flightmode
               if(flightmode == 'hover'):
                  w.rotor_aero_eta[iseg]  = emp_rotors.FM
               else:
                  w.rotor_aero_eta[iseg]  = prop.eta

#====================================================================
# see if this wing group has a specified lift fraction
# if not, get it from the other wing group
# in any aircraft, there can only be 2 wing groups (or less)
#====================================================================
         
            try:
               w.lift_frac = Wing['liftfraction']
            except:
               for j in range(ngroups):
                  if(j != i):
                     k3                       = 'group' + str(j)
                     wing.groups[j].lift_frac = aircraft['wing'][k3]['liftfraction']
                     w.lift_frac              = 1.0 - wing.groups[j].lift_frac 

#====================================================================
# Set rotor properties based on 
# 1. prescribed radius or disk loading
# 2. prescribed tip speed or blade loading
#====================================================================

      ngroups           = rotor.ngroups
      for i in range(ngroups):
         k2             = 'group'+str(i)
         Rotor          = aircraft['rotor'][k2]       # from classified combination list
         r              = rotor.groups[i]             # actual class instance used for calculations
         if r.set_radius:
            r.radius      = Rotor['radius']
         if r.set_DL:
            r.diskloading = Rotor['DL']*lb2kg*grav/(f2m*f2m)

#====================================================================
# do nothing for span-driven rotor 
#====================================================================

#====================================================================
# Tip speed or blade loading set
#====================================================================

         if r.set_Vtip:
            r.tipspeed    = Rotor['Vtip'] 

         if r.set_BL:
            r.ctsigma     = Rotor['ctsigma']

         if r.set_sigma:
            r.solidity    = Rotor['solidity']

#====================================================================
# blade aspect ratio, rotor solidity and # blades
#====================================================================

#         r.aspectratio    = Rotor['blade_aspect_ratio']
         r.nblade         = Rotor['Nb']
         r.nu_beta        = Rotor['flap_freq']
         r.RPM_ratio      = Rotor['cruise_rpm_ratio']

#====================================================================
# empirical properties: may change later for perturbation studies,
# hence assignment is performed here..
#====================================================================

         r.ipf          = emp_rotors.induced_power_factor 
         r.cd0          = emp_rotors.cd0
         r.hvr_dwld     = emp_rotors.hover_dwld_factor
         r.kint         = emp_rotors.kint 
         r.fm           = emp_rotors.FM
         try:
            r.thrust_share = emp_rotors.hover_thrust
         except:
            r.thrust_share = 'wing_based'
            
#====================================================================
# body properties: set for "aircraft" class
#====================================================================

      Body              = empirical.Aerodynamics.Body 
      aircraft['fdrag'] = Body.flat_plate_factor

#====================================================================
# For tilt-rotor: intialize weight at tip of fixed wing to 20% of 
# take-off weight + 2 kg ; because why not
#====================================================================

      if aid == 2:
         self.wingtip_mass = 0.2*self.massTakeoff + 2.0e0

      return None 

#======================================================================
# END OF FILE
#======================================================================