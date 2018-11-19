#
# HYDRA CLASS INTERFACE
#
#####################################################################

#====================================================================
# standard Python modules
#====================================================================
#import matplotlib.pyplot as plt

import sys,os,copy,string,types,itertools
import numpy as np
import sys, os 
#from set_user_data     import *
from extract_plot_data import write_to_file, write_header

# does not write pyc
#sys.dont_write_bytecode=True
from   mpi4py import MPI
comm     = MPI.COMM_WORLD
rank     = comm.Get_rank()
nprocs   = comm.Get_size()

#====================================================================
# main class definition
#====================================================================

import loops,set_inputs
import FEA_class
import update_airframe_weight
import march_one_step
import iterate_design
import run_hydra
import engines
import empty_weight
import logdata
import cruise
import extract_summary_data
import hover
import idle
import flat_plate_buildup
import pbmodel, bemt_performance
from segment               import all_segments 
from component_classes     import prop_sizing, constants, wings 
from rotor_class           import rotors, motors
from dict2obj              import obj

#====================================================================
# HYDRA MASTER CLASS
#====================================================================

class hydraInterface(loops                  . _loops,
                     set_inputs             . _set_inputs,
                     FEA_class              . _airframeconfig,
                     update_airframe_weight . _update_airframe_weight,
                     march_one_step         . _march_one_step,
                     iterate_design         . _iterate_design,
                     engines                . _engines,
                     empty_weight           . _empty_weight,
                     logdata                . _logdata,
                     run_hydra              . _run_hydra,
                     cruise                 . _cruisemodel,
                     extract_summary_data   . _postprocessor, 
                     hover                  . _hovermodel,
                     idle                   . _idlemodel,
                     bemt_performance       . _bemt_model,
                     flat_plate_buildup     . _find_f):
   """
   hydraDriver class to provide Python API calls

   @author: Bharath Govindarajan, Ananth Sridharan
   """

#====================================================================
# class constructor: only input is a bunch of dictionaries
#====================================================================

   def __init__(self,all_dict):

      self.all_dict = all_dict
      ad            = all_dict

#====================================================================
# initialize some numbers that are remembered for output
#====================================================================

      self.mass_battery    = 0.0e0
      self.f_plate         = 0.0e0 
      self.maxfuelflowrate = 0.0e0 
      self.maxtorque       = 0.0e0 
      self.Volume          = 0.0e0 

#====================================================================
# initialize tip weight for tilt rotor wing
#====================================================================

      self.wingtip_mass    = 0.0e0 

#====================================================================
# interpret empirical parameter data and mission information
#====================================================================
      
      self.emp_data     = obj(ad['empirical'])
      self.constraints  = obj(ad['sizing']['Constraints'])

#====================================================================
# Find payload and common equipment required
#====================================================================

      payload           = self.find_payload(ad['aircraft'])
      seats_etc         = self.find_common_equip(ad['aircraft'])
      ad['aircraft']['common_equipment']  = seats_etc

      self.mission      = all_segments(ad['mission'], payload)

#====================================================================
# Powerplant type
#====================================================================

      self.etype        = self.all_dict['aircraft']['engineType']

#initialize classes for rotor and wing sizing information
      self.constants    = constants()
      # print 'initializing wing',all_dict['sizing']['Wings']
      self.prop         = prop_sizing(ad['aircraft'])
      self.motor        = motors(all_dict['empirical']['Motors'], self.mission.nseg)
      self.rotor        = rotors(all_dict['sizing']['Rotors'], self.mission.nseg)
      self.wing         = wings(all_dict['sizing']['Wings'], self.mission.nseg)

#====================================================================
# Create connections to remember rotor -> wing and rotor -> motor maps
#====================================================================

      self.map_all(all_dict['config'])
      
#====================================================================
# find number of tilting rotors, edgewise rotors and cruise propellers
#====================================================================

      for i in range(self.wing.ngroups):

         rgid           = self.wing.groups[i].rotor_group_id 
         group          = self.rotor.groups[rgid]

         nrotors        = self.wing.groups[i].nrotors
         if group.type == 'tilting':
            self.rotor.ntilt     = self.rotor.ntilt + nrotors
         elif group.type == 'lift':
            self.rotor.nlift     = self.rotor.nlift + nrotors 
         elif group.type == 'cruise':
            self.rotor.ncruise   = self.rotor.ncruise + nrotors
         else:
            print('CRITICAL ERROR: unknown rotor type - ',group.type)
            quit('valid types are tilting, lift or cruise')

#====================================================================
# perform some pre-processing
#====================================================================
      
      self.bemt_mode    = 0
      self.first_init()      
      self.kick_off()

#====================================================================
# remember output file name
#====================================================================

      if(rank == 0):
         self.fname     = 'summary.dat'
      else:
         self.fname     = 'summary_p' + str(rank) + '.dat'

#====================================================================
# function to find mass of seats, HUD, heater elements
# roll up is 14 kg / passenger + whatever user gives as fixed input
#====================================================================

   def find_common_equip(self, data):
      if('pax_count' in data):
         npax           = max(0,data['pax_count'])
      else:
         npax           = 0

      common_pax        = data['common_per_pax']
      mass_common       = common_pax*npax + data['common_equipment']  

      return mass_common 

#====================================================================
# function to find payload from a given passenger count
# other option is that its explicitly specified in kilograms, in which
# case the passenger count -> payload map is used and added!
#====================================================================

   def find_payload(self, data):

      payload     = 0.0
      if('mass_payload' in data):
         payload  = payload + data['mass_payload']            # mass, kg 

      if('pax_count' in data):
         npax        = data['pax_count']
      else:
         npax        = 0

#======================
# 1st passenger 150 kg
#======================

      if npax > 0:
         payload  = payload + 150.0

#======================
# 2nd, 3rd: 125 kg each
#======================

         n           = min(npax-1,2)         # number of pax at 125 kg each
         n           = max(n,0)
         payload     = payload + 125*float(n) 

#=========================
# 4th onwards: 120 kg each
#=========================

         n           = max(npax-3,0) 
         payload     = payload + 120*float(n)

      return payload
#====================================================================
# interpret entries in input dictionaries to sizing variables
#====================================================================

   def kick_off(self):

#====================================================================
# initialize bemt and fea memory
#====================================================================

      self.initialize_bemt()

      if (self.all_dict['sizing']['ifea'] == 1):
         self.initializeFEA()

      return None
      
#====================================================================
# setup parameter sweeps
#====================================================================

   def setupLoops(self):

      lists, inp_h      = self.setup_loops() 

      flat              = [[(k,v) for v in vs] for k, vs in sorted(lists.items())]     # alphabetical sorting
      all_combinations  = [dict(items) for items in itertools.product(*flat)]

#====================================================================
# find average # of cases per processor
#====================================================================
   
      ncases            = len(all_combinations)
      ncase_avg         = int(ncases/nprocs)
#====================================================================
# Workload subdivision done here
#====================================================================

      for icom,com in enumerate(all_combinations):
         all_combinations[icom]['icom'] = icom 
         all_combinations[icom]['rank'] = nprocs-1   # default rank is last process

#====================================================================
# loop over each process rank and see if case ID can be assigned 
# to this process
#====================================================================

         for prank in range(nprocs):
            id_start          = ncase_avg*prank                       # starting ID
            id_end            = id_start + ncase_avg - 1
            if (icom >= id_start and icom <= id_end):
               all_combinations[icom]['rank']    = prank 

#====================================================================
#debug prints
#      comm.Barrier()
#      if(rank == 0):
#         for icom,com in enumerate(all_combinations):
#            print('combination #',icom+1,' assigned to rank # ',all_combinations[icom]['rank'])
#      quit()
#end debug prints
#====================================================================

#====================================================================
# write header for summary.dat
#====================================================================

      here              = os.getcwd()
      fname             = self.fname
      f                 = open(self.fname,'w')
      f.close()

      if(rank == 0):
         with open(fname,'w') as f:
            write_header(f, inp_h)
      return all_combinations,ncases

#====================================================================
# log data regarding a design
#====================================================================

   def writeLogData(self, icom):
      
      filename = './output/logs/log'+str(icom)+'.yaml'

      self.writelogdata(filename)

#====================================================================
# update combination dictionary based on internal constraints that 
# change the original inputs 
#====================================================================

#====================================================================
# write output to file
#====================================================================

   def write_summary(self, com):

      with open(self.fname,'a') as f:
         write_to_file(com, self.summary, f)
      
#====================================================================
# clean up memory: right now, only FEA has allocatable memory,??
# what about old sizing code?
#====================================================================

   def cleanup(self):
      if (self.all_dict['sizing']['ifea'] == 1):
         self.FEA_clean()
      return None

#====================================================================
# get weight, power, fuel, radius and wing span
#====================================================================

   def get_essentials(self):

#====================================================================
# run code and extract quantities of interest -> dictionary
#====================================================================

      self.first_init()
      self.set_and_run()
      array                = {}
      array['Weight']      = self.massTakeoff            # kg
      array['Power']       = self.p_ins                  # kW, installed power

      if(self.mass_fuel > 1e-3):
         array['Fuel']     = self.mass_fuel              # kgs

      if(self.mass_battery > 1e-3):
         array['Battery']  = self.mass_battery           # kgs
      wing                 = self.wing 
      ngroups              = wing.ngroups
      for i in range(ngroups):
         group             = wing.groups[i]
         rid               = group.rotor_group_id

         array['Radius']   = rotor.groups[rid].radius           # meters
         array['Span0']    = group.span    # meters
      array['Payload']     = self.mission.payload        # kg
      array['Cost']        = self.costs.variable_costs   # USD/hr
      array['valid']       = self.valid                  # valid design or not?
      return array

#====================================================================
# check if sizing iterations have converged
#====================================================================
      
   def check_converged(self, error, itercount, icontinue):

#====================================================================
# Initialize data
#====================================================================
   
      constraint        = self.constraints
      wing              = self.wing 
      rotor             = self.rotor 
      aid               = self.all_dict['aircraft']['aircraftID']
      geom              = self.emp_data.Geometry
      self.valid        = True
      self.errmsg       = 'valid design'
      iterCountMax      = 100
      errorLimit        = 1.e-3           
      Rmax              = constraint.max_rotor_radius
      Wmax              = constraint.max_gtow
      BLmax             = constraint.max_ct_sigma      # blade loading limits
         
      if (error < errorLimit): 
         converged = True
      else:
         converged = False 

#======================================================================
# Blade loading within bounds
#======================================================================

      ngr               = rotor.ngroups
      ctsig_c           = -1.0
      for i in range(ngr):
         group          = rotor.groups[i]
         ctsigma        = group.ctsigma
         ctsig_c        = max(ctsig_c, ctsigma-BLmax)

#====================================================================
# rotor radius > constraint
#====================================================================

         if (group.radius >= Rmax):
            converged         = True            
            self.errmsg       = '   -- warning: max rotor radius exceeded constraint \n' + \
                               str(group.radius) + ' > ' + str(Rmax)
            self.valid        = False
            converged         = True 

      if ctsig_c < 0:
         self.valid  = True  
      else:
         self.valid  = False
         self.errmsg = 'Hover blade loading too high: CT/sigma = ' + str(round(rotor.groups[0].ctsigma,3))

#====================================================================
# too many iterations
#====================================================================

      if(itercount >= iterCountMax):
         self.valid     = False 
         self.errmsg    = '   -- warning: not converged after ' + str(iterCountMax) + ' iterations'
         converged      = True 
      
#====================================================================
# too heavy
#====================================================================

#      print(self.massTakeoff, Wmax)
      if (self.massTakeoff >= Wmax):
#         print('anytime now!!')
         converged            = True 
         self.valid           = False
         self.errmsg          = ' --- warning: weight limit exceeded \n' + \
                                str(self.massTakeoff) + ' > ' + str(Wmax)
         
#====================================================================
# FEA, blade weight model did not converge/find a good solution
#====================================================================

      if icontinue != 1:
         self.valid        = False
         converged         = True 
         self.errmsg       = 'something went wrong with high fidelity FE models'

#====================================================================
# For tilt-rotor/tilt-wing: check wing span and radius ratios
#====================================================================

      if aid == 2:

#====================================================================
# check that wing span is not too large
#====================================================================

         bfus           = geom.fuselage_width
         size           = 1.0 + geom.clearance*0.5     # clearance for rotor on each side
         for i in range(wing.ngroups):
            # print('this is group # ',i)
            group       = wing.groups[i]
            nr          = group.nrotors/group.nwings
            multiplier  = size*float(nr) - 1.0

            rgid        = group.rotor_group_id
            radius      = rotor.groups[rgid].radius
            # print('HI',group.span-bfus,rotor.radius*6,multiplier)
            if  group.span > constraint.max_rotor_radius:
               self.valid     = False 
               self.errmsg    = 'wing span = ' + str(round(group.span,3)) + \
                                '> max size allowed: ' + str(constraint.max_rotor_radius)

#====================================================================
# check that rotor fits on wing
# (span - fuselage)/2 = available length along which rotors can be placed
#====================================================================
         
            Rmax     = (group.span - bfus)*0.5/multiplier
            # print('multiplier is ',multiplier)
            if radius/Rmax > 1.001:
               self.valid     = False 
               self.errmsg    =  ' radius ' + str(round(radius,4)) +  ' m > allowed to fit on wing: span = ' + str(round(group.span,2)) + ' m; ' +\
                                 'available half span = ' + str(round((group.span-bfus)/2.0,4)) + ' m'  + \
                                 ' Span covered by rotors is ' + str(round(multiplier*radius,4)) +' m'
            # print('multiplier is ' + str(multiplier) + '\n' + \
                             # 'wing span is ' + str(round(group.span,2)))
#            print(rotor.radius, Rmax)
#            print(self.errmsg)
#====================================================================
# check that rotor fits on smaller wing
#====================================================================

      return converged


#====================================================================
# mapping from wing group(s) -> rotor group 
#and rotor group(s) -> motor group
#====================================================================

   def map_all(self, config_dict):

#====================================================================
# rotor performance/sizing: need to know mapping from wings to rotors
#====================================================================

      Wing_rotor_map    = config_dict['Wings']
      
      for wing_group_name,rotor_group_name in Wing_rotor_map.items():
         wing_group_id  = find_group_id(wing_group_name, self.wing)
         rotor_group_id = find_group_id(rotor_group_name, self.rotor)

         self.wing.groups[wing_group_id].rotor_group_id = rotor_group_id
         self.rotor.groups[rotor_group_id].wing_group_ids.append(wing_group_id)

#====================================================================
# motor sizing: need to know mapping from rotors to motors
#====================================================================

      Rotor_motor_map   = config_dict['Rotors']
      for rotor_group_name,motor_group_name in Rotor_motor_map.items():
         rotor_group_id  = find_group_id(rotor_group_name, self.rotor)
         motor_group_id  = find_group_id(motor_group_name, self.motor)

         self.rotor.groups[rotor_group_id].motor_group_id = motor_group_id
         self.motor.groups[motor_group_id].rotor_group_ids.append(rotor_group_id)

#====================================================================
# Count how many rotors there are in a single rotor group
#====================================================================
      
      for i in range(self.rotor.ngroups):
         wing_group_ids     = self.rotor.groups[i].wing_group_ids
      
         nrotors            = 0
         for wgid in wing_group_ids:
            group           = self.wing.groups[wgid]
            nrotors         = nrotors + group.nrotors

#find total # of rotors of this group type, in the aircraft
         self.rotor.groups[i].nrotors  = nrotors 

#====================================================================
# Count how many motors are there in a single motor group with the 
# following logic:
# wing group has rotor groups, which are powered by motor groups
# backtrack as follows:
#====================================================================

#====================================================================
# Loop over motor groups
#====================================================================

      Ntotal      = 0

      for i in range(self.motor.ngroups):
         mgroup   = self.motor.groups[i]

         nmotors  = 0 

#====================================================================
# Find which rotor groups are associated with these motors
#====================================================================

         rotor_group_ids    = mgroup.rotor_group_ids

#====================================================================
# Loop over rotor groups and find wings groups where they are mounted
#====================================================================

         for rgid in rotor_group_ids: 
            nmotors         = nmotors + self.rotor.groups[rgid].nrotors

#====================================================================
# remember values in a class
#====================================================================

         mgroup.nmotors = nmotors
         Ntotal         = Ntotal + nmotors
         
      self.rotor.nrotors   = Ntotal
      self.motor.nmotors   = Ntotal 

#====================================================================
# Print values to screen for a sanity check
#====================================================================

      if(rank == 0):
         print('hello there: here is a summary of rotors, wings and motors')
         print('\n------------------------------------------------')
         print('            Rotor groups and associations     ')
         print('------------------------------------------------')
         for i in range(self.rotor.ngroups):
            group       = self.rotor.groups[i]
            print ('%20s %2d\t | \t %s' % ('group number: ',i, group.key))
            print (' number of rotors   :    | \t ',group.nrotors)            
            print (' linked motor groups:    | \t ',group.motor_group_id)
            print( ' linked wing  groups:    | \t ',group.wing_group_ids)

         print('\n------------------------------------------------')
         print('            Wing groups and associations     ')
         print('------------------------------------------------')
         for i in range(self.wing.ngroups):
            group       = self.wing.groups[i]
            print ('%20s %2d\t | \t %s' % ('group number: ',i, group.key))
            print (' linked rotor groups:    | \t ',group.rotor_group_id)
            print (' number of wings    :    | \t ',group.nwings)

         print('\n------------------------------------------------')
         print('            Motor groups and associations     ')
         print('------------------------------------------------')
         for i in range(self.motor.ngroups):
            group       = self.motor.groups[i]
            print ('%20s %2d\t | \t %s' % ('group number',i, group.key))
            print (' linked rotor groups:    | \t ',group.rotor_group_ids)
            print (' number of motors   :    | \t ',group.nmotors)
         print('\n')
      comm.Barrier()

#====================================================================
# Given a dictionary with a collection of groups (rotors, motors or 
# wings, this function finds the dictionary key of the group that 
# matches a string value contained under the value "key")
#====================================================================

def find_group_id(target_key, group_class):
   ngroups     = group_class.ngroups
   for i in range(ngroups):
      key    = group_class.groups[i].key 
      if(key == target_key):
         print('match found: group # ',i,' has key ', key)
         return i
         break
   quit('ERROR: COULD NOT MAP IT! CONFIGURATION NOT VALID')