#
# HYDRA MASTER CLASS 
#
#====================================================================
# standard Python modules
#====================================================================
#import matplotlib.pyplot as plt

import sys,os,copy,string,types,itertools
import numpy as np
import sys, os 
#from set_user_data     import *
from extract_plot_data import write_to_file, write_header
from footprint import footprint_calc
# does not write pyc
#sys.dont_write_bytecode=True
try:
   from   mpi4py import MPI
   comm     = MPI.COMM_WORLD
   rank     = comm.Get_rank()
   nprocs   = comm.Get_size()
except:
   rank     = 0 
   nprocs   = 1
   pass 

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
from component_classes     import prop_sizing, constants
from fuselage_class        import fuselage
from rotor_class           import rotors
from motor_class           import motors 
from engine_class          import all_engines
from powerplant_class      import powerplant
from wing_class            import wings
from dict2obj              import obj
from transmission_class    import transmission 

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

   def __init__(self, all_dict, MPI_found):

      self.all_dict        = all_dict
      ad                   = all_dict
      self.MPI_found       = MPI_found
      
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
# initialize classes for rotor -> must be present!
#====================================================================

      self.rotor        = rotors(ad['sizing']['Rotors'], self.mission.nseg)

#====================================================================
# initialize fuselage class
#====================================================================

      self.fuselage     = fuselage(ad['sizing']['Fuselage'],self.mission.nseg)

#====================================================================
#initialize classes for constants and propellers
#====================================================================

      self.constants    = constants()
      self.prop         = prop_sizing(ad['aircraft'])

#====================================================================
# Powerplant class: producing mechanical or electrical energy
#====================================================================
      
      self.powerplant   = powerplant(ad, self.mission.nseg)

#====================================================================
# transmission class: to send energy from powerplant to rotor
#====================================================================

      self.transmission = transmission(ad, self.mission.nseg, self.rotor.ngroups)

#====================================================================
# initialize motor classes as necessary
#====================================================================

#      if(self.transmission.motors_present):
#         self.motor     = motors(ad['empirical']['Motors'], self.mission.nseg)
#      else:
#         self.motor     = motors()

#====================================================================
# wing design parameters: initialize only when wing present
#====================================================================

      if('Wings' in all_dict['sizing']):
         self.wing      = wings(all_dict['sizing']['Wings'], self.mission.nseg)
      else:
         self.wing      = wings()

#====================================================================
# Create connections to remember rotor -> wing, rotor -> motor and 
# fuselage -> rotor maps
#====================================================================

      self.map_all(all_dict['config'])
      
#====================================================================
# find number of tilting rotors, edgewise rotors and cruise propellers
#====================================================================

      if(bool(self.wing)):
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
               print('rotor group id is ',rgid)
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
      self.run_hydra()
      array                = np.zeros(8)
      if self.valid:
         array[0]          = 1 
      else:
         array[0]          = 0
      array[1]             = self.massTakeoff            # kg
      array[2]             = self.p_ins                  # kW, installed power

      if(self.powerplant.mass_fuel > 1e-3):
         array[3]          = self.powerplant.mass_fuel   # kgs

      if(self.mass_battery > 1e-3):
         array[4]          = self.powerplant.mass_battery# kgs

      array[5]             = self.mission.payload        # kg
      array[6]             = self.costs.variable_costs   # USD/hr

      array[7]             = self.footprint

      return array

#====================================================================
# function to calculate, update and store footprint
#====================================================================

   def footprint_update(self):
      wing                 = self.wing 
      rotor                = self.rotor
      ngroups              = wing.ngroups
      footprint            = 0.0
      spans                = []
      for i in range(ngroups):
         group             = wing.groups[i]
         rid               = group.rotor_group_id
         for i in range(group.nwings):
            spans.append(group.span)

      if(ngroups == 0):
         spans.append(rotor.groups[0].radius*2)
      l_fus                = self.geometry.fuselage_length
      self.footprint       = footprint_calc(spans, l_fus)

      return None

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
      geom              = self.geometry
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
# weight margin or payload negative
#====================================================================

#      if self.massEmptyGroup['margin'] < 0 or self.msision.payload < 0:
#         self.valid        = False
#         self.errmsg       = 'margin negative, or payload negative'

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
# footprint constraint
#====================================================================

            if  self.footprint > constraint.max_rotor_radius:
               self.valid     = False 
               self.errmsg    = 'total footprint violated!: 2R = ' + str(round(self.footprint,3)) + ' m'

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
# mapping from wing group(s)/fuselage -> rotor group --> Xmsn --> powerplant
#====================================================================

   def map_all(self, config_dict):

#====================================================================
# rotor performance/sizing: need to know mapping from wings to rotors
# Presently, only one rotor group allowed per wing
#====================================================================

      Nrotors              = 0
      if(self.wing.ngroups > 0 and 'Wings' in config_dict):

         Wing_rotor_map    = config_dict['Wings']

         for wing_group_name,rotor_group_name in Wing_rotor_map.items():

            wing_group_id  = find_group_id(wing_group_name, self.wing)
            rotor_group_id = find_group_id(rotor_group_name, self.rotor)

            wgroup         = self.wing.groups[wing_group_id]

            wgroup.rotor_group_id = rotor_group_id
            self.rotor.groups[rotor_group_id].wing_group_ids.append(wing_group_id)

# find total # rotors in system: add rotors from all wings in this group
            Nrotors        = Nrotors + wgroup.nrotors

#====================================================================
# fuselage-mounted rotors: edgewise rotors (SMR/Coax)
# Presently, only one rotor group allowed per fuselage
#====================================================================

      if(self.fuselage.nrotors > 0):
         rotor_group_name              = config_dict['Fuselage']
         rotor_group_id                = find_group_id(rotor_group_name, self.rotor)

         self.fuselage.rotor_group_id  = rotor_group_id
         rgroup                        = self.rotor.groups[rotor_group_id]
         rgroup.type                   = 'lift'
         rgroup.fuse_group_id          = 0

         Nrotors                       = Nrotors + self.fuselage.nrotors

         print('assuming fuselage-mounted rotors are edgewise.. set attributes based on inputs')

#====================================================================
# remember total number of rotors
#====================================================================
      
      self.rotor.nrotors               = Nrotors 

#====================================================================
# Create rotor -> transmission mapping
# loop over all rotor groups
#====================================================================

      RXMap       = config_dict['Rotors']

      for rgname,xgname in RXMap.items():
         rgid     = find_group_id(rgname, self.rotor)
         xgid     = find_group_id(xgname, self.transmission)

# one rotor can be powered by a single transmission group
         self.rotor.groups[rgid].xmsn_group_id = xgid

# remember input power type for this rotor group
         self.rotor.groups[rgid].input_power = self.transmission.groups[xgid].type

# one transmission can transmit power to multiple rotor groups         
         self.transmission.groups[xgid].rotor_group_ids.append(rotor_group_id)

#====================================================================
# Create transmission -> powerplant maps
#====================================================================

      XPmap       = config_dict['Transmission']
      Nplant      = 0
      for xgname,xgdata in XPmap.items():

         pgname   = xgdata['Powerplant']
         xgid     = find_group_id(xgname, self.transmission)
         pgid     = find_group_id(pgname, self.powerplant)
         prat     = xgdata['PowerFraction']
# one transmission can be connected to multiple powerplants
         tgroup   = self.transmission.groups[xgid]
         pgroup   = self.powerplant.groups[pgid]

         tgroup.powerplant_group_ids.append(pgid)
         tgroup.powerplant_power_rat.append(prat)

         if(tgroup.type == 'electric' and pgroup.type in ['turboshaft','piston_engine']):
            quit('add logic for including generator after engine')
            
# one powerplant can be connected to multiple transmissions
         pgroup.xmsn_group_ids.append(xgid)

         Nplant   = Nplant + pgroup.npowerplant

# remember total number of engines
      self.powerplant.nengines = Nplant

#====================================================================
# Count how many rotors there are in a single rotor group
#====================================================================
      
      if(self.wing.ngroups > 0):
         for i in range(self.rotor.ngroups):
            wing_group_ids     = self.rotor.groups[i].wing_group_ids
      
            nrotors            = 0
            for wgid in wing_group_ids:
               group           = self.wing.groups[wgid]
               nrotors         = nrotors + group.nrotors

#find total # of rotors of this group type, in the aircraft
            self.rotor.groups[i].nrotors  = nrotors 

#====================================================================
# count how many rotors there are in the fuselage
#====================================================================

      if(self.fuselage.nrotors > 0):
         rgid                             = self.fuselage.rotor_group_id
         self.rotor.groups[rgid].nrotors  = self.fuselage.nrotors 

#====================================================================
# Count how many motors there are in a single motor group with the 
# following logic: wing group has rotor groups, which are powered by motor groups
# backtrack as follows:
#====================================================================

#====================================================================
# Loop over transmission groups
#====================================================================

      for i in range(self.transmission.ngroups):
         group    = self.transmission.groups[i]
         if(group.type == 'electric'):

#====================================================================
# Find which rotor groups are associated with this transmission
#====================================================================

            rotor_group_ids    = group.rotor_group_ids

#====================================================================
# Loop over rotor sets belonging to this transmission group
# count how many drive motors are needed
#====================================================================

            nmotors            = 0
            for rgid in rotor_group_ids: 
               nmotors         = nmotors + self.rotor.groups[rgid].nrotors

#====================================================================
# remember values in a class
#====================================================================

            group.nmotors  = nmotors

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
            print ('%20s %2d\t | \t  %s' % ('group number: ',i, group.key))
            print ('  number of rotors:      | \t ',group.nrotors)            

            print (' linked xmsn groups:     | \t ',group.xmsn_group_id)

         if(self.wing.ngroups > 0):

            print( ' linked wing  groups:    | \t ',group.wing_group_ids)
            print('\n------------------------------------------------')
            print('            Wing groups and associations     ')
            print('------------------------------------------------')
            for i in range(self.wing.ngroups):
               group       = self.wing.groups[i]
               print ('%20s %2d\t | \t  %s' % ('group number: ',i, group.key))
               print (' linked rotor groups:    | \t ',self.rotor.groups[group.rotor_group_id].key)
               print (' number of wings    :    | \t ',group.nwings)

         if(self.transmission.ngroups > 0):
            print('\n------------------------------------------------')
            print('        Transmission groups and associations     ')
            print('------------------------------------------------')
            for i in range(self.transmission.ngroups):
               group       = self.transmission.groups[i]
               print ('%20s %2d\t | \t %s' % ('group number',i, group.key))
               keys        = []
               for gid in group.rotor_group_ids:
                  keys.append(self.rotor.groups[gid].key)
               print (' linked rotor groups:    | \t',keys)
               if(group.nmotors > 0):
                  print (' number of motors   :    | \t ',group.nmotors)

         if(self.powerplant.ngroups > 0):
            print('\n------------------------------------------------')
            print('        Powerplant groups and associations     ')
            print('------------------------------------------------')
            for i in range(self.powerplant.ngroups):
               group       = self.powerplant.groups[i]
               print ('%20s %2d\t | \t %s' % ('group number',i, group.key))
               keys        = [] 
               for gid in group.xmsn_group_ids:
                  keys.append(self.transmission.groups[gid].key)
               print (' linked xmsn groups:     | \t',keys)

            print('\n')
      if(self.MPI_found):
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
   print('target key is ', target_key)
   quit('ERROR: COULD NOT MAP IT! CONFIGURATION NOT VALID')