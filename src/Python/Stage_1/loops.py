#=======================================================================
#
# Setup parameter space to sweep over
#
#=======================================================================
import os,sys
import numpy as np 

try:
   from mpi4py import MPI 
   comm  = MPI.COMM_WORLD
   rank  = comm.Get_rank()
except:
   rank  = 0
   
from data_processing   import run_read
class _loops: 

# ######################################################################
#
# setup_loops
#
# ######################################################################

   def setup_loops(self): 
   
      case_sen_dict = self.all_dict['sizing']
      sizing_dict   = {}
      for key in sorted(case_sen_dict):
         sizing_dict[key.lower()] = case_sen_dict[key]

#=======================================================================
# screen print
#=======================================================================

      if(rank == 0):
         print('\n=========================|========================')
         print('   Input parameters\t | \t values            ')
         print( '=========================|========================')

         for key,value in sorted(sizing_dict.items()):
            if type(value) is dict:
               print (key)
               for k,v in sorted(value.items()):
                  if type(v) is dict:
                     print ('\t ',k)
                     for k2,v2 in sorted(v.items()):
                        print ('%20s \t | \t %s' % (k2, str(v2)))
                  else:
                     print ('%20s \t | \t %s' % (k, str(v)))
            else:
               print ('%20s \t | \t %s' % (key, str(value)))
         print( '==================================================')
      
      if(self.MPI_found):
         comm.Barrier()

#=======================================================================
# Extract entries and add to a list of lists
#=======================================================================

      all_lists   = {}

#=======================================================================
# wing design variables
#find # wing groups. Members within a group are identical   
#=======================================================================

      if 'wings' in sizing_dict:
#loop over wing groups, find design parameters
         Wings       = sizing_dict['wings']
         for igroup in range(self.wing.ngroups):
            grp            = self.wing.groups[igroup]
            key            = grp.key            # how to access this object through dict addressing
            svars          = Wings[key]
            st             = 'Wing_group' + str(igroup)
            for key,value in svars.items():

#=======================================================================
#loop over design parameters that are not rotor-related
#=======================================================================

               st1         = st.lower() + '_' + key 
               all_lists[st1] = value

#=======================================================================
# Rotor design variables
#=======================================================================

      Rotors       = sizing_dict['rotors']
      keys_to_pop  = []
      for igroup in range(self.rotor.ngroups):
         grp            = self.rotor.groups[igroup]
         key            = grp.key            # how to access this object through dict addressing
         svars          = Rotors[key]
         st             = 'Rotor_group' + str(igroup)


         for key,value in svars.items():

#=======================================================================
# for span-driven rotor sizing, remove disk loading and radius from 
#=======================================================================

            if grp.span_driven and (key == 'radius' or key == 'DL' or key == 'disk_loading'):
               print('AHA! found you! I will not include this value in sizing loops: ',key)
            else:
               st1         = st.lower() + '_' + key 
               all_lists[st1] = value
               # print(st1,value)
#=======================================================================
#loop over design parameters pertaining to rotors
# the dictionary containing these details is now called "value"
#=======================================================================

#=======================================================================
# add default values if not found, and throw error if inputs n/a
#=======================================================================

         if 'Nb' not in svars:
            print('could not find #blades: defaulting to 3')
            Nb             = [3]
            st1            = st.lower() + '_' + 'number_blades'
            all_lists[st1] = Nb  

         if 'flap_freq' not in svars:
            print('fl_freq not found for this rotor set: defaulting to flap frequency of 1.1/rev')
            fl_freq     = [1.1]
            st1            = st.lower() + '_' + 'flap_freq'
            all_lists[st1] = fl_freq

         if 'radius' not in svars and 'DL' not in svars:
            pref           = st.lower()
            print(st.lower(),' has neither radius nor disk loading specified: ')
            print('span driven option chosen to size this rotor group')

#=======================================================================
# cruise rpm ratio 
#=======================================================================

         if 'cruise_rpm_ratio' not in svars:
            print('defaulting to cruise RPM ratio = 50%')
            rpm_ratio   = [0.5]
            st1            = st.lower() + '_' + 'cruise_rpm_ratio'
            all_lists[st1] = rpm_ratio

#=======================================================================
# tip speed 
#=======================================================================

         if 'Vtip' not in svars:
            st1            = st.lower() + '_' + 'rotor_Vtip'
            Vtip           = [170]
            print('Defaulting to 170 m/s for hover tip speed')
            all_lists[st1] = Vtip

#=======================================================================
# rotor solidity and blade loading both not found
#=======================================================================
            
         if ('solidity' not in svars) and ('ctsigma' not in svars):
            print('solidity and blade loading both not found: assigning defaults')
            print('defaulting to rotor solidity of 0.1')
            solidity       = [0.1]
            st1            = st.lower() + '_' + 'solidity'
            all_lists[st1] = solidity

#=======================================================================
# rotor solidity and blade loading both found
#=======================================================================

         elif ('solidity' in svars) and ('ctsigma' in svars):
            print('cannot specify both solidity and blade loading: pick one!')
            quit('CRITICAL PROGRAM ERROR: TERMINATING')

#=======================================================================
# only solidity specified: blade loading calculated
#=======================================================================

         else:
            print('only one of solidity OR blade loading specified: using it to calculate the other')

#=======================================================================
# first column: design id #
#=======================================================================

      strings = 'Design ID #'

#=======================================================================
# assemble header string
#=======================================================================

      string1              = '{:^15}'.format(strings)
# End of operations
      return all_lists, string1

# ######################################################################
#
# Modify case with itertools generated loop
#
# ######################################################################

   def modify_case(self, com):
 
# alter the aircraft dictionary directly based
# on the particular combination of inputs
# that is incoming through 'com' --> combination of DVs

      aircraft_dict  = self.all_dict['aircraft']

#========================================================================
# Wing design variables
#========================================================================

      aircraft_dict['wing']         = {}
      aircraft_dict['rotor']        = {}

      if(bool(self.wing)):
         for i in range(self.wing.ngroups):
            aircraft_dict['wing']['group'+str(i)]          = {} 

      for key,lst in com.items():
         if key.startswith('wing_'):
            kspl                       = key.split('_')
            k0                         = kspl[0]            # wingX
            k1                         = kspl[1]            # groupY 
            replace_str                = k0+'_'+k1+'_'
            k                          = key.replace(replace_str,"")
            aircraft_dict[k0][k1][k]   = lst
#========================================================================
# For rotor, nest properties in a dictionary under wing dictionary
#========================================================================

      for i in range(self.rotor.ngroups):
         aircraft_dict['rotor']['group'+str(i)]          = {} 

      for key,lst in com.items():
         if key.startswith('rotor_'):
            kspl                       = key.split('_')
            k0                         = kspl[0]            # rotorX
            k1                         = kspl[1]            # groupY 
            replace_str                = k0+'_'+k1+'_'
            k                          = key.replace(replace_str,"")
            aircraft_dict[k0][k1][k]   = lst

#========================================================================
# calculate blade aspect ratio from solidity
#========================================================================

      for i in range(self.rotor.ngroups):
         rotor                       = aircraft_dict['rotor']['group'+str(i)]
#         rotor['blade_aspect_ratio'] = rotor['Nb']/np.pi/rotor['solidity']

#========================================================================
# set sizing option switches and remember in class group
#========================================================================

         class_group                 = self.rotor.groups[i]
         if('radius' in rotor):
            class_group.set_radius   = True 
         elif('DL' in rotor):
            class_group.set_DL       = True 
         else:
            class_group.span_driven  = True 

         if('Vtip' in rotor):
            class_group.set_Vtip     = True 
         
         if('ctsigma' in rotor):
            class_group.set_BL       = True 
         
         if('solidity' in rotor):
            class_group.set_sigma    = True 

#========================================================================
# create dictionary for data for one parametric combination
#========================================================================

      self.current_dict = {}
      self.current_dict = {'weight_hist' : [] }

      self.com          = com          # remember original combination in self
      return None
   
