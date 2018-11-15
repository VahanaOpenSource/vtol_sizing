import numpy, sys
from conversions import *
#=========================================================================   
# Write file header information
#=========================================================================   

def write_header(f, inp_h):
#replaced by input header
#   string1              = '{:^15} {:^15} {:^15}'.format('DL (lb/sq.ft)','Blade Solidity','# blades')
#   string2              = '{:^15} {:^15} {:^15}'.format('wing AR','Wing_LF' ,'Tip Speed')
#   string3              = '{:^15} {:^15} {:^15}'.format('Rpm ratio','   Flap frequency','Wing Cl target')
  
   string4              = '{:^15} {:^15} {:^15}'.format('    Weight (kg)', '   Power (kW)', 'Fuel/Batt (kg)')
   string5              = '{:^15}'.format('  Empty Wt (kg)')
   string7              = '{:^15} {:^15} {:^15}'.format('    Payload (kg)', 'Op Cost ($/hr)', 'Valid design')

   string               = inp_h + string4 + string5 + string7
   string               = string + '\n'
   f.write(string)

   return None

#=========================================================================   
# Write run information from one point to file
#=========================================================================   

def write_to_file(com, data1, f):


   icom                 = com['icom']
   Wt                   = data1['Wt']
   Power                = data1['Power']
   Fuel                 = data1['Fuel']
   Empty                = data1['Empty']
   payload              = data1['payload']
   op_cost              = data1['op_cost']
   valid                = data1['valid']

   string0              = '{:15d}'.format(icom)
   string1              = '{:^15.2f} {:^15.2f} {:^15.3f}'.format(Wt, Power, Fuel)
   string2              = '{:^15.2f} {:^15.2f} {:^15.5f}'.format(Empty, payload, op_cost) # Wlift)
   string3              = '{:^15d}'.format(valid)   

   line                 = string0 + string1 + string2 + string3 + '\n'
   f.write(line)
   
   return None

#=========================================================================   
# Break line of data into two dictionaries (input and output )
#=========================================================================   

def create_subset(array):

   subset                  = {};    data = {}

#=========================================================================   
# "Input"
#=========================================================================   

   subset['icom']          = int(array[0]);      iout = 1
   data['Wt']              = float(array[iout]); iout = iout + 1     # col index 11
   data['Power']           = float(array[iout]); iout = iout + 1     # col index 12
   data['Fuel']            = float(array[iout]); iout = iout + 1     # col index 17
   data['Empty']           = float(array[iout]); iout = iout + 1     # col index 18
   data['payload']         = float(array[iout]); iout = iout + 1     # col index 20
   data['op_cost']         = float(array[iout]); iout = iout + 1     # col index 21
   data['valid']           =   int(array[iout]); iout = iout + 1     # col index 22

   return subset, data

#=========================================================================   
# write data from FEA. Use sparingly. Only if you want to analyse 
# an isolated case.
#=========================================================================   
