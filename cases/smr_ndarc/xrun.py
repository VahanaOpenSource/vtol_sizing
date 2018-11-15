# ===================================================================
# Driver file for hydra sizing
#
# Bharath Govindarajan, Ananth Sridharan
# ===================================================================

print ''
print '     _  _    _  _    ____    ____     __  '
print '    / )( \  ( \/ )  (    \  (  _ \   /__\ '
print '    ) __ (   )  /    ) D (   )   /  /    \\'
print '    \_)(_/  (__/    (____/  (__\_)  \_/\_/'
print ''

print ''
print ' -------------------------------------------'
print ' Check input.yaml and set directory paths'
print ' -------------------------------------------'
print ''

# ===================================================================
# Import built-in modules
# ===================================================================

import os, sys, time, shutil
import copy
import itertools
import yaml
import matplotlib.pyplot as plt

sys.dont_write_bytecode=True

# ===================================================================
# Set path and import routines
# ===================================================================

#path_loc = os.path.normpath(os.path.join(here,'../Postprocessing'))
#sys.path.insert(0,path_loc)
sys.path.append('../../src/Python/Stage_0')
sys.path.append('../../src/Python/Stage_1')
sys.path.append('../../src/Python/Stage_1/afdd')
sys.path.append('../../src/Python/Stage_1/engines')
sys.path.append('../../src/Python/Stage_3')
sys.path.append('./Postprocessing')

from hydraInterface    import hydraInterface

# ===================================================================
# clear log data 
# ===================================================================

log_dir   = './output/logs/'
try:
  shutil.rmtree(log_dir)
except:
  pass

os.makedirs(log_dir)

# ===================================================================
# parse through the inputs and create the subset data structures
# ===================================================================

# ===================================================================
# load default file
# ===================================================================

with open("defaults.yaml","r") as f:
   fdefault = yaml.load(f)

# ===================================================================
# update with user-defined sizing file
# ===================================================================

with open("input.yaml","r") as f:
   fyaml = yaml.load(f)

fyaml['Sizing']   . update(fdefault['Sizing'])
fyaml['Aircraft'] . update(fdefault['Aircraft'])

sizing_dict      = fyaml['Sizing']
mission_dict     = fyaml['Mission']
aircraft_dict    = fyaml['Aircraft']
paths_dict       = fyaml['Paths']

empirical_dict   = fdefault['Empirical']

all_dict = { 'paths'     : paths_dict,     \
             'empirical' : empirical_dict, \
             'mission'   : mission_dict,   \
             'aircraft'  : aircraft_dict,  \
             'sizing'    : sizing_dict }

# ===================================================================
# create instance of hydraInterface class
# ===================================================================

x = hydraInterface(all_dict)

# ===================================================================
# create loops based on input.yaml 
# ===================================================================

all_combinations, ncases = x.setupLoops()

# ===================================================================
# loop through all possible combinations
# ===================================================================

t1 =  time.time()
for icom,com in enumerate(all_combinations):
   
# ===================================================================
# print statement is continued
# ===================================================================

   print ("\n - Case: %7d of %7d ..... " %(icom+1,ncases)),

# ===================================================================
#  'case' dictionaries adjusted here for next design
# ===================================================================

   x.modifyCase(com)

# ===================================================================
# run the case and save data to file 
# ===================================================================

   x.setAndRun(com)

# ===================================================================
# write out log file
# ===================================================================
   
   x.writeLogData(icom, x.getEmptyWeightGroup() )

#perf=extract_performance('./output/log.yaml')
#perf.power_curve()
#perf.energy_consumption()
#perf.hover_endurance()
#perf.cruise_performance()
#
#plt.show()

# ===================================================================
# final wrap-up for program termination
# ===================================================================

# ===================================================================
# print message to screen with diagnostics
# ===================================================================

print 'number of valid designs are ',x.valid_count
t2 =  time.time()
print 'time taken is ',t2-t1

# ===================================================================
# deallocate memory and terminate program
# ===================================================================

x.cleanup()

# ===================================================================
# END OF FILE
# ===================================================================
