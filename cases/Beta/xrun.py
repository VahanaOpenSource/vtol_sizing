# ===================================================================
# Import built-in modules
# ===================================================================

from __future__ import print_function
import os, sys, time, shutil, numpy
import copy
import itertools
import yaml
import matplotlib.pyplot as plt

from mpi4py import MPI 
comm      = MPI.COMM_WORLD
rank      = comm.Get_rank()
nprocs    = comm.Get_size()

# ===================================================================
# Driver file for hydra sizing
#
# Bharath Govindarajan, Ananth Sridharan
# ===================================================================

#sys.dont_write_bytecode=True

# ===================================================================
# Set path and import routines
# ===================================================================

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
if(rank == 0):
  try:
    shutil.rmtree(log_dir)
  except:
    pass
  os.makedirs(log_dir)

comm.Barrier()

# ===================================================================
# load empirical parameters and inputs
# ===================================================================

with open("defaults.yaml","r") as f:
   fdefault = yaml.load(f)

with open("input.yaml","r") as f:
   fyaml = yaml.load(f)

# ===================================================================
# Get super-dictionary
# ===================================================================

fyaml['Sizing']. update(fdefault['Sizing'])

sizing_dict      = fyaml['Sizing']
mission_dict     = fyaml['Mission']
aircraft_dict    = fyaml['Aircraft']
configuration    = fyaml['Configuration']
empirical_dict   = fdefault['Empirical']
ops_dict         = fdefault['Operations']
acq_dict         = fdefault['Acquisition']
redund           = fdefault['Redundancies']

#beta_factors     = fdefault['Beta_factors']

all_dict = { 'empirical' : empirical_dict,  \
             'mission'   : mission_dict,    \
             'aircraft'  : aircraft_dict,   \
             'sizing'    : sizing_dict,     \
             'operations': ops_dict,        \
             'purchase'  : acq_dict,        \
             'redund'    : redund,          \
             'config'    : configuration}

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

valid_count     = 0
t1              =  time.time()
for icom,com in enumerate(all_combinations):
   
# ===================================================================
# print statement is continued
# ===================================================================

# ===================================================================
#  'case' dictionaries adjusted here for next design
# ===================================================================

  if(rank == com['rank']):

# ===================================================================
# run the case and save data to file 
# ===================================================================

    x.modify_case(com)
    x.run_hydra()
    x.write_summary(com)
    print ("\n - Rank %3d, Case: %7d of %7d ..... " %(com['rank'],icom,ncases),x.errmsg)
    if x.valid:
      x.writeLogData(icom)
      valid_count     = valid_count + 1
      x.errmsg        = x.errmsg + 'valid design!'

# ===================================================================
# final wrap-up for program termination
# ===================================================================

# ===================================================================
# print message to screen with diagnostics
# ===================================================================

comm.Barrier()
t2 =  time.time()
if(rank == 0):
  print ('\n Computation time taken is %6.4f seconds' %(t2-t1))

print ('process # ',rank,': number of valid designs are ',valid_count)

#====================================================================
# loop through all output files: do with processor rank 0
#====================================================================

if(rank == 0 and nprocs > 1):
  with open(x.fname,'a') as f:
      for i in range(1,nprocs):
#          print('getting data from processor rank # ',i)
          fname   = 'summary_p' + str(i) + '.dat'
          with open(fname,'r') as origin:
            for line in origin:
                f.write(line)
          os.remove(fname)
  t3        = time.time()
  print ('File reconstruction time taken is ',t3-t2,'seconds')

# ===================================================================
# deallocate memory and terminate program
# ===================================================================

x.cleanup()

# ===================================================================
# END OF FILE
# ===================================================================