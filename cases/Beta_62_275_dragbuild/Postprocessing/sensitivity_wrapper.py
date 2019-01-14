#====================================================================
# Main function
#====================================================================

#====================================================================
# import modules used
#====================================================================

from __future__ import print_function
import sys, os  
sys.path.append(os.getcwd())
sys.path.append('../../src/Python/Postprocessing/')
sys.path.append('../../src/Python/Stage_0/')
sys.path.append('../../src/Python/Stage_1/')
sys.path.append('../../src/Python/Stage_3/')
sys.path.append('../../src/Python/Stage_1/afdd/')
sys.path.append('../../src/Python/Stage_1/engines/')

from sensitivities import dv_sensitivities
import pickle

#====================================================================
# error traps
#====================================================================

if len(sys.argv) != 2:
   print(" -- Error: performance.py requires a design ID")
   print(" -- Usage: 'python performance.py <ID> '")

#====================================================================
# get design ID and build file name, write status message
#====================================================================

n        = int(sys.argv[1])

DIR      = 'output/logs/'
fname    = DIR + 'log'+str(n)+'.txt'

with open(fname,'rb') as f:
  design    = pickle.load(f)

print(" -- Chosen design ID",n)

#====================================================================
# Initialize class, perform perturbation on wing/rotor design vars
#====================================================================

sensit 	= dv_sensitivities(design, n)
sensit.make_plots()
