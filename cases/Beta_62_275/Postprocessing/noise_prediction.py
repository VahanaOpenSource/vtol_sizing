#====================================================================
# Driver program
#====================================================================

import sys,os

sys.path.append('../../src/Python/Postprocessing/')
sys.path.append('../../src/Python/Stage_0/')
sys.path.append('../../src/Python/Stage_1/')
sys.path.append('../../src/Python/Stage_1/afdd/')
sys.path.append('../../src/Python/Stage_1/engines/')
sys.path.append('../../src/Python/Stage_3/')
sys.path.append(os.getcwd())

from noise_generator import noise 

#====================================================================
# get design ID
#====================================================================

if len(sys.argv) < 2:
   print(" -- Error: noise_prediction.py requires a design ID")
   print(" -- Usage: 'python noise_prediction.py <ID> '")
   quit()

else:
   n           = sys.argv[1]
   filename    = 'output/logs/log'+n+'.txt'
   footprint   = noise(filename)

#====================================================================
# draw noise footprint
#====================================================================
