#====================================================================
# quantify sensitivity of design to various uncertainties in modeling
#====================================================================

import sys,os
sys.path.append('../../src/Python/Postprocessing/')
sys.path.append('../../src/Python/Stage_0/')
sys.path.append('../../src/Python/Stage_1/')
sys.path.append('../../src/Python/Stage_1/afdd/')
sys.path.append('../../src/Python/Stage_1/engines/')
sys.path.append('../../src/Python/Stage_3/')
sys.path.append(os.getcwd())

if len(sys.argv) == 1:
   print(" Performing gradient-based optimization, multiple initial conditions")
   method	= 'gradient_based'
else:
   print(" Performing single differential evolution optimization")
   method	= 'diffl_evo'

from optimize import optimize

optimize(os.getcwd(), method)

