#====================================================================
# Python program to perform parametric perturbations
#====================================================================
import matplotlib 
import os, numpy, sys
import pickle
#====================================================================
# define locations and import modules
#====================================================================

sys.path.append('../../src/Python/Stage_0')
sys.path.append('../../src/Python/Stage_1')
sys.path.append('../../src/Python/Stage_1/afdd')
sys.path.append('../../src/Python/Stage_1/engines')
sys.path.append('../../src/Python/Stage_3')

import bemt

#====================================================================
# set fonts for plotting
#====================================================================

import matplotlib.pyplot as plt
font={'weight' : 'normal',
      'size'   : 16}
matplotlib.rc('font',**font)

#======================================================================
# populate inputs to bemt
#======================================================================

interface               = bemt.bemt_interface
inputs                  = interface.inputs
inputs.nb               = 2                      # number of blades
inputs.radius           = 0.3e0                  # radius, m
inputs.solidity         = 0.12e0                 # blade solidity
inputs.omegah           = 1500.0*numpy.pi/30.0   # hover rotor speed, rad/s
inputs.use_tables       = 1                      # 1 for True, 0 for False
inputs.cl_alpha         = 5.73e0                 # 1 for True, 0 for False
inputs.cdo              = 0.028e0                # 1 for True, 0 for False
inputs.sweep_rpm        = 0                      # 1 = test all cruise RPM
inputs.rpm_ratio        = 1.0                    # cruise to hover RPM ratio
inputs.mode_operation   = 1                      # point evaluation of "best" design, fixed RPM

#======================================================================
# set design 
#======================================================================

design                  = bemt.bemt_interface.best_design
design.valid            = 1
design.taper            = 3.0
design.thx              = 18.0
design.thtip            = 24.0
design.x                = 0.6 
design.omega            = inputs.omegah 

#======================================================================
# set flight condition details
#======================================================================

inputs.n_flight         = 1
flt                     = interface.flt
print(bemt.bemt_interface.inputs)

#======================================================================
# 5 min hover segment, MSL
#======================================================================

flt[1].time             =  300.0e0          # segment time, seconds
flt[1].vc               =    0.0e0          # climb velocity, m/s
flt[1].thrust           =   58.0e0          # thrust, N
flt[1].rho              =    1.225e0        # density, kg/cu.m
flt[1].altitude         =    0.00e0         # altitude in meters

#======================================================================
# 30 min cruise segment, MSL
#======================================================================

flt[0].time             = 1800.0e0          # 30 min segment
flt[0].vc               =   18.0e0          # climb velocity, m/s
flt[0].thrust           =   11.0e0          # thrust, N
flt[0].rho              =    1.225e0        # density, kg/cu.m
flt[0].altitude         =    0.00e0         # altitude in meters

#======================================================================
# call bemt analysis: calculate performance
#======================================================================

bemt.setup_bemt()

#======================================================================
# screen output for results
#======================================================================

for i in range(inputs.n_flight):
   print ('================================')
   print ('Flight condition # %3d' % (i+1))
   print ('================================')
   print ('Root collective  is %12.3f degrees' % flt[i].coll)
   print ('Rotor efficiency is %12.3f' % (flt[i].efficiency))
   print ('Power required   is %12.3f watts' % flt[i].power)
   print (' ')

#======================================================================
# print results to file
#======================================================================

br 		= interface.best_rotor
rotor 	= {'nb'			: br.nb,
		   'radius'		: br.radius, 
		   'solidity'	: br.solidity,
		   'omega'		: br.omega, 
		   'th0'		: br.th0,
		   'power_watts': br.power, 
		   'rcout' 		: br.rcout,
		   'nondimlspan': br.r, 
		   'chord_m' 	: br.chord, 
		   'twist_deg' 	: br.twist, 
		   'alpha_deg' 	: br.alpha,
		   'inflow' 	: br.inflow, 
		   'dctdr' 		: br.dctdr,
		   'dcpdr' 		: br.dcpdr,
		   'eta' 		: br.eta}

with open("details.pkl","wb") as f:
	pickle.dump(rotor,f,protocol=2)

with open("details.pkl","rb") as f:
	a 	= pickle.load(f)

	for k,v in a.items():
		print(k,v)