#=========================================================================   
# Python function to call HYDRA in a loop and do pre/postprocessing
#=========================================================================   

from performance_sweep import *  
#=========================================================================
# Import built-in modules
#=========================================================================

import shutil, itertools, os, copy, numpy

#=========================================================================
# Copy template inputs to program inputs directory
#=========================================================================

code_dir             = '../../../PrasadUM/'
tar_dir              = os.path.join(code_dir,'Inputs/')
src_dir              = os.path.join(code_dir,'Inputs_samples/HYDRA_v2/')

if os.path.exists(tar_dir):
   shutil.rmtree(tar_dir)
shutil.copytree(src_dir, tar_dir)

#=========================================================================   
# Choose global properties
#=========================================================================      

Nr                      = 1               # number of rotors
Nb                      = 4               # blades per rotor
DL                      = 12              # disk loading,          lb/sq. ft
GTOW                    = 24938           # gross take-off weight, lbs
Vtip                    = 787.4           # hover tip speed,       ft/s
wing_LF                 = 0.4             # wing lift fraction
sigma                   = 0.09            # rotor solidity
lift_offset             = 0.5             # target lift offset
Mtip                    = 0.6             # adv. tip Mach # in cruise
Vcruise                 = 240             # cruise speed,          knots
density_alt             = 1829            # density altitude,      meters
blade_mass              = 198.3           # blade mass,            kilograms
nu_beta                 = 1.2             # hover flap frequency,  per rev
shaft_tilt              =-0.63            # rotor tilt angle, deg (+ve fwd)

#=========================================================================      
# Compute global properties
#=========================================================================      

Target_T                = GTOW * (1.0 - wing_LF) / float(Nr)# thrust, lbs
Target_T                = Target_T - Nb*blade_mass*2.2         # subtract blade weight
Disk_area               = 1.15*GTOW /  DL / float(Nr)             # area, sq.ft
Radius                  = numpy.sqrt(Disk_area/numpy.pi)          # Radius, ft

#=========================================================================      
# Overwrite
#=========================================================================      
#print Radius 
#quit()
#Radius                  = 18.543
#Target_T                = 14260
#Radius                  = 27.7            # radius, feet
#print Radius
#quit()
#=========================================================================      
# Derived quantities: DO NOT TOUCH
#=========================================================================      

Omega_hover             = Vtip/Radius
Vinf                    = Vcruise * 1.853 * 5/18                  # in m/s
Omega                   = (Mtip * 340.44 - Vinf) / (0.3048 * Radius)  # cruise rotor speed, rad/s
Target_L                = Target_T * lift_offset * Radius         # Roll moment, ft-lb
cbar                    = sigma * numpy.pi * Radius / float(Nb)   # mean chord, ft

mu                      = Vinf/(Omega*Radius*0.3048)

#=========================================================================      
# Convert radius to meters
#=========================================================================      

Radius                  = Radius * 0.3048 

#=========================================================================      
# Store CSD wind-tunnel trim parameters in a dictionary
#=========================================================================      

CSD_params              = {'Thrust' : Target_T/2.2*9.8,                    \
                                'L' : Target_L/2.2*9.8*0.3048,             \
                           'Radius' : Radius * 0.3048        ,             \
                                'M' : 0.0,                                 \
                             'Vinf' : Vinf,                                \
                            'alpha' : shaft_tilt*numpy.pi/180.0,           \
                               'Nb' : Nb ,                                 \
                             'cbar' : cbar,                                \
                       'density_alt': density_alt,                         \
                              'mbl' : blade_mass,                          \
                          'Omega_h' : Omega_hover,                         \
                          'Omega_c' : Omega,                               \
                          'nu_beta' : nu_beta}

#=========================================================================      
# Define parameter space for blade geometry
#=========================================================================      

#x        = [0.1]#[0.4,  0.5, 0.6, 0.7];
#thx      = [-1]#[-7.5, -5,  -2.5]
#thtip    = [-10]#[-12.5,-10, -7.5, -5]
#y        = [0.5,0.6,0.7];  
#ctip     = [1.0,0.9, 0.8, 0.7, 0.6];
#croot    = [1.0,  1.1, 1.2];

#=========================================================================      
#Baseline: 
#=========================================================================      

#BILINEAR TWIST, NO TAPER
#x        = [ 0.3,0.4,0.5,0.6,0.7];
#thx      = [-5.0,-4,-3,-2,-1];
#thtip    = [-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2]#[-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0];
#y        = [0.5]
#croot    = [1.0]
#ctip     = [1.0]

#BEST: dynamic inflow
x        = [0.3];
thx      = [-3];
thtip    = [-3];
y        = [0.5]
croot    = [1.0]
ctip     = [1.0]

#=========================================================================      
# DO NOT TOUCH: dimensionalize chord and generate list combinations
#=========================================================================      

croot_list              = [cbar*c for c in croot]
ctip_list               = [cbar*c for c in ctip]

lists                   =  [x,thx,thtip,y,croot_list,ctip_list]
all_combinations        = list(itertools.product(*lists))
ncases                  = len(all_combinations)

#=========================================================================      
# Initialize output file with header
#=========================================================================      

f = open('../../Outputs/Stage_1/performance_summary.dat','w') 
write_header(f)
f.close()

#=========================================================================   
# Initialize sweep
#=========================================================================   

for icom,com in enumerate(all_combinations):

   print 'running case # ',icom+1, 'of', ncases

#=========================================================================   
# run the code and remember outputs
#=========================================================================   

   CSD_outputs     = CSD_caller(com, CSD_params)

#=========================================================================   
# Write information to file/screen/whatever
#=========================================================================   

   f = open('../../Outputs/Stage_1/performance_summary.dat','a')
   write_to_file(com, CSD_outputs, f)
   f.close()

#print 'number of valid designs is ',valid_count
