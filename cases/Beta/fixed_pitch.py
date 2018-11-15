#====================================================================
# Python program to perform parametric perturbations
#====================================================================
import matplotlib 
matplotlib.use('Agg')

#====================================================================
# define locations and import modules
#====================================================================
import sys
sys.path.append('../../src/Python/Stage_0')
sys.path.append('../../src/Python/Stage_1')
sys.path.append('../../src/Python/Stage_1/afdd')
sys.path.append('../../src/Python/Stage_1/engines')
sys.path.append('../../src/Python/Stage_3')

import os, numpy, math, pickle
import matplotlib.pyplot as plt
from hydraInterface    import hydraInterface
import bemt

#====================================================================
# set fonts for plotting
#====================================================================

font={'weight' : 'normal',
      'size'   : 16}
matplotlib.rc('font',**font)

#====================================================================
# quantify sensitivity of design to various uncertainties in modeling
#====================================================================

from sys import argv

if len(argv) != 2:
   print(" -- Error: fixed_pitch.py requires a design ID")
   print(" -- Usage: 'python fixed_pitch.py <ID> '")

n        = int(argv[1])
print (" -- Chosen design ID",n)
DIR      = 'output/logs/'; fname    = DIR+'log'+str(n)+'.txt'

#====================================================================
# get design ID, determine file name to read data from
#====================================================================

print ('loading HYDRA class from pickle output..',)
with open(fname,'rb') as f:
   design = pickle.load(f)
print ('success!')
design.kick_off()
baseline          = design.get_essentials()
nvars             = len(baseline)
design.bemt_mode  = 1          # fix design, calculate performance
dvars             = baseline.keys()
print ('restored baseline with variable-pitch, variable RPM configuration')
print ('generating fixed-pitch, variable RPM performance curves')
#print baseline

#====================================================================
# run BEMT analysis in fixed-pitch mode  
#====================================================================

design.bemt_mode = 2 
design.bemt_model()

interface               = bemt.bemt_interface
N                       = design.mission.nseg
flt                     = interface.flt
inputs                  = interface.inputs

for i in range(N):
   n0                   = flt[i].nvalid_rpms[0]
   T_fixedth            = flt[i].t_rpmsweep[:n0,0]
   P_fixedth            = flt[i].p_rpmsweep[:n0,0]

   n1                   = flt[i].nvalid_rpms[1]
   T_varth              = flt[i].t_rpmsweep[:n1,1]
   P_varth              = flt[i].p_rpmsweep[:n1,1]

   fname                = 'flight_condition_' + str(i+1) + '.png'

   plt.figure(i)

   l1    = 'fixed coll = ' + str(round(flt[0].coll,1))
   l2    = 'vary  coll = ' + str(round(flt[i].coll,1))

   plt.plot(T_fixedth,P_fixedth,'x',label=l1)
   plt.plot(T_varth,P_varth,'-',label=l2)
   plt.plot(flt[i].thrust,flt[i].power,'ro',markerfacecolor='None',markersize=10,markeredgewidth=2.5,label='design point')
   plt.xlabel('Thrust (N)')
   plt.ylabel('Power (watts)')
   title_str            = 'Flight condition # ' + str(i+1) + ': V = ' + str(round(flt[i].vc,1)) + ' m/s'
   plt.title(title_str)
   plt.legend()
   plt.ylim(ymin=0,ymax=flt[i].power*1.5)
   plt.xlim(xmax=flt[i].thrust*1.5)
   plt.grid(linestyle='--')
   plt.tight_layout(pad=0.25)
   plt.savefig(fname,format='png',dpi=150)
   plt.close()
#   print ('target thrust, power are %6.3f %6.3f' %(flt[i].thrust, flt[i].power/746))
#   for j in xrange(n0):
#      print ('%3d %6.2f %6.2f %6.2f' % (j,Omega[j]/inputs.omegah*100,Trange[j],Prange[j]/746))

   
