#====================================================================
# Python program to perform parametric perturbations
#====================================================================

#====================================================================
# define locations and import modules
#====================================================================
import sys
sys.path.append('../../src/Python/Stage_0')
sys.path.append('../../src/Python/Stage_1')
sys.path.append('../../src/Python/Stage_1/afdd')
sys.path.append('../../src/Python/Stage_1/engines')
sys.path.append('../../src/Python/Stage_3')

import os, numpy, math, pickle, matplotlib
import matplotlib.pyplot as plt
from hydraInterface    import hydraInterface

#====================================================================
# set fonts for plotting
#====================================================================

font={'weight' : 'bold',
      'size'   : 18}
matplotlib.rc('font',**font)

#====================================================================
# quantify sensitivity of design to various uncertainties in modeling
#====================================================================

from sys import argv

if len(argv) != 2:
   print " -- Error: perturbation.py requires a design ID"
   print " -- Usage: 'python performance.py <ID> '"

n        = int(argv[1])
print " -- Chosen design ID",n
DIR      = 'output/logs/'; fname    = DIR+'log'+str(n)+'.txt'

#====================================================================
# get design ID, determine file name to read data from
#====================================================================

print 'loading HYDRA class from pickle output..',
with open(fname) as f:
   design = pickle.load(f)
print 'success!'
design.kick_off()
baseline          = design.get_essentials()
nvars             = len(baseline)
design.bemt_mode  = 1          # fix design, calculate performance
dvars             = baseline.keys()
print 'got baseline'

#====================================================================
# loop over and perturb empirical parameters, store percent change of
# Weight
# Power
# Fuel
# Radius
# Wing span 
#====================================================================

factor   = [0.9, 0.95,1.01,1.05,1.1]        # test
npert    = len(factor)

#first perturb aerodynamics
Aerod    = design.emp_data.Aerodynamics.__dict__
Weight   = design.emp_data.Tech_factors.__dict__

params   = [Aerod,Weight]

for param in params:
   out      = {}

#loop over all aero parameter "groups" (see defaults.yaml for structure)
   for group,v in param.iteritems():
      print 'key is ',group
      val         = v.__dict__ 
      out[group]  = {}

#within each sub group, loop over empirical parameters
      for key,value in val.iteritems():
         original    = value

#initialize storage of sensitivities
         out[group][key] = numpy.zeros((npert,nvars))

#loop over perturbation magnitudes
         for fcount,f in enumerate(factor):
            val[key] = original*f
            print 'setting ',group,"['",key,"'] to ",val[key]
            print 'redoing sizing'
            new      = design.get_essentials()
            per      = []

#loop over design variables of interest
            for i in xrange(nvars):
               dkey  = dvars[i]  
               diff  = numpy.round((new[dkey] - baseline[dkey])/baseline[dkey]*100,4)
               out[group][key][fcount,i] = diff

         print 'resetting ',group,"['",key,"'] to ",original
         val[key]    = original 

#check if group had any effect whatsoeve
         print factor
         print out[group][key]
         if abs(numpy.amax(out[group][key])) <= 1e-3:      # less than 0.1%: no effect 
            del out[group][key]
            print 'deleting out[',"'",group,"']['","'",key,"']"
         print dvars
      # rinp = raw_input('OK?')

   factor   = numpy.asarray(factor)
   xvals    = factor*100 - 100
   markers  = ['','o','^','s','>']
#=============================================================================
# Create plots
#=============================================================================
   for k,v in out.iteritems():        # loop over groups
      for key,value in v.iteritems():  # loop over parameters

         prefix = k + '_' + key  
         plt.figure(1)
         count = 0
         for idnum,dpar in enumerate(dvars):
            yval  = value[:,idnum]
            if abs(numpy.amax(yval)) > 1.e-3 and min(numpy.absolute(yval)) > 0.1:
               plt.plot(xvals,yval,marker=markers[count],label=dpar,linewidth=2.2)
               count = count + 1 
   
         if count > 0:
            print prefix
            plt.xlabel('% change in ' + prefix.replace('_',' ').lower())
            plt.ylabel('% change in design')
            plt.legend()
            xmin  = min(numpy.amin(xvals),numpy.amin(yval))
            xmax  = max(numpy.amax(xvals),numpy.amax(yval))
            plt.grid(linestyle='--',color='0.85')
            plt.xlim([xmin,xmax])
            plt.ylim([xmin,xmax])
            plt.tight_layout(pad=0.2)
            plt.savefig(prefix+'.png')
            plt.close()      
#=============================================================================
# repeat process for weights
#=============================================================================

