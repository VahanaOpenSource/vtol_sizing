#====================================================================
# Python program to perform parametric perturbations
#====================================================================

from __future__ import print_function
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

from matplotlib.backends.backend_pdf import PdfPages

#====================================================================
# set fonts for plotting
#====================================================================

from matplotlib import rc

rc('text',usetex=True)
font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 18}
rc('font', **font)

#====================================================================
# quantify sensitivity of design to various uncertainties in modeling
#====================================================================

from sys import argv

if len(argv) != 2:
   print(" -- Error: perturbation.py requires a design ID")
   print(" -- Usage: 'python performance.py <ID> '")

n        = int(argv[1])
print(" -- Chosen design ID",n)
DIR      = 'output/logs/'; fname    = DIR+'log'+str(n)+'.txt'

#====================================================================
# get design ID, determine file name to read data from
#====================================================================

print('loading HYDRA class from pickle output..',)
with open(fname,'rb') as f:
   design = pickle.load(f)
print('success!')
design.kick_off()
baseline          = design.get_essentials()
nvars             = len(baseline)
design.bemt_mode  = 1          # fix design, calculate performance
dvars             = list(baseline.keys())
# print('got baseline',dvars)

#====================================================================
# loop over and perturb empirical parameters, store percent change of
# some big-picture variables
#====================================================================

factor   = [0.9, 0.95,1.01,1.05,1.1]        # test
npert    = len(factor)

#first perturb aerodynamics
Aerod    = design.emp_data.Aerodynamics.__dict__
Weight   = design.emp_data.Tech_factors.__dict__
Battery  = design.emp_data.Battery.__dict__
Costs    = design.all_dict['purchase']
params   = [Aerod,Weight,Battery,Costs]

fname    = 'sensitivities_design_' + str(n) + '.pdf'

#====================================================================
# Open PDF to write images to
#====================================================================

with PdfPages(fname) as pdf:

#====================================================================
# loop over perturbation parameters
#====================================================================

 for param in params:
   out      = {}

#====================================================================
#loop over "groups" within parameters
#====================================================================

   for group,v in param.items():

#====================================================================
# if it is a class object, continue
#====================================================================

      if isinstance(v,type(Weight['Weight_scaling'])):
         val         = v.__dict__ 
      elif isinstance(v,dict):
         val         = v 
      else:
         print('ignoring object',group)
         continue
      out[group]  = {}

#      print(val.keys())
#within each sub group, loop over design parameter assumed for sizing
      for key,value in val.items():

         original    = value

#initialize storage of sensitivities
         out[group][key] = numpy.zeros((npert,nvars))

#loop over perturbation magnitudes
         for fcount,f in enumerate(factor):
            val[key] = original*f
#            print('setting ',group,"['",key,"'] to ",val[key])
#            print('redoing sizing')
            new      = design.get_essentials()
            per      = []

#loop over design variables of interest
            for i in range(nvars):
               dkey  = dvars[i]  
               # print(dkey,baseline[dkey])
               diff  = numpy.round((new[dkey] - baseline[dkey])/baseline[dkey]*100,4)
               out[group][key][fcount,i] = diff

#         print('resetting ',group,"['",key,"'] to ",original)
         val[key]    = original 

#check if group had any effect whatsoeve
         if abs(numpy.amax(out[group][key])) <= 1e-3:      # less than 0.1%: no effect 
            del out[group][key]
            print('deleting out[',"'",group,"']['","'",key,"']")
#         print(dvars)
      # rinp = raw_input('OK?')

   factor   = numpy.asarray(factor)
   xvals    = factor*100 - 100
   markers  = ['','o','s','','','']
   styles   = ['-','--','',':',':','-']
#=============================================================================
# Create plots
#=============================================================================
   for k,v in out.items():        # loop over groups
      for key,value in v.items():  # loop over parameters

         prefix = key + '_(' + k + ')'
         plt.figure(1)
         count = 0
         for idnum,dpar in enumerate(dvars):
            yval  = value[:,idnum]
            if abs(numpy.amax(yval)) > 1.e-3 and min(numpy.absolute(yval)) > 0.05:
               plt.plot(xvals,yval,marker=markers[count],label=dpar,linewidth=2.2,linestyle=styles[count])
               count = count + 1 
   
         if count > 0:
            print(prefix)
            plt.xlabel('\% change in ' + prefix.replace('_',' ').lower())
            plt.ylabel('\% change in design')
            title_str   = 'Effect of ' + k +' ' + key 
            plt.title(title_str.replace('_',' '))
            plt.legend()
            xmin  = min(numpy.amin(xvals),numpy.amin(yval))
            xmax  = max(numpy.amax(xvals),numpy.amax(yval))
            plt.grid(linestyle='--',color='0.85')
            plt.xlim([xmin,xmax])
            plt.ylim([xmin,xmax])
            plt.tight_layout(pad=0.2)
            pdf.savefig()
            plt.close()      
#=============================================================================
# repeat process for weights
#=============================================================================
 
