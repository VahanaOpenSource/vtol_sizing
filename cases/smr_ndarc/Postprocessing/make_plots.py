#=========================================================================
# python script to read UH60 wind-tunnel test data from Ames
#=========================================================================

import numpy, copy
import sys, os
from pylab import plt, gca
import matplotlib
import pylab

home           = os.path.expanduser('~')

#=========================================================================
# add directories to path
#=========================================================================

sys.path.insert(0,home+'/Dropbox/py/packages/plotting/')

#=========================================================================
# Import user-defined modules/routines
#=========================================================================

from save_png  import save_png
from mod_ticks import mod_ticks
from plot_fn   import modify_label
#=========================================================================
# Intialize data "super dictionary"
#=========================================================================

all_data        = {}

#=========================================================================
# data set 
#=========================================================================

all_data['185'] = numpy.loadtxt('dash185.dat',delimiter=',',skiprows=1)
all_data['230'] = numpy.loadtxt('dash230.dat',delimiter=',',skiprows=1)

with open('dash185.dat','r') as f:
    header      = f.next().split(',')

ncols           = len(header)-1

xlbl            = 'Time on station (hours)'

#=========================================================================
# create plots
#=========================================================================

plt.rc('text', usetex=True)

font = {'family' : 'Times',
      'weight' : 'bold',
      'size'   : 17}
matplotlib.rc('font', **font)

#=========================================================================
# teeter angle
#=========================================================================

for i in xrange(1,len(header)-1):
    plt.figure(i)
    for k,v in all_data.iteritems():    
        plt.plot(v[:,0],v[:,i],linewidth=2.4,marker='o')

    plt.xlabel(modify_label(xlbl))
    print modify_label(header[i].strip())
    plt.ylabel(modify_label(header[i]))
 #   a = gca()
 #   mod_ticks(gca())
    plt.grid(True,linestyle='--')
    pylab.tight_layout(pad=1.0)
    filename        = header[i].split('(')[0].strip()
    print filename
    save_png(filename)
    plt.close()

