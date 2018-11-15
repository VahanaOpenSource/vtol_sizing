# ===================================================================
# Post processing file for hydra
#
# Bharath Govindarajan, Ananth Sridharan
# ===================================================================


print ''
print '     _  _    _  _    ____    ____     __  '
print '    / )( \  ( \/ )  (    \  (  _ \   /__\ '
print '    ) __ (   )  /    ) D (   )   /  /    \\'
print '    \_)(_/  (__/    (____/  (__\_)  \_/\_/'
print ''

print ''
print '-------------------------------------------'
print '                 POST-PROCESS              '
print '-------------------------------------------'
print ''

# ===================================================================
# Import built-in modules
# ===================================================================

import os,sys
import copy
import itertools
import yaml
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

sys.dont_write_bytecode=True

font={'weight' : 'normal',
      'size'   : 16}
matplotlib.rc('font',**font)

# ===================================================================
# user placed constraints
# ===================================================================

ctsigma_max   = 0.15
max_dimension = 10 # ft


# ===================================================================
# find total number of files
# ===================================================================

DIR = '../output/logs/'

nfiles = len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR,name))])

print ' -- Total number of files: ',nfiles

# ===================================================================
# allocate space
# ===================================================================

index = np.zeros(nfiles,dtype=int)
gtow  = np.zeros(nfiles,dtype=float)
power = np.zeros(nfiles,dtype=float) 

# ===================================================================
# loop through the files
# ===================================================================

kk=0
for n in range(nfiles):
   fname=DIR+'log'+str(n)+'.yaml'

# ===================================================================
   # load yaml file
# ===================================================================

   with open(fname) as f:
      fyaml=yaml.load(f)
   
# ===================================================================
# proceed only if design is valid
# ===================================================================

   if fyaml['valid'][0] == 1:

      ctsigma  = fyaml['rotor']['ct_sigma'][0] 
      radius   = fyaml['rotor']['radius'][0]  # ft
      wingspan = fyaml['wing']['span'][0]     # ft

# ===================================================================
# check if contraints met
# ===================================================================

      if ( ctsigma      < ctsigma_max   and 
           4.445*radius < max_dimension and
           wingspan     < max_dimension):

         gtow[kk]   = fyaml['vehicle']['gtow'][0]
         power[kk]  = fyaml['vehicle']['power_installed'][0]
         index[kk]  = n
         kk+=1

nvalid = kk

print ' -- Total number of valid designs: ',nvalid

# ===================================================================
# sort based on [user choice]
# ===================================================================

myList     = gtow[:nvalid]

sortID     = sorted(range(len(myList)),key=lambda x:myList[x])
sortedList = np.sort(myList) 

index0   = [x for (y,x) in sorted(zip(myList[:nvalid],index[:nvalid]))]
gtow0    = [x for (y,x) in sorted(zip(myList[:nvalid],gtow[:nvalid]))]
power0   = [x for (y,x) in sorted(zip(myList[:nvalid],power[:nvalid]))]

# ===================================================================
# some bar charts for comparison 
# ===================================================================

nshow = min(10,nvalid)
print ' -- Showing only top ',nshow,'designs'

# ===================================================================
# show GTOW
# ===================================================================

fig   = plt.figure()
ax    = fig.add_subplot(111)
ypos  = np.arange(nshow)
ax.barh(ypos,gtow0[:nshow],color='r')
ax.set_xlabel('GTOW, lb')
ax.set_ylabel('Design ID Number')
ax.set_yticks(ypos)
ax.set_yticklabels(index0)
ax.xaxis.grid(True,linestyle='--')
plt.gca().invert_yaxis()

# ===================================================================
# show Installed power
# ===================================================================

fig   = plt.figure()
ax    = fig.add_subplot(111)
ypos  = np.arange(nshow)
ax.barh(ypos,power0[:nshow],color='b')
ax.set_xlabel('Installed power, hp')
ax.set_ylabel('Design ID Number')
ax.set_yticks(ypos)
ax.set_yticklabels(index0)
ax.xaxis.grid(True,linestyle='--')
plt.gca().invert_yaxis()

plt.show()

# ===================================================================
# END OF FILE
# ===================================================================
