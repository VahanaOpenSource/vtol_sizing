import sys,os,os.path
import numpy as np
import yaml
import matplotlib
import matplotlib.pyplot as plt

font={'weight' : 'normal',
      'size'   : 14}
matplotlib.rc('font',**font)


# find number of files
# path joining version for other paths
DIR = '../output/logs/'
xaxis = 'disk_loading' # wing_ar or solidity or disk_loading or tip_speed or number_blades


nfiles  = len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])

print ' -- Total number of files',nfiles

# allocate
xval   = np.zeros(nfiles,dtype=float)
gtow   = np.zeros(nfiles,dtype=float)
power  = np.zeros(nfiles,dtype=float)
wfus   = np.zeros(nfiles,dtype=float)
wrotor = np.zeros(nfiles,dtype=float)
ctsig  = np.zeros(nfiles,dtype=float)

# loop through the file
kk=0
for n in range(nfiles):
   if n < 10:
      fname=DIR+'log00'+str(n)+'.yaml'
   elif n<100:
      fname=DIR+'log0'+str(n)+'.yaml'
   elif n<1000:
      fname=DIR+'log'+str(n)+'.yaml'


   # load yaml file
   with open(fname) as f:
      fyaml=yaml.load(f)
   
   # proceed only if valid
   if fyaml['valid'][0] == 1:

      if xaxis=='wing_ar':
         xval[kk] = fyaml['wing']['aspect_ratio'][0]
      else:
         xval[kk]   = fyaml['rotor'][xaxis][0]

      gtow[kk]   = fyaml['vehicle']['gtow'][0]
      power[kk]  = fyaml['vehicle']['power_installed'][0]
      wfus[kk]   = fyaml['vehicle']['empty_weight']['fuselage'][0]
      wrotor[kk] = fyaml['vehicle']['empty_weight']['rotor']['total'][0]
      ctsig[kk]  = fyaml['rotor']['ct_sigma'][0]
      kk+=1



# sort array
xval0   = np.sort(xval[0:kk])
gtow0   = [x for (y,x) in sorted(zip(xval[0:kk],gtow[0:kk]))]
power0  = [x for (y,x) in sorted(zip(xval[0:kk],power[0:kk]))]
wfus0   = [x for (y,x) in sorted(zip(xval[0:kk],wfus[0:kk]))]
wrotor0 = [x for (y,x) in sorted(zip(xval[0:kk],wrotor[0:kk]))]
ctsig0  = [x for (y,x) in sorted(zip(xval[0:kk],ctsig[0:kk]))]


# plot figure
if xaxis=='solidity':
   xlabel = 'Rotor solidity'
elif xaxis=='disk_loading':
   xlabel = 'Disk loading, lb/ft$^2$'
elif xaxis=='tip_speed':
   xlabel = 'Tip speed, ft/s'
elif xaxis=='number_blades':
   xlabel = 'Number of blades'
elif xaxis=='wing_ar':
   xlabel = 'Wing aspect ratio'
else:
   print '\n -- what is the xlabel?\n\n'
   sys.exit(1)


plt.figure()
plt.plot(xval0,gtow0,'o-',color='C3',label='GTOW')
plt.xlabel(xlabel)
plt.ylabel('GTOW, lb')
plt.legend()
plt.grid(linestyle='--')

plt.figure()
plt.plot(xval0,wfus0,'o-',color='C3',label='Wt Fuselage')
plt.xlabel(xlabel)
plt.ylabel('Airframe weight, lb')
plt.legend()
plt.grid(linestyle='--')

plt.figure()
plt.plot(xval0,wrotor0,'o-',color='C3',label='Wt Rotor')
plt.xlabel(xlabel)
plt.ylabel('Rotor weight, lb')
plt.legend()
plt.grid(linestyle='--')

plt.figure()
plt.plot(xval0,power0,'o-',color='C3',label='Power')
plt.xlabel(xlabel)
plt.ylabel('Installed power, hp')
plt.legend()
plt.grid(linestyle='--')

plt.figure()
plt.plot(xval0,ctsig0,'o-',color='C3',label='Hover $C_T/\sigma$')
plt.xlabel(xlabel)
plt.ylabel('Hover $C_T/\sigma$')
plt.legend()
plt.grid(linestyle='--')




plt.show()

