import sys,os,os.path
import numpy as np
import yaml
import matplotlib
import matplotlib.pyplot as plt

font={'weight' : 'normal',
      'size'   : 22}
matplotlib.rc('font',**font)


# find number of files
# path joining version for other paths
DIR = '../output/logs/'
xaxis = 'tip_speed' # solidity or disk_loading or tip_speed or number_blades


nfiles  = len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])

print ' -- Total number of files',nfiles

# allocate
xval   = np.zeros(nfiles,dtype=float)
gtow   = np.zeros(nfiles,dtype=float)
power  = np.zeros(nfiles,dtype=float)
wfus   = np.zeros(nfiles,dtype=float)
wrotor = np.zeros(nfiles,dtype=float)
ctsig  = np.zeros(nfiles,dtype=float)
wspan  = np.zeros(nfiles,dtype=float)
wingcl = np.zeros(nfiles,dtype=float)
wingar = np.zeros(nfiles,dtype=float)

# loop through the file
kk=0
for n in range(nfiles):
   fname=DIR+'log'+str(n)+'.yaml'


   # load yaml file
   with open(fname) as f:
      fyaml=yaml.load(f)
   
   # proceed only if valid
   if fyaml['valid'][0] == 1:
      xval[kk]   = fyaml['rotor'][xaxis][0]
      gtow[kk]   = fyaml['vehicle']['gtow'][0]
      power[kk]  = fyaml['vehicle']['power_installed'][0]
      wfus[kk]   = fyaml['vehicle']['empty_weight']['fuselage'][0]
      wrotor[kk] = fyaml['vehicle']['empty_weight']['rotor']['total'][0]
      wspan[kk]  = fyaml['wing']['span'][0]
      wingcl[kk] = fyaml['wing']['cl'][0]
      wingar[kk] = fyaml['wing']['aspect_ratio'][0]
      kk+=1


nvalid=kk

unique_wingcl = np.unique(wingcl[:nvalid])
unique_wingar = np.unique(wingar[:nvalid])

nwingcl = len(unique_wingcl)
nwingar = len(unique_wingar)

#
# plot direction 1
#

clines_gtow1=np.zeros([nwingcl,nwingar],dtype=float)
clines_span1=np.zeros([nwingcl,nwingar],dtype=float)
ndatapts1=np.zeros(nwingcl,dtype=int)

for j in range(nwingcl):
   wcl=unique_wingcl[j]

   kk=0
   for k in range(nvalid):
      if wcl==wingcl[k]:
         clines_gtow1[j,kk]=gtow[k]
         clines_span1[j,kk]=wspan[k]
         kk+=1
         
   
   ndatapts1[j]=kk

#
# plot direction 2
#

clines_gtow2=np.zeros([nwingcl,nwingar],dtype=float)
clines_span2=np.zeros([nwingcl,nwingar],dtype=float)
ndatapts2=np.zeros(nwingcl,dtype=int)

for j in range(nwingar):
   war=unique_wingar[j]

   kk=0
   for k in range(nvalid):
      if war==wingar[k]:
         clines_gtow2[j,kk]=gtow[k]
         clines_span2[j,kk]=wspan[k]
         kk+=1
         
   
   ndatapts2[j]=kk

plt.figure()

for j in range(nwingcl):
   plt.plot(clines_span1[j,:ndatapts1[j]],clines_gtow1[j,:ndatapts1[j]],'o-',linewidth=4.0,label='Wing $C_L$ %2.1f' % unique_wingcl[j])


for j in range(nwingar):
   plt.plot(clines_span2[j,:ndatapts2[j]],clines_gtow2[j,:ndatapts2[j]],'o-',linewidth=4.0,label='Wing AR %2.1f' % unique_wingar[j])

#plt.plot(clines_span2[1,:ndatapts2[1]],clines_gtow2[1,:ndatapts2[1]],'o-',color='C5',linewidth=4.0,label='Wing AR %2.1f' % unique_wingar[j])
#plt.plot(clines_span1[1,:ndatapts1[1]],clines_gtow1[1,:ndatapts1[1]],'o-',color='C1',linewidth=4.0,label='Wing $C_L$ %2.1f' % unique_wingcl[j])

#plt.plot(wspan[:nvalid],gtow[:nvalid],'o',color='C3')
plt.xlabel('Wing span, ft')
plt.ylabel('GTOW, lb')
plt.grid(linestyle='--')
plt.legend()
plt.xlim([2,10])
plt.ylim([45,70])
plt.show()
