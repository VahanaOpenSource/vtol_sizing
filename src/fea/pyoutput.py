import matplotlib.pyplot as plt
import numpy as np

f = open('./outputs/optimize.dat','r')

line = f.next()
line = f.next()
line    = line.strip() # gets rid of the \n at the end
columns = line.split()
#
nlines=int(columns[0])
#
fos        = np.zeros(nlines)
weight     = np.zeros(nlines)
deflection = np.zeros(nlines)
itcount    = range(1,nlines+1)
#
for j in range(nlines):
   line = f.next()
   line    = line.strip() # gets rid of the \n at the end
   columns = line.split()

   fos[j]        = columns[1]
   weight[j]     = columns[2]
   deflection[j] = columns[3]

f.close()
#
#
#

fig,ax1=plt.subplots()

ax2 = ax1.twinx()
ax1.plot(itcount, fos, 'ro-',label='FOS')
ax1.axhline(y=1.5,color='black',linestyle='--')
ax2.plot(itcount, deflection, 'bo-',label='Deflection')
ax2.axhline(y=0.1,color='black',linestyle='--')

ax1.set_xlabel('Iteration count')
ax1.set_ylabel('Factor of safety')
ax2.set_ylabel('Deflection, m')


lines1,label1 = ax1.get_legend_handles_labels()
lines2,label2 = ax2.get_legend_handles_labels()

ax2.legend(lines1+lines2,label1+label2,loc=0)
plt.grid()

fig,ax1=plt.subplots()

ax1.plot(itcount, weight, 'ro-')
ax1.set_xlabel('Iteration count')
ax1.set_ylabel('Airframe weight, kg')
ax2.set_ylabel('Deflection, m')

plt.grid()
plt.show()
