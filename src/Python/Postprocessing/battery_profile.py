#=================================================================
# file to visualize battery usage over mission duration
#
#=================================================================

import yaml,sys,os,numpy
import matplotlib
import matplotlib.pyplot as plt

#==========================
# use these lines for latex
#==========================
from matplotlib import rc

rc('text',usetex=True)
font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 12}

rc('font', **font)
#==========================

#=====================================
# use this line to save images to PDF
#=====================================

from matplotlib.backends.backend_pdf import PdfPages

#=====================================

class battery_profile:
   
# ===================================================================
# initialize and load yaml file
# ===================================================================

   def __init__(self,fname,n):

      print(" -- Chosen design ID",n)

      try:
         with open(fname) as f:
            fyaml=yaml.load(f)
      except:
         print('looked for file',fname)
         print('CRITICAL ERROR: file does not exist')
         quit()

# ===================================================================
# basic dictionary pointers and parameters
# ===================================================================

      self.mission      = fyaml['Mission']
      self.vehicle      = fyaml['vehicle']
      self.battery      = fyaml['Battery']
      self.id           = n 
#====================================================================
# This function creates a line plot for the battery assumed in sizing
#====================================================================

   def current_battery(self): 

#====================================================================
# define ideal battery min/max limits
#====================================================================

      CAmin       = 5.0 
      CAmax       = 100.0 
      CArange     = CAmax-CAmin;

#====================================================================
#loop over all segments (including reserve)
#====================================================================

      m           = self.mission
      b           = self.battery
      soh         = b['state_of_health']
      E_init      = soh*100 #*b['rated_capacity']
      E_final     = b['depth_discharge']#*b['rated_capacity']

      nseg        = m['nsegments']
      dE          = numpy.zeros(nseg) 
      Estatus     = numpy.zeros(nseg+1); Estatus[0] = E_init
      time        = numpy.zeros(nseg+1)

      plt.figure(1)
      plt.subplot(211)
      plt.title('\\textbf{Perfect battery assumptions}')
      plt.ylim([0,100])
      plt.xlabel('\\textbf{Mission time (min)}')
      plt.ylabel('\\textbf{Battery charge (\%)}')
      plt.subplot(212)
      plt.title('\\textbf{Vahana battery assumptions}')
      plt.ylim([0,100])
      plt.ylabel('\\textbf{Battery charge (\%)}')
      plt.xlabel('\\textbf{Mission time (min)}')
      plt.tight_layout(pad=1)

#====================================================================
# find total time of mission
#====================================================================

      #      plt.grid(linestyle=':')
      dist        = 0.0
      for i in range(nseg):

         P              = m['battery_power_draw'][i]           # kW
         t              = m['time'][i]                         # hrs
         V              = m['cruise_speed'][i]*1.853/3.6       # in m/s

         dE[i]          = P*t/60.0/b['rated_capacity']*100.0   # % of rated energy used
         time[i+1]      = t + time[i]                    # time counter since start of mission
         Estatus[i+1]   = Estatus[i] - dE[i]
         typ            = m['segment_type'][i].rstrip()
         print(i,P)
         if(typ == 'all'): 
            style       = '-bo'
         elif(typ == 'reserve'):
            style       = '--b'
            reserve_pow = P 
            reserve_V   = V
            reserve_d   = V*t*60*0.001        # in km
#            print('reserve power is ',reserve_pow)
#====================================================================
# now find extra energy assumed by ***UNNAMED DIVISION***
#====================================================================

      E_used            = numpy.sum(dE)                   # % energy used
      if(E_used < 100.0):
         delE           = 100.0 - E_used
         print('Percentage of energy actually used  = ',100 -delE, ' %')

         delE           = CArange - E_used
         str1           = str(round(delE,2))
         print('Ideal battery: assumes ' + str1 + '% more usable energy for sizing')

         deltaEnergy    = delE*b['rated_capacity']*0.01     # convert unused percent to actual kWh value
         dt             = deltaEnergy/reserve_pow*60        # additional reserve time 
         d              = dt*reserve_V*60/1000              # extra reserve distance
         # print('extra time is ',dt,reserve_V,d);quit()
      else:
         print('warning: energy used is ',round(E_used,0), '% available is ', b['rated_capacity'], ' kWh')
#====================================================================
# set x limits on plot
#====================================================================
      
      # print(time[-1]+1+dt)
      plt.subplot(211)
      print(dt)
      plt.xlim([time[0], time[-1]+1+dt])
      plt.subplot(212)
      plt.xlim([time[0], time[-1]+1+dt])

#====================================================================
# find type of segment
#====================================================================

      for i in range(nseg):

         P              = m['battery_power_draw'][i]           # kW
         t              = m['time'][i]                         # hrs

         typ            = m['segment_type'][i].rstrip()
         V              = m['cruise_speed'][i]*1.853/3.6       # in m/s

         dtemp          = V*t*60*0.001                         # in km
         dist           = dist + dtemp        

         slope       = numpy.arctan2(dE[i],t)*180.0/numpy.pi 
         x0          = time[i] + 0.3*t
         y0          = Estatus[i] - 0.3*dE[i] +2
         l2          = numpy.asarray([x0,y0])
         trans_angle = plt.gca().transData.transform_angles(numpy.array((slope,)),
                                                l2.reshape((1, 2)))[0]
         if(typ == 'all'): 
            lstyle      = '-'
            lcolor      = 'blue'
            marker      = 'o'

            if(V != 0):
               str1        = '\\textbf{range = ' + str(round(dist,1)) + ' km}'
               # print('annotating at',x0,y0)
               plt.annotate(str1,xy=(x0,y0),xytext=(x0,y0),rotation=-trans_angle,color=lcolor)
         elif(typ == 'reserve'):
            lstyle      = ':'
            lcolor      = 'green'
            marker      = ''
            str1        = '\\textbf{\\qquad' + str(round(reserve_d,0)) + ' km reserve}'
            x0          = time[i]
            y0          = Estatus[i]+2
            l2          = numpy.asarray([x0,y0])

            plt.annotate(str1,xy=(x0,y0),xytext=(x0,y0),rotation=-trans_angle,color=lcolor)
         else:
            quit('unknown type of mission segment: can be all or reserve')


         # print(Estatus[i],Estatus[i+1])
         if(i == nseg-1):
            lstr        = 'Vahana \\quad    : range = ' + str(round(dist-reserve_d,1)) + ' km'
            plt.plot((time[i],time[i+1]),(Estatus[i],Estatus[i+1]),marker=marker,linestyle=lstyle,color=lcolor,linewidth=2.4,label=lstr)
         else:
            plt.plot((time[i],time[i+1]),(Estatus[i],Estatus[i+1]),marker=marker,linestyle=lstyle,color=lcolor,linewidth=2.4)

      plt.plot((0,time[-1]),(Estatus[-1],Estatus[-1]),':b',linewidth=1.0,alpha=0.5)

#====================================================================
# now draw **UNNAMED** published profile with 100% -> 5% battery
#====================================================================
      
      Estatus[0]        = 100.0
      time[0]           = 0.0
      dist              = 0.0
      for i in range(nseg):
         V              = m['cruise_speed'][i]
         t              = m['time'][i]                         # hrs         
         typ, style     = self.segment_type(i)
         # print(typ)
         if(typ == 'reserve'):
            t           = t + dt 
            dE[i]       = dE[i] + delE
            # print('YOOO',delE*b['rated_capacity']/100/(dt/60))
         time[i+1]      = t + time[i]                    # time counter since start of mission
         Estatus[i+1]   = Estatus[i] - dE[i]
         dist           = dist + V*1.853*5.0/18.0*t*60*0.001        # in km

      for i in range(nseg):

         typ, style     = self.segment_type(i)

         if(i == nseg-1):
            lstr        = 'Ideal battery: range = ' + str(round(dist,1)) + ' km'
            plt.subplot(211)
            plt.plot((time[i],time[i+1]),(Estatus[i],Estatus[i+1]),style,linewidth=2.4,)
            plt.subplot(212)
            plt.plot((time[i],time[i+1]),(Estatus[i],Estatus[i+1]),style,linewidth=2.4,alpha=0.3)
         else:
            plt.subplot(211)
            plt.plot((time[i],time[i+1]),(Estatus[i],Estatus[i+1]),style,linewidth=2.4)
            plt.subplot(212)
            plt.plot((time[i],time[i+1]),(Estatus[i],Estatus[i+1]),style,linewidth=2.4,alpha=0.3)

         t              = time[i+1] - time[i]
         delE           = Estatus[i+1] - Estatus[i]
         if(typ == 'reserve'):
            reserve_d   = reserve_V*1.853*5.0/18.0*t*60*0.001        # in km
            str1        = '\\textbf{range = ' + str(round(dist,1)) + ' km}'
            # print(time[i],t,Estatus[i],delE)

            slope       = numpy.arctan2(dE[i],t)*180.0/numpy.pi 
            x0          = time[i] + 0.3*t
            y0          = Estatus[i] - 0.3*dE[i] +2
            l2          = numpy.asarray([x0,y0])
            trans_angle = plt.gca().transData.transform_angles(numpy.array((slope,)),
                                                   l2.reshape((1, 2)))[0]

            plt.subplot(211)
            plt.annotate(str1,xy=(time[i]+0.1*t,Estatus[i]),xytext=(time[i]+0.1*t,Estatus[i]),rotation=-trans_angle,color='r')
            plt.subplot(212)
            plt.annotate(str1,xy=(time[i]+0.1*t,Estatus[i]),xytext=(time[i]+0.1*t,Estatus[i]),rotation=-trans_angle,color='r')

      plt.plot((0,time[-1]),(Estatus[-1],Estatus[-1]),':r',linewidth=1.0,alpha=0.3)
      # plt.subplot(211)
      # plt.legend()
#      plt.subplot(212)
#      plt.legend()

#====================================================================
# show/save profile
#====================================================================
      
      self.pdf.savefig()

#====================================================================
# function to identify segment type
#====================================================================

   def segment_type(self,segid):

      m              = self.mission
      typ            = m['segment_type'][segid].rstrip()
      if(typ == 'all'):
         style       = '-r'
      elif(typ == 'reserve'):
         style       = '-r'
      else:
         quit('unknown segment type')

      return typ, style

#====================================================================
# driver function to render plots
#====================================================================

   def render_plots(self):
      fname       = 'battery_design_'+str(self.id)+'.pdf'
      with PdfPages(fname) as pdf:

         self.pdf  = pdf 
         self.current_battery()
         print('saved battery draw profile to '+fname)