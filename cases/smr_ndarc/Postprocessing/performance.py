#
# file to extract performance
#
# power curve
# payload range diagram
# payload endurance diagram

import yaml
import sys,os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

kts2fps = 1.68781
fps2kts = 1.0/kts2fps
fps2nm  = 0.000164579
rho     = 0.002378 # slugs/ft3

font={'weight' : 'normal',
      'size'   : 18}
matplotlib.rc('font',**font)


class extract_performance:
   

   def __init__(self,filename):

      if not os.path.isfile(filename):
         print '\n\nFile:',filename,'does not exist.\n\n'
         sys.exit(1)

      with open(filename) as f:
         fyaml=yaml.load(f)

      # 
      # basic dictionary pointers and parameters
      #
      self.rotor     = fyaml['rotor']
      self.wing      = fyaml['wing']
      self.vehicle   = fyaml['vehicle']
      self.mission   = fyaml['mission']
      self.sfc_h     = self.mission['sfc_segment'][0]
      self.sfc_c     = self.mission['sfc_segment'][1]
      self.gtow      = self.vehicle['gtow'][0]
      self.f_factor  = self.vehicle['f_factor'][0]

      self.Vcruise   = max(self.mission['cruise_speed'])

      self.Vthresh   = 0.7*self.Vcruise*kts2fps

      # rotor parameters
      self.rotor_radius = self.rotor['radius'][0]
      self.rotor_area   = np.pi * (self.rotor_radius)**2.0
      self.sigma        = self.rotor['solidity'][0]
      self.cd0          = self.rotor['cd0'][0]
      self.ipf          = self.rotor['ipf'][0]
      self.hvrdwld      = self.rotor['hvr_dwld'][0]
      self.vtip         = self.rotor['tip_speed'][0]
      self.rpm_ratio    = self.rotor['cruise_rpm_ratio'][0]
      #quit('ok?')
      # wing parameters
      self.wing_area = self.wing['span'][0] * self.wing['chord'][0]
      self.oswald    = self.wing['oswald'][0]
      self.wing_ar   = self.wing['aspect_ratio'][0]
      self.cd0w      = self.wing['cd0'][0]
      self.FM        = self.rotor['hover_FM'][0]
      self.eta_xmsn  = self.rotor['eta_xmsn'][0]

# available power at rotor is stored in a parameter called "pwr_installed" (units are hp)
      self.pwr_installed = self.eta_xmsn*self.vehicle['power_installed'][0]
      # print 'INSTALLED POWER',self.pwr_installed
      # print self.pwr_installed
      # quit('hel')
      self.total_energy  = self.vehicle['total_energy'][0] # hp
      # non-dimensional terms
      self.non_dim_thrust = rho * self.rotor_area * self.vtip**2
      self.non_dim_power  = self.non_dim_thrust * self.vtip

      self.nrotor = self.rotor['nrotor'][0]
      self.nwing  = self.wing['nwing'][0]

      self.BIGVALUE = 1e15

      # read airfoil data
      self.get_airfoil_data()

# ===================================================================
# get airfoil data from NACA 0012 tables (for testing purposes)
# ===================================================================

   def get_airfoil_data(self):
      
      kk=0
      #with open('./Postprocessing/table.data','r') as affile:
      with open('table.data','r') as affile:
         for lid,line in enumerate(affile):
            if lid==0:
               continue
      
            temp=line.rstrip('\n').split()
      
            if lid==1:
               naoa=int(temp[0])
               self.aoadata=np.zeros(naoa,dtype=float)
               self.cldata=np.zeros(naoa,dtype=float)
               self.cddata=np.zeros(naoa,dtype=float)
      
            if lid>1:
               self.aoadata[kk]=float(temp[0])
               self.cldata[kk]=float(temp[1])
               self.cddata[kk]=float(temp[2])
               kk+=1
      
      
      if False:
         plt.figure(100)
         plt.subplot(211)
         plt.plot(self.aoadata,self.cldata,'o-')

         plt.subplot(212)
         plt.plot(self.aoadata,self.cddata,'o-')

         plt.show()

      self.clmax = max(self.cldata)

      self.vstall = np.sqrt( ( 2.0 * self.gtow/self.nwing)/(self.clmax*rho*self.wing_area))



# ===================================================================
# solve for the power at a given trim state
# ===================================================================

   def trim_solve(self,vinf,weight,alpha_guess): # vinf in ft/s

      alpha_old = alpha_guess #0.0 # initial guess
      
      converged=False
      
      maxcount=5000
      
      error=100.0
      kk = 0


      if abs(vinf) < 1.e-3:
         converged=True
      
      #
      # trim convergence loop
      #
      while (not converged):

         # incoming weight does not take into account hover download factor
         # vary from 'ALL' at alpha=0, to '0' at alpha=pi/2
         weight_temp = weight*(1.0 + self.hvrdwld*(1.0 - 2.0*alpha_old/np.pi))
      

         # aoa of wing is the complimentary angle of the shaft-tilt
         aoa_wing = 0.5*np.pi - alpha_old

         # get cl and cd from tables
         cl_wing   = np.interp(aoa_wing/np.pi*180.0, self.aoadata, self.cldata)
         cd0_wing  = np.interp(aoa_wing/np.pi*180.0, self.aoadata, self.cddata)
         
         cd_wing = (cd0_wing + (cl_wing**2.0) /(np.pi * self.wing_ar * self.oswald))

         # lift and drag of wing

         lift_wing = 0.5 * cl_wing * rho * vinf*vinf * self.wing_area
         drag_wing = 0.5 * cd_wing * rho * vinf*vinf * self.wing_area

         # drag of fuselage
         drag_fus  = 0.5 * rho * self.f_factor * vinf*vinf

         # r1 and r2
         r1 = self.nwing*drag_wing + drag_fus
         r2 = weight_temp - self.nwing*lift_wing

         # update aoa
         alpha_new = math.atan2(r1,r2)

         # error
         #if (abs(alpha_old) > 1.e-3):
         error = abs(alpha_new - alpha_old) #/ abs(alpha_old)

         if error < 0.01/180.0*np.pi:
            converged = True

         # update aoa
         omega=0.99
         alpha_old = omega * alpha_old + (1.0-omega) * alpha_new

         kk+=1

         if kk==maxcount:
            break
      
      # -- end while loop

      weight_temp = weight*(1.0 + self.hvrdwld*(1.0 - 2.0*alpha_old/np.pi))

      
      aoa_wing  = 0.5*np.pi - alpha_old
      cl_wing   = np.interp(aoa_wing/np.pi*180.0, self.aoadata, self.cldata)
      cd0_wing  = np.interp(aoa_wing/np.pi*180.0, self.aoadata, self.cddata)
      cd_wing   = (cd0_wing + cl_wing**2/(np.pi * self.wing_ar * self.oswald))
      lift_wing = 0.5 * cl_wing * rho * vinf*vinf * self.wing_area
      drag_wing = 0.5 * cd_wing * rho * vinf*vinf * self.wing_area
      drag_fus  = 0.5 * rho * self.f_factor * vinf*vinf

      r1 = self.nwing*drag_wing + drag_fus
      r2 = weight_temp - self.nwing*lift_wing

      thrust = np.sqrt(r1*r1+r2*r2)
      if (not converged):
         print ' -- trim solver failed at speed ',vinf*fps2kts,'knots'
         #sys.exit(1)
      
      #print '-------------------------------------------'
      #print 'r1,r2     ',r1,r2
      #print 'thrust    ', thrust
      #print 'winglift  ', lift_wing
      #print 'drag_wing ',drag_wing
      #print 'drag_fus  ', drag_fus
      #print 'aoa       ',alpha_old/np.pi * 180.0
      #print 'gtow      ',self.gtow
      #print 't_para    ',thrust*math.sin(alpha_old)
      #print 't_perp    ',thrust*math.cos(alpha_old)

      return thrust, (self.nwing *lift_wing), alpha_old


# ===================================================================
# find rotor power (assuming rotor_alpha = 0)
# ===================================================================

   def find_rotor_power(self,total_thrust,vinf):

      
# ===================================================================
# low speed: full rotor RPM
# ===================================================================

      if vinf < self.Vthresh:
         factor   = 1.0 
      elif vinf <= self.Vcruise*kts2fps:
         factor   = 1.0 - (1.0-self.rpm_ratio)*(vinf-self.Vthresh)/(self.Vcruise*kts2fps-self.Vthresh)
      else:
         factor   = self.rpm_ratio 

#      print vinf, self.Vthresh, self.Vcruise*kts2fps,factor 
      # induced power
      thrust      = total_thrust / self.nrotor
      ct          = thrust/(self.non_dim_thrust*factor*factor)
      lh          = np.sqrt(0.5*ct)
      mu          = vinf/(self.vtip*factor)
      mu_lh       = mu/lh
      l_lh        = np.sqrt( -0.5*mu_lh**2 + np.sqrt(1.0+0.25*mu_lh**4)) 
      cpi         = self.ipf*ct*(l_lh*lh)
      pwr_induced = self.nrotor * cpi * self.non_dim_power / 550.0*factor*factor*factor
      
      # rotor profile power
      alpha       = 0.5*np.pi
      mu_rotor    = np.multiply(vinf,np.cos(alpha)/(self.vtip*factor))
      cp0         = 0.125 * self.sigma * self.cd0 * ( 1.0 + 4.65 * mu_rotor**2)
      pwr_profile = self.nrotor*cp0*self.non_dim_power / 550.0 * factor*factor*factor # convert to HP

      return pwr_induced,pwr_profile

# ===================================================================
# wing related powers
# ===================================================================
   
   def find_wing_power(self,total_lift,vinf, aoa_wing):


      if self.nwing==0 or abs(total_lift < 1.e-5):
         return 0,0,0
      
      cl   = np.interp(aoa_wing/np.pi*180.0, self.aoadata, self.cldata)
      cd0w = np.interp(aoa_wing/np.pi*180.0, self.aoadata, self.cddata)
#
# wing induced power
#
#lift = total_lift / self.nwing
#cl   = lift/(0.5*rho*self.wing_area*vinf*vinf)
      cdi  = cl* cl /(np.pi * self.wing_ar * self.oswald)
      
      pwr_induced = self.nwing * cdi * 0.5 * rho * vinf**3 * self.wing_area / 550.0
# 
# wing profile power
#
      pwr_profile = self.nwing * self.cd0w * 0.5 * rho * vinf**3 * self.wing_area / 550.0
# 
# wing drag
#
      wing_drag   = (pwr_induced+pwr_profile) * 550.0 /vinf

      return pwr_induced,pwr_profile,wing_drag

# ===================================================================
# fuselage parasitic power
# ===================================================================

   def find_parasitic_power(self,vinf,alphashaft):

      mu             = vinf/self.vtip
      cpp            = 0.5 * self.f_factor / self.rotor_area * mu**3
      pwr_parasitic  = cpp*self.non_dim_power / 550.0 # convert to HP

      return pwr_parasitic

# ===================================================================
# power curve
# ===================================================================
   
   def power_curve(self):
      
      #
      # forward flight speed-range
      # 
      vmin              = 0.0
      vmax              = 1.2*self.Vcruise * kts2fps # knots to ft/s

      if vmax < 100*kts2fps:
         vmax           = 100*kts2fps
      nvinf             = 100
      vinf_range        = np.linspace(vmin,vmax,nvinf)
      #vinf_range=np.array([19.9,19.99]) * kts2fps
      #nvinf=2


      gtow              = self.vehicle['gtow'][0]           # lbs
      payload           = self.vehicle['payload'][0]        # lbs

      weight_array      = np.asarray([gtow, gtow-payload])
      nwts              = np.shape(weight_array)[0]

# allocate memory for power components (energy method)
      pwr_rotor_induced = np.zeros((nvinf, nwts))
      pwr_rotor_profile = np.zeros((nvinf, nwts))
      pwr_wing_induced  = np.zeros((nvinf, nwts))
      pwr_wing_profile  = np.zeros((nvinf, nwts))
      pwr_parasitic     = np.zeros((nvinf, nwts))
      alpha_shaft_tilt  = np.zeros((nvinf, nwts))

      for iw in xrange(nwts):
         weight         = weight_array[iw]
         for j in range(nvinf):
            vinf        = vinf_range[j]

            alpha_guess = 0.0 if j==0 else alpha_shaft_tilt[j-1,iw]

#            print vinf, weight, alpha_guess
            [rotor_thrust,wing_lift,alphashaft] = self.trim_solve(vinf,weight,alpha_guess)

            alpha_shaft_tilt[j,iw]  = alphashaft

            wing_aoa = 0.5*np.pi - alphashaft

            [pind1,ppro1]           = self.find_rotor_power(rotor_thrust,vinf)
            [pind2,ppro2,wdrag]     = self.find_wing_power(wing_lift,vinf,wing_aoa)
            pwr_parasitic[j,iw]     = self.find_parasitic_power(vinf,alphashaft)

#            print vinf*fps2kts,pind1+ppro1+pind2+ppro2,pwr_parasitic[j,iw]
            pwr_rotor_induced[j,iw] = pind1
            pwr_rotor_profile[j,iw] = ppro1

            pwr_wing_induced[j,iw]  = pind2
            pwr_wing_profile[j,iw]  = ppro2

      # -- end j loop

      # 
      # total power
      # 
      pwr_total=(pwr_rotor_induced +  
                 pwr_rotor_profile + 
                 pwr_wing_induced  + 
                 pwr_wing_profile  +
                 pwr_parasitic)

      plt.figure()
      plt.plot([vmin,vmax],[self.pwr_installed,self.pwr_installed],'k--',linewidth=3.0)
#      plt.plot(vinf_range*fps2kts,pwr_parasitic,linewidth=3.0,label='parasitic',linestyle='--')
#      plt.plot(vinf_range*fps2kts,pwr_rotor_induced,linewidth=3.0,label='rotor induced',linestyle='--')
#      plt.plot(vinf_range*fps2kts,pwr_rotor_profile,linewidth=3.0,label='rotor profile',linestyle='--')
#      if self.nwing>0:
#         plt.plot(vinf_range*fps2kts,pwr_wing_induced,linewidth=3.0,label='wing induced',linestyle='--')
#         plt.plot(vinf_range*fps2kts,pwr_wing_profile,linewidth=3.0,label='wing profile',linestyle='--')
#      plt.plot(vinf_range*fps2kts,pwr_total,linewidth=3.0,color='k',label='total')
      plt.plot(vinf_range*fps2kts,pwr_total,linewidth=3.0)
      plt.xlabel('Cruise speed, knots')
      plt.ylabel('Power required, hp')
      plt.grid(linestyle='--')
      plt.xlim([0.0,vmax/kts2fps])
      plt.ylim([0.0,1.2*self.pwr_installed])
      plt.legend()
      plt.tight_layout()
      plt.savefig('power.png',format='png',dpi=150)
      #plt.title('Power curve')
      plt.figure()
      plt.plot(vinf_range*fps2kts,alpha_shaft_tilt/np.pi*180.0,'o-')
      plt.xlabel('Cruise speed, knots')
      plt.ylabel('Shaft tilt, degrees')
      plt.grid(linestyle='--')
      plt.xlim([0.0,vmax/kts2fps])
      plt.ylim([0.0,90])
      plt.tight_layout()
      plt.savefig('alpha.png',format='png',dpi=150)

# ===================================================================
# plot energy consumed with time
# ===================================================================

   def energy_consumption(self):

      nsegments       = self.mission['nsegments'][0]

      energy_consumed = np.zeros(nsegments,dtype='float')


      for n in range(nsegments):
         energy_segment = self.mission['energy_segment'][n]
         time_segment   = self.mission['time'][0] / 60.0 # convert from min to hr

         if n==0:
            energy_consumed[n]=energy_segment
         else:
            energy_consumed[n]=energy_consumed[n-1]+energy_segment

      
      seg_range=np.arange(nsegments)+1

      plt.figure()
      plt.plot([0,nsegments],[self.total_energy,self.total_energy],'k--',linewidth=2,label='Total energy')
      plt.bar(seg_range,energy_consumed,color='C2')
      plt.xticks(seg_range)
      plt.xlabel('Segment ID')
      plt.ylabel('Energy, hp-hr')
      plt.grid(linestyle='--')
      plt.legend()
      plt.title('Energy consumption')


# ===================================================================
# plot hover endurance  with payload
# ===================================================================

   def hover_endurance(self):

# weight minus the payload = fuel weight
      mfuel          = self.vehicle['fuel'][0]
      mpay           = self.vehicle['payload'][0]
      empty_wt       = self.gtow - mpay - mfuel
      wt             = self.gtow - mpay            # lbs; empty + fuel
      max_fuel       = mfuel*1.35 
      max_payload    = mpay*1.35

#      print self.gtow,empty_wt,mpay,mfuel
      N              = 40
      payload        = np.linspace(0.7*mpay,mpay+mfuel,N)    # fuel: from zero to full fuel 
      fuel           = np.zeros_like(payload)         # payload, lbs
      thover         = np.zeros_like(payload)         # hover time

# ===================================================================
# calculate payload assuming payload + fuel = constant, ie GTOW = constant
# sfc          = lb/hp-hr 
# mass fuel    = lbs 
# power        = hp
#
#
# sfc          = (mass fuel)/(power*time)
#
# so time      = fuel / (power * sfc)
#
# ===================================================================

      for ifuel in xrange(N):
         fuel[ifuel]    = mpay + mfuel - payload[ifuel]

# ===================================================================
# limit payload and fuel to 135% of design values
# ===================================================================

         if fuel[ifuel] > max_fuel:
            fuel[ifuel] = max_fuel

         if payload[ifuel] > max_payload:
            payload[ifuel] = max_payload

         mpayload       = payload[ifuel]

# ===================================================================
# calculate thrust, power and endurance
# ===================================================================

         thrust         = (empty_wt + payload[ifuel] + fuel[ifuel])*(1.0+self.hvrdwld)
         thrust         = thrust/self.nrotor
         power          = thrust**(1.5)/np.sqrt(2.0*rho*self.rotor_area)/self.FM/550.0
         power          = power*self.nrotor

         thover[ifuel]  = 0.0
         if (power < self.pwr_installed):
            thover[ifuel]= fuel[ifuel]/(power*self.sfc_h) * 60.0 # convert from hrs to mins

# ===================================================================
# make a plot
# ===================================================================
         
      plt.figure()
      plt.plot(thover, payload,'o-')
      plt.ylabel('Payload, lb')
      plt.xlabel('Time in hover, min')
      plt.grid(linestyle='--')
      plt.title('Payload-Endurance [Hover]')
      plt.tight_layout()
      plt.savefig('hover_endurance.png',format='png',dpi=150)
#====================================================================
# plot payload-endurance and payload-range with speed in cruise
#
# payload_range contains the payload range data
# vinf_ff contains ff speed data
#====================================================================
   
   def cruise_performance(self):

#====================================================================
# find empty weight, payload and fuel weights at design point
#====================================================================

      mfuel          = self.vehicle['fuel'][0]
      mpay           = self.vehicle['payload'][0]
      empty_wt       = self.gtow - mpay - mfuel
      wt             = self.gtow - mpay            # lbs; empty + fuel
      max_fuel       = mfuel*1.15 
      max_payload    = mpay*1.15
      N              = 40
      payload        = np.linspace(0.7*mpay,mpay+mfuel,N)    # fuel: from zero to full fuel 
      fuel           = np.zeros_like(payload)         # payload, lbs
      thover         = np.zeros_like(payload)         # hover time
      
# range of vinf
      V1          = self.Vcruise - 40
      V2          = self.Vcruise + 20
      nV          = 5
      vinf_range  = np.linspace(V1,V2,nV)
      vinf_range  = vinf_range*kts2fps
      nvinf       = len(vinf_range)

# allocate
      mtx_endu    = np.zeros([N,nvinf],dtype='float')
      mtx_range   = np.zeros([N,nvinf],dtype='float')

# loops...
      for j in range(nvinf):
         vinf = vinf_range[j]
         for k in xrange(N):

# ===================================================================
# limit payload and fuel to 135% of design values
# ===================================================================

            fuel[k]    = mpay + mfuel - payload[k]
            if fuel[k] > max_fuel:
               fuel[k] = max_fuel

            if payload[k] > max_payload:
               payload[k] = max_payload

# ===================================================================
# calculate thrust, power and endurance
# ===================================================================

            vehicle_weight = empty_wt + payload[k] + fuel[k]

            alpha_guess = 0.0 

            [rotor_thrust,wing_lift,alphashaft] = self.trim_solve(vinf,vehicle_weight,alpha_guess)

            wing_aoa    = 0.5*np.pi - alphashaft

            [pind1,ppro1]        = self.find_rotor_power(rotor_thrust,vinf)
            [pind2,ppro2,wdrag]  = self.find_wing_power(wing_lift,vinf,wing_aoa)
            para                 = self.find_parasitic_power(vinf,alphashaft)

            power = pind1+ppro1+pind2+ppro2+para
            
            if (power < self.pwr_installed):
               
               t              = fuel[k]/(power*self.sfc_c) * 60.0 # convert from hrs to mins
               mtx_endu[k,j]  = t                        # time in mins
               mtx_range[k,j] = vinf * t * fps2nm * 60   # nautical miles

      
      # plot thingies
      plt.figure()
      for j in range(nvinf):
         plt.plot(mtx_endu[:,j],payload,'o-',label=('%4.0f kts' % (vinf_range[j]*fps2kts)))
      plt.xlabel('Time in cruise, min')
      plt.ylabel('Payload, lb')
      plt.grid(linestyle='--')
      plt.title('Payload-Endurance [Cruise]')
      plt.legend(loc='best',fontsize=15)
      plt.tight_layout()
      plt.savefig('cruise_endurance.png',format='png',dpi=150)
      plt.figure()
      for j in range(nvinf):
         plt.plot(mtx_range[:,j],payload,'o-',label=('%4.0f kts' % (vinf_range[j]*fps2kts)))
      plt.xlabel('Range, nm')
      plt.ylabel('Payload, lb')
      plt.grid(linestyle='--')
      plt.title('Payload-Range [Cruise]')
      plt.legend(loc='best',fontsize=15)
      plt.tight_layout()
      plt.savefig('cruise_range.png',format='png',dpi=150)

# ===================================================================
# Main function
# ===================================================================

from sys import argv

if len(argv) != 2:
   print " -- Error: performance.py requires a design ID"
   print " -- Usage: 'python performance.py <ID> '"

# get design ID
n = int(argv[1])

DIR = '../output/logs/'
#if n < 10:
#   fname=DIR+'log00'+str(n)+'.yaml'
#elif n<100:
#   fname=DIR+'log0'+str(n)+'.yaml'
#elif n<1000:
fname=DIR+'log'+str(n)+'.yaml'

print " -- Chosen design ID",n


perf=extract_performance(fname)
perf.power_curve()
perf.cruise_performance()
perf.hover_endurance()
#plt.show()

