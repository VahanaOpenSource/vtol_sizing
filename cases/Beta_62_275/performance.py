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
from numpy import cos,sin,tan,sqrt

kts2fps = 1.68781
fps2kts = 1.0/kts2fps
fps2nm  = 0.000164579
r2d     = 180.0/np.pi
d2r     = 1.0/r2d 
font    = {'weight' : 'normal',
           'size'   : 18}
matplotlib.rc('font',**font)

class extract_performance:
   
# ===================================================================
# initialize 
# ===================================================================

   def __init__(self,filename):

      if not os.path.isfile(filename):
         print '\n\nFile:',filename,'does not exist.\n\n'
         sys.exit(1)

      with open(filename) as f:
         fyaml=yaml.load(f)

# ===================================================================
# basic dictionary pointers and parameters
# ===================================================================

      self.rotor        = fyaml['rotor']
      self.wing         = fyaml['wing']['set0']
      self.vehicle      = fyaml['vehicle']
      self.mission      = fyaml['mission']
      self.sfc_h        = self.mission['sfc_segment'][0]
      self.sfc_c        = self.mission['sfc_segment'][1]
      self.gtow         = self.vehicle['gtow'][0]
      self.f_factor     = self.vehicle['f_factor'][0]
      self.Vcruise      = max(self.mission['cruise_speed'])       # knots  
      self.Vthresh      = 0.65*self.Vcruise                        # knots 
      self.Vthresh2     = self.Vcruise-10 

# ===================================================================
# rotor parameters
# ===================================================================

      self.rotor_radius = self.rotor['radius'][0]
      self.rotor_area   = np.pi * (self.rotor_radius)**2.0
      self.sigma        = self.rotor['solidity'][0]
      self.cd0          = self.rotor['cd0'][0]
      self.ipf          = self.rotor['ipf'][0]
      self.hvrdwld      = self.rotor['hvr_dwld'][0]
      self.vtip         = self.rotor['tip_speed'][0]
      self.rpm_ratio    = self.rotor['cruise_rpm_ratio'][0]

# ===================================================================
# wing parameters
# ===================================================================

      self.nwing        = self.wing['nwing'][0]
      self.wing_area    = self.wing['span'][0] * self.wing['chord'][0]
      self.oswald       = self.wing['oswald'][0]
      self.wing_ar      = self.wing['aspect_ratio'][0]
      self.cd0w         = self.wing['cd0'][0]
      self.FM           = self.rotor['hover_FM'][0]
      self.eta_xmsn     = self.rotor['eta_xmsn'][0]

# ===================================================================
# get density in slug/cu.ft
# ===================================================================

      self.all_rho      = np.asarray(self.mission['density'])*0.002378/1.2256

# ===================================================================
# propeller and wing 
# ===================================================================

      self.prop_present = False 
      self.wing_present = False 
      try:
         self.prop         = fyaml['prop']
         self.eta_p        = self.prop['eta'][0]
         self.prop_present = True
      except:
         pass 
         print 'propeller not present in system'
   
      if self.wing_area * self.nwing > 0:
         self.wing_present = True 
         ARe               = self.wing_ar * self.oswald
         wing_Cla          = ARe/(ARe+2.0)*2*np.pi          # lift-curve slope, /rad
         self.wing_aoa     = self.wing['cl'][0]/wing_Cla    # in /rad
         self.wing_cla     = wing_Cla                       # /rad
         self.alpha_max    = 0.85/self.wing_cla             # max angle of attack, rad
         self.alpha_max    = self.alpha_max * r2d           # max angle of attack, deg

# ===================================================================
# configuration
# ===================================================================

      self.aircraftID   = self.vehicle['aircraftID'][0]
      if self.aircraftID == 2:            # tilt something
         if self.hvrdwld < 0.1:
            self.tilt_type = 'tilt_wing'
         else: 
            self.tilt_type = 'tilt_rotor'
      elif self.aircraftID == 5:          # tail-sitter
         self.tilt_type    = 'tilt_body'
      else:                               # regular or compound
         self.tilt_type = 'none'

# ===================================================================
# available power at rotor is stored in a parameter called "pwr_installed" (units are hp)
# ===================================================================

      self.pwr_installed   = self.eta_xmsn*self.vehicle['power_installed'][0]
      self.total_energy    = self.vehicle['total_energy'][0] # hp

      self.nrotor          = self.rotor['nrotor'][0]
      self.BIGVALUE        = 1e15

      self.get_airfoil_data()

# ===================================================================
# get airfoil data from NACA 0012 tables (for testing purposes)
# ===================================================================

   def get_airfoil_data(self):
      
      kk = 0
      with open('Postprocessing/table.data','r') as affile:
         for lid,line in enumerate(affile):
            if lid == 0:
               continue      
            temp = line.rstrip('\n').split()
      
            if lid == 1:
               naoa              = int(temp[0])
               self.aoadata      = np.zeros(naoa,dtype=float)
               self.cldata       = np.zeros(naoa,dtype=float)
               self.cddata       = np.zeros(naoa,dtype=float)
      
            if lid > 1:
               self.aoadata[kk]  = float(temp[0])
               self.cldata[kk]   = float(temp[1])
               self.cddata[kk]   = float(temp[2])
               kk += 1

      self.clmax     = max(self.cldata)

#====================================================================
# find wing loads
# alpha_w = wing angle of incidence, deg
#====================================================================

   def wing_loads(self, alpha_w, vinf, rho):

#dynamic pressure
      q              = 0.5*rho*vinf*vinf 

#wing properties
      qS             = q * self.wing_area
      nwing          = self.nwing

# get cl and cd from tables, implement 3d correction
      cl_wing        = np.interp(alpha_w, self.aoadata, self.cldata)
      cd0_wing       = np.interp(alpha_w, self.aoadata, self.cddata)

#implement 3d correction for wing 
      ARe            = self.wing_ar * self.oswald
      K              = 1.0/(np.pi*ARe)
      cl_wing        = cl_wing  * ARe/(2.0+ARe)         
      cd_wing        = cd0_wing + cl_wing*cl_wing *K

#====================================================================
# lift and drag of wing
#====================================================================

      Lw             = cl_wing * qS * nwing 
      Dw             = cd_wing * qS * nwing 

      return Lw, Dw 

#====================================================================
# This function calculates body pitch attitudes at a flight condition
# vinf in ft/s
# alpha_r = rotor angle of attack
# alpha_w = wing angle of attack
#====================================================================

   def trim_solve(self, vinf, weight, rho, alpha_guess): 

      converged      = False
      maxcount       = 5000
      error          = 100.0
      kk             = 0

      q              = 0.5*rho*vinf*vinf 

#====================================================================
# incoming weight does not take into account hover download factor
# add extra thrust up to mu = 0.1, linearly vary down from mu=0.05 
#====================================================================

      mu                = vinf/self.vtip
      if mu < 0.05:
         over_T         = self.hvrdwld
      elif mu < 0.1:
         over_T         = self.hvrdwld*(0.1-mu)/0.05
      else:
         over_T         = 0.0

#      print over_T
#total body forces
      Fzreq             = weight*(1.0+over_T)
      Fxreq             = q * self.f_factor 

      if self.prop_present:
         soln_type         = 'analytical'     # propeller: keep vehicle flat
      else:                                # 
         if self.wing_present:
            if self.tilt_type == 'tilt_rotor':
               soln_type   = 'analytical'
            else:
               soln_type   = 'iterative'
         else:
            soln_type      = 'analytical'

         option      = 'iterate'    
#         option      = 'sweep'    

#====================================================================
# Extract wing properties if present
#====================================================================

      if self.wing_present:
         if vinf > 0:
            CLreq          = Fzreq/(q*self.wing_area*self.nwing)
            alpha_w        = CLreq/(self.wing_cla)*r2d

            if alpha_w > self.alpha_max:
               alpha_w     = self.alpha_max 
            CLw            = self.wing_cla * (alpha_w * d2r)

            ARe            = self.wing_ar * self.oswald
            K              = 1.0/(np.pi*ARe)
            CDw            = self.cd0w + K*CLw*CLw 

            Lw             = q*self.wing_area*CLw 
            Dw             = q*self.wing_area*CDw  

         else:
            alpha_w        = 0.0
            Lw             = 0.0
            Dw             = 0.0
         # print vinf, Fzreq, alpha_w, Lw
      else:
         Lw             = 0.0
         Dw             = 0.0      
      
#====================================================================
# Analytical solutions
# we know the wing loads and fuselage drag/gravity, so find rotor 
# thrust, shaft tilt and calculate power. easy peasy
#====================================================================

      if soln_type == 'analytical':
         Fxreq       = Fxreq + Dw
         Fzreq       = Fzreq - Lw
         soln_found  = True 
         if self.prop_present:
            Pp       = Fxreq * vinf/(self.eta_p*550.0)
            Fxreq    = 0.0 
            T        = Fzreq
            alpha    = 0.0
         else:           
            T        = sqrt(Fxreq*Fxreq + Fzreq*Fzreq)
            alpha    = math.atan2(Fxreq,Fzreq)         # in radians
            Pp       = 0.0


#         print vinf/kts2fps, T, mu*tan(alpha)

         Pr          = self.rotor_power(T, vinf, rho, alpha)

         P           = Pr + Pp

#         print ('%8.3f %8.3f %8.3f' % (vinf/kts2fps, Fxreq, Dw))

#====================================================================
# sweep all angles of attack: tilt-wing or tilt-body
# also account for downwash from rotor on wing changing the aoa
# take simple uniform inflow blowing over whole wing with extra flow
# oriented along the shaft, pointing down/back
#====================================================================

      else:

         alpha       = alpha_guess

         if option == 'sweep':
            err         = 200.0
            iters       = 0
            nal         = 96
            alpha_arr   = np.linspace(-5,90,nal)         # -5 to 90, steps of 1 deg
            dal_flap    = np.asarray([0])
            nflap       = len(dal_flap)
            res         = np.zeros_like(alpha_arr)           # -5 to 90, steps of 1 deg
            soln_found  = False 
            # print 'trying at airspeed = ',vinf/kts2fps,alpha_guess
            conv        = 0.2
            ial         = -1
            for ial in xrange(nal):       # loop over shaft tilt

               alpha    = alpha_arr[ial]
               alpha_w0 = 90.0 - alpha_arr[ial] + self.wing_aoa*r2d      

   #loop over flap settings at this wing tilt angle
               for iflap in xrange(nflap):
                  alpha_w  = alpha_w0 + dal_flap[iflap]
                  Lw, Dw   = self.wing_loads(alpha_w, vinf, rho)

   # find force along X, Z reqd from rotor
                  FxR      = Fxreq + Dw
                  FzR      = Fzreq - Lw

   # update rotor TPP orientation angle
                  al_new   = math.atan2(FxR, FzR)
                  dres     = alpha - al_new*r2d


                  if iflap == 0:
                     minres  = dres
                     minflap = dal_flap[iflap]

                  if abs(dres) <= minres:
                     minres  = dres
                     minflap = dal_flap[iflap]

               res[ial]    = minres

   # y axis: residual (alpha reqd - alpha assumed)
   # x axis: alpha assumed
   # when residual function crosses y=0, its indicated by product of residuals 
   # at left and right ends of interval being non-positive
   # then find the x value at which y = 0 

               print 'rotor angle wrt vertical is ',al_new
               print ial,res[ial]
#               temp = raw_input('yoyo')

               if ial > 0:
                  if res[ial]*res[ial-1] <= 0:
                     soln_found = True 
                     xl         = alpha_arr[ial-1];  yl = res[ial-1]
                     xr         = alpha_arr[ial]  ;  yr = res[ial]
   #find slope of line
                     slope      = (yr-yl) / (xr-xl)

   # line eqn is (y-y0) = slope*(x-x0)
   # at y = 0, x = y0/slope + x0
                     if slope == 0:          # cannot happen, but what the hell..
                        alpha   = xr 
                     else:
                        alpha   = yr/slope + xr

                     print 'V = ',vinf/kts2fps,'found solution at alpha = ', alpha
                     break 

   #with the right shaft tilt, find rotor power
            if soln_found:
               alpha_w  = 90.0 - alpha + self.wing_aoa*r2d + minflap
               Lw, Dw   = self.wing_loads(alpha_w, vinf, rho)

   # find force along X, Z reqd from rotor
               FxR      = Fxreq + Dw
               FzR      = Fzreq - Lw
               T           = sqrt(FxR*FxR + FzR*FzR)
               P           = self.rotor_power(T, vinf, rho, alpha*d2r)
            else:
               P           = 0.0
               print 'WARNING: DID NOT FIND A SOLUTION: ARE YOU KIDDING ME'
               print 'V = ',vinf/kts2fps

#====================================================================
# iterative solution for aoa
#====================================================================

         else: #option == 'iterate':
            err         = 200.0
            soln_found  = False 
            alpha       = alpha_guess

            while not soln_found:
               alpha_w0 = 90.0 - alpha + self.wing_aoa*r2d      
#loop over flap settings at this wing tilt angle
               Lw, Dw   = self.wing_loads(alpha_w, vinf, rho)

# find force along X, Z reqd from rotor
               FxR      = Fxreq + Dw
               FzR      = Fzreq - Lw

# update rotor TPP orientation angle
               al_new   = math.atan2(FxR, FzR)*r2d
               alpha    = 0.95*al_new + 0.05*alpha 

               if abs(al_new - alpha) < 1e-4:
                  soln_found = True 

            T           = sqrt(FxR*FxR + FzR*FzR)
            P           = self.rotor_power(T, vinf, rho, alpha*d2r)

# print FxR, FzR,al_new,err
# print vinf,P, alpha
# x = raw_input('ok?')

      alpha = alpha*d2r
      return P, alpha, soln_found

#====================================================================
# inflow convergence using fixed-point iterations
# mu = advance ratio, alpha = fwd tilt in radians, CT = thrust coeff
#====================================================================

   def converge_inflow(self, mu, mutanal, CT):

#      mutanal  = mu*tan(alpha)
      err      = 20
      musq     = mu*mu
      lam      = 0.05               # initial guess
      while (err > 1e-5):
         lam2  = mutanal + CT*0.5/sqrt(musq + lam*lam) 
         err   = abs(lam2-lam)
         lam   = lam2*0.95 + lam*0.05

      return lam

#====================================================================
# find rotor power (assuming rotor_alpha = 0)
# vinf = ft/s, alpha = fwd tpp tilt (radians), total_thrust (lbs)
#====================================================================

   def rotor_power(self, total_thrust, vinf, rho, alpha):

#====================================================================
# low speed: full rotor RPM
#====================================================================

      A           = self.rotor_area                  # disk area, sq.ft
      V1          = self.Vthresh*kts2fps 
      V2          = self.Vthresh2*kts2fps
#      print self.Vthresh, self.Vcruise
#      quit()
      if vinf < V1:
         factor   = 1.0 
      elif vinf <= V2:
         factor   = 1.0 - (1.0-self.rpm_ratio)*(vinf-V1)/(V2-V1)
      else:
         factor   = self.rpm_ratio 
#      print 'airspeed is ',vinf/kts2fps 
#      yemp = raw_input('yo?') 
#      print 'rotor rpm ratio is ',factor 
      vtip        = self.vtip*factor                 # Vtip at this vinf
      mu          = vinf*cos(alpha)/vtip             # advance ratio
      Tnd         = rho*A*vtip*vtip                  # normalized thrust
      Pnd         = rho*A*vtip*vtip*vtip/550.0       # to convert Cp to hp
      CT          = total_thrust/(self.nrotor*Tnd)   # thrust coeff.

#converge inflow (total inflow)
      mutanal     = vinf*sin(alpha)/vtip 
      lam         = self.converge_inflow(mu, mutanal, CT) 

      ald         = alpha*r2d
      if ald < 10.0:
         ipf      = self.ipf 
      elif ald > 90:
         ipf      = 1.0       # no induced losses
      else:                   # between 10 and 100
         ipf      = 1.0 + (self.ipf-1)*(90.0-ald)/80.0
         # print (100-ald)/90.0        
#calculate rotor induced and profile power 
      cpi         = self.ipf*CT*lam*self.nrotor
      Pi          = cpi*Pnd
      cp0         = 0.125*self.sigma*self.cd0*(1.0+4.65*mu*mu)
      Po          = self.nrotor*cp0*Pnd 
      Ptotal      = Pi + Po
      # print total_thrust/2/A
      # print 'FM = ',(CT*sqrt(0.5*CT))/(lam*CT + cp0)
      # print ('%8.3f %8.3f %8.5f %8.5f %8.2f' % (vinf/kts2fps, rho, CT*sqrt(0.5*CT)*ipf*Pnd*self.nrotor, Pi, Po))
      return Ptotal

# ===================================================================
# power curve
# ===================================================================
   
   def power_curve(self, flt_id):
      
#====================================================================
# forward flight speed-range
#====================================================================

      vmin              = 0.0
      vmax              = 1.05*self.Vcruise * kts2fps # knots to ft/s
#      vmax              = 200 * kts2fps # knots to ft/s
      # print 'VMAX = ',self.Vcruise 
      # x = raw_input('yoe?')
      if vmax < 80*kts2fps:
         vmax           = 80*kts2fps
      nvinf             = 80
      vinf_range        = np.linspace(vmin,vmax,nvinf)

      gtow              = self.vehicle['gtow'][0]              # lbs
      payload           = self.vehicle['payload'][0]           # lbs

      weights           = np.asarray(self.mission['segment_wt'])   # in lbs
      nflt              = len(flt_id)
      indices           = flt_id 
#      tags              = ['take-off','loiter']
      tags              = []
#====================================================================
# allocate memory for power components (energy method)
#====================================================================

      Power             = np.zeros((nvinf, nflt))
      alpha_s           = np.zeros((nvinf, nflt))
      LbyD              = np.zeros((nvinf, nflt))
      xaxis             = np.zeros((nvinf, nflt))

#====================================================================
# loop over densities, weights and speeds
#====================================================================

      for iw in xrange(nflt):
         idx            = indices[iw]
         rho            = self.all_rho[idx]       # density in slug/cu.ft
         weight         = weights[idx]
         vmin           = 0.0 
#         vmax           = 1.1*self.mission['cruise_speed'][idx]*kts2fps
#         if vmax == 0.0:
#            vmax        = 1.1* self.Vcruise * kts2fps 
         print vmax / kts2fps 
         xaxis[:,iw]    = np.linspace(vmin, vmax, nvinf)

         for j in range(nvinf):
            vinf        = xaxis[j,iw]

#====================================================================
# find trim solution for a given flight condition
#====================================================================

            if j == 0:
               alpha_guess    = 0.0 
            else:
               alpha_guess    = alpha_s[j-1,iw]
            P,als,valid       = self.trim_solve(vinf, weight, rho, alpha_guess)
            if not valid:
               P              = Power[j-1,iw]
               als            = alpha_s[j-1,iw]
               xaxis[j,iw]    = xaxis[j-1,iw]
            Power[j,iw]       = P
            LbyD[j,iw]        = weight*vinf/(P*550.0)
            alpha_s[j,iw]     = als*r2d 

#====================================================================
# find cruise point, mark it with red circle
#====================================================================

         Pc,alc,valid         = self.trim_solve(self.Vcruise*kts2fps, weight, rho, alpha_guess)

#====================================================================
# update plots
#====================================================================

         plt.figure(2)
         plt.plot(xaxis[:,iw]*fps2kts,Power[:,iw],linewidth=2.4)#label=tags[iw])
         if iw == 0:
            print 'hello',self.Vcruise, Pc
            plt.plot(self.Vcruise,Pc,'ro',markersize=12,markerfacecolor='None',markeredgewidth=2.2)
   
         plt.figure(3)
         plt.plot(xaxis[:,iw]*fps2kts,alpha_s[:,iw],linewidth=2.4)#label=tags[iw])
         if iw == 0:
            plt.plot(self.Vcruise,als*r2d,'ro',markersize=12,markerfacecolor='None',markeredgewidth=2.2)

         plt.figure(4)
         plt.plot(xaxis[:,iw]*fps2kts, LbyD[:,iw],linewidth=2.4)#label=tags[iw]
         if iw == 0:
            plt.plot(self.Vcruise,weight*self.Vcruise*kts2fps/Pc/550.0,'ro',markerfacecolor='None',markersize=12,markeredgewidth=2.2)

#====================================================================
# total power
#====================================================================

      plt.figure(2)
      plt.plot([vmin/kts2fps,vmax/kts2fps],[self.pwr_installed,self.pwr_installed],'k-.',linewidth=3.0)
      plt.xlabel('Cruise speed, knots')
      plt.ylabel('Power required, hp')
      plt.legend(loc='best')
      plt.grid(linestyle='--')
      xmax     = round(vmax/kts2fps,-1)
      plt.xlim([0.0,xmax])
      plt.ylim([0.0,1.2*self.pwr_installed])
      plt.legend()
      plt.tight_layout()
      plt.savefig('power_curve.png',format='png',dpi=150)
      #plt.title('Power curve')
      plt.figure(3)
      plt.xlabel('Cruise speed, knots')
      plt.ylabel('Shaft tilt, degrees')
      plt.grid(linestyle='--')
      plt.xlim([0.0,xmax])
      plt.ylim([0.0,100])
      plt.tight_layout()
      plt.savefig('alpha.png',format='png',dpi=150)

#====================================================================
# L/D plot: W*V/P
#====================================================================

      plt.figure(4)
#      for iw in xrange(nflt):
#         idx            = indices[iw]
#         rho            = self.all_rho[idx]       # density in slug/cu.ft
#         weight         = weights[idx]
      plt.xlabel('Cruise speed, knots')
      plt.ylabel('Lift to Drag ratio')
      plt.grid(linestyle='--')
      # print xmax, 
      # plt.xlim([0.0,xmax])
      plt.ylim([0.0,np.ceil(np.amax(LbyD))])
      plt.legend(loc='best')
      plt.tight_layout()
      plt.savefig('LbyD.png',format='png',dpi=150)

#      quit('OK POWER CURVE?')
#====================================================================
# plot energy consumed with time
#====================================================================

   def energy_consumption(self):

      nsegments         = self.mission['nsegments'][0]
      energy_consumed   = np.zeros(nsegments,dtype='float')

      for n in range(nsegments):
         energy_segment = self.mission['energy_segment'][n]
         time_segment   = self.mission['time'][0] / 60.0 # convert from min to hr

         if n == 0:
            energy_consumed[n] = energy_segment
         else:
            energy_consumed[n] = energy_consumed[n-1]+energy_segment

      seg_range         = np.arange(nsegments)+1

      plt.figure()
      plt.plot([0,nsegments],[self.total_energy,self.total_energy],'k--',linewidth=2.5,label='Total energy')
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

      rho            = self.all_rho[0]
# weight minus the payload = fuel weight
      try:
         mfuel       = self.vehicle['fuel'][0]
      except:
         mfuel       = self.vehicle['battery'][0]

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

         weight               = (empty_wt + payload[ifuel] + fuel[ifuel])#*(1.0+self.hvrdwld)
         power,alpha_s,valid  = self.trim_solve(0.0, self.gtow, rho, 0.0)

         thover[ifuel]        = 0.0
         if (power < self.pwr_installed):
            thover[ifuel]     = fuel[ifuel]/(power*self.sfc_h) * 60.0 # convert from hrs to mins

# ===================================================================
# reference condition: design point
# ===================================================================

      power,alpha_s,valid     = self.trim_solve(0.0, self.gtow, rho, 0.0)
      tref                    = mfuel/(power*self.sfc_h)*60.0 

# ===================================================================
# make a plot
# ===================================================================
         
      plt.figure()
      plt.plot(thover, payload,'-',markevery=3,linewidth=2.2)
      plt.plot(tref,mpay,'ro',markersize=12,markerfacecolor='None',markeredgewidth=2.2)

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

      try:
         mfuel       = self.vehicle['fuel'][0]
      except:
         mfuel       = self.vehicle['battery'][0] 
      mpay           = self.vehicle['payload'][0]
      empty_wt       = self.gtow - mpay - mfuel
      wt             = self.gtow - mpay            # lbs; empty + fuel
      max_fuel       = mfuel*1.35 
      max_payload    = mpay*1.35
      N              = 40
      payload        = np.linspace(0.7*mpay,mpay+mfuel,N)    # fuel: from zero to full fuel 
      fuel           = np.zeros_like(payload)         # payload, lbs
      thover         = np.zeros_like(payload)         # hover time
      
#====================================================================
# range of vinf
#====================================================================

      V1          = self.Vcruise - 10
      V2          = self.Vcruise + 10
      nV          = 3
      vinf_range  = np.linspace(V1,V2,nV)
      vinf_range  = vinf_range*kts2fps
      styles      = ['-','-','--']
      markers     = ['','s','o']

#====================================================================
# allocate
#====================================================================

      mtx_endu    = np.zeros([N,nV],dtype='float')
      mtx_range   = np.zeros([N,nV],dtype='float')

      for j in range(nV):        # different airspeeds
         vinf = vinf_range[j]
         for k in xrange(N):        # different payloads

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

            weight               = empty_wt + payload[k] + fuel[k]
            alpha_guess          = 0.0 
            rho                  = self.all_rho[1]
            power,alpha_s,valid  = self.trim_solve(vinf, weight, rho, alpha_guess)

            if (power < self.pwr_installed):
               
               t                 = fuel[k]/(power*self.sfc_c) * 60.0 # convert from hrs to mins
               mtx_endu[k,j]     = t                        # time in mins
               mtx_range[k,j]    = vinf * t * fps2nm * 60   # nautical miles

#====================================================================
# also find baseline condition for pure cruise
#====================================================================
      
      vinf                       = self.Vcruise*kts2fps
      power,alpha_s,valid        = self.trim_solve(vinf, self.gtow, rho, alpha_guess)

      t                          = mfuel/(power*self.sfc_c) * 60.0 # convert from hrs to mins
      end_ref                    = t                        # time in mins
      range_ref                  = vinf * t * fps2nm * 60   # nautical miles

#====================================================================
# create plots
#====================================================================

      plt.figure()
      for j in range(nV):
         plt.plot(mtx_endu[:,j],payload,'o-',label=('%4.0f kts' % (vinf_range[j]*fps2kts)),linestyle=styles[j],marker=markers[j],linewidth=2.2,markevery=3)
      plt.plot(end_ref,mpay,'ro',markersize=12,markerfacecolor='None',markeredgewidth=2.2)
      plt.xlabel('Time in cruise, min')
      plt.ylabel('Payload, lb')
      plt.grid(linestyle='--')
      plt.title('Payload-Endurance [Cruise]')
      plt.legend(loc='best',fontsize=15)
      plt.tight_layout()
      plt.savefig('cruise_endurance.png',format='png',dpi=150)
      plt.figure()
      for j in range(nV):
         plt.plot(mtx_range[:,j],payload,'o-',label=('%4.0f kts' % (vinf_range[j]*fps2kts)),linestyle=styles[j],marker=markers[j],linewidth=2.2,markevery=3)
      plt.plot(range_ref,mpay,'ro',markersize=12,markerfacecolor='None',markeredgewidth=2.2)
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

if len(argv) < 2:
   print " -- Error: performance.py requires a design ID"
   print " -- Usage: 'python performance.py <ID> '"

#====================================================================
# get design ID
#====================================================================

n           = int(argv[1])
DIR         = 'output/logs/'

if len(argv) == 3:
   id_flt   = argv[2]
   indices  = id_flt.split(',')
   for ix in xrange(len(indices)):
      indices[ix] = int(indices[ix])
else:
   indices   = [0]

#====================================================================
#if n < 10:
#   fname=DIR+'log00'+str(n)+'.yaml'
#elif n<100:
#   fname=DIR+'log0'+str(n)+'.yaml'
#elif n<1000:
#====================================================================

fname = DIR+'log'+str(n)+'.yaml'
print " -- Chosen design ID",n
perf  = extract_performance(fname)
perf.power_curve(indices)
perf.cruise_performance()
perf.hover_endurance()