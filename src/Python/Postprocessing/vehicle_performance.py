#====================================================================
#
# file to extract performance
#
# power curve
# payload range diagram
# payload endurance diagram
#====================================================================

import yaml, sys, os, matplotlib, numpy
import matplotlib.pyplot as plt
from numpy import pi
#====================================================================
# unit conversions
#====================================================================

grav    = 9.81             # m/s/s
kts2mps = 0.5147           
r2d     = 180.0/numpy.pi 
d2r     = numpy.pi/180.0

#==========================
# use these lines for latex
#==========================

from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages

rc('text',usetex=True)
font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 13}

rc('font', **font)

#====================================================================
# Class definition
#====================================================================

class vehicle_performance:

#====================================================================
# initialization: load data
#====================================================================
   
   def __init__(self, design, table_data):

#====================================================================
# determine tilt type if relevant
#====================================================================
      
      adict                = design.all_dict['aircraft']
      aID                  = adict['aircraftID']

      self.wing_present    = bool(design.wing)
      if aID == 2:            # tilt something
         if design.rotor.groups[0].hvr_dwld < 0.1:
            self.tilt_type = 'tilt_wing'
         else: 
            self.tilt_type = 'tilt_rotor'
      elif aID == 5:          # tail-sitter
         self.tilt_type    = 'tilt_body'
      else:                               # regular or compound
         self.tilt_type    = 'none'
 
      self.design          = design
      self.rho             = design.mission.segment[0].rho 

      self.Vcruise         = 0.0
      mission              = design.mission
      self.tcruise         = 0.0
      self.thover          = 0.0
      for i in range(mission.nseg):
         segment           = mission.segment[i]
         self.Vcruise      = max(self.Vcruise,segment.cruisespeed*kts2mps)
         if(segment.flightmode == 'cruise'):
            self.tcruise   = self.tcruise + segment.time          # cruise time in min
         elif(segment.flightmode == 'hover'):
            self.thover    = self.thover  + segment.time
      self.Vthresh         = 0.7*self.Vcruise
      self.Vthresh2        = self.Vcruise - 10

#====================================================================
# vehicle global parameters
#====================================================================

      self.gtow            = design.massTakeoff*grav 
      self.f_factor        = design.f_plate

#====================================================================
# Installed power line
#====================================================================

      if design.etype == 'electric_motor':
         Eins                 = design.engine.Eins
         self.pwr_installed   = 6*Eins                  # 6C is max "installed" power for battery
      else:
         self.pwr_installed   = design.p_ins

#====================================================================
# read airfoil tables from input file (0012 right now)
#====================================================================

      self.get_airfoil_data(table_data)

# ===================================================================
# get airfoil data from NACA 0012 tables (for testing purposes)
# ===================================================================

   def get_airfoil_data(self, table_data):
      
      kk=0
      with open(table_data,'r') as affile:
         for lid,line in enumerate(affile):
            if lid==0:
               continue
      
            temp=line.rstrip('\n').split()
      
            if lid==1:
               naoa=int(temp[0])
               self.aoadata=numpy.zeros(naoa,dtype=float)
               self.cldata=numpy.zeros(naoa,dtype=float)
               self.cddata=numpy.zeros(naoa,dtype=float)
      
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

      self.clmax  = 1.0        

      design      = self.design

      if(bool(design.wing)):
         self.vstall = design.wing.groups[0].stall_speed
         print('static stall speed is ', self.vstall,' m/s')

#====================================================================
# find wing loads
# alpha_w = wing angle of incidence, deg
#====================================================================

   def wing_loads(self, alpha_w, vinf, rho, group):

#dynamic pressure
      q              = 0.5*rho*vinf*vinf 

#wing properties
      qS             = q * group.area
      nwing          = group.nwings

# get cl and cd from tables, implement 3d correction
      cl_wing        = numpy.interp(alpha_w, self.aoadata, self.cldata)
      cd0_wing       = numpy.interp(alpha_w, self.aoadata, self.cddata)

#implement 3d correction for wing 
      ARe            = group.aspectratio*group.oswald
      K              = 1.0/(numpy.pi*ARe)
      cl_wing        = cl_wing  * ARe/(2.0+ARe)         
      cd_wing        = cd0_wing + cl_wing*cl_wing *K

#====================================================================
# lift and drag of wing
#====================================================================

      Lw             = cl_wing * qS * nwing 
      Dw             = cd_wing * qS * nwing *1.05
      # print('wing CL',cl_wing,Lw)

      return Lw, Dw 

#====================================================================
# This function calculates body pitch attitudes at a flight condition
# vinf in ft/s
# alpha_r = rotor angle of attack
# alpha_w = wing angle of attack
#====================================================================

   def trim_solve(self, vinf, weight, rho, alpha_guess, dalpha_max=1e5): 

      converged      = False
      maxcount       = 5000
      error          = 100.0
      kk             = 0
      q              = 0.5*rho*vinf*vinf 

#====================================================================
# incoming weight does not take into account hover download factor
# add extra thrust up to mu = 0.1, linearly vary down from mu=0.05 
#====================================================================

      design            = self.design 
      rotor             = design.rotor.groups[0]
      Vtip              = rotor.tipspeed        # m/s, hover tip speed
      mu                = vinf/Vtip
      dwld              = rotor.hvr_dwld
      if mu < 0.05:
         over_T         = dwld
      elif mu < 0.1:
         over_T         = dwld*(0.1-mu)/0.05
      else:
         over_T         = 0.0

      Fzreq             = weight*(1.0+over_T)
      Fxreq             = q * self.f_factor 

      adict             = self.design.all_dict['aircraft']
      if adict['npropeller'] > 0:
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

#====================================================================
# Extract wing properties if present
#====================================================================

      Lw             = 0.0; Dw = 0.0
      wing           = self.design.wing
      if self.wing_present:
         if vinf > 0:
            for igroup in range(wing.ngroups):
               group       =wing.groups[igroup]
               CLreq       = Fzreq*group.lift_frac/(q*group.area*group.nwings)
               ARe         = group.aspectratio*group.oswald
               CLa         = 0.11*ARe/(ARe+2.0)
               alpha_w     = CLreq/CLa                # AoA, deg 
               if alpha_w > 12.0:                     # stall AoA, deg
                  alpha_w  = 12.0 
               CLw         = CLa * alpha_w 
               K           = 1.0/(pi*ARe)
               CDw         = (group.cd0 + K*CLw*CLw)*1.05
               Lw          = Lw + q*group.area*CLw*group.nwings
            Dw             = Dw + q*group.area*CDw*group.nwings

         else:
            alpha_w        = 0.0
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
            T        = numpy.sqrt(Fxreq*Fxreq + Fzreq*Fzreq)
            alpha    = numpy.arctan2(Fxreq,Fzreq)         # in radians
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
            alpha_arr   = numpy.linspace(-5,90,nal)         # -5 to 90, steps of 1 deg
            dal_flap    = numpy.asarray([0])
            nflap       = len(dal_flap)
            res         = numpy.zeros_like(alpha_arr)           # -5 to 90, steps of 1 deg
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
                  al_new   = math.arctan2(FxR, FzR)
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

               print('rotor angle wrt vertical is ',al_new)
               print(ial,res[ial])
#               temp = raw_inumpyut('yoyo')

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

                     print('V = ',vinf/kts2fps,'found solution at alpha = ', alpha)
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
               print('WARNING: DID NOT FIND A SOLUTION: ARE YOU KIDDING ME')
               print('V = ',vinf/kts2fps)

#====================================================================
# iterative solution for aoa
#====================================================================

         else: #option == 'iterate':
            err         = 200.0
            soln_found  = False 
            alpha       = alpha_guess
            while not soln_found:
               Lw       = 0.0; Dw = 0.0
               for gname,group in self.design.wing.groups.items():
                  ARe      = group.aspectratio*group.oswald
                  CLa      = 0.11*ARe/(2.0+ARe)    # lift curve slope / deg
                  alpha_m  = group.cl/CLa                                       # wing mount angle, deg
                  alpha_w0 = 90.0 - alpha + alpha_m
#loop over flap settings at this wing tilt angle
                  Lwg, Dwg = self.wing_loads(alpha_w, vinf, rho, group)
                  Lw       = Lwg + Lw
                  Dw       = Dwg + Dw
                  # print(alpha,alpha_w0,Lwg,Fzreq)
                  # print(alpha_w0,Lwg,self.gtow)
# find force along X, Z reqd from rotor
               FxR      = Fxreq + Dw
               FzR      = Fzreq - Lw
               # x1  = input('yes?')
# update rotor TPP orientation angle
               al_new      = numpy.arctan2(FxR, FzR)*r2d
               alpha_prev  = alpha
               alpha    = 0.95*al_new + 0.05*alpha 

               if abs(alpha - alpha_guess) > dalpha_max:
                  print(alpha_guess,dalpha_max,alpha)
                  max_allow = alpha_guess + numpy.sign(alpha-alpha_guess)*dalpha_max
                  al_print  = alpha
                  alpha     = max_allow
                  # print(alpha_guess,'capping AoA variation: wanted',al_print,' allowed: ',max_allow)
               if abs(alpha_prev - alpha) < 1e-4:
                  soln_found = True 

            T           = numpy.sqrt(FxR*FxR + FzR*FzR)
            P           = self.rotor_power(T, vinf, rho, alpha*d2r)

# print FxR, FzR,al_new,err
# print vinf,P, alpha
# x = raw_inumpyut('ok?')

      alpha    = alpha*d2r
      LbyD     = weight*vinf/P*0.001

      P        = P/self.eta_xmsn(vinf,0)          # get power draw at source
      return P, alpha, soln_found, LbyD

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
         lam2  = mutanal + CT*0.5/numpy.sqrt(musq + lam*lam) 
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

      rotor       = self.design.rotor

      V1          = self.Vthresh
      V2          = self.Vthresh2

      for i in range(rotor.ngroups):
         group       = rotor.groups[0]
         A           = group.area
         nrotor      = group.nrotors
         rpm_ratio   = group.RPM_ratio

#      if vinf < V1:
#         factor   = 1.0 
#      elif vinf <= V2:
#         factor   = 1.0 - (1.0-rpm_ratio)*(vinf-V1)/(V2-V1)
#      else:
         factor      = rpm_ratio 
#      print 'airspeed is ',vinf/kts2fps 
#      yemp = raw_inumpyut('yo?') 
#      print 'rotor rpm ratio is ',factor 
         converged   = False 
         sigma       = group.solidity
         while (not converged):
            vtip           = group.tipspeed*factor   # Vtip at this vinf
            Tnd            = rho*A*vtip*vtip                  # normalized thrust
            CT             = total_thrust/(nrotor*Tnd)        # thrust coeff.
            CTsig          = CT/sigma 
            if(CTsig > 0.14):
               factor      = numpy.sqrt(CTsig/0.14)*factor*1.01 
            else:
               converged   = True
         Pnd         = rho*A*vtip*vtip*vtip/1000.0      # to convert Cp to kW
         mu          = vinf*numpy.cos(alpha)/vtip       # advance ratio

#converge inflow (total inflow)
         mutanal     = vinf*numpy.sin(alpha)/vtip 
         lam         = self.converge_inflow(mu, mutanal, CT) 
         ald         = alpha*r2d
         cd0         = group.cd0#*(90.0-ald)/80.0
         ripf        = group.ipf
         cipf        = 1.0/self.design.prop.eta
         if ald < 10.0:
            ipf      = ripf
         elif ald > 85:
            ipf      = cipf       # no induced losses
         else:                   # between 10 and 100
            ipf      = ripf + (cipf-ripf)*(85.0-ald)/75.0
#calculate rotor induced and profile power 
         cpi         = ipf*CT*lam*nrotor
         Pi          = cpi*Pnd
         cp0         = 0.125*sigma*cd0*(1.0+4.65*mu*mu)
         Po          = nrotor*cp0*Pnd 
         Ptotal      = Pi + Po
      # print total_thrust/2/A
      # print 'FM = ',(CT*sqrt(0.5*CT))/(lam*CT + cp0)
      # print ('%8.3f %8.3f %8.5f %8.5f %8.2f' % (vinf/kts2fps, rho, CT*sqrt(0.5*CT)*ipf*Pnd*self.nrotor, Pi, Po))
      # print(alpha*r2d,total_thrust, vinf, Ptotal)
      
      return Ptotal

# ===================================================================
# transmission efficiency calculation
# ===================================================================
   
   def eta_xmsn(self, vinf, imotor):

      V1       = self.Vthresh
      V2       = self.Vthresh2

      eta_h    = self.design.motor.groups[imotor].hover_efficiency
      eta_c    = self.design.motor.groups[imotor].cruise_efficiency
      if vinf < V1:
         eta   = eta_h 
      elif vinf <= V2:
         eta   = eta_h + (eta_c-eta_h)/(V2-V1)*(vinf-V1)
      else:
         eta   = eta_c 

      return eta 
# ===================================================================
# power curve
# ===================================================================
   
   def power_curve(self):
      
#====================================================================
# forward flight speed-range
#====================================================================

      vmin              = 0.0
      vmax              = round(1.2*self.Vcruise,-1)
      nvinf             = 50

      vinf_range        = numpy.linspace(vmin,vmax,nvinf)

      gtow              = self.gtow
      nflt              = 2
      mf, mp, me        = self.get_wts()        # fuel, payload and empty mass
      weights           = [me+mp+mf, me+mf]     # laden, unladen masses

#====================================================================
# allocate memory for power components (energy method)
#====================================================================

      Power             = numpy.zeros((nvinf, nflt))
      alpha_s           = numpy.zeros((nvinf, nflt))
      LbyD              = numpy.zeros((nvinf, nflt))
      xaxis             = numpy.zeros((nvinf, nflt))

#====================================================================
# loop over loading conditions
#====================================================================

      for iw in range(nflt):
         rho            = self.rho
         weight         = weights[iw]*grav
         vmin           = 0.0 
         xaxis[:,iw]    = numpy.linspace(vmin, vmax, nvinf)

         for j in range(nvinf):
            vinf        = xaxis[j,iw]

#====================================================================
# find trim solution for a given flight condition
#====================================================================

            if j == 0:
               alpha_guess    = 0.0 
            else:
               alpha_guess    = alpha_s[j-1,iw]

            dalpha            = 7.5#*numpy.pi/180.0
            P,als,valid,LtoD  = self.trim_solve(vinf, weight, self.rho, alpha_guess)

#====================================================================
# If trim didnt converge, use previous point
#====================================================================

            if not valid:
               P              = Power[j-1,iw]
               als            = alpha_s[j-1,iw]
               xaxis[j,iw]    = xaxis[j-1,iw]
            Power[j,iw]       = P
            LbyD[j,iw]        = LtoD
            alpha_s[j,iw]     = als*r2d 

#====================================================================
# find cruise point, mark it with red circle
#====================================================================

         Pc,alc,valid,LtoD    = self.trim_solve(self.Vcruise, weight, rho, alpha_guess)

#====================================================================
# update plots
#====================================================================

         plt.figure(2)
         if(iw == 0):
            label             = 'with payload'
            color             = 'b'
            style             = '-'
            alpha             = 1.0
         else:
            label             = 'no payload'
            color             = 'r'
            style             = '--'
            alpha             = 0.4
         plt.plot(xaxis[:,iw],Power[:,iw],linestyle=style,linewidth=2.4,label=label,color=color,alpha=alpha)
         if iw == 0:
            plt.plot(self.Vcruise,Pc,'ro',markersize=12,markerfacecolor='None',markeredgewidth=2.2)
   
         plt.figure(3)
         plt.plot(xaxis[:,iw],alpha_s[:,iw],linestyle=style,linewidth=2.4,label=label,color=color,alpha=alpha)
         if iw == 0:
            plt.plot(self.Vcruise,als*r2d,'ro',markersize=12,markerfacecolor='None',markeredgewidth=2.2)

         plt.figure(4)
         plt.plot(xaxis[:,iw], LbyD[:,iw],linestyle=style,linewidth=2.4,label=label,color=color,alpha=alpha)
         if iw == 0:
            plt.plot(self.Vcruise,LtoD,'ro',markerfacecolor='None',markersize=12,markeredgewidth=2.2)

#====================================================================
# total power
#====================================================================

      plt.figure(2)
      if(self.design.etype == 'electric_motor'):
         lab   = '6 C'
      else:
         lab   = 'installed'

      plt.plot([vmin,vmax],[self.pwr_installed,self.pwr_installed],'k--',linewidth=3.0,label=lab)
      plt.xlabel('\\textbf{Cruise speed (m/s)}')
      plt.ylabel('\\textbf{Battery power required (kW)}')
      plt.grid(linestyle=':')
      xmax     = round(vmax,-1)
      plt.xlim([0.0,xmax])

      plt.ylim([0.0,1.2*self.pwr_installed])
      plt.legend(loc='best')
      plt.tight_layout()
      # plt.savefig('power_curve.png',format='png',dpi=150)
      self.pdf.savefig()
      plt.close()

      plt.figure(3)
      plt.xlabel('\\textbf{Cruise speed (m/s)}')
      plt.ylabel('\\textbf{Wing tilt (deg)}')
      plt.grid(linestyle=':')
      # plt.xlim([0.0,xmax])
      plt.ylim([0.0,100])
      plt.legend(loc='best')
      plt.tight_layout()
      # plt.savefig('alpha.png',format='png',dpi=150)
      self.pdf.savefig()
      plt.close()
#====================================================================
# L/D plot: W*V/P
#====================================================================

      plt.figure(4)
#      for iw in xrange(nflt):
#         idx            = indices[iw]
#         rho            = self.all_rho[idx]       # density in slug/cu.ft
#         weight         = weights[idx]
      plt.xlabel('\\textbf{Cruise speed (m/s)}')
      plt.ylabel('\\textbf{Lift to Drag ratio WV/P$_\\textrm{rotor}$}')
      plt.grid(linestyle=':')
      plt.legend(loc='best')

      # print xmax, 
      # plt.xlim([0.0,xmax])
      plt.ylim([0.0,numpy.ceil(numpy.amax(LbyD*1.1))])
      plt.tight_layout()
      # plt.savefig('LbyD.png',format='png',dpi=150)
      self.pdf.savefig()
      plt.close()
# ===================================================================
# return fuel, payload and empty mass
# ===================================================================

   def get_wts(self):

      design         = self.design
      mfuel          = design.mass_battery
      mpay           = design.mission.payload
      empty_wt       = design.massempty

      return mfuel, mpay, empty_wt

# ===================================================================
# return energy capacity, SOH and DOD
# ===================================================================

   def get_energy(self):

      design                  = self.design 
      Pack                    = design.engine.Pack
      Ebatt                   = design.engine.Eins        # kW-hr
      SOH                     = Pack.SOH
      DOD                     = Pack.DOD_min

      return Ebatt, SOH, DOD

# ===================================================================
# plot hover endurance  with payload
# ===================================================================

   def hover_endurance(self):

      mfuel, mpay, empty_wt   = self.get_wts()
      Ebatt, SOH , DOD        = self.get_energy()

      max_fuel       = mfuel*1.35
      max_payload    = mpay*1.35

#      print self.gtow,empty_wt,mpay,mfuel
      N              = 40
      fuel           = numpy.linspace(0,mpay+mfuel,N)       # fuel: from zero to full supply
      payload        = numpy.zeros_like(fuel)               # payload, kg
      thover         = numpy.zeros((N,2))                   # hover  time

      tcruise        = self.tcruise
      Pc,alc,valid,l = self.trim_solve(self.Vcruise, self.gtow, self.rho, 90.0)
      Ph,alh,valid,l = self.trim_solve(0.0         , self.gtow, self.rho,  0.0)
      Ecruise        = Pc*tcruise/60.0
      Eactual        = Ebatt*(SOH-DOD) - Ecruise 
      tref           = Eactual/Ph*60

#====================================================================
# calculate payload assuming payload + fuel = constant, ie GTOW = constant
# mass fuel    = kgs
# power        = kw
# sfc          = ("mass" fuel)/(power*time)
#
# so time      = fuel / (power * sfc)
#
#====================================================================

      for i in range(N):
         payload[i]    = mpay + mfuel - fuel[i]

#====================================================================
# limit payload and fuel to 135% of design values
#====================================================================

         if fuel[i] > max_fuel:
            fuel[i] = max_fuel

         if payload[i] > max_payload:
            payload[i] = max_payload

#====================================================================
# calculate thrust, power and endurance in hover and cruise
#====================================================================
         
         weight         = (empty_wt + payload[i] + fuel[i])*grav
         Ph,alh,valid,l = self.trim_solve(0.0,          weight, self.rho,  0.0)
         Pc,alc,valid,l = self.trim_solve(self.Vcruise, weight, self.rho, 90.0)
         Ecruise        = Pc*tcruise/60.0
         phover         = Ph

         Eactual        = Ebatt*(SOH-DOD) - Ecruise 
         Eideal         = Ebatt           - Ecruise
         Eibatt         = Eideal *fuel[i]/mfuel
         Eabatt         = Eactual*fuel[i]/mfuel

         # print(payload[i],fuel[i],thrust)
         if (phover < self.pwr_installed and Eibatt >0 and Eabatt > 0):

            thover[i,0] = Eibatt/phover * 60.0 # convert from hrs to mins
            thover[i,1] = Eabatt/phover * 60.0 # convert from hrs to mins

# ===================================================================
# make a plot
# ===================================================================
         
      plt.figure()
      plt.plot(thover[:,0], payload,'-r',label='Ideal battery',linewidth=2.4,alpha=0.4)
      plt.plot(thover[:,1], payload,'-o',label='Vahana assumptions',linewidth=2.4,markevery=4)
      plt.plot(tref, mpay,'o',markersize=16,markerfacecolor='none',markeredgewidth=2.5,markeredgecolor='r')
      plt.ylabel('\\textbf{Payload (kg)}')
      plt.xlabel('\\textbf{Time in hover (min)}')
      plt.xlim([0,numpy.max(thover[:,0]*1.1)])
      plt.ylim([0,round(numpy.max(payload)*1.1,-2)])
      plt.grid(linestyle=':')
      plt.title('\\textbf{Payload-Endurance in hover, fixed cruise distance}')
      plt.legend(loc='best')
      plt.tight_layout()
      # plt.savefig('hover_endurance.png',format='png',dpi=150)
      self.pdf.savefig()
      plt.close()
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

      mfuel, mpay, empty_wt   = self.get_wts()
      Ebatt, SOH , DOD        = self.get_energy()

      max_fuel       = mfuel*1.35
      max_payload    = mpay*1.35

      N              = 40

#====================================================================
# available energy for cruise 
#====================================================================

      fuel           = numpy.linspace(0,mpay+mfuel,N)       # fuel: from zero to full supply
      payload        = numpy.zeros_like(fuel)               # payload, kg

      Pc,alc,valid,l = self.trim_solve(self.Vcruise, self.gtow, self.rho, 90.0)
      Ph,alh,valid,l = self.trim_solve(0.0         , self.gtow, self.rho,  0.0)

      thover         = self.thover 
      Ehover         = Ph*thover/60.0
      Eactual        = Ebatt*(SOH-DOD) - Ehover
      tref           = Eactual/Pc*60.0
      dV             = 7.5
      all_V          = [self.Vcruise, self.Vcruise-dV, self.Vcruise+dV]
      nV             = len(all_V)
      tcruise        = numpy.zeros((N,2,nV))                   # cruise time
      # print(tref);quit()
#====================================================================
# calculate payload assuming payload + fuel = constant, ie GTOW = constant
# mass fuel    = kgs
# power        = kw
# sfc          = ("mass" fuel)/(power*time)
#
# so time      = fuel / (power * sfc)
#
#====================================================================

   
      for i in range(N):
         payload[i]    = mpay + mfuel - fuel[i]

#====================================================================
# limit payload and fuel to 135% of design values
#====================================================================

         if fuel[i] > max_fuel:
            fuel[i] = max_fuel

         if payload[i] > max_payload:
            payload[i] = max_payload

#====================================================================
# calculate thrust, power and endurance 
#====================================================================

         weight         = (empty_wt + payload[i] + fuel[i])*grav
         Ph,alh,valid,l = self.trim_solve(0.0,          weight, self.rho,  0.0)
         Ehover         = Ph*thover/60.0
         Eactual        = Ebatt*(SOH-DOD) - Ehover 
         Eideal         = Ebatt           - Ehover
         Eibatt         = Eideal *fuel[i]/mfuel
         Eabatt         = Eactual*fuel[i]/mfuel         

         for j in range(nV):
            V                 = all_V[j]
            Pc,alc,valid,l    = self.trim_solve(V, weight, self.rho, 90.0)

            if (Pc < self.pwr_installed and Eibatt > 0 and Eabatt >0):
               tcruise[i,0,j] = Eibatt/Pc * 60.0 # convert from hrs to mins
               tcruise[i,1,j] = Eabatt/Pc * 60.0 # convert from hrs to mins

#====================================================================
# payload-endurance
#====================================================================
      
      plt.figure()
      for j in range(nV):
         V     = all_V[j]
         plt.plot(tcruise[:,1,j], payload,'-o',label=str(int(round(V,0))) + ' m/s',linewidth=2.0)

      for j in range(nV):
         V     = all_V[j]
         if(V == self.Vcruise):

            plt.plot(tcruise[:,0,j], payload,'-r',label=str(int(round(self.Vcruise,0))) + ' m/s, ideal battery',linewidth=2.4,alpha=0.4)
   
      plt.plot(tref, mpay,'o',markeredgecolor='r',markersize=16,markerfacecolor='none',markeredgewidth=2.5)
      plt.ylabel('\\textbf{Payload (kg)}')
      plt.xlabel('\\textbf{Time in cruise (min)}')
      plt.xlim([0,numpy.max(tcruise[:,0,1]*1.1)])
      plt.ylim([0,round(numpy.max(payload)*1.1,-2)])
      plt.grid(linestyle=':')
      plt.title('\\textbf{Payload-Endurance in Cruise}')
      plt.legend(loc='best')
      plt.tight_layout()
      # plt.savefig('cruise_endurance.png',format='png',dpi=150)
      self.pdf.savefig()
      plt.close()
#====================================================================
# payload-range
#====================================================================

      Range    = numpy.zeros_like(tcruise)
      plt.figure()
      for j in range(nV):
         V     = all_V[j]
         plt.plot(tcruise[:,1,j]*0.06*V, payload,'-o',label=str(int(round(V,0))) + ' m/s',linewidth=2.0)
         Range[:,1,j] = tcruise[:,1,j]*0.06*V
      for j in range(nV):
         V     = all_V[j]
         if(V == self.Vcruise):
            plt.plot(tcruise[:,0,j]*0.06*V, payload,'-r',label=str(int(round(self.Vcruise,0))) + ' m/s, ideal battery',linewidth=2.4,alpha=0.4)
            Range[:,0,j] = tcruise[:,0,j]*0.06*V

      # plt.figure()
      # plt.plot(Range[:,0], payload,'-r',label='Ideal battery assumptions',linewidth=2.4,alpha=0.4)
      # plt.plot(Range[:,1], payload,'-o',label='Vahana: ',linewidth=2.4)
      plt.plot(tref*self.Vcruise*0.06, mpay,'o',markeredgecolor='r',markersize=16,markerfacecolor='none',markeredgewidth=2.5)
      plt.ylabel('\\textbf{Payload (kg)}')
      plt.xlabel('\\textbf{Cruise range (km)}')
      plt.xlim([0,numpy.max(Range[:,0]*1.1)])
      plt.ylim([0,round(numpy.max(payload)*1.1,-2)])
      plt.grid(linestyle=':')
      plt.title('\\textbf{Payload-Range [3 min hover]}')
      plt.legend(loc='best')
      plt.tight_layout()
      # plt.savefig('cruise_range.png',format='png',dpi=150)
      self.pdf.savefig()
      plt.close()
#====================================================================
# driver function to launch plots
#====================================================================

   def make_plots(self, design_id):

      fname          = 'performance_design_' + str(design_id) + '.pdf'
      with PdfPages(fname) as pdf:
         self.pdf    = pdf 
         self.power_curve()
         self.hover_endurance()
         self.cruise_performance()
      print('saved results to file ',fname )