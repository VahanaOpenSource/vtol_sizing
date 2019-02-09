#====================================================================
#
# file to extract performance
#
# power curve
# payload range diagram
# payload endurance diagram
#====================================================================

import matplotlib, numpy
from wing_carpets_v2 import wing_carpets as wc 
from line_sensitivities import line_sensitivities as ls
import copy
#====================================================================
# unit conversions
#====================================================================

grav    = 9.81             # m/s/s

#==========================
# use these lines for latex
#==========================

from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

rc('text',usetex=True)
font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 13}

rc('font', **font)

#====================================================================
# Class definition
#====================================================================

class dv_sensitivities:

#====================================================================
# initialization: load data
#====================================================================
   
   def __init__(self, design, design_id):
      print('ok')
      self.design_id  = design_id
      self.design     = design 
      self.fname      = 'sensitivities_design_' + str(design_id) + '.pdf'

#====================================================================
# sensitivities for wing design
#====================================================================

   def rotor_sensitivities(self):

     design           = self.design 
     rotor            = design.rotor 
     vtip_min         = 130.0;  vtip_max  = 180.0; dvtip  = 10.0 
     nvtip            = int((vtip_max-vtip_min)/dvtip)+1
     vtip_all         = numpy.linspace(vtip_min,vtip_max,nvtip)

     BL_min           = 0.08;  BL_max  = 0.135; dBL  = 0.005
     nBL              = int((BL_max-BL_min)/dBL)+1
     BL_all           = numpy.linspace(BL_min,BL_max,nBL)

#====================================================================
# converge for baseline design, get array of quantities to plot
#====================================================================

     baseline         = design.get_essentials()

     nqtys            = len(baseline)
     nextra           = 2
     adict            = design.all_dict['aircraft']
     baseline_adict   = copy.deepcopy(adict)

#====================================================================
# loop over rotor groups
#====================================================================

     for i in range(rotor.ngroups):
        group         = rotor.groups[i]
        all_data      = numpy.zeros((nvtip,nBL,nqtys+nextra))
        inp_group     = adict['rotor']['group'+str(i)]
        xref          = group.diskloading/9.81*2.2*0.3048**2
        yref          = group.solidity
        # print(xref,yref,baseline)
        for iBL,BL in enumerate(BL_all):
           for ivtip,vtip in enumerate(vtip_all):

              inp_group['ctsigma']  = BL 
              inp_group['Vtip']     = vtip

#====================================================================
# converge to design, save data to matrix
#====================================================================

              all_data[ivtip,iBL,:nqtys]    = design.get_essentials()

              # print(vtip,BL,group.solidity,group.diskloading)
#====================================================================
# calculate span of wing group being perturbed, as well as wing loading
#====================================================================
              
              all_data[ivtip,iBL,nqtys]   = group.diskloading
              all_data[ivtip,iBL,nqtys+1] = group.solidity

#====================================================================
# draw carpet
#====================================================================

        design.all_dict['aircraft']     = baseline_adict 

        Mass        = all_data[:,:,1]
        Power       = all_data[:,:,2]
        Fuel        = all_data[:,:,3]
        Battery     = all_data[:,:,4]
        Cost        = all_data[:,:,6]
        valid       = all_data[:,:,0]
        DL          = all_data[:,:,nqtys]/9.81*2.2*0.3048**2
        Sigma       = all_data[:,:,nqtys+1]
        cstr        = 'rotor (group ' + str(i) + ') design variables'

        ystr        = 'Solidity'
        xstr        = 'Disk loading (lb/sq.ft)'
        astr        = ' m/s'
        bstr        = 'C$_T$/$\\sigma$='

        wc(vtip_all, BL_all, DL, Sigma, Mass,    xref,yref, baseline[1], valid, 'Mass (kg)'   , cstr,xstr,ystr,astr,bstr)
        self.pdf.savefig(); plt.close()
        wc(vtip_all, BL_all, DL, Sigma, Power,   xref,yref, baseline[2], valid, 'Power (kW)'  , cstr,xstr,ystr,astr,bstr)
        self.pdf.savefig(); plt.close()
        wc(vtip_all, BL_all, DL, Sigma, Battery, xref,yref, baseline[4], valid, 'Battery (kg)', cstr,xstr,ystr,astr,bstr)
        self.pdf.savefig(); plt.close()
        wc(vtip_all, BL_all, DL, Sigma, Cost,    xref,yref, baseline[6], valid, '(USD/hr)'    , cstr,xstr,ystr,astr,bstr)
        self.pdf.savefig(); plt.close()

#====================================================================
# sensitivities for wing design
#====================================================================

   def wing_sensitivities(self):

     design           = self.design 
     wing             = design.wing 
     AR_min           = 4.0;  AR_max  = 12.0; dAR  = 1.0 
     nAR              = int((AR_max-AR_min)/dAR)+1
     AR_all           = numpy.linspace(AR_min,AR_max,nAR)

     CL_min           = 0.2;  CL_max  = 0.8; dCL  = 0.1
     nCL              = int((CL_max-CL_min)/dCL)+1
     CL_all           = numpy.linspace(CL_min,CL_max,nCL)

#====================================================================
# converge for baseline design, get array of quantities to plot
#====================================================================

     baseline         = design.get_essentials()

     nqtys            = len(baseline)
     nextra           = 2
     adict            = design.all_dict['aircraft']
     baseline_adict   = copy.deepcopy(adict)

#====================================================================
# loop over wing groups
#====================================================================

     for i in range(wing.ngroups):
        group         = wing.groups[i]
        all_data      = numpy.zeros((nAR,nCL,nqtys+nextra))
        inp_group     = adict['wing']['group'+str(i)]
        xref          = group.span 
        yref          = baseline[1]*grav*group.lift_frac/(group.area*group.nwings)
        for iCL,CL in enumerate(CL_all):
           for iAR,AR in enumerate(AR_all):

              inp_group['aspectratio']  = AR 
              inp_group['cl']           = CL 

#====================================================================
# converge to design, save data to matrix
#====================================================================

              all_data[iAR,iCL,:nqtys]    = design.get_essentials()

#====================================================================
# calculate span of wing group being perturbed, as well as wing loading
#====================================================================
              
              all_data[iAR,iCL,nqtys]     = group.span
              mass                        = all_data[iAR,iCL,1]
              all_data[iAR,iCL,nqtys+1]   = mass*grav*group.lift_frac/(group.nwings*group.area)      # N/sq.m

#====================================================================
# draw carpet
#====================================================================

        design.all_dict['aircraft']     = baseline_adict 

        Mass        = all_data[:,:,1]
        Power       = all_data[:,:,2]
        Fuel        = all_data[:,:,3]
        Battery     = all_data[:,:,4]
        Cost        = all_data[:,:,6]
        valid       = all_data[:,:,0]
        Span        = all_data[:,:,nqtys]
        WL          = all_data[:,:,nqtys+1]
        prefix      = 'wing (group ' + str(i) 

        ystr        = 'Wing loading (N/m$^2$)'
        xstr        = 'Wing span (m)'
        cstr        = prefix + ') design variables'
        astr        = 'AR='
        bstr        = 'C$_L$='        

        wc(AR_all, CL_all, Span, WL, Mass,    xref,yref, baseline[1], valid, 'Mass (kg)'   , cstr,xstr,ystr,astr,bstr)
        self.pdf.savefig(); plt.close()
        wc(AR_all, CL_all, Span, WL, Power,   xref,yref, baseline[2], valid, 'Power (kW)'  , cstr,xstr,ystr,astr,bstr)
        self.pdf.savefig(); plt.close()
        wc(AR_all, CL_all, Span, WL, Battery, xref,yref, baseline[4], valid, 'Battery (kg)', cstr,xstr,ystr,astr,bstr)
        self.pdf.savefig(); plt.close()
        wc(AR_all, CL_all, Span, WL, Cost,    xref,yref, baseline[6], valid, '(USD/hr)'    , cstr,xstr,ystr,astr,bstr)
        self.pdf.savefig(); plt.close()

#====================================================================
# show sensitivity to cruise airspeed (line plot)
#====================================================================

   def airspeed_sensitivity(self):

      design           = self.design
      mission          = design.mission 
      nseg             = mission.nseg

#====================================================================
# identify cruise segment for sizing
#====================================================================

      seg_ids          = []
      speeds           = []
      for i in range(nseg):
         segment       = mission.segment[i]
         if(segment.type == 'all' and segment.flightmode == 'cruise'):
            seg_ids.append(i)                       # segment id
            speeds.append(segment.cruisespeed)      # in knots

#====================================================================
# converge for baseline design, get array of quantities to plot
#====================================================================

      baseline         = design.get_essentials()
      nqtys            = len(baseline)
      adict            = design.all_dict['aircraft']
      baseline_adict   = copy.deepcopy(adict)
      xref             = speeds[0]

#====================================================================
# setup airspeed array, initialize output data storage
#====================================================================

      nspd             = 11
      ratios           = numpy.linspace(0.75,1.25,nspd)
      all_data         = numpy.zeros((nspd,nqtys))

#====================================================================
# loop over different airspeeds, set segment properties and calculate
# new vehicle after sizing
#====================================================================

      for i in range(nspd):
         for seg_id,speed in zip(seg_ids,speeds):
            mission.segment[seg_id].cruisespeed = speed*ratios[i]

         all_data[i,:] = design.get_essentials()

#         for seg_id,speed in zip(seg_ids,speeds):
#            print(mission.segment[seg_id].cruisespeed,mission.segment[seg_id].time)

#====================================================================
# reset state to baseline: airspeeds
#====================================================================

      for seg_id,speed in zip(seg_ids,speeds):
         mission.segment[seg_id].cruisespeed = speed

#====================================================================
# plot sensitivities
#====================================================================

      Mass        = all_data[:,1]
      Power       = all_data[:,2]
      Fuel        = all_data[:,3]
      Battery     = all_data[:,4]
      Cost        = all_data[:,6]
      valid       = all_data[:,0]

      title       = 'Effect of cruise airspeed'
      xstr        = 'Cruise airspeed (knots)'

#line plot 1: mass and power
      ls(ratios*speeds[0], Mass, Power, xref,baseline[1], baseline[2],valid, xstr, 'Mass (kg)', 'Power (kW)',title)
      self.pdf.savefig(); plt.close()
#line plot 2: battery mass and cost
      ls(ratios*speeds[0], Battery, Cost, xref,baseline[4], baseline[6],valid, xstr, 'Battery (kg)', 'Cost (USD/hr)',title)
      self.pdf.savefig(); plt.close()

#====================================================================
# driver function to launch plots
#====================================================================

   def make_plots(self):

      fname          = self.fname
      with PdfPages(fname) as pdf:
         self.pdf    = pdf 
         self.airspeed_sensitivity()
         self.wing_sensitivities();print('completed wing perturbation plots')
         self.rotor_sensitivities();print('completed rotor DV perturbations, plots')
      print('saved results to file ',fname )