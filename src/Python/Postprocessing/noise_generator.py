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
import matplotlib        as mpl

#=====================================
# use this line to save images to PDF
#=====================================

from matplotlib.backends.backend_pdf import PdfPages

#====================================================================
# unit conversions
#====================================================================

# kts2fps = 1.68781
# fps2kts = 1.0/kts2fps
# fps2nm  = 0.000164579
# rho     = 0.002378         # slugs/ft3
# grav    = 32.18            # ft/s/s

#si units
grav    = 9.81             # m/s/s
# rho     = 1.2256           # kg/cu.m
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
# extra paths
#====================================================================

import os, numpy, math, pickle, matplotlib, sys, time

sys.path.append('../Stage_0')
sys.path.append('../Stage_1')
sys.path.append('../Stage_1/afdd')
sys.path.append('../Stage_1/engines')
sys.path.append('../Stage_3')

from hydraInterface     import hydraInterface

#====================================================================
# Class definition
#====================================================================

def noise(filename):

   footprint      = {}
   
#====================================================================
# initialization: load data
#====================================================================
   
   if not os.path.isfile(filename):
      print('\n\nFile:',filename,'does not exist.\n\n')
      sys.exit(1)

   with open(filename,'rb') as f:
      design    = pickle.load(f)

#====================================================================
# get ratio of hover to cruise RPM
#====================================================================

   rotor        = design.rotor
   R            = rotor.radius
   Omega_ratio  = design.all_dict['aircraft']['cruise_rpm_ratio']

#====================================================================
# Generate rotor layout
#====================================================================

   l            = 5.25            # longitudinal spacing of wings 
   h            = 1.5*R           #     vertical spacing of wings
   altitude     = 250.0*0.3048    # altitude meters

#====================================================================
# loop over wings, get rotors per wing
#====================================================================

   nr           = []
   geom         = design.emp_data.Geometry
   bf           = geom.fuselage_width
   clearance    = geom.clearance*R     # clearance between rotor tip and another object
   x            = []        # location of rotor in body X axis, m
   y            = []        # location of rotor in body Y axis, m
   z            = []        # location of rotor in body Z axis, m

   wings        = design.wing.groups
   ngroups      = design.wing.ngroups

   l_shaft      = R

   for key,group in wings.items():
      for i in range(group.nwings):
          nr.append(int(group.nrotors))
          l_shaft   = min(l_shaft, group.chord)
          # print(group.span)
   # print('shaft length',l_shaft)
#====================================================================
# loop over rotor groups
#====================================================================
   
   for ir,nr_wing in enumerate(nr):

#====================================================================
# place first set at x=0, z=0
#====================================================================

      if(ir == 0):
         xr    = 0.5*l; zr    = 0.5*h

#====================================================================
# place second set at x=-l, z=-h (back, up above CG)
#====================================================================

      else:
         xr    =-0.5*l  ; zr =-0.5*h

#====================================================================
# place rotors spanwise along the wing
#====================================================================

      for j in range(int(nr_wing/2)):

#====================================================================
# start inside and go out to wing tip;
# these are motor mount locations!
#====================================================================

        yr    = bf*0.5+(clearance+2*R)*(j)+clearance+R

        x.append(xr); y.append(yr) ; z.append(zr)
        x.append(xr); y.append(-yr); z.append(zr)

   print('   Hub X,   Hub Y,  Hub Z')
   print('=========================')
   for xr,yr,zr in zip(x,y,z):
      print('  {:6.3f}   {:6.3f}  {:6.3f} '.format(xr,yr,zr))

#====================================================================
# calculate shaft tilt in hover vs. cruise
#====================================================================

   alpha_s_h   = numpy.ones(len(x))*numpy.pi*0.5    # 90 deg pitch up for hover orientation
   alpha_s_c   = numpy.zeros(len(x))

#====================================================================
# get distance and angle at observer from rotor axis for different points 
# generate noise footprint for hover
#====================================================================

   Sh, thetah,obs,hubs = observer_posn(x,y,z,alpha_s_h,l_shaft,altitude,R)
   noise_h             = noise_footprint(Sh, thetah, design.rotor, obs)

#====================================================================
# draw rotor discs in space and noise levels too
#====================================================================

   draw_plot(R, hubs, alpha_s_h, noise_h, obs)

#====================================================================
# cruise condition: noise prediction
#====================================================================

#   Sc, thetac  = observer_posn(x,y,z,alpha_s_c,l_shaft,altitude,R)
#   noise_c     = noise_footprint(Sc, thetac)

#====================================================================
# return noise footprint
#====================================================================

   return footprint

#====================================================================
# function to calculate observer position and rotor positions wrt CG
# of vehicle; observer positions are computed on a rectangular grid at
# points at a given altitude in meters
#
# xr,yr,zr: locations of the wing points where rotors are mounted, body axes, meters
# alpha_s : shaft tilt angles of the rotors, radians, zero = airplane mode, positive nose up
# l_shaft : shaft length, meters; distance from rotor tilt axis and rotor hub
#
# returns:
#       S: distance from rotor hub to observer, meters
#   theta: angle between rotor Z-axis (direction of thrust) and line joining
#          rotor hub to observer
#     obs: observer locations on the ground, (3xnobsxnobs) array
#    hubs: X,Y,Z (body axes) location of rotor hubs on vehicle
#====================================================================

def observer_posn(xr,yr,zr,alpha_s,l_shaft,altitude,radius):

  xMax    = 150           # meters away from CG, measured in horizontal plane along X,Y axes
  dX      = 10            # grid size on ground, meters
  n       = 2*int((xMax/dX)+1)
  line    = numpy.linspace(-xMax,xMax,n)
  obs     = numpy.zeros((3,n,n))
  z       = altitude*numpy.ones((n,n))

  nr      = len(xr)

#====================================================================
# generate xy grid on ground
#====================================================================

  for i in range(n):
    for j in range(n):
      obs[0,i,j]  = line[i]
      obs[1,i,j]  = line[j]
      obs[2,i,j]  = altitude 

#====================================================================
# now generate positions of the rotor hubs in space based on tilt 
# angles alpha_s
#====================================================================

  hubs   = numpy.zeros((3,nr))
  vecs   = numpy.zeros((3,nr))
  for ir in range(nr):
     hubs[1,ir]   = yr[ir]       # y position is unaltered
     hubs[0,ir]   = xr[ir] + l_shaft*numpy.cos(alpha_s[ir])
     hubs[2,ir]   = zr[ir] - l_shaft*numpy.sin(alpha_s[ir]) 

#====================================================================
# generate normal vectors: remember that this normal is in body frame!
#====================================================================

     vecs[0,ir]   = -numpy.cos(alpha_s[ir])
     vecs[2,ir]   = -numpy.sin(alpha_s[ir])

#====================================================================
# for each (X,Y) point on the ground, calculate distance from each 
# rotor hub to that point
#====================================================================

  S       = numpy.zeros((n,n,nr))
  theta   = numpy.zeros_like(S)
  
  for i in range(n):
    for j in range(n):
      point           = obs[:,i,j]

      for ir in range(nr):
        rhub          = hubs[:,ir]
        dX            = point[0] - rhub[0]       # vector X component from hub to observer
        dY            = point[1] - rhub[1]       # vector Y component from hub to observer
        dZ            = point[2] - rhub[2]       # vector Z component from hub to observer
        S[i,j,ir]     = numpy.sqrt(dX*dX + dY*dY + dZ*dZ)
        dX            = dX/S[i,j,ir]
        dY            = dY/S[i,j,ir]
        dZ            = dZ/S[i,j,ir]
        theta[i,j,ir] = numpy.arccos(vecs[0,ir]*dX + vecs[1,ir]*dY + vecs[2,ir]*dZ)
        # print(dX,dY,dZ,vecs[:,ir])
        # print(theta[i,j,ir])
        # x = input('ok theta?')
#====================================================================
# end of operations
#====================================================================

  return S, theta, obs, hubs

#====================================================================
# function to draw rotor discs in space
# radius    = rotor radius, meters
# hubs      = x,y,z locations of hubs in space, body axes, size (3 x nrotors)
# alpha_s   = shaft tilt angle wrt horizontal, radians, +ve pitch up
# noise     = SPL array from all rotors, (nobs x nobs) array
# obs       = observer locations, (3 x nobs x nobs) array 
#====================================================================

def draw_plot(radius, hubs, alpha_s, noise, obs):
  fig   = plt.figure()
  ncir         = 41
  theta        = numpy.linspace(0.0,2*numpy.pi,ncir)
  costh        = numpy.cos(theta)
  sinth        = numpy.sin(theta)
  zcir         = radius*costh
  ycir         = radius*sinth
  yc2          = ycir
  nr           = numpy.shape(hubs)[1]
  plot_3d      = False 

  fname        = 'noise_profile_hover.pdf'
  with PdfPages(fname) as pdf:

    if(plot_3d):
      ax    = fig.add_subplot(111, projection='3d')
      ax.scatter(-hubs[0,:],hubs[1,:],-hubs[2,:],'bs')
    else:
      plt.plot(-hubs[0,:],hubs[1,:],'bs')

# draw rotor discs
    for ir in range(nr):
      xc2          = zcir*numpy.sin(alpha_s[ir])
      zc2          = zcir*numpy.cos(alpha_s[ir])

# 3d: all coords, 2d only top view
      if(plot_3d):
        ax.plot(xc2-hubs[0,ir],yc2+hubs[1,ir],zc2-hubs[2,ir],'-k')
      else:
        plt.plot(xc2-hubs[0,ir],yc2+hubs[1,ir],'-k')

#====================================
# draw bounding box to get equal axes 
#====================================

    if(plot_3d):
      X         = -hubs[0,:]      # from FD body axes to plot axes
      Y         =  hubs[1,:]      # obtained with 180 deg rotation about Y  
      Z         = -hubs[2,:]      #
      delta     = radius 
      max_range = numpy.array([X.max()-X.min()+radius, Y.max()-Y.min()+radius, Z.max()-Z.min()+radius]).max()
      Xb = 0.5*max_range*numpy.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
      Yb = 0.5*max_range*numpy.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
      Zb = 0.5*max_range*numpy.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
      for xb, yb, zb in zip(Xb, Yb, Zb):
       ax.plot([xb], [yb], [zb], 'w')

#====================================
# draw contour map
#====================================

    else:
      plt.axis('equal')
      # print(numpy.shape(noise))
      plt.contourf(obs[0,:,:],obs[1,:,:],noise,20,cmap='jet')
      norm  = mpl.colors.Normalize(vmin=45,vmax=75)
      plt.colorbar(norm=norm)
      maxN  = numpy.max(noise)
      flatn = numpy.ndarray.flatten(noise)
      imaxN = numpy.where(flatn >= maxN-0.5)

#====================================
# find locations where noise is within 2dB of this maximum
#====================================

      xmax  = numpy.ndarray.flatten(obs[0,:,:])[imaxN]
      ymax  = numpy.ndarray.flatten(obs[1,:,:])[imaxN]
      r     = numpy.sqrt(numpy.square(xmax) + numpy.square(ymax))
      rmin  = numpy.min(r)
      rmax  = numpy.max(r)

      xmin  = rmin*costh; xmax  = rmax*costh
      ymin  = rmin*sinth; ymax  = rmax*sinth

#====================================
# plot
#====================================

      plt.plot(xmax,ymax,'--k',markersize=4,markerfacecolor='None')
      plt.plot(xmin,ymin,'--k',markersize=4,markerfacecolor='None')
      plt.xlabel('Distance (m)')
      plt.ylabel('Distance (m)')
      title = 'Hover altitude = ' + str(int(obs[2,0,0]/0.3048)) + ' feet: sound pressure level (dB)'
      plt.title(title)
    # plt.clim(45,75)
#====================================
# show the plot
#====================================

    #plt.show()
    pdf.savefig()
  return

#====================================================================
# function to predict noise variation over a square on the ground
# S           = distance from rotor hub to observer locations,        (nobs x nobs x nrotors)
# theta       = angle between thrust direction and vector to observer (nobs x nobs x nrotors)
# rotor       = pointer to rotor class entry of hydra design output
#====================================================================

def noise_footprint(S_all, theta_all, rotor, obs):

  nobs        = numpy.shape(S_all)[0]
  nr          = numpy.shape(S_all)[2]
  dnoise      = numpy.zeros(nr)
  noise       = numpy.zeros((nobs,nobs))

  S           = S_all[0,0,0]
  theta       = theta_all[0,0,0]
  # print(obs[:,0,0],S,theta)
  # print(rotor.fanNoise(S, theta))
  # quit()
  # print(numpy.shape(S_all)),nr
  for i in range(nobs):
    for j in range(nobs):
      for ir in range(nr):
        S           =     S_all[i,j,ir]
        theta       = theta_all[i,j,ir]
        dnoise[ir]  = rotor.fanNoise(S,theta)
        noise[i,j]  = noise[i,j] + (10**(dnoise[ir]*0.1))
        # print(i,j,S,theta,dnoise[ir])
        # x1=input('ok?')
      noise[i,j]    = 10*numpy.log10(noise[i,j])
  return noise 
