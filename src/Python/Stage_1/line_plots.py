#=======================================================================
# Function to make line plot
#=======================================================================
from   pylab         import *
from   matplotlib    import pyplot
import numpy
import math, numpy, os
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties

#=========================================================================
# function to plot multiple data sets stored in dictionaries and save pngs
#=========================================================================
   
def plot_wrapper(plot_data, keys, labels, save_folder):

   marker_size    = 5

#=========================================================================
# Create folder if required
#=========================================================================
   
   if not os.path.exists(save_folder):
       os.makedirs(save_folder)

#=========================================================================
# loop over all plots to create
#=========================================================================

   ikey           = 0
#   print 'printing keys'
#   for key,data in plot_data.iteritems():
#      print 'data set name is ',key
#      for k2,d2 in data.iteritems():
#         print 'data subset name is ',k2

   lsty_defaults  = ['-','-','--','--']     ; ilin = -1
   lcol_defaults  = ['k','r','b','0.7']   ; icol = -1
   mtyp_defaults  = ['','s','','o']      ; imar = -1

   for key,ylabel in zip(keys,labels):
      ikey        = ikey + 1
      print 'looking for ylabel id ',ylabel

      ilin        = -1; icol = -1; imar = -1
#=========================================================================
# loop over all data sets
#=========================================================================
      
      plot_erase   = True
      for setname,dataset in sorted(plot_data.iteritems()):

#=========================================================================
# extract x axis and plot options
#=========================================================================

#         print setname, dataset
#=========================================================================
# if key exists in dataset, plot it
#=========================================================================
   
         if key in dataset:

            x         = dataset['x']
            xlbl      = dataset['xlbl']
            try:
              lsty      = dataset['lstyle']
            except:
              # print ilin
              ilin      = ilin + 1; lsty = lsty_defaults[ilin]
              if (ilin >= len(lsty_defaults)):
                ilin    = -1

            try:
              lcol      = dataset['lcolor']
            except:
              icol      = icol + 1; lcol = lcol_defaults[icol]
              print setname,lcol_defaults[icol]
              if (icol >= len(lcol_defaults)):
                icol    = -1



            try:
              mtyp      = dataset['mtype']
            except:
              imar      = imar + 1; mtyp = mtyp_defaults[imar]
              if (imar >= len(mtyp_defaults)):
                imar    = -1

            try:
              tag    = dataset['tag']
            except:
              tag    = ''

            try:
              plot_type = dataset['plot_type'].lower()
            except:
              plot_type = 'line'
              # print 'WARNING: PLOT TYPE NOT FOUND; DEFAULTING TO LINE'

#=======================================================================
# Line/Bar graph: call the plot wrapper
#=======================================================================

            if plot_type in ['bar','line']:

#               print 'adding data for ',key
               y           = dataset[key]  
               fid         = ikey      

               a, ax,f     = make_plot(x, y, ikey, lcol, mtyp, lsty,
                     2.5, marker_size, xlbl, ylabel, plot_erase, tag, plot_type)
               plot_erase  = False      # stop erasing figure when adding lines

#=========================================================================
# Contour plot: extract additional information and call the wrapper
#=========================================================================

            elif plot_type == 'contour':
               y           = dataset['y']         # y data must be given like x
               zmin        = dataset['zmin']      # zmin and zmax given like x
               zmax        = dataset['zmax']
               z           = dataset['key']       # z: level data is in the key
               title_str   = key      
               make_contour(x, y, z, fid, title_str, zmin, zmax)

#=======================================================================
# Other plot types: not yet initialized
#=======================================================================

            else:
               print "allowed types are 'line','bar' or 'contour'"
               print 'you may also use lower/upper case letters'
               print 'unfortunately, user specified this:',plot_type
               quit('ERROR: UNKNOWN PLOT TYPE: PROGRAM EXITING')

#=======================================================================
# save the figure and close it
#=======================================================================

      img_name      = save_folder + os.sep + key
      print 'saving figure ',img_name
      try:
        save_png(img_name)
      except:
        save_eps(img_name)
      plt.close(fid)
      
   return None

#=========================================================================
# function to plot multiple data sets 
#=========================================================================

def make_plot(x, y,  fid, line_color, marker_type, line_style,             \
                   line_width,  marker_size, xlbl, ylbl, clear_fig, tag,   \
                   plot_type):

#=======================================================================
# Set default font sizes for plotting
#=======================================================================

   plt.rc('text', usetex=True)

   fs       = 25
   font = {'family' : 'Times',
           'weight' : 'bold',
           'size'   : fs}

   matplotlib.rc('font', **font)

   fontP    = FontProperties()
   fontP.set_family('monospace')

#=======================================================================
# determine marker interval
#=======================================================================
   
   if line_style == '' or line_style == ' ':
      mark_int = 1
   else:
      if len(x) <= 6:
         mark_int = 1
      else:
         mark_int = math.ceil(float(len(x))/6.0)

   mark_int = 1
#   print 'marker interval is ',mark_int
  
#=======================================================================
# draw a figure
#=======================================================================

   fig         = figure(fid)      

   if clear_fig:
      clf()

   ax          = fig.add_subplot(111)

   if plot_type == 'line':
      pyplot.plot(x,y,linewidth=line_width,color = line_color, linestyle = line_style,
                                   marker = marker_type, markersize = marker_size,
                                   markeredgewidth = 1.0,markeredgecolor = line_color,
                                   markevery = mark_int)
   elif plot_type == 'bar':
#      x = scipy.arange(4)
#      y = scipy.array([4,7,6,5])
#      f = figure()
#      ax = f.add_axes([0.1, 0.1, 0.8, 0.8])
      ax.bar(x, y, align='center', color = line_color)
      ax.set_xticks(x)
      ax.set_xticklabels(xlbl,fontsize=fs, fontweight='bold')
#      ax.bar(x, y,color=line_color,align='center')
#      ax.set_xticks(x)
#      ax.set_xticklabels(['1','2','3'])
      #ax.set_xticklabels(['Aye', 'Bee', 'Cee', 'Dee'])
      #ax.xticks(x, xlbl, rotation='vertical')
      print 'hello'
   else:
      print 'UNKNOWN PLOT OPTION: MUST BE line OR bar'
      print 'CRITICAL ERROR : PROGRAM TERMINATING'
      quit()

#=======================================================================
#if plotting azimuth for rotorcraft, use [0,360] xticks
#=======================================================================

   if plot_type == 'line':  
     if(xlbl[0:7].lower() == 'azimuth'):
        if(min(x) >= -2.0 and max(x) <= 362.0):
           pyplot.xlim(0.0,360.0)   
           xticks = numpy.arange(0,450,90)
           plt.xticks(xticks)
#           print 'modifying ticks for azimuth'
  #      else:
  #         print 'WARNING: TICKS NOT MODIFIED FOR AZIM'
  #         print xlbl, x
#     else:
#        print xlbl
#   plt.xticks([1,2,3,4])

#for fan plots
#     pyplot.ylim(0,10)
#=======================================================================
# create x and y labels, set ylabel rotation to 0
#=======================================================================

   if plot_type == 'line':
     if(xlbl != ''):
        plt.xlabel(xlbl,fontsize=fs, fontweight='bold')
        ax.xaxis.set_label_coords(0.5 ,-0.075)

   if(ylbl != ''):
      plt.ylabel(ylbl,fontsize=fs, fontweight='bold', rotation=0, ha='left')
   
#=======================================================================
# move the labels to non-intrusive locations on plot
#=======================================================================

   ax.yaxis.set_label_coords(0.0,1.05)
#   ax.yaxis.set_label_position("right")

 #  ax.set_xticks(np.arange(0,6,1))
 #  ax.set_yticks(np.arange(0,6,1))
 #  label = ax.set_xlabel('xlabel', ha='left', va = 'top', )#fontsize = 9)

#   ax.set_xmargin(0.01)
#=======================================================================
# make x,y grid values bold     
#=======================================================================

   a           = gca()
   xticks      = a.get_xticks()
   yticks      = a.get_yticks()

#=======================================================================
# for power plots, use zero as lower limit
#=======================================================================

   if(ylbl.lower().startswith('power')):
      pyplot.ylim(0,float(yticks[-1]))

   yticks      = a.get_yticks()
   xtick2      = []
   ytick2      = []
   if plot_type == 'line':
      for tick in xticks:
        xtick2.append(format_number_AS(tick))
      a.set_xticklabels(xtick2, font)
  
   for tick in yticks:
      ytick2.append(format_number_AS(tick))

   a.set_yticklabels(ytick2, font)

#=======================================================================
# Plot annotation and bounds
#=======================================================================

   sizes = fig.get_size_inches() 

#=======================================================================
#Add plot annotation if requested
#=======================================================================

   if len(tag)>0:
 
#=======================================================================
#pad the string with leading spaces
#=======================================================================

      pad = '' ; spc = ' '
      for i in range(-2,len(tag)):
         pad = pad + spc
#            tag = pad+tag        
                 
#=======================================================================
#Draw an aesthetically pleasing arrow (30 or 45 deg line)
#=======================================================================

      ymin  = min(y) ; ymax = max(y) ; xmin = min(x); xmax = max(x);
      dy    = 0.1*(ymax-ymin) ; dx   = 0;
 
      plt.annotate(tag, color=line_color, xy=(x[-1],y[-1]), 
            xytext=(x[len(x)-1], y[len(y)-1]), fontsize=18, fontweight='bold')
#            arrowprops=dict(facecolor=clr, edgecolor = clr, 
#            arrowstyle="-|>",linewidth=0.004) )
      
      fig.set_size_inches(sizes[0]+1,sizes[1]+0.2)

   else:
      fig.set_size_inches(sizes[0],sizes[1]+0.2)

#=======================================================================
# turn on the grid      
#=======================================================================

   grid(True)

#=======================================================================
# Remove the "bounding box" lines to the top and right:
#=======================================================================

   for side in ['right','top']:
      ax.spines[side].set_visible(False)

#=======================================================================
# end of operations: return axes handle for additional postprocessing
#=======================================================================

   return a,ax,font

import os

#=============================================================================
# Save png
#=============================================================================

def save_png(filename):

  plt.tight_layout()
  try:
    plt.tight_layout(pad=0.1)
  except:
    pass
  try:
    os.remove(filename)
  except:
    pass

  plt.savefig(filename+'.png',format='png',dpi=300)
#  except:
#      plt.savefig(filename+'.png',format='png',dpi=300)
#  
  return 


#=============================================================================
# save as eps
#=============================================================================

def save_eps(filename):

  plt.tight_layout()
  try:
      plt.tight_layout(pad=0.1)
  except:
      pass
  try:
    os.remove(filename)
  except:
    pass

  plt.savefig(filename+'.jpg',dpi=150)
#  except:
#      plt.savefig(filename+'.png',format='png',dpi=300)
#  
  return 

import decimal
import random

#=============================================================================
# number formatting function by AS to remove trailing zeros
#=============================================================================

def format_number_AS(in_num):

    try:
        dec = decimal.Decimal(in_num)
    except:
        print 'neither integer nor float: what do you want to format? '
        return num
        #return 'bad programmer: BAD! enter a number only'

#    if len(str(num)) < 4:
#        return str(num)
    
    num             = float(in_num)
    outstr          = str(num)
    int_part        = int(num)
    dec_part        = num - float(int_part)
    nmax            = 6

    if dec_part == 0.0:
        return str(int_part)
    else:
#find how many digits exist in integer part
#        print 'original number', num
#        print 'integer part', int_part
#        print 'decimal part ',dec_part
        int_str  = str(int_part)

        if dec_part < 0.0 and int_part == 0:
            int_str = '-0'

        dec_part = round(dec_part,nmax)
#        print 'integer string',int_str
        fl_part  = str(dec_part).split('.')
        if len(fl_part) > 1:
            fl_part = fl_part[1]
        else:
            fl_part = fl_part[0]
        nint     = len(int_str)
#        print 'floating point string',fl_part
#have at most six characters in xlabel
        if nint >= nmax:
            return int_str
#            print 'returning integer string'
        else:
            for i in xrange(len(fl_part),5):
                fl_part = fl_part + '0'
            out_str = int_str + '.' + fl_part[0:nmax-1-nint]
            out_str = out_str.rstrip('0')
            out_str = out_str.rstrip('.')
        
    return outstr
#======================================================================
# Python function that creates a contour plot with auto label generation
#
#  Dependencies: pylab, numpy, matplotlib, AS custom number formatting
#
#
#======================================================================

#AS Jan 25th 2016 @ 35k ft: halfway from SFO to SEA
#Imports moved to top of function and condensed: 
#import pylab as pl
#import numpy as np
#from pylab import *
#from matplotlib import pyplot
#from line_plots import format_number_AS   
#import matplotlib.cm as cm
#======================================================================
#  Inputs:
#======================================================================
#        (a) x, y       : Cartesian coordinates of points where data is available
#        (b) z          : Value of data. X,Y,Z must be of the same size, 2d numpy arrays
#        (c) fid        : figure id number. Mostly a dummy number since we usually
#                         do not overlay contours
#        (d) title_str  : File name to use to save file as png
#        (e) zmin       : lower limit cutoff
#        (f) zmax       : upper limit cutoff. The data is saturated at
#                         zmin if it is lower than zmin, and similarly
#                         capped at zmax if it is larger than zmax
#
#======================================================================
# Outputs : None
#======================================================================

# def make_contour(x,y,z,fid,title_str,zmin,zmax):

# #=======================================================================
# # Data preprocessing
# #=======================================================================

#    zabs  = max(abs(zmin),abs(zmax))

# #=======================================================================
# # If either one is zero, or all values are +ve, or all values are -ve
# #=======================================================================

#    if zmin*zmax >= 0:
#       zmin2       = zmin
#       zmax2       = zmax 
#       map_style   = 'Greys_r'

# #=======================================================================
# # If some values are on one side of zero and some on the other
# #=======================================================================

#    else:
#       zmin2       =-zabs
#       zmax2       = zabs
#       map_style   = 'Greys'#'seismic'

# #   map_style      = 'grays'

#    zmin2 = 0
#    zmax2 = 0.1
# #=======================================================================
# # Saturate z to avoid blank regions in contour plots
# #=======================================================================

#    zmin2       = round(zmin2,2)
#    zmax2       = round(zmax2,2)
#    z           = numpy.clip(z,zmin2,zmax2)

# #=======================================================================
# # Set ticks 
# #=======================================================================
   
#    ticks = numpy.linspace(zmin2,zmax2,5)
#    level = numpy.linspace(zmin2,zmax2,num=51)

# #=======================================================================
# # Set fonts
# #=======================================================================

#    plt.rc('text', usetex=True)
#    font = {'family' : 'Times',
#            'weight' : 'bold',
#            'size'   : 20}
#    matplotlib.rc('font', **font)

#    fig   = figure(fid)
#    plt.clf()

# #=======================================================================
# # Call contour plot maker
# #=======================================================================

#    plt.contourf(x,y,z,levels=level,cmap=map_style)
   
#    cbar  = plt.colorbar()
# #   tick2 = []
# #   for tick in ticks:
# #      tick2.append(round(tick,4))

#    cbar.set_ticks(ticks)
#    cbar.set_ticklabels(ticks)
   
# #===========================================================================
# # Format axes
# #===========================================================================

#    a           = gca()
#    xticks      = a.get_xticks()
#    yticks      = a.get_yticks()
#    xtick2      = []
#    ytick2      = []
#    for tick in xticks:
#       xtick2.append(format_number_AS(tick))
#    for tick in yticks:
#       ytick2.append(format_number_AS(tick))
#    a.set_xticklabels(xtick2, font)
#    a.set_yticklabels(ytick2, font)

#    return None
