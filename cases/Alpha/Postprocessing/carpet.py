#=======================================================================
# python script to create contour plots of resulting design aspects
# with variation of rotor and wing parameters
#=======================================================================

import plotly, numpy, sys
sys.path.append('../../src/Python/Postprocessing/')
sys.path.append('../../src/Python/Stage_1/')

from draw_carpet  import draw_a_carpet as dac 
from wing_carpets import wing_carpets as wc 
from hydra_summary_map  import hydra_summary_map as hsmap

#=======================================================================
# define file name
#=======================================================================

fname             = 'summary.dat'
opt_name          = 'best_design.dat'

#=======================================================================
# read files for anchor point information and all data
#=======================================================================

with open(opt_name,'r') as f:
    header      = f.readline()
    line        = next(f)
    vals        = line.split()

for ival,val in enumerate(vals):
    vals[ival]  = float(val)

data            = numpy.loadtxt(fname,skiprows=1)
ndat            = numpy.shape(data)[0]
valid           = data[:,-1]

maps            = hsmap()

#=======================================================================
# read all rows in fname
#=======================================================================

with open(fname,'r') as f:
   lines        = f.readlines()

#=======================================================================
# determine which columns to freeze
#=======================================================================

ab_rotor        = [maps['DL'],maps['sigma'],maps['vtip']]            # a,b axes for rotor design variables    
xy_rotor        = [maps['DL'],maps['CTsigma']]          # x,y axes for rotor design 
ab_wings        = [maps['wing_AR'],maps['wing_CL'],maps['vtip']]     # a,b axes for wing  design variables
xy_wings        = [maps['wing_b'],maps['wing_CL']]      # x,y axes for wing  design

#=======================================================================
# find axes names
#=======================================================================

ref_rotor       = []
ref_wings       = []
for param in ab_rotor:
   ref_rotor.append(vals[param])

for param in ab_wings:
   ref_wings.append(vals[param])

print(ref_rotor,ref_wings)
#=======================================================================
# isolate vehicle designs for fixed wing  designs 
# isolate vehicle designs for fixed rotor designs 
#=======================================================================

target_cols     = ab_wings
target_vals     = ref_wings
rotor_array     = []; rotor_rows    = [];
wings_array     = []; wing_rows     = [];

for idat in range(ndat):               # loop over all rows
   line        = data[idat,:]

   valid_line  = True 
   for imatch,match in zip(ab_wings,ref_wings):       # loop over all target cols
      # print(imatch,match)
      if line[imatch] != match:
         valid_line = False; break 

   if valid_line:
      rotor_array.append(line); rotor_rows.append(idat)

   valid_line  = True 
   for imatch,match in zip(ab_rotor,ref_rotor):       # loop over all target cols
      if line[imatch] != match:
         valid_line = False; break 

   if valid_line:# and line[17] < 0.14:
      wings_array.append(line); wing_rows.append(idat)

#=======================================================================
# convert to numpy arrays and save to disk
#=======================================================================

dat1            = numpy.asarray(rotor_array)
dat2            = numpy.asarray(wings_array)
valid1          = dat1[:,maps['valid']]
valid2          = dat2[:,maps['valid']]

with open('rotor_summary.dat','w') as r:
   r.write(header)
   for row in rotor_rows:
      r.write(lines[row+1])

with open('wings_summary.dat','w') as w:
   w.write(header)
   for row in wing_rows:
      w.write(lines[row+1])

#=======================================================================
# get parameters to plot
#=======================================================================

DL              = dat1[:,maps['DL']]
sigma           = dat1[:,maps['sigma']]
Wt              = dat1[:,maps['Weight']]
Power           = dat1[:,maps['Power']]
Radius          = dat1[:,maps['Radius']]
Span            = dat1[:,maps['wing_b']]
Fuel            = dat1[:,maps['Fuel']]
CTsigma         = dat1[:,maps['CTsigma']]
Cost            = dat1[:,maps['Cost']]

#=======================================================================
# draw variation of vehicle metrics with rotor design parameters
#=======================================================================

ystr = '<b> Hover blade loading (CT/sigma) <b>'
xstr = '<b> Disk loading (lb/sq.ft) <b>'

dac(DL,sigma,Wt   ,DL,CTsigma,'weight1'   ,valid1, 'Take-off mass (kg)',xstr, ystr)
dac(DL,sigma,Power,DL,CTsigma,'power1'    ,valid1, 'Power (kW)',xstr, ystr)
dac(DL,sigma,Fuel ,DL,CTsigma,'battery1'  ,valid1, 'Battery (kg)',xstr, ystr)
dac(DL,sigma,Cost ,DL,CTsigma,'cost1'     ,valid1, 'Cost (USD/hr)',xstr, ystr)

#=======================================================================
# draw variation of vehicle metrics with wing design parameters
#=======================================================================

DL              = dat2[:,maps['DL']]
sigma           = dat2[:,maps['sigma']]
Wt              = dat2[:,maps['Weight']]
Power           = dat2[:,maps['Power']]
Radius          = dat2[:,maps['Radius']]
Span            = dat2[:,maps['wing_b']]
Fuel            = dat2[:,maps['Fuel']]
Cost            = dat1[:,maps['Cost']]
CTsigma         = dat2[:,maps['CTsigma']]
AR              = dat2[:,maps['wing_AR']]
CL              = dat2[:,maps['wing_CL']]

Area            = numpy.multiply(numpy.reciprocal(AR),numpy.square(Span))
WL              = numpy.multiply(numpy.reciprocal(Area),Wt)*2.2*0.3048**2       # wing loading

ystr = '<b> Wing loading (lb/sq.ft) </b>'
xstr = '<b>Wing span (m) </b>'

wc(AR, CL, Wt   , Span, WL, 'weight2'  ,valid2, 'Take-off mass (kg)',xstr,ystr)
wc(AR, CL, Power, Span, WL, 'power2'   ,valid2, 'Power (kW)',xstr,ystr)
wc(AR, CL, Fuel , Span, WL, 'battery2' ,valid2, 'Battery (kg)',xstr,ystr)
wc(AR, CL, Cost , Span, WL, 'cost2'    ,valid2, 'Cost (USD/hr)',xstr,ystr)