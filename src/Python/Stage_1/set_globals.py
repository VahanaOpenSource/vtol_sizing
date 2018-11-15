#=========================================================================
#=========================================================================
# Declare global variables: import this module and access dictionary
# modules
#=========================================================================
#=========================================================================
import numpy 

def init_globals(xaxis, zaxis):

#=========================================================================
# set parameters
#=========================================================================

   global_set                    = {}

#=======================================================================
# define x axis
#=======================================================================

   xstrings                      = ['Speed', 'Payload', 'Disk Loading', 'Rotor Solidity']

   units                         = ['(knots)', '(kg)', '(lb/sq.ft)','']

   zstrings                      = ['V = ', 'Payload = ','DL = ','\sigma = ']

   data_keys                     = ['speed_array','payload_array','DL_array','sigma_array']

   subset_keys                   = ['speed', 'payload','disk loading', 'solidity']

#=======================================================================
# Find x axis value
#=======================================================================

   for val, key, unit, zstring, skey in zip(xstrings, data_keys, units, zstrings, subset_keys):
      LHS                        = xaxis.lower().replace(" ","")
      RHS                        = val.lower().replace(" ","").lower()
      if LHS in RHS or RHS in LHS:
         
         xstr                    = val + ' ' + unit
         global_set['xkey']      = key
         global_set['xsubkey']   = skey
         print 'x axis string is ',xaxis
         print 'data key name is ',key 
         print 'units are        ',unit 
   
#=======================================================================
# Find z axis value
#=======================================================================

      LHS                        = zaxis.lower().replace(" ","")
      if LHS in RHS or RHS in LHS:
         global_set['zkey']      = key
         global_set['zsubkey']   = skey
         global_set['zprefix']   = zstring
         print 'z axis string is ',zstrings[xstrings.index(val)]
         print 'data key name is ',key 

#=======================================================================
# Create a dictionary for data storage/processing
#=======================================================================

   global_set['blade_AR_array']  = []
   global_set['DL_array']        = []
   global_set['Nb']              = 0
   global_set['payload_array']   = []
   global_set['sigma_array']     = []
   global_set['speed_array']     = []
   global_set['xstr']            = xstr
   global_set['zstr_array']      = []
   global_set['strings']         = xstrings

#=========================================================================
# End of operations
#=========================================================================

   return global_set

#=========================================================================
#=========================================================================
# Function to process additional information and set plot arrays
#=========================================================================
#=========================================================================

def mod_globals(global_set, a1, a2, a3, a4, a5):

#=========================================================================
# Remember data in global dictionary
#=========================================================================

   global_set['DL_array']        = a1
   global_set['blade_AR_array']  = a2
   global_set['speed_array']     = a3
   global_set['payload_array']   = a4
   global_set['Nb']              = a5

#=======================================================================
# Compute solidity and remember
#=======================================================================
   
   Nb                      = global_set['Nb']
   blade_AR_array          = global_set['blade_AR_array']
   speed_array             = global_set['speed_array']
   payload_array           = global_set['payload_array']
   DL_array                = global_set['DL_array']
   zprefix                 = global_set['zprefix']
   for blade_AR in global_set['blade_AR_array']:
      this_sigma           = float(Nb)/numpy.pi/blade_AR
      global_set['sigma_array'].append(this_sigma)
      global_set['zstr_array'].append(zprefix + str(this_sigma))
      
   sigma_array             = global_set['sigma_array']

#=======================================================================
# check inputs: no more than two-parameter sweeps allowed
#=======================================================================

   all_arrays                 = [   speed_array, payload_array, DL_array,     \
                                 blade_AR_array]
   plot_arrays                = [   speed_array, payload_array, DL_array,     \
                                    sigma_array]
   global_set['all_arrays']   = all_arrays
   global_set['plot_arrays']  = plot_arrays
   num_loops                  = 0
   for l in all_arrays:
      if len(l) > 1:
         num_loops            = num_loops + 1

   if num_loops > 2:   
      print 'CRITICAL ERROR: CANNOT SWEEP OVER > 2 PARAMETERS!'
      print 'PROGRAM WILL EXIT'
      quit()

#=========================================================================
# create subset data for each individual run
#=========================================================================

   subset                  = {'speed':0,  'payload':0.0, 'disk loading':0.0, 
                                 'Nb':0, 'solidity':0.0}

   global_set['subset']    = subset

#=========================================================================
# End of operations
#=========================================================================

   return None

#=======================================================================
#=======================================================================
# Set plot properties for each line
#=======================================================================
#=======================================================================

def set_plot_prop(data, iz, global_set):

#=======================================================================
# Define line styles
#=======================================================================

   zval                    = global_set['zval']

   line_colors             = ['r','b','k','m','firebrick','purple']
   line_styles             = ['-','--','-','--','-','--']
   marker_styles           = ['','','o','s','>','^']

   data['x']               = global_set[global_set['xkey']]
   data['xlbl']            = global_set['xstr']
   izplot                  = iz
   data['lcolor']          = line_colors[iz]
   data['lstyle']          = line_styles[iz]
   data['mtype']           = marker_styles[iz]
   ztruncate               = numpy.around(zval,3)
   data['tag']             = global_set['zstr_array'][iz]
   print data['tag']

   return None

#=======================================================================
#=======================================================================
# Set design parameters for single run of design code
#=======================================================================
#=======================================================================

def set_conditions(global_set):

#=======================================================================
# default values
#=======================================================================

   subset                  = global_set['subset']
   subset['speed']         = global_set['speed_array'][0]
   subset['payload']       = global_set['payload_array'][0]
   subset['disk loading']  = global_set['DL_array'][0]
   subset['solidity']      = global_set['sigma_array'][0]
   subset['Nb']            = global_set['Nb']

   xskey                   = global_set['xsubkey']
   zskey                   = global_set['zsubkey']

   subset[xskey]           = global_set['xval']
   subset[zskey]           = global_set['zval']

   return subset