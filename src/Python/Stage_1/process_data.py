#====================================================================
# Program to postprocess design sweep
#====================================================================

import numpy, copy
import os, shutil, sys

#====================================================================
# insert path
#====================================================================
   
sys.path.insert(0,'/home/yosemite/Dropbox/py/packages/plotting/')

#====================================================================
# import modules
#====================================================================

from line_plots         import plot_wrapper
from extract_plot_data  import *
from run_sim            import create_dir
from data_processing    import append_plot_data, final_preparations

#====================================================================
# Set name of file
#====================================================================

#path           = 'sh60/'
#path           = 'asym/'
#path           = os.getcwd()
#path           = 'coax/'
#path           = 'results/smr/V160/Nb6/'
#FOR SMR: 160
#path_list      = ['results/smr/V160/Nb3/','results/smr/V160/Nb4/',   \
#                  'results/smr/V160/Nb5/','results/smr/V160/Nb6/']

#FOR COAX: simple vs. CSD, 220 knots
#path_list      = ['results/coax/V220/Nb4/',                          \
#                  'results/CSD/coax/V220/Nb4/0deg_twist/',           \
#                  'results/CSD/coax/V220/Nb4/5deg_twist/',           \
#                  'results/CSD/coax/V220/Nb4/10deg_twist/',          \
#                  'results/CSD/coax/V220/Nb4/12p5deg_twist/']

#FOR COAX vs ASYM 200
#path_list      = ['results/coax/V200/Nb4/','results/asym/V200/Nb4/']

#FOR COAX vs ASYM 220
#path_list      = ['results/coax/V220/Nb4/','results/asym/V220/Nb4/']

#FOR ASYM:simple vs CSD 240
#path_list      = [#'results/asym/V240/Nb4/',                    \
#                  'results/CSD/asym/V240/Nb4/0deg_twist/',     \
#                  'results/CSD/asym/V240/Nb4/5deg_twist/',     \
#                  'results/CSD/asym/V240/Nb4/10deg_twist/']

#path_list      = ['smr_nub/Nb4/','coax_nub/Nb4/']
#path_list       = ['../../Outputs/Stage_1/ERF_data/half_wing_flap/']
path_list       = ['../../Outputs/Stage_1/']
#FOR coax: simple vs CSD 240
#path_list      = ['results/coax/V240/Nb4/',                    \
#                  'results/CSD/coax/V240/Nb4/0deg_twist/',     \
#                  'results/CSD/coax/V240/Nb4/5deg_twist/',     \
#                  'results/CSD/coax/V240/Nb4/10deg_twist/']

#path_list       = ['results/coax/V240/Nb3/','results/coax/V240/Nb4/',
#                  'results/coax/V240/Nb5/','results/coax/V240/Nb6/']
#COAX VS ASYM
#path_list      = ['results/coax/V240/Nb4/','results/asym/V240/Nb4/']
#path_list      = ['results/coax/V240/Nb3/','results/coax/V240/Nb4/',   \
#                  'results/coax/V240/Nb5/','results/coax/V240/Nb6/']

#FOR OTHERS
#path_list      = ['results/coax/V200/Nb4/','results/coax/V240/Nb4/', \
#                  'results/coax/V220/Nb4/','results/coax/V230/Nb4/', \
#                  'results/coax/V210/Nb4/']#,   \

#====================================================================
# Define and create folder for output images
#====================================================================

images_path    = '../../Outputs/Stage_1/images/'
all_data       = {}
create_dir(images_path)

#====================================================================
# loop over plot data paths and create sets
#====================================================================

for path in path_list: 
   k, yl       = append_plot_data(path, all_data)

#====================================================================
# Remove internal dictionary counter 
#====================================================================

final_preparations(all_data)

quit()
#====================================================================
# Call plot wrapper function for each parameter perturbation, all sets
#====================================================================

for key,data in all_data.iteritems():
   folder_name       = images_path + 'sweep_' + str(key)
   plot_wrapper(data, k, yl, folder_name)

# #====================================================================
# # For bar graphs
# #====================================================================

# big_data             = {}; keys     =  []
# sub_data             = {}; labels   =  []

# #====================================================================
# # Essentials for plotting
# #====================================================================

# sub_data['x']        = numpy.linspace(1,len(unique_rowid),len(unique_rowid))#unique_ids
# sub_data['xlbl']     = 'Case number'
# sub_data['lstyle']   = '-'
# sub_data['lcolor']   = 'r'
# sub_data['mtype']    = 'o'
# sub_data['tag']      = ''
# sub_data['plot_type']= 'bar'

# #====================================================================
# # Actual data values to be plotted
# #====================================================================

# sub_data['GTOW']     = all_vals[unique_rowid,7] ; keys.append('GTOW')   ; labels.append('Take-off Weight (lb)')
# sub_data['Power']    = all_vals[unique_rowid,8] ; keys.append('Power')  ; labels.append('Power (hp)')
# sub_data['Radius']   = all_vals[unique_rowid,9] ; keys.append('Radius') ; labels.append('Rotor Radius (ft)')
# sub_data['Fuel']     = all_vals[unique_rowid,13]; keys.append('Fuel')   ; labels.append('Fuel Weight (lb)')
# sub_data['Empty']    = all_vals[unique_rowid,14]; keys.append('Empty')  ; labels.append('Empty Weight (lb)')

# #====================================================================
# # Place
# #====================================================================

# big_data['SMR']      = sub_data

# #====================================================================
# # Call plot wrapper
# #====================================================================

# plot_wrapper(big_data, keys, labels, images_path.rstrip('/'))

#====================================================================
# End of operations
#====================================================================
