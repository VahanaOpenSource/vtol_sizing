#====================================================================
# Program to postprocess design sweep over rotor geometry
#====================================================================

import numpy, copy
import os, shutil, sys
#from line_plots         import plot_wrapper
#from extract_plot_data  import *

home                 	= os.path.expanduser('~') 
add_path 				= home + '/Dropbox/py/packages/batchrun/'
sys.path.insert(0,add_path)

from run_sim            import create_dir
from data_processing    import append_plot_data_geom, final_preparations

#====================================================================
# Set name of file
#====================================================================

#FOR ASYM:simple vs CSD 240
#path_list      = ['performance/twist_sweep/']
path_list 		= ['../../Outputs/Stage_1/']
#====================================================================
# Define and create folder for output images
#====================================================================

images_path    = 'geom_sweep_images/'
all_data       = {}
create_dir(images_path)

#====================================================================
# loop over plot data paths and create sets
#====================================================================

for path in path_list: 
   k, yl       = append_plot_data_geom(path, all_data)

#====================================================================
# Remove internal dictionary counter 
#====================================================================

final_preparations(all_data)
quit()
#====================================================================
# Call plot wrapper function for each parameter perturbation, all sets
#====================================================================

#for key,data in all_data.iteritems():
#   folder_name       = images_path + 'sweep_' + str(key)
#   plot_wrapper(data, k, yl, folder_name)