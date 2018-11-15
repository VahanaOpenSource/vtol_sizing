#====================================================================
# Program to postprocess design sweep
#====================================================================

import numpy, copy
import os, shutil, sys

#====================================================================
# insert path
#====================================================================
   
#sys.path.insert(0,'/home/yosemite/Dropbox/py/packages/plotting/')
sys.path.insert(0,'../../src/Python/Stage_1/')

#====================================================================
# import modules
#====================================================================

from data_processing    import append_plot_data, final_preparations

#====================================================================
# Set name of file
#====================================================================

path_list       = ['./']

#====================================================================
# Define and create folder for output images
#====================================================================

all_data       = {}

#icol           =  9             # rank by GTOW
#icol           = 10             # rank by Power installed
#icol           = 11             # rank by smallest rotor
#icol           = 16             # rank by fuel/battery weight 
#icol           = 17             # rank by empty weight
icol           =-19            # rank by max payload
#icol           = 20            # rank by operating cost/hr
#====================================================================
#16 = fuel/battery, 
#10 = GTOW, 
#11 = Power, 
#12 = Radius
#====================================================================

constraints     = {'b': 14.0, 'DL': 20, 'CTsigma': 0.14}#, 'wingAR':8}

#====================================================================
# loop over plot data paths and create sets
#====================================================================

for path in path_list: 
   k, yl       = append_plot_data(path, all_data, icol, constraints)

#====================================================================
# Remove internal dictionary counter 
#====================================================================

final_preparations(all_data)