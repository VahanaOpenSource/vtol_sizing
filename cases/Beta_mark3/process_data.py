#====================================================================
# Program to postprocess design sweep
#====================================================================

import numpy, copy
import os, shutil, sys
from sys import argv

#====================================================================
# insert path
#====================================================================
   
#sys.path.insert(0,'/home/yosemite/Dropbox/py/packages/plotting/')
sys.path.insert(0,'../../src/Python/Stage_1/')

#====================================================================
# import modules
#====================================================================

from extract_plot_data  import *
from run_sim            import create_dir
from data_processing    import append_plot_data, final_preparations

#====================================================================
# Set name of file
#====================================================================

path_list       = ['./']

#====================================================================
# Define and create folder for output images
#====================================================================

all_data       = {}
     
#====================================================================
#icol           =  1             # rank by GTOW
#icol           =  2             # rank by Power installed
#icol           =  3             # rank by fuel/battery weight 
#icol           =  4             # rank by empty weight
#icol           = -5             # rank by max payload
icol           =   6             # rank by operating cost/hr
#====================================================================

constraints     = {'DL': 20, 'CTsigma': 0.141}

if len(argv) < 2:
    file_in         = 'summary.dat'
    file_out        = 'best_design.dat'
elif(len(argv) == 2):
    if(argv[1] == 'optim'):
        file_in     = 'optim_summary.dat'
        file_out    = 'optim_ranked.dat'
        print('ranking optmized designs')
else:
    quit('ERROR: usage is python process_data.py or python process_data.py optim')

#====================================================================
# loop over plot data paths and create sets
#====================================================================

for path in path_list: 
   k, yl       = append_plot_data(path, file_in, file_out, all_data, icol, constraints)

#====================================================================
# Remove internal dictionary counter 
#====================================================================

final_preparations(all_data)