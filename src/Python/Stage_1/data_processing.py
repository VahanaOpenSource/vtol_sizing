#====================================================================
# Read data from one run and return list
#====================================================================
import copy, os, numpy, sys
from extract_plot_data import write_header, write_to_file, create_subset
from locator            import *

def run_read(path,filename):


#====================================================================
# Get file name
#====================================================================

   fname          = os.path.join(path,filename)

#====================================================================
# Open and find number of rows and columns
#====================================================================

   nlines         = -1
   with open(fname,'r') as f:
      for lid,line in enumerate(f):

         nlines   = nlines + 1

#====================================================================
# get header line
#====================================================================

         if lid == 0:
            header = line
         else:
            if lid == 1:
               temp  = line.rstrip('\n').split()
               ncols = len(temp)

#====================================================================
# Initialize array and list
#====================================================================

   all_vals    = numpy.zeros((nlines,ncols+1))     # +1 cols for buffer
   list_vals   = []

#====================================================================
# Open file and load information into return data types
#====================================================================

   with open(fname,'r') as f:
      for lid,line in enumerate(f):
         if lid > 0:
            temp        = line.rstrip('\n').lstrip(' ').rstrip(' ').split()
            this_list   = []
            for ival,val in enumerate(temp):
               this_list.append(float(val))
               all_vals[lid-1,ival] = float(val)

            list_vals.append(this_list)

#====================================================================
# End of operations
#====================================================================

   return list_vals, all_vals, ncols, header

#====================================================================
# Function to check if a design is valid
#====================================================================

def check_valid(row_data, constraints):

#====================================================================
# max limits: CT sigma, disk loading and span
#====================================================================

   try:
      CTsig_max   = constraints['CTsigma']
   except:  
      CTsig_max   = 0.16

   try:
      DL_max      = constraints['DL']
   except:
      DL_max      = 1e6

   try:
      bmax        = constraints['b']
   except:
      bmax        = 1e10

   try:
      ARmax       = constraints['wingAR']
   except:
      ARmax       = 40

#====================================================================
# interpret row data
#====================================================================

   iout           = 0
   GTOW           = row_data[iout+1]
   Power          = row_data[iout+2]
   valid          = row_data[iout+7]

   if valid == 1:
      row_valid   = True 
   else:
      row_valid   = False 

   return row_valid

#====================================================================
# Function to return row indices for best designs
#====================================================================
from performance_sweep import write_header as wh, write_to_file as wtf

def find_best(path, f_out, all_vals, nbest, icol, h, constraints):

#====================================================================
# Find indices for top "nbest" entries according to various criteria
#====================================================================

#      id_sort     = list(numpy.argsort(all_vals[:,icol])[::-1])
#   else:
   if icol < 0:          # payload: find max 
      reverse  = True 
      icol     = abs(icol)
   else:
      reverse  = False 

   id_sort     = list(numpy.argsort(all_vals[:,icol]))

   if reverse:          # payload: find max 
      id_sort  = id_sort[::-1]

#   id_GTOW        = list(numpy.argsort(all_vals[:,7] ))#[0:nbest]
#   id_Power       = list(numpy.argsort(all_vals[:,8] ))#[0:nbest]
#   id_Radius      = list(numpy.argsort(all_vals[:,9] ))#[0:nbest]
#   id_Fuel        = list(numpy.argsort(all_vals[:,13]))
   
#   print id_Fuel
#   id_Fuel        =  id_Fuel[0:nbest]
#====================================================================
# Filter out designs with CT/sigma > some value
#====================================================================

   best           = []
   
#====================================================================
# Min fuel weight
#====================================================================

#   print numpy.shape(all_vals)
   ibest          = -1
   for irow in id_sort:

      valid_row   = check_valid(all_vals[irow,:], constraints) 
      if valid_row:
         best.append(irow)

         if ibest == -1:
            ibest = irow

#====================================================================
# Take first "nbest" designs
#====================================================================

   nrows          = numpy.shape(best)[0]
   nmax           = min(nrows, nbest)
   best           = best[0:nmax]

#====================================================================
# Write best designs to file
#====================================================================

   best_designs_file    =  os.path.join(path,f_out) #'best_design.dat'

   with open(best_designs_file,'w') as f:

      f.write('%s' % h)      
      case, data     = create_subset(all_vals[ibest,:])

      for row in best:

         case, data     = create_subset(all_vals[row,:])
         write_to_file(case, data, f)

#====================================================================
# End of operations
#====================================================================

   return best, ibest

#====================================================================
# Function to get unique data rows after applying a filter
#====================================================================

def get_unique_rows(path, all_vals, all_rows, indices, uids, h):

#====================================================================
# Filter out columns based on index list
#====================================================================

   big_list             = []
   for rowid in all_rows:
      this_list         = []
      for colid in indices:
         this_list.append(all_vals[rowid,colid])
      big_list.append(tuple(this_list))

#====================================================================
# Find unique rows
#====================================================================

   unique_list          = set(big_list)

#====================================================================
# Map rows of unique_list back to big_list
#====================================================================

   unique_rowid         = []
   for row in unique_list:
      unique_rowid.append(big_list.index(row))

#====================================================================
# Write data to file
#====================================================================

   if len(all_vals[0,:]) > 9:                
      reduced_data_file    =  os.path.join(path,'reduced_data.dat')
      with open(reduced_data_file,'w') as f:

         f.write('%s' % h)
         for idx in unique_rowid:
         
            case, data     = create_subset(all_vals[idx,:])
            write_to_file(case, data, f)

#====================================================================
# Get parameter lists for making 1d plots
#====================================================================

   parameter_lists      = []
   for icol in indices:
      this_list         = numpy.unique(all_vals[:,icol])
      parameter_lists.append(this_list)

#====================================================================
# truncate parameter lists: temp only!
#====================================================================

#   parameter_lists[0]   = parameter_lists[0][0:6]

#====================================================================
# End of operations
#====================================================================

   return  unique_list, unique_rowid, parameter_lists

#====================================================================
# Function to get row indices for single parameter perturbations 
#====================================================================

def get_rowids(indices,  template_row, parameter_lists, unique_list, \
               unique_rowid, all_rows, nmatch):

#====================================================================
# Separate data into lists for single-parameter perturbations
#====================================================================

   rowids               = []

   for icol in indices:                   # loop over all parameters

#====================================================================
# Initialize match set
#====================================================================

      data_set          = template_row[indices]
      rowid             = []

#====================================================================
# search through big_list 
#====================================================================

      for param_val in parameter_lists[icol]:
         data_set[icol] = param_val                # data set to match 

         for irow,row in enumerate(unique_list):   # loop over unique_list
            if row[0:nmatch] == tuple(data_set):  # if right entry 
               uniquelist_id  = irow
               biglist_id     = unique_rowid[irow]
               globalid       = all_rows[biglist_id]
               rowid.append(globalid)      

      rowids.append(rowid)

   return rowids

#====================================================================
# Pack 1d parameter sweeps into dictionaries for easy processing
#====================================================================

def pack_data(all_vals, indices, rowids, parameter_lists, all_data):

#====================================================================
# Initialize subsets of all_data
#====================================================================

   param_names                = ['DL','Solidity','Nb','Wing_AR','Wing_Loading','Vtip','Mtip','nu_beta']
   if not all_data:
      for param in param_names:
         all_data[param]         = {}
      all_data['counter']        = 0
   else:
      all_data['counter']       += 1

   iset                          = all_data['counter']
   sub_key                       = str(iset)

#====================================================================
# Define array of line styles and colors
#====================================================================

   line_styles                = ['-', '--', '-', '--', '-', '--']
   marker_styles              = ['o',  's', '^',  '>', '+', 'x' ]
   line_colors                = ['r',  'k', 'b',  'm', 'firebrick', 'purple']

#====================================================================
# Set plot properties
#====================================================================

   local_superset             = {}
   data                       = {}
   data['lstyle']             =   line_styles[iset]
   data['lcolor']             =   line_colors[iset]
   data['mtype']              = marker_styles[iset]
   data['tag']                = ''
   data['plot_type']          = 'line'
   xlabels                    = [' DL (lb/sq.ft)','Rotor Solidity',                 \
                                 ' Number of Blades','Wing Aspect Ratio',           \
                                 ' Wing lift fraction','Tip Speed (m/s)','Adv. Tip Mach Number','Flap Frequency']

   keys                       = ['GTOW','Power','Radius','Fuel','Empty','Ctsigma','RotorLO']
   ylabels                    = ['Take-off Weight (lb)','Power (hp)', 'Rotor Radius (ft)', \
                                 'Fuel Wt (lb)', 'Empty Weight (lb)','Blade loading','Lift Offset']

#   for key, value in all_data.iteritems():
#      print key
#====================================================================
# Pack data in appropriate sets/subsets
#====================================================================

   for icol, rowlist, xvals, param in zip(indices, rowids, parameter_lists, param_names):
      this_data                  = copy.copy(data)
      this_data['x']             = xvals
      xlbl                       = xlabels[icol]
      this_data['xlbl']          = xlbl

      iout                       = 1
      this_data['GTOW']          = all_vals[rowlist,iout]
      this_data['Power']         = all_vals[rowlist,iout+1]
      this_data['Fuel']          = all_vals[rowlist,iout+2]
      this_data['Empty']         = all_vals[rowlist,iout+3]

      all_data[param][sub_key]   = copy.deepcopy(this_data)
#      local_superset[icol]       = 
#====================================================================
# Store in superset
#====================================================================

#      all_data[iset]             = copy.deepcopy(this_data)
#      param_sweep_superset[icol] = copy.deepcopy(this_data)

#====================================================================
# End of operations
#====================================================================

   return keys, ylabels

#====================================================================
# Pack 1d parameter sweeps into dictionaries for easy processing
#====================================================================

def pack_data_geom(all_vals, indices, rowids, parameter_lists, all_data):

#====================================================================
# Initialize subsets of all_data
#====================================================================

   param_names                = ['x','thx','th1','y','croot','ctip']
   if not all_data:
      for param in param_names:
         all_data[param]         = {}
      all_data['counter']        = 0
   else:
      all_data['counter']       += 1

   iset                          = all_data['counter']
   sub_key                       = str(iset)

#====================================================================
# Define array of line styles and colors
#====================================================================

   line_styles                = ['-',  '-', '-', '--', '--', '--']
   marker_styles              = ['o',  's', '^',  '>', '+', 'x' ]
   line_colors                = ['r',  'k', 'b',  'r', 'k', 'b']

#====================================================================
# Set plot properties
#====================================================================

   local_superset             = {}
   data                       = {}
   data['lstyle']             =   line_styles
   data['lcolor']             =   line_colors
   data['mtype']              = marker_styles
   data['tag']                = ''
   data['plot_type']          = 'line'
   xlabels                    = [' Bilinear Twist Junction','Twist at Junction (deg)',   \
                                 ' Tip Twist (deg)','Bilinear Taper Junction',           \
                                 ' Chord at Junction (ft)','Tip Chord (ft)']

   keys                       = [['RotorPower','FxhubPower','TotalPower']]
   ylabels                    = ['Power (hp)']

#====================================================================
# Pack data in appropriate sets/subsets
#====================================================================

   for icol, rowlist, xvals, param in zip(indices, rowids, parameter_lists, param_names):
      this_data                  = copy.copy(data)
      this_data['x']             = xvals
      xlbl                       = xlabels[icol]
      this_data['xlbl']          = xlbl

      this_data['TotalPower']    = all_vals[rowlist,8]
      this_data['RotorPower']    = all_vals[rowlist,6]
      this_data['FxhubPower']    = all_vals[rowlist,8] - all_vals[rowlist,6]

      all_data[param][sub_key]   = copy.deepcopy(this_data)

#====================================================================
# Store in superset
#====================================================================

#      all_data[iset]             = copy.deepcopy(this_data)
#      param_sweep_superset[icol] = copy.deepcopy(this_data)

#====================================================================
# End of operations
#====================================================================

   return keys, ylabels

#====================================================================
# function to prepare all_data dictionary to plot
#====================================================================

def final_preparations(all_data):

   all_data.pop('counter',None)

#====================================================================
# End of operations
#====================================================================

   return None

#====================================================================
# Function to load run data into plottable dictionaries
# icol = column index to sort by
#====================================================================

def append_plot_data(path, f_in, f_out, all_data, icol, constraints):

#====================================================================
# Load data into a list
#====================================================================

   filename                      = f_in
   list_vals, all_vals, ncols, h = run_read(path,filename)

#====================================================================
# Get the best values (ranked increasing order of values)
#====================================================================

   nbest                = 1000
   best, ibest          = find_best(path, f_out, all_vals, nbest, icol, h, constraints)

#====================================================================
# Choose columns used to filter out unique data 
#====================================================================

   nmatch               = 0
   indices              = range(0,nmatch)
   uids                 = range(0,ncols);# del uids[6]

#====================================================================
# Choose design with minimum fuel requirement
#====================================================================

   template_row         = all_vals[ibest,:]

#====================================================================
# Get row indices for one-parameter perturbations around reqd. row
#====================================================================

   all_rows             = locator(all_vals, template_row, indices)

#====================================================================
# get list with filter applied for columns
#====================================================================

   unique_list, unique_rowid, parameter_lists =                      \
         get_unique_rows(path, all_vals, all_rows, indices, uids, h)

#====================================================================
# Get ids of row for single parameter perturbations
#====================================================================

   rowids               = get_rowids(indices, template_row,          \
                              parameter_lists, unique_list,          \
                              unique_rowid, all_rows, nmatch)

#====================================================================
# Prepare 1d parameter sweep plot data
#====================================================================

   k, yl                = pack_data(all_vals, indices, rowids,       \
                                      parameter_lists, all_data)

#====================================================================
# End of operations
#====================================================================

   return k, yl
