#====================================================================
# Function that finds "neighboring" data sets in all_data (2d array)
# that are "similar" to "template_data" list along column indices 
# given by "template_cols"
#
# Returns: 2d ARRAY
# value: row index in original data set 
# dimension 1: for multiple values of specified parameter
# dimension 2: each parameter specified by col id in template_cols
#====================================================================
import numpy

def locator(all_data, template_data, template_cols):

#====================================================================
# loop over big data set
#====================================================================
   
   nrows, ncols      = numpy.shape(all_data)       # size of data set
   ntcols            = len(template_cols)       # number of template cols

   rows              = []

#   print 'template columns are ',template_cols, nrows
#   print template_data

#====================================================================
# Loop over all rows
#====================================================================

   for irow in range(0,nrows):

#====================================================================
# Check if any element in this row and column in template_cols 
# matches the corresponding value in template_data
#====================================================================

      nmatch         = 0
      for icol in template_cols:
         imiss       = 0
         if abs(all_data[irow,icol] - template_data[icol]) < 1e-5:
            nmatch  += 1
         else:
            imiss    = icol 

#====================================================================
# If reqd. # dimensions match, remember this row index
#====================================================================

      if nmatch >= ntcols-1:
         rows.append(irow)

#====================================================================
# End of operations
#====================================================================

   return rows

#====================================================================
# Extract unique row ids corresponding using a certain col set ("cols")
# from a 2d array "all_vals" based on rows with index "all_row_ids" 
#====================================================================

def unique_rows(cols, all_vals, all_row_ids):


   unique_ids     = []
   unique_data    = []

#====================================================================
# Loop over all rows 
#====================================================================

   for idx in all_row_ids:
      this_list   = all_vals[idx,cols]

#====================================================================
# Check if a duplicate row exists
#====================================================================

      found_dupl  = False
      for l in unique_data:
         if set(this_list).intersection(l) == set(this_list):
            found_dupl = True
#            print 'FOUND DUPLICATE'
            break

#====================================================================
# If duplicate not found, add this row to the list of unique rows
#====================================================================

      if not found_dupl:
         unique_data.append(this_list)
         unique_ids.append(idx)

#====================================================================
# End of operations
#====================================================================

   return unique_ids, unique_data