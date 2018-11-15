#=======================================================================
#  This function copies the contents of "srcfile" to "tarfile". Also,
#     the function searches for lines containing characters similar to 
#     those in the list "input_strings", and replaces it with a line
#     containing the corresponding entries of the list "value" 
#
#      The purpose is to preserve "most" of the default inputs but
#  selectively modify certain others to run parametric sweeps/test cases
#=======================================================================

import re
import copy

def input_swap(srcfile,tarfile,input_strings,values):

#=======================================================================
#                 Begin executable code
#=======================================================================

      source                  = open(srcfile, 'r')
      target                  = open(tarfile, 'w')
      blanks                  = '      '

#=======================================================================
# print source file name, strings to identify and values to enter
#=======================================================================
      
#      print 'source file is ',source.name
#      print input_strings,values
      
#=======================================================================
#pad input string with wildcard characters to use regular expressions
#=======================================================================

      strng_list              = copy.deepcopy(input_strings)
      for i in range(0,len(input_strings)):
            strng_list[i] = '(.*)'+input_strings[i].lower()+'(.*)'
      
#=======================================================================
#check each line in source file
#=======================================================================

      for line in source:
            replace_flag      = False

            temp              = line.rstrip('\n').lower()

#=======================================================================
#loop over all key strings and identify whether it needs to be replaced
#=======================================================================

            for i in range(0,len(strng_list)):
                  if re.match(strng_list[i], line.lower()): #OLD LINE
#                  if re.search(strng_list[i], line.lower()): #NEW LINE
#                  if temp.startswith(strng_list[i]):
                        replace_flag = True
                        replace_id   = i                  

#=======================================================================
#initialize to default
#=======================================================================

            i = -1

#=======================================================================
#if the replacement event was triggered, inject custom line
#=======================================================================

            if replace_flag:
                  i = replace_id
                  print >> target, str(values[i]),blanks,blanks,blanks#,'TEST CASE',input_strings[i]
            else:
                  print >> target, line,        # lines in files have a default \n at the end, though invisible

#=======================================================================
#close the files        
#=======================================================================

      source.close()
      target.close()
      
#=======================================================================
#return statement
#=======================================================================

      return None
