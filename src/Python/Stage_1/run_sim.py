#=======================================================================
#       user-defined python functions by Ananth Sridharan
#  These functions are used to post-process outputs of ORCA/PrasadUM
#=======================================================================

#=======================================================================
#            Import required modules for operations
#=======================================================================

import os
import subprocess
import sys
import shutil
#import pywake

#=======================================================================
#       Function to copy all files from one folder to another
#=======================================================================

def copy_files(source,target):

#define a default 'ok' status for exit unless input folder doesnt exist
      errcode     = ''

#if directory exists, flush all files
      src_exists  = os.path.isdir(source)
      tar_exists  = os.path.isdir(target)
      
#check if directories exist
      if not src_exists:
            errcode = errcode + 'source folder not found:' + source
      if not tar_exists:
            errcode = errcode + 'target folder not found:' + target
      if src_exists and tar_exists:
            errcode = 'ok'

#copy all files from input path to temp folder
      if errcode == 'ok':

#this is linux system command to copy files recursively
            os.system('cp -r '+source+'* '+target)

#this is python command to copy files recursively
#remove target destination : copytree requires that destination not exist
#            shutil.rmtree(target)
#            shutil.copytree(source,target)

      else:
            sys.exit(errcode)
            
      return errcode
      
#=======================================================================
#    Function to create a directory in current folder
#=======================================================================

def create_dir(temp_folder):

#get current directory
      current_dir = os.getcwd()
      
#define a temp directory location to backup inputs
      temp_dir    = current_dir + '/' + temp_folder

#check if directory exists
      dir_exists  = os.path.isdir(temp_dir)

#remove the directory and recreate it [this is a problem with rmtree: it
#deletes the directory along with its files, hence we have to recreate  
      if not dir_exists:
            os.mkdir(temp_dir)
            
      return None
      
#=======================================================================
#           Function to execute a compiled fortran/C program
#=======================================================================

def run_code(exec_path, exec_name):
      
#=======================================================================
#check if directory exists      
#=======================================================================

      dir_exists  = os.path.isdir(exec_path)
      errcode     = 'ok'

      if(dir_exists):

#=======================================================================
#check if executable exists
#=======================================================================

            file_exists = os.path.isfile(exec_path+exec_name)

#=======================================================================
#if file exists, "cd" to that location, execute code and "cd" back
#=======================================================================

            if(file_exists):
                  #print 'changing current working directory to ',exec_path

#=======================================================================
#fortran version
#=======================================================================

                  current_dir = os.getcwd()
                  os.chdir(exec_path)
                  #print os.getcwd()

                  #print 'calling Fortran code'
                  subprocess.call('./'+exec_name)

#=======================================================================
#python module version
#=======================================================================

#                  print 'calling pysadum'
#                  pywake.pysadum()
                  
#=======================================================================

                  #print 'reverting to original working directory',current_dir
                  os.chdir(current_dir)
            else:
                  errcode = 'directory exists but executable not found'
                  print('ERROR: directory exists; program not found')
      else:
            errcode = 'directory not found'
            print ('ERROR: directory not found')
      
#end of operations      
      return errcode
