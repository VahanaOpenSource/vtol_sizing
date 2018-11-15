# ===================================================================
#
# set fortran data from python/yaml inputs
#
# ===================================================================

import sys,os
import numpy as np
import yaml

import fea

class _airframeconfig:

   def initializeFEA(self):

      layout   = fea.layout
      material = fea.material
      current  = fea.currentvalues

      ndof_beam             = layout.ndof_beam
      self.nmatprop         = material.nmatprop
      current.maxdeflection = -1.0

# ===================================================================
# read from config.inp
# ===================================================================

      f = open("config.inp","r")
      
      line    = f.next()
      #
      line    = f.next()
      line    = line.strip()
      columns = line.split()

      nnodes = int(columns[0])
      nelem  = int(columns[1])
      #
      memberID  = np.zeros(nelem)
      #
      # node position
      #
      line    = f.next()
      for j in range(nnodes):
         line    = f.next()
      line    = f.next()
      #
      # element connectivity
      #
      for j in range(nelem):
         line    = f.next()
         line    = line.strip()
         columns = line.split()

         memberID[j]=int(columns[2])

      line    = f.next()
      line    = f.next()
      line    = line.strip()
      columns = line.split()

      nloading  = int(columns[0])
      maxnforce = int(columns[1])
      #
      # loop through the number of loading conditions
      #
      nforce  = np.zeros(nloading)
      for n in range(nloading):
         line    = f.next()
         line    = f.next()
         line    = line.strip()
         columns = line.split()

         nforce[n] = int(columns[0])
         #
         # force/moment boundary conditions
         #
         for j in range(int(nforce[n])):
            line    = f.next()
      #
      # read displacement condition
      #
      line    = f.next()
      line    = f.next()
      line    = line.strip()
      columns = line.split()

      ndisp  = int(columns[0])

      for j in range(ndisp):
         line    = f.next()
         line    = line.strip()
      #
      # read material properties
      #
      self.matprop = np.zeros(self.nmatprop)
      line    = f.next()
      for j in range(self.nmatprop):
         line    = f.next()
         line    = line.strip()
         columns = line.split()
         self.matprop[j] = float(columns[0])

      line    = f.next()
      line    = f.next()
      line    = line.strip()
      columns = line.split()
      self.cstype = str(columns[0])

      f.close()

# ===================================================================
# create Fortran link to FEA
# ===================================================================

      layout.nnodes    = nnodes
      layout.nelem     = nelem
      layout.nloading  = nloading
      layout.maxnforce = maxnforce
      layout.ndisp     = ndisp

      #
      # allocate fortran arrays
      #
      nrotor = self.all_dict['aircraft']['nrotor']
      
      layout.initializelayout()

      # copy data into fortran arrays
      for j in range(nelem):
         layout.memberid[j]=memberID[j]
      
      #
      # fea preprocess
      #
      fea.preprocess()

      current.initializecvals()


      #
      # initialize material database
      #
      material.initializematerial(nelem)
      #
      # copy to Fortran
      #
      for j in range(self.nmatprop):
         material.matprop[j][:] = self.matprop[j]

      material.cstype = self.cstype

      return None

# ###################################################################
# 
# resetMaterialProperties
#
# ###################################################################

   def resetMaterialProperties(self):

      material = fea.material

      for j in range(self.nmatprop):
         material.matprop[j][:] = self.matprop[j]

      material.cstype = self.cstype

      return None

# ###################################################################
# 
# setupFEA routine
#
# ###################################################################

   def setupFEA(self):

      layout   = fea.layout
      material = fea.material

      ndof_beam = layout.ndof_beam
      nmatprop  = material.nmatprop
      
# ===================================================================
# read from config.inp
# ===================================================================

      f = open("config.inp","r")
      
      line    = f.next()
      #
      line    = f.next()
      line    = line.strip()
      columns = line.split()

      nnodes = int(columns[0])
      nelem  = int(columns[1])
      
      rxyz      = np.zeros([nnodes,3])
      conn      = np.zeros([nelem,2])
      memberID  = np.zeros(nelem)
      #
      # node position
      #
      line    = f.next()
      for j in range(nnodes):
         line    = f.next()
         line    = line.strip()
         columns = line.split()

         rxyz[j][0] = columns[0]
         rxyz[j][1] = columns[1]
         rxyz[j][2] = columns[2]
         
      line    = f.next()
      #
      # element connectivity
      #
      for j in range(nelem):
         line    = f.next()
         line    = line.strip()
         columns = line.split()

         conn[j][0]  = int(columns[0])
         conn[j][1]  = int(columns[1])
         memberID[j] = int(columns[2])

      line    = f.next()
      line    = f.next()
      line    = line.strip()
      columns = line.split()

      nloading  = int(columns[0])
      maxnforce = int(columns[1])

      nforce         = np.zeros(nloading)
      extforce       = np.zeros([nloading,maxnforce])
      extforce_dofid = np.zeros([nloading,maxnforce])
      #
      # loop through the number of loading conditions
      #
      for n in range(nloading):
         line    = f.next()
         line    = f.next()
         line    = line.strip()
         columns = line.split()

         nforce[n] = int(columns[0])

         extforce[n][:] = 0.0
         #
         # force/moment boundary conditions
         #
         for j in range(int(nforce[n])):
            line    = f.next()
            line    = line.strip()
            columns = line.split()

            nodeid = int(columns[0])
            dofid  = int(columns[1])
            dummy  = float(columns[2])
            
            extforce[n,j] = extforce[n,j] + dummy

            extforce_dofid[n,j] = ndof_beam * (nodeid-1) + dofid
      #
      # read displacement condition
      #
      line    = f.next()
      line    = f.next()
      line    = line.strip()
      columns = line.split()

      ndisp  = int(columns[0])
      
      extdisp       = np.zeros(ndisp)
      extdisp_dofid = np.zeros(ndisp)

      for j in range(ndisp):
         line    = f.next()
         line    = line.strip()
         columns = line.split()

         nodeid     = int(columns[0])
         dofid      = int(columns[1])
         extdisp[j] = float(columns[2])

         extdisp_dofid[j] = ndof_beam * (nodeid-1) + dofid
      #
      # read material properties
      #
      matprop = np.zeros(nmatprop)
      line    = f.next()
      for j in range(nmatprop):
         line    = f.next()
         line    = line.strip()
         columns = line.split()
         matprop[j] = float(columns[0])

      line    = f.next()
      line    = f.next()
      line    = line.strip()
      columns = line.split()
      cstype = str(columns[0])

      f.close()

# ===================================================================
# create Fortran link to FEA
# ===================================================================

      layout.nnodes    = nnodes
      layout.nelem     = nelem
      layout.nloading  = nloading
      layout.maxnforce = maxnforce
      layout.ndisp     = ndisp
      #
      # allocate fortran arrays
      #
      nrotor = self.all_dict['aircraft']['nrotor']
      
      #print 'nelem ... ',nelem
      #layout.initializelayout()
      #return None
      #
      # copy data into Fortran arrays
      #
      for j in range(nnodes):
         layout.nodepos[0][j] = rxyz[j][0]
         layout.nodepos[1][j] = rxyz[j][1]
         layout.nodepos[2][j] = rxyz[j][2]
      
      for j in range(nelem):
         layout.conn[0][j]  = conn[j][0]
         layout.conn[1][j]  = conn[j][1]
         layout.memberid[j] = memberID[j]




      for n in range(nloading):

         layout.nforce[n] = int(nforce[n])

         layout.extforce[n][:] = 0.0
        
         for j in range(int(nforce[n])):
            
            layout.extforce[n,j]       = extforce[n,j]
            layout.extforce_dofid[n,j] = extforce_dofid[n,j]

      for j in range(ndisp):
         layout.extdisp[j]       = extdisp[j] 
         layout.extdisp_dofid[j] = extdisp_dofid[j]
      ##
      ## initialize material database
      ##
      #material.initializematerial(nelem)
      ##
      ## copy to Fortran
      ##
      #for j in range(nmatprop):
      #   material.matprop[j][:] = matprop[j]

      #material.cstype = cstype

# ===================================================================
# fea preprocess
# ===================================================================
      
      #fea.preprocess()
      
      return None

# ###################################################################
# 
# print fea output
#
# ###################################################################

   def writeFEAoutput(self):

      fea.writefeaoutputs()

      return None

# ===================================================================
# END OF FILE
# ===================================================================
