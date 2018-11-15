# ===================================================================
#
# set fortran data from python/yaml inputs
#
# ===================================================================

import sys,os
import numpy as np
import yaml

import fea
layout                  = fea.layout
material                = fea.material
current                 = fea.currentvalues
ndof_beam               = layout.ndof_beam

from   mpi4py import MPI
comm  = MPI.COMM_WORLD
rank  = comm.Get_rank()

class _airframeconfig:

#====================================================================
# deallocate all memory
#====================================================================

   def FEA_clean(self):
      print('deallocating fortran variables for FEA')
      fea.deallocate_variables()

#====================================================================
# function to read configuration and initialize FEA for airframe
# this function is run ONCE when HYDRA is initialized
#====================================================================

   def initializeFEA(self):
      
      self.FEA_clean()
      
#====================================================================
# set pointers to fortran memory
#====================================================================

      self.nmatprop           = material.nmatprop
      current.maxdeflection   = -1.0
      iodata                  = fea.fea_types.iodata

      comm.Barrier()
#====================================================================
# read details from config.inp
#====================================================================
      
      with open("config.inp","r") as f:

#====================================================================
# skip two lines
#====================================================================
      
         line     = next(f)
         line     = next(f)

#====================================================================
# third line has # of nodes, and # of elements: these do not change! 
#====================================================================

         line     = line.strip()
         columns  = line.split()

         nnodes   = int(columns[0])
         nelem    = int(columns[1])

#====================================================================
# node positions
#====================================================================

         rxyz        = np.zeros((3,nnodes))

         line        = next(f)
         for j in range(nnodes):
            line        = next(f).rstrip('\n').split()
            rxyz[0,j]   = float(line[0])
            rxyz[1,j]   = float(line[1])
            rxyz[2,j]   = float(line[2])

         line        = next(f)

#====================================================================
# element connectivity: "member" is a set of elements 
# remember list of member ID #s
#====================================================================

         conn        = np.zeros((nelem,2),dtype=int)
         memberID    = np.zeros(nelem,dtype=int)

         for j in range(nelem):
            line     = next(f)
            line     = line.strip()
            columns  = line.split()

            conn[j,0]   = int(columns[0])
            conn[j,1]   = int(columns[1])
            memberID[j] = int(columns[2])

#====================================================================
# next line is dummy, and the line after that contains 
# number of flight conditions, and max # loads applied per condition
#====================================================================

         line    = next(f)
         line    = next(f)
         line    = line.strip()
         columns = line.split()

         nloading  = int(columns[0])
         maxnforce = int(columns[1]);  

#====================================================================
# loop through the flight conditions and read external forcing values
#====================================================================

         forceid        = {}
         nforce         = np.zeros(nloading,dtype=int)      # number of loads/flight condition        
         for n in range(nloading):                          # loop over loads
            line        = next(f)                          # skip dummy line
            line        = next(f)
            line        = line.strip()
            columns     = line.split()

            nforce[n]   = int(columns[0])       # find number of forces in this flight condition

            nf          = nforce[n]
            Fdofs       = np.zeros((2,nf),dtype=int)

            for j in range(int(nforce[n])):                    # loop over # forces
               line        = next(f).rstrip('\n').split()     # skip it 
               Fdofs[0,j]  = int(line[0])                   # node                 
               Fdofs[1,j]  = int(line[1])                   # DOF 

            forceid[n]   = Fdofs 

#====================================================================
# read displacement boundary conditions
#====================================================================

         line        = next(f)                 # skip dummy row
         line        = next(f)
         line        = line.strip()
         columns     = line.split()
         ndisp       = int(columns[0])          # number of BC to apply

#remember boundary nodes and constrained DOF
         bc_dof      = np.zeros(ndisp,dtype=int)
         bc_node     = np.zeros(ndisp,dtype=int)

         for j in range(ndisp):              # skip actual nodes that are constrained
            line        = next(f)               # not relevant for allocation
            line        = line.rstrip('\n').strip().split()
            bc_node[j]  = int(line[0]) 
            bc_dof[j]   = int(line[1])

#====================================================================
# read properties for the construction material (single construction)
# store in "matprop" entry under the class
#====================================================================

         self.matprop         = np.zeros(self.nmatprop)     # initialize

         line = next(f)
         for j in range(self.nmatprop):                     # read 
            line              = next(f)
            line              = line.strip()
            columns           = line.split()
            self.matprop[j]   = float(columns[0])

#====================================================================
# read cross-section type 
#====================================================================

         line        = next(f)              # skip dummy line
         line        = next(f)              # value is here
         line        = line.strip()
         columns     = line.split()
         self.cstype = str(columns[0])

#====================================================================
# read reference lengths for nondiml model of airframe
#====================================================================

         line        = next(f)              # skip dummy line
         line        = next(f).rstrip('\n').split()

         ref_R       = float(line[0])  
         ref_b       = float(line[1])  

#==========================================================================
# read wing nodes that must be scaled by span (everything else scales by R)
#==========================================================================

         line        = next(f)              # skip dummy line
         line        = next(f).rstrip('\n')

         n_orignodes = int(line)
         orig_nodes  = iodata.orig_nodes        # from fortran memory
         for j in range(n_orignodes):
            orig_nodes[j]  = next(f).rstrip('\n')

#==========================================================================
# read wing nodes that must be interpolated between other points
#==========================================================================

         line              = next(f)
         line              = next(f).rstrip('\n')
         n_wingnodes       = int(line)
         iw_nodes          = iodata.iw_nodes
         wtw_nodes         = iodata.wtw_nodes

         for j in range(n_wingnodes):
            line           = next(f).rstrip('\n').split()
            iw_nodes[0,j]  = int(line[0])          # actual wing node
            iw_nodes[1,j]  = int(line[1])          # left lookup node
            iw_nodes[2,j]  = int(line[2])          # right lookup node

            wtw_nodes[0,j] = float(line[3])        # left interp weight
            wtw_nodes[1,j] = float(line[4])        # right interp weight

#====================================================================
# read rotor nodes count and actual # of rotor nodes
#====================================================================
      
         line              = next(f)        # skip dummy line
         line              = next(f).rstrip('\n').split()
         nrotors           = int(line[0])
         iodata.id_thrust  = int(line[1])
         iodata.id_torque  = int(line[2])

         iodata.nrotors = nrotors
         line              = next(f)        # skip dummy line
         for i in range(nrotors):
            line                 = next(f).rstrip('\n').split()
            iodata.id_rotor[i]   = int(line[0])
            iodata.dirn[i]       = float(line[1])

#====================================================================
# read wing related parameters (# wings, # elements per wing and 
#element IDs in that wing
#====================================================================
      
         line              = next(f)                 # skip dummy line
         line              = next(f).rstrip('\n')    # number of wings
         nwings            = int(line)

         iodata.nwings     = nwings          # store in data structure

#============================
# loop over wings
#============================

         for i in range(nwings):
            line                 = next(f)        # skip dummy line
            line                 = next(f).rstrip('\n')

#============================
# set # wing elements 
#============================

            nwelem               = int(line)
            iodata.nwing_elem[i] = nwelem 

#============================
# read wing element IDs
#============================
            
            line                 = next(f)        # skip dummy line
            for j in range(nwelem):
               line                = next(f).rstrip('\n')
               iodata.iw_elem[j,i] = int(line)

#====================================================================
# store details in fortran memory
#====================================================================

      layout.nnodes    = nnodes
      layout.nelem     = nelem
      layout.nloading  = nloading
      layout.maxnforce = maxnforce
      layout.ndisp     = ndisp

#====================================================================
# allocate fortran memory based on array dimensions in input file
#====================================================================

      layout.initializelayout()

#====================================================================
# copy data into fortran arrays
#====================================================================

      for j in range(nelem):
         layout.memberid[j]=memberID[j]
      
#====================================================================
# RUN fea preprocessor and initialize memory
#====================================================================

      fea.preprocess()
      current.initializecvals()

#====================================================================
# initialize material database
#====================================================================

      material.initializematerial(nelem)

#====================================================================
# copy to Fortran
#====================================================================

      for j in range(self.nmatprop):
         material.matprop[j,:] = self.matprop[j]

      material.cstype = self.cstype

#====================================================================
# dictionary with configuration details
#====================================================================

      self.fea_dict    = {'nnodes': nnodes,   'nelem'    : nelem, 
                            'rxyz': rxyz,     'memberID' : memberID ,
                        'nloading': nloading, 'maxnforce': maxnforce, 
                          'nforce': nforce,       'ndisp': ndisp}

#====================================================================
# set reference radius and single-number values
#====================================================================

      iodata.ref_r         = ref_R        # reference radius 
      iodata.ref_b         = ref_b        # reference span 

      iodata.norig_nodes   = n_orignodes  # number of nodes to scale by span
      iodata.nwing_nodes   = n_wingnodes  # number of nodes to interpolate

#====================================================================
# transfer normalized coordinates to fortran (one-time cost)
#====================================================================

      fea.transfer_norm_xyz(nnodes, rxyz)

#====================================================================
# connections between nodes and member IDs
#====================================================================
   
      connf      = layout.conn
      memid      = layout.memberid
      for j in range(nelem):
         connf[0,j]  = conn[j,0]
         connf[1,j]  = conn[j,1]
         memid[j]    = memberID[j]   

#====================================================================
# move displacement boundary conditions to fortran memory
#====================================================================

      ndof_beam      = layout.ndof_beam
      bc_id          = layout.extdisp_dofid
      ext_disp       = layout.extdisp
      for j in range(ndisp):
         nodeid      = bc_node[j]
         dofid       = bc_dof[j]
         bc_id[j]    = ndof_beam * (nodeid-1) + dofid
         ext_disp[j] = 0.0e0
#         print ('constrained node %d with DOF %d' % (nodeid,dofid))

#====================================================================
# set force boundary condition info
#====================================================================

      force_dof      = layout.extforce_dofid
#      print type(force_dof)
      for n in range(nloading):
         nloads            = int(nforce[n])
         layout.nforce[n]  = nloads 
         for j in range(nloads):
            nodeid         = forceid[n][0,j]
            dofid          = forceid[n][1,j]

            force_dof[n,j] = ndof_beam * (nodeid-1) + dofid

#====================================================================
# set boundary condition indices for forcing and displacement
#====================================================================
   
      fea.assign_indices()

#====================================================================
# end of operations
#====================================================================

      return None

# ###################################################################
# 
# resetMaterialProperties for all elements!!!
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
# print fea output
#
# ###################################################################

   def writeFEAoutput(self):

      fea.writefeaoutputs()

      return None

#====================================================================
# method to transfer reference data to fortran memory 
#====================================================================

#   def 
#====================================================================
# function to scale finite element model based on loads and 
# other design parameters
#====================================================================

   def scale_dimensions(self, radius):


      # self.fea_dict    = {'nnodes': nnodes,   'nelem'    : nelem, 
      #                       'rxyz': rxyz,     'memberID' : memberID ,
      #                   'nloading': nloading, 'maxnforce': maxnforce, 
      #                     'nforce': nforce,       'ndisp': ndisp}


      return None
# ===================================================================
# END OF FILE
# ===================================================================