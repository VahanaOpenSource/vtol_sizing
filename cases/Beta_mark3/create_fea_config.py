import os,sys,copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import numpy as np
from matplotlib import cm, colors, patches
#
radius = 1.0 

iplot = False

from   mpi4py import MPI
comm  = MPI.COMM_WORLD
rank  = comm.Get_rank()

#==============================================================================
# function to create user airframe model for FEA and weight calculations for 
# load-bearing structural members 
#==============================================================================
# two inputs: parameters of the aircraft (params) and a dictionary with 
# relevant numbers passed to actual model generator from this program, which 
# is a geometry generator
#
#
#   Aircraft parameters are :
# radius,     # rotor radius [m]
# wingspan,   # wing span [m]
# thrust,     # thrust in [N] (per rotor)
# hvr_torque, # torque in [Nm] (per rotor)
# ff_torque,  # torque in [Nm] (per rotor)
# drag,       # drag in forward flight [N] (per rotor) 
# winglift):  # lift on each wing (per wing)
#==============================================================================

def createAirframeConfig(params, print_flag):

#==============================================================================
# unpack parameters from sizing (used to scale model + apply loads)
#==============================================================================

   radius     = params['R']
   wingspan   = params['b']
   thrust     = params['T']
   hvr_torque = params['Qh']
   ff_torque  = params['Qc']
   drag       = params['D']
   winglift   = params['L']

#==============================================================================
# screen print
#==============================================================================
  
   if False:
      print( '\n')
      print( '----------------------------------')
      print( '    IN CREATE AIRFRAME CONFIG     ')
      print( '----------------------------------')
      print( ' radius     [m]    %9.3f' % ( radius     ))
      print( ' wingspan   [m]    %9.3f' % ( wingspan   ))
      print( ' thrust     [N]    %9.3f' % ( thrust     ))
      print( ' hvr_torque [Nm]   %9.3f' % ( hvr_torque ))
      print( ' ff_torque  [Nm]   %9.3f' % ( ff_torque  ))
      print( ' drag       [N]    %9.3f' % ( drag       ))
      print( ' winglift   [N]    %9.3f' % ( winglift   ))
      print( '----------------------------------')
      print( '\n')

#==============================================================================
# begin executable code
#==============================================================================

   BASE      = 1
   nbeamelem = 3
   ihover    = True
   iforward  = True

   #get normalized half span
   span = wingspan * 0.5 * 0.1905/radius

   # scale the rotor mount point location with the wing span
   # performed to prevent excessive bending moments at the root
   rmnt = 0.2329 

#==============================================================================
#changing the wing span to be at least rotor diameter + clearance
#==============================================================================

   if (span < rmnt):                  # nondiml span > 0.23 R
      span            = 1.1 * rmnt 
      wingspan        = radius*span / (0.5*0.1905)
      params['b']     = wingspan
#      print 'warning!!!! span is changing'
#==============================================================================
#
#==============================================================================

   #
   # position of the points
   #
   rxyz0 = [ [  -span,    rmnt,  0.0  ],#  1 wing tip top left  
             [  -rmnt,    rmnt,  0.0  ],#  2 rotor mount
             [   rmnt,    rmnt,  0.0  ],#  3 rotor mount
             [   span,    rmnt,  0.0  ],#  4 wing tip top right
             [ 0.0   ,  0.0   ,  0.0  ],#  5 center of vehicle
             [  -span,   -rmnt,  0.0  ],#  6 wing tip bottom left 
             [  -rmnt,   -rmnt,  0.0  ],#  7 rotor mount
             [   rmnt,   -rmnt,  0.0  ],#  8 rotor mount
             [   span,   -rmnt,  0.0  ],#  9 wing tip bottom right
             [-0.0900,  0.0900,  0.0  ],# 10
             [ 0.0900,  0.0900,  0.0  ],# 11
             [-0.0900, -0.0900,  0.0  ],# 12
             [ 0.0900, -0.0900,  0.0  ],# 13 
             [-0.0900,  0.0900, -0.21 ],# 14
             [ 0.0900,  0.0900, -0.21 ],# 15
             [-0.0900, -0.0900, -0.21 ],# 16
             [ 0.0900, -0.0900, -0.21 ],# 17
             [  -rmnt,    rmnt,  0.14 ],# 18 rotor hub 1
             [   rmnt,    rmnt,  0.14 ],# 19 rotor hub 2
             [  -rmnt,   -rmnt,  0.14 ],# 20 rotor hub 3
             [   rmnt,   -rmnt,  0.14 ] # 21 rotor hub 4
           ]

   # divide by radius of aierz (0.1905 m) and multiply by 
   # actual radius
   rxyz0 = np.multiply(rxyz0,radius/0.1905)
   #
   # element connectivity
   #
   conn0 = np.array([[  1,  2, 2],  #  1    => wing
                     [  2,  3, 0],  #  2
                     [  3,  4, 2],  #  3    => wing 
                     [  6,  7, 2],  #  4    => wing
                     [  7,  8, 0],  #  5
                     [  8,  9, 2],  #  6    => wing
                     [  2, 10, 2],  #  7
                     [ 10,  5, 2],  #  8
                     [  5, 13, 2],  #  9
                     [ 13,  8, 2],  # 10
                     [  3, 11, 2],  # 11
                     [ 11,  5, 2],  # 12
                     [  5, 12, 2],  # 13
                     [ 12,  7, 2],  # 14
                     [ 10, 14, 0],  # 15
                     [ 11, 15, 0],  # 16
                     [ 12, 16, 0],  # 17
                     [ 13, 17, 0],  # 18
                     [ 18,  2, 2],  # 19
                     [ 19,  3, 2],  # 20
                     [ 20,  7, 2],  # 21
                     [ 21,  8, 2],  # 22
                     [ 14, 15, 2],  # 23
                     [ 15, 17, 2],  # 24
                     [ 17, 16, 2],  # 25
                     [ 16, 14, 2]]) # 26

#=======================================================================
#which nodes are on wing tips
#=======================================================================

   wing_id    = [1,2]
   nwings     = len(wing_id)
   orig_nodes = [[1,4],[6,9]]            # original nodes + intermediate points

#=======================================================================
# other nodes that are also on wing (rotor mounts)
#=======================================================================

   xn_nodes     = [[2,3], [7,8]]

   rotor_nodes  = [18, 19, 20, 21]
   dirn         = [1,  -1, -1,  1]
   thrust_id    = 3
   torque_id    = 6

#=======================================================================
# additional nodes generated to be on the wing
#=======================================================================

   wing_nodes = []                      # only original nodes

#=======================================================================
# list of beam elements that are mapped to a fixed wing
#=======================================================================

   wing_elem  = {}

#=======================================================================
# interpolation weights and IDs of source nodes for generated elements
#=======================================================================

   node_wts   = []                      # interpolation weights 
   node_ids   = []                      # interpolation indices

#=======================================================================
# types of different loadings to design for (hover, forward flight)
#=======================================================================

   if(ihover==True and iforward==False):                  
      nloadings = 1
   elif(ihover==False and iforward==True):
      nloadings = 1
   elif(ihover==True and iforward==True):
      nloadings = 2
   else:
      print ('Incorrect loading condition...',ihover,iforward)
      sys.exit(1)

#=======================================================================
# force/moment loading on the structure
# array with 3 columns per entry (nodeid, dofid, value of force/moment)
#forceid  = np.array([[1,0,0]])
#=======================================================================

   allforce = []
 
   forceOut = {}

#=======================================================================
# hover
#=======================================================================

   if(ihover):
      forceid = np.array([[18, 3,  thrust],
                          [18, 6,  hvr_torque],
                          [19, 3,  thrust],
                          [19, 6, -hvr_torque],
                          [20, 3,  thrust],
                          [20, 6, -hvr_torque],
                          [21, 3,  thrust],
                          [21, 6,  hvr_torque]])
      
#=======================================================================
# append loading information
#=======================================================================

      allforce.append(forceid)
   
#=======================================================================
# forward flight
#=======================================================================

   if(iforward):
      forceid = np.array([[18, 3, drag], 
                          [18, 6, ff_torque], 
                          [19, 3, drag], 
                          [19, 6,-ff_torque], 
                          [20, 3, drag], 
                          [20, 6,-ff_torque], 
                          [21, 3, drag], 
                          [21, 6, ff_torque]])
      
#=======================================================================
# append loading information
#=======================================================================

      allforce.append(forceid)

#=======================================================================
# distributed loads (member ID, dofid, load (N/m))
#=======================================================================

   alldistload = []
   
#=======================================================================
# hover
#=======================================================================

   if(ihover):
      distloadid = np.array([])
      
#=======================================================================
# append distributed load
#=======================================================================

      alldistload.append(distloadid)
   
#=======================================================================
# forward flight
#=======================================================================

   if(iforward):   
      distloadid = np.array([
                             [[1,2,3],2,winglift], # simulating lift from wing
                             [[4,5,6],2,winglift]
                           ])
      
#=======================================================================
# append distributed load
#=======================================================================

      alldistload.append(distloadid)

#=======================================================================
# array with 3 columns per entry (nodeid, dofid, value of displacement)
#=======================================================================

   dispid = np.array([[5, 1, 0.0],
                      [5, 2, 0.0],
                      [5, 3, 0.0],
                      [5, 4, 0.0],
                      [5, 5, 0.0],
                      [5, 6, 0.0]])

#=======================================================================
#=======================================================================
# the rest of the code takes the above inputs and creates
# the input file necessary for the fea code
# no user intervention required
#=======================================================================
#=======================================================================
   
# ===================================================================
# given the above rxyz and conn
# create the node and conn array assuming nbeamelem
# ===================================================================
   
   conn0   -= BASE # convert to base 0 from base 1
   
   nnodes0  = len(rxyz0)
   nelem0   = len(conn0)
   
   nnodes   = nnodes0 + nelem0*(nbeamelem-1)
   nelem    = nelem0 * nbeamelem
   
   rxyz     = np.zeros([3,nnodes])
   conn     = np.zeros([nelem,2])
   memberid = np.zeros(nelem)

   for i in wing_id:
     wing_elem[i] = []

#=======================================================================
# copy positions
#=======================================================================

   for j in range(nnodes0):
      rxyz[:,j] = rxyz0[j][:]

#=======================================================================
# create member to node array
#=======================================================================

   maxconn      = 10
   nmembers     = len(conn0)
   member2node  = np.empty([nmembers,maxconn])
   member2node.fill(-1)
   nmember2node = np.zeros(nmembers)

#=======================================================================
# create the position and conn arrays 
#=======================================================================

   inode        = nnodes0
   ielem        = 0
  
   for j in range(nelem0):

#==================================================================
# find nodes at the ends of this member
#==================================================================

      n1        = conn0[j][0]
      n2        = conn0[j][1]

#==================================================================
# flag to indicate whether the end nodes are wing tips or not
#==================================================================

      interp_node   = False 
      iw            = 0
      for iwing in wing_id:
        test1       = (n1+1) in orig_nodes[iwing-1]   # true if n1 is wing tip
        test2       = (n2+1) in orig_nodes[iwing-1]   # true if n2 is wing tip
        test3       = (n1+1) in xn_nodes[iwing-1]     # true if n1 is rotor mount
        test4       = (n2+1) in xn_nodes[iwing-1]     # true if n2 is rotor mount

#==================================================================
# if only one node in a member is at the wing tip, then the other 
# is an intermediate node. We test for this condition using XOR 
# the Python operator for XOR (exclusive OR) is the carat sign (^)
# if the interp_node flag is triggered, then intermediate nodes 
# between the wing tip and other locations may be generated 
# by the code. Those nodes are stored separately and their positions
# are interpolated in Fortran based on the end conditions.
#==================================================================

        if (test1 ^ test2):
            iw          = iwing       # wing id for wing tip
            interp_node = True

#==================================================================
# flag to test whether the end nodes are on the wing 
#  if the end nodes are on the wing (rotor mounts or wing tips)
# then any elements generated for this "member" are wing elements
# for which distributed lift is applied in the FE model. We want 
# to remember a list of all wing elements and their associated 
# wing identification numbers, wing_id.
#
# Coding wise, we test whether the end nodes (n1 and n2) for this
# "member" correspond to wing tips (test1 for n1, and test2 for n2)
# or rotor mounts (test3 for n1, and test4 for n2)
# If both conditions are satisfied, then its a wing member and all 
# elements in this member are wing elements
#==================================================================

#==================================================================
# check if the end node "n1" is at wing tip or intersection 
# i.e. whether n1 is in orig_nodes or xn_nodes
# if either is true, n1 is on a wing 
# similarly for n2
#==================================================================

        test3   = test1 or test3      # true if n1 is on a lifting wing 
        test4   = test2 or test4      # true if n2 is on a lifting wing

#=======================================================================
# check if both nodes are on wings; if so, its a wing member
#=======================================================================

        is_wing = test3 and test4     # if wing member, is_wing is true

        if is_wing:
          iw      = iwing
          break
#==================================================================
# end locations of the member
#==================================================================

      r1        = rxyz0[n1][:]
      r2        = rxyz0[n2][:]
   
      if(conn0[j][2] == -1):        # no extra subdivisions needed
         nn = nbeamelem - 1
      else:
         nn = conn0[j][2]
   
#==================================================================
# Case with 1 beam element
#==================================================================

      if(nn==0):
         
         memid           = j
   
         conn[ielem][0]  = int(n1)        # left end  node
         conn[ielem][1]  = int(n2)        # right end node
         memberid[ielem] = memid

#==================================================================
# add conn[ielem][0] to member2node
#==================================================================

         if (conn[ielem][0] in member2node[memid][:]) == False:
            member2node[memid][nmember2node[memid]] = conn[ielem][0]
            nmember2node[memid] += 1
   
         # add conn[ielem][1] to member2node
         if (conn[ielem][1] in member2node[memid][:]) == False:
            member2node[memid][nmember2node[memid]] = conn[ielem][1]
            nmember2node[memid] += 1

#==================================================================
# test for whether its a wing element; if so, remember ID number
#==================================================================
   
         if is_wing:
            wing_elem[iw].append(ielem+1)   # fortran indexing

#==================================================================
# increment element counter
#==================================================================
         
         ielem += 1
   

#==================================================================
# more than 1 beam element
# add intermediate nodes to the model
#==================================================================

      else:

#==================================================================
# generate the points (0 to nn)
#==================================================================

         for k in range(nn):
            rxyz[0,inode] = r1[0] + (r2[0]-r1[0])/(nn+1) * (k+1)
            rxyz[1,inode] = r1[1] + (r2[1]-r1[1])/(nn+1) * (k+1)
            rxyz[2,inode] = r1[2] + (r2[2]-r1[2])/(nn+1) * (k+1)

#==================================================================
# node indices
#==================================================================
            
            nn1   = inode - 1
            nn2   = inode + 0
    
#==================================================================
#interpolation weights to determine intermediate point location
#==================================================================

            w2    = 1.0/(nn+1.0)*float(k+1)
            w1    = 1.0 - w2
            if(k==0):
               nn1 = n1

#==================================================================
# remember connection information
#==================================================================
   
            conn[ielem][0]  = int(nn1)
            conn[ielem][1]  = int(nn2)
            memberid[ielem] = j
            
            memid           = int(memberid[ielem])
   
#==================================================================
# add conn[ielem][0] to member2node
#==================================================================

            if (conn[ielem][0] in member2node[memid][:]) == False:
               member2node[memid][int(nmember2node[memid])] = conn[ielem][0]
               nmember2node[memid] += 1
   
#==================================================================
# add conn[ielem][1] to member2node
#==================================================================

            if (conn[ielem][1] in member2node[memid][:]) == False:
               member2node[memid][int(nmember2node[memid])] = conn[ielem][1]
               nmember2node[memid] += 1

#==================================================================
# tag rows which correspond to fixed wing members
# clearly, we're generating a new node so add it to wing_nodes
# store interpolation indices and weights too..
#==================================================================

            if interp_node:         
               wing_nodes.append(inode+1)
               node_wts.append([w1,w2])
               node_ids.append([n1+1,n2+1])

#==================================================================
# check if its a wing element and remember
#==================================================================

            if is_wing:
              wing_elem[iw].append(ielem+1)
              #print nn1+1,nn2+1

#              print 'adding element to wing',n1+1,n2+1,'element ',ielem,iw
#            else:
#                print 'not a wing element',n1+1,n2+1,'element',ielem

#==================================================================
# increment element counter
#==================================================================

            ielem += 1
            inode += 1
   
#==================================================================
# remember connection between elements
#==================================================================
   
         conn[ielem][0]  = int(inode-1)
         conn[ielem][1]  = int(n2)
         memberid[ielem] = j
   
         memid = int(memberid[ielem])
   
#==================================================================
# add conn[ielem][0] to member2node
#==================================================================

         if (conn[ielem][0] in member2node[memid][:]) == False:
            member2node[memid][int(nmember2node[memid])] = conn[ielem][0]
            nmember2node[memid] += 1
   
         # add conn[ielem][1] to member2node
         if (conn[ielem][1] in member2node[memid][:]) == False:
            member2node[memid][int(nmember2node[memid])] = conn[ielem][1]
            nmember2node[memid] += 1

#==================================================================
# check if its a wing element and remember
#==================================================================

         if is_wing:
            wing_elem[iw].append(ielem+1)
            #print inode,n2+1

#            print 'adding element to wing',n1+1,n2+1,'element ',ielem,iw
#         else:
#            print 'not a wing element',n1+1,n2+1,'element',ielem

#==================================================================
# add one more to element counter
#==================================================================
   
         ielem += 1

#==================================================================
# end of loop over prescribed connections 
#==================================================================

#==================================================================
# loop over the number of loading conditins
#==================================================================

   maxnforce      = 0
   nforces        = []
     
   for n in range(nloadings):
   
#      print 'n is ',n 
      forceid    = allforce[n]
      distloadid = alldistload[n]
   
      #
      # logic to find the nodes to apply the distributed loads
      # 
      # loop through the number of distributed loads
      ndistload = len(distloadid)

#      print 'n is ',n,' at the top-ish'

      for j in range(ndistload):
      
         member_length = 0.0
         nmember_load  = len(distloadid[j][0])
         nodelist      = np.zeros(0,dtype=int)

#         print 'n is ',n,' at the top-middle'
      
         # loop through the members for each distributed load
         for k in range(nmember_load):
            memid = int(distloadid[j][0][k]) - BASE
      
            # add length to member
            n1 = conn0[memid][0]
            n2 = conn0[memid][1]
            member_length += np.sqrt(  (rxyz[0,n1]-rxyz[0,n2])*(rxyz[0,n1]-rxyz[0,n2])  \
                                     + (rxyz[1,n1]-rxyz[1,n2])*(rxyz[1,n1]-rxyz[1,n2])  \
                                     + (rxyz[2,n1]-rxyz[2,n2])*(rxyz[2,n1]-rxyz[2,n2]))
      
            # find a list of all nodes belonging to the members
            for node in member2node[memid][:]:
               node = int(node)
               if (node >= 0 and ((node in nodelist)==False)):
                  nodelist = np.append(nodelist,node)
      
         # apply the distributed load (defined as force/length)
         # to the nodes in nodelist. Must multiply by the length
         # of the member to ensure compatible dimensions
         #print 'nodelist',nodelist
         #for node in nodelist:
         #   print 'node',node,'position',rxyz[node,:]
         #
         # currently nodes are arranged in ascending order of 'x'
         # which is convienient..
         #
         nnodes = len(nodelist)
         midpoints = np.zeros([nnodes-1,3])
         for kk in range(nnodes-1):
            n1 = nodelist[kk]
            n2 = nodelist[kk+1]
            midpoints[kk,:] = 0.5*(rxyz[:,n1] + rxyz[:,n2])
         
         integ_limits = np.zeros([nnodes,2])

         #print '\n\nmidpoints',midpoints

         #
         # find limits of integration
         #
         for kk in range(nnodes):

            if kk==0:
               n1 = nodelist[kk]
               xstart = rxyz[:,n1]
            else:
               xstart = midpoints[kk-1,:]
            #
            if kk == nnodes-1:
               n1 = nodelist[kk]
               xend = rxyz[:,n2]
            else:
               xend   = midpoints[kk,:]
            
            integ_limits[kk,0] = xstart[0]
            integ_limits[kk,1] = xend[0]
         #print '\n\ninteglimits',integ_limits
         #
         # find area under curve
         #
#         print 'n is ',n,' at the middle'
         npts = 10
         totarea = 0
         area_node = np.zeros(nnodes)
         for kk in range(nnodes):
            xstart = integ_limits[kk,0]
            xend   = integ_limits[kk,1]
            
            deltax = (xend-xstart)/npts
            
            # rectangular rule
            area = 0
            for n1 in range(npts):
               x1 = xstart + n1 * deltax
               x2 = x1 + deltax
               xmid = 0.5*(x1+x2)

#original: put EDL
#               area += np.sqrt(1.0 - (2*xmid/wingspan)**2) * deltax

#change: put UDL (AS)               
               area += deltax*np.pi/4.0
            area_node[kk] = area
            totarea += area
         #
         # find L0
         #
         L0 = 4*winglift/(wingspan*np.pi)     
 
         #
         dofid            = distloadid[j][1]
         force_per_length = distloadid[j][2]
               
         # elliptical distribution of force
         kk = -1
         temp =0
         for node in nodelist:
            kk+=1
            force_in_node = L0 * area_node[kk]
            forceid = np.vstack((forceid,(int(node)+BASE,int(dofid),force_in_node)))
            temp += force_in_node

#      print 'n is ',n,' at the bottom'
      forceOut[str(n)] = forceid

      nforce      = len(forceid)
      nforces.append(nforce) 
      maxnforce   = max(maxnforce,nforce)
   
   nforces = np.asarray(nforces)

# ===================================================================
# 
# ===================================================================

# ===================================================================
# write out to file
# ===================================================================
   
   #convert back to base 1
   conn     += BASE
   memberid += BASE
   
   nnodes         = inode
   nelem          = ielem
   #     
   ndisp = len(dispid)

#   quit('all ok?')
#==================================================================================
# store values in a dictionaries
#==================================================================================

   if(print_flag == 1 and rank == 0):
     with open('config.inp','w') as f:
     
       f.write("# nnodes, nelem\n")
       f.write(" %d %d\n" % (nnodes,nelem) )
       
       f.write("# node positions\n")
       for j in range(nnodes):
          f.write("%20.12e %20.12e %20.12e\n" % (rxyz[0,j], rxyz[1,j],rxyz[2,j]) )
       
       f.write("# element connectivity\n")
       for j in range(nelem):
          f.write("%d %d %d\n" % (conn[j][0], conn[j][1], memberid[j] ))
       
       f.write("# number of loading conditions, maxnforce (for allocation)\n")
       f.write("%d %d\n" % (nloadings,maxnforce))
       
       for n in range(nloadings):
          
          forceid = forceOut[str(n)]
       
          nforce = len(forceid)
          f.write("# force/moment boundary condition (nforce \\ nodeid dof id || 1-Fx 2-Fy 3-Fz 4-Mx 5-My 6-Mz)\n")
          f.write(" %d\n" % nforce)
          for j in range(nforce):
             f.write(" %d %d %20.12e\n" % (forceid[j][0],forceid[j][1],forceid[j][2]))
            
       f.write("# displacement boundary condition (ndisp \\ nodeid, dof_id, displacement || 1-u 2-v 3-w 4-thx 5-thy 6-thz)\n")
       f.write(" %d\n" % ndisp)
       for j in range(ndisp):
          f.write(" %d %d %e\n" % (dispid[j][0],dispid[j][1],dispid[j][2]))
       
       
       #
       # write material properties
       #
       f.write("# element properties (Al 6063-T6)\n")
       f.write("68.9e9     #  1 - elastic modulus\n")
       f.write("0.33       #  2 - poissons ratio\n")
       f.write("2700.0     #  3 - density of material\n")
       f.write("3.3333e-5  #  4 - cross-sectional area\n")
       f.write("9.2592e-11 #  5 - Iy\n")
       f.write("9.2592e-11 #  6 - Iz\n")
       f.write("1.8518e-10 #  7 - J (polar moment)\n")
       f.write("0.0        #  8 - ky (shear factor)\n") 
       f.write("0.0        #  9 - kz (shear factor)\n")
       f.write("2.8867e-3  # 10 - max y dimension\n")
       f.write("2.8867e-3  # 11 - max z dimension\n")
       f.write("276.e6     # 12 - Yield stress\n")

#==================================================================================
# write cs type
#==================================================================================

       f.write("#   cross section type\n")
       #f.write("solidsquare\n")
       f.write("hollowcircle\n")
       f.write("#  reference radius, span\n")
       f.write("%20.12e %20.12e\n" % (radius,wingspan))

#==================================================================================
#write original node ids to scale by wing span
#==================================================================================

       f.write("#  fixed wing node ids: to scale by span\n")
       orig_nodes = np.asarray(orig_nodes).flatten()
       N_orig     = np.shape(orig_nodes)[0]

       f.write("%4d\n" % len(orig_nodes)) 
       for j in orig_nodes: 
          f.write("%6d \n" % j)

#==================================================================================
#rest of nodes are to be interpolated between radius-scaled and span-scaled
#order is : node id, 2 nodes its interpolated from, interpolation weights
#for original wing nodes (wing tips), its referenced to its own nondiml position
#==================================================================================

       f.write("#  interpolation node ids: to correct behavior in future\n")
       f.write("%4d\n" % len(wing_nodes)) 
       for ij,j in enumerate(wing_nodes):
          il  = node_ids[ij][0]
          ir  = node_ids[ij][1]
          wl  = node_wts[ij][0]
          wr  = node_wts[ij][1]
          f.write("%6d %6d %6d %20.12e %20.12e \n" % (j,il,ir,wl,wr))

#==================================================================================
#rest of nodes are to be interpolated between radius-scaled and span-scaled
#order is : node id, 2 nodes its interpolated from, interpolation weights
#for original wing nodes (wing tips), its referenced to its own nondiml position
#==================================================================================

       rotor_nodes   = np.asarray(rotor_nodes).flatten()
       f.write("#Number of rotors, DOF id for thrust and torque\n")
       f.write("%6d %6d %6d\n " % (len(rotor_nodes), thrust_id, torque_id))
       f.write("#Rotor nodes and direction of rotation\n")
       nhub       = len(rotor_nodes)
       for i in range(nhub):
          f.write("%6d %6.2f \n" % (rotor_nodes[i],float(dirn[i])))

#==================================================================================
#rest of nodes are to be interpolated between radius-scaled and span-scaled
#order is : node id, 2 nodes its interpolated from, interpolation weights
#for original wing nodes (wing tips), its referenced to its own nondiml position
#==================================================================================

       f.write("#Number of wings\n")
       nwing    = len(wing_id)
       f.write("%6d \n" % nwing)

       for j in wing_id:
          temp  = wing_elem[j]
          f.write("# number of elements in wing %d \n" % j)
          f.write("%6d \n" % len(temp))
          f.write("# element IDs for wing %d \n" % j)
          for i in temp:
            f.write("%6d \n" % i)

       f.close()
       print ('wrote config.inp file for airframe structure')
     if iplot:
       plt.show()
   # comm.Barrier()
 
#==================================================================================
# new style: pass through dictionary
#==================================================================================

   matprop        = np.zeros(12) # SI units for mat properties
   matprop[0]     = 68.9e9          # Young's modulus 
   matprop[1]     =  0.33e0         # Poisson's ratio 
   matprop[2]     = 2700.0          # Density 
   matprop[3]     = 3.33e-5         # CS area, initial 
   matprop[4]     = 9.259e-11       # Iy, initial 
   matprop[5]     = 9.259e-11       # Iz, initial 
   matprop[6]     = 2.259e-10       # Ix, initial 
   matprop[9]     = 2.8867e-3       # max y dimension 
   matprop[10]    = 2.8867e-3       # max z dimension 
   matprop[11]    = 276.e6          # Yield stress for Von Mises 

   info           = {'nnodes': nnodes, 'nelem': nelem, 'rxyz': rxyz,  'conn': conn, 
                     'memid' : memberid, 'nloading':nloadings, 'maxnforce':maxnforce,
                     'nforces': nforces, 'forceid': forceOut, 'ndisp':ndisp, 'dispid': dispid, 
                     'matprop': matprop, 'cstype': 'hollowcircle', 'print': False} 

#   print wing_nodes,orig_nodes   
   return info     

# ===================================================================
# END OF FILE
# ===================================================================


# ===================================================================
# USE LINES BELOW TO GENERATE A TEMPLATE
# ===================================================================

R           = 1.0
b           = 4.0

params      = {'R':R, 'b': b, 'T':100.0, 'Qh':0.0, 
               'Qc':0.0, 'D': 10.0, 'L':100.0}

createAirframeConfig(params, True)
# ===================================================================
# END OF FILE
# ===================================================================

