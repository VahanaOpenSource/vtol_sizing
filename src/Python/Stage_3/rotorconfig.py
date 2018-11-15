import sys
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

#from matplotlib.patches import Ellipse, Arc
from matplotlib import cm, colors, patches
#
class Arrow3D(FancyArrowPatch):

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

class _rotorconfig:

   def arrangeRotors(self,rotor,nr,rb,rs):
              
      if(nr != rotor['nrotor']):
         print 'Incorrect number of rotor'
         sys.exit(1)
   
      # central location
      center = [0.0,0.0]
      
      nrotor   = rotor['nrotor']
      ncluster = rotor['ncluster']
      
      fig = plt.figure()
      ax  = plt.subplot(projection='3d')
      
      ntheta = 36
      theta  = np.linspace(0,2.0*np.pi,num=ntheta)
      for j in range(nrotor):
         xpos   = rotor['rotor_position'][j][0]
         ypos   = rotor['rotor_position'][j][1]
         zpos   = 0.0
         radius = rotor['radius'][j]
      
         x = xpos + radius * np.cos(theta)
         y = ypos + radius * np.sin(theta)
         z = zpos + np.zeros(ntheta)
      
         # plot circle
         plt.plot(x,y,z,'r-')
      
         # draw thrust vector
         dzpos=0.015
         a = Arrow3D([xpos, xpos], [ypos, ypos], [zpos, zpos+dzpos], mutation_scale=20,
                     lw=1, arrowstyle="-|>", color="b")
         ax.add_artist(a)
      
      # draw payload vector
      dzpos=0.015
      xpos=center[0]; ypos=center[0]; zpos=0
      a = Arrow3D([xpos, xpos], [ypos, ypos], [zpos, zpos-dzpos], mutation_scale=20,
                  lw=1, arrowstyle="-|>", color="green")
      ax.add_artist(a)
      
      # plot lines
      for j in range(nrotor):
         clusterID = rotor['clusterID'][j]
         xrotor    = rotor['rotor_position'][j][0]
         yrotor    = rotor['rotor_position'][j][1]
         xclus     = rotor['cluster_position'][clusterID][0]
         yclus     = rotor['cluster_position'][clusterID][1]
         
         # plot center to cluster center
         plt.plot([center[0],xclus],[center[1],yclus],linestyle='--',color='black',linewidth=2,marker='o')
      
         # plot from cluster center to rotor center
         plt.plot([xclus,xrotor],[yclus,yrotor],linestyle='--',color='black',linewidth=2,marker='o')
      
      # figure label 
      plt.xlim(-7.5*rb,7.5*rb)
      plt.ylim(-5*rb,5*rb)
      plt.gca().set_aspect('equal',adjustable='box')
      plt.grid()
      ax.set_xlabel('x/R')
      ax.set_ylabel('y/R')
      ax.set_zlabel('z/R')
      
      plt.draw()
      
      # ===================================================================
      # Output to rotor.inp (needed by FEA code)
      # ===================================================================
      
      BASE = 1
      
      nBeamElem = 5 # number of elements for each connection
      
      
      nelem = ncluster * nBeamElem + nrotor * nBeamElem
      
      nnodes = ncluster * (nBeamElem+1) - ncluster + 1 \
             + nrotor   * (nBeamElem+1) - nrotor
      
      rxyz = np.zeros([nnodes,3])
      conn = np.zeros([nelem,2])
      
      # special node indices
      idxCluster = np.zeros(ncluster+1)
      for j in range(ncluster):
         idxCluster[j+1] = int(nBeamElem*(j+1) + BASE)
      
      idxRotor = np.zeros(nrotor)
      for j in range(nrotor):
         idxRotor[j] = nBeamElem*ncluster  + nBeamElem*(j+1) + BASE
      
      # ===================================================================
      # coordinate data
      # ===================================================================
      
      # connections from center to cluster centers
      kk1 = 0
      posEnd1 = [center[0], center[1], 0.0]
      for j in range(ncluster):
         xclus = rotor['cluster_position'][j][0]
         yclus = rotor['cluster_position'][j][1]
         posEnd2 = [xclus,yclus,0.0]
      
         for k in range(nBeamElem):
            
            if(j>0 and k==0):
               continue;
      
            rxyz[kk1][0] = posEnd1[0] + (posEnd2[0]-posEnd1[0])*k/(nBeamElem)
            rxyz[kk1][1] = posEnd1[1] + (posEnd2[1]-posEnd1[1])*k/(nBeamElem)
            kk1=kk1+1
         
         # add the last point
         rxyz[kk1][0] = posEnd2[0]
         rxyz[kk1][1] = posEnd2[1]
         kk1 = kk1 + 1
      
      # connections from cluster center to rotors
      for j in range(nrotor):
         clusterID = rotor['clusterID'][j]
         xclus     = rotor['cluster_position'][clusterID][0]
         yclus     = rotor['cluster_position'][clusterID][1]
         posEnd1   = [xclus,yclus,0.0]
      
         posEnd2[0] = rotor['rotor_position'][j][0]
         posEnd2[1] = rotor['rotor_position'][j][1]
         posEnd2[2] = 0.0
      
         for k in range(nBeamElem):
            
            if(k==0):
               continue
      
            rxyz[kk1][0] = posEnd1[0] + (posEnd2[0]-posEnd1[0])*k/(nBeamElem)
            rxyz[kk1][1] = posEnd1[1] + (posEnd2[1]-posEnd1[1])*k/(nBeamElem)
            kk1=kk1+1
         
         # add the last point
         rxyz[kk1][0] = posEnd2[0]
         rxyz[kk1][1] = posEnd2[1]
         kk1 = kk1 + 1
      
      # ===================================================================
      # connectivity data
      # ===================================================================
      
      kk2 = 0
      kk3 = 0 
      memberID = np.zeros(nelem)
      # for center -> cluster
      for j in range(ncluster):
         for k in range(nBeamElem):
      
            n1 = nBeamElem * j  + k + BASE
            n2 = n1 + 1
            
            if (k==0):
               n1 = 1
      
            if(k==nBeamElem-1):
               n2 = idxCluster[j+1]
        
            conn[kk2][0]  = int(n1)
            conn[kk2][1]  = int(n2)
            memberID[kk2] = int(kk3 + BASE)
            kk2 = kk2 + 1
         
         kk3 = kk3 + 1
      
      
      # for cluster center to rotor
      for j in range(nrotor):
         
         clusterID = rotor['clusterID'][j]
      
         for k in range(nBeamElem):
      
            n1 = (nBeamElem*ncluster) + nBeamElem * j  + k + BASE
            n2 = n1 + 1
            
            if (k==0):
               n1 = idxCluster[clusterID+1]
      
            if(k==nBeamElem-1):
               n2 = idxRotor[j]
        
            conn[kk2][0]  = int(n1)
            conn[kk2][1]  = int(n2)
            memberID[kk2] = int(kk3 + BASE)
            kk2 = kk2 + 1
         
         kk3 = kk3 + 1
      
      #
      f=open('config.inp','w')
      
      f.write("# nnodes, nelem\n")
      f.write(" %d %d\n" % (nnodes,nelem) )
      f.write("# node positions\n")
      for j in range(nnodes):
         f.write("%e %e %e\n" % (rxyz[j][0], rxyz[j][1],rxyz[j][2]) )
      f.write("# element connectivity\n")
      for j in range(nelem):
         f.write("%d %d %d\n" % (conn[j][0], conn[j][1], memberID[j] ))
      
      f.close()
      
      plt.show()
   
# # ===================================================================
# # my rotor
# 
# 
# nr     = 4
# radius = 0.914 #m
# 
# rb     = 1 * radius
# rs     = 0.5*rb
# 
# 
# 
# arrangeRotors(nr,rb,rs)


# ===================================================================
# END OF FILE
# ===================================================================

