# Analyze hover margins

# Tilt angles to test

#==========================
# use these lines for latex
#==========================

from matplotlib import pyplot as plt 
from matplotlib import rc
import numpy
rc('text',usetex=True)
font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 12}

rc('font', **font)

from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D

pdf         = PdfPages('motor_out.pdf')
plotOn      = True

#0 to 360 to make a circle representing a rotor 
th          = numpy.linspace(0,2*numpy.pi,60)
costh       = numpy.cos(th)
sinth       = numpy.sin(th)

#==========================

import numpy,time 
from scipy.optimize import linprog
from setup_lp import setup_lp, get_rotor_dirn
from run_lp import run_lp
d2r         = numpy.pi/180.0
g           = 9.81      # accln due to gravity, m/s^2

#=================================================================================
# constants and options
#=================================================================================

dMzdTau     = 0.5       # Ratio of torque to aero torque (wing recovers some swirl)
m           = 725       # vehicle mass, kg
W           = m*g 
rho         = 1.225     # density, kg/cu.m
nrotor      = 8

#=================================================================================
# rotor radius, shaft mounting angles wrt wing, mount locations 
#=================================================================================

r           = 0.75      # rotor radius, meters
radius      = numpy.asarray([r,r,r,r,r,r,r,r])
alpha_s     = numpy.asarray([3,5,5,3,3,1,1,3]).astype(float)
y_rotor     = numpy.asarray([2.9,1.23,-1.23,-2.9,2.9,1.23,-1.23,-2.9])
x_rotor     = numpy.asarray([2.125,2.125,2.125,2.125,-2.125,-2.125,-2.125,-2.125])

#=================================================================================
# Estimate effect of control surfaces to set max yaw moment
#=================================================================================

dClda       = 0.0016      # moment per degree of deflection from AVL
dClde       = 0.0006      # moment per degree of deflection from AVL
T           = m * g / 8 # Average motor thrust at hover
v           = numpy.sqrt(T / (2 * rho * numpy.pi * r*r)) # Induced velocity, m/s
MzSurface   = 0.5 * rho * v*v * 8.8113 * 6.1882 * (dClde + dClda) * 20 # => ~100 Nm yawing moment available at hover with 20 deg deflections

#=================================================================================
# Setup optimization problem: bounds, inequality and equality constraints
#=================================================================================

ndof        = 4
ncontrols   = nrotor+1
linp_dict   = setup_lp(nrotor,MzSurface,m*g)

#=================================================================================
# Define +1 as right handed rotation about thrust direction
#=================================================================================

kkrange     = [0,1,2,3,4,5,6]
als         = numpy.zeros_like(alpha_s)

#=================================================================================
# Vehicle properties
#=================================================================================

dTdTau      = 0.01409146 / (r * 0.00162954) # Thrust per torque = Ct / r*Cq

dTdTau_arr  = []
for i in range(nrotor):
    dTdTau_arr.append(dTdTau)

dTdTau_arr  = numpy.asarray(dTdTau_arr)

#=================================================================================
# wing tilt angles
#=================================================================================

fwdTiltTest = [-8,-4,-2,0,2,4,8] 
aftTiltTest = [-8,-4,-2,0,2,4,8]

#=================================================================================
# initialize output arrays
#=================================================================================

nfwd        = len(fwdTiltTest)
naft        = len(aftTiltTest)
maxThrust   = numpy.zeros((nfwd,naft,len(kkrange)))
maxMz       = numpy.zeros_like(maxThrust)

#=================================================================================
# loop over all wing tilt angles
#=================================================================================

for ii in range(nfwd):                  # loop over fore wing tilt angles

    fwdTilt = fwdTiltTest[ii]

    for jj in range(naft):              # loop over  aft wing tilt angles

        aftTilt = aftTiltTest[jj]

#set resultant shaft angles wrt body -z axis
        for i in range(nrotor):
            if(i <= 3):
                als[i]  = (alpha_s[i] + fwdTilt)*d2r
            else:
                als[i]  = (alpha_s[i] + aftTilt)*d2r

        for ik,kk in enumerate(kkrange):              # loop over rotor rotation directions

            kk      = kk + 1
            rotDir  = get_rotor_dirn(kk)

            x,maxthr,meantau,maxQ   =  run_lp(als, rotDir, dTdTau_arr, dMzdTau, ndof, ncontrols,  \
                                              x_rotor, y_rotor, MzSurface, linp_dict)

#=================================================================================
# save maximum torque ratio and max yawing moment 
#=================================================================================

            maxThrust[ii,jj,ik]   = maxthr / meantau
            maxMz[ii,jj,ik]       = abs(maxQ)

#=================================================================================
# create plots
#=================================================================================
 
            if plotOn and fwdTilt == 2 and aftTilt == -8:
                xPos =-y_rotor
                yPos = x_rotor

                for i in range(ncontrols):
                    plt.figure(1)
                    for j in range(nrotor):

                        if(i == 0 or j != (i-1)):
                            alpha   = 1.0    
                        else:
                            alpha   = 0.4
                        if rotDir[j] > 0:
                            plt.plot(xPos[j]+r*costh,yPos[j]+r*sinth,'b',linestyle='--',alpha=alpha)
                            plt.plot(xPos[j]+r*costh[0],yPos[j]+r*sinth[0],'b^',markersize=8,alpha=alpha)
                        else:
                            plt.plot(xPos[j]+r*costh,yPos[j]+r*sinth,'k',alpha=alpha)
                            plt.plot(xPos[j]+r*costh[0],yPos[j]+r*sinth[0],'kv',markersize=8,alpha=alpha)
                        if(j == i-1 and i > 0):
                            plt.plot(xPos[j],yPos[j],'rx',markersize=32,linewidth=3)
                        plt.text(xPos[j],yPos[j],str(round(x[ncontrols*i+j]/meantau,3)),    \
                                 horizontalalignment='center',verticalalignment='center',fontsize=10)

                    plt.grid('True',linestyle=':')
                    if(i == 0):
                        title_str   = 'Layout ' + str(kk) + ': Nominal case: ' + ' Surface Mz = ' + str(round(x[9*(i+1)-1],3)) + ' Nm'
                    else:
                        title_str   = 'Layout ' + str(kk) + ': Rotor '+str(i) + ' Out' + ' Surface Mz = ' + str(round(x[9*(i+1)-1],3)) + ' Nm'
                    plt.title(title_str)
                    plt.xlim([min(xPos)-2,max(xPos)+2]); plt.ylim([min(yPos)-2,max(yPos)+2])
                    plt.axis('equal')
                    plt.axis('off')                    ; plt.tight_layout()
                    pdf.savefig()                      ; plt.close()


#=================================================================================
# plots
#=================================================================================

if nfwd > 1 and naft > 1:

#=================================================================================
# Thrust ratio vs. wing tilt: wireplots
#=================================================================================

    fig     = plt.figure()
    ax      = fig.gca(projection = '3d')
    plt.title('Thrust Ratio vs. Wing Tilt')
        
    X,Y     = numpy.meshgrid(fwdTiltTest,aftTiltTest)
    for i,kk in enumerate(kkrange):
        ax.plot_wireframe(X,Y,maxThrust[:,:,i].T,color='C'+str(i))

    ax.set_xlabel('Fwd Tilt [deg]')
    ax.set_ylabel('Aft Tilt [deg]')
    ax.set_zlabel('Motor Size')
    ax.view_init(12, 45)
    pdf.savefig()
        # plt.show()
    plt.close()

#=================================================================================
# Thrust ratio vs. wing tilt: wireplots
#=================================================================================

    fig     = plt.figure()
    ax      = fig.gca(projection = '3d')

    plt.title('Max. yawing moment vs. Wing Tilt')

    for i,kk in enumerate(kkrange):
        ax.plot_wireframe(X,Y,maxMz[:,:,i].T,color='C'+str(i))

    ax.set_xlabel('Fwd Tilt [deg]')
    ax.set_ylabel('Aft Tilt [deg]')
    ax.set_zlabel('Max Yaw Moment [Nm]')
    ax.view_init(elev=10, azim=45)
    pdf.savefig()
    # plt.show()
    plt.close()
    

#=================================================================================
# Thrust ratio vs. yawing moment
#=================================================================================

    plt.figure()
    plt.title('Thrust Ratio vs. Yawing Moment')

    markers     = ['o','^','s']
    colors      = ['2','0','1']
    icolor      =  -1
    imark       =  0
    for ik,kk in enumerate(kkrange):

        icolor   = icolor + 1
        for ift,ft in enumerate(fwdTiltTest):
            for iat,at in enumerate(aftTiltTest):
                if(ift == 0 and iat == 0):
                    plt.plot(maxThrust[ift,iat,ik],maxMz[ift,iat,ik],'o',color='C'+str(colors[icolor]),marker=markers[imark],markersize=6,label='Layout ' + str(ik+1),markerfacecolor=None)
                else:
                    plt.plot(maxThrust[ift,iat,ik],maxMz[ift,iat,ik],'o',color='C'+str(colors[icolor]),marker=markers[imark],markersize=6,markerfacecolor=None)

        if(icolor == 2):
            icolor  = -1
            imark   = imark+1

    plt.grid(True,linestyle=':')
    plt.xlabel('Motor Thrust Ratio');
    plt.ylabel('Max Yaw Moment [Nm]')
    plt.legend(loc='best')
    pdf.savefig()
    plt.close()
    # % Design Point
    # % fwd 2, aft -8, (old fwd 0, aft -6)
    # %     i = find(fwdTiltTest == 0);
    # %     j = find(aftTiltTest == -6);
    # i = find(fwdTiltTest == 2);
    # j = find(aftTiltTest == -8);
    # plot(maxThrust(i,j,3),maxMz(i,j,3), 'bx','MarkerSize',10)


pdf.close()