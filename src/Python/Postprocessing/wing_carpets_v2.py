import matplotlib.pyplot as plt
import numpy 
import matplotlib
def wing_carpets(A,B,X,Y,Z,xref,yref,zref,valid,caxis,cstr,xstr,ystr,astr,bstr):


#====================================================================
# Wings; use left halign, else use right for rotors
#====================================================================

    if(astr.startswith('AR')):
        halign = 'left'
        valign = 'top'
        offset = -1
        angle  = -30
        asuf   = ''
        rem    = 0
    else:
        halign = 'right'
        valign = 'baseline'
        offset = 0
        angle  = 0
        asuf   = astr 
        astr   = ''
        rem    = 1
#====================================================================
# set custom color map scale
#====================================================================

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["blue","white","red"])

    zmax    = min(numpy.amax(Z),zref+0.35*abs(zref)); zmax   = round(zmax,0)
    zmin    = max(numpy.amin(Z),zref-0.35*abs(zref)); zmin   = round(zmin,0)
    levels  = numpy.linspace(zmin,zmax,41)
    norm    = plt.Normalize(zmin,zmax)
#    print('min/max limits on contour plot are ',zmin,zmax)
#====================================================================
# create contour plot
#====================================================================

    plt.contourf(X, Y, Z, vmin=zmin,vmax=zmax,levels=levels,cmap=cmap)
    cbar    = plt.colorbar()

#====================================================================
# add labels and title
#====================================================================

    plt.xlabel('\\textbf{'+xstr+'}'); plt.ylabel('\\textbf{'+ystr+'}')
    plt.title('\\textbf{'+cstr+'}')
    # plt.show()
    cbar.ax.set_title('\\textbf{'+caxis+'}',horizontalalignment='left')

#====================================================================
# draw lines of constant A,B; we have x,y coords in (X,Y); 
# A, B arrays are input, and Z has data in matrix form 
# according to rows of A variations, and columns of B variations
# connect lines of X,Y
#====================================================================
    
    for ib,b in enumerate(B):
        X_dash  = X[:,ib]
        Y_dash  = Y[:,ib]
        if(ib % 2 == 0):
            plt.plot(X_dash,Y_dash,'--',linewidth=0.7,color='0.5')
            plt.text(X_dash[offset],Y_dash[offset],'\\textbf{'+bstr+str(round(b,2))+'}',color='0.5',fontsize=10)
        else:
            plt.plot(X_dash,Y_dash,'--',linewidth=0.5,color='0.5')

#====================================================================
# draw lines for A variation
#====================================================================

    for ia,a in enumerate(A):
        X_dash  = X[ia,:]
        Y_dash  = Y[ia,:]
        if(ia % 2 == rem):
            plt.plot(X_dash,Y_dash,'--',color='0.5',linewidth=0.7)
            plt.text(X_dash[0],Y_dash[0],'\\textbf{'+astr+str(int(round(a,0)))+asuf+'}',rotation=angle,\
                color='0.5',horizontalalignment=halign,verticalalignment=valign,fontsize=10)
        else:
            plt.plot(X_dash,Y_dash,'--',color='0.5',linewidth=0.5)

#====================================================================
# put a red "x" on invalid designs
#====================================================================

    nx,ny       = numpy.shape(valid)
    ind_i       = []
    ind_j       = []
    for i in range(nx):
        for j in range(ny):
            if(valid[i,j] == 0):
                ind_i.append(i)
                ind_j.append(j)

    Xv          = X[ind_i,ind_j]
    Yv          = Y[ind_i,ind_j]
    plt.plot(Xv,Yv,'rx',markersize=4)
    plt.plot(Xv,Yv,'o',markersize=8,markerfacecolor='none',color='0.75')

    plt.plot(xref,yref,'--o',color='orange',linewidth=2.4,markersize=12,markerfacecolor='none',markeredgewidth=2.4)

#    plt.plot(Xv,Yv,'.',markersize=2)
#====================================================================
# set plot axes limits
#====================================================================

    xmin        = numpy.amin(X);
    xmax        = numpy.amax(X);
    ymin        = numpy.amin(Y);
    ymax        = numpy.amax(Y); 

    if(astr.startswith('AR')):
        xmax = xmax + 0.25*abs(xmax)
        ymax = ymax + 0.3*abs(ymax)
        ymin = ymin - 0.5*abs(ymin)
        xmin = xmin - 0.25*abs(xmin)
    else:
        xmax = xmax + 0.1*abs(xmax-xmin)
        xmin = xmin - 0.1*abs(xmax-xmin)
        ymax = ymax + 0.1*abs(ymax-ymin)
        ymin = ymin - 0.1*abs(ymax-ymin)
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    plt.tight_layout(pad=1.1)
#    plt.show()        
#====================================================================
# save file
#====================================================================

    return None