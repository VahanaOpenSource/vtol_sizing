import matplotlib.pyplot as plt
import numpy 
import matplotlib
def line_sensitivities(x,y,y2,xref,yref,yref2,valid,xstr,ystr,ystr2,title):

#====================================================================
# create line plot
#====================================================================

    fig, ax1 = plt.subplots()

    ax1.plot(x, y ,'-',linewidth=2.0)

#====================================================================
# add labels and title
#====================================================================

    ax1.set_xlabel('\\textbf{'+xstr +'}')
    ax1.set_ylabel('\\textbf{'+ystr +'}',color='C0')
    plt.title('\\textbf{'+title+'}')
    ax1.tick_params('y', colors='C0')

#====================================================================
# plot second set
#====================================================================

    ax2      = ax1.twinx()
    ax2.plot(x, y2,'-',linewidth=2.0,color='r')
    ax2.set_ylabel('\\textbf{'+ystr2+'}',color='r')
    ax2.tick_params('y', colors='r')

#====================================================================
# put a red "x" on invalid designs
#====================================================================

    v   = [i for i, x in enumerate(valid) if x == 0]
    Xv  = x[v]
    Yv  = y[v]
    Yv2 = y2[v]

    ax1.plot(Xv,Yv,'rx',markersize=4)
    ax1.plot(Xv,Yv,'o',markersize=8,markerfacecolor='none',color='0.75')

    ax2.plot(Xv,Yv2,'rx',markersize=4)
    ax2.plot(Xv,Yv2,'o',markersize=8,markerfacecolor='none',color='0.75')

    ax1.plot(xref,yref,'--o',color='C0',linewidth=2.4,markersize=12,markerfacecolor='none',markeredgewidth=2.4)
    ax1.plot(xref,yref,'.',color='C0',linewidth=2.4,markersize=4,markeredgewidth=2.4)

    ax2.plot(xref,yref2,'--o',color='r',linewidth=2.4,markersize=12,markerfacecolor='none',markeredgewidth=2.4)
    ax2.plot(xref,yref2,'.',color='r',linewidth=2.4,markersize=4,markeredgewidth=2.4)

#====================================================================
# draw a vertical black line at airspeed of interest
#====================================================================

    ylims   = ax1.get_ylim()
    xlims   = (xref,xref)
    ax1.plot(xlims,ylims,'--k')
    ax1.text(xref,ylims[1],'design point',rotation=0,verticalalignment='bottom')

#====================================================================
# set plot axes limits
#====================================================================

    xmin    = numpy.amin(x)
    xmax    = numpy.amax(x)
    ymin    = numpy.amin(y)
    ymax    = numpy.amax(y) 

    dx      = 0.1*abs(xmax-xmin)
    xmax    = xmax + 0.1*abs(xmax-xmin)
    xmin    = xmin - 0.1*abs(xmax-xmin)
    ymax    = ymax + 0.1*abs(ymax-ymin)
    ymin    = ymin - 0.1*abs(ymax-ymin)

    plt.xlim([xmin,xmax])
#    plt.ylim([ymin,ymax])
    plt.tight_layout(pad=1.1)
    # plt.show()        
#====================================================================
# save file
#====================================================================

    return None