import plotly.graph_objs as go
#import plotly.plotly as py

import plotly.offline as offline
from selenium import webdriver
# from matplotlib import rc 
# #==========================
# # use these lines for latex
# #==========================

# rc('text',usetex=True)
# font = {'family' : 'serif',
#         'weight' : 'bold',
#         'size'   : 18}

# rc('font', **font)

def draw_a_carpet(A,B,Z,X,Y,save_file,valid,cstr,xstr,ystr):

    minval          = round(min(Z),-1)
    maxval          = round(max(Z),-1)
    wtstep          = (maxval-minval)/10

    print('range to plot is ',minval,maxval)
#define the font
    font1 =dict(
        family='Helvetica',
        size=18
    )

#draw color contour
    trace1 = go.Contourcarpet(
        a = A,
        b = B,
        z = Z,
        autocontour = False,
        contours = dict(
            start = minval,
            end = maxval,
            size = wtstep
        ),
        line = dict(
            width = 0.5,
            smoothing = 0
        )
        , colorbar = dict(
        #     len = 0.4,
             x = 0.9,
        title=cstr,
        tickfont=font1,
        titlefont=dict(size=16)
        )
    )


#define carpet location
    trace2 = go.Carpet(
        a = A,
        b = B,
        x = X,
        y = Y,
        aaxis =        dict(
           # tickprefix = 'DL = ',   
           smoothing = 0,
           nticks   = 4,
           type = 'linear',
           showticklabels='none',
#           categoryorder='category '
       ),
        baxis = dict(
#            tickprefix  = 'sigma=',
#            smoothing   = 0,
#            tickmode    = 'array',
#            tickvals    = [min(B),max(B)],
#            tickfont    = font1,
           showticklabels='none',
#            nticks      = 3,
        )
    )

#=========================================================
# highlight valid designs
#=========================================================

    v   = [i for i, x in enumerate(valid) if x == 1]
    x   = X[v]
    y   = Y[v]
    trace3 =     go.Scatter(
            mode        = 'markers',
            x           = x,
            y           = y,
            opacity     = 0.2,
            marker  = dict(symbol='circle-dot',size = 10,line = dict(color = 'black', width = 2),color='white'),
    showlegend = False    )

#=========================================================
# highlight min and max
#=========================================================

    Zlist       = list(Z[v])
    highest     = max(Zlist)
    lowest      = min(Zlist)
    ihigh       = [ list(Z).index(highest) ]
    ilow        = [ list(Z).index(lowest)  ]
    print(ihigh,ilow,highest,lowest)
    trace4 =     go.Scatter(
            mode    = 'markers',
            x       = X[ilow], 
            y       = Y[ilow],
            marker  = dict(symbol='circle-dot',size = 12,line = dict(color = 'green', width = 4),color='white'),    showlegend = False    )

    trace5 =     go.Scatter(
            mode    = 'markers',
            x       = X[ihigh], 
            y       = Y[ihigh],
            marker  = dict(symbol='circle-dot',size = 12,line = dict(color = 'red', width = 4),color='white'),    showlegend = False    )

#=========================================================
# call plotly interface
#=========================================================

    data = [trace1, trace2,trace3,trace4,trace5]
    #data = [trace1]#, trace2]

    layout = go.Layout(
        margin = dict(
            t = 30,
            r = 60,
            b = 60,
            l = 70
        ),
        yaxis = dict(
#            range = [0.0,0.18],
            title = ystr,
            titlefont = font1,
            tickfont  = font1
        ),
        xaxis = dict(
#            range = [min(trace2.x)-2,max(trace2.x)+2],
            title = xstr,
            titlefont = font1,
            tickfont  = font1
        )
    )

#    try:
#    fig = go.Figure(data = data, layout = layout)
#        py.image.save_as(fig, filename=save_file+'.png')
    offline.plot({'data':data,'layout':layout},image='svg', auto_open=False, image_width=1000, image_height=500,show_link=False)
    try:
        driver = webdriver.PhantomJS(executable_path="/usr/local/bin/phantomjs")
    except:
        driver = webdriver.PhantomJS(executable_path="/usr/bin/phantomjs")
    driver.set_window_size(1000, 500)
    driver.get('temp-plot.html')
    driver.save_screenshot(save_file + '.png')
    print('saving file',save_file+'.png')

#    except:
#        py.sign_in('nenanth', 'eLGDoeHWkICx2mqrdpUu') # Replace the username, and API key with your credentials.
#        fig = go.Figure(data = data, layout = layout)
#        py.image.save_as(fig, filename=save_file+'.png')

    return None