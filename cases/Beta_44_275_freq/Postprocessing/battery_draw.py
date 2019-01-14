#====================================================================
# Driver program
#====================================================================
import sys 

sys.path.append('../../src/Python/Postprocessing/')

from battery_profile import battery_profile 

#====================================================================
# get design ID
#====================================================================

if len(sys.argv) < 2:
   print(" -- Error: battery_draw.py requires a design ID")
   print(" -- Usage: 'python battery_draw.py <ID> <optim>'")
   quit()

if len(sys.argv) > 2:
    if(sys.argv[2] == 'optim'):
        print('postprocessing optimized results')
        DIR         = 'output/logs/optim_'
    else:
        quit("third argument must be 'optim'")
else:
    print('postprocessing parameter sweep results')
    DIR         = 'output/logs/'

n           = sys.argv[1]
file        = ''
fname       = DIR+'log'+str(n)+'.yaml'

design      = battery_profile(fname,n)
design.render_plots()

#====================================================================
# draw acquisition pie chart
#====================================================================

