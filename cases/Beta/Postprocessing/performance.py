#====================================================================
# Main function
#====================================================================

#====================================================================
# import modules used
#====================================================================

from __future__ import print_function
import sys; sys.path.append('../../src/Python/Postprocessing/')

from vehicle_performance import vehicle_performance

#====================================================================
# error traps
#====================================================================

if len(sys.argv) != 2:
   print(" -- Error: performance.py requires a design ID")
   print(" -- Usage: 'python performance.py <ID> '")

#====================================================================
# get design ID and build file name, write status message
#====================================================================

n        = int(sys.argv[1])

DIR      = 'output/logs/'
fname    = DIR + 'log'+str(n)+'.yaml'

print(" -- Chosen design ID",n)

#====================================================================
# Initialize class and load data
#====================================================================

perf 	 = vehicle_performance(fname, n, 'Postprocessing/table.data')
perf.make_plots()

