#====================================================================
# Driver program
#====================================================================
import sys 

sys.path.append('../../src/Python/Postprocessing/')

from pie_generator import breakdown 

#====================================================================
# get design ID
#====================================================================

if len(sys.argv) < 2:
   print(" -- Error: piethon.py requires a design ID")
   print(" -- Usage: 'python piethon.py <ID> '")
   quit()

else:
   n           = sys.argv[1]
   design      = breakdown(n)

#====================================================================
# draw pie charts
#====================================================================

   design.render_plots()


