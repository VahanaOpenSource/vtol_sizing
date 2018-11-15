#
# file to visualize cost breakdown
#
# acquisition costs
# annual (fixed) costs
# operating costs (fixed and variable)

import yaml,sys,os
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
#==========================
# use these lines for latex
#==========================

rc('text',usetex=True)
font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 18}

rc('font', **font)
#==========================

#=====================================
# use this line to save images to PDF
#=====================================

from matplotlib.backends.backend_pdf import PdfPages

#=====================================

class breakdown:
   
# ===================================================================
# initialize and load yaml file
# ===================================================================

   def __init__(self,n):

      DIR         = 'output/logs/'
      fname       = DIR+'log'+str(n)+'.yaml'
      print(" -- Chosen design ID",n)

      try:
         with open(fname) as f:
            fyaml=yaml.load(f)
      except:
         print('looked for file',fname)
         print('CRITICAL ERROR: file does not exist')
         quit()

# ===================================================================
# basic dictionary pointers and parameters
# ===================================================================

      self.costs        = fyaml['Costs']
      self.id           = n 
      self.weights      = fyaml['Weights']
      self.drags        = fyaml['Parasitic_drag']

#====================================================================
# This function creates a pie chart for acquisition cost breakdown
#====================================================================

   def var_cost_pie(self, angle): 

      var               = self.costs['variable_cost_breakdown']
      var_total         = self.costs['Variable_operating_costs']
      labels            = []
      values            = []
      expl              = []
      expl_val          = 0.05

#====================================================================
#loop over all elements in acquisition cost
#====================================================================

      for key in var.keys():#[-1:0:-1]:

         value          = var[key]
         percent        = value[1]
         labels.append(key.replace('_',' ') + ': ' + str(round(value[1],1)) + '%')
         values.append(value[1])
         expl.append(expl_val)
      
      fname             = 'var_cost_design_' + str(self.id) + '.png'

#====================================================================
# calculate miscellaneous percentages
# plot secondary pie if applicable
#====================================================================

      Total             = float(var_total[0])

#====================================================================
# draw primary pie
#====================================================================

      t1                = 'Operating cost breakdown'
      data              = {'values':values,'labels':labels,'expl':expl}

      self.draw_pie(data,fname,expl_val,t1,'\\textdollar'+str(int(Total))+ '/hr', angle)

#====================================================================
# This function creates a pie chart for operating cost breakdown
#====================================================================

   def fixed_cost_pie(self, angle): 

      fix               = self.costs['fixed_cost_breakdown']
      fix_total         = self.costs['Fixed_operating_costs']
      labels            = []
      values            = []
      expl              = []
      expl_val          = 0.025

#====================================================================
#loop over all elements in acquisition cost
#====================================================================

      for key in fix.keys():

         value          = fix[key]
         percent        = value[1]
         labels.append(key.replace('_',' ') + ': ' + str(round(value[1],1)) + '%')
         values.append(value[1])
         expl.append(expl_val)
      
      fname             = 'fixed_cost_design_' + str(self.id) + '.png'
      title             = 'fixed cost breakdown'

#====================================================================
# calculate miscellaneous percentages
# plot secondary pie if applicable
#====================================================================

      Total             = float(fix_total[0])

#====================================================================
# draw primary pie
#====================================================================

      t1                = 'Annual fixed cost breakdown'
      data              = {'values':values,'labels':labels,'expl':expl}
      self.draw_pie(data,fname,expl_val,t1,'\\textdollar' + str(int(Total)), angle)      

#====================================================================
# This function creates a pie chart for weight breakdown
#====================================================================

   def weights_pie(self, angle): 

      ewb               = self.weights['empty_weight']

      f                 = float(self.weights['battery'][0])
      p                 = float(self.weights['payload'][0])
      ewp1              = ewb['total']
      e                 = float(ewp1[0])
      ewp               = float(ewp1[1])
      W                 = e + f + p

      labels            = []
      values            = []
      expl_val          = 0.02

#====================================================================
#loop over all elements in acquisition cost
#====================================================================

      for key,val in ewb.items():
         if key != 'total':
            if(isinstance(val,dict)):
               v     = val['total']
            else:
               v     = val 
            percen_e = v[1]/ewp*100.0
            labels.append(key.replace('_',' ') + ': ' + str(int(round(v[0],1))) + ' kg')
            values.append(percen_e)
      
#====================================================================
# draw empty weight breakdown pie
#====================================================================

      fname             = 'weights_empty_design_' + str(self.id) + '.png'
      t1                = 'Empty weight breakdown'
      data              = {'values':values,'labels':labels}
      self.draw_pie(data,fname,expl_val,t1,str(int(round(e,0))) + ' kg', angle)

#====================================================================
# draw big-picture weight breakdown pie
#====================================================================

      values            = [] 
      labels            = []
      self.weights['empty_weight']['total'] = ewp1
      for key,val in self.weights.items():
         if(isinstance(val,dict)):
            v     = val['total']
         else:
            v     = val 
         labels.append(key.replace('_',' ') + ': ' + str(int(round(v[0],1))) + ' kg')
         values.append(v[0])

      fname             = 'weights_design_' + str(self.id) + '.png'
      t1                = 'Vehicle weight breakdown'
      data              = {'values':values,'labels':labels}
      self.draw_pie(data,fname,expl_val,t1, str(int(round(W,0))) + ' kg', angle)      

#====================================================================
# This function creates a pie chart for flat-plate area breakdown
#====================================================================

   def flat_plate_breakdown(self,angle): 

      f                 = self.drags
      ftot              = float(f['flat_plate_area'])

      try:
         f_details      = f['flat_plate_breakdown']
      except:
         print('could not find flat-plate area; parametric model may have been used')
         return 
      labels            = []
      values            = []
      expl_val          = 0.1

#====================================================================
#loop over all elements in acquisition cost
#====================================================================

      for key,v in f_details.items():
         percent        = v[0]
         labels.append(key.replace('_',' ') + ': ' + str(round(percent,1)) + '%')
         values.append(percent)
      
#====================================================================
# draw empty weight breakdown pie
#====================================================================

      fname             = 'parasitic_drag_design_' + str(self.id) + '.png'
      t1                = 'Equivalent flat-plate area = ' + str(round(ftot,3)) + \
                          'sq.m: breakdown (wings not included)'
      data              = {'values':values,'labels':labels}
      self.draw_pie(data,fname,expl_val,t1,angle)      

#====================================================================
# This function creates a pie chart for acquisition cost breakdown
#====================================================================

   def acqusition_pie(self,angle): 

      acq               = self.costs['acquisition_cost_breakdown']
      acq_total         = self.costs['Frame_acquisition']
      labels            = []
      values            = []
      expl              = []
      misc_val          = [] 
      misc_key          = []
      misc_expl         = []
      misc_total        = 0.0
      expl_val          = 0.3

#====================================================================
#loop over all elements in acquisition cost
#====================================================================

      for key in acq.keys():

         value          = acq[key]
#====================================================================
# if % of acquisition cost < 2%, groups separately
#====================================================================

         percent        = value[1]
         if(percent < 2.5):
            misc_key.append(key)
            misc_val.append(value[1])
            misc_expl.append(expl_val)
            misc_total += value[1]           # total % contribution from miscellaneous components 
         else:
            labels.append(key.replace('_',' ') + ': ' + str(round(value[1],1)) + '%')
            values.append(value[1])
            expl.append(expl_val)
      
      if(misc_total >0):
         labels.append('misc: ' + str(round(misc_total,2)) + '%')
         values.append(misc_total)
         expl.append(expl_val)

      fname             = 'acquisition_design_' + str(self.id) + '.png'
      fname2            = 'acquisition_miscel_' + str(self.id) + '.png'

#====================================================================
# draw primary pie
#====================================================================

      Total             = float(self.costs['Frame_acquisition'][0])*1e6
      t1                = 'Acquisition cost breakdown'
      data              = {'values':values,'labels':labels,'expl':expl}
      self.draw_pie(data,fname,0.2,t1,'\\textdollar'+str(int(Total)),angle)      

#====================================================================
# calculate miscellaneous percentages
# plot secondary pie if applicable
#====================================================================

      if(misc_total > 0):
         misc_key2        = []
         for il,l in enumerate(misc_key):
            misc_val[il]   = misc_val[il]/misc_total*100.0
            misc_key2.append(l.replace('_',' ') + ': ' + str(round(misc_val[il],2)) + '%')

         t2                = 'Miscelleneous acquisition cost breakdown'
         data              = {'values':misc_val,'labels':misc_key2,'expl':misc_expl}
         misc_total        = misc_total*Total/100.0
         self.draw_pie(data,fname2,0.02,t2,'\\textdollar'+str(int(round(misc_total,0))),angle)

#====================================================================
# This function creates a pie chart from a dictionary and saves it 
# to a png file of given name
#====================================================================

   def draw_pie(self,inp, savefile, expl,titlestr,cstr,angle): 

#====================================================
      fig, ax        = plt.subplots(figsize=(14, 9), subplot_kw=dict(aspect="equal"))

      recipe         = []
      for label in inp['labels']:
         recipe.append('\\textbf{' + label.replace('%','\%') + '}')

      data           = inp['values']
      exp            = np.ones(len(inp['labels']))*expl

      wedges, texts  = ax.pie(data, wedgeprops=dict(width=0.5), startangle= angle, explode=exp)

      bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
      kw = dict(xycoords='data', textcoords='data', arrowprops=dict(arrowstyle="-"),
             bbox=bbox_props, zorder=0, va="center")

      for i, p in enumerate(wedges):
         ang = (p.theta2 - p.theta1)/2. + p.theta1
         y = np.sin(np.deg2rad(ang))
         x = np.cos(np.deg2rad(ang))
         horizontalalignment  = {-1: "right", 1: "left"}[int(np.sign(x))]
         connectionstyle      = "angle,angleA=0,angleB={}".format(ang)
         kw["arrowprops"].update({"connectionstyle": connectionstyle})
         ax.annotate(recipe[i], xy=(x, y), xytext=(1.45*np.sign(x), 1.4*y),
                    horizontalalignment=horizontalalignment, **kw)

      clen  = len(cstr)
      
      if(cstr.startswith('\\')):
         dx    = clen*0.5*0.03   
      else:
         dx    = clen*0.5*0.07
      cstr  = '\\textbf{' + cstr.replace('%','\%') + '}'
      
      ax.annotate(cstr, xy=(-dx, 0), xytext=(-dx, 0),fontsize=27)

#for latex
      t     = '\\textbf{\\underline{' + titlestr + '}}'
      t     = t.replace('%','\%')
      ax.set_title(t,y = 1.08,fontweight='bold',fontsize=25)

#      plt.tight_layout(pad=1)
#for PNG files
#      plt.savefig(savefile)
#for PDF file
      self.pdf.savefig()

      return None

#====================================================================
# driver function to render plots
#====================================================================

   def render_plots(self, angle=40):
      fname       = 'costs_design_'+str(self.id)+'.pdf'
      with PdfPages(fname) as pdf:

         self.pdf  = pdf 
         self.weights_pie(angle)
         self.var_cost_pie(0)
         self.fixed_cost_pie(0)
         self.acqusition_pie(0)
         self.flat_plate_breakdown(0)
         print('completed pie chart generation; saved to ',fname)