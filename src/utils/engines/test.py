# python code to test engine library

import os,sys
import yaml
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

# ==================================================================
# load turboshaft data
# ==================================================================

class turboshaft:

# ==================================================================
# class initialization
# ==================================================================

   def __init__(self):
      with open("turboshaft.yaml","r") as f:
         fturbo = yaml.load(f)
       
      engine = fturbo['Turboshaft']
      
      # total number of engines
      neng = len(engine)
   
      self.sfc_to   = np.zeros(neng)
      self.power_to = np.zeros(neng)
      self.weight   = np.zeros(neng)
      #
      # create plot of sfc_to vs p_to
      #
      kk = 0
      for j in range(neng):
         if j < 10:
            eng = 'eng00' + str(j)
         elif j < 100:
            eng = 'eng0' + str(j)
         else:
            eng = 'eng' + str(j)
      
         powerTmp  = engine[eng]['power_to']
         sfcTmp    = engine[eng]['sfc_to']
      
         if powerTmp > 0.0 and sfcTmp > 0.0:
            self.power_to[kk] = powerTmp
            self.sfc_to[kk]   = sfcTmp
            kk = kk + 1
      
      self.sfc_params = scipy.optimize.curve_fit(lambda t, a,b: a+b*np.log(t), self.power_to, self.sfc_to)

      #
      # create plot of weight vs p_to
      #
      kk = 0
      for j in range(neng):
         if j < 10:
            eng = 'eng00' + str(j)
         elif j < 100:
            eng = 'eng0' + str(j)
         else:
            eng = 'eng' + str(j)
      
         powerTmp  = engine[eng]['power_to']
         weightTmp = engine[eng]['weight']
      
         if powerTmp > 0.0 and weightTmp > 0.0:
            self.power_to[kk] = powerTmp
            self.weight[kk]   = weightTmp
            kk = kk + 1
      
      self.weight_params = scipy.optimize.curve_fit(lambda t, a, b, c : a+b*t+c*t*t, self.power_to, self.weight)

# ==================================================================
# get SFC for a given power
# ==================================================================

   def getEngineSFC(self,power):
      
      sfc = self.sfc_params[0][0] + self.sfc_params[0][1] * np.log(power)

      return sfc

# ==================================================================
# get weight for a given power
# ==================================================================

   def getEngineWeight(self,power):
      
      weight = self.weight_params[0][0] + self.weight_params[0][1] * power + self.weight_params[0][2]*power*power

      return weight

# ==================================================================
# get weight for a given power
# ==================================================================

   def plotTrends(self):

      len_power_range = 50
      power_range  = np.linspace(100.0, 12000.0, num=len_power_range)
      weight_range = np.zeros(len_power_range)
      for j in range(len_power_range):
         a = self.weight_params[0][0]
         b = self.weight_params[0][1]
         c = self.weight_params[0][2]
      
         weight_range[j] = a + b * power_range[j] + c * power_range[j] * power_range[j]
      
      plt.plot(self.power_to,self.weight,linestyle='None',marker='x')
      plt.plot(power_range,weight_range,color='red',linestyle='--')
      plt.ylabel('Weight, lb')
      
      #plt.plot(power_to,sfc_to,linestyle='None',marker='o')
      #plt.plot(power_range,sfc_range,color='red',linestyle='--')
      #plt.ylabel('SFC, lb/hp-hr')

      plt.xlabel('Takeoff power, hp')
      plt.grid()
      plt.show()
   



x = turboshaft()
x.plotTrends()
