#
# electric motors; assume one motor per rotor.
#
import numpy 
class battery_pack:

#====================================================================
# class instance initialization
# bunch of efficiencies
#====================================================================

   def __init__(self, Battery, nseg):

      Cell              = Battery['Cell'] 
      Pack              = Battery['Pack']

      self.energy_buf   = Pack['SOH'] - Pack['DOD_min']                 # buffer  fraction inverse for SOH,remaining reserve adjustment
      self.IF           = 1.0/Pack['integ_fac']
      self.VIF          = 1.0/Pack['vol_fac']
      self.cell_V       = Cell['volume']                                # in liters
      self.en_density   = Cell['energy_vol']                            # Watt-hours/liter
      self.dTdt         = numpy.asarray([0.172702, 0.126487, 0.0342976])   # temperature-Crating curve fit parameters
      self.sp_energy    = Cell['sp_energy']
      self.en_scaling   = 1.0/(self.energy_buf)                         # includes energy buffer, pack integration factor
      self.Tmax         = Cell['Tmax']                                  # max cell temperature
      self.E_operations = 0.0
      self.Pack         = Pack
      self.C_rating     = numpy.zeros(nseg)
      self.cell_temp    = numpy.zeros(nseg)

#====================================================================
# Force sizing battery by energy only (match GB's model)
#====================================================================

      try:
         self.en_sizing = (Battery['Force_sizing'] == 'energy')
      except:
         self.en_sizing = False 

#====================================================================
# this function calculates the battery weight for an electric motor
# based drive system for rotors/propellers, given some relevant inputs
# in a dictionary
#====================================================================

   def getWeight(self, inputs):

      battery_fac    = inputs['tech_bat']
      mission        = inputs['mission']

      # print(mission.__dict__);quit()
#====================================================================
# calculate energy capacity to install in battery including buffers
#====================================================================

#====================================================================
# energy-based prediction of battery weight
#====================================================================
      
      max_C          = 0.0
      Power_draw     = inputs['P_rated']
      Pmax           = numpy.amax(Power_draw)*0.001  # max power required from battery, kW
      Ereq           = inputs['E_stored']                    # in kWh
      E_segment      = inputs['segment_energy']
      Etotal         = Ereq*self.en_scaling

#====================================================================
# track energy used in normal mission profile separately
#====================================================================
   
      E_operations   = 0.0

      for i in range(mission.nseg):
         segment                 = mission.segment[i]
         self.C_rating[i]        = Power_draw[i]/Etotal*0.001     # convert power to kW, use to calculate C-rating

#====================================================================
# perform piecewise C-rating corrections for energy capacity
# and recalculate new C-ratings
#====================================================================

      Etotal                     = 0.0
      for i in range(mission.nseg):
         segment                 = mission.segment[i]
         correction              = 1.0 #+ 0.0542*segment.c_rating            # only apply 50% of correction (because we'll overshoot again)
         Etotal                  = Etotal       + correction*E_segment[i]  # kW-hr
         E_operations            = E_operations + correction*E_segment[i]*segment.op_cost_factor

      self.E_operations          = E_operations                                   # energy used for normal mission
      Etotal                     = Etotal * self.en_scaling

#====================================================================
# update C-ratings, calculate cell mass
#====================================================================

      for i in range(mission.nseg):
         segment                 = mission.segment[i]
         self.C_rating[i]        = Power_draw[i]/Etotal*0.001
         max_C                   = max(max_C,self.C_rating[i])

      self.max_c_rating          = max_C 
      m_battery                  = self.energy_to_mass(Etotal)

#====================================================================
# Step 7 in spreadsheet: check for max C-rating < 2
# if yes, size by energy requirements; already done!
# if max C-rating > 2, size by power  requirements, proceed further
#====================================================================

      if(max_C > 2 and not self.en_sizing):

#====================================================================
# Step 8: find average voltage during discharge
#====================================================================

         Iavg        = 900             # current draw at max power, Amps
         Vbatt       = Pmax*1000/Iavg  # nominal battery potential, V
         Vcell       = 3.6             # nominal cell    potential, V

#====================================================================
# Step 9: set # of serial cells based on battery V reqd and cell V
#====================================================================

         xS          = int(Vbatt/Vcell) + 1 

#====================================================================
# Step 10: calculate max cell current from max C-rating and cell cap
#====================================================================

         cell_maxC   = 20.0/3                     # max C-rating
         cell_Ah     = 3.0*self.energy_buf        # cell Amp-hr rating, adjusted for SOH
         Imax_cell   = cell_maxC*cell_Ah

#====================================================================
# Step 11: find # of parallel cells required
#====================================================================
   
         yP          = int(Iavg/Imax_cell)+1

#====================================================================
# Step 12: calculate battery energy in cells
#====================================================================

         cell_E      = xS*yP*Vcell*cell_Ah*0.001      # in kWh

         ratio       = cell_E/Etotal                  # ratio of avail. to reqd. energy
#         print('stored energy is ',cell_E, Etotal)
#====================================================================
# Step 13: increase parallel cell count to match energy req.
#====================================================================

         if ratio < 1:                                # less than reqd
            yP_new      = int(yP/ratio)+1
            # print(ratio,yP,'-->',yP_new)
            Etotal      = cell_E*yP_new/yP 
            yP          = yP_new 
            m_battery   = self.energy_to_mass(Etotal)
#            print('ratio is ',ratio,'adjusted')
         else:
            Etotal      = cell_E
#         print('cell energy build-up is ',xS*yP*cell_Ah*Vcell*0.001)
#         print('cell volume build up is ',xS*yP*0.01708*self.VIF,'liters')
#         print('cell volume from energy=',Etotal*1000/self.en_density*self.VIF/self.energy_buf,'liters')
         
#====================================================================
# Step 14: calculate cell C-rating for each phase and temp ramp rate
#====================================================================
         
         Tinit                   = mission.segment[0].temperature
         safe_T                  = False                 # get it? SAFETY = SAFE_T!
         dTtotal_max             = self.Tmax - Tinit     # max allowed temp increase

         while not safe_T:
            twavg_C                 = 0         
            T                       = Tinit
            Etotal                  = Vbatt*yP*cell_Ah*0.001
            for i in range(mission.nseg):
               segment              = mission.segment[i]
               C                    = Power_draw[i]/Etotal

#====================================================================
# Step 15: calcuate cell temperature rate of increase @ present C-rating 
#          and also final cell temperature at end of mission
#====================================================================

               dTdt                 = C*(self.dTdt[0] + C*(self.dTdt[1] + C*self.dTdt[2]))
               T                    = T + dTdt*segment.time
               self.cell_temp[i]    = T               # in Celcius
            dTtotal                 = T - Tinit 
            # print('temp,',T)
#====================================================================
# Step 16: check if final cell temperature < 70deg C; increase yP
# Step 17: increase yP if required 
# Step 18: done above, find battery energy for new yP value
#====================================================================
            
            if dTtotal > dTtotal_max: 
               dyP                  = int(min(yP*(dTtotal/dTtotal_max - 1.0),10))+1
               yP                   = yP + dyP
            else:
               safe_T               = True
               Etotal               = xS*yP*Vcell*cell_Ah*0.001      # in kWh

#====================================================================
# Step 19: find total battery energy: done above as Etotal
#====================================================================

#====================================================================
# Step 20: find cell mass and
# Step 21: apply pack integration factor for casing, etc.
#====================================================================

         m_battery                  = self.energy_to_mass(Etotal)

#====================================================================
# Step 22: volume of all cells ; done based on energy + buffer
#====================================================================
         
#         V_battery                  = xS*yP*0.01708*self.VIF

#====================================================================
# End of operations ; calculate mass and volume
# assume technology factor applies to both battery mass and volume
#====================================================================

      m_battery         = m_battery * battery_fac
      V_battery         = Etotal*1000.0*self.VIF/(self.en_density)*battery_fac

# remember relevant values in this object
      self.Eins         = Etotal                    # total Energy, kWh
      self.Vbatt        = V_battery*0.001           # battery volume, cu.m
      self.m_batt       = m_battery 

# return a dictionary 
      out_dict          = {'total': 0.0}
      return out_dict

#====================================================================
# calculate battery mass from rated cell specific energy (Wh/kg) and
# pack integ. factor
#====================================================================

   def energy_to_mass(self, Etotal):
      return (Etotal * 1000/self.sp_energy*self.IF)

#====================================================================
# end of operations
#====================================================================