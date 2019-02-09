#====================================================================
# obtain the weight of a turboshaft engine, given its max installed 
# power "P" in kilowatts. output is a dictionary 
#====================================================================

def turbo_engine(P):       # Power in Kwatts, mass_fuel in kgs

#====================================================================
# get engine weight based on curve fits
#====================================================================

  lb2kg       = 1.0/2.2                     # conversion from lb to kg
  P_hp        = P/0.746                     # power (installed) in Hp      
  w_engine    = 7.3874* (P_hp**0.5519)      # engine weight, lbs

  m_engine    = w_engine*lb2kg

#====================================================================
# use total fuel weight to calculate fuel system handling weight 
#====================================================================
      
#  w_fuel      = mass_fuel * 2.2               # in lbs
#  w_tank      = 454.6*(w_fuel/6.5)**(-0.0566) # fuel system wt, lbs
#  m_tank      = w_tank*lb2kg                  # fuel sys mass, kg
#  m_tank      = mass_fuel*0.5           

  return m_engine # mass, kg 

#====================================================================
# get SFC given flight and atmospheric conditions
# powerReq = engine output required, Pmax = max installed
#====================================================================

def getSFC(powerReq, Pmax, delta, theta):  

#====================================================================
# calculate sfc base value based on data fits
# base sfc changes with power output of engine - less efficient at lower range
#====================================================================
      
  if powerReq < 1.e-3:
     sfc_base    = 0.0e0
  else:
     sfc_base    = 1.5089*((Pmax)**(-0.161e0))  # fuel consumption, lb/hp-hr

     if sfc_base > 0.82:     # curve not valid for SFC > 0.65: truncate
        sfc_base = 0.82

#====================================================================
# get SFC scaling when operating at non-optimal (sub-max) conditions
#====================================================================
      
  x          = powerReq/Pmax                      # power ratio
  SFCratio   = 0.9526*(x**(-0.256))               # sfc scaling with power rating

  if SFCratio < 1.0:
     SFCratio    = 1.0

  sfc_base   = sfc_base*SFCratio

#====================================================================
# scalable engine data (ESF = engine scaling factor) to correct for 
# sfc variation with altitude and temperature for a given power 
# output - is it double counting?
#====================================================================

#VT eqn
  sfc_corr   = 1.0 -0.299*(delta-1) + 1.253*(theta-1)

#MG eqn
#  ESF        = (1.0 - KT*(theta-1.0))*(1.0 + KD*(delta-1.0))
#  sfc_corr   = (-0.00932*ESF*ESF + 0.865*ESF + 0.4450)/(ESF+0.3010)

#====================================================================
# after temperature and altitude corrections
#====================================================================

  sfc        = sfc_base#*sfc_corr
#  sfc        = sfc_base 
  
#====================================================================
# converting from lb/hp-hr to kg/kW-hr
#====================================================================

  sfc        = sfc*0.45359/0.7457

  return sfc
