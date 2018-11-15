#====================================================================
# obtain the weight of a turboshaft engine, given its max installed 
# power "P" in kilowatts. output is a dictionary 
#====================================================================

def piston_engine(P):       # Power in Kwatts, mass_fuel in kgs

#====================================================================
# get engine weight based on curve fits
#====================================================================

  lb2kg       = 1.0/2.2                     # conversion from lb to kg
  P_hp        = P/0.746                     # power (installed) in Hp      

  if P_hp <= 4:
    w_engine  = 14.0                        # up to 5 Hp: 14 lb EL-005 engine
  elif P_hp <= 20:
    w_engine  = 14.0/4.0*P_hp 
  else:
    w_engine  = 2.3668* (P_hp**0.9155)      # engine weight, lbs

  m_engine    = w_engine*lb2kg

#====================================================================
# use total fuel weight to calculate fuel system handling weight 
#====================================================================
      
  return m_engine # weight  dictionary [kg]

#====================================================================
# get SFC given flight and atmospheric conditions
# powerReq = engine output required, Pmax = max installed; units=kW
#====================================================================

def getSFC(theta, delta, powerReq, Pmax, KT, KD):  

#====================================================================
# calculate sfc base value based on data fits
# base sfc changes with power output of engine - less efficient at lower range
#====================================================================
  
  P_hp          = Pmax/0.746
  
  if P_hp <= 4.e0:
     sfc_base   = 0.42e0
  elif P_hp < 56.e0:
     sfc_base   = -0.0046*P_hp + 0.5935
  else:
     sfc_base   = 0.5185*(P_hp**(-0.09717e0))  # fuel consumption, lb/hp-hr


  if sfc_base < 0.3e0:
    sfc_base    = 0.3e0

  if sfc_base > 1.5e0: 
    sfc_base  = 1.5e0

#====================================================================
# get SFC scaling when operating at non-optimal (sub-max) conditions
#====================================================================
      
  x          = powerReq/Pmax                      # power ratio
  SFCratio   = 0.9526*(x**(-0.256))               # sfc scaling with power rating
#
  if SFCratio < 1.0:
     SFCratio    = 1.0

  sfc_base   = sfc_base*SFCratio

#====================================================================
# scalable engine data (ESF = engine scaling factor) to correct for 
# sfc variation with altitude and temperature for a given power 
# output - is it double counting?
#====================================================================

#  ESF        = (1.0 - KT*(theta-1.0))*(1.0 + KD*(delta-1.0))
#  sfc_corr   = (-0.00932*ESF*ESF + 0.865*ESF + 0.4450)/(ESF+0.3010)

#====================================================================
# after temperature and altitude corrections
#====================================================================

#  sfc        = sfc_corr*sfc_base
  
  sfc        = sfc_base 
  
#====================================================================
# converting from lb/hp-hr to kg/kW-hr
#====================================================================

  sfc        = sfc*0.45359/0.7457

  return sfc
