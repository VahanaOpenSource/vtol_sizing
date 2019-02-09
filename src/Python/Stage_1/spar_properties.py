
def spar_properties(material):

  spar              = {}
  
# source: http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=mtp641
  if(material == 'titanium'):
    spar['rho']     = 4500.0e0                 # material density, kg/cu.m
    spar['sigma_y'] =  880.0e6                 # Yield stress,      Pa
    spar['E']       =  114.0e9                 # Young's modulus,   Pa     
    spar['G']       =   44.0e9                 # Shear modulus,     Pa 
    spar['tau_y']   =  760.0e6                 # ult shear stress,  Pa 

#source: http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma6061t6
  elif(material == 'aluminum'):
    spar['rho']     = 2700.0e0                 # density,         kg/cu.m
    spar['sigma_y'] =  276.0e6                 # yield stress,    Pa
    spar['E']       =   71.0e9                 # Young's modulus, Pa     
    spar['G']       =   27.3e9                 # Shear modulus,   Pa 
    spar['tau_y']   =  207.0e6                 # yield shear strs,Pa

#source: http://www.acpsales.com/upload/Mechanical-Properties-of-Carbon-Fiber-Composite-Materials.pdf
  elif(material == 'isolated_uniaxial_carbon'):
    spar['rho']     = 1660.0e0                 # density,         kg/cu.m
    # spar['sigma_y'] = 1200.0e6                 # ult. T/C stress, Pa 
    spar['sigma_y'] =  450.0e6                 # from Zach uni carbon below
    spar['E']       =  122.0e9                 # Young's modulus, Pa     
    spar['E2']      =   10.0e9
    spar['G']       =    5.0e9                 # Shear modulus,   Pa 
    spar['tau_y']   =   70.0e9                 # ult shear stress Pa
  elif(material == 'isolated_090_carbon'):
    spar['rho']     = 1600.0e0                 # density         ,kg/cu.m
    spar['sigma_y'] =  570.0e6                 # T/C limit stress,Pa
    spar['E']       =   70.0e9                 # Young's modulus, Pa     
    spar['E2']      =   70.0e9
    spar['G']       =    5.0e9                 # Shear modulus,   Pa 
    # spar['tau_y']   =   90.0e6                 # Shear stress,    Pa
    spar['tau_y']   =   47.0e6                 # zach bi carbon
  elif(material == 'isolated_pm45_carbon'):
    spar['rho']     = 1600.0e0                 # density, kg/cu.m
    spar['sigma_y'] =  110.0e6                 # ultimate stress in Tension/compression
#    spar['sigma_y'] =  275.0e6                 # from zach_bi_carbon
    spar['E']       =   17.0e9                 # E; Young's modulus, Pa
    spar['E2']      =   17.0e9                 # E; Young's modulus, Pa
    spar['G']       =   47.0e9                 # G; Shear modulus,   Pa (high-modulus version)
    spar['tau_y']   =  210.0e6                 # ult. shear stress,  Pa

#ZL values from aerodesigntool repo
  elif(material == 'uniaxial_carbon'):
    spar['rho']     = 1660.0e0                 # 
    spar['sigma_y'] =  450.0e6
    spar['E']       =  122.0e9                 # Young's modulus, Pa     
    spar['G']       =    5.0e9                 # Shear modulus, Pa 
  elif(material == '090_carbon'):
    spar['rho']     = 1660.0e0                 # 
    spar['sigma_y'] =  275.0e6                 # max tensile stress, Pa
    spar['tau_y']   =   47.0e6                 # max yield stress,   Pa
    spar['E']       =   70.0e9                 # Young's modulus,    Pa     
    spar['G']       =   33.3e9                 # Shear modulus,      Pa

  else:
    quit('critical error: unknown spar material: must be titanium or aluminum')

  return spar