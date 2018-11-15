
def spar_properties(material):

  spar              = {}
  
  if(material == 'titanium'):
    spar['rho']     = 4500.0e0                 # material density, kg/cu.m
    spar['sigma_y'] =  880.0e6                 # allowed stress in tension, Pa
    spar['E']       =  117.0e9                 # Young's modulus, Pa     
    spar['G']       =   43.0e9                 # Shear modulus, Pa 
  elif(material == 'aluminum'):
    spar['rho']     = 2700.0e0                 # 
    spar['sigma_y'] =  350.0e6
    spar['E']       =   71.0e9                 # Young's modulus, Pa     
    spar['G']       =   27.3e9                 # Shear modulus, Pa 
  elif(material == 'uniaxial_carbon'):
    spar['rho']     = 1660.0e0                 # 
    spar['sigma_y'] =  450.0e6
    spar['E']       =  122.0e9                 # Young's modulus, Pa     
    spar['G']       =   27.3e9                 # Shear modulus, Pa 
  elif(material == '090_carbon'):
    spar['rho']     = 1600.0e0                 # 
    spar['sigma_y'] =  275.0e6
    spar['E']       =   70.0e9                 # Young's modulus, Pa     
    spar['G']       =   47.3e9                 # Shear modulus, Pa 
  else:
    quit('critical error: unknown spar material: must be titanium or aluminum')


  return spar