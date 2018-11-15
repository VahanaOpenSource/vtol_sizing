#====================================================================
#
# Nacelles and air induction group
#
# engine-exhaust system, engine accessories, engine support
# use for piston engine or turboshaft engine
#====================================================================

#====================================================================
# empirical parameters
#====================================================================

#f_lub    = 1.4799 # 1.0 if lubrication is included in engine weight
#k_0exh   = 0.0    # ??
#k_1exh   = 0.0    # ??
f_lub    = 1.0 
f_airind = 0.3
f_pylon  = 0.00     # verify with MGB
#s_nac    = 20.0     # because included in other wetted area (huh? AS)

from conversions import *
#====================================================================
# begin routine
#====================================================================

def engine_accessories(vehicle_parameters):

#====================================================================
# unpack input dictionary
#====================================================================

   nengine         = vehicle_parameters['nengine']
   gtow            = vehicle_parameters['gtow']
   power_installed = vehicle_parameters['pwr_installed']   
   powerplant_tot  = vehicle_parameters['engine_wt']
   fac             = vehicle_parameters['tech_factors'].powerplant
   s_nac           = 50.0/1600.0*power_installed 
   wght_eng        = 0.667*powerplant_tot/float(nengine)    # 2/3 of engine+accs+exhaust sys

#====================================================================
# engine accessories: assume its included in engine weight model
#====================================================================

#   wght_acc        = 2.088 * f_lub * (wght_eng/nengine)**0.5919 * nengine**0.7858

#====================================================================
# weight of engine support structure, engine cowling, pylon 
# support structure and air induction components
#====================================================================

   wght_supt   = ( 0.0412    * (1-f_airind) * 
                   wght_eng**1.1433 * nengine**1.3762 )
   
   wght_cowl   = 0.2315*s_nac**1.3476  
#   print 'check engine cowling weight',wght_cowl,' lbs'
   
   wght_pylon  = f_pylon * gtow
   
   wght_airind = 0.0412 * f_airind * wght_eng**1.1433 * nengine**1.3762

#====================================================================
# apply tech factors, find total weight,
# pack breakdown into a dictionary and write to disk
#====================================================================
   
   wght_supt   = wght_supt  *fac 
   wght_cowl   = wght_cowl  *fac 
   wght_pylon  = wght_pylon *fac 
   wght_airind = wght_airind*fac 
   
   total       = wght_supt + wght_cowl + wght_pylon + wght_airind #+ wght_acc

   engine_acc  = {'engine_support': wght_supt *lb2kg, 'cowling':wght_cowl  *lb2kg, 
                  'pylon_support' : wght_pylon*lb2kg, 'air_ind':wght_airind*lb2kg,
                  'total': total*lb2kg}
    
   return engine_acc



   #nengine = aircraft.engine.number
   #P       = aircraft.engine.ratedP
   #W       = aircraft.engine.weight

   #wght_eng = nengine * W

   #wght_exh = nengine * (k_0exh + k_1exh*P)

   #try:
   #    wght_acc = 2.088 * f_lub * (wght_eng / nengine)**0.5919 * nengine**0.7858
   #except:
   #    print "Warning: error in engine acc weight ", wght_eng, nengine
   #    wght_acc = 0.0

   #if(aircraft.engine.__class__.__name__ == 'Hybrid_Special414'):
   #    wght_acc = 0.0


   #if(not bfile is None):
   #    bfile.write("---- ENGINE ----\n")
   #    bfile.write("wght_eng   %9.2f\n"%(wght_eng ))
   #    bfile.write("wght_exh   %9.2f\n"%(wght_exh   ))
   #    bfile.write("wght_acc   %9.2f\n"%(wght_acc  )) 

   #return wght_eng + wght_exh + wght_acc
