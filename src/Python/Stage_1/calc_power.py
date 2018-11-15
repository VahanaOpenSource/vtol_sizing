#=========================================================================   
# Python function to call PRASADUM for power calculations for asymm.
# compound with wing, and store results in an output dictionary
#
#  Inputs: 
#        (a) : max airspeed, knots
#        (b) : dictionary with aircraft properties [outputs of driver_sweep]
#        (c) : location to save output files and images
#
#=========================================================================   
import math, numpy
import shutil, os, sys
from run_CSD               import run_CSD 
from set_PRASADUM_inputs   import set_csd_rotor
from csd_rotor             import set_CSD_inputs
from performance_defaults  import *

def calc_power(Aircraft, code_dir, save_path):

#=========================================================================   
# Copy inputs from template directory to input directory of CSD solver
#=========================================================================   

   tar_dir              = os.path.join(code_dir,'Inputs/')
   output_dir           = os.path.join(code_dir,'Outputs/')
   src_dir              = os.path.join(code_dir,'Inputs_samples/Rotor_only/HYDRA_v2/')

   if os.path.exists(tar_dir):
      shutil.rmtree(tar_dir)
   shutil.copytree(src_dir, tar_dir)

#=========================================================================   
# parse the information and create derived quantities 
#=========================================================================   

   Flight      = Aircraft['Flight']
   Vmax        = Flight['V']
   alt         = Flight['alt']

   Airframe    = Aircraft['Airframe']
   Blade       = Aircraft['Blade']
   Rotor       = Aircraft['Rotor']
   Wing        = Aircraft['Wing']

   Wing        =  wing_calcs(Wing, Airframe, Flight)
   Rotor       = rotor_calcs(Rotor, Flight)
   Blade       = blade_calcs(Blade, Rotor)

#=========================================================================   
# Flush output folder (create if reqd)
#=========================================================================   

   output_dict = {}

   if os.path.exists(save_path):
      shutil.rmtree(save_path)
   os.makedirs(save_path)

#=========================================================================   
# Easy to use shortcuts
#=========================================================================   

   eta_p             = 0.85
   rcout             = 0.15

   Rotor             = Aircraft['Rotor']
   Radius            = Rotor['Radius']
   Mtip              = Rotor['Mtip']
   Nb                = Rotor['Nb']

   dirs              = {'prasadum': code_dir }

#=========================================================================   
# Get constants 
#=========================================================================   

   Weight            = Airframe['Wt']     # lbs
   atype             = Airframe['atype']  # configuration description
   CDo               = Wing['CDo']
   K                 = Wing['K']
   Swing             = Wing['Area']    # sq ft
   CLw               = Wing['CL']
   Span              = Wing['Span']

#=========================================================================   
# Compute total flat-plate area (all sources)
#=========================================================================   

   fwing             = 2*Swing*CDo

#=========================================================================   
# Create airspeed loop
#=========================================================================   

   Varray            = [0]
   if Vmax <= 0:
      keep_extending = False
   else:
      keep_extending = True


   while keep_extending:
      Varray.append(Varray[-1]+10)
      if Varray[-1] >= Vmax:
         keep_extending = False

   #Varray           = [Vmax]
   print 'Airspeed loop (knots) is \n',Varray

#=========================================================================   
# Rotor RPM schedule
#=========================================================================   

   Vlim           = Mtip*340.44       # in m/s
   Vtip_h         = Rotor['Vtiph']    # in m/s

#=========================================================================   
# Write header
#=========================================================================   

   f                 =  open(save_path+'Powercurve.dat','w')
   f.write('#     Speed(knots)   Wing_Lift(lb)   Total_Drag(lb)   Hub_Drag(lb)    Prop_SHP     Rotor_SHP       lift_offset       th0           th1c           th1s \n' )

#=========================================================================   
# Calculate drag from wing and airframe
#=========================================================================   

   f_airframe        = Aircraft['Airframe']['f']
   f_total           = f_airframe + fwing
   rho               = Flight['rho']            # slug/cu.ft

#=========================================================================   
# find RPM in high speed cruise at Vmax
#=========================================================================   

   kts2mps           = 1.853*5.0/18.0          # multiply by this number to convert knots to m/s
   VSI               = Vmax*kts2mps            # cruise speed, m/s
   hover_omega       = Vtip_h/Rotor['Radius']  # rad/s

   if Vtip_h + VSI <= Vlim:
      cruise_omega   = hover_omega
   else:
      Vtip           = Mtip * 340.44 - VSI            # m/s
      cruise_omega   = Vtip / Radius                  # cruise rotor speed, rad/s

#=========================================================================   
# find delta Omega (rad/s) between hover and vcruise
# also find the break points to change RPM
#=========================================================================   

   domega            = hover_omega - cruise_omega
   Vbreak            = (Vlim - Vtip_h)/kts2mps

   if Vbreak < 160:
      Vbreak         = 160
   dV                = Vmax - Vbreak 

#=========================================================================   
# Loop over airspeeds
#=========================================================================   

   # print domega, hover_omega;quit()
   for Vkts in Varray:
      VSI            = Vkts * 1.853  * 5/18        # m/s
      Vinf           = VSI / 0.3048                # ft/s

#=========================================================================   
# initialize output list
#=========================================================================   

      output_list    = [Vkts]

#=========================================================================   
# Create rotor data input file from template
# NOTE: RPM CHANGES WITH AIRSPEED; HAS TO BE HERE!
#=========================================================================   

      if Vkts <= Vbreak:
         Rotor['omega'] = hover_omega
      else:
         Rotor['omega'] = hover_omega - domega*(Vkts-Vbreak)/dV
#=========================================================================   
# Create CSD input file for rotor and blade properties
#=========================================================================   
   
      set_CSD_inputs(Rotor, Blade, tar_dir, src_dir)

#=========================================================================   
# Calc. download factor
#=========================================================================   

      if Vkts <= 50:
         dwf         = 1 - 0.1*(Vkts-50)/50
      else:
         dwf         = 1.0

#=========================================================================   
# Wing and body aerodynamic loads
# WIng lift coefficient is constant assumption: same angle of attack when
# there is no interference. AT low speeds, dyn. pressure is so low it doesnt 
# matter anyway.. lift is quadratic with Vinf
#=========================================================================   

      dyn_pr         = 0.5*rho * Vinf*Vinf
      Drag           = dyn_pr  * f_total                 #drag, lbs
      L_wing         = dyn_pr  * Swing * CLw             #lift, lbs

#=========================================================================   
# Output lists entry # 2: Wing lift
#=========================================================================   

      output_list.append(L_wing)

#=========================================================================   
# Compute CSD targets: SI units.. don't ask
#=========================================================================   

      T_target       = dwf*(Weight - L_wing) /  2.2 * 9.8      # Newtons
      T_target       = T_target/Rotor['NR']
      if atype == 'symmetric':
         L_target    = 0.0
      elif atype == 'coaxial' or atype == 'asymmetric':
         L_target    = Rotor['loff']*T_target*Rotor['Radius']*Vkts/Vmax   # linear variation with airspeed
      else:
         print 'valid aircraft types are coaxial, asymmetric or symmetric'
         quit('unknown configuration: cannot calculate power')

      shaft_tilt     = 0.0

      L_target       = 0.0          # FOR SYMMETRIC WING..

#=========================================================================   
# Run CSD code and convert outputs to FPS and Hp
#=========================================================================   

#      print 'weight',Weight, 'wing lift',L_wing
#      print 'inputs are',T_target,L_target,VSI
      info, Power, Fx = run_CSD(T_target, L_target, VSI, 0.0, alt, dirs)

      Fx             = Fx     / 9.8 * 2.2          # lbs
      Power          = Power  / 746                # Hp

#multiply by #rotors
      Fx             = Fx     * Rotor['NR']
      Power          = Power  * Rotor['NR']
#=========================================================================   
# Save CSD outputs in subfolder
#=========================================================================   

      subfolder_name = 'V' + str(Vkts)
      path           = os.path.join(save_path,subfolder_name)

      shutil.copytree(output_dir, path)

#=========================================================================   
# Update total drag
#=========================================================================   
   
      Drag           = Drag + Fx
#      print Drag, Fx;quit()
#=========================================================================   
# Output lists entry # 3, 4: Total Drag and hub X-force
#=========================================================================   

      output_list.append(Drag)
      output_list.append(Fx)

#=========================================================================   
# Compute propeller power and add rotor power
#=========================================================================   

      Prop_power     = Drag * Vinf / eta_p / 550.00         # in hp
      Power          = Prop_power + Power                   # add rotor power

      lift_offset    = L_target/Radius/T_target/0.3048      # nondiml

#=========================================================================   
# Output list entry # 5, 6, 7: Propeller and rotor shaft power, lift offset
#=========================================================================   

      output_list.append(Prop_power)
      output_list.append(Power - Prop_power)
      output_list.append(lift_offset)

#=========================================================================   
# Read trim_rotor_controls.dat 
#=========================================================================   

      controls_file     = os.path.join(path,'trim_rotor_controls_deg.dat')         
      controls          = numpy.loadtxt(controls_file)

#=========================================================================   
# Output list entry # 8, 9, 10: Rotor coll, th1c, th1s
#=========================================================================   

      output_list.append(controls[1])
      output_list.append(controls[2])
      output_list.append(controls[3])

#=========================================================================   
# write to file
#=========================================================================   

      for item in output_list:
         f.write("%15.4f" % item)
      f.write("\n")

   f.close()

#=========================================================================   
# End of operations
#=========================================================================   

   return output_dict
