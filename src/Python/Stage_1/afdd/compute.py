import alighting
import compute
import conversions
import drive
import empennage
import engine
import engine_support
import flight_controls
import fuel_system
import fuselage
import hydraulics
import icing
import rotor

def check():
    print "checked!"

# # ==============================================================
# #
# #                  [OTHER SYSTEMS AND EQUIPMENT]
# #
# # ==============================================================

# # Mass of common equipment
# mass_common =  50 + ... # flight control, cockpit controls
#                50 + ... # automatic flight control system (fly-by-wire)
#               100 + ... # auxillary power group
#               100 + ... # instruments group
#                 0 + ... # hydraulic equipment group
#               400 + ... # electrical group
#               300 + ... # avionics group
#               200 + ... # furnishings group
#                 0 + ... # environmental control group
#                 0 + ... # anti-icing group
#               100 + ... # load handling group (internal)
#               100 + ... # load handling group (external)
#                50      # unusuable fluids
           
# # Convert masses to (kg) ---------------------------------------
# mass_common = mass_common*lb2kg
           
# # ==============================================================
# #
# #                  [EVALUATE TOTAL EMPTY MASS]
# #
# # ==============================================================

# mass_empty_array = [mass_tiltrotor_wing_group mass_rotor_group...
# 	mass_empennage_group mass_fuselage_group mass_alighting_gear_group...
# 	mass_engine_group mass_airinduction_group mass_fuel_system_group...
# 	mass_drive_system_group mass_flight_control_group ...
# 	mass_hydraulic_group mass_anti_icing_group mass_common]
	
# mass_empty = sum(mass_empty_array)


# # ##############################################################
# #                    --- END OF FILE ---
# # ##############################################################
