#====================================================================
# Wing defaults
#====================================================================
import numpy 
def wing_calcs(Wing, Airframe,Flight):

    AR                  = Wing['AR']  
    Wing['CDo']         = 0.012                                # profile drag coef
    Wing['K']           = 1/numpy.pi/AR/0.85                   # K in K CL^2
    Wing['CL']          = numpy.sqrt(Wing['CDo']/Wing['K'])    # CL

    rho                 = Flight['rho']                        # in lb/cu.ft
    Vcruise             = Flight['V']                          # in knots
    cruise_lift         = Airframe['Wt']*Wing['fw']            # lbs
    cruise_speed        = Vcruise*1.853*5.0/18.0/0.3048        # ft/s
    fsdyn               = 0.5*rho*cruise_speed**2              # dyn. pressure
    Wing['Area']        = cruise_lift/fsdyn/Wing['CL']         # sq.ft
    Wing['Span']        = numpy.sqrt(Wing['Area'] * AR)        # ft

    return Wing

#====================================================================
# Rotor property calculations
#====================================================================

def rotor_calcs(Rotor, Flight):

#tip speed in cruise for edgewise rotor
    Mtip                 = Rotor['Mtip']
    Vcruise              = Flight['V']                                # knots
    Vtip                 = Mtip * 340.44 - Vcruise * 1.853 * 5/18     # tip speed, m/s
    Vtip                 = Vtip / 0.3048                              # tip speed, ft/s
    Rotor['omega']       = Vtip / Rotor['Rft']                        # cruise rotor speed, rad/s
    Rotor['Radius']      = Rotor['Rft']*0.3048                        # in meters
    Rotor['omega_hover'] = Rotor['Vtiph']/Rotor['Radius']

    return Rotor

#====================================================================
# Blade property calculations
#====================================================================

def blade_calcs(Blade, Rotor):

#mean chord
    crd                  = Rotor['sigma']*Rotor['Radius']*numpy.pi/Rotor['Nb']   # calculate mean chord
    crd0                 = crd
    crd1                 = crd

    Blade['root_chord']  = crd0      # in m 
    Blade['tip_chord']   = crd1      # in m 

    return Blade 