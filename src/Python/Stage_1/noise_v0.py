#===============================================================================
# Compute noise from an individual fan/propeller/rotor
#===============================================================================
# Based on propeller noise from "A Review of Aerodynamic
# Noise from Propellers, Rotors, and Lift Fans" JPL Tech report 32-1462
# Inputs
#  T - thrust [N]
#  P - power [kW]
#  omega - fan rotational speed [rad/s]
#  geom - fan geometry structure
#  air - air property structure
#  S - distance from observer location to noise source [m]
#  theta - angle from thrust direction to observer location [rad]
#
# Outputs:
#   spl - Sound pressure level in dB
#
#===============================================================================

# Constants
NtoLb = 0.224809           # Newtons to lb
WtoHp = 0.00134102         # Watts to hp
f2m   = 0.3048             # feet to meters
Pref  = 0.0002             # dynes/cm^2
k     = 6.1e-27            # proportionality constant

def fanNoise(T,P,omega,geom,air,S,theta):

#===============================================================================
# Compute thrust coefficient
#===============================================================================

    Ct      = T ./ (air.rho * pi * geom.R^4 * omega.^2)

#===============================================================================
# Average blade CL
#===============================================================================

    CL      = 6 * Ct / geom.solidity

#===============================================================================
# Thrust in pounds
#===============================================================================

    T_lb    = T * NtoLb

#===============================================================================
# Parameters
#===============================================================================

    VTip    = omega * geom.R        # tip speed, m/s
    S_ft    = S / fttom             # Distance to source in ft
    R_ft    = geom.R / fttom        # Radius in ft
    A_ft2   = pi * R_ft*R_ft        # Disc area in ft^2
    P_Hp    = P * WtoHp             # Power in Hp
    M_t     = VTip / air.a          # tip mach number
    A_b     = A_ft2 * geom.solidity # blade area, sq. ft

#===============================================================================
# Compute rotational noise
#===============================================================================

    m_max   = 10                    # Maximum harmonic number
    p_m     = numpy.zeros(m_max)
    for m in xrange(1,m_max+1):
        p_m(m) = 169.3 * m * geom.B * R_ft * M_t ./ (S_ft * A_ft2) .* ...
            (0.76 * P_Hp / M_t^2 - T_lb * cos(theta)) .* ...
            besselj(m * geom.B, 0.8 * M_t * m * geom.B * sin(theta))

#===============================================================================
# Compute total RMS sound pressure level in dynes/cm^2 = 0.1 Pa
#===============================================================================

    p = sqrt(sum(p_m.^2))

#===============================================================================
# Convert to SPL in dB's
#===============================================================================

    SPL_rotational = 20 * log10(p / Pref)

#===============================================================================
# Vortex Noise (pg. 11)
# SPL = 10 * log(k * A_b * V_0.7^6 / 10^-16) + 20 * log(CL/0.4) [dB at 300 ft]
# where
# k = constant of proportionality = 6.1e-27
# A_b = propeller blade area [ft^2]
# V_0.7 = velocity at 0.7 * radius
#===============================================================================

    SPL_vortex  = 10 * log10(k * A_b * (0.7 * VTip).^6 / 1e-16) + 20 * log10(CL / 0.4)
# Total noise
# Adding atmospheric attenuation of 6 dB for every doubling of distance
# from reference
    spl = 10 * log10( 10**(SPL_rotational / 10) + 10**(SPL_vortex / 10)) - 6 * log2(S_ft/300)

    return spl
