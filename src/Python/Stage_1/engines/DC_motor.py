#====================================================================
# Brushless DC motor weight sizing: depends only on rated power
# For small scale motors,
# weight [g] = 0.2513 * max_power [W]
# For larger scale, a more gentle curve fit is applied
# note:valid up to 10 kW
#
# Each motor/generator requires a speed controller, which is essentially
# an AC/DC rectifier that also regulates the rotational speed. 
#
# The speed controller weight depends on the current, and is sized here
# assuming a 12 Volt setting. High tension lines are also possible but 
# present an operational safety issue, hence not considered in this version.
# Therefore, the motor sizing is only a function of the rated power
# and not anything else
#====================================================================

#====================================================================
# Input:  P = rated power (kW) - how much power the motor is handling
# Output: m =       MASS  (Kg)  
#====================================================================

def DC_motor(P):

#====================================================================
# small scale fit
#====================================================================

    if P <= 10.0:         
        m_motor = 0.2513 * P
        
#====================================================================
# electronic speed controller weight based on the current drawn and
# is assumed to be a linear relationship from the 8 lb QBiT design
# at small scale; full scale fits include the esc!
#====================================================================

        V           = 12.0
        I           = P / V            # in Kilo-Amperes, i.e. 10^3 * A

#====================================================================
# 8lb QBT draws 7.25 [A] and weighs 32 [g], i.e. 32 g per 7.25 Amps
# so 
# ESC weight in kg = 32/7.25 * (current in Kilo-Amperes)
#====================================================================

        m_esc      = I * 32.0 / 7.25

        m_total    = m_motor  + m_esc

#====================================================================
# full-scale data: the speed controller is included in the data fit
# fit is power in kW vs. mass in kg 
#====================================================================

    else:
        m_total    = 0.8373*P**0.726
        #print('rated power in kW is ',P)
    return m_total 