from DC_motor import DC_motor 

#====================================================================
# mechanical transmission weight: depends on Power
# it consists of only one component: the transmission!
#
# input:   P in kilowatts 
# output:  transmission mass in kgs
#====================================================================

def mechanical_transmission(P_engine):

    lb2kg       = 1/2.2                       
    P_hp        = P_engine/0.746              # power rating in Hp
    w_trans     = 0.4857*(P_hp**0.9472)       # in lbs
    mass        = w_trans*lb2kg               # in kg

    return mass 

#====================================================================
# electric transmission for fuel-based powerplants: depends on Power
# it consists of two parts
# (a) generator => converts engine power to electricity
# (b) motor     => converts  electricity to rotor power
# Note: wire weight is not counted here!
# input:  P in kilowatts, where P is power output from powerplant
#====================================================================

def electric_transmission(P):

    eta_gen     = 0.93            # generator efficiency

#====================================================================
# find generator weight (including its speed controller)
#====================================================================

    P_gen       = P               # generator rated power, kW
    m_gen       = DC_motor(P_gen)

#====================================================================
# find motor weight including its speed controller
#====================================================================

    P_motor     = P_gen*eta_gen   # motor power in kW
    m_motor     = DC_motor(P_motor)

    mass        = m_gen + m_motor 

    return mass