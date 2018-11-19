#====================================================================
# function to calculate area based weight of a fixed wing
# includes skin only; add extra for ribs by boosting min gauge
# to 1.5mm
#====================================================================

def skin_mass(chord, rho, span):

	t 			= min(0.005*chord,0.0015)
	perimeter 	= 2.1*chord
	skin_vol 	= perimeter*t*span
	mass_skin 	= rho*skin_vol 

	return mass_skin