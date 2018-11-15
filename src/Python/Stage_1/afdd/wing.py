#===============================================================
#            [TILT-ROTOR WING GROUP]
#===============================================================
# Section 29-1.1 NDARC v1_11, pg 245
# Original formulation from Chappell and Peyran
#===============================================================

from conversions import *

#===============================================================
# Empirical parameters
#===============================================================

#===============================================================
#Note: assuming 0/90 Carbon: G_tb same as 27 GPa, E_tb same as 70 GPa
#===============================================================

G_tb        	= 4e6*144.0 		# torque box  shear modulus, lb/sq.ft
E_tb        	= 10e6*144.0 		# torque box Youngs modulus, lb/sq.ft
E_sp        	= E_tb 				# wing spar modulus, lb/sq.ft
dens_tb       	= 3.1045        	# density, torsion box  (slug/ft3) ~ carbon fiber
dens_spar     	= 3.1045        	# density, spar         (slug/ft3) ~ carbon fiber
e_tb          	= 0.8           	# structural effeciency factor for torsion box
e_sp          	= 0.8           	# structural effeciency factor for spar
ct            	= 0.75          	# weight correction for spar taper (equiv. stiffness)
chord_tb_wing 	= 0.6 				# ratio of torque box chord to wing chord
omega_B       	= 0.5 				# wing flap bending freq normalized by rotor speed, /rev
omega_C       	= 0.8 				# wing chordwise bending freq. normalized by rotor speed, /rev
omega_T 		= 0.9 				# wing torsion freq. normalized by rotor speed, /rev
wing_tc       	= 0.23 				# wing thickness to chord ratio

#===============================================================
# other parameters
#===============================================================

loading_fair  	= 2.0				# Unit weight of LE/TE fairings (lb/sq.ft)
loading_flap  	= 3.0 				# Unit weight of control surface (lb/sq.ft)
f_fit         	= 0.12 				# weight fraction of the wing
f_fold        	= 0.0				# wing fold/tilt

#===============================================================
# Begin execution
#===============================================================

def tiltrotor_wing_wt(vparams):
	W     		= vparams['gtow'] 			# take-off weight, lbs
	Omega 		= vparams['Omega']			# rotor speed 	 , rad/s
	R 			= vparams['radius']			# rotor radius   , ft
	wtip 		= vparams['wing_tip_wt']	# wing tip weight, lbs
	T 			= vparams['design_T']		# design thrust  , lbs
	tilt_wing 	= vparams['tilt_wing'] 		# is it a tilting wing with tilt-rotor?
	rotor_wt    = vparams['wt_rotors'] 	    # in lbs; weight of all rotors
	fac 		= vparams['tech_factors'].wing 
	wings		= vparams['wing']
	nwing_grp  	= wings.ngroups

#===============================================================
# derived quantities
#===============================================================

	all_mass 	= {}
	grand_tally = 0.0
	for i in range(nwing_grp):
		wing 	= wings.groups[i]
		nwings  = wing.nwings

		b_w 	= wing.span*m2f
		S_w 	= wing.area*m2f*m2f
		c 		= wing.chord*m2f

		S_fair      = 0.10*S_w 					# fairing area, sq. ft
		S_flap      = 0.25*S_w 					# flap    area, sq. ft
		r_pylon     = 0.3*R  					# pylon rad. gyration, ft (to find pitch inertia of tip mass)
		mtip 		= wtip/32.2 				# tip mass, slugs

		f_tip       = wtip/W 					# tip weight ratio
		f_mode      = 1.0 - f_tip 				# mode shape correction based on tip weight ratio

		t_w         = wing_tc*c  				# wing thickness, ft
		c_tb        = chord_tb_wing*c 			# torque box chord, ft

		wtb       	= chord_tb_wing 			# short name for tb chord/wing chord ratio
		tauw 		= wing_tc 					# short name for wing t/c ratio

#===============================================================
# Form factors for wing sections
# B  = beam bending, C = chord bending, T = torsion, 
# VH = vertical/horizontal spar cap bending
#===============================================================

		F_B       	= 0.073*sin(2*pi*(tauw-0.151)/0.1365) +  0.14598*tauw + 	  \
					  0.610*sin(2*pi*(wtb+0.08)/2.1560)   - (0.4126-1.6309*tauw)* \
					  (wtb-0.131) + 0.0081

		F_C       	= 0.640424*wtb*wtb - 0.89717*wtb + 0.4615*tauw+0.655317

		F_T       	= ((0.27-tauw)/0.12)*0.12739*(-0.96+sqrt(3.32 + 94.6788*wtb - \
					   (wtb/0.08344)**2)) - 2.7545*wtb*wtb + 5.1799*wtb - 0.2683
			
		F_VH   	 	= 0.25*sin(5.236*wtb) + 0.325

#===============================================================
# Find wing torsion stiffness from assumed torsion frequency
# and find torque box CS area from stiffness, and finally mass
# Note: NDARC uses b_w - w_attach as effective length for torsion
# however, we're using 90% of wing span b_w 
#===============================================================

		GJ        	= (omega_T*Omega)**2*0.25*(0.9*b_w)*mtip*r_pylon*r_pylon

		A_tb      	= 4*GJ/(G_tb*F_T*t_w*t_w)

		mass_box  	= A_tb*dens_tb*b_w/e_tb 		# in slugs
		wt_box 		= mass_box * 32.2 			    # in lbs 

#===============================================================
# Find wing bending stiffnesses from assumed mode frequencies 
# then calculate the torque box, spar bending stiffnesses
# proceed to area calculation, then volume and finally mass
#===============================================================

		EI_C      	= (omega_C*Omega)**2*(1/48.0)*b_w**3*mtip*f_mode
		EI_B      	= (omega_B*Omega)**2*(1/48.0)*b_w**3*mtip*f_mode

		EI_Ctb    	= E_tb*F_C*A_tb*c_tb*c_tb*0.25

		EI_Csp    	= EI_C - EI_Ctb 
		EI_Csp = 0 if(EI_Csp < 0) else EI_Csp
		A_Csp       = EI_Csp/(E_sp*c_tb*c_tb*0.25)

		EI_Btb      = E_tb* F_B * A_tb * t_w*t_w * 0.25
		EI_VH     	= E_sp*F_VH * A_Csp* t_w*t_w * 0.25

		EI_Bsp    	= EI_B - EI_Btb - EI_VH
		EI_Bsp = 0 if(EI_Bsp < 0) else EI_Bsp

		A_Bsp     	= EI_Bsp/(E_sp*t_w*t_w*0.25)

		A_sp      	= A_Csp + A_Bsp#if(A_sp<0) A_sp=0end
		mass_spar 	= ct*A_sp*dens_spar*b_w/e_sp 		# slugs
		mass_prim 	= (mass_box + mass_spar) 			# slugs 

		wt_spar 	= mass_spar * 32.2  				# lbs
		wt_prim 	= mass_prim * 32.2  				# lbs 

#===============================================================
# Aerodynamic surfaces: fairings and flaps
#===============================================================

		wt_fair 	= S_fair*loading_fair 		# lbs
		wt_flap 	= S_flap*loading_flap 		# lbs

		mass_flap 	= wt_flap / 32.2 			# slugs 
		mass_fair 	= wt_fair / 32.2 			# slugs 

#===============================================================
# fitting: in NDARC, its f_fit*T, here its a different expression
#===============================================================

		mass_fit  	= f_fit/(1.0-f_fit)*(mass_prim+mass_fair+mass_flap)
		wt_fit 		= mass_fit * 32.2 			# lbs

#===============================================================
# folding mechanisms: fraction of everything else
#===============================================================

		mass_fold 	= f_fold*(mass_prim+mass_fair+mass_flap+mass_fit)
		wt_fold 	= mass_fold*32.2 			# lbs 

#===============================================================
# add tilt mechanism weights for tilt-rotor/tilt wing: extra
# 15% (or 10% for tilt-wing) of all pieces that "tilt"
#===============================================================

		if tilt_wing:
			f_tilt   = 0.13
			wt_tilt  = f_tilt*(wt_prim+wt_fair+wt_flap+wt_fit+wtip)
		else: 							# tilt-rotor: 15% of tilting components
			f_tilt   = 0.13
			wt_tilt  = f_tilt*wtip

#===============================================================
# apply tech factor, calculate total weight in lbs
#===============================================================

		wt_prim 	= wt_prim*fac*nwings
		wt_fair 	= wt_fair*fac*nwings
		wt_flap 	= wt_flap*fac*nwings
		wt_fit  	= wt_fit *fac*nwings
		wt_fold 	= wt_fold*fac*nwings
		wt_tilt 	= wt_tilt*fac*nwings

		total       = wt_prim + wt_fair + wt_flap + wt_fit + wt_fold + \
					  wt_tilt

		grand_tally = grand_tally + total 	
#===============================================================
# tiltrotor wing mass (kg) dictionary with breakdown and total
#===============================================================

		key 			= 'group' + str(i)
		all_mass[key] 	= {'primary': wt_prim*lb2kg, 'fairing': wt_fair*lb2kg, 
					         'flaps': wt_flap*lb2kg, 'fitting':  wt_fit*lb2kg,
					          'tilt': wt_tilt*lb2kg,
					       'folding': wt_fold*lb2kg, 'total'  :   total*lb2kg}


	all_mass['total'] = grand_tally*lb2kg
	return all_mass