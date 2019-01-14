import numpy
from FEM import *
from numpy import pi,log
E 		= 122.0e9 					# Young's modulus, Nm^2
rho 	= 1650.0					# density of CF, meters
tbyc 	= 0.168

#======================================================================
# theory solution for natural frequency
# fn = input frequency in Hz
#======================================================================

def theory_soln(r1,taper,L,MbyA,Masses,xposn,fn):

#	LHS, RHS, dRHS 	= coefficients(r1,taper,L,Masses,xposn)
	LHS, RHS, dRHS 	= coefficients_simple(r1,taper,L,Masses,xposn)
	omegan 	= fn*2*pi
	beta 	= MbyA/tbyc
	coef 	= E*pi*0.5/(omegan*omegan)*LHS-rho*pi*RHS 
	b 		= beta*RHS + dRHS
	t 		= b/coef
	# print("{:20.12f} {:20.12f} {:20.12f}".format(taper,RHS/dRHS,rho*pi*RHS/(E*pi*0.5/(omegan*omegan)*LHS)))
#calculate spar mass from thickness
	rbar 	= (1+taper)*0.5*r1
	mspar 	= 2*pi*rbar*t*rho 		# spar 
	Mspar 	= mspar*L 

	return t, Mspar

#======================================================================
# deflection coeffs
#======================================================================

#======================================================================
# function to extract RHS coefficients
#======================================================================

#def RHS(r1,taper):

#======================================================================
# LHS, RHS coefficients
#======================================================================

def coefficients(r1,taper,L,Masses,xposn):
	r2 		= taper*r1
	rprime 	= (r2-r1)/L 
	ratio 	= r1/r2 
	a 		= 1-ratio*ratio 
	b 		= 1-ratio 
	c 		= 1/(rprime**3)*log(1.0/ratio)
	rat 	= L/r1
	LHS 	= a*0.5*(rat*rat/rprime+2*rat/rprime**2+1.0/rprime**3) + 	\
			  b*(-2/rprime**3-2*rat/rprime**2) + 						\
			  c 

#======================================================================
# calculate "C" coefficients
#======================================================================

	cm2 	= -L/(2*rprime)-r1/(2*rprime**2)
	cm1 	= 1/(rprime**2)
	c0 		= L/(2*rprime*r1*r1) - 1.0/(2.0*r1*rprime**2)

#======================================================================
# calculate "B" coefficients
#======================================================================
	
	b1 		= c0 
	b0 		= cm2/(r1*rprime) - cm1*log(r1)/rprime 
	bm1 	= -cm2/rprime 
	bl 		= cm1/rprime 

#======================================================================
# calculate RHS coefficients (KE)
#======================================================================
	
	I 		= b0*b0*(r1*L+rprime*L*L*0.5)
	II 		= b1*b1*(r1*L**3/3.0+rprime*L**4/4.0)
	III 	= bl*bl/(4.0*rprime)*(r2*r2*(2.0*log(r2)**2-2*log(r2)+1) - \
								  r1*r1*(2.0*log(r1)**2-2*log(r1)+1))
	IV 		= bm1*bm1/rprime*log(r2/r1)
	V 		= 2*b0*b1*(r1*L*L*0.5+rprime*L**3/3.0)
	VI 		= 2*b0*bl/4/rprime*(r2**2*(2*log(r2)-1) - r1**2*(2*log(r1)-1))
	VII 	= 2*b0*bm1*L 
	VIII 	= 2*b1*bl/rprime**2*(r2**3/9.0*(3*log(r2)-1)-r1**3/9.0*(3*log(r1)-1) -r1*0.25*( r2**2*(2*log(r2)-1) - r1**2*(2*log(r1)-1) ))
	IX 		= b1*bm1*L*L 
	X 		= 2*bl*bm1/rprime*(r2*(log(r2)-1) - r1*(log(r1)-1) )

	RHS 	= I + II + III + IV + V + VI + VII + VIII + IX + X 

	dRHS 	= 0.0

	for Mk,xk in zip(Masses,xposn):
		w 		= beam_shape(b0,b1,bm1,bl,xk,r1,rprime)
		dRHS 	= dRHS + Mk*0.5*w*w

#	print(RHS,dRHS,RHS/dRHS)
	return LHS, RHS, dRHS

#======================================================================
# LHS evaluation: terms A,B,C
#======================================================================

def get_ABC(taper):

	A 		= (1+taper)**4/taper**2*0.5
	B 		= (1+taper)**3/(taper-1)**2*(0.5-1/taper+0.5/taper**2)*(-2.0)
	C 		=((1+taper)/(taper-1))**3*(log(taper)-1.5-0.5/taper**2+2.0/taper)

	return A+B+C


#======================================================================
# RHS coefficients: poly in logr1
#======================================================================

def RHS_coeffs(taper):

	b0bar, brbar, b1bar, blbar, bmbar = bbars(taper)

#first term
	I_C 	= b0bar*b0bar*0.5
#	I_L 	= 2*b0bar*brbar*0.5
#	I_Q 	= brbar*brbar*0.5

#second term
	II_C 	= b1bar*b1bar*(3*taper+1)/(12*(taper+1))*0.25

#third term
	const 	= blbar*blbar/(taper*taper-1)*0.25
	III_C 	= const*(taper*taper-1+2*taper*taper*(log(taper)**2-log(taper)))
#	III_L 	= const*(4*taper*taper*log(taper) - 2*taper*taper + 2)
#	III_Q 	= const*(2*taper*taper - 2)

#fourth term 
	IV_C 	= bmbar*bmbar*(taper+1)/(taper-1)*log(taper)*0.25

#fifth term 
	const 	= b1bar*(2*taper+1)/(6*(1+taper))
	V_C 	= const*b0bar 
#	V_L 	= const*brbar

#sixth term 
	const 	= taper*taper/(taper*taper-1)*log(taper)-0.5 
	VI_C 	= blbar*b0bar*const 
	# VI_L 	= blbar*(b0bar + brbar*const) 
#	VI_Q 	= blbar*brbar 

#seventh term 
	VII_C 	= bmbar*b0bar 
	# VII_L 	= bmbar*brbar 

#eighth term 
	const 	= b1bar*blbar/(taper-1)/(taper*taper-1)
	taper3 	= taper*taper*taper 
	logt 	= log(taper)
	VIII_C 	= const*(taper3*logt/3.0 - taper*taper*logt*0.5 + \
			  (taper*taper-1)*0.25 - (taper3-1)/9.0)
	# VIII_L 	= const*((taper3-1)/3.0 - (taper*taper-1)/2.0)

#ninth term 
	IX_C 	= b1bar*bmbar*0.25 

#tenth term 
	const 	= bmbar*blbar 
	X_C 	= const*(taper*log(taper)/(taper-1) - 1)
	# X_L 	= const

	RHS_C   = I_C + II_C + III_C + IV_C + V_C + VI_C + VII_C + VIII_C + IX_C + X_C
	# RHS_L   = I_L +      + III_L +        V_L + VI_L + VII_L + VIII_L + 	   X_L
#	RHS_Q   = I_Q + 	   III_Q +        	    VI_Q 

#	print('hi',taper,I_L,III_L,V_L,VI_L,VII_L,VIII_L,X_L)
	return RHS_C

#======================================================================
# LHS, RHS coefficients
#======================================================================

def coefficients_simple(r1,taper,L,Masses,xposn):

	r2 		= taper*r1
	rprime 	= (r2-r1)/L 
	ratio 	= r1/r2 

#======================================================================
# alternate formulation of LHS terms
#======================================================================

	cbar 	= r1/tbyc*(1+taper)				# mean chord
	AR 		= 2*L/cbar 						# wing aspect ratio 
	K 		= 2.0*tbyc/AR
	K3inv 	= 1.0/(K*K*K)
	taper2 	= taper*taper 
	taper3 	= taper*taper2

#cubic approx. for LHS terms sum
	ABC2	= 1.626090567083019 - 0.482597591190902*taper + \
			  2.119679981957482*taper2 - 0.595753854491954*taper3

	LHS 	= ABC2*K3inv

#======================================================================
# alternate formulation for deflection coefficients
#======================================================================

	logr1 	= log(r1)
	
#======================================================================
# v2.0 calculations for individual terms of RHS
#======================================================================

	scale 	= L*L*K3inv/(K*K)

#	RHS_C = RHS_coeffs(taper)
	RHS_C 	= 0.369176557170034*taper3 + 0.151779690454288*taper2 + \
  			  0.269000169945435*taper  + 0.048225660574691

	RHS_v3 	= RHS_C*scale

#======================================================================
# accummulate RHS
#======================================================================

	dRHS 	= 0.0
	b0bar, brbar, b1bar, blbar, bmbar = bbars(taper)
	b0 		= (b0bar+brbar*logr1)*K3inv
	b1 		= b1bar*0.5/L*K3inv 
	bm 		= bmbar*L*K*0.5*K3inv 
	bl 		= blbar*K3inv
	for Mk,xk in zip(Masses,xposn):

		r 		= r1 + rprime*xk
		w 		= b0 + b1*xk + bm/r+bl*log(r)
		dRHS 	= dRHS + Mk*0.5*w*w

	return LHS, RHS_v3, dRHS

#======================================================================
# bbar coefficients for RHS calculations
#======================================================================

def bbars(taper):
	temp 	= ((taper+1)/(taper-1))
	temp 	= temp*temp*temp
	b1bar 	= temp*(taper-1.0)*(taper-2.0)
	bmbar  	= taper*temp/(taper+1)#(taper+1)**2/(taper-1)**3
	blbar 	= temp 				  #((taper+1)/(taper-1))**3
	b0bar 	= -temp*0.5*taper
	brbar 	= -temp
	return b0bar, brbar, b1bar, blbar, bmbar 

#======================================================================
# calculate extra kinetic energy terms due to lumped masses
#======================================================================

def beam_shape(b0,b1,bm1,bl,xk,rroot,rprime):

	r 	= rroot+rprime*xk
	w 	= b0+b1*xk+bm1/r+bl*log(r)

	return w 

#======================================================================
# beam shape functions
#======================================================================

def get_H(x,le):
	H 		= numpy.zeros(4)
	Hpp		= numpy.zeros(4)
	linv 	= 1.0/le 
	lins 	= linv*linv 

	xbar 	= x*linv 

	xbarsq 	= xbar*xbar 
	xbarcu 	= xbar*xbarsq 
	H[0]    = 1.0 - 3*xbarsq + 2*xbarcu
	H[1] 	= le*(xbarcu-2*xbarsq+xbar)
	H[2] 	= 1.0 - H[0]
	H[3] 	= le*(xbarcu-xbarsq)

	Hpp[0]  = (12.0*xbar-6.0)*lins
	Hpp[1] 	= ( 6.0*xbar-4.0)*linv
	Hpp[2] 	= -Hpp[0]
	Hpp[3] 	= ( 6.0*xbar-2.0)*linv 

	return H, Hpp

#======================================================================
# function that calculates radius of cross-section
#======================================================================

def get_r(rroot,rprime,x):
	r 		= rroot + rprime*x
	return r 

#======================================================================
# function calculate structural mass and bending stiffness for 
# a given location
#======================================================================

def get_mEI(x,rroot,rprime,t,MbyA=0.0):

	r 		= get_r(rroot,rprime,x)
	mNS 	= MbyA*2.0/tbyc*r
	m 		= 2*pi*rho*r*t + mNS
	EI 		= pi*r*r*r*t*E
	return m, EI 

#======================================================================
# python function to calculate element mass and stiffness matrices for
# vertical bending for a tapered beam
#======================================================================

from gaussian import gaussian
def taper_element(le,t,rroot,L,xroot=0.0,taper=1.0,MbyA=0.0):

	M 		= numpy.zeros((4,4))
	K 		= numpy.zeros((4,4))
	NG 		= 6
	xg, wg 	= gaussian(NG)
	rtip 	= rroot*taper
	rprime 	= (rtip-rroot)/L
	x 		= xg*le 

	for ig in range(NG):
		H,Hpp 	= get_H(x[ig],le)

		xbeam 	= xroot + x[ig]
		m,EI 	= get_mEI(xbeam,rroot,rprime,t,MbyA)
		dx 		= wg[ig]*le
		for i in range(4):
			for j in range(4):
				M[i,j] 	= M[i,j] + dx*H[i]*H[j]*m
				K[i,j]  = K[i,j] + dx*Hpp[i]*Hpp[j]*EI 
	return M, K 
	
#======================================================================
# get FEM based natural frequency solution for a tapered beam
#======================================================================

def taper_beam(MbyA, L, rroot, t, taper, Masses, xposn):

#======================================================================
# generate element lengths and tag elements with tip masses
#======================================================================

	le,elem_id 	= generate_elements(xposn,L,0.05*L)
	Me 			= []
	Ke 			= [] 
	ie 			= -1
	im 			= 0
	le 			= numpy.asarray(le)
	ne 			= len(le)
	xleft 		= numpy.zeros(ne)
	suml 		= 0.0
	for i in range(ne):
		xleft[i] 	= suml 
		suml 		= suml + le[i]

#======================================================================
# loop over elements
#======================================================================

	for l in le:
		ie 		= ie+1

#======================================================================
# BUild mass, stiffness matrices
#======================================================================

		M,K 	= taper_element(l,t,rroot,L,xleft[ie],taper,MbyA)

#======================================================================
# add tip mass effect if this element is tagged as one with a lumped mass
#======================================================================

		if(ie in elem_id):
			M[2,2] 	= M[2,2] + Masses[im]
			im 		= im + 1

#======================================================================
# append to list
#======================================================================

		Me.append(M)
		Ke.append(K)

#======================================================================
# Assemble and solve
#======================================================================

	fn 			= nat_freq(Me,Ke)[0]/(2.0*numpy.pi) 		# nat frequency, Hz
	return fn