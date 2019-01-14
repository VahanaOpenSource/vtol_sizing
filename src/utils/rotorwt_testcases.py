from __future__ import print_function
import sys
sys.path.append('../Python/Stage_1')
from blade_wt_modelv2 import *

#============================================================
# UH-60 blade properties
#============================================================
R1 	 = {'chord': 0.53,
		'R'    : 8.17,
		'Omega': 28,
		'nz'   : 3.8,
		'Fz'   : 71272.0,
		'nb'   : 4,
		'nr'   : 1, 
		'mat'  : 'titanium',
		'betap': 11} 			# blade hinge removes root bending moment

UH60,mt = blade_wt_modelv2(R1, 1.0)
print(UH60['blades']/R1['nb']*2.2, 'lbs')

#============================================================
# Bo-105 blade properties
#============================================================

R1 	 = {'chord': 0.27,
		'R'    : 4.9,
		'Omega': 44.4,
		'nz'   : 3.8,
		'Fz'   : 24500.0,
		'nb'   : 4,
		'nr'   : 1, 
		'mat'  : 'uniaxial_carbon',
		'betap': 7} 					# extra coning due to normal load factor

Bo105,mt = blade_wt_modelv2(R1, 1.0)
print(Bo105['blades']/R1['nb'],' kgs')

#============================================================
# Vahana
#============================================================

R1 	 = {'chord': 0.082,
		'R'    : 0.75,
		'Omega': 300,
		'nz'   : 3.8,
		'Fz'   : 1530,
		'nb'   : 3,
		'nr'   : 1, 
		'mat'  : 'uniaxial_carbon'} 			# precone removes moment

Alpha,mt = blade_wt_modelv2(R1, 1.0)
print(Alpha['blades']/R1['nb'],' kgs')
