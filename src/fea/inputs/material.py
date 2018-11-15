#
# calculator for cross-sectional properties
# 
# x is the axial direction and y,z and the cross-sectional axes
#

import numpy as np

# ===================================================================
# solid circle
# ===================================================================

density  = 2700       # SI [kg/m3] Aluminum
length   = 1          # SI [m]
outerDia = 2*7.89206e-2 #10 * 1.e-2 # SI [m]
innerDia = 0.85*outerDia # SI [m]

#innerDia =  5 * 1.e-2 # SI [m]

r2 = 0.5*outerDia
r1 = 0.5*innerDia

area = np.pi*(r2*r2-r1*r1)

Ix = np.pi/2*(r2**4 - r1**4)
Iy = 0.5*Ix
Iz = Iy

J = np.pi/2*(r2**4 - r1**4)

print ' - Density         %10.4e kg/m3' % density
print ' - Length          %10.4e m'     % length
print ' - Outer diamater  %10.4e m'     % outerDia
print ' - Inner diamater  %10.4e m'     % innerDia
print ' - Area            %10.4e m2'    % area
print ' - Iy              %10.4e kg-m2' % Iy
print ' - Iz              %10.4e kg-m2' % Iz
#print ' - Ix              %10.4e kg-m2' % Ix
print ' - J               %10.4e kg-m2' % J
print ' - Max y dimension %10.4e m'     % r2
print ' - Max z dimension %10.4e m'     % r2






