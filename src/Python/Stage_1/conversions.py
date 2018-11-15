from numpy import pi, sqrt, cos, sin

grav    = 9.80665;   # gravitation constant (m/s2)
rho0    = 1.2256;    # standard air density (kg/m3)
f2m     = 0.3048;    # conversion from feet to meter
m2f     = 1/f2m;     # conversion from meter to feet
hp2kw   = 0.7457;    # conversion from hp to KW
kw2hp   = 1/hp2kw;   # conversion from KW to hp
kts2mps = 0.5144;    # conversion from knots to m/s
mps2kts = 1/kts2mps; # conversion from m/s to knots
lb2kg   = 0.4536;    # conversion from lb to kg
kg2lb   = 1/lb2kg;   # conversion from kg to lb
mps2kph = 3.6;       # conversion from m/s to km/hr
hr2sec  = 3600;      # conversion from hour to sec
sec2hr  = 1/hr2sec;  # conversion from sec to hour
rpm2rps = 2*pi/60;   # conversion from RPM to rad/sec
in2m    = 0.0254;    # conversion from inch to m
m2in    = 1/in2m;    # conversion from m to inch
Nm2ftlb = 0.7365;    # conversion from N-m to ft-lb
N2lb    = 0.22481;   # conversion from Newton to lb
lb2N    = 1.0/N2lb   # conversion from lb to Newtons