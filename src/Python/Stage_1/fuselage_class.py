import numpy 

Nchannels   = 3               # number of parallel channels 
c1          = 0.17            # kg/m of signal wire ==> based on ALpha
c2          = 0.0057          # kg/kW/m of power wire, for (+) and (-)

#====================================================================
# Fuselage class
#====================================================================

class fuselage:

   def __init__(self, fuselage_dict, nseg):

      """
      function to initialize constants for the fuselage 
      Inputs:
      1. fuselage_dict: dictionary containing details: # rotors and lift fraction
      2. emp_dict:      dictionary containing geometry parameters for sizing
      3. nseg         : number of mission segments
      """

      self.rotor_group_id     = -1
      self.nrotors            = fuselage_dict['nrotors'][0]
      self.lift_frac          = fuselage_dict['liftfraction'][0]

      self.nrotor_groups      = 0
      self.rotor_thrust       = numpy.zeros(nseg)
      self.rotor_power        = numpy.zeros(nseg)
      self.motor_power        = numpy.zeros(nseg)
      self.rotor_torque       = numpy.zeros(nseg)
      self.rotor_aero_eta     = 0.75*numpy.ones(nseg)

# power carried by a drivetrain
      self.drivetrain_power   = numpy.zeros(nseg)

# masses and offsets of fuselage-mounted rotors 
      self.ycoords            = 0.0
      self.masses             = 0.0

#====================================================================
# geometry class that remembers where rotors and wings are 
# mainly used to find weight of wires
#====================================================================

class geometry:
   def __init__(self, geom_dict, nrotors, nwings):
      """
      function to initialize vehicle geometry-related parameters
      Inputs:
      1. geom_dict:  dictionary containing fuselage length and width, meters
                     (a) 

      2. nrotors  : total number of rotors in the system 
      3. nwings   : total number of wings  in the system

      """

      self.fuselage_length    = geom_dict['fuselage_length']
      self.fuselage_width     = geom_dict['fuselage_width']
      self.clearance          = geom_dict['clearance']

# rotor hub and wing root locations wrt cg
# X,Y,Z location used to calculate wire lengths for electric transmissions
      self.rotor_xyz          = numpy.zeros((3,nrotors))

# Z location of rotors wrt battery = fuselage height (worst case scenario for all rotors on high wings)
# here, assume height = width
      self.rotor_xyz[2,:]     = self.fuselage_width

# X location of rotors wrt battery = half of fuselage length
      self.rotor_xyz[0,:]     = self.fuselage_length*0.5

# group ID numbers for fuselages and wings: -1 => not attached, anything else => group number of component 
# that this rotor is mounted on 

      self.fuse_grp_id        =-numpy.ones(nrotors)
      self.wing_grp_id        =-numpy.ones(nrotors)
      self.rotor_grp_id       =-numpy.ones(nrotors)
      self.xmsn_grp_id        =-numpy.ones(nrotors)
      self.wire_length        = numpy.zeros(nrotors)

# only lateral (Y) coordinate decided by span-driven sizing for wing-mounted rotors, etc 

#====================================================================
# function to set rotor locations for wing- and fuselage-mounted rotors
#====================================================================

   def rotor_placement(self, rotor, wing, fuselage, transmission):      
      
      """
      This function determines the distance between individual rotor hubs and primary structures
      on which the rotors are mounted; also finds masses of signal wires/power cables
      For fuselage-mounted multi-rotors, the distance is 1.2 radii 
      For     wing-mounted multi-rotors, the distance is set by clearances
      
      This third coordinate is use to calculate wire length for elctric transmissions
      Inputs:
      1. rotor          : class containing details of various rotor groups 
      2. wing           : class containing details of various  wing groups
      3. fuselage       : class containing details of the fuselage 
      4. transmission   : class containing details of the various transmission groups

      Outputs:
      1. mass_sw        : mass of signal wires, kg 
      2. mass_pw        : mass of power  wires, kg
      """

      for rgid,rgroup in rotor.groups.items():

         xgid     = rgroup.xmsn_group_id
         irotor   = -1
         l_signal = self.fuselage_length

# for electric transmissions, proceed
         if(transmission.groups[xgid].type == 'electric'):

# if the rotor group is present on the fuselage, set coordinate to 1.2R for those rotors
# and tag it as being on the fuselage
            if(rgroup.fuse_group_id != -1):

               nrotors                        = fuselage.nrotors 
               for i in range(nrotors):
                  irotor                      = irotor + 1            # rotor number   
                  offset                      = 1.2*rgroup.radius
                  self.rotor_xyz[1,irotor]    = offset
                  self.fuse_grp_id[irotor]    = 0
                  self.xmsn_grp_id[irotor]    = xgid
                  self.rotor_grp_id[irotor]   = rgid

# remember mount offset and assembly mass in array
                  fuselage.ycoords            = offset
                  fuselage.masses             = rgroup.mass_assembly

# update signal wire length
                  self.wire_length[irotor]    = offset 
                  
# if the rotor group is also present on a wing, keep going
            wgids       = rgroup.wing_group_ids
            if(bool(wgids)):
               for wgid in wgids:
                  ycoords     = []
                  masses      = []
                  wgroup      = wing.groups[wgid]
                  nrotors     = int(wgroup.nrotors/wgroup.nwings/2)   # must be even number

# overhang length ahead of wing = one radius + 30% wing chord
                  L_mount     = rgroup.radius + wgroup.chord*0.3

# loop over duplicate wings in this group
                  for j in range(wgroup.nwings):

# assign coordinates for wing-mounted rotors
                     offset      = 0.5*self.fuselage_width + self.clearance*rgroup.radius

# add signal wire length from wing tip to wing tip
                     l_signal    = l_signal + wgroup.span
# rotors on one half-wing
                     for i in range(nrotors):
                        irotor                 = irotor + 1
                        y_rotor                = offset+rgroup.radius

# set coordinate and remember that this rotor is located on a wing
                        self.rotor_xyz[1,irotor]  = y_rotor
                        self.wing_grp_id[irotor]  = wgid
                        self.xmsn_grp_id[irotor]  = xgid
                        self.rotor_grp_id[irotor] = rgid
                        self.wire_length[irotor]  = y_rotor + L_mount

# add mirror images: rotors on other half-wing
                        self.rotor_xyz[1,irotor+nrotors]  = y_rotor
                        self.rotor_grp_id[irotor+nrotors] = xgid
                        self.xmsn_grp_id[irotor+nrotors]  = xgid
                        self.wing_grp_id[irotor+nrotors]  = wgid
                        self.wire_length[irotor+nrotors]  = y_rotor + L_mount

# remember array of rotor assembly coordinates
                        if(j == 0):
                           wgroup.ycoords[i]    = y_rotor
                           wgroup.masses[i]     = rgroup.mass_assembly
                           wgroup.rotor_T[i]    = numpy.amax(wgroup.rotor_thrust)
# update offset parameter for next rotor hub
                        offset                  = offset + rgroup.radius*(2.0 + self.clearance)

# skip "nrotor" indices: mirror images
                     irotor                  = irotor + nrotors

#====================================================================
# calculate weight of power wires and signal wires
#====================================================================
      
# masses in the fuselage
      P           = 0.0
      l_signal    = l_signal + numpy.sum(self.wire_length)
      mass_sw     = c1*Nchannels*l_signal
      Pdsum       = 0.0

# masses along mount structures 
      for i in range(rotor.nrotors):
         rgid     = self.rotor_grp_id[i]
         xgid     = self.xmsn_grp_id[i]
         tgroup   = transmission.groups[xgid]

# calculate rated line power for electric transmissions
         if(tgroup.type == 'electric'):
               P1 = transmission.groups[xgid].line_rated_power
               P  = P + P1

# find sigma power x length over all power calbes
         Pdsum    = Pdsum + self.wire_length[i]*P1

# add up power cables in fuselage
      Pdsum       = Pdsum + (self.fuselage_length+2*self.fuselage_width)*P      

# calculate mass of power cables
      mass_pw        = c2*Pdsum
      
      wires          = {'signal_wires': mass_sw, 'power_cables': mass_pw}
      return wires, P  