! ###################################################################
!
! Contains the routines and definitions related to the
! standard international localAtmosphere
!
! ###################################################################
#define type(x) TYPE(x), target

module isa
 
   implicit none

   real(kind=8), parameter :: gamma  = 1.4d0
   real(kind=8), parameter :: R_air  = 287.058d0
   real(kind=8), parameter :: rhoSL  = 1.2256d0
   real(kind=8), parameter :: tempSL = 273.15d0

   type :: atm
      real(kind=8) :: temp
      real(kind=8) :: delta
      real(kind=8) :: pres
      real(kind=8) :: sigma
      real(kind=8) :: altDens
      real(kind=8) :: rho
      real(kind=8) :: theta
      real(kind=8) :: speedSound
   end type! atm


contains
   
   !
   ! Note that altitude is in [meters] and deltaTempISA is in [degree C]
   !
   subroutine calculateDensityAltitude( deltaTempISA, altitude, localAtm)

      !
      ! variable declarations
      !
      real(kind=8), intent(in)  :: deltaTempIsa, altitude
      type(atm),    intent(out) :: localAtm

      real(kind=8), parameter   :: m2f = 3.28083989501312d0 
      real(kind=8), parameter   :: f2m = 0.3048d0
!
!
! compute T at altitude (according to ISA conditions) in deg C
      localAtm%temp  = 15.d0-0.001981d0*altitude*m2f 
      
! true temperature at altitude
      localAtm%temp  = localAtm%temp + deltaTempISA

! compute p/p0 (delta)
      localAtm%delta = (1.d0-6.876e-6*altitude*m2f)**5.265d0
      
! compute pressure altitude (m)
      localAtm%pres  = (518.4d0/0.00357d0*(1.d0- (localAtm%delta**0.1903d0) ))*f2m 
      
! compute rho/rho0 (sigma)
      localAtm%sigma = 288.16d0/(localAtm%temp + 273.16d0)* &
      	(1.d0 - 0.001981d0 * localAtm%pres * m2f/288.16d0)**5.256d0 
      
! compute height for this density ratio, aka, dens altitude
      localAtm%altDens  = (518.4d0/0.00357d0*(1.d0-localAtm%sigma**0.235d0))*f2m 
      
! compute equivalent air density based on density altitude
      localAtm%rho   = rhoSL*localAtm%sigma 
      
! Ratio of temperatures
      localAtm%theta = (localAtm%temp + tempSL)/(15.d0 + tempSL) 
      
! speed of sound at this altitude
      localAtm%speedSound = sqrt(gamma * R_air*(localAtm%temp + tempSL))  ! SI (m/s)

   end subroutine

end module isa
! ###################################################################
! END OF FILE
! ###################################################################
