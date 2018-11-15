!====================================================================
! Fortran subroutine to calculate rotor power from momentum theory 
! severely limited.. fail, even
!====================================================================

   subroutine rotorPower(forwardSpeed, alpha, thrust, OmegaR, rho,   &
                         R, sigma, cd0, powerTotal)

!   use currentvalues
   implicit none

!====================================================================
! Input/output
!====================================================================

   real(kind=8), intent(in)  :: forwardSpeed, alpha, thrust, OmegaR,    &
                                rho, R, sigma, cd0
   real(kind=8), intent(out) :: powerTotal

!====================================================================
! Local variables
!====================================================================

   integer                   :: flag
   real(kind=8), parameter   :: ipf   = 1.15
   real(kind=8), parameter   :: pi    = 4.0d0*atan(1.0d0)
   real(kind=8)              :: lambda0,lambda,err,lambdaInduced
   real(kind=8)              :: powerInduced,ct,fProfile,mu
   real(kind=8)              :: powerProfile,cdeff,  ipfcr, ipf2
   real(kind=8)              :: area

!====================================================================
! Begin executable code
!====================================================================

   area     = pi*R*R 

!====================================================================
! compute ct and mu
!====================================================================

   ct       = thrust/(rho * area * OmegaR**2)
   mu       = forwardSpeed*cos(alpha)/OmegaR

!====================================================================
! Induced power required (fixed point iteration)
!====================================================================

   lambda0  = sqrt(0.5*ct)
   flag     = 1
   do while (flag == 1)
      lambda = mu*tan(alpha) + 0.5*ct/sqrt(lambda0**2 + mu**2)
      err    = abs(lambda-lambda0)/lambda

      if (err < 0.0005) flag = 0
      
      lambda0 = lambda0*0.05d0+lambda*0.95d0
   end do

!====================================================================
! compute induced inflow through disk
!====================================================================

   lambdaInduced = lambda - mu*tan(alpha)

!====================================================================
!Compute induced power factor (Johnson NDARC paper: curve fit by AS)
!====================================================================

   if (mu .le. 0.2) then               !linear region
      ipfcr       = 1.125d0 + mu
   elseif (mu .le. 0.5) then          !polyfit
      ipfcr       = 1.125d0 + 0.2d0 + 0.9*(mu-0.2d0) + 30.0d0*(mu-0.2d0)**2.0
   else 
      ipfcr       = 4.3d0       ! cutoff limiter for power factor
   end if 

!====================================================================
! for negative induced inflow (e.g. descent), set ind. power
! factor to 1.0, i.e. no power "gains"!
!====================================================================

   if (lambdaInduced .le. 1.d-3) then 
      ipfcr       = 1.0d0/ipfcr           ! absorb only 1/K of power in auto-rotation
   end if 

!====================================================================
! evaluate induced power
!====================================================================

   powerInduced  = ipfcr*thrust*lambdaInduced*OmegaR

!====================================================================
!"Propulsive" power (also has induced power factor)
!      write(*,*) ipfcr
!====================================================================

   if(alpha .lt. 0.d0) then 
      ipf2        = 1.0d0/ipf
   else 
      ipf2        = ipf 
   end if
!   write(*,*) powerInduced
   powerInduced  = powerInduced + ipf2*thrust*mu*OmegaR*tan(alpha)
    
!====================================================================
! THIS IS ONLY POWER NEEEDED TO SPIN THE BLADES: HUB X FORCE IS DONE 
! SEPARATELY
!====================================================================

   cdeff          = cd0*(1.d0 + 1.55d0*mu*mu)
   if (mu .gt. 0.3d0) then
       cdeff      = cdeff + (mu - 0.3d0)/0.38d0*1.25d0*0.03d0   ! extrapolation
   end if 

   fProfile      = sigma*cdeff/8.d0
   powerProfile  = rho * area * OmegaR**3 * fProfile

!====================================================================
! total power required for the rotor
!====================================================================

!   write(*,*) powerInduced,powerProfile
   powerTotal    = powerInduced + powerProfile

!====================================================================
! End of operations
!====================================================================

   end subroutine rotorPower