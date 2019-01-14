!===========================================================================
!> This routine obtains a converged inflow solution for a given pitch 
!> angle and flt condition
!===========================================================================

subroutine converge_inflow(flt , Rotor, Thrust, Power, converged)

    use bemt_data_structures
    use bemt_types
    implicit none

!===========================================================================
! inputs
!===========================================================================

    type(rotor_def)  , intent(inout)   :: Rotor
    type(flt_def)    , intent(in)      :: flt 

!===========================================================================
! Outputs
!===========================================================================

    real(kind=rdp), intent(out)     :: Thrust, Power
    logical,        intent(out)     :: converged

!===========================================================================
! Local variables
!===========================================================================

    integer                         :: Nb, counter, k
    real(kind=rdp), dimension(NSeg) :: pitch, lam, lami, sigma, lam_new, Cfx, Cfz
    real(kind=rdp)                  :: R, omega, Vtip, A, lamc, CT, rho, a1, b1, c1, d1
    real(kind=rdp)                  :: err, F, Cl, Cd, temp, determ, phi, dr
    real(kind=rdp)                  :: CL_alpha, alpha, all_err(nseg)
    real(kind=rdp)                  :: nu, Reynolds, kelvin, spd_sound, Mach
    real(kind=rdp)                  :: relax1, relax2, cphi, sphi, dlam, coef(Nseg)
    real(kind=rdp)                  :: vres(Nseg), ainv, nuinv, adeg, s1, al1, Cl1

!===========================================================================
! begin executable code
!===========================================================================

    CL_alpha        = Rotor % CL_alpha
    Nb              = Rotor % Nb
    R               = Rotor % radius
    omega           = Rotor % omega 
    rho             = flt  % rho 
    sigma           = real(Nb) * Rotor % chord / (pi * R)

    kelvin          = 288.0d0 - flt % altitude * 0.0065d0
    spd_sound       = sqrt(1.4d0*287.05d0*kelvin)
    nu              = 1.458D-6*sqrt(kelvin)/(1.0D+0+110.0D0/kelvin)
    nu              = nu / flt % rho 
    ainv            = one/spd_sound
    nuinv           = one/nu
!    nu              = 1.827d-5/flt % rho 

!===========================================================================
! linear lifting limit
!===========================================================================

    al1             = 10.d0                 ! linear limit
    al1             = al1*d2r
    Cl1             = CL_alpha*al1          ! lift at al1 
    s1              = -(Cl1)/(0.5d0*pi-al1) ! CL_alpha in stall, per radian

!===========================================================================
! derived quantities
!===========================================================================

    A               = pi * R * R                            ! disk area
    Vtip            = omega * R                             ! tip speed 
    CT              = flt  % Thrust / (rho * A * Vtip * Vtip) ! thrust coeff target
    lamc            = flt  % Vc / Vtip                        ! climb ratio

    relax1          = 0.6d0
    relax2          = 0.5d0

!===========================================================================
! convert collective + twist to radians
!===========================================================================
    
    pitch           = (Rotor % th0 + Rotor % twist) * d2r

    do k = 1, Nseg 
        pitch(k)    = (Rotor % th0 + Rotor % twist(k) )*d2r
        coef(k)     =-real(Nb) * half * (one - Rotor%r(k)) / Rotor % r(k)
        Vres(k)     = sqrt(Rotor % r(k) * Rotor % r(k) + lamc*lamc)
    end do 

!===========================================================================
!First generate an approximate solution using uniform inflow initial guess
!===========================================================================

    lam             = 0.05d0 + lamc

    counter         = 0 
    err             = one; all_err = one

    do while (err .ge. 1.d-2 .and. counter .le. 10) 
        counter          = counter +1

        err              = zero 
!!$OMP PARALLEL DO DEFAULT SHARED PRIVATE(k,phi,f,temp,determ)
        do k = 1, Nseg

!===========================================================================
! if this station is converged, dont perform iterations
!===========================================================================

            if(all_err(k) .le. 1.d-4) then 

            else
                lami(k)      = lam(k) - lamc             ! induced inflow
                phi          = lami(k)/Rotor % r(k)      ! inflow angle phi

!===========================================================================
!Iteration expression for total inflow ratio distribution along the span
!===========================================================================
    
                temp         = sigma(k)* CL_alpha * 0.125d0
                determ       = temp * pitch(k) * Rotor % r(k) +                   &
                               (temp*half - lamc*half)**2
    
                lam_new(k)   = lamc*half - temp*half +                          &
                              sign(one,determ)*sqrt(abs(determ))

                dlam         = lam_new(k) - lam(k)
                all_err(k)   = abs(dlam)
                lam(k)       = lam(k) + relax1*dlam
                lami(k)      = lam(k) - lamc
            end if 

            err          = err + all_err(k)

        end do              
    end do

!===========================================================================
!Use approximate solution as a guess for the "right" solution (large angles, tables etc)
!===========================================================================

    err         =  one; all_err = one
    counter     =  0
    do while (err .ge. 1e-3 .and. counter .le. 50) 

        counter    = counter + 1    
        err        = zero 

!===========================================================================
! loop over stations
!===========================================================================

        do k = 1, Nseg

!===========================================================================
! if this station is converged, do nothing
!===========================================================================

            if (all_err(k) .le. 1.d-4) then 

            else
                phi        = atan(lam(k)/Rotor % r(k)) !inflow angle phi

!===========================================================================
! calculate tip loss factor for outboard stations
!===========================================================================
                
                temp       = coef(k)/phi
                if(temp .le. -5.d0) then
                    F      = one 
                else 
                    f          = exp(temp)
                    if (f .ge. 0.99d0) f = 0.99d0
                    F          = twobypi*acos(f)    !prandtl tip loss factor
                end if 

!===========================================================================
! find total velocity, mach number, Reynolds number and angle of attack
!===========================================================================

!                temp        = Rotor % r(k) * Rotor % r(k) + lam(k) * lam(k)
!                temp        = sqrt(temp)            ! nondiml resultant vel 
                temp        = vres(k)*Vtip
                Mach        = temp * ainv
                Reynolds    = temp*Rotor % chord(k)*nuinv 
                alpha       = (pitch(k) - phi)
                adeg        = alpha*r2d 

                Rotor % alpha(k)    = adeg 

!===========================================================================
! Find Cl, Cd from tables
!===========================================================================

                if(Rotor % use_tables) then        
                    call table_lookup( Rotor % r(k), adeg, Mach, Reynolds, Cl, Cd)
                else
                    Cl         = CL_alpha * (pitch(k) - phi)
                    Cd         = Rotor % Cdo + 0.65d0*alpha*alpha
                end if 

!===========================================================================
! get vertical and in-plane loading components
!===========================================================================

                temp       = one/vres(k) 
                sphi       = lam(k)*temp 
                cphi       = Rotor % r(k)*temp
                Cfz(k)     = Cl*cphi - Cd*sphi
                Cfx(k)     = Cd*cphi + Cl*sphi

!===========================================================================
!setting up quadratic in lam_i
!===========================================================================

                a1         = one - sigma(k) * Cfz(k) * 0.125d0 / (F*Rotor % r(k))
                b1         = lamc * (two*a1 - one)
                temp       = vres(k)*vres(k)!Rotor % r(k)* Rotor % r(k) + lamc * lamc
                c1         = temp * (a1-one)
                d1         = b1*b1 - 4.0d0*a1*c1

!===========================================================================
!solving for the physical root
!===========================================================================

                lami(k)    = half/a1 * (-b1+sqrt(abs(d1))*sign(one,d1))

                lam_new(k) = lamc + lami(k)
                dlam       = lam_new(k) - lam(k)
                lam(k)     = lam(k) + relax2*dlam

                Rotor % inflow(k)   = lam(k)


                all_err(k) = abs(dlam)

            end if 
            err        = err + all_err(k)

!===========================================================================
!relaxation-type iteration
!===========================================================================

        end do 
    end do

!===========================================================================
! get thrust and power
!===========================================================================

    dr          = Rotor % r(2) - Rotor % r(1)

    Thrust      = zero; Power  = zero
    do k = 1, Nseg
        temp    = Rotor % r(k) * Rotor % r(k) + lam(k) * lam(k)     ! nondiml dyn pressure * 2

        Thrust  = Thrust + half * sigma(k) * temp * Cfz(k) * dr
        Power   = Power  + half * sigma(k) * temp * Cfx(k) * dr * Rotor % r(k)

        Rotor % dCTdr(k) = half * sigma(k) * temp * Cfz(k)
        Rotor % dCPdr(k) = half * sigma(k) * temp * Cfx(k) * Rotor % r(k)

    end do 
    
    Thrust      = Thrust * rho * A * Vtip * Vtip 
    Power       = Power  * rho * A * Vtip * Vtip * Vtip 

!===========================================================================
! if error is still significant, return the defaults
!===========================================================================

    converged   = .true.
    if(err .ge. 0.01d0) then
        Power       = 1.d16                       ! some value
        converged   = .false.
    end if 
    
!===========================================================================
! end of operations
!===========================================================================

return
end subroutine 
