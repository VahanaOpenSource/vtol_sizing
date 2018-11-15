!====================================================================
! some universal constants
!====================================================================

module precision
    integer, parameter          :: rdp      = 8       ! double precision
    integer, parameter          :: Nseg     = 40      ! spanwise segments
    integer, parameter          :: mxafoil  = 1       ! airfoil in a blade
    integer, parameter          :: NRPM     = 24      ! number of RPMs for RPM sweep (operation mode = 2)
    integer, parameter          :: ntaper   = 3      ! 1 to 3 
    integer, parameter          :: nthx     = 5      ! 8 to 20
    integer, parameter          :: nthtip   = 7      ! 10 to 35
    integer, parameter          :: nx       = 3      ! bilinear twist locations
    integer, parameter          :: mxOmega  = 6      ! 0.5 to 1 , steps of 0.1
    integer, parameter          :: nfc      = 8     ! max flight conditions

    real(kind=rdp), parameter   ::  ZERO    = 0.0d0
    real(kind=rdp), parameter   ::  one     = 1.0d0
    real(kind=rdp), parameter   ::  HALF    = 0.5d0
    real(kind=rdp), parameter   ::  TWO     = 2.0d0
    real(kind=rdp), parameter   ::  THREE   = 3.0d0
    real(kind=rdp), parameter   ::  pi      = 4.0d0*atan(1.0d0)
    real(kind=rdp), parameter   ::  d2r     = pi/180.0d0
    real(kind=rdp), parameter   ::  r2d     = 180.0d0/pi
    real(kind=rdp), parameter   ::  twobypi = two/pi

end module precision
