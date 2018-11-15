#define type(x) Type(X), target
module airfoil_design_types

!=======================================================================
! material properties
! E       = Young's modulus 
! G       = Shear   modulus 
! sigma_y = yield stress 
! rho     = density
!=======================================================================

   type :: Material_Type_Def
    real(kind=8)    ::       E, G, sigma_y, rho

   end type 

!=======================================================================
! big-picture stuff (SI units)
! rotor chord, radius, rotor speed, load factor, vertical thrust, torque
!=======================================================================

   type :: Rotor_Info_Def


      real(kind=8) :: chord, R, Omega, nz, Fz, Mz 

   end type 

!=======================================================================
! spar properties
! t  = thickness, b = width, xNA = neutral axis location
! xc = chordwise location, E = YOung's modulus, G = shear modulus 
! sigma_y = yield stress, rho = density, m = mass/span in kg/m
!=======================================================================

   type :: Spar_Info_Def


      real(kind=8) ::    t = 2.0e-3
      real(kind=8) ::  xNA = 0.25e0 
      real(kind=8) ::   xc = 0.25e0
      real(kind=8) ::    m = 0.0e0

      Type(material_type_def) :: Mat

   end type 

!=======================================================================
! skin properties
! t = thickness, m_LE = leading edge weight/span
!=======================================================================

   type :: Skin_Info_Def

      real(kind=8) ::    t = 1.e-6
      real(kind=8) :: m_LE = 0.e0

   end type !Skin_Info_Def

!=======================================================================
! cross-section properties (skin inertias) normalized by chord and 
! skin thickness wherever applicable
! note: 1 layer uniform thickness assumed
!=======================================================================

   type :: Contour_Prop_Def

      real(kind=8) :: A_skin    = 2.0392e0
      real(kind=8) :: Izzskin   = 0.2952e0
      real(kind=8) :: Iyyskin   = 0.00398e0 
      real(kind=8) :: Atotal    = 0.0822e0 
      real(kind=8) :: mi2a_fill = 0.006993e0
      real(kind=8) :: dsbyt     = 2.0392e0 
      real(kind=8) :: mi1a_fill = 0.03456e0
      real(kind=8) :: mi1a_skin = 1.0052e0

   end type

!=======================================================================
! section properties
! m_total: mass/sapn, J_total: polar M.I of CS, xcg: chordwise cg loc 
! m_LE: leading edge mass/span, spar_SF: spar safety factor 
! skin_SF: skin safety factor 
!=======================================================================

   type :: Cross_Section_Def

    real(kind=8)    :: m, Io, xcg, m_LE, spar_SF, skin_SF, na
   
   end type 

end module airfoil_design_types
