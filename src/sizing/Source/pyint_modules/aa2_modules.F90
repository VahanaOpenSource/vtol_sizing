#define type(x) TYPE(x), target

!=======================================================================
! blade cross-section design
!=======================================================================

module CrossSectionData

   use airfoil_design_types
!=======================================================================
! spar materials
!=======================================================================

    type(material_type_def) :: Titanium, Aluminum, Carbon_uni, Carbon_090 

    integer                 :: material_init = 0 

!=======================================================================
! routines to perform initialization
!=======================================================================

    contains
        subroutine init_materials() 
        implicit none
            Titanium % rho       = 4500.0e0                 ! material density, kg/cu.m
            Titanium % sigma_y   =  880.0e6                 ! allowed stress in tension, Pa
            Titanium % E         =  117.0e9                 ! Young's modulus, Pa     
            Titanium % G         =   43.0e9                 ! Shear modulus, Pa 

            Aluminum % rho       = 2700.0e0                 ! 
            Aluminum % sigma_y   =  200.0e6
            Aluminum % E         =   71.0e9                 ! Young's modulus, Pa     
            Aluminum % G         =   27.3e9                 ! Shear modulus, Pa 

            Carbon_uni % rho     = 1600.0e0                 ! 
            Carbon_uni % sigma_y = 2400.0e6
            Carbon_uni % E       =  122.0e9                 ! Young's modulus, Pa     
            Carbon_uni % G       =   27.3e9                 ! Shear modulus, Pa 

            Carbon_090 % rho     = 1600.0e0                 ! 
            Carbon_090 % sigma_y =  600.0e6
            Carbon_090 % E       =   70.0e9                 ! Young's modulus, Pa     
            Carbon_090 % G       =   27.3e9                 ! Shear modulus, Pa 

            material_init        = 1                        ! trip flag
        end subroutine


! ###################################################################
! END OF FILE
! ###################################################################
end module CrossSectionData
