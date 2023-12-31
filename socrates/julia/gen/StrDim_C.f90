! Autogenerated from Fortran file src/radiance_core/def_dimen.F90
! svn revision 1226

module StrDim_C

use, intrinsic :: iso_c_binding

use UTILITIES_CF
use def_dimen

implicit none

private

contains

    function create_StrDim() result(dimen_Cptr) & 
        bind(C, NAME='PS_create_StrDim')

        implicit none
        type(c_ptr) :: dimen_Cptr

        TYPE(StrDim), pointer :: dimen

        ALLOCATE(dimen)

        dimen_Cptr = c_loc(dimen)
    end function create_StrDim

    subroutine delete_StrDim(dimen_Cptr) &
        bind(C, NAME='PS_delete_StrDim')

        implicit none
        type(c_ptr), value, intent(in) :: dimen_Cptr

        TYPE(StrDim), pointer :: dimen

        call C_F_POINTER(dimen_Cptr, dimen)

        DEALLOCATE(dimen)

    end subroutine delete_StrDim

    function StrDim_get_integer(dimen_Cptr, field_Cstr, val) result(field_OK) &
        bind(C, NAME='PS_StrDim_get_integer')

        implicit none
        type(c_ptr), value, intent(in)          :: dimen_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        integer(c_int), intent(out) :: val

        ! local variables
        TYPE(StrDim), pointer                  :: dimen
        character(len=:), allocatable           :: field

        ! convert types
        call C_F_POINTER(dimen_Cptr, dimen)
        field = c_to_f_string(field_Cstr)

        field_OK = .TRUE.
    
        if (field == 'nd_profile') then
            val = dimen%nd_profile
        else if (field == 'nd_layer') then
            val = dimen%nd_layer
        else if (field == 'nd_layer_clr') then
            val = dimen%nd_layer_clr
        else if (field == 'id_cloud_top') then
            val = dimen%id_cloud_top
        else if (field == 'nd_2sg_profile') then
            val = dimen%nd_2sg_profile
        else if (field == 'nd_flux_profile') then
            val = dimen%nd_flux_profile
        else if (field == 'nd_radiance_profile') then
            val = dimen%nd_radiance_profile
        else if (field == 'nd_j_profile') then
            val = dimen%nd_j_profile
        else if (field == 'nd_channel') then
            val = dimen%nd_channel
        else if (field == 'nd_max_order') then
            val = dimen%nd_max_order
        else if (field == 'nd_direction') then
            val = dimen%nd_direction
        else if (field == 'nd_viewing_level') then
            val = dimen%nd_viewing_level
        else if (field == 'nd_sph_coeff') then
            val = dimen%nd_sph_coeff
        else if (field == 'nd_brdf_basis_fnc') then
            val = dimen%nd_brdf_basis_fnc
        else if (field == 'nd_brdf_trunc') then
            val = dimen%nd_brdf_trunc
        else if (field == 'nd_tile_type') then
            val = dimen%nd_tile_type
        else if (field == 'nd_aerosol_mode') then
            val = dimen%nd_aerosol_mode
        else if (field == 'nd_profile_aerosol_prsc') then
            val = dimen%nd_profile_aerosol_prsc
        else if (field == 'nd_profile_cloud_prsc') then
            val = dimen%nd_profile_cloud_prsc
        else if (field == 'nd_opt_level_aerosol_prsc') then
            val = dimen%nd_opt_level_aerosol_prsc
        else if (field == 'nd_opt_level_cloud_prsc') then
            val = dimen%nd_opt_level_cloud_prsc
        else if (field == 'nd_phf_term_aerosol_prsc') then
            val = dimen%nd_phf_term_aerosol_prsc
        else if (field == 'nd_phf_term_cloud_prsc') then
            val = dimen%nd_phf_term_cloud_prsc
        else if (field == 'nd_cloud_component') then
            val = dimen%nd_cloud_component
        else if (field == 'nd_cloud_type') then
            val = dimen%nd_cloud_type
        else if (field == 'nd_cloud_representation') then
            val = dimen%nd_cloud_representation
        else if (field == 'nd_column') then
            val = dimen%nd_column
        else if (field == 'nd_subcol_gen') then
            val = dimen%nd_subcol_gen
        else if (field == 'nd_subcol_req') then
            val = dimen%nd_subcol_req
        else if (field == 'nd_overlap_coeff') then
            val = dimen%nd_overlap_coeff
        else if (field == 'nd_source_coeff') then
            val = dimen%nd_source_coeff
        else if (field == 'nd_region') then
            val = dimen%nd_region
        else if (field == 'nd_point_tile') then
            val = dimen%nd_point_tile
        else if (field == 'nd_tile') then
            val = dimen%nd_tile
        else
            field_OK = .FALSE.
        end if

    end function StrDim_get_integer

    function StrDim_set_integer(dimen_Cptr, field_Cstr, val) result(field_OK) &
        bind(C, NAME='PS_StrDim_set_integer')

        implicit none
        type(c_ptr), value, intent(in)          :: dimen_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        integer(c_int), value, intent(in) :: val

        ! local variables
        TYPE(StrDim), pointer                  :: dimen
        character(len=:), allocatable           :: field

        ! convert types
        call C_F_POINTER(dimen_Cptr, dimen)
        field = c_to_f_string(field_Cstr)

        field_OK = .TRUE.
    
        if (field == 'nd_profile') then
            dimen%nd_profile = val
        else if (field == 'nd_layer') then
            dimen%nd_layer = val
        else if (field == 'nd_layer_clr') then
            dimen%nd_layer_clr = val
        else if (field == 'id_cloud_top') then
            dimen%id_cloud_top = val
        else if (field == 'nd_2sg_profile') then
            dimen%nd_2sg_profile = val
        else if (field == 'nd_flux_profile') then
            dimen%nd_flux_profile = val
        else if (field == 'nd_radiance_profile') then
            dimen%nd_radiance_profile = val
        else if (field == 'nd_j_profile') then
            dimen%nd_j_profile = val
        else if (field == 'nd_channel') then
            dimen%nd_channel = val
        else if (field == 'nd_max_order') then
            dimen%nd_max_order = val
        else if (field == 'nd_direction') then
            dimen%nd_direction = val
        else if (field == 'nd_viewing_level') then
            dimen%nd_viewing_level = val
        else if (field == 'nd_sph_coeff') then
            dimen%nd_sph_coeff = val
        else if (field == 'nd_brdf_basis_fnc') then
            dimen%nd_brdf_basis_fnc = val
        else if (field == 'nd_brdf_trunc') then
            dimen%nd_brdf_trunc = val
        else if (field == 'nd_tile_type') then
            dimen%nd_tile_type = val
        else if (field == 'nd_aerosol_mode') then
            dimen%nd_aerosol_mode = val
        else if (field == 'nd_profile_aerosol_prsc') then
            dimen%nd_profile_aerosol_prsc = val
        else if (field == 'nd_profile_cloud_prsc') then
            dimen%nd_profile_cloud_prsc = val
        else if (field == 'nd_opt_level_aerosol_prsc') then
            dimen%nd_opt_level_aerosol_prsc = val
        else if (field == 'nd_opt_level_cloud_prsc') then
            dimen%nd_opt_level_cloud_prsc = val
        else if (field == 'nd_phf_term_aerosol_prsc') then
            dimen%nd_phf_term_aerosol_prsc = val
        else if (field == 'nd_phf_term_cloud_prsc') then
            dimen%nd_phf_term_cloud_prsc = val
        else if (field == 'nd_cloud_component') then
            dimen%nd_cloud_component = val
        else if (field == 'nd_cloud_type') then
            dimen%nd_cloud_type = val
        else if (field == 'nd_cloud_representation') then
            dimen%nd_cloud_representation = val
        else if (field == 'nd_column') then
            dimen%nd_column = val
        else if (field == 'nd_subcol_gen') then
            dimen%nd_subcol_gen = val
        else if (field == 'nd_subcol_req') then
            dimen%nd_subcol_req = val
        else if (field == 'nd_overlap_coeff') then
            dimen%nd_overlap_coeff = val
        else if (field == 'nd_source_coeff') then
            dimen%nd_source_coeff = val
        else if (field == 'nd_region') then
            dimen%nd_region = val
        else if (field == 'nd_point_tile') then
            dimen%nd_point_tile = val
        else if (field == 'nd_tile') then
            dimen%nd_tile = val
        else
            field_OK = .FALSE.
        end if

    end function StrDim_set_integer

end module StrDim_C
