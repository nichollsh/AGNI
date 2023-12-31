! Autogenerated from Fortran file src/radiance_core/def_out.F90
! svn revision 1226

module StrOut_C

use, intrinsic :: iso_c_binding

use UTILITIES_CF
use def_out

implicit none

private

contains

    function create_StrOut() result(radout_Cptr) & 
        bind(C, NAME='PS_create_StrOut')

        implicit none
        type(c_ptr) :: radout_Cptr

        TYPE(StrOut), pointer :: radout

        ALLOCATE(radout)

        radout_Cptr = c_loc(radout)
    end function create_StrOut

    subroutine delete_StrOut(radout_Cptr) &
        bind(C, NAME='PS_delete_StrOut')

        implicit none
        type(c_ptr), value, intent(in) :: radout_Cptr

        TYPE(StrOut), pointer :: radout

        call C_F_POINTER(radout_Cptr, radout)

        DEALLOCATE(radout)

    end subroutine delete_StrOut

    function StrOut_get_real_array(radout_Cptr, field_Cstr, loc, dims_C, ndim, lbounds_C) result(field_OK) &
        bind(C, NAME='PS_StrOut_get_real_array')

        implicit none
        type(c_ptr), value, intent(in)          :: radout_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        type(c_ptr), intent(out)                :: loc
        type(c_ptr), value, intent(in)          :: dims_C
        integer(c_int), intent(inout)           :: ndim
        type(c_ptr), value, intent(in)          :: lbounds_C

        ! local variables
        TYPE(StrOut), pointer             :: radout
        character(len=:), allocatable           :: field
        integer(c_int), pointer, dimension(:)   :: dims, lbounds    

        ! convert types
        call C_F_POINTER(radout_Cptr, radout)
        field = c_to_f_string(field_Cstr)
        call C_F_POINTER(dims_C, dims, [ndim])
        call C_F_POINTER(lbounds_C, lbounds, [ndim])

        field_OK = .TRUE.
    
        if (field == 'flux_direct') then
            if (allocated(radout%flux_direct)) then
                loc = C_LOC(radout%flux_direct)
                ndim = size(shape(radout%flux_direct))
                dims(1:ndim) = shape(radout%flux_direct)
                lbounds(1:ndim) = lbound(radout%flux_direct)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_down') then
            if (allocated(radout%flux_down)) then
                loc = C_LOC(radout%flux_down)
                ndim = size(shape(radout%flux_down))
                dims(1:ndim) = shape(radout%flux_down)
                lbounds(1:ndim) = lbound(radout%flux_down)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_up') then
            if (allocated(radout%flux_up)) then
                loc = C_LOC(radout%flux_up)
                ndim = size(shape(radout%flux_up))
                dims(1:ndim) = shape(radout%flux_up)
                lbounds(1:ndim) = lbound(radout%flux_up)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_div') then
            if (allocated(radout%flux_div)) then
                loc = C_LOC(radout%flux_div)
                ndim = size(shape(radout%flux_div))
                dims(1:ndim) = shape(radout%flux_div)
                lbounds(1:ndim) = lbound(radout%flux_div)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_direct_clear') then
            if (allocated(radout%flux_direct_clear)) then
                loc = C_LOC(radout%flux_direct_clear)
                ndim = size(shape(radout%flux_direct_clear))
                dims(1:ndim) = shape(radout%flux_direct_clear)
                lbounds(1:ndim) = lbound(radout%flux_direct_clear)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_down_clear') then
            if (allocated(radout%flux_down_clear)) then
                loc = C_LOC(radout%flux_down_clear)
                ndim = size(shape(radout%flux_down_clear))
                dims(1:ndim) = shape(radout%flux_down_clear)
                lbounds(1:ndim) = lbound(radout%flux_down_clear)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_up_clear') then
            if (allocated(radout%flux_up_clear)) then
                loc = C_LOC(radout%flux_up_clear)
                ndim = size(shape(radout%flux_up_clear))
                dims(1:ndim) = shape(radout%flux_up_clear)
                lbounds(1:ndim) = lbound(radout%flux_up_clear)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_div_clear') then
            if (allocated(radout%flux_div_clear)) then
                loc = C_LOC(radout%flux_div_clear)
                ndim = size(shape(radout%flux_div_clear))
                dims(1:ndim) = shape(radout%flux_div_clear)
                lbounds(1:ndim) = lbound(radout%flux_div_clear)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_direct_div') then
            if (allocated(radout%flux_direct_div)) then
                loc = C_LOC(radout%flux_direct_div)
                ndim = size(shape(radout%flux_direct_div))
                dims(1:ndim) = shape(radout%flux_direct_div)
                lbounds(1:ndim) = lbound(radout%flux_direct_div)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_direct_sph') then
            if (allocated(radout%flux_direct_sph)) then
                loc = C_LOC(radout%flux_direct_sph)
                ndim = size(shape(radout%flux_direct_sph))
                dims(1:ndim) = shape(radout%flux_direct_sph)
                lbounds(1:ndim) = lbound(radout%flux_direct_sph)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_direct_clear_div') then
            if (allocated(radout%flux_direct_clear_div)) then
                loc = C_LOC(radout%flux_direct_clear_div)
                ndim = size(shape(radout%flux_direct_clear_div))
                dims(1:ndim) = shape(radout%flux_direct_clear_div)
                lbounds(1:ndim) = lbound(radout%flux_direct_clear_div)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_direct_clear_sph') then
            if (allocated(radout%flux_direct_clear_sph)) then
                loc = C_LOC(radout%flux_direct_clear_sph)
                ndim = size(shape(radout%flux_direct_clear_sph))
                dims(1:ndim) = shape(radout%flux_direct_clear_sph)
                lbounds(1:ndim) = lbound(radout%flux_direct_clear_sph)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'radiance') then
            if (allocated(radout%radiance)) then
                loc = C_LOC(radout%radiance)
                ndim = size(shape(radout%radiance))
                dims(1:ndim) = shape(radout%radiance)
                lbounds(1:ndim) = lbound(radout%radiance)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'photolysis') then
            if (allocated(radout%photolysis)) then
                loc = C_LOC(radout%photolysis)
                ndim = size(shape(radout%photolysis))
                dims(1:ndim) = shape(radout%photolysis)
                lbounds(1:ndim) = lbound(radout%photolysis)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'solar_tail_flux') then
            if (allocated(radout%solar_tail_flux)) then
                loc = C_LOC(radout%solar_tail_flux)
                ndim = size(shape(radout%solar_tail_flux))
                dims(1:ndim) = shape(radout%solar_tail_flux)
                lbounds(1:ndim) = lbound(radout%solar_tail_flux)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'contrib_funci') then
            if (allocated(radout%contrib_funci)) then
                loc = C_LOC(radout%contrib_funci)
                ndim = size(shape(radout%contrib_funci))
                dims(1:ndim) = shape(radout%contrib_funci)
                lbounds(1:ndim) = lbound(radout%contrib_funci)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'contrib_funcf') then
            if (allocated(radout%contrib_funcf)) then
                loc = C_LOC(radout%contrib_funcf)
                ndim = size(shape(radout%contrib_funcf))
                dims(1:ndim) = shape(radout%contrib_funcf)
                lbounds(1:ndim) = lbound(radout%contrib_funcf)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_up_tile') then
            if (allocated(radout%flux_up_tile)) then
                loc = C_LOC(radout%flux_up_tile)
                ndim = size(shape(radout%flux_up_tile))
                dims(1:ndim) = shape(radout%flux_up_tile)
                lbounds(1:ndim) = lbound(radout%flux_up_tile)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_up_blue_tile') then
            if (allocated(radout%flux_up_blue_tile)) then
                loc = C_LOC(radout%flux_up_blue_tile)
                ndim = size(shape(radout%flux_up_blue_tile))
                dims(1:ndim) = shape(radout%flux_up_blue_tile)
                lbounds(1:ndim) = lbound(radout%flux_up_blue_tile)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_direct_blue_surf') then
            if (allocated(radout%flux_direct_blue_surf)) then
                loc = C_LOC(radout%flux_direct_blue_surf)
                ndim = size(shape(radout%flux_direct_blue_surf))
                dims(1:ndim) = shape(radout%flux_direct_blue_surf)
                lbounds(1:ndim) = lbound(radout%flux_direct_blue_surf)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_down_blue_surf') then
            if (allocated(radout%flux_down_blue_surf)) then
                loc = C_LOC(radout%flux_down_blue_surf)
                ndim = size(shape(radout%flux_down_blue_surf))
                dims(1:ndim) = shape(radout%flux_down_blue_surf)
                lbounds(1:ndim) = lbound(radout%flux_down_blue_surf)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_up_blue_surf') then
            if (allocated(radout%flux_up_blue_surf)) then
                loc = C_LOC(radout%flux_up_blue_surf)
                ndim = size(shape(radout%flux_up_blue_surf))
                dims(1:ndim) = shape(radout%flux_up_blue_surf)
                lbounds(1:ndim) = lbound(radout%flux_up_blue_surf)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_direct_band') then
            if (allocated(radout%flux_direct_band)) then
                loc = C_LOC(radout%flux_direct_band)
                ndim = size(shape(radout%flux_direct_band))
                dims(1:ndim) = shape(radout%flux_direct_band)
                lbounds(1:ndim) = lbound(radout%flux_direct_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_direct_div_band') then
            if (allocated(radout%flux_direct_div_band)) then
                loc = C_LOC(radout%flux_direct_div_band)
                ndim = size(shape(radout%flux_direct_div_band))
                dims(1:ndim) = shape(radout%flux_direct_div_band)
                lbounds(1:ndim) = lbound(radout%flux_direct_div_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_direct_sph_band') then
            if (allocated(radout%flux_direct_sph_band)) then
                loc = C_LOC(radout%flux_direct_sph_band)
                ndim = size(shape(radout%flux_direct_sph_band))
                dims(1:ndim) = shape(radout%flux_direct_sph_band)
                lbounds(1:ndim) = lbound(radout%flux_direct_sph_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_down_band') then
            if (allocated(radout%flux_down_band)) then
                loc = C_LOC(radout%flux_down_band)
                ndim = size(shape(radout%flux_down_band))
                dims(1:ndim) = shape(radout%flux_down_band)
                lbounds(1:ndim) = lbound(radout%flux_down_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_up_band') then
            if (allocated(radout%flux_up_band)) then
                loc = C_LOC(radout%flux_up_band)
                ndim = size(shape(radout%flux_up_band))
                dims(1:ndim) = shape(radout%flux_up_band)
                lbounds(1:ndim) = lbound(radout%flux_up_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_div_band') then
            if (allocated(radout%flux_div_band)) then
                loc = C_LOC(radout%flux_div_band)
                ndim = size(shape(radout%flux_div_band))
                dims(1:ndim) = shape(radout%flux_div_band)
                lbounds(1:ndim) = lbound(radout%flux_div_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_direct_clear_band') then
            if (allocated(radout%flux_direct_clear_band)) then
                loc = C_LOC(radout%flux_direct_clear_band)
                ndim = size(shape(radout%flux_direct_clear_band))
                dims(1:ndim) = shape(radout%flux_direct_clear_band)
                lbounds(1:ndim) = lbound(radout%flux_direct_clear_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_direct_clear_div_band') then
            if (allocated(radout%flux_direct_clear_div_band)) then
                loc = C_LOC(radout%flux_direct_clear_div_band)
                ndim = size(shape(radout%flux_direct_clear_div_band))
                dims(1:ndim) = shape(radout%flux_direct_clear_div_band)
                lbounds(1:ndim) = lbound(radout%flux_direct_clear_div_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_direct_clear_sph_band') then
            if (allocated(radout%flux_direct_clear_sph_band)) then
                loc = C_LOC(radout%flux_direct_clear_sph_band)
                ndim = size(shape(radout%flux_direct_clear_sph_band))
                dims(1:ndim) = shape(radout%flux_direct_clear_sph_band)
                lbounds(1:ndim) = lbound(radout%flux_direct_clear_sph_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_down_clear_band') then
            if (allocated(radout%flux_down_clear_band)) then
                loc = C_LOC(radout%flux_down_clear_band)
                ndim = size(shape(radout%flux_down_clear_band))
                dims(1:ndim) = shape(radout%flux_down_clear_band)
                lbounds(1:ndim) = lbound(radout%flux_down_clear_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_up_clear_band') then
            if (allocated(radout%flux_up_clear_band)) then
                loc = C_LOC(radout%flux_up_clear_band)
                ndim = size(shape(radout%flux_up_clear_band))
                dims(1:ndim) = shape(radout%flux_up_clear_band)
                lbounds(1:ndim) = lbound(radout%flux_up_clear_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'flux_div_clear_band') then
            if (allocated(radout%flux_div_clear_band)) then
                loc = C_LOC(radout%flux_div_clear_band)
                ndim = size(shape(radout%flux_div_clear_band))
                dims(1:ndim) = shape(radout%flux_div_clear_band)
                lbounds(1:ndim) = lbound(radout%flux_div_clear_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'contrib_funci_band') then
            if (allocated(radout%contrib_funci_band)) then
                loc = C_LOC(radout%contrib_funci_band)
                ndim = size(shape(radout%contrib_funci_band))
                dims(1:ndim) = shape(radout%contrib_funci_band)
                lbounds(1:ndim) = lbound(radout%contrib_funci_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'contrib_funcf_band') then
            if (allocated(radout%contrib_funcf_band)) then
                loc = C_LOC(radout%contrib_funcf_band)
                ndim = size(shape(radout%contrib_funcf_band))
                dims(1:ndim) = shape(radout%contrib_funcf_band)
                lbounds(1:ndim) = lbound(radout%contrib_funcf_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'actinic_flux') then
            if (allocated(radout%actinic_flux)) then
                loc = C_LOC(radout%actinic_flux)
                ndim = size(shape(radout%actinic_flux))
                dims(1:ndim) = shape(radout%actinic_flux)
                lbounds(1:ndim) = lbound(radout%actinic_flux)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'actinic_flux_clear') then
            if (allocated(radout%actinic_flux_clear)) then
                loc = C_LOC(radout%actinic_flux_clear)
                ndim = size(shape(radout%actinic_flux_clear))
                dims(1:ndim) = shape(radout%actinic_flux_clear)
                lbounds(1:ndim) = lbound(radout%actinic_flux_clear)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'actinic_flux_band') then
            if (allocated(radout%actinic_flux_band)) then
                loc = C_LOC(radout%actinic_flux_band)
                ndim = size(shape(radout%actinic_flux_band))
                dims(1:ndim) = shape(radout%actinic_flux_band)
                lbounds(1:ndim) = lbound(radout%actinic_flux_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'actinic_flux_clear_band') then
            if (allocated(radout%actinic_flux_clear_band)) then
                loc = C_LOC(radout%actinic_flux_clear_band)
                ndim = size(shape(radout%actinic_flux_clear_band))
                dims(1:ndim) = shape(radout%actinic_flux_clear_band)
                lbounds(1:ndim) = lbound(radout%actinic_flux_clear_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'photolysis_rate') then
            if (allocated(radout%photolysis_rate)) then
                loc = C_LOC(radout%photolysis_rate)
                ndim = size(shape(radout%photolysis_rate))
                dims(1:ndim) = shape(radout%photolysis_rate)
                lbounds(1:ndim) = lbound(radout%photolysis_rate)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'photolysis_rate_clear') then
            if (allocated(radout%photolysis_rate_clear)) then
                loc = C_LOC(radout%photolysis_rate_clear)
                ndim = size(shape(radout%photolysis_rate_clear))
                dims(1:ndim) = shape(radout%photolysis_rate_clear)
                lbounds(1:ndim) = lbound(radout%photolysis_rate_clear)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'photolysis_div') then
            if (allocated(radout%photolysis_div)) then
                loc = C_LOC(radout%photolysis_div)
                ndim = size(shape(radout%photolysis_div))
                dims(1:ndim) = shape(radout%photolysis_div)
                lbounds(1:ndim) = lbound(radout%photolysis_div)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'photolysis_div_clear') then
            if (allocated(radout%photolysis_div_clear)) then
                loc = C_LOC(radout%photolysis_div_clear)
                ndim = size(shape(radout%photolysis_div_clear))
                dims(1:ndim) = shape(radout%photolysis_div_clear)
                lbounds(1:ndim) = lbound(radout%photolysis_div_clear)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'tot_cloud_cover') then
            if (allocated(radout%tot_cloud_cover)) then
                loc = C_LOC(radout%tot_cloud_cover)
                ndim = size(shape(radout%tot_cloud_cover))
                dims(1:ndim) = shape(radout%tot_cloud_cover)
                lbounds(1:ndim) = lbound(radout%tot_cloud_cover)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'cloud_absorptivity') then
            if (allocated(radout%cloud_absorptivity)) then
                loc = C_LOC(radout%cloud_absorptivity)
                ndim = size(shape(radout%cloud_absorptivity))
                dims(1:ndim) = shape(radout%cloud_absorptivity)
                lbounds(1:ndim) = lbound(radout%cloud_absorptivity)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'cloud_weight_absorptivity') then
            if (allocated(radout%cloud_weight_absorptivity)) then
                loc = C_LOC(radout%cloud_weight_absorptivity)
                ndim = size(shape(radout%cloud_weight_absorptivity))
                dims(1:ndim) = shape(radout%cloud_weight_absorptivity)
                lbounds(1:ndim) = lbound(radout%cloud_weight_absorptivity)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'ls_cloud_absorptivity') then
            if (allocated(radout%ls_cloud_absorptivity)) then
                loc = C_LOC(radout%ls_cloud_absorptivity)
                ndim = size(shape(radout%ls_cloud_absorptivity))
                dims(1:ndim) = shape(radout%ls_cloud_absorptivity)
                lbounds(1:ndim) = lbound(radout%ls_cloud_absorptivity)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'ls_cloud_weight_absorptivity') then
            if (allocated(radout%ls_cloud_weight_absorptivity)) then
                loc = C_LOC(radout%ls_cloud_weight_absorptivity)
                ndim = size(shape(radout%ls_cloud_weight_absorptivity))
                dims(1:ndim) = shape(radout%ls_cloud_weight_absorptivity)
                lbounds(1:ndim) = lbound(radout%ls_cloud_weight_absorptivity)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'cnv_cloud_absorptivity') then
            if (allocated(radout%cnv_cloud_absorptivity)) then
                loc = C_LOC(radout%cnv_cloud_absorptivity)
                ndim = size(shape(radout%cnv_cloud_absorptivity))
                dims(1:ndim) = shape(radout%cnv_cloud_absorptivity)
                lbounds(1:ndim) = lbound(radout%cnv_cloud_absorptivity)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'cnv_cloud_weight_absorptivity') then
            if (allocated(radout%cnv_cloud_weight_absorptivity)) then
                loc = C_LOC(radout%cnv_cloud_weight_absorptivity)
                ndim = size(shape(radout%cnv_cloud_weight_absorptivity))
                dims(1:ndim) = shape(radout%cnv_cloud_weight_absorptivity)
                lbounds(1:ndim) = lbound(radout%cnv_cloud_weight_absorptivity)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'cloud_extinction') then
            if (allocated(radout%cloud_extinction)) then
                loc = C_LOC(radout%cloud_extinction)
                ndim = size(shape(radout%cloud_extinction))
                dims(1:ndim) = shape(radout%cloud_extinction)
                lbounds(1:ndim) = lbound(radout%cloud_extinction)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'cloud_weight_extinction') then
            if (allocated(radout%cloud_weight_extinction)) then
                loc = C_LOC(radout%cloud_weight_extinction)
                ndim = size(shape(radout%cloud_weight_extinction))
                dims(1:ndim) = shape(radout%cloud_weight_extinction)
                lbounds(1:ndim) = lbound(radout%cloud_weight_extinction)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'ls_cloud_extinction') then
            if (allocated(radout%ls_cloud_extinction)) then
                loc = C_LOC(radout%ls_cloud_extinction)
                ndim = size(shape(radout%ls_cloud_extinction))
                dims(1:ndim) = shape(radout%ls_cloud_extinction)
                lbounds(1:ndim) = lbound(radout%ls_cloud_extinction)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'ls_cloud_weight_extinction') then
            if (allocated(radout%ls_cloud_weight_extinction)) then
                loc = C_LOC(radout%ls_cloud_weight_extinction)
                ndim = size(shape(radout%ls_cloud_weight_extinction))
                dims(1:ndim) = shape(radout%ls_cloud_weight_extinction)
                lbounds(1:ndim) = lbound(radout%ls_cloud_weight_extinction)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'cnv_cloud_extinction') then
            if (allocated(radout%cnv_cloud_extinction)) then
                loc = C_LOC(radout%cnv_cloud_extinction)
                ndim = size(shape(radout%cnv_cloud_extinction))
                dims(1:ndim) = shape(radout%cnv_cloud_extinction)
                lbounds(1:ndim) = lbound(radout%cnv_cloud_extinction)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'cnv_cloud_weight_extinction') then
            if (allocated(radout%cnv_cloud_weight_extinction)) then
                loc = C_LOC(radout%cnv_cloud_weight_extinction)
                ndim = size(shape(radout%cnv_cloud_weight_extinction))
                dims(1:ndim) = shape(radout%cnv_cloud_weight_extinction)
                lbounds(1:ndim) = lbound(radout%cnv_cloud_weight_extinction)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'aerosol_absorption_band') then
            if (allocated(radout%aerosol_absorption_band)) then
                loc = C_LOC(radout%aerosol_absorption_band)
                ndim = size(shape(radout%aerosol_absorption_band))
                dims(1:ndim) = shape(radout%aerosol_absorption_band)
                lbounds(1:ndim) = lbound(radout%aerosol_absorption_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'aerosol_scattering_band') then
            if (allocated(radout%aerosol_scattering_band)) then
                loc = C_LOC(radout%aerosol_scattering_band)
                ndim = size(shape(radout%aerosol_scattering_band))
                dims(1:ndim) = shape(radout%aerosol_scattering_band)
                lbounds(1:ndim) = lbound(radout%aerosol_scattering_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'aerosol_asymmetry_band') then
            if (allocated(radout%aerosol_asymmetry_band)) then
                loc = C_LOC(radout%aerosol_asymmetry_band)
                ndim = size(shape(radout%aerosol_asymmetry_band))
                dims(1:ndim) = shape(radout%aerosol_asymmetry_band)
                lbounds(1:ndim) = lbound(radout%aerosol_asymmetry_band)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'spherical_path') then
            if (allocated(radout%spherical_path)) then
                loc = C_LOC(radout%spherical_path)
                ndim = size(shape(radout%spherical_path))
                dims(1:ndim) = shape(radout%spherical_path)
                lbounds(1:ndim) = lbound(radout%spherical_path)
            else
                loc = C_NULL_PTR
            end if

        else
            field_OK = .FALSE.
        end if

    end function StrOut_get_real_array

end module StrOut_C
