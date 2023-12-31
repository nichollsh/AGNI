! Autogenerated from Fortran file src/radiance_core/def_aer.F90
! svn revision 1226

module StrAer_C

use, intrinsic :: iso_c_binding

use UTILITIES_CF
use def_aer

implicit none

private

contains

    function create_StrAer() result(aer_Cptr) & 
        bind(C, NAME='PS_create_StrAer')

        implicit none
        type(c_ptr) :: aer_Cptr

        TYPE(StrAer), pointer :: aer

        ALLOCATE(aer)

        aer_Cptr = c_loc(aer)
    end function create_StrAer

    subroutine delete_StrAer(aer_Cptr) &
        bind(C, NAME='PS_delete_StrAer')

        implicit none
        type(c_ptr), value, intent(in) :: aer_Cptr

        TYPE(StrAer), pointer :: aer

        call C_F_POINTER(aer_Cptr, aer)

        DEALLOCATE(aer)

    end subroutine delete_StrAer

    function StrAer_get_integer(aer_Cptr, field_Cstr, val) result(field_OK) &
        bind(C, NAME='PS_StrAer_get_integer')

        implicit none
        type(c_ptr), value, intent(in)          :: aer_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        integer(c_int), intent(out) :: val

        ! local variables
        TYPE(StrAer), pointer                  :: aer
        character(len=:), allocatable           :: field

        ! convert types
        call C_F_POINTER(aer_Cptr, aer)
        field = c_to_f_string(field_Cstr)

        field_OK = .TRUE.
    
        if (field == 'n_mode') then
            val = aer%n_mode
        else
            field_OK = .FALSE.
        end if

    end function StrAer_get_integer

    function StrAer_set_integer(aer_Cptr, field_Cstr, val) result(field_OK) &
        bind(C, NAME='PS_StrAer_set_integer')

        implicit none
        type(c_ptr), value, intent(in)          :: aer_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        integer(c_int), value, intent(in) :: val

        ! local variables
        TYPE(StrAer), pointer                  :: aer
        character(len=:), allocatable           :: field

        ! convert types
        call C_F_POINTER(aer_Cptr, aer)
        field = c_to_f_string(field_Cstr)

        field_OK = .TRUE.
    
        if (field == 'n_mode') then
            aer%n_mode = val
        else
            field_OK = .FALSE.
        end if

    end function StrAer_set_integer

    function StrAer_get_integer_array(aer_Cptr, field_Cstr, loc, dims_C, ndim, lbounds_C) result(field_OK) &
        bind(C, NAME='PS_StrAer_get_integer_array')

        implicit none
        type(c_ptr), value, intent(in)          :: aer_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        type(c_ptr), intent(out)                :: loc
        type(c_ptr), value, intent(in)          :: dims_C
        integer(c_int), intent(inout)           :: ndim
        type(c_ptr), value, intent(in)          :: lbounds_C

        ! local variables
        TYPE(StrAer), pointer             :: aer
        character(len=:), allocatable           :: field
        integer(c_int), pointer, dimension(:)   :: dims, lbounds    

        ! convert types
        call C_F_POINTER(aer_Cptr, aer)
        field = c_to_f_string(field_Cstr)
        call C_F_POINTER(dims_C, dims, [ndim])
        call C_F_POINTER(lbounds_C, lbounds, [ndim])

        field_OK = .TRUE.
    
        if (field == 'mr_type_index') then
            if (allocated(aer%mr_type_index)) then
                loc = C_LOC(aer%mr_type_index)
                ndim = size(shape(aer%mr_type_index))
                dims(1:ndim) = shape(aer%mr_type_index)
                lbounds(1:ndim) = lbound(aer%mr_type_index)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'mr_source') then
            if (allocated(aer%mr_source)) then
                loc = C_LOC(aer%mr_source)
                ndim = size(shape(aer%mr_source))
                dims(1:ndim) = shape(aer%mr_source)
                lbounds(1:ndim) = lbound(aer%mr_source)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'n_opt_level_prsc') then
            if (allocated(aer%n_opt_level_prsc)) then
                loc = C_LOC(aer%n_opt_level_prsc)
                ndim = size(shape(aer%n_opt_level_prsc))
                dims(1:ndim) = shape(aer%n_opt_level_prsc)
                lbounds(1:ndim) = lbound(aer%n_opt_level_prsc)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'n_phase_term_prsc') then
            if (allocated(aer%n_phase_term_prsc)) then
                loc = C_LOC(aer%n_phase_term_prsc)
                ndim = size(shape(aer%n_phase_term_prsc))
                dims(1:ndim) = shape(aer%n_phase_term_prsc)
                lbounds(1:ndim) = lbound(aer%n_phase_term_prsc)
            else
                loc = C_NULL_PTR
            end if

        else
            field_OK = .FALSE.
        end if

    end function StrAer_get_integer_array

    function StrAer_get_real_array(aer_Cptr, field_Cstr, loc, dims_C, ndim, lbounds_C) result(field_OK) &
        bind(C, NAME='PS_StrAer_get_real_array')

        implicit none
        type(c_ptr), value, intent(in)          :: aer_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        type(c_ptr), intent(out)                :: loc
        type(c_ptr), value, intent(in)          :: dims_C
        integer(c_int), intent(inout)           :: ndim
        type(c_ptr), value, intent(in)          :: lbounds_C

        ! local variables
        TYPE(StrAer), pointer             :: aer
        character(len=:), allocatable           :: field
        integer(c_int), pointer, dimension(:)   :: dims, lbounds    

        ! convert types
        call C_F_POINTER(aer_Cptr, aer)
        field = c_to_f_string(field_Cstr)
        call C_F_POINTER(dims_C, dims, [ndim])
        call C_F_POINTER(lbounds_C, lbounds, [ndim])

        field_OK = .TRUE.
    
        if (field == 'mix_ratio') then
            if (allocated(aer%mix_ratio)) then
                loc = C_LOC(aer%mix_ratio)
                ndim = size(shape(aer%mix_ratio))
                dims(1:ndim) = shape(aer%mix_ratio)
                lbounds(1:ndim) = lbound(aer%mix_ratio)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'mean_rel_humidity') then
            if (allocated(aer%mean_rel_humidity)) then
                loc = C_LOC(aer%mean_rel_humidity)
                ndim = size(shape(aer%mean_rel_humidity))
                dims(1:ndim) = shape(aer%mean_rel_humidity)
                lbounds(1:ndim) = lbound(aer%mean_rel_humidity)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'mode_mix_ratio') then
            if (allocated(aer%mode_mix_ratio)) then
                loc = C_LOC(aer%mode_mix_ratio)
                ndim = size(shape(aer%mode_mix_ratio))
                dims(1:ndim) = shape(aer%mode_mix_ratio)
                lbounds(1:ndim) = lbound(aer%mode_mix_ratio)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'mode_absorption') then
            if (allocated(aer%mode_absorption)) then
                loc = C_LOC(aer%mode_absorption)
                ndim = size(shape(aer%mode_absorption))
                dims(1:ndim) = shape(aer%mode_absorption)
                lbounds(1:ndim) = lbound(aer%mode_absorption)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'mode_scattering') then
            if (allocated(aer%mode_scattering)) then
                loc = C_LOC(aer%mode_scattering)
                ndim = size(shape(aer%mode_scattering))
                dims(1:ndim) = shape(aer%mode_scattering)
                lbounds(1:ndim) = lbound(aer%mode_scattering)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'mode_asymmetry') then
            if (allocated(aer%mode_asymmetry)) then
                loc = C_LOC(aer%mode_asymmetry)
                ndim = size(shape(aer%mode_asymmetry))
                dims(1:ndim) = shape(aer%mode_asymmetry)
                lbounds(1:ndim) = lbound(aer%mode_asymmetry)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'pressure_prsc') then
            if (allocated(aer%pressure_prsc)) then
                loc = C_LOC(aer%pressure_prsc)
                ndim = size(shape(aer%pressure_prsc))
                dims(1:ndim) = shape(aer%pressure_prsc)
                lbounds(1:ndim) = lbound(aer%pressure_prsc)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'absorption_prsc') then
            if (allocated(aer%absorption_prsc)) then
                loc = C_LOC(aer%absorption_prsc)
                ndim = size(shape(aer%absorption_prsc))
                dims(1:ndim) = shape(aer%absorption_prsc)
                lbounds(1:ndim) = lbound(aer%absorption_prsc)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'scattering_prsc') then
            if (allocated(aer%scattering_prsc)) then
                loc = C_LOC(aer%scattering_prsc)
                ndim = size(shape(aer%scattering_prsc))
                dims(1:ndim) = shape(aer%scattering_prsc)
                lbounds(1:ndim) = lbound(aer%scattering_prsc)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'phase_fnc_prsc') then
            if (allocated(aer%phase_fnc_prsc)) then
                loc = C_LOC(aer%phase_fnc_prsc)
                ndim = size(shape(aer%phase_fnc_prsc))
                dims(1:ndim) = shape(aer%phase_fnc_prsc)
                lbounds(1:ndim) = lbound(aer%phase_fnc_prsc)
            else
                loc = C_NULL_PTR
            end if

        else
            field_OK = .FALSE.
        end if

    end function StrAer_get_real_array

end module StrAer_C
