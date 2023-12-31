! Autogenerated from Fortran file src/radiance_core/def_atm.F90
! svn revision 1226

module StrAtm_C

use, intrinsic :: iso_c_binding

use UTILITIES_CF
use def_atm

implicit none

private

contains

    function create_StrAtm() result(atm_Cptr) & 
        bind(C, NAME='PS_create_StrAtm')

        implicit none
        type(c_ptr) :: atm_Cptr

        TYPE(StrAtm), pointer :: atm

        ALLOCATE(atm)

        atm_Cptr = c_loc(atm)
    end function create_StrAtm

    subroutine delete_StrAtm(atm_Cptr) &
        bind(C, NAME='PS_delete_StrAtm')

        implicit none
        type(c_ptr), value, intent(in) :: atm_Cptr

        TYPE(StrAtm), pointer :: atm

        call C_F_POINTER(atm_Cptr, atm)

        DEALLOCATE(atm)

    end subroutine delete_StrAtm

    function StrAtm_get_integer(atm_Cptr, field_Cstr, val) result(field_OK) &
        bind(C, NAME='PS_StrAtm_get_integer')

        implicit none
        type(c_ptr), value, intent(in)          :: atm_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        integer(c_int), intent(out) :: val

        ! local variables
        TYPE(StrAtm), pointer                  :: atm
        character(len=:), allocatable           :: field

        ! convert types
        call C_F_POINTER(atm_Cptr, atm)
        field = c_to_f_string(field_Cstr)

        field_OK = .TRUE.
    
        if (field == 'n_profile') then
            val = atm%n_profile
        else if (field == 'n_layer') then
            val = atm%n_layer
        else if (field == 'n_direction') then
            val = atm%n_direction
        else if (field == 'n_viewing_level') then
            val = atm%n_viewing_level
        else
            field_OK = .FALSE.
        end if

    end function StrAtm_get_integer

    function StrAtm_set_integer(atm_Cptr, field_Cstr, val) result(field_OK) &
        bind(C, NAME='PS_StrAtm_set_integer')

        implicit none
        type(c_ptr), value, intent(in)          :: atm_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        integer(c_int), value, intent(in) :: val

        ! local variables
        TYPE(StrAtm), pointer                  :: atm
        character(len=:), allocatable           :: field

        ! convert types
        call C_F_POINTER(atm_Cptr, atm)
        field = c_to_f_string(field_Cstr)

        field_OK = .TRUE.
    
        if (field == 'n_profile') then
            atm%n_profile = val
        else if (field == 'n_layer') then
            atm%n_layer = val
        else if (field == 'n_direction') then
            atm%n_direction = val
        else if (field == 'n_viewing_level') then
            atm%n_viewing_level = val
        else
            field_OK = .FALSE.
        end if

    end function StrAtm_set_integer

    function StrAtm_get_real_array(atm_Cptr, field_Cstr, loc, dims_C, ndim, lbounds_C) result(field_OK) &
        bind(C, NAME='PS_StrAtm_get_real_array')

        implicit none
        type(c_ptr), value, intent(in)          :: atm_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        type(c_ptr), intent(out)                :: loc
        type(c_ptr), value, intent(in)          :: dims_C
        integer(c_int), intent(inout)           :: ndim
        type(c_ptr), value, intent(in)          :: lbounds_C

        ! local variables
        TYPE(StrAtm), pointer             :: atm
        character(len=:), allocatable           :: field
        integer(c_int), pointer, dimension(:)   :: dims, lbounds    

        ! convert types
        call C_F_POINTER(atm_Cptr, atm)
        field = c_to_f_string(field_Cstr)
        call C_F_POINTER(dims_C, dims, [ndim])
        call C_F_POINTER(lbounds_C, lbounds, [ndim])

        field_OK = .TRUE.
    
        if (field == 'lon') then
            if (allocated(atm%lon)) then
                loc = C_LOC(atm%lon)
                ndim = size(shape(atm%lon))
                dims(1:ndim) = shape(atm%lon)
                lbounds(1:ndim) = lbound(atm%lon)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'lat') then
            if (allocated(atm%lat)) then
                loc = C_LOC(atm%lat)
                ndim = size(shape(atm%lat))
                dims(1:ndim) = shape(atm%lat)
                lbounds(1:ndim) = lbound(atm%lat)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'direction') then
            if (allocated(atm%direction)) then
                loc = C_LOC(atm%direction)
                ndim = size(shape(atm%direction))
                dims(1:ndim) = shape(atm%direction)
                lbounds(1:ndim) = lbound(atm%direction)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'viewing_level') then
            if (allocated(atm%viewing_level)) then
                loc = C_LOC(atm%viewing_level)
                ndim = size(shape(atm%viewing_level))
                dims(1:ndim) = shape(atm%viewing_level)
                lbounds(1:ndim) = lbound(atm%viewing_level)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'mass') then
            if (allocated(atm%mass)) then
                loc = C_LOC(atm%mass)
                ndim = size(shape(atm%mass))
                dims(1:ndim) = shape(atm%mass)
                lbounds(1:ndim) = lbound(atm%mass)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'density') then
            if (allocated(atm%density)) then
                loc = C_LOC(atm%density)
                ndim = size(shape(atm%density))
                dims(1:ndim) = shape(atm%density)
                lbounds(1:ndim) = lbound(atm%density)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'p') then
            if (allocated(atm%p)) then
                loc = C_LOC(atm%p)
                ndim = size(shape(atm%p))
                dims(1:ndim) = shape(atm%p)
                lbounds(1:ndim) = lbound(atm%p)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'p_level') then
            if (allocated(atm%p_level)) then
                loc = C_LOC(atm%p_level)
                ndim = size(shape(atm%p_level))
                dims(1:ndim) = shape(atm%p_level)
                lbounds(1:ndim) = lbound(atm%p_level)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 't') then
            if (allocated(atm%t)) then
                loc = C_LOC(atm%t)
                ndim = size(shape(atm%t))
                dims(1:ndim) = shape(atm%t)
                lbounds(1:ndim) = lbound(atm%t)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 't_level') then
            if (allocated(atm%t_level)) then
                loc = C_LOC(atm%t_level)
                ndim = size(shape(atm%t_level))
                dims(1:ndim) = shape(atm%t_level)
                lbounds(1:ndim) = lbound(atm%t_level)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'r_layer') then
            if (allocated(atm%r_layer)) then
                loc = C_LOC(atm%r_layer)
                ndim = size(shape(atm%r_layer))
                dims(1:ndim) = shape(atm%r_layer)
                lbounds(1:ndim) = lbound(atm%r_layer)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'r_level') then
            if (allocated(atm%r_level)) then
                loc = C_LOC(atm%r_level)
                ndim = size(shape(atm%r_level))
                dims(1:ndim) = shape(atm%r_level)
                lbounds(1:ndim) = lbound(atm%r_level)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'gas_mix_ratio') then
            if (allocated(atm%gas_mix_ratio)) then
                loc = C_LOC(atm%gas_mix_ratio)
                ndim = size(shape(atm%gas_mix_ratio))
                dims(1:ndim) = shape(atm%gas_mix_ratio)
                lbounds(1:ndim) = lbound(atm%gas_mix_ratio)
            else
                loc = C_NULL_PTR
            end if

        else
            field_OK = .FALSE.
        end if

    end function StrAtm_get_real_array

end module StrAtm_C
