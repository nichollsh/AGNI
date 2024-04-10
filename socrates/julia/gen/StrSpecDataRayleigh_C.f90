! Autogenerated from Fortran file src/radiance_core/def_spectrum.F90
! svn revision 1226

module StrSpecDataRayleigh_C

use, intrinsic :: iso_c_binding

use UTILITIES_CF
use def_spectrum

implicit none

private

contains

    function StrSpecDataRayleigh_get_integer(spectrum_Cptr, field_Cstr, val) result(field_OK) &
        bind(C, NAME='PS_StrSpecDataRayleigh_get_integer')

        implicit none
        type(c_ptr), value, intent(in)          :: spectrum_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        integer(c_int), intent(out) :: val

        ! local variables
        TYPE(StrSpecData), pointer                  :: spectrum
        character(len=:), allocatable           :: field

        ! convert types
        call C_F_POINTER(spectrum_Cptr, spectrum)
        field = c_to_f_string(field_Cstr)

        field_OK = .TRUE.
    
        if (field == 'i_rayleigh_scheme') then
            val = spectrum%Rayleigh%i_rayleigh_scheme
        else if (field == 'n_gas_rayleigh') then
            val = spectrum%Rayleigh%n_gas_rayleigh
        else
            field_OK = .FALSE.
        end if

    end function StrSpecDataRayleigh_get_integer

    function StrSpecDataRayleigh_set_integer(spectrum_Cptr, field_Cstr, val) result(field_OK) &
        bind(C, NAME='PS_StrSpecDataRayleigh_set_integer')

        implicit none
        type(c_ptr), value, intent(in)          :: spectrum_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        integer(c_int), value, intent(in) :: val

        ! local variables
        TYPE(StrSpecData), pointer                  :: spectrum
        character(len=:), allocatable           :: field

        ! convert types
        call C_F_POINTER(spectrum_Cptr, spectrum)
        field = c_to_f_string(field_Cstr)

        field_OK = .TRUE.
    
        if (field == 'i_rayleigh_scheme') then
            spectrum%Rayleigh%i_rayleigh_scheme = val
        else if (field == 'n_gas_rayleigh') then
            spectrum%Rayleigh%n_gas_rayleigh = val
        else
            field_OK = .FALSE.
        end if

    end function StrSpecDataRayleigh_set_integer

    function StrSpecDataRayleigh_get_integer_array(spectrum_Cptr, field_Cstr, loc, dims_C, ndim, lbounds_C) result(field_OK) &
        bind(C, NAME='PS_StrSpecDataRayleigh_get_integer_array')

        implicit none
        type(c_ptr), value, intent(in)          :: spectrum_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        type(c_ptr), intent(out)                :: loc
        type(c_ptr), value, intent(in)          :: dims_C
        integer(c_int), intent(inout)           :: ndim
        type(c_ptr), value, intent(in)          :: lbounds_C

        ! local variables
        TYPE(StrSpecData), pointer             :: spectrum
        character(len=:), allocatable           :: field
        integer(c_int), pointer, dimension(:)   :: dims, lbounds    

        ! convert types
        call C_F_POINTER(spectrum_Cptr, spectrum)
        field = c_to_f_string(field_Cstr)
        call C_F_POINTER(dims_C, dims, [ndim])
        call C_F_POINTER(lbounds_C, lbounds, [ndim])

        field_OK = .TRUE.
    
        if (field == 'index_rayleigh') then
            if (allocated(spectrum%Rayleigh%index_rayleigh)) then
                loc = C_LOC(spectrum%Rayleigh%index_rayleigh)
                ndim = size(shape(spectrum%Rayleigh%index_rayleigh))
                dims(1:ndim) = shape(spectrum%Rayleigh%index_rayleigh)
                lbounds(1:ndim) = lbound(spectrum%Rayleigh%index_rayleigh)
            else
                loc = C_NULL_PTR
            end if

        else
            field_OK = .FALSE.
        end if

    end function StrSpecDataRayleigh_get_integer_array

    function StrSpecDataRayleigh_get_real_array(spectrum_Cptr, field_Cstr, loc, dims_C, ndim, lbounds_C) result(field_OK) &
        bind(C, NAME='PS_StrSpecDataRayleigh_get_real_array')

        implicit none
        type(c_ptr), value, intent(in)          :: spectrum_Cptr
        type(c_ptr), value, intent(in)          :: field_Cstr
        logical(c_bool)                         :: field_OK
        type(c_ptr), intent(out)                :: loc
        type(c_ptr), value, intent(in)          :: dims_C
        integer(c_int), intent(inout)           :: ndim
        type(c_ptr), value, intent(in)          :: lbounds_C

        ! local variables
        TYPE(StrSpecData), pointer             :: spectrum
        character(len=:), allocatable           :: field
        integer(c_int), pointer, dimension(:)   :: dims, lbounds    

        ! convert types
        call C_F_POINTER(spectrum_Cptr, spectrum)
        field = c_to_f_string(field_Cstr)
        call C_F_POINTER(dims_C, dims, [ndim])
        call C_F_POINTER(lbounds_C, lbounds, [ndim])

        field_OK = .TRUE.
    
        if (field == 'rayleigh_coeff') then
            if (allocated(spectrum%Rayleigh%rayleigh_coeff)) then
                loc = C_LOC(spectrum%Rayleigh%rayleigh_coeff)
                ndim = size(shape(spectrum%Rayleigh%rayleigh_coeff))
                dims(1:ndim) = shape(spectrum%Rayleigh%rayleigh_coeff)
                lbounds(1:ndim) = lbound(spectrum%Rayleigh%rayleigh_coeff)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'rayleigh_coeff_gas') then
            if (allocated(spectrum%Rayleigh%rayleigh_coeff_gas)) then
                loc = C_LOC(spectrum%Rayleigh%rayleigh_coeff_gas)
                ndim = size(shape(spectrum%Rayleigh%rayleigh_coeff_gas))
                dims(1:ndim) = shape(spectrum%Rayleigh%rayleigh_coeff_gas)
                lbounds(1:ndim) = lbound(spectrum%Rayleigh%rayleigh_coeff_gas)
            else
                loc = C_NULL_PTR
            end if

        else
            field_OK = .FALSE.
        end if

    end function StrSpecDataRayleigh_get_real_array

end module StrSpecDataRayleigh_C