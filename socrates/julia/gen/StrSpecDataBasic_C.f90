! Autogenerated from Fortran file src/radiance_core/def_spectrum.F90
! svn revision 1226

module StrSpecDataBasic_C

use, intrinsic :: iso_c_binding

use UTILITIES_CF
use def_spectrum

implicit none

private

contains

    function StrSpecDataBasic_get_integer(spectrum_Cptr, field_Cstr, val) result(field_OK) &
        bind(C, NAME='PS_StrSpecDataBasic_get_integer')

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
    
        if (field == 'n_band') then
            val = spectrum%Basic%n_band
        else
            field_OK = .FALSE.
        end if

    end function StrSpecDataBasic_get_integer

    function StrSpecDataBasic_set_integer(spectrum_Cptr, field_Cstr, val) result(field_OK) &
        bind(C, NAME='PS_StrSpecDataBasic_set_integer')

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
    
        if (field == 'n_band') then
            spectrum%Basic%n_band = val
        else
            field_OK = .FALSE.
        end if

    end function StrSpecDataBasic_set_integer

    function StrSpecDataBasic_get_logical_array(spectrum_Cptr, field_Cstr, loc, dims_C, ndim, lbounds_C) result(field_OK) &
        bind(C, NAME='PS_StrSpecDataBasic_get_logical_array')

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
    
        if (field == 'l_present') then
            if (allocated(spectrum%Basic%l_present)) then
                loc = C_LOC(spectrum%Basic%l_present)
                ndim = size(shape(spectrum%Basic%l_present))
                dims(1:ndim) = shape(spectrum%Basic%l_present)
                lbounds(1:ndim) = lbound(spectrum%Basic%l_present)
            else
                loc = C_NULL_PTR
            end if

        else
            field_OK = .FALSE.
        end if

    end function StrSpecDataBasic_get_logical_array

    function StrSpecDataBasic_get_integer_array(spectrum_Cptr, field_Cstr, loc, dims_C, ndim, lbounds_C) result(field_OK) &
        bind(C, NAME='PS_StrSpecDataBasic_get_integer_array')

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
    
        if (field == 'n_band_exclude') then
            if (allocated(spectrum%Basic%n_band_exclude)) then
                loc = C_LOC(spectrum%Basic%n_band_exclude)
                ndim = size(shape(spectrum%Basic%n_band_exclude))
                dims(1:ndim) = shape(spectrum%Basic%n_band_exclude)
                lbounds(1:ndim) = lbound(spectrum%Basic%n_band_exclude)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'index_exclude') then
            if (allocated(spectrum%Basic%index_exclude)) then
                loc = C_LOC(spectrum%Basic%index_exclude)
                ndim = size(shape(spectrum%Basic%index_exclude))
                dims(1:ndim) = shape(spectrum%Basic%index_exclude)
                lbounds(1:ndim) = lbound(spectrum%Basic%index_exclude)
            else
                loc = C_NULL_PTR
            end if

        else
            field_OK = .FALSE.
        end if

    end function StrSpecDataBasic_get_integer_array

    function StrSpecDataBasic_get_real_array(spectrum_Cptr, field_Cstr, loc, dims_C, ndim, lbounds_C) result(field_OK) &
        bind(C, NAME='PS_StrSpecDataBasic_get_real_array')

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
    
        if (field == 'wavelength_long') then
            if (allocated(spectrum%Basic%wavelength_long)) then
                loc = C_LOC(spectrum%Basic%wavelength_long)
                ndim = size(shape(spectrum%Basic%wavelength_long))
                dims(1:ndim) = shape(spectrum%Basic%wavelength_long)
                lbounds(1:ndim) = lbound(spectrum%Basic%wavelength_long)
            else
                loc = C_NULL_PTR
            end if

        else if (field == 'wavelength_short') then
            if (allocated(spectrum%Basic%wavelength_short)) then
                loc = C_LOC(spectrum%Basic%wavelength_short)
                ndim = size(shape(spectrum%Basic%wavelength_short))
                dims(1:ndim) = shape(spectrum%Basic%wavelength_short)
                lbounds(1:ndim) = lbound(spectrum%Basic%wavelength_short)
            else
                loc = C_NULL_PTR
            end if

        else
            field_OK = .FALSE.
        end if

    end function StrSpecDataBasic_get_real_array

end module StrSpecDataBasic_C
