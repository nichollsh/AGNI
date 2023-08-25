! C to F utility functions for string conversion etc

module UTILITIES_CF

use, intrinsic :: iso_c_binding

implicit none

private


public :: C_strlen, c_to_f_string, copy_f_to_c_string, copy_c_to_f_string
public :: test_double_val, test_double_ref

interface
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! C utility functions for internal use
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Interface to C strlen
    ! see eg http://fortranwiki.org/fortran/show/c_interface_module
    ! Fortran doesn't appear to have a way to say C_char_ptr
    function C_strlen(C_string_ptr) result(length) bind(C, name='strlen')
        use, intrinsic :: iso_c_binding
        implicit none
        type(c_ptr), value, intent(in) :: C_string_ptr
        integer(c_size_t) :: length
    end function C_strlen 
end interface

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Test functions for argument passing etc
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! minimal function to test argument passing from/to C via Julia
    function test_double_val(din) result(dout) bind(C, NAME='PS_test_double_val')
        implicit none
        real(c_double), value :: din
        real(c_double) :: dout       

        dout = din + 1.0

    end function test_double_val

    ! minimal function to test argument passing from/to C via Julia
    ! If din is passed as a C null pointer, Fortran sees this as a missing optional argument
    function test_double_ref(din) result(dout) bind(C, NAME='PS_test_double_ref')
        implicit none
        real(c_double), optional :: din
        real(c_double) :: dout       

        if (PRESENT(din)) then
            dout = din + 1.0
        else
            dout = -1.0
        endif

    end function test_double_ref

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Utility functions for internal use
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! see eg http://fortranwiki.org/fortran/show/c_interface_module
    ! convert a null-terminated C string C_string (as a c_ptr) to a fortran string F_string
    ! returns empty string if C_string is a null pointer
    function c_to_f_string(string_C) result(string_F)
        implicit none   
        character(len=:),allocatable                        :: string_F
        type(c_ptr), intent(in)                             :: string_C
        ! local variables
        character(len=1,kind=c_char), dimension(:), pointer :: p_chars
        integer(c_size_t)                                   :: length, i
        logical(c_bool)                                     :: length_OK
    
        if (.not. c_associated(string_C)) then
            string_F = ''
        else                   
            allocate(character(len=C_strlen(string_C)) :: string_F)
            length_OK = copy_c_to_f_string(string_F, string_C)
        end if
            
    end function c_to_f_string

    ! Copy C string into supplied Fortran string
    ! returns .FALSE. if Fortran string is too short
    function copy_c_to_f_string(string_F, string_C) result(length_OK)
        implicit none
        logical(c_bool)                                     :: length_OK
        character(len=*), intent(out)                       :: string_F
        type(c_ptr), intent(in)                             :: string_C
        ! local variables
        character(len=1,kind=c_char), dimension(:), pointer :: p_chars
        integer(c_size_t)                                   :: length, i
            
        length = min(C_strlen(string_C), len(string_F))
        if (length .LE. len(string_F)) then
            length_OK = .TRUE.
        else
            length_OK = .FALSE.
        end if

        string_F = ''
               
        call c_f_pointer(string_C, p_chars, [C_strlen(string_C)])

        forall (i=1:length)
            string_F(i:i) = p_chars(i)
        end forall
                    
    end function copy_c_to_f_string

    ! copy a Fortran string into supplied C buffer
    ! returns .FALSE. if buffer is too short
    function copy_f_to_c_string(string_C, string_C_len, string_F) result(length_OK)
        implicit none
        logical(c_bool)                 :: length_OK        
        type(c_ptr), intent(in)         :: string_C
        integer(c_size_t), intent(in)   :: string_C_len  ! Max string length,
                                                       ! INCLUDING THE TERMINAL NUL
        character(len=*), intent(in)    :: string_F
        ! local variables
        character(len=1,kind=C_char), dimension(:), pointer :: p_chars
        integer :: i, strlen
      
        strlen = min(len(string_F), string_C_len-1)
        if (strlen .EQ. len(string_F)) then
            length_OK = .TRUE.
        else
            length_OK = .FALSE.
        end if

        call C_F_pointer(string_C, p_chars, [strlen+1])
        forall (i=1:strlen)
          p_chars(i) = string_F(i:i)
        end forall
        p_chars(strlen+1) = C_NULL_char
    end function copy_f_to_c_string

  end module UTILITIES_CF