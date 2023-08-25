! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the variables in the Socrates atm type
!
!------------------------------------------------------------------------------
module socrates_set_atm
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_ATM'
contains

subroutine set_atm(atm, dimen, spectrum, n_profile, n_layer, &
  profile_list, n_layer_stride, n_level_stride, &
  p_layer, t_layer, mass, density, p_level, t_level, r_layer, r_level, &
  p_layer_1d, t_layer_1d, mass_1d, density_1d, p_level_1d, t_level_1d, &
  r_layer_1d, r_level_1d, &
  h2o, co2, o3, n2o, ch4, o2, so2, n2, cfc11, cfc12, cfc113, hcfc22, hfc134a, &
  h2o_1d, co2_1d, o3_1d, n2o_1d, ch4_1d, o2_1d, so2_1d, n2_1d, cfc11_1d, &
  cfc12_1d, cfc113_1d, hcfc22_1d, hfc134a_1d, &
  h2o_mix_ratio, co2_mix_ratio, o3_mix_ratio, n2o_mix_ratio, ch4_mix_ratio, &
  o2_mix_ratio, so2_mix_ratio, n2_mix_ratio, cfc11_mix_ratio, cfc12_mix_ratio, &
  cfc113_mix_ratio, hcfc22_mix_ratio, hfc134a_mix_ratio, &
  l_h2o_well_mixed, l_co2_well_mixed, l_o3_well_mixed, l_n2o_well_mixed, &
  l_ch4_well_mixed, l_o2_well_mixed, l_so2_well_mixed, l_n2_well_mixed, &
  l_cfc11_well_mixed, l_cfc12_well_mixed, l_cfc113_well_mixed, &
  l_hcfc22_well_mixed, l_hfc134a_well_mixed, &
  l_invert, l_profile_last, l_debug, i_profile_debug)

use def_atm,      only: StrAtm, allocate_atm
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use realtype_rd,  only: RealK, RealExt
use gas_list_pcf, only: ip_h2o, ip_co2, ip_o3, ip_n2o, ip_ch4, ip_o2, ip_so2, &
  ip_n2, ip_cfc11, ip_cfc12, ip_cfc113, ip_hcfc22, ip_hfc134a

implicit none


! Atmospheric properties:
type(StrAtm),      intent(out) :: atm

! Dimensions:
type(StrDim),      intent(in)  :: dimen

! Spectral data:
type(StrSpecData), intent(in)  :: spectrum

integer, intent(in) :: n_profile
!   Number of atmospheric profiles for radiation calculations
integer, intent(in) :: n_layer
!   Number of atmospheric layers for radiation calculations
integer, intent(in), optional :: profile_list(:)
!   List of profiles to use from input fields
integer, intent(in), optional :: n_layer_stride
!   Number of layers in input 1d arrays
integer, intent(in), optional :: n_level_stride
!   Number of levels in input 1d arrays

real(RealExt), intent(in), optional :: p_layer(:, :), p_layer_1d(:)
!   Pressure at layer centres
real(RealExt), intent(in), optional :: t_layer(:, :), t_layer_1d(:)
!   Temperature at layer centres
real(RealExt), intent(in), optional :: mass(:, :), mass_1d(:)
!   Mass of layer (kg m-2)
real(RealExt), intent(in), optional :: density(:, :), density_1d(:)
!   Density of layer (kg m-3)
real(RealExt), intent(in), optional :: p_level(:, :), p_level_1d(:)
!   Pressure at layer boundaries
real(RealExt), intent(in), optional :: t_level(:, :), t_level_1d(:)
!   Temperature at layer boundaries
real(RealExt), intent(in), optional :: r_layer(:, :), r_layer_1d(:)
!   Radius (height from centre of planet) at layer centres
real(RealExt), intent(in), optional :: r_level(:, :), r_level_1d(:)
!   Radius (height from centre of planet) at layer boundaries

real(RealExt), intent(in), dimension(:, :), optional :: &
  h2o, co2, o3, n2o, ch4, o2, so2, n2, cfc11, cfc12, cfc113, hcfc22, hfc134a
!   Full field mass mixing ratios

real(RealExt), intent(in), dimension(:), optional :: &
  h2o_1d, co2_1d, o3_1d, n2o_1d, ch4_1d, o2_1d, so2_1d, n2_1d, cfc11_1d, &
  cfc12_1d, cfc113_1d, hcfc22_1d, hfc134a_1d
!   1d mass mixing ratios

real(RealExt), intent(in), optional :: &
  h2o_mix_ratio, co2_mix_ratio, o3_mix_ratio, n2o_mix_ratio, ch4_mix_ratio, &
  o2_mix_ratio, so2_mix_ratio, n2_mix_ratio, cfc11_mix_ratio, cfc12_mix_ratio, &
  cfc113_mix_ratio, hcfc22_mix_ratio, hfc134a_mix_ratio
!   Well mixed mass mixing ratios

logical, intent(in), optional :: &
  l_h2o_well_mixed, l_co2_well_mixed, l_o3_well_mixed, l_n2o_well_mixed, &
  l_ch4_well_mixed, l_o2_well_mixed, l_so2_well_mixed, l_n2_well_mixed, &
  l_cfc11_well_mixed, l_cfc12_well_mixed, l_cfc113_well_mixed, &
  l_hcfc22_well_mixed, l_hfc134a_well_mixed
!   Flag to use the well mixed ratios

logical, intent(in), optional :: l_invert
!   Flag to invert fields in the vertical
logical, intent(in), optional :: l_profile_last
!   Loop over profiles is last in input fields

logical, intent(in), optional :: l_debug
integer, intent(in), optional :: i_profile_debug
!   Options for outputting debugging information

! Local variables.
integer :: i, ii, l, ll, i_gas, list(n_profile)
integer :: layer_offset, level_offset, stride_layer, stride_level
logical :: l_last

call allocate_atm(atm, dimen, spectrum)

! Set up atmosphere grid
atm%n_profile = n_profile
atm%n_layer   = n_layer

if (present(profile_list)) then
  list = profile_list(1:n_profile)
else
  do l=1, n_profile
    list(l) = l
  end do
end if

layer_offset = 0
level_offset = 0
if (present(l_invert)) then
  if (l_invert) then
    ! The layer is indexed using an inverted loop counter
    layer_offset = n_layer + 1
    level_offset = n_layer
  end if
end if

if (present(l_profile_last)) then
  l_last = l_profile_last
else
  l_last = .false.
end if

! Set the number of layers and levels in the 1d arrays
if (present(n_layer_stride)) then
  stride_layer = n_layer_stride
else
  stride_layer = n_layer
end if
if (present(n_level_stride)) then
  stride_level = n_level_stride
else
  stride_level = n_layer + 1
end if

! Set the pressures, temperatures, masses (per square metre) and densities
call set_layer_field(atm%p, p_layer, p_layer_1d)
call set_layer_field(atm%t, t_layer, t_layer_1d)
call set_layer_field(atm%mass, mass, mass_1d)
call set_layer_field(atm%density, density, density_1d)
call set_layer_field(atm%r_layer, r_layer, r_layer_1d)

call set_level_field(atm%p_level, p_level, p_level_1d)
call set_level_field(atm%t_level, t_level, t_level_1d)
call set_level_field(atm%r_level, r_level, r_level_1d)

! Set gas mass mixing ratios
do i_gas=1, spectrum%gas%n_absorb
  select case(spectrum%gas%type_absorb(i_gas))
  case(ip_h2o)
    call set_gas_mix_ratio(h2o, h2o_1d, h2o_mix_ratio, l_h2o_well_mixed)
  case(ip_co2)
    call set_gas_mix_ratio(co2, co2_1d, co2_mix_ratio, l_co2_well_mixed)
  case(ip_o3)
    call set_gas_mix_ratio(o3, o3_1d, o3_mix_ratio, l_o3_well_mixed)
  case(ip_n2o)
    call set_gas_mix_ratio(n2o, n2o_1d, n2o_mix_ratio, l_n2o_well_mixed)
  case(ip_ch4)
    call set_gas_mix_ratio(ch4, ch4_1d, ch4_mix_ratio, l_ch4_well_mixed)
  case(ip_o2)
    call set_gas_mix_ratio(o2, o2_1d, o2_mix_ratio, l_o2_well_mixed)
  case(ip_so2)
    call set_gas_mix_ratio(so2, so2_1d, so2_mix_ratio, l_so2_well_mixed)
  case(ip_cfc11)
    call set_gas_mix_ratio(cfc11, cfc11_1d, cfc11_mix_ratio, l_cfc11_well_mixed)
  case(ip_cfc12)
    call set_gas_mix_ratio(cfc12, cfc12_1d, cfc12_mix_ratio, l_cfc12_well_mixed)
  case(ip_cfc113)
    call set_gas_mix_ratio(cfc113, cfc113_1d, cfc113_mix_ratio, &
                           l_cfc113_well_mixed)
  case(ip_hcfc22)
    call set_gas_mix_ratio(hcfc22, hcfc22_1d, hcfc22_mix_ratio, &
                           l_hcfc22_well_mixed)
  case(ip_hfc134a)
    call set_gas_mix_ratio(hfc134a, hfc134a_1d, hfc134a_mix_ratio, &
                           l_hfc134a_well_mixed)
  case default
    do i=1, n_layer
      do l=1, n_profile
        atm%gas_mix_ratio(l, i, i_gas) = 0.0_RealK
      end do
    end do
  end select
end do

if (present(l_debug)) then
  if (l_debug) then
    if (present(i_profile_debug)) then
      l = i_profile_debug
    else
      l = 1
    end if
    write(1000+l,'(A)') 'PRESS(PA) TEMP(K) MASS(KGM-2) DENSITY(KGM-3) GASES'
    do i=1, n_layer
      write(1000+l,'(7(1pe16.8))') &
        atm%p(l, i), atm%t(l, i), atm%mass(l, i), atm%density(l, i), &
        atm%gas_mix_ratio(l, i, 1:MIN(spectrum%gas%n_absorb,3))
    end do
    write(1000+l,'(A)') 'PLEV(PA) TEMP(K)'
    do i=0, n_layer
      write(1000+l,'(2(1pe16.8))') atm%p_level(l, i), atm%t_level(l, i)
    end do
  end if
end if

contains

  subroutine set_layer_field(out_field, full_field, oned_field)
    implicit none

    real(RealK), intent(inout) :: out_field(:, :)
!     Output field
    real(RealExt), intent(in), optional :: full_field(:, :)
!     Full field variable
    real(RealExt), intent(in), optional :: oned_field(:)
!     One-dimensional variable

    if (present(full_field)) then
      if (l_last) then
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            out_field(l, i) = real(full_field(ii, list(l)), RealK)
          end do
        end do
      else
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            out_field(l, i) = real(full_field(list(l), ii), RealK)
          end do
        end do
      end if
    else if (present(oned_field)) then
      if (l_last) then
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            ll = stride_layer*(list(l)-1) + ii
            out_field(l, i) = real(oned_field(ll), RealK)
          end do
        end do
      else
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            ll = n_profile*(ii-1) + list(l)
            out_field(l, i) = real(oned_field(ll), RealK)
          end do
        end do
      end if
    end if
  end subroutine set_layer_field


  subroutine set_level_field(out_field, full_field, oned_field)
    implicit none

    real(RealK), intent(inout) :: out_field(:, 0:)
!     Output field
    real(RealExt), intent(in), optional :: full_field(:, :)
!     Full field variable
    real(RealExt), intent(in), optional :: oned_field(:)
!     One-dimensional variable

    if (present(full_field)) then
      if (l_last) then
        do i=0, n_layer
          ii = abs(level_offset-i) + 1
          do l=1, n_profile
            out_field(l, i) = real(full_field(ii, list(l)), RealK)
          end do
        end do
      else
        do i=0, n_layer
          ii = abs(level_offset-i) + 1
          do l=1, n_profile
            out_field(l, i) = real(full_field(list(l), ii), RealK)
          end do
        end do
      end if
    else if (present(oned_field)) then
      if (l_last) then
        do i=0, n_layer
          ii = abs(level_offset-i) + 1
          do l=1, n_profile
            ll = stride_level*(list(l)-1) + ii
            out_field(l, i) = real(oned_field(ll), RealK)
          end do
        end do
      else
        do i=0, n_layer
          ii = abs(level_offset-i) + 1
          do l=1, n_profile
            ll = n_profile*(ii-1) + list(l)
            out_field(l, i) = real(oned_field(ll), RealK)
          end do
        end do
      end if
    end if
  end subroutine set_level_field


  subroutine set_gas_mix_ratio(full_field, oned_field, mix_ratio, l_well_mixed)
    implicit none

    real(RealExt), intent(in), optional :: full_field(:, :)
!     Full field mass mixing ratio
    real(RealExt), intent(in), optional :: oned_field(:)
!     One-dimensional mass mixing ratio
    real(RealExt), intent(in), optional :: mix_ratio
!     Well mixed mass mixing ratio
    logical, intent(in), optional :: l_well_mixed
!     Flag to use the well mixed ratio

    integer :: i_select
    integer, parameter :: ip_zero = 0
    integer, parameter :: ip_full_field = 1
    integer, parameter :: ip_oned_field = 2
    integer, parameter :: ip_well_mixed = 3

    if (present(full_field)) then
      i_select=ip_full_field
    else if (present(oned_field)) then
      i_select=ip_oned_field
    else if (present(mix_ratio)) then
      i_select=ip_well_mixed
    else
      i_select=ip_zero
    end if
    if ( present(mix_ratio) .and. present(l_well_mixed) ) then
      if (l_well_mixed) then
        i_select=ip_well_mixed
      end if
    end if

    select case(i_select)
    case(ip_full_field)
      if (l_last) then
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            atm%gas_mix_ratio(l, i, i_gas) &
              = max(real(full_field(ii, list(l)), RealK), 0.0_RealK)
          end do
        end do
      else
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            atm%gas_mix_ratio(l, i, i_gas) &
              = max(real(full_field(list(l), ii), RealK), 0.0_RealK)
          end do
        end do
      end if
    case(ip_oned_field)
      if (l_last) then
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            ll = stride_layer*(list(l)-1) + ii
            atm%gas_mix_ratio(l, i, i_gas) &
              = max(real(oned_field(ll), RealK), 0.0_RealK)
          end do
        end do
      else
        do i=1, n_layer
          ii = abs(layer_offset-i)
          do l=1, n_profile
            ll = n_profile*(ii-1) + list(l)
            atm%gas_mix_ratio(l, i, i_gas) &
              = max(real(oned_field(ll), RealK), 0.0_RealK)
          end do
        end do
      end if
    case(ip_well_mixed)
      do i=1, n_layer
        do l=1, n_profile
          atm%gas_mix_ratio(l, i, i_gas) = real(mix_ratio, RealK)
        end do
      end do
    case(ip_zero)
      do i=1, n_layer
        do l=1, n_profile
          atm%gas_mix_ratio(l, i, i_gas) = 0.0_RealK
        end do
      end do
    end select
  end subroutine set_gas_mix_ratio

end subroutine set_atm
end module socrates_set_atm
