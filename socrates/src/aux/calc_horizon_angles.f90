! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Calculate ancillary horizon angles (currently dummy values)

program calc_horizon_angles

  use netcdf
  use realtype_rd, only: RealK
  use rad_ccf, only: pi
  implicit none

! Orography ancillary copied from:
!   '/data/users/lfric/data/ancils/basic-gal/Quagga/C192/n512e_l70/
!     orography/gmted_ramp2/qrparm.orog.ugrid.nc'
  character(*), parameter :: file_orog='qrparm.orog.ugrid.nc'
  character(*), parameter :: dim_name='dim0'
  character(*), parameter :: orog_name='surface_altitude'
  character(*), parameter :: hor_ang_name='horizon_angle'
  character(*), parameter :: hor_asp_name='horizon_aspect'

  integer :: n_horiz_ang = 16
  integer :: n_horiz_layer = 1
  integer :: num_h_ang, num_h_asp

  integer :: i, k, l

  integer :: status, ncid, dimid1, dimid2, dimid3, varid, n_profile
  real(RealK), allocatable :: horizon_angle(:, :)
  real(RealK), allocatable :: horizon_aspect(:, :)

! Open orography ancillary file read/write
  call nf(nf90_open(trim(file_orog), NF90_WRITE, ncid))

! Get number of points
  call nf(nf90_inq_dimid(ncid, dim_name, dimid1))
  call nf(nf90_inquire_dimension(ncid, dimid1, len=n_profile))

! See if horizon angles are already present
  status = nf90_inq_varid(ncid, hor_ang_name, varid)
  if (status == NF90_NOERR) then
    call nf(nf90_close(ncid))
    stop 'Horizon angles already present: will not be overwritten'
  end if

! If not, add them to the file
  num_h_ang = n_horiz_ang*n_horiz_layer
  num_h_asp = n_horiz_ang
  allocate(horizon_angle(n_profile, num_h_ang))
  allocate(horizon_aspect(n_profile, num_h_asp))

! Set dummy values
  ! Set a flat horizon
  horizon_angle = pi/2.0_RealK
  do k = 1, n_horiz_ang
    do l = 1, n_profile
      ! Equal spacing of aspect angles around horizon
      horizon_aspect(l, k) = 2.0_RealK * pi * real(k-1, RealK) &
                           / real(n_horiz_ang, RealK)
    end do
  end do
  
! Add fields to the orog ancillary file
  call nf(nf90_def_dim(ncid, 'num_h_ang', num_h_ang, dimid2))
  call nf(nf90_def_dim(ncid, 'num_h_asp', num_h_asp, dimid3))

  call nf(nf90_def_var(ncid, 'horizon_angle', NF90_DOUBLE, &
    (/ dimid1, dimid2 /), varid))
  call nf(nf90_put_att(ncid, varid, 'standard_name', 'horizon_angle'))
  call nf(nf90_put_att(ncid, varid, 'long_name', 'horizon_angle'))
  call nf(nf90_put_att(ncid, varid, 'units', 'radian'))
  call nf(nf90_put_att(ncid, varid, 'location', 'face'))
  call nf(nf90_put_att(ncid, varid, 'mesh', 'dynamics'))
  call nf(nf90_put_att(ncid, varid, 'online_operation', 'once'))
  call nf(nf90_put_att(ncid, varid, 'coordinates', 'dynamics_face_x dynamics_face_y'))
  call nf(nf90_enddef(ncid))
  call nf(nf90_put_var(ncid, varid, horizon_angle))

  call nf(nf90_redef(ncid))
  call nf(nf90_def_var(ncid, 'horizon_aspect', NF90_DOUBLE, &
    (/ dimid1, dimid3 /), varid))
  call nf(nf90_put_att(ncid, varid, 'standard_name', 'horizon_aspect'))
  call nf(nf90_put_att(ncid, varid, 'long_name', 'horizon_aspect'))
  call nf(nf90_put_att(ncid, varid, 'units', 'radian'))
  call nf(nf90_put_att(ncid, varid, 'location', 'face'))
  call nf(nf90_put_att(ncid, varid, 'mesh', 'dynamics'))
  call nf(nf90_put_att(ncid, varid, 'online_operation', 'once'))
  call nf(nf90_put_att(ncid, varid, 'coordinates', 'dynamics_face_x dynamics_face_y'))
  call nf(nf90_enddef(ncid))
  call nf(nf90_put_var(ncid, varid, horizon_aspect))

  deallocate(horizon_aspect)
  deallocate(horizon_angle)
  call nf(nf90_close(ncid))

contains

  subroutine nf(status)
    use netcdf
    integer, intent(in):: status
    if (status /= NF90_NOERR) then
       write(*,*) 'netCDF-ERROR: ',nf90_strerror(status)
       stop 'Stopped while performing netCDF operation'
    end if
  end subroutine nf

end program calc_horizon_angles
