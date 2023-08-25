! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read the external data for the CKD foreign continuum.
!
SUBROUTINE set_extern_ckd_frn_data &
!
(ierr)
!
! Method:
!   The name of the file is requested. The file is opened
!   and read. Directives encoutered are processed as the
!   data are read.
!
!
! Modules used:
  USE realtype_rd
  USE def_std_io_icf
  USE ckd_extern_data
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
!
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
!
!
! Local variables.
  CHARACTER (LEN=80) :: line
!   Line read from input file
  INTEGER :: iu_ckd_data
!   Unit number for reading CKD data
  INTEGER :: ios
!   I/O error flag
  INTEGER :: k
!   Loop variable
!
! Subroutines called:
  EXTERNAL &
      open_file_in, &
      get_free_unit 
!
!
!
!
  CALL get_free_unit(ierr, iu_ckd_data)
!
!
!
  CALL open_file_in(ierr, iu_ckd_data, &
    'Enter the name of the file containing the foreign CKD data.')
!
  DO 
    READ(iu_ckd_data, "(a)", IOSTAT=ios) line
    IF (ios /= 0) THEN
      WRITE(iu_err, "(/a)") "*** Error reading foreign continuum."
      ierr = i_err_fatal
      RETURN
    ENDIF
    IF (line(1:11) == "*BEGIN_DATA") EXIT
  ENDDO
!
! Set up the size of the table
  READ(iu_ckd_data, "(3(/, T38, E12.5))") &
    c_foreign_h2o_296%table_start, &
    c_foreign_h2o_296%table_end, &
    c_foreign_h2o_296%table_inc
  READ(iu_ckd_data, "(T38, i5)") c_foreign_h2o_296%n_freq
  ALLOCATE(c_foreign_h2o_296%c(c_foreign_h2o_296%n_freq))
  READ(iu_ckd_data, "(//, 4(3X, 1E12.5))", IOSTAT=ios) &
    (c_foreign_h2o_296%c(k), k=1, c_foreign_h2o_296%n_freq)
  IF (ios /= 0) THEN
    WRITE(iu_err, "(/a)") "*** Error reading foreign continuum."
    ierr = i_err_fatal
    RETURN
  ENDIF
!
  CLOSE(iu_ckd_data)
!
!
!
  RETURN
END SUBROUTINE set_extern_ckd_frn_data
