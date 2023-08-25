! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to calculate the solar flux at the top of the atmosphere.
!
! Method:
!	A solar constant and a day of the year are supplied. The
!	solar flux at the top of the atmosphere is calculated 
!	accordingly using the same algorithm as in the Unified Model.
!
!- ---------------------------------------------------------------------
      PROGRAM solar_inc
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE rad_ccf, ONLY: pi, eccentricity, length_year, day_perihelion
      USE def_std_io_icf
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!
!     Declaration of variables.
      INTEGER
     &    ierr
!           Error flag
      INTEGER
     &    ios
!           I/O error flag
      REAL  (RealK) ::
     &    day
!           Day of year
     &  , solar_constant
!           User's mean value of solar constant
     &  , solar_constant_local
!           Local solar constant
      REAL  (RealK) ::
     &    m
!           Mean anomaly
     &  , v
!           True anomaly
!
      data ierr/i_normal/
!
!
!     Obtain day of year and solar constant.
      WRITE(iu_stdout, '(/a)')
     &  'enter the day of the year and the mean solar constant.'
1     read(iu_stdin, *, iostat=ios) day, solar_constant
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ illegal input: please re-enter.'
        goto 1
      ENDIF
!
      m=2.0_RealK*pi*(day-day_perihelion)/length_year
      v=m+(2.0_RealK*eccentricity
     &  -0.25_RealK*eccentricity**3)*sin(m)
     &  +0.8_RealK*eccentricity**2*sin(2.0_RealK*m)
     &  +(13.0_RealK*eccentricity**3
     &  /12.0_RealK)*sin(3.0_RealK*m)
      solar_constant_local=solar_constant*(
     &   (1.0_RealK+0.5_RealK*eccentricity**2)
     &   *(1.0_RealK+eccentricity*cos(v))
     &   /(1.0_RealK-eccentricity**2))**2
!
      WRITE(iu_stdout, '(/a, 1pe12.5)') 
     &  'solar irradiance at the toa = ', solar_constant_local
!
!
!
      STOP
      END
