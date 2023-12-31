! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History:
! 05/01/15  Dustin Swales - Original version
! 04/04/18  Rodrigo Guzman- Added CALIOP-like Ground LIDar routines (GLID)
! 10/04/18  Rodrigo Guzman- Added ATLID-like (EarthCare) lidar routines (ATLID)
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module cosp_optics
  USE COSP_KINDS, ONLY: wp,dp
  USE COSP_MATH_CONSTANTS,  ONLY: pi
  USE COSP_PHYS_CONSTANTS,  ONLY: rholiq,km,rd,grav
  USE MOD_MODIS_SIM,        ONLY: get_g_nir,get_ssa_nir,phaseIsLiquid,phaseIsIce
  implicit none
  
  real(wp),parameter ::        & !
       ice_density   = 0.93_wp   ! Ice density used in MODIS phase partitioning

  interface cosp_simulator_optics
     module procedure cosp_simulator_optics2D, cosp_simulator_optics3D
  end interface cosp_simulator_optics
  
contains
  ! ##########################################################################
  !                          COSP_SIMULATOR_OPTICS
  !
  ! Used by: ISCCP, MISR and MODIS simulators
  ! ##########################################################################
  subroutine cosp_simulator_optics2D(dim1,dim2,dim3,flag,varIN1,varIN2,varOUT)
    ! INPUTS
    integer,intent(in) :: &
         dim1,   & ! Dimension 1 extent (Horizontal)
         dim2,   & ! Dimension 2 extent (Subcolumn)
         dim3      ! Dimension 3 extent (Vertical)
    real(wp),intent(in),dimension(dim1,dim2,dim3) :: &
         flag      ! Logical to determine the of merge var1IN and var2IN
    real(wp),intent(in),dimension(dim1,     dim3) :: &
         varIN1, & ! Input field 1
         varIN2    ! Input field 2
    ! OUTPUTS
    real(wp),intent(out),dimension(dim1,dim2,dim3) :: &
         varOUT    ! Merged output field
    ! LOCAL VARIABLES
    integer :: j
    
    varOUT(1:dim1,1:dim2,1:dim3) = 0._wp
    do j=1,dim2
       where(flag(:,j,:) .eq. 1)
          varOUT(:,j,:) = varIN2
       endwhere
       where(flag(:,j,:) .eq. 2)
          varOUT(:,j,:) = varIN1
       endwhere
    enddo
  end subroutine cosp_simulator_optics2D
  subroutine cosp_simulator_optics3D(dim1,dim2,dim3,flag,varIN1,varIN2,varOUT)
    ! INPUTS
    integer,intent(in) :: &
         dim1,   & ! Dimension 1 extent (Horizontal)
         dim2,   & ! Dimension 2 extent (Subcolumn)
         dim3      ! Dimension 3 extent (Vertical)
    real(wp),intent(in),dimension(dim1,dim2,dim3) :: &
         flag      ! Logical to determine the of merge var1IN and var2IN
    real(wp),intent(in),dimension(dim1,dim2,dim3) :: &
         varIN1, & ! Input field 1
         varIN2    ! Input field 2
    ! OUTPUTS
    real(wp),intent(out),dimension(dim1,dim2,dim3) :: &
         varOUT    ! Merged output field
    
    varOUT(1:dim1,1:dim2,1:dim3) = 0._wp
   where(flag(:,:,:) .eq. 1)
       varOUT(:,:,:) = varIN2
    endwhere
    where(flag(:,:,:) .eq. 2)
       varOUT(:,:,:) = varIN1
    endwhere
    
  end subroutine cosp_simulator_optics3D
  
  ! ##############################################################################
  !                           MODIS_OPTICS_PARTITION
  !
  ! For the MODIS simulator, there are times when only a sinlge optical depth
  ! profile, cloud-ice and cloud-water are provided. In this case, the optical
  ! depth is partitioned by phase.
  ! ##############################################################################
  subroutine MODIS_OPTICS_PARTITION(npoints,nlev,ncolumns,cloudWater,cloudIce,waterSize, &
                                    iceSize,tau,tauL,tauI)
    ! INPUTS
    INTEGER,intent(in) :: &
         npoints,   & ! Number of horizontal gridpoints
         nlev,      & ! Number of levels
         ncolumns     ! Number of subcolumns
    REAL(wp),intent(in),dimension(npoints,nlev,ncolumns) :: &
         cloudWater, & ! Subcolumn cloud water content
         cloudIce,   & ! Subcolumn cloud ice content
         waterSize,  & ! Subcolumn cloud water effective radius
         iceSize,    & ! Subcolumn cloud ice effective radius
         tau           ! Optical thickness
    
    ! OUTPUTS
    real(wp),intent(out),dimension(npoints,nlev,ncolumns) :: &
         tauL,       & ! Partitioned liquid optical thickness.
         tauI          ! Partitioned ice optical thickness.
    ! LOCAL VARIABLES
    real(wp),dimension(nlev,ncolumns) :: fracL
    integer                           :: i
    
    
    do i=1,npoints
       where(cloudIce(i,:, :) <= 0.) 
          fracL(:, :) = 1._wp
       elsewhere
          where (cloudWater(i,:, :) <= 0.) 
             fracL(:, :) = 0._wp
          elsewhere 
             ! Geometic optics limit - tau as LWP/re  (proportional to LWC/re) 
             fracL(:, :) = (cloudWater(i,:, :)/waterSize(i,:, :)) / &
                  (cloudWater(i,:, :)/waterSize(i,:, :) + cloudIce(i,:, :)/(ice_density * iceSize(i,:, :)) ) 
          end where
       end where
       tauL(i,:, :) = fracL(:, :) * tau(i,:, :) 
       tauI(i,:, :) = tau(i,:, :) - tauL(i,:, :)
    enddo
    
  end subroutine MODIS_OPTICS_PARTITION
  ! ########################################################################################
  !                                   MODIS_OPTICS
  ! 
  ! ########################################################################################
  subroutine modis_optics(nPoints,nLevels,nSubCols,tauLIQ,sizeLIQ,tauICE,sizeICE,fracLIQ, g, w0)
    ! INPUTS
    integer, intent(in)                                      :: nPoints,nLevels,nSubCols
    real(wp),intent(in),dimension(nPoints,nSubCols,nLevels)  :: tauLIQ, sizeLIQ, tauICE, sizeICE
    ! OUTPUTS
    real(wp),intent(out),dimension(nPoints,nSubCols,nLevels) :: g,w0,fracLIQ
    ! LOCAL VARIABLES
    real(wp), dimension(nLevels)            :: water_g, water_w0, ice_g, ice_w0,tau
    integer :: i,j
    
    ! Initialize
    g(1:nPoints,1:nSubCols,1:nLevels)  = 0._wp
    w0(1:nPoints,1:nSubCols,1:nLevels) = 0._wp
    
    do j =1,nPoints
       do i=1,nSubCols
          water_g(1:nLevels)  = get_g_nir(  phaseIsLiquid, sizeLIQ(j,i,1:nLevels)) 
          water_w0(1:nLevels) = get_ssa_nir(phaseIsLiquid, sizeLIQ(j,i,1:nLevels))
          ice_g(1:nLevels)    = get_g_nir(  phaseIsIce,    sizeICE(j,i,1:nLevels))
          ice_w0(1:nLevels)   = get_ssa_nir(phaseIsIce,    sizeICE(j,i,1:nLevels))
          
          ! Combine ice and water optical properties
          tau(1:nLevels) = tauICE(j,i,1:nLevels) + tauLIQ(j,i,1:nLevels) 
          where (tau(1:nLevels) > 0) 
             w0(j,i,1:nLevels) = (tauLIQ(j,i,1:nLevels)*water_w0(1:nLevels) + tauICE(j,i,1:nLevels) *ice_w0(1:nLevels)) / &
                                  (tau(1:nLevels))
             g(j,i,1:nLevels) = (tauLIQ(j,i,1:nLevels)*water_g(1:nLevels)*water_w0(1:nLevels) + tauICE(j,i,1:nLevels) * &
                  ice_g(1:nLevels) * ice_w0(1:nLevels)) / (w0(j,i,1:nLevels) * tau(1:nLevels))
          end where
       enddo
    enddo
    
    ! Compute the total optical thickness and the proportion due to liquid in each cell
    do i=1,npoints
       where(tauLIQ(i,1:nSubCols,1:nLevels) + tauICE(i,1:nSubCols,1:nLevels) > 0.) 
          fracLIQ(i,1:nSubCols,1:nLevels) = tauLIQ(i,1:nSubCols,1:nLevels)/ &
               (tauLIQ(i,1:nSubCols,1:nLevels) + tauICE(i,1:nSubCols,1:nLevels))
       elsewhere
          fracLIQ(i,1:nSubCols,1:nLevels) = 0._wp
       end  where
    enddo
    
  end subroutine modis_optics
  
  ! ######################################################################################
  ! SUBROUTINE lidar_optics
  ! ######################################################################################
  subroutine lidar_optics(npoints, ncolumns, nlev, npart, ice_type, lidar_freq, lground, &
       q_lsliq, q_lsice, q_cvliq, q_cvice, ls_radliq, ls_radice, cv_radliq, cv_radice,   &
       pres, presf, temp, beta_mol, betatot, tau_mol, tautot, tautot_S_liq, tautot_S_ice,&
       betatot_ice, betatot_liq, tautot_ice, tautot_liq)

    ! ####################################################################################
    ! NOTE: Using "grav" from cosp_constants.f90, instead of grav=9.81, introduces
    ! changes of up to 2% in atb532 adn 0.003% in parasolRefl and lidarBetaMol532. 
    ! This also results in  small changes in the joint-histogram, cfadLidarsr532.
    ! ####################################################################################
    
    ! INPUTS
    INTEGER,intent(in) :: & 
         npoints,      & ! Number of gridpoints
         ncolumns,     & ! Number of subcolumns
         nlev,         & ! Number of levels
         npart,        & ! Number of cloud meteors (stratiform_liq, stratiform_ice, conv_liq, conv_ice). 
         ice_type,     & ! Ice particle shape hypothesis (0 for spheres, 1 for non-spherical)
         lidar_freq      ! Lidar frequency (nm). Use to change between lidar platforms
    logical,intent(in) :: &
         lground         ! True for ground-based lidar
    REAL(WP),intent(in),dimension(npoints,nlev) :: &
         temp,         & ! Temperature of layer k
         pres,         & ! Pressure at full levels
         ls_radliq,    & ! Effective radius of LS liquid particles (meters)
         ls_radice,    & ! Effective radius of LS ice particles (meters)
         cv_radliq,    & ! Effective radius of CONV liquid particles (meters)
         cv_radice       ! Effective radius of CONV ice particles (meters)
    REAL(WP),intent(in),dimension(npoints,ncolumns,nlev) :: &
         q_lsliq,      & ! LS sub-column liquid water mixing ratio (kg/kg)
         q_lsice,      & ! LS sub-column ice water mixing ratio (kg/kg)
         q_cvliq,      & ! CONV sub-column liquid water mixing ratio (kg/kg)
         q_cvice         ! CONV sub-column ice water mixing ratio (kg/kg)
    REAL(WP),intent(in),dimension(npoints,nlev+1) :: &
         presf           ! Pressure at half levels
    
    ! OUTPUTS
    REAL(WP),intent(out),dimension(npoints,ncolumns,nlev)       :: &
         betatot,        & ! 
         tautot            ! Optical thickess integrated from top
    REAL(WP),intent(out),dimension(npoints,nlev) :: &
         beta_mol,       & ! Molecular backscatter coefficient
         tau_mol           ! Molecular optical depth
    ! OUTPUTS (optional)
    REAL(WP),optional,intent(out),dimension(npoints,ncolumns) :: &
         tautot_S_liq,   & ! TOA optical depth for liquid
         tautot_S_ice      ! TOA optical depth for ice
    REAL(WP),optional,intent(out),dimension(npoints,ncolumns,nlev)       :: &
         betatot_ice,    & ! Backscatter coefficient for ice particles
         betatot_liq,    & ! Backscatter coefficient for liquid particles
         tautot_ice,     & ! Total optical thickness of ice
         tautot_liq        ! Total optical thickness of liq
    
    ! LOCAL VARIABLES
    REAL(WP),dimension(npart)              :: rhopart
    REAL(WP),dimension(npart,5)            :: polpart 
    REAL(WP),dimension(npoints,nlev)       :: rhoair,alpha_mol
    REAL(WP),dimension(npoints,nlev+1)     :: zheight          
    REAL(WP),dimension(npoints,nlev,npart) :: rad_part,kp_part,qpart,alpha_part,tau_part
    real(wp)                               :: Cmol,rdiffm,x
    logical                                :: lparasol,lphaseoptics
    INTEGER                                :: i,j,k,icol,zi,zf,zinc,zoffset
    
    ! Local data
    REAL(WP),PARAMETER :: rhoice     = 0.5e+03    ! Density of ice (kg/m3) 
    REAL(WP),PARAMETER :: Cmol_532nm = 6.2446e-32 ! Wavelength dependent
    REAL(WP),PARAMETER :: Cmol_355nm = 3.2662e-31! Wavelength dependent
    REAL(WP),PARAMETER :: rdiffm_532nm = 0.7_wp     ! Multiple scattering correction parameter
    REAL(WP),PARAMETER :: rdiffm_355nm = 0.6_wp     ! Multiple scattering correction parameter
    REAL(WP),PARAMETER :: Qscat      = 2.0_wp     ! Particle scattering efficiency at 532 nm
    ! Local indicies for large-scale and convective ice and liquid 
    INTEGER,PARAMETER  :: INDX_LSLIQ  = 1
    INTEGER,PARAMETER  :: INDX_LSICE  = 2
    INTEGER,PARAMETER  :: INDX_CVLIQ  = 3
    INTEGER,PARAMETER  :: INDX_CVICE  = 4
    
    ! Polarized optics parameterization
    ! Polynomial coefficients for spherical liq/ice particles derived from Mie theory.
    ! Polynomial coefficients for non spherical particles derived from a composite of
    ! Ray-tracing theory for large particles (e.g. Noel et al., Appl. Opt., 2001)
    ! and FDTD theory for very small particles (Yang et al., JQSRT, 2003).
    ! We repeat the same coefficients for LS and CONV cloud to make code more readable
    REAL(WP),PARAMETER,dimension(5) :: &
         polpartCVLIQ  = (/ 2.6980e-8_wp,  -3.7701e-6_wp,  1.6594e-4_wp,    -0.0024_wp,    0.0626_wp/), &
         polpartLSLIQ  = (/ 2.6980e-8_wp,  -3.7701e-6_wp,  1.6594e-4_wp,    -0.0024_wp,    0.0626_wp/), &
         polpartCVICE0 = (/-1.0176e-8_wp,   1.7615e-6_wp, -1.0480e-4_wp,     0.0019_wp,    0.0460_wp/), &
         polpartLSICE0 = (/-1.0176e-8_wp,   1.7615e-6_wp, -1.0480e-4_wp,     0.0019_wp,    0.0460_wp/), &
         polpartCVICE1 = (/ 1.3615e-8_wp, -2.04206e-6_wp, 7.51799e-5_wp, 0.00078213_wp, 0.0182131_wp/), &
         polpartLSICE1 = (/ 1.3615e-8_wp, -2.04206e-6_wp, 7.51799e-5_wp, 0.00078213_wp, 0.0182131_wp/)
    ! ##############################################################################

    ! Which LIDAR frequency are we using?
    if (lidar_freq .eq. 355) then
       Cmol   = Cmol_355nm
       rdiffm = rdiffm_355nm
    endif
    if (lidar_freq .eq. 532) then
       Cmol   = Cmol_532nm
       rdiffm = rdiffm_532nm
    endif
    
    ! Do we need to generate optical inputs for Parasol simulator?
    lparasol = .false.
    if (present(tautot_S_liq) .and. present(tautot_S_ice)) lparasol = .true.
    
    ! Are optical-depths and backscatter coefficients for ice and liquid requested?
    lphaseoptics=.false.
    if (present(betatot_ice) .and. present(betatot_liq) .and. present(tautot_liq) .and. &
         present(tautot_ice)) lphaseoptics=.true.

    ! Is this lidar spaceborne (default) or ground-based (lground=.true.)?
    zi   = 2
    zf   = nlev
    zinc = 1
    zoffset = -1
    if (lground) then
       zi   = nlev-1
       zf   = 1
       zinc = -1
       zoffset = 1
    endif
    
    ! Liquid/ice particles
    rhopart(INDX_LSLIQ)  = rholiq
    rhopart(INDX_LSICE)  = rhoice
    rhopart(INDX_CVLIQ)  = rholiq
    rhopart(INDX_CVICE)  = rhoice
    
    ! LS and CONV Liquid water coefficients
    polpart(INDX_LSLIQ,1:5)  = polpartLSLIQ
    polpart(INDX_CVLIQ,1:5)  = polpartCVLIQ
    
    ! LS and CONV Ice water coefficients
    if (ice_type .eq. 0) then
       polpart(INDX_LSICE,1:5) = polpartLSICE0
       polpart(INDX_CVICE,1:5) = polpartCVICE0
    endif
    if (ice_type .eq. 1) then
       polpart(INDX_LSICE,1:5) = polpartLSICE1
       polpart(INDX_CVICE,1:5) = polpartCVICE1
    endif
    
    ! Effective radius particles:
    rad_part(1:npoints,1:nlev,INDX_LSLIQ)  = ls_radliq(1:npoints,1:nlev)
    rad_part(1:npoints,1:nlev,INDX_LSICE)  = ls_radice(1:npoints,1:nlev)
    rad_part(1:npoints,1:nlev,INDX_CVLIQ)  = cv_radliq(1:npoints,1:nlev)
    rad_part(1:npoints,1:nlev,INDX_CVICE)  = cv_radice(1:npoints,1:nlev)    
    rad_part(1:npoints,1:nlev,1:npart)     = MAX(rad_part(1:npoints,1:nlev,1:npart),0._wp)
    rad_part(1:npoints,1:nlev,1:npart)     = MIN(rad_part(1:npoints,1:nlev,1:npart),70.0e-6_wp)
    
    ! Density (clear-sky air)
    rhoair(1:npoints,1:nlev) = pres(1:npoints,1:nlev)/(rd*temp(1:npoints,1:nlev))
    
    ! Altitude at half pressure levels:
    zheight(1:npoints,nlev+1) = 0._wp
    do k=nlev,1,-1
       zheight(1:npoints,k) = zheight(1:npoints,k+1) &
            -(presf(1:npoints,k)-presf(1:npoints,k+1))/(rhoair(1:npoints,k)*grav)
    enddo
    
    ! ##############################################################################
    ! *) Molecular alpha, beta and optical thickness
    ! ##############################################################################
    
    beta_mol(1:npoints,1:nlev)  = pres(1:npoints,1:nlev)/km/temp(1:npoints,1:nlev)*Cmol
    alpha_mol(1:npoints,1:nlev) = 8._wp*pi/3._wp * beta_mol(1:npoints,1:nlev)
    
    ! Optical thickness of each layer (molecular)  
    tau_mol(1:npoints,1:nlev) = alpha_mol(1:npoints,1:nlev)*(zheight(1:npoints,1:nlev)-&
         zheight(1:npoints,2:nlev+1))
             
    ! Optical thickness from TOA to layer k (molecular)
    DO k = zi,zf,zinc
       tau_mol(1:npoints,k) = tau_mol(1:npoints,k) + tau_mol(1:npoints,k+zoffset)
    ENDDO    

    do k = 1,ncolumns
      betatot(1:npoints,k,1:nlev) = beta_mol(1:npoints,1:nlev)
      tautot (1:npoints,k,1:nlev) = tau_mol (1:npoints,1:nlev)
    enddo
    if (lphaseoptics) then
       betatot_liq(1:npoints,1:ncolumns,1:nlev) = betatot(1:npoints,1:ncolumns,1:nlev)
       betatot_ice(1:npoints,1:ncolumns,1:nlev) = betatot(1:npoints,1:ncolumns,1:nlev)
       tautot_liq (1:npoints,1:ncolumns,1:nlev) = tautot(1:npoints,1:ncolumns,1:nlev)
       tautot_ice (1:npoints,1:ncolumns,1:nlev) = tautot(1:npoints,1:ncolumns,1:nlev)      
    endif

    ! ##############################################################################
    ! *) Particles alpha, beta and optical thickness
    ! ##############################################################################
    ! Polynomials kp_lidar derived from Mie theory
    do i = 1, npart
      do j = 1, nlev
        do k = 1, npoints
          if (rad_part(k,j,i) .gt. 0._wp) then
            x = rad_part(k,j,i)*1.0e6_wp
            kp_part(k,j,i) = &
               polpart(i,1)*x**4 &
               + polpart(i,2)*x**3 &
               + polpart(i,3)*x**2 &
               + polpart(i,4)*x &
               + polpart(i,5)
          else
            kp_part(k,j,i) = 0._wp
          endif
        enddo    
      enddo    
    enddo    

    ! Initialize (if necessary)
    if (lparasol) then
       tautot_S_liq(1:npoints,1:ncolumns) = 0._wp
       tautot_S_ice(1:npoints,1:ncolumns) = 0._wp
    endif

    ! Loop over all subcolumns
    do icol=1,ncolumns
       ! ##############################################################################
       ! Mixing ratio particles in each subcolum
       ! ##############################################################################
       qpart(1:npoints,1:nlev,INDX_LSLIQ) =  q_lsliq(1:npoints,icol,1:nlev)
       qpart(1:npoints,1:nlev,INDX_LSICE)  = q_lsice(1:npoints,icol,1:nlev)
       qpart(1:npoints,1:nlev,INDX_CVLIQ)  = q_cvliq(1:npoints,icol,1:nlev)
       qpart(1:npoints,1:nlev,INDX_CVICE)  = q_cvice(1:npoints,icol,1:nlev)

       ! ##############################################################################
       ! Alpha and optical thickness (particles)
       ! ##############################################################################
       ! Alpha of particles in each subcolumn:
       do i = 1, npart
          where (rad_part(1:npoints,1:nlev,i) .gt. 0.0)
             alpha_part(1:npoints,1:nlev,i) = 3._wp/4._wp * Qscat &
                  * rhoair(1:npoints,1:nlev) * qpart(1:npoints,1:nlev,i) &
                  / (rhopart(i) * rad_part(1:npoints,1:nlev,i) )
          elsewhere
             alpha_part(1:npoints,1:nlev,i) = 0._wp
          endwhere
       enddo
       
       ! Optical thicknes
       tau_part(1:npoints,1:nlev,1:npart) = rdiffm * alpha_part(1:npoints,1:nlev,1:npart)
       do i = 1, npart
          ! Optical thickness of each layer (particles)
          tau_part(1:npoints,1:nlev,i) = tau_part(1:npoints,1:nlev,i) &
               & * (zheight(1:npoints,1:nlev)-zheight(1:npoints,2:nlev+1) )
          ! Optical thickness from TOA to layer k (particles)
          do k=zi,zf,zinc
             tau_part(1:npoints,k,i) = tau_part(1:npoints,k,i) + tau_part(1:npoints,k+zoffset,i)
          enddo
       enddo
 
       ! ##############################################################################
       ! Beta and optical thickness (total=molecular + particules)
       ! ##############################################################################
       
       DO i = 1, npart
          betatot(1:npoints,icol,1:nlev) = betatot(1:npoints,icol,1:nlev) + &
               kp_part(1:npoints,1:nlev,i)*alpha_part(1:npoints,1:nlev,i)
          tautot(1:npoints,icol,1:nlev) = tautot(1:npoints,icol,1:nlev)  + &
              tau_part(1:npoints,1:nlev,i)
       ENDDO
       
       ! ##############################################################################
       ! Beta and optical thickness (liquid/ice)
       ! ##############################################################################
       if (lphaseoptics) then
          ! Ice
          betatot_ice(1:npoints,icol,1:nlev) = betatot_ice(1:npoints,icol,1:nlev)+ &
               kp_part(1:npoints,1:nlev,INDX_LSICE)*alpha_part(1:npoints,1:nlev,INDX_LSICE)+ &
               kp_part(1:npoints,1:nlev,INDX_CVICE)*alpha_part(1:npoints,1:nlev,INDX_CVICE)
          tautot_ice(1:npoints,icol,1:nlev) = tautot_ice(1:npoints,icol,1:nlev)  + &
               tau_part(1:npoints,1:nlev,INDX_LSICE) + &
               tau_part(1:npoints,1:nlev,INDX_CVICE)
          
          ! Liquid
          betatot_liq(1:npoints,icol,1:nlev) = betatot_liq(1:npoints,icol,1:nlev)+ &
               kp_part(1:npoints,1:nlev,INDX_LSLIQ)*alpha_part(1:npoints,1:nlev,INDX_LSLIQ)+ &
               kp_part(1:npoints,1:nlev,INDX_CVLIQ)*alpha_part(1:npoints,1:nlev,INDX_CVLIQ)
          tautot_liq(1:npoints,icol,1:nlev) = tautot_liq(1:npoints,icol,1:nlev)  + &
               tau_part(1:npoints,1:nlev,INDX_LSLIQ) + &
               tau_part(1:npoints,1:nlev,INDX_CVLIQ)
       endif

       ! ##############################################################################    
       ! Optical depths used by the PARASOL simulator
       ! ##############################################################################             
       if (lparasol) then
          tautot_S_liq(:,icol) = tau_part(:,nlev,1)+tau_part(:,nlev,3)
          tautot_S_ice(:,icol) = tau_part(:,nlev,2)+tau_part(:,nlev,4)              
       endif
    enddo
    
  end subroutine lidar_optics

end module cosp_optics
