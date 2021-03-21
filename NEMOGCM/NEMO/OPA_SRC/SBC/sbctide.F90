MODULE sbctide
  !!=================================================================================
  !!                       ***  MODULE  sbctide  ***
  !! Initialization of tidal forcing
  !! History :  9.0  !  07  (O. Le Galloudec)  Original code
  !!=================================================================================
  !! * Modules used
  USE oce             ! ocean dynamics and tracers variables
  USE dom_oce         ! ocean space and time domain
  USE in_out_manager  ! I/O units
  USE ioipsl          ! NetCDF IPSL library
  USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
  USE phycst
  USE daymod
  USE dynspg_oce
  USE tideini
  USE iom

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: pot_astro

#if defined key_tide

  LOGICAL, PUBLIC, PARAMETER ::   lk_tide  = .TRUE.
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: amp_pot,phi_pot
  !!---------------------------------------------------------------------------------
  !!   OPA 9.0 , LODYC-IPSL  (2003)
  !!---------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE sbc_tide ( kt )
    !!----------------------------------------------------------------------
    !!                 ***  ROUTINE sbc_tide  ***
    !!----------------------------------------------------------------------      
    !! * Arguments
    INTEGER, INTENT( in ) ::   kt     ! ocean time-step
    !!----------------------------------------------------------------------

    IF ( kt == nit000 .AND. .NOT. lk_dynspg_ts )  CALL ctl_stop( 'STOP', 'sbc_tide : tidal potential use only with time splitting' )

    IF ( nsec_day == NINT(0.5 * rdttra(1)) ) THEN
      !
      kt_tide = kt

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_tide : (re)Initialization of the tidal potential at kt=',kt
         WRITE(numout,*) '~~~~~~~ '
      ENDIF

      IF(lwp) THEN
         IF ( kt == nit000 ) WRITE(numout,*) 'Apply astronomical potential : ln_tide_pot =', ln_tide_pot
         CALL flush(numout)
      ENDIF

      IF ( kt == nit000 ) ALLOCATE(amp_pot(jpi,jpj,nb_harmo))
      IF ( kt == nit000 ) ALLOCATE(phi_pot(jpi,jpj,nb_harmo))
      IF ( kt == nit000 ) ALLOCATE(pot_astro(jpi,jpj))

      amp_pot(:,:,:) = 0.e0
      phi_pot(:,:,:) = 0.e0
      pot_astro(:,:) = 0.e0

      IF ( ln_tide_pot ) CALL tide_init_potential
      !
    ENDIF

  END SUBROUTINE sbc_tide

  SUBROUTINE tide_init_potential
    !!----------------------------------------------------------------------
    !!                 ***  ROUTINE tide_init_potential  ***
    !!----------------------------------------------------------------------
    !! * Local declarations
    INTEGER  :: ji,jj,jk
    REAL(wp) :: zcons,ztmp1,ztmp2,zlat,zlon


    DO jk=1,nb_harmo
       zcons=0.7*Wave(ntide(jk))%equitide*ftide(jk)
       do ji=1,jpi
          do jj=1,jpj
             ztmp1 = amp_pot(ji,jj,jk)*COS(phi_pot(ji,jj,jk))
             ztmp2 = -amp_pot(ji,jj,jk)*SIN(phi_pot(ji,jj,jk))
             zlat = gphit(ji,jj)*rad !! latitude en radian
             zlon = glamt(ji,jj)*rad !! longitude en radian
             ! le potentiel est composé des effets des astres:
             IF (Wave(ntide(jk))%nutide .EQ.1) THEN
                ztmp1= ztmp1 + zcons*(SIN(2.*zlat))*COS(v0tide(jk)+utide(jk)+Wave(ntide(jk))%nutide*zlon)
                ztmp2= ztmp2 - zcons*(SIN(2.*zlat))*SIN(v0tide(jk)+utide(jk)+Wave(ntide(jk))%nutide*zlon)
             ENDIF
             IF (Wave(ntide(jk))%nutide.EQ.2) THEN
                ztmp1= ztmp1 + zcons*(COS(zlat)**2)*COS(v0tide(jk)+utide(jk)+Wave(ntide(jk))%nutide*zlon)
                ztmp2= ztmp2 - zcons*(COS(zlat)**2)*SIN(v0tide(jk)+utide(jk)+Wave(ntide(jk))%nutide*zlon)
             ENDIF
             amp_pot(ji,jj,jk)=SQRT(ztmp1**2+ztmp2**2)
             phi_pot(ji,jj,jk)=ATAN2(-ztmp2/MAX(1.E-10,SQRT(ztmp1**2+ztmp2**2)),ztmp1/MAX(1.E-10,SQRT(ztmp1**2+ztmp2**2)))
          enddo
       enddo
    END DO

  END SUBROUTINE tide_init_potential

#else
  !!----------------------------------------------------------------------
  !!   Default case :   Empty module
  !!----------------------------------------------------------------------
  LOGICAL, PUBLIC, PARAMETER ::   lk_tide = .FALSE.
CONTAINS
  SUBROUTINE sbc_tide( kt )      ! Empty routine
    INTEGER         , INTENT(in) ::   kt         ! ocean time-step
    WRITE(*,*) 'sbc_tide: You should not have seen this print! error?', kt
  END SUBROUTINE sbc_tide
#endif
  !!======================================================================

END MODULE sbctide
