MODULE tideini
  !!=================================================================================
  !!                       ***  MODULE  tideini  ***
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
  USE tide_mod
  USE iom

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::  &
       omega_tide,  &
       v0tide,      &
       utide,       &
       ftide

  LOGICAL, PUBLIC :: ln_tide_pot = .false., ln_tide_ramp = .false.
  REAL(wp), PUBLIC :: rdttideramp 
  INTEGER, PUBLIC :: nb_harmo
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: ntide
  INTEGER, PUBLIC :: kt_tide

  !!---------------------------------------------------------------------------------
  !!   OPA 9.0 , LODYC-IPSL  (2003)
  !!---------------------------------------------------------------------------------

CONTAINS
   
  SUBROUTINE tide_init ( kt )
    !!----------------------------------------------------------------------
    !!                 ***  ROUTINE tide_init  ***
    !!----------------------------------------------------------------------      
    !! * Local declarations
    INTEGER  :: ji, jk
    INTEGER, INTENT( in ) ::   kt     ! ocean time-step
    CHARACTER(LEN=4), DIMENSION(jpmax_harmo) :: clname
    !
    NAMELIST/nam_tide/ln_tide_pot, ln_tide_ramp, rdttideramp, clname
    !!----------------------------------------------------------------------

    IF ( kt == nit000 ) THEN
       !
       IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'tide_init : Initialization of the tidal components'
          WRITE(numout,*) '~~~~~~~~~ '
       ENDIF
       !
       CALL tide_init_Wave
       !
       clname(:)=''
       !
       ! Read Namelist nam_tide
       REWIND ( numnam )
       READ   ( numnam, nam_tide )
       !
       nb_harmo=0
       DO jk=1,jpmax_harmo
          DO ji=1,jpmax_harmo
             IF(TRIM(clname(jk)) .eq. Wave(ji)%cname_tide) THEN
                nb_harmo=nb_harmo+1
             ENDIF
          END DO
       ENDDO
       !
       IF(lwp) THEN
          WRITE(numout,*) '        Namelist nam_tide'
          WRITE(numout,*) '        nb_harmo    = ', nb_harmo
          WRITE(numout,*) '        ln_tide_ramp = ', ln_tide_ramp 
          WRITE(numout,*) '        rdttideramp = ', rdttideramp
          IF (ln_tide_ramp.AND.((nitend-nit000+1)*rdt/rday < rdttideramp)) &
          & CALL ctl_stop('rdttideramp must be lower than run duration')
          IF (ln_tide_ramp.AND.(rdttideramp<0.)) &
          & CALL ctl_stop('rdttideramp must be positive')
          CALL flush(numout)
       ENDIF
       !
       ALLOCATE(ntide(nb_harmo))
       DO jk=1,nb_harmo
          DO ji=1,jpmax_harmo
             IF (TRIM(clname(jk)) .eq. Wave(ji)%cname_tide) THEN
                ntide(jk) = ji
                EXIT
             END IF
          END DO
       END DO
       !
       ALLOCATE(omega_tide(nb_harmo))
       ALLOCATE(v0tide    (nb_harmo))
       ALLOCATE(utide     (nb_harmo))
       ALLOCATE(ftide     (nb_harmo))
       kt_tide = kt
       !
    ENDIF

    IF ( nsec_day == NINT(0.5 * rdttra(1)) ) THEN
       !
       IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'tide_ini : Update of the tidal components at kt=',kt
          WRITE(numout,*) '~~~~~~~~ '
       ENDIF
       CALL tide_harmo(omega_tide, v0tide, utide, ftide, ntide, nb_harmo)
       DO jk =1,nb_harmo
         IF(lwp) WRITE(numout,*) Wave(ntide(jk))%cname_tide,utide(jk),ftide(jk),v0tide(jk),omega_tide(jk)
         call flush(numout)
       END DO
       !
       kt_tide = kt
       !
    ENDIF

  END SUBROUTINE tide_init
   
END MODULE tideini
