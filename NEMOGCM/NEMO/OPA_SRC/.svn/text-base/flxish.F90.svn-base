MODULE flxish
  !!==========================================================================
  !!                 *** MODULE flxish  ***
  !! add net heat and fresh water flux from ice shelf melting
  !! into the adjacent ocean using the parameterisation by
  !! Beckmann and Goosse (2003), "A parameterization of ice shelf-ocean
  !!     interaction for climate models", Ocean Modelling 5(2003) 157-170.
  !!  (hereafter BG)
  !!==========================================================================
#if defined key_iceshelf || defined key_esopa
   !!----------------------------------------------------------------------
   !!              key_iceshelf                             ICE SHELF PARAM
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   flx_ish      : routine called from step
   !!----------------------------------------------------------------------
  !! * Modules used
  USE iom
  USE oce            ! ocean parameters
  USE dom_oce        ! ocean space and time domain
  USE phycst         ! physcial parameters
  USE par_kind       ! kind parameters
  USE par_oce        ! set the ocean parameters
  USE lbclnk         ! lateral bondary condition
  USE in_out_manager !
  USE wrk_nemo       ! for dynamical allocation
  USE lib_mpp        ! distributed memory computing library
  USE timing         ! Timing
   
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC flx_ish     !      Called by step.F90

  INTEGER, PARAMETER  :: jpish = 1000  ! max number of ice-shelf points
  INTEGER , PUBLIC    :: nishg         ! Total number of gridpoints for ice shelf 

  REAL(wp)   ::      zrhois   = 920._wp        !: density of the ice shelf
  REAL(wp)   ::      zlis     = 3.34e+5_wp     !: latent heat of fusion of the ice shelf
  REAL(wp)   ::      zgamat   = 1.e-4_wp       !: ice shelf-ocean heat exchange coefficient
  
  ! public in order to be able to output then 
  REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:,:) ::  qish            !: net heat flux from ice shelf
  REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:,:) ::  qfmass          !: fresh water flux 
  REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)   ::  qish2d          
  REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION (:,:)   ::  qfmass2d

  INTEGER,  ALLOCATABLE, SAVE, DIMENSION(:)   :: iishf         !: ice shelf location index in i direction
  INTEGER,  ALLOCATABLE, SAVE, DIMENSION(:)   :: ijshf         !:     ...  location index in j direction
  INTEGER,  ALLOCATABLE, SAVE, DIMENSION(:)   :: ikshf         !:     ...      base index in k direction
  INTEGER,  ALLOCATABLE, SAVE, DIMENSION(:)   :: ishsca        !:     ...  Leff need for parametrization
  REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: rmskisdta   ! mask iceshelf
  REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: rLefisdta   ! mask iceshelf

  REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   :: shfdepth     !:depth of ice shelf base
  REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:) :: sume3t       !: sum of the scale factor in depth between 200m and 600m

  LOGICAL, PUBLIC, PARAMETER ::  lk_iceshelf  = .TRUE.   !: 

  !! * Substitutions
#  include "domzgr_substitute.h90"
  !!----------------------------
  
CONTAINS
  INTEGER FUNCTION flx_ish_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION flx_ish_alloc  ***
      !!----------------------------------------------------------------------
      flx_ish_alloc = 0       ! set to zero if no array to be allocated
      IF( .NOT. ALLOCATED( qish ) ) THEN
         ALLOCATE( qish(jpi,jpj,jpk), qfmass(jpi,jpj,jpk), qish2d(jpi,jpj), qfmass2d(jpi,jpj), &
               &    iishf(jpish), ijshf(jpish), ikshf(jpish), ishsca(jpish),                   &
               &    shfdepth(jpish), sume3t(jpi,jpj),                                          &
               &    rmskisdta(jpi,jpj), rLefisdta(jpi,jpj), STAT= flx_ish_alloc )
         !
         IF( lk_mpp             )   CALL mpp_sum ( flx_ish_alloc )
         IF( flx_ish_alloc /= 0 )   CALL ctl_warn('flx_ish_alloc: failed to allocate arrays.')
         !
      ENDIF
  END FUNCTION flx_ish_alloc

  SUBROUTINE flx_ish(kt)
    !!----------------------------------------------------------------------
    !!                       ***  ROUTINE flx_ish  ***
    !!
    !! ** Purpose :   Add heat and fresh water fluxes due to ice shelf melting
    !!
    !! History :
    !!      !  06-02  (C. Wang) Original code
    !!----------------------------------------------------------------------

    INTEGER, INTENT ( in ) :: kt

    INTEGER :: ji, jj, jk, jish  !temporary integer
    INTEGER :: ijkmin, ikbot
    INTEGER :: ii, ij, ik 
    INTEGER :: inum

    REAL(wp) ::   rmin, zhdept   ! temporary scalars
    
    REAL(wp), POINTER, DIMENSION(:,:)   :: zt_sum      ! sum of the temperature between 200m and 600m
    REAL(wp), POINTER, DIMENSION(:,:)   :: zt_ave      ! averaged temperature between 200m and 600m

    REAL(wp), POINTER, DIMENSION(:,:,:) :: zta         ! temporary variable
    REAL(wp), POINTER, DIMENSION(:,:,:) :: zt_frz      ! freezing point temperature at depth z
    
    CHARACTER(LEN=256)  :: cnameis                     ! name of iceshelf file
    CHARACTER (LEN=32)  :: cvarmask                    ! variable name for iceshelf mask
    CHARACTER (LEN=32)  :: cvarLeff                    ! variable name for efficient Length scale
    !!----------------------------------------------------------------------
    IF ( nn_timing == 1 ) CALL timing_start('flx_ish')
     !
    CALL wrk_alloc( jpi,jpj, zt_sum, zt_ave    )
    CALL wrk_alloc( jpi,jpj, jpk, zta, zt_frz  )

    !! Initialize arrays to 0 (each step)
    zt_frz(:,:,:) = 0.e0_wp
    zt_ave(:,:)   = 0.e0_wp
    
    zt_sum(:,:)   = 0.e0_wp
    zta(:,:,:)    = 0.e0_wp

    ! Read Data and save some integral values
    IF( kt == nit000) THEN
       IF ( lwp ) WRITE(numout,*) 'flxish : reading data ... '
       IF ( lwp ) WRITE(numout,*) '~~~~~~ '
       !
       IF ( flx_ish_alloc() /= 0 )  CALL ctl_stop( 'STOP', 'flx_ish : unable to allocate arrays' )

       cnameis = 'mask_iceshelf.nc'     ! ice shelf mask file
       cvarmask  = 'zish'               ! variable name for ice shelf mask
       cvarLeff  = 'Leff'               ! variable name for Efficient Length scale

       CALL iom_open( cnameis, inum )
       CALL iom_get( inum, jpdom_data, cvarmask, rmskisdta, 1)
       CALL iom_get( inum, jpdom_data, cvarLeff, rLefisdta, 1)
       CALL iom_close(inum)

       nishg = 0
       rmin  = 999999999.
       sume3t(:,:) = 0.e0_wp
       ikshf = 99999

       DO ji=1,jpi
          DO jj=1,jpj
             ! IS is detected when next test is true
             ikbot  = mbkt(ji,jj)
             zhdept = fsdept( ji,jj, ikbot)
             IF ( ( rmskisdta(ji,jj) > 50. ) .AND. ( zhdept >  rmskisdta(ji,jj) ) ) THEN

                nishg=nishg+1

                iishf(nishg) = mig(ji)    ! global position of ice shelf point I
                ijshf(nishg) = mjg(jj)    ! global position of ice shelf point J
                shfdepth(nishg) = rmskisdta(ji,jj)  ! 
                ishsca(nishg)  = rLefisdta(ji,jj)*1000

                DO jk=1,jpk
                   IF ((fsdept(ji,jj,jk) <=  600.) .AND. (fsdept(ji,jj,jk) >  200.)) THEN
                      sume3t(ji,jj) = sume3t(ji,jj)+fse3t(ji,jj,jk)*tmask(ji,jj,jk) ! total thickness use in mean value of T
                   END IF

                   IF (ABS(fsdept(ji,jj,jk)-shfdepth(nishg)) <=  rmin) THEN
                      rmin = ABS(fsdept(ji,jj,jk)-shfdepth(nishg))
                      IF ( fsdept(ji,jj,jk) <= shfdepth(nishg) ) ijkmin= jk+1
                      IF ( fsdept(ji,jj,jk) >= shfdepth(nishg) ) ijkmin= jk
                   END IF
                END DO
                ikshf(nishg) = ijkmin

                IF (shfdepth(nishg) > 1000) THEN
                   PRINT *, 'PB ds le mask ice shelf'
                   nstop = nstop + 1                   
                END IF
                IF ((ikshf(nishg) > 40 .OR. ikshf(nishg) < 4) .AND. (tmask(ji,jj,ikshf(nishg)) ==1 ) ) THEN
                   PRINT *, 'PB ds la recherche du niveau de la base de l isf'
                   nstop = nstop + 1
                END IF
                IF (ikshf(nishg) /=  0) PRINT *, 'Niveau concerne pour l iceshelf', ikshf(nishg),'ji',mig(ji),'jj',mjg(jj)
                qish2d(ji,jj)   = -99999999999999999.
                qfmass2d(ji,jj) = -9999999999999999.
             END IF
          END DO
       END DO

       IF ( lwp ) WRITE(numout,*) 'flxish : reading data done '
       IF ( lwp ) WRITE(numout,*) '         Number of ice-shelf grid points ', nishg
       IF ( lwp ) WRITE(numout,*) '~~~~~~ '

       IF( lwp ) WRITE(numout,*)
       IF( lwp ) WRITE(numout,*) 'flx_ish: heat flux of the ice shelf'
       IF( lwp ) WRITE(numout,*) '~~~~~~~~~'
    ENDIF

    
!   qish(:,:,:) = 0.e0_wp
!   qfmass(:,:,:) = 0.e0_wp

!   qish2d(:,:) = 0.e0_wp
!   qfmass2d(:,:) = 0.e0_wp

    ! This test is false only in the very first time step of a run (JMM ???- Initialy build to skip 1rst year of run )
    ! IF (kt >= INT(86400*365/rdt)  ) THEN        !  for the first year, Ice shelf melting not added 
    IF (kt >=  0) THEN       

    DO jish = 1, nishg
       ii = iishf(jish)
       ij = ijshf(jish)
       ik = ikshf(jish)
       DO jj = mj0(ij), mj1(ij)
          DO ji = mi0(ii), mi1(ii)
             !IF (sume3t(ji,jj) .NE. 0) THEN
                ! freezing point temperature  at ice shelf base BG eq. 2 (JMM sign pb ??? +7.64e-4 !!!)
                ! after verif with UNESCO, wrong sign in BG eq. 2
                zt_frz(ji,jj,ik) = 0.0939-0.057 * tsb(ji,jj,ik,jp_sal) - 7.64e-4 * shfdepth(jish)
    ! 3. -----------the average temperature between 200m and 600m ---------------------
             DO jk = 1, jpk
                IF ((fsdept(ji,jj,jk) <= 600.) .AND. (fsdept(ji,jj,jk) >= 200.)) THEN
                   zt_sum(ji,jj) = zt_sum(ji,jj) + tsb(ji,jj,jk,jp_tem) * fse3t(ji,jj,jk) * tmask(ji,jj,jk)  ! sum temp
                ENDIF               
             ENDDO
             zt_ave(ji,jj) = zt_sum(ji,jj)/sume3t(ji,jj) ! calcul mean value
    
    ! 4. ------------Net heat flux and fresh water flux due to the ice shelf
                !IF ((fsdept(ji,jj,ikshf(jish)) .LE. 600.) .AND. (fsdept(ji,jj,ikshf(jish)) .GE. 200.)) THEN
             zta(ji,jj,ik) = zgamat * ishsca(jish) * (zt_ave(ji,jj)-zt_frz(ji,jj,ik)) * tmask(ji,jj,ik)
             ! For those corresponding to zonal boundary    
             IF ( kt==nit000 )  PRINT *, zt_ave(ji,jj), zta(ji,jj,ik), tmask(ji,jj,ik), zt_frz(ji,jj,ik)
             qish(ji,jj,ik) = rau0 * rcp * e1t(ji,jj) * zta(ji,jj,ik) !heat flux

             tsa(ji,jj,ik,jp_tem) = tsa(ji,jj,ik,jp_tem) -     & 
                      &                  zta(ji,jj,ik)/fse3t(ji,jj,ik)/e2t(ji,jj) !add to temperature trend
             qish2d(ji,jj)= qish(ji,jj,ik) ! for output
             
             !reflechir sur ice shelf vertical et horizontal (avec e1 et e2) ....                  
             qfmass(ji,jj,ik) = qish(ji,jj,ik) / zrhois / zlis            !fresh water flux m3/s                  
             tsa(ji,jj,ik,jp_sal) = tsa(ji,jj,ik,jp_sal) -     &
                      &                  qfmass(ji,jj,ik)*tsb(ji,jj,ik,jp_sal)/fse3t(ji,jj,ik)/e1t(ji,jj)/e2t(ji,jj) ! salinity trend
             qfmass2d(ji,jj) = qfmass(ji,jj,ik) ! for output
             !add to salinity trend
             !ENDIF
             !ENDIF
             !ENDDO
             ENDDO
          ENDDO
       ENDDO

       CALL lbc_lnk(qish(:,:,:), 'T', 1.)
       CALL lbc_lnk(qfmass(:,:,:), 'T', 1.)
    END IF

    IF (MOD(kt,50)==0) THEN
       IF ( lwp ) WRITE(numout,*) 'flxish : '
       IF ( lwp ) WRITE(numout,*) '~~~~~~ '
       IF ( lwp ) WRITE(numout,*) ' max abs(qish)   = ', MAXVAL(ABS(qish2d))
       IF ( lwp ) WRITE(numout,*) ' max abs(qfmass) = ', MAXVAL(ABS(qish2d))
    END IF
    CALL wrk_dealloc( jpi,jpj, zt_sum, zt_ave     )
    CALL wrk_dealloc( jpi,jpj, jpk, zta, zt_frz   )
    !
    IF( nn_timing == 1 )  CALL timing_stop('flx_ish')

  END SUBROUTINE flx_ish
#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                                        NO ICE SHELF
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_iceshelf  = .FALSE.  
CONTAINS
   SUBROUTINE flx_ish( kt )          ! Empty routine
      WRITE(*,*) 'flx_ish: You should not have seen this print! error?', kt
   END SUBROUTINE flx_ish
#endif

  !!======================================================================
END MODULE flxish
