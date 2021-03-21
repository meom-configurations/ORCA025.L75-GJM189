MODULE limsbc_2
   !!======================================================================
   !!                       ***  MODULE limsbc_2   ***
   !! LIM-2 :   updates the fluxes at the ocean surface with ice-ocean fluxes
   !!======================================================================
   !! History :  LIM  ! 2000-01 (H. Goosse) Original code
   !!            1.0  ! 2002-07 (C. Ethe, G. Madec) re-writing F90
   !!            3.0  ! 2006-07 (G. Madec) surface module
   !!            3.3  ! 2009-05 (G. Garric, C. Bricaud) addition of the lim2_evp case
   !!             -   ! 2010-11 (G. Madec) ice-ocean stress computed at each ocean time-step
   !!           3.3.1 ! 2011-01 (A. R. Porter, STFC Daresbury) dynamical allocation
   !!            3.5  ! 2012-11 ((G. Madec, Y. Aksenov, A. Coward) salt and heat fluxes associated with e-p
   !!----------------------------------------------------------------------
#if defined key_lim2
   !!----------------------------------------------------------------------
   !!   'key_lim2'                                    LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_sbc_alloc_2 : allocate the limsbc arrays
   !!   lim_sbc_init    : initialisation
   !!   lim_sbc_flx_2   : update mass, heat and salt fluxes at the ocean surface
   !!   lim_sbc_tau_2   : update i- and j-stresses, and its modulus at the ocean surface
   !!----------------------------------------------------------------------
   USE par_oce          ! ocean parameters
   USE phycst           ! physical constants
   USE dom_oce          ! ocean domain
   USE dom_ice_2        ! LIM-2: ice domain
   USE ice_2            ! LIM-2: ice variables
   USE sbc_ice          ! surface boundary condition: ice
   USE sbc_oce          ! surface boundary condition: ocean
   USE sbccpl
   USE cpl_oasis3, ONLY : lk_cpl
   USE oce       , ONLY : sshn, sshb, snwice_mass, snwice_mass_b, snwice_fmass 
   USE albedo           ! albedo parameters
   USE lbclnk           ! ocean lateral boundary condition - MPP exchanges
   USE lib_mpp          ! MPP library
   USE wrk_nemo         ! work arrays
   USE in_out_manager   ! I/O manager
   USE diaar5, ONLY :   lk_diaar5
   USE iom              ! I/O library
   USE prtctl           ! Print control
   USE lib_fortran      ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_sbc_init_2     ! called by ice_init_2
   PUBLIC   lim_sbc_flx_2      ! called by sbc_ice_lim_2
   PUBLIC   lim_sbc_tau_2      ! called by sbc_ice_lim_2

   REAL(wp)  ::   r1_rdtice            ! = 1. / rdt_ice 
   REAL(wp)  ::   epsi16 = 1.e-16_wp   ! constant values
   REAL(wp)  ::   rzero  = 0._wp       !     -      -
   REAL(wp)  ::   rone   = 1._wp       !     -      -
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   soce_0, sice_0   ! constant SSS and ice salinity used in levitating sea-ice case
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   utau_oce, vtau_oce   ! air-ocean surface i- & j-stress              [N/m2]
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tmod_io              ! modulus of the ice-ocean relative velocity   [m/s]

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/LIM2 4.0 , UCL - NEMO Consortium (2011)
   !! $Id: limsbc_2.F90 3905 2013-05-24 13:40:39Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION lim_sbc_alloc_2()
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE lim_sbc_alloc_2 ***
      !!-------------------------------------------------------------------
      ALLOCATE( soce_0(jpi,jpj) , utau_oce(jpi,jpj) ,                       &
         &      sice_0(jpi,jpj) , vtau_oce(jpi,jpj) , tmod_io(jpi,jpj), STAT=lim_sbc_alloc_2)
         !
      IF( lk_mpp               )   CALL mpp_sum( lim_sbc_alloc_2 )
      IF( lim_sbc_alloc_2 /= 0 )   CALL ctl_warn('lim_sbc_alloc_2: failed to allocate arrays.')
      !
   END FUNCTION lim_sbc_alloc_2


   SUBROUTINE lim_sbc_flx_2( kt )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE lim_sbc_2 ***
      !!  
      !! ** Purpose :   Update surface ocean boundary condition over areas
      !!              that are at least partially covered by sea-ice
      !!         
      !! ** Action  : - comput. of the momentum, heat and freshwater/salt
      !!                fluxes at the ice-ocean interface.
      !!              - Update the fluxes provided to the ocean
      !!     
      !! ** Outputs : - qsr     : sea heat flux    : solar 
      !!              - qns     : sea heat flux    : non solar (including heat content of the mass flux)
      !!              - emp     : freshwater budget: mass flux 
      !!              - sfx     : freshwater budget: salt flux due to Freezing/Melting
      !!              - utau    : sea surface i-stress (ocean referential)
      !!              - vtau    : sea surface j-stress (ocean referential)
      !!              - fr_i    : ice fraction
      !!              - tn_ice  : sea-ice surface temperature
      !!              - alb_ice : sea-ice alberdo (lk_cpl=T)
      !!
      !! References : Goosse, H. et al. 1996, Bul. Soc. Roy. Sc. Liege, 65, 87-90.
      !!              Tartinville et al. 2001 Ocean Modelling, 3, 95-108.
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! number of iteration
      !!
      INTEGER  ::   ji, jj   ! dummy loop indices
      INTEGER  ::   ii0, ii1, ij0, ij1         ! local integers
      INTEGER  ::   ifvt, i1mfr, idfr, iflt    !   -       -
      INTEGER  ::   ial, iadv, ifral, ifrdv    !   -       -
      REAL(wp) ::   zqsr,     zqns,   zfmm     ! local scalars
      REAL(wp) ::   zinda,    zfsalt, zemp     !   -      -
      REAL(wp) ::   zemp_snw, zqhc,   zcd      !   -      -
      REAL(wp) ::   zswitch                    !   -      -
      REAL(wp), POINTER, DIMENSION(:,:)   ::   zqnsoce       ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zalb, zalbp   ! 2D/3D workspace
      !!---------------------------------------------------------------------
     
      CALL wrk_alloc( jpi, jpj, zqnsoce )
      CALL wrk_alloc( jpi, jpj, 1, zalb, zalbp )

      SELECT CASE( nn_ice_embd )                 ! levitating or embedded sea-ice option
        CASE( 0    )   ;   zswitch = 1           ! (0) standard levitating sea-ice : salt exchange only
        CASE( 1, 2 )   ;   zswitch = 0           ! (1) levitating sea-ice: salt and volume exchange but no pressure effect
                                                 ! (2) embedded sea-ice : salt and volume fluxes and pressure
      END SELECT                                 !    

      !------------------------------------------!
      !      heat flux at the ocean surface      !
      !------------------------------------------!

      zqnsoce(:,:) = qns(:,:)
      DO jj = 1, jpj
         DO ji = 1, jpi
            zinda   = 1.0   - MAX( rzero , SIGN( rone, - ( 1.0 - pfrld(ji,jj) )   ) )
            ifvt    = zinda * MAX( rzero , SIGN( rone,  - phicif(ji,jj)           ) )
            i1mfr   = 1.0   - MAX( rzero , SIGN( rone, - ( 1.0 - frld(ji,jj) )    ) )
            idfr    = 1.0   - MAX( rzero , SIGN( rone, frld(ji,jj) - pfrld(ji,jj) ) )
            iflt    = zinda  * (1 - i1mfr) * (1 - ifvt )
            ial     = ifvt   * i1mfr + ( 1 - ifvt ) * idfr
            iadv    = ( 1  - i1mfr ) * zinda
            ifral   = ( 1  - i1mfr * ( 1 - ial ) )   
            ifrdv   = ( 1  - ifral * ( 1 - ial ) ) * iadv 

!!$            attempt to explain the tricky flags set above....
!!$            zinda   = 1.0 - AINT( pfrld(ji,jj) )                   ! = 0. if ice-free ocean else 1. (after ice adv, but before ice thermo)
!!$            i1mfr   = 1.0 - AINT(  frld(ji,jj) )                   ! = 0. if ice-free ocean else 1. (after ice thermo)
!!$
!!$            IF( phicif(ji,jj) <= 0. ) THEN   ;   ifvt = zinda      ! = zinda if previous thermodynamic step overmelted the ice???
!!$            ELSE                             ;   ifvt = 0.         ! 
!!$            ENDIF
!!$
!!$            IF( frld(ji,jj) >= pfrld(ji,jj) ) THEN   ;   idfr = 0.  !   = 0. if lead fraction increases due to ice thermodynamics
!!$            ELSE                                     ;   idfr = 1.   
!!$            ENDIF
!!$
!!$            iflt    = zinda  * (1 - i1mfr) * (1 - ifvt )    !   = 1. if ice (not only snow) at previous time and ice-free ocean currently
!!$
!!$            ial     = ifvt   * i1mfr    +    ( 1 - ifvt ) * idfr
!!$                    = i1mfr if ifvt = 1 i.e. 
!!$                    = idfr  if ifvt = 0
!!$!                 snow no ice   ice         ice or nothing  lead fraction increases
!!$!                 at previous   now           at previous
!!$!                -> ice area increases  ???         -> ice area decreases ???
!!$
!!$            iadv    = ( 1  - i1mfr ) * zinda
!!$!                     pure ocean      ice at
!!$!                     at current      previous
!!$!                        -> = 1. if ice disapear between previous and current
!!$
!!$            ifral   = ( 1  - i1mfr * ( 1 - ial ) )  
!!$!                            ice at     ???
!!$!                            current         
!!$!                         -> ???
!!$
!!$            ifrdv   = ( 1  - ifral * ( 1 - ial ) ) * iadv
!!$!                                                    ice disapear
!!$
!!$

            !   computation the solar flux at ocean surface
#if defined key_coupled 
            zqsr = qsr_tot(ji,jj) + ( fstric(ji,jj) - qsr_ice(ji,jj,1) ) * ( 1.0 - pfrld(ji,jj) )
#else
            zqsr = pfrld(ji,jj) * qsr(ji,jj)  + ( 1.  - pfrld(ji,jj) ) * fstric(ji,jj)
#endif            
            !  computation the non solar heat flux at ocean surface
            zqns    =  - ( 1. - thcm(ji,jj) ) * zqsr                                              &   ! part of the solar energy used in leads
               &       + iflt    * ( fscmbq(ji,jj) + ffltbif(ji,jj) )                             &
               &       + ifral   * ( ial * qcmif(ji,jj) + (1 - ial) * qldif(ji,jj) ) * r1_rdtice  &
               &       + ifrdv   * (       qfvbq(ji,jj) +             qdtcn(ji,jj) ) * r1_rdtice 

            fsbbq(ji,jj) = ( 1.0 - ( ifvt + iflt ) ) * fscmbq(ji,jj)     ! store residual heat flux (to put into the ocean at the next time-step)
            zqhc = ( rdq_snw(ji,jj)                                     &
                 & + rdq_ice(ji,jj) * ( 1.- zswitch) ) * r1_rdtice       ! heat flux due to snow ( & ice heat content, 
            !                                                            !           if ice/ocean mass exchange active) 
            qsr  (ji,jj) = zqsr                                          ! solar heat flux 
            qns  (ji,jj) = zqns - fdtcn(ji,jj) + zqhc                    ! non solar heat flux 
            !
            !                          !------------------------------------------!
            !                          !  mass and salt flux at the ocean surface !
            !                          !------------------------------------------!
            !
            ! mass flux at the ocean-atmosphere interface (open ocean fraction = leads area)
#if defined key_coupled
            !                                                  ! coupled mode: 
            zemp = + emp_tot(ji,jj)                            &     ! net mass flux over the grid cell (ice+ocean area)
               &   - emp_ice(ji,jj) * ( 1. - pfrld(ji,jj) )          ! minus the mass flux intercepted by sea-ice
#else
            !                                                  ! forced  mode: 
            zemp = + emp(ji,jj)     *         frld(ji,jj)      &     ! mass flux over open ocean fraction 
               &   - tprecip(ji,jj) * ( 1. -  frld(ji,jj) )    &     ! liquid precip. over ice reaches directly the ocean
               &   + sprecip(ji,jj) * ( 1. - pfrld(ji,jj) )          ! snow is intercepted by sea-ice (previous frld)
#endif            
            !
            ! mass flux at the ocean/ice interface (sea ice fraction)
            zemp_snw = rdm_snw(ji,jj) * r1_rdtice                    ! snow melting = pure water that enters the ocean
            zfmm     = rdm_ice(ji,jj) * r1_rdtice                    ! Freezing minus Melting (F-M)

            fmmflx(ji,jj) = zfmm                                     ! F/M mass flux save at least for biogeochemical model

            ! salt flux at the ice/ocean interface (sea ice fraction) [PSU*kg/m2/s]
            zfsalt = - sice_0(ji,jj) * zfmm                          ! F-M salt exchange
            zcd    =   soce_0(ji,jj) * zfmm                          ! concentration/dilution term due to F-M
            !
            ! salt flux only       : add concentration dilution term in salt flux  and no  F-M term in volume flux
            ! salt and mass fluxes : non concentration dilution term in salt flux  and add F-M term in volume flux
            sfx (ji,jj) = zfsalt +                  zswitch  * zcd   ! salt flux (+ C/D if no ice/ocean mass exchange)
            emp (ji,jj) = zemp   + zemp_snw + ( 1.- zswitch) * zfmm  ! mass flux (+ F/M mass flux if ice/ocean mass exchange)
            !
         END DO
      END DO
      !                                !------------------------------------------!
      !                                !    mass of snow and ice per unit area    !
      !                                !------------------------------------------!
      IF( nn_ice_embd /= 0 ) THEN      ! embedded sea-ice (mass required)
         snwice_mass_b(:,:) = snwice_mass(:,:)                  ! save mass from the previous ice time step
         !                                                      ! new mass per unit area
         snwice_mass  (:,:) = tms(:,:) * ( rhosn * hsnif(:,:) + rhoic * hicif(:,:)  ) * ( 1.0 - frld(:,:) )
         !                                                      ! time evolution of snow+ice mass
         snwice_fmass (:,:) = ( snwice_mass(:,:) - snwice_mass_b(:,:) ) / rdt_ice
      ENDIF

      CALL iom_put( 'hflx_ice_cea', - fdtcn(:,:) )      
      CALL iom_put( 'qns_io_cea', qns(:,:) - zqnsoce(:,:) * pfrld(:,:) )      
      CALL iom_put( 'qsr_io_cea', fstric(:,:) * (1.e0 - pfrld(:,:)) )

      IF( lk_diaar5 ) THEN       ! AR5 diagnostics
         CALL iom_put( 'isnwmlt_cea'  ,                 rdm_snw(:,:) * r1_rdtice )
         CALL iom_put( 'fsal_virt_cea',   soce_0(:,:) * rdm_ice(:,:) * r1_rdtice )
         CALL iom_put( 'fsal_real_cea', - sice_0(:,:) * rdm_ice(:,:) * r1_rdtice )
      ENDIF

      !-----------------------------------------------!
      !   Coupling variables                          !
      !-----------------------------------------------!

#if defined key_coupled
      tn_ice(:,:,1) = sist(:,:)          ! sea-ice surface temperature       
      ht_i(:,:,1) = hicif(:,:)
      ht_s(:,:,1) = hsnif(:,:)
      a_i(:,:,1) = fr_i(:,:)
      !                                  ! Computation of snow/ice and ocean albedo
      CALL albedo_ice( tn_ice, ht_i, ht_s, zalbp, zalb )
      alb_ice(:,:,1) =  0.5 * ( zalbp(:,:,1) + zalb (:,:,1) )   ! Ice albedo (mean clear and overcast skys)
      CALL iom_put( "icealb_cea", alb_ice(:,:,1) * fr_i(:,:) )  ! ice albedo
#endif

      IF(ln_ctl) THEN            ! control print
         CALL prt_ctl(tab2d_1=qsr   , clinfo1=' lim_sbc: qsr    : ', tab2d_2=qns   , clinfo2=' qns     : ')
         CALL prt_ctl(tab2d_1=emp   , clinfo1=' lim_sbc: emp    : ', tab2d_2=sfx   , clinfo2=' sfx     : ')
         CALL prt_ctl(tab2d_1=utau  , clinfo1=' lim_sbc: utau   : ', mask1=umask,   &
            &         tab2d_2=vtau  , clinfo2=' vtau    : '        , mask2=vmask )
         CALL prt_ctl(tab2d_1=fr_i  , clinfo1=' lim_sbc: fr_i   : ', tab2d_2=tn_ice(:,:,1), clinfo2=' tn_ice  : ')
      ENDIF 
      !
      CALL wrk_dealloc( jpi, jpj, zqnsoce )
      CALL wrk_dealloc( jpi, jpj, 1, zalb, zalbp )
      !
   END SUBROUTINE lim_sbc_flx_2


   SUBROUTINE lim_sbc_tau_2( kt , pu_oce, pv_oce )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE lim_sbc_tau_2 ***
      !!  
      !! ** Purpose : Update the ocean surface stresses due to the ice
      !!         
      !! ** Action  : * at each ice time step (every nn_fsbc time step):
      !!                - compute the modulus of ice-ocean relative velocity 
      !!                  at T-point (C-grid) or I-point (B-grid)
      !!                      tmod_io = rhoco * | U_ice-U_oce |
      !!                - update the modulus of stress at ocean surface
      !!                      taum = frld * taum + (1-frld) * tmod_io * | U_ice-U_oce |
      !!              * at each ocean time step (each kt): 
      !!                  compute linearized ice-ocean stresses as
      !!                      Utau = tmod_io * | U_ice - pU_oce |
      !!                using instantaneous current ocean velocity (usually before)
      !!
      !!    NB: - the averaging operator used depends on the ice dynamics grid (cp_ice_msh='I' or 'C')
      !!        - ice-ocean rotation angle only allowed in cp_ice_msh='I' case
      !!        - here we make an approximation: taum is only computed every ice time step
      !!          This avoids mutiple average to pass from T -> U,V grids and next from U,V grids 
      !!          to T grid. taum is used in TKE and GLS, which should not be too sensitive to this approximation...
      !!
      !! ** Outputs : - utau, vtau : surface ocean i- and j-stress (u- & v-pts) updated with ice-ocean fluxes
      !!              - taum       : modulus of the surface ocean stress (T-point) updated with ice-ocean fluxes
      !!---------------------------------------------------------------------
      INTEGER ,                     INTENT(in) ::   kt               ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pu_oce, pv_oce   ! surface ocean currents
      !!
      INTEGER  ::   ji, jj   ! dummy loop indices
      REAL(wp) ::   zfrldu, zat_u, zu_i, zutau_ice, zu_t, zmodt   ! local scalar
      REAL(wp) ::   zfrldv, zat_v, zv_i, zvtau_ice, zv_t, zmodi   !   -      -
      REAL(wp) ::   zsang, zumt                                   !    -         -
      REAL(wp), POINTER, DIMENSION(:,:) ::   ztio_u, ztio_v   ! ocean stress below sea-ice
      !!---------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi, jpj, ztio_u, ztio_v )
      !
      SELECT CASE( cp_ice_msh )     
      !                             !-----------------------!
      CASE( 'I' )                   !  B-grid ice dynamics  !   I-point (i.e. F-point with sea-ice indexation)
         !                          !--=--------------------!
         !
         IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN     !==  Ice time-step only  ==! (i.e. surface module time-step)
!CDIR NOVERRCHK
            DO jj = 1, jpj                               !* modulus of ice-ocean relative velocity at I-point
!CDIR NOVERRCHK
               DO ji = 1, jpi
                  zu_i  = u_ice(ji,jj) - u_oce(ji,jj)                   ! ice-ocean relative velocity at I-point
                  zv_i  = v_ice(ji,jj) - v_oce(ji,jj)
                  tmod_io(ji,jj) = SQRT( zu_i * zu_i + zv_i * zv_i )    ! modulus of this velocity (at I-point)
               END DO
            END DO
!CDIR NOVERRCHK
            DO jj = 1, jpjm1                             !* update the modulus of stress at ocean surface (T-point)
!CDIR NOVERRCHK
               DO ji = 1, jpim1   ! NO vector opt.
                  !                                               ! modulus of U_ice-U_oce at T-point
                  zumt  = 0.25_wp * (  tmod_io(ji+1,jj) + tmod_io(ji  ,jj  )    &   
                     &               + tmod_io(ji,jj+1) + tmod_io(ji+1,jj+1)  )
                  !                                               ! update the modulus of stress at ocean surface
                  taum(ji,jj) = frld(ji,jj) * taum(ji,jj) + ( 1._wp - frld(ji,jj) ) * rhoco * zumt * zumt
               END DO
            END DO
            CALL lbc_lnk( taum, 'T', 1. )
            !
            utau_oce(:,:) = utau(:,:)                    !* save the air-ocean stresses at ice time-step
            vtau_oce(:,:) = vtau(:,:)
            !
         ENDIF
         !
         !                                        !==  at each ocean time-step  ==!
         !
         !                                               !* ice/ocean stress WITH a ice-ocean rotation angle at I-point
         DO jj = 2, jpj
            zsang  = SIGN( 1._wp, gphif(1,jj) ) * sangvg          ! change the cosine angle sign in the SH 
            DO ji = 2, jpi    ! NO vect. opt. possible
               ! ... ice-ocean relative velocity at I-point using instantaneous surface ocean current at u- & v-pts
               zu_i = u_ice(ji,jj) - 0.5_wp * ( pu_oce(ji-1,jj  ) + pu_oce(ji-1,jj-1) ) * tmu(ji,jj)
               zv_i = v_ice(ji,jj) - 0.5_wp * ( pv_oce(ji  ,jj-1) + pv_oce(ji-1,jj-1) ) * tmu(ji,jj)
               ! ... components of stress with a ice-ocean rotation angle 
               zmodi = rhoco * tmod_io(ji,jj)                     
               ztio_u(ji,jj) = zmodi * ( cangvg * zu_i - zsang * zv_i )
               ztio_v(ji,jj) = zmodi * ( cangvg * zv_i + zsang * zu_i )
            END DO
         END DO
         !                                               !* surface ocean stresses at u- and v-points
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! NO vector opt.
               !                                   ! ice-ocean stress at U and V-points  (from I-point values)
               zutau_ice  = 0.5_wp * ( ztio_u(ji+1,jj) + ztio_u(ji+1,jj+1) )
               zvtau_ice  = 0.5_wp * ( ztio_v(ji,jj+1) + ztio_v(ji+1,jj+1) )
               !                                   ! open-ocean (lead) fraction at U- & V-points (from T-point values)
               zfrldu = 0.5_wp * ( frld(ji,jj) + frld(ji+1,jj) )
               zfrldv = 0.5_wp * ( frld(ji,jj) + frld(ji,jj+1) )
               !                                   ! update the surface ocean stress (ice-cover wheighted)
               utau(ji,jj) = zfrldu * utau_oce(ji,jj) + ( 1._wp - zfrldu ) * zutau_ice
               vtau(ji,jj) = zfrldv * vtau_oce(ji,jj) + ( 1._wp - zfrldv ) * zvtau_ice
            END DO
         END DO
         CALL lbc_lnk( utau, 'U', -1. )   ;   CALL lbc_lnk( vtau, 'V', -1. )     ! lateral boundary condition
         !
         !
         !                          !-----------------------!
      CASE( 'C' )                   !  C-grid ice dynamics  !   U & V-points (same as in the ocean)
         !                          !--=--------------------!
         !
         IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN     !==  Ice time-step only  ==! (i.e. surface module time-step)
!CDIR NOVERRCHK
            DO jj = 2, jpjm1                          !* modulus of the ice-ocean velocity at T-point
!CDIR NOVERRCHK
               DO ji = fs_2, fs_jpim1
                  zu_t  = u_ice(ji,jj) + u_ice(ji-1,jj) - u_oce(ji,jj) - u_oce(ji-1,jj)   ! 2*(U_ice-U_oce) at T-point
                  zv_t  = v_ice(ji,jj) + v_ice(ji,jj-1) - v_oce(ji,jj) - v_oce(ji,jj-1)      
                  zmodt =  0.25_wp * (  zu_t * zu_t + zv_t * zv_t  )                      ! |U_ice-U_oce|^2
                  !                                               ! update the modulus of stress at ocean surface
                  taum   (ji,jj) = frld(ji,jj) * taum(ji,jj) + ( 1._wp - frld(ji,jj) ) * rhoco * zmodt
                  tmod_io(ji,jj) = SQRT( zmodt ) * rhoco          ! rhoco*|Uice-Uoce|
               END DO
            END DO
            CALL lbc_lnk( taum, 'T', 1. )   ;   CALL lbc_lnk( tmod_io, 'T', 1. )
            !
            utau_oce(:,:) = utau(:,:)                 !* save the air-ocean stresses at ice time-step
            vtau_oce(:,:) = vtau(:,:)
            !
         ENDIF
         !
         !                                        !==  at each ocean time-step  ==!
         !
         DO jj = 2, jpjm1                             !* ice stress over ocean WITHOUT a ice-ocean rotation angle
            DO ji = fs_2, fs_jpim1
               !                                            ! ocean area at u- & v-points
               zfrldu  = 0.5_wp * ( frld(ji,jj) + frld(ji+1,jj) )
               zfrldv  = 0.5_wp * ( frld(ji,jj) + frld(ji,jj+1) )
               !                                            ! quadratic drag formulation without rotation
               !                                            ! using instantaneous surface ocean current
               zutau_ice = 0.5 * ( tmod_io(ji,jj) + tmod_io(ji+1,jj) ) * ( u_ice(ji,jj) - pu_oce(ji,jj) )
               zvtau_ice = 0.5 * ( tmod_io(ji,jj) + tmod_io(ji,jj+1) ) * ( v_ice(ji,jj) - pv_oce(ji,jj) )
               !                                            ! update the surface ocean stress (ice-cover wheighted)
               utau(ji,jj) = zfrldu * utau_oce(ji,jj) + ( 1._wp - zfrldu ) * zutau_ice
               vtau(ji,jj) = zfrldv * vtau_oce(ji,jj) + ( 1._wp - zfrldv ) * zvtau_ice
            END DO
         END DO
         CALL lbc_lnk( utau, 'U', -1. )   ;   CALL lbc_lnk( vtau, 'V', -1. )   ! lateral boundary condition
         !
      END SELECT

      IF(ln_ctl)   CALL prt_ctl( tab2d_1=utau, clinfo1=' lim_sbc: utau   : ', mask1=umask,   &
         &                       tab2d_2=vtau, clinfo2=' vtau    : '        , mask2=vmask )
      !  
      CALL wrk_dealloc( jpi, jpj, ztio_u, ztio_v )
      !
   END SUBROUTINE lim_sbc_tau_2


   SUBROUTINE lim_sbc_init_2
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE lim_sbc_init  ***
      !!             
      !! ** Purpose : Preparation of the file ice_evolu for the output of
      !!      the temporal evolution of key variables
      !!
      !! ** input   : Namelist namicedia
      !!-------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'lim_sbc_init_2 : LIM-2 sea-ice - surface boundary condition'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   '

      !                                      ! allocate lim_sbc arrays
      IF( lim_sbc_alloc_2() /= 0 )   CALL ctl_stop( 'STOP', 'lim_sbc_flx_2 : unable to allocate arrays' )
      !
      r1_rdtice = 1._wp / rdt_ice
      !
      soce_0(:,:) = soce                     ! constant SSS and ice salinity used in levitating sea-ice case
      sice_0(:,:) = sice
      !
      IF( cp_cfg == "orca" ) THEN            ! decrease ocean & ice reference salinities in the Baltic sea 
         WHERE( 14._wp <= glamt(:,:) .AND. glamt(:,:) <= 32._wp .AND.   &
            &   54._wp <= gphit(:,:) .AND. gphit(:,:) <= 66._wp         ) 
            soce_0(:,:) = 4._wp
            sice_0(:,:) = 2._wp
         END WHERE
      ENDIF
      !                                      ! embedded sea ice
      IF( nn_ice_embd /= 0 ) THEN            ! mass exchanges between ice and ocean (case 1 or 2) set the snow+ice mass
         snwice_mass  (:,:) = tms(:,:) * ( rhosn * hsnif(:,:) + rhoic * hicif(:,:)  ) * ( 1.0 - frld(:,:) )
         snwice_mass_b(:,:) = snwice_mass(:,:)
      ELSE
         snwice_mass  (:,:) = 0.e0           ! no mass exchanges
         snwice_mass_b(:,:) = 0.e0           ! no mass exchanges
      ENDIF
      IF( nn_ice_embd == 2 .AND.          &  ! full embedment (case 2) & no restart : 
         &   .NOT.ln_rstart ) THEN           ! deplete the initial ssh below sea-ice area
         sshn(:,:) = sshn(:,:) - snwice_mass(:,:) * r1_rau0
         sshb(:,:) = sshb(:,:) - snwice_mass(:,:) * r1_rau0
      ENDIF
      !
   END SUBROUTINE lim_sbc_init_2

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module        NO LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
#endif 

   !!======================================================================
END MODULE limsbc_2
