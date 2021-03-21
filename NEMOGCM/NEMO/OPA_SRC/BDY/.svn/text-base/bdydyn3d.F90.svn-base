MODULE bdydyn3d
   !!======================================================================
   !!                       ***  MODULE  bdydyn3d  ***
   !! Unstructured Open Boundary Cond. :   Flow relaxation scheme on baroclinic velocities
   !!======================================================================
   !! History :  3.4  !  2011     (D. Storkey) new module as part of BDY rewrite 
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) Optimization of BDY communications
   !!----------------------------------------------------------------------
#if defined key_bdy 
   !!----------------------------------------------------------------------
   !!   'key_bdy' :                    Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------
   !!   bdy_dyn3d        : apply open boundary conditions to baroclinic velocities
   !!   bdy_dyn3d_frs    : apply Flow Relaxation Scheme
   !!----------------------------------------------------------------------
   USE timing          ! Timing
   USE wrk_nemo        ! Memory Allocation
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE bdy_oce         ! ocean open boundary conditions
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  !
   Use phycst

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dyn3d     ! routine called by bdy_dyn
   PUBLIC   bdy_dyn3d_dmp ! routine called by step

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdydyn.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_dyn3d( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d  ***
      !!
      !! ** Purpose : - Apply open boundary conditions for baroclinic velocities
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) :: kt     ! Main time step counter
      !!
      INTEGER               :: ib_bdy ! loop index
      !!

      DO ib_bdy=1, nb_bdy

!!$         IF ( using Orlanski radiation conditions ) THEN 
!!$            CALL bdy_rad( kt,  bdyidx(ib_bdy) )
!!$         ENDIF

         SELECT CASE( nn_dyn3d(ib_bdy) )
         CASE(jp_none)
            CYCLE
         CASE(jp_frs)
            CALL bdy_dyn3d_frs( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy )
         CASE(2)
            CALL bdy_dyn3d_spe( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy )
         CASE(3)
            CALL bdy_dyn3d_zro( idx_bdy(ib_bdy), dta_bdy(ib_bdy), kt, ib_bdy )
         CASE DEFAULT
            CALL ctl_stop( 'bdy_dyn3d : unrecognised option for open boundaries for baroclinic velocities' )
         END SELECT
      ENDDO

   END SUBROUTINE bdy_dyn3d

   SUBROUTINE bdy_dyn3d_spe( idx, dta, kt , ib_bdy )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_spe  ***
      !!
      !! ** Purpose : - Apply a specified value for baroclinic velocities
      !!                at open boundaries.
      !!
      !!----------------------------------------------------------------------
      INTEGER                     ::   kt
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      !!
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_spe')
      !
      igrd = 2                      ! Relaxation of zonal velocity
      DO jb = 1, idx%nblenrim(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            ua(ii,ij,jk) = dta%u3d(jb,jk) * umask(ii,ij,jk)
         END DO
      END DO
      !
      igrd = 3                      ! Relaxation of meridional velocity
      DO jb = 1, idx%nblenrim(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            va(ii,ij,jk) = dta%v3d(jb,jk) * vmask(ii,ij,jk)
         END DO
      END DO
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )   ;   CALL lbc_bdy_lnk( va, 'V', -1.,ib_bdy )   ! Boundary points should be updated
      !
      IF( kt .eq. nit000 ) CLOSE( unit = 102 )

      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_spe')

   END SUBROUTINE bdy_dyn3d_spe

   SUBROUTINE bdy_dyn3d_zro( idx, dta, kt, ib_bdy )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_zro  ***
      !!
      !! ** Purpose : - baroclinic velocities = 0. at open boundaries.
      !!
      !!----------------------------------------------------------------------
      INTEGER                     ::   kt
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      !!
      INTEGER  ::   ib, ik         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd, zcoef   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_zro')
      !
      igrd = 2                       ! Everything is at T-points here
      DO ib = 1, idx%nblenrim(igrd)
         ii = idx%nbi(ib,igrd)
         ij = idx%nbj(ib,igrd)
         DO ik = 1, jpkm1
            ua(ii,ij,ik) = 0._wp
         END DO
      END DO

      igrd = 3                       ! Everything is at T-points here
      DO ib = 1, idx%nblenrim(igrd)
         ii = idx%nbi(ib,igrd)
         ij = idx%nbj(ib,igrd)
         DO ik = 1, jpkm1
            va(ii,ij,ik) = 0._wp
         END DO
      END DO
      !
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )   ;   CALL lbc_bdy_lnk( va, 'V', -1.,ib_bdy )   ! Boundary points should be updated
      !
      IF( kt .eq. nit000 ) CLOSE( unit = 102 )

      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_zro')

   END SUBROUTINE bdy_dyn3d_zro

   SUBROUTINE bdy_dyn3d_frs( idx, dta, kt, ib_bdy )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_frs  ***
      !!
      !! ** Purpose : - Apply the Flow Relaxation Scheme for baroclinic velocities
      !!                at open boundaries.
      !!
      !! References :- Engedahl H., 1995: Use of the flow relaxation scheme in 
      !!               a three-dimensional baroclinic ocean model with realistic
      !!               topography. Tellus, 365-382.
      !!----------------------------------------------------------------------
      INTEGER                     ::   kt
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      !!
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_frs')
      !
      igrd = 2                      ! Relaxation of zonal velocity
      DO jb = 1, idx%nblen(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            zwgt = idx%nbw(jb,igrd)
            ua(ii,ij,jk) = ( ua(ii,ij,jk) + zwgt * ( dta%u3d(jb,jk) - ua(ii,ij,jk) ) ) * umask(ii,ij,jk)
         END DO
      END DO
      !
      igrd = 3                      ! Relaxation of meridional velocity
      DO jb = 1, idx%nblen(igrd)
         DO jk = 1, jpkm1
            ii   = idx%nbi(jb,igrd)
            ij   = idx%nbj(jb,igrd)
            zwgt = idx%nbw(jb,igrd)
            va(ii,ij,jk) = ( va(ii,ij,jk) + zwgt * ( dta%v3d(jb,jk) - va(ii,ij,jk) ) ) * vmask(ii,ij,jk)
         END DO
      END DO 
      CALL lbc_bdy_lnk( ua, 'U', -1., ib_bdy )   ;   CALL lbc_bdy_lnk( va, 'V', -1.,ib_bdy )   ! Boundary points should be updated
      !
      IF( kt .eq. nit000 ) CLOSE( unit = 102 )

      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_frs')

   END SUBROUTINE bdy_dyn3d_frs

   SUBROUTINE bdy_dyn3d_dmp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn3d_dmp  ***
      !!
      !! ** Purpose : Apply damping for baroclinic velocities at open boundaries.
      !!
      !!----------------------------------------------------------------------
      INTEGER                     ::   kt
      !!
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      INTEGER  ::  ib_bdy          ! loop index
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 ) CALL timing_start('bdy_dyn3d_dmp')
      !
      !-------------------------------------------------------
      ! Remove barotropic part from before velocity
      !-------------------------------------------------------
      CALL wrk_alloc(jpi,jpj,pu2d,pv2d) 

      pu2d(:,:) = 0.e0
      pv2d(:,:) = 0.e0

      DO jk = 1, jpkm1
#if defined key_vvl
         pu2d(:,:) = pu2d(:,:) + fse3u_b(:,:,jk)* ub(:,:,jk)   *umask(:,:,jk) 
         pv2d(:,:) = pv2d(:,:) + fse3v_b(:,:,jk)* vb(:,:,jk)   *vmask(:,:,jk)
#else
         pu2d(:,:) = pu2d(:,:) + fse3u_0(:,:,jk) * ub(:,:,jk)  * umask(:,:,jk)
         pv2d(:,:) = pv2d(:,:) + fse3v_0(:,:,jk) * vb(:,:,jk)  * vmask(:,:,jk)
#endif
      END DO

      IF( lk_vvl ) THEN
         pu2d(:,:) = pu2d(:,:) * umask(:,:,1) / ( hu_0(:,:) + sshu_b(:,:) + 1._wp - umask(:,:,1) )
         pv2d(:,:) = pv2d(:,:) * vmask(:,:,1) / ( hv_0(:,:) + sshv_b(:,:) + 1._wp - vmask(:,:,1) )
      ELSE
         pu2d(:,:) = pv2d(:,:) * hur(:,:)
         pv2d(:,:) = pu2d(:,:) * hvr(:,:)
      ENDIF

      DO ib_bdy=1, nb_bdy
         IF ( ln_dyn3d_dmp(ib_bdy).and.nn_dyn3d(ib_bdy).gt.0 ) THEN
            igrd = 2                      ! Relaxation of zonal velocity
            DO jb = 1, idx_bdy(ib_bdy)%nblen(igrd)
               ii   = idx_bdy(ib_bdy)%nbi(jb,igrd)
               ij   = idx_bdy(ib_bdy)%nbj(jb,igrd)
               zwgt = idx_bdy(ib_bdy)%nbd(jb,igrd)
               DO jk = 1, jpkm1
                  ua(ii,ij,jk) = ( ua(ii,ij,jk) + zwgt * ( dta_bdy(ib_bdy)%u3d(jb,jk) - &
                                   ub(ii,ij,jk) + pu2d(ii,ij)) ) * umask(ii,ij,jk)
               END DO
            END DO
            !
            igrd = 3                      ! Relaxation of meridional velocity
            DO jb = 1, idx_bdy(ib_bdy)%nblen(igrd)
               ii   = idx_bdy(ib_bdy)%nbi(jb,igrd)
               ij   = idx_bdy(ib_bdy)%nbj(jb,igrd)
               zwgt = idx_bdy(ib_bdy)%nbd(jb,igrd)
               DO jk = 1, jpkm1
                  va(ii,ij,jk) = ( va(ii,ij,jk) + zwgt * ( dta_bdy(ib_bdy)%v3d(jb,jk) -  &
                                   vb(ii,ij,jk) + pv2d(ii,ij)) ) * vmask(ii,ij,jk)
               END DO
            END DO
         ENDIF
      ENDDO
      !
      CALL wrk_dealloc(jpi,jpj,pu2d,pv2d) 
      !
      CALL lbc_lnk( ua, 'U', -1. )   ;   CALL lbc_lnk( va, 'V', -1. )   ! Boundary points should be updated
      !
      IF( nn_timing == 1 ) CALL timing_stop('bdy_dyn3d_dmp')

   END SUBROUTINE bdy_dyn3d_dmp

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_dyn3d( kt )      ! Empty routine
      WRITE(*,*) 'bdy_dyn3d: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dyn3d

   SUBROUTINE bdy_dyn3d_dmp( kt )      ! Empty routine
      WRITE(*,*) 'bdy_dyn3d_dmp: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dyn3d_dmp

#endif

   !!======================================================================
END MODULE bdydyn3d
