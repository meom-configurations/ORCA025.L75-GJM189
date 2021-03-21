MODULE dynkeg
   !!======================================================================
   !!                       ***  MODULE  dynkeg  ***
   !! Ocean dynamics:  kinetic energy gradient trend
   !!======================================================================
   !! History :  1.0  !  87-09  (P. Andrich, m.-a. Foujols)  Original code
   !!            7.0  !  97-05  (G. Madec)  Split dynber into dynkeg and dynhpg
   !!            9.0  !  02-07  (G. Madec)  F90: Free form and module
   !!----------------------------------------------------------------------
   
   !!----------------------------------------------------------------------
   !!   dyn_keg      : update the momentum trend with the horizontal tke
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE prtctl          ! Print control
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_keg    ! routine called by step module
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynkeg.F90 3231 2011-12-21 09:11:11Z smasson $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_keg( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_keg  ***
      !!
      !! ** Purpose :   Compute the now momentum trend due to the horizontal
      !!      gradient of the horizontal kinetic energy and add it to the 
      !!      general momentum trend.
      !!
      !! ** Method  :   Compute the now horizontal kinetic energy 
      !!         zhke = 1/2 [ mi-1( un^2 ) + mj-1( vn^2 ) ]
      !!      Take its horizontal gradient and add it to the general momentum
      !!      trend (ua,va).
      !!         ua = ua - 1/e1u di[ zhke ]
      !!         va = va - 1/e2v dj[ zhke ]
      !!
      !! ** Action : - Update the (ua, va) with the hor. ke gradient trend
      !!             - save this trends (l_trddyn=T) for post-processing
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   iim1, iip1, ijm1, ijp1
      REAL(wp) ::   zu, zv       ! temporary scalars
      REAL(wp) ::   zcoef, zcoef2, zcoef3        ! temporary scalars
      REAL(wp) ::   zuij, zuijm1, zuijp1,  zuibarj2, zuim1j, zuim1jm1, zuim1jp1, zuim1barj2
      REAL(wp) ::   zvij, zvim1j, zvip1j,  zvjbari2, zvijm1, zvim1jm1, zvip1jm1, zvjm1bari2
!JM     REAL(wp), POINTER, DIMENSION(:,:,:) :: zhke
      REAL(wp), POINTER, DIMENSION(:,:) :: zhke
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztrdu, ztrdv 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_keg')
      !
!JM      CALL wrk_alloc( jpi, jpj, jpk, zhke )
      CALL wrk_alloc( jpi, jpj, zhke )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_keg : kinetic energy gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      IF( l_trddyn ) THEN           ! Save ua and va trends
         CALL wrk_alloc( jpi,jpj,jpk, ztrdu, ztrdv )
         ztrdu(:,:,:) = ua(:,:,:) 
         ztrdv(:,:,:) = va(:,:,:) 
      ENDIF
      
      zcoef=0.5 *1./6_wp  ! 1/12
      zcoef2 = 1./6_wp   ! 1/6
      zcoef3 = 1./48_wp
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 2, jpjm1         ! Horizontal kinetic energy at T-point
            ijm1=jj-1
            ijp1=jj+1
            DO ji = fs_2, jpim1   ! vector opt.
               iim1=ji-1
               iip1=ji+1
! Hollingsworth correction 
!following A2000
!              zu = zcoef * &
!                  & ( ( 2 * un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk) &
!                  &    + 0.25 * (un(ji-1,jj-1  ,jk) + un(ji-1,jj+1  ,jk)) &
!                  &           * (un(ji-1,jj-1  ,jk) + un(ji-1,jj+1  ,jk)) ) &
!                  &  + ( 2 * un(ji,jj  ,jk) * un(ji,jj  ,jk) &
!                  &    + 0.25 * (un(ji,jj-1  ,jk) + un(ji,jj+1  ,jk)) &
!                  &           * (un(ji,jj-1  ,jk) + un(ji,jj+1  ,jk)) ) )
!              zv =  zcoef * &
!                  & ( ( 2 * vn(ji,jj-1  ,jk) * vn(ji,jj-1  ,jk) &
!                  &    + 0.25 * (vn(ji-1,jj-1  ,jk) + vn(ji+1,jj-1  ,jk)) &
!                  &           * (vn(ji-1,jj-1  ,jk) + vn(ji+1,jj-1  ,jk)) ) &
!                  &  + ( 2 * vn(ji,jj  ,jk) * vn(ji,jj  ,jk) &
!                  &    + 0.25 * (vn(ji-1,jj  ,jk) + vn(ji+1,jj  ,jk)) &
!                  &           * (vn(ji-1,jj  ,jk) + vn(ji+1,jj  ,jk)) ) )

               zu = zcoef * &
                   & ( ( 2 * un(iim1,jj  ,jk) * un(iim1,jj  ,jk) &
                   &    + 0.25 * (un(iim1,ijm1  ,jk) + un(iim1,ijp1  ,jk)) &
                   &           * (un(iim1,ijm1  ,jk) + un(iim1,ijp1  ,jk)) ) &
                   &  + ( 2 * un(ji,jj  ,jk) * un(ji,jj  ,jk) &
                   &    + 0.25 * (un(ji,ijm1  ,jk) + un(ji,ijp1  ,jk)) &
                   &           * (un(ji,ijm1  ,jk) + un(ji,ijp1  ,jk)) ) )
               zv =  zcoef * &
                   & ( ( 2 * vn(ji,ijm1  ,jk) * vn(ji,ijm1  ,jk) &
                   &    + 0.25 * (vn(iim1,ijm1  ,jk) + vn(iip1,ijm1  ,jk)) &
                   &           * (vn(iim1,ijm1  ,jk) + vn(iip1,ijm1  ,jk)) ) &
                   &  + ( 2 * vn(ji,jj  ,jk) * vn(ji,jj  ,jk) &
                   &    + 0.25 * (vn(iim1,jj  ,jk) + vn(iip1,jj  ,jk)) &
                   &           * (vn(iim1,jj  ,jk) + vn(iip1,jj  ,jk)) ) )


! JMM try to optimize because dynkeg is now very expensive !
! 
!              zuim1j     = un(iim1, jj    ,jk)
!              zuim1jm1   = un(iim1, ijm1  ,jk)
!              zuim1jp1   = un(iim1, ijp1  ,jk)
!              zuim1barj2 = zuim1jm1 + zuim1jp1

!              zuij       = un(ji,   jj    ,jk)
!              zuijm1     = un(ji,   ijm1  ,jk)
!              zuijp1     = un(ji,   ijp1  ,jk)
!              zuibarj2   = zuijm1 + zuijp1

!              zvijm1     = vn(ji,   ijm1  ,jk)
!              zvim1jm1   = vn(iim1, ijm1  ,jk)
!              zvip1jm1   = vn(iip1, ijm1  ,jk)
!              zvjm1bari2 = zvim1jm1 + zvip1jm1

!              zvij       = vn(ji,   jj    ,jk)
!              zvim1j     = vn(iim1, jj    ,jk)
!              zvip1j     = vn(iip1, jj    ,jk)
!              zvjbari2   = zvim1j + zvip1j

!              zu = zcoef2 * ( zuim1j*zuim1j + zuij  *zuij) + zcoef3 *( zuim1barj2 * zuim1barj2 + zuibarj2 * zuibarj2  )
!              zv = zcoef2 * ( zvijm1*zvijm1 + zvij  *zvij) + zcoef3 *( zvjm1bari2 * zvjm1bari2 + zvjbari2 * zvjbari2  )

!JM               zhke(ji,jj,jk) = zv + zu
               zhke(ji,jj) = zv + zu
            END DO  
         END DO 
!
         CALL lbc_lnk( zhke, 'T', 1. )
! 
         DO jj = 2, jpjm1       ! add the gradient of kinetic energy to the general momentum trends
            DO ji = fs_2, fs_jpim1   ! vector opt.
!JM               ua(ji,jj,jk) = ua(ji,jj,jk) - ( zhke(ji+1,jj  ,jk) - zhke(ji,jj,jk) ) / e1u(ji,jj)
!JM               va(ji,jj,jk) = va(ji,jj,jk) - ( zhke(ji  ,jj+1,jk) - zhke(ji,jj,jk) ) / e2v(ji,jj)
               ua(ji,jj,jk) = ua(ji,jj,jk) - ( zhke(ji+1,jj  ) - zhke(ji,jj) ) / e1u(ji,jj)
               va(ji,jj,jk) = va(ji,jj,jk) - ( zhke(ji  ,jj+1) - zhke(ji,jj) ) / e2v(ji,jj)
            END DO 
         END DO
!!gm idea to be tested  ==>>   is it faster on scalar computers ?
!         DO jj = 2, jpjm1       ! add the gradient of kinetic energy to the general momentum trends
!            DO ji = fs_2, fs_jpim1   ! vector opt.
!               ua(ji,jj,jk) = ua(ji,jj,jk) - 0.25 * ( + un(ji+1,jj  ,jk) * un(ji+1,jj  ,jk)   &
!                  &                                   + vn(ji+1,jj-1,jk) * vn(ji+1,jj-1,jk)   &
!                  &                                   + vn(ji+1,jj  ,jk) * vn(ji+1,jj  ,jk)   &
!                  !
!                  &                                   - un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)   &
!                  &                                   - vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)   &
!                  &                                   - vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)   ) / e1u(ji,jj)
!                  !
!               va(ji,jj,jk) = va(ji,jj,jk) - 0.25 * (   un(ji-1,jj+1,jk) * un(ji-1,jj+1,jk)   &
!                  &                                   + un(ji  ,jj+1,jk) * un(ji  ,jj+1,jk)   &
!                  &                                   + vn(ji  ,jj+1,jk) * vn(ji  ,jj+1,jk)   &
!                  !
!                  &                                   - un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)   &
!                  &                                   - un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)   &
!                  &                                   - vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)   ) / e2v(ji,jj)
!            END DO 
!         END DO
!!gm en idea            <<==
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      IF( l_trddyn ) THEN      ! save the Kinetic Energy trends for diagnostic
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_mod( ztrdu, ztrdv, jpdyn_trd_keg, 'DYN', kt )
         CALL wrk_dealloc( jpi,jpj,jpk, ztrdu, ztrdv )
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' keg  - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
!JM      CALL wrk_dealloc( jpi, jpj, jpk, zhke )
      CALL wrk_dealloc( jpi, jpj, zhke )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_keg')
      !
   END SUBROUTINE dyn_keg

   !!======================================================================
END MODULE dynkeg
