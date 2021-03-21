MODULE dynldf_bilap
   !!======================================================================
   !!                     ***  MODULE  dynldf_bilap  ***
   !! Ocean dynamics:  lateral viscosity trend
   !!======================================================================
   !! History :  OPA  ! 1990-09  (G. Madec)  Original code
   !!            4.0  ! 1993-03  (M. Guyon)  symetrical conditions (M. Guyon)
   !!            6.0  ! 1996-01  (G. Madec)  statement function for e3
   !!            8.0  ! 1997-07  (G. Madec)  lbc calls
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90: Free form and module
   !!            2.0  ! 2004-08  (C. Talandier) New trends organization
   !!            2.0  ! 2007-12  (G. Hervieux) No slip accurate version Gaelle Hervieux)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf_bilap : update the momentum trend with the lateral diffusion
   !!                   using an iso-level bilaplacian operator
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ldfdyn_oce      ! ocean dynamics: lateral physics
   USE in_out_manager  ! I/O manager
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_ldf_bilap   ! called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldfdyn_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynldf_bilap.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_ldf_bilap( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_bilap  ***
      !!
      !! ** Purpose :   Compute the before trend of the lateral momentum
      !!      diffusion and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The before horizontal momentum diffusion trend is a 
      !!      bi-harmonic operator (bilaplacian type) which separates the
      !!      divergent and rotational parts of the flow.
      !!      Its horizontal components are computed as follow:
      !!      laplacian:
      !!          zlu = 1/e1u di[ hdivb ] - 1/(e2u*e3u) dj-1[ e3f rotb ]
      !!          zlv = 1/e2v dj[ hdivb ] + 1/(e1v*e3v) di-1[ e3f rotb ]
      !!      third derivative:
      !!       * multiply by the eddy viscosity coef. at u-, v-point, resp.
      !!          zlu = ahmu * zlu
      !!          zlv = ahmv * zlv
      !!       * curl and divergence of the laplacian
      !!          zuf = 1/(e1f*e2f) ( di[e2v zlv] - dj[e1u zlu] )
      !!          zut = 1/(e1t*e2t*e3t) ( di[e2u*e3u zlu] + dj[e1v*e3v zlv] )
      !!      bilaplacian:
      !!              diffu = 1/e1u di[ zut ] - 1/(e2u*e3u) dj-1[ e3f zuf ]
      !!              diffv = 1/e2v dj[ zut ] + 1/(e1v*e3v) di-1[ e3f zuf ]
      !!      If ln_sco=F and ln_zps=F, the vertical scale factors in the
      !!      rotational part of the diffusion are simplified
      !!      Add this before trend to the general trend (ua,va):
      !!            (ua,va) = (ua,va) + (diffu,diffv)
      !!      'key_trddyn' defined: the two components of the horizontal
      !!                               diffusion trend are saved.
      !!
      !! ** Action : - Update (ua,va) with the before iso-level biharmonic
      !!               mixing trend.
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk                  ! dummy loop indices
      REAL(wp) ::   zua, zva, zbt, ze2u, ze2v   ! temporary scalar
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zcu, zcv
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zuf, zut, zlu, zlv
!{ DRAKKAR : for no slip accurate G. Hervieux version
      REAL(wp), POINTER, DIMENSION(:,:  ) ::  zrotb  ! temporary workspace
!}
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_ldf_bilap')
      !
      CALL wrk_alloc( jpi, jpj,      zcu, zcv, zrotb    )
      CALL wrk_alloc( jpi, jpj, jpk, zuf, zut, zlu, zlv ) 
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_ldf_bilap : iso-level bilaplacian operator'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF

!!bug gm this should be enough
!!$      zuf(:,:,jpk) = 0.e0
!!$      zut(:,:,jpk) = 0.e0
!!$      zlu(:,:,jpk) = 0.e0
!!$      zlv(:,:,jpk) = 0.e0
      zuf(:,:,:) = 0._wp
      zut(:,:,:) = 0._wp
      zlu(:,:,:) = 0._wp
      zlv(:,:,:) = 0._wp

#if defined key_noslip_accurate
     zrotb(:,:) = 0.e0
#endif

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
#if defined key_noslip_accurate
         zcu(:,:) = ub(:,:,jk)
         zcv(:,:) = vb(:,:,jk)
         CALL  nsacurl(kt, jk, zcu, zcv, zrotb)
#else
         zrotb(:,:)=rotb(:,:,jk)
#endif
         ! Laplacian
         ! ---------

         IF( ln_sco .OR. ln_zps ) THEN   ! s-coordinate or z-coordinate with partial steps
            zuf(:,:,jk) = zrotb(:,:) * fse3f(:,:,jk)
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zlu(ji,jj,jk) = - ( zuf(ji,jj,jk) - zuf(ji,jj-1,jk) ) / ( e2u(ji,jj) * fse3u(ji,jj,jk) )   &
                     &         + ( hdivb(ji+1,jj,jk) - hdivb(ji,jj,jk) ) / e1u(ji,jj)
   
                  zlv(ji,jj,jk) = + ( zuf(ji,jj,jk) - zuf(ji-1,jj,jk) ) / ( e1v(ji,jj) * fse3v(ji,jj,jk) )   &
                     &         + ( hdivb(ji,jj+1,jk) - hdivb(ji,jj,jk) ) / e2v(ji,jj)
               END DO
            END DO
         ELSE                            ! z-coordinate - full step
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zlu(ji,jj,jk) = - ( zrotb (ji  ,jj) - zrotb (ji,jj-1) ) / e2u(ji,jj)   &
                     &         + ( hdivb(ji+1,jj,jk) - hdivb(ji,jj  ,jk) ) / e1u(ji,jj)
   
                  zlv(ji,jj,jk) = + ( zrotb (ji,jj  ) - zrotb (ji-1,jj) ) / e1v(ji,jj)   &
                     &         + ( hdivb(ji,jj+1,jk) - hdivb(ji  ,jj,jk) ) / e2v(ji,jj)
               END DO  
            END DO  
         ENDIF
      END DO
      CALL lbc_lnk( zlu, 'U', -1. )   ;   CALL lbc_lnk( zlv, 'V', -1. )   ! Boundary conditions

         
      DO jk = 1, jpkm1
   
         ! Third derivative
         ! ----------------
         
         ! Multiply by the eddy viscosity coef. (at u- and v-points)
         zlu(:,:,jk) = zlu(:,:,jk) * ( fsahmu(:,:,jk) * (1-nkahm_smag) + nkahm_smag)

         zlv(:,:,jk) = zlv(:,:,jk) * ( fsahmv(:,:,jk) * (1-nkahm_smag) + nkahm_smag)

#if defined key_noslip_accurate
         CALL nsacurl(kt,jk,zlu(:,:,jk),zlv(:,:,jk),zrotb)
         zuf (:,:,jk) = zrotb(:,:) * fse3f(:,:,jk)
#else
         ! Contravariant "laplacian"
         zcu(:,:) = e1u(:,:) * zlu(:,:,jk)
         zcv(:,:) = e2v(:,:) * zlv(:,:,jk)
         
         ! Laplacian curl ( * e3f if s-coordinates or z-coordinate with partial steps)
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zuf(ji,jj,jk) = fmask(ji,jj,jk) * (  zcv(ji+1,jj  ) - zcv(ji,jj)      &
                  &                            - zcu(ji  ,jj+1) + zcu(ji,jj)  )   &
                  &       * fse3f(ji,jj,jk) / ( e1f(ji,jj)*e2f(ji,jj) )
            END DO  
         END DO  
#endif

         ! Laplacian Horizontal fluxes
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zlu(ji,jj,jk) = e2u(ji,jj) * fse3u(ji,jj,jk) * zlu(ji,jj,jk)
               zlv(ji,jj,jk) = e1v(ji,jj) * fse3v(ji,jj,jk) * zlv(ji,jj,jk)
            END DO
         END DO

         ! Laplacian divergence
         DO jj = 2, jpj
            DO ji = fs_2, jpi   ! vector opt.
               zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk)
               zut(ji,jj,jk) = (  zlu(ji,jj,jk) - zlu(ji-1,jj  ,jk)   &
                  &             + zlv(ji,jj,jk) - zlv(ji  ,jj-1,jk) ) / zbt
            END DO
         END DO
      END DO


      ! boundary conditions on the laplacian curl and div (zuf,zut)
!!bug gm no need to do this 2 following lbc...
      CALL lbc_lnk( zuf, 'F', 1. )
      CALL lbc_lnk( zut, 'T', 1. )

      DO jk = 1, jpkm1      
   
         ! Bilaplacian
         ! -----------

         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ze2u = e2u(ji,jj) * fse3u(ji,jj,jk)
               ze2v = e1v(ji,jj) * fse3v(ji,jj,jk)
               ! horizontal biharmonic diffusive trends
               zua = - ( zuf(ji  ,jj,jk) - zuf(ji,jj-1,jk) ) / ze2u   &
                  &  + ( zut(ji+1,jj,jk) - zut(ji,jj  ,jk) ) / e1u(ji,jj)

               zva = + ( zuf(ji,jj  ,jk) - zuf(ji-1,jj,jk) ) / ze2v   &
                  &  + ( zut(ji,jj+1,jk) - zut(ji  ,jj,jk) ) / e2v(ji,jj)
               ! add it to the general momentum trends
               ua(ji,jj,jk) = ua(ji,jj,jk) + zua * ( fsahmu(ji,jj,jk)*nkahm_smag +(1 -nkahm_smag ))
               va(ji,jj,jk) = va(ji,jj,jk) + zva * ( fsahmv(ji,jj,jk)*nkahm_smag +(1 -nkahm_smag ))
            END DO
         END DO

         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      CALL wrk_dealloc( jpi, jpj,      zcu, zcv, zrotb    )
      CALL wrk_dealloc( jpi, jpj, jpk, zuf, zut, zlu, zlv ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_ldf_bilap')
      !
   END SUBROUTINE dyn_ldf_bilap

#if defined key_noslip_accurate
   SUBROUTINE nsacurl( kt, kk, piu, piv, protb  )
      !!----------------------------------------------------------------------
      !!                     *** SUBROUTINE nsacurl ***
      !!
      !! ** Purpose :  laplacien
      !!
      !!
      !! History :  
      !!   NEMO  2.0  ! 2007-12  (G. Hervieux) No slip accurate version Gaelle Hervieux)
      !!         3.4  ! 2012-01  (J.M. Molines) Cleaning ...
      !!----------------------------------------------------------------------
      INTEGER,                     INTENT(in)  :: kt,kk   ! level index
      REAL(wp), DIMENSION(jpi,jpj),INTENT(in)  :: piu,piv ! array in
      REAL(wp), DIMENSION(jpi,jpj),INTENT(out) :: protb   ! array out
      !!
      INTEGER :: ji, jj, jl    ! dummy loop indices
      INTEGER :: ii, ij        ! temporary integer
      INTEGER :: ijt, iju      ! temporary integer   
      REAL(wp), DIMENSION(:,:), POINTER :: zphiu, zphiv 
      REAL(wp), DIMENSION(:,:), POINTER :: zphip1u, zphip1v 
      REAL(wp), DIMENSION(:,:), POINTER :: zphip2u, zphip2v 
      REAL(wp), DIMENSION(:,:), POINTER :: zphim1u, zphim1v 
      REAL(wp), DIMENSION(:,:), POINTER :: zphim2u, zphim2v 
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('nsacurl')
      !
      CALL wrk_alloc( jpi,jpj, zphiu, zphiv, zphip1u, zphip1v, zphip2u, zphip2v , &
              &                              zphim1u, zphim1v, zphim2u, zphim2v   )

      IF( kt == nit000 .AND. kk == 1 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'NSA curl (gh) '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~' 
         IF(lwp) WRITE(numout,*) 
      ENDIF
      
       protb (:,:) = 0.e0 
       zphiu(:,:)=e1u(:,:)*piu(:,:)
       zphiv(:,:)=e2v(:,:)*piv(:,:)

        DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               protb(ji,jj) = ( zphiv(ji+1,jj) - zphiv(ji,jj)   &
                  &           - zphiu(ji,jj+1) + zphiu(ji,jj) ) &
                  &           * fmask(ji,jj,kk) / ( e1f(ji,jj) * e2f(ji,jj) )
            END DO
        END DO  

       ! use auxiliary arrays for mpp computation with jpreci=1 !
       DO jj=1,jpj-1
        zphip1u(:,jj) = zphiu (:,jj+1)
       END DO
       CALL lbc_lnk ( zphip1u , 'U' , -1. )
      
       DO ji=1,jpi-1
        zphip1v(ji,:) = zphiv (ji+1, :)
       END DO
       CALL lbc_lnk (zphip1v ,'V' ,-1. )

       DO jj=2,jpj
        zphim1u(:,jj) = zphiu (:,jj-1)
       END DO
       CALL lbc_lnk (zphim1u,'U',-1. )

       DO ji=2,jpi
        zphim1v(ji,:) = zphiv (ji-1, :)
       END DO
       CALL lbc_lnk (zphim1v,'V',-1. )

       DO jj=1,jpj-1
        zphip2u(:,jj) = zphip1u (:,jj+1)
       END DO
       CALL lbc_lnk (zphip2u,'U',-1. )
       
       DO ji=1,jpi-1
        zphip2v(ji,:) = zphip1v (ji+1, :)
       END DO
       CALL lbc_lnk (zphip2v,'V',-1. )

       DO jj=2,jpj
        zphim2u(:,jj) = zphim1u (:,jj-1)
       END DO
       CALL lbc_lnk (zphim2u,'U',-1. )

       DO ji=2,jpi
        zphim2v(ji,:) = zphim1v (ji-1, :)
       END DO
       CALL lbc_lnk (zphim2v,'V',-1. )

        ! south 
        DO jl = 1, npcoa(3,kk)  
            ji = nicoa(jl,3,kk)
            jj = njcoa(jl,3,kk)
            protb(ji,jj) = -( 4*zphip1u(ji,jj)  - zphip2u(ji,jj) + 0.2*zphip2u(ji,jj+1) )   &
               &/ ( e1f(ji  ,jj  ) * e2f(ji  ,jj  )   )        * fmask(ji,jj,kk) * 0.5    
        END DO 

        ! north 
        DO jl = 1, npcoa(4,kk)  
            ji = nicoa(jl,4,kk)
            jj = njcoa(jl,4,kk)
            protb(ji,jj) = -( -4*zphiu( ji, jj) + zphim1u(ji,jj) - 0.2*zphim2u( ji ,jj ) )  &
              &/ ( e1f(ji  ,jj  ) * e2f(ji  ,jj  ) ) * fmask(ji,jj,kk) * 0.5  
        END DO
       
        ! west
        DO jl = 1, npcoa(1,kk)  
            ji = nicoa(jl,1,kk)
            jj = njcoa(jl,1,kk)     
            protb(ji,jj) = (  4*zphip1v (ji,jj) - zphip2v(ji, jj) +  0.2*zphip2v(ji+1,jj)  ) &
              &  /    ( e1f(ji  ,jj  ) * e2f(ji  ,jj) ) * fmask(ji,jj,kk) * 0.5
        END DO        

        ! east 
        DO jl = 1, npcoa(2,kk)  
            ji = nicoa(jl,2,kk)
            jj = njcoa(jl,2,kk)
            protb(ji,jj) =(- 4* zphiv(ji,jj) + zphim1v(ji,jj)  - 0.2* zphim2v(ji,jj )  )   &
              &/    ( e1f(ji  ,jj  ) * e2f(ji  ,jj  )   ) * fmask(ji,jj,kk) * 0.5 
        END DO

        !! west-south 
        DO jl = 1, npcoa(5,kk)  
            ji = nicoa(jl,5,kk)
            jj = njcoa(jl,5,kk)   
            protb(ji,jj) = ((  4*zphip1v(ji,jj) - zphip2v (ji,jj) +   0.2*zphip2v(ji+1,jj)  ) &
              &  -  (  4*zphip1u(ji,jj) - zphip2u(ji,jj)  +   0.2*zphip2u(ji,jj+1)  ))&
              &  /      ( e1f(ji  ,jj  ) * e2f(ji  ,jj  )   ) * fmask(ji,jj,kk) * 0.5 
        END DO
 
        ! west-north
        DO jl = 1, npcoa(6,kk)  
            ji = nicoa(jl,6,kk)
            jj = njcoa(jl,6,kk)  
            protb(ji,jj) = ((   4*zphip1v(ji,jj) -  zphip2v(ji,jj) + 0.2*zphip2v(ji+1,jj)  )    &
              &  - ( - 4*zphiu  (ji,jj) +  zphim1u(ji,jj) - 0.2*zphim2u(ji  ,jj)  ))   &
              &  / (     e1f(ji  ,jj  ) * e2f(ji  ,jj  )   ) * fmask(ji,jj,kk) * 0.5
        END DO

        ! east-south 
        DO jl = 1, npcoa(7,kk)  
            ji = nicoa(jl,7,kk)
            jj = njcoa(jl,7,kk)
            protb(ji,jj) =((- 4*zphiv  (ji,jj) +   zphim1v(ji,jj) -  0.2* zphim2v(ji,jj  )  )   &
              &     - (  4*zphip1u(ji,jj) -   zphip2u(ji,jj) +  0.2* zphip2u(ji,jj+1)  ))&
              &/     ( e1f(ji  ,jj  ) * e2f(ji  ,jj  )  ) * fmask(ji,jj,kk) * 0.5
        END DO

        !! east-north
        DO jl = 1, npcoa(8,kk)  
            ji = nicoa(jl,8,kk)
            jj = njcoa(jl,8,kk)  
            protb(ji,jj) =((- 4*zphiv(ji,jj) + zphim1v(ji,jj) - 0.2*zphim2v(ji,jj)   )  &
              &     - (- 4*zphiu(ji,jj) + zphim1u(ji,jj) - 0.2*zphim2u(ji,jj) ))   &
              &     / (     e1f(ji  ,jj  ) * e2f(ji  ,jj  )   ) * fmask(ji,jj,kk) * 0.5
        END DO
      CALL lbc_lnk( protb , 'F', -1. )     ! F-point,    sign change 

      CALL wrk_dealloc( jpi,jpj, zphiu, zphiv, zphip1u, zphip1v, zphip2u, zphip2v , &
              &                                zphim1u, zphim1v, zphim2u, zphim2v   )
      !
      IF( nn_timing == 1 )  CALL timing_stop('nsacurl')
   END SUBROUTINE nsacurl
#else

   SUBROUTINE nsacurl( kt, kk, piu, piv, protb  )
      !!----------------------------------------------------------------------
      !!                     *** SUBROUTINE nsacurl ***
      !!
      !! ** Purpose :  dummy routine 
      !!
      !!----------------------------------------------------------------------
      INTEGER,  INTENT(in)  :: kt,kk   ! level index
      REAL(wp), INTENT(in)  :: piu,piv ! array in
      REAL(wp), INTENT(out) :: protb   ! array out
   PRINT *, ' ERROR : you must no go this way '

   END SUBROUTINE nsacurl
#endif
         
   !!======================================================================
END MODULE dynldf_bilap
