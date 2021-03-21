#define SPONGE && define SPONGE_TOP

Module agrif_opa_sponge
#if defined key_agrif  && ! defined key_offline
   USE par_oce
   USE oce
   USE dom_oce
   USE in_out_manager
   USE agrif_oce
   USE wrk_nemo  

   IMPLICIT NONE
   PRIVATE

   PUBLIC Agrif_Sponge, Agrif_Sponge_Tra, Agrif_Sponge_Dyn, interptsn, interpun, interpvn

  !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/NST 3.3 , NEMO Consortium (2010)
   !! $Id: agrif_opa_sponge.F90 4009 2013-08-28 14:18:17Z cbricaud $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   CONTAINS

   SUBROUTINE Agrif_Sponge_Tra
      !!---------------------------------------------
      !!   *** ROUTINE Agrif_Sponge_Tra ***
      !!---------------------------------------------
      !!
      INTEGER :: ji,jj,jk,jn
      REAL(wp) :: timecoeff
      REAL(wp) :: ztsa, zabe1, zabe2, zbtr
      REAL(wp), POINTER, DIMENSION(:,:    ) :: ztu, ztv
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: ztab
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: tsbdiff

#if defined SPONGE
      CALL wrk_alloc( jpi, jpj, ztu, ztv )
      CALL wrk_alloc( jpi, jpj, jpk, jpts, ztab, tsbdiff  )

      timecoeff = REAL(Agrif_NbStepint(),wp)/Agrif_rhot()

      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = .TRUE.
      ztab = 0.e0
      CALL Agrif_Bc_Variable(ztab, tsa_id,calledweight=timecoeff,procname=interptsn)
      Agrif_UseSpecialValue = .FALSE.

      tsbdiff(:,:,:,:) = tsb(:,:,:,:) - ztab(:,:,:,:)

      CALL Agrif_Sponge

      DO jn = 1, jpts
         DO jk = 1, jpkm1
            !
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zabe1 = umask(ji,jj,jk) * spe1ur(ji,jj) * fse3u(ji,jj,jk)
                  zabe2 = vmask(ji,jj,jk) * spe2vr(ji,jj) * fse3v(ji,jj,jk)
                  ztu(ji,jj) = zabe1 * ( tsbdiff(ji+1,jj  ,jk,jn) - tsbdiff(ji,jj,jk,jn) )
                  ztv(ji,jj) = zabe2 * ( tsbdiff(ji  ,jj+1,jk,jn) - tsbdiff(ji,jj,jk,jn) )
               ENDDO
            ENDDO

            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zbtr = spbtr2(ji,jj) / fse3t(ji,jj,jk)
                  ! horizontal diffusive trends
                  ztsa = zbtr * (  ztu(ji,jj) - ztu(ji-1,jj  )   &
                  &              + ztv(ji,jj) - ztv(ji  ,jj-1)  )
                  ! add it to the general tracer trends
                  tsa(ji,jj,jk,jn) = tsa(ji,jj,jk,jn) + ztsa
               END DO
            END DO
            !
         ENDDO
      ENDDO

      CALL wrk_dealloc( jpi, jpj, ztu, ztv )
      CALL wrk_dealloc( jpi, jpj, jpk, jpts, ztab, tsbdiff  )
#endif

   END SUBROUTINE Agrif_Sponge_Tra

   SUBROUTINE Agrif_Sponge_dyn
      !!---------------------------------------------
      !!   *** ROUTINE Agrif_Sponge_dyn ***
      !!---------------------------------------------
      !!
      INTEGER :: ji,jj,jk
      REAL(wp) :: timecoeff
      REAL(wp) :: ze2u, ze1v, zua, zva, zbtr
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ubdiff, vbdiff
      REAL(wp), POINTER, DIMENSION(:,:,:) :: rotdiff, hdivdiff
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztab

#if defined SPONGE
      CALL wrk_alloc( jpi, jpj, jpk, ztab, ubdiff, vbdiff, rotdiff, hdivdiff )

      timecoeff = REAL(Agrif_NbStepint(),wp)/Agrif_rhot()

      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = ln_spc_dyn
      ztab = 0.e0
      CALL Agrif_Bc_Variable(ztab, ua_id,calledweight=timecoeff,procname=interpun)
      Agrif_UseSpecialValue = .FALSE.

      ubdiff(:,:,:) = ( ub(:,:,:) - ztab(:,:,:) ) * umask(:,:,:)

      ztab = 0.e0
      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = ln_spc_dyn
      CALL Agrif_Bc_Variable(ztab, va_id,calledweight=timecoeff,procname=interpvn)
      Agrif_UseSpecialValue = .FALSE.

      vbdiff(:,:,:) = ( vb(:,:,:) - ztab(:,:,:) ) * vmask(:,:,:)

      CALL Agrif_Sponge

      DO jk = 1,jpkm1
         ubdiff(:,:,jk) = ubdiff(:,:,jk) * spe1ur2(:,:)
         vbdiff(:,:,jk) = vbdiff(:,:,jk) * spe2vr2(:,:)
      ENDDO
      
      hdivdiff = 0.
      rotdiff = 0.

      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         !                                             ! --------
         ! Horizontal divergence                       !   div
         !                                             ! --------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zbtr = spbtr2(ji,jj) / fse3t(ji,jj,jk)
               hdivdiff(ji,jj,jk) =  (  e2u(ji  ,jj  ) * fse3u(ji  ,jj  ,jk) * ubdiff(ji  ,jj  ,jk)     &
                  &                   - e2u(ji-1,jj  ) * fse3u(ji-1,jj  ,jk) * ubdiff(ji-1,jj  ,jk)     &
                  &                   + e1v(ji  ,jj  ) * fse3v(ji  ,jj  ,jk) * vbdiff(ji  ,jj  ,jk)     &
                  &                   - e1v(ji  ,jj-1) * fse3v(ji  ,jj-1,jk) * vbdiff(ji  ,jj-1,jk)  ) * zbtr
            END DO
         END DO

         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               zbtr = spbtr3(ji,jj) * fse3f(ji,jj,jk)
               rotdiff(ji,jj,jk) = (  e2v(ji+1,jj  ) * vbdiff(ji+1,jj  ,jk) - e2v(ji,jj) * vbdiff(ji,jj,jk)    &
                  &                 - e1u(ji  ,jj+1) * ubdiff(ji  ,jj+1,jk) + e1u(ji,jj) * ubdiff(ji,jj,jk)  ) &
                  &               * fmask(ji,jj,jk) * zbtr
            END DO
         END DO

      ENDDO

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ! horizontal diffusive trends
               zua = - ( rotdiff (ji  ,jj,jk) - rotdiff (ji,jj-1,jk) ) / ( e2u(ji,jj) * fse3u(ji,jj,jk) )   &
                     + ( hdivdiff(ji+1,jj,jk) - hdivdiff(ji,jj  ,jk) ) / e1u(ji,jj)

               zva = + ( rotdiff (ji,jj  ,jk) - rotdiff (ji-1,jj,jk) ) / ( e1v(ji,jj) * fse3v(ji,jj,jk) )   &
                     + ( hdivdiff(ji,jj+1,jk) - hdivdiff(ji  ,jj,jk) ) / e2v(ji,jj)
               ! add it to the general momentum trends
               ua(ji,jj,jk) = ua(ji,jj,jk) + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zva
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj, jpk, ztab, ubdiff, vbdiff, rotdiff, hdivdiff )
#endif

   END SUBROUTINE Agrif_Sponge_dyn

   SUBROUTINE Agrif_Sponge
      !!---------------------------------------------
      !!   *** ROUTINE  Agrif_Sponge ***
      !!---------------------------------------------
      INTEGER  :: ji,jj,jk
      INTEGER  :: ispongearea, ilci, ilcj
      LOGICAL  :: ll_spdone
      REAL(wp) :: z1spongearea, zramp
      REAL(wp), POINTER, DIMENSION(:,:) :: ztabramp

#if defined SPONGE || defined SPONGE_TOP
      ll_spdone=.TRUE.
      IF (( .NOT. spongedoneT ).OR.( .NOT. spongedoneU )) THEN
         ! Define ramp from boundaries towards domain interior
         ! at T-points
         ! Store it in ztabramp
         ll_spdone=.FALSE.

         CALL wrk_alloc( jpi, jpj, ztabramp )

         ispongearea  = 2 + 2 * Agrif_irhox()
         ilci = nlci - ispongearea
         ilcj = nlcj - ispongearea 
         z1spongearea = 1._wp / REAL( ispongearea - 2 )
         spbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )

         ztabramp(:,:) = 0.

         IF( (nbondi == -1) .OR. (nbondi == 2) ) THEN
            DO jj = 1, jpj
               IF ( umask(2,jj,1) == 1._wp ) THEN
                 DO ji = 2, ispongearea                  
                    ztabramp(ji,jj) = ( ispongearea-ji ) * z1spongearea
                 END DO
               ENDIF
            ENDDO
         ENDIF

         IF( (nbondi == 1) .OR. (nbondi == 2) ) THEN
            DO jj = 1, jpj
               IF ( umask(nlci-2,jj,1) == 1._wp ) THEN
                  DO ji = ilci+1,nlci-1
                     zramp = (ji - (ilci+1) ) * z1spongearea
                     ztabramp(ji,jj) = MAX( ztabramp(ji,jj), zramp )
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

         IF( (nbondj == -1) .OR. (nbondj == 2) ) THEN
            DO ji = 1, jpi
               IF ( vmask(ji,2,1) == 1._wp ) THEN
                  DO jj = 2, ispongearea
                     zramp = ( ispongearea-jj ) * z1spongearea
                     ztabramp(ji,jj) = MAX( ztabramp(ji,jj), zramp )
                  END DO
               ENDIF
            ENDDO
         ENDIF

         IF( (nbondj == 1) .OR. (nbondj == 2) ) THEN
            DO ji = 1, jpi
               IF ( vmask(ji,nlcj-2,1) == 1._wp ) THEN
                  DO jj = ilcj+1,nlcj-1
                     zramp = (jj - (ilcj+1) ) * z1spongearea
                     ztabramp(ji,jj) = MAX( ztabramp(ji,jj), zramp )
                  END DO
               ENDIF
            ENDDO
         ENDIF

      ENDIF

      ! Tracers
      IF( .NOT. spongedoneT ) THEN
         spe1ur(:,:) = 0.
         spe2vr(:,:) = 0.

         IF( (nbondi == -1) .OR. (nbondi == 2) ) THEN
            spe1ur(2:ispongearea-1,:       ) = visc_tra                                        &
               &                             *    0.5 * (  ztabramp(2:ispongearea-1,:      )   &
               &                                         + ztabramp(3:ispongearea  ,:      ) ) &
               &                             * e2u(2:ispongearea-1,:) / e1u(2:ispongearea-1,:)

            spe2vr(2:ispongearea  ,1:jpjm1 ) = visc_tra                                        &
               &                             *    0.5 * (  ztabramp(2:ispongearea  ,1:jpjm1)   &
               &                                         + ztabramp(2:ispongearea,2  :jpj  ) ) &
               &                             * e1v(2:ispongearea,1:jpjm1) / e2v(2:ispongearea,1:jpjm1)
         ENDIF

         IF( (nbondi == 1) .OR. (nbondi == 2) ) THEN
            spe1ur(ilci+1:nlci-2,:        ) = visc_tra                                   &
               &                            * 0.5 * (  ztabramp(ilci+1:nlci-2,:      )   & 
               &                                     + ztabramp(ilci+2:nlci-1,:      ) ) &
               &                            * e2u(ilci+1:nlci-2,:) / e1u(ilci+1:nlci-2,:)

            spe2vr(ilci+1:nlci-1,1:jpjm1  )  = visc_tra                                  &
               &                            * 0.5 * (  ztabramp(ilci+1:nlci-1,1:jpjm1)   & 
               &                                     + ztabramp(ilci+1:nlci-1,2:jpj  ) ) & 
               &                            * e1v(ilci+1:nlci-1,1:jpjm1) / e2v(ilci+1:nlci-1,1:jpjm1)
         ENDIF

         IF( (nbondj == -1) .OR. (nbondj == 2) ) THEN
            spe1ur(1:jpim1,2:ispongearea  ) = visc_tra                                     &
               &                            * 0.5 * (  ztabramp(1:jpim1,2:ispongearea  )   & 
               &                                     + ztabramp(2:jpi  ,2:ispongearea  ) ) &
               &                            * e2u(1:jpim1,2:ispongearea) / e1u(1:jpim1,2:ispongearea)
   
            spe2vr(:      ,2:ispongearea-1) = visc_tra                                     &
               &                            * 0.5 * (  ztabramp(:      ,2:ispongearea-1)   &
               &                                     + ztabramp(:      ,3:ispongearea  ) ) &
               &                            * e1v(:,2:ispongearea-1) / e2v(:,2:ispongearea-1)
         ENDIF

         IF( (nbondj == 1) .OR. (nbondj == 2) ) THEN
            spe1ur(1:jpim1,ilcj+1:nlcj-1) = visc_tra                                   &
               &                          * 0.5 * (  ztabramp(1:jpim1,ilcj+1:nlcj-1)   &
               &                                   + ztabramp(2:jpi  ,ilcj+1:nlcj-1) ) &
               &                                * e2u(1:jpim1,ilcj+1:nlcj-1) / e1u(1:jpim1,ilcj+1:nlcj-1)

            spe2vr(:      ,ilcj+1:nlcj-2) = visc_tra                                   &
               &                          * 0.5 * (  ztabramp(:      ,ilcj+1:nlcj-2)   &
               &                                   + ztabramp(:      ,ilcj+2:nlcj-1) ) &
               &                                * e1v(:,ilcj+1:nlcj-2) / e2v(:,ilcj+1:nlcj-2)
         ENDIF
         spongedoneT = .TRUE.
      ENDIF

      ! Dynamics
      IF( .NOT. spongedoneU ) THEN
         spe1ur2(:,:) = 0.
         spe2vr2(:,:) = 0.

         IF( (nbondi == -1) .OR. (nbondi == 2) ) THEN
            spe1ur2(2:ispongearea-1,:      ) = visc_dyn                                   &
               &                             * 0.5 * (  ztabramp(2:ispongearea-1,:      ) &
               &                                      + ztabramp(3:ispongearea  ,:      ) )
            spe2vr2(2:ispongearea  ,1:jpjm1) = visc_dyn                                   &
               &                             * 0.5 * (  ztabramp(2:ispongearea  ,1:jpjm1) &
               &                                      + ztabramp(2:ispongearea  ,2:jpj  ) ) 
         ENDIF

         IF( (nbondi == 1) .OR. (nbondi == 2) ) THEN
            spe1ur2(ilci+1:nlci-2  ,:      ) = visc_dyn                                   &
               &                             * 0.5 * (  ztabramp(ilci+1:nlci-2, :       ) &
               &                                      + ztabramp(ilci+2:nlci-1, :       ) )                      
            spe2vr2(ilci+1:nlci-1  ,1:jpjm1) = visc_dyn                                   &
               &                             * 0.5 * (  ztabramp(ilci+1:nlci-1,1:jpjm1  ) &
               &                                      + ztabramp(ilci+1:nlci-1,2:jpj    ) ) 
         ENDIF

         IF( (nbondj == -1) .OR. (nbondj == 2) ) THEN
            spe1ur2(1:jpim1,2:ispongearea  ) = visc_dyn                                   &  
               &                             * 0.5 * (  ztabramp(1:jpim1,2:ispongearea  ) &
               &                                      + ztabramp(2:jpi  ,2:ispongearea  ) ) 
            spe2vr2(:      ,2:ispongearea-1) = visc_dyn                                   &
               &                             * 0.5 * (  ztabramp(:      ,2:ispongearea-1) &
               &                                      + ztabramp(:      ,3:ispongearea  ) )
         ENDIF

         IF( (nbondj == 1) .OR. (nbondj == 2) ) THEN
            spe1ur2(1:jpim1,ilcj+1:nlcj-1  ) = visc_dyn                                   &
               &                             * 0.5 * (  ztabramp(1:jpim1,ilcj+1:nlcj-1  ) &
               &                                      + ztabramp(2:jpi  ,ilcj+1:nlcj-1  ) ) 
            spe2vr2(:      ,ilcj+1:nlcj-2  ) = visc_dyn                                   &
               &                             * 0.5 * (  ztabramp(:      ,ilcj+1:nlcj-2  ) &
               &                                      + ztabramp(:      ,ilcj+2:nlcj-1  ) )
         ENDIF
         spongedoneU = .TRUE.
         spbtr3(:,:) = 1. / ( e1f(:,:) * e2f(:,:) )
      ENDIF
      !
      IF (.NOT.ll_spdone) CALL wrk_dealloc( jpi, jpj, ztabramp )
      !
#endif

   END SUBROUTINE Agrif_Sponge

   SUBROUTINE interptsn(tabres,i1,i2,j1,j2,k1,k2,n1,n2)
      !!---------------------------------------------
      !!   *** ROUTINE interptsn ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2,n1,n2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) :: tabres

      tabres(i1:i2,j1:j2,k1:k2,n1:n2) = tsn(i1:i2,j1:j2,k1:k2,n1:n2)

   END SUBROUTINE interptsn

   SUBROUTINE interpun(tabres,i1,i2,j1,j2,k1,k2)
      !!---------------------------------------------
      !!   *** ROUTINE interpun ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres

      tabres(i1:i2,j1:j2,k1:k2) = un(i1:i2,j1:j2,k1:k2)

   END SUBROUTINE interpun

   SUBROUTINE interpvn(tabres,i1,i2,j1,j2,k1,k2)
      !!---------------------------------------------
      !!   *** ROUTINE interpvn ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres

      tabres(i1:i2,j1:j2,k1:k2) = vn(i1:i2,j1:j2,k1:k2)

   END SUBROUTINE interpvn

#else
CONTAINS

   SUBROUTINE agrif_opa_sponge_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_OPA_sponge_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_opa_sponge : You should not have seen this print! error?'
   END SUBROUTINE agrif_opa_sponge_empty
#endif

END MODULE agrif_opa_sponge
