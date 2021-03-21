#define SPONGE_TOP

Module agrif_top_sponge
#if defined key_agrif && defined key_top
   USE par_oce
   USE oce
   USE dom_oce
   USE in_out_manager
   USE agrif_oce
   USE agrif_opa_sponge
   USE trc
   USE lib_mpp
   USE wrk_nemo  

   IMPLICIT NONE
   PRIVATE

   PUBLIC Agrif_Sponge_Trc, interptrn

  !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/NST 3.3 , NEMO Consortium (2010)
   !! $Id: agrif_top_sponge.F90 3680 2012-11-27 14:42:24Z rblod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   CONTAINS

   SUBROUTINE Agrif_Sponge_Trc
      !!---------------------------------------------
      !!   *** ROUTINE Agrif_Sponge_Trc ***
      !!---------------------------------------------
      !! 
      INTEGER :: ji,jj,jk,jn
      REAL(wp) :: timecoeff
      REAL(wp) :: ztra, zabe1, zabe2, zbtr
      REAL(wp), POINTER, DIMENSION(:,:) :: ztru, ztrv
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: ztabr
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: trbdiff

#if defined SPONGE_TOP
      CALL wrk_alloc( jpi, jpj, ztru, ztrv )
      CALL wrk_alloc( jpi, jpj, jpk, jptra, ztabr, trbdiff )

      timecoeff = REAL(Agrif_NbStepint(),wp)/Agrif_rhot()

      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = .TRUE.
      ztabr = 0.e0
      CALL Agrif_Bc_Variable(ztabr, tra_id,calledweight=timecoeff,procname=interptrn)
      Agrif_UseSpecialValue = .FALSE.

      trbdiff(:,:,:,:) = trb(:,:,:,:) - ztabr(:,:,:,:)

      CALL Agrif_sponge

      DO jn = 1, jptra
         DO jk = 1, jpkm1
            !
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zabe1 = umask(ji,jj,jk) * spe1ur(ji,jj) * fse3u(ji,jj,jk)
                  zabe2 = vmask(ji,jj,jk) * spe2vr(ji,jj) * fse3v(ji,jj,jk)
                  ztru(ji,jj) = zabe1 * ( trbdiff(ji+1,jj  ,jk,jn) - trbdiff(ji,jj,jk,jn) )
                  ztrv(ji,jj) = zabe2 * ( trbdiff(ji  ,jj+1,jk,jn) - trbdiff(ji,jj,jk,jn) )
               ENDDO
            ENDDO

            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  zbtr = spbtr2(ji,jj) / fse3t(ji,jj,jk)
                  ! horizontal diffusive trends
                  ztra = zbtr * ( ztru(ji,jj) - ztru(ji-1,jj) + ztrv(ji,jj) - ztrv(ji,jj-1)  )
                  ! add it to the general tracer trends
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
               END DO
            END DO
            !
         ENDDO
      ENDDO
 
      CALL wrk_dealloc( jpi, jpj, ztru, ztrv )
      CALL wrk_dealloc( jpi, jpj, jpk, jptra, trbdiff, ztabr )

#endif

   END SUBROUTINE Agrif_Sponge_Trc

   SUBROUTINE interptrn(tabres,i1,i2,j1,j2,k1,k2,n1,n2)
      !!---------------------------------------------
      !!   *** ROUTINE interptn ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2,n1,n2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) :: tabres
      !
      tabres(i1:i2,j1:j2,k1:k2,n1:n2) = trn(i1:i2,j1:j2,k1:k2,n1:n2)

   END SUBROUTINE interptrn

#else
CONTAINS

   SUBROUTINE agrif_top_sponge_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_top_sponge_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_top_sponge : You should not have seen this print! error?'
   END SUBROUTINE agrif_top_sponge_empty
#endif

END MODULE agrif_top_sponge
