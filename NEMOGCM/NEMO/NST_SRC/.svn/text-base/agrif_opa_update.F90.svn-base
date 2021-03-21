#define TWO_WAY

MODULE agrif_opa_update
#if defined key_agrif  && ! defined key_offline
   USE par_oce
   USE oce
   USE dom_oce
   USE agrif_oce
   USE in_out_manager  ! I/O manager
   USE lib_mpp
   USE wrk_nemo  

   IMPLICIT NONE
   PRIVATE

   PUBLIC Agrif_Update_Tra, Agrif_Update_Dyn

   INTEGER, PUBLIC :: nbcline = 0

   !!----------------------------------------------------------------------
   !! NEMO/NST 3.3 , NEMO Consortium (2010)
   !! $Id: agrif_opa_update.F90 3918 2013-06-13 10:50:37Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE Agrif_Update_Tra( kt )
      !!---------------------------------------------
      !!   *** ROUTINE Agrif_Update_Tra ***
      !!---------------------------------------------
      !!
      INTEGER, INTENT(in) :: kt
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: ztab

       
      IF((Agrif_NbStepint() .NE. (Agrif_irhot()-1)).AND.(kt /= 0)) RETURN
#if defined TWO_WAY
      CALL wrk_alloc( jpi, jpj, jpk, jpts, ztab )

      Agrif_UseSpecialValueInUpdate = .TRUE.
      Agrif_SpecialValueFineGrid = 0.

      IF (MOD(nbcline,nbclineupdate) == 0) THEN
         CALL Agrif_Update_Variable(ztab,tsn_id, procname=updateTS)
      ELSE
         CALL Agrif_Update_Variable(ztab,tsn_id,locupdate=(/0,2/), procname=updateTS)
      ENDIF

      Agrif_UseSpecialValueInUpdate = .FALSE.

      CALL wrk_dealloc( jpi, jpj, jpk, jpts, ztab )
#endif

   END SUBROUTINE Agrif_Update_Tra

   SUBROUTINE Agrif_Update_Dyn( kt )
      !!---------------------------------------------
      !!   *** ROUTINE Agrif_Update_Dyn ***
      !!---------------------------------------------
      !!
      INTEGER, INTENT(in) :: kt
      REAL(wp), POINTER, DIMENSION(:,:) :: ztab2d
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztab


      IF ((Agrif_NbStepint() .NE. (Agrif_irhot()-1)).AND.(kt /= 0)) Return
#if defined TWO_WAY
      CALL wrk_alloc( jpi, jpj,      ztab2d )
      CALL wrk_alloc( jpi, jpj, jpk, ztab   )

      IF (mod(nbcline,nbclineupdate) == 0) THEN
         CALL Agrif_Update_Variable(ztab,un_id,procname = updateU)
         CALL Agrif_Update_Variable(ztab,vn_id,procname = updateV)
      ELSE
         CALL Agrif_Update_Variable(ztab,un_id,locupdate=(/0,1/),procname = updateU)
         CALL Agrif_Update_Variable(ztab,vn_id,locupdate=(/0,1/),procname = updateV)         
      ENDIF

      CALL Agrif_Update_Variable(ztab2d,e1u_id,procname = updateU2d)
      CALL Agrif_Update_Variable(ztab2d,e2v_id,procname = updateV2d)  

      nbcline = nbcline + 1

      Agrif_UseSpecialValueInUpdate = ln_spc_dyn
      Agrif_SpecialValueFineGrid = 0.
      CALL Agrif_Update_Variable(ztab2d,sshn_id,procname = updateSSH)
      Agrif_UseSpecialValueInUpdate = .FALSE.

      CALL wrk_dealloc( jpi, jpj,      ztab2d )
      CALL wrk_dealloc( jpi, jpj, jpk, ztab   )

!Done in step
!      CALL Agrif_ChildGrid_To_ParentGrid()
!      CALL recompute_diags( kt )
!      CALL Agrif_ParentGrid_To_ChildGrid()

#endif

   END SUBROUTINE Agrif_Update_Dyn

   SUBROUTINE recompute_diags( kt )
      !!---------------------------------------------
      !!   *** ROUTINE recompute_diags ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: kt

   END SUBROUTINE recompute_diags

   SUBROUTINE updateTS( tabres, i1, i2, j1, j2, k1, k2, n1, n2, before )
      !!---------------------------------------------
      !!           *** ROUTINE updateT ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"

      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2,n1,n2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) :: tabres
      LOGICAL, iNTENT(in) :: before

      INTEGER :: ji,jj,jk,jn

      IF (before) THEN
         DO jn = n1,n2
            DO jk=k1,k2
               DO jj=j1,j2
                  DO ji=i1,i2
                     tabres(ji,jj,jk,jn) = tsn(ji,jj,jk,jn)
                  END DO
               END DO
            END DO
         END DO
      ELSE
         DO jn = n1,n2
            DO jk=k1,k2
               DO jj=j1,j2
                  DO ji=i1,i2
                     IF( tabres(ji,jj,jk,jn) .NE. 0. ) THEN 
                         tsn(ji,jj,jk,jn) = tabres(ji,jj,jk,jn) * tmask(ji,jj,jk)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      ENDIF

   END SUBROUTINE updateTS

   SUBROUTINE updateu( tabres, i1, i2, j1, j2, k1, k2, before )
      !!---------------------------------------------
      !!           *** ROUTINE updateu ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"

      INTEGER, INTENT(in) :: i1, i2, j1, j2, k1, k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before

      INTEGER :: ji, jj, jk
      REAL(wp) :: zrhoy

      IF (before) THEN
         zrhoy = Agrif_Rhoy()
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  tabres(ji,jj,jk) = e2u(ji,jj) * un(ji,jj,jk)
                  tabres(ji,jj,jk) = tabres(ji,jj,jk) * fse3u(ji,jj,jk)
               END DO
            END DO
         END DO
         tabres = zrhoy * tabres
      ELSE
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  un(ji,jj,jk) = tabres(ji,jj,jk) / (e2u(ji,jj))
                  un(ji,jj,jk) = un(ji,jj,jk) * umask(ji,jj,jk)
                  un(ji,jj,jk) = un(ji,jj,jk) / fse3u(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

   END SUBROUTINE updateu

   SUBROUTINE updatev( tabres, i1, i2, j1, j2, k1, k2, before )
      !!---------------------------------------------
      !!           *** ROUTINE updatev ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"

      INTEGER :: i1,i2,j1,j2,k1,k2
      INTEGER :: ji,jj,jk
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2) :: tabres
      LOGICAL :: before

      REAL(wp) :: zrhox

      IF (before) THEN
         zrhox = Agrif_Rhox()
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  tabres(ji,jj,jk) = e1v(ji,jj) * vn(ji,jj,jk)
                  tabres(ji,jj,jk) = tabres(ji,jj,jk) * fse3v(ji,jj,jk)
               END DO
            END DO
         END DO
         tabres = zrhox * tabres
      ELSE
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  vn(ji,jj,jk) = tabres(ji,jj,jk) / (e1v(ji,jj))
                  vn(ji,jj,jk) = vn(ji,jj,jk) * vmask(ji,jj,jk)
                  vn(ji,jj,jk) = vn(ji,jj,jk) / fse3v(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

   END SUBROUTINE updatev

   SUBROUTINE updateu2d( tabres, i1, i2, j1, j2, before )
      !!---------------------------------------------
      !!          *** ROUTINE updateu2d ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"

      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before

      INTEGER :: ji, jj, jk
      REAL(wp) :: zrhoy
      REAL(wp) :: zhinv

      IF (before) THEN
         zrhoy = Agrif_Rhoy()
         DO jk = 1,jpkm1
            DO jj=j1,j2
               DO ji=i1,i2
                  tabres(ji,jj) = tabres(ji,jj) + fse3u(ji,jj,jk) * un(ji,jj,jk)
               END DO
            END DO
         END DO
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj) = tabres(ji,jj) * e2u(ji,jj)
            END DO
         END DO
         tabres = zrhoy * tabres
      ELSE
         DO jj=j1,j2
            DO ji=i1,i2
               IF(umask(ji,jj,1) .NE. 0.) THEN             
                  spgu(ji,jj) = 0.e0
                  DO jk=1,jpk
                     spgu(ji,jj) = spgu(ji,jj) + fse3u(ji,jj,jk) * un(ji,jj,jk)
                  END DO
                  spgu(ji,jj) = spgu(ji,jj) * e2u(ji,jj)
                  zhinv = (tabres(ji,jj)-spgu(ji,jj))/(hu(ji,jj)*e2u(ji,jj))
                  Do jk=1,jpk              
                     un(ji,jj,jk) = un(ji,jj,jk) + zhinv
                     un(ji,jj,jk) = un(ji,jj,jk) * umask(ji,jj,jk)            
                  END DO
               ENDIF
            END DO
         END DO
      ENDIF

   END SUBROUTINE updateu2d

   SUBROUTINE updatev2d( tabres, i1, i2, j1, j2, before )
      !!---------------------------------------------
      !!          *** ROUTINE updatev2d ***
      !!---------------------------------------------

      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before

      INTEGER :: ji, jj, jk
      REAL(wp) :: zrhox
      REAL(wp) :: zhinv

      IF (before) THEN
         zrhox = Agrif_Rhox()
         tabres = 0.e0
         DO jk = 1,jpkm1
            DO jj=j1,j2
               DO ji=i1,i2
                  tabres(ji,jj) = tabres(ji,jj) + fse3v(ji,jj,jk) * vn(ji,jj,jk)
               END DO
            END DO
         END DO
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj) = tabres(ji,jj) * e1v(ji,jj)
            END DO
         END DO
         tabres = zrhox * tabres
      ELSE
         DO jj=j1,j2
            DO ji=i1,i2
               IF(vmask(ji,jj,1) .NE. 0.) THEN             
                  spgv(ji,jj) = 0.
                  DO jk=1,jpk
                     spgv(ji,jj) = spgv(ji,jj) + fse3v(ji,jj,jk) * vn(ji,jj,jk)
                  END DO
                  spgv(ji,jj) = spgv(ji,jj) * e1v(ji,jj)
                  zhinv = (tabres(ji,jj)-spgv(ji,jj))/(hv(ji,jj)*e1v(ji,jj))
                  DO jk=1,jpk             
                     vn(ji,jj,jk) = vn(ji,jj,jk) + zhinv
                     vn(ji,jj,jk) = vn(ji,jj,jk) * vmask(ji,jj,jk)
                  END DO
               ENDIF
            END DO
         END DO
      ENDIF

   END SUBROUTINE updatev2d

   SUBROUTINE updateSSH( tabres, i1, i2, j1, j2, before )
      !!---------------------------------------------
      !!          *** ROUTINE updateSSH ***
      !!---------------------------------------------
#  include "domzgr_substitute.h90"

      INTEGER, INTENT(in) :: i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before

      INTEGER :: ji, jj
      REAL(wp) :: zrhox, zrhoy

      IF (before) THEN
         zrhox = Agrif_Rhox()
         zrhoy = Agrif_Rhoy()
         DO jj=j1,j2
            DO ji=i1,i2
               tabres(ji,jj) = e1t(ji,jj) * e2t(ji,jj) * sshn(ji,jj)
            END DO
         END DO
         tabres = zrhox * zrhoy * tabres
      ELSE
         DO jj=j1,j2
            DO ji=i1,i2
               sshn(ji,jj) = tabres(ji,jj) / (e1t(ji,jj) * e2t(ji,jj))
               sshn(ji,jj) = sshn(ji,jj) * tmask(ji,jj,1)
            END DO
         END DO
      ENDIF

   END SUBROUTINE updateSSH

#else
CONTAINS
   SUBROUTINE agrif_opa_update_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_opa_update_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_opa_update : You should not have seen this print! error?'
   END SUBROUTINE agrif_opa_update_empty
#endif
END MODULE agrif_opa_update
