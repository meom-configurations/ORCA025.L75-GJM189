MODULE p2zsed
   !!======================================================================
   !!                         ***  MODULE p2zsed  ***
   !! TOP :   PISCES Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :    -   !  1995-06 (M. Levy)  original code
   !!              -   !  2000-12 (E. Kestenare)  clean up
   !!             2.0  !  2007-12  (C. Deltel, G. Madec)  F90 + simplifications
   !!----------------------------------------------------------------------
#if defined key_pisces_reduced
   !!----------------------------------------------------------------------
   !!   'key_pisces_reduced'                                     LOBSTER bio-model
   !!----------------------------------------------------------------------
   !!   p2z_sed        :  Compute loss of organic matter in the sediments
   !!----------------------------------------------------------------------
   USE oce_trc         !
   USE trc
   USE sms_pisces
   USE lbclnk
   USE trdmod_oce
   USE trdmod_trc
   USE iom
   USE prtctl_trc      ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p2z_sed         ! called in ???
   PUBLIC   p2z_sed_init    ! called in ???

   REAL(wp), PUBLIC ::   sedlam     = 3.86e-7   !: time coefficient of POC remineralization in sediments
   REAL(wp), PUBLIC ::   sedlostpoc = 0.        ! mass of POC lost in sediments 
   REAL(wp), PUBLIC ::   vsed       = 3.47e-5   ! detritus sedimentation speed [m/s] 
   REAL(wp), PUBLIC ::   xhr        = -0.858    ! coeff for martin''s remineralisation profile

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p2z_sed.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p2z_sed( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_sed  ***
      !!
      !! ** Purpose :   compute the now trend due to the vertical sedimentation of
      !!              detritus and add it to the general trend of detritus equations
      !!
      !! ** Method  :   this ROUTINE compute not exactly the advection but the
      !!              transport term, i.e.  dz(wt) and dz(ws)., dz(wtr)
      !!              using an upstream scheme
      !!              the now vertical advection of tracers is given by:
      !!                      dz(trn wn) = 1/bt dk+1( e1t e2t vsed (trn) )
      !!              add this trend now to the general trend of tracer (ta,sa,tra):
      !!                             tra = tra + dz(trn wn)
      !!        
      !!              IF 'key_diabio' is defined, the now vertical advection
      !!              trend of passive tracers is saved for futher diagnostics.
      !!---------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      INTEGER  ::   ji, jj, jk, jl, ierr
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zw2d
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zwork, ztra, ztrbio
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p2z_sed')
      !
      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' p2z_sed: LOBSTER sedimentation'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~'
      ENDIF

      ! Allocate temporary workspace
                       CALL wrk_alloc( jpi, jpj, jpk, zwork, ztra )
      IF( ln_diatrc )  CALL wrk_alloc( jpi, jpj, zw2d )
      IF( l_trdtrc ) THEN
         CALL wrk_alloc( jpi, jpj, jpk, ztrbio )
         ztrbio(:,:,:) = tra(:,:,:,jpdet)
      ENDIF

      ! sedimentation of detritus  : upstream scheme
      ! --------------------------------------------

      ! for detritus sedimentation only - jpdet
      zwork(:,:,1  ) = 0.e0      ! surface value set to zero
      zwork(:,:,jpk) = 0.e0      ! bottom value  set to zero

      ! tracer flux at w-point: we use -vsed (downward flux)  with simplification : no e1*e2
      DO jk = 2, jpkm1
         zwork(:,:,jk) = -vsed * trn(:,:,jk-1,jpdet)
      END DO

      ! tracer flux divergence at t-point added to the general trend
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztra(ji,jj,jk)  = - ( zwork(ji,jj,jk) - zwork(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
               tra(ji,jj,jk,jpdet) = tra(ji,jj,jk,jpdet) + ztra(ji,jj,jk) 
            END DO
         END DO
      END DO

      IF( ln_diatrc ) THEN 
         DO jk = 1, jpkm1
            zw2d(:,:) = zw2d(:,:) + ztra(:,:,jk) * fse3t(:,:,jk) * 86400.
         END DO
         IF( lk_iomput )  THEN
           CALL iom_put( "TDETSED", zw2d )
         ELSE
           trc2d(:,:,jp_pcs0_2d + 7) = zw2d(:,:)
         ENDIF
         CALL wrk_dealloc( jpi, jpj, zw2d )
      ENDIF
      !
      IF( ln_diabio )  trbio(:,:,:,jp_pcs0_trd + 7) = ztra(:,:,:)
      CALL wrk_dealloc( jpi, jpj, jpk, zwork, ztra )
      !
      IF( l_trdtrc ) THEN
         ztrbio(:,:,:) = tra(:,:,:,jpdet) - ztrbio(:,:,:)
         jl = jp_pcs0_trd + 7
         CALL trd_mod_trc( ztrbio, jl, kt )   ! handle the trend
         CALL wrk_dealloc( jpi, jpj, jpk, ztrbio )
      ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sed')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p2z_sed')
      !
   END SUBROUTINE p2z_sed

   SUBROUTINE p2z_sed_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p2z_sed_init  ***
      !!
      !! ** Purpose :   Parameters from aphotic layers to sediment
      !!
      !! ** Method  :   Read the namlobsed namelist and check the parameters
      !!
      !!----------------------------------------------------------------------
      NAMELIST/namlobsed/ sedlam, sedlostpoc, vsed, xhr

      REWIND( numnatp )
      READ  ( numnatp, namlobsed )

      IF(lwp) THEN
          WRITE(numout,*) ' Namelist namlobsed'
          WRITE(numout,*) '    time coeff of POC in sediments                       sedlam    =', sedlam
          WRITE(numout,*) '    Sediment geol loss for POC                           sedlostpoc=', sedlostpoc
          WRITE(numout,*) '    detritus sedimentation speed                         vsed      =', 86400 * vsed  , ' d'
          WRITE(numout,*) '    coeff for martin''s remineralistion                  xhr       =', xhr
          WRITE(numout,*) ' '
      ENDIF
      !
   END SUBROUTINE p2z_sed_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p2z_sed( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p2z_sed: You should not have seen this print! error?', kt
   END SUBROUTINE p2z_sed
#endif 

   !!======================================================================
END MODULE  p2zsed
