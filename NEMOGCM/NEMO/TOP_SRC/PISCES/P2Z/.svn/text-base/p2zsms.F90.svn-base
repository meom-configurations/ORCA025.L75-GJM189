MODULE p2zsms
   !!======================================================================
   !!                         ***  MODULE p2zsms  ***
   !! TOP :   Time loop of LOBSTER model
   !!======================================================================
   !! History :   1.0  !            M. Levy
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
#if defined key_pisces_reduced
   !!----------------------------------------------------------------------
   !!   'key_pisces_reduced'                              LOBSTER bio-model
   !!----------------------------------------------------------------------
   !!   p2zsms        :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc          !
   USE trc
   USE p2zbio
   USE p2zopt
   USE p2zsed
   USE p2zexp
   USE trdmod_oce
   USE trdmod_trc_oce
   USE trdmod_trc
   USE trdmld_trc

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p2z_sms    ! called in p2zsms.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p2zsms.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p2z_sms( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_sms  ***
      !!
      !! ** Purpose :  Managment of the call to Biological sources and sinks 
      !!               routines of LOBSTER bio-model 
      !!
      !! ** Method  : - ???
      !! --------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !
      INTEGER :: jn
      !! --------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p2z_sms')
      !
      CALL p2z_opt( kt )      ! optical model
      CALL p2z_bio( kt )      ! biological model
      CALL p2z_sed( kt )      ! sedimentation model
      CALL p2z_exp( kt )      ! export

      IF( l_trdtrc ) THEN
         DO jn = jp_pcs0, jp_pcs1
           CALL trd_mod_trc( tra(:,:,:,jn), jn, jptra_trd_sms, kt )   ! save trends
         END DO
      END IF

      IF( lk_trdmld_trc )  CALL trd_mld_bio( kt )   ! trends: Mixed-layer
      !
      IF( nn_timing == 1 )  CALL timing_stop('p2z_sms')
      !
   END SUBROUTINE p2z_sms

#else
   !!======================================================================
   !!  Dummy module :                                     No passive tracer
   !!======================================================================
CONTAINS
   SUBROUTINE p2z_sms( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p2z_sms: You should not have seen this print! error?', kt
   END SUBROUTINE p2z_sms
#endif 

   !!======================================================================
END MODULE  p2zsms
