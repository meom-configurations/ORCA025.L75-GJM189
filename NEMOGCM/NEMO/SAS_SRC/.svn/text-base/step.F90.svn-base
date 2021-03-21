MODULE step
   !!======================================================================
   !!                       ***  MODULE step  ***
   !! Time-stepping    : manager of the ocean, tracer and ice time stepping
   !!                    version for standalone surface scheme
   !!======================================================================
   !! History :  OPA  !  1991-03  (G. Madec)  Original code
   !!             .   !    .                                                     
   !!             .   !    .                                                     
   !!   NEMO     3.5  !  2012-03  (S. Alderson)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   stp             : OPA system time-stepping
   !!----------------------------------------------------------------------
   USE oce              ! ocean dynamics and tracers variables
   USE dom_oce          ! ocean space and time domain variables 
   USE in_out_manager   ! I/O manager
   USE iom              !
   USE lbclnk
#if defined key_iomput
   USE xios
#endif

   USE daymod           ! calendar                         (day     routine)

   USE sbcmod           ! surface boundary condition       (sbc     routine)
   USE sbcrnf           ! surface boundary condition: runoff variables

   USE eosbn2           ! equation of state                (eos_bn2 routine)

   USE diawri           ! Standard run outputs             (dia_wri routine)
   USE stpctl           ! time stepping control            (stp_ctl routine)
   USE prtctl           ! Print control                    (prt_ctl routine)

   USE timing           ! Timing            

   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp   ! called by opa.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "zdfddm_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: step.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_agrif
   SUBROUTINE stp( )
      INTEGER             ::   kstp   ! ocean time-step index
#else
   SUBROUTINE stp( kstp )
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index
#endif
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp  ***
      !!                      
      !! ** Purpose : - Time stepping of SBC (surface boundary)
      !! 
      !! ** Method  : -1- Update forcings and data  
      !!              -2- Outputs and diagnostics
      !!----------------------------------------------------------------------
      INTEGER ::   indic    ! error indicator if < 0
      !! ---------------------------------------------------------------------

#if defined key_agrif
      kstp = nit000 + Agrif_Nb_Step()
# if defined key_iomput
      IF( Agrif_Nbstepint() == 0 )   CALL iom_swap
# endif   
#endif   
      IF( kstp == nit000 )   CALL iom_init            ! iom_put initialization (must be done after nemo_init for AGRIF+XIOS+OASIS)
      IF( kstp /= nit000 )   CALL day( kstp )             ! Calendar (day was already called at nit000 in day_init)
                             CALL iom_setkt( kstp )       ! say to iom that we are at time step kstp

                             CALL sbc    ( kstp )         ! Sea Boundary Condition (including sea-ice)

                             CALL dia_wri( kstp )         ! ocean model: outputs

                             indic = 0                    ! although indic is not changed in stp_ctl
                                                          ! need to keep the same interface 
                             CALL stp_ctl( kstp, indic )
#if defined key_iomput
      IF( kstp == nitend )   CALL xios_context_finalize() ! needed for XIOS+AGRIF
#endif
      !
      IF( nn_timing == 1 .AND.  kstp == nit000  )   CALL timing_reset
      !
   END SUBROUTINE stp

   !!======================================================================
END MODULE step
