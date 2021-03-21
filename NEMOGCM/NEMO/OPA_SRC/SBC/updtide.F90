MODULE updtide
  !!=================================================================================
  !!                       ***  MODULE  updtide  ***
  !! Initialization of tidal forcing
  !! History :  9.0  !  07  (O. Le Galloudec)  Original code
  !!=================================================================================
#if defined key_tide
  !! * Modules used
  USE oce             ! ocean dynamics and tracers variables
  USE dom_oce         ! ocean space and time domain
  USE in_out_manager  ! I/O units
  USE phycst
  USE sbctide
  USE dynspg_oce
  USE tideini, ONLY: ln_tide_ramp, rdttideramp

  IMPLICIT NONE
  PUBLIC

  !! * Routine accessibility
  PUBLIC upd_tide
  !!---------------------------------------------------------------------------------
  !!   OPA 9.0 , LODYC-IPSL  (2003)
  !!---------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE upd_tide (kt,kit)
    !!----------------------------------------------------------------------
    !!                 ***  ROUTINE upd_tide  ***
    !!----------------------------------------------------------------------      
    !! * Local declarations

    INTEGER, INTENT( in ) ::   kt,kit      ! ocean time-step index
    INTEGER  :: ji,jj,jk
    REAL (wp) :: zramp
    REAL (wp), DIMENSION(nb_harmo) :: zwt 
    !...............................................................................

    pot_astro(:,:)=0.e0
    zramp = 1.e0

    IF (lk_dynspg_ts) THEN
       zwt(:) = omega_tide(:)* ((kt-kt_tide)*rdt + kit*(rdt/REAL(nn_baro,wp)))
       IF (ln_tide_ramp) THEN
          zramp = MIN(MAX( ((kt-nit000)*rdt + kit*(rdt/REAL(nn_baro,wp)))/(rdttideramp*rday),0.),1.)
       ENDIF
    ELSE
       zwt(:) = omega_tide(:)*(kt-kt_tide)*rdt
       IF (ln_tide_ramp) THEN
          zramp = MIN(MAX( ((kt-nit000)*rdt)/(rdttideramp*rday),0.),1.) 
       ENDIF  
    ENDIF

    do jk=1,nb_harmo
       do ji=1,jpi
          do jj=1,jpj
             pot_astro(ji,jj)=pot_astro(ji,jj) + zramp*(amp_pot(ji,jj,jk)*COS(zwt(jk)+phi_pot(ji,jj,jk)))      
          enddo
       enddo
    enddo

  END SUBROUTINE upd_tide

#else
  !!----------------------------------------------------------------------
  !!   Dummy module :                                        NO TIDE
  !!----------------------------------------------------------------------
CONTAINS
  SUBROUTINE upd_tide( kt,kit )          ! Empty routine
    INTEGER,INTENT (IN) :: kt, kit
    WRITE(*,*) 'upd_tide: You should not have seen this print! error?', kt
  END SUBROUTINE upd_tide

#endif

  !!======================================================================

END MODULE updtide
