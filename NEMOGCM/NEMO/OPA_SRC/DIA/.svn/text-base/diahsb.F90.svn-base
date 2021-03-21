MODULE diahsb
   !!======================================================================
   !!                       ***  MODULE  diahsb  ***
   !! Ocean diagnostics: Heat, salt and volume budgets
   !!======================================================================
   !! History :  3.3  ! 2010-09  (M. Leclair)  Original code 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! surface thermohaline fluxes
   USE in_out_manager  ! I/O manager
   USE domvvl          ! vertical scale factors
   USE traqsr          ! penetrative solar radiation
   USE trabbc          ! bottom boundary condition 
   USE lib_mpp         ! distributed memory computing library
   USE trabbc          ! bottom boundary condition
   USE obc_par         ! (for lk_obc)
   USE bdy_par         ! (for lk_bdy)
   USE timing          ! preformance summary
   USE lib_fortran
   USE sbcrnf

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_hsb        ! routine called by step.F90
   PUBLIC   dia_hsb_init   ! routine called by opa.F90

   LOGICAL, PUBLIC ::   ln_diahsb  = .FALSE.   !: check the heat and salt budgets

   INTEGER                                 ::   numhsb                           !
   REAL(dp)                                ::   surf_tot   , vol_tot             !
   REAL(dp)                                ::   frc_t      , frc_s     , frc_v   ! global forcing trends
   REAL(dp)                                ::   frc_wn_t      , frc_wn_s ! global forcing trends
   REAL(dp)                                ::   fact1                            ! conversion factors
   REAL(dp)                                ::   fact21    , fact22               !     -         -
   REAL(dp)                                ::   fact31    , fact32               !     -         -
   REAL(dp), DIMENSION(:,:)  , ALLOCATABLE ::   surf      , ssh_ini              !
   REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::   hc_loc_ini, sc_loc_ini, e3t_ini  !
   REAL(dp), DIMENSION(:,:)  , ALLOCATABLE ::   ssh_hc_loc_ini, ssh_sc_loc_ini

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: diahsb.F90 4071 2013-10-17 12:35:17Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dia_hsb( kt )
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hsb  ***
      !!     
      !! ** Purpose: Compute the ocean global heat content, salt content and volume conservation
      !!	
      !! ** Method : - Compute the deviation of heat content, salt content and volume
      !!	            at the current time step from their values at nit000
      !!	            - Compute the contribution of forcing and remove it from these deviations
      !!
      !! ** Action : Write the results in the 'heat_salt_volume_budgets.txt' ASCII file
      !!---------------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER    ::   jk                          ! dummy loop indice
      REAL(dp)   ::   zdiff_hc    , zdiff_sc      ! heat and salt content variations
      REAL(dp)   ::   zdiff_hc1   , zdiff_sc1     ! heat and salt content variations of ssh
      REAL(dp)   ::   zdiff_v1    , zdiff_v2      ! volume variation
      REAL(dp)   ::   zerr_hc1    , zerr_sc1      ! Non conservation due to free surface
      REAL(dp)   ::   z1_rau0                     ! local scalars
      REAL(dp)   ::   zdeltat                     !    -     -
      REAL(dp)   ::   z_frc_trd_t , z_frc_trd_s   !    -     -
      REAL(dp)   ::   z_frc_trd_v                 !    -     -
      REAL(dp)   ::   z_wn_trd_t , z_wn_trd_s   !    -     -
      REAL(dp)   ::   z_ssh_hc , z_ssh_sc   !    -     -
      !!---------------------------------------------------------------------------
      IF( nn_timing == 1 )   CALL timing_start('dia_hsb')

      ! ------------------------- !
      ! 1 - Trends due to forcing !
      ! ------------------------- !
      z1_rau0 = 1.e0 / rau0
      z_frc_trd_v = z1_rau0 * glob_sum( - ( emp(:,:) - rnf(:,:) ) * surf(:,:) )     ! volume fluxes
      z_frc_trd_t =           glob_sum( sbc_tsc(:,:,jp_tem) * surf(:,:) )     ! heat fluxes
      z_frc_trd_s =           glob_sum( sbc_tsc(:,:,jp_sal) * surf(:,:) )     ! salt fluxes
      ! Add runoff heat & salt input
      IF( ln_rnf    )   z_frc_trd_t = z_frc_trd_t + glob_sum( rnf_tsc(:,:,jp_tem) * surf(:,:) )
      IF( ln_rnf_sal)   z_frc_trd_s = z_frc_trd_s + glob_sum( rnf_tsc(:,:,jp_sal) * surf(:,:) )
      ! Add penetrative solar radiation
      IF( ln_traqsr )   z_frc_trd_t = z_frc_trd_t + r1_rau0_rcp * glob_sum( qsr     (:,:) * surf(:,:) )
      ! Add geothermal heat flux
      IF( ln_trabbc )   z_frc_trd_t = z_frc_trd_t +  glob_sum( qgh_trd0(:,:) * surf(:,:) )
      IF( .NOT. lk_vvl ) THEN
         z_wn_trd_t = - glob_sum( surf(:,:) * wn(:,:,1) * tsb(:,:,1,jp_tem) )
         z_wn_trd_s = - glob_sum( surf(:,:) * wn(:,:,1) * tsb(:,:,1,jp_sal) )
      ENDIF

      frc_v = frc_v + z_frc_trd_v * rdt
      frc_t = frc_t + z_frc_trd_t * rdt
      frc_s = frc_s + z_frc_trd_s * rdt
      !                                          ! Advection flux through fixed surface (z=0)
      IF( .NOT. lk_vvl ) THEN
         frc_wn_t = frc_wn_t + z_wn_trd_t * rdt
         frc_wn_s = frc_wn_s + z_wn_trd_s * rdt
      ENDIF

      ! ----------------------- !
      ! 2 -  Content variations !
      ! ----------------------- !
      zdiff_v2 = 0.d0
      zdiff_hc = 0.d0
      zdiff_sc = 0.d0

      ! volume variation (calculated with ssh)
      zdiff_v1 = glob_sum( surf(:,:) * ( sshn(:,:) - ssh_ini(:,:) ) )

      ! heat & salt content variation (associated with ssh)
      IF( .NOT. lk_vvl ) THEN
         z_ssh_hc = glob_sum( surf(:,:) * ( tsn(:,:,1,jp_tem) * sshn(:,:) - ssh_hc_loc_ini(:,:) ) )
         z_ssh_sc = glob_sum( surf(:,:) * ( tsn(:,:,1,jp_sal) * sshn(:,:) - ssh_sc_loc_ini(:,:) ) )
      ENDIF

      DO jk = 1, jpkm1
        ! volume variation (calculated with scale factors)
         zdiff_v2 = zdiff_v2 + glob_sum( surf(:,:) * tmask(:,:,jk)   &
            &                       * ( fse3t_n(:,:,jk)         &
            &                           - e3t_ini(:,:,jk) ) )
         ! heat content variation
         zdiff_hc = zdiff_hc + glob_sum( surf(:,:) * tmask(:,:,jk)          &
            &                       * ( fse3t_n(:,:,jk) * tsn(:,:,jk,jp_tem)   &
            &                           - hc_loc_ini(:,:,jk) ) )
         ! salt content variation
         zdiff_sc = zdiff_sc + glob_sum( surf(:,:) * tmask(:,:,jk)          &
            &                       * ( fse3t_n(:,:,jk) * tsn(:,:,jk,jp_sal)   &
            &                           - sc_loc_ini(:,:,jk) ) )
      ENDDO

      ! Substract forcing from heat content, salt content and volume variations
      zdiff_v1 = zdiff_v1 - frc_v
      IF( lk_vvl )   zdiff_v2 = zdiff_v2 - frc_v
      zdiff_hc = zdiff_hc - frc_t
      zdiff_sc = zdiff_sc - frc_s
      IF( .NOT. lk_vvl ) THEN
         zdiff_hc1 = zdiff_hc + z_ssh_hc 
         zdiff_sc1 = zdiff_sc + z_ssh_sc
         zerr_hc1  = z_ssh_hc - frc_wn_t
         zerr_sc1  = z_ssh_sc - frc_wn_s
      ENDIF
      
      ! ----------------------- !
      ! 3 - Diagnostics writing !
      ! ----------------------- !
      zdeltat  = 1.e0 / ( ( kt - nit000 + 1 ) * rdt )
      IF( lk_vvl ) THEN
         WRITE(numhsb , 9020) kt , zdiff_hc / vol_tot , zdiff_hc * fact1  * zdeltat,                                &
            &                      zdiff_sc / vol_tot , zdiff_sc * fact21 * zdeltat, zdiff_sc * fact22 * zdeltat,   &
            &                      zdiff_v1           , zdiff_v1 * fact31 * zdeltat, zdiff_v1 * fact32 * zdeltat,   &
            &                      zdiff_v2           , zdiff_v2 * fact31 * zdeltat, zdiff_v2 * fact32 * zdeltat
      ELSE
         WRITE(numhsb , 9030) kt , zdiff_hc1 / vol_tot , zdiff_hc1 * fact1  * zdeltat,                                &
            &                      zdiff_sc1 / vol_tot , zdiff_sc1 * fact21 * zdeltat, zdiff_sc1 * fact22 * zdeltat,   &
            &                      zdiff_v1            , zdiff_v1  * fact31 * zdeltat, zdiff_v1  * fact32 * zdeltat,   &
            &                      zerr_hc1 / vol_tot  , zerr_sc1 / vol_tot
      ENDIF

      IF ( kt == nitend ) CLOSE( numhsb )

      IF( nn_timing == 1 )   CALL timing_stop('dia_hsb')

9020  FORMAT(I5,11D15.7)
9030  FORMAT(I5,10D15.7)
      !
   END SUBROUTINE dia_hsb


   SUBROUTINE dia_hsb_init
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hsb  ***
      !!     
      !! ** Purpose: Initialization for the heat salt volume budgets
      !!	
      !! ** Method : Compute initial heat content, salt content and volume
      !!
      !! ** Action : - Compute initial heat content, salt content and volume
      !!             - Initialize forcing trends
      !!             - Compute coefficients for conversion
      !!---------------------------------------------------------------------------
      CHARACTER (len=32) ::   cl_name  ! output file name
      INTEGER            ::   jk       ! dummy loop indice
      INTEGER            ::   ierror   ! local integer
      !!
      NAMELIST/namhsb/ ln_diahsb
      !!----------------------------------------------------------------------
      !
      REWIND ( numnam )              ! Read Namelist namhsb 
      READ   ( numnam, namhsb )
      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dia_hsb_init : check the heat and salt budgets'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namhsb : set hsb parameters'
         WRITE(numout,*) '      Switch for hsb diagnostic (T) or not (F)  ln_diahsb  = ', ln_diahsb
      ENDIF

      IF( .NOT. ln_diahsb )   RETURN
      IF( .NOT. lk_mpp_rep ) &
        CALL ctl_stop (' Your global mpp_sum if performed in single precision - 64 bits -', &
             &         ' whereas the global sum to be precise must be done in double precision ',&
             &         ' please add key_mpp_rep')

      ! ------------------- !
      ! 1 - Allocate memory !
      ! ------------------- !
      ALLOCATE( hc_loc_ini(jpi,jpj,jpk), sc_loc_ini(jpi,jpj,jpk), &
         &      ssh_hc_loc_ini(jpi,jpj), ssh_sc_loc_ini(jpi,jpj), &
         &      e3t_ini(jpi,jpj,jpk)                            , &
         &      surf(jpi,jpj),  ssh_ini(jpi,jpj), STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_hsb: unable to allocate hc_loc_ini' )   ;   RETURN
      ENDIF

      ! ----------------------------------------------- !
      ! 2 - Time independant variables and file opening !
      ! ----------------------------------------------- !
      WRITE(numout,*) "dia_hsb: heat salt volume budgets activated"
      WRITE(numout,*) "~~~~~~~  output written in the 'heat_salt_volume_budgets.txt' ASCII file"
      IF( lk_obc .or. lk_bdy ) THEN
         CALL ctl_warn( 'dia_hsb does not take open boundary fluxes into account' )         
      ENDIF
      cl_name    = 'heat_salt_volume_budgets.txt'                         ! name of output file
      surf(:,:) = e1t(:,:) * e2t(:,:) * tmask(:,:,1) * tmask_i(:,:)      ! masked surface grid cell area
      surf_tot  = glob_sum( surf(:,:) )                                       ! total ocean surface area
      vol_tot   = 0.d0                                                   ! total ocean volume
      DO jk = 1, jpkm1
         vol_tot  = vol_tot + glob_sum( surf(:,:) * tmask(:,:,jk)     &
            &                         * fse3t_n(:,:,jk)         )
      END DO

      CALL ctl_opn( numhsb , cl_name , 'UNKNOWN' , 'FORMATTED' , 'SEQUENTIAL' , 1 , numout , lwp , 1 )
      IF( lk_vvl ) THEN
         !                   12345678901234567890123456789012345678901234567890123456789012345678901234567890 -> 80
         WRITE( numhsb, 9010 ) "kt   |     heat content budget     |            salt content budget             ",   &
            !                                                   123456789012345678901234567890123456789012345 -> 45
            &                                                  "|            volume budget (ssh)             ",   &
            !                                                   678901234567890123456789012345678901234567890 -> 45
            &                                                  "|            volume budget (e3t)             "
         WRITE( numhsb, 9010 ) "     |      [C]         [W/m2]     |     [psu]        [mmm/s]          [SV]     ",   &
            &                                                  "|     [m3]         [mmm/s]          [SV]     ",   &
            &                                                  "|     [m3]         [mmm/s]          [SV]     "
      ELSE
         !                   12345678901234567890123456789012345678901234567890123456789012345678901234567890 -> 80
         WRITE( numhsb, 9011 ) "kt   |     heat content budget     |            salt content budget             ",   &
            !                                                   123456789012345678901234567890123456789012345 -> 45
            &                                                  "|            volume budget (ssh)             ",   &
            !                                                   678901234567890123456789012345678901234567890 -> 45
            &                                                  "|  Non conservation due to free surface      "
         WRITE( numhsb, 9011 ) "     |      [C]         [W/m2]     |     [psu]        [mmm/s]          [SV]     ",   &
            &                                                  "|     [m3]         [mmm/s]          [SV]     ",   &
            &                                                  "|  [heat - C]     [salt - psu]                "
      ENDIF
      ! --------------- !
      ! 3 - Conversions ! (factors will be multiplied by duration afterwards)
      ! --------------- !

      ! heat content variation   =>   equivalent heat flux:
      fact1  = rau0 * rcp / surf_tot                                         ! [C*m3]   ->  [W/m2]
      ! salt content variation   =>   equivalent EMP and equivalent "flow": 
      fact21 = 1.e3  / ( soce * surf_tot )                                   ! [psu*m3] ->  [mm/s]
      fact22 = 1.e-6 / soce                                                  ! [psu*m3] ->  [Sv]
      ! volume variation         =>   equivalent EMP and equivalent "flow":
      fact31 = 1.e3  / surf_tot                                              ! [m3]     ->  [mm/s]
      fact32 = 1.e-6                                                         ! [m3]     ->  [SV]

      ! ---------------------------------- !
      ! 4 - initial conservation variables !
      ! ---------------------------------- !
      ssh_ini(:,:) = sshn(:,:)                                       ! initial ssh
      DO jk = 1, jpk
         e3t_ini   (:,:,jk) = fse3t_n(:,:,jk)                        ! initial vertical scale factors
         hc_loc_ini(:,:,jk) = tsn(:,:,jk,jp_tem) * fse3t_n(:,:,jk)   ! initial heat content
         sc_loc_ini(:,:,jk) = tsn(:,:,jk,jp_sal) * fse3t_n(:,:,jk)   ! initial salt content
      END DO
      frc_v = 0.d0                                           ! volume       trend due to forcing
      frc_t = 0.d0                                           ! heat content   -    -   -    -   
      frc_s = 0.d0                                           ! salt content   -    -   -    -         
      IF( .NOT. lk_vvl ) THEN
         ssh_hc_loc_ini(:,:) = tsn(:,:,1,jp_tem) * ssh_ini(:,:)   ! initial heat content associated with ssh
         ssh_sc_loc_ini(:,:) = tsn(:,:,1,jp_sal) * ssh_ini(:,:)   ! initial salt content associated with ssh
         frc_wn_t = 0.d0
         frc_wn_s = 0.d0
      ENDIF
      !
9010  FORMAT(A80,A45,A45)
9011  FORMAT(A80,A45,A45)
      !
   END SUBROUTINE dia_hsb_init

   !!======================================================================
END MODULE diahsb
