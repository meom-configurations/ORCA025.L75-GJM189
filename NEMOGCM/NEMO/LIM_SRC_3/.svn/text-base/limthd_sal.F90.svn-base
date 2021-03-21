MODULE limthd_sal
   !!======================================================================
   !!                       ***  MODULE limthd_sal ***
   !! LIM-3 sea-ice :  computation of salinity variations in the ice
   !!======================================================================
   !! History :   -   ! 2003-05 (M. Vancoppenolle) UCL-ASTR first coding for LIM3-1D
   !!            3.0  ! 2005-12 (M. Vancoppenolle) adapted to the 3-D version
   !!            4.0  ! 2011-02 (G. Madec) dynamical allocation
   !!---------------------------------------------------------------------
#if defined key_lim3
   !!----------------------------------------------------------------------
   !!   'key_lim3'                                      LIM-3 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_thd_sal   : salinity variations in the ice
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE phycst         ! physical constants (ocean directory)
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE ice            ! LIM variables
   USE par_ice        ! LIM parameters
   USE thd_ice        ! LIM thermodynamics
   USE limvar         ! LIM variables
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! work arrays
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_thd_sal        ! called by limthd module
   PUBLIC   lim_thd_sal_init   ! called by iceini module

   !!----------------------------------------------------------------------
   !! NEMO/LIM3 3.4 , UCL - NEMO Consortium (2011)
   !! $Id: limthd_sal.F90 3625 2012-11-21 13:19:18Z acc $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_thd_sal( kideb, kiut )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE lim_thd_sal  ***    
      !!   
      !! ** Purpose :   computes new salinities in the ice
      !!
      !! ** Method  :  3 possibilities
      !!               -> num_sal = 1 -> Sice = cst    [ice salinity constant in both time & space] 
      !!               -> num_sal = 2 -> Sice = S(z,t) [Vancoppenolle et al. 2005]
      !!               -> num_sal = 3 -> Sice = S(z)   [multiyear ice]
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kideb, kiut   ! thickness category index
      !
      INTEGER  ::   ji, jk     ! dummy loop indices 
      REAL(wp) ::   zsold, iflush, iaccrbo, igravdr, isnowic, i_ice_switch,  ztmelts   ! local scalars
      REAL(wp) ::   zaaa, zbbb, zccc, zdiscrim   ! local scalars
      REAL(wp), POINTER, DIMENSION(:) ::   ze_init, zhiold, zsiold
      !!---------------------------------------------------------------------

      CALL wrk_alloc( jpij, ze_init, zhiold, zsiold )

      !------------------------------------------------------------------------------|
      ! 1) Constant salinity, constant in time                                       |
      !------------------------------------------------------------------------------|
!!gm comment: if num_sal = 1 s_i_new, s_i_b and sm_i_b can be set to bulk_sal one for all in the initialisation phase !!
!!gm           ===>>>   simplification of almost all test on num_sal value
      IF(  num_sal == 1  ) THEN
            s_i_b  (kideb:kiut,1:nlay_i) =  bulk_sal
            sm_i_b (kideb:kiut)          =  bulk_sal 
            s_i_new(kideb:kiut)          =  bulk_sal
      ENDIF

      !------------------------------------------------------------------------------|
      !  Module 2 : Constant salinity varying in time                                |
      !------------------------------------------------------------------------------|

      IF(  num_sal == 2  ) THEN

         !---------------------------------
         ! Thickness at previous time step
         !---------------------------------
         DO ji = kideb, kiut
            zhiold(ji) = ht_i_b(ji) - dh_i_bott(ji) - dh_snowice(ji) - dh_i_surf(ji)
         END DO

         !---------------------
         ! Global heat content
         !---------------------
         ze_init(:)  =  0._wp
         DO jk = 1, nlay_i
            DO ji = kideb, kiut
               ze_init(ji) = ze_init(ji) + q_i_b(ji,jk) * ht_i_b(ji) / nlay_i
            END DO
         END DO

         DO ji = kideb, kiut
            !
            ! Switches 
            !----------
            iflush       =         MAX( 0._wp , SIGN( 1.0 , t_su_b(ji) - rtt )        )    ! =1 if summer 
            igravdr      =         MAX( 0._wp , SIGN( 1.0 , t_bo_b(ji) - t_su_b(ji) ) )    ! =1 if t_su < t_bo
            iaccrbo      =         MAX( 0._wp , SIGN( 1.0 , dh_i_bott(ji) )           )    ! =1 if bottom accretion
            i_ice_switch = 1._wp - MAX ( 0._wp , SIGN( 1._wp , - ht_i_b(ji) + 1.e-2 ) )
            isnowic      = 1._wp - MAX ( 0._wp , SIGN( 1._wp , - dh_snowice(ji) ) ) * i_ice_switch   ! =1 if snow ice formation

            !---------------------
            ! Salinity tendencies
            !---------------------
            !                                   ! drainage by gravity drainage
            dsm_i_gd_1d(ji) = - igravdr * MAX( sm_i_b(ji) - sal_G , 0._wp ) / time_G * rdt_ice 
            !                                   ! drainage by flushing  
            dsm_i_fl_1d(ji) = - iflush  * MAX( sm_i_b(ji) - sal_F , 0._wp ) / time_F * rdt_ice

            !-----------------
            ! Update salinity   
            !-----------------
            ! only drainage terms ( gravity drainage and flushing )
            ! snow ice / bottom sources are added in lim_thd_ent to conserve energy
            zsiold(ji) = sm_i_b(ji)
            sm_i_b(ji) = sm_i_b(ji) + dsm_i_fl_1d(ji) + dsm_i_gd_1d(ji)

            ! if no ice, salinity = 0.1
            i_ice_switch = 1._wp - MAX ( 0._wp, SIGN( 1._wp , - ht_i_b(ji) ) )
            sm_i_b(ji)   = i_ice_switch * sm_i_b(ji) + s_i_min * ( 1._wp - i_ice_switch )
         END DO ! ji

         CALL lim_var_salprof1d( kideb, kiut )         ! Salinity profile


         !----------------------------
         ! Heat flux - brine drainage
         !----------------------------

         DO ji = kideb, kiut
!!gm useless
            ! iflush  : 1 if summer 
            iflush  =  MAX( 0._wp , SIGN( 1._wp , t_su_b(ji) - rtt        )  ) 
            ! igravdr : 1 if t_su lt t_bo
            igravdr =  MAX( 0._wp , SIGN( 1._wp , t_bo_b(ji) - t_su_b(ji) )  ) 
            ! iaccrbo : 1 if bottom accretion
            iaccrbo =  MAX( 0._wp , SIGN( 1._wp , dh_i_bott(ji)           )  )
!!gm end useless
            !
            fhbri_1d(ji) = 0._wp
         END DO ! ji

         !----------------------------
         ! Salt flux - brine drainage
         !----------------------------
         DO ji = kideb, kiut
            i_ice_switch = 1._wp - MAX(  0._wp, SIGN( 1._wp , - ht_i_b(ji) )  )
            sfx_bri_1d(ji) = sfx_bri_1d(ji) - i_ice_switch * rhoic * a_i_b(ji) * ht_i_b(ji)         &
               &           * ( MAX( dsm_i_gd_1d(ji) + dsm_i_fl_1d(ji) , sm_i_b(ji) - zsiold(ji) ) ) * r1_rdtice
         END DO

         ! Only necessary for conservation check since salinity is modified
         !--------------------
         ! Temperature update
         !--------------------
         DO jk = 1, nlay_i
            DO ji = kideb, kiut
               ztmelts    =  -tmut*s_i_b(ji,jk) + rtt
               !Conversion q(S,T) -> T (second order equation)
               zaaa         =  cpic
               zbbb         =  ( rcp - cpic ) * ( ztmelts - rtt ) + q_i_b(ji,jk) / rhoic - lfus
               zccc         =  lfus * ( ztmelts - rtt )
               zdiscrim     =  SQRT(  MAX( zbbb*zbbb - 4.0*zaaa*zccc, 0._wp )  )
               t_i_b(ji,jk) =  rtt - ( zbbb + zdiscrim ) / ( 2.0 *zaaa )
            END DO
         END DO
         !
      ENDIF 

      !------------------------------------------------------------------------------|
      !  Module 3 : Profile of salinity, constant in time                            |
      !------------------------------------------------------------------------------|

      IF(  num_sal == 3  )   CALL lim_var_salprof1d( kideb, kiut )


      !------------------------------------------------------------------------------|
      ! 5) Computation of salt flux due to Bottom growth
      !------------------------------------------------------------------------------|
      ! note: s_i_new = bulk_sal in constant salinity case
      DO ji = kideb, kiut
         sfx_thd_1d(ji) = sfx_thd_1d(ji) - s_i_new(ji) * rhoic * a_i_b(ji) * MAX( dh_i_bott(ji) , 0._wp ) * r1_rdtice
      END DO
      !
      CALL wrk_dealloc( jpij, ze_init, zhiold, zsiold )
      !
   END SUBROUTINE lim_thd_sal


   SUBROUTINE lim_thd_sal_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE lim_thd_sal_init  ***
      !!
      !! ** Purpose :   initialization of ice salinity parameters
      !!
      !! ** Method  :   Read the namicesal namelist and check the parameter
      !!              values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namicesal
      !!-------------------------------------------------------------------
      NAMELIST/namicesal/ num_sal, bulk_sal, sal_G, time_G, sal_F, time_F,   &
         &                s_i_max, s_i_min, s_i_0, s_i_1
      !!-------------------------------------------------------------------
      !
      REWIND( numnam_ice )                   ! Read Namelist namicesal
      READ  ( numnam_ice  , namicesal )
      !
      IF(lwp) THEN                           ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'lim_thd_sal_init : Ice parameters for salinity '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) ' switch for salinity num_sal        : ', num_sal
         WRITE(numout,*) ' bulk salinity value if num_sal = 1 : ', bulk_sal
         WRITE(numout,*) ' restoring salinity for GD          : ', sal_G
         WRITE(numout,*) ' restoring time for GD              : ', time_G
         WRITE(numout,*) ' restoring salinity for flushing    : ', sal_F
         WRITE(numout,*) ' restoring time for flushing        : ', time_F
         WRITE(numout,*) ' Maximum tolerated ice salinity     : ', s_i_max
         WRITE(numout,*) ' Minimum tolerated ice salinity     : ', s_i_min
         WRITE(numout,*) ' 1st salinity for salinity profile  : ', s_i_0
         WRITE(numout,*) ' 2nd salinity for salinity profile  : ', s_i_1
      ENDIF
      !
   END SUBROUTINE lim_thd_sal_init

#else
   !!----------------------------------------------------------------------
   !!   Default option         Dummy Module          No LIM-3 sea-ice model
   !!----------------------------------------------------------------------
#endif
   !!======================================================================
END MODULE limthd_sal
