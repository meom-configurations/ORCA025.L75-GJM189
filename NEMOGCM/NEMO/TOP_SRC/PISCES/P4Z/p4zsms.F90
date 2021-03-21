MODULE p4zsms
   !!======================================================================
   !!                         ***  MODULE p4zsms  ***
   !! TOP :   PISCES Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4zsms        :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE trcdta
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zbio          !  Biological model
   USE p4zche          !  Chemical model
   USE p4zlys          !  Calcite saturation
   USE p4zflx          !  Gas exchange
   USE p4zsbc          !  External source of nutrients
   USE p4zsed          !  Sedimentation
   USE p4zint          !  time interpolation
   USE iom             !  I/O manager
   USE trdmod_oce      !  Ocean trends variables
   USE trdmod_trc      !  TOP trends variables
   USE sedmodel        !  Sediment model
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sms_init    ! called in p4zsms.F90
   PUBLIC   p4z_sms    ! called in p4zsms.F90

   REAL(wp) :: alkbudget, no3budget, silbudget, ferbudget
   INTEGER ::  numco2, numnut  !: logical unit for co2 budget

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zsms.F90 3320 2012-03-05 16:37:52Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_sms( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sms  ***
      !!
      !! ** Purpose :   Managment of the call to Biological sources and sinks 
      !!              routines of PISCES bio-model
      !!
      !! ** Method  : - at each new day ...
      !!              - several calls of bio and sed ???
      !!              - ...
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      INTEGER ::   jnt, jn, jl
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:,:)  :: ztrdpis
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sms')
      !
      IF( l_trdtrc )  THEN
         CALL wrk_alloc( jpi, jpj, jpk, jp_pisces, ztrdpis ) 
         DO jn = 1, jp_pisces
            jl = jn + jp_pcs0 - 1
            ztrdpis(:,:,:,jn) = trn(:,:,:,jl)
         ENDDO
      ENDIF
      !
      IF( kt == nittrc000 ) THEN
        !
        CALL p4z_che                              ! initialize the chemical constants
        !
        IF( .NOT. ln_rsttr ) THEN  ;   CALL p4z_ph_ini   !  set PH at kt=nit000 
        ELSE                       ;   CALL p4z_rst( nittrc000, 'READ' )  !* read or initialize all required fields 
        ENDIF
        !
      ENDIF

      IF( ln_pisdmp .AND. MOD( kt - nn_dttrc, nn_pisdmp ) == 0 )   CALL p4z_dmp( kt )      ! Relaxation of some tracers
      !
      IF( ndayflxtr /= nday_year ) THEN      ! New days
         !
         ndayflxtr = nday_year

         IF(lwp) write(numout,*)
         IF(lwp) write(numout,*) ' New chemical constants and various rates for biogeochemistry at new day : ', nday_year
         IF(lwp) write(numout,*) '~~~~~~'

         CALL p4z_che              ! computation of chemical constants
         CALL p4z_int( kt )        ! computation of various rates for biogeochemistry
         !
      ENDIF

      IF( ll_sbc ) CALL p4z_sbc( kt )   ! external sources of nutrients 

      DO jnt = 1, nrdttrc          ! Potential time splitting if requested
         !
         CALL p4z_bio (kt, jnt)    ! Biology
         CALL p4z_sed (kt, jnt)    ! Sedimentation
         !
         DO jn = jp_pcs0, jp_pcs1
            trb(:,:,:,jn) = trn(:,:,:,jn)
         ENDDO
         !
      END DO

      IF( l_trdtrc )  THEN
         DO jn = 1, jp_pisces
            jl = jn + jp_pcs0 - 1
            ztrdpis(:,:,:,jn) = ( ztrdpis(:,:,:,jn) - trn(:,:,:,jl) ) * rfact2r
         ENDDO
      ENDIF

      CALL p4z_lys( kt )             ! Compute CaCO3 saturation
      CALL p4z_flx( kt )             ! Compute surface fluxes

      DO jn = jp_pcs0, jp_pcs1
        CALL lbc_lnk( trn(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk( trb(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk( tra(:,:,:,jn), 'T', 1. )
      END DO
      !
      IF( lk_sed ) THEN 
         !
         CALL sed_model( kt )     !  Main program of Sediment model
         !
         DO jn = jp_pcs0, jp_pcs1
           CALL lbc_lnk( trn(:,:,:,jn), 'T', 1. )
         END DO
         !
      ENDIF
      !
      IF( lrst_trc )  CALL p4z_rst( kt, 'WRITE' )  !* Write PISCES informations in restart file 
      !
      IF( l_trdtrc ) THEN
         DO jn = 1, jp_pisces
            jl = jn + jp_pcs0 - 1
             ztrdpis(:,:,:,jn) = ztrdpis(:,:,:,jn) + tra(:,:,:,jl)
             CALL trd_mod_trc( ztrdpis(:,:,:,jn), jn, jptra_trd_sms, kt )   ! save trends
          END DO
          CALL wrk_dealloc( jpi, jpj, jpk, jp_pisces, ztrdpis ) 
      END IF
      !
      CALL p4z_chk_mass( kt ) ! Mass conservation checking

      IF( nn_timing == 1 )  CALL timing_stop('p4z_sms')
      !
      !
   END SUBROUTINE p4z_sms

   SUBROUTINE p4z_sms_init
      !!----------------------------------------------------------------------
      !!                     ***  p4z_sms_init  ***  
      !!
      !! ** Purpose :   read PISCES namelist
      !!
      !! ** input   :   file 'namelist.trc.s' containing the following
      !!             namelist: natext, natbio, natsms
      !!                       natkriest ("key_kriest")
      !!----------------------------------------------------------------------
      NAMELIST/nampisbio/ nrdttrc, wsbio, xkmort, ferat3, wsbio2, niter1max, niter2max
#if defined key_kriest
      NAMELIST/nampiskrp/ xkr_eta, xkr_zeta, xkr_ncontent, xkr_mass_min, xkr_mass_max
#endif
      NAMELIST/nampisdmp/ ln_pisdmp, nn_pisdmp
      NAMELIST/nampismass/ ln_check_mass
      !!----------------------------------------------------------------------


      REWIND( numnatp )
      READ  ( numnatp, nampisbio )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' Namelist : nampisbio'
         WRITE(numout,*) '    frequence pour la biologie                nrdttrc   =', nrdttrc
         WRITE(numout,*) '    POC sinking speed                         wsbio     =', wsbio
         WRITE(numout,*) '    half saturation constant for mortality    xkmort    =', xkmort
         WRITE(numout,*) '    Fe/C in zooplankton                       ferat3    =', ferat3
         WRITE(numout,*) '    Big particles sinking speed               wsbio2    =', wsbio2
         WRITE(numout,*) '    Maximum number of iterations for POC      niter1max =', niter1max
         WRITE(numout,*) '    Maximum number of iterations for GOC      niter2max =', niter2max
      ENDIF

#if defined key_kriest

      !                               ! nampiskrp : kriest parameters
      !                               ! -----------------------------
      xkr_eta      = 0.62
      xkr_zeta     = 1.62
      xkr_ncontent = 5.7E-6
      xkr_mass_min = 0.0002
      xkr_mass_max = 1.

      REWIND( numnatp )                     ! read natkriest
      READ  ( numnatp, nampiskrp )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : nampiskrp'
         WRITE(numout,*) '    Sinking  exponent                        xkr_eta      = ', xkr_eta
         WRITE(numout,*) '    N content exponent                       xkr_zeta     = ', xkr_zeta
         WRITE(numout,*) '    N content factor                         xkr_ncontent = ', xkr_ncontent
         WRITE(numout,*) '    Minimum mass for Aggregates              xkr_mass_min = ', xkr_mass_min
         WRITE(numout,*) '    Maximum mass for Aggregates              xkr_mass_max = ', xkr_mass_max
         WRITE(numout,*)
     ENDIF


     ! Computation of some variables
     xkr_massp = xkr_ncontent * 7.625 * xkr_mass_min**xkr_zeta

#endif

      ln_pisdmp = .true.
      nn_pisdmp = 1

      REWIND( numnatp )
      READ  ( numnatp, nampisdmp )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : nampisdmp'
         WRITE(numout,*) '    Relaxation of tracer to glodap mean value             ln_pisdmp      =', ln_pisdmp
         WRITE(numout,*) '    Frequency of Relaxation                               nn_pisdmp      =', nn_pisdmp
         WRITE(numout,*) ' '
      ENDIF

      ln_check_mass = .false.
      REWIND( numnatp )       
      READ  ( numnatp, nampismass )
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameter for mass conservation checking'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Flag to check mass conservation of NO3/Si/TALK ln_check_mass = ', ln_check_mass
      ENDIF

   END SUBROUTINE p4z_sms_init

   SUBROUTINE p4z_ph_ini
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE p4z_ini_ph  ***
      !!
      !!  ** Purpose : Initialization of chemical variables of the carbon cycle
      !!---------------------------------------------------------------------
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zcaralk, zbicarb, zco3
      REAL(wp) ::  ztmas, ztmas1
      !!---------------------------------------------------------------------

      ! Set PH from  total alkalinity, borat (???), akb3 (???) and ak23 (???)
      ! --------------------------------------------------------
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztmas   = tmask(ji,jj,jk)
               ztmas1  = 1. - tmask(ji,jj,jk)
               zcaralk = trn(ji,jj,jk,jptal) - borat(ji,jj,jk) / (  1. + 1.E-8 / ( rtrn + akb3(ji,jj,jk) )  )
               zco3    = ( zcaralk - trn(ji,jj,jk,jpdic) ) * ztmas + 0.5e-3 * ztmas1
               zbicarb = ( 2. * trn(ji,jj,jk,jpdic) - zcaralk )
               hi(ji,jj,jk) = ( ak23(ji,jj,jk) * zbicarb / zco3 ) * ztmas + 1.e-9 * ztmas1
            END DO
         END DO
     END DO
     !
   END SUBROUTINE p4z_ph_ini

   SUBROUTINE p4z_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE p4z_rst  ***
      !!
      !!  ** Purpose : Read or write variables in restart file:
      !!
      !!  WRITE(READ) mode:
      !!       kt        : number of time step since the begining of the experiment at the
      !!                   end of the current(previous) run
      !!---------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      !
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zcaralk, zbicarb, zco3
      REAL(wp) ::  ztmas, ztmas1
      !!---------------------------------------------------------------------

      IF( TRIM(cdrw) == 'READ' ) THEN
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' p4z_rst : Read specific variables from pisces model '
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
         ! 
         IF( iom_varid( numrtr, 'PH', ldstop = .FALSE. ) > 0 ) THEN
            CALL iom_get( numrtr, jpdom_autoglo, 'PH' , hi(:,:,:)  )
         ELSE
!            hi(:,:,:) = 1.e-9 
            CALL p4z_ph_ini
         ENDIF
         CALL iom_get( numrtr, jpdom_autoglo, 'Silicalim', xksi(:,:) )
         IF( iom_varid( numrtr, 'Silicamax', ldstop = .FALSE. ) > 0 ) THEN
            CALL iom_get( numrtr, jpdom_autoglo, 'Silicamax' , xksimax(:,:)  )
         ELSE
            xksimax(:,:) = xksi(:,:)
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN
         IF( kt == nitrst ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'p4z_rst : write pisces restart file  kt =', kt
            IF(lwp) WRITE(numout,*) '~~~~~~~'
         ENDIF
         CALL iom_rstput( kt, nitrst, numrtw, 'PH', hi(:,:,:) )
         CALL iom_rstput( kt, nitrst, numrtw, 'Silicalim', xksi(:,:) )
         CALL iom_rstput( kt, nitrst, numrtw, 'Silicamax', xksimax(:,:) )
      ENDIF
      !
   END SUBROUTINE p4z_rst

   SUBROUTINE p4z_dmp( kt )
      !!----------------------------------------------------------------------
      !!                    ***  p4z_dmp  ***
      !!
      !! ** purpose  : Relaxation of some tracers
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT( in )  ::     kt ! time step
      !
      REAL(wp) ::  alkmean = 2426.     ! mean value of alkalinity ( Glodap ; for Goyet 2391. )
      REAL(wp) ::  po4mean = 2.165     ! mean value of phosphates
      REAL(wp) ::  no3mean = 30.90     ! mean value of nitrate
      REAL(wp) ::  silmean = 91.51     ! mean value of silicate
      !
      REAL(wp) :: zarea, zalksum, zpo4sum, zno3sum, zsilsum
      !!---------------------------------------------------------------------


      IF(lwp)  WRITE(numout,*)
      IF(lwp)  WRITE(numout,*) ' p4z_dmp : Restoring of nutrients at time-step kt = ', kt
      IF(lwp)  WRITE(numout,*)

      IF( cp_cfg == "orca" .AND. .NOT. lk_c1d ) THEN      ! ORCA configuration (not 1D) !
         !                                                    ! --------------------------- !
         ! set total alkalinity, phosphate, nitrate & silicate
         zarea          = 1._wp / glob_sum( cvol(:,:,:) ) * 1e6              

         zalksum = glob_sum( trn(:,:,:,jptal) * cvol(:,:,:)  ) * zarea
         zpo4sum = glob_sum( trn(:,:,:,jppo4) * cvol(:,:,:)  ) * zarea * po4r
         zno3sum = glob_sum( trn(:,:,:,jpno3) * cvol(:,:,:)  ) * zarea * rno3
         zsilsum = glob_sum( trn(:,:,:,jpsil) * cvol(:,:,:)  ) * zarea
 
         IF(lwp) WRITE(numout,*) '       TALK mean : ', zalksum
         trn(:,:,:,jptal) = trn(:,:,:,jptal) * alkmean / zalksum

         IF(lwp) WRITE(numout,*) '       PO4  mean : ', zpo4sum
         trn(:,:,:,jppo4) = trn(:,:,:,jppo4) * po4mean / zpo4sum

         IF(lwp) WRITE(numout,*) '       NO3  mean : ', zno3sum
         trn(:,:,:,jpno3) = trn(:,:,:,jpno3) * no3mean / zno3sum

         IF(lwp) WRITE(numout,*) '       SiO3 mean : ', zsilsum
         trn(:,:,:,jpsil) = MIN( 400.e-6,trn(:,:,:,jpsil) * silmean / zsilsum )
         !
      ENDIF

   END SUBROUTINE p4z_dmp


   SUBROUTINE p4z_chk_mass( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_chk_mass  ***
      !!
      !! ** Purpose :  Mass conservation check 
      !!
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      
      !!
      !!---------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN 
         IF( ln_check_mass .AND. lwp) THEN      !   Open budget file of NO3, ALK, Si, Fer
            CALL ctl_opn( numco2, 'carbon.budget'  , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
            CALL ctl_opn( numnut, 'nutrient.budget', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
         ENDIF
      ENDIF

      IF( ln_check_mass .AND. kt == nitend ) THEN      !   Compute the budget of NO3, ALK, Si, Fer
         no3budget = glob_sum( (   trn(:,:,:,jpno3) + trn(:,:,:,jpnh4)  &
            &                    + trn(:,:,:,jpphy) + trn(:,:,:,jpdia)  &
            &                    + trn(:,:,:,jpzoo) + trn(:,:,:,jpmes)  &
            &                    + trn(:,:,:,jppoc)                     &
#if ! defined key_kriest
            &                    + trn(:,:,:,jpgoc)                     &
#endif
            &                    + trn(:,:,:,jpdoc)                     ) * cvol(:,:,:)  ) 
         ! 
         silbudget = glob_sum( (   trn(:,:,:,jpsil) + trn(:,:,:,jpgsi)  &
            &                    + trn(:,:,:,jpdsi)                     ) * cvol(:,:,:)  )
         ! 
         alkbudget = glob_sum( (   trn(:,:,:,jpno3) * rno3              &
            &                    + trn(:,:,:,jptal)                     &
            &                    + trn(:,:,:,jpcal) * 2.                ) * cvol(:,:,:)  )
         ! 
         ferbudget = glob_sum( (   trn(:,:,:,jpfer) + trn(:,:,:,jpnfe)  &
            &                    + trn(:,:,:,jpdfe)                     &
#if ! defined key_kriest
            &                    + trn(:,:,:,jpbfe)                     &
#endif
            &                    + trn(:,:,:,jpsfe)                     &
            &                    + trn(:,:,:,jpzoo) * ferat3            &
            &                    + trn(:,:,:,jpmes) * ferat3            ) * cvol(:,:,:)  )

         !
         t_atm_co2_flx  = t_atm_co2_flx / glob_sum( e1e2t(:,:) )
         t_oce_co2_flx  = t_oce_co2_flx         * 12. / 1.e15 * (-1 )
         tpp            = tpp           * 1000. * 12. / 1.E15
         t_oce_co2_exp  = t_oce_co2_exp * 1000. * 12. / 1.E15
         !
         no3budget = no3budget / areatot
         silbudget = silbudget / areatot
         alkbudget = alkbudget / areatot
         ferbudget = ferbudget / areatot
         !
         IF(lwp) THEN
            WRITE(numco2,9000) ndastp, t_atm_co2_flx, t_oce_co2_flx, tpp, t_oce_co2_exp
            WRITE(numnut,9500) ndastp, alkbudget, no3budget, silbudget, ferbudget
         ENDIF
         !
      ENDIF
       !
 9000  FORMAT(i8,f10.5,e18.10,f10.5,f10.5)
 9500  FORMAT(i8,4e18.10)     
       !
   END SUBROUTINE p4z_chk_mass

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sms( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p4z_sms: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_sms
#endif 

   !!======================================================================
END MODULE p4zsms 
