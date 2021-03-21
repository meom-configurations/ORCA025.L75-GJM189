MODULE p4zsed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :   PISCES Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06 (C. Ethe) USE of fldread
   !!             3.5  !  2012-07 (O. Aumont) improvment of river input of nutrients 
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_sed        :  Compute loss of organic matter in the sediments
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zopt          !  optical model
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zrem          !  Remineralisation of organic matter
   USE p4zsbc          !  External source of nutrients 
   USE p4zint          !  interpolation and computation of various fields
   USE iom             !  I/O manager
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sed   

   !! * Module variables
   REAL(wp) :: ryyss                    !: number of seconds per year 
   REAL(wp) :: r1_ryyss                 !: inverse of ryyss
   REAL(wp) :: rmtss                    !: number of seconds per month
   REAL(wp) :: r1_rday                  !: inverse of rday

   INTEGER ::  numnit  


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Header:$ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sed( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed  ***
      !!
      !! ** Purpose :   Compute loss of organic matter in the sediments. This
      !!              is by no way a sediment model. The loss is simply 
      !!              computed to balance the inout from rivers and dust
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER  ::   ji, jj, jk, ikt
#if ! defined key_sed
      REAL(wp) ::   zsumsedsi, zsumsedpo4, zsumsedcal
      REAL(wp) ::   zrivalk, zrivsil, zrivno3
#endif
      REAL(wp) ::  zwflux, zfminus, zfplus
      REAL(wp) ::  zlim, zfact, zfactcal
      REAL(wp) ::  zo2, zno3, zflx, zpdenit, z1pdenit, zdenitt, zolimit
      REAL(wp) ::  zsiloss, zcaloss, zwsbio3, zwsbio4, zwscal, zdep, zwstpoc
      REAL(wp) ::  ztrfer, ztrpo4, zwdust
      REAL(wp) ::  zrdenittot, zsdenittot, znitrpottot
      !
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zpdep, zsidep, zwork1, zwork2, zwork3, zwork4
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zdenit2d, zironice
      REAL(wp), POINTER, DIMENSION(:,:,:) :: znitrpot, zirondep
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sed')
      !
      IF( kt == nittrc000 .AND. jnt == 1 )  THEN
         ryyss    = nyear_len(1) * rday    ! number of seconds per year and per month
         rmtss    = ryyss / raamo
         r1_rday  = 1. / rday
         r1_ryyss = 1. / ryyss
         IF( ln_check_mass .AND. lwp)  &
           &  CALL ctl_opn( numnit, 'nitrogen.budget', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
      ENDIF
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, zdenit2d, zwork1, zwork2, zwork3, zwork4 )
      CALL wrk_alloc( jpi, jpj, jpk, znitrpot )
      zdenit2d(:,:) = 0.e0

      ! Iron input/uptake due to sea ice : Crude parameterization based on Lancelot et al.
      ! ----------------------------------------------------
      IF( ln_ironice ) THEN  
         !                                              
         CALL wrk_alloc( jpi, jpj, zironice )
         !                                              
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep    = rfact2 / fse3t(ji,jj,1)
               zwflux  = fmmflx(ji,jj) / 1000._wp
               zfminus = MIN( 0._wp, -zwflux ) * trn(ji,jj,1,jpfer) * zdep
               zfplus  = MAX( 0._wp, -zwflux ) * icefeinput * zdep
               zironice(ji,jj) =  zfplus + zfminus
            END DO
         END DO
         !
         trn(:,:,1,jpfer) = trn(:,:,1,jpfer) + zironice(:,:) 
         !                                              
         IF( ln_diatrc .AND. lk_iomput .AND. jnt == nrdttrc )   &
            &   CALL iom_put( "Ironice", zironice(:,:) * 1.e+3 * rfact2r * fse3t(:,:,1) * tmask(:,:,1) ) ! iron flux from ice
         CALL wrk_dealloc( jpi, jpj, zironice )
         !                                              
      ENDIF

      ! Add the external input of nutrients from dust deposition
      ! ----------------------------------------------------------
      IF( ln_dust ) THEN
         !                                              
         CALL wrk_alloc( jpi, jpj,      zpdep, zsidep )
         CALL wrk_alloc( jpi, jpj, jpk, zirondep      )
         !                                              ! Iron and Si deposition at the surface
         IF( ln_solub ) THEN
            zirondep(:,:,1) = solub(:,:) * dust(:,:) * rfact2 / fse3t(:,:,1) / ( 55.85 * rmtss ) + 3.e-10 * r1_ryyss 
         ELSE
            zirondep(:,:,1) = dustsolub  * dust(:,:) * rfact2 / fse3t(:,:,1) / ( 55.85 * rmtss ) + 3.e-10 * r1_ryyss 
         ENDIF
         zsidep(:,:) = 8.8 * 0.075 * dust(:,:) * rfact2 / fse3t(:,:,1) / ( 28.1  * rmtss )
         zpdep (:,:) = 0.1 * 0.021 * dust(:,:) * rfact2 / fse3t(:,:,1) / ( 31.   * rmtss ) / po4r 
         !                                              ! Iron solubilization of particles in the water column
         zwdust = 0.005 / ( wdust * 55.85 * 30.42 ) / ( 45. * rday ) 
         DO jk = 2, jpkm1
            zirondep(:,:,jk) = dust(:,:) * zwdust * rfact2 * EXP( -fsdept(:,:,jk) / 1000. )
         END DO
         !                                              ! Iron solubilization of particles in the water column
         trn(:,:,1,jppo4) = trn(:,:,1,jppo4) + zpdep   (:,:)
         trn(:,:,1,jpsil) = trn(:,:,1,jpsil) + zsidep  (:,:)
         trn(:,:,:,jpfer) = trn(:,:,:,jpfer) + zirondep(:,:,:) 
         !                                              
         IF( ln_diatrc ) THEN
            zfact = 1.e+3 * rfact2r
            IF( lk_iomput ) THEN
               IF( jnt == nrdttrc ) THEN
                  CALL iom_put( "Irondep", zirondep(:,:,1) * zfact * fse3t(:,:,1) * tmask(:,:,1) ) ! surface downward dust depo of iron
                  CALL iom_put( "pdust"  , dust(:,:) / ( wdust * 30.42 * 0.035 )  * tmask(:,:,1) ) ! dust concentration at surface
               ENDIF
            ELSE
               trc2d(:,:,jp_pcs0_2d + 11) = zirondep(:,:,1) * zfact * fse3t(:,:,1) * tmask(:,:,1)
            ENDIF
         ENDIF
         CALL wrk_dealloc( jpi, jpj,      zpdep, zsidep )
         CALL wrk_dealloc( jpi, jpj, jpk, zirondep      )
         !                                              
      ENDIF
     
      ! Add the external input of nutrients from river
      ! ----------------------------------------------------------
      IF( ln_river ) THEN
         trn(:,:,1,jppo4) = trn(:,:,1,jppo4) + rivdip(:,:) * rfact2
         trn(:,:,1,jpno3) = trn(:,:,1,jpno3) + rivdin(:,:) * rfact2
         trn(:,:,1,jpfer) = trn(:,:,1,jpfer) + rivdic(:,:) * 5.e-5 * rfact2
         trn(:,:,1,jpsil) = trn(:,:,1,jpsil) + rivdsi(:,:) * rfact2
         trn(:,:,1,jpdic) = trn(:,:,1,jpdic) + rivdic(:,:) * rfact2
         trn(:,:,1,jptal) = trn(:,:,1,jptal) + ( rivalk(:,:) - rno3 * rivdin(:,:) ) * rfact2
      ENDIF
      
      ! Add the external input of nutrients from nitrogen deposition
      ! ----------------------------------------------------------
      IF( ln_ndepo ) THEN
         trn(:,:,1,jpno3) = trn(:,:,1,jpno3) + nitdep(:,:) * rfact2
         trn(:,:,1,jptal) = trn(:,:,1,jptal) - rno3 * nitdep(:,:) * rfact2
      ENDIF

      ! Add the external input of iron from sediment mobilization
      ! ------------------------------------------------------
      IF( ln_ironsed ) THEN
         trn(:,:,:,jpfer) = trn(:,:,:,jpfer) + ironsed(:,:,:) * rfact2
         !
         IF( ln_diatrc .AND. lk_iomput .AND. jnt == nrdttrc )   &
            &   CALL iom_put( "Ironsed", ironsed(:,:,:) * 1.e+3 * tmask(:,:,:) ) ! iron inputs from sediments
      ENDIF

      ! Add the external input of iron from hydrothermal vents
      ! ------------------------------------------------------
      IF( ln_hydrofe ) THEN
         trn(:,:,:,jpfer) = trn(:,:,:,jpfer) + hydrofe(:,:,:) * rfact2
         !
         IF( ln_diatrc .AND. lk_iomput .AND. jnt == nrdttrc )   &
            &   CALL iom_put( "HYDR", hydrofe(:,:,:) * 1.e+3 * tmask(:,:,:) ) ! hydrothermal iron input
      ENDIF

#if ! defined key_sed
      ! Computation of the sediment denitrification proportion: The metamodel from midlleburg (2006) is being used
      ! -------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
           IF( tmask(ji,jj,1) == 1 ) THEN
              ikt = mbkt(ji,jj)
# if defined key_kriest
              zflx =    trn(ji,jj,ikt,jppoc) * wsbio3(ji,jj,ikt)    * 1E3 * 1E6 / 1E4
# else
              zflx = (  trn(ji,jj,ikt,jpgoc) * wsbio4(ji,jj,ikt)   &
                &     + trn(ji,jj,ikt,jppoc) * wsbio3(ji,jj,ikt) )  * 1E3 * 1E6 / 1E4
#endif
              zflx  = LOG10( MAX( 1E-3, zflx ) )
              zo2   = LOG10( MAX( 10. , trn(ji,jj,ikt,jpoxy) * 1E6 ) )
              zno3  = LOG10( MAX( 1.  , trn(ji,jj,ikt,jpno3) * 1E6 * rno3 ) )
              zdep  = LOG10( fsdepw(ji,jj,ikt+1) )
              zdenit2d(ji,jj) = -2.2567 - 1.185 * zflx - 0.221 * zflx**2 - 0.3995 * zno3 * zo2 + 1.25 * zno3    &
              &                + 0.4721 * zo2 - 0.0996 * zdep + 0.4256 * zflx * zo2
              zdenit2d(ji,jj) = 10.0**( zdenit2d(ji,jj) )
           ENDIF
         END DO
      END DO 
      ! Loss of biogenic silicon, Caco3 organic carbon in the sediments. 
      ! First, the total loss is computed.
      ! The factor for calcite comes from the alkalinity effect
      ! -------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt = mbkt(ji,jj) 
# if defined key_kriest
            zwork1(ji,jj) = trn(ji,jj,ikt,jpgsi) * wscal (ji,jj,ikt)
            zwork2(ji,jj) = trn(ji,jj,ikt,jppoc) * wsbio3(ji,jj,ikt)
# else
            zwork1(ji,jj) = trn(ji,jj,ikt,jpgsi) * wsbio4(ji,jj,ikt)
            zwork2(ji,jj) = trn(ji,jj,ikt,jpgoc) * wsbio4(ji,jj,ikt) + trn(ji,jj,ikt,jppoc) * wsbio3(ji,jj,ikt) 
# endif
            ! For calcite, burial efficiency is made a function of saturation
            zfactcal      = MIN( excess(ji,jj,ikt), 0.2 )
            zfactcal      = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )
            zwork3(ji,jj) = trn(ji,jj,ikt,jpcal) * wscal (ji,jj,ikt) * 2.e0 * zfactcal
         END DO
      END DO
      zsumsedsi  = glob_sum( zwork1(:,:) * e1e2t(:,:) ) * r1_rday
      zsumsedpo4 = glob_sum( zwork2(:,:) * e1e2t(:,:) ) * r1_rday
      zsumsedcal = glob_sum( zwork3(:,:) * e1e2t(:,:) ) * r1_rday
#endif

      ! This loss is scaled at each bottom grid cell for equilibrating the total budget of silica in the ocean.
      ! Thus, the amount of silica lost in the sediments equal the supply at the surface (dust+rivers)
      ! ------------------------------------------------------
#if ! defined key_sed
      zrivsil =  1._wp - ( sumdepsi + rivdsiinput * r1_ryyss ) / zsumsedsi
      zrivno3 =  1._wp - ( rivdininput * r1_ryyss ) / zsumsedpo4
#endif

      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstep / fse3t(ji,jj,ikt)
            zwsbio4 = wsbio4(ji,jj,ikt) * zdep
            zwscal  = wscal (ji,jj,ikt) * zdep
# if defined key_kriest
            zsiloss = trn(ji,jj,ikt,jpgsi) * zwsbio4
# else
            zsiloss = trn(ji,jj,ikt,jpgsi) * zwscal
# endif
            zcaloss = trn(ji,jj,ikt,jpcal) * zwscal
            !
            trn(ji,jj,ikt,jpgsi) = trn(ji,jj,ikt,jpgsi) - zsiloss
            trn(ji,jj,ikt,jpcal) = trn(ji,jj,ikt,jpcal) - zcaloss
#if ! defined key_sed
            trn(ji,jj,ikt,jpsil) = trn(ji,jj,ikt,jpsil) + zsiloss * zrivsil 
            zfactcal = MIN( excess(ji,jj,ikt), 0.2 )
            zfactcal = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )
            zrivalk  =  1._wp - ( rivalkinput * r1_ryyss ) * zfactcal / zsumsedcal 
            trn(ji,jj,ikt,jptal) =  trn(ji,jj,ikt,jptal) + zcaloss * zrivalk * 2.0
            trn(ji,jj,ikt,jpdic) =  trn(ji,jj,ikt,jpdic) + zcaloss * zrivalk
#endif
         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt     = mbkt(ji,jj)
            zdep    = xstep / fse3t(ji,jj,ikt)
            zwsbio4 = wsbio4(ji,jj,ikt) * zdep
            zwsbio3 = wsbio3(ji,jj,ikt) * zdep
# if ! defined key_kriest
            trn(ji,jj,ikt,jpgoc) = trn(ji,jj,ikt,jpgoc) - trn(ji,jj,ikt,jpgoc) * zwsbio4
            trn(ji,jj,ikt,jppoc) = trn(ji,jj,ikt,jppoc) - trn(ji,jj,ikt,jppoc) * zwsbio3
            trn(ji,jj,ikt,jpbfe) = trn(ji,jj,ikt,jpbfe) - trn(ji,jj,ikt,jpbfe) * zwsbio4
            trn(ji,jj,ikt,jpsfe) = trn(ji,jj,ikt,jpsfe) - trn(ji,jj,ikt,jpsfe) * zwsbio3
            zwstpoc =  trn(ji,jj,ikt,jpgoc) * zwsbio4 + trn(ji,jj,ikt,jppoc) * zwsbio3 
# else
            trn(ji,jj,ikt,jpnum) = trn(ji,jj,ikt,jpnum) - trn(ji,jj,ikt,jpnum) * zwsbio4
            trn(ji,jj,ikt,jppoc) = trn(ji,jj,ikt,jppoc) - trn(ji,jj,ikt,jppoc) * zwsbio3
            trn(ji,jj,ikt,jpsfe) = trn(ji,jj,ikt,jpsfe) - trn(ji,jj,ikt,jpsfe) * zwsbio3
            zwstpoc = trn(ji,jj,ikt,jppoc) * zwsbio3 
# endif

#if ! defined key_sed
            ! The 0.5 factor in zpdenit and zdenitt is to avoid negative NO3 concentration after both denitrification
            ! in the sediments and just above the sediments. Not very clever, but simpliest option.
            zpdenit  = MIN( 0.5 * ( trn(ji,jj,ikt,jpno3) - rtrn ) / rdenit, zdenit2d(ji,jj) * zwstpoc * zrivno3 )
            z1pdenit = zwstpoc * zrivno3 - zpdenit
            zolimit = MIN( ( trn(ji,jj,ikt,jpoxy) - rtrn ) / o2ut, z1pdenit * ( 1.- nitrfac(ji,jj,ikt) ) )
            zdenitt = MIN(  0.5 * ( trn(ji,jj,ikt,jpno3) - rtrn ) / rdenit, z1pdenit * nitrfac(ji,jj,ikt) )
            trn(ji,jj,ikt,jpdoc) = trn(ji,jj,ikt,jpdoc) + z1pdenit - zolimit - zdenitt
            trn(ji,jj,ikt,jppo4) = trn(ji,jj,ikt,jppo4) + zpdenit + zolimit + zdenitt
            trn(ji,jj,ikt,jpnh4) = trn(ji,jj,ikt,jpnh4) + zpdenit + zolimit + zdenitt
            trn(ji,jj,ikt,jpno3) = trn(ji,jj,ikt,jpno3) - rdenit * (zpdenit + zdenitt)
            trn(ji,jj,ikt,jpoxy) = trn(ji,jj,ikt,jpoxy) - zolimit * o2ut
            trn(ji,jj,ikt,jptal) = trn(ji,jj,ikt,jptal) + rno3 * (zolimit + (1.+rdenit) * (zpdenit + zdenitt) )
            trn(ji,jj,ikt,jpdic) = trn(ji,jj,ikt,jpdic) + zpdenit + zolimit + zdenitt
            zwork4(ji,jj) = rdenit * zpdenit * fse3t(ji,jj,ikt)
#endif
         END DO
      END DO

      ! Nitrogen fixation process
      !-----------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               !                                       ! Potential nitrogen fixation dependant on temperature and iron
               zlim = ( 1.- xnanono3(ji,jj,jk) - xnanonh4(ji,jj,jk) )
               IF( zlim <= 0.2 )   zlim = 0.01
#if defined key_degrad
               zfact = zlim * rfact2 * facvol(ji,jj,jk)
#else
               zfact = zlim * rfact2
#endif
               ztrfer = biron(ji,jj,jk)       / ( concfediaz + biron(ji,jj,jk)       )
               ztrpo4 = trn  (ji,jj,jk,jppo4) / ( concnnh4   + trn  (ji,jj,jk,jppo4) ) 
               znitrpot(ji,jj,jk) =  MAX( 0.e0, ( 0.6 * tgfunc(ji,jj,jk) - 2.15 ) * r1_rday )   &
                 &         *  zfact * MIN( ztrfer, ztrpo4 ) * ( 1.- EXP( -etot(ji,jj,jk) / diazolight ) )
            END DO
         END DO
      END DO
 
      IF( ln_check_mass ) THEN
         ! The total gain from nitrogen fixation is scaled to balance the loss by denitrification
         ! -------------------------------------------------------------
         zrdenittot   = glob_sum ( denitr(:,:,:) * rdenit * xnegtr(:,:,:) * cvol(:,:,:) )
         zsdenittot   = glob_sum ( zwork4(:,:)   * e1e2t(:,:) )
         znitrpottot  = glob_sum ( znitrpot(:,:,:)                        * cvol(:,:,:) )
         IF( kt == nitend .AND. jnt == nrdttrc ) THEN
            zfact = 1.e+3 * rfact2r * rno3 * 365. * 86400. * 14. / 1e12
            IF(lwp) WRITE(numnit,9100) ndastp, znitrpottot * nitrfix * zfact, zrdenittot * zfact , zsdenittot * zfact
         ENDIF
       ENDIF

      ! Nitrogen change due to nitrogen fixation
      ! ----------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! zfact = znitrpot(ji,jj,jk) * ( zrdenittot + zsdenittot ) / znitrpottot
               zfact = znitrpot(ji,jj,jk) * nitrfix
               trn(ji,jj,jk,jpnh4) = trn(ji,jj,jk,jpnh4) +             zfact
               trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal) + rno3      * zfact
               trn(ji,jj,jk,jpoxy) = trn(ji,jj,jk,jpoxy) + o2nit     * zfact 
               trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4) + 30. / 46. * zfact
           END DO
         END DO 
      END DO
      !
      IF( ln_diatrc ) THEN
         zfact = 1.e+3 * rfact2r
         IF( lk_iomput ) THEN
            IF( jnt == nrdttrc ) THEN
               CALL iom_put( "Nfix"  , znitrpot(:,:,:) * nitrfix * rno3 * zfact * tmask(:,:,:) )  ! nitrogen fixation 
               CALL iom_put( "Sdenit", zwork4(:,:)               * rno3 * zfact * tmask(:,:,1) )  ! Nitrate reduction in the sediments
            ENDIF
         ELSE
            trc2d(:,:,jp_pcs0_2d + 12) = znitrpot(:,:,1) * nitrfix * zfact * fse3t(:,:,1) * tmask(:,:,1)
         ENDIF
      ENDIF
      !
      IF(ln_ctl) THEN  ! print mean trends (USEd for debugging)
         WRITE(charout, fmt="('sed ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=trn, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zdenit2d, zwork1, zwork2, zwork3, zwork4 )
      CALL wrk_dealloc( jpi, jpj, jpk, znitrpot )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sed')
      !
 9100  FORMAT(i8,3f10.5)
      !
   END SUBROUTINE p4z_sed

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sed                         ! Empty routine
   END SUBROUTINE p4z_sed
#endif 

   !!======================================================================
END MODULE  p4zsed
