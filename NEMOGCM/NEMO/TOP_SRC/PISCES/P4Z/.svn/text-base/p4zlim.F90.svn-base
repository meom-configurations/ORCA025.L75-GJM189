MODULE p4zlim
   !!======================================================================
   !!                         ***  MODULE p4zlim  ***
   !! TOP :   PISCES 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-04  (O. Aumont, C. Ethe) Limitation for iron modelled in quota 
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_lim        :   Compute the nutrients limitation terms 
   !!   p4z_lim_init   :   Read the namelist 
   !!----------------------------------------------------------------------
   USE oce_trc         ! Shared ocean-passive tracers variables
   USE trc             ! Tracers defined
   USE sms_pisces      ! PISCES variables
   USE p4zopt          ! Optical
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC p4z_lim    
   PUBLIC p4z_lim_init    

   !! * Shared module variables
   REAL(wp), PUBLIC ::  concnno3  = 1.e-6_wp      !:  NO3, PO4 half saturation   
   REAL(wp), PUBLIC ::  concdno3  = 2.e-6_wp      !:  Phosphate half saturation for diatoms  
   REAL(wp), PUBLIC ::  concnnh4  = 1.e-7_wp      !:  NH4 half saturation for phyto  
   REAL(wp), PUBLIC ::  concdnh4  = 2.e-7_wp      !:  NH4 half saturation for diatoms
   REAL(wp), PUBLIC ::  concnfer  = 1.e-9_wp      !:  Iron half saturation for nanophyto 
   REAL(wp), PUBLIC ::  concdfer  = 2.e-9_wp      !:  Iron half saturation for diatoms  
   REAL(wp), PUBLIC ::  concbno3  = 2.5e-7_wp      !:  NO3 half saturation  for bacteria 
   REAL(wp), PUBLIC ::  concbnh4  = 2.5e-8_wp      !:  NH4 half saturation for bacteria
   REAL(wp), PUBLIC ::  xsizedia  = 5.e-7_wp      !:  Minimum size criteria for diatoms
   REAL(wp), PUBLIC ::  xsizephy  = 1.e-6_wp      !:  Minimum size criteria for nanophyto
   REAL(wp), PUBLIC ::  xsizern   = 3.0_wp        !:  Size ratio for nanophytoplankton
   REAL(wp), PUBLIC ::  xsizerd   = 4.0_wp        !:  Size ratio for diatoms
   REAL(wp), PUBLIC ::  xksi1     = 2.E-6_wp      !:  half saturation constant for Si uptake 
   REAL(wp), PUBLIC ::  xksi2     = 3.33e-6_wp    !:  half saturation constant for Si/C 
   REAL(wp), PUBLIC ::  xkdoc     = 417.e-6_wp    !:  2nd half-sat. of DOC remineralization  
   REAL(wp), PUBLIC ::  concbfe   = 1.E-11_wp     !:  Fe half saturation for bacteria 
   REAL(wp), PUBLIC ::  qnfelim   = 7.E-6_wp      !:  optimal Fe quota for nanophyto
   REAL(wp), PUBLIC ::  qdfelim   = 7.E-6_wp      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  caco3r    = 0.16_wp       !:  mean rainratio 

   ! Coefficient for iron limitation
   REAL(wp) ::  xcoef1   = 0.0016  / 55.85  
   REAL(wp) ::  xcoef2   = 1.21E-5 * 14. / 55.85 / 7.625 * 0.5 * 1.5
   REAL(wp) ::  xcoef3   = 1.15E-4 * 14. / 55.85 / 7.625 * 0.5 
   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zlim.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_lim( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_lim  ***
      !!
      !! ** Purpose :   Compute the co-limitations by the various nutrients
      !!              for the various phytoplankton species
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in)  :: kt, jnt
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zlim1, zlim2, zlim3, zlim4, zno3, zferlim
      REAL(wp) ::   zconcd, zconcd2, zconcn, zconcn2
      REAL(wp) ::   z1_trndia, z1_trnphy, ztem1, ztem2, zetot1, zetot2
      REAL(wp) ::   zdenom, zratio, zironmin
      REAL(wp) ::   zconc1d, zconc1dnh4, zconc0n, zconc0nnh4   
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_lim')
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               
               ! Tuning of the iron concentration to a minimum level that is set to the detection limit
               !-------------------------------------
               zno3    = trn(ji,jj,jk,jpno3) / 40.e-6
               zferlim = MAX( 3e-11 * zno3 * zno3, 5e-12 )
               zferlim = MIN( zferlim, 5e-11 )
               trn(ji,jj,jk,jpfer) = MAX( trn(ji,jj,jk,jpfer), zferlim )

               ! Computation of a variable Ks for iron on diatoms taking into account
               ! that increasing biomass is made of generally bigger cells
               !------------------------------------------------
               zconcd   = MAX( 0.e0 , trn(ji,jj,jk,jpdia) - xsizedia )
               zconcd2  = trn(ji,jj,jk,jpdia) - zconcd
               zconcn   = MAX( 0.e0 , trn(ji,jj,jk,jpphy) - xsizephy )
               zconcn2  = trn(ji,jj,jk,jpphy) - zconcn
               z1_trnphy   = 1. / ( trn(ji,jj,jk,jpphy) + rtrn )
               z1_trndia   = 1. / ( trn(ji,jj,jk,jpdia) + rtrn )

               concdfe(ji,jj,jk) = MAX( concdfer, ( zconcd2 * concdfer + concdfer * xsizerd * zconcd ) * z1_trndia )
               zconc1d           = MAX( concdno3, ( zconcd2 * concdno3 + concdno3 * xsizerd * zconcd ) * z1_trndia )
               zconc1dnh4        = MAX( concdnh4, ( zconcd2 * concdnh4 + concdnh4 * xsizerd * zconcd ) * z1_trndia )

               concnfe(ji,jj,jk) = MAX( concnfer, ( zconcn2 * concnfer + concnfer * xsizern * zconcn ) * z1_trnphy )
               zconc0n           = MAX( concnno3, ( zconcn2 * concnno3 + concnno3 * xsizern * zconcn ) * z1_trnphy )
               zconc0nnh4        = MAX( concnnh4, ( zconcn2 * concnnh4 + concnnh4 * xsizern * zconcn ) * z1_trnphy )

               ! Michaelis-Menten Limitation term for nutrients Small bacteria
               ! -------------------------------------------------------------
               zdenom = 1. /  ( concbno3 * concbnh4 + concbnh4 * trn(ji,jj,jk,jpno3) + concbno3 * trn(ji,jj,jk,jpnh4) )
               xnanono3(ji,jj,jk) = trn(ji,jj,jk,jpno3) * concbnh4 * zdenom
               xnanonh4(ji,jj,jk) = trn(ji,jj,jk,jpnh4) * concbno3 * zdenom
               !
               zlim1    = xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk)
               zlim2    = trn(ji,jj,jk,jppo4) / ( trn(ji,jj,jk,jppo4) + concbnh4 )
               zlim3    = trn(ji,jj,jk,jpfer) / ( concbfe + trn(ji,jj,jk,jpfer) )
               zlim4    = trn(ji,jj,jk,jpdoc) / ( xkdoc   + trn(ji,jj,jk,jpdoc) )
               xlimbacl(ji,jj,jk) = MIN( zlim1, zlim2, zlim3 )
               xlimbac (ji,jj,jk) = MIN( zlim1, zlim2, zlim3 ) * zlim4

               ! Michaelis-Menten Limitation term for nutrients Small flagellates
               ! -----------------------------------------------
               zdenom = 1. /  ( zconc0n * zconc0nnh4 + zconc0nnh4 * trn(ji,jj,jk,jpno3) + zconc0n * trn(ji,jj,jk,jpnh4) )
               xnanono3(ji,jj,jk) = trn(ji,jj,jk,jpno3) * zconc0nnh4 * zdenom
               xnanonh4(ji,jj,jk) = trn(ji,jj,jk,jpnh4) * zconc0n    * zdenom
               !
               zlim1    = xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk)
               zlim2    = trn(ji,jj,jk,jppo4) / ( trn(ji,jj,jk,jppo4) + zconc0nnh4 )
               zratio   = trn(ji,jj,jk,jpnfe) * z1_trnphy 
               zironmin = xcoef1 * trn(ji,jj,jk,jpnch) * z1_trnphy + xcoef2 * zlim1 + xcoef3 * xnanono3(ji,jj,jk)
               zlim3    = MAX( 0.,( zratio - zironmin ) / qnfelim )
               xlimnfe(ji,jj,jk) = MIN( 1., zlim3 )
               xlimphy(ji,jj,jk) = MIN( zlim1, zlim2, zlim3 )
               !
               !   Michaelis-Menten Limitation term for nutrients Diatoms
               !   ----------------------------------------------
               zdenom   = 1. / ( zconc1d * zconc1dnh4 + zconc1dnh4 * trn(ji,jj,jk,jpno3) + zconc1d * trn(ji,jj,jk,jpnh4) )
               xdiatno3(ji,jj,jk) = trn(ji,jj,jk,jpno3) * zconc1dnh4 * zdenom
               xdiatnh4(ji,jj,jk) = trn(ji,jj,jk,jpnh4) * zconc1d    * zdenom
               !
               zlim1    = xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk)
               zlim2    = trn(ji,jj,jk,jppo4) / ( trn(ji,jj,jk,jppo4) + zconc1dnh4  )
               zlim3    = trn(ji,jj,jk,jpsil) / ( trn(ji,jj,jk,jpsil) + xksi(ji,jj) )
               zratio   = trn(ji,jj,jk,jpdfe) * z1_trndia
               zironmin = xcoef1 * trn(ji,jj,jk,jpdch) * z1_trndia + xcoef2 * zlim1 + xcoef3 * xdiatno3(ji,jj,jk)
               zlim4    = MAX( 0., ( zratio - zironmin ) / qdfelim )
               xlimdfe(ji,jj,jk) = MIN( 1., zlim4 )
               xlimdia(ji,jj,jk) = MIN( zlim1, zlim2, zlim3, zlim4 )
               xlimsi(ji,jj,jk)  = MIN( zlim1, zlim2, zlim4 )
           END DO
         END DO
      END DO

      ! Compute the fraction of nanophytoplankton that is made of calcifiers
      ! --------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zlim1 =  ( trn(ji,jj,jk,jpno3) * concnnh4 + trn(ji,jj,jk,jpnh4) * concnno3 )    &
                  &   / ( concnno3 * concnnh4 + concnnh4 * trn(ji,jj,jk,jpno3) + concnno3 * trn(ji,jj,jk,jpnh4) ) 
               zlim2  = trn(ji,jj,jk,jppo4) / ( trn(ji,jj,jk,jppo4) + concnnh4 )
               zlim3  = trn(ji,jj,jk,jpfer) / ( trn(ji,jj,jk,jpfer) +  5.E-11   )
               ztem1  = MAX( 0., tsn(ji,jj,jk,jp_tem) )
               ztem2  = tsn(ji,jj,jk,jp_tem) - 10.
               zetot1 = MAX( 0., etot(ji,jj,jk) - 1.) / ( 4. + etot(ji,jj,jk) ) 
               zetot2 = 30. / ( 30. + etot(ji,jj,jk) ) 

               xfracal(ji,jj,jk) = caco3r * MIN( zlim1, zlim2, zlim3 )                  &
                  &                       * ztem1 / ( 0.1 + ztem1 )                     &
                  &                       * MAX( 1., trn(ji,jj,jk,jpphy) * 1.e6 / 2. )  &
                  &                       * zetot1 * zetot2               &
                  &                       * ( 1. + EXP(-ztem2 * ztem2 / 25. ) )         &
                  &                       * MIN( 1., 50. / ( hmld(ji,jj) + rtrn ) )
               xfracal(ji,jj,jk) = MIN( 0.8 , xfracal(ji,jj,jk) )
               xfracal(ji,jj,jk) = MAX( 0.02, xfracal(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
      IF( ln_diatrc .AND. lk_iomput .AND. jnt == nrdttrc ) THEN        ! save output diagnostics
        !
        CALL iom_put( "xfracal", xfracal(:,:,:) * tmask(:,:,:) )  ! euphotic layer deptht
        CALL iom_put( "LNnut"  , xlimphy(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        CALL iom_put( "LDnut"  , xlimdia(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        CALL iom_put( "LNFe"   , xlimnfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        CALL iom_put( "LDFe"   , xlimdfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        !
      ENDIF

      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_lim')
      !
   END SUBROUTINE p4z_lim

   SUBROUTINE p4z_lim_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_lim_init  ***
      !!
      !! ** Purpose :   Initialization of nutrient limitation parameters
      !!
      !! ** Method  :   Read the nampislim namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampislim
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampislim/ concnno3, concdno3, concnnh4, concdnh4, concnfer, concdfer, concbfe,   &
         &                concbno3, concbnh4, xsizedia, xsizephy, xsizern, xsizerd,          & 
         &                xksi1, xksi2, xkdoc, qnfelim, qdfelim, caco3r

      REWIND( numnatp )                     ! read numnat
      READ  ( numnatp, nampislim )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for nutrient limitations, nampislim'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    mean rainratio                           caco3r    = ', caco3r
         WRITE(numout,*) '    NO3 half saturation of nanophyto         concnno3  = ', concnno3
         WRITE(numout,*) '    NO3 half saturation of diatoms           concdno3  = ', concdno3
         WRITE(numout,*) '    NH4 half saturation for phyto            concnnh4  = ', concnnh4
         WRITE(numout,*) '    NH4 half saturation for diatoms          concdnh4  = ', concdnh4
         WRITE(numout,*) '    half saturation constant for Si uptake   xksi1     = ', xksi1
         WRITE(numout,*) '    half saturation constant for Si/C        xksi2     = ', xksi2
         WRITE(numout,*) '    half-sat. of DOC remineralization        xkdoc     = ', xkdoc
         WRITE(numout,*) '    Iron half saturation for nanophyto       concnfer  = ', concnfer
         WRITE(numout,*) '    Iron half saturation for diatoms         concdfer  = ', concdfer
         WRITE(numout,*) '    size ratio for nanophytoplankton         xsizern   = ', xsizern
         WRITE(numout,*) '    size ratio for diatoms                   xsizerd   = ', xsizerd
         WRITE(numout,*) '    NO3 half saturation of bacteria          concbno3  = ', concbno3
         WRITE(numout,*) '    NH4 half saturation for bacteria         concbnh4  = ', concbnh4
         WRITE(numout,*) '    Minimum size criteria for diatoms        xsizedia  = ', xsizedia
         WRITE(numout,*) '    Minimum size criteria for nanophyto      xsizephy  = ', xsizephy
         WRITE(numout,*) '    Fe half saturation for bacteria          concbfe   = ', concbfe
         WRITE(numout,*) '    optimal Fe quota for nano.               qnfelim   = ', qnfelim
         WRITE(numout,*) '    Optimal Fe quota for diatoms             qdfelim   = ', qdfelim
      ENDIF

   END SUBROUTINE p4z_lim_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_lim                   ! Empty routine
   END SUBROUTINE p4z_lim
#endif 

   !!======================================================================
END MODULE  p4zlim
