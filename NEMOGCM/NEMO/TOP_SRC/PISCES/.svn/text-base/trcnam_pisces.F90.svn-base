MODULE trcnam_pisces
   !!======================================================================
   !!                      ***  MODULE trcnam_pisces  ***
   !! TOP :   initialisation of some run parameters for PISCES bio-model
   !!======================================================================
   !! History :    -   !  1999-10 (M.A. Foujols, M. Levy) original code
   !!              -   !  2000-01 (L. Bopp) hamocc3, p3zd
   !!             1.0  !  2003-08 (C. Ethe)  module F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) from trcnam.pisces.h90
   !!----------------------------------------------------------------------
#if defined key_pisces || defined key_pisces_reduced
   !!----------------------------------------------------------------------
   !!   'key_pisces'   :                                   PISCES bio-model
   !!----------------------------------------------------------------------
   !! trc_nam_pisces       : PISCES model namelist read
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE sms_pisces      ! sms trends
   USE trdmod_trc_oce
   USE iom             ! I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_pisces   ! called by trcnam.F90 module


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_pisces.F90 3680 2012-11-27 14:42:24Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_pisces
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_pisces  ***  
      !!
      !! ** Purpose :   read PISCES namelist
      !!
      !! ** input   :   file 'namelist.trc.sms' containing the following
      !!             namelist: natext, natbio, natsms
      !!                       natkriest ("key_kriest")
      !!----------------------------------------------------------------------
      !!
      INTEGER :: jl, jn
      TYPE(DIAG), DIMENSION(jp_pisces_2d)  :: pisdia2d
      TYPE(DIAG), DIMENSION(jp_pisces_3d)  :: pisdia3d
      TYPE(DIAG), DIMENSION(jp_pisces_trd) :: pisdiabio
      CHARACTER(LEN=20)   ::   clname
      !!
      NAMELIST/nampisdia/ pisdia3d, pisdia2d     ! additional diagnostics
#if defined key_pisces_reduced
      NAMELIST/nampisdbi/ pisdiabio
#endif

      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      clname = 'namelist_pisces'
#if defined key_pisces
      IF(lwp) WRITE(numout,*) ' trc_nam_pisces : read PISCES namelist'
#else
      IF(lwp) WRITE(numout,*) ' trc_nam_pisces : read LOBSTER namelist'
#endif
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      CALL ctl_opn( numnatp, TRIM( clname ), 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )

      !
      IF( .NOT.lk_iomput .AND. ln_diatrc ) THEN
         !
         ! Namelist nampisdia
         ! -------------------
         DO jl = 1, jp_pisces_2d
            WRITE(pisdia2d(jl)%sname,'("2D_",I1)') jl                      ! short name
            WRITE(pisdia2d(jl)%lname,'("2D DIAGNOSTIC NUMBER ",I2)') jl    ! long name
            pisdia2d(jl)%units = ' '                                       ! units
         END DO
         !                                 ! 3D output arrays
         DO jl = 1, jp_pisces_3d
            WRITE(pisdia3d(jl)%sname,'("3D_",I1)') jl                      ! short name
            WRITE(pisdia3d(jl)%lname,'("3D DIAGNOSTIC NUMBER ",I2)') jl    ! long name
            pisdia3d(jl)%units = ' '                                       ! units
         END DO

         REWIND( numnatp )               ! 
         READ  ( numnatp, nampisdia )

         DO jl = 1, jp_pisces_2d
            jn = jp_pcs0_2d + jl - 1
            ctrc2d(jn) = pisdia2d(jl)%sname
            ctrc2l(jn) = pisdia2d(jl)%lname
            ctrc2u(jn) = pisdia2d(jl)%units
         END DO

         DO jl = 1, jp_pisces_3d
            jn = jp_pcs0_3d + jl - 1
            ctrc3d(jn) = pisdia3d(jl)%sname
            ctrc3l(jn) = pisdia3d(jl)%lname
            ctrc3u(jn) = pisdia3d(jl)%units
         END DO

         IF(lwp) THEN                   ! control print
            WRITE(numout,*)
            WRITE(numout,*) ' Namelist : natadd'
            DO jl = 1, jp_pisces_3d
               jn = jp_pcs0_3d + jl - 1
               WRITE(numout,*) '  3d diag nb : ', jn, '    short name : ', ctrc3d(jn), &
                 &             '  long name  : ', ctrc3l(jn), '   unit : ', ctrc3u(jn)
            END DO
            WRITE(numout,*) ' '

            DO jl = 1, jp_pisces_2d
               jn = jp_pcs0_2d + jl - 1
               WRITE(numout,*) '  2d diag nb : ', jn, '    short name : ', ctrc2d(jn), &
                 &             '  long name  : ', ctrc2l(jn), '   unit : ', ctrc2u(jn)
            END DO
            WRITE(numout,*) ' '
         ENDIF
         !
      ENDIF

#if defined key_pisces_reduced

      IF( ( .NOT.lk_iomput .AND. ln_diabio ) .OR. lk_trdmld_trc ) THEN
         !
         ! Namelist nampisdbi
         ! -------------------
         DO jl = 1, jp_pisces_trd
            IF(     jl <  10 ) THEN   ;   WRITE (pisdiabio(jl)%sname,'("BIO_",I1)') jl      ! short name
            ELSEIF (jl < 100 ) THEN   ;   WRITE (pisdiabio(jl)%sname,'("BIO_",I2)') jl
            ELSE                      ;   WRITE (pisdiabio(jl)%sname,'("BIO_",I3)') jl
            ENDIF
            WRITE(pisdiabio(jl)%lname,'("BIOLOGICAL TREND NUMBER ",I2)') jl                 ! long name
            pisdiabio(jl)%units = 'mmoleN/m3/s '                                            ! units
         END DO

         REWIND( numnatp )
         READ  ( numnatp, nampisdbi )

         DO jl = 1, jp_pisces_trd
            jn = jp_pcs0_trd + jl - 1
            ctrbio(jl) = pisdiabio(jl)%sname
            ctrbil(jl) = pisdiabio(jl)%lname
            ctrbiu(jl) = pisdiabio(jl)%units
         END DO

         IF(lwp) THEN                   ! control print
            WRITE(numout,*)
            WRITE(numout,*) ' Namelist : nampisdbi'
            DO jl = 1, jp_pisces_trd
               jn = jp_pcs0_trd + jl - 1
               WRITE(numout,*) '  biological trend No : ', jn, '    short name : ', ctrbio(jn), &
                 &             '  long name  : ', ctrbio(jn), '   unit : ', ctrbio(jn)
            END DO
            WRITE(numout,*) ' '
         END IF
         !
      END IF

#endif

   END SUBROUTINE trc_nam_pisces

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No PISCES bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_pisces                      ! Empty routine
   END  SUBROUTINE  trc_nam_pisces
#endif  

   !!======================================================================
END MODULE trcnam_pisces
