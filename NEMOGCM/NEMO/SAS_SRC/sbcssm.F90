MODULE sbcssm
   !!======================================================================
   !!                       ***  MODULE  sbcssm  ***
   !! Off-line : interpolation of the physical fields
   !!======================================================================
   !! History : 
   !!   NEMO         3.4  ! 2012-03 First version by S. Alderson 
   !!                     !         Heavily derived from Christian's dtadyn routine
   !!                     !         in OFF_SRC
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_ssm_init : initialization, namelist read, and SAVEs control
   !!   sbc_ssm      : Interpolation of the fields
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE c1d             ! 1D configuration: lk_c1d
   USE dom_oce         ! ocean domain: variables
   USE zdf_oce         ! ocean vertical physics: variables
   USE sbc_oce         ! surface module: variables
   USE phycst          ! physical constants
   USE eosbn2          ! equation of state - Brunt Vaisala frequency
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE zpshde          ! z-coord. with partial steps: horizontal derivatives
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE lib_mpp         ! distributed memory computing library
   USE prtctl          ! print control
   USE fldread         ! read input fields 
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ssm_init   ! called by sbc_init
   PUBLIC   sbc_ssm        ! called by sbc

   CHARACTER(len=100)   ::   cn_dir     = './'    !: Root directory for location of ssm files
   LOGICAL              ::   ln_3d_uv   = .true.  !: specify whether input velocity data is 3D
   INTEGER  , SAVE      ::   nfld_3d
   INTEGER  , SAVE      ::   nfld_2d

   INTEGER  , PARAMETER ::   jpfld_3d = 4   ! maximum number of files to read
   INTEGER  , PARAMETER ::   jpfld_2d = 1   ! maximum number of files to read
   INTEGER  , SAVE      ::   jf_tem         ! index of temperature
   INTEGER  , SAVE      ::   jf_sal         ! index of salinity
   INTEGER  , SAVE      ::   jf_usp         ! index of u velocity component
   INTEGER  , SAVE      ::   jf_vsp         ! index of v velocity component
   INTEGER  , SAVE      ::   jf_ssh         ! index of sea surface height

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_ssm_3d  ! structure of input fields (file information, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_ssm_2d  ! structure of input fields (file information, fields read)

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OFF 3.3 , NEMO Consortium (2010)
   !! $Id: sbcssm.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_ssm( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ssm  ***
      !!
      !! ** Purpose :  Prepares dynamics and physics fields from a NEMO run
      !!               for an off-line simulation using surface processes only
      !!
      !! ** Method : calculates the position of data 
      !!             - interpolates data if needed
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) ::   ztinta     ! ratio applied to after  records when doing time interpolation
      REAL(wp) ::   ztintb     ! ratio applied to before records when doing time interpolation
      !!----------------------------------------------------------------------
      
      !
      IF( nn_timing == 1 )  CALL timing_start( 'sbc_ssm')

      IF( nfld_3d > 0 ) CALL fld_read( kt, 1, sf_ssm_3d )      !==   read data at kt time step   ==!
      IF( nfld_2d > 0 ) CALL fld_read( kt, 1, sf_ssm_2d )      !==   read data at kt time step   ==!
      ! 
      IF( ln_3d_uv ) THEN
         ssu_m(:,:) = sf_ssm_3d(jf_usp)%fnow(:,:,1) * umask(:,:,1)    ! u-velocity
         ssv_m(:,:) = sf_ssm_3d(jf_vsp)%fnow(:,:,1) * vmask(:,:,1)    ! v-velocity 
      ELSE
         ssu_m(:,:) = sf_ssm_2d(jf_usp)%fnow(:,:,1) * umask(:,:,1)    ! u-velocity
         ssv_m(:,:) = sf_ssm_2d(jf_vsp)%fnow(:,:,1) * vmask(:,:,1)    ! v-velocity 
      ENDIF
      !
      sst_m(:,:) = sf_ssm_2d(jf_tem)%fnow(:,:,1) * tmask(:,:,1)    ! temperature
      sss_m(:,:) = sf_ssm_2d(jf_sal)%fnow(:,:,1) * tmask(:,:,1)    ! salinity
      ssh_m(:,:) = sf_ssm_2d(jf_ssh)%fnow(:,:,1) * tmask(:,:,1)    ! sea surface height
      !
      tsn(:,:,1,jp_tem) = sst_m(:,:)
      tsn(:,:,1,jp_sal) = sss_m(:,:)
      IF ( nn_ice == 1 ) THEN
         tsb(:,:,1,jp_tem) = sst_m(:,:)
         tsb(:,:,1,jp_sal) = sss_m(:,:)
      ENDIF
      ub (:,:,1       ) = ssu_m(:,:)
      vb (:,:,1       ) = ssv_m(:,:)

      IF(ln_ctl) THEN                  ! print control
         CALL prt_ctl(tab2d_1=sst_m, clinfo1=' sst_m   - : ', mask1=tmask, ovlap=1   )
         CALL prt_ctl(tab2d_1=sss_m, clinfo1=' sss_m   - : ', mask1=tmask, ovlap=1   )
         CALL prt_ctl(tab2d_1=ssu_m, clinfo1=' ssu_m   - : ', mask1=umask, ovlap=1   )
         CALL prt_ctl(tab2d_1=ssv_m, clinfo1=' ssv_m   - : ', mask1=vmask, ovlap=1   )
         CALL prt_ctl(tab2d_1=ssh_m, clinfo1=' ssh_m   - : ', mask1=tmask, ovlap=1   )
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'sbc_ssm')
      !
   END SUBROUTINE sbc_ssm


   SUBROUTINE sbc_ssm_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ssm_init  ***
      !!
      !! ** Purpose :   Initialisation of the dynamical data     
      !! ** Method  : - read the data namsbc_ssm namelist
      !!
      !! ** Action  : - read parameters
      !!----------------------------------------------------------------------
      INTEGER  :: ierr, ierr0, ierr1, ierr2, ierr3   ! return error code
      INTEGER  :: ifpr                               ! dummy loop indice
      INTEGER  :: inum, idv, idimv, jpm              ! local integer
      !!
      CHARACTER(len=100)                     ::  cn_dir       ! Root directory for location of core files
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::  slf_3d       ! array of namelist information on the fields to read
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::  slf_2d       ! array of namelist information on the fields to read
      TYPE(FLD_N) :: sn_tem, sn_sal                     ! information about the fields to be read
      TYPE(FLD_N) :: sn_usp, sn_vsp, sn_ssh
      !
      NAMELIST/namsbc_sas/cn_dir, ln_3d_uv, sn_tem, sn_sal, sn_usp, sn_vsp, sn_ssh

      !!----------------------------------------------------------------------
      !                                   ! ============
      !                                   !   Namelist
      !                                   ! ============
      ! (NB: frequency positive => hours, negative => months)
      !                !   file      ! frequency !  variable  ! time intep !  clim  ! 'yearly' or ! weights  ! rotation   !
      !                !   name      !  (hours)  !   name     !   (T/F)    !  (T/F) !  'monthly'  ! filename ! pairs      !
      sn_usp  = FLD_N( 'ssm_grid_U' ,    120    , 'vozocrtx' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_vsp  = FLD_N( 'ssm_grid_V' ,    120    , 'vomecrty' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_tem  = FLD_N( 'ssm_grid_T' ,    120    , 'sosstsst' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_sal  = FLD_N( 'ssm_grid_T' ,    120    , 'sosaline' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      sn_ssh  = FLD_N( 'ssm_grid_T' ,    120    , 'sossheig' ,  .true.    , .true. ,   'yearly'  , ''       , ''         )
      !
      REWIND( numnam )                          ! read in namlist namsbc_ssm
      READ  ( numnam, namsbc_sas )
      !                                         ! store namelist information in an array
      !                                         ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_sas : standalone surface scheme '
         WRITE(numout,*) '~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namsbc_sas'
         WRITE(numout,*)
      ENDIF
      
      !
      !! switch off stuff that isn't sensible with a standalone module
      !! note that we need sbc_ssm called first in sbc
      !
      IF( ln_cpl ) THEN
         IF( lwp ) WRITE(numout,*) 'Coupled mode not sensible with StandAlone Surface scheme'
         ln_cpl = .FALSE.
      ENDIF
      IF( ln_apr_dyn ) THEN
         IF( lwp ) WRITE(numout,*) 'No atmospheric gradient needed with StandAlone Surface scheme'
         ln_apr_dyn = .FALSE.
      ENDIF
      IF( ln_dm2dc ) THEN
         IF( lwp ) WRITE(numout,*) 'No diurnal cycle needed with StandAlone Surface scheme'
         ln_dm2dc = .FALSE.
      ENDIF
      IF( ln_rnf ) THEN
         IF( lwp ) WRITE(numout,*) 'No runoff needed with StandAlone Surface scheme'
         ln_rnf = .FALSE.
      ENDIF
      IF( ln_ssr ) THEN
         IF( lwp ) WRITE(numout,*) 'No surface relaxation needed with StandAlone Surface scheme'
         ln_ssr = .FALSE.
      ENDIF
      IF( nn_fwb > 0 ) THEN
         IF( lwp ) WRITE(numout,*) 'No freshwater budget adjustment needed with StandAlone Surface scheme'
         nn_fwb = 0
      ENDIF
      IF( nn_closea > 0 ) THEN
         IF( lwp ) WRITE(numout,*) 'No closed seas adjustment needed with StandAlone Surface scheme'
         nn_closea = 0
      ENDIF

      ! 
      !! following code is a bit messy, but distinguishes between when u,v are 3d arrays and
      !! when we have other 3d arrays that we need to read in
      !! so if a new field is added i.e. jf_new, just give it the next integer in sequence
      !! for the corresponding dimension (currently if ln_3d_uv is true, 4 for 2d and 3 for 3d,
      !! alternatively if ln_3d_uv is false, 6 for 2d and 1 for 3d), reset nfld_3d, nfld_2d,
      !! and the rest of the logic should still work
      !
      jf_tem = 1 ; jf_sal = 2 ; jf_ssh = 3
      !
      IF( ln_3d_uv ) THEN
         jf_usp = 1 ; jf_vsp = 2
         nfld_3d  = 2
         nfld_2d  = 3
      ELSE
         jf_usp = 4 ; jf_vsp = 5
         nfld_3d  = 0
         nfld_2d  = 5
      ENDIF

      IF( nfld_3d > 0 ) THEN
         ALLOCATE( slf_3d(nfld_3d), STAT=ierr )         ! set slf structure
         IF( ierr > 0 ) THEN
            CALL ctl_stop( 'sbc_ssm_init: unable to allocate slf 3d structure' )   ;   RETURN
         ENDIF
         IF( ln_3d_uv ) THEN
            slf_3d(jf_usp) = sn_usp
            slf_3d(jf_vsp) = sn_vsp
         ENDIF
      ENDIF

      IF( nfld_2d > 0 ) THEN
         ALLOCATE( slf_2d(nfld_2d), STAT=ierr )         ! set slf structure
         IF( ierr > 0 ) THEN
            CALL ctl_stop( 'sbc_ssm_init: unable to allocate slf 2d structure' )   ;   RETURN
         ENDIF
         slf_2d(jf_tem) = sn_tem ; slf_2d(jf_sal) = sn_sal ; slf_2d(jf_ssh) = sn_ssh
         IF( .NOT. ln_3d_uv ) THEN
            slf_2d(jf_usp) = sn_usp ; slf_2d(jf_vsp) = sn_vsp
         ENDIF
      ENDIF
      !
      IF( nfld_3d > 0 ) THEN
         ALLOCATE( sf_ssm_3d(nfld_3d), STAT=ierr )         ! set sf structure
         IF( ierr > 0 ) THEN
            CALL ctl_stop( 'sbc_ssm_init: unable to allocate sf structure' )   ;   RETURN
         ENDIF
         DO ifpr = 1, nfld_3d
                                       ALLOCATE( sf_ssm_3d(ifpr)%fnow(jpi,jpj,jpk)    , STAT=ierr0 )
            IF( slf_3d(ifpr)%ln_tint ) ALLOCATE( sf_ssm_3d(ifpr)%fdta(jpi,jpj,jpk,2)  , STAT=ierr1 )
            IF( ierr0 + ierr1 > 0 ) THEN
               CALL ctl_stop( 'sbc_ssm_init : unable to allocate sf_ssm_3d array structure' )   ;   RETURN
            ENDIF
         END DO
         !                                         ! fill sf with slf_i and control print
         CALL fld_fill( sf_ssm_3d, slf_3d, cn_dir, 'sbc_ssm_init', '3D Data in file', 'namsbc_ssm' )
      ENDIF

      IF( nfld_2d > 0 ) THEN
         ALLOCATE( sf_ssm_2d(nfld_2d), STAT=ierr )         ! set sf structure
         IF( ierr > 0 ) THEN
            CALL ctl_stop( 'sbc_ssm_init: unable to allocate sf 2d structure' )   ;   RETURN
         ENDIF
         DO ifpr = 1, nfld_2d
                                       ALLOCATE( sf_ssm_2d(ifpr)%fnow(jpi,jpj,1)    , STAT=ierr0 )
            IF( slf_2d(ifpr)%ln_tint ) ALLOCATE( sf_ssm_2d(ifpr)%fdta(jpi,jpj,1,2)  , STAT=ierr1 )
            IF( ierr0 + ierr1 > 0 ) THEN
               CALL ctl_stop( 'sbc_ssm_init : unable to allocate sf_ssm_2d array structure' )   ;   RETURN
            ENDIF
         END DO
         !
         CALL fld_fill( sf_ssm_2d, slf_2d, cn_dir, 'sbc_ssm_init', '2D Data in file', 'namsbc_ssm' )
      ENDIF
      !
      ! lim code currently uses surface temperature and salinity in tsn array for initialisation
      ! and ub, vb arrays in ice dynamics
      ! so allocate enough of arrays to use
      !
      ierr3 = 0
      jpm = MAX(jp_tem, jp_sal)
      ALLOCATE( tsn(jpi,jpj,1,jpm), STAT=ierr0 )
      ALLOCATE( ub(jpi,jpj,1)     , STAT=ierr1 )
      ALLOCATE( vb(jpi,jpj,1)     , STAT=ierr2 )
      IF ( nn_ice == 1 ) ALLOCATE( tsb(jpi,jpj,1,jpm), STAT=ierr3 )
      ierr = ierr0 + ierr1 + ierr2 + ierr3
      IF( ierr > 0 ) THEN
         CALL ctl_stop('sbc_ssm_init: unable to allocate surface arrays')
      ENDIF
      !
      ! finally tidy up

      IF( nfld_3d > 0 ) DEALLOCATE( slf_3d, STAT=ierr )
      IF( nfld_2d > 0 ) DEALLOCATE( slf_2d, STAT=ierr )
      !
   END SUBROUTINE sbc_ssm_init

   !!======================================================================
END MODULE sbcssm
