MODULE iom
   !!=====================================================================
   !!                    ***  MODULE  iom ***
   !! Input/Output manager :  Library to read input files
   !!====================================================================
   !! History :  2.0  ! 2005-12  (J. Belier) Original code
   !!            2.0  ! 2006-02  (S. Masson) Adaptation to NEMO
   !!            3.0  ! 2007-07  (D. Storkey) Changes to iom_gettime
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie and G. Reffray)  add C1D case  
   !!--------------------------------------------------------------------

   !!--------------------------------------------------------------------
   !!   iom_open       : open a file read only
   !!   iom_close      : close a file or all files opened by iom
   !!   iom_get        : read a field (interfaced to several routines)
   !!   iom_gettime    : read the time axis cdvar in the file
   !!   iom_varid      : get the id of a variable in a file
   !!   iom_rstput     : write a field in a restart file (interfaced to several routines)
   !!--------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE c1d             ! 1D vertical configuration
   USE flo_oce         ! floats module declarations
   USE lbclnk          ! lateal boundary condition / mpp exchanges
   USE iom_def         ! iom variables definitions
   USE iom_ioipsl      ! NetCDF format with IOIPSL library
   USE iom_nf90        ! NetCDF format with native NetCDF library
   USE iom_rstdimg     ! restarts access direct format "dimg" style...
   USE in_out_manager  ! I/O manager
   USE lib_mpp           ! MPP library
#if defined key_iomput
   USE sbc_oce, ONLY :   nn_fsbc                ! ocean space and time domain
   USE trc_oce, ONLY :   nn_dttrc               ! frequency of step on passive tracers
   USE icb_oce, ONLY :   nclasses, class_num    ! iceberg classes
   USE domngb          ! ocean space and time domain
   USE phycst          ! physical constants
   USE dianam          ! build name of file
   USE xios
# endif
   USE ioipsl, ONLY :  ju2ymds    ! for calendar

   IMPLICIT NONE
   PUBLIC   !   must be public to be able to access iom_def through iom
   
#if defined key_iomput
   LOGICAL, PUBLIC, PARAMETER ::   lk_iomput = .TRUE.        !: iom_put flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_iomput = .FALSE.       !: iom_put flag
#endif
   PUBLIC iom_init, iom_swap, iom_open, iom_close, iom_setkt, iom_varid, iom_get, iom_gettime, iom_rstput, iom_put
   PUBLIC iom_getatt

   PRIVATE iom_rp0d, iom_rp1d, iom_rp2d, iom_rp3d
   PRIVATE iom_g0d, iom_g1d, iom_g2d, iom_g3d, iom_get_123d
   PRIVATE iom_p1d, iom_p2d, iom_p3d
#if defined key_iomput
   PRIVATE iom_set_domain_attr, iom_set_axis_attr, iom_set_field_attr, iom_set_file_attr, iom_get_file_attr, iom_set_grid_attr
   PRIVATE set_grid, set_scalar, set_xmlatt, set_mooring, iom_update_file_name, iom_sdate
# endif

   INTERFACE iom_get
      MODULE PROCEDURE iom_g0d, iom_g1d, iom_g2d, iom_g3d
   END INTERFACE
   INTERFACE iom_getatt
      MODULE PROCEDURE iom_g0d_intatt
   END INTERFACE
   INTERFACE iom_rstput
      MODULE PROCEDURE iom_rp0d, iom_rp1d, iom_rp2d, iom_rp3d
   END INTERFACE
  INTERFACE iom_put
     MODULE PROCEDURE iom_p0d, iom_p1d, iom_p2d, iom_p3d
  END INTERFACE
#if defined key_iomput
   INTERFACE iom_setkt
      MODULE PROCEDURE xios_update_calendar
   END INTERFACE
# endif

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: iom.F90 4191 2013-11-14 10:05:09Z cetlod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE iom_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE   ***
      !!
      !! ** Purpose :   
      !!
      !!----------------------------------------------------------------------
#if defined key_iomput
      TYPE(xios_time)   :: dtime    = xios_time(0, 0, 0, 0, 0, 0)
      CHARACTER(len=19) :: cldate 
      CHARACTER(len=10) :: clname
      INTEGER           ::   ji
      !!----------------------------------------------------------------------

      clname = "nemo"
      IF( TRIM(Agrif_CFixed()) /= '0' )   clname = TRIM(Agrif_CFixed())//"_"//TRIM(clname)
      CALL xios_context_initialize(TRIM(clname), mpi_comm_opa)
      CALL iom_swap

      ! calendar parameters
      SELECT CASE ( nleapy )        ! Choose calendar for IOIPSL
      CASE ( 1)   ;   CALL xios_set_context_attr(TRIM(clname), calendar_type= "Gregorian")
      CASE ( 0)   ;   CALL xios_set_context_attr(TRIM(clname), calendar_type= "NoLeap")
      CASE (30)   ;   CALL xios_set_context_attr(TRIM(clname), calendar_type= "D360")
      END SELECT
      WRITE(cldate,"(i4.4,'-',i2.2,'-',i2.2,' 00:00:00')") nyear,nmonth,nday 
      CALL xios_set_context_attr(TRIM(clname), start_date=cldate )

      ! horizontal grid definition
      CALL set_scalar
      CALL set_grid( "T", glamt, gphit ) 
      CALL set_grid( "U", glamu, gphiu )
      CALL set_grid( "V", glamv, gphiv )
      CALL set_grid( "W", glamt, gphit )

      ! vertical grid definition
      CALL iom_set_axis_attr( "deptht", gdept_0 )
      CALL iom_set_axis_attr( "depthu", gdept_0 )
      CALL iom_set_axis_attr( "depthv", gdept_0 )
      CALL iom_set_axis_attr( "depthw", gdepw_0 )
# if defined key_floats
      CALL iom_set_axis_attr( "nfloat", (/ (REAL(ji,wp), ji=1,nfloat) /) )
# endif
      CALL iom_set_axis_attr( "icbcla", class_num )
      
      ! automatic definitions of some of the xml attributs
      CALL set_xmlatt

      ! end file definition
      dtime%second = rdt
      CALL xios_set_timestep(dtime)
      CALL xios_close_context_definition()
      
      CALL xios_update_calendar(0)
#endif
      
   END SUBROUTINE iom_init


   SUBROUTINE iom_swap
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_swap  ***
      !!
      !! ** Purpose :  swap context between different agrif grid for xmlio_server
      !!---------------------------------------------------------------------
#if defined key_iomput
      TYPE(xios_context) :: nemo_hdl

     IF( TRIM(Agrif_CFixed()) == '0' ) THEN
        CALL xios_get_handle("nemo",nemo_hdl)
     ELSE
        CALL xios_get_handle(TRIM(Agrif_CFixed())//"_nemo",nemo_hdl)
     ENDIF
     CALL xios_set_current_context(nemo_hdl)

#endif
   END SUBROUTINE iom_swap


   SUBROUTINE iom_open( cdname, kiomid, ldwrt, kdom, kiolib, ldstop, ldiof )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_open  ***
      !!
      !! ** Purpose :  open an input file (return 0 if not found)
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   )           ::   cdname   ! File name
      INTEGER         , INTENT(  out)           ::   kiomid   ! iom identifier of the opened file
      LOGICAL         , INTENT(in   ), OPTIONAL ::   ldwrt    ! open in write modeb          (default = .FALSE.)
      INTEGER         , INTENT(in   ), OPTIONAL ::   kdom     ! Type of domain to be written (default = jpdom_local_noovlap)
      INTEGER         , INTENT(in   ), OPTIONAL ::   kiolib   ! library used to open the file (default = jpnf90) 
      LOGICAL         , INTENT(in   ), OPTIONAL ::   ldstop   ! stop if open to read a non-existing file (default = .TRUE.)
      LOGICAL         , INTENT(in   ), OPTIONAL ::   ldiof    ! Interp On the Fly, needed for AGRIF (default = .FALSE.)

      CHARACTER(LEN=256)    ::   clname    ! the name of the file based on cdname [[+clcpu]+clcpu]
      CHARACTER(LEN=256)    ::   cltmpn    ! tempory name to store clname (in writting mode)
      CHARACTER(LEN=10)     ::   clsuffix  ! ".nc" or ".dimg"
      CHARACTER(LEN=15)     ::   clcpu     ! the cpu number (max jpmax_digits digits)
      CHARACTER(LEN=256)    ::   clinfo    ! info character
      LOGICAL               ::   llok      ! check the existence 
      LOGICAL               ::   llwrt     ! local definition of ldwrt
      LOGICAL               ::   llnoov    ! local definition to read overlap
      LOGICAL               ::   llstop    ! local definition of ldstop
      LOGICAL               ::   lliof     ! local definition of ldiof
      INTEGER               ::   iolib     ! library do we use to open the file
      INTEGER               ::   icnt      ! counter for digits in clcpu (max = jpmax_digits)
      INTEGER               ::   iln, ils  ! lengths of character
      INTEGER               ::   idom      ! type of domain
      INTEGER               ::   istop     ! 
      INTEGER, DIMENSION(2,5) ::   idompar ! domain parameters: 
      ! local number of points for x,y dimensions
      ! position of first local point for x,y dimensions
      ! position of last local point for x,y dimensions
      ! start halo size for x,y dimensions
      ! end halo size for x,y dimensions
      !---------------------------------------------------------------------
      ! Initializations and control
      ! =============
      kiomid = -1
      clinfo = '                    iom_open ~~~  '
      istop = nstop
      ! if iom_open is called for the first time: initialize iom_file(:)%nfid to 0
      ! (could be done when defining iom_file in f95 but not in f90)
      IF( Agrif_Root() ) THEN
         IF( iom_open_init == 0 ) THEN
            iom_file(:)%nfid = 0
            iom_open_init = 1
         ENDIF
      ENDIF
      ! do we read or write the file?
      IF( PRESENT(ldwrt) ) THEN   ;   llwrt = ldwrt
      ELSE                        ;   llwrt = .FALSE.
      ENDIF
      ! do we call ctl_stop if we try to open a non-existing file in read mode?
      IF( PRESENT(ldstop) ) THEN   ;   llstop = ldstop
      ELSE                         ;   llstop = .TRUE.
      ENDIF
      ! what library do we use to open the file?
      IF( PRESENT(kiolib) ) THEN   ;   iolib = kiolib
      ELSE                         ;   iolib = jpnf90
      ENDIF
      ! are we using interpolation on the fly?
      IF( PRESENT(ldiof) ) THEN   ;   lliof = ldiof
      ELSE                        ;   lliof = .FALSE.
      ENDIF
      ! do we read the overlap 
      ! ugly patch SM+JMM+RB to overwrite global definition in some cases
      llnoov = (jpni * jpnj ) == jpnij .AND. .NOT. lk_agrif 
      ! create the file name by added, if needed, TRIM(Agrif_CFixed()) and TRIM(clsuffix)
      ! =============
      clname   = trim(cdname)
      IF ( .NOT. Agrif_Root() .AND. .NOT. lliof ) THEN
         iln    = INDEX(clname,'/') 
         cltmpn = clname(1:iln)
         clname = clname(iln+1:LEN_TRIM(clname))
         clname=TRIM(cltmpn)//TRIM(Agrif_CFixed())//'_'//TRIM(clname)
      ENDIF
      ! which suffix should we use?
      SELECT CASE (iolib)
      CASE (jpioipsl ) ;   clsuffix = '.nc'
      CASE (jpnf90   ) ;   clsuffix = '.nc'
      CASE (jprstdimg) ;   clsuffix = '.dimg'
      CASE DEFAULT     ;   clsuffix = ''
         CALL ctl_stop( TRIM(clinfo), 'accepted IO library are only jpioipsl, jpnf90 and jprstdimg' )
      END SELECT
      ! Add the suffix if needed
      iln = LEN_TRIM(clname)
      ils = LEN_TRIM(clsuffix)
      IF( iln <= ils .OR. INDEX( TRIM(clname), TRIM(clsuffix), back = .TRUE. ) /= iln - ils + 1 )   &
         &   clname = TRIM(clname)//TRIM(clsuffix)
      cltmpn = clname   ! store this name
      ! try to find if the file to be opened already exist
      ! =============
      INQUIRE( FILE = clname, EXIST = llok )
      IF( .NOT.llok ) THEN
         ! we try to add the cpu number to the name
         IF( iolib == jprstdimg ) THEN   ;   WRITE(clcpu,*) narea
         ELSE                            ;   WRITE(clcpu,*) narea-1
         ENDIF
         clcpu  = TRIM(ADJUSTL(clcpu))
         iln = INDEX(clname,TRIM(clsuffix), back = .TRUE.)
         clname = clname(1:iln-1)//'_'//TRIM(clcpu)//TRIM(clsuffix)
         icnt = 0
         INQUIRE( FILE = clname, EXIST = llok ) 
         ! we try different formats for the cpu number by adding 0
         DO WHILE( .NOT.llok .AND. icnt < jpmax_digits )
            clcpu  = "0"//trim(clcpu)
            clname = clname(1:iln-1)//'_'//TRIM(clcpu)//TRIM(clsuffix)
            INQUIRE( FILE = clname, EXIST = llok )
            icnt = icnt + 1
         END DO
      ENDIF
      IF( llwrt ) THEN
         ! check the domain definition
! JMM + SM: ugly patch before getting the new version of lib_mpp)
!         idom = jpdom_local_noovlap   ! default definition
         IF( llnoov ) THEN   ;   idom = jpdom_local_noovlap   ! default definition
         ELSE                ;   idom = jpdom_local_full      ! default definition
         ENDIF
         IF( PRESENT(kdom) )   idom = kdom
         ! create the domain informations
         ! =============
         SELECT CASE (idom)
         CASE (jpdom_local_full)
            idompar(:,1) = (/ jpi             , jpj              /)
            idompar(:,2) = (/ nimpp           , njmpp            /)
            idompar(:,3) = (/ nimpp + jpi - 1 , njmpp + jpj - 1  /)
            idompar(:,4) = (/ nldi - 1        , nldj - 1         /)
            idompar(:,5) = (/ jpi - nlei      , jpj - nlej       /)
         CASE (jpdom_local_noextra)
            idompar(:,1) = (/ nlci            , nlcj             /)
            idompar(:,2) = (/ nimpp           , njmpp            /)
            idompar(:,3) = (/ nimpp + nlci - 1, njmpp + nlcj - 1 /)
            idompar(:,4) = (/ nldi - 1        , nldj - 1         /)
            idompar(:,5) = (/ nlci - nlei     , nlcj - nlej      /)
         CASE (jpdom_local_noovlap)
            idompar(:,1) = (/ nlei  - nldi + 1, nlej  - nldj + 1 /)
            idompar(:,2) = (/ nimpp + nldi - 1, njmpp + nldj - 1 /)
            idompar(:,3) = (/ nimpp + nlei - 1, njmpp + nlej - 1 /)
            idompar(:,4) = (/ 0               , 0                /)
            idompar(:,5) = (/ 0               , 0                /)
         CASE DEFAULT
            CALL ctl_stop( TRIM(clinfo), 'wrong value of kdom, only jpdom_local* cases are accepted' )
         END SELECT
      ENDIF
      ! Open the NetCDF or RSTDIMG file
      ! =============
      ! do we have some free file identifier?
      IF( MINVAL(iom_file(:)%nfid) /= 0 )   &
         &   CALL ctl_stop( TRIM(clinfo), 'No more free file identifier', 'increase jpmax_files in iom_def' )
      ! if no file was found...
      IF( .NOT. llok ) THEN
         IF( .NOT. llwrt ) THEN   ! we are in read mode 
            IF( llstop ) THEN   ;   CALL ctl_stop( TRIM(clinfo), 'File '//TRIM(cltmpn)//'* not found' )
            ELSE                ;   istop = nstop + 1   ! make sure that istop /= nstop so we don't open the file
            ENDIF
         ELSE                     ! we are in write mode so we 
            clname = cltmpn       ! get back the file name without the cpu number
         ENDIF
      ELSE
         IF( llwrt .AND. .NOT. ln_clobber ) THEN   ! we stop as we want to write in a new file 
            CALL ctl_stop( TRIM(clinfo), 'We want to write in a new file but '//TRIM(clname)//' already exists...' )
            istop = nstop + 1                      ! make sure that istop /= nstop so we don't open the file
         ENDIF
      ENDIF
      IF( istop == nstop ) THEN   ! no error within this routine
         SELECT CASE (iolib)
         CASE (jpioipsl )   ;   CALL iom_ioipsl_open(  clname, kiomid, llwrt, llok, idompar )
         CASE (jpnf90   )   ;   CALL iom_nf90_open(    clname, kiomid, llwrt, llok, idompar )
         CASE (jprstdimg)   ;   CALL iom_rstdimg_open( clname, kiomid, llwrt, llok, idompar )
         CASE DEFAULT
            CALL ctl_stop( TRIM(clinfo)//' accepted IO library are only jpioipsl, jpnf90 and jprstdimg' )
         END SELECT
      ENDIF
      !
   END SUBROUTINE iom_open


   SUBROUTINE iom_close( kiomid )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_close  ***
      !!
      !! ** Purpose : close an input file, or all files opened by iom
      !!--------------------------------------------------------------------
      INTEGER, INTENT(inout), OPTIONAL ::   kiomid   ! iom identifier of the file to be closed
      !                                              ! return 0 when file is properly closed
      !                                              ! No argument: all files opened by iom are closed

      INTEGER ::   jf         ! dummy loop indices
      INTEGER ::   i_s, i_e   ! temporary integer
      CHARACTER(LEN=100)    ::   clinfo    ! info character
      !---------------------------------------------------------------------
      !
      clinfo = '                    iom_close ~~~  '
      IF( PRESENT(kiomid) ) THEN
         i_s = kiomid
         i_e = kiomid
      ELSE
         i_s = 1
         i_e = jpmax_files
      ENDIF

      IF( i_s > 0 ) THEN
         DO jf = i_s, i_e
            IF( iom_file(jf)%nfid > 0 ) THEN
               SELECT CASE (iom_file(jf)%iolib)
               CASE (jpioipsl )   ;   CALL iom_ioipsl_close(  jf )
               CASE (jpnf90   )   ;   CALL iom_nf90_close(    jf )
               CASE (jprstdimg)   ;   CALL iom_rstdimg_close( jf )
               CASE DEFAULT
                  CALL ctl_stop( TRIM(clinfo)//' accepted IO library are only jpioipsl, jpnf90 and jprstdimg' )
               END SELECT
               iom_file(jf)%nfid       = 0          ! free the id 
               IF( PRESENT(kiomid) )   kiomid = 0   ! return 0 as id to specify that the file was closed
               IF(lwp) WRITE(numout,*) TRIM(clinfo)//' close file: '//TRIM(iom_file(jf)%name)//' ok'
            ELSEIF( PRESENT(kiomid) ) THEN
               WRITE(ctmp1,*) '--->',  kiomid
               CALL ctl_stop( TRIM(clinfo)//' Invalid file identifier', ctmp1 )
            ENDIF
         END DO
      ENDIF
      !    
   END SUBROUTINE iom_close


   FUNCTION iom_varid ( kiomid, cdvar, kdimsz, ldstop )  
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  iom_varid  ***
      !!
      !! ** Purpose : get the id of a variable in a file (return 0 if not found)
      !!-----------------------------------------------------------------------
      INTEGER              , INTENT(in   )           ::   kiomid   ! file Identifier
      CHARACTER(len=*)     , INTENT(in   )           ::   cdvar    ! name of the variable
      INTEGER, DIMENSION(:), INTENT(  out), OPTIONAL ::   kdimsz   ! size of the dimensions
      LOGICAL              , INTENT(in   ), OPTIONAL ::   ldstop   ! stop if looking for non-existing variable (default = .TRUE.)
      !
      INTEGER                        ::   iom_varid, iiv, i_nvd
      LOGICAL                        ::   ll_fnd
      CHARACTER(LEN=100)             ::   clinfo                   ! info character
      LOGICAL                        ::   llstop                   ! local definition of ldstop
      !!-----------------------------------------------------------------------
      iom_varid = 0                         ! default definition
      ! do we call ctl_stop if we look for non-existing variable?
      IF( PRESENT(ldstop) ) THEN   ;   llstop = ldstop
      ELSE                         ;   llstop = .TRUE.
      ENDIF
      !
      IF( kiomid > 0 ) THEN
         clinfo = 'iom_varid, file: '//trim(iom_file(kiomid)%name)//', var: '//trim(cdvar)
         IF( iom_file(kiomid)%nfid == 0 ) THEN 
            CALL ctl_stop( trim(clinfo), 'the file is not open' )
         ELSE
            ll_fnd  = .FALSE.
            iiv = 0
            !
            DO WHILE ( .NOT.ll_fnd .AND. iiv < iom_file(kiomid)%nvars )
               iiv = iiv + 1
               ll_fnd  = ( TRIM(cdvar) == TRIM(iom_file(kiomid)%cn_var(iiv)) )
            END DO
            !
            IF( .NOT.ll_fnd ) THEN
               iiv = iiv + 1
               IF( iiv <= jpmax_vars ) THEN
                  SELECT CASE (iom_file(kiomid)%iolib)
                  CASE (jpioipsl )   ;   iom_varid = iom_ioipsl_varid( kiomid, cdvar, iiv, kdimsz )
                  CASE (jpnf90   )   ;   iom_varid = iom_nf90_varid  ( kiomid, cdvar, iiv, kdimsz )
                  CASE (jprstdimg)   ;   iom_varid = -1   ! all variables are listed in iom_file
                  CASE DEFAULT   
                     CALL ctl_stop( TRIM(clinfo)//' accepted IO library are only jpioipsl, jpnf90 and jprstdimg' )
                  END SELECT
               ELSE
                  CALL ctl_stop( trim(clinfo), 'Too many variables in the file '//iom_file(kiomid)%name,   &
                        &                         'increase the parameter jpmax_vars')
               ENDIF
               IF( llstop .AND. iom_varid == -1 )   CALL ctl_stop( TRIM(clinfo)//' not found' ) 
            ELSE
               iom_varid = iiv
               IF( PRESENT(kdimsz) ) THEN 
                  i_nvd = iom_file(kiomid)%ndims(iiv)
                  IF( i_nvd == size(kdimsz) ) THEN
                     kdimsz(:) = iom_file(kiomid)%dimsz(1:i_nvd,iiv)
                  ELSE
                     WRITE(ctmp1,*) i_nvd, size(kdimsz)
                     CALL ctl_stop( trim(clinfo), 'error in kdimsz size'//trim(ctmp1) )
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      !
   END FUNCTION iom_varid


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_get
   !!----------------------------------------------------------------------
   SUBROUTINE iom_g0d( kiomid, cdvar, pvar )
      INTEGER         , INTENT(in   )                 ::   kiomid    ! Identifier of the file
      CHARACTER(len=*), INTENT(in   )                 ::   cdvar     ! Name of the variable
      REAL(wp)        , INTENT(  out)                 ::   pvar      ! read field
      !
      INTEGER               :: idvar   ! variable id
      !
      IF( kiomid > 0 ) THEN
         idvar = iom_varid( kiomid, cdvar )
         IF( iom_file(kiomid)%nfid > 0 .AND. idvar > 0 ) THEN
            SELECT CASE (iom_file(kiomid)%iolib)
            CASE (jpioipsl )   ;   CALL iom_ioipsl_get(  kiomid, idvar, pvar )
            CASE (jpnf90   )   ;   CALL iom_nf90_get(    kiomid, idvar, pvar )
            CASE (jprstdimg)   ;   CALL iom_rstdimg_get( kiomid, idvar, pvar )
            CASE DEFAULT    
               CALL ctl_stop( 'iom_g0d: accepted IO library are only jpioipsl, jpnf90 and jprstdimg' )
            END SELECT
         ENDIF
      ENDIF
   END SUBROUTINE iom_g0d

   SUBROUTINE iom_g1d( kiomid, kdom, cdvar, pvar, ktime, kstart, kcount )
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(wp)        , INTENT(  out), DIMENSION(:)           ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      INTEGER         , INTENT(in   ), DIMENSION(1), OPTIONAL ::   kstart    ! start axis position of the reading 
      INTEGER         , INTENT(in   ), DIMENSION(1), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom       , cdvar        , pv_r1d=pvar,   &
              &                                                     ktime=ktime, kstart=kstart, kcount=kcount )
      ENDIF
   END SUBROUTINE iom_g1d

   SUBROUTINE iom_g2d( kiomid, kdom, cdvar, pvar, ktime, kstart, kcount )
      INTEGER         , INTENT(in   )                           ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                           ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                           ::   cdvar     ! Name of the variable
      REAL(wp)        , INTENT(  out), DIMENSION(:,:)           ::   pvar      ! read field
      INTEGER         , INTENT(in   )                , OPTIONAL ::   ktime     ! record number
      INTEGER         , INTENT(in   ), DIMENSION(2)  , OPTIONAL ::   kstart    ! start axis position of the reading 
      INTEGER         , INTENT(in   ), DIMENSION(2)  , OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom       , cdvar        , pv_r2d=pvar,   &
              &                                                     ktime=ktime, kstart=kstart, kcount=kcount )
      ENDIF
   END SUBROUTINE iom_g2d

   SUBROUTINE iom_g3d( kiomid, kdom, cdvar, pvar, ktime, kstart, kcount )
      INTEGER         , INTENT(in   )                             ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                             ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                             ::   cdvar     ! Name of the variable
      REAL(wp)        , INTENT(  out), DIMENSION(:,:,:)           ::   pvar      ! read field
      INTEGER         , INTENT(in   )                  , OPTIONAL ::   ktime     ! record number
      INTEGER         , INTENT(in   ), DIMENSION(3)    , OPTIONAL ::   kstart    ! start axis position of the reading 
      INTEGER         , INTENT(in   ), DIMENSION(3)    , OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom       , cdvar        , pv_r3d=pvar,   &
              &                                                     ktime=ktime, kstart=kstart, kcount=kcount )
      ENDIF
   END SUBROUTINE iom_g3d
   !!----------------------------------------------------------------------

   SUBROUTINE iom_get_123d( kiomid, kdom  , cdvar ,   &
         &                  pv_r1d, pv_r2d, pv_r3d,   &
         &                  ktime , kstart, kcount  )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_get_123d  ***
      !!
      !! ** Purpose : read a 1D/2D/3D variable
      !!
      !! ** Method : read ONE record at each CALL
      !!-----------------------------------------------------------------------
      INTEGER                    , INTENT(in   )           ::   kiomid     ! Identifier of the file
      INTEGER                    , INTENT(in   )           ::   kdom       ! Type of domain to be read
      CHARACTER(len=*)           , INTENT(in   )           ::   cdvar      ! Name of the variable
      REAL(wp), DIMENSION(:)     , INTENT(  out), OPTIONAL ::   pv_r1d     ! read field (1D case)
      REAL(wp), DIMENSION(:,:)   , INTENT(  out), OPTIONAL ::   pv_r2d     ! read field (2D case)
      REAL(wp), DIMENSION(:,:,:) , INTENT(  out), OPTIONAL ::   pv_r3d     ! read field (3D case)
      INTEGER                    , INTENT(in   ), OPTIONAL ::   ktime      ! record number
      INTEGER , DIMENSION(:)     , INTENT(in   ), OPTIONAL ::   kstart     ! start position of the reading in each axis 
      INTEGER , DIMENSION(:)     , INTENT(in   ), OPTIONAL ::   kcount     ! number of points to be read in each axis
      !
      LOGICAL                        ::   llnoov      ! local definition to read overlap
      INTEGER                        ::   jl          ! loop on number of dimension 
      INTEGER                        ::   idom        ! type of domain
      INTEGER                        ::   idvar       ! id of the variable
      INTEGER                        ::   inbdim      ! number of dimensions of the variable
      INTEGER                        ::   idmspc      ! number of spatial dimensions 
      INTEGER                        ::   itime       ! record number
      INTEGER                        ::   istop       ! temporary value of nstop
      INTEGER                        ::   ix1, ix2, iy1, iy2   ! subdomain indexes
      INTEGER                        ::   ji, jj      ! loop counters
      INTEGER                        ::   irankpv       ! 
      INTEGER                        ::   ind1, ind2  ! substring index
      INTEGER, DIMENSION(jpmax_dims) ::   istart      ! starting point to read for each axis
      INTEGER, DIMENSION(jpmax_dims) ::   icnt        ! number of value to read along each axis 
      INTEGER, DIMENSION(jpmax_dims) ::   idimsz      ! size of the dimensions of the variable
      INTEGER, DIMENSION(jpmax_dims) ::   ishape      ! size of the dimensions of the variable
      REAL(wp)                       ::   zscf, zofs  ! sacle_factor and add_offset
      INTEGER                        ::   itmp        ! temporary integer
      CHARACTER(LEN=256)             ::   clinfo      ! info character
      CHARACTER(LEN=256)             ::   clname      ! file name
      CHARACTER(LEN=1)               ::   clrankpv, cldmspc      ! 
      !---------------------------------------------------------------------
      !
      clname = iom_file(kiomid)%name   !   esier to read
      clinfo = '          iom_get_123d, file: '//trim(clname)//', var: '//trim(cdvar)
      ! local definition of the domain ?
      idom = kdom
      ! do we read the overlap 
      ! ugly patch SM+JMM+RB to overwrite global definition in some cases
      llnoov = (jpni * jpnj ) == jpnij .AND. .NOT. lk_agrif 
      ! check kcount and kstart optionals parameters...
      IF( PRESENT(kcount) .AND. (.NOT. PRESENT(kstart)) ) CALL ctl_stop(trim(clinfo), 'kcount present needs kstart present')
      IF( PRESENT(kstart) .AND. (.NOT. PRESENT(kcount)) ) CALL ctl_stop(trim(clinfo), 'kstart present needs kcount present')
      IF( PRESENT(kstart) .AND. idom /= jpdom_unknown   ) CALL ctl_stop(trim(clinfo), 'kstart present needs kdom = jpdom_unknown')

      ! Search for the variable in the data base (eventually actualize data)
      istop = nstop
      idvar = iom_varid( kiomid, cdvar )
      !
      IF( idvar > 0 ) THEN
         ! to write iom_file(kiomid)%dimsz in a shorter way !
         idimsz(:) = iom_file(kiomid)%dimsz(:, idvar) 
         inbdim = iom_file(kiomid)%ndims(idvar)            ! number of dimensions in the file
         idmspc = inbdim                                   ! number of spatial dimensions in the file
         IF( iom_file(kiomid)%luld(idvar) )   idmspc = inbdim - 1
         IF( idmspc > 3 )   CALL ctl_stop(trim(clinfo), 'the file has more than 3 spatial dimensions this case is not coded...') 
         !
         ! update idom definition...
         ! Identify the domain in case of jpdom_auto(glo/dta) definition
         IF( idom == jpdom_autoglo .OR. idom == jpdom_autodta ) THEN            
            IF( idom == jpdom_autoglo ) THEN   ;   idom = jpdom_global 
            ELSE                               ;   idom = jpdom_data
            ENDIF
            ind1 = INDEX( clname, '_', back = .TRUE. ) + 1
            ind2 = INDEX( clname, '.', back = .TRUE. ) - 1
            IF( ind2 > ind1 ) THEN   ;   IF( VERIFY( clname(ind1:ind2), '0123456789' ) == 0 )   idom = jpdom_local   ;   ENDIF
         ENDIF
         ! Identify the domain in case of jpdom_local definition
         IF( idom == jpdom_local ) THEN
            IF(     idimsz(1) == jpi               .AND. idimsz(2) == jpj               ) THEN   ;   idom = jpdom_local_full
            ELSEIF( idimsz(1) == nlci              .AND. idimsz(2) == nlcj              ) THEN   ;   idom = jpdom_local_noextra
            ELSEIF( idimsz(1) == (nlei - nldi + 1) .AND. idimsz(2) == (nlej - nldj + 1) ) THEN   ;   idom = jpdom_local_noovlap
            ELSE   ;   CALL ctl_stop( trim(clinfo), 'impossible to identify the local domain' )
            ENDIF
         ENDIF
         !
         ! check the consistency between input array and data rank in the file
         !
         ! initializations
         itime = 1
         IF( PRESENT(ktime) ) itime = ktime

         irankpv = 1 * COUNT( (/PRESENT(pv_r1d)/) ) + 2 * COUNT( (/PRESENT(pv_r2d)/) ) + 3 * COUNT( (/PRESENT(pv_r3d)/) )
         WRITE(clrankpv, fmt='(i1)') irankpv
         WRITE(cldmspc , fmt='(i1)') idmspc
         !
         IF(     idmspc <  irankpv ) THEN 
            CALL ctl_stop( TRIM(clinfo), 'The file has only '//cldmspc//' spatial dimension',   &
               &                         'it is impossible to read a '//clrankpv//'D array from this file...' )
         ELSEIF( idmspc == irankpv ) THEN
            IF( PRESENT(pv_r1d) .AND. idom /= jpdom_unknown )   &
               &   CALL ctl_stop( TRIM(clinfo), 'case not coded...You must use jpdom_unknown' )
         ELSEIF( idmspc >  irankpv ) THEN
               IF( PRESENT(pv_r2d) .AND. itime == 1 .AND. idimsz(3) == 1 .AND. idmspc == 3 ) THEN
                  CALL ctl_warn( trim(clinfo), '2D array but 3 spatial dimensions for the data...'              ,   &
                        &         'As the size of the z dimension is 1 and as we try to read the first record, ',   &
                        &         'we accept this case, even if there is a possible mix-up between z and time dimension' )   
                  idmspc = idmspc - 1
               ELSE
                  CALL ctl_stop( TRIM(clinfo), 'To keep iom lisibility, when reading a '//clrankpv//'D array,'         ,   &
                     &                         'we do not accept data with more than '//cldmspc//' spatial dimension',   &
                     &                         'Use ncwa -a to suppress the unnecessary dimensions' )
               ENDIF
         ENDIF

         !
         ! definition of istart and icnt
         !
         icnt  (:) = 1
         istart(:) = 1
         istart(idmspc+1) = itime

         IF(              PRESENT(kstart)       ) THEN ; istart(1:idmspc) = kstart(1:idmspc) ; icnt(1:idmspc) = kcount(1:idmspc)
         ELSE
            IF(           idom == jpdom_unknown ) THEN                                       ; icnt(1:idmspc) = idimsz(1:idmspc)
            ELSE 
               IF( .NOT. PRESENT(pv_r1d) ) THEN   !   not a 1D array
                  IF(     idom == jpdom_data    ) THEN ; istart(1:2) = (/ mig(1), mjg(1) /)  ! icnt(1:2) done bellow
                  ELSEIF( idom == jpdom_global  ) THEN ; istart(1:2) = (/ nimpp , njmpp  /)  ! icnt(1:2) done bellow
                  ENDIF
                  ! we do not read the overlap                     -> we start to read at nldi, nldj
! JMM + SM: ugly patch before getting the new version of lib_mpp)
!                  IF( idom /= jpdom_local_noovlap )   istart(1:2) = istart(1:2) + (/ nldi - 1, nldj - 1 /)
                  IF( llnoov .AND. idom /= jpdom_local_noovlap ) istart(1:2) = istart(1:2) + (/ nldi - 1, nldj - 1 /)
                  ! we do not read the overlap and the extra-halos -> from nldi to nlei and from nldj to nlej 
! JMM + SM: ugly patch before getting the new version of lib_mpp)
!                  icnt(1:2) = (/ nlei - nldi + 1, nlej - nldj + 1 /)
                  IF( llnoov ) THEN   ;   icnt(1:2) = (/ nlei - nldi + 1, nlej - nldj + 1 /)
                  ELSE                ;   icnt(1:2) = (/ nlci           , nlcj            /)
                  ENDIF
                  IF( PRESENT(pv_r3d) ) THEN
                     IF( idom == jpdom_data ) THEN   ; icnt(3) = jpkdta
                     ELSE                            ; icnt(3) = jpk
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         ! check that istart and icnt can be used with this file
         !-
         DO jl = 1, jpmax_dims
            itmp = istart(jl)+icnt(jl)-1
            IF( itmp > idimsz(jl) .AND. idimsz(jl) /= 0 ) THEN
               WRITE( ctmp1, FMT="('(istart(', i1, ') + icnt(', i1, ') - 1) = ', i5)" ) jl, jl, itmp
               WRITE( ctmp2, FMT="(' is larger than idimsz(', i1,') = ', i5)"         ) jl, idimsz(jl)
               CALL ctl_stop( trim(clinfo), 'start and count too big regarding to the size of the data, ', ctmp1, ctmp2 )     
            ENDIF
         END DO

         ! check that icnt matches the input array
         !-     
         IF( idom == jpdom_unknown ) THEN
            IF( irankpv == 1 )        ishape(1:1) = SHAPE(pv_r1d)
            IF( irankpv == 2 )        ishape(1:2) = SHAPE(pv_r2d)
            IF( irankpv == 3 )        ishape(1:3) = SHAPE(pv_r3d)
            ctmp1 = 'd'
         ELSE
            IF( irankpv == 2 ) THEN
! JMM + SM: ugly patch before getting the new version of lib_mpp)
!               ishape(1:2) = SHAPE(pv_r2d(nldi:nlei,nldj:nlej  ))   ;   ctmp1 = 'd(nldi:nlei,nldj:nlej)'
               IF( llnoov ) THEN ; ishape(1:2)=SHAPE(pv_r2d(nldi:nlei,nldj:nlej  )) ; ctmp1='d(nldi:nlei,nldj:nlej)'
               ELSE              ; ishape(1:2)=SHAPE(pv_r2d(1   :nlci,1   :nlcj  )) ; ctmp1='d(1:nlci,1:nlcj)'
               ENDIF
            ENDIF
            IF( irankpv == 3 ) THEN 
! JMM + SM: ugly patch before getting the new version of lib_mpp)
!               ishape(1:3) = SHAPE(pv_r3d(nldi:nlei,nldj:nlej,:))   ;   ctmp1 = 'd(nldi:nlei,nldj:nlej,:)'
               IF( llnoov ) THEN ; ishape(1:3)=SHAPE(pv_r3d(nldi:nlei,nldj:nlej,:)) ; ctmp1='d(nldi:nlei,nldj:nlej,:)'
               ELSE              ; ishape(1:3)=SHAPE(pv_r3d(1   :nlci,1   :nlcj,:)) ; ctmp1='d(1:nlci,1:nlcj,:)'
               ENDIF
            ENDIF
         ENDIF
         
         DO jl = 1, irankpv
            WRITE( ctmp2, FMT="(', ', i1,'): ', i5,' /= icnt(', i1,'):', i5)" ) jl, ishape(jl), jl, icnt(jl)
            IF( ishape(jl) /= icnt(jl) )   CALL ctl_stop( TRIM(clinfo), 'size(pv_r'//clrankpv//TRIM(ctmp1)//TRIM(ctmp2) )
         END DO

      ENDIF

      ! read the data
      !-     
      IF( idvar > 0 .AND. istop == nstop ) THEN   ! no additional errors until this point...
         !
         ! find the right index of the array to be read
! JMM + SM: ugly patch before getting the new version of lib_mpp)
!         IF( idom /= jpdom_unknown ) THEN   ;   ix1 = nldi   ;   ix2 = nlei      ;   iy1 = nldj   ;   iy2 = nlej
!         ELSE                               ;   ix1 = 1      ;   ix2 = icnt(1)   ;   iy1 = 1      ;   iy2 = icnt(2)
!         ENDIF
         IF( llnoov ) THEN
            IF( idom /= jpdom_unknown ) THEN   ;   ix1 = nldi   ;   ix2 = nlei      ;   iy1 = nldj   ;   iy2 = nlej
            ELSE                               ;   ix1 = 1      ;   ix2 = icnt(1)   ;   iy1 = 1      ;   iy2 = icnt(2)
            ENDIF
         ELSE
            IF( idom /= jpdom_unknown ) THEN   ;   ix1 = 1      ;   ix2 = nlci      ;   iy1 = 1      ;   iy2 = nlcj
            ELSE                               ;   ix1 = 1      ;   ix2 = icnt(1)   ;   iy1 = 1      ;   iy2 = icnt(2)
            ENDIF
         ENDIF
      
         SELECT CASE (iom_file(kiomid)%iolib)
         CASE (jpioipsl )   ;   CALL iom_ioipsl_get(  kiomid, idvar, inbdim, istart, icnt, ix1, ix2, iy1, iy2,   &
            &                                         pv_r1d, pv_r2d, pv_r3d )
         CASE (jpnf90   )   ;   CALL iom_nf90_get(    kiomid, idvar, inbdim, istart, icnt, ix1, ix2, iy1, iy2,   &
            &                                         pv_r1d, pv_r2d, pv_r3d )
         CASE (jprstdimg)   ;   CALL iom_rstdimg_get( kiomid, idom, idvar, ix1, ix2, iy1, iy2,   &
            &                                         pv_r1d, pv_r2d, pv_r3d )
         CASE DEFAULT    
            CALL ctl_stop( TRIM(clinfo)//' accepted IO library are only jpioipsl, jpnf90 and jprstdimg' )
         END SELECT

         IF( istop == nstop ) THEN   ! no additional errors until this point...
            IF(lwp) WRITE(numout,"(10x,' read ',a,' (rec: ',i4,') in ',a,' ok')") TRIM(cdvar), itime, TRIM(iom_file(kiomid)%name)
          
            !--- overlap areas and extra hallows (mpp)
            IF(     PRESENT(pv_r2d) .AND. idom /= jpdom_unknown ) THEN
               CALL lbc_lnk( pv_r2d,'Z',-999.,'no0' )
            ELSEIF( PRESENT(pv_r3d) .AND. idom /= jpdom_unknown ) THEN
               ! this if could be simplified with the new lbc_lnk that works with any size of the 3rd dimension
               IF( icnt(3) == jpk ) THEN
                  CALL lbc_lnk( pv_r3d,'Z',-999.,'no0' )
               ELSE   ! put some arbitrary value (a call to lbc_lnk will be done later...)
                  DO jj = nlcj+1, jpj   ;   pv_r3d(1:nlci, jj, :) = pv_r3d(1:nlci, nlej, :)   ;   END DO
                  DO ji = nlci+1, jpi   ;   pv_r3d(ji    , : , :) = pv_r3d(nlei  , :   , :)   ;   END DO
               ENDIF
            ENDIF
            
            ! C1D case : always call lbc_lnk to replicate the central value over the whole 3X3 domain
            IF( lk_c1d .AND. PRESENT(pv_r2d) )   CALL lbc_lnk( pv_r2d,'Z',1. )
            IF( lk_c1d .AND. PRESENT(pv_r3d) )   CALL lbc_lnk( pv_r3d,'Z',1. )
    
            !--- Apply scale_factor and offset
            zscf = iom_file(kiomid)%scf(idvar)      ! scale factor
            zofs = iom_file(kiomid)%ofs(idvar)      ! offset
            IF(     PRESENT(pv_r1d) ) THEN
               IF( zscf /= 1. )   pv_r1d(:) = pv_r1d(:) * zscf 
               IF( zofs /= 0. )   pv_r1d(:) = pv_r1d(:) + zofs
            ELSEIF( PRESENT(pv_r2d) ) THEN
!CDIR COLLAPSE
               IF( zscf /= 1.)   pv_r2d(:,:) = pv_r2d(:,:) * zscf
!CDIR COLLAPSE
               IF( zofs /= 0.)   pv_r2d(:,:) = pv_r2d(:,:) + zofs
            ELSEIF( PRESENT(pv_r3d) ) THEN
!CDIR COLLAPSE
               IF( zscf /= 1.)   pv_r3d(:,:,:) = pv_r3d(:,:,:) * zscf
!CDIR COLLAPSE
               IF( zofs /= 0.)   pv_r3d(:,:,:) = pv_r3d(:,:,:) + zofs
            ENDIF
            !
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE iom_get_123d


   SUBROUTINE iom_gettime( kiomid, ptime, cdvar, kntime, cdunits, cdcalendar )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE iom_gettime  ***
      !!
      !! ** Purpose : read the time axis cdvar in the file 
      !!--------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   kiomid     ! file Identifier
      REAL(wp), DIMENSION(:)     , INTENT(  out) ::   ptime      ! the time axis
      CHARACTER(len=*), OPTIONAL , INTENT(in   ) ::   cdvar      ! time axis name
      INTEGER         , OPTIONAL , INTENT(  out) ::   kntime     ! number of times in file
      CHARACTER(len=*), OPTIONAL , INTENT(  out) ::   cdunits    ! units attribute of time coordinate
      CHARACTER(len=*), OPTIONAL , INTENT(  out) ::   cdcalendar ! calendar attribute of 
      !
      INTEGER, DIMENSION(1) :: kdimsz
      INTEGER            ::   idvar    ! id of the variable
      CHARACTER(LEN=32)  ::   tname    ! local name of time coordinate
      CHARACTER(LEN=100) ::   clinfo   ! info character
      !---------------------------------------------------------------------
      !
      IF ( PRESENT(cdvar) ) THEN
         tname = cdvar
      ELSE
         tname = iom_file(kiomid)%uldname
      ENDIF
      IF( kiomid > 0 ) THEN
         clinfo = 'iom_gettime, file: '//trim(iom_file(kiomid)%name)//', var: '//trim(tname)
         IF ( PRESENT(kntime) ) THEN
            idvar  = iom_varid( kiomid, tname, kdimsz = kdimsz )
            kntime = kdimsz(1)
         ELSE
            idvar = iom_varid( kiomid, tname )
         ENDIF
         !
         ptime(:) = 0. ! default definition
         IF( idvar > 0 ) THEN
            IF( iom_file(kiomid)%ndims(idvar) == 1 ) THEN
               IF( iom_file(kiomid)%luld(idvar) ) THEN
                  IF( iom_file(kiomid)%dimsz(1,idvar) <= size(ptime) ) THEN
                     SELECT CASE (iom_file(kiomid)%iolib)
                     CASE (jpioipsl )   ;   CALL iom_ioipsl_gettime( kiomid, idvar, ptime, cdunits, cdcalendar )
                     CASE (jpnf90   )   ;   CALL iom_nf90_gettime(   kiomid, idvar, ptime, cdunits, cdcalendar )
                     CASE (jprstdimg)   ;   CALL ctl_stop( TRIM(clinfo)//' case IO library == jprstdimg not coded...' )
                     CASE DEFAULT    
                        CALL ctl_stop( TRIM(clinfo)//' accepted IO library are only jpioipsl, jpnf90 and jprstdimg' )
                     END SELECT
                  ELSE
                     WRITE(ctmp1,*) 'error with the size of ptime ',size(ptime),iom_file(kiomid)%dimsz(1,idvar)
                     CALL ctl_stop( trim(clinfo), trim(ctmp1) )
                  ENDIF
               ELSE
                  CALL ctl_stop( trim(clinfo), 'variable dimension is not unlimited... use iom_get' )
               ENDIF
            ELSE
               CALL ctl_stop( trim(clinfo), 'the variable has more than 1 dimension' )
            ENDIF
         ELSE
            CALL ctl_stop( trim(clinfo), 'variable not found in '//iom_file(kiomid)%name )
         ENDIF
      ENDIF
      !
   END SUBROUTINE iom_gettime


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_getatt
   !!----------------------------------------------------------------------
   SUBROUTINE iom_g0d_intatt( kiomid, cdatt, pvar )
      INTEGER         , INTENT(in   )                 ::   kiomid    ! Identifier of the file
      CHARACTER(len=*), INTENT(in   )                 ::   cdatt     ! Name of the attribute
      INTEGER         , INTENT(  out)                 ::   pvar      ! read field
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) THEN
            SELECT CASE (iom_file(kiomid)%iolib)
            CASE (jpioipsl )   ;   CALL ctl_stop('iom_getatt: only nf90 available')
            CASE (jpnf90   )   ;   CALL iom_nf90_getatt( kiomid, cdatt, pvar )
            CASE (jprstdimg)   ;   CALL ctl_stop('iom_getatt: only nf90 available')
            CASE DEFAULT    
               CALL ctl_stop( 'iom_g0d_att: accepted IO library are only jpioipsl, jpnf90 and jprstdimg' )
            END SELECT
         ENDIF
      ENDIF
   END SUBROUTINE iom_g0d_intatt


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_rstput
   !!----------------------------------------------------------------------
   SUBROUTINE iom_rp0d( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file 
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(wp)        , INTENT(in)                         ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      INTEGER :: ivid   ! variable id
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) THEN
            ivid = iom_varid( kiomid, cdvar, ldstop = .FALSE. )
            SELECT CASE (iom_file(kiomid)%iolib)
            CASE (jpioipsl )   ;   CALL iom_ioipsl_rstput( kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r0d = pvar )
            CASE (jpnf90   )   ;   CALL iom_nf90_rstput(   kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r0d = pvar )
            CASE (jprstdimg)   ;   IF( kt == kwrite )    CALL iom_rstdimg_rstput( kiomid, cdvar, ivid, pvar )
            CASE DEFAULT     
               CALL ctl_stop( 'iom_rp0d: accepted IO library are only jpioipsl, jpnf90 and jprstdimg' )
            END SELECT
         ENDIF
      ENDIF
   END SUBROUTINE iom_rp0d

   SUBROUTINE iom_rp1d( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file 
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(wp)        , INTENT(in), DIMENSION(          :) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      INTEGER :: ivid   ! variable id
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) THEN
            ivid = iom_varid( kiomid, cdvar, ldstop = .FALSE. )
            SELECT CASE (iom_file(kiomid)%iolib)
            CASE (jpioipsl )   ;   CALL iom_ioipsl_rstput( kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r1d = pvar )
            CASE (jpnf90   )   ;   CALL iom_nf90_rstput(   kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r1d = pvar )
            CASE (jprstdimg)   ;   IF( kt == kwrite )    CALL iom_rstdimg_rstput( kiomid, cdvar, ivid, pv_r1d = pvar )
            CASE DEFAULT     
               CALL ctl_stop( 'iom_rp1d: accepted IO library are only jpioipsl, jpnf90 and jprstdimg' )
            END SELECT
         ENDIF
      ENDIF
   END SUBROUTINE iom_rp1d

   SUBROUTINE iom_rp2d( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file 
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(wp)        , INTENT(in), DIMENSION(:,    :    ) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      INTEGER :: ivid   ! variable id
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) THEN
            ivid = iom_varid( kiomid, cdvar, ldstop = .FALSE. )
            SELECT CASE (iom_file(kiomid)%iolib)
            CASE (jpioipsl )   ;   CALL iom_ioipsl_rstput( kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r2d = pvar )
            CASE (jpnf90   )   ;   CALL iom_nf90_rstput(   kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r2d = pvar )
            CASE (jprstdimg)   ;   IF( kt == kwrite )   CALL iom_rstdimg_rstput( kiomid, cdvar, ivid, pv_r2d = pvar ) 
            CASE DEFAULT     
               CALL ctl_stop( 'iom_rp2d: accepted IO library are only jpioipsl, jpnf90 and jprstdimg' )
            END SELECT
         ENDIF
      ENDIF
   END SUBROUTINE iom_rp2d

   SUBROUTINE iom_rp3d( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file 
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(wp)        , INTENT(in),       DIMENSION(:,:,:) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      INTEGER :: ivid   ! variable id
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) THEN
            ivid = iom_varid( kiomid, cdvar, ldstop = .FALSE. )
            SELECT CASE (iom_file(kiomid)%iolib)
            CASE (jpioipsl )   ;   CALL iom_ioipsl_rstput( kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r3d = pvar )
            CASE (jpnf90   )   ;   CALL iom_nf90_rstput(   kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r3d = pvar )
            CASE (jprstdimg)   ;   IF( kt == kwrite )   CALL iom_rstdimg_rstput( kiomid, cdvar, ivid, pv_r3d = pvar )
            CASE DEFAULT     
               CALL ctl_stop( 'iom_rp3d: accepted IO library are only jpioipsl and jprstdimg' )
            END SELECT
         ENDIF
      ENDIF
   END SUBROUTINE iom_rp3d


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_put
   !!----------------------------------------------------------------------
   SUBROUTINE iom_p0d( cdname, pfield0d )
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      REAL(wp)        , INTENT(in) ::   pfield0d
#if defined key_iomput
      CALL xios_send_field(cdname, (/pfield0d/))
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield0d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p0d

   SUBROUTINE iom_p1d( cdname, pfield1d )
      CHARACTER(LEN=*)          , INTENT(in) ::   cdname
      REAL(wp),     DIMENSION(:), INTENT(in) ::   pfield1d
#if defined key_iomput
      CALL xios_send_field( cdname, RESHAPE( (/pfield1d/), (/1,1,SIZE(pfield1d)/) ) )
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield1d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p1d

   SUBROUTINE iom_p2d( cdname, pfield2d )
      CHARACTER(LEN=*)            , INTENT(in) ::   cdname
      REAL(wp),     DIMENSION(:,:), INTENT(in) ::   pfield2d
#if defined key_iomput
      CALL xios_send_field(cdname, pfield2d)
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield2d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p2d

   SUBROUTINE iom_p3d( cdname, pfield3d )
      CHARACTER(LEN=*)                , INTENT(in) ::   cdname
      REAL(wp),       DIMENSION(:,:,:), INTENT(in) ::   pfield3d
#if defined key_iomput
      CALL xios_send_field(cdname, pfield3d)
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield3d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p3d
   !!----------------------------------------------------------------------

#if defined key_iomput

   SUBROUTINE iom_set_domain_attr( cdid, ni_glo, nj_glo, ibegin, jbegin, ni, nj, zoom_ibegin, zoom_jbegin, zoom_ni, zoom_nj,   &
      &                                    data_dim, data_ibegin, data_ni, data_jbegin, data_nj, lonvalue, latvalue, mask )
      CHARACTER(LEN=*)                 , INTENT(in) ::   cdid
      INTEGER                , OPTIONAL, INTENT(in) ::   ni_glo, nj_glo, ibegin, jbegin, ni, nj
      INTEGER                , OPTIONAL, INTENT(in) ::   data_dim, data_ibegin, data_ni, data_jbegin, data_nj
      INTEGER                , OPTIONAL, INTENT(in) ::   zoom_ibegin, zoom_jbegin, zoom_ni, zoom_nj
      REAL(wp), DIMENSION(:) , OPTIONAL, INTENT(in) ::   lonvalue, latvalue
      LOGICAL, DIMENSION(:,:), OPTIONAL, INTENT(in) ::   mask

      IF ( xios_is_valid_domain     (cdid) ) THEN
         CALL xios_set_domain_attr     ( cdid, ni_glo=ni_glo, nj_glo=nj_glo, ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,   &
            &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni, data_jbegin=data_jbegin, data_nj=data_nj ,   &
            &    zoom_ibegin=zoom_ibegin, zoom_jbegin=zoom_jbegin, zoom_ni=zoom_ni, zoom_nj=zoom_nj,                       &
            &    lonvalue=lonvalue, latvalue=latvalue,mask=mask )
      ENDIF

      IF ( xios_is_valid_domaingroup(cdid) ) THEN
         CALL xios_set_domaingroup_attr( cdid, ni_glo=ni_glo, nj_glo=nj_glo, ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,   &
            &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni, data_jbegin=data_jbegin, data_nj=data_nj ,   &
            &    zoom_ibegin=zoom_ibegin, zoom_jbegin=zoom_jbegin, zoom_ni=zoom_ni, zoom_nj=zoom_nj,                       &
            &    lonvalue=lonvalue, latvalue=latvalue,mask=mask )
      ENDIF
      CALL xios_solve_inheritance()

   END SUBROUTINE iom_set_domain_attr


   SUBROUTINE iom_set_axis_attr( cdid, paxis )
      CHARACTER(LEN=*)      , INTENT(in) ::   cdid
      REAL(wp), DIMENSION(:), INTENT(in) ::   paxis
      IF ( xios_is_valid_axis     (cdid) )   CALL xios_set_axis_attr     ( cdid, size=size(paxis),value=paxis )
      IF ( xios_is_valid_axisgroup(cdid) )   CALL xios_set_axisgroup_attr( cdid, size=size(paxis),value=paxis )
      CALL xios_solve_inheritance()
   END SUBROUTINE iom_set_axis_attr


   SUBROUTINE iom_set_field_attr( cdid, freq_op, freq_offset )
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   freq_op
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   freq_offset
      IF ( xios_is_valid_field     (cdid) )   CALL xios_set_field_attr     ( cdid, freq_op=freq_op, freq_offset=freq_offset )
      IF ( xios_is_valid_fieldgroup(cdid) )   CALL xios_set_fieldgroup_attr( cdid, freq_op=freq_op, freq_offset=freq_offset )
      CALL xios_solve_inheritance()
   END SUBROUTINE iom_set_field_attr


   SUBROUTINE iom_set_file_attr( cdid, name, name_suffix )
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   name, name_suffix
      IF ( xios_is_valid_file     (cdid) )   CALL xios_set_file_attr     ( cdid, name=name, name_suffix=name_suffix )
      IF ( xios_is_valid_filegroup(cdid) )   CALL xios_set_filegroup_attr( cdid, name=name, name_suffix=name_suffix )
      CALL xios_solve_inheritance()
   END SUBROUTINE iom_set_file_attr


   SUBROUTINE iom_get_file_attr( cdid, name, name_suffix, output_freq )
      CHARACTER(LEN=*)          , INTENT(in ) ::   cdid
      CHARACTER(LEN=*),OPTIONAL , INTENT(out) ::   name, name_suffix, output_freq
      LOGICAL                                 ::   llexist1,llexist2,llexist3
      !---------------------------------------------------------------------
      IF( PRESENT( name        ) )   name = ''          ! default values
      IF( PRESENT( name_suffix ) )   name_suffix = ''
      IF( PRESENT( output_freq ) )   output_freq = ''
      IF ( xios_is_valid_file     (cdid) ) THEN
         CALL xios_solve_inheritance()
         CALL xios_is_defined_file_attr     ( cdid, name = llexist1, name_suffix = llexist2, output_freq = llexist3)
         IF(llexist1)   CALL xios_get_file_attr     ( cdid, name = name )
         IF(llexist2)   CALL xios_get_file_attr     ( cdid, name_suffix = name_suffix )
         IF(llexist3)   CALL xios_get_file_attr     ( cdid, output_freq = output_freq )
      ENDIF
      IF ( xios_is_valid_filegroup(cdid) ) THEN
         CALL xios_solve_inheritance()
         CALL xios_is_defined_filegroup_attr( cdid, name = llexist1, name_suffix = llexist2, output_freq = llexist3)
         IF(llexist1)   CALL xios_get_filegroup_attr( cdid, name = name )
         IF(llexist2)   CALL xios_get_filegroup_attr( cdid, name_suffix = name_suffix )
         IF(llexist3)   CALL xios_get_filegroup_attr( cdid, output_freq = output_freq )
      ENDIF
   END SUBROUTINE iom_get_file_attr


   SUBROUTINE iom_set_grid_attr( cdid, mask )
      CHARACTER(LEN=*)                   , INTENT(in) ::   cdid
      LOGICAL, DIMENSION(:,:,:), OPTIONAL, INTENT(in) ::   mask
      IF ( xios_is_valid_grid     (cdid) )   CALL xios_set_grid_attr     ( cdid, mask=mask )
      IF ( xios_is_valid_gridgroup(cdid) )   CALL xios_set_gridgroup_attr( cdid, mask=mask )
      CALL xios_solve_inheritance()
   END SUBROUTINE iom_set_grid_attr


   SUBROUTINE set_grid( cdgrd, plon, plat )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_grid  ***
      !!
      !! ** Purpose :   define horizontal grids
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=1)            , INTENT(in) ::   cdgrd
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   plon
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   plat
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zmask
      INTEGER  :: ni,nj
      
      ni=nlei-nldi+1 ; nj=nlej-nldj+1

      CALL iom_set_domain_attr("grid_"//cdgrd, ni_glo=jpiglo, nj_glo=jpjglo, ibegin=nimpp+nldi-1, jbegin=njmpp+nldj-1, ni=ni, nj=nj)
      CALL iom_set_domain_attr("grid_"//cdgrd, data_dim=2, data_ibegin = 1-nldi, data_ni = jpi, data_jbegin = 1-nldj, data_nj = jpj)
      CALL iom_set_domain_attr("grid_"//cdgrd, lonvalue = RESHAPE(plon(nldi:nlei, nldj:nlej),(/ ni*nj /)),   &
         &                                     latvalue = RESHAPE(plat(nldi:nlei, nldj:nlej),(/ ni*nj /)))  

      IF ( ln_mskland ) THEN
         ! mask land points, keep values on coast line -> specific mask for U, V and W points
         SELECT CASE ( cdgrd )
         CASE('T')   ;   zmask(:,:,:)       = tmask(:,:,:)
         CASE('U')   ;   zmask(2:jpim1,:,:) = tmask(2:jpim1,:,:) + tmask(3:jpi,:,:)   ;   CALL lbc_lnk( zmask, 'U', 1. )
         CASE('V')   ;   zmask(:,2:jpjm1,:) = tmask(:,2:jpjm1,:) + tmask(:,3:jpi,:)   ;   CALL lbc_lnk( zmask, 'V', 1. )
         CASE('W')   ;   zmask(:,:,2:jpk  ) = tmask(:,:,1:jpkm1) + tmask(:,:,2:jpk)   ;   zmask(:,:,1) = tmask(:,:,1)
         END SELECT
         !
         CALL iom_set_domain_attr( "grid_"//cdgrd       , mask = RESHAPE(zmask(nldi:nlei,nldj:nlej,1),(/ni,nj    /)) /= 0. )
         CALL iom_set_grid_attr  ( "grid_"//cdgrd//"_3D", mask = RESHAPE(zmask(nldi:nlei,nldj:nlej,:),(/ni,nj,jpk/)) /= 0. )
      ENDIF
      
   END SUBROUTINE set_grid


   SUBROUTINE set_scalar
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_scalar  ***
      !!
      !! ** Purpose :   define fake grids for scalar point
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(1,1) ::   zz = 1.
      !!----------------------------------------------------------------------
      CALL iom_set_domain_attr('scalarpoint', ni_glo=jpnij, nj_glo=1, ibegin=narea, jbegin=1, ni=1, nj=1)
      CALL iom_set_domain_attr('scalarpoint', data_dim=1)
      CALL iom_set_domain_attr('scalarpoint', lonvalue=(/ zz /), latvalue=(/ zz /))

   END SUBROUTINE set_scalar


   SUBROUTINE set_xmlatt
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_xmlatt  ***
      !!
      !! ** Purpose :   automatic definitions of some of the xml attributs...
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=1),DIMENSION( 3) ::   clgrd                    ! suffix name
      CHARACTER(len=256)             ::   clsuff                   ! suffix name
      CHARACTER(len=1)               ::   cl1                      ! 1 character
      CHARACTER(len=2)               ::   cl2                      ! 1 character
      INTEGER                        ::   ji, jg                   ! loop counters
      INTEGER                        ::   ix, iy                   ! i-,j- index
      REAL(wp)        ,DIMENSION(11) ::   zlontao                  ! longitudes of tao    moorings
      REAL(wp)        ,DIMENSION( 7) ::   zlattao                  ! latitudes  of tao    moorings
      REAL(wp)        ,DIMENSION( 4) ::   zlonrama                 ! longitudes of rama   moorings
      REAL(wp)        ,DIMENSION(11) ::   zlatrama                 ! latitudes  of rama   moorings
      REAL(wp)        ,DIMENSION( 3) ::   zlonpira                 ! longitudes of pirata moorings
      REAL(wp)        ,DIMENSION( 9) ::   zlatpira                 ! latitudes  of pirata moorings
      !!----------------------------------------------------------------------
      ! 
      ! frequency of the call of iom_put (attribut: freq_op)
      WRITE(cl1,'(i1)')        1   ;   CALL iom_set_field_attr('field_definition', freq_op = cl1//'ts', freq_offset='0ts')
      WRITE(cl1,'(i1)')  nn_fsbc   ;   CALL iom_set_field_attr('SBC'             , freq_op = cl1//'ts', freq_offset='0ts')
      WRITE(cl1,'(i1)') nn_dttrc   ;   CALL iom_set_field_attr('ptrc_T'          , freq_op = cl1//'ts', freq_offset='0ts')
      WRITE(cl1,'(i1)') nn_dttrc   ;   CALL iom_set_field_attr('diad_T'          , freq_op = cl1//'ts', freq_offset='0ts')
       
      ! output file names (attribut: name)
      DO ji = 1, 9
         WRITE(cl1,'(i1)') ji 
         CALL iom_update_file_name('file'//cl1)
      END DO
      DO ji = 1, 99
         WRITE(cl2,'(i2.2)') ji 
         CALL iom_update_file_name('file'//cl2)
      END DO

      ! Zooms...
      clgrd = (/ 'T', 'U', 'W' /) 
      DO jg = 1, SIZE(clgrd)                                                                   ! grid type
         cl1 = clgrd(jg)
         ! Equatorial section (attributs: jbegin, ni, name_suffix)
         CALL dom_ngb( 0., 0., ix, iy, cl1 )
         CALL iom_set_domain_attr ('Eq'//cl1, zoom_jbegin=iy, zoom_ni=jpiglo)
         CALL iom_get_file_attr   ('Eq'//cl1, name_suffix = clsuff             )
         CALL iom_set_file_attr   ('Eq'//cl1, name_suffix = TRIM(clsuff)//'_Eq')
         CALL iom_update_file_name('Eq'//cl1)
      END DO
      ! TAO moorings (attributs: ibegin, jbegin, name_suffix)
      zlontao = (/ 137.0, 147.0, 156.0, 165.0, -180.0, -170.0, -155.0, -140.0, -125.0, -110.0, -95.0 /)
      zlattao = (/  -8.0,  -5.0,  -2.0,   0.0,    2.0,    5.0,    8.0 /)
      CALL set_mooring( zlontao, zlattao )
      ! RAMA moorings (attributs: ibegin, jbegin, name_suffix)
      zlonrama = (/  55.0,  67.0, 80.5, 90.0 /)
      zlatrama = (/ -16.0, -12.0, -8.0, -4.0, -1.5, 0.0, 1.5, 4.0, 8.0, 12.0, 15.0 /)
      CALL set_mooring( zlonrama, zlatrama )
      ! PIRATA moorings (attributs: ibegin, jbegin, name_suffix)
      zlonpira = (/ -38.0, -23.0, -10.0 /)
      zlatpira = (/ -19.0, -14.0,  -8.0, 0.0, 4.0, 8.0, 12.0, 15.0, 20.0 /)
      CALL set_mooring( zlonpira, zlatpira )
      
   END SUBROUTINE set_xmlatt


   SUBROUTINE set_mooring( plon, plat)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_mooring  ***
      !!
      !! ** Purpose :   automatic definitions of moorings xml attributs...
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(in) ::  plon, plat           ! longitudes/latitudes oft the mooring
      !
!!$      CHARACTER(len=1),DIMENSION(4) ::   clgrd = (/ 'T', 'U', 'V', 'W' /)   ! suffix name
      CHARACTER(len=1),DIMENSION(1) ::   clgrd = (/ 'T' /)        ! suffix name
      CHARACTER(len=256)            ::   clname                   ! file name
      CHARACTER(len=256)            ::   clsuff                   ! suffix name
      CHARACTER(len=1)              ::   cl1                      ! 1 character
      CHARACTER(len=6)              ::   clon,clat                ! name of longitude, latitude
      INTEGER                       ::   ji, jj, jg               ! loop counters
      INTEGER                       ::   ix, iy                   ! i-,j- index
      REAL(wp)                      ::   zlon, zlat
      !!----------------------------------------------------------------------
      DO jg = 1, SIZE(clgrd)
         cl1 = clgrd(jg)
         DO ji = 1, SIZE(plon)
            DO jj = 1, SIZE(plat)
               zlon = plon(ji)
               zlat = plat(jj)
               ! modifications for RAMA moorings
               IF( zlon ==  67. .AND. zlat ==  15. )   zlon =  65.
               IF( zlon ==  90. .AND. zlat <=  -4. )   zlon =  95.
               IF( zlon ==  95. .AND. zlat ==  -4. )   zlat =  -5.
               ! modifications for PIRATA moorings
               IF( zlon == -38. .AND. zlat == -19. )   zlon = -34.
               IF( zlon == -38. .AND. zlat == -14. )   zlon = -32.
               IF( zlon == -38. .AND. zlat ==  -8. )   zlon = -30.
               IF( zlon == -38. .AND. zlat ==   0. )   zlon = -35.
               IF( zlon == -23. .AND. zlat ==  20. )   zlat =  21.
               IF( zlon == -10. .AND. zlat == -14. )   zlat = -10.
               IF( zlon == -10. .AND. zlat ==  -8. )   zlat =  -6.
               IF( zlon == -10. .AND. zlat ==   4. ) THEN   ;   zlon = 0.   ;   zlat = 0.   ;   ENDIF
               CALL dom_ngb( zlon, zlat, ix, iy, cl1 )
               IF( zlon >= 0. ) THEN  
                  IF( zlon == REAL(NINT(zlon), wp) ) THEN   ;   WRITE(clon, '(i3,  a)') NINT( zlon), 'e'
                  ELSE                                      ;   WRITE(clon, '(f5.1,a)')       zlon , 'e'
                  ENDIF
               ELSE             
                  IF( zlon == REAL(NINT(zlon), wp) ) THEN   ;   WRITE(clon, '(i3,  a)') NINT(-zlon), 'w'
                  ELSE                                      ;   WRITE(clon, '(f5.1,a)')      -zlon , 'w'
                  ENDIF
               ENDIF
               IF( zlat >= 0. ) THEN  
                  IF( zlat == REAL(NINT(zlat), wp) ) THEN   ;   WRITE(clat, '(i2,  a)') NINT( zlat), 'n'
                  ELSE                                      ;   WRITE(clat, '(f4.1,a)')       zlat , 'n'
                  ENDIF
               ELSE             
                  IF( zlat == REAL(NINT(zlat), wp) ) THEN   ;   WRITE(clat, '(i2,  a)') NINT(-zlat), 's'
                  ELSE                                      ;   WRITE(clat, '(f4.1,a)')      -zlat , 's'
                  ENDIF
               ENDIF
               clname = TRIM(ADJUSTL(clat))//TRIM(ADJUSTL(clon))
               CALL iom_set_domain_attr (TRIM(clname)//cl1, zoom_ibegin= ix, zoom_jbegin= iy)
               CALL iom_get_file_attr   (TRIM(clname)//cl1, name_suffix = clsuff                         )
               CALL iom_set_file_attr   (TRIM(clname)//cl1, name_suffix = TRIM(clsuff)//'_'//TRIM(clname))
               CALL iom_update_file_name(TRIM(clname)//cl1)
            END DO
         END DO
      END DO
      
   END SUBROUTINE set_mooring

   
   SUBROUTINE iom_update_file_name( cdid )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_update_file_name  ***
      !!
      !! ** Purpose :   
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      !
      CHARACTER(LEN=256) ::   clname
      CHARACTER(LEN=20)  ::   clfreq
      CHARACTER(LEN=20)  ::   cldate
      INTEGER            ::   idx
      INTEGER            ::   jn
      INTEGER            ::   itrlen
      INTEGER            ::   iyear, imonth, iday, isec
      REAL(wp)           ::   zsec
      LOGICAL            ::   llexist
      !!----------------------------------------------------------------------

      DO jn = 1,2

         IF( jn == 1 )   CALL iom_get_file_attr( cdid, name        = clname, output_freq = clfreq )
         IF( jn == 2 )   CALL iom_get_file_attr( cdid, name_suffix = clname )

         IF ( TRIM(clname) /= '' ) THEN 

            idx = INDEX(clname,'@expname@') + INDEX(clname,'@EXPNAME@')
            DO WHILE ( idx /= 0 ) 
               clname = clname(1:idx-1)//TRIM(cexper)//clname(idx+9:LEN_TRIM(clname))
               idx = INDEX(clname,'@expname@') + INDEX(clname,'@EXPNAME@')
            END DO

            idx = INDEX(clname,'@freq@') + INDEX(clname,'@FREQ@')
            DO WHILE ( idx /= 0 ) 
               IF ( TRIM(clfreq) /= '' ) THEN
                  itrlen = LEN_TRIM(clfreq)
                  IF ( clfreq(itrlen-1:itrlen) == 'mo' ) clfreq = clfreq(1:itrlen-1)
                  clname = clname(1:idx-1)//TRIM(clfreq)//clname(idx+6:LEN_TRIM(clname))
               ELSE
                  CALL ctl_stop('error in the name of file id '//TRIM(cdid),   &
                     & ' attribute output_freq is undefined -> cannot replace @freq@ in '//TRIM(clname) )
               ENDIF
               idx = INDEX(clname,'@freq@') + INDEX(clname,'@FREQ@')
            END DO

            idx = INDEX(clname,'@startdate@') + INDEX(clname,'@STARTDATE@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( fjulday - rdttra(1) / rday )
               clname = clname(1:idx-1)//TRIM(cldate)//clname(idx+11:LEN_TRIM(clname))
               idx = INDEX(clname,'@startdate@') + INDEX(clname,'@STARTDATE@')
            END DO

            idx = INDEX(clname,'@startdatefull@') + INDEX(clname,'@STARTDATEFULL@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( fjulday - rdttra(1) / rday, ldfull = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)//clname(idx+15:LEN_TRIM(clname))
               idx = INDEX(clname,'@startdatefull@') + INDEX(clname,'@STARTDATEFULL@')
            END DO

            idx = INDEX(clname,'@enddate@') + INDEX(clname,'@ENDDATE@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( fjulday + rdttra(1) / rday * REAL( nitend - nit000, wp ), ld24 = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)//clname(idx+9:LEN_TRIM(clname))
               idx = INDEX(clname,'@enddate@') + INDEX(clname,'@ENDDATE@')
            END DO

            idx = INDEX(clname,'@enddatefull@') + INDEX(clname,'@ENDDATEFULL@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( fjulday + rdttra(1) / rday * REAL( nitend - nit000, wp ), ld24 = .TRUE., ldfull = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)//clname(idx+13:LEN_TRIM(clname))
               idx = INDEX(clname,'@enddatefull@') + INDEX(clname,'@ENDDATEFULL@')
            END DO

            IF( jn == 1 )   CALL iom_set_file_attr( cdid, name        = clname )
            IF( jn == 2 )   CALL iom_set_file_attr( cdid, name_suffix = clname )

         ENDIF

      END DO

   END SUBROUTINE iom_update_file_name


   FUNCTION iom_sdate( pjday, ld24, ldfull )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_sdate  ***
      !!
      !! ** Purpose :   send back the date corresponding to the given julian day
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )           ::   pjday         ! julian day
      LOGICAL , INTENT(in   ), OPTIONAL ::   ld24          ! true to force 24:00 instead of 00:00
      LOGICAL , INTENT(in   ), OPTIONAL ::   ldfull        ! true to get the compleate date: yyyymmdd_hh:mm:ss
      !
      CHARACTER(LEN=20) ::   iom_sdate
      CHARACTER(LEN=50) ::   clfmt                         !  format used to write the date 
      INTEGER           ::   iyear, imonth, iday, ihour, iminute, isec
      REAL(wp)          ::   zsec
      LOGICAL           ::   ll24, llfull
      !
      IF( PRESENT(ld24) ) THEN   ;   ll24 = ld24
      ELSE                       ;   ll24 = .FALSE.
      ENDIF

      IF( PRESENT(ldfull) ) THEN   ;   llfull = ldfull
      ELSE                         ;   llfull = .FALSE.
      ENDIF

      CALL ju2ymds( pjday, iyear, imonth, iday, zsec )
      isec = NINT(zsec)

      IF ( ll24 .AND. isec == 0 ) THEN   ! 00:00 of the next day -> move to 24:00 of the current day
         CALL ju2ymds( pjday - 1., iyear, imonth, iday, zsec )
         isec = 86400
      ENDIF

      IF( iyear < 10000 ) THEN   ;   clfmt = "i4.4,2i2.2"                ! format used to write the date 
      ELSE                       ;   WRITE(clfmt, "('i',i1,',2i2.2')") INT(LOG10(REAL(iyear,wp))) + 1
      ENDIF
      
      IF( llfull ) THEN 
         clfmt = TRIM(clfmt)//",'_',i2.2,':',i2.2,':',i2.2"
         ihour   = isec / 3600
         isec    = MOD(isec, 3600)
         iminute = isec / 60
         isec    = MOD(isec, 60)
         WRITE(iom_sdate, '('//TRIM(clfmt)//')') iyear, imonth, iday, ihour, iminute, isec    ! date of the end of run
      ELSE
         WRITE(iom_sdate, '('//TRIM(clfmt)//')') iyear, imonth, iday                          ! date of the end of run
      ENDIF

   END FUNCTION iom_sdate

#else

   SUBROUTINE iom_setkt( kt )
      INTEGER, INTENT(in   )::   kt 
      IF( .FALSE. )   WRITE(numout,*) kt   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_setkt

#endif
   
   !!======================================================================
END MODULE iom
