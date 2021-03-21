MODULE obcdta
   !!==============================================================================
   !!                            ***  MODULE obcdta  ***
   !! Open boundary data : read the data for the open boundaries.
   !!==============================================================================
   !! History :  OPA  ! 1998-05 (J.M. Molines) Original code
   !!            8.5  ! 2002-10 (C. Talandier, A-M. Treguier) Free surface, F90
   !!   NEMO     1.0  ! 2004-06 (F. Durand, A-M. Treguier) Netcdf BC files on input
   !!            3.0  ! 2007-2008 (C. Langlais, P. Mathiot, J.M. Molines) high frequency boundaries data
   !!------------------------------------------------------------------------------
#if defined key_obc
   !!------------------------------------------------------------------------------
   !!   'key_obc'         :                                Open Boundary Conditions
   !!------------------------------------------------------------------------------
   !!   obc_dta           : read u, v, t, s data along each open boundary
   !!------------------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE phycst          ! physical constants
   USE obc_par         ! ocean open boundary conditions
   USE obc_oce         ! ocean open boundary conditions
   USE in_out_manager  ! I/O logical units
   USE lib_mpp         ! distributed memory computing
   USE dynspg_oce      ! ocean: surface pressure gradient
   USE ioipsl          ! now only for  ymds2ju function 
   USE iom             ! 

   IMPLICIT NONE
   PRIVATE

   PUBLIC   obc_dta         ! routine  called by step.F90
   PUBLIC   obc_dta_bt      ! routine  called by dynspg_ts.F90
   PUBLIC   obc_dta_alloc   ! function called by obcini.F90

   !$AGRIF_DO_NOT_TREAT
   REAL(wp),  DIMENSION(2)              ::   zjcnes_obc   ! 
   REAL(wp),  DIMENSION(:), ALLOCATABLE ::   ztcobc
   !$AGRIF_END_DO_NOT_TREAT
   REAL(wp) :: rdt_obc
   REAL(wp) :: zjcnes
   INTEGER :: imm0, iyy0, idd0, iyy, imm, idd
   INTEGER :: nt_a=2, nt_b=1, itobc, ndate0_cnes, nday_year0
   INTEGER ::  itobce, itobcw, itobcs, itobcn, itobc_b  ! number of time steps in OBC files

   INTEGER ::   ntobc        ! where we are in the obc file
   INTEGER ::   ntobc_b      ! first record used
   INTEGER ::   ntobc_a      ! second record used

   CHARACTER (len=40) ::   cl_obc_eTS, cl_obc_eU   ! name of data files
   CHARACTER (len=40) ::   cl_obc_wTS, cl_obc_wU   !   -       -
   CHARACTER (len=40) ::   cl_obc_nTS, cl_obc_nV   !   -       -
   CHARACTER (len=40) ::   cl_obc_sTS, cl_obc_sV   !   -       -

   ! bt arrays for interpolating time dependent data on the boundaries
   INTEGER ::   nt_m=0, ntobc_m
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ubtedta, vbtedta, sshedta      ! East
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ubtwdta, vbtwdta, sshwdta	   ! West
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ubtndta, vbtndta, sshndta	   ! North
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ubtsdta, vbtsdta, sshsdta	   ! South
   ! arrays used for interpolating time dependent data on the boundaries
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: uedta, vedta, tedta, sedta    ! East
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: uwdta, vwdta, twdta, swdta    ! West
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: undta, vndta, tndta, sndta    ! North
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: usdta, vsdta, tsdta, ssdta    ! South

   ! Masks set to .TRUE. after successful allocation below
   LOGICAL , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ltemsk, luemsk, lvemsk  ! boolean msks
   LOGICAL , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ltwmsk, luwmsk, lvwmsk  ! used for outliers
   LOGICAL , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ltnmsk, lunmsk, lvnmsk  ! checks
   LOGICAL , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ltsmsk, lusmsk, lvsmsk

   !! * Substitutions
#  include "obc_vectopt_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obcdta.F90 3565 2012-11-15 18:05:12Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION obc_dta_alloc()
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE obc_dta_alloc  ***
      !!-------------------------------------------------------------------
      INTEGER :: ierr(2)
      !!-------------------------------------------------------------------
# if defined key_dynspg_ts
      ALLOCATE(   &     ! time-splitting : 0:jptobc
         ! bt arrays for interpolating time dependent data on the boundaries
         &      ubtedta  (jpj,0:jptobc) , vbtedta  (jpj,0:jptobc) , sshedta  (jpj,0:jptobc) ,    &
         &      ubtwdta  (jpj,0:jptobc) , vbtwdta  (jpj,0:jptobc) , sshwdta  (jpj,0:jptobc) ,    &
         &      ubtndta  (jpi,0:jptobc) , vbtndta  (jpi,0:jptobc) , sshndta  (jpi,0:jptobc) ,    &
         &      ubtsdta  (jpi,0:jptobc) , vbtsdta  (jpi,0:jptobc) , sshsdta  (jpi,0:jptobc) ,    &
         ! arrays used for interpolating time dependent data on the boundaries
         &      uedta(jpj,jpk,0:jptobc) , vedta(jpj,jpk,0:jptobc)                           ,     &
         &      tedta(jpj,jpk,0:jptobc) , sedta(jpj,jpk,0:jptobc)                           ,     &
         &      uwdta(jpj,jpk,0:jptobc) , vwdta(jpj,jpk,0:jptobc)                           ,     &
         &      twdta(jpj,jpk,0:jptobc) , swdta(jpj,jpk,0:jptobc)                           ,     &
         &      undta(jpi,jpk,0:jptobc) , vndta(jpi,jpk,0:jptobc)                           ,     &
         &      tndta(jpi,jpk,0:jptobc) , sndta(jpi,jpk,0:jptobc)                           ,     &
         &      usdta(jpi,jpk,0:jptobc) , vsdta(jpi,jpk,0:jptobc)                           ,     &
         &      tsdta(jpi,jpk,0:jptobc) , ssdta(jpi,jpk,0:jptobc)                           , STAT=ierr(1) )
# else
      ALLOCATE(   &     ! no time splitting : 1:jptobc
         ! bt arrays for interpolating time dependent data on the boundaries
         &      ubtedta  (jpj,jptobc) , vbtedta  (jpj,jptobc) , sshedta  (jpj,jptobc)  ,     &
         &      ubtwdta  (jpj,jptobc) , vbtwdta  (jpj,jptobc) , sshwdta  (jpj,jptobc)  ,     &
         &      ubtndta  (jpi,jptobc) , vbtndta  (jpi,jptobc) , sshndta  (jpi,jptobc)  ,     &
         &      ubtsdta  (jpi,jptobc) , vbtsdta  (jpi,jptobc) , sshsdta  (jpi,jptobc)  ,     &
         ! arrays used for interpolating time dependent data on the boundaries
         &      uedta(jpj,jpk,jptobc) , vedta(jpj,jpk,jptobc)                          ,     &
         &      tedta(jpj,jpk,jptobc) , sedta(jpj,jpk,jptobc)                          ,     &
         &      uwdta(jpj,jpk,jptobc) , vwdta(jpj,jpk,jptobc)                          ,     &
         &      twdta(jpj,jpk,jptobc) , swdta(jpj,jpk,jptobc)                          ,     &
         &      undta(jpi,jpk,jptobc) , vndta(jpi,jpk,jptobc)                          ,     &
         &      tndta(jpi,jpk,jptobc) , sndta(jpi,jpk,jptobc)                          ,     &
         &      usdta(jpi,jpk,jptobc) , vsdta(jpi,jpk,jptobc)                          ,     &
         &      tsdta(jpi,jpk,jptobc) , ssdta(jpi,jpk,jptobc)                          , STAT=ierr(1) )
# endif

      ALLOCATE( ltemsk(jpj,jpk) , luemsk(jpj,jpk) , lvemsk(jpj,jpk) ,     &
         &      ltwmsk(jpj,jpk) , luwmsk(jpj,jpk) , lvwmsk(jpj,jpk) ,     &
         &      ltnmsk(jpi,jpk) , lunmsk(jpi,jpk) , lvnmsk(jpi,jpk) ,     &
         &      ltsmsk(jpi,jpk) , lusmsk(jpi,jpk) , lvsmsk(jpi,jpk) , STAT=ierr(2) )

      obc_dta_alloc = MAXVAL( ierr )
      IF( lk_mpp )   CALL mpp_sum( obc_dta_alloc )

      IF( obc_dta_alloc == 0 )  THEN         ! Initialise mask values following successful allocation
         !      east            !          west            !          north           !          south           !
         ltemsk(:,:) = .TRUE.   ;   ltwmsk(:,:) = .TRUE.   ;   ltnmsk(:,:) = .TRUE.   ;   ltsmsk(:,:) = .TRUE.
         luemsk(:,:) = .TRUE.   ;   luwmsk(:,:) = .TRUE.   ;   lunmsk(:,:) = .TRUE.   ;   lusmsk(:,:) = .TRUE.
         lvemsk(:,:) = .TRUE.   ;   lvwmsk(:,:) = .TRUE.   ;   lvnmsk(:,:) = .TRUE.   ;   lvsmsk(:,:) = .TRUE.
      END IF
      !
   END FUNCTION obc_dta_alloc


   SUBROUTINE obc_dta( kt )
      !!---------------------------------------------------------------------------
      !!                      ***  SUBROUTINE obc_dta  ***
      !!                    
      !! ** Purpose :   Find the climatological  boundary arrays for the specified date, 
      !!                The boundary arrays are netcdf files. Three possible cases: 
      !!                - one time frame only in the file (time dimension = 1).
      !!                in that case the boundary data does not change in time.
      !!                - many time frames. In that case,  if we have 12 frames
      !!                we assume monthly fields. 
      !!                Else, we assume that time_counter is in seconds 
      !!                since the beginning of either the current year or a reference
      !!                year given in the namelist.
      !!                (no check is done so far but one would have to check the "unit"
      !!                 attribute of variable time_counter).
      !!
      !!---------------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
      !
      INTEGER, SAVE :: immfile, iyyfile                     !
      INTEGER :: nt              !  record indices (incrementation)
      REAL(wp) ::   zsec, zxy, znum, zden ! time interpolation weight
      !!---------------------------------------------------------------------------

      ! 0.  initialisation :
      ! --------------------
      IF ( kt == nit000  )  CALL obc_dta_ini ( kt )
      IF ( nobc_dta == 0 )  RETURN   ! already done in obc_dta_ini
      IF ( itobc == 1    )  RETURN   ! case of only one time frame in file done in obc_dta_ini

      ! in the following code, we assume that obc data are read from files, with more than 1 time frame in it

      iyyfile=iyy ; immfile = 00  ! set component of the current file name
      IF ( cffile /= 'annual') immfile = imm   ! 
      IF ( ln_obc_clim       ) iyyfile = 0000  ! assume that climatological files are labeled y0000

      ! 1. Synchronize time of run with time of data files
      !---------------------------------------------------
      ! nday_year is the day number in the current year ( 1 for 01/01 )
      zsec=MOD( (kt-nit000)*rdt - (nday_year - nday_year0 )*rday, rday ) ! number of seconds in the current day
      IF (ln_obc_clim)  THEN 
         zjcnes = nday_year - 1  + zsec/rday
      ELSE
         zjcnes = zjcnes + rdt/rday
      ENDIF

      ! look for 'before' record number in the current file
      ntobc = nrecbef ()  ! this function return the record number for 'before', relative to zjcnes

      IF (MOD(kt-1,10)==0) THEN
         IF (lwp) WRITE(numout,*) 'kt= ',kt,' zjcnes =', zjcnes,' ndastp =',ndastp, 'mm =',imm 
      END IF

      ! 2. read a new data if necessary 
      !--------------------------------
      IF ( ntobc /= ntobc_b ) THEN
         ! we need to read the 'after' record
         ! swap working index:
# if defined key_dynspg_ts
         nt=nt_m ; nt_m=nt_b ; nt_b=nt
# endif
         nt=nt_b ; nt_b=nt_a ; nt_a=nt
         ntobc_b = ntobc

         ! new record number :
         ntobc_a = ntobc_a + 1 

         ! all tricky things related to record number, changing files etc... are managed by obc_read

         CALL obc_read (kt, nt_a, ntobc_a, iyyfile, immfile )

         ! update zjcnes_obc
# if defined key_dynspg_ts
         ntobc_m=mod(ntobc_b-2+itobc,itobc)+1
         zjcnes_obc(nt_m)= ztcobc(ntobc_m)
# endif
         zjcnes_obc(nt_b)= ztcobc(ntobc_b)
         zjcnes_obc(nt_a)= ztcobc(ntobc_a)
      ENDIF

      ! 3.   interpolation at each time step
      ! ------------------------------------
      IF( ln_obc_clim) THEN
         znum= MOD(zjcnes           - zjcnes_obc(nt_b), REAL(nyear_len(1),wp) )
         IF( znum < 0 ) znum = znum + REAL(nyear_len(1),wp)
         zden= MOD(zjcnes_obc(nt_a) - zjcnes_obc(nt_b), REAL(nyear_len(1),wp) ) 
         IF( zden < 0 ) zden = zden + REAL(nyear_len(1),wp)
      ELSE
         znum= zjcnes           - zjcnes_obc(nt_b)
         zden= zjcnes_obc(nt_a) - zjcnes_obc(nt_b)
      ENDIF
      zxy = znum / zden

      IF( lp_obc_east ) THEN
         !  fills sfoe, tfoe, ufoe ,vfoe
         sfoe(:,:) = zxy * sedta (:,:,nt_a) + (1. - zxy)*sedta(:,:,nt_b)
         tfoe(:,:) = zxy * tedta (:,:,nt_a) + (1. - zxy)*tedta(:,:,nt_b)
         ufoe(:,:) = zxy * uedta (:,:,nt_a) + (1. - zxy)*uedta(:,:,nt_b)
         vfoe(:,:) = zxy * vedta (:,:,nt_a) + (1. - zxy)*vedta(:,:,nt_b)
      ENDIF

      IF( lp_obc_west) THEN
         !  fills sfow, tfow, ufow ,vfow
         sfow(:,:) = zxy * swdta (:,:,nt_a) + (1. - zxy)*swdta(:,:,nt_b)
         tfow(:,:) = zxy * twdta (:,:,nt_a) + (1. - zxy)*twdta(:,:,nt_b)
         ufow(:,:) = zxy * uwdta (:,:,nt_a) + (1. - zxy)*uwdta(:,:,nt_b)
         vfow(:,:) = zxy * vwdta (:,:,nt_a) + (1. - zxy)*vwdta(:,:,nt_b)
      ENDIF

      IF( lp_obc_north) THEN
         !  fills sfon, tfon, ufon ,vfon
         sfon(:,:) = zxy * sndta (:,:,nt_a) + (1. - zxy)*sndta(:,:,nt_b)
         tfon(:,:) = zxy * tndta (:,:,nt_a) + (1. - zxy)*tndta(:,:,nt_b)
         ufon(:,:) = zxy * undta (:,:,nt_a) + (1. - zxy)*undta(:,:,nt_b)
         vfon(:,:) = zxy * vndta (:,:,nt_a) + (1. - zxy)*vndta(:,:,nt_b)
      ENDIF

      IF( lp_obc_south) THEN
         !  fills sfos, tfos, ufos ,vfos
         sfos(:,:) = zxy * ssdta (:,:,nt_a) + (1. - zxy)*ssdta(:,:,nt_b)
         tfos(:,:) = zxy * tsdta (:,:,nt_a) + (1. - zxy)*tsdta(:,:,nt_b)
         ufos(:,:) = zxy * usdta (:,:,nt_a) + (1. - zxy)*usdta(:,:,nt_b)
         vfos(:,:) = zxy * vsdta (:,:,nt_a) + (1. - zxy)*vsdta(:,:,nt_b)
      ENDIF
   END SUBROUTINE obc_dta


   SUBROUTINE obc_dta_ini( kt )
      !!-----------------------------------------------------------------------------
      !!                       ***  SUBROUTINE obc_dta_ini  ***
      !!
      !! ** Purpose :   When obc_dta first call, realize some data initialization
      !!----------------------------------------------------------------------------
      INTEGER, INTENT(in)  :: kt      ! ocean time-step index
      !
      INTEGER ::   ji, jj   ! dummy loop indices
      INTEGER, SAVE :: immfile, iyyfile                     !

      ! variables for the julian day calculation
      INTEGER :: iyear, imonth, iday
      REAL(wp) :: zsec , zjulian, zjuliancnes   

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*)  'obc_dta : find boundary data'
      IF(lwp) WRITE(numout,*)  '~~~~~~~'
      IF (lwp) THEN
         IF ( nobc_dta == 0 ) THEN 
            WRITE(numout,*)  '          OBC data taken from initial conditions.'
         ELSE      
            WRITE(numout,*)  '          OBC data taken from netcdf files.'
         ENDIF
      ENDIF
      nday_year0 = nday_year  ! to remember the day when kt=nit000

      sedta(:,:,:) = 0.e0 ; tedta(:,:,:) = 0.e0 ; uedta(:,:,:) = 0.e0 ; vedta(:,:,:) = 0.e0 ! East
      swdta(:,:,:) = 0.e0 ; twdta(:,:,:) = 0.e0 ; uwdta(:,:,:) = 0.e0 ; vwdta(:,:,:) = 0.e0 ! West
      sndta(:,:,:) = 0.e0 ; tndta(:,:,:) = 0.e0 ; undta(:,:,:) = 0.e0 ; vndta(:,:,:) = 0.e0 ! North
      ssdta(:,:,:) = 0.e0 ; tsdta(:,:,:) = 0.e0 ; usdta(:,:,:) = 0.e0 ; vsdta(:,:,:) = 0.e0 ! South

      sfoe(:,:) = 0.e0  ; tfoe(:,:) = 0.e0 ; ufoe(:,:) = 0.e0 ; vfoe(:,:) = 0.e0   ! East
      sfow(:,:) = 0.e0  ; tfow(:,:) = 0.e0 ; ufow(:,:) = 0.e0 ; vfow(:,:) = 0.e0   ! West
      sfon(:,:) = 0.e0  ; tfon(:,:) = 0.e0 ; ufon(:,:) = 0.e0 ; vfon(:,:) = 0.e0   ! North
      sfos(:,:) = 0.e0  ; tfos(:,:) = 0.e0 ; ufos(:,:) = 0.e0 ; vfos(:,:) = 0.e0   ! South

      IF (nobc_dta == 0 ) THEN   ! boundary data are the initial data of this run (set only at nit000)
         IF (lp_obc_east) THEN  ! East
            DO ji = nie0 , nie1    
               sfoe(nje0:nje1,:) = temsk(nje0:nje1,:) * tsn(ji+1 , nje0:nje1 , :,jp_sal) * tmask(ji+1,nje0:nje1 , :)
               tfoe(nje0:nje1,:) = temsk(nje0:nje1,:) * tsn(ji+1 , nje0:nje1 , :,jp_tem) * tmask(ji+1,nje0:nje1 , :)
               ufoe(nje0:nje1,:) = uemsk(nje0:nje1,:) * un (ji   , nje0:nje1 , :)        * umask(ji,  nje0:nje1 , :)
               vfoe(nje0:nje1,:) = vemsk(nje0:nje1,:) * vn (ji+1 , nje0:nje1 , :)        * vmask(ji+1,nje0:nje1 , :)
            END DO
         ENDIF

         IF (lp_obc_west) THEN  ! West
            DO ji = niw0 , niw1    
               sfow(njw0:njw1,:) = twmsk(njw0:njw1,:) * tsn(ji , njw0:njw1 , :,jp_sal) * tmask(ji , njw0:njw1 , :)
               tfow(njw0:njw1,:) = twmsk(njw0:njw1,:) * tsn(ji , njw0:njw1 , :,jp_tem) * tmask(ji , njw0:njw1 , :)
               ufow(njw0:njw1,:) = uwmsk(njw0:njw1,:) * un (ji , njw0:njw1 , :)        * umask(ji , njw0:njw1 , :)
               vfow(njw0:njw1,:) = vwmsk(njw0:njw1,:) * vn (ji , njw0:njw1 , :)        * vmask(ji , njw0:njw1 , :)
            END DO
         ENDIF

         IF (lp_obc_north) THEN ! North
            DO jj = njn0 , njn1
               sfon(nin0:nin1,:) = tnmsk(nin0:nin1,:) * tsn(nin0:nin1 , jj+1 , :,jp_sal) * tmask(nin0:nin1 , jj+1 , :)
               tfon(nin0:nin1,:) = tnmsk(nin0:nin1,:) * tsn(nin0:nin1 , jj+1 , :,jp_tem) * tmask(nin0:nin1 , jj+1 , :)
               ufon(nin0:nin1,:) = unmsk(nin0:nin1,:) * un (nin0:nin1 , jj+1 , :)        * umask(nin0:nin1 , jj+1 , :)
               vfon(nin0:nin1,:) = vnmsk(nin0:nin1,:) * vn (nin0:nin1 , jj   , :)        * vmask(nin0:nin1 , jj   , :)
            END DO
         ENDIF

         IF (lp_obc_south) THEN ! South
            DO jj = njs0 , njs1
               sfos(nis0:nis1,:) = tsmsk(nis0:nis1,:) * tsn(nis0:nis1 , jj , :,jp_sal) * tmask(nis0:nis1 , jj , :)
               tfos(nis0:nis1,:) = tsmsk(nis0:nis1,:) * tsn(nis0:nis1 , jj , :,jp_tem) * tmask(nis0:nis1 , jj , :)
               ufos(nis0:nis1,:) = usmsk(nis0:nis1,:) * un (nis0:nis1 , jj , :)        * umask(nis0:nis1 , jj , :)
               vfos(nis0:nis1,:) = vsmsk(nis0:nis1,:) * vn (nis0:nis1 , jj , :)        * vmask(nis0:nis1 , jj , :)
            END DO
         ENDIF
         RETURN  ! exit the routine all is done
      ENDIF  ! nobc_dta = 0 

!!!! In the following OBC data are read from files.
      ! all logical-mask are initialzed to true when declared
      WHERE ( temsk == 0 ) ltemsk=.FALSE. 
      WHERE ( uemsk == 0 ) luemsk=.FALSE. 
      WHERE ( vemsk == 0 ) lvemsk=.FALSE. 

      WHERE ( twmsk == 0 ) ltwmsk=.FALSE. 
      WHERE ( uwmsk == 0 ) luwmsk=.FALSE. 
      WHERE ( vwmsk == 0 ) lvwmsk=.FALSE. 

      WHERE ( tnmsk == 0 ) ltnmsk=.FALSE. 
      WHERE ( unmsk == 0 ) lunmsk=.FALSE. 
      WHERE ( vnmsk == 0 ) lvnmsk=.FALSE. 

      WHERE ( tsmsk == 0 ) ltsmsk=.FALSE. 
      WHERE ( usmsk == 0 ) lusmsk=.FALSE. 
      WHERE ( vsmsk == 0 ) lvsmsk=.FALSE. 

      iyear=1950;  imonth=01; iday=01;  zsec=0. 
      ! zjuliancnes : julian day corresonding  to  01/01/1950
      CALL ymds2ju(iyear, imonth, iday,zsec , zjuliancnes)

      !current year and curent month 
      iyy=INT(ndastp/10000) ; imm=INT((ndastp -iyy*10000)/100) ; idd=(ndastp-iyy*10000-imm*100)
      IF (iyy <  1900)  iyy = iyy+1900  ! always assume that years are on 4 digits.
      CALL ymds2ju(iyy, imm, idd ,zsec , zjulian)
      ndate0_cnes = zjulian - zjuliancnes   ! jcnes day when call to obc_dta_ini

      iyyfile=iyy ; immfile=0  ! set component of the current file name
      IF ( cffile /= 'annual') immfile=imm
      IF ( ln_obc_clim) iyyfile = 0  ! assume that climatological files are labeled y0000

      CALL obc_dta_chktime ( iyyfile, immfile )

      IF ( itobc == 1 ) THEN 
         ! in this case we will provide boundary data only once.
         nt_a=1 ; ntobc_a=1
         CALL obc_read (nit000, nt_a, ntobc_a, iyyfile, immfile) 
         IF( lp_obc_east ) THEN
            !  fills sfoe, tfoe, ufoe ,vfoe
            sfoe(:,:) =  sedta (:,:,1) ; tfoe(:,:) =  tedta (:,:,1)
            ufoe(:,:) =  uedta (:,:,1) ; vfoe(:,:) =  vedta (:,:,1)
         ENDIF

         IF( lp_obc_west) THEN
            !  fills sfow, tfow, ufow ,vfow
            sfow(:,:) =  swdta (:,:,1) ; tfow(:,:) =  twdta (:,:,1)
            ufow(:,:) =  uwdta (:,:,1) ; vfow(:,:) =  vwdta (:,:,1)
         ENDIF

         IF( lp_obc_north) THEN
            !  fills sfon, tfon, ufon ,vfon
            sfon(:,:) =  sndta (:,:,1) ; tfon(:,:) =  tndta (:,:,1)
            ufon(:,:) =  undta (:,:,1) ; vfon(:,:) =  vndta (:,:,1)
         ENDIF

         IF( lp_obc_south) THEN
            !  fills sfos, tfos, ufos ,vfos
            sfos(:,:) =  ssdta (:,:,1) ; tfos(:,:) =  tsdta (:,:,1)
            ufos(:,:) =  usdta (:,:,1) ; vfos(:,:) =  vsdta (:,:,1)
         ENDIF
         RETURN  ! we go out of obc_dta_ini -------------------------------------->>>>>
      ENDIF

      ! nday_year is the day number in the current year ( 1 for 01/01 )
      ! we suppose that we always start from the begining of a day
      !    zsec=MOD( (kt-nit000)*rdt - (nday_year - nday_year0 )*rday, rday ) ! number of seconds in the current day
      zsec=0.e0  ! here, kt=nit000, nday_year = ndat_year0 

      IF (ln_obc_clim)  THEN 
         zjcnes = nday_year - 1  + zsec/rday  ! for clim file time is in days in a year
      ELSE
         zjcnes = ndate0_cnes + (nday_year - nday_year0 ) + zsec/rday
      ENDIF

      ! look for 'before' record number in the current file
      ntobc = nrecbef ()

      IF (lwp) WRITE(numout,*) 'obc files frequency :',cffile
      IF (lwp) WRITE(numout,*) ' zjcnes0 =',zjcnes,' ndastp0 =',ndastp
      IF (lwp) WRITE(numout,*) ' annee0 ',iyy,' month0 ', imm,' day0 ', idd
      IF (lwp) WRITE(numout,*) 'first file open :',cl_obc_nTS

      ! record initialisation
      !--------------------
      nt_b = 1 ; nt_a = 2

      ntobc_a = ntobc + 1
      ntobc_b = ntobc

      CALL obc_read (kt, nt_b, ntobc_b, iyyfile, immfile)  ! read 'before' fields
      CALL obc_read (kt, nt_a, ntobc_a, iyyfile, immfile)  ! read 'after' fields

      ! additional frame in case of time-splitting
# if defined key_dynspg_ts
      nt_m = 0
      ntobc_m=mod(ntobc_b-2+itobc,itobc)+1
      zjcnes_obc(nt_m)= ztcobc(ntobc_m) ! FDbug has not checked that this is correct!!
      IF (ln_rstart) THEN
         CALL obc_read (kt, nt_m, ntobc_m, iyyfile, immfile)  ! read 'after' fields
      ENDIF
# endif

      zjcnes_obc(nt_b)= ztcobc(ntobc_b)
      zjcnes_obc(nt_a)= ztcobc(ntobc_a)
      ! 
   END SUBROUTINE obc_dta_ini


   SUBROUTINE obc_dta_chktime (kyyfile, kmmfile)
      !
      ! check the number of time steps in the files and read ztcobc 
      !
      ! * Arguments
      INTEGER, INTENT(in) :: kyyfile, kmmfile
      ! * local variables
      INTEGER :: istop       ! error control
      INTEGER :: ji          ! dummy loop index

      INTEGER ::  idvar, id_e, id_w, id_n, id_s       ! file identifiers
      INTEGER, DIMENSION(1)  :: itmp
      CHARACTER(LEN=25) :: cl_vname

      ntobc_a = 0; itobce =0 ; itobcw = 0; itobcn = 0; itobcs = 0
      ! build file name
      WRITE(cl_obc_eTS ,'("obc_east_TS_y",i4.4,"m",i2.2,".nc")'  ) kyyfile,kmmfile
      WRITE(cl_obc_wTS ,'("obc_west_TS_y",i4.4,"m",i2.2,".nc")'  ) kyyfile,kmmfile
      WRITE(cl_obc_nTS ,'("obc_north_TS_y",i4.4,"m",i2.2,".nc")' ) kyyfile,kmmfile
      WRITE(cl_obc_sTS ,'("obc_south_TS_y",i4.4,"m",i2.2,".nc")' ) kyyfile,kmmfile

      cl_vname = 'time_counter'
      IF ( lp_obc_east ) THEN
         CALL iom_open ( cl_obc_eTS , id_e )
         idvar = iom_varid( id_e, cl_vname, kdimsz = itmp ); itobce=itmp(1)
      ENDIF
      IF ( lp_obc_west ) THEN
         CALL iom_open ( cl_obc_wTS , id_w )
         idvar = iom_varid( id_w, cl_vname, kdimsz = itmp ) ; itobcw=itmp(1)
      ENDIF
      IF ( lp_obc_north ) THEN
         CALL iom_open ( cl_obc_nTS , id_n )
         idvar = iom_varid( id_n, cl_vname, kdimsz = itmp ) ; itobcn=itmp(1)
      ENDIF
      IF ( lp_obc_south ) THEN
         CALL iom_open ( cl_obc_sTS , id_s )
         idvar = iom_varid( id_s, cl_vname, kdimsz = itmp ) ; itobcs=itmp(1)
      ENDIF

      itobc = MAX( itobce, itobcw, itobcn, itobcs )
      istop = 0
      IF ( lp_obc_east  .AND. itobce /= itobc ) istop = istop+1 
      IF ( lp_obc_west  .AND. itobcw /= itobc ) istop = istop+1      
      IF ( lp_obc_north .AND. itobcn /= itobc ) istop = istop+1
      IF ( lp_obc_south .AND. itobcs /= itobc ) istop = istop+1 
      nstop = nstop + istop

      IF ( istop /=  0 )  THEN
         WRITE(ctmp1,*) ' east, west, north, south: ', itobce, itobcw, itobcn, itobcs
         CALL ctl_stop( 'obcdta : all files must have the same number of time steps', ctmp1 )
      ENDIF

      IF ( itobc == 1 ) THEN 
         IF (lwp) THEN
            WRITE(numout,*) ' obcdta found one time step only in the OBC files'
            IF (ln_obc_clim) THEN
               ! OK no problem
            ELSE
               ln_obc_clim=.true.
               WRITE(numout,*) ' we force ln_obc_clim to T'
            ENDIF
         ENDIF
      ELSE
         IF ( ALLOCATED(ztcobc) ) DEALLOCATE ( ztcobc )
         ALLOCATE (ztcobc(itobc))
         DO ji=1,1   ! use a dummy loop to read ztcobc only once
            IF ( lp_obc_east ) THEN
               CALL iom_gettime ( id_e, ztcobc, cl_vname ) ; CALL iom_close (id_e) ; EXIT
            ENDIF
            IF ( lp_obc_west ) THEN
               CALL iom_gettime ( id_w, ztcobc, cl_vname ) ; CALL iom_close (id_w) ; EXIT
            ENDIF
            IF ( lp_obc_north ) THEN
               CALL iom_gettime ( id_n, ztcobc, cl_vname ) ; CALL iom_close (id_n) ; EXIT
            ENDIF
            IF ( lp_obc_south ) THEN
               CALL iom_gettime ( id_s, ztcobc, cl_vname ) ; CALL iom_close (id_s) ; EXIT
            ENDIF
         END DO
         rdt_obc = ztcobc(2)-ztcobc(1)  !  just an information, not used for any computation
         IF (lwp) WRITE(numout,*) ' obcdta found', itobc,' time steps in the OBC files'
         IF (lwp) WRITE(numout,*) ' time step of obc data :', rdt_obc,' days'            
      ENDIF
      zjcnes = zjcnes - rdt/rday  ! trick : zcnes is always incremented by rdt/rday in obc_dta!
   END SUBROUTINE obc_dta_chktime

# if defined key_dynspg_ts || defined key_dynspg_exp
   SUBROUTINE obc_dta_bt( kt, kbt )
      !!---------------------------------------------------------------------------
      !!                      ***  SUBROUTINE obc_dta  ***
      !!
      !! ** Purpose :   time interpolation of barotropic data for time-splitting scheme
      !!                Data at the boundary must be in m2/s 
      !!
      !! History :  9.0  !  05-11 (V. garnier) Original code
      !!---------------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
      INTEGER, INTENT( in ) ::   kbt         ! barotropic ocean time-step index
      !
      INTEGER ::   ji, jj  ! dummy loop indices
      INTEGER ::   i15
      INTEGER ::   itobcm, itobcp
      REAL(wp) ::  zxy
      INTEGER ::   isrel           ! number of seconds since 1/1/1992
      !!---------------------------------------------------------------------------

      ! 1.   First call: check time frames available in files.
      ! -------------------------------------------------------

      IF( kt == nit000 ) THEN

         ! 1.1  Barotropic tangential velocities set to zero
         ! -------------------------------------------------
         IF( lp_obc_east  ) vbtfoe(:) = 0.e0
         IF( lp_obc_west  ) vbtfow(:) = 0.e0
         IF( lp_obc_south ) ubtfos(:) = 0.e0
         IF( lp_obc_north ) ubtfon(:) = 0.e0

         ! 1.2  Sea surface height and normal barotropic velocities set to zero
         !                               or initial conditions if nobc_dta == 0
         ! --------------------------------------------------------------------

         IF( lp_obc_east ) THEN
            ! initialisation to zero
            sshedta(:,:) = 0.e0
            ubtedta(:,:) = 0.e0
            vbtedta(:,:) = 0.e0 ! tangential component
            !                                        ! ================== !
            IF( nobc_dta == 0 )   THEN               ! initial state used !
               !                                     ! ================== !
               !  Fills sedta, tedta, uedta (global arrays)
               !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
               DO ji = nie0, nie1
                  DO jj = 1, jpj
                     sshedta(jj,1) = sshn(ji+1,jj) * tmask(ji+1,jj,1)
                  END DO
               END DO
            ENDIF
         ENDIF

         IF( lp_obc_west) THEN
            ! initialisation to zero
            sshwdta(:,:) = 0.e0
            ubtwdta(:,:) = 0.e0
            vbtwdta(:,:) = 0.e0 ! tangential component
            !                                        ! ================== !
            IF( nobc_dta == 0 )   THEN               ! initial state used !
               !                                     ! ================== !
               !  Fills swdta, twdta, uwdta (global arrays)
               !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
               DO ji = niw0, niw1
                  DO jj = 1, jpj
                     sshwdta(jj,1) = sshn(ji,jj) * tmask(ji,jj,1)
                  END DO
               END DO
            ENDIF
         ENDIF

         IF( lp_obc_north) THEN
            ! initialisation to zero
            sshndta(:,:) = 0.e0
            ubtndta(:,:) = 0.e0 ! tangential component
            vbtndta(:,:) = 0.e0
            !                                        ! ================== !
            IF( nobc_dta == 0 ) THEN                 ! initial state used !
               !                                     ! ================== !
               !  Fills sndta, tndta, vndta (global arrays)
               !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
               DO jj = njn0, njn1
                  DO ji = 1, jpi
                     sshndta(ji,1) = sshn(ji,jj+1) * tmask(ji,jj+1,1)
                  END DO
               END DO
            ENDIF
         ENDIF

         IF( lp_obc_south) THEN
            ! initialisation to zero
            sshsdta(:,:) = 0.e0
            ubtsdta(:,:) = 0.e0 ! tangential component
            vbtsdta(:,:) = 0.e0
            !                                        ! ================== !
            IF( nobc_dta == 0 )   THEN               ! initial state used !
               !                                     ! ================== !
               !  Fills ssdta, tsdta, vsdta (global arrays)
               !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
               DO jj = njs0, njs1
                  DO ji = 1, jpi
                     sshsdta(ji,1) = sshn(ji,jj) * tmask(ji,jj,1)
                  END DO
               END DO
            ENDIF
         ENDIF

         IF( nobc_dta == 0 ) CALL obc_depth_average(1)   ! depth averaged velocity from the OBC depth-dependent frames

      ENDIF        !       END kt == nit000

      !!------------------------------------------------------------------------------------
      ! 2.      Initialize the time we are at. Does this every time the routine is called,
      !         excepted when nobc_dta = 0
      !

      ! 3.  Call at every time step : Linear interpolation of BCs to current time step
      ! ----------------------------------------------------------------------

      IF( lk_dynspg_ts ) THEN
         isrel = (kt-1)*rdt + kbt*(rdt/REAL(nn_baro,wp))
      ELSE IF( lk_dynspg_exp ) THEN
         isrel=kt*rdt
      ENDIF

      itobcm = nt_b
      itobcp = nt_a
      IF( itobc == 1 .OR. nobc_dta == 0 ) THEN
         zxy = 0.e0
         itobcm = 1
         itobcp = 1
      ELSE IF( itobc == 12 ) THEN
         i15   = nday / 16
         zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
      ELSE
         zxy = (zjcnes_obc(nt_a)-FLOAT(isrel)) / (zjcnes_obc(nt_a)-zjcnes_obc(nt_b))
         IF( zxy < 0. ) THEN   ! case of extrapolation, switch to old time frames
            itobcm = nt_m
            itobcp = nt_b
            zxy = (zjcnes_obc(nt_b)-FLOAT(isrel)) / (zjcnes_obc(nt_b)-zjcnes_obc(nt_m))
         ENDIF
      ENDIF

      IF( lp_obc_east ) THEN           !  fills sshfoe, ubtfoe (local to each processor)
         DO jj = 1, jpj
            sshfoe(jj) = zxy * sshedta(jj,itobcp) + (1.-zxy) * sshedta(jj,itobcm)
            ubtfoe(jj) = zxy * ubtedta(jj,itobcp) + (1.-zxy) * ubtedta(jj,itobcm)
            vbtfoe(jj) = zxy * vbtedta(jj,itobcp) + (1.-zxy) * vbtedta(jj,itobcm)
         END DO
      ENDIF

      IF( lp_obc_west) THEN            !  fills sshfow, ubtfow (local to each processor)
         DO jj = 1, jpj
            sshfow(jj) = zxy * sshwdta(jj,itobcp) + (1.-zxy) * sshwdta(jj,itobcm)
            ubtfow(jj) = zxy * ubtwdta(jj,itobcp) + (1.-zxy) * ubtwdta(jj,itobcm)
            vbtfow(jj) = zxy * vbtwdta(jj,itobcp) + (1.-zxy) * vbtwdta(jj,itobcm)
         END DO
      ENDIF

      IF( lp_obc_north) THEN           !  fills sshfon, vbtfon (local to each processor)
         DO ji = 1, jpi
            sshfon(ji) = zxy * sshndta(ji,itobcp) + (1.-zxy) * sshndta(ji,itobcm)
            ubtfon(ji) = zxy * ubtndta(ji,itobcp) + (1.-zxy) * ubtndta(ji,itobcm)
            vbtfon(ji) = zxy * vbtndta(ji,itobcp) + (1.-zxy) * vbtndta(ji,itobcm)
         END DO
      ENDIF

      IF( lp_obc_south) THEN           !  fills sshfos, vbtfos (local to each processor)
         DO ji = 1, jpi
            sshfos(ji) = zxy * sshsdta(ji,itobcp) + (1.-zxy) * sshsdta(ji,itobcm)
            ubtfos(ji) = zxy * ubtsdta(ji,itobcp) + (1.-zxy) * ubtsdta(ji,itobcm)
            vbtfos(ji) = zxy * vbtsdta(ji,itobcp) + (1.-zxy) * vbtsdta(ji,itobcm)
         END DO
      ENDIF

   END SUBROUTINE obc_dta_bt

# else
   !!-----------------------------------------------------------------------------
   !!   Default option
   !!-----------------------------------------------------------------------------
   SUBROUTINE obc_dta_bt ( kt, kbt )       ! Empty routine
      !! * Arguments
      INTEGER,INTENT(in) :: kt
      INTEGER, INTENT( in ) ::   kbt         ! barotropic ocean time-step index
      WRITE(*,*) 'obc_dta_bt: You should not have seen this print! error?', kt
      WRITE(*,*) 'obc_dta_bt: You should not have seen this print! error?', kbt
   END SUBROUTINE obc_dta_bt
# endif

   SUBROUTINE obc_read (kt, nt_x, ntobc_x, iyy, imm)
      !!-------------------------------------------------------------------------
      !!                      ***  ROUTINE obc_read  ***
      !!
      !! ** Purpose :  Read the boundary data in files identified by iyy and imm
      !!               According to the validated open boundaries, return the 
      !!               following arrays :
      !!                sedta, tedta : East OBC salinity and temperature
      !!                uedta, vedta :   "   "  u and v velocity component      
      !!
      !!                swdta, twdta : West OBC salinity and temperature
      !!                uwdta, vwdta :   "   "  u and v velocity component      
      !!
      !!                sndta, tndta : North OBC salinity and temperature
      !!                undta, vndta :   "   "  u and v velocity component      
      !!
      !!                ssdta, tsdta : South OBC salinity and temperature
      !!                usdta, vsdta :   "   "  u and v velocity component      
      !!
      !! ** Method  :  These fields are read in the record ntobc_x of the files.
      !!               The number of records is already known. If  ntobc_x is greater
      !!               than the number of record, this routine will look for next file,
      !!               updating the indices (case of inter-annual obcs) or loop at the
      !!               begining in case of climatological file (ln_obc_clim = true ).
      !! -------------------------------------------------------------------------
      !! History:     !  2005  ( P. Mathiot, C. Langlais ) Original code
      !!              !  2008  ( J,M, Molines ) Use IOM and cleaning
      !!--------------------------------------------------------------------------

      ! * Arguments
      INTEGER, INTENT( in ) :: kt, nt_x
      INTEGER, INTENT( inout ) :: ntobc_x , iyy, imm      ! yes ! inout !

      ! * Local variables
      CHARACTER (len=40) :: &    ! file names
         cl_obc_eTS   , cl_obc_eU,  cl_obc_eV,&
         cl_obc_wTS   , cl_obc_wU,  cl_obc_wV,&
         cl_obc_nTS   , cl_obc_nU,  cl_obc_nV,&
         cl_obc_sTS   , cl_obc_sU,  cl_obc_sV

      INTEGER :: ikprint
      REAL(wp) :: zmin, zmax   ! control of boundary values

      !IOM stuff
      INTEGER :: id_e, id_w, id_n, id_s
      INTEGER, DIMENSION(2) :: istart, icount

      !--------------------------------------------------------------------------
      IF ( ntobc_x > itobc ) THEN
         IF (ln_obc_clim) THEN  ! just loop on the same file
            ntobc_x = 1 
         ELSE
            ! need to change file : it is always for an 'after' data
            IF ( cffile == 'annual' ) THEN ! go to next year file
               iyy = iyy + 1
            ELSE IF ( cffile =='monthly' ) THEN  ! go to next month file
               imm = imm + 1 
               IF ( imm == 13 ) THEN 
                  imm = 1 ; iyy = iyy + 1
               ENDIF
            ELSE
               ctmp1='obcread : this type of obc file is not supported :( '
               ctmp2=TRIM(cffile)
               CALL ctl_stop (ctmp1, ctmp2)
               ! cffile should be either annual or monthly ...
            ENDIF
            ! as the file is changed, need to update itobc etc ...
            CALL obc_dta_chktime (iyy,imm)
            ntobc_x = nrecbef() + 1 ! remember : this case occur for an after data
         ENDIF
      ENDIF

      IF( lp_obc_east ) THEN 
         ! ... Read datafile and set temperature, salinity and normal velocity
         ! ... initialise the sedta, tedta, uedta arrays
         WRITE(cl_obc_eTS ,'("obc_east_TS_y"  ,i4.4,"m",i2.2,".nc")' ) iyy,imm
         WRITE(cl_obc_eU  ,'("obc_east_U_y"   ,i4.4,"m",i2.2,".nc")' ) iyy,imm
         WRITE(cl_obc_eV  ,'("obc_east_V_y"   ,i4.4,"m",i2.2,".nc")' ) iyy,imm
         ! JMM this may change depending on the obc data format ...
         istart(:)=(/nje0+njmpp-1,1/) ; icount(:)=(/nje1-nje0 +1,jpk/)
         IF (lwp) WRITE(numout,*) 'read data in :', TRIM(cl_obc_eTS)
         IF (nje1 >= nje0 ) THEN
            CALL iom_open ( cl_obc_eTS , id_e )
            CALL iom_get ( id_e, jpdom_unknown, 'votemper', tedta(nje0:nje1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
            CALL iom_get ( id_e, jpdom_unknown, 'vosaline', sedta(nje0:nje1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
# if defined key_dynspg_ts || defined key_dynspg_exp
            CALL iom_get ( id_e, jpdom_unknown, 'vossurfh', sshedta(nje0:nje1,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
# endif
            CALL iom_close (id_e)
            !
            CALL iom_open ( cl_obc_eU , id_e )
            CALL iom_get  ( id_e, jpdom_unknown, 'vozocrtx', uedta(nje0:nje1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
            CALL iom_close ( id_e )
            !
            CALL iom_open ( cl_obc_eV , id_e )
            CALL iom_get ( id_e, jpdom_unknown, 'vomecrty', vedta(nje0:nje1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
            CALL iom_close ( id_e )

            ! mask the boundary values
            tedta(:,:,nt_x) = tedta(:,:,nt_x)*temsk(:,:) ;  sedta(:,:,nt_x) = sedta(:,:,nt_x)*temsk(:,:)
            uedta(:,:,nt_x) = uedta(:,:,nt_x)*uemsk(:,:) ;  vedta(:,:,nt_x) = vedta(:,:,nt_x)*vemsk(:,:)

            ! check any outliers 
            zmin=MINVAL( sedta(:,:,nt_x), mask=ltemsk ) ; zmax=MAXVAL(sedta(:,:,nt_x), mask=ltemsk)
            IF (  zmin < 5 .OR. zmax > 50)   THEN
               CALL ctl_stop('Error in sedta',' routine obcdta')
            ENDIF
            zmin=MINVAL( tedta(:,:,nt_x), mask=ltemsk ) ; zmax=MAXVAL(tedta(:,:,nt_x), mask=ltemsk)
            IF (  zmin < -10. .OR. zmax > 40)   THEN
               CALL ctl_stop('Error in tedta',' routine obcdta')
            ENDIF
            zmin=MINVAL( uedta(:,:,nt_x), mask=luemsk ) ; zmax=MAXVAL(uedta(:,:,nt_x), mask=luemsk)
            IF (  zmin < -5. .OR. zmax > 5.)   THEN
               CALL ctl_stop('Error in uedta',' routine obcdta')
            ENDIF
            zmin=MINVAL( vedta(:,:,nt_x), mask=lvemsk ) ; zmax=MAXVAL(vedta(:,:,nt_x), mask=lvemsk)
            IF (  zmin < -5. .OR. zmax > 5.)   THEN
               CALL ctl_stop('Error in vedta',' routine obcdta')
            ENDIF

            !               Usually printout is done only once at kt = nit000, unless nprint (namelist) > 1      
            IF ( lwp .AND.  ( kt == nit000 .OR. nprint /= 0 )  ) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' Read East OBC data records ', ntobc_x
               ikprint = jpj/20 +1
               WRITE(numout,*) ' Temperature  record 1 - printout every 3 level'
               CALL prihre( tedta(:,:,nt_x), jpj, jpk, 1, jpj, ikprint, jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Salinity  record 1 - printout every 3 level'
               CALL prihre( sedta(:,:,nt_x), jpj, jpk, 1, jpj, ikprint, jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Normal velocity U  record 1  - printout every 3 level'
               CALL prihre( uedta(:,:,nt_x), jpj, jpk, 1, jpj, ikprint, jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Tangential velocity V  record 1  - printout every 3 level'
               CALL prihre( vedta(:,:,nt_x), jpj, jpk, 1, jpj, ikprint, jpk, 1, -3, 1., numout )
            ENDIF
         ENDIF
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( lp_obc_west ) THEN
         ! ... Read datafile and set temperature, salinity and normal velocity
         ! ... initialise the swdta, twdta, uwdta arrays
         WRITE(cl_obc_wTS ,'("obc_west_TS_y"  ,i4.4,"m",i2.2,".nc")' ) iyy,imm
         WRITE(cl_obc_wU  ,'("obc_west_U_y"   ,i4.4,"m",i2.2,".nc")' ) iyy,imm
         WRITE(cl_obc_wV  ,'("obc_west_V_y"   ,i4.4,"m",i2.2,".nc")' ) iyy,imm
         istart(:)=(/njw0+njmpp-1,1/) ; icount(:)=(/njw1-njw0 +1,jpk/)
         IF (lwp) WRITE(numout,*) 'read data in :', TRIM(cl_obc_wTS)

         IF ( njw1 >= njw0 ) THEN
            CALL iom_open ( cl_obc_wTS , id_w )
            CALL iom_get ( id_w, jpdom_unknown, 'votemper', twdta(njw0:njw1,:,nt_x), & 
               &               ktime=ntobc_x , kstart=istart, kcount= icount )

            CALL iom_get ( id_w, jpdom_unknown, 'vosaline', swdta(njw0:njw1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount)
# if defined key_dynspg_ts || defined key_dynspg_exp
            CALL iom_get ( id_w, jpdom_unknown, 'vossurfh', sshwdta(njw0:njw1,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
# endif
            CALL iom_close (id_w)
            !
            CALL iom_open ( cl_obc_wU , id_w )
            CALL iom_get  ( id_w, jpdom_unknown, 'vozocrtx', uwdta(njw0:njw1,:,nt_x),&
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
            CALL iom_close ( id_w )
            !
            CALL iom_open ( cl_obc_wV , id_w )
            CALL iom_get ( id_w, jpdom_unknown, 'vomecrty', vwdta(njw0:njw1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
            CALL iom_close ( id_w )

            ! mask the boundary values
            twdta(:,:,nt_x) = twdta(:,:,nt_x)*twmsk(:,:) ;  swdta(:,:,nt_x) = swdta(:,:,nt_x)*twmsk(:,:)
            uwdta(:,:,nt_x) = uwdta(:,:,nt_x)*uwmsk(:,:) ;  vwdta(:,:,nt_x) = vwdta(:,:,nt_x)*vwmsk(:,:)

            ! check any outliers
            zmin=MINVAL( swdta(:,:,nt_x), mask=ltwmsk ) ; zmax=MAXVAL(swdta(:,:,nt_x), mask=ltwmsk)
            IF (  zmin < 5 .OR. zmax > 50)   THEN
               CALL ctl_stop('Error in swdta',' routine obcdta')
            ENDIF
            zmin=MINVAL( twdta(:,:,nt_x), mask=ltwmsk ) ; zmax=MAXVAL(twdta(:,:,nt_x), mask=ltwmsk)
            IF (  zmin < -10. .OR. zmax > 40)   THEN
               CALL ctl_stop('Error in twdta',' routine obcdta')
            ENDIF
            zmin=MINVAL( uwdta(:,:,nt_x), mask=luwmsk ) ; zmax=MAXVAL(uwdta(:,:,nt_x), mask=luwmsk)
            IF (  zmin < -5. .OR. zmax > 5.)   THEN
               CALL ctl_stop('Error in uwdta',' routine obcdta')
            ENDIF
            zmin=MINVAL( vwdta(:,:,nt_x), mask=lvwmsk ) ; zmax=MAXVAL(vwdta(:,:,nt_x), mask=lvwmsk)
            IF (  zmin < -5. .OR. zmax > 5.)   THEN
               CALL ctl_stop('Error in vwdta',' routine obcdta')
            ENDIF


            IF ( lwp .AND.  ( kt == nit000 .OR. nprint /= 0 )  ) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' Read West OBC data records ', ntobc_x
               ikprint = jpj/20 +1
               WRITE(numout,*) ' Temperature  record 1 - printout every 3 level'
               CALL prihre( twdta(:,:,nt_x), jpj, jpk, 1, jpj, ikprint, jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Salinity  record 1 - printout every 3 level'
               CALL prihre( swdta(:,:,nt_x),jpj,jpk, 1, jpj, ikprint,   jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Normal velocity U  record 1  - printout every 3 level'
               CALL prihre( uwdta(:,:,nt_x), jpj, jpk, 1, jpj, ikprint, jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Tangential velocity V  record 1  - printout every 3 level'
               CALL prihre( vwdta(:,:,nt_x), jpj, jpk, 1, jpj, ikprint, jpk, 1, -3, 1., numout )
            ENDIF
         END IF
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF( lp_obc_north) THEN
         WRITE(cl_obc_nTS ,'("obc_north_TS_y" ,i4.4,"m",i2.2,".nc")' ) iyy,imm
         WRITE(cl_obc_nV  ,'("obc_north_V_y"  ,i4.4,"m",i2.2,".nc")' ) iyy,imm
         WRITE(cl_obc_nU  ,'("obc_north_U_y"  ,i4.4,"m",i2.2,".nc")' ) iyy,imm
         istart(:)=(/nin0+nimpp-1,1/) ; icount(:)=(/nin1-nin0 +1,jpk/)
         IF (lwp) WRITE(numout,*) 'read data in :', TRIM(cl_obc_nTS)
         IF ( nin1 >= nin0 ) THEN
            CALL iom_open ( cl_obc_nTS , id_n )
            CALL iom_get ( id_n, jpdom_unknown, 'votemper', tndta(nin0:nin1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
            CALL iom_get ( id_n, jpdom_unknown, 'vosaline', sndta(nin0:nin1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
# if defined key_dynspg_ts || defined key_dynspg_exp
            CALL iom_get ( id_n, jpdom_unknown, 'vossurfh', sshndta(nin0:nin1,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
# endif
            CALL iom_close (id_n)
            !
            CALL iom_open ( cl_obc_nU , id_n )
            CALL iom_get  ( id_n, jpdom_unknown, 'vozocrtx', undta(nin0:nin1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
            CALL iom_close ( id_n )
            !
            CALL iom_open ( cl_obc_nV , id_n )
            CALL iom_get  ( id_n, jpdom_unknown, 'vomecrty', vndta(nin0:nin1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
            CALL iom_close ( id_n )

            ! mask the boundary values
            tndta(:,:,nt_x) = tndta(:,:,nt_x)*tnmsk(:,:) ;  sndta(:,:,nt_x) = sndta(:,:,nt_x)*tnmsk(:,:)
            undta(:,:,nt_x) = undta(:,:,nt_x)*unmsk(:,:) ;  vndta(:,:,nt_x) = vndta(:,:,nt_x)*vnmsk(:,:)

            ! check any outliers
            zmin=MINVAL( sndta(:,:,nt_x), mask=ltnmsk ) ; zmax=MAXVAL(sndta(:,:,nt_x), mask=ltnmsk)
            IF (  zmin < 5 .OR. zmax > 50)   THEN
               CALL ctl_stop('Error in sndta',' routine obcdta')
            ENDIF
            zmin=MINVAL( tndta(:,:,nt_x), mask=ltnmsk ) ; zmax=MAXVAL(tndta(:,:,nt_x), mask=ltnmsk)
            IF (  zmin < -10. .OR. zmax > 40)   THEN
               CALL ctl_stop('Error in tndta',' routine obcdta')
            ENDIF
            zmin=MINVAL( undta(:,:,nt_x), mask=lunmsk ) ; zmax=MAXVAL(undta(:,:,nt_x), mask=lunmsk)
            IF (  zmin < -5. .OR. zmax > 5.)   THEN
               CALL ctl_stop('Error in undta',' routine obcdta')
            ENDIF
            zmin=MINVAL( vndta(:,:,nt_x), mask=lvnmsk ) ; zmax=MAXVAL(vndta(:,:,nt_x), mask=lvnmsk)
            IF (  zmin < -5. .OR. zmax > 5.)   THEN
               CALL ctl_stop('Error in vndta',' routine obcdta')
            ENDIF

            IF ( lwp .AND.  ( kt == nit000 .OR. nprint /= 0 )  ) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' Read North OBC data records ', ntobc_x
               ikprint = jpi/20 +1
               WRITE(numout,*) ' Temperature  record 1 - printout every 3 level'
               CALL prihre( tndta(:,:,nt_x), jpi, jpk, 1, jpi, ikprint, jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Salinity  record 1 - printout every 3 level'
               CALL prihre( sndta(:,:,nt_x), jpi, jpk, 1, jpi, ikprint, jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Normal velocity V  record 1  - printout every 3 level'
               CALL prihre( vndta(:,:,nt_x), jpi, jpk, 1, jpi, ikprint, jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Tangential  velocity U  record 1  - printout every 3 level'
               CALL prihre( undta(:,:,nt_x), jpi, jpk, 1, jpi, ikprint, jpk, 1, -3, 1., numout )
            ENDIF
         ENDIF
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF( lp_obc_south) THEN 
         WRITE(cl_obc_sTS ,'("obc_south_TS_y" ,i4.4,"m",i2.2,".nc")' ) iyy,imm
         WRITE(cl_obc_sV  ,'("obc_south_V_y"  ,i4.4,"m",i2.2,".nc")' ) iyy,imm
         WRITE(cl_obc_sU  ,'("obc_south_U_y"  ,i4.4,"m",i2.2,".nc")' ) iyy,imm
         istart(:)=(/nis0+nimpp-1,1/) ; icount(:)=(/nis1-nis0 +1,jpk/)
         IF (lwp) WRITE(numout,*) 'read data in :', TRIM(cl_obc_sTS)
         IF ( nis1 >= nis0 ) THEN 
            CALL iom_open ( cl_obc_sTS , id_s )
            CALL iom_get ( id_s, jpdom_unknown, 'votemper', tsdta(nis0:nis1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
            CALL iom_get ( id_s, jpdom_unknown, 'vosaline', ssdta(nis0:nis1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
# if defined key_dynspg_ts || defined key_dynspg_exp
            CALL iom_get ( id_s, jpdom_unknown, 'vossurfh', sshsdta(nis0:nis1,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
# endif
            CALL iom_close (id_s)
            !
            CALL iom_open ( cl_obc_sU , id_s )
            CALL iom_get  ( id_s, jpdom_unknown, 'vozocrtx', usdta(nis0:nis1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
            CALL iom_close ( id_s )
            !
            CALL iom_open ( cl_obc_sV , id_s )
            CALL iom_get  ( id_s, jpdom_unknown, 'vomecrty', vsdta(nis0:nis1,:,nt_x), &
               &               ktime=ntobc_x , kstart=istart, kcount= icount )
            CALL iom_close ( id_s )

            ! mask the boundary values
            tsdta(:,:,nt_x) = tsdta(:,:,nt_x)*tsmsk(:,:) ;  ssdta(:,:,nt_x) = ssdta(:,:,nt_x)*tsmsk(:,:)
            usdta(:,:,nt_x) = usdta(:,:,nt_x)*usmsk(:,:) ;  vsdta(:,:,nt_x) = vsdta(:,:,nt_x)*vsmsk(:,:)

            ! check any outliers
            zmin=MINVAL( ssdta(:,:,nt_x), mask=ltsmsk ) ; zmax=MAXVAL(ssdta(:,:,nt_x), mask=ltsmsk)
            IF (  zmin < 5 .OR. zmax > 50)   THEN
               CALL ctl_stop('Error in ssdta',' routine obcdta')
            ENDIF
            zmin=MINVAL( tsdta(:,:,nt_x), mask=ltsmsk ) ; zmax=MAXVAL(tsdta(:,:,nt_x), mask=ltsmsk)
            IF (  zmin < -10. .OR. zmax > 40)   THEN
               CALL ctl_stop('Error in tsdta',' routine obcdta')
            ENDIF
            zmin=MINVAL( usdta(:,:,nt_x), mask=lusmsk ) ; zmax=MAXVAL(usdta(:,:,nt_x), mask=lusmsk)
            IF (  zmin < -5. .OR. zmax > 5.)   THEN
               CALL ctl_stop('Error in usdta',' routine obcdta')
            ENDIF
            zmin=MINVAL( vsdta(:,:,nt_x), mask=lvsmsk ) ; zmax=MAXVAL(vsdta(:,:,nt_x), mask=lvsmsk)
            IF (  zmin < -5. .OR. zmax > 5.)   THEN
               CALL ctl_stop('Error in vsdta',' routine obcdta')
            ENDIF

            IF ( lwp .AND.  ( kt == nit000 .OR. nprint /= 0 )  ) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' Read South OBC data records ', ntobc_x
               ikprint = jpi/20 +1
               WRITE(numout,*) ' Temperature  record 1 - printout every 3 level'
               CALL prihre( tsdta(:,:,nt_x), jpi, jpk, 1, jpi, ikprint, jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Salinity  record 1 - printout every 3 level'
               CALL prihre( ssdta(:,:,nt_x), jpi, jpk, 1, jpi, ikprint, jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Normal velocity V  record 1  - printout every 3 level'
               CALL prihre( vsdta(:,:,nt_x), jpi, jpk, 1, jpi, ikprint, jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Tangential velocity U  record 1  - printout every 3 level'
               CALL prihre( usdta(:,:,nt_x), jpi, jpk, 1, jpi, ikprint, jpk, 1, -3, 1., numout )
            ENDIF
         ENDIF
      ENDIF

# if defined key_dynspg_ts || defined key_dynspg_exp
      CALL obc_depth_average(nt_x)   ! computation of depth-averaged velocity
# endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   END SUBROUTINE obc_read


   INTEGER FUNCTION nrecbef()
      !!-----------------------------------------------------------------------
      !!                     ***    FUNCTION nrecbef   ***
      !!
      !!  Purpose : - provide the before record number in files, with respect to zjcnes
      !!
      !!    History : 2008-04 : ( J.M. Molines ) Original code
      !!-----------------------------------------------------------------------

      INTEGER :: it , idum

      idum = itobc
      DO it =1, itobc
         IF ( ztcobc(it) > zjcnes ) THEN ;  idum = it - 1 ; EXIT ;  ENDIF
         ENDDO
         ! idum can be 0 (climato, before first record)
         IF ( idum == 0 ) THEN
            IF ( ln_obc_clim ) THEN
               idum = itobc
            ELSE
               ctmp1='obc_dta: find ntobc == 0 for  non climatological file '
               ctmp2='consider adding a first record in your data file '
               CALL ctl_stop(ctmp1, ctmp2)
            ENDIF
         ENDIF
         ! idum can be itobc ( zjcnes > ztcobc (itobc) )
         !  This is not a problem ...
         nrecbef = idum

      END FUNCTION nrecbef


      SUBROUTINE obc_depth_average(nt_x)
         !!-----------------------------------------------------------------------
         !!                     ***    ROUTINE obc_depth_average   ***
         !!
         !!  Purpose : - compute the depth-averaged velocity from depth-dependent OBC frames
         !!
         !!    History : 2009-01 : ( Fred Dupont ) Original code
         !!-----------------------------------------------------------------------

         ! * Arguments
         INTEGER, INTENT( in ) :: nt_x

         ! * Local variables
         INTEGER :: ji, jj, jk


         IF( lp_obc_east ) THEN
            ! initialisation to zero
            ubtedta(:,nt_x) = 0.e0
            vbtedta(:,nt_x) = 0.e0
            DO ji = nie0, nie1
               DO jj = 1, jpj
                  DO jk = 1, jpkm1
                     ubtedta(jj,nt_x) = ubtedta(jj,nt_x) + uedta(jj,jk,nt_x)*fse3u(ji,jj,jk)
                     vbtedta(jj,nt_x) = vbtedta(jj,nt_x) + vedta(jj,jk,nt_x)*fse3v(ji+1,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF

         IF( lp_obc_west) THEN
            ! initialisation to zero
            ubtwdta(:,nt_x) = 0.e0
            vbtwdta(:,nt_x) = 0.e0
            DO ji = niw0, niw1
               DO jj = 1, jpj
                  DO jk = 1, jpkm1
                     ubtwdta(jj,nt_x) = ubtwdta(jj,nt_x) + uwdta(jj,jk,nt_x)*fse3u(ji,jj,jk)
                     vbtwdta(jj,nt_x) = vbtwdta(jj,nt_x) + vwdta(jj,jk,nt_x)*fse3v(ji,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF

         IF( lp_obc_north) THEN
            ! initialisation to zero
            ubtndta(:,nt_x) = 0.e0
            vbtndta(:,nt_x) = 0.e0
            DO jj = njn0, njn1
               DO ji = 1, jpi
                  DO jk = 1, jpkm1
                     ubtndta(ji,nt_x) = ubtndta(ji,nt_x) + undta(ji,jk,nt_x)*fse3u(ji,jj+1,jk)
                     vbtndta(ji,nt_x) = vbtndta(ji,nt_x) + vndta(ji,jk,nt_x)*fse3v(ji,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF

         IF( lp_obc_south) THEN
            ! initialisation to zero
            ubtsdta(:,nt_x) = 0.e0
            vbtsdta(:,nt_x) = 0.e0
            DO jj = njs0, njs1
               DO ji = nis0, nis1
                  DO jk = 1, jpkm1
                     ubtsdta(ji,nt_x) = ubtsdta(ji,nt_x) + usdta(ji,jk,nt_x)*fse3u(ji,jj,jk)
                     vbtsdta(ji,nt_x) = vbtsdta(ji,nt_x) + vsdta(ji,jk,nt_x)*fse3v(ji,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF

      END SUBROUTINE obc_depth_average

#else
      !!------------------------------------------------------------------------------
      !!   default option:           Dummy module          NO Open Boundary Conditions
      !!------------------------------------------------------------------------------
   CONTAINS
      SUBROUTINE obc_dta( kt )             ! Dummy routine
         INTEGER, INTENT (in) :: kt
         WRITE(*,*) 'obc_dta: You should not have seen this print! error?', kt
      END SUBROUTINE obc_dta
      !!-----------------------------------------------------------------------------
      !!   Default option
      !!-----------------------------------------------------------------------------
      SUBROUTINE obc_dta_bt ( kt, kbt )     ! Empty routine
         INTEGER,INTENT(in) :: kt
         INTEGER, INTENT( in ) ::   kbt     ! barotropic ocean time-step index
         WRITE(*,*) 'obc_dta_bt: You should not have seen this print! error?', kt
         WRITE(*,*) 'obc_dta_bt: You should not have seen this print! error?', kbt
      END SUBROUTINE obc_dta_bt
#endif
   !!==============================================================================
   END MODULE obcdta
