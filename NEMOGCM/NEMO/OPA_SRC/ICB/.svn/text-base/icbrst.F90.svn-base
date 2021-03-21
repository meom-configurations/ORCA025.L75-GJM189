MODULE icbrst

   !!======================================================================
   !!                       ***  MODULE  icbrst  ***
   !! Ocean physics:  read and write iceberg restart files
   !!======================================================================
   !! History : 3.3.1 !  2010-01  (Martin&Adcroft) Original code
   !!            -    !  2011-03  (Madec)          Part conversion to NEMO form
   !!            -    !                            Removal of mapping from another grid
   !!            -    !  2011-04  (Alderson)       Split into separate modules
   !!            -    !  2011-04  (Alderson)       Restore restart routine
   !!            -    !                            Currently needs a fixed processor
   !!            -    !                            layout between restarts
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   icb_rst_read    : read restart file
   !!   icb_rst_write   : write restart file
   !!----------------------------------------------------------------------
   USE par_oce        ! NEMO parameters
   USE dom_oce        ! NEMO domain
   USE in_out_manager ! NEMO IO routines
   USE lib_mpp        ! NEMO MPI library, lk_mpp in particular
   USE netcdf         ! netcdf routines for IO
   USE icb_oce        ! define iceberg arrays
   USE icbutl         ! iceberg utility routines

   IMPLICIT NONE
   PRIVATE

   PUBLIC   icb_rst_read    ! routine called in icbini.F90 module
   PUBLIC   icb_rst_write   ! routine called in icbstp.F90 module
   
   INTEGER ::   nlonid, nlatid, nxid, nyid, nuvelid, nvvelid
   INTEGER ::   nmassid, nthicknessid, nwidthid, nlengthid
   INTEGER ::   nyearid, ndayid
   INTEGER ::   nscaling_id, nmass_of_bits_id, nheat_density_id, numberid
   INTEGER ::   nsiceid, nsheatid, ncalvid, ncalvhid, nkountid
   INTEGER ::   nret, ncid, nc_dim
   
   INTEGER,  DIMENSION(3)                  :: nstrt3, nlngth3

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2011)
   !! $Id:$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE icb_rst_read()
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE icb_rst_read  ***
      !!
      !! ** Purpose :   read a iceberg restart file
      !!      NB: for this version, we just read back in the restart for this processor
      !!      so we cannot change the processor layout currently with iceberg code
      !!----------------------------------------------------------------------
      INTEGER                      ::   idim, ivar, iatt
      INTEGER                      ::   jn, iunlim_dim, ibergs_in_file
      INTEGER                      ::   iclass
      INTEGER, DIMENSION(1)        ::   istrt, ilngth, idata
      INTEGER, DIMENSION(2)        ::   istrt2, ilngth2
      INTEGER, DIMENSION(nkounts)  ::   idata2
      REAL(wp), DIMENSION(1)       ::   zdata                                         ! need 1d array to read in with
                                                                                            ! start and count arrays
      LOGICAL                      ::   ll_found_restart
      CHARACTER(len=80)            ::   cl_filename
      CHARACTER(len=NF90_MAX_NAME) ::   cl_dname
      TYPE(iceberg)                ::   localberg ! NOT a pointer but an actual local variable
      TYPE(point)                  ::   localpt   ! NOT a pointer but an actual local variable
      !!----------------------------------------------------------------------

      ! Find a restart file
      cl_filename = ' '
      IF ( lk_mpp ) THEN
         cl_filename = ' '
         WRITE( cl_filename, '("restart_icebergs_",I4.4,".nc")' ) narea-1
         INQUIRE( file=TRIM(cl_filename), exist=ll_found_restart )
      ELSE
         cl_filename = 'restart_icebergs.nc'
         INQUIRE( file=TRIM(cl_filename), exist=ll_found_restart )
      ENDIF

      IF ( .NOT. ll_found_restart) THEN                     ! only do the following if a file was found
         CALL ctl_stop('icebergs: no restart file found')
      ENDIF

      IF (nn_verbose_level >= 0 .AND. lwp)  &
         WRITE(numout,'(2a)') 'icebergs, read_restart_bergs: found restart file = ',TRIM(cl_filename)

      nret = NF90_OPEN(TRIM(cl_filename), NF90_NOWRITE, ncid)
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, read_restart_bergs: nf_open failed')

      nret = nf90_inquire(ncid, idim, ivar, iatt, iunlim_dim)
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, read_restart_bergs: nf_inquire failed')

      IF( iunlim_dim .NE. -1) THEN

         nret = nf90_inquire_dimension(ncid, iunlim_dim, cl_dname, ibergs_in_file)
         IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, read_restart_bergs: nf_inq_dimlen failed')

         nret = NF90_INQ_VARID(ncid, 'number', numberid)
         nret = NF90_INQ_VARID(ncid, 'mass_scaling', nscaling_id)
         nret = NF90_INQ_VARID(ncid, 'xi', nxid)
         nret = NF90_INQ_VARID(ncid, 'yj', nyid)
         nret = NF90_INQ_VARID(ncid, 'lon', nlonid)
         nret = NF90_INQ_VARID(ncid, 'lat', nlatid)
         nret = NF90_INQ_VARID(ncid, 'uvel', nuvelid)
         nret = NF90_INQ_VARID(ncid, 'vvel', nvvelid)
         nret = NF90_INQ_VARID(ncid, 'mass', nmassid)
         nret = NF90_INQ_VARID(ncid, 'thickness', nthicknessid)
         nret = NF90_INQ_VARID(ncid, 'width', nwidthid)
         nret = NF90_INQ_VARID(ncid, 'length', nlengthid)
         nret = NF90_INQ_VARID(ncid, 'year', nyearid)
         nret = NF90_INQ_VARID(ncid, 'day', ndayid)
         nret = NF90_INQ_VARID(ncid, 'mass_of_bits', nmass_of_bits_id)
         nret = NF90_INQ_VARID(ncid, 'heat_density', nheat_density_id)

         ilngth(1) = 1
         istrt2(1) = 1
         ilngth2(1) = nkounts
         ilngth2(2) = 1
         DO jn=1, ibergs_in_file

            istrt(1) = jn
            istrt2(2) = jn

            nret = NF90_GET_VAR(ncid, numberid, idata2, istrt2, ilngth2 )
            localberg%number(:) = idata2(:)

            nret = NF90_GET_VAR(ncid, nscaling_id, zdata, istrt, ilngth )
            localberg%mass_scaling = zdata(1)

            nret = NF90_GET_VAR(ncid, nlonid, zdata, istrt, ilngth)
            localpt%lon = zdata(1)
            nret = NF90_GET_VAR(ncid, nlatid, zdata, istrt, ilngth)
            localpt%lat = zdata(1)
            IF (nn_verbose_level >= 2 .AND. lwp) THEN
               WRITE(numout,'(a,i5,a,2f10.4,a,i5)') 'icebergs, read_restart_bergs: berg ',jn,' is at ', &
                                              localpt%lon,localpt%lat,' on PE ',narea-1
            ENDIF
            nret = NF90_GET_VAR(ncid, nxid, zdata, istrt, ilngth)
            localpt%xi = zdata(1)
            nret = NF90_GET_VAR(ncid, nyid, zdata, istrt, ilngth)
            localpt%yj = zdata(1)
            nret = NF90_GET_VAR(ncid, nuvelid, zdata, istrt, ilngth )
            localpt%uvel = zdata(1)
            nret = NF90_GET_VAR(ncid, nvvelid, zdata, istrt, ilngth )
            localpt%vvel = zdata(1)
            nret = NF90_GET_VAR(ncid, nmassid, zdata, istrt, ilngth )
            localpt%mass = zdata(1)
            nret = NF90_GET_VAR(ncid, nthicknessid, zdata, istrt, ilngth )
            localpt%thickness = zdata(1)
            nret = NF90_GET_VAR(ncid, nwidthid, zdata, istrt, ilngth )
            localpt%width = zdata(1)
            nret = NF90_GET_VAR(ncid, nlengthid, zdata, istrt, ilngth )
            localpt%length = zdata(1)
            nret = NF90_GET_VAR(ncid, nyearid, idata, istrt, ilngth )
            localpt%year = idata(1)
            nret = NF90_GET_VAR(ncid, ndayid, zdata, istrt, ilngth )
            localpt%day = zdata(1)
            nret = NF90_GET_VAR(ncid, nmass_of_bits_id, zdata, istrt, ilngth )
            localpt%mass_of_bits = zdata(1)
            nret = NF90_GET_VAR(ncid, nheat_density_id, zdata, istrt, ilngth )
            localpt%heat_density = zdata(1)
            !
            CALL icb_utl_add( localberg, localpt )
         END DO
         !
      ENDIF

      nret = NF90_INQ_DIMID( ncid, 'c', nc_dim )
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, read_restart: nf_inq_dimid c failed')

      nret = NF90_INQUIRE_DIMENSION( ncid, nc_dim, cl_dname, iclass )
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, read_restart: nf_inquire_dimension failed')

      nret = NF90_INQ_VARID(ncid, 'kount'       , nkountid)
      nret = NF90_INQ_VARID(ncid, 'calving'     , ncalvid)
      nret = NF90_INQ_VARID(ncid, 'calving_hflx', ncalvhid)
      nret = NF90_INQ_VARID(ncid, 'stored_ice'  , nsiceid)
      nret = NF90_INQ_VARID(ncid, 'stored_heat' , nsheatid)

      nstrt3(1) = 1
      nstrt3(2) = 1
      nlngth3(1) = jpi
      nlngth3(2) = jpj
      nlngth3(3) = 1

      DO jn = 1, iclass
         nstrt3(3) = jn
         nret      = NF90_GET_VAR( ncid, nsiceid , griddata, nstrt3, nlngth3 )
         berg_grid%stored_ice(:,:,jn) = griddata(:,:,1)
      END DO

      nret = NF90_GET_VAR( ncid, ncalvid , src_calving          (:,:) )
      nret = NF90_GET_VAR( ncid, ncalvhid, src_calving_hflx     (:,:) )
      nret = NF90_GET_VAR( ncid, nsheatid, berg_grid%stored_heat(:,:) )
      nret = NF90_GET_VAR( ncid, nkountid, idata2(:) )
      num_bergs(:) = idata2(:)

      ! Finish up
      nret = NF90_CLOSE(ncid)
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, read_restart: nf_close failed')

      ! Sanity check
      jn = icb_utl_count()
      IF (nn_verbose_level >= 0)   &
         WRITE(numout,'(2(a,i5))') 'icebergs, read_restart_bergs: # bergs =',jn,' on PE',narea-1
      IF( lk_mpp ) THEN
         CALL mpp_sum(ibergs_in_file)
         CALL mpp_sum(jn)
      ENDIF
      IF(lwp)   WRITE(numout,'(a,i5,a,i5,a)') 'icebergs, read_restart_bergs: there were',ibergs_in_file,   &
         &                                    ' bergs in the restart file and', jn,' bergs have been read'
      !
      IF( lwp .and. nn_verbose_level >= 0)  WRITE(numout,'(a)') 'icebergs, read_restart_bergs: completed'
      !
   END SUBROUTINE icb_rst_read


   SUBROUTINE icb_rst_write( kt )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE icb_rst_write  ***
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) :: kt
      !
      INTEGER ::   jn   ! dummy loop index
      INTEGER ::   ix_dim, iy_dim, ik_dim, in_dim
      CHARACTER(len=80)      :: cl_filename
      TYPE(iceberg), POINTER :: this
      TYPE(point)  , POINTER :: pt
      !!----------------------------------------------------------------------

      IF( lk_mpp ) THEN
         WRITE(cl_filename,'("icebergs_",I8.8,"_restart_",I4.4,".nc")') kt, narea-1
      ELSE
         WRITE(cl_filename,'("icebergs_",I8.8,"_restart.nc")') kt
      ENDIF
      IF (nn_verbose_level >= 0) WRITE(numout,'(2a)') 'icebergs, write_restart: creating ',TRIM(cl_filename)

      nret = NF90_CREATE(TRIM(cl_filename), NF90_CLOBBER, ncid)
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_create failed')

      ! Dimensions
      nret = NF90_DEF_DIM(ncid, 'x', jpi, ix_dim)
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_def_dim x failed')

      nret = NF90_DEF_DIM(ncid, 'y', jpj, iy_dim)
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_def_dim y failed')

      nret = NF90_DEF_DIM(ncid, 'c', nclasses, nc_dim)
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_def_dim c failed')

      nret = NF90_DEF_DIM(ncid, 'k', nkounts, ik_dim)
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_def_dim k failed')

      IF (associated(first_berg)) then
         nret = NF90_DEF_DIM(ncid, 'n', NF90_UNLIMITED, in_dim)
         IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_def_dim n failed')
      ENDIF

      ! Variables
      nret = NF90_DEF_VAR(ncid, 'kount'       , NF90_INT   , (/ ik_dim /), nkountid)
      nret = NF90_DEF_VAR(ncid, 'calving'     , NF90_DOUBLE, (/ ix_dim, iy_dim /), ncalvid)
      nret = NF90_DEF_VAR(ncid, 'calving_hflx', NF90_DOUBLE, (/ ix_dim, iy_dim /), ncalvhid)
      nret = NF90_DEF_VAR(ncid, 'stored_ice'  , NF90_DOUBLE, (/ ix_dim, iy_dim, nc_dim /), nsiceid)
      nret = NF90_DEF_VAR(ncid, 'stored_heat' , NF90_DOUBLE, (/ ix_dim, iy_dim /), nsheatid)

      ! Attributes
      nret = NF90_PUT_ATT(ncid, ncalvid , 'long_name', 'iceberg calving')
      nret = NF90_PUT_ATT(ncid, ncalvid , 'units', 'some')
      nret = NF90_PUT_ATT(ncid, ncalvhid, 'long_name', 'heat flux associated with iceberg calving')
      nret = NF90_PUT_ATT(ncid, ncalvhid, 'units', 'some')
      nret = NF90_PUT_ATT(ncid, nsiceid , 'long_name', 'stored ice used to calve icebergs')
      nret = NF90_PUT_ATT(ncid, nsiceid , 'units', 'kg/s')
      nret = NF90_PUT_ATT(ncid, nsheatid, 'long_name', 'heat in stored ice used to calve icebergs')
      nret = NF90_PUT_ATT(ncid, nsheatid, 'units', 'J/kg/s')

      IF ( ASSOCIATED(first_berg) ) THEN

         ! Only add berg variables for this PE if we have anything to say

         ! Variables
         nret = NF90_DEF_VAR(ncid, 'lon', NF90_DOUBLE, in_dim, nlonid)
         nret = NF90_DEF_VAR(ncid, 'lat', NF90_DOUBLE, in_dim, nlatid)
         nret = NF90_DEF_VAR(ncid, 'xi', NF90_DOUBLE, in_dim, nxid)
         nret = NF90_DEF_VAR(ncid, 'yj', NF90_DOUBLE, in_dim, nyid)
         nret = NF90_DEF_VAR(ncid, 'uvel', NF90_DOUBLE, in_dim, nuvelid)
         nret = NF90_DEF_VAR(ncid, 'vvel', NF90_DOUBLE, in_dim, nvvelid)
         nret = NF90_DEF_VAR(ncid, 'mass', NF90_DOUBLE, in_dim, nmassid)
         nret = NF90_DEF_VAR(ncid, 'thickness', NF90_DOUBLE, in_dim, nthicknessid)
         nret = NF90_DEF_VAR(ncid, 'width', NF90_DOUBLE, in_dim, nwidthid)
         nret = NF90_DEF_VAR(ncid, 'length', NF90_DOUBLE, in_dim, nlengthid)
         nret = NF90_DEF_VAR(ncid, 'number', NF90_INT, (/ik_dim,in_dim/), numberid)
         nret = NF90_DEF_VAR(ncid, 'year', NF90_INT, in_dim, nyearid)
         nret = NF90_DEF_VAR(ncid, 'day', NF90_DOUBLE, in_dim, ndayid)
         nret = NF90_DEF_VAR(ncid, 'mass_scaling', NF90_DOUBLE, in_dim, nscaling_id)
         nret = NF90_DEF_VAR(ncid, 'mass_of_bits', NF90_DOUBLE, in_dim, nmass_of_bits_id)
         nret = NF90_DEF_VAR(ncid, 'heat_density', NF90_DOUBLE, in_dim, nheat_density_id)

         ! Attributes
         nret = NF90_PUT_ATT(ncid, nlonid, 'long_name', 'longitude')
         nret = NF90_PUT_ATT(ncid, nlonid, 'units', 'degrees_E')
         nret = NF90_PUT_ATT(ncid, nlatid, 'long_name', 'latitude')
         nret = NF90_PUT_ATT(ncid, nlatid, 'units', 'degrees_N')
         nret = NF90_PUT_ATT(ncid, nxid, 'long_name', 'x grid box position')
         nret = NF90_PUT_ATT(ncid, nxid, 'units', 'fractional')
         nret = NF90_PUT_ATT(ncid, nyid, 'long_name', 'y grid box position')
         nret = NF90_PUT_ATT(ncid, nyid, 'units', 'fractional')
         nret = NF90_PUT_ATT(ncid, nuvelid, 'long_name', 'zonal velocity')
         nret = NF90_PUT_ATT(ncid, nuvelid, 'units', 'm/s')
         nret = NF90_PUT_ATT(ncid, nvvelid, 'long_name', 'meridional velocity')
         nret = NF90_PUT_ATT(ncid, nvvelid, 'units', 'm/s')
         nret = NF90_PUT_ATT(ncid, nmassid, 'long_name', 'mass')
         nret = NF90_PUT_ATT(ncid, nmassid, 'units', 'kg')
         nret = NF90_PUT_ATT(ncid, nthicknessid, 'long_name', 'thickness')
         nret = NF90_PUT_ATT(ncid, nthicknessid, 'units', 'm')
         nret = NF90_PUT_ATT(ncid, nwidthid, 'long_name', 'width')
         nret = NF90_PUT_ATT(ncid, nwidthid, 'units', 'm')
         nret = NF90_PUT_ATT(ncid, nlengthid, 'long_name', 'length')
         nret = NF90_PUT_ATT(ncid, nlengthid, 'units', 'm')
         nret = NF90_PUT_ATT(ncid, numberid, 'long_name', 'iceberg number on this processor')
         nret = NF90_PUT_ATT(ncid, numberid, 'units', 'count')
         nret = NF90_PUT_ATT(ncid, nyearid, 'long_name', 'calendar year of calving event')
         nret = NF90_PUT_ATT(ncid, nyearid, 'units', 'years')
         nret = NF90_PUT_ATT(ncid, ndayid, 'long_name', 'year day of calving event')
         nret = NF90_PUT_ATT(ncid, ndayid, 'units', 'days')
         nret = NF90_PUT_ATT(ncid, nscaling_id, 'long_name', 'scaling factor for mass of calving berg')
         nret = NF90_PUT_ATT(ncid, nscaling_id, 'units', 'none')
         nret = NF90_PUT_ATT(ncid, nmass_of_bits_id, 'long_name', 'mass of bergy bits')
         nret = NF90_PUT_ATT(ncid, nmass_of_bits_id, 'units', 'kg')
         nret = NF90_PUT_ATT(ncid, nheat_density_id, 'long_name', 'heat density')
         nret = NF90_PUT_ATT(ncid, nheat_density_id, 'units', 'J/kg')

      ENDIF ! associated(first_berg)

      ! End define mode
      nret = NF90_ENDDEF(ncid)

      ! --------------------------------
      ! now write some data

      nstrt3(1) = 1
      nstrt3(2) = 1
      nlngth3(1) = jpi
      nlngth3(2) = jpj
      nlngth3(3) = 1

      DO jn=1,nclasses
         griddata(:,:,1) = berg_grid%stored_ice(:,:,jn)
         nstrt3(3) = jn
         nret = NF90_PUT_VAR( ncid, nsiceid, griddata, nstrt3, nlngth3 )
         IF (nret .ne. NF90_NOERR) THEN
            IF( lwp ) WRITE(numout,*) TRIM(NF90_STRERROR( nret ))
            CALL ctl_stop('icebergs, write_restart: nf_put_var stored_ice failed')
         ENDIF
      ENDDO
      IF( lwp ) WRITE(numout,*) 'file: ',TRIM(cl_filename),' var: stored_ice  written'

      nret = NF90_PUT_VAR( ncid, nkountid, num_bergs(:) )
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_put_var kount failed')

      nret = NF90_PUT_VAR( ncid, nsheatid, berg_grid%stored_heat(:,:) )
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_put_var stored_heat failed')
      IF( lwp ) WRITE(numout,*) 'file: ',TRIM(cl_filename),' var: stored_heat written'

      nret = NF90_PUT_VAR( ncid, ncalvid , src_calving(:,:) )
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_put_var calving failed')
      nret = NF90_PUT_VAR( ncid, ncalvhid, src_calving_hflx(:,:) )
      IF (nret .ne. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_put_var calving_hflx failed')
      IF( lwp ) WRITE(numout,*) 'file: ',TRIM(cl_filename),' var: calving written'

      IF ( ASSOCIATED(first_berg) ) THEN

         ! Write variables
         ! just write out the current point of the trajectory

         this => first_berg
         jn = 0
         DO WHILE (ASSOCIATED(this))
            pt => this%current_point
            jn=jn+1

            nret = NF90_PUT_VAR(ncid, numberid, this%number, (/1,jn/), (/nkounts,1/) )
            nret = NF90_PUT_VAR(ncid, nscaling_id, this%mass_scaling, (/ jn /) )

            nret = NF90_PUT_VAR(ncid, nlonid, pt%lon, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, nlatid, pt%lat, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, nxid, pt%xi, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, nyid, pt%yj, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, nuvelid, pt%uvel, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, nvvelid, pt%vvel, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, nmassid, pt%mass, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, nthicknessid, pt%thickness, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, nwidthid, pt%width, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, nlengthid, pt%length, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, nyearid, pt%year, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, ndayid, pt%day, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, nmass_of_bits_id, pt%mass_of_bits, (/ jn /) )
            nret = NF90_PUT_VAR(ncid, nheat_density_id, pt%heat_density, (/ jn /) )

            this=>this%next
         END DO
         !
      ENDIF ! associated(first_berg)

      ! Finish up
      nret = NF90_CLOSE(ncid)
      IF (nret /= NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_close failed')
      !
   END SUBROUTINE icb_rst_write
   !
END MODULE icbrst
