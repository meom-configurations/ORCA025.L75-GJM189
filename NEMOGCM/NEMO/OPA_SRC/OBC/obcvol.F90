MODULE obcvol
   !!=================================================================================
   !!                       ***  MODULE  obcvol  ***
   !! Ocean dynamic :  Volume constraint when OBC and Free surface are used
   !!=================================================================================
#if   defined key_obc   &&   ! defined key_vvl
   !!---------------------------------------------------------------------------------
   !!   'key_obc'               and   NOT                 open boundary conditions
   !!   'key_vvl'                                         constant volume free surface
   !!---------------------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain 
   USE sbc_oce         ! ocean surface boundary conditions
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions
   USE lib_mpp         ! for mppsum
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC obc_vol        ! routine called by dynspg_flt

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "obc_vectopt_loop_substitute.h90"
   !!---------------------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obcvol.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!---------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE obc_vol ( kt )
      !!------------------------------------------------------------------------------
      !!                      ***  ROUTINE obcvol  ***
      !!
      !! ** Purpose : 
      !!      This routine is called in dynspg_flt to control 
      !!      the volume of the system. A correction velocity is calculated
      !!      to correct the total transport through the OBC. 
      !!      The total depth used is constant (H0) to be consistent with the 
      !!      linear free surface coded in OPA 8.2
      !!
      !! ** Method :  
      !!      The correction velocity (zubtpecor here) is defined calculating
      !!      the total transport through all open boundaries (trans_obc) minus
      !!      the cumulate E-P flux (zCflxemp) divided by the total lateral 
      !!      surface (obcsurftot) of these OBC. 
      !!
      !!      zubtpecor = [trans_obc - zCflxemp ]*(1./obcsurftot)
      !!
      !!      with zCflxemp => sum of (Evaporation minus Precipitation)
      !!                       over all the domain in m3/s at each time step.
      !!
      !!      zCflxemp < 0 when precipitation dominate
      !!      zCflxemp > 0 when evaporation dominate
      !!
      !!      There are 2 options (user's desiderata): 
      !!
      !!         1/ The volume changes according to E-P, this is the default
      !!            option. In this case the cumulate E-P flux are setting to
      !!            zero (zCflxemp=0) to calculate the correction velocity. So
      !!            it will only balance the flux through open boundaries.
      !!            (set volemp to 0 in tne namelist for this option)
      !!
      !!         2/ The volume is constant even with E-P flux. In this case
      !!            the correction velocity must balance both the flux 
      !!            through open boundaries and the ones through the free
      !!            surface. 
      !!            (set volemp to 1 in tne namelist for this option)
      !!
      !! History :
      !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Original code
      !!----------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk
      REAL(wp) ::   zubtpecor
      REAL(wp) ::   zCflxemp
      REAL(wp) ::   ztransw, ztranse, ztransn, ztranss, ztranst
      !!-----------------------------------------------------------------------------

      IF( kt == nit000 ) THEN 
         IF(lwp) WRITE(numout,*)'        '
         IF(lwp) WRITE(numout,*)'obc_vol : Correction of velocities along OBC'
         IF(lwp) WRITE(numout,*)'~~~~~~~'
         IF(lwp) WRITE(numout,*)'        '
      END IF 

      ! 1. Calculate the cumulate surface Flux zCflxemp (m3/s) over all the domain.
      ! ---------------------------------------------------------------------------

      zCflxemp = SUM ( ( emp(:,:)-rnf(:,:) )*obctmsk(:,:)* e1t(:,:) * e2t(:,:)  / rau0 ) 

      IF( lk_mpp )   CALL mpp_sum( zCflxemp )   ! sum over the global domain

      ! 2. Barotropic velocity for each open boundary
      ! ---------------------------------------------

      zubtpecor = 0.e0

      ! ... East open boundary
      IF( lp_obc_east ) THEN                      ! ... Total transport through the East OBC
         DO ji = fs_nie0, fs_nie1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  zubtpecor = zubtpecor - ua(ji,jj,jk)*e2u(ji,jj)*fse3u(ji,jj,jk) * &
             &     uemsk(jj,jk)*MAX(obctmsk(ji,jj),obctmsk(ji+1,jj) )
               END DO
            END DO
         END DO
      END IF 

      ! ... West open boundary
      IF( lp_obc_west ) THEN                      ! ... Total transport through the West OBC
         DO ji = fs_niw0, fs_niw1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  zubtpecor = zubtpecor + ua(ji,jj,jk)*e2u(ji,jj)*fse3u(ji,jj,jk) * &
             &    uwmsk(jj,jk) *MAX(obctmsk(ji,jj),obctmsk(ji+1,jj) )
               END DO
            END DO
         END DO
       ENDIF

      ! ... North open boundary
      IF( lp_obc_north ) THEN                     ! ... Total transport through the North OBC
         DO jj = fs_njn0, fs_njn1 ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  zubtpecor = zubtpecor - va(ji,jj,jk)*e1v(ji,jj)*fse3v(ji,jj,jk) * &
             &    vnmsk(ji,jk) * MAX(obctmsk(ji,jj),obctmsk(ji,jj+1) )
               END DO
            END DO
         END DO
       ENDIF

      ! ... South open boundary
      IF( lp_obc_south ) THEN                     ! ... Total transport through the South OBC
         DO jj = fs_njs0, fs_njs1 ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  zubtpecor = zubtpecor + va(ji,jj,jk)*e1v(ji,jj)*fse3v(ji,jj,jk) * &
             &    vsmsk(ji,jk) * MAX(obctmsk(ji,jj),obctmsk(ji,jj+1) )
               END DO
            END DO
         END DO
       ENDIF

      IF( lk_mpp )   CALL mpp_sum( zubtpecor )   ! sum over the global domain


      ! 3. The normal velocity correction
      ! ---------------------------------
      IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
         IF(lwp) WRITE(numout,*)'        '
         IF(lwp) WRITE(numout,*)'obc_vol : time step :', kt
         IF(lwp) WRITE(numout,*)'~~~~~~~ '
         IF(lwp) WRITE(numout,*)'          cumulate flux EMP :', zCflxemp,' (m3/s)'
         IF(lwp) WRITE(numout,*)'          lateral transport :',zubtpecor,'(m3/s)'
         IF(lwp) WRITE(numout,*)'          net inflow        :',zubtpecor-zCflxemp,'(m3/s)'
      ENDIF

      zubtpecor = (zubtpecor - zCflxemp*volemp)*(1./obcsurftot)

      IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
         IF(lwp) WRITE(numout,*)'          total lateral surface of OBC :',obcsurftot,'(m2)'
         IF(lwp) WRITE(numout,*)'          correction velocity zubtpecor :',zubtpecor,'(m/s)'
         IF(lwp) WRITE(numout,*)'        '
      END IF 

      ! 4. Correction of the total velocity on each open 
      !    boundary to respect the mass flux conservation
      ! -------------------------------------------------

      ztranse = 0.e0   ; ztransw = 0.e0 ; ztransn = 0.e0 ; ztranss = 0.e0
      ztranst = 0.e0  ! total

      IF( lp_obc_west ) THEN
         ! ... correction of the west velocity
         DO ji = fs_niw0, fs_niw1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  ua(ji,jj,jk) = ua(ji,jj,jk) - zubtpecor*uwmsk(jj,jk)
                  ztransw= ztransw + ua(ji,jj,jk)*fse3u(ji,jj,jk)*e2u(ji,jj)*uwmsk(jj,jk) * &
             &    MAX(obctmsk(ji,jj),obctmsk(ji+1,jj) )
               END DO
            END DO
         END DO

         IF( lk_mpp )   CALL mpp_sum( ztransw )   ! sum over the global domain

         IF( lwp .AND. MOD( kt, nwrite ) == 0)  WRITE(numout,*)'          West OB transport ztransw :', ztransw,'(m3/s)'
      END IF 

      IF( lp_obc_east ) THEN

         ! ... correction of the east velocity
         DO ji = fs_nie0, fs_nie1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  ua(ji,jj,jk) = ua(ji,jj,jk) + zubtpecor*uemsk(jj,jk)
                  ztranse= ztranse + ua(ji,jj,jk)*fse3u(ji,jj,jk)*e2u(ji,jj)*uemsk(jj,jk) * &
            &     MAX(obctmsk(ji,jj),obctmsk(ji+1,jj) )
               END DO
            END DO
         END DO

         IF( lk_mpp )   CALL mpp_sum( ztranse )   ! sum over the global domain

         IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
            IF(lwp) WRITE(numout,*)'          East OB transport ztranse :', ztranse,'(m3/s)'
         END IF 

      END IF 
      !{ DRAKKAR : special treatment for NATL025
      IF ( cp_cfg == 'natl' .AND. jp_cfg == 025 ) THEN
         ! no correction through the northern boundary in NATL025 ( see also obcini for obcsurftot computation)
      ELSE
         IF( lp_obc_north ) THEN
   
            ! ... correction of the north velocity
            DO jj = fs_njn0, fs_njn1 ! Vector opt.
               DO jk = 1, jpkm1
                  DO ji =  1, jpi
                     va(ji,jj,jk) = va(ji,jj,jk) + zubtpecor*vnmsk(ji,jk)
                     ztransn= ztransn + va(ji,jj,jk)*fse3v(ji,jj,jk)*e1v(ji,jj)*vnmsk(ji,jk) * &
              &      MAX(obctmsk(ji,jj),obctmsk(ji,jj+1) )
                  END DO
               END DO
            END DO
            IF( lk_mpp )   CALL mpp_sum( ztransn )   ! sum over the global domain
   
            IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
               IF(lwp) WRITE(numout,*)'          North OB transport ztransn :', ztransn,'(m3/s)'
            END IF 

         END IF 
      END IF
      !}

      IF( lp_obc_south ) THEN

         ! ... correction of the south velocity
         DO jj = fs_njs0, fs_njs1 ! Vector opt.
            DO jk = 1, jpkm1
               DO ji =  1, jpi
                  va(ji,jj,jk) = va(ji,jj,jk) - zubtpecor*vsmsk(ji,jk)
                  ztranss= ztranss + va(ji,jj,jk)*fse3v(ji,jj,jk)*e1v(ji,jj)*vsmsk(ji,jk) * &
            &     MAX(obctmsk(ji,jj),obctmsk(ji,jj+1) )
               END DO
            END DO
         END DO
         IF( lk_mpp )   CALL mpp_sum( ztranss )   ! sum over the global domain

         IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
            IF(lwp) WRITE(numout,*)'          South OB transport ztranss :', ztranss,'(m3/s)'
         END IF 

      END IF 

      ! 5. Check the cumulate transport through OBC
      !    once barotropic velocities corrected
      ! -------------------------------------------


      IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
         ztranst = ztransw - ztranse + ztranss - ztransn
         IF(lwp) WRITE(numout,*)'        '
         IF(lwp) WRITE(numout,*)'          Cumulate transport ztranst =', ztranst,'(m3/s)'
         IF(lwp) WRITE(numout,*)'          Balance  =', ztranst - zCflxemp ,'(m3/s)'
         IF(lwp) WRITE(numout,*)'        '
      END IF 

   END SUBROUTINE obc_vol

#else
   !!---------------------------------------------------------------------------------
   !!  Default option :                                                   Empty module
   !!---------------------------------------------------------------------------------
CONTAINS
   SUBROUTINE obc_vol        ! Empty routine
   END SUBROUTINE obc_vol
#endif

   !!=================================================================================
END MODULE obcvol
