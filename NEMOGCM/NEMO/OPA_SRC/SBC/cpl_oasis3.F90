MODULE cpl_oasis3
   !!======================================================================
   !!                    ***  MODULE cpl_oasis  ***
   !! Coupled O/A : coupled ocean-atmosphere case using OASIS3 V. prism_2_4
   !!               special case: NEMO OPA/LIM coupled to ECHAM5
   !!=====================================================================
   !! History :   
   !!   9.0  !  04-06  (R. Redler, NEC Laboratories Europe, Germany) Original code
   !!   " "  !  04-11  (R. Redler, NEC Laboratories Europe; N. Keenlyside, W. Park, IFM-GEOMAR, Germany) revision
   !!   " "  !  04-11  (V. Gayler, MPI M&D) Grid writing
   !!   " "  !  05-08  (R. Redler, W. Park) frld initialization, paral(2) revision
   !!   " "  !  05-09  (R. Redler) extended to allow for communication over root only
   !!   " "  !  06-01  (W. Park) modification of physical part
   !!   " "  !  06-02  (R. Redler, W. Park) buffer array fix for root exchange
   !!   3.4  !  11-11  (C. Harris) Changes to allow mutiple category fields
   !!----------------------------------------------------------------------
#if defined key_oasis3
   !!----------------------------------------------------------------------
   !!   'key_oasis3'                    coupled Ocean/Atmosphere via OASIS3
   !!----------------------------------------------------------------------
   !!   cpl_prism_init     : initialization of coupled mode communication
   !!   cpl_prism_define   : definition of grid and fields
   !!   cpl_prism_snd     : snd out fields in coupled mode
   !!   cpl_prism_rcv     : receive fields in coupled mode
   !!   cpl_prism_finalize : finalize the coupled mode communication
   !!----------------------------------------------------------------------
   USE mod_prism_proto              ! OASIS3 prism module
   USE mod_prism_def_partition_proto! OASIS3 prism module for partitioning
   USE mod_prism_put_proto          ! OASIS3 prism module for snding
   USE mod_prism_get_proto          ! OASIS3 prism module for receiving
   USE mod_comprism_proto           ! OASIS3 prism module to get coupling frequency
   USE par_oce                      ! ocean parameters
   USE dom_oce                      ! ocean space and time domain
   USE in_out_manager               ! I/O manager
   USE lbclnk                       ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   cpl_prism_init
   PUBLIC   cpl_prism_define
   PUBLIC   cpl_prism_snd
   PUBLIC   cpl_prism_rcv
   PUBLIC   cpl_prism_freq
   PUBLIC   cpl_prism_finalize

   LOGICAL, PUBLIC, PARAMETER ::   lk_cpl = .TRUE.   !: coupled flag
   INTEGER, PUBLIC            ::   OASIS_Rcv  = 1    !: return code if received field
   INTEGER, PUBLIC            ::   OASIS_idle = 0    !: return code if nothing done by oasis
   INTEGER                    ::   ncomp_id          ! id returned by prism_init_comp
   INTEGER                    ::   nerror            ! return error code

   INTEGER, PARAMETER ::   nmaxfld=40    ! Maximum number of coupling fields
   
   TYPE, PUBLIC ::   FLD_CPL               !: Type for coupling field information
      LOGICAL               ::   laction   ! To be coupled or not
      CHARACTER(len = 8)    ::   clname    ! Name of the coupling field   
      CHARACTER(len = 1)    ::   clgrid    ! Grid type  
      REAL(wp)              ::   nsgn      ! Control of the sign change
      INTEGER, DIMENSION(9) ::   nid       ! Id of the field (no more than 9 categories)
      INTEGER               ::   nct       ! Number of categories in field
   END TYPE FLD_CPL

   TYPE(FLD_CPL), DIMENSION(nmaxfld), PUBLIC ::   srcv, ssnd   !: Coupling fields

   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   exfld   ! Temporary buffer for receiving

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: cpl_oasis3.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE cpl_prism_init( kl_comm )
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_prism_init  ***
      !!
      !! ** Purpose :   Initialize coupled mode communication for ocean
      !!    exchange between AGCM, OGCM and COUPLER. (OASIS3 software)
      !!
      !! ** Method  :   OASIS3 MPI communication 
      !!--------------------------------------------------------------------
      INTEGER, INTENT(out) ::   kl_comm   ! local communicator of the model
      !!--------------------------------------------------------------------

      ! WARNING: No write in numout in this routine
      !============================================

      !------------------------------------------------------------------
      ! 1st Initialize the PRISM system for the application
      !------------------------------------------------------------------
      CALL prism_init_comp_proto ( ncomp_id, 'oceanx', nerror )
      IF ( nerror /= PRISM_Ok ) &
         CALL prism_abort_proto (ncomp_id, 'cpl_prism_init', 'Failure in prism_init_comp_proto')

      !------------------------------------------------------------------
      ! 3rd Get an MPI communicator for OPA local communication
      !------------------------------------------------------------------

      CALL prism_get_localcomm_proto ( kl_comm, nerror )
      IF ( nerror /= PRISM_Ok ) &
         CALL prism_abort_proto (ncomp_id, 'cpl_prism_init','Failure in prism_get_localcomm_proto' )
      !
   END SUBROUTINE cpl_prism_init


   SUBROUTINE cpl_prism_define( krcv, ksnd )
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_prism_define  ***
      !!
      !! ** Purpose :   Define grid and field information for ocean
      !!    exchange between AGCM, OGCM and COUPLER. (OASIS3 software)
      !!
      !! ** Method  :   OASIS3 MPI communication 
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   krcv, ksnd     ! Number of received and sent coupling fields
      !
      INTEGER :: id_part
      INTEGER :: paral(5)       ! OASIS3 box partition
      INTEGER :: ishape(2,2)    ! shape of arrays passed to PSMILe
      INTEGER :: ji,jc          ! local loop indicees
      CHARACTER(LEN=8) :: zclname
      !!--------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'cpl_prism_define : initialization in coupled ocean/atmosphere case'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)

      !
      ! ... Define the shape for the area that excludes the halo
      !     For serial configuration (key_mpp_mpi not being active)
      !     nl* is set to the global values 1 and jp*glo.
      !
      ishape(:,1) = (/ 1, nlei-nldi+1 /)
      ishape(:,2) = (/ 1, nlej-nldj+1 /)
      !
      ! ... Allocate memory for data exchange
      !
      ALLOCATE(exfld(nlei-nldi+1, nlej-nldj+1), stat = nerror)
      IF( nerror > 0 ) THEN
         CALL prism_abort_proto ( ncomp_id, 'cpl_prism_define', 'Failure in allocating exfld')   ;   RETURN
      ENDIF
      !
      ! -----------------------------------------------------------------
      ! ... Define the partition 
      ! -----------------------------------------------------------------
      
      paral(1) = 2                                              ! box partitioning
      paral(2) = jpiglo * (nldj-1+njmpp-1) + (nldi-1+nimpp-1)   ! NEMO lower left corner global offset    
      paral(3) = nlei-nldi+1                                    ! local extent in i 
      paral(4) = nlej-nldj+1                                    ! local extent in j
      paral(5) = jpiglo                                         ! global extent in x
      
      IF( ln_ctl ) THEN
         WRITE(numout,*) ' multiexchg: paral (1:5)', paral
         WRITE(numout,*) ' multiexchg: jpi, jpj =', jpi, jpj
         WRITE(numout,*) ' multiexchg: nldi, nlei, nimpp =', nldi, nlei, nimpp
         WRITE(numout,*) ' multiexchg: nldj, nlej, njmpp =', nldj, nlej, njmpp
      ENDIF
      
      CALL prism_def_partition_proto ( id_part, paral, nerror )
      !
      ! ... Announce send variables. 
      !
      DO ji = 1, ksnd
         IF ( ssnd(ji)%laction ) THEN 
            DO jc = 1, ssnd(ji)%nct
               IF ( ssnd(ji)%nct .gt. 1 ) THEN
                  WRITE(zclname,'( a7, i1)') ssnd(ji)%clname,jc
               ELSE
                  zclname=ssnd(ji)%clname
               ENDIF
               WRITE(numout,*) "Define",ji,jc,zclname," for",PRISM_Out
               CALL prism_def_var_proto (ssnd(ji)%nid(jc), zclname, id_part, (/ 2, 0/),   &
                    PRISM_Out, ishape, PRISM_REAL, nerror)
               IF ( nerror /= PRISM_Ok ) THEN
                  WRITE(numout,*) 'Failed to define transient ', ji, TRIM(zclname)
                  CALL prism_abort_proto ( ssnd(ji)%nid(jc), 'cpl_prism_define', 'Failure in prism_def_var')
               ENDIF
            END DO
         ENDIF
      END DO
      !
      ! ... Announce received variables. 
      !
      DO ji = 1, krcv
         IF ( srcv(ji)%laction ) THEN 
            DO jc = 1, srcv(ji)%nct
               IF ( srcv(ji)%nct .gt. 1 ) THEN
                  WRITE(zclname,'( a7, i1)') srcv(ji)%clname,jc
               ELSE
                  zclname=srcv(ji)%clname
               ENDIF
               WRITE(numout,*) "Define",ji,jc,zclname," for",PRISM_In
               CALL prism_def_var_proto ( srcv(ji)%nid(jc), zclname, id_part, (/ 2, 0/),   &
                    &                      PRISM_In    , ishape   , PRISM_REAL, nerror)
               IF ( nerror /= PRISM_Ok ) THEN
                  WRITE(numout,*) 'Failed to define transient ', ji, TRIM(zclname)
                  CALL prism_abort_proto ( srcv(ji)%nid(jc), 'cpl_prism_define', 'Failure in prism_def_var')
               ENDIF
            END DO
         ENDIF
      END DO
      
      !------------------------------------------------------------------
      ! End of definition phase
      !------------------------------------------------------------------
      
      CALL prism_enddef_proto(nerror)
      IF( nerror /= PRISM_Ok )   CALL prism_abort_proto ( ncomp_id, 'cpl_prism_define', 'Failure in prism_enddef')
      !
   END SUBROUTINE cpl_prism_define
   
   
   SUBROUTINE cpl_prism_snd( kid, kstep, pdata, kinfo )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_prism_snd  ***
      !!
      !! ** Purpose : - At each coupling time-step,this routine sends fields
      !!      like sst or ice cover to the coupler or remote application.
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   kid       ! variable index in the array
      INTEGER                   , INTENT(  out) ::   kinfo     ! OASIS3 info argument
      INTEGER                   , INTENT(in   ) ::   kstep     ! ocean time-step in seconds
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdata
      !!
      INTEGER                                   ::   jc        ! local loop index
      !!--------------------------------------------------------------------
      !
      ! snd data to OASIS3
      !
      DO jc = 1, ssnd(kid)%nct

         CALL prism_put_proto ( ssnd(kid)%nid(jc), kstep, pdata(nldi:nlei, nldj:nlej,jc), kinfo )
         
         IF ( ln_ctl ) THEN        
            IF ( kinfo == PRISM_Sent     .OR. kinfo == PRISM_ToRest .OR.   &
                 & kinfo == PRISM_SentOut  .OR. kinfo == PRISM_ToRestOut ) THEN
               WRITE(numout,*) '****************'
               WRITE(numout,*) 'prism_put_proto: Outgoing ', ssnd(kid)%clname
               WRITE(numout,*) 'prism_put_proto: ivarid ', ssnd(kid)%nid(jc)
               WRITE(numout,*) 'prism_put_proto:  kstep ', kstep
               WRITE(numout,*) 'prism_put_proto:   info ', kinfo
               WRITE(numout,*) '     - Minimum value is ', MINVAL(pdata(:,:,jc))
               WRITE(numout,*) '     - Maximum value is ', MAXVAL(pdata(:,:,jc))
               WRITE(numout,*) '     -     Sum value is ', SUM(pdata(:,:,jc))
               WRITE(numout,*) '****************'
            ENDIF
         ENDIF

      ENDDO
      !
    END SUBROUTINE cpl_prism_snd


   SUBROUTINE cpl_prism_rcv( kid, kstep, pdata, kinfo )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_prism_rcv  ***
      !!
      !! ** Purpose : - At each coupling time-step,this routine receives fields
      !!      like stresses and fluxes from the coupler or remote application.
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   kid       ! variable index in the array
      INTEGER                   , INTENT(in   ) ::   kstep     ! ocean time-step in seconds
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pdata     ! IN to keep the value if nothing is done
      INTEGER                   , INTENT(  out) ::   kinfo     ! OASIS3 info argument
      !!
      INTEGER                                   ::   jc        ! local loop index
      LOGICAL                                   ::   llaction
      !!--------------------------------------------------------------------
      !
      ! receive local data from OASIS3 on every process
      !
      DO jc = 1, srcv(kid)%nct

         CALL prism_get_proto ( srcv(kid)%nid(jc), kstep, exfld, kinfo )         
         
         llaction = .false.
         IF( kinfo == PRISM_Recvd   .OR. kinfo == PRISM_FromRest .OR.   &
              kinfo == PRISM_RecvOut .OR. kinfo == PRISM_FromRestOut )   llaction = .TRUE.
         
         IF ( ln_ctl )   WRITE(numout,*) "llaction, kinfo, kstep, ivarid: " , llaction, kinfo, kstep, srcv(kid)%nid(jc)
         
         IF ( llaction ) THEN
            
            kinfo = OASIS_Rcv
            pdata(nldi:nlei, nldj:nlej,jc) = exfld(:,:)
            
            !--- Fill the overlap areas and extra hallows (mpp)
            !--- check periodicity conditions (all cases)
            CALL lbc_lnk( pdata(:,:,jc), srcv(kid)%clgrid, srcv(kid)%nsgn )   
            
            IF ( ln_ctl ) THEN        
               WRITE(numout,*) '****************'
               WRITE(numout,*) 'prism_get_proto: Incoming ', srcv(kid)%clname
               WRITE(numout,*) 'prism_get_proto: ivarid '  , srcv(kid)%nid(jc)
               WRITE(numout,*) 'prism_get_proto:   kstep', kstep
               WRITE(numout,*) 'prism_get_proto:   info ', kinfo
               WRITE(numout,*) '     - Minimum value is ', MINVAL(pdata(:,:,jc))
               WRITE(numout,*) '     - Maximum value is ', MAXVAL(pdata(:,:,jc))
               WRITE(numout,*) '     -     Sum value is ', SUM(pdata(:,:,jc))
               WRITE(numout,*) '****************'
            ENDIF
            
         ELSE
            kinfo = OASIS_idle     
         ENDIF
         
      ENDDO
      !
   END SUBROUTINE cpl_prism_rcv


   INTEGER FUNCTION cpl_prism_freq( kid )  
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_prism_freq  ***
      !!
      !! ** Purpose : - send back the coupling frequency for a particular field
      !!----------------------------------------------------------------------
      INTEGER,INTENT(in) ::   kid   ! variable index 
      !!----------------------------------------------------------------------
      cpl_prism_freq = ig_def_freq( kid )
      !
   END FUNCTION cpl_prism_freq


   SUBROUTINE cpl_prism_finalize
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_prism_finalize  ***
      !!
      !! ** Purpose : - Finalizes the coupling. If MPI_init has not been
      !!      called explicitly before cpl_prism_init it will also close
      !!      MPI communication.
      !!----------------------------------------------------------------------
      !
      DEALLOCATE( exfld )
      CALL prism_terminate_proto( nerror )         
      !
   END SUBROUTINE cpl_prism_finalize

#else
   !!----------------------------------------------------------------------
   !!   Default case          Dummy module          Forced Ocean/Atmosphere
   !!----------------------------------------------------------------------
   USE in_out_manager               ! I/O manager
   LOGICAL, PUBLIC, PARAMETER :: lk_cpl = .FALSE.   !: coupled flag
   PUBLIC cpl_prism_init
   PUBLIC cpl_prism_finalize
CONTAINS
   SUBROUTINE cpl_prism_init (kl_comm) 
      INTEGER, INTENT(out)   :: kl_comm       ! local communicator of the model
      kl_comm = -1
      WRITE(numout,*) 'cpl_prism_init: Error you sould not be there...'
   END SUBROUTINE cpl_prism_init
   SUBROUTINE cpl_prism_finalize
      WRITE(numout,*) 'cpl_prism_finalize: Error you sould not be there...'
   END SUBROUTINE cpl_prism_finalize
#endif

   !!=====================================================================
END MODULE cpl_oasis3
