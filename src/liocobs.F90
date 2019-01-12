! Copyright: CNRS - Université de Grenoble
!
! Contributors : Jean-Michel Brankart, Charles-Emmanuel Testut, Laurent Parent,
!                Emmanuel Cosme, Claire Lauvernet, Frédéric Castruccio
!
! Jean-Michel.Brankart@hmg.inpg.fr
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ---                                                           ---
! ---                  LIOCOBS.F90                                ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 07-11  (J.M. Brankart)                     ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE evalhdrfilecobs  : Read individual '.cobs' file header
! --- SUBROUTINE readvalfilecobs  : Read individual '.cobs' file
! ---                               (observed values only)
! --- SUBROUTINE readcfgfilecobs  : Read individual '.cobs' file
! ---                               (observation configuration only)
! --- SUBROUTINE writehdrfilecobs : Write individual '.cobs' file header
! --- SUBROUTINE writefilecobs    : Write individual '.cobs' file
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liocobs
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrfilecobs,readvalfilecobs,readcfgfilecobs
      PUBLIC writehdrfilecobs,writefilecobs

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrfilecobs (kfninobs,indobs,jpoloc,jpitploc)
!---------------------------------------------------------------------
!
!  Purpose : Read individual '.cobs' file header
!  -------
!  Method : Read header of NetCDF observation file
!  ------
!  Input :  kfninobs : filename
!  -----    indobs   : observed variable index
!
!  Output : jpoloc   : number of observations in file
!  ------   jpitploc : number of interpolation point in obs operator
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      INTEGER, intent(in) :: indobs
      INTEGER, intent(out) :: jpoloc,jpitploc
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER, dimension(1:nbvar) :: var_dim1, var_nbr1
      CHARACTER(len=varlg), dimension(1:nbvar) :: var_nam1
      INTEGER :: jdta, inddta, dtaend1, varlg1
!
      CHARACTER(len=word80) :: kform
      INTEGER :: ierr, idf, idobs, iditp, iddta, idnam
      INTEGER :: idvnam, idvdim, idvnbr
      INTEGER, dimension(2) :: vstart,vcount
      LOGICAL :: filexists
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./evalhdrobs/evalhdrfilecobs'
         WRITE(numout,*) '    ==> READING file ',kfninobs(1:lenv(kfninobs))
      ENDIF
!
! -1.- Open NetCDF observation file
! ---------------------------------
!
      INQUIRE (FILE=kfninobs,EXIST=filexists)
      IF (.NOT.filexists) THEN
        jpoloc = 0 ; jpitploc = 0
        RETURN
      ENDIF
!
      ierr = NF90_OPEN(kfninobs,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
!
! -2.- Read and check observation file dimensions
! -----------------------------------------------
!
      ierr = NF90_INQ_DIMID(idf,'obsidx',idobs)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'itpidx',iditp)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'dta',iddta)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'namelength',idnam)
      IF (ierr.NE.0) GOTO 103
!
      ierr = NF90_INQUIRE_DIMENSION(idf,idobs,len=jpoloc)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,iditp,len=jpitploc)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,iddta,len=dtaend1)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,idnam,len=varlg1)
      IF (ierr.NE.0) GOTO 103
!
      IF (varlg1.NE.varlg) GOTO 104
      IF (dtaend1.NE.dtaend) GOTO 104
!
! -3.- Read and check observation file header
! -------------------------------------------
!
      ierr = NF90_INQ_VARID(idf,'name',idvnam)
      IF (ierr.NE.0) GOTO 105
      ierr = NF90_INQ_VARID(idf,'dim',idvdim)
      IF (ierr.NE.0) GOTO 105
      ierr = NF90_INQ_VARID(idf,'nbr',idvnbr)
      IF (ierr.NE.0) GOTO 105
!
      vstart=(/1,1/) ; vcount=(/varlg,dtaend/)
      ierr = NF90_GET_VAR(idf,idvnam,var_nam1,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
      vstart=(/1,1/) ; vcount=(/dtaend,1/)
      ierr = NF90_GET_VAR(idf,idvdim,var_dim1,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
      ierr = NF90_GET_VAR(idf,idvnbr,var_nbr1,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
!
      DO jdta = 1,dtaend
         inddta=dta_ord(jdta)
         IF (var_nam1(jdta).NE.dta_nam(inddta)) GOTO 104
         IF (var_dim1(jdta).NE.dta_dim(inddta)) GOTO 104
         IF (var_nbr1(jdta).NE.dta_nbr(inddta)) GOTO 104
      ENDDO
!
! -4.- Close observation file
! ---------------------------
!
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
!
! -5.- Control print
! ------------------
!
      IF (nprint.GE.3) THEN
         kform='(8x,a,i5)'
         WRITE(numout,kform) '- Number of observations in file: ',jpoloc
         WRITE(numout,kform) '- Number of interpolation points: ',jpitploc
      ENDIF
!
! -5.- Coherence test
! -------------------
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'liocobs','evalhdrfilecobs')
 1001 CALL printerror2(0,1001,3,'liocobs','evalhdrfilecobs')
!
 101  WRITE (texterror,*) 'Observation file does not exist',kfninobs
      CALL printerror2(0,101,3,'liocobs','evalhdrfilecobs',comment=texterror)
 102  WRITE (texterror,*) 'Bad NetCDF file: ',kfninobs
      CALL printerror2(0,102,3,'liocobs','evalhdrfilecobs',comment=texterror)
 103  WRITE (texterror,*) 'Error reading dimensions in file: ',kfninobs
      CALL printerror2(0,103,3,'liocobs','evalhdrfilecobs',comment=texterror)
 104  WRITE (texterror,*) 'Obs file inconsistent with config: ',kfninobs
      CALL printerror2(0,104,3,'liocobs','evalhdrfilecobs',comment=texterror)
 105  WRITE (texterror,*) 'Error reading header in file: ',kfninobs
      CALL printerror2(0,105,3,'liocobs','evalhdrfilecobs',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfninobs
      CALL printerror2(0,107,3,'liocobs','evalhdrfilecobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readvalfilecobs(kfninobs,kvecto,kjobs)
!---------------------------------------------------------------------
!
!  Purpose : Read individual '.cobs' file (observed values only)
!  -------
!  Method : Read NetCDF observation file
!  ------
!  Input :  kfninobs    : filename
!  -----    kjobs       : observed variable index
!
!  Output : kvecto      : vector of observed values in file
!  ------ 
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      BIGREAL, dimension(:), intent(out) :: kvecto
      INTEGER, intent(in) :: kjobs
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jposize,jpitpsize,indobs,inddbs
      CHARACTER(len=word80) :: kform
      BIGREAL4, dimension(:), allocatable :: ptabo
      BIGREAL :: sxyo_moy,sxyo_ect
      INTEGER :: jpolocobs,jpitplocobs
!
      INTEGER :: ierr, idf, idobs, iditp, iddta, idnam, idvobs
      INTEGER, dimension(2) :: vstart,vcount
      LOGICAL :: filexists
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readvalobs/readvalfilecobs'
         WRITE(numout,*) '    ==> READING file ',kfninobs(1:lenv(kfninobs))
      ENDIF
!
! Set input array sizes and observation indices
      jposize = size(kvecto,1)
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
      jpitpsize=obs_itp(indobs,inddbs)
!
! -1.- Open NetCDF observation file
! ---------------------------------
!
      INQUIRE (FILE=kfninobs,EXIST=filexists)
      IF (.NOT.filexists) GOTO 101
!
      ierr = NF90_OPEN(kfninobs,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
!
! -2.- Read and check observation file dimensions
! -----------------------------------------------
! Get dimensions ids
      ierr = NF90_INQ_DIMID(idf,'obsidx',idobs)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'itpidx',iditp)
      IF (ierr.NE.0) GOTO 103
!
! Read dimensions
      ierr = NF90_INQUIRE_DIMENSION(idf,idobs,len=jpolocobs)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,iditp,len=jpitplocobs)
      IF (ierr.NE.0) GOTO 103
!
! Control print
      IF (nprint.GE.4) THEN
         kform='(8x,a,i5)'
         WRITE(numout,kform) '- Number of observations in file: ',jpolocobs
         WRITE(numout,kform) '- Number of interpolation points: ',jpitplocobs
      ENDIF
!
! Coherence test
      IF (jpolocobs.NE.jposize) GOTO 104
      IF (jpitplocobs.NE.jpitpsize) GOTO 104
!
! -3.- Read observation values
! ----------------------------
!
      vstart=(/1,1/) ; vcount=(/jposize,jpitpsize/)
      allocate ( ptabo(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
!
! Read observation values
      ierr = NF90_INQ_VARID(idf,obs_nam(indobs,inddbs),idvobs)
      IF (ierr.NE.0) GOTO 105
      ierr = NF90_GET_VAR(idf,idvobs,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
!
! Center/reduce input data if required
      IF (lmoyect) THEN
         sxyo_moy=FREAL4(obs_moy(indobs,inddbs))
         sxyo_ect=FREAL4(obs_ect(indobs,inddbs))
         kvecto(:) = (FREAL(ptabo(:))-sxyo_moy)/sxyo_ect
      ELSE
         kvecto(:) = FREAL(ptabo(:))
      ENDIF
!
! -4.- Close observation file
! ---------------------------
!
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
!
! --- deallocation ptabo
      IF (allocated(ptabo)) deallocate (ptabo)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'liocobs','readvalfilecobs')
 1001 CALL printerror2(0,1001,3,'liocobs','readvalfilecobs')
!
 101  WRITE (texterror,*) 'Observation file does not exist',kfninobs
      CALL printerror2(0,101,3,'liocobs','readvalfilecobs',comment=texterror)
 102  WRITE (texterror,*) 'Bad NetCDF file: ',kfninobs
      CALL printerror2(0,102,3,'liocobs','readvalfilecobs',comment=texterror)
 103  WRITE (texterror,*) 'Error reading dimensions in file: ',kfninobs
      CALL printerror2(0,103,3,'liocobs','readvalfilecobs',comment=texterror)
 104  WRITE (texterror,*) 'Obs file inconsistent with config: ',kfninobs
      CALL printerror2(0,104,3,'liocobs','readvalfilecobs',comment=texterror)
 105  WRITE (texterror,*) 'Error reading obs values in file: ',kfninobs
      CALL printerror2(0,105,3,'liocobs','readvalfilecobs',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfninobs
      CALL printerror2(0,107,3,'liocobs','readvalfilecobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcfgfilecobs(kfninobs,kjobs,kflagcfg, &
     &                           kvectorms,kgridijkobs,kposcoefobs)
!---------------------------------------------------------------------
!
!  Purpose : Read individual '.cobs' file (observation configuration only)
!  -------
!  Method : Read NetCDF observation file
!  ------
!  Input :  kfninobs    : filename
!  -----    kjobs       : observed variable index
!
!  Output : kgridijkobs : observation location (x,y,z)
!  ------   kposcoefobs : observation operator (interpolation points
!                         and interpolation coefficients)
!           kvectorms   : associated error value (obsolete)
!           kflagcfg    : what element of the configuration to read
!                         (1=kvectorms, 2=kgridijkobs, 3=kposcoefobs)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      INTEGER, intent(in) :: kjobs,kflagcfg
      BIGREAL, dimension(:), optional, intent(out) :: kvectorms
      TYPE (type_gridijk), dimension(:), optional, intent(out) ::  &
     &     kgridijkobs
      TYPE (type_poscoef), dimension(:,:), optional, intent(out) ::  &
     &     kposcoefobs
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jposize,jpitpsize,indobs,inddbs
      CHARACTER(len=word80) :: kform
      BIGREAL4, dimension(:), allocatable :: ptabo
      INTEGER :: jpolocobs,jpitplocobs,jitp
!
      INTEGER :: ierr, idf, idobs, iditp, iddta, idnam
      INTEGER :: idvlon, idvlat, idvdep, idvpos, idvcoe
      INTEGER, dimension(2) :: vstart,vcount
      LOGICAL :: filexists
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readcfgobs/readcfgfilecobs'
         WRITE(numout,*) '    ==> READING file ',kfninobs(1:lenv(kfninobs))
      ENDIF
!
! Set input array sizes and observation indices
! Check coherence of input array sizes
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
      SELECT CASE (kflagcfg)
      CASE(1)
         IF (.NOT.(present(kvectorms))) GOTO 1000
         jposize = size(kvectorms,1)
         jpitpsize = obs_itp(indobs,inddbs)
      CASE(2)
         IF (.NOT.(present(kgridijkobs))) GOTO 1000
         jposize = size(kgridijkobs,1)
         jpitpsize = obs_itp(indobs,inddbs)
      CASE(3)
         IF (.NOT.(present(kposcoefobs))) GOTO 1000
         jposize = size(kposcoefobs,1)
         jpitpsize = size(kposcoefobs,2)
         IF (jpitpsize.LT.obs_itp(indobs,inddbs)) GOTO 102
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -1.- Open NetCDF observation file
! ---------------------------------
!
      INQUIRE (FILE=kfninobs,EXIST=filexists)
      IF (.NOT.filexists) GOTO 101
!
      ierr = NF90_OPEN(kfninobs,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
!
! -2.- Read and check observation file dimensions
! -----------------------------------------------
! Get dimensions ids
      ierr = NF90_INQ_DIMID(idf,'obsidx',idobs)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'itpidx',iditp)
      IF (ierr.NE.0) GOTO 103
!
! Read dimensions
      ierr = NF90_INQUIRE_DIMENSION(idf,idobs,len=jpolocobs)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,iditp,len=jpitplocobs)
      IF (ierr.NE.0) GOTO 103
!
! Control print
      IF (nprint.GE.4) THEN
         kform='(8x,a,i5)'
         WRITE(numout,kform) '- Number of observations in file: ',jpolocobs
         WRITE(numout,kform) '- Number of interpolation points: ',jpitplocobs
      ENDIF
!
! Coherence test
      IF (jpolocobs.NE.jposize) GOTO 104
      IF (jpitplocobs.NE.jpitpsize) GOTO 104
!
! -3.- Read observation file
! --------------------------
!
      vstart=(/1,1/) ; vcount=(/jposize,jpitpsize/)
      allocate ( ptabo(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
!
! Read associated error value (obsolete->zero)
      IF (kflagcfg.EQ.1) THEN
        kvectorms(:) = FREAL(0.0)
      ENDIF
!
! Read observation longitudes
      IF (kflagcfg.EQ.2) THEN
         ierr = NF90_INQ_VARID(idf,'lon',idvlon)
         IF (ierr.NE.0) GOTO 105
         ierr = NF90_GET_VAR(idf,idvlon,ptabo,start=vstart,count=vcount)
         IF (ierr.NE.0) GOTO 105
         kgridijkobs(:)%longi = FREAL(ptabo(:))
      ENDIF
!
! Read observation latitudes
      IF (kflagcfg.EQ.2) THEN
         ierr = NF90_INQ_VARID(idf,'lat',idvlat)
         IF (ierr.NE.0) GOTO 105
         ierr = NF90_GET_VAR(idf,idvlat,ptabo,start=vstart,count=vcount)
         IF (ierr.NE.0) GOTO 105
         kgridijkobs(:)%latj = FREAL(ptabo(:))
      ENDIF
!
! Read observation depths
      IF (kflagcfg.EQ.2) THEN
         ierr = NF90_INQ_VARID(idf,'depth',idvdep)
         IF (ierr.NE.0) GOTO 105
         ierr = NF90_GET_VAR(idf,idvdep,ptabo,start=vstart,count=vcount)
         IF (ierr.NE.0) GOTO 105
         kgridijkobs(:)%levk = FREAL(ptabo(:))
      ENDIF
!
! Read observation operator
      IF (kflagcfg.EQ.3) THEN
         DO jitp=1,jpitpsize
            vstart=(/1,jitp/) ; vcount=(/jposize,1/)
!
            ierr = NF90_INQ_VARID(idf,'itp_pos',idvpos)
            IF (ierr.NE.0) GOTO 105
            ierr = NF90_GET_VAR(idf,idvpos,ptabo,start=vstart,count=vcount)
            IF (ierr.NE.0) GOTO 105
            kposcoefobs(:,jitp)%pos = FREAL(ptabo(:))
!
            ierr = NF90_INQ_VARID(idf,'itp_coef',idvcoe)
            IF (ierr.NE.0) GOTO 105
            ierr = NF90_GET_VAR(idf,idvcoe,ptabo,start=vstart,count=vcount)
            IF (ierr.NE.0) GOTO 105
            kposcoefobs(:,jitp)%coef = FREAL(ptabo(:))
!     
         ENDDO
      ENDIF
!
! -4.- Close observation file
! ---------------------------
!
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
!
! --- deallocation
      IF (allocated(ptabo)) deallocate(ptabo)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'liocobs','readcfgfilecobs')
 1001 CALL printerror2(0,1001,3,'liocobs','readcfgfilecobs')
!
 101  WRITE (texterror,*) 'Observation file does not exist',kfninobs
      CALL printerror2(0,101,3,'liocobs','readcfgfilecobs',comment=texterror)
 102  WRITE (texterror,*) 'Bad NetCDF file: ',kfninobs
      CALL printerror2(0,102,3,'liocobs','readcfgfilecobs',comment=texterror)
 103  WRITE (texterror,*) 'Error reading dimensions in file: ',kfninobs
      CALL printerror2(0,103,3,'liocobs','readcfgfilecobs',comment=texterror)
 104  WRITE (texterror,*) 'Obs file inconsistent with config: ',kfninobs
      CALL printerror2(0,104,3,'liocobs','readcfgfilecobs',comment=texterror)
 105  WRITE (texterror,*) 'Error reading obs config in file: ',kfninobs
      CALL printerror2(0,105,3,'liocobs','readcfgfilecobs',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfninobs
      CALL printerror2(0,107,3,'liocobs','readcfgfilecobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writehdrfilecobs(kfnoutobs,kjobs)
!---------------------------------------------------------------------
!
!  Purpose : Write empty '.cobs' file (i.e. no file)
!  -------
!  Method : Do nothing (check that the file does not exist)
!  ------
!  Input :  kfnoutobs : filename
!  -----    kjobs    : observed variable index
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutobs
      INTEGER, intent(in) :: kjobs
      LOGICAL :: filexists
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writeobs/writehdrfilecobs'
         WRITE(numout,*) '    ==> WRITING file ',kfnoutobs(1:lenv(kfnoutobs))
      ENDIF
!
! -1.- Error if NetCDF observation file already exists
! ----------------------------------------------------
!
      INQUIRE (FILE=kfnoutobs,EXIST=filexists)
      IF (filexists) GOTO 101
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'liocobs','writehdrfilecobs')
!
 101  WRITE (texterror,*) 'Observation file already exists=',kfnoutobs
      CALL printerror2(0,101,3,'liocobs','writehdrfilecobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writefilecobs(kfnoutobs,kvecto,kvectorms,kgridijkobs, &
     &                         kposcoefobs,kjobs)
!---------------------------------------------------------------------
!
!  Purpose : Write individual '.cobs' file
!  -------
!  Method : Write direct access binary file
!  ------
!  Input : kfnoutobs   : filename
!  -----   kjobs       : observed variable index
!          kvecto      : vector of observed values in file
!          kvectorms   : associated error value (obsolete)
!          kgridijkobs : observation location (x,y,z)
!          kposcoefobs : observation operator (interpolation points
!                         and interpolation coefficients)
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutobs
      BIGREAL, dimension(:), intent(in) :: kvecto,kvectorms
      TYPE (type_gridijk), dimension(:), intent(in)  :: kgridijkobs
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
      INTEGER, intent(in) :: kjobs
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jposize,jpitpsize,indobs,inddbs
      INTEGER :: jitp
      BIGREAL4, dimension(:), allocatable :: ptabo
      BIGREAL4 :: sxyo_moy,sxyo_ect
!
      INTEGER, dimension(1:nbvar) :: var_dim1, var_nbr1
      CHARACTER(len=varlg), dimension(1:nbvar) :: var_nam1
      INTEGER :: jdta, inddta
!
      INTEGER :: ierr, idf, idobs, iditp, iddta, idnam
      INTEGER :: idvnam, idvdim, idvnbr, idvobs
      INTEGER :: idvlon, idvlat, idvdep, idvpos, idvcoe
      INTEGER, dimension(2) :: vstart,vcount
      LOGICAL :: filexists
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writeobs/writehdrfilecobs'
         WRITE(numout,*) '    ==> WRITING file ',kfnoutobs(1:lenv(kfnoutobs))
      ENDIF
!
! Set input array sizes and observation indices
      jposize = size(kvecto,1)
      jpitpsize = size(kposcoefobs,2)
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
!
! Check coherence of input array sizes
      IF (jposize.NE.size(kvectorms,1)) GOTO 102
      IF (jposize.NE.size(kgridijkobs,1)) GOTO 102
      IF (jposize.NE.size(kposcoefobs,1)) GOTO 102
      IF (jpitpsize.NE.obs_itp(indobs,inddbs)) GOTO 102
!
! -1.- Create NetCDF observation file
! -----------------------------------
!
      INQUIRE (FILE=kfnoutobs,EXIST=filexists)
      IF (filexists) GOTO 101
!
      ierr = NF90_CREATE(kfnoutobs,NF90_CLOBBER,idf)
      IF (ierr.NE.0) GOTO 102
!
! Create dimensions
      ierr = NF90_DEF_DIM(idf,'obsidx',jposize,idobs)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_DEF_DIM(idf,'itpidx',jpitpsize,iditp)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_DEF_DIM(idf,'dta',dtaend,iddta)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_DEF_DIM(idf,'namelength',varlg,idnam)
      IF (ierr.NE.0) GOTO 103
!
! Create variables
      ierr = NF90_DEF_VAR(idf,'name', &
     &                    NF90_CHAR,(/idnam,iddta/),idvnam)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'dim', &
     &                    NF90_INT,(/iddta/),idvdim)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'nbr', &
     &                    NF90_INT,(/iddta/),idvnbr)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'lon', &
     &                    NF90_FLOAT,(/idobs/),idvlon)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'lat', &
     &                    NF90_FLOAT,(/idobs/),idvlat)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'depth', &
     &                    NF90_FLOAT,(/idobs/),idvdep)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,obs_nam(indobs,inddbs), &
     &                    NF90_FLOAT,(/idobs/),idvobs)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'itp_pos', &
     &                    NF90_FLOAT,(/idobs,iditp/),idvpos)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'itp_coef', &
     &                    NF90_FLOAT,(/idobs,iditp/),idvcoe)
      IF (ierr.NE.0) GOTO 104
!
! Terminate NetCDF definition mode
      ierr = NF90_ENDDEF(idf)
      IF (ierr.NE.0) GOTO 102
!
! -2.- Write observation file header
! ----------------------------------
!
!     ierr = NF90_CLOSE(idf)
      DO jdta = 1,dtaend
         inddta=dta_ord(jdta)
         var_nam1(jdta) = dta_nam(inddta)
         var_dim1(jdta) = dta_dim(inddta)
         var_nbr1(jdta) = dta_nbr(inddta)
      ENDDO
!
      vstart=(/1,1/) ; vcount=(/varlg,dtaend/)
      ierr = NF90_PUT_VAR(idf,idvnam,values=var_nam1,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
      vstart=(/1,1/) ; vcount=(/dtaend,1/)
      ierr = NF90_PUT_VAR(idf,idvdim,values=var_dim1,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
      ierr = NF90_PUT_VAR(idf,idvnbr,values=var_nbr1,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
!
! -3.- Write observation file
! ---------------------------
!
      vstart=(/1,1/) ; vcount=(/jposize,jpitpsize/)
      allocate ( ptabo(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
!
! Write observation longitudes
      ptabo(:) = FREAL4(kgridijkobs(:)%longi)
      ierr = NF90_PUT_VAR(idf,idvlon,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
!
! Write observation latitudes
      ptabo(:) = FREAL4(kgridijkobs(:)%latj)
      ierr = NF90_PUT_VAR(idf,idvlat,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
!
! Write observation depths
      ptabo(:) = FREAL4(kgridijkobs(:)%levk)
      ierr = NF90_PUT_VAR(idf,idvdep,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
!
! Write observation values
      IF (lmoyect) THEN
         sxyo_moy=FREAL4(obs_moy(indobs,inddbs))
         sxyo_ect=FREAL4(obs_ect(indobs,inddbs))
         ptabo(:) = FREAL4(kvecto(:))*sxyo_ect+sxyo_moy
      ELSE
         ptabo(:) = FREAL4(kvecto(:))
      ENDIF
      ierr = NF90_PUT_VAR(idf,idvobs,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
!
! Write observation operator
      DO jitp=1,obs_itp(indobs,inddbs)
         vstart=(/1,jitp/) ; vcount=(/jposize,1/)
!
         ptabo(:) = FREAL4(kposcoefobs(:,jitp)%pos)
         ierr = NF90_PUT_VAR(idf,idvpos,ptabo,start=vstart,count=vcount)
         IF (ierr.NE.0) GOTO 106
!
         ptabo(:) = FREAL4(kposcoefobs(:,jitp)%coef)
         ierr = NF90_PUT_VAR(idf,idvcoe,ptabo,start=vstart,count=vcount)
         IF (ierr.NE.0) GOTO 106
      ENDDO
!
! -4.- Close observation file
! ---------------------------
!
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'liocobs','writefilecobs')
 1001 CALL printerror2(0,1001,3,'liocobs','writefilecobs')
!
 101  WRITE (texterror,*) 'Observation file already exists: ',kfnoutobs
      CALL printerror2(0,101,3,'liocobs','writefilecobs',comment=texterror)
 102  WRITE (texterror,*) 'Error creating Netcdf file: ',kfnoutobs
      CALL printerror2(0,102,1,'liocobs','writefilecobs',comment=texterror)
 103  WRITE (texterror,*) 'Error creating dimension in Netcdf file: ',kfnoutobs
      CALL printerror2(0,103,1,'liocobs','writefilecobs',comment=texterror)
 104  WRITE (texterror,*) 'Error creating variable in Netcdf file: ',kfnoutobs
      CALL printerror2(0,104,1,'liocobs','writefilecobs',comment=texterror)
 105  WRITE (texterror,*) 'Error writing header variable in file: ',kfnoutobs
      CALL printerror2(0,105,1,'liocobs','writefilecobs',comment=texterror)
 106  WRITE (texterror,*) 'Error writing obs variable in file: ',kfnoutobs
      CALL printerror2(0,106,1,'liocobs','writefilecobs',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfnoutobs
      CALL printerror2(0,107,1,'liocobs','writefilecobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liocobs
