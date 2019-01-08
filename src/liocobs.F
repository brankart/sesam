C Copyright: CNRS - Université de Grenoble
C
C Contributors : Jean-Michel Brankart, Charles-Emmanuel Testut, Laurent Parent,
C                Emmanuel Cosme, Claire Lauvernet, Frédéric Castruccio
C
C Jean-Michel.Brankart@hmg.inpg.fr
C
C This software is governed by the CeCILL license under French law and
C abiding by the rules of distribution of free software.  You can  use,
C modify and/ or redistribute the software under the terms of the CeCILL
C license as circulated by CEA, CNRS and INRIA at the following URL
C "http://www.cecill.info".
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C ---                                                           ---
C ---                  LIOCOBS.F                                ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 07-11  (J.M. Brankart)                     ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE evalhdrfilecobs  : Read individual '.cobs' file header
C --- SUBROUTINE readvalfilecobs  : Read individual '.cobs' file
C ---                               (observed values only)
C --- SUBROUTINE readcfgfilecobs  : Read individual '.cobs' file
C ---                               (observation configuration only)
C --- SUBROUTINE writehdrfilecobs : Write individual '.cobs' file header
C --- SUBROUTINE writefilecobs    : Write individual '.cobs' file
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liocobs
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrfilecobs,readvalfilecobs,readcfgfilecobs
      PUBLIC writehdrfilecobs,writefilecobs

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrfilecobs (kfninobs,indobs,jpoloc,jpitploc)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read individual '.cobs' file header
CCC  -------
CCC  Method : Read header of NetCDF observation file
CCC  ------
CCC  Input :  kfninobs : filename
CCC  -----    indobs   : observed variable index
CCC
CCC  Output : jpoloc   : number of observations in file
CCC  ------   jpitploc : number of interpolation point in obs operator
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
C----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      INTEGER, intent(in) :: indobs
      INTEGER, intent(out) :: jpoloc,jpitploc
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER, dimension(1:nbvar) :: var_dim1, var_nbr1
      CHARACTER(len=varlg), dimension(1:nbvar) :: var_nam1
      INTEGER :: jdta, inddta, dtaend1, varlg1
C
      CHARACTER(len=word80) :: kform
      INTEGER :: ierr, idf, idobs, iditp, iddta, idnam
      INTEGER :: idvnam, idvdim, idvnbr
      INTEGER, dimension(2) :: vstart,vcount
      LOGICAL :: filexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./evalhdrobs/evalhdrfilecobs'
         WRITE(numout,*) '    ==> READING file ',kfninobs(1:lenv(kfninobs))
      ENDIF
C
C -1.- Open NetCDF observation file
C ---------------------------------
C
      INQUIRE (FILE=kfninobs,EXIST=filexists)
      IF (.NOT.filexists) THEN
        jpoloc = 0 ; jpitploc = 0
        RETURN
      ENDIF
C
      ierr = NF90_OPEN(kfninobs,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
C
C -2.- Read and check observation file dimensions
C -----------------------------------------------
C
      ierr = NF90_INQ_DIMID(idf,'obsidx',idobs)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'itpidx',iditp)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'dta',iddta)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'namelength',idnam)
      IF (ierr.NE.0) GOTO 103
C
      ierr = NF90_INQUIRE_DIMENSION(idf,idobs,len=jpoloc)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,iditp,len=jpitploc)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,iddta,len=dtaend1)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,idnam,len=varlg1)
      IF (ierr.NE.0) GOTO 103
C
      IF (varlg1.NE.varlg) GOTO 104
      IF (dtaend1.NE.dtaend) GOTO 104
C
C -3.- Read and check observation file header
C -------------------------------------------
C
      ierr = NF90_INQ_VARID(idf,'name',idvnam)
      IF (ierr.NE.0) GOTO 105
      ierr = NF90_INQ_VARID(idf,'dim',idvdim)
      IF (ierr.NE.0) GOTO 105
      ierr = NF90_INQ_VARID(idf,'nbr',idvnbr)
      IF (ierr.NE.0) GOTO 105
C
      vstart=(/1,1/) ; vcount=(/varlg,dtaend/)
      ierr = NF90_GET_VAR(idf,idvnam,var_nam1,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
      vstart=(/1,1/) ; vcount=(/dtaend,1/)
      ierr = NF90_GET_VAR(idf,idvdim,var_dim1,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
      ierr = NF90_GET_VAR(idf,idvnbr,var_nbr1,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
C
      DO jdta = 1,dtaend
         inddta=dta_ord(jdta)
         IF (var_nam1(jdta).NE.dta_nam(inddta)) GOTO 104
         IF (var_dim1(jdta).NE.dta_dim(inddta)) GOTO 104
         IF (var_nbr1(jdta).NE.dta_nbr(inddta)) GOTO 104
      ENDDO
C
C -4.- Close observation file
C ---------------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
C -5.- Control print
C ------------------
C
      IF (nprint.GE.3) THEN
         kform='(8x,a,i5)'
         WRITE(numout,kform) '- Number of observations in file: ',jpoloc
         WRITE(numout,kform) '- Number of interpolation points: ',jpitploc
      ENDIF
C
C -5.- Coherence test
C -------------------
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocobs','evalhdrfilecobs')
 1001 CALL printerror2(0,1001,3,'liocobs','evalhdrfilecobs')
C
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
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readvalfilecobs(kfninobs,kvecto,kjobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read individual '.cobs' file (observed values only)
CCC  -------
CCC  Method : Read NetCDF observation file
CCC  ------
CCC  Input :  kfninobs    : filename
CCC  -----    kjobs       : observed variable index
CCC
CCC  Output : kvecto      : vector of observed values in file
CCC  ------ 
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      BIGREAL, dimension(:), intent(out) :: kvecto
      INTEGER, intent(in) :: kjobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jposize,jpitpsize,indobs,inddbs
      CHARACTER(len=word80) :: kform
      BIGREAL4, dimension(:), allocatable :: ptabo
      BIGREAL :: sxyo_moy,sxyo_ect
      INTEGER :: jpolocobs,jpitplocobs
C
      INTEGER :: ierr, idf, idobs, iditp, iddta, idnam, idvobs
      INTEGER, dimension(2) :: vstart,vcount
      LOGICAL :: filexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readvalobs/readvalfilecobs'
         WRITE(numout,*) '    ==> READING file ',kfninobs(1:lenv(kfninobs))
      ENDIF
C
C Set input array sizes and observation indices
      jposize = size(kvecto,1)
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
      jpitpsize=obs_itp(indobs,inddbs)
C
C -1.- Open NetCDF observation file
C ---------------------------------
C
      INQUIRE (FILE=kfninobs,EXIST=filexists)
      IF (.NOT.filexists) GOTO 101
C
      ierr = NF90_OPEN(kfninobs,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
C
C -2.- Read and check observation file dimensions
C -----------------------------------------------
C Get dimensions ids
      ierr = NF90_INQ_DIMID(idf,'obsidx',idobs)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'itpidx',iditp)
      IF (ierr.NE.0) GOTO 103
C
C Read dimensions
      ierr = NF90_INQUIRE_DIMENSION(idf,idobs,len=jpolocobs)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,iditp,len=jpitplocobs)
      IF (ierr.NE.0) GOTO 103
C
C Control print
      IF (nprint.GE.4) THEN
         kform='(8x,a,i5)'
         WRITE(numout,kform) '- Number of observations in file: ',jpolocobs
         WRITE(numout,kform) '- Number of interpolation points: ',jpitplocobs
      ENDIF
C
C Coherence test
      IF (jpolocobs.NE.jposize) GOTO 104
      IF (jpitplocobs.NE.jpitpsize) GOTO 104
C
C -3.- Read observation values
C ----------------------------
C
      vstart=(/1,1/) ; vcount=(/jposize,jpitpsize/)
      allocate ( ptabo(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
C Read observation values
      ierr = NF90_INQ_VARID(idf,obs_nam(indobs,inddbs),idvobs)
      IF (ierr.NE.0) GOTO 105
      ierr = NF90_GET_VAR(idf,idvobs,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
C
C Center/reduce input data if required
      IF (lmoyect) THEN
         sxyo_moy=FREAL4(obs_moy(indobs,inddbs))
         sxyo_ect=FREAL4(obs_ect(indobs,inddbs))
         kvecto(:) = (FREAL(ptabo(:))-sxyo_moy)/sxyo_ect
      ELSE
         kvecto(:) = FREAL(ptabo(:))
      ENDIF
C
C -4.- Close observation file
C ---------------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
C --- deallocation ptabo
      IF (allocated(ptabo)) deallocate (ptabo)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocobs','readvalfilecobs')
 1001 CALL printerror2(0,1001,3,'liocobs','readvalfilecobs')
C
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
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcfgfilecobs(kfninobs,kjobs,kflagcfg,
     $                           kvectorms,kgridijkobs,kposcoefobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read individual '.cobs' file (observation configuration only)
CCC  -------
CCC  Method : Read NetCDF observation file
CCC  ------
CCC  Input :  kfninobs    : filename
CCC  -----    kjobs       : observed variable index
CCC
CCC  Output : kgridijkobs : observation location (x,y,z)
CCC  ------   kposcoefobs : observation operator (interpolation points
CCC                         and interpolation coefficients)
CCC           kvectorms   : associated error value (obsolete)
CCC           kflagcfg    : what element of the configuration to read
CCC                         (1=kvectorms, 2=kgridijkobs, 3=kposcoefobs)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      INTEGER, intent(in) :: kjobs,kflagcfg
      BIGREAL, dimension(:), optional, intent(out) :: kvectorms
      TYPE (type_gridijk), dimension(:), optional, intent(out) :: 
     $     kgridijkobs
      TYPE (type_poscoef), dimension(:,:), optional, intent(out) :: 
     $     kposcoefobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jposize,jpitpsize,indobs,inddbs
      CHARACTER(len=word80) :: kform
      BIGREAL4, dimension(:), allocatable :: ptabo
      INTEGER :: jpolocobs,jpitplocobs,jitp
C
      INTEGER :: ierr, idf, idobs, iditp, iddta, idnam
      INTEGER :: idvlon, idvlat, idvdep, idvpos, idvcoe
      INTEGER, dimension(2) :: vstart,vcount
      LOGICAL :: filexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readcfgobs/readcfgfilecobs'
         WRITE(numout,*) '    ==> READING file ',kfninobs(1:lenv(kfninobs))
      ENDIF
C
C Set input array sizes and observation indices
C Check coherence of input array sizes
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
C
C -1.- Open NetCDF observation file
C ---------------------------------
C
      INQUIRE (FILE=kfninobs,EXIST=filexists)
      IF (.NOT.filexists) GOTO 101
C
      ierr = NF90_OPEN(kfninobs,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
C
C -2.- Read and check observation file dimensions
C -----------------------------------------------
C Get dimensions ids
      ierr = NF90_INQ_DIMID(idf,'obsidx',idobs)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'itpidx',iditp)
      IF (ierr.NE.0) GOTO 103
C
C Read dimensions
      ierr = NF90_INQUIRE_DIMENSION(idf,idobs,len=jpolocobs)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,iditp,len=jpitplocobs)
      IF (ierr.NE.0) GOTO 103
C
C Control print
      IF (nprint.GE.4) THEN
         kform='(8x,a,i5)'
         WRITE(numout,kform) '- Number of observations in file: ',jpolocobs
         WRITE(numout,kform) '- Number of interpolation points: ',jpitplocobs
      ENDIF
C
C Coherence test
      IF (jpolocobs.NE.jposize) GOTO 104
      IF (jpitplocobs.NE.jpitpsize) GOTO 104
C
C -3.- Read observation file
C --------------------------
C
      vstart=(/1,1/) ; vcount=(/jposize,jpitpsize/)
      allocate ( ptabo(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
C Read associated error value (obsolete->zero)
      IF (kflagcfg.EQ.1) THEN
        kvectorms(:) = FREAL(0.0)
      ENDIF
C
C Read observation longitudes
      IF (kflagcfg.EQ.2) THEN
         ierr = NF90_INQ_VARID(idf,'lon',idvlon)
         IF (ierr.NE.0) GOTO 105
         ierr = NF90_GET_VAR(idf,idvlon,ptabo,start=vstart,count=vcount)
         IF (ierr.NE.0) GOTO 105
         kgridijkobs(:)%longi = FREAL(ptabo(:))
      ENDIF
C
C Read observation latitudes
      IF (kflagcfg.EQ.2) THEN
         ierr = NF90_INQ_VARID(idf,'lat',idvlat)
         IF (ierr.NE.0) GOTO 105
         ierr = NF90_GET_VAR(idf,idvlat,ptabo,start=vstart,count=vcount)
         IF (ierr.NE.0) GOTO 105
         kgridijkobs(:)%latj = FREAL(ptabo(:))
      ENDIF
C
C Read observation depths
      IF (kflagcfg.EQ.2) THEN
         ierr = NF90_INQ_VARID(idf,'depth',idvdep)
         IF (ierr.NE.0) GOTO 105
         ierr = NF90_GET_VAR(idf,idvdep,ptabo,start=vstart,count=vcount)
         IF (ierr.NE.0) GOTO 105
         kgridijkobs(:)%levk = FREAL(ptabo(:))
      ENDIF
C
C Read observation operator
      IF (kflagcfg.EQ.3) THEN
         DO jitp=1,jpitpsize
            vstart=(/1,jitp/) ; vcount=(/jposize,1/)
C
            ierr = NF90_INQ_VARID(idf,'itp_pos',idvpos)
            IF (ierr.NE.0) GOTO 105
            ierr = NF90_GET_VAR(idf,idvpos,ptabo,start=vstart,count=vcount)
            IF (ierr.NE.0) GOTO 105
            kposcoefobs(:,jitp)%pos = FREAL(ptabo(:))
C
            ierr = NF90_INQ_VARID(idf,'itp_coef',idvcoe)
            IF (ierr.NE.0) GOTO 105
            ierr = NF90_GET_VAR(idf,idvcoe,ptabo,start=vstart,count=vcount)
            IF (ierr.NE.0) GOTO 105
            kposcoefobs(:,jitp)%coef = FREAL(ptabo(:))
C     
         ENDDO
      ENDIF
C
C -4.- Close observation file
C ---------------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
C --- deallocation
      IF (allocated(ptabo)) deallocate(ptabo)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocobs','readcfgfilecobs')
 1001 CALL printerror2(0,1001,3,'liocobs','readcfgfilecobs')
C
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
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writehdrfilecobs(kfnoutobs,kjobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write empty '.cobs' file (i.e. no file)
CCC  -------
CCC  Method : Do nothing (check that the file does not exist)
CCC  ------
CCC  Input :  kfnoutobs : filename
CCC  -----    kjobs    : observed variable index
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnoutobs
      INTEGER, intent(in) :: kjobs
      LOGICAL :: filexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writeobs/writehdrfilecobs'
         WRITE(numout,*) '    ==> WRITING file ',kfnoutobs(1:lenv(kfnoutobs))
      ENDIF
C
C -1.- Error if NetCDF observation file already exists
C ----------------------------------------------------
C
      INQUIRE (FILE=kfnoutobs,EXIST=filexists)
      IF (filexists) GOTO 101
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocobs','writehdrfilecobs')
C
 101  WRITE (texterror,*) 'Observation file already exists=',kfnoutobs
      CALL printerror2(0,101,3,'liocobs','writehdrfilecobs',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writefilecobs(kfnoutobs,kvecto,kvectorms,kgridijkobs,
     $                         kposcoefobs,kjobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write individual '.cobs' file
CCC  -------
CCC  Method : Write direct access binary file
CCC  ------
CCC  Input : kfnoutobs   : filename
CCC  -----   kjobs       : observed variable index
CCC          kvecto      : vector of observed values in file
CCC          kvectorms   : associated error value (obsolete)
CCC          kgridijkobs : observation location (x,y,z)
CCC          kposcoefobs : observation operator (interpolation points
CCC                         and interpolation coefficients)
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnoutobs
      BIGREAL, dimension(:), intent(in) :: kvecto,kvectorms
      TYPE (type_gridijk), dimension(:), intent(in)  :: kgridijkobs
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
      INTEGER, intent(in) :: kjobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jposize,jpitpsize,indobs,inddbs
      INTEGER :: jitp
      BIGREAL4, dimension(:), allocatable :: ptabo
      BIGREAL4 :: sxyo_moy,sxyo_ect
C
      INTEGER, dimension(1:nbvar) :: var_dim1, var_nbr1
      CHARACTER(len=varlg), dimension(1:nbvar) :: var_nam1
      INTEGER :: jdta, inddta
C
      INTEGER :: ierr, idf, idobs, iditp, iddta, idnam
      INTEGER :: idvnam, idvdim, idvnbr, idvobs
      INTEGER :: idvlon, idvlat, idvdep, idvpos, idvcoe
      INTEGER, dimension(2) :: vstart,vcount
      LOGICAL :: filexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writeobs/writehdrfilecobs'
         WRITE(numout,*) '    ==> WRITING file ',kfnoutobs(1:lenv(kfnoutobs))
      ENDIF
C
C Set input array sizes and observation indices
      jposize = size(kvecto,1)
      jpitpsize = size(kposcoefobs,2)
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
C
C Check coherence of input array sizes
      IF (jposize.NE.size(kvectorms,1)) GOTO 102
      IF (jposize.NE.size(kgridijkobs,1)) GOTO 102
      IF (jposize.NE.size(kposcoefobs,1)) GOTO 102
      IF (jpitpsize.NE.obs_itp(indobs,inddbs)) GOTO 102
C
C -1.- Create NetCDF observation file
C -----------------------------------
C
      INQUIRE (FILE=kfnoutobs,EXIST=filexists)
      IF (filexists) GOTO 101
C
      ierr = NF90_CREATE(kfnoutobs,NF90_CLOBBER,idf)
      IF (ierr.NE.0) GOTO 102
C
C Create dimensions
      ierr = NF90_DEF_DIM(idf,'obsidx',jposize,idobs)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_DEF_DIM(idf,'itpidx',jpitpsize,iditp)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_DEF_DIM(idf,'dta',dtaend,iddta)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_DEF_DIM(idf,'namelength',varlg,idnam)
      IF (ierr.NE.0) GOTO 103
C
C Create variables
      ierr = NF90_DEF_VAR(idf,'name',
     $                    NF90_CHAR,(/idnam,iddta/),idvnam)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'dim',
     $                    NF90_INT,(/iddta/),idvdim)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'nbr',
     $                    NF90_INT,(/iddta/),idvnbr)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'lon',
     $                    NF90_FLOAT,(/idobs/),idvlon)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'lat',
     $                    NF90_FLOAT,(/idobs/),idvlat)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'depth',
     $                    NF90_FLOAT,(/idobs/),idvdep)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,obs_nam(indobs,inddbs),
     $                    NF90_FLOAT,(/idobs/),idvobs)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'itp_pos',
     $                    NF90_FLOAT,(/idobs,iditp/),idvpos)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_DEF_VAR(idf,'itp_coef',
     $                    NF90_FLOAT,(/idobs,iditp/),idvcoe)
      IF (ierr.NE.0) GOTO 104
C
C Terminate NetCDF definition mode
      ierr = NF90_ENDDEF(idf)
      IF (ierr.NE.0) GOTO 102
C
C -2.- Write observation file header
C ----------------------------------
C
C     ierr = NF90_CLOSE(idf)
      DO jdta = 1,dtaend
         inddta=dta_ord(jdta)
         var_nam1(jdta) = dta_nam(inddta)
         var_dim1(jdta) = dta_dim(inddta)
         var_nbr1(jdta) = dta_nbr(inddta)
      ENDDO
C
      vstart=(/1,1/) ; vcount=(/varlg,dtaend/)
      ierr = NF90_PUT_VAR(idf,idvnam,values=var_nam1,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
      vstart=(/1,1/) ; vcount=(/dtaend,1/)
      ierr = NF90_PUT_VAR(idf,idvdim,values=var_dim1,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
      ierr = NF90_PUT_VAR(idf,idvnbr,values=var_nbr1,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 105
C
C -3.- Write observation file
C ---------------------------
C
      vstart=(/1,1/) ; vcount=(/jposize,jpitpsize/)
      allocate ( ptabo(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
C Write observation longitudes
      ptabo(:) = FREAL4(kgridijkobs(:)%longi)
      ierr = NF90_PUT_VAR(idf,idvlon,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
C
C Write observation latitudes
      ptabo(:) = FREAL4(kgridijkobs(:)%latj)
      ierr = NF90_PUT_VAR(idf,idvlat,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
C
C Write observation depths
      ptabo(:) = FREAL4(kgridijkobs(:)%levk)
      ierr = NF90_PUT_VAR(idf,idvdep,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
C
C Write observation values
      IF (lmoyect) THEN
         sxyo_moy=FREAL4(obs_moy(indobs,inddbs))
         sxyo_ect=FREAL4(obs_ect(indobs,inddbs))
         ptabo(:) = FREAL4(kvecto(:))*sxyo_ect+sxyo_moy
      ELSE
         ptabo(:) = FREAL4(kvecto(:))
      ENDIF
      ierr = NF90_PUT_VAR(idf,idvobs,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
C
C Write observation operator
      DO jitp=1,obs_itp(indobs,inddbs)
         vstart=(/1,jitp/) ; vcount=(/jposize,1/)
C
         ptabo(:) = FREAL4(kposcoefobs(:,jitp)%pos)
         ierr = NF90_PUT_VAR(idf,idvpos,ptabo,start=vstart,count=vcount)
         IF (ierr.NE.0) GOTO 106
C
         ptabo(:) = FREAL4(kposcoefobs(:,jitp)%coef)
         ierr = NF90_PUT_VAR(idf,idvcoe,ptabo,start=vstart,count=vcount)
         IF (ierr.NE.0) GOTO 106
      ENDDO
C
C -4.- Close observation file
C ---------------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocobs','writefilecobs')
 1001 CALL printerror2(0,1001,3,'liocobs','writefilecobs')
C
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
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liocobs
