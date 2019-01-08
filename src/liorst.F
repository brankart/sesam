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
C ---                   LIORST.F                                ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12  (C.E. Testut)                       ---
C --- modification : 99-05  (C.E. Testut)                       ---
C --- modification : 01-06  (C.E. Testut)                       ---
C --- modification : 03-04  (J.M. Brankart)                     ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE  readrst  : Read Vx or Vy vector from user defined format
C --- SUBROUTINE  writerst : Write Vx vector in user defined format
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liorst
      use mod_main
      use utilvct
      use utilfiles
      IMPLICIT NONE
      PRIVATE

      PUBLIC readrst,writerst

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readrst(kfninsxy,kvects,klectinfo,kflagxyo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read Vx or Vy vector from user defined format
CCC  -------
CCC  Method : Loop on variable fields defined in SESAM configuration,
CCC  ------   read them in input file, and store them in Vx or Vy vector
CCC
CCC  Input :  kfninsxy  : filename
CCC  -----    klectinfo : read or not header of input files
CCC           kflagxyo  : vector type (Vx/Vy/Vo)
CCC
CCC  Output : kvects : 1D vector object (Vx or Vy)
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : jpxend, jpyend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfninsxy
      BIGREAL, dimension(:), intent(out) :: kvects
      LOGICAL, intent(in) :: klectinfo
      INTEGER, intent(in) :: kflagxyo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=bgword) :: kform
      CHARACTER(len=varlg), dimension(1:nbvar) :: sxy_nam
      INTEGER :: jpsend,sxyend
      INTEGER, dimension(1:nbvar) :: sxy_ord,sxy_dim,
     $     sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt,sxy_nbr,sxy_ind
      INTEGER :: jsxy,indsxy,somsxynbr,sompartsxynbr
      INTEGER :: nsdeb, nsfin, nstot
      INTEGER :: allocok,jpisize,jpjsize,jpksize,jpssize
      BIGREAL8, allocatable, dimension(:,:) :: ptabij
      INTEGER :: jinf,ji,jj,jk,jt,js
      INTEGER :: nrec,iblocrst
CC----------------------------------------------------------------------
      jpssize=size(kvects,1)
CCC---------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.3) THEN
         WRITE(numout,*) '*** ROUTINE : ../readvar/readrst'
         WRITE(numout,*) '    ==> READING file ',
     $            kfninsxy(1:lenv(kfninsxy))
      ENDIF
C
C Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE(1)
C --- var
         jpsend  = jpxend
         sxyend  = varend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=var_ord(jsxy)
            indsxy= sxy_ord(jsxy)
C
            sxy_nam(indsxy)=var_nam(indsxy)
            sxy_dim(indsxy)=var_dim(indsxy)
            sxy_jpi(indsxy)=var_jpi(indsxy)
            sxy_jpj(indsxy)=var_jpj(indsxy)
            sxy_jpk(indsxy)=var_jpk(indsxy)
            sxy_jpt(indsxy)=var_jpt(indsxy)
            sxy_nbr(indsxy)=var_nbr(indsxy)
            sxy_ind(indsxy)=var_ind(indsxy)
C
         ENDDO
      CASE(2)
C --- dta
         jpsend  = jpyend
         sxyend  = dtaend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=dta_ord(jsxy)
            indsxy= sxy_ord(jsxy)
C
            sxy_nam(indsxy)=dta_nam(indsxy)
            sxy_dim(indsxy)=dta_dim(indsxy)
            sxy_jpi(indsxy)=dta_jpi(indsxy)
            sxy_jpj(indsxy)=dta_jpj(indsxy)
            sxy_jpk(indsxy)=dta_jpk(indsxy)
            sxy_jpt(indsxy)=dta_jpt(indsxy)
            sxy_nbr(indsxy)=dta_nbr(indsxy)
            sxy_ind(indsxy)=dta_ind(indsxy)
C
         ENDDO
      CASE(3)
         GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Set array size and allocate 2D array to read from file
      jpisize=1
      jpjsize=1
      jpksize=size(mask,3)
      DO jsxy = 1,sxyend
         indsxy=sxy_ord(jsxy)
         jpisize=MAX(sxy_jpi(indsxy),jpisize)
         jpjsize=MAX(sxy_jpj(indsxy),jpjsize)
      ENDDO
C
      allocate (ptabij(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabij(:,:) = FREAL8(0.0)
C
C -1.- Open file
C --------------
C Open user defined .rst file format
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C INSERT INSTRUCTIONS TO OPEN .rst FILES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Example: direct access binary file:
C     irecl =
C     CALL openfile(numfil,kfninsxy,clold,clunf,cldir,irecl)
C
C -2.- Read file header
C ---------------------
C Read header of user defined .rst file format
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C INSERT INSTRUCTIONS TO READ .rst FILE HEADER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Example: direct access binary file (1st record):
C     READ( UNIT=numfil ,REC=1, ERR=101, IOSTAT=iost ) ...
C
C Store information in memory (useful to write restart file)
C (Define variables in 'mod_cfgxyo.F' to store required information)
      IF (klectinfo) THEN
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C INSERT INSTRUCTIONS TO STORE .rst HEADER INFORMATIONS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      ENDIF
C
C -3.- Read 2D variable fields
C ----------------------------
C
      nsdeb = 1
      nstot = 0
      nsfin = 1
C
      DO jsxy = 1,sxyend
         indsxy=sxy_ord(jsxy)
C
C Read variable number 'indsxy' from file
C [Variable name: sxy_nam(indsxy)]
C and store it in output vector
         IF (nsdeb.NE.sxy_ind(indsxy)) GOTO 1000
         IF (sxy_nbr(indsxy).EQ.0) GOTO 1000
         nstot = sxy_nbr(indsxy)
         nsfin = nsdeb - 1 + nstot
         IF (nsfin.GT.size(kvects,1)) GOTO 1000
C
         js=nsdeb
         jt=1
         somsxynbr = 0
         DO jk=1,sxy_jpk(indsxy)
            sompartsxynbr=0
            IF (js.LT.jpssize) THEN
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C INSERT INSTRUCTIONS TO READ LAYER jk IN .rst file
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Example: direct access binary file:
C   nrec =
C   READ(UNIT=numfil,REC=nrec+jk-1,ERR=101,IOSTAT=iost) ptabij(:,:)
               CALL mk8vct(kvects(js:nsfin),ptabij(:,:),jk,jt,
     $              jsxy,sompartsxynbr,kflagxyo)
C Use 'mk8vct' to read real(kind=8) type
C and 'mk4vct' to read real(kind=4) type
            ENDIF
            js=js+sompartsxynbr
            somsxynbr= somsxynbr+sompartsxynbr
         END DO
C
C Control print
         IF (nprint.GE.2) THEN
            WRITE(numout,*) ' '
            SELECT CASE(kflagxyo)
            CASE(1)
C --- var
               WRITE(numout,*) '    ==> STORING variable ',
     $              var_nam(indsxy)(1:lenv(var_nam(indsxy))),
     $              ' in Vx vector'
            CASE(2)
C --- dta
               WRITE(numout,*) '    ==> STORING variable ',
     $              dta_nam(indsxy)(1:lenv(dta_nam(indsxy))),
     $              ' in Vy vector'
            CASE DEFAULT
               GOTO 1000
            END SELECT
            kform='(8x,a,i9)'
            WRITE(numout,kform) '- Initial index : ',nsdeb
            WRITE(numout,kform) '- Final index   : ',nsfin
            WRITE(numout,kform) '- Segment size  : ',nsfin-nsdeb+1
         ENDIF 
         IF (nstot.NE.somsxynbr) GOTO 1000
         nsdeb = nsfin+1
C
      ENDDO
C
      jpsend = nsfin
C
C Check size of vector object
C
      SELECT CASE(kflagxyo)
      CASE(1)
C --- var
         IF (jpsend.NE.jpxend) GOTO 102
      CASE(2)
C --- dta
         IF (jpsend.NE.jpyend) GOTO 103
      CASE(3)
         GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -3.- Close file
C ---------------
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C INSERT INSTRUCTIONS TO CLOSE .rst FILES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Example: direct access binary file:
C     CLOSE (UNIT=numfil)
C
C --- deallocation
      IF (allocated(ptabij)) deallocate(ptabij)
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'liorst','readrst')
 1001 CALL printerror2(0,1001,3,'liorst','readrst')
C
 101  WRITE (texterror,*) 'Error reading .rst file, iost=',iost
      CALL printerror2(0,101,3,'liorst','readrst',comment=texterror)
 102  WRITE (texterror,*) 'incoherence between jpxend =',jpsend,
     $     ' and jpxend from mask =',jpxend
      CALL printerror2(0,102,3,'liorst','readrst',comment=texterror)
 103  WRITE (texterror,*) 'incoherence between jpyend =',jpsend,
     $     ' and jpyend from mask =',jpyend
      CALL printerror2(0,103,3,'liorst','readrst',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writerst(kfnoutrst,kvectxin)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write Vx vector in user defined format
CCC  -------
CCC  Method : Loop on variable fields defined in SESAM configuration,
CCC  ------   get them from Vx vector, and write them in output file
CCC
CCC  Input :  kfnoutrst : filename
CCC  -----    kvectxin  : 1D vector object (Vx)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnoutrst
      BIGREAL, dimension(:), intent(in) :: kvectxin
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpisize,jpjsize,jpksize
      BIGREAL8, allocatable, dimension(:,:) :: ptabij
      BIGREAL8 :: spval
      INTEGER :: jvar,ji,jj,jk,jt,jx
      INTEGER :: ndeb, indvar, nrec, iblocrst
      INTEGER :: kvardim,somvarnbr,sompartvarnbr,flagxyo
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writevar/writerst'
         WRITE(numout,*) '    ==> WRITING file ',
     $            kfnoutrst(1:lenv(kfnoutrst))
      ENDIF
C
C Set array size and allocate 2D array to write in file
      jpksize=size(mask,3)
      jpisize=1
      jpjsize=1
      DO jvar = 1,varend
         indvar=var_ord(jvar)
         jpisize=MAX(var_jpi(indvar),jpisize)
         jpjsize=MAX(var_jpj(indvar),jpjsize)
      ENDDO
C
      allocate (ptabij(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabij(:,:) = FREAL8(0.0)
      spval = FREAL8(0.0)
C
C -1.- Open file
C --------------
C Open user defined .rst file format
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C INSERT INSTRUCTIONS TO OPEN .rst FILES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Example: direct access binary file:
C     irecl =
C     CALL openfile(numfil,kfnoutrst,clold,clunf,cldir,irecl)
C
C -2.- Write file header
C-----------------------
C Get file header information (from memory if previously read)
C (from variables defined in 'mod_cfgxyo.F', see readrst.F)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C INSERT INSTRUCTIONS TO GET .rst HEADER INFORMATIONS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Write header of user defined .rst file format
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C INSERT INSTRUCTIONS TO WRITE .rst FILE HEADER
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Example: direct access binary file (1st record):
C     WRITE( UNIT=numfil ,REC=1, ERR=101, IOSTAT=iost ) ...
C
C -3.- Write 2D variable fields
C -----------------------------
C
      flagxyo=1
      ndeb = 1
      DO jvar = 1,varend
         indvar=var_ord(jvar)
C
C Get variable 'indvar' from input vector and write it in file
C [Variable name: var_nam(indvar)]
         jt=1
         jx=ndeb
         somvarnbr=0
         DO jk=1,var_jpk(indvar)
            ptabij(:,:) = FREAL8(0.0)
            sompartvarnbr=0
            IF (jx.LT.size(kvectxin,1)) THEN
C Use 'unmk8vct' to write real(kind=8) type
C and 'unmk4vct' to write real(kind=4) type
               CALL unmk8vct(kvectxin(jx:),ptabij(:,:),jk,jt,
     $              jvar,sompartvarnbr,spval,flagxyo)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C INSERT INSTRUCTIONS TO WRITE LAYER jk IN .rst file
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Example: direct access binary file:
C   nrec =
C   WRITE(UNIT=numfil,REC=nrec+jk-1,ERR=101,IOSTAT=iost) ptabij(:,:)
            ENDIF
            jx=jx+sompartvarnbr
            somvarnbr= somvarnbr+sompartvarnbr
         ENDDO
C
         IF (var_nbr(indvar).NE.somvarnbr) GOTO 1000
         ndeb = ndeb + var_nbr(indvar)
C
      ENDDO
C
C -3.- Close file
C ---------------
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C INSERT INSTRUCTIONS TO CLOSE .rst FILES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Example: direct access binary file:
C     CLOSE (UNIT=numfil)
C
C --- deallocation
      IF (allocated(ptabij)) deallocate(ptabij)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liorst','writerst')
 1001 CALL printerror2(0,1001,3,'liorst','writerst')
C
 101  WRITE (texterror,*) 'Error writing .rst file, iost=',iost
      CALL printerror2(0,101,3,'liorst','writerst',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liorst
