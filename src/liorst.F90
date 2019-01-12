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
! ---                   LIORST.F90                                ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12  (C.E. Testut)                       ---
! --- modification : 99-05  (C.E. Testut)                       ---
! --- modification : 01-06  (C.E. Testut)                       ---
! --- modification : 03-04  (J.M. Brankart)                     ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  readrst  : Read Vx or Vy vector from user defined format
! --- SUBROUTINE  writerst : Write Vx vector in user defined format
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liorst
      use mod_main
      use utilvct
      use utilfiles
      IMPLICIT NONE
      PRIVATE

      PUBLIC readrst,writerst

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readrst(kfninsxy,kvects,klectinfo,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Read Vx or Vy vector from user defined format
!  -------
!  Method : Loop on variable fields defined in SESAM configuration,
!  ------   read them in input file, and store them in Vx or Vy vector
!
!  Input :  kfninsxy  : filename
!  -----    klectinfo : read or not header of input files
!           kflagxyo  : vector type (Vx/Vy/Vo)
!
!  Output : kvects : 1D vector object (Vx or Vy)
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : jpxend, jpyend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninsxy
      BIGREAL, dimension(:), intent(out) :: kvects
      LOGICAL, intent(in) :: klectinfo
      INTEGER, intent(in) :: kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: kform
      CHARACTER(len=varlg), dimension(1:nbvar) :: sxy_nam
      INTEGER :: jpsend,sxyend
      INTEGER, dimension(1:nbvar) :: sxy_ord,sxy_dim, &
     &     sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt,sxy_nbr,sxy_ind
      INTEGER :: jsxy,indsxy,somsxynbr,sompartsxynbr
      INTEGER :: nsdeb, nsfin, nstot
      INTEGER :: allocok,jpisize,jpjsize,jpksize,jpssize
      BIGREAL8, allocatable, dimension(:,:) :: ptabij
      INTEGER :: jinf,ji,jj,jk,jt,js
      INTEGER :: nrec,iblocrst
!----------------------------------------------------------------------
      jpssize=size(kvects,1)
!---------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.3) THEN
         WRITE(numout,*) '*** ROUTINE : ../readvar/readrst'
         WRITE(numout,*) '    ==> READING file ', &
     &            kfninsxy(1:lenv(kfninsxy))
      ENDIF
!
! Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         jpsend  = jpxend
         sxyend  = varend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=var_ord(jsxy)
            indsxy= sxy_ord(jsxy)
!
            sxy_nam(indsxy)=var_nam(indsxy)
            sxy_dim(indsxy)=var_dim(indsxy)
            sxy_jpi(indsxy)=var_jpi(indsxy)
            sxy_jpj(indsxy)=var_jpj(indsxy)
            sxy_jpk(indsxy)=var_jpk(indsxy)
            sxy_jpt(indsxy)=var_jpt(indsxy)
            sxy_nbr(indsxy)=var_nbr(indsxy)
            sxy_ind(indsxy)=var_ind(indsxy)
!
         ENDDO
      CASE(2)
! --- dta
         jpsend  = jpyend
         sxyend  = dtaend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=dta_ord(jsxy)
            indsxy= sxy_ord(jsxy)
!
            sxy_nam(indsxy)=dta_nam(indsxy)
            sxy_dim(indsxy)=dta_dim(indsxy)
            sxy_jpi(indsxy)=dta_jpi(indsxy)
            sxy_jpj(indsxy)=dta_jpj(indsxy)
            sxy_jpk(indsxy)=dta_jpk(indsxy)
            sxy_jpt(indsxy)=dta_jpt(indsxy)
            sxy_nbr(indsxy)=dta_nbr(indsxy)
            sxy_ind(indsxy)=dta_ind(indsxy)
!
         ENDDO
      CASE(3)
         GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Set array size and allocate 2D array to read from file
      jpisize=1
      jpjsize=1
      jpksize=size(mask,3)
      DO jsxy = 1,sxyend
         indsxy=sxy_ord(jsxy)
         jpisize=MAX(sxy_jpi(indsxy),jpisize)
         jpjsize=MAX(sxy_jpj(indsxy),jpjsize)
      ENDDO
!
      allocate (ptabij(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabij(:,:) = FREAL8(0.0)
!
! -1.- Open file
! --------------
! Open user defined .rst file format
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INSERT INSTRUCTIONS TO OPEN .rst FILES
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Example: direct access binary file:
!     irecl =
!     CALL openfile(numfil,kfninsxy,clold,clunf,cldir,irecl)
!
! -2.- Read file header
! ---------------------
! Read header of user defined .rst file format
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INSERT INSTRUCTIONS TO READ .rst FILE HEADER
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Example: direct access binary file (1st record):
!     READ( UNIT=numfil ,REC=1, ERR=101, IOSTAT=iost ) ...
!
! Store information in memory (useful to write restart file)
! (Define variables in 'mod_cfgxyo.F90' to store required information)
      IF (klectinfo) THEN
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INSERT INSTRUCTIONS TO STORE .rst HEADER INFORMATIONS
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      ENDIF
!
! -3.- Read 2D variable fields
! ----------------------------
!
      nsdeb = 1
      nstot = 0
      nsfin = 1
!
      DO jsxy = 1,sxyend
         indsxy=sxy_ord(jsxy)
!
! Read variable number 'indsxy' from file
! [Variable name: sxy_nam(indsxy)]
! and store it in output vector
         IF (nsdeb.NE.sxy_ind(indsxy)) GOTO 1000
         IF (sxy_nbr(indsxy).EQ.0) GOTO 1000
         nstot = sxy_nbr(indsxy)
         nsfin = nsdeb - 1 + nstot
         IF (nsfin.GT.size(kvects,1)) GOTO 1000
!
         js=nsdeb
         jt=1
         somsxynbr = 0
         DO jk=1,sxy_jpk(indsxy)
            sompartsxynbr=0
            IF (js.LT.jpssize) THEN
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INSERT INSTRUCTIONS TO READ LAYER jk IN .rst file
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Example: direct access binary file:
!   nrec =
!   READ(UNIT=numfil,REC=nrec+jk-1,ERR=101,IOSTAT=iost) ptabij(:,:)
               CALL mk8vct(kvects(js:nsfin),ptabij(:,:),jk,jt, &
     &              jsxy,sompartsxynbr,kflagxyo)
! Use 'mk8vct' to read real(kind=8) type
! and 'mk4vct' to read real(kind=4) type
            ENDIF
            js=js+sompartsxynbr
            somsxynbr= somsxynbr+sompartsxynbr
         END DO
!
! Control print
         IF (nprint.GE.2) THEN
            WRITE(numout,*) ' '
            SELECT CASE(kflagxyo)
            CASE(1)
! --- var
               WRITE(numout,*) '    ==> STORING variable ', &
     &              var_nam(indsxy)(1:lenv(var_nam(indsxy))), &
     &              ' in Vx vector'
            CASE(2)
! --- dta
               WRITE(numout,*) '    ==> STORING variable ', &
     &              dta_nam(indsxy)(1:lenv(dta_nam(indsxy))), &
     &              ' in Vy vector'
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
!
      ENDDO
!
      jpsend = nsfin
!
! Check size of vector object
!
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         IF (jpsend.NE.jpxend) GOTO 102
      CASE(2)
! --- dta
         IF (jpsend.NE.jpyend) GOTO 103
      CASE(3)
         GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -3.- Close file
! ---------------
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INSERT INSTRUCTIONS TO CLOSE .rst FILES
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Example: direct access binary file:
!     CLOSE (UNIT=numfil)
!
! --- deallocation
      IF (allocated(ptabij)) deallocate(ptabij)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'liorst','readrst')
 1001 CALL printerror2(0,1001,3,'liorst','readrst')
!
 101  WRITE (texterror,*) 'Error reading .rst file, iost=',iost
      CALL printerror2(0,101,3,'liorst','readrst',comment=texterror)
 102  WRITE (texterror,*) 'incoherence between jpxend =',jpsend, &
     &     ' and jpxend from mask =',jpxend
      CALL printerror2(0,102,3,'liorst','readrst',comment=texterror)
 103  WRITE (texterror,*) 'incoherence between jpyend =',jpsend, &
     &     ' and jpyend from mask =',jpyend
      CALL printerror2(0,103,3,'liorst','readrst',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writerst(kfnoutrst,kvectxin)
!---------------------------------------------------------------------
!
!  Purpose : Write Vx vector in user defined format
!  -------
!  Method : Loop on variable fields defined in SESAM configuration,
!  ------   get them from Vx vector, and write them in output file
!
!  Input :  kfnoutrst : filename
!  -----    kvectxin  : 1D vector object (Vx)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutrst
      BIGREAL, dimension(:), intent(in) :: kvectxin
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jpisize,jpjsize,jpksize
      BIGREAL8, allocatable, dimension(:,:) :: ptabij
      BIGREAL8 :: spval
      INTEGER :: jvar,ji,jj,jk,jt,jx
      INTEGER :: ndeb, indvar, nrec, iblocrst
      INTEGER :: kvardim,somvarnbr,sompartvarnbr,flagxyo
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writevar/writerst'
         WRITE(numout,*) '    ==> WRITING file ', &
     &            kfnoutrst(1:lenv(kfnoutrst))
      ENDIF
!
! Set array size and allocate 2D array to write in file
      jpksize=size(mask,3)
      jpisize=1
      jpjsize=1
      DO jvar = 1,varend
         indvar=var_ord(jvar)
         jpisize=MAX(var_jpi(indvar),jpisize)
         jpjsize=MAX(var_jpj(indvar),jpjsize)
      ENDDO
!
      allocate (ptabij(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabij(:,:) = FREAL8(0.0)
      spval = FREAL8(0.0)
!
! -1.- Open file
! --------------
! Open user defined .rst file format
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INSERT INSTRUCTIONS TO OPEN .rst FILES
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Example: direct access binary file:
!     irecl =
!     CALL openfile(numfil,kfnoutrst,clold,clunf,cldir,irecl)
!
! -2.- Write file header
!-----------------------
! Get file header information (from memory if previously read)
! (from variables defined in 'mod_cfgxyo.F90', see readrst.F)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INSERT INSTRUCTIONS TO GET .rst HEADER INFORMATIONS
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Write header of user defined .rst file format
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INSERT INSTRUCTIONS TO WRITE .rst FILE HEADER
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Example: direct access binary file (1st record):
!     WRITE( UNIT=numfil ,REC=1, ERR=101, IOSTAT=iost ) ...
!
! -3.- Write 2D variable fields
! -----------------------------
!
      flagxyo=1
      ndeb = 1
      DO jvar = 1,varend
         indvar=var_ord(jvar)
!
! Get variable 'indvar' from input vector and write it in file
! [Variable name: var_nam(indvar)]
         jt=1
         jx=ndeb
         somvarnbr=0
         DO jk=1,var_jpk(indvar)
            ptabij(:,:) = FREAL8(0.0)
            sompartvarnbr=0
            IF (jx.LT.size(kvectxin,1)) THEN
! Use 'unmk8vct' to write real(kind=8) type
! and 'unmk4vct' to write real(kind=4) type
               CALL unmk8vct(kvectxin(jx:),ptabij(:,:),jk,jt, &
     &              jvar,sompartvarnbr,spval,flagxyo)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INSERT INSTRUCTIONS TO WRITE LAYER jk IN .rst file
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Example: direct access binary file:
!   nrec =
!   WRITE(UNIT=numfil,REC=nrec+jk-1,ERR=101,IOSTAT=iost) ptabij(:,:)
            ENDIF
            jx=jx+sompartvarnbr
            somvarnbr= somvarnbr+sompartvarnbr
         ENDDO
!
         IF (var_nbr(indvar).NE.somvarnbr) GOTO 1000
         ndeb = ndeb + var_nbr(indvar)
!
      ENDDO
!
! -3.- Close file
! ---------------
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INSERT INSTRUCTIONS TO CLOSE .rst FILES
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Example: direct access binary file:
!     CLOSE (UNIT=numfil)
!
! --- deallocation
      IF (allocated(ptabij)) deallocate(ptabij)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'liorst','writerst')
 1001 CALL printerror2(0,1001,3,'liorst','writerst')
!
 101  WRITE (texterror,*) 'Error writing .rst file, iost=',iost
      CALL printerror2(0,101,3,'liorst','writerst',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liorst
