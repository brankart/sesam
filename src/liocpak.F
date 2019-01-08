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
C ---                  LIOCPAK.F                                ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-10  (C.E. Testut)                       ---
C --- modification : 01-06  (C.E. Testut)                       ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE  readcpak      : Read Vx or Vy vector from 'cpak' file
C --- SUBROUTINE  readpartcpak  : Read Vx segment from 'cpak' file
C --- SUBROUTINE  writecpak     : Write Vx vector in 'cpak' file
C --- SUBROUTINE  writepartcpak : Write Vx segment in 'cpak' file
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liocpak
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC readcpak,readpartcpak,writecpak,writepartcpak

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcpak(kfninsxy,kvects,klectinfo,kflagxyo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read Vx or Vy vector from 'cpak' file
CCC  -------
CCC  Method : Check cpak file compatibility and
CCC  ------   read required information from file
CCC
CCC  Input :  kfninsxy  : filename
CCC  -----    klectinfo : read or not header of 'cpak' file
CCC           kflagxyo  : vector type (Vx,Vy)
CCC  Output : kvects    : 1D vector (Vx or Vy)
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : jpxend, jpyend
      use utilcdfpak
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
      CHARACTER(len=word80) :: cdrec1
      INTEGER :: jpxend1, varend1, varlg1
      INTEGER, dimension(1:nbvar) :: var_ord1,var_dim1,var_nbr1
      CHARACTER(len=varlg), dimension(1:nbvar) :: var_nam1
      BIGREAL4, dimension(1:nbvar) :: var_moy1,var_ect1
      BIGREAL :: sxy_moy,sxy_ect
      INTEGER :: allocok,jpssize
      BIGREAL4, allocatable, dimension(:) :: ptabx
      LOGICAL :: incompatible,sameend,sameord,affect,transfert
      INTEGER :: jvar,jvar1,indvar,jdta,inddta,jrec
      INTEGER :: jindy,jindybeg,jindyend,jindx,jindxbeg,jindxend
      INTEGER :: js,jsdeb,jsfin,lgcdrec,lgvarnam
      LOGICAL :: ltest
CC----------------------------------------------------------------------
C
C Check input arguments
      jpssize = size(kvects,1)
      IF ((kflagxyo.EQ.1).AND.(jpssize.NE.jpxend)) GOTO 104
      IF ((kflagxyo.EQ.2).AND.(jpssize.NE.jpyend)) GOTO 104
      IF (kflagxyo.EQ.3) GOTO 106
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readvar/readcpak'
         WRITE(numout,*) '    ==> READING file ',kfninsxy(1:lenv(kfninsxy))
      ENDIF
C
C -1.- Read 'cpak' file dimensions
C --------------------------------
C
      CALL cdfrdimpak(kfninsxy,jpxend1,varend1,varlg1,cdrec1)
C
C -2.- Read 'cpak' file header
C ----------------------------
C
      CALL cdfrhdrpak(kfninsxy,var_nam1,var_dim1,
     $     var_nbr1,var_moy1,var_ect1)
C
C -3.- Check 'cpak' file compatibilty
C -----------------------------------
C
C Compute order of variables fields in 'cpak' file
      DO jvar1=1,varend1
         affect = .FALSE.
         LOOP1 : DO indvar=1,nbvar 
            IF (var_nam1(jvar1).EQ.var_nam(indvar)) THEN
               var_ord1(jvar1)=indvar
               affect = .TRUE.
               EXIT LOOP1
            ENDIF
         ENDDO LOOP1  
         IF (.NOT.(affect)) GOTO 103
      ENDDO
      var_ord1((varend1+1):nbvar) = 0
C
C Check if 'cpak' file configuration is identical to SESAM configuration
      sameend = .TRUE.
      sameend = ((varend1.EQ.varend).AND.sameend)
      sameend = ((jpxend1.EQ.jpxend).AND.sameend)
      sameend = ((varlg1.LE.varlg).AND.sameend)
      sameord = .TRUE.
      DO jvar1=1,varend1
         sameord = ((var_ord1(jvar1).EQ.var_ord(jvar1)).AND.sameord)
      ENDDO
C
      incompatible = (.NOT.(sameord.AND.sameend))
      DO jvar=1,varend
         indvar=var_ord(jvar)
         incompatible = ((var_nam(indvar).NE.var_nam1(jvar))
     $        .OR.incompatible)
         incompatible = ((var_dim(indvar).NE.var_dim1(jvar))
     $        .OR.incompatible)
         incompatible = ((var_nbr(indvar).NE.var_nbr1(jvar))
     $        .OR.incompatible)
      ENDDO
C
      IF (incompatible) GOTO 102
C
C -4.- Take decision for centering/reducing input variables
C ---------------------------------------------------------
C
      transfert=.FALSE.
      DO jvar=1,varend
         indvar=var_ord(jvar)
         transfert = (
     $       ((FREAL4(var_moy(indvar)).NE.var_moy1(jvar)).AND.(lmoyect))
     $       .OR. ((var_moy1(jvar).NE.(FREAL4(0.0))).AND.(.NOT.lmoyect))
     $       .OR. transfert )
         transfert = (
     $       ((FREAL4(var_ect(indvar)).NE.var_ect1(jvar)).AND.(lmoyect))
     $       .OR. ((var_ect1(jvar).NE.(FREAL4(1.0))).AND.(.NOT.lmoyect))
     $       .OR. transfert )
      ENDDO
C
C -5.- Read Vx or Vy vector from cpak file
C ----------------------------------------
C
      SELECT CASE (kflagxyo)
C
C Read Vx vector
      CASE(1)
C
C Allocate Vx vector working array (kr4 real kind)
         allocate (ptabx(1:jpxend1), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptabx(:) = FREAL4(0.0)
C
C Read Vx vector
         CALL cdfrpak(kfninsxy,ptabx,1,jpxend1)
C
C Center/reduce input data
         IF (transfert) THEN
            DO jvar=1,varend
               indvar=var_ord(jvar)
               jsdeb=var_ind(indvar)
               jsfin=var_ind(indvar)+var_nbr(indvar)-1
               IF (lmoyect) THEN
                  sxy_ect=FREAL(var_ect1(jvar))/var_ect(indvar)
                  sxy_moy=(FREAL(var_moy1(jvar))-var_moy(indvar))
     $                 /var_ect(indvar)
               ELSE
                  sxy_ect=FREAL(var_ect1(jvar))
                  sxy_moy=FREAL(var_moy1(jvar))
               ENDIF
               kvects(jsdeb:jsfin)=FREAL(ptabx(jsdeb:jsfin))
     $              *sxy_ect+sxy_moy
           ENDDO
         ELSE
            kvects(:)=FREAL(ptabx(:))
         ENDIF
C
C Read Vy vector
      CASE(2)
C
C Allocate Vx vector working array (kr4 real kind)
         allocate ( ptabx(1:jpyend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptabx(:) = FREAL4(0.0)
C
C Fill progressively Vy vector
         jindybeg=1
         jindyend=jpyend
         IF (tabindxtoy(jindyend).GT.jpxend1) GOTO 1000
C
         jindy=1
         LOOP2 : DO jindxbeg=tabindxtoy(jindybeg),
     $                       tabindxtoy(jindyend),jpyend
            jindxend=MIN(tabindxtoy(jindyend),jindxbeg-1+jpyend)
            IF ( (tabindxtoy(jindy).LE.jindxend)
     $           .AND.(tabindxtoy(jindy).GE.jindxbeg) ) THEN
C
               CALL cdfrpak(kfninsxy,ptabx(1),jindxbeg,jindxend)
C
               jindx=tabindxtoy(jindy)-jindxbeg+1
               ltest=((tabindxtoy(jindy)-jindxbeg+1).LE.size(ptabx,1))
               DO WHILE ( (jindy.LE.jindyend)
     $              .AND.(ltest) )
                  jindx=tabindxtoy(jindy)-jindxbeg+1
                  kvects(jindy)=FREAL(ptabx(jindx))
                  jindy=jindy+1
                  IF (jindy.LE.jindyend) ltest =
     $                 ((tabindxtoy(jindy)-jindxbeg+1).LE.size(ptabx,1))
               ENDDO
C
            ENDIF
            IF (jindy.GT.jindyend) EXIT LOOP2
         ENDDO LOOP2
C
C Center/reduce input data
         IF (transfert) THEN
            DO jdta=1,dtaend
               inddta=dta_ord(jdta)
               jvar=1
               indvar = var_ord(jvar)
               DO WHILE ((inddta.NE.indvar).AND.(jvar.LT.varend))
                  jvar=jvar+1
                  indvar = var_ord(jvar)
               ENDDO  
               IF (inddta.NE.indvar) GOTO 1000
               jsdeb=dta_ind(inddta)
               jsfin=dta_ind(inddta)+dta_nbr(inddta)-1
               IF (lmoyect) THEN
                  sxy_ect=FREAL(var_ect1(jvar))/dta_ect(inddta)
                  sxy_moy=(FREAL(var_moy1(jvar))-dta_moy(inddta))
     $                 /dta_ect(inddta)
               ELSE
                  sxy_ect=FREAL(var_ect1(jvar))
                  sxy_moy=FREAL(var_moy1(jvar))
               ENDIF
               kvects(jsdeb:jsfin)=kvects(jsdeb:jsfin)*sxy_ect+sxy_moy
            ENDDO
         ENDIF
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C --- deallocation
      IF (allocated(ptabx)) deallocate(ptabx)
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'liocpak','readcpak')
 1001 CALL printerror2(0,1001,3,'liocpak','readcpak')
C
 102  WRITE (texterror,*) 'Parameters in .cpak file header',
     $     ' incompatible with SESAM configuration'
      CALL printerror2(0,102,3,'liocpak','readcpak',comment=texterror)
 103  WRITE (texterror,*) 'Variable names in .cpak file',
     $     ' incompatible with SESAM configuration'
      CALL printerror2(0,103,3,'liocpak','readcpak',comment=texterror)
 104  WRITE (texterror,*) 'Invalid input array size'
      CALL printerror2(0,104,1,'liocpak','readcpak',comment=texterror)
 106  WRITE (texterror,*) 'Cannot read Vo vector from cpak file'
      CALL printerror2(0,106,1,'liocpak','readcpak',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readpartcpak(kfninsxy,kvects,kjns,klectinfo,kflagxyo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read Vx segment from 'cpak' file
CCC  -------
CCC  Method : Check cpak file compatibility and
CCC  ------   read required segment from file
CCC
CCC  Input :  kfninsxy  : filename
CCC  -----    klectinfo : read or not header of 'cpak' file
CCC           kflagxyo  : vector type (Vx,Vy)
CCC           kjns      : index of Vx segment to read
CCC  Output : kvects    : 1D vector (Vx segment)
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : jpxend, jpx, arraynx_jindxbeg,
     $     arraynx_jpindxend
      use utilcdfpak
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfninsxy
      BIGREAL, dimension(:), intent(out) :: kvects
      INTEGER, intent(in) :: kjns
      LOGICAL, intent(in) :: klectinfo
      INTEGER, intent(in) :: kflagxyo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=word80) :: cdrec1
      INTEGER :: jpxend1
      INTEGER :: varend1
      INTEGER :: varlg1
      INTEGER, dimension(1:nbvar) :: var_ord1,var_dim1,var_nbr1
      CHARACTER(len=varlg), dimension(1:nbvar) :: var_nam1
      BIGREAL4, dimension(1:nbvar) :: var_moy1,var_ect1
      BIGREAL :: sxy_moy,sxy_ect
      INTEGER :: allocok,jpssize
      BIGREAL4, allocatable, dimension(:) :: ptabx
      LOGICAL :: incompatible,sameend,sameord,affect,transfert
      INTEGER :: jvar,jvar1,indvar,jdta,inddta,jrec
      INTEGER :: jindy,jindybeg,jindyend,jindx,jindxbeg,jindxend
      INTEGER :: js,jsdeb,jsfin,lgcdrec,lgvarnam
CC----------------------------------------------------------------------
C
C Check input arguments
      jpssize = size(kvects,1)
      IF ((kflagxyo.EQ.1).AND.(jpssize.NE.jpx)) GOTO 104
      IF (kflagxyo.EQ.2) GOTO 105
      IF (kflagxyo.EQ.3) GOTO 106
C
C Control print
      IF (kjns.EQ.1) THEN
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readvar/readpartcpak'
         WRITE(numout,*) '    ==> READING file ',kfninsxy(1:lenv(kfninsxy))
      ENDIF
      ENDIF
C
C -1.- Read 'cpak' file header
C ----------------------------
C
      CALL cdfrhdrpak(kfninsxy,var_nam1,var_dim1,
     $        var_nbr1,var_moy1,var_ect1)
C
      IF (kjns.EQ.1) THEN
C
C -2.- Read 'cpak' file dimensions
C --------------------------------
C
         CALL cdfrdimpak(kfninsxy,jpxend1,
     $        varend1,varlg1,cdrec1)
C
C -3.- Check 'cpak' file compatibilty
C -----------------------------------
C
C Compute order of variables fields in 'cpak' file
         DO jvar1=1,varend1
            affect = .FALSE.
            LOOP1 : DO indvar=1,nbvar 
               IF (var_nam1(jvar1).EQ.var_nam(indvar)) THEN
                  var_ord1(jvar1)=indvar
                  affect = .TRUE.
                  EXIT LOOP1
               ENDIF
            ENDDO LOOP1  
            IF (.NOT.(affect)) GOTO 103
         ENDDO
         var_ord1((varend1+1):nbvar) = 0
C
C Check if 'cpak' file configuration is identical to SESAM configuration
         sameend = .TRUE.
         sameend = ((varend1.EQ.varend).AND.sameend)
         sameend = ((jpxend1.EQ.jpxend).AND.sameend)
         sameend = ((varlg1.LE.varlg).AND.sameend)
         sameord = .TRUE.
         DO jvar1=1,varend1
            sameord = ((var_ord1(jvar1).EQ.var_ord(jvar1)).AND.sameord)
         ENDDO
C     
         incompatible = (.NOT.(sameord.AND.sameend))
         DO jvar=1,varend
            indvar=var_ord(jvar)
            incompatible = ((var_nam(indvar).NE.var_nam1(jvar))
     $           .OR.incompatible)
            incompatible = ((var_dim(indvar).NE.var_dim1(jvar))
     $           .OR.incompatible)
            incompatible = ((var_nbr(indvar).NE.var_nbr1(jvar))
     $           .OR.incompatible)
         ENDDO
C
         IF (incompatible) GOTO 102
      ENDIF
C
C -4.- Take decision for centering/reducing input variables
C ---------------------------------------------------------
C
      transfert=.FALSE.
      DO jvar=1,varend
         indvar=var_ord(jvar)
         transfert = (
     $       ((FREAL4(var_moy(indvar)).NE.var_moy1(jvar)).AND.(lmoyect))
     $       .OR. ((var_moy1(jvar).NE.(FREAL4(0.0))).AND.(.NOT.lmoyect))
     $       .OR. transfert )
         transfert = (
     $       ((FREAL4(var_ect(indvar)).NE.var_ect1(jvar)).AND.(lmoyect))
     $       .OR. ((var_ect1(jvar).NE.(FREAL4(1.0))).AND.(.NOT.lmoyect))
     $       .OR. transfert )
      ENDDO
C
C -5.- Read Vx segment from cpak file
C -----------------------------------
C
C Allocate Vx vector working array (kr4 real kind)
      allocate ( ptabx(1:jpx), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabx(:) = FREAL4(0.0)
C
C Read Vx segment
      jindxbeg=arraynx_jindxbeg(kjns)
      jindxend=arraynx_jindxbeg(kjns)-1+arraynx_jpindxend(kjns)
      CALL cdfrpak(kfninsxy,ptabx,jindxbeg,jindxend)
C
C Center/reduce input data
      IF (transfert) THEN
         DO jvar=1,varend
            indvar=var_ord(jvar)
            jsdeb=var_ind(indvar)-jindxbeg+1
            jsfin=var_ind(indvar)+var_nbr(indvar)-1-jindxbeg+1
            IF (lmoyect) THEN
               sxy_ect=FREAL(var_ect1(jvar))/var_ect(indvar)
               sxy_moy=(FREAL(var_moy1(jvar))-var_moy(indvar))
     $              /var_ect(indvar)
            ELSE
               sxy_ect=FREAL(var_ect1(jvar))
               sxy_moy=FREAL(var_moy1(jvar))
            ENDIF
            IF ((jsfin.GE.1).AND.(jsdeb.LE.jpx)) THEN
               kvects(MAX(jsdeb,1):MIN(jsfin,jpx))=
     $              FREAL(ptabx(MAX(jsdeb,1):MIN(jsfin,jpx)))
     $              *sxy_ect+sxy_moy
            ENDIF
         ENDDO
      ELSE
         kvects(:)=FREAL(ptabx(:))
      ENDIF
C
C --- deallocation
      IF (allocated(ptabx)) deallocate(ptabx)
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'liocpak','readpartcpak')
 1001 CALL printerror2(0,1001,3,'liocpak','readpartcpak')
C
 102  WRITE (texterror,*) 'Parameters in .cpak file header',
     $     ' incompatible with SESAM configuration'
      CALL printerror2(0,102,3,'liocpak','readpartcpak',comment=texterror)
 103  WRITE (texterror,*) 'Variable names in .cpak file',
     $     ' incompatible with SESAM configuration'
      CALL printerror2(0,103,3,'liocpak','readpartcpak',comment=texterror)
 104  WRITE (texterror,*) 'Invalid input array size'
      CALL printerror2(0,104,1,'liocpak','readpartcpak',comment=texterror)
 105  WRITE (texterror,*) 'Cannot read segment of Vy vector from cpak file'
      CALL printerror2(0,105,1,'liocpak','readpartcpak',comment=texterror)
 106  WRITE (texterror,*) 'Cannot read Vo vector from cpak file'
      CALL printerror2(0,106,1,'liocpak','readpartcpak',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writecpak(kfnoutcpak,kvectxin)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write Vx vector in 'cpak' file
CCC  -------
CCC  Method : Check if cpak file already exists
CCC  ------   Create it if necessary
CCC           Write or overwrite Vx vector in file
CCC
CCC  Input : kfnoutcpak : filename
CCC  -----   kvectxin   : 1D vector (Vx)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend
      use utilcdfpak
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnoutcpak
      BIGREAL, dimension(:), intent(in) :: kvectxin
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=word80) :: cdrec1,cdrec2,text
      INTEGER :: jpxend1,varend1,varlg1
      INTEGER :: jpxend2,varend2,varlg2
      INTEGER, dimension(1:nbvar) :: var_dim1, var_nbr1
      CHARACTER(len=varlg), dimension(1:nbvar) :: var_nam1
      BIGREAL4, dimension(1:nbvar) :: var_moy1,var_ect1
      INTEGER allocok,jpxsize
#if ! defined ALLREAL8
      BIGREAL4, allocatable, dimension(:) :: ptabx
#endif
      INTEGER :: jx,jvar,jvar1,indvar,indvar1
      LOGICAL :: existence
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writevar/writecpak'
         WRITE(numout,*) '    ==> WRITING file ',kfnoutcpak(1:lenv(kfnoutcpak))
      ENDIF
C
C Check vector size
      jpxsize = size(kvectxin,1)
      IF (jpxsize.NE.jpxend) GOTO 1000
C
C -1.- Set cpak header variables
C ------------------------------
C
      jpxend1  = jpxsize
      varend1  = varend
      varlg1   = varlg
C
      DO jvar1 = 1,varend
         indvar=var_ord(jvar1)
         var_nam1(jvar1) = var_nam(indvar)
         var_dim1(jvar1) = var_dim(indvar)
         var_nbr1(jvar1) = var_nbr(indvar)
      ENDDO
C
      IF (lmoyect) THEN
         DO jvar1 = 1,varend
            indvar=var_ord(jvar1)
            var_moy1(jvar1) = FREAL4(var_moy(indvar))
            var_ect1(jvar1) = FREAL4(var_ect(indvar))
         ENDDO
      ELSE
         DO jvar1 = 1,varend
            indvar=var_ord(jvar1)
            var_moy1(jvar1) = FREAL4(0.0)
            var_ect1(jvar1) = FREAL4(1.0)
         ENDDO
      ENDIF
C
      WRITE (cdrec1,'(a)') 'Vx =>'
      DO jvar1=1,varend1
         text=cdrec1
         WRITE (cdrec1,'(a,a)') text(1:lenv(text)),
     $        var_nam1(jvar1)(1:lenv(var_nam1(jvar1)))  
      ENDDO
C
C -2.- Check or write 'cpak' file dimensions
C ------------------------------------------
C
C Check if 'cpak' file already exists
      INQUIRE (FILE=kfnoutcpak,EXIST=existence)
C
      IF (existence) THEN
         CALL cdfrdimpak(kfnoutcpak,jpxend2,varend2,varlg2,cdrec2)
         IF (jpxend2.NE.jpxend1) GOTO 103
         IF (varend2.NE.varend1) GOTO 103
         IF (varlg2.NE.varlg1) GOTO 103
      ELSE
         CALL cdfwdimpak(kfnoutcpak,jpxend1,varend1,varlg1,cdrec1)
      ENDIF
C
C -3.- Write 'cpak' file header
C ------------------------------
C
      CALL cdfwhdrpak(kfnoutcpak,var_nam1,var_dim1,var_nbr1,var_moy1,var_ect1)
C
C -4.- Write Vx vector in cpak file
C ---------------------------------
C
#if ! defined ALLREAL8
      allocate (ptabx(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabx(:) = FREAL4(kvectxin(:))
      CALL cdfwpak(kfnoutcpak,ptabx(1),1,jpxend1)
      IF (allocated(ptabx)) deallocate (ptabx)
#else
      CALL cdfwpak(kfnoutcpak,kvectxin(1),1,jpxend1)
#endif
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'liocpak','writecpak')
 1001 CALL printerror2(0,1001,3,'liocpak','writecpak')
C
 103  WRITE (texterror,*) 'inconsistent dimensions in output cpak file:',
     $                        kfnoutcpak(1:lenv(kfnoutcpak))
      CALL printerror2(0,103,3,'liocpak','writecpak',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writepartcpak(kfnoutcpak,kvectxin,kjnx)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write Vx segment in 'cpak' file
CCC  -------
CCC  Method :
CCC  ------
CCC  Input : kfnoutcpak : filename
CCC  -----   kvectxin   : 1D vector (Vx segment)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend, jpx, arraynx_jindxbeg,
     $     arraynx_jpindxend
      use utilcdfpak
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnoutcpak
      BIGREAL, dimension(:), intent(in) :: kvectxin
      INTEGER , intent(in) :: kjnx
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=word80) :: cdrec1,cdrec2,text
      INTEGER :: jpxend1,varend1,varlg1
      INTEGER :: jpxend2,varend2,varlg2
      INTEGER, dimension(1:nbvar) :: var_dim1, var_nbr1
      CHARACTER(len=varlg), dimension(1:nbvar) :: var_nam1
      BIGREAL4, dimension(1:nbvar) :: var_moy1,var_ect1
      INTEGER allocok,jpxsize
#if ! defined ALLREAL8
      BIGREAL4, allocatable, dimension(:) :: ptabx
#endif
      INTEGER :: jindxbeg,jindxend
      INTEGER :: jx,jvar,jvar1,indvar,indvar1
      LOGICAL :: existence
CC----------------------------------------------------------------------
C
C Control print
      IF (kjnx.EQ.1) THEN
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writevar/writecpak'
         WRITE(numout,*) '    ==> WRITING file ',kfnoutcpak(1:lenv(kfnoutcpak))
      ENDIF
      ENDIF
C
C Check vector size
      jpxsize = size(kvectxin,1)
      IF (jpxsize.NE.jpx) GOTO 1000
C
C Check and write header of cpak file only for 1st segment
      IF (kjnx.EQ.1) THEN
C
C -1.- Set cpak header variables
C ------------------------------
C
         jpxend1  = jpxend
         varend1  = varend
         varlg1   = varlg
         DO jvar1 = 1,varend
            indvar=var_ord(jvar1)
            var_nam1(jvar1) = var_nam(indvar)
            var_dim1(jvar1) = var_dim(indvar)
            var_nbr1(jvar1) = var_nbr(indvar)
         ENDDO
C
         IF (lmoyect) THEN
            DO jvar1 = 1,varend
               indvar=var_ord(jvar1)
               var_moy1(jvar1) = FREAL4(var_moy(indvar))
               var_ect1(jvar1) = FREAL4(var_ect(indvar))
            ENDDO
         ELSE
            DO jvar1 = 1,varend
               indvar=var_ord(jvar1)
               var_moy1(jvar1) = FREAL4(0.0)
               var_ect1(jvar1) = FREAL4(1.0)
            ENDDO
         ENDIF
C
         WRITE (cdrec1,'(a)') 'Vx =>'
         DO jvar1=1,varend1
            text=cdrec1
            WRITE (cdrec1,'(a,a)') text(1:lenv(text)),
     $           var_nam1(jvar1)(1:lenv(var_nam1(jvar1)))  
         ENDDO
C
C -2.- Check or write 'cpak' file dimensions
C ------------------------------------------
C
         INQUIRE (FILE=kfnoutcpak,EXIST=existence)
C
         IF (existence) THEN
            CALL cdfrdimpak(kfnoutcpak,jpxend2,varend2,varlg2,cdrec2)
            IF (jpxend2.NE.jpxend1) GOTO 103
            IF (varend2.NE.varend1) GOTO 103
            IF (varlg2.NE.varlg1) GOTO 103
         ELSE
            CALL cdfwdimpak(kfnoutcpak,jpxend1,varend1,varlg1,cdrec1)
         ENDIF
C
C -3.- Write 'cpak' file header
C -----------------------------
C
         CALL cdfwhdrpak(kfnoutcpak,var_nam1,var_dim1,var_nbr1,var_moy1,var_ect1)
C
      ENDIF
C
C -4.- Write Vx segment in cpak file
C ----------------------------------
C
      jindxbeg=arraynx_jindxbeg(kjnx)
      jindxend=arraynx_jindxbeg(kjnx)-1+arraynx_jpindxend(kjnx)
C
#if ! defined ALLREAL8
      allocate ( ptabx(1:jpxsize), stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      ptabx(:) = FREAL4(kvectxin(:))
      CALL cdfwpak(kfnoutcpak,ptabx(1),jindxbeg,jindxend)
      IF (allocated(ptabx)) deallocate (ptabx)
#else
      CALL cdfwpak(kfnoutcpak,kvectxin(1),jindxbeg,jindxend)
#endif
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'liocpak','writepartcpak')
 1001 CALL printerror2(0,1001,3,'liocpak','writeparcpak')
C
 103  WRITE (texterror,*) 'inconsistent dimensions in output cpak file:',
     $                        kfnoutcpak(1:lenv(kfnoutcpak))
      CALL printerror2(0,103,3,'liocpak','writeparcpak',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liocpak

