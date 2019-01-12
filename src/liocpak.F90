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
! ---                  LIOCPAK.F90                                ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-10  (C.E. Testut)                       ---
! --- modification : 01-06  (C.E. Testut)                       ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  readcpak      : Read Vx or Vy vector from 'cpak' file
! --- SUBROUTINE  readpartcpak  : Read Vx segment from 'cpak' file
! --- SUBROUTINE  writecpak     : Write Vx vector in 'cpak' file
! --- SUBROUTINE  writepartcpak : Write Vx segment in 'cpak' file
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liocpak
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC readcpak,readpartcpak,writecpak,writepartcpak

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcpak(kfninsxy,kvects,klectinfo,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Read Vx or Vy vector from 'cpak' file
!  -------
!  Method : Check cpak file compatibility and
!  ------   read required information from file
!
!  Input :  kfninsxy  : filename
!  -----    klectinfo : read or not header of 'cpak' file
!           kflagxyo  : vector type (Vx,Vy)
!  Output : kvects    : 1D vector (Vx or Vy)
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : jpxend, jpyend
      use utilcdfpak
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
!----------------------------------------------------------------------
!
! Check input arguments
      jpssize = size(kvects,1)
      IF ((kflagxyo.EQ.1).AND.(jpssize.NE.jpxend)) GOTO 104
      IF ((kflagxyo.EQ.2).AND.(jpssize.NE.jpyend)) GOTO 104
      IF (kflagxyo.EQ.3) GOTO 106
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readvar/readcpak'
         WRITE(numout,*) '    ==> READING file ',kfninsxy(1:lenv(kfninsxy))
      ENDIF
!
! -1.- Read 'cpak' file dimensions
! --------------------------------
!
      CALL cdfrdimpak(kfninsxy,jpxend1,varend1,varlg1,cdrec1)
!
! -2.- Read 'cpak' file header
! ----------------------------
!
      CALL cdfrhdrpak(kfninsxy,var_nam1,var_dim1, &
     &     var_nbr1,var_moy1,var_ect1)
!
! -3.- Check 'cpak' file compatibilty
! -----------------------------------
!
! Compute order of variables fields in 'cpak' file
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
!
! Check if 'cpak' file configuration is identical to SESAM configuration
      sameend = .TRUE.
      sameend = ((varend1.EQ.varend).AND.sameend)
      sameend = ((jpxend1.EQ.jpxend).AND.sameend)
      sameend = ((varlg1.LE.varlg).AND.sameend)
      sameord = .TRUE.
      DO jvar1=1,varend1
         sameord = ((var_ord1(jvar1).EQ.var_ord(jvar1)).AND.sameord)
      ENDDO
!
      incompatible = (.NOT.(sameord.AND.sameend))
      DO jvar=1,varend
         indvar=var_ord(jvar)
         incompatible = ((var_nam(indvar).NE.var_nam1(jvar)) &
     &        .OR.incompatible)
         incompatible = ((var_dim(indvar).NE.var_dim1(jvar)) &
     &        .OR.incompatible)
         incompatible = ((var_nbr(indvar).NE.var_nbr1(jvar)) &
     &        .OR.incompatible)
      ENDDO
!
      IF (incompatible) GOTO 102
!
! -4.- Take decision for centering/reducing input variables
! ---------------------------------------------------------
!
      transfert=.FALSE.
      DO jvar=1,varend
         indvar=var_ord(jvar)
         transfert = ( &
     &       ((FREAL4(var_moy(indvar)).NE.var_moy1(jvar)).AND.(lmoyect)) &
     &       .OR. ((var_moy1(jvar).NE.(FREAL4(0.0))).AND.(.NOT.lmoyect)) &
     &       .OR. transfert )
         transfert = ( &
     &       ((FREAL4(var_ect(indvar)).NE.var_ect1(jvar)).AND.(lmoyect)) &
     &       .OR. ((var_ect1(jvar).NE.(FREAL4(1.0))).AND.(.NOT.lmoyect)) &
     &       .OR. transfert )
      ENDDO
!
! -5.- Read Vx or Vy vector from cpak file
! ----------------------------------------
!
      SELECT CASE (kflagxyo)
!
! Read Vx vector
      CASE(1)
!
! Allocate Vx vector working array (kr4 real kind)
         allocate (ptabx(1:jpxend1), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptabx(:) = FREAL4(0.0)
!
! Read Vx vector
         CALL cdfrpak(kfninsxy,ptabx,1,jpxend1)
!
! Center/reduce input data
         IF (transfert) THEN
            DO jvar=1,varend
               indvar=var_ord(jvar)
               jsdeb=var_ind(indvar)
               jsfin=var_ind(indvar)+var_nbr(indvar)-1
               IF (lmoyect) THEN
                  sxy_ect=FREAL(var_ect1(jvar))/var_ect(indvar)
                  sxy_moy=(FREAL(var_moy1(jvar))-var_moy(indvar)) &
     &                 /var_ect(indvar)
               ELSE
                  sxy_ect=FREAL(var_ect1(jvar))
                  sxy_moy=FREAL(var_moy1(jvar))
               ENDIF
               kvects(jsdeb:jsfin)=FREAL(ptabx(jsdeb:jsfin)) &
     &              *sxy_ect+sxy_moy
           ENDDO
         ELSE
            kvects(:)=FREAL(ptabx(:))
         ENDIF
!
! Read Vy vector
      CASE(2)
!
! Allocate Vx vector working array (kr4 real kind)
         allocate ( ptabx(1:jpyend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptabx(:) = FREAL4(0.0)
!
! Fill progressively Vy vector
         jindybeg=1
         jindyend=jpyend
         IF (tabindxtoy(jindyend).GT.jpxend1) GOTO 1000
!
         jindy=1
         LOOP2 : DO jindxbeg=tabindxtoy(jindybeg), &
     &                       tabindxtoy(jindyend),jpyend
            jindxend=MIN(tabindxtoy(jindyend),jindxbeg-1+jpyend)
            IF ( (tabindxtoy(jindy).LE.jindxend) &
     &           .AND.(tabindxtoy(jindy).GE.jindxbeg) ) THEN
!
               CALL cdfrpak(kfninsxy,ptabx(1),jindxbeg,jindxend)
!
               jindx=tabindxtoy(jindy)-jindxbeg+1
               ltest=((tabindxtoy(jindy)-jindxbeg+1).LE.size(ptabx,1))
               DO WHILE ( (jindy.LE.jindyend) &
     &              .AND.(ltest) )
                  jindx=tabindxtoy(jindy)-jindxbeg+1
                  kvects(jindy)=FREAL(ptabx(jindx))
                  jindy=jindy+1
                  IF (jindy.LE.jindyend) ltest = &
     &                 ((tabindxtoy(jindy)-jindxbeg+1).LE.size(ptabx,1))
               ENDDO
!
            ENDIF
            IF (jindy.GT.jindyend) EXIT LOOP2
         ENDDO LOOP2
!
! Center/reduce input data
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
                  sxy_moy=(FREAL(var_moy1(jvar))-dta_moy(inddta)) &
     &                 /dta_ect(inddta)
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
!
! --- deallocation
      IF (allocated(ptabx)) deallocate(ptabx)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'liocpak','readcpak')
 1001 CALL printerror2(0,1001,3,'liocpak','readcpak')
!
 102  WRITE (texterror,*) 'Parameters in .cpak file header', &
     &     ' incompatible with SESAM configuration'
      CALL printerror2(0,102,3,'liocpak','readcpak',comment=texterror)
 103  WRITE (texterror,*) 'Variable names in .cpak file', &
     &     ' incompatible with SESAM configuration'
      CALL printerror2(0,103,3,'liocpak','readcpak',comment=texterror)
 104  WRITE (texterror,*) 'Invalid input array size'
      CALL printerror2(0,104,1,'liocpak','readcpak',comment=texterror)
 106  WRITE (texterror,*) 'Cannot read Vo vector from cpak file'
      CALL printerror2(0,106,1,'liocpak','readcpak',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readpartcpak(kfninsxy,kvects,kjns,klectinfo,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Read Vx segment from 'cpak' file
!  -------
!  Method : Check cpak file compatibility and
!  ------   read required segment from file
!
!  Input :  kfninsxy  : filename
!  -----    klectinfo : read or not header of 'cpak' file
!           kflagxyo  : vector type (Vx,Vy)
!           kjns      : index of Vx segment to read
!  Output : kvects    : 1D vector (Vx segment)
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : jpxend, jpx, arraynx_jindxbeg, &
     &     arraynx_jpindxend
      use utilcdfpak
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninsxy
      BIGREAL, dimension(:), intent(out) :: kvects
      INTEGER, intent(in) :: kjns
      LOGICAL, intent(in) :: klectinfo
      INTEGER, intent(in) :: kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
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
!----------------------------------------------------------------------
!
! Check input arguments
      jpssize = size(kvects,1)
      IF ((kflagxyo.EQ.1).AND.(jpssize.NE.jpx)) GOTO 104
      IF (kflagxyo.EQ.2) GOTO 105
      IF (kflagxyo.EQ.3) GOTO 106
!
! Control print
      IF (kjns.EQ.1) THEN
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readvar/readpartcpak'
         WRITE(numout,*) '    ==> READING file ',kfninsxy(1:lenv(kfninsxy))
      ENDIF
      ENDIF
!
! -1.- Read 'cpak' file header
! ----------------------------
!
      CALL cdfrhdrpak(kfninsxy,var_nam1,var_dim1, &
     &        var_nbr1,var_moy1,var_ect1)
!
      IF (kjns.EQ.1) THEN
!
! -2.- Read 'cpak' file dimensions
! --------------------------------
!
         CALL cdfrdimpak(kfninsxy,jpxend1, &
     &        varend1,varlg1,cdrec1)
!
! -3.- Check 'cpak' file compatibilty
! -----------------------------------
!
! Compute order of variables fields in 'cpak' file
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
!
! Check if 'cpak' file configuration is identical to SESAM configuration
         sameend = .TRUE.
         sameend = ((varend1.EQ.varend).AND.sameend)
         sameend = ((jpxend1.EQ.jpxend).AND.sameend)
         sameend = ((varlg1.LE.varlg).AND.sameend)
         sameord = .TRUE.
         DO jvar1=1,varend1
            sameord = ((var_ord1(jvar1).EQ.var_ord(jvar1)).AND.sameord)
         ENDDO
!     
         incompatible = (.NOT.(sameord.AND.sameend))
         DO jvar=1,varend
            indvar=var_ord(jvar)
            incompatible = ((var_nam(indvar).NE.var_nam1(jvar)) &
     &           .OR.incompatible)
            incompatible = ((var_dim(indvar).NE.var_dim1(jvar)) &
     &           .OR.incompatible)
            incompatible = ((var_nbr(indvar).NE.var_nbr1(jvar)) &
     &           .OR.incompatible)
         ENDDO
!
         IF (incompatible) GOTO 102
      ENDIF
!
! -4.- Take decision for centering/reducing input variables
! ---------------------------------------------------------
!
      transfert=.FALSE.
      DO jvar=1,varend
         indvar=var_ord(jvar)
         transfert = ( &
     &       ((FREAL4(var_moy(indvar)).NE.var_moy1(jvar)).AND.(lmoyect)) &
     &       .OR. ((var_moy1(jvar).NE.(FREAL4(0.0))).AND.(.NOT.lmoyect)) &
     &       .OR. transfert )
         transfert = ( &
     &       ((FREAL4(var_ect(indvar)).NE.var_ect1(jvar)).AND.(lmoyect)) &
     &       .OR. ((var_ect1(jvar).NE.(FREAL4(1.0))).AND.(.NOT.lmoyect)) &
     &       .OR. transfert )
      ENDDO
!
! -5.- Read Vx segment from cpak file
! -----------------------------------
!
! Allocate Vx vector working array (kr4 real kind)
      allocate ( ptabx(1:jpx), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabx(:) = FREAL4(0.0)
!
! Read Vx segment
      jindxbeg=arraynx_jindxbeg(kjns)
      jindxend=arraynx_jindxbeg(kjns)-1+arraynx_jpindxend(kjns)
      CALL cdfrpak(kfninsxy,ptabx,jindxbeg,jindxend)
!
! Center/reduce input data
      IF (transfert) THEN
         DO jvar=1,varend
            indvar=var_ord(jvar)
            jsdeb=var_ind(indvar)-jindxbeg+1
            jsfin=var_ind(indvar)+var_nbr(indvar)-1-jindxbeg+1
            IF (lmoyect) THEN
               sxy_ect=FREAL(var_ect1(jvar))/var_ect(indvar)
               sxy_moy=(FREAL(var_moy1(jvar))-var_moy(indvar)) &
     &              /var_ect(indvar)
            ELSE
               sxy_ect=FREAL(var_ect1(jvar))
               sxy_moy=FREAL(var_moy1(jvar))
            ENDIF
            IF ((jsfin.GE.1).AND.(jsdeb.LE.jpx)) THEN
               kvects(MAX(jsdeb,1):MIN(jsfin,jpx))= &
     &              FREAL(ptabx(MAX(jsdeb,1):MIN(jsfin,jpx))) &
     &              *sxy_ect+sxy_moy
            ENDIF
         ENDDO
      ELSE
         kvects(:)=FREAL(ptabx(:))
      ENDIF
!
! --- deallocation
      IF (allocated(ptabx)) deallocate(ptabx)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'liocpak','readpartcpak')
 1001 CALL printerror2(0,1001,3,'liocpak','readpartcpak')
!
 102  WRITE (texterror,*) 'Parameters in .cpak file header', &
     &     ' incompatible with SESAM configuration'
      CALL printerror2(0,102,3,'liocpak','readpartcpak',comment=texterror)
 103  WRITE (texterror,*) 'Variable names in .cpak file', &
     &     ' incompatible with SESAM configuration'
      CALL printerror2(0,103,3,'liocpak','readpartcpak',comment=texterror)
 104  WRITE (texterror,*) 'Invalid input array size'
      CALL printerror2(0,104,1,'liocpak','readpartcpak',comment=texterror)
 105  WRITE (texterror,*) 'Cannot read segment of Vy vector from cpak file'
      CALL printerror2(0,105,1,'liocpak','readpartcpak',comment=texterror)
 106  WRITE (texterror,*) 'Cannot read Vo vector from cpak file'
      CALL printerror2(0,106,1,'liocpak','readpartcpak',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writecpak(kfnoutcpak,kvectxin)
!---------------------------------------------------------------------
!
!  Purpose : Write Vx vector in 'cpak' file
!  -------
!  Method : Check if cpak file already exists
!  ------   Create it if necessary
!           Write or overwrite Vx vector in file
!
!  Input : kfnoutcpak : filename
!  -----   kvectxin   : 1D vector (Vx)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend
      use utilcdfpak
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutcpak
      BIGREAL, dimension(:), intent(in) :: kvectxin
!----------------------------------------------------------------------
! local declarations
! ==================
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
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writevar/writecpak'
         WRITE(numout,*) '    ==> WRITING file ',kfnoutcpak(1:lenv(kfnoutcpak))
      ENDIF
!
! Check vector size
      jpxsize = size(kvectxin,1)
      IF (jpxsize.NE.jpxend) GOTO 1000
!
! -1.- Set cpak header variables
! ------------------------------
!
      jpxend1  = jpxsize
      varend1  = varend
      varlg1   = varlg
!
      DO jvar1 = 1,varend
         indvar=var_ord(jvar1)
         var_nam1(jvar1) = var_nam(indvar)
         var_dim1(jvar1) = var_dim(indvar)
         var_nbr1(jvar1) = var_nbr(indvar)
      ENDDO
!
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
!
      WRITE (cdrec1,'(a)') 'Vx =>'
      DO jvar1=1,varend1
         text=cdrec1
         WRITE (cdrec1,'(a,a)') text(1:lenv(text)), &
     &        var_nam1(jvar1)(1:lenv(var_nam1(jvar1)))  
      ENDDO
!
! -2.- Check or write 'cpak' file dimensions
! ------------------------------------------
!
! Check if 'cpak' file already exists
      INQUIRE (FILE=kfnoutcpak,EXIST=existence)
!
      IF (existence) THEN
         CALL cdfrdimpak(kfnoutcpak,jpxend2,varend2,varlg2,cdrec2)
         IF (jpxend2.NE.jpxend1) GOTO 103
         IF (varend2.NE.varend1) GOTO 103
         IF (varlg2.NE.varlg1) GOTO 103
      ELSE
         CALL cdfwdimpak(kfnoutcpak,jpxend1,varend1,varlg1,cdrec1)
      ENDIF
!
! -3.- Write 'cpak' file header
! ------------------------------
!
      CALL cdfwhdrpak(kfnoutcpak,var_nam1,var_dim1,var_nbr1,var_moy1,var_ect1)
!
! -4.- Write Vx vector in cpak file
! ---------------------------------
!
#if ! defined ALLREAL8
      allocate (ptabx(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabx(:) = FREAL4(kvectxin(:))
      CALL cdfwpak(kfnoutcpak,ptabx(1),1,jpxend1)
      IF (allocated(ptabx)) deallocate (ptabx)
#else
      CALL cdfwpak(kfnoutcpak,kvectxin(1),1,jpxend1)
#endif
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'liocpak','writecpak')
 1001 CALL printerror2(0,1001,3,'liocpak','writecpak')
!
 103  WRITE (texterror,*) 'inconsistent dimensions in output cpak file:', &
     &                        kfnoutcpak(1:lenv(kfnoutcpak))
      CALL printerror2(0,103,3,'liocpak','writecpak',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writepartcpak(kfnoutcpak,kvectxin,kjnx)
!---------------------------------------------------------------------
!
!  Purpose : Write Vx segment in 'cpak' file
!  -------
!  Method :
!  ------
!  Input : kfnoutcpak : filename
!  -----   kvectxin   : 1D vector (Vx segment)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend, jpx, arraynx_jindxbeg, &
     &     arraynx_jpindxend
      use utilcdfpak
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutcpak
      BIGREAL, dimension(:), intent(in) :: kvectxin
      INTEGER , intent(in) :: kjnx
!----------------------------------------------------------------------
! local declarations
! ==================
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
!----------------------------------------------------------------------
!
! Control print
      IF (kjnx.EQ.1) THEN
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writevar/writecpak'
         WRITE(numout,*) '    ==> WRITING file ',kfnoutcpak(1:lenv(kfnoutcpak))
      ENDIF
      ENDIF
!
! Check vector size
      jpxsize = size(kvectxin,1)
      IF (jpxsize.NE.jpx) GOTO 1000
!
! Check and write header of cpak file only for 1st segment
      IF (kjnx.EQ.1) THEN
!
! -1.- Set cpak header variables
! ------------------------------
!
         jpxend1  = jpxend
         varend1  = varend
         varlg1   = varlg
         DO jvar1 = 1,varend
            indvar=var_ord(jvar1)
            var_nam1(jvar1) = var_nam(indvar)
            var_dim1(jvar1) = var_dim(indvar)
            var_nbr1(jvar1) = var_nbr(indvar)
         ENDDO
!
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
!
         WRITE (cdrec1,'(a)') 'Vx =>'
         DO jvar1=1,varend1
            text=cdrec1
            WRITE (cdrec1,'(a,a)') text(1:lenv(text)), &
     &           var_nam1(jvar1)(1:lenv(var_nam1(jvar1)))  
         ENDDO
!
! -2.- Check or write 'cpak' file dimensions
! ------------------------------------------
!
         INQUIRE (FILE=kfnoutcpak,EXIST=existence)
!
         IF (existence) THEN
            CALL cdfrdimpak(kfnoutcpak,jpxend2,varend2,varlg2,cdrec2)
            IF (jpxend2.NE.jpxend1) GOTO 103
            IF (varend2.NE.varend1) GOTO 103
            IF (varlg2.NE.varlg1) GOTO 103
         ELSE
            CALL cdfwdimpak(kfnoutcpak,jpxend1,varend1,varlg1,cdrec1)
         ENDIF
!
! -3.- Write 'cpak' file header
! -----------------------------
!
         CALL cdfwhdrpak(kfnoutcpak,var_nam1,var_dim1,var_nbr1,var_moy1,var_ect1)
!
      ENDIF
!
! -4.- Write Vx segment in cpak file
! ----------------------------------
!
      jindxbeg=arraynx_jindxbeg(kjnx)
      jindxend=arraynx_jindxbeg(kjnx)-1+arraynx_jpindxend(kjnx)
!
#if ! defined ALLREAL8
      allocate ( ptabx(1:jpxsize), stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      ptabx(:) = FREAL4(kvectxin(:))
      CALL cdfwpak(kfnoutcpak,ptabx(1),jindxbeg,jindxend)
      IF (allocated(ptabx)) deallocate (ptabx)
#else
      CALL cdfwpak(kfnoutcpak,kvectxin(1),jindxbeg,jindxend)
#endif
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'liocpak','writepartcpak')
 1001 CALL printerror2(0,1001,3,'liocpak','writeparcpak')
!
 103  WRITE (texterror,*) 'inconsistent dimensions in output cpak file:', &
     &                        kfnoutcpak(1:lenv(kfnoutcpak))
      CALL printerror2(0,103,3,'liocpak','writeparcpak',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liocpak

