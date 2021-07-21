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
! ---                    UTILVALID.F90                            ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12 (C.E. Testut)                        ---
! --- modification : 99-11 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- modification : 07-11 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE fildirbas    : Set ensemble member filenames
! --- SUBROUTINE fildirnam    : Set ensemble directory name
! --- FUNCTION existfile      : Check existence of all files 
! ---                           consituting a SESAM object
! --- FUNCTION validmod       : Check validity of module name
! --- FUNCTION validswi       : Check switch validity
! --- FUNCTION validextbas    : Check validity of covariance
! ---                           directory 
! --- FUNCTION validextdtabas : Check validity of Cy directory
! --- FUNCTION validextobsbas : Check validity of Co directory
! --- FUNCTION validextvarbas : Check validity of Cx directory
! --- FUNCTION validextzonbas : Check validity of Cz directory
! --- FUNCTION validextdbs    : Check validity of Io filename
! --- FUNCTION validextdta    : Check validity of Vy filename
! --- FUNCTION validextobs    : Check validity of Vo filename
! --- FUNCTION validextvar    : Check validity of Vx filename
! --- FUNCTION validextzon    : Check validity of Vz filename
! --- FUNCTION validmsk       : Check validity of Vx and Vy masks
! ---                           Check if Vy is included in Vx
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilvalid
      IMPLICIT NONE
      PRIVATE

      PUBLIC fildirbas,fildirnam
      PUBLIC existfile,validmod,validswi,validextbas
      PUBLIC validextdtabas,validextobsbas,validextvarbas
      PUBLIC validextzonbas,validextdbs,validextdta,validextobs
      PUBLIC validextvar,validextzon,validmsk

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE fildirbas (filenamout,dirnamin,kjprout,knumjr,kserie)
!---------------------------------------------------------------------
!
!  Purpose : Create error mode filenames in covariance directory
!  -------
!  Method : Analyse and check covariance directory name
!  ------   Compute error mode filename
!
!  Input : dirnamin   : covariance input directory
!  -----   knumjr     : index of required error mode
!          kserie     : type of files in covariance directory
!                       0=information file series 
!                            *knumjr=1 : valb.txt
!                            *knumjr=2 : valp.txt
!                            *knumjr=3 : vectp.txt
!                       1=mode file series
!                       2=amplitude file series
!  Output : filenamout : filename for required error mode
!  ------   kjprout    : rank(+1) of covariance matrix
!
!---------------------------------------------------------------------
! module
! ======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ==================
      CHARACTER(len=*), intent(out) :: filenamout
      CHARACTER(len=*), intent(in) :: dirnamin
      INTEGER, intent(out) :: kjprout
      INTEGER, intent(in) :: knumjr,kserie
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jext,jextbas,lnumber
      INTEGER :: lextbas,lext,lnum,jtype,xpos
      CHARACTER(len=bgword) :: extbas,ext,num
      INTEGER :: ldir,ldirnum,ldirnumext,ldirnumextbas
      CHARACTER(len=bgword) :: dir,dirnum,dirnumext,dirnumextbas,cform
      INTEGER :: jvar,indvar,jdta,inddta,jobs,indobs,inddbs
      LOGICAL :: existence,existence1,extunit
      CHARACTER(len=1) :: chdigit
!----------------------------------------------------------------------
!
! -1.- Analyse and check directory name
! -------------------------------------
!
! Check validity and type of input directory name
      IF (.NOT.(validextbas(dirnamin))) GOTO 101
!
      jtype = 1
      IF (validextdtabas(dirnamin)) THEN
         jtype=3
      ELSEIF (validextobsbas(dirnamin)) THEN
         jtype=4
      ELSEIF (validextvarbas(dirnamin)) THEN
         jtype=5
      ELSEIF (validextzonbas(dirnamin)) THEN
         jtype=6
      ELSE
         GOTO 101
      ENDIF
!
      ldirnumextbas=lenv(dirnamin)
      WRITE (dirnumextbas,'(a)') dirnamin(1:ldirnumextbas)
!
! Get extension of covariance directory (e.g. .bas)
      jextbas=indext(dirnumextbas,extbastab,nbextbas)
      lextbas=lenv(extbastab(jextbas))
      extbas=extbastab(jextbas)(1:lextbas)
      ldirnumext=ldirnumextbas-lextbas
! Get covariance directory name without type extension (e.g. .bas)
      WRITE (dirnumext,'(a)') dirnumextbas(1:ldirnumext)
!
! Get extension of object filenames (e.g. .dta, .cdf))
      SELECT CASE(jtype)
      CASE(3)
         jext=indext(dirnumext,extdtatab,nbextdta)
         lext=lenv(extdtatab(jext))
         ext=extdtatab(jext)(1:lext)
         extunit=extdtaunit(jext)
      CASE(4)
         jext=indext(dirnumext,extobstab,nbextobs)
         lext=lenv(extobstab(jext))
         ext=extobstab(jext)(1:lext)
         extunit=extobsunit(jext)
      CASE(5)
         jext=indext(dirnumext,extvartab,nbextvar)
         lext=lenv(extvartab(jext))
         ext=extvartab(jext)(1:lext)
         extunit=extvarunit(jext)
      CASE(6)
         jext=indext(dirnumext,extzontab,nbextzon)
         lext=lenv(extzontab(jext))
         ext=extzontab(jext)(1:lext)
         extunit=extzonunit(jext)
      CASE DEFAULT
         GOTO 1000
      END SELECT
! Get covariance directory name without file extension (e.g. .cdf.bas)
      ldirnum=ldirnumext-lext
      WRITE (dirnum,'(a)') dirnumext(1:ldirnum)
!
! Get rank from covariance directory name
      lnumber=nbdigits
      IF (lnumber.GT.9) GOTO 102
      IF (lnumber.LT.1) GOTO 103
      WRITE (chdigit,'(i1.1)') lnumber
      ldir=ldirnum-lnumber
      WRITE (dir,'(a)') dirnum(1:ldir)
      kjprout=mkint( dirnum(ldir+1:ldirnum) )
!
! -1.- Compute error mode filename
! --------------------------------
!
      IF (knumjr.GT.10**nbdigits-1) GOTO 1000
!
      SELECT CASE (knumjr)
      CASE (-3)
          WRITE (filenamout,'(a)') dir(1:ldir)
      CASE (-2)
          WRITE (filenamout,'(a)') dirnum(1:ldirnum)
      CASE (-1)
          WRITE (filenamout,'(a)') dirnumext(1:ldirnumext)
      CASE (0:999999999)
         SELECT CASE (kserie)
         CASE (0)
            SELECT CASE(knumjr)
            CASE (0)
               WRITE (filenamout,'(a)') dirnumextbas(1:ldirnumextbas)
            CASE (1)
! base type information file
               filenamout='valb.txt'
            CASE (2)
! eigenvalue information file
               filenamout='valp.txt'
            CASE (3)
! eigenvector information file
               filenamout='vctp.txt'
            CASE DEFAULT
               GOTO 1000
            END SELECT
         CASE (1)
! mode files series
            IF (extunit) THEN
               cform='("vct",i'//chdigit//'.'//chdigit//',a)'
               WRITE (filenamout,cform) &
     &                 knumjr,ext(1:lext)
            ELSE
               cform='("vct",a1,i'//chdigit//'.'//chdigit//',a)'
               WRITE (filenamout,cform) &
     &                 etoile,knumjr,ext(1:lext)
            ENDIF
         CASE (2)
! amplitude files series
            cform='("apl",i'//chdigit//'.'//chdigit//',".txt")'
            WRITE (filenamout,cform) knumjr
         CASE DEFAULT
            GOTO 1000
         END SELECT
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilbas','fildirbas')
!
 101  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,101,3,'utilbas','fildirbas',comment=texterror)
 102  WRITE (texterror,*) 'Nb of digits in cov. directory name > 9'
      CALL printerror2(0,102,3,'utilbas','fildirbas',comment=texterror)
 103  WRITE (texterror,*) 'Nb of digits in cov. directory name < 1'
      CALL printerror2(0,103,3,'utilbas','fildirbas',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE fildirnam (dirnamout,dirnamin,kjscl)
!---------------------------------------------------------------------
!
!  Purpose : Set directory name for one of the scale in multiple scale ensemble
!  -------
!  Method : Replace star charcater '$' by scale index (kjscl)
!  ------
!
!  Input : dirnamin   : name of ensemble directory (with star)
!  -----   kjscl      : scale index
!
!  Output : dirnamout  : actual name of the directory for the required scale
!  ------
!---------------------------------------------------------------------
! module
! ======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ==================
      CHARACTER(len=*), intent(out) :: dirnamout
      CHARACTER(len=*), intent(in) :: dirnamin
      INTEGER, intent(in) :: kjscl
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpos, jstar, ldirnam
      CHARACTER(len=1) :: chdigit
!----------------------------------------------------------------------
!
! -1.- Analyse and check directory name
! -------------------------------------
!
      ldirnam=lenv(dirnamin)
!
! Get index of star charcater '$' in input directory name
      DO jpos=1,ldirnam
        IF (dirnamin(jpos:jpos).EQ.'@') THEN
          jstar = jpos
          EXIT
        ENDIF
        IF (jpos.EQ.ldirnam) GOTO 101
      ENDDO
!
! Insert scale index instead of the star character
      IF (kjscl.GT.9) GOTO 102
      IF (kjscl.LT.1) GOTO 103
      WRITE (chdigit,'(i1.1)') kjscl
      
      dirnamout=dirnamin(1:jstar-1)//chdigit//dirnamin(jstar+1:ldirnam)

      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilbas','fildirnam')
!
 101  WRITE (texterror,*) 'No $ character in input ensemble name'
      CALL printerror2(0,101,3,'utilbas','fildirnam',comment=texterror)
 102  WRITE (texterror,*) 'Index of scale > 9'
      CALL printerror2(0,102,3,'utilbas','fildirnam',comment=texterror)
 103  WRITE (texterror,*) 'Index of scale < 1'
      CALL printerror2(0,103,3,'utilbas','fildirnam',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION existfile (textin,kflagxyoz)
!
! --- Module declaration
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
! --- Variable declaration
      CHARACTER(len=*), intent(in) :: textin
      INTEGER, intent(in) :: kflagxyoz
      CHARACTER(len=bgword) :: fnamein
      LOGICAL :: existence,existence1
      INTEGER :: jext,xpos,jvar,indvar,jdta,inddta, &
     &     jobs,indobs,inddbs
!
      existence=.TRUE.
      existence1=.FALSE.
!
! Check existence of all files consituting a SESAM object
!
      SELECT CASE(kflagxyoz)
! ==> Vx object (var)
      CASE (1)
         jext=indext(textin,extvartab,nbextvar)
         IF (extvarunit(jext)) THEN
            INQUIRE (FILE=textin,EXIST=existence1)
            existence=(existence1.AND.existence)
         ELSE
            xpos=posit(textin,etoile)
            IF (xpos.LT.1) GOTO 1000
            DO jvar=1,varend
               indvar=var_ord(jvar)
               WRITE(fnamein,'(A,A,A)')  &
     &              textin(1:(xpos-1)), &
     &              varinam(indvar)(1:lenv(varinam(indvar))), &
     &              textin((xpos+1):lenv(textin))
               existence1=.FALSE.
               INQUIRE (FILE=fnamein,EXIST=existence1)
               existence=(existence1.AND.existence)
            ENDDO
         ENDIF
! ==> Vy object (dta)
      CASE (2)
         jext=indext(textin,extdtatab,nbextdta)
         IF (extdtaunit(jext)) THEN
            INQUIRE (FILE=textin,EXIST=existence1)
            existence=(existence1.AND.existence)
         ELSE
            xpos=posit(textin,etoile)
            IF (xpos.LT.1) GOTO 1000
            DO jdta=1,dtaend
               inddta=dta_ord(jdta)
               WRITE(fnamein,'(A,A,A)')  &
     &              textin(1:(xpos-1)), &
     &              dtainam(inddta)(1:lenv(dtainam(inddta))), &
     &              textin((xpos+1):lenv(textin))
               existence1=.FALSE.
               INQUIRE (FILE=fnamein,EXIST=existence1)
               existence=(existence1.AND.existence)
            ENDDO
         ENDIF 
! ==> Vo object (obs)
      CASE (3)
         jext=indext(textin,extobstab,nbextobs)
         IF (extobsunit(jext)) THEN
            INQUIRE (FILE=textin,EXIST=existence1)
            existence=(existence1.AND.existence)
         ELSE
            xpos=posit(textin,etoile)
            IF (xpos.LT.1) GOTO 1000
            DO jobs=1,obsend
               indobs=obs_ord(jobs)
               inddbs=obsnord(jobs)
               WRITE(fnamein,'(A,A,A)')  &
     &              textin(1:(xpos-1)), &
     &              obsinam(indobs,inddbs) &
     &              (1:lenv(obsinam(indobs,inddbs))), &
     &              textin((xpos+1):lenv(textin))
               existence1=.FALSE.
               INQUIRE (FILE=fnamein,EXIST=existence1)
               existence=(existence1.AND.existence)
            ENDDO
         ENDIF 
! ==> Vz object (zon)
      CASE (4)
         jext=indext(textin,extzontab,nbextzon)
         IF (extzonunit(jext)) THEN
            INQUIRE (FILE=textin,EXIST=existence1)
            existence=(existence1.AND.existence)
        ELSE
            xpos=posit(textin,etoile)
            IF (xpos.LT.1) GOTO 1000
            DO jdta=1,dtaend
               inddta=dta_ord(jdta)
               WRITE(fnamein,'(A,A,A)')  &
     &              textin(1:(xpos-1)), &
     &              dtainam(inddta)(1:lenv(dtainam(inddta))), &
     &              textin((xpos+1):lenv(textin))
               existence1=.FALSE.
               INQUIRE (FILE=fnamein,EXIST=existence1)
               existence=(existence1.AND.existence)
            ENDDO
         ENDIF
      CASE DEFAULT 
         GOTO 1000
      END SELECT
! 
      existfile=existence
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilvalid','existfile')
!
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validmod (textin)
!
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jmod
!
      validmod=.FALSE.
      DO jmod = 1,nbmod
         IF (modbool(jmod)) THEN
            validmod = ((textin.EQ.modtab(jmod)).OR.validmod)
         ENDIF
      ENDDO
!
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validswi (textin)
!
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jarg
!
      validswi=.FALSE.
      DO jarg = 1,nbarg
         IF (swibool(jarg)) THEN
            validswi = ((textin.EQ.switab(jarg)).OR.validswi)
         ENDIF
      ENDDO         
!
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextbas (textin)
!
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
!
      validextbas=validextdtabas(textin) &
     &             .OR.validextobsbas(textin) &
     &             .OR.validextvarbas(textin) &
     &             .OR.validextzonbas(textin)
!
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextzonbas (textin)
!
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!
      CHARACTER(len=*), intent(in) :: textin
      LOGICAL :: result
      INTEGER :: jprbas,jextbas,jextext, &
     &     lextbas,ltextin,lextext,lnumber
!
      jextbas=indext(textin,extbastab,nbextbas)
      SELECT CASE (jextbas)
      CASE (1:nbextbas)
         result = .TRUE.
         lextbas=lenv(extbastab(jextbas))
         ltextin=lenv(textin)
         IF (indext(textin(1:(ltextin-lextbas)),extzontab,nbextzon) &
     &        .NE.0) THEN
            jextext=indext(textin(1:(ltextin-lextbas)), &
     &           extzontab,nbextzon)
            lextext=lenv(extzontab(jextext))
            lnumber=nbdigits
            IF ((ltextin-lextbas-lextext-lnumber+1).GT.1) THEN
               jprbas=mkint( textin( (ltextin-lextbas-lextext-lnumber+1) &
     &              :(ltextin-lextbas-lextext)  ) )
               result = .TRUE.
            ELSE
               result = .FALSE.
            ENDIF
         ELSE
            result =.FALSE.
         ENDIF
         result = result.AND.extbasbool(jextbas)
      CASE DEFAULT
         result=.FALSE.            
      END SELECT
!
      validextzonbas=result
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextvarbas (textin)
!
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      LOGICAL :: result
      INTEGER :: jprbas,jextbas,jextext, &
     &     lextbas,ltextin,lextext,lnumber
!
      jextbas=indext(textin,extbastab,nbextbas)
      SELECT CASE (jextbas)
      CASE (1:nbextbas)
         result = .TRUE.
         lextbas=lenv(extbastab(jextbas))
         ltextin=lenv(textin)
         IF (indext(textin(1:(ltextin-lextbas)),extvartab,nbextvar) &
     &        .NE.0) THEN
            jextext=indext(textin(1:(ltextin-lextbas)), &
     &           extvartab,nbextvar)
            lextext=lenv(extvartab(jextext))
            lnumber=nbdigits
            IF ((ltextin-lextbas-lextext-lnumber+1).GT.1) THEN
               jprbas=mkint( textin( (ltextin-lextbas-lextext-lnumber+1) &
     &              :(ltextin-lextbas-lextext)  ) )
               result = .TRUE.
            ELSE
               result = .FALSE.
            ENDIF
         ELSE
            result =.FALSE.
         ENDIF
         result = result.AND.extbasbool(jextbas)           
      CASE DEFAULT
         result=.FALSE.            
      END SELECT
!
      validextvarbas=result
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextobsbas (textin)
!
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      LOGICAL :: result
      INTEGER :: jprbas,jextbas,jextext, &
     &     lextbas,ltextin,lextext,lnumber
!
      jextbas=indext(textin,extbastab,nbextbas)
      SELECT CASE (jextbas)
      CASE (1:nbextbas)
         result = .TRUE.
         lextbas=lenv(extbastab(jextbas))
         ltextin=lenv(textin)
         IF (indext(textin(1:(ltextin-lextbas)),extobstab,nbextobs) &
     &        .NE.0) THEN
            jextext=indext(textin(1:(ltextin-lextbas)), &
     &           extobstab,nbextobs)
            lextext=lenv(extobstab(jextext))
            lnumber=nbdigits
            IF ((ltextin-lextbas-lextext-lnumber+1).GT.1) THEN
               jprbas=mkint( textin( (ltextin-lextbas-lextext-lnumber+1) &
     &              :(ltextin-lextbas-lextext)  ) )
               result = .TRUE.
            ELSE
               result = .FALSE.
            ENDIF
         ELSE
            result =.FALSE.
         ENDIF 
         result = result.AND.extbasbool(jextbas)                    
      CASE DEFAULT
         result=.FALSE.            
      END SELECT
!
      validextobsbas=result
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextdtabas (textin)
!
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      LOGICAL :: result
      INTEGER :: jprbas,jextbas,jextext, &
     &     lextbas,ltextin,lextext,lnumber
!
      jextbas=indext(textin,extbastab,nbextbas)
      SELECT CASE (jextbas)
      CASE (1:nbextbas)
         result = .TRUE.
         lextbas=lenv(extbastab(jextbas))
         ltextin=lenv(textin)
         IF (indext(textin(1:(ltextin-lextbas)),extdtatab,nbextdta) &
     &        .NE.0) THEN
            jextext=indext(textin(1:(ltextin-lextbas)), &
     &           extdtatab,nbextdta)
            lextext=lenv(extdtatab(jextext))
            lnumber=nbdigits
            IF ((ltextin-lextbas-lextext-lnumber+1).GT.1) THEN
               jprbas=mkint( textin( (ltextin-lextbas-lextext-lnumber+1) &
     &              :(ltextin-lextbas-lextext)  ) )
               result = .TRUE.
            ELSE
               result = .FALSE.
            ENDIF
         ELSE
            result =.FALSE.
         ENDIF   
         result = result.AND.extbasbool(jextbas)                       
      CASE DEFAULT
         result=.FALSE.            
      END SELECT
!
      validextdtabas=result
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextdbs (textin)
!
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jextdbs
!
      jextdbs=indext(textin,extdbstab,nbextdbs)
      SELECT CASE (jextdbs)
      CASE (1:nbextdbs)
         validextdbs = extdbsbool(jextdbs)
      CASE DEFAULT
         validextdbs = .FALSE.            
      END SELECT
!
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextdta (textin)
!
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jextdta,spos
!
      jextdta=indext(textin,extdtatab,nbextdta)
      SELECT CASE (jextdta)
      CASE (1:nbextdta)
         IF (extdtaunit(jextdta)) THEN
            validextdta = extdtabool(jextdta)
         ELSE
            spos=posit(textin,etoile)
            validextdta =(((spos.GT.1).AND.(spos.LT.lenv(textin))) &
     &           .AND.(extdtabool(jextdta)))
         ENDIF
      CASE DEFAULT
         validextdta = .FALSE.            
      END SELECT
!     
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextobs (textin)
!
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jextobs,spos
!
      jextobs=indext(textin,extobstab,nbextobs)
      SELECT CASE (jextobs)
      CASE (1:nbextobs)
         IF (extobsunit(jextobs)) THEN
            validextobs = extobsbool(jextobs)
         ELSE
            spos=posit(textin,etoile)
            validextobs = (((spos.GT.1).AND.(spos.LT.lenv(textin))) &
     &           .AND.(extobsbool(jextobs)))
         ENDIF
      CASE DEFAULT
         validextobs = .FALSE.            
      END SELECT
!
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextvar (textin)
!
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jextvar,spos
!
      jextvar=indext(textin,extvartab,nbextvar)
      SELECT CASE (jextvar)
      CASE (1:nbextvar)
         IF (extvarunit(jextvar)) THEN
            validextvar = extvarbool(jextvar)
         ELSE
            spos=posit(textin,etoile)
            validextvar = (((spos.GT.1).AND.(spos.LT.lenv(textin))) &
     &           .AND.(extvarbool(jextvar)))
         ENDIF
      CASE DEFAULT
         validextvar = .FALSE.            
      END SELECT
!
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextzon (textin)
!
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jextzon,spos
!
      jextzon=indext(textin,extzontab,nbextzon)
      SELECT CASE (jextzon)
      CASE (1:nbextzon)
         IF (extzonunit(jextzon)) THEN
            validextzon = extzonbool(jextzon)
         ELSE           
            spos=posit(textin,etoile)
            validextzon = (((spos.GT.1).AND.(spos.LT.lenv(textin))) &
     &           .AND.(extzonbool(jextzon)))
         ENDIF
      CASE DEFAULT
         validextzon = .FALSE.            
      END SELECT
!
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validmsk (jvar,jdta)
!
! --- Module declaration
      use mod_main
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
! --- Variable declaration
      INTEGER, intent(in) :: jvar,jdta
      INTEGER :: indvar,inddta,indvarmsk,inddtamsk
      LOGICAL :: result
      INTEGER :: ji, jj, jk, jt
      INTEGER :: jpifin, jpjfin, jpkfin, jptfin
!
      indvar=var_ord(jvar)
      indvarmsk=jvar-1
      inddta=dta_ord(jdta)
      inddtamsk=jdta-1+varend
! Check variables dimensions
      result=(var_dim(indvar).GE.dta_dim(inddta))
      result=(var_jpi(indvar).GE.dta_jpi(inddta)).AND.result
      IF (dta_dim(inddta).GE.2) THEN
         result=(var_jpj(indvar).GE.dta_jpj(inddta)).AND.result
         IF (dta_dim(inddta).GE.3) THEN
            result=(var_jpk(indvar).GE.dta_jpk(inddta)).AND.result
            IF (dta_dim(inddta).GE.4) THEN
               result=(var_jpt(indvar).GE.dta_jpt(inddta)).AND.result
            ENDIF
         ENDIF
      ENDIF
! Set variable array dimensions
      SELECT CASE (dta_dim(inddta))
      CASE (1)
! ==>  1D
         jpifin = dta_jpi(inddta)
         jpjfin = 1
         jpkfin = 1
         jptfin = 1
      CASE (2)
! ==>  2D
         jpifin = dta_jpi(inddta)
         jpjfin = dta_jpj(inddta)
         jpkfin = 1
         jptfin = 1
      CASE (3)
! ==>  3D
         jpifin = dta_jpi(inddta)
         jpjfin = dta_jpj(inddta)
         jpkfin = dta_jpk(inddta)     
         jptfin = 1
      CASE (4)
! ==>  4D
         jpifin = dta_jpi(inddta)
         jpjfin = dta_jpj(inddta)
         jpkfin = dta_jpk(inddta)  
         jptfin = dta_jpt(inddta)
      CASE DEFAULT
! ==> ERROR
         GOTO 1000
      END SELECT
! Check if Vy vector is included in Vx vector
      jt=1
      DO WHILE ((jt.LE.jptfin).AND.result)
         jk=1
         DO WHILE ((jk.LE.jpkfin).AND.result)
            jj=1
            DO WHILE ((jj.LE.jpjfin).AND.result)
               ji=1
               DO WHILE ((ji.LE.jpifin).AND.result)
                  IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) THEN
                     result=((IBITS(mask(ji,jj,jk,jt),indvarmsk,1).NE.0) &
     &                          .AND.result)
                  ENDIF
                  ji=ji+1
               ENDDO         
               jj=jj+1
            ENDDO                 
            jk=jk+1
         ENDDO         
         jt=jt+1
      ENDDO         
!
      validmsk=result
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,2,'utilvalid','validmsk')
!
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilvalid
