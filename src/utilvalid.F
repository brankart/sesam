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
C ---                    UTILVALID.F                            ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12 (C.E. Testut)                        ---
C --- modification : 99-11 (C.E. Testut)                        ---
C --- modification : 01-06 (C.E. Testut)                        ---
C --- modification : 07-11 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE fildirbas    : Create error mode filenames
C --- FUNCTION existfile      : Check existence of all files 
C ---                           consituting a SESAM object
C --- FUNCTION validmod       : Check validity of module name
C --- FUNCTION validswi       : Check switch validity
C --- FUNCTION validextbas    : Check validity of covariance
C ---                           directory 
C --- FUNCTION validextdtabas : Check validity of Cy directory
C --- FUNCTION validextobsbas : Check validity of Co directory
C --- FUNCTION validextvarbas : Check validity of Cx directory
C --- FUNCTION validextzonbas : Check validity of Cz directory
C --- FUNCTION validextdbs    : Check validity of Io filename
C --- FUNCTION validextdta    : Check validity of Vy filename
C --- FUNCTION validextobs    : Check validity of Vo filename
C --- FUNCTION validextvar    : Check validity of Vx filename
C --- FUNCTION validextzon    : Check validity of Vz filename
C --- FUNCTION validmsk       : Check validity of Vx and Vy masks
C ---                           Check if Vy is included in Vx
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilvalid
      IMPLICIT NONE
      PRIVATE

      PUBLIC fildirbas,existfile,validmod,validswi,validextbas
      PUBLIC validextdtabas,validextobsbas,validextvarbas
      PUBLIC validextzonbas,validextdbs,validextdta,validextobs
      PUBLIC validextvar,validextzon,validmsk

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE fildirbas (filenamout,dirnamin,kjprout,knumjr,kserie)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Create error mode filenames in covariance directory
CCC  -------
CCC  Method : Analyse and check covariance directory name
CCC  ------   Compute error mode filename
CCC
CCC  Input : dirnamin   : covariance input directory
CCC  -----   knumjr     : index of required error mode
CCC          kserie     : type of files in covariance directory
CCC                       0=information file series 
CCC                            *knumjr=1 : valb.txt
CCC                            *knumjr=2 : valp.txt
CCC                            *knumjr=3 : vectp.txt
CCC                       1=mode file series
CCC                       2=amplitude file series
CCC  Output : filenamout : filename for required error mode
CCC  ------   kjprout    : rank(+1) of covariance matrix
CCC
CCC---------------------------------------------------------------------
CC module
CC ======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ==================
      CHARACTER(len=*), intent(out) :: filenamout
      CHARACTER(len=*), intent(in) :: dirnamin
      INTEGER, intent(out) :: kjprout
      INTEGER, intent(in) :: knumjr,kserie
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jext,jextbas,lnumber
      INTEGER :: lextbas,lext,lnum,jtype,xpos
      CHARACTER(len=bgword) :: extbas,ext,num
      INTEGER :: ldir,ldirnum,ldirnumext,ldirnumextbas
      CHARACTER(len=bgword) :: dir,dirnum,dirnumext,dirnumextbas,cform
      INTEGER :: jvar,indvar,jdta,inddta,jobs,indobs,inddbs
      LOGICAL :: existence,existence1,extunit
      CHARACTER(len=1) :: chdigit
CC----------------------------------------------------------------------
C
C -1.- Analyse and check directory name
C -------------------------------------
C
C Check validity and type of input directory name
      IF (.NOT.(validextbas(dirnamin))) GOTO 101
C
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
C
      ldirnumextbas=lenv(dirnamin)
      WRITE (dirnumextbas,'(a)') dirnamin(1:ldirnumextbas)
C
C Get extension of covariance directory (e.g. .bas)
      jextbas=indext(dirnumextbas,extbastab,nbextbas)
      lextbas=lenv(extbastab(jextbas))
      extbas=extbastab(jextbas)(1:lextbas)
      ldirnumext=ldirnumextbas-lextbas
C Get covariance directory name without type extension (e.g. .bas)
      WRITE (dirnumext,'(a)') dirnumextbas(1:ldirnumext)
C
C Get extension of object filenames (e.g. .dta, .cdf))
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
C Get covariance directory name without file extension (e.g. .cdf.bas)
      ldirnum=ldirnumext-lext
      WRITE (dirnum,'(a)') dirnumext(1:ldirnum)
C
C Get rank from covariance directory name
      lnumber=nbdigits
      IF (lnumber.GT.9) GOTO 102
      IF (lnumber.LT.1) GOTO 103
      WRITE (chdigit,'(i1.1)') lnumber
      ldir=ldirnum-lnumber
      WRITE (dir,'(a)') dirnum(1:ldir)
      kjprout=mkint( dirnum(ldir+1:ldirnum) )
C
C -1.- Compute error mode filename
C --------------------------------
C
      IF (knumjr.GT.10**nbdigits-1) GOTO 1000
C
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
C base type information file
               filenamout='valb.txt'
            CASE (2)
C eigenvalue information file
               filenamout='valp.txt'
            CASE (3)
C eigenvector information file
               filenamout='vctp.txt'
            CASE DEFAULT
               GOTO 1000
            END SELECT
         CASE (1)
C mode files series
            IF (extunit) THEN
               cform='("vct",i'//chdigit//'.'//chdigit//',a)'
               WRITE (filenamout,cform)
     $                 knumjr,ext(1:lext)
            ELSE
               cform='("vct",a1,i'//chdigit//'.'//chdigit//',a)'
               WRITE (filenamout,cform)
     $                 etoile,knumjr,ext(1:lext)
            ENDIF
         CASE (2)
C amplitude files series
            cform='("apl",i'//chdigit//'.'//chdigit//',".txt")'
            WRITE (filenamout,cform) knumjr
         CASE DEFAULT
            GOTO 1000
         END SELECT
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilbas','fildirbas')
C
 101  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,101,3,'utilbas','fildirbas',comment=texterror)
 102  WRITE (texterror,*) 'Nb of digits in cov. directory name > 9'
      CALL printerror2(0,102,3,'utilbas','fildirbas',comment=texterror)
 103  WRITE (texterror,*) 'Nb of digits in cov. directory name < 1'
      CALL printerror2(0,103,3,'utilbas','fildirbas',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION existfile (textin,kflagxyoz)
C
C --- Module declaration
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
C --- Variable declaration
      CHARACTER(len=*), intent(in) :: textin
      INTEGER, intent(in) :: kflagxyoz
      CHARACTER(len=bgword) :: fnamein
      LOGICAL :: existence,existence1
      INTEGER :: jext,xpos,jvar,indvar,jdta,inddta,
     $     jobs,indobs,inddbs
C
      existence=.TRUE.
      existence1=.FALSE.
C
C Check existence of all files consituting a SESAM object
C
      SELECT CASE(kflagxyoz)
C ==> Vx object (var)
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
               WRITE(fnamein,'(A,A,A)') 
     $              textin(1:(xpos-1)),
     $              varinam(indvar)(1:lenv(varinam(indvar))),
     $              textin((xpos+1):lenv(textin))
               existence1=.FALSE.
               INQUIRE (FILE=fnamein,EXIST=existence1)
               existence=(existence1.AND.existence)
            ENDDO
         ENDIF
C ==> Vy object (dta)
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
               WRITE(fnamein,'(A,A,A)') 
     $              textin(1:(xpos-1)),
     $              dtainam(inddta)(1:lenv(dtainam(inddta))),
     $              textin((xpos+1):lenv(textin))
               existence1=.FALSE.
               INQUIRE (FILE=fnamein,EXIST=existence1)
               existence=(existence1.AND.existence)
            ENDDO
         ENDIF 
C ==> Vo object (obs)
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
               WRITE(fnamein,'(A,A,A)') 
     $              textin(1:(xpos-1)),
     $              obsinam(indobs,inddbs)
     $              (1:lenv(obsinam(indobs,inddbs))),
     $              textin((xpos+1):lenv(textin))
               existence1=.FALSE.
               INQUIRE (FILE=fnamein,EXIST=existence1)
               existence=(existence1.AND.existence)
            ENDDO
         ENDIF 
C ==> Vz object (zon)
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
               WRITE(fnamein,'(A,A,A)') 
     $              textin(1:(xpos-1)),
     $              dtainam(inddta)(1:lenv(dtainam(inddta))),
     $              textin((xpos+1):lenv(textin))
               existence1=.FALSE.
               INQUIRE (FILE=fnamein,EXIST=existence1)
               existence=(existence1.AND.existence)
            ENDDO
         ENDIF
      CASE DEFAULT 
         GOTO 1000
      END SELECT
C 
      existfile=existence
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilvalid','existfile')
C
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validmod (textin)
C
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jmod
C
      validmod=.FALSE.
      DO jmod = 1,nbmod
         IF (modbool(jmod)) THEN
            validmod = ((textin.EQ.modtab(jmod)).OR.validmod)
         ENDIF
      ENDDO
C
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validswi (textin)
C
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jarg
C
      validswi=.FALSE.
      DO jarg = 1,nbarg
         IF (swibool(jarg)) THEN
            validswi = ((textin.EQ.switab(jarg)).OR.validswi)
         ENDIF
      ENDDO         
C
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextbas (textin)
C
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
C
      validextbas=validextdtabas(textin)
     $             .OR.validextobsbas(textin)
     $             .OR.validextvarbas(textin)
     $             .OR.validextzonbas(textin)
C
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextzonbas (textin)
C
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
C
      CHARACTER(len=*), intent(in) :: textin
      LOGICAL :: result
      INTEGER :: jprbas,jextbas,jextext,
     $     lextbas,ltextin,lextext,lnumber
C
      jextbas=indext(textin,extbastab,nbextbas)
      SELECT CASE (jextbas)
      CASE (1:nbextbas)
         result = .TRUE.
         lextbas=lenv(extbastab(jextbas))
         ltextin=lenv(textin)
         IF (indext(textin(1:(ltextin-lextbas)),extzontab,nbextzon)
     $        .NE.0) THEN
            jextext=indext(textin(1:(ltextin-lextbas)),
     $           extzontab,nbextzon)
            lextext=lenv(extzontab(jextext))
            lnumber=nbdigits
            IF ((ltextin-lextbas-lextext-lnumber+1).GT.1) THEN
               jprbas=mkint( textin( (ltextin-lextbas-lextext-lnumber+1)
     $              :(ltextin-lextbas-lextext)  ) )
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
C
      validextzonbas=result
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextvarbas (textin)
C
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      LOGICAL :: result
      INTEGER :: jprbas,jextbas,jextext,
     $     lextbas,ltextin,lextext,lnumber
C
      jextbas=indext(textin,extbastab,nbextbas)
      SELECT CASE (jextbas)
      CASE (1:nbextbas)
         result = .TRUE.
         lextbas=lenv(extbastab(jextbas))
         ltextin=lenv(textin)
         IF (indext(textin(1:(ltextin-lextbas)),extvartab,nbextvar)
     $        .NE.0) THEN
            jextext=indext(textin(1:(ltextin-lextbas)),
     $           extvartab,nbextvar)
            lextext=lenv(extvartab(jextext))
            lnumber=nbdigits
            IF ((ltextin-lextbas-lextext-lnumber+1).GT.1) THEN
               jprbas=mkint( textin( (ltextin-lextbas-lextext-lnumber+1)
     $              :(ltextin-lextbas-lextext)  ) )
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
C
      validextvarbas=result
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextobsbas (textin)
C
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      LOGICAL :: result
      INTEGER :: jprbas,jextbas,jextext,
     $     lextbas,ltextin,lextext,lnumber
C
      jextbas=indext(textin,extbastab,nbextbas)
      SELECT CASE (jextbas)
      CASE (1:nbextbas)
         result = .TRUE.
         lextbas=lenv(extbastab(jextbas))
         ltextin=lenv(textin)
         IF (indext(textin(1:(ltextin-lextbas)),extobstab,nbextobs)
     $        .NE.0) THEN
            jextext=indext(textin(1:(ltextin-lextbas)),
     $           extobstab,nbextobs)
            lextext=lenv(extobstab(jextext))
            lnumber=nbdigits
            IF ((ltextin-lextbas-lextext-lnumber+1).GT.1) THEN
               jprbas=mkint( textin( (ltextin-lextbas-lextext-lnumber+1)
     $              :(ltextin-lextbas-lextext)  ) )
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
C
      validextobsbas=result
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextdtabas (textin)
C
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      LOGICAL :: result
      INTEGER :: jprbas,jextbas,jextext,
     $     lextbas,ltextin,lextext,lnumber
C
      jextbas=indext(textin,extbastab,nbextbas)
      SELECT CASE (jextbas)
      CASE (1:nbextbas)
         result = .TRUE.
         lextbas=lenv(extbastab(jextbas))
         ltextin=lenv(textin)
         IF (indext(textin(1:(ltextin-lextbas)),extdtatab,nbextdta)
     $        .NE.0) THEN
            jextext=indext(textin(1:(ltextin-lextbas)),
     $           extdtatab,nbextdta)
            lextext=lenv(extdtatab(jextext))
            lnumber=nbdigits
            IF ((ltextin-lextbas-lextext-lnumber+1).GT.1) THEN
               jprbas=mkint( textin( (ltextin-lextbas-lextext-lnumber+1)
     $              :(ltextin-lextbas-lextext)  ) )
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
C
      validextdtabas=result
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextdbs (textin)
C
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jextdbs
C
      jextdbs=indext(textin,extdbstab,nbextdbs)
      SELECT CASE (jextdbs)
      CASE (1:nbextdbs)
         validextdbs = extdbsbool(jextdbs)
      CASE DEFAULT
         validextdbs = .FALSE.            
      END SELECT
C
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextdta (textin)
C
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jextdta,spos
C
      jextdta=indext(textin,extdtatab,nbextdta)
      SELECT CASE (jextdta)
      CASE (1:nbextdta)
         IF (extdtaunit(jextdta)) THEN
            validextdta = extdtabool(jextdta)
         ELSE
            spos=posit(textin,etoile)
            validextdta =(((spos.GT.1).AND.(spos.LT.lenv(textin)))
     $           .AND.(extdtabool(jextdta)))
         ENDIF
      CASE DEFAULT
         validextdta = .FALSE.            
      END SELECT
C     
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextobs (textin)
C
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jextobs,spos
C
      jextobs=indext(textin,extobstab,nbextobs)
      SELECT CASE (jextobs)
      CASE (1:nbextobs)
         IF (extobsunit(jextobs)) THEN
            validextobs = extobsbool(jextobs)
         ELSE
            spos=posit(textin,etoile)
            validextobs = (((spos.GT.1).AND.(spos.LT.lenv(textin)))
     $           .AND.(extobsbool(jextobs)))
         ENDIF
      CASE DEFAULT
         validextobs = .FALSE.            
      END SELECT
C
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextvar (textin)
C
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jextvar,spos
C
      jextvar=indext(textin,extvartab,nbextvar)
      SELECT CASE (jextvar)
      CASE (1:nbextvar)
         IF (extvarunit(jextvar)) THEN
            validextvar = extvarbool(jextvar)
         ELSE
            spos=posit(textin,etoile)
            validextvar = (((spos.GT.1).AND.(spos.LT.lenv(textin)))
     $           .AND.(extvarbool(jextvar)))
         ENDIF
      CASE DEFAULT
         validextvar = .FALSE.            
      END SELECT
C
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validextzon (textin)
C
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: textin
      INTEGER :: jextzon,spos
C
      jextzon=indext(textin,extzontab,nbextzon)
      SELECT CASE (jextzon)
      CASE (1:nbextzon)
         IF (extzonunit(jextzon)) THEN
            validextzon = extzonbool(jextzon)
         ELSE           
            spos=posit(textin,etoile)
            validextzon = (((spos.GT.1).AND.(spos.LT.lenv(textin)))
     $           .AND.(extzonbool(jextzon)))
         ENDIF
      CASE DEFAULT
         validextzon = .FALSE.            
      END SELECT
C
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      LOGICAL FUNCTION validmsk (jvar,jdta)
C
C --- Module declaration
      use mod_main
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
C --- Variable declaration
      INTEGER, intent(in) :: jvar,jdta
      INTEGER :: indvar,inddta,indvarmsk,inddtamsk
      LOGICAL :: result
      INTEGER :: ji, jj, jk, jt
      INTEGER :: jpifin, jpjfin, jpkfin, jptfin
C
      indvar=var_ord(jvar)
      indvarmsk=jvar-1
      inddta=dta_ord(jdta)
      inddtamsk=jdta-1+varend
C Check variables dimensions
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
C Set variable array dimensions
      SELECT CASE (dta_dim(inddta))
      CASE (1)
C ==>  1D
         jpifin = dta_jpi(inddta)
         jpjfin = 1
         jpkfin = 1
         jptfin = 1
      CASE (2)
C ==>  2D
         jpifin = dta_jpi(inddta)
         jpjfin = dta_jpj(inddta)
         jpkfin = 1
         jptfin = 1
      CASE (3)
C ==>  3D
         jpifin = dta_jpi(inddta)
         jpjfin = dta_jpj(inddta)
         jpkfin = dta_jpk(inddta)     
         jptfin = 1
      CASE (4)
C ==>  4D
         jpifin = dta_jpi(inddta)
         jpjfin = dta_jpj(inddta)
         jpkfin = dta_jpk(inddta)  
         jptfin = dta_jpt(inddta)
      CASE DEFAULT
C ==> ERROR
         GOTO 1000
      END SELECT
C Check if Vy vector is included in Vx vector
      jt=1
      DO WHILE ((jt.LE.jptfin).AND.result)
         jk=1
         DO WHILE ((jk.LE.jpkfin).AND.result)
            jj=1
            DO WHILE ((jj.LE.jpjfin).AND.result)
               ji=1
               DO WHILE ((ji.LE.jpifin).AND.result)
                  IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) THEN
                     result=((IBITS(mask(ji,jj,jk,jt),indvarmsk,1).NE.0)
     $                          .AND.result)
                  ENDIF
                  ji=ji+1
               ENDDO         
               jj=jj+1
            ENDDO                 
            jk=jk+1
         ENDDO         
         jt=jt+1
      ENDDO         
C
      validmsk=result
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,2,'utilvalid','validmsk')
C
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilvalid
