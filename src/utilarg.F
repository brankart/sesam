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
C ---                    UTILARG.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12 (C.E. Testut)                        ---
C --- modification : 99-05 (C.E. Testut)                        ---
C --- modification : 03-02 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE affectarg : Check argument validity
C --- SUBROUTINE affectargext : Check file extension validity
C --- SUBROUTINE infoarg : Print information about current action
C --- SUBROUTINE infoargext : Print information about arguments
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilarg
      use mod_main
      use utilvalid
      IMPLICIT NONE
      PRIVATE

      PUBLIC affectarg,affectargext,infoarg,infoargext

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE affectarg (wrd,jmod,jcas,swibool1,modbool1,extbasbool1,
     $     extdbsbool1,extdtabool1,extobsbool1,extvarbool1,extzonbool1)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Check argument validity
CCC  -------
CCC  Method : To each switch correspond particular argument checkings
CCC  ------   Specific error messages are displayed
CCC
CCC  Input : module, switch, argument, and tables of validity
CCC  -----
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: wrd
      INTEGER, intent(in) :: jmod,jcas
      LOGICAL, dimension(:), intent(in) :: swibool1, modbool1
      LOGICAL, dimension(0:), intent(in) :: 
     $     extbasbool1, extdbsbool1, extdtabool1,
     $     extobsbool1, extvarbool1, extzonbool1
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=bgword) :: extension
      INTEGER :: jhelp,jhelp1,jext,ltext
      INTEGER :: jvar,indvar,jobs,indobs,inddbs,flagext
      LOGICAL :: autorise,inbool
CC----------------------------------------------------------------------
C
      SELECT CASE (jcas)
      CASE (1)
C ---  -mode
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largmod) GOTO 102
         IF (.NOT.(validmod(wrd))) GOTO 103
         IF (wrd.NE.modtab(jmod)) GOTO 103
         IF (.NOT.(modbool1(jmod))) GOTO 1000
         argmod = wrd
         largmod=.TRUE.
      CASE (2)
C --- -help (option)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (larghelp) GOTO 102
         jhelp=0
         DO jhelp=1,nbhelp
            IF (wrd.EQ.helptab(jhelp)) THEN
               jhelp1=jhelp
            ENDIF
         ENDDO
         IF (jhelp1.EQ.0) GOTO 103
         arghelp = wrd
         larghelp=.TRUE.
      CASE (3)
C --- -list (option)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (larglist) GOTO 102
         arglist = wrd
         larglist=.TRUE.
      CASE (4)
C --- -varmsk (option)
         IF (wrd.NE.'nomask') THEN
            flagext=5
            inbool=.TRUE.
            CALL affectargext (wrd,jmod,argvarmsk,largvarmsk,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
            SELECT CASE(flagext)
            CASE(5)
               jext=indext(wrd,extvartab,nbextvar)
               IF (.NOT.(extvarbool1(jext))) GOTO 106
               extension=extvartab(jext)
               autorise=((jext.EQ.1).OR.(jext.EQ.2).OR.(jext.EQ.3))
            CASE DEFAULT
               GOTO 1000
            END SELECT
            IF (.NOT.(autorise)) GOTO 106
         ELSE
            IF (largvarmsk) GOTO 102
            argvarmsk = wrd
            largvarmsk=.TRUE.
         ENDIF
      CASE (5)
C --- -dtamsk (option)
         IF (wrd.NE.'nomask') THEN
            flagext=3
            inbool=.TRUE.
            CALL affectargext (wrd,jmod,argdtamsk,largdtamsk,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
            SELECT CASE(flagext)
            CASE(3)
               jext=indext(wrd,extdtatab,nbextdta)
               IF (.NOT.(extdtabool1(jext))) GOTO 106
               extension=extdtatab(jext)
               autorise=((jext.EQ.1).OR.(jext.EQ.2))
            CASE DEFAULT
               GOTO 1000
            END SELECT
            IF (.NOT.(autorise)) GOTO 106
         ELSE
            IF (largdtamsk) GOTO 102
            argdtamsk = wrd
            largdtamsk=.TRUE.
         ENDIF
      CASE (6)
C --- -weight (option)
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argweight,largweight,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.(
     $        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
      CASE (7)
C --- -oestd (option)
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argoestd,largoestd,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.(
     $        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
      CASE (8)
C --- -bias (option)
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argbias,largbias,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.(
     $        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
      CASE (9)
C --- -outinfo (option)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largoutinfo) GOTO 102
         argoutinfo = wrd
         largoutinfo=.TRUE.
      CASE (10)
C --- -inxbas
         flagext=1
         IF (.NOT.validextvarbas(wrd)) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginxbas,larginxbas,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existbas=.TRUE.
      CASE (11)
C --- -indbs
         flagext=2
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argindbs,largindbs,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.2).AND.(.NOT.(
     $        extdbsbool1(indext(wrd,extdbstab,nbextdbs))))) GOTO 106
         existdbs=.TRUE.
      CASE (12)
C --- -indta
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argindta,largindta,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (13)
C --- -inobs
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginobs,larginobs,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.(
     $        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existobs=.TRUE.
      CASE (14)
C --- -invar
         flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginvar,larginvar,
     $        flagext, switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (15)
C --- -inxbasref
         flagext=1
         IF (.NOT.validextvarbas(wrd)) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginxbasref,larginxbasref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existbas=.TRUE.
      CASE (16)
C --- -indbsref
         flagext=2
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argindbsref,largindbsref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.2).AND.(.NOT.(
     $        extdbsbool1(indext(wrd,extdbstab,nbextdbs))))) GOTO 106
         existdbs=.TRUE.
      CASE (17)
C --- -indtaref
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argindtaref,largindtaref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (18)
C --- -inobsref
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginobsref,larginobsref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.(
     $        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existobs=.TRUE.
      CASE (19)
C --- -invarref
         flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginvarref,larginvarref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (20)
C --- -outxbas
         flagext=1
         IF (.NOT.validextvarbas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutxbas,largoutxbas,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existbas=.TRUE.
      CASE (21)
C --- -outdta
         flagext=3
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutdta,largoutdta,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         existdta=.TRUE.
      CASE (22)
C --- -outobs
         flagext=4
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutobs,largoutobs,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.4).AND.(.NOT.(
     $        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
C         existobs=.TRUE.     
      CASE (23)
C --- -outvar
         flagext=5
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutvar,largoutvar,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (24)
C --- -outxbasref
         flagext=1
         IF (.NOT.validextvarbas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutxbasref,largoutxbasref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existbas=.TRUE.
      CASE (25)
C --- -outdtaref
         flagext=3
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutdtaref,largoutdtaref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         existdta=.TRUE.
      CASE (26)
C --- -outobsref
         flagext=4
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutobsref,largoutobsref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.4).AND.(.NOT.(
     $        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
C         existobs=.TRUE.     
      CASE (27)
C --- -outvarref
         flagext=5
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutvarref,largoutvarref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (30)
C --- -reducedta
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argreducedta,largreducedta,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (31)
C --- -outparadap
         argoutparadap=wrd
         largoutparadap=.TRUE.
      CASE (32)
C --- -outrz
C        jext=indext(wrd,extdtatab,nbextdta)
C        IF (jext.EQ.0) GOTO 106
C        IF (.NOT.(extdtabool1(jext))) GOTO 106
C        autorise=((jext.EQ.1).OR.(jext.EQ.2).OR.(jext.EQ.3))
C        IF (.NOT.(autorise)) GOTO 106
         argoutrz=wrd
         largoutrz=.TRUE.
      CASE (33)
C --- -typeoper
         argtypeoper=wrd
         largtypeoper=.TRUE.
      CASE (34)
C --- -typedtadiag
         jext=indext(wrd,extdtatab,nbextdta)
         IF (jext.EQ.0) GOTO 106
         IF (.NOT.(extdtabool1(jext))) GOTO 106
         autorise=((jext.EQ.1).OR.(jext.EQ.2).OR.(jext.EQ.3))
         IF (.NOT.(autorise)) GOTO 106
         argtypedtadiag=wrd
         largtypedtadiag=.TRUE.
      CASE (35)
C --- -diffobsref
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argdiffobsref,largdiffobsref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.(
     $        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
      CASE (36)
C --- -diffdtaref
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argdiffdtaref,largdiffdtaref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (37)
C --- -diffvarref
         flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argdiffvarref,largdiffvarref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (38)
C --- -diffobsorg
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argdiffobsorg,largdiffobsorg,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.(
     $        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
      CASE (39)
C --- -diffdtaorg
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argdiffdtaorg,largdiffdtaorg,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (40)
C --- -diffvarorg
         flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argdiffvarorg,largdiffvarorg,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (41)
C --- -inerrdta
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginerrdta,larginerrdta,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (42)
C --- -outerrdta
         flagext=3
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argouterrdta,largouterrdta,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         existdta=.TRUE.
      CASE (43)
C --- -instddta
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginstddta,larginstddta,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (44)
C --- -outstddta
         flagext=3
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutstddta,largoutstddta,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         existdta=.TRUE.
      CASE (45)
C --- -inybasref
         flagext=1
         IF ((.NOT.validextvarbas(wrd))
     $        .AND.(.NOT.validextdtabas(wrd))) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginybasref,larginybasref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existybas=.TRUE.
      CASE (46)
C --- -outybasref
         flagext=1
         IF (.NOT.validextdtabas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutybasref,largoutybasref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existybas=.TRUE.
      CASE (47)
C --- -outbiasdta
         flagext=3
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutbiasdta,largoutbiasdta,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         existdta=.TRUE.
      CASE (48)
C --- -affectobs
         IF (.NOT.(swibool1(jcas))) GOTO 101
         LOOPAFFECT : DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            ltext=lenv(obs_nam(indobs,inddbs))
            IF ((wrd(1:ltext).EQ.obs_nam(indobs,inddbs)(1:ltext))
     $              .AND.((wrd((ltext+1):(ltext+1))).EQ.(' '))) THEN
               argaffectobs=wrd 
               largaffectobs= .TRUE.
               EXIT LOOPAFFECT
            ENDIF            
         ENDDO LOOPAFFECT
         IF (.NOT.largaffectobs) GOTO 103
      CASE (49)
C --- -nullobs
         IF (.NOT.(swibool1(jcas))) GOTO 101
         LOOPNULL : DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            ltext=lenv(obs_nam(indobs,inddbs))
            IF (((wrd(1:ltext).EQ.obs_nam(indobs,inddbs)(1:ltext))
     $           .AND.((wrd((ltext+1):(ltext+1))).EQ.(' ')))
     $           .OR.(wrd(1:4).EQ.'ALL ')) THEN
               argnullobs=wrd 
               largnullobs= .TRUE.
               EXIT LOOPNULL
            ENDIF            
         ENDDO LOOPNULL
         IF (.NOT.(largnullobs)) GOTO 103
      CASE (50)
C --- -inpartobs
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginpartobs,larginpartobs,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.(
     $        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.(
     $        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         IF (flagext.EQ.4) existobs=.TRUE.
      CASE (51)
C --- -outdiaghst
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largoutdiaghst) GOTO 102
         argoutdiaghst = wrd
         largoutdiaghst=.TRUE.
      CASE (52)
C --- -configobs
         flagext=4
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argconfigobs,largconfigobs,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.4).AND.(.NOT.(
     $        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         existobs=.TRUE.     
      CASE (53)
C --- -fixjpx (option 1)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largfixjpx) GOTO 102
         argfixjpx = wrd
         largfixjpx=.TRUE.
      CASE (54)
C --- -reducevar
         flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argreducevar,largreducevar,
     $        flagext, switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (55)
C --- -scale
         IF (.NOT.(swibool1(jcas))) GOTO 101
         argscale = wrd
         largscale=.TRUE.
      CASE (56)
C --- -biasdbs
         flagext=2
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argbiasdbs,largbiasdbs,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.2).AND.(.NOT.(
     $        extdbsbool1(indext(wrd,extdbstab,nbextdbs))))) GOTO 106
         IF (indext(argbiasdbs,extdbstab,nbextdbs).NE.2) GOTO 106
         existdbs=.TRUE.
      CASE (57)
C --- -inpartvar
         flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginpartvar,larginpartvar,
     $        flagext, switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (58)
C --- -outpartvar
         flagext=5
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutpartvar,largoutpartvar,
     $        flagext, switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.(
     $        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
      CASE (59)
C --- -inzon
         flagext=6
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginzon,larginzon,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.6).AND.(.NOT.(
     $        extzonbool1(indext(wrd,extzontab,nbextzon))))) GOTO 106
         existzon=.TRUE.
      CASE (60)
C --- -outzon
         flagext=6
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutzon,largoutzon,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.6).AND.(.NOT.(
     $        extzonbool1(indext(wrd,extzontab,nbextzon))))) GOTO 106
      CASE (61)
C --- -oecorrel
         flagext=1
         IF (.NOT.validextzonbas(wrd)) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argoecorrel,largoecorrel,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existzbas=.TRUE.
      CASE (62)
C --- -fecorrel
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largfecorrel) GOTO 102
         argfecorrel = wrd
         largfecorrel=.TRUE.
      CASE (63)
C --- -zonindex
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largzonindex) GOTO 102
         argzonindex = wrd
         largzonindex=.TRUE.
      CASE (64)
C --- -incfg
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largincfg) GOTO 102
         argincfg = wrd
         largincfg=.TRUE.
      CASE (65)
C --- -inybas
         flagext=1
         IF ((.NOT.validextvarbas(wrd))
     $        .AND.(.NOT.validextdtabas(wrd))) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginybas,larginybas,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existybas=.TRUE.
      CASE (66)
C --- -ouytbas
         flagext=1
         IF (.NOT.validextdtabas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutybas,largoutybas,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existybas=.TRUE.
      CASE (67)
C --- -inobas
         flagext=1
         IF ((.NOT.validextvarbas(wrd))
     $        .AND.(.NOT.validextdtabas(wrd))
     $        .AND.(.NOT.validextobsbas(wrd))) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginobas,larginobas,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existobas=.TRUE.
      CASE (68)
C --- -outobas
         flagext=1
         IF (.NOT.validextobsbas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutobas,largoutobas,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existobas=.TRUE.
      CASE (69)
C --- -inzbas
         flagext=1
         IF (.NOT.validextzonbas(wrd)) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginzbas,larginzbas,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existzbas=.TRUE.
      CASE (70)
C --- -outzbas
         flagext=1
         IF (.NOT.validextzonbas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutzbas,largoutzbas,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existzbas=.TRUE.
      CASE (71)
C --- -configzon
         flagext=6
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argconfigzon,largconfigzon,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.6).AND.(.NOT.(
     $        extzonbool1(indext(wrd,extzontab,nbextzon))))) GOTO 106
         existzon=.TRUE.     
      CASE (72)
C --- -coefrmax
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largcoefrmax) GOTO 102
         argcoefrmax = wrd
         largcoefrmax=.TRUE.
      CASE (73)
C --- -disable
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largdisable) GOTO 102
         argdisable = wrd
         largdisable=.TRUE.
      CASE (74)
C --- -inobasref
         flagext=1
         IF ((.NOT.validextvarbas(wrd))
     $        .AND.(.NOT.validextdtabas(wrd))
     $        .AND.(.NOT.validextobsbas(wrd))) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginobasref,larginobasref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existobas=.TRUE.
      CASE (75)
C --- -fixjpz (option 1)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largfixjpz) GOTO 102
         argfixjpz = wrd
         largfixjpz=.TRUE.
      CASE (76)
C --- -inparadap (option 2)
         arginparadap = wrd
         larginparadap=.TRUE.
      CASE (77)
C --- -fixjpu (option 1)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largfixjpu) GOTO 102
         argfixjpu = wrd
         largfixjpu=.TRUE.
      CASE (78)
C --- -inptzon
         flagext=6
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginptzon,larginptzon,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.6).AND.(.NOT.(
     $        extzonbool1(indext(wrd,extzontab,nbextzon))))) GOTO 106
         existzon=.TRUE.
      CASE (79)
C --- -outptzon
         flagext=6
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutptzon,largoutptzon,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.6).AND.(.NOT.(
     $        extzonbool1(indext(wrd,extzontab,nbextzon))))) GOTO 106
      CASE (80)
C --- -action (option 1)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largaction) GOTO 102
         argaction = wrd
         largaction=.TRUE.
      CASE (81)
C --- -inzbasref
         flagext=1
         IF (.NOT.validextzonbas(wrd)) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginzbasref,larginzbasref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existzbas=.TRUE.
      CASE (82)
C --- -inzonref
         flagext=6
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginzonref,larginzonref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.6).AND.(.NOT.(
     $        extzonbool1(indext(wrd,extzontab,nbextzon))))) GOTO 106
         existzon=.TRUE.
      CASE (83)
C --- -inoptcfg
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (larginoptcfg) GOTO 102
         arginoptcfg = wrd
         larginoptcfg=.TRUE.
      CASE (84)
C --- -connect
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largconnect) GOTO 102
         argconnect = wrd
         largconnect=.TRUE.
      CASE (85)
C --- -incstr
         flagext=1
         IF (.NOT.validextvarbas(wrd)) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argincstr,largincstr,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existcstr=.TRUE.
      CASE (86)
C --- -outobasref
         flagext=1
         IF (.NOT.validextobsbas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutobasref,largoutobasref,
     $        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.(
     $        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existobas=.TRUE.
      CASE (87)
C --- -insmocfg
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (larginsmocfg) GOTO 102
         arginsmocfg = wrd
         larginsmocfg=.TRUE.
      CASE (88)
C --- -outsmocfg
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largoutsmocfg) GOTO 102
         argoutsmocfg = wrd
         largoutsmocfg=.TRUE.
      CASE (89)
C --- -inrz
         arginrz=wrd
         larginrz=.TRUE.
      CASE DEFAULT
         GOTO 1000
      END SELECT  
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(jmod,1000,1,'readarg','affectarg')
C
 101   WRITE (texterror,*) 'Switch ',
     $     switab(jcas)(1:lenv(switab(jcas))),
     $     ' is not available for this module (see -help swi)'
      CALL printerror2(jmod,101,3,'readarg','affectarg',
     $      comment=texterror)
 102   WRITE (texterror,*) 'The same switch ',
     $     switab(jcas)(1:lenv(switab(jcas))),
     $     ' cannot be used more than once'
      CALL printerror2(jmod,102,3,'readarg','affectarg',
     $      comment=texterror)
 103   WRITE (texterror,*) 'Argument ',
     $     switab(jcas)(1:lenv(switab(jcas))),
     $     ' is not valid'
      CALL printerror2(jmod,103,3,'readarg','affectarg',
     $      comment=texterror)
 106   WRITE (texterror,*) 'File extension ',extension(1:lenv(extension)),
     $     ' is not available',
     $     ' with switch ',switab(jcas)(1:lenv(switab(jcas))),
     $     ' (see -help ext)'
      CALL printerror2(jmod,106,3,'readarg','affectarg',
     $      comment=texterror)

      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE affectargext (wrd,jmod,argument,largument,
     $     flagext,swi,lswi,inbool)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Check validity of file extensions
CCC  -------
CCC  Method : To each object type correspond particular checkings
CCC  ------   Specific error messages are displayed
CCC
CCC  Input : module, switch, argument, object type (flagext),
CCC  ------  input or output object (inbool)
CCC
CCC  Output : activate argument in SESAM: argument,largument
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: wrd,swi
      CHARACTER(len=*), intent(out) :: argument
      INTEGER, intent(in) :: flagext,jmod
      LOGICAL, intent(inout) :: largument
      LOGICAL, intent(in) :: lswi,inbool
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=bgword) :: extension
      INTEGER :: jext
CC----------------------------------------------------------------------
C
      IF (largument) GOTO 102
      SELECT CASE (flagext)
      CASE (1)
C --- bas
         IF ( (.NOT.(validextvarbas(wrd)))
     $        .AND.(.NOT.(validextdtabas(wrd)))
     $        .AND.(.NOT.(validextobsbas(wrd)))
     $        .AND.(.NOT.(validextzonbas(wrd)))
     $        ) GOTO 103
         jext=indext(wrd,extbastab,nbextbas)
         extension=extbastab(jext)
         IF (.NOT.(extbasbool(jext))) GOTO 104
         IF (((extbasmem(jext)/10).EQ.0).AND.(inbool)) GOTO 106
         IF (((MOD(extbasmem(jext),10)).EQ.0).AND.(.NOT.inbool)) GOTO 107
      CASE (2)
C --- dbs
         IF (.NOT.(validextdbs(wrd))) GOTO 103
         jext=indext(wrd,extdbstab,nbextdbs)
         extension=extdbstab(jext)
         IF (.NOT.(extdbsbool(jext))) GOTO 104
         IF (((extdbsmem(jext)/10).EQ.0).AND.(inbool)) GOTO 106
         IF (((MOD(extdbsmem(jext),10)).EQ.0).AND.(.NOT.inbool)) GOTO 107
      CASE (3)
C --- dta
         IF (.NOT.(validextdta(wrd))) GOTO 103
         jext=indext(wrd,extdtatab,nbextdta)
         extension=extdtatab(jext)
         IF (.NOT.(extdtabool(jext))) GOTO 104
         IF (((extdtamem(jext)/10).EQ.0).AND.(inbool)) GOTO 106
         IF (((MOD(extdtamem(jext),10)).EQ.0).AND.(.NOT.inbool)) GOTO 107
      CASE (4)
C --- obs
         IF (.NOT.(validextobs(wrd))) GOTO 103
         jext=indext(wrd,extobstab,nbextobs)
         extension=extobstab(jext)
         IF (.NOT.(extobsbool(jext))) GOTO 104
         IF (((extobsmem(jext)/10).EQ.0).AND.(inbool)) GOTO 106
         IF (((MOD(extobsmem(jext),10)).EQ.0).AND.(.NOT.inbool)) GOTO 107
      CASE (5)
C --- var
         IF (.NOT.(validextvar(wrd))) GOTO 103
         jext=indext(wrd,extvartab,nbextvar)
         extension=extvartab(jext)
         IF (.NOT.(extvarbool(jext))) GOTO 104
         IF (((extvarmem(jext)/10).EQ.0).AND.(inbool)) GOTO 106
         IF (((MOD(extvarmem(jext),10)).EQ.0).AND.(.NOT.inbool)) GOTO 107
      CASE (6)
C --- zon
         IF (.NOT.(validextzon(wrd))) GOTO 103
         jext=indext(wrd,extzontab,nbextzon)
         extension=extzontab(jext)
         IF (.NOT.(extzonbool(jext))) GOTO 104
         IF (((extzonmem(jext)/10).EQ.0).AND.(inbool)) GOTO 106
         IF (((MOD(extzonmem(jext),10)).EQ.0).AND.(.NOT.inbool)) GOTO 107
C
      CASE DEFAULT
         GOTO 1000
      END SELECT
      argument = wrd
      largument=.TRUE.
      RETURN
C
C --- error management
C
 1000 CALL printerror2(jmod,1000,1,'readarg','affectargext')
C
 102  WRITE (texterror,*) 'The same switch ',
     $     swi(1:lenv(swi)),
     $     ' cannot be used more than once'
      CALL printerror2(jmod,102,3,'readarg','affectargext',
     $      comment=texterror)
 103   WRITE (texterror,*) 'Argument ',
     $      wrd(1:lenv(wrd)),
     $      ' is not valid'
      CALL printerror2(jmod,103,3,'readarg','affectargext',
     $      comment=texterror)
 104  WRITE (texterror,*) 'Extension ',extension(1:lenv(extension)),
     $     ' is not available with switch ',swi(1:lenv(swi)),
     $     ' (see -help ext)'
      CALL printerror2(jmod,104,3,'readarg','affectargext',
     $      comment=texterror)
 106  WRITE (texterror,*) 'Extension ',extension(1:lenv(extension)),
     $     ' is not valid with switch ',swi(1:lenv(swi)),
     $     ' for input files (see -help ext)' 
      CALL printerror2(jmod,106,3,'readarg','affectargext',
     $      comment=texterror)
 107  WRITE (texterror,*) 'Extension ',extension(1:lenv(extension)),
     $     ' is not valid with switch ',swi(1:lenv(swi)),
     $     ' for output files (see -help ext)' 
      CALL printerror2(jmod,107,3,'readarg','affectargext',
     $      comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE infoarg (jmod,jaction,swibool1,modbool1,extbasbool1,
     $     extdbsbool1,extdtabool1,extobsbool1,extvarbool1,extzonbool1)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Print information about current action
CCC  -------
CCC  Method : Loop over all user switches and write
CCC  ------   information about them
CCC
CCC  Input : mode, action and tables of usage
CCC  -----
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      INTEGER, intent(in) :: jmod,jaction
      LOGICAL, dimension(:), intent(in) :: swibool1, modbool1
      LOGICAL, dimension(0:), intent(in) :: 
     $     extbasbool1, extdbsbool1, extdtabool1,
     $     extobsbool1, extvarbool1, extzonbool1
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=bgword) :: argument
      INTEGER :: jcas
      INTEGER :: jext,xypos,flagext
      INTEGER :: jvar,indvar,jdta,inddta
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         DO jcas=1,nbarg
            jext = 0
            xypos= 0
            SELECT CASE (jcas)
            CASE (1)
C --- -mode
               IF (largmod) THEN
                  WRITE(numout,*) ' -----------------------------'
                  WRITE(numout,*) ' module = ',
     $                 modtab(jmod)(1:lenv(modtab(jmod)))
                  WRITE(numout,*) ' action = ',jaction
                  WRITE(numout,*) ' -----------------------------'
              ENDIF
            CASE (2)
C --- -help
               IF (larghelp) THEN
                  WRITE(numout,*) ' -----------------------------'
                  WRITE(numout,*) ' help  = ',
     $                 arghelp(1:lenv(arghelp))
                  WRITE(numout,*) ' -----------------------------'
               ENDIF
            CASE (3)
C --- -list (option)
               IF (larglist) THEN
                  WRITE(numout,*) ' namelist (option name) = ',
     $                 arglist(1:lenv(arglist))
               ELSE
                  WRITE(numout,*) ' namelist (default name)  = ',
     $                 fnamlist(1:lenv(fnamlist))
               ENDIF
            CASE (4)
C --- -varmsk (option)
               IF (largvarmsk) THEN
                  argument=argvarmsk
                  WRITE(numout,*) ' Vx mask (option name) ='
                  IF (argvarmsk(1:lenv(argvarmsk)).EQ.'nomask') THEN
                     WRITE(numout,*) '      no mask read !'
                  ELSE
                     IF (.NOT.(validextvar(argvarmsk))) GOTO 101
                     jext=indext(argvarmsk,extvartab,nbextvar)
                     IF (.NOT.(extvarunit(jext))) THEN
                        xypos=posit(argvarmsk,etoile)
                        DO jvar=1,varend
                           indvar=var_ord(jvar)
                           WRITE (clname,'(6x,a,a,a)') argvarmsk(1:(xypos-1)),
     $                       varinam(indvar)(1:lenv(varinam(indvar))),
     $                       argvarmsk((xypos+1):lenv(argvarmsk))
                           IF (jvar.EQ.1) THEN
                              WRITE(numout,*) '      mask(',jvar,'<=>',
     $                          var_nam(indvar)(1:lenv(var_nam(indvar)))
     $                          ,')  = ',clname(1:lenv(clname))
                           ELSE
                              WRITE(numout,*) '       "" (',jvar,'<=>',
     $                          var_nam(indvar)(1:lenv(var_nam(indvar)))
     $                          ,')  = ',clname(1:lenv(clname))
                           ENDIF
                        ENDDO
                     ELSE
                        WRITE(numout,*) '      mask  = ',
     $                    argvarmsk(1:lenv(argvarmsk))
                     ENDIF
                  ENDIF
               ELSE
                  WRITE(numout,*) ' Vx mask (default name) ='
                  DO jvar=1,varend
                     indvar=var_ord(jvar)
                     IF (jvar.EQ.1) THEN
                        WRITE(numout,*) '      mask(',jvar,'<=>',
     $                       var_nam(indvar)(1:lenv(var_nam(indvar))),
     $                       ')  = ',
     $                       varfmsk(indvar)(1:lenv(varfmsk(indvar)))
                     ELSE
                        WRITE(numout,*) '       "" (',jvar,'<=>',
     $                       var_nam(indvar)(1:lenv(var_nam(indvar))),
     $                       ')  = ',
     $                       varfmsk(indvar)(1:lenv(varfmsk(indvar)))
                        ENDIF
                     ENDDO
               ENDIF
            CASE (5)
C --- -dtamsk (option)
               IF (largdtamsk) THEN
                  argument=argdtamsk
                  WRITE(numout,*) ' Vy mask (option name) ='
                  IF (argdtamsk(1:lenv(argdtamsk)).EQ.'nomask') THEN
                     WRITE(numout,*) '      no mask read !'
                  ELSE
                     IF (.NOT.(validextdta(argdtamsk))) GOTO 101
                     jext=indext(argdtamsk,extdtatab,nbextdta)
                     IF (.NOT.(extdtaunit(jext))) THEN
                        xypos=posit(argdtamsk,etoile)
                        DO jdta=1,dtaend
                           inddta=dta_ord(jdta)
                           WRITE (clname,'(6x,a,a,a)') argdtamsk(1:(xypos-1)),
     $                       dtainam(inddta)(1:lenv(dtainam(inddta))),
     $                       argdtamsk((xypos+1):lenv(argdtamsk))
                           IF (jdta.EQ.1) THEN
                              WRITE(numout,*) '      mask(',jdta,'<=>',
     $                          dta_nam(inddta)(1:lenv(dta_nam(inddta)))
     $                          ,')  = ',clname(1:lenv(clname))
                           ELSE
                              WRITE(numout,*) '       "" (',jdta,'<=>',
     $                          dta_nam(inddta)(1:lenv(dta_nam(inddta)))
     $                          ,')  = ',clname(1:lenv(clname))
                           ENDIF
                        ENDDO
                     ELSE
                        WRITE(numout,*) '      mask  = ',
     $                    argdtamsk(1:lenv(argdtamsk))
                     ENDIF
                  ENDIF
               ELSE
                  WRITE(numout,*) ' Vy mask (default name) ='
                  DO jdta=1,dtaend
                     inddta=dta_ord(jdta)
                     IF (jdta.EQ.1) THEN
                        WRITE(numout,*) '      mask(',jdta,'<=>',
     $                       dta_nam(inddta)(1:lenv(dta_nam(inddta))),
     $                       ')  = ',
     $                       dtafmsk(inddta)(1:lenv(dtafmsk(inddta)))
                     ELSE
                        WRITE(numout,*) '       "" (',jdta,'<=>',
     $                       dta_nam(inddta)(1:lenv(dta_nam(inddta))),
     $                       ')  = ',
     $                       dtafmsk(inddta)(1:lenv(dtafmsk(inddta)))
                        ENDIF
                     ENDDO
               ENDIF
            CASE (6)
C --- -weight (option)
               IF (largweight) THEN
                  flagext=5
                  IF (validextobs(argweight)) flagext=4
                  IF (validextdta(argweight)) flagext=3
                  CALL infoargext(argweight,'argweight',flagext,jmod)
                  IF (flagext.EQ.5) THEN
                    IF (extvarmem(indext(argweight,extvartab,nbextvar))/10
     $                   .GT.nallmem) GOTO 102
                  ENDIF
               ENDIF
            CASE (7)
C --- -oestd (option)
               IF (largoestd) THEN
                  flagext=5
                  IF (validextobs(argoestd)) flagext=4
                  IF (validextdta(argoestd)) flagext=3
                  CALL infoargext(argoestd,'argoestd',flagext,jmod)
                  IF (flagext.EQ.5) THEN
                    IF (extvarmem(indext(argoestd,extvartab,nbextvar))/10
     $                   .GT.nallmem) GOTO 102
                  ENDIF
               ENDIF
            CASE (8)
C --- -bias (option)
               IF (largbias) THEN
                  flagext=5
                  IF (validextobs(argbias)) flagext=4
                  IF (validextdta(argbias)) flagext=3
                  CALL infoargext(argbias,'argbias',flagext,jmod)
                  IF (flagext.EQ.5) THEN
                    IF (extvarmem(indext(argbias,extvartab,nbextvar))/10
     $                   .GT.nallmem) GOTO 102
                  ENDIF
               ENDIF
            CASE (9)
C --- -outinfo (option)
               IF (largoutinfo) THEN
                  WRITE(numout,*) ' listing file (option name) = ',
     $                 argoutinfo(1:lenv(argoutinfo))
               ENDIF
            CASE (10)
C --- -inxbas
               IF (larginxbas) THEN
                  flagext=1
                  CALL infoargext(arginxbas,'arginxbas',flagext,jmod)
                  IF (extvarmem(indext(arginxbas(1:lenv(arginxbas)-
     $                 lenv(extbastab(indext(arginxbas,extbastab,nbextbas)))),
     $                 extvartab, nbextvar ))/10.GT.nallmem) GOTO 103
               ENDIF
            CASE (11)
C --- -indbs
               IF (largindbs) THEN
                  flagext=2
                  CALL infoargext(argindbs,'argindbs',flagext,jmod)
               ENDIF
            CASE (12)
C --- -indta
               IF (largindta) THEN
                  flagext=3
                  IF (validextvar(argindta)) flagext=5
                  CALL infoargext(argindta,'argindta',flagext,jmod)
               ENDIF
            CASE (13)
C --- -inobs
               IF (larginobs) THEN
                  flagext=4
                  IF (validextdta(arginobs)) flagext=3
                  IF (validextvar(arginobs)) flagext=5
                  CALL infoargext(arginobs,'arginobs',flagext,jmod)
               ENDIF
            CASE (14)
C --- -invar
               IF (larginvar) THEN
                  flagext=5
                  CALL infoargext(arginvar,'arginvar',flagext,jmod)
                  IF (extvarmem(indext(arginvar,extvartab,nbextvar))/10
     $                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (15)
C --- -inxbasref
               IF (larginxbasref) THEN
                  flagext=1
                  CALL infoargext(arginxbasref,'arginxbasref',flagext,
     $                 jmod)
                  IF (extvarmem(indext(arginxbasref(1:lenv(arginxbasref)-
     $                 lenv(extbastab(indext(arginxbasref,extbastab,nbextbas)))),
     $                 extvartab, nbextvar ))/10.GT.nallmem) GOTO 103
               ENDIF
            CASE (16)
C --- -indbsref
               IF (largindbsref) THEN
                  flagext=2
                  CALL infoargext(argindbsref,'argindbsref',flagext,jmod)
               ENDIF
            CASE (17)
C --- -indtaref
               IF (largindtaref) THEN
                  flagext=3
                  IF (validextvar(argindtaref)) flagext=5
                  CALL infoargext(argindtaref,'argindtaref',flagext,jmod)
               ENDIF
            CASE (18)
C --- -inobsref
               IF (larginobsref) THEN
                  flagext=4
                  IF (validextdta(arginobsref)) flagext=3
                  IF (validextvar(arginobsref)) flagext=5
                  CALL infoargext(arginobsref,'arginobsref',flagext,jmod)
               ENDIF
            CASE (19)
C --- -invarref
               IF (larginvarref) THEN
                  flagext=5
                  CALL infoargext(arginvarref,'arginvarref',flagext,jmod)
                  IF (extvarmem(indext(arginvarref,extvartab,nbextvar))/10
     $                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (20)
C --- -outxbas
               IF (largoutxbas) THEN
                  flagext=1
                  CALL infoargext(argoutxbas,'argoutxbas',flagext,jmod)
                  IF (extvarmem(indext(argoutxbas(1:lenv(argoutxbas)-
     $                 lenv(extbastab(indext(argoutxbas,extbastab,nbextbas)))),
     $                 extvartab, nbextvar ))/10.GT.nallmem) GOTO 103
               ENDIF
            CASE (21)
C --- -outdta
               IF (largoutdta) THEN
                  flagext=3
                  CALL infoargext(argoutdta,'argoutdta',flagext,jmod)
               ENDIF
            CASE (22)
C --- -outobs
               IF (largoutobs) THEN
                  flagext=4
                  CALL infoargext(argoutobs,'argoutobs',flagext,jmod)
               ENDIF
            CASE (23)
C --- -outvar
               IF (largoutvar) THEN
                  flagext=5
                  CALL infoargext(argoutvar,'argoutvar',flagext,jmod)
                  IF (extvarmem(indext(argoutvar,extvartab,nbextvar))/10
     $                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (24)
C --- -outxbasref
               IF (largoutxbasref) THEN
                  flagext=1
                  CALL infoargext(argoutxbasref,'argoutxbasref',flagext,jmod)
                  IF (extvarmem(indext(argoutxbasref(1:lenv(argoutxbasref)-
     $                 lenv(extbastab(indext(argoutxbasref,extbastab,nbextbas)))),
     $                 extvartab, nbextvar ))/10.GT.nallmem) GOTO 103
               ENDIF
            CASE (25)
C --- -outdtaref
               IF (largoutdtaref) THEN
                  flagext=3
                  CALL infoargext(argoutdtaref,'argoutdtaref',flagext,jmod)
               ENDIF
            CASE (26)
C --- -outobsref
               IF (largoutobsref) THEN
                  flagext=4
                  CALL infoargext(argoutobsref,'argoutobsref',flagext,jmod)
               ENDIF
            CASE (27)
C --- -outvarref
               IF (largoutvarref) THEN
                  flagext=5
                  CALL infoargext(argoutvarref,'argoutvarref',flagext,jmod)
                  IF (extvarmem(indext(argoutvarref,extvartab,nbextvar))/10
     $                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (30)
C --- -reducedta
               IF (largreducedta) THEN
                  flagext=3
                  IF (validextvar(argreducedta)) flagext=5
                  CALL infoargext(argreducedta,'argreducedta',
     $                 flagext,jmod)
               ENDIF
            CASE (31)
C --- -outparadap
              IF (largoutparadap) THEN
                  WRITE(numout,*) ' output adaptive parameter file  = ',
     $                 argoutparadap(1:lenv(argoutparadap))
              ENDIF
            CASE (32)
C --- -outrz
              IF (largoutrz) THEN
                  WRITE(numout,*) ' reduced space file  = ',
     $                 argoutrz(1:lenv(argoutrz))
              ENDIF
            CASE (33)
C --- -typeoper
              IF (largtypeoper) THEN
                  WRITE(numout,*) ' operation type = ',
     $                 argtypeoper(1:lenv(argtypeoper))
              ENDIF
            CASE (34)
C --- -typedtadiag
              IF (largtypedtadiag) THEN
                  WRITE(numout,*) ' type of Vy diagnostic files  = ',
     $                 argtypedtadiag(1:lenv(argtypedtadiag))
              ENDIF
            CASE (35)
C --- -diffobsref
               IF (largdiffobsref) THEN
                  flagext=4
                  IF (validextvar(argdiffobsref)) flagext=5
                  IF (validextdta(argdiffobsref)) flagext=3
                  CALL infoargext(argdiffobsref,'argdiffobsref',flagext,jmod)
               ENDIF
            CASE (36)
C --- -diffdtaref
               IF (largdiffdtaref) THEN
                  flagext=3
                  IF (validextvar(argdiffdtaref)) flagext=5
                  CALL infoargext(argdiffdtaref,'argdiffdtaref',flagext,jmod)
               ENDIF
            CASE (37)
C --- -diffvarref
               IF (largdiffvarref) THEN
                  flagext=5
                  CALL infoargext(argdiffvarref,'argdiffvarref',flagext,jmod)
               ENDIF
            CASE (38)
C --- -diffobsorg
               IF (largdiffobsorg) THEN
                  flagext=4
                  IF (validextvar(argdiffobsorg)) flagext=5
                  IF (validextdta(argdiffobsorg)) flagext=3
                  CALL infoargext(argdiffobsorg,'argdiffobsorg',flagext,jmod)
               ENDIF
            CASE (39)
C --- -diffdtaorg
               IF (largdiffdtaorg) THEN
                  flagext=3
                  IF (validextvar(argdiffdtaorg)) flagext=5
                  CALL infoargext(argdiffdtaorg,'argdiffdtaorg',flagext,jmod)
               ENDIF
            CASE (40)
C --- -diffvarorg
               IF (largdiffvarorg) THEN
                  flagext=5
                  CALL infoargext(argdiffvarorg,'argdiffvarorg',flagext,jmod)
               ENDIF
            CASE (41)
C --- -inerrdta
               IF (larginerrdta) THEN
                  flagext=3
                  IF (validextvar(arginerrdta)) flagext=5
                  CALL infoargext(arginerrdta,'arginerrdta',flagext,jmod)
               ENDIF
            CASE (42)
C --- -outerrdta
               IF (largouterrdta) THEN
                  flagext=3
                  CALL infoargext(argouterrdta,'argouterrdta',flagext,jmod)
               ENDIF
            CASE (43)
C --- -instddta
               IF (larginstddta) THEN
                  flagext=3
                  IF (validextvar(arginstddta)) flagext=5
                  CALL infoargext(arginstddta,'arginstddta',flagext,jmod)
               ENDIF
            CASE (44)
C --- -outstddta
               IF (largoutstddta) THEN
                  flagext=3
                  CALL infoargext(argoutstddta,'argoutstddta',flagext,jmod)
               ENDIF
            CASE (45)
C --- -inybasref
               IF (larginybasref) THEN
                  flagext=1
                  CALL infoargext(arginybasref,'arginybasref',flagext,jmod)
               ENDIF
            CASE (46)
C --- -outybasref
               IF (largoutybasref) THEN
                  flagext=1
                  CALL infoargext(argoutybasref,'argoutybasref',flagext,jmod)
               ENDIF
            CASE (47)
C --- -outbiasdta
               IF (largoutbiasdta) THEN
                  flagext=3
                  CALL infoargext(argoutbiasdta,'argoutbiasdta',flagext,jmod)
               ENDIF
            CASE (48)
C --- -affectobs
               IF (largaffectobs) THEN
                  WRITE(numout,*) ' observation variable (affectobs) = ',
     $                 argaffectobs(1:lenv(argaffectobs))
              ENDIF
           CASE (49)
C --- -nullobs
               IF (largnullobs) THEN
                  WRITE(numout,*) ' observation variable (nullobs) = ',
     $                 argnullobs(1:lenv(argnullobs))
              ENDIF
            CASE (50)
C --- -inpartobs (option 2)
               IF (larginpartobs) THEN
                  WRITE(numout,*) ' obs partition file (inpartobs) = ',
     $                 arginpartobs(1:lenv(arginpartobs))
               ENDIF
            CASE (51)
C --- -argoutdiaghst (option 2)
               IF (largoutdiaghst) THEN
                  WRITE(numout,*) ' output diag file (outdiaghst) = ',
     $                 argoutdiaghst(1:lenv(argoutdiaghst))
               ENDIF
            CASE (52)
C --- -configobs
               IF (largconfigobs) THEN
                  flagext=4
                  CALL infoargext(argconfigobs,'argconfigobs',flagext,jmod)
               ENDIF
            CASE (53)
C --- -fixjpx (option 1)
               IF (largfixjpx) THEN
                  WRITE(numout,*) ' Vx block size (fixjpx) = ',
     $                 argfixjpx(1:lenv(argfixjpx))
               ENDIF
            CASE (54)
C --- -reducevar
               IF (largreducevar) THEN
                  flagext=5
                  CALL infoargext(argreducevar,'argreducevar',flagext,jmod)
                  IF (extvarmem(indext(argreducevar,extvartab,nbextvar))/10
     $                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (55)
C --- -scale (option 2)
               IF (largscale) THEN
                  WRITE(numout,*) ' observation scaling factor (scale) = ',
     $                 argscale(1:lenv(argscale))
               ENDIF
            CASE (56)
C --- -biasdbs (option 2)
               IF (largbiasdbs) THEN
                  flagext=2
                  CALL infoargext(argbiasdbs,'argbiasdbs',flagext,jmod)
               ENDIF
            CASE (57)
C --- -inpartvar
               IF (larginpartvar) THEN
                  flagext=5
                  CALL infoargext(arginpartvar,'arginpartvar',flagext,jmod)
                  IF (extvarmem(indext(arginpartvar,extvartab,nbextvar))/10
     $                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (58)
C --- -outpartvar
               IF (largoutpartvar) THEN
                  flagext=5
                  CALL infoargext(argoutpartvar,'argoutpartvar',flagext,jmod)
                  IF (extvarmem(indext(argoutpartvar,extvartab,nbextvar))/10
     $                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (59)
C --- -inzon
               IF (larginzon) THEN
                  flagext=6
                  CALL infoargext(arginzon,'arginzon',flagext,jmod)
               ENDIF
            CASE (60)
C --- -outzon
               IF (largoutzon) THEN
                  flagext=6
                  CALL infoargext(argoutzon,'argoutzon',flagext,jmod)
               ENDIF
            CASE (61)
C --- -oecorrel
               IF (largoecorrel) THEN
                  WRITE(numout,*) ' observation error correlation matrix = ',
     $                 argoecorrel(1:lenv(argoecorrel))
               ENDIF
            CASE (62)
C --- -fecorrel
               IF (largfecorrel) THEN
                  WRITE(numout,*) ' forecast error correlation matrix = ',
     $                 argfecorrel(1:lenv(argfecorrel))
               ENDIF
            CASE (63)
C --- -zonindex
               IF (largzonindex) THEN
                  WRITE(numout,*) ' index in local data section (zonindex) = ',
     $                 argzonindex(1:lenv(argzonindex))
               ENDIF
            CASE (64)
C --- -incfg
               IF (largincfg) THEN
                  WRITE(numout,*) ' input configuration file (incfg) = ',
     $                 argincfg(1:lenv(argincfg))
               ENDIF
            CASE (65)
C --- -inybas
               IF (larginybas) THEN
                  flagext=1
                  CALL infoargext(arginybas,'arginybas',flagext,jmod)
               ENDIF
            CASE (66)
C --- -outybas
               IF (largoutybas) THEN
                  flagext=1
                  CALL infoargext(argoutybas,'argoutybas',flagext,jmod)
               ENDIF
            CASE (67)
C --- -inobas
               IF (larginobas) THEN
                  flagext=1
                  CALL infoargext(arginobas,'arginobas',flagext,jmod)
               ENDIF
            CASE (68)
C --- -outobas
               IF (largoutobas) THEN
                  flagext=1
                  CALL infoargext(argoutobas,'argoutobas',flagext,jmod)
               ENDIF
            CASE (69)
C --- -inzbas
               IF (larginzbas) THEN
                  flagext=1
                  CALL infoargext(arginzbas,'arginzbas',flagext,jmod)
               ENDIF
            CASE (70)
C --- -outzbas
               IF (largoutzbas) THEN
                  flagext=1
                  CALL infoargext(argoutzbas,'argoutzbas',flagext,jmod)
               ENDIF
            CASE (71)
C --- -configzon
               IF (largconfigzon) THEN
                  flagext=6
                  CALL infoargext(argconfigzon,'argconfigzon',flagext,jmod)
               ENDIF
            CASE (72)
C --- -coefrmax
               IF (largcoefrmax) THEN
                  WRITE(numout,*) ' maximum ROA coefficient (coefrmax) = ',
     $                 argcoefrmax(1:lenv(argcoefrmax))
               ENDIF
            CASE (73)
C --- -disable
               IF (largdisable) THEN
                  WRITE(numout,*) ' disable option = ',
     $                 argdisable(1:lenv(argdisable))
               ENDIF
            CASE (74)
C --- -inobasref
               IF (larginobasref) THEN
                  flagext=1
                  CALL infoargext(arginobasref,'arginobasref',flagext,jmod)
               ENDIF
            CASE (75)
C --- -fixjpz (option 1)
               IF (largfixjpz) THEN
                  WRITE(numout,*) ' maximum number of local ',
     $                 'data sections in memory (fixjpz) = ',
     $                 argfixjpz(1:lenv(argfixjpz))
               ENDIF
            CASE (76)
C --- -inparadap (option 2)
               IF (larginparadap) THEN
                  WRITE(numout,*) ' input adaptive parameter file  = ',
     $                 arginparadap(1:lenv(arginparadap))
              ENDIF
            CASE (77)
C --- -fixjpu (option 1)
               IF (largfixjpu) THEN
                  WRITE(numout,*) ' initial local observation vector ',
     $                 'size (fixjpu)  = ',argfixjpu(1:lenv(argfixjpu))
               ENDIF
            CASE (78)
C --- -inptzon
               IF (larginptzon) THEN
                  flagext=6
                  CALL infoargext(arginptzon,'arginptzon',flagext,jmod)
               ENDIF
            CASE (79)
C --- -outptzon
               IF (largoutptzon) THEN
                  flagext=6
                  CALL infoargext(argoutptzon,'argoutptzon',flagext,jmod)
               ENDIF
            CASE (80)
C --- -action (option 1)
               IF (largaction) THEN
                  WRITE(numout,*) ' action number = ',
     $                 argaction(1:lenv(argaction))
               ENDIF
            CASE (81)
C --- -inzbasref
               IF (larginzbasref) THEN
                  flagext=1
                  CALL infoargext(arginzbasref,'arginzbasref',flagext,jmod)
               ENDIF
            CASE (82)
C --- -inzonref
               IF (larginzonref) THEN
                  flagext=6
                  CALL infoargext(arginzonref,'arginzonref',flagext,jmod)
               ENDIF
            CASE (83)
C --- -inoptcfg
               IF (larginoptcfg) THEN
                  WRITE(numout,*) ' optional config file (inoptcfg) = ',
     $                 arginoptcfg(1:lenv(arginoptcfg))
               ENDIF
            CASE (84)
C --- -connect
               IF (largconnect) THEN
                  WRITE(numout,*) ' connection file = ',
     $                 argconnect(1:lenv(argconnect))
               ENDIF
            CASE (85)
C --- -incstr
               IF (largincstr) THEN
                  WRITE(numout,*) ' constraint directory = ',
     $                 argincstr(1:lenv(argincstr))
               ENDIF
            CASE (87)
C --- -insmocfg
               IF (larginsmocfg) THEN
                  WRITE(numout,*) ' optional config file (insmocfg) = ',
     $                 arginsmocfg(1:lenv(arginsmocfg))
               ENDIF
            CASE (88)
C --- -outsmocfg
               IF (largoutsmocfg) THEN
                  WRITE(numout,*) ' optional config file (outsmocfg) = ',
     $                 argoutsmocfg(1:lenv(argoutsmocfg))
               ENDIF
            CASE (89)
C --- -inrz
               IF (larginrz) THEN
                  WRITE(numout,*) ' reduced space file  = ',
     $                 arginrz(1:lenv(arginrz))
              ENDIF
            END SELECT  
         ENDDO
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(jmod,1000,1,'readarg','infoarg')
C
 101  WRITE (texterror,*) 'Invalid mask extensiona'
      CALL printerror2(jmod,101,3,'readarg','infoarg',comment=texterror)
 102  WRITE (texterror,*) 'Cannot read Vx vector block by block',
     $         ' with such file format'
      CALL printerror2(jmod,102,3,'readarg','infoarg',comment=texterror)
 103  WRITE (texterror,*) 'Cannot read Cx covariance block by block',
     $         ' with such file format'
      CALL printerror2(jmod,103,3,'readarg','infoarg',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE infoargext (argument,text,flagext,jmod)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Print information about argument of given extension
CCC  -------
CCC  Method : 
CCC  ------   
CCC  Input : argument, argument name, type of argument, module
CCC  -----
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: argument,text
      INTEGER, intent(in) :: flagext,jmod
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jext,xypos
      INTEGER :: jvar,indvar,jdta,inddta
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         SELECT CASE (flagext)
         CASE (1)
C --- bas
            IF ( (.NOT.(validextvarbas(argument)))
     $        .AND.(.NOT.(validextdtabas(argument)))
     $        .AND.(.NOT.(validextobsbas(argument)))
     $        .AND.(.NOT.(validextzonbas(argument)))
     $         ) GOTO 101
            jext=indext(argument,extbastab,nbextbas)
            WRITE(numout,*) ' ',text(1:lenv(text)),' ='
            WRITE(numout,*) '      argument  = ',argument(1:lenv(argument))
            WRITE(numout,*) '      jextbas   = ',jext
            WRITE(numout,*) '      extension = ',extbastab(jext)
     $           (1:lenv(extbastab(jext)))
         CASE (2)
C --- dbs
            IF (.NOT.(validextdbs(argument))) GOTO 101
            jext=indext(argument,extdbstab,nbextdbs)
            WRITE(numout,*) ' ',text(1:lenv(text)),' ='
            WRITE(numout,*) '      argument  = ',argument(1:lenv(argument))
            WRITE(numout,*) '      jextdbs   = ',jext
            WRITE(numout,*) '      extension = ',extdbstab(jext)
     $           (1:lenv(extdbstab(jext)))
         CASE (3)
C --- dta
            IF (.NOT.(validextdta(argument))) GOTO 101
            jext=indext(argument,extdtatab,nbextdta)
            WRITE(numout,*) ' ',text(1:lenv(text)),' ='
            IF (.NOT.(extdtaunit(jext))) THEN
               xypos=posit(argument,etoile)
               DO jdta=1,dtaend
                  inddta=dta_ord(jdta)
                  WRITE (clname,'(6x,a,a,a)') argument(1:(xypos-1)),
     $                 dtainam(inddta)(1:lenv(dtainam(inddta))),
     $                 argument((xypos+1):lenv(argument))
                  IF (jdta.EQ.1) THEN
                     WRITE(numout,*) '      argument(',jdta,'<=>'
     $                 ,dta_nam(inddta)(1:lenv(dta_nam(inddta)))
     $                 ,')  = ',clname(1:lenv(clname))
                  ELSE
                     WRITE(numout,*) '         ""   (',jdta,'<=>'
     $                 ,dta_nam(inddta)(1:lenv(dta_nam(inddta)))
     $                 ,')  = ',clname(1:lenv(clname))
                  ENDIF
               ENDDO
            ELSE
               WRITE(numout,*) '      argument  = ',
     $              argument(1:lenv(argument))
            ENDIF
            WRITE(numout,*) '      jextdta   = ',jext
            WRITE(numout,*) '      extension = ',extdtatab(jext)
     $           (1:lenv(extdtatab(jext)))
         CASE (4)
C --- obs
            IF (.NOT.(validextobs(argument))) GOTO 101
            jext=indext(argument,extobstab,nbextobs)
            WRITE(numout,*) ' ',text(1:lenv(text)),' ='
            WRITE(numout,*) '      argument  = ',
     $              argument(1:lenv(argument))
            WRITE(numout,*) '      jextobs   = ',jext
            WRITE(numout,*) '      extension = ',extobstab(jext)
     $           (1:lenv(extobstab(jext)))
         CASE (5)
C --- var
            IF (.NOT.(validextvar(argument))) GOTO 101
            jext=indext(argument,extvartab,nbextvar)
            WRITE(numout,*) ' ',text(1:lenv(text)),' ='
            IF (.NOT.(extvarunit(jext))) THEN
               xypos=posit(argument,etoile)
               DO jvar=1,varend
                  indvar=var_ord(jvar)
                  WRITE (clname,'(6x,a,a,a)') argument(1:(xypos-1)),
     $                 varinam(indvar)(1:lenv(varinam(indvar))),
     $                 argument((xypos+1):lenv(argument))
                  IF (jvar.EQ.1) THEN
                     WRITE(numout,*) '      argument(',jvar,'<=>'
     $                 ,var_nam(indvar)(1:lenv(var_nam(indvar)))
     $                 ,')  = ',clname(1:lenv(clname))
                  ELSE
                     WRITE(numout,*) '         ""   (',jvar,'<=>'
     $                 ,var_nam(indvar)(1:lenv(var_nam(indvar)))
     $                 ,')  = ',clname(1:lenv(clname))
                  ENDIF
               ENDDO
            ELSE
               WRITE(numout,*) '      argument  = ',
     $              argument(1:lenv(argument))
            ENDIF
            WRITE(numout,*) '      jextvar   = ',jext
            WRITE(numout,*) '      extension = ',extvartab(jext)
     $           (1:lenv(extvartab(jext)))
         CASE (6)
C --- zon
            IF (.NOT.(validextzon(argument))) GOTO 101
            jext=indext(argument,extzontab,nbextzon)
            WRITE(numout,*) ' ',text(1:lenv(text)),' ='
            IF (.NOT.(extzonunit(jext))) THEN
               xypos=posit(argument,etoile)
               DO jdta=1,dtaend
                  inddta=dta_ord(jdta)
                  WRITE (clname,'(6x,a,a,a)') argument(1:(xypos-1)),
     $                 dtainam(inddta)(1:lenv(dtainam(inddta))),
     $                 argument((xypos+1):lenv(argument))
                  IF (jdta.EQ.1) THEN
                     WRITE(numout,*) '      argument(',jdta,'<=>'
     $                 ,dta_nam(inddta)(1:lenv(dta_nam(inddta)))
     $                 ,')  = ',clname(1:lenv(clname))
                  ELSE
                     WRITE(numout,*) '         ""   (',jdta,'<=>'
     $                 ,dta_nam(inddta)(1:lenv(dta_nam(inddta)))
     $                 ,')  = ',clname(1:lenv(clname))
                  ENDIF
               ENDDO
            ELSE
               WRITE(numout,*) '      argument  = ',
     $              argument(1:lenv(argument))
            ENDIF
            WRITE(numout,*) '      jextzon   = ',jext
            WRITE(numout,*) '      extension = ',extzontab(jext)
     $           (1:lenv(extzontab(jext)))
C
         CASE DEFAULT
            GOTO 1000
         END SELECT  
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(jmod,1000,1,'readarg','infoargext')
C
 101  WRITE (texterror,*) 'argument ',argument(1:lenv(argument)),
     $     ' is not valid'
      CALL printerror2(jmod,101,3,'readarg','infoargext',
     $      comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilarg
