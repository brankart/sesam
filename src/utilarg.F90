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
! ---                    UTILARG.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12 (C.E. Testut)                        ---
! --- modification : 99-05 (C.E. Testut)                        ---
! --- modification : 03-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE affectarg : Check argument validity
! --- SUBROUTINE affectargext : Check file extension validity
! --- SUBROUTINE infoarg : Print information about current action
! --- SUBROUTINE infoargext : Print information about arguments
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilarg
      use mod_main
      use utilvalid
      IMPLICIT NONE
      PRIVATE

      PUBLIC affectarg,affectargext,infoarg,infoargext

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE affectarg (wrd,jmod,jcas,swibool1,modbool1,extbasbool1, &
     &     extdbsbool1,extdtabool1,extobsbool1,extvarbool1,extzonbool1)
!---------------------------------------------------------------------
!
!  Purpose : Check argument validity
!  -------
!  Method : To each switch correspond particular argument checkings
!  ------   Specific error messages are displayed
!
!  Input : module, switch, argument, and tables of validity
!  -----
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: wrd
      INTEGER, intent(in) :: jmod,jcas
      LOGICAL, dimension(:), intent(in) :: swibool1, modbool1
      LOGICAL, dimension(0:), intent(in) ::  &
     &     extbasbool1, extdbsbool1, extdtabool1, &
     &     extobsbool1, extvarbool1, extzonbool1
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: extension
      INTEGER :: jhelp,jhelp1,jext,ltext
      INTEGER :: jvar,indvar,jobs,indobs,inddbs,flagext
      LOGICAL :: autorise,inbool
!----------------------------------------------------------------------
!
      SELECT CASE (jcas)
      CASE (1)
! ---  -mode
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largmod) GOTO 102
         IF (.NOT.(validmod(wrd))) GOTO 103
         IF (wrd.NE.modtab(jmod)) GOTO 103
         IF (.NOT.(modbool1(jmod))) GOTO 1000
         argmod = wrd
         largmod=.TRUE.
      CASE (2)
! --- -help (option)
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
! --- -list (option)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (larglist) GOTO 102
         arglist = wrd
         larglist=.TRUE.
      CASE (4)
! --- -varmsk (option)
         IF (wrd.NE.'nomask') THEN
            flagext=5
            inbool=.TRUE.
            CALL affectargext (wrd,jmod,argvarmsk,largvarmsk, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
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
! --- -dtamsk (option)
         IF (wrd.NE.'nomask') THEN
            flagext=3
            inbool=.TRUE.
            CALL affectargext (wrd,jmod,argdtamsk,largdtamsk, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
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
! --- -weight (option)
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argweight,largweight, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.( &
     &        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
      CASE (7)
! --- -oestd (option)
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argoestd,largoestd, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.( &
     &        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
      CASE (8)
! --- -bias (option)
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argbias,largbias, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.( &
     &        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
      CASE (9)
! --- -outinfo (option)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largoutinfo) GOTO 102
         argoutinfo = wrd
         largoutinfo=.TRUE.
      CASE (10)
! --- -inxbas
         flagext=1
         IF (.NOT.validextvarbas(wrd)) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginxbas,larginxbas, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existbas=.TRUE.
      CASE (11)
! --- -indbs
         flagext=2
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argindbs,largindbs, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.2).AND.(.NOT.( &
     &        extdbsbool1(indext(wrd,extdbstab,nbextdbs))))) GOTO 106
         existdbs=.TRUE.
      CASE (12)
! --- -indta
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argindta,largindta, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (13)
! --- -inobs
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginobs,larginobs, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.( &
     &        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existobs=.TRUE.
      CASE (14)
! --- -invar
         flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginvar,larginvar, &
     &        flagext, switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (15)
! --- -inxbasref
         flagext=1
         IF (.NOT.validextvarbas(wrd)) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginxbasref,larginxbasref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existbas=.TRUE.
      CASE (16)
! --- -indbsref
         flagext=2
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argindbsref,largindbsref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.2).AND.(.NOT.( &
     &        extdbsbool1(indext(wrd,extdbstab,nbextdbs))))) GOTO 106
         existdbs=.TRUE.
      CASE (17)
! --- -indtaref
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argindtaref,largindtaref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (18)
! --- -inobsref
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginobsref,larginobsref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.( &
     &        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existobs=.TRUE.
      CASE (19)
! --- -invarref
         flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginvarref,larginvarref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (20)
! --- -outxbas
         flagext=1
         IF (.NOT.validextvarbas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutxbas,largoutxbas, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existbas=.TRUE.
      CASE (21)
! --- -outdta
         flagext=3
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutdta,largoutdta, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         existdta=.TRUE.
      CASE (22)
! --- -outobs
         flagext=4
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutobs,largoutobs, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.4).AND.(.NOT.( &
     &        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
!         existobs=.TRUE.     
      CASE (23)
! --- -outvar
         flagext=5
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutvar,largoutvar, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (24)
! --- -outxbasref
         flagext=1
         IF (.NOT.validextvarbas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutxbasref,largoutxbasref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existbas=.TRUE.
      CASE (25)
! --- -outdtaref
         flagext=3
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutdtaref,largoutdtaref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         existdta=.TRUE.
      CASE (26)
! --- -outobsref
         flagext=4
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutobsref,largoutobsref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.4).AND.(.NOT.( &
     &        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
!         existobs=.TRUE.     
      CASE (27)
! --- -outvarref
         flagext=5
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutvarref,largoutvarref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (30)
! --- -reducedta
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argreducedta,largreducedta, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (31)
! --- -outparadap
         argoutparadap=wrd
         largoutparadap=.TRUE.
      CASE (32)
! --- -outrz
!        jext=indext(wrd,extdtatab,nbextdta)
!        IF (jext.EQ.0) GOTO 106
!        IF (.NOT.(extdtabool1(jext))) GOTO 106
!        autorise=((jext.EQ.1).OR.(jext.EQ.2).OR.(jext.EQ.3))
!        IF (.NOT.(autorise)) GOTO 106
         argoutrz=wrd
         largoutrz=.TRUE.
      CASE (33)
! --- -typeoper
         argtypeoper=wrd
         largtypeoper=.TRUE.
      CASE (34)
! --- -typedtadiag
         jext=indext(wrd,extdtatab,nbextdta)
         IF (jext.EQ.0) GOTO 106
         IF (.NOT.(extdtabool1(jext))) GOTO 106
         autorise=((jext.EQ.1).OR.(jext.EQ.2).OR.(jext.EQ.3))
         IF (.NOT.(autorise)) GOTO 106
         argtypedtadiag=wrd
         largtypedtadiag=.TRUE.
      CASE (35)
! --- -diffobsref
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argdiffobsref,largdiffobsref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.( &
     &        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
      CASE (36)
! --- -diffdtaref
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argdiffdtaref,largdiffdtaref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (37)
! --- -diffvarref
         flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argdiffvarref,largdiffvarref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (38)
! --- -diffobsorg
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argdiffobsorg,largdiffobsorg, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.( &
     &        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
      CASE (39)
! --- -diffdtaorg
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argdiffdtaorg,largdiffdtaorg, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (40)
! --- -diffvarorg
         flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argdiffvarorg,largdiffvarorg, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (41)
! --- -inerrdta
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginerrdta,larginerrdta, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (42)
! --- -outerrdta
         flagext=3
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argouterrdta,largouterrdta, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         existdta=.TRUE.
      CASE (43)
! --- -instddta
         flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginstddta,larginstddta, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existdta=.TRUE.
      CASE (44)
! --- -outstddta
         flagext=3
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutstddta,largoutstddta, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         existdta=.TRUE.
      CASE (45)
! --- -inybasref
         flagext=1
         IF ((.NOT.validextvarbas(wrd)) &
     &        .AND.(.NOT.validextdtabas(wrd))) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginybasref,larginybasref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existybas=.TRUE.
      CASE (46)
! --- -outybasref
         flagext=1
         IF (.NOT.validextdtabas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutybasref,largoutybasref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existybas=.TRUE.
      CASE (47)
! --- -outbiasdta
         flagext=3
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutbiasdta,largoutbiasdta, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         existdta=.TRUE.
      CASE (48)
! --- -affectobs
         IF (.NOT.(swibool1(jcas))) GOTO 101
         LOOPAFFECT : DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            ltext=lenv(obs_nam(indobs,inddbs))
            IF ((wrd(1:ltext).EQ.obs_nam(indobs,inddbs)(1:ltext)) &
     &              .AND.((wrd((ltext+1):(ltext+1))).EQ.(' '))) THEN
               argaffectobs=wrd 
               largaffectobs= .TRUE.
               EXIT LOOPAFFECT
            ENDIF            
         ENDDO LOOPAFFECT
         IF (.NOT.largaffectobs) GOTO 103
      CASE (49)
! --- -nullobs
         IF (.NOT.(swibool1(jcas))) GOTO 101
         LOOPNULL : DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            ltext=lenv(obs_nam(indobs,inddbs))
            IF (((wrd(1:ltext).EQ.obs_nam(indobs,inddbs)(1:ltext)) &
     &           .AND.((wrd((ltext+1):(ltext+1))).EQ.(' '))) &
     &           .OR.(wrd(1:4).EQ.'ALL ')) THEN
               argnullobs=wrd 
               largnullobs= .TRUE.
               EXIT LOOPNULL
            ENDIF            
         ENDDO LOOPNULL
         IF (.NOT.(largnullobs)) GOTO 103
      CASE (50)
! --- -inpartobs
         flagext=4
         IF (validextdta(wrd)) flagext=3
         IF (validextvar(wrd)) flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginpartobs,larginpartobs, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.3).AND.(.NOT.( &
     &        extdtabool1(indext(wrd,extdtatab,nbextdta))))) GOTO 106
         IF ((flagext.EQ.4).AND.(.NOT.( &
     &        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         IF (flagext.EQ.4) existobs=.TRUE.
      CASE (51)
! --- -anamorphosis
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (larganamorphosis) GOTO 102
         arganamorphosis = wrd
         larganamorphosis=.TRUE.
      CASE (52)
! --- -configobs
         flagext=4
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argconfigobs,largconfigobs, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.4).AND.(.NOT.( &
     &        extobsbool1(indext(wrd,extobstab,nbextobs))))) GOTO 106
         existobs=.TRUE.     
      CASE (53)
! --- -fixjpx (option 1)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largfixjpx) GOTO 102
         argfixjpx = wrd
         largfixjpx=.TRUE.
      CASE (54)
! --- -reducevar
         flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argreducevar,largreducevar, &
     &        flagext, switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (55)
! --- -scale
         IF (.NOT.(swibool1(jcas))) GOTO 101
         argscale = wrd
         largscale=.TRUE.
      CASE (56)
! --- -biasdbs
         flagext=2
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argbiasdbs,largbiasdbs, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.2).AND.(.NOT.( &
     &        extdbsbool1(indext(wrd,extdbstab,nbextdbs))))) GOTO 106
         IF (indext(argbiasdbs,extdbstab,nbextdbs).NE.2) GOTO 106
         existdbs=.TRUE.
      CASE (57)
! --- -inpartvar
         flagext=5
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginpartvar,larginpartvar, &
     &        flagext, switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
         existvar=.TRUE.
      CASE (58)
! --- -outpartvar
         flagext=5
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutpartvar,largoutpartvar, &
     &        flagext, switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.5).AND.(.NOT.( &
     &        extvarbool1(indext(wrd,extvartab,nbextvar))))) GOTO 106
      CASE (59)
! --- -inzon
         flagext=6
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginzon,larginzon, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.6).AND.(.NOT.( &
     &        extzonbool1(indext(wrd,extzontab,nbextzon))))) GOTO 106
         existzon=.TRUE.
      CASE (60)
! --- -outzon
         flagext=6
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutzon,largoutzon, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.6).AND.(.NOT.( &
     &        extzonbool1(indext(wrd,extzontab,nbextzon))))) GOTO 106
      CASE (61)
! --- -oecorrel
         flagext=1
         IF (.NOT.validextzonbas(wrd)) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argoecorrel,largoecorrel, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existzbas=.TRUE.
      CASE (62)
! --- -fecorrel
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largfecorrel) GOTO 102
         argfecorrel = wrd
         largfecorrel=.TRUE.
      CASE (63)
! --- -zonindex
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largzonindex) GOTO 102
         argzonindex = wrd
         largzonindex=.TRUE.
      CASE (64)
! --- -incfg
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largincfg) GOTO 102
         argincfg = wrd
         largincfg=.TRUE.
      CASE (65)
! --- -inybas
         flagext=1
         IF ((.NOT.validextvarbas(wrd)) &
     &        .AND.(.NOT.validextdtabas(wrd))) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginybas,larginybas, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existybas=.TRUE.
      CASE (66)
! --- -ouytbas
         flagext=1
         IF (.NOT.validextdtabas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutybas,largoutybas, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existybas=.TRUE.
      CASE (67)
! --- -inobas
         flagext=1
         IF ((.NOT.validextvarbas(wrd)) &
     &        .AND.(.NOT.validextdtabas(wrd)) &
     &        .AND.(.NOT.validextobsbas(wrd))) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginobas,larginobas, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existobas=.TRUE.
      CASE (68)
! --- -outobas
         flagext=1
         IF (.NOT.validextobsbas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutobas,largoutobas, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existobas=.TRUE.
      CASE (69)
! --- -inzbas
         flagext=1
         IF (.NOT.validextzonbas(wrd)) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginzbas,larginzbas, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existzbas=.TRUE.
      CASE (70)
! --- -outzbas
         flagext=1
         IF (.NOT.validextzonbas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutzbas,largoutzbas, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existzbas=.TRUE.
      CASE (71)
! --- -configzon
         flagext=6
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argconfigzon,largconfigzon, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.6).AND.(.NOT.( &
     &        extzonbool1(indext(wrd,extzontab,nbextzon))))) GOTO 106
         existzon=.TRUE.     
      CASE (72)
! --- -coefrmax
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largcoefrmax) GOTO 102
         argcoefrmax = wrd
         largcoefrmax=.TRUE.
      CASE (73)
! --- -disable
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largdisable) GOTO 102
         argdisable = wrd
         largdisable=.TRUE.
      CASE (74)
! --- -inobasref
         flagext=1
         IF ((.NOT.validextvarbas(wrd)) &
     &        .AND.(.NOT.validextdtabas(wrd)) &
     &        .AND.(.NOT.validextobsbas(wrd))) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginobasref,larginobasref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existobas=.TRUE.
      CASE (75)
! --- -fixjpz (option 1)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largfixjpz) GOTO 102
         argfixjpz = wrd
         largfixjpz=.TRUE.
      CASE (76)
! --- -inparadap (option 2)
         arginparadap = wrd
         larginparadap=.TRUE.
      CASE (77)
! --- -fixjpu (option 1)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largfixjpu) GOTO 102
         argfixjpu = wrd
         largfixjpu=.TRUE.
      CASE (78)
! --- -inptzon
         flagext=6
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginptzon,larginptzon, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.6).AND.(.NOT.( &
     &        extzonbool1(indext(wrd,extzontab,nbextzon))))) GOTO 106
         existzon=.TRUE.
      CASE (79)
! --- -outptzon
         flagext=6
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutptzon,largoutptzon, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.6).AND.(.NOT.( &
     &        extzonbool1(indext(wrd,extzontab,nbextzon))))) GOTO 106
      CASE (80)
! --- -action (option 1)
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largaction) GOTO 102
         argaction = wrd
         largaction=.TRUE.
      CASE (81)
! --- -inzbasref
         flagext=1
         IF (.NOT.validextzonbas(wrd)) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginzbasref,larginzbasref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existzbas=.TRUE.
      CASE (82)
! --- -inzonref
         flagext=6
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,arginzonref,larginzonref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.6).AND.(.NOT.( &
     &        extzonbool1(indext(wrd,extzontab,nbextzon))))) GOTO 106
         existzon=.TRUE.
      CASE (83)
! --- -inoptcfg
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (larginoptcfg) GOTO 102
         arginoptcfg = wrd
         larginoptcfg=.TRUE.
      CASE (84)
! --- -connect
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largconnect) GOTO 102
         argconnect = wrd
         largconnect=.TRUE.
      CASE (85)
! --- -incstr
         flagext=1
         IF (.NOT.validextvarbas(wrd)) GOTO 106
         inbool=.TRUE.
         CALL affectargext (wrd,jmod,argincstr,largincstr, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existcstr=.TRUE.
      CASE (86)
! --- -outobasref
         flagext=1
         IF (.NOT.validextobsbas(wrd)) GOTO 106
         inbool=.FALSE.
         CALL affectargext (wrd,jmod,argoutobasref,largoutobasref, &
     &        flagext,switab(jcas),swibool1(jcas),inbool)
         IF ((flagext.EQ.1).AND.(.NOT.( &
     &        extbasbool1(indext(wrd,extbastab,nbextbas))))) GOTO 106
         existobas=.TRUE.
      CASE (87)
! --- -insmocfg
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (larginsmocfg) GOTO 102
         arginsmocfg = wrd
         larginsmocfg=.TRUE.
      CASE (88)
! --- -outsmocfg
         IF (.NOT.(swibool1(jcas))) GOTO 101
         IF (largoutsmocfg) GOTO 102
         argoutsmocfg = wrd
         largoutsmocfg=.TRUE.
      CASE (89)
! --- -inrz
         arginrz=wrd
         larginrz=.TRUE.
      CASE (90)
! --- -iterate
         IF (largiterate) GOTO 102
         argiterate = wrd
         largiterate = .TRUE.
      CASE DEFAULT
         GOTO 1000
      END SELECT  
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(jmod,1000,1,'readarg','affectarg')
!
 101   WRITE (texterror,*) 'Switch ', &
     &     switab(jcas)(1:lenv(switab(jcas))), &
     &     ' is not available for this module (see -help swi)'
      CALL printerror2(jmod,101,3,'readarg','affectarg', &
     &      comment=texterror)
 102   WRITE (texterror,*) 'The same switch ', &
     &     switab(jcas)(1:lenv(switab(jcas))), &
     &     ' cannot be used more than once'
      CALL printerror2(jmod,102,3,'readarg','affectarg', &
     &      comment=texterror)
 103   WRITE (texterror,*) 'Argument ', &
     &     switab(jcas)(1:lenv(switab(jcas))), &
     &     ' is not valid'
      CALL printerror2(jmod,103,3,'readarg','affectarg', &
     &      comment=texterror)
 106   WRITE (texterror,*) 'File extension ',extension(1:lenv(extension)), &
     &     ' is not available', &
     &     ' with switch ',switab(jcas)(1:lenv(switab(jcas))), &
     &     ' (see -help ext)'
      CALL printerror2(jmod,106,3,'readarg','affectarg', &
     &      comment=texterror)

      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE affectargext (wrd,jmod,argument,largument, &
     &     flagext,swi,lswi,inbool)
!---------------------------------------------------------------------
!
!  Purpose : Check validity of file extensions
!  -------
!  Method : To each object type correspond particular checkings
!  ------   Specific error messages are displayed
!
!  Input : module, switch, argument, object type (flagext),
!  ------  input or output object (inbool)
!
!  Output : activate argument in SESAM: argument,largument
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: wrd,swi
      CHARACTER(len=*), intent(out) :: argument
      INTEGER, intent(in) :: flagext,jmod
      LOGICAL, intent(inout) :: largument
      LOGICAL, intent(in) :: lswi,inbool
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: extension
      INTEGER :: jext
!----------------------------------------------------------------------
!
      IF (largument) GOTO 102
      SELECT CASE (flagext)
      CASE (1)
! --- bas
         IF ( (.NOT.(validextvarbas(wrd))) &
     &        .AND.(.NOT.(validextdtabas(wrd))) &
     &        .AND.(.NOT.(validextobsbas(wrd))) &
     &        .AND.(.NOT.(validextzonbas(wrd))) &
     &        ) GOTO 103
         jext=indext(wrd,extbastab,nbextbas)
         extension=extbastab(jext)
         IF (.NOT.(extbasbool(jext))) GOTO 104
         IF (((extbasmem(jext)/10).EQ.0).AND.(inbool)) GOTO 106
         IF (((MOD(extbasmem(jext),10)).EQ.0).AND.(.NOT.inbool)) GOTO 107
      CASE (2)
! --- dbs
         IF (.NOT.(validextdbs(wrd))) GOTO 103
         jext=indext(wrd,extdbstab,nbextdbs)
         extension=extdbstab(jext)
         IF (.NOT.(extdbsbool(jext))) GOTO 104
         IF (((extdbsmem(jext)/10).EQ.0).AND.(inbool)) GOTO 106
         IF (((MOD(extdbsmem(jext),10)).EQ.0).AND.(.NOT.inbool)) GOTO 107
      CASE (3)
! --- dta
         IF (.NOT.(validextdta(wrd))) GOTO 103
         jext=indext(wrd,extdtatab,nbextdta)
         extension=extdtatab(jext)
         IF (.NOT.(extdtabool(jext))) GOTO 104
         IF (((extdtamem(jext)/10).EQ.0).AND.(inbool)) GOTO 106
         IF (((MOD(extdtamem(jext),10)).EQ.0).AND.(.NOT.inbool)) GOTO 107
      CASE (4)
! --- obs
         IF (.NOT.(validextobs(wrd))) GOTO 103
         jext=indext(wrd,extobstab,nbextobs)
         extension=extobstab(jext)
         IF (.NOT.(extobsbool(jext))) GOTO 104
         IF (((extobsmem(jext)/10).EQ.0).AND.(inbool)) GOTO 106
         IF (((MOD(extobsmem(jext),10)).EQ.0).AND.(.NOT.inbool)) GOTO 107
      CASE (5)
! --- var
         IF (.NOT.(validextvar(wrd))) GOTO 103
         jext=indext(wrd,extvartab,nbextvar)
         extension=extvartab(jext)
         IF (.NOT.(extvarbool(jext))) GOTO 104
         IF (((extvarmem(jext)/10).EQ.0).AND.(inbool)) GOTO 106
         IF (((MOD(extvarmem(jext),10)).EQ.0).AND.(.NOT.inbool)) GOTO 107
      CASE (6)
! --- zon
         IF (.NOT.(validextzon(wrd))) GOTO 103
         jext=indext(wrd,extzontab,nbextzon)
         extension=extzontab(jext)
         IF (.NOT.(extzonbool(jext))) GOTO 104
         IF (((extzonmem(jext)/10).EQ.0).AND.(inbool)) GOTO 106
         IF (((MOD(extzonmem(jext),10)).EQ.0).AND.(.NOT.inbool)) GOTO 107
!
      CASE DEFAULT
         GOTO 1000
      END SELECT
      argument = wrd
      largument=.TRUE.
      RETURN
!
! --- error management
!
 1000 CALL printerror2(jmod,1000,1,'readarg','affectargext')
!
 102  WRITE (texterror,*) 'The same switch ', &
     &     swi(1:lenv(swi)), &
     &     ' cannot be used more than once'
      CALL printerror2(jmod,102,3,'readarg','affectargext', &
     &      comment=texterror)
 103   WRITE (texterror,*) 'Argument ', &
     &      wrd(1:lenv(wrd)), &
     &      ' is not valid'
      CALL printerror2(jmod,103,3,'readarg','affectargext', &
     &      comment=texterror)
 104  WRITE (texterror,*) 'Extension ',extension(1:lenv(extension)), &
     &     ' is not available with switch ',swi(1:lenv(swi)), &
     &     ' (see -help ext)'
      CALL printerror2(jmod,104,3,'readarg','affectargext', &
     &      comment=texterror)
 106  WRITE (texterror,*) 'Extension ',extension(1:lenv(extension)), &
     &     ' is not valid with switch ',swi(1:lenv(swi)), &
     &     ' for input files (see -help ext)' 
      CALL printerror2(jmod,106,3,'readarg','affectargext', &
     &      comment=texterror)
 107  WRITE (texterror,*) 'Extension ',extension(1:lenv(extension)), &
     &     ' is not valid with switch ',swi(1:lenv(swi)), &
     &     ' for output files (see -help ext)' 
      CALL printerror2(jmod,107,3,'readarg','affectargext', &
     &      comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE infoarg (jmod,jaction,swibool1,modbool1,extbasbool1, &
     &     extdbsbool1,extdtabool1,extobsbool1,extvarbool1,extzonbool1)
!---------------------------------------------------------------------
!
!  Purpose : Print information about current action
!  -------
!  Method : Loop over all user switches and write
!  ------   information about them
!
!  Input : mode, action and tables of usage
!  -----
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, intent(in) :: jmod,jaction
      LOGICAL, dimension(:), intent(in) :: swibool1, modbool1
      LOGICAL, dimension(0:), intent(in) ::  &
     &     extbasbool1, extdbsbool1, extdtabool1, &
     &     extobsbool1, extvarbool1, extzonbool1
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: argument
      INTEGER :: jcas
      INTEGER :: jext,xypos,flagext
      INTEGER :: jvar,indvar,jdta,inddta
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         DO jcas=1,nbarg
            jext = 0
            xypos= 0
            SELECT CASE (jcas)
            CASE (1)
! --- -mode
               IF (largmod) THEN
                  WRITE(numout,*) ' -----------------------------'
                  WRITE(numout,*) ' module = ', &
     &                 modtab(jmod)(1:lenv(modtab(jmod)))
                  WRITE(numout,*) ' action = ',jaction
                  WRITE(numout,*) ' -----------------------------'
              ENDIF
            CASE (2)
! --- -help
               IF (larghelp) THEN
                  WRITE(numout,*) ' -----------------------------'
                  WRITE(numout,*) ' help  = ', &
     &                 arghelp(1:lenv(arghelp))
                  WRITE(numout,*) ' -----------------------------'
               ENDIF
            CASE (3)
! --- -list (option)
               IF (larglist) THEN
                  WRITE(numout,*) ' namelist (option name) = ', &
     &                 arglist(1:lenv(arglist))
               ELSE
                  WRITE(numout,*) ' namelist (default name)  = ', &
     &                 fnamlist(1:lenv(fnamlist))
               ENDIF
            CASE (4)
! --- -varmsk (option)
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
                           WRITE (clname,'(6x,a,a,a)') argvarmsk(1:(xypos-1)), &
     &                       varinam(indvar)(1:lenv(varinam(indvar))), &
     &                       argvarmsk((xypos+1):lenv(argvarmsk))
                           IF (jvar.EQ.1) THEN
                              WRITE(numout,*) '      mask(',jvar,'<=>', &
     &                          var_nam(indvar)(1:lenv(var_nam(indvar))) &
     &                          ,')  = ',clname(1:lenv(clname))
                           ELSE
                              WRITE(numout,*) '       "" (',jvar,'<=>', &
     &                          var_nam(indvar)(1:lenv(var_nam(indvar))) &
     &                          ,')  = ',clname(1:lenv(clname))
                           ENDIF
                        ENDDO
                     ELSE
                        WRITE(numout,*) '      mask  = ', &
     &                    argvarmsk(1:lenv(argvarmsk))
                     ENDIF
                  ENDIF
               ELSE
                  WRITE(numout,*) ' Vx mask (default name) ='
                  DO jvar=1,varend
                     indvar=var_ord(jvar)
                     IF (jvar.EQ.1) THEN
                        WRITE(numout,*) '      mask(',jvar,'<=>', &
     &                       var_nam(indvar)(1:lenv(var_nam(indvar))), &
     &                       ')  = ', &
     &                       varfmsk(indvar)(1:lenv(varfmsk(indvar)))
                     ELSE
                        WRITE(numout,*) '       "" (',jvar,'<=>', &
     &                       var_nam(indvar)(1:lenv(var_nam(indvar))), &
     &                       ')  = ', &
     &                       varfmsk(indvar)(1:lenv(varfmsk(indvar)))
                        ENDIF
                     ENDDO
               ENDIF
            CASE (5)
! --- -dtamsk (option)
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
                           WRITE (clname,'(6x,a,a,a)') argdtamsk(1:(xypos-1)), &
     &                       dtainam(inddta)(1:lenv(dtainam(inddta))), &
     &                       argdtamsk((xypos+1):lenv(argdtamsk))
                           IF (jdta.EQ.1) THEN
                              WRITE(numout,*) '      mask(',jdta,'<=>', &
     &                          dta_nam(inddta)(1:lenv(dta_nam(inddta))) &
     &                          ,')  = ',clname(1:lenv(clname))
                           ELSE
                              WRITE(numout,*) '       "" (',jdta,'<=>', &
     &                          dta_nam(inddta)(1:lenv(dta_nam(inddta))) &
     &                          ,')  = ',clname(1:lenv(clname))
                           ENDIF
                        ENDDO
                     ELSE
                        WRITE(numout,*) '      mask  = ', &
     &                    argdtamsk(1:lenv(argdtamsk))
                     ENDIF
                  ENDIF
               ELSE
                  WRITE(numout,*) ' Vy mask (default name) ='
                  DO jdta=1,dtaend
                     inddta=dta_ord(jdta)
                     IF (jdta.EQ.1) THEN
                        WRITE(numout,*) '      mask(',jdta,'<=>', &
     &                       dta_nam(inddta)(1:lenv(dta_nam(inddta))), &
     &                       ')  = ', &
     &                       dtafmsk(inddta)(1:lenv(dtafmsk(inddta)))
                     ELSE
                        WRITE(numout,*) '       "" (',jdta,'<=>', &
     &                       dta_nam(inddta)(1:lenv(dta_nam(inddta))), &
     &                       ')  = ', &
     &                       dtafmsk(inddta)(1:lenv(dtafmsk(inddta)))
                        ENDIF
                     ENDDO
               ENDIF
            CASE (6)
! --- -weight (option)
               IF (largweight) THEN
                  flagext=5
                  IF (validextobs(argweight)) flagext=4
                  IF (validextdta(argweight)) flagext=3
                  CALL infoargext(argweight,'argweight',flagext,jmod)
                  IF (flagext.EQ.5) THEN
                    IF (extvarmem(indext(argweight,extvartab,nbextvar))/10 &
     &                   .GT.nallmem) GOTO 102
                  ENDIF
               ENDIF
            CASE (7)
! --- -oestd (option)
               IF (largoestd) THEN
                  flagext=5
                  IF (validextobs(argoestd)) flagext=4
                  IF (validextdta(argoestd)) flagext=3
                  CALL infoargext(argoestd,'argoestd',flagext,jmod)
                  IF (flagext.EQ.5) THEN
                    IF (extvarmem(indext(argoestd,extvartab,nbextvar))/10 &
     &                   .GT.nallmem) GOTO 102
                  ENDIF
               ENDIF
            CASE (8)
! --- -bias (option)
               IF (largbias) THEN
                  flagext=5
                  IF (validextobs(argbias)) flagext=4
                  IF (validextdta(argbias)) flagext=3
                  CALL infoargext(argbias,'argbias',flagext,jmod)
                  IF (flagext.EQ.5) THEN
                    IF (extvarmem(indext(argbias,extvartab,nbextvar))/10 &
     &                   .GT.nallmem) GOTO 102
                  ENDIF
               ENDIF
            CASE (9)
! --- -outinfo (option)
               IF (largoutinfo) THEN
                  WRITE(numout,*) ' listing file (option name) = ', &
     &                 argoutinfo(1:lenv(argoutinfo))
               ENDIF
            CASE (10)
! --- -inxbas
               IF (larginxbas) THEN
                  flagext=1
                  CALL infoargext(arginxbas,'arginxbas',flagext,jmod)
                  IF (extvarmem(indext(arginxbas(1:lenv(arginxbas)- &
     &                 lenv(extbastab(indext(arginxbas,extbastab,nbextbas)))), &
     &                 extvartab, nbextvar ))/10.GT.nallmem) GOTO 103
               ENDIF
            CASE (11)
! --- -indbs
               IF (largindbs) THEN
                  flagext=2
                  CALL infoargext(argindbs,'argindbs',flagext,jmod)
               ENDIF
            CASE (12)
! --- -indta
               IF (largindta) THEN
                  flagext=3
                  IF (validextvar(argindta)) flagext=5
                  CALL infoargext(argindta,'argindta',flagext,jmod)
               ENDIF
            CASE (13)
! --- -inobs
               IF (larginobs) THEN
                  flagext=4
                  IF (validextdta(arginobs)) flagext=3
                  IF (validextvar(arginobs)) flagext=5
                  CALL infoargext(arginobs,'arginobs',flagext,jmod)
               ENDIF
            CASE (14)
! --- -invar
               IF (larginvar) THEN
                  flagext=5
                  CALL infoargext(arginvar,'arginvar',flagext,jmod)
                  IF (extvarmem(indext(arginvar,extvartab,nbextvar))/10 &
     &                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (15)
! --- -inxbasref
               IF (larginxbasref) THEN
                  flagext=1
                  CALL infoargext(arginxbasref,'arginxbasref',flagext, &
     &                 jmod)
                  IF (extvarmem(indext(arginxbasref(1:lenv(arginxbasref)- &
     &                 lenv(extbastab(indext(arginxbasref,extbastab,nbextbas)))), &
     &                 extvartab, nbextvar ))/10.GT.nallmem) GOTO 103
               ENDIF
            CASE (16)
! --- -indbsref
               IF (largindbsref) THEN
                  flagext=2
                  CALL infoargext(argindbsref,'argindbsref',flagext,jmod)
               ENDIF
            CASE (17)
! --- -indtaref
               IF (largindtaref) THEN
                  flagext=3
                  IF (validextvar(argindtaref)) flagext=5
                  CALL infoargext(argindtaref,'argindtaref',flagext,jmod)
               ENDIF
            CASE (18)
! --- -inobsref
               IF (larginobsref) THEN
                  flagext=4
                  IF (validextdta(arginobsref)) flagext=3
                  IF (validextvar(arginobsref)) flagext=5
                  CALL infoargext(arginobsref,'arginobsref',flagext,jmod)
               ENDIF
            CASE (19)
! --- -invarref
               IF (larginvarref) THEN
                  flagext=5
                  CALL infoargext(arginvarref,'arginvarref',flagext,jmod)
                  IF (extvarmem(indext(arginvarref,extvartab,nbextvar))/10 &
     &                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (20)
! --- -outxbas
               IF (largoutxbas) THEN
                  flagext=1
                  CALL infoargext(argoutxbas,'argoutxbas',flagext,jmod)
                  IF (extvarmem(indext(argoutxbas(1:lenv(argoutxbas)- &
     &                 lenv(extbastab(indext(argoutxbas,extbastab,nbextbas)))), &
     &                 extvartab, nbextvar ))/10.GT.nallmem) GOTO 103
               ENDIF
            CASE (21)
! --- -outdta
               IF (largoutdta) THEN
                  flagext=3
                  CALL infoargext(argoutdta,'argoutdta',flagext,jmod)
               ENDIF
            CASE (22)
! --- -outobs
               IF (largoutobs) THEN
                  flagext=4
                  CALL infoargext(argoutobs,'argoutobs',flagext,jmod)
               ENDIF
            CASE (23)
! --- -outvar
               IF (largoutvar) THEN
                  flagext=5
                  CALL infoargext(argoutvar,'argoutvar',flagext,jmod)
                  IF (extvarmem(indext(argoutvar,extvartab,nbextvar))/10 &
     &                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (24)
! --- -outxbasref
               IF (largoutxbasref) THEN
                  flagext=1
                  CALL infoargext(argoutxbasref,'argoutxbasref',flagext,jmod)
                  IF (extvarmem(indext(argoutxbasref(1:lenv(argoutxbasref)- &
     &                 lenv(extbastab(indext(argoutxbasref,extbastab,nbextbas)))), &
     &                 extvartab, nbextvar ))/10.GT.nallmem) GOTO 103
               ENDIF
            CASE (25)
! --- -outdtaref
               IF (largoutdtaref) THEN
                  flagext=3
                  CALL infoargext(argoutdtaref,'argoutdtaref',flagext,jmod)
               ENDIF
            CASE (26)
! --- -outobsref
               IF (largoutobsref) THEN
                  flagext=4
                  CALL infoargext(argoutobsref,'argoutobsref',flagext,jmod)
               ENDIF
            CASE (27)
! --- -outvarref
               IF (largoutvarref) THEN
                  flagext=5
                  CALL infoargext(argoutvarref,'argoutvarref',flagext,jmod)
                  IF (extvarmem(indext(argoutvarref,extvartab,nbextvar))/10 &
     &                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (30)
! --- -reducedta
               IF (largreducedta) THEN
                  flagext=3
                  IF (validextvar(argreducedta)) flagext=5
                  CALL infoargext(argreducedta,'argreducedta', &
     &                 flagext,jmod)
               ENDIF
            CASE (31)
! --- -outparadap
              IF (largoutparadap) THEN
                  WRITE(numout,*) ' output adaptive parameter file  = ', &
     &                 argoutparadap(1:lenv(argoutparadap))
              ENDIF
            CASE (32)
! --- -outrz
              IF (largoutrz) THEN
                  WRITE(numout,*) ' reduced space file  = ', &
     &                 argoutrz(1:lenv(argoutrz))
              ENDIF
            CASE (33)
! --- -typeoper
              IF (largtypeoper) THEN
                  WRITE(numout,*) ' operation type = ', &
     &                 argtypeoper(1:lenv(argtypeoper))
              ENDIF
            CASE (34)
! --- -typedtadiag
              IF (largtypedtadiag) THEN
                  WRITE(numout,*) ' type of Vy diagnostic files  = ', &
     &                 argtypedtadiag(1:lenv(argtypedtadiag))
              ENDIF
            CASE (35)
! --- -diffobsref
               IF (largdiffobsref) THEN
                  flagext=4
                  IF (validextvar(argdiffobsref)) flagext=5
                  IF (validextdta(argdiffobsref)) flagext=3
                  CALL infoargext(argdiffobsref,'argdiffobsref',flagext,jmod)
               ENDIF
            CASE (36)
! --- -diffdtaref
               IF (largdiffdtaref) THEN
                  flagext=3
                  IF (validextvar(argdiffdtaref)) flagext=5
                  CALL infoargext(argdiffdtaref,'argdiffdtaref',flagext,jmod)
               ENDIF
            CASE (37)
! --- -diffvarref
               IF (largdiffvarref) THEN
                  flagext=5
                  CALL infoargext(argdiffvarref,'argdiffvarref',flagext,jmod)
               ENDIF
            CASE (38)
! --- -diffobsorg
               IF (largdiffobsorg) THEN
                  flagext=4
                  IF (validextvar(argdiffobsorg)) flagext=5
                  IF (validextdta(argdiffobsorg)) flagext=3
                  CALL infoargext(argdiffobsorg,'argdiffobsorg',flagext,jmod)
               ENDIF
            CASE (39)
! --- -diffdtaorg
               IF (largdiffdtaorg) THEN
                  flagext=3
                  IF (validextvar(argdiffdtaorg)) flagext=5
                  CALL infoargext(argdiffdtaorg,'argdiffdtaorg',flagext,jmod)
               ENDIF
            CASE (40)
! --- -diffvarorg
               IF (largdiffvarorg) THEN
                  flagext=5
                  CALL infoargext(argdiffvarorg,'argdiffvarorg',flagext,jmod)
               ENDIF
            CASE (41)
! --- -inerrdta
               IF (larginerrdta) THEN
                  flagext=3
                  IF (validextvar(arginerrdta)) flagext=5
                  CALL infoargext(arginerrdta,'arginerrdta',flagext,jmod)
               ENDIF
            CASE (42)
! --- -outerrdta
               IF (largouterrdta) THEN
                  flagext=3
                  CALL infoargext(argouterrdta,'argouterrdta',flagext,jmod)
               ENDIF
            CASE (43)
! --- -instddta
               IF (larginstddta) THEN
                  flagext=3
                  IF (validextvar(arginstddta)) flagext=5
                  CALL infoargext(arginstddta,'arginstddta',flagext,jmod)
               ENDIF
            CASE (44)
! --- -outstddta
               IF (largoutstddta) THEN
                  flagext=3
                  CALL infoargext(argoutstddta,'argoutstddta',flagext,jmod)
               ENDIF
            CASE (45)
! --- -inybasref
               IF (larginybasref) THEN
                  flagext=1
                  CALL infoargext(arginybasref,'arginybasref',flagext,jmod)
               ENDIF
            CASE (46)
! --- -outybasref
               IF (largoutybasref) THEN
                  flagext=1
                  CALL infoargext(argoutybasref,'argoutybasref',flagext,jmod)
               ENDIF
            CASE (47)
! --- -outbiasdta
               IF (largoutbiasdta) THEN
                  flagext=3
                  CALL infoargext(argoutbiasdta,'argoutbiasdta',flagext,jmod)
               ENDIF
            CASE (48)
! --- -affectobs
               IF (largaffectobs) THEN
                  WRITE(numout,*) ' observation variable (affectobs) = ', &
     &                 argaffectobs(1:lenv(argaffectobs))
              ENDIF
           CASE (49)
! --- -nullobs
               IF (largnullobs) THEN
                  WRITE(numout,*) ' observation variable (nullobs) = ', &
     &                 argnullobs(1:lenv(argnullobs))
              ENDIF
            CASE (50)
! --- -inpartobs (option 2)
               IF (larginpartobs) THEN
                  WRITE(numout,*) ' obs partition file (inpartobs) = ', &
     &                 arginpartobs(1:lenv(arginpartobs))
               ENDIF
            CASE (51)
! --- -arganamorphosis (option 2)
               IF (larganamorphosis) THEN
                  WRITE(numout,*) ' anamorphosis transformation = ', &
     &                 arganamorphosis(1:lenv(arganamorphosis))
               ENDIF
            CASE (52)
! --- -configobs
               IF (largconfigobs) THEN
                  flagext=4
                  CALL infoargext(argconfigobs,'argconfigobs',flagext,jmod)
               ENDIF
            CASE (53)
! --- -fixjpx (option 1)
               IF (largfixjpx) THEN
                  WRITE(numout,*) ' Vx block size (fixjpx) = ', &
     &                 argfixjpx(1:lenv(argfixjpx))
               ENDIF
            CASE (54)
! --- -reducevar
               IF (largreducevar) THEN
                  flagext=5
                  CALL infoargext(argreducevar,'argreducevar',flagext,jmod)
                  IF (extvarmem(indext(argreducevar,extvartab,nbextvar))/10 &
     &                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (55)
! --- -scale (option 2)
               IF (largscale) THEN
                  WRITE(numout,*) ' observation scaling factor (scale) = ', &
     &                 argscale(1:lenv(argscale))
               ENDIF
            CASE (56)
! --- -biasdbs (option 2)
               IF (largbiasdbs) THEN
                  flagext=2
                  CALL infoargext(argbiasdbs,'argbiasdbs',flagext,jmod)
               ENDIF
            CASE (57)
! --- -inpartvar
               IF (larginpartvar) THEN
                  flagext=5
                  CALL infoargext(arginpartvar,'arginpartvar',flagext,jmod)
                  IF (extvarmem(indext(arginpartvar,extvartab,nbextvar))/10 &
     &                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (58)
! --- -outpartvar
               IF (largoutpartvar) THEN
                  flagext=5
                  CALL infoargext(argoutpartvar,'argoutpartvar',flagext,jmod)
                  IF (extvarmem(indext(argoutpartvar,extvartab,nbextvar))/10 &
     &                   .GT.nallmem) GOTO 102
               ENDIF
            CASE (59)
! --- -inzon
               IF (larginzon) THEN
                  flagext=6
                  CALL infoargext(arginzon,'arginzon',flagext,jmod)
               ENDIF
            CASE (60)
! --- -outzon
               IF (largoutzon) THEN
                  flagext=6
                  CALL infoargext(argoutzon,'argoutzon',flagext,jmod)
               ENDIF
            CASE (61)
! --- -oecorrel
               IF (largoecorrel) THEN
                  WRITE(numout,*) ' observation error correlation matrix = ', &
     &                 argoecorrel(1:lenv(argoecorrel))
               ENDIF
            CASE (62)
! --- -fecorrel
               IF (largfecorrel) THEN
                  WRITE(numout,*) ' forecast error correlation matrix = ', &
     &                 argfecorrel(1:lenv(argfecorrel))
               ENDIF
            CASE (63)
! --- -zonindex
               IF (largzonindex) THEN
                  WRITE(numout,*) ' index in local data section (zonindex) = ', &
     &                 argzonindex(1:lenv(argzonindex))
               ENDIF
            CASE (64)
! --- -incfg
               IF (largincfg) THEN
                  WRITE(numout,*) ' input configuration file (incfg) = ', &
     &                 argincfg(1:lenv(argincfg))
               ENDIF
            CASE (65)
! --- -inybas
               IF (larginybas) THEN
                  flagext=1
                  CALL infoargext(arginybas,'arginybas',flagext,jmod)
               ENDIF
            CASE (66)
! --- -outybas
               IF (largoutybas) THEN
                  flagext=1
                  CALL infoargext(argoutybas,'argoutybas',flagext,jmod)
               ENDIF
            CASE (67)
! --- -inobas
               IF (larginobas) THEN
                  flagext=1
                  CALL infoargext(arginobas,'arginobas',flagext,jmod)
               ENDIF
            CASE (68)
! --- -outobas
               IF (largoutobas) THEN
                  flagext=1
                  CALL infoargext(argoutobas,'argoutobas',flagext,jmod)
               ENDIF
            CASE (69)
! --- -inzbas
               IF (larginzbas) THEN
                  flagext=1
                  CALL infoargext(arginzbas,'arginzbas',flagext,jmod)
               ENDIF
            CASE (70)
! --- -outzbas
               IF (largoutzbas) THEN
                  flagext=1
                  CALL infoargext(argoutzbas,'argoutzbas',flagext,jmod)
               ENDIF
            CASE (71)
! --- -configzon
               IF (largconfigzon) THEN
                  flagext=6
                  CALL infoargext(argconfigzon,'argconfigzon',flagext,jmod)
               ENDIF
            CASE (72)
! --- -coefrmax
               IF (largcoefrmax) THEN
                  WRITE(numout,*) ' maximum ROA coefficient (coefrmax) = ', &
     &                 argcoefrmax(1:lenv(argcoefrmax))
               ENDIF
            CASE (73)
! --- -disable
               IF (largdisable) THEN
                  WRITE(numout,*) ' disable option = ', &
     &                 argdisable(1:lenv(argdisable))
               ENDIF
            CASE (74)
! --- -inobasref
               IF (larginobasref) THEN
                  flagext=1
                  CALL infoargext(arginobasref,'arginobasref',flagext,jmod)
               ENDIF
            CASE (75)
! --- -fixjpz (option 1)
               IF (largfixjpz) THEN
                  WRITE(numout,*) ' maximum number of local ', &
     &                 'data sections in memory (fixjpz) = ', &
     &                 argfixjpz(1:lenv(argfixjpz))
               ENDIF
            CASE (76)
! --- -inparadap (option 2)
               IF (larginparadap) THEN
                  WRITE(numout,*) ' input adaptive parameter file  = ', &
     &                 arginparadap(1:lenv(arginparadap))
              ENDIF
            CASE (77)
! --- -fixjpu (option 1)
               IF (largfixjpu) THEN
                  WRITE(numout,*) ' initial local observation vector ', &
     &                 'size (fixjpu)  = ',argfixjpu(1:lenv(argfixjpu))
               ENDIF
            CASE (78)
! --- -inptzon
               IF (larginptzon) THEN
                  flagext=6
                  CALL infoargext(arginptzon,'arginptzon',flagext,jmod)
               ENDIF
            CASE (79)
! --- -outptzon
               IF (largoutptzon) THEN
                  flagext=6
                  CALL infoargext(argoutptzon,'argoutptzon',flagext,jmod)
               ENDIF
            CASE (80)
! --- -action (option 1)
               IF (largaction) THEN
                  WRITE(numout,*) ' action number = ', &
     &                 argaction(1:lenv(argaction))
               ENDIF
            CASE (81)
! --- -inzbasref
               IF (larginzbasref) THEN
                  flagext=1
                  CALL infoargext(arginzbasref,'arginzbasref',flagext,jmod)
               ENDIF
            CASE (82)
! --- -inzonref
               IF (larginzonref) THEN
                  flagext=6
                  CALL infoargext(arginzonref,'arginzonref',flagext,jmod)
               ENDIF
            CASE (83)
! --- -inoptcfg
               IF (larginoptcfg) THEN
                  WRITE(numout,*) ' optional config file (inoptcfg) = ', &
     &                 arginoptcfg(1:lenv(arginoptcfg))
               ENDIF
            CASE (84)
! --- -connect
               IF (largconnect) THEN
                  WRITE(numout,*) ' connection file = ', &
     &                 argconnect(1:lenv(argconnect))
               ENDIF
            CASE (85)
! --- -incstr
               IF (largincstr) THEN
                  WRITE(numout,*) ' constraint directory = ', &
     &                 argincstr(1:lenv(argincstr))
               ENDIF
            CASE (87)
! --- -insmocfg
               IF (larginsmocfg) THEN
                  WRITE(numout,*) ' optional config file (insmocfg) = ', &
     &                 arginsmocfg(1:lenv(arginsmocfg))
               ENDIF
            CASE (88)
! --- -outsmocfg
               IF (largoutsmocfg) THEN
                  WRITE(numout,*) ' optional config file (outsmocfg) = ', &
     &                 argoutsmocfg(1:lenv(argoutsmocfg))
               ENDIF
            CASE (89)
! --- -inrz
               IF (larginrz) THEN
                  WRITE(numout,*) ' reduced space file  = ', &
     &                 arginrz(1:lenv(arginrz))
               ENDIF
            CASE (90)
! --- -iterate
               IF (largiterate) THEN
                  WRITE(numout,*) ' number of MCMC iterations  = ', &
     &                 argiterate(1:lenv(argiterate))
               ENDIF

            END SELECT  
         ENDDO
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(jmod,1000,1,'readarg','infoarg')
!
 101  WRITE (texterror,*) 'Invalid mask extensiona'
      CALL printerror2(jmod,101,3,'readarg','infoarg',comment=texterror)
 102  WRITE (texterror,*) 'Cannot read Vx vector block by block', &
     &         ' with such file format'
      CALL printerror2(jmod,102,3,'readarg','infoarg',comment=texterror)
 103  WRITE (texterror,*) 'Cannot read Cx covariance block by block', &
     &         ' with such file format'
      CALL printerror2(jmod,103,3,'readarg','infoarg',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE infoargext (argument,text,flagext,jmod)
!---------------------------------------------------------------------
!
!  Purpose : Print information about argument of given extension
!  -------
!  Method : 
!  ------   
!  Input : argument, argument name, type of argument, module
!  -----
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: argument,text
      INTEGER, intent(in) :: flagext,jmod
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jext,xypos
      INTEGER :: jvar,indvar,jdta,inddta
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         SELECT CASE (flagext)
         CASE (1)
! --- bas
            IF ( (.NOT.(validextvarbas(argument))) &
     &        .AND.(.NOT.(validextdtabas(argument))) &
     &        .AND.(.NOT.(validextobsbas(argument))) &
     &        .AND.(.NOT.(validextzonbas(argument))) &
     &         ) GOTO 101
            jext=indext(argument,extbastab,nbextbas)
            WRITE(numout,*) ' ',text(1:lenv(text)),' ='
            WRITE(numout,*) '      argument  = ',argument(1:lenv(argument))
            WRITE(numout,*) '      jextbas   = ',jext
            WRITE(numout,*) '      extension = ',extbastab(jext) &
     &           (1:lenv(extbastab(jext)))
         CASE (2)
! --- dbs
            IF (.NOT.(validextdbs(argument))) GOTO 101
            jext=indext(argument,extdbstab,nbextdbs)
            WRITE(numout,*) ' ',text(1:lenv(text)),' ='
            WRITE(numout,*) '      argument  = ',argument(1:lenv(argument))
            WRITE(numout,*) '      jextdbs   = ',jext
            WRITE(numout,*) '      extension = ',extdbstab(jext) &
     &           (1:lenv(extdbstab(jext)))
         CASE (3)
! --- dta
            IF (.NOT.(validextdta(argument))) GOTO 101
            jext=indext(argument,extdtatab,nbextdta)
            WRITE(numout,*) ' ',text(1:lenv(text)),' ='
            IF (.NOT.(extdtaunit(jext))) THEN
               xypos=posit(argument,etoile)
               DO jdta=1,dtaend
                  inddta=dta_ord(jdta)
                  WRITE (clname,'(6x,a,a,a)') argument(1:(xypos-1)), &
     &                 dtainam(inddta)(1:lenv(dtainam(inddta))), &
     &                 argument((xypos+1):lenv(argument))
                  IF (jdta.EQ.1) THEN
                     WRITE(numout,*) '      argument(',jdta,'<=>' &
     &                 ,dta_nam(inddta)(1:lenv(dta_nam(inddta))) &
     &                 ,')  = ',clname(1:lenv(clname))
                  ELSE
                     WRITE(numout,*) '         ""   (',jdta,'<=>' &
     &                 ,dta_nam(inddta)(1:lenv(dta_nam(inddta))) &
     &                 ,')  = ',clname(1:lenv(clname))
                  ENDIF
               ENDDO
            ELSE
               WRITE(numout,*) '      argument  = ', &
     &              argument(1:lenv(argument))
            ENDIF
            WRITE(numout,*) '      jextdta   = ',jext
            WRITE(numout,*) '      extension = ',extdtatab(jext) &
     &           (1:lenv(extdtatab(jext)))
         CASE (4)
! --- obs
            IF (.NOT.(validextobs(argument))) GOTO 101
            jext=indext(argument,extobstab,nbextobs)
            WRITE(numout,*) ' ',text(1:lenv(text)),' ='
            WRITE(numout,*) '      argument  = ', &
     &              argument(1:lenv(argument))
            WRITE(numout,*) '      jextobs   = ',jext
            WRITE(numout,*) '      extension = ',extobstab(jext) &
     &           (1:lenv(extobstab(jext)))
         CASE (5)
! --- var
            IF (.NOT.(validextvar(argument))) GOTO 101
            jext=indext(argument,extvartab,nbextvar)
            WRITE(numout,*) ' ',text(1:lenv(text)),' ='
            IF (.NOT.(extvarunit(jext))) THEN
               xypos=posit(argument,etoile)
               DO jvar=1,varend
                  indvar=var_ord(jvar)
                  WRITE (clname,'(6x,a,a,a)') argument(1:(xypos-1)), &
     &                 varinam(indvar)(1:lenv(varinam(indvar))), &
     &                 argument((xypos+1):lenv(argument))
                  IF (jvar.EQ.1) THEN
                     WRITE(numout,*) '      argument(',jvar,'<=>' &
     &                 ,var_nam(indvar)(1:lenv(var_nam(indvar))) &
     &                 ,')  = ',clname(1:lenv(clname))
                  ELSE
                     WRITE(numout,*) '         ""   (',jvar,'<=>' &
     &                 ,var_nam(indvar)(1:lenv(var_nam(indvar))) &
     &                 ,')  = ',clname(1:lenv(clname))
                  ENDIF
               ENDDO
            ELSE
               WRITE(numout,*) '      argument  = ', &
     &              argument(1:lenv(argument))
            ENDIF
            WRITE(numout,*) '      jextvar   = ',jext
            WRITE(numout,*) '      extension = ',extvartab(jext) &
     &           (1:lenv(extvartab(jext)))
         CASE (6)
! --- zon
            IF (.NOT.(validextzon(argument))) GOTO 101
            jext=indext(argument,extzontab,nbextzon)
            WRITE(numout,*) ' ',text(1:lenv(text)),' ='
            IF (.NOT.(extzonunit(jext))) THEN
               xypos=posit(argument,etoile)
               DO jdta=1,dtaend
                  inddta=dta_ord(jdta)
                  WRITE (clname,'(6x,a,a,a)') argument(1:(xypos-1)), &
     &                 dtainam(inddta)(1:lenv(dtainam(inddta))), &
     &                 argument((xypos+1):lenv(argument))
                  IF (jdta.EQ.1) THEN
                     WRITE(numout,*) '      argument(',jdta,'<=>' &
     &                 ,dta_nam(inddta)(1:lenv(dta_nam(inddta))) &
     &                 ,')  = ',clname(1:lenv(clname))
                  ELSE
                     WRITE(numout,*) '         ""   (',jdta,'<=>' &
     &                 ,dta_nam(inddta)(1:lenv(dta_nam(inddta))) &
     &                 ,')  = ',clname(1:lenv(clname))
                  ENDIF
               ENDDO
            ELSE
               WRITE(numout,*) '      argument  = ', &
     &              argument(1:lenv(argument))
            ENDIF
            WRITE(numout,*) '      jextzon   = ',jext
            WRITE(numout,*) '      extension = ',extzontab(jext) &
     &           (1:lenv(extzontab(jext)))
!
         CASE DEFAULT
            GOTO 1000
         END SELECT  
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(jmod,1000,1,'readarg','infoargext')
!
 101  WRITE (texterror,*) 'argument ',argument(1:lenv(argument)), &
     &     ' is not valid'
      CALL printerror2(jmod,101,3,'readarg','infoargext', &
     &      comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilarg
