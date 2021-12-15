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
! ---                    UTILSESAM                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12 (C.E. Testut)                        ---
! --- modification : 99-05 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SESAM general utilities
! ---
! --- FUNCTION readnextline : read next uncommented line from text file
! --- FUNCTION lenv : get length of character variable
! --- FUNCTION indext : get index of file extension
! --- FUNCTION posit : get position of specified character in string
! --- FUNCTION mkint : read integer on character string
! --- SUBROUTINE PRINTERROR2 : error display (2nd version)
! --- SUBROUTINE PRINTCOMMENT : print general error comments
! --- SUBROUTINE SHELLORDER : execute shell script
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      CHARACTER(len=2048) FUNCTION readnextline(knumfil, &
     &     ktextexclusion) 
!
! --- Module declaration
      use mod_main , only : lenv,printerror2
      use mod_cfgxyo , only : texterror
      IMPLICIT NONE
! --- Variable declaration
      INTEGER, intent(in) :: knumfil
      CHARACTER(len=*), intent(in) :: ktextexclusion
      INTEGER :: lentxt,badrec
      CHARACTER(len=2048) :: text,filenam
! -----------------------------------------------------------------
      lentxt=lenv(ktextexclusion)
      text=ktextexclusion
!
! Read next line (not beginning with ktextexclusion)
      DO WHILE (text(1:lentxt).EQ.ktextexclusion(1:lentxt))
         READ(UNIT=knumfil,FMT='(a)',ERR=101) text
      ENDDO
!
      readnextline=text
!
      RETURN
!
! --- error management
!
 101  INQUIRE(UNIT=knumfil,NAME=filenam,NEXTREC=badrec)
      WRITE (texterror,*) 'error reading file ',filenam(1:lenv(filenam)), &
     &     ', record ',badrec
      CALL printerror2(0,101,3,'utilsesam','readnextline', &
     &     comment=texterror)
!
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER FUNCTION lenv (ctext) 
!
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: ctext
!
      lenv = len_trim(ctext)
!
      RETURN
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER FUNCTION indext (ctext,ctab,kjpext)
!
! --- Module declaration
      use mod_main , only : lenv,printerror2
      use mod_cfgxyo , only : texterror
      IMPLICIT NONE
! --- Variable declaration 
      INTEGER, intent(in) :: kjpext
      CHARACTER(len=*), intent(in) :: ctext
      CHARACTER(len=*), dimension(:), intent(in) :: ctab
      INTEGER :: kjext, kindext
      INTEGER :: finctext, deb, finctab
      LOGICAL :: lfound
!
      finctext=lenv(ctext)
      kindext=0
      lfound=.FALSE.
!
      IF (finctext.NE.0) THEN
         DO kjext=1,kjpext
            finctab=lenv(ctab(kjext))
            deb=finctext-finctab+1
            IF (deb.LT.1) GOTO 101
            IF (ctext(deb:finctext).EQ.ctab(kjext)(1:finctab)) THEN
               IF (lfound) THEN
                  GOTO 102
               ELSE
                  kindext=kjext
                  lfound=.TRUE.
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
      indext=kindext
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilsesam','indext')
!
 101  WRITE (texterror,*) 'Filename shorter than extension'
      CALL printerror2(0,101,3,'utilsesam','indext', &
     &     comment=texterror)
 102  WRITE (texterror,*) 'Internal error: bad extension definition '
      CALL printerror2(0,102,3,'utilsesam','indext', &
     &     comment=texterror)
!
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER FUNCTION posit (ctext,cwrd)
!
! --- result:
! posit =  n : position of first character of cwrd in ctext
! posit =  0 : no localisation
! posit = -1 : several localisations
!
! --- Module declaration
      use mod_main , only : lenv,printerror2
      use mod_cfgxyo , only : texterror
      IMPLICIT NONE
! --- Variable declaration 
      CHARACTER(len=*), intent(in) :: cwrd, ctext
      INTEGER :: kposit, i0, i1, istart
      LOGICAL :: localise
      INTEGER :: iwrd, itext
!
      localise=.FALSE.
      kposit=0
      iwrd  = lenv (cwrd)
      itext = lenv(ctext)
!
      IF (itext.GE.iwrd) THEN
         istart=itext-iwrd+1
         charloop: DO i0=istart,1,-1
             i1=i0+iwrd-1
             IF (ctext(i0:i1).EQ.cwrd(1:iwrd)) THEN
                IF (localise) THEN
                   kposit=-1
                   EXIT charloop
                ELSE
                   kposit = i0
                ENDIF
             ENDIF
          ENDDO charloop
      ENDIF
!
      posit=kposit
!
      RETURN
!
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER FUNCTION mkint (ctext) 
!
      use mod_main , only : lenv,printerror2
      use mod_cfgxyo , only : texterror
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: ctext
!
      READ(ctext,*,ERR=101) mkint
!
      RETURN
!
! --- error management
!
 101  WRITE (texterror,*) 'Invalid number:',ctext(1:lenv(ctext))
      CALL printerror2(0,101,1,'utilsesam','mkint',comment=texterror)
!
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE PRINTERROR2 (sesamhelp,sesamiost,errortype, &
     &     namefile,nameroutine,comment)
!
! --- Module declaration
      use mod_main , only : lenv,printcomment
      use arghelp
      IMPLICIT NONE
! --- Variable declaration 
      INTEGER, intent(in) :: sesamhelp,sesamiost
      CHARACTER(len=*), intent(in) :: namefile,nameroutine
      INTEGER, intent(in) :: errortype
      CHARACTER(len=*), optional, intent(in) :: comment
!
      INTEGER :: jbar,lcomment
      INTEGER, parameter :: barlength=72
      CHARACTER(len=barlength) :: text,phrase,ctext0,ctext2,ctext3, &
     &     ctext4,ctext5,ctext60,ctext61,ctext62,ctext63,ctext7
      CHARACTER(len=1), parameter :: frame='#'
!
! Build error text to display
!
! --- make ctext0 (#######)
      ctext0=REPEAT(frame,barlength)
! --- make ctext2 (error label: )
      WRITE (text,'(A1,1X,"error label :",1X,I4)') frame(1:1), &
     &        sesamiost
      WRITE (ctext2,'(A,1X,A1)')  &
     &     text(1:barlength-2),frame(1:1)
! --- make ctext3 (source file: )
      WRITE (text,'(A1,1X,"source file :",1X,A)') frame(1:1), &
     &     namefile(1:min(barlength-2-2,lenv(namefile)))
      WRITE (ctext3,'(A,1X,A1)')  &
     &     text(1:barlength-2),frame(1:1)
! --- make ctext4 (subprogram : )
      WRITE (text,'(A1,1X,"subprogram  :",1X,A)') frame(1:1), &
     &     nameroutine(1:min(barlength-2-2,lenv(nameroutine)))
      WRITE (ctext4,'(A,1X,A1)')  &
     &     text(1:barlength-2),frame(1:1)
! --- make ctext5 (error type : )
      SELECT CASE (errortype)
      CASE (1)
         ctext5='program error'
      CASE (2)
         ctext5='configuration error'
      CASE (3)
         ctext5='user error'
      CASE DEFAULT
         ctext5='unknown'
      END SELECT
      WRITE (text,'(A1,1X,"error type  :",1X,A)') frame(1:1), &
     &     ctext5(1:min(barlength-2-2,lenv(ctext5)))
      WRITE (ctext5,'(A,1X,A1)')  &
     &     text(1:barlength-2),frame(1:1)
! --- make ctext6 (comment:    )
      ctext60=''
      ctext61=''
      ctext62=''
      ctext63=''
      IF (.NOT.(present(comment))) THEN
         WRITE (text,'(A1,1X,"description :")') frame(1:1)
         WRITE (ctext60,'(A,1X,A1)') text(1:barlength-2),frame(1:1)
         CALL PRINTCOMMENT(sesamiost,ctext61,ctext62,ctext63)
         WRITE (text,'(A1,3X,A)') frame(1:1), &
     &        ctext61(1:min(barlength-2-4,lenv(ctext61)))
         WRITE (ctext61,'(A,1X,A1)')  text(1:barlength-2),frame(1:1)
         WRITE (text,'(A1,3X,A)') frame(1:1), &
     &        ctext62(1:min(barlength-2-4,lenv(ctext62)))
         WRITE (ctext62,'(A,1X,A1)')  text(1:barlength-2),frame(1:1)
         WRITE (text,'(A1,3X,A)') frame(1:1), &
     &        ctext63(1:min(barlength-2-4,lenv(ctext63)))
         WRITE (ctext63,'(A,1X,A1)') text(1:barlength-2),frame(1:1)
      ELSE
         WRITE (text,'(A1,1X,"description :")') frame(1:1)
         WRITE (ctext60,'(A,1X,A1)') text(1:barlength-2),frame(1:1)
         lcomment=lenv(comment)
         WRITE (text,'(A1,3X,A)') frame(1:1), &
     &        comment(1:min(barlength-2-4,lcomment))
         WRITE (ctext61,'(A,1X,A1)') text(1:barlength-2),frame(1:1)
         IF (lcomment.GT.(barlength-2-4)) THEN
            WRITE (text,'(A1,3X,A)') frame(1:1), &
     &           comment((barlength-2-4+1): &
     &           min(2*(barlength-2-4),lcomment))
            WRITE (ctext62,'(A,1X,A1)') text(1:barlength-2),frame(1:1)
            IF (lcomment.GT.(2*(barlength-2-4))) THEN
               WRITE (text,'(A1,3X,A)') frame(1:1), &
     &              comment(2*(barlength-2-4+1): &
     &              min(3*(barlength-2-4),lcomment))
               WRITE (ctext63,'(A,1X,A1)') text(1:barlength-2), &
     &              frame(1:1)
            ELSE
               text=''
               WRITE (ctext63,'(A1,1X,A,1X,A1)') frame(1:1), &
     &              text(1:barlength-4),frame(1:1)
            ENDIF
         ELSE
            text=''
            WRITE (ctext62,'(A1,1X,A,1X,A1)') frame(1:1), &
     &           text(1:barlength-4),frame(1:1)
            ctext63=ctext62
         ENDIF
      ENDIF
! --- make ctext7 (debugging: )
      SELECT CASE (errortype)
      CASE (1)
         ctext7='contact SESAM wizard or debug yourself'
      CASE (2)
         ctext7='correct SESAM configuration'
      CASE (3)
         ctext7='play again'
      CASE DEFAULT
         ctext7='play again'
      END SELECT
      WRITE (text,'(A1,1X,"debbugging :",1X,A)') frame(1:1), &
     &     ctext7(1:min(barlength-2-2,lenv(ctext7)))
      WRITE (ctext7,'(A,1X,A1)') text(1:barlength-2),frame(1:1)
!
! Display error text
!
      PRINT *,' '
      PRINT *,ctext0(1:lenv(ctext0))
      PRINT *,ctext2(1:lenv(ctext2))
      PRINT *,ctext3(1:lenv(ctext3))
      PRINT *,ctext4(1:lenv(ctext4))
      PRINT *,ctext5(1:lenv(ctext5))
      PRINT *,ctext60(1:lenv(ctext60))
      PRINT *,ctext61(1:lenv(ctext61))
      PRINT *,ctext62(1:lenv(ctext62))
      PRINT *,ctext63(1:lenv(ctext63))
      PRINT *,ctext7(1:lenv(ctext7))
      PRINT *,ctext0(1:lenv(ctext0))
      PRINT *,' '
!
! Display appropriate online help
!
      CALL readarghelp (sesamhelp)
!
      RETURN
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE SHELLORDER (ctext)
!
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: ctext
! 
#if defined _CRAYC90 || defined _CRAYT3E
      CALL ishell(ctext)
#else
      CALL system(ctext)
#endif
      RETURN
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE PRINTCOMMENT(sesamiost,ctext1,ctext2,ctext3)
!
! --- Module declaration
      use mod_main , only : lenv
      IMPLICIT NONE
! --- Variable declaration 
      INTEGER, intent(in) :: sesamiost
      CHARACTER(len=*), intent(out) :: ctext1,ctext2,ctext3
      INTEGER, parameter :: barlength=72
!
      ctext1=''
      ctext2=''
      ctext3=''
!
      SELECT CASE (sesamiost)
      CASE (1000)
         ctext1='SESAM internal inconsistency'
         ctext2='illicit operation'
      CASE (1001)
         ctext1='dynamic allocation error'
         ctext2='not enough memory for this configuration'
         ctext3='use SESAM options: fixjpx, fixjpz'
      CASE (1002)
         ctext1='incompatibility with SESAM settings'
         ctext2='(in defcst.switch.h)'
      CASE (1003)
         ctext1='bad SESAM parallelization settings'
         ctext2='(usually: too many processors for this task)'
      CASE (1004)
         ctext1='incompatible options used'
      CASE DEFAULT
         ctext1='unknown error label'
         ctext2='upgrade SESAM error management'
      END SELECT
!
      ctext1=ctext1(1:min(barlength-2-4,lenv(ctext1)))
      ctext2=ctext2(1:min(barlength-2-4,lenv(ctext2)))
      ctext3=ctext3(1:min(barlength-2-4,lenv(ctext3)))
!     
      RETURN
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
