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
C ---                    ARGSESAM.F                             ---
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
C --- SUBROUTINE readarg : Interpret user commandline
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE argsesam
      use mod_main
      use argmode
      use arghelp
      use utilvalid
      IMPLICIT NONE
      PRIVATE

      PUBLIC readarg

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readarg (execaction)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Interpret user command to know which module and
CCC  -------   which action to execute
CCC
CCC  Method : Check commandline validity and determine if the user
CCC  ------   wants help or wants to run a SESAM module.
CCC           Further analyse commandline [readarghelp, readargmode].
CCC
CCC  Input : user commandline (routine argument) as
CCC  ------  a list of switch and arguments
CCC
CCC  Output : nmode, naction (in mod_cfgxyo.F)
CCC  ------
CCC---------------------------------------------------------------------
CC module
CC ======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      TYPE (type_swiarg), dimension(1:nbargmax), intent(in) :: 
     $     execaction
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=swilg+1) :: swi
      CHARACTER(len=bgword) :: wrd
      INTEGER :: jmod,jhelp,jloop
      INTEGER :: jargmod,jarghelp,indmod,indhelp
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/readarg :'
         WRITE(numout,*) '         interpret user commandline'
      ENDIF
C
      indmod  = -1
      indhelp = 0
      jargmod = 1
      jarghelp= 2
C
C -1.- Loop on user commandline (switches and arguments)
C ------------------------------------------------------
C Check commandline validity and determine if the user
C wants help or wants to run a SESAM module
C
      argloop: DO jloop = 1, nbargmax
C
         swi=execaction(jloop)%swi
         wrd=execaction(jloop)%arg
C ==> End of loop
         IF (swi.EQ.'  ') THEN 
            IF (jloop.EQ.1) GOTO 107
            EXIT argloop
         ENDIF
C ==> Check switch validity
         IF (swi((swilg+1):(swilg+1)).NE.' ') GOTO 101
         IF ((.NOT.(validswi(swi))).AND.
     $        (swi(1:lenv(swi)).NE.switab(jarghelp)(1:2))) 
     $        GOTO 101
C ==> Check if argument is not empty
         IF ((wrd.EQ.'  ').AND.
     $        (swi(1:lenv(swi)).NE.switab(jarghelp)(1:2)).AND.
     $        (swi.NE.switab(jarghelp))) GOTO 102
C ==> Search for switch '-help', and determine which
C        help information must be displayed
         IF ((swi.EQ.switab(jarghelp)).OR.
     $        (swi(1:lenv(swi)).EQ.switab(jarghelp)(1:2))) THEN
            IF (wrd.EQ.'  ') GOTO 107
            DO jhelp=1,nbhelp
               IF ((helpbool(jhelp)).AND.(wrd.EQ.helptab(jhelp))) THEN
                  IF (indmod.EQ.0) GOTO 106
                  indmod  = 0
                  indhelp = -jhelp
               ENDIF
            ENDDO
            IF (indhelp.EQ.0) THEN
               DO jmod=1,nbmod
                  IF ((modbool(jmod)).AND.(wrd.EQ.modtab(jmod))) THEN
                     indmod  = 0
                     indhelp = jmod
                  ENDIF
               ENDDO
            ENDIF
            IF (indhelp.EQ.0) THEN
               indmod  = 0
               indhelp = (nbmod+1)
            ELSE
               larghelp=.TRUE.
               arghelp = wrd
            ENDIF
         ENDIF
C ==> Search for switch '-mode', and determine which SESAM module to run
         IF (swi.EQ.switab(jargmod)) THEN
            IF (.NOT.(validmod(wrd))) GOTO 103
            DO jmod=1,nbmod
               IF ((modbool(jmod)).AND.(wrd.EQ.modtab(jmod))) THEN
                  IF (indmod.NE.-1) GOTO 105
                  indmod = jmod
               ENDIF
            ENDDO
         ENDIF
C
      END DO argloop
C
C -2.- Analyse SESAM help or SESAM module commandline
C ---------------------------------------------------
C
      SELECT CASE (indmod)
      CASE (0)
C ===> SESAM help
         CALL readarghelp(indhelp)
      CASE (1:nbmod)
C ===> SESAM module
         CALL readargmode(indmod,execaction)
      CASE DEFAULT
         GOTO 104
      END SELECT
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'readarg','readarg')
C
 101  WRITE (texterror,*) 'Bad switch: ',swi(1:lenv(swi))
      CALL printerror2(0,101,3,'readarg','readarg',
     $      comment=texterror)
 102  WRITE (texterror,*) 'Switch without argument: ',
     $                                   swi(1:lenv(swi))
      CALL printerror2(0,102,3,'readarg','readarg',
     $      comment=texterror)
 103  WRITE (texterror,*) 'Invalid SESAM mode: ',
     $                                   wrd(1:lenv(wrd))
      CALL printerror2(0,103,3,'readarg','readarg',
     $      comment=texterror)
 104  WRITE (texterror,*) 'Switch -mode or -help is needed'
      CALL printerror2(0,104,3,'readarg','readarg',
     $      comment=texterror)
 105  WRITE (texterror,*) 'Switch -mode cannot be used more than once'
      CALL printerror2(0,105,3,'readarg','readarg',
     $      comment=texterror)
 106  WRITE (texterror,*) 'Switch -help cannot be used more than once'
      CALL printerror2(0,106,3,'readarg','readarg',
     $      comment=texterror)
 107  PRINT *, '******************************************'
      PRINT *, '*        Welcome in SESAM program        *'
      PRINT *, '*          For help use -help <>         *'
      PRINT *, '******************************************'
      CALL readarghelp ((nbmod+1))
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE argsesam
