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
C ---                    MODCORR.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2004-06 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE modcorr
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modcorr
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute correlation coefficients
CCC  -------
CCC  Method : Identify action and perform required operations
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use algocorr
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: actioncorr, numfila, jxyo, flagxyo, jtype
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
      INTEGER :: ji1,jj1,jk1,indvar1,jvar1,inddta1,jdta1
      INTEGER :: indvarmsk,inddtamsk,indvar,inddta,jx,jy
      INTEGER :: ji,jj,jk,jt,jvar,jdta
      CHARACTER(len=varlg) :: dta_nam1, var_nam1
      LOGICAL :: found
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: CORR  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Define operations to perform as a function of action index
C ---------------------------------------------------------------
C
      actioncorr=naction
C
C -2.- Read index of variable from configuration file
C ---------------------------------------------------
C
      textexclusion='#'
      numfila=10
      CALL openfile(numfila,argincfg)
      text=readnextline(numfila,textexclusion)
      CLOSE(numfila)
C
      SELECT CASE (actioncorr)
      CASE (1)
C
         found=.FALSE.
         READ(text,*,ERR=103) var_nam1,ji1,jj1,jk1,jtype
         DO jvar1 = 1,varend
            indvar1 = var_ord(jvar1)
            IF (var_nam1(1:lenv(var_nam1))
     $          .EQ.var_nam(indvar1)(1:lenv(var_nam(indvar1))) ) THEN
               indvar=indvar1
               jvar=jvar1
               found=.TRUE.
            ENDIF
         ENDDO 
         IF (.NOT.found) GOTO 101
C
         found=.FALSE.
         indvarmsk=jvar-1
         jx=var_ind(indvar)
         DO jt=1,var_jpt(indvar)
         DO jk=1,var_jpk(indvar)
         DO jj=1,var_jpj(indvar)
         DO ji=1,var_jpi(indvar)
            IF (IBITS(mask(ji,jj,jk,jt),indvarmsk,1).NE.0) THEN
              IF ((ji.EQ.ji1).AND.(jj.EQ.jj1).AND.(jk.EQ.jk1)) THEN
                jxyo=jx
                found=.TRUE.
              ENDIF  
              jx=jx+1
            ENDIF  
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         IF (.NOT.found) GOTO 102
      CASE (2)
C
         found=.FALSE.
         READ(text,*,ERR=103) dta_nam1,ji1,jj1,jk1
         DO jdta1 = 1,dtaend
            inddta1 = dta_ord(jdta1)
            IF (dta_nam1(1:lenv(dta_nam1))
     $          .EQ.dta_nam(inddta1)(1:lenv(dta_nam(inddta1))) ) THEN
               inddta=inddta1
               jdta=jdta1
               found=.TRUE.
            ENDIF
         ENDDO 
         IF (.NOT.found) GOTO 101
C
         found=.FALSE.
         inddtamsk=jdta-1+varend
         jy=dta_ind(inddta)
         DO jt=1,var_jpt(indvar)
         DO jk=1,dta_jpk(inddta)
         DO jj=1,dta_jpj(inddta)
         DO ji=1,dta_jpi(inddta)
            IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) THEN
              IF ((ji.EQ.ji1).AND.(jj.EQ.jj1).AND.(jk.EQ.jk1)) THEN
                jxyo=jy
                found=.TRUE.
              ENDIF  
              jy=jy+1
            ENDIF  
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         IF (.NOT.found) GOTO 102
      CASE (3)
         READ(text,*,ERR=103) jxyo
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -3.- Perform required action
C ----------------------------
C
      print *, 'jxyo',jxyo
      SELECT CASE (actioncorr)
      CASE (1)
C Action: -inxbas *.var.bas -outvar *.var -incfg *.cfg
         flagxyo = 1
         CALL calccorr(arginxbas,argoutvar,jxyo,jtype,flagxyo,argconfigobs)
      CASE (2)
C Action: -inybas *.var.bas -outdta *.dta -incfg *.cfg
         flagxyo = 2
         CALL calccorr(arginybas,argoutdta,jxyo,jtype,flagxyo,argconfigobs)
      CASE (3)
C Action: -inobas *.var.bas -outobs *.obs -incfg *.cfg -configobs *.obs
         flagxyo = 3
         CALL calccorr(arginobas,argoutobs,jxyo,jtype,flagxyo,argconfigobs)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module CORR    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'modcorr','modcorr')
C
 101  WRITE (texterror,*) 'invalid variable name'
      CALL printerror2(0,101,3,'modcorr','modcorr',comment=texterror)
 102  WRITE (texterror,*) 'invalid grid location'
      CALL printerror2(0,102,3,'modcorr','modcorr',comment=texterror)
 103  WRITE (texterror,*) 'bad input configuration file'
      CALL printerror2(0,103,3,'modcorr','modcorr',comment=texterror)
C
      END
