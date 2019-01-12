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
! ---                    MODCORR.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2004-06 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE modcorr
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modcorr
!---------------------------------------------------------------------
!
!  Purpose : Compute correlation coefficients
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use algocorr
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actioncorr, numfila, jxyo, flagxyo, jtype
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
      INTEGER :: ji1,jj1,jk1,indvar1,jvar1,inddta1,jdta1
      INTEGER :: indvarmsk,inddtamsk,indvar,inddta,jx,jy
      INTEGER :: ji,jj,jk,jt,jvar,jdta
      CHARACTER(len=varlg) :: dta_nam1, var_nam1
      LOGICAL :: found
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: CORR  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
! ---------------------------------------------------------------
!
      actioncorr=naction
!
! -2.- Read index of variable from configuration file
! ---------------------------------------------------
!
      textexclusion='#'
      numfila=10
      CALL openfile(numfila,argincfg)
      text=readnextline(numfila,textexclusion)
      CLOSE(numfila)
!
      SELECT CASE (actioncorr)
      CASE (1)
!
         found=.FALSE.
         READ(text,*,ERR=103) var_nam1,ji1,jj1,jk1,jtype
         DO jvar1 = 1,varend
            indvar1 = var_ord(jvar1)
            IF (var_nam1(1:lenv(var_nam1)) &
     &          .EQ.var_nam(indvar1)(1:lenv(var_nam(indvar1))) ) THEN
               indvar=indvar1
               jvar=jvar1
               found=.TRUE.
            ENDIF
         ENDDO 
         IF (.NOT.found) GOTO 101
!
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
!
         found=.FALSE.
         READ(text,*,ERR=103) dta_nam1,ji1,jj1,jk1
         DO jdta1 = 1,dtaend
            inddta1 = dta_ord(jdta1)
            IF (dta_nam1(1:lenv(dta_nam1)) &
     &          .EQ.dta_nam(inddta1)(1:lenv(dta_nam(inddta1))) ) THEN
               inddta=inddta1
               jdta=jdta1
               found=.TRUE.
            ENDIF
         ENDDO 
         IF (.NOT.found) GOTO 101
!
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
!
! -3.- Perform required action
! ----------------------------
!
      print *, 'jxyo',jxyo
      SELECT CASE (actioncorr)
      CASE (1)
! Action: -inxbas *.var.bas -outvar *.var -incfg *.cfg
         flagxyo = 1
         CALL calccorr(arginxbas,argoutvar,jxyo,jtype,flagxyo,argconfigobs)
      CASE (2)
! Action: -inybas *.var.bas -outdta *.dta -incfg *.cfg
         flagxyo = 2
         CALL calccorr(arginybas,argoutdta,jxyo,jtype,flagxyo,argconfigobs)
      CASE (3)
! Action: -inobas *.var.bas -outobs *.obs -incfg *.cfg -configobs *.obs
         flagxyo = 3
         CALL calccorr(arginobas,argoutobs,jxyo,jtype,flagxyo,argconfigobs)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module CORR    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modcorr','modcorr')
!
 101  WRITE (texterror,*) 'invalid variable name'
      CALL printerror2(0,101,3,'modcorr','modcorr',comment=texterror)
 102  WRITE (texterror,*) 'invalid grid location'
      CALL printerror2(0,102,3,'modcorr','modcorr',comment=texterror)
 103  WRITE (texterror,*) 'bad input configuration file'
      CALL printerror2(0,103,3,'modcorr','modcorr',comment=texterror)
!
      END
