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
! ---                    DEFCST.F90                               ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12 (C.E. Testut)                        ---
! --- modification : 00-02 (C.E. Testut)                        ---
! --- modification : 03-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  defcst : Define SESAM configuration
! --- SUBROUTINE  defswi : Define SESAM structure
! --- SUBROUTINE  defctl : Load user configuration
! --- SUBROUTINE  defchk : Check user configuration
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE defcst 
!---------------------------------------------------------------------
!
!  Purpose : Define SESAM constants, load and check user configuration
!  -------
!  Method : 1) Define SESAM constants (list of possible actions,
!  ------         list of possible switches, list of possible object
!                 formats, help information, ...)   [defswi]
!           2) Load SESAM user configuration (defined by the user
!                 in 'defcst.control.h').           [defctl]
!           3) Check user configuration             [defchk]
!
!---------------------------------------------------------------------
! module
! ======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!
      IF ((nprint.GE.0).AND.(jproc.EQ.0)) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/defcst :'
         WRITE(numout,*) '         define SESAM configuration'
      ENDIF
!
! -1.- Define SESAM constants (switches, modules, actions,...)
! ------------------------------------------------------------
!
      CALL defswi
!
! -2.- Load SESAM user configuration
! ----------------------------------
!
      CALL defctl
!
! -3.- Check user configuration
! -----------------------------
!
      CALL defchk
!
      RETURN
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE defswi
!---------------------------------------------------------------------
!
!  Purpose : Define SESAM constants (list of possible actions,
!  -------   list of possible switches, list of possible object
!            formats, help information, ...)
!
!  Method : Include file "defcst.switch.h"
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jarg,jmod,jaction,jorder,jhelp
      INTEGER :: jextbas,jextdbs,jextdta,jextobs,jextvar,jextzon
!----------------------------------------------------------------------
!
      IF ((nprint.GE.0).AND.(jproc.EQ.0)) THEN
         WRITE(numout,*) '*** sesam/defcst/defswi :'
         WRITE(numout,*) '         define SESAM structure'
      ENDIF
!
! -1.- Include file with SESAM structure definition
! -------------------------------------------------
!
#include "defcst.switch.h90"
!
      RETURN
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE defctl
!---------------------------------------------------------------------
!
!  Purpose : Load SESAM user configuration (definition of SESAM objects:
!  -------   list of variables, masks, observations structure,...)
!
!  Method : Read SESAM configuration file
!  ------
!---------------------------------------------------------------------
! module
! ======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: indvar,inddta,indobs,inddbs
      INTEGER :: jvar,jdta,jobs,jndbs
!----------------------------------------------------------------------
!
      IF ((nprint.GE.0).AND.(jproc.EQ.0)) THEN
         WRITE(numout,*) '*** ROUTINE sesam/defcst/defctl :'
         WRITE(numout,*) '         load user configuration'
      ENDIF
!
! -1.- Default SESAM configuration
! --------------------------------
!
#include "defcst.control.h90"
!
! -2.- Read user defined SESAM configuration file
! -----------------------------------------------
!
      CALL readlist
!
      RETURN
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE defchk 
!---------------------------------------------------------------------
!
!  Purpose : Check user configuration
!  -------
!
!  Method : Compute number and order of variables in SESAM objects
!  ------   Pre-initialize size and position of variables in SESAM objects
!           Check variable and mask dimensions
!           Write information about user configuration
!
!---------------------------------------------------------------------
! module
! ======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpdbsend, &
     &     jprend,jpoend,jpitpend,jpxend,jpx,jpyend
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: indvar,inddta,indobs,inddbs
      INTEGER :: jvar,jdta,jobs,jndbs
      INTEGER :: jdta1,jobs1,jinf
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** sesam/defcst/defchk :'
         WRITE(numout,*) '         check user configuration'
      ENDIF
!
! -1.- Compute number and order of variables in SESAM objects
! -----------------------------------------------------------
!  Vx: varend, var_ord; Vy: dtaend, dta_ord; Vo: obsend, obs_ord
!  This is computed from user defined 'varend' and 'var_ord'
      jdta1=1
      DO jvar=1,varend
         indvar=var_ord(jvar)
         IF (dta_act(indvar)) THEN
            dta_ord(jdta1)=indvar
            jdta1=jdta1+1
         ENDIF
      ENDDO  
      dtaend=jdta1-1
      DO jvar=jdta1,nbvar
         dta_ord(jvar)= 0
      ENDDO
!
      obsend=0
      jobs1=1
      DO jdta=1,dtaend
         inddta=dta_ord(jdta)
         IF (.NOT.(dta_act(inddta))) GOTO 1000
         DO jobs=jobs1,(jobs1+obsndbs(inddta)-1)
            obs_ord(jobs)=inddta
            obsnord(jobs)=jobs-jobs1+1
         ENDDO
         jobs1=jobs1+obsndbs(inddta)
         obsend=obsend+obsndbs(inddta)
      ENDDO
      DO jobs=obsend+1,nbobs
         obs_ord(jobs)=0
         obsnord(jobs)=0
      ENDDO
!
! -2.- Pre-initialize size and position of variables in SESAM objects
! -------------------------------------------------------------------
!  Vx: var_nbr, var_ind; Vy: dta_nbr, dta_ind; Vo: obs_nbr, obs_ind
!  xxx_nbr: size of the variable; xxx_ind: index of the variable in SESAM Object
!  Check variable and mask dimensions
!
! ==> A Variable configuration
      DO jvar=1,varend
         indvar=var_ord(jvar)
         vardmsk(indvar)=var_dim(indvar)
         IF (var_dim(indvar).NE.vardmsk(indvar)) GOTO 101
         var_nbr(indvar)= 1
         var_ind(indvar)= 1
      ENDDO
! ==> B Data section configuration
      DO jdta=1,dtaend
         inddta=dta_ord(jdta)
         dtadmsk(inddta)=dta_dim(inddta)
         IF (dta_dim(inddta).LT.dtadmsk(inddta)) GOTO 101
         IF (dta_dim(inddta).GT.var_dim(inddta)) GOTO 101
         dta_nbr(inddta)= 1
         dta_ind(inddta)= 1
      ENDDO
! ==> C Observation configuration
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         IF (obs_dim(indobs,inddbs).GT.dta_dim(indobs)) GOTO 101
         obs_nbr(indobs,inddbs)= 1
         obs_ind(indobs,inddbs)= 1
      ENDDO
!
! ==> Check mask dimensions
      DO jvar=1,varend
         indvar=var_ord(jvar)
         IF (var_dim(indvar).GT.vardmsk(indvar)) THEN
            PRINT *,'----WARNING----'
            PRINT *,var_nam(indvar)(1:lenv(var_nam(indvar))), &
     &           ' : mask and variable dimensions are different'
         ENDIF
      ENDDO
      DO jdta=1,dtaend
         inddta=dta_ord(jdta)
         IF (dta_dim(inddta).GT.dtadmsk(inddta)) THEN
            PRINT *,'----WARNING----'
            PRINT *,dta_nam(inddta)(1:lenv(dta_nam(inddta))), &
     &           ' : mask and variable dimensions are different'
         ENDIF
      ENDDO
!
! -2.- Pre-initialize size of SESAM objects
! -----------------------------------------
!  Vx: jpxend; Vy: jpyend; Vo: jpoend
!  Rank of covariance matrix: jprend
!  Number of interpolation points in observation operator: jpitpend
!  Size of database vector: jpdbsend
      jpxend=0
      DO jvar=1,varend
         jpxend=jpxend+var_nbr(var_ord(jvar))
      ENDDO
      jpx=jpxend
!
      jpyend=0
      DO jdta=1,dtaend
         jpyend=jpyend+dta_nbr(dta_ord(jdta))
      ENDDO
!
      jpoend=0
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         jpoend=jpoend+obs_nbr(indobs,inddbs)
      ENDDO
!
      jpitpend=0
      jpdbsend=jpoend
      jprend=0
!
! -3.- Pre-initialize object existence flag
! -----------------------------------------
      existbas=.FALSE.
      existdbs=.FALSE.
      existdta=.FALSE.
      existobs=.FALSE.
      existvar=.FALSE.
!
! -4.- Write information about user configuration
! -----------------------------------------------
      IF(nprint.GE.1) THEN
         WRITE(numout,*) '... Vx object configuration:'
         WRITE(numout,*) ' varend = ',varend
         WRITE(numout,12) '-'
         WRITE(numout,10) 'ord','ind','var','dim','nmax','moy','ect'
         WRITE(numout,12) '-'
         DO jvar=1,varend
            indvar=var_ord(jvar)
            WRITE(numout,11) jvar,var_ord(jvar),var_nam(indvar), &
     &           var_dim(indvar),var_nbr(indvar), &
     &           var_moy(indvar),var_ect(indvar)
         ENDDO
         WRITE(numout,12) '-'
         WRITE(numout,*) '... Vy object configuration:'
         WRITE(numout,*) ' dtaend = ',dtaend
         WRITE(numout,12) '-'
         WRITE(numout,10) 'ord','ind','dta','dim','nmax','moy','ect'
         WRITE(numout,12) '-'
         DO jdta=1,dtaend
            inddta=dta_ord(jdta)
            WRITE(numout,11) jdta,dta_ord(jdta), &
     &           dta_nam(inddta), &
     &           dta_dim(inddta),dta_nbr(inddta), &
     &           dta_moy(inddta),dta_ect(inddta)
         ENDDO
         WRITE(numout,12) '-'
         WRITE(numout,*) '... Vo object configuration:'
         WRITE(numout,*) ' obsend = ',obsend
         WRITE(numout,12) '-'
         WRITE(numout,13) 'obs','dta','dbs','name','dim','size','ref','std'
         WRITE(numout,12) '-'
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            WRITE(numout,14) jobs,obs_ord(jobs),obsnord(jobs), &
     &           obs_nam(indobs,inddbs), &
     &           obs_dim(indobs,inddbs),obs_nbr(indobs,inddbs), &
     &           obs_moy(indobs,inddbs),obs_ect(indobs,inddbs)
         ENDDO
         WRITE(numout,12) '-'
      ENDIF
!
      RETURN
!
! --- format definitions
!
 10   FORMAT ("|",7(2X,A4,2X,"|"))
 11   FORMAT ("|",2(I5,3X,"|"),3X,A5,"|",2X,I3,"D",2X,"|", &
     &     I7,1X,"|",2(E7.1,1X,"|"))
 12   FORMAT (A1,72("-"))
 13   FORMAT ("|",8(2X,A4,2X,"|"))
 14   FORMAT ("|",3(I5,3X,"|"),3X,A5,"|",2X,I3,"D",2X,"|", &
     &     I7,1X,"|",2(E7.1,1X,"|"))
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'defcst','defchk')
!
 101  WRITE (texterror,*) 'incoherent variable dimensions in user configuration'
      CALL printerror2(0,101,3,'defcst','defchk',comment=texterror)
!
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
