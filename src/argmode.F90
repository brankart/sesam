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
! ---                  ARGMODE.F90                                ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-05 (C.E. Testut)                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE readargmode
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE argmode
      use mod_main
      use mod_cfgxyo
      use utilarg
      use utilvalid
      IMPLICIT NONE
      PRIVATE

      PUBLIC readargmode

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readargmode(jmod,execaction)
!---------------------------------------------------------------------
!
!
!  Purpose : Interpret user command to know which module and
!  -------   which action to execute
!
!  Method : Check all switch arguments
!  ------   Determine which action to perform from list of switches
!           Print information
!
!  Input : user commandline
!  -----
!  Output : nmode, naction (in mod_cfgxyo.F90)
!  ------
!---------------------------------------------------------------------
! modules
! =======
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, intent(in) :: jmod
      TYPE (type_swiarg), dimension(1:nbargmax), intent(in) ::  &
     &     execaction
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=swilg+1) :: swi
      CHARACTER(len=bgword) :: wrd
      INTEGER :: jarg,jcas,jloop,indmod
      INTEGER :: jaction,jorder
      INTEGER :: jextbas,jextdbs,jextdta,jextobs,jextvar,jextzon
      INTEGER :: jext,xypos
      LOGICAL :: sortie,incompatible,laction,lorder
      INTEGER :: jvar,indvar,jdta,inddta,jnbswiopt
      LOGICAL :: warning
      INTEGER :: ierror
      CHARACTER(len=bgword) :: text1,text2
      LOGICAL, dimension(1:nbarg) :: swibool1,swibool2
      LOGICAL, dimension(1:nbmod) :: modbool1
      LOGICAL, dimension(0:nbextbas)  :: extbasbool1
      LOGICAL, dimension(0:nbextdbs)  :: extdbsbool1
      LOGICAL, dimension(0:nbextdta)  :: extdtabool1
      LOGICAL, dimension(0:nbextobs)  :: extobsbool1
      LOGICAL, dimension(0:nbextvar)  :: extvarbool1
      LOGICAL, dimension(0:nbextzon)  :: extzonbool1
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/readarg/readargmode :'
         WRITE(numout,*) '         determine SESAM action from commandline'
      ENDIF
!
! -1.- Initialisation
! -------------------
!
! -1.1- Build table (swibool1) of switches
!       that are possible in module jmod
!
! optional switches
      DO jarg=1,nbarg
         swibool1(jarg)= .FALSE.
         swibool2(jarg)= .FALSE.
         IF (swiopt(jarg).GE.1) THEN
            swibool1(jarg) = swibool(jarg)
         ENDIF
      ENDDO
! switch -mode 
      swibool1(1) = swibool(1)
!
! required switches
      DO jaction=1,nbaction
         IF  (tabvalorder(jmod,jaction).NE.0) THEN
            IF (nborder.LT.(tabvalorder(jmod,jaction)/10)) GOTO 107
            DO jorder=1,((tabvalorder(jmod,jaction)/10))
               jarg=taborder(jmod,jaction,jorder)
               IF (jarg.GT.nbarg) GOTO 106
               swibool1(jarg) = swibool(jarg)
            ENDDO
         ENDIF
      ENDDO
!
! -1.2- Build table (modbool1) in which current
!       module is activated
!
      modbool1(1:nbmod)= .FALSE.
      modbool1(jmod)= modbool(jmod)
      nmode=jmod
!
! -1.3- Copy tables of possible file extensions
!
      extbasbool1(0:nbextbas)= extbasbool(0:nbextbas)
      extdbsbool1(0:nbextdbs)= extdbsbool(0:nbextdbs)
      extdtabool1(0:nbextdta)= extdtabool(0:nbextdta)
      extobsbool1(0:nbextobs)= extobsbool(0:nbextobs)
      extvarbool1(0:nbextvar)= extvarbool(0:nbextvar)
      extzonbool1(0:nbextzon)= extzonbool(0:nbextzon)
!
! -2.- Check all switch arguments
! -------------------------------
      argloop: DO jloop = 1, nbargmax
!
         swi=execaction(jloop)%swi
         wrd=execaction(jloop)%arg
! ==> End of loop
         IF (swi.EQ.'  ') THEN 
            EXIT argloop
         ENDIF
! ==> Write switch and argument
         IF (nprint.GE.1) THEN
            WRITE(numout,*) ' switch/argument: ', swi(1:lenv(swi)), &
     &                                       ' ', wrd(1:lenv(wrd))
         ENDIF
! ==> Check switch validity
         IF (.NOT.(validswi(swi))) GOTO 101
! ==> Check switch validity for this module
         jcas=0
         DO jarg=1,nbarg
            IF (swibool(jarg)) THEN
               IF (swi.EQ.switab(jarg)) THEN
                  IF (swibool1(jarg)) THEN
                     jcas=jarg
                  ELSE
                     GOTO 103
                  ENDIF
               ENDIF                
            ENDIF                
         ENDDO
         IF (jcas.EQ.0) GOTO 1000          
!
! ==> Check argument validity
         CALL affectarg(wrd,jmod,jcas,swibool1(:),modbool1(:), &
     &        extbasbool1(:),extdbsbool1(:),extdtabool1(:), &
     &        extobsbool1(:),extvarbool1(:),extzonbool1(:))
! ==> Set swibool2 to .TRUE. if switch is present in the commandline
         swibool2(jcas)= .TRUE.
!
      END DO argloop
!
! -3.- Determine which action to perform from list of arguments
! -------------------------------------------------------------
! Check if list of arguments corresponds to existing action
!
      incompatible = .FALSE.
      incompatible = ( .NOT.(largmod) .OR. incompatible )
!
      naction=0
      actionloop : DO jaction=1,nbaction
         IF (tabvalorder(jmod,jaction).NE.0) THEN
            laction=.TRUE.
!  Check if all switches of this action are in commandline
            DO jorder=1,((tabvalorder(jmod,jaction)/10))
               jarg=taborder(jmod,jaction,jorder)
               laction=( laction .AND. swibool2(jarg) )
            ENDDO
!  Check if all switches of commandline are in this action
            IF (laction) THEN
               DO jarg=3,nbarg
                  IF (swiopt(jarg).EQ.0) THEN
                     lorder=.FALSE.
                     DO jorder=1,((tabvalorder(jmod,jaction)/10))
                        lorder= ( lorder .OR.  &
     &                       (jarg.EQ.taborder(jmod,jaction,jorder)) )
                     ENDDO
                     IF (.NOT.(lorder)) THEN
                        laction=( laction .AND. .NOT.(swibool2(jarg)) ) 
                     ENDIF   
                  ENDIF              
               ENDDO
            ENDIF
! Set swibool2 to false for all switches of this action
            IF (laction) THEN
               DO jorder=1,((tabvalorder(jmod,jaction)/10))
                  jarg=taborder(jmod,jaction,jorder)
                  swibool2(jarg)=.FALSE.
               ENDDO
! Set naction to jaction and exit loop
               naction=jaction
            ENDIF
         ENDIF
         IF (naction.NE.0) THEN
            EXIT actionloop
         ENDIF
      ENDDO actionloop
!
      incompatible = ( (naction.EQ.0) .OR.incompatible )
      IF (incompatible) GOTO 105
!
! -3.- Check optional switches
! ----------------------------
!
      DO jarg=3,nbarg
         IF ((swibool2(jarg)).AND.(swiopt(jarg).EQ.0)) GOTO 105
         IF ((swibool2(jarg)).AND.(swiopt(jarg).EQ.2)) THEN
            warning=.TRUE.
            DO jnbswiopt=1,nbswiopt
               warning=(warning.AND.(tabnumswiopt(jmod,jnbswiopt).NE.jarg))
            ENDDO
            IF (warning) THEN
               print *,'WARNING: optional switch', &
     &              switab(jarg)(1:lenv(switab(jarg))), &
     &              'is not active in this module'
            ENDIF
         ENDIF
      ENDDO
!
! -4.- Set memory treatment of Vx and Cx vectors
! ----------------------------------------------
! nallmem = 2 : Vx vectors are loaded in memory block by block
! nallmem = 3 : Vx vectors are fully loaded in memory 
! 
      IF (largfixjpx) THEN
         nallmem=2
      ELSE
         nallmem=3
      ENDIF       
!
! -5.- Print information
! ----------------------
! And check if Vx and Cx file formats are compatible
! with memory treatment of Vx and Cx vectors
! 
      CALL infoarg(nmode,naction,swibool1(:),modbool1(:), &
     &        extbasbool1(:),extdbsbool1(:),extdtabool1(:), &
     &        extobsbool1(:),extvarbool1(:),extzonbool1(:))
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'readargmode','readargmode')
 1001 CALL printerror2(0,1001,3,'readargmode','readargmode')
!
 101  WRITE (texterror,*) 'Bad switch: ',swi(1:lenv(swi))
      CALL printerror2(0,101,3,'readargmode','readargmode', &
     &     comment=texterror)
 103  WRITE (texterror,*) 'Switch ', &
     &     swi(1:lenv(swi)),' is invalid for this module (see -help swi)'
      CALL printerror2(0,103,3,'readargmode','readargmode', &
     &     comment=texterror)
 105  WRITE (texterror,*) 'Incorrect list of switches for this module', &
     &     ' (see -help swi)'
      CALL printerror2(0,105,3,'readargmode','readargmode', &
     &     comment=texterror)
 106  WRITE (texterror,*) 'Invalid switch address jarg=',jarg, &
     &     'jmod=',jmod,'jaction=',jaction,'jorder=',jorder
      CALL printerror2(0,106,3,'readargmode','readargmode', &
     &     comment=texterror)
 107  WRITE (texterror,*) 'Insufficient array size: nborder=',nborder
      CALL printerror2(0,107,3,'readargmode','readargmode', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE argmode
