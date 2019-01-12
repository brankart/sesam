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
! ---                    MODGLBEOF.F90                            ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 98-06 (C.E. Testut)                        ---
! --- modification : 99-11 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- modification : 03-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE modglbeof
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modglbeof
!---------------------------------------------------------------------
!
!  Purpose : Compute EOF decomposition (global/local/bubble)
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use algobase
      use utilvalid
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actioneof,flaganlxyo,flagcovxyo,flagincovxyoz, &
     &     flagbicovoz,flaggloloc,ios,jchar
      CHARACTER(len=bgword) :: inbasxyo,outbasxyo,incovxyoz,outcovoz
      CHARACTER(len=bgword) :: inpartxyo,inzon,configo,textdisable
      CHARACTER(len=bgword) :: printinbasxyo,printoutbasxyo,printconfigo, &
     &     printincovxyoz,printoutcovoz,printinpartxyo,printinzon, &
     &     vectindep,printlamddabasxyoz
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         SELECT CASE (nmode)
         CASE (10)
           WRITE(numout,*) '&  Running SESAM module: GEOF  &'
         CASE (11)
           WRITE(numout,*) '&  Running SESAM module: LEOF  &'
         CASE (12)
           WRITE(numout,*) '&  Running SESAM module: BEOF  &'
         END SELECT
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of
!      module and action index
!----------------------------------------------------------------
! Set EOF particular action index
      SELECT CASE (nmode)
      CASE (10)
! --- GEOF
         SELECT CASE (naction)
         CASE (1)
            actioneof=1
         CASE (2)
            actioneof=2
         CASE (3)
            actioneof=5
         CASE (4)
            actioneof=6
         CASE (5)
            actioneof=12
         CASE (6)
            actioneof=15
         CASE (7)
            actioneof=16
         CASE (8)
            actioneof=22
         CASE DEFAULT
            GOTO 1000
         END SELECT
      CASE (11)
! --- LEOF
         SELECT CASE (naction)
         CASE (1)
            actioneof=3
         CASE (2)
            actioneof=7
         CASE (3)
            actioneof=13
         CASE (4)
            actioneof=17
         CASE DEFAULT
            GOTO 1000
         END SELECT
      CASE (12)
! --- BEOF
         SELECT CASE (naction)
         CASE (1)
            actioneof=4
         CASE (2)
            actioneof=9
         CASE (3)
            actioneof=10
         CASE (4)
            actioneof=11
         CASE (5)
            actioneof=14
         CASE (6)
            actioneof=19
         CASE (7)
            actioneof=20
         CASE (8)
            actioneof=21
         CASE DEFAULT
            GOTO 1000
         END SELECT
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Define actions to perform
!
! -A- flaganlxyo
! = 1 => EOF output is in Vx space
! = 2 => EOF output is in Vy space
! = 3 => EOF output is in Vo space
      SELECT CASE (actioneof)
      CASE (1,2,3,4,5,6,7,8,9,10,11)
         flaganlxyo=1
      CASE (12,13,14,15,16,17,18,19,20,21)
         flaganlxyo=2
      CASE (22,23,24)
         flaganlxyo=3
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -B- flagcovxyo
! = 1 => Metric is in Vx space
! = 2 => Metric is in Vy space
! = 3 => Metric is in Vo space
      flagcovxyo=0
      SELECT CASE (actioneof)
      CASE (1)
         flagcovxyo=1
      CASE (2,3,4,10,12,13,14,20)
         flagcovxyo=2
      CASE (5,6,7,8,9,11,15,16,17,18,19,21,22,23,24)
         flagcovxyo=3
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -C- flagincovyoz
! = 1 => Input covariance is available in Cx space only
! = 2 => Input covariance is available in Cy space only
! = 3 => Input covariance is also available in Co space
! = 4 => Input covariance is also available in Cz space
      SELECT CASE (actioneof)
      CASE (1)
         flagincovxyoz=1
      CASE (2,3,4,12,13,14)
         flagincovxyoz=2
      CASE (5,6,7,8,9,15,16,17,18,19,22,23)
         flagincovxyoz=3
      CASE (10,11,20,21,24)
         flagincovxyoz=4
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -D- flagbicovoz
! = 0 => Do not output EOF decomposition in Cz space
! = 3 => Output EOF decomposition in Co space as well
! = 4 => Output EOF decomposition in Cz space as well
      SELECT CASE (actioneof)
      CASE (1,2,3,5,7,12,13,15,17,22,23)
         flagbicovoz=0
      CASE (6,8,16,18)
         flagbicovoz=3
      CASE (4,9,10,11,14,19,20,21,24)
         flagbicovoz=4
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -E- flaggloloc
! = 0 => Global EOF decomposition
! = 2 => Local EOF decomposition
! = 3 => Local EOF decomposition
! = 4 => Bubble EOF decomposition
      flaggloloc=0
      SELECT CASE (actioneof)
      CASE (1,2,5,6,12,15,16,22)
         flaggloloc=0
      CASE (3,13)
         flaggloloc=2
      CASE (7,8,17,18,23)
         flaggloloc=3
      CASE (4,9,10,11,14,19,20,21,24)
         flaggloloc=4
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -2.- Define input and output file or directory names
!      Define symbols for printing EOF algorithm
! ---------------------------------------------------------------------
! inbasxyo  : Input covariance matrix (or ensemble of states)
! outbasxyo : Output reduced order covariance matrix (EOF decomposition)
! printinbasxyo  : Symbol for input covariance matrix
! printoutbasxyo : Symbol for output reduced order covariance matrix
      SELECT CASE (flaganlxyo)
      CASE (1)
         inbasxyo=arginxbas
         outbasxyo=argoutxbas
         printinbasxyo='(A)'
         printoutbasxyo='(S)'
      CASE (2)
         inbasxyo=arginybas
         outbasxyo=argoutybas
         printinbasxyo='(HyxA)'
         printoutbasxyo='(HyxS)'
      CASE (3)
         inbasxyo=arginobas
         outbasxyo=argoutobas
         printinbasxyo='(HoxA)'
         printoutbasxyo='(HoxS)'
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! inyo    : Observation vector
! configo : Observation operator
! printconfigo : Symbol for observation operator
      SELECT CASE (flagcovxyo)
      CASE (1)
         configo=' '
         printconfigo=' '
      CASE (2)
         configo=' '
         printconfigo='Hyx'
      CASE (3)
         configo=argconfigobs
         printconfigo='Hox'
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! incovyoz : Covariance matrix from which computing input
!            covariance in observation space
      SELECT CASE (flagincovxyoz)
      CASE (1)
         incovxyoz=arginxbas
      CASE (2)
         incovxyoz=arginybas
      CASE (3)
         incovxyoz=arginobas
      CASE (4)
         incovxyoz=arginzbas
         printconfigo='Hzx'
      CASE DEFAULT
         GOTO 1000
      END SELECT
      WRITE (printincovxyoz,'("(",A,"A)")')  &
     &     printconfigo(1:lenv(printconfigo))
      WRITE (printlamddabasxyoz,'("(",A,"S)")')  &
     &     printconfigo(1:lenv(printconfigo))
!
! inpartxyo : Partition in subsystems
! inzon     : Influence bubble for all subsystems
! printinpartxyo : Symbol for partition in subsystems
! printinzon     : Symbol for influence bubbles
      SELECT CASE (flaggloloc)
      CASE (0)
         inpartxyo=''
         inzon=''
         printinpartxyo=''
         IF (largweight) THEN 
            printinzon=' (W)'
         ELSE
            printinzon=''
         ENDIF
      CASE (2,3,4)
         inpartxyo=arginpartvar
         inzon=arginzon
         printinpartxyo='Pt'
         printinzon=' (Wz)'
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! outcovoz : Reduced order error covariance (EOF decomposition)
!            if computed in Vo or Vz spaces
! printoutcovoz : Symbol for reduced order error covariance
!                 in Vo or Vz spaces
      SELECT CASE (flagbicovoz)
      CASE (0)
         outcovoz=''
         printoutcovoz=''
      CASE (3)
         outcovoz=argoutobas
      CASE (4)
         outcovoz=argoutzbas
      CASE DEFAULT
         GOTO 1000
      END SELECT
      WRITE (printoutcovoz,'("(",A,"S)")')  &
     &     printconfigo(1:lenv(printconfigo))
!
      vectindep='(1/(jpr-1))'
!
! -2.- Perform additional checkings
! ---------------------------------
! Check argweight
      IF (largweight) THEN
      IF (.NOT.(validextvar(argweight) &
     &     .OR.((flagcovxyo.GE.2).AND.(validextdta(argweight))) &
     &     .OR.((flagcovxyo.GE.3).AND.(validextobs(argweight))) )) &
     &     GOTO 112
      ENDIF
! Check argoestd
      IF (largoestd) THEN
      IF (.NOT.(validextvar(argoestd) &
     &     .OR.((flagcovxyo.GE.2).AND.(validextdta(argoestd))) &
     &     .OR.((flagcovxyo.GE.3).AND.(validextobs(argoestd))) ))  &
     &     GOTO 112
      ENDIF
! Check argreducevar
      IF (largreducevar) THEN
      IF (.NOT.( &
     &     ((flaganlxyo.GE.1).AND.(validextvar(argreducevar))) &
     &     .OR.((flaganlxyo.GE.2).AND.(validextdta(argreducevar))) ))  &
     &     GOTO 112
      ENDIF
! Check argdisable
      IF (largdisable) THEN
         DO jchar=1,5
            IF ((argdisable(jchar:jchar).NE.'T') &
     &           .AND.(argdisable(jchar:jchar).NE.'F')) GOTO 113
         ENDDO
      ENDIF
!
      textdisable='TT'
      IF (largdisable) THEN
         textdisable=argdisable(2:2)
      ENDIF
!
! -3.- Print short description of EOF algorithm
! ---------------------------------------------
!
      IF (nprint.GE.1) THEN
      print '(A)',' ------------------------------'
      print '(A,A,A)',' --------- ',modtab(nmode)(1:lenv(modtab(nmode))), &
     &     ' mode ----------'
      print '(A,I2,A)',' ---- action number ',naction, &
     &     ' --------'
      print '(A)',' ------------------------------'
!
      print '(A,I2)',' flaganlxyo = ',flaganlxyo
      print *,'  ',printinbasxyo(1:lenv(printinbasxyo)), &
     &     ' = ',inbasxyo(1:lenv(inbasxyo))
      print *,'  ',printoutbasxyo(1:lenv(printoutbasxyo)), &
     &     ' = ',outbasxyo(1:lenv(outbasxyo))
!
      print '(A,I2)',' flagcovxyo = ',flagcovxyo
      IF (flagcovxyo.GE.3) print *,'  ', &
     &     printconfigo(1:lenv(printconfigo)), &
     &        ' = ',configo(1:lenv(configo))
!
      print '(A,I2)',' flagincovxyoz = ',flagincovxyoz
      print *,'  ',printincovxyoz(1:lenv(printincovxyoz)), &
     &        ' = ',incovxyoz(1:lenv(incovxyoz))
!
      print '(A,I2)',' flagbicovooz = ',flagbicovoz
      IF (flagbicovoz.NE.0) print *,'  ', &
     &     printoutcovoz(1:lenv(printoutcovoz)), &
     &        ' = ',outcovoz(1:lenv(outcovoz))
!
      print '(A,I2)',' flaggloloc = ',flaggloloc
      IF (flaggloloc.NE.0) THEN 
         print *,'  ',printinpartxyo(1:lenv(printinpartxyo)), &
     &        ' = ',inpartxyo(1:lenv(inpartxyo))
         print *,'  ',printinzon(2:lenv(printinzon)), &
     &        ' = ',inzon(1:lenv(inzon))
      ELSE
         IF (largweight) print *,'  ',printinzon(2:lenv(printinzon)), &
     &        ' = ',argweight(1:lenv(argweight))
      ENDIF
!
      print *,' textdisable=',textdisable(1:lenv(textdisable))
!
      print '(A,I2,A)',' ---- Algorithme numero ',naction,' ----'
      print *,' cov = U^ lambda U = ',vectindep(1:lenv(vectindep)),' ', &
     &     printincovxyoz(1:lenv(printincovxyoz)),'^', &
     &     printinzon(1:lenv(printinzon)),' ', &
     &     printincovxyoz(1:lenv(printincovxyoz))
      print *,' Pf = ', &
     &     printoutbasxyo(1:lenv(printoutbasxyo)),' ', &
     &     printoutbasxyo(1:lenv(printoutbasxyo)),'^ = ', &
     &     vectindep(1:lenv(vectindep)),' ', &
     &     printinbasxyo(1:lenv(printinbasxyo)),' U U^ ', &
     &     printinbasxyo(1:lenv(printinbasxyo)),'^'
      IF (flagbicovoz.NE.0) print *,' (', &
     &     printconfigo(1:lenv(printconfigo)),')(Pf)(', &
     &     printconfigo(1:lenv(printconfigo)),')^ = (', &
     &     printconfigo(1:lenv(printconfigo)), &
     &     printoutbasxyo(1:lenv(printoutbasxyo)),') (', &
     &     printconfigo(1:lenv(printconfigo)), &
     &     printoutbasxyo(1:lenv(printoutbasxyo)),')^ = ', &
     &     vectindep(1:lenv(vectindep)),' (', &
     &     printconfigo(1:lenv(printconfigo)), &
     &     printinbasxyo(1:lenv(printinbasxyo)),') U U^ (', &
     &     printconfigo(1:lenv(printconfigo)), &
     &     printinbasxyo(1:lenv(printinbasxyo)),')^'
      print *,' lambda = ', &
     &     printlamddabasxyoz(1:lenv(printlamddabasxyoz)),'^', &
     &     printinzon(1:lenv(printinzon)),' ', &
     &     printlamddabasxyoz(1:lenv(printlamddabasxyoz))
      print '(A)',' ------------------------------'
      ENDIF
!
! -4.- Perform required action (EOF decomposition)
! ------------------------------------------------
!
      CALL calceofs(flaganlxyo,inbasxyo,outbasxyo, &
     &     flagcovxyo, &
     &     flagincovxyoz,incovxyoz, &
     &     flagbicovoz,outcovoz, &
     &     flaggloloc,inpartxyo,inzon,configo, &
     &     textdisable)
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         SELECT CASE (nmode)
         CASE (10)
           WRITE(numout,*) '&  End of SESAM module GEOF    &'
         CASE (11)
           WRITE(numout,*) '&  End of SESAM module LEOF    &'
         CASE (12)
           WRITE(numout,*) '&  End of SESAM module BEOF    &'
         END SELECT
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modglbeof','modglbeof')
!
 112  WRITE (texterror,*) 'Optional argument incompatible with', &
     &     'required action'
      CALL printerror2(0,112,3,'modglbeof','modglbeof', &
     &      comment=texterror)
 113  WRITE (texterror,*) 'Invalid argument: ', &
     &     argdisable(1:lenv(argdisable))
      CALL printerror2(0,113,3,'modglbeof','modglbeof', &
     &      comment=texterror)
!
      END
