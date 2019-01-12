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
! ---                    MODGLBROA.F90                            ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 98-06 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- modification : 03-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE modglbroa
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modglbroa
!---------------------------------------------------------------------
!
!  Purpose : Reduced order analysis (global, local or bubble)
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use algoroa
      use utilvalid
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actionroa,flaganlxyo,flagcovyo,flagrefyo, &
     &     flagincovyoz,flagbicovoz,flaggloloc,ios,jchar
      CHARACTER(len=bgword) :: outxyo,inrefxyo,inbasxyo,outbasxyo
      CHARACTER(len=bgword) :: inyo,inrefyo,incovyoz,outcovoz
      CHARACTER(len=bgword) :: inpartxyo,inzon,configo,textdisable
      CHARACTER(len=bgword) :: printoutxyo,printinrefxyo,printinbasxyo, &
     &     printoutbasxyo,printinyo,printconfigo,printinrefyo, &
     &     printincovyoz,printinpartxyo,printinzon,printoutcovoz, &
     &     printbias,printinptzon
      BIGREAL :: coefrmax
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         SELECT CASE (nmode)
         CASE (14)
           WRITE(numout,*) '&  Running SESAM module: GROA  &'
         CASE (15)
           WRITE(numout,*) '&  Running SESAM module: LROA  &'
         CASE (16)
           WRITE(numout,*) '&  Running SESAM module: BROA  &'
         END SELECT
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of
!      module and action index
!----------------------------------------------------------------
! Set ROA particular action index
      SELECT CASE (nmode)
      CASE (14)
! --- GROA
         SELECT CASE (naction)
         CASE (1)
            actionroa=1
         CASE (2)
            actionroa=4
         CASE (3)
            actionroa=5
         CASE (4)
            actionroa=6
         CASE (5)
            actionroa=11
         CASE (6)
            actionroa=14
         CASE (7)
            actionroa=15
         CASE (8)
            actionroa=16
         CASE (9)
            actionroa=21
         CASE DEFAULT
            GOTO 1000
         END SELECT
      CASE (15)
! --- LROA
         SELECT CASE (naction)
         CASE (1)
            actionroa=2
         CASE (2)
            actionroa=7
         CASE (3)
            actionroa=8
         CASE (4)
            actionroa=9
         CASE (5)
            actionroa=12
         CASE (6)
            actionroa=17
         CASE (7)
            actionroa=18
         CASE (8)
            actionroa=19
         CASE (9)
            actionroa=22
         CASE DEFAULT
            GOTO 1000
         END SELECT
      CASE (16)
! --- BROA
         SELECT CASE (naction)
         CASE (1)
            actionroa=3
         CASE (2)
            actionroa=10
         CASE (3)
            actionroa=13
         CASE (4)
            actionroa=20
         CASE (5)
            actionroa=23
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
! = 1 => Analysis output is in Vx space
! = 2 => Analysis output is in Vy space
! = 3 => Analysis output is in Vo space
      SELECT CASE (actionroa)
      CASE (1,2,3,4,5,6,7,8,9,10)
         flaganlxyo=1
      CASE (11,12,13,14,15,16,17,18,19,20)
         flaganlxyo=2
      CASE (21,22,23)
         flaganlxyo=3
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -B- flagcovxy
! = 2 => Observation vector is in Vy space
! = 3 => Observation vector is in Vo space
      SELECT CASE (actionroa)
      CASE (1,2,3,11,12,13)
         flagcovyo=2
      CASE (4,5,6,7,8,9,10,14,15,16,17,18,19,20,21,22,23)
         flagcovyo=3
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -C- flagrefyo
! = 1 => Background (forecast) is available in Vx space only
! = 2 => Background (forecast) is available in Vy space only
! = 3 => Background (forecast) is also available in Vo space
      SELECT CASE (actionroa)
      CASE (1,2,3,4,7)
         flagrefyo=1
      CASE (11,12,13,14,17)
         flagrefyo=2
      CASE (5,6,8,9,10,15,16,18,19,20,21,22,23)
         flagrefyo=3
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -D- flagincovyoz
! = 1 => Background (forecast) error covariance is available in Cx space only
! = 2 => Background (forecast) error covariance is available in Cy space only
! = 3 => Background (forecast) error covariance is also available in Co space
! = 4 => Background (forecast) error covariance is also available in Cz space
      SELECT CASE (actionroa)
      CASE (1,2,4,5,7,8)
         flagincovyoz=1
      CASE (11,12,14,15,17,18)
         flagincovyoz=2
      CASE (6,9,16,19,21,22)
         flagincovyoz=3
      CASE (3,10,13,20,23)
         flagincovyoz=4
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -E- flagbicovoz
! = 0 => Do not output analysis error covariance in Cz space
! = 3 => Output analysis error covariance in Co space as well
! = 4 => Output analysis error covariance in Cz space as well
      flagbicovoz=0
      SELECT CASE (actionroa)
      CASE (1,2,4,5,7,8,11,12,14,15,17,18,21,22)
         flagbicovoz=0
      CASE (6,9,16,19)
         flagbicovoz=3
      CASE (3,10,13,20,23)
         flagbicovoz=4
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -F- flaggloloc
! = 0 => Global reduced order analysis
! = 2 => Local reduced order analysis
! = 3 => Local reduced order analysis
! = 4 => Bubble reduced order analysis
      flaggloloc=0
      SELECT CASE (actionroa)
      CASE (1,4,5,6,11,14,15,16,21)
         flaggloloc=0
      CASE (2,12)
         flaggloloc=2
      CASE (7,8,9,17,18,19,22)
         flaggloloc=3
      CASE (3,10,13,20,23)
         flaggloloc=4
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -2.- Define input and output file or directory names
!      Define symbols for printing ROA algorithm
! -----------------------------------------------------
! outxy    : Analysis output
! inrefxy  : Background (forecast)
! inxybas  : Background (forecast) error covariance
! outxybas : Analysis error covariance
! printoutxyo    : Symbol for analysis output
! printinrefxy   : Symbol for background (forecast)
! printinbasxyo  : Symbol for background (forecast) error covariance
! printoutbasxyo : Symbol for analysis error covariance
      SELECT CASE (flaganlxyo)
      CASE (1)
         outxyo=argoutvar
         inrefxyo=arginvarref
         inbasxyo=arginxbas
         outbasxyo=argoutxbas
         printoutxyo='Xa'
         printinrefxyo='Xf'
         printinbasxyo='(Sf)'
         printoutbasxyo='(Sa)'
      CASE (2)
         outxyo=argoutdta
         inrefxyo=argindtaref
         inbasxyo=arginybas
         outbasxyo=argoutybas
         printoutxyo='Ya'
         printinrefxyo='Yf'
         printinbasxyo='(HyxSf)'
         printoutbasxyo='(HyxSa)'
      CASE (3)
         outxyo=argoutobs
         inrefxyo=arginobsref
         inbasxyo=arginobas
         outbasxyo=argoutobas
         printoutxyo='Oa'
         printinrefxyo='Of'
         printinbasxyo='(HoxSf)'
         printoutbasxyo='(HoxSa)'
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! inyo    : Observation vector
! configo : Observation operator
! printinyo    : Symbol for observation vector
! printconfigo : Symbol for observation operator
      SELECT CASE (flagcovyo)
      CASE (2)
         inyo=argindta
         configo=' '
         printinyo='Yt'
         printconfigo='Hyx'
      CASE (3)
         inyo=arginobs
         configo=argconfigobs
         printinyo='Ot'
         printconfigo='Hox'
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! inrefyo  : Vector from which computing background (forecast)
!            in observation space
! printinrefyo : Symbol for background (forecast) in observation space
      SELECT CASE (flagrefyo)
      CASE (1)
         inrefyo=arginvarref
         WRITE (printinrefyo,'("(",A,"Xf)")')  &
     &        printconfigo(1:lenv(printconfigo))
      CASE (2)
         inrefyo=argindtaref
         WRITE (printinrefyo,'("(",A,"Yf)")')  &
     &        printconfigo(1:lenv(printconfigo))
      CASE (3)
         inrefyo=arginobsref
         WRITE (printinrefyo,'("(",A,"Of)")')  &
     &        printconfigo(1:lenv(printconfigo))

      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! incovyoz : Covariance matrix from which computing background
!            (forecast) error covariance in observation space
      SELECT CASE (flagincovyoz)
      CASE (1)
         incovyoz=arginxbas
      CASE (2)
         incovyoz=arginybas
      CASE (3)
         incovyoz=arginobas
      CASE (4)
         incovyoz=arginzbas
         printconfigo='Hzx'
      CASE DEFAULT
         GOTO 1000
      END SELECT
      WRITE(printincovyoz,'("(",A,"Sf)")')  &
     &     printconfigo(1:lenv(printconfigo))
!
! inpartxyo : Partition in subsystems
! inzon     : Influence bubble for all subsystems
! printinpartxyo : Symbol for partition in subsystems
! printinzon     : Symbol for influence bubbles
! printinptzon   : Symbol for pointers on subsystems
      SELECT CASE (flaggloloc)
      CASE (0)
         inpartxyo=''
         inzon=''
         printinpartxyo=''
         printinzon=''
         printinptzon=''
      CASE (2,3,4)
         inpartxyo=arginpartvar
         inzon=arginzon
         printinpartxyo='Part'
         printinzon=' Wz '
         printinptzon=' Pt '
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! outcovoz : Analysis error covariance
!            if computed in Vo or Vz spaces
! printoutcovoz : Symbol for analysis error covariance
!                 in Vo or Vz spaces
      SELECT CASE (flagbicovoz)
      CASE (0)
         outcovoz=''
         printoutcovoz=''
      CASE (3)
         outcovoz=argoutobas
         WRITE(printoutcovoz,'("(",A,"Sf)")')  &
     &        printconfigo(1:lenv(printconfigo))
      CASE (4)
         outcovoz=argoutzbas
         WRITE(printoutcovoz,'("(",A,"Sf)")') &
     &        printconfigo(1:lenv(printconfigo))
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! printbias : Symbol for constant field added to observations
      printbias=''
      IF (largbias) printbias=' + bias'
!
! -2.- Perform additional checkings
! ---------------------------------
! Check argcoefrmax
      IF (largcoefrmax) THEN
         READ(argcoefrmax,*,IOSTAT=ios) coefrmax
         IF (ios.NE.0) GOTO 114
      ENDIF
! Check argdisable
      IF (largdisable) THEN
         DO jchar=1,5
            IF ((argdisable(jchar:jchar).NE.'T') &
     &           .AND.(argdisable(jchar:jchar).NE.'F')) GOTO 113
         ENDDO
      ENDIF
! Check argweight
      IF (largweight) THEN
        IF (.NOT.(validextvar(argweight) &
     &     .OR.((flagcovyo.GE.2).AND.(validextdta(argweight))) &
     &     .OR.((flagcovyo.GE.3).AND.(validextobs(argweight))) ))  &
     &     GOTO 112
      ENDIF
! Check argoestd
      IF (largoestd) THEN
        IF (.NOT.(validextvar(argoestd) &
     &     .OR.((flagcovyo.GE.2).AND.(validextdta(argoestd))) &
     &     .OR.((flagcovyo.GE.3).AND.(validextobs(argoestd))) ))  &
     &     GOTO 112
      ENDIF
! Check argbias
      IF (largbias) THEN 
         IF (.NOT.(validextvar(argbias) &
     &     .OR.((flagcovyo.GE.2).AND.(validextdta(argbias))) &
     &     .OR.((flagcovyo.GE.3).AND.(validextobs(argbias))) ))  &
     &     GOTO 112
      ENDIF
! Check argreducevar
      IF (largreducevar) THEN
         IF (.NOT.( &
     &     ((flaganlxyo.GE.1).AND.(validextvar(argreducevar))) &
     &     .OR.((flaganlxyo.GE.2).AND.(validextdta(argreducevar))) ))  &
     &     GOTO 1000
      ENDIF
! 
      textdisable='TTTTT'
      IF (largdisable) THEN
         textdisable=argdisable(1:5)
      ENDIF
!
! -3.- Print short description of ROA algorithm
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
      print '(A,I2)', ' flaganlxyo =',flaganlxyo
      print *,'  ',printoutxyo(1:lenv(printoutxyo)), &
     &     ' = ',outxyo(1:lenv(outxyo))
      print *,'  ',printinrefxyo(1:lenv(printinrefxyo)), &
     &     ' = ',inrefxyo(1:lenv(inrefxyo))
      print *,'  ',printinbasxyo(1:lenv(printinbasxyo)), &
     &     ' = ',inbasxyo(1:lenv(inbasxyo))
      print *,'  ',printoutbasxyo(1:lenv(printoutbasxyo)), &
     &     ' = ',outbasxyo(1:lenv(outbasxyo))
!
      print '(A,I2)',' flagcovyo = ',flagcovyo
      print *,'  ',printinyo(1:lenv(printinyo)), &
     &        ' = ',inyo(1:lenv(inyo))
      IF (flagcovyo.GE.3) THEN
         print *,'  ',printconfigo(1:lenv(printconfigo)), &
     &        ' = ',configo(1:lenv(configo))
      ENDIF
!
      print '(A,I2)',' flagrefyo = ',flagrefyo
      print *,'  ',printinrefyo(1:lenv(printinrefyo)), &
     &        ' =',inrefyo(1:lenv(inrefyo))
!
      print '(A,I2)',' flagincovyoz = ',flagincovyoz
      print *,'  ',printincovyoz(1:lenv(printincovyoz)), &
     &        ' = ',incovyoz(1:lenv(incovyoz))
!
      print '(A,I2)',' flagbicvooz = ',flagbicovoz
      IF (flagbicovoz.NE.0)  print *,'  ', &
     &     printoutcovoz(1:lenv(printoutcovoz)), &
     &        ' = ',outcovoz(1:lenv(outcovoz))
!
      print '(A,I2)',' flaggloloc = ',flaggloloc
      IF (flaggloloc.NE.0) THEN 
         print *,'  ', &
     &     printinpartxyo(1:lenv(printinpartxyo)), &
     &        ' = ',inpartxyo(1:lenv(inpartxyo))
         print *,'  ',printinzon(2:lenv(printinzon)), &
     &        ' = ',inzon(1:lenv(inzon))
         IF (larginptzon) print *,'  ', &
     &        printinptzon(1:lenv(printinptzon)), &
     &        ' = ',arginptzon(1:lenv(arginptzon))
      ENDIF
!
      print *,' textdisable=',textdisable(1:lenv(textdisable))
!
      print '(A,I2,A)',' ---- Algorithm number ',naction,' ----'
      print *,' innov = { ',printinyo(1:lenv(printinyo)),' - ', &
     &     printinrefyo(1:lenv(printinrefyo)), &
     &     printbias(1:lenv(printbias)),' }'
      print *,' gamma = [ I + ',printincovyoz(1:lenv(printincovyoz)), &
     &     printinzon(1:lenv(printinzon)),'.Rd', &
     &     printinzon(1:lenv(printinzon)),'.', &
     &     printincovyoz(1:lenv(printincovyoz)),' ]'
      print *,' coefr = gamma',printinzon(1:lenv(printinzon)),'.Rd', &
     &     printinzon(1:lenv(printinzon)),'.innov'
      print *,' ',printoutxyo(1:lenv(printoutxyo)),' = ', &
     &     printinrefxyo(1:lenv(printinrefxyo)),' + ', &
     &     printinbasxyo(1:lenv(printinbasxyo)),'.coefr'
      print *,' ',printoutbasxyo(1:lenv(printoutbasxyo)),' = ', &
     &     printinbasxyo(1:lenv(printinbasxyo)),'.gamma'
      IF (flagbicovoz.NE.0) print *,' ', &
     &     printoutcovoz(1:lenv(printoutcovoz)),' = ', &
     &     printincovyoz(1:lenv(printincovyoz)),'.gamma'
      print '(A)',' ------------------------------'
      ENDIF
!
! -4.- Perform required action (reduced order analysis)
! -----------------------------------------------------
!
      CALL calcroa(flaganlxyo,outxyo,inrefxyo,inbasxyo,outbasxyo, &
     &              flagcovyo,inyo,inrefyo, &
     &              flagincovyoz,incovyoz, &
     &              flagbicovoz,outcovoz, &
     &              flaggloloc,inpartxyo,inzon,configo, &
     &              textdisable)
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         SELECT CASE (nmode)
         CASE (14)
           WRITE(numout,*) '&  End of SESAM module GROA    &'
         CASE (15)
           WRITE(numout,*) '&  End of SESAM module LROA    &'
         CASE (16)
           WRITE(numout,*) '&  End of SESAM module BROA    &'
         END SELECT
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modglbroa','modglbroa')
 1001 CALL printerror2(0,1001,3,'modglbroa','modglbroa')
!
 112  WRITE (texterror,*) 'Optional argument incompatible with', &
     &     'required action'
      CALL printerror2(0,112,3,'modglbroa','modglbroa', &
     &      comment=texterror)
 113  WRITE (texterror,*) 'Invalid argument: ', &
     &     argdisable(1:lenv(argdisable))
      CALL printerror2(0,113,3,'modglbroa','modglbroa', &
     &      comment=texterror)
 114  WRITE (texterror,*) 'Invalid argument: ', &
     &     argcoefrmax(1:lenv(argcoefrmax))
      CALL printerror2(0,114,3,'modglbroa','modglbroa', &
     &      comment=texterror)
!
      END
