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
! ---                    HIOCFG.F90                               ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 01-06 (C.E. Testut)                        ---
! --- revised      : 03-03 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE evalhdrcfgoper  : Read header of configuration file
! ---                              (mode oper)
! --- SUBROUTINE readcfgoper     : Read configuration file (mode oper)
! ---
! --- SUBROUTINE evalhdrcfgarea  : Read header of configuration file
! ---                              (mode diff)
! --- SUBROUTINE readcfgarea     : Read configuration file (mode diff)
! --- SUBROUTINE evalhdrcfgsmo   : Read header of smoother config files
! --- SUBROUTINE readcfgsmo      : Read smoother configuration files
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE hiocfg
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrcfgoper,readcfgoper,evalhdrcfgarea
      PUBLIC readcfgarea,evalhdrcfgsmo,readcfgsmo

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcfgoper(kargcfgoper,kjpgroup)
!---------------------------------------------------------------------
!
!  Purpose : Read header of configuration file (mode oper)
!  -------
!  Method : Read first uncommented record of ASCII configuration file
!  ------
!  Input :  kargcfgoper : filename
!  -----
!  Output : kjpgroup : Number of filenames in configuration file
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kargcfgoper
      INTEGER, intent(out) :: kjpgroup
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: nbtotarea, delta,ios
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modoper/algooper/evalhdrcfgoper :'
         WRITE(numout,*) '    ==> READING cfg file ',kargcfgoper(1:lenv(kargcfgoper))
      ENDIF
!
! Open configuration file
      CALL openfile(numfil,kargcfgoper)
!
! Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios,ERR=101) kjpgroup
!
      CLOSE (UNIT=numfil)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiocfg','evalhdrcfgoper')
 1001 CALL printerror2(0,1001,3,'hiocfg','evalhdrcfgoper')
!
 101  WRITE (texterror,*) 'Error reading .cfg file'
      CALL printerror2(0,101,3,'hiocfg','evalhdrcfgoper',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcfgoper(kargcfgoper,knam_groupoper)
!---------------------------------------------------------------------
!
!  Purpose : Read configuration file (mode oper)
!  -------
!  Method : Read uncommented records of ASCII configuration file
!  ------
!  Input :  kargcfgoper : filename
!  -----
!  Output : knam_groupoper : List of filenames in configuration files
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kargcfgoper
      CHARACTER(len=*), dimension(:), intent(out) :: knam_groupoper
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpgroup,jgroup,ios
      CHARACTER(len=hgword) :: text
      CHARACTER(len=bgword) :: grounam,kform
      CHARACTER(len=1) :: textexclusion
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modoper/algooper/readcfgoper :'
         WRITE(numout,*) '    ==> READING cfg file ',kargcfgoper(1:lenv(kargcfgoper))
      ENDIF
!
! Open configuration file
      CALL openfile(numfil,kargcfgoper)
!
! Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios) jpgroup
!
      IF ((jpgroup).NE.size(knam_groupoper,1)) GOTO 1000
!
! Read configuration file
      DO jgroup=1,jpgroup
         knam_groupoper(jgroup)=readnextline(numfil,textexclusion)
      ENDDO
!
! Control print
      IF (nprint.GE.2) THEN
         kform='(8x,a,i4)' 
         WRITE(numout,kform)  '- Number of filenames: ',jpgroup
         kform='(8x,a)' 
         WRITE(numout,kform)  '- List of filenames:'
         kform='(15x,a)' 
         DO jgroup=1,jpgroup
            WRITE(numout,kform)  &
     &           knam_groupoper(jgroup)(1:lenv(knam_groupoper(jgroup)))
         ENDDO
      ENDIF

      CLOSE (UNIT=numfil)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiocfg','readcfgoper')
 1001 CALL printerror2(0,1001,3,'hiocfg','readcfgoper')
!
 101  WRITE (texterror,*) 'Error reading .cfg file'
      CALL printerror2(0,101,3,'hiocfg','readcfgoper',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcfgarea(kargcfgarea,kjpgroup,kjparea)
!---------------------------------------------------------------------
!
!  Purpose : Read header of configuration file (mode diff)
!  -------
!  Method : Read first uncommented record of contour ASCII file
!  ------
!  Input :  kargcfgarea : filename
!  -----
!  Output : kjpgroup : Number of regions (group of areas)
!  ------   kjparea  : Maximum number of areas (element of a partition)
!                      per region
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kargcfgarea
      INTEGER, intent(out) :: kjpgroup,kjparea
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: nbtotarea, delta,ios
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/moddiff/evalhdrcfgarea :'
         WRITE(numout,*) '    ==> READING cfg file ',kargcfgarea(1:lenv(kargcfgarea))
      ENDIF
!
! Open configuration file
      CALL openfile(numfil,kargcfgarea)
!
! Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios) kjpgroup, kjparea, nbtotarea, delta
!
! If error, read simple version of configuration file
      IF ((nbtotarea.LT.1).OR.(ios.NE.0)) THEN
         READ(text,FMT=*,ERR=101) kjpgroup, kjparea
      ELSE
         kjpgroup=kjpgroup*nbtotarea
         kjparea=kjparea*nbtotarea
      ENDIF
!
      CLOSE (UNIT=numfil)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiocfg','evalhdrcfgarea')
!
 101  WRITE (texterror,*) 'Error reading .cfg file'
      CALL printerror2(0,101,3,'hiocfg','evalhdrcfgarea',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcfgarea(kargcfgarea,knam_grouparea, &
     &     ktab_nbarea,ktab_grouparea)
!---------------------------------------------------------------------
!
!  Purpose : Read configuration file (mode diff)
!  -------
!  Method :
!  ------
!  Input :  kargcfgarea : filename
!  -----
!  Output : knam_grouparea : List of diagnostic region names
!  ------   ktab_nbarea    : Number of areas for each region
!           ktab_grouparea : List of areas for each region
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ==================
      CHARACTER(len=*), intent(in) :: kargcfgarea
      CHARACTER(len=*), dimension(:), intent(out) :: knam_grouparea
      INTEGER, dimension(:), intent(out) :: ktab_nbarea
      INTEGER, dimension(:,:), intent(out) :: ktab_grouparea
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpgroup,jparea,jgroup,nbtotarea,delta,jtotarea,ios
      LOGICAL :: ltotarea
      CHARACTER(len=hgword) :: text
      CHARACTER(len=bgword) :: grounam,kform
      CHARACTER(len=1) :: textexclusion
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/moddiff/readcfgarea :'
         WRITE(numout,*) '    ==> READING cfg file ',kargcfgarea(1:lenv(kargcfgarea))
      ENDIF
!
! Open configuration file
      CALL openfile(numfil,kargcfgarea)
!
! Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios) jpgroup, jparea, nbtotarea, delta
!
! If error, read simple version of configuration file
      IF ((nbtotarea.LT.1).OR.(ios.NE.0)) THEN
         READ(text,FMT=*,ERR=101) jpgroup, jparea
         nbtotarea=1
         delta=0
         ltotarea=.FALSE.
      ELSE
         ltotarea=.TRUE.
      ENDIF
!
! Check size of input arrays
      IF ((jpgroup*nbtotarea).NE.size(knam_grouparea,1)) GOTO 1000
      IF ((jpgroup*nbtotarea).NE.size(ktab_nbarea,1)) GOTO 1000
      IF ((jpgroup*nbtotarea).NE.size(ktab_grouparea,1)) GOTO 1000
      IF ((jparea*nbtotarea).NE.size(ktab_grouparea,2)) GOTO 1000
!
! Loop over diagnostic regions
      DO jgroup=1,jpgroup
! --- Read diagnostic region name
         text=readnextline(numfil,textexclusion)
         READ(text,FMT=*,ERR=101) knam_grouparea(jgroup)
! --- Read number of areas (element of a partition) for this region
         text=readnextline(numfil,textexclusion)
         READ(text,FMT=*,ERR=101)  &
     &        ktab_nbarea(jgroup)
         IF (ktab_nbarea(jgroup).GT.jparea) GOTO 103
         IF (ktab_nbarea(jgroup).LE.0) GOTO 103
! --- Read list of area indices for this region
         text=readnextline(numfil,textexclusion)
         READ(text,FMT=*,ERR=101)  &
     &        ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))
      ENDDO
!
! Automatically generate additional regions (advanced option)
      DO jtotarea=2,nbtotarea
         DO jgroup=1,jpgroup
            WRITE(text,'(I1,A10)') jtotarea,knam_grouparea(jgroup)
            knam_grouparea(jgroup+(jtotarea-1)*jpgroup) = text
            ktab_nbarea(jgroup+(jtotarea-1)*jpgroup)= &
     &           ktab_nbarea(jgroup)*jtotarea
            ktab_grouparea((jgroup+(jtotarea-1)*jpgroup), &
     &           1:(ktab_nbarea(jgroup+(jtotarea-2)*jpgroup))) = &
     &           ktab_grouparea((jgroup+(jtotarea-2)*jpgroup), &
     &           1:(ktab_nbarea(jgroup+(jtotarea-2)*jpgroup)))
            ktab_grouparea(jgroup+(jtotarea-1)*jpgroup, &
     &           (1+ktab_nbarea(jgroup+ &
     &           (jtotarea-2)*jpgroup)):(ktab_nbarea(jgroup+ &
     &           (jtotarea-1)*jpgroup))) = ktab_grouparea(jgroup, &
     &           1:(ktab_nbarea(jgroup)))+delta*(jtotarea-1)
         ENDDO
      ENDDO
!
! Control print
      IF (nprint.GE.2) THEN
         kform='(8x,a,i4)'
         WRITE(numout,kform) 'Number of regions: ',jpgroup
         WRITE(numout,kform) 'Maximum number of areas per region: ',jparea
         kform='(8x,a)'
         WRITE(numout,kform)  '- List of regions and list of areas:'
         DO jgroup=1,jpgroup*nbtotarea
            WRITE(numout,*) knam_grouparea(jgroup)(1:lenv(knam_grouparea(jgroup)))
            WRITE(numout,*) ktab_nbarea(jgroup)
            WRITE(numout,*) ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))
         ENDDO
      ENDIF

      CLOSE (UNIT=numfil)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiocfg','readcfgarea')
 1001 CALL printerror2(0,1001,3,'hiocfg','readcfgarea')
!
 101  WRITE (texterror,*) 'Error reading .cfg file'
      CALL printerror2(0,101,3,'hiocfg','readcfgarea',comment=texterror)
 103  WRITE (texterror,*) 'Incoherent parameters in file :', &
     &             kargcfgarea(1:lenv(kargcfgarea))
      CALL printerror2(0,103,3,'hiocfg','readcfgarea',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcfgsmo(karginsmocfg,kargoutsmocfg,kjsmo)
!---------------------------------------------------------------------
!
!  Purpose : Read header of smoother configuration file
!            (mode groa and lroa, options insmocfg and outsmocfg)
!  -------
!  Method : Read first uncommented record of ASCII file
!  ------
!  Input :  karginsmocfg,kargoutsmocfg : filenames
!  -----
!  Output : kjsmo : Number of retrospective analysis steps
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: karginsmocfg,kargoutsmocfg
      INTEGER, intent(out) :: kjsmo
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: kjsmoin,kjsmoout,ios
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) &
     &  '*** ROUTINE : sesam/modglbroa/algoroa/evalhdrcfgsmo :'
         WRITE(numout,*) &
     &  '    ==> READING cfg file ',karginsmocfg(1:lenv(karginsmocfg))
         WRITE(numout,*) &
     &  '    ==> READING cfg file ',kargoutsmocfg(1:lenv(kargoutsmocfg))
      ENDIF
!
! Open configuration file (in)
      CALL openfile(numfil,karginsmocfg)
! Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios) kjsmoin
      CLOSE (UNIT=numfil)
!
! Open configuration file (out)
      CALL openfile(numfil,kargoutsmocfg)
! Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios) kjsmoout
      CLOSE (UNIT=numfil)
!
! Check errors
      IF ((kjsmoin.NE.kjsmoout).OR.(ios.NE.0)) THEN
         print*, 'jsmoin, jsmoout :', kjsmoin,kjsmoout
         GOTO 101
      ENDIF
!
      kjsmo=kjsmoin
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiocfg','evalhdrcfgsmo')
!
 101  WRITE (texterror,*) 'Error reading .cfg file'
      CALL printerror2(0,101,3,'hiocfg','evalhdrcfgsmo',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcfgsmo(karginsmocfg,kargoutsmocfg,kjsmo, &
     &           kvctnaminsmo,kvctnamoutsmo,kdirnaminsmo,kdirnamoutsmo)
!---------------------------------------------------------------------
!
!  Purpose : Read smoother configuration files
!            (mode groa and lroa, options insmocfg and outsmocfg)
!  -------
!  Method : Read uncommented record of ASCII file
!  ------
!  Input :  karginsmocfg,kargoutsmocfg : filenames
!           kjsmo : Number of retrospective analysis steps
!  -----
!  Output : kvctnaminsmo : input analysis vector filename
!           kvctnamoutsmo: output analysis vector filename
!           kdirnaminsmo : input analysis error modes filename
!           kdirnamoutsmo: output analysis error modes filename
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: karginsmocfg,kargoutsmocfg
      INTEGER, intent(in) :: kjsmo
      CHARACTER(len=*), dimension(kjsmo), intent(out) :: &
     &       kvctnaminsmo,kvctnamoutsmo,kdirnaminsmo,kdirnamoutsmo
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: kjsmoin,kjsmoout,ios,js
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) &
     &  '*** ROUTINE : sesam/modglbroa/algoroa/readcfgsmo :'
         WRITE(numout,*) &
     &  '    ==> READING cfg file ',karginsmocfg(1:lenv(karginsmocfg))
         WRITE(numout,*) &
     &  '    ==> READING cfg file ',kargoutsmocfg(1:lenv(kargoutsmocfg))
      ENDIF
!
! Open configuration file (in)
      CALL openfile(numfil,karginsmocfg)
! Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios) kjsmoin
      DO js=1,kjsmo
         text=readnextline(numfil,textexclusion)
         READ(text,FMT=*,IOSTAT=ios) kvctnaminsmo(js)
         text=readnextline(numfil,textexclusion)
         READ(text,FMT=*,IOSTAT=ios) kdirnaminsmo(js)
      ENDDO
      CLOSE (UNIT=numfil)
!
! Open configuration file (out)
      CALL openfile(numfil,kargoutsmocfg)
! Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios) kjsmoout
      DO js=1,kjsmo
         text=readnextline(numfil,textexclusion)
         READ(text,FMT=*,IOSTAT=ios) kvctnamoutsmo(js)
         text=readnextline(numfil,textexclusion)
         READ(text,FMT=*,IOSTAT=ios) kdirnamoutsmo(js)
      ENDDO
      CLOSE (UNIT=numfil)
!
! Check errors
      IF ((kjsmoin.NE.kjsmoout).OR.(ios.NE.0)) THEN
         GOTO 101
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiocfg','readcfgsmo')
!
 101  WRITE (texterror,*) 'Error reading .cfg file'
      CALL printerror2(0,101,3,'hiocfg','readcfgsmo',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hiocfg
