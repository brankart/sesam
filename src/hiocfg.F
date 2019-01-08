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
C ---                    HIOCFG.F                               ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 01-06 (C.E. Testut)                        ---
C --- revised      : 03-03 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE evalhdrcfgoper  : Read header of configuration file
C ---                              (mode oper)
C --- SUBROUTINE readcfgoper     : Read configuration file (mode oper)
C ---
C --- SUBROUTINE evalhdrcfgarea  : Read header of configuration file
C ---                              (mode diff)
C --- SUBROUTINE readcfgarea     : Read configuration file (mode diff)
C --- SUBROUTINE evalhdrcfgsmo   : Read header of smoother config files
C --- SUBROUTINE readcfgsmo      : Read smoother configuration files
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE hiocfg
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrcfgoper,readcfgoper,evalhdrcfgarea
      PUBLIC readcfgarea,evalhdrcfgsmo,readcfgsmo

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcfgoper(kargcfgoper,kjpgroup)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read header of configuration file (mode oper)
CCC  -------
CCC  Method : Read first uncommented record of ASCII configuration file
CCC  ------
CCC  Input :  kargcfgoper : filename
CCC  -----
CCC  Output : kjpgroup : Number of filenames in configuration file
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kargcfgoper
      INTEGER, intent(out) :: kjpgroup
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: nbtotarea, delta,ios
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modoper/algooper/evalhdrcfgoper :'
         WRITE(numout,*) '    ==> READING cfg file ',kargcfgoper(1:lenv(kargcfgoper))
      ENDIF
C
C Open configuration file
      CALL openfile(numfil,kargcfgoper)
C
C Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios,ERR=101) kjpgroup
C
      CLOSE (UNIT=numfil)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiocfg','evalhdrcfgoper')
 1001 CALL printerror2(0,1001,3,'hiocfg','evalhdrcfgoper')
C
 101  WRITE (texterror,*) 'Error reading .cfg file'
      CALL printerror2(0,101,3,'hiocfg','evalhdrcfgoper',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcfgoper(kargcfgoper,knam_groupoper)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read configuration file (mode oper)
CCC  -------
CCC  Method : Read uncommented records of ASCII configuration file
CCC  ------
CCC  Input :  kargcfgoper : filename
CCC  -----
CCC  Output : knam_groupoper : List of filenames in configuration files
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kargcfgoper
      CHARACTER(len=*), dimension(:), intent(out) :: knam_groupoper
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jpgroup,jgroup,ios
      CHARACTER(len=hgword) :: text
      CHARACTER(len=bgword) :: grounam,kform
      CHARACTER(len=1) :: textexclusion
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modoper/algooper/readcfgoper :'
         WRITE(numout,*) '    ==> READING cfg file ',kargcfgoper(1:lenv(kargcfgoper))
      ENDIF
C
C Open configuration file
      CALL openfile(numfil,kargcfgoper)
C
C Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios) jpgroup
C
      IF ((jpgroup).NE.size(knam_groupoper,1)) GOTO 1000
C
C Read configuration file
      DO jgroup=1,jpgroup
         knam_groupoper(jgroup)=readnextline(numfil,textexclusion)
      ENDDO
C
C Control print
      IF (nprint.GE.2) THEN
         kform='(8x,a,i4)' 
         WRITE(numout,kform)  '- Number of filenames: ',jpgroup
         kform='(8x,a)' 
         WRITE(numout,kform)  '- List of filenames:'
         kform='(15x,a)' 
         DO jgroup=1,jpgroup
            WRITE(numout,kform) 
     $           knam_groupoper(jgroup)(1:lenv(knam_groupoper(jgroup)))
         ENDDO
      ENDIF

      CLOSE (UNIT=numfil)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiocfg','readcfgoper')
 1001 CALL printerror2(0,1001,3,'hiocfg','readcfgoper')
C
 101  WRITE (texterror,*) 'Error reading .cfg file'
      CALL printerror2(0,101,3,'hiocfg','readcfgoper',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcfgarea(kargcfgarea,kjpgroup,kjparea)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read header of configuration file (mode diff)
CCC  -------
CCC  Method : Read first uncommented record of contour ASCII file
CCC  ------
CCC  Input :  kargcfgarea : filename
CCC  -----
CCC  Output : kjpgroup : Number of regions (group of areas)
CCC  ------   kjparea  : Maximum number of areas (element of a partition)
CCC                      per region
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kargcfgarea
      INTEGER, intent(out) :: kjpgroup,kjparea
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: nbtotarea, delta,ios
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/moddiff/evalhdrcfgarea :'
         WRITE(numout,*) '    ==> READING cfg file ',kargcfgarea(1:lenv(kargcfgarea))
      ENDIF
C
C Open configuration file
      CALL openfile(numfil,kargcfgarea)
C
C Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios) kjpgroup, kjparea, nbtotarea, delta
C
C If error, read simple version of configuration file
      IF ((nbtotarea.LT.1).OR.(ios.NE.0)) THEN
         READ(text,FMT=*,ERR=101) kjpgroup, kjparea
      ELSE
         kjpgroup=kjpgroup*nbtotarea
         kjparea=kjparea*nbtotarea
      ENDIF
C
      CLOSE (UNIT=numfil)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiocfg','evalhdrcfgarea')
C
 101  WRITE (texterror,*) 'Error reading .cfg file'
      CALL printerror2(0,101,3,'hiocfg','evalhdrcfgarea',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcfgarea(kargcfgarea,knam_grouparea,
     $     ktab_nbarea,ktab_grouparea)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read configuration file (mode diff)
CCC  -------
CCC  Method :
CCC  ------
CCC  Input :  kargcfgarea : filename
CCC  -----
CCC  Output : knam_grouparea : List of diagnostic region names
CCC  ------   ktab_nbarea    : Number of areas for each region
CCC           ktab_grouparea : List of areas for each region
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ==================
      CHARACTER(len=*), intent(in) :: kargcfgarea
      CHARACTER(len=*), dimension(:), intent(out) :: knam_grouparea
      INTEGER, dimension(:), intent(out) :: ktab_nbarea
      INTEGER, dimension(:,:), intent(out) :: ktab_grouparea
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jpgroup,jparea,jgroup,nbtotarea,delta,jtotarea,ios
      LOGICAL :: ltotarea
      CHARACTER(len=hgword) :: text
      CHARACTER(len=bgword) :: grounam,kform
      CHARACTER(len=1) :: textexclusion
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/moddiff/readcfgarea :'
         WRITE(numout,*) '    ==> READING cfg file ',kargcfgarea(1:lenv(kargcfgarea))
      ENDIF
C
C Open configuration file
      CALL openfile(numfil,kargcfgarea)
C
C Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios) jpgroup, jparea, nbtotarea, delta
C
C If error, read simple version of configuration file
      IF ((nbtotarea.LT.1).OR.(ios.NE.0)) THEN
         READ(text,FMT=*,ERR=101) jpgroup, jparea
         nbtotarea=1
         delta=0
         ltotarea=.FALSE.
      ELSE
         ltotarea=.TRUE.
      ENDIF
C
C Check size of input arrays
      IF ((jpgroup*nbtotarea).NE.size(knam_grouparea,1)) GOTO 1000
      IF ((jpgroup*nbtotarea).NE.size(ktab_nbarea,1)) GOTO 1000
      IF ((jpgroup*nbtotarea).NE.size(ktab_grouparea,1)) GOTO 1000
      IF ((jparea*nbtotarea).NE.size(ktab_grouparea,2)) GOTO 1000
C
C Loop over diagnostic regions
      DO jgroup=1,jpgroup
C --- Read diagnostic region name
         text=readnextline(numfil,textexclusion)
         READ(text,FMT=*,ERR=101) knam_grouparea(jgroup)
C --- Read number of areas (element of a partition) for this region
         text=readnextline(numfil,textexclusion)
         READ(text,FMT=*,ERR=101) 
     $        ktab_nbarea(jgroup)
         IF (ktab_nbarea(jgroup).GT.jparea) GOTO 103
         IF (ktab_nbarea(jgroup).LE.0) GOTO 103
C --- Read list of area indices for this region
         text=readnextline(numfil,textexclusion)
         READ(text,FMT=*,ERR=101) 
     $        ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))
      ENDDO
C
C Automatically generate additional regions (advanced option)
      DO jtotarea=2,nbtotarea
         DO jgroup=1,jpgroup
            WRITE(text,'(I1,A10)') jtotarea,knam_grouparea(jgroup)
            knam_grouparea(jgroup+(jtotarea-1)*jpgroup) = text
            ktab_nbarea(jgroup+(jtotarea-1)*jpgroup)=
     $           ktab_nbarea(jgroup)*jtotarea
            ktab_grouparea((jgroup+(jtotarea-1)*jpgroup),
     $           1:(ktab_nbarea(jgroup+(jtotarea-2)*jpgroup))) =
     $           ktab_grouparea((jgroup+(jtotarea-2)*jpgroup),
     $           1:(ktab_nbarea(jgroup+(jtotarea-2)*jpgroup)))
            ktab_grouparea(jgroup+(jtotarea-1)*jpgroup,
     $           (1+ktab_nbarea(jgroup+
     $           (jtotarea-2)*jpgroup)):(ktab_nbarea(jgroup+
     $           (jtotarea-1)*jpgroup))) = ktab_grouparea(jgroup,
     $           1:(ktab_nbarea(jgroup)))+delta*(jtotarea-1)
         ENDDO
      ENDDO
C
C Control print
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
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiocfg','readcfgarea')
 1001 CALL printerror2(0,1001,3,'hiocfg','readcfgarea')
C
 101  WRITE (texterror,*) 'Error reading .cfg file'
      CALL printerror2(0,101,3,'hiocfg','readcfgarea',comment=texterror)
 103  WRITE (texterror,*) 'Incoherent parameters in file :',
     $             kargcfgarea(1:lenv(kargcfgarea))
      CALL printerror2(0,103,3,'hiocfg','readcfgarea',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcfgsmo(karginsmocfg,kargoutsmocfg,kjsmo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read header of smoother configuration file
CCC            (mode groa and lroa, options insmocfg and outsmocfg)
CCC  -------
CCC  Method : Read first uncommented record of ASCII file
CCC  ------
CCC  Input :  karginsmocfg,kargoutsmocfg : filenames
CCC  -----
CCC  Output : kjsmo : Number of retrospective analysis steps
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: karginsmocfg,kargoutsmocfg
      INTEGER, intent(out) :: kjsmo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: kjsmoin,kjsmoout,ios
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
     $  '*** ROUTINE : sesam/modglbroa/algoroa/evalhdrcfgsmo :'
         WRITE(numout,*)
     $  '    ==> READING cfg file ',karginsmocfg(1:lenv(karginsmocfg))
         WRITE(numout,*)
     $  '    ==> READING cfg file ',kargoutsmocfg(1:lenv(kargoutsmocfg))
      ENDIF
C
C Open configuration file (in)
      CALL openfile(numfil,karginsmocfg)
C Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios) kjsmoin
      CLOSE (UNIT=numfil)
C
C Open configuration file (out)
      CALL openfile(numfil,kargoutsmocfg)
C Read header of configuration file
      textexclusion='#'
      text=readnextline(numfil,textexclusion)
      READ(text,FMT=*,IOSTAT=ios) kjsmoout
      CLOSE (UNIT=numfil)
C
C Check errors
      IF ((kjsmoin.NE.kjsmoout).OR.(ios.NE.0)) THEN
         print*, 'jsmoin, jsmoout :', kjsmoin,kjsmoout
         GOTO 101
      ENDIF
C
      kjsmo=kjsmoin
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiocfg','evalhdrcfgsmo')
C
 101  WRITE (texterror,*) 'Error reading .cfg file'
      CALL printerror2(0,101,3,'hiocfg','evalhdrcfgsmo',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcfgsmo(karginsmocfg,kargoutsmocfg,kjsmo,
     $           kvctnaminsmo,kvctnamoutsmo,kdirnaminsmo,kdirnamoutsmo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read smoother configuration files
CCC            (mode groa and lroa, options insmocfg and outsmocfg)
CCC  -------
CCC  Method : Read uncommented record of ASCII file
CCC  ------
CCC  Input :  karginsmocfg,kargoutsmocfg : filenames
CCC           kjsmo : Number of retrospective analysis steps
CCC  -----
CCC  Output : kvctnaminsmo : input analysis vector filename
CCC           kvctnamoutsmo: output analysis vector filename
CCC           kdirnaminsmo : input analysis error modes filename
CCC           kdirnamoutsmo: output analysis error modes filename
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: karginsmocfg,kargoutsmocfg
      INTEGER, intent(in) :: kjsmo
      CHARACTER(len=*), dimension(kjsmo), intent(out) ::
     $       kvctnaminsmo,kvctnamoutsmo,kdirnaminsmo,kdirnamoutsmo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: kjsmoin,kjsmoout,ios,js
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
     $  '*** ROUTINE : sesam/modglbroa/algoroa/readcfgsmo :'
         WRITE(numout,*)
     $  '    ==> READING cfg file ',karginsmocfg(1:lenv(karginsmocfg))
         WRITE(numout,*)
     $  '    ==> READING cfg file ',kargoutsmocfg(1:lenv(kargoutsmocfg))
      ENDIF
C
C Open configuration file (in)
      CALL openfile(numfil,karginsmocfg)
C Read header of configuration file
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
C
C Open configuration file (out)
      CALL openfile(numfil,kargoutsmocfg)
C Read header of configuration file
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
C
C Check errors
      IF ((kjsmoin.NE.kjsmoout).OR.(ios.NE.0)) THEN
         GOTO 101
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiocfg','readcfgsmo')
C
 101  WRITE (texterror,*) 'Error reading .cfg file'
      CALL printerror2(0,101,3,'hiocfg','readcfgsmo',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hiocfg
