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
C ---                    MODDIFF.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12 (C.E. Testut)                        ---
C --- modification : 99-11 (C.E. Testut)                        ---
C --- modification : 01-06 (C.E. Testut)                        ---
C --- modification : 03-02 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE  moddiff
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE moddiff
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute RMS differences between vector objects
CCC  -------
CCC  Method : Identify action and perform required operations
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use algodiff
      use hiocfg
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: actiondiff,flaganlxyo,flaginpart,flaginreforg
      INTEGER :: jpgroup,jparea
      LOGICAL :: linorgxyo
      CHARACTER(len=bgword) :: inxyo,inrefxyo,inorgxyo
      CHARACTER(len=bgword), dimension(:), allocatable :: nam_grouparea
      INTEGER, dimension(:), allocatable :: tab_nbarea
      INTEGER, dimension(:,:), allocatable :: tab_grouparea
      INTEGER :: allocok
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: DIFF  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Define operations to perform as a function of action index
C ---------------------------------------------------------------
C
      actiondiff=naction
C
C -A- flaganlxyo
C = 1 => operate on var files
C = 2 => operate on dta files
C = 3 => operate on zon files
      flaganlxyo=0
      SELECT CASE (actiondiff)
      CASE (1,4,7,10,13)
         flaganlxyo=1
      CASE (2,5,8,11,14)
         flaganlxyo=2
      CASE (3,6,9,12,15)
         flaganlxyo=3
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -B- flaginreforg
C = 1 => operate with 1 input vector
C = 2 => operate with 2 input vectors
C = 3 => operate with 3 input vectors
      flaginreforg=0
      SELECT CASE (actiondiff)
      CASE (13,14,15)
         flaginreforg=1
      CASE (1,2,3,7,8,9)
         flaginreforg=2
      CASE (4,5,6,10,11,12)
         flaginreforg=3
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -C- flaginpart
C = 0 => compute global (or level by level) RMS differences
C = 1 => compute regional RMS differences
      flaginpart=0
      SELECT CASE (actiondiff)
      CASE (1,2,3,4,5,6)
         flaginpart=0
      CASE (7,8,9,10,11,12,13,14,15)
         flaginpart=1
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C  -2.- Define input and output file names
C ----------------------------------------
C inxyo,inrefxyo,inorgxyo,linorgxyo
C
      inxyo=''
      inrefxyo=''
      inorgxyo=''
      linorgxyo=.FALSE.
      SELECT CASE (flaganlxyo)
      CASE (1)
         inxyo=arginvar
         IF (flaginreforg.GE.2) inrefxyo=argdiffvarref
         IF (flaginreforg.GE.3) THEN
            linorgxyo=.TRUE.
            inorgxyo=argdiffvarorg
         ENDIF
      CASE (2)
         inxyo=argindta
         IF (flaginreforg.GE.2) inrefxyo=argdiffdtaref
         IF (flaginreforg.GE.3) THEN
            linorgxyo=.TRUE.
            inorgxyo=argdiffdtaorg
         ENDIF
      CASE (3)
         inxyo=arginobs
         IF (flaginreforg.GE.2) inrefxyo=argdiffobsref
         IF (flaginreforg.GE.3) THEN
            linorgxyo=.TRUE.
            inorgxyo=argdiffobsorg
         ENDIF
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C  -3.- Read definition of regions (if flaginpart=1)
C --------------------------------------------------
C
      SELECT CASE (flaginpart)
      CASE (0)
C Do nothing
      CASE (1)
C Read header of region configuration file
         CALL evalhdrcfgarea(argincfg,jpgroup,jparea)
C
C --- allocate nam_grouparea
         allocate ( nam_grouparea(1:jpgroup), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         nam_grouparea(:) = 'nam'
C --- allocate tab_nbarea
         allocate ( tab_nbarea(1:jpgroup), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         tab_nbarea(:) = 0
C --- allocate tab_grouparea
         allocate ( tab_grouparea(1:jpgroup,1:jparea), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         tab_grouparea(:,:) = 0
C
C Read region configuration file
         CALL readcfgarea(argincfg,nam_grouparea(:),
     $        tab_nbarea(:),tab_grouparea(:,:))
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -4.- Perform required action
C ----------------------------
C
      SELECT CASE (flaginpart)
      CASE (0)
         CALL algodiffbylev(inxyo,inrefxyo,
     $        inorgxyo,linorgxyo,
     $        flaganlxyo)
      CASE (1)
         IF (flaginreforg.GE.2) THEN
            CALL algodiffbyzon(inxyo,inrefxyo,
     $           inorgxyo,linorgxyo,
     $           flaganlxyo,nam_grouparea(:),
     $           tab_nbarea(:),tab_grouparea(:,:))
         ELSE
            CALL algomeanstdbyzon(inxyo,flaganlxyo,
     $           nam_grouparea(:),
     $           tab_nbarea(:),tab_grouparea(:,:))
         ENDIF
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module DIFF    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'moddiff','moddiff')
 1001 CALL printerror2(0,1001,3,'moddiff','moddiff')
C
      END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
