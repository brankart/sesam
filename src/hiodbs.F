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
C ---                  HIODBS.F                                 ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 01-06  (C.E. Testut)                       ---
C --- modification : 03-03  (J.M. Brankart)                     ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE evalhdrdbs
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE hiodbs
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrdbs

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrdbs (kfnindbs,jpdbsout)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read header of observation database
CCC  -------
CCC  Method : Call appropriate routine to read the right 'dbs' format
CCC  ------
CCC  Input :  kfnindbs : database filename
CCC  -----
CCC  Output : jpdbsout : database size
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use lioadbs
      use lioncdbs
      use liocrg
      use utilvalid
      IMPLICIT NONE
C----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnindbs
      INTEGER, intent(out) :: jpdbsout
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jextdbs,jobs,indobs,inddbs,ltext
      LOGICAL :: found
CC----------------------------------------------------------------------
C
C Check validity of input file
      IF (.NOT.(validextdbs(kfnindbs))) GOTO 101
      jextdbs=indext(kfnindbs,extdbstab,nbextdbs)
C
      found=.FALSE.
      LOOP1 : DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         ltext=lenv(obs_nam(indobs,inddbs))
         IF ((argaffectobs(1:ltext).EQ.obs_nam(indobs,inddbs)(1:ltext))
     $        .AND.((argaffectobs((ltext+1):(ltext+1))).EQ.(' '))) THEN
            found=.TRUE.
            EXIT LOOP1
         ENDIF
      ENDDO LOOP1
      IF (.NOT.found) GOTO 102
C
C Select the right 'dbs' file format
      SELECT CASE (jextdbs)
      CASE (2)
         CALL evalhdrcrg (kfnindbs,jpdbsout)
      CASE (3)
         CALL evalhdrncdbs (kfnindbs,jobs,jpdbsout)
      CASE (4)
         CALL evalhdradbs (kfnindbs,jpdbsout)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiodbs','evalhdrdbs')
C
 101  WRITE (texterror,*) 'Invalid file extension: ',
     $            kfnindbs(1:lenv(kfnindbs))
      CALL printerror2(0,101,3,'hiodbs','evalhdrdbs',comment=texterror)
 102  WRITE (texterror,*) 'Bad observation name'
      CALL printerror2(0,102,3,'hiodbs','evalhdrdbs',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hiodbs
