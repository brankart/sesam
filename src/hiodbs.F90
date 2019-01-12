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
! ---                  HIODBS.F90                                 ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 01-06  (C.E. Testut)                       ---
! --- modification : 03-03  (J.M. Brankart)                     ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE evalhdrdbs
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE hiodbs
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrdbs

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrdbs (kfnindbs,jpdbsout)
!---------------------------------------------------------------------
!
!  Purpose : Read header of observation database
!  -------
!  Method : Call appropriate routine to read the right 'dbs' format
!  ------
!  Input :  kfnindbs : database filename
!  -----
!  Output : jpdbsout : database size
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use lioadbs
      use lioncdbs
      use liocrg
      use utilvalid
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnindbs
      INTEGER, intent(out) :: jpdbsout
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jextdbs,jobs,indobs,inddbs,ltext
      LOGICAL :: found
!----------------------------------------------------------------------
!
! Check validity of input file
      IF (.NOT.(validextdbs(kfnindbs))) GOTO 101
      jextdbs=indext(kfnindbs,extdbstab,nbextdbs)
!
      found=.FALSE.
      LOOP1 : DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         ltext=lenv(obs_nam(indobs,inddbs))
         IF ((argaffectobs(1:ltext).EQ.obs_nam(indobs,inddbs)(1:ltext)) &
     &        .AND.((argaffectobs((ltext+1):(ltext+1))).EQ.(' '))) THEN
            found=.TRUE.
            EXIT LOOP1
         ENDIF
      ENDDO LOOP1
      IF (.NOT.found) GOTO 102
!
! Select the right 'dbs' file format
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
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiodbs','evalhdrdbs')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &            kfnindbs(1:lenv(kfnindbs))
      CALL printerror2(0,101,3,'hiodbs','evalhdrdbs',comment=texterror)
 102  WRITE (texterror,*) 'Bad observation name'
      CALL printerror2(0,102,3,'hiodbs','evalhdrdbs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hiodbs
