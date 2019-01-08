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
C ---                  HIOCNT.F                                 ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 01-06  (J.M. Brankart)                     ---
C --- revised      : 03-03  (J.M. Brankart)                     ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE evalhdrcnt : Read contour file header
C --- SUBROUTINE readcnt    : Read contour file
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE hiocnt
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrcnt,readcnt

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcnt (kfilecnt)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read contour file header
CCC  ------
CCC  Method : Read first uncommented record of contour ASCII file
CCC  ------
CCC  Input :  kfilecnt  : filename
CCC  -----
CCC  Output : Common variables (described in module mod_cont):
CCC  ------   jpc, jppend, jplayend
CCC
CCC---------------------------------------------------------------------
C modules
      use mod_main
      use mod_cfgxyo
      use mod_cont
      use utilfiles
      IMPLICIT NONE
C
C header declarations
      CHARACTER(len=*), intent(in) :: kfilecnt
C
C local declarations
      CHARACTER(len=bgword) :: fline
C
C Open contour file
      CALL openfile(10,kfilecnt)
C
C Read contour file header (1st uncommented line)
C
      fline='#'
      DO WHILE (fline(1:1).EQ.'#')
         READ(10,'(a)') fline
      ENDDO
      READ(fline,*) jpc, jppend, jplayend
C
      CLOSE(10)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiocnt','evalhdrcnt')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcnt (kfilecnt,ktypcnt)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read contour file
CCC  -------
CCC  Method : Loop over contours and read each contour definition
CCC  ------
CCC  Input :  kfilecnt : contour filename
CCC  -----    ktypcnt  : contour type
CCC
CCC  Output : Common variables (in module mod_cont):
CCC  ------   jpp, jplay, contij, contidx, contlevmin, contlevmax
CCC
CCC---------------------------------------------------------------------
C modules
      use mod_main
      use mod_cfgxyo
      use mod_cont
      use utilfiles
      IMPLICIT NONE
C
C header declarations
      CHARACTER(len=*), intent(in) :: kfilecnt
      INTEGER, intent(in) :: ktypcnt
CCC---------------------------------------------------------------------
C local declarations
      INTEGER :: jc, jp, jlay
      CHARACTER(len=bgword) :: fline
CCC---------------------------------------------------------------------
C
C Open contour file
      CALL openfile(10,kfilecnt)
C
C Read contour file header (1st uncommented line)
C
      fline='#'
      DO WHILE (fline(1:1).EQ.'#')
         READ(10,'(a)') fline
      ENDDO
      READ(fline,*) jpc, jppend, jplayend
C
C Loop over contours
      DO jc = 1,jpc
C
C Read number of vertices defining each contour
         fline='#'
         DO WHILE (fline(1:1).EQ.'#')
            READ(10,'(a)') fline
         ENDDO
         READ(fline,*) jpp(jc)
C
C Read position of vertices defining each contour
         DO jp = 1,jpp(jc)
            fline='#'
            DO WHILE (fline(1:1).EQ.'#')
               READ(10,'(a)') fline
            ENDDO
            READ(fline,*) contij(jp,jc)%longi, contij(jp,jc)%latj
         ENDDO
C
C Read number of horizontal slices in cylinders
C defined by the polygonal contours
         fline='#'
         DO WHILE (fline(1:1).EQ.'#')
            READ(10,'(a)') fline
         ENDDO
         READ(fline,*) jplay(jc)
C
C Read minimum and maximum level defining each slice
C Read parameters of linear function (if ktypcnt=1)
C defined inside polygonal slice
         DO jlay = 1,jplay(jc)
            fline='#'
            DO WHILE (fline(1:1).EQ.'#')
               READ(10,'(a)') fline
            ENDDO
C
            SELECT CASE(ktypcnt)
            CASE(0)
               READ(fline,*) contidx(jlay,jc), 
     $                   contlevmin(jlay,jc), contlevmax(jlay,jc)
            CASE(1)
               READ(fline,*) contidx(jlay,jc), 
     $                   contlevmin(jlay,jc), contlevmax(jlay,jc),
     $                   contf0(jlay,jc), contfx(jlay,jc),
     $                   contfy(jlay,jc), contfz(jlay,jc)
            CASE DEFAULT
               GOTO 1000
            END SELECT
C
         ENDDO
      ENDDO
C
      CLOSE(10)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiocnt','readcnt')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hiocnt
