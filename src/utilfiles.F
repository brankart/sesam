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
C ---                    UTILFILES.F                             ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 98-06 (C.E. Testut)                        ---
C --- modification : 01-06 (C.E. Testut)                        ---
C --- modification : 03-03 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE openfile
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilfiles
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
      PRIVATE

      PUBLIC openfile

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE openfile(kunit,kfile,kstatus,kform,kaccess,krecl)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Open file and check if required file is available
CCC  -------   in the right access mode
CCC
CCC  Method : Check possible errors and open file
CCC  ------
CCC  Input : kunit   : FORTRAN unit to open
CCC  -----   kfile   : file to open
CCC          kstatus : status (OLD, NEW or UNKNOWN)
CCC          kform   : format specification (FORMATTED or UNFORMATTED)
CCC          kaccess : access specification (SEQUENTIAL or DIRECT)
CCC          krecl   : record length (for direct access)
CCC---------------------------------------------------------------------
CC modules
CC =======
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfile
      INTEGER, intent(in) :: kunit
      CHARACTER(len=*), optional, intent(in) :: kstatus,kaccess,kform
      INTEGER, optional, intent(in) :: krecl
CC----------------------------------------------------------------------
CC Local declarations
CC ==================
      LOGICAL :: existence
      CHARACTER(len=bgword) :: kstatus1,kaccess1,kform1
CC----------------------------------------------------------------------
C
C Set default values for optional arguments
C -----------------------------------------
C
      kstatus1='OLD'
      kform1='FORMATTED'
      kaccess1='SEQUENTIAL'
      IF (present(kstatus)) kstatus1 = kstatus
      IF (present(kform))   kform1   = kform
      IF (present(kaccess)) kaccess1 = kaccess
C
C Interpret optional arguments and check file
C -------------------------------------------
C
      INQUIRE (UNIT=kunit,OPENED=existence,ERR=101,IOSTAT=iost)
      IF (existence) GOTO 106
C
      INQUIRE (FILE=kfile,EXIST=existence,ERR=101,IOSTAT=iost)
      SELECT CASE(kstatus1(1:1))
      CASE ('o','O')
         kstatus1='OLD'
         IF (.NOT.existence) GOTO 107
      CASE ('n','N')
         kstatus1='NEW'
         IF (existence) GOTO 108
      CASE ('u','U')
         kstatus1='UNKNOWN'
      CASE DEFAULT
         GOTO 102
      END SELECT
C
      SELECT CASE(kform1(1:1))
      CASE ('f','F')
         kform1='FORMATTED'
      CASE ('u','U')
         kform1='UNFORMATTED'
      CASE DEFAULT
         GOTO 103
      END SELECT
C
      SELECT CASE(kaccess1(1:1))
      CASE ('s','S')
         kaccess1='SEQUENTIAL'
      CASE ('d','D')
         kaccess1='DIRECT'
         IF (kform1(1:1).EQ.'F') GOTO 109
         IF (.NOT.present(krecl)) GOTO 105
      CASE DEFAULT
         GOTO 104
      END SELECT
C
C Open file
C ---------
C
      IF ( kaccess1(1:6).EQ.'DIRECT' )  THEN
         OPEN( UNIT=kunit, FILE=kfile, FORM=kform1, ACCESS=kaccess1,
     $        STATUS=kstatus1, RECL=krecl, ERR=101, IOSTAT=iost )
      ELSE
         OPEN( UNIT=kunit, FILE=kfile, FORM=kform1, ACCESS=kaccess1,
     $        STATUS=kstatus1, ERR=101, IOSTAT=iost)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'openfile','openfile')
C
 101  IF (nprint.GE.1) THEN
         WRITE(numout,*) ' ERROR : cannot open file: ', 
     $                                 kfile(1:lenv(kfile))
         WRITE(numout,*) '   + unit   = ',kunit
         IF ( kaccess(1:6).EQ.'DIRECT' )  THEN
           WRITE(numout,*) '   + recl   = ',krecl
         ENDIF
         WRITE(numout,*) '   + status = ',kstatus1(1:lenv(kstatus1))
         WRITE(numout,*) '   + form   = ',kform1(1:lenv(kform1))
         WRITE(numout,*) '   + access = ',kaccess1(1:lenv(kaccess1))
         WRITE(numout,*) '   + iostat = ',iost
      ENDIF
C
      IF ( kaccess(1:6).EQ.'DIRECT' ) THEN
         WRITE (texterror,'(3A,I3,A,I5,7A,I4,A)') 'cannot open file : ',
     $        kfile(1:lenv(kfile)),
     $        ' ; (unit,recl,status,form,access,iostat) = (',
     $        kunit,',',krecl,',',kstatus1(1:lenv(kstatus1)),',',
     $        kform1(1:lenv(kform1)),',',
     $        kaccess1(1:lenv(kaccess1)),',',iost,') '
      ELSE
         WRITE (texterror,'(3A,I3,7A,I4,A)') 'cannot open file : ',
     $        kfile(1:lenv(kfile)),
     $        ' ; (unit,status,form,access,iostat) = (',
     $        kunit,',',kstatus1(1:lenv(kstatus1)),',',
     $        kform1(1:lenv(kform1)),',',
     $        kaccess1(1:lenv(kaccess1)),',',iost,') '
      ENDIF
      CALL printerror2(0,101,3,'openfile','openfile',comment=texterror)
 102  WRITE (texterror,*) 'Invalid open status: ',kstatus1
      CALL printerror2(0,102,1,'openfile','openfile',comment=texterror)
 103  WRITE (texterror,*) 'Invalid open format: ',kform1
      CALL printerror2(0,103,1,'openfile','openfile',comment=texterror)
 104  WRITE (texterror,*) 'Invalid open access: ',kaccess1
      CALL printerror2(0,104,1,'openfile','openfile',comment=texterror)
 105  WRITE (texterror,*) 'Record length required to open file in direct access'
      CALL printerror2(0,105,1,'openfile','openfile',comment=texterror)
 106  WRITE (texterror,*) 'FORTRAN unit already connected: ',kunit
      CALL printerror2(0,106,1,'openfile','openfile',comment=texterror)
 107  WRITE (texterror,*) 'File does not exist: ',kfile
      CALL printerror2(0,107,3,'openfile','openfile',comment=texterror)
 108  WRITE (texterror,*) 'File already exists: ',kfile
      CALL printerror2(0,108,3,'openfile','openfile',comment=texterror)
 109  WRITE (texterror,*) 'Cannot open formatted file in direct access'
      CALL printerror2(0,109,1,'openfile','openfile',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilfiles
