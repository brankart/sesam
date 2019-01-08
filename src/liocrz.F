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
C ---                  LIOCRZ.F                                ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 08-10  (J.M. Brankart)                     ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE evalhdrcrz  : Read '.crz' file header
C --- SUBROUTINE readveccrz  : Read vector in '.crz' file
C --- SUBROUTINE readmat2crz  : Read 2D matrix in '.crz' file
C --- SUBROUTINE readmat3crz : Read 3D matrix in '.crz' file
C --- SUBROUTINE readmat4crz : Read 4D matrix in '.crz' file
C --- SUBROUTINE writehdrcrz : Write '.crz' file header
C --- SUBROUTINE writeveccrz : Write vector in '.crz' file
C --- SUBROUTINE writemat2crz : Write 2D matrix in '.crz' file
C --- SUBROUTINE writemat3crz: Write 3D matrix in '.crz' file
C --- SUBROUTINE writemat4crz: Write 4D matrix in '.crz' file
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liocrz
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrcrz,readveccrz,readmat2crz,readmat3crz,readmat4crz
      PUBLIC writehdrcrz,writeveccrz,writemat2crz,writemat3crz,writemat4crz

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcrz(kfnincrz,jpr,jpz,jpiobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read '.crz' file header
CCC  -------
CCC  Method : Read header of NetCDF file
CCC  ------
CCC  Input :  kfnincrz : filename
CCC  -----
CCC  Output : jpr   : dimension of the reduced space
CCC  ------   jpz   : number of local subsystems
CCC           jpiobs: number of segments in the observation vector
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
C----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrz
      INTEGER, intent(out) :: jpr,jpz,jpiobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=word80) :: kform
      INTEGER :: ierr, idf, idr1, idr2, idz, idiobs
      INTEGER :: jpr1, jpr2
      LOGICAL :: filexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./evalhdrcrz'
         WRITE(numout,*) '    ==> READING file ',
     $                           kfnincrz(1:lenv(kfnincrz))
      ENDIF
C
C -1.- Open NetCDF file
C ---------------------
C
      INQUIRE (FILE=kfnincrz,EXIST=filexists)
      IF (.NOT.filexists) GOTO 101
C
      ierr = NF90_OPEN(kfnincrz,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
C
C -2.- Read and check file dimensions
C -----------------------------------
C
      ierr = NF90_INQ_DIMID(idf,'r1',idr1)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'r2',idr2)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'z',idz)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQ_DIMID(idf,'iobs',idiobs)
      IF (ierr.NE.0) GOTO 103
C
      ierr = NF90_INQUIRE_DIMENSION(idf,idr1,len=jpr1)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,idr2,len=jpr2)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=jpz)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,idiobs,len=jpiobs)
      IF (ierr.NE.0) GOTO 103
C
      IF (jpr1.NE.jpr2) GOTO 104
      jpr = jpr1
C
C -3.- Close NetCDF file
C ----------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
C -4.- Control print
C ------------------
C
      IF (nprint.GE.3) THEN
         kform='(8x,a,i5)'
         WRITE(numout,kform) '- Dimension of the reduced space: ',jpr
         WRITE(numout,kform) '- Number of local subsystems: ',jpz
         WRITE(numout,kform) '- Number of observation segments: ',jpiobs
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrz','evalhdrcrz')
 1001 CALL printerror2(0,1001,3,'liocrz','evalhdrcrz')
C
 101  WRITE (texterror,*) '.crz file does not exist',kfnincrz
      CALL printerror2(0,101,3,'liocrz','evalhdrcrz',comment=texterror)
 102  WRITE (texterror,*) 'Bad NetCDF file: ',kfnincrz
      CALL printerror2(0,102,3,'liocrz','evalhdrcrz',comment=texterror)
 103  WRITE (texterror,*) 'Error reading dimensions in file: ',kfnincrz
      CALL printerror2(0,103,3,'liocrz','evalhdrcrz',comment=texterror)
 104  WRITE (texterror,*) 'Inconsistent .crz file: ',kfnincrz
      CALL printerror2(0,104,3,'liocrz','evalhdrcrz',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfnincrz
      CALL printerror2(0,107,3,'liocrz','evalhdrcrz',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readveccrz(kfnincrz,kvectnam,kvect)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read vector in '.crz' file
CCC  -------
CCC  Method : Read vector in NetCDF file
CCC  ------
CCC  Input :  kfnincrz    : filename
CCC  -----    kvectnam    : variable name
CCC
CCC  Output : kvect       : vector
CCC  ------ 
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrz, kvectnam
      BIGREAL, dimension(:), intent(out) :: kvect
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpssize
      CHARACTER(len=word80) :: kform
      BIGREAL4, dimension(:), allocatable :: ptabo
C
      INTEGER :: ierr, idf, idv, idd, ndims
      INTEGER, dimension(1) :: vstart,vcount,idims
      LOGICAL :: filexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readveccrz'
         WRITE(numout,*) '    ==> READING file ',
     $                        kfnincrz(1:lenv(kfnincrz))
      ENDIF
C
C Set output vector size
      jpssize = size(kvect,1)
C
C -1.- Open NetCDF file
C ---------------------
C
      INQUIRE (FILE=kfnincrz,EXIST=filexists)
      IF (.NOT.filexists) GOTO 101
C
      ierr = NF90_OPEN(kfnincrz,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
C
C -2.- Read and check vector dimension
C ------------------------------------
C Get variable id
      ierr = NF90_INQ_VARID(idf,kvectnam,idv)
      IF (ierr.NE.0) GOTO 104
C Get variable number dimensions
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,ndims=ndims)
      IF (ierr.NE.0) GOTO 104
C Check if it is a 1D vector
      IF (ndims.NE.1) GOTO 104
C Get dimension id
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,dimids=idims)
      IF (ierr.NE.0) GOTO 104
C Get vector size
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(1),len=vcount(1))
      IF (ierr.NE.0) GOTO 103
C Check if it is equal to output vector size
      IF (vcount(1).NE.jpssize) GOTO 105
C
C -3.- Read vector from file
C --------------------------
C
      vstart=1
      allocate ( ptabo(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
      ierr = NF90_GET_VAR(idf,idv,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
C
      kvect(:) = FREAL(ptabo(:))
C
C --- deallocation ptabo
      IF (allocated(ptabo)) deallocate (ptabo)
C
C -4.- Close NetCDF file
C ----------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrz','readveccrz')
 1001 CALL printerror2(0,1001,3,'liocrz','readveccrz')
C
 101  WRITE (texterror,*) 'File does not exist',kfnincrz
      CALL printerror2(0,101,3,'liocrz','readveccrz',comment=texterror)
 102  WRITE (texterror,*) 'Bad NetCDF file: ',kfnincrz
      CALL printerror2(0,102,3,'liocrz','readveccrz',comment=texterror)
 103  WRITE (texterror,*) 'Error reading dimension in file: ',kfnincrz
      CALL printerror2(0,103,3,'liocrz','readveccrz',comment=texterror)
 104  WRITE (texterror,*) 'Error reading variable in file: ',kfnincrz
      CALL printerror2(0,104,3,'liocrz','readveccrz',comment=texterror)
 105  WRITE (texterror,*) 'Inconsistent dimensions in file: ',kfnincrz
      CALL printerror2(0,105,3,'liocrz','readveccrz',comment=texterror)
 106  WRITE (texterror,*) 'Error reading vector in file: ',kfnincrz
      CALL printerror2(0,106,3,'liocrz','readveccrz',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfnincrz
      CALL printerror2(0,107,3,'liocrz','readveccrz',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readmat2crz(kfnincrz,kmatnam,kmat)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read 2D matrix in '.crz' file
CCC  -------
CCC  Method : Read matrix in NetCDF file
CCC  ------
CCC  Input :  kfnincrz    : filename
CCC  -----    kmatnam     : variable name
CCC
CCC  Output : kmat        : matrix
CCC  ------   
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrz, kmatnam
      BIGREAL, dimension(:,:), intent(out) :: kmat
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpssize1,jpssize2
      CHARACTER(len=word80) :: kform
      BIGREAL4, dimension(:,:), allocatable :: ptabo
C
      INTEGER :: ierr, idf, idv, ndims
      INTEGER, dimension(2) :: vstart,vcount,idims
      LOGICAL :: filexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readmat2crz'
         WRITE(numout,*) '    ==> READING file ',
     $                        kfnincrz(1:lenv(kfnincrz))
      ENDIF
C
C Set output matrix size
      jpssize1 = size(kmat,1)
      jpssize2 = size(kmat,2)
C
C -1.- Open NetCDF file
C ---------------------
C
      INQUIRE (FILE=kfnincrz,EXIST=filexists)
      IF (.NOT.filexists) GOTO 101
C
      ierr = NF90_OPEN(kfnincrz,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
C
C -2.- Read and check matrix dimension
C ------------------------------------
C Get variable id
      ierr = NF90_INQ_VARID(idf,kmatnam,idv)
      IF (ierr.NE.0) GOTO 104
C Get variable number dimensions
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,ndims=ndims)
      IF (ierr.NE.0) GOTO 104
C Check if it is a 2D matrix
      IF (ndims.NE.2) GOTO 104
C Get dimension id
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,dimids=idims)
      IF (ierr.NE.0) GOTO 104
C Get matrix size
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(1),len=vcount(1))
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(2),len=vcount(2))
      IF (ierr.NE.0) GOTO 103
C Check if it is equal to output matrix size
      IF (vcount(1).NE.jpssize1) GOTO 105
      IF (vcount(2).NE.jpssize2) GOTO 105
C
C -3.- Read matrix file
C ---------------------
C
      vstart=(/1,1/)
      allocate ( ptabo(1:jpssize1,1:jpssize2), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
      ierr = NF90_GET_VAR(idf,idv,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
C
      kmat(:,:) = FREAL(ptabo(:,:))
C
C --- deallocation
      IF (allocated(ptabo)) deallocate(ptabo)
C
C -4.- Close observation file
C ---------------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrz','readmat2crz')
 1001 CALL printerror2(0,1001,3,'liocrz','readmat2crz')
C
 101  WRITE (texterror,*) 'File does not exist',kfnincrz
      CALL printerror2(0,101,3,'liocrz','readmat2crz',comment=texterror)
 102  WRITE (texterror,*) 'Bad NetCDF file: ',kfnincrz
      CALL printerror2(0,102,3,'liocrz','readmat2crz',comment=texterror)
 103  WRITE (texterror,*) 'Error reading dimension in file: ',kfnincrz
      CALL printerror2(0,103,3,'liocrz','readmat2crz',comment=texterror)
 104  WRITE (texterror,*) 'Error reading variable in file: ',kfnincrz
      CALL printerror2(0,104,3,'liocrz','readmat2crz',comment=texterror)
 105  WRITE (texterror,*) 'Inconsistent dimensions in file: ',kfnincrz
      CALL printerror2(0,105,3,'liocrz','readmat2crz',comment=texterror)
 106  WRITE (texterror,*) 'Error reading matrix in file: ',kfnincrz
      CALL printerror2(0,106,3,'liocrz','readmat2crz',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfnincrz
      CALL printerror2(0,107,3,'liocrz','readmat2crz',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readmat3crz(kfnincrz,kmatnam,kmat)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read 3D matrix in '.crz' file
CCC  -------
CCC  Method : Read matrix in NetCDF file
CCC  ------
CCC  Input :  kfnincrz    : filename
CCC  -----    kmatnam    : variable name
CCC
CCC  Output : kmat        : matrix
CCC  ------   
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrz, kmatnam
      BIGREAL, dimension(:,:,:), intent(out) :: kmat
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpssize1,jpssize2,jpssize3
      CHARACTER(len=word80) :: kform
      BIGREAL4, dimension(:,:,:), allocatable :: ptabo
C
      INTEGER :: ierr, idf, idv, ndims
      INTEGER, dimension(3) :: vstart,vcount,idims
      LOGICAL :: filexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readmat3crz'
         WRITE(numout,*) '    ==> READING file ',
     $                        kfnincrz(1:lenv(kfnincrz))
      ENDIF
C
C Set output matrix size
      jpssize1 = size(kmat,1)
      jpssize2 = size(kmat,2)
      jpssize3 = size(kmat,3)
C
C -1.- Open NetCDF file
C ---------------------
C
      INQUIRE (FILE=kfnincrz,EXIST=filexists)
      IF (.NOT.filexists) GOTO 101
C
      ierr = NF90_OPEN(kfnincrz,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
C
C -2.- Read and check matrix dimension
C ------------------------------------
C Get variable id
      ierr = NF90_INQ_VARID(idf,kmatnam,idv)
      IF (ierr.NE.0) GOTO 104
C Get variable number dimensions
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,ndims=ndims)
      IF (ierr.NE.0) GOTO 104
C Check if it is a 3D matrix
      IF (ndims.NE.3) GOTO 104
C Get dimension id
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,dimids=idims)
      IF (ierr.NE.0) GOTO 104
C Get matrix size
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(1),len=vcount(1))
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(2),len=vcount(2))
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(3),len=vcount(3))
      IF (ierr.NE.0) GOTO 103
C Check if it is equal to output matrix size
      IF (vcount(1).NE.jpssize1) GOTO 105
      IF (vcount(2).NE.jpssize2) GOTO 105
      IF (vcount(3).NE.jpssize3) GOTO 105
C
C -3.- Read matrix file
C ---------------------
C
      vstart=(/1,1,1/)
      allocate ( ptabo(1:jpssize1,1:jpssize2,1:jpssize3), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
      ierr = NF90_GET_VAR(idf,idv,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
C
      kmat(:,:,:) = FREAL(ptabo(:,:,:))
C
C --- deallocation
      IF (allocated(ptabo)) deallocate(ptabo)
C
C -4.- Close observation file
C ---------------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrz','readmat3crz')
 1001 CALL printerror2(0,1001,3,'liocrz','readmat3crz')
C
 101  WRITE (texterror,*) 'File does not exist',kfnincrz
      CALL printerror2(0,101,3,'liocrz','readmat3crz',comment=texterror)
 102  WRITE (texterror,*) 'Bad NetCDF file: ',kfnincrz
      CALL printerror2(0,102,3,'liocrz','readmat3crz',comment=texterror)
 103  WRITE (texterror,*) 'Error reading dimension in file: ',kfnincrz
      CALL printerror2(0,103,3,'liocrz','readmat3crz',comment=texterror)
 104  WRITE (texterror,*) 'Error reading variable in file: ',kfnincrz
      CALL printerror2(0,104,3,'liocrz','readmat3crz',comment=texterror)
 105  WRITE (texterror,*) 'Inconsistent dimensions in file: ',kfnincrz
      CALL printerror2(0,105,3,'liocrz','readmat3crz',comment=texterror)
 106  WRITE (texterror,*) 'Error reading matrix in file: ',kfnincrz
      CALL printerror2(0,106,3,'liocrz','readmat3crz',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfnincrz
      CALL printerror2(0,107,3,'liocrz','readmat3crz',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readmat4crz(kfnincrz,kmatnam,kmat)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read 4D matrix in '.crz' file
CCC  -------
CCC  Method : Read matrix in NetCDF file
CCC  ------
CCC  Input :  kfnincrz    : filename
CCC  -----    kmatnam    : variable name
CCC
CCC  Output : kmat        : matrix
CCC  ------   
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrz, kmatnam
      BIGREAL, dimension(:,:,:,:), intent(out) :: kmat
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpssize1,jpssize2,jpssize3,jpssize4
      CHARACTER(len=word80) :: kform
      BIGREAL4, dimension(:,:,:,:), allocatable :: ptabo
C
      INTEGER :: ierr, idf, idv, ndims
      INTEGER, dimension(4) :: vstart,vcount,idims
      LOGICAL :: filexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readmat4crz'
         WRITE(numout,*) '    ==> READING file ',
     $                        kfnincrz(1:lenv(kfnincrz))
      ENDIF
C
C Set output matrix size
      jpssize1 = size(kmat,1)
      jpssize2 = size(kmat,2)
      jpssize3 = size(kmat,3)
      jpssize4 = size(kmat,4)
C
C -1.- Open NetCDF file
C ---------------------
C
      INQUIRE (FILE=kfnincrz,EXIST=filexists)
      IF (.NOT.filexists) GOTO 101
C
      ierr = NF90_OPEN(kfnincrz,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
C
C -2.- Read and check matrix dimension
C ------------------------------------
C Get variable id
      ierr = NF90_INQ_VARID(idf,kmatnam,idv)
      IF (ierr.NE.0) GOTO 104
C Get variable number dimensions
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,ndims=ndims)
      IF (ierr.NE.0) GOTO 104
C Check if it is a 4D matrix
      IF (ndims.NE.4) GOTO 104
C Get dimension id
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,dimids=idims)
      IF (ierr.NE.0) GOTO 104
C Get matrix size
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(1),len=vcount(1))
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(2),len=vcount(2))
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(3),len=vcount(3))
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(4),len=vcount(4))
      IF (ierr.NE.0) GOTO 103
C Check if it is equal to output matrix size
      IF (vcount(1).NE.jpssize1) GOTO 105
      IF (vcount(2).NE.jpssize2) GOTO 105
      IF (vcount(3).NE.jpssize3) GOTO 105
      IF (vcount(4).NE.jpssize4) GOTO 105
C
C -3.- Read matrix file
C ---------------------
C
      vstart=(/1,1,1,1/)
      allocate ( ptabo(1:jpssize1,1:jpssize2,1:jpssize3,1:jpssize4),
     $           stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
      ierr = NF90_GET_VAR(idf,idv,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
C
      kmat(:,:,:,:) = FREAL(ptabo(:,:,:,:))
C
C --- deallocation
      IF (allocated(ptabo)) deallocate(ptabo)
C
C -4.- Close observation file
C ---------------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrz','readmat4crz')
 1001 CALL printerror2(0,1001,3,'liocrz','readmat4crz')
C
 101  WRITE (texterror,*) 'File does not exist',kfnincrz
      CALL printerror2(0,101,3,'liocrz','readmat4crz',comment=texterror)
 102  WRITE (texterror,*) 'Bad NetCDF file: ',kfnincrz
      CALL printerror2(0,102,3,'liocrz','readmat4crz',comment=texterror)
 103  WRITE (texterror,*) 'Error reading dimension in file: ',kfnincrz
      CALL printerror2(0,103,3,'liocrz','readmat4crz',comment=texterror)
 104  WRITE (texterror,*) 'Error reading variable in file: ',kfnincrz
      CALL printerror2(0,104,3,'liocrz','readmat4crz',comment=texterror)
 105  WRITE (texterror,*) 'Inconsistent dimensions in file: ',kfnincrz
      CALL printerror2(0,105,3,'liocrz','readmat4crz',comment=texterror)
 106  WRITE (texterror,*) 'Error reading matrix in file: ',kfnincrz
      CALL printerror2(0,106,3,'liocrz','readmat4crz',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfnincrz
      CALL printerror2(0,107,3,'liocrz','readmat4crz',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writehdrcrz(kfnincrz,jpr,jpz,jpiobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write .crz' file header
CCC  -------
CCC  Method : Write NetCDF dimensions
CCC  ------
CCC  Input :  kfnincrz : filename
CCC  -----    jpr   : dimension of the reduced space
CCC           jpz   : number of local subsystems
CCC           jpiobs: number of segments in the observation vector
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrz
      INTEGER, intent(in) :: jpr, jpz, jpiobs
      LOGICAL :: filexists
      INTEGER :: ierr, idf, idr1, idr2, idz, idiobs
      INTEGER :: jpr1, jpr2, jpz1, jpiobs1
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writehdrcrz'
         WRITE(numout,*) '    ==> WRITING file ',
     $                           kfnincrz(1:lenv(kfnincrz))
      ENDIF
C
C -1.- Check If file exists
C -------------------------
C
      INQUIRE (FILE=kfnincrz,EXIST=filexists)
      IF (.NOT.filexists) THEN
C
C -2.- If file does not exist, create the file and dimensions
C -----------------------------------------------------------
C
        ierr = NF90_CREATE(kfnincrz,NF90_CLOBBER,idf)
        IF (ierr.NE.0) GOTO 101
        ierr = NF90_DEF_DIM(idf,'r1',jpr,idr1)
        IF (ierr.NE.0) GOTO 105
        ierr = NF90_DEF_DIM(idf,'r2',jpr,idr2)
        IF (ierr.NE.0) GOTO 105
        ierr = NF90_DEF_DIM(idf,'z',jpz,idz)
        IF (ierr.NE.0) GOTO 105
        ierr = NF90_DEF_DIM(idf,'iobs',jpiobs,idiobs)
        IF (ierr.NE.0) GOTO 105
        ierr = NF90_ENDDEF(idf)
        IF (ierr.NE.0) GOTO 101
C
      ELSE
C
C -3.- Else, read dimensions and check that they are consistent
C -------------------------------------------------------------
C
        ierr = NF90_OPEN(kfnincrz,NF90_NOWRITE,idf)
        IF (ierr.NE.0) GOTO 102

        ierr = NF90_INQ_DIMID(idf,'r1',idr1)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_INQ_DIMID(idf,'r2',idr2)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_INQ_DIMID(idf,'z',idz)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_INQ_DIMID(idf,'iobs',idiobs)
        IF (ierr.NE.0) GOTO 103
C
        ierr = NF90_INQUIRE_DIMENSION(idf,idr1,len=jpr1)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_INQUIRE_DIMENSION(idf,idr2,len=jpr2)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=jpz1)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=jpiobs1)
        IF (ierr.NE.0) GOTO 103
C
        IF (jpr1.NE.jpr) GOTO 104
        IF (jpr2.NE.jpr) GOTO 104
        IF (jpz1.NE.jpz) GOTO 104
        IF (jpiobs1.NE.jpiobs) GOTO 104
C
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrz','writehdrcrz')
C
 101  WRITE (texterror,*) 'Error creating Netcdf file: ',kfnincrz
      CALL printerror2(0,101,1,'liocrz','writehdrcrz',comment=texterror)
 102  WRITE (texterror,*) 'Bad NetCDF file: ',kfnincrz
      CALL printerror2(0,102,3,'liocrz','writehdrcrz',comment=texterror)
 103  WRITE (texterror,*) 'Error reading dimensions in file: ',kfnincrz
      CALL printerror2(0,103,3,'liocrz','writehdrcrz',comment=texterror)
 104  WRITE (texterror,*) 'Inconsistent .crz file: ',kfnincrz
      CALL printerror2(0,104,3,'liocrz','writehdrcrz',comment=texterror)
 105  WRITE (texterror,*) 'Error creating Netcdf dimension: ',kfnincrz
      CALL printerror2(0,103,1,'liocrz','writehdrcrz',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeveccrz(kfnincrz,kvectnam,kvect,ktyp)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write vector in '.crz' file
CCC  -------
CCC  Method : Write vector in NetCDF file
CCC  ------
CCC  Input : kfnincrz    : filename
CCC  -----   kvectnam    : variable name
CCC          kvect       : vector
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrz,kvectnam,ktyp
      BIGREAL, dimension(:), intent(in) :: kvect
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok, jpssize
      BIGREAL4, dimension(:), allocatable :: ptabo
C
      INTEGER :: ierr, idf, idv
      INTEGER, dimension(1) :: vstart, vcount, idims, idims1
      LOGICAL :: varexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writeveccrz'
         WRITE(numout,*) '    ==> WRITING file ',
     $                           kfnincrz(1:lenv(kfnincrz))
      ENDIF
C
C Set output vector size
      jpssize = size(kvect,1)
C
C -1.- Set variable dimension list
C --------------------------------
C
      ierr = NF90_OPEN(kfnincrz,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
C
      SELECT CASE (ktyp)
      CASE('r')
        ierr = NF90_INQ_DIMID(idf,'r1',idims(1))
        IF (ierr.NE.0) GOTO 101
      CASE('z')
        ierr = NF90_INQ_DIMID(idf,'z',idims(1))
        IF (ierr.NE.0) GOTO 101
      CASE DEFAULT
        GOTO 1000
      END SELECT
C
C -2.- Get dimensions ids and check that they are consistent
C ----------------------------------------------------------
C
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(1),len=vcount(1))
      IF (ierr.NE.0) GOTO 101
      IF (vcount(1).NE.jpssize) GOTO 108
C
C -3.- Check if the variable exists in the NetCDF file
C ----------------------------------------------------
C
      ierr = NF90_OPEN(kfnincrz,NF90_WRITE,idf)
      IF (ierr.NE.0) GOTO 102
      ierr = NF90_INQ_VARID(idf,kvectnam,idv)
      varexists = (ierr.EQ.0)
C
      IF (.NOT.varexists) THEN
C
C -4.- If variable does not exist, create it
C ------------------------------------------
C
        ierr = NF90_REDEF(idf)
        IF (ierr.NE.0) GOTO 104
        ierr = NF90_DEF_VAR(idf,kvectnam,NF90_FLOAT,idims,idv)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_ENDDEF(idf)
        IF (ierr.NE.0) GOTO 104
C
      ENDIF
C
C -5.- Write vector in file
C -------------------------
C
      vstart=1
      allocate ( ptabo(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
      ptabo(:) = FREAL4(kvect(:))
      ierr = NF90_PUT_VAR(idf,idv,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
C
C -6.- Close observation file
C ---------------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrz','writeveccrz')
 1001 CALL printerror2(0,1001,3,'liocrz','writeveccrz')
C
 101  WRITE (texterror,*) 'Incorrect Netcdf dimension: ',kfnincrz
      CALL printerror2(0,101,3,'liocrz','writeveccrz',comment=texterror)
 102  WRITE (texterror,*) 'Error opening Netcdf file: ',kfnincrz
      CALL printerror2(0,102,1,'liocrz','writeveccrz',comment=texterror)
 103  WRITE (texterror,*) 'Error creating Netcdf variable: ',kfnincrz
      CALL printerror2(0,103,1,'liocrz','writeveccrz',comment=texterror)
 104  WRITE (texterror,*) 'Error with Netcdf file: ',kfnincrz
      CALL printerror2(0,104,1,'liocrz','writeveccrz',comment=texterror)
 105  WRITE (texterror,*) 'Incorrect Netcdf variable: ',kfnincrz
      CALL printerror2(0,105,1,'liocrz','writeveccrz',comment=texterror)
 106  WRITE (texterror,*) 'Error writing vector in file: ',kfnincrz
      CALL printerror2(0,106,1,'liocrz','writeveccrz',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfnincrz
      CALL printerror2(0,107,1,'liocrz','writeveccrz',comment=texterror)
 108  WRITE (texterror,*) 'Inconsistent dimensions in file: ',kfnincrz
      CALL printerror2(0,108,3,'liocrz','writeveccrz',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writemat2crz(kfnincrz,kmatnam,kmat,ktyp)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write 2D matrix in '.crz' file
CCC  -------
CCC  Method : Write matrix in NetCDF file
CCC  ------
CCC  Input : kfnincrz   : filename
CCC  -----   kmatnam    : variable name
CCC          kmat       : matrix
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrz,kmatnam,ktyp
      BIGREAL, dimension(:,:), intent(in) :: kmat
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok, jpssize1, jpssize2
      BIGREAL4, dimension(:,:), allocatable :: ptabo
C
      INTEGER :: ierr, idf, idv
      INTEGER, dimension(2) :: vstart, vcount, idims, idims1
      LOGICAL :: varexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writemat2crz'
         WRITE(numout,*) '    ==> WRITING file ',
     $                           kfnincrz(1:lenv(kfnincrz))
      ENDIF
C
C Set output matrix size
      jpssize1 = size(kmat,1)
      jpssize2 = size(kmat,2)
C
C -1.- Set variable dimension list
C --------------------------------
C
      ierr = NF90_OPEN(kfnincrz,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
C
      SELECT CASE (ktyp)
      CASE('rr')
        ierr = NF90_INQ_DIMID(idf,'r1',idims(1))
        IF (ierr.NE.0) GOTO 101
        ierr = NF90_INQ_DIMID(idf,'r2',idims(2))
        IF (ierr.NE.0) GOTO 101
      CASE('rz')
        ierr = NF90_INQ_DIMID(idf,'r1',idims(1))
        IF (ierr.NE.0) GOTO 101
        ierr = NF90_INQ_DIMID(idf,'z',idims(2))
        IF (ierr.NE.0) GOTO 101
      CASE('zi')
        ierr = NF90_INQ_DIMID(idf,'z',idims(1))
        IF (ierr.NE.0) GOTO 101
        ierr = NF90_INQ_DIMID(idf,'iobs',idims(2))
        IF (ierr.NE.0) GOTO 101
      CASE DEFAULT
        GOTO 1000
      END SELECT
C
C -2.- Get dimensions sizes and check that they are consistent
C ------------------------------------------------------------
C
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(1),len=vcount(1))
      IF (ierr.NE.0) GOTO 101
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(2),len=vcount(2))
      IF (ierr.NE.0) GOTO 101
      IF (vcount(1).NE.jpssize1) GOTO 108
      IF (vcount(2).NE.jpssize2) GOTO 108
C
C -3.- Check if the variable exists in the NetCDF file
C ----------------------------------------------------
C
      ierr = NF90_OPEN(kfnincrz,NF90_WRITE,idf)
      IF (ierr.NE.0) GOTO 102
      ierr = NF90_INQ_VARID(idf,kmatnam,idv)
      varexists = (ierr.EQ.0)
C
      IF (.NOT.varexists) THEN
C
C -4.- If variable does not exist, create it
C ------------------------------------------
C
        ierr = NF90_REDEF(idf)
        IF (ierr.NE.0) GOTO 104
        ierr = NF90_DEF_VAR(idf,kmatnam,NF90_FLOAT,idims,idv)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_ENDDEF(idf)
        IF (ierr.NE.0) GOTO 104
C
      ENDIF
C
C -5.- Write matrix in file
C -------------------------
C
      vstart=(/1,1/)
      allocate ( ptabo(1:jpssize1,1:jpssize2), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
      ptabo(:,:) = FREAL4(kmat(:,:))
      ierr = NF90_PUT_VAR(idf,idv,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
C
C -6.- Close observation file
C ---------------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrz','writemat2crz')
 1001 CALL printerror2(0,1001,3,'liocrz','writemat2crz')
C
 101  WRITE (texterror,*) 'Incorrect Netcdf dimension: ',kfnincrz
      CALL printerror2(0,101,3,'liocrz','writemat2crz',comment=texterror)
 102  WRITE (texterror,*) 'Error opening Netcdf file: ',kfnincrz
      CALL printerror2(0,102,1,'liocrz','writemat2crz',comment=texterror)
 103  WRITE (texterror,*) 'Error creating Netcdf variable: ',kfnincrz
      CALL printerror2(0,103,1,'liocrz','writemat2crz',comment=texterror)
 104  WRITE (texterror,*) 'Error with Netcdf file: ',kfnincrz
      CALL printerror2(0,104,1,'liocrz','writemat2crz',comment=texterror)
 105  WRITE (texterror,*) 'Incorrect Netcdf variable: ',kfnincrz
      CALL printerror2(0,105,1,'liocrz','writemat2crz',comment=texterror)
 106  WRITE (texterror,*) 'Error writing vector in file: ',kfnincrz
      CALL printerror2(0,106,1,'liocrz','writemat2crz',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfnincrz
      CALL printerror2(0,107,1,'liocrz','writemat2crz',comment=texterror)
 108  WRITE (texterror,*) 'Inconsistent dimensions in file: ',kfnincrz
      CALL printerror2(0,108,3,'liocrz','writemat2crz',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writemat3crz(kfnincrz,kmatnam,kmat,ktyp)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write 3D matrix in '.crz' file
CCC  -------
CCC  Method : Write matrix in NetCDF file
CCC  ------
CCC  Input : kfnincrz   : filename
CCC  -----   kmatnam    : variable name
CCC          kmat       : matrix
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrz,kmatnam,ktyp
      BIGREAL, dimension(:,:,:), intent(in) :: kmat
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok, jpssize1, jpssize2, jpssize3
      BIGREAL4, dimension(:,:,:), allocatable :: ptabo
C
      INTEGER :: ierr, idf, idv
      INTEGER, dimension(3) :: vstart, vcount, idims, idims1
      LOGICAL :: varexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writemat3crz'
         WRITE(numout,*) '    ==> WRITING file ',
     $                           kfnincrz(1:lenv(kfnincrz))
      ENDIF
C
C Set output matrix size
      jpssize1 = size(kmat,1)
      jpssize2 = size(kmat,2)
      jpssize3 = size(kmat,3)
C
C -1.- Set variable dimension list
C --------------------------------
C
      ierr = NF90_OPEN(kfnincrz,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
C
      SELECT CASE (ktyp)
      CASE('rrz')
        ierr = NF90_INQ_DIMID(idf,'r1',idims(1))
        IF (ierr.NE.0) GOTO 101
        ierr = NF90_INQ_DIMID(idf,'r2',idims(2))
        IF (ierr.NE.0) GOTO 101
        ierr = NF90_INQ_DIMID(idf,'z',idims(3))
        IF (ierr.NE.0) GOTO 101
      CASE('rzi')
        ierr = NF90_INQ_DIMID(idf,'r1',idims(1))
        IF (ierr.NE.0) GOTO 101
        ierr = NF90_INQ_DIMID(idf,'z',idims(2))
        IF (ierr.NE.0) GOTO 101
        ierr = NF90_INQ_DIMID(idf,'iobs',idims(3))
        IF (ierr.NE.0) GOTO 101
      CASE DEFAULT
        GOTO 1000
      END SELECT
C
C -2.- Get dimensions sizes and check that they are consistent
C ------------------------------------------------------------
C
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(1),len=vcount(1))
      IF (ierr.NE.0) GOTO 101
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(2),len=vcount(2))
      IF (ierr.NE.0) GOTO 101
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(3),len=vcount(3))
      IF (ierr.NE.0) GOTO 101
      IF (vcount(1).NE.jpssize1) GOTO 108
      IF (vcount(2).NE.jpssize2) GOTO 108
      IF (vcount(3).NE.jpssize3) GOTO 108
C
C -3.- Check if the variable exists in the NetCDF file
C ----------------------------------------------------
C
      ierr = NF90_OPEN(kfnincrz,NF90_WRITE,idf)
      IF (ierr.NE.0) GOTO 102
      ierr = NF90_INQ_VARID(idf,kmatnam,idv)
      varexists = (ierr.EQ.0)
C
      IF (.NOT.varexists) THEN
C
C -4.- If variable does not exist, create it
C ------------------------------------------
C
        ierr = NF90_REDEF(idf)
        IF (ierr.NE.0) GOTO 104
        ierr = NF90_DEF_VAR(idf,kmatnam,NF90_FLOAT,idims,idv)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_ENDDEF(idf)
        IF (ierr.NE.0) GOTO 104
C
      ENDIF
C
C -5.- Write matrix in file
C -------------------------
C
      vstart=(/1,1,1/)
      allocate ( ptabo(1:jpssize1,1:jpssize2,1:jpssize3), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
      ptabo(:,:,:) = FREAL4(kmat(:,:,:))
      ierr = NF90_PUT_VAR(idf,idv,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
C
C -6.- Close observation file
C ---------------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrz','writemat3crz')
 1001 CALL printerror2(0,1001,3,'liocrz','writemat3crz')
C
 101  WRITE (texterror,*) 'Incorrect Netcdf dimension: ',kfnincrz
      CALL printerror2(0,101,3,'liocrz','writemat3crz',comment=texterror)
 102  WRITE (texterror,*) 'Error opening Netcdf file: ',kfnincrz
      CALL printerror2(0,102,1,'liocrz','writemat3crz',comment=texterror)
 103  WRITE (texterror,*) 'Error creating Netcdf variable: ',kfnincrz
      CALL printerror2(0,103,1,'liocrz','writemat3crz',comment=texterror)
 104  WRITE (texterror,*) 'Error with Netcdf file: ',kfnincrz
      CALL printerror2(0,104,1,'liocrz','writemat3crz',comment=texterror)
 105  WRITE (texterror,*) 'Incorrect Netcdf variable: ',kfnincrz
      CALL printerror2(0,105,1,'liocrz','writemat3crz',comment=texterror)
 106  WRITE (texterror,*) 'Error writing vector in file: ',kfnincrz
      CALL printerror2(0,106,1,'liocrz','writemat3crz',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfnincrz
      CALL printerror2(0,107,1,'liocrz','writemat3crz',comment=texterror)
 108  WRITE (texterror,*) 'Inconsistent dimensions in file: ',kfnincrz
      CALL printerror2(0,108,3,'liocrz','writemat3crz',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writemat4crz(kfnincrz,kmatnam,kmat,ktyp)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write 4D matrix in '.crz' file
CCC  -------
CCC  Method : Write matrix in NetCDF file
CCC  ------
CCC  Input : kfnincrz   : filename
CCC  -----   kmatnam    : variable name
CCC          kmat       : matrix
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrz,kmatnam,ktyp
      BIGREAL, dimension(:,:,:,:), intent(in) :: kmat
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok, jpssize1, jpssize2, jpssize3, jpssize4
      BIGREAL4, dimension(:,:,:,:), allocatable :: ptabo
C
      INTEGER :: ierr, idf, idv
      INTEGER, dimension(4) :: vstart, vcount, idims, idims1
      LOGICAL :: varexists
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writemat4crz'
         WRITE(numout,*) '    ==> WRITING file ',
     $                           kfnincrz(1:lenv(kfnincrz))
      ENDIF
C
C Set output matrix size
      jpssize1 = size(kmat,1)
      jpssize2 = size(kmat,2)
      jpssize3 = size(kmat,3)
      jpssize4 = size(kmat,4)
C
C -1.- Set variable dimension list
C --------------------------------
C
      ierr = NF90_OPEN(kfnincrz,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 102
C
      SELECT CASE (ktyp)
      CASE('rrzi')
        ierr = NF90_INQ_DIMID(idf,'r1',idims(1))
        IF (ierr.NE.0) GOTO 101
        ierr = NF90_INQ_DIMID(idf,'r2',idims(2))
        IF (ierr.NE.0) GOTO 101
        ierr = NF90_INQ_DIMID(idf,'z',idims(3))
        IF (ierr.NE.0) GOTO 101
        ierr = NF90_INQ_DIMID(idf,'iobs',idims(4))
        IF (ierr.NE.0) GOTO 101
      CASE DEFAULT
        GOTO 1000
      END SELECT
C
C -2.- Get dimensions sizes and check that they are consistent
C ------------------------------------------------------------
C
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(1),len=vcount(1))
      IF (ierr.NE.0) GOTO 101
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(2),len=vcount(2))
      IF (ierr.NE.0) GOTO 101
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(3),len=vcount(3))
      IF (ierr.NE.0) GOTO 101
      ierr = NF90_INQUIRE_DIMENSION(idf,idims(4),len=vcount(4))
      IF (ierr.NE.0) GOTO 101
      IF (vcount(1).NE.jpssize1) GOTO 108
      IF (vcount(2).NE.jpssize2) GOTO 108
      IF (vcount(3).NE.jpssize3) GOTO 108
      IF (vcount(4).NE.jpssize4) GOTO 108
C
C -3.- Check if the variable exists in the NetCDF file
C ----------------------------------------------------
C
      ierr = NF90_OPEN(kfnincrz,NF90_WRITE,idf)
      IF (ierr.NE.0) GOTO 102
      ierr = NF90_INQ_VARID(idf,kmatnam,idv)
      varexists = (ierr.EQ.0)
C
      IF (.NOT.varexists) THEN
C
C -4.- If variable does not exist, create it
C ------------------------------------------
C
        ierr = NF90_REDEF(idf)
        IF (ierr.NE.0) GOTO 104
        ierr = NF90_DEF_VAR(idf,kmatnam,NF90_FLOAT,idims,idv)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_ENDDEF(idf)
        IF (ierr.NE.0) GOTO 104
C
      ENDIF
C
C -5.- Write matrix in file
C -------------------------
C
      vstart=(/1,1,1,1/)
      allocate ( ptabo(1:jpssize1,1:jpssize2,1:jpssize3,1:jpssize4),
     $           stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
      ptabo(:,:,:,:) = FREAL4(kmat(:,:,:,:))
      ierr = NF90_PUT_VAR(idf,idv,ptabo,start=vstart,count=vcount)
      IF (ierr.NE.0) GOTO 106
C
C -6.- Close observation file
C ---------------------------
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
      RETURN
C
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrz','writemat4crz')
 1001 CALL printerror2(0,1001,3,'liocrz','writemat4crz')
C
 101  WRITE (texterror,*) 'Incorrect Netcdf dimension: ',kfnincrz
      CALL printerror2(0,101,3,'liocrz','writemat4crz',comment=texterror)
 102  WRITE (texterror,*) 'Error opening Netcdf file: ',kfnincrz
      CALL printerror2(0,102,1,'liocrz','writemat4crz',comment=texterror)
 103  WRITE (texterror,*) 'Error creating Netcdf variable: ',kfnincrz
      CALL printerror2(0,103,1,'liocrz','writemat4crz',comment=texterror)
 104  WRITE (texterror,*) 'Error with Netcdf file: ',kfnincrz
      CALL printerror2(0,104,1,'liocrz','writemat4crz',comment=texterror)
 105  WRITE (texterror,*) 'Incorrect Netcdf variable: ',kfnincrz
      CALL printerror2(0,105,1,'liocrz','writemat4crz',comment=texterror)
 106  WRITE (texterror,*) 'Error writing vector in file: ',kfnincrz
      CALL printerror2(0,106,1,'liocrz','writemat4crz',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfnincrz
      CALL printerror2(0,107,1,'liocrz','writemat4crz',comment=texterror)
 108  WRITE (texterror,*) 'Inconsistent dimensions in file: ',kfnincrz
      CALL printerror2(0,108,3,'liocrz','writemat4crz',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liocrz
