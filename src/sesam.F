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
C ---                  SESAM.F                                  ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12  (C.E. Testut)                       ---
C --- modification : 99-05  (C.E. Testut)                       ---
C --- modification : 01-06  (C.E. Testut)                       ---
C --- modification : 03-01  (J.M. Brankart)                     ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE sesam
C --- 
C -----------------------------------------------------------------
CCC*********************************************************************
CCC*********************************************************************
CCC
CCC                 Institut de Mecanique de Grenoble
CCC
CCC                          L E G I / I M G
CCC       Laboratoire des Ecoulements Geophysiques et Industriels
CCC                           Equipe MEOM
CCC               Modelisation des Ecoulements Oceaniques
CCC                    de Moyenne et grande Echelle
CCC  
CCC*********************************************************************
CCC 
CCC       
CCC                            S E S A M :
CCC
CCC       An integrated System of Sequential assimilation modules
CCC   
CCC   
CCC Authors  :  Charles-Emmanuel Testut
CCC             Jean-Michel Brankart
CCC             Laurent Parent
CCC  
CCC Original version : December 1997
CCC
CCC Version SESAM 3.2 : January 2003
CCC Version SESAM 4.1 : March 2009
CCC Version SESAM 4.3 : May 2011
CCC
CCC*********************************************************************
      SUBROUTINE sesam(textcommand)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Integrated System of Sequential assimilation modules
CCC  -------
CCC  Method :   1) Define SESAM constants (list of possible actions,
CCC  ------           list of possible switches, list of possible object
CCC                   formats, help information, ...)
CCC             2) Load SESAM configuration (defined by the user
CCC                   in 'defcst.control.h').
CCC             3) Read and interpret user command to decide which
CCC                   module to invoke and which action to perform
CCC             4) Read SESAM namelist and SESAM mask files.
CCC             5) Check consistency of SESAM parameters
CCC             6) Run required SESAM module
CCC   
CCC  Input : User command, provided either as an argument of the
CCC  -----      SESAM routine, or as commandline arguments.
CCC
CCC  Output : no
CCC  ------
CCC  References : SESAM 3.2 reference manual.
CCC  ----------
CCC---------------------------------------------------------------------
CCC                                      
CCC  main SESAM routine calling graph
CCC  --------------------------------
CCC   sesam                 main SESAM routine
CCC   ==>  defcst           define SESAM constants and load user configuration
CCC   ==>  readarg          interpret user command to know which module
CCC                            to invoke and which action to perform
CCC   ==>  readmsk          read SESAM mask files
CCC   ==>  evalconfig       set up SESAM internal configuration
CCC   ==>  modintf (01)     module 01: interface an object from a file
CCC                            format to another
CCC   ==>  modcorr (02)     module 02: compute correlation coefficients
CCC   ==>  modtgop (03)     module 03: truncated Gaussian module
CCC   ==>  modobsv (04)     module 04: observation management
CCC   ==>  modfilt (05)     module 05: low-pass filter
CCC   ==>  modadap (06)     module 06: adaptive parameter estimation
CCC   ==>  moddiff (07)     module 07: compute RMS difference between
CCC                            vector objects
CCC   ==>  modoerr (08)     module 08: observation error management
CCC   ==>  modoper (09)     module 09: arithmetic operation
CCC   ==>  modgroa (10)     module 10: global EOF decomposition
CCC   ==>  modlroa (11)     module 11: local EOF decomposition
CCC   ==>  modbroa (12)     module 12: bubble EOF decomposition
CCC   ==>  modzone (13)     module 13: local data section management
CCC   ==>  modgeof (14)     module 14: global reduced order analysis
CCC   ==>  modleof (15)     module 15: local reduced order analysis
CCC   ==>  modbeof (16)     module 16: bubble reduced order analysis
CCC   ==>  modgreg (17)     module 17: global linear regression
CCC   ==>  modlreg (18)     module 18: local linear regression
CCC   ==>  modbreg (19)     module 19: bubble linear regression
CCC   ==>  modvari (20)     module 20: compute variances
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mpitime
      use argsesam
      use hiomsk
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*) :: textcommand
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=13) :: cmd
      BIGREAL, parameter :: xprec=1.0_kr
      BIGREAL8, parameter :: xprec8=1.0_kr8
      BIGREAL4, parameter :: xprec4=1.0_kr4
      TYPE (type_swiarg), dimension(1:nbargmax) :: execaction
      INTEGER :: jloop,jindbeg,jindend
#if defined GETARG
      EXTERNAL getarg
#endif
CC----------------------------------------------------------------------
C
      nprint = 0
C
C - Initialisation of the parallel code
#if defined MPI
      CALL mpi_init(mpi_code)
      CALL mpi_comm_size(mpi_comm_world,jpproc,mpi_code)
      CALL mpi_comm_rank(mpi_comm_world,jproc,mpi_code)
      call MPI_TIMER(0)
#endif
C
C -0.- Get user command from input
C --------------------------------
C
      execaction(:) = type_swiarg(' ',' ')
      IF (VERIFY(textcommand,' ').NE.0) THEN
C  get user command from routine argument
         jindbeg=VERIFY(textcommand(jindend+1:),' ')
         jindend=SCAN(textcommand(jindbeg:),' ')-1
         argloop1: DO jloop = 1, nbargmax
            execaction(jloop)%swi=textcommand(jindbeg:jindend) 
            jindbeg=VERIFY(textcommand(jindend+1:),' ')
            IF (jindbeg.EQ.0) EXIT argloop1
            jindbeg=jindbeg+jindend
            jindend=SCAN(textcommand(jindbeg:),' ')-1
            jindend=jindend+jindbeg
            execaction(jloop)%arg=textcommand(jindbeg:jindend)   
            jindbeg=VERIFY(textcommand(jindend+1:),' ')
            IF (jindbeg.EQ.0) EXIT argloop1
            jindbeg=jindbeg+jindend
            jindend=SCAN(textcommand(jindbeg:),' ')-1
            jindend=jindend+jindbeg
         END DO argloop1  
      ELSE
C  get user command from commandline arguments
         argloop: DO jloop = 1, nbargmax
#if defined GETARG
           CALL getarg ((jloop*2-1),execaction(jloop)%swi)
           CALL getarg ((jloop*2),execaction(jloop)%arg)
#else
           CALL get_command_argument ((jloop*2-1),execaction(jloop)%swi)
           CALL get_command_argument ((jloop*2),execaction(jloop)%arg)
#endif
           IF (execaction(jloop)%swi.EQ.' ') EXIT argloop
         END DO argloop
      ENDIF
C
C -1.- Open SESAM listing file and SESAM configuration file
C ---------------------------------------------------------
C
      DO jloop = 1, nbargmax
         IF (execaction(jloop)%swi(1:MAX(1,lenv(execaction(jloop)%swi)))
     $        .EQ.'-outinfo') fnamout=execaction(jloop)%arg
      ENDDO
      CALL openfile(numout,fnamout,kstatus=clunk)
C
      DO jloop = 1, nbargmax
         IF (execaction(jloop)%swi(1:MAX(1,lenv(execaction(jloop)%swi)))
     $        .EQ.'-config')    fnamlist=execaction(jloop)%arg
      ENDDO
      CALL openfile(numnam,fnamlist,kstatus=clunk)
C
      IF ((nprint.GE.0).AND.(jproc.EQ.0)) THEN
         WRITE(numout,*)
         WRITE(numout,*) repeat('@',59)
         WRITE(numout,*) '@'
         WRITE(numout,*) '@             S E S A M   Software'
         WRITE(numout,*) '@'
         WRITE(numout,*) '@   SystEm of Sequential Assimilation Modules'
         WRITE(numout,*) '@ '
         WRITE(numout,*) '@'
         WRITE(numout,*) '@ Equipe MEOM'
         WRITE(numout,*) '@ LEGI/IMG'
         WRITE(numout,*) '@ Grenoble'
         WRITE(numout,*) '@'
         WRITE(numout,*) '@ version ',version(:lenv(version))
         WRITE(numout,*) '@ (',version_year(:lenv(version_year)),')'
         WRITE(numout,*) '@'
         WRITE(numout,*) repeat('@',59)
      ENDIF
C
      IF ((nprint.GE.0).AND.(jproc.EQ.0)) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Check SESAM real variable types:'
         WRITE(numout,*) ' ================================'
         WRITE(numout,*) ' Precision kr :'
         WRITE(numout,*) ' --------------'
         WRITE(numout,*) ' kr              =',kr
         WRITE(numout,*) ' digits(xprec)   =',digits(xprec)
         WRITE(numout,*) ' epsilon(xprec)  =',epsilon(xprec)
         WRITE(numout,*) ' huge(xprec)     =',huge(xprec)
         WRITE(numout,*) ' precision(xprec)=',precision(xprec)
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Precision kr8:'
         WRITE(numout,*) ' --------------'
         WRITE(numout,*) ' kr8              =',kr8
         WRITE(numout,*) ' digits(xprec8)   =',digits(xprec8)
         WRITE(numout,*) ' epsilon(xprec8)  =',epsilon(xprec8)
         WRITE(numout,*) ' huge(xprec8)     =',huge(xprec8)
         WRITE(numout,*) ' precision(xprec8)=',precision(xprec8)
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Precision kr4:'
         WRITE(numout,*) ' --------------'
         WRITE(numout,*) ' kr4              =',kr4
         WRITE(numout,*) ' digits(xprec4)   =',digits(xprec4)
         WRITE(numout,*) ' epsilon(xprec4)  =',epsilon(xprec4)
         WRITE(numout,*) ' huge(xprec4)     =',huge(xprec4)
         WRITE(numout,*) ' precision(xprec4)=',precision(xprec4)
         WRITE(numout,*)
      ENDIF 
C
C -2.- Initialisation of SESAM parameters :
C -----------------------------------------
C
C -2.1- Define SESAM constants
C
      CALL defcst
C
C -2.2- Interpret user command to know which module to invoke
C
      CALL readarg(execaction)
C
C -2.3- Read SESAM masks
C
      CALL readmsk
C
C -2.4- Set up SESAM internal configuration
C
      CALL evalconfig
C
C -3.- Run appropriate SESAM module
C ---------------------------------
C
      SELECT CASE (nmode)
      CASE (1)
C -3.1- INTF mode : interface an object from a file format to another
C
         CALL modintf
C
      CASE (2)
C -3.2- ---- mode : compute correlation coefficients
C
         CALL modcorr
C
      CASE (3)
C -3.3- ---- mode : truncated Gaussian operations
C
         CALL modtgop
C
      CASE (4)
C -3.4- OBSV mode : observation management
C
         CALL modobsv
C
      CASE (5)
C -3.5- FILT mode : low-pass filter
C      
         CALL modfilt
C
      CASE (6)
C -3.6- ADAP mode : adaptive parameter estimation
C
         CALL modadap
C
      CASE (7)
C -3.7- DIFF mode : compute RMS difference between vector objects
C
         CALL moddiff
C
      CASE (8)
C -3.8- OERR mode : observation error management
C
         CALL modoerr
C
      CASE (9)
C -3.9- OPER mode : arithmetic operation
C
         CALL modoper
C
      CASE (10)
C -3.10- GEOF mode : global EOF decomposition
C
         CALL modglbeof
C
      CASE (11)
C -3.11- LEOF mode : local EOF decomposition
C
         CALL modglbeof
C
      CASE (12)
C -3.12- BEOF mode : bubble EOF decomposition
C
         CALL modglbeof
C
      CASE (13)
C -3.13- ZONE mode : local data section management
C
         CALL modzone
C
      CASE (14)
C -3.14- GROA mode : global reduced order analysis
C
         CALL modglbroa
C
      CASE (15)
C -3.15- LROA mode : local reduced order analysis
C
         CALL modglbroa
C
      CASE (16)
C -3.16- BROA mode : bubble reduced order analysis
C
         CALL modglbroa
C
      CASE (17)
C -3.17- GREG mode : global linear regression
C
         CALL modglbreg
C
      CASE (18)
C -3.18- LREG mode : local linear regression
C
         CALL modglbreg
C
      CASE (19)
C -3.19- BREG mode : bubble linear regression
C
         CALL modglbreg
C
      CASE (20)
C -3.20- VARI mode : variance computation
C
         CALL modvari
C
      CASE (21)
C -3.21- ANAM mode : anamorphosis
C
         CALL modanam
C
      CASE (22)
C -3.22- SCOR mode : probabilistic scores (CRPS,...)
C
         CALL modscor
C
      CASE (23)
C -3.23- SPCT mode : spectral transformation (spherical harmonics)
C
         CALL modspct
C
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -4.- Create empty 'SESAMOK' file 
C --------------------------------
C (to prove the job has terminated successfully)
C
      IF ((nprint.GE.0).AND.(jproc.EQ.0)) THEN
         cmd='touch SESAMOK'
         CALL shellorder(cmd)
      ENDIF
C
C -5.- This is the End
C --------------------
C
      IF ((nprint.GE.0).AND.(jproc.EQ.0)) THEN
         WRITE(numout,*)
         WRITE(numout,*) '@'
         WRITE(numout,*) '@ End of ',version(:lenv(version))
         WRITE(numout,*) repeat('@',59)
      ENDIF
C
#if defined MPI
      call MPI_TIMER(1)
      CALL mpi_finalize(mpi_code)
#endif
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'sesam','sesam')
C
      END
