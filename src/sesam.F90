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
! ---                  SESAM.F90                                  ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12  (C.E. Testut)                       ---
! --- modification : 99-05  (C.E. Testut)                       ---
! --- modification : 01-06  (C.E. Testut)                       ---
! --- modification : 03-01  (J.M. Brankart)                     ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE sesam
! --- 
! -----------------------------------------------------------------
!*********************************************************************
!*********************************************************************
!
!                 Institut de Mecanique de Grenoble
!
!                          L E G I / I M G
!       Laboratoire des Ecoulements Geophysiques et Industriels
!                           Equipe MEOM
!               Modelisation des Ecoulements Oceaniques
!                    de Moyenne et grande Echelle
!  
!*********************************************************************
! 
!       
!                            S E S A M :
!
!       An integrated System of Sequential assimilation modules
!   
!   
! Authors  :  Charles-Emmanuel Testut
!             Jean-Michel Brankart
!             Laurent Parent
!  
! Original version : December 1997
!
! Version SESAM 3.2 : January 2003
! Version SESAM 4.1 : March 2009
! Version SESAM 4.3 : May 2011
!
!*********************************************************************
      SUBROUTINE sesam(textcommand)
!---------------------------------------------------------------------
!
!  Purpose : Integrated System of Sequential assimilation modules
!  -------
!  Method :   1) Define SESAM constants (list of possible actions,
!  ------           list of possible switches, list of possible object
!                   formats, help information, ...)
!             2) Load SESAM configuration (defined by the user
!                   in 'defcst.control.h').
!             3) Read and interpret user command to decide which
!                   module to invoke and which action to perform
!             4) Read SESAM namelist and SESAM mask files.
!             5) Check consistency of SESAM parameters
!             6) Run required SESAM module
!   
!  Input : User command, provided either as an argument of the
!  -----      SESAM routine, or as commandline arguments.
!
!  Output : no
!  ------
!  References : SESAM 3.2 reference manual.
!  ----------
!---------------------------------------------------------------------
!                                      
!  main SESAM routine calling graph
!  --------------------------------
!   sesam                 main SESAM routine
!   ==>  defcst           define SESAM constants and load user configuration
!   ==>  readarg          interpret user command to know which module
!                            to invoke and which action to perform
!   ==>  readmsk          read SESAM mask files
!   ==>  evalconfig       set up SESAM internal configuration
!   ==>  modintf (01)     module 01: interface an object from a file
!                            format to another
!   ==>  modcorr (02)     module 02: compute correlation coefficients
!   ==>  modtgop (03)     module 03: truncated Gaussian module
!   ==>  modobsv (04)     module 04: observation management
!   ==>  modfilt (05)     module 05: low-pass filter
!   ==>  modadap (06)     module 06: adaptive parameter estimation
!   ==>  moddiff (07)     module 07: compute RMS difference between
!                            vector objects
!   ==>  modoerr (08)     module 08: observation error management
!   ==>  modoper (09)     module 09: arithmetic operation
!   ==>  modgroa (10)     module 10: global EOF decomposition
!   ==>  modlroa (11)     module 11: local EOF decomposition
!   ==>  modbroa (12)     module 12: bubble EOF decomposition
!   ==>  modzone (13)     module 13: local data section management
!   ==>  modgeof (14)     module 14: global reduced order analysis
!   ==>  modleof (15)     module 15: local reduced order analysis
!   ==>  modbeof (16)     module 16: bubble reduced order analysis
!   ==>  modgreg (17)     module 17: global linear regression
!   ==>  modlreg (18)     module 18: local linear regression
!   ==>  modbreg (19)     module 19: bubble linear regression
!   ==>  modvari (20)     module 20: compute variances
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mpitime
      use argsesam
      use hiomsk
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*) :: textcommand
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=13) :: cmd
      BIGREAL, parameter :: xprec=1.0_kr
      BIGREAL8, parameter :: xprec8=1.0_kr8
      BIGREAL4, parameter :: xprec4=1.0_kr4
      TYPE (type_swiarg), dimension(1:nbargmax) :: execaction
      INTEGER :: jloop,jindbeg,jindend
#if defined GETARG
      EXTERNAL getarg
#endif
!----------------------------------------------------------------------
!
      nprint = 0
!
! - Initialisation of the parallel code
#if defined MPI
      CALL mpi_init(mpi_code)
      CALL mpi_comm_size(mpi_comm_world,jpproc,mpi_code)
      CALL mpi_comm_rank(mpi_comm_world,jproc,mpi_code)
      call MPI_TIMER(0)
#endif
!
! -0.- Get user command from input
! --------------------------------
!
      execaction(:) = type_swiarg(' ',' ')
      IF (VERIFY(textcommand,' ').NE.0) THEN
!  get user command from routine argument
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
!  get user command from commandline arguments
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
!
! -1.- Open SESAM listing file and SESAM configuration file
! ---------------------------------------------------------
!
      DO jloop = 1, nbargmax
         IF (execaction(jloop)%swi(1:MAX(1,lenv(execaction(jloop)%swi))) &
     &        .EQ.'-outinfo') fnamout=execaction(jloop)%arg
      ENDDO
      CALL openfile(numout,fnamout,kstatus=clunk)
!
      DO jloop = 1, nbargmax
         IF (execaction(jloop)%swi(1:MAX(1,lenv(execaction(jloop)%swi))) &
     &        .EQ.'-config')    fnamlist=execaction(jloop)%arg
      ENDDO
      CALL openfile(numnam,fnamlist,kstatus=clunk)
!
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
!
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
!
! -2.- Initialisation of SESAM parameters :
! -----------------------------------------
!
! -2.1- Define SESAM constants
!
      CALL defcst
!
! -2.2- Interpret user command to know which module to invoke
!
      CALL readarg(execaction)
!
! -2.3- Read SESAM masks
!
      CALL readmsk
!
! -2.4- Set up SESAM internal configuration
!
      CALL evalconfig
!
! -3.- Run appropriate SESAM module
! ---------------------------------
!
      SELECT CASE (nmode)
      CASE (1)
! -3.1- INTF mode : interface an object from a file format to another
!
         CALL modintf
!
      CASE (2)
! -3.2- ---- mode : compute correlation coefficients
!
         CALL modcorr
!
      CASE (3)
! -3.3- ---- mode : truncated Gaussian operations
!
         CALL modtgop
!
      CASE (4)
! -3.4- OBSV mode : observation management
!
         CALL modobsv
!
      CASE (5)
! -3.5- FILT mode : low-pass filter
!      
         CALL modfilt
!
      CASE (6)
! -3.6- ADAP mode : adaptive parameter estimation
!
         CALL modadap
!
      CASE (7)
! -3.7- DIFF mode : compute RMS difference between vector objects
!
         CALL moddiff
!
      CASE (8)
! -3.8- OERR mode : observation error management
!
         CALL modoerr
!
      CASE (9)
! -3.9- OPER mode : arithmetic operation
!
         CALL modoper
!
      CASE (10)
! -3.10- GEOF mode : global EOF decomposition
!
         CALL modglbeof
!
      CASE (11)
! -3.11- LEOF mode : local EOF decomposition
!
         CALL modglbeof
!
      CASE (12)
! -3.12- BEOF mode : bubble EOF decomposition
!
         CALL modglbeof
!
      CASE (13)
! -3.13- ZONE mode : local data section management
!
         CALL modzone
!
      CASE (14)
! -3.14- GROA mode : global reduced order analysis
!
         CALL modglbroa
!
      CASE (15)
! -3.15- LROA mode : local reduced order analysis
!
         CALL modglbroa
!
      CASE (16)
! -3.16- BROA mode : bubble reduced order analysis
!
         CALL modglbroa
!
      CASE (17)
! -3.17- GREG mode : global linear regression
!
         CALL modglbreg
!
      CASE (18)
! -3.18- LREG mode : local linear regression
!
         CALL modglbreg
!
      CASE (19)
! -3.19- BREG mode : bubble linear regression
!
         CALL modglbreg
!
      CASE (20)
! -3.20- VARI mode : variance computation
!
         CALL modvari
!
      CASE (21)
! -3.21- ANAM mode : anamorphosis
!
         CALL modanam
!
      CASE (22)
! -3.22- SCOR mode : probabilistic scores (CRPS,...)
!
         CALL modscor
!
      CASE (23)
! -3.23- SPCT mode : spectral transformation (spherical harmonics)
!
         CALL modspct
!
      CASE (24)
! -3.24- MCMC mode : MCMC sampler
!
         CALL modmcmc
!
      CASE (25)
! -3.24- RANK mode : computation of ranks of data within an ensemble
!
         CALL modrank
!
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -4.- Create empty 'SESAMOK' file 
! --------------------------------
! (to prove the job has terminated successfully)
!
      IF ((nprint.GE.0).AND.(jproc.EQ.0)) THEN
         cmd='touch SESAMOK'
         CALL shellorder(cmd)
      ENDIF
!
! -5.- This is the End
! --------------------
!
      IF ((nprint.GE.0).AND.(jproc.EQ.0)) THEN
         WRITE(numout,*)
         WRITE(numout,*) '@'
         WRITE(numout,*) '@ End of ',version(:lenv(version))
         WRITE(numout,*) repeat('@',59)
      ENDIF
!
#if defined MPI
      call MPI_TIMER(1)
      CALL mpi_finalize(mpi_code)
#endif
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'sesam','sesam')
!
      END
