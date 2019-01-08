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
CC---------------------------------------------------------------------
CC
CC                        PARAM.CONTROL.H
CC
CC---------------------------------------------------------------------
CC
CC  PURPOSE :
CC  ---------
CC  Include file with SESAM main parameters
CC  (size of SESAM configuration arrays, kind of real
CC  variables, i/o options, SESAM version,...).
CC
CC  MODIFICATIONS :
CC  ---------------
CC  original     : 97-12 (C.E. Testut)
CC  modification : 99-10 (C.E. Testut)
CC  modification : 03-02 (J.M. Brankart)
CC---------------------------------------------------------------------
CC A/ User configuration size parameters
CC---------------------------------------------------------------------
C Maximum number of variables in state vector object
      INTEGER, parameter :: nbvar=32
C Maximum number of nested grids
      INTEGER, parameter :: nbnest=3
C Maximum number of observation data sets for each observed variable
      INTEGER, parameter :: jpndbs=9
C Maximum number of observation data sets
      INTEGER, parameter :: nbobs=nbvar*jpndbs
CC
CC---------------------------------------------------------------------
CC B/ SESAM configuration size parameters
CC---------------------------------------------------------------------
CC
C Length of most character variables
      INTEGER, parameter :: word80=80,bgword=120,hgword=2048
C Length of specific character variables
C (switches, module names, help keywords, file extensions, variable names)
      INTEGER, parameter :: swilg=12,modlg=5,helplg=6,extlg=6,varlg=5
C
C Number of possible commandline switches
      INTEGER, parameter :: nbarg=89
C Number of modules, of help keywords
      INTEGER, parameter :: nbmod=23,nbhelp=3
C Number of optional switches
      INTEGER, parameter :: nbswiopt=14
C Maximum number of actions per module
      INTEGER, parameter :: nbaction=30
C Maximum number of required switches per action
      INTEGER, parameter :: nborder=12
C Maximum number of argument in commandline
      INTEGER, parameter :: nbargmax = ( nbswiopt + 5 + 1 + nborder)
C Maximum number of possible file formats (bas, dbs)
      INTEGER, parameter :: nbext=2,nbextbas=1,nbextdbs=4
C Maximum number of possible file formats (dta, obs, var, zon)
      INTEGER, parameter :: nbextdta=5,nbextobs=3,nbextvar=8,nbextzon=2
C Maximum number of variables in user defined '.rst' files
      INTEGER, parameter :: nbinf=10
C Number of possible space to work in: (x, y, o)
      INTEGER, parameter :: nbxyo=3
C Maximum number of z objects that is possible to load simultaneously
      INTEGER, parameter :: nbztyp=4
C     INTEGER, parameter :: jpndbfr=1
CC
CC---------------------------------------------------------------------
CC C/ SESAM file units and file names: default values
CC---------------------------------------------------------------------
CC
      INTEGER, parameter :: numout=1
      CHARACTER(len=bgword) :: fnamout='sesam.output'
      CHARACTER(len=bgword) :: fnamlist='sesamlist'
      INTEGER, parameter :: numnam=2,numfil=3,numfil1=10,numfil2=20
      CHARACTER(len=1) :: etoile='#'
CC
CC---------------------------------------------------------------------
CC D/ Parameters defining kind of real variables
CC---------------------------------------------------------------------
CC
#if defined KR4EQ4
      INTEGER, parameter :: kr4=4
#elif defined KR4EQ8
      INTEGER, parameter :: kr4=8
#else
      INTEGER, parameter :: kr4=4
#endif
#if defined KR8EQ8
      INTEGER, parameter :: kr8=8
#else
      INTEGER, parameter :: kr8=8
#endif
#if defined KREQ4
      INTEGER, parameter :: kr=4
#elif defined KREQ8
      INTEGER, parameter :: kr=8
#elif defined KREQ16
      INTEGER, parameter :: kr=16
#else
      INTEGER, parameter :: kr=8
#endif
CC
CC---------------------------------------------------------------------
CC E/ Parameters for opening direct access files
CC---------------------------------------------------------------------
CC
#if defined _NEC
      INTEGER, parameter :: ibloc=1
#else
      INTEGER, parameter :: ibloc=512
#endif
      INTEGER, parameter :: jpbyt8=8,jpbyt4=4,jpmsk=2
CC
CC---------------------------------------------------------------------
CC F/ SESAM version
CC---------------------------------------------------------------------
CC
      CHARACTER(len=bgword) :: version='SESAM 4.3', version_year='2012'
CC---------------------------------------------------------------------
