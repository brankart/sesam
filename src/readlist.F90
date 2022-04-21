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
! ---                  READLIST.F90                               ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 07-11  (J.M. Brankart)                     ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE readlist
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readlist 
!---------------------------------------------------------------------
!
! Purpose : Read SESAM configuration file
! -------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use algospct, only : loc_time_scale, loc_radius_in_deg
      use utilconstraint, only : dyn_constraint, dyn_constraint_std
      use ensdam_mcmc_update
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: indvar, jndbs, integ, eof
      CHARACTER(len=bgword) :: line, key, attr, word
      BIGREAL :: bigreal
      BIGREAL4 :: bigreal4
      LOGICAL :: boolean
      LOGICAL, dimension(1:nbvar) :: lvarinam, lvaronam, lvardmsk
      LOGICAL, dimension(1:nbvar) :: lvarfmsk, lvarvmsk, lvaremsk
      LOGICAL, dimension(1:nbvar) :: lvarmsea, lvarifil, lvarofil
      LOGICAL, dimension(1:nbvar) :: lvarxdim, lvarydim, lvarzdim
      LOGICAL, dimension(1:nbvar) :: lvartdim, lvarxnam, lvarynam
      LOGICAL, dimension(1:nbvar) :: lvarznam, lvartnam
      LOGICAL, dimension(1:nbvar) :: lvarfgrd, lvaregrd, lvarngrd
      LOGICAL, dimension(1:nbvar) :: ldta_nam, ldta_dim
      LOGICAL, dimension(1:nbvar) :: ldta_moy, ldta_ect
      LOGICAL, dimension(1:nbvar) :: ldtainam, ldtaonam, ldtadmsk
      LOGICAL, dimension(1:nbvar) :: ldtafmsk, ldtavmsk, ldtaemsk
      LOGICAL, dimension(1:nbvar) :: ldtamsea, ldtaifil, ldtaofil
      LOGICAL, dimension(1:nbvar) :: ldtafgrd, ldtaegrd, ldtangrd
      LOGICAL, dimension(1:nbvar) :: ldtaflev, ldtaelev, ldtaxdim
      LOGICAL, dimension(1:nbvar) :: ldtaydim, ldtazdim, ldtatdim
      LOGICAL, dimension(1:nbvar) :: ldtaxnam, ldtaynam, ldtaznam
      LOGICAL, dimension(1:nbvar) :: ldtatnam
      LOGICAL, dimension(1:nbvar,1:jpndbs) :: lobs_dim, lobs_moy, lobs_ect
      LOGICAL, dimension(1:nbvar,1:jpndbs) :: lobsinam, lobsonam, lobs_rms
      LOGICAL, dimension(1:nbvar,1:jpndbs) :: lobsifil
!----------------------------------------------------------------------
! Initialize logical arrays to TRUE (=parameters not read)
      lvarinam(1:nbvar) = .TRUE.
      lvarifil(1:nbvar) = .TRUE.
      lvaronam(1:nbvar) = .TRUE.
      lvarofil(1:nbvar) = .TRUE.
      lvarxdim(1:nbvar) = .TRUE.
      lvarydim(1:nbvar) = .TRUE.
      lvarzdim(1:nbvar) = .TRUE.
      lvartdim(1:nbvar) = .TRUE.
      lvardmsk(1:nbvar) = .TRUE.
      lvarfmsk(1:nbvar) = .TRUE.
      lvaremsk(1:nbvar) = .TRUE.
      lvarvmsk(1:nbvar) = .TRUE.
      lvarmsea(1:nbvar) = .TRUE.
      lvarfgrd(1:nbvar) = .TRUE.
      lvaregrd(1:nbvar) = .TRUE.
      lvarngrd(1:nbvar) = .TRUE.
      lvarxnam(1:nbvar) = .TRUE.
      lvarynam(1:nbvar) = .TRUE.
      lvarznam(1:nbvar) = .TRUE.
      lvartnam(1:nbvar) = .TRUE.
!
      ldta_nam(1:nbvar) = .TRUE.
      ldta_dim(1:nbvar) = .TRUE.
      ldta_moy(1:nbvar) = .TRUE.
      ldta_ect(1:nbvar) = .TRUE.
      ldtainam(1:nbvar) = .TRUE.
      ldtaifil(1:nbvar) = .TRUE.
      ldtaonam(1:nbvar) = .TRUE.
      ldtaofil(1:nbvar) = .TRUE.
      ldtaxdim(1:nbvar) = .TRUE.
      ldtaydim(1:nbvar) = .TRUE.
      ldtazdim(1:nbvar) = .TRUE.
      ldtatdim(1:nbvar) = .TRUE.
      ldtaxdim(1:nbvar) = .TRUE.
      ldtaydim(1:nbvar) = .TRUE.
      ldtazdim(1:nbvar) = .TRUE.
      ldtatdim(1:nbvar) = .TRUE.
      ldtadmsk(1:nbvar) = .TRUE.
      ldtafmsk(1:nbvar) = .TRUE.
      ldtaemsk(1:nbvar) = .TRUE.
      ldtavmsk(1:nbvar) = .TRUE.
      ldtamsea(1:nbvar) = .TRUE.
      ldtafgrd(1:nbvar) = .TRUE.
      ldtaegrd(1:nbvar) = .TRUE.
      ldtangrd(1:nbvar) = .TRUE.
      ldtaxnam(1:nbvar) = .TRUE.
      ldtaynam(1:nbvar) = .TRUE.
      ldtaznam(1:nbvar) = .TRUE.
      ldtatnam(1:nbvar) = .TRUE.
      ldtaflev(1:nbvar) = .TRUE.
      ldtaelev(1:nbvar) = .TRUE.
!
      lobs_dim(1:nbvar,1:jpndbs) = .TRUE.
      lobs_moy(1:nbvar,1:jpndbs) = .TRUE.
      lobs_ect(1:nbvar,1:jpndbs) = .TRUE.
      lobsinam(1:nbvar,1:jpndbs) = .TRUE.
      lobsonam(1:nbvar,1:jpndbs) = .TRUE.
      lobs_rms(1:nbvar,1:jpndbs) = .TRUE.
      lobsifil(1:nbvar,1:jpndbs) = .TRUE.
!
      IF ((nprint.GE.0).AND.(jproc.EQ.0)) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/defctl/readlist :'
         WRITE(numout,*) '         read SESAM configuration file'
      ENDIF
!
! -1.- Loop on uncommented line of configuration files
! ----------------------------------------------------
!
      loopconfig : DO
!
        line='#'
        DO WHILE (line(1:1)=='#'.or.LEN_TRIM(line)==0)
          READ(numnam,*,IOSTAT=eof) line
          if (eof < 0 ) EXIT loopconfig
        ENDDO
!
        CALL readconfigline(line,key,indvar,jndbs,attr)
        IF (indvar.GT.nbvar) GOTO 103
!
        SELECT CASE(key)
!
! Get global SESAM parameters
        CASE('NPRINT')
          READ(attr,*,ERR=102) nprint
          IF (jproc.NE.0) nprint = 0
        CASE('TRADITIONAL')
          READ(attr,*,ERR=102) traditional
        CASE('LMOYECT')
          READ(attr,*,ERR=102) lmoyect
        CASE('FACTOUBLI')
          READ(attr,*,ERR=102) factoubli
        CASE('FORGEXP')
          READ(attr,*,ERR=102) forgexp
        CASE('BAS_DIGITS')
          READ(attr,*,ERR=102) nbdigits
        CASE('DISABLE')
          READ(attr,*,ERR=102) argdisable
!
! Observation error type
        CASE('OBSERROR_CDF')
          READ(attr,*,ERR=102) obserror_type_sesam
        CASE('OECORRELTYP')
          READ(attr,*,ERR=102) oecorreltyp
          IF ((oecorreltyp.LE.0).OR.(oecorreltyp.GT.2)) GOTO 106
!
! Parameterization of MCMC sampler
        CASE('MCMC_NSCALE')
          READ(attr,*,ERR=102) jpscl
          IF ((jpscl.LT.1).OR.(jpscl.GT.nbscl)) GOTO 106
        CASE('MCMC_SCALE_MULTIPLICITY')
          IF (indvar.LE.0) GOTO 108
          IF (indvar.GT.nbscl) GOTO 109
          READ(attr,*,ERR=102) scl_mult(indvar)
        CASE('MCMC_CONTROL_PRINT')
          READ(attr,*,ERR=102) mcmc_control_print
        CASE('MCMC_CONVERGENCE_CHECK')
          READ(attr,*,ERR=102) mcmc_convergence_check
        CASE('MCMC_CONVERGENCE_STOP')
          READ(attr,*,ERR=102) mcmc_convergence_stop
        CASE('MCMC_ZERO_START')
          READ(attr,*,ERR=102) mcmc_zero_start
!
! Parameterization of dynamical constraint
        CASE('DYN_CONSTRAINT')
          READ(attr,*,ERR=102) dyn_constraint
        CASE('DYN_CONSTRAINT_STD')
          READ(attr,*,ERR=102) dyn_constraint_std
!
! Parameterization of localization radius fro random sampling
        CASE('LOC_TIME_SCALE')
          READ(attr,*,ERR=102) loc_time_scale
        CASE('LOC_RADIUS_IN_DEG')
          READ(attr,*,ERR=102) loc_radius_in_deg
!
! Parameterization of regression in SPCT module
        CASE('REGR_TYPE')
          READ(attr,*,ERR=102) regr_type_sesam
        CASE('REGR_MAXITER')
          READ(attr,*,ERR=102) regr_maxiter_sesam
          IF (regr_maxiter_sesam.LE.0) GOTO 106
        CASE('REGR_MAXBLOC')
          READ(attr,*,ERR=102) regr_maxbloc_sesam
          IF (regr_maxbloc_sesam.LE.0) GOTO 106
        CASE('REGR_OVERLAP')
          READ(attr,*,ERR=102) regr_overlap_sesam
          IF (regr_overlap_sesam.LE.0) GOTO 106
        CASE('REGR_EPSILON')
          READ(attr,*,ERR=102) regr_epsilon_sesam
          IF (regr_epsilon_sesam.LE.0.0) GOTO 106
        CASE('REGR_RHO')
          READ(attr,*,ERR=102) regr_rho_sesam
          IF ((regr_rho_sesam.LT.0.0).OR.(regr_rho_sesam.GT.1.0)) GOTO 106
        CASE('SPECIAL_VALUE')
          READ(attr,*,ERR=102) special_value
!
! ==> A State vector configuration
!
! ==> A.1. Configuration of the variable
        CASE('VAREND')
          READ(attr,*,ERR=102) varend
          IF (varend.GT.nbvar) GOTO 103
        CASE('VAR_NAM')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) var_nam(indvar)
        CASE('VAR_DIM')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) var_dim(indvar)
        CASE('VARNEST')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) varnest(indvar)
        CASE('VAR_MOY')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) var_moy(indvar)
        CASE('VAR_ECT')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) var_ect(indvar)
! ==> A.2. Configuration of the input/ouput files
        CASE('VARINAM')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) varinam(indvar)
          lvarinam(indvar) = .FALSE.
        CASE('VARIFIL')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) varifil(indvar)
          lvarifil(indvar) = .FALSE.
        CASE('VARIPOS')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) varipos(indvar)
        CASE('VARONAM')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) varonam(indvar)
          lvaronam(indvar) = .FALSE.
        CASE('VAROFIL')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) varofil(indvar)
          lvarofil(indvar) = .FALSE.
        CASE('VAROPOS')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) varopos(indvar)
        CASE('VARXDIM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            varxdim(indvar) = word
            lvarxdim(indvar) = .FALSE.
          ELSE
            WHERE (lvarxdim(1:nbvar)) varxdim(:) = word
          ENDIF
        CASE('VARYDIM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            varydim(indvar) = word
            lvarydim(indvar) = .FALSE.
          ELSE
            WHERE (lvarydim(1:nbvar)) varydim(:) = word
          ENDIF
        CASE('VARZDIM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            varzdim(indvar) = word
            lvarzdim(indvar) = .FALSE.
          ELSE
            WHERE (lvarzdim(1:nbvar)) varzdim(:) = word
          ENDIF
        CASE('VARTDIM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            vartdim(indvar) = word
            lvartdim(indvar) = .FALSE.
          ELSE
            WHERE (lvartdim(1:nbvar)) vartdim(:) = word
          ENDIF
! ==> A.3. Configuration of the mask file
        CASE('VARFMSK')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            varfmsk(indvar) = word
            lvarfmsk(indvar) = .FALSE.
          ELSE
            WHERE (lvarfmsk(1:nbvar)) varfmsk(:) = word
          ENDIF
        CASE('VAREMSK')
          READ(attr,*,ERR=102) integ
          IF (indvar.NE.0) THEN
            varemsk(indvar) = integ
            lvaremsk(indvar) = .FALSE.
          ELSE
            WHERE (lvaremsk(1:nbvar)) varemsk(:) = integ
          ENDIF
        CASE('VARDMSK')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) vardmsk(indvar)
          lvardmsk(indvar) = .FALSE.
        CASE('VARPMSK')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) varpmsk(indvar)
        CASE('VARVMSK')
          READ(attr,*,ERR=102) bigreal4
          IF (indvar.NE.0) THEN
            varvmsk(indvar) = bigreal4
            lvarvmsk(indvar) = .FALSE.
          ELSE
            WHERE (lvarvmsk(1:nbvar)) varvmsk(:) = bigreal4
          ENDIF
        CASE('VARMSEA')
          READ(attr,*,ERR=102) boolean
          IF (indvar.NE.0) THEN
            varmsea(indvar) = boolean
            lvarmsea(indvar) = .FALSE.
          ELSE
            WHERE (lvarmsea(1:nbvar)) varmsea(:) = boolean
          ENDIF
! ==> A.4. Configuration of grid files
        CASE('VARFGRD')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            varfgrd(indvar) = word
            lvarfgrd(indvar) = .FALSE.
          ELSE
            WHERE (lvarfgrd(1:nbvar)) varfgrd(:) = word
          ENDIF
        CASE('VAREGRD')
          READ(attr,*,ERR=102) integ
          IF (indvar.NE.0) THEN
            varegrd(indvar) = integ
            lvaregrd(indvar) = .FALSE.
          ELSE
            WHERE (lvaregrd(1:nbvar)) varegrd(:) = integ
          ENDIF
        CASE('VARNGRD')
          READ(attr,*,ERR=102) integ
          IF (indvar.NE.0) THEN
            varngrd(indvar) = integ
            lvarngrd(indvar) = .FALSE.
          ELSE
            WHERE (lvarngrd(1:nbvar)) varngrd(:) = integ
          ENDIF
        CASE('VARXNAM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            varxnam(indvar) = word
            lvarxnam(indvar) = .FALSE.
          ELSE
            WHERE (lvarxnam(1:nbvar)) varxnam(:) = word
          ENDIF
        CASE('VARYNAM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            varynam(indvar) = word
            lvarynam(indvar) = .FALSE.
          ELSE
            WHERE (lvarynam(1:nbvar)) varynam(:) = word
          ENDIF
        CASE('VARZNAM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            varznam(indvar) = word
            lvarznam(indvar) = .FALSE.
          ELSE
            WHERE (lvarznam(1:nbvar)) varznam(:) = word
          ENDIF
        CASE('VARTNAM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            vartnam(indvar) = word
            lvartnam(indvar) = .FALSE.
          ELSE
            WHERE (lvartnam(1:nbvar)) vartnam(:) = word
          ENDIF
!
! ==> B Data section vector configuration
!
! ==> B.0. Does the data section include this variable?
        CASE('DTA_ACT')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) dta_act(indvar)
! ==> B.1. Configuration of the data section variable
        CASE('DTA_NAM')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) dta_nam(indvar)
          ldta_nam(indvar) = .FALSE.
        CASE('DTA_DIM')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) dta_dim(indvar)
          ldta_dim(indvar) = .FALSE.
        CASE('DTA_MOY')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) var_moy(indvar)
          ldta_moy(indvar) = .FALSE.
        CASE('DTA_ECT')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) var_ect(indvar)
          ldta_ect(indvar) = .FALSE.
! ==> B.2. Configuration of the input/ouput files
        CASE('DTAINAM')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) dtainam(indvar)
          ldtainam(indvar) = .FALSE.
        CASE('DTAIFIL')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) dtaifil(indvar)
          ldtaifil(indvar) = .FALSE.
        CASE('DTAIPOS')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) dtaipos(indvar)
        CASE('DTAONAM')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) dtaonam(indvar)
          ldtaonam(indvar) = .FALSE.
        CASE('DTAOFIL')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) dtaofil(indvar)
          ldtaofil(indvar) = .FALSE.
        CASE('DTAOPOS')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) dtaopos(indvar)
        CASE('DTAXDIM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            dtaxdim(indvar) = word
            ldtaxdim(indvar) = .FALSE.
          ELSE
            WHERE (ldtaxdim(1:nbvar)) dtaxdim(:) = word
            ldtaxdim(:) = .FALSE.
          ENDIF
        CASE('DTAYDIM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            dtaydim(indvar) = word
            ldtaydim(indvar) = .FALSE.
          ELSE
            WHERE (ldtaydim(1:nbvar)) dtaydim(:) = word
            ldtaydim(:) = .FALSE.
          ENDIF
        CASE('DTAZDIM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            dtazdim(indvar) = word
            ldtazdim(indvar) = .FALSE.
          ELSE
            WHERE (ldtazdim(1:nbvar)) dtazdim(:) = word
            ldtazdim(:) = .FALSE.
          ENDIF
        CASE('DTATDIM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            dtatdim(indvar) = word
            ldtatdim(indvar) = .FALSE.
          ELSE
            WHERE (ldtatdim(1:nbvar)) dtatdim(:) = word
            ldtatdim(:) = .FALSE.
          ENDIF
! ==> B.3. Configuration of the mask file
        CASE('DTAFMSK')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            dtafmsk(indvar) = word
            ldtafmsk(indvar) = .FALSE.
          ELSE
            WHERE (ldtafmsk(1:nbvar)) dtafmsk(:) = word
          ENDIF
        CASE('DTAEMSK')
          READ(attr,*,ERR=102) integ
          IF (indvar.NE.0) THEN
            dtaemsk(indvar) = integ
            ldtaemsk(indvar) = .FALSE.
          ELSE
            WHERE (ldtaemsk(1:nbvar)) dtaemsk(:) = integ
          ENDIF
        CASE('DTADMSK')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) dtadmsk(indvar)
          ldtadmsk(indvar) = .FALSE.
        CASE('DTAPMSK')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) dtapmsk(indvar)
        CASE('DTAVMSK')
          READ(attr,*,ERR=102) bigreal4
          IF (indvar.NE.0) THEN
            dtavmsk(indvar) = bigreal4
            ldtavmsk(indvar) = .FALSE.
          ELSE
            WHERE (ldtavmsk(1:nbvar)) dtavmsk(:) = bigreal4
          ENDIF
        CASE('DTAMSEA')
          READ(attr,*,ERR=102) boolean
          IF (indvar.NE.0) THEN
            dtamsea(indvar) = boolean
            ldtamsea(indvar) = .FALSE.
          ELSE
            WHERE (ldtamsea(1:nbvar)) dtamsea(:) = boolean
          ENDIF
! ==> B.4. Observation error standard deviation
        CASE('DTA_RMS')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) dta_rms(indvar)
! ==> B.5. Configuration of grid and level files
        CASE('DTAFGRD')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            dtafgrd(indvar) = word
            ldtafgrd(indvar) = .FALSE.
          ELSE
            WHERE (ldtafgrd(1:nbvar)) dtafgrd(:) = word
          ENDIF
        CASE('DTAEGRD')
          READ(attr,*,ERR=102) integ
          IF (indvar.NE.0) THEN
            dtaegrd(indvar) = integ
            ldtaegrd(indvar) = .FALSE.
          ELSE
            WHERE (ldtaegrd(1:nbvar)) dtaegrd(:) = integ
          ENDIF
        CASE('DTANGRD')
          READ(attr,*,ERR=102) integ
          IF (indvar.NE.0) THEN
            dtangrd(indvar) = integ
            ldtangrd(indvar) = .FALSE.
          ELSE
            WHERE (ldtangrd(1:nbvar)) dtangrd(:) = integ
          ENDIF
        CASE('DTAXNAM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            dtaxnam(indvar) = word
            ldtaxnam(indvar) = .FALSE.
          ELSE
            WHERE (ldtaxnam(1:nbvar)) dtaxnam(:) = word
          ENDIF
        CASE('DTAYNAM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            dtaynam(indvar) = word
            ldtaynam(indvar) = .FALSE.
          ELSE
            WHERE (ldtaynam(1:nbvar)) dtaynam(:) = word
          ENDIF
        CASE('DTAZNAM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            dtaznam(indvar) = word
            ldtaznam(indvar) = .FALSE.
          ELSE
            WHERE (ldtaznam(1:nbvar)) dtaznam(:) = word
          ENDIF
        CASE('DTATNAM')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            dtatnam(indvar) = word
            ldtatnam(indvar) = .FALSE.
          ELSE
            WHERE (ldtatnam(1:nbvar)) dtatnam(:) = word
          ENDIF
        CASE('DTAFLEV')
          READ(attr,*,ERR=102) word
          IF (indvar.NE.0) THEN
            dtaflev(indvar) = word
            ldtaflev(indvar) = .FALSE.
          ELSE
            WHERE (ldtaflev(1:nbvar)) dtaflev(:) = word
          ENDIF
        CASE('DTAELEV')
          READ(attr,*,ERR=102) integ
          IF (indvar.NE.0) THEN
            dtaelev(indvar) = integ
            ldtaelev(indvar) = .FALSE.
          ELSE
            WHERE (ldtaelev(1:nbvar)) dtaelev(:) = integ
          ENDIF
!
! ==> C Observation vector configuration
!
        CASE('OBSNDBS')
          IF (indvar.EQ.0) GOTO 104
          READ(attr,*,ERR=102) obsndbs(indvar)
! ==> C.1. Configuration of observation variables
        CASE('OBS_NAM')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_nam(indvar,jndbs)
        CASE('OBS_DIM')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_dim(indvar,jndbs)
          lobs_dim(indvar,jndbs) = .FALSE.
        CASE('OBS_MOY')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_moy(indvar,jndbs)
          lobs_moy(indvar,jndbs) = .FALSE.
        CASE('OBS_ECT')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_ect(indvar,jndbs)
          lobs_ect(indvar,jndbs) = .FALSE.
! ==> C.2. Configuration of the input/ouput files
        CASE('OBSINAM')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsinam(indvar,jndbs)
          lobsinam(indvar,jndbs) = .FALSE.
        CASE('OBSONAM')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsonam(indvar,jndbs)
          lobsonam(indvar,jndbs) = .FALSE.
! ==> C.3. Observation error standard deviation
        CASE('OBS_RMS')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_rms(indvar,jndbs)
          lobs_rms(indvar,jndbs) = .FALSE.
! ==> C.4. Configuration of database files (.ncdbs)
        CASE('OBSIFIL')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsifil(indvar,jndbs)
          lobsifil(indvar,jndbs) = .FALSE.
        CASE('OBSXDIM')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsxdim(indvar,jndbs)
        CASE('OBSXNAM')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsxnam(indvar,jndbs)
        CASE('OBSYDIM')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsydim(indvar,jndbs)
        CASE('OBSYNAM')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsynam(indvar,jndbs)
        CASE('OBSZDIM')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obszdim(indvar,jndbs)
        CASE('OBSZNAM')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsznam(indvar,jndbs)
        CASE('OBSTDIM')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obstdim(indvar,jndbs)
        CASE('OBSTNAM')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obstnam(indvar,jndbs)
! ==> C.5. Configuration of the extraction from the database
        CASE('OBSEXCL')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsexcl(indvar,jndbs)
        CASE('OBS_MIN')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_min(indvar,jndbs)
        CASE('OBS_MAX')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_max(indvar,jndbs)
!
        CASE('OBSIMIN')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsimin(indvar,jndbs)
        CASE('OBSIMAX')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsimax(indvar,jndbs)
        CASE('OBSJMIN')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsjmin(indvar,jndbs)
        CASE('OBSJMAX')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obsjmax(indvar,jndbs)
        CASE('OBSKMIN')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obskmin(indvar,jndbs)
        CASE('OBSKMAX')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obskmax(indvar,jndbs)
        CASE('OBSTMIN')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obstmin(indvar,jndbs)
        CASE('OBSTMAX')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obstmax(indvar,jndbs)
!
        CASE('OBS_LON_DEF')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_lon_def(indvar,jndbs)
        CASE('OBS_LAT_DEF')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_lat_def(indvar,jndbs)
        CASE('OBS_DEP_DEF')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_dep_def(indvar,jndbs)
        CASE('OBS_TIM_DEF')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_tim_def(indvar,jndbs)
!
        CASE('OBS_LON_MIN')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_lon_min(indvar,jndbs)
        CASE('OBS_LON_MAX')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_lon_max(indvar,jndbs)
        CASE('OBS_LAT_MIN')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_lat_min(indvar,jndbs)
        CASE('OBS_LAT_MAX')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_lat_max(indvar,jndbs)
        CASE('OBS_DEP_MIN')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_dep_min(indvar,jndbs)
        CASE('OBS_DEP_MAX')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_dep_max(indvar,jndbs)
        CASE('OBS_TIM_MIN')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_tim_min(indvar,jndbs)
        CASE('OBS_TIM_MAX')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_tim_max(indvar,jndbs)
!
        CASE('OBS_SIZ_MAX')
          IF (indvar.EQ.0) GOTO 104
          IF (jndbs.EQ.0) GOTO 105
          READ(attr,*,ERR=102) obs_siz_max(indvar,jndbs)
!
! ==> D Nesting configuration
        CASE('NESTEND')
          READ(attr,*,ERR=102) nestend
          IF ((nestend.LT.0).OR.(nestend.GT.nbnest)) GOTO 107
        CASE('NESTXORI')
          IF ((indvar.LT.0).OR.(indvar.GT.nbnest)) GOTO 107
          READ(attr,*,ERR=102) nestxori(indvar)
        CASE('NESTYORI')
          IF ((indvar.LT.0).OR.(indvar.GT.nbnest)) GOTO 107
          READ(attr,*,ERR=102) nestyori(indvar)
        CASE('NESTXRES')
          IF ((indvar.LT.0).OR.(indvar.GT.nbnest)) GOTO 107
          READ(attr,*,ERR=102) nestxres(indvar)
        CASE('NESTYRES')
          IF ((indvar.LT.0).OR.(indvar.GT.nbnest)) GOTO 107
          READ(attr,*,ERR=102) nestyres(indvar)
!
        CASE DEFAULT
          GOTO 101
        END SELECT
!
      ENDDO loopconfig
!
! -2.- Put the default values to undefined parameters
! ---------------------------------------------------
!
      WHERE (lvarinam(:)) varinam(:) = var_nam(:)
      WHERE (lvarifil(:)) varifil(:) = var_nam(:)
      WHERE (lvaronam(:)) varonam(:) = var_nam(:)
      WHERE (lvarofil(:)) varofil(:) = var_nam(:)
      WHERE (lvardmsk(:)) vardmsk(:) = var_dim(:)
!
      WHERE (ldta_nam(:)) dta_nam(:) = var_nam(:)
      WHERE (ldta_dim(:)) dta_dim(:) = var_dim(:)
      WHERE (ldta_moy(:)) dta_moy(:) = var_moy(:)
      WHERE (ldta_ect(:)) dta_moy(:) = var_ect(:)
!
      WHERE (ldtainam(:)) dtainam(:) = dta_nam(:)
      WHERE (ldtaifil(:)) dtaifil(:) = dta_nam(:)
      WHERE (ldtaonam(:)) dtaonam(:) = dta_nam(:)
      WHERE (ldtaofil(:)) dtaofil(:) = dta_nam(:)
      WHERE (ldtaxdim(:)) dtaxdim(:) = varxdim(:)
      WHERE (ldtaydim(:)) dtaydim(:) = varydim(:)
      WHERE (ldtazdim(:)) dtazdim(:) = varzdim(:)
      WHERE (ldtatdim(:)) dtatdim(:) = vartdim(:)
      WHERE (ldtadmsk(:)) dtadmsk(:) = dta_dim(:)
!
      WHERE (lobsinam(:,:)) obsinam(:,:) = obs_nam(:,:)
      WHERE (lobsonam(:,:)) obsonam(:,:) = obs_nam(:,:)
      WHERE (lobsifil(:,:)) obsifil(:,:) = obs_nam(:,:)
!
      DO jndbs=1,jpndbs
        WHERE (lobs_dim(:,jndbs)) obs_dim(:,jndbs) = dta_dim(:)
        WHERE (lobs_moy(:,jndbs)) obs_moy(:,jndbs) = dta_moy(:)
        WHERE (lobs_ect(:,jndbs)) obs_ect(:,jndbs) = dta_ect(:)
        WHERE (lobs_rms(:,jndbs)) obs_rms(:,jndbs) = dta_rms(:)
      ENDDO
!
! set order of variables in memory
      DO indvar=1,varend
         var_ord(indvar)=indvar
      ENDDO
!
! -3.- Check global SESAM parameters
! ----------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '    SESAM global parameters :'
         WRITE(numout,*) '     lmoyect       = ',lmoyect
         WRITE(numout,*) '     nprint        = ',nprint
         WRITE(numout,*) '     factoubli     = ',factoubli
      ENDIF
!
      IF (nprint.LT.0) THEN
         print *,'WARNING: Bad nprint value in namelist'
         print *,'         nprint set to 0'
         nprint  = 0
      ENDIF
      IF (nprint.GT.4) THEN
         print *,'WARNING: Bad nprint value in namelist'
         print *,'         nprint set to 3'
         nprint  = 4
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'readlist','readlist')
!
 101  WRITE (texterror,*) 'Bad keyword in configuration file: ',key
      CALL printerror2(0,101,3,'readlist','readlist',comment=texterror)
 102  WRITE (texterror,*) 'Bad attribute in configuration file: ',attr
      CALL printerror2(0,102,3,'readlist','readlist',comment=texterror)
 103  WRITE (texterror,*) 'Too large var index in config file'
      CALL printerror2(0,103,3,'readlist','readlist',comment=texterror)
 104  WRITE (texterror,*) 'Missing var index in config file key: ',key
      CALL printerror2(0,104,3,'readlist','readlist',comment=texterror)
 105  WRITE (texterror,*) 'Missing dbs index in config file key: ',key
      CALL printerror2(0,105,3,'readlist','readlist',comment=texterror)
 106  WRITE (texterror,*) 'Bad parameter value in sesamlist: ',key
      CALL printerror2(0,106,3,'readlist','readlist',comment=texterror)
 107  WRITE (texterror,*) 'Bad nesting level in config file key: ',key
      CALL printerror2(0,107,3,'readlist','readlist',comment=texterror)
 108  WRITE (texterror,*) 'Missing scale index in config file key: ',key
      CALL printerror2(0,108,3,'readlist','readlist',comment=texterror)
 109  WRITE (texterror,*) 'Bad scale index in config file key: ',key
      CALL printerror2(0,109,3,'readlist','readlist',comment=texterror)
!
      END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readconfigline(line,key,indvar,jndbs,attr)
!
! --- Module declaration
      use mod_main , only : printerror2
      use mod_cfgxyo , only : texterror
      IMPLICIT NONE
! --- Variable declaration
      CHARACTER(len=*), intent(inout) :: line
      CHARACTER(len=*), intent(out) :: key, attr
      INTEGER, intent(out) :: indvar, jndbs
!
      INTEGER :: idx, len_line
! -----------------------------------------------------------------
      indvar=0 ; jndbs=0
!
! Eliminate comments after the # symbol
      key=line
      len_line=len_trim(line)
      DO idx=1,len_line
        IF (line(idx:idx)=='#') THEN
          key=line(1:idx-1)
          EXIT
        ENDIF
      ENDDO
!
! Separate key (before =) and attribute (after =)
      line=key
      len_line=len_trim(line)
      DO idx=1,len_line
        IF (line(idx:idx)=='=') THEN
          key=line(1:idx-1)
          attr=line(idx+1:)
          EXIT
        ENDIF
        IF (idx==len_line) GOTO 101
      ENDDO
!
! Extract observation database index from key (if any)
      line=key
      len_line=len_trim(line)
      DO idx=1,len_line
        IF (line(idx:idx)==':') THEN
          key=line(1:idx-1)
          READ(line(idx+1:),*,ERR=102) jndbs
          EXIT
        ENDIF
      ENDDO
!
! Extract variable index from key (if any)
      line=key
      len_line=len_trim(line)
      DO idx=1,len_line
        IF (line(idx:idx)=='-') THEN
          key=line(1:idx-1)
          READ(line(idx+1:),*,ERR=102) indvar
          EXIT
        ENDIF
      ENDDO
!
      RETURN
!
! --- error management
!
 101  WRITE (texterror,*) 'Line without "=" in configuration file'
      CALL printerror2(0,101,3,'readlist','readnextline',comment=texterror)
 102  WRITE (texterror,*) 'Bad var index in config file keyword: ',key
      CALL printerror2(0,102,3,'readlist','readnextline',comment=texterror)
!
      END
