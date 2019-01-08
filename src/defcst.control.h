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
CC                         DEFCST.CONTROL.H
CC
CC---------------------------------------------------------------------
CC
CC  Purpose :
CC  ---------
CC  Include file with SESAM user configuration parameters
CC  For each of the variable fields (indvar=1,nbvar) in
CC  the user defined state vector, 3 blocks (A, B and C) must exist
CC  defining respectively the configuration of the Vx, Vy
CC  and Vo objects.
CC 
CC  At the end of the file, a last block specify
CC  the number of variables from the list to take
CC  into account (varend), and for each of these
CC  varend variables, a loop defines their index
CC  in the overall list previously defined.
CC  This method allows to make an exhautsive lists of all
CC  variables of the system, and to change easily SESAM configuration,
CC  by selecting some variables from the list.
CC   
CC  Modifications :
CC  -------------
CC  original     : 97-12 (C.E. Testut)
CC  modification : 99-04 (C.E. Testut)
CC  modification : 03-02 (J.M. Brankart)
CC  modification : 07-11 (J.M. Brankart)
CC
C---------------------------------------------------------------------
C Description of the configuration variables
C---------------------------------------------------------------------
C Block A: Vx object
C---------------------------------------------------------------------
C var_nam : Name of the state variable field (a maximum of 5 characters).
C var_dim : Number of dimensions of the state variable field (1, 2, 3 or 4).
C var_ord : Order of that variable in memory
C var_moy : Offset used for centering the variables (if lmoyect=.TRUE.).
C           It is substracted from the data after they are read and
C           added to the data before they are written.
C var_ect : Scale factor used for reducing the variables (if lmoyect=.TRUE.).
C           The data are divided by this factor after they are read and
C           multiplied by this factor before they are written.
C varinam : Character string inserted in the input file name
C           (in place of the `#' character), when the Vx object
C           is splitted in an ensemble of files according to the
C           variable field (the default for the 'bimg' and 'dimg' formats).
C varipos : Two digit integer specifying (i) the number of variable fields
C           in the file 'varinam' and (ii) the index
C           of the 'var_nam' variable in this file.
C varonam : Character string inserted in the output file name
C           (in place of the `#' character), when the Vx object
C           is splitted in an ensemble of files according to the
C           variable field (the default for the 'bimg' and 'dimg' formats).
C varopos : Two digit integer specifying (i) the number of variable fields
C           in the file 'varonam' and (ii) the index
C           of the 'var_nam' variable in this file.
C varfmsk : Name of the mask file for this variable field.
C           (It must exist when running SESAM!)
C varemsk : Type of the varfmsk mask file (1 for bimg, 2 for dimg,
C           3 for cdf).
C vardmsk : Dimension of the mask (usually equal to var_dim).
C varpmsk : Two digit integer specifying (i) the number of variable fields
C           in the file varfmsk and (ii) the index
C           of the var_nam variable in this file.
C varvmsk : Special value of the mask file.
C varmsea : Set to .TRUE. if the special value means inside the system
C           and to .FALSE. if the special value means outside the system.
C varfgrd : Name of the grid file for this state variable field.
C varegrd : Type of the 'varfgrd' grid file (1 for bimg, 2 for dimg,
C           3 for cdf, 4 for nc).
C varngrd : Type of 2D grid (1 for regular grids, 2 for semi-regular grids,
C           3 for irregular grids).
C---------------------------------------------------------------------
C Block B: Vy object
C---------------------------------------------------------------------
C dta_act : Set to .TRUE. if the data section object contains
C           a portion of the corresponding var_nam} variable field.
C           If it is set to .FALSE., the var_nam variable field
C           is outside of the data section object.
C dta_nam : Name of the data section variable field (a maximum of 5 characters).
C dta_dim : Number of dimensions of the data section variable field
C           (dta_dim <= var_dim).
C           (These dimensions are always the dta_dim first dimensions
C           of the state variable field.)
C dta_moy : Offset used for centering the variables (if lmoyect=.TRUE.).
C           It is substracted from the data after they are read and
C           added to the data before they are written.
C dta_ect : Scale factor used for reducing the variables (if lmoyect=.TRUE.).
C           The data are divided by this factor after they are read and
C           multiplied by this factor before they are written.
C dtainam : Character string inserted in the input file name
C           (in place of the `#' character), when the Vy object
C           is splitted in an ensemble of files according to the
C           variable field (the default for the 'bdta' and 'dta' formats).
C dtaipos : Two digit integer specifying (i) the number of variable fields
C           in the file 'dtainam' and (ii) the index
C           of the 'dta_nam' variable in this file.
C dtaonam : Character string inserted in the output file name
C           (in place of the `#' character), when the Vy object
C           is splitted in an ensemble of files according to the
C           variable field (the default for the 'bdta' and 'dta' formats).
C dtaopos : Two digit integer specifying (i) the number of variable fields
C           in the file 'dtaonam' and (ii) the index
C           of the 'dta_nam' variable in this file.
C dtafmsk : Name of the mask file for this variable field.
C           (It must exist when running SESAM!)
C           Using this mask, the user is free to include any subset
C           of the state variable field in the data section object.
C dtaemsk : Type of the 'varfmsk' mask file (1 for bdta, 2 for dta,
C           3 for cdta).
C dtadmsk : Dimension of the mask (usually equal to dta_dim).
C dtapmsk : Two digit integer specifying (i) the number of variable fields
C           in the file 'dtafmsk' and (ii) the index
C           of the 'dta_nam' variable in this file.
C dtavmsk : Special value of the mask file.
C dtamsea : Set to .TRUE. if the special value means inside the system
C           and to .FALSE. if the special value means outside the system.
C dta_rms : Default observation error standard deviation.
C           (This can be changed interactively.)
C dtafgrd : Name of the grid file for this data section variable field.
C dtaegrd : Type of the 'dtafgrd' grid file (1 for bdta, 2 for dta,
C           3 for cdta, 4 for ncdta).
C dtangrd : Type of 2D grid (1 for regular grids, 2 for semi-regular grids,
C           3 for irregular grids).
C dtaflev : Name of the level file for this data section variable field.
C dtaelev : Type of the 'varflev' level file (1 for bdta, 2 for dta,
C           3 for cdta).
C---------------------------------------------------------------------
C Block C: Vo object
C---------------------------------------------------------------------
C obsndbs : Number of observed property associated to the variable field.
C obs_nam : Name of the observed properties (a maximum of 5 characters).
C obs_dim : Dimension of the space in which the observations are localized
C           (obs_dim <= dta_dim).
C           It is used for the generation of the observation operator.
C obs_rms : Default observation error standard deviation.
C           (This can be changed interactively.)
C obs_moy : Offset used for centering the variables (if lmoyect=.TRUE.).
C           It is substracted from the data after they are read and
C           added to the data before they are written.
C           (In this example, it is set equal to var_moy.)
C obs_ect : Scale factor used for reducing the variables (if lmoyect=.TRUE.).
C           The data are divided by this factor after they are read and
C           multiplied by this factor before they are written.
C           (In this example, it is set equal to var_ect.)
C obsinam : Character string inserted in the input file name
C           (in place of the `#' character), when the Vo object
C           is splitted in an ensemble of files according to the
C           observed property (like in the obs format).
C obsonam : Character string inserted in the output file name
C           (in place of the `#' character), when the Vo object
C           is splitted in an ensemble of files according to the
C           observed property (like in the obs format).
C---------------------------------------------------------------------
C Nesting configuration
C---------------------------------------------------------------------
C nestend : number of nested grids
C nestxori, nestyori : coordinates of the nested grid south-west corner in grid 0
C nestxres, nestyres : nested grid resolution factor
C varnest : grid level for each Vx variable
C---------------------------------------------------------------------
C Following is the default SESAM configuration
C that can be changed directly here or dynamically
C to the configuration file
C---------------------------------------------------------------------
C
      DO indvar=1,nbvar
C ----------------------------------------------------------------------
C
C ==> A State vector configuration
C
C ==> A.1. Configuration of the variable
         WRITE(var_nam(indvar),'("VAR",i2.2)') indvar
         var_ord(indvar)= indvar
         var_dim(indvar)= 3
         varnest(indvar)= 0
         var_moy(indvar)= FREAL(0.0)
         var_ect(indvar)= FREAL(1.0)
C ==> A.2. Configuration of the input/ouput files
         varinam(indvar)= var_nam(indvar)
         varifil(indvar)= var_nam(indvar)
         varipos(indvar)= 0101
         varonam(indvar)= var_nam(indvar)
         varofil(indvar)= var_nam(indvar)
         varopos(indvar)= 0101
         varxdim(indvar)= 'lon'
         varydim(indvar)= 'lat'
         varzdim(indvar)= 'depth'
         vartdim(indvar)= 'time'
C ==> A.3. Configuration of the mask file
         varfmsk(indvar)= 'maskvar.cdf'
         varemsk(indvar)= 3
         vardmsk(indvar)= var_dim(indvar)
         varpmsk(indvar)= 11
         varvmsk(indvar)= FREAL(-9999.0)
         varmsea(indvar)= .FALSE.
C ==> A.4. Configuration of grid files
         varfgrd(indvar)= 'geom.cdf'
         varegrd(indvar)= 3
         varngrd(indvar)= 2
         varxnam(indvar)= 'lonxy'
         varynam(indvar)= 'latxy'
         varznam(indvar)= 'depth'
         vartnam(indvar)= 'time'
C
C ==> B Data section vector configuration
C
C ==> B.0. Does the data section include this variable?
         dta_act(indvar)= .FALSE.
C ==> B.1. Configuration of the data section variable
         dta_nam(indvar)= var_nam(indvar)
         dta_dim(indvar)= var_dim(indvar)
         dta_moy(indvar)= var_moy(indvar)
         dta_ect(indvar)= var_ect(indvar)
C ==> B.2. Configuration of the input/ouput files
         dtainam(indvar)= dta_nam(indvar)
         dtaifil(indvar)= dta_nam(indvar)
         dtaipos(indvar)= 11
         dtaonam(indvar)= dta_nam(indvar)
         dtaofil(indvar)= dta_nam(indvar)
         dtaopos(indvar)= 11
         dtaxdim(indvar)= varxdim(indvar)
         dtaydim(indvar)= varydim(indvar)
         dtazdim(indvar)= varzdim(indvar)
         dtatdim(indvar)= vartdim(indvar)
C ==> B.3. Configuration of the mask file
         dtafmsk(indvar)= 'maskdta.cdta'
         dtaemsk(indvar)= 3
         dtapmsk(indvar)= 11
         dtadmsk(indvar)= dta_dim(indvar)
         dtavmsk(indvar)= FREAL(-9999.0)
         dtamsea(indvar)= .FALSE.
C ==> B.4. Observation error standard deviation
         dta_rms(indvar)= FREAL(1.0)
C ==> B.5. Configuration of grid and level files
         dtafgrd(indvar)= 'geom.cdf'
         dtaegrd(indvar)= 3
         dtangrd(indvar)= 2
         dtaxnam(indvar)= 'lonxy'
         dtaynam(indvar)= 'latxy'
         dtaznam(indvar)= 'depth'
         dtatnam(indvar)= 'time'
         dtaflev(indvar)= 'lev.cdf'
         dtaelev(indvar)= 3
C
C ==> C Observation vector configuration
C
         obsndbs(indvar)=1 
C ==> C.1. Configuration of observation variables
         DO jndbs=1,jpndbs
           WRITE(obs_nam(indvar,jndbs),'("O",i2.2,"_",i1.1)') indvar,jndbs
           obs_dim(indvar,jndbs)= dta_dim(indvar)
           obs_moy(indvar,jndbs)= var_moy(indvar)
           obs_ect(indvar,jndbs)= var_ect(indvar)
C ==> C.2. Configuration of the input/ouput files
           obsinam(indvar,jndbs)= obs_nam(indvar,jndbs)
           obsonam(indvar,jndbs)= obs_nam(indvar,jndbs)
C ==> C.3. Observation error standard deviation
           obs_rms(indvar,jndbs)= dta_rms(indvar)
C ==> C.4. Configuration of the database files (.ncdbs)
           obsifil(indvar,jndbs)= dtaifil(indvar)
           obsxdim(indvar,jndbs)= 'lon'
           obsxnam(indvar,jndbs)= 'lon'
           obsydim(indvar,jndbs)= 'lat'
           obsynam(indvar,jndbs)= 'lat'
           obszdim(indvar,jndbs)= 'depth'
           obsznam(indvar,jndbs)= 'depth'
           obstdim(indvar,jndbs)= 'time'
           obstnam(indvar,jndbs)= 'time'
C ==> C.5. Configuration of the extraction from the database
           obsexcl(indvar,jndbs)= -HUGE(obsexcl(indvar,jndbs))
           obs_min(indvar,jndbs)= -HUGE(obs_min(indvar,jndbs))
           obs_max(indvar,jndbs)= HUGE(obs_max(indvar,jndbs))
           obsimin(indvar,jndbs)= 1
           obsimax(indvar,jndbs)= HUGE(obsimax(indvar,jndbs))
           obsjmin(indvar,jndbs)= 1
           obsjmax(indvar,jndbs)= HUGE(obsjmax(indvar,jndbs))
           obskmin(indvar,jndbs)= 1
           obskmax(indvar,jndbs)= HUGE(obskmax(indvar,jndbs))
           obstmin(indvar,jndbs)= 1
           obstmax(indvar,jndbs)= HUGE(obstmax(indvar,jndbs))
           obs_lon_min(indvar,jndbs)= -HUGE(obs_lon_min(indvar,jndbs))
           obs_lon_max(indvar,jndbs)= HUGE(obs_lon_max(indvar,jndbs))
           obs_lat_min(indvar,jndbs)= -HUGE(obs_lat_min(indvar,jndbs))
           obs_lat_max(indvar,jndbs)= HUGE(obs_lat_max(indvar,jndbs))
           obs_dep_min(indvar,jndbs)= -HUGE(obs_dep_min(indvar,jndbs))
           obs_dep_max(indvar,jndbs)= HUGE(obs_dep_max(indvar,jndbs))
           obs_tim_min(indvar,jndbs)= -HUGE(obs_tim_min(indvar,jndbs))
           obs_tim_max(indvar,jndbs)= HUGE(obs_tim_max(indvar,jndbs))
           obs_siz_max(indvar,jndbs)= 0
         ENDDO
C ----------------------------------------------------------------------
      ENDDO
C ----------------------------------------------------------------------
C Number of variables effectively taken into account
      varend = 1
C
C Next variables are excluded from SESAM vector objects
      DO indvar=varend+1,nbvar
         var_ord(indvar)=0
      ENDDO
C
C Nesting configuration
      nestend = 0
      nestxori(0:nbnest)=1
      nestyori(0:nbnest)=1
      nestxres(0:nbnest)=1
      nestyres(0:nbnest)=1
C
C Only first variable is susceptible to be observed
      dta_act(1)= .TRUE.
C ----------------------------------------------------------------------
