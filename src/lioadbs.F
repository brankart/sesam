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
C ---                   LIOADBS.F                               ---
C ---                                                           ---
C --- original     : 2002-02 (J.M. Brankart)                    ---
C --- modification : 2003-03 (J.M. Brankart)                    ---
C --- modification : 2010-04 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE  readadbs    : Read ASCII observation database file
C --- SUBROUTINE  evalhdradbs : Read header of ASCII database file
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE lioadbs
      use mod_main
      use utilfiles
      IMPLICIT NONE
      PRIVATE

      PUBLIC readadbs,evalhdradbs

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readadbs (kfnindbs,kvectdbs,kjobs,
     $                     spvaldbsout,kgridij,kgridijk,kgrid)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read ASCII observation database file
CCC  -------
CCC  Method :  Read data extraction criterion from configuration file
CCC  ------    Read and select observations from ASCII database records
CCC
CCC  Input :   kfnindbs    : filename
CCC  -----     dbsobsnam   : variable name (not used)
CCC            spvaldbsout : special value
CCC
CCC  Output :  kvectdbs    : observation values
CCC  ------    kgridij     : observation location (2D case)
CCC            kgridijk    : observation location (3D case)
CCC            kgrid       : observation location (4D case)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnindbs
      BIGREAL, dimension(:), intent(out) :: kvectdbs
      INTEGER, intent(in) :: kjobs
      BIGREAL, intent(out) :: spvaldbsout
      TYPE (type_gridij), optional, dimension (:) :: kgridij
      TYPE (type_gridijk), optional, dimension (:) :: kgridijk
      TYPE (type_grid4d), optional, dimension (:) :: kgrid
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER, parameter :: numfila=99
      INTEGER :: jdbs,jpdbs,jpdbssize,jo,jpo,dtadim,indobs,inddbs
      BIGREAL :: tim_min,tim_max,lon_min,lon_max,dbs_min,dbs_max
      BIGREAL :: lat_min,lat_max,dep_min,dep_max,dbsexcl
      BIGREAL :: olon,olat,odep,otim,oval
      LOGICAL :: inside
      CHARACTER(len=1) :: textexclusion
      CHARACTER(len=word80) :: kform
      CHARACTER(len=hgword) :: line
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
     $    '*** ROUTINE : /sesam/modobsv/mkdbstoobs/mkadbstoobs/readadbs'
         WRITE(numout,*) '    ==> READING file ',kfnindbs(1:lenv(kfnindbs))
      ENDIF
C
C Get database array size
      jpdbssize=size(kvectdbs,1)
C
C Get observation and database index
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
C
C Check grid type
      IF (present(kgridij)) THEN
        IF (present(kgridijk)) GOTO 1000
        IF (present(kgrid)) GOTO 1000
        IF (jpdbssize.NE.size(kgridij,1)) GOTO 102
        dtadim = 2
      ENDIF
C
      IF (present(kgridijk)) THEN
        IF (present(kgridij)) GOTO 1000
        IF (present(kgrid)) GOTO 1000
        IF (jpdbssize.NE.size(kgridijk,1)) GOTO 102
        dtadim = 3
      ENDIF
C
      IF (present(kgrid)) THEN
        IF (present(kgridij)) GOTO 1000
        IF (present(kgridijk)) GOTO 1000
        IF (jpdbssize.NE.size(kgrid,1)) GOTO 102
        dtadim = 4
      ENDIF
C
C -1.- Get data extraction parameters
C -----------------------------------
C dbsexcl   :  observation exclusion value
C dbs_min   :  minimum observation value
C dbs_max   :  maximum observation value
C lon_min   :  minimum observation longitude
C lon_max   :  maximum observation longitude
C lat_min   :  minimum observation latitude
C lat_max   :  maximum observation latitude
C dep_min   :  minimum observation depth
C dep_max   :  maximum observation depth
C tim_min   :  minimum observation time
C tim_max   :  maximum observation time
C
      dbsexcl=obsexcl(indobs,inddbs)
      spvaldbsout = FREAL(dbsexcl)
      dbs_min=obs_min(indobs,inddbs)
      dbs_max=obs_max(indobs,inddbs)
      lon_min=obs_lon_min(indobs,inddbs)
      lon_max=obs_lon_max(indobs,inddbs)
      lat_min=obs_lat_min(indobs,inddbs)
      lat_max=obs_lat_max(indobs,inddbs)
      dep_min=obs_dep_min(indobs,inddbs)
      dep_max=obs_dep_max(indobs,inddbs)
      tim_min=obs_tim_min(indobs,inddbs)
      tim_max=obs_tim_max(indobs,inddbs)
C
C -2.- Read database file header
C ------------------------------
C
      textexclusion='#'
      CALL openfile(numfila,kfnindbs)
      line=readnextline(numfila,textexclusion)
      READ(line,*) jpdbs
      IF (jpdbs.NE.jpdbssize) GOTO 1000
C
C -3.- Read observations
C ----------------------
C Initialize with special values
      kvectdbs(:) = spvaldbsout
C
      SELECT CASE(dtadim)
      CASE(2)
         kgridij(:)%longi = spvaldbsout
         kgridij(:)%latj = spvaldbsout
      CASE(3)
         kgridijk(:)%longi = spvaldbsout
         kgridijk(:)%latj = spvaldbsout
         kgridijk(:)%levk = spvaldbsout
      CASE(4)
         kgrid(:)%lon = spvaldbsout
         kgrid(:)%lat = spvaldbsout
         kgrid(:)%dep = spvaldbsout
         kgrid(:)%tim = spvaldbsout
      END SELECT
C
C Loop over uncommented ASCII file records
      jpo = 0
      DO jdbs=1,jpdbs
         line=readnextline(numfila,textexclusion)
         READ(line,*) olon,olat,odep,otim,oval
C
C Check if observation is inside prescribed limits
         inside=.TRUE.
         inside=inside.AND.(otim.GE.tim_min)
         inside=inside.AND.(otim.LE.tim_max)
         inside=inside.AND.(odep.GE.dep_min)
         inside=inside.AND.(odep.LE.dep_max)
         inside=inside.AND.(olon.GE.lon_min)
         inside=inside.AND.(olon.LE.lon_max)
         inside=inside.AND.(olat.GE.lat_min)
         inside=inside.AND.(olat.LE.lat_max)
C
         inside=inside.AND.(oval.GE.dbs_min)
         inside=inside.AND.(oval.LE.dbs_max)
         inside=inside.AND.(oval.NE.dbsexcl)
C
         IF (inside) THEN
C
            IF (dtangrd(indobs).EQ.4) THEN
              DO WHILE (olon.GE.FREAL4(442.))
                olon = olon - FREAL4(360)
              ENDDO
              DO WHILE (olon.LT.FREAL4(78.))
                olon = olon + FREAL4(360)
              ENDDO
            ENDIF
C
C Store observation in output arrays
            jpo = jpo + 1
            IF (jpo.GT.jpdbssize) GOTO 1000
            kvectdbs(jpo) = oval
C
            SELECT CASE(dtadim)
            CASE(2)
               kgridij(jpo)%longi = olon
               kgridij(jpo)%latj = olat
            CASE(3)
               kgridijk(jpo)%longi = olon
               kgridijk(jpo)%latj = olat
               kgridijk(jpo)%levk = odep
            CASE(4)
               kgrid(jpo)%lon = olon
               kgrid(jpo)%lat = olat
               kgrid(jpo)%dep = odep
               kgrid(jpo)%tim = otim
            END SELECT
C
         ENDIF
C
      ENDDO
C
      CLOSE(numfila)
C
C -4.- Control print
C ------------------
C
      IF (nprint.GE.4) THEN
         WRITE(numout,*) ' number of observations: ',jpo
         WRITE(numout,*) ' sample of observations:'
         SELECT CASE(dtadim)
         CASE(2)
           WRITE(numout,*) ' index lon lat value'
           kform='(i9,3e12.3)'
         CASE(3)
           WRITE(numout,*) ' index lon lat depth value'
           kform='(i9,4e12.3)'
         CASE(4)
           WRITE(numout,*) ' index lon lat depth timr value'
           kform='(i9,5e12.3)'
         END SELECT
         DO jo=1,jpo,MAX(1,jpo/25)
            SELECT CASE(dtadim)
            CASE(2)
               WRITE(numout,kform) jo, kgridij(jo)%longi,
     $            kgridij(jo)%latj, kvectdbs(jo)
            CASE(3)
               WRITE(numout,kform) jo, kgridijk(jo)%longi,
     $            kgridijk(jo)%latj, kgridijk(jpo)%levk,
     $            kvectdbs(jo)
            CASE(4)
               WRITE(numout,kform) jo, kgrid(jo)%lon,
     $            kgrid(jo)%lat, kgrid(jpo)%dep,
     $            kgrid(jo)%tim, kvectdbs(jo)
            END SELECT
         ENDDO
      ENDIF
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'lioadbs','readadbs')
 1001 CALL printerror2(0,1001,3,'lioadbs','readadbs')
C
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,1,'lioadbs','readadbs',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdradbs (kfnindbs,jpdbsout)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read header of ASCII observation database file
CCC  -------
CCC  Method : Read database size from 1st uncommented record of ASCII file
CCC  ------
CCC  Input :  kfnindbs : filename
CCC  -----
CCC  Output : jpdbsout : database size
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnindbs
      INTEGER, intent(out) :: jpdbsout
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=word80) :: kform
      CHARACTER(len=hgword) :: line
      CHARACTER(len=1), parameter :: textexclusion='#'
      INTEGER, parameter :: numfila=99
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : /sesam/evalconfig/./evalhdradbs'
         WRITE(numout,*) '    ==> READING file ',kfnindbs(1:lenv(kfnindbs))
      ENDIF
C
C -1.- Read header of .adbs file
C ------------------------------
C
      CALL openfile(numfila,kfnindbs)
      line=readnextline(numfila,textexclusion)
      READ(line,*) jpdbsout
      CLOSE(numfila)
C
C -2.- Control print
C ------------------
C
      IF (nprint.GE.2) THEN
         kform='(8x,a,i10)'
         WRITE(numout,kform) '- Database size: ',jpdbsout
      ENDIF
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'lioadbs','evalhdradbs')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE lioadbs
