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
C ---                    UTILOBSV.F                             ---
C ---                                                           ---
C --- modification : 99-11 (J.M. Brankart)                      ---
C --- modification : 00-03 (J.M. Brankart)                      ---
C --- modification : 07-11 (J.M. Brankart)                      ---
C ---                                                          ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE crglocate
C --- SUBROUTINE crginterp
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilobs
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC crglocate,crginterp

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE crglocate(jicrg,jjcrg,ricrg,rjcrg,kgridijkobs,
     $                     kloncrgbias,klatcrgbias)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : locate observation 2D location in regular hor. grid
CCC  -------
CCC  Method :  dichotomy in longitude and latitude 1D arrays
CCC  ------
CCC  Input :   kgridijkobs  : observation location
CCC  -----     kloncrgbias  : longitude grid
CCC            klatcrgbias  : latitude grid
CCC  Output :  jicrg, jjcrg : SW grid location
CCC  ------    ricrg, rjcrg : postion in the grid cell (0-1)
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      INTEGER, intent(out) :: jicrg,jjcrg
      BIGREAL, intent(out) :: ricrg,rjcrg
      TYPE (type_gridijk), intent(in) :: kgridijkobs
      BIGREAL, dimension(:), intent(in) :: kloncrgbias,klatcrgbias
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL :: lon,lat
      INTEGER :: jpicrg,jpjcrg,il,ih,im
CC----------------------------------------------------------------------
C
C Get size of input arrays
      jpicrg=size(kloncrgbias,1)
      jpjcrg=size(klatcrgbias,1)
C
C Get observation location
      lon = kgridijkobs%longi
      lat = kgridijkobs%latj
C
C Dichotomy in longitude
      IF ( (lon.GE.kloncrgbias(1)) .AND. (lon.LE.kloncrgbias(jpicrg)) ) THEN
        il = 1
        ih = jpicrg
        DO WHILE (il.LT.ih-1)
          im = (il+ih)/2
          IF (lon.GT.kloncrgbias(im)) THEN
            il = im
          ELSE
            ih = im
          ENDIF
        ENDDO
        jicrg = il
        ricrg = ( lon - kloncrgbias(il) ) / ( kloncrgbias(il+1) - kloncrgbias(il) )
      ELSE
        jicrg = 0
        ricrg = FREAL(0.0)
      ENDIF
C
C Dichotomy in latitude
      IF ( (lat.GE.klatcrgbias(1)) .AND. (lat.LE.klatcrgbias(jpjcrg)) ) THEN
        il = 1
        ih = jpjcrg
        DO WHILE (il.LT.ih-1)
          im = (il+ih)/2
          IF (lat.GT.klatcrgbias(im)) THEN
            il = im
          ELSE
            ih = im
          ENDIF
        ENDDO
        jjcrg = il
        rjcrg = ( lat - klatcrgbias(il) ) / ( klatcrgbias(il+1) - klatcrgbias(il) )
      ELSE
        jjcrg = 0
        rjcrg = FREAL(0.0)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilobsv','crglocate')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE crginterp(kobsbias,kmatcrg,kspvalcrg,
     $                     kjicrg,kjjcrg,kricrg,krjcrg)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Interpolate in 2D observation array
CCC  -------
CCC  Method :  Bilinear interpolation
CCC  ------
CCC  Input :   kjicrg, kjjcrg : SW grid location
CCC  -----     kricrg, krjcrg : postion in the grid cell (0-1)
CCC            kspvalcrg      : special value in observation array
CCC            kmatcrg        : 2D observation array
CCC  Output :  kobsbias       : interpolated observation value
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      INTEGER, intent(in) :: kjicrg,kjjcrg
      BIGREAL, intent(in) :: kricrg,krjcrg
      BIGREAL, intent(in) :: kspvalcrg
      BIGREAL, intent(out) :: kobsbias
      BIGREAL, dimension(:,:), intent(in) :: kmatcrg
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL :: lon,lat,w00,w10,w01,w11
      INTEGER :: ji,jj,jpicrg,jpjcrg
      LOGICAL :: bad
CC----------------------------------------------------------------------
C
C Get size of input arrays
      jpicrg=size(kmatcrg,1)
      jpjcrg=size(kmatcrg,2)
CC----------------------------------------------------------------------
C
C Check if all 4 values of the grid cell are not special value
      bad = .FALSE.
      bad = bad .OR. (kjicrg.LE.0)
      bad = bad .OR. (kjicrg.GT.jpicrg)
      bad = bad .OR. (kjjcrg.LE.0)
      bad = bad .OR. (kjjcrg.GT.jpjcrg)
      IF (.not.bad) THEN
         bad = bad .OR. (kmatcrg(kjicrg,kjjcrg).EQ.kspvalcrg)
         bad = bad .OR. (kmatcrg(kjicrg+1,kjjcrg).EQ.kspvalcrg)
         bad = bad .OR. (kmatcrg(kjicrg,kjjcrg+1).EQ.kspvalcrg)
         bad = bad .OR. (kmatcrg(kjicrg+1,kjjcrg+1).EQ.kspvalcrg)
      ENDIF
C
C Bilinear interpolation
      IF (bad) THEN
         kobsbias = kspvalcrg
      ELSE
         w11 = kricrg * krjcrg
         w10 = kricrg * (FREAL(1.0) - krjcrg)
         w01 = (FREAL(1.0) - kricrg) * krjcrg
         w00 = (FREAL(1.0) - kricrg) * (FREAL(1.0) - krjcrg)
         kobsbias = FREAL(0.0)
         kobsbias = kobsbias + kmatcrg(kjicrg,kjjcrg) * w00
         kobsbias = kobsbias + kmatcrg(kjicrg+1,kjjcrg) * w10
         kobsbias = kobsbias + kmatcrg(kjicrg,kjjcrg+1) * w01
         kobsbias = kobsbias + kmatcrg(kjicrg+1,kjjcrg+1) * w11
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilobsv','crginterp')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilobs
