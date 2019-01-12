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
! ---                    UTILOBSV.F90                             ---
! ---                                                           ---
! --- modification : 99-11 (J.M. Brankart)                      ---
! --- modification : 00-03 (J.M. Brankart)                      ---
! --- modification : 07-11 (J.M. Brankart)                      ---
! ---                                                          ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE crglocate
! --- SUBROUTINE crginterp
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilobs
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC crglocate,crginterp

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE crglocate(jicrg,jjcrg,ricrg,rjcrg,kgridijkobs, &
     &                     kloncrgbias,klatcrgbias)
!---------------------------------------------------------------------
!
!  Purpose : locate observation 2D location in regular hor. grid
!  -------
!  Method :  dichotomy in longitude and latitude 1D arrays
!  ------
!  Input :   kgridijkobs  : observation location
!  -----     kloncrgbias  : longitude grid
!            klatcrgbias  : latitude grid
!  Output :  jicrg, jjcrg : SW grid location
!  ------    ricrg, rjcrg : postion in the grid cell (0-1)
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, intent(out) :: jicrg,jjcrg
      BIGREAL, intent(out) :: ricrg,rjcrg
      TYPE (type_gridijk), intent(in) :: kgridijkobs
      BIGREAL, dimension(:), intent(in) :: kloncrgbias,klatcrgbias
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL :: lon,lat
      INTEGER :: jpicrg,jpjcrg,il,ih,im
!----------------------------------------------------------------------
!
! Get size of input arrays
      jpicrg=size(kloncrgbias,1)
      jpjcrg=size(klatcrgbias,1)
!
! Get observation location
      lon = kgridijkobs%longi
      lat = kgridijkobs%latj
!
! Dichotomy in longitude
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
!
! Dichotomy in latitude
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
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilobsv','crglocate')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE crginterp(kobsbias,kmatcrg,kspvalcrg, &
     &                     kjicrg,kjjcrg,kricrg,krjcrg)
!---------------------------------------------------------------------
!
!  Purpose : Interpolate in 2D observation array
!  -------
!  Method :  Bilinear interpolation
!  ------
!  Input :   kjicrg, kjjcrg : SW grid location
!  -----     kricrg, krjcrg : postion in the grid cell (0-1)
!            kspvalcrg      : special value in observation array
!            kmatcrg        : 2D observation array
!  Output :  kobsbias       : interpolated observation value
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, intent(in) :: kjicrg,kjjcrg
      BIGREAL, intent(in) :: kricrg,krjcrg
      BIGREAL, intent(in) :: kspvalcrg
      BIGREAL, intent(out) :: kobsbias
      BIGREAL, dimension(:,:), intent(in) :: kmatcrg
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL :: lon,lat,w00,w10,w01,w11
      INTEGER :: ji,jj,jpicrg,jpjcrg
      LOGICAL :: bad
!----------------------------------------------------------------------
!
! Get size of input arrays
      jpicrg=size(kmatcrg,1)
      jpjcrg=size(kmatcrg,2)
!----------------------------------------------------------------------
!
! Check if all 4 values of the grid cell are not special value
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
!
! Bilinear interpolation
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
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilobsv','crginterp')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilobs
