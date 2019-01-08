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
C ---                   ALGODIFF.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-11 (C.E. Testut)                        ---
C --- modification : 01-06 (C.E. Testut)                        ---
C --- modification : 06-03 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE algodiffbyzon
C --- SUBROUTINE algodiffbylev
C --- SUBROUTINE algomeanstdbyzon
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algodiff
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC algodiffbyzon,algodiffbylev,algomeanstdbyzon

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algodiffbyzon(karginxyo,kargdiffxyoref,
     $     kargdiffxyoorg,klargdiffxyoorg,kflaganlxyo,knam_grouparea,
     $        ktab_nbarea,ktab_grouparea)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Algorithm for difference calcul on state :
CCC  -------
CCC  Method :
CCC  ------
CCC  Input :
CCC  -----
CCC  Output :
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : spvalvar,spvaldta,spvalobs,
     $     jpx,jpxend,jpo,jpitpend,jpy,poscoefobs,
     $     arraynx_jindxbeg,arraynx_jpindxend
      use hioxyo
      use utilvalid
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      INTEGER, intent(in) :: kflaganlxyo
      LOGICAL, intent(in) :: klargdiffxyoorg
      CHARACTER(len=*), intent(in) :: karginxyo,kargdiffxyoref,
     $     kargdiffxyoorg
      CHARACTER(len=*), dimension(:), intent(in) :: knam_grouparea
      INTEGER, dimension(:), intent(in) :: ktab_nbarea
      INTEGER, dimension(:,:), intent(in) :: ktab_grouparea
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:), allocatable, save :: vectsref
      BIGREAL, dimension(:), allocatable, save :: vectszon
      BIGREAL, dimension(:), allocatable, save :: vectyzon
C
      INTEGER :: allocok,jpssize,jpysize,jpitpsize
      INTEGER :: jpgroupsize,jpareasize,jgroup,jarea
      LOGICAL :: lectinfo,lmoyectold
      INTEGER :: jnxyo,flagxyo1,flagcfg
      INTEGER :: jsend,jindsbeg,jindsend,kjzondeb,kjzonfin
      INTEGER :: xyoend,jxyo,indxyo,inddbs
      INTEGER :: js,jsdeb,jsfin,jind,kjzon
      INTEGER :: tabxyoend,tabxyoend1,jtabxyo
      INTEGER, dimension(1:1) :: jitp
      INTEGER, dimension(:), allocatable :: xyo_ind,xyo_nbr
      INTEGER, dimension(:,:), allocatable :: xyo_nbrgrp
      BIGREAL, dimension(:), allocatable :: xyo_ect
      BIGREAL, dimension(:,:), allocatable :: xyo_RMS_AR
      CHARACTER(len=varlg), dimension(:), allocatable :: xyo_nam
      INTEGER, dimension (:,:), allocatable :: tabnbrxyo
      BIGREAL, dimension(:,:), allocatable :: RMS_ARxyo
      BIGREAL, dimension(:), allocatable :: vcts_sr
      BIGREAL :: spval1,spval
      INTEGER :: labtitr,labnbre
      CHARACTER(len=bgword) :: ligne,lignbre,lignbre1,
     $     methode,fhisto1,fhisto2
      CHARACTER(len=bgword) :: myformat
      INTEGER, parameter  :: zonkfreq = 1
CC----------------------------------------------------------------------
C
      jpysize=jpy
      jpgroupsize=size(ktab_grouparea,1)
      jpareasize=size(ktab_grouparea,2)
      SELECT CASE (kflaganlxyo)
      CASE (1)
         jpssize=jpx
         jpitpsize=1
         xyoend  = varend
      CASE (2)
         jpssize=jpy
         jpitpsize=1
         xyoend  = dtaend
      CASE (3)
         jpssize=jpo
         jpitpsize=jpitpend
         xyoend  = obsend
      CASE DEFAULT
         GOTO 1000
      END SELECT
C --- allocation xyo_nbr
      allocate ( xyo_nbr(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_nbr(:) = 0
C --- allocation xyo_nbrgrp
      allocate ( xyo_nbrgrp(1:xyoend,1:(jpgroupsize+1)), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_nbrgrp(:,:) = 0
C --- allocation xyo_ind
      allocate ( xyo_ind(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_ind(:) = 0
C --- allocation xyo_ect
      allocate ( xyo_ect(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_ect(:) = FREAL(1.0)
C --- allocation xyo_nam
      allocate ( xyo_nam(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_nam(:) = ' '
CCC---------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine algodiffbyzon &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*)
      ENDIF
C
C -0.- Initialisation :
C ---------------------
C
      IF ((ktab_nbarea(1).GT.0).AND.
     $     (ktab_nbarea(1).LE.jpareasize)) THEN
         kjzondeb = ktab_grouparea(1,1)
         kjzonfin = ktab_grouparea(1,1)
         DO kjzon=2,ktab_nbarea(1)
            kjzondeb = MIN(kjzondeb,ktab_grouparea(1,kjzon))
            kjzonfin = MAX(kjzonfin,ktab_grouparea(1,kjzon))
         ENDDO
      ELSE
         kjzondeb=1
         kjzonfin=1
      ENDIF
C
      DO jgroup=2,jpgroupsize
         IF ((ktab_nbarea(jgroup).GT.0).AND.
     $        (ktab_nbarea(jgroup).LE.jpareasize)) THEN
         DO kjzon=1,ktab_nbarea(jgroup)
            kjzondeb = MIN(kjzondeb,ktab_grouparea(jgroup,kjzon))
            kjzonfin  =MAX(kjzonfin,ktab_grouparea(jgroup,kjzon))
         ENDDO
         ENDIF
      ENDDO
C
      spval1=spvalvar
C --- allocation poscoefobs
      allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      SELECT CASE (kflaganlxyo)
      CASE (1)
         DO js=1,jpssize
            poscoefobs(js,:) = type_poscoef(js,FREAL(1.0))
         ENDDO
      CASE (2)
         DO js=1,jpssize
            poscoefobs(js,:) = type_poscoef(js,FREAL(1.0))
         ENDDO
      CASE (3)
         poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
C --- reading config.obs
         flagcfg=3
         CALL readcfgobs (argconfigobs,flagcfg,
     $        kposcoefobs=poscoefobs(:,:))
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C --- allocation vects
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
C --- allocation vectsref
      allocate ( vectsref(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectsref(:) = FREAL(0.0)
C --- allocation vectszon
      allocate ( vectszon(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectszon(:) = FREAL(0.0)
C
      SELECT CASE (kflaganlxyo)
      CASE (1,2)
C --- nothing
      CASE (3)
C --- allocation vectyzon
         allocate ( vectyzon(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectyzon(:) = FREAL(0.0)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      SELECT CASE (kflaganlxyo)
      CASE(1)
         DO jxyo = 1,xyoend
            indxyo= var_ord(jxyo)
            xyo_nam(jxyo)=var_nam(indxyo)
            xyo_ind(jxyo)=var_ind(indxyo)
            xyo_nbr(jxyo)=var_nbr(indxyo)
            xyo_ect(jxyo)=var_ect(indxyo)
         ENDDO
      CASE(2)
         DO jxyo = 1,xyoend
            indxyo=dta_ord(jxyo)
            xyo_nam(jxyo)=dta_nam(indxyo)
            xyo_ind(jxyo)=dta_ind(indxyo)
            xyo_nbr(jxyo)=dta_nbr(indxyo)
            xyo_ect(jxyo)=dta_ect(indxyo)
         ENDDO
      CASE(3)
         DO jxyo = 1,xyoend
            indxyo=obs_ord(jxyo)
            inddbs=obsnord(jxyo)
            xyo_nam(jxyo)=obs_nam(indxyo,inddbs)
            xyo_ind(jxyo)=obs_ind(indxyo,inddbs)
            xyo_nbr(jxyo)=obs_nbr(indxyo,inddbs)
            xyo_ect(jxyo)=obs_ect(indxyo,inddbs)
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      tabxyoend=(kjzonfin-kjzondeb+1)
      tabxyoend1=tabxyoend+1
C
C -1.- Initialisation :
C ---------------------
C
C --- allocation tabnbrxyo
      allocate ( tabnbrxyo(0:tabxyoend1,1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabnbrxyo(:,:) = 0
C --- allocation RMS_ARxyo
      allocate ( RMS_ARxyo(0:tabxyoend1,1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      RMS_ARxyo(:,:) = FREAL(0.0)
C --- allocation xyo_RMS_AR
      allocate ( xyo_RMS_AR(1:xyoend,1:(1+jpgroupsize)), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_RMS_AR(:,:) = FREAL(0.0)
C
C -2.- Run the DIFF analysis :
C ----------------------------
C
      DO jnxyo=1,limjpnxyo(kflaganlxyo)
         IF (kflaganlxyo.EQ.1) THEN
            jindsbeg=arraynx_jindxbeg(jnxyo)
            jindsend=arraynx_jindxbeg(jnxyo)-1+arraynx_jpindxend(jnxyo)
            jsend=arraynx_jpindxend(jnxyo)
         ELSE
            jindsbeg=1
            jindsend=jpssize
            jsend=jpssize
         ENDIF
C --- reading object diffxyoref
         lectinfo=.TRUE.
         CALL readxyo(kargdiffxyoref,vectsref(:),
     $        jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
         IF (validextvar(kargdiffxyoref)) THEN
            spval=spvalvar
         ELSEIF (validextdta(kargdiffxyoref)) THEN
            spval=spvaldta
         ELSEIF (validextobs(kargdiffxyoref)) THEN
            spval=spvalobs
         ENDIF
         IF (ANY(vectsref(:).EQ.spval)) GOTO 101
         lectinfo=.FALSE.
         IF (largbias) THEN
            print *,'activate bias : diffref + bias'
            CALL readxyo(argbias,vects(:),
     $           jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:)) 
            IF (validextvar(argbias)) THEN
               spval=spvalvar
            ELSEIF (validextdta(argbias)) THEN
               spval=spvaldta
            ELSEIF (validextobs(argbias)) THEN
               spval=spvalobs
            ENDIF
            IF (ANY(vects(:).EQ.spval)) GOTO 102
            vectsref(:) = vectsref(:) + vects(:)
         ENDIF
C --- reading object inxyo
         CALL readxyo(karginxyo,vects(:),
     $     jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:)) 
         IF (validextvar(karginxyo)) THEN
            spval=spvalvar
         ELSEIF (validextdta(karginxyo)) THEN
            spval=spvaldta
         ELSEIF (validextobs(karginxyo)) THEN
            spval=spvalobs
         ENDIF
         IF (ANY(vects(:).EQ.spval)) GOTO 103
C --- reading object inpartvar
         lmoyectold=lmoyect
         lmoyect=.FALSE.
         SELECT CASE (kflaganlxyo)
         CASE (1,2)
            CALL readxyo(arginpartvar,vectszon(:),
     $           jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:)) 
         CASE (3)
            flagxyo1=2
            CALL readxyo(arginpartvar,vectyzon(:),
     $           jnxyo,lectinfo,flagxyo1,poscoefobs(:,:)) 
C           vectszon(:jsend) = vectyzon(poscoefobs(:jsend,1)%pos)
            DO js=1,jsend
              jitp(1:1)=MAXLOC(poscoefobs(js,1:jpitpend)%coef)
              vectszon(js) = vectyzon(poscoefobs(js,jitp(1))%pos)
            ENDDO
         CASE DEFAULT
            GOTO 1000
         END SELECT
         lmoyect=lmoyectold
         vectszon(1:jsend)=FREAL(NINT(vectszon(1:jsend)))
         DO js=1,jsend
            vectszon(js)=vectszon(js)-FREAL(kjzondeb-1)
            vectszon(js)=MAX(FREAL(0.0),vectszon(js))
            vectszon(js)=MIN(vectszon(js),FREAL(kjzonfin-kjzondeb+1+1))
         ENDDO
#if defined _NEC
CCC!CDIR NOVECTOR
         IF (ALL(NINT(vectszon(1:jsend)).GT.tabxyoend1)) GOTO 1000
CCC         IF (INT(MAXVAL(vectszon(1:jsend))).GT.tabxyoend1) GOTO 1000
CCC!CDIR NOVECTOR
         IF (ALL(NINT(vectszon(1:jsend)).LT.0)) GOTO 1000
CCC         IF (INT(MINVAL(vectszon(1:jsend))).LT.0) GOTO 1000
#else
         IF (INT(MAXVAL(vectszon(1:jsend))).GT.tabxyoend1) GOTO 1000
         IF (INT(MINVAL(vectszon(1:jsend))).LT.0) GOTO 1000
#endif
C
C --- allocation vcts_sr
         allocate ( vcts_sr(1:jsend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vcts_sr(:) = FREAL(0.0)
         vcts_sr(1:jsend) = vects(1:jsend)-vectsref(1:jsend)
C
         DO jxyo=1,xyoend
            jsdeb = MAX(1,xyo_ind(jxyo)-(jindsbeg-1))
            jsfin = MIN(jsend,xyo_nbr(jxyo)+xyo_ind(jxyo)-1-(jindsbeg-1))
            DO js=jsdeb,jsfin
               RMS_ARxyo(INT(vectszon(js)),jxyo) = 
     $              RMS_ARxyo(INT(vectszon(js)),jxyo) +
     $              vcts_sr(js)*vcts_sr(js)
            ENDDO
C
            DO js=jsdeb,jsfin
               tabnbrxyo(INT(vectszon(js)),jxyo) = 
     $              tabnbrxyo(INT(vectszon(js)),jxyo) + 1
            ENDDO
         ENDDO
C
         IF (allocated(vcts_sr)) deallocate(vcts_sr)
C
      ENDDO
C
      DO jxyo=1,xyoend
      DO jgroup=1,jpgroupsize
      DO jarea=1,ktab_nbarea(jgroup)
         jtabxyo=ktab_grouparea(jgroup,jarea)-(kjzondeb-1)
         IF (tabnbrxyo(jtabxyo,jxyo).GT.0) THEN
            xyo_RMS_AR(jxyo,jgroup+1) =
     $           xyo_RMS_AR(jxyo,jgroup+1)+RMS_ARxyo(jtabxyo,jxyo)
            xyo_nbrgrp(jxyo,jgroup+1) = 
     $           xyo_nbrgrp(jxyo,jgroup+1) + tabnbrxyo(jtabxyo,jxyo)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
#if defined _NEC
!CDIR NOVECTOR
      DO jxyo=1,xyoend
!CDIR NOVECTOR
         DO jtabxyo=1,tabxyoend
#else
      DO jxyo=1,xyoend
         DO jtabxyo=1,tabxyoend
#endif
         IF (tabnbrxyo(jtabxyo,jxyo).GT.0) THEN
            xyo_RMS_AR(jxyo,1) =
     $           xyo_RMS_AR(jxyo,1)+RMS_ARxyo(jtabxyo,jxyo)
            xyo_nbrgrp(jxyo,1) = 
     $           xyo_nbrgrp(jxyo,1) + tabnbrxyo(jtabxyo,jxyo)
            RMS_ARxyo(jtabxyo,jxyo) = FREAL(SQRT(RMS_ARxyo(jtabxyo,jxyo) /
     $           tabnbrxyo(jtabxyo,jxyo)))
            IF ((lmoyect)
     $           .AND.(xyo_ect(jxyo).NE.FREAL(1.0))) THEN
               RMS_ARxyo(jtabxyo,jxyo) = RMS_ARxyo(jtabxyo,jxyo) *
     $              xyo_ect(jxyo)
            ENDIF
         ELSE
            RMS_ARxyo(jtabxyo,jxyo)   = FREAL(0.0)
         ENDIF
      ENDDO
      ENDDO
C
      DO jxyo =1,xyoend
         xyo_RMS_AR(jxyo,1) = xyo_RMS_AR(jxyo,1) +
     $        RMS_ARxyo(0,jxyo) + RMS_ARxyo(tabxyoend1,jxyo)
         xyo_nbrgrp(jxyo,1) = xyo_nbrgrp(jxyo,1) +
     $        tabnbrxyo(0,jxyo) + tabnbrxyo(tabxyoend1,jxyo)
      ENDDO
C
      DO jind =1,jpgroupsize+1
      DO jxyo =1,xyoend
         IF (xyo_nbrgrp(jxyo,jind).NE.FREAL(0.0)) THEN
            xyo_RMS_AR(jxyo,jind) = 
     $           FREAL(SQRT(xyo_RMS_AR(jxyo,jind)/xyo_nbrgrp(jxyo,jind)))
            IF ((lmoyect).AND.(xyo_ect(jxyo).NE.FREAL(1.0))) THEN
               xyo_RMS_AR(jxyo,jind) = xyo_RMS_AR(jxyo,jind) *
     $              xyo_ect(jxyo)
            ENDIF
         ELSE
            xyo_RMS_AR(jxyo,jind) = FREAL(0.0)
         ENDIF
      ENDDO
      ENDDO
C
C -3.- Print-Control :
C --------------------
C
      WRITE(ligne,'(A,A)') 
     $        ' ------------------------------------',
     $        '-------------------------'
      ASSIGN 11 to labtitr
      ASSIGN 20 to labnbre
C
      PRINT *,'  -------------------  Array of ERRORS'
     $       ,'  -------------------'
      PRINT '(A)', ligne(1:lenv(ligne))
      PRINT labtitr, 'Var','Method','Nbr pts'
     $        ,' RMS[AR] ',' AR_R(%) '
      PRINT '(A)', ligne(1:lenv(ligne))
      IF (nprint.GE.0) THEN
         WRITE(numout,'(A)') '   --------  Array of ERRORS  --------'
         WRITE(numout,'(A)') ligne(1:lenv(ligne))
         WRITE(numout,labtitr) 'Var','Methode','Nbr pts'
     $           ,' RMS[AR] ',' AR_R(%) '
         WRITE(numout,'(A)') ligne(1:lenv(ligne))
      ENDIF
C
      DO jxyo=1,xyoend
C jind=1
         jind=1
         methode='-ALL AREAS-'
         WRITE(lignbre,labnbre) xyo_nam(jxyo),methode,
     $        xyo_nbrgrp(jxyo,jind),xyo_RMS_AR(jxyo,jind)
         IF (nprint.GE.0) WRITE(numout,'(A)') lignbre(1:lenv(lignbre))
         PRINT '(A)',lignbre(1:lenv(lignbre))
C 
         DO jgroup=1,jpgroupsize
            methode=knam_grouparea(jgroup)(1:11)
            WRITE(lignbre,labnbre) xyo_nam(jxyo),methode,
     $           xyo_nbrgrp(jxyo,jgroup+1),xyo_RMS_AR(jxyo,jgroup+1)
            SELECT CASE (ktab_nbarea(jgroup))
            CASE (7:)
               IF (nprint.GE.0) WRITE(numout,'(A)') lignbre(1:lenv(lignbre))
               PRINT '(A)',lignbre(1:lenv(lignbre))
            CASE (2:6)
               WRITE(myformat,'(A5,I1,A17)') '(A,A,',
     $              ktab_nbarea(jgroup)-1,'(I3.3,1X),I3.3,A)'
C            print *,myformat
               IF (nprint.GE.0) WRITE(numout,myformat) lignbre(1:lenv(lignbre)),'   areas=(',
     $              (ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))),')'
               PRINT myformat,lignbre(1:lenv(lignbre)),'   areas=(',
     $              (ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))),')'
            CASE (1)
               WRITE(myformat,'(A)') '(A,A,I3.3,A)'
C            print *,myformat
               IF (nprint.GE.0) WRITE(numout,myformat) lignbre(1:lenv(lignbre)),'   areas=(',
     $              (ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))),')'
               PRINT myformat,lignbre(1:lenv(lignbre)),'   areas=(',
     $              (ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))),')'
            CASE DEFAULT
               GOTO 1000
            END SELECT
         ENDDO
C
         DO jtabxyo=1,tabxyoend
            WRITE(methode,'("jzon",I3.3,"-",I3.3)') 
     $           jtabxyo+(kjzondeb-1),
     $           jtabxyo+(kjzondeb-1)
            WRITE(lignbre,labnbre) xyo_nam(jxyo),methode,
     $           tabnbrxyo(jtabxyo,jxyo),RMS_ARxyo(jtabxyo,jxyo)
            IF (nprint.GE.0) WRITE(numout,'(A)') lignbre(1:lenv(lignbre))
C            IF (MOD(jtabxyo,zonkfreq).EQ.0) 
C     $           PRINT '(A)',lignbre(1:lenv(lignbre))
         ENDDO
         PRINT '(A)', ligne(1:lenv(ligne))
         WRITE(numout,'(A)') ligne(1:lenv(ligne))
      ENDDO
C
C --- deallocate
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(vectsref)) deallocate(vectsref)
      IF (allocated(vectszon)) deallocate(vectszon)
      IF (allocated(vectyzon)) deallocate(vectyzon)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(xyo_ind)) deallocate(xyo_ind)
      IF (allocated(xyo_nbr)) deallocate(xyo_nbr)
      IF (allocated(xyo_ect)) deallocate(xyo_ect)
      IF (allocated(xyo_nam)) deallocate(xyo_nam)
      IF (allocated(tabnbrxyo)) deallocate(tabnbrxyo)
      IF (allocated(RMS_ARxyo)) deallocate(RMS_ARxyo)
C
      RETURN
C
C --- format definition
C
C --- titles
 10   FORMAT(3X,A3,3X,"|",3X,A7,3X,"|",1X,A7,1(1X,"|",2X,A9,1X))
 11   FORMAT(3X,A3,3X,"|",3X,A7,3X,"|",1X,A7,2(1X,"|",2X,A9,1X))
 12   FORMAT(3X,A3,3X,"|",3X,A7,3X,"|",1X,A7,3(1X,"|",2X,A9,1X))
C --- values  
 20   FORMAT(2X,A5,2X,"|",1X,A11,3X,I7.7,3X,1(1PE11.5E2,3X))
 21   FORMAT(2X,A5,2X,"|",1X,A11,3X,I7.7,3X,2(1PE11.5E2,3X))
 22   FORMAT(2X,A5,2X,"|",1X,A11,3X,I7.7,3X,3(1PE11.5E2,3X))
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algodiff','algodiffbyzon')
 1001 CALL printerror2(0,1001,3,'algodiff','algodiffbyzon')
C
 101  WRITE (texterror,*) 'spval exist in the file ',
     $     kargdiffxyoref(1:lenv(kargdiffxyoref)),
     $     ' => interf. it in .obs'
      CALL printerror2(0,101,3,'algodiff','algodiffbyzon',
     $     comment=texterror)
 102  WRITE (texterror,*) 'spval exist in the file ',
     $     argbias(1:lenv(argbias)),
     $     ' => interf. it in .obs'
      CALL printerror2(0,102,3,'algodiff','algodiffbyzon',
     $     comment=texterror)
 103  WRITE (texterror,*) 'spval exist in the file ',
     $     karginxyo(1:lenv(karginxyo)),
     $     ' => interf. it in .obs'
      CALL printerror2(0,103,3,'algodiff','algodiffbyzon',
     $     comment=texterror)
 104  WRITE (texterror,*) 'spval exist in the file ',
     $     kargdiffxyoorg(1:lenv(kargdiffxyoorg)),
     $     ' => interf. it in .obs'
      CALL printerror2(0,104,3,'algodiff','algodiffbyzon',
     $     comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algodiffbylev(karginxyo,kargdiffxyoref,
     $     kargdiffxyoorg,klargdiffxyoorg,kflaganlxyo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Algorithm for difference calcul on state :
CCC  -------
CCC  Method :
CCC  ------
CCC  Input :
CCC  -----
CCC  Output :
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : spvalvar,spvaldta,spvalobs,
     $     jpx,jpxend,jpo,jpitpend,jpy,poscoefobs,
     $     arraynx_jindxbeg,arraynx_jpindxend
      use hioxyo
      use utilvalid
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      INTEGER, intent(in) :: kflaganlxyo
      LOGICAL, intent(in) :: klargdiffxyoorg
      CHARACTER(len=*), intent(in) :: karginxyo,kargdiffxyoref,
     $     kargdiffxyoorg
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:), allocatable, save :: vectsref
      BIGREAL, dimension(:), allocatable, save :: vectsorg
C
      INTEGER :: allocok,jpssize,jpysize,jpitpsize
      LOGICAL :: lectinfo
      INTEGER :: jnxyo,flagxyo1,flagcfg
      INTEGER :: jlevdeb,jlevfin,jsend,jindsbeg,jindsend
      INTEGER :: xyoend,jxyo,indxyo,inddbs,indxyomsk,indxyomskdeb
      INTEGER :: ji,jj,jk,jt,js,jsdeb,jsfin,jind
      INTEGER :: tabxyoend,jtabxyo,jtabxyo1,jtabxyo2,
     $     jtabxyodeb,jtabxyofin
      INTEGER, dimension(:), allocatable :: xyo_dim,
     $     xyo_ind,xyo_jpi,xyo_jpj,xyo_jpk,xyo_jpt
      INTEGER, dimension(:), allocatable :: xyo_nbr,xyo_nbr1
      BIGREAL, dimension(:), allocatable :: xyo_ect
      BIGREAL, dimension(:), allocatable :: xyo_RMS_AR,
     $     xyo_rap_AR_R,xyo_rap_AR_LR
      CHARACTER(len=varlg), dimension(:), allocatable :: xyo_nam
      INTEGER, dimension (:), allocatable :: tabindxyo,
     $     tabnbrxyo,tabnbrxyo1,tabjxyoxyo,tabjkxyo
      BIGREAL, dimension(:), allocatable :: RMS_ARxyo,rap_AR_Rxyo,
     $     rap_AR_LRxyo
      BIGREAL, dimension(:), allocatable :: vcts_sr,vctso_sr
      BIGREAL :: spval1,spval
      INTEGER :: labtitr,labnbre
      CHARACTER(len=bgword) :: ligne,lignbre,lignbre1,
     $     methode,fhisto1,fhisto2
      INTEGER, parameter  :: zonkfreq = 1
CC----------------------------------------------------------------------
C
      SELECT CASE (kflaganlxyo)
      CASE (1)
         jpssize=jpx
         jpitpsize=1
         xyoend  = varend
      CASE (2)
         jpssize=jpy
         jpitpsize=1
         xyoend  = dtaend
      CASE (3)
         jpssize=jpo
         jpitpsize=jpitpend
         xyoend  = obsend
      CASE DEFAULT
         GOTO 1000
      END SELECT
C --- allocation xyo_dim
      allocate ( xyo_dim(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_dim(:) = 0
C --- allocation xyo_ind
      allocate ( xyo_ind(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_ind(:) = 0
C --- allocation xyo_nbr
      allocate ( xyo_nbr(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_nbr(:) = 0
C --- allocation xyo_nbr1
      allocate ( xyo_nbr1(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_nbr1(:) = 0
C --- allocation xyo_jpi
      allocate ( xyo_jpi(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_jpi(:) = 0
C --- allocation xyo_jpj
      allocate ( xyo_jpj(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_jpj(:) = 0
C --- allocation xyo_jpk
      allocate ( xyo_jpk(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_jpk(:) = 0
C --- allocation xyo_jpt
      allocate ( xyo_jpt(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_jpt(:) = 0
C --- allocation xyo_ect
      allocate ( xyo_ect(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_ect(:) = FREAL(1.0)
C --- allocation xyo_nam
      allocate ( xyo_nam(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_nam(:) = ' '
CCC---------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine algodiffbylev &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*)
      ENDIF
C
C -0.- Initialisation :
C ---------------------
C
      spval1=spvalvar
C --- allocation poscoefobs
      allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      SELECT CASE (kflaganlxyo)
      CASE (1,4)
         DO js=1,jpssize
            poscoefobs(js,:) = type_poscoef(js,FREAL(1.0))
         ENDDO
      CASE (2,5)
         DO js=1,jpssize
            poscoefobs(js,:) = type_poscoef(js,FREAL(1.0))
         ENDDO
      CASE (3,6)
         poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
C --- reading config.obs
         flagcfg=3
         CALL readcfgobs (argconfigobs,flagcfg,
     $        kposcoefobs=poscoefobs(:,:))
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C --- allocation vects
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
C --- allocation vectsref
      allocate ( vectsref(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectsref(:) = FREAL(0.0)
C --- allocation vectsorg
      IF (klargdiffxyoorg) THEN
         allocate ( vectsorg(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectsorg(:) = FREAL(0.0)
      ENDIF
C
      SELECT CASE (kflaganlxyo)
      CASE(1)
         DO jxyo = 1,xyoend
            indxyo= var_ord(jxyo)
            xyo_nam(jxyo)=var_nam(indxyo)
            xyo_dim(jxyo)=var_dim(indxyo)
            xyo_ind(jxyo)=var_ind(indxyo)
            xyo_nbr(jxyo)=var_nbr(indxyo)
            xyo_jpi(jxyo)=var_jpi(indxyo)
            xyo_jpj(jxyo)=var_jpj(indxyo)
            xyo_jpk(jxyo)=var_jpk(indxyo)
            xyo_jpt(jxyo)=var_jpt(indxyo)
            xyo_ect(jxyo)=var_ect(indxyo)
         ENDDO
      CASE(2)
         DO jxyo = 1,xyoend
            indxyo=dta_ord(jxyo)
            xyo_nam(jxyo)=dta_nam(indxyo)
            xyo_dim(jxyo)=dta_dim(indxyo)
            xyo_ind(jxyo)=dta_ind(indxyo)
            xyo_nbr(jxyo)=dta_nbr(indxyo)
            xyo_jpi(jxyo)=dta_jpi(indxyo)
            xyo_jpj(jxyo)=dta_jpj(indxyo)
            xyo_jpk(jxyo)=dta_jpk(indxyo)
            xyo_jpt(jxyo)=dta_jpt(indxyo)
            xyo_ect(jxyo)=dta_ect(indxyo)
         ENDDO
      CASE(3)
         DO jxyo = 1,xyoend
            indxyo=obs_ord(jxyo)
            inddbs=obsnord(jxyo)
            xyo_nam(jxyo)=obs_nam(indxyo,inddbs)
            xyo_dim(jxyo)=obs_dim(indxyo,inddbs)
            xyo_ind(jxyo)=obs_ind(indxyo,inddbs)
            xyo_nbr(jxyo)=obs_nbr(indxyo,inddbs)
            xyo_jpi(jxyo)=dta_jpi(indxyo)
            xyo_jpj(jxyo)=dta_jpj(indxyo)
            xyo_jpk(jxyo)=dta_jpk(indxyo)
            xyo_jpt(jxyo)=dta_jpt(indxyo)
            xyo_ect(jxyo)=obs_ect(indxyo,inddbs)
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C      print *,kflaganlxyo
      SELECT CASE (kflaganlxyo)
      CASE (1)
         indxyomskdeb=0
         tabxyoend=0
         DO jxyo = 1,xyoend
            tabxyoend=tabxyoend+xyo_jpk(jxyo)*xyo_jpt(jxyo)
         ENDDO
      CASE (2)
         indxyomskdeb=varend
         tabxyoend=0
         DO jxyo = 1,xyoend
            tabxyoend=tabxyoend+xyo_jpk(jxyo)*xyo_jpt(jxyo)
         ENDDO
      CASE (3)
         indxyomskdeb=varend
         tabxyoend=xyoend
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -1.- Initialisation :
C ---------------------
C
C --- allocation tabindxyo
      allocate ( tabindxyo(1:tabxyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabindxyo(:) = 0
C --- allocation tabnbrxyo
      allocate ( tabnbrxyo(1:tabxyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabnbrxyo(:) = 0
C --- allocation tabnbrxyo1
      allocate ( tabnbrxyo1(1:tabxyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabnbrxyo1(:) = 0
C --- allocation tabjxyoxyo
      allocate ( tabjxyoxyo(1:tabxyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabjxyoxyo(:) = 0
C --- allocation tabjkxyo
      allocate ( tabjkxyo(1:tabxyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabjkxyo(:) = 0
C --- allocation RMS_ARxyo
      allocate ( RMS_ARxyo(1:tabxyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      RMS_ARxyo(:) = FREAL(0.0)
C --- allocation rap_AR_Rxyo
      allocate ( rap_AR_Rxyo(1:tabxyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      rap_AR_Rxyo(:) = FREAL(0.0)
      IF (klargdiffxyoorg) THEN
C --- allocation rap_AR_LRxyo
         allocate ( rap_AR_LRxyo(1:tabxyoend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         rap_AR_LRxyo(:) = FREAL(0.0)
      ENDIF
C --- allocation xyo_RMS_AR
      allocate ( xyo_RMS_AR(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_RMS_AR(:) = FREAL(0.0)
C --- allocation xyo_rap_AR_R
      allocate ( xyo_rap_AR_R(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_rap_AR_R(:) = FREAL(0.0)
      IF (klargdiffxyoorg) THEN
C --- allocation xyo_rap_AR_LR
         allocate ( xyo_rap_AR_LR(1:xyoend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         xyo_rap_AR_LR(:) = FREAL(0.0)
      ENDIF
C
C      print *,kflaganlxyo
      SELECT CASE (kflaganlxyo)
      CASE(1,2)
         jtabxyo=1
         js=1
         DO jxyo = 1,xyoend
            indxyomsk=indxyomskdeb+jxyo-1
            DO jt=1,xyo_jpt(jxyo)
            DO jk=1,xyo_jpk(jxyo)
               tabindxyo(jtabxyo)=js
               tabjxyoxyo(jtabxyo)=jxyo
               tabjkxyo(jtabxyo)=jk
               DO jj=1,xyo_jpj(jxyo)
               DO ji=1,xyo_jpi(jxyo)
                  js=js+ABS(IBITS(mask(ji,jj,jk,jt),indxyomsk,1))
               ENDDO
               ENDDO
               tabnbrxyo(jtabxyo)=js-tabindxyo(jtabxyo)
               jtabxyo=jtabxyo+1
            ENDDO
            ENDDO
         ENDDO
         IF ((kflaganlxyo.EQ.1).AND.(nallmem.LT.3)) THEN
            IF (jpxend.NE.(js-1)) THEN
               print *,'bug',jpxend,js-1
               GOTO 1000
            ENDIF
         ELSE
            IF (jpssize.NE.(js-1)) THEN
               print *,'bug',jpssize,js-1
               GOTO 1000
            ENDIF
         ENDIF
      CASE(3)
         tabindxyo(:) = xyo_ind(:)
         tabnbrxyo(1:tabxyoend) = xyo_nbr(1:tabxyoend)
         DO jxyo=1,tabxyoend
            tabjxyoxyo(jxyo) = jxyo
         ENDDO
         tabjkxyo(:) = 1
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -2.- Run the DIFF analysis :
C ----------------------------
C
      DO jnxyo=1,limjpnxyo(kflaganlxyo)
         IF (kflaganlxyo.EQ.1) THEN
            jindsbeg=arraynx_jindxbeg(jnxyo)
            jindsend=arraynx_jindxbeg(jnxyo)-1+arraynx_jpindxend(jnxyo)
            jsend=arraynx_jpindxend(jnxyo)
         ELSE
            jindsbeg=1
            jindsend=jpssize
            jsend=jpssize
         ENDIF
C --- reading object diffxyoref
         lectinfo=.TRUE.
         CALL readxyo(kargdiffxyoref,vectsref(:),
     $        jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
         IF (validextvar(kargdiffxyoref)) THEN
            spval=spvalvar
         ELSEIF (validextdta(kargdiffxyoref)) THEN
            spval=spvaldta
         ELSEIF (validextobs(kargdiffxyoref)) THEN
            spval=spvalobs
         ENDIF
         IF (ANY(vectsref(:).EQ.spval)) GOTO 101
         lectinfo=.FALSE.
         IF (largbias) THEN
            print *,'activate bias : diffref + bias'
            CALL readxyo(argbias,vects(:),
     $           jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:)) 
            IF (validextvar(argbias)) THEN
               spval=spvalvar
            ELSEIF (validextdta(argbias)) THEN
               spval=spvaldta
            ELSEIF (validextobs(argbias)) THEN
               spval=spvalobs
            ENDIF
            IF (ANY(vects(:).EQ.spval)) GOTO 102
            vectsref(:) = vectsref(:) + vects(:)
         ENDIF
C --- reading object inxyo
         CALL readxyo(karginxyo,vects(:),
     $     jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
         IF (validextvar(karginxyo)) THEN
            spval=spvalvar
         ELSEIF (validextdta(karginxyo)) THEN
            spval=spvaldta
         ELSEIF (validextobs(karginxyo)) THEN
            spval=spvalobs
         ENDIF
         IF (ANY(vects(:).EQ.spval)) GOTO 103
C --- reading object diffxyoorg
         IF (klargdiffxyoorg) THEN
            CALL readxyo(kargdiffxyoorg,vectsorg(:),
     $           jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
            IF (validextvar(kargdiffxyoorg)) THEN
               spval=spvalvar
            ELSEIF (validextdta(kargdiffxyoorg)) THEN
               spval=spvaldta
            ELSEIF (validextobs(kargdiffxyoorg)) THEN
               spval=spvalobs
            ENDIF
            IF (ANY(vectsorg(:).EQ.spval)) GOTO 104
         ENDIF
C
C --- allocation vcts_sr
         allocate ( vcts_sr(1:jsend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vcts_sr(:) = FREAL(0.0)
         IF (klargdiffxyoorg) THEN
C --- allocation vctso_sr
            allocate ( vctso_sr(1:jsend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            vctso_sr(:) = FREAL(0.0)
         ENDIF
C
         vcts_sr(1:jsend) = vects(1:jsend)-vectsref(1:jsend)
         IF (klargdiffxyoorg) vctso_sr(1:jsend) = 
     $        vectsorg(1:jsend)-vectsref(1:jsend)
C
         DO jtabxyo=1,tabxyoend
            jsdeb = MAX(1,tabindxyo(jtabxyo)-(jindsbeg-1))
            jsfin = MIN(jsend,tabindxyo(jtabxyo)+tabnbrxyo(jtabxyo)-1-(jindsbeg-1))
            IF (jsfin.GE.jsdeb) THEN
               RMS_ARxyo(jtabxyo) = RMS_ARxyo(jtabxyo) +
     $              FREAL(DOT_PRODUCT(vcts_sr(jsdeb:jsfin),
     $              vcts_sr(jsdeb:jsfin)))
               rap_AR_Rxyo(jtabxyo) = rap_AR_Rxyo(jtabxyo) +
     $              FREAL(DOT_PRODUCT(vectsref(jsdeb:jsfin),
     $              vectsref(jsdeb:jsfin)))
               IF (klargdiffxyoorg) rap_AR_LRxyo(jtabxyo) = 
     $              rap_AR_LRxyo(jtabxyo) +
     $              FREAL(DOT_PRODUCT(vctso_sr(jsdeb:jsfin),
     $              vctso_sr(jsdeb:jsfin)))
               tabnbrxyo1(jtabxyo) = tabnbrxyo1(jtabxyo) +
     $              jsfin-jsdeb+1
            ENDIF
         ENDDO
C
         IF (allocated(vcts_sr)) deallocate(vcts_sr)
         IF (allocated(vctso_sr)) deallocate(vctso_sr)
C
      ENDDO
C
C      print *,'beforeany',kflaganlxyo
      IF (ANY(tabnbrxyo1(:).NE.tabnbrxyo(:))) THEN
         print *,'bugANY'
         GOTO 1000
      ENDIF
C
      DO jtabxyo=1,tabxyoend
         IF (tabnbrxyo(jtabxyo).GT.0) THEN
            xyo_RMS_AR(tabjxyoxyo(jtabxyo)) =
     $           xyo_RMS_AR(tabjxyoxyo(jtabxyo)) +
     $           RMS_ARxyo(jtabxyo)
            xyo_rap_AR_R(tabjxyoxyo(jtabxyo)) =
     $           xyo_rap_AR_R(tabjxyoxyo(jtabxyo)) +
     $           rap_AR_Rxyo(jtabxyo)
            IF (klargdiffxyoorg) THEN
               xyo_rap_AR_LR(tabjxyoxyo(jtabxyo)) =
     $              xyo_rap_AR_LR(tabjxyoxyo(jtabxyo)) +
     $              rap_AR_LRxyo(jtabxyo)
            ENDIF
            xyo_nbr1(tabjxyoxyo(jtabxyo)) = 
     $           xyo_nbr1(tabjxyoxyo(jtabxyo)) +
     $           tabnbrxyo(jtabxyo)
            IF (klargdiffxyoorg) THEN
               IF (rap_AR_LRxyo(jtabxyo).NE.FREAL(0.0)) THEN
                  rap_AR_LRxyo(jtabxyo) = RMS_ARxyo(jtabxyo) /
     $                 rap_AR_LRxyo(jtabxyo)
                  rap_AR_LRxyo(jtabxyo) =
     $                 FREAL(SQRT(rap_AR_LRxyo(jtabxyo))*100.0)
               ENDIF
            ENDIF
            IF (rap_AR_Rxyo(jtabxyo).NE.FREAL(0.0)) THEN
               rap_AR_Rxyo(jtabxyo) = RMS_ARxyo(jtabxyo) /
     $              rap_AR_Rxyo(jtabxyo)
               rap_AR_Rxyo(jtabxyo) = 
     $              FREAL(SQRT(rap_AR_Rxyo(jtabxyo))*100.0)
            ENDIF
            RMS_ARxyo(jtabxyo) = FREAL(SQRT(RMS_ARxyo(jtabxyo) /
     $           tabnbrxyo(jtabxyo)))
            IF ((lmoyect)
     $           .AND.(xyo_ect(tabjxyoxyo(jtabxyo)).NE.FREAL(1.0))) THEN
               RMS_ARxyo(jtabxyo)  =RMS_ARxyo(jtabxyo) *
     $              xyo_ect(tabjxyoxyo(jtabxyo))
            ENDIF
         ELSE
            RMS_ARxyo(jtabxyo)   = FREAL(0.0)
            rap_AR_Rxyo(jtabxyo) = spval1
            IF (klargdiffxyoorg) rap_AR_LRxyo(jtabxyo)= spval1
         ENDIF
      ENDDO
C
C      print *,'beforeany',kflaganlxyo
      IF (ANY(xyo_nbr1(:).NE.xyo_nbr(:))) THEN
         print *,'bugANY2'
         GOTO 1000
      ENDIF
C
      DO jxyo =1,xyoend
         IF (klargdiffxyoorg) THEN
            IF (xyo_rap_AR_LR(jxyo).NE.FREAL(0.0)) THEN
               xyo_rap_AR_LR(jxyo) = xyo_RMS_AR(jxyo) / 
     $              xyo_rap_AR_LR(jxyo)
               xyo_rap_AR_LR(jxyo) = 
     $              FREAL(SQRT(xyo_rap_AR_LR(jxyo))*100.0)
            ELSE
               xyo_rap_AR_LR(jxyo) = spval1
            ENDIF
         ENDIF
         IF (xyo_rap_AR_R(jxyo).NE.FREAL(0.0)) THEN
            xyo_rap_AR_R(jxyo) = xyo_RMS_AR(jxyo) /
     $           xyo_rap_AR_R(jxyo)
            xyo_rap_AR_R(jxyo) = 
     $           FREAL(SQRT(xyo_rap_AR_R(jxyo))*100.0)
         ELSE
            xyo_rap_AR_R(jxyo) = spval1
         ENDIF
         xyo_RMS_AR(jxyo) = 
     $        FREAL(SQRT(xyo_RMS_AR(jxyo)/xyo_nbr(jxyo)))
         IF ((lmoyect).AND.(xyo_ect(jxyo).NE.FREAL(1.0))) THEN
            xyo_RMS_AR(jxyo) = xyo_RMS_AR(jxyo) *
     $           xyo_ect(jxyo)
         ENDIF
      ENDDO
C
C -3.- Print-Control :
C --------------------
C
C
      IF (klargdiffxyoorg) THEN
         WRITE(ligne,'(A,A)')
     $        ' ------------------------------------',
     $        '---------------------------------------'
         ASSIGN 12 to labtitr
      ELSE
         WRITE(ligne,'(A,A)')
     $        ' ------------------------------------',
     $        '-------------------------'
         ASSIGN 11 to labtitr
      ENDIF
      ASSIGN 20 to labnbre
C
      PRINT *,'  -------------------  Array of ERRORS'
     $       ,'  -------------------'
      PRINT '(A)', ligne(1:lenv(ligne))
      IF (klargdiffxyoorg) THEN
         PRINT labtitr, 'Var','Method ','Nbr pts'
     $        ,' RMS[AR] ',' AR_R(%) ',' AR_LR(%)'
      ELSE
         PRINT labtitr, 'Var','Method ','Nbr pts'
     $        ,' RMS[AR] ',' AR_R(%) '
      ENDIF
      PRINT '(A)', ligne(1:lenv(ligne))
      IF (nprint.GE.0) THEN
         WRITE(numout,'(A)') '   ---------  Array of ERRORS  ---------'
         WRITE(numout,'(A)') ligne(1:lenv(ligne))
         IF (klargdiffxyoorg) THEN
            WRITE(numout,labtitr) 'Var','Method','Nbr pts'
     $        ,' RMS[AR] ',' AR_R(%) ',' AR_LR(%)'
         ELSE
            WRITE(numout,labtitr) 'Var','Method','Nbr pts'
     $           ,' RMS[AR] ',' AR_R(%) '
         ENDIF
         WRITE(numout,'(A)') ligne(1:lenv(ligne))
      ENDIF
C
      jtabxyo=1
      DO jxyo=1,xyoend
         IF (xyo_dim(jxyo).GE.3) THEN
            WRITE(methode,'("zonk",I3.3,"-",I3.3)') 
     $           1,xyo_jpk(jxyo)
         ELSE
            WRITE(methode,'("zonk",I3.3,"-",I3.3)') 
     $           0,0
         ENDIF
         WRITE(lignbre,labnbre) xyo_nam(jxyo),methode,
     $        xyo_nbr(jxyo),xyo_RMS_AR(jxyo)
         IF (xyo_rap_AR_R(jxyo).NE.spval1) THEN
            WRITE(lignbre1,'(A,3X,F11.7)') lignbre(1:lenv(lignbre))
     $           ,xyo_rap_AR_R(jxyo)
            lignbre=lignbre1
         ELSE
            WRITE(lignbre1,'(A,3X,A11)') lignbre(1:lenv(lignbre))
     $           ,'- div / 0 -'
            lignbre=lignbre1
         ENDIF
         IF (klargdiffxyoorg) THEN
            IF (xyo_rap_AR_LR(jxyo).NE.spval1) THEN
               WRITE(lignbre1,'(A,3X,F11.7)') lignbre(1:lenv(lignbre))
     $              ,xyo_rap_AR_LR(jxyo)
               lignbre=lignbre1
            ELSE
               WRITE(lignbre1,'(A,3X,A11)') lignbre(1:lenv(lignbre))
     $              ,'- div / 0 -'
               lignbre=lignbre1
            ENDIF
         ENDIF
         IF (nprint.GE.0) WRITE(numout,'(A)') lignbre(1:lenv(lignbre))
         PRINT '(A)',lignbre(1:lenv(lignbre))
         IF (  ((kflaganlxyo.NE.3).AND.(xyo_dim(jxyo).GE.3)) ) THEN
            jtabxyo2=MIN(tabxyoend,(jtabxyo+xyo_jpk(jxyo)*xyo_jpt(jxyo)-1))
            DO jtabxyo1=jtabxyo,jtabxyo2
               WRITE(methode,'("zonk",I3.3,"-",I3.3)') 
     $                 (jtabxyo1-jtabxyo+1),(jtabxyo1-jtabxyo+1)
               WRITE(lignbre,labnbre) xyo_nam(jxyo),methode,
     $              tabnbrxyo(jtabxyo1),RMS_ARxyo(jtabxyo1)
               IF (rap_AR_Rxyo(jtabxyo1).NE.spval1) THEN
                  WRITE(lignbre1,'(A,3X,F11.7)') lignbre(1:lenv(lignbre))
     $                 ,rap_AR_Rxyo(jtabxyo1)
                  lignbre=lignbre1
               ELSE
                  WRITE(lignbre1,'(A,3X,A11)') lignbre(1:lenv(lignbre))
     $                 ,'- div / 0 -'
                  lignbre=lignbre1
               ENDIF
               IF (klargdiffxyoorg) THEN
                  IF (rap_AR_LRxyo(jtabxyo1).NE.spval1) THEN
                     WRITE(lignbre1,'(A,3X,F11.7)') lignbre(1:lenv(lignbre))
     $                    ,rap_AR_LRxyo(jtabxyo1)
                     lignbre=lignbre1
                  ELSE
                     WRITE(lignbre1,'(A,3X,A11)') lignbre(1:lenv(lignbre))
     $                    ,'- div / 0 -'
                     lignbre=lignbre1
                 ENDIF
               ENDIF
               IF (nprint.GE.0) WRITE(numout,'(A)') lignbre(1:lenv(lignbre))
               IF (MOD(jtabxyo1-jtabxyo,zonkfreq).EQ.0) 
     $              PRINT '(A)',lignbre(1:lenv(lignbre))
            ENDDO
         ENDIF
         jtabxyo=jtabxyo+xyo_jpk(jxyo)*xyo_jpt(jxyo)
         PRINT '(A)', ligne(1:lenv(ligne))
         WRITE(numout,'(A)') ligne(1:lenv(ligne))
      ENDDO
C
C --- deallocate
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(vectsref)) deallocate(vectsref)
      IF (allocated(vectsorg)) deallocate(vectsorg)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(xyo_dim)) deallocate(xyo_dim)
      IF (allocated(xyo_ind)) deallocate(xyo_ind)
      IF (allocated(xyo_nbr)) deallocate(xyo_nbr)
      IF (allocated(xyo_nbr1)) deallocate(xyo_nbr1)
      IF (allocated(xyo_jpi)) deallocate(xyo_jpi)
      IF (allocated(xyo_jpj)) deallocate(xyo_jpj)
      IF (allocated(xyo_jpk)) deallocate(xyo_jpk)
      IF (allocated(xyo_jpt)) deallocate(xyo_jpt)
      IF (allocated(xyo_ect)) deallocate(xyo_ect)
      IF (allocated(xyo_nam)) deallocate(xyo_nam)
      IF (allocated(tabindxyo)) deallocate(tabindxyo)
      IF (allocated(tabnbrxyo)) deallocate(tabnbrxyo)
      IF (allocated(tabnbrxyo1)) deallocate(tabnbrxyo1)
      IF (allocated(tabjxyoxyo)) deallocate(tabjxyoxyo)
      IF (allocated(tabjkxyo)) deallocate(tabjkxyo)
      IF (allocated(RMS_ARxyo)) deallocate(RMS_ARxyo)
      IF (allocated(rap_AR_Rxyo)) deallocate(rap_AR_Rxyo)
      IF (allocated(rap_AR_LRxyo)) deallocate(rap_AR_LRxyo)
C
      RETURN
C
C --- format definitions
C
C --- titles
 10   FORMAT(3X,A3,3X,"|",3X,A7,3X,"|",1X,A7,1(1X,"|",2X,A9,1X))
 11   FORMAT(3X,A3,3X,"|",3X,A7,3X,"|",1X,A7,2(1X,"|",2X,A9,1X))
 12   FORMAT(3X,A3,3X,"|",3X,A7,3X,"|",1X,A7,3(1X,"|",2X,A9,1X))
C --- values
 20   FORMAT(2X,A5,2X,"|",1X,A11,3X,I7.7,3X,1(1PE11.5E2,3X))
 21   FORMAT(2X,A5,2X,"|",1X,A11,3X,I7.7,3X,2(1PE11.5E2,3X))
 22   FORMAT(2X,A5,2X,"|",1X,A11,3X,I7.7,3X,3(1PE11.5E2,3X))
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algodiff','algodiffbylev')
 1001 CALL printerror2(0,1001,3,'algodiff','algodiffbylev')
C
 101  WRITE (texterror,*) 'spval exist in the file ',
     $     kargdiffxyoref(1:lenv(kargdiffxyoref)),
     $     ' => interf. it in .obs'
      CALL printerror2(0,101,3,'algodiff','algodiffbylev',
     $     comment=texterror)
 102  WRITE (texterror,*) 'spval exist in the file ',
     $     argbias(1:lenv(argbias)),
     $     ' => interf. it in .obs'
      CALL printerror2(0,102,3,'algodiff','algodiffbylev',
     $     comment=texterror)
 103  WRITE (texterror,*) 'spval exist in the file ',
     $     karginxyo(1:lenv(karginxyo)),
     $     ' => interf. it in .obs'
      CALL printerror2(0,103,3,'algodiff','algodiffbylev',
     $     comment=texterror)
 104  WRITE (texterror,*) 'spval exist in the file ',
     $     kargdiffxyoorg(1:lenv(kargdiffxyoorg)),
     $     ' => interf. it in .obs'
      CALL printerror2(0,104,3,'algodiff','algodiffbylev',
     $     comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algomeanstdbyzon(karginxyo,kflaganlxyo,knam_grouparea,
     $        ktab_nbarea,ktab_grouparea)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Algorithm for difference calcul on state :
CCC  -------
CCC  Method :
CCC  ------
CCC  Input :
CCC  -----
CCC  Output :
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : spvalvar,spvaldta,spvalobs,
     $     jpx,jpxend,jpo,jpitpend,jpy,poscoefobs,
     $     arraynx_jindxbeg,arraynx_jpindxend
      use hioxyo
      use utilvalid
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      INTEGER, intent(in) :: kflaganlxyo
      CHARACTER(len=*), intent(in) :: karginxyo
      CHARACTER(len=*), dimension(:), intent(in) :: knam_grouparea
      INTEGER, dimension(:), intent(in) :: ktab_nbarea
      INTEGER, dimension(:,:), intent(in) :: ktab_grouparea
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:), allocatable, save :: vectszon
      BIGREAL, dimension(:), allocatable, save :: vectyzon
C
      INTEGER :: allocok,jpssize,jpysize,jpitpsize
      INTEGER :: jpgroupsize,jpareasize,jgroup,jarea
      LOGICAL :: lectinfo,lmoyectold
      INTEGER :: jnxyo,flagxyo1,flagcfg
      INTEGER :: jsend,jindsbeg,jindsend,kjzondeb,kjzonfin
      INTEGER :: xyoend,jxyo,indxyo,inddbs
      INTEGER :: js,jsdeb,jsfin,jind,kjzon
      INTEGER, dimension(1:1) :: jitp
      INTEGER :: tabxyoend,tabxyoend1,jtabxyo
      INTEGER, dimension(:), allocatable :: xyo_ind,xyo_nbr
      INTEGER, dimension(:,:), allocatable :: xyo_nbrgrp
      BIGREAL, dimension(:), allocatable :: xyo_ect
      BIGREAL, dimension(:,:), allocatable :: xyo_zonmean,
     $     xyo_zonstd
      CHARACTER(len=varlg), dimension(:), allocatable :: xyo_nam
      INTEGER, dimension (:,:), allocatable :: tabnbrxyo
      BIGREAL, dimension(:,:), allocatable :: zonmeanxyo,zonstdxyo
      BIGREAL :: spval1,spval
      INTEGER :: labtitr,labnbre
      CHARACTER(len=bgword) :: ligne,lignbre,lignbre1,
     $     methode,fhisto1,fhisto2
      CHARACTER(len=bgword) :: myformat
      INTEGER, parameter  :: zonkfreq = 1
CC----------------------------------------------------------------------
C
      jpysize=jpy
      jpgroupsize=size(ktab_grouparea,1)
      jpareasize=size(ktab_grouparea,2)
      SELECT CASE (kflaganlxyo)
      CASE (1)
         jpssize=jpx
         jpitpsize=1
         xyoend  = varend
      CASE (2)
         jpssize=jpy
         jpitpsize=1
         xyoend  = dtaend
      CASE (3)
         jpssize=jpo
         jpitpsize=jpitpend
         xyoend  = obsend
      CASE DEFAULT
         GOTO 1000
      END SELECT
C --- allocation xyo_nbr
      allocate ( xyo_nbr(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_nbr(:) = 0
C --- allocation xyo_nbrgrp
      allocate ( xyo_nbrgrp(1:xyoend,1:(jpgroupsize+1)), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_nbrgrp(:,:) = 0
C --- allocation xyo_ind
      allocate ( xyo_ind(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_ind(:) = 0
C --- allocation xyo_ect
      allocate ( xyo_ect(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_ect(:) = FREAL(1.0)
C --- allocation xyo_nam
      allocate ( xyo_nam(1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_nam(:) = ' '
CCC---------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine algomeanstdbyzon &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*)
      ENDIF
C
C -0.- Initialisation :
C ---------------------
C
      IF ((ktab_nbarea(1).GT.0).AND.
     $     (ktab_nbarea(1).LE.jpareasize)) THEN
         kjzondeb = ktab_grouparea(1,1)
         kjzonfin = ktab_grouparea(1,1)
         DO kjzon=2,ktab_nbarea(1)
            kjzondeb = MIN(kjzondeb,ktab_grouparea(1,kjzon))
            kjzonfin = MAX(kjzonfin,ktab_grouparea(1,kjzon))
         ENDDO
      ELSE
         kjzondeb=1
         kjzonfin=1
      ENDIF
      DO jgroup=2,jpgroupsize
         IF ((ktab_nbarea(jgroup).GT.0).AND.
     $        (ktab_nbarea(jgroup).LE.jpareasize)) THEN
         DO kjzon=1,ktab_nbarea(jgroup)
            kjzondeb = MIN(kjzondeb,ktab_grouparea(jgroup,kjzon))
            kjzonfin = MAX(kjzonfin,ktab_grouparea(jgroup,kjzon))
         ENDDO
         ENDIF
      ENDDO
C BUG UQBAR vectorization
C      IF ((ktab_nbarea(1).GT.0).AND.
C     $     (ktab_nbarea(1).LE.jpareasize)) THEN
C         kjzondeb = MINVAL(ktab_grouparea(1,1:ktab_nbarea(1)))
C         kjzonfin = MAXVAL(ktab_grouparea(1,1:ktab_nbarea(1)))
C      ELSE
C         kjzondeb=1
C         kjzonfin=1
C      ENDIF
C      DO jgroup=2,jpgroupsize
C            IF ((ktab_nbarea(jgroup).GT.0).AND.
C     $        (ktab_nbarea(jgroup).LE.jpareasize)) THEN
C            kjzondeb = MIN(kjzondeb,
C     $           MINVAL(ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))))
C            kjzonfin  =MAX(kjzonfin,
C     $           MAXVAL(ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))))
C         ENDIF
C      ENDDO
      spval1=spvalvar
C --- allocation poscoefobs
      allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      SELECT CASE (kflaganlxyo)
      CASE (1)
         DO js=1,jpssize
            poscoefobs(js,:) = type_poscoef(js,FREAL(1.0))
         ENDDO
      CASE (2)
         DO js=1,jpssize
            poscoefobs(js,:) = type_poscoef(js,FREAL(1.0))
         ENDDO
      CASE (3)
         poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
C --- reading config.obs
         flagcfg=3
         CALL readcfgobs (argconfigobs,flagcfg,
     $        kposcoefobs=poscoefobs(:,:))
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C --- allocation vects
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
C
C --- allocation vectszon
      allocate ( vectszon(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectszon(:) = FREAL(0.0)
C
      SELECT CASE (kflaganlxyo)
      CASE (1,2)
C --- nothing
      CASE (3)
C --- allocation vectyzon
         allocate ( vectyzon(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectyzon(:) = FREAL(0.0)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      SELECT CASE (kflaganlxyo)
      CASE(1)
         DO jxyo = 1,xyoend
            indxyo= var_ord(jxyo)
            xyo_nam(jxyo)=var_nam(indxyo)
            xyo_ind(jxyo)=var_ind(indxyo)
            xyo_nbr(jxyo)=var_nbr(indxyo)
            xyo_ect(jxyo)=var_ect(indxyo)
         ENDDO
      CASE(2)
         DO jxyo = 1,xyoend
            indxyo=dta_ord(jxyo)
            xyo_nam(jxyo)=dta_nam(indxyo)
            xyo_ind(jxyo)=dta_ind(indxyo)
            xyo_nbr(jxyo)=dta_nbr(indxyo)
            xyo_ect(jxyo)=dta_ect(indxyo)
         ENDDO
      CASE(3)
         DO jxyo = 1,xyoend
            indxyo=obs_ord(jxyo)
            inddbs=obsnord(jxyo)
            xyo_nam(jxyo)=obs_nam(indxyo,inddbs)
            xyo_ind(jxyo)=obs_ind(indxyo,inddbs)
            xyo_nbr(jxyo)=obs_nbr(indxyo,inddbs)
            xyo_ect(jxyo)=obs_ect(indxyo,inddbs)
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      tabxyoend=(kjzonfin-kjzondeb+1)
      tabxyoend1=tabxyoend+1
C
C -1.- Initialisation :
C ---------------------
C
C --- allocation tabnbrxyo
      allocate ( tabnbrxyo(0:tabxyoend1,1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabnbrxyo(:,:) = 0
C --- allocation zonmeanxyo
      allocate ( zonmeanxyo(0:tabxyoend1,1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      zonmeanxyo(:,:) = FREAL(0.0)
C --- allocation zonstdxyo
      allocate ( zonstdxyo(0:tabxyoend1,1:xyoend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      zonstdxyo(:,:) = FREAL(0.0)
C --- allocation xyo_zonmean
      allocate ( xyo_zonmean(1:xyoend,1:(1+jpgroupsize)), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_zonmean(:,:) = FREAL(0.0)
C --- allocation xyo_zonstd
      allocate ( xyo_zonstd(1:xyoend,1:(1+jpgroupsize)), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      xyo_zonstd(:,:) = FREAL(0.0)
C
C -2.- Run the DIFF analysis :
C ----------------------------
C
      DO jnxyo=1,limjpnxyo(kflaganlxyo)
         IF (kflaganlxyo.EQ.1) THEN
            jindsbeg=arraynx_jindxbeg(jnxyo)
            jindsend=arraynx_jindxbeg(jnxyo)-1+arraynx_jpindxend(jnxyo)
            jsend=arraynx_jpindxend(jnxyo)
         ELSE
            jindsbeg=1
            jindsend=jpssize
            jsend=jpssize
         ENDIF
C --- reading object inxyo
         CALL readxyo(karginxyo,vects(:),
     $     jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:)) 
         IF (validextvar(karginxyo)) THEN
            spval=spvalvar
         ELSEIF (validextdta(karginxyo)) THEN
            spval=spvaldta
         ELSEIF (validextobs(karginxyo)) THEN
            spval=spvalobs
         ENDIF
         IF (ANY(vects(:).EQ.spval)) GOTO 101
C --- reading object inpartvar
         lmoyectold=lmoyect
         lmoyect=.FALSE.
         SELECT CASE (kflaganlxyo)
         CASE (1,2)
            CALL readxyo(arginpartvar,vectszon(:),
     $           jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:)) 
         CASE (3)
            flagxyo1=2
            CALL readxyo(arginpartvar,vectyzon(:),
     $           jnxyo,lectinfo,flagxyo1,poscoefobs(:,:)) 
C           vectszon(:jsend) = vectyzon(poscoefobs(:jsend,1)%pos)
            DO js=1,jsend
              jitp(1:1)=MAXLOC(poscoefobs(js,1:jpitpend)%coef)
              vectszon(js) = vectyzon(poscoefobs(js,jitp(1))%pos)
            ENDDO
         CASE DEFAULT
            GOTO 1000
         END SELECT
         lmoyect=lmoyectold
         vectszon(1:jsend)=FREAL(NINT(vectszon(1:jsend)))
         DO js=1,jsend
            vectszon(js)=vectszon(js)-FREAL(kjzondeb-1)
            vectszon(js)=MAX(FREAL(0.0),vectszon(js))
            vectszon(js)=MIN(vectszon(js),FREAL(kjzonfin-kjzondeb+1+1))
         ENDDO
#if defined _NEC
CCC!CDIR NOVECTOR
         IF (ALL(NINT(vectszon(1:jsend)).GT.tabxyoend1)) GOTO 1000
CCC         IF (INT(MAXVAL(vectszon(1:jsend))).GT.tabxyoend1) GOTO 1000
CCC!CDIR NOVECTOR
         IF (ALL(NINT(vectszon(1:jsend)).LT.0)) GOTO 1000
CCC         IF (INT(MINVAL(vectszon(1:jsend))).LT.0) GOTO 1000
#else
         IF (INT(MAXVAL(vectszon(1:jsend))).GT.tabxyoend1) GOTO 1000
         IF (INT(MINVAL(vectszon(1:jsend))).LT.0) GOTO 1000
#endif
C
         DO jxyo=1,xyoend
            jsdeb = MAX(1,xyo_ind(jxyo)-(jindsbeg-1))
            jsfin = MIN(jsend,xyo_nbr(jxyo)+xyo_ind(jxyo)-1-(jindsbeg-1))
            DO js=jsdeb,jsfin
C 1=>X
               zonmeanxyo(INT(vectszon(js)),jxyo) = 
     $              zonmeanxyo(INT(vectszon(js)),jxyo) +
     $              vects(js)
C 2=>X*X
               zonstdxyo(INT(vectszon(js)),jxyo) = 
     $              zonstdxyo(INT(vectszon(js)),jxyo) +
     $              vects(js)*vects(js)
            ENDDO
C
            DO js=jsdeb,jsfin
               tabnbrxyo(INT(vectszon(js)),jxyo) = 
     $              tabnbrxyo(INT(vectszon(js)),jxyo) + 1
            ENDDO
         ENDDO
C
      ENDDO
C
      DO jxyo=1,xyoend
      DO jgroup=1,jpgroupsize
      DO jarea=1,ktab_nbarea(jgroup)
         jtabxyo=ktab_grouparea(jgroup,jarea)-(kjzondeb-1)
         IF (tabnbrxyo(jtabxyo,jxyo).GT.0) THEN
C
            xyo_zonmean(jxyo,jgroup+1) =
     $           xyo_zonmean(jxyo,jgroup+1)+zonmeanxyo(jtabxyo,jxyo)
            xyo_zonstd(jxyo,jgroup+1) =
     $           xyo_zonstd(jxyo,jgroup+1)+zonstdxyo(jtabxyo,jxyo)
C
            xyo_nbrgrp(jxyo,jgroup+1) = 
     $           xyo_nbrgrp(jxyo,jgroup+1) + tabnbrxyo(jtabxyo,jxyo)
C
         ENDIF
      ENDDO
      ENDDO
      ENDDO
#if defined _NEC
!CDIR NOVECTOR
      DO jxyo=1,xyoend
!CDIR NOVECTOR
         DO jtabxyo=1,tabxyoend
#else
      DO jxyo=1,xyoend
         DO jtabxyo=1,tabxyoend
#endif
         IF (tabnbrxyo(jtabxyo,jxyo).GT.0) THEN
            xyo_zonmean(jxyo,1) =
     $           xyo_zonmean(jxyo,1)+zonmeanxyo(jtabxyo,jxyo)
            xyo_zonstd(jxyo,1) =
     $           xyo_zonstd(jxyo,1)+zonstdxyo(jtabxyo,jxyo)
            xyo_nbrgrp(jxyo,1) = 
     $           xyo_nbrgrp(jxyo,1) + tabnbrxyo(jtabxyo,jxyo)

            zonmeanxyo(jtabxyo,jxyo) = zonmeanxyo(jtabxyo,jxyo) /
     $           FREAL(tabnbrxyo(jtabxyo,jxyo))
            zonstdxyo(jtabxyo,jxyo) = SQRT(zonstdxyo(jtabxyo,jxyo)
C     $           * zonstdxyo(jtabxyo,jxyo)
     $           / FREAL(tabnbrxyo(jtabxyo,jxyo))
     $           - zonmeanxyo(jtabxyo,jxyo)
     $           * zonmeanxyo(jtabxyo,jxyo))
         ELSE
            zonmeanxyo(jtabxyo,jxyo)   = spval1
            zonstdxyo(jtabxyo,jxyo) = spval1
         ENDIF
      ENDDO
      ENDDO
C
      DO jxyo =1,xyoend
         xyo_zonmean(jxyo,1) = xyo_zonmean(jxyo,1) +
     $        zonmeanxyo(0,jxyo) + zonmeanxyo(tabxyoend1,jxyo)
         xyo_zonstd(jxyo,1) = xyo_zonstd(jxyo,1) +
     $        zonstdxyo(0,jxyo) + zonstdxyo(tabxyoend1,jxyo)
         xyo_nbrgrp(jxyo,1) = xyo_nbrgrp(jxyo,1) +
     $        tabnbrxyo(0,jxyo) + tabnbrxyo(tabxyoend1,jxyo)
      ENDDO
C
      DO jind =1,jpgroupsize+1
      DO jxyo =1,xyoend
         IF (xyo_nbrgrp(jxyo,jind).NE.FREAL(0.0)) THEN
            xyo_zonmean(jxyo,jind) = 
     $           xyo_zonmean(jxyo,jind)/FREAL(xyo_nbrgrp(jxyo,jind))
            xyo_zonstd(jxyo,jind) = SQRT(xyo_zonstd(jxyo,jind)
C     $           * xyo_zonstd(jxyo,jind)
     $           / FREAL(xyo_nbrgrp(jxyo,jind))
     $           - xyo_zonmean(jxyo,jind)
     $           * xyo_zonmean(jxyo,jind))
         ELSE
            xyo_zonmean(jxyo,jind) = spval1
            xyo_zonstd(jxyo,jind) = spval1
         ENDIF
      ENDDO
      ENDDO
C
C -3.- Print-Control :
C --------------------
C
      WRITE(ligne,'(A,A)') 
     $     ' ------------------------------------',
     $     '-------------------------'
      ASSIGN 11 to labtitr
      ASSIGN 20 to labnbre
C
      PRINT *,'  -------------------  Array of ERRORS'
     $       ,'  -------------------'
      PRINT '(A)', ligne(1:lenv(ligne))
      PRINT labtitr, 'Var','Method','Nbr pts'
     $        ,'  MEAN   ','   STD   '
      PRINT '(A)', ligne(1:lenv(ligne))
      IF (nprint.GE.0) THEN
         WRITE(numout,'(A)') '   --------  Array of VALUE   --------'
         WRITE(numout,'(A)') ligne(1:lenv(ligne))
         WRITE(numout,labtitr) 'Var','Methode','Nbr pts'
     $        ,'  MEAN   ','   STD   '
         WRITE(numout,'(A)') ligne(1:lenv(ligne))
      ENDIF
C
      DO jxyo=1,xyoend
C jind=1
         jind=1
         methode='-ALL AREAS-'
         WRITE(lignbre,labnbre) xyo_nam(jxyo),methode,
     $        xyo_nbrgrp(jxyo,jind),xyo_zonmean(jxyo,jind)
         IF (xyo_zonstd(jxyo,jind).NE.spval1) THEN
            WRITE(lignbre1,'(A,3X,F11.7)') lignbre(1:lenv(lignbre))
     $           ,xyo_zonstd(jxyo,jind)
            lignbre=lignbre1
         ELSE
            WRITE(lignbre1,'(A,3X,A11)') lignbre(1:lenv(lignbre))
     $           ,'- div / 0 -'
            lignbre=lignbre1
         ENDIF
         IF (nprint.GE.0) WRITE(numout,'(A)') lignbre(1:lenv(lignbre))
         PRINT '(A)',lignbre(1:lenv(lignbre))
C 
         DO jgroup=1,jpgroupsize
            methode=knam_grouparea(jgroup)(1:11)
            WRITE(lignbre,labnbre) xyo_nam(jxyo),methode,
     $           xyo_nbrgrp(jxyo,jgroup+1),xyo_zonmean(jxyo,jgroup+1)
            IF (xyo_zonstd(jxyo,jgroup+1).NE.spval1) THEN
               WRITE(lignbre1,'(A,3X,F11.7)') lignbre(1:lenv(lignbre))
     $              ,xyo_zonstd(jxyo,jgroup+1)
               lignbre=lignbre1
            ELSE
               WRITE(lignbre1,'(A,3X,A11)') lignbre(1:lenv(lignbre))
     $              ,'- div / 0 -'
               lignbre=lignbre1
            ENDIF
            SELECT CASE (ktab_nbarea(jgroup))
            CASE (7:)
               IF (nprint.GE.0) WRITE(numout,'(A)') lignbre(1:lenv(lignbre))
               PRINT '(A)',lignbre(1:lenv(lignbre))
            CASE (2:6)
               WRITE(myformat,'(A5,I1,A17)') '(A,A,',
     $              ktab_nbarea(jgroup)-1,'(I3.3,1X),I3.3,A)'
               IF (nprint.GE.0) WRITE(numout,myformat) lignbre(1:lenv(lignbre)),'   areas=(',
     $              (ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))),')'
               PRINT myformat,lignbre(1:lenv(lignbre)),'   areas=(',
     $              (ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))),')'
            CASE (1)
               WRITE(myformat,'(A)') '(A,A,I3.3,A)'
               IF (nprint.GE.0) WRITE(numout,myformat) lignbre(1:lenv(lignbre)),'   areas=(',
     $              (ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))),')'
               PRINT myformat,lignbre(1:lenv(lignbre)),'   areas=(',
     $              (ktab_grouparea(jgroup,1:ktab_nbarea(jgroup))),')'
            CASE DEFAULT
               GOTO 1000
            END SELECT
         ENDDO
C
         DO jtabxyo=1,tabxyoend
            WRITE(methode,'("jzon",I3.3,"-",I3.3)') 
     $           jtabxyo+(kjzondeb-1),
     $           jtabxyo+(kjzondeb-1)
            WRITE(lignbre,labnbre) xyo_nam(jxyo),methode,
     $           tabnbrxyo(jtabxyo,jxyo),zonmeanxyo(jtabxyo,jxyo)
            IF (zonstdxyo(jtabxyo,jxyo).NE.spval1) THEN
               WRITE(lignbre1,'(A,3X,F11.7)') lignbre(1:lenv(lignbre))
     $              ,zonstdxyo(jtabxyo,jxyo)
               lignbre=lignbre1
            ELSE
               WRITE(lignbre1,'(A,3X,A11)') lignbre(1:lenv(lignbre))
     $              ,'- div / 0 -'
               lignbre=lignbre1
            ENDIF
            IF (nprint.GE.0) WRITE(numout,'(A)') lignbre(1:lenv(lignbre))
         ENDDO
         PRINT '(A)', ligne(1:lenv(ligne))
         WRITE(numout,'(A)') ligne(1:lenv(ligne))
      ENDDO
C
C --- deallocate
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(vectszon)) deallocate(vectszon)
      IF (allocated(vectyzon)) deallocate(vectyzon)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(xyo_ind)) deallocate(xyo_ind)
      IF (allocated(xyo_nbr)) deallocate(xyo_nbr)
      IF (allocated(xyo_ect)) deallocate(xyo_ect)
      IF (allocated(xyo_nam)) deallocate(xyo_nam)
      IF (allocated(tabnbrxyo)) deallocate(tabnbrxyo)
      IF (allocated(zonmeanxyo)) deallocate(zonmeanxyo)
      IF (allocated(zonstdxyo)) deallocate(zonstdxyo)
C
      RETURN
C
C --- format definitions
C
C --- titles
 10   FORMAT(3X,A3,3X,"|",3X,A7,3X,"|",1X,A7,1(1X,"|",2X,A9,1X))
 11   FORMAT(3X,A3,3X,"|",3X,A7,3X,"|",1X,A7,2(1X,"|",2X,A9,1X))
 12   FORMAT(3X,A3,3X,"|",3X,A7,3X,"|",1X,A7,3(1X,"|",2X,A9,1X))
C --- values  
 20   FORMAT(2X,A5,2X,"|",1X,A11,3X,I7.7,3X,1(1PE11.4E2,3X))
 21   FORMAT(2X,A5,2X,"|",1X,A11,3X,I7.7,3X,2(1PE11.4E2,3X))
 22   FORMAT(2X,A5,2X,"|",1X,A11,3X,I7.7,3X,3(1PE11.4E2,3X))
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algodiff','algomeanstdbyzon')
 1001 CALL printerror2(0,1001,3,'algodiff','algomeanstdbyzon')
C
 101  WRITE (texterror,*) 'spval exist in the file ',
     $     karginxyo(1:lenv(karginxyo)),
     $     ' => interf. it in .obs'
      CALL printerror2(0,101,3,'algodiff','algomeanstdbyzon',
     $     comment=texterror)
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algodiff
