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
C ---                   MKPARTVAR.F                             ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2000-02 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C ---
C --- SUBROUTINE partvar
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkpartvar
      use mod_main
      use utilzone
      IMPLICIT NONE
      PRIVATE

      PUBLIC partvar

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE partvar(kargincfg,kargoutpartvar)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Generate a partition from a set of contours
CCC  -------        
CCC  Method :
CCC  ------
CCC  Input :           : no
CCC  -----
CCC  Output :          : no
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_coord
      use mod_cont
      use mod_spacexyo , only : jpx, jpxend
      use hioxyo
      use hiocnt
      use hiogrd
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kargincfg,kargoutpartvar
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vectx
C
      INTEGER :: allocok,jpxsize,cpt
      INTEGER :: ji,jj,jk,jt,jx,jvar,jvar1,jc,jp,jnxyo,jlay,typcnt
      INTEGER :: varecnt1, flagxy
      INTEGER :: indvar1, indvar, indvarmsk
      INTEGER, dimension(:), allocatable  :: varecnt
      CHARACTER(len=bgword) :: varfcnt1
      CHARACTER(len=bgword), dimension(:), allocatable  :: varfcnt
      CHARACTER(len=varlg) :: var_nam1
      LOGICAL :: lmoyectold
      TYPE (type_gridij) :: grdpt
      BIGREAL :: levpt
CC----------------------------------------------------------------------
C
      jpxsize=jpx
C --- allocation vectx
      allocate ( vectx(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectx(:) = FREAL(0.0)
C --- allocation varecnt
      allocate ( varecnt(1:nbvar), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      varecnt(:) = 0
C --- allocation varfcnt
      allocate ( varfcnt(1:nbvar), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      varfcnt(:) = ''
C ----------------------------------------------------------------------
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine mkpartvar &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Initialisation
C -------------------
C
      CALL openfile(10,kargincfg)
C
C --- read grid description
C
      DO jvar1 = 1,varend
         READ(10,*) var_nam1, varfcnt1, varecnt1
         DO jvar = 1,varend
            indvar = var_ord(jvar)
            IF (var_nam1(1:lenv(var_nam1))
     $          .EQ.var_nam(indvar)(1:lenv(var_nam(indvar))) ) THEN
               varfcnt(indvar) = varfcnt1
               varecnt(indvar) = varecnt1
            ENDIF
         ENDDO
      ENDDO
      CLOSE(10)
C
C -2.- Generate the partition :
C -----------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'MKPARTVAR : ',
     $         'Generating the partiton..'
      ENDIF
C
C --- Loop on variables
C
      cpt = 0
      jx = 0
      DO jvar = 1,varend
         indvar = var_ord(jvar)
         indvarmsk = jvar - 1
C
C --- allocate and read contour file
C
         CALL evalhdrcnt(varfcnt(indvar))
C
C --- allocation jpp
         allocate ( jpp(1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         jpp(:) = 0
C --- allocation contij
         allocate ( contij(1:jppend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contij(:,:)%longi = FREAL(0.0)
         contij(:,:)%latj = FREAL(0.0)
C --- allocation jplay
         allocate ( jplay(1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         jplay(:) = 0
C --- allocation contlevmin
         allocate ( contlevmin(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contlevmin(:,:) = FREAL(0.0)
C --- allocation contlevmax
         allocate ( contlevmax(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contlevmax(:,:) = FREAL(0.0)
C --- allocation contidx
         allocate ( contidx(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contidx(:,:) = 0
C
         typcnt = 0
         CALL readcnt(varfcnt(indvar),typcnt)
C
         IF (MOD(varecnt(indvar),10).EQ.0) THEN
C
C --- allocate and read grid file
C
            IF (varngrd(indvar).LE.2) THEN
C --- allocation longi
               allocate ( longi(1:var_jpi(indvar)) , stat = allocok )
               IF (allocok.GT.0) GOTO 1001
               longi(:) = FREAL(0.0)
C --- allocation latj
               allocate ( latj(1:var_jpj(indvar)) , stat = allocok )
               IF (allocok.GT.0) GOTO 1001
               latj(:) = FREAL(0.0)
            ELSE
C --- allocation gridij
               allocate ( gridij(1:var_jpi(indvar),1:var_jpj(indvar)) ,
     $              stat = allocok )
               IF (allocok.GT.0) GOTO 1001
               gridij(:,:) = type_gridij(FREAL(0.0),FREAL(0.0))
            ENDIF
C
            flagxy=1
            CALL readgrd(flagxy,indvar)
C
         ENDIF
C
C --- Loop on grid indices
C
         DO jt=1,var_jpt(indvar)
         DO jk=1,var_jpk(indvar)
            DO jj=1,var_jpj(indvar)
            DO ji=1,var_jpi(indvar)
               IF (IBITS(mask(ji,jj,jk,jt),indvarmsk,1).NE.0) THEN
                  jx = jx + 1
C
                  IF (MOD(varecnt(indvar),10).EQ.0) THEN
                     IF (varngrd(indvar).LE.2) THEN
                        grdpt%longi = longi(ji)
                        grdpt%latj  = latj(jj)
                     ELSE
                        grdpt%longi = gridij(ji,jj)%longi
                        grdpt%latj  = gridij(ji,jj)%latj
                     ENDIF
                  ELSE
                     grdpt%longi = FREAL(ji)
                     grdpt%latj  = FREAL(jj)
                  ENDIF
C
                  IF (varecnt(indvar)/10.EQ.0) THEN
                     levpt = var_lev(jk,indvar)
                  ELSE
                     levpt = FREAL(jk)
                  ENDIF
C
                  DO jc = 1,jpc
                     IF (inarea(grdpt,jc)) THEN
                        DO jlay = 1,jplay(jc)
                           IF ( ( levpt.GE.contlevmin(jlay,jc) ) .AND.
     $                        ( levpt.LE.contlevmax(jlay,jc) ) ) THEN
                              cpt = cpt + 1
                              vectx(jx) = contidx(jlay,jc)
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            ENDDO
         ENDDO
         ENDDO
C
C --- deallocate grid and contour arrays
C
         IF (allocated(longi)) deallocate(longi)
         IF (allocated(latj)) deallocate(latj)
         IF (allocated(gridij)) deallocate(gridij)
         IF (allocated(contij)) deallocate(contij)
         IF (allocated(contlevmin)) deallocate(contlevmin)
         IF (allocated(contlevmax)) deallocate(contlevmax)
         IF (allocated(jpp)) deallocate(jpp)
         IF (allocated(jplay)) deallocate(jplay)
         IF (allocated(contidx)) deallocate(contidx)
C
      ENDDO
C
      IF (jx.NE.jpxend) GOTO 1000
      print *, 'Number of points in the partition : ',
     $                 cpt,'/',jpx
C
C -3.- Writing the partition
C --------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'MKPARTVAR : ',
     $             'Writing the partition... '
      ENDIF
C
      lmoyectold=lmoyect
      lmoyect=.FALSE.
      jnxyo = 1
      CALL writevar (kargoutpartvar,vectx(:),jnxyo)
      lmoyect=lmoyectold
C
C --- deallocate
      IF (allocated(vectx)) deallocate(vectx)
      IF (allocated(varecnt)) deallocate(varecnt)
      IF (allocated(varfcnt)) deallocate(varfcnt)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'mkpartvar','partvar')
 1001 CALL printerror2(0,1001,3,'mkpartvar','partvar')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkpartvar
