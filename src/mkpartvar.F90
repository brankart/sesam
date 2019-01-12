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
! ---                   MKPARTVAR.F90                             ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2000-02 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! ---
! --- SUBROUTINE partvar
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkpartvar
      use mod_main
      use utilzone
      IMPLICIT NONE
      PRIVATE

      PUBLIC partvar

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE partvar(kargincfg,kargoutpartvar)
!---------------------------------------------------------------------
!
!  Purpose : Generate a partition from a set of contours
!  -------        
!  Method :
!  ------
!  Input :           : no
!  -----
!  Output :          : no
!  ------
!---------------------------------------------------------------------
! modules
! =======
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
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kargincfg,kargoutpartvar
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vectx
!
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
!----------------------------------------------------------------------
!
      jpxsize=jpx
! --- allocation vectx
      allocate ( vectx(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectx(:) = FREAL(0.0)
! --- allocation varecnt
      allocate ( varecnt(1:nbvar), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      varecnt(:) = 0
! --- allocation varfcnt
      allocate ( varfcnt(1:nbvar), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      varfcnt(:) = ''
! ----------------------------------------------------------------------
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine mkpartvar &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Initialisation
! -------------------
!
      CALL openfile(10,kargincfg)
!
! --- read grid description
!
      DO jvar1 = 1,varend
         READ(10,*) var_nam1, varfcnt1, varecnt1
         DO jvar = 1,varend
            indvar = var_ord(jvar)
            IF (var_nam1(1:lenv(var_nam1)) &
     &          .EQ.var_nam(indvar)(1:lenv(var_nam(indvar))) ) THEN
               varfcnt(indvar) = varfcnt1
               varecnt(indvar) = varecnt1
            ENDIF
         ENDDO
      ENDDO
      CLOSE(10)
!
! -2.- Generate the partition :
! -----------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'MKPARTVAR : ', &
     &         'Generating the partiton..'
      ENDIF
!
! --- Loop on variables
!
      cpt = 0
      jx = 0
      DO jvar = 1,varend
         indvar = var_ord(jvar)
         indvarmsk = jvar - 1
!
! --- allocate and read contour file
!
         CALL evalhdrcnt(varfcnt(indvar))
!
! --- allocation jpp
         allocate ( jpp(1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         jpp(:) = 0
! --- allocation contij
         allocate ( contij(1:jppend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contij(:,:)%longi = FREAL(0.0)
         contij(:,:)%latj = FREAL(0.0)
! --- allocation jplay
         allocate ( jplay(1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         jplay(:) = 0
! --- allocation contlevmin
         allocate ( contlevmin(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contlevmin(:,:) = FREAL(0.0)
! --- allocation contlevmax
         allocate ( contlevmax(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contlevmax(:,:) = FREAL(0.0)
! --- allocation contidx
         allocate ( contidx(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contidx(:,:) = 0
!
         typcnt = 0
         CALL readcnt(varfcnt(indvar),typcnt)
!
         IF (MOD(varecnt(indvar),10).EQ.0) THEN
!
! --- allocate and read grid file
!
            IF (varngrd(indvar).LE.2) THEN
! --- allocation longi
               allocate ( longi(1:var_jpi(indvar)) , stat = allocok )
               IF (allocok.GT.0) GOTO 1001
               longi(:) = FREAL(0.0)
! --- allocation latj
               allocate ( latj(1:var_jpj(indvar)) , stat = allocok )
               IF (allocok.GT.0) GOTO 1001
               latj(:) = FREAL(0.0)
            ELSE
! --- allocation gridij
               allocate ( gridij(1:var_jpi(indvar),1:var_jpj(indvar)) , &
     &              stat = allocok )
               IF (allocok.GT.0) GOTO 1001
               gridij(:,:) = type_gridij(FREAL(0.0),FREAL(0.0))
            ENDIF
!
            flagxy=1
            CALL readgrd(flagxy,indvar)
!
         ENDIF
!
! --- Loop on grid indices
!
         DO jt=1,var_jpt(indvar)
         DO jk=1,var_jpk(indvar)
            DO jj=1,var_jpj(indvar)
            DO ji=1,var_jpi(indvar)
               IF (IBITS(mask(ji,jj,jk,jt),indvarmsk,1).NE.0) THEN
                  jx = jx + 1
!
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
!
                  IF (varecnt(indvar)/10.EQ.0) THEN
                     levpt = var_lev(jk,indvar)
                  ELSE
                     levpt = FREAL(jk)
                  ENDIF
!
                  DO jc = 1,jpc
                     IF (inarea(grdpt,jc)) THEN
                        DO jlay = 1,jplay(jc)
                           IF ( ( levpt.GE.contlevmin(jlay,jc) ) .AND. &
     &                        ( levpt.LE.contlevmax(jlay,jc) ) ) THEN
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
!
! --- deallocate grid and contour arrays
!
         IF (allocated(longi)) deallocate(longi)
         IF (allocated(latj)) deallocate(latj)
         IF (allocated(gridij)) deallocate(gridij)
         IF (allocated(contij)) deallocate(contij)
         IF (allocated(contlevmin)) deallocate(contlevmin)
         IF (allocated(contlevmax)) deallocate(contlevmax)
         IF (allocated(jpp)) deallocate(jpp)
         IF (allocated(jplay)) deallocate(jplay)
         IF (allocated(contidx)) deallocate(contidx)
!
      ENDDO
!
      IF (jx.NE.jpxend) GOTO 1000
      print *, 'Number of points in the partition : ', &
     &                 cpt,'/',jpx
!
! -3.- Writing the partition
! --------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'MKPARTVAR : ', &
     &             'Writing the partition... '
      ENDIF
!
      lmoyectold=lmoyect
      lmoyect=.FALSE.
      jnxyo = 1
      CALL writevar (kargoutpartvar,vectx(:),jnxyo)
      lmoyect=lmoyectold
!
! --- deallocate
      IF (allocated(vectx)) deallocate(vectx)
      IF (allocated(varecnt)) deallocate(varecnt)
      IF (allocated(varfcnt)) deallocate(varfcnt)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mkpartvar','partvar')
 1001 CALL printerror2(0,1001,3,'mkpartvar','partvar')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkpartvar
