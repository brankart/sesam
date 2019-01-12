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
! ---                  HIOXYO.F90                                 ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12  (C.E. Testut)                       ---
! --- modification : 99-05  (C.E. Testut)                       ---
! --- modification : 01-06  (C.E. Testut)                       ---
! --- modification : 03-02  (J.M. Brankart)                     ---
! --- modification : 07-11  (J.M. Brankart)                     ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  readxyo    : Read Vx, Vy or Vo vector from input files
! --- SUBROUTINE  writeyo    : Write Vy or Vo vector to output files
! --- 
! --- SUBROUTINE  readfil    : Read Vx or Vy vector from 'var' or 'dta'
! ---                          input files (structured formats)
! --- SUBROUTINE  writenfil  : Write Vx or Vy vector in 'var' or 'dta'
! ---                          output files (structured formats)
! --- 
! --- SUBROUTINE  readvar    : Read Vx or Vy vector from 'var' file
! --- SUBROUTINE  writevar   : Write Vx vector in 'var' file
! --- 
! --- SUBROUTINE  readdta    : Read Vy vector from 'dta' file
! --- SUBROUTINE  writedta   : Write Vy vector in 'dta' file
! --- 
! --- SUBROUTINE  evalhdrobs : Read 'obs' file header
! --- SUBROUTINE  readvalobs : Read Vo vector from 'obs' file
! ---                          (observed values only)
! --- SUBROUTINE  readcfgobs : Read Vo vector from 'obs' file
! ---                          (observation configuration only)
! --- SUBROUTINE  writeobs   : Write Vo vector in 'obs' file
! --- SUBROUTINE  writesingleobs : Write segment of Vo vector in 'obs' file
! ---                          (corresponding to one observation database)
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE hioxyo
      use mod_main
      use liobimg
      use liodimg
      use liocdf
      use lionc
      use liocpak
      use lioobs
      use liocobs
      use liorst
      use utilvalid
      IMPLICIT NONE
      PRIVATE

      PUBLIC readxyo,writeyo,readfil,writenfil
      PUBLIC readvar,writevar,readdta,writedta
      PUBLIC evalhdrobs,readvalobs,readcfgobs
      PUBLIC writeobs,writesingleobs

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readxyo(fninxyo,kvects,kjns,klectinfo, &
     &     kflagxyo,kposcoefobs) 
!---------------------------------------------------------------------
!
!  Purpose : Read Vx, Vy or Vo vector from input files
!  -------
!  Method : Call appropriate routine to read the right vector object
!  ------
!  Input :     fninxyo   : filename
!  -----       kjns      : block index
!              klectinfo : read or not header of input files
!              kflagxyo  : vector type (1=Vx,2=Vy,3=Vo)
!              kposcoefobs : observation operator (for Vo objects only)
!  Output :    kvects    : 1D vector object (Vx, Vy or Vo)
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpyend
      use utilmkh
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*) :: fninxyo
      BIGREAL, dimension(:), intent(out) ::  kvects
      INTEGER, intent(in) :: kjns
      LOGICAL, intent(in) :: klectinfo
      INTEGER, intent(in) :: kflagxyo
      TYPE (type_poscoef), dimension(:,:), optional, intent(in) :: &
     &     kposcoefobs
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable :: vecty
!
      INTEGER :: allocok,jpssize
      INTEGER :: flagxyo,flagext
!----------------------------------------------------------------------
!
! Check validity of input filename
      jpssize = size(kvects,1)
      SELECT CASE (kflagxyo)
      CASE(1)
         IF (.NOT.(validextvar(fninxyo))) GOTO 101
      CASE(2)
         IF ((.NOT.(validextvar(fninxyo))) &
     &        .AND.(.NOT.(validextdta(fninxyo)))) GOTO 101
      CASE(3)
         IF (.NOT.(present(kposcoefobs))) GOTO 1000
         IF (jpssize.NE.size(kposcoefobs,1)) GOTO 1000 
         IF ((.NOT.(validextvar(fninxyo))) &
     &        .AND.(.NOT.(validextdta(fninxyo))) &
     &        .AND.(.NOT.(validextobs(fninxyo)))) GOTO 101
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Set file type (var/dta/obs)
      flagext=0
      IF (validextdta(fninxyo)) flagext=3
      IF (validextobs(fninxyo)) flagext=4
      IF (validextvar(fninxyo)) flagext=5
!
! Allocate temporary Vy vector (to read Vo from var or dta file)
      IF (kflagxyo.EQ.3) THEN
         IF ( (flagext.EQ.3) .OR. (flagext.EQ.5) ) THEN
            allocate ( vecty(1:jpyend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            vecty(:) = FREAL(0.0)
         ENDIF
      ENDIF
!
! Use appropriate routine to read var, dta or obs file
! ----------------------------------------------------
!
      SELECT CASE (kflagxyo)
      CASE(1,2)
         SELECT CASE (flagext)
         CASE(3)
! Read Vy vector from dta file
            CALL readdta(fninxyo,kvects,klectinfo)
         CASE(5)
! Read Vx or Vy vector from var file
            CALL readvar(fninxyo,kvects,kjns,klectinfo,kflagxyo) 
         CASE DEFAULT
            GOTO 1000
         END SELECT
      CASE(3)
         SELECT CASE (flagext)
         CASE(3)
! Read Vo vector from dta file
            CALL readdta(fninxyo,vecty,klectinfo)
            CALL mkhytoo(vecty(:),kvects(:),kposcoefobs(:,:))
         CASE(5)
! Read Vo vector from var file
            flagxyo=2
            CALL readvar(fninxyo,vecty,kjns,klectinfo,flagxyo)
            CALL mkhytoo(vecty(:),kvects(:),kposcoefobs(:,:))
         CASE(4)
! Read Vo vector from obs file
            CALL readvalobs(fninxyo,kvects(:))
         CASE DEFAULT
            GOTO 1000
         END SELECT
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (allocated(vecty)) deallocate(vecty)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'hioxyo','readxyo')
 1001 CALL printerror2(0,1001,3,'hioxyo','readxyo')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &                       fninxyo(1:lenv(fninxyo))
      CALL printerror2(0,101,3,'hioxyo','readxyo',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeyo(kfnoutyo,kvects,kflagxyo,kvectsrms, &
     &                   kgridijkobs,kposcoefobs) 
!---------------------------------------------------------------------
!
!  Purpose : Write Vy or Vo vector to output files
!  -------
!  Method : Call appropriate routine to write the right vector object
!  ------
!  Input :  kfnoutyo : output filename
!  -----    kvects : 1D vector object (Vy or Vo) to write
!           kflagxyo : vector type (2=Vy,3=Vo)
!           kvectsrms : observation error (for Vo vectors only) [obsolete]
!           kgridijkobs : observation location (for Vo vectors only)
!           kposcoefobs : observation operator (for Vo vectors only)
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutyo
      INTEGER, intent(in) :: kflagxyo
      BIGREAL, dimension(:), intent(in) :: kvects
      BIGREAL, dimension(:), optional, intent(in) :: kvectsrms
      TYPE (type_gridijk), dimension(:), optional,  intent(in)  ::  &
     &     kgridijkobs
      TYPE (type_poscoef), dimension(:,:), optional,  intent(in) ::  &
     &     kposcoefobs
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpssize
!----------------------------------------------------------------------
!
! Check validity of output filename
      jpssize = size(kvects,1)
      SELECT CASE (kflagxyo)
      CASE (2)
         IF (.NOT.(validextdta(kfnoutyo))) GOTO 101
      CASE (3)
         IF (.NOT.(validextobs(kfnoutyo))) GOTO 101
         IF (.NOT.(present(kvectsrms))) GOTO 1000
         IF (.NOT.(present(kgridijkobs))) GOTO 1000
         IF (.NOT.(present(kposcoefobs))) GOTO 1000
         IF (jpssize.NE.size(kvectsrms,1))  GOTO 1000
         IF (jpssize.NE.size(kgridijkobs,1))  GOTO 1000
         IF (jpssize.NE.size(kposcoefobs,1))  GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Use appropriate routine to write dta or obs file
! ------------------------------------------------
!
      SELECT CASE (kflagxyo)
      CASE (1)
         GOTO 1000
      CASE (2)
         CALL writedta(kfnoutyo,kvects) 
      CASE (3)
         CALL writeobs(kfnoutyo,kvects,kvectsrms,kgridijkobs, &
     &     kposcoefobs)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hioxyo','writeyo')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &                       kfnoutyo(1:lenv(kfnoutyo))
      CALL printerror2(0,101,3,'hixyo','writeyo',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readfil(kfninsxy,kvects,klectinfo,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Read Vx or Vy vector from 'var' or 'dta'
!  -------   input files (structured formats)
!
!  Method : Read all necessary files of any structured format
!  ------   (BIMG, DIMG or NETCDF)
!
!  Input  : kfninsxy : filename(s)
!  -----    klectinfo : read or not header of input files
!           kflagxyo : vector type (Vx/Vy/Vo)
!  Output : kvects : 1D vector object (Vx or Vy)
!  ------
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend,jpyend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninsxy
      BIGREAL, dimension(:), intent(out) :: kvects
      LOGICAL, intent(in) :: klectinfo
      INTEGER, intent(in) :: kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpsend,sxyend
      INTEGER, dimension(1:nbvar) :: sxy_ord,sxy_dim, &
     &     sxy_nbr,sxy_ind
      LOGICAL :: extsxyunit
      BIGREAL, dimension(1:nbvar) :: sxy_moy,sxy_ect
      INTEGER :: jsxy,indsxy,spos,somsxynbr
      INTEGER :: nsdeb, nsfin, nstot
      INTEGER :: jextsxy
      CHARACTER(len=bgword) :: fname, kform
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.3) THEN
         WRITE(numout,*) '*** ROUTINE : ./readfil'
         WRITE(numout,*) '    ==> READING file ', &
     &            kfninsxy(1:lenv(kfninsxy))
      ENDIF
!
! Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         jpsend  = jpxend
         sxyend  = varend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=var_ord(jsxy)
            indsxy= sxy_ord(jsxy)
            sxy_dim(indsxy)=var_dim(indsxy)
            sxy_nbr(indsxy)=var_nbr(indsxy)
            sxy_moy(indsxy)=var_moy(indsxy)
            sxy_ect(indsxy)=var_ect(indsxy)
            sxy_ind(indsxy)=var_ind(indsxy)
         ENDDO
      CASE(2)
! --- dta
         jpsend  = jpyend
         sxyend  = dtaend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=dta_ord(jsxy)
            indsxy= sxy_ord(jsxy)
            sxy_dim(indsxy)=dta_dim(indsxy)
            sxy_nbr(indsxy)=dta_nbr(indsxy)
            sxy_moy(indsxy)=dta_moy(indsxy)
            sxy_ect(indsxy)=dta_ect(indsxy)
            sxy_ind(indsxy)=dta_ind(indsxy)
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Get file extension
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         IF (validextvar(kfninsxy)) THEN
            jextsxy = indext(kfninsxy,extvartab,nbextvar)
            extsxyunit=extvarunit(jextsxy)
         ELSE
            GOTO 102
         ENDIF
      CASE(2)
! --- dta
         IF (validextdta(kfninsxy)) THEN
            jextsxy = indext(kfninsxy,extdtatab,nbextdta)
            extsxyunit=extdtaunit(jextsxy)
         ELSEIF (validextvar(kfninsxy)) THEN
            jextsxy=indext(kfninsxy,extvartab,nbextvar)
            extsxyunit=extvarunit(jextsxy)
         ELSE
            GOTO 102
         ENDIF
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Get position of variable name in file name
      IF (.NOT.(extsxyunit)) THEN
         spos = posit(kfninsxy,etoile)
         IF (spos.LT.1) GOTO 102
      ELSE
         spos = 1
      ENDIF
!
! -1.- Read all variable field in SESAM object
! --------------------------------------------
!
      nsdeb = 1
      nstot = 0
      nsfin = 0
      DO jsxy = 1,sxyend
         indsxy=sxy_ord(jsxy)
!
! Construct filename
         SELECT CASE(kflagxyo)
         CASE(1)
! --- var
            IF (extsxyunit) THEN
               WRITE (fname,'(a)') kfninsxy(1:lenv(kfninsxy))
            ELSE
               WRITE (fname,'(a,a,a)') kfninsxy(1:(spos-1)), &
     &              varinam(indsxy)(1:lenv(varinam(indsxy))), &
     &              kfninsxy((spos+1):lenv(kfninsxy))
            ENDIF
         CASE(2)
! --- dta
            IF (extsxyunit) THEN
               WRITE (fname,'(a)') kfninsxy(1:lenv(kfninsxy))
            ELSE
               IF (validextvar(kfninsxy)) WRITE (fname,'(a,a,a)')  &
     &              kfninsxy(1:(spos-1)), &
     &              varinam(indsxy)(1:lenv(varinam(indsxy))), &
     &              kfninsxy((spos+1):lenv(kfninsxy))
               IF (validextdta(kfninsxy)) WRITE (fname,'(a,a,a)')  &
     &              kfninsxy(1:(spos-1)), &
     &              dtainam(indsxy)(1:lenv(dtainam(indsxy))), &
     &              kfninsxy((spos+1):lenv(kfninsxy))
            ENDIF
         CASE DEFAULT
            GOTO 1000
         END SELECT
!
! Check if we are at the right place in 1D vector
         IF (nsdeb.NE.sxy_ind(indsxy)) GOTO 1000
         IF (sxy_nbr(indsxy).EQ.0) GOTO 1000
         nstot = sxy_nbr(indsxy)
         nsfin = nsdeb - 1 + nstot
         IF (nsfin.GT.size(kvects,1)) GOTO 1000
!
! Read variable field (using the right format)
         somsxynbr = 0
         SELECT CASE (jextsxy)
         CASE (1)
            CALL readbimg(fname,kvects(nsdeb:nsfin), &
     &           jsxy,somsxynbr,kflagxyo)
         CASE (2)
            CALL readdimg(fname,kvects(nsdeb:nsfin), &
     &           jsxy,somsxynbr,kflagxyo)
         CASE (3)
            CALL readcdf(fname,kvects(nsdeb:nsfin), &
     &           jsxy,somsxynbr,kflagxyo)
         CASE (4)
            CALL readnc(fname,kvects(nsdeb:nsfin), &
     &           jsxy,somsxynbr,kflagxyo)
         CASE DEFAULT
            GOTO 1000
         END SELECT
!
! Control print
         IF (nprint.GE.3) THEN
            SELECT CASE(kflagxyo)
            CASE(1)
! --- var
               WRITE(numout,*) '    ==> STORING variable ', &
     &              var_nam(indsxy)(1:lenv(var_nam(indsxy))), &
     &              ' in Vx vector'
            CASE(2)
! --- dta
               WRITE(numout,*) '    ==> STORING variable ', &
     &              dta_nam(indsxy)(1:lenv(dta_nam(indsxy))), &
     &              ' in Vy vector'
            CASE DEFAULT
               GOTO 1000
            END SELECT
            kform='(8x,a,i9)'
            WRITE(numout,kform) '- Initial index : ',nsdeb
            WRITE(numout,kform) '- Final index   : ', nsfin 
            WRITE(numout,kform) '- Segment size  : ', nsfin-nsdeb+1 
         ENDIF
         IF (nstot.NE.somsxynbr) GOTO 1000
         nsdeb = nsfin+1
!
      ENDDO
!
      jpsend = nsfin
!
! Check size of vector object
!
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         IF (jpsend.NE.jpxend) GOTO 104
      CASE(2)
! --- dta
         IF (jpsend.NE.jpyend) GOTO 105
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'hioxyo','readfil')
 1001 CALL printerror2(0,1001,3,'hioxyo','readfil')
!
 102  WRITE (texterror,*) 'bad input file name:', &
     &                 kfninsxy(1:lenv(kfninsxy))
      CALL printerror2(0,102,3,'hioxyo','readfil',comment=texterror)
 104  WRITE (texterror,*) 'incoherence between jpxend =',jpsend, &
     &     ' and jpxend from mask =',jpxend
      CALL printerror2(0,104,3,'hioxyo','readfil',comment=texterror)
 105  WRITE (texterror,*) 'incoherence between jpyend =',jpsend, &
     &     ' and jpyend from mask =',jpyend
      CALL printerror2(0,105,3,'hioxyo','readfil',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writenfil(kfnoutsxy,kvects,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Write Vx or Vy vector in 'var' or 'dta'
!  -------   output files (structured formats)
!
!  Method : Write all necessary files of any structured format
!  ------   (BIMG, DIMG or NETCDF)
!
!  Input :  kfnoutsxy: filename(s)
!  -----    kvects : 1D vector object (Vx or Vy)
!           kflagxyo : vector type (Vx,Vy,Vo)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutsxy
      BIGREAL, dimension(:), intent(in) :: kvects
      INTEGER, intent(in) :: kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: sxyend
      INTEGER, dimension(1:nbvar) :: sxy_ord,sxy_dim, &
     &     sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt,sxy_nbr,sxy_ind,sxyopos
      CHARACTER(len=varlg), dimension(1:nbvar) :: sxy_nam
      CHARACTER(len=bgword), dimension(1:nbvar) :: sxyonam
      LOGICAL, dimension(1:nbvar) :: sxy_towrite
      LOGICAL :: extsxyunit
      INTEGER :: jsxy,jsxy1,indsxy,indsxy1
      INTEGER :: spos
      CHARACTER(len=3) :: sxy
      INTEGER :: jtxy,indtxy,txyend
      INTEGER, dimension(1:nbvar) :: txy_indvct,txy_indmsk
      INTEGER :: jextsxy
      CHARACTER(len=bgword) :: fname, kform
      LOGICAL :: done
      CHARACTER(len=word80)  ::cdrec1,cdrec2,cdrec3,cdrec4
      BIGREAL :: idast
!---------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./writenfil'
         WRITE(numout,*) '    ==> WRITING file ', &
     &            kfnoutsxy(1:lenv(kfnoutsxy))
      ENDIF
!
! Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         sxy     = 'var'
         sxyend  = varend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=var_ord(jsxy)
            indsxy= sxy_ord(jsxy)
!
            sxy_nam(indsxy)=var_nam(indsxy)
            sxy_dim(indsxy)=var_dim(indsxy)
            sxy_jpi(indsxy)=var_jpi(indsxy)
            sxy_jpj(indsxy)=var_jpj(indsxy)
            sxy_jpk(indsxy)=var_jpk(indsxy)
            sxy_jpt(indsxy)=var_jpt(indsxy)
!
            sxy_nbr(indsxy)=var_nbr(indsxy)
            sxy_ind(indsxy)=var_ind(indsxy)
            sxyonam(indsxy)=varonam(indsxy)
            sxyopos(indsxy)=varopos(indsxy)
         ENDDO
      CASE(2)
! --- dta
         sxy     = 'dta'
         sxyend  = dtaend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=dta_ord(jsxy)
            indsxy= sxy_ord(jsxy)
!
            sxy_nam(indsxy)=dta_nam(indsxy)
            sxy_dim(indsxy)=dta_dim(indsxy)
            sxy_jpi(indsxy)=dta_jpi(indsxy)
            sxy_jpj(indsxy)=dta_jpj(indsxy)
            sxy_jpk(indsxy)=dta_jpk(indsxy)
            sxy_jpt(indsxy)=dta_jpt(indsxy)
!
            sxy_nbr(indsxy)=dta_nbr(indsxy)
            sxy_ind(indsxy)=dta_ind(indsxy)
            sxyonam(indsxy)=dtaonam(indsxy)
            sxyopos(indsxy)=dtaopos(indsxy)
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Get file extension
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         IF (.NOT.(validextvar(kfnoutsxy))) GOTO 102
         jextsxy = indext(kfnoutsxy,extvartab,nbextvar)
         extsxyunit=extvarunit(jextsxy)
      CASE(2)
! --- dta
         IF (.NOT.(validextdta(kfnoutsxy))) GOTO 102
         jextsxy = indext(kfnoutsxy,extdtatab,nbextdta)
         extsxyunit=extdtaunit(jextsxy)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Get position of variable name in file name
      IF (.NOT.(extsxyunit)) THEN
         spos = posit(kfnoutsxy,etoile)
         IF (spos.LT.1) GOTO 102
      ELSE
         spos = 1
      ENDIF
!
! Check coherence between array sxy_nbr and sxy_ind
      indsxy=sxy_ord(1)
      IF (sxy_ind(indsxy).NE.1) GOTO 104
      DO jsxy=2,sxyend
         indsxy1=sxy_ord(jsxy-1)
         indsxy=sxy_ord(jsxy)
         IF ((sxy_ind(indsxy)).NE.(sxy_ind(indsxy1)+sxy_nbr(indsxy1))) &
     &        GOTO 104
      ENDDO
!
! Initialize table of variable remaining to write
      DO jsxy=1,nbvar
         sxy_towrite(jsxy)=.FALSE.
      ENDDO
      DO jsxy=1,sxyend
         indsxy=sxy_ord(jsxy)
         sxy_towrite(indsxy)=.TRUE.         
      ENDDO
!
! Control print
      IF (nprint.GE.3) THEN
         DO jsxy=1,sxyend
            indsxy=sxy_ord(jsxy)
!
            SELECT CASE(kflagxyo)
            CASE(1)
! --- var
               WRITE(numout,*) '    ==> EXTRACTING variable ', &
     &              sxy_nam(indsxy)(1:lenv(sxy_nam(indsxy))), &
     &              ' from Vx vector'
            CASE(2)
! --- dta
               WRITE(numout,*) '    ==> EXTRACTING variable ', &
     &              sxy_nam(indsxy)(1:lenv(sxy_nam(indsxy))), &
     &              ' from Vy vector'
            CASE DEFAULT
               GOTO 1000
            END SELECT
!
            kform='(8x,a,i9)'
            WRITE(numout,kform) '- Initial index : ', &
     &                  sxy_ind(indsxy)
            WRITE(numout,kform) '- Final index   : ', &
     &                  sxy_ind(indsxy)+sxy_nbr(indsxy)-1
            WRITE(numout,kform) '- Segment size  : ', &
     &                  sxy_nbr(indsxy)
         ENDDO
      ENDIF
!
! -1.- Loop on variables to write
! -------------------------------
!
      done=.FALSE.
      jsxy=1
      DO WHILE ((.NOT.(done)).AND.(jsxy.LE.sxyend))
!
! Initialize tables of variables to write in the same file
         DO jsxy1=1,nbvar
            txy_indvct(jsxy1) = 0
            txy_indmsk(jsxy1)  = 0
         ENDDO
!
! Find first variable not yet written
         indsxy=sxy_ord(jsxy)
         DO WHILE ((.NOT.(sxy_towrite(indsxy))).AND.(jsxy.LE.sxyend))
            jsxy=jsxy+1
            indsxy=sxy_ord(jsxy)
         ENDDO
!
! -2.- Build tables of variables to write in the same file
! --------------------------------------------------------
! indtxy = index of variable in file
! ! order of variables in file is indifferent in NetCDF files (jextsxy=3)
         jtxy   = 1
         IF (jextsxy.NE.3) THEN
            indtxy = mod(sxyopos(indsxy),100)
         ELSE
            indtxy = jtxy
         ENDIF
         txy_indvct(indtxy) = sxy_ind(indsxy)
         txy_indmsk(indtxy) = jsxy 
         sxy_towrite(indsxy)  = .FALSE.
         DO jsxy1=jsxy+1,sxyend
            indsxy1=sxy_ord(jsxy1)
            IF ((sxyonam(indsxy).EQ.sxyonam(indsxy1)).OR.extsxyunit) THEN
               IF (sxy_jpi(indsxy).NE.sxy_jpi(indsxy1)) GOTO 104
               IF ((sxy_jpj(indsxy).NE.sxy_jpj(indsxy1)) &
     &              .AND.(sxy_dim(indsxy).GE.2)) GOTO 104
               IF ((jextsxy.NE.3).AND.(jextsxy.NE.4)) THEN
                  IF (sxy_dim(indsxy).NE.sxy_dim(indsxy1)) GOTO 104
                  IF ((sxy_jpk(indsxy).NE.sxy_jpk(indsxy1)) &
     &                 .AND.(sxy_dim(indsxy).GE.3)) GOTO 104
                  IF ((sxy_jpt(indsxy).NE.sxy_jpt(indsxy1)) &
     &                 .AND.(sxy_dim(indsxy).GE.4))  GOTO 104
               ELSE
                  IF (sxy_dim(indsxy).GT.4) GOTO 104
               ENDIF
               jtxy=jtxy+1
               IF (jextsxy.NE.3) THEN
                  indtxy = mod(sxyopos(indsxy1),100)
               ELSE
                  indtxy = jtxy
               ENDIF
               txy_indvct(indtxy) = sxy_ind(indsxy1)
               txy_indmsk(indtxy) = jsxy1 
               sxy_towrite(indsxy1) = .FALSE.
            ENDIF
         ENDDO
!
! Check internal coherence of the table
         txyend=jtxy
         IF (.NOT.(extsxyunit)) THEN
            DO jtxy=1,txyend
               jsxy1=txy_indmsk(jtxy)
               indsxy1=sxy_ord(jsxy1)
               IF (txyend.NE.(sxyopos(indsxy1)/100)) GOTO 104   
            ENDDO
         ENDIF
!
! Set parameters to write in file headers
         IF (extsxyunit) THEN
            WRITE (fname,'(a)') kfnoutsxy(1:lenv(kfnoutsxy))
            WRITE (cdrec1,'(a)') 'SESAM output data'
         ELSE
            WRITE (fname,'(a,a,a)') kfnoutsxy(1:(spos-1)), &
     &           sxyonam(indsxy)(1:lenv(sxyonam(indsxy))), &
     &           kfnoutsxy((spos+1):lenv(kfnoutsxy))
            WRITE (cdrec1,'("Variables",1X,a)') &
     &           sxyonam(indsxy)(2:(lenv(sxyonam(indsxy))-1))
         ENDIF
         idast=FREAL(0.)
         WRITE (cdrec2,'(a)') ' '
         WRITE (cdrec3,'(a)') ' '
         WRITE (cdrec4,'(a)') ' '
!
! -3.- Control print
! ------------------
!
         IF (nprint.GE.3) THEN
            WRITE(numout,*) '    ==> WRITING file ', &
     &           fname(1:lenv(fname))
            kform='(8x,2a)'
            WRITE(numout,kform) '- Title : ',cdrec1
            kform='(8x,a,i3)'
            WRITE(numout,kform) '- Number of variables :',txyend
            DO jtxy=1,txyend
               jsxy1=txy_indmsk(jtxy)
               indsxy1=sxy_ord(jsxy1)
               kform='(12x,a,i3)'
               WRITE(numout,kform) '+ variable index     : ',jsxy1
               kform='(12x,2a)'
               WRITE(numout,kform) '+ variable name      : ', &
     &                         sxy_nam(indsxy1)
               kform='(12x,a,i3,a)'
               WRITE(numout,kform) '+ variable dimension :', &
     &                         sxy_dim(indsxy1) ,'D'
               IF (jextsxy.NE.3) THEN
                  kform='(12x,a,i3)'
                  WRITE (numout,kform) ' + position in file: ', &
     &                         mod(sxyopos(indsxy1),100)
               ENDIF
            ENDDO
         ENDIF
!
! Write ensemble of variable fields (in the right format)
         SELECT CASE (jextsxy)
         CASE (1)
            CALL writenbimg(fname,cdrec1,cdrec2,cdrec3,cdrec4,idast, &
     &                  kvects(:),txy_indvct,txy_indmsk,txyend,kflagxyo)
         CASE (2)
            CALL writendimg(fname,cdrec1,idast,kvects(:), &
     &                  txy_indvct,txy_indmsk,txyend,kflagxyo)
         CASE (3)
            CALL writencdf(fname,cdrec1,idast,kvects(:), &
     &                  txy_indvct,txy_indmsk,txyend,kflagxyo)
         CASE (4)
            CALL writenc(fname,kvects(:), &
     &                   txy_indvct,txy_indmsk,txyend,kflagxyo)
         CASE DEFAULT
            GOTO 1000
         END SELECT
!
! -4.- Exit variable loop if all variable are written
! ---------------------------------------------------
!
         done=.TRUE.
         DO jsxy1=1,sxyend
            indsxy1=sxy_ord(jsxy1)
            IF (sxy_towrite(indsxy1)) THEN
               done=.FALSE.
            ENDIF
         ENDDO
      ENDDO
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'hioxyo','writenfil')
 1001 CALL printerror2(0,1001,3,'hioxyo','writenfil')
!
 102  WRITE (texterror,*) 'bad output file name:', &
     &               kfnoutsxy(1:lenv(kfnoutsxy))
      CALL printerror2(0,102,3,'hioxyo','writenfil',comment=texterror)
 104  WRITE (texterror,*) 'incoherent SESAM configuration'
      CALL printerror2(0,104,3,'hioxyo','writenfil',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readvar(kfninvar,kvects,kjns,klectinfo,kflagxyo) 
!---------------------------------------------------------------------
!
!  Purpose : Read Vx or Vy vector from 'var' file
!  -------
!  Method : Call appropriate routine to read the right 'var' format
!  ------
!  Input :     kfninvar  : filename
!  -----       kjns      : block index
!              klectinfo : read or not header of input files
!              kflagxyo  : vector type (Vx/Vy/Vo)
!  Output :    kvects    : 1D vector object (Vx or Vy)
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpx,jpxend,jpnxend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*) :: kfninvar
      BIGREAL, dimension(:), intent(out) ::  kvects
      INTEGER, intent(in) :: kjns
      LOGICAL, intent(in) :: klectinfo
      INTEGER, intent(in) :: kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jextvar
!----------------------------------------------------------------------
      IF (.NOT.(validextvar(kfninvar))) GOTO 101
      jextvar = indext(kfninvar,extvartab,nbextvar)
      IF ((kflagxyo.NE.1).AND.(kflagxyo.NE.2)) GOTO 1000
!
! Read 'var' file
! ---------------
!
      IF ((kflagxyo.EQ.1).AND.(nallmem.EQ.2)) THEN
!
! Read current block of Vx vector
!
         IF (jpnxend.EQ.1) GOTO 1000
         IF (jpx.EQ.1) GOTO 1000     
         IF (jpx.EQ.jpxend) GOTO 1000     
         IF ((extvarmem(indext(kfninvar,extvartab,nbextvar))/10) &
     &    .GT.nallmem) GOTO 104
!
         SELECT CASE (jextvar)
         CASE (5)
! --- Format cpak
            CALL readpartcpak(kfninvar,kvects(:),kjns,klectinfo,kflagxyo)
         CASE (1,2,3,4,8)
! --- Formats bimg,dimg,cdf,rst
            GOTO 104
         CASE DEFAULT
            GOTO 1002
         END SELECT
!
      ELSE
!
! Read full Vx or Vy vector
!
         SELECT CASE (jextvar)
         CASE (1,2,3,4)
! --- Formats bimg,dimg,cdf
            CALL readfil(kfninvar,kvects(:),klectinfo,kflagxyo)
         CASE (5)
! --- Format cpak
            CALL readcpak(kfninvar,kvects(:),klectinfo,kflagxyo)
         CASE (8)
! --- Format rst
            CALL readrst(kfninvar,kvects(:),klectinfo,kflagxyo)
         CASE DEFAULT
            GOTO 1002
         END SELECT
!
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hioxyo','readvar')
 1002 CALL printerror2(0,1002,1,'hioxyo','readvar')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &        kfninvar(1:lenv(kfninvar))
      CALL printerror2(0,101,3,'hixyo','readvar',comment=texterror)
 104  WRITE (texterror,*) 'Cannot read Vx vector block by block', &
     &     ' with such file format: ',kfninvar(1:lenv(kfninvar))
      CALL printerror2(0,104,3,'hixyo','readvar',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writevar(kfnoutvar,kvectx,kjnx) 
!---------------------------------------------------------------------
!
!  Purpose : Write Vx vector in 'var' file
!  -------
!  Method : Call appropriate routine to write the right 'var' format
!  ------
!  Input :   kfnoutvar  : filename
!  -----     kvectx     : 1D vector object (Vx) to write
!            kjnx       : block index
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpx, jpxend, jpnxend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutvar
      BIGREAL, dimension(:), intent(in) :: kvectx
      INTEGER, intent(in) :: kjnx
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jextvar,flagxyo,jjproc
!----------------------------------------------------------------------
      IF (.NOT.(validextvar(kfnoutvar))) GOTO 101
      jextvar = indext(kfnoutvar,extvartab,nbextvar)
      flagxyo=1
!
! Write 'var' file
! ----------------
!
      IF ((flagxyo.EQ.1).AND.(nallmem.EQ.2)) THEN
!
! Write current block of Vx vector
!
         IF (jpnxend.EQ.1) GOTO 1000
         IF (jpx.EQ.1) GOTO 1000
         IF (jpx.EQ.jpxend) GOTO 1000
         IF ((extvarmem(indext(kfnoutvar,extvartab,nbextvar))/10) &
     &    .GT.nallmem) GOTO 104
!
         SELECT CASE (jextvar)
         CASE (5)
! --- Format cpak
            DO jjproc = 0,jpproc-1
              IF (jproc.EQ.jjproc) THEN
                CALL writepartcpak (kfnoutvar,kvectx(:),kjnx)
              ENDIF
#if defined MPI
              CALL MPI_BARRIER(MPI_COMM_WORLD,mpi_code)
#endif
            ENDDO
         CASE (1,2,3,4,8)
! --- Formats bimg,dimg,cdf,rst
            GOTO 104
         CASE DEFAULT
            GOTO 1002
         END SELECT
!
      ELSE
!
! Write full Vx vector
!
         SELECT CASE (jextvar)
         CASE (1,2,3,4)
! --- Formats bimg,dimg,cdf
            CALL writenfil(kfnoutvar,kvectx(:),flagxyo)
         CASE (5)
! --- Format cpak
            CALL writecpak(kfnoutvar,kvectx(:))
         CASE (8)
! --- Format rst
            CALL writerst(kfnoutvar,kvectx(:))
         CASE DEFAULT
            GOTO 1002
         END SELECT
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hioxxyo','writevar')
 1002 CALL printerror2(0,1002,1,'hioxxyo','writevar')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &     kfnoutvar(1:lenv(kfnoutvar))
      CALL printerror2(0,101,3,'hioxxyo','writevar',comment=texterror)
 104  WRITE (texterror,*) 'Cannot write Vx vector block by block', &
     &     ' with such file format: ',kfnoutvar(1:lenv(kfnoutvar))
      CALL printerror2(0,104,3,'hioxxyo','writevar',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readdta(fnindta,kvecty,klectinfo) 
!---------------------------------------------------------------------
!
!  Purpose : Read Vy vector from 'dta' file
!  -------
!  Method : Call appropriate routine to read the right 'dta' format
!  ------
!  Input :   fnindta   : filename
!  -----     klectinfo : read or not header of input files
!  Output :  kvecty    : 1D vector object (Vy)
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: fnindta
      BIGREAL, dimension(:), intent(out) :: kvecty
      LOGICAL, intent(in) :: klectinfo
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jextdta,flagxyo
!----------------------------------------------------------------------
      IF (.NOT.(validextdta(fnindta))) GOTO 101
      flagxyo=2
      jextdta=indext(fnindta,extdtatab,nbextdta)
!
! Read Vy vector
! --------------
!
      SELECT CASE (jextdta)
      CASE (1,2,3,4)
! Formats bdta,dta,cdta
         CALL readfil(fnindta,kvecty(:),klectinfo,flagxyo)
      CASE DEFAULT
         GOTO 1002
      END SELECT
!
      RETURN
!
! --- error management
!
 1002 CALL printerror2(0,1002,1,'hioxyo','readdta')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &       fnindta(1:lenv(fnindta))
      CALL printerror2(0,101,3,'hixyo','readdta',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writedta(kfnoutdta,kvecty) 
!---------------------------------------------------------------------
!
!  Purpose : Write Vy vector in 'dta' file
!  -------
!  Method : Call appropriate routine to write the right 'dta' format
!  ------
!  Input :  kfnoutdta  : filename
!  -----    kvecty     : 1D vector object (Vy) to write
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutdta
      BIGREAL, dimension(:), intent(in) :: kvecty
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jextdta,flagxyo
!----------------------------------------------------------------------
      IF (.NOT.(validextdta(kfnoutdta))) GOTO 101
      jextdta = indext(kfnoutdta,extdtatab,nbextdta)
      flagxyo=2
!
! Write Vy vector
! ---------------
!
      SELECT CASE (jextdta)
      CASE (1,2,3,4)
! Formats bdta,dta,cdta
         CALL writenfil(kfnoutdta,kvecty(:),flagxyo)
      CASE DEFAULT
         GOTO 1002
      END SELECT
!
      RETURN
!
! --- error management
!
 1002 CALL printerror2(0,1002,1,'hioxyo','writedta')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &             kfnoutdta(1:lenv(kfnoutdta))
       CALL printerror2(0,101,3,'hioxyo','writedta',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrobs (kfninobs,jpoendout,jpitpendout)
!---------------------------------------------------------------------
!
!  Purpose : Read observation file header
!  --------
!  Method : Loop on observation databases and call appropriate
!  ------   routine to read individual observation files
!
!  Input :  kfninobs : filename
!  -----
!  Output : jpoendout   : number of observation values
!  ------   jpitpendout : number of interpolation points
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      INTEGER, intent(out) :: jpoendout,jpitpendout
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: fname
      INTEGER :: jextobs,jpoloc,jpitploc,jobs,indobs,inddbs,spos
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : ./evalhdrobs :'
         WRITE(numout,*) '         read observation file header'
      ENDIF
!
! Select the right 'obs' file format
      IF (.NOT.(validextobs(kfninobs))) GOTO 101
      jextobs=indext(kfninobs,extobstab,nbextobs)
      SELECT CASE (jextobs)
      CASE (1,2)
! --- obs format
         jpoloc=0
         jpitploc=0
!
! Loop over observation databases
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            IF (extobsunit(jextobs)) GOTO 1002
!
! Build individual observation filename
            spos = posit(kfninobs,etoile)
            WRITE (fname,'(a,a,a)') kfninobs(1:(spos-1)), &
     &           obsinam(indobs,inddbs) &
     &           (1:lenv(obsinam(indobs,inddbs))), &
     &           kfninobs((spos+1):lenv(kfninobs))
!
! Read header of each individual observation file
            SELECT CASE (jextobs)
            CASE (1)
              CALL evalhdrfileobs (fname,indobs,jpoloc,jpitploc)
            CASE (2)
              CALL evalhdrfilecobs(fname,indobs,jpoloc,jpitploc)
            END SELECT
            obs_nbr(indobs,inddbs)=jpoloc
            obs_itp(indobs,inddbs)=jpitploc
         ENDDO
!
! Compute header values for full observation vector
! jpoendout   = sum of individual file contributions
! jpitpendout = max of individual file contributions
         indobs=obs_ord(1)
         inddbs=obsnord(1)
         obs_ind(indobs,inddbs)=1
         jpoendout=obs_nbr(indobs,inddbs)
         jpitpendout=obs_itp(indobs,inddbs)
         DO jobs=2,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            obs_ind(indobs,inddbs)=jpoendout+1
            jpoendout=obs_nbr(indobs,inddbs)+jpoendout
            jpitpendout=MAX(jpitpendout,obs_itp(indobs,inddbs))
         ENDDO
      CASE DEFAULT
         GOTO 1002
      END SELECT
!
      RETURN
!
! --- error management
!
 1002 CALL printerror2(0,1002,1,'hioxyo','evalhdrobs')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &          kfninobs(1:lenv(kfninobs))
      CALL printerror2(0,101,3,'hioxyo','evalhdrobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readvalobs(kfninobs,kvecto)
!---------------------------------------------------------------------
!
!  Purpose : Read Vo vector from 'obs' file (observed value only)
!  --------
!  Method : Loop on observation databases and call appropriate
!  ------   routine to read individual observation files
!
!  Input :  kfninobs  : filename
!  -----
!  Output : kvecto    : 1D vector object (Vo)
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      BIGREAL, dimension(:), intent(out) :: kvecto
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: fname
      INTEGER :: jposize
      INTEGER :: jextobs
      INTEGER :: jobs,indobs,inddbs,jdim,jrec,jobsdeb
      INTEGER :: jodebloc,jofinloc,jpoloc,jpoendloc
      INTEGER :: jo,jo1,jodeb,jofin,spos
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readvalobs :'
         WRITE(numout,*) '         read observation values'
      ENDIF
! Get size of input array
      jposize = size(kvecto,1)
! Check file extension
      IF (.NOT.(validextobs(kfninobs))) GOTO 101
      jextobs=indext(kfninobs,extobstab,nbextobs)
!
! Select the right 'obs' file format
      SELECT CASE (jextobs)
      CASE (1,2)
! --- obs format
         jodebloc=0
         jofinloc=0
         jpoloc=0
         jpoendloc=0
!
! Loop over observation databases
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            IF (extobsunit(jextobs)) GOTO 1002
!
! Build individual observation filename
            spos = posit(kfninobs,etoile)
            WRITE (fname,'(a,a,a)') kfninobs(1:(spos-1)), &
     &           obsinam(indobs,inddbs) &
     &           (1:lenv(obsinam(indobs,inddbs))), &
     &           kfninobs((spos+1):lenv(kfninobs))
            jodebloc=obs_ind(indobs,inddbs)
            jofinloc=obs_ind(indobs,inddbs)-1+obs_nbr(indobs,inddbs)
            jpoloc=obs_nbr(indobs,inddbs)
!
! Read each individual observation file
            IF (jofinloc.GE.jodebloc) THEN
               SELECT CASE (jextobs)
               CASE (1)
                 CALL readvalfileobs (fname,kvecto(jodebloc:jofinloc),jobs)
               CASE (2)
                 CALL readvalfilecobs(fname,kvecto(jodebloc:jofinloc),jobs)
               END SELECT
            ENDIF
            jpoendloc=jpoloc+jpoendloc
         ENDDO
      CASE DEFAULT
         GOTO 1002
      END SELECT
!
! Check coherence with array sizes
      IF (jpoendloc.NE.jposize) GOTO 102
!
! Control Print
      IF (nprint.GE.4) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            WRITE(numout,*) ' ==> Observation database : ', &
     &           obs_nam(indobs,inddbs)(1:lenv(obs_nam(indobs,inddbs)))
            IF (obs_nbr(indobs,inddbs).GE.1) THEN
               jodeb=obs_ind(indobs,inddbs)
               jofin=obs_ind(indobs,inddbs)-1+obs_nbr(indobs,inddbs)
               jofin=MIN(jofin,jposize)
               DO jo=jodeb,jofin,10
                  WRITE(numout,'(a,10(1X,i5.5))')' ind |', &
     &                 (jo1,jo1=jo,MIN(jo+9,jofin))
                  WRITE(numout,'(a,10(1X,F5.1))')' val |', &
     &                 kvecto(jo:MIN(jo+9,jofin))
               ENDDO
            ELSE
               WRITE(numout,*) ' empty observation database'
            ENDIF
         ENDDO
      ENDIF
!
      RETURN
!
! --- error management
!
 1002 CALL printerror2(0,1002,1,'hioxyo','readvalobs')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &          kfninobs(1:lenv(kfninobs))
      CALL printerror2(0,101,3,'hioxyo','readvalobs',comment=texterror)
 102  WRITE (texterror,*) 'Incoherence between input arrays'
      CALL printerror2(0,102,3,'hioxyo','readvalobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcfgobs(kfninobs,kflagcfg,kvectorms, &
     &     kgridijkobs,kposcoefobs)
!---------------------------------------------------------------------
!
!  Purpose : Read Vo vector from 'obs' file
!  -------   (observation configuration only)
!
!  Method : Loop on observation databases and call appropriate
!  ------   routine to read individual observation files
!
!  Input :  kfninobs    : filename
!  -----
!  Output : kflagcfg    : what element of the configuration to read
!  ------                 (1=kvectorms, 2=kgridijkobs, 3=kposcoefobs)
!                         (obsolete parameter)
!           kvectorms   : associated error value (obsolete parameter)
!           kgridijkobs : observation location (x,y,z)
!           kposcoefobs : observation operator (interpolation points
!                         and interpolation coefficients)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      INTEGER, intent(in) :: kflagcfg
      BIGREAL, optional, dimension(:), intent(out) :: kvectorms
      TYPE (type_gridijk), optional, dimension(:), intent(out) ::  &
     &     kgridijkobs
      TYPE (type_poscoef), optional, dimension(:,:), intent(out) ::  &
     &     kposcoefobs
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: fname
      INTEGER :: jposize,jpitpsize
      INTEGER :: jextobs
      INTEGER :: jobs,indobs,inddbs,jdim,jrec,jitp,jobsdeb
      INTEGER :: jodebloc,jofinloc,jpoloc,jpitploc,jpoendloc,jpitpendloc
      INTEGER :: jo,jo1,jodeb,jofin,spos
!----------------------------------------------------------------------
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readcfgobs :'
         WRITE(numout,*) '         read observation configuration'
      ENDIF
!
! Select element of the configuration to read 
! Set vector sizes
      jpitpsize = 0
      IF (present(kvectorms)) THEN
        jposize = size(kvectorms,1)
        DO jobs=1,obsend
           indobs=obs_ord(jobs)
           inddbs=obsnord(jobs)
           jpitpsize = MAX(jpitpsize,obs_itp(indobs,inddbs))
        ENDDO
      ENDIF
      IF (present(kgridijkobs)) THEN
        jposize = size(kgridijkobs,1)
        DO jobs=1,obsend
           indobs=obs_ord(jobs)
           inddbs=obsnord(jobs)
           jpitpsize = MAX(jpitpsize,obs_itp(indobs,inddbs))
        ENDDO
      ENDIF
      IF (present(kposcoefobs)) THEN
        jposize = size(kposcoefobs,1)
        jpitpsize = size(kposcoefobs,2)
        DO jobs=1,obsend
           indobs=obs_ord(jobs)
           inddbs=obsnord(jobs)
           IF (jpitpsize.LT.obs_itp(indobs,inddbs)) GOTO 102
        ENDDO
      ENDIF
!
! Check file extension
      IF (.NOT.(validextobs(kfninobs))) GOTO 101
      jextobs=indext(kfninobs,extobstab,nbextobs)
!
! Select the right 'obs' file format
      SELECT CASE (jextobs)
      CASE (1,2)
! --- obs format
         jodebloc=0
         jofinloc=0
         jpoloc=0
         jpoendloc=0
         jpitploc=0
         jpitpendloc=0
!
! Loop over observation databases
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            IF (extobsunit(jextobs)) GOTO 1002
!
! Build individual observation filename
            spos = posit(kfninobs,etoile)
            WRITE (fname,'(a,a,a)') kfninobs(1:(spos-1)), &
     &           obsinam(indobs,inddbs) &
     &           (1:lenv(obsinam(indobs,inddbs))), &
     &           kfninobs((spos+1):lenv(kfninobs))
            jodebloc=obs_ind(indobs,inddbs)
            jofinloc=obs_ind(indobs,inddbs)-1+obs_nbr(indobs,inddbs)
            jpoloc=obs_nbr(indobs,inddbs)
            jpitploc=obs_itp(indobs,inddbs)
!
! Read each individual observation file
            IF (jofinloc.GE.jodebloc) THEN
               IF (present(kvectorms)) kvectorms(jodebloc:jofinloc)=FREAL(0.0)
               SELECT CASE (jextobs)
               CASE (1)
                  IF (present(kgridijkobs)) CALL readcfgfileobs (fname, &
     &               jobs,2,kgridijkobs=kgridijkobs(jodebloc:jofinloc))
                  IF (present(kposcoefobs)) CALL readcfgfileobs (fname, &
     &               jobs,3,kposcoefobs=kposcoefobs(jodebloc:jofinloc,:jpitploc))
               CASE (2)
                  IF (present(kgridijkobs)) CALL readcfgfilecobs(fname, &
     &               jobs,2,kgridijkobs=kgridijkobs(jodebloc:jofinloc))
                  IF (present(kposcoefobs)) CALL readcfgfilecobs(fname, &
     &               jobs,3,kposcoefobs=kposcoefobs(jodebloc:jofinloc,:jpitploc))
               END SELECT
            ENDIF
            jpoendloc=jpoloc+jpoendloc
            jpitpendloc=MAX(jpitploc,jpitpendloc)
         ENDDO
      CASE DEFAULT
         GOTO 1002
      END SELECT
!
! Check coherence with array sizes
      IF (jpoendloc.NE.jposize) GOTO 102
      IF (present(kposcoefobs).AND.(jpitpendloc.NE.jpitpsize)) GOTO 102
!
! Control Print
      IF (nprint.GE.3) THEN
         WRITE(numout,30) '-'
         WRITE(numout,40) 'obs','dta','dbs','name', &
     &                    'dim','size','ref','std'
         WRITE(numout,30) '-'
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            WRITE(numout,50) jobs,obs_ord(jobs),obsnord(jobs), &
     &           obs_nam(indobs,inddbs), &
     &           obs_dim(indobs,inddbs),obs_nbr(indobs,inddbs), &
     &           obs_moy(indobs,inddbs),obs_ect(indobs,inddbs)
         ENDDO
         WRITE(numout,30) '-'
      ENDIF
!
      IF (nprint.GE.4) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            WRITE(numout,*) ' '
            WRITE(numout,*) ' ==> Observation database : ', &
     &           obs_nam(indobs,inddbs)(1:lenv(obs_nam(indobs,inddbs)))
            IF (obs_nbr(indobs,inddbs).GE.1) THEN
               jodeb=obs_ind(indobs,inddbs)
               jofin=obs_ind(indobs,inddbs)-1+obs_nbr(indobs,inddbs)
               jofin=MIN(jofin,jposize)
               DO jo=jodeb,jofin,10
                  WRITE(numout,'(a,10(1X,i5.5))')' ind |', &
     &                 (jo1,jo1=jo,MIN(jo+9,jofin))
                  IF (present(kgridijkobs)) THEN
                     WRITE(numout,'(a,10(1X,F5.1))')' lon |', &
     &                    kgridijkobs(jo:MIN(jo+9,jofin))%longi
                     WRITE(numout,'(a,10(1X,F5.1))')' lat |', &
     &                    kgridijkobs(jo:MIN(jo+9,jofin))%latj
                     WRITE(numout,'(a,10(1X,F5.1))')' lev |', &
     &                    kgridijkobs(jo:MIN(jo+9,jofin))%levk
                  ENDIF
                  IF (present(kposcoefobs)) THEN
                     DO jitp=1,obs_itp(indobs,inddbs)
                        WRITE(numout,'(a,i1.1,a,10(1X,i5.5))') ' pos', &
     &                      jitp,'|',kposcoefobs(jo:MIN(jo+9,jofin),jitp)%pos
                        WRITE(numout,'(a,i1.1,a,10(1X,F5.1))') ' cof', &
     &                      jitp,'|',kposcoefobs(jo:MIN(jo+9,jofin),jitp)%coef
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               WRITE(numout,*) ' empty observation database'
            ENDIF
         ENDDO
      ENDIF
!
      RETURN
!
! --- format definitions
!
 10   FORMAT ("|",7(2X,A4,2X,"|"))
 20   FORMAT ("|",2(I5,3X,"|"),3X,A5,"|",2X,I3,"D",2X,"|", &
     &   I7,1X,"|",2(E8.1,1X,"|"))
 30   FORMAT (A1,64("-"))
 40   FORMAT ("|",8(2X,A4,2X,"|"))
 50   FORMAT ("|",3(I5,3X,"|"),3X,A5,"|",2X,I3,"D",2X,"|", &
     &   I7,1X,"|",2(E8.1,1X,"|"))
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hioxyo','readcfgobs')
 1002 CALL printerror2(0,1002,1,'hioxyo','readcfgobs')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &          kfninobs(1:lenv(kfninobs))
      CALL printerror2(0,101,3,'hioxyo','readcfgobs',comment=texterror)
 102  WRITE (texterror,*) 'Incoherence between input arrays'
      CALL printerror2(0,102,3,'hioxyo','readcfgobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeobs(kfnoutobs,kvecto,kvectorms,kgridijkobs, &
     &     kposcoefobs)
!---------------------------------------------------------------------
!
!  Purpose : Write Vo vector in 'obs' file
!  -------
!  Method : Loop on observation databases and call appropriate
!  ------   routine to write individual observation files
!
!  Input :  kfnoutobs  : filename  
!  -----    kvecto     : 1D vector object (Vo)
!           kvectorms   : associated error value [obsolete]
!           kgridijkobs : observation location (x,y,z)
!           kposcoefobs : observation operator (interpolation points
!                         and interpolation coefficients)
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutobs
      BIGREAL, dimension(:), intent(in) :: kvecto,kvectorms
      TYPE (type_gridijk), dimension(:), intent(in)  :: kgridijkobs
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: fname
      INTEGER :: jextobs
      INTEGER :: jposize,jpitpsize,jobs,indobs,inddbs
      INTEGER :: jodebloc,jofinloc,jpoloc,jpitploc
      INTEGER :: jdim,jrec,jinterp,jobsdeb,spos
      INTEGER :: jo,jo1,jodeb,jofin,jitp,jpoendloc,jpitpendloc
!----------------------------------------------------------------------
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./writeobs :'
         WRITE(numout,*) '         write observation object'
      ENDIF
! Get size of input arrays
      jposize = size(kvecto,1)
      jpitpsize = size(kposcoefobs,2)
! Test coherence of input arrays
      IF (jposize.NE.size(kvectorms,1)) GOTO 1000
      IF (jposize.NE.size(kposcoefobs,1)) GOTO 1000
! Check file extension
      IF (.NOT.(validextobs(kfnoutobs))) GOTO 101
      jextobs=indext(kfnoutobs,extobstab,nbextobs)
!
! Select the right 'obs' file format
      SELECT CASE (jextobs)
      CASE (1,2)
! --- obs format
         jodebloc=0
         jpoloc=0
         jpitploc=0
         jpoendloc=0
         jpitpendloc=0
!
! Loop over observation databases
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            IF (extobsunit(jextobs)) GOTO 1002
!
! Build individual observation filename
            spos = posit(kfnoutobs,etoile)
            WRITE (fname,'(a,a,a)') kfnoutobs(1:(spos-1)), &
     &           obsonam(indobs,inddbs)(1:lenv(obsonam(indobs,inddbs))), &
     &           kfnoutobs((spos+1):lenv(kfnoutobs))
            jodebloc=obs_ind(indobs,inddbs)
            jofinloc=obs_ind(indobs,inddbs)-1+obs_nbr(indobs,inddbs)
            jpoloc=obs_nbr(indobs,inddbs)
            jpitploc=obs_itp(indobs,inddbs)
!
! Write each individual observation file
            IF (jpoloc.NE.(jofinloc-jodebloc+1)) GOTO 1000
            IF (jofinloc.LT.jodebloc) THEN
              SELECT CASE (jextobs)
              CASE (1)
                CALL writehdrfileobs (fname,jobs)
              CASE (2)
                CALL writehdrfilecobs(fname,jobs)
              END SELECT
            ELSE
              SELECT CASE (jextobs)
              CASE (1)
                CALL writefileobs (fname, &
     &               kvecto(jodebloc:jofinloc), &
     &               kvectorms(jodebloc:jofinloc), &
     &               kgridijkobs(jodebloc:jofinloc), &
     &               kposcoefobs(jodebloc:jofinloc,:jpitploc),jobs)
              CASE (2)
                CALL writefilecobs(fname, &
     &               kvecto(jodebloc:jofinloc), &
     &               kvectorms(jodebloc:jofinloc), &
     &               kgridijkobs(jodebloc:jofinloc), &
     &               kposcoefobs(jodebloc:jofinloc,:jpitploc),jobs)
              END SELECT
            ENDIF
            jpoendloc=jpoloc+jpoendloc
            jpitpendloc=MAX(jpitploc,jpitpendloc)
         ENDDO
      CASE DEFAULT
         GOTO 1002
      END SELECT
!
! Check coherence with array sizes
      IF ((jpoendloc.NE.jposize).OR.(jpitpendloc.NE.jpitpsize)) GOTO 1000
!
! Control Print
      IF (nprint.GE.4) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            WRITE(numout,*) ' ==> Observation database : ', &
     &           obs_nam(indobs,inddbs)(1:lenv(obs_nam(indobs,inddbs)))
            IF (obs_nbr(indobs,inddbs).GE.1) THEN
               jodeb=obs_ind(indobs,inddbs)
               jofin=obs_ind(indobs,inddbs)-1+obs_nbr(indobs,inddbs)
               jofin=MIN(jofin,jposize)
               DO jo=jodeb,jofin,10
                  WRITE(numout,'(a,10(1X,i5.5))')' ind |', &
     &                 (jo1,jo1=jo,MIN(jo+9,jofin))
                  WRITE(numout,'(a,10(1X,F5.1))')' lon |', &
     &                 kgridijkobs(jo:MIN(jo+9,jofin))%longi
                  WRITE(numout,'(a,10(1X,F5.1))')' lat |', &
     &                 kgridijkobs(jo:MIN(jo+9,jofin))%latj
                  WRITE(numout,'(a,10(1X,F5.1))')' lev |', &
     &                 kgridijkobs(jo:MIN(jo+9,jofin))%levk
                  WRITE(numout,'(a,10(1X,F5.1))')' val |', &
     &                 kvecto(jo:MIN(jo+9,jofin))
                  WRITE(numout,'(a,10(1X,F5.2))')' rms |', &
     &                 kvectorms(jo:MIN(jo+9,jofin))
                  DO jitp=1,obs_itp(indobs,inddbs)
                     WRITE(numout,'(a,i1.1,a,10(1X,i5.5))') ' pos', &
     &                    jitp,'|',kposcoefobs(jo:MIN(jo+9,jofin),jitp)%pos
                     WRITE(numout,'(a,i1.1,a,10(1X,F5.1))') ' cof', &
     &                    jitp,'|',kposcoefobs(jo:MIN(jo+9,jofin),jitp)%coef
                  ENDDO
               ENDDO
            ELSE
               WRITE(numout,*) ' empty observation database'
            ENDIF
         ENDDO
      ENDIF
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hioxyo','writeobs')
 1002 CALL printerror2(0,1002,1,'hioxyo','writeobs')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &       kfnoutobs(1:lenv(kfnoutobs))
      CALL printerror2(0,101,3,'hioxyo','writeobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writesingleobs(kfnoutobs, &
     &     kvecto,kvectorms,kgridijkobs,kposcoefobs,kjobs)
!---------------------------------------------------------------------
!
!  Purpose : Write segment of Vo vector in 'obs' file
!  -------   (corresponding to one observation database)
!
!  Method : Select observation database and call appropriate
!  ------   routine to write corresponding observation file
!
!  Input : kfnoutobs   : filename  
!  -----   kvecto      : 1D vector object (Vo)
!          kvectorms   : associated error value (obsolete)
!          kgridijkobs : observation location (x,y,z)
!          kposcoefobs : observation operator (interpolation points
!                         and interpolation coefficients)
!          kjobs       : observation index
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutobs
      BIGREAL, dimension(:), intent(in) :: kvecto,kvectorms
      TYPE (type_gridijk), dimension(:), intent(in)  :: kgridijkobs
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
      INTEGER, intent(in) :: kjobs
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: fname
      INTEGER :: jextobs
      INTEGER :: jposizeloc,jpitpsizeloc,indobs,inddbs
      INTEGER :: jodebloc,jofinloc,jpoloc,jpitploc
      INTEGER :: jdim,jrec,jinterp,jobsdeb,spos
      INTEGER :: jo,jo1,jodeb,jofin,jitp,jpoendloc,jpitpendloc
!----------------------------------------------------------------------
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./writesingleobs :'
         WRITE(numout,*) '         write single observation file'
      ENDIF
! Get size of input arrays
      jposizeloc = size(kvecto,1)
      jpitpsizeloc = size(kposcoefobs,2)
! Test coherence of input arrays
      IF (jposizeloc.NE.size(kvectorms,1)) GOTO 1000
      IF (jposizeloc.NE.size(kposcoefobs,1)) GOTO 1000
! Check file extension
      IF (.NOT.(validextobs(kfnoutobs))) GOTO 101
      jextobs=indext(kfnoutobs,extobstab,nbextobs)
!
! Build individual observation filename
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
      jodebloc=1
      jofinloc=jposizeloc
      jpitploc=jpitpsizeloc 
      IF (extobsunit(jextobs)) THEN
         WRITE (fname,'(a)') kfnoutobs
      ELSE
         spos=posit(kfnoutobs,etoile)
         WRITE (fname,'(a,a,a)') kfnoutobs(1:(spos-1)), &
     &        obsinam(indobs,inddbs) &
     &        (1:lenv(obsinam(indobs,inddbs))), &
     &        kfnoutobs((spos+1):lenv(kfnoutobs))
      ENDIF 
!
! Select the right 'obs' file format
      SELECT CASE (jextobs)
      CASE (1)
! --- obs format
         CALL writefileobs (fname,kvecto(jodebloc:jofinloc), &
     &        kvectorms(jodebloc:jofinloc), &
     &        kgridijkobs(jodebloc:jofinloc), &
     &        kposcoefobs(jodebloc:jofinloc,:jpitploc),kjobs)
      CASE (2)
! --- cobs format
         CALL writefilecobs(fname,kvecto(jodebloc:jofinloc), &
     &        kvectorms(jodebloc:jofinloc), &
     &        kgridijkobs(jodebloc:jofinloc), &
     &        kposcoefobs(jodebloc:jofinloc,:jpitploc),kjobs)
      CASE DEFAULT
         GOTO 1002
      END SELECT
!
! Control Print
      IF (nprint.GE.4) THEN
         WRITE(numout,*) ' ==> Observation database : ', &
     &        obs_nam(indobs,inddbs)(1:lenv(obs_nam(indobs,inddbs)))
         IF (obs_nbr(indobs,inddbs).GE.1) THEN
            jodeb=obs_ind(indobs,inddbs)
            jofin=obs_ind(indobs,inddbs)-1+obs_nbr(indobs,inddbs)
            jofin=MIN(jofin,jposizeloc)
            DO jo=jodeb,jofin,10
               WRITE(numout,'(a,10(1X,i5.5))')' ind |', &
     &              (jo1,jo1=jo,MIN(jo+9,jofin))
               WRITE(numout,'(a,10(1X,F5.1))')' lon |', &
     &              kgridijkobs(jo:MIN(jo+9,jofin))%longi
               WRITE(numout,'(a,10(1X,F5.1))')' lat |', &
     &              kgridijkobs(jo:MIN(jo+9,jofin))%latj
               WRITE(numout,'(a,10(1X,F5.1))')' lev |', &
     &              kgridijkobs(jo:MIN(jo+9,jofin))%levk
               WRITE(numout,'(a,10(1X,F5.1))')' val |', &
     &              kvecto(jo:MIN(jo+9,jofin))
               WRITE(numout,'(a,10(1X,F5.2))')' rms |', &
     &              kvectorms(jo:MIN(jo+9,jofin))
               DO jitp=1,obs_itp(indobs,inddbs)
                  WRITE(numout,'(a,i1.1,a,10(1X,i5.5))') ' pos', &
     &                 jitp,'|',kposcoefobs(jo:MIN(jo+9,jofin),jitp)%pos
                  WRITE(numout,'(a,i1.1,a,10(1X,F5.1))') ' cof', &
     &                 jitp,'|',kposcoefobs(jo:MIN(jo+9,jofin),jitp)%coef
               ENDDO
            ENDDO
         ELSE
            WRITE(numout,*) ' empty observation database'
         ENDIF
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hioxyo','writesingleobs')
 1002 CALL printerror2(0,1002,1,'hioxyo','writesingleobs')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &       kfnoutobs(1:lenv(kfnoutobs))
      CALL printerror2(0,101,3,'hioxyo','writesingleobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hioxyo
