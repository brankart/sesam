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
! ---                   HIOMSK.F90                                ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12  (C.E. Testut)                       ---
! --- modification : 99-05  (C.E. Testut)                       ---
! --- modification : 01-06  (C.E. Testut)                       ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  readmsk : Read SESAM configuration masks
! --- SUBROUTINE  evalhdrmsk : Read mask file header
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE hiomsk
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC readmsk,evalhdrmsk

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readmsk
!---------------------------------------------------------------------
!
!  Purpose : Read SESAM configuration masks
!  -------
!  Method : Call specific routines as a function
!  ------   of mask file format (bimg,dimg,cdf)
!
!---------------------------------------------------------------------
! module
! ======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use liobimg
      use liodimg
      use liocdf
      use lionc
      use utilvalid
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: fname
      INTEGER :: allocok,ji,jj,jk,jt
      INTEGER :: jext,spos,kjpi,kjpj,kjpk,kjpt
      INTEGER :: jvar,indvar,jdta,inddta,flagxyo
      LOGICAL :: warning
      LOGICAL :: nodtamsk=.FALSE.,novarmsk=.FALSE.
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/readmsk :'
         WRITE(numout,*) '         read SESAM configuration masks'
      ENDIF
!
      jext=1
!
! Do not use default Vx mask files if optional switch -varmsk is present
!
      IF (largvarmsk) THEN
         IF (argvarmsk(1:lenv(argvarmsk)).EQ.'nomask') THEN
            novarmsk=.TRUE.
         ELSE
            DO jvar =1,varend
               indvar = var_ord(jvar)
               IF (.NOT.(validextvar(argvarmsk))) GOTO 1000
               spos = posit(argvarmsk,etoile)
               jext = indext(argvarmsk,extvartab,nbextvar)
               IF ((jext.NE.1).AND.(jext.NE.2).AND.(jext.NE.3)) GOTO 101
               WRITE (fname,'(a,a,a)') argvarmsk(1:(spos-1)), &
     &                         varinam(indvar)(1:lenv(varinam(indvar))), &
     &                         argvarmsk((spos+1):lenv(argvarmsk))
               varfmsk(indvar)= fname(1:lenv(fname))
               varemsk(indvar)= jext
               varmsea(indvar)= .FALSE.
            ENDDO         
         ENDIF
      ENDIF
!
! Do not use default Vy mask files if optional switch -dtamsk is present
!
      IF (largdtamsk) THEN
         IF (argdtamsk(1:lenv(argdtamsk)).EQ.'nomask') THEN
            nodtamsk=.TRUE.
         ELSE
            DO jdta =1,dtaend
               inddta = dta_ord(jdta)
               IF (.NOT.(validextdta(argdtamsk))) GOTO 1000
               spos = posit(argdtamsk,etoile)
               jext = indext(argdtamsk,extdtatab,nbextdta)
               IF ((jext.NE.1).AND.(jext.NE.2).AND.(jext.NE.3)) GOTO 102
               WRITE (fname,'(a,a,a)') argdtamsk(1:(spos-1)), &
     &                         dtainam(inddta)(1:lenv(dtainam(inddta))), &
     &                         argdtamsk((spos+1):lenv(argdtamsk))
               dtafmsk(inddta)= fname(1:lenv(fname))
               dtaemsk(inddta)= jext
               dtamsea(inddta)= .FALSE.
            ENDDO         
         ENDIF
      ENDIF
!
! -1.- Read size of configuration variables (jpi,jpj,jpk,jpt)
! -----------------------------------------------------------
!   var_jpi()  ----> Vx object, size of 1st dimension (i.e. longitude)
!   var_jpj()  ----> Vx object, size of 2nd dimension (i.e. latitude)
!   var_jpk()  ----> Vx object, size of 3rd dimension (i.e. depth)
!   var_jpt()  ----> Vx object, size of 4th dimension (i.e. time)
!   dta_jpi()  ----> Vy object, size of 1st dimension (i.e. longitude)
!   dta_jpj()  ----> Vy object, size of 2nd dimension (i.e. latitude)
!   dta_jpk()  ----> Vy object, size of 3rd dimension (i.e. depth)
!   dta_jpt()  ----> Vy object, size of 4th dimension (i.e. time)
!
! -1.1- Vx dimensions
!
      flagxyo = 1
      DO jvar=1,varend
         indvar=var_ord(jvar)
         kjpi=1
         kjpj=1
         kjpk=1
         kjpt=1
         CALL evalhdrmsk(indvar,kjpi,kjpj,kjpk,kjpt,flagxyo)
!
         var_jpi(indvar)=kjpi
         var_jpj(indvar)=kjpj
         var_jpk(indvar)=kjpk
         var_jpt(indvar)=kjpt
         var_nbr(indvar)= var_jpi(indvar)*var_jpj(indvar) &
     &        *var_jpk(indvar)*var_jpt(indvar)
      ENDDO
!
! -1.2- Vy dimensions
!
      flagxyo = 2
      DO jdta=1,dtaend
         inddta=dta_ord(jdta)
         kjpi=1
         kjpj=1
         kjpk=1
         kjpt=1
         CALL evalhdrmsk(inddta,kjpi,kjpj,kjpk,kjpt,flagxyo)
!
         dta_jpi(inddta)=kjpi
         dta_jpj(inddta)=kjpj
         dta_jpk(inddta)=kjpk
         dta_jpt(inddta)=kjpt
         dta_nbr(inddta)= dta_jpi(inddta)*dta_jpj(inddta) &
     &        *dta_jpk(inddta)*dta_jpt(inddta)
      ENDDO
!
! -1.3- Coherence test
!
      DO jdta=1,dtaend
         inddta=dta_ord(jdta)
         IF (dta_jpi(inddta).GT.var_jpi(inddta)) GOTO 103
         IF (dta_jpj(inddta).GT.var_jpj(inddta)) GOTO 103
         IF (dta_jpk(inddta).GT.var_jpk(inddta)) GOTO 103
         IF (dta_jpt(inddta).GT.var_jpt(inddta)) GOTO 103
      ENDDO
!
! -1.4- Print information about variables dimensions
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*)  'Variable dimensions in Vx object (var)'
         WRITE(numout,30) '-'
         WRITE(numout,10) 'ord','ind','var','jpi','jpj','jpk','jpt'
         WRITE(numout,30) '-'
         DO jvar=1,varend
            indvar=var_ord(jvar)
            WRITE(numout,20) jvar,var_ord(jvar), &
     &           var_nam(indvar), &
     &           var_jpi(indvar),var_jpj(indvar), &
     &           var_jpk(indvar),var_jpt(indvar)
         ENDDO
         WRITE(numout,30) '-'
         WRITE(numout,*)  'Variable dimensions in Vy object (dta)'
         WRITE(numout,30) '-'
         WRITE(numout,10) 'ord','ind','dta','jpi','jpj','jpk','jpt'
         WRITE(numout,30) '-'
         DO jdta=1,dtaend
            inddta=dta_ord(jdta)
            WRITE(numout,20) jdta,dta_ord(jdta), &
     &           dta_nam(inddta), &
     &           dta_jpi(inddta),dta_jpj(inddta), &
     &           dta_jpk(inddta),dta_jpt(inddta)
         ENDDO
         WRITE(numout,30) '-'
      ENDIF
!
! -1.5- Warning : Vx and Vy of different dimensions
!
      warning = .FALSE.
      DO jdta=1,dtaend
         inddta=dta_ord(jdta)
         warning= ( (warning) .AND. &
     &        ((dta_jpi(inddta).NE.var_jpi(inddta)) &
     &        .AND.(dta_dim(inddta).GE.1))  )
         warning= ( (warning) .AND. &
     &        ((dta_jpj(inddta).NE.var_jpj(inddta)) &
     &        .AND.(dta_dim(inddta).GE.2))  )
         warning= ( (warning) .AND. &
     &        ((dta_jpk(inddta).NE.var_jpk(inddta)) &
     &        .AND.(dta_dim(inddta).GE.3))  )
         warning= ( (warning) .AND. &
     &        ((dta_jpt(inddta).NE.var_jpt(inddta)) &
     &        .AND.(dta_dim(inddta).GE.4))  )
         IF (warning) THEN
            PRINT *,'----WARNING----'
            PRINT *,dta_nam(inddta)(1:lenv(dta_nam(inddta))) &
     &           ,'(jpidta,jpjdta,jpkdta,jptdta)', &
     &           'different from ', &
     &           var_nam(inddta)(1:lenv(var_nam(inddta))) &
     &           ,'(jpivar,jpjvar,jpkvar,jptvar)'
         ENDIF
      ENDDO
!
! -1.6- Define mask dimensions and allocate mask arrays
!
      kjpi=1
      kjpj=1
      kjpk=1
      kjpt=1
      DO jvar =1,varend
         indvar = var_ord(jvar)
         IF (vardmsk(indvar).GE.1) kjpi=MAX(kjpi,var_jpi(indvar))
         IF (vardmsk(indvar).GE.2) kjpj=MAX(kjpj,var_jpj(indvar))
         IF (vardmsk(indvar).GE.3) kjpk=MAX(kjpk,var_jpk(indvar))
         IF (vardmsk(indvar).GE.4) kjpt=MAX(kjpt,var_jpt(indvar))
      ENDDO  
! 
      IF (nprint.GE.1) THEN
         WRITE(numout,*) 'Mask dimensions :'
         WRITE(numout,*) '     jpi=',kjpi
         WRITE(numout,*) '     jpj=',kjpj
         WRITE(numout,*) '     jpk=',kjpk
         WRITE(numout,*) '     jpt=',kjpt
      ENDIF
!
! --- allocation mask
      allocate ( mask(1:kjpi,1:kjpj,1:kjpk,1:kjpt), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      mask(:,:,:,:) = 0
! --- allocation var_lev
      allocate ( var_lev(1:kjpk,1:nbvar), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      var_lev(:,:) = FREAL(0.0)
!
      DO jvar=1,varend
         indvar=var_ord(jvar)
         DO jk=1,var_jpk(indvar)
            var_lev(jk,indvar)=FREAL(jk)
         ENDDO
      ENDDO
!
! -2.- Load configuration masks
! -----------------------------
!
! -2.1- read Vx configuration mask (var)
!
      IF (.NOT.novarmsk) THEN
!
         flagxyo=1
         DO jvar =1,varend
            indvar = var_ord(jvar)
            jext = varemsk(indvar)
            SELECT CASE (jext)
            CASE (1)
               CALL readmskbimg(jvar,flagxyo)
            CASE (2)
               GOTO 104
            CASE (3)
               CALL readmskcdf(jvar,flagxyo)
            CASE (4)
               CALL readmsknc(jvar,flagxyo)
            CASE DEFAULT
               GOTO 1000
            END SELECT
         ENDDO         
!
      ELSE
!
         DO jvar =1,varend
            mask(:,:,:,:) = mask(:,:,:,:) + (2**(jvar-1))
         ENDDO
!
      ENDIF
!
! -2.1- read Vy configuration mask (dta)
!
      IF (.NOT.nodtamsk) THEN
!
         flagxyo=2
         DO jdta =1,dtaend
            inddta = dta_ord(jdta)
            jext = dtaemsk(inddta)
            SELECT CASE (jext)
            CASE (1)
               CALL readmskbimg(jdta,flagxyo)
            CASE (2)
               GOTO 104
            CASE (3)
               CALL readmskcdf(jdta,flagxyo)
            CASE (4)
               CALL readmsknc(jdta,flagxyo)
            CASE DEFAULT
               GOTO 1000
            END SELECT
            jvar=1
            indvar = var_ord(jvar)
            DO WHILE ((indvar.NE.inddta).AND.(jvar.LT.varend))
               jvar=jvar+1
               indvar = var_ord(jvar)
            ENDDO  
            IF (indvar.NE.inddta) GOTO 1000
            IF (.NOT.(validmsk(jvar,jdta))) GOTO 103
         ENDDO
!
      ELSE
!
         DO jdta =1,dtaend
            mask(:,:,:,:) = mask(:,:,:,:) + (2**(jdta-1+varend))
         ENDDO
!
      ENDIF
!
      RETURN
!
! --- format definitions
!
 10   FORMAT ("|",7(2X,A4,2X,"|"))
 20   FORMAT ("|",2(I5,3X,"|"),3X,A5,"|",4(2X,I4,2X,"|"))
 30   FORMAT (A1,64("-"))
 40   FORMAT ("|",8(2X,A4,2X,"|"))
 50   FORMAT ("|",3(I5,3X,"|"),3X,A5,"|",2X,I3,"D",2X,"|",I7,1X,"|",2(E7.1,1X,"|"))
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiomsk','readmsk')
 1001 CALL printerror2(0,1001,3,'hiomsk','readmsk')
!
 101  WRITE (texterror,*) 'Bad Vx (var) mask extension : ', &
     &     extvartab(jext)
      CALL printerror2(0,101,3,'hiomsk','readmsk',comment=texterror)
 102  WRITE (texterror,*) 'Bad Vy (dta) mask extension : ', &
     &     extdtatab(jext)
      CALL printerror2(0,102,3,'hiomsk','readmsk',comment=texterror)
 103  WRITE (texterror,*) 'Bad mask configuration'
      CALL printerror2(0,103,3,'hiomsk','readmsk',comment=texterror)
 104  WRITE (texterror,*) 'Invalid mask file format: dimg or dta'
      CALL printerror2(0,104,3,'hiomsk','readmsk',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrmsk(kindxy,kjpi,kjpj,kjpk,kjpt,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Read mask file header
!  -------
!  Method : Call specific routine as a function
!  ------   of mask file format (bimg,dimg,cdf)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use liobimg
      use liodimg
      use liocdf
      use lionc
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, intent(out) :: kjpi,kjpj,kjpk,kjpt
      INTEGER, intent(in) :: kflagxyo,kindxy
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: fname
      INTEGER :: jext,sxydmsk
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/readmsk/evalhdrmsk :'
         WRITE(numout,*) '         read mask file header'
      ENDIF
!
! Set mask filename, mask filetype, mask dimensions
!
      SELECT CASE (kflagxyo)
      CASE(1)
         jext = varemsk(kindxy)
         fname = varfmsk(kindxy)
         sxydmsk = vardmsk(kindxy)
      CASE(2)
         jext = dtaemsk(kindxy)
         fname = dtafmsk(kindxy)
         sxydmsk = dtadmsk(kindxy)
      CASE(3)
         GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Read mask dimension sizes
!
      kjpi=1
      kjpj=1
      kjpk=1
      kjpt=1
      SELECT CASE (jext)
      CASE (1)
         CALL evalhdrmskbimg(fname,kjpi,kjpj,kjpk,kjpt)
      CASE (2)
         CALL evalhdrmskdimg(fname,kjpi,kjpj,kjpk,kjpt)
      CASE (3)
         CALL evalhdrmskcdf(fname,kjpi,kjpj,kjpk,kjpt)
      CASE (4)
         CALL evalhdrmsknc(fname,kjpi,kjpj,kjpk,kjpt,kindxy,kflagxyo)
      CASE DEFAULT
         GOTO 1000
      END SELECT
      IF (sxydmsk.LT.1) kjpi=1
      IF (sxydmsk.LT.2) kjpj=1
      IF (sxydmsk.LT.3) kjpk=1
      IF (sxydmsk.LT.4) kjpt=1
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiomsk','evalhdrmsk')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hiomsk
