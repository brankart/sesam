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
! ---                   EVALCONFIG.F90                            ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 98-06  ( C.E. Testut)                      ---
! --- modification : 99-05 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- modification : 03-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  evalconfig
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalconfig
!---------------------------------------------------------------------
!
!  Purpose : Set up and check SESAM internal configuration
!  -------
!  Method : Configure Vx and Vy vectors
!  ------   Compute table of pointers from Vy vector to Vx vector
!           Configure Vo vectors
!           Configure Io objects
!           Configure Vz vectors
!           Configure SESAM covariance matrices
!           Configure SESAM utilities to save computer memory
!           Configure Vx blocks
!           Configure grid connections
!           Print information about SESAM internal configuration
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : jpdbs,jpdbsend,jpmend,jpsmplend, &
     &     jprend,jpo,jpoend,jpitp,jpitpend,jpy,jpyend,jpz,jpperc, &
     &     jpx,jpxend,jpnxend,arraynx_jindxbeg, &
     &     arraynx_jpindxend,arraynx_jindvarbeg, &
     &     arraynx_jindkbeg,arraynx_jindtbeg
      use hioxyo
      use hiozon
      use hiodbs
      use utilvalid
      use mkconnect
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: fname
      INTEGER :: serie,numjr,ios,allocok
      INTEGER :: jx,jy,indvarmsk,inddtamsk,somtot,sompart
      INTEGER :: jvar,indvar,jdta,inddta,jobs,indobs,inddbs
      INTEGER :: ji,jj,jk,jt,jarg,jpisize,jpjsize,jpksize,jptsize
      INTEGER :: nbtabnx,jtabnx,jnx
      INTEGER, dimension(:), allocatable :: tabnx_jindxbeg, &
     &     tabnx_jpindxend,tabnx_jindvarbeg, &
     &     tabnx_jindkbeg,tabnx_jindtbeg
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/evalconfig :'
         WRITE(numout,*) ' setup and check SESAM internal configuration'
      ENDIF
!
! -0.- Initialization
! -------------------
! Initialize size and position of variables in SESAM objects
! Vx: var_nbr, var_ind; Vy: dta_nbr, dta_ind; Vo: obs_nbr, obs_ind
      var_nbr(:) = 1
      var_ind(:) = 1
      dta_nbr(:) = 1
      dta_ind(:) = 1
!
      DO jobs = 1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         obs_nbr(indobs,inddbs)=1
         obs_ind(indobs,inddbs)=1
      ENDDO
!
! Check if Vy mask (dta) is included in Vx mask (var)
      DO jdta =1,dtaend
         inddta = dta_ord(jdta)
         jvar=1
         indvar = var_ord(jvar)
         DO WHILE ((indvar.NE.inddta).AND.(jvar.LT.varend))
            jvar=jvar+1
            indvar = var_ord(jvar)
         ENDDO  
         IF (indvar.NE.inddta) GOTO 1000
         IF (.NOT.(validmsk(jvar,jdta))) GOTO 101
      ENDDO
! Compute number of 2D slices in Vx masks
      nbtabnx=0
      DO jvar =1,varend
         indvar = var_ord(jvar)
         SELECT CASE (var_dim(indvar))
         CASE (1,2)
            nbtabnx=nbtabnx+1
         CASE (3)
            nbtabnx=nbtabnx+var_jpk(indvar)
         CASE (4)
            nbtabnx=nbtabnx+var_jpk(indvar)*var_jpt(indvar)
         CASE DEFAULT
            GOTO 1000
         END SELECT               
      ENDDO
!
! Arrays pointing on 2D slices in Vx and Vy objects
!
! --- allocation tabnx_jindxbeg
      allocate ( tabnx_jindxbeg(1:nbtabnx) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabnx_jindxbeg(:) = 0
! --- allocation tabnx_jpindxend
      allocate ( tabnx_jpindxend(1:nbtabnx) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabnx_jpindxend(:) = 0
! --- allocation tabnx_jindvarbeg
      allocate ( tabnx_jindvarbeg(1:nbtabnx) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabnx_jindvarbeg(:) = 0
! --- allocation tabnx_jindkbeg
      allocate ( tabnx_jindkbeg(1:nbtabnx) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabnx_jindkbeg(:) = 0
! --- allocation tabnx_jindtbeg
      allocate ( tabnx_jindtbeg(1:nbtabnx) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabnx_jindtbeg(:) = 0
!
! -1.- Configuration of Vx and Vy vectors
! ---------------------------------------
! Compute size and position of variables in SESAM Vx and Vy vectors:
!        var_ind,var_nbr,dta_ind,dta_nbr
! Compute size of SESAM Vx and Vy vectors:
!        jpxend,jpyend
!
      jpxend=1
      jpyend=1
!
      jtabnx=1
      jx=1
      jy=1
      DO jvar=1,varend
         indvar=var_ord(jvar)
         indvarmsk=jvar-1
         var_ind(indvar)=jx
         IF (dta_act(indvar)) THEN
            jdta=1
            inddta = dta_ord(jdta)
            DO WHILE ((inddta.NE.indvar).AND.(jdta.LT.dtaend))
               jdta=jdta+1
               inddta = dta_ord(jdta)
            ENDDO  
            IF (inddta.NE.indvar) GOTO 1000
            inddtamsk=jdta-1+varend
            dta_ind(inddta)=jy
         ELSE
            inddtamsk=indvarmsk
         ENDIF
         IF (dta_act(indvar)) THEN
            DO jt=1,var_jpt(indvar)
            DO jk=1,var_jpk(indvar)
!
            IF (jtabnx.GT.nbtabnx) GOTO 1000
            tabnx_jindxbeg(jtabnx)=jx
!
            DO jj=1,var_jpj(indvar)
            DO ji=1,var_jpi(indvar)
               jy=jy+ABS(IBITS(mask(ji,jj,jk,jt),indvarmsk,1))* &
     &               ABS(IBITS(mask(ji,jj,jk,jt),inddtamsk,1))
               jx=jx+ABS(IBITS(mask(ji,jj,jk,jt),indvarmsk,1))
            ENDDO
            ENDDO
!
            tabnx_jpindxend(jtabnx)=jx-tabnx_jindxbeg(jtabnx)
            tabnx_jindvarbeg(jtabnx)=indvar
            tabnx_jindkbeg(jtabnx)=jk
            tabnx_jindtbeg(jtabnx)=jt
            jtabnx=jtabnx+1
!
            ENDDO
            ENDDO
         ELSE
            DO jt=1,var_jpt(indvar)
            DO jk=1,var_jpk(indvar)
!
            IF (jtabnx.GT.nbtabnx) GOTO 1000
            tabnx_jindxbeg(jtabnx)=jx
!
            DO jj=1,var_jpj(indvar)
            DO ji=1,var_jpi(indvar)
               jx=jx+ABS(IBITS(mask(ji,jj,jk,jt),indvarmsk,1))
            ENDDO
            ENDDO
!
            tabnx_jpindxend(jtabnx)=jx-tabnx_jindxbeg(jtabnx)
            tabnx_jindvarbeg(jtabnx)=indvar
            tabnx_jindkbeg(jtabnx)=jk
            tabnx_jindtbeg(jtabnx)=jt
            jtabnx=jtabnx+1
!
            ENDDO
            ENDDO
         ENDIF
         var_nbr(indvar)=jx-var_ind(indvar)
         IF (dta_act(indvar)) THEN
            dta_nbr(indvar)=jy-dta_ind(indvar)
         ENDIF
         IF (var_nbr(indvar).LT.1) GOTO 107
         IF ((dta_act(indvar)).AND.(dta_nbr(indvar).LT.1)) GOTO 108
      ENDDO
      jpxend=jx-1
      jpyend=jy-1
!
! Print information about SESAM Vx and Vy objects
!
      IF (nprint.GE.3) THEN
         WRITE (numout,*) ' '
         WRITE (numout,*) ' Indices of 2D slices in Vx object:'
         WRITE (numout,*) ' ----------------------------------'
         WRITE (numout,*) ' slice start_index val_number var_index', &
     &                    ' level_index time_index var_name'
         DO jtabnx=1,nbtabnx
            WRITE (numout,'(1x,i4,5(1x,i9),8x,a5)') jtabnx, &
     &               tabnx_jindxbeg(jtabnx), &
     &               tabnx_jpindxend(jtabnx), &
     &               tabnx_jindvarbeg(jtabnx), &
     &               tabnx_jindkbeg(jtabnx), &
     &               tabnx_jindtbeg(jtabnx), &
     &               var_nam(tabnx_jindvarbeg(jtabnx))
         ENDDO
      ENDIF
      IF (SUM(tabnx_jpindxend(:)).NE.jpxend) GOTO 1000
!
! Compute table of pointers from Vy vector to Vx vector
! (telling which element of Vx vector corresponds to
! each element of Vy vector)
!
! --- allocation tabindxtoy
      allocate ( tabindxtoy(1:jpyend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabindxtoy(:) = 0
!
      jx=1
      jy=1
      DO jvar=1,varend
         indvar=var_ord(jvar)
         indvarmsk=jvar-1
         var_ind(indvar)=jx
         IF (dta_act(indvar)) THEN
            jdta=1
            inddta = dta_ord(jdta)
            DO WHILE ((inddta.NE.indvar).AND.(jdta.LT.dtaend))
               jdta=jdta+1
               inddta = dta_ord(jdta)
            ENDDO  
            IF (inddta.NE.indvar) GOTO 1000
            inddtamsk=jdta-1+varend
            dta_ind(inddta)=jy
         ELSE
            inddtamsk=indvarmsk
         ENDIF
         IF (dta_act(indvar)) THEN
            DO jt=1,var_jpt(indvar)
            DO jk=1,var_jpk(indvar)
            DO jj=1,var_jpj(indvar)
            DO ji=1,var_jpi(indvar)
               IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) THEN
                  tabindxtoy(jy)=jx
                  jy=jy+ABS(IBITS(mask(ji,jj,jk,jt),indvarmsk,1))
               ENDIF
               jx=jx+ABS(IBITS(mask(ji,jj,jk,jt),indvarmsk,1))
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ELSE
            DO jt=1,var_jpt(indvar)
            DO jk=1,var_jpk(indvar)
            DO jj=1,var_jpj(indvar)
            DO ji=1,var_jpi(indvar)
               jx=jx+ABS(IBITS(mask(ji,jj,jk,jt),indvarmsk,1))
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF
         var_nbr(indvar)=jx-var_ind(indvar)
         IF (dta_act(indvar)) THEN
            dta_nbr(inddta)=jy-dta_ind(inddta)
         ENDIF
      ENDDO
      IF (jpxend.NE.(jx-1)) GOTO 1000
      IF (jpyend.NE.(jy-1)) GOTO 1000
!
! -2.- Configuration of Vo vectors
! --------------------------------
! Set up and check observation space
! jpoend,jpitpend,obs_nbr,obs_ind,obs_itp
!
      jpoend=0
      jpitpend=0
      IF (existobs) THEN
         IF (largconfigobs) THEN
            IF (.NOT.(validextobs(argconfigobs))) GOTO 1000
            CALL evalhdrobs (argconfigobs,jpoend,jpitpend)
            IF ((jpoend.EQ.0).AND.(nmode.NE.0)) GOTO 109
            SELECT CASE (nmode)
            CASE (0,1,4,5,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24)
! -3.0- Help  config
! -3.1- Mode intf : 
! -3.4- Mode obsv : 
! -3.5- Mode filt : 
! -3.7- Mode diff : 
! -3.8- Mode oerr : 
! -3.9- Mode oper :
! -3.10- Mode geof : 
! -3.11- Mode leof :
! -3.12- Mode beof :
! -3.14- Mode groa : 
! -3.15- Mode lroa : 
! -3.16- Mode broa : 
! -3.17- Mode greg : 
! -3.18- Mode lreg :
! -3.19- Mode breg :
! -3.20- Mode vari :
! -3.24- Mode mcmc :
            CASE (2,3,6,13)
! -3.2- Mode corr : 
! -3.3- Mode tgop : 
! -3.6- Mode adap :
! -3.13- Mode zone :
               GOTO 1000
            CASE DEFAULT
               GOTO 1000
            END SELECT
         ENDIF
      ELSE
         jpoend=1
         jpitpend=1
      ENDIF
      IF ((jpoend.EQ.0).AND.(nmode.NE.0)) GOTO 1000
!
! -3.- Configuration of Io objects
! --------------------------------
! jpdbsend,jpdbs
!
      jpdbsend=0
      IF (existdbs) THEN
         IF (largindbs) THEN
            IF (.NOT.(validextdbs(argindbs))) GOTO 1000
            CALL evalhdrdbs (argindbs,jpdbsend)
            SELECT CASE (nmode)
            CASE (0,4)
! -3.0- Help  config
! -3.4- Mode obsv :
            CASE (1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24)
! -3.1- Mode intf :
! -3.2- Mode corr :
! -3.3- Mode ---- :
! -3.5- Mode filt :
! -3.6- Mode adap :
! -3.7- Mode diff :
! -3.8- Mode oerr :
! -3.9- Mode oper :
! -3.10- Mode geof : 
! -3.11- Mode leof :
! -3.12- Mode beof :
! -3.13- Mode zone :
! -3.14- Mode groa : 
! -3.15- Mode lroa : 
! -3.16- Mode broa : 
! -3.17- Mode greg : 
! -3.18- Mode lreg :
! -3.19- Mode breg :
! -3.20- Mode vari :
! -3.24- Mode mcmc :
               GOTO 1000
            CASE DEFAULT
               GOTO 1000
            END SELECT
         ENDIF
      ENDIF
!
! -3.- Configuration of Vz objects
! --------------------------------
! zon_jpi,zon_jpj,zon_jpk,zon_jpt,jpz,jpbub
!
      zon_jpi=0
      zon_jpj=0
      zon_jpk=0
      zon_jpt=0
      jpz=0
      jpbub=0
      IF (existzon) THEN
         IF ((largconfigzon).OR.(larginzon)) THEN
            IF (largconfigzon) THEN
               IF (.NOT.(validextzon(argconfigzon))) GOTO 1000
               CALL evalhdrzon (argconfigzon,zon_jpi,zon_jpj,zon_jpk, &
     &              zon_jpt,jpbub,jpz)
            ELSEIF (larginzon) THEN
               IF (.NOT.(validextzon(arginzon))) GOTO 1000
               CALL evalhdrzon (arginzon,zon_jpi,zon_jpj,zon_jpk, &
     &              zon_jpt,jpbub,jpz)
            ENDIF
            SELECT CASE (nmode)
            CASE (0,1,3,5,9,10,11,12,13,14,15,16,17,18,19)
! -3.0- Help  config
! -3.1- Mode intf :
! -3.3- Mode tgop :
! -3.5- Mode filt :
! -3.9- Mode oper :
! -3.10- Mode geof : 
! -3.11- Mode leof :
! -3.12- Mode beof :
! -3.13- Mode zone :
! -3.14- Mode groa : 
! -3.15- Mode lroa : 
! -3.16- Mode broa : 
! -3.17- Mode greg : 
! -3.18- Mode lreg :
! -3.19- Mode breg :
            CASE (2,4,6,7,8,20,21,23,24)
! -3.2- Mode corr :
! -3.4- Mode obsv :
! -3.6- Mode adap :
! -3.7- Mode diff :
! -3.8- Mode oerr :
! -3.10- Mode vari :
! -3.10- Mode mcmc :
               GOTO 1000
            CASE DEFAULT
               GOTO 1000
            END SELECT
      ENDIF
      ELSE
         zon_jpi=1
         zon_jpj=1
         zon_jpk=1
         zon_jpt=1
         jpz=1
         jpbub=1
      ENDIF
      IF (jpz.EQ.0) GOTO 1000
      IF (jpbub.EQ.0) GOTO 1000
      IF (zon_jpi.EQ.0) GOTO 1000
      IF (zon_jpj.EQ.0) GOTO 1000
      IF (zon_jpk.EQ.0) GOTO 1000
      IF (zon_jpt.EQ.0) GOTO 1000
!
! -4.- Print information about SESAM objects configuration
! --------------------------------------------------------
!     
      IF (nprint.GE.1) THEN
         WRITE (numout,*)
         WRITE (numout,*) ' Vx configuration (var)'
         WRITE (numout,*) ' ----------------------'
         WRITE (numout,10)  &
     &        'Variables',' var_ind ',' var_nbr ',' ratio(%)'
         somtot=0
         DO jvar=1,varend
            indvar=var_ord(jvar)
            sompart=var_jpi(indvar)*var_jpj(indvar)* &
     &           var_jpk(indvar)*var_jpt(indvar)
            somtot=somtot+sompart
            WRITE (numout,11) var_nam(indvar), &
     &           var_ind(indvar),var_nbr(indvar), &
     &       INT(FREAL(var_nbr(indvar)*100)/ FREAL(sompart))
         ENDDO
         WRITE (numout,*) ' Vx size (jpxend) = ',jpxend
         WRITE (numout,*) ' valid points ratio (%) = ', &
     &        INT(FREAL(jpxend*100)/ FREAL(somtot))
!
         WRITE (numout,*)
         WRITE (numout,*) ' Vy configuration (dta)'
         WRITE (numout,*) ' ----------------------'
         WRITE (numout,10)  &
     &        'Variables',' dta_ind ',' dta_nbr ',' ratio(%)'
         somtot=0
         DO jdta=1,dtaend
            inddta=dta_ord(jdta)
            sompart=dta_jpi(inddta)*dta_jpj(inddta)* &
     &           dta_jpk(inddta)*dta_jpt(inddta)
            somtot=somtot+sompart
            WRITE (numout,11) dta_nam(inddta), &
     &              dta_ind(inddta),dta_nbr(inddta), &
     &       INT(FREAL(dta_nbr(inddta)*100)/ FREAL(sompart))
         ENDDO
         WRITE (numout,*) ' Vy size (jpyend) = ',jpyend
         WRITE (numout,*) ' valid points ratio (%) = ', &
     &        INT(FREAL(jpyend*100)/ FREAL(somtot))
         IF (existobs) THEN
            WRITE (numout,*)
            WRITE (numout,*) ' Vo configuration (obs)'
            WRITE (numout,*) ' ----------------------'
            WRITE (numout,12) 'Variables',' obs_ind ',' obs_nbr ', &
     &           ' jpitploc',' ratio(%)'
            DO jobs=1,obsend
               indobs=obs_ord(jobs)
               inddbs=obsnord(jobs)
               WRITE (numout,13) obs_nam(indobs,inddbs), &
     &              obs_ind(indobs,inddbs),obs_nbr(indobs,inddbs), &
     &              obs_itp(indobs,inddbs), &
     &              INT(FREAL(obs_nbr(indobs,inddbs)*100)/FREAL(jpoend))
            ENDDO
         WRITE (numout,*) ' Vo size (jpoend) = ',jpoend
         WRITE (numout,*) ' Nbr of interpolation points (jpitpend) = ',jpitpend
         ENDIF
         IF (existdbs) THEN
            WRITE (numout,*)
            WRITE (numout,*) ' Io configuration (dbs)'
            WRITE (numout,*) ' ----------------------'
            WRITE (numout,18) '  obs_nam  ',' jpdbsend  '
            WRITE (numout,19) argaffectobs(1:lenv(argaffectobs)),jpdbsend
         ENDIF
         IF (existzon) THEN
            WRITE (numout,*)
            WRITE (numout,*) ' Vz configuration (zon)'
            WRITE (numout,*) ' ----------------------'
            WRITE (numout,20) ' zon_jpi ',' zon_jpj ',' zon_jpk ', &
     &           ' zon_jpt ','   jpz   ','  jpbub  '
            WRITE (numout,21) zon_jpi,zon_jpj,zon_jpk,zon_jpt, &
     &           jpz,jpbub
         ENDIF
      ENDIF
!
! -5.- Configuration of SESAM covariance matrices
! -----------------------------------------------
! set 'jprend' : covariance reduced rank
! check if rank of covariance matrices are all coherent
!
      jprend=0
      IF ((existbas).OR.(existybas).OR.(existobas).OR.(existzbas)) THEN
         serie=0
         numjr=1
         SELECT CASE (nmode)
         CASE (0,1,2,4,8,10,11,12,14,15,16,17,18,19,20)
! -5.0- Help Config
! -5.1- Mode intf :
! -5.2- Mode corr :
! -5.4- Mode obsv :
! -5.8- Mode oerr :
! -5.10- Mode geof : 
! -5.11- Mode leof :
! -5.12- Mode beof :
! -5.14- Mode groa : 
! -5.15- Mode lroa : 
! -5.16- Mode broa : 
! -5.17- Mode greg : 
! -5.18- Mode lreg :
! -5.19- Mode breg :
! -5.20- Mode vari :
            IF (larginxbas) THEN
               IF (.NOT.(validextvarbas(arginxbas))) GOTO 1000
               CALL fildirbas (fname,arginxbas,jprend,numjr,serie)
            ELSEIF (larginybas) THEN
               IF ((.NOT.validextvarbas(arginybas)) &
     &              .AND.(.NOT.validextdtabas(arginybas))) GOTO 1000
               CALL fildirbas (fname,arginybas,jprend,numjr,serie)
            ELSEIF (larginobas) THEN
               IF ((.NOT.validextvarbas(arginobas)) &
     &              .AND.(.NOT.validextdtabas(arginobas)) &
     &              .AND.(.NOT.validextobsbas(arginobas))) GOTO 1000
               CALL fildirbas (fname,arginobas,jprend,numjr,serie)
            ELSEIF (larginzbas) THEN
               IF (.NOT.(validextzonbas(arginzbas))) GOTO 1000
               CALL fildirbas (fname,arginzbas,jprend,numjr,serie)
            ELSE
               GOTO 1000
            ENDIF
         CASE (3)
! -5.3- Mode tgop :
            IF (larginxbas) THEN
               IF (.NOT.(validextvarbas(arginxbas))) GOTO 1000
               CALL fildirbas (fname,arginxbas,jprend,numjr,serie)
            ELSE
               GOTO 1000
            ENDIF
            IF (largoutxbas) THEN
               IF (.NOT.(validextvarbas(argoutxbas))) GOTO 1000
               CALL fildirbas (fname,argoutxbas,jpsmplend,numjr,serie)
            ENDIF
         CASE (5,6,7,9)
! -5.5- Mode filt :
! -5.6- Mode adap :
! -5.7- Mode diff :
! -5.9- Mode oper :
            GOTO 103
         CASE (13)
! -5.13- Mode zone :
            IF (largoutzbas) THEN
               IF (.NOT.(validextzonbas(argoutzbas))) GOTO 1000
               CALL fildirbas (fname,argoutzbas,jprend,numjr,serie)
            ELSE
               GOTO 1000
            ENDIF
! -5.21- Mode anam :
         CASE (21)
            IF (larginxbas) THEN
               IF (.NOT.(validextvarbas(arginxbas))) GOTO 1000
               CALL fildirbas (fname,arginxbas,jprend,numjr,serie)
            ELSEIF (larginybas) THEN
               IF ((.NOT.validextvarbas(arginybas)) &
     &              .AND.(.NOT.validextdtabas(arginybas))) GOTO 1000
               CALL fildirbas (fname,arginybas,jprend,numjr,serie)
            ELSEIF (larginobas) THEN
               IF ((.NOT.validextvarbas(arginobas)) &
     &              .AND.(.NOT.validextdtabas(arginobas)) &
     &              .AND.(.NOT.validextobsbas(arginobas))) GOTO 1000
               CALL fildirbas (fname,arginobas,jprend,numjr,serie)
            ELSE
               jprend=1
            ENDIF
!
            IF (larginxbasref) THEN
               IF (.NOT.(validextvarbas(arginxbasref))) GOTO 1000
               CALL fildirbas (fname,arginxbasref,jpperc,numjr,serie)
            ELSEIF (larginybasref) THEN
               IF ((.NOT.validextvarbas(arginybasref)) &
     &              .AND.(.NOT.validextdtabas(arginybasref))) GOTO 1000
               CALL fildirbas (fname,arginybasref,jpperc,numjr,serie)
            ELSEIF (larginobasref) THEN
               IF ((.NOT.validextvarbas(arginobasref)) &
     &              .AND.(.NOT.validextdtabas(arginobasref)) &
     &              .AND.(.NOT.validextobsbas(arginobasref))) GOTO 1000
               CALL fildirbas (fname,arginobasref,jpperc,numjr,serie)
            ELSEIF (largoutxbasref) THEN
               IF (.NOT.(validextvarbas(argoutxbasref))) GOTO 1000
               CALL fildirbas (fname,argoutxbasref,jpperc,numjr,serie)
            ELSEIF (largoutybasref) THEN
               IF (.NOT.validextdtabas(argoutybasref)) GOTO 1000
               CALL fildirbas (fname,argoutybasref,jpperc,numjr,serie)
            ELSEIF (largoutobasref) THEN
               IF (.NOT.validextobsbas(argoutobasref)) GOTO 1000
               CALL fildirbas (fname,argoutobasref,jpperc,numjr,serie)
            ELSEIF (largoutybas) THEN
               IF (.NOT.(validextdtabas(argoutybas))) GOTO 1000
               CALL fildirbas (fname,argoutybas,jpsmplend,numjr,serie)
            ELSEIF (largoutobas) THEN
               IF (.NOT.(validextobsbas(argoutobas))) GOTO 1000
               CALL fildirbas (fname,argoutobas,jpsmplend,numjr,serie)
            ELSE
               GOTO 1000
            ENDIF
! -5.22- Mode scor :
         CASE (22)
            IF (larginxbas) THEN
               IF (.NOT.(validextvarbas(arginxbas))) GOTO 1000
               CALL fildirbas (fname,arginxbas,jprend,numjr,serie)
            ELSEIF (larginybas) THEN
               IF ((.NOT.validextvarbas(arginybas)) &
     &              .AND.(.NOT.validextdtabas(arginybas))) GOTO 1000
               CALL fildirbas (fname,arginybas,jprend,numjr,serie)
            ELSEIF (larginobas) THEN
               IF ((.NOT.validextvarbas(arginobas)) &
     &              .AND.(.NOT.validextdtabas(arginobas)) &
     &              .AND.(.NOT.validextobsbas(arginobas))) GOTO 1000
               CALL fildirbas (fname,arginobas,jprend,numjr,serie)
            ELSE
               jprend=1
            ENDIF
! -5.23- Mode spct :
         CASE (23)
            IF (larginxbas) THEN
               IF (.NOT.(validextvarbas(arginxbas))) GOTO 1000
               CALL fildirbas (fname,arginxbas,jprend,numjr,serie)
            ELSEIF (larginybas) THEN
               IF ((.NOT.validextvarbas(arginybas)) &
     &              .AND.(.NOT.validextdtabas(arginybas))) GOTO 1000
               CALL fildirbas (fname,arginybas,jprend,numjr,serie)
            ELSEIF (larginobas) THEN
               IF ((.NOT.validextvarbas(arginobas)) &
     &              .AND.(.NOT.validextdtabas(arginobas)) &
     &              .AND.(.NOT.validextobsbas(arginobas))) GOTO 1000
               CALL fildirbas (fname,arginobas,jprend,numjr,serie)
            ELSE
               jprend=1
            ENDIF
! -5.24- Mode mcmc :
         CASE (24)
            IF (larginxbas) THEN
               IF (.NOT.(validextvarbas(arginxbas))) GOTO 1000
               CALL fildirbas (fname,arginxbas,jprend,numjr,serie)
            ELSEIF (larginybas) THEN
               IF ((.NOT.validextvarbas(arginybas)) &
     &              .AND.(.NOT.validextdtabas(arginybas))) GOTO 1000
               CALL fildirbas (fname,arginybas,jprend,numjr,serie)
            ELSEIF (larginobas) THEN
               IF ((.NOT.validextvarbas(arginobas)) &
     &              .AND.(.NOT.validextdtabas(arginobas)) &
     &              .AND.(.NOT.validextobsbas(arginobas))) GOTO 1000
               CALL fildirbas (fname,arginobas,jprend,numjr,serie)
            ELSE
               jprend=1
            ENDIF

            IF (larganamorphosis) THEN
               IF ((.NOT.validextvarbas(arganamorphosis)) &
     &              .AND.(.NOT.validextdtabas(arganamorphosis)) &
     &              .AND.(.NOT.validextobsbas(arganamorphosis))) GOTO 1000
               CALL fildirbas (fname,arganamorphosis,jpperc,numjr,serie)
            ELSE
               jpperc=1
            ENDIF

            IF (largoutxbas) THEN
               IF (.NOT.(validextvarbas(argoutxbas))) GOTO 1000
               CALL fildirbas (fname,argoutxbas,jpsmplend,numjr,serie)
            ELSEIF (largoutybas) THEN
               IF (.NOT.(validextdtabas(argoutybas))) GOTO 1000
               CALL fildirbas (fname,argoutybas,jpsmplend,numjr,serie)
            ELSEIF (largoutobas) THEN
               IF (.NOT.(validextobsbas(argoutobas))) GOTO 1000
               CALL fildirbas (fname,argoutobas,jpsmplend,numjr,serie)
            ELSE
               GOTO 1000
            ENDIF

            IF (largiterate) THEN
               READ(argiterate,*,IOSTAT=ios) maxiter
               IF (ios.NE.0) GOTO 111
            ENDIF
!
         CASE DEFAULT
            GOTO 1000
         END SELECT
      ELSE
         jprend=1
      ENDIF
      IF (jprend.EQ.0) GOTO 1000
!
! -6.- Configuration of SESAM covariance matrices
! -----------------------------------------------
! set 'jpmend' : number of constraints
! check if rank of covariance matrices are all coherent
!
      jpmend=0
      IF ((existcstr)) THEN
         serie=0
         numjr=1
         SELECT CASE (nmode)
         CASE (3)
            IF (largincstr) THEN
               IF (.NOT.(validextvarbas(argincstr))) GOTO 1000
               CALL fildirbas (fname,argincstr,jpmend,numjr,serie)
            ELSE
               GOTO 1000
            ENDIF
         CASE DEFAULT
            GOTO 1000
         END SELECT
      ELSE
         jpmend=1
      ENDIF
      IF (jpmend.LT.0) GOTO 1000
!
! -7.- Configuration of SESAM utilities to save computer memory
! -------------------------------------------------------------
! optional switches: fixjpx, fixjpz, fixjpu
!
      IF (largfixjpx) THEN
         READ(argfixjpx,'(I8)',IOSTAT=ios) jpfixjpx
         IF (ios.NE.0) THEN
            PRINT *,'WARNING: invalid fixjpx argument'
            PRINT *,'         a division in two segments is applied'
            jpfixjpx = jpxend-1
         ENDIF
         IF (jpfixjpx.GE.jpxend) THEN
            PRINT *,'WARNING: fixjpx larger than state vector'
            PRINT *,'         a division in two segments is applied'
            jpfixjpx = jpxend-1
         ENDIF
! Compute the smallest jpfixjpx with the same number of segments
! and make the number of segments a multiple of the number of processors
         jpfixjpx = ( jpxend - 1 ) / jpfixjpx + 1
         jpfixjpx = ( jpfixjpx - 1 ) / jpproc + 1
         jpfixjpx = jpfixjpx * jpproc
         jpfixjpx = ( jpxend - 1 ) / jpfixjpx + 1
      ELSEIF (jpproc.GT.1) THEN
         jpfixjpx = ( jpxend - 1 ) / jpproc + 1
      ELSE
         jpfixjpx = jpxend
      ENDIF
!
      IF (largfixjpz) THEN
         READ(argfixjpz,'(I8)',IOSTAT=ios) jpfixjpz
         IF (ios.NE.0) THEN
            PRINT *,'WARNING: invalid fixjpz argument'
            jpfixjpz = jpz
            largfixjpz = .FALSE.
         ENDIF
         IF (existzon) THEN
            jpfixjpz = MIN(jpfixjpz,jpz)
         ELSE
            PRINT *,'WARNING: ignored fixjpz switch'
            largfixjpz = .FALSE.
         ENDIF
      ELSE
         jpfixjpz = jpz
      ENDIF
!
      IF (largfixjpu) THEN
         READ(argfixjpu,'(I8)',IOSTAT=ios) jpfixjpu
         IF (ios.NE.0) THEN
            PRINT *,'WARNING: invalid fixjpu argument'
            jpfixjpu = jpy
            largfixjpu = .FALSE.
         ENDIF
      ELSE
         jpfixjpu = jpy
      ENDIF
!
! -8.- Vx block configuration
! ---------------------------
! Set size of Vx block in memory (jpx)
! and number of blocks (jpnxend)
! arraynx_... ,nallmem :
!
! nallmem = 2 : Vx vectors are loaded in memory block by block
! nallmem = 3 : Vx vectors are fully loaded in memory 
!
      IF (largfixjpx.OR.(jpproc.GT.1)) THEN
         nallmem=2
         jpx = jpfixjpx
      ELSE
         nallmem=3
         jpx = jpxend
      ENDIF       
!
! initialize arrays defining blocks of Vx vectors
      SELECT CASE (nallmem)
      CASE (2)
         IF (jpx.EQ.jpxend) GOTO 1000
! set number of blocks
         jpnxend = INT((jpxend-1)/jpx)+1
! allocate arrays
         allocate ( arraynx_jindxbeg(1:jpnxend) , stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         arraynx_jindxbeg(:) = 0
         allocate ( arraynx_jpindxend(1:jpnxend) , stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         arraynx_jpindxend(:) = 0
         allocate ( arraynx_jindvarbeg(1:jpnxend) , stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         arraynx_jindvarbeg(:) = 0
         allocate ( arraynx_jindkbeg(1:jpnxend) , stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         arraynx_jindkbeg(:) = 0
         allocate ( arraynx_jindtbeg(1:jpnxend) , stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         arraynx_jindtbeg(:) = 0
!
         arraynx_jindxbeg(1)=1
         DO jnx=2,jpnxend
            arraynx_jindxbeg(jnx)=arraynx_jindxbeg(jnx-1)+jpx
         ENDDO
!
         DO jnx=1,jpnxend-1
            arraynx_jpindxend(jnx)=jpx
         ENDDO
         arraynx_jpindxend(jpnxend)=jpxend-(jpnxend-1)*jpx
!
         jtabnx=1
         DO jnx=1,jpnxend
            DO WHILE ((jtabnx.LE.nbtabnx).AND.(.NOT. &
     &           ((arraynx_jindxbeg(jnx).GE.tabnx_jindxbeg(jtabnx)) &
     &           .AND.(arraynx_jindxbeg(jnx).LE. &
     &         (tabnx_jindxbeg(jtabnx))+tabnx_jpindxend(jtabnx)-1))))
               jtabnx=jtabnx+1
            ENDDO
            arraynx_jindvarbeg(jnx)=tabnx_jindvarbeg(jtabnx)
            arraynx_jindkbeg(jnx)=tabnx_jindkbeg(jtabnx)
            arraynx_jindtbeg(jnx)=tabnx_jindtbeg(jtabnx)
         ENDDO
! Print information about Vx blocks structure
         IF (nprint.GE.2) THEN
            WRITE (numout,*)
            WRITE (numout,*) ' Vx object block structure'
            WRITE (numout,*) ' -------------------------'
            WRITE (numout,*) ' block Vx_start Vx_end var k t'
            DO jnx=1,jpnxend
               WRITE (numout,*) ' ',jnx,'  ',arraynx_jindxbeg(jnx), &
     &              '  ',arraynx_jpindxend(jnx), &
     &              '  ',arraynx_jindvarbeg(jnx), &
     &              '  ',arraynx_jindkbeg(jnx), &
     &              '  ',arraynx_jindtbeg(jnx)
            ENDDO
         ENDIF
         IF (SUM(arraynx_jpindxend(:)).NE.jpxend) GOTO 1000
      CASE (3)
         IF (jpx.NE.jpxend) GOTO 1000
         jpnxend = 1
! allocate arrays
         allocate ( arraynx_jindxbeg(1:jpnxend) , stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         arraynx_jindxbeg(1) = 1
         allocate ( arraynx_jpindxend(1:jpnxend) , stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         arraynx_jpindxend(1) = jpx
         allocate ( arraynx_jindvarbeg(1:jpnxend) , stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         arraynx_jindvarbeg(1) = var_ord(1)
         allocate ( arraynx_jindkbeg(1:jpnxend) , stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         arraynx_jindkbeg(:) = 1
         allocate ( arraynx_jindtbeg(1:jpnxend) , stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         arraynx_jindtbeg(:) = 1
         IF (nprint.GE.3) THEN
            WRITE (numout,*)
            WRITE (numout,*) ' Vx object block structure'
            WRITE (numout,*) ' -------------------------'
            WRITE (numout,*) ' block Vx_start Vx_end var k t'
            DO jnx=1,jpnxend
               WRITE (numout,'(6(1x,i9))') jnx, &
     &                    arraynx_jindxbeg(jnx), &
     &                    arraynx_jpindxend(jnx), &
     &                    arraynx_jindvarbeg(jnx), &
     &                    arraynx_jindkbeg(jnx), &
     &                    arraynx_jindtbeg(jnx)
            ENDDO
         ENDIF
         IF (SUM(arraynx_jpindxend(:)).NE.jpxend) GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      jpy=jpyend
      jpo=jpoend
      jpitp=jpitpend
      jpdbs=jpdbsend
      limjpnxyo(1:3) = (/ jpnxend,1,1 /)
!
! Print information about Vx blocks size
      IF (nprint.GE.1) THEN
         WRITE (numout,*)
         WRITE (numout,*)  ' Vx object block configuration'
         WRITE (numout,*)  ' -----------------------------'
         WRITE (numout,*)  ' size of Vx blocks in memory = ',jpx
         WRITE (numout,*)  ' number of blocks = ',jpnxend
      ENDIF
!
! -9.0- Initialize arrays containing grid connections
! ---------------------------------------------------
!
      IF (largconnect) CALL connect_init
!
      RETURN
!
! --- format definitions
!
 10   FORMAT(2X,A9,3(1X,"|",2X,A9,1X))
 11   FORMAT(5X,A3,4X,"|",1X,3(I11,3X))
 12   FORMAT(2X,A9,4(1X,"|",2X,A9,1X))
 13   FORMAT(5X,A3,4X,"|",1X,4(I11,3X))
 14   FORMAT(A5,8X,A)
 15   FORMAT(A1,A3,A1,3X,A,1X,A1,A,A1)
 16   FORMAT(2X,A9,3(1X,"|",2X,A9,1X))
 17   FORMAT(5X,A3,4X,"|",1X,3(I11,3X))
 18   FORMAT(1X,A11,3(1X,"|",1X,A11))
 19   FORMAT(4X,A5,4X,"|",1X,3(I11,3X))
 20   FORMAT(1X,A11,5(1X,"|",1X,A11))
 21   FORMAT(2X,I9,5(1X,"|",2X,I9,1X))
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'evalconfig','evalconfig')
 1001 CALL printerror2(0,1001,3,'evalconfig','evalconfig')
!
 101  WRITE (texterror,*) 'Error Vy mask (dta) not included', &
     &                    ' in Vx mask (var)'
      CALL printerror2(0,101,3,'evalconfig','evalconfig', &
     &     comment=texterror)
 103  WRITE (texterror,*) 'Cannot use covariance matrices', &
     &                    ' with this module'
      CALL printerror2(0,103,1,'evalconfig','evalconfig', &
     &     comment=texterror)
 107  WRITE (texterror,*) 'Bad Vx mask (var): variable ', &
     &     var_nam(indvar)(1:lenv(var_nam(indvar))), &
     &     ', size= ',var_nbr(indvar)
      CALL printerror2(0,107,3,'evalconfig','evalconfig', &
     &     comment=texterror)
 108  WRITE (texterror,*) 'Bad Vy mask (dta): variable ', &
     &     dta_nam(indvar)(1:lenv(dta_nam(indvar))), &
     &     ', size= ',dta_nbr(indvar)
      CALL printerror2(0,108,3,'evalconfig','evalconfig', &
     &     comment=texterror)
 109  WRITE (texterror,*) 'Error: empty observation vector', &
     &     ' (jpoend=',jpoend,',jpitpend=',jpitpend,')'
      CALL printerror2(0,109,3,'evalconfig','evalconfig', &
     &     comment=texterror)
 110  WRITE (texterror,*) 'The two input ensembles', &
     &                    ' must have the same size'
      CALL printerror2(0,110,3,'evalconfig','evalconfig', &
     &     comment=texterror)
 111  WRITE (texterror,*) 'Bad number of iterations'
      CALL printerror2(0,111,3,'evalconfig','evalconfig', &
     &     comment=texterror)
!
      END
