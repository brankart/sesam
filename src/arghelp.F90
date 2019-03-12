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
! ---                  ARGHELP.F90                                ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-05 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        --- 
! --- modification : 03-02 (J.M. Brankart)                      --- 
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE readarghelp : Display online help
! --- FUNCTION ETAT : Convert 'TRUE/FALSE' into 'ON/OFF'
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE arghelp
      use mod_main
      use utilarg
      use utilvalid
      IMPLICIT NONE
      PRIVATE

      PUBLIC readarghelp

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readarghelp (numhelp)
!---------------------------------------------------------------------
!
!  Purpose : Analyse help commandline and display
!  -------   SESAM online help
!
!  Input : help type
!  -----
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpoend,jpitpend,jpdbsend, &
     &     jpx,jpxend,jpnxend,jpyend,jprend,jpz
      use hiomsk
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, intent(in) :: numhelp
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=swilg+1) :: swi
      CHARACTER(len=bgword) :: wrd,progname
      INTEGER :: jmod,jaction,jorder,jarg,jporder,nbprintorder
      INTEGER :: jactiondeb,jactionfin,jactionone,pos,posfin,lengthmax
      INTEGER :: jcas,jloop,indmod,jtype,jtypebeg,jtypeend,jtypenum
      INTEGER :: jextbas,jextdbs,jextdta,jextobs,jextvar,jextzon
      INTEGER :: jvar,indvar,jdta,inddta,jobs,indobs,inddbs,jnbswiopt
      CHARACTER(len=255) :: text1,text11,text2,text21
      CHARACTER(len=bgword) :: textmodel,textdate,textjpr,textfvar
      LOGICAL :: etatbool
      LOGICAL, dimension(1:nbarg) :: swibool1,swibool2
      LOGICAL, dimension(1:nbmod) :: modbool1
      LOGICAL, dimension(0:nbextbas) :: extbasbool1
      LOGICAL, dimension(0:nbextdbs) :: extdbsbool1
      LOGICAL, dimension(0:nbextdta) :: extdtabool1
      LOGICAL, dimension(0:nbextobs) :: extobsbool1
      LOGICAL, dimension(0:nbextvar) :: extvarbool1
      LOGICAL, dimension(0:nbextzon) :: extzonbool1
      INTEGER ierror
#if defined GETARG
      EXTERNAL getarg
#endif
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/readarg/readarghelp :'
         WRITE(numout,*) '        determine SESAM help from commandline'
      ENDIF
!
      PRINT *, ' '
      SELECT CASE (numhelp)
      CASE (:-(nbhelp+1),0,(nbmod+2):)
!
! These cases should be impossible
         DO jarg=2,2
            PRINT *, ' ',switab(jarg)(1:lenv(switab(jarg))), &
     &           ' <',swihelp(jarg)(1:lenv(swihelp(jarg))), &
     &           '>  ==> possible help switches'
         ENDDO
         PRINT *, ' '
         PRINT *, 'play again'
!
      CASE (-1,1:(nbmod+1))
!
! -1.- Display general help information
! -------------------------------------
! -help <all,ext,mode> -action <num>
!
! Build table (swibool1) of switches
! that are possible in this help context
         DO jarg=1,nbarg
            swibool1(jarg)= .FALSE.
            swibool2(jarg)= .FALSE.
         ENDDO
! ==> -help
         swibool1(2) = swibool(2)
         larghelp=.FALSE.
! ==> -action (optional)
         swibool1(80) = swibool(80)

!
! Read user commandline (again!), check switch and argument validity
         argloop1: DO jloop = 1, nbargmax
!
#if defined GETARG
            CALL getarg ((jloop*2-1),swi)
            CALL getarg ((jloop*2),wrd)
#else 
            CALL get_command_argument ((jloop*2-1),swi)
            CALL get_command_argument ((jloop*2),wrd)
#endif
! ==> End of loop
            IF (swi.EQ.'  ') THEN 
               EXIT argloop1
            ENDIF
! ==> Write switch and argument
            IF (nprint.GE.1) THEN
               WRITE(numout,*) ' switch/argument: ', swi(1:lenv(swi)), &
     &                                          ' ', wrd(1:lenv(wrd))
            ENDIF
! ==> Check switch validity
            IF (swi((swilg+1):(swilg+1)).NE.' ') EXIT argloop1
            IF ((.NOT.(validswi(swi))).AND. &
     &           (swi(1:lenv(swi)).NE.switab(2)(1:2))) THEN 
               EXIT argloop1
            ENDIF
            IF ((wrd.EQ.'  ').AND. &
     &           (swi(1:lenv(swi)).NE.switab(2)(1:2)).AND. &
     &           (swi.NE.switab(2))) THEN 
               EXIT argloop1
            ENDIF
! ==> Check if all switches are valid in help context
            jcas=0
            DO jarg=1,nbarg
               IF (swibool(jarg)) THEN
                  IF ((swi.EQ.switab(jarg)).OR. &
     &                 ((jarg.EQ.2).AND. &
     &                 (swi(1:lenv(swi)).EQ.switab(2)(1:2)))) THEN
                     IF (swibool1(jarg)) THEN
                        jcas=jarg
                     ELSE
                        EXIT argloop1
                     ENDIF
                  ENDIF                
               ENDIF                
            ENDDO
            IF (jcas.EQ.0) EXIT argloop1
            IF (jcas.NE.2) THEN
!  ==> Check argument validity
               CALL affectarg(wrd,jmod,jcas,swibool1,modbool1, &
     &           extbasbool1,extdbsbool1,extdtabool1, &
     &           extobsbool1,extvarbool1,extzonbool1)
! ==> Set swibool2 to .TRUE. if switch is present in the commandline
               swibool2(jcas)= .TRUE.
            ENDIF
!     
         END DO argloop1
!
! Display help on general optional switches
         IF (.NOT.largaction) THEN
            DO jarg=2,2
               PRINT *, ' ',switab(jarg)(1:lenv(switab(jarg))), &
     &           ' <',swihelp(jarg)(1:lenv(swihelp(jarg))), &
     &           '>  ==> possible help switches'
            ENDDO
            PRINT *, ' '
            PRINT 1050,'state','general optional switches :'
            DO jarg=1,nbarg
               IF (swiopt(jarg).EQ.1) then
                  PRINT 1060,'(',ETAT(swibool(jarg)),')',                   &
     &              switab(jarg)(1:lenv(switab(jarg))), &
     &              '<',swihelp(jarg)(1:lenv(swihelp(jarg))),'>'
               ENDIF
            ENDDO
            PRINT *, ' '
            SELECT CASE (numhelp)
            CASE (-1,nbmod+1)
               PRINT *, ' listing of SESAM modules :'
            CASE (1:nbmod)
               PRINT *, ' list of action for this module :'
            CASE DEFAULT
               PRINT *, ' program error :'
            END SELECT
         ENDIF
!         
! Display information on SESAM modules
         jtypebeg=1
         jtypeend=3
         DO jtype=jtypebeg,jtypeend
            jtypenum=1
            IF ((numhelp.EQ.-1).OR.(numhelp.EQ.(nbmod+1))) THEN
               PRINT *, ' '
               SELECT CASE (jtype)
               CASE (1)
                  PRINT *, ' analysis modes :'
               CASE (2)
                  PRINT *, ' diagnostic modes :'
               CASE (3)
                  PRINT *, ' utility modes :'
               CASE DEFAULT
                  PRINT *, ' program error :'
               END SELECT
            ENDIF
            DO jmod=1,nbmod
               IF ((modtype(jmod).EQ.jtype).AND.modbool(jmod).AND. &
     &             (((numhelp.EQ.jmod).OR. &
     &             (numhelp.EQ.(nbmod+1)).OR.(numhelp.EQ.-1)) &
     &             .AND.(tabvalorder(jmod,1).NE.0))) THEN
                 IF (numhelp.EQ.jmod) THEN
                    PRINT '(a,a,a,a)', ' -mode ', &
     &                modtab(jmod)(1:lenv(modtab(jmod))),' :' &
     &                   ,modhelp(jmod)(1:lenv(modhelp(jmod)))
                 ELSE
                    PRINT '(I2,a,a,a,a)', jtypenum,'  =>: -mode ', &
     &                modtab(jmod)(1:lenv(modtab(jmod))),' :' &
     &                   ,modhelp(jmod)(1:lenv(modhelp(jmod)))
                    jtypenum=1+jtypenum
                 ENDIF
               ENDIF
!
! Display information on specified SESAM module
               IF ((modtype(jmod).EQ.jtype).AND.modbool(jmod).AND. &
     &             (((numhelp.EQ.jmod).OR.(numhelp.EQ.-1)) &
     &             .AND.(tabvalorder(jmod,1).NE.0))) THEN
                 IF (largaction) PRINT *,' '
                 PRINT 1030,'action','state','rule','list of required switches'
                 IF (largaction) THEN
                    READ(argaction,*) jactiondeb
                    jactionfin=jactiondeb
                    nbprintorder=1
                 ELSE
                    jactiondeb=1
                    jactionfin=nbaction
                    nbprintorder=3
                 ENDIF               
!
! Display information on all SESAM actions corresponding to this SESAM module
                 DO jaction=jactiondeb,jactionfin
                   IF  (tabvalorder(jmod,jaction).NE.0) THEN
                     etatbool=.TRUE.
                     DO jporder=1,((((tabvalorder(jmod,jaction)/10)-1)/nbprintorder)+1)
                     text1=' '
                       DO jorder=1+(jporder-1)*nbprintorder, &
     &                       MIN((tabvalorder(jmod,jaction)/10),jporder*nbprintorder)
                         text11=text1
                         jarg=taborder(jmod,jaction,jorder)
                         WRITE (text1,'(A,1X,A,1X,A1,A,A1)')  &
     &                        text11(1:lenv(text11)), &
     &                        switab(jarg)(1:lenv(switab(jarg))),'<', &
     &                        swihelp(jarg)(1:lenv(swihelp(jarg))),'>'    
                         etatbool=(etatbool.OR.swibool(jarg))
                       ENDDO
                       etatbool=(etatbool.AND. &
     &                    (MOD(tabvalorder(jmod,jaction),10).NE.0))
                       IF (etatbool) THEN
                         IF (jporder.EQ.1) THEN
                            PRINT 1040,jaction,'(', &
     &                          ETAT(etatbool) &
     &                          ,')',(MOD(tabvalorder(jmod,jaction),10)) &
     &                          ,text1(1:lenv(text1)) 
                         ELSE                 
                            PRINT 1045,text1(1:lenv(text1))
                         ENDIF
                       ENDIF
                     ENDDO
                   ENDIF
                 ENDDO
!
! Display information on optional switches particular to this SESAM module
                 LOOP1 : DO jnbswiopt=1,nbswiopt
                    IF (tabnumswiopt(jmod,jnbswiopt).NE.0) then
                       IF (largaction) PRINT *,' '
                       PRINT 1050,'state','particular optional switches :'
                       EXIT LOOP1
                    ENDIF
                 ENDDO LOOP1
                 DO jnbswiopt=1,nbswiopt
                    IF (tabnumswiopt(jmod,jnbswiopt).NE.0) then
                       PRINT 1060,'(', &
     &                    ETAT(swibool(tabnumswiopt(jmod,jnbswiopt))),')',                   &
     &                    switab(tabnumswiopt(jmod,jnbswiopt)) &
     &                    (1:lenv(switab(tabnumswiopt(jmod,jnbswiopt)))), &
     &                    '<',swihelp(tabnumswiopt(jmod,jnbswiopt)) &
     &                    (1:lenv(swihelp(tabnumswiopt(jmod,jnbswiopt)))),'>'
                    ENDIF
                 ENDDO
                 IF ((.NOT.largaction).AND.(numhelp.EQ.jmod)) PRINT *, &
     &              '   for further details, type -action <numaction>'
!
! Display information on one specific action of this SESAM module
                 IF (largaction) THEN
                   lengthmax=60
                   IF (lenv(tabhelporder(jmod,jactiondeb)).GT.0) THEN
                     PRINT *,' '
                     pos=1
                     DO WHILE (pos.LE.lenv(tabhelporder(jmod,jactiondeb)))
                       posfin=MIN(pos-1+lengthmax, &
     &                       lenv(tabhelporder(jmod,jactiondeb)))
                       IF (posfin.NE.lenv(tabhelporder(jmod,jactiondeb)))  &
     &                        THEN
                         DO WHILE (.NOT.( &
     &                     (posfin.EQ.lenv(tabhelporder(jmod,jactiondeb))) &
     &                     .OR.tabhelporder(jmod,jactiondeb) &
     &                     (posfin+1:posfin+1).EQ.' '))
                           posfin=posfin-1
                           IF (posfin.EQ.pos) THEN 
                             posfin=lenv(tabhelporder(jmod,jactiondeb))
                           ENDIF
                         ENDDO
                       ENDIF
                       IF (pos.EQ.1) THEN
                         PRINT 1046,'  TASK : ', &
     &                     tabhelporder(jmod,jactiondeb)(pos:posfin)
                       ELSE
                         PRINT 1047,tabhelporder(jmod,jactiondeb)(pos:posfin)
                       ENDIF
                       pos=posfin+2
                     ENDDO
                   ENDIF
                 ENDIF
!
                 PRINT *, ' '
               ENDIF
            ENDDO
         ENDDO
!
      CASE (-2)
!
         DO jarg=2,2
            PRINT *, ' ',switab(jarg)(1:lenv(switab(jarg))), &
     &           ' <',swihelp(jarg)(1:lenv(swihelp(jarg))), &
     &           '>  ==> possible help switches'
         ENDDO
         PRINT *, ' '
!
! -2.- Display help information on file extensions
! ------------------------------------------------
         textmodel= '[model]'
         textdate = '[date]'
         textjpr  = '[NBR]'
         textfvar = '[extxyoz]'
         PRINT *, ' 0  => unavailable format'
         PRINT *, ' 1  => available format'
         PRINT *, ' 2  => available format if all vectors ', &
     &        'cannot fit in memory'
         PRINT *, ' 3  => available format if all vectors ', &
     &        'can fit in memory'
!
! -2.1- Observation vector format
         PRINT *, ' '
         PRINT 1010,'state','in','out', &
     &        'Observation vector format (internal convention : o,obs)'
         DO jextobs=1,nbextobs
            IF (extobstab(jextobs).NE.'-----') THEN
               IF (extobsunit(jextobs)) THEN
                  WRITE(text1,*) textmodel(1:lenv(textmodel)), &
     &                 textdate(1:lenv(textdate)), &
     &                 extobstab(jextobs)
               ELSE
                  WRITE(text1,*) textmodel(1:lenv(textmodel)),'#', &
     &                 textdate(1:lenv(textdate)), &
     &                 extobstab(jextobs)
               ENDIF
               PRINT 1020,'(',ETAT(extobsbool(jextobs)),')', &
     &              (extobsmem(jextobs)/10), &
     &              (MOD(extobsmem(jextobs),10)), &
     &              text1(1:lenv(text1))
            ENDIF
         ENDDO
!     
! -2.2- Observation database format :
         PRINT *, ' '
         PRINT 1010,'state','in','out', &
     &        'Observation database format (internal convention : dbs)'
         DO jextdbs=1,nbextdbs
            IF (extdbstab(jextdbs).NE.'-----') THEN
               WRITE(text1,*) textmodel(1:lenv(textmodel)), &
     &              textdate(1:lenv(textdate)), &
     &              extdbstab(jextdbs)
               PRINT 1020,'(',ETAT(extdbsbool(jextdbs)),')', &
     &              (extdbsmem(jextdbs)/10), &
     &              (MOD(extdbsmem(jextdbs),10)), &
     &              text1(1:lenv(text1))
            ENDIF
         ENDDO
!
! -2.3- data section vector format
         PRINT *, ' '
         PRINT 1010,'state','in','out', &
     &        'Data section vector format (internal convention : y,dta)'
         DO jextdta=1,nbextdta
            IF (extdtatab(jextdta).NE.'-----') THEN
               IF (extdtaunit(jextdta)) THEN
                  WRITE(text1,*) textmodel(1:lenv(textmodel)), &
     &                 textdate(1:lenv(textdate)), &
     &                 extdtatab(jextdta)
               ELSE
                  WRITE(text1,*) textmodel(1:lenv(textmodel)),'#', &
     &                 textdate(1:lenv(textdate)), &
     &                 extdtatab(jextdta)
               ENDIF
               PRINT 1020,'(',ETAT(extdtabool(jextdta)),')', &
     &              (extdtamem(jextdta)/10), &
     &              (MOD(extdtamem(jextdta),10)), &
     &              text1(1:lenv(text1))
            ENDIF
         ENDDO
!
! -2.4- State vector format
         PRINT *, ' '
         PRINT 1010,'state','in','out', &
     &        'State vector format (internal convention : x,var)'
         DO jextvar=1,nbextvar
            IF (extvartab(jextvar).NE.'-----') THEN
               IF (extvarunit(jextvar)) THEN
                  WRITE(text1,*) textmodel(1:lenv(textmodel)), &
     &                 textdate(1:lenv(textdate)), &
     &                 extvartab(jextvar)
               ELSE
                  WRITE(text1,*) textmodel(1:lenv(textmodel)),'#', &
     &                 textdate(1:lenv(textdate)), &
     &                 extvartab(jextvar)
               ENDIF
               PRINT 1020,'(',ETAT(extvarbool(jextvar)),')', &
     &              (extvarmem(jextvar)/10), &
     &              (MOD(extvarmem(jextvar),10)), &
     &              text1(1:lenv(text1))
            ENDIF
         ENDDO
!
! -2.5- Local data section vector format
         PRINT *, ' '
         PRINT 1010,'state','in','out', &
     &        'Local data section vector format (internal convention : z,zon)'
         DO jextzon=1,nbextzon
            IF (extzontab(jextzon).NE.'-----') THEN
               IF (extzonunit(jextzon)) THEN
                  WRITE(text1,*) textmodel(1:lenv(textmodel)), &
     &                 textdate(1:lenv(textdate)), &
     &                 extzontab(jextzon)
               ELSE
                  WRITE(text1,*) textmodel(1:lenv(textmodel)),'#', &
     &                 textdate(1:lenv(textdate)), &
     &                 extzontab(jextzon)
               ENDIF
               PRINT 1020,'(',ETAT(extzonbool(jextzon)),')', &
     &              (extzonmem(jextzon)/10), &
     &              (MOD(extzonmem(jextzon),10)), &
     &              text1(1:lenv(text1))
            ENDIF
         ENDDO
!
! -2.6- Covariance directory format :
         PRINT *, ' '
         PRINT 1010,'state','in','out', &
     &        'Covariance directory format (internal convention : bas)'
         DO jextbas=1,nbextbas
            IF (extbastab(jextbas).NE.'-----') THEN
               WRITE(text1,*) textmodel(1:lenv(textmodel)), &
     &              textdate(1:lenv(textdate)), &
     &              textjpr(1:lenv(textjpr)),'.', &
     &              textfvar(1:lenv(textfvar)), &
     &              extbastab(jextbas)
               PRINT 1020,'(',ETAT(extbasbool(jextbas)),')', &
     &              (extbasmem(jextbas)/10), &
     &              (MOD(extbasmem(jextbas),10)), &
     &              text1(1:lenv(text1))
            ENDIF
         ENDDO
         PRINT *, ' '
         PRINT *, ' ',textjpr(1:lenv(textjpr)), &
     &        ' = rank(+1) of covariance matrix   ex: 010'
         PRINT *, ' ',textfvar(1:lenv(textfvar)), &
     &        ' = file extension (from var,dta,obs,zon types)'
         PRINT *, ' '
!
      CASE (-3)
!
         DO jarg=2,2
            PRINT *, ' ',switab(jarg)(1:lenv(switab(jarg))), &
     &           ' <',swihelp(jarg)(1:lenv(swihelp(jarg))), &
     &           '>  ==> possible help switches'
         ENDDO
         PRINT *, ' '
!
! -3.- Display help information on SESAM configuration
! ----------------------------------------------------
! -help <config> -varmsk -dtamsk -weight -bias -list
!
! Build table (swibool1) of switches
! that are possible in this help context
         DO jarg=1,nbarg
            swibool1(jarg)= .FALSE.
            swibool2(jarg)= .FALSE.
         ENDDO
! ===> -help
         swibool1(2) = swibool(2)
         larghelp=.FALSE.
! ===> -varmsk (optional)
         swibool1(4) = swibool(4)
! ===> -dtamsk (optional)
         swibool1(5) = swibool(5)
! ===> -memory (optional)
         swibool1(9) = swibool(9)
! ===> -inxbas (optional)
         swibool1(10) = swibool(10)
! ===> -configobs (optional)
         swibool1(52) = swibool(52)
! ===> -fixjpx (option)
         swibool1(53) = swibool(53)
! ===> -configzon (optional)
         swibool1(71) = swibool(71)
!
! Build empty table (modbool1) of module
         modbool1(1:nbmod)= .FALSE.
         nmode=0
!
! Copy tables of possible file extensions
         extbasbool1(0:nbextbas)= extbasbool(0:nbextbas)
         extdbsbool1(0:nbextdbs)= extdbsbool(0:nbextdbs)
         extdtabool1(0:nbextdta)= extdtabool(0:nbextdta)
         extobsbool1(0:nbextobs)= extobsbool(0:nbextobs)
         extvarbool1(0:nbextvar)= extvarbool(0:nbextvar)
         extzonbool1(0:nbextzon)= extzonbool(0:nbextzon)
!
! Read user commandline (again!), check switch and argument validity
         argloop2: DO jloop = 1, nbargmax
!
#if defined GETARG
            CALL getarg ((jloop*2-1),swi)
            CALL getarg ((jloop*2),wrd)
#else 
            CALL get_command_argument ((jloop*2-1),swi)
            CALL get_command_argument ((jloop*2),wrd)
#endif
! ==> End of loop
            IF (swi.EQ.'  ') THEN 
               EXIT argloop2
            ENDIF
! ==> Write switch and argument
            IF (nprint.GE.1) THEN
               WRITE(numout,*) ' switch/argument: ', swi(1:lenv(swi)), &
     &                                               wrd(1:lenv(wrd))
            ENDIF
! ==> Check switch validity
            IF (swi((swilg+1):(swilg+1)).NE.' ') EXIT argloop2
            IF (wrd.EQ.'  ') EXIT argloop2
            IF (.NOT.(validswi(swi))) EXIT argloop2
! ==> Check if all switches are valid in help context
            jcas=0
            DO jarg=1,nbarg
               IF (swibool(jarg)) THEN
                  IF ((swi.EQ.switab(jarg)).OR.((jarg.EQ.2).AND.(swi.EQ.switab(2)(1:2)))) THEN
                     IF (swibool1(jarg)) THEN
                        jcas=jarg
                     ELSE
                        EXIT argloop2
                     ENDIF
                  ENDIF                
               ENDIF                
            ENDDO
!
!  ==> Check argument validity
            IF (jcas.EQ.0) EXIT argloop2
            CALL affectarg(wrd,jmod,jcas,swibool1,modbool1, &
     &           extbasbool1,extdbsbool1,extdtabool1, &
     &           extobsbool1,extvarbool1,extzonbool1)
!
! ==> Set swibool2 to .TRUE. if switch is present in the commandline
            swibool2(jcas)= .TRUE.
!     
         END DO argloop2
!
! Read SESAM mask
         CALL readmsk
!
! Evaluate SESAM configuration
         CALL evalconfig
!
! Display help on optional switches in <-help config> context
         PRINT 1050,'state','optional switches :'
         jarg=4
         PRINT 1060,'(',ETAT(swibool(jarg)),')',                   &
     &        switab(jarg)(1:lenv(switab(jarg))), &
     &        '<',swihelp(jarg)(1:lenv(swihelp(jarg))),'>'
         jarg=5
         PRINT 1060,'(',ETAT(swibool(jarg)),')',                   &
     &        switab(jarg)(1:lenv(switab(jarg))), &
     &        '<',swihelp(jarg)(1:lenv(swihelp(jarg))),'>'
         jarg=9
         PRINT 1060,'(',ETAT(swibool(jarg)),')',                   &
     &        switab(jarg)(1:lenv(switab(jarg))), &
     &        '<',swihelp(jarg)(1:lenv(swihelp(jarg))),'>'
         jarg=10
         PRINT 1060,'(',ETAT(swibool(jarg)),')',                   &
     &        switab(jarg)(1:lenv(switab(jarg))), &
     &        '<',swihelp(jarg)(1:lenv(swihelp(jarg))),'>'
         jarg=52
         PRINT 1060,'(',ETAT(swibool(jarg)),')',                   &
     &        switab(jarg)(1:lenv(switab(jarg))), &
     &        '<',swihelp(jarg)(1:lenv(swihelp(jarg))),'>'
         jarg=71
         PRINT 1060,'(',ETAT(swibool(jarg)),')',                   &
     &        switab(jarg)(1:lenv(switab(jarg))), &
     &        '<',swihelp(jarg)(1:lenv(swihelp(jarg))),'>'
!
! Display information on current SESAM configuration
!
! --- display Vx object configuration
         PRINT *, ' '
         PRINT *, ' Vx object CONFIGURATION (var):'
         PRINT *, ' ------------------------------'
         PRINT 1070, &
     &        'Variables',' var_ind ',' var_nbr ','  jpxend '
         DO jvar=1,varend
            indvar=var_ord(jvar)
            PRINT 1080, var_nam(indvar), &
     &           var_ind(indvar),var_nbr(indvar),jpxend
         ENDDO
! --- display Vy object configuration
         PRINT *, ' '
         PRINT *, ' Vy object CONFIGURATION (dta):'
         PRINT *, ' ------------------------------'
         PRINT *, ' '
         PRINT 1070, &
     &        'Variables',' dta_ind ',' dta_nbr ','  jpyend '
         DO jdta=1,dtaend
            inddta=dta_ord(jdta)
            PRINT 1080, dta_nam(inddta), &
     &              dta_ind(inddta),dta_nbr(inddta),jpyend
         ENDDO
! --- display Vz object configuration
         IF (existzon) THEN
            PRINT *, ' '
            PRINT *, ' Vz object CONFIGURATION (zon):'
            PRINT *, ' ------------------------------'
            PRINT *, ' '
            PRINT 1130,' zon_jpi ',' zon_jpj ',' zon_jpk ',' zon_jpt ', &
     &           '   jpz   ','  jpbub  '
            PRINT 1140,zon_jpi,zon_jpj,zon_jpk,zon_jpt, &
     &           jpz,jpbub
            PRINT *, ' '
            PRINT *, ' '
         ENDIF
! --- display Vo object configuration
         IF (existobs) THEN
            PRINT *, ' '
            PRINT *, ' Vo object CONFIGURATION (obs):'
            PRINT *, ' ------------------------------'
            PRINT *, ' '
            PRINT 1090,'Variables',' obs_ind ',' obs_nbr ', &
     &           '  jpoend ',' jpitpend'
            DO jobs=1,obsend
               indobs=obs_ord(jobs)
               inddbs=obsnord(jobs)
               PRINT 1100, obs_nam(indobs,inddbs), &
     &              obs_ind(indobs,inddbs),obs_nbr(indobs,inddbs), &
     &              jpoend,obs_itp(indobs,inddbs)
            ENDDO
         ENDIF
! --- display Io object configuration
         IF (existdbs) THEN
            PRINT *, ' '
            PRINT *, ' Io object CONFIGURATION (dbs):'
            PRINT *, ' ------------------------------'
            PRINT *, ' '
            PRINT 1090, '  obs_nam  ',' jpdbsend  '
            IF (largaffectobs) THEN
               PRINT 1100, argaffectobs(1:lenv(argaffectobs)),jpdbsend
            ELSE
               PRINT 1100, 'undefined observation name',jpdbsend
            ENDIF
         ENDIF
! --- display Vx blocks configuration
         PRINT *, ' '
         PRINT *, ' Vx object block CONFIGURATION :'
         PRINT *, ' -------------------------------'
         PRINT *, ' '
         PRINT *, ' size of Vx blocks in memory = ',jpx
         PRINT *, ' number of blocks = ',jpnxend
         PRINT *, ' '
      CASE DEFAULT
         PRINT *, 'Help not available on that subject'
      END SELECT
!
      text1=version
      text2=version_year
      PRINT *, 'Program ',text1(1:lenv(text1)), &
     &         '  year ',text2(1:lenv(text2))
      PRINT *, ' '
      STOP 'End of SESAM online help'
!
! --- format definitions
!
 1010 FORMAT(A5,1X,A5,2X,A5,8X,A)
 1020 FORMAT(A1,A3,A1,5X,I1,5X,I1,5X,A)
 1030 FORMAT(A6,1X,A5,1X,A4,8X,A)
 1040 FORMAT(2X,I2,2X,A1,A3,A1,3X,I1,3X,A)
 1045 FORMAT(18X,A)
 1046 FORMAT(A,A)
 1047 FORMAT(10X,A)
 1050 FORMAT(A5,8X,A)
 1060 FORMAT(A1,A3,A1,3X,A,1X,A1,A,A1)
 1070 FORMAT(2X,A9,3(1X,"|",2X,A9,1X))
 1080 FORMAT(5X,A3,4X,"|",1X,3(I11,3X))
 1090 FORMAT(2X,A9,1X,"|",2X,A11,1X)
 1100 FORMAT(5X,A5,2X,"|",1X,I11,3X)
 1110 FORMAT(1X,A11,4(1X,"|",1X,A11))
 1120 FORMAT(A12,1X,"|",1X,4(I11,3X))
 1130 FORMAT(2X,A9,5(1X,"|",2X,A9,1X))
 1140 FORMAT(2X,I9,5(1X,"|",2X,I9,1X))
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      CHARACTER(len=3) FUNCTION ETAT (state)
!
      IMPLICIT NONE
      LOGICAL, intent(in) :: state
!
      IF (state) THEN
         ETAT=' ON'
      ELSE
         ETAT='OFF'
      ENDIF
!
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE arghelp
