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
CC---------------------------------------------------------------------
CC
CC                         DEFCST.SWITCH.H
CC
CC---------------------------------------------------------------------
CC
CC  Purpose :
CC  ---------
CC  Include file with SESAM structure parameters
CC  (definition of commandline switches, definition of SESAM
CC  modules, list of possible actions for every modules,
CC  list of SESAM objects, description of available file
CC  formats for every objects, parameters with SESAM online help, ...)
CC   
CC  Modifications :
CC  -------------
CC  original     : 97-12 (C.E. Testut)
CC  modification : 01-06 (C.E. Testut)
CC  modification : 03-02 (J.M. Brankart)
CC
CC---------------------------------------------------------------------
C A/ Define all commandline switches available in SESAM
C 1. switab = switch name
C 2. swihelp = help about the nature of the argument of the switch
C 3. swibool = switch is active (TRUE or FALSE)
C 4. swiopt = switch is optional (0=false, 1=general optional switch,
C                   2=optional switch particular to a few moodules).
C 5. arg<switchname> = initialize the argument of the switch (to empty value)
C 6. larg<switchname> = argument of the switch is not empty (TRUE or FALSE)
C                       (initialized to FALSE).
      DO jarg=1,nbarg+1
         SELECT CASE (jarg)
         CASE (1)
            switab(jarg) = '-mode'
            WRITE (swihelp(jarg),'(A,A,A)') 
     $           '(groa,lroa,broa,filt,greg,lreg,breg),',
     $           '(adap,diff,geof,leof,beof),',
     $           '(intf,obsv,oerr,oper,zone)'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argmod       = ''
            largmod      = .FALSE.
         CASE (2)
            switab(jarg) = '-help'
            swihelp(jarg)= 'mode,ext,config,all'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arghelp      = ''
            larghelp     = .FALSE.
         CASE (3)
            switab(jarg) = '-config'
            swihelp(jarg)= 'SESAM configuration'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 1
            arglist      = ''
            larglist     = .FALSE.
         CASE (4)
            switab(jarg) = '-varmsk'
            swihelp(jarg)= 'mask.[bimg,dimg,cdf]'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 1
            argvarmsk    = ''
            largvarmsk   = .FALSE.
         CASE (5)
            switab(jarg) = '-dtamsk'
            swihelp(jarg)= 'mask.[bdta,dta,cdta]'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 1
            argdtamsk    = ''
            largdtamsk   = .FALSE.
         CASE (6)
            switab(jarg) = '-weight'
            swihelp(jarg)= 'file_xyo'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            argweight    = ''
            largweight   = .FALSE.
         CASE (7)
            switab(jarg) = '-oestd'
            swihelp(jarg)= 'file_xyo'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            argoestd     = ''
            largoestd    = .FALSE.
         CASE (8)
            switab(jarg) = '-bias'
            swihelp(jarg)= 'file_xyo'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            argbias      = ''
            largbias     = .FALSE.
         CASE (9)
            switab(jarg) = '-outinfo'
            swihelp(jarg)= 'listing file'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 1
            argoutinfo   = ''
            largoutinfo  = .FALSE.
         CASE (10)
            switab(jarg) = '-inxbas'
            swihelp(jarg)= 'dir_xbas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginxbas    = ''
            larginxbas   = .FALSE.
         CASE (11)
            switab(jarg) = '-indbs'
            swihelp(jarg)= 'file_dbs'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argindbs     = ''
            largindbs    = .FALSE.
         CASE (12)
            switab(jarg) = '-indta'
            swihelp(jarg)= 'file_xy'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argindta     = ''
            largindta    = .FALSE.
         CASE (13)
            switab(jarg) = '-inobs'
            swihelp(jarg)= 'file_xyo'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginobs     = ''
            larginobs    = .FALSE.
         CASE (14)
            switab(jarg) = '-invar'
            swihelp(jarg)= 'file_x'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginvar     = ''
            larginvar    = .FALSE.
         CASE (15)
            switab(jarg) = '-inxbasref'
            swihelp(jarg)= 'dir_xbas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginxbasref = ''
            larginxbasref= .FALSE.
         CASE (16)
            switab(jarg) = '-indbsref'
            swihelp(jarg)= 'file_dbs'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argindbsref  = ''
            largindbsref = .FALSE.
         CASE (17)
            switab(jarg) = '-indtaref'
            swihelp(jarg)= 'file_xy'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argindtaref  = ''
            largindtaref = .FALSE.
         CASE (18)
            switab(jarg) = '-inobsref'
            swihelp(jarg)= 'file_xyo'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginobsref  = ''
            larginobsref = .FALSE.
         CASE (19)
            switab(jarg) = '-invarref'
            swihelp(jarg)= 'file_x'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginvarref  = ''
            larginvarref = .FALSE.
         CASE (20)
            switab(jarg) = '-outxbas'
            swihelp(jarg)= 'dir_xbas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutxbas   = ''
            largoutxbas  = .FALSE.
         CASE (21)
            switab(jarg) = '-outdta'
            swihelp(jarg)= 'file_y'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutdta    = ''
            largoutdta   = .FALSE.
         CASE (22)
            switab(jarg) = '-outobs'
            swihelp(jarg)= 'file_o'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutobs    = ''
            largoutobs   = .FALSE.
         CASE (23)
            switab(jarg) = '-outvar'
            swihelp(jarg)= 'file_x'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutvar    = ''
            largoutvar   = .FALSE.
         CASE (24)
            switab(jarg) = '-outxbasref'
            swihelp(jarg)= 'dir_xbas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutxbasref= ''
            largoutxbasref= .FALSE.
         CASE (25)
            switab(jarg) = '-outdtaref'
            swihelp(jarg)= 'file_y'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutdtaref = ''
            largoutdtaref= .FALSE.
         CASE (26)
            switab(jarg) = '-outobsref'
            swihelp(jarg)= 'file_o'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutobsref = ''
            largoutobsref= .FALSE.
         CASE (27)
            switab(jarg) = '-outvarref'
            swihelp(jarg)= 'file_x'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutvarref = ''
            largoutvarref= .FALSE.
         CASE (30)
            switab(jarg) = '-reducedta'
            swihelp(jarg)= 'file_xy'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            argreducedta = ''
            largreducedta= .FALSE.
         CASE (31)
            switab(jarg) = '-outparadap'
            swihelp(jarg)= 'adaptive parameters'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutparadap= ''
            largoutparadap= .FALSE.
         CASE (32)
            switab(jarg) = '-outrz'
            swihelp(jarg)= 'file_rz'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            argoutrz      = ''
            largoutrz    = .FALSE.
         CASE (33)
            switab(jarg) = '-typeoper'
            swihelp(jarg)= 'operation type'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argtypeoper  = ''
            largtypeoper = .FALSE.
         CASE (34)
            switab(jarg) = '-typedtadiag'
            swihelp(jarg)= 'dta_type diag'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 1
            argtypedtadiag= ''
            largtypedtadiag= .FALSE.
         CASE (35)
            switab(jarg) = '-diffobsref'
            swihelp(jarg)= 'file_o'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argdiffobsref= ''
            largdiffobsref= .FALSE.
         CASE (36)
            switab(jarg) = '-diffdtaref'
            swihelp(jarg)= 'file_xy'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argdiffdtaref= ''
            largdiffdtaref= .FALSE.
        CASE (37)
            switab(jarg) = '-diffvarref'
            swihelp(jarg)= 'file_x'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argdiffvarref= ''
            largdiffvarref= .FALSE.
         CASE (38)
            switab(jarg) = '-diffobsorg'
            swihelp(jarg)= 'file_o'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argdiffobsorg= ''
            largdiffobsorg= .FALSE.
         CASE (39)
            switab(jarg) = '-diffdtaorg'
            swihelp(jarg)= 'file_xy'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argdiffdtaorg= ''
            largdiffdtaorg= .FALSE.
         CASE (40)
            switab(jarg) = '-diffvarorg'
            swihelp(jarg)= 'file_x'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argdiffvarorg= ''
            largdiffvarorg= .FALSE.
         CASE (41)
            switab(jarg) = '-inerrdta'
            swihelp(jarg)= 'file_y'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginerrdta  = ''
            larginerrdta = .FALSE.
         CASE (42)
            switab(jarg) = '-outerrdta'
            swihelp(jarg)= 'file_y'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argouterrdta = ''
            largouterrdta= .FALSE.
         CASE (43)
            switab(jarg) = '-instddta'
            swihelp(jarg)= 'file_xy'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginstddta  = ''
            larginstddta = .FALSE.
         CASE (44)
            switab(jarg) = '-outstddta'
            swihelp(jarg)= 'file_y'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutstddta = ''
            largoutstddta= .FALSE.
         CASE (45)
            switab(jarg) = '-inybasref'
            swihelp(jarg)= 'dir_xybas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginybasref = ''
            larginybasref= .FALSE.
         CASE (46)
            switab(jarg) = '-outybasref'
            swihelp(jarg)= 'dir_ybas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutybasref= ''
            largoutybasref= .FALSE.
        CASE (47)
            switab(jarg) = '-outbiasdta'
            swihelp(jarg)= 'file_y'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutbiasdta= ''
            largoutbiasdta= .FALSE.
        CASE (48)
            switab(jarg) = '-affectobs'
            swihelp(jarg)= 'obs_nam'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argaffectobs = ''
            largaffectobs= .FALSE.
        CASE (49)
            switab(jarg) = '-nullobs'
            swihelp(jarg)= 'ALL,obs_nam'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argnullobs   = ''
            largnullobs  = .FALSE.
         CASE (50)
            switab(jarg) = '-inpartobs'
            swihelp(jarg)= 'file_xyo'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            arginpartobs = ''
            larginpartobs= .FALSE.
         CASE (51)
            switab(jarg) = '-outdiaghst'
            swihelp(jarg)= 'file_hst'
            swibool(jarg)= .FALSE.
            swiopt(jarg) = 2
            argoutdiaghst= ''
            largoutdiaghst= .FALSE.
         CASE (52)
            switab(jarg) = '-configobs'
            swihelp(jarg)= 'file_o'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argconfigobs = ''
            largconfigobs= .FALSE.
         CASE (53)
            switab(jarg) = '-fixjpx'
            swihelp(jarg)= 'block size'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 1
            argfixjpx    = ''
            largfixjpx= .FALSE.
         CASE (54)
            switab(jarg) = '-reducevar'
            swihelp(jarg)= 'file_x'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            argreducevar = ''
            largreducevar= .FALSE.
         CASE (55)
            switab(jarg) = '-scale'
            swihelp(jarg)= 'scale'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            argscale     = ''
            largscale    = .FALSE.
         CASE (56)
            switab(jarg) = '-biasdbs'
            swihelp(jarg)= 'dbs_file'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            argbiasdbs   = ''
            largbiasdbs  = .FALSE.
         CASE (57)
            switab(jarg) = '-inpartvar'
            swihelp(jarg)= 'file_x'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginpartvar = ''
            larginpartvar= .FALSE.
         CASE (58)
            switab(jarg) = '-outpartvar'
            swihelp(jarg)= 'file_x'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutpartvar = ''
            largoutpartvar= .FALSE.
         CASE (59)
            switab(jarg) = '-inzon'
            swihelp(jarg)= 'file_z'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginzon     = ''
            larginzon    = .FALSE.
         CASE (60)
            switab(jarg) = '-outzon'
            swihelp(jarg)= 'file_z'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutzon    = ''
            largoutzon   = .FALSE.
         CASE (61)
            switab(jarg) = '-oecorrel'
            swihelp(jarg)= 'dir_zbas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            argoecorrel  = ''
            largoecorrel = .FALSE.
         CASE (62)
            switab(jarg) = '-fecorrel'
            swihelp(jarg)= 'dir_zbas'
            swibool(jarg)= .FALSE.
            swiopt(jarg) = 2
            argfecorrel  = ''
            largfecorrel = .FALSE.
         CASE (63)
            switab(jarg) = '-zonindex'
            swihelp(jarg)= 'bubble_index'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argzonindex  = ''
            largzonindex = .FALSE.
         CASE (64)
            switab(jarg) = '-incfg'
            swihelp(jarg)= 'file_cfg'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argincfg     = ''
            largincfg    = .FALSE.
         CASE (65)
            switab(jarg) = '-inybas'
            swihelp(jarg)= 'dir_xybas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginybas    = ''
            larginybas   = .FALSE.
         CASE (66)
            switab(jarg) = '-outybas'
            swihelp(jarg)= 'dir_ybas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutybas   = ''
            largoutybas  = .FALSE.
         CASE (67)
            switab(jarg) = '-inobas'
            swihelp(jarg)= 'dir_xyobas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginobas    = ''
            larginobas   = .FALSE.
         CASE (68)
            switab(jarg) = '-outobas'
            swihelp(jarg)= 'dir_obas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutobas   = ''
            largoutobas  = .FALSE.
         CASE (69)
            switab(jarg) = '-inzbas'
            swihelp(jarg)= 'dir_zbas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginzbas    = ''
            larginzbas   = .FALSE.
         CASE (70)
            switab(jarg) = '-outzbas'
            swihelp(jarg)= 'dir_zbas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutzbas   = ''
            largoutzbas  = .FALSE.
         CASE (71)
            switab(jarg) = '-configzon'
            swihelp(jarg)= 'file_z'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argconfigzon = ''
            largconfigzon= .FALSE.
         CASE (72)
            switab(jarg) = '-coefrmax'
            swihelp(jarg)= 'real_number'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            argcoefrmax  = ''
            largcoefrmax = .FALSE.
         CASE (73)
            switab(jarg) = '-disable'
            swihelp(jarg)= 'string'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            argdisable   = ''
            largdisable  = .FALSE.
         CASE (74)
            switab(jarg) = '-inobasref'
            swihelp(jarg)= 'dir_xyobas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginobasref = ''
            larginobasref= .FALSE.
         CASE (75)
            switab(jarg) = '-fixjpz'
            swihelp(jarg)= 'integer'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 1
            argfixjpz    = ''
            largfixjpz= .FALSE.
         CASE (76)
            switab(jarg) = '-inparadap'
            swihelp(jarg)= 'adaptive parameters'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            arginparadap = ''
            larginparadap= .FALSE.
         CASE (77)
            switab(jarg) = '-fixjpu'
            swihelp(jarg)= 'integer'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 1
            argfixjpu    = ''
            largfixjpu= .FALSE.
         CASE (78)
            switab(jarg) = '-inptzon'
            swihelp(jarg)= 'file_z'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            arginptzon   = ''
            larginptzon  = .FALSE.
         CASE (79)
            switab(jarg) = '-outptzon'
            swihelp(jarg)= 'file_z'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutptzon  = ''
            largoutptzon = .FALSE.
         CASE (80)
            switab(jarg) = '-action'
            swihelp(jarg)= 'action_number'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 1
            argaction    = ''
            largaction   = .FALSE.
         CASE (81)
            switab(jarg) = '-inzbasref'
            swihelp(jarg)= 'dir_zbas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginzbasref = ''
            larginzbasref= .FALSE.
         CASE (82)
            switab(jarg) = '-inzonref'
            swihelp(jarg)= 'file_z'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginzonref  = ''
            larginzonref = .FALSE.
         CASE (83)
            switab(jarg) = '-inoptcfg'
            swihelp(jarg)= 'file_cfg'
            swibool(jarg)= .TRUE. 
            swiopt(jarg) = 2
            arginoptcfg  = ''
            larginoptcfg = .FALSE.
         CASE (84)
            switab(jarg) = '-connect'
            swihelp(jarg)= 'file_cfg'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 1
            arginoptcfg  = ''
            larginoptcfg = .FALSE.
         CASE (85)
            switab(jarg) = '-incstr'
            swihelp(jarg)= 'dir_xbas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argincstr = ''
            largincstr= .FALSE.
         CASE (86)
            switab(jarg) = '-outobasref'
            swihelp(jarg)= 'dir_obas'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            argoutobasref= ''
            largoutobasref= .FALSE.
         CASE (87)
            switab(jarg) = '-insmocfg'
            swihelp(jarg)= 'file_cfg'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            arginsmocfg  = ''
            larginsmocfg = .FALSE.
         CASE (88)
            switab(jarg) = '-outsmocfg'
            swihelp(jarg)= 'file_cfg'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 2
            argoutsmocfg  = ''
            largoutsmocfg = .FALSE.
         CASE (89)
            switab(jarg) = '-inrz'
            swihelp(jarg)= 'file_rz'
            swibool(jarg)= .TRUE.
            swiopt(jarg) = 0
            arginrz       = ''
            larginrz     = .FALSE.
         CASE (nbarg+1)
            IF (jarg.NE.90) STOP 'defcst.switch.h : bad nbarg parameter'
         CASE DEFAULT
            switab(jarg) = '--------'
            swihelp(jarg)= '--------'
            swibool(jarg)= .FALSE.
            swiopt(jarg) = 0
         END SELECT  
      ENDDO
C
C B/ Define all SESAM modules
C 1. modtab = module name
C 2. modhelp = help about the nature of the module
C 3. modbool = module is active (TRUE or FALSE)
C 4. modtype = type of module (1=analysis modules, 2=diagnostic modules,
C                   3=utility modules)
      DO jmod=1,nbmod+1
         SELECT CASE (jmod)
         CASE (1)
            modtab(jmod) = 'intf'
            modhelp(jmod)= 
     $           'Interface an object from a file format to another'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 3
         CASE (2)
            modtab(jmod) = 'corr'
            modhelp(jmod)= 
     $           'Compute correlations from reduced order cov. matrix'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 2
         CASE (3)
            modtab(jmod) = 'tgop'
            modhelp(jmod)=
     $           'Operations on truncated Gaussian pdf'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 1
         CASE (4)
            modtab(jmod) = 'obsv'
            modhelp(jmod)= 'Observation management'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 3
         CASE (5)
            modtab(jmod) = 'filt'
            modhelp(jmod)= 'Low-pass filter'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 1
         CASE (6)
            modtab(jmod) = 'adap'
            modhelp(jmod)= 'Adaptive parameter estimation'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 2
         CASE (7)
            modtab(jmod) = 'diff'
            modhelp(jmod)= 
     $           'Compute RMS difference between vector objects'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 2
         CASE (8)
            modtab(jmod) = 'oerr'
            modhelp(jmod)= 'Observation error management'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 3
         CASE (9)
            modtab(jmod) = 'oper'
            modhelp(jmod)= 'Arithmetic operations'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 3
         CASE (10)
            modtab(jmod) = 'geof'
            modhelp(jmod)= 'Global EOF decomposition'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 2
         CASE (11)
            modtab(jmod) = 'leof'
            modhelp(jmod)= 'Local EOF decomposition'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 2
         CASE (12)
            modtab(jmod) = 'beof'
            modhelp(jmod)= 'Bubble EOF decomposition'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 2
         CASE (13)
            modtab(jmod) = 'zone'
            modhelp(jmod)= 'Local data section management'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 3
         CASE (14)
            modtab(jmod) = 'groa'
            modhelp(jmod)= 'Global reduced order analysis'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 1
         CASE (15)
            modtab(jmod) = 'lroa'
            modhelp(jmod)= 'Local reduced order analysis'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 1
         CASE (16)
            modtab(jmod) = 'broa'
            modhelp(jmod)= 'Bubble reduced order analysis'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 1
         CASE (17)
            modtab(jmod) = 'greg'
            modhelp(jmod)= 'Global linear regression'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 1
         CASE (18)
            modtab(jmod) = 'lreg'
            modhelp(jmod)= 'Local linear regression'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 1
         CASE (19)
            modtab(jmod) = 'breg'
            modhelp(jmod)= 'Bubble linear regression'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 1
         CASE (20)
            modtab(jmod) = 'vari'
            modhelp(jmod)=
     $           'Compute variances from reduced order cov. matrix'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 2
         CASE (21)
            modtab(jmod) = 'anam'
            modhelp(jmod)=
     $           'Anamorphosis operations'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 2
         CASE (22)
            modtab(jmod) = 'scor'
            modhelp(jmod)=
     $           'Probabilistic scores (CRPS,...)'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 2
         CASE (23)
            modtab(jmod) = 'spct'
            modhelp(jmod)=
     $           'Spectral transformation'
            modbool(jmod)= .TRUE.
            modtype(jmod)= 2
         CASE (nbmod+1)
            IF (jmod.NE.24) STOP 'defcst.switch.h : bad nbmod parameter'
         CASE DEFAULT
            modtab(jmod) = '----'
            modhelp(jmod)= ''
            modbool(jmod)= .FALSE.
            modtype(jmod)= 3
         END SELECT         
      ENDDO
C
C C/ Define SESAM help keywords
C 1. helptab = help keyword name
C 2. helpbool = help keyword is active (TRUE or FALSE)
      DO jhelp=1,nbhelp+1
         SELECT CASE (jhelp)
         CASE (1)
            helptab(jhelp) = 'all'
            helpbool(jhelp)= .TRUE.
         CASE (2)
            helptab(jhelp) = 'ext'
            helpbool(jhelp)= .TRUE.
         CASE (3)
            helptab(jhelp) = 'config'
            helpbool(jhelp)= .TRUE.
         CASE (nbhelp+1)
            IF (jhelp.NE.4) STOP 'defcst.switch.h : bad nbhelp parameter'
         CASE DEFAULT
            helptab(jhelp) = '----'
            helpbool(jhelp)= .FALSE.
         END SELECT         
      ENDDO
C
C D/ Define SESAM covariance directory extensions
C 1. extbastab = covariance directory extension name
C 2. extbasbool = covariance directory extension is active (TRUE or FALSE)
C 3. extbasmem = type of format (unit=as output format, ten=as input format):
C                0=unavailable, 1=available, >1=available with restrictions
      DO jextbas=0,nbextbas+1
         SELECT CASE (jextbas)
         CASE (0)
            extbasbool(jextbas)= .FALSE.
         CASE (1)
            extbastab(jextbas) = '.bas'
            extbasbool(jextbas)= .TRUE.
            extbasmem(jextbas)= 11
         CASE (nbextbas+1)
            IF (jextbas.NE.2) STOP 'defcst.switch.h : bad nbextbas parameter'
         CASE DEFAULT
            extbastab(jextbas) = '-----'
            extbasbool(jextbas)= .FALSE.
            extbasmem(jextbas)= 11
         END SELECT         
      ENDDO
C
C E/ Define SESAM observation database file extensions
C 1. extdbstab = observation database file extension name
C 2. extdbsbool = observation database file extension is active (TRUE or FALSE)
C 3. extdbsmem = type of format (unit=as output format, ten=as input format):
C                0=unavailable, 1=available, >1=available with restrictions
      DO jextdbs=0,nbextdbs+1
         SELECT CASE (jextdbs)
         CASE (0)
            extdbsbool(jextdbs)= .FALSE.
         CASE (1)
            extdbstab(jextdbs) = '.tra'
            extdbsbool(jextdbs)= .FALSE.
            extdbsmem(jextdbs)= 10
         CASE (2)
            extdbstab(jextdbs) = '.crg'
            extdbsbool(jextdbs)= .TRUE.
            extdbsmem(jextdbs)= 10
         CASE (3)
            extdbstab(jextdbs) = '.ncdbs'
            extdbsbool(jextdbs)= .TRUE.
            extdbsmem(jextdbs)= 10
         CASE (4)
            extdbstab(jextdbs) = '.adbs'
            extdbsbool(jextdbs)= .TRUE.
            extdbsmem(jextdbs)= 10
         CASE (nbextdbs+1)
            IF (jextdbs.NE.5) STOP 'defcst.switch.h : bad nbextdbs parameter'
         CASE DEFAULT
            extdbstab(jextdbs) = '-----'
            extdbsbool(jextdbs)= .FALSE.
            extdbsmem(jextdbs)= 10
         END SELECT         
      ENDDO
C
C F/ Define SESAM dta file extensions
C 1. extdtatab = dta file extension name
C 2. extdtabool = dta file extension is active (TRUE or FALSE)
C 3. extdtamem = type of format (unit=as output format, ten=as input format):
C                0=unavailable, 1=available, >1=available with restrictions
C 4. extdtaunit = different variables are in the same file (TRUE or FALSE)
      DO jextdta=0,nbextdta+1
         SELECT CASE (jextdta)
         CASE (0)
            extdtabool(jextdta)= .FALSE.
         CASE (1)
            extdtatab(jextdta) = '.bdta'
            extdtabool(jextdta)= .TRUE.
            extdtamem(jextdta) = 11
            extdtaunit(jextdta)= .FALSE.
         CASE (2)
            extdtatab(jextdta) = '.dta'
            extdtabool(jextdta)= .TRUE.
            extdtamem(jextdta) = 11
            extdtaunit(jextdta)= .FALSE.
         CASE (3)
            extdtatab(jextdta) = '.cdta'
            extdtabool(jextdta)= .TRUE.
            extdtamem(jextdta) = 11
            extdtaunit(jextdta)= .TRUE.
         CASE (4)
            extdtatab(jextdta) = '.ncdta'
            extdtabool(jextdta)= .TRUE.
            extdtamem(jextdta) = 11
            extdtaunit(jextdta)= .FALSE.
         CASE (5)
            extdtatab(jextdta) = '.yplm'
            extdtabool(jextdta)= .FALSE.
            extdtamem(jextdta) = 11
            extdtaunit(jextdta)= .FALSE.
         CASE (nbextdta+1)
            IF (jextdta.NE.6) STOP 'defcst.switch.h : bad nbextdta parameter'
         CASE DEFAULT
            extdtatab(jextdta) = '-----'
            extdtabool(jextdta)= .FALSE.
            extdtamem(jextdta) = 11
            extdtaunit(jextdta)= .TRUE.
         END SELECT         
      ENDDO
C
C G/ Define SESAM obs file extensions
C 1. extobstab = obs file extension name
C 2. extobsbool = obs file extension is active (TRUE or FALSE)
C 3. extobsmem = type of format (unit=as output format, ten=as input format):
C                0=unavailable, 1=available, >1=available with restrictions
C 4. extobsunit = different variables are in the same file (TRUE or FALSE)
      DO jextobs=0,nbextobs+1
         SELECT CASE (jextobs)
         CASE (0)
            extobsbool(jextobs)= .FALSE.
         CASE (1)
            extobstab(jextobs) = '.obs'
            extobsbool(jextobs)= .TRUE.
            extobsmem(jextobs)= 11
            extobsunit(jextobs)= .FALSE.
         CASE (2)
            extobstab(jextobs) = '.cobs'
            extobsbool(jextobs)= .TRUE.
            extobsmem(jextobs)= 11
            extobsunit(jextobs)= .FALSE.
         CASE (3)
            extobstab(jextobs) = '.oplm'
            extobsbool(jextobs)= .FALSE.
            extobsmem(jextobs)= 11
            extobsunit(jextobs)= .FALSE.
         CASE (nbextobs+1)
            IF (jextobs.NE.4) STOP 'defcst.switch.h : bad nbextobs parameter'
         CASE DEFAULT
            extobstab(jextobs) = '-----'
            extobsbool(jextobs)= .FALSE.
            extobsmem(jextobs)= 11
            extobsunit(jextobs)= .TRUE.
         END SELECT         
      ENDDO
C
C H/ Define SESAM var file extensions
C 1. extvartab = var file extension name
C 2. extvarbool = var file extension is active (TRUE or FALSE)
C 3. extvarmem = type of format (unit=as output format, ten=as input format):
C                0=unavailable, 1=available, >1=available with restrictions
C 4. extvarunit = different variables are in the same file (TRUE or FALSE)
      DO jextvar=1,nbextvar+1
         SELECT CASE (jextvar)
         CASE (0)
            extvarbool(jextvar)= .FALSE.
         CASE (1)
            extvartab(jextvar) = '.bimg'
            extvarbool(jextvar)= .TRUE.
            extvarmem(jextvar) = 33
            extvarunit(jextvar)= .FALSE.
         CASE (2)
            extvartab(jextvar) = '.dimg'
            extvarbool(jextvar)= .TRUE.
            extvarmem(jextvar) = 33
            extvarunit(jextvar)= .FALSE.
         CASE (3)
            extvartab(jextvar) = '.cdf'
            extvarbool(jextvar)= .TRUE.
            extvarmem(jextvar) = 33
            extvarunit(jextvar)= .TRUE.
         CASE (4)
            extvartab(jextvar) = '.nc'
            extvarbool(jextvar)= .TRUE.
            extvarmem(jextvar) = 33
            extvarunit(jextvar)= .FALSE.
         CASE (5)
            extvartab(jextvar) = '.cpak'
            extvarbool(jextvar)= .TRUE.
            extvarmem(jextvar) = 11
            extvarunit(jextvar)= .TRUE.
         CASE (6)
            extvartab(jextvar) = '.dpak'
            extvarbool(jextvar)= .FALSE.
            extvarmem(jextvar) = 33
            extvarunit(jextvar)= .TRUE.
         CASE (7)
            extvartab(jextvar) = '.mpak'
            extvarbool(jextvar)= .FALSE.
            extvarmem(jextvar) = 33
            extvarunit(jextvar)= .TRUE.
         CASE (8)
            extvartab(jextvar) = '.rst'
            extvarbool(jextvar)= .TRUE.
            extvarmem(jextvar) = 33
            extvarunit(jextvar)= .TRUE.
         CASE (nbextvar+1)
            IF (jextvar.NE.9) STOP 'defcst.switch.h : bad nbextvar parameter'
         CASE DEFAULT
            extvartab(jextvar) = '-----'
            extvarbool(jextvar)= .FALSE.
            extvarmem(jextvar) = 33
            extvarunit(jextvar)= .TRUE.
         END SELECT         
      ENDDO
C
C I/ Define SESAM zon file extensions
C 1. extzontab = zon file extension name
C 2. extzonbool = zon file extension is active (TRUE or FALSE)
C 3. extzonmem = type of format (unit=as output format, ten=as input format):
C                0=unavailable, 1=available, >1=available with restrictions
C 4. extzonunit = different variables are in the same file (TRUE or FALSE)
      DO jextzon=0,nbextzon+1
         SELECT CASE (jextzon)
         CASE (0)
            extzonbool(jextzon)= .FALSE.
         CASE (1)
            extzontab(jextzon) = '.czon'
            extzonbool(jextzon)= .TRUE.
            extzonmem(jextzon)= 11
            extzonunit(jextzon)= .TRUE.
         CASE (2)
            extzontab(jextzon) = '.zplm'
            extzonbool(jextzon)= .FALSE.
            extzonmem(jextzon)= 11
            extzonunit(jextzon)= .TRUE.
         CASE (nbextzon+1)
            IF (jextzon.NE.3) STOP 'defcst.switch.h : bad nbextzon parameter'
         CASE DEFAULT
            extzontab(jextzon) = '-----'
            extzonbool(jextzon)= .FALSE.
            extzonmem(jextzon)= 11
            extzonunit(jextzon)= .TRUE.
         END SELECT         
      ENDDO
C
C J/ Define all possible actions for every SESAM modules
C 1. tabvalorder(jmod,jaction)  = Number of switches for action 'jaction' (x10),
C                                 and type of action (unit): obsolete feature
C 2. tabhelporder(jmod,jaction) = Help for action 'jaction' of module 'jmod'
C 3. taborder(jmod,jaction,:)   = List of switch indices defining action 'jaction'
C                                    of module 'jmod'
C 4. tabnumswiopt(jmod,:)   = List of optional switch indices
C                                    particular to module 'jmod'
C
C -- Initialisation :
      taborder(:,:,:)=0
      tabvalorder(:,:)=0
      tabhelporder(:,:)=''
      tabnumswiopt(:,:)=0
C
C -I.- INTF module : Interface an object from a file format to another
C --------------------------------------------------------------------
      jmod=1
C -I.1- =>: -invar <file_x> -outvar <file_x>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Interface a state vector',
     $    ' (invar) to another (outvar)'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=23
C -I.2- =>: -invar <file_x> -outdta <file_y>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Extract a data section',
     $    ' (outdta) from a state vector (invar)'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=21
C -I.3- =>: -invar <file_x> -outobs <file_o> -configobs <file_o>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Extract the observation',
     $     ' vector (outobs) from a state vector (invar)',
     $     ' according to the Hxo operator (configobs)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=22
      taborder(jmod,jaction,3)=52
C -I.4- =>: -invar <file_x> -outzon <file_z> -configzon <file_z>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Interface a state vector',
     $    ' to a local data section vector (outzon)',
     $    ' according to the Hxz operator (configzon)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=60
      taborder(jmod,jaction,3)=71
C -I.5- =>: -indta <file_xy> -outdta <file_y>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Interface a',
     $    ' data section vector (indta) to another (outdta)'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=21
C -I.6- =>: -indta <file_xy> -outobs <file_o> -configobs <file_o>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Extract the observation',
     $     ' vector (outobs) from a data section vector (indta)',
     $     ' according to the Hxo operator (configobs)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=22
      taborder(jmod,jaction,3)=52
C -I.7- =>: -indta <file_xy> -outzon <file_z> -configzon <file_z>
      jaction=7
      WRITE(tabhelporder(jmod,jaction),*) 'Interface a',
     $    ' data section vector (indta) to a',
     $    ' local data section vector (outzon)',
     $    ' according to the Hxz operator (configzon)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=60
      taborder(jmod,jaction,3)=71
C -I.8- =>: -inobs <file_xyo> -outobs <file_o> -configobs <file_o>
      jaction=8
      WRITE(tabhelporder(jmod,jaction),*) 'Interface a',
     $    ' observation vector (inobs) to another (outobs)',
     $     ' according to the Ho operator (configobs)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=22
      taborder(jmod,jaction,3)=52
C -I.9- =>: -inzon <file_z> -outzon <file_z> -configzon <file_z>
      jaction=9
      WRITE(tabhelporder(jmod,jaction),*) 'Interface a',
     $    ' local data section vector (inzon) to another (outzon)',
     $     ' according to the Hz operator (configzon)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=59
      taborder(jmod,jaction,2)=60
      taborder(jmod,jaction,3)=71
C -I.10- =>: -inxbas <file_xbas> -outxbas <file_xbas>
      jaction=10
      WRITE(tabhelporder(jmod,jaction),*) 'Interface an ensemble of',
     $    ' state vectors (inxbas) to another (outxbas)'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=20
C -I.11- =>: -inxbas <file_xbas> -outybas <file_ybas>
      jaction=11
      WRITE(tabhelporder(jmod,jaction),*) 'Interface an ensemble of',
     $    ' state vectors (inxbas) to an ensemble of',
     $    ' data section vectors (outybas)'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=66
C -I.12- =>: -inxbas <file_xbas> -outobas <file_obas> -configobs <file_o>
      jaction=12
      WRITE(tabhelporder(jmod,jaction),*) 'Interface an ensemble of',
     $    ' state vectors (inxbas) to an ensemble of',
     $    ' observation vectors (outobas)',
     $    ' according to the Hxo operator (configobs)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=68
      taborder(jmod,jaction,3)=52
C -I.13- =>: -inxbas <file_xbas> -outzbas <file_zbas> -configzon <file_z>
      jaction=13
      WRITE(tabhelporder(jmod,jaction),*) 'Interface an ensemble of',
     $    ' state vectors (inxbas) to an ensemble of',
     $    ' local data section vectors (outzbas)',
     $    ' according to the Hxz operator (configzon)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=70
      taborder(jmod,jaction,3)=71
C -I.14- =>: -inybas <file_xybas> -outybas <file_ybas>
      jaction=14
      WRITE(tabhelporder(jmod,jaction),*) 'Interface an ensemble of',
     $    ' data section vectors (inybas) to another (outybas)'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=66
C -I.15- =>: -inybas <file_xybas> -outobas <file_obas> -configobs <file_o>
      jaction=15
      WRITE(tabhelporder(jmod,jaction),*) 'Interface an ensemble of',
     $    ' data section vectors (inybas) to an ensemble of',
     $    ' observation section vectors (outobas)',
     $    ' according to the Hyo operator (configobs)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=68
      taborder(jmod,jaction,3)=52
C -I.16- =>: -inybas <file_xybas> -outzbas <file_zbas> -configzon <file_z>
      jaction=16
      WRITE(tabhelporder(jmod,jaction),*) 'Interface an ensemble of',
     $    ' data section vectors (inybas) to an ensemble of',
     $    ' local data section vectors (outzbas)',
     $    ' according to the Hz operator (configzon)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=70
      taborder(jmod,jaction,3)=71
C -I.17- =>: -inobas <file_xyobas> -outobas <file_obas> -configobs <file_o>
      jaction=17
      WRITE(tabhelporder(jmod,jaction),*) 'Interface an ensemble of',
     $    ' observation vectors (inobas) to another',
     $    ' (outobas)',
     $    ' according to the Hyo operator (configobs)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=68
      taborder(jmod,jaction,3)=52
C -I.18- =>: -inzbas <file_zbas> -outzbas <file_zbas> -configzon <file_z>
      jaction=18
      WRITE(tabhelporder(jmod,jaction),*) 'Interface an ensemble of',
     $    ' local data section vectors (inzbas) to another',
     $    ' (outzbas)',
     $    ' according to the Hz operator (configzon)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=69
      taborder(jmod,jaction,2)=70
      taborder(jmod,jaction,3)=71
C
C -II.- CORR module : Compute correlation coefficients
C ----------------------------------------------------
      jmod=2
C -II.1- =>: -inxbas <file_xbas> -outvar <file_x> -incfg <file_cfg>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Compute correlations',
     $   ' coefficients (outvar) from state reduced order',
     $   ' covariance matrix (inxbas)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=23
      taborder(jmod,jaction,3)=64
C -II.2- =>: -inybas <file_ybas> -outdta <file_y> -incfg <file_cfg>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Compute correlations',
     $   ' coefficients (outdta) from state reduced order',
     $   ' covariance matrix (inybas)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=21
      taborder(jmod,jaction,3)=64
C -II.3- =>: -inobas <file_obas> -outobs <file_o> -incfg <file_cfg> -configobs <file_o>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Compute correlations',
     $   ' coefficients (outobs) from state reduced order',
     $   ' covariance matrix (inobas)'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=22
      taborder(jmod,jaction,3)=64
      taborder(jmod,jaction,4)=52
C
C -III.- TGOP module : Truncated Gaussian operations
C --------------------------------------------------
      jmod=3
C -III.1- =>: -outvar <file_x> -invar <file_x> -inxbas <file_xbas> -incstr <file_xbas>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Compute maximum likelihood or variance minimizing estimator',
     $    ' of TG pdf from the TG parameters'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=14
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=85
C -III.2- =>: -outvar <file_x> -invar <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -incstr <file_xbas>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Compute maximum likelihood or variance minimizing estimator',
     $    ' of TG pdf (and associated error covariance) from the TG parameters'
      tabvalorder(jmod,jaction)=52 
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=14
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=20
      taborder(jmod,jaction,5)=85
C -III.3- =>: -invar <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -incstr <file_xbas>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Compute sample of a TG pdf from the TG parameters'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=20
      taborder(jmod,jaction,4)=85
C -III.4- =>: -outvar <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -incstr <file_xbas>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Compute maximum likelihood estimator',
     $    ' of the TG parameters from a sample'
      tabvalorder(jmod,jaction)=42 
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=20
      taborder(jmod,jaction,4)=85
C -III.5- =>: -outvar <file_x> -inxbas <file_xbas> -inxbasref <file_xbas>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Compare two samples using',
     $    ' Kolmogorov-Smirnov test'
      tabvalorder(jmod,jaction)=32 
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=15
C -III.6- =>: -outvar <file_x> -invar <file_x> -inxbas <file_xbas> -incstr <file_xbas> -inpartvar <file_x>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Compute local maximum likelihood',
     $    ' or variance minimizing estimator of TG pdf from the TG parameters'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=14
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=85
      taborder(jmod,jaction,5)=57
C
C -IV.- OBSV module : Observation management 
C ------------------------------------------
      jmod=4
C --- optional switches -reducedta (30) -scale(55) -biasdbs(56) -inoptcfg(83)
      tabnumswiopt(jmod,1)=30
      tabnumswiopt(jmod,2)=55
      tabnumswiopt(jmod,3)=56
      tabnumswiopt(jmod,4)=83
C -IV.1- =>: -indbs <file_dbs> -outobs <file_o> -affectobs <obs_nam>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Extract an',
     $   ' observation vector (outobs) from a database (indbs)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=11
      taborder(jmod,jaction,2)=22
      taborder(jmod,jaction,3)=48
C -IV.2- =>: -outobs <file_o> -nullobs <obs_nam>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Create an',
     $   ' empty observation vector'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=22
      taborder(jmod,jaction,2)=49
C -IV.3- =>: -inobs <file_o> -outdta <file_y> -configobs <file_o>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Make an approximate',
     $   ' conversion from an observation vector (inobs)',
     $   ' to a data section vector (outdta)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=21
      taborder(jmod,jaction,3)=52
C -IV.4- =>: -indta <file_y> -outobs <file_o>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Convert a',
     $   ' data section vector (indta) into',
     $   ' observation vector (outobs)'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=22
C -IV.5- =>: -indta <file_y> -outobs <file_o> -configobs <file_o>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Extract the',
     $   ' observation equivalent (outobs) from',
     $   ' a data section vector (indta)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=22
      taborder(jmod,jaction,3)=52
C
C -V.- FILT module : Low-pass filter
C ----------------------------------
C
      jmod=5
C --- optional switches -inptzon(78)
      tabnumswiopt(jmod,1)=78
C -V.1- =>: -incfg <incfg.cfg> -outzon <file_z>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Generate a',
     $   ' homogeneous filter (outzon)'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=64
      taborder(jmod,jaction,2)=60
C -V.2- =>: -incfg <incfg.cfg> -indta <file_xy> -outzon <file_z>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Generate an',
     $   ' inhomogeneous filter (outzon)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=64
      taborder(jmod,jaction,2)=12
      taborder(jmod,jaction,3)=60
C -V.3- =>: -indta <file_xy> -inzon <file_z> -outdta <file_y>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Apply the',
     $   ' filter (inzon) on a data section vector',
     $   ' (indta) using a Kyy gain'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=59
      taborder(jmod,jaction,3)=21
C -V.4- =>: -inobs <file_xy> -inzon <file_z> -outdta <file_y> -configobs <file_o>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Apply the',
     $   ' filter (inzon) on an observation vector',
     $   ' (inobs) using a Kyo gain'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=59
      taborder(jmod,jaction,3)=21
      taborder(jmod,jaction,4)=52
C
C -VI.- ADAP module : Adaptive parameter estimation
C -------------------------------------------------
C
      jmod=6
C --- optional switches -disable(73) -inparadap(76) -outrz(32) -inpartobs(50)
      tabnumswiopt(jmod,1)=73
      tabnumswiopt(jmod,2)=76
      tabnumswiopt(jmod,3)=32
      tabnumswiopt(jmod,4)=50
C -VI.1- =>: -incfg <file_cfg> -outparadap <file_crz>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Compute optimal',
     $   ' parameter estimates'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=64
      taborder(jmod,jaction,2)=31
C
C -VII.- DIFF module : Compute RMS difference between vector objects
C ------------------------------------------------------------------
C
      jmod=7
C --- optional switches -outdiaghst(51) -bias(8)
      tabnumswiopt(jmod,1)=51
      tabnumswiopt(jmod,2)=8
C -VII.1- =>: -invar <file_x> -diffvarref <file_x>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Compute global',
     $   ' RMS difference by level and variable between',
     $   ' two state vectors'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=37
C -VII.2- =>: -indta <file_y> -diffdtaref <file_y>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Compute global',
     $   ' RMS difference by level and variable between',
     $   ' two data section vectors'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=36
C -VII.3- =>: -inobs <file_o> -diffobsref <file_o> -configobs <file_o>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Compute global',
     $   ' RMS difference by level and variable between',
     $   ' two observation vectors'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=35
      taborder(jmod,jaction,3)=52
C -VII.4- =>: -invar <file_x> -diffvarref <file_x> -diffvarorg <file_x>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Compute global',
     $   ' RMS difference by level and variable between',
     $   ' two state vectors'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=37
      taborder(jmod,jaction,3)=40
C -VII.5- =>: -indta <file_y> -diffdtaref <file_y> -diffdtaorg <file_y>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Compute global',
     $   ' RMS difference by level and variable between',
     $   ' two data section vectors'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=36
      taborder(jmod,jaction,3)=39
C -VII.6- =>: -inobs <file_o> -diffobsref <file_o> -diffobsorg <file_o> -configobs <file_o>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Compute global',
     $   ' RMS difference by level and variable between',
     $   ' two observation vectors'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=35
      taborder(jmod,jaction,3)=38
      taborder(jmod,jaction,4)=52
C -VII.7- =>: -invar <file_x> -diffvarref <file_x> -inpartvar <file_x> -incfg <incfg.cfg>
      jaction=7
      WRITE(tabhelporder(jmod,jaction),*) 'Compute regional',
     $   ' RMS difference by variables between',
     $   ' two states vectors',
     $   ' using a partition in subsystems (inpartvar,incfg)'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=37
      taborder(jmod,jaction,3)=57
      taborder(jmod,jaction,4)=64
C -VII.8- =>: -indta <file_y> -diffdtaref <file_y> -inpartvar <file_x> -incfg <incfg.cfg>
      jaction=8
      WRITE(tabhelporder(jmod,jaction),*) 'Compute regional',
     $   ' RMS difference by variables between',
     $   ' two data section vectors',
     $   ' using a partition in subsystems (inpartvar,incfg)'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=36
      taborder(jmod,jaction,3)=57
      taborder(jmod,jaction,4)=64
C -VII.9- =>: -inobs <file_o> -diffobsref <file_o> -configobs <file_o> -inpartvar <file_x> -incfg <incfg.cfg>
      jaction=9
      WRITE(tabhelporder(jmod,jaction),*) 'Compute regional',
     $   ' RMS difference by variables between',
     $   ' two observation vectors',
     $   ' using a partition in subsystems (inpartvar,incfg)'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=35
      taborder(jmod,jaction,3)=52
      taborder(jmod,jaction,4)=57
      taborder(jmod,jaction,5)=64
C -VII.10- =>: -invar <file_x> -diffvarref <file_x> -diffvarorg <file_x> -inpartvar <file_x> -incfg <incfg.cfg>
      jaction=10
      WRITE(tabhelporder(jmod,jaction),*) 'Compute regional',
     $   ' RMS difference by variables between',
     $   ' two state vectors',
     $   ' using a partition in subsystems (inpartvar,incfg)'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=37
      taborder(jmod,jaction,3)=40
      taborder(jmod,jaction,4)=57
      taborder(jmod,jaction,5)=64
C -VII.11- =>: -indta <file_y> -diffdtaref <file_y> -diffdtaorg <file_y> -inpartvar <file_x> -incfg <incfg.cfg>
      jaction=11
      WRITE(tabhelporder(jmod,jaction),*) 'Compute regional',
     $   ' RMS difference by variables between',
     $   ' two data section vectors',
     $   ' using a partition in subsystems (inpartvar,incfg)'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=36
      taborder(jmod,jaction,3)=39
      taborder(jmod,jaction,4)=57
      taborder(jmod,jaction,5)=64
C -VII.12- =>: -inobs <file_o> -diffobsref <file_o> -diffobsorg <file_o> -configobs <file_o> -inpartvar <file_x> -incfg <incfg.cfg>
      jaction=12
      WRITE(tabhelporder(jmod,jaction),*) 'Compute regional',
     $   ' RMS difference by variables between',
     $   ' two observation vectors',
     $   ' using a partition in subsystems (inpartvar,incfg)'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=35
      taborder(jmod,jaction,3)=38
      taborder(jmod,jaction,4)=52
      taborder(jmod,jaction,5)=57
      taborder(jmod,jaction,6)=64
C -VII.13- =>: -invar <file_x> -inpartvar <file_x> -incfg <incfg.cfg>
      jaction=13
      WRITE(tabhelporder(jmod,jaction),*) 'Compute regional',
     $   ' mean and std by variable for',
     $   ' a state vector',
     $   ' using a partition in subsystems (inpartvar,incfg)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=57
      taborder(jmod,jaction,3)=64
C -VII.14- =>: -indta <file_y> -inpartvar <file_x> -incfg <incfg.cfg>
      jaction=14
      WRITE(tabhelporder(jmod,jaction),*) 'Compute regional',
     $   ' mean and std by variable for',
     $   ' a data section vector',
     $   ' using a partition in subsystems (inpartvar,incfg)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=57
      taborder(jmod,jaction,3)=64
C -VII.15- =>: -inobs <file_o> -configobs <file_o> -inpartvar <file_x> -incfg <incfg.cfg>
      jaction=15
      WRITE(tabhelporder(jmod,jaction),*) 'Compute regional',
     $   ' mean and std by variable for',
     $   ' an observation vector',
     $   ' using a partition in subsystems (inpartvar,incfg)'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=52
      taborder(jmod,jaction,3)=57
      taborder(jmod,jaction,4)=64
C
C -VIII.- OERR module : Observation error management
C --------------------------------------------------
C
      jmod=8
C --- optional switches -weight(6) -oestd(7) -bias(8) -oecorrel(61)
      tabnumswiopt(jmod,1)=6
      tabnumswiopt(jmod,2)=7
      tabnumswiopt(jmod,3)=8
      tabnumswiopt(jmod,4)=61
C -VIII.1- =>: -outerrdta <file_y> -incfg <incfg.cfg>
      jaction=1
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=42
      taborder(jmod,jaction,2)=64
C -VIII.2- =>: -inerrdta <file_y> -outerrdta <file_y> -incfg <incfg.cfg>
      jaction=2
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=41
      taborder(jmod,jaction,2)=42
      taborder(jmod,jaction,3)=64
C -VIII.3- =>: -instddta <file_xy> -indta <file_xy> -outstddta <file_y> -typeoper <operation>
      jaction=3
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=43
      taborder(jmod,jaction,2)=12
      taborder(jmod,jaction,3)=44
      taborder(jmod,jaction,4)=33
C -VIII.4- =>: -instddta <file_xy> -inobs <file_xyo> -configobs <file_o> -outstddta <file_y> -typeoper <operation> 
      jaction=4
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=43
      taborder(jmod,jaction,2)=13
      taborder(jmod,jaction,3)=52
      taborder(jmod,jaction,4)=44
      taborder(jmod,jaction,5)=33
C -VIII.5- =>: -inobas <file_xyobas> -outerrdta <file_y>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Compute',
     $   ' the representativeness error (outerrdta) from',
     $   ' an ensemble of observation vectors (inobas)'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=42
C -VIII.6- =>: -inobas <file_xyobas> -inybas <file_xybas> -outerrdta <file_y> -outstddta <file_y>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Compute',
     $   ' the correlation (outerrdta) and',
     $   ' RMS difference (outstddta) between',
     $   ' an ensemble of data section vectors (inybas) and',
     $   ' an ensemble of observation vectors (inobas)'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=42
      taborder(jmod,jaction,4)=44
C -VIII.7- =>: -inobas <file_xyobas> -inybas <file_xybas> -inybasref <file_xybas> -outerrdta <file_y> -outstddta <file_y>
      jaction=7
      tabvalorder(jmod,jaction)=50
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=45
      taborder(jmod,jaction,4)=42
      taborder(jmod,jaction,5)=44
C -VIII.8- =>: -inobas <file_xyobas> -inybas <file_xybas> -outbiasdta <file_y>
      jaction=8
      WRITE(tabhelporder(jmod,jaction),*) 'Compute the mean difference',
     $   ' between an ensemble of data section vectors (inybas) and',
     $   ' an ensemble of observation vectors (inobas)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=47
C -VIII.9- =>: -inobas <file_xyobas> -outdtaref <file_y> -outstddta <file_y>
      jaction=9
      WRITE(tabhelporder(jmod,jaction),*) 'Compute',
     $   ' the mean (outdtaref) and the standard deviation (outstddta) of',
     $   ' an ensemble of observation vectors (inobas)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=25
      taborder(jmod,jaction,3)=44
C
C -IX.- OPER module : Arithmetic operation
C ----------------------------------------
C
      jmod=9
C -IX.1- =>: -incfg <file_cfg> -outvar <file_x> -typeoper <operation>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Compute a',
     $   ' state vector (outvar) from parameter file',
     $   ' by specified arithmetic operation:',
     $   ' list of operations = (mean)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=64
      taborder(jmod,jaction,2)=23
      taborder(jmod,jaction,3)=33
C -IX.2- =>: -incfg <file_cfg> -outdta <file_y> -typeoper <operation>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Compute a',
     $   ' data section vector (outdta) from parameter file',
     $   ' by specified arithmetic operation:',
     $   ' list of operations = (mean)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=64
      taborder(jmod,jaction,2)=21
      taborder(jmod,jaction,3)=33
C -IX.3- =>: -incfg <file_cfg> -outobs <file_o> -configobs <file_o> -typeoper <operation>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Compute an',
     $   ' observation vector (outobs) from parameter file',
     $   ' by specified arithmetic operation:',
     $   ' list of operations = (mean)'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=64
      taborder(jmod,jaction,2)=22
      taborder(jmod,jaction,3)=33
      taborder(jmod,jaction,4)=52
C -IX.4- =>: -incfg <file_cfg> -outzon <file_z> -configzon <file_z> -typeoper <operation>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Compute a',
     $   ' local data section vector (outzon) from parameter file',
     $   ' by specified arithmetic operation:',
     $   ' list of operations = ()'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=64
      taborder(jmod,jaction,2)=60
      taborder(jmod,jaction,3)=33
      taborder(jmod,jaction,4)=71
C -IX.5- =>: -invar <file_x> -outvar <file_x> -typeoper <operation>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Compute a',
     $   ' state vector (outvar) by arithmetic operation on',
     $   ' input state vector (invar):',
     $   ' list of operations = (+_a,-_a,x_a,/_a,inv_a,pow_a,cst_a,csp_a,',
     $   'pow2,sqrt,abs,opp,inv)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=23
      taborder(jmod,jaction,3)=33
C -IX.6- =>: -indta <file_xy> -outdta <file_y> -typeoper <operation>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Compute a',
     $   ' data section vector (outdta)',
     $   ' by arithmetic operation on',
     $   ' input data section vector (indta):',
     $   ' list of operations = (+_a,-_a,x_a,/_a,inv_a,pow_a,cst_a,csp_a,',
     $   'pow2,sqrt,abs,opp,inv)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=21
      taborder(jmod,jaction,3)=33
C -IX.7- =>: -inobs <file_xyo> -outobs <file_o> -configobs <file_o> -typeoper <operation>
      jaction=7
      WRITE(tabhelporder(jmod,jaction),*) 'Compute an',
     $   ' observation vector (outobs)',
     $   ' by arithmetic operation on',
     $   ' input observation vector (inobs):',
     $   ' list of operations = (+_a,-_a,x_a,/_a,inv_a,pow_a,cst_a,csp_a,',
     $   'pow2,sqrt,abs,opp,inv)'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=22
      taborder(jmod,jaction,3)=33
      taborder(jmod,jaction,4)=52
C -IX.8- =>: -inzon <file_z> -outzon <file_z> -configzon <file_z> -typeoper <operation>
      jaction=8
      WRITE(tabhelporder(jmod,jaction),*) 'Compute a',
     $   ' local data section vector (outzon)',
     $   ' by arithmetic operation on',
     $   ' input local data section vector (inzon):',
     $   ' list of operations = (+_a,-_a,x_a,/_a,inv_a,pow_a,cst_a,csp_a,',
     $   'pow2,sqrt,abs,opp,inv)'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=59
      taborder(jmod,jaction,2)=60
      taborder(jmod,jaction,3)=33
      taborder(jmod,jaction,4)=71
C -IX.9- =>: -invar <file_x> -invarref <file_x> -outvar <file_x> -typeoper <operation>
      jaction=9
      WRITE(tabhelporder(jmod,jaction),*) 'Compute a',
     $   ' state vector (outvar)',
     $   ' by arithmetic operation on',
     $   ' two input state vectors, a=invar and b=invarref, ',
     $   ' list of operations = (rstopa,+,-,min,max,x,/,',
     $   '/sqrt,sqrt<a2-b2>,replace_spval)'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=19
      taborder(jmod,jaction,3)=23
      taborder(jmod,jaction,4)=33
C -IX.10- =>: -indta <file_xy> -indtaref <file_xy> -outdta <file_y> -typeoper <operation>
      jaction=10
      WRITE(tabhelporder(jmod,jaction),*) 'Compute a',
     $   ' data section vector (outdta)',
     $   ' by arithmetic operation on',
     $   ' two input data section vectors, a=indta and b=indtaref, ',
     $   ' list of operations = (+,-,min,max,x,/,',
     $   '/sqrt,sqrt<a2-b2>,replace_spval)'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=17
      taborder(jmod,jaction,3)=21
      taborder(jmod,jaction,4)=33
C -IX.11- =>: -inobs <file_xyo>-inobsref <file_xy> -outobs <file_o> -configobs <file_o> -typeoper <operation>
      jaction=11
      WRITE(tabhelporder(jmod,jaction),*) 'Compute an',
     $   ' observation vector (outobs)',
     $   ' by arithmetic operation on',
     $   ' two input observation vectors, a=inobs and b=inobsref, ',
     $   ' list of operations = (+,-,min,max,x,/,',
     $   '/sqrt,sqrt<a2-b2>,replace_spval)'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=18
      taborder(jmod,jaction,3)=22
      taborder(jmod,jaction,4)=33
      taborder(jmod,jaction,5)=52
C -IX.12- =>: -inzon <file_z>-inzonref <file_z> -outzon <file_z> -configzon <file_z> -typeoper <operation>
      jaction=12
      WRITE(tabhelporder(jmod,jaction),*) 'Compute a',
     $   ' local data section vector (outzon)',
     $   ' by arithmetic operation on',
     $   ' two input observation vectors, a=inzon and b=inzonref, ',
     $   ' list of operations = (+,-,min,max,x,/,',
     $   '/sqrt,sqrt<a2-b2>)'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=59
      taborder(jmod,jaction,2)=82
      taborder(jmod,jaction,3)=60
      taborder(jmod,jaction,4)=33
      taborder(jmod,jaction,5)=71
C
C -X.- GEOF module : Global EOF decomposition
C -------------------------------------------
      jmod=10
C --- optional switches -weight(6) -oestd(7) -reducevar(54) -disable(73)
      tabnumswiopt(jmod,1)=6
      tabnumswiopt(jmod,2)=7
      tabnumswiopt(jmod,3)=54
      tabnumswiopt(jmod,4)=73
C -X.1- =>: -inxbas <file_xbas> -outxbas <file_xbas>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Global Eof analysis of an',
     $    ' ensemble of state vectors (inxbas,jprmax).',
     $    ' The action writes the Eof vectors (outxbas).'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=20
C -X.2- =>: -inxbas <file_xbas> -outxbas <file_xbas> -inybas <file_xybas> 
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Global Eof analysis of an',
     $    ' ensemble of state vectors (inxbas,jprmax)',
     $    ' with the eigenspace defined by a',
     $    ' data section ensemble (inybas,jprmax)',
     $    ' The action writes the Eof vectors (outxbas).'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=20
      taborder(jmod,jaction,3)=65
C -X.3- =>: -inxbas <file_xbas> -outxbas <file_xbas> -inobas <file_xyobas> -configobs <file_o>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Global Eof analysis of an',
     $    ' ensemble of state vectors (inxbas,jprmax)',
     $    ' with the eigenspace defined by an',
     $    ' observation ensemble (inobas,jprmax)',
     $    ' The action writes the Eof vectors (outxbas).'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=20
      taborder(jmod,jaction,3)=67
      taborder(jmod,jaction,4)=52
C -X.4- =>: -inxbas <file_xbas> -outxbas <file_xbas> -inobas <file_xyobas> -configobs <file_o> -outobas <file_obas>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Global Eof analysis of an',
     $     ' ensemble of state vectors (inxbas,jprmax)',
     $     ' with the eigenspace defined by an',
     $     ' observation ensemble (inobas,jprmax)',
     $     ' The action writes the Eof vectors in state space (outxbas)',
     $     ' and in observation space (outobas)'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=20
      taborder(jmod,jaction,3)=67
      taborder(jmod,jaction,4)=52
      taborder(jmod,jaction,5)=68
C -X.5- =>: -inybas <file_xybas> -outybas <file_ybas>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Global Eof analysis of an',
     $     ' ensemble of data section vectors (inybas,jprmax).',
     $     ' The action writes the Eof vectors (outybas).'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=66
C -X.6- =>: -inybas <file_xybas> -outybas <file_ybas> -inobas <file_xyobas> -configobs <file_o>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Global Eof analysis of an',
     $     ' ensemble of data section vectors (inybas,jprmax)',
     $     ' with the eigenspace',
     $     ' defined by an observation ensemble (inobas,jprmax).',
     $     ' The action writes the Eof vectors (outybas).'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=66
      taborder(jmod,jaction,3)=67
      taborder(jmod,jaction,4)=52
C -X.7- =>: -inybas <file_xybas> -outybas <file_ybas> -inobas <file_xyobas> -configobs <file_o> -outobas <file_obas>
      jaction=7
      WRITE(tabhelporder(jmod,jaction),*) 'Global Eof analysis of an',
     $     ' ensemble of data section vectors (inybas,jprmax)',
     $     ' with the eigenspace',
     $     ' defined by an observation ensemble (inobas,jprmax).',
     $     ' The action writes the Eof vectors in data section space (outybas)',
     $     ' and in observation space (outobas)'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=66
      taborder(jmod,jaction,3)=67
      taborder(jmod,jaction,4)=52
      taborder(jmod,jaction,5)=68
C -X.8- =>: -inobas <file_xyobas> -outobas <file_obas> -configobs <file_o>
      jaction=8
      WRITE(tabhelporder(jmod,jaction),*) 'Global Eof analysis of an',
     $     ' ensemble of observation vectors (inobas,jprmax).',
     $     ' The action writes the Eof vectors (outobas).'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=68
      taborder(jmod,jaction,3)=52
C
C -XI.- LEOF module : Local EOF decomposition
C -------------------------------------------
      jmod=11
C --- optional switches -weight(6) -oestd(7) -reducevar(54) -oecorrel(61) -disable(73)
      tabnumswiopt(jmod,1)=6
      tabnumswiopt(jmod,2)=7
      tabnumswiopt(jmod,3)=54
      tabnumswiopt(jmod,4)=61
      tabnumswiopt(jmod,5)=73
C -XI.1- =>: -inxbas <file_xbas> -outxbas <file_xbas> -inybas <file_xybas> -inpartvar <file_x> -inzon <file_z>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Local Eof analysis of an',
     $    ' ensemble of state vectors (inobas,jprmax).'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=20
      taborder(jmod,jaction,3)=65
      taborder(jmod,jaction,4)=57
      taborder(jmod,jaction,5)=59
C -XI.2- =>: -inxbas <file_xbas> -outxbas <file_xbas> -inobas <file_xyobas> -configobs <file_o> -inpartvar <file_x> -inzon <file_z>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Local Eof analysis of an',
     $    ' ensemble of state vectors with the eigenspace defined by an',
     $    ' observation ensemble (inobas).'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=20
      taborder(jmod,jaction,3)=67
      taborder(jmod,jaction,4)=52
      taborder(jmod,jaction,5)=57
      taborder(jmod,jaction,6)=59
C -XI.- =>: -inxbas <file_xbas> -outxbas <file_xbas> -inobas <file_xyobas> -configobs <file_o> -outobas <file_obas> -inpartvar <file_x> -inzon <file_z>
C     jaction=
C     WRITE(tabhelporder(jmod,jaction),*) 'Local Eof analysis of an',
C    $    ' ensemble of state vectors with the eigenspace defined by an',
C    $    ' observation ensemble (inobas).',
C    $    ' The action writes the Eof vectors in state space (outxbas)',
C    $    ' and in observation space (outobas)'
C     tabvalorder(jmod,jaction)=72
C     taborder(jmod,jaction,1)=10
C     taborder(jmod,jaction,2)=20
C     taborder(jmod,jaction,3)=67
C     taborder(jmod,jaction,4)=52
C     taborder(jmod,jaction,5)=68
C     taborder(jmod,jaction,6)=57
C     taborder(jmod,jaction,7)=59
C -XI.3- =>: -inybas <file_xybas> -outybas <file_ybas> -inpartvar <file_x> -inzon <file_z>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Local Eof analysis of an',
     $   ' ensemble of data section vectors.'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=66
      taborder(jmod,jaction,3)=57
      taborder(jmod,jaction,4)=59
C -XI.4- =>: -inybas <file_xybas> -outybas <file_ybas> -inobas <file_xyobas> -configobs <file_o> -inpartvar <file_x> -inzon <file_z>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Local Eof analysis of an',
     $   ' ensemble of data section vectors with the eigenspace',
     $   ' defined by an observation ensemble (inobas).'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=66
      taborder(jmod,jaction,3)=67
      taborder(jmod,jaction,4)=52
      taborder(jmod,jaction,5)=57
      taborder(jmod,jaction,6)=59
C -XI.- =>: -inybas <file_xybas> -outybas <file_ybas> -inobas <file_xyobas> -configobs <file_o> -outobas <file_obas> -inpartvar <file_x> -inzon <file_z>
C      jaction=
C      WRITE(tabhelporder(jmod,jaction),*) 'Local Eof analysis of an',
C     $   ' ensemble of data section vectors with the eigenspace',
C     $   ' defined by an observation ensemble (inobas).'
C      tabvalorder(jmod,jaction)=72
C      taborder(jmod,jaction,1)=65
C      taborder(jmod,jaction,2)=66
C      taborder(jmod,jaction,3)=67
C      taborder(jmod,jaction,4)=52
C      taborder(jmod,jaction,5)=68
C      taborder(jmod,jaction,6)=57
C      taborder(jmod,jaction,7)=59
C -XI.- =>: -inobas <file_xyobas> -outobas <file_obas> -configobs <file_o> -inpartvar <file_x> -inzon <file_z>
C      jaction=
C      WRITE(tabhelporder(jmod,jaction),*) 'Local Eof analysis of an',
C     $   ' ensemble of observation vectors.'
C      tabvalorder(jmod,jaction)=52
C      taborder(jmod,jaction,1)=67
C      taborder(jmod,jaction,2)=68
C      taborder(jmod,jaction,3)=52
C      taborder(jmod,jaction,4)=57
C      taborder(jmod,jaction,5)=59
C
C -XII.- BEOF module : Bubble EOF decomposition
C ---------------------------------------------
      jmod=12
C --- optional switches -weight(6) -oestd(7) -reducevar(54) -oecorrel(61) -disable(73)
      tabnumswiopt(jmod,1)=6
      tabnumswiopt(jmod,2)=7
      tabnumswiopt(jmod,3)=54
      tabnumswiopt(jmod,4)=61
      tabnumswiopt(jmod,5)=73
C -XII.1- =>: -inxbas <file_xbas> -outxbas <file_xbas> -inybas <file_xybas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble Eof analysis of an',
     $    ' ensemble of state vectors.'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=20
      taborder(jmod,jaction,3)=65
      taborder(jmod,jaction,4)=70
      taborder(jmod,jaction,5)=57
      taborder(jmod,jaction,6)=59
C -XII.2- =>: -inxbas <file_xbas> -outxbas <file_xbas> -inobas <file_xyobas> -configobs <file_o> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble Eof analysis of an',
     $    ' ensemble of state vectors with the eigenspace defined by an',
     $    ' observation ensemble (inobas).'
      tabvalorder(jmod,jaction)=72
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=20
      taborder(jmod,jaction,3)=67
      taborder(jmod,jaction,4)=52
      taborder(jmod,jaction,5)=70
      taborder(jmod,jaction,6)=57
      taborder(jmod,jaction,7)=59
C -XII.3- =>: -inxbas <file_xbas> -outxbas <file_xbas> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble Eof analysis of an',
     $    ' ensemble of state vectors with the eigenspace defined by a',
     $    ' local data section ensemble (inzbas).'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=20
      taborder(jmod,jaction,3)=69
      taborder(jmod,jaction,4)=70
      taborder(jmod,jaction,5)=57
      taborder(jmod,jaction,6)=59
C -XII.4- =>: -inxbas <file_xbas> -outxbas <file_xbas> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z> -configobs <file_o>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble Eof analysis of an',
     $    ' ensemble of state vectors with the eigenspace defined by a',
     $    ' local data section ensemble (inzbas)',
     $    ' at observation locations (configobs).'
      tabvalorder(jmod,jaction)=72
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=20
      taborder(jmod,jaction,3)=69
      taborder(jmod,jaction,4)=70
      taborder(jmod,jaction,5)=57
      taborder(jmod,jaction,6)=59
      taborder(jmod,jaction,7)=52
C -XII.5- =>: -inybas <file_xybas> -outybas <file_ybas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble Eof analysis of an',
     $    ' ensemble of data section vectors.'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=66
      taborder(jmod,jaction,3)=70
      taborder(jmod,jaction,4)=57
      taborder(jmod,jaction,5)=59
C -XII.6- =>: -inybas <file_xybas> -outybas <file_ybas> -inobas <file_xyobas> -configobs <file_o> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble Eof analysis of an',
     $    ' ensemble of data section vectors with the eigenspace',
     $    ' defined by an observation ensemble (inobas).'
      tabvalorder(jmod,jaction)=72
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=66
      taborder(jmod,jaction,3)=67
      taborder(jmod,jaction,4)=52
      taborder(jmod,jaction,5)=70
      taborder(jmod,jaction,6)=57
      taborder(jmod,jaction,7)=59
C -XII.7- =>: -inybas <file_xybas> -outybas <file_ybas> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=7
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble Eof analysis of an',
     $    ' ensemble of data section vectors with the eigenspace',
     $    ' defined by a local data section ensemble (inzbas).'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=66
      taborder(jmod,jaction,3)=69
      taborder(jmod,jaction,4)=70
      taborder(jmod,jaction,5)=57
      taborder(jmod,jaction,6)=59
C -XII.8- =>: -inybas <file_xybas> -outybas <file_ybas> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z> -configobs <file_o>
      jaction=8
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble Eof analysis of an',
     $    ' ensemble of data section vectors with the eigenspace',
     $    ' defined by a local data section ensemble (inzbas) at',
     $    ' observation locations (configobs).'
      tabvalorder(jmod,jaction)=72
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=66
      taborder(jmod,jaction,3)=69
      taborder(jmod,jaction,4)=70
      taborder(jmod,jaction,5)=57
      taborder(jmod,jaction,6)=59
      taborder(jmod,jaction,7)=52
C -XII.- =>: -inobas <file_xyobas> -outobas <file_obas> -configobs <file_o> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
C      jaction=
C      WRITE(tabhelporder(jmod,jaction),*) 'Bubble Eof analysis of an',
C     $    ' ensemble of observation vectors.'
C      tabvalorder(jmod,jaction)=72
C      taborder(jmod,jaction,1)=67
C      taborder(jmod,jaction,2)=68
C      taborder(jmod,jaction,3)=52
C      taborder(jmod,jaction,4)=69
C      taborder(jmod,jaction,5)=70
C      taborder(jmod,jaction,6)=57
C      taborder(jmod,jaction,7)=59
C
C -XIII.- ZONE module : Local data section management
C ---------------------------------------------------
C
      jmod=13
C --- optional switches -reducevar(54)
      tabnumswiopt(jmod,1)=54
C -XIII.1- =>: -outdta <file_y> -inzon <file_z> -zonindex <jzonindex>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Extract in a dta file',
     $    ' a local data section related to particular index.'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=21
      taborder(jmod,jaction,2)=59
      taborder(jmod,jaction,3)=63
C -XIII.2- =>: -outpartvar <file_x> -outzon <file_z> -incfg <incfg.cfg>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Generate a partition and',
     $    ' the corresponding influence bubbles from user parameter file.'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=58
      taborder(jmod,jaction,2)=60
      taborder(jmod,jaction,3)=64
C -XIII.3- =>: -inzon <file_z> -incfg <incfg.cfg> -outzbas <file_zbas>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Generate local correlation',
     $    ' matrices (outzbas) from user parameter file (incfg)',
     $    ' and predefined local data section configuration (inzon).'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=59
      taborder(jmod,jaction,2)=64
      taborder(jmod,jaction,3)=70
C -XIII.4- =>: -incfg <incfg.cfg> -outpartvar <file_x>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Generate a partition from a',
     $    ' set of contours (incfg).'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=64
      taborder(jmod,jaction,2)=58
C -XIII.5- =>: -incfg <incfg.cfg> -outvar <file_x>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Build a state vector',
     $    ' (outvar) using a set of elementary functions (incfg).'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=64
      taborder(jmod,jaction,2)=23
C -XIII.6- =>: -incfg <incfg.cfg> -outzon <file_z>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Build a set of influence',
     $    ' bubbles (outzon) from a set of contours (incfg).'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=64
      taborder(jmod,jaction,2)=60
C -XIII.7- =>: -inzon <file_z> -outptzon <file_z>
      jaction=7
      WRITE(tabhelporder(jmod,jaction),*) 'Generate the pointers',
     $    ' (outptzon) corresponding to a local data section'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=59
      taborder(jmod,jaction,2)=79
C -XIII.8- =>: -outvar <file_x> -inrz <file_rz> -inpartvar <file_x> -incfg <incfg.cfg>
      jaction=8
      WRITE(tabhelporder(jmod,jaction),*) 'Extract in a var file',
     $    ' a parameter stored in rz file.'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=89
      taborder(jmod,jaction,3)=57
      taborder(jmod,jaction,4)=64
C
C -XIV.- GROA module : Global reduced order analysis 
C --------------------------------------------------
C
      jmod=14
C --- optional switches -weight(6) -oestd(7) -bias(8) -reducevar(54) -coefrmax(72) -disable(73) -inparadap(76) -outrz(32) -inpartobs(50) -insmocfg(87) -outsmocfg(88)
      tabnumswiopt(jmod,1)=6
      tabnumswiopt(jmod,2)=7
      tabnumswiopt(jmod,3)=8
      tabnumswiopt(jmod,4)=54
      tabnumswiopt(jmod,5)=72
      tabnumswiopt(jmod,6)=73
      tabnumswiopt(jmod,7)=76
      tabnumswiopt(jmod,8)=32
      tabnumswiopt(jmod,9)=50
      tabnumswiopt(jmod,10)=87
      tabnumswiopt(jmod,11)=88
C -XIV.1- =>: -outvar <file_x> -invarref <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -indta <file_y>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Global reduced order',
     $    ' analysis using a global Kxy gain'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=19
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=20
      taborder(jmod,jaction,5)=12
C -XIV.2- =>: -outvar <file_x> -invarref <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -inobs <file_yo> -configobs <file_o>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Global reduced order',
     $    ' analysis using a global Kxo gain'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=19
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=20
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
C -XIV.3- =>: -outvar <file_x> -invarref <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -inobs <file_yo> -configobs <file_o> -inobsref <file_xyo>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Global reduced order',
     $    ' analysis using a global Kxo gain with a user-defined',
     $    ' observation equivalent for the first guess (inobsref)'
      tabvalorder(jmod,jaction)=72
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=19
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=20
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=18
C -XIV.4- =>: -outvar <file_x> -invarref <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -inobs <file_yo> -configobs <file_o> -inobsref <file_xyo> -inobas <file_xyobas> -outobas <file_obas>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Global reduced order',
     $    ' analysis using a global Kxo gain with a user-defined',
     $    ' observation equivalent for the first guess (inobsref) and',
     $    ' the background error covariance (inobas)'
      tabvalorder(jmod,jaction)=92
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=19
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=20
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=18
      taborder(jmod,jaction,8)=67
      taborder(jmod,jaction,9)=68
C -XIV.5- =>: -outdta <file_y> -indtaref <file_x> -inybas <file_xybas> -outybas <file_ybas> -indta <file_y>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Global reduced order',
     $    ' analysis using a global Kyy gain'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=21
      taborder(jmod,jaction,2)=17
      taborder(jmod,jaction,3)=65
      taborder(jmod,jaction,4)=66
      taborder(jmod,jaction,5)=12
C -XIV.6- =>: -outdta <file_y> -indtaref <file_x> -inybas <file_xybas> -outybas <file_ybas> -inobs <file_yo> -configobs <file_o>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Global reduced order',
     $    ' analysis using a global Kyo gain'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=21
      taborder(jmod,jaction,2)=17
      taborder(jmod,jaction,3)=65
      taborder(jmod,jaction,4)=66
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
C -XIV.7- => -outdta <file_y> -indtaref <file_x> -inybas <file_xybas> -outybas <file_ybas> -inobs <file_yo> -configobs <file_o> -inobsref <file_xyo>
      jaction=7
      WRITE(tabhelporder(jmod,jaction),*) 'Global reduced order',
     $    ' analysis using a global Kyo gain with a user-defined',
     $    ' observation equivalent for the first guess (inobsref)'
      tabvalorder(jmod,jaction)=72
      taborder(jmod,jaction,1)=21
      taborder(jmod,jaction,2)=17
      taborder(jmod,jaction,3)=65
      taborder(jmod,jaction,4)=66
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=18
C -XIV.8- =>: -outdta <file_y> -indtaref <file_x> -inybas <file_xybas> -outybas <file_ybas> -inobs <file_yo> -configobs <file_o> -inobsref <file_xyo> -inobas <file_xyobas> -outobas <file_obas>
      jaction=8
      WRITE(tabhelporder(jmod,jaction),*) 'Global reduced order',
     $    ' analysis using a global Kyo gain with a user-defined',
     $    ' observation equivalent for the first guess (inobsref) and',
     $    ' the background error covariance (inobas)'
      tabvalorder(jmod,jaction)=92
      taborder(jmod,jaction,1)=21
      taborder(jmod,jaction,2)=17
      taborder(jmod,jaction,3)=65
      taborder(jmod,jaction,4)=66
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=18
      taborder(jmod,jaction,8)=67
      taborder(jmod,jaction,9)=68
C -XIV.9- =>: -outobs <file_o> -inobsref <file_xyo> -inobas <file_xyobas> -outobas <file_obas> -inobs <file_yo> -configobs <file_o>
      jaction=9
      WRITE(tabhelporder(jmod,jaction),*) 'Global reduced order',
     $    ' analysis using a global Koo gain'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=22
      taborder(jmod,jaction,2)=18
      taborder(jmod,jaction,3)=67
      taborder(jmod,jaction,4)=68
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
C
C -XV.- LROA module : Local reduced order analysis
C ------------------------------------------------
C
      jmod=15
C --- optional switches -weight(6) -oestd(7) -bias(8) -reducevar(54) -oecorrel(61) -coefrmax(72) -disable(73) -inparadap(76) -inptzon(78) -outrz(32) -insmocfg(87) -outsmocfg(88) -typedtadiag(34)
      tabnumswiopt(jmod,1)=6
      tabnumswiopt(jmod,2)=7
      tabnumswiopt(jmod,3)=8
      tabnumswiopt(jmod,4)=54
      tabnumswiopt(jmod,5)=61
      tabnumswiopt(jmod,6)=72
      tabnumswiopt(jmod,7)=73
      tabnumswiopt(jmod,8)=76
      tabnumswiopt(jmod,9)=78
      tabnumswiopt(jmod,10)=32
      tabnumswiopt(jmod,11)=50
      tabnumswiopt(jmod,12)=87
      tabnumswiopt(jmod,13)=88
      tabnumswiopt(jmod,14)=34
C -XV.1- =>: -outvar <file_x> -invarref <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -indta <file_y> -inpartvar <file_x> -inzon <file_z>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Local reduced order',
     $    ' analysis using a local Kxy gain'
      tabvalorder(jmod,jaction)=72
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=19
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=20
      taborder(jmod,jaction,5)=12
      taborder(jmod,jaction,6)=57
      taborder(jmod,jaction,7)=59
C -XV.- =>: -outvar <file_x> -invarref <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -inobs <file_yo> -configobs <file_o> -inpartvar <file_x> -inzon <file_z>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Local reduced order',
     $    ' analysis using a local Kxo gain'
      tabvalorder(jmod,jaction)=82
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=19
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=20
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=57
      taborder(jmod,jaction,8)=59
C -XV.3- =>: -outvar <file_x> -invarref <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -inobs <file_yo> -configobs <file_o> -inobsref <file_xyo> -inpartvar <file_x> -inzon <file_z>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Local reduced order',
     $    ' analysis using a local Kxo gain with a user-defined',
     $    ' observation equivalent for the first guess (inobsref)'
      tabvalorder(jmod,jaction)=92
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=19
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=20
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=18
      taborder(jmod,jaction,8)=57
      taborder(jmod,jaction,9)=59
C -XV.4- =>: -outvar <file_x> -invarref <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -inobs <file_yo> -configobs <file_o> -inobsref <file_xyo> -inobas <file_xyobas> -outobas <file_obas> -inpartvar <file_x> -inzon <file_z>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Local reduced order',
     $    ' analysis using a local Kxo gain with a user-defined',
     $    ' observation equivalent for the first guess (inobsref) and',
     $    ' the background error covariance (inobas)'
      tabvalorder(jmod,jaction)=112
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=19
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=20
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=18
      taborder(jmod,jaction,8)=67
      taborder(jmod,jaction,9)=68
      taborder(jmod,jaction,10)=57
      taborder(jmod,jaction,11)=59
C -XV.5- =>: -outdta <file_y> -indtaref <file_x> -inybas <file_xybas> -outybas <file_ybas> -indta <file_y> -inpartvar <file_x> -inzon <file_z>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Local reduced order',
     $    ' analysis using a local Kyy gain'
      tabvalorder(jmod,jaction)=72
      taborder(jmod,jaction,1)=21
      taborder(jmod,jaction,2)=17
      taborder(jmod,jaction,3)=65
      taborder(jmod,jaction,4)=66
      taborder(jmod,jaction,5)=12
      taborder(jmod,jaction,6)=57
      taborder(jmod,jaction,7)=59
C -XV.6- =>: -outdta <file_y> -indtaref <file_x> -inybas <file_xybas> -outybas <file_ybas> -inobs <file_yo> -configobs <file_o> -inpartvar <file_x> -inzon <file_z>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Local reduced order',
     $    ' analysis using a local Kyo gain'
      tabvalorder(jmod,jaction)=82
      taborder(jmod,jaction,1)=21
      taborder(jmod,jaction,2)=17
      taborder(jmod,jaction,3)=65
      taborder(jmod,jaction,4)=66
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=57
      taborder(jmod,jaction,8)=59
C -XV.7- =>: -outdta <file_y> -indtaref <file_x> -inybas <file_xybas> -outybas <file_ybas> -inobs <file_yo> -configobs <file_o> -inobsref <file_xyo> -inpartvar <file_x> -inzon <file_z>
      jaction=7
      WRITE(tabhelporder(jmod,jaction),*) 'Local reduced order',
     $    ' analysis using a local Kyo gain with a user-defined',
     $    ' observation equivalent for the first guess (inobsref)'
      tabvalorder(jmod,jaction)=92
      taborder(jmod,jaction,1)=21
      taborder(jmod,jaction,2)=17
      taborder(jmod,jaction,3)=65
      taborder(jmod,jaction,4)=66
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=18
      taborder(jmod,jaction,8)=57
      taborder(jmod,jaction,9)=59
C -XV.8- =>: -outdta <file_y> -indtaref <file_x> -inybas <file_xybas> -outybas <file_ybas> -inobs <file_yo> -configobs <file_o> -inobsref <file_xyo> -inobas <file_xyobas> -outobas <file_obas> -inpartvar <file_x> -inzon <file_z>
      jaction=8
      WRITE(tabhelporder(jmod,jaction),*) 'Local reduced order',
     $    ' analysis using a local Kyo gain with a user-defined',
     $    ' observation equivalent for the first guess (inobsref) and',
     $    ' the background error covariance (inobas)'
      tabvalorder(jmod,jaction)=112
      taborder(jmod,jaction,1)=21
      taborder(jmod,jaction,2)=17
      taborder(jmod,jaction,3)=65
      taborder(jmod,jaction,4)=66
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=18
      taborder(jmod,jaction,8)=67
      taborder(jmod,jaction,9)=68
      taborder(jmod,jaction,10)=57
      taborder(jmod,jaction,11)=59
C -XV.9- =>: -outobs <file_o> -inobsref <file_xyo> -inobas <file_xyobas> -outobas <file_obas> -inobs <file_yo> -configobs <file_o> -inpartvar <file_x> -inzon <file_z>
      jaction=9
      WRITE(tabhelporder(jmod,jaction),*) 'Local reduced order',
     $    ' analysis using a local Koo gain'
      tabvalorder(jmod,jaction)=82
      taborder(jmod,jaction,1)=22
      taborder(jmod,jaction,2)=18
      taborder(jmod,jaction,3)=67
      taborder(jmod,jaction,4)=68
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=57
      taborder(jmod,jaction,8)=59
C
C -XVI.- BROA module : Bubble reduced order analysis
C --------------------------------------------------
C
      jmod=16
C --- optional switches -weight(6) -oestd(7) -bias(8) -reducevar(54) -oecorrel(61) -coefrmax(72) -disable(73) -inparadap(76) -inptzon(78) -outrz(32)
      tabnumswiopt(jmod,1)=6
      tabnumswiopt(jmod,2)=7
      tabnumswiopt(jmod,3)=8
      tabnumswiopt(jmod,4)=54
      tabnumswiopt(jmod,5)=61
      tabnumswiopt(jmod,6)=72
      tabnumswiopt(jmod,7)=73
      tabnumswiopt(jmod,8)=76
      tabnumswiopt(jmod,9)=78
      tabnumswiopt(jmod,10)=32
      tabnumswiopt(jmod,11)=50
C -XVI.1- =>: -outvar <file_x> -invarref <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -indta <file_y> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble reduced order',
     $    ' analysis using a local Kxy gain with a',
     $    ' background error covariance by bubble (inzbas)',
     $    ' and subdomain (inxbas)'
      tabvalorder(jmod,jaction)=92
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=19
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=20
      taborder(jmod,jaction,5)=12
      taborder(jmod,jaction,6)=69
      taborder(jmod,jaction,7)=70
      taborder(jmod,jaction,8)=57
      taborder(jmod,jaction,9)=59
C -XVI.2- =>: -outvar <file_x> -invarref <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -inobs <file_yo> -configobs <file_o> -inobsref <file_xyo> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble reduced order',
     $    ' analysis using a local Kxo gain with a',
     $    ' background error covariance by bubble (inzbas)',
     $    ' and subdomain (inxbas) and with a user-defined',
     $    ' observation equivalent for the first guess (inobsref)'
      tabvalorder(jmod,jaction)=112
      taborder(jmod,jaction,1)=23
      taborder(jmod,jaction,2)=19
      taborder(jmod,jaction,3)=10
      taborder(jmod,jaction,4)=20
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=18
      taborder(jmod,jaction,8)=69
      taborder(jmod,jaction,9)=70
      taborder(jmod,jaction,10)=57
      taborder(jmod,jaction,11)=59
C -XVI.3- =>: -outdta <file_y> -indtaref <file_x> -inybas <file_xybas> -outybas <file_ybas> -indta <file_y> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble reduced order',
     $    ' analysis using a local Kyy gain with a',
     $    ' background error covariance by bubble (inzbas)',
     $    ' and subdomain (inybas)'
      tabvalorder(jmod,jaction)=92
      taborder(jmod,jaction,1)=21
      taborder(jmod,jaction,2)=17
      taborder(jmod,jaction,3)=65
      taborder(jmod,jaction,4)=66
      taborder(jmod,jaction,5)=12
      taborder(jmod,jaction,6)=69
      taborder(jmod,jaction,7)=70
      taborder(jmod,jaction,8)=57
      taborder(jmod,jaction,9)=59
C -XVI.4- =>: -outdta <file_y> -indtaref <file_x> -inybas <file_xybas> -outybas <file_ybas> -inobs <file_yo> -configobs <file_o> -inobsref <file_xyo> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble reduced order',
     $    ' analysis using a local Kyo gain with a',
     $    ' background error covariance by bubble (inzbas)',
     $    ' and subdomain (inybas) and with a user-defined',
     $    ' observation equivalent for the first guess (inobsref)'
      tabvalorder(jmod,jaction)=112
      taborder(jmod,jaction,1)=21
      taborder(jmod,jaction,2)=17
      taborder(jmod,jaction,3)=65
      taborder(jmod,jaction,4)=66
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=18
      taborder(jmod,jaction,8)=69
      taborder(jmod,jaction,9)=70
      taborder(jmod,jaction,10)=57
      taborder(jmod,jaction,11)=59
C -XVI.5- =>: -outobs <file_o> -inobsref <file_xyo> -inobas <file_xyobas> -outobas <file_obas> -inobs <file_yo> -configobs <file_o> -inzbas <file_zbas> -outzbas <file_zbas>  -inpartvar <file_x> -inzon <file_z>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble reduced order',
     $    ' analysis using a local Koo gain with a',
     $    ' background error covariance by bubble (inzbas)',
     $    ' and subdomain (inobas) and with a user-defined',
     $    ' observation equivalent for the first guess (inobsref)'
      tabvalorder(jmod,jaction)=102
      taborder(jmod,jaction,1)=22
      taborder(jmod,jaction,2)=18
      taborder(jmod,jaction,3)=67
      taborder(jmod,jaction,4)=68
      taborder(jmod,jaction,5)=13
      taborder(jmod,jaction,6)=52
      taborder(jmod,jaction,7)=69
      taborder(jmod,jaction,8)=70
      taborder(jmod,jaction,9)=57
      taborder(jmod,jaction,10)=59
C
C -XVII.- GREG module : Global linear regression
C ----------------------------------------------
      jmod=17
C --- optional switches -weight(6) -oestd(7) -reducevar(54) -disable(73)
      tabnumswiopt(jmod,1)=6
      tabnumswiopt(jmod,2)=7
      tabnumswiopt(jmod,3)=54
      tabnumswiopt(jmod,4)=73
C -XVII.1- =>: -inxbasref <file_xbas> -inxbas <file_xbas> -outxbas <file_xbas>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Global linear regression of',
     $    ' an ensemble of state vectors (inxbasref,jpd) using an',
     $    ' ensemble of state vectors modes (inxbas,jpr).',
     $    ' The action writes the projected',
     $    ' state vectors ensemble (outxbas,jpd).'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=15
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=20
C -XVII.2- =>: -inybasref <file_xybas> -inxbas <file_xbas> -outxbas <file_xbas> -inybas <file_xybas> 
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Global linear regression of',
     $    ' a data section ensemble (inybasref,jpd) using an',
     $    ' ensemble of state vectors modes (inxbas,jpr)',
     $    ' with the eigenspace defined by a',
     $    ' data section modes ensemble (inybas,jpr).',
     $    ' The action writes the projected',
     $    ' state vectors ensemble (outxbas,jpd).'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=45
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=20
      taborder(jmod,jaction,4)=65
C -XVII.3- =>: -inobasref <file_xyobas> -inxbas <file_xbas> -outxbas <file_xbas> -inobas <file_xyobas> -configobs <file_o>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Global linear regression of',
     $    ' an observation ensemble (inobasref,jpd) using an',
     $    ' ensemble of state vectors modes (inxbas,jpr)',
     $    ' with the eigenspace defined by an',
     $    ' observation modes ensemble (inobas,jpr).',
     $    ' The action writes the projected',
     $    ' state vectors ensemble (outxbas,jpd).'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=74
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=20
      taborder(jmod,jaction,4)=67
      taborder(jmod,jaction,5)=52
C -XVII.4- =>: -inobasref <file_xyobas> -inxbas <file_xbas> -outxbas <file_xbas> -inobas <file_xyobas> -configobs <file_o> -outobas <file_obas>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Global linear regression of',
     $    ' an observation ensemble (inobasref,jpd) using an',
     $    ' ensemble of state vectors modes (inxbas,jpr)',
     $    ' with the eigenspace defined by an',
     $    ' observation modes ensemble (inobas,jpr).',
     $    ' The action writes the projected',
     $    ' observation ensemble (outobas,jpd)',
     $    ' and the projected',
     $    ' state vectors ensemble (outxbas,jpd).'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=74
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=20
      taborder(jmod,jaction,4)=67
      taborder(jmod,jaction,5)=52
      taborder(jmod,jaction,6)=68
C -XVII.5- =>: -inybasref <file_xybas> -inybas <file_xybas> -outybas <file_ybas>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Global linear regression of',
     $    ' a data section ensemble (inybasref,jpd) using an',
     $    ' ensemble of data section vectors modes (inybas,jpr).',
     $    ' The action writes the projected',
     $    ' data section vectors ensemble (outybas,jpd).'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=45
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=66
C -XVII.6- =>: -inobasref <file_xyobas> -inybas <file_xybas> -outybas <file_ybas> -inobas <file_xyobas> -configobs <file_o>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Global linear regression of',
     $    ' an observation ensemble (inobasref,jpd) using an',
     $    ' ensemble of data section vectors modes (inybas,jpr)',
     $    ' with the eigenspace defined by an',
     $    ' observation modes ensemble (inobas,jpr).',
     $    ' The action writes the projected',
     $    ' data section vectors ensemble (outybas,jpd).'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=74
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=66
      taborder(jmod,jaction,4)=67
      taborder(jmod,jaction,5)=52
C -XVII.7- =>: -inobasref <file_xyobas> -inybas <file_xybas> -outybas <file_ybas> -inobas <file_xyobas> -configobs <file_o> -outobas <file_obas>
      jaction=7
      WRITE(tabhelporder(jmod,jaction),*) 'Global linear regression of',
     $    ' an observation ensemble (inobasref,jpd) using an',
     $    ' ensemble of data section vectors modes (inybas,jpr)',
     $    ' with the eigenspace defined by an',
     $    ' observation modes ensemble (inobas,jpr).',
     $    ' The action writes the projected',
     $    ' data section vectors ensemble (outybas,jpd)',
     $    ' and the projected',
     $    ' observation vectors ensemble (outobas,jpd).'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=74
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=66
      taborder(jmod,jaction,4)=67
      taborder(jmod,jaction,5)=52
      taborder(jmod,jaction,6)=68
C -XVII.8- =>: -inobasref <file_xyobas> -inobas <file_xyobas> -outobas <file_obas> -configobs <file_o>
      jaction=8
      WRITE(tabhelporder(jmod,jaction),*) 'Global linear regression of',
     $    ' an observation ensemble (inobasref,jpd) using an',
     $    ' ensemble of observation vectors modes (inobas,jpr).',
     $    ' The action writes the projected',
     $    ' observation vectors ensemble (outobas,jpd).'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=74
      taborder(jmod,jaction,2)=67
      taborder(jmod,jaction,3)=68
      taborder(jmod,jaction,4)=52
C
C -XVIII.- LREG module : Local linear regression
C ----------------------------------------------
      jmod=18
C --- optional switches -weight(6) -oestd(7) -reducevar(54) -oecorrel(61) -disable(73)
      tabnumswiopt(jmod,1)=6
      tabnumswiopt(jmod,2)=7
      tabnumswiopt(jmod,3)=54
      tabnumswiopt(jmod,4)=61
      tabnumswiopt(jmod,5)=73
C -XVIII.1- =>: -inybasref <file_xybas> -inxbas <file_xbas> -outxbas <file_xbas> -inybas <file_xybas> -inpartvar <file_x> -inzon <file_z>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Local linear regression',
     $    ' (bubble by bubble) of',
     $    ' an ensemble of data section vectors',
     $    ' (inybasref,jpd) using an',
     $    ' ensemble of state vectors modes (inxbas,jpr)',
     $    ' with the eigenspace defined by a',
     $    ' data section modes ensemble (inybas,jpr).',
     $    ' The action writes the projected',
     $    ' state vectors ensemble (outxbas,jpd).'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=45
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=20
      taborder(jmod,jaction,4)=65
      taborder(jmod,jaction,5)=57
      taborder(jmod,jaction,6)=59
C -XVIII.2- =>: -inobasref <file_xyobas> -inxbas <file_xbas> -outxbas <file_xbas> -inobas <file_xyobas> -configobs <file_o> -inpartvar <file_x> -inzon <file_z>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Local linear regression',
     $    ' (bubble by bubble) of',
     $    ' an ensemble of observation vectors',
     $    ' (inobasref,jpd) using an',
     $    ' ensemble of state vectors modes (inxbas,jpr)',
     $    ' with the eigenspace defined by an',
     $    ' observation modes ensemble (inobas,jpr).',
     $    ' The action writes the projected',
     $    ' state vectors ensemble (outxbas,jpd).'
      tabvalorder(jmod,jaction)=72
      taborder(jmod,jaction,1)=74
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=20
      taborder(jmod,jaction,4)=67
      taborder(jmod,jaction,5)=52
      taborder(jmod,jaction,6)=57
      taborder(jmod,jaction,7)=59
C -XVIII.- =>: -inobasref <file_xyobas> -inxbas <file_xbas> -outxbas <file_xbas> -inobas <file_xyobas> -configobs <file_o> -outobas <file_obas> -inpartvar <file_x> -inzon <file_z>
C      jaction=
C      tabvalorder(jmod,jaction)=82
C      taborder(jmod,jaction,1)=74
C      taborder(jmod,jaction,2)=10
C      taborder(jmod,jaction,3)=20
C      taborder(jmod,jaction,4)=67
C      taborder(jmod,jaction,5)=52
C      taborder(jmod,jaction,6)=68
C      taborder(jmod,jaction,7)=57
C      taborder(jmod,jaction,8)=59
C -XVIII.3- =>: -inybasref <file_xybas> -inybas <file_xybas> -outybas <file_ybas> -inpartvar <file_x> -inzon <file_z>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Local linear regression',
     $    ' (bubble by bubble) of',
     $    ' an ensemble of data section vectors',
     $    ' (inybasref,jpd) using an',
     $    ' ensemble of data section vectors modes (inybas,jpr).',
     $    ' The action writes the projected',
     $    ' data section vectors ensemble (outybas,jpd).'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=45
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=66
      taborder(jmod,jaction,4)=57
      taborder(jmod,jaction,5)=59
C -XVIII.4- =>: -inobasref <file_xyobas> -inybas <file_xybas> -outybas <file_ybas> -inobas <file_xyobas> -configobs <file_o> -inpartvar <file_x> -inzon <file_z>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Local linear regression',
     $    ' (bubble by bubble) of',
     $    ' an ensemble of observation vectors',
     $    ' (inobasref,jpd) using an',
     $    ' ensemble of data section vectors modes (inybas,jpr)',
     $    ' with the eigenspace defined by an',
     $    ' observation modes ensemble (inobas,jpr).',
     $    ' The action writes the projected',
     $    ' data section vectors ensemble (outybas,jpd)'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=74
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=66
      taborder(jmod,jaction,4)=67
      taborder(jmod,jaction,5)=52
      taborder(jmod,jaction,6)=57
      taborder(jmod,jaction,7)=59
C -XVIII.- =>: -inobasref <file_xyobas> -inybas <file_xybas> -outybas <file_ybas> -inobas <file_xyobas> -configobs <file_o> -outobas <file_obas> -inpartvar <file_x> -inzon <file_z>
C      jaction=
C      tabvalorder(jmod,jaction)=82
C      taborder(jmod,jaction,1)=74
C      taborder(jmod,jaction,2)=65
C      taborder(jmod,jaction,3)=66
C      taborder(jmod,jaction,4)=67
C      taborder(jmod,jaction,5)=52
C      taborder(jmod,jaction,6)=68
C      taborder(jmod,jaction,7)=57
C      taborder(jmod,jaction,8)=59
C -XVIII.- =>: -inobasref <file_xyobas> -inobas <file_xyobas> -outobas <file_obas> -configobs <file_o> -inpartvar <file_x> -inzon <file_z>
C      jaction=
C      tabvalorder(jmod,jaction)=62
C      taborder(jmod,jaction,1)=74
C      taborder(jmod,jaction,2)=67
C      taborder(jmod,jaction,3)=68
C      taborder(jmod,jaction,4)=52
C      taborder(jmod,jaction,5)=57
C      taborder(jmod,jaction,6)=59
C
C -XIX.- BREG module : Bubble linear regression
C ---------------------------------------------
      jmod=19
C --- optional switches -weight(6) -oestd(7) -reducevar(54) -oecorrel(61) -disable(73)
      tabnumswiopt(jmod,1)=6
      tabnumswiopt(jmod,2)=7
      tabnumswiopt(jmod,3)=54
      tabnumswiopt(jmod,4)=61
      tabnumswiopt(jmod,5)=73
C -XIX.1- =>: -inybasref <file_xybas> -inxbas <file_xbas> -outxbas <file_xbas> -inybas <file_xybas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble linear regression',
     $    ' (bubble by bubble) of',
     $    ' an ensemble of data section vectors (inybasref,jpd) using an',
     $    ' ensemble of state vectors modes (inxbas,jpr)',
     $    ' with the eigenspace defined by a',
     $    ' data section modes ensemble (inybas,jpr).',
     $    ' The action writes the projected',
     $    ' state vectors ensemble (outxbas,jpd),',
     $    ' and the projected',
     $    ' local data section vectors ensemble (outzbas,jpd).'
      tabvalorder(jmod,jaction)=72
      taborder(jmod,jaction,1)=45
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=20
      taborder(jmod,jaction,4)=65
      taborder(jmod,jaction,5)=70
      taborder(jmod,jaction,6)=57
      taborder(jmod,jaction,7)=59
C -XIX.2- =>: -inobasref <file_xyobas> -inxbas <file_xbas> -outxbas <file_xbas> -inobas <file_xyobas> -configobs <file_o> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble linear regression',
     $    ' (bubble by bubble) of',
     $    ' an ensemble of observation vectors (inobasref,jpd) using an',
     $    ' ensemble of state vectors modes (inxbas,jpr)',
     $    ' with the eigenspace defined by an',
     $    ' observation modes ensemble (inobas,jpr).',
     $    ' The action writes the projected',
     $    ' state vectors ensemble (outxbas,jpd),',
     $    ' and the projected',
     $    ' local data section vectors ensemble (outzbas,jpd).'
      tabvalorder(jmod,jaction)=82
      taborder(jmod,jaction,1)=74
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=20
      taborder(jmod,jaction,4)=67
      taborder(jmod,jaction,5)=52
      taborder(jmod,jaction,6)=70
      taborder(jmod,jaction,7)=57
      taborder(jmod,jaction,8)=59
C -XIX.3- =>: -inzbasref <file_zbas> -inxbas <file_xbas> -outxbas <file_xbas> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble linear regression',
     $    ' (bubble by bubble) of',
     $    ' an ensemble of local data section vectors (inzbasref,jpd)',
     $    ' using an',
     $    ' ensemble of state vectors modes (inxbas,jpr)',
     $    ' with the eigenspace defined by a',
     $    ' local data section modes ensemble (inzbas,jpr).',
     $    ' The action writes the projected',
     $    ' state vectors ensemble (outxbas,jpd),',
     $    ' and the projected',
     $    ' local data section vectors ensemble (outzbas,jpd).'
      tabvalorder(jmod,jaction)=72
      taborder(jmod,jaction,1)=81
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=20
      taborder(jmod,jaction,4)=69
      taborder(jmod,jaction,5)=70
      taborder(jmod,jaction,6)=57
      taborder(jmod,jaction,7)=59
C -XIX.4- =>: -inzbasref <file_zbas> -inxbas <file_xbas> -outxbas <file_xbas> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z> -configobs <file_o>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble linear regression',
     $    ' (bubble by bubble) of',
     $    ' an ensemble of local data section vectors (inzbasref,jpd)',
     $    ' at observation locations (configobs) using an',
     $    ' ensemble of state vectors modes (inxbas,jpr)',
     $    ' with the eigenspace defined by a',
     $    ' local data section modes ensemble (inzbas,jpr)',
     $    ' at observation locations (configobs).',
     $    ' The action writes the projected',
     $    ' state vectors ensemble (outxbas,jpd),',
     $    ' and the projected',
     $    ' local data section vectors ensemble (outzbas,jpd).'
      tabvalorder(jmod,jaction)=82
      taborder(jmod,jaction,1)=81
      taborder(jmod,jaction,2)=10
      taborder(jmod,jaction,3)=20
      taborder(jmod,jaction,4)=69
      taborder(jmod,jaction,5)=70
      taborder(jmod,jaction,6)=57
      taborder(jmod,jaction,7)=59
      taborder(jmod,jaction,8)=52
C -XIX.5- =>: -inybasref <file_xybas> -inybas <file_xybas> -outybas <file_ybas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble linear regression',
     $    ' (bubble by bubble) of',
     $    ' an ensemble of data section vectors (inybasref,jpd)',
     $    ' using an',
     $    ' ensemble of data section vectors modes (inybas,jpr).',
     $    ' The action writes the projected',
     $    ' data section vectors ensemble (outybas,jpd),',
     $    ' and the projected',
     $    ' local data section vectors ensemble (outzbas,jpd).'
      tabvalorder(jmod,jaction)=62
      taborder(jmod,jaction,1)=45
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=66
      taborder(jmod,jaction,4)=70
      taborder(jmod,jaction,5)=57
      taborder(jmod,jaction,6)=59
C -XIX.6- =>: -inobasref <file_xyobas> -inybas <file_xybas> -outybas <file_ybas> -inobas <file_xyobas> -configobs <file_o> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble linear regression',
     $    ' (bubble by bubble) of',
     $    ' an ensemble of data section vectors (inobasref,jpd)',
     $    ' using an',
     $    ' ensemble of data section vectors modes (inybas,jpr)',
     $    ' with the eigenspace defined by an',
     $    ' observation modes ensemble (inobas,jpr).',
     $    ' The action writes the projected',
     $    ' data section vectors ensemble (outybas,jpd),',
     $    ' and the projected',
     $    ' local data section vectors ensemble (outzbas,jpd).'
      tabvalorder(jmod,jaction)=82
      taborder(jmod,jaction,1)=74
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=66
      taborder(jmod,jaction,4)=67
      taborder(jmod,jaction,5)=52
      taborder(jmod,jaction,6)=70
      taborder(jmod,jaction,7)=57
      taborder(jmod,jaction,8)=59
C -XIX.7- =>: -inybasref <file_xybas> -inybas <file_xybas> -outybas <file_ybas> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
      jaction=7
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble linear regression',
     $    ' (bubble by bubble) of',
     $    ' an ensemble of local data section vectors (inybasref,jpd)',
     $    ' using an',
     $    ' ensemble of data section vectors modes (inybas,jpr)',
     $    ' with the eigenspace defined by an',
     $    ' local data section modes ensemble (inzbas,jpr).',
     $    ' The action writes the projected',
     $    ' data section vectors ensemble (outybas,jpd),',
     $    ' and the projected',
     $    ' local data section vectors ensemble (outzbas,jpd).'
      tabvalorder(jmod,jaction)=72
      taborder(jmod,jaction,1)=45
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=66
      taborder(jmod,jaction,4)=69
      taborder(jmod,jaction,5)=70
      taborder(jmod,jaction,6)=57
      taborder(jmod,jaction,7)=59
C -XIX.8- =>: -inzbasref <file_zbas> -inybas <file_xybas> -outybas <file_ybas> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z> -configobs <file_o>
      jaction=8
      WRITE(tabhelporder(jmod,jaction),*) 'Bubble linear regression',
     $    ' (bubble by bubble) of',
     $    ' an ensemble of local data section vectors (inzbasref,jpd)',
     $    ' at observation locations (configobs) using an',
     $    ' ensemble of data section vectors modes (inybas,jpr)',
     $    ' with the eigenspace defined by an',
     $    ' local data section modes ensemble (inzbas,jpr)',
     $    ' at observation locations (configobs).',
     $    ' The action writes the projected',
     $    ' data section vectors ensemble (outybas,jpd),',
     $    ' and the projected',
     $    ' local data section vectors ensemble (outzbas,jpd).'
      tabvalorder(jmod,jaction)=82
      taborder(jmod,jaction,1)=81
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=66
      taborder(jmod,jaction,4)=69
      taborder(jmod,jaction,5)=70
      taborder(jmod,jaction,6)=57
      taborder(jmod,jaction,7)=59
      taborder(jmod,jaction,8)=52
C -XIX.- =>: -inzbasref <file_zbas> -inobas <file_xyobas> -outobas <file_obas> -configobs <file_o> -inzbas <file_zbas> -outzbas <file_zbas> -inpartvar <file_x> -inzon <file_z>
C      jaction=
C      tabvalorder(jmod,jaction)=72
C      taborder(jmod,jaction,1)=81
C      taborder(jmod,jaction,2)=67
C      taborder(jmod,jaction,3)=68
C      taborder(jmod,jaction,4)=52
C      taborder(jmod,jaction,5)=69
C      taborder(jmod,jaction,6)=70
C      taborder(jmod,jaction,7)=57
C      taborder(jmod,jaction,8)=59
C
C -XX.- VARI module : Compute variances
C -------------------------------------
      jmod=20
C -II.1- =>: -inxbas <file_xbas> -outvar <file_x>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Compute variances',
     $   ' (outvar) from state reduced order',
     $   ' covariance matrix (inxbas)'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=23
C -II.2- =>: -inybas <file_ybas> -outdta <file_y>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Compute variances',
     $   ' (outdta) from data section reduced order',
     $   ' covariance matrix (inybas)'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=21
C -II.3- =>: -inobas <file_obas> -outobs <file_o> -configobs <file_o>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Compute variances',
     $   ' (outobs) from observation reduced order',
     $   ' covariance matrix (inobas)'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=22
      taborder(jmod,jaction,3)=52
C
C -XXI.- ANAM module : Anamorphosis operations
C --------------------------------------------
      jmod=21
C --- optional switches -oestd(7)
      tabnumswiopt(jmod,1)=7
C -XXI.1- =>: -inxbas <file_xbas> -outxbasref <file_xbas>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Compute Cx ensemble percentiles'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=24
C -XXI.2- =>: -inybas <file_xybas> -outybasref <file_ybas>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Compute Cy ensemble percentiles'
      tabvalorder(jmod,jaction)=22
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=46
C -XXI.3- =>: -inobas <file_xyobas> -outobasref <file_obas> -configobs <file_o>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Compute Co ensemble percentiles'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=86
      taborder(jmod,jaction,3)=52
C -XXI.4- =>: -invar <file_x> -inxbasref <file_xbas> -outvar <file_x> -typeoper <operation>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Anamorphosis of state vector'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=15
      taborder(jmod,jaction,3)=23
      taborder(jmod,jaction,4)=33
C -XXI.5- =>: -indta <file_xy> -inybasref <file_xybas> -outdta <file_y> -typeoper <operation>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Anamorphosis of data section vector'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=45
      taborder(jmod,jaction,3)=21
      taborder(jmod,jaction,4)=33
C -XXI.6- =>: -inobs <file_x> -inobasref <file_xbas> -outobs <file_o> -typeoper <operation> -configobs <file_o>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Anamorphosis of observation vector'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=74
      taborder(jmod,jaction,3)=22
      taborder(jmod,jaction,4)=33
      taborder(jmod,jaction,5)=52
C -XXI.7- =>: -inxbas <file_xbas> -inxbasref <file_xbas> -outxbas <file_xbas> -typeoper <operation>
      jaction=7
      WRITE(tabhelporder(jmod,jaction),*) 'Anamorphosis of state ensemble'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=15
      taborder(jmod,jaction,3)=20
      taborder(jmod,jaction,4)=33
C -XXI.8- =>: -inybas <file_xybas> -inybasref <file_xybas> -outybas <file_ybas> -typeoper <operation>
      jaction=8
      WRITE(tabhelporder(jmod,jaction),*) 'Anamorphosis of data section ensemble'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=45
      taborder(jmod,jaction,3)=66
      taborder(jmod,jaction,4)=33
C -XXI.9- =>: -inobas <file_xyobas> -inobasref <file_xyobas> -outobas <file_obas> -typeoper <operation> -configobs <file_o>
      jaction=9
      WRITE(tabhelporder(jmod,jaction),*) 'Anamorphosis of observation ensemble'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=74
      taborder(jmod,jaction,3)=68
      taborder(jmod,jaction,4)=33
      taborder(jmod,jaction,5)=52
C -XXI.10- =>: -indta <file_xy> -inybas <file_xybas> -outybas <file_xybas>
      jaction=10
      WRITE(tabhelporder(jmod,jaction),*) 'Ensemble anamorphosis of data section vector'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=65
      taborder(jmod,jaction,3)=66
C -XXI.11- =>: -inobs <file_xyo> -inobas <file_xyobas> -outobas <file_obas> -configobs <file_o>
      jaction=11
      WRITE(tabhelporder(jmod,jaction),*) 'Ensemble anamorphosis of observation vector'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=67
      taborder(jmod,jaction,3)=68
      taborder(jmod,jaction,4)=52
C
C -XXII.- SCOR module : Probabilistic scores
C ------------------------------------------
      jmod=22
C -XXII.1- =>: -inxbas <file_xbas> -invar <file_x> -typeoper <operation> -inpartvar <file_x>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Compute Cx probabilistic score'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=14
      taborder(jmod,jaction,3)=33
      taborder(jmod,jaction,4)=57
C -XXII.2- =>: -inybas <file_xybas> -indta <file_xy> -typeoper <operation> -inpartvar <file_x>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Compute Cy probabilistic score'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=12
      taborder(jmod,jaction,3)=33
      taborder(jmod,jaction,4)=57
C -XXII.3- =>: -inobas <file_xyobas> -inobs <file_xyo> -configobs <file_o> -typeoper <operation> -inpartvar <file_x>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Compute Co probabilistic score'
      tabvalorder(jmod,jaction)=52
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=13
      taborder(jmod,jaction,3)=52
      taborder(jmod,jaction,4)=33
      taborder(jmod,jaction,5)=57
C
C -XXIII.- SPCT module : Spectral transformation
C ----------------------------------------------
      jmod=23
C -XXIII.1- =>: -invar <file_x> -outvar <file_x> -typeoper <operation>
      jaction=1
      WRITE(tabhelporder(jmod,jaction),*) 'Spectrum of state vector'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=14
      taborder(jmod,jaction,2)=23
      taborder(jmod,jaction,3)=33
C -XXIII.2- =>: -indta <file_xy> -outdta <file_y> -typeoper <operation>
      jaction=2
      WRITE(tabhelporder(jmod,jaction),*) 'Spectrum of data section vector'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=12
      taborder(jmod,jaction,2)=21
      taborder(jmod,jaction,3)=33
C -XXIII.3- =>: -inobs <file_xyo> -configobs <file_o> -outdta <file_y> -typeoper <operation>
      jaction=3
      WRITE(tabhelporder(jmod,jaction),*) 'Spectrum from observation vector'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=13
      taborder(jmod,jaction,2)=52
      taborder(jmod,jaction,3)=21
      taborder(jmod,jaction,4)=33
C -XXIII.4- =>: -inxbas <file_xbas> -outxbas <file_xbas> -typeoper <operation>
      jaction=4
      WRITE(tabhelporder(jmod,jaction),*) 'Spectrum of state ensemble'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=10
      taborder(jmod,jaction,2)=20
      taborder(jmod,jaction,3)=33
C -XXIII.5- =>: -inybas <file_xybas> -outybas <file_ybas> -typeoper <operation>
      jaction=5
      WRITE(tabhelporder(jmod,jaction),*) 'Spectrum of data section ensemble'
      tabvalorder(jmod,jaction)=32
      taborder(jmod,jaction,1)=65
      taborder(jmod,jaction,2)=66
      taborder(jmod,jaction,3)=33
C -XXIII.6- =>: -inobas <file_xyobas> -configobs <file_o> -outybas <file_ybas> -typeoper <operation>
      jaction=6
      WRITE(tabhelporder(jmod,jaction),*) 'Spectrum of observation ensemble'
      tabvalorder(jmod,jaction)=42
      taborder(jmod,jaction,1)=67
      taborder(jmod,jaction,2)=52
      taborder(jmod,jaction,3)=66
      taborder(jmod,jaction,4)=33
