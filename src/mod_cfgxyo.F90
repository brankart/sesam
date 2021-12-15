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
!---------------------------------------------------------------------
!
!                        MODULE MOD_CFGXYO
!
!---------------------------------------------------------------------
!
!  Purpose : Module with SESAM configuration arrays
!  -------   (internal configuration defined in defcst.switch.h, and
!            user configuration defined in defcst.control.h)
!
!  original     : 99/05 (C.E. Testut)
!  modification : 03/02 (J.M. Brankart)
!---------------------------------------------------------------------
#include "config.main.h90"
!
      MODULE mod_cfgxyo
!
      use mod_main
!
      IMPLICIT NONE
!
! -1.- Arrays with SESAM keywords (switches, modules, extensions,...)
! -------------------------------------------------------------------
!
      INTEGER, save :: nbdigits=4
      INTEGER, save :: nallmem,jpfixjpx,jpfixjpz,jpfixjpu,maxiter
      INTEGER, dimension(1:nbxyo), save :: limjpnxyo
      INTEGER, dimension(1:nbztyp), save :: limjpnz,jpnbub,jpbubend
      LOGICAL :: lsplitobs, lsplitstate
!
! Arrays with SESAM switches
! switab  : switch name
! swihelp : help about the nature of the argument of the switch
! swibool : switch is active (TRUE or FALSE)
! swiopt  : switch is optional (0=false, 1=general optional switch,
!                   2=optional switch particular to a few moodules).
      CHARACTER(len=swilg), dimension(1:nbarg), save :: switab
      CHARACTER(len=bgword), dimension(1:nbarg), save :: swihelp
      LOGICAL, dimension(1:nbarg), save :: swibool
      INTEGER, dimension(1:nbarg), save :: swiopt
!
! Arrays with SESAM modules
! modtab  : module name
! modhelp : help about the nature of the module
! modbool : module is active (TRUE or FALSE)
! modtype : type of module (1=analysis modules, 2=diagnostic modules,
!                  3=utility modules)
      CHARACTER(len=modlg), dimension(1:nbmod), save :: modtab
      CHARACTER(len=bgword), dimension(1:nbmod), save :: modhelp
      LOGICAL, dimension(1:nbmod), save :: modbool
      INTEGER, dimension(1:nbmod), save :: modtype
!
! Arrays with SESAM help keywords
! helptab  : help keyword name
! helpbool : help keyword is active (TRUE or FALSE)
      CHARACTER(len=helplg), dimension(1:nbhelp), save :: helptab
      LOGICAL, dimension(1:nbhelp), save :: helpbool
!
! Arrays with SESAM covariance directory extensions
! extbastab  : covariance directory extension name
! extbasbool : covariance directory extension is active (TRUE or FALSE)
! extbasmem  : type of format (unit=as output format, ten=as input format):
!                0=unavailable, 1=available, >1=available with restrictions
      CHARACTER(len=extlg), dimension(1:nbextbas), save :: extbastab
      LOGICAL, dimension(0:nbextbas), save :: extbasbool
      INTEGER, dimension(1:nbextbas), save :: extbasmem
!
! Arrays with SESAM observation database file extensions
! extdbstab  : observation database file extension name
! extdbsbool : observation database file extension is active (TRUE or FALSE)
! extdbsmem  : type of format (unit=as output format, ten=as input format):
!                0=unavailable, 1=available, >1=available with restrictions
      CHARACTER(len=extlg), dimension(1:nbextdbs), save :: extdbstab
      LOGICAL, dimension(0:nbextdbs), save :: extdbsbool
      INTEGER, dimension(1:nbextdbs), save :: extdbsmem
!
! Arrays with SESAM dta file extensions
! extdtatab  : dta file extension name
! extdtabool : dta file extension is active (TRUE or FALSE)
! extdtamem  : type of format (unit=as output format, ten=as input format):
!                0=unavailable, 1=available, >1=available with restrictions
! extdtaunit : different variables are in the same file (TRUE or FALSE)
      CHARACTER(len=extlg), dimension(1:nbextdta), save :: extdtatab
      LOGICAL, dimension(0:nbextdta), save :: extdtabool
      LOGICAL, dimension(1:nbextdta), save :: extdtaunit
      INTEGER, dimension(1:nbextdta), save :: extdtamem
!
! Arrays with SESAM obs file extensions
! extobstab  : obs file extension name
! extobsbool : obs file extension is active (TRUE or FALSE)
! extobsmem  : type of format (unit=as output format, ten=as input format):
!               0=unavailable, 1=available, >1=available with restrictions
! extobsunit : different variables are in the same file (TRUE or FALSE)
      CHARACTER(len=extlg), dimension(1:nbextobs), save :: extobstab
      LOGICAL, dimension(0:nbextobs), save :: extobsbool
      LOGICAL, dimension(1:nbextobs), save :: extobsunit
      INTEGER, dimension(1:nbextobs), save :: extobsmem
!
! Arrays with SESAM var file extensions
! extvartab  : var file extension name
! extvarbool : var file extension is active (TRUE or FALSE)
! extvarmem  : type of format (unit=as output format, ten=as input format):
!                0=unavailable, 1=available, >1=available with restrictions
! extvarunit : different variables are in the same file (TRUE or FALSE)
      CHARACTER(len=extlg), dimension(1:nbextvar), save :: extvartab
      LOGICAL, dimension(0:nbextvar), save :: extvarbool
      LOGICAL, dimension(1:nbextvar), save :: extvarunit
      INTEGER, dimension(1:nbextvar), save :: extvarmem
!
! Arrays with SESAM zon file extensions
! extzontab  : zon file extension name
! extzonbool : zon file extension is active (TRUE or FALSE)
! extzonmem  : type of format (unit=as output format, ten=as input format):
!                0=unavailable, 1=available, >1=available with restrictions
! extzonunit : different variables are in the same file (TRUE or FALSE)
      CHARACTER(len=extlg), dimension(1:nbextzon), save :: extzontab
      LOGICAL, dimension(0:nbextzon), save :: extzonbool
      LOGICAL, dimension(1:nbextzon), save :: extzonunit
      INTEGER, dimension(1:nbextzon), save :: extzonmem
!
! Arrays defining all possible actions for every SESAM modules
! tabvalorder  : Number of switches for action (x10),
!                  and type of action (unit): obsolete feature
! tabhelporder : Help for every actions of each module
! taborder     : List of switch indices defining the action
! tabnumswiopt : List of optional switch indices particular to each module
      INTEGER, dimension(1:nbmod,1:nbaction,1:nborder), save :: taborder
      INTEGER, dimension(1:nbmod,1:nbaction), save :: tabvalorder
      CHARACTER(len=hgword), dimension(1:nbmod,1:nbaction), save ::  &
     &     tabhelporder
      INTEGER, dimension(1:nbmod,1:nbswiopt), save :: tabnumswiopt
!
! -2.- Specific variables for every SESAM arguments
! -------------------------------------------------
! arg<switchname>  : argument of switch <switchname>
! larg<switchname> : argument of the switch exists (TRUE or FALSE)
! --> mode switch
      CHARACTER(len=bgword), save :: argmod
      LOGICAL, save :: largmod
! --> optional switches:
! -help -list -varmsk -dtamsk -weight -oestd -bias -outinfo -fixjpx -fixjpu
! -reducevar -reducedta -scale -biasdbs -oecorrel -fecorrel -coefrmax -disable
! -inrz -outrz -outparadap -fixjpz -inpartobs -anamorphosis -inparadap -action
      CHARACTER(len=bgword), save :: arghelp,arglist,argvarmsk, &
     &     argdtamsk,argweight,argoestd, &
     &     argbias,argoutinfo,argfixjpx,argreducevar,argreducedta, &
     &     argscale,argbiasdbs,argoecorrel,argfecorrel, &
     &     argcoefrmax,argdisable,argoutparadap,arginrz, &
     &     argoutrz,arginpartobs,arganamorphosis,argfixjpz, &
     &     arginparadap,argfixjpu,argaction
      LOGICAL, save :: larghelp,larglist,largvarmsk,largdtamsk, &
     &     largweight,largoestd, &
     &     largbias,largoutinfo,largfixjpx,largreducevar,largreducedta, &
     &     largscale,largbiasdbs,largoecorrel,largfecorrel, &
     &     largcoefrmax,largdisable,largoutparadap,larginrz, &
     &     largoutrz,larginpartobs,larganamorphosis,largfixjpz, &
     &     larginparadap,largfixjpu,largaction
! --> required switches: in<...>
      CHARACTER(len=bgword), save :: arginxbas,argindbs,argindta, &
     &     arginobs,arginvar,arginxbasref,argindbsref,argindtaref, &
     &     arginobsref,arginvarref,arginpartvar,arginzon, &
     &     arginzonref, &
     &     arginybas,arginobas,arginzbas,arginptzon, &
     &     arginybasref,arginobasref,arginzbasref, argincstr
      LOGICAL, save :: larginxbas,largindbs,largindta,larginobs, &
     &     larginvar,larginxbasref,largindbsref,largindtaref, &
     &     larginobsref,larginvarref,larginpartvar,larginzon, &
     &     larginzonref, &
     &     larginybas,larginobas,larginzbas,larginptzon, &
     &     larginybasref,larginobasref,larginzbasref, largincstr
! --> required switches: out<...>
      CHARACTER(len=bgword), save :: argoutxbas,argoutdta,argoutobs, &
     &     argoutvar,argoutxbasref,argoutdtaref,argoutobsref, &
     &     argoutvarref,argoutpartvar,argoutzon, &
     &     argoutybas,argoutobas,argoutzbas,argoutptzon
      LOGICAL, save :: largoutxbas,largoutdta,largoutobs,largoutvar, &
     &     largoutxbasref,largoutdtaref,largoutobsref,largoutvarref, &
     &     largoutpartvar,largoutzon, &
     &     largoutybas,largoutobas,largoutzbas,largoutptzon
! --> other required switches:
      CHARACTER(len=bgword), save :: argtypeoper,argtypedtadiag, &
     &     argdiffobsref,argdiffdtaref,argdiffvarref, &
     &     argdiffobsorg,argdiffdtaorg,argdiffvarorg, &
     &     arginerrdta,argouterrdta,arginstddta,argoutstddta, &
     &     argoutobasref,argoutybasref,argoutbiasdta, &
     &     argaffectobs,argnullobs,argconfigobs,argconfigzon, &
     &     argzonindex,argincfg,arginoptcfg,argconnect, &
     &     arginsmocfg,argoutsmocfg,argiterate
      LOGICAL, save :: largtypeoper,largtypedtadiag, &
     &     largdiffobsref,largdiffdtaref,largdiffvarref, &
     &     largdiffobsorg,largdiffdtaorg,largdiffvarorg, &
     &     larginerrdta,largouterrdta,larginstddta,largoutstddta, &
     &     largoutobasref,largoutybasref,largoutbiasdta, &
     &     largaffectobs,largnullobs,largconfigobs,largconfigzon, &
     &     largzonindex,largincfg,larginoptcfg,largconnect, &
     &     larginsmocfg,largoutsmocfg,largiterate
!
! -3.- Arrays with SESAM user configuration
! -----------------------------------------
! See defcst.control.h for a description of the variables
! that are not described here
!
! Vx object configuration
! varend  : number of variable fields in Vx configuration
! var_ind : pointer of variable field in SESAM Vx vector
! var_nbr : size of variable field in SESAM Vx vector
      INTEGER, save :: varend
      INTEGER, dimension(1:nbvar), save :: var_ord,var_dim, &
     &      var_jpi,var_jpj,var_jpk,var_jpt,var_nbr,var_ind, &
     &      varipos,varopos,varemsk,varpmsk,vardmsk, &
     &      varegrd,varngrd
      BIGREAL, dimension(1:nbvar), save :: var_moy,var_ect
      BIGREAL4, dimension(1:nbvar), save :: varvmsk
      BIGREAL, dimension(:,:), allocatable, save :: var_lev
      CHARACTER(len=varlg), dimension(1:nbvar), save :: var_nam
      CHARACTER(len=bgword), dimension(1:nbvar), save :: varfmsk, &
     &       varinam,varonam,varifil,varofil,varxdim,varydim, &
     &       varzdim,vartdim,varfgrd, &
     &       varxnam,varynam,varznam,vartnam
      LOGICAL, dimension(1:nbvar), save :: varmsea
! Vy object configuration
! dtaend : number of variable fields in Vx configuration
! dta_ind : pointer of variable field in SESAM Vy vector
! dta_nbr : size of variable field in SESAM Vy vector
      INTEGER, save :: dtaend
      INTEGER, dimension(1:nbvar), save :: dta_ord,dta_dim, &
     &      dta_jpi,dta_jpj,dta_jpk,dta_jpt,dta_nbr,dta_ind, &
     &      dtaipos,dtaopos,dtaemsk,dtapmsk,dtadmsk, &
     &      dtaegrd,dtangrd,dtaelev
      BIGREAL, dimension(1:nbvar), save :: dta_moy,dta_ect,dta_rms
      BIGREAL4, dimension(1:nbvar), save :: dtavmsk
      CHARACTER(len=varlg), dimension(1:nbvar), save :: dta_nam
      CHARACTER(len=bgword), dimension(1:nbvar), save :: dtainam, &
     &     dtaonam,dtaifil,dtaofil,dtafmsk,dtafgrd,dtaflev, &
     &     dtaxdim,dtaydim,dtazdim,dtatdim, &
     &     dtaxnam,dtaynam,dtaznam,dtatnam
      LOGICAL, dimension(1:nbvar), save :: dta_act,dtamsea
! Vo object configuration
! obsend : number of observations in Vx configuration
! obs_ind : pointer of observation property in SESAM Vo vector
! obs_nbr : size of observation property in SESAM Vo vector
      INTEGER, save :: obsend
      INTEGER, dimension(1:nbvar), save :: obsndbs
      INTEGER, dimension(1:nbobs), save :: obs_ord,obsnord
      INTEGER, dimension(1:nbvar,1:jpndbs), save :: obs_dim, &
     &     obs_nbr,obs_ind,obs_itp
      BIGREAL, dimension(1:nbvar,1:jpndbs), save :: obs_moy, &
     &     obs_ect,obs_rms
      CHARACTER(len=varlg), dimension(1:nbvar,1:jpndbs), save :: obs_nam
      CHARACTER(len=bgword), dimension(1:nbvar,1:jpndbs), save :: &
     &     obsinam,obsonam,obsifil,obsxdim,obsxnam,obsydim,obsynam, &
     &     obszdim,obsznam,obstdim,obstnam
      INTEGER, dimension(1:nbvar,1:jpndbs), save :: obsimin,obsimax, &
     &     obsjmin,obsjmax,obskmin,obskmax,obstmin,obstmax,obs_siz_max
      BIGREAL4, dimension(1:nbvar,1:jpndbs), save :: obsexcl
      BIGREAL4, dimension(1:nbvar,1:jpndbs), save :: obs_min,obs_max, &
     &     obs_lon_min,obs_lon_max,obs_lat_min,obs_lat_max, &
     &     obs_dep_min,obs_dep_max,obs_tim_min,obs_tim_max, &
     &     obs_lon_def,obs_lat_def,obs_dep_def,obs_tim_def
! Vz object configuration
! zon_jpi,zon_jpj,zon_jpk,zon_jpt : size of local data sections
! jpbub : number of different local data sections
      INTEGER, save :: zon_jpi,zon_jpj,zon_jpk,zon_jpt,jpbub
! Nesting configuration
! nestend : number of nested grids
! nestx1, nesty1 : coordinates of the nested grid south-west corner in grid 0
! nestxres, nestyres : nested grid resolution factor
! varnest : grid level for each Vx variable
      INTEGER, save :: nestend
      INTEGER, dimension(0:nbnest), save :: nestxori, nestyori
      INTEGER, dimension(0:nbnest), save :: nestxres, nestyres
      INTEGER, dimension(1:nbvar), save :: varnest
!
! -4.- SESAM namelist variables
! -----------------------------
! ==> SESAM global parameters
      LOGICAL, save :: lmoyect = .FALSE.
      LOGICAL, save :: traditional = .FALSE.
      INTEGER, save :: nprint = 0
! ==> Parameters for roa modules
      BIGREAL, save :: factoubli = 1.0
      BIGREAL, save :: forgexp = 1.0
! ==> Parameters for adap module
      INTEGER, save :: oecorreltyp = 1
! ==> Parameters for aanm module
      CHARACTER(len=80), save :: obserror_type_sesam = 'gaussian'
      BIGREAL, save :: special_value = -999.
! ==> Parameters for spct module
      CHARACTER(len=80) :: regr_type_sesam='local' ! regression type
      INTEGER, SAVE :: regr_maxiter_sesam=50   ! maximum number of iteration
      INTEGER, SAVE :: regr_maxbloc_sesam=1    ! maximum size of local blocks
      INTEGER, SAVE :: regr_overlap_sesam=1    ! overlapping of local blocks
      REAL(KIND=8), SAVE :: regr_epsilon_sesam=0.01 ! relative difference for convergence test
      REAL(KIND=8), SAVE :: regr_rho_sesam=1.0 ! weight according to signal std (rho=1) or not (rho=0)
! ==> Parameters for mcmc module
      INTEGER, SAVE :: jpscl=2                       ! number of scales in multiple scale ensemble
      INTEGER, dimension(1:nbscl), save :: scl_mult  ! multiplicity of each scale in Schur products
!
! -5.- Miscellaneous variables
! ----------------------------
! SESAM module and action index
      INTEGER, save :: nmode,naction
! Parameters for i/o operations
      CHARACTER(len=bgword), parameter :: clold='OLD'
      CHARACTER(len=bgword), parameter :: clnew='NEW'
      CHARACTER(len=bgword), parameter :: clunk='UNKNOWN'
      CHARACTER(len=bgword), parameter :: clfor='FORMATTED'
      CHARACTER(len=bgword), parameter :: clunf='UNFORMATTED'
      CHARACTER(len=bgword), parameter :: clseq='SEQUENTIAL'
      CHARACTER(len=bgword), parameter :: cldir='DIRECT'
! Variables for i/o operations
      CHARACTER(len=bgword), save :: clname
      INTEGER, save :: irecl, iost
      INTEGER, save :: wribefore=0
! Error message
      CHARACTER(len=4*bgword), save :: texterror
! Existence of various objects in current action
      LOGICAL, save :: existbas,existybas,existobas,existzbas, &
     &     existdbs,existdta,existobs,existvar,existzon,existcstr
!
      END MODULE mod_cfgxyo
!----------------------------------------------------------------------
