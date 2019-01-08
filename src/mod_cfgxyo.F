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
CCC---------------------------------------------------------------------
CCC
CCC                        MODULE MOD_CFGXYO
CCC
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Module with SESAM configuration arrays
CCC  -------   (internal configuration defined in defcst.switch.h, and
CCC            user configuration defined in defcst.control.h)
CCC
CCC  original     : 99/05 (C.E. Testut)
CCC  modification : 03/02 (J.M. Brankart)
CCC---------------------------------------------------------------------
#include "config.main.h"
C
      MODULE mod_cfgxyo
C
      use mod_main
C
      IMPLICIT NONE
C
C -1.- Arrays with SESAM keywords (switches, modules, extensions,...)
C -------------------------------------------------------------------
C
      INTEGER, save :: nbdigits=4
      INTEGER, save :: nallmem,jpfixjpx,jpfixjpz,jpfixjpu
      INTEGER, dimension(1:nbxyo), save :: limjpnxyo
      INTEGER, dimension(1:nbztyp), save :: limjpnz,jpnbub,jpbubend
C
C Arrays with SESAM switches
C switab  : switch name
C swihelp : help about the nature of the argument of the switch
C swibool : switch is active (TRUE or FALSE)
C swiopt  : switch is optional (0=false, 1=general optional switch,
C                   2=optional switch particular to a few moodules).
      CHARACTER(len=swilg), dimension(1:nbarg), save :: switab
      CHARACTER(len=bgword), dimension(1:nbarg), save :: swihelp
      LOGICAL, dimension(1:nbarg), save :: swibool
      INTEGER, dimension(1:nbarg), save :: swiopt
C
C Arrays with SESAM modules
C modtab  : module name
C modhelp : help about the nature of the module
C modbool : module is active (TRUE or FALSE)
C modtype : type of module (1=analysis modules, 2=diagnostic modules,
C                  3=utility modules)
      CHARACTER(len=modlg), dimension(1:nbmod), save :: modtab
      CHARACTER(len=bgword), dimension(1:nbmod), save :: modhelp
      LOGICAL, dimension(1:nbmod), save :: modbool
      INTEGER, dimension(1:nbmod), save :: modtype
C
C Arrays with SESAM help keywords
C helptab  : help keyword name
C helpbool : help keyword is active (TRUE or FALSE)
      CHARACTER(len=helplg), dimension(1:nbhelp), save :: helptab
      LOGICAL, dimension(1:nbhelp), save :: helpbool
C
C Arrays with SESAM covariance directory extensions
C extbastab  : covariance directory extension name
C extbasbool : covariance directory extension is active (TRUE or FALSE)
C extbasmem  : type of format (unit=as output format, ten=as input format):
C                0=unavailable, 1=available, >1=available with restrictions
      CHARACTER(len=extlg), dimension(1:nbextbas), save :: extbastab
      LOGICAL, dimension(0:nbextbas), save :: extbasbool
      INTEGER, dimension(1:nbextbas), save :: extbasmem
C
C Arrays with SESAM observation database file extensions
C extdbstab  : observation database file extension name
C extdbsbool : observation database file extension is active (TRUE or FALSE)
C extdbsmem  : type of format (unit=as output format, ten=as input format):
C                0=unavailable, 1=available, >1=available with restrictions
      CHARACTER(len=extlg), dimension(1:nbextdbs), save :: extdbstab
      LOGICAL, dimension(0:nbextdbs), save :: extdbsbool
      INTEGER, dimension(1:nbextdbs), save :: extdbsmem
C
C Arrays with SESAM dta file extensions
C extdtatab  : dta file extension name
C extdtabool : dta file extension is active (TRUE or FALSE)
C extdtamem  : type of format (unit=as output format, ten=as input format):
C                0=unavailable, 1=available, >1=available with restrictions
C extdtaunit : different variables are in the same file (TRUE or FALSE)
      CHARACTER(len=extlg), dimension(1:nbextdta), save :: extdtatab
      LOGICAL, dimension(0:nbextdta), save :: extdtabool
      LOGICAL, dimension(1:nbextdta), save :: extdtaunit
      INTEGER, dimension(1:nbextdta), save :: extdtamem
C
C Arrays with SESAM obs file extensions
C extobstab  : obs file extension name
C extobsbool : obs file extension is active (TRUE or FALSE)
C extobsmem  : type of format (unit=as output format, ten=as input format):
C               0=unavailable, 1=available, >1=available with restrictions
C extobsunit : different variables are in the same file (TRUE or FALSE)
      CHARACTER(len=extlg), dimension(1:nbextobs), save :: extobstab
      LOGICAL, dimension(0:nbextobs), save :: extobsbool
      LOGICAL, dimension(1:nbextobs), save :: extobsunit
      INTEGER, dimension(1:nbextobs), save :: extobsmem
C
C Arrays with SESAM var file extensions
C extvartab  : var file extension name
C extvarbool : var file extension is active (TRUE or FALSE)
C extvarmem  : type of format (unit=as output format, ten=as input format):
C                0=unavailable, 1=available, >1=available with restrictions
C extvarunit : different variables are in the same file (TRUE or FALSE)
      CHARACTER(len=extlg), dimension(1:nbextvar), save :: extvartab
      LOGICAL, dimension(0:nbextvar), save :: extvarbool
      LOGICAL, dimension(1:nbextvar), save :: extvarunit
      INTEGER, dimension(1:nbextvar), save :: extvarmem
C
C Arrays with SESAM zon file extensions
C extzontab  : zon file extension name
C extzonbool : zon file extension is active (TRUE or FALSE)
C extzonmem  : type of format (unit=as output format, ten=as input format):
C                0=unavailable, 1=available, >1=available with restrictions
C extzonunit : different variables are in the same file (TRUE or FALSE)
      CHARACTER(len=extlg), dimension(1:nbextzon), save :: extzontab
      LOGICAL, dimension(0:nbextzon), save :: extzonbool
      LOGICAL, dimension(1:nbextzon), save :: extzonunit
      INTEGER, dimension(1:nbextzon), save :: extzonmem
C
C Arrays defining all possible actions for every SESAM modules
C tabvalorder  : Number of switches for action (x10),
C                  and type of action (unit): obsolete feature
C tabhelporder : Help for every actions of each module
C taborder     : List of switch indices defining the action
C tabnumswiopt : List of optional switch indices particular to each module
      INTEGER, dimension(1:nbmod,1:nbaction,1:nborder), save :: taborder
      INTEGER, dimension(1:nbmod,1:nbaction), save :: tabvalorder
      CHARACTER(len=hgword), dimension(1:nbmod,1:nbaction), save :: 
     $     tabhelporder
      INTEGER, dimension(1:nbmod,1:nbswiopt), save :: tabnumswiopt
C
C -2.- Specific variables for every SESAM arguments
C -------------------------------------------------
C arg<switchname>  : argument of switch <switchname>
C larg<switchname> : argument of the switch exists (TRUE or FALSE)
C --> mode switch
      CHARACTER(len=bgword), save :: argmod
      LOGICAL, save :: largmod
C --> optional switches:
C -help -list -varmsk -dtamsk -weight -oestd -bias -outinfo -fixjpx -fixjpu
C -reducevar -reducedta -scale -biasdbs -oecorrel -fecorrel -coefrmax -disable
C -inrz -outrz -outparadap -fixjpz -inpartobs -outdiaghst -inparadap -action
      CHARACTER(len=bgword), save :: arghelp,arglist,argvarmsk,
     $     argdtamsk,argweight,argoestd,
     $     argbias,argoutinfo,argfixjpx,argreducevar,argreducedta,
     $     argscale,argbiasdbs,argoecorrel,argfecorrel,
     $     argcoefrmax,argdisable,argoutparadap,arginrz,
     $     argoutrz,arginpartobs,argoutdiaghst,argfixjpz,
     $     arginparadap,argfixjpu,argaction
      LOGICAL, save :: larghelp,larglist,largvarmsk,largdtamsk,
     $     largweight,largoestd,
     $     largbias,largoutinfo,largfixjpx,largreducevar,largreducedta,
     $     largscale,largbiasdbs,largoecorrel,largfecorrel,
     $     largcoefrmax,largdisable,largoutparadap,larginrz,
     $     largoutrz,larginpartobs,largoutdiaghst,largfixjpz,
     $     larginparadap,largfixjpu,largaction
C --> required switches: in<...>
      CHARACTER(len=bgword), save :: arginxbas,argindbs,argindta,
     $     arginobs,arginvar,arginxbasref,argindbsref,argindtaref,
     $     arginobsref,arginvarref,arginpartvar,arginzon,
     $     arginzonref,
     $     arginybas,arginobas,arginzbas,arginptzon,
     $     arginybasref,arginobasref,arginzbasref, argincstr
      LOGICAL, save :: larginxbas,largindbs,largindta,larginobs,
     $     larginvar,larginxbasref,largindbsref,largindtaref,
     $     larginobsref,larginvarref,larginpartvar,larginzon,
     $     larginzonref,
     $     larginybas,larginobas,larginzbas,larginptzon,
     $     larginybasref,larginobasref,larginzbasref, largincstr
C --> required switches: out<...>
      CHARACTER(len=bgword), save :: argoutxbas,argoutdta,argoutobs,
     $     argoutvar,argoutxbasref,argoutdtaref,argoutobsref,
     $     argoutvarref,argoutpartvar,argoutzon,
     $     argoutybas,argoutobas,argoutzbas,argoutptzon
      LOGICAL, save :: largoutxbas,largoutdta,largoutobs,largoutvar,
     $     largoutxbasref,largoutdtaref,largoutobsref,largoutvarref,
     $     largoutpartvar,largoutzon,
     $     largoutybas,largoutobas,largoutzbas,largoutptzon
C --> other required switches:
      CHARACTER(len=bgword), save :: argtypeoper,argtypedtadiag,
     $     argdiffobsref,argdiffdtaref,argdiffvarref,
     $     argdiffobsorg,argdiffdtaorg,argdiffvarorg,
     $     arginerrdta,argouterrdta,arginstddta,argoutstddta,
     $     argoutobasref,argoutybasref,argoutbiasdta,
     $     argaffectobs,argnullobs,argconfigobs,argconfigzon,
     $     argzonindex,argincfg,arginoptcfg,argconnect,
     $     arginsmocfg,argoutsmocfg
      LOGICAL, save :: largtypeoper,largtypedtadiag,
     $     largdiffobsref,largdiffdtaref,largdiffvarref,
     $     largdiffobsorg,largdiffdtaorg,largdiffvarorg,
     $     larginerrdta,largouterrdta,larginstddta,largoutstddta,
     $     largoutobasref,largoutybasref,largoutbiasdta,
     $     largaffectobs,largnullobs,largconfigobs,largconfigzon,
     $     largzonindex,largincfg,larginoptcfg,largconnect,
     $     larginsmocfg,largoutsmocfg
C
C -3.- Arrays with SESAM user configuration
C -----------------------------------------
C See defcst.control.h for a description of the variables
C that are not described here
C
C Vx object configuration
C varend  : number of variable fields in Vx configuration
C var_ind : pointer of variable field in SESAM Vx vector
C var_nbr : size of variable field in SESAM Vx vector
      INTEGER, save :: varend
      INTEGER, dimension(1:nbvar), save :: var_ord,var_dim,
     $      var_jpi,var_jpj,var_jpk,var_jpt,var_nbr,var_ind,
     $      varipos,varopos,varemsk,varpmsk,vardmsk,
     $      varegrd,varngrd
      BIGREAL, dimension(1:nbvar), save :: var_moy,var_ect
      BIGREAL4, dimension(1:nbvar), save :: varvmsk
      BIGREAL, dimension(:,:), allocatable, save :: var_lev
      CHARACTER(len=varlg), dimension(1:nbvar), save :: var_nam
      CHARACTER(len=bgword), dimension(1:nbvar), save :: varfmsk,
     $       varinam,varonam,varifil,varofil,varxdim,varydim,
     $       varzdim,vartdim,varfgrd,
     $       varxnam,varynam,varznam,vartnam
      LOGICAL, dimension(1:nbvar), save :: varmsea
C Vy object configuration
C dtaend : number of variable fields in Vx configuration
C dta_ind : pointer of variable field in SESAM Vy vector
C dta_nbr : size of variable field in SESAM Vy vector
      INTEGER, save :: dtaend
      INTEGER, dimension(1:nbvar), save :: dta_ord,dta_dim,
     $      dta_jpi,dta_jpj,dta_jpk,dta_jpt,dta_nbr,dta_ind,
     $      dtaipos,dtaopos,dtaemsk,dtapmsk,dtadmsk,
     $      dtaegrd,dtangrd,dtaelev
      BIGREAL, dimension(1:nbvar), save :: dta_moy,dta_ect,dta_rms
      BIGREAL4, dimension(1:nbvar), save :: dtavmsk
      CHARACTER(len=varlg), dimension(1:nbvar), save :: dta_nam
      CHARACTER(len=bgword), dimension(1:nbvar), save :: dtainam,
     $     dtaonam,dtaifil,dtaofil,dtafmsk,dtafgrd,dtaflev,
     $     dtaxdim,dtaydim,dtazdim,dtatdim,
     $     dtaxnam,dtaynam,dtaznam,dtatnam
      LOGICAL, dimension(1:nbvar), save :: dta_act,dtamsea
C Vo object configuration
C obsend : number of observations in Vx configuration
C obs_ind : pointer of observation property in SESAM Vo vector
C obs_nbr : size of observation property in SESAM Vo vector
      INTEGER, save :: obsend
      INTEGER, dimension(1:nbvar), save :: obsndbs
      INTEGER, dimension(1:nbobs), save :: obs_ord,obsnord
      INTEGER, dimension(1:nbvar,1:jpndbs), save :: obs_dim,
     $     obs_nbr,obs_ind,obs_itp
      BIGREAL, dimension(1:nbvar,1:jpndbs), save :: obs_moy,
     $     obs_ect,obs_rms
      CHARACTER(len=varlg), dimension(1:nbvar,1:jpndbs), save :: obs_nam
      CHARACTER(len=bgword), dimension(1:nbvar,1:jpndbs), save ::
     $     obsinam,obsonam,obsifil,obsxdim,obsxnam,obsydim,obsynam,
     $     obszdim,obsznam,obstdim,obstnam
      INTEGER, dimension(1:nbvar,1:jpndbs), save :: obsimin,obsimax,
     $     obsjmin,obsjmax,obskmin,obskmax,obstmin,obstmax,obs_siz_max
      BIGREAL4, dimension(1:nbvar,1:jpndbs), save :: obsexcl
      BIGREAL4, dimension(1:nbvar,1:jpndbs), save :: obs_min,obs_max,
     $     obs_lon_min,obs_lon_max,obs_lat_min,obs_lat_max,
     $     obs_dep_min,obs_dep_max,obs_tim_min,obs_tim_max,
     $     obs_lon_def,obs_lat_def,obs_dep_def,obs_tim_def
C Vz object configuration
C zon_jpi,zon_jpj,zon_jpk,zon_jpt : size of local data sections
C jpbub : number of different local data sections
      INTEGER, save :: zon_jpi,zon_jpj,zon_jpk,zon_jpt,jpbub
C Nesting configuration
C nestend : number of nested grids
C nestx1, nesty1 : coordinates of the nested grid south-west corner in grid 0
C nestxres, nestyres : nested grid resolution factor
C varnest : grid level for each Vx variable
      INTEGER, save :: nestend
      INTEGER, dimension(0:nbnest), save :: nestxori, nestyori
      INTEGER, dimension(0:nbnest), save :: nestxres, nestyres
      INTEGER, dimension(1:nbvar), save :: varnest
C
C -4.- SESAM namelist variables
C -----------------------------
C ==> SESAM global parameters
      LOGICAL, save :: lmoyect = .FALSE.
      LOGICAL, save :: traditional = .FALSE.
      INTEGER, save :: nprint = 0
C ==> Parameters for roa modules
      BIGREAL, save :: factoubli = 1.0
      BIGREAL, save :: forgexp = 1.0
C ==> Parameters for adap module
      INTEGER, save :: oecorreltyp = 1
C ==> Parameters for aanm module
      CHARACTER(len=80), save :: obserror_cdf = 'gaussian'
C ==> Parameters for spct module
      CHARACTER(len=80) :: regr_type='local' ! regression type
      INTEGER, SAVE :: regr_maxiter=50   ! maximum number of iteration
      INTEGER, SAVE :: regr_maxbloc=1    ! maximum size of local blocks
      INTEGER, SAVE :: regr_overlap=1    ! overlapping of local blocks
      REAL(KIND=8), SAVE :: regr_epsilon=0.01 ! relative difference for convergence test
      REAL(KIND=8), SAVE :: regr_rho=1.0 ! weight according to signal std (rho=1) or not (rho=0)
C
C -5.- Miscellaneous variables
C ----------------------------
C SESAM module and action index
      INTEGER, save :: nmode,naction
C Parameters for i/o operations
      CHARACTER(len=bgword), parameter :: clold='OLD'
      CHARACTER(len=bgword), parameter :: clnew='NEW'
      CHARACTER(len=bgword), parameter :: clunk='UNKNOWN'
      CHARACTER(len=bgword), parameter :: clfor='FORMATTED'
      CHARACTER(len=bgword), parameter :: clunf='UNFORMATTED'
      CHARACTER(len=bgword), parameter :: clseq='SEQUENTIAL'
      CHARACTER(len=bgword), parameter :: cldir='DIRECT'
C Variables for i/o operations
      CHARACTER(len=bgword), save :: clname
      INTEGER, save :: irecl, iost
      INTEGER, save :: wribefore=0
C Error message
      CHARACTER(len=4*bgword), save :: texterror
C Existence of various objects in current action
      LOGICAL, save :: existbas,existybas,existobas,existzbas,
     $     existdbs,existdta,existobs,existvar,existzon,existcstr
C
      END MODULE mod_cfgxyo
CC----------------------------------------------------------------------
