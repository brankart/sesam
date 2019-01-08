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
C---------------------------------------------------------------------
C
C                        CONFIG.MAIN
C
C---------------------------------------------------------------------
C
C  Purpose :
C  ---------
C     Definition of keywords specifying the type
C     of SESAM 'real' and 'integer' variables.
C     This file must be included in all SESAM FORTRAN files.
C   
C  Modifications :
C  -------------
C     original     : 97-12 (C.E. Testut)
C     modification : 99-12 (C.E. Testut)
C     modification : 03-02 (J.M. Brankart)
C
C---------------------------------------------------------------------
C List of keywords defined in this include file
C
C BIGREAL  : kind of real for SESAM internal computations
C BIGREAL8 : kind of real for SESAM i/o (8-byte real)
C BIGREAL4 : kind of real for SESAM i/o (4-byte real)
C FREAL    : function transforming integer or real variable into BIGREAL
C FREAL8   : function transforming integer or real variable into BIGREAL8
C FREAL4   : function transforming integer or real variable into BIGREAL4
C
C The definitions can be modified by precompilation options
C like 'cpp_real4' or 'cpp_real16', or as a function of the computer.
C---------------------------------------------------------------------
# undef cpp_real16
# undef cpp_real4
C---------------------------------------------------------------------
C List of keywords defining the type of computer:
C 
C cray     :  all CRAY computers
C _CRAYT3E :  CRAY T3E
C _CRAYC90 :  CRAY C90
C __uxpv__ :  Fujitsu
C SX       :  NEC SX-5
C hpux     :  HP K200
C
C---------------------------------------------------------------------
C Definition of keywords specifying the kind of real variables
C (Values for kr, kr4, kr8 are defined in 'mod_main.F' as
C a function of KR<y>EQ<nn> keywords defined in the following).
#define BIGREAL8 REAL(KIND=kr8)
#define BIGREAL4 REAL(KIND=kr4)
#define BIGREAL REAL(KIND=kr)
C Defintions for NEC SX-5 computers
#if defined SX
#  define _NEC
#  define FREAL8 dble
#  define FREAL4 real
#  define KR8EQ8
#  define KR4EQ4
#  define SMALLINTGR2 INTEGER(KIND=4)
#  if defined cpp_real4
#     define FREAL real
#     define KREQ4
#  elif defined cpp_real16
#     define FREAL qfloat
#     define KREQ16
#  else
#     define FREAL dble
#     define KREQ8
#  endif
C Defintions for other computers
C They work for Linux, IBM SP3 and SP4, Fujistsu, HP workstations.
#else
#  define FREAL8 dble
#  define FREAL4 real
#  define KR8EQ8
#  define KR4EQ4
#  define SMALLINTGR2 INTEGER*4
#  if defined cpp_real4
#     define FREAL real
#     define KREQ4
#  elif defined cpp_real16
#     define FREAL qfloat
#     define KREQ16
#  else
#     define FREAL dble
#     define KREQ8
#  endif
#endif
