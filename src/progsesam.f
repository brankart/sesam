C Copyright: CNRS - Université de Grenoble
C
C Contributors : Jean-Michel Brankart, Charles-Emmanuel Testut, Laurent Parent,
C                Emmanuel Cosme, Claire Lauvernet, Frédéric Castruccio
C
C Jean-Michel.Brankart@hmg.inpg.fr
C
C -----------------------------------------------------------------------------
C
C This software is a computer program whose purpose is to perform
C the various basic operations that are required
C in sequential/stochastic data assimilation systems.
C
C This software is governed by the CeCILL license under French law and
C abiding by the rules of distribution of free software.  You can  use,
C modify and/ or redistribute the software under the terms of the CeCILL
C license as circulated by CEA, CNRS and INRIA at the following URL
C "http://www.cecill.info".
C
C As a counterpart to the access to the source code and  rights to copy,
C modify and redistribute granted by the license, users are provided only
C with a limited warranty  and the software's author,  the holder of the
C economic rights,  and the successive licensors  have only  limited
C liability.
C
C In this respect, the user's attention is drawn to the risks associated
C with loading,  using,  modifying and/or developing or reproducing the
C software by the user in light of its specific status of free software,
C that may mean  that it is complicated to manipulate,  and  that  also
C therefore means  that it is reserved for developers  and  experienced
C professionals having in-depth computer knowledge. Users are therefore
C encouraged to load and test the software's suitability as regards their
C requirements in conditions enabling the security of their systems and/or
C data to be ensured and,  more generally, to use and operate it in the
C same conditions as regards security.
C
C The fact that you are presently reading this means that you have had
C knowledge of the CeCILL license and that you accept its terms.
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C ---                                                           ---
C ---                    PROGSESAM.F                            ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12 (C.E. Testut)                        ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C --- 
C --- PROGRAM progsesam
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PROGRAM progsesam
CCC----------------------------------------------------------------
CCC
CCC   PURPOSE :
CCC   ---------
CCC       Main program calling the SESAM routine
CCC       with the action to perform in argument.
CCC       If the argument is empty, the SESAM action
CCC       will be read by SESAM as the commandline arguments.
CCC
CCC----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(len=1200) :: textcommand
C
      textcommand=''
C
      CALL sesam (textcommand)
C
      END    
