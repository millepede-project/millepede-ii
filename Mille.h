#ifndef MILLE_H
#define MILLE_H

/** \file
 *  Define class Mille.
 *
 *  \author Gero Flucke, University Hamburg, 2006
 *
 *  \copyright
 *  Copyright (c) 2009 - 2015 Deutsches Elektronen-Synchroton,
 *  Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
 *  This library is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version. \n\n
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details. \n\n
 *  You should have received a copy of the GNU Library General Public
 *  License along with this program (see the file COPYING.LIB for more
 *  details); if not, write to the Free Software Foundation, Inc.,
 *  675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <fstream>
#include <vector>

#ifdef SUPPORT_ROOTIO 
// fwd declare ROOT classes if needed
class TTree;
class TFile; 
#endif 


/**
 * \class Mille
 *
 *  Class to write a C binary (cf. below) file of a given name and to fill it
 *  with information used as input to **pede**.
 *  Use its member functions \c mille(), \c special(), \c kill() and \c end()
 *  as you would use the fortran \ref mille.f90 "MILLE"
 *  and its entry points \c MILLSP, \c KILLE and \c ENDLE.
 *
 *  For debugging purposes constructor flags enable switching to text output and/or
 *  to write also derivatives and labels which are ==0.
 *  But note that **pede** will not be able to read text output and has not been tested with
 *  derivatives/labels ==0.
 *
 *  author    : Gero Flucke
 *  date      : October 2006
 *  $Revision: 1.3 $
 *  $Date: 2007/04/16 17:47:38 $
 *  (last update by $Author: flucke $)
 */

/// Class to write C binary file.
class Mille 
{
 public:

  enum class OutputMode{
    Cbinary=1,
    TextFile=2,
    RootBinary=3,
    _nOutModes
  };

  Mille(const char *outFileName, OutputMode out= OutputMode::Cbinary, bool writeZero = false);
  ~Mille();

  void mille(int NLC, const float *derLc, int NGL, const float *derGl,
	     const int *label, float rMeas, float sigma);
  void special(int nSpecial, const float *floatings, const int *integers);
  void kill();
  void end();

#ifdef SUPPORT_ROOTIO 
  void toROOT(); 
#endif 
 private:
  void newSet();
  bool checkBufferSize(int nLocal, int nGlobal);

  void openC(const std::string  & fname); 
  void openRoot(const std::string & fname); 

  void writeC(); 
  void writeText(); 
  void writeRoot(); 

  void closeC(); 
  void closeRoot(); 

  std::ofstream myOutFile; ///< C-binary for output
  OutputMode outputMode_ = OutputMode::Cbinary;         ///< output mode
  bool myWriteZero;        ///< if true also write out derivatives/labels ==0
  /// buffer size for ints and floats
  std::vector<int>   myBufferInt;   ///< to collect labels etc.
  std::vector<float> myBufferFloat; ///< to collect derivatives etc.
  int   myBufferPos; ///< position in buffer
  bool  myHasSpecial; ///< if true, special(..) already called for this record
  /// largest label allowed: 2^31 - 1
  enum {myMaxLabel = (0xFFFFFFFF - (1 << 31))};
#ifdef SUPPORT_ROOTIO 
  TFile* rootOutFile_ = nullptr; 
  TTree* outTree_ = nullptr; 
  std::vector<double> myBufferDummyDouble; ///< dummy
#endif 
};
#endif
