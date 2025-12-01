
/** \file
 *  Create Millepede-II C-binary record.
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

/**
 *  author    : Gero Flucke, University Hamburg, 2006
 *  date      : October 2006
 *  $Revision: 1.3 $
 *  $Date: 2007/04/16 17:47:38 $
 *  (last update by $Author: flucke $)
 */

#include "Mille.h"

#include <fstream>
#include <iostream>

#ifdef SUPPORT_ROOTIO 
#include "TTree.h"
#include "TFile.h" 
#endif 

//___________________________________________________________________________

/// Opens outFileName (by default as binary file).
/**
 * \param[in] outFileName  file name
 * \param[in] asBinary     flag for binary
 * \param[in] writeZero    flag for keeping of zeros
 */
Mille::Mille(const char *outFileName, OutputMode out, bool writeZero) : 
  outputMode_(out), myWriteZero(writeZero), myBufferPos(-1), myHasSpecial(false)
{
  if (outputMode_ != OutputMode::RootBinary) openC(outFileName);
  else openRoot(outFileName); 
}

//___________________________________________________________________________
/// Closes file.
Mille::~Mille()
{
  if (outputMode_ != OutputMode::RootBinary){
    closeC();
  }
  else closeRoot(); 
}

//___________________________________________________________________________
/// Add measurement to buffer.
/**
 * \param[in]    NLC    number of local derivatives
 * \param[in]    derLc  local derivatives
 * \param[in]    NGL    number of global derivatives
 * \param[in]    derGl  global derivatives
 * \param[in]    label  global labels
 * \param[in]    rMeas  measurement (residuum)
 * \param[in]    sigma  error
 */
void Mille::mille(int NLC, const float *derLc,
		  int NGL, const float *derGl, const int *label,
		  float rMeas, float sigma)
{
  if (sigma <= 0.) return;
  if (myBufferInt.empty()) this->newSet(); // start, e.g. new track

  // reserve sufficient memory to avoid frequent reallocation
  myBufferFloat.reserve(NGL+NLC+2 + myBufferFloat.size()); 
  myBufferInt.reserve(NGL+NLC+2 + myBufferInt.size()); 

  // first store measurement
  myBufferFloat.push_back(rMeas);
  myBufferInt.push_back(0);

  // store local derivatives and local 'lables' 1,...,NLC
  for (int i = 0; i < NLC; ++i) {
    if (derLc[i] || myWriteZero) { // by default store only non-zero derivatives
      myBufferFloat.push_back(derLc[i]);
      myBufferInt.push_back(i+1);
    }
  }

  // store uncertainty of measurement in between locals and globals
  myBufferFloat.push_back(sigma);
  myBufferInt.push_back(0);

  // store global derivatives and their labels
  for (int i = 0; i < NGL; ++i) {
    if (derGl[i] || myWriteZero) { // by default store only non-zero derivatives
      if ((label[i] > 0 || myWriteZero) && label[i] <= myMaxLabel) { // and for valid labels
        myBufferFloat.push_back(derGl[i]); // global derivatives
        myBufferInt.push_back(label[i]);  // index of global parameter
      } else {
	std::cerr << "Mille::mille: Invalid label " << label[i] 
		  << " <= 0 or > " << myMaxLabel << std::endl; 
      }
    }
  }
}

//___________________________________________________________________________
/// Add special data to buffer.
/**
 * \param[in]    nSpecial   number of floats/ints
 * \param[in]    floatings  floats
 * \param[in]    integers   ints
 */
void Mille::special(int nSpecial, const float *floatings, const int *integers)
{
  if (nSpecial == 0) return;
  if (myBufferInt.empty()) this->newSet(); // start, e.g. new track
  if (myHasSpecial) {
    std::cerr << "Mille::special: Special values already stored for this record."
	      << std::endl; 
    return;
  }
  myHasSpecial = true; // after newSet() (Note: MILLSP sets to buffer position...)

  //  myBufferFloat[.]  | myBufferInt[.]
  // ------------------------------------
  //      0.0           |      0
  //  -float(nSpecial)  |      0
  //  The above indicates special data, following are nSpecial floating and nSpecial integer data.

  // zeros 
  myBufferFloat.push_back( 0.);
  myBufferInt  .push_back( 0);

  // nSpecial and zero
  myBufferFloat.push_back( nSpecial);
  myBufferInt  .push_back( 0);

  for (int i = 0; i < nSpecial; ++i) {
    myBufferFloat.push_back( floatings[i]);
    myBufferInt  .push_back( integers[i]);
  }
}

//___________________________________________________________________________
/// Reset buffers, i.e. kill derivatives accumulated for current set.
void Mille::kill()
{
  myBufferFloat.clear();
  myBufferInt.clear();
}

//___________________________________________________________________________
/// Write buffer (set of derivatives with same local parameters) to file.
void Mille::end()
{
  if (myBufferInt.size() > 1) { // only if anything stored...
    if (outputMode_ == OutputMode::Cbinary) {
      writeC(); 
    } else  if (outputMode_ == OutputMode::TextFile){
      writeText(); 
    }
    else{
      writeRoot(); 
    }
  }
  kill(); // reset for next entry 
} 

//___________________________________________________________________________
/// Initialize for new set of locals, e.g. new track.
void Mille::newSet()
{
  myHasSpecial = false;
  myBufferFloat.push_back(0);
  myBufferInt.push_back(0); 
}

//___________________________________________________________________________
/// Enough space for next nLocal + nGlobal derivatives incl. measurement?
/**
 * \param[in]   nLocal  number of local derivatives
 * \param[in]   nGlobal number of global derivatives
 * \return      true if sufficient space available (else false)
 */
bool Mille::checkBufferSize(int nLocal, int nGlobal)
{
  // if (myBufferPos + nLocal + nGlobal + 2 >= myBufferSize) {
  //   ++(myBufferInt[0]); // increase error count
  //   std::cerr << "Mille::checkBufferSize: Buffer too short (" 
	//       << myBufferSize << "),"
	//       << "\n need space for nLocal (" << nLocal<< ")"
	//       << "/nGlobal (" << nGlobal << ") local/global derivatives, " 
	//       << myBufferPos + 1 << " already stored!"
	//       << std::endl;
  //   return false;
  // } else {
    return true;
  // }
}

void Mille::openC(const std::string  & fname){
  myOutFile = std::ofstream(fname, (outputMode_ == OutputMode::Cbinary ? (std::ios::binary | std::ios::out) : std::ios::out));
  // Instead myBufferPos(-1), myHasSpecial(false) and the following two lines
  // we could call newSet() and kill()...
  if (!myOutFile.is_open()) {
    std::cerr << "Mille::openC: Could not open " << fname 
	      << " as output file." << std::endl;
  }
} 
void Mille::openRoot(const std::string & fname){
  #ifdef SUPPORT_ROOTIO
  rootOutFile_ = TFile::Open(fname.c_str(),"RECREATE"); 
  if (!rootOutFile_ || !rootOutFile_->IsOpen()){
    std::cerr << "Mille::openRoot: Could not open " << fname 
	      << " as output file." << std::endl;
  }
  outTree_ = new TTree("MilleRecords","MilleRecords");
  outTree_->SetDirectory(rootOutFile_); 
  outTree_->Branch("floats",&myBufferFloat); 
  outTree_->Branch("ints",&myBufferInt); 
  outTree_->Branch("doubles",&myBufferDummyDouble); 
  #else 
  std::cerr << " Error: Requesting ROOT output, but Mille was not built with ROOT support. Will not do anything. "<<std::endl; 
  #endif  
  
} 

void Mille::writeC(){
    const int numWordsToWrite = 2 * myBufferFloat.size(); 
    myOutFile.write(reinterpret_cast<const char*>(&numWordsToWrite), 
        sizeof(numWordsToWrite));
    myOutFile.write(reinterpret_cast<char*>(myBufferFloat.data()), 
        (myBufferFloat.size()) * sizeof(myBufferFloat[0]));
    myOutFile.write(reinterpret_cast<char*>(myBufferInt.data()), 
        (myBufferFloat.size()) * sizeof(myBufferInt[0]));
  
} 
void Mille::writeText(){
    const int numWordsToWrite = 2 * myBufferFloat.size(); 
    myOutFile << numWordsToWrite << "\n";
    for (int i = 0; i < myBufferPos+1; ++i) {
      myOutFile << myBufferFloat[i] << " ";
    }
    myOutFile << "\n";
    
    for (int i = 0; i < myBufferPos+1; ++i) {
      myOutFile << myBufferInt[i] << " ";
    }
    myOutFile << "\n";
  
} 
void Mille::writeRoot(){
#ifdef SUPPORT_ROOTIO
  outTree_->Fill();
#endif 
} 

void Mille::closeC(){
  myOutFile.close();
} 
void Mille::closeRoot(){
  #ifdef SUPPORT_ROOTIO
  rootOutFile_->Write(); 
  rootOutFile_->Close();
  #endif  
} 
