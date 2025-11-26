
/** \file
 *  Read from ROOT files - implementation mimicking the readc interface.
 *
 * \author Max Goblirsch-Kolb, DESY (2025)
 *
 *  \copyright
 *  Copyright (c) 2025 Deutsches Elektronen-Synchroton,
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

#include "readc.h"
#include <stdio.h>
#include <vector>
#ifdef USE_ZLIB
#include <zlib.h>
#endif
#include <string>
#include <memory>
#include <iostream>
#include "TTree.h"
#include "TFile.h"

class fileHandler{
    public:
    fileHandler( const std::string & f=""){
        file_.reset(TFile::Open(f.c_str(),"READ")); 
		if (!file_ || !file_->IsOpen()){
            throw std::runtime_error("could not open file "+f);
		}
    };
	bool connect(){
		file_->GetObject<TTree>("MilleRecords",tree_);
		if (!tree_){
            throw std::runtime_error("could not read tree from file "+std::string(file_->GetName()));
			return false; 
		}
		else {
			tree_->SetBranchAddress("doubles",&doubles_); 
			tree_->SetBranchAddress("ints",	  &ints_); 
			tree_->SetCacheSize(100000000);	// 100 MB
			tree_->AddBranchToCache("*",true); 
			tree_->StopCacheLearningPhase();
			nentries_ = tree_->GetEntries(); 
			ints_ = nullptr;
			doubles_ = nullptr; 
		}	
		return true; 
	}
	int readNext(double* targetDouble, int* targetInt, int maxBufferSize){
		if (entryNumber_ >= nentries_){
			// EOF 
			return 0; 
		}
		int got = tree_->GetEntry(entryNumber_++);
		if (!got) return 0; 
		int entrySize = ints_->size();
		if (entrySize < maxBufferSize){
			std::memcpy(targetInt, 	  ints_->data(),   entrySize*sizeof(ints_->front()   )); 
			std::memcpy(targetDouble, doubles_->data(),entrySize*sizeof(doubles_->front())); 
		}
		// buffer too small - return required size 
		else{
			return -1 * entrySize; 
		}
		return entrySize; 
	}
	void rewind(){
		entryNumber_ = 0;
	}
	void close(){
		file_->Close(); 
		tree_ = nullptr;
	}
	~fileHandler(){
		std::cout <<__LINE__<<std::endl;
		file_.release(); 
		std::cout <<__LINE__<<std::endl;
		std::cout <<__LINE__<<std::endl;
		doubles_ = nullptr;
		std::cout <<__LINE__<<std::endl;
		ints_ = nullptr;
		std::cout <<__LINE__<<std::endl;
	}

    private:
	std::unique_ptr<TFile> file_ = nullptr;
	TTree* tree_ = nullptr;
	Long64_t entryNumber_ = 0; 
	Long64_t nentries_ = 0; 
	std::vector<double>* doubles_ = nullptr; 
	std::vector<int>* ints_ = nullptr; 
};

/* ________ global variables used for file handling __________ */

std::vector<std::unique_ptr<fileHandler>> files_;   ///< pointer to list of pointers to opened binary files

/*______________________________________________________________*/
/// Initialises the 'global' variables used for file handling.
/**
 * \param[in]  nFiles  Maximal number of C binary files to use.
 */
void initc(int nFiles) {
	files_.reserve(nFiles); 
}

/*______________________________________________________________*/
/// Open file.
void openc(const char *fileName, int lengthFileName, int nFileIn, int *errorFlag)
/**
 * \param[in]  fileName  File name
 * \param[in]  lengthFileName  Length of file name
 * \param[in]  nFileIn  File number (1 .. maxNumFiles) or <=0 for next one
 * \param[out] errorFlag error flag:
 *      * 0: if file opened and OK,
 *      * 1: if too many files open,
 *      * 2: if file could not be opened
 *      * 3: if file opened, but with error (can that happen?)
 */
{

	std::string fname = fileName;
	fname = fname.substr(0,lengthFileName); 
	fname[lengthFileName] = '\0';
	int fileIndex = nFileIn - 1; /* index of specific file */
	if (fileIndex < 0) fileIndex = files_.size(); /* next one */
	if (files_.size() <= fileIndex) files_.resize(fileIndex+1); 
	files_.at(fileIndex) = std::move(std::make_unique<fileHandler>(fname.c_str())); 
	files_.at(fileIndex)->connect();
}

/*______________________________________________________________*/
/// Close file.
/**
 * \param[in]  nFileIn  File number (1 .. maxNumFiles)
 */
void closec(int nFileIn) {
	files_.at(nFileIn-1)->close();
}

/*______________________________________________________________*/
/// Rewind file.
/**
 * \param[in]  nFileIn  File number (1 .. maxNumFiles)
 */
void resetc(int nFileIn) {
	files_.at(nFileIn-1)->rewind(); 
}

/*______________________________________________________________*/
/// Read record from file.
/**
 * \param[out]    bufferDouble  read buffer for doubles
 * \param[out]    bufferFloat   read buffer for floats
 * \param[out]    bufferInt     read buffer for integers
 * \param[in,out] lengthBuffers in: buffer length, out: number of floats/ints in records
 *                              (> buffer size: record skipped)
 * \param[in]     nFileIn       File number (1 .. maxNumFiles)
 * \param[out]    errorFlag     error flag:
 *      *  -1: pointer to a buffer or lengthBuffers are null
 *      *  -2: problem reading record length
 *      *  -4: given buffers too short for record
 *      *  -8: problem with stream or EOF reading floats
 *      * -16: problem with stream or EOF reading ints
 *      * -32: problem with stream or EOF reading doubles
 *      *  =0: reached end of file (or read empty record?!)
 *      *  =4: found floats
 *      *  =8: found doubles
 */
void readc(double *bufferDouble, float *bufferFloat, int *bufferInt,
		int *lengthBuffers, int nFileIn, int *errorFlag) {
	/* No return value since to be called as subroutine from fortran,
	 negative *errorFlag are errors, otherwise fine.

	 *nFileIn: number of the file the record is read from,
	 starting from 1 (not 0)
	 */
	if (nFileIn > files_.size()){
		std::cerr <<" Requested file " <<nFileIn<<" out of " <<files_.size()<<std::endl; 
		return; 
	}
	auto fH = files_.at(nFileIn-1).get(); 
	*errorFlag = 0;
	if (!bufferFloat || !bufferInt || !lengthBuffers) {
		*errorFlag = -1;
		return;
	}
	int nRead = fH->readNext(bufferDouble, bufferInt, *lengthBuffers); 
	if (nRead < 0){
		printf("readC: given buffers too short (%d, need > %d)\n", *lengthBuffers,
		 -nRead); 
		/* skip floats */
		*errorFlag = -4;
		*lengthBuffers = -nRead;
		return;
	} 
	else if (nRead == 0){
		// done! 
		return;
	}
	else {
		*lengthBuffers = nRead;
	}
	*errorFlag = 8;
}
