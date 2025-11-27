
/** \file
*  Read from ROOT files - implementation of the readc interface.
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

#include "readROOT.h"
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

/// @brief Helper class to wrap the process of interacting with a ROOT TFile object
/// Internally contains a TFile to access one input file, and a TTree to read the 
/// contained information.
/// Provides an interface adapted to the 'readc' module, 
/// with functionality to open / close / rewind and to incrementally access
/// the next record.  
class fileHandler{
    public:
    /// @brief constructor
    /// @param f name of the ROOT file to load 
    /// Actual file access is delayed to the 
    /// connect() call. 
    fileHandler( const std::string & f){
        fname_ = f; 
    };
    /// @brief Open a ROOT file. 
    /// Load the contained Tree and configure it for reading. 
    /// Enable a TTreeCache to optimise network I/O performance.
    /// @return A bool, indicating the success of the operation. 
    bool connect(){
        // Access the ROOT file. 
        file_.reset(TFile::Open(fname_.c_str(),"READ")); 
        if (!file_ || !file_->IsOpen()){
            std::cerr << "could not open file "<<fname_<<std::endl;
            return false; 
        }
        // Load the tree
        file_->GetObject<TTree>("MilleRecords",tree_);
        if (!tree_){
            std::cerr << "could not read the tree 'MilleRecords' from file "<<fname_<<std::endl;
            return false; 
        }
        // connect our branches 
        tree_->SetBranchAddress("doubles",&doubles_); 
        tree_->SetBranchAddress("ints",	  &ints_); 
        tree_->SetBranchAddress("floats",	  &floats_); 
        // configure TTreeCache 
        tree_->SetCacheSize(100000000);	// 100 MB
        tree_->AddBranchToCache("*",true); 
        tree_->StopCacheLearningPhase();
        // check number of entries 
        nentries_ = tree_->GetEntries();
        // ensure valid status of internal pointers 
        ints_ = nullptr;
        doubles_ = nullptr; 
        floats_ = nullptr; 
        return true; 
    }
    /// @brief Read the next record from the connected tree. 
    /// Will populate the passed buffers. 
    /// Following the readc interface, the double buffer
    /// Will always be populated, even when reading floats. 
    /// @param targetDouble read buffer for doubles
    /// @param targetFloat read buffer for floats 
    /// @param targetInt read buffer for integers (= labels) 
    /// @param maxBufferSize buffer size 
    /// @param foundDoubles will be set to true if double-precision content was detected, else to false. 
    /// @return returns the number of read entries (recordLength / 2) in the given record. If the buffer size is exceeded,
    /// the negative size that would have been needed will be returned.  
    /// A zero return indicates EOF was reached. 
    int readNext(double* targetDouble, float* targetFloat, int* targetInt, int maxBufferSize, bool & foundDoubles){
        if (entryNumber_ >= nentries_){
            // We are done reading the tree - return 0 to indicate EOF 
            return 0; 
        }
        // read the next entry, increment internal entry counter
        int got = tree_->GetEntry(entryNumber_++);
        // if we read nothing, we are at the end of the file - return EOF
        if (!got) return 0; 
        // check the length of the label array for this record 
        int entrySize = ints_->size();
        // if we fit into the buffer, copy the information into the buffers  
        if (entrySize < maxBufferSize){
            // label info
            std::memcpy(targetInt, 	  ints_->data(),   entrySize*sizeof(ints_->front()   )); 
            // double-precision float info 
            if (doubles_->size() > 0){
                foundDoubles = true; 
                std::memcpy(targetDouble, doubles_->data(),entrySize*sizeof(doubles_->front())); 
            }
            // Single precision float info - here we also copy the information into the double buffer. 
            // Pede internally always accesses real numbers via the double buffer. 
            else{
                foundDoubles = false;
                std::memcpy(targetFloat, floats_->data(),entrySize*sizeof(floats_->front())); 
                size_t iD = 0; 
                // double-convert and return in double array 
                for (auto & f : *floats_) targetDouble[iD++] = f; 
            }
        }
        // buffer too small - return the required size 
        else{
            return -1 * entrySize; 
        }
        return entrySize; 
    }
    /// @brief Rewind the file to the first entry. 
    void rewind(){
        entryNumber_ = 0;
    }
    /// @brief Clos the file and reset our tree pointer
    void close(){
        file_->Close(); 
        tree_ = nullptr;
    }
    /// @brief Destructor. 
    ~fileHandler(){
        // Release the unique_ptr to prevent a ROOT related double-delete 
        file_.release(); 
        doubles_ = nullptr;
        ints_ = nullptr;
    }
    
    private:
    std::string fname_ = "";    ///< name of the ROOT input file 
    std::unique_ptr<TFile> file_ = nullptr; ///< TFile object to read input file 
    TTree* tree_ = nullptr; ///< TTree object to interface to the stored information
    Long64_t entryNumber_ = 0;  ///< Current entry number - this is the next one to be read
    Long64_t nentries_ = 0; ///< Total number of entries in the tree
    std::vector<double>* doubles_ = nullptr; ///< array to access double-precision reals in the tree
    std::vector<float>* floats_ = nullptr; ///< array to access single-precision reals in the tree
    std::vector<int>* ints_ = nullptr; ///< array to access parameter labels in the tree
};

/* ________ global variables used for file handling __________ */

std::vector<std::unique_ptr<fileHandler>> g_rootFileHandlers;   ///< pointer to list of pointers to opened binary files
unsigned int g_maxNumFiles;      ///< max number of files
unsigned int g_numAllFiles;      ///< number of opened files

/*______________________________________________________________*/
/// Initialises the 'global' variables used for file handling.
void initroot(int nFiles) {
    g_rootFileHandlers.resize(nFiles); 
    g_maxNumFiles = nFiles;
    g_numAllFiles = 0; 
}

/*______________________________________________________________*/
/// Open file.
void openroot(const char *fileName, int lengthFileName, int nFileIn, int *errorFlag)
{
    if (!errorFlag)
        return; /* 'printout' error? */
    int fileIndex = nFileIn - 1; /* index of specific file */

    if (fileIndex < 0) fileIndex = g_numAllFiles; /* next one */

    // invalid file index 
    if (fileIndex >= g_maxNumFiles) {
        *errorFlag = 1;
        return; 
    }
    // some sanitation of the input file name
    std::string fname = fileName;
    fname = fname.substr(0,lengthFileName); 
    fname[lengthFileName] = '\0';
    // now access the file 
    g_rootFileHandlers.at(fileIndex) = std::move(std::make_unique<fileHandler>(fname.c_str())); 
    if (!g_rootFileHandlers.at(fileIndex)->connect()){
        // failed to connect
        *errorFlag = 2; 
        return; 
    }
    // all good!
    ++g_numAllFiles; 
    *errorFlag = 0; 
}

/*______________________________________________________________*/
/// Close file.
void closeroot(int nFileIn) {
    int fileIndex = nFileIn - 1; /* index of current file */
    if (fileIndex < 0)
        return; /* no file opened at all... */
    g_rootFileHandlers.at(fileIndex)->close();
}

/*______________________________________________________________*/
/// Rewind file.
void resetroot(int nFileIn) {
    int fileIndex = nFileIn - 1; /* index of current file */
    if (fileIndex < 0)
        return; /* no file opened at all... */
    g_rootFileHandlers.at(fileIndex)->rewind();
}

/**______________________________________________________________*/
// Read record from file.

void readroot(double *bufferDouble, float *bufferFloat, int *bufferInt,
    int *lengthBuffers, int nFileIn, int *errorFlag) {
        /* No return value since to be called as subroutine from fortran,
        negative *errorFlag are errors, otherwise fine.
        
        *nFileIn: number of the file the record is read from,
        starting from 1 (not 0)
        */
	if (!errorFlag)
		return;
	*errorFlag = 0;
    // check for valid buffers. 
    if (!bufferFloat || !bufferInt || !lengthBuffers) {
        *errorFlag = -1;
        return;
    }
	int fileIndex = nFileIn - 1; /* index of current file */
	if (fileIndex < 0)
		return; /* no file opened at all... */
    // check for a valid file index
    if (fileIndex >= g_rootFileHandlers.size()){
        std::cerr <<" Requested file " <<nFileIn<<" out of " <<g_rootFileHandlers.size()<<std::endl;
		*errorFlag = -3; 
        return; 
    }
    /// check if the file handler is well-defined
    if (!g_rootFileHandlers.at(fileIndex)){
        std::cerr <<" File " <<nFileIn<<" out of " <<g_rootFileHandlers.size()<<" is in invalid state"<<std::endl;
		*errorFlag = -3; 
        return; 
    }
    // access the requested handler 
    auto fH = g_rootFileHandlers.at(fileIndex).get(); 
    // read the next record via ROOT's TTree I/O 
    bool foundDouble = false; 
    int nRead = fH->readNext(bufferDouble, bufferFloat, bufferInt, *lengthBuffers,foundDouble); 
    // negative nRead: buffer too small. |nRead| is the size we would need. 
    if (nRead < 0){
        printf("readC: given buffers too short (%d, need > %d)\n", *lengthBuffers,
            -nRead); 
            *errorFlag = -4;
            *lengthBuffers = -nRead;
            return;
    } 
    else if (nRead == 0){
        // done with the file. 
        *errorFlag = 0; 
        return;
    }
    else {
        // successful read. 
        *lengthBuffers = nRead;
    }
    // set the error flag to positive value to indicate a succesful read of either double or float precision
    *errorFlag = (foundDouble ? 8 : 4);
}
