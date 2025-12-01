/**
 * \file Common pieces for I/O unit tests. 
    \author Maximilian Goblirsch-Kolb (DESY) 
    \copyright
    Copyright (c) 2025 Deutsches Elektronen-Synchroton,
    Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU Library General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. \n\n
    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Library General Public License for more details. \n\n
    You should have received a copy of the GNU Library General Public
    License along with this program (see the file COPYING.LIB for more
    details); if not, write to the Free Software Foundation, Inc.,
    675 Mass Ave, Cambridge, MA 02139, USA.
 */

#pragma once 

#include <iostream> 
#include <vector> 
#include <random> 
#include <functional> 
#include "Mille.h"

/// @brief Write a dummy binary with randomised content, and return information on the expected content 
/// @param [in] nTracks number of records (= tracks) to generate 
/// @param [in] fname name of the binary to write
/// @param [in] mode File type to write 
/// @param [out] expectedFloat populate with the expected content of the float array - outer vector represents records
/// @param [out] expectedInt populate with the expected content of the int array - outer vector represents records
/// @return size of the largest record - defines the needed buffer length. 
size_t writeBinary(size_t nTracks, const std::string & fname, Mille::OutputMode mode, std::vector<std::vector<float>> & expectedFloat, std::vector<std::vector<int>> & expectedInt){

    // seed the random engine with a constant to get reproducible behaviour 
    std::default_random_engine rdm(42); 
    std::uniform_int_distribution intDistro(1,255); // for labels
    std::uniform_int_distribution nMeasurements(1,6); // for multiplicity
    std::uniform_real_distribution floatDistro(-1.,1.); // for derivatives and residuals 

    Mille mille(fname.c_str(),mode); 
    std::vector<float> localDerivatives;
    std::vector<float> globalDerivatives;
    std::vector<int> localLabels;
    std::vector<int> globalLabels;
    size_t nMax = 0; 
    for (size_t iTrack = 0; iTrack < nTracks; ++iTrack){
        // track header: 0 - 0 (nErrors)
        expectedFloat.push_back(std::vector<float>{0.}); 
        expectedInt.push_back(std::vector<int>{0}); 
        auto & currentInts = expectedInt.back(); 
        auto & currentFloats = expectedFloat.back(); 
        
        // random set of measurements on track
        int nMeas = nMeasurements(rdm);

        // assume same n(locals) for each measurement on one track
        int nLoc = nMeasurements(rdm);

        for (size_t meas = 0; meas < nMeas; ++meas){
            // also randomise derivatives
            int nGlob = nMeasurements(rdm);
            // clear the local storage 
            localDerivatives.clear();
            globalDerivatives.clear();
            localLabels.clear();
            globalLabels.clear();

            // roll local derivatives
            for (int loc = 0; loc < nLoc; ++loc){
                // make ~half of the local derivatives zero 
                if (floatDistro(rdm) > 0){
                    localDerivatives.push_back(floatDistro(rdm));
                    localLabels.push_back(loc+1);
                }
                else{
                    localDerivatives.push_back(0); 
                }
            }
            // roll local derivatives
            for (int glob = 0; glob < nGlob; ++glob){
                globalDerivatives.push_back(floatDistro(rdm));
                globalLabels.push_back(intDistro(rdm));
            }
            // roll measurement outcome 
            float residual = floatDistro(rdm); 
            float residualError = std::abs(floatDistro(rdm)); 

            // push into expected structure 

            // residual 
            currentInts.push_back(0);
            currentFloats.push_back(residual);
            
            // local derivatives
            for (float d : localDerivatives){
                if (d!=0) currentFloats.push_back(d); 
            }
            currentInts.insert(     currentInts.end(),  localLabels.begin(),        localLabels.end()); 
            
            // measurement error 
            currentInts.push_back(0);
            currentFloats.push_back(residualError);
            
            // global derivatives
            currentInts.insert(     currentInts.end(),  globalLabels.begin(),       globalLabels.end()); 
            currentFloats.insert(   currentFloats.end(),globalDerivatives.begin(),  globalDerivatives.end()); 

            // now also store in Mille 
            mille.mille(nLoc,localDerivatives.data(),nGlob,globalDerivatives.data(),globalLabels.data(),residual,residualError); 

        }
        // check max record size
        if (currentFloats.size() >= nMax) nMax = currentFloats.size(); 
        // track done
        mille.end();
    }
    return nMax; 
}

/// @brief Helper function to test reading back binaries and compare them to a reference.
/// The I/O functions are passed as function pointers, allowing to test different
/// implementations of the same interface. This would be so much nicer 
/// using C++ polymorphism, but the methods have to be called from fortran,
/// so we prefer different symbol names. 
/// @param nTracks  Number of tracks used for the test (= expectation of records to find)
/// @param binaryName  Name of the binary to read
/// @param bufferSize Buffer size - can be set using the return of writeBinary 
/// @param expectedFloat - expected float buffer structure for each record. 
/// @param expectedInt - expected int buffer structure for each record 
/// @param initfunc - initialisation function to call (e.g. initc / initroot) 
/// @param openfunc - file open function to call (e.g. openc/openroot) 
/// @param readfunc - reading function to call (e.g. readc/readroot) 
/// @param rewindfunc - rewind function to call (e.g. resetc/resetroot) 
/// @return status code - 0 is success, non-zero integers denote failure modes.
int testReadBack(size_t nTracks, 
                 const std::string & binaryName, 
                 int bufferSize, 
                 const std::vector<std::vector<float>> & expectedFloat, 
                 const std::vector<std::vector<int>> & expectedInt,
                 std::function<void(int)> initfunc,
                 std::function<void(const char *, int , int , int *)> openfunc,
                 std::function<void(double *, float *, int *, int *, int , int *)> readfunc,
                 std::function<void(int)> rewindfunc)

                 {  
    // book read buffer
    std::vector<double> bufferDouble(bufferSize,0);
    std::vector<float> bufferFloat(bufferSize,0);
    std::vector<int> bufferInt(bufferSize,0); 
    // init the reader 
    initfunc(1);
    // open the file, check if successful 
    int errorFlag = 0; 
    openfunc(binaryName.c_str(), binaryName.size(), 1, &errorFlag); 
    if (errorFlag){
        std::cerr << "Open failed with status "<<errorFlag<<std::endl;
        return 1; 
    }
    const float tolerance = std::numeric_limits<float>::epsilon();
    int bufSize = bufferSize; 
    // read tracks sequentially 
    for (int k =0; k < nTracks; ++k){
        bufSize = bufferSize; 
        // call read method, check return 
        errorFlag = 0; 
        readfunc(bufferDouble.data(),bufferFloat.data(), bufferInt.data(), &bufSize, 1, &errorFlag); 
        if (errorFlag != 4){
            std::cerr << "Read did not deliver intended error flag - expected 4, got "<<errorFlag<<std::endl; 
            return 2; 
        }
        // retrieve the expected result 
        const std::vector<float> & expFloat = expectedFloat.at(k);
        const std::vector<int> & expInt = expectedInt.at(k);
        // first check size of returned arrays - compatible with expectation? 
        if (bufSize != expInt.size()){
            std::cerr <<" Read did not read the expected Int entries - expected "<<expInt.size()<<", got "<<bufSize<<std::endl; 
            return 3;
        }
        if (bufSize != expFloat.size()){
            std::cerr <<" Read did not read the expected Float entries - expected "<<expFloat.size()<<", got "<<bufSize<<std::endl; 
            return 3;
        }
        // now also validate the content. Can't wait for std::ranges::enumarate!!  
        bool contentOK = true; 
        for (size_t k = 0; k < expInt.size(); ++k){
            if (expInt.at(k) != bufferInt.at(k)){
                std::cerr << "Int buffer mismatch at position "<<k<<" - expect "<<expInt.at(k)<<", got "<<bufferInt.at(k)<<std::endl; 
                contentOK=false; 
            }
            if (std::abs(expFloat.at(k) - bufferDouble.at(k))>tolerance){
                std::cerr << "Int buffer mismatch at position "<<k<<" - expect "<<expFloat.at(k)<<", got "<<bufferDouble.at(k)<<std::endl; 
                contentOK=false; 
            }
        }
        if (!contentOK){
            std::cerr <<" Read file content not compatible with expectation. See above for mismatches. "<<std::endl;
            return 4; 
        }
    }
    // make sure reading beyond the number of written entries results in a null exit code. 
    errorFlag = 0; 
    bufSize = bufferSize; 
    readfunc(bufferDouble.data(),bufferFloat.data(), bufferInt.data(), &bufSize, 1, &errorFlag); 
    if (errorFlag != 0){
        std::cerr << " Read did not return error code 0 at EOF, instead got "<<errorFlag<<std::endl;
        return 5;  
    }

    // rewind the file
    rewindfunc(1);  
    errorFlag = 0; 
    bufSize = bufferSize; 
    readfunc(bufferDouble.data(),bufferFloat.data(), bufferInt.data(), &bufSize, 1, &errorFlag); 
    if (errorFlag != 4){
        std::cerr << " Read after rewind did not return expected code 4, instead got "<<errorFlag<<std::endl;
        return 6;  
    }

    // one more check if we are back to entry 1
    if (bufSize != expectedInt.at(0).size()){
        std::cerr <<" Read after rewind did not read the expected Int entries - expected "<<expectedInt.at(0).size()<<", got "<<bufSize<<std::endl; 
        return 7;
    }
    std::cout <<" All records read back as expected!" <<std::endl; 

    return 0; 

}
