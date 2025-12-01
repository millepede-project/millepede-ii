
/**
    \file Unit test for the readc module.
    Define a known (non-physical) set of inputs,
    dump them using Mille as C binaries, and read them back
    using readc. 
    \author Maximilian Goblirsch-Kolb (DESY) 
    \author Claus Kleinwort, DESY (maintenance and developement)
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

#include "readc.h"
#include "testMilleIObase.h"


int main(int, char**){

    // generate a set of random inputs
    const size_t nTracks = 5; 
    const std::string binaryName = "TestCbinary.dat";

    // store the expected output 
    std::vector<std::vector<float>> expectedFloat{};
    std::vector<std::vector<int>> expectedInt{};
    
    // write the dummy binary 
    int readBufferSize = writeBinary(nTracks, binaryName, Mille::OutputMode::Cbinary, expectedFloat, expectedInt); 

    int ret = testReadBack(nTracks, binaryName, readBufferSize, expectedFloat, expectedInt,
                            initc, openc, readc, resetc); 

    if(!ret) std::cout << " Readback test on C binary OK!"<<std::endl;
    
    return ret; 
}
