#include "readc.h"

#include "TTree.h"
#include "TFile.h"

#include <iostream> 
#include <memory> 
#include <vector> 

int main(int argc, char** argv){

    if (argc < 3){
        std::cerr << "Usage: cToRoot <target ROOT file> <one or many mille binaries>"<<std::endl;
        return 1; 
    }
    const std::string outFile = argv[1]; 
    auto fout = std::make_unique<TFile>(outFile.c_str(), "RECREATE"); 
    std::vector<std::string> inMilleFiles{};
    inMilleFiles.reserve(argc - 2); 
    for (size_t i = 2; i < argc; ++i) inMilleFiles.push_back(argv[i]); 
    initc(inMilleFiles.size()); 
    std::vector<double> doubleBuffer (1e6,0);
    std::vector<float> floatBuffer (1e6,0);
    std::vector<int> intBuffer (1e6,0);

    auto tree = std::make_unique<TTree>("MilleRecords","MilleRecords"); 
    tree->SetDirectory(fout.get()); 
    std::vector<double> doubles{};
    std::vector<float> floats{};
    std::vector<int> ints{}; 
    doubles.reserve(10000); 
    floats.reserve(10000); 
    ints.reserve(10000); 
    int floatingPointFormat = 0; 
    double residual; 
    double resError; 
    int trackIndex = 0; 
    tree->Branch("ints",&ints); 
    tree->Branch("floats",&floats); 
    tree->Branch("doubles",&doubles); 
    tree->Branch("floatingPointFormat",&floatingPointFormat); 
    tree->SetAutoFlush(100); 

    // read from C-files
    int err = 0; 
    for (size_t iFile = 0; iFile < inMilleFiles.size(); ++iFile){
        openc(inMilleFiles.at(iFile).c_str(),inMilleFiles.at(iFile).size(), iFile+1, &err);
        if (err){
            std::cerr <<" OPEN ERROR: "<<err<<std::endl;    
            continue; 
        }
        err = 4; 
        while(err > 0){
            floatingPointFormat = 0; 
            doubles.clear();
            floats.clear(); 
            ints.clear();
            int lBuffers = doubleBuffer.size(); 
            readc(doubleBuffer.data() , floatBuffer.data() , intBuffer.data() , &lBuffers, iFile+1, &err); 
            if (err < 0){
                std::cerr <<" READ ERROR: "<<err<<std::endl;    
                break; 
            }
            if (err == 0){
                // EOF - can stop 
                std::cout << " We are done! " <<std::endl;
                break; 
            }
            if (err == 4){
                floatingPointFormat = 4; 
                floats.assign(floatBuffer.begin(),floatBuffer.begin()+lBuffers); 
            }
            else if (err == 8){
                floatingPointFormat = 8; 
                doubles.assign(doubleBuffer.begin(),doubleBuffer.begin()+lBuffers); 
            }
            ints.assign(intBuffer.begin(),intBuffer.begin()+lBuffers); 
            tree->Fill(); 
            if (++trackIndex % 1000 == 0) std::cout << " done with record "<<trackIndex<<std::endl; 
        }
        std::cout << " ==> Done with "<<inMilleFiles[iFile]<<std::endl;
    }
    fout->Write();
    // dump to ROOT 
    return 0;
}
