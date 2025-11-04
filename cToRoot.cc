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
    std::cout << "hello, world "<<std::endl;
    initc(inMilleFiles.size()); 
    std::vector<double> doubleBuffer (1e7,0);
    std::vector<float> floatBuffer (1e7,0);
    std::vector<int> intBuffer (1e7,0);

    auto tree = std::make_unique<TTree>("MilleRecords","MilleRecords"); 
    std::vector<double> globalDeriv;
    std::vector<double> localDeriv;
    std::vector<double> specialFloats;
    std::vector<int> globalLabels; 
    std::vector<int> localLabels; 
    std::vector<int> specialInts; 
    double residual; 
    double resError; 
    int trackIndex = 0; 
    tree->Branch("globalDerivatives",&globalDeriv); 
    tree->Branch("localDerivatives",&localDeriv); 
    tree->Branch("specialFloats",&specialFloats); 
    tree->Branch("specialInts",&specialInts); 
    tree->Branch("globalLabels",&globalLabels); 
    tree->Branch("localLabels",&localLabels); 
    tree->Branch("residual",&residual); 
    tree->Branch("resError",&resError); 
    tree->Branch("trackIndex",&trackIndex); 

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
            int lBuffers = doubleBuffer.size(); 
            readc(doubleBuffer.data() , floatBuffer.data() , intBuffer.data() , &lBuffers, iFile+1, &err); 
            if (err < 0){
                std::cerr <<" READ ERROR: "<<err<<std::endl;    
                break; 
            }
            residual = floatBuffer[1]; 
            globalDeriv.clear(); 
            localDeriv.clear(); 
            globalLabels.clear(); 
            localLabels.clear(); 
            bool foundRes = false; 
            bool foundSigma = false; 
            bool foundSpecial = false; 
            for (int iIndex = 1; iIndex < lBuffers && trackIndex < 3; ++iIndex){
                std::cout << iIndex <<"  "<<intBuffer[iIndex]<<"    "<<floatBuffer[iIndex]<<std::endl;
                if (intBuffer[iIndex] == 0){
                    if (!foundRes){
                        residual = floatBuffer[iIndex]; 
                        foundRes = true; 
                        std::cout <<" << first entry >>"<<std::endl; 
                    }
                    else if (floatBuffer[iIndex] == 0 && floatBuffer[iIndex+1] < 0 && std::abs(float(int(floatBuffer[iIndex+1])) - floatBuffer[iIndex+1]) < 1e-4  && !foundSpecial){
                        std::cout <<" << special block >>"<<std::endl; 
                        foundSpecial = true; 
                        foundSigma = true; 
                    }
                    else if (!foundSigma){
                        std::cout <<" << sigma block >>"<<std::endl; 
                        resError = floatBuffer[iIndex]; 
                        foundSigma = true; 
                    }
                    else {
                        std::cout <<" << next entry >>"<<std::endl; 
                        foundRes = true;
                        foundSigma = false; 
                        foundSpecial = false; 
                        tree->Fill(); 
                        globalDeriv.clear(); 
                        localDeriv.clear(); 
                        globalLabels.clear(); 
                        localLabels.clear(); 
                        specialFloats.clear();
                        specialInts.clear();
                        residual = floatBuffer[iIndex]; 
                    }
                }
                else{
                    if (foundSpecial){
                        specialInts.push_back(intBuffer[iIndex]); 
                        specialFloats.push_back(floatBuffer[iIndex]); 
                    }
                    else if (foundSigma){
                        globalLabels.push_back(intBuffer[iIndex]); 
                        globalDeriv.push_back(floatBuffer[iIndex]); 
                    }
                    else{
                        localLabels.push_back(intBuffer[iIndex]); 
                        localDeriv.push_back(floatBuffer[iIndex]); 
                    }
                }
                std::cout << " [ res "<<residual <<" sig "<<resError<<" ix "<<trackIndex<<" #li "<< localLabels.size()<<" #gi "<<globalLabels.size()<<" #ld "<< localDeriv.size()<<" #gd "<<globalDeriv.size() <<" #sp "<<specialInts.size()<<std::endl; 
            }
            ++trackIndex;
            if (trackIndex % 100000 == 0) std::cout << " done with track "<<trackIndex<<std::endl; 
        }
        std::cout << " ==> Done with "<<inMilleFiles[iFile]<<std::endl;
    }
    fout->Write();
    // dump to ROOT 
    return 0;
}
