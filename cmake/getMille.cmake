# Get an external Mille installation 

include(FetchContent)

# Flag to adapt config file creation later on 
set(MILLE_AS_SUBDIR True) 

set(Recommended_Mille_Version V01-00-02)

FetchContent_Declare(Mille 
    URL https://gitlab.desy.de/millepede/mille/-/archive/${Recommended_Mille_Version}/mille-${Recommended_Mille_Version}.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP True
    SOURCE_SUBDIR DummyForCMake
)

FetchContent_MakeAvailable(Mille)
add_subdirectory(${mille_SOURCE_DIR} ${mille_BINARY_DIR})
set(Mille_MODULE_DIRS ${mille_BINARY_DIR}/modules)
