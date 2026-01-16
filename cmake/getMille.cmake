# Get an external Mille installation 

include(FetchContent)

set(MILLE_TARGET_VERSION main) 
# Flag to adapt config file creation later on 
set(MILLE_AS_SUBDIR True) 


FetchContent_Declare(Mille 
    GIT_REPOSITORY https://gitlab.desy.de/millepede/mille.git
    GIT_TAG ${MILLE_TARGET_VERSION}
)

FetchContent_GetProperties(Mille)
if(NOT Mille_POPULATED)
  FetchContent_Populate(Mille)
  add_subdirectory(${mille_SOURCE_DIR} ${mille_BINARY_DIR})
  set(Mille_MODULE_DIRS ${mille_BINARY_DIR}/modules)
endif()
