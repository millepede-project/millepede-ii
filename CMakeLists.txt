CMAKE_MINIMUM_REQUIRED(VERSION 3.26.5)

PROJECT(millepedeII LANGUAGES Fortran C CXX)

set(CMAKE_Fortran_COMPILER gcc)

# check for ROOT 
find_package(ROOT COMPONENTS Core Hist Graf Gpad RIO QUIET )

# Arguments for CMake
option(DEBUG "Debug" OFF)
option(SUPPORT_READ_C "read C files" ON)
option(SUPPORT_ZLIB "support Z lib for file compression" ON)
option(SUPPORT_C_RFIO "Use C RFIO" OFF)
option(SUPPORT_OPENMP "Support OPENMP" ON)
option(LAPACK_OPENBLAS "Use LAPACK with OpenBLAS" ON)
option(LAPACK_MKL "Use LAPACK with MKL. Define MKL_DIR if installed in a nonstandard location." OFF)
option(PARDISO "Support Intel oneMKL PARDISO. Define MKL_DIR if installed in a nonstandard location." OFF)
option(SCOREP "Add SCORE-P instrumentation. Define SCOREP_DIR if installed in a nonstandard location." OFF)

# Check dependencies
if (SUPPORT_ZLIB)
    find_package(ZLIB  REQUIRED)
endif()

if (SUPPORT_C_RFIO)
    find_package(SHIFT  REQUIRED)
endif()

set(useLP_OPENBLAS false)
if (LAPACK_OPENBLAS)
    find_package( LAPACK REQUIRED )
    set(useLP_OPENBLAS true)
endif()

if (LAPACK_MKL)
    find_package( MKL REQUIRED )
    set(useLP_OPENBLAS false)
endif()



# Destination targets for build results
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)


if (DEBUG)
    message("Setting flags for debug build")    
    set(CMAKE_Fortran_FLAGS "-Wall -fautomatic -fno-backslash -O3  -g -fcheck=all -fbacktrace -cpp -ffpe-trap=zero,overflow,underflow")
    set(CMAKE_C_FLAGS "-Wall -Wno-unused-function -O3")
else()
    set(CMAKE_Fortran_FLAGS "-Wall -fautomatic -fno-backslash -O3 -cpp")
    set(CMAKE_C_FLAGS "-Wall -Wno-unused-function -O3")
endif()

# add our main targets 

# ===== Mille ======= 
add_library(Mille SHARED)
target_sources(Mille 
            PRIVATE Mille.cc 
            PUBLIC 
                FILE_SET HEADERS 
                FILES Mille.h)

# ===== Pede ======== 
add_executable(pede)
target_sources(pede
  PUBLIC
    mpdef.f90 mpdalc.f90 mpmon.f90 mpbits.f90 mpqldec.f90 mptest1.f90 mptest2.f90 mille.f90 mpnum.f90 mptext.f90 mphistab.f90 
	minresDataModule.f90 minresModule.f90 minresqlpDataModule.f90 minresqlpBlasModule.f90 minresqlpModule.f90 
    randoms.f90 vertpr.f90 linesrch.f90 Dbandmatrix.f90 mpmod.f90 pede.f90
)
target_link_libraries(pede PUBLIC gfortran m )

# Special setup implementing build options

# ===== LAPACK with MKL ===== 
if (LAPACK_MKL)
    message ("Building with LAPACK support via MKL")
    target_include_directories(pede PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(pede PUBLIC mkl_rt)
    target_compile_definitions(pede PUBLIC LAPACK64=\"MKL\")
    target_link_directories(pede PUBLIC ${MKL_ROOT}/lib/${MKL_ARCH})
    if (PARDISO)
        message ("Using Intel oneMKL PARDISO")
        target_compile_definitions(pede PUBLIC PARDISO)
        target_include_directories(pede PUBLIC ${MKL_ROOT}/include)
    endif()
endif()
# ===== LAPACK with OpenMP ===== 
if (useLP_OPENBLAS)
    message ("Building with LAPACK support via OpenBLAS")
    target_link_libraries(pede PUBLIC openblasp64)
    target_compile_definitions(pede PUBLIC LAPACK64=\"OPENBLAS\")
endif()
# ===== Reading of C files ===== 
if (SUPPORT_READ_C)
    message ("Building with C-file support")
    add_library(readc )
    target_sources(readc PUBLIC readc.c)
    target_compile_definitions(pede PUBLIC READ_C_FILES)
    target_link_libraries(pede PUBLIC readc)
    # ===== LibZ for file compression ===== 
    if (SUPPORT_ZLIB)
        message ("Building with libz support")
        target_compile_definitions(readc PUBLIC USE_ZLIB)
        target_link_libraries(readc PUBLIC z)
    endif()
endif()
# ==== OpenMP ==== 
if (SUPPORT_OPENMP)
    message ("Building with openMP support")
    target_compile_options(pede PUBLIC -fopenmp)
    target_link_libraries(pede PUBLIC gomp)
endif()

# install headers for Mille 
install(
  TARGETS Mille
  LIBRARY DESTINATION lib
  FILE_SET HEADERS 
  DESTINATION include
)
# add SCORE-P instrumentation 
if (SCOREP)
    if (NOT DEFINED SCOREP_DIR)
        message ( FATAL_ERROR "Please set SCOREP_DIR to the base path of your SCORE-P installation to enable SCORE-P profiling. ")  
    else()
        message ("Adding instrumentation for SCORE-P profiling") 
        set(CMAKE_C_COMPILER "${SCOREP_DIR}/bin/scorep-gcc")
        set(CMAKE_Fortran_COMPILER "${SCOREP_DIR}/bin/scorep-gcc")
        target_include_directories(pede PUBLIC ${SCOREP_DIR}/include)
    endif()
endif()


# install pede executable 
install(
  TARGETS pede
  DESTINATION bin
)

# Standalone histogrammer (if we have ROOT)
if (${ROOT_FOUND})
    set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
    message ("Found ROOT - will build histogram reader") 
    add_executable(readPedeHists readPedeHists.cpp)
    target_link_libraries(readPedeHists PRIVATE ${ROOT_LIBRARIES})
    target_include_directories(readPedeHists PRIVATE ${ROOT_INCLUDE_DIRS})
    install (
        TARGETS readPedeHists
        DESTINATION bin
    )
endif()

# finally, prepare a setup script for the user 
# to get pede into their environment 
configure_file(
    ${millepedeII_SOURCE_DIR}/cmake/setup.sh.cmake.in
    ${millepedeII_BINARY_DIR}/mp2setup.sh @ONLY)

install(FILES ${millepedeII_BINARY_DIR}/mp2setup.sh
         DESTINATION . )

if (SCOREP AND DEFINED SCOREP_DIR)
    message("=============================
Please make sure to 
    export SCOREP_WRAPPER_INSTRUMENTER_FLAGS=\"--pomp --user --thread=omp\"
before invoking make with SCORE-P instrumentation
=============================" )
endif()
