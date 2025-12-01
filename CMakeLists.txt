CMAKE_MINIMUM_REQUIRED(VERSION 3.26.5)


# discover the compilers from PATH (when system install is overridden)
find_program(gf_loc gfortran)
find_program(gcc_loc gcc)
find_program(gpp_loc g++)

set(CMAKE_Fortran_COMPILER ${gf_loc})
set(CMAKE_C_COMPILER ${gcc_loc})
set(CMAKE_CXX_COMPILER ${gpp_loc})


PROJECT(millepedeII LANGUAGES Fortran C CXX)
include(CTest)

# Arguments for CMake
option(DEBUG "Debug" OFF)
option(SUPPORT_READ_C "read C files" ON)
option(SUPPORT_ROOT "build ROOT file I/O and root plotting scripts" ON)
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

set(build_ROOT_comps false)
if(SUPPORT_ROOT)
    # check for ROOT 
    if (DEFINED ENV{ROOTSYS})
        # follow ROOTSYS if it has been defined 
        set (ROOT_DIR $ENV{ROOTSYS}/cmake) 
        message ("Taking ROOT from ROOTSYS: $ENV{ROOTSYS}")
    endif()
    find_package(ROOT COMPONENTS Core Hist Graf Gpad RIO Net QUIET PATHS ${ROOT_DIR} )
    if (${ROOT_FOUND})
        set(build_ROOT_comps true)
    endif()
endif()

set(build_CTOROOT false)
if (build_ROOT_comps AND SUPPORT_READ_C)
    set(build_CTOROOT true)
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
    randoms.f90  vertpr.f90 linesrch.f90 Dbandmatrix.f90 mpmod.f90 mppar.f90 mplapack.f90 pede.f90
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
    add_library(readc SHARED)
    target_sources(readc 
                PRIVATE readc.c
                PUBLIC 
                    FILE_SET HEADERS 
                    FILES readc.h)
    target_link_libraries(pede PUBLIC readc)
    target_compile_definitions(pede PUBLIC READ_C_FILES)
    # target_link_libraries(pede PUBLIC readc)
    # ===== LibZ for file compression ===== 
    if (SUPPORT_ZLIB)
        message ("Building with libz support")
        target_compile_definitions(readc PUBLIC USE_ZLIB)
        target_link_libraries(readc PUBLIC z)
    endif()
    install(
        TARGETS readc
        LIBRARY DESTINATION lib
        FILE_SET HEADERS 
        DESTINATION include
    )
endif()

if (build_ROOT_comps)
    message ("Found ROOT - will build histogram reader and ROOT I/O lib") 

    set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})

    # Read ROOT files
    add_library(readROOT SHARED)
    target_sources(readROOT 
                PRIVATE readROOT.cc 
                PUBLIC 
                    FILE_SET HEADERS 
                    FILES readROOT.h)
    target_link_libraries(readROOT PUBLIC ${ROOT_LIBRARIES})
    target_include_directories(readROOT PUBLIC ${ROOT_INCLUDE_DIRS})
    target_compile_definitions(pede PUBLIC READ_ROOT_FILES)
    target_link_libraries(pede PUBLIC readROOT)
    install(
        TARGETS readROOT
        LIBRARY DESTINATION lib
        FILE_SET HEADERS 
        DESTINATION include
    )

    # enable ROOT support in Mille
    target_compile_definitions(Mille PUBLIC SUPPORT_ROOTIO)
    target_link_libraries(Mille PUBLIC ${ROOT_LIBRARIES})
    target_include_directories(Mille PUBLIC ${ROOT_INCLUDE_DIRS})

    # histogrammer 
    add_executable(readPedeHists readPedeHists.cpp)
    target_link_libraries(readPedeHists PRIVATE ${ROOT_LIBRARIES})
    target_include_directories(readPedeHists PRIVATE ${ROOT_INCLUDE_DIRS})
    install (
        TARGETS readPedeHists
        DESTINATION bin
    )

    # C-to-ROOT conversion 
    if (build_CTOROOT)
        add_executable(cToRoot cToRoot.cc)
        target_link_libraries(cToRoot PUBLIC readc ${ROOT_LIBRARIES})
        target_include_directories(cToRoot PRIVATE ${ROOT_INCLUDE_DIRS})
        install (
            TARGETS cToRoot
            DESTINATION bin
        )
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
    message ("Adding instrumentation for SCORE-P profiling") 
    # first, we try to find a system-wide installation. 
    # Score-P does not come with a CMake config, so look for the main 
    # program 
    find_program (
          found_SCOREP
          NAMES scorep 
          PATHS ENV PATH 
          ) 
    # if we find no installation on PATH, we assume there is a custom setup. 
    # In this case, the user should supply "SCOREP_DIR" to specify where to look. 
    if (NOT found_SCOREP) 
        if (NOT DEFINED SCOREP_DIR)
            message ( FATAL_ERROR "Please set SCOREP_DIR to the base path of your SCORE-P installation to enable SCORE-P profiling. ")  
        else()
            set(CMAKE_C_COMPILER "${SCOREP_DIR}/bin/scorep-gcc")
            set(CMAKE_Fortran_COMPILER "${SCOREP_DIR}/bin/scorep-gcc")
            target_include_directories(pede PUBLIC ${SCOREP_DIR}/include)
        endif()
    else()
        set(CMAKE_C_COMPILER "scorep-gcc")
        set(CMAKE_Fortran_COMPILER "scorep-gcc")
    endif()
endif()


# install pede executable 
install(
  TARGETS pede
  DESTINATION bin
)
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


# Unit tests 
subdirs(test)
