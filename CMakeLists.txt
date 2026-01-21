CMAKE_MINIMUM_REQUIRED(VERSION 3.26.5)


# discover the compilers from PATH (when system install is overridden)
find_program(gf_loc gfortran)
find_program(gcc_loc gcc)
find_program(gpp_loc g++)

set(CMAKE_Fortran_COMPILER ${gf_loc})
set(CMAKE_C_COMPILER ${gcc_loc})
set(CMAKE_CXX_COMPILER ${gpp_loc})


PROJECT(millepedeII LANGUAGES Fortran C CXX)
set(millepedeII_VERSION 05.00.0)
include(CTest)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE "${CMAKE_BINARY_DIR}/MillePedeInstall")
  set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT False)
endif()


find_package(Mille QUIET )
if (${Mille_FOUND})
    message ("Found existing Mille installation")
else()
    message ("No existing Mille installation found - will clone and build the recommended version")
    include ("cmake/getMille.cmake")
endif()

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
    find_package(ROOT CONFIG QUIET )
    if (${ROOT_FOUND})
        message("Building ROOT related components")
        set(build_ROOT_comps true)
    endif()
endif()


if (LAPACK_MKL)
    find_package( MKL REQUIRED )
    set(useLP_OPENBLAS false)
endif()



# Destination targets for build results
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)



if (DEBUG)
    message("Setting flags for debug build")    
    set(CMAKE_Fortran_FLAGS "-Wall -fautomatic -fno-backslash -O3  -g -fcheck=all -fbacktrace -cpp -ffpe-trap=zero,overflow,underflow")
    set(CMAKE_C_FLAGS "-Wall -Wno-unused-function -O3")
else()
    set(CMAKE_Fortran_FLAGS "-Wall -fautomatic -fno-backslash -O3 -cpp")
    set(CMAKE_C_FLAGS "-Wall -Wno-unused-function -O3")
endif()

# add our main targets 

# ===== Pede ======== 
add_executable(pede)
target_sources(pede
  PRIVATE
    src/mpdef.f90 src/mpdalc.f90 src/mpmon.f90 src/mpbits.f90 src/mpqldec.f90 src/mptest1.f90 src/mptest2.f90 src/mille.f90 src/mpnum.f90 src/mptext.f90 src/mphistab.f90 
	src/minresDataModule.f90 src/minresModule.f90 src/minresqlpDataModule.f90 src/minresqlpBlasModule.f90 src/minresqlpModule.f90 
    src/randoms.f90  src/vertpr.f90 src/linesrch.f90 src/Dbandmatrix.f90 src/mpmod.f90 src/mppar.f90 src/mplapack.f90 src/pede.f90
)
target_link_libraries(pede PUBLIC gfortran m MilleC)
add_dependencies(pede  MilleF)
target_link_libraries(pede INTERFACE MilleF)
message ("Will include ${Mille_MODULE_DIRS}")
target_include_directories(pede PRIVATE ${Mille_MODULE_DIRS})

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

if (build_ROOT_comps)
    message ("Found ROOT - will build histogram reader and ROOT I/O lib") 

    set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
    # histogrammer 
    add_executable(readPedeHists)
    target_sources(readPedeHists PRIVATE util/readPedeHists.cpp)
    target_link_libraries(readPedeHists PUBLIC ROOT::Core ROOT::Tree ROOT::RIO ROOT::Graf ROOT::Gpad ROOT::Hist)
    install (
        TARGETS readPedeHists
        EXPORT millepedeII-targets
        DESTINATION bin
    )

endif()

# ==== OpenMP ==== 
if (SUPPORT_OPENMP)
    message ("Building with openMP support")
    target_compile_options(pede PUBLIC -fopenmp)
    target_link_libraries(pede PUBLIC gomp)
endif()


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
  EXPORT millepedeII-targets
  DESTINATION bin
)

# install python modules 
install (FILES python/comparePedeResults.py  python/__init__.py python/compareResults.py python/tinypede.py python/readMilleBinary.py
        DESTINATION python/Millepede)


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

install (EXPORT millepedeII-targets
        DESTINATION cmake 
)

# add a target to generate API documentation with Doxygen
FIND_PACKAGE(Doxygen)

if(DOXYGEN_FOUND)

    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)

    ADD_CUSTOM_TARGET(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )

ELSE()
    MESSAGE( STATUS "Doxygen not installed or found -- API documentation can't be generatred!" )

ENDIF(DOXYGEN_FOUND)

# Create and install package configuration and version files.
configure_file(
${millepedeII_SOURCE_DIR}/cmake/millepedeII-config.cmake.in
${millepedeII_BINARY_DIR}/cmake/millepedeIIConfig.cmake @ONLY)

configure_file(
${millepedeII_SOURCE_DIR}/cmake/millepedeII-config-version.cmake.in
${millepedeII_BINARY_DIR}/cmake/millepedeIIConfig-version.cmake @ONLY)

install(FILES ${millepedeII_BINARY_DIR}/cmake/millepedeIIConfig.cmake
        ${millepedeII_BINARY_DIR}/cmake/millepedeIIConfig-version.cmake
        DESTINATION cmake/)


# Unit tests 
subdirs(test)
