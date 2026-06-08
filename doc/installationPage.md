 \page installation_page Installation and compilation

## Requirements

### Fortran standards

Since version V04-13-00 at least fortran 2003 is required. The code complies with fortran 2023 (`gcc -std=f2023 -fall-intrinsics`).

### Mille

The package depends on [Mille]((https://gitlab.desy.de/millepede/Mille) to read input binaries.
It will search for an existing installation using `MILLE_DIR` and the `CMAKE_PREFIX_PATH`. 
If no installation is found, the recommended release will be automatically downloaded and compiled during the build of Millepede-II. 

### Optional dependencies

Optionally, the package can be compiled to use 
- [OpenMP](https://www.openmp.org/) for parallelisation, 
- [OpenBLAS](http://www.openmathlib.org/OpenBLAS/) or [Intel oneAPI MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) to use LAPACK routines for optimised solving, 
- [Intel MKL PARDISO](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-0/pardiso.html) to further improve the solution of sparse systems, 
- [ROOT](https://root.cern/) to visualise monitoring histograms and use ROOT-based I/O 

Instrumentation for profiling using [SCORE-P](https://www.vi-hps.org/projects/score-p/overview/overview.html) is implemented and can also be enabled during the build if desired. 

For the LAPACK libraries, 64-bit integer interface are required. 
By default, the build routine tries to find LAPACK on the system and will enable the corresponding features if successful.
If multiple LAPACK installations are found, by default MKL is preferred over OpenBLAS which takes priority over all others. 
See the configuration options below for how to enforce a strict requirement of a certain LAPACK flavour. 

### Documentation

To build the documentation, [doxygen](https://www.doxygen.nl/) version 1.9.8 or later is required to correctly generate internal links between the pages. 
This is for example available in the `gitlab.desy.de:5555/millepede/millepede-ii/mp2dev:latest` docker image.


## Installation

 To install **Millepede** (on a linux system):

### Download the sources 
Download the software package from the DESY \c gitlab server to
    \a target directory, e.g. (shallow clone):

         git clone --depth 1 --branch V05-01-05 \
             https://gitlab.desy.de/millepede/millepede-ii.git target

Then compile using the following steps:
- Create a `build` folder, outside the source folder
- navigate into `build` and invoke `cmake <path_to_source>`. Here you can set additional options (see next sub-section) to configure the build.  
- then compile the sources by calling `make` followed by `make install`. 
- If you wish to compile the documentation, additionally do `make doc`. 

#### Configuring the CMake build
You can pass additional options / flags to CMake: 
- `-DCMAKE_INSTALL_PREFIX=<...>`: Set the install location for the binaries. The default is a subfolder `MillePedeInstall` within your build folder. 
- `-DDEBUG=on`: Enable debug build 
- `-SUPPORT_LAPACK=off`: Disable (default: on) LAPACK support. 
- `-DLAPACK_OPENBLAS=on`: Require (default: off) the use of OpenBLAS for LAPACK. Will throw a fatal error if OpenBLAS is not found. 
- `-DLAPACK_MKL=on`: Require (default: off) the use of Intel MKL for LAPACK. Priority over OpenBLAS if that flag is set as well. If you have a nonstandard installation, set `MKL_DIR` to point to your MKL installation (the standard intel oneAPI setup scripts will do this for you).  Will throw a fatal error if MKL is not found. 
- `-DPARDISO=on`: Enable (default: off) support for Intel PARDISO. Depends on MKL (and will enable the corresponding hard requirement automatically).  
- `-DSCOREP=on`: Enable (default: off) instrumentation for profiling with SCORE_P. If you have a nonstandard installation, set `SCOREP_DIR` to point to your Score-P installation. 
- `-DSUPPORT_OPENMP=off`: Disable (default: on) support for parallelisation with OpenMP.
- `-DSUPPORT_ROOT=off`: Disable (default: on) support for ROOT for file access and histogramming, if ROOT is found. 


##### Before V05-01-00 
Up to and including V05-00-00, the build options were slightly different: 
- `-DCMAKE_INSTALL_PREFIX=<...>`: Set the install location for the binaries. The default is a subfolder `MillePedeInstall` within your build folder. 
- `-DDEBUG=on`: Enable debug build 
- `-DLAPACK_OPENBLAS=off`: Disable (default: on) LAPACK support using OpenBLAS 
- `-DLAPACK_MKL=on`: Enable (default: off) LAPACK support using Intel MKL. Disables OpenBLAS if set. If you have a nonstandard installation, set `MKL_DIR` to point to your MKL installation (the standard intel oneAPI setup scripts will do this for you). 
- `-DPARDISO=on`: Enable (default: off) support for Intel PARDISO. 
- `-DSCOREP=on`: Enable (default: off) instrumentation for profiling with SCORE_P. If you have a nonstandard installation, set `SCOREP_DIR` to point to your Score-P installation. 
- `-DSUPPORT_OPENMP=off`: Disable (default: on) support for parallelisation with OpenMP.
- `-DSUPPORT_ROOT=off`: Disable (default: on) support for ROOT for file access and histogramming, if ROOT is found. 

#### Bonus: Setup script
The build will place a `mp2setup.sh` script in your chosen installation folder.
If this is not a folder on the default path, you can call this script to place `pede` and the standalone programs on your `PATH`, and to make the libraries visible. 
This removes the need to point by hand to the `pede` location. 
If the installation automatically downloaded `Mille` for you, the install script will also handle the corresponding components.  



### Testing your build

 Optionally check the installation by running the simple test case:

         ./pede -t

    This will create (and use) the necessary text and binary files.

 Alternatively tarballs can be found [here](http://www.desy.de/~kleinwrt/MP2/tar). 

## Alternative: Docker images 

- A pre-made docker image containing the standard (LAPACK-OpenBLAS, OpenMP) version of `pede` built on a minimal `EL9` stack is available under the DESY docker image registry: 
```
docker pull gitlab.desy.de:5555/millepede/millepede-ii/mp2runtime:latest
```  
This allows to run the program in a containerised environment if desired. 

- A second image, `gitlab.desy.de:5555/millepede/millepede-ii/mp2dev_base:latest` provides a development environment for building and testing millepede-II in its default configuration. 
