 \page installation_page Installation and compilation

## Requirements

### Fortran standards

Since version V04-13-00 at least fortran 2003 is required. The code complies with fortran 2023 (`gcc -std=f2023 -fall-intrinsics`).

### Optional dependencies

Optionally, the package can be compiled to use 
- [OpenMP](https://www.openmp.org/) for parallelisation, 
- [OpenBLAS](http://www.openmathlib.org/OpenBLAS/) or [Intel oneAPI MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) to use LAPACK routines for optimised solving, 
- [Intel MKL PARDISO](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-0/pardiso.html) to further improve the solution of sparse systems, 
- [ROOT](https://root.cern/) to visualise monitoring histograms and use ROOT-based I/O 

Instrumentation for profiling using [SCORE-P](https://www.vi-hps.org/projects/score-p/overview/overview.html) is implemented and can also be enabled during the build if desired. 

### Documentation

To build the documentation, [doxygen](https://www.doxygen.nl/) version 1.9.8 or later is required to correctly generate internal links between the pages. 
This is for example available in the `gitlab.desy.de:5555/millepede/millepede-ii/mp2dev:latest` docker image.


## Installation

 To install **Millepede** (on a linux system):

### Download the sources 
Download the software package from the DESY \c gitlab server to
    \a target directory, e.g. (shallow clone):

         git clone --depth 1 --branch V04-19-03 \
             https://gitlab.desy.de/millepede/millepede-ii.git target

Then compile using either of the following options: 

 ### Classical Makefile
Directly invoke `make` to create the **Pede** executable (in \a target directory):
```
            make pede
```
   In this workflow, you can configure the build by directly editing the `Makefile` 

### CMake-based build
Alternatively, you can build the package using `CMake`.
In this setup, automatic dependency detection will be attempted. 
- Create a `build` folder, outside the source folder
- navigate into `build` and invoke `cmake <path_to_source>`   
- then compile the sources by calling `make` followed by `make install`. 

#### Configuring the CMake build
You can pass additional options / flags to CMake: 
- `-DCMAKE_INSTALL_PREFIX=<...>`: Set the install location for the binaries. This is needed if you have no root permissions, as the default is a system-wide installation. 
- `-DDEBUG=on`: Enable debug build 
- `-DLAPACK_OPENBLAS=off`: Disable (default: on) LAPACK support using OpenBLAS 
- `-DLAPACK_MKL=on`: Enable (default: off) LAPACK support using Intel MKL. Disables OpenBLAS if set. If you have a nonstandard installation, set `MKL_DIR` to point to your MKL installation (the standard intel oneAPI setup scripts will do this for you). 
- `-DPARDISO=on`: Enable (default: off) support for Intel PARDISO. 
- `-DSCOREP=on`: Enable (default: off) instrumentation for profiling with SCORE_P. If you have a nonstandard installation, set `SCOREP_DIR` to point to your Score-P installation. 
- `-DSUPPORT_OPENMP=off`: Disable (default: on) support for parallelisation with OpenMP.
- `-DSUPPORT_ZLIB=off`: Disable (default: on) support for lib-z compression for C-binaries. 
- `-DSUPPORT_READ_C=off`: Disable (default: on) support for reading C-files. 
- `-DSUPPORT_C_RFIO=on`: Enable (default: off) use of C RFIO for binary reading. 
- `-DSUPPORT_ROOT=off`: Disable (default: on) support for ROOT for file access and histogramming, if ROOT is found. 

#### Bonus: Setup script
The build will place a `mp2setup.sh` script in your chosen installation folder.
If this is not a folder on the default path, you can call this script to place `pede` and the standalone programs on your `PATH`, and to make the libraries visible. 
This removes the need to point by hand to the `pede` location. 



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
