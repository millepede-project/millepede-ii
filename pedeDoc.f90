
!> \file
!! Millepede II program, documentation.
!!
!! \author Volker Blobel, University Hamburg, 2005-2009 (initial Fortran77 version)
!! \author Gero Flucke, University Hamburg (support of C-type binary files)
!! \author Claus Kleinwort, DESY (maintenance and developement)
!!
!! \copyright
!! Copyright (c) 2009 - 2024 Deutsches Elektronen-Synchroton,
!! Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
!! This library is free software; you can redistribute it and/or modify
!! it under the terms of the GNU Library General Public License as
!! published by the Free Software Foundation; either version 2 of the
!! License, or (at your option) any later version. \n\n
!! This library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Library General Public License for more details. \n\n
!! You should have received a copy of the GNU Library General Public
!! License along with this program (see the file COPYING.LIB for more
!! details); if not, write to the Free Software Foundation, Inc.,
!! 675 Mass Ave, Cambridge, MA 02139, USA.
!!

!> \mainpage Overview
!!
!! \section intro_sec Introduction
!! In certain least squares fit problems with a very large number of parameters
!! the set of parameters can be divided into two classes, global and local parameters.
!! Local parameters are those parameters which are present only in subsets of the
!! data. Detector alignment and calibration based on track fits is one of the problems,
!! where the interest is only in optimal values of the global parameters, the
!! alignment parameters. The method, called Millepede, to solve the linear least
!! squares problem with a simultaneous fit of all global and local parameters,
!! irrespectively of the number of local parameters, is described in the draft manual.
!! (Correlated measurements need to be transformed into independent measurements
!! by diagonalization of their covariance matrix.)
!!
!! The Millepede method and the initial implementation has been
!! developed by [V. Blobel](http://www.desy.de/~blobel) from he University of Hamburg.
!! Meanwhile the code is maintained at DESY by the statistics tools group of the
!! analysis center of the Helmholtz [Terascale](https://terascale.de) alliance
!! using [GitLab](https://about.gitlab.com) ([code and wiki](https://gitlab.desy.de)).
!!
!! The Millepede II software is provided by DESY under the terms of the
!! [LGPLv2 license](http://www.gnu.org/licenses/old-licenses/lgpl-2.0-standalone.html).
!!
!! \section install_sec Installation
!! To install **Millepede** (on a linux system):
!! 1. Download the software package from the DESY \c gitlab server to
!!    \a target directory, e.g. (shallow clone):
!!
!!         git clone --depth 1 --branch V04-17-07 \
!!             https://gitlab.desy.de/claus.kleinwort/millepede-ii.git target
!!
!! 2. Create **Pede** executable (in \a target directory):
!!
!!         make pede
!!
!! 3. Optionally check the installation by running the simple test case:
!!
!!         ./pede -t
!!
!!    This will create (and use) the necessary text and binary files.
!!
!! Alternatively tarballs can be found [here](http://www.desy.de/~kleinwrt/MP2/tar).
!!
#include "pedeChanges.f90"
!!
!! \section tools_sec Tools
!! The subdirectory \c tools contains some useful scripts:
!! * \c readMilleBinary.py: Python script to read binary files and print
!!   records in text form.
!! * \c compareResults.py: Python3 script to compare result files (<tt>millepede.res</tt>).
!! * \c readPedeHists.C: ROOT script to read and convert the **Millepede**
!!   histogram file <tt>millepede.his</tt>.
!! * \c lapack: Test programs to print LAPACK (library) configuration (MKL, OpenBLAS).
!! * \c tinypede.py: **Pede** implementation in python3 with basic functionality for
!!   illustration or testing with small problems.
!!
!! \section julia_sec Julia
!! The subdirectory \c julia contains some of the above tools reimplemented in [Julia](https://julialang.org):
!! * \c readMilleBinary.jl
!! * \c tinypede.jl
!!
!! \section details_sec Details
!!
!! Detailed information is available at:
!!
!! \subpage draftman_page
!!
!! \subpage changes_page
!!
!! \subpage option_page
!!
!! \subpage exit_code_page
!!
!! \subpage troubleshooting_page
!!
!! \subpage test_brlf_page "Example"
!!
!! \section Contact
!!
!! For information exchange the **Millepede** mailing list
!! anacentre-millepede2@desy.de should be used.
!!
!! \section Legacy
!! The subdirectory \c legacy contains the original \ref millepede1.f90
!! "Millepede-I" implementation from Volker
!! Blobel (2000) and a Millepede-I to Millepede-II \ref mp1to2.f90 "interface"
!! (creating Millepede-II input files from Millepede-I calls, developed for COMPASS at CERN).
!!
!! \section ref_sec References
!!
!! 1. A New Method for the High-Precision Alignment of Track Detectors,
!!    Volker Blobel and Claus Kleinwort, Proceedings of the Conference on
!!    Adcanced Statistical Techniques in Particle Physics, Durham, 18 - 22 March 2002,
!!    Report DESY 02-077 (June 2002) and
!!    [hep-ex/0208021](http://arxiv.org/abs/hep-ex/0208021)
!! 2.  Alignment Algorithms, V. Blobel,
!!    [Proceedings](http://cdsweb.cern.ch/search?p=reportnumber%3ACERN-2007-004)
!!    of the LHC Detector Alignment Workshop, September 4 - 6 2006, CERN
!! 3. Software alignment for Tracking Detectors, V. Blobel,
!!    NIM A, 566 (2006), pp. 5-13,
!!    [doi:10.1016/j.nima.2006.05.157](http://dx.doi.org/10.1016/j.nima.2006.05.157)
!! 4. A new fast track-fit algorithm based on broken lines, V. Blobel,
!!    NIM A, 566 (2006), pp. 14-17,
!!    [doi:10.1016/j.nima.2006.05.156](http://dx.doi.org/10.1016/j.nima.2006.05.156)
!! 5. Millepede 2009, V. Blobel,
!!    [Contribution](https://indico.cern.ch/conferenceOtherViews.py?view=standard&confId=50502)
!!    to the 3rd LHC Detector Alignment Workshop, June 15 - 16 2009, CERN
!! 6. General Broken Lines as advanced track fitting method, C. Kleinwort,
!!    NIM A, 673 (2012), pp. 107-110,
!!    [doi:10.1016/j.nima.2012.01.024](http://dx.doi.org/10.1016/j.nima.2012.01.024)
!! 7. Volker Blobel und Erich Lohrmann, Statistische und numerische Methoden der
!!    Datenanalyse, Teubner Studienb&uuml;cher, B.G. Teubner, Stuttgart, 1998.
!!    [Online-Ausgabe](http://www.desy.de/~blobel/eBuch.pdf).
!! 8. [Systems Optimization Laboratory](http://web.stanford.edu/group/SOL/software/minres),
!!    Stanford University;\n
!!    C. C. Paige and M. A. Saunders (1975),
!!    Solution of sparse indefinite systems of linear equations,
!!    SIAM J. Numer. Anal. 12(4), pp. 617-629.
!! 9. [Systems Optimization Laboratory](http://web.stanford.edu/group/SOL/software/minresqlp),
!!    Stanford University;\n
!!    Sou-Cheng Choi, Christopher Paige, and Michael Saunders,
!!    MINRES-QLP: A Krylov subspace method for indefinite or singular
!!    symmetric systems, SIAM Journal of Scientific Computing 33:4, 1810-1836, 2011,
!!    [doi:10.1137/100787921](http://dx.doi.org/10.1137/100787921)
!!

!> \page changes_page Major changes
!! Major changes with respect to the \ref draftman_page "draft manual".
!! \tableofcontents
!!
!! \section ch-entries-ext Extended entries cut
!! The \ref cmd-entries "entries" cut has now three arguments:
!! 1. \ref mpmod::mreqenf "mreqenf": Minimum required (number of) \ref an-entries "entries"
!!    in binary files for global parameters to be *variable*.
!! 2. \ref mpmod::mreqena "mreqena": Minimum required (number of) entries
!!    (with record/track) \ref localfit-rejection "accepted" after the local (track) fit for
!!    variable global parameters. Due to the tigthening of the \ref cmd-chisqcut "chisqcut" during
!!    the internal pede iterations the number of *accepted* entries can get far lower than the
!!    number of entries in the binary files(&ge; **mreqenf**). If the number of entries gets lower
!!    than the number of variable global parameters determined from the related measurements
!!    there will be no more proper solution (for this part of the global linear equation
!!    system). Therefore *severe* warnings are issued
!!    in case the number of accepted entries is below **mreqena** for any variable
!!    global parameter in any iteration.
!!    As the gobal matrix is usually constructed only in the first iteration it probably will
!!    not get singular in this case, but the solutions for the affected parameters are not
!!    independent anymore.
!! 3. **iarg3** (\ref mpmod::iteren "iteren"=mreqenf*iarg3): Iteration of entries cut.
!!    It can be problematic if different global parameters related to the same measurements
!!    get different (number of) entries. If some are below and some are above **mreqenf**
!!    some are fixed and some are variable which may not make sense. An example is the
!!    rigid body alignment of a planar silicon pixel sensor in the local (uvw) system. For the
!!    two indendent measurements (in u and v) the corresponding alignment offsets (du, dv)
!!    appear each only for *one* measurement. The offset perpendicular to the sensor plane (dw)
!!    and the rotations affect *both* measurement. Therefore the offsets (du, dv) have a factor
!!    two smaller number of entries as the other four. To handle this situation the entries cut
!!    (on the binary files) can be iterated: If in a measurement some global parameters are fixed
!!    (due to entries < **mreqenf**) for the other parameters the cut value is increased by
!!    the constant factor **iarg3**.
!!    For the axample above a factor two would be appropriate (all six parameters fixed now).
!!
!!    Another option to avoid this problem is to not count the entries on the level of equations
!!    (measurements) but on the level of records (tracks): \ref cmd-countrecords "countrecords"
!!
!!    A third option is to produce the binary files without zero supression (of global derivatives).
!! \section ch-methods Solution methods
!! The following methods to obtain the solution \f$\Vek{x}\f$ from a
!! linear equation system \f$\Vek{A}\cdot\Vek{x}=\Vek{b}\f$ are implemented:
!! \subsection ch-inv Inversion
!! The solution and the covariance matrix \f$\Vek{A}^{-1}\f$ are obtained by
!! \ref an-inv "inversion" of \f$\Vek{A}\f$.
!! Available are the value, error and global correlation for all global parameters.
!! The matrix inversion \ref sqminl "routine" has been \ref ch-openmp "parallelized"
!! and can be used for up to several 10000 parameters.
!! \subsection ch-diag Diagonalization
!! The solution and the covariance matrix \f$\Vek{A}^{-1}\f$ are obtained by
!! \ref an-diag "diagonalization" of \f$\Vek{A}\f$.
!! Available are the value, error, global correlation and
!! eigenvalue (and eigenvector) for all global parameters.
!! \subsection ch-minres Minimal Residual Method (MINRES)
!! The solution is obtained by minimizing \f$\Vert\Vek{A}\cdot\Vek{x}-\Vek{b}\Vert_2\f$
!! iteratively and is only approximate. \ref minresmodule::minres "MINRES"  [\ref ref_sec "ref 8"] is a special case of the
!! generalized minimal residual method (\ref an-gmres "GMRES") for symmetric matrices.
!! Preconditioning with a band matrix of zero or finite
!! \ref mpmod::mbandw "bandwidth" is possible.
!! Individual columns \f$\Vek{c_i}\f$ of the covariance matrix can be calculated by
!! solving \f$\Vek{A}\cdot\Vek{c}_i=\Vek{1}_i\f$ where \f$\Vek{1}_i\f$ is the i-th
!! column on the unit matrix.
!! The most time consuming part (\ref avprod "product" matrix times vector per iteration)
!! has been \ref ch-openmp "parallelized".
!! Available are the value for all (and optionally error, global correlation
!! for few) global parameters.
!! \subsection ch-minresqlp Advanced Minimal Residual Method (MINRES-QLP)
!! The \ref minresqlpmodule::minresqlp "MINRES-QLP" implementation [\ref ref_sec "ref 9"]
!! is a MINRES evolution with improved norm estimates and stopping conditions
!! (leading potentially to different numbers of internal iterations).
!! Internally it uses QLP instead of the QR factorization in
!! MINRES which should be numerically superior and allows to find for
!! singular systems the minimal length (pseudo-inverse) solution.
!!
!! The default behavior is to start (the internal iterations) with QR factorization
!! and to switch to QLP if the (estimated) matrix condition exceeds
!! \ref cmd-mrestranscond "mrtcnd". Pure QR or QLP factorization can be enforced
!! by \ref cmd-mresmode "mrmode".
!!
!! \subsection ch-elim-const Elimination of constraints
!! As alternative to the Lagrange multiplier method the solution by elimination
!! has been added for problems with linear equality constraints.
!! A \ref mpqldec::qldec "QL factorization" (with Householder reflections) of the
!! transposed constraints matrix is used to transform to an unconstrained problem.
!! For sparse matrix storage the sparsity of the global matrix is preserved.
!!
!! \subsection ch-mchdec Decomposition
!! The solution is obtained by a root-free Cholesky decomposition (LDLt).
!! The covarinance matrix is *not* being calclulated.
!! The method ia about a factor 2-3 faster than inversion (according to several tests).
!! It is restricted to solution by elimination for problems with linear equality constraints
!! (requires positive definite matrix). Supports block matrix storage for disjoint parameter blocks.
!!
!! \subsection ch-lapack LAPACK
!! For these optional methods the solution is obtained by matrix factorization from an external LAPACK library.
!! There exist (open or proprietary) implementations heavily optimized for specific hardware
!! (and partially multi-threaded)
!! which could easily be an order of magnitude faster than e.g. the custom code used for
!! \ref ch-mchdec "decomposition".
!! Tested has been the [Intel MKL](https://software.intel.com/en-us/intel-mkl) and
!! [OpenBLAS](https://www.openblas.net).
!!
!! The symmetric global matrix can be stored in **packed** triangular (similar to **full** for inversion etc)
!! or **unpacked** quadaratic \ref mpmod::matsto "shape".
!! For unpacked storage the \ref ch-elim-const "elimination of constraints" is
!! implemented with LAPACK (DGEQLF, DORMQL) as alternative.
!!
!! For positive definite matrices (elimination of constraints) a Cholesky factorization is used (DPPTRF/DPOTRF) and
!! for indefinite matrices (Lagrange multipliers for constraints) a Bunch-Kaufman factorization (DSPTRF/DSYTRF).
!! With the option \ref cmd-lapackerr "LAPACKwitherrors" the inverse matrix is calculated from the
!! factorization to provide parameter errors (DPPTRI/DPOTRI, DSPTRI/DSYTRI, potentially slow, single-threaded).
!!
!! \subsection ch-pardiso PARDISO
!! This optional method is using the [Intel oneMKL PARDISO](
!! https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-2/onemkl-pardiso-parallel-direct-sparse-solver-iface.html)
!! shared-memory multiprocessing parallel direct sparse solver. This aims at use cases with really sparse global
!! matrices (e.g. less than 10% nonzero elements).
!!
!! The upper triangle of the (symmetric) global matrix is stored in the 3-array variation of the compressed sparse row format (CSR3).
!! It can be checked for a (fixed size) \b block structure to reduce the storage size (BSR3).
!! (This could be further improved by reordering - not implemeted.)
!!
!! For problems with linear equality constraints the usage of Lagrange multipliers is enforced.
!!
!! For the internal steering parameters \c IPARM(64) the default values are used.
!! They can be modified with the \ref cmd-pardiso "pardiso" keyword
!! (specifying pairs of indices (1-64) and values):
!!
!!          pardiso
!!          ...      ...
!!          index    value
!!          ...      ...
!!
!! Related options:
!!   - \ref cmd-blocksize-pd "blocksizePARDISO", e.g. combinations of factors 2 and 3:
!!
!!           blocksizePARDISO 2 3 4 6 9 12
!!
!!   - \ref cmd-debug-pd "debugPARDISO".
!!
!! \section ch-regul Regularization
!! Optionally a term \f$\tau\cdot\Vert\Vek{x}\Vert\f$ can be added to the objective function
!! (to be minimized) where \f$\Vek{x}\f$ is the vector of global parameters
!! weighted with the inverse of their individual pre-sigma values.
!!
!! \section ch-locfit Local fit
!! In case the \ref par-locfitv "local fit" is a track fit with proper description of multiple
!! scattering in the detector material additional local parameters have to be introduced
!! for each scatterer and solution by *inversion* can get time consuming
!! (~ \f$n_{lp}^3\f$ for \f$n_{lp}\f$ local parameters). For trajectories based on
!! **broken lines** [\ref ref_sec "ref 4,6"] the corresponding matrix \f$\Vek{\Gamma}\f$
!! has a bordered band structure (\f$\Gamma_{ij}=0\f$ for \f$\min(i,j)>b\f$
!! (border size) and \f$|i-j|>m\f$ (bandwidth)). With
!! <i>root-free Cholesky decomposition</i> the time for the solution is linear
!! and for the calculation of \f$\Gamma^{-1}\f$
!! (needed for the construction of the global matrix) quadratic in \f$n_{lp}\f$.
!! The condition of the diagonal matrix from the decomposition (of the band part)
!! can be used to reject ill conditioned cases (see \ref cmd-maxlocalcond).
!! For each local fit the structure of \f$\Vek{\Gamma}\f$ is checked and the faster
!! solution method selected automatically.
!!
!! \section ch-openmp Parallelization
!! The code has been largely parallelized using [OpenMP&tm;](www.openmp.org).
!! This includes the reading of binary files, the local fits, the construction of the
!! sparsity structure and filling of the global matrix and the global fit
!! (except by diagonalization). The number of threads is set by the command
!! \ref cmd-threads.
!!
!! The OpenMP environment can be displayed with:
!!
!!       export OMP_DISPLAY_ENV=TRUE
!!
!! \b Caching. The records are read in blocks into a *read cache* and processed from
!! there in parallel, each record by a single thread. For the filling of the global
!! matrix the (zero-compressed) update matrices (\f$\Vek{\D C}_1+\Vek{\D C}_2\f$ from
!! equations \ref eq-c1 "(15)", \ref eq-c2 "(16)")
!! produced by each local fit are collected in a
!! *write cache*. After processing the block of records this is used to update
!! the global matrix in parallel, each row by a single thread.
!! The total cache size can be changed by the command \ref cmd-cache.
!!
!! \section ch-sparsemat Sparse matrix storage
!! \subsection ch-compression Compression
!! In sparse storage mode for each row the list of column indices (and values) for the
!! non-zero elements are stored. With compression regions of continous column indices
!! are represented by the first index and their number (packed into a single 32bit
!! integer). Compression is selected by the command \ref cmd-compress.
!! In addition rare elements can be neglected (,histogrammed) or stored in single instead
!! of double precision according to the \ref cmd-pairentries command.
!! \subsection ch-pargroup Parameter groups
!! With the implementation of parameter groups (sets of adjacent global parameters (labels)
!! appearing in the binary files *always* together) the sparsity structure addresses
!! not single values anymore but block matrices with sizes according to the contributing
!! parameter groups. Groups with adjacent column ranges are combined into a single block matrix.
!! The offset and first column index in a (block) row is stored.
!!
!! \section ch-gzip Gzipped C binary files
!! The [zlib](zlib.net) can be used to directly read *gzipped* C binary files.
!! In this case reading with multiple threads
!! (each file by single thread) can speed up the decompression.
!!
!! \section ch-weightedfiles Weighted binary files
!! A constant weight can be applied to all records in a binary file:
!!
!!         <file name> -- <weight> ! <comment>
!!
!! \section ch-transf Transformation from FORTRAN77 to Fortran90
!! The **Millepede** source code has been formally transformed from <i>fixed form</i>
!! FORTRAN77 to <i>free form</i> Fortran90 (using TO_F90 by Alan Miller)
!! and (most parts) modernized:
!! - <tt>IMPLICIT NONE</tt> everywhere. Unused variables removed.
!! - \c COMMON blocks replaced by \c MODULEs.
!! - Backward \c GOTOs replaced by proper \c DO loops.
!! - \c INTENT (input/output) of arguments described.
!! - Code documented with doxygen.
!!
!! Unused parts of the code (like the interactive mode) have been removed.
!! The reference compiler for the Fortran90 version is gcc-4.6.2 (gcc-4.4.4 works too).
!!
!! \section ch-memmanage Memory management
!! The memory management for dynamic data structures (matrices, vectors, ..)
!! has been changed from a \ref an-dynal "subdivided" *static* \c COMMON block to
!! *dynamic* (\c ALLOCATABLE) Fortran90 arrays. One **Pede** executable is now
!! sufficient for all application sizes.
!!
!! \section ch-readbuf Read buffer size
!! In the \ref sssec-loop1 "first loop" over all binary files a preset
!! \ref mpmod::ndimbuf "read buffer size" is used. Too large records are skipped,
!! but the maximal record length is still being updated. If any records had to be skipped
!! the read buffer size is afterwards adjusted according to the maximal record length
!! and the first loop is repeated.
!!
!! \section ch-numbin Number of binary files
!! The number of binary files has no hard-coded limit anymore, but is calculated from
!! the steering file and resources (file names, descriptors, ..)
!! are allocated dynamically. Some resources may be limited by the system.
!!
!! \section ch-checkinput Check input mode
!! The **Pede** executable can be run in special modes to only check the input
!! without producing a solution. This can be useful in case of problems.
!! Checked are the global parameters and their (equality) constraints from the steering
!! and binary files. Debug information is written to the standard output and the
!! <tt>millepede.res</tt> text file.
!!
!! 1. Basic checks: Command '\ref cmd-checkinput 1' or command line option '\ref opt-c'
!!
!!   To standard output:
!!   - For each constraint the number of all and of variable global parameters involved.
!!     *Empty* constraints without variable parameters should be removed.
!!   - Constraint group definitions (first and last (sorted (by global labels)) constraint,
!!     first and last global label). Are disjoint (in global labels) sets of constraints.
!!   - Constraint block definitions (first and last constraint group, first and last global label).
!!     Are non overlapping (in global labels) sets of constraints.
!!
!!   To text file:
!!   - For each global parameter the label, value, presigma, number of entries, constraint
!!     group (0 for none) and status (variable or fixed).
!!   - For the first parameter of a parameter group (consecutive parameter (labels) appearing
!!     always togther) the line starts with '<tt>!></tt>'.
!!   - For each constraint group the number of contributing constraints, number of entries
!!     and first and last global label.
!!
!! 2. More checks: Command '\ref cmd-checkinput 2' or command line option '\ref opt-C'
!!
!!   Additionally to standard output:
!!   - For each **sorted** constraint the index by appearance (in text file) with number
!!     of first line and first and last global label.
!!   - For each constraint group the number of constraints and the rank of the corresponding
!!     block of the product matrix. Must be equal!
!!
!!   Additionally to text file:
!!   - For each global parameter the appearance (in binary files) statistics
!!     (first file and record number, last file and record number and number of files)
!!     and number of paired global parameters (appearing together).
!!   - For each constraint group the label ranges of (not contributing) paired global parameters.
!!   - For each constraint group the appearance statistics (using the contributing global labels).
!!
!! \section ch-parcom Global parameter comments
!! For each global parameter a \ref mpdef::itemclen "brief" **comment** can be defined to annotate the <tt>millepede.res</tt> text file.
!! The syntax is similar as for the "Parameter" or "Constraint" keywords:
!!
!!         Comment
!!         ...     ...
!!         label   text
!!         ...     ...
!!


!> \page option_page List of options and commands
!!
!! \tableofcontents
!!
!! \section sec-opt Command line options:
!! \subsection opt-t1 -t
!! Create text and binary files for \ref mptest1.f90 "wire chamber" test case, set
!! \ref mpmod::ictest "ictest" to 1.
!! \subsection opt-t2 -t=track-model
!! Create text and binary files for \ref mptest2.f90 "silicon strip tracker" test case
!! using \a track-models with different accounting for multiple scattering, set
!! \ref mpmod::ictest "ictest" to 2..6.
!! \subsection opt-s -s
!! Solution is not iterated.
!! Automatically switched on in case of rank deficits for constraints.
!! \subsection opt-f -f
!! Force iterating of solution (in case of rank deficits for constraints).
!! \subsection opt-c -c
!! Check input (binary files, constraints). No solution is determined. (\ref mpmod::icheck "icheck"=1)
!! \subsection opt-C -C
!! Check input (binary files, constraints, appearance). No solution is determined. (\ref mpmod::icheck "icheck"=2)
!!
!! \section sec-cmd Steering file commands:
!! In general the commands are defined by a single line:
!!
!!         keyword   number1  number2  ...
!!
!! For those specifying \ref sssec-parinf "properties" of the global parameters
!! (\a keyword = \c parameter, \c constraint or \c measurement (or \c comment))
!! for each involved global parameter (identified by a \ref an-glolab "label")
!! one additional line follows:
!!
!!         label     number1  number2  ...
!!
!! Default values for the numerical arguments are shown in
!! the command descriptions in '[]'. Missing arguments without default
!! values have no effect.
!!
!! \subsection cmd-bandwidth bandwidth
!! Set band width \ref mpmod::mbandw "mbandw" for
!! \ref minresmodule::minres "MINRES" preconditioner to \a number1 [0]
!! and additional flag \ref mpmod::lprecm "lprecm" to \a number2 [0].
!! \subsection cmd-blocksize-pd blocksizePARDISO
!! Enable for \ref ch-pardiso "PARDISO" storage in BSR3 format
!! (three-array variation of the block compressed sparse row format).
!! Provide (one to ten, increasingly ordered) candidate block sizes. The block size with the smallest
!! memory footprint (of the global matrix in BSR3 format) is selected.
!! \subsection cmd-cache cache
!! Set (read+write) cache size \ref mpmod::ncache "ncache" to \a number1.
!! Define cache size and average fill level.
!! \subsection cmd-cfiles Cfiles
!! Following binaries are C files.
!! \subsection cmd-checkinput checkinput
!! Set check input flag \ref mpmod::icheck "icheck" to \a number1 [1].
!! Similar to \ref opt-c "-c" or \ref opt-C "-C".
!! For mpmod::icheck "icheck" >0 no solution is performed but input statistics is checked in detail.
!! With mpmod::icheck "icheck" >1 the appearance range (first/last file,record and number of files)
!! of global parameters is determined too.
!! \subsection cmd-checkpargroups checkparametergroups
!! Set flag \ref mpmod::ichkpg "ichkpg" to 1 (true) to enable the checking of
!! (the rank of) \ref ch-pargroup "parameter groups".
!! \subsection cmd-chisqcut chisqcut
!! For local fit \ref an-chisq "setChi^2" cut \ref mpmod::chicut "chicut" to \a number1 [1.],
!! \ref mpmod::chirem "chirem" to \a number2 [1.].
!! \subsection cmd-comment comment
!! Define \ref ch-parcom "comments" for global parameters.
!! \subsection cmd-compress compress
!! Obsolete. Compression is default.
!! \subsection cmd-closeandreopen closeandreopen
!! Set flag \ref mpmod::keepopen "keepOpen" to zero to enable closing and reopening of binary files
!! to limit the number of concurrently open files.
!! \subsection cmd-constraint constraint
!! Define \ref sssec_consinf "constraints" for global parameters.
!! \subsection cmd-countrecords countrecords
!! Set flag \ref mpmod::mcount "mcount" to 1 (true) to enable parameter counting om record level.
!! \subsection cmd-debug debug
!! Set number of records with debug printout \ref mpmod::mdebug "mdebug" to
!! \a number1 [3], number of measurements with printout \ref mpmod::mdebg2 "mdebg2" to \a number2.
!! \subsection cmd-debug-pd debugPARDISO
!! Set \ref ch-pardiso "PARDISO" message level \ref mpmod::ipddbg "ipddbg" to
!! \a number1 [0]. Dump steering array \c IPARM(64) for value > 0.
!! \subsection cmd-dwfractioncut dwfractioncut
!! Set \ref an-dwcut "down-weighting fraction" cut \ref mpmod::dwcut "dwcut"
!! to \a number1 (max. 0.5).
!! \subsection cmd-outlierfracwarnthreshold outlierfracwarnthreshold
!! Set \ref an-outlierfracwarnthreshold "warning level for fraction of large chi2 entries" \ref mpmod::warnthresholdchi2 "warnthresholdchi2"
!! to \a number1 [1]. Value is expected to be in percent (1.0 parsed as 1%).
!! \subsection cmd-entries entries
!! Set \ref an-entries "entries" cuts for variable global parameter
!! \ref mpmod::mreqenf "mreqenf" to \a number1 [25],
!! \ref mpmod::mreqena "mreqena" to \a number2 [10] and
!! \ref mpmod::iteren "iteren" to the product of \a number1 and \a number3 [0].
!! \subsection cmd-errlabels errlabels
!! Define (up to 100 in total) global labels \a number1 .. \a numberN
!! for which the parameter errors are calculated for method MINRES too
!! (by \ref solglo "solving" \f$\Vek{C}\cdot\Vek{x}_i = \Vek{b}^i, b^i_j = \delta_{ij} \f$).
!! \subsection cmd-force force
!! Set force (iterations) flag \ref mpmod::iforce "iforce" to 1 (true).
!! Same as \ref opt-f "-f".
!! \subsection cmd-fortranfiles fortranfiles
!! Following binaries are Fortran files.
!! \subsection cmd-globalcorr globalcorr
!! Set flag \ref mpmod::igcorr "igcorr" for output of global correlations to 1 (true).
!! \subsection cmd-histprint histprint
!! Set flag \ref mpmod::nhistp "nhistp" for \ref an-histpr "histogram printout"
!! to 1 (true).
!! \subsection cmd-hugecut hugecut
!! For local fit set Chi^2 cut \ref mpmod::chhuge "chhuge"
!! for \ref sssec-outlierdeb "unreasonable data" to \a number1 [1.].
!! \subsection cmd-iterateentries iterateentries
!! Set maximum value \ref mpmod::iteren "iteren" for iteration of entries cut to
!! \a number1 [maxint]. Can alternatively be set by the \ref cmd-entries command.
!! For parameters with less entries the cut will be iterated ignoring measurements with
!! at least one parameter below \ref mpmod::mreqenf "mreqenf".
!! \subsection cmd-lapackerr lapackwitherrors
!! Set flag \ref mpmod::ilperr "ilperr" for calculation of inverse matrix
!! for parameter errors by \ref ch-lapack "LAPACK" to 1 (true).
!! \subsection cmd-linesearch linesearch
!! The mode \ref mpmod::lsearch "lsearch" of the \ref par-linesearch "line search"
!! to improve the solution is set to \a number1.
!! \subsection cmd-localfit localfit
!! For local fit set number of iterations \ref mpmod::lfitnp "lfitnp"
!! with calculation of pulls to \a number1, flag \ref mpmod::lfitbb "lfitbb"
!! for auto-detection of bordered band matrices to \a number2.
!! \subsection cmd-matiter matiter
!! Set number of iterations \ref mpmod::matrit "matrit" with (re)calcuation of
!! global matrix to \a number1.
!! \subsection cmd-matmoni matmoni
!! Set record interval \ref mpmod::matmon "matmon" for monitoring of (sparse) matrix
!! construction to \a number1.
!! \subsection cmd-maxlocalcond maxlocalcond
!! Set maximal Log10(condition) of decomposition of band part for local fit
!! \ref mpmod::cndlmx "cndlmx" to \a number1. Records with larger condition will be rejected.
!! \subsection cmd-maxrecord maxrecord
!! Set record limit \ref mpmod::mxrec "mxrec" to \a number1.
!! \subsection cmd-measurement measurement
!! Define (additional) \ref sssec_gpm "measurements" for global parameters.
!! \subsection cmd-memorydebug memorydebug
!! Set debug flag \ref mpmod::memdbg "memdbg" for memory management
!! to \a number1 [1].
!! \subsection cmd-method method
!! Has special format:
!!
!!         method   name     number1  number2
!!
!! Set \ref ch-methods "solution method" \ref mpmod::metsol "metsol" and
!! storage mode \ref mpmod::matsto "matsto" according to \a name,
!! (\c inversion : (1,1), \c diagonalization : (2,1),
!! \c decomposition : (3,1),
!! \c fullMINRES : (4,1) or \c sparseMINRES : (4,2),
!! \c fullMINRES-QLP : (5,1) or \c sparseMINRES-QLP : (5,2),
!! \c fullLAPACK factorization : (7,1), \c unpackedLAPACK factorization : (8,0)),
!! \c sparsePARDISO factorization : (9,3)
!! (minimum) number of iterations \ref mpmod::mitera "mitera" to \a number1,
!! convergence limit \ref mpmod::dflim "dflim" to \a number2.
!!
!! \c Inversion and \c diagonalization provide in addition to the solution the parameter errors
!! (from the diagonal of the inverted global matrix). Solutions with \c MINRES are only approximate.
!! \subsection cmd-monres monitorresiduals
!! Set flag \ref mpmod::imonit "imonit" for monitoring of residuals to \a number1 [3]
!! and increase number of bins (of size 0.1) for internal storage to \a number2 [100].
!! Monitoring mode \ref mpmod::imonmd "imonmd" is 0.
!! \subsection cmd-monpull monitorpulls
!! Set flag \ref mpmod::imonit "imonit" for monitoring of pulls to \a number1 [3]
!! and increase number of bins (of size 0.1) for internal storage to \a number2 [100].
!! Monitoring mode \ref mpmod::imonmd "imonmd" is 1.
!! \subsection cmd-monpgs monitorprogress
!! For progress monitoring set for repetition rate \c nrep the start value \ref mpmod::monpg1 "monpg1"
!! to \a number1 [1] and maximum increase \ref mpmod::monpg2 "monpg2" to \a number2 [1024].
!! Monitored are operations (inversion, decomposition, similarity) on the global and the constraints matrices.
!! If the (outermost loop) index is greater equal \c nrep the index is printed and \c nrep updated
!! (+ min(\c nrep, \c monpg2)).
!! \subsection cmd-mresmode mresmode
!! Set \ref minresqlpmodule::minresqlp "MINRES-QLP" factorization mode
!!  \ref mpmod::mrmode "mrmode" to \a number1.
!! \subsection cmd-mrestranscond mrestranscond
!! Set \ref minresqlpmodule::minresqlp "MINRES-QLP" transition (matrix) condition
!!  \ref mpmod::mrtcnd "mrtcnd" to \a number1.
!! \subsection cmd-mrestol mrestol
!! Set tolerance criterion \ref mpmod::mrestl "mrestl" for \ref minresmodule::minres "MINRES"
!! to \a number1 (\f$10^{-10}\f$ .. \f$10^{-4}\f$).
!! \subsection cmd-nofeasiblestart nofeasiblestart
!! Set flag \ref mpmod::nofeas "nofeas" for \ref an-nofeas "skipping"
!! making parameters feasible to \a number1 [1].
!! \subsection cmd-outlierdownweighting outlierdownweighting
!! For local fit set number of \ref sssec-outlow "outlier"
!! \ref an-downw "down-weighting" iterations
!! \ref mpmod::lhuber "lhuber" to \a number1.
!! \subsection cmd-pairentries pairentries
!! Set entries cut for variable global parameter pairs \ref mpmod::mreqpe "mreqpe"
!! to \a number1, histogram upper bound \ref mpmod::mhispe "mhispe" for pairs
!! to \a number2 (<1: no histogramming), upper bound \ref mpmod::msngpe "msngpe"
!! for pair entries with single precision storage
!! to \a number3.
!! \subsection cmd-parameter parameter
!! Define \ref sssec-parinf "initial value, pre-sigma" for global parameters.
!! \subsection cmd-pardiso pardiso
!! Modify for \ref ch-pardiso "PARDISO" the internal steering parameters.
!! \subsection cmd-postprocessing postprocessing
!! Define post processing *string*. Will be executed by system at end of **pede**.
!! \subsection cmd-presigma presigma
!! Set default pre-sigma \ref mpmod::regpre "regpre" to \a number1 [1].
!! \subsection cmd-print print
!! Set print level \ref mpmod::mprint "mprint" to \a number1 [1].
!! \subsection cmd-printcounts printcounts
!! Set flag \ref mpmod::ipcntr "ipcntr" to \a number1 [1].
!! The counters for the global parameters from the accepted local fits (=1)
!! or from the binary files (>1) will be printed in the result file. Alternatively
!! the counters for zero global derivatives from the binary files (<0) can be selected.
!! \subsection cmd-printrecord printrecord
!! \ref an-recpri "Record" numbers with printout.
!! \subsection cmd-pullrange pullrange
!! Set (symmetric) range \ref mpmod::prange "prange" for histograms
!! of pulls, normalized residuals to \a number1 (=0: auto-ranging).
!! \subsection cmd-readerroraseof readerroraseof
!! Set flag \ref mpmod::ireeof "ireeof" to 1 (true) to treat read errors for binary files
!! as end-of-file instead of aborting.
!! \subsection cmd-regularisation regularisation
!! Set flag \ref mpmod::nregul "nregul" for regularization to 1 (true),
!! regularization parameter \ref mpmod::regula "regula" to \a number2,
!! default pre-sigma \ref mpmod::regpre "regpre" to \a number3.
!! \subsection cmd-regularization regularization
!! Set flag \ref mpmod::nregul "nregul" for regularization to 1 (true),
!! regularization parameter \ref mpmod::regula "regula" to \a number2,
!! default pre-sigma \ref mpmod::regpre "regpre" to \a number3.
!! \subsection cmd-resolveredundancycons resolveredundancycons
!! Set flag \ref mpmod::irslvrc "irslvrc" to 1 (true).
!! Redundancy constraints will be resolved 
!! (parameters appearing in constraints will be fixed, constraints skipped).
!! \subsection cmd-scaleerrors scaleerrors
!! Set measurement scaling factors \ref mpmod::dscerr "dscerr"
!! to \a number1 [1.] and \a number2 [\a number1].
!! First value is for "global" measurements (with global derivatives),
!! second for "local" measurements (without global derivatives).
!! \subsection cmd-skipemptycons skipemptycons
!! Set flag \ref mpmod::iskpec "iskpec" to 1 (true).
!! Empty constraints (without variable parameters) will be skipped.
!! \subsection cmd-subito subito
!! Set subito (no iterations) flag \ref mpmod::isubit "isubit" to 1 (true).
!! Same as \ref opt-s "-s".
!! \subsection cmd-threads threads
!! Set number \ref mpmod::mthrd "mthrd" of OpenMP&tm; threads for processing
!! to \a number1,
!! number \ref mpmod::mthrdr "mthrdr" of threads for reading
!! binary files to \a number2 [\a number1].
!! \subsection cmd-weightedcons weightedcons
!! Set flag \ref mpmod::iwcons "iwcons" to \a number1 [1].
!! Implements \ref sssec_consinf "weighted constraints" for global parameters.
!! \subsection cmd-withelim withelimination
!! Set flag \ref mpmod::icelim "icelim" to 1 (true).
!! Selects solution by elimination for linear equality constraints.
!! \subsection cmd-withlapackelim withlapackelimination
!! Set flag \ref mpmod::icelim "icelim" to 2 (LAPACK).
!! Selects solution by elimination for linear equality constraints with LAPACK.
!! Only available for unpacked LAPACK!
!! \subsection cmd-withmult withmultipliers
!! Set flag \ref mpmod::icelim "icelim" to 0 (false).
!! Selects solution by Lagrange multipliers for linear equality constraints.
!! \subsection cmd-wolfe wolfe
!! For strong Wolfe condition in \ref par-linesearch "line search"
!! set parameter \ref mpmod::wolfc1 "wolfc1" to \a number1, \ref mpmod::wolfc2
!! "wolfc2" to \a number2.

!> \page exit_code_page List of exit codes
!! The exit code and message of the **Pede** executable can be found in the
!! file <tt>millepede.end</tt> :
!!    + <b>-1</b>   Still running or crashed
!!    + **00**   Ended normally
!!    + **01**   Ended with warnings (bad measurements)
!!    + **02**   Ended with severe warnings (insufficient measurements)
!!    + **03**   Ended with severe warnings (bad global matrix)
!!    + **04**   Ended with severe warnings (bad binary file(s))
!!    + **05**   Ended without solution (empty constraints)
!!    + **10**   Aborted, no steering file
!!    + **11**   Aborted, open error for steering file
!!    + **12**   Aborted, second text file in command line
!!    + **13**   Aborted, unknown keywords in steering file
!!    + **14**   Aborted, no binary files
!!    + **15**   Aborted, open error(s) for binary files
!!    + **16**   Aborted, open error(s) for text files
!!    + **17**   Aborted, file name too long
!!    + **18**   Aborted, read error(s) for binary files
!!    + **19**   Aborted, binary file(s) modified
!!    + **20**   Aborted, bad binary records
!!    + **21**   Aborted, no labels/parameters defined
!!    + **22**   Aborted, no variable global parameters
!!    + **23**   Aborted, bad matrix index
!!    + **24**   Aborted, vector/matrix size mismatch
!!    + **25**   Aborted, result vector contains NaNs
!!    + **26**   Aborted, too many rejects
!!    + **27**   Aborted, singular QL decomposition of constraints matrix
!!    + **28**   Aborted, no local parameters
!!    + **29**   Aborted, factorization of global matrix failed
!!    + **30**   Aborted, memory allocation failed
!!    + **31**   Aborted, memory deallocation failed
!!    + **32**   Aborted, iteration limit reached in diagonalization
!!    + **33**   Aborted, stack overflow in quicksort
!!    + **34**   Aborted, pattern string too long - obsolete
!!    + **35**   Aborted, mismatch of number of global parameters
!!    + **36**   Aborted, unselected solution method requested
!!    + **40**   Aborted, other errors

!> \page troubleshooting_page Troubleshooting
!!
!! \section ch-experience General experience:
!!
!! In the ideal case the track and geometry models are correct and complete,
!! the measurement errors are gaussian and described perfectly
!! and all dependencies are linear. With this (probably unrealistic)
!! requirements pede should end *normally* (exit code 0 
!! (in <tt>millepede.end</tt>)).
!!
!! In the real world measurement errors can be non gaussian (e.g.
!! binary readout or multiple scattering (tails)) or dependencies
!! non linear (e.g. curvature or rotations). Therefore pede usually
!! ends with *warnings* (exit code 1).
!!
!! In case of problems with the linear equation system (e.g. too few measurements,
!! linear dependent constraints or global matrix not positive definite) pede
!! still tries to get a solution and ends with *severe warnings*
!! (exit code 2). The solution is very probable not reliable and should be
!! 'handled with care'. The pede input (steering and binary files) should be
!! \ref ch-checkinput "checked" thoroughly.
!!
!! \section ch-aborts Aborts
!! In case of severe (usually technical) problems pede aborts (exit code >= 10).
!! It may be helpful to \ref ch-checkinput "check" the pede input
!! (steering and binary files).
!!
!! For 'too many rejects' (exit code 26) usually the measurement errors
!! are largely underestimated. It may help to inflate those in pede
!! (command \ref cmd-scaleerrors) or to skip the internal iterations
!! (command \ref cmd-subito or command line option \ref opt-s)
!! to get *some* solution.
!!
!! \section ch-segfault Segmentation fault
!! In case of segmentation faults the OpenMP stack size may be to small (default 8M?).
!! Increased by defining the shell variable <tt>OMP_STACKSIZE</tt>:
!!
!!        export OMP_STACKSIZE=32M
!!
!! For further debugging one should add the compiler options
!! \"<tt>-g -fcheck=all -fbacktrace</tt>\"
!! to the \c F_FLAGS in the \c Makefile and recompile.
