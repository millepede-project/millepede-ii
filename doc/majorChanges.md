\page changes_page Major changes

Major changes with respect to the \ref draftman_page "draft manual".
\tableofcontents

# External Mille library
  Starting with version 05-00-00, the built-in components for reading `Mille` binaries have been replaced by an external package.
  This can read a range of formats,
  1. Uncompressed C-binaries - expected file extension .dat
  2. Compressed C-binaries - expected file extension .dat.gz
  3. ROOT binaries - expected file extension .root 
  4. Plain text binaries - expected file extension .csv

  The internal file format used in the external package is fully compatible with the existing one, all existing binaries should still work. 
  


# Extended entries cut
 The [entries](option_page.html#entries) cut has now three arguments:
 1. \ref mpmod::mreqenf "mreqenf": Minimum required (number of) \ref an-entries "entries"
    in binary files for global parameters to be *variable*.
 2. \ref mpmod::mreqena "mreqena": Minimum required (number of) entries
    (with record/track) \ref localfit-rejection "accepted" after the local (track) fit for
    variable global parameters. Due to the tigthening of the [chisqcut](option_page.html#chisqcut) during
    the internal pede iterations the number of *accepted* entries can get far lower than the
    number of entries in the binary files(&ge; **mreqenf**). If the number of entries gets lower
    than the number of variable global parameters determined from the related measurements
    there will be no more proper solution (for this part of the global linear equation
    system). Therefore *severe* warnings are issued
    in case the number of accepted entries is below **mreqena** for any variable
    global parameter in any iteration.
    As the gobal matrix is usually constructed only in the first iteration it probably will
    not get singular in this case, but the solutions for the affected parameters are not
    independent anymore.
 3. **iarg3** (\ref mpmod::iteren "iteren"=mreqenf*iarg3): Iteration of entries cut.
    It can be problematic if different global parameters related to the same measurements
    get different (number of) entries. If some are below and some are above **mreqenf**
    some are fixed and some are variable which may not make sense. An example is the
    rigid body alignment of a planar silicon pixel sensor in the local (uvw) system. For the
    two indendent measurements (in u and v) the corresponding alignment offsets (du, dv)
    appear each only for *one* measurement. The offset perpendicular to the sensor plane (dw)
    and the rotations affect *both* measurement. Therefore the offsets (du, dv) have a factor
    two smaller number of entries as the other four. To handle this situation the entries cut
    (on the binary files) can be iterated: If in a measurement some global parameters are fixed
    (due to entries < **mreqenf**) for the other parameters the cut value is increased by
    the constant factor **iarg3**.
    For the axample above a factor two would be appropriate (all six parameters fixed now).

    Another option to avoid this problem is to not count the entries on the level of equations
    (measurements) but on the level of records (tracks): [countrecords](option_page.html#countrecords).

    A third option is to produce the binary files without zero supression (of global derivatives).

 # Solution methods
 The following methods to obtain the solution \f$\vec{x}\f$ from a
 linear equation system \f$\mathbf{A}\cdot\vec{x}=\vec{b}\f$ are implemented:
 ## Inversion
 The solution and the covariance matrix \f$\mathbf{A}^{-1}\f$ are obtained by
 \ref an-inv "inversion" of \f$\mathbf{A}\f$.
 Available are the value, error and global correlation for all global parameters.
 The matrix inversion \ref sqminl "routine" has been [parallelized](#parallelization) 
 and can be used for up to several 10000 parameters.
 ## Diagonalization
 The solution and the covariance matrix \f$\mathbf{A}^{-1}\f$ are obtained by
 \ref an-diag "diagonalization" of \f$\mathbf{A}\f$.
 Available are the value, error, global correlation and
 eigenvalue (and eigenvector) for all global parameters.
 ## Minimal Residual Method (MINRES)
 The solution is obtained by minimizing \f$\Vert\mathbf{A}\cdot\vec{x}-\vec{b}\Vert_2\f$
 iteratively and is only approximate. \ref minresmodule::minres "MINRES"  [ref 8](index.html#references) is a special case of the
 generalized minimal residual method (\ref an-gmres "GMRES") for symmetric matrices.
 Preconditioning with a band matrix of zero or finite
 \ref mpmod::mbandw "bandwidth" is possible.
 Individual columns \f$\vec{c_i}\f$ of the covariance matrix can be calculated by
 solving \f$\mathbf{A}\cdot\vec{c}_i=\vec{1}_i\f$ where \f$\vec{1}_i\f$ is the i-th
 column on the unit matrix.
 The most time consuming part (\ref avprod "product" matrix times vector per iteration)
 has been [parallelized](#parallelization).
 Available are the value for all (and optionally error, global correlation
 for few) global parameters.
 ## Advanced Minimal Residual Method (MINRES-QLP)
 The \ref minresqlpmodule::minresqlp "MINRES-QLP" implementation [ref 9](index.html#references)
 is a MINRES evolution with improved norm estimates and stopping conditions
 (leading potentially to different numbers of internal iterations).
 Internally it uses QLP instead of the QR factorization in
 MINRES which should be numerically superior and allows to find for
 singular systems the minimal length (pseudo-inverse) solution.

 The default behavior is to start (the internal iterations) with QR factorization
 and to switch to QLP if the (estimated) matrix condition exceeds
 [mrestranscond](option_page.html#mrestranscond) ("mrtcnd" in sources). Pure QR or QLP factorization can be enforced
 by [mresmode](option_page.html#mresmode) ("mrmode" in sources).

 ## Elimination of constraints
 As alternative to the Lagrange multiplier method the solution by elimination
 has been added for problems with linear equality constraints.
 A \ref mpqldec::qldec "QL factorization" (with Householder reflections) of the
 transposed constraints matrix is used to transform to an unconstrained problem.
 For sparse matrix storage the sparsity of the global matrix is preserved.

 ## Decomposition
 The solution is obtained by a root-free Cholesky decomposition (LDLt).
 The covarinance matrix is *not* being calclulated.
 The method ia about a factor 2-3 faster than inversion (according to several tests).
 It is restricted to solution by elimination for problems with linear equality constraints
 (requires positive definite matrix). Supports block matrix storage for disjoint parameter blocks.

 ## LAPACK
 For these optional methods the solution is obtained by matrix factorization from an external LAPACK library.
 There exist (open or proprietary) implementations heavily optimized for specific hardware
 (and partially multi-threaded)
 which could easily be an order of magnitude faster than e.g. the custom code used for
 [decomposition](#decomposition).
 Tested has been the [Intel MKL](https://software.intel.com/en-us/intel-mkl) and
 [OpenBLAS](https://www.openblas.net).

 The symmetric global matrix can be stored in **packed** triangular (similar to **full** for inversion etc)
 or **unpacked** quadaratic \ref mpmod::matsto "shape".
 For unpacked storage the ["elimination of constraints"](#elimination-of-constraints) is
 implemented with LAPACK (DGEQLF, DORMQL) as alternative.

 For positive definite matrices (elimination of constraints) a Cholesky factorization is used (DPPTRF/DPOTRF) and
 for indefinite matrices (Lagrange multipliers for constraints) a Bunch-Kaufman factorization (DSPTRF/DSYTRF).
 With the option [lapackwitherrors](option_page.html#lapackwitherrors) the inverse matrix is calculated from the  factorization to provide parameter errors (DPPTRI/DPOTRI, DSPTRI/DSYTRI, potentially slow, single-threaded).

 ## PARDISO
 This optional method is using the [Intel oneMKL PARDISO](
 https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-2/onemkl-pardiso-parallel-direct-sparse-solver-iface.html)
 shared-memory multiprocessing parallel direct sparse solver. This aims at use cases with really sparse global
 matrices (e.g. less than 10% nonzero elements).

 The upper triangle of the (symmetric) global matrix is stored in the 3-array variation of the compressed sparse row format (CSR3).
 It can be checked for a (fixed size) \b block structure to reduce the storage size (BSR3).
 (This could be further improved by reordering - not implemeted.)

 For problems with linear equality constraints the usage of Lagrange multipliers is enforced.

 For the internal steering parameters \c IPARM(64) the default values are used.
 They can be modified with the [pardiso](option_page.html#pardiso) keyword
 (specifying pairs of indices (1-64) and values):

          pardiso
          ...      ...
          index    value
          ...      ...

 Related options:
   - \ref [blocksizePARDISO](option_page.html#blocksizepardiso), e.g. combinations of factors 2 and 3:

           blocksizePARDISO 2 3 4 6 9 12

   - [debug](option_page.html#debugpardiso).

 # Regularization
 Optionally a term \f$\tau\cdot\Vert\vec{x}\Vert\f$ can be added to the objective function
 (to be minimized) where \f$\vec{x}\f$ is the vector of global parameters
 weighted with the inverse of their individual pre-sigma values.

 # Local fit
 In case the \ref par-locfitv "local fit" is a track fit with proper description of multiple
 scattering in the detector material additional local parameters have to be introduced
 for each scatterer and solution by *inversion* can get time consuming
 (~ \f$n_{lp}^3\f$ for \f$n_{lp}\f$ local parameters). For trajectories based on
 **broken lines** [refs 4,6](index.html#references) the corresponding matrix \f$\mathbf{\Gamma}\f$
 has a bordered band structure (\f$\Gamma_{ij}=0\f$ for \f$\min(i,j)>b\f$
 (border size) and \f$|i-j|>m\f$ (bandwidth)). With
 <i>root-free Cholesky decomposition</i> the time for the solution is linear
 and for the calculation of \f$\Gamma^{-1}\f$
 (needed for the construction of the global matrix) quadratic in \f$n_{lp}\f$.
 The condition of the diagonal matrix from the decomposition (of the band part)
 can be used to reject ill conditioned cases (see [maxlocalcond](option_page.html#maxlocalcond)).
 For each local fit the structure of \f$\mathbf{\Gamma}\f$ is checked and the faster
 solution method selected automatically.

 # Parallelization
 The code has been largely parallelized using [OpenMP&tm;](http://www.openmp.org).
 This includes the reading of binary files, the local fits, the construction of the
 sparsity structure and filling of the global matrix and the global fit
 (except by diagonalization). The number of threads is set by the command
 [threads](option_page.html#threads).
 The OpenMP environment can be displayed with:

       export OMP_DISPLAY_ENV=TRUE

 \b Caching. The records are read in blocks into a *read cache* and processed from
 there in parallel, each record by a single thread. For the filling of the global
 matrix the (zero-compressed) update matrices (\f$\mathbf{\Delta C}_1+\mathbf{\Delta C}_2\f$ from
 equations \ref eq-c1 "(15)", \ref eq-c2 "(16)")
 produced by each local fit are collected in a
 *write cache*. After processing the block of records this is used to update
 the global matrix in parallel, each row by a single thread.
 The total cache size can be changed by the command [cache](#option_page.html#cache).



 # Sparse matrix storage

 ## Compression
 In sparse storage mode for each row the list of column indices (and values) for the
 non-zero elements are stored. With compression regions of continous column indices
 are represented by the first index and their number (packed into a single 32bit
 integer). Compression is selected by the command [compress](option_page.html#compress).
 In addition rare elements can be neglected (,histogrammed) or stored in single instead
 of double precision according to the [pairentries](option_page.html#pairentries) command.
 ## Parameter groups
 With the implementation of parameter groups (sets of adjacent global parameters (labels)
 appearing in the binary files *always* together) the sparsity structure addresses
 not single values anymore but block matrices with sizes according to the contributing
 parameter groups. Groups with adjacent column ranges are combined into a single block matrix.
 The offset and first column index in a (block) row is stored.

 # Gzipped C binary files
 The [zlib](zlib.net) can be used to directly read *gzipped* C binary files.
 In this case reading with multiple threads
 (each file by single thread) can speed up the decompression.

 # Weighted binary files
 A constant weight can be applied to all records in a binary file:

         <file name> -- <weight> ! <comment>

 # Transformation from FORTRAN77 to Fortran90
 The **Millepede** source code has been formally transformed from <i>fixed form</i>
 FORTRAN77 to <i>free form</i> Fortran90 (using TO_F90 by Alan Miller)
 and (most parts) modernized:
 - <tt>IMPLICIT NONE</tt> everywhere. Unused variables removed.
 - \c COMMON blocks replaced by \c MODULEs.
 - Backward \c GOTOs replaced by proper \c DO loops.
 - \c INTENT (input/output) of arguments described.
 - Code documented with doxygen.

 Unused parts of the code (like the interactive mode) have been removed.
 The reference compiler for the Fortran90 version is gcc-4.6.2 (gcc-4.4.4 works too).

 # Memory management
 The memory management for dynamic data structures (matrices, vectors, ..)
 has been changed from a \ref an-dynal "subdivided" *static* \c COMMON block to
 *dynamic* (\c ALLOCATABLE) Fortran90 arrays. One **Pede** executable is now
 sufficient for all application sizes.

 # Read buffer size
 In the \ref sssec-loop1 "first loop" over all binary files a preset
 \ref mpmod::ndimbuf "read buffer size" is used. Too large records are skipped,
 but the maximal record length is still being updated. If any records had to be skipped
 the read buffer size is afterwards adjusted according to the maximal record length
 and the first loop is repeated.

 # Number of binary files
 The number of binary files has no hard-coded limit anymore, but is calculated from
 the steering file and resources (file names, descriptors, ..)
 are allocated dynamically. Some resources may be limited by the system.

 # Check input mode
 The **Pede** executable can be run in special modes to only check the input
 without producing a solution. This can be useful in case of problems.
 Checked are the global parameters and their (equality) constraints from the steering
 and binary files. Debug information is written to the standard output and the
 <tt>millepede.res</tt> text file.

 1. Basic checks: Command '[checkinput](option_page.html#checkinput) 1' or command line option [-c](option_page.html#-c)

   To standard output:
   - For each constraint the number of all and of variable global parameters involved.
     *Empty* constraints without variable parameters should be removed.
   - Constraint group definitions (first and last (sorted (by global labels)) constraint,
     first and last global label). Are disjoint (in global labels) sets of constraints.
   - Constraint block definitions (first and last constraint group, first and last global label).
     Are non overlapping (in global labels) sets of constraints.

   To text file:
   - For each global parameter the label, value, presigma, number of entries, constraint
     group (0 for none) and status (variable or fixed).
   - For the first parameter of a parameter group (consecutive parameter (labels) appearing
     always togther) the line starts with '<tt></tt>'.
   - For each constraint group the number of contributing constraints, number of entries
     and first and last global label.

 2. More checks: Command '[checkinput](option_page.html#checkinput) 2' or command line option [-C](option_page.html#-c)

   Additionally to standard output:
   - For each **sorted** constraint the index by appearance (in text file) with number
     of first line and first and last global label.
   - For each constraint group the number of constraints and the rank of the corresponding
     block of the product matrix. Must be equal!

   Additionally to text file:
   - For each global parameter the appearance (in binary files) statistics
     (first file and record number, last file and record number and number of files)
     and number of paired global parameters (appearing together).
   - For each constraint group the label ranges of (not contributing) paired global parameters.
   - For each constraint group the appearance statistics (using the contributing global labels).

 # Global parameter comments
 For each global parameter a \ref mpdef::itemclen "brief" **comment** can be defined to annotate the <tt>millepede.res</tt> text file.
 The syntax is similar as for the "Parameter" or "Constraint" keywords:

         Comment
         ...     ...
         label   text
         ...     ...
