\page release_page List of stable releases
 
## Version history

### Version V05-00-00 (new)

Retire Makefile build, adopt external Mille package

Clone Id replaced by direct version printout starting with 05-00-00.

### Version V04-19-03 

Fix an FPE due to incomplete array initialisation - no impact on result, only relevant for debug builds.

(Clone Id: c5abe0f4ddabe767dedc2f9c03d53cd6607db438)

### Version V04-19-02

Bugfix to cToRoot program for large input binaries. 

(Clone Id: 6b6285a829566d03769f4dd301c5422ddb23d86e)

### Version V04-19-01

Changes to CMake configuration - 'millepedeII' is now discoverable by find_package. Use location-indepedent ROOT linking 

(Clone Id: 6b6285a829566d03769f4dd301c5422ddb23d86e)

### Version V04-19-00 

Option to store and read binaries in ROOT format, C-to-root binary conversion, unit tests. 

(Clone Id: 6b6285a829566d03769f4dd301c5422ddb23d86e)

### Version V04-18-00

New CMake-based build, revised documentation, option to set number of records in internal test. 

(Clone Id: 685557b74073e7efc2b280d883640d7f58fc4ebe)

### Version V04-17-07

Configurable warning threshold for outlier fraction (`outlierfracwarnthreshold`)

(Clone Id: b15d9f1893ea8ea05892c71644543c861e35bd00)

### Version V04-17-06

Some optimisations for `julia/tinypede.jl` tool.

(Clone Id: 8eeee5b5d97c2f73818cac1d7c8a50bd0befa14f)

### Version V04-17-05

Some tools have been reimplemented in [Julia](https://julialang.org) (`julia/readMilleBinary.jl`, `julia/tinypede.jl`).

(Clone Id: a5226a0c0002cfcb2639c6a9b17286b8d1b9f8e3)

### Version V04-17-04

Proper abort (message) in case of wrong binary file type (Fortranfiles or Cfiles) added. A wrong file type typically leads to bad records with read errors or invalid record size.

(Clone Id: d281f2231ebd7fabc50925038cb0bd05e9f5bff2)

### Version V04-17-03

Implement post processing at end of **pede**. The post processing string is executed by the (operating) system (e.g. option "`postprocessing pwd; ls -l millepede.??? ! simple check`").

(Clone Id: dd0c569a1aafb6f6eb2f26b9b9537a685639ca25)

### Version V04-17-01

For instrumentation and profiling with [Score-P](http:score-p.org) switched from POMP (`scorep --pomp`) (based on OPARI2, to be superseded by OMPT (OpenMP 5.0)) to user (`scorep --user`) regions.

(Clone Id: 2abb4a313debe0a8077a24e1c46a3c0875f9041b)

### Version V04-17-00

Modernisation of development environment (from EL7 to EL9 (gcc11.4), from ompP to [Score-P](http:score-p.org) for profiling).

(Clone Id: 75827788c0ff437e6144a75b1027a9baea8a0ce5)

### Version V04-16-05

For local fits with bordered-band matrix structure (e.g. using (general) broken lines) now use the condition of the diagonal matrix from the root-free Cholesky decomposition (of the band part) to optionally reject 'ill conditioned' records. With the new option `maxlocalcond` an upper limit on Log10(condition) can be set.

(Clone Id: cd0c67b54a617e1263fab81c9ad3c282df8a9e80)

### Version V04-16-04

In first loop over binary files for counting and grouping of global labels now ignore those with zero global derivative. (Previously zero-suppressed (global derivatives in the) binary files had been assumed.)

(Clone Id: c5a7342b3793f36f30ad6658bc1f72bc74250cf7)

### Version V04-16-03

Added `legacy` folder (Millepede-I code and documentation).

(Clone Id: cfed4d61561fe1111c04b75d70c482d786ba4fff)

### Version V04-16-02

Added `tinypede.py` to tools: **Pede** implementation in python3 with basic functionality for illustration or testing with small problems.

(Clone Id: 33bbf5bd3de60676558a758745bccee1006d0a75)

### Version V04-16-01

Counters scaling with the number of records are now long (64bit) integers.

(Clone Id: 5056d5ae06dcee7bffa50aeeb65cbf3c3645125b)

### Version V04-16-00

Summation of global Chi2 and NDF revised (no more 32bit variables, update per record). Before in rare cases integer overflows lead to negative Chi2 values.

(Clone Id: c6be33cb87ef25e4226d5088ec53ebf6d39f4712)

### Version V04-15-00

Fix severe problem with external measurements depending on multiple global parameters (was broken before).

(Clone Id: c19e0780f000076328c33b39811fa9ce02890377)

### Version V04-14-00

New optional (solution) method `sparsePARDISO` using the [Intel oneMKL PARDISO](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-2/onemkl-pardiso-parallel-direct-sparse-solver-iface.html) shared-memory multiprocessing parallel direct sparse solver added. This aims at use cases with really sparse global matrices (e.g. less than 10% nonzero elements).

(Clone Id: a02374cf8daa2b6d49a3ad18b1b1537e00828738)

### Version V04-13-06

Fix line search parameters for LAPACK.

(Clone Id: 39f9b5d09202485a84b06607c7c0bf468a32869a)

### Version V04-13-05

Fix problem for block diagonal global matrix (keeping single block).

(Clone Id: 2124f90ac72e5e156faa22990e65f20002203779)

### Version V04-13-04

Fix problem with monitoring of residuals. Calculate _skyline_ fraction for sparse matrices.

(Clone Id: 686292e46bfdda2b2b6d721969a8d471462099b5)

### Version V04-13-03

Several new features:

- Check for _redundancy_ constraints: Constraints from a 1-to-1 (hierarchy) relation describing a linear transformation between two equivalent sets of global parameters. Optionally with `resolveredundancycons` resolve those: Fix contributing global parameters and skip constraint (to save resources and burden on numerics).
- Optionally check with `checkparametergroups` (the rank of) the global parameter groups to verify the linear independency of the global derivatives (in a group).
- Optionally add a brief `comment` to global parameters in the steering to annotate the results file.

(Clone Id: d36231056d77020a30ba0bd90e988374eaee657a)

### Version V04-13-02

Cleanup and fixes (constraint elimination with LAPACK using constraint groups, proper termination of file names for 'iso_c_binding').

(Clone Id: 9bc61c3279cd02bdab89d764c971d7e16dbae12e)

### Version V04-13-01

Tool `readMilleBinary.py` is now compatible with python3, updated CLI.

(Clone Id: 648815480555bcf98cd72903fc89ecefa9af9be2)

### Version V04-13-00

Fortran/C interoperability uses now 'iso_c_binding' (fortran 2003) instead of 'cfortran.h'.

(Clone Id: c102eccc9e1e85a68259a259882bc15e22a0ca59)

### Version V04-12-03

Fix global parameter errors for solution by diagonalization using elimination of constraints.

(Clone Id: 56e1ade063eecf11086597550d153cd1c9cd9000)

### Version V04-12-02

Force check input mode (2) in case of empty constraints (no variable parameters). No solution will be calculated.

(Clone Id: 610e2b012f1a4e0da4505318b3495b869e40780a)

### Version V04-12-01

Detailed example (-t=BRLF) included in doxygen documentation. Python2 script tools/compareResults.py added.

(Clone Id: cc523e423547ec1e335b86a2ead55e4b2e6b3040)

### Version V04-12-00

More code modernisation to comply the with fortran standard 2018. Still some GNU fortran extensions are used.

(Clone Id: 389408e8154aac37a7e22c031482e2f68e0e7620)

### Version V04-11-04

Fix (uninitialsed values) and cleanup for internal silicon strip tracker example.

(Clone Id: 72d776418dbf55914ee2c77faa40dd1403b41f70)

### Version V04-11-03

Cleanup, (rare) problem with construction of parameter groups fixed (avoiding aborts with exit code 35).

(Clone Id: 7dbd0f53285a84a16e56ec2315d424d04c88d681)

### Version V04-11-02

Cleanup, test programs to print LAPACK (library) configuration added to tools directory.

(Clone Id: f48bce3c294217cac74be41d6ecd4da10d15519e)

### Version V04-11-01

Constraint groups implemented, documentation (check input, trouble shooting) added.

(Clone Id: c17943c47abedc4cdae3634611ae3b3b2f65f3eb)

### Version V04-11-00

Migration to GitLab.

(`git clone --depth 1 --branch V04-11-00` Id: 15b0d4c1bfffa9656e136d720d05bbc8751e233b)

### Version V04-10-01

Fortran code modernized further (assumed-size array arguments replaced). Additional cases of out of range read access fixed.

### Version V04-10-00

The code quality has been improved. Compilation with '-fcheck=all’ has been used to check for out of index range array access (few cases for read access found and fixed) and ancient Fortran statements have been replaced. For one user an EQUIVALENCE statement in the internal histogramming produced (with gcc10) memory corruption with NANs leading to an indefinite loop.

### Version V04-09-05

Fixed uninitialized values in QL decomposition of constraints matrix (introduced in V04-09-04).

### Version V04-09-04

Exploit sparsity and decomposition of constraints matrix into disjoint blocks for all solution methods (e.g. QL decomposition, MINRES preconditioner).

### Version V04-09-03

Memory allocation for internal histograms is now dynamic in terms of number of binary files.

### Version V04-09-02

The warning messages for severe warnings (exit code 2) have been improved.

### Version V04-09-01

All indices for constraint matrices have been upgraded to 64 bit. Previously a segmentation fault happened for (number of parameters) \* (number of constraints) \>= 2^31.

### Version V04-09-00

Optional usage of external linear algebra libraries to speed up solution. Implemented has been LAPACK and tested with the [Intel MKL](https://software.intel.com/en-us/intel-mkl) and the [OpenBLAS](https://www.openblas.net) library.

Two new solution methods based on matrix factorisation are available optionally: `fullLAPACK` and `unpackedLAPACK`. They differ in the LAPACK storage mode (packed vs unpacked).

### Version V04-08-02

New command `monitorprogress` to monitor progress in operations on global and constraints matrices with additional printout into the log file.

### Version V04-08-01

Fixing potential integer overflows for elimination of constraints. This happens for all solution methods except MINRES if the square of the number of (variable global) parameters is larger than 2^31.

### Version V04-08-00

New solution method (Cholesky) `decomposition` added. It's faster than `inversion` but does not calculates parameter errors.

### Version V04-07-04

Bug fix (integer overflow) and little speedup (\~15%) for solution of large problems by inversion.

### Version V04-07-03

The solution by inversion with Lagrange multipliers has been fixed. (It was broken since V04-06-00.)

### Version V04-07-02

More robust sorting of constraints for splitting into disjoint blocks.

### Version V04-07-01

The counting of the appearance of global parameters in the binary files can now be done on record (e.g. track) level instead of equation (e.g. measurement) level. This is enabled with the new command `countrecords` and makes the iteration of the first data loop (by `iterateentries`) obsolete.

### Version V04-07-00

Implementation of parameter groups (sets of adjacent global parameters (labels) appearing in the binary files always together). Used to speed up construction of global matrix. Similarity operations are now aware of sparse (rectangular) matrices.

Solution with method `sparseMINRES` got considerably faster.

### Version V04-06-01

Modifications for compilation with PGI compiler (`make -f Makefile_pgi`).

### Version V04-06-00

Checking global parameters for disjoint blocks now. In case of solution by inversion (optionally with constraints handled by elimination) switch to block diagonal storage mode (to save memory and cpu time).

### Version V04-05-04

Method `sparseMINRES` fine again.

### Version V04-05-03 (bad)

Improved stability. The line search is now skipped in case the current solution is numerically too close to the minimum. Method `sparseMINRES` broken due to bad back port from development version.

### Version V04-05-02

Two (very old) cases of access to arrays outside of the defined bounds (reported by gcc9) have been fixed.

### Version V04-05-01

Fix for modified text parser (concerning leading blanks in lines with only numbers).

### Version V04-05-00

Update of (approximate) string matching for keyword detection. Matching is now symmetric in pattern and text. Previously e.g. a binary file with the letters from '`Cfiles`' in the name in that order was treated as that keyword and not as a binary file.

### Version V04-04-01

Cleanup of operations (open, close, rewind) on binary files. New command `closeandreopen` to enable closing and reopening of binary files to limit the number of concurrently open files. The modification dates of the files are monitored to ensure data integrity.

### Version V04-04-00

Constraints are now sorted and split into disjoint blocks to speed up calculation of rank and QL decomposition by block matrix algebra. This works best if the label sets of the involved alignable objects are disjoint too.

### Version V04-03-10

Some minor fixes, additional level of detail (appearance range of global parameters in binary files) for checkinput mode.

### Version V04-03-09

Fixes for special data blocks and not decreasing function value.

### Version V04-03-08

More debug output in case of problems reading the binary Cfiles.

### Version V04-03-07

Fix for downweighting correction (Cauchy) of global chi2.

### Version V04-03-06

Fix for abort on singular QL decomposition of constraints matrix.

### Version V04-03-05

Monitoring of pulls (`monitorpulls`) and scaling of measurement errors (`scaleerrors`) from binary files added. Abort now in case of singular QL decomposition of constraints matrix (solution by elimination).

### Version V04-03-04

Bug fix for monitoring of residuals.

### Version V04-03-03

Contains few fixes, mainly improved response to numerical problem leading before to "Function not decreasing".

### Version V04-03-02

The monitoring of the residuals per local fit cycle has been implemented and is selected by `monitorresiduals`. The normalized residuals are grouped by the first global label and the median and the RMS (from the median of the absolute deviations) per group are written to the new text file `millepede.mon`.

### Version V04-03-01

Contains fixes and improvements for MINRES preconditioning:

- Preconditioning with a band matrix works now for solution by elimination of constraints too.
- Preconditioning with a skyline matrix added to avoid the potential rank deficits of a band matrix (by reduction to the diagonal).

### Version V04-03-00

Contains some new features:

- Solution by elimination for problems with linear equality constraints has been implemented (as default, new command `withelimination`) in addition to the Lagrange multiplier method (new command `withmultipliers`).
- The skipping of _empty_ constraints (without any variable parameters) can be enabled with the new command `skipemptycons`.
- With new command `checkinput` detailed checks of input data (binary files, constraints) are performed, but no solution will be determined. Some input statistics is available in the output file `millepede.res`.
- The iteration of the entries cut has been implemented with the new command `iterateentries`. In the second iteration measurements with any parameters fixed by the previous entries cut are skipped. Useful if for some measurements different parameters have different number of entries.

### Version V04-02-03

Contains some new features:

- Dynamic entries check added. Due to rejected local fits degrees of freedoms may get lost. The entries from accepted local fits are counted and checked for each internal iteration.
- Weighted constraints implemented with new option `weightedcons`.

### Version V04-02-02

Contains some new features:

- Printout of global parameter counts from binary files added to result file with new option `printcounts`.
- Aborting in case of (C) read errors.
- Storage of values read from text files as doubles implemented.

### Version V04-02-01

Integer overflow in auto detection of bordered band matrices for local fits with huge number of parameters fixed.

### Version V04-02-00

The reading of records from C binary files containing doubles (indicated by negative record length) instead of floats has been implemented. The complete local fit is now in double precision.

### Version V04-01-04

Integer overflow in auto detection of bordered band matrices for local fits with huge number of parameters fixed.

### Version V04-01-03

Fixed V04-01-02.

### Version V04-01-02

Has bug, crashing with gcc-4.8.

Contains some minor fixes, e.g:

- Matrix inversion for local fit uses tolerance as specified by compiler instead of fixed value.
- Python script (`tools/readMilleBinary.py`) to print binary files accepts now command line arguments (file name, number of records).

### Version V04-01-01

The exit code and message of the **pede** executable are written to the file `millepede.end`:

- \-1: Still running or crashed\<
- 00: Ended normally
- 01: Ended with warnings (bad measurements)
- 02: Ended with severe warnings (bad global matrix)
- 10-17: Aborted, problem with steering
- 20-26: Aborted, problem with binary data
- 30-34: Aborted, internal problem

### Version V04-01-00

Solution method [MINRES-QLP](http://www.stanford.edu/group/SOL/software/minresqlp/) added. This advanced MINRES version is selected by '`method sparseMINRES-QLP ..`' or '`method fullMINRES-QLP ..`'.

### Version V04-00-04

Fixes for integrity check of binary files (program abort) and reading of MINRES tolerance parameter (`mrestol` set to zero).

### Version V04-00-03

Added checks of integrity of binary files. The integers are checked for the consistency of labels and markers (zeros) and the floats for NaNs (not a number).

### Version V04-00-02

Allow parallel make and other minor changes to Makefile, e.g.remove option for LARGE integers to be INTEGER\*4 to allow compilation on 32-bit systems (but medium size problem now require a bit more memory).

### Version V04-00-01

The _line search_ to optimize the solution is now configurable. The new option `linesearch` allows to select the iterations (solutions) with line search: \>2: all, =2: all with (next) Chi2 cut scaling factor equal to 1., =1: last, \<1: none. The default value is '2', the previous behaviour corresponds to '3'.

Minor fixes. Option '`maxrecord`' working again with OpenMP.

### Version V04-00-00

First "fortran90" version. The source code has been transformed from fixed form (punch card image) Fortran77 to free form fortran90. Major changes:

- Use of dynamic memory allocation in fortran90 instead of fixed size COMMON block for dynamic data structures (matrices, vectors). Only one **Pede** executable is sufficient for all application sizes. This requires gcc version 4.4 or higher.
- Documentation with _doxygen_. This includes the draft manual from V. Blobel (2007), the major changes with respect to that and the list of **Pede** options and commands. It is available in [html](http://www.desy.de/~kleinwrt/MP2/doc/html/index.html).

---

### Version V03-04-05

Bug fixes:

- Linesearch: Use proper size of vectors (NVGB instead of NAGB).
- Measurements defined in steering files: Label not replaced anymore by index during construction of sparsity structure of global matrix.

### Version V03-04-04

- Adjustments/comments to Makefile for more recent gcc (e.g. 4.6.X).
- New option '`memorydebug`' (1: just printing, \>1: more debug).
- Little bugfix, cleanup of algebra routines.
- Add plot "\<ndf\> vs file".

### Version V03-04-03

Minor bug fix to tools/readPedeHists.C to read pede histogram output.

### Version V03-04-02

Bug fixes: (global matrix) update cache buffer overrun from rejected records; integer overflow for final rejection ratio.

### Version V03-04-01

Compilation problem for gcc version 4.4/4.5 avoided for MINRES pre-conditioner. New features:

- Method MINRES makes use of pre-conditioners for error calculation (option '`errlabels`') too.
- Method INVERSION calculates global correlations for option '`globalcorr`'.

### Version V03-04-00

New features:

- Methods INVERSION and fullMINRES parallelized.
- Interfaced to [GeneralBrokenLines](https://gitlab.desy.de/claus.kleinwort/general-broken-lines/-/wikis/home) (external track seed).
- For insufficient constraints switch to SUBITO (-s) mode to avoid dead loops (unless option FORCE (-f) specified).

### Version V03-03-00

Cleanup of printout and new features:

- For method MINRES new option '`errlabels <lab1> .. <labn>`' to define up to 100 global parameters for which errors are calculated.
- More parts of the code (like the local fit) are parallelized with OpenMP. With '`threads <nproc> <nread>`' the number of threads for processing (nproc) and reading binary files (nread) in parallel can be specified.
- With Makefile option 'SUPPORT_ZLIB = yes' gzipped binary files can be read directly.
- Weights have been implemented for the records. They are defined per binary file with the file name '\<filename\> -- \<weight\>'.

### Version V03-02-02

Methods inversion and diagonalization working again!

### Version V03-02-01 (MINRES only)

Methods inversion and diagonalization broken!

- Bug fix: nsingle=2^n-1 now possible for option '`pairentries`'.
- New option '`maxrecord <record>`' to set record limit for input from binary files.
- New option '`matmoni <record>`' to monitor construction of sparse global matrix. Every \<record\> records the number of off-diagonal elements, the compression ratio and the memory consumption is printed.
- Cleanup: For OpenMP the maximal number of threads is set at run time. For the final Chi2/ndf the number of degrees of freedom (ndf) is corrected for the number of global parameters.

### Version V03-02-00 (MINRES only)

Methods inversion and diagonalization broken!

- Optimization of memory usage for sparse storage of global matrix. With option '`compress`' the sparsity information is compressed and off-diagonal elements with single entries are stored in single precision (lossless as input from binary file is single precision only). With '`pairentries <ncut> <nhist> <nsingle>`' off-diagonals with up to \<nsingle\> entries are stored in single precision (unintended feature: for nsingle=2^n-1 and nhist\<=nsingle all off-diagonals are in single precision).
- Parallel processing for solution with '`method sparseMINRES`'. Most time is spent in AVPROD (multiplication of global matrix with vector). This has been parallelized with [OpenMP™](http://openmp.org/wp). The number of threads are specified with '`threads <nthread>`'.
- The default of the internal convergence parameter RTOL for MINRES has been changed from 10^-8 to 10^-6. It can be set by '`mrestol <rtol>`'.

### Version V03-01-00

Large Integers are used for internal (FORTRAN) memory management. The Makefile allows to select 32 bit (LARGE_SIZE=4) or 64 bit (LARGE_SIZE=8) integers.

### Version V03-00-02 (Last 32 bit recommended version)

- (V04-00-02 should work on 32-bit systems as well)
- several bug fixes
- cut '`pairentries <ncut> <nhist>`' implemented. Rare global parameter pairs (off-diagonal elements of global matrix) with less than \<ncut\> entries are neglected. Optional \<nhist\> defines upper bound for histogramming.

### Version V03-00-01

- fix that additional text steering files stated after the first binary file are silently ignored (bug present since V02-00-00)
- small fix for auto detection of bordered band matrices
- '`chisqcut f1 f2`' with f1\<1.0 will disable Chi2 cut (including Huge) in FIRST iteration
- '`subito`' option reimplemented (lost since V02-00-00)
- '`entries n`' checks for n\>0 (n=max(1,n))
- for `debug`' some more statistics for fixed parameters is calculated (fraction of equations, records with fixed parameters)

### Version V03-00-00

- local (track) fit improvements:
  + fast methods for band matrix structure coming from broken lines track fit
  + fix single precision leak in local fit
  + reject fits with rank deficit
- new inline simulation
  + simple silicon strip tracker
  + multiple scattering in detector material
  + track models with different accounting of multiple scattering (-t=SL0/SLE/BP/BRLC/BRLF, details in mptst2.F)
- resurrected Eigenvalue output lost in V02-00-00
- histogram output:
  + separate residuals for measurements with and without global parameters (e.g. for multiple scattering treatment via pseudo measurements)
  + add histograms for pulls of measurements
- detect NAN to avoid endless loops
- tool readMilleBinary.py to deal also with Fortran binaries
- improved some text output

---

### Version V02-00-01

- Fix histogram output (histograms partly empty due to bug introduced when cleaning compiler warnings).
- Moved macro readPedeHists.C to new 'tools' subdirectory, added python script readMilleBinary.py to dump binary files.

### Version V02-00-00

Recent development has been presented at the [3rd LHC Detector Alignment Workshop in June 2009](http://indico.cern.ch/conferenceDisplay.py?confId=50502) by Prof. Blobel. The main new features/changes are:

- Solution methods:
  + iterative solution method named correctly MINRES, not GMRES,
  + remaining methods are '`sparseMINRES`', '`fullMINRES`', '`inversion`' and '`diagonalisation`', others have been removed.
- New preconditioning for MINRES method:
  + not limiting the number of parameters to 46340,
  + preconditioning always applied (except if '`bandwidth 0`' specified),
  + recommendation: do _not_ specify bandwidth option (uses bandwith 1),
  + if bandwith option explicitely requested, another method will be used.
- Fix problem in hash values for label pairs that was limiting to 46340 parameters.
- The line search algorithm that is applied after calculating the step direction in a millepede iteration looks for slope = 0, not a function value minimum. This achieves better stability with outlier rejection.
- Speed up of LOOP2 for sparse matrix storage.
- Further new options:
  + '`matiter <number>`': Force re-calculation of the matrix in the first \<number\> iterations - slower, but more robust in case of large number of rejects due to too many outliers.
  + '`presigma <number>`': Assign pre-sigma \<number\> to all labels where no explicit presigma is set.
  + '`regularisation <number>`: Regularisation with τ = \<number\>, i.e. add τ ∙ ‖x‖ to the minimised function where x is the vector of alignment parameters, weighted with the inverse of their individual presigma; this should reduce the influence of weak modes, but needs careful studies.
  + '`regularisation <number1> <number2>`: Same as having both '`regularisation <number1>`' and '`presigma <number2>`'.
  + '`hugecut <number>`': to determine huge outlier rejection cut that is always active, irrespective of chisqcut settings, default stays at 50, minimum value is 1.0.
- Steering files now must end with '`.txt`'.
- The subversion repository includes the macro readPedeHists.C to convert the pede histogram file millepede.his into ROOT format and/or produce a postscript file from it.

The manual still needs to be updated.
