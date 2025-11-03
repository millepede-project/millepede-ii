
!> \file
!! Millepede II program, inline changelog for doxy page.
!!
!! \section news_sec News
!! * 131008: New solution method \ref ch-minresqlp "MINRES-QLP"
!! [\ref ref_sec "ref 9"] implemented.
!! * 140226: Reading of C binary files containing *doubles* implemented.
!! * 141020: Storage of values read from text files as *doubles* implemented.
!! * 141125: Dynamic entries (from accepted local fits) check implemented.
!!   (Rejection of local fits may lead to the loss of degrees of freedom.)
!!   Printout of global parameter counters with new command \ref cmd-printcounts.
!! * 141126: Weighted constraints implemented (with new command \ref cmd-weightedcons).
!! * 150210: Solution by elimination for problems with linear equality constraints
!!   has been implemented (as default, new command \ref cmd-withelim) in addition to the
!!   Lagrange multiplier method (new command \ref cmd-withmult).
!! * 150218: Skipping *empty* constraints (without variable parameters).
!!   With new command \ref cmd-checkinput detailed check of input data (binary files,
!!   constraints) is performed, but no solution will be determined.
!!   Some input statistics is available in the output file <tt>millepede.res</tt>.
!! * 150226: Iteration of entries cut with new command \ref cmd-iterateentries.
!!   In the second iteration measurements with any parameters fixed by the 
!!   previous entries cut are skipped. Useful if parameters of measurements have
!!   different number of entries. 
!! * 150420: Skipping of empty constraints has to be enabled by new command \ref
!!   cmd-skipemptycons.
!! * 150901: Preconditioning for MINRES with skyline matrix (avoiding rank deficits 
!!   of band matrix) added (selected by second argument in \ref cmd-bandwidth >0).
!! * 150925: Monitoring of residuals per local fit cycle is selected by \ref cmd-monres.
!!   The normalized residuals are grouped by the first global label and the median 
!!   and the RMS (from the median of the absolute deviations) per group are
!!   written to <tt>millepede.mon</tt>.
!! * 170502: Monitoring of pulls per local fit cycle is selected by \ref cmd-monpull.
!!   The scaling of measurement errors is enabled by \ref cmd-scaleerrors.
!!   Pede will abort now for constraints with a singular QL decomposition
!!   of the constraints matrix (solution by elemination).
!!   This problem is usually caused by *empty* constraints (see \ref
!!   cmd-skipemptycons).
!! * 170831: More debug information for problems with reading Cfiles. Don't stop
!!   after read error for \ref cmd-checkinput mode.
!! * 180525: Some fixes: Proper handling of special (debug) data blocks in binary
!!   files, proper exit code (3) for 'function not decreasing'.
!! * 180815: Some minor fixes, additional level of detail (appearance range of global
!!   parameters in binary files) for \ref cmd-checkinput mode.
!! * 190319: Constraints are now sorted and split into disjoint blocks to speed up
!!   calculation of rank and QL decomposition by block matrix algebra.
!!   This works best if the label sets of the involved alignable objects are disjoint too.
!! * 190412: Cleanup of operations (open, close, rewind) on binary files. New command
!!   \ref cmd-closeandreopen to enable closing and reopening of binary files
!!   to limit the number of concurrently open files. The modification dates of the
!!   files are monitored to ensure data integrity.
!! * 190430: Update of (approximate) string matching for keyword detection.
!!   Matching is now symmetric in pattern and text. Previously e.g. a binary file
!!   with the letters from '<tt>Cfiles</tt>' in the name in that order was
!!   treated as that keyword and not as a binary file.
!! * 191004: Checking global parameters for disjoint blocks. In case of solution by
!!   inversion (optionally with constraints handled by elimination) switch to
!!   \ref mpmod::npblck "block diagonal" storage mode.
!! * 200429: Modifications for compilation with PGI compiler (make -f Makefile_pgi).
!! * 200701: Implementation of \ref ch-pargroup "parameter groups" (sets of adjacent global parameters (labels)
!!   appearing in the binary files *always* together). Used to speed up construction of global matrix.
!!   Similarity operations are now aware of sparse (rectangular) matrices.
!! * 200716: The counting of the appearance of global parameters in the binary files can
!!   now be done on record (e.g. track) level instead of equation (e.g. measurement) level.
!!   This is enabled with the new command \ref cmd-countrecords and makes the iteration
!!   of the first data loop (by \ref cmd-iterateentries) obsolete.
!! * 201027: New solution method \ref ch-mchdec "decomposition" implemented.
!! * 201214: New command \ref cmd-monpgs to monitor progress in operations
!!   on global and constraints matrices.
!! * 210301: New solution methods \ref ch-lapack "fullLAPACK" and \ref ch-lapack "unpackedLAPACK"
!!   (matrix factorization) based on [LAPACK](http://www.netlib.org/lapack/) can be included optionally
!!   (at compile time, <tt>-DLAPACK64=..</tt>).
!! * 210728: Exploit decomposition of constraints matrix into disjoint blocks for all
!!   solution methods (e.g. QL decomposition, MINRES preconditioner).
!! * 211008: Fortran code modernized (EQUIVALENCE and ENTRY statements replaced) and checked
!!   (compiling with '-fcheck=all').
!! * 211022: Fortran code modernized further (assumed-size array arguments replaced).
!! * 211101: Migration from DESY \c svn to \c gitlab server (includes wiki).
!! * 211210: Further exploration of sparsity of constraints matrix (C). First the constraints are now
!!   split into disjoint *groups* for on optimized check of the rank of the product matrix (C*C^t).
!!   The groups are then combined into non overlapping *blocks* for an efficient QL decompsition in
!!   case of elimination of constraints. For solution by "unpackedLAPACK" the QL decomposition is
!!   now using the internal sparsity-aware code as default. To use LAPACK routines for this the
!!   new command \ref cmd-withlapackelim has to be used.
!! * 211222: Constraints groups included in \ref cmd-checkinput. Documentation for
!!   \ref ch-checkinput and \ref troubleshooting_page added.
!! * 220616: Cleanup, test programs to print LAPACK (library) configuration added to \c tools directory.
!! * 220817: Cleanup, (rare) problem with construction of \ref ch-pargroup "parameter groups" fixed
!!   (avoiding aborts with exit code 35).
!! * 221010: Fix (uninitialsed values) and cleanup for internal silicon strip tracker example.
!! * 221017: More code modernisation to comply the with fortran standard 2018
!!   (<tt>gcc11 -std=f2018 -fall-intrinsics</tt>). Still some GNU fortran extensions are used:
!!   <tt>etime, fdate, getarg, getenv, iargc, stat, system, time</tt>
!! * 221122: Cleanup and documentation/exercises for \ref test_brlf_page "example"
!!   (internal test case -t=BRLF).
!! * 221212: Force \ref ch-checkinput "check input mode" (2) in case of accepted *empty* constraints
!!   (no variable parameters). No solution will be calculated.
!! * 230201: Fix global parameter errors for solution by diagonalization using elimination of constraints.
!! * 230321: Fortran/C interoperability uses now 'iso_c_binding' (fortran 2003) instead of 'cfortran.h'.
!! * 230322: Tool \c readMilleBinary.py now compatible with python3, updated CLI.
!! * 230502: Cleanup and fixes (constraint elimination with LAPACK using constraint groups,
!!   proper termination of file names for 'iso_c_binding')
!! * 230515: Check for redundancy constraints: Constraint groups defining linear transformation between two groups
!!   of equivalent global parameters. With the new command \ref cmd-resolveredundancycons they can be resolved
!!   (to save resources and burden on numerics).
!!   Optionally add (brief) \ref  ch-parcom "comments" for global parameters to annotate the results file.
!! * 230516: New command \ref cmd-checkpargroups to check (the rank (linear independency of
!!   global derivatives) for) (global) \ref ch-pargroup "parameter groups".
!! * 230617: Fix problem with monitoring of residuals. Calculate *skyline* fraction for sparse matrices.
!! * 230822: Fix problem for block diagonal global matrix (keeping single block).
!! * 231020: Define proper \ref par-linesearch "line search" parameters for LAPACK too.
!! * 231218: New optional method \ref ch-pardiso "sparsePARDISO"
!!   using the Intel oneMKL PARDISO solver for sparse matrices.
!! * 240214: Fix severe problem with external measurements depending on multiple global parameters.
!! * 240227: Quick fix for possible integer overflow in summing up global Chi2 (ADDSUM). Needs careful revision.
!! * 240229: Summation of global Chi2 and NDF revised (no more 32bit variables, update per record).
!! * 240412: Counters scaling with the number of records are now long (64bit) integers.
!! * 240429: Added \c tinypede.py to tools.
!! * 240502: Added \ref Legacy "legacy" (Millepede-I) folder.
!! * 240610: In first loop over binary files for counting and grouping of global labels
!!   ignore those with zero global derivative. (Option \ref cmd-printcounts "printcounts -1" will print their numbers.)
!! * 240626: For local fits with bordered-band matrix structure use condition of diagonal matrix from
!!   root-free Cholesky decomposition (of band part) to optionally reject records (see \ref cmd-maxlocalcond
!!   and internal histogram 16).
!! * 240708: Code complies the with fortran standard 2023
!!   (<tt>gcc14 -std=f2023 -fall-intrinsics</tt>). Still GNU fortran extensions are used.
!! * 240716: Modernisation of development environment (from EL7 to EL9 (gcc11.4),
!!   from ompP to [Score-P](http:score-p.org) for profiling).
!! * 240731: For instrumentation and profiling with [Score-P](http:score-p.org) switched
!!   from POMP (<tt>scorep --pomp</tt>) (based on OPARI2, to be superseded by OMPT (OpenMP 5.0))
!!   to user (<tt>scorep --user</tt>) regions.
!! * 240918: Allow for post processing of results (option \ref cmd-postprocessing).
!! * 241205: Proper abort (message) in case of wrong binary file type (Fortranfiles or Cfiles).
!! * 250306: Some tools have been reimplemented in \ref julia_sec.
!! * 250819: Some optimizations for \c tinypede.jl tool.
!! * 251013: New command \ref cmd-outlierfracwarnthreshold "outlierfracwarnthreshold" to set warning threshold for outlier fraction. readMilleBinary to py3, small cosmetic changes to readc.c 
