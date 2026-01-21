 \page option_page List of options and commands

 \tableofcontents

 # Command line options:
 ## -t
 Create text and binary files for \ref mptest1.f90 "wire chamber" test case, set
 \ref mpmod::ictest "ictest" to 1.
 ## -t=track-model
 Create text and binary files for \ref mptest2.f90 "silicon strip tracker" test case
 using \a track-models with different accounting for multiple scattering, set
 \ref mpmod::ictest "ictest" to 2..6.
 ## -s
 Solution is not iterated.
 Automatically switched on in case of rank deficits for constraints.
 ## -f
 Force iterating of solution (in case of rank deficits for constraints).
 ## -c
 Check input (binary files, constraints). No solution is determined. (\ref mpmod::icheck "icheck"=1)
 ## -C
 Check input (binary files, constraints, appearance). No solution is determined. (\ref mpmod::icheck "icheck"=2)
 ## -n=[Nevents]
 If running with one of the [-t](#-t) options, generate `Nevents` test records (default: 10000).  

 # Steering file commands:
 In general the commands are defined by a single line:

         keyword   number1  number2  ...

 For those specifying \ref sssec-parinf "properties" of the global parameters
 (\a keyword = \c parameter, \c constraint or \c measurement (or \c comment))
 for each involved global parameter (identified by a \ref an-glolab "label")
 one additional line follows:

         label     number1  number2  ...

 Default values for the numerical arguments are shown in
 the command descriptions in '[]'. Missing arguments without default
 values have no effect.

 ## bandwidth
 Set band width \ref mpmod::mbandw "mbandw" for
 \ref minresmodule::minres "MINRES" preconditioner to \a number1 [0]
 and additional flag \ref mpmod::lprecm "lprecm" to \a number2 [0].
 ## blocksizePARDISO
 Enable for \ref ch-pardiso "PARDISO" storage in BSR3 format
 (three-array variation of the block compressed sparse row format).
 Provide (one to ten, increasingly ordered) candidate block sizes. The block size with the smallest
 memory footprint (of the global matrix in BSR3 format) is selected.
 ## cache
 Set (read+write) cache size \ref mpmod::ncache "ncache" to \a number1.
 Define cache size and average fill level.
 ## Cfiles
 Following binaries are C files - this refers to any format supported by the external Mille (anything not Fortran-binary).
 ## checkinput
 Set check input flag \ref mpmod::icheck "icheck" to \a number1 [1].
 Similar to \ref opt-c "-c" or \ref opt-C "-C".
 For mpmod::icheck "icheck" >0 no solution is performed but input statistics is checked in detail.
 With mpmod::icheck "icheck" >1 the appearance range (first/last file,record and number of files)
 of global parameters is determined too.
 ## cmd-checkpargroups checkparametergroups
 Set flag \ref mpmod::ichkpg "ichkpg" to 1 (true) to enable the checking of
 (the rank of) \ref ch-pargroup "parameter groups".
 ## chisqcut
 For local fit \ref an-chisq "setChi^2" cut \ref mpmod::chicut "chicut" to \a number1 [1.],
 \ref mpmod::chirem "chirem" to \a number2 [1.].
 ## comment
 Define \ref ch-parcom "comments" for global parameters.
 ## compress
 Obsolete. Compression is default.
 ## closeandreopen
 Set flag \ref mpmod::keepopen "keepOpen" to zero to enable closing and reopening of binary files
 to limit the number of concurrently open files.
 ## constraint
 Define \ref sssec_consinf "constraints" for global parameters.
 ## countrecords
 Set flag \ref mpmod::mcount "mcount" to 1 (true) to enable parameter counting om record level.
 ## debug
 Set number of records with debug printout \ref mpmod::mdebug "mdebug" to
 \a number1 [3], number of measurements with printout \ref mpmod::mdebg2 "mdebg2" to \a number2.
 ## debugPARDISO
 Set \ref ch-pardiso "PARDISO" message level \ref mpmod::ipddbg "ipddbg" to
 \a number1 [0]. Dump steering array \c IPARM(64) for value > 0.
 ##  dwfractioncut
 Set \ref an-dwcut "down-weighting fraction" cut \ref mpmod::dwcut "dwcut"
 to \a number1 (max. 0.5).
 ## outlierfracwarnthreshold
 Set \ref an-outlierfracwarnthreshold "warning level for fraction of large chi2 entries" \ref mpmod::warnthresholdchi2 "warnthresholdchi2"
 to \a number1 [1]. Value is expected to be in percent (1.0 parsed as 1%).
 ## entries
 Set \ref an-entries "entries" cuts for variable global parameter
 \ref mpmod::mreqenf "mreqenf" to \a number1 [25],
 \ref mpmod::mreqena "mreqena" to \a number2 [10] and
 \ref mpmod::iteren "iteren" to the product of \a number1 and \a number3 [0].
 ## errlabels
 Define (up to 100 in total) global labels \a number1 .. \a numberN
 for which the parameter errors are calculated for method MINRES too
 (by \ref solglo "solving" \f$\mathbf{C}\cdot\vec{x}_i = \vec{b}^i, b^i_j = \delta_{ij} \f$).
 ## force
 Set force (iterations) flag \ref mpmod::iforce "iforce" to 1 (true).
 Same as [-f](#-f).
 ## fortranfiles
 Following binaries are Fortran files.
 ## globalcorr
 Set flag \ref mpmod::igcorr "igcorr" for output of global correlations to 1 (true).
 ## histprint
 Set flag \ref mpmod::nhistp "nhistp" for \ref an-histpr "histogram printout"
 to 1 (true).
 ## hugecut
 For local fit set Chi^2 cut \ref mpmod::chhuge "chhuge"
 for \ref sssec-outlierdeb "unreasonable data" to \a number1 [1.].
 ## iterateentries
 Set maximum value \ref mpmod::iteren "iteren" for iteration of entries cut to
 \a number1 [maxint]. Can alternatively be set by the [entries](#entries) command.
 For parameters with less entries the cut will be iterated ignoring measurements with
 at least one parameter below \ref mpmod::mreqenf "mreqenf".
 ## lapackwitherrors
 Set flag \ref mpmod::ilperr "ilperr" for calculation of inverse matrix
 for parameter errors by \ref ch-lapack "LAPACK" to 1 (true).
 ## linesearch
 The mode \ref mpmod::lsearch "lsearch" of the \ref par-linesearch "line search"
 to improve the solution is set to \a number1.
 ## localfit
 For local fit set number of iterations \ref mpmod::lfitnp "lfitnp"
 with calculation of pulls to \a number1, flag \ref mpmod::lfitbb "lfitbb"
 for auto-detection of bordered band matrices to \a number2.
 ## matiter
 Set number of iterations \ref mpmod::matrit "matrit" with (re)calcuation of
 global matrix to \a number1.
 ## matmoni
 Set record interval \ref mpmod::matmon "matmon" for monitoring of (sparse) matrix
 construction to \a number1.
 ## maxlocalcond
 Set maximal Log10(condition) of decomposition of band part for local fit
 \ref mpmod::cndlmx "cndlmx" to \a number1. Records with larger condition will be rejected.
 ## maxrecord
 Set record limit \ref mpmod::mxrec "mxrec" to \a number1.
 ## measurement
 Define (additional) \ref sssec_gpm "measurements" for global parameters.
 ## memorydebug
 Set debug flag \ref mpmod::memdbg "memdbg" for memory management
 to \a number1 [1].
 ## method
 Has special format:

         method   name     number1  number2

 Set \ref ch-methods "solution method" \ref mpmod::metsol "metsol" and
 storage mode \ref mpmod::matsto "matsto" according to \a name,
 (\c inversion : (1,1), \c diagonalization : (2,1),
 \c decomposition : (3,1),
 \c fullMINRES : (4,1) or \c sparseMINRES : (4,2),
 \c fullMINRES-QLP : (5,1) or \c sparseMINRES-QLP : (5,2),
 \c fullLAPACK factorization : (7,1), \c unpackedLAPACK factorization : (8,0)),
 \c sparsePARDISO factorization : (9,3)
 (minimum) number of iterations \ref mpmod::mitera "mitera" to \a number1,
 convergence limit \ref mpmod::dflim "dflim" to \a number2.

 \c Inversion and \c diagonalization provide in addition to the solution the parameter errors
 (from the diagonal of the inverted global matrix). Solutions with \c MINRES are only approximate.

 ## monitorresiduals
 Set flag \ref mpmod::imonit "imonit" for monitoring of residuals to \a number1 [3]
 and increase number of bins (of size 0.1) for internal storage to \a number2 [100].
 Monitoring mode \ref mpmod::imonmd "imonmd" is 0.
 ## monitorpulls
 Set flag \ref mpmod::imonit "imonit" for monitoring of pulls to \a number1 [3]
 and increase number of bins (of size 0.1) for internal storage to \a number2 [100].
 Monitoring mode \ref mpmod::imonmd "imonmd" is 1.
 ## monitorprogress
 For progress monitoring set for repetition rate \c nrep the start value \ref mpmod::monpg1 "monpg1"
 to \a number1 [1] and maximum increase \ref mpmod::monpg2 "monpg2" to \a number2 [1024].
 Monitored are operations (inversion, decomposition, similarity) on the global and the constraints matrices.
 If the (outermost loop) index is greater equal \c nrep the index is printed and \c nrep updated
 (+ min(\c nrep, \c monpg2)).
 ## mresmode
 Set \ref minresqlpmodule::minresqlp "MINRES-QLP" factorization mode
  \ref mpmod::mrmode "mrmode" to \a number1.
 ## mrestranscond
 Set \ref minresqlpmodule::minresqlp "MINRES-QLP" transition (matrix) condition
  \ref mpmod::mrtcnd "mrtcnd" to \a number1.
 ## mrestol
 Set tolerance criterion \ref mpmod::mrestl "mrestl" for \ref minresmodule::minres "MINRES"
 to \a number1 (\f$10^{-10}\f$ .. \f$10^{-4}\f$).
 ## nofeasiblestart
 Set flag \ref mpmod::nofeas "nofeas" for \ref an-nofeas "skipping"
 making parameters feasible to \a number1 [1].
 ## outlierdownweighting
 For local fit set number of \ref sssec-outlow "outlier"
 \ref an-downw "down-weighting" iterations
 \ref mpmod::lhuber "lhuber" to \a number1.
 ## pairentries
 Set entries cut for variable global parameter pairs \ref mpmod::mreqpe "mreqpe"
 to \a number1, histogram upper bound \ref mpmod::mhispe "mhispe" for pairs
 to \a number2 (<1: no histogramming), upper bound \ref mpmod::msngpe "msngpe"
 for pair entries with single precision storage
 to \a number3.
 ## parameter
 Define \ref sssec-parinf "initial value, pre-sigma" for global parameters.
 ## pardiso
 Modify for [PARDISO](changes_page.html#pardiso) the internal steering parameters.
 ## postprocessing
 Define post processing *string*. Will be executed by system at end of **pede**.
 ## presigma
 Set default pre-sigma \ref mpmod::defaultPreSigma "defaultPreSigma" to \a number1 [1].
 ## print
 Set print level \ref mpmod::mprint "mprint" to \a number1 [1].
 ## printcounts
 Set flag \ref mpmod::ipcntr "ipcntr" to \a number1 [1].
 The counters for the global parameters from the accepted local fits (=1)
 or from the binary files (>1) will be printed in the result file. Alternatively
 the counters for zero global derivatives from the binary files (<0) can be selected.
 ## printrecord
 \ref an-recpri "Record" numbers with printout.
 ## pullrange
 Set (symmetric) range \ref mpmod::prange "prange" for histograms
 of pulls, normalized residuals to \a number1 (=0: auto-ranging).
 ## readerroraseof
 Set flag \ref mpmod::ireeof "ireeof" to 1 (true) to treat read errors for binary files
 as end-of-file instead of aborting.
 ## regularisation / regularization
 Set flag \ref mpmod::nregul "nregul" for regularization to 1 (true),
 regularization parameter \ref mpmod::regula "regula" to \a number2,
 default pre-sigma \ref mpmod::defaultPreSigma "defaultPreSigma" to \a number3.
 ## resolveredundancycons
 Set flag \ref mpmod::irslvrc "irslvrc" to 1 (true).
 Redundancy constraints will be resolved 
 (parameters appearing in constraints will be fixed, constraints skipped).
 ## scaleerrors
 Set measurement scaling factors \ref mpmod::errorScaleFactor "errorScaleFactor"
 to \a number1 [1.] and \a number2 [\a number1].
 First value is for "global" measurements (with global derivatives),
 second for "local" measurements (without global derivatives).
 ## skipemptycons
 Set flag \ref mpmod::iskpec "iskpec" to 1 (true).
 Empty constraints (without variable parameters) will be skipped.
 ## subito
 Set subito (no iterations) flag \ref mpmod::isubit "isubit" to 1 (true).
 Same as [-s](#-s).
 ## threads
 Set number \ref mpmod::nOMPThreads "nOMPThreads" of OpenMP&tm; threads for processing
 to \a number1,
 number \ref mpmod::numberOfReadingThreads "numberOfReadingThreads" of threads for reading
 binary files to \a number2 [\a number1].
 ## weightedcons
 Set flag \ref mpmod::iwcons "iwcons" to \a number1 [1].
 Implements \ref sssec_consinf "weighted constraints" for global parameters.
 ## withelimination
 Set flag \ref mpmod::icelim "icelim" to 1 (true).
 Selects solution by elimination for linear equality constraints.
 ## withlapackelimination
 Set flag \ref mpmod::icelim "icelim" to 2 (LAPACK).
 Selects solution by elimination for linear equality constraints with LAPACK.
 Only available for unpacked LAPACK!
 ## withmultipliers
 Set flag \ref mpmod::icelim "icelim" to 0 (false).
 Selects solution by Lagrange multipliers for linear equality constraints.
 ## wolfe
 For strong Wolfe condition in \ref par-linesearch "line search"
 set parameter \ref mpmod::wolfc1 "wolfc1" to \a number1, \ref mpmod::wolfc2
 "wolfc2" to \a number2.
