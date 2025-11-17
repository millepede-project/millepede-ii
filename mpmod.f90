!> \file
!! Data structures.
!!
!! \author Claus Kleinwort, DESY, 2012 (Claus.Kleinwort@desy.de)
!!
!! \copyright
!! Copyright (c) 2012 - 2024 Deutsches Elektronen-Synchroton,
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

!> Parameters, variables, dynamic arrays.
!!
!! For parameters which can be set from command line or
!! steering files more details are available in: \ref option_page.

MODULE mpmod
    USE mpdef
    IMPLICIT NONE
    SAVE
    ! steering parameters
    INTEGER(mpi) :: ictest=0  !< test mode '-t'
    INTEGER(mpi) :: ntestEvents=10000  !< test stats '-n'
    INTEGER(mpi) :: metsol=0  !< solution method (1: inversion, 2: diagonalization, 3: decomposition, 4: MINRES, 5: \ref minresqlpmodule::minresqlp "MINRES-QLP", 7/8: LAPACK, 9: (Intel oneMKL) PARDISO)
    INTEGER(mpi) :: matsto=2  !< (global) matrix storage mode (0: unpacked, 1: full = packed, 2: sparse(custom), 3: sparse(CSR3, BSR3))
    INTEGER(mpi) :: matbsz=1  !< (global) matrix (fixed) block size, only used for BSR3 storage mode (Intel oneMKL PARDISO)
    INTEGER(mpi) :: mprint=1  !< print flag (0: minimal, 1: normal, >1: more)
    INTEGER(mpi) :: mdebug=0  !< debug flag (number of records to print)
    INTEGER(mpi) :: mdebg2=10 !< number of measurements for record debug printout
    INTEGER(mpi) :: mreqenf=25 !< required number of entries (for variable global parameter from binary Files)
    INTEGER(mpi) :: mreqena=10 !< required number of entries (for variable global parameter from Accepted local fits)
    INTEGER(mpi) :: mitera=1  !< number of iterations
    INTEGER(mpi) :: nloopn=0  !< number of data reading, fitting loops
    INTEGER(mpi) :: mbandw=0  !< band width of preconditioner matrix
    INTEGER(mpi) :: lprecm=0  !< additional flag for preconditioner (band) matrix (>0: preserve rank by skyline matrix)
    INTEGER(mpi) :: lunkno=0  !< flag for unkown keywords
    INTEGER(mpi) :: lhuber=0  !< Huber down-weighting flag
    REAL(mps)    :: chicut=0.0  !< cut in terms of 3-sigma cut, first iteration
    REAL(mps)    :: chirem=0.0  !< cut in terms of 3-sigma cut, other iterations, approaching 1.
    REAL(mps)    :: chhuge=50.0 !< cut in terms of 3-sigma for unreasonable data, all iterations
    INTEGER(mpl) :: nrecpr=0  !< record number with printout
    INTEGER(mpl) :: nrecp2=0  !< record number with printout
    INTEGER(mpl) :: nrec1 =0  !< record number with largest residual
    INTEGER(mpl) :: nrec2 =0  !< record number with largest chi^2/Ndf
    REAL(mps)    :: value1=0.0!< largest residual
    REAL(mps)    :: value2=0.0!< largest chi^2/Ndf
    REAL(mps)    :: dwcut=0.0 !< down-weight fraction cut
    REAL(mps)    :: warnThresholdChi2=100  !< warning threshold for large chi2 tracks, internally stored in units of 1e-4 (100 = 1%)
    INTEGER(mpi) :: isubit=0  !< subito flag '-s'
    REAL(mps)    :: wolfc1=0.0!< C_1 of strong Wolfe condition
    REAL(mps)    :: wolfc2=0.0!< C_2 of strong Wolfe condition
    REAL(mpd) :: mrestl=1.0E-06 !< tolerance criterion for MINRES-QLP
    REAL(mpd) :: mrtcnd=1.0E+07 !< transition (QR -> QLP) (matrix) condition for MINRES-QLP
    INTEGER(mpi) :: mrmode=0  !< MINRES-QLP mode (0: QR+QLP, 1: only QR, 2: only QLP factorization)
    INTEGER(mpi) :: nofeas=0  !< flag for skipping making parameters feasible
    INTEGER(mpi) :: nhistp=0  !< flag for histogram printout
    REAL(mps)    :: delfun=0.0!< expected function change
    REAL(mps)    :: actfun=0.0!< actual function change
    REAL(mps)    :: angras=0.0!< angle between gradient and search direction
    INTEGER(mpi) :: iterat=0  !< iterations in solution
    INTEGER(mpi) :: nregul=0  !< regularization flag
    REAL(mps)    :: regula=1.0!< regularization parameter, add regula * norm(global par.) to objective function
    REAL(mps)    :: defaultPreSigma=0.0!< default presigma
    INTEGER(mpi) :: matrit=0  !< matrix calculation up to iteration MATRIT
    INTEGER(mpi) :: icalcm=0  !< calculation mode (for \ref xloopn "XLOOPN") , >0: calculate matrix
    INTEGER(mpi), DIMENSION(2) :: nbndr =0  !< number of records with bordered band matrix for local fit (upper/left, lower/right)
    INTEGER(mpi) :: nbdrx =0  !< max border size for local fit
    INTEGER(mpi) :: nbndx =0  !< max band width for local fit
    INTEGER(mpl) :: nrecer=0  !< record with error (rank deficit or Not-a-Number) for printout
    INTEGER(mpl) :: nrec3 = huge(nrec3) !< (1.) record number with error
    INTEGER(mpi) :: mreqpe=1  !< min number of pair entries
    INTEGER(mpi) :: mhispe=0  !< upper bound for pair entry histogrammimg
    INTEGER(mpi) :: msngpe=-1 !< upper bound for pair entry single precision storage
    INTEGER(mpi) :: mextnd=0  !< flag for extended storage (both 'halves' of sym. mat. for improved access patterns)
    INTEGER(mpi) :: nOMPThreads =1  !< number of (OpenMP) threads
    INTEGER(mpi) :: mxrec =0  !< max number of records
    INTEGER(mpi) :: matmon=0  !< record interval for monitoring of (sparse) matrix construction
    INTEGER(mpi) :: lfitnp=huge(lfitnp) !< local fit: number of iteration to calculate pulls
    INTEGER(mpi) :: lfitbb=1  !< local fit: check for bordered band matrix (if >0)
    INTEGER(mpi) :: mnrsel=0  !< number of MINRES error labels in LBMNRS (calc err, corr with SOLGLO)
    INTEGER(mpi) :: cacheBufferSize=-1 !< buffer size for caching (default 100MB per thread)
    REAL(mps), DIMENSION(3) :: fcache = (/ 0.8,  0., 0. /) !< read cache, average fill level; write cache; dynamic size
    INTEGER(mpi) :: numberOfReadingThreads=1  !< number of threads for reading binary files
    INTEGER(mpi) :: mnrsit=0  !< total number of MINRES internal iterations
    INTEGER(mpi) :: iforce=0  !< switch to SUBITO for (global) rank defects if zero
    INTEGER(mpi) :: igcorr=0  !< flag for output of global correlations for inversion, =0: none
    INTEGER(mpi) :: memdbg=0  !< debug flag for memory management
    REAL(mps)    :: prange=0.0!< range (-PRANGE..PRANGE) for histograms of pulls, norm. residuals
    INTEGER(mpi) :: lsearch=2 !< iterations (solutions) with line search:
                         !! >2: all, =2: all with (next) Chi2 cut scaling factor =1., =1: last, <1: none
    INTEGER(mpi) :: ipcntr=0  !< flag for output of global parameter counts (entries), =0: none, =1: local fits, >1: binary files
    INTEGER(mpi) :: iwcons=0  !< flag for weighting of constraints (>0: weighting with \ref globalparlabelcounter "globalParLabelCounter", else: none)
    INTEGER(mpi) :: icelim=1  !< flag for using elimination (instead of multipliers) for constraints
    INTEGER(mpi) :: icheck=0  !< flag for checking input only (no solution determined)
    INTEGER(mpi) :: ichkpg=0  !< flag for checking (rank of) parameter groups
    INTEGER(mpi) :: iteren=0  !< entries cut is iterated for parameters with less entries (if > \ref mreqenf)
    INTEGER(mpi) :: iskpec=0  !< flag for skipping empty constraints (no variable parameters)
    INTEGER(mpi) :: irslvrc=0 !< flag for resolving redundancy constraints (two equivalent parameter groups)
    INTEGER(mpi) :: imonit=0  !< flag for monitoring residuals per local fit cycle (=0: none, <0: all, bit 0: first, bit 1: last)
    INTEGER(mpi) :: measBins=100 !< number of bins per measurement for monitoring
    INTEGER(mpi) :: imonmd=0  !< monitoring mode: 0:residuals (normalized to average error), 1:pulls
    INTEGER(mpi) :: scaleErrors=0  !< flag for scaling of errors
    REAL(mpd), DIMENSION(2) :: errorScaleFactor = (/ 1.0, 1.0 /) !< scaling factors for errors of 'global' and 'local' measurement
    INTEGER(mpi) :: keepOpen=1 !< flag for keeping binary files open
    INTEGER(mpi) :: ireeof=0 !< flag for treating (binary file) read errors as end-of-file
    INTEGER(mpi) :: mcount=0 !< flag for grouping and counting global parameters on equlation (0) or record (1) level
    INTEGER(mpi) :: monpg1=0 !< progress monitoring, repetition rate start value
    INTEGER(mpi) :: monpg2=0 !< progress monitoring, repetition rate max increase
#ifdef LAPACK64
    INTEGER(mpi) :: ilperr=0 !< flag to calculate parameter errors with LAPACK
#ifdef PARDISO
    INTEGER(mpi) :: ipddbg=0 !< flag for debugging Intel oneMKL PARDISO
#endif
#endif
    REAL(mps)    :: cndlmx=100.0!< cut on log10(condition of band part) of local (bordered-band matrix) fit

    ! variables
    INTEGER(mpi) :: lunmon !< unit for monitoring output file
    INTEGER(mpi) :: lunlog !< unit for logfile
    INTEGER(mpi) :: lvllog !< log level
    INTEGER(mpi) :: nTotalGlobalPar !< total number of global parameters
    INTEGER(mpi) :: nVarGlobalPar !< number of variable global parameters
    INTEGER(mpi) :: nAllActivePar !< number of all parameters (var. global par. + Lagrange mult.)
    INTEGER(mpi) :: nFitPar !< number of fit parameters
    INTEGER(mpi) :: nConstraints !< number of constraints
    INTEGER(mpi) :: nEmptyConstraints !< number of empty constraints (no variable parameters)
    INTEGER(mpi) :: nTotalParGroups !< number of parameter groups
    INTEGER(mpi) :: nVarParGroups !< number of variable parameter groups
    INTEGER(mpi) :: nAllActiveParGroups !< number of all parameter groups (variable + Lagrange mult.)
    INTEGER(mpi) :: nParBlocks !< number of (disjoint) parameter blocks (>1: block diagonal storage)
    INTEGER(mpi) :: nConstraintGroups !< number of (disjoint) constraint groups
    INTEGER(mpi) :: nConstraintBlocks !< number of (non overlapping) constraint blocks
    INTEGER(mpl) :: mszcon !< (integrated block) matrix size for constraint matrix
    INTEGER(mpl) :: mszprd !< (integrated block) matrix size for (constraint) product matrix
    INTEGER(mpi), DIMENSION(3) :: nprecond !< number of constraints (blocks), matrix size for preconditioner
    INTEGER(mpl) :: mszpcc !< (integrated block) matrix size for constraint matrix for preconditioner
    INTEGER(mpi) :: nGlobalMaxPerRec !< max number of global paramters per record
    INTEGER(mpi) :: nLocalMaxPerRec !< max number of local paramters per record
    INTEGER(mpi) :: nEqMaxPerRec !< max number of equations (measurements) per record
    INTEGER(mpl) :: nrec  !< number of records read
    INTEGER(mpl) :: nrecd !< number of records read containing doubles
    INTEGER(mpl) :: nbEqRead  !< number of equations (measurements) read
    INTEGER(mpl) :: nbEqWithGlobPar  !< number of equations read with global parameters
    INTEGER(mpl) :: nbGlobDeriv  !< number of global derivatives read
    INTEGER(mpl) :: nbZeroGlobDeriv  !< number of zero global derivatives read
    REAL(mps)    :: dflim !< convergence limit
    INTEGER(mpl), DIMENSION(6) :: nrejec !< rejected records
    REAL(mps), DIMENSION(0:8) :: times !< cpu time counters
    REAL(mps)    :: stepl !< step length (line search)
    CHARACTER (LEN=74) :: textl !< name of current MP 'module' (step)
    LOGICAL :: newite !< flag for new iteration
    INTEGER(mpl) :: ndfsum !< sum(ndf)
    INTEGER(mpi) :: iitera !< MINRES iterations
    INTEGER(mpi) :: istopa !< MINRES istop (convergence)
    INTEGER(mpi) :: lsinfo !< line search: returned information
    REAL         :: rstart !< cpu start time for solution iterations
    REAL(mps)    :: deltim !< cpu time difference
    INTEGER(mpi) :: npresg !< number of pre-sigmas
    INTEGER(mpl) :: nrecal !< number of records
    INTEGER(mpi) :: ndefec=0 !< rank deficit for global matrix (from inversion)
    INTEGER(mpi) :: ndefpg=0 !< number of parameter groups with rank deficit (from inversion)
    INTEGER(mpi) :: nmiss1=0 !< rank deficit for constraints
    INTEGER(mpi) :: nalow=0 !< (sum of) global parameters with too few accepted entries
    INTEGER(mpi) :: nxlow=0 !< (max of) global parameters with too few accepted entries for icalcm=1
    INTEGER(mpi) :: lcalcm !< last calclation mode
    INTEGER(mpi) :: nspc=1 !< number of precision for sparse global matrix (1=D, 2=D+F)
    INTEGER(mpi) :: nencdb !< encoding info (number bits for column counter)
    INTEGER(mpi) :: numMeas !< number of measurement groups for monitoring
    REAL(mpd), PARAMETER :: measBinSize=0.1 !< bins size for monitoring 
    INTEGER(mpi), DIMENSION(100) :: lbmnrs !< MINRES error labels
    REAL(mpd) :: fvalue !< function value (chi2 sum) solution
    REAL(mpd) :: flines !< function value line search
    REAL(mpd) :: sumndf !< weighted sum(ndf)
    INTEGER(mpi) :: nrderr=0 !< number of binary files with read errors
    INTEGER(mpi) :: mpdbsz=0 !< PARDISO, number of block sizes to be tried (by PBSBITS)
    INTEGER(mpi), DIMENSION(10) :: ipdbsz !< PARDISO, list of block sizes to be tried (by PBSBITS)

    ! each loop
    INTEGER(mpi) :: nbReadRecords     !< number of buffers (records) in (read) block
    INTEGER(mpi) :: numBlocks         !< number of (read) blocks
    INTEGER(mpl) :: sumRecords        !< sum of records
    INTEGER(mpl) :: skippedRecords    !< number of skipped records (buffer too small)
    INTEGER(mpi) :: minRecordsInBlock !< min. records in block
    INTEGER(mpi) :: maxRecordsInBlock !< max. records in block
    INTEGER(mpi) :: lenGlobalVec      !< length of global vector 'b' (A*x=b)
    ! dynamic arrays
    !======================================================
    ! global parameters
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: globalParameter !< global parameters (start values + sum(x_i))
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: globalParCopy !< copy of global parameters
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: globalCorrections !< correction x_i (from A*x_i=b_i in iteration i)
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: globalParStart     !< start value for global parameters
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: globalParPreSigma  !< pre-sigma for global parameters
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: globalParPreWeight !< weight from pre-sigma
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: globalParCons !< global parameters (number of) constraints
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: globalParComments !< global parameters comments
    ! global matrix, vector
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: globalMatD !< global matrix 'A' (double, full or sparse)
    REAL(mps), DIMENSION(:), ALLOCATABLE :: globalMatF !< global matrix 'A' (float part for compressed sparse)
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: globalVector !< global vector 'x' (in A*x=b)
    INTEGER(mpl), DIMENSION(:), ALLOCATABLE :: globalRowOffsets !< row offsets for full or unpacked matrix
    INTEGER(mpl), DIMENSION(:), ALLOCATABLE :: globalCounter !< global counter (entries in 'x')
    ! AVPROD (A*x=b) by MINRES
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: vecXav !< vector x for AVPROD (A*x=b)
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: vecBav !< vector b for AVPROD (A*x=b)
    ! preconditioning
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: matPreCond !< preconditioner matrix (band and other parts)
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: indPreCond !< preconditioner pointer array
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: blockPreCond !< preconditioner (constraint) blocks
    INTEGER(mpl), DIMENSION(:), ALLOCATABLE :: offPreCond !< preconditioner (block matrix) offsets
    ! accurate sumation (per thread)
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: globalChi2SumD    !< fractional part of Chi2 sum
    INTEGER(mpl), DIMENSION(:), ALLOCATABLE :: globalChi2SumI !< integer part of Chi2 sum
    INTEGER(mpl), DIMENSION(:), ALLOCATABLE :: globalNdfSum   !< NDF sum
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: globalNdfSumW     !< weighted NDF sum
    ! auxiliary vectors
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: workspaceD !< (general) workspace (D)
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: workspaceDiag !< diagonal of global matrix (for global corr.)
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: workspaceRow !< (pivot) row of global matrix (for global corr.)
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: workspaceLinesearch !< workspace line search
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: workspaceDiagonalization !< workspace diag.
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: workspaceEigenValues !< workspace eigen values
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: workspaceEigenVectors !< workspace eigen vectors
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: workspaceI !< (general) workspace (I)
#ifdef LAPACK64
    ! LAPACK
    INTEGER(mpl):: lplwrk=1 !< length of LAPACK WORK array
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: lapackQL !< LAPACK QL (QL decomp.)
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: lapackTAU !< LAPACK TAU (QL decomp.)
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: lapackWORK !< LAPACK work array
    INTEGER(mpl), DIMENSION(:), ALLOCATABLE :: lapackIPIV !< LAPACK IPIV (pivot)
#endif
    ! constraint matrix, residuals
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: matConsProduct !< product matrix of constraints
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: vecConsResiduals !< residuals of constraints
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: vecConsSolution !< solution for constraint elimination
    ! constraint sorting, blocks
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: vecConsStart !< start of constraint in listConstraints (unsorted input)
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: matConsRanges !< parameter ranges for constraints
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: matConsSort !< keys and index for sorting
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: matConsGroups !< start of constraint groups, parameter range 
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: vecConsGroupCounts !< counter for constraint groups 
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: matConsGroupStats !< statistics for constraint groups 
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: matConsBlocks !< start of constraint blocks, parameter range 
    ! monitoring of input residuals
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: measIndex !< mapping of 1. global label to measurement index
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: measHists !< measurement histograms (100 bins per thread)
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: measRes !< average measurement error   
    ! global parameter mapping
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: globalParLabelIndex !< global parameters label, total -> var. index, group infos.
                                                                    !!
                                                                    !! 1 : global label \n
                                                                    !! 2 : index in variable parameter list (negative if fixed)\n
                                                                    !! 3 : label of par. group start \n
                                                                    !! 4 : index of par. group 
    INTEGER(mpl), DIMENSION(:), ALLOCATABLE :: globalParLabelCounter !< global parameters label counters
    INTEGER(mpl), DIMENSION(:), ALLOCATABLE :: globalParLabelZeros   !< global parameters label with zero derivative counters
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: globalParHashTable    !< global parameters hash table
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: globalParVarToTotal   !< global parameters variable -> total index
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: globalAllParToGroup   !< all parameters variable -> group index
    INTEGER(mpi), DIMENSION(-8:0) :: globalParHeader = 0 !< global parameters (mapping) header
                                                    !!
                                                    !!  0: length of labels/indices; \n
                                                    !! -1: number of stored items; \n
                                                    !! -2: =0 during build-up; \n
                                                    !! -3: next number; \n
                                                    !! -4: (largest) prime number (< length); \n
                                                    !! -5: number of overflows; \n
                                                    !! -6: nr of variable parameters; \n
                                                    !! -7: call counter for build-up;
                                                    !! -8: number of sorted items;
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: globalTotIndexGroups   !< 1. (total) index, size per group
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: globalAllIndexGroups   !< 1. (all variable) index per group

    ! row, column information for custom sparse matrix (MINRES)
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: sparseMatrixColumns   !< (compressed) list of columns for sparse matrix
    INTEGER(mpl), DIMENSION(:,:), ALLOCATABLE :: sparseMatrixOffsets !< row offsets for column list, sparse matrix elements
    ! row, column information for CSR3 sparse matrix (PARDISO)
    INTEGER(mpl), DIMENSION(:), ALLOCATABLE :: csr3RowOffsets !< row offsets for column list
    INTEGER(mpl), DIMENSION(:), ALLOCATABLE :: csr3ColumnList !< list of columns for sparse matrix
    ! read buffer
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: readBufferInfo !< buffer management (per thread)
                                                            !! 1: file index \n
                                                            !! 2: active thread number \n 
                                                            !! 3: record to read \n
                                                            !! 4: record counter  \n 
                                                            !! 5: offset \n 
                                                            !! 6: nb records cont. doubles

    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: readBufferPointer !< pointer to used buffers. Stores pointer (=offset) to the first int/float read from file for a given record index
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: readBufferDataI !< integer data. In addition to file content, we inject two header words ahead of each record:
    !! 1. local rec. number 
    !! 2. pointer to *last* entry of record) \n
    !! 
    !! Pointers in readBufferPointer point to second field of header, 
    !! before the data content
    REAL(mpr4), DIMENSION(:), ALLOCATABLE :: readBufferDataF !< float data
    REAL(mpr8), DIMENSION(:), ALLOCATABLE :: readBufferDataD !< double data. In addition to file content, we inject two header words ahead of each record:
    !! 1. file number of this record (as double...)
    !! 2. weight assigned to the input file
    !! 
    !! Pointers in readBufferPointer point to second field of header, 
    !! before the data content

    ! global parameter usage in record
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: globalIndexUsage !< indices of global par in record
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: backIndexUsage   !< list of global par in record
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: appearanceCounter !< appearance statistics for global par (first/last file,record). 5 sequential entries per par.
    !< structure: 
    !< 1) first file with record
    !< 2) local record nb within file (1) 
    !< 3) last file with record
    !< 4) local record nb within file (3) 
    !< 5) total appearance counter of parameter 

    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: pairCounter !< number of paired parameters (in equations)
    ! global parameter usage from all records
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: globalIndexRanges   !< global par ranges
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: matParBlockOffsets   !< global par block offsets (parameter, constraint blocks)
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: vecParBlockConOffsets   !< global par block (constraint) offsets 
    ! local fit
    REAL(mpd), DIMENSION(:), ALLOCATABLE::blvec  !< local fit vector 'b' (in A*x=b), replaced by 'x'
    REAL(mpd), DIMENSION(:), ALLOCATABLE::clmat  !< local fit matrix 'A' (in A*x=b)
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE:: ibandh !< local fit 'band width histogram' (band size autodetection)
    ! scratch arrays for local fit
    REAL(mpd), DIMENSION(:), ALLOCATABLE::vbnd !< local fit band part of 'A'
    REAL(mpd), DIMENSION(:), ALLOCATABLE::vbdr !< local fit border part of 'A'
    REAL(mpd), DIMENSION(:), ALLOCATABLE::aux  !< local fit 'solutions for border rows'
    REAL(mpd), DIMENSION(:), ALLOCATABLE::vbk  !< local fit 'matrix for border solution'
    REAL(mpd), DIMENSION(:), ALLOCATABLE::vzru !< local fit 'border solution'
    REAL(mpd), DIMENSION(:), ALLOCATABLE::scdiag !< local fit workspace (D)
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE:: scflag         !< local fit workspace (I)
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: localEquations !< indices (decodeNextMeasurement) for local equations (measurements)
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: localCorrections !< local fit corrections (to residuals)
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: localGlobalMatrix !< matrix correlating local and global par, content
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: localGlobalMap !< matrix correlating local and global par, map (counts)
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: localGlobalStructure !< matrix correlating local and global par, (sparsity) structure
    ! update of global matrix
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: writeBufferInfo  !< write buffer management (per thread)
    REAL(mps), DIMENSION(:,:), ALLOCATABLE :: writeBufferData     !< write buffer data (largest residual, Chi2/ndf, per thread)
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: writeBufferIndices !< write buffer for indices
    REAL(mpd), DIMENSION(:), ALLOCATABLE :: writeBufferUpdates !< write buffer for update matrices
    INTEGER(mpi), DIMENSION(-6:6) :: writeBufferHeader = 0 !< write buffer header (-6..-1: updates, 1..6: indices)
                                                      !!
                                                      !! +/-1: buffer size (words) per thread; \n
                                                      !! +/-2: min number of free words; \n
                                                      !! +/-3: number of buffer flushes; \n
                                                      !! +/-4: number of buffer overruns; \n
                                                      !! +/-5: average fill level; \n
                                                      !! +/-6: peak fill level;
    !> list items from steering file
    INTEGER(mpi) :: lenParameters=0   !< length of list of parameters from steering file
    TYPE(listItem), DIMENSION(:), ALLOCATABLE :: listParameters   !< list of parameters from steering file
    INTEGER(mpi) :: lenPresigmas=0    !< length of list of pre-sigmas from steering file
    TYPE(listItem), DIMENSION(:), ALLOCATABLE :: listPreSigmas    !< list of pre-sgmas from steering file
    INTEGER(mpi) :: lenConstraints=0  !< length of list of constraints from steering file
    TYPE(listItem), DIMENSION(:), ALLOCATABLE :: listConstraints  !< list of constraints from steering file
    INTEGER(mpi) :: numMeasurements=0 !< number of (external) measurements from steering file
    INTEGER(mpi) :: lenMeasurements=0 !< length of list of (external) measurements from steering file
    TYPE(listItem), DIMENSION(:), ALLOCATABLE :: listMeasurements !< list of (external) measurements from steering file
    INTEGER(mpi) :: lenComments=0  !< length of list of (global parameter) comments from steering file
    TYPE(listItemC), DIMENSION(:), ALLOCATABLE :: listComments  !< list of comments from steering file
#ifdef LAPACK64
#ifdef PARDISO
    INTEGER(mpi) :: lenPARDISO=0  !< length of list of Intel oneMKL PARDISO parameters (indices 1..64)
    TYPE(listItemI), DIMENSION(:), ALLOCATABLE :: listPARDISO !< list of Intel oneMKL PARDISO parameters
    INTEGER(mpl) :: ipdmem=0 !< memory (kB) used by Intel oneMKL PARDISO
#endif
#endif
    INTEGER(mpi) :: lenPostProc=0  !< length of post processing string
    CHARACTER (LEN=1024) :: cPostProc !< post processing string
    !======================================================
    ! file information
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: mfd   !< file mode: cbinary =1, text =2, fbinary=3
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: lfd   !< length of file name
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: nfd   !< index (line) in (steering) file
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: recordNbInFile !< (1,.)=  number of records in file, (2,..)= file order
    INTEGER(mpl), DIMENSION(:), ALLOCATABLE :: integratedRecordNb   !< file: integrated record numbers (=offset)
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: jfd   !< file: number of accepted records
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: dfd   !< file: ndf sum
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: maxRecPerFile   !< file: max. record size
    REAL(mps), DIMENSION(:), ALLOCATABLE :: cfd      !< file: chi2 sum
    REAL(mps), DIMENSION(:), ALLOCATABLE :: ofd      !< file: option
    REAL(mps), DIMENSION(:), ALLOCATABLE :: fileLevelWeight      !< binary file: weight
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: sfd !< offset (1,..), length (2,..) of binary file name in tfd
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: yfd   !< binary file: modification date
    CHARACTER (LEN=1024) :: filnam !< name of steering file
    INTEGER(mpi) :: nfnam  !< length of sterring file name
    CHARACTER, DIMENSION(:), ALLOCATABLE :: tfd !< file names (concatenation)
    INTEGER(mpi) :: ifile  !< current file (index)
    INTEGER(mpi) :: nfiles !< number of files
    INTEGER(mpi) :: nBinaryFiles  !< number of binary files
    INTEGER(mpi) :: nfilf  !< number of Fortran binary files
    INTEGER(mpi) :: nfilc  !< number of C binary files
    INTEGER(mpi) :: nfilw  !< number of weighted binary files
    INTEGER(mpi) :: readBufferSize=10000 !< default read buffer size (I/F words, half record length)

END MODULE mpmod
