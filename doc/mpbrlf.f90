!> \file
!! -t=BRLF test case.

! No code - only documentation for doxygen.

!> \page test_brlf_page Internal silicon tracker test case (-t=BRLF)
!!
!! \author C. Kleinwort, DESY, 2022
!!
!! \remark
!! <details>
!!   <summary>Click (on triangle):</summary>
!!   for (un)folding details.
!! </details>
!!
!! \tableofcontents
!!
!! \section test_brlf_setup Setup
!! \subsection test_brlf_setup_det Detector
!! Test case with a simple silicon strip tracker with ten equidistant, parallel planes in a beam telescope like
!! configuration without magnetic field (\ref mptest2.f90 "implementation" &rarr; More).
!! The strip sensor have a size of 2 cm in X and 4 cm in Y and a thickness corresponding to
!! 2% of a radiation length. The planes are 10 cm apart and oriented
!! perpendicular to the Z-axis (as average beam direction). Each plane (1..10) has a layer of 50 sensors
!! (0..49, 10 columns(X), 5 rows(Y)) measuring in the X-direction (with 20 &micro;m resolution).
!! Every third plane (1,4,7,10) has an additional &plusmn;5&deg; stereo layer.
!! Ten thousand tracks are simulated (with \f$\log(p)\f$ uniform in 10..100 GeV/c).
!!
!! \subsection test_brlf_setup_track Track model
!! The *broken lines* track model is used to describe the multiple scattering in the silicon sensors [\ref ref_sec "ref 4,6"].
!! For the seed trajectory the true track parameters at the first layer (before any scattering) are used.
!! Other models with different handling of multiple scattering are available too (see \ref mptest2.f90 "implementation"  &rarr; More).
!!
!! \subsection test_brlf_setup_align Alignment
!! As geometry reference the average corrections for some columns of modules in the first and last plane are
!! forced to zero (using (four) \ref ssec-lincon "linear equality constraints") to avoid correlations with
!! the four degrees of freedom provided by the track model.
!!
!! There are 700 alignment parameters in total:
!!   - 500 (all planes) in X-direction (global label: <tt>50*\#plane+\#sensor</tt>)
!!   - 200 (stereo planes) in Y-direction (global label: <tt>50*\#plane+\#sensor+1000</tt>)
!!
!! \subsection test_brlf_setup_misalign Misalignment
!! All modules (except those appearing in the constraints) are mislaligned with a RMS of 100 &micro;m in X and Y.
!!
!! \subsection test_brlf_setup_size Size of alignment problem
!! The number of global parameters (here 700) and the number of constraints (\c NCGB, here 4)
!! determine the (maximal) size of the global matrix. The actual (total) number of global parameters
!! (identified by their global labels) appearing in the binary files (\c NTGB) may be smaller.
!! In addition global parameters may have only few *entries* (number of appearances) producing possibly inaccurate results.
!! An \ref ch-entries-ext "entries" cut is available to fix parameters with too few entries.
!! Only the (\c NVGB) variable global parameters are used by pede. The number of fit parameters (\c NFGB) is
!! \c NVGB-NCGB for elimination of constraints or \c NVGB+NCGB for Lagrange multipliers.
!!
!! \section test_brlf_prep Preparations
!! 1. Download the software package from the DESY \c gitlab server to
!!    \a target directory, e.g. (shallow clone):
!!
!!         git clone --depth 1 --branch V04-12-01 \
!!             https://gitlab.desy.de/claus.kleinwort/millepede-ii.git target
!!
!! 2. Create **Pede** executable (in \a target directory):
!!
!!         make pede
!!
!! 3. Optionally cleanup:
!!
!!         rm mp2???.???
!!
!! \section test_brlf_create Create pede files
!!
!!
!!         ./pede -t=BRLF > pede.dump
!!
!! 1. Creates the pede input binary (<tt>mp2tst.bin</tt>) and \ref ssec-textfiles "text files":
!!
!!         cksum mp2???.???
!!         2762434670      448 mp2con.txt
!!         3457733667     1971 mp2str.txt
!!         2736700272 14919648 mp2tst.bin
!!
!!    - <details><summary><a href="../examples/brlf/run0/mp2con.txt">constraints</a> file
!!      <tt>mp2con.txt</tt></summary><pre>
!! Constraint  0.0
!!      54   1.00000
!!      64   1.00000
!!      74   1.00000
!!      84   1.00000
!!      94   1.00000
!! Constraint  0.0
!!    1054   1.00000
!!    1064   1.00000
!!    1074   1.00000
!!    1084   1.00000
!!    1094   1.00000
!! Constraint  0.0
!!     504   1.00000
!!     514   1.00000
!!     524   1.00000
!!     534   1.00000
!!     544   1.00000
!! Constraint  0.0
!!    1504   1.00000
!!    1514   1.00000
!!    1524   1.00000
!!    1534   1.00000
!!    1544   1.00000
!!      </pre></details>
!!    - <details><summary><a href="../examples/brlf/run0/mp2str.txt">steering</a> file
!!      <tt>mp2str.txt</tt></summary><pre>
!!*        *** Default test steering file ***
!!
!!*            Additinal files
!!fortranfiles ! following bin files are fortran
!!mp2con.txt   ! constraints text file
!!mp2tst.bin   ! binary data file
!!Cfiles       ! following bin files are Cfiles
!!
!!*            Selection of variable global parameters
!!*entries  10 ! lower   limit on number of entries/parameter
!!*entries  25 ! default limit on number of entries/parameter
!!*entries  50 ! higher  limit on number of entries/parameter
!!
!!*            Initial parameter values, presigmas
!!*Parameter        ! define parameter attributes (start  of list)
!!*205  -0.01   0.  ! start value
!!*206   0.01   0.  ! start value
!!*215   0.01  -1.  ! fix parameter at value
!!*216   0.    -1.  ! fix parameter at value
!!
!!*            Handling of outliers, tails etc
!!*hugecut 50.0            !cut factor in iteration 0
!!*chisqcut 30.0 6.0       ! cut factor in iterations 1 and 2
!!*outlierdownweighting  2 ! number of internal iterations (> 1)
!!*dwfractioncut      0.2  ! 0 < value < 0.5
!!*presigma           0.01 ! default value for presigma
!!*regularisation 1.0      ! regularisation factor
!!*regularisation 1.0 0.01 ! regularisation factor, pre-sigma
!!
!!*            Solution methods
!!method fullMINRES       3 0.01 ! (approximate) MINRES, full storage
!!method sparseMINRES     3 0.01 ! (approximate) MINRES, sparse storage
!!*mrestol      1.0D-8           ! epsilon for MINRES convergence
!!*bandwidth 0                   ! width of MINRES precond. band matrix
!!method diagonalization 3 0.001 ! diagonalization (-> millepede.eve)
!!method decomposition   3 0.001 ! Cholesky decomposition
!!method inversion       3 0.001 ! Gauss matrix inversion
!!* last method is applied
!!
!!*            Additional output. monitoring
!!printcounts           ! print number of entries
!!*monitorresiduals     ! poor man's DMR (-> millepede.mon)
!!*printrecord   1  2   ! debug printout for records
!!*printrecord  -1 -1   ! debug printout for bad data records
!!
!!end ! optional for end-of-data
!!      </pre></details>
!!
!!    <details><summary>
!!    With the python2 script readMilleBinary.py binary files can be read and \ref mpmille "records" printed in text form:
!!    </summary>
!!    The *broken lines* track model defines at each (scattering) layer offsets in
!!    X (odd local labels) and Y-direction (even local labels)
!!    as fit parameeters. The measurement in each layer depends ("local derivatives") on a
!!    linear combination of the offsets in that layer depending on the orientation in the XY-plane
!!    (0.&deg;, &plusmn;5&deg;).
!!    <pre>
!!Detected Fortran binary file
!! === NR  1 185
!! -g- meas.  1 73 1 1 -0.00794364139438 0.00200000009499
!! local   [1]
!! local   [1.0]
!! global  [73]
!! global  [1.0]
!! -g- meas.  2 73 2 2 -0.0057476432994 0.00200000009499
!! local   [3, 4]
!! local   [0.9961847066879272, 0.08726999908685684]
!! global  [73, 1073]
!! global  [0.9961847066879272, 0.08726999908685684]
!! ...
!! -l- meas.  15 1 3 0 0.0 2.02223127417e-05
!! local   [1, 3, 5]
!! local   [2.0, -2.1052632331848145, 0.10526315867900848]
!! -l- meas.  16 2 3 0 0.0 2.02223127417e-05
!! local   [2, 4, 6]
!! local   [2.0, -2.1052632331848145, 0.10526315867900848]
!! ...
!!    </pre>
!!    Kinks constructed from three adjacent layers are describing the multiple scattering in the
!!    middle layer and are connecting the layers (in the XZ and YZ planes).
!!    The local derivatives depend on the distances between the layers.
!!    </details>
!!
!!         python2 tools/readMilleBinary.py mp2tst.bin -n 3 > mp2binary.txt
!!
!!    The first three tracks are printed into <a href="../examples/brlf/mp2binary.txt"><tt>mp2binary.txt</tt></a>.
!! 2. Runs pede and creates output files (<tt>millepede.???</tt>):
!!    + <a href="../examples/brlf/run0/.pede.dump.txt"><tt>pede.dump</tt></a>:
!!      Redirected standard output, large overlap with <tt>millepede.log</tt>.
!!      Information from the preparation and execution of the \ref ssec-globalfit "global fit". E.g.:
!!      - <details><summary>Pede start: Version, configuration and start time.</summary><pre>
!! (\$Id: cc523e423547ec1e335b86a2ead55e4b2e6b3040 \$)
!! using OpenMP (TM)
!! compiled with gcc 12.2.0
!!
!!   <  Millepede II-P starting ... Tue Nov 22 10:20:58 2022
!!                                  HOSTNAME
!!        </pre></details>
!!
!!               head -8 pede.dump
!!
!!      - <details><summary>End of preparation: Number of global parameters (total, variable, fit),
!!        constraints (number, rank, ..), solution method, storage mode and size of global matrix.
!!        </summary><pre>
!!     NTGB =       700 = total number of parameters
!!                      (all parameters, appearing in binary files)
!!     NVGB =       683 = number of variable parameters
!!                      (appearing in fit matrix/vectors)
!!     NAGB =       683 = number of all parameters
!!                      (including Lagrange multiplier or reduced)
!!   NTPGRP =       700 = total number of parameter groups
!!   NVPGRP =       683 = number of variable parameter groups
!!     NFGB =       679 = number of fit parameters
!!     NOFF =    232903 = max number of off-diagonal elements
!!     NCGB =         4 = number of constraints
!!    NAGBN =        18 = max number of global parameters in an event
!!    NALCN =        28 = max number of local parameters in an event
!!    NAEQN =        38 = max number of equations in an event
!!   NCACHE =  25000000 = number of words for caching
!!                      cache splitting   56.7 %   3.0 %  40.2 %
!!
!!
!! Solution method and matrix-storage mode:
!!      METSOL = 1:  matrix inversion
!!                   with           3  iterations
!!      MATSTO = 1:  full symmetric matrix, (n*n+n)/2 elements
!! Convergence assumed, if expected dF <  0.1000E-02
!! Constraints handled by elimination
!!
!! Rank of product matrix of constraints is           4  for           4  constraint equations
!!
!! QL decomposition of constraints matrix
!!    largest  |eigenvalue| of L:   -2.2360679774997898
!!    smallest |eigenvalue| of L:   -2.2360679774997898
!!
!! Size of global matrix:           1  MB
!!
!! _______________________________LOOP2-end__________________________________
!!        </pre></details>
!!
!!               grep -B33 "LOOP2-end" pede.dump
!!
!!        For sparse storage information about the sparsity structure is added (use <tt>-B48</tt>).
!!      - <details><summary>End of global fit:
!!        Quality: final Sum(Chi^2)/Sum(Ndf) (should be close to 1.),
!!        number of rejected records/tracks (if any, should be <1%).
!!        </summary><pre>
!!                                                                       00000   S
!!  2  3 0.10259E+06 0.28E-14                0  0           0.0       0  00000 F
!!
!! Sum(Chi^2)/Sum(Ndf) =   102589.35168164221
!!                     / (       99992 -         679 )
!!                     =   1.0329901592101960
!!
!!
!! _____________________________Iteration-end________________________________
!!        </pre></details>
!!
!!               grep -A5 -B3 "/Sum(Ndf)" pede.dump
!!
!!      - <details><summary>Pede end: End time and maximum memory allocation.
!!        </summary><pre>
!!
!!   <  Millepede II-P ending   ... Tue Nov 22 10:20:58 2022
!!
!!      Peak dynamic memory allocation:    0.102192 GB
!!
!!        </pre></details>
!!
!!               tail -5 pede.dump
!!
!!      - In this special case (creation of input files) the simulated and fitted
!!        <a href="../examples/brlf/run0/pede.dump.special.txt">(mis)alignment parameters</a>
!!        (and their differences) are printed close to the end of the file.
!!
!!               grep -A708 "Misalignment test Si tracker" pede.dump
!!
!!        From this a steering file (fragment)
!!        <a href="../examples/brlf/mp2param.txt"><tt>mp2param.txt</tt></a>
!!        defining the simulated values
!!        as start values ("starting from truth") can be created, e.g.:
!!
!!               grep -A708 ... | awk 'BEGIN {print "Parameter"} NR>9 {print $2,$3,0.}' > mp2param.txt
!!
!!    + <a href="../examples/brlf/run0/.millepede.end.txt"><tt>millepede.end</tt></a>:
!!      Single line text file with \ref exit_code_page "exit code" and message.
!!      Be careful with results for exit code > 1.
!!    + <a href="../examples/brlf/run0/.millepede.log.txt"><tt>millepede.log</tt></a>:
!!      Log file, large overlap with <tt>pede.dump</tt>.
!!    + <a href="../examples/brlf/run0/.millepede.his.txt"><tt>millepede.his</tt></a>:
!!      Internal histogram file. \ref readPedeHists.C "Script" for viewing with ROOT available.
!!    + <a href="../examples/brlf/run0/.millepede.res.txt"><tt>millepede.res</tt></a>:
!!      Result file with alignment corrections from global fit.
!!      <details><summary>Contains one line for each global parameter with three
!!      \ref sssec-parinf "input values" and for variable parameters up to three output values:
!!      </summary><pre>
!! Parameter   ! first 3 elements per line are significant (if used as input)
!!        50     0.12387E-01    0.0000       0.12387E-01   0.38268E-03         185
!!        51    -0.10100E-02    0.0000      -0.10100E-02   0.31170E-03         350
!!        52     0.35959E-02    0.0000       0.35959E-02   0.28972E-03         366
!!        53     0.10121E-01    0.0000       0.10121E-01   0.28810E-03         331
!! ...
!!       150      0.0000        0.0000
!! ...
!!      </pre></details>
!!      - Label.
!!      - Total value (start value (usually 0.), updated from fit).
!!      - Pre-sigma (usually zero or negative to fix parameter).
!!      - Correction value from fit.
!!      - Error of correction for some \ref ch-methods "solution methods" (e.g. \ref ch-inv "inversion").
!!      - Number of entries (only with option \ref cmd-printcounts).
!!
!!    (\ref sssec-outfiles in draft manual.)
!!
!! Save steering and output files:
!!
!!        mkdir run0; cp mp2???.txt run0; cp pede.dump run0; cp millepede.??? run0
!!
!! \section test_brlf_exercises Exercises
!! \subsection test_brlf_check_input 1. Check input
!! Run pede in special mode to \ref ch-checkinput "check input" only (no solution calculated):
!!
!!        ./pede -C mp2str.txt > pede.dump
!!
!! The standard output in <a href="../examples/brlf/run1/.pede.dump.txt"><tt>pede.dump</tt></a>
!! contains now a lot of information about the (linear equality) constraints:
!! - <details><summary>After \"\c LOOP2\" start: For each constraint (in order of appearance in steering files)
!!   the number of all and of variable global parameters involved.
!!   </summary><pre>
!! _________________________________LOOP2____________________________________
!!
!! constraint     1 :         5 parameters,        5 variable
!! constraint     2 :         5 parameters,        5 variable
!! constraint     3 :         5 parameters,        5 variable
!! constraint     4 :         5 parameters,        5 variable
!! PRPCON:           4  constraints accepted
!!   </pre>
!!   **Empty** constraints (without variable global parameters) will cause singularities.
!!   </details>
!!
!!        grep -A6 '_LOOP2_' pede.dump
!!
!! - <details><summary>Around \"<tt>PRPCON: constraints split</tt>\": The constraints
!!   are sorted by (global) label range, split into disjoint groups and
!!   combined into non overlapping blocks.
!!   </summary>
!!   Definition of **sorted** constraints (with reference to steering file
!!   (index, first line number)), **groups** (with first/last constraint)
!!   and **blocks** (with first/last group) and label ranges:
!!   <pre>
!!  Cons. sorted       index file: index, first line first label  last label
!!  Cons. group        index first cons.  last cons. first label  last label
!!  Cons. block        index first group  last group first label  last label
!!  Cons. sorted           1           1           1          54          94
!!  Cons. group            1           1           1          54          94
!!  Cons. block            1           1           1          54          94
!!  Cons. sorted           2           3          13         504         544
!!  Cons. group            2           2           2         504         544
!!  Cons. block            2           2           2         504         544
!!  Cons. sorted           3           2           7        1054        1094
!!  Cons. group            3           3           3        1054        1094
!!  Cons. block            3           3           3        1054        1094
!!  Cons. sorted           4           4          19        1504        1544
!!  Cons. group            4           4           4        1504        1544
!!  Cons. block            4           4           4        1504        1544
!!
!! PRPCON: constraints split into            4 (disjoint) groups,
!!         groups  combined  into            4 (non overlapping) blocks
!!         max group size (cons., par.)            1          41
!!         max block size (cons., par.)            1          41
!!         total block matrix sizes                       164                    4
!!   </pre></details>
!!
!!        grep -A4 -B16 'PRPCON: constraints' pede.dump
!!
!! - <details><summary>Before \"<tt>Rank of product</tt>\": For each constraint group the size
!!   and the rank of the product matrix (\f$\Vek{A} \Vek{A}\trans\f$
!!   in equation \ref eq-parsol "(11)") to verify linear independence.
!!   </summary><pre>
!!  Constraint group, \#con, rank           1           1           1
!!  Constraint group, \#con, rank           2           1           1
!!  Constraint group, \#con, rank           3           1           1
!!  Constraint group, \#con, rank           4           1           1
!!
!! Rank of product matrix of constraints is           4  for           4  constraint equations
!!   </pre></details>
!!
!!        grep  -B6 'Rank of product' pede.dump
!!
!! The results file <a href="../examples/brlf/run1/.millepede.res.txt"><tt>millepede.res</tt></a>
!! contains now several different sections:
!! - <details><summary>**Global parameter entries**. One line for each global parameter with the
!!   three \ref sssec-parinf "input values", the number of entries from the binary files, the
!!   index of the constraint group the parameter is part of (0 for none) and the parameter status
!!   (variable or fixed).
!!   </summary>
!!   For the first parameter of a parameter group (consecutive parameter (labels) appearing
!!   always togther) the line starts with '!>'.
!!   <pre>
!! ! === global parameters ===
!! ! fixed-1: by pre-sigma, -2: by entries cut, -3: by iterated entries cut
!! !      Label       Value     Pre-sigma         Entries Cons. group  Status
!! !>        50      0.0000        0.0000             185           0  variable
!! !>        51      0.0000        0.0000             350           0  variable
!! !>        52      0.0000        0.0000             366           0  variable
!! !>        53      0.0000        0.0000             331           0  variable
!! !>        54      0.0000        0.0000             357           1  variable
!! !>        55      0.0000        0.0000             373           0  variable
!! !>        56      0.0000        0.0000             376           0  variable
!! ...
!! !>       150      0.0000        0.0000              18           0  fixed-2
!! ...
!!   </pre></details>
!!
!!        grep -e"\! " -e "\!>" millepede.res
!!
!! - <details><summary>**Global parameter appearance**. One line for each global parameter
!!   with the file and record number of the first and last appearance in the binary
!!   files, the number of binary files with appearance and the number of
!!   paired global parameters (appearing together in measurements).
!!   </summary><pre>
!! !.
!! !.Appearance statistics
!! !.     Label  First file and record  Last file and record   \#files  \#paired-par
!! !.        50          1        169          1       9931          1          1
!! !.        51          1         63          1       9976          1          1
!! !.        52          1         15          1       9914          1          1
!! !.        53          1         64          1       9992          1          1
!! !.        54          1        132          1       9968          1          1
!! !.        55          1         22          1       9984          1          1
!! !.        56          1          2          1       9984          1          1
!! ...
!! !.       100          1        550          1       9931          1          0
!! !.       101          1         63          1       9925          1          0
!! !.       102          1         58          1       9976          1          0
!! ...
!!   </pre>
!!   The first plane (labels 50..99) contains a stereo layer. There each alignment
!!   offset in X or Y is paired with one in the other (Y or X) direction.
!!   </details>
!!
!!        grep -e"\!\." millepede.res
!!
!! - <details><summary>**Constraint group entries**. One line for each constraint
!!   group (described by contributing (variable global) parameters)
!!   with number of constraints, number of entries and label range.
!!   Additional lines with paired label ranges.
!!   </summary><pre>
!! * === constraint groups ===
!! *  Group  \#Cons.     Entries First label  Last label   Paired label range
!! *      1       1        2231          54          94
!! *:                                                        1054 ..        1054
!! *:                                                        1064 ..        1064
!! *:                                                        1074 ..        1074
!! *:                                                        1084 ..        1084
!! *:                                                        1094 ..        1094
!! *      2       1        2049         504         544
!! *:                                                        1504 ..        1504
!! *:                                                        1514 ..        1514
!! *:                                                        1524 ..        1524
!! *:                                                        1534 ..        1534
!! *:                                                        1544 ..        1544
!! ...
!!   </pre></details>
!!
!!        grep -e"* " -e"*:" millepede.res
!!
!! - <details><summary>**Constraint group appearance**. One line for each constraint
!!   group (described by contributing (variable global) parameters)
!!   with the file and record number of the first and last appearance in the binary
!!   files and the number of binary files with appearance.
!!   </summary><pre>
!! *.
!! *.Appearance statistics
!! *.     Group  First file and record  Last file and record   #files
!! *.         1          1         13          1       9998          1
!! *.         2          1          5          1       9963          1
!! *.         3          1         13          1       9998          1
!! *.         4          1          5          1       9963          1
!!   </pre></details>
!!
!!        grep -e"*\." millepede.res
!!
!!
!! Save steering and output files:
!!
!!        mkdir run1; cp mp2???.txt run1; cp pede.dump run1; cp millepede.??? run1
!!
!! \subsection test_brlf_check_reject 2. Reject bad tracks
!! Depending on the quality of the \ref par-locfitv "local (track) fit" very bad
!! tracks are \ref localfit-rejection "rejected"
!! (\f$ \chi^2 > f_i\cdot \chi^2_\mathrm{cut} \f$) automatically and bad tracks can be
!! removed by the user with the \"<tt>\ref cmd-chisqcut "chisqcut f1 f2"</tt>\" option. The
!! scaling factor \f$f_i\f$ should start with a large value \f$f_1\f$ to allow for initial
!! misalignment and be further reduced (to 1.) in the internal \ref par-iter "iterations"
!! of the \ref sssec-glofit "global fit".
!!
!! \subsubsection test_brlf_modify_steering Modify steering:
!! <details><summary>Uncomment \"<tt>chisqcut</tt>\" option in the steering file.
!! </summary><pre>
!! chisqcut 30.0 6.0        ! cut factor in iterations 1 and 2
!! </pre></details>
!!
!! \subsubsection test_brlf_rerun_pede Rerun pede:
!!
!!        ./pede mp2str.txt > pede.dump
!!
!! \subsubsection test_brlf_check_output Check output:
!! <details><summary>Watch the pede exit code in <tt>millepede.end</tt>.
!! </summary><pre>
!!     2   Ended with severe warnings (insufficient measurements)
!! </pre></details>
!!
!! The severe warnings indicate a major problem. More details are
!! in <a href="../examples/brlf/run2/.pede.dump.txt"><tt>pede.dump</tt></a>.
!! - <details><summary> During iterations:
!!   </summary><pre>
!!  3  4 0.10086E+06 0.17E+02                0  0           1.6      33  00000 F
!!                                                                       00000   S
!!
!!  ... warning ...
!!  global parameters with too few (< MREQENA) accepted entries:            2
!!  minimum entries:            7  for label          209
!!
!!  4  5 0.10013E+06 0.53E+01 0.781 0.062    0  0 -1 1.000  1.0      79  00000 F
!! ...
!!   </pre>
!!   After rejection of bad tracks there are only 7 entries left for the global
!!   parameter with the label "209". This is less than the minumum required by
!!   <tt>MREQENA</tt> (default 10, second parameter of \c entries cut).
!!   </details>
!!
!!        grep -C3 '... warning ...' pede.dump
!!
!! - <details><summary> After iterations:
!!   </summary><pre>
!!       ngWarningWarningWarningWarningWarningWarningWarningWarningW
!!       gWarningWarningWarningWarningWarningWarningWarningWarningWa
!!
!!         Possible bad elements =           6  in global vector
!!         (too few accepted entries)
!!         (indicated in millepede.res by counts<0)
!!         => please check mille data and ENTRIES cut
!!
!!       WarningWarningWarningWarningWarningWarningWarningWarningWar
!!   </pre>
!!   During the iterations six times a global parameter had too few accepted entries.
!!   All those parameters are flagged in \c millepede.res by an negative
!!   count value (-(\#entries+1)).
!!   </details>
!!
!!        grep -C4 'too few accepted' pede.dump
!!
!!
!! Get affected global parameter from result file
!! <a href="../examples/brlf/run2/.millepede.res.txt"><tt>millepede.res</tt></a>.
!! <details><summary>Select lines with negative count value, e.g.:
!! </summary><pre>
!!       209     0.76209E-02    0.0000       0.76209E-02   0.54649E-03          -8
!!       249    -0.11705E-01    0.0000      -0.11705E-01   0.51044E-03          -9
!! </pre>
!! The global parameters with labels 209/249 have finally less than 10 entries.
!! </details>
!!
!!        awk 'NF==6 && $6<0' millepede.res
!!
!! \subsubsection test_brlf_recheck_input Recheck input:
!! - <details><summary>Get entries in binary files for affected global parameters
!!   from \ref test_brlf_check_input "check input" results, e.g.:
!!   </summary><pre>
!! !>       209      0.0000        0.0000              35           0  variable
!! !>       249      0.0000        0.0000              40           0  variable
!!   </pre>
!!   The global parameters with labels 209/249 have initialy more than 25 entries.
!!   </details>
!!
!!        awk '($1=="!>" || $1=="! ") && ($2==209 || $2==249)' run1/millepede.res
!!
!! - <details><summary>Get number of paired parameters for affected global parameters
!!   from \ref test_brlf_check_input "check input" results, e.g.:
!!   </summary><pre>
!! !.       209          1         39          1       9867          1          1
!! !.       249          1        637          1       9547          1          1
!!   </pre>
!!   The global parameters with labels 209/249 describe offsets in a stereo layer.
!!   Therefore the X and Y offsets for each sensor are paired (appearing together
!!   for each measurement). The number of paired parameters (last column) is 1.
!!   </details>
!!
!!        awk '$1=="!." && ($2==209 || $2==249)' run1/millepede.res
!!
!! \subsubsection test_brlf_resolve_options Options to resolve severe warnings
!! Options to resolve the severe warnings for problematic labels 209/249:
!! 1. Increase the statistics (more/larger binary files) to get the minimal
!!    number of accepted entries (7) above the
!!    \ref ch-entries-ext "entries" cut value **mreqena** (10).
!! 2. Decrease the \ref ch-entries-ext "entries" cut value **mreqena** (10) below
!!    the minimal number of accepted entries (7). It must still be larger than
!!    the maximal number of paired parameters (1).
!! 3. Increase the \ref ch-entries-ext "entries" cut value **mreqenf** (25) above
!!    the maximum number of entries (40) to fix those parameters.
!! 4. Fix those parameters in the steering file:
!!
!!            Parameter
!!            209  0.  -1.
!!            249  0.  -1.
!!
!!
!! Save steering and output files:
!!
!!        mkdir run2; cp mp2???.txt run2; cp pede.dump run2; cp millepede.??? run2
!!
!! \subsection test_brlf_check_reject2_acc 3. Reject bad tracks with decreased (accepted) entries cut (Option 2)
!!
!! \subsubsection test_brlf_modify_steering2_acc Modify steering further:
!! <details><summary>Uncomment \"<tt>entries 25</tt>\" option in the steering file and
!! add **mreqena** with value 5.
!! </summary><pre>
!! entries  25 5 ! default  limit on number of entries/parameter, reduced accepted entries
!! </pre></details>
!!
!! \subsubsection test_brlf_rerun_pede2_acc Rerun pede:
!!
!!        ./pede mp2str.txt > pede.dump
!!
!! \subsubsection test_brlf_check_output2_acc Check output:
!! <details><summary>Watch the pede exit code in <tt>millepede.end</tt>.
!! </summary><pre>
!!     0   Ended normally
!! </pre></details>
!!
!! No warnings are reported. More details are
!! in <a href="../examples/brlf/run3/.pede.dump.txt"><tt>pede.dump</tt></a>.
!!
!! <details><summary>Quality of global fit.
!! </summary><pre>
!! Data rejected in last loop:
!!               0  (rank deficit/NaN)            0  (Ndf=0)              0  (huge)             79  (large)
!!
!! Sum(Chi^2)/Sum(Ndf) =   100125.94658222546
!!                     / (       99992 -         679 )
!!                     =   1.0081857015921929
!!
!!
!! _____________________________Iteration-end________________________________
!! </pre>
!! The final Sum(Chi^2)/Sum(Ndf) is very close to 1, the rejection rate
!! less than 1% (79/10000) and 679 (of maximal 696) parameters have been fitted.
!! </details>
!!
!!        grep -A5 -B3 "/Sum(Ndf)" pede.dump
!!
!!
!! The quality is fine.
!!
!! Save steering and output files:
!!
!!        mkdir run3; cp mp2???.txt run3; cp pede.dump run3; cp millepede.??? run3
!!
!! \subsection test_brlf_check_reject2 4. Reject bad tracks with increased entries cut (Option 3)
!!
!! \subsubsection test_brlf_modify_steering2 Modify steering further:
!! <details><summary>Comment \"<tt>entries 25</tt>\" again and 
!! uncomment \"<tt>entries 50</tt>\" option in the steering file.
!! </summary><pre>
!! *entries  25 5 ! default  limit on number of entries/parameter, reduced accepted entries
!! entries  50 ! higher  limit on number of entries/parameter
!! </pre></details>
!!
!! \subsubsection test_brlf_rerun_pede2 Rerun pede:
!!
!!        ./pede mp2str.txt > pede.dump
!!
!! \subsubsection test_brlf_check_output2 Check output:
!! <details><summary>Watch the pede exit code in <tt>millepede.end</tt>.
!! </summary><pre>
!!     1   Ended with warnings (bad measurements)
!! </pre></details>
!!
!! The warnings indicate some problem. More details are
!! in <a href="../examples/brlf/run4/.pede.dump.txt"><tt>pede.dump</tt></a>.
!! - <details><summary> Quality of global fit:
!!   </summary><pre>
!! Data rejected in last loop:
!!               0  (rank deficit/NaN)            0  (Ndf=0)              0  (huge)            276  (large)
!!
!! Sum(Chi^2)/Sum(Ndf) =   104495.47920560479
!!                     / (       99992 -         645 )
!!                     =   1.0518231975359578
!!
!!
!!
!!   </pre>
!!   The final Sum(Chi^2)/Sum(Ndf) is close to 1, but the rejection rate
!!   is allmost 3% (276/10000) and only 645 (of maximal 696) parameters have been fitted.
!!   </details>
!!
!!        grep -A5 -B3 "/Sum(Ndf)" pede.dump
!!
!! - <details><summary> After iterations:
!!   </summary><pre>
!!       ngWarningWarningWarningWarningWarningWarningWarningWarningW
!!       gWarningWarningWarningWarningWarningWarningWarningWarningWa
!!
!!         Fraction of rejects =   2.76 %  (should be far below 1 %)
!!         => please provide correct mille data
!!
!!       WarningWarningWarningWarningWarningWarningWarningWarningWar
!!   </pre></details>
!!
!!        grep -C3 'Fraction of rejects' pede.dump
!!
!!
!! The result is not optimal. Too many globals parameters have been fixed
!! at the **misaligned** values spoiling the local/track fit.
!!
!! Save steering and output files:
!!
!!        mkdir run4; cp mp2???.txt run4; cp pede.dump run4; cp millepede.??? run4
!!
!! \subsection test_brlf_check_reject2_truth 5. Reject bad tracks with increased entries cut (from truth)
!!
!! \subsubsection test_brlf_modify_steering2_truth Modify steering further:
!! <details><summary>Add text file \"\c mp2param.txt\" with start values to the steering file
!! (e.g. after constraints file).
!! </summary><pre>
!! mp2con.txt   ! constraints text file
!! mp2param.txt ! start values (from truth)
!! mp2tst.bin   ! binary data file
!! </pre></details>
!!
!! This defines a different starting geometry (in this case the "truth").
!!
!! \subsubsection test_brlf_rerun_pede2_truth Rerun pede:
!!
!!        ./pede mp2str.txt > pede.dump
!!
!! \subsubsection test_brlf_check_output2_truth Check output:
!! <details><summary>Watch the pede exit code in <tt>millepede.end</tt>.
!! </summary><pre>
!!     0   Ended normally
!! </pre></details>
!!
!! No warnings are reported. More details are
!! in <a href="../examples/brlf/run5/.pede.dump.txt"><tt>pede.dump</tt></a>.
!!
!! <details><summary>Quality of global fit.
!! </summary><pre>
!! Data rejected in last loop:
!!               0  (rank deficit/NaN)            0  (Ndf=0)              0  (huge)             24  (large)
!!
!! Sum(Chi^2)/Sum(Ndf) =   99050.817948768963
!!                     / (       99992 -         645 )
!!                     =  0.99701871167492695
!!
!!
!! _____________________________Iteration-end________________________________
!! </pre>
!! The final Sum(Chi^2)/Sum(Ndf) is very close to 1 and the rejection rate less
!! than 1% (24/10000) despite only 645 (of maximal 696) parameters have been fitted.
!! </details>
!!
!!        grep -A5 -B3 "/Sum(Ndf)" pede.dump
!!
!!
!! The quality is fine. Starting values matter for fixed parameters.
!!
!! Save steering and output files:
!!
!!        mkdir run5; cp mp2???.txt run5; cp pede.dump run5; cp millepede.??? run5
!!
!! \subsection test_brlf_compare 6. Compare millepede results
!! With the python2 script compareResults.py two millepede result text files
!! (\c millepede.res) can be compared. As only the first three columns
!! are evaluated a steering file (fragment) defining start values can be used too, e.g:
!!
!!        python2 compareResults.py mp2param.txt run4/millepede.res
!!
!! <details><summary>This compares the results from example 4 with the "truth".
!! </summary><pre>
!! output file         mp2compare.txt
!! total parameters    700
!! matched parameters  700
!! mean1, mean2        -0.000189471428571 -0.000542676439286
!! rms1, rms2          0.00994723953843 0.00994581992618
!! correlation         0.945525591556
!! </pre>
!! For both sets of parameters the mean and RMS and the correlation of the sets is calculated.
!! </details>
!!
!! A text file
!! (default <a href="../examples/brlf/mp2compare.txt"><tt>mp2compare.txt</tt></a>)
!! is created.
!! <details><summary>
!! It contains the parameter labels and values and can be directly used
!! as ROOT input (with \"<tt>tree->ReadFile()</tt>\"):
!! </summary><pre>
!! label/I:value_1/F:value_2/F
!!        50      0.01339     0.012557
!!        51     -0.00037  -0.00082443
!!        52      0.00403    0.0037749
!!        53      0.01061     0.010278
!!        54           -0   0.00020854
!!        55     -0.00712   -0.0075308
!! ...
!!      1545     -0.01218    -0.010553
!!      1546      -0.0029    0.0032485
!!      1547     -0.01605   -0.0082543
!!      1548     -0.01081   -0.0027929
!!      1549      0.00236     0.008585
!! </pre>
!!    For all labels appearing in boths files ("matched"): label value1 value2.
!!    For all labels appearing only in one file: -label, value, value
!! </details>
!!
!! \section test_brlf_other Other examples
!!
!! Other examples
!! (<a href="https://www.desy.de/~kleinwrt/GBL/doc/cpp/html/">C++</a>,
!! <a href="https://www.desy.de/~kleinwrt/GBL/doc/python/html/">Python2</a> or
!! <a href="https://www.desy.de/~kleinwrt/GBL/doc/python3/html/">Python3</a>)
!! are included in the
!! <a href="https://gitlab.desy.de/claus.kleinwort/general-broken-lines/-/wikis/home">
!! GeneralBrokenLines package</a> to explore track fitting and track based alignment with MP2.
!! Those will produce binary files of \ref sssec-fileinf "type" <tt>Cfiles</tt>.
