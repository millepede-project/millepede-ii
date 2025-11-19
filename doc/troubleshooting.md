 \page troubleshooting_page Troubleshooting

 \section ch-experience General experience:

 In the ideal case the track and geometry models are correct and complete,
 the measurement errors are gaussian and described perfectly
 and all dependencies are linear. With this (probably unrealistic)
 requirements pede should end *normally* (exit code 0 
 (in <tt>millepede.end</tt>)).

 In the real world measurement errors can be non gaussian (e.g.
 binary readout or multiple scattering (tails)) or dependencies
 non linear (e.g. curvature or rotations). Therefore pede usually
 ends with *warnings* (exit code 1).

 In case of problems with the linear equation system (e.g. too few measurements,
 linear dependent constraints or global matrix not positive definite) pede
 still tries to get a solution and ends with *severe warnings*
 (exit code 2). The solution is very probable not reliable and should be
 'handled with care'. The pede input (steering and binary files) should be
 \ref ch-checkinput "checked" thoroughly.

 \section ch-aborts Aborts
 In case of severe (usually technical) problems pede aborts (exit code >= 10).
 It may be helpful to \ref ch-checkinput "check" the pede input
 (steering and binary files).

 For 'too many rejects' (exit code 26) usually the measurement errors
 are largely underestimated. It may help to inflate those in pede
 (command [scaleerrors](option_page.html#scaleerrors)) or to skip the internal iterations
 (command [subito](option_page.html#subito) or command line option [-s](#-s))
 to get *some* solution.

 \section ch-segfault Segmentation fault
 In case of segmentation faults the OpenMP stack size may be to small (default 8M?).
 Increased by defining the shell variable <tt>OMP_STACKSIZE</tt>:

        export OMP_STACKSIZE=32M

 For further debugging one should add the compiler options
 \"<tt>-g -fcheck=all -fbacktrace</tt>\"
 to the \c F_FLAGS in the \c Makefile and recompile.
