
 \file
 Millepede II program, documentation.

 \author Volker Blobel, University Hamburg, 2005-2009 (initial Fortran77 version)
 \author Gero Flucke, University Hamburg (support of C-type binary files)
 \author Claus Kleinwort, DESY (maintenance and developement)
 \author Maximilian Golbirsch-Kolb, DESY (maintenance and developement)

 \copyright
 Copyright (c) 2009 - 2025 Deutsches Elektronen-Synchroton,
 Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Library General Public License as
 published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version. \n\n
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Library General Public License for more details. \n\n
 You should have received a copy of the GNU Library General Public
 License along with this program (see the file COPYING.LIB for more
 details); if not, write to the Free Software Foundation, Inc.,
 675 Mass Ave, Cambridge, MA 02139, USA.


 \mainpage Overview

# Introduction
 In certain least squares fit problems with a very large number of parameters
 the set of parameters can be divided into two classes, global and local parameters.
 Local parameters are those parameters which are present only in subsets of the
 data. Detector alignment and calibration based on track fits is one of the problems,
 where the interest is only in optimal values of the global parameters, the
 alignment parameters. The method, called Millepede, to solve the linear least
 squares problem with a simultaneous fit of all global and local parameters,
 irrespectively of the number of local parameters, is described in the draft manual.
 (Correlated measurements need to be transformed into independent measurements
 by diagonalization of their covariance matrix.)

 The Millepede method and the initial implementation has been
 developed by [V. Blobel](http://www.desy.de/~blobel) from he University of Hamburg.
 Meanwhile the code is maintained at DESY by the statistics tools group of the
 analysis center of the Helmholtz [Terascale](https://terascale.de) alliance
 using [GitLab](https://about.gitlab.com) ([code and wiki](https://gitlab.desy.de)).

 The Millepede II software is provided by DESY under the terms of the
 [LGPLv2 license](http://www.gnu.org/licenses/old-licenses/lgpl-2.0-standalone.html).

# Installation
 To install **Millepede** (on a linux system):
 1. Download the software package from the DESY \c gitlab server to
    \a target directory, e.g. (shallow clone):

         git clone --depth 1 --branch V04-17-07 \
             https://gitlab.desy.de/millepede/millepede-ii.git target

 2. Create **Pede** executable (in \a target directory):

         make pede

 3. Optionally check the installation by running the simple test case:

         ./pede -t

    This will create (and use) the necessary text and binary files.

 Alternatively tarballs can be found [here](http://www.desy.de/~kleinwrt/MP2/tar).

# News

\subpage changelog_page

# Tools
 The subdirectory \c tools contains some useful scripts:
 * \c readMilleBinary.py: Python script to read binary files and print
   records in text form.
 * \c compareResults.py: Python3 script to compare result files (<tt>millepede.res</tt>).
 * \c readPedeHists.C: ROOT script to read and convert the **Millepede**
   histogram file <tt>millepede.his</tt>.
 * \c lapack: Test programs to print LAPACK (library) configuration (MKL, OpenBLAS).
 * \c tinypede.py: **Pede** implementation in python3 with basic functionality for
   illustration or testing with small problems.

# Julia
 The subdirectory \c julia contains some of the above tools reimplemented in [Julia](https://julialang.org):
 * \c readMilleBinary.jl
 * \c tinypede.jl

 # Details

 Detailed information is available at:

 \subpage draftman_page

 \subpage changes_page - documentation of new features added after the draft manual.

 \subpage option_page 

 \subpage exit_code_page

 \subpage troubleshooting_page

 \subpage test_brlf_page "Example"

 # Contact

 For information exchange the **Millepede** mailing list
 anacentre-millepede2@desy.de should be used.

# Legacy
 The subdirectory \c legacy contains the original \ref millepede1.f90
 "Millepede-I" implementation from Volker
 Blobel (2000) and a Millepede-I to Millepede-II \ref mp1to2.f90 "interface"
 (creating Millepede-II input files from Millepede-I calls, developed for COMPASS at CERN).

# References

 1. A New Method for the High-Precision Alignment of Track Detectors,
    Volker Blobel and Claus Kleinwort, Proceedings of the Conference on
    Adcanced Statistical Techniques in Particle Physics, Durham, 18 - 22 March 2002,
    Report DESY 02-077 (June 2002) and
    [hep-ex/0208021](http://arxiv.org/abs/hep-ex/0208021)
 2.  Alignment Algorithms, V. Blobel,
    [Proceedings](http://cdsweb.cern.ch/search?p=reportnumber%3ACERN-2007-004)
    of the LHC Detector Alignment Workshop, September 4 - 6 2006, CERN
 3. Software alignment for Tracking Detectors, V. Blobel,
    NIM A, 566 (2006), pp. 5-13,
    [doi:10.1016/j.nima.2006.05.157](http://dx.doi.org/10.1016/j.nima.2006.05.157)
 4. A new fast track-fit algorithm based on broken lines, V. Blobel,
    NIM A, 566 (2006), pp. 14-17,
    [doi:10.1016/j.nima.2006.05.156](http://dx.doi.org/10.1016/j.nima.2006.05.156)
 5. Millepede 2009, V. Blobel,
    [Contribution](https://indico.cern.ch/conferenceOtherViews.py?view=standard&confId=50502)
    to the 3rd LHC Detector Alignment Workshop, June 15 - 16 2009, CERN
 6. General Broken Lines as advanced track fitting method, C. Kleinwort,
    NIM A, 673 (2012), pp. 107-110,
    [doi:10.1016/j.nima.2012.01.024](http://dx.doi.org/10.1016/j.nima.2012.01.024)
 7. Volker Blobel und Erich Lohrmann, Statistische und numerische Methoden der
    Datenanalyse, Teubner Studienb&uuml;cher, B.G. Teubner, Stuttgart, 1998.
    [Online-Ausgabe](http://www.desy.de/~blobel/eBuch.pdf).
 8. [Systems Optimization Laboratory](http://web.stanford.edu/group/SOL/software/minres),
    Stanford University;\n
    C. C. Paige and M. A. Saunders (1975),
    Solution of sparse indefinite systems of linear equations,
    SIAM J. Numer. Anal. 12(4), pp. 617-629.
 9. [Systems Optimization Laboratory](http://web.stanford.edu/group/SOL/software/minresqlp),
    Stanford University;\n
    Sou-Cheng Choi, Christopher Paige, and Michael Saunders,
    MINRES-QLP: A Krylov subspace method for indefinite or singular
    symmetric systems, SIAM Journal of Scientific Computing 33:4, 1810-1836, 2011,
    [doi:10.1137/100787921](http://dx.doi.org/10.1137/100787921)


 \subpage exit_code_page
 
 \subpage troubleshooting_page
