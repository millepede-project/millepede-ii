# Millepede-II

![MP2 logo](./mp2-logo.png)

Millepede II is a package for linear least squares fits with a large number of parameters. Developed for the alignment and calibration of tracking detectors.

Please check out the [Package Documentation](http://millepede.pages.desy.de/millepede-ii) for detailed instructions on how to set up and use Millepede-II.

## Availability 

The source code of Millepede II is provided by DESY under the terms of the [LGPLv2 license](http://www.gnu.org/licenses/old-licenses/lgpl-2.0-standalone.html).

The latest release recommended for use is `V04-17-07`. 

The recommended version is publicly available via the git command line:

`git clone --depth 1 --branch V04-17-07  https://gitlab.desy.de/millepede/millepede-ii.git MillepedeII`

For development, please clone the repository:

`git clone https://gitlab.desy.de/millepede/millepede-ii.git MillepedeII`

Please see the [documentation](https://millepede.pages.desy.de/millepede-ii/installation_page.html) for compilation instructions.  

In addition, docker images are available under `gitlab.desy.de:5555/millepede/millepede-ii/mp2runtime` 


- The software can be freely used for research and education. We expect that all publications describing work using this software quote at least one reference (see [here](http://www.desy.de/~blobel/mptalks.html) or [references](references)).

- Disclaimer: This software is provided without any expressed or implied warranty. In particular there is no warranty of any kind concerning the fitness of this software for any particular purpose.

## Requirements

### Fortran standards

Since version V04-13-00 at least fortran 2003 is required. The code complies with fortran 2023 (`gcc -std=f2023 -fall-intrinsics`).

### Optional dependencies

Optionally, the package can be compiled to use 
- [OpenMP](https://www.openmp.org/) for parallelisation, 
- [OpenBLAS](http://www.openmathlib.org/OpenBLAS/) or [Intel oneAPI MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) to use LAPACK routines for optimised solving, 
- [Intel MKL PARDISO](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-0/pardiso.html) to further improve the solution of sparse systems, 
- [ROOT](https://root.cern/) to visualise monitoring histograms

Instrumentation for profiling using [SCORE-P](https://www.vi-hps.org/projects/score-p/overview/overview.html) is implemented and can also be enabled during the build if desired. 

### Documentation

The documentation needs `doxygen` version 1.9.8 or later to correctly generate internal links between the pages.  

## History 

Millepede II is the successor of Millepede, a package for linear least squares fits with a large number of parameters. 

Both have been developed as experiment independent programs by Prof. V. Blobel (Univ. Hamburg) primary for the alignment and calibration of tracking detectors. Documentation and source code of Millepede and of the first version of Millepede II are available on [his web page](http://www.desy.de/~blobel/mptalks.html) (or [mirror](http://www.desy.de/~sschmitt/blobel/mptalks.html)). 

Meanwhile, the Statistics Tools group of the Analysis Centre of the Helmholtz [Terascale](http://terascale.de) Alliance took over the maintenance of Millepede II. It's listed in the [HIFIS spotlights](https://www.hifis.net/spotlights/millepede2).
The first Analysis Centre release, version V02-00-01, contains V. Blobel's 2009 development and small changes/fixes coming from CMS experience by G. Flucke.

## Recent development 
 Further development, mainly by C. Kleinwort for usage in CMS, adds
- auto-detection of bordered band matrices in the local fit to speed up processing if e.g. fed with tracks fitted using the broken lines technique (cf. [NIM A566:14-17,2006](http://www-library.desy.de/cgi-bin/spiface/find/hep/www?j=NUIMA,A566,14) by V. Blobel and [GeneralBrokenLines](https://gitlab.desy.de/claus.kleinwort/general-broken-lines/-/wikis/home)
- the possibility to ignore or compress rare off-diagonal elements of the global matrix (to reduce space requirements),
- shared-memory parallelisation of most parts of the code (track fits, MINRES, inversion) using [OpenMP™](https://www.openmp.org/)
- the possibility to read gzipped C/C++ binaries
- further fixes/improvements not only from CMS experience.

Using these features allows CMS silicon tracker alignment to determine 200 000 parameters in one go fitting more than 20 million tracks, including 3.6 million cosmic ray tracks and about 375 thousand muon pairs from Z decays reparametrised as one fit object with a vertex constraint and adding a virtual Z mass measurement (cf. CMS-CR-2011-323). The computing requirements on an Intel™ Xeon™ L5520 with 2.27 GHz using eight threads are 44.5 h CPU (but only 10 h wall clock time) and less than 32 GB of RAM.

## References
- A New Method for the High-Precision Alignment of Track Detectors, Volker Blobel and 
  Claus Kleinwort, Proceedings of the Conference on Adcanced Statistical Techniques in 
  Particle Physics, Durham, 18 - 22 March 2002, Report DESY 02-077 (June 2002) and 
  [hep-ex/0208021](http://arxiv.org/abs/hep-ex/0208021)
- Alignment Algorithms, V. Blobel, 
  [Proceedings](http://cdsweb.cern.ch/search?p=reportnumber%3ACERN-2007-004) of the 
  LHC Detector Alignment Workshop, September 4 - 6 2006, CERN
- Software alignment for Tracking Detectors, V. Blobel, NIM A, 566 (2006), pp. 5-13, 
  [doi:10.1016/j.nima.2006.05.157](http://dx.doi.org/10.1016/j.nima.2006.05.157)
