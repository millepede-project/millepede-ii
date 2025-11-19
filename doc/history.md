\page history_page History

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

