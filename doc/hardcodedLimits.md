 \page limits_page Hardcoded limits

## Hardcoded limits

### Since V04-00-00 (fortran90)
For the input to pede the resources (memory, handles for C binary files, ..) are allocated 
dynamically and may be limited by the system. The are no explicit limits on the number of 
global parameters or number of (C) binary files.

For the internal monitoring some arrays still have fixed size:
- Maximum number of points for xy-scatter data (e.g. used for chi^2/ndf per file plots): 
  200 (500 since rev92) (`JFLC(5,IG) = 200` in 'mphistab.F', take care that 
  `NARR = NUMGXY * "desiredValue""` - if more points requested, points are averaged).
  
Since version V04-09-03 the memory allocation for internal histograms is now dynamic 
in terms of number of binary files.

For the steering text files the maximal length of a line is 1024 characters. This limits 
the length of file names too.

### Before V04-00-00 (Fortran77)
Since old versions of pede are Fortran77, several array lengths have a fixed size, e.g.
- Maximum number of files (binary and text files): 500 (`PARAMETER (MFILES=500`) in 
  'mpinds.inc'), maximum number of C-binaries is 490 (`#define MAXNUMFILES 490` in 
  'readc.c')
- Maximum number of characters of file names: 500 (`CHARACTER*500 TFD(MFILES)` in 
  'mpinds.inc')
- Maximum number of points for xy-scatter data (e.g. used for chi^2/ndf per file plots):
  200 (500 since rev92) (`JFLC(5,IG) = 200` in 'mphistab.F', take care that 
  `NARR = NUMGXY * "desiredValue"` - if more points requested, points are averaged).

The given values refer to tag V03-04-04.
