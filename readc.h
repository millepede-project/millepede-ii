/** \file
 *  Read from C/C++ binary files - headers.
 *
 * \author Gero Flucke, University Hamburg, 2006
 * \author Claus Kleinwort, DESY (maintenance and developement)
 *
 *  \copyright
 *  Copyright (c) 2009 - 2023 Deutsches Elektronen-Synchroton,
 *  Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
 *  This library is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version. \n\n
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details. \n\n
 *  You should have received a copy of the GNU Library General Public
 *  License along with this program (see the file COPYING.LIB for more
 *  details); if not, write to the Free Software Foundation, Inc.,
 *  675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  C-methods to handle input of C/C++ binary files as input for
 *  the fortran **pede** program (see \ref peread).
 *  \deprecated
 *  This includes macros utilising \c cfortran.h to allow direct callability
 *  from fortran.
 *
 *  \c initc() has to be called once in the beginning,
 *  followed by one or several calls to \c openc() to open one or several files.
 *  \c readc() is then called to read the records sequentially. \c resetc()
 *  allows to rewind files.
 *
 *  If compiled with preprocessor macro \c USE_SHIFT_RFIO, uses \c libRFIO,
 *  i.e. includes \c shift.h instead of \c stdio.h
 *
 *  If compiled with preprocessor macro \c USE_ZLIB, uses \c libz,
 *  enables direct reading of gzipped files.
 *
 *  Written by Gero Flucke (gero.flucke@cern.ch) in 2006/7
 *  - update on July 14th, 2008
 *  - update on October 29th, 2008: return for file number in \c readC()
 *
 *  Major updates on April 24th, 2012 by C.Kleinwort:
 *  - skip records larger than buffer size (to determine max record length)
 *  - dynamic allocation of file pointer list (no hard-coded max number of files)
 *
 *  Major update on February 26th, 2014 by C.Kleinwort:
 *  - implement reading of records containing doubles (instead of floats)
 *    indicated by negative record length.
 *
 *  Major update on April 10th, 2019 by C.Kleinwort:
 *  - Option to close and reopen files 
 *
 *  Last major update on March 21th, 2023 by C.Kleinwort:
 *  - Fortran/C interoperability uses now 'iso_c_binding' (fortran 2003) instead of 'cfortran.h'
 *    (Proper string termination (file names) by thomas.white@desy.de)
 */
/*______________________________________________________________*/
/// Initialises the 'global' variables used for file handling.
#ifdef __cplusplus
extern "C"
{
#endif
/**
 * \param[in]  nFiles  Maximal number of C binary files to use.
 */
 void initc(int nFiles);

/*______________________________________________________________*/
/// Open file.
/**
 * \param[in]  fileName  File name
 * \param[in]  lengthFileName  Length of file name
 * \param[in]  nFileIn  File number (1 .. maxNumFiles) or <=0 for next one
 * \param[out] errorFlag error flag:
 *      * 0: if file opened and OK,
 *      * 1: if too many files open,
 *      * 2: if file could not be opened
 *      * 3: if file opened, but with error (can that happen?)
 */
 void openc(const char *fileName, int lengthFileName, int nFileIn, int *errorFlag); 
/*______________________________________________________________*/
/// Close file.
/**
 * \param[in]  nFileIn  File number (1 .. maxNumFiles)
 */
void closec(int nFileIn);
/*______________________________________________________________*/
/// Rewind file.
/**
 * \param[in]  nFileIn  File number (1 .. maxNumFiles)
 */
void resetc(int nFileIn);
/*______________________________________________________________*/
/// Read record from file.
/**
 * \param[out]    bufferDouble  read buffer for doubles
 * \param[out]    bufferFloat   read buffer for floats
 * \param[out]    bufferInt     read buffer for integers
 * \param[in,out] lengthBuffers in: buffer length, out: number of floats/ints in records
 *                              (> buffer size: record skipped)
 * \param[in]     nFileIn       File number (1 .. maxNumFiles)
 * \param[out]    errorFlag     error flag:
 *      *  -1: pointer to a buffer or lengthBuffers are null
 *      *  -2: problem reading record length
 *      *  -4: given buffers too short for record
 *      *  -8: problem with stream or EOF reading floats
 *      * -16: problem with stream or EOF reading ints
 *      * -32: problem with stream or EOF reading doubles
 *      *  =0: reached end of file (or read empty record?!)
 *      *  =4: found floats
 *      *  =8: found doubles
 */
/* No return value since to be called as subroutine from fortran,
    negative *errorFlag are errors, otherwise fine.

    *nFileIn: number of the file the record is read from,
    starting from 1 (not 0)
    */
void readc(double *bufferDouble, float *bufferFloat, int *bufferInt,
		int *lengthBuffers, int nFileIn, int *errorFlag);
#ifdef __cplusplus
}
#endif
