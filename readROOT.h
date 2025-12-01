/** \file
 *  Read from ROOT files - headers.
 *
 * \author Maximilian Goblirsch-Kolb, DESY 
 *
 *  \copyright
 *  Copyright (c) 2025 Deutsches Elektronen-Synchroton,
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
 *  C-methods to handle input of ROOT files as input for
 *  the fortran **pede** program (see \ref peread).
 * 
 *  Follows the interface conventions of the \ref readc module:  
 *  \c initroot() is called once in the beggining, calls to \c openroot() 
 *  open input files, and subsequent calls to \c readroot() read one record
 *  at a time. \c resetroot() rewinds the files to their first entry, and 
 *  \c closeroot closes the files. 
 *
 */
/*______________________________________________________________*/
/// Avoid C++ name mangling for fortran interoperability 
#ifdef __cplusplus
extern "C"
{
#endif

/**
 * \param[in]  nFiles  Maximal number of ROOT files to use.
 */
 void initroot(int nFiles);

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
 void openroot(const char *fileName, int lengthFileName, int nFileIn, int *errorFlag); 
/*______________________________________________________________*/
/// Close file.
/**
 * \param[in]  nFileIn  File number (1 .. maxNumFiles)
 */
void closeroot(int nFileIn);
/*______________________________________________________________*/
/// Rewind file.
/**
 * \param[in]  nFileIn  File number (1 .. maxNumFiles)
 */
void resetroot(int nFileIn);
/*______________________________________________________________*/
/// Read record from file.
/**
 * \param[out]    bufferDouble  read buffer for doubles
 * \param[out]    bufferFloat   read buffer for floats. Not used, but included for interface consistency with \ref readc. 
 * \param[out]    bufferInt     read buffer for integers
 * \param[in,out] lengthBuffers in: buffer length, out: number of floats/ints in records
 *                              (> buffer size: record skipped)
 * \param[in]     nFileIn       File number (1 .. maxNumFiles)
 * \param[out]    errorFlag     error flag:
 *      *  -1: pointer to a buffer or lengthBuffers are null
 *      *  -2: problem reading record length
 *      *  -3: Invalid input file 
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
void readroot(double *bufferDouble, float *bufferFloat, int *bufferInt,
		int *lengthBuffers, int nFileIn, int *errorFlag);
#ifdef __cplusplus
}
#endif
