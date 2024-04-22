#!/usr/bin/env python

## \file
# Read millepede binary file and print records
#
# \author Claus Kleinwort, DESY, 2009-2022 (Claus.Kleinwort@desy.de)
# \author Thomas Eichlersmith, Univ. Minnesota, 2023 (eichl008@umn.edu), major upgrade for python3 and CLI
#
#  \copyright
#  Copyright (c) 2009 - 2023 Deutsches Elektronen-Synchroton,
#  Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
#  This library is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Library General Public License as
#  published by the Free Software Foundation; either version 2 of the
#  License, or (at your option) any later version. \n\n
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Library General Public License for more details. \n\n
#  You should have received a copy of the GNU Library General Public
#  License along with this program (see the file COPYING.LIB for more
#  details); if not, write to the Free Software Foundation, Inc.,
#  675 Mass Ave, Cambridge, MA 02139, USA.
#
# Use the `-h` or `--help` command options to print out the full list of options.
# Here are a few helpful examples for what this script is helpful for.
#
# Select a Few Tracks:
#
# You can look at specific tracks by changing the number of records to print
# and the number of records to skip before starting to print. Perhaps, you know
# that the 100'th track is "bad" for some reason, then you could look at it and
# its neighbors with
#
#   ./readMilleBinary.py --num_records 3 --skip-records 99 input.bin
#
# Check Format:
#
# You can check the format of a mille data file by running over all the records
# and keeping the printout quiet.
#
#   ./readMilleBinary.py --num-records -1 --quiet input.bin
#
# If the format is correct, it will just printout the number of records and
# have an exit status of 0. If it is incorrect, it will print an error and
# have an non-zero exist status. *Warning*: The binary format cannot _by itself_
# distinguish between a file that has been broken at a record boundary
# and one that has successfully been closed correctly. For this reason, more 
# thorough testing will be needed by the user if they wish to eliminate this 
# possibility.
#
# Description of the output from readMilleBinary.py
#    -  Records (tracks) start with \c '===' followed by record number and length 
#       (<0 for binary files containing doubles)
#    -  Measurements: A measurement with global derivatives is called a 'global measurement', 
#       otherwise 'local measurement'. Usually the real measurements from the detectors are 'global'
#       ones and virtual measurements e.g. to describe multiple scattering are 'local'.
#    -  'Global' measurements start with \c '-g-' followed by measurement number, first global label,
#       number of local and global derivatives, measurement value and error. The next lines contain 
#       local and global labels (array('i')) and derivatives (array('f') or array('d')).
#    -  'Local' measurements start with \c '-l-' followed by measurement number, first local label, 
#       number of local and global derivatives, measurement value and error. The next lines contain
#       local labels (array('i')) and derivatives (array('f') or array('d')).
#
# ---
#
# Tested with SL4, SL5, SL6
#
# Used containers to manually test Python 2.7.18 and Python 3.11.2
#
#   - Install [distrobox](https://distrobox.privatedns.org/)
#     - source code: https://github.com/89luca89/distrobox
#   - Construct a distrobox for the desired python version
#       `distrobox create -i python:2 -n py2`
#   - Enter it to do test runs
#       `distrobox enter py2`
#
# One can check how this performs against a corrupted mille data file
# by using `dd` to intentionally only copy the first N bytes of a file.
#
#   dd count=10 if=uncorrupted.bin of=corrupted.bin
#
# will work if 'uncorrupted.bin' is greater than 512*10 bytes large.

# in Python2.7, we can use the beta version of the print function
#    imports from __future__ need to happen before any other code
from __future__ import print_function

import sys

# CLI module distributed with Python
import argparse
# packing/unpacking structured binary data with Python
import struct

parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='read a mille binary file and print its data'
        )
parser.add_argument('--type', choices=['c', 'fortran', 'autodetect'], default='autodetect',
                    help='type of binary file that will be read')
parser.add_argument('filename', help='binary file to read')
parser.add_argument('-n', '--num-records', type=int, default=10,
                    help='Number of records (i.e. tracks) to print to terminal. Continue until the end of the file if negative.')
parser.add_argument('-s', '--skip-records', type=int, default=0,
                    help='number of records (tracks) to skip before starting to print.')
parser.add_argument('--min-val', type=float,
                    help='minimum value to print derivatives')
parser.add_argument('--quiet', action='store_true',
                    help='do not print any information to terminal')

arg = parser.parse_args()

# the argument parser makes sure the user provides one of the elements 
#   in the choices list, so if we get here, we can assume that arg.type
#   is one of them.
Cfiles = -1
if arg.type == 'c':
    CFiles = 1
elif arg.type == 'fortran':
    CFiles = 0
else:
    # need to auto-detect
    f = open(arg.filename, "rb")
    header_words = struct.unpack('ii', f.read(8))
    f.close()
    Cfiles = 1  # C
    if header_words[0] == 4 * (header_words[1] + 1): 
        Cfiles = 0  # Fortran
        print("Detected Fortran binary file")

# read file
f = open(arg.filename, "rb")


def unpack(typechar, number=1):
    """unpack a certain number of the input type

    We read from the file stored in the variable `f`.

    Raises
    ------
    EOFError
        if there is no data returned when the read is done
    ValueError
        if data is returned but the length does not match
        the amount of data requested

    Arguments
    ---------
    typechar : str
        single-character name of type to unpack: 'i', 'f', or 'd'
    number : int, optional
        number of that type to unpack, default one 

    Returns
    -------
    array of input length and type, filled with data unpacked
    from file
    """
    # 'i' and 'f' are 4 bytes, 'd' is 8 bytes
    bytesper = 4
    if typechar == 'd':
        bytesper = 8

    total_bytes = bytesper * number
    bin_data = f.read(total_bytes)
    if len(bin_data) != total_bytes:
        if len(bin_data) == 0:
            raise EOFError()
        else:
            raise ValueError('Requested %d bytes but only got %d from the file.' % (total_bytes, len(bin_data)))

    # struct knows that 'i' is 4 bytes, 'f' is 4 bytes, and 'd' is 8 bytes
    # https://docs.python.org/2.7/library/struct.html#format-characters
    return struct.unpack(typechar * number, bin_data)


nrec = 0
try:
    while (nrec < arg.num_records + arg.skip_records) or (arg.num_records < 0):
# read 1 record
        nr = 0
        if (Cfiles == 0):
            lenf = struct.unpack('i', f.read(4))

        try:
            length = unpack('i')
        except EOFError:
            # EOF is allowed on first word of format,
            #   that would mean we are done reading
            break

        # using bit-shifting instead of division since
        #   integer-division was promoted to its own operator
        #   in Python3, luckily shifting by 1 is the same as
        #   integer division by 2
        nr = abs(length[0] >> 1)
        nrec += 1

        floattype = 'f'
        if length[0] < 0:
            floattype = 'd'

        # read read nr floats and then nr integers
        glder = unpack(floattype, nr)
        inder = unpack('i', nr)

        if (Cfiles == 0):
            lenf = unpack('i')

        if (nrec <= arg.skip_records):  # must be after last fromfile
            continue

        if arg.quiet:
            continue

        print(" === NR ", nrec, length[0] / 2)

        # no details, only header
        if arg.num_records < 0:
            continue

        i = 0
        nh = 0
        ja = 0
        jb = 0
        jsp = 0
        nsp = 0
        while (i < (nr - 1)):
            i += 1
            while (i < nr) and (inder[i] != 0): i += 1
            ja = i
            i += 1
            while (i < nr) and (inder[i] != 0): i += 1
            jb = i
            i += 1
            # special data ?
            if (ja + 1 == jb) and (glder[jb] < 0.):
                jsp = jb
                nsp = int(-glder[jb])
                i += nsp - 1
                print(' ### spec. ', nsp, inder[jsp + 1:i + 1], glder[jsp + 1:i + 1])
                continue
            while (i < nr) and (inder[i] != 0): i += 1
            i -= 1
            nh += 1
            if (jb < i):
# measurement with global derivatives
                print(' -g- meas. ', nh, inder[jb + 1], jb - ja - 1, i - jb, glder[ja], glder[jb])
            else:
# measurement without global derivatives
                print(' -l- meas. ', nh, inder[ja + 1], jb - ja - 1, i - jb, glder[ja], glder[jb])
            if (ja + 1 < jb):
                lab = []
                val = []
                for k in range(ja + 1, jb):
                    if arg.min_val is None:
                        lab.append(inder[k])
                        val.append(glder[k])
                    elif abs(glder[k]) >= arg.min_val:
                        lab.append(inder[k])
                        val.append(glder[k])
                print(" local  ", lab)
                print(" local  ", val)
            if (jb + 1 < i + 1):
                lab = []
                val = []
                for k in range(jb + 1, i + 1):
                    if arg.min_val is None:
                        lab.append(inder[k])
                        val.append(glder[k]) 
                    elif abs(glder[k]) >= arg.min_val:
                        lab.append(inder[k])
                        val.append(glder[k])
                print(" global ", lab)
                print(" global ", val)

except EOFError:
    if (nr > 0):
        print(" >>> error: end of file before end of record", nrec)
        sys.exit(1)
except ValueError as e:
    print(" >>> error: unable to unpack values before end of record", nrec, *(e.args))
    sys.exit(2)

print(" end of file after", nrec, "records")
f.close()
