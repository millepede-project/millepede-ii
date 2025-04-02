#= (doxygen like documentation style)
## \file
# Read millepede binary file and print records
#
# \author Claus Kleinwort, DESY, 2025 (Claus.Kleinwort@desy.de), julia version
# \author Claus Kleinwort, DESY, 2009-2022 (Claus.Kleinwort@desy.de), python version
# \author Thomas Eichlersmith, Univ. Minnesota, 2023 (eichl008@umn.edu), major upgrade for python3 and CLI
#
#  \copyright
#  Copyright (c) 2009 - 2025 Deutsches Elektronen-Synchroton,
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
#
# One can check how this performs against a corrupted mille data file
# by using `dd` to intentionally only copy the first N bytes of a file.
#
#   dd count=10 if=uncorrupted.bin of=corrupted.bin
#
# will work if 'uncorrupted.bin' is greater than 512*10 bytes large.
=#

using ArgParse
# to read gzipped binary files:
using GZip
# replace 'open()' by 'GZip.open()' 

# define command line arguments
s = ArgParseSettings()
@add_arg_table s begin
	"--num-records", "-n"
	help = "number of records (i.e. tracks) to print to terminal. Continue until the end of the file if negative."
	arg_type = Int
	default = 10
	"--skip-records", "-s"
	help = "number of records (i.e. tracks) to skip before starting to print."
	arg_type = Int
	default = 0
	"--min-val"
	help = "minimum value to print derivatives."
	arg_type = Float32
	"--quiet"
	help = "do not print any information to terminal"
	action = :store_true
	"--type"
	help = "type of binary file to read (\"c\" or \"fortran\" )"
	arg_type = String
	range_tester = x -> x in ["autodetect", "c", "fortran"]
	default = "autodetect"
	"filename"
	help = "name of binary file to read"
	arg_type = String
	required = true
end
parsed_args = parse_args(ARGS, s)

# get file type
isCfile = true
if parsed_args["type"] == "c"
	isCFile = true
elseif parsed_args["type"] == "fortran"
	isCFile = false
else
	# need to auto-detect
	io = GZip.open(parsed_args["filename"], "r")
	n1 = read(io, Int32) # fortran: total record length (bytes)
	n2 = read(io, Int32) # fortran: MP2 record length (words)
	close(io)
	isCfile = true  # C
	if n1 == 4 * (n2 + 1)
		isCfile = false  # Fortran
		println("Detected Fortran binary file")
	end
end

# read file
io = GZip.open(parsed_args["filename"], "r")
maxRec = parsed_args["num-records"]
skipRec = parsed_args["skip-records"]
nRec = 0

while !eof(io) && (nRec < maxRec + skipRec || maxRec < 0)
	global nRec += 1
	# fortran: skip total record length
	if !isCfile
		reclen = read(io, Int32)
	end

	nr = read(io, Int32)
	nw = div(abs(nr), 2)
	# read read nw floats and then nw integers
	if nr > 0
		glder = Vector{Float32}(undef, nw)
	else
		glder = Vector{Float64}(undef, nw)
	end
	read!(io, glder)
	inder = Vector{Int32}(undef, nw)
	read!(io, inder)

	# fortran: skip total record length
	if !isCfile
		reclen = read(io, Int32)
	end

	if nRec <= skipRec
		continue
	end

	if parsed_args["quiet"]
		continue
	end

	println(" === NR $nRec $nw")

	# no details, only header
	if maxRec < 0
		continue
	end

	# decode record
	i = 1
	nh = 0
	ja = 1
	jb = 1
	jsp = 1
	nsp = 0
	while i < nw - 1
		i += 1
		while (i <= nw) && (inder[i] != 0)
			i += 1
		end
		ja = i
		i += 1
		while (i <= nw) && (inder[i] != 0)
			i += 1
		end
		jb = i
		i += 1
		# special data ?
		if (ja + 1 == jb) && (glder[jb] < 0.0)
			jsp = jb
			nsp = int(-glder[jb])
			i += nsp - 1
			println(" ### spec. $nsp $(inder[jsp + 1:i + 1]) $(glder[jsp + 1:i + 1])")
			continue
		end
		while (i <= nw) && (inder[i] != 0)
			i += 1
		end
		i -= 1
		nh += 1
		if (jb < i)
			# measurement with global derivatives
			println(" -g- meas. $nh $(inder[jb + 1]) $(jb - ja - 1) $(i - jb) $(glder[ja]) $(glder[jb])")
		else
			# measurement without global derivatives
			println(" -l- meas. $nh $(inder[ja + 1]) $(jb - ja - 1) $(i - jb) $(glder[ja]) $(glder[jb])")
		end
		if ja + 1 < jb
			lab = Int32[]
			val = typeof(glder[1])[]
			for k ∈ ja+1:jb-1
				if parsed_args["min-val"] === nothing
					push!(lab, inder[k])
					push!(val, glder[k])
				elseif abs(glder[k]) > parsed_args["min-val"]
					push!(lab, inder[k])
					push!(val, glder[k])
				end
			end
			println(" local  ", lab)
			println(" local  ", val)
		end
		if jb + 1 < i + 1
			lab = Int32[]
			val = typeof(glder[1])[]
			for k ∈ jb+1:i
				if parsed_args["min-val"] === nothing
					push!(lab, inder[k])
					push!(val, glder[k])
				elseif abs(glder[k]) > parsed_args["min-val"]
					push!(lab, inder[k])
					push!(val, glder[k])
				end
			end
			println(" global  ", lab)
			println(" global  ", val)
		end
	end
end

# summary
if maxRec < 0 || nRec < maxRec
	println(" end of file after $nRec records")
end
