"""
Read (and decode) MP2 binary record

### Arguments
- `io`: file to read from
- `isCFile::bool`: is C type file (byte stream)

### Returns 
- list of measurements (indices), labels and derivatives
"""
function readRecord(io, isCfile)
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

	# decode record
	# - step 1: count measurements
	numMeas = 0
	# internal indices
	i = 1
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
		numMeas += 1
		i -= 1
	end
	measurements = Vector{Tuple{Int32, Int32, Int32}}(undef, numMeas)
	# - step 2: fill measurements (indices)
	numMeas = 0
	# internal indices
	i = 1
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
		numMeas += 1
		#println("numMeas $numMeas $ja $jb $i")
		measurements[numMeas] = (ja, jb, i)
		i -= 1
	end
	return measurements, inder, glder
end
