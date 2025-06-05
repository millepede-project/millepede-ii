"""
Read (and decode) MP2 binary record, 'special data (blocks)' NOT supported!

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
	#  ignore error counter
	inder[1] = -1
	#  find zero indices
	ind0 = findall(iszero, inder)
	#  add begin of next record
	push!(ind0, nw + 1)
	#  define measurements
	measurements = [(ind0[i], ind0[i+1], ind0[i+2]) for i in 1:2:length(ind0)-1]
	return measurements, inder, glder
end
