#= (doxygen like documentation style)
## \file
# Tiny PEDE implementation
#
# \author Claus Kleinwort, DESY, 2025 (Claus.Kleinwort@desy.de),
#
#  \copyright
#  Copyright (c) 2025 Deutsches Elektronen-Synchroton,
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
=#

"""
tiny PEDE implementation
"""
module Pede

using Distributions

# Read MP2 record
include("readRecord.jl")

# MP-II initialisation
include("pedeInit.jl")

# MP-II steering
include("pedeSteer.jl")

"""
Define (variable) global parameters

Read all binary files to get list of all global parameters.
Apply (minimum number of) entries cut and list of fixed parameters
to get list of free/variable parameters.
"""
function defineParameters()
	# loop over records    
	for (fileName, fileType, maxRec) in binaryFiles
		if fileType == 'F'
			println(" Reading Fortran binary file ", fileName)
		else
			println(" Reading binary file ", fileName)
		end
		io = open(fileName, "r")
		nRec = 0
		while !eof(io) && (nRec < maxRec || maxRec < 0)
			nRec += 1
			meas, inder, glder = readRecord(io, fileType != 'F')
			# count global labels
			for indices in meas
				for i in indices[2]+1:indices[3]-1
					if !haskey(parCounters, inder[i])
						parCounters[inder[i]] = 0
					end
					parCounters[inder[i]] += 1
				end
			end
		end
		println("   records found ", nRec)
		close(io)
		Pede.numRec += nRec
	end

	# fix parameters
	for label in fixedParLabels
		if haskey(parCounters, label)
			parCounters[label] *= -1
		end
	end
	# get (indices for) variable parameters
	index = 0
	for (label, counts) in sort(collect(parCounters))
		if counts >= entriesCut
			index += 1
			parIndices[label] = index
		end
	end
	println(" Found in binary files")
	println("   total number of records    ", numRec)
	println("   total number of parameters ", length(parCounters))
	println("   minimum number of entries  ", entriesCut)
	println("   number of free parameters  ", length(parIndices))
	println()
end

"""
Add constraints 

Add (linear equality) constraints from text file.
Accept only constraints with at least one variable parameter.
Assume zero r.h.s. (C*p=0). 

"""
function addConstraints()
	# loop over records    
	for fileName in constraintFiles
		nlines = 0
		nelem = 0
		ncons = 0
		println(" Reading constraint file ", fileName)
		io = open(fileName, "r")
		for line in readlines(io)
			nlines += 1
			fields = split(line, " ", keepempty = false)
			# to few fields ?
			if length(fields) < 2
				continue
			end
			# new constraint ?
			if lowercase(fields[1]) == "constraint"
				ncons += 1
				if nelem > 0
					Pede.numCons += 1
				end
				nelem = 0
			else
				# pair(label, value)
				label = parse(Int32, fields[1])
				value = parse(Float64, fields[2])
				# variable parameter ?
				if haskey(parIndices, label)
					nelem = nelem + 1
					push!(consElements, (numCons + 1, parIndices[label], value))
				end
			end
		end
		if nelem > 0
			Pede.numCons += 1
		end
		println("   number of constr. in file  ", ncons)
		println("   number of accepted constr. ", numCons)
		close(io)
	end
end

"""
Construct linear equation system

Perform local fits, reject bad records (local fit failed) and
records with too large Chi2(ndf) and build up global matrix and vector.
"""
function constructGlobalEquationSystem()
	println()
	println(" Constructing linear equations system")
	println("   local chi2 scaling factor  ", chi2Factor)
	# reset
	fill!(globalMatrix, 0.0)
	fill!(globalVector, 0.0)
	"Number "
	Pede.numAccepted = 0
	Pede.numRejected = 0
	Pede.numBad = 0
	Pede.sumNdf = 0
	Pede.sumChi2 = 0.0
	# loop over records
	for (fileName, fileType, maxRec) in binaryFiles
		io = open(fileName, "r")
		nRec = 0
		while !eof(io) && (nRec < maxRec || maxRec < 0)
			nRec += 1
			meas, inder, glder = readRecord(io, fileType != 'F')
			ndf, Chi2 = localFit(meas, inder, glder, true)
			# record accepted?
			if ndf > 0
				Pede.numAccepted += 1
				Pede.sumNdf += ndf
				Pede.sumChi2 += Chi2
			end
		end
	end
	# record statistics
	println("   number of rejected records ", numRejected)
	println("   number of bad records      ", numBad)
	# matrix information
	matSize = numVarPar + numCons
	nonZero = count(x -> x != 0.0, globalMatrix)
	println("   size of global matrix      ", matSize)
	println("   fraction of elements <> 0. ", trunc(100.0 * nonZero / matSize^2), '%')
	# add Lagrange multipliers
	for (iCons, label, value) in consElements
		globalMatrix[numVarPar+iCons, label] = value
		globalMatrix[label, numVarPar+iCons] = value
	end
end

"""
Perform local fit

Build and solve local linear equation system for a (MP2) record.
Optionally reject records depending on quality of local fit.
Update objective function (sum(ndf), sum(Chi2)) from accepted local fits
and (optionally) the global vector and matrix. 

### Arguments
- `record`: list of measurements (indices), labels and derivatives
- `updateGobalMatrix::bool`: update global matrix

"""
function localFit(measurements, indices, derivatives, updateGlobalMatrix)
	# number of local parameters
	numLocal = 0
	# number of measurements
	numMeas = 0
	# 1st loop over measurements: counting
	for (indMeas, indErr, indEnd) in measurements
		numLocal = max(numLocal, indices[indErr-1])
		numMeas += 1
	end
	#println("$numMeas $numLocal")
	#=
	Normal equaltions: D'*W*D = D'*W*m
	 - D: local derivatives
	 - W: weight/precision matrix of measurements
	 - m: measurements

	Measurements are independent, W is a diagonal matrix. Splitting into 2 parts (W=S'*S)
	leads to simplified normal equation: (S*D)'*(S*D) = (S*D)'*(S*m) 
	=#
	# measurements and local derivatives (scaled by error)
	vecMeasS = zeros(numMeas)
	locDerS = zeros(numMeas, numLocal)
	# 2nd loop over measurements: local deivatives
	iMeas = 0
	for (indMeas, indErr, indEnd) in measurements
		aMeas = derivatives[indMeas]
		aErr = derivatives[indErr]
		iMeas += 1
		vecMeasS[iMeas] = aMeas / aErr
		for i in indMeas+1:indErr-1
			locDerS[iMeas, indices[i]] = derivatives[i] / aErr
		end
	end
	localMatrix = locDerS' * locDerS
	localVector = locDerS' * vecMeasS
	#println(localMatrix)
	# local fit
	try
		# covariance matrix
		localCov = inv(localMatrix)
		# solution
		localSol = localCov * localVector
		# (scaled) residuals
		locResidualsS = vecMeasS - locDerS * localSol
		# chi2, ndf
		locNdf = numMeas - numLocal
		locChi2 = locResidualsS' * locResidualsS
		#println("local fit, ndf $locNdf Chi2 $locChi2")
		prob = 1.0 - cdf(Chisq(locNdf), locChi2 / chi2Factor)
		# reject (at 3 sigma)
		if prob < 0.0027
			Pede.numRejected += 1
			return -1, 0.0
		end
		# update global matrix?
		if updateGlobalMatrix
			# 3rd loop over measurements: get active global parameters
			# number of active variable global parameters
			numActive = 0
			for (indMeas, indErr, indEnd) in measurements
				for i in indErr+1:indEnd-1
					index = get(parIndices, indices[i], 0)
					# replace label by index in place
					indices[i] = index
					# variable?
					if index > 0
						if activeVarParIndices[index] == 0
							numActive += 1
							activeVarParList[numActive] = index
							activeVarParIndices[index] = numActive
						end
					end
				end
			end
			# sort active parameters
			activeGloPar = sort(activeVarParList[1:numActive])
			for i in 1:numActive
				activeVarParIndices[activeGloPar[i]] = i
			end
			#println("numActive $numActive :", activeGloPar)
			# 4th loop over measurements
			# global derivatives - only active parmeters!
			gloDerS = zeros(numMeas, numActive)
			iMeas = 0
			for (indMeas, indErr, indEnd) in measurements
				aErr = derivatives[indErr]
				iMeas += 1
				for i in indErr+1:indEnd-1
					index = indices[i]
					if index > 0
						relIndex = activeVarParIndices[index]
						gloDerS[iMeas, relIndex] = derivatives[i] / aErr
					end
				end
			end
			# compressed update matrix, vector (global parameter part)
			compressedUpdateMatrix = gloDerS' * gloDerS
			compressedUpdateVector = gloDerS' * vecMeasS
			# account for correlations (via local parameters)
			matGloLoc = gloDerS' * locDerS
			compressedUpdateMatrix -= matGloLoc * (localCov * matGloLoc')
			compressedUpdateVector -= matGloLoc * localSol
			# apply update ("expand")
			for i in 1:numActive
				index = activeGloPar[i]
				globalVector[index] += compressedUpdateVector[i]
				for j in 1:numActive
					globalMatrix[index, activeGloPar[j]] += compressedUpdateMatrix[i, j]
				end
			end
			# cleanup  (active global parameters)
			for i in 1:numActive
				activeVarParIndices[activeGloPar[i]] = 0
			end
		end

		return locNdf, locChi2

		# local fit failed
	catch e
		println("Inversion of local matrix failed: $e")
		println(" numMeas $numMeas numLocal $numLocal")
		Pede.numBad += 1
		return -2, 0.0
	end
end

"""
Solve linear equation system

Solve linear equation system (global fit) and print results.
"""
function solveGlobalEquationSystem()
	println()
	println(" Solving linear equations system")
	try
		# global covarinace matrix
		globalCov = inv(globalMatrix)
		# global solution
		globalSol = globalCov * globalVector
		# Chi2 reduction
		deltaChi2 = globalVector' * globalSol

		println("   sum(ndf)                 ", sumNdf - numVarPar)
		println("   initial sum(chi2)        ", sumChi2)
		println("   chi2 reduction by fit    ", deltaChi2)
		println("   final sum(chi2)          ", sumChi2 - deltaChi2)
		println("   final sum(chi2)/sum(ndf) ", (sumChi2 - deltaChi2) / (sumNdf - numVarPar))
		println()

		# write results to text file
		io = open("tinypede.res", "w")
		println(io, " Solution")
		println(io, "   parameter label, #entries, correction, error")

		for (label, counts) in sort(collect(parCounters))
			if haskey(parIndices, label)
				index = parIndices[label]
				println(io, "$label $counts $(globalSol[index]) $(sqrt(globalCov[index,index]))")
			else
				println(io, "$label $counts fixed")
			end
		end
		close(io)

		# global fit failed	
	catch e
		println("Inversion of global matrix failed: $e")
		println("   >>> improve input data <<<")
		return
	end
end

"""
Pede step

Implements only very basic **PEDE** functionality. 
 - meant as demonstration (of MP2 math/method)
 - plain julia arrays
 - optimisations:
   + exploit sparsitity of global derivatives (in single record)

Performs local and global fits. Features: 
 - fix parameters
 - apply (linear equality) constraints (using Lagrange multipliers)

"""
function pede()
	#=
	initialise
	=#
	startTime = time()
	println("\n TinyPede - a simple PEDE implementation in julia\n")
	#=
		get variable global parameters
	=#
	defineParameters()
	Pede.numVarPar = length(parIndices)
	Pede.activeVarParList = zeros(Int32, Pede.numVarPar)
	Pede.activeVarParIndices = zeros(Int32, Pede.numVarPar)
	println("time elapsed $(time()-startTime)")

	#=
		add constraints
	=#
	addConstraints()
	#println("consElements $consElements")
	println("time elapsed $(time()-startTime)")

	#=
		construct global (linear) equation system, Lagrange multipliers
	=#
	Pede.globalMatrix = Matrix{Float64}(undef, numVarPar + numCons, numVarPar + numCons)
	Pede.globalVector = Vector{Float64}(undef, numVarPar + numCons)
	constructGlobalEquationSystem()
	println("time elapsed $(time()-startTime)")

	#=
		solve global (linear) equation system
	=#
	solveGlobalEquationSystem()
	println("time elapsed $(time()-startTime)")
end

pede()
end
