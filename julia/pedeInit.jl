#
# MP2 initialization
#
"Number of records"
numRec = 0
"Number of variable global parameters"
numVarPar = 0
"Number of constraints"
numCons = 0
"Number of accepted records"
numAccepted = 0
"Number of rejected recods"
numRejected = 0
"Number of bad records (fit failed)"
numBad = 0
"Sum(Ndf)"
sumNdf = 0
"Sum(Chi2)"
sumChi2 = 0.0
"Global parameters (label->counts)"
parCounters = Dict{Int32, Int}()
"Variable global parameters (label->index)"
parIndices = Dict{Int32, Int32}()
"Constraints elements (~ sparse constraints matrix)"
consElements = Tuple{Int32, Int32, Float64}[]
"Global matrix"
globalMatrix = Float64[;;]
"Global matrix (symmetric view)"
globalMatrixSym = Symmetric(Float64[;;])
"Global vector"
globalVector = Float64[;]
"List of active variable global parameters in record"
activeVarParList = Int32[;]
"Counters of active variable global parameters in record"
activeVarParIndices = Int32[;]
