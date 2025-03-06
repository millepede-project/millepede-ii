#
# MP2-II steering
#
"List of binary files (nam, type, max. records)"
binaryFiles = [("mp2tst.bin", 'F', 10000)] # fortran file (from "pede -t"), max. 10000 records
"Entries cut for variable parameters (min. number of appearances in binay files)"
entriesCut = 25
"List of fixed parameters (labels)"
fixedParLabels = [] #[56, 1056]
"List of text files with constraints"
constraintFiles = ["mp2con.txt"]
#constraintFiles = String[] # no constraints
" Chi2 scaling factor"
chi2Factor = 50.0
