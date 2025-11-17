 \page exit_code_page List of exit codes
 The exit code and message of the **Pede** executable can be found in the
 file <tt>millepede.end</tt> :
    + <b>-1</b>   Still running or crashed
    + **00**   Ended normally
    + **01**   Ended with warnings (bad measurements)
    + **02**   Ended with severe warnings (insufficient measurements)
    + **03**   Ended with severe warnings (bad global matrix)
    + **04**   Ended with severe warnings (bad binary file(s))
    + **05**   Ended without solution (empty constraints)
    + **10**   Aborted, no steering file
    + **11**   Aborted, open error for steering file
    + **12**   Aborted, second text file in command line
    + **13**   Aborted, unknown keywords in steering file
    + **14**   Aborted, no binary files
    + **15**   Aborted, open error(s) for binary files
    + **16**   Aborted, open error(s) for text files
    + **17**   Aborted, file name too long
    + **18**   Aborted, read error(s) for binary files
    + **19**   Aborted, binary file(s) modified
    + **20**   Aborted, bad binary records
    + **21**   Aborted, no labels/parameters defined
    + **22**   Aborted, no variable global parameters
    + **23**   Aborted, bad matrix index
    + **24**   Aborted, vector/matrix size mismatch
    + **25**   Aborted, result vector contains NaNs
    + **26**   Aborted, too many rejects
    + **27**   Aborted, singular QL decomposition of constraints matrix
    + **28**   Aborted, no local parameters
    + **29**   Aborted, factorization of global matrix failed
    + **30**   Aborted, memory allocation failed
    + **31**   Aborted, memory deallocation failed
    + **32**   Aborted, iteration limit reached in diagonalization
    + **33**   Aborted, stack overflow in quicksort
    + **34**   Aborted, pattern string too long - obsolete
    + **35**   Aborted, mismatch of number of global parameters
    + **36**   Aborted, unselected solution method requested
    + **40**   Aborted, other errors
