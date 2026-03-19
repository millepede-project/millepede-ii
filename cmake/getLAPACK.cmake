# set which LAPACK versions we allow
if (LAPACK_MKL)
set (__LP_ALLOW_MKL True)
set (__LP_ALLOW_BLAS False) 
set (__LP_ALLOW_OTHER False)
elseif (LAPACK_OPENBLAS)
set (__LP_ALLOW_MKL False)
set (__LP_ALLOW_BLAS True) 
set (__LP_ALLOW_OTHER False)
else() 
set (__LP_ALLOW_MKL True)
set (__LP_ALLOW_BLAS True) 
set (__LP_ALLOW_OTHER True)
endif()

# always require int64 interface. 
set(BLA_SIZEOF_INTEGER 8) 

# MKL prioritised over OPENBLAS
if (LAPACK_MKL) 
    set (LAPACK_OPENBLAS False)
endif()

# PARDISO requires MKL as an explicit dependency
if (PARDISO)
    find_package( MKL REQUIRED )
endif()

# Now check possible LAPACK sources by descending priority 

# first: MKL 
if (${__LP_ALLOW_MKL})
    set (BLA_VENDOR Intel10_64_dyn)
    find_package(LAPACK QUIET )
    if (LAPACK_FOUND) 
        set (LAPACK_FLAVOUR "MKL")
    endif()
endif() 

# then, OpenBLAS
if ( (NOT LAPACK_FOUND) AND ${__LP_ALLOW_BLAS}) 
    set (BLA_VENDOR OpenBLAS)
    find_package(LAPACK QUIET )
    if (LAPACK_FOUND) 
        set (LAPACK_FLAVOUR "OPENBLAS")
    endif() 
endif() 

# Finally, anything else we can find 
if (NOT LAPACK_FOUND AND ${__LP_ALLOW_OTHER})
    unset(BLA_VENDOR) 
    find_package(LAPACK QUIET )
    if (LAPACK_FOUND) 
        set (LAPACK_FLAVOUR "OTHER")
    endif() 
endif()

if (LAPACK_FOUND) 
message("Found LAPACK: ${LAPACK_LIBRARIES}")
endif()


