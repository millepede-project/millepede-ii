!> \file
!! Bit field counters.
!!
!! \author Volker Blobel, University Hamburg, 2005-2009 (initial Fortran77 version)
!! \author Claus Kleinwort, DESY (maintenance and developement)
!!
!! \copyright
!! Copyright (c) 2009 - 2024 Deutsches Elektronen-Synchroton,
!! Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
!! This library is free software; you can redistribute it and/or modify
!! it under the terms of the GNU Library General Public License as
!! published by the Free Software Foundation; either version 2 of the
!! License, or (at your option) any later version. \n\n
!! This library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Library General Public License for more details. \n\n
!! You should have received a copy of the GNU Library General Public
!! License along with this program (see the file COPYING.LIB for more
!! details); if not, write to the Free Software Foundation, Inc.,
!! 675 Mass Ave, Cambridge, MA 02139, USA.
!!
!! Count pairs of global parameters for sparse storage of global matrix,
!! apply pair entries cut and build (compressed) sparsity structure (row offsets, column lists).
!!
!! In sparse storage mode for each row the list of column indices with non zero elements
!! (and those elements) are stored. With compression this list is represented by the
!! first column and their number for continous regions (encoded in single INTEGER(mpi) words).
!! Rare elements may be stored in single precision.
!!
!! An additional bit map is used to monitor the parameter pairs for measurements (or 'equations').
!!
!! To make use of the Intel oneMKL PARDISO solver for sparse systems the 'CSR3'
!! (or 'BSR3') format (using the 'upper' triangle) is needed.
!! For this a new set of routines have been added (starting with 'P').
!! In this case the bit array has two parts: A triangular one for pairs of parameter groups
!! and a rectangular one for Lagrange multipliers (parameter constraint pairs).

#ifdef SCOREP_USER_ENABLE
#include "scorep/SCOREP_User.inc"
#endif

!> Bit field data.
MODULE mpbits
    USE mpdef
    IMPLICIT NONE

    INTEGER(mpl) :: ndimb  !< dimension for bit (field) array
    INTEGER(mpl) :: ndimb2 !< dimension for bit map
    INTEGER(mpi) :: n      !< matrix size (counters, sparse, triangular part)
    INTEGER(mpi) :: nar    !< additional rows (counters, sparse, rectangular part)
    INTEGER(mpi) :: nac    !< additional columns (counters, sparse, rectangular part)
    INTEGER(mpi) :: n2     !< matrix size (map)
    INTEGER(mpi) :: ibfw   !< bit field width
    INTEGER(mpi) :: ireqpe !< min number of pair entries
    INTEGER(mpi) :: isngpe !< upper bound for pair entry single precision storage
    INTEGER(mpi) :: iextnd !< flag for extended storage (both 'halves' of sym. mat. for improved access patterns)
    INTEGER(mpi) :: nspc   !< number of precision for sparse global matrix (1=D, 2=D+f)
    INTEGER(mpi) :: mxcnt  !< max value for bit field counters
    INTEGER(mpi) :: nthrd  !< number of threads
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: bitFieldCounters !< fit field counters for global parameters pairs (tracks)
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: bitMap !< fit field map for global parameters pairs (measurements)
    INTEGER(mpi), PARAMETER :: bs = BIT_SIZE(1_mpi)  !< number of bits in INTEGER(mpi)

END MODULE mpbits

!> Fill bit fields (counters, triangular part).
!!
!! \param [in]    im     first index
!! \param [in]    jm     second index
!! \param [in]    inc    increment (usually 1)
!!
SUBROUTINE inbits(im,jm,inc)        ! include element (I,J)
    USE mpbits
    IMPLICIT NONE

    INTEGER(mpi), INTENT(IN) :: im
    INTEGER(mpi), INTENT(IN) :: jm
    INTEGER(mpi), INTENT(IN) :: inc

    INTEGER(mpl) :: l
    INTEGER(mpi) :: i
    INTEGER(mpi) :: j
    INTEGER(mpi) :: noffj
    INTEGER(mpi) :: nout
    INTEGER(mpi) :: m
    INTEGER(mpi) :: icount
    INTEGER(mpi) :: jcount
    INTEGER(mpl) :: noffi

    ! diagonal included now !
    !IF(im == jm) RETURN  ! diagonal
    j=MIN(im,jm)
    i=MAX(im,jm)
    IF(j <= 0) RETURN    ! out low
    IF(i > n) RETURN    ! out high
    IF(iextnd >= 0) THEN
        ! lower triangle, row 'i', col 'j'
        noffi=INT(i-1,mpl)*INT(i,mpl)*INT(ibfw,mpl)/2 ! row offset for J=1
        noffj=(j-1)*ibfw
        l=noffi/bs+i+noffj/bs ! row offset + column offset
        !     add I instead of 1 to keep bit maps of different rows in different words (openMP !)
    ELSE
        ! upper triangle, row 'j', col 'i'
        noffi=INT(j-1,mpl)*INT(2*n+2-j,mpl)/2 ! row offset for I=J
        noffj=i-j
        l=noffi/bs+j+noffj/bs ! row offset + column offset
        !     add J instead of 1 to keep bit maps of different rows in different words (openMP !)
    ENDIF
    m=MOD(noffj,bs)
    IF (ibfw <= 1) THEN
        bitFieldCounters(l)=ibset(bitFieldCounters(l),m)
    ELSE
        !        get counter from bit field
        icount=0
        nout=m+ibfw-bs ! number of bits outside word
        IF (nout <= 0) THEN
            ! inside single word
            CALL mvbits(bitFieldCounters(l),m,ibfw,icount,0)
        ELSE
            ! spread over two words
            CALL mvbits(bitFieldCounters(l),m,ibfw-nout,icount,0)
            CALL mvbits(bitFieldCounters(l+1),0,nout,icount,ibfw-nout)
        ENDIF
        !        increment
        jcount=icount
        icount=MIN(icount+inc,mxcnt)
        !        store counter into bit field
        IF (icount /= jcount) THEN
            IF (nout <= 0) THEN
                ! inside single word
                CALL mvbits(icount,0,ibfw,bitFieldCounters(l),m)
            ELSE
                ! spread over two words
                CALL mvbits(icount,0,ibfw-nout,bitFieldCounters(l),m)
                CALL mvbits(icount,ibfw-nout,nout,bitFieldCounters(l+1),0)
            ENDIF
        END IF
    END IF
    RETURN

END SUBROUTINE inbits

!> Fill bit fields (counters, rectangular part).
!!
!! \param [in]    i     row (parameter)
!! \param [in]    j     column (constraint)
!!
SUBROUTINE irbits(i,j)        ! include element (I,J)
    USE mpbits
    IMPLICIT NONE

    INTEGER(mpi), INTENT(IN) :: i
    INTEGER(mpi), INTENT(IN) :: j

    INTEGER(mpl) :: l
    INTEGER(mpi) :: noffj
    INTEGER(mpi) :: m
    INTEGER(mpl) :: noffi

    IF(i > nar.OR.j > nac) RETURN    ! out high
    ! upper triangle, rectangular part
    noffi=ndimb+INT(i-1,mpl)*INT(nac/bs+1,mpl)+1 ! row offset for J=1
    noffj=j-1
    l=noffi+noffj/bs ! row offset + column offset
    m=MOD(noffj,bs)
    bitFieldCounters(l)=ibset(bitFieldCounters(l),m)
    RETURN

END SUBROUTINE irbits

!> Calculate bit (field) array size, encoding.
!!
!! \param [in]    in        matrix size
!! \param [in]    jreqpe    min number of pair entries
!! \param [in]    jhispe    upper bound for pair entry histogrammimg
!! \param [in]    jsngpe    upper bound for pair entry single precision storage
!! \param [in]    jextnd    flag for extended storage
!! \param [out]   idimb     dimension for bit (field) array
!! \param [out]   ispc      number of precision for sparse global matrix
!!
SUBROUTINE clbits(in,jreqpe,jhispe,jsngpe,jextnd,idimb,ispc)
    USE mpbits
    USE mpdalc
    IMPLICIT NONE

    INTEGER(mpi), INTENT(IN) :: in
    INTEGER(mpi), INTENT(IN) :: jreqpe
    INTEGER(mpi), INTENT(IN) :: jhispe
    INTEGER(mpi), INTENT(IN) :: jsngpe
    INTEGER(mpi), INTENT(IN) :: jextnd
    INTEGER(mpl), INTENT(OUT) :: idimb
    INTEGER(mpi), INTENT(OUT) :: ispc

    INTEGER(mpl) :: noffd
    INTEGER(mpi) :: i
    INTEGER(mpi) :: icount
    INTEGER(mpi) :: mb
    !$    INTEGER(mpi) :: OMP_GET_MAX_THREADS
    ! save input parameter 
    n=in
    nar=0
    nac=0
    ireqpe=jreqpe
    isngpe=jsngpe
    iextnd=max(0,jextnd)
    ! number of precision types (D, F)
    ispc=1
    ! bit field size
    ibfw=1 ! number of bits needed to count up to ICOUNT
    mxcnt=1
    IF (jextnd >= 0) THEN
        ! optional larger bit fields for lower triangle
        icount=MAX(jsngpe+1,jhispe)
        icount=MAX(jreqpe,icount)
        DO i=1,30
            IF (icount > mxcnt) THEN
                ibfw=ibfw+1
                mxcnt=mxcnt*2+1
            END IF
        END DO
        IF (jsngpe>0) ispc=2
    END IF
    ! bit field array size
    noffd=INT(n,mpl)*INT(n+1,mpl)*INT(ibfw,mpl)/2
    ndimb=noffd/bs+n
    idimb=ndimb
    nspc=ispc
    mb=INT(4.0E-6*REAL(ndimb,mps),mpi)
    WRITE(*,*) ' '
    WRITE(*,*) 'CLBITS: symmetric matrix of dimension',n
    WRITE(*,*) 'CLBITS: off-diagonal elements',noffd
    IF (mb > 0) THEN
        WRITE(*,*) 'CLBITS: dimension of bit-array',ndimb , '(',mb,'MB)'
    ELSE
        WRITE(*,*) 'CLBITS: dimension of bit-array',ndimb , '(< 1 MB)'
    END IF
    CALL mpalloc(bitFieldCounters,ndimb,'INBITS: bit storage')
    bitFieldCounters=0
    nthrd=1
    !$ NTHRD=OMP_GET_MAX_THREADS()
    RETURN
END SUBROUTINE clbits

!> Calculate bit field array size (PARDISO).
!!
!! Count pairs of parameter groups in triangular part,
!! pairs of constraints and parameters in rectangular part.
!!
!! \param [in]    in        matrix size (triangular part, number of parameter groups)
!! \param [in]    inar      rows for rectangular part (number of parameters)
!! \param [in]    inac      cols for rectangular part (number of constraints)
!! \param [out]   idimb     dimension for bit (field) array
!!
SUBROUTINE plbits(in,inar,inac,idimb)
    USE mpbits
    USE mpdalc
    IMPLICIT NONE

    INTEGER(mpi), INTENT(IN) :: in
    INTEGER(mpi), INTENT(IN) :: inar
    INTEGER(mpi), INTENT(IN) :: inac
    INTEGER(mpl), INTENT(OUT) :: idimb

    INTEGER(mpl) :: noffd
    INTEGER(mpi) :: mb
    !$    INTEGER(mpi) :: OMP_GET_MAX_THREADS
    ! save input parameter
    n=in
    nar=inar
    nac=inac
    ireqpe=1
    isngpe=0
    iextnd=-1 ! upper triangle
    ibfw=1
    ! bit field array size
    noffd=INT(n,mpl)*INT(n+1,mpl)/2
    ndimb=noffd/bs+n
    idimb=ndimb+INT(nar,mpl)*INT(nac/bs+1,mpl)
    nspc=1
    mb=INT(4.0E-6*REAL(idimb,mps),mpi)
    WRITE(*,*) ' '
    WRITE(*,*) 'PLBITS: symmetric matrix of dimension',n, '(',nar, nac,')'
    WRITE(*,*) 'PLBITS: off-diagonal elements',noffd, '(',INT(nar,mpl)*INT(nac,mpl),')'
    IF (mb > 0) THEN
        WRITE(*,*) 'PLBITS: dimension of bit-array',idimb , '(',mb,'MB)'
    ELSE
        WRITE(*,*) 'PLBITS: dimension of bit-array',idimb , '(< 1 MB)'
    END IF
    CALL mpalloc(bitFieldCounters,idimb,'INBITS: bit storage')
    bitFieldCounters=0
    nthrd=1
    !$ NTHRD=OMP_GET_MAX_THREADS()
    RETURN
END SUBROUTINE plbits

!> Analyze bit fields.
!!
!! \param [in]     npgrp   parameter groups
!! \param [out]    ndims   (1): (reduced) size of bit array; (2): size of column lists;
!!                         (3/4): number of (double/single precision) off diagonal elements;
!! \param[out]     nsparr  row offsets
!! \param[in]      ihst    >0: histogram number
!!
SUBROUTINE ndbits(npgrp,ndims,nsparr,ihst)
    USE mpbits
    USE mpdalc
    IMPLICIT NONE

    INTEGER(mpi), DIMENSION(:), INTENT(IN) :: npgrp
    INTEGER(mpl), DIMENSION(4), INTENT(OUT) :: ndims
    INTEGER(mpl), DIMENSION(:,:), INTENT(OUT) :: nsparr
    INTEGER(mpi), INTENT(IN) :: ihst

    INTEGER(mpi) :: nwcp(0:1)
    INTEGER(mpi) :: irgn(2)
    INTEGER(mpi) :: inr(2)
    INTEGER(mpi) :: ichunk
    INTEGER(mpi) :: i
    INTEGER(mpi) :: j
    INTEGER(mpi) :: jcol
    INTEGER(mpi) :: last
    INTEGER(mpi) :: lrgn
    INTEGER(mpi) :: next
    INTEGER(mpi) :: icp
    INTEGER(mpi) :: mm
    INTEGER(mpi) :: jp
    INTEGER(mpi) :: n1
    INTEGER(mpi) :: nd
    INTEGER(mpi) :: ib
    INTEGER(mpi) :: ir
    INTEGER(mpi) :: icount
    INTEGER(mpi) :: iproc
    INTEGER(mpi) :: iword
    INTEGER(mpi) :: mb
    INTEGER(mpl) :: ll
    INTEGER(mpl) :: lb
    INTEGER(mpl) :: nin
    INTEGER(mpl) :: npar
    INTEGER(mpl) :: ntot
    INTEGER(mpl) :: nskyln
    INTEGER(mpl) :: noffi
    INTEGER(mpl) :: noffj
    REAL(mps) :: cpr
    REAL(mps) :: fracu
    REAL(mps) :: fracs
    REAL(mps) :: fracz
    LOGICAL :: btest
    !$    INTEGER(mpi) :: OMP_GET_THREAD_NUM
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: lastRowInCol

    ll=INT(n,mpl)*INT(nthrd,mpl)
    CALL mpalloc(lastRowInCol,ll,'NDBITS: last (non zero) row in col.')
    lastRowInCol=0

    nd=npgrp(n+1)-npgrp(1) ! number of diagonal elements

    ndims(1)=ndimb
    ndims(2)=0
    ndims(3)=0
    ndims(4)=0
    ntot=0
    ll=0
    lb=0
    ichunk=MIN((n+nthrd-1)/nthrd/32+1,256)
    ! reduce bit field counters to (precision type) bits, analyze precision type bit fields ('1st half' (j<=i))
    ! parallelize row loop
    ! private copy of NTOT for each thread, combined at end, init with 0.
    !$OMP  PARALLEL DO &
    !$OMP  PRIVATE(I,NOFFI,LL,MM,LB,MB,IWORD,IPROC,J,ICOUNT,IB,INR,IRGN,LAST,LRGN,NEXT,JP,IR,NPAR) &
    !$OMP  REDUCTION(+:NTOT) &
    !$OMP  SCHEDULE(DYNAMIC,ICHUNK)
    DO i=1,n          
        noffi=INT(i-1,mpl)*INT(i,mpl)*INT(ibfw,mpl)/2
        ll=noffi/bs+i
        mm=0
        lb=ll
        mb=0
        iword=0 ! reset temporary bit fields
        iproc=0
        !$ IPROC=OMP_GET_THREAD_NUM()         ! thread number
        inr(1)=0
        inr(2)=0
        irgn(1)=1 ! 'end marker' region
        irgn(2)=1
        last=0
        lrgn=0
        npar=0

        DO j=1,i         !  loop until diagonal element                
            ! get (pair) counter
            icount=0
            next=0
            DO ib=0,ibfw-1
                IF (btest(bitFieldCounters(ll),mm)) icount=ibset(icount,ib)
                mm=mm+1
                IF (mm >= bs) THEN
                    ll=ll+1
                    mm=mm-bs
                END IF
            END DO

            IF (icount > 0) THEN
                npar=npar+npgrp(j+1)-npgrp(j)
                IF (iproc == 0.AND.ihst > 0) CALL hmpent(ihst,REAL(icount,mps))
            END IF

            ! keep pair ?
            IF (icount >= ireqpe) THEN
                next=1 ! double
                IF (icount <= isngpe) next=2 ! single
                iword=ibset(iword,mb+next-1)
                inr(next)=inr(next)+npgrp(j+1)-npgrp(j) ! number of parameters
                IF (next /= last) THEN
                    irgn(next)=irgn(next)+1
                END IF
                jcol=j+iproc*n
                lastRowInCol(jcol)=max(lastRowInCol(jcol),i)
            END IF
            last=next
            ! save condensed bitfield
            mb=mb+nspc
            IF (mb >= bs) THEN
                bitFieldCounters(lb)=iword ! store
                iword=0
                lb=lb+1
                mb=mb-bs
            END IF
        END DO
        bitFieldCounters(lb)=iword ! store
        ntot=ntot+npar*(npgrp(i+1)-npgrp(i))

        ! save row statistics
        ir=i+1
        DO jp=1,nspc
            nsparr(1,ir)=irgn(jp)    ! number of regions per row and precision
            nsparr(2,ir)=inr(jp)     ! number of columns per row and precision (groups)
            ir=ir+n+1
        END DO
        
    END DO
    !$OMP END PARALLEL DO

    ! analyze precision type bit fields for extended storage, check for row compression

    ! parallelize row loop
    ! private copy of NDIMS for each thread, combined at end, init with 0.
    !$OMP  PARALLEL DO &
    !$OMP  PRIVATE(I,NOFFI,NOFFJ,LL,MM,INR,IRGN,LAST,LRGN,J,NEXT,ICP,NWCP,JP,IR,IB) &
    !$OMP  REDUCTION(+:NDIMS) &
    !$OMP  SCHEDULE(DYNAMIC,ICHUNK)
    DO i=1,n
        ! restore row statistics
        ir=i+1
        DO jp=1,nspc
            irgn(jp)=INT(nsparr(1,ir),mpi)    ! number of regions per row and precision
            inr(jp)=INT(nsparr(2,ir),mpi)     ! number of columns per row and precision (groups)
            ir=ir+n+1
        END DO

        ! analyze precision type bit fields for extended storage ('2nd half' (j>i) too) ?
        IF (iextnd > 0) THEN

            noffj=(i-1)*nspc
            mm=INT(MOD(noffj,INT(bs,mpl)),mpi)

            last=0
            lrgn=0

            ! remaining columns
            DO j=i+1, n
                ! index for pair (J,I)
                noffi=INT(j-1,mpl)*INT(j,mpl)*INT(ibfw,mpl)/2 ! for I=1
                ll=noffi/bs+j+noffj/bs ! row offset + column offset

                ! get precision type
                next=0
                DO ib=0,nspc-1
                    IF (btest(bitFieldCounters(ll),mm+ib)) next=ibset(next,ib)
                END DO

                ! keep pair ?
                IF (next > 0) THEN
                    inr(next)=inr(next)+npgrp(j+1)-npgrp(j) ! number of parameters
                    IF (next /= last) THEN
                        irgn(next)=irgn(next)+1
                    END IF
                END IF
                last=next
            END DO
        END IF

        ! row statistics, compression
        ir=i+1
        DO jp=1,nspc
            icp=0
            nwcp(0)=inr(jp)                       ! list of column indices (default)
            IF (inr(jp) > 0) THEN
                nwcp(1)=irgn(jp)*2                ! list of regions (group starts and offsets)
                ! compress row ?
                IF ((nwcp(1) < nwcp(0)).OR.iextnd > 0) THEN
                    icp=1
                END IF
                ! total space
                ndims(2)   =ndims(2)   +nwcp(icp)
                ndims(jp+2)=ndims(jp+2)+nwcp(0)*(npgrp(i+1)-npgrp(i))
            END IF
            ! per row and precision
            nsparr(1,ir)=nwcp(icp)
            nsparr(2,ir)=nwcp(0)*(npgrp(i+1)-npgrp(i))
            ir=ir+n+1
        END DO
    END DO
    !$OMP END PARALLEL DO

    ! sum up, fill row offsets
    lb=0
    n1=0
    ll=nd
    DO jp=1,nspc
        DO i=1,n
            n1=n1+1
            nsparr(1,n1)=lb
            nsparr(2,n1)=ll
            lb=lb+nsparr(1,n1+1)
            ll=ll+nsparr(2,n1+1)
        END DO
        n1=n1+1
        nsparr(1,n1)=lb
        nsparr(2,n1)=ll
        ll=0
    END DO

    ! look at skyline
    ! first combine threads
    ll=n
    DO j=2,nthrd
        DO i=1,n
            ll=ll+1
            lastRowInCol(i)=max(lastRowInCol(i),lastRowInCol(ll))
        END DO
    END DO
    ! sumup
    nskyln=0
    DO i=1,n
        npar=npgrp(lastRowInCol(i)+1)-npgrp(i)
        nskyln=nskyln+npar*(npgrp(i+1)-npgrp(i))
    END DO
    ! cleanup
    CALL mpdealloc(lastRowInCol)

    nin=ndims(3)+ndims(4)
    fracz=200.0*REAL(ntot,mps)/REAL(nd,mps)/REAL(nd-1,mps)
    fracu=200.0*REAL(nin,mps)/REAL(nd,mps)/REAL(nd-1,mps)
    fracs=200.0*REAL(nskyln,mps)/REAL(nd,mps)/REAL(nd+1,mps)
    WRITE(*,*) ' '
    WRITE(*,*) 'NDBITS: number of diagonal elements',nd
    WRITE(*,*) 'NDBITS: number of used off-diagonal elements',nin
    WRITE(*,1000) 'fraction of non-zero off-diagonal elements', fracz
    WRITE(*,1000) 'fraction of used off-diagonal elements', fracu
    cpr=100.0*REAL(mpi*ndims(2)+mpd*ndims(3)+mps*ndims(4),mps)/REAL((mpd+mpi)*nin,mps)
    WRITE(*,1000) 'compression ratio for off-diagonal elements', cpr
    WRITE(*,1000) 'fraction inside skyline ', fracs
1000 FORMAT(' NDBITS: ',a,f6.2,' %')
    RETURN
END SUBROUTINE ndbits

!> Analyze bit fields.
!!
!! Check block structure (for PARDISO BSR3 storage)
!!
!! \param[in]      npgrp   parameter groups
!! \param[in]      ibsize  block size
!! \param[out]     nnzero  number of non-zero elements
!! \param[out]     nblock  number of blocks used
!! \param[out]     nbkrow  number of (column) blocks in row (blocks)
!!
SUBROUTINE pbsbits(npgrp,ibsize,nnzero,nblock,nbkrow)
    USE mpbits
    USE mpdalc
    IMPLICIT NONE

    INTEGER(mpi), DIMENSION(:), INTENT(IN) :: npgrp
    INTEGER(mpi), INTENT(IN) :: ibsize
    INTEGER(mpl), INTENT(OUT) :: nnzero
    INTEGER(mpl), INTENT(OUT) :: nblock
    INTEGER(mpi), DIMENSION(:),INTENT(OUT) :: nbkrow

    INTEGER(mpi) :: ichunk
    INTEGER(mpi) :: i
    INTEGER(mpi) :: ib
    INTEGER(mpi) :: igrpf
    INTEGER(mpi) :: igrpl
    INTEGER(mpi) :: iproc
    INTEGER(mpi) :: ioffb
    INTEGER(mpi) :: ioffg
    INTEGER(mpi) :: j
    INTEGER(mpi) :: mb
    INTEGER(mpi) :: mbt
    INTEGER(mpi) :: mm
    INTEGER(mpi) :: nd
    INTEGER(mpi) :: ngrp
    INTEGER(mpi) :: ir
    INTEGER(mpi) :: irfrst
    INTEGER(mpi) :: irlast
    INTEGER(mpi) :: jb
    INTEGER(mpi) :: jc
    INTEGER(mpl) :: length
    INTEGER(mpl) :: ll
    INTEGER(mpl) :: noffi
    INTEGER(mpl), PARAMETER :: two=2
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: rowBlocksToGroups
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: blockCounter
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: groupList
    !$    INTEGER(mpi) :: OMP_GET_THREAD_NUM

    LOGICAL :: btest

    nnzero=0
    nblock=0
    nbkrow=0

    !$POMP INST BEGIN(pbsbits)
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_BY_NAME_BEGIN("pbsbits", SCOREP_USER_REGION_TYPE_COMMON)
#endif
    nd=npgrp(n+1)-npgrp(1) ! number of diagonal elements
    mb=(nd+nac-1)/ibsize+1 ! max. number of blocks per row/column
    mbt=(nd-1)/ibsize+1 ! max. number of blocks in triangular part
    length=INT(mb,mpl)*INT(nthrd,mpl)
    CALL mpalloc(blockCounter,length,'PBBITS: block counter') 
    length=INT(n,mpl)*INT(nthrd,mpl)
    CALL mpalloc(groupList,length,'PBBITS: group list')

    ! mapping row blocks to parameters groups
    length=INT(mbt,mpl)
    CALL mpalloc(rowBlocksToGroups,two,length,'mapping row blocks to par. groups (I)')
    rowBlocksToGroups(:,:)=0
    igrpf=1 ! first group of block
    igrpl=1 ! last  group of block
    ir=1    ! first row   of block
    DO i=1,mbt
        DO WHILE (igrpf < n .AND. npgrp(igrpf+1) <= ir)
            igrpf=igrpf+1
        END DO
        rowBlocksToGroups(1,i)=igrpf
        ir=ir+ibsize
        DO WHILE (igrpl < n .AND. npgrp(igrpl+1) < ir)
            igrpl=igrpl+1
        END DO
        rowBlocksToGroups(2,i)=igrpl
    END DO

    ll=0
    ichunk=MIN((mbt+nthrd-1)/nthrd/32+1,256)
    ! analyze bit fields ('upper half' (j>=i))
    ! parallelize row loop
    !$OMP  PARALLEL DO &
    !$OMP  PRIVATE(IPROC,IOFFB,IOFFG,I,NOFFI,LL,MM,NGRP,J,IR,IRFRST,IRLAST) &
    !$OMP  REDUCTION(+:NNZERO,NBLOCK) &
    !$OMP  SCHEDULE(DYNAMIC,ICHUNK)
    DO ib=1,mbt
        irfrst=ibsize*(ib-1)+1  ! first row in block
        irlast=ibsize*ib        ! last row in block
        iproc=0
        !$ IPROC=OMP_GET_THREAD_NUM()         ! thread number
        ioffb=iproc*mb
        ioffg=iproc*n
        blockCounter(ioffb+1:ioffb+mb)=0
        DO i=rowBlocksToGroups(1,ib),rowBlocksToGroups(2,ib)
            noffi=INT(i-1,mpl)*INT(2*n+2-i,mpl)/2 ! row offset for I=J
            ll=noffi/bs+i
            mm=0
            ngrp=0
            DO j=i,n         !  loop from diagonal element
                IF (btest(bitFieldCounters(ll),mm)) THEN ! found group
                    ngrp=ngrp+1
                    groupList(ioffg+ngrp)=j
                END IF
                mm=mm+1
                IF (mm >= bs) THEN
                    ll=ll+1
                    mm=mm-bs
                END IF
            END DO
            ! analyze rows (in overlap of group 'i 'and block 'ib')
            DO ir=max(irfrst,npgrp(i)),min(irlast,npgrp(i+1)-1)
                ! triangular part (parameter groups)
                DO j=1,ngrp
                    ! column loop
                    DO jc=max(ir,npgrp(groupList(ioffg+j))),npgrp(groupList(ioffg+j)+1)-1
                        ! block number
                        jb=(jc-1)/ibsize+1
                        blockCounter(ioffb+jb)=blockCounter(ioffb+jb)+1
                    END DO
                END DO
                ! rectangular part (parameters, constraints)
                noffi=ndimb+INT(ir-1,mpl)*INT(nac/bs+1,mpl)+1 ! row offset for J=1
                ll=noffi ! row offset
                mm=0
                DO j=1,nac
                    IF (btest(bitFieldCounters(ll),mm)) THEN ! found group
                        ! block number
                        jb=(nd+j-1)/ibsize+1
                        blockCounter(ioffb+jb)=blockCounter(ioffb+jb)+1
                    END IF
                    mm=mm+1
                    IF (mm >= bs) THEN
                        ll=ll+1
                        mm=mm-bs
                    END IF
                END DO
            END DO
        END DO
        ! end of 'row' block
        DO j=1,mb
            IF (blockCounter(ioffb+j) > 0) THEN
                nnzero=nnzero+blockCounter(ioffb+j)
                nblock=nblock+1
                nbkrow(ib)=nbkrow(ib)+1
            ENDIF
        END DO
    END DO
    !$OMP END PARALLEL DO

    ! 'empty' diagonal elements needed too
    DO ib=mbt+1,mb
        nnzero=nnzero+ibsize
        nblock=nblock+1
        nbkrow(ib)=1
    END DO

    ! cleanup
    CALL mpdealloc(groupList)
    CALL mpdealloc(blockCounter)
    CALL mpdealloc(rowBlocksToGroups)

#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_BY_NAME_END("pbsbits")
#endif
    !$POMP INST END(pbsbits)
    WRITE(*,*) ' '
    WRITE(*,*) 'PBSBITS: number of used elements', nnzero
    WRITE(*,1000) 'fraction of used elements', 200.0*REAL(nnzero,mps)/REAL(nd+nac,mps)/REAL(nd+nac+1,mps)
    WRITE(*,*) 'PBSBITS: block size', ibsize
    WRITE(*,*) 'PBSBITS: number of (used) blocks', nblock
    WRITE(*,1000) 'fraction of used storage ', 100.0*REAL(ibsize*ibsize+1,mps)*REAL(nblock,mps)/REAL(2*nnzero,mps)
1000 FORMAT(' PBSBITS: ',a,f7.2,' %')
    RETURN
END SUBROUTINE pbsbits

!> Analyze bit fields.
!!
!! Calculate BSR3 column list (for PARDISO)
!!
!! \param[in]      npgrp   parameter groups
!! \param[in]      ibsize  block size
!! \param[in]      nsparr  row offsets
!! \param[out]     nsparc  column list
!!
SUBROUTINE pblbits(npgrp,ibsize,nsparr,nsparc)
    USE mpbits
    USE mpdalc
    IMPLICIT NONE

    INTEGER(mpi), DIMENSION(:), INTENT(IN) :: npgrp
    INTEGER(mpi), INTENT(IN) :: ibsize
    INTEGER(mpl), DIMENSION(:), INTENT(IN) :: nsparr
    INTEGER(mpl), DIMENSION(:), INTENT(OUT) :: nsparc

    INTEGER(mpi) :: ichunk
    INTEGER(mpi) :: i
    INTEGER(mpi) :: ib
    INTEGER(mpi) :: igrpf
    INTEGER(mpi) :: igrpl
    INTEGER(mpi) :: iproc
    INTEGER(mpi) :: ioffb
    INTEGER(mpi) :: ioffg
    INTEGER(mpi) :: j
    INTEGER(mpi) :: mb
    INTEGER(mpi) :: mbt
    INTEGER(mpi) :: mm
    INTEGER(mpi) :: nd
    INTEGER(mpi) :: ngrp
    INTEGER(mpi) :: ir
    INTEGER(mpi) :: irfrst
    INTEGER(mpi) :: irlast
    INTEGER(mpi) :: jb
    INTEGER(mpi) :: jc
    INTEGER(mpl) :: kk
    INTEGER(mpl) :: length
    INTEGER(mpl) :: ll
    INTEGER(mpl) :: noffi
    INTEGER(mpl), PARAMETER :: two=2
    INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: rowBlocksToGroups
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: blockCounter
    INTEGER(mpi), DIMENSION(:), ALLOCATABLE :: groupList
    !$    INTEGER(mpi) :: OMP_GET_THREAD_NUM

    LOGICAL :: btest

    !$POMP INST BEGIN(pblbits)
#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_BY_NAME_BEGIN("pblbits", SCOREP_USER_REGION_TYPE_COMMON)
#endif
    nd=npgrp(n+1)-npgrp(1) ! number of diagonal elements
    mb=(nd+nac-1)/ibsize+1 ! max. number of blocks per row/column
    mbt=(nd-1)/ibsize+1 ! max. number of blocks in triangular part
    length=INT(mb,mpl)*INT(nthrd,mpl)
    CALL mpalloc(blockCounter,length,'PBBITS: block counter')
    length=INT(n,mpl)*INT(nthrd,mpl)
    CALL mpalloc(groupList,length,'PBBITS: group list')

    ! mapping row blocks to parameters groups
    length=INT(mbt,mpl)
    CALL mpalloc(rowBlocksToGroups,two,length,'mapping row blocks to par. groups (I)')
    rowBlocksToGroups(:,:)=0
    igrpf=1 ! first group of block
    igrpl=1 ! last  group of block
    ir=1    ! first row   of block
    DO i=1,mbt
        DO WHILE (igrpf < n .AND. npgrp(igrpf+1) <= ir)
            igrpf=igrpf+1
        END DO
        rowBlocksToGroups(1,i)=igrpf
        ir=ir+ibsize
        DO WHILE (igrpl < n .AND. npgrp(igrpl+1) < ir)
            igrpl=igrpl+1
        END DO
        rowBlocksToGroups(2,i)=igrpl
    END DO

    ll=0
    ichunk=MIN((mbt+nthrd-1)/nthrd/32+1,256)
    ! analyze bit fields ('upper half' (j>=i))
    ! parallelize row loop
    !$OMP  PARALLEL DO &
    !$OMP  PRIVATE(IPROC,IOFFB,IOFFG,I,NOFFI,KK,LL,MM,NGRP,J,IR,IRFRST,IRLAST) &
    !$OMP  SCHEDULE(DYNAMIC,ICHUNK)
    DO ib=1,mbt
        irfrst=ibsize*(ib-1)+1  ! first row in block
        irlast=ibsize*ib        ! last row in block
        iproc=0
        !$ IPROC=OMP_GET_THREAD_NUM()         ! thread number
        ioffb=iproc*mb
        ioffg=iproc*n
        blockCounter(ioffb+1:ioffb+mb)=0
        DO i=rowBlocksToGroups(1,ib),rowBlocksToGroups(2,ib)
            noffi=INT(i-1,mpl)*INT(2*n+2-i,mpl)/2 ! row offset for I=J
            ll=noffi/bs+i
            mm=0
            ngrp=0
            DO j=i,n         !  loop from diagonal element
                IF (btest(bitFieldCounters(ll),mm)) THEN ! found group
                    ngrp=ngrp+1
                    groupList(ioffg+ngrp)=j
                END IF
                mm=mm+1
                IF (mm >= bs) THEN
                    ll=ll+1
                    mm=mm-bs
                END IF
            END DO
            ! analyze rows (in overlap of group 'i 'and block 'ib')
            DO ir=max(irfrst,npgrp(i)),min(irlast,npgrp(i+1)-1)
                ! triangular part (parameter groups)
                DO j=1,ngrp
                    ! column loop
                    DO jc=max(ir,npgrp(groupList(ioffg+j))),npgrp(groupList(ioffg+j)+1)-1
                        ! block number
                        jb=(jc-1)/ibsize+1
                        blockCounter(ioffb+jb)=blockCounter(ioffb+jb)+1
                    END DO
                END DO
                ! rectangular part (parameters, constraints)
                noffi=ndimb+INT(ir-1,mpl)*INT(nac/bs+1,mpl)+1 ! row offset for J=1
                ll=noffi ! row offset
                mm=0
                DO j=1,nac
                    IF (btest(bitFieldCounters(ll),mm)) THEN ! found group
                        ! block number
                        jb=(nd+j-1)/ibsize+1
                        blockCounter(ioffb+jb)=blockCounter(ioffb+jb)+1
                    END IF
                    mm=mm+1
                    IF (mm >= bs) THEN
                        ll=ll+1
                        mm=mm-bs
                    END IF
                END DO
            END DO
        END DO
        ! end of 'row' block
        kk=nsparr(ib) ! offset for row block
        DO j=1,mb
            IF (blockCounter(ioffb+j) > 0) THEN
                ! store used column block indices
                nsparc(kk)=j
                kk=kk+1
            ENDIF
        END DO
    END DO
    !$OMP END PARALLEL DO

    ! 'empty' diagonal elements needed too
    DO ib=mbt+1,mb
        kk=nsparr(ib)
        nsparc(kk)=ib
    END DO

    ! cleanup
    CALL mpdealloc(groupList)
    CALL mpdealloc(blockCounter)
    CALL mpdealloc(rowBlocksToGroups)

#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_BY_NAME_END("pblbits")
#endif
    !$POMP INST END(pblbits)
    WRITE(*,*) ' '
    WRITE(*,*) 'PBLBITS: column list constructed ',nsparr(mb+1)-nsparr(1), ' words'
    CALL mpdealloc(bitFieldCounters)
    RETURN
END SUBROUTINE pblbits


!> Analyze bit fields.
!!
!! Calculate CSR3 row offsets (for PARDISO)
!!
!! \param[in]      npgrp   parameter groups
!! \param[out]     nsparr  row offsets (preset with number of Lagrange multipliers (at row+1))
!!
SUBROUTINE prbits(npgrp,nsparr)
    USE mpbits
    USE mpdalc
    IMPLICIT NONE

    INTEGER(mpi), DIMENSION(:), INTENT(IN) :: npgrp
    INTEGER(mpl), DIMENSION(:), INTENT(OUT) :: nsparr

    INTEGER(mpi) :: ichunk
    INTEGER(mpi) :: i
    INTEGER(mpi) :: j
    INTEGER(mpi) :: mm
    INTEGER(mpi) :: nd
    INTEGER(mpi) :: ir
    INTEGER(mpl) :: ll
    INTEGER(mpl) :: npar
    INTEGER(mpl) :: nparc
    INTEGER(mpl) :: ntot
    INTEGER(mpl) :: noffi

    LOGICAL :: btest

    nd=npgrp(n+1)-npgrp(1) ! number of diagonal elements

    ntot=0
    ll=0
    ichunk=MIN((n+nthrd-1)/nthrd/32+1,256)
    ! analyze bit fields ('upper half' (j>=i))
    ! parallelize row loop
    ! private copy of NTOT for each thread, combined at end, init with 0.
    !$OMP  PARALLEL DO &
    !$OMP  PRIVATE(I,NOFFI,LL,MM,J,IR,NPAR,NPARC) &
    !$OMP  SCHEDULE(DYNAMIC,ICHUNK)
    DO i=1,n
        noffi=INT(i-1,mpl)*INT(2*n+2-i,mpl)/2 ! row offset for I=J
        ll=noffi/bs+i
        mm=0
        npar=0
        ! triangular part (parameter groups)
        DO j=i,n         !  loop from diagonal element
            IF (btest(bitFieldCounters(ll),mm)) npar=npar+npgrp(j+1)-npgrp(j) ! number of parameters
            mm=mm+1
            IF (mm >= bs) THEN
                ll=ll+1
                mm=mm-bs
            END IF
        END DO

        ! save number of parameters for row(s)
        DO ir=npgrp(i),npgrp(i+1)-1
            ! rectangular part (parameters, constraints)
            noffi=ndimb+INT(ir-1,mpl)*INT(nac/bs+1,mpl)+1 ! row offset for J=1
            ll=noffi ! row offset
            mm=0
            nparc=0
            DO j=1,nac
                IF (btest(bitFieldCounters(ll),mm)) THEN ! found parameter/constraint combination
                    nparc=nparc+1
                END IF
                mm=mm+1
                IF (mm >= bs) THEN
                    ll=ll+1
                    mm=mm-bs
                END IF
            END DO
            nsparr(ir+1)=npar+nparc
            npar=npar-1
        END DO

    END DO
    !$OMP END PARALLEL DO

    ! sum up, fill row offsets
    nsparr(1)=1
    DO i=1,nd
        nsparr(i+1)=nsparr(i+1)+nsparr(i)
    END DO
    ! 'empty' diagonal elements needed too
    DO i=nd+1,nd+nac
        nsparr(i+1)=nsparr(i)+1
    END DO
    ntot=nsparr(nd+nac+1)-nsparr(1)

    WRITE(*,*) ' '
    WRITE(*,*) 'PRBITS: number of diagonal elements',nd+nac
    WRITE(*,*) 'PRBITS: number of used elements',ntot
    WRITE(*,1000) 'fraction of used elements', 200.0*REAL(ntot,mps)/REAL(nd+nac,mps)/REAL(nd+nac+1,mps)
1000 FORMAT(' PRBITS: ',a,f6.2,' %')
    RETURN
END SUBROUTINE prbits

!> Analyze bit fields.
!!
!! Calculate CSR3 column list (for PARDISO)
!!
!! \param[in]      npgrp   parameter groups
!! \param[in]      nsparr  row offsets
!! \param[out]     nsparc  column list
!!
SUBROUTINE pcbits(npgrp,nsparr,nsparc)
    USE mpbits
    USE mpdalc
    IMPLICIT NONE

    INTEGER(mpi), DIMENSION(:), INTENT(IN) :: npgrp
    INTEGER(mpl), DIMENSION(:), INTENT(IN) :: nsparr
    INTEGER(mpl), DIMENSION(:), INTENT(OUT) :: nsparc

    INTEGER(mpi) :: ichunk
    INTEGER(mpi) :: i
    INTEGER(mpi) :: j
    INTEGER(mpi) :: mm
    INTEGER(mpi) :: nd
    INTEGER(mpi) :: ic
    INTEGER(mpi) :: ir
    INTEGER(mpl) :: kk
    INTEGER(mpl) :: ll

    INTEGER(mpl) :: noffi
    INTEGER(mpl) :: noffr

    LOGICAL :: btest

    nd=npgrp(n+1)-npgrp(1) ! number of diagonal elements

    ll=0
    ichunk=MIN((n+nthrd-1)/nthrd/32+1,256)
    ! analyze bit fields ('upper half' (j>=i))
    ! parallelize row loop
    ! private copy of NTOT for each thread, combined at end, init with 0.
    !$OMP  PARALLEL DO &
    !$OMP  PRIVATE(I,NOFFI,LL,MM,J,IR,KK,IC) &
    !$OMP  SCHEDULE(DYNAMIC,ICHUNK)
    DO i=1,n
        noffi=INT(i-1,mpl)*INT(2*n+2-i,mpl)/2 ! row offset for I=J
        ! fill column list for row(s)
        DO ir=npgrp(i),npgrp(i+1)-1
            ! triangular part (parameter groups)
            ll=noffi/bs+i
            mm=0
            kk=nsparr(ir)
            DO j=i,n         !  loop from diagonal element
                IF (btest(bitFieldCounters(ll),mm)) THEN
                    ! 'all' columns in group
                    DO ic=max(ir,npgrp(j)),npgrp(j+1)-1
                        nsparc(kk)=ic
                        kk=kk+1
                    END DO
                ENDIF
                mm=mm+1
                IF (mm >= bs) THEN
                    ll=ll+1
                    mm=mm-bs
                END IF
            END DO
            ! rectangular part (parameters, constraints)
            noffr=ndimb+INT(ir-1,mpl)*INT(nac/bs+1,mpl)+1 ! row offset for J=1
            ll=noffr ! row offset
            mm=0
            DO j=1,nac
                IF (btest(bitFieldCounters(ll),mm)) THEN ! found parameter/constraint combination
                    nsparc(kk)=nd+j
                    kk=kk+1
                END IF
                mm=mm+1
                IF (mm >= bs) THEN
                    ll=ll+1
                    mm=mm-bs
                END IF
            END DO
        END DO

    END DO
    !$OMP END PARALLEL DO

    ! 'empty' diagonal elements needed too
    DO ir=nd+1,nd+nac
        kk=nsparr(ir)
        nsparc(kk)=ir
    END DO

    WRITE(*,*) ' '
    WRITE(*,*) 'PCBITS: column list constructed ',nsparr(nd+nac+1)-nsparr(1), ' words'
    CALL mpdealloc(bitFieldCounters)
    RETURN
END SUBROUTINE pcbits

!> Check sparsity of matrix.
!!
!! \param [in]     npgrp   parameter groups
!! \param [out]    ndims   (1): number of non-zero elements; (2): size of column lists;
!!                         (3/4): number of (double/single precision) off diagonal elements;
!!
SUBROUTINE ckbits(npgrp,ndims)
    USE mpbits
    IMPLICIT NONE

    INTEGER(mpi), DIMENSION(:), INTENT(IN) :: npgrp
    INTEGER(mpl), DIMENSION(4), INTENT(OUT) :: ndims

    INTEGER(mpi) :: nwcp(0:1)
    INTEGER(mpi) :: irgn(2)
    INTEGER(mpi) :: inr(2)
    INTEGER(mpl) :: ll
    INTEGER(mpl) :: noffi
    INTEGER(mpi) :: i
    INTEGER(mpi) :: j
    INTEGER(mpi) :: last
    INTEGER(mpi) :: lrgn
    INTEGER(mpi) :: next
    INTEGER(mpi) :: icp
    INTEGER(mpi) :: ib
    INTEGER(mpi) :: icount
    INTEGER(mpi) :: kbfw
    INTEGER(mpi) :: jp
    INTEGER(mpi) :: mm
    LOGICAL :: btest

    DO i=1,4
        ndims(i)=0
    END DO
    kbfw=1
    IF (ibfw > 1) kbfw=2
    ll=0

    DO i=1,n
        noffi=INT(i-1,mpl)*INT(i,mpl)*INT(ibfw,mpl)/2
        ll=noffi/bs+i
        mm=0
        inr(1)=0
        inr(2)=0
        irgn(1)=1
        irgn(2)=1
        last=0
        lrgn=0
        DO j=1,i
            icount=0
            next=0
            DO ib=0,ibfw-1
                IF (btest(bitFieldCounters(ll),mm)) icount=ibset(icount,ib)
                mm=mm+1
                IF (mm >= bs) THEN
                    ll=ll+1
                    mm=mm-bs
                END IF
            END DO

            IF (icount > 0) ndims(1)=ndims(1)+1
            !           keep pair ?
            IF (icount >= ireqpe) THEN
                next=1 ! double
                IF (icount <= isngpe) next=2 ! single
                inr(next)=inr(next)+npgrp(j+1)-npgrp(j) ! number of parameters
                IF (next /= last) THEN
                    irgn(next)=irgn(next)+1
                END IF
            END IF
            last=next
        END DO

        DO jp=1,kbfw
            IF (inr(jp) > 0) THEN
                icp=0
                nwcp(0)=inr(jp)         ! list of column indices (default)
                nwcp(1)=irgn(jp)*2      ! list of regions (group starts and offsets)
                ! compress row ?
                IF ((nwcp(1) < nwcp(0)).OR.iextnd > 0) THEN
                    icp=1
                END IF
                ! total space
                ndims(2)   =ndims(2)   +nwcp(icp)
                ndims(jp+2)=ndims(jp+2)+nwcp(0)*(npgrp(i+1)-npgrp(i))
            END IF
        END DO

    END DO

    RETURN
END SUBROUTINE ckbits

!> Create sparsity information.
!!
!! \param [in]     npgrp   parameter groups
!! \param[in ]    nsparr  row offsets
!! \param[out]    nsparc  column indices
!!
SUBROUTINE spbits(npgrp,nsparr,nsparc)               ! collect elements
    USE mpbits
    USE mpdalc
    IMPLICIT NONE

    INTEGER(mpi), DIMENSION(:), INTENT(IN) :: npgrp
    INTEGER(mpl), DIMENSION(:,:), INTENT(IN) :: nsparr
    INTEGER(mpi), DIMENSION(:), INTENT(OUT) :: nsparc

    INTEGER(mpl) :: kl
    INTEGER(mpl) :: l
    INTEGER(mpl) :: ll
    INTEGER(mpl) :: l1
    INTEGER(mpl) :: n1
    INTEGER(mpl) :: ndiff
    INTEGER(mpl) :: noffi
    INTEGER(mpl) :: noffj
    INTEGER(mpi) :: i
    INTEGER(mpi) :: j
    INTEGER(mpi) :: jb
    INTEGER(mpi) :: k
    INTEGER(mpi) :: m
    INTEGER(mpi) :: ichunk
    INTEGER(mpi) :: next
    INTEGER(mpi) :: last

    LOGICAL :: btest

    ichunk=MIN((n+nthrd-1)/nthrd/32+1,256)
    DO jb=0,nspc-1
        ! parallelize row loop
        !$OMP  PARALLEL DO &
        !$OMP  PRIVATE(I,N1,NOFFI,NOFFJ,L,M,KL,L1,NDIFF,LAST,LL,J,NEXT) &
        !$OMP  SCHEDULE(DYNAMIC,ICHUNK)
        DO i=1,n
            n1=i+jb*(n+1)
            noffi=INT(i-1,mpl)*INT(i,mpl)*INT(ibfw,mpl)/2
            l=noffi/bs+i
            m=jb
            kl=nsparr(1,n1)    ! pointer to row in NSPARC
            l1=nsparr(2,n1)    ! pointer to row in sparse matrix
            ndiff=(nsparr(1,n1+1)-kl)*(npgrp(i+1)-npgrp(i))-(nsparr(2,n1+1)-l1) ! 0 for no compression
            ll=l1
            last=0
            
            DO j=1,i         !  loop until diagonal element                
                next=0
                IF(bitFieldCounters(l) /= 0) THEN
                    IF(btest(bitFieldCounters(l),m)) THEN
                        IF (ndiff == 0) THEN
                            DO k=npgrp(j),npgrp(j+1)-1
                                kl=kl+1
                                nsparc(kl)=k ! column index
                            END DO
                        ELSE
                            next=1
                            IF (last == 0) THEN
                                kl=kl+1
                                nsparc(kl)=INT(ll-l1,mpi)
                                kl=kl+1
                                nsparc(kl)=j
                            END IF
                        END IF
                        ll=ll+(npgrp(j+1)-npgrp(j))
                    END IF
                END IF
                last=next
                m=m+nspc
                IF (m >= bs) THEN
                    m=m-bs
                    l=l+1
                END IF
            END DO

            ! extended storage ('2nd half' too) ? 
            IF (iextnd > 0) THEN
                noffj=(i-1)*nspc
                m=INT(MOD(noffj,INT(bs,mpl)),mpi)+jb
                last=0
                ! remaining columns
                DO j=i+1, n
                    ! index for pair (J,I)
                    noffi=INT(j-1,mpl)*INT(j,mpl)*INT(ibfw,mpl)/2 ! for I=1
                    l=noffi/bs+j+noffj/bs ! row offset + column offset
                    next=0
                    IF(btest(bitFieldCounters(l),m)) THEN
                        IF (ndiff == 0) THEN
                            DO k=npgrp(j),npgrp(j+1)-1
                                kl=kl+1
                                nsparc(kl)=k ! column index
                            END DO
                        ELSE
                            next=1
                            IF (last == 0) THEN
                                kl=kl+1
                                nsparc(kl)=INT(ll-l1,mpi)
                                kl=kl+1
                                nsparc(kl)=j
                            END IF
                        END IF
                        ll=ll+(npgrp(j+1)-npgrp(j))
                    END IF
                    last=next

                END DO
            END IF
            ! next offset, end marker
            IF (ndiff /= 0) THEN
                kl=kl+1
                nsparc(kl)=INT(ll-l1,mpi)
                kl=kl+1
                nsparc(kl)=n+1
            END IF

        END DO
    !$OMP END PARALLEL DO
    END DO 
        
    n1=(n+1)*nspc
    WRITE(*,*) ' '
    WRITE(*,*) 'SPBITS: sparse structure constructed ',nsparr(1,n1), ' words'
    WRITE(*,*) 'SPBITS: dimension parameter of matrix',nsparr(2,1)
    IF (ibfw <= 1) THEN
        WRITE(*,*) 'SPBITS: index of last used location',nsparr(2,n1)
    ELSE
        WRITE(*,*) 'SPBITS: index of last used double',nsparr(2,n1/2)
        WRITE(*,*) 'SPBITS: index of last used single',nsparr(2,n1)
    END IF
    CALL mpdealloc(bitFieldCounters)
    RETURN
END SUBROUTINE spbits


!> Clear (additional) bit map.
!!
!! \param [in]    in        matrix size
!
SUBROUTINE clbmap(in)
    USE mpbits
    USE mpdalc
    IMPLICIT NONE

    INTEGER(mpi), INTENT(IN) :: in

    INTEGER(mpl) :: noffd
    INTEGER(mpi) :: mb

    ! save input parameter 
    n2=in    
    ! bit field array size
    noffd=INT(n2,mpl)*INT(n2-1,mpl)/2
    ndimb2=noffd/bs+n2
    mb=INT(4.0E-6*REAL(ndimb2,mps),mpi)
    WRITE(*,*) ' '
    IF (mb > 0) THEN
        WRITE(*,*) 'CLBMAP: dimension of bit-map',ndimb2 , '(',mb,'MB)'
    ELSE
        WRITE(*,*) 'CLBMAP: dimension of bit-map',ndimb2 , '(< 1 MB)'
    END IF
    CALL mpalloc(bitMap,ndimb2,'INBMAP: bit storage')
    bitMap=0
    RETURN
END SUBROUTINE clbmap    

!> Fill bit map.
!!
!! \param [in]    im     first index
!! \param [in]    jm     second index
!!
SUBROUTINE inbmap(im,jm)        ! include element (I,J)
    USE mpbits
    IMPLICIT NONE

    INTEGER(mpi), INTENT(IN) :: im
    INTEGER(mpi), INTENT(IN) :: jm

    INTEGER(mpl) :: l
    INTEGER(mpi) :: i
    INTEGER(mpi) :: j
    INTEGER(mpi) :: noffj
    INTEGER(mpl) :: noffi
    INTEGER(mpi) :: m

    IF(im == jm) RETURN  ! diagonal
    j=MIN(im,jm)
    i=MAX(im,jm)
    IF(j <= 0) RETURN    ! out low
    IF(i > n2) RETURN    ! out high
    noffi=INT(i-1,mpl)*INT(i-2,mpl)/2 ! for J=1
    noffj=(j-1)
    l=noffi/bs+i+noffj/bs ! row offset + column offset
    !     add I instead of 1 to keep bit maps of different rows in different words (openMP !)
    m=MOD(noffj,bs)
    bitMap(l)=ibset(bitMap(l),m)
    RETURN
END SUBROUTINE inbmap 

!> Get pairs (statistic) from map.
!!
!! \param [in]     ngroup  number of parameter groups
!! \param [in]     npgrp   parameter groups
!! \param [out]    npair   number of paired parameters
!!
SUBROUTINE gpbmap(ngroup,npgrp,npair)
    USE mpbits
    IMPLICIT NONE

    INTEGER(mpi), INTENT(IN) :: ngroup
    INTEGER(mpi), DIMENSION(:,:), INTENT(IN) :: npgrp
    INTEGER(mpi), DIMENSION(:), INTENT(OUT) :: npair

    INTEGER(mpl) :: l
    INTEGER(mpl) :: noffi
    INTEGER(mpi) :: i
    INTEGER(mpi) :: j
    INTEGER(mpi) :: m
    LOGICAL :: btest

    npair(1:ngroup)=0
    l=0

    DO i=1,ngroup
        npair(i)=npair(i)+npgrp(2,i)-1 ! from own group
        noffi=INT(i-1,mpl)*INT(i-2,mpl)/2
        l=noffi/bs+i
        m=0
        DO j=1,i-1
            IF (btest(bitMap(l),m)) THEN
                ! from other group
                npair(i)=npair(i)+npgrp(2,j)
                npair(j)=npair(j)+npgrp(2,i)
            END IF
            m=m+1
            IF (m >= bs) THEN
                l=l+1
                m=m-bs
            END IF
        END DO
    END DO

    RETURN
END SUBROUTINE gpbmap

!> Get paired (parameter) groups from map.
!!
!! \param [in]     ipgrp   parameter group
!! \param [out]    npair   number of paired parameters
!! \param [in]     npgrp   paired parameter groups (for ipgrp)
!!
SUBROUTINE ggbmap(ipgrp,npair,npgrp)
    USE mpbits
    IMPLICIT NONE

    INTEGER(mpi), INTENT(IN) :: ipgrp
    INTEGER(mpi), INTENT(OUT) :: npair
    INTEGER(mpi), DIMENSION(:), INTENT(OUT) :: npgrp

    INTEGER(mpl) :: l
    INTEGER(mpl) :: noffi
    INTEGER(mpi) :: noffj
    INTEGER(mpi) :: i
    INTEGER(mpi) :: j
    LOGICAL :: btest

    npair=0

    i=ipgrp
    noffi=INT(i-1,mpl)*INT(i-2,mpl)/2 ! for J=1    
    l=noffi/bs+i! row offset
    !     add I instead of 1 to keep bit maps of different rows in different words (openMP !)
    DO j=1,ipgrp-1
        noffj=j-1
        IF (btest(bitMap(l+noffj/bs),MOD(noffj,bs))) THEN
            npair=npair+1
            npgrp(npair)=j
        END IF 
    END DO
 
    noffj=ipgrp-1
    DO i=ipgrp+1,n2
        noffi=INT(i-1,mpl)*INT(i-2,mpl)/2 ! for J=1    
        l=noffi/bs+i ! row offset   
        IF (btest(bitMap(l+noffj/bs),MOD(noffj,bs))) THEN
            npair=npair+1
            npgrp(npair)=i
        END IF 
    END DO

    RETURN
END SUBROUTINE ggbmap
