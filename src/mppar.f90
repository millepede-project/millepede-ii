!> \file
!! Parameter management

!! \author Maximilian Goblirsch-Kolb, DESY, 2025 (goblirsc@cern.ch)
!!
!! \copyright
!! Copyright (c) 2025 Deutsches Elektronen-Synchroton,
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

!> Functions to manage parameters
!!
MODULE mppar    
    USE mpmod
    USE mpdalc
    IMPLICIT NONE

    contains 
        
    
    !> largest prime number < N.
    !!
    !! \param [in] n N
    !! \return largest prime number < N

    INTEGER(mpi) FUNCTION iprime(n)
        USE mpdef

        IMPLICIT NONE
        INTEGER(mpi), INTENT(IN) :: n
        INTEGER(mpi) :: nprime
        INTEGER(mpi) :: nsqrt
        INTEGER(mpi) :: i
        !     ...
        SAVE
        nprime=n                               ! max number
        IF(MOD(nprime,2) == 0) nprime=nprime+1 ! ... odd number
        outer: DO
            nprime=nprime-2                        ! next lower odd number
            nsqrt=INT(SQRT(REAL(nprime,mps)),mpi)
            DO i=3,nsqrt,2                         !
                IF(i*(nprime/i) == nprime) CYCLE outer   ! test prime number
            END DO
            EXIT outer ! found
        END DO outer
        iprime=nprime
    END FUNCTION iprime
    !> Translate labels to indices (for global parameters).
    !!
    !! Functions indexForGlobalLabel and subroutine UPONE are
    !! used to collect items, i.e. labels, and to order and translate them.
    !!
    !! In the first phase items are collected and stored by calling
    !! <tt>IRES=indexForGlobalLabel(ITEM)</tt>.
    !!
    !! At the first entry the two sub-arrays "a" (globalParLabelIndex)
    !! and "b" (globalParHashTable) of length 2N
    !! are generated with a start length for N=128 entries.
    !! In array "a" two words are reserved for each item: (ITEM, count).
    !! The function indexForGlobalLabel(ITEM) returns the number of the item.
    !! At each entry the argument is compared with the already stored items,
    !! new items are stored. Search
    !! for entries is done using hash-indices, stored in sub-array "b".
    !! The initial hash-index is
    !!
    !!        j = 1 + mod(ITEM, n_prime) + N
    !!
    !! where n_prime is the largest prime number less than N.
    !! At each entry the count is increased by one. If N items are stored,
    !! the size of the sub-arrays is increased by calling
    !! <tt>CALL UPONE</tt>.
    !!
    !! \param[in] item  label
    !! \return index

    !> Make usable (sort items and redefine hash indices).
    SUBROUTINE useone
        USE mpmod

        IMPLICIT NONE
        INTEGER(mpi) :: i
        INTEGER(mpi) :: j
        INTEGER(mpi) :: k
        SAVE
        !     ...
        IF (globalParHeader(-1) > globalParHeader(-8)) THEN
            CALL sort22l(globalParLabelIndex,globalParLabelCounter,globalParHeader(-1)) ! sort items
            ! redefine hash
            globalParHashTable = 0
            outer: DO i=1,globalParHeader(-1)
                j=1+MOD(globalParLabelIndex(1,i),globalParHeader(-4))+globalParHeader(-0)
                inner: DO
                    k=j
                    j=globalParHashTable(k)
                    IF(j == 0) EXIT inner    ! unused hash code
                    IF(j == i) CYCLE outer ! found
                ENDDO inner
                globalParHashTable(k)=i
            END DO outer
            globalParHeader(-8)=globalParHeader(-1)
        END IF
    END SUBROUTINE useone                  ! make usable


    INTEGER(mpi) FUNCTION indexOfGlobalLabel(item)             ! translate 1-D identifier to nrs
        USE mpmod
        USE mpdalc

        IMPLICIT NONE
        INTEGER(mpi), INTENT(IN) :: item
        INTEGER(mpi) :: j
        INTEGER(mpi) :: k
        INTEGER(mpl) :: length
        INTEGER(mpl), PARAMETER :: four = 4

        indexOfGlobalLabel=0
        IF(item <= 0) RETURN
        IF(globalParHeader(-1) == 0) THEN
            length=128                   ! initial number
            CALL mpalloc(globalParLabelIndex,four,length,'indexOfGlobalLabel: label & index')
            CALL mpalloc(globalParLabelCounter,length,'indexOfGlobalLabel: counter') ! updated in pargrp
            CALL mpalloc(globalParHashTable,2*length,'indexOfGlobalLabel: hash pointer')
            globalParHashTable = 0
            globalParHeader(-0)=INT(length,mpi)       ! length of labels/indices
            globalParHeader(-1)=0                 ! number of stored items
            globalParHeader(-2)=0                 ! =0 during build-up
            globalParHeader(-3)=INT(length,mpi)       ! next number
            globalParHeader(-4)=iprime(globalParHeader(-0))    ! prime number
            globalParHeader(-5)=0                 ! number of overflows
            globalParHeader(-6)=0                 ! nr of variable parameters
            globalParHeader(-8)=0                 ! number of sorted items
        END IF
        outer: DO
            j=1+MOD(item,globalParHeader(-4))+globalParHeader(-0)
            inner: DO ! normal case: find item
                k=j
                j=globalParHashTable(k)
                IF(j == 0) EXIT inner    ! unused hash code
                IF(item == globalParLabelIndex(1,j)) EXIT outer ! found
            END DO inner
            ! not found
            IF(globalParHeader(-1) == globalParHeader(-0).OR.globalParHeader(-2) /= 0) THEN
                globalParHeader(-5)=globalParHeader(-5)+1 ! overflow
                j=0
                RETURN
            END IF
            globalParHeader(-1)=globalParHeader(-1)+1      ! increase number of elements
            globalParHeader(-3)=globalParHeader(-1)
            j=globalParHeader(-1)
            globalParHashTable(k)=j                ! hash index
            globalParLabelIndex(1,j)=item          ! add new item
            globalParLabelIndex(2,j)=0             ! reset index (for variable par.)
            globalParLabelIndex(3,j)=0             ! reset group info (first label)
            globalParLabelIndex(4,j)=0             ! reset group info (group index)
            globalParLabelCounter(j)=0             ! reset (long) counter
            IF(globalParHeader(-1) /= globalParHeader(-0)) EXIT outer
            ! update with larger dimension and redefine index
            globalParHeader(-3)=globalParHeader(-3)*2
            CALL upone
            IF (lvllog > 1) WRITE(lunlog,*) 'indexOfGlobalLabel: array increased to',  &
                globalParHeader(-3),' words'
        END DO outer

        ! counting now in pargrp
        !IF(globalParHeader(-2) == 0) THEN
        !    globalParLabelIndex(2,j)=globalParLabelIndex(2,j)+1 ! increase counter
        !    globalParHeader(-7)=globalParHeader(-7)+1
        !END IF
        indexOfGlobalLabel=j
    END FUNCTION indexOfGlobalLabel


    !> Update, redefine hash indices.
    SUBROUTINE upone
        USE mpmod
        USE mpdalc

        IMPLICIT NONE
        INTEGER(mpi) :: i
        INTEGER(mpi) :: j
        INTEGER(mpi) :: k
        INTEGER(mpi) :: nused
        LOGICAL :: finalUpdate
        INTEGER(mpl) :: oldLength
        INTEGER(mpl) :: newLength
        INTEGER(mpl), PARAMETER :: four = 4
        INTEGER(mpi), DIMENSION(:,:), ALLOCATABLE :: tempArr
        INTEGER(mpl), DIMENSION(:), ALLOCATABLE :: tempVec
        SAVE
        !     ...
        finalUpdate=(globalParHeader(-3) == globalParHeader(-1))
        IF(finalUpdate) THEN ! final (cleanup) call
            IF (globalParHeader(-1) > globalParHeader(-8)) THEN
                CALL sort22l(globalParLabelIndex,globalParLabelCounter,globalParHeader(-1)) ! sort items
                globalParHeader(-8)=globalParHeader(-1)
            END IF
        END IF
        ! save old LabelIndex
        nused = globalParHeader(-1)
        oldLength = globalParHeader(-0)
        CALL mpalloc(tempArr,four,oldLength,'indexOfGlobalLabel: temp array')
        tempArr(:,1:nused)=globalParLabelIndex(:,1:nused)
        CALL mpalloc(tempVec,oldLength,'indexOfGlobalLabel: temp vector')
        tempVec(1:nused)=globalParLabelCounter(1:nused)
        CALL mpdealloc(globalParLabelIndex)
        CALL mpdealloc(globalParLabelCounter)
        CALL mpdealloc(globalParHashTable)
        ! create new LabelIndex
        newLength = globalParHeader(-3)
        CALL mpalloc(globalParLabelIndex,four,newLength,'indexOfGlobalLabel: label & index')
        CALL mpalloc(globalParLabelCounter,newLength,'indexOfGlobalLabel: counter')
        CALL mpalloc(globalParHashTable,2*newLength,'indexOfGlobalLabel: hash pointer')
        globalParHashTable = 0
        globalParLabelIndex(:,1:nused) = tempArr(:,1:nused) ! copy back saved content
        globalParLabelCounter(1:nused) = tempVec(1:nused)   ! copy back saved content
        CALL mpdealloc(tempVec)
        CALL mpdealloc(tempArr)
        globalParHeader(-0)=INT(newLength,mpi)   ! length of labels/indices
        globalParHeader(-3)=globalParHeader(-1)
        globalParHeader(-4)=iprime(globalParHeader(-0))          ! prime number < LNDA
        ! redefine hash
        outer: DO i=1,globalParHeader(-1)
            j=1+MOD(globalParLabelIndex(1,i),globalParHeader(-4))+globalParHeader(-0)
            inner: DO
                k=j
                j=globalParHashTable(k)
                IF(j == 0) EXIT inner    ! unused hash code
                IF(j == i) CYCLE outer ! found
            ENDDO inner
            globalParHashTable(k)=i
        END DO outer
        IF(.NOT.finalUpdate) RETURN

        globalParHeader(-2)=1       ! set flag to inhibit further updates
        IF (lvllog > 1) THEN
            WRITE(lunlog,*) ' '
            WRITE(lunlog,*) 'indexOfGlobalLabel: array reduced to',newLength,' words'
            WRITE(lunlog,*) 'indexOfGlobalLabel:',globalParHeader(-1),' items stored.'
        END IF
    END SUBROUTINE upone                  !

    !> Prepare records.
    !!
    !! For global parameters replace label by index (<tt>indexForGlobalLabel</tt>).
    !!
    !! \param[in] mode <=0: build index table (indexForGlobalLabel) for global variables; \n
    !!                 >0: use index table, can be parallelized, optional scale errors 
    !!
    SUBROUTINE prepareRecords(mode)
        USE mpmod

        IMPLICIT NONE

        INTEGER(mpi), INTENT(IN) :: mode

        INTEGER(mpi) :: currRecord
        INTEGER(mpi) :: ichunk
        INTEGER(mpi) :: lastGlobal
        INTEGER(mpi) :: globalIndex
        INTEGER(mpi) :: j
        INTEGER(mpi) :: startLocal
        INTEGER(mpi) :: startGlobal
        INTEGER(mpi) :: startSpecial
        INTEGER(mpi) :: endOfentry
        INTEGER(mpi), PARAMETER :: maxbad = 100 ! max number of bad records with print out
        INTEGER(mpi) :: nbad
        INTEGER(mpi) :: nerr

        ! behaviour for *using* index table
        IF (mode > 0) THEN
#ifdef __PGIC__
            ! to prevent "PGF90-F-0000-Internal compiler error. Could not locate uplevel instance for stblock"
            ichunk=256
#else    
            ichunk=MIN((nbReadRecords+nOMPThreads-1)/nOMPThreads/32+1,256)
#endif
            ! parallelize record loop
            !$OMP  PARALLEL DO &
            !$OMP   DEFAULT(PRIVATE) &
            !$OMP   SHARED(nbReadRecords,readBufferPointer,readBufferDataI,readBufferDataD,ICHUNK,scaleErrors,errorScaleFactor) &
            !$OMP   SCHEDULE(DYNAMIC,ICHUNK)
            DO currRecord=1,nbReadRecords ! buffer for current record
                ! init iterator for last global derivative on measurement
                lastGlobal=readBufferPointer(currRecord)+1
                ! pick up the final entry of our record
                endOfentry=readBufferDataI(readBufferPointer(currRecord))
                DO ! loop over measurements
                    ! identify data blocks within entry
                    CALL decodeNextMeasurement(endOfentry,lastGlobal,startLocal,startGlobal,startSpecial)
                    ! startLocal: Points to the " 0    residual" entry 
                    ! startGlobal: Points to the " 0    sigma"  entry
                    ! lastGlobal: points to the final global derivative of this measurement
                    IF(startGlobal == 0) EXIT
                    ! loop over global labels
                    DO j=1,lastGlobal-startGlobal    
                        ! overwrite global labels with corresponding internal indices
                        readBufferDataI(startGlobal+j)=indexOfGlobalLabel( readBufferDataI(startGlobal+j) ) ! translate to index
                    END DO
                    ! scale error ?
                    IF (scaleErrors > 0) THEN
                        IF (startGlobal < lastGlobal) THEN
                            readBufferDataD(startGlobal) = readBufferDataD(startGlobal) * errorScaleFactor(1) ! 'global' measurement
                        ELSE
                            readBufferDataD(startGlobal) = readBufferDataD(startGlobal) * errorScaleFactor(2) ! 'local' measurement
                        END IF 
                    END IF    
                END DO
            END DO
            !$OMP  END PARALLEL DO
        END IF

        !$POMP INST BEGIN(peprep)
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_REGION_BY_NAME_BEGIN("UR_peprep", SCOREP_USER_REGION_TYPE_COMMON)
#endif
        ! behaviour for *building* index table
        IF (mode <= 0) THEN
            nbad=0
            DO currRecord=1,nbReadRecords ! buffer for current record
                CALL pechk(currRecord,nerr)
                IF(nerr > 0) THEN
                    nbad=nbad+1
                    IF(nbad >= maxbad) EXIT
                ELSE
                    lastGlobal=readBufferPointer(currRecord)+1
                    endOfentry=readBufferDataI(readBufferPointer(currRecord))
                    DO ! loop over measurements within record
                        CALL decodeNextMeasurement(endOfentry,lastGlobal,startLocal,startGlobal,startSpecial)
                        ! this changes lastGlobal to point to the last global derivative
                        ! of the current measurement.
                        ! It will give us the next measurement within the record.
                        IF(startGlobal == 0) EXIT
                        nbEqRead=nbEqRead+1

                        ! no global derivatives on this measurement
                        IF(startGlobal == lastGlobal) CYCLE  
                        
                        nbEqWithGlobPar=nbEqWithGlobPar+1
                        nbGlobDeriv=nbGlobDeriv+(lastGlobal-startGlobal)
                        DO j=1,lastGlobal-startGlobal
                            ! this will also register the labels in our map,
                            ! the first time they are seen. 
                            globalIndex=indexOfGlobalLabel( readBufferDataI(startGlobal+j) ) ! generate index
                        END DO
                    END DO
                END IF
            END DO
            IF(nbad > 0) THEN
                CALL peend(20,'Aborted, bad binary records')
                STOP 'PEREAD: stopping due to bad records'
            END IF
        END IF
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_REGION_BY_NAME_END("UR_peprep")
#endif
        !$POMP INST END(peprep)

    END SUBROUTINE prepareRecords

end module mppar
