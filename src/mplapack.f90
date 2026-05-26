!> \file
!! Lapack-specific subroutines.
!!
!! \author Claus Kleinwort, DESY (original implementation)
!! \author Maximilian Goblirsch-Kolb DESY (refactoring)
!!
!! \copyright
!! Copyright (c) 2009 - 2025 Deutsches Elektronen-Synchroton,
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

#ifdef PARDISO
    INCLUDE 'mkl_pardiso.f90'
#endif

MODULE mplapack    
    USE mpmod
    USE mpdalc
    IMPLICIT NONE

#ifdef LAPACK64
    contains 
    !> Solution by factorization.
    !!
    !! Using LAPACK routines, packed storage (DyPTRF, DyPTRS)
    !! Solve A*x=b with A=LL^t positive definite (Cholesky, elimination)
    !! or with A=LDL^t indefinite (Bunch-Kaufman, Lagrange multipliers)

    SUBROUTINE mdptrf

        IMPLICIT NONE
        INTEGER(mpi) :: i
        INTEGER(mpi) :: ib
        INTEGER(mpi) :: icoff
        INTEGER(mpi) :: ipoff
        INTEGER(mpi) :: j
        INTEGER(mpi) :: lun
        INTEGER(mpi) :: ncon
        INTEGER(mpi) :: nfit
        INTEGER(mpi) :: npar
        INTEGER(mpl) :: imoff
        INTEGER(mpl) :: ioff1
        INTEGER(mpi) :: infolp
        REAL(mpd) :: matij
            
        EXTERNAL avprds

        SAVE
        !     ...
        lun=lunlog                       ! log file

        IF(icalcm == 1) THEN
            IF(ilperr == 1) THEN
                ! save diagonal (for global correlation)
                DO i=1,nAllActivePar
                    workspaceDiag(i)=matij(i,i)
                END DO
            END IF
            ! use elimination for constraints ?
            IF(nFitPar < nVarGlobalPar) THEN
                ! monitor progress
                IF(monpg1 > 0) THEN
                    WRITE(lunlog,*) 'Shrinkage of global matrix (A->Q^t*A*Q)'
                    CALL monini(lunlog,monpg1,monpg2)
                END IF
                CALL qlssq(avprds,globalMatD,size(globalMatD,kind=mpl),globalRowOffsets,.true.) ! Q^t*A*Q
                IF(monpg1 > 0) CALL monend()
            END IF
        END IF

        ! loop over blocks (multiple blocks only with elimination !)
        DO ib=1,nParBlocks
            ipoff=matParBlockOffsets(1,ib)          ! parameter offset for block
            npar=matParBlockOffsets(1,ib+1)-ipoff   ! number of parameters in block
            icoff=vecParBlockConOffsets(ib)         ! constraint offset for block
            ncon=vecParBlockConOffsets(ib+1)-icoff  ! number of constraints in block
            imoff=globalRowOffsets(ipoff+1)+ipoff   ! block offset in global matrix
            nfit=npar+ncon; IF (icelim > 0) nfit=npar-ncon ! number of fit parameters in block
            ! use elimination for constraints ?
            IF(nfit < npar) THEN
                CALL qlsetb(ib)
                ! solve L^t*y=d by backward substitution
                CALL qlbsub(vecConsResiduals(icoff+1:),vecConsSolution)
                ! transform, reduce rhs
                CALL qlmlq(globalCorrections(ipoff+1:),1,.true.) ! Q^t*b
                ! correction from eliminated part
                DO i=1,nfit
                    DO j=1,ncon
                        ioff1=globalRowOffsets(nfit+j+ipoff)+i+ipoff ! local (nfit+j,i)
                        globalCorrections(i+ipoff)=globalCorrections(i+ipoff)-globalMatD(ioff1)*vecConsSolution(j)
                    END DO
                END DO
            END IF

            IF(icalcm == 1) THEN
                ! multipliers?
                IF (nfit > npar) THEN
                    ! monitor progress
                    IF(monpg1 > 0) THEN
                        WRITE(lunlog,*) 'Factorization of global matrix (A->L*D*L^t)'
                        CALL monini(lunlog,monpg1,monpg2)
                    END IF
                    !$POMP INST BEGIN(dsptrf)
#ifdef SCOREP_USER_ENABLE
                    SCOREP_USER_REGION_BY_NAME_BEGIN("UR_dsptrf", SCOREP_USER_REGION_TYPE_COMMON)
#endif
                    CALL dsptrf('U',INT(nfit,mpl),globalMatD(imoff+1:),lapackIPIV(ipoff+1:),infolp)
#ifdef SCOREP_USER_ENABLE
                    SCOREP_USER_REGION_BY_NAME_END("UR_dsptrf")
#endif
                    !$POMP INST END(dsptrf)
                    IF(monpg1 > 0) CALL monend()
                ELSE
                    ! monitor progress
                    IF(monpg1 > 0) THEN
                        WRITE(lunlog,*) 'Factorization of global matrix (A->L*L^t)'
                        CALL monini(lunlog,monpg1,monpg2)
                    END IF
                    !$POMP INST BEGIN(dpptrf)
#ifdef SCOREP_USER_ENABLE
                    SCOREP_USER_REGION_BY_NAME_BEGIN("UR_dpptrf", SCOREP_USER_REGION_TYPE_COMMON)
#endif
                    CALL dpptrf('U',INT(nfit,mpl),globalMatD(imoff+1:),infolp)
#ifdef SCOREP_USER_ENABLE
                    SCOREP_USER_REGION_BY_NAME_END("UR_dpptrf")
#endif
                    !$POMP INST END(dpptrf)
                    IF(monpg1 > 0) CALL monend()
                ENDIF
                ! check result
                IF(infolp==0) THEN
                    IF(nParBlocks == 1) THEN
                        WRITE(lun,*) 'No rank defect of the symmetric matrix'
                    ELSE
                        WRITE(lun,*) 'No rank defect of the symmetric block', ib, ' of size', npar
                    END IF
                ELSE
                    ndefec=ndefec+1 ! (lower limit of) rank defect
                    WRITE(*,*)   'Warning: factorization of the symmetric',nfit,  &
                        '-by-',nfit,' failed at index ', infolp
                    WRITE(lun,*)   'Warning: factorization of the symmetric',nfit,  &
                        '-by-',nfit,' failed at index ', infolp
                    CALL peend(29,'Aborted, factorization of global matrix failed')
                    STOP 'mdptrf: bad matrix'
                END IF
            END IF
            ! backward/forward substitution
            ! multipliers?
            IF (nfit > npar) THEN
                CALL dsptrs('U',INT(nfit,mpl),1_mpl,globalMatD(imoff+1:),lapackIPIV(ipoff+1:),&
                    globalCorrections(ipoff+1:),INT(nfit,mpl),infolp)
                IF(infolp /= 0) PRINT *, ' DSPTRS failed: ', infolp
            ELSE
                CALL dpptrs('U',INT(nfit,mpl),1_mpl,globalMatD(imoff+1:),&
                    globalCorrections(ipoff+1:),INT(nfit,mpl),infolp)
                IF(infolp /= 0) PRINT *, ' DPPTRS failed: ', infolp
            ENDIF

            !use elimination for constraints ?
            IF(nfit < npar) THEN
                ! extend, transform back solution
                globalCorrections(nfit+1+ipoff:npar+ipoff)=vecConsSolution(1:ncon)
                CALL qlmlq(globalCorrections(ipoff+1:),1,.false.) ! Q*x
            END IF
        END DO
        
    END SUBROUTINE mdptrf

    !> Solution by factorization.
    !!
    !! Using LAPACK routines, unpacked storage (DPOTRF/DSYTRF, DPOTRS/DSYTRS)
    !! Solve A*x=b with A=LL^t positive definite (Cholesky, elimination)
    !! or with A=LDL^t indefinite (Bunch-Kaufman, Lagrange multipliers)

    SUBROUTINE mdutrf

        IMPLICIT NONE
        INTEGER(mpi) :: i
        INTEGER(mpi) :: ib
        INTEGER(mpi) :: icoff
        INTEGER(mpi) :: ipoff
        INTEGER(mpi) :: j
        INTEGER(mpi) :: lun
        INTEGER(mpi) :: ncon
        INTEGER(mpi) :: nfit
        INTEGER(mpi) :: npar
        INTEGER(mpl) :: imoff
        INTEGER(mpl) :: ioff1
        INTEGER(mpl) :: iloff
        INTEGER(mpi) :: infolp
        
        REAL(mpd) :: matij
            
        EXTERNAL avprds

        SAVE
        !     ...
        lun=lunlog                       ! log file

        IF(icalcm == 1) THEN
            IF(ilperr == 1) THEN
                ! save diagonal (for global correlation)
                DO i=1,nAllActivePar
                    workspaceDiag(i)=matij(i,i)
                END DO
            END IF
            ! use elimination for constraints ?
            IF(nFitPar < nVarGlobalPar) THEN
                ! monitor progress
                IF(monpg1 > 0) THEN
                    WRITE(lunlog,*) 'Shrinkage of global matrix (A->Q^t*A*Q)'
                    CALL monini(lunlog,monpg1,monpg2)
                END IF
                IF (icelim > 1) THEN
                    CALL lpavat(.true.)
                ELSE
                    CALL qlssq(avprds,globalMatD,size(globalMatD,kind=mpl),globalRowOffsets,.true.) ! Q^t*A*Q            
                END IF    
                IF(monpg1 > 0) CALL monend()
            END IF
        END IF

        ! loop over blocks (multiple blocks only with elimination !)
        iloff=0 ! offset of L in lapackQL
        DO ib=1,nParBlocks
            ipoff=matParBlockOffsets(1,ib)          ! parameter offset for block
            npar=matParBlockOffsets(1,ib+1)-ipoff   ! number of parameters in block
            icoff=vecParBlockConOffsets(ib)         ! constraint offset for block
            ncon=vecParBlockConOffsets(ib+1)-icoff  ! number of constraints in block
            imoff=globalRowOffsets(ipoff+1)+ipoff   ! block offset in global matrix
            nfit=npar+ncon; IF (icelim > 0) nfit=npar-ncon ! number of fit parameters in block
            ! use elimination for constraints ?
            IF(nfit < npar) THEN
                IF (icelim > 1) THEN
                    ! solve L^t*y=d by backward substitution
                    vecConsSolution(1:ncon)=vecConsResiduals(icoff+1:icoff+ncon)
                    CALL dtrtrs('L','T','N',INT(ncon,mpl),1_mpl,lapackQL(iloff+npar-ncon+1:),INT(npar,mpl),&
                        vecConsSolution,INT(ncon,mpl),infolp)
                    IF(infolp /= 0) PRINT *, ' DTRTRS failed: ', infolp
                    ! transform, reduce rhs, Q^t*b
                    CALL dormql('L','T',INT(npar,mpl),1_mpl,INT(ncon,mpl),lapackQL(iloff+1:),INT(npar,mpl),&
                        lapackTAU(icoff+1:),globalCorrections(ipoff+1:),INT(npar,mpl),lapackWORK,lplwrk,infolp)
                    IF(infolp /= 0) PRINT *, ' DORMQL failed: ', infolp
                ELSE
                    CALL qlsetb(ib)
                    ! solve L^t*y=d by backward substitution
                    CALL qlbsub(vecConsResiduals(icoff+1:),vecConsSolution)
                    ! transform, reduce rhs
                    CALL qlmlq(globalCorrections(ipoff+1:),1,.true.) ! Q^t*b   
                END IF
                ! correction from eliminated part
                DO i=1,nfit
                    DO j=1,ncon
                        ioff1=globalRowOffsets(nfit+j+ipoff)+i+ipoff ! local (nfit+j,i)
                        globalCorrections(i+ipoff)=globalCorrections(i+ipoff)-globalMatD(ioff1)*vecConsSolution(j)
                    END DO
                END DO
            END IF

            IF(icalcm == 1) THEN
                ! multipliers?
                IF (nfit > npar) THEN
                    ! monitor progress
                    IF(monpg1 > 0) THEN
                        WRITE(lunlog,*) 'Factorization of global matrix (A->L*D*L^t)'
                        CALL monini(lunlog,monpg1,monpg2)
                    END IF
                    !$POMP INST BEGIN(dsytrf)
#ifdef SCOREP_USER_ENABLE
                    SCOREP_USER_REGION_BY_NAME_BEGIN("UR_dsytrf", SCOREP_USER_REGION_TYPE_COMMON)
#endif
                    CALL dsytrf('U',INT(nfit,mpl),globalMatD(imoff+1:),INT(nfit,mpl),&
                        lapackIPIV(ipoff+1:),lapackWORK,lplwrk,infolp)
#ifdef SCOREP_USER_ENABLE
                    SCOREP_USER_REGION_BY_NAME_END("UR_dsytrf")
#endif
                    !$POMP INST END(dsytrf)
                    IF(monpg1 > 0) CALL monend()
                ELSE
                    ! monitor progress
                    IF(monpg1 > 0) THEN
                        WRITE(lunlog,*) 'Factorization of global matrix (A->L*L^t)'
                        CALL monini(lunlog,monpg1,monpg2)
                    END IF
                    !$POMP INST BEGIN(dpotrf)
#ifdef SCOREP_USER_ENABLE
                    SCOREP_USER_REGION_BY_NAME_BEGIN("UR_dpotrf", SCOREP_USER_REGION_TYPE_COMMON)
#endif
                    CALL dpotrf('U',INT(nfit,mpl),globalMatD(imoff+1:),INT(npar,mpl),infolp)
#ifdef SCOREP_USER_ENABLE
                    SCOREP_USER_REGION_BY_NAME_END("UR_dpotrf")
#endif
                    !$POMP INST END(dpotrf)
                    IF(monpg1 > 0) CALL monend()
                ENDIF
                ! check result
                IF(infolp==0) THEN
                    IF(nParBlocks == 1) THEN
                        WRITE(lun,*) 'No rank defect of the symmetric matrix'
                    ELSE
                        WRITE(lun,*) 'No rank defect of the symmetric block', ib, ' of size', npar
                    END IF
                ELSE
                    ndefec=ndefec+1 ! (lower limit of) rank defect
                    WRITE(*,*)   'Warning: factorization of the symmetric',nfit,  &
                        '-by-',nfit,' failed at index ', infolp
                    WRITE(lun,*)   'Warning: factorization of the symmetric',nfit,  &
                        '-by-',nfit,' failed at index ', infolp
                    CALL peend(29,'Aborted, factorization of global matrix failed')
                    STOP 'mdutrf: bad matrix'
                END IF
            END IF
            ! backward/forward substitution
            ! multipliers?
            IF (nfit > npar) THEN
                CALL dsytrs('U',INT(nfit,mpl),1_mpl,globalMatD(imoff+1:),INT(nfit,mpl),&
                    lapackIPIV(ipoff+1:),globalCorrections(ipoff+1:),INT(nfit,mpl),infolp)
                IF(infolp /= 0) PRINT *, ' DSYTRS failed: ', infolp
            ELSE
                CALL dpotrs('U',INT(nfit,mpl),1_mpl,globalMatD(imoff+1:),INT(npar,mpl),&
                    globalCorrections(ipoff+1:),INT(npar,mpl),infolp)
                IF(infolp /= 0) PRINT *, ' DPOTRS failed: ', infolp
            ENDIF

            !use elimination for constraints ?
            IF(nfit < npar) THEN
                IF (icelim > 1) THEN
                    ! correction from eliminated part
                    globalCorrections(nfit+1+ipoff:npar+ipoff)=vecConsSolution(1:ncon)
                    ! extend, transform back solution, Q*x
                    CALL dormql('L','N',INT(npar,mpl),1_mpl,INT(ncon,mpl),lapackQL(iloff+1:),INT(npar,mpl),&
                        lapackTAU(icoff+1:),globalCorrections(ipoff+1:),INT(npar,mpl),lapackWORK,lplwrk,infolp)
                    IF(infolp /= 0) PRINT *, ' DORMQL failed: ', infolp
                ELSE
                    ! extend, transform back solution
                    globalCorrections(nfit+1+ipoff:npar+ipoff)=vecConsSolution(1:ncon)
                    CALL qlmlq(globalCorrections(ipoff+1:),1,.false.) ! Q*x            
                END IF    
            END IF
            iloff=iloff+INT(npar,mpl)*INT(ncon,mpl)
        END DO
        
    END SUBROUTINE mdutrf

    !> QL decomposition.
    !!
    !! QL decomposition of constraints matrix for solution by elimination
    !! for unpacked storage using LAPACK.
    !! Optionally split into disjoint blocks.
    !!
    !! \param [in]     a  packed constraint matrix
    !! \param [out]    emin  eigenvalue with smallest absolute value
    !! \param [out]    emax  eigenvalue with largest absolute value
    !!
    SUBROUTINE lpqldec(a,emin,emax)

        IMPLICIT NONE
        INTEGER(mpi) :: ib
        INTEGER(mpi) :: icb
        INTEGER(mpi) :: icboff
        INTEGER(mpi) :: icblst
        INTEGER(mpi) :: icoff
        INTEGER(mpi) :: icfrst
        INTEGER(mpi) :: iclast
        INTEGER(mpi) :: ipfrst
        INTEGER(mpi) :: iplast
        INTEGER(mpi) :: ipoff
        INTEGER(mpi) :: i
        INTEGER(mpi) :: j
        INTEGER(mpi) :: ncon
        INTEGER(mpi) :: npar
        INTEGER(mpi) :: npb
        INTEGER(mpl) :: imoff
        INTEGER(mpl) :: iloff
        INTEGER(mpi) :: infolp
        INTEGER :: nbopt, ILAENV
        
        REAL(mpd), INTENT(IN) :: a(mszcon)
        REAL(mpd), INTENT(OUT)         :: emin
        REAL(mpd), INTENT(OUT)         :: emax
        SAVE

        PRINT *
        ! loop over blocks (multiple blocks only with elimination !)
        iloff=0 ! size of unpacked constraint matrix
        DO ib=1,nParBlocks
            ipoff=matParBlockOffsets(1,ib)          ! parameter offset for block
            npar=matParBlockOffsets(1,ib+1)-ipoff   ! number of parameters in block
            icoff=vecParBlockConOffsets(ib)         ! constraint offset for block
            ncon=vecParBlockConOffsets(ib+1)-icoff  ! number of constraints in block
            iloff=iloff+INT(npar,mpl)*INT(ncon,mpl)
        END DO
        ! allocate
        CALL mpalloc(lapackQL, iloff, 'LAPACK QL (QL decomp.) ')
        lapackQL=0.
        iloff=nConstraints
        CALL mpalloc(lapackTAU, iloff, 'LAPACK TAU (QL decomp.) ')
        ! fill
        iloff=0 ! offset of unpacked constraint matrix block
        imoff=0 ! offset of packed constraint matrix block
        DO ib=1,nParBlocks
            ipoff=matParBlockOffsets(1,ib)          ! parameter offset for block
            npar=matParBlockOffsets(1,ib+1)-ipoff   ! number of parameters in block
            icoff=vecParBlockConOffsets(ib)         ! constraint offset for block
            ncon=vecParBlockConOffsets(ib+1)-icoff  ! number of constraints in block
            IF(ncon <= 0) CYCLE
            ! block with constraints
            icboff=matParBlockOffsets(2,ib) ! constraint block offset
            icblst=matParBlockOffsets(2,ib+1) ! constraint block offset
            DO icb=icboff+1,icboff+icblst
                icfrst=matConsBlocks(1,icb)       ! first constraint in block
                iclast=matConsBlocks(1,icb+1)-1   ! last constraint in block
                DO j=icfrst,iclast
                    ipfrst=matConsRanges(3,j)-ipoff ! first (rel.) parameter
                    iplast=matConsRanges(4,j)-ipoff ! last  (rel.) parameters
                    npb=iplast-ipfrst+1
                    lapackQL(iloff+ipfrst:iloff+iplast)=a(imoff+1:imoff+npb)
                    imoff=imoff+npb
                    iloff=iloff+npar
                END DO
            END DO
        END DO
        ! decompose
        iloff=0 ! offset of unpacked constraint matrix block
        emax=-1.
        emin=1.
        DO ib=1,nParBlocks
            ipoff=matParBlockOffsets(1,ib)          ! parameter offset for block
            npar=matParBlockOffsets(1,ib+1)-ipoff   ! number of parameters in block
            icoff=vecParBlockConOffsets(ib)         ! constraint offset for block
            ncon=vecParBlockConOffsets(ib+1)-icoff  ! number of constraints in block
            IF(ncon <= 0) CYCLE
            ! block with constraints
            nbopt = ilaenv( 1_mpl, 'DGEQLF', '', int(npar,mpl), int(ncon,mpl), int(npar,mpl), -1_mpl ) ! optimal block size
            PRINT *, 'LAPACK optimal block size for DGEQLF:', nbopt
            lplwrk=int(ncon,mpl)*int(nbopt,mpl)
            CALL mpalloc(lapackWORK, lplwrk,'LAPACK WORK array (d)')
            !$POMP INST BEGIN(dgeqlf)
#ifdef SCOREP_USER_ENABLE
            SCOREP_USER_REGION_BY_NAME_BEGIN("UR_dgeqlf", SCOREP_USER_REGION_TYPE_COMMON)
#endif
            CALL dgeqlf(INT(npar,mpl),INT(ncon,mpl),lapackQL(iloff+1:),INT(npar,mpl),&
                lapackTAU(icoff+1:),lapackWORK,lplwrk,infolp)
            IF(infolp /= 0) PRINT *, ' DGEQLF failed: ', infolp
#ifdef SCOREP_USER_ENABLE
            SCOREP_USER_REGION_BY_NAME_END("UR_dgeqlf")
#endif
            !$POMP INST END(dgeqlf)
            CALL mpdealloc(lapackwork)
            iloff=iloff+INT(npar,mpl)*INT(ncon,mpl)
            ! get min/max diaginal element of L
            imoff=iloff
            IF(emax < emin) THEN
                emax=lapackQL(imoff)
                emin=emax
            END IF
            DO i=1,ncon
                IF (ABS(emax) < ABS(lapackQL(imoff))) emax=lapackQL(imoff)
                IF (ABS(emin) > ABS(lapackQL(imoff))) emin=lapackQL(imoff)
                imoff=imoff-npar-1
            END DO
        END DO
        PRINT *
    END SUBROUTINE lpqldec

    !> Similarity transformation by Q(t).
    !!
    !! Similarity transformation for global matrix by Q from QL decomposition
    !! for unpacked storage using LAPACK.
    !!
    !! Global matrix A is replaced by Q*A*Q^t (t=false) or Q^t*A*Q (t=true)
    !!
    !! \param [in]     t        use transposed of Q
    !!
    SUBROUTINE lpavat(t)

        IMPLICIT NONE
        INTEGER(mpi) :: i
        INTEGER(mpi) :: ib
        INTEGER(mpi) :: icoff
        INTEGER(mpi) :: ipoff
        INTEGER(mpi) :: j
        INTEGER(mpi) :: ncon
        INTEGER(mpi) :: npar
        INTEGER(mpl) :: imoff
        INTEGER(mpl) :: iloff
        INTEGER(mpi) :: infolp
        CHARACTER (LEN=1) :: transr, transl
                    
        LOGICAL, INTENT(IN) :: t
        SAVE

        IF (t) THEN     ! Q^t*A*Q
            transr='N'
            transl='T'
        ELSE            ! Q*A*Q^t
            transr='T'
            transl='N'
        ENDIF
            
        ! loop over blocks (multiple blocks only with elimination !)
        iloff=0 ! offset of L in lapackQL
        DO ib=1,nParBlocks
            ipoff=matParBlockOffsets(1,ib)          ! parameter offset for block
            npar=matParBlockOffsets(1,ib+1)-ipoff   ! number of parameters in block
            icoff=vecParBlockConOffsets(ib)         ! constraint offset for block
            ncon=vecParBlockConOffsets(ib+1)-icoff  ! number of constraints in block
            imoff=globalRowOffsets(ipoff+1)+ipoff   ! block offset in global matrix
            IF(ncon <= 0 ) CYCLE
            
            !$POMP INST BEGIN(dormql)
#ifdef SCOREP_USER_ENABLE
            SCOREP_USER_REGION_BY_NAME_BEGIN("UR_dormql", SCOREP_USER_REGION_TYPE_COMMON)
#endif
            ! expand matrix (copy lower to upper triangle)
            ! parallelize row loop
            ! slot of 32 'I' for next idle thread
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(J) &
            !$OMP SCHEDULE(DYNAMIC,32)
            DO i=ipoff+1,ipoff+npar
                DO j=ipoff+1,i-1
                    globalMatD(globalRowOffsets(j)+i)=globalMatD(globalRowOffsets(i)+j)
                ENDDO
            ENDDO
            ! A*Q
            CALL dormql('R',transr,int(npar,mpl),int(npar,mpl),int(ncon,mpl),lapackQL(iloff+1:),&
                INT(npar,mpl),lapackTAU(icoff+1:),globalMatD(imoff+1:),int(npar,mpl),&
                lapackWORK,lplwrk,infolp)
            IF(infolp /= 0) PRINT *, ' DORMQL failed: ', infolp
            ! Q^t*(A*Q)
            CALL dormql('L',transl,int(npar,mpl),int(npar,mpl),int(ncon,mpl),lapackQL(iloff+1:),&
                INT(npar,mpl),lapackTAU(icoff+1:),globalMatD(imoff+1:),int(npar,mpl),&
                lapackWORK,lplwrk,infolp)
            IF(infolp /= 0) PRINT *, ' DORMQL failed: ', infolp
#ifdef SCOREP_USER_ENABLE
            SCOREP_USER_REGION_BY_NAME_END("UR_dormql")
#endif
            !$POMP INST END(dormql)

            iloff=iloff+INT(npar,mpl)*INT(ncon,mpl)
        END DO

    END SUBROUTINE lpavat

#ifdef PARDISO
    !===============================================================================
    ! Copyright 2004-2022 Intel Corporation.
    !
    ! This software and the related documents are Intel copyrighted  materials,  and
    ! your use of  them is  governed by the  express license  under which  they were
    ! provided to you (License).  Unless the License provides otherwise, you may not
    ! use, modify, copy, publish, distribute,  disclose or transmit this software or
    ! the related documents without Intel''s prior written permission.
    !
    ! This software and the related documents  are provided as  is,  with no express
    ! or implied  warranties,  other  than those  that are  expressly stated  in the
    ! License.
    !===============================================================================
    !
    !   Content : Intel(R) oneAPI Math Kernel Library (oneMKL) PARDISO Fortran-90
    !             use case
    !
    !*******************************************************************************

    !> Solution with Intel(R) oneAPI Math Kernel Library (oneMKL) PARDISO
    !!
    !! Sparse matrix
    !!
    SUBROUTINE mspardiso
        USE mkl_pardiso
        IMPLICIT NONE

        !.. Internal solver memory pointer
        TYPE(MKL_PARDISO_HANDLE) :: pt(64)    ! Handle to internal data structure
        !.. All other variables
        INTEGER(mpl), PARAMETER :: maxfct =1  ! Max. number of factors with identical sparsity structure kept in memory
        INTEGER(mpl), PARAMETER :: mnum = 1   ! Actual factor to use
        INTEGER(mpl), PARAMETER :: nrhs = 1   ! Number of right hand sides

        INTEGER(mpl) :: mtype                  ! Matrix type (symmetric, pos. def.: 2, indef.: -2)
        INTEGER(mpl) :: phase                  ! Solver phase(s) to be executed
        INTEGER(mpl) :: error                  ! Error code
        INTEGER(mpl) :: msglvl                 ! Message level

        INTEGER(mpi) :: i
        INTEGER(mpl) :: ij
        INTEGER(mpl) :: idum(1)
        INTEGER(mpi) :: lun
        INTEGER(mpl) :: length
        INTEGER(mpi) :: nfill
        INTEGER(mpi) :: npdblk
        REAL(mpd) :: adum(1)
        REAL(mpd) :: ddum(1)

        INTEGER(mpl) :: iparm(64)
        REAL(mpd), ALLOCATABLE :: b( : )       ! Right hand side (of equations system)
        REAL(mpd), ALLOCATABLE :: x( : )       ! Solution (of equations system)
        SAVE

        lun=lunlog                       ! log file

        error  = 0 ! initialize error flag
        msglvl = ipddbg ! print statistical information
        npdblk=(nFitPar-1)/matbsz+1 ! number of row blocks

        IF(icalcm == 1) THEN
            mtype = 2                    ! positive definite symmetric matrix
            IF (nFitPar > nVarGlobalPar)  mtype = -2 ! indefinte symmetric matrix (Lagrange multipliers)

            !$POMP INST BEGIN(mspd00)
#ifdef SCOREP_USER_ENABLE
            SCOREP_USER_REGION_BY_NAME_BEGIN("UR_mspd00", SCOREP_USER_REGION_TYPE_COMMON)
#endif
            WRITE(*,*)
            WRITE(*,*) 'MSPARDISO: number of non-zero elements = ', csr3RowOffsets(npdblk+1)-csr3RowOffsets(1)
            ! fill up last block?
            nfill = npdblk*matbsz-nFitPar
            IF (nfill > 0) THEN
                WRITE(*,*) 'MSPARDISO: number of rows to fill up   = ', nfill
                ! end of last block
                ij = (csr3RowOffsets(npdblk+1)-csr3RowOffsets(1))*INT(matbsz,mpl)*INT(matbsz,mpl)
                DO i=1,nfill
                    globalMatD(ij) = 1.0_mpd
                    ij = ij-matbsz-1 ! back one row and one column in last block
                END DO
            END IF

            ! close previous PARADISO run
            IF (ipdmem > 0) THEN
                !.. Termination and release of memory
                phase = -1 ! release internal memory
                CALL pardiso_64(pt, maxfct, mnum, mtype, phase, INT(npdblk,mpl), adum, idum, idum, &
                    idum, nrhs, iparm, msglvl, ddum, ddum, error)
                IF (error /= 0) THEN
                    WRITE(lun,*) 'The following ERROR was detected: ', error
                    WRITE(*,'(A,2I10)') ' PARDISO release failed (phase, error): ', phase, error
                    IF (ipddbg == 0) WRITE(*,*) '        rerun with "debugPARDISO" for more info'
                    CALL peend(40,'Aborted, other error: PARDISO release')
                    STOP 'MSPARDISO: stopping due to error in PARDISO release'
                END IF
                ipdmem=0
            END IF

            !..
            !.. Set up PARDISO control parameter
            !..
            iparm=0 ! using defaults
            iparm(2)  =  2 ! fill-in reordering from METIS
            iparm(10) =  8 ! perturb the pivot elements with 1E-8
            iparm(18) = -1 ! Output: number of nonzeros in the factor LU
            iparm(19) = -1 ! Output: Mflops for LU factorization
            iparm(21) =  1 ! pivoting for symmetric indefinite matrices
            DO i=1, lenPARDISO
                iparm(listPARDISO(i)%label)=listPARDISO(i)%ivalue
            END DO
            IF (iparm(1) == 0) WRITE(lun,*) 'PARDISO using defaults '
            IF (iparm(43) /= 0) THEN
                WRITE(lun,*) 'PARDISO: computation of the diagonal of inverse matrix not implemented !'
                iparm(43) =  0 ! no computation of the diagonal of inverse matrix
            END IF

            ! necessary for the FIRST call of the PARDISO solver.
            DO i = 1, 64
                pt(i)%DUMMY =  0
            END DO
#ifdef SCOREP_USER_ENABLE
            SCOREP_USER_REGION_BY_NAME_END("UR_mspd00")
#endif
            !$POMP INST END(mspd00)
        END IF

        IF(icalcm == 1) THEN
            ! monitor progress
            IF(monpg1 > 0) THEN
                WRITE(lunlog,*) 'Decomposition of global matrix (A->L*D*L^t)'
                CALL monini(lunlog,monpg1,monpg2)
            END IF
            ! decompose and solve
            !.. Reordering and Symbolic Factorization, This step also allocates
            ! all memory that is necessary for the factorization
            !$POMP INST BEGIN(mspd11)
#ifdef SCOREP_USER_ENABLE
            SCOREP_USER_REGION_BY_NAME_BEGIN("UR_mspd11", SCOREP_USER_REGION_TYPE_COMMON)
#endif
            phase = 11 ! only reordering and symbolic factorization
            IF (matbsz > 1) THEN
                iparm(1) = 1 ! non default setting
                iparm(37) = matbsz ! using BSR3 instead of CSR3
            END IF
            IF (ipddbg > 0) THEN
                DO i=1,64
                    WRITE(lun,*) ' iparm(',i,') =', iparm(i)
                END DO
            END IF
            CALL pardiso_64(pt, maxfct, mnum, mtype, phase, INT(npdblk,mpl), globalMatD, csr3RowOffsets, csr3ColumnList, &
                idum, nrhs, iparm, msglvl, ddum, ddum, error)
#ifdef SCOREP_USER_ENABLE
            SCOREP_USER_REGION_BY_NAME_END("UR_mspd11")
#endif
            !$POMP INST END(mspd11)
            WRITE(lun,*) 'PARDISO reordering completed ... '
            WRITE(lun,*) 'PARDISO peak memory required (KB)', iparm(15)
            IF (ipddbg > 0) THEN
                DO i=1,64
                    WRITE(lun,*) ' iparm(',i,') =', iparm(i)
                END DO
            END IF
            IF (error /= 0) THEN
                WRITE(lun,*) 'The following ERROR was detected: ', error
                WRITE(*,'(A,2I10)') ' PARDISO decomposition failed (phase, error): ', phase, error
                IF (ipddbg == 0) WRITE(*,*) '        rerun with "debugPARDISO" for more info'
                CALL peend(40,'Aborted, other error: PARDISO reordering')
                STOP 'MSPARDISO: stopping due to error in PARDISO reordering'
            END IF
            IF (iparm(60) == 0) THEN
                ipdmem=ipdmem+max(iparm(15),iparm(16))+iparm(17) ! in core
            ELSE
                ipdmem=ipdmem+max(iparm(15),iparm(16))+iparm(63) ! out of core
            END IF
            WRITE(lun,*) 'Size (KB) of allocated memory  = ',ipdmem
            WRITE(lun,*) 'Number of nonzeros in factors  = ',iparm(18)
            WRITE(lun,*) 'Number of factorization MFLOPS = ',iparm(19)

            !.. Factorization.
            !$POMP INST BEGIN(mspd22)
#ifdef SCOREP_USER_ENABLE
            SCOREP_USER_REGION_BY_NAME_BEGIN("UR_mspd22", SCOREP_USER_REGION_TYPE_COMMON)
#endif
            phase = 22 ! only factorization
            CALL pardiso_64(pt, maxfct, mnum, mtype, phase, INT(npdblk,mpl), globalMatD, csr3RowOffsets, csr3ColumnList, &
                idum, nrhs, iparm, msglvl, ddum, ddum, error)
#ifdef SCOREP_USER_ENABLE
            SCOREP_USER_REGION_BY_NAME_END("UR_mspd22")
#endif
            !$POMP INST END(mspd22)
            WRITE(lun,*) 'PARDISO factorization completed ... '
            IF (ipddbg > 0) THEN
                DO i=1,64
                    WRITE(lun,*) ' iparm(',i,') =', iparm(i)
                END DO
            END IF
            IF (error /= 0) THEN
                WRITE(lun,*) 'The following ERROR was detected: ', error
                WRITE(*,'(A,2I10)') ' PARDISO decomposition failed (phase, error): ', phase, error
                IF (ipddbg == 0) WRITE(*,*) '        rerun with "debugPARDISO" for more info'
                CALL peend(40,'Aborted, other error: PARDISO factorization')
                STOP 'MSPARDISO: stopping due to error in PARDISO factorization'
            ENDIF
            IF (mtype < 0) THEN
                IF (iparm(14) > 0) &
                    WRITE(lun,*) 'Number of perturbed pivots     = ',iparm(14)
                WRITE(lun,*) 'Number of positive eigenvalues = ',iparm(22)-nfill
                WRITE(lun,*) 'Number of negative eigenvalues = ',iparm(23)
            ELSE IF (iparm(30) > 0) THEN
                WRITE(lun,*) 'Equation with bad pivot (<=0.) = ',iparm(30)
            END IF

            IF (monpg1 > 0) CALL monend()
        END IF

        ! backward/forward substitution
        !.. Back substitution and iterative refinement
        length=nFitPar+nfill
        CALL mpalloc(b,length,' PARDISO r.h.s')
        CALL mpalloc(x,length,' PARDISO solution')
        b(:nFitPar) = globalCorrections
        !$POMP INST BEGIN(mspd33)
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_REGION_BY_NAME_BEGIN("UR_mspd33", SCOREP_USER_REGION_TYPE_COMMON)
#endif
        iparm(6) = 0 ! don''t update r.h.s. with solution
        phase = 33 ! only solving
        CALL pardiso_64(pt, maxfct, mnum, mtype, phase, INT(npdblk,mpl), globalMatD, csr3RowOffsets, csr3ColumnList, &
            idum, nrhs, iparm, msglvl, b, x, error)
#ifdef SCOREP_USER_ENABLE
        SCOREP_USER_REGION_BY_NAME_END("UR_mspd33")
#endif
        !$POMP INST END(mspd33)
        globalCorrections = x(:nFitPar)
        CALL mpdealloc(x)
        CALL mpdealloc(b)
        WRITE(lun,*) 'PARDISO solve completed ... '
        IF (error /= 0) THEN
            WRITE(lun,*) 'The following ERROR was detected: ', error
            WRITE(*,'(A,2I10)') ' PARDISO decomposition failed (phase, error): ', phase, error
            IF (ipddbg == 0) WRITE(*,*) '        rerun with "debugPARDISO" for more info'
            CALL peend(40,'Aborted, other error: PARDISO solve')
            STOP 'MSPARDISO: stopping due to error in PARDISO solve'
        ENDIF

    END SUBROUTINE mspardiso
#endif
#endif
end module
