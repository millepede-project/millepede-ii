! Code converted using TO_F90 by Alan Miller
! Date: 2024-04-25  Time: 15:10:45

!> \file
!! Millepede-I to Millepede-II interface.
!!
!! \author Claus Kleinwort, DESY, 2018-2024 (Claus.Kleinwort@desy.de)
!!
!! \copyright
!! Copyright (c) 2024 Deutsches Elektronen-Synchroton,
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
!! Produce MP2 steering, constraints and binary file from MP1 calls.
!!
!! \verbatim
!! Implemented are:
!!     INITGL      initialization
!!     PARGLO       optional: initialize parameters with nonzero values
!!     PARSIG       optional: define sigma for single parameter
!!     INITCS       optional: constraints
!!     EQULOC      equations for local fit
!!     FITLOC      local parameter fit
!!     FITGLO      final global parameter fit
!!
!! Not implemented/dummy are:
!!     INITUN       optional: unit for iterations
!!     ERRPAR       optional: get parameter errors
!!     CORPAR       optional: get parameter correlations
!!     PRTGLO       optional: print results
!! \endverbatim

!> Initialization of package.
!!
!!     NAGB = number of global parameters
!!     DERGB(1) ... DERGB(NAGB) = derivatives w.r.t. global parameters
!!     NALC = number of local parameters (maximum)
!!     DERLC(1) ... DERLC(NALC) = derivatives w.r.t. local parameters
SUBROUTINE initgl(nagbar,nalcar,nstd,iprlim)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters

    INTEGER, INTENT(IN)                      :: nagbar
    INTEGER, INTENT(IN)                      :: nalcar
    INTEGER, INTENT(IN OUT)                  :: nstd
    INTEGER, INTENT(IN)                      :: iprlim
    INTEGER, PARAMETER :: mglobl=1400
    INTEGER, PARAMETER :: mlocal=10
    INTEGER, PARAMETER :: nstore=10000
    INTEGER, PARAMETER :: mcs=10
    INTEGER, PARAMETER :: mgl=mglobl+mcs
    !     derived parameters
    INTEGER, PARAMETER :: msymgb  =(mglobl*mglobl+mglobl)/2
    INTEGER, PARAMETER :: msym    =(mgl*mgl+mgl)/2
    INTEGER, PARAMETER :: msymlc=(mlocal*mlocal+mlocal)/2
    INTEGER, PARAMETER :: mrecta= mglobl*mlocal
    INTEGER, PARAMETER :: mglocs= mglobl*mcs
    INTEGER, PARAMETER :: msymcs= (mcs*mcs+mcs)/2
    DOUBLE PRECISION :: cgmat,clmat,clcmat,bgvec,blvec,  &
        corrm,corrv,summ,diag,scdiag,pparm,dparm,adercs, arhs
    LOGICAL :: scflag
    COMMON/lsqred/cgmat(msym),clmat(msymlc),clcmat(mrecta),  &
        diag(mgl),bgvec(mgl),blvec(mlocal),  &
        corrm(msymgb),corrv(mglobl),psigm(mglobl),  &
        pparm(mglobl),adercs(mglocs),arhs(mcs), dparm(mglobl),  &
        scdiag(mglobl),summ,scflag(mglobl),  &
        indgb(mglobl),indlc(mlocal),loctot,locrej,  &
        nagb,nalc,nsum,nhist,mhist(51),khist(51),lhist(51),  &
        nst,nfl,indst(nstore),arest(nstore),itert,lunit,ncs,  &
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim
    !     ------------------------------------------------------------------
    INTEGER :: ndr(7)
    DATA ndr/1,2,5,10,20,50,100/
    !     ...
    icnlim=iprlim
    IF(icnlim >= 0) WRITE(*,199)
199 FORMAT( ' MP-I to MP-II interface: create steering, binary files'/  &
        '                                                       '/  &
        '              *   o   o                        o       '/  &
        '                  o   o                        o       '/  &
        '   o ooooo    o   o   o    oo   ooo    oo    ooo   oo  '/  &
        '    o  o  o   o   o   o   o  o  o  o  o  o  o  o  o  o '/  &
        '    o  o  o   o   o   o   oooo  o  o  oooo  o  o  oooo '/  &
        '    o  o  o   o   o   o   o     ooo   o     o  o  o    '/  &
        '    o  o  o   o   oo  oo   oo   o      oo    ooo   oo  starting'/  &
        '                                o                      ')
    lunit=51 ! unit for binary file
    ncs  =0  ! number of constraints
    nagb=nagbar
    nalc=nalcar
    IF(icnlim >= 0) THEN
        WRITE(*,*) '                               '
        WRITE(*,*) 'Number of global parameters       ',nagb
        WRITE(*,*) 'Number of local parameters        ',nalc
    END IF
    IF(nstd /= 3) THEN
        WRITE(*,*) '                               '
        WRITE(*,*) 'MP-II uses nstd=3 istead of       ',nstd
    END IF
    WRITE(*,*) '                               '

    IF(nagb > mglobl.OR.nalc > mlocal) THEN
        WRITE(*,*) 'Too many parameter - STOP'
        STOP
    END IF
    !     reset input for global variables
    DO i=1,nagb
        pparm(i)=0.0   ! previous values of parameters set to zero
        psigm(i)=-1.0  ! no sigma defined for parameter I
    END DO
    !     open binary file
    OPEN(UNIT=lunit,ACCESS='SEQUENTIAL',FORM='UNFORMATTED', FILE='mp2tst.bin')
    !     output buffer
    nst=1
    indst(nst)=0
    arest(nst)=0.
END SUBROUTINE initgl

!> Initialize global parameters.
!!
!!     optional: initialize parameters with nonzero values
SUBROUTINE parglo(par)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters

    REAL, INTENT(IN)                         :: par(*)
    INTEGER, PARAMETER :: mglobl=1400
    INTEGER, PARAMETER :: mlocal=10
    INTEGER, PARAMETER :: nstore=10000
    INTEGER, PARAMETER :: mcs=10
    INTEGER, PARAMETER :: mgl=mglobl+mcs
    !     derived parameters
    INTEGER, PARAMETER :: msymgb  =(mglobl*mglobl+mglobl)/2
    INTEGER, PARAMETER :: msym    =(mgl*mgl+mgl)/2
    INTEGER, PARAMETER :: msymlc=(mlocal*mlocal+mlocal)/2
    INTEGER, PARAMETER :: mrecta= mglobl*mlocal
    INTEGER, PARAMETER :: mglocs= mglobl*mcs
    INTEGER, PARAMETER :: msymcs= (mcs*mcs+mcs)/2
    DOUBLE PRECISION :: cgmat,clmat,clcmat,bgvec,blvec,  &
        corrm,corrv,summ,diag,scdiag,pparm,dparm,adercs, arhs
    LOGICAL :: scflag
    COMMON/lsqred/cgmat(msym),clmat(msymlc),clcmat(mrecta),  &
        diag(mgl),bgvec(mgl),blvec(mlocal),  &
        corrm(msymgb),corrv(mglobl),psigm(mglobl),  &
        pparm(mglobl),adercs(mglocs),arhs(mcs), dparm(mglobl),  &
        scdiag(mglobl),summ,scflag(mglobl),  &
        indgb(mglobl),indlc(mlocal),loctot,locrej,  &
        nagb,nalc,nsum,nhist,mhist(51),khist(51),lhist(51),  &
        nst,nfl,indst(nstore),arest(nstore),itert,lunit,ncs,  &
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim
    !     ------------------------------------------------------------------


    DO i=1,nagb
        pparm(i)=par(i)
    END DO
END SUBROUTINE parglo

!> Define sigma for single parameter (optional).
SUBROUTINE parsig(INDEX,sigma)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters

    INTEGER, INTENT(IN)                      :: INDEX
    REAL, INTENT(IN)                         :: sigma
    INTEGER, PARAMETER :: mglobl=1400
    INTEGER, PARAMETER :: mlocal=10
    INTEGER, PARAMETER :: nstore=10000
    INTEGER, PARAMETER :: mcs=10
    INTEGER, PARAMETER :: mgl=mglobl+mcs
    !     derived parameters
    INTEGER, PARAMETER :: msymgb  =(mglobl*mglobl+mglobl)/2
    INTEGER, PARAMETER :: msym    =(mgl*mgl+mgl)/2
    INTEGER, PARAMETER :: msymlc=(mlocal*mlocal+mlocal)/2
    INTEGER, PARAMETER :: mrecta= mglobl*mlocal
    INTEGER, PARAMETER :: mglocs= mglobl*mcs
    INTEGER, PARAMETER :: msymcs= (mcs*mcs+mcs)/2
    DOUBLE PRECISION :: cgmat,clmat,clcmat,bgvec,blvec,  &
        corrm,corrv,summ,diag,scdiag,pparm,dparm,adercs, arhs
    LOGICAL :: scflag
    COMMON/lsqred/cgmat(msym),clmat(msymlc),clcmat(mrecta),  &
        diag(mgl),bgvec(mgl),blvec(mlocal),  &
        corrm(msymgb),corrv(mglobl),psigm(mglobl),  &
        pparm(mglobl),adercs(mglocs),arhs(mcs), dparm(mglobl),  &
        scdiag(mglobl),summ,scflag(mglobl),  &
        indgb(mglobl),indlc(mlocal),loctot,locrej,  &
        nagb,nalc,nsum,nhist,mhist(51),khist(51),lhist(51),  &
        nst,nfl,indst(nstore),arest(nstore),itert,lunit,ncs,  &
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim
    !     ------------------------------------------------------------------
    IF(INDEX < 1.OR.INDEX > nagb) RETURN
    IF(sigma < 0.0) RETURN
    psigm(INDEX)=sigma
END SUBROUTINE parsig

!> Set nonlinear flag for single parameter (optional).
SUBROUTINE nonlin(INDEX)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters

    INTEGER, INTENT(IN)                      :: INDEX
    INTEGER, PARAMETER :: mglobl=1400
    INTEGER, PARAMETER :: mlocal=10
    INTEGER, PARAMETER :: nstore=10000
    INTEGER, PARAMETER :: mcs=10
    INTEGER, PARAMETER :: mgl=mglobl+mcs
    !     derived parameters
    INTEGER, PARAMETER :: msymgb  =(mglobl*mglobl+mglobl)/2
    INTEGER, PARAMETER :: msym    =(mgl*mgl+mgl)/2
    INTEGER, PARAMETER :: msymlc=(mlocal*mlocal+mlocal)/2
    INTEGER, PARAMETER :: mrecta= mglobl*mlocal
    INTEGER, PARAMETER :: mglocs= mglobl*mcs
    INTEGER, PARAMETER :: msymcs= (mcs*mcs+mcs)/2
    DOUBLE PRECISION :: cgmat,clmat,clcmat,bgvec,blvec,  &
        corrm,corrv,summ,diag,scdiag,pparm,dparm,adercs, arhs
    LOGICAL :: scflag
    COMMON/lsqred/cgmat(msym),clmat(msymlc),clcmat(mrecta),  &
        diag(mgl),bgvec(mgl),blvec(mlocal),  &
        corrm(msymgb),corrv(mglobl),psigm(mglobl),  &
        pparm(mglobl),adercs(mglocs),arhs(mcs), dparm(mglobl),  &
        scdiag(mglobl),summ,scflag(mglobl),  &
        indgb(mglobl),indlc(mlocal),loctot,locrej,  &
        nagb,nalc,nsum,nhist,mhist(51),khist(51),lhist(51),  &
        nst,nfl,indst(nstore),arest(nstore),itert,lunit,ncs,  &
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim
    !     ------------------------------------------------------------------
    IF(INDEX < 1.OR.INDEX > nagb) RETURN
    nlnpa(INDEX)=1
END SUBROUTINE nonlin

!> Define unit for iterations (optional) - dummy.
SUBROUTINE initun(lun,cutfac)
    WRITE(*,*) ' INITUN is dummy !', lun, cutfac
END SUBROUTINE initun

!> Add constraint (optional).
SUBROUTINE constf(dercs,rhs)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters

    REAL, INTENT(IN)                         :: dercs(*)
    REAL, INTENT(IN)                         :: rhs
    INTEGER, PARAMETER :: mglobl=1400
    INTEGER, PARAMETER :: mlocal=10
    INTEGER, PARAMETER :: nstore=10000
    INTEGER, PARAMETER :: mcs=10
    INTEGER, PARAMETER :: mgl=mglobl+mcs
    !     derived parameters
    INTEGER, PARAMETER :: msymgb  =(mglobl*mglobl+mglobl)/2
    INTEGER, PARAMETER :: msym    =(mgl*mgl+mgl)/2
    INTEGER, PARAMETER :: msymlc=(mlocal*mlocal+mlocal)/2
    INTEGER, PARAMETER :: mrecta= mglobl*mlocal
    INTEGER, PARAMETER :: mglocs= mglobl*mcs
    INTEGER, PARAMETER :: msymcs= (mcs*mcs+mcs)/2
    DOUBLE PRECISION :: cgmat,clmat,clcmat,bgvec,blvec,  &
        corrm,corrv,summ,diag,scdiag,pparm,dparm,adercs, arhs
    LOGICAL :: scflag
    COMMON/lsqred/cgmat(msym),clmat(msymlc),clcmat(mrecta),  &
        diag(mgl),bgvec(mgl),blvec(mlocal),  &
        corrm(msymgb),corrv(mglobl),psigm(mglobl),  &
        pparm(mglobl),adercs(mglocs),arhs(mcs), dparm(mglobl),  &
        scdiag(mglobl),summ,scflag(mglobl),  &
        indgb(mglobl),indlc(mlocal),loctot,locrej,  &
        nagb,nalc,nsum,nhist,mhist(51),khist(51),lhist(51),  &
        nst,nfl,indst(nstore),arest(nstore),itert,lunit,ncs,  &
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim
    !     ------------------------------------------------------------------


    IF(ncs >= mcs) STOP '<INITCS> too many constraints'
    DO i=1,nagb
        adercs(nagb*ncs+i)=dercs(i)
    END DO
    ncs=ncs+1
    arhs(ncs)=rhs
END SUBROUTINE constf

!> Add single equation with its derivatives.
!!
!!     DERGB(1) ... DERGB(NAGB) = derivatives w.r.t. global parameters
!!     DERLC(1) ... DERLC(NALC) = derivatives w.r.t. local parameters
!!     RMEAS       = measured value
!!     SIGMA       = standard deviation
!!     (WGHT       = weight = 1/SIGMA**2)
SUBROUTINE equloc(dergb,derlc,rrmeas,sigma)

    REAL, INTENT(OUT)                        :: dergb(*)
    REAL, INTENT(OUT)                        :: derlc(*)
    REAL, INTENT(IN)                         :: rrmeas
    REAL, INTENT(IN)                         :: sigma
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    INTEGER, PARAMETER :: mglobl=1400
    INTEGER, PARAMETER :: mlocal=10
    INTEGER, PARAMETER :: nstore=10000
    INTEGER, PARAMETER :: mcs=10
    INTEGER, PARAMETER :: mgl=mglobl+mcs
    !     derived parameters
    INTEGER, PARAMETER :: msymgb  =(mglobl*mglobl+mglobl)/2
    INTEGER, PARAMETER :: msym    =(mgl*mgl+mgl)/2
    INTEGER, PARAMETER :: msymlc=(mlocal*mlocal+mlocal)/2
    INTEGER, PARAMETER :: mrecta= mglobl*mlocal
    INTEGER, PARAMETER :: mglocs= mglobl*mcs
    INTEGER, PARAMETER :: msymcs= (mcs*mcs+mcs)/2
    DOUBLE PRECISION :: cgmat,clmat,clcmat,bgvec,blvec,  &
        corrm,corrv,summ,diag,scdiag,pparm,dparm,adercs, arhs
    LOGICAL :: scflag
    COMMON/lsqred/cgmat(msym),clmat(msymlc),clcmat(mrecta),  &
        diag(mgl),bgvec(mgl),blvec(mlocal),  &
        corrm(msymgb),corrv(mglobl),psigm(mglobl),  &
        pparm(mglobl),adercs(mglocs),arhs(mcs), dparm(mglobl),  &
        scdiag(mglobl),summ,scflag(mglobl),  &
        indgb(mglobl),indlc(mlocal),loctot,locrej,  &
        nagb,nalc,nsum,nhist,mhist(51),khist(51),lhist(51),  &
        nst,nfl,indst(nstore),arest(nstore),itert,lunit,ncs,  &
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim
    !     ------------------------------------------------------------------

    !     ...
    rmeas=rrmeas
    IF(sigma <= 0.0) THEN
        DO i=1,nalc           ! local parameters
            derlc(i)=0.0         !  reset
        END DO
        DO i=1,nagb           ! global parameters
            dergb(i)=0.0         !  reset
        END DO
        RETURN
    END IF
    nonzer=0
    ialc=0
    iblc=-1
    DO i=1,nalc              ! count number of local parameters
        IF(derlc(i) /= 0.0) THEN
            nonzer=nonzer+1
            IF(ialc == 0) ialc=i ! first and last index
            iblc=i
        END IF
    END DO
    iagb=0
    ibgb=-1
    DO i=1,nagb              ! ... plus global parameters
        IF(dergb(i) /= 0.0) THEN
            nonzer=nonzer+1
            IF(iagb == 0) iagb=i ! first and last index
            ibgb=i
        END IF
    END DO
    IF(nst+nonzer+2 >= nstore) THEN
        nfl=1   ! set overflow flag
        RETURN  ! ignore data
    END IF
    nst=nst+1
    indst(nst)=0
    arest(nst)=rmeas
    DO i=ialc,iblc           ! local parameters
        IF(derlc(i) /= 0.0) THEN
            nst=nst+1
            indst(nst)=i         ! store index ...
            arest(nst)=derlc(i)  ! ... and value of nonzero derivative
            derlc(i)=0.0         ! reset
        END IF
    END DO
    nst=nst+1
    indst(nst)=0
    arest(nst)=sigma
    DO i=iagb,ibgb           ! global parameters
        IF(dergb(i) /= 0.0) THEN
            nst=nst+1
            indst(nst)=i         ! store index ...
            arest(nst)=dergb(i)  ! ... and value of nonzero derivative
            dergb(i)=0.0         ! reset
        END IF
    END DO
END SUBROUTINE equloc

!> Reset derivatives.
!!
!!     DERGB(1) ... DERGB(NAGB) = derivatives w.r.t. global parameters
!!     DERLC(1) ... DERLC(NALC) = derivatives w.r.t. local parameters
SUBROUTINE zerloc(dergb,derlc)

    REAL, INTENT(OUT)                        :: dergb(*)
    REAL, INTENT(OUT)                        :: derlc(*)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    INTEGER, PARAMETER :: mglobl=1400
    INTEGER, PARAMETER :: mlocal=10
    INTEGER, PARAMETER :: nstore=10000
    INTEGER, PARAMETER :: mcs=10
    INTEGER, PARAMETER :: mgl=mglobl+mcs
    !     derived parameters
    INTEGER, PARAMETER :: msymgb  =(mglobl*mglobl+mglobl)/2
    INTEGER, PARAMETER :: msym    =(mgl*mgl+mgl)/2
    INTEGER, PARAMETER :: msymlc=(mlocal*mlocal+mlocal)/2
    INTEGER, PARAMETER :: mrecta= mglobl*mlocal
    INTEGER, PARAMETER :: mglocs= mglobl*mcs
    INTEGER, PARAMETER :: msymcs= (mcs*mcs+mcs)/2
    DOUBLE PRECISION :: cgmat,clmat,clcmat,bgvec,blvec,  &
        corrm,corrv,summ,diag,scdiag,pparm,dparm,adercs, arhs
    LOGICAL :: scflag
    COMMON/lsqred/cgmat(msym),clmat(msymlc),clcmat(mrecta),  &
        diag(mgl),bgvec(mgl),blvec(mlocal),  &
        corrm(msymgb),corrv(mglobl),psigm(mglobl),  &
        pparm(mglobl),adercs(mglocs),arhs(mcs), dparm(mglobl),  &
        scdiag(mglobl),summ,scflag(mglobl),  &
        indgb(mglobl),indlc(mlocal),loctot,locrej,  &
        nagb,nalc,nsum,nhist,mhist(51),khist(51),lhist(51),  &
        nst,nfl,indst(nstore),arest(nstore),itert,lunit,ncs,  &
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim
    !     ------------------------------------------------------------------


    DO i=1,nalc           ! local parameters
        derlc(i)=0.0         !  reset
    END DO
    DO i=1,nagb           ! global parameters
        dergb(i)=0.0         !  reset
    END DO
END SUBROUTINE zerloc

!> Fit after end of local block.
SUBROUTINE fitloc
    !     faster(?) version
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    INTEGER, PARAMETER :: mglobl=1400
    INTEGER, PARAMETER :: mlocal=10
    INTEGER, PARAMETER :: nstore=10000
    INTEGER, PARAMETER :: mcs=10
    INTEGER, PARAMETER :: mgl=mglobl+mcs
    !     derived parameters
    INTEGER, PARAMETER :: msymgb  =(mglobl*mglobl+mglobl)/2
    INTEGER, PARAMETER :: msym    =(mgl*mgl+mgl)/2
    INTEGER, PARAMETER :: msymlc=(mlocal*mlocal+mlocal)/2
    INTEGER, PARAMETER :: mrecta= mglobl*mlocal
    INTEGER, PARAMETER :: mglocs= mglobl*mcs
    INTEGER, PARAMETER :: msymcs= (mcs*mcs+mcs)/2
    DOUBLE PRECISION :: cgmat,clmat,clcmat,bgvec,blvec,  &
        corrm,corrv,summ,diag,scdiag,pparm,dparm,adercs, arhs
    LOGICAL :: scflag
    COMMON/lsqred/cgmat(msym),clmat(msymlc),clcmat(mrecta),  &
        diag(mgl),bgvec(mgl),blvec(mlocal),  &
        corrm(msymgb),corrv(mglobl),psigm(mglobl),  &
        pparm(mglobl),adercs(mglocs),arhs(mcs), dparm(mglobl),  &
        scdiag(mglobl),summ,scflag(mglobl),  &
        indgb(mglobl),indlc(mlocal),loctot,locrej,  &
        nagb,nalc,nsum,nhist,mhist(51),khist(51),lhist(51),  &
        nst,nfl,indst(nstore),arest(nstore),itert,lunit,ncs,  &
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim
    !     ------------------------------------------------------------------


    IF(nst > 1) THEN ! write to binary file
        WRITE(lunit) nst*2,(arest(i),i=1,nst),(indst(i),i=1,nst)
    END IF

    ENTRY killoc
    IF(itert <= 1) THEN
        !        histogram of used store space
        IF(nfl == 0) THEN
            ibin=INT(1.0+50.0*FLOAT(nst)/FLOAT(nstore))
            ibin=MIN(ibin,50)
            khist(ibin)=khist(ibin)+1
        ELSE
            khist(51)=khist(51)+1
        END IF
    END IF
    nst=1      ! reset counter
    nfl=0      ! reset overflow flag
END SUBROUTINE fitloc

!> Final global fit.
SUBROUTINE fitglo(par)

    REAL, INTENT(IN OUT)                     :: par(*)

    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    INTEGER, PARAMETER :: mglobl=1400
    INTEGER, PARAMETER :: mlocal=10
    INTEGER, PARAMETER :: nstore=10000
    INTEGER, PARAMETER :: mcs=10
    INTEGER, PARAMETER :: mgl=mglobl+mcs
    !     derived parameters
    INTEGER, PARAMETER :: msymgb  =(mglobl*mglobl+mglobl)/2
    INTEGER, PARAMETER :: msym    =(mgl*mgl+mgl)/2
    INTEGER, PARAMETER :: msymlc=(mlocal*mlocal+mlocal)/2
    INTEGER, PARAMETER :: mrecta= mglobl*mlocal
    INTEGER, PARAMETER :: mglocs= mglobl*mcs
    INTEGER, PARAMETER :: msymcs= (mcs*mcs+mcs)/2
    DOUBLE PRECISION :: cgmat,clmat,clcmat,bgvec,blvec,  &
        corrm,corrv,summ,diag,scdiag,pparm,dparm,adercs, arhs
    LOGICAL :: scflag
    COMMON/lsqred/cgmat(msym),clmat(msymlc),clcmat(mrecta),  &
        diag(mgl),bgvec(mgl),blvec(mlocal),  &
        corrm(msymgb),corrv(mglobl),psigm(mglobl),  &
        pparm(mglobl),adercs(mglocs),arhs(mcs), dparm(mglobl),  &
        scdiag(mglobl),summ,scflag(mglobl),  &
        indgb(mglobl),indlc(mlocal),loctot,locrej,  &
        nagb,nalc,nsum,nhist,mhist(51),khist(51),lhist(51),  &
        nst,nfl,indst(nstore),arest(nstore),itert,lunit,ncs,  &
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim
    !     ------------------------------------------------------------------
    !     close binary file
    CLOSE(lunit)
    !     count parameter settings
    npar=0
    DO i=1,nagb           ! global parameters
        IF (pparm(i) /= 0.0.OR.psigm(i) >= 0.0) npar=npar+1
        par(i)=0.
    END DO
    !     create steering file
    luns=7
    OPEN(UNIT=luns,ACCESS='SEQUENTIAL',FORM='FORMATTED', FILE='mp2str.txt')
    WRITE(luns,101) '*            Default test steering file'
    WRITE(luns,101) 'fortranfiles ! following bin files are fortran'
    WRITE(luns,101) 'mp2tst.bin   ! binary data file'
    WRITE(luns,101) ' '
    IF(ncs > 0) THEN
        WRITE(luns,101) 'mp2con.txt   ! constraints text file '
        WRITE(luns,101) ' '
    END IF
    IF(npar > 0) THEN
        WRITE(luns,101) 'Parameter ! with start values or pre-sigma'
        DO i=1,nagb           ! global parameters
            IF (pparm(i) /= 0.0.OR.psigm(i) >= 0.0) THEN
                presig=psigm(i)
                !     parameter fixed by PRESIG = 0 (MP1) -> <0 (MP2)
                IF (presig == 0.0) presig=-1.0
                WRITE(luns,103) i, pparm(i), presig
            END IF
        END DO
        WRITE(luns,101) ' '
    END IF
    WRITE(luns,101) 'chisqcut 30.0 6.0 ! Chi2/ndf cut factors'
    WRITE(luns,101) 'method inversion 3 0.001   ! Gauss matrix inversion'
    WRITE(luns,101) ' '
    WRITE(luns,101) 'end ! optional for end-of-data'
    CLOSE(luns)
    !     create constraints file
    IF (ncs > 0) THEN
        lunc=9
        OPEN(UNIT=lunc,ACCESS='SEQUENTIAL',FORM='FORMATTED', FILE='mp2con.txt')
        DO i=0, ncs-1
            WRITE(lunc,101) 'Constraint  0.0'
            DO j=1,nagb
                IF (adercs(nagb*i+j) /= 0.0) WRITE(lunc,102) j, adercs(nagb*i+j)
            END DO
        END DO
        CLOSE(lunc)
    END IF
    WRITE(*,199)
199 FORMAT( '                                                       '/  &
        '              *   o   o                        o       '/  &
        '                  o   o                        o       '/  &
        '   o ooooo    o   o   o    oo   ooo    oo    ooo   oo  '/  &
        '    o  o  o   o   o   o   o  o  o  o  o  o  o  o  o  o '/  &
        '    o  o  o   o   o   o   oooo  o  o  oooo  o  o  oooo '/  &
        '    o  o  o   o   o   o   o     ooo   o     o  o  o    '/  &
        '    o  o  o   o   oo  oo   oo   o      oo    ooo   oo   ending.'/  &
        '                                o                      '/  &
        '                                                       '/  &
        ' MP-I to MP-II interface: create steering, binary files')
101 FORMAT(a)
102 FORMAT(i8,f10.5)
103 FORMAT(i8,2F10.5)

END SUBROUTINE fitglo

!> Return error for parameter I - dummy.
FUNCTION errpar(i)
    errpar=0.0
    WRITE(*,*) ' ERRPAR is dummy !', i
END FUNCTION errpar

!> Return correlation between parameters I and J - dummy.
FUNCTION corpar(i,j)
    corpar=0.0
    WRITE(*,*) ' CORPAR is dummy !', i, j
END FUNCTION corpar

!> Print result on file - dummy.
SUBROUTINE prtglo(lun)
    WRITE(*,*) ' PRTGLO is dummy !', lun
END SUBROUTINE prtglo
