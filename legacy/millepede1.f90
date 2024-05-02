! Code converted using TO_F90 by Alan Miller
! Date: 2024-04-29  Time: 13:19:49

!> \file
!! Millepede(-I) subroutines.
!!
!! \author Volker Blobel, University Hamburg, 2000 (original Fortran77 code)
!! \author Claus Kleinwort, DESY, 2024 (reformatting)
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
!! \verbatim
!!
!!                Millepede - Linear Least Squares
!!                ================================
!!    A Least Squares Method for Detector Alignment - Fortran code
!!
!!     TESTMP      short test program for detector aligment
!!                 with GENER (generator) + ZRAND, ZNORM (random gen.)
!!
!!     The execution of the test program needs a MAIN program:
!!                 CALL TESTMP(0)
!!                 CALL TESTMP(1)
!!                 END
!!
!!     INITGL      initialization
!!     PARGLO       optional: initialize parameters with nonzero values
!!     PARSIG       optional: define sigma for single parameter
!!     INITUN       optional: unit for iterations
!!     CONSTF       optional: constraints
!!     EQULOC      equations for local fit
!!     ZERLOC
!!     FITLOC      local parameter fit (+entry KILLOC)
!!     FITGLO      final global parameter fit
!!     ERRPAR      parameter errors
!!     CORPAR      parameter correlations
!!     PRTGLO      print result on file
!!
!!     Special matrix subprograms (all in Double Precision):
!!        SPMINV   matrix inversion + solution
!!        SPAVAT   special double matrix product
!!        SPAX     product matrix times vector
!!
!!     PXHIST      histogram printing
!!     CHFMT       formatting of real numbers
!!     CHINDL      limit of chi**2/nd
!!     FITMUT/FITMIN vector input/ouput for parallel processing
!! \endverbatim

!> Test of millepede.
!!
!!     IARG = 0   test with constraint
!!     IARG = 1   test without constraint
SUBROUTINE testmp(iarg)       ! test program IARG = 0 or 1

    INTEGER, INTENT(IN)                      :: iarg
    REAL :: dergb(10),derlc(2),par(10)
    INTEGER, PARAMETER :: nplan=10
    REAL :: x(nplan),y(nplan),sigma(nplan),bias(nplan),heff(nplan)
    REAL :: frot(nplan)
    DATA x/5.0,9.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0/
    DATA sigma/2*0.0020,8*0.0300/
    DATA bias /0.00,-0.04,0.15,0.03,-0.075,0.045,0.035,-0.08,0.09,-0.05/
    DATA heff/4*0.9,0.50,5*0.9/
    !     ...
    ncases=1000
    CALL initgl(10,2,3,1)        ! define dimension parameters
    CALL parsig(1,0.0)           ! parameter 1 is fixed
    IF(iarg == 0) THEN
        frot(1)=0.0
        DO i=2,10
            frot(i)=1.0/(x(i)-x(1))
        END DO
        CALL constf(frot,0.0)     ! constraint: total rotation zero
    ELSE
        DO i=3,10
            CALL parsig(i,0.03)      ! parameters 3...8
        END DO
    END IF
    CALL initun(11,10.0)         ! option iterations
    CALL zerloc(dergb,derlc)     ! initialization of arrays to zero
    !     -------------- loop
    DO nc=1,ncases
        !      generate straight line parameters
        CALL gener(y,a,b,x,sigma,bias,heff)
        DO i=1,nplan
            !       calibrate  Y(I) = A + B * X(I)
            IF(y(i) /= 0.0) THEN
                derlc(1)= 1.0
                derlc(2)= x(i)
                dergb(i)= 1.0
                CALL equloc(dergb,derlc,y(i),sigma(i)) ! single measurement
            END IF
        END DO
        CALL fitloc
    END DO ! end loop cases
!     -------------- loop end
CONTINUE
CALL fitglo(par)             ! final solution
END SUBROUTINE testmp

!> Generate Y values.
SUBROUTINE gener(y,a,b,x,sigma,bias,heff)

    INTEGER, PARAMETER :: nplan=10

    REAL, INTENT(OUT)                        :: y(*)
    REAL, INTENT(OUT)                        :: a
    REAL, INTENT(OUT)                        :: b
    REAL, INTENT(IN)                         :: x(nplan)
    REAL, INTENT(IN)                         :: sigma(nplan)
    REAL, INTENT(IN)                         :: bias(nplan)
    REAL, INTENT(IN)                         :: heff(nplan)

    !     ...
    !     generate straight line parameters
    a=20.0*zrand()-10.0        ! A =   -10 ...  +10   uniform
    b= 0.3*znorm()             ! B =  -0.3 ...  +0.3  gaussian
    DO i=1,nplan  ! planes
        IF(zrand() < heff(i)) THEN
            y(i)=a+b*x(i) + bias(i)+sigma(i)*znorm()
        ELSE
            y(i)=0.0  ! unmeasured
        END IF
    END DO
END SUBROUTINE gener

!> Return random number U(0,1).
FUNCTION zrand()
    !     (simple generator, showing principle)
    PARAMETER (ia=205,ic=29573,im=139968)
    DATA last/4711/
    last=MOD(ia*last+ic,im)
    IF(last == 0) last=MOD(ia*last+ic,im)
    zrand=FLOAT(last)/FLOAT(im)
END

!> Return random number N(0,1).
FUNCTION znorm()
    !     (simple generator, showing principle)
    PARAMETER (ia=205,ic=29573,im=139968)
    DATA last/4711/
    znorm=-6.0
    DO i=1,12
        last=MOD(ia*last+ic,im)
        znorm=znorm+FLOAT(last)/FLOAT(im)
    END DO
END

!> Initialization of package.
!!
!!     NAGB = number of global parameters
!!     DERGB(1) ... DERGB(NAGB) = derivatives w.r.t. global parameters
!!     NALC = number of local parameters (maximum)
!!     DERLC(1) ... DERLC(NALC) = derivatives w.r.t. local parameters
SUBROUTINE initgl(nagbar,nalcar,nstd,iprlim)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    INTEGER :: ndr(7)
    DATA ndr/1,2,5,10,20,50,100/
    !     ...
    icnlim=iprlim
    IF(icnlim >= 0) WRITE(*,199)
199 FORMAT( '                                                       '/  &
        '              *   o   o                        o       '/  &
        '                  o   o                        o       '/  &
        '   o ooooo    o   o   o    oo   ooo    oo    ooo   oo  '/  &
        '    o  o  o   o   o   o   o  o  o  o  o  o  o  o  o  o '/  &
        '    o  o  o   o   o   o   oooo  o  o  oooo  o  o  oooo '/  &
        '    o  o  o   o   o   o   o     ooo   o     o  o  o    '/  &
        '    o  o  o   o   oo  oo   oo   o      oo    ooo   oo  starting'/  &
        '                                o                      ')
    itert=0  ! reset iteration counter/flag
    lunit=0  ! unit for scratch storage
    cfact=1.0! factor for cut during iterations
    ncs  =0  ! number of constraints
    icnpr=0  ! number of printouts
    loctot=0 ! total number of local fits
    locrej=0 ! number of rejected fits
    nstdev=MAX(0,MIN(nstd,3))  ! 0 ... 4
    cfactr=1.0
    nagb=nagbar
    nalc=nalcar
    IF(icnlim >= 0) THEN
        WRITE(*,*) '                               '
        WRITE(*,*) 'Number of global parameters       ',nagb
        WRITE(*,*) 'Number of local parameters        ',nalc
        WRITE(*,*) '   Number of standard deviations  ',nstdev
        WRITE(*,*) '   Number of test printouts       ',icnlim
    END IF
    IF(nagb > mglobl.OR.nalc > mlocal) THEN
        WRITE(*,*) 'Too many parameter - STOP'
        STOP
    END IF
    IF(nstdev /= 0.AND.icnlim >= 0) THEN
        WRITE(*,*) 'Final cut corresponds to',nstdev, ' standard deviations.'
        WRITE(*,*) 'The actual cuts are made in Chisquare/Ndf at:'
        DO i=1,7
            WRITE(*,101) ndr(i),chindl(nstdev,ndr(i))
101         FORMAT(20X,'Ndf =',i4,'   Limit =',f8.2)
        END DO
    END IF
    !     reset matrices for global variables
    DO i=1,nagb
        bgvec(i)=0.0
        pparm(i)=0.0   ! previous values of parameters set to zero
        dparm(i)=0.0   ! corrections of parameters zero
        psigm(i)=-1.0  ! no sigma defined for parameter I
        nlnpa(i)=0     ! linear parameter by default
    END DO
    DO i=1,(nagb*nagb+nagb)/2
        cgmat(i)=0.0
    END DO
    !     reset histogram
    nhist=0
    DO i=1,51
        mhist(i)=0
        khist(i)=0
        lhist(i)=0
    END DO
    !     reset matrices for local variables
    summ=0.0D0
    nsum=0
    DO i=1,nalc
        blvec(i)=0.0
    END DO
    DO i=1,(nalc*nalc+nalc)/2
        clmat(i)=0.0
    END DO
    nst=0        ! reset counter for derivatives
    nfl=0        ! reset flag for derivative storage
END SUBROUTINE initgl

!> Initialize global parameters.
!!
!!     optional: initialize parameters with nonzero values
SUBROUTINE parglo(par)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    REAL :: par(*)
    DO i=1,nagb
        pparm(i)=par(i)
    END DO
END SUBROUTINE parglo

!> Define sigma for single parameter (optional).
SUBROUTINE parsig(INDEX,sigma)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    IF(INDEX < 1.OR.INDEX > nagb) RETURN
    IF(sigma < 0.0) RETURN
    psigm(INDEX)=sigma
END SUBROUTINE parsig

!> Set nonlinear flag for single parameter (optional).
SUBROUTINE nonlin(INDEX)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    IF(INDEX < 1.OR.INDEX > nagb) RETURN
    nlnpa(INDEX)=1
END SUBROUTINE nonlin

!> Define unit for iterations (optional).
SUBROUTINE initun(lun,cutfac)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    LOGICAL :: opn
    !     test file for existence
    INQUIRE(UNIT=lun,OPENED=opn,IOSTAT=ios)
    IF(ios /= 0) THEN
        STOP '<INITUN: Inquire error>'
    END IF
    IOSTAT=0
    IF(.NOT.opn) THEN
        OPEN(UNIT=lun,FORM='UNFORMATTED',STATUS='SCRATCH',IOSTAT=ios)
        IF(ios /= 0) THEN
            STOP '<INITUN: Open error>'
        END IF
        IF(icnlim >= 0) WRITE(*,*) 'Scratch file opened'
    END IF
    lunit=lun
    cfactr=MAX(1.0,cutfac)     ! > 1.0
    IF(icnlim >= 0) WRITE(*,*) 'Initial cut factor is',cfactr
    itert=1  ! iteration 1 is first iteration
END SUBROUTINE initun

!> Add constraint (optional).
SUBROUTINE constf(dercs,rhs)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    REAL :: dercs(*)
    IF(ncs >= mcs) STOP '<INITCS> too many constraints'
    DO i=1,nagb
        adercs(nagb*ncs+i)=dercs(i)
    END DO
    ncs=ncs+1
    arhs(ncs)=rhs
    WRITE(*,*) 'Number of constraints increased to ',ncs
END SUBROUTINE constf

!> Add single equation with its derivatives.
!!
!!     DERGB(1) ... DERGB(NAGB) = derivatives w.r.t. global parameters
!!     DERLC(1) ... DERLC(NALC) = derivatives w.r.t. local parameters
!!     RMEAS       = measured value
!!     SIGMA       = standard deviation
!!     (WGHT       = weight = 1/SIGMA**2)
SUBROUTINE equloc(dergb,derlc,rrmeas,sigma)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    REAL :: dergb(*),derlc(*)
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
    wght=1.0/sigma**2
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
    arest(nst)=wght
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
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    REAL :: dergb(*),derlc(*)
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
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    icnpr=icnpr+1
    incrp=0
    IF(itert == 1) THEN
        IF(nst > 0) THEN ! write to scratch file
            WRITE(lunit) nst,(indst(i),i=1,nst),(arest(i),i=1,nst)
        END IF
    END IF
    !     reset matrices for local variables ooooooooooooooooooooooooooooooo
    summ=0.0D0
    nsum=0
    DO i=1,nalc
        blvec(i)=0.0
    END DO
    DO i=1,(nalc*nalc+nalc)/2
        clmat(i)=0.0
    END DO
    !     reset pointer matrix for mixed variables
    DO i=1,nagb
        indnz(i)=0
    END DO
    nagbn=0          ! actual number of global parameters
    !     normal equations - symmetric matrix for local parameters
    IF(nst <= 1) GO TO 100
    ist=0
10  ja =0 ! new equation
    jb =0
20  ist=ist+1
    IF(ist > nst.OR.indst(ist) == 0) THEN
        IF(ja == 0) THEN
            ja=ist             ! first zero: measured value
        ELSE IF(jb == 0) THEN
            jb=ist             ! second zero: weight
        ELSE
            ist=ist-1          ! end of equation
            rmeas=arest(ja)          ! use the data
            !           subtract global ... from measured value
            DO j=1,ist-jb
                ij=indst(jb+j)           ! subtract from measured value ...
                IF(nlnpa(ij) == 0) THEN  ! ... for linear parameter
                    rmeas=rmeas-arest(jb+j)*REAL(pparm(ij)+dparm(ij))
                ELSE                     ! ... for nonlinear parameter
                    rmeas=rmeas-arest(jb+j)*REAL(dparm(ij))
                END IF
            END DO
            !           end-subtract
            wght =arest(jb)          ! ... and the weight
            DO j=1,jb-ja-1           ! number of derivatives
                ij=indst(ja+j)
                blvec(ij)=blvec(ij)+wght*rmeas*arest(ja+j)
                DO k=1,j
                    ik=indst(ja+k)
                    jk=(ij*ij-ij)/2+ik
                    clmat(jk)=clmat(jk)+wght*arest(ja+j)*arest(ja+k)
                END DO
            END DO
    
            IF(ist == nst) GO TO 21
            GO TO 10 ! next equation
        END IF
    END IF
    IF(ist <= nst) GO TO 20
    !     end of data - determine local fit parameters
21  CALL spminv(clmat,blvec,nalc,nrank,scdiag,scflag)
    ! inversion and solution
    IF(icnpr <= icnlim) THEN
        WRITE(*,*) '                                                  '
        WRITE(*,*) '__________________________________________________'
        WRITE(*,*) 'Printout of local fit (FITLOC) with rank=',nrank
        WRITE(*,*) '   Result of local fit: (Index/Parameter/error)'
        WRITE(*,103) (j,blvec(j),SQRT(clmat((j*j+j)/2)),j=1,jb-ja-1)
103     FORMAT(2(i6,2G12.4))
    END IF
    !     calculate residuals
    summ=0.0D0
    nsum=0
    !     second loop for residual calculation -----------------------------
    ist=0
30  ja =0 ! new equation
    jb =0
40  ist=ist+1
    IF(ist > nst.OR.indst(ist) == 0) THEN
        IF(ja == 0) THEN
            ja=ist             ! first zero: measured value
        ELSE IF(jb == 0) THEN
            jb=ist             ! second zero: weight
        ELSE
            ist=ist-1          ! end of equation
            IF(icnpr <= icnlim.AND.incrp <= 11) THEN
                nderlc=jb-ja-1
                ndergl=ist-jb
                incrp=incrp+1
                WRITE(*,*) '            '
                WRITE(*,*) incrp,'. equation: ',  &
                    'measured value ',arest(ja),' +-',1.0/SQRT(arest(jb))
                IF(incrp <= 11) THEN
                    WRITE(*,*) ' number of derivates (global, local):',ndergl,',',nderlc
                    WRITE(*,*) '  Global derivatives are: (index/derivative/parvalue)'
                    WRITE(*,101) (indst(jb+j),arest(jb+j), pparm(indst(jb+j)),j=1,ist-jb)
101                 FORMAT(2(i6,2G12.4))
                    WRITE(*,*) '  Local derivatives are: (index/derivative)'
                    WRITE(*,102) (indst(ja+j),arest(ja+j),j=1,jb-ja-1)
102                 FORMAT(3(i6,g14.6))
                END IF
                IF(incrp == 11) THEN
                    WRITE(*,*)  '...   (+ further equations)'
                END IF
            END IF
            rmeas=arest(ja)          ! use the data
            !           subtract global ... from measured value
            DO j=1,ist-jb
                ij=indst(jb+j)           ! subtract from measured value ...
                IF(nlnpa(ij) == 0) THEN  ! ... for linear parameter
                    rmeas=rmeas-arest(jb+j)*REAL(pparm(ij)+dparm(ij))
                ELSE                     ! ... for nonlinear parameter
                    rmeas=rmeas-arest(jb+j)*REAL(dparm(ij))
                END IF
            END DO
            !           end-subtract
            wght =arest(jb)          ! ... and the weight
            DO j=1,jb-ja-1           ! number of derivatives
                ij=indst(ja+j)
                rmeas=rmeas-arest(ja+j)*REAL(blvec(ij))
            END DO
            IF(icnpr <= icnlim) THEN
                WRITE(*,104) wght*rmeas**2,rmeas
104             FORMAT(' Chi square contribution=',f8.2,  &
                    5X,g12.4,'= residuum')
            END IF
            !           residual histogram
            ihbin=INT(1.0+5.0*(rmeas*SQRT(wght)+5.0))
            IF(ihbin < 1) ihbin=51
            IF(ihbin > 50) ihbin=51
            lhist(ihbin)=lhist(ihbin)+1  ! residial histogram
            summ=summ+wght*rmeas**2
            nsum=nsum+1              ! ... and count equation
            IF(ist == nst) GO TO 41
            GO TO 30 ! next equation
        END IF
    END IF
    IF(ist <= nst) GO TO 40
41  ndf=nsum-nrank
    IF(icnpr <= icnlim) THEN
        WRITE(*,*) 'Entry number',icnpr
        WRITE(*,*) 'Final chi square, degrees of freedom',summ,ndf
    END IF
    rms=0.0
    IF(ndf > 0) THEN             ! histogram entry
        rms=REAL(summ)/REAL(ndf)
        nhist=nhist+1
        ihbin=INT(1+5.0*rms)                 ! bins width 0.2
        IF(ihbin < 1) ihbin= 1
        IF(ihbin > 50) ihbin=51
        mhist(ihbin)=mhist(ihbin)+1
    END IF
    loctot=loctot+1
    !     make eventual cut
    IF(nstdev /= 0.AND.ndf > 0) THEN
        cutval=chindl(nstdev,ndf)*cfactr
        IF(icnpr <= icnlim) THEN
            WRITE(*,*) 'Reject if Chisq/Ndf=',rms,' > ',cutval
        END IF
        IF(rms > cutval) THEN
            locrej=locrej+1
            GO TO 100
        END IF
    END IF
    !     third loop for global parameters ---------------------------------
    ist=0
50  ja =0 ! new equation
    jb =0
60  ist=ist+1
    IF(ist > nst.OR.indst(ist) == 0) THEN
        IF(ja == 0) THEN
            ja=ist             ! first zero: measured value
        ELSE IF(jb == 0) THEN
            jb=ist             ! second zero: weight
        ELSE
            ist=ist-1          ! end of equation
            rmeas=arest(ja)          ! use the data
            !           subtract global ... from measured value
            DO j=1,ist-jb
                ij=indst(jb+j)           ! subtract from measured value ...
                IF(nlnpa(ij) == 0) THEN  ! ... for linear parameter
                    rmeas=rmeas-arest(jb+j)*REAL(pparm(ij)+dparm(ij))
                ELSE                     ! ... for nonlinear parameter
                    rmeas=rmeas-arest(jb+j)*REAL(dparm(ij))
                END IF
            END DO
            !           end-subtract
            wght =arest(jb)          ! ... and the weight
            !           normal equations - symmetric matrix for global parameters
            DO j=1,ist-jb
                ij=indst(jb+j)
                bgvec(ij)=bgvec(ij)+wght*rmeas*arest(jb+j)
                DO k=1,j
                    ik=indst(jb+k)
                    jk=(ij*ij-ij)/2+ik
                    cgmat(jk)=cgmat(jk)+wght*arest(jb+j)*arest(jb+k)
                END DO
            END DO
            !           normal equations - rectangular matrix for global/local pars
            DO j=1,ist-jb
                ij=indst(jb+j)   ! index of global variable
                !            new code start
                ijn=indnz(ij)    ! get index of index
                IF(ijn == 0) THEN
                    !               new global variable - initialize matrix row
                    DO k=1,nalc
                        clcmat(nagbn*nalc+k)=0.0
                    END DO
                    nagbn=nagbn+1
                    indnz(ij)=nagbn     ! insert pointer
                    indbk(nagbn)=ij     ! ... and pointer back
                    ijn=nagbn
                END IF
                !            new code end
                DO k=1,jb-ja-1
                    ik=indst(ja+k)
                    !             JK=IK+(IJ-1)*NALC    ! old code
                    jk=ik+(ijn-1)*nalc   ! new code
                    clcmat(jk)=clcmat(jk)+wght*arest(jb+j)*arest(ja+k) !<=
                END DO
            END DO
            IF(ist == nst) GO TO 70
            GO TO 50 ! next equation
        END IF
    END IF
    IF(ist <= nst) GO TO 60
    !     -------------------------------------------------------------
    !     update global matrices
70  CALL spavat(clmat,clcmat,corrm,nalc,nagbn) ! correction matrix
    CALL spax(clcmat,blvec,corrv,nagbn,nalc)   ! correction vector
    ijn=0
    DO in=1,nagbn
        i=indbk(in)      ! get pointer back to global index
        bgvec(i)=bgvec(i)-corrv(in)
        DO jn=1,in
            j=indbk(jn)
            ijn=ijn+1
            IF(i >= j) THEN
                ij=j+(i*i-i)/2
            ELSE
                ij=i+(j*j-j)/2
            END IF
            cgmat(ij)=cgmat(ij)-corrm(ijn)
        END DO
    END DO
    ENTRY killoc
100 IF(itert <= 1) THEN
        !        histogram of used store space
        IF(nfl == 0) THEN
            ibin=INT(1.0+50.0*FLOAT(nst)/FLOAT(nstore))
            ibin=MIN(ibin,50)
            khist(ibin)=khist(ibin)+1
        ELSE
            khist(51)=khist(51)+1
        END IF
    END IF
    nst=0      ! reset counter
    nfl=0      ! reset overflow flag
END SUBROUTINE fitloc

!> Final global fit.
SUBROUTINE fitglo(par)
    REAL :: par(*)
    CHARACTER (LEN=2) :: patext
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    DOUBLE PRECISION :: dsum
    IF(itert <= 1) itelim=10 ! maximum number of iterations
    IF(itert /= 0) THEN
        lthist=0
        DO i=1,51
            lthist=lthist+lhist(i)
        END DO
        IF(icnlim >= 0) THEN
            WRITE(*,*) '  Initial residual histogram:',  &
                '  Total',lthist,' entries with',lhist(51), ' overflows'
            CALL pxhist(lhist,50,-5.0,+5.0) ! histogram printout
        END IF
    END IF
    khtot=0
    DO i=1,51
        khtot=khtot+khist(i)
    END DO
    IF(icnlim >= 0) THEN
        WRITE(*,*) ' '
        WRITE(*,*) 'Histogram of used local fit storage:',  &
            ' total',khtot,' local fits with',khist(51), ' overflows'
    END IF
    IF(khist(51) /= 0.AND.icnlim >= 0) THEN
        WRITE(*,*) 'Parameter NSTORE to small! Change and rerun!'
    END IF
    IF(icnlim >= 0) THEN
        CALL pxhist(khist,50,0.0,FLOAT(nstore)) ! histogram printout
    END IF
10  IF(icnlim >= 0) THEN
        WRITE(*,*) '    '
        WRITE(*,*) '... making global fit ...'
    END IF
    nn=0
    ii=0          ! modify matrix acccording to PSIGM value
    DO i=1,nagb
        ll=ii
        ii=ii+i
        IF(psigm(i) == 0.0) THEN
            nn=nn+1        ! count number of fixed parameters
            DO j=1,nagb
                ll=ll+1
                IF(j > i) ll=ll+j-2
                cgmat(ll)=0.0 ! reset row and column for parameter
            END DO
        ELSE IF(psigm(i) > 0.0) THEN
            cgmat(ii)=cgmat(ii)+1.0/psigm(i)**2 ! add to diagonal
        END IF
    END DO
    ii=0                ! start saving diagonal elements
    DO i=1,nagb
        ii=ii+i
        diag(i)=cgmat(ii)  ! save original diagonal elements
    END DO
    nvar=nagb
    WRITE(*,*) 'Number of constraints ',ncs
    IF(ncs /= 0) THEN   ! add constraints
        ii=(nagb*nagb+nagb)/2
        DO i=1,ncs  ! loop on constraints
            dsum=arhs(i)
            DO j=1,nagb
                cgmat(ii+j)=adercs(nagb*(i-1)+j)
                dsum=dsum-adercs(nagb*(i-1)+j)*(pparm(j)+dparm(j))
            END DO
            DO j=1,i
                cgmat(ii+nagb+j)=0.0
            END DO
            nvar=nvar+1
            ii  =ii+nvar
            bgvec(nvar)=dsum
        END DO
    END IF
    !     =================================================
    CALL spminv(cgmat,bgvec,nvar,nrank,scdiag,scflag)!matrix inversion
    ndefec=nvar-nrank-nn  ! rank defect
    !     =====================================================
    DO i=1,nagb
        dparm(i)=dparm(i)+bgvec(i)   ! accumulate corrections
    END DO
    CALL prtglo(66)
    IF(icnlim >= 0) THEN
        WRITE(*,*) 'The rank defect of the symmetric',nvar,'-by-',nvar,  &
            ' matrix is ',ndefec,' (should be zero).'
    END IF
    IF(itert == 0.OR.nstdev == 0.OR.itert >= itelim) GO TO 90
    !     iterations
    IF(icnlim >= 0) THEN
        WRITE(*,*) '                   '
        WRITE(*,*) '  Total',loctot,' local fits,',locrej,' rejected.'
        WRITE(*,*) '  Histogram of Chisq/Ndf:',  &
            '  Total',nhist,' entries with',mhist(51), ' overflows'
        CALL pxhist(mhist,50,0.0,10.0) ! histogram printout
    END IF
    !     reset histogram
    nhist=0
    DO i=1,51
        mhist(i)=0
        lhist(i)=0
    END DO
    itert=itert+1
    loctot=0
    locrej=0
    IF(cfactr /= 1.0) THEN
        cfactr=SQRT(cfactr)
        IF(cfactr < 1.2) THEN
            cfactr=1.0
            itelim=itert+1
        END IF
    END IF
    IF(icnlim >= 0) WRITE(*,107) itert,cfactr
107 FORMAT(' Iteration',i3,' with cut factor=',f6.2)
    !     reset matrices for global variables
    DO i=1,nagb
        bgvec(i)=0.0
    END DO
    DO i=1,(nagb*nagb+nagb)/2
        cgmat(i)=0.0
    END DO
    REWIND lunit
20  READ(lunit,END=10) nst,(indst(i),i=1,nst),(arest(i),i=1,nst)
    CALL fitloc
    GO TO 20
    !     ==================================================================
90  IF(icnlim >= 0) THEN
        WRITE(*,*) '                         '
        WRITE(*,*) '         Result of fit for global parameters'
        WRITE(*,*) '         ==================================='
        WRITE(*,101)
    END IF
    ii=0
    DO i=1,nagb
        ii=ii+i
        ERR=REAL(SQRT(ABS(cgmat(ii))))
        IF(cgmat((i*i+i)/2) < 0.0) ERR=-ERR
        gcor=0.0
        IF(cgmat(ii)*diag(i) > 0.0) THEN
            !         global correlation
            gcor=REAL(SQRT(ABS(1.0-1.0/(cgmat(ii)*diag(i)))))
        END IF
        IF(i <= 25.OR.nagb-i <= 25) THEN
            patext='  '
            IF(nlnpa(i) /= 0) patext='nl'
            IF(icnlim >= 0) WRITE(*,102) i,patext,  &
                pparm(i),pparm(i)+dparm(i),dparm(i),bgvec(i),ERR,gcor
        END IF
        par(i)=REAL(pparm(i)+dparm(i))    ! copy of result to array in argument
    END DO
    DO i=1,ncs                   ! constraints
        IF(i == 1.AND.icnlim >= 0) WRITE(*,*) '                      '
        dsum=0.0
        DO j=1,nagb
            dsum=dsum+adercs(nagb*(i-1)+j)*(pparm(j)+dparm(j))
        END DO
        IF(icnlim >= 0) WRITE(*,106) i,dsum,arhs(i),dsum-arhs(i)
    END DO
    IF(icnlim >= 0) THEN
        WRITE(*,*) '                   '
        WRITE(*,*) '  Total',loctot,' local fits,',locrej,' rejected.'
        WRITE(*,*) '  Histogram of RMS:',  &
            '  Total',nhist,' entries with',mhist(51), ' overflows'
        CALL pxhist(mhist,50,0.0,10.0) ! histogram printout
    END IF
    lthist=0
    DO i=1,51
        lthist=lthist+lhist(i)
    END DO
    IF(icnlim >= 0) THEN
        WRITE(*,*) '  Residual histogram:',  &
            '  Total',lthist,' entries with',lhist(51), ' overflows'
        CALL pxhist(lhist,50,-5.0,+5.0) ! histogram printout
        WRITE(*,199)
    END IF
199 FORMAT( '                                                       '/  &
        '              *   o   o                        o       '/  &
        '                  o   o                        o       '/  &
        '   o ooooo    o   o   o    oo   ooo    oo    ooo   oo  '/  &
        '    o  o  o   o   o   o   o  o  o  o  o  o  o  o  o  o '/  &
        '    o  o  o   o   o   o   oooo  o  o  oooo  o  o  oooo '/  &
        '    o  o  o   o   o   o   o     ooo   o     o  o  o    '/  &
        '    o  o  o   o   oo  oo   oo   o      oo    ooo   oo   ending.'/  &
        '                                o                      ')
101 FORMAT(1X,'   I         initial       final      differ',  &
        '     lastcor    Error glcor'/  &
        1X,' ---     ----------- ----------- -----------',  &
        ' ----------- -------- -----')
102 FORMAT(1X,i4,1X,a2,1X,4F12.5,f9.5,f6.3)
106 FORMAT(' Constraint',i2,'     Sum - RHS =',g12.5,' -',g12.5, ' = ',g12.5)
END SUBROUTINE fitglo

!> Return error for parameter I.
FUNCTION errpar(i)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    errpar=0.0
    IF(i <= 0.OR.i > nagb) RETURN
    ii=(i*i+i)/2
    errpar=REAL(SQRT(ABS(cgmat(ii))))
    IF(cgmat(ii) < 0.0) errpar=errpar
END FUNCTION errpar

!> Return correlation between parameters I and J.
FUNCTION corpar(i,j)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    corpar=0.0
    IF(i <= 0.OR.i > nagb) RETURN
    IF(j <= 0.OR.j > nagb) RETURN
    IF(i == j)              RETURN
    ii=(i*i+i)/2
    jj=(j*j+j)/2
    k=MAX(i,j)
    ij=(k*k-k)/2+MIN(i,j)
    ERR=REAL(SQRT(ABS(cgmat(ii)*cgmat(jj))))
    IF(ERR /= 0.0) corpar=REAL(cgmat(ij))/ERR
END FUNCTION corpar

!> Print result on file.
SUBROUTINE prtglo(lun)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    CHARACTER (LEN=2) :: patext
    !     ...
    lup=lun
    IF(lup == 0) lup=6
    WRITE(lup,*) '         Result of fit for global parameters'
    WRITE(lup,*) '         ==================================='
    WRITE(lup,101)
    ii=0
    DO i=1,nagb
        ii=ii+i
        ERR=REAL(SQRT(ABS(cgmat(ii))))
        IF(cgmat(ii) < 0.0) ERR=-ERR
        gcor=0.0
        IF(cgmat(ii)*diag(i) > 0.0) THEN
            !         global correlation
            gcor=REAL(SQRT(ABS(1.0-1.0/(cgmat(ii)*diag(i)))))
        END IF
        patext='  '
        IF(nlnpa(i) /= 0) patext='nl'
        WRITE(lup,102) i,patext,  &
            pparm(i),pparm(i)+dparm(i),dparm(i),bgvec(i),ERR,gcor
    END DO
101 FORMAT(1X,'   I         initial       final      differ',  &
        '     lastcor    Error glcor'/  &
        1X,' ---     ----------- ----------- -----------',  &
        ' ----------- -------- -----')
102 FORMAT(1X,i4,1X,a2,1X,4F12.5,f9.5,f6.3)
END SUBROUTINE prtglo


!> Obtain solution of a system of linear equations with symmetric matrix and the inverse.
!!
!!                    - - -
!!        CALL SPMINV(V,B,N,NRANK,...,...)      solve  V * X = B
!!                    - -   -----
!!
!!           V = symmetric N-by-N matrix in symmetric storage mode
!!               V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, . . .
!!               replaced by inverse matrix
!!           B = N-vector, replaced by solution vector
!!
!!     DIAG(N) =  double precision scratch array
!!     FLAG(N) =  logical scratch array
!!
!!     Method of solution is by elimination selecting the  pivot  on  the
!!     diagonal each stage. The rank of the matrix is returned in  NRANK.
!!     For NRANK ne N, all remaining  rows  and  cols  of  the  resulting
!!     matrix V and the corresponding elements of  B  are  set  to  zero.
SUBROUTINE spminv(v,b,n,nrank,diag,flag)

    DOUBLE PRECISION :: v(*),b(n),diag(n),vkk,vjk,eps
    LOGICAL :: flag(*)
    PARAMETER     (eps=1.0D-10)
    !     ...
    DO i=1,n
        flag(i)=.true.             ! reset flags
        diag(i)=ABS(v((i*i+i)/2))  ! save abs of diagonal elements
    END DO
    nrank=0
    DO i=1,n                    ! start of loop
        k  =0
        jj =0
        kk =0
        vkk=0.0D0
        DO j=1,n                   ! search for pivot
            jj=jj+j
            IF(flag(j)) THEN          ! not used so far
                IF(ABS(v(jj)) > MAX(ABS(vkk),eps*diag(j))) THEN
                    vkk=v(jj)           ! pivot (candidate)
                    k  =j               ! index of pivot
                    kk =jj              ! index of diagonal element
                END IF
            END IF
        END DO
        IF(k /= 0) THEN            ! pivot found
            nrank=nrank+1           ! increase rank and ...
            flag(k)=.false.         ! ... reset flag
            vkk    =1.0/vkk
            v(kk)  =-vkk
            b(k)   =b(k)*vkk
            jk     =kk-k
            jl     =0
            DO j=1,n                ! elimination
                IF(j == k) THEN
                    jk=kk
                    jl=jl+j
                ELSE
                    IF(j < k) THEN
                        jk=jk+1
                    ELSE
                        jk=jk+j-1
                    END IF
                    vjk  =v(jk)
                    v(jk)=vkk*vjk
                    b(j) =b(j)-b(k)*vjk
                    lk   =kk-k
                    DO l=1,j
                        jl=jl+1
                        IF(l == k) THEN
                            lk=kk
                        ELSE
                            IF(l < k) THEN
                                lk=lk+1
                            ELSE
                                lk=lk+l-1
                            END IF
                            v(jl)=v(jl)-v(lk)*vjk
                        END IF
                    END DO
                END IF
            END DO
        ELSE
            DO k=1,n
                IF(flag(k)) THEN
                    b(k)=0.0D0       ! clear vector element
                    DO j=1,k
                        IF(flag(j)) v((k*k-k)/2+j)=0.0D0  ! clear matrix row/col
                    END DO
                END IF
            END DO
            GO TO 10
        END IF
    END DO             ! end of loop
    10   DO ij=1,(n*n+n)/2
        v(ij)=-v(ij)      ! finally reverse sign of all matrix elements
    END DO
END SUBROUTINE spminv

!> Similarity operation A*V*A^t.
!!
!!     multiply symmetric N-by-N matrix from the left with general M-by-N
!!     matrix and from the right with the transposed of the same  general
!!     matrix  to  form  symmetric  M-by-M   matrix   (used   for   error
!!     propagation).
!!
!!                    - -   - -                                   T
!!        CALL SPAVAT(V,A,W,N,M)         W   =   A   *   V   *   A
!!                        -             M*M     M*N     N*N     N*M
!!
!!        where V = symmetric N-by-N matrix
!!              A = general N-by-M matrix
!!              W = symmetric M-by-M matrix
SUBROUTINE spavat(v,a,w,n,m)

    DOUBLE PRECISION :: v,a,w,cik
    DIMENSION v(*),a(*),w(*)
    !     ...
    DO i=1,(m*m+m)/2
        w(i)=0.0                ! reset output matrix
    END DO
    il=-n
    ijs=0
    DO i=1,m                 ! do I
        ijs=ijs+i-1             !
        il=il+n                 !
        lkl=0                   !
        DO k=1,n                !   do K
            cik=0.0D0              !
            lkl=lkl+k-1            !
            lk=lkl                 !
            DO l=1,k               !     do L
                lk=lk+1               !     .
                cik=cik+a(il+l)*v(lk) !     .
            END DO                 !     end do L
            DO l=k+1,n             !     do L
                lk=lk+l-1             !     .
                cik=cik+a(il+l)*v(lk) !     .
            END DO                 !     end do L
            jk=k                   !
            ij=ijs                 !
            DO j=1,i               !     do J
                ij=ij+1               !     .
                w(ij)=w(ij)+cik*a(jk) !     .
                jk=jk+n               !     .
            END DO                 !     end do J
        END DO                  !   end do K
    END DO                   ! end do I
END SUBROUTINE spavat

!> Multiply general M-by-N matrix A and N-vector X.
!!
!!                   - -   - -
!!        CALL  SPAX(A,X,Y,M,N)          Y   :=   A   *    X
!!                       -               M       M*N       N
!!
!!        where A = general M-by-N matrix (A11 A12 ... A1N  A21 A22 ...)
!!              X = N vector
!!              Y = M vector
SUBROUTINE spax(a,x,y,m,n)

    DOUBLE PRECISION :: a(*),x(*),y(*)
    !     ...
    ij=0
    DO i=1,m
        y(i)=0.0D0
        DO j=1,n
            ij=ij+1
            y(i)=y(i)+a(ij)*x(j)
        END DO
    END DO
END SUBROUTINE spax

!> Print X histogram.
SUBROUTINE pxhist(inc,n,xa,xb)
    PARAMETER (maxpl=70)
    INTEGER :: inc(n),num(maxpl)
    PARAMETER (np=maxpl/10+1)
    REAL :: p(np)
    EQUIVALENCE (nval,fval)
    CHARACTER (LEN=1) :: text(10)*130  ! 10 rows of X's
    CHARACTER (LEN=4) :: echar
    CHARACTER (LEN=8) :: xchar(maxpl)
    CHARACTER (LEN=1) :: chn(0:10)*1
    DATA chn/'0','1','2','3','4','5','6','7','8','9','0'/
    !     ...
    IF(n < 10) RETURN
    ntyp=0
    DO i=1,n
        IF(inc(i) /= 0) THEN
            nval=inc(i)
            IF(ABS(fval) < 1.0E-30) THEN
                ntyp=ntyp+1  ! integer
            ELSE
                ntyp=ntyp-1  ! floating point
            END IF
        END IF
    END DO
    m=1        ! M bins combined to 1
10  nred=n/m
    IF(nred > maxpl) THEN
        m=m+m   ! M = power of 2
        GO TO 10
    END IF
    nred=n/m   ! reduced number of bins
    nk=0
    IF(ntyp >= 0) THEN               ! integer
        DO i=1,n,m ! add M bins together
            mum=0
            DO k=i,MIN(i+m-1,n)
                mum=mum+inc(k)
            END DO
            nk=nk+1
            num(nk)=mum  ! copy to array NUM(.)
        END DO
    ELSE                            ! floating point
        DO i=1,n,m ! add M bins together
            sum=0
            DO k=i,MIN(i+m-1,n)
                nval=inc(k)
                sum=sum+fval
            END DO
            nk=nk+1
            num(nk)=INT(sum+0.5)  ! copy to array NUM(.)
        END DO
    END IF
    inmax=1        ! find maximum bin
    DO i=1,nk
        IF(num(i) > num(inmax)) inmax=i
    END DO
    inh=num(inmax) ! maximum bin content
    idiv=1+(inh-1)/10     ! X equivalent
    nr=inh/idiv    ! number of X lines
    DO l=1,nr
        text(l)=' '   ! blank text line
    END DO
    DO k=1,nred
        lr=num(k)/idiv
        IF(lr /= 0) THEN
            DO l=1,lr
                text(l)(k:k)='X'
            END DO
        ELSE
            IF(num(k) /= 0) text(1)(i:i)='.'
        END IF
    END DO
    DO l=nr,1,-1  ! print X's
        WRITE(*,103) text(l)(1:nred)
    END DO
    n10=1+(nred-1)/10
    text(1)=' '
    DO i=1,n10
        iup=MIN(10*i,nred)
        text(1)(10*i-9:iup)='----+----+'
        IF(iup == 10*i) text(1)(iup:iup)=chn(i)
    END DO
    WRITE(*,103) text(1)(1:nred)
    nap=nred/10+1
    IF(nap*10-10 > nred) nap=nap-1
    DO l=1,nap
        p(l)=xa+FLOAT(10*l-10)*(xb-xa)/FLOAT(nred)
    END DO
    CALL chfmt(p,nap,xchar,echar)
    WRITE(*,104) (xchar(l),l=1,nap)
    IF(echar /= ' ') WRITE(*,105) (echar,l=1,nap)
    WRITE(*,103) ' '
    nt=0
    DO k=1,10
        text(k)=' '
        DO i=1,nred
            IF(num(i) /= 0) THEN
                nt=k
                l=MOD(num(i),10)
                num(i)=num(i)/10
                text(k)(i:i)=chn(l)
            END IF
        END DO
    END DO
    DO l=nt,1,-1
        WRITE(*,103) text(l)(1:nred)
    END DO
    WRITE(*,103) ' '
103 FORMAT(5X,a)
104 FORMAT(12(a8,2X)/)
105 FORMAT(3X,12(a4,6X))
END SUBROUTINE pxhist

!> Prepare printout of array of real numbers as character strings.
!!
!!                  - -
!!       CALL CHFMT(X,N,XCHAR,ECHAR)
!!                      ----- -----
!!     where X( )     = array of n real values
!!           XCHAR( ) = array of n character*8 variables
!!           ECHAR    = character*4 variable
!!
!!     CHFMT converts an array of  n  real  values  into n character*  8
!!     variables (containing the values as text) and  a  common  exponent
!!     for printing. unneccessary zeros are suppressed.
!!
!!
!!     example: x(1)=1200.0, x(2)=1700.0 with n=2 are converted to
!!               xchar(1)='  1.2   ', xchar(2)='  1.7   ', echar='e 03'
SUBROUTINE chfmt(x,n,xchar,echar)

    REAL :: x(*)
    CHARACTER (LEN=1) :: ch(10)
    CHARACTER (LEN=4) :: echar
    CHARACTER (LEN=8) :: xchar(n),fxf,sxf,null
    DATA null/'00000000'/,ch/'0','1','2','3','4','5','6','7','8','9'/
    !     ...
    !     determine factor fc, so that fc*xmax is 5 digit number
    ip=0
    xm=0.0
    DO i=1,n
        xm=AMAX1(xm,ABS(x(i)))
    END DO
    IF(xm /= 0.0) THEN
        jp=104-IFIX(ALOG10(ABS(xm))+100.04)
        fc=10.0**jp
    ELSE
        jp=5
        fc=1.0
    END IF
    !     store digits as characters and find jm = first nonzero digit
    im=6
    DO i=1,n
        fxf=null
        ij=INT(fc*ABS(x(i))+0.5)
        jm=6
        IF(ij /= 0) THEN
            DO j=1,5
                jn=MOD(ij,10)
                ij=ij/10
                IF(jn /= 0.AND.jm == 6) jm=j
                fxf(j:j)=ch(jn+1)
            END DO
            im=MIN(im,jm)
        END IF
        xchar(i)=fxf
    END DO
    jm=im
    !     determine exponent as a multiple of 3
32  IF(jp < 1) THEN
        jp=jp+3
        ip=ip+3
        GO TO 32
    END IF
34  IF(jp > jm+4.OR.jp >= 8) THEN
        jp=jp-3
        ip=ip-3
        GO TO 34
    END IF
    !     loop to convert to print format
    ja=MIN(jm,jp)
    jb=MAX(6,jp+1)
    DO  i=1,n
        fxf=xchar(i)
        sxf=' '
        ib=7+(jb-ja)/2
        DO  j=ja,jb
            IF(fxf(j:j) /= ch(1)) GO TO 70
            IF(j > jp+1) GO TO 50
            IF(fxf /= null.OR.j >= jp) GO TO 70
            ib=ib-1
            CYCLE
            50 DO k=j,jb
                IF(fxf(k:k) /= ch(1)) GO TO 70
            END DO
            CYCLE
            !     insert digit
70          ib=ib-1
            sxf(ib:ib)=fxf(j:j)
            IF(j == jp) THEN
                !     insert decimal dot
                ib=ib-1
                sxf(ib:ib)='.'
            END IF
        END DO
        !     insert - sign
        IF(x(i) < 0.0) sxf(ib-1:ib-1)='-'
        xchar(i)=sxf
    END DO
    !     prepare print format for exponent
    echar=' '
    IF(ip /= 0) THEN
        echar='E 0 '
        IF(ip <= 0) THEN
            echar(2:2)='-'
            ip=IABS(ip)
        END IF
        j=MOD(ip,10)
        echar(4:4)=ch(j+1)
        ip=(ip-j)/10
        IF(ip /= 0) THEN
            j=MOD(ip,10)
            echar(3:3)=ch(j+1)
        END IF
    END IF
END SUBROUTINE chfmt

!> Return limit in chi^2/ND for N sigmas (N=1, 2 or 3).
FUNCTION chindl(n,nd)
    REAL :: pn(3),sn(3),table(30,3)
    DATA pn/0.31731,0.0455002785,2.69985E-3/         ! probabilities
    DATA sn/0.47523,1.690140,2.782170/
    DATA table/ 1.0000, 1.1479, 1.1753, 1.1798, 1.1775, 1.1730, 1.1680, 1.1630,  &
        1.1581, 1.1536, 1.1493, 1.1454, 1.1417, 1.1383, 1.1351, 1.1321,  &
        1.1293, 1.1266, 1.1242, 1.1218, 1.1196, 1.1175, 1.1155, 1.1136,  &
        1.1119, 1.1101, 1.1085, 1.1070, 1.1055, 1.1040,  &
        4.0000, 3.0900, 2.6750, 2.4290, 2.2628, 2.1415, 2.0481, 1.9736,  &
        1.9124, 1.8610, 1.8171, 1.7791, 1.7457, 1.7161, 1.6897, 1.6658,  &
        1.6442, 1.6246, 1.6065, 1.5899, 1.5745, 1.5603, 1.5470, 1.5346,  &
        1.5230, 1.5120, 1.5017, 1.4920, 1.4829, 1.4742,  &
        9.0000, 5.9146, 4.7184, 4.0628, 3.6410, 3.3436, 3.1209, 2.9468,  &
        2.8063, 2.6902, 2.5922, 2.5082, 2.4352, 2.3711, 2.3143, 2.2635,  &
        2.2178, 2.1764, 2.1386, 2.1040, 2.0722, 2.0428, 2.0155, 1.9901,  &
        1.9665, 1.9443, 1.9235, 1.9040, 1.8855, 1.8681/
    !     ...
    IF(nd < 1) THEN
        chindl=0.0
    ELSE
        m=MAX(1,MIN(n,3))         ! 1, 2 or 3 sigmas
        IF(nd <= 30) THEN
            chindl=table(nd,m)     ! from table
        ELSE                      ! approximation for ND > 30
            chindl=(sn(m)+SQRT(FLOAT(nd+nd-1)))**2/FLOAT(nd+nd)
        END IF
    END IF
END FUNCTION chindl

!> Get matrix information out.
SUBROUTINE fitmut(nvec,vec)
    !     ------------------------------------------------------------------
    !     Basic dimension parameters
    PARAMETER (mglobl=1400,mlocal=10,nstore=10000,mcs=10) ! dimensions
    PARAMETER (mgl=mglobl+mcs)                     ! derived parameter
    !     derived parameters
    PARAMETER (msymgb  =(mglobl*mglobl+mglobl)/2, msym    =(mgl*mgl+mgl)/2,  &
        msymlc=(mlocal*mlocal+mlocal)/2, mrecta= mglobl*mlocal,  &
        mglocs= mglobl*mcs, msymcs= (mcs*mcs+mcs)/2 )
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
        nlnpa(mglobl),nstdev,cfactr,icnpr,icnlim, indnz(mglobl),indbk(mglobl)
    !     ------------------------------------------------------------------
    REAL :: vec(*)
    !     ...
    nvec=(nagb*nagb+nagb)/2+nagb
    DO i=1,nvec-nagb
        vec(i)=REAL(cgmat(i))
    END DO
    DO i=1,nagb
        vec(nvec-nagb+i)=REAL(bgvec(i))
    END DO
    WRITE(*,*) 'FITMUT ',nvec,(vec(i),i=1,nvec)
    RETURN
    ENTRY fitmin(nvec,vec)
    !     insert information
    IF(nvec /= (nagb*nagb+nagb)/2+nagb) THEN
        WRITE(*,*) ' Wrong dimensions in FITMUT/FITMIN'
        WRITE(*,*) ' Argument NVEC =',nvec
        WRITE(*,*) ' Expected NVEC =',(nagb*nagb+nagb)/2+nagb
        STOP 'FITMIN'
    END IF
    DO i=1,nvec-nagb
        cgmat(i)=vec(i)
    END DO
    DO i=1,nagb
        bgvec(i)=vec(nvec-nagb+i)
    END DO
END SUBROUTINE fitmut











