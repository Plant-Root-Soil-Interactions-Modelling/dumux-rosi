! ==============================================================================
! Source file WATER    |||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
!> main subroutine WATER calculates water flow in the soil 
      SUBROUTINE WATER(t,dt,dtOpt,tOld,ReDo,IAD,IADN,IADD)
      USE ParamData, ONLY: lChem,maxbnd,maxplant,iter,iter_tot
      USE PlntData, ONLY:TotSur
      USE tmctrl, ONlY: dtMin
      USE GridData, ONLY : epslonPH,epslonR,epslonWC,epslonS,RelEps,factorRelEps,&
           nPt,Axy,Dxy,Exy,sinkOld,sink,itMax,betac,wn,betac,dxgrid,dygrid,dzgrid,nElm,&
           subN,deter,sink,sink_cube,sink_cubeOld,HElm,HElmOld,betac_cube,betac_cube2
      USE MatData
      USE DoussanMat, ONLY : PHr,PHrOld,sinkRold,SinkR,PHrtemp,PHstemp,delta2, &
           switchcriterion,savelast,stressBC,PHs,delta2old,hx_min,nplant
      USE RootData, ONLY :nrec,lCou,lDou,lFed,ldJvL,ltemp
      USE WatFun
      USE SolData
      USE RhizoData
      IMPLICIT NONE

      REAL(sp),INTENT(in):: told
      REAL(sp),INTENT(inout)::dt,t,dtopt
      LOGICAL,INTENT(out) :: ReDo

      REAL(sp):: tol_R,tol_S,tol_PH
      REAL(sp):: deltaPHr,deltaPHint,deltaS,BCr
      REAL(dp):: deltaPHs, maxPhrtot, maxPhr(1:maxplant)
      !REAL(sp) :: t0,t1!A1(maxbnd,maxnod),B1(maxnod)
      INTEGER(ap) :: IAD(maxbnd,nPt),IADN(nPt),IADD(nPt),iter_root,iter_soil
      INTEGER(ap) :: ipl,BCtp,i
      LOGICAL :: ItCrit,stressBC_old,lOrt,Celia=.false.
      REAL(dp), ALLOCATABLE,DIMENSION (:) ::sol
      ALLOCATE(sol(nPt),B(nPt))
      !> \param t current time
	  !> \param dt time step size
	  !> \param dtOpt optimal time step size calculated in the time step control TmCont
	  !> \param tOld time before last time step
	  !> \param ReDo logical to rerun the model from a given time step
	  !> \param IAD adjacency matrix (nodal connections)
	  !> \param IADN number of adjacent nodes in IAD (self-inclusive)
	  !> \param IADD position of diagonal in adjacency matrix
!----------------------- initialisation ------------------------------------
      iter_root=0
      ReDo=.FALSE.
      theta_old=theta
      IF (lDou) THEN
         DO ipl=1,nplant
            PHrOld(1:nrec+1,ipl)=PHr(1:nrec+1,ipl) !changed MJ june09
            sinkRold(0:nrec,ipl)=sinkR(0:nrec,ipl)
         ENDDO
      ENDIF
      iter_soil = iter
      tol_PH=epslonPH
      tol_R=epslonR
      tol_S=epslonS
      iter = 0
! ----------------------- Start of FEM iteration loop -------------------------
    250 CONTINUE
      
      IF (lDou) stressBC_old=stressBC
! NonEquilibrium in the rhizosphere
      IF (lRhizo .AND. RhizoModel .EQ. 3) THEN ! Only Relevant for the Dynamic Scenario 
        DO i=1,nPt
            IF (Rnorm(i) .NE. 1.) THEN
                CALL thtNonEq(hnew(i), dt, i)
            ELSE
                thetaNonEq(i) = 0.0
                thetaTot(i) = theta(i)
                hEqRhizo(i) = hNew(i)
                tauTht(i) = 1.0
            ENDIF
        ENDDO
      ENDIF
      CALL SetMat(1)
      IF (Celia) THEN
        ! PRINT *,'celia'!
         CALL ResetCelia(dt,IAD,IADN,IADD)
      ELSE
         CALL Reset(dt,IAD,IADN,IADD)
      ENDIF
      
      CALL Dirich(IADD)
      CALL SolveIlut(sol)
! test for convergence:
      ItCrit=.FALSE.!true if no convergence for soil or root
      hTemp=hNew
      hNew=sol
     
      iter=iter+1
      iter_tot=iter_tot+1

      IF (RelEps) THEN
         tol_PH=epslonPH
         IF (epslonPH .LT. ABS(MAXVAL(hNew(1:nPt)/factorRelEps))) tol_PH= ABS(MAXVAL(hNew(1:nPt)/factorRelEps))
      ENDIF
      deltaPHs=MAXVAL(ABS(hNew(1:nPt)-hTemp(1:nPt)))
!print *, 'errPH=',deltaPHs
      IF (deltaPHs.GT.tol_PH) THEN !.OR.deltath.gt.eps_WC concelled 24.06.09
         WRITE (*,'(''x'',$)')
         ItCrit=.TRUE.
      ENDIF
      IF (.NOT.(itCrit)) THEN !in case there is convergence for soil (itCrit==.false.)
         sink_cubeOld=sink_cube
         IF (lDou) THEN
            DO ipl=1,nplant
               PHrTemp(1:nrec+1,ipl)=PHr(1:nrec+1,ipl)
            ENDDO
            CALL SolveRoot(t,dt,.FALSE.,iter_root)
!relative or absolute error
            IF (RelEps) THEN 
               DO ipl=1,nplant
                  maxPhr(ipl)=MAXVAL(PHr(1:nrec+1,1:nplant)/factorRelEps)
               ENDDO
               maxPhrtot=MAXVAL(maxPhr(1:nplant))
               tol_R=epslonR
               IF (epslonR .LT. ABS(maxPhrtot)) tol_R= ABS(maxPhrtot) ! we take the largest tolerance
            ENDIF

!check for all nodes if h_value is lower than relative or absolute critetrium
            maxPhr=0
            DO ipl=1,nplant
               maxPhr(ipl)=MAXVAL(ABS(PHr(1:nrec+1,ipl)-PHrTemp(1:nrec+1,ipl)))
            ENDDO
            deltaPHr=MAXVAL(maxPhr(1:nplant))
            deltaS=MAXVAL(ABS(sink_cubeOld-sink_cube))
            IF ((deltaPHr.GT.tol_R).OR.(deltaS.GT.tol_S)) THEN
               WRITE(*,'(''+'',$)')
               ItCrit=.TRUE.! no convergence for root-> rerun soil and root.
            ENDIF
         ELSEIF (lFed.OR.lCou) THEN
            IF (ldJvL) THEN
               HElmOld=HElm
            ENDIF
            CALL SetBCroot(t,BCr,BCtp)
            CALL SetSnk(t,BCr,BCtp)
            IF (ldJvL) THEN
               deltaPHint=MAXVAL(ABS(HElmOld-HElm))
               IF (deltaPHint.GT.tol_R) THEN
                  WRITE(*,'(''*'',$)')
                  ItCrit=.TRUE.! no convergence for soil-root interface water potential-> rerun soil and root.
               ENDIF
            ENDIF
            deltaS=MAXVAL(ABS(sink_cubeOld-sink_cube))
            IF (deltaS.GT.tol_S) THEN
               WRITE(*,'(''+'',$)')
               ItCrit=.TRUE.! no convergence for root-> rerun soil and root.
            ENDIF
         ENDIF
      ENDIF

      IF ((savelast) .AND. (.NOT.(iTCrit)) ) RETURN

      IF (ItCrit) THEN
         IF (iter.LT.itMax) THEN
            GOTO 250
         ELSE IF (dt.LE.dtMin) THEN
            STOP 'Cannot converge, have to stop here. Sorry.'
            !WRITE (*,'(''Cannot converge, have to stop here. Sorry.'')')
            !savelast=.TRUE.
            !IF (savelast) RETURN
         ELSE
            WRITE (*,'(''!'',$)')
            WRITE (*,'('' Iter. tot.: '',I5,$)') iter_tot
            hNew =hOld
            hTemp=hOld
         
            dt=MAX(dt/3._sp,dtMin) !changed javaux july 09
            dtOpt=dt
            t=tOld+dt
            IF (lDou) THEN
               stressBC=stressBC_old
               PHr(1:nrec+1,:)=PHrOld(1:nrec+1,:)
               !print*,'hier water',phr(1:2,1)
            ENDIF
            ReDo=.TRUE.
            GOTO 200
         ENDIF
     ELSE
!iteration is OK just make a last update of soil WC& PH with the right sink term
     ENDIF
     CALL WatInf(dt)
     CALL SetMat(1)
    
! --------------------- End of FEM iteration loop -----------------------------
200  DEALLOCATE(sol,B)
     RETURN
    
     END SUBROUTINE WATER
!*******************************************************************************
!> solve the linear equation system Ax=b - preconditioner and iterative solver can be chosen here    
	SUBROUTINE SOLVEILUT(sol)
      USE MatData
      USE GridData, ONLY : nPt

      INTEGER(ap) ::lfil,maxits=1000,ierr
      REAL(dp) :: droptol
      INTEGER(ap) ::ipar(16),nwk
      REAL(dp):: fpar(16)
      REAL(dp), INTENT(out) ::sol(nPt)

      INTEGER(ap) , ALLOCATABLE, DIMENSION(:) :: jwork,ju
      REAL(dp), ALLOCATABLE,DIMENSION (:) ::work,xran,ww
      REAL(dp), ALLOCATABLE,DIMENSION (:,:) ::vv
      EXTERNAL cgnr, fom, ilut,ilu0,ulid
      EXTERNAL cg,bcg,dbcg,bcgstab,tfqmr,gmres,fgmres,dqgmres
    
      ALLOCATE (jwork(2*nPt))
      ALLOCATE (xran(nPt))
      ALLOCATE (work(nPt+1))
      ALLOCATE (vv(nPt,20))
      ALLOCATE (ww(nPt*40))
      ALLOCATE (ju(nPt))
      ALLOCATE (alu(numnz))
      ALLOCATE (jlu(numnz))
	  
	  !> \param sol solution vector of the linear equation system -> pressure head

      !set the parameters for the iterative solvers
      ipar(2) = 1      !1 is left preconditioning,2 right,3 both, 0 = no precon
      ipar(3) = 1      !1 -> convergence test based on residualnorms; 2->based on change in approximate solution
      ipar(4) = nPt*40      !# elements in array w
      ipar(5) = 16      !size of krylov subspace
      ipar(6) = maxits      !max number of matrix-vector operations
      fpar(1) = 1.0E-8   !relative tolerance
      fpar(2) = 1.0E-12  !absolute tolerance
      lfil = 1
      droptol =1.0E-2
      nwk = 2*numNZ

      !preconditioner
      !IF (MOD(time_step,10).EQ.0) THEN
         CALL ilut (nPt,A_sparse,JCOL,IROW,lfil,droptol,alu,jlu,ju,nwk,work,jwork,ierr)
      !END IF
      time_step = time_step +1
      !solve
      xran=0.0 !initial guess
      CALL runrc(nPt,B,sol,ipar,fpar,ww,xran,A_sparse,JCOL,IROW,alu,jlu,ju,bcgstab)

      DEALLOCATE (jwork)
      DEALLOCATE (xran)
      DEALLOCATE (work)
      DEALLOCATE (vv)
      DEALLOCATE (ww)
      DEALLOCATE (ju)
      DEALLOCATE (alu)
      DEALLOCATE (jlu)

     END SUBROUTINE SOLVEILUT
!**************************************************************************
!>     the actual tester. It starts the iterative linear system solvers
!>     with a initial guess suppied by the user.
!>
!>     The structure {au, jau, ju} is assumed to have the output from
!>     the ILU* routines in ilut.f.   
SUBROUTINE runrc(n,rhs,sol,ipar,fpar,wk,guess,a,ja,ia,&
    au,jau,ju,solver)
      IMPLICIT NONE
      INTEGER n,ipar(16),ia(n+1),ja(*),ju(*),jau(*)
      REAL*8 fpar(16),rhs(n),sol(n),guess(n),wk(*),a(*),au(*)
      EXTERNAL solver

!c-----------------------------------------------------------------------
!c     local variables
      INTEGER i, its
      REAL*8 res, dnrm2
!c     external dtime
      EXTERNAL dnrm2
      SAVE its,res

!ipar(2) can be 0, 1, 2, please don´t use 3
!      IF (ipar(2).GT.2) THEN
!         PRINT *, 'I can not do both left and right preconditioning.'
!         RETURN
!      ENDIF

!normal execution

      its = 0
      res = 0.0D0
      DO i = 1, n
         sol(i) = guess(i)
      ENDDO
      ipar(1) = 0

  10  CALL solver(n,rhs,sol,ipar,fpar,wk)

!output the residuals

      IF (ipar(7).NE.its) THEN
         its = ipar(7)
      ENDIF
      res = fpar(5)
!c
      IF (ipar(1).EQ.1) THEN
         CALL amux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
         GOTO 10
      ELSE IF (ipar(1).EQ.2) THEN
         CALL atmux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
         GOTO 10
      ELSE IF (ipar(1).EQ.3 .OR. ipar(1).EQ.5) THEN
         CALL lusol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
         GOTO 10
      ELSE IF (ipar(1).EQ.4 .OR. ipar(1).EQ.6) THEN
         CALL lutsol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
         GOTO 10
      ELSE IF (ipar(1).LE.0) THEN
         IF (ipar(1).EQ.0) THEN
 !            WRITE (*,'(''*'',$)')
            !Iterative solver has satisfied convergence test.
         ELSE IF (ipar(1).EQ.-1) THEN
            PRINT *, 'Iterative solver has iterated too many times.'
         ELSE IF (ipar(1).EQ.-2) THEN
            PRINT *, 'Iterative solver was not given enough work space.'
            PRINT *, 'The work space should at least have ', ipar(4),&
           ' elements.'
         ELSE IF (ipar(1).EQ.-3) THEN
            PRINT *, 'Iterative sovler is facing a break-down.'
         ELSE
            PRINT *, 'Iterative solver terminated. code =', ipar(1)
         ENDIF
      ENDIF
      RETURN
      END
!c-----end-of-runrc
  FUNCTION distdot(n,x,ix,y,iy)
      INTEGER n, ix, iy
      REAL*8 distdot, x(*), y(*), ddot
      EXTERNAL ddot
      distdot = ddot(n,x,ix,y,iy)
      RETURN
    END FUNCTION distdot
!*********************************************************************
!> sets the boundary conditions in the matrix and in the right-hand side vector 
!> of the linear equation system
    SUBROUTINE SetBC(t)
      USE typedef
      USE BoundData
      USE GridData
      USE SolData, ONLY: hnew,Kode
      USE CumData
      IMPLICIT NONE
      REAL(sp) :: h(maxnod),t,head,Qtot,wtot,m,b
      REAL(sp) , DIMENSION (:) :: volflw(1+(homogene-1)*mxBcCh),irrigflw(nBcPts)
      INTEGER(ap) :: i,jj,n
	  !> \param t current simulation time
	  
      Q = 0._sp
! calculate current value for Q-BC´s, if any:
      IF (nQbcCh.EQ.0) GOTO 25
      i=0
   22 i=i+1
      IF(Kode(iBCPt(i)).EQ.-1) THEN
         IF (t.GE.tQbcCh(nQbcCh)) THEN
            IF (homogene.EQ.1) THEN	
               volflw=Qbc(nQbcCh,1)
            ELSE
               volflw(1:mxBcCh)=Qbc(nQbcCh,:)
            ENDIF
         ELSE
            IF (t.GT.tQbcCh(i)) GOTO 22
            IF (homogene.EQ.1) THEN	
               volflw=Qbc(i-1,1)+(t-tQbcCh(i-1))/(tQbcCh(i)-tQbcCh(i-1))*(Qbc(i,1)-Qbc(i-1,1))
            ELSE
               volflw(1:mxBcCh)=Qbc(i-1,:)+(t-tQbcCh(i-1))/(tQbcCh(i)-tQbcCh(i-1))*(Qbc(i,:)-Qbc(i-1,:))
            ENDIF
         ENDIF
      ENDIF
! calculate current value for I-BC´s, if any:			
   25 IF (nIbcCh.EQ.0) GOTO 30
      irrigflw(:)=0
      ! loop over irrigators
      n=0
      DO i=1,nBcPts
         IF (Kode(iBcPt(i)).EQ.-3) THEN
            n=n+1      !irrigator number
            IF(t.GE.tIBcCh(n,nIBcCh)) THEN
               irrigflw(i)=Ibc(n,nIBcCh)
            ELSE
               jj=nIbcCh
    27         jj=jj-1
               IF((t.LT.tIBcCh(n,jj)).AND.(jj.GT.1)) GOTO 27
               m=(Ibc(n,jj+1)-Ibc(n,jj))/(tIBcCh(n,jj+1)-tIBcCh(n,jj))
               b=Ibc(n,jj)-tIBcCh(n,jj)*m
               irrigflw(i)=m*t+b
            ENDIF
         ENDIF
      ENDDO
! calculate current value for h-BC´s, if any:
   30 IF (nhbcCh.EQ.0) GOTO 40
      i=1
   33 i=i+1
      IF (i.GT.nhbcCh) THEN
         head=hbc(nhbcCh)
      ELSE
         IF (t.GT.thbcCh(i)) GOTO 33
         head=hbc(i-1)+(t-thbcCh(i-1))/(thbcCh(i)-thbcCh(i-1))*(hbc(i)-hbc(i-1))
      ENDIF
!calculate QBC*surface
   40 Qtot=0
      wtot=0
      DO i=1,nBCPts
         IF ((Kode(iBCPt(i)).EQ.+1).OR.(Kode(iBCPt(i)).EQ.+2)) THEN
            hnew(iBCPt(i))=head
         ELSEIF (Kode(iBCPt(i)).EQ.-1) THEN
            IF (homogene.EQ.1) THEN						
	        Q(iBCPt(i))=volflw(1)*Width(iBCPt(i))
            ELSE
               Q(iBCPt(i))=volflw(i)*width(iBCPt(i)) 
            ENDIF
	     Qtot=Qtot+Q(iBCPt(i))
            wtot=wtot+Width(iBCPt(i))
         ELSEIF (Kode(iBCPt(i)).EQ.-3) THEN
            Q(iBCPt(i))=irrigflw(i)*dxgrid*dygrid
            Qtot=Qtot+Q(iBCPt(i))
            wtot=wtot+Width(i)
         ENDIF   
      END DO
! calculate current value for solute transport boundary conditions, if any:
      IF (nCBnd1.EQ.0) GOTO 50
      i=0
   55 i=i+1
      IF (i.GE.nCBnd1) THEN
         cBound(1)=CBnd1(nCBnd1)
      ELSE
         IF (t.GT.tCBnd1(i)) GOTO 55
         cBound(1)=CBnd1(i-1)+(t-tCBnd1(i-1))/(tCBnd1(i)-tCBnd1(i-1))*(CBnd1(i)-CBnd1(i-1))
      ENDIF
   50 IF (nCBnd2.EQ.0) GOTO 60
      i=0
   66 i=i+1
      IF (i.GE.nCBnd2) THEN
         cBound(2)=CBnd2(nCBnd2)
      ELSE
         IF (t.GT.tCBnd2(i)) GOTO 66
         cBound(2)=CBnd2(i-1)+(t-tCBnd2(i-1))/(tCBnd2(i)-tCBnd2(i-1))*(CBnd2(i)-CBnd2(i-1))
      ENDIF
   60 RETURN
    END SUBROUTINE SetBC
!*********************************************************************************
!> calculates the width at the boundary faces similar to the SWMS_3D model    
	SUBROUTINE CalcWidthSWMS_new
      USE GridData
      USE DomData
      USE BoundData, ONLY: xqmax1,xqmin2,qfun
      USE SolData, ONLY: MatNum
      IMPLICIT NONE
      INTEGER(ap):: iB,ix,iy,iCord,iDiv,iRest,i,bord1,bord2
      REAL (sp):: WidthX,WidthY
      iB=nPt-(nx*ny)
      i=0
      DO iy=1,ny
         DO ix=1,nx
            i=i+1
            iB=iB+1
            IF(qfun.EQ.5) THEN  !for cylindrical/arbitrary setup              
                  !X direction
                  IF(ix.NE.1 .AND. ix.NE.nx) THEN
                     IF(MatNum(i-1).EQ.1 .AND. MatNum(i+1).EQ.1) THEN
                        WidthX =  xCol(ix+1)-xCol(ix-1)
                     ELSE IF(MatNum(i-1).NE.1) THEN !left border
                        WidthX = xCol(ix+1)-xCol(ix)
                     ELSE IF(MatNum(i+1).NE.1) THEN !right border
                        WidthX = xCol(ix)-xCol(ix-1)
                     END IF
                  ELSE IF(ix.EQ.1) THEN
                     WidthX = xCol(ix+1)-xCol(1)			
                  ELSE IF(ix.EQ.nx) THEN
                     WidthX = xCol(nx)-xCol(nx-1)
                  END IF
                  !Y direction
                  IF(iy.NE.1 .AND. iy.NE.ny) THEN
                     IF(MatNum(i-ny).EQ.1 .AND. MatNum(i+ny).EQ.1) THEN
                        WidthY =  yCol(iy+1)-yCol(iy-1)
                     ELSE IF(MatNum(i-ny).NE.1) THEN !lower border
                        WidthY = yCol(iy+1)-xCol(iy)
                     ELSE IF(MatNum(i+ny).NE.1) THEN !upper border
                        WidthY = yCol(iy)-yCol(iy-1)
                     END IF
                  ELSE IF(iy.EQ.1) THEN
                     WidthY = yCol(iy+1)-yCol(1)			
                  ELSE IF(iy.EQ.ny) THEN
                     WidthY = yCol(ny)-yCol(ny-1)
                  END IF

            ELSE               
               IF(ix.NE.1.AND.ix.NE.nx) THEN
                  IF(qfun.EQ.4) THEN  !for split setup introducing a second and third border at xqmax1,xqmin2
                     bord1=(xqmax1-xmin)/dxgrid
                     bord2=(xqmin2-xmin)/dxgrid+1
                     IF(ix.NE.bord1.AND.ix.NE.bord2) THEN
                        WidthX=(xCol(ix+1)-xCol(ix-1))
                     ELSE IF(ix.EQ.bord1) THEN
                        WidthX=xcol(bord1)-xcol(bord1-1)
                     ELSE IF(ix.EQ.bord2) THEN
                        WidthX=xcol(bord2+1)-xcol(bord2)
                     END IF
                  ELSE
                     WidthX=(xCol(ix+1)-xCol(ix-1))
                  END IF
               ELSE IF(ix.EQ.1) THEN
                  WidthX=(xCol(ix+1)-xCol(1))
                  IF(continu) WidthX=(xCol(2)-xCol(nx)+nx*dxGrid)				
               ELSE IF(ix.EQ.nx) THEN
                  WidthX=(xCol(nx)-xCol(nx-1))
                  IF(continu) WidthX=(xCol(1)-xCol(nx-1)+nx*dxGrid)
               END IF
               IF(iy.NE.1.AND.iy.NE.ny) THEN
                  WidthY=(yCol(iy+1)-yCol(iy-1))
               ELSE IF(iy.EQ.1) THEN
                  WidthY=(yCol(iy+1)-yCol(1))
                  IF(continu) WidthY=(yCol(2)-yCol(ny)+ny*dyGrid)
               ELSE IF(iy.EQ.ny) THEN
                  WidthY=(yCol(ny)-yCol(ny-1))
                  IF(continu) WidthY=(yCol(1)-yCol(ny-1)+ny*dyGrid)
               END IF
            END IF
            iCord=ix+iy
            iRest=MOD(iCord,2)
            iDiv=3
            IF(iRest.EQ.1) iDiv=6 
            Width(i)=WidthX*WidthY/iDiv !top surface
            iCord=ix+iy+nz
            iRest=MOD(iCord,2)
            iDiv=3
            IF(iRest.EQ.0) iDiv=6
            Width(iB)= WidthX*WidthY/iDiv !bottom surface
         END DO
      END DO
      RETURN
    END SUBROUTINE CalcWidthSWMS_new
!****************************************************************
!> calculates the soil conductivity - the used model is  defined in the input files    
	 SUBROUTINE SetMat(K)

      USE GridData, ONLY : nPt,Axy,Bxy,Dxy,Exy
      USE SolData
      USE WatFun
      USE RhizoData
      IMPLICIT NONE

      REAL(sp):: ci,cpi,Ti!,t0,t1
      REAL(sp):: sl,s2,him,hi2,hi1
      INTEGER(ap):: it,i,K,m
      IF (lRhizo) THEN
        CALL SetMatRhizo(K)
      ELSE
        !$OMP PARALLEL DO private(M,hi1,hi2,hiM,iT,Sl,ci,Cpi,Ti)
              DO i=1,nPt
        ! calculate nodal conductivity values:
      
                 IF(K.EQ.1) conO(i)=con(i)
                 M=MatNum(i)
                 hi1=hTemp(i)/Axy(i) !javaux
                 hi2=hNew(i)/Axy(i)!javaux
                 hiM=0.1_dp*hi1+0.9_dp*hi2
                 IF (soiltab) THEN	
                    ci=FKP_soiltab(hiM,M)
                 ELSEIF (lTab.AND.(hiM.GE.hTab(nTab)).AND.(hiM.LE.hTab(1))) THEN
                    iT=INT((LOG10(-hiM)-alh1)/dlh)+1
                    Sl=(ConTab(iT+1,M)-ConTab(iT,M))/(hTab(iT+1)-hTab(iT))
                    ci=ConTab(iT,M)+Sl*(hiM-hTab(iT))
                 ELSE
                    ci=FKP(hiM,par(:,M),i)
                 ENDIF
                 con(i)=ci*Bxy(i)
                 IF(K.EQ.0) conO(i)=con(i)

        ! calculate nodal capacity values:
                 hi1=hOld(i)/Axy(i)
                 hi2=hNew(i)/Axy(i)
                 hiM=(hi2+hi1)/2
                 IF (soiltab) THEN	
                    Cpi=FCP_soiltab(hiM,M)
                    Ti=Fth_soiltab(hiM,M)
                 ELSEIF (lTab.AND.(hiM.GE.hTab(nTab)).AND.(hiM.LE.hTab(1))) THEN
                    iT=INT((LOG10(-hiM)-alh1)/dlh)+1
                    Sl =(CapTab(iT+1,M)-CapTab(iT,M))/(hTab(iT+1)-hTab(iT))
                    S2 =(TheTab(iT+1,M)-TheTab(iT,M))/(hTab(iT+1)-hTab(iT))
                    Cpi=CapTab(iT,M)+Sl*(hiM-hTab(iT))
                    Ti=TheTab(iT,M)+S2*(hiM-hTab(iT))
                 ELSE
                    Cpi=FCP(hiM,par(:,M))
                    Ti=Fth(hiM,par(:,M),i)
                 ENDIF
                 Cap(i)=Cpi*Dxy(i)/Axy(i)
                 IF (soiltab) THEN
                    theta(i)=TheTab(nTab,M)*Exy(i)+(Ti-TheTab(nTab,M))*Dxy(i)
                 ELSE
                    theta(i)=par(2,M)*Exy(i)+(Ti-par(2,M))*Dxy(i)
                 ENDIF
         ENDDO
      ENDIF
!$OMP END PARALLEL DO
!call cpu_time(t1)

      RETURN
      END SUBROUTINE SetMat !MJ09
!************************************************************************************
!> SetMat for the Rhizosphere
!> Calculate soil conductivity and capacity
   SUBROUTINE SetMatRhizo(K)
        USE GridData, ONLY : nPt
        USE SolDAta
        USE WatFun
        USE RhizoData
        IMPLICIT NONE

        REAL(sp)    :: him,hi2,hi1
        INTEGER(ap) :: i,K,M
        M = MatNum(1)
        Do i=1,nPt
                !! In this section the conductivity is calculated
                IF(K .EQ. 1) conO(i)=con(i)
                hi1 = hTemp(i)
                hi2 = hNew(i)
                hiM = 0.1_dp*hi1 + 0.9_dp*hi2
                theta(i) = Fth(hiM,par(:,M),i)
                ! Calculate the initial thetaTot water contet
                IF (K.EQ. 0) CALL IniRhizoTheta(hiM,i)
                ! If rhizoStatic, calculate the ThetaTot
                IF (K .EQ. 1 .AND. RhizoModel .EQ. 2) Call IniRhizoTheta(hiM,i)
                con(i) = FKP(hiM,par(:,M),i)
                IF (K .EQ. 0) conO(i)=con(i)
                
                !! In this section the capacity is calculated
                !! For the capacity hiM is calculated differently. 
                hi1=hOld(i)
                hiM = (hi2+hi1)/2
                ! Capacity
                Cap(i) = FCP(hiM,par(:,M))
                ! The water content is recalculated based on the hiM
                theta(i) = Fth(hiM,par(:,M),i)
                IF (K .EQ. 0) CALL IniRhizoTheta(hiM,i)
                IF (K .EQ. 1 .AND. RhizoModel .NE. 1) CALL IniRhizoTheta(hiM,i)
        ENDDO
        RETURN
  END SUBROUTINE SetMatRhizo
      

!************************************************************************************
!> calculate the geometry, shape_functions etc needed for the finite element method   
!> takes care of splitting the cube into 5 tetrahedrals 
   SUBROUTINE CalcGeom
      USE GridData
      USE CumData
      USE DomData
      USE MatData
      IMPLICIT NONE

      INTEGER(ap) :: iE,iSE,i,j,k,l,ii,jj
      REAL(sp) :: xGridi,xGridj,xGridk,xGridl,yGridi,yGridj,yGridk,yGridl
      REAL(sp) :: Cxx,Czz,Cyy,Cxy,Cxz,Cyz,VE,ai(1:4)
      INTEGER(ap) :: ex,ey

! Loop on elements:
      ex=0
      ey=1
      DO 200 iE=1,nElm  ! loop over all cubes
         Cxx=ConAxx(iE)
         Cyy=ConAyy(iE)
         Czz=ConAzz(iE)
         Cxy=ConAxy(iE)
         Cxz=ConAxz(iE)
         Cyz=ConAyz(iE)
         ex=ex+1							!Count the position of the element in the grid (Couvreur dec 2009)
         IF (ex.GT.nex) THEN
            ex=1
            ey=ey+1
            IF (ey.GT.ney) ey=1
         ENDIF
         ! Loop on subelements: 5 tetrahedrals per soil cube
         DO 120 iSE=1,5     
            i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
            j=elmnod(iL(2,iSE,subN(iE)),iE)
            k=elmnod(iL(3,iSE,subN(iE)),iE)
            l=elmnod(iL(4,iSE,subN(iE)),iE)
            !if master node then coefficient calculated as normal; secondly signs are needed for orientation
            IF (((ex.EQ.nex).OR.(ey.EQ.ney)).AND.continu) THEN					!"Bridge sub-element?" (Couvreur dec 2009)
               xGridi=xGrid(i)
               xGridj=xGrid(j)
               xGridk=xGrid(k)
               xGridl=xGrid(l)
               yGridi=yGrid(i)
               yGridj=yGrid(j)
               yGridk=yGrid(k)
               yGridl=yGrid(l)
               IF (ex.EQ.nex) THEN
                  SELECT CASE(subN(iE))
                  CASE(1) 
                     SELECT CASE(iSE)
                     CASE (1,3) 
                        xGridj = xGridj + nex*dxgrid
                     CASE (2,4)
                        xGridj = xGridj + nex*dxgrid
                        xGridk = xGridk + nex*dxgrid
                        xGridl = xGridl + nex*dxgrid                           
                     CASE (5)
                        xGridj = xGridj + nex*dxgrid
                        xGridl = xGridl + nex*dxgrid
                     END SELECT
                  CASE (2)
                     SELECT CASE(iSE)
                     CASE (1,3) 
                        xGridj = xGridj + nex*dxgrid
                     CASE (2) 
                        xGridj = xGridj + nex*dxgrid
                        xGridk = xGridk + nex*dxgrid
                        xGridl = xGridl + nex*dxgrid                           
                     CASE (4)
                        xGridk = xGridk + nex*dxgrid
                        xGridj = xGridj + nex*dxgrid  
                     CASE (5) 
                        xGridj = xGridj + nex*dxgrid
                        xGridi = xGridi + nex*dxgrid
                        xGridl = xGridl + nex*dxgrid  
                     END SELECT
                  END SELECT
               END IF
               
               IF (ey.EQ.ney) THEN
                  SELECT CASE (subN(iE))
                  CASE (1) 
                     SELECT CASE (iSE)
                     CASE (1) 
                        yGridj = yGridj + ney*dygrid
                        yGridk = yGridk + ney*dygrid
                        yGridl = yGridl + ney*dygrid
                     CASE (2,3) 
                        yGridk = yGridk + ney*dygrid
                     CASE (4) 
                        yGridi = yGridi + ney*dygrid
                        yGridk = yGridk + ney*dygrid
                        yGridl = yGridl + ney*dygrid
                     CASE (5) 
                        yGridk = yGridk + ney*dygrid
                        yGridl = yGridl + ney*dygrid
                     END SELECT
                  CASE(2) 
                     SELECT CASE(iSE)
                     CASE(1,2) 
                        yGridk = yGridk + ney*dygrid  
                     CASE (3,5)
                        yGridk = yGridk + ney*dygrid
                        yGridj = yGridj + ney*dygrid
                        yGridl = yGridl + ney*dygrid
                     CASE (4)
                        yGridk = yGridk + ney*dygrid
                        yGridl = yGridl + ney*dygrid
                     END SELECT
                  END SELECT
               END IF
               
               bi(1,iSE,iE)=-(yGridk-yGridj)*(zGrid(l)-zGrid(j))+(yGridl-yGridj)*(zGrid(k)-zGrid(j))			!bi, ci and di are now globals 
               bi(2,iSE,iE)=+(yGridl-yGridk)*(zGrid(i)-zGrid(k))-(yGridi-yGridk)*(zGrid(l)-zGrid(k))
               bi(3,iSE,iE)=-(yGridi-yGridl)*(zGrid(j)-zGrid(l))+(yGridj-yGridl)*(zGrid(i)-zGrid(l))
               bi(4,iSE,iE)=+(yGridj-yGridi)*(zGrid(k)-zGrid(i))-(yGridk-yGridi)*(zGrid(j)-zGrid(i))
               ci(1,iSE,iE)=+(xGridk-xGridj)*(zGrid(l)-zGrid(j))-(xGridl-xGridj)*(zGrid(k)-zGrid(j))
               ci(2,iSE,iE)=-(xGridl-xGridk)*(zGrid(i)-zGrid(k))+(xGridi-xGridk)*(zGrid(l)-zGrid(k))
               ci(3,iSE,iE)=+(xGridi-xGridl)*(zGrid(j)-zGrid(l))-(xGridj-xGridl)*(zGrid(i)-zGrid(l))
               ci(4,iSE,iE)=-(xGridj-xGridi)*(zGrid(k)-zGrid(i))+(xGridk-xGridi)*(zGrid(j)-zGrid(i))
               di(1,iSE,iE)=-(xGridk-xGridj)*(yGridl-yGridj)+(xGridl-xGridj)*(yGridk-yGridj)
               di(2,iSE,iE)=+(xGridl-xGridk)*(yGridi-yGridk)-(xGridi-xGridk)*(yGridl-yGridk)
               di(3,iSE,iE)=-(xGridi-xGridl)*(yGridj-yGridl)+(xGridj-xGridl)*(yGridi-yGridl)
               di(4,iSE,iE)=+(xGridj-xGridi)*(yGridk-yGridi)-(xGridk-xGridi)*(yGridj-yGridi)
               ! coefficient a-> shape_function = a + bx +cy +dz
               ai(1)=xGridj*yGridk*zGrid(l) + xGridk*yGridl*zGrid(j) + xGridl*yGridj*zGrid(k) - &
                    xGridl*yGridk*zGrid(j) - xGridj*yGridl*zGrid(k) - xGridk*yGridj*zGrid(l)
               ai(2)=xGridi*yGridl*zGrid(k) + xGridk*yGridi*zGrid(l) + xGridl*yGridk*zGrid(i) - &
                    xGridi*yGridk*zGrid(l) - xGridk*yGridl*zGrid(i) - xGridl*yGridi*zGrid(k)
               ai(3)=xGridi*yGridj*zGrid(l) + xGridj*yGridl*zGrid(i) + xGridl*yGridi*zGrid(j) - &
                    xGridi*yGridl*zGrid(j) - xGridj*yGridi*zGrid(l) - xGridl*yGridj*zGrid(i)
               ai(4)=- xGridi*yGridj*zGrid(k) - xGridj*yGridk*zGrid(i) - xGridk*yGridi*zGrid(j) + &
                    xGridi*yGridk*zGrid(j) + xGridj*yGridi*zGrid(k) + xGridk*yGridj*zGrid(i)
               
               Deter(iSE,iE)=(xGridl-xGridi)*bi(4,iSE,iE)+(yGridl-yGridi)*ci(4,iSE,iE)+(zGrid(l)-zGrid(i))*di(4,iSE,iE)	!Deter, Ax, Ay, Az and B1fact are now globals (Couvreur mar 2010)
               Ax(:,iSE,iE)=Cxx*bi(:,iSE,iE)+Cxy*ci(:,iSE,iE)+Cxz*di(:,iSE,iE)
               Ay(:,iSE,iE)=Cxy*bi(:,iSE,iE)+Cyy*ci(:,iSE,iE)+Cyz*di(:,iSE,iE)
               Az(:,iSE,iE)=Cxz*bi(:,iSE,iE)+Cyz*ci(:,iSE,iE)+Czz*di(:,iSE,iE)
            ELSE
               bi(1,iSE,iE)=-(yGrid(k)-yGrid(j))*(zGrid(l)-zGrid(j))+(yGrid(l)-yGrid(j))*(zGrid(k)-zGrid(j))
               bi(2,iSE,iE)=+(yGrid(l)-yGrid(k))*(zGrid(i)-zGrid(k))-(yGrid(i)-yGrid(k))*(zGrid(l)-zGrid(k))
               bi(3,iSE,iE)=-(yGrid(i)-yGrid(l))*(zGrid(j)-zGrid(l))+(yGrid(j)-yGrid(l))*(zGrid(i)-zGrid(l))
               bi(4,iSE,iE)=+(yGrid(j)-yGrid(i))*(zGrid(k)-zGrid(i))-(yGrid(k)-yGrid(i))*(zGrid(j)-zGrid(i))
               ci(1,iSE,iE)=+(xGrid(k)-xGrid(j))*(zGrid(l)-zGrid(j))-(xGrid(l)-xGrid(j))*(zGrid(k)-zGrid(j))
               ci(2,iSE,iE)=-(xGrid(l)-xGrid(k))*(zGrid(i)-zGrid(k))+(xGrid(i)-xGrid(k))*(zGrid(l)-zGrid(k))
               ci(3,iSE,iE)=+(xGrid(i)-xGrid(l))*(zGrid(j)-zGrid(l))-(xGrid(j)-xGrid(l))*(zGrid(i)-zGrid(l))
               ci(4,iSE,iE)=-(xGrid(j)-xGrid(i))*(zGrid(k)-zGrid(i))+(xGrid(k)-xGrid(i))*(zGrid(j)-zGrid(i))
               di(1,iSE,iE)=-(xGrid(k)-xGrid(j))*(yGrid(l)-yGrid(j))+(xGrid(l)-xGrid(j))*(yGrid(k)-yGrid(j))
               di(2,iSE,iE)=+(xGrid(l)-xGrid(k))*(yGrid(i)-yGrid(k))-(xGrid(i)-xGrid(k))*(yGrid(l)-yGrid(k))
               di(3,iSE,iE)=-(xGrid(i)-xGrid(l))*(yGrid(j)-yGrid(l))+(xGrid(j)-xGrid(l))*(yGrid(i)-yGrid(l))
               di(4,iSE,iE)=+(xGrid(j)-xGrid(i))*(yGrid(k)-yGrid(i))-(xGrid(k)-xGrid(i))*(yGrid(j)-yGrid(i))
               ! coefficient a-> shape_function = a + bx +cy +dz
               ai(1)=xGrid(j)*yGrid(k)*zGrid(l) + xGrid(k)*yGrid(l)*zGrid(j) + xGrid(l)*yGrid(j)*zGrid(k) - &
                    xGrid(l)*yGrid(k)*zGrid(j) - xGrid(j)*yGrid(l)*zGrid(k) - xGrid(k)*yGrid(j)*zGrid(l)
               ai(2)=xGrid(i)*yGrid(l)*zGrid(k) + xGrid(k)*yGrid(i)*zGrid(l) + xGrid(l)*yGrid(k)*zGrid(i) - &
                    xGrid(i)*yGrid(k)*zGrid(l) - xGrid(k)*yGrid(l)*zGrid(i) - xGrid(l)*yGrid(i)*zGrid(k)
               ai(3)=xGrid(i)*yGrid(j)*zGrid(l) + xGrid(j)*yGrid(l)*zGrid(i) + xGrid(l)*yGrid(i)*zGrid(j) - &
                    xGrid(i)*yGrid(l)*zGrid(j) - xGrid(j)*yGrid(i)*zGrid(l) - xGrid(l)*yGrid(j)*zGrid(i)
               ai(4)=- xGrid(i)*yGrid(j)*zGrid(k) - xGrid(j)*yGrid(k)*zGrid(i) - xGrid(k)*yGrid(i)*zGrid(j) + &
                    xGrid(i)*yGrid(k)*zGrid(j) + xGrid(j)*yGrid(i)*zGrid(k) + xGrid(k)*yGrid(j)*zGrid(i)

               Deter(iSE,iE)=(xGrid(l)-xGrid(i))*bi(4,iSE,iE)+(yGrid(l)-yGrid(i))*ci(4,iSE,iE)+(zGrid(l)-zGrid(i))*di(4,iSE,iE)
               Ax(:,iSE,iE)=Cxx*bi(:,iSE,iE)+Cxy*ci(:,iSE,iE)+Cxz*di(:,iSE,iE)
               Ay(:,iSE,iE)=Cxy*bi(:,iSE,iE)+Cyy*ci(:,iSE,iE)+Cyz*di(:,iSE,iE)
               Az(:,iSE,iE)=Cxz*bi(:,iSE,iE)+Cyz*ci(:,iSE,iE)+Czz*di(:,iSE,iE)
            ENDIF
            VE=ABS(Deter(iSE,iE))/6.
            IF (Deter(iSE,iE).EQ.0 .OR. VE.EQ.0) THEN
               WRITE(*,*)'Ve=0',VE,'deter=',Deter(iSE,iE),'element nr.',iE
               PRINT*,'coo i',xgrid(i),ygrid(i),zgrid(i)
               PRINT*,'coo j',xgrid(j),ygrid(j),zgrid(j)
               PRINT*,'coo k',xgrid(k),ygrid(k),zgrid(k)
               PRINT*,'coo l',xgrid(l),ygrid(l),zgrid(l)
               READ(*,*)
            ENDIF

            DO 110 ii=1,4
               B1fact(ii,iSE,iE)=Cxz*bi(ii,iSE,iE)+Cyz*ci(ii,iSE,iE)+Czz*di(ii,iSE,iE)
               DO 100 jj=1,4
                  E(jj,ii,iSE,iE) = Cxx*bi(ii,iSE,iE)*bi(jj,iSE,iE)+Cyy*ci(ii,iSE,iE)*ci(jj,iSE,iE)+Czz*di(ii,iSE,iE)*di(jj,iSE,iE) &
                  +Cxy*(bi(ii,iSE,iE)*ci(jj,iSE,iE)+bi(jj,iSE,iE)*ci(ii,iSE,iE))+Cxz*(bi(ii,iSE,iE)*di(jj,iSE,iE)+bi(jj,iSE,iE)*di(ii,iSE,iE)) &
                  +Cyz*(ci(ii,iSE,iE)*di(jj,iSE,iE)+ci(jj,iSE,iE)*di(ii,iSE,iE))
               100  CONTINUE
            110  CONTINUE
         120  CONTINUE
      200  CONTINUE
      RETURN
      END SUBROUTINE CalcGeom
!****************************************************************************
!> assambles the matrix and the right-hand side of the linear equation system
!> based on the standard galerkin finite element method
    SUBROUTINE Reset(dt,IAD,IADN,IADD)
      USE ParamData, ONLY: maxbnd
      USE GridData, ONLY: nPt,nElm,rootSkOld,subN,elmnod,iL,deter,sink,sink_cube,b1fact, iadd_temp,e,width
      USE PlntData, ONLY: Tpot,Tact
      USE CumData
      USE SolData
      USE DomData
      USE MatData
      USE BoundData,ONLY: iBCPt
      USE RootData,ONLY: lDou,lCou
      USE tmctrl, ONLY: tlevel_soil
      USE RhizoData, ONLY: Rnorm, tauTht, hEqRhizo, RhizoModel, lRhizo
      IMPLICIT NONE

      INTEGER(ap) :: ind,ind1
      INTEGER(ap) :: iSE,i,j,k,l,iE,n,ii,jj,kk,iG,jG,kkk!, kk2, NiG, DDiG
      REAL(sp) :: cone,cape

      REAL(sp) :: fmul,bmul
      REAL(sp) :: amul,QN
      REAL(sp) :: dt
      REAL(sp) :: SinkE
      REAL(sp) :: VE
      INTEGER(ap) :: IAD(maxbnd,nPt),IADN(nPt),IADD(nPt)
      REAL(sp), ALLOCATABLE,DIMENSION (:) :: B1,F,DS
      ALLOCATE (B1(nPt),F(nPt),DS(nPt))
	  !> \param dt time step size
	  !> \param IAD adjacency matrix (nodal connections)
	  !> \param IADN number of adjacent nodes in IAD (self-inclusive)
	  !> \param IADD position of diagonal in adjacency matrix
      a_sparse(1:irow(npt+1)-1) = 0.0d0
      B(1:npt)=0.0
      B1=0.0
      DS=0.0
      F=0.0
      RootSkold=0.0

      !save matrix index for sparse format in the first step -> tlevel_soil = false
      k= 0
      IF(.NOT.tlevel_soil) THEN
         DO  iE=1,nElm  ! cubes
            DO  iSE=1,5! Loop on subelements: 5 tetrahedrals per soil cube
               DO  ii=1,4
                  iG=elmnod(iL(ii,iSE,subN(iE)),iE)
                  !NiG = IADN(iG)
                  !DDiG = IADD(iG)
                  DO  jj=1,4
                     jG=elmnod(iL(jj,iSE,subN(iE)),iE)
                     !CALL Find2(iG,jG,kk,nPt,maxbnd,IAD,NiG,DDiG)
                     CALL FIND (iG,jG,kk,nPt,maxbnd,IAD,IADN)
                     k = k +1
                     iadd_temp(k) = kk
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         tlevel_soil = .TRUE.
      ENDIF
     
      ! Loop on elements:
      kkk = 0 !index for matrix index: iadd_temp
      DO 200 iE=1,nElm  ! cubes
         ! Loop on subelements: 5 tetrahedrals per soil cube
         DO 120 iSE=1,5
            i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tetrahedrals which is equal to a cubic.
            j=elmnod(iL(2,iSE,subN(iE)),iE)
            k=elmnod(iL(3,iSE,subN(iE)),iE)
            l=elmnod(iL(4,iSE,subN(iE)),iE)
            CapE=(Cap(i)+Cap(j)+Cap(k)+Cap(l))/4.
            ConE=(Con(i)+Con(j)+Con(k)+Con(l))/4.
            VE=ABS(Deter(iSE,iE))/6.			!Deter and B1fact are now globals 
            AMul=ConE/VE/36.		![/L²/T]
            BMul=ConE/6.			![L/T]
            FMul=VE/20.			![L³]
            SinkE=FMul*5*sink_cube(iE)	![L³/T]
            DS(i)=DS(i)+SinkE		![L³/T]
            DS(j)=DS(j)+SinkE
            DS(k)=DS(k)+SinkE
            DS(l)=DS(l)+SinkE
            !SinkE=(sink(i)+sink(j)+sink(k)+sink(l))/4	!Previous way to proceed ("sink node")
            !DS(i)=DS(i)+FMul*(4*SinkE+sink(i))
            !DS(j)=DS(j)+FMul*(4*SinkE+sink(j))
            !DS(k)=DS(k)+FMul*(4*SinkE+sink(k))
            !DS(l)=DS(l)+FMul*(4*SinkE+sink(l))
            F(i)=F(i)+FMul*(4*CapE+Cap(i)) !*5	![L³/P] or [L²]
            F(j)=F(j)+FMul*(4*CapE+Cap(j)) !*5
            F(k)=F(k)+FMul*(4*CapE+Cap(k)) !*5
            F(l)=F(l)+FMul*(4*CapE+Cap(l)) !*5
            DO 110 ii=1,4
               !NiG = IADN(iG)
               !DDiG = IADD(iG)
               iG=elmnod(iL(ii,iSE,subN(iE)),iE)
               !F(iG)=F(iG)+FMul*(4*CapE+Cap(iG)) !*5
               B1(iG)=B1(iG)+BMul*B1fact(ii,iSE,iE)	![L³/T]	(B1fact [L²])
               DO 100 jj=1,4
                  jG=elmnod(iL(jj,iSE,subN(iE)),iE)
                  kkk = kkk+1
                  !call FIND (iG,jG,kk,nPt,maxbnd,IAD,IADN)
                  !normal value added to array A_sparse (CSR format)
                  A_sparse(IROW(iG)-1+iadd_temp(kkk))=A_sparse(IROW(iG)-1+iadd_temp(kkk))+AMul*E(jj,ii,iSE,iE)              
                  !A_sparse(IROW(iG)-1+kk)=A_sparse(IROW(iG)-1+kk)+AMul*E(jj,ii,iSE,iE)                  
               100  CONTINUE
            110  CONTINUE
         120  CONTINUE
      200  CONTINUE

      DO 240 N=1,nPt
         IF (Kode(N).GT.0) THEN !hbc
            QN=B1(N)+DS(N)+F(N)*(theta(N)-theta_old(N))/dt ![L³/T]	Tom had deleted the term +F(N) etc... same than SWMS
            ind=IROW(N)
            ind1=IROW(N+1)-1
            DO j=ind,ind1
               QN=QN+A_sparse(j)*hNew(JCOL(j))
               !if a row is part of a bound.cond. then get the values in A_sparse and multiply by the corresponding column values in hNew
            ENDDO
            Q(N)=QN
         ELSEIF (Kode(N).EQ.0) THEN !no BC defined
            Q(N)=0
!         elseif (Kode(N).EQ.-1) then !rain (homogeneous or not)
         ELSEIF (Kode(N).EQ.-2) THEN !free drainage
            Q(N)=-con(N)*Width(N)
         ENDIF
         ! Construction of effective matrix:
         ! use row index - 1 to get the number of previous entries, add up the entry number of the diagonal value to find the location of the diagonal value in A_sparse
         
         IF(.NOT. lRhizo .OR. RhizoModel .EQ. 1 .OR. RhizoModel .EQ. 2) THEN
             ! A matrix when non equilibrium is not considered.
             A_sparse(IROW(N)-1+IADD(N)) = A_sparse(IROW(N)-1+IADD(N)) +  F(N)/dt !CSR format
         ELSE
            ! A Matrix when non equilibrium is considered.
             A_sparse(IROW(N)-1+IADD(N)) = A_sparse(IROW(N)-1+IADD(N)) + Rnorm(N)*F(N)/dt + (1. - Rnorm(N))*FMul/tauTht(N) !CSR format
         ENDIF
      240 CONTINUE
         ! Complete construction of RHS vector:

         IF(.NOT. lRhizo .OR. RhizoModel .EQ. 1 .OR. RhizoModel .EQ. 2) THEN
            B=Q-B1+hold*F/dt-DS
        ELSE
            B = Q - B1 + Rnorm*hold*F/dt-DS + (1. - Rnorm)*hEqRhizo*FMul/tauTht
        ENDIF
         
         DEALLOCATE (B1,F,DS)
         RETURN
       END SUBROUTINE Reset
!****************************************************************************
!> assambles the matrix and the right-hand side of the linear equation system
!> based on the standard galerkin finite element method
!> used the mass lumping algorith from Celia et al 1990 - not clear if this still works
    SUBROUTINE ResetCelia(dt,IAD,IADN,IADD)
      USE ParamData, ONLY: maxbnd
      USE GridData, ONLY: nPt,nElm,rootSkOld,subN,elmnod,iL,deter,sink,sink_cube,b1fact, iadd_temp,e,width
      USE PlntData, ONLY: Tpot,Tact
      USE CumData
      USE SolData
      USE DomData
      USE MatData
      USE WatFun
      USE BoundData,ONLY: iBCPt
      USE tmctrl, ONLY: tlevel_soil
      IMPLICIT NONE

      INTEGER(ap) :: ind,ind1
      INTEGER(ap) :: iSE,i,j,k,l,iE,n,ii,jj,kk,iG,jG,kkk
      REAL(sp) :: cone,cape
      REAL(sp) :: fmul,bmul
      REAL(sp) :: amul,QN,VE
      REAL(sp) :: dt
      REAL(sp) :: SinkE
      !REAL(sp) :: vol_elem_sink
      INTEGER(ap) :: IAD(maxbnd,nPt),IADN(nPt),IADD(nPt)
      REAL(sp), ALLOCATABLE,DIMENSION (:) ::B1, F ,DS
      ALLOCATE (B1(nPt),F(nPt),DS(nPt))
      !> \param dt time step size
	  !> \param IAD adjacency matrix (nodal connections)
	  !> \param IADN number of adjacent nodes in IAD (self-inclusive)
	  !> \param IADD position of diagonal in adjacency matrix
      a_sparse(1:irow(npt+1)-1) = 0.0d0
      B(1:npt)=0.0
      B1=0.0
      DS=0.0
      F=0.0
      RootSkold=0.0

      !save matrix index for sparse format in the first step -> tlevel_soil = false
      k= 0
      IF(.NOT.tlevel_soil) THEN
         DO  iE=1,nElm  ! cubes
            DO  iSE=1,5! Loop on subelements: 5 tetrahedrals per soil cube
               DO  ii=1,4
                  iG=elmnod(iL(ii,iSE,subN(iE)),iE)
                  !NiG = IADN(iG)
                  !DDiG = IADD(iG)
                  DO  jj=1,4
                     jG=elmnod(iL(jj,iSE,subN(iE)),iE)
                     !CALL Find2(iG,jG,kk,nPt,maxbnd,IAD,NiG,DDiG)
                     CALL FIND (iG,jG,kk,nPt,maxbnd,IAD,IADN)
                     k = k +1
                     iadd_temp(k) = kk
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         tlevel_soil = .TRUE.
      ENDIF
      
      kkk=0
      ! Loop on elements:
      DO 200 iE=1,nElm  !cubes
         ! Loop on subelements: 5 tetrahedrals per half soil cube
         DO 120 iSE=1,5
            i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
            j=elmnod(iL(2,iSE,subN(iE)),iE)
            k=elmnod(iL(3,iSE,subN(iE)),iE)
            l=elmnod(iL(4,iSE,subN(iE)),iE)
            CapE=(Cap(i)+Cap(j)+Cap(k)+Cap(l))/4
            ConE=(Con(i)+Con(j)+Con(k)+Con(l))/4
            VE=ABS(Deter(iSE,iE))/6.			!Deter and B1fact are now globals (Couvreur mar 2010)
            AMul=ConE/VE/36.
            BMul=ConE/6.
            FMul=VE/20.
            SinkE=(sink(i)+sink(j)+sink(k)+sink(l))/4
            !DS(i)=DS(i)+FMul*(4*SinkE+sink(i))
            !DS(j)=DS(j)+FMul*(4*SinkE+sink(j))
            !DS(k)=DS(k)+FMul*(4*SinkE+sink(k))
            !DS(l)=DS(l)+FMul*(4*SinkE+sink(l))
            DS(i)=DS(i)+FMul*5*sink_cube(iE)
            DS(j)=DS(j)+FMul*5*sink_cube(iE)
            DS(k)=DS(k)+FMul*5*sink_cube(iE)
            DS(l)=DS(l)+FMul*5*sink_cube(iE)
!flow leaving soil is average sink times their volume  over all the elements
            RootSkold=RootSkold+VE*SinkE
            !RootSkold=RootSkold+VE*sink_cube(iE)
            IF (Tpot.NE.0) THEN!
                  !              write(*,*)'Rootsk=',RootSk,'SE',SinkE
                   !             write(*,*)'sink per node=',sink(k),sink(j),sink(l),sink(i)
              ! vol_elem_sink = vol_elem_sink + VE
            ENDIF
            DO 110 ii=1,4
               !NiG = IADN(iG)
               !DDiG = IADD(iG)
               iG=elmnod(iL(ii,iSE,subN(iE)),iE)
               F(iG)=F(iG)+FMul*5!(4*CapE+Cap(iG)) !*5
               B1(iG)=B1(iG)+BMul*B1fact(ii,iSE,iE)
               DO 100 jj=1,4
                  jG=elmnod(iL(jj,iSE,subN(iE)),iE)
                  kkk = kkk+1
                  !call FIND (iG,jG,kk,nPt,maxbnd,IAD,IADN)
                  !normal value added to array A_sparse (CSR format)
                  A_sparse(IROW(iG)-1+iadd_temp(kkk))=A_sparse(IROW(iG)-1+iadd_temp(kkk))+AMul*E(jj,ii,iSE,iE)              
                  !A_sparse(IROW(iG)-1+kk)=A_sparse(IROW(iG)-1+kk)+AMul*E(jj,ii,iSE,iE)                  
               100  CONTINUE
            110  CONTINUE
          120  CONTINUE
      200  CONTINUE
       
      DO 240 N=1,nPt         
         IF (Kode(N).GT.0) THEN !hbc
            QN=B1(N)+DS(N)+F(N)*(theta(N)-theta_old(N))/dt !Tom had deleted the term +F(N) etc... same than SWMS
            ind=IROW(N)
            ind1=IROW(N+1)-1
            DO j=ind,ind1
               QN=QN+A_sparse(j)*hNew(JCOL(j))
               !if a row is part of a bound.cond. then get the values in A_sparse and multiply by the corresponding column values in hNew
            ENDDO
            Q(N)=QN
         ELSEIF (Kode(N).EQ.0) THEN !no BC defined
            Q(N)=0
!         elseif (Kode(N).EQ.-1) then !rain (homogeneous or not)
         ELSEIF (Kode(N).EQ.-2) THEN !free drainage
            Q(N)=-con(N)*Width(N)
            !print *,  'width: ', Width(N)
         ENDIF
      240 CONTINUE
      
! Complete construction of RHS vector and form effective matrix:
! use row index - 1 to get the number of previous entries, add up the entry number of the diagonal value to find the location of the diagonal value in A_sparse
      DO 310 i=1,nPt
         B(i)=F(i)*Cap(i)*hNew(i)/dt-F(i)*(theta(i)-theta_old(i))/dt+ Q(i)-B1(i)-DS(i) ! =Q(i)-B1(i)+hold(i)*F(i)/dt-DS(i)
         A_sparse(IROW(i)-1+IADD(i)) = A_sparse(IROW(i)-1+IADD(i)) + F(i)*Cap(i)/dt  !+ F(i)/dt !CSR format
      310 CONTINUE
         
      DEALLOCATE (B1,F,DS)
      RETURN
   END SUBROUTINE ResetCelia
!***********************************************************************
!> adding Dirichlet boundary condition to the sparse matrix 
 SUBROUTINE Dirich(IADD)
     USE GridData, ONLY:
     USE MatData
     USE SolData, ONLY: Kode, hnew
     IMPLICIT NONE
 
     INTEGER(ap) :: n,IADD(nPt),ind
	 !> \param IADD position of diagonal in adjacency matrix
!usage of compressed sparse row format
! IA(1:nPt+1) represent the rows + 1. The values in IA correspond to the positions in JA
! (column indices) and AA (nonzero values). The indicators ind and ind1 can be defined
! which show for one row the indices for the arrays JA and AA.
!
! example (from Youcef Saad: numerical methods for large eigenvalue problems 1992)
!
! AA = 1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12.
! JA = 1 4 1 2 4 1 3 4 5 3 4 5
! IA = 1 3 6 10 12 13
!
! A = |   1. 0. 0.  2.  0.   |
!     |   3. 4. 0.  5.  0.   |
!     |   6. 0. 7.  8.  9.   |
!     |   0. 0. 10. 11. 0.   |
!     |   0. 0. 0.  0.  12.  |
     DO 70 N=1,nPt
        IF (KODE(N).LT.1) GOTO 70
        !ind is row index; to find where diagonal position in A_sparse is
        !located take row index -1 (gives number of previous entries in
        !A_sparse and add up the number which corresponds to the diagonal value
        !, i.e., IADD(N)
        ind=IROW(N)
        A_sparse(ind-1+IADD(N))=10.d30
        B(N)=10.d30*hNew(N)
70      CONTINUE
        RETURN
END 
!*****************************************************************************
!> time step control: calculates the next time step size for  
!> water, solute and root growth
      SUBROUTINE TmCont(iter,t,dt,dtOpt,tCallR,tFEMRoo,dtMaxC,tcBCr,tProf,tProbe)
      USE typedef
      USE tmctrl, ONlY: dtMax,dtMin,facdec,facinc,tmax
      USE Rootdata, ONLY :lno_RWU,lno_root_growth
      IMPLICIT NONE
      REAL(sp), INTENT(out) :: dt
      REAL(sp), INTENT(in) :: t,dtMaxC,tCallR,tFEMRoo,tcBCr,tProf,tProbe
      REAL(sp), INTENT(inout) :: dtOpt
      REAL(sp)::tfix
      INTEGER(ap) :: iter
	  !> \param iter
	  !> \param t current simulation time
 	  !> \param dt current time step size
  	  !> \param dtOpt optimal time step size calculated in this subroutine
  	  !> \param tCallR
	  !> \param tFEMRoo output time for the next FEM/SOIL and ROOT output
 	  !> \param dtMaxC
  	  !> \param tcBCr
  	  !> \param tProf
	  !> \param tProbe
      dtMax=MIN(dtMax,dtMaxC)
      tFix=MIN(tcallr,tFEMRoo,tMax,tcBCr,tProf,tProbe)
      IF (iter.LE.3.AND.(tFix-t).GE.FacInc*dtOpt)  dtOpt=MIN(dtMax,FacInc*dtOpt)!3
      IF (iter.GE.7)  dtOpt=MAX(dtMin,FacDec*dtOpt) !7
      dt=MIN(dtOpt,tFix-t)
      dt=MIN((tFix-t)/ANINT((tFix-t)/dt),dtMax)
      IF(tFix-t.NE.dt.AND.dt.GT.(tFix-t)/2._dp) dt=(tFix-t)/2._dp
      PRINT *,iter,'dt=',dt

      RETURN
      END SUBROUTINE TmCont
!****************************************************************
      SUBROUTINE WatInf(dt)
      USE GridData, ONLY: dxGrid,dyGrid,dzGrid,nElm,nPt,RootSk,sink,sink_cube,Vn
      USE CumData
      USE RootData, ONLY: lFed, lCou
      USE SolData, ONLY: Kode
      IMPLICIT NONE
      REAL(sp):: vMean(3),dt
      INTEGER(ap)::i,j
      vMean=0.0_sp

      DO i=1,nPt
         j=iabs(Kode(i))
         wCumA=wCumA+ABS(Q(i))*dt
         IF (j.NE.0) THEN
            vMean(j)=vMean(j)-Q(i)
         ENDIF
      END DO
      IF (lCou.OR.lFed) THEN
         RootSk=0.0_dp
         DO i=1,nElm
            RootSk=RootSk+sink_cube(i)*dxGrid*dyGrid*dzGrid
         ENDDO
      ENDIF
      wCumA=wCumA+(RootSk*dt)!abs(RootSk*dt)
      CumRt=CumRt+RootSk*dt
      VolSink = RootSk*dt
      wCumT=CumRt
      CumQ(1)=CumQ(1)+(vMean(1)+vMean(3))*dt
      CumQ(2)=CumQ(2)+vMean(2)*dt
      wCumT=wCumT+CumQ(1)+CumQ(2)
      RETURN
    END SUBROUTINE WatInf
!****************************************************************
!< SUBROUTINE TO CALCULATE THE INITIAL WATER CONTENT IN THE RHIZOSPEHRE
  SUBROUTINE IniRhizoTheta(h0,nodeID)
     USE Typedef
     USE RhizoData, ONLY: bulkPara, StaticRhizoPara, thetaTot, Rnorm
     USE solData,   ONLY: theta
     USE GridData,  ONLY: nPt
     IMPLICIT NONE
     REAL(sp) :: h0, thetaEqIni=0.
     REAL(dp) :: lambda_s, hcr_s, thtR, thtS
     INTEGER  :: nodeID

     thtR = bulkPara(1); thtS = bulkPara(2)
     lambda_s = StaticRhizoPara(1)
     hcr_s    = StaticRhizoPara(2)
    
     IF (Rnorm(nodeID) .LT. 1.) THEN
         IF (h0 .LT. hcr_s) THEN
             thetaEqIni = thtR +(thtS-thtR)*(hcr_s/h0)**lambda_s
         ELSE
             thetaEqIni = thtS
         ENDIF
         thetaTot(nodeID) = Rnorm(nodeID)*theta(nodeID) + (1. - Rnorm(nodeID))*thetaEqIni
     ELSE
         thetaTot(nodeID) = theta(nodeID)
     ENDIF
     RETURN
  END SUBROUTINE IniRhizoTheta
!***************************************************************
!< SUBROUTINE TO CALCULATE heq (equilibrium pressure head in the rhizopshere
  SUBROUTINE hEq(thtIn, hIn, i)
    ! This subroutine calculate hEq base on theta and location (i.e. cTot)
    USE Typedef
    USE RhizoData, ONLY: bulkPara, RhizoPara, hEqRhizo, cTot_r
    IMPLICIT NONE
    REAL(dp) :: thtR, thtS, hcr, lambda, cw, rhow, rhob, beta, omega
    REAL(sp) :: thtIn, hIn, h1
    INTEGER  :: i

    thtR = bulkPara(1); thtS = bulkPara(2); lambda = bulkPara(3); hcr = bulkPara(4)
    omega = RhizoPara(1); beta = RhizoPara(2); rhob = RhizoPara(8); rhow = RhizoPara(9)
    
    IF (thtIn .LT. thtS) THEN
        h1 = hcr/((thtIn-thtR)/(thtS-thtR))**(1./lambda)
    ELSE
        h1 = hIn
    ENDIF
    cw = cTot_r(i)*rhob/(rhow*thtIn)
    hEqRhizo(i) = h1 - omega*cw**beta
    RETURN
  END SUBROUTINE hEq
!**************************************************************
!< SUBROUTINE TO CALCULATE thetaNonEquilibrium
  SUBROUTINE thtNonEq(hNew, dt, i)
      USE Typedef
      USE SolData, ONLY: theta, theta_old
      USE RhizoData, ONLY: bulkPara, RhizoPara, thetaNonEq, thetaNonEqOld, thetaTot, hEqRhizo, thetaTot, Rnorm, tauTht
      IMPLICIT NONE
      REAL(sp) :: hNew, dt
      REAL(sp) :: tolerance, diffTht
      REAL(dp) :: gamm, tau0, thtS
      INTEGER  :: i, maxIter, iter
      ! Initiate parameters
      tau0 = RhizoPara(6)
      gamm = RhizoPara(7)
      thtS = bulkPara(2)
      maxIter = 500
      tolerance = 1e-4
      iter = 0
      diffTht = 1.0
      thetaNonEq(i) = theta(i)*0.5
      ! The main loop that calculate nonequilibrium parameters: 
      ! At each iteration the the equilibrium presure head at the rhizosphere 
      ! is calculate. Following, tau, the mucilage water content and 
      ! the total water content  are calculated. 
      ! The loop is ended when difference in the local water content between two successive
      ! iteration is lower than the tolerance or the iteration exceed the maximum
      ! allowablw iteration. 
      DO WHILE (ABS(diffTht) .GT. tolerance .AND. iter .LT. maxIter)
        CALL hEq(thetaNonEq(i), hnew, i) ! Calculating the new equilibrium pressure head
        diffTht = thetaTot(i)
        tauTht(i) = thetaNonEq(i)**(-gamm)*tau0
        ! d(theta)\dt = 1/tau * h-heq ! 
        thetaNonEq(i) = thetaNonEqOld(i) + (1./tauTht(i))*(hnew - hEqRhizo(i))*dt
        IF (thetaNonEq(i) .GE. thtS) THEN
            thetaNonEq(i) = thtS
            thetaTot(i) = Rnorm(i)*theta(i) + (1. - Rnorm(i))*thetaNonEq(i)
            GOTO 251
        ENDIF
        thetaTot(i) = Rnorm(i)*theta(i) + (1. - Rnorm(i))*thetaNonEq(i)
        diffTht = diffTht - thetaTot(i) ! difference between thetaTot within two iterations.
        iter = iter + 1
        IF (iter .GE. maxIter .AND. ABS(diffTht) .GT. tolerance) THEN
            WRITE(*,*) '!!!-!!!'
            WRITE(*,*) 'rhizoERR = ' , diffTht
            GOTO 251
        ENDIF
      ENDDO
    251 CONTINUE
!    WRITE(*,*) iter, diffTht, thetaTot(i), thetaNonEq(i), theta(i)
  END SUBROUTINE thtNonEq





               
      



    

