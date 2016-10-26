! Makefile .   To switch OpenMp on: pgf90 -mp -o executable file.f90
!                                   export OMP_NUM_THREADS= 2
!                                   ./executable
! To use openmp -> switch module OMP_LIB on!!! gfortran compiler can't read OMP_LIB
! for gfortran only new compiler (v3) supports openmp
!******************************************************************************
!    Mathieu Javaux- Tom Schroeder- Jan Vanderborght Fz-Juelich (Germany)	 *
!										        *
!                                                                             *
!                                                                             *
!    R-SWMS_3D :An algorithm for three-dimensional, simultaneous modeling     *
!    of root growth, transient soil water flow, solute transport, and root    *
!    water and solute uptake.						        *
!    coupled with RootTyp (Pages et al.)                                      *
!                                                                             *
!                        Version 3.2.3.           August09 		        *
!                                                                             *
!    Program Version for VAX - FORTRAN                                        *
!                                                                             *
!   based on the code of                                                      *
!   Francesca Somma & Volker Clausnitzer, University of California, Davis     *
!       (1998)                                                                *
!******************************************************************************
       MODULE ParamData
       Use typedef
       IMPLICIT NONE
       INTEGER(sp),PARAMETER :: maxnod=50000,mxpnts=100,nTab=100,maxplant=10
       INTEGER(sp),PARAMETER :: maxbdr=2500,mxBcCh=5000,maxIrrig=10
       INTEGER(sp),PARAMETER :: mxtime=1000,mxdpth=30,maxelm=90000
       INTEGER(sp),PARAMETER :: maxgrw=100000,maxrec=1000000,maxest=100
       INTEGER(sp),PARAMETER :: maxemg=10,maxord=6,maxmat=10,maxbnd=2500
       INTEGER(sp),PARAMETER :: maxobs=30
       integer(dp),parameter:: maxAMG=5000000,maxamg2=5000000 !array size for sparse A matrix and sparse row index matrix
       REAL(sp), Parameter :: pi=3.141592!65358979323
       END MODULE ParamData
!************************************************************************
       MODULE PlntData
! data needed for plant
      Use typedef
      USE ParamData
      IMPLICIT NONE
! root characteristics and Doussan input variables
      INTEGER(sp) :: nTpot,ntTpLA,nfTpLA,ncTpLA,ntLA,ntRSR,nscRSR
      INTEGER(sp) :: ntW,nsfW,nscW,ncnc,ns,nBCr,nsfRSR,typeBCr(mxBcCh)
      REAL(sp) ::tTpot(mxBcCh),Tpotc(mxBcCh),SpWgt,TotSur
      REAL(sp) :: tTpLA(mxBcCh),TpLAc(mxBcCh)
      REAL(sp) :: sfTpLA(mxBcCh),fTpLAc(mxBcCh),scTpLA(mxBcCh)
      REAL(sp) :: tLA(mxBcCh),LAc(mxBcCh),cTpLAc(mxBcCh)
      REAL(sp) :: h50,p50,p1,p2,CMm,VMax,fk,xin,h0,h1,h2,h3
      REAL(sp) :: tRSR(mxBcCh),RSRc(mxBcCh)
      REAL(sp) :: sfRSR(mxBcCh),fRSRc(mxBcCh),scRSR(mxBcCh)
      REAL(sp) :: sc(mxBcCh),rsc(mxBcCh),cncp(mxBcCh),rscnc(mxBcCh)
      REAL(sp) :: tW(mxBcCh),Wc(mxBcCh),cRSRc(mxBcCh)
      REAL(sp) :: sfW(mxBcCh),fWc(mxBcCh),scW(mxBcCh),cWc(mxBcCh)
      REAL(sp) :: tBCr(mxBcCh),BCroot(mxBcCh)!Doussan in BC
      REAL(dp) :: Tpot,Tact(maxplant)
      END MODULE PlntData
!************************************************************************
      MODULE RootData
! data needed for root
      Use typedef
      Use paramData
      IMPLICIT NONE
      INTEGER(sp) :: nbr,nrec,ngrow,naxemg,nvch(maxord),nMPLch(maxord)
      REAL(sp) :: xplant(maxplant),yplant(maxplant)
      INTEGER(sp) :: ordseg(maxrec),ibrseg(maxrec),iaxis(maxgrw),num_seg(maxgrw),br_rec(maxgrw)
      INTEGER(sp) :: nUrf(maxord)
      INTEGER :: irnd
      INTEGER(sp) :: irecpr(maxrec),irecsg(maxgrw),nnewax(maxemg)
      INTEGER(sp) :: nestbl(maxgrw),ordgrw(maxgrw),ibrgrw(maxgrw)
      REAL(dp) :: seglen(maxrec),segsur(maxrec),segmas(maxrec),segdiam(maxrec),brlgth(maxgrw)
      REAL(sp) :: xs(maxrec),ys(maxrec),zs(maxrec),timorg(maxrec)
      REAL(sp) :: xg(maxgrw),yg(maxgrw),zg(maxgrw),brlmax(maxord)
      REAL(sp) :: ovrlen(maxgrw),agevch(maxord,mxpnts)
      REAL(sp) :: timest(maxgrw,maxest),vch(maxord,mxpnts)
      REAL(sp) :: sMPLch(maxord,mxpnts),MPLch(maxord,mxpnts)
      REAL(sp) :: age(maxord,mxpnts),Urf(maxord,mxpnts)
      REAL(sp) :: tnewax(maxemg),strsen(maxord),rdmang(maxord),ran(maxemg)
      REAL(sp) :: brspac(maxord-1),brnang(maxord-1),dtbrch(maxord-1)
      LOGICAL stopgr(maxgrw),toosml(maxrec),no_root_gwth,connex(maxgrw),rrt
      END MODULE RootData
!************************************************************************
      MODULE GridData
      USE ParamData
      Use typedef
      IMPLICIT NONE
      REAL(sp), allocatable, dimension(:) :: xgrid,ygrid,zgrid,Axy,Bxy,Dxy,Exy,betaw,betac,Width
      REAL(dp), allocatable,dimension(:) :: Wn,sink,csink,sinkOld
      REAL(dp), allocatable,dimension(:,:) :: VE
      REAL(dp), allocatable,dimension(:,:,:) :: B1fact
      REAL(dp), allocatable,dimension(:,:,:,:) :: E
      INTEGER(sp), allocatable,dimension (:,:) :: elmnod
      REAL(sp) :: dxgrid,dygrid,dzgrid,epslonPH,epslonWC,epslonR,factorRelEps,epslonS
      INTEGER(sp) ::nPt,nElm,nel,nBand,nBCPts,nex,ney,nez,itMax,itMaxRoot,nx,ny,nz
      CHARACTER Headln*80,LnUnit*5,TmUnit*5,MsUnit*5,CnUnit*5
      REAL(sp) ::xCol(100),yCol(100)
      REAL(dp) :: RootSk,checkSink,RootSkOld
      logical :: RelEps,continu
      contains
        SUBROUTINE IniGrid
	    IMPLICIT NONE
	    allocate (xgrid(1:nPt))
           allocate (ygrid(1:nPt))
           allocate (zgrid(1:nPt))
           allocate (Axy(1:nPt))
           allocate (Bxy(1:nPt))
           allocate (Dxy(1:nPt))
           allocate (Exy(1:nPt))
           allocate (betaw(1:nPt))
           allocate (betac(1:nPt))
           allocate (sink(1:nPt))
           allocate (sinkOld(1:nPt))
           allocate (csink(1:nPt))
           allocate (Wn(1:nPt))
           allocate (Width(1:nPt))
           allocate (elmnod(1:nElm,1:6))
           allocate (VE(1:nElm,1:3))
           allocate (B1fact(1:nElm,1:3,1:4))
           allocate (E(1:nElm,1:3,1:4,1:4))
	 END SUBROUTINE IniGrid
      END MODULE GridData
!*************************************************************************
       MODULE DoussanMat
       Use typedef
       Use paramData
       Use GridData
       USE RootData
       Use SparseMatrix !multiple roots
       IMPLICIT NONE

! Doussan matrices
       REAL(dp), allocatable, dimension (:,:) :: PHr_sub,PH_root_sub,axialRootFlow
       REAL(sp), allocatable, dimension (:,:) :: k_ave,B,Qd,Qi,Q_bc
       REAL(dp), allocatable, dimension (:) :: Lr,Khr,GH,qroot,l_seg
       REAL(dp), allocatable, dimension (:) :: Phi_mat,h_mat2,Joutr
       REAL(dp), allocatable, dimension (:) :: delta2,delta2old,sinkRtemp
       REAL(dp), allocatable, dimension (:) :: tempSinkR,curr_BCr,BCr_usr
       REAL(dp), allocatable, dimension (:,:) :: PHs,PHr,SinkROld,PHrOld,PHrTemp,PHsTemp,sinkR
       REAL(dp) ::KhRoot(1:3,mxBcCh),LrRoot(1:3,mxBcCh),Jintot,sinktot
       REAL(sp) :: ageKh(1:3,mxBcCh),ageLr(1:3,mxBcCh),hx_min,stresfun,stresval1,stresval2
       REAL(sp), allocatable, dimension (:,:,:) :: PH_micro2
       REAL(sp), allocatable, dimension (:,:,:,:) :: w_dis,cent,cp_mean,Intc
       REAL(dp), allocatable, dimension (:,:,:) :: w_sub,l_sub,sum_dis
       INTEGER(sp) :: nBCn,nKh(1:3),nLr(1:3),nmax,nrecOld,isubmax=10,nplant
       INTEGER(dp) :: nLibr=100e3,indexValue=1000,switchcriterion=1,n=50
       INTEGER(dp), allocatable, dimension(:) :: counter2,BCtp_usr,curr_BCtp
       INTEGER(sp), allocatable, dimension(:) :: nBC_irecn,nBC_iprvn,no_voxels !ija,sa,
       INTEGER(sp), allocatable, dimension(:,:) :: numNodes_voxel,voxel_no,nsub
       INTEGER(sp), allocatable, dimension (:,:,:) :: voxel_node
       INTEGER(sp), allocatable, dimension (:,:,:,:) :: loc_Q,transroot
       integer(sp)::count_nodes
       LOGICAL :: loop1=.true., stressBC, ave, old,oldT,eqDis,ItCrit_root=.true.,once=.true.
       LOGICAL :: moment=.false.,savelast=.false.,switchSolve,tcheck=.false.
       LOGICAL :: tcheck2=.false.,ana_aan
       Type(SparseMatrixType), allocatable, dimension (:) :: plantmatrix !multiple roots

       contains
          SUBROUTINE IniMat
	   use RootData, only :nrec
	   IMPLICIT NONE
	   IF (loop1) THEN
	      loop1=.false.
	   else
             deallocate(delta2)
             deallocate(delta2old)
             deallocate(numNodes_voxel)
             deallocate(w_dis)
             deallocate(cent)
             deallocate(cp_mean)
             deallocate(Intc)
             deallocate(tempSinkR)
             deallocate(nsub)
             deallocate(w_sub)
             deallocate(l_sub)
             deallocate(l_seg)
             deallocate(GH)
             deallocate(Joutr)
             deallocate(axialRootFlow)
             deallocate(sinkR)
             deallocate(sinkRtemp)
             deallocate(Lr)
             deallocate(Khr)
             deallocate(PHs)
             deallocate(PHr)
             deallocate(PHrOld)
             deallocate(SinkROld)
             deallocate(PHrTemp)
             deallocate(PHr_sub)
             deallocate(PHsTemp)
             deallocate(PH_root_sub)
             deallocate(Phi_mat)
             deallocate(h_mat2)
             deallocate(loc_Q)
	      deallocate(Qi)
	      deallocate(Q_bc)
             deallocate(Qd)
             deallocate(sum_dis)
             deallocate(PH_micro2)
             deallocate (voxel_no)
             deallocate (no_voxels)
             deallocate(voxel_node)
	      deallocate(transroot)
             if (ave) then
                deallocate(counter2)
             endif
	   endif
          if (ana_aan) then
             if (ave) then
                allocate(counter2(1))
                counter2=0.
             endif
             allocate(PH_micro2(1,1,1)) !1:500 -> length r is lower than 500, so no. of compartments lower
             PH_micro2=0
             allocate (Phi_mat(1))
             Phi_mat = 0._sp
             allocate (h_mat2(1))
             h_mat2 = 0._sp
          else
             if (ave) then
                allocate(counter2(1:Nelm))
                counter2=0.
             endif
             allocate(PH_micro2(0:nrec,1:500,1:isubmax)) !1:500 -> length r is lower than 500, so no. of compartments lower
             PH_micro2=0
             allocate(numNodes_voxel(0:nrec,1:isubmax))
             numNodes_voxel=0
             allocate (Phi_mat(1:nLibr-1))
             Phi_mat = 0._sp
             allocate (h_mat2(1:nLibr-1))
             h_mat2 = 0._sp
          endif
          if (.not.(old)) then
             if ((ave) .or. (eqdis)) then
                allocate (voxel_no(0:nrec,1:isubmax))
                voxel_no=0
                allocate (no_voxels(1:nElm))
                no_voxels=0
                allocate (voxel_node(1:nElm,1:2,indexValue))
                voxel_node=0
                allocate(numNodes_voxel(0:nrec,1:isubmax))
                numNodes_voxel=0
             else
                allocate(numNodes_voxel(1,1))
                numNodes_voxel=0
                allocate (voxel_no(1,1))
                voxel_no=0
                allocate (no_voxels(1))
                no_voxels=0
                allocate (voxel_node(1,1,1))
                voxel_node=0
             endif
          else
             allocate(numNodes_voxel(1,1))
             numNodes_voxel=0
             allocate (voxel_no(1,1))
             voxel_no=0
             allocate (no_voxels(1))
             no_voxels=0
             allocate (voxel_node(1,1,1))
             voxel_node=0
          endif
          allocate(delta2(1:nrec+1))
          delta2=0
          allocate(delta2old(1:nrec+1))
          delta2old=0
          allocate(tempSinkR(1:nrec+1))
          tempSinkR=0
	   allocate (nsub(0:nrec+ngrow,1:nplant))
	   nsub=0._dp
	   allocate (l_seg(0:nrec))
	   l_seg=0._dp
	   allocate (w_sub(0:nrec,1:isubmax,1:nplant))
	   w_sub=0._dp
	   allocate (l_sub(0:nrec,1:isubmax,1:nplant))
	   l_sub=0._dp
          allocate (PHr_sub(0:nrec,1:isubmax))
	   PHr_sub=0._dp
          allocate (PH_root_sub(0:nrec,1:isubmax))
	   PH_root_sub=0._dp
	   allocate (GH(0:nrec))
	   GH=0._dp
	   allocate (sinkR(0:nrec,1:nplant))
	   sinkR=0._dp
          allocate (BCtp_usr(1:nplant))
          BCtp_usr=0
          allocate (curr_BCtp(1:nplant))
          curr_BCtp=0
          allocate (curr_BCr(1:nplant))
          curr_BCr=0
          allocate (BCr_usr(1:nplant))
          BCr_usr=0
          allocate (sinkRtemp(0:nrec))
	   sinkRtemp=0._dp
          allocate (Joutr(0:nrec))
          Joutr=0._dp
          allocate (axialRootFlow(0:nrec,1:nplant))
          axialRootFlow=0._dp
	   allocate (Lr(0:nrec))
	   Lr=0._dp
	   allocate (Khr(0:nrec))
	   Khr=0._dp
	   allocate (PHs(0:nrec,1:nplant))
	   PHs=0._dp
	   allocate (PHr(1:nrec+1,1:nplant))
          PHr=0._dp
          allocate (PHrOld(1:nrec+1,1:nplant))
          PHrOld=0._dp
          allocate (SinkROld(0:nrec,1:nplant))
          SinkROld=0._dp
          allocate (PHrTemp(1:nrec+1,1:nplant))
          PHrTemp=0._dp
          allocate (PHsTemp(0:nrec,1:nplant))
          PHsTemp=0._dp
	   allocate (transroot(0:nrec+ngrow,1:2,1:isubmax,1:nplant))
	   transroot=0
	   allocate (loc_Q(0:nrec,1:8,1:isubmax,1:nplant))
	   loc_Q=0
	   allocate (w_dis(0:nrec,1:8,1:isubmax,1:nplant))
	   w_dis=0._dp
          allocate (cent(0:nrec,1:3,1:isubmax,1:nplant))
          cent=0._dp
          allocate (cp_mean(0:nrec,1:3,1:isubmax,1:nplant))
          cp_mean=0._dp
          allocate (Intc(0:nrec+ngrow,1:3,1:isubmax,1:nplant))
          Intc=0._dp
	   allocate (Qi(0:nrec,1:nplant))
	   Qi=0._dp
	   allocate (Q_bc(0:nrec,1:nplant))
	   Q_bc=0._dp
	   allocate (Qd(0:nrec,1:nplant))
	   Qd=0._dp
	   allocate (sum_dis(0:nrec,1:isubmax,1:nplant))
	   sum_dis=0._dp
          allocate (k_ave(0:nrec,1:isubmax))
	   k_ave=0._dp
          allocate (B(0:nrec,1:isubmax))
	   B=0._dp
          allocate (plantmatrix(1:nplant))
	   END SUBROUTINE IniMat
       END MODULE DoussanMat
!*******************************************************************
 Module NumericalRecipes
      USE Typedef
       
      Contains
       
      SUBROUTINE linbcg(n,b,x,itol,tol,itmax,iter,err,ipl)
! solve inverse matrix with preconditionned biconjugate gradient method      
! code taken from Numerical recipes in fortran 77, p. 79
      USE SparseMatrix
      USE DoussanMat, only: plantmatrix
      IMPLICIT NONE
      
      INTEGER(sp), intent(in) :: itmax,itol,n,ipl
      INTEGER(sp), intent(out) :: iter
      REAL(dp), intent(inout):: tol,b(:),x(:)
      REAL (dp), intent(out) ::err
      INTEGER (sp) :: j
      REAL(dp) :: ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm,EPS
      REAL(dp) :: p(n),pp(n),r(n),rr(n),z(n),zz(n)
      PARAMETER (EPS=1.E-14_dp)
      
      iter=0     
      call SM_multiply(plantmatrix(ipl), x, r)
      do 11 j=1,n
         r(j)=b(j)-r(j)
         rr(j)=r(j)
11    continue
      znrm=1._dp
      if (itol.eq.1) then
         bnrm=snrm(n,b,itol)
      else if (itol.eq.2) then
         call SM_divide_by_diagonal(plantmatrix(ipl), b,z)
         bnrm=snrm(n,z,itol)
      else if (itol.eq.3.or.itol.eq.4) then
         call SM_divide_by_diagonal(plantmatrix(ipl), b,z)
         bnrm=snrm(n,z,itol)
         call SM_divide_by_diagonal(plantmatrix(ipl), r,z)
         znrm=snrm(n,z,itol)
      else
	  print *,'illegal itol in linbcg'
      endif
      call SM_divide_by_diagonal(plantmatrix(ipl), r,z)     

100   if (iter.le.itmax) then
         iter=iter+1
         zm1nrm=znrm
         call SM_divide_by_diagonal(plantmatrix(ipl), rr,zz)
         bknum=0._dp
         do 12 j=1,n
            bknum=bknum+z(j)*rr(j)
12       continue
         if(iter.eq.1) then
            do 13 j=1,n
               p(j)=z(j)
               pp(j)=zz(j)
13          continue
         else
            bk=bknum/bkden
            do 14 j=1,n
               p(j)=bk*p(j)+z(j)
               pp(j)=bk*pp(j)+zz(j)
14          continue
         endif
         bkden=bknum
!        call atimes(n,p,z,0)
         call SM_multiply(plantmatrix(ipl), p, z)
	  akden=0._dp
         do 15 j=1,n
            akden=akden+z(j)*pp(j)
15       continue
         ak=bknum/akden
!         call atimes(n,pp,zz,1)

         call SM_multiply_transpose(plantmatrix(ipl), pp, zz)
         do 16 j=1,n
            x(j)=x(j)+ak*p(j)
            r(j)=r(j)-ak*z(j)
            rr(j)=rr(j)-ak*zz(j)
16       continue  
!print *, '3496',x(1:5),'ak',ak,'z',z(1:5),'zz',zz(1:5),'p',p(1:5)
         call SM_divide_by_diagonal(plantmatrix(ipl), r,z)
         if(itol.eq.1.or.itol.eq.2)then
            znrm=1._dp
            err=snrm(n,r,itol)/bnrm
         else if(itol.eq.3.or.itol.eq.4)then
            znrm=snrm(n,z,itol)
            if(abs(zm1nrm-znrm).gt.EPS*znrm) then
               dxnrm=abs(ak)*snrm(n,p,itol)
               err=znrm/abs(zm1nrm-znrm)*dxnrm
            else
               err=znrm/bnrm
               goto 100
            endif
            xnrm=snrm(n,x,itol)
            if(err.le.0.5_dp*xnrm) then
               err=err/xnrm
            else
               err=znrm/bnrm
               goto 100
            endif
         endif
         if(err.gt.tol) then
            goto 100
         endif
      endif	
      return      
      END SUBROUTINE linbcg
!  (C) Copr. 1986-92 Numerical Recipes Software '%12'%*Sim+).
!***************************************************************************
      FUNCTION snrm(n,sx,itol)
! compute one or two norms of a vector sx(1:n), as signaled by itol. used by linbcg
! from Numerical Recipes in F., p.81      
      INTEGER(sp):: n,itol,i,isamax
      REAL(dp):: sx(:),snrm
      
      if (itol.le.3)then
         snrm=0._dp
         do 11 i=1,n
            snrm=snrm+sx(i)**2
11       continue
         snrm=sqrt(snrm)
      else
         isamax=1
         do 12 i=1,n
            if(abs(sx(i)).gt.abs(sx(isamax))) isamax=i
12       continue
         snrm=abs(sx(isamax))
      endif
      return
      END FUNCTION snrm

    End Module NumericalRecipes
!========================================================================
      MODULE tmctrl
      Use typedef
      IMPLICIT NONE
      INTEGER(sp), PARAMETER :: mxoFEM=50,mxoRoo=50,mxProf=1000000
      REAL(sp) :: tOuFEM(mxoFEM),tOuRoo(mxoRoo),touProf(mxProf),touProbe(mxProf)
      REAL(sp) :: dtroot,dtMin,dtMax,FacInc,FacDec,tmax,t_begin,dtProf,dtProbe
      INTEGER(sp) :: nOuFEM,nOuRoo,nouProf,nouProbe
      END MODULE tmctrl
!*******************************************************************
      MODULE ConData
      USE typedef
      USE ParamData
      IMPLICIT NONE
! conpar,conimp
      REAL(sp):: impc(maxnod),coptma,cmax,cmin,coptmi
      END MODULE ConData
!********************************************************************
      MODULE BoundData
!data for boundary conditions
!COMMON /bndary/      /cBound/
      Use typedef
      USE ParamData
      IMPLICIT NONE
      REAL(sp) :: cBound(2),CBnd2(mxBcCh)
      REAL(sp) :: tQbcCh(mxBcCh),tIbcCh(mxBcCh),Ibc(mxBcCh),thbcCh(mxBcCh),hbc(mxBcCh)
      REAL(sp) :: tCBnd1(mxBcCh),CBnd1(mxBcCh),tCBnd2(mxBcCh)
      REAL(sp), dimension (:,:) :: Qbc(mxBcCh,mxBcCh)
      INTEGER(sp) :: iBCPt(maxbdr+maxIrrig),nQbcCh,nIbcCh,nCBnd1,nhbcCh,nCBnd2,homogene
      END MODULE BoundData
!********************************************************************
      MODULE TempData
      USE ParamData
      USE typedef
      IMPLICIT NONE
!data for temperature
!common/tempar//TempTm/ /tmptre/
      REAL(sp):: time(mxtime),depth(mxdpth),temtim(mxtime,mxdpth)
      REAL(sp):: tem(maxnod),impt(maxnod)
      REAL(sp):: tmin,topt,tmax,trange,tmid,expo
      INTEGER(sp):: nt,nz_temp
      END MODULE TempData
!********************************************************************
      MODULE CumData
      USE Typedef
      USE Paramdata
      IMPLICIT NONE
!common flow,cmconc,error      tensor, list
      REAL(sp) :: Qc(maxnod),Q(maxnod)
      REAL(dp) :: cumCh0=0._dp,cumChr=0._dp,cumCh1=0._dp,CumRt=0._dp
      REAL(dp),dimension(2):: CumQ=(/0._dp,0._dp/),ChemS=(/0._dp,0._dp/)
      REAL(dp) :: wCumA=0._dp,wCumT=0._dp,cCumA,cCumT,wVolI,cVolI,wBalR
      REAL(dp) :: VolSink=0._dp, VolQ=0._dp, WatVolOld=0._dp
      REAL(sp) :: WatIn(maxelm),SolIn(maxelm)
      REAL(sp) :: ConAxx(maxelm),ConAyy(maxelm),ConAzz(maxelm)
      REAL(sp):: ConAxy(maxelm),ConAxz(maxelm),ConAyz(maxelm)
      INTEGER(sp) :: ListNE(maxnod)
      END MODULE CumData
!********************************************************************
      MODULE DomData
!common Domlim
       USE TypeDef
       USE ParamData
       IMPLICIT NONE
       REAL(sp):: xmin,xmax,ymin,ymax,zmin,zmax
       END MODULE DomData
!********************************************************************
      MODULE GeoData
      USE TYPEDEF
! common  /axspar//latpar/
      IMPLICIT NONE
      INTEGER(sp), PARAMETER :: maxaxs=100,mxpnts=100
      REAL(sp) :: geoaxs(maxaxs),angaxs(maxaxs,mxpnts)
      REAL(sp) :: tempax(maxaxs,mxpnts)
      REAL(sp):: anglat(mxpnts),templt(mxpnts),geolat
      INTEGER(sp) ::nangax(maxaxs),nanglt
      END MODULE GeoData
!********************************************************************
       MODULE MatData
       USE TypeDef
       USE GridData, only :nband,nPt 
       USE ParamData
       IMPLICIT NONE
       integer(sp) :: NumNZ
       REAL(sp), allocatable, dimension (:,:) :: A
      ! REAL(sp):: A(maxbnd,maxnod)
       real(dp):: A_sparse(maxamg2)
       real(dp):: B(maxnod)
       integer(sp):: IROW(maxamg2),JCOL(maxamg2)
       REAL(sp) :: As(maxbnd,maxnod),Bs(maxnod)
    !   LOGICAL :: loop1=.true.
!
 !      contains
  !           SUBROUTINE IniMatData
!	             allocate(A(1:nband,1:nPt))
 !            END SUBROUTINE IniMatData

	   END MODULE MatData
!********************************************************************
       MODULE ObsData
       USE TypeDef
       USE Paramdata
	   !output variables for observation probes
	   ! npr: number of probes
          ! Pr: Probe ID
	   ! Pt: plane direction (perp to X (1), to Y (2) or to Z (3)) or probe given by the user (4))
	   ! CrP: crossing point or number of user nodes (if Pt=4)
	   ! VarP: variable (WC=1, PH=2, both=3)
	   ! distrP: 1=all, 2=average
       IMPLICIT NONE
	   INTEGER(sp) :: npr, Pt(maxobs),nodebyPr(maxobs),nprof,Pr(maxobs)
	   INTEGER(sp) :: NodePr(maxobs,1000) !1000=max number of node for a given plane
	   INTEGER(sp) :: VarP(maxobs),distrP(maxobs),Pl(maxobs),varProf(maxobs)
	   REAL(dp) ::CrP(maxobs)
	   END MODULE ObsData
!********************************************************************
     MODULE SolData
       USE Typedef
       USE ParamData
       IMPLICIT NONE
       REAL(sp) :: ChPar(10,maxmat),par(11,maxmat),epsi,Peclet,Courant
       REAL(sp) :: PeCr,Fc(maxnod),Gc(maxnod)
       INTEGER(sp) ::NLevel,nMat
       REAL(sp) :: Vx(maxnod),Vy(maxnod),Vz(maxnod)
       REAL(sp) :: Dispxx(maxnod),Dispyy(maxnod),Dispzz(maxnod)
       REAL(sp) :: Dispxy(maxnod),Dispxz(maxnod),Dispyz(maxnod)
       REAL(sp):: conO(maxnod),con(maxnod),cap(maxnod)
      END MODULE SolData
!********************************************************************
       MODULE StrData
       USE TypeDef
       USE ParamData
       IMPLICIT NONE
!strpar and strgth
       REAL (sp):: s(maxnod),imps(maxnod),simp,refgrd
       END MODULE StrData
!********************************************************************
      MODULE WatFun
      USE Typedef
      USE ParamData
      IMPLICIT NONE
!COMMON FILE:Table
      REAL(sp) :: hTab(nTab),ConTab(nTab,maxmat),CapTab(nTab,maxmat)
      REAL(sp) :: TheTab(nTab,maxmat)
      REAL(sp) ::alh1,dlh
      CONTAINS

! current local theta (volumetric soil water content) values:
   REAL(sp) FUNCTION Fth(h,Parloc)
      REAL(dp) :: thr,ths,THT
      REAL(sp),intent(in) :: h,Parloc(:)

      thr=Parloc(2)
      ths=Parloc(3)
      IF (h.LT.0.0_dp) THEN
         THT=FThNrm(h,Parloc)
         Fth=thr+THT*(ths-thr)
      ELSE
         Fth=ths
      ENDIF
      RETURN
      END FUNCTION Fth
!********************************************************************
! current local THETA (0...1 normalized theta) values:
! updated Auguts 09 MJ

      REAL(sp) FUNCTION FThNrm(h,Parloc)
      REAL(dp) :: a,n,m,a2,n2, m2,w1,w2
      REAL(sp),intent(in) :: h,Parloc(:)
      INTEGER :: MatMod

      a=Parloc(4)
      n=Parloc(5)
      m=1._dp-1._dp/n
      MatMod=Parloc(1)

      SELECT CASE(MatMod)
         CASE (1) ! Mualem Van Genuchten
            IF (h.LT.0.0_dp) THEN
               FthNrm=(1._dp+(a*ABS(h))**n)**(-m)
            ELSE
               FthNrm=1.0_dp
            ENDIF
         CASE (5) !Dual porosity model
            w2=parloc(8)
            a2=parloc(9)
            n2=parloc(10)
            w1=1._dp-w2
            m2=1._dp-1._dp/n2
            IF (h.LT.0.0_dp) THEN
               FthNrm=w1*(1._dp+(a*ABS(h))**n)**(-m) + w2*(1._dp+(a2*ABS(h))**n2)**(-m2)
            ELSE
               FthNrm=1.0_dp
         ENDIF
      END SELECT
      RETURN
      END FUNCTION FThNrm
!********************************************************************
! current local conductivity values:
! updated Auguts 09 MJ

  REAL(sp) FUNCTION FKP(h,Parloc)
      REAL(dp) ::  Ks,THT,lambda,a,a2,n,n2,m,m2,w1,w2,Sv1,Sv2,rNumer,rDenom
      REAL(sp),intent(in) :: h,Parloc(:)
      INTEGER :: MatMod
      
      a=Parloc(4)
      n=Parloc(5)
      Ks=Parloc(6)
      lambda=Parloc(7)
      m=1._dp-1._dp/n

      MatMod=Parloc(1)
      SELECT CASE(MatMod)
         CASE (1) ! Mualem Van Genuchten
            IF (h.LT.0.0_dp) THEN
               THT=(1._dp+(a*ABS(h))**n)**(-m)
               FKP=Ks*(THT**lambda)*(1._dp-(1._dp-THT**(1._dp/m))**m)**2
            ELSE
               FKP=Ks
            ENDIF
         CASE (5) !Dual porosity model
            w2=parloc(8)
            a2=parloc(9)
            n2=parloc(10)
            w1=1.0-w2
            m2=1._dp-1._dp/n2
            IF (h.LT.0.0_dp) THEN
               THT=w1*(1._dp+(a*ABS(h))**n)**(-m)+w2*(1._dp+(a2*ABS(h))**n2)**(-m2)
               Sv1=(a*ABS(h))**(n-1)
               Sv2=(a2*ABS(h))**(n2-1)
               rNumer=w1*a*(1._dp-Sv1*(1._dp+(a*ABS(h))**n)**(-m)) +&
                  w2*a2*(1._dp-Sv2*(1._dp+(a2*ABS(h))**n2)**(-m2)) 
               rDenom=w1*a+w2*a2
               FKP=Ks*(THT**lambda)*(rNumer/rDenom)**2
            ELSE
               FKP=Ks
            ENDIF            
      END SELECT
      RETURN
      END FUNCTION FKP
!********************************************************************
! current local soil water capacity values:
! updated Auguts 09 MJ

   REAL(sp) FUNCTION FCP(h,Parloc)
      REAL(dp) :: thr,ths,a,n,m,C2a,C2b,m2,w1,W2,a2,n2
      INTEGER :: MatMod
      REAL(sp), Intent(in) :: h,Parloc(:)
      
      MatMod=Parloc(1)
      thr=parloc(2)
      ths=parloc(3)
      a=parloc(4)
      n=parloc(5)
      m=1._dp-1._dp/n
      SELECT CASE(MatMod)
         CASE (1) ! Mualem Van Genuchten
            IF (h.LT.0.0_dp) THEN
               FCP=(ths-thr) *(a*n*m*((a*ABS(h))**(n-1._dp)))/((1._dp+(a*ABS(h))**n)**(m+1._dp))
            ELSE
               FCP=0.0_dp
            ENDIF
         CASE (5) !Dual porosity model
            w2=parloc(8)
            a2=parloc(9)
            n2=parloc(10)
            w1=1._dp-w2
            m2=1._dp-1._dp/n2
            IF (h.LT.0.0_dp) THEN
               C2a=(ths-thr) *(a*n*m*((a*ABS(h))**(n-1._dp)))/((1._dp+(a*ABS(h))**n)**(m+1._dp))
               C2b=(ths-thr) *(a2*n2*m2*((a2*ABS(h))**(n2-1._dp)))/((1._dp+(a2*ABS(h))**n2)**(m2+1._dp))
               FCP= w1*C2a+w2*C2b
            ELSE
               FCP=0.0_dp
            ENDIF            
      END SELECT
      RETURN
      END FUNCTION FCP
      END MODULE WatFun

!====================================================================
!main program
!====================================================================
      USE typedef
      USE CumData
      USE PlntData
      USE tmctrl
      USE GridData
      USE RootData
      USE MatData
      USE DoussanMat
      USE SparseMatrix
      USE soldata
      USE ObsData
      USE DomData
      USE BoundData, only : tQbcCh
      IMPLICIT NONE

!variable declaration
      LOGICAL temp,lChem,ReDo,toxi,it1,lDou,lOrt
      REAL(sp) :: dtMaxC=1.E+30_dp,delbm,w,LA,LAmsh,mroot,mshoot
      REAL(sp) ::theta_old(maxnod),tcBCr,tProf,tProbe
      REAL(sp) :: tcallr,dmroot,told,dtold,tfem,troo,rs,concrs
      REAL (dp) ::t0,t1!,delta2
      REAL(sp) :: hOld(maxnod),hTemp(maxnod),hNew(maxnod),theta(maxnod)
      REAL(dp) ::Conc(maxnod)
      REAL(sp) :: random,dt,tpulse,t,savg,cavg,dtopt,rsr,dr,dsh,grwfac, maxdelta2
      INTEGER(sp) :: icount,iter,i,level,naxes,norder,kOuRoo=0,koufem=0,kouprof=0,kouprobe=0
      INTEGER(sp) :: iter_tot=0,kBCr=0,ind,daytime(8),ip,ipl
      INTEGER(sp) :: Kode(maxnod),KodCB(maxbdr+maxIrrig),MatNum(maxnod),kaxemg=0
      INTEGER(sp):: IAD(maxbnd,maxnod),IADN(maxnod),IADD(maxnod),iter_root
      logical :: root_growth,itCrit,stressBC_old,ObsOK,profOK,no_root
      CHARACTER form*3,file*8

      common /root_grow/ root_growth
! initiate random number generator:
!      random=secnds(.0_dp)
!      irnd=2*INT(random)+1
      call random_number(random)
      irnd=2*INT(random)+1
      call date_and_time(values=daytime)
!      irnd=2*(daytime(5)*3600+daytime(6)*60+daytime(7))+1
 print *,'need random number generator',irnd
!      print *,'1',irnd
!      irnd=2*JINT(random)+1
!      print *,'2',irnd
!      CALL random_seed(put=1)
!      CALL random_generator(harvest=rd)
!      irnd=rd(1)
!      print *,'2',irnd
!open simulation summary file
      OPEN (UNIT=15,FILE='out/simul_sum.out',STATUS='UNKNOWN')
      WRITE(15,112)
112   FORMAT('+++++ SIMULATION INPUT SUMMARY +++++ ')
      CLOSE(15)
      ObsOK=.FALSE.

!-----------------------------------INPUT---------------------------------------
! get input for the specific problem (domain, BC, IC, control parameters):
      CALL Applic(hOld,hNew,hTemp,dt,tPulse,Kode,KodCB,MatNum,Conc,lOrt,ObsOK,profOK,no_root)
! calculate Wn (new Javaux)
      CALL CalcWnodes
! get solute transport information
      CALL ChemIn(lChem)
! get soil material input and set up K,C-table for interpolation:
      CALL SoilIn
      IF (.not.(no_root)) THEN
! initial root system and growth parameters:
         CALL RootIn(naxes,norder,t,mroot,mshoot,LA,sAvg,cAvg,toxi)
! get Doussan model input information
         CALL DouIn(lDou)
! get plant parameter input:
         CALL PlntIn(level,lChem,toxi,lDou,t)
! soil temperature over depth / in time:
         CALL TempIn(temp)

!--------------------------------- OUTPUT FILES ------------------------------------
! open log file:
         WRITE (file,'(A7)')'out/log'
         DO ipl=1,nplant
            WRITE (file(8:8),'(I1)') ipl
            OPEN (UNIT=10,FILE=file,STATUS='UNKNOWN')
            WRITE (10,90)
90          FORMAT('  Time       Tpot       Tact   soil_lim_nodes    PH_collar  grwfac     sAvg       cAvg       mShoot     mRoot')
            CLOSE (10)
         ENDDO
      ELSE
         lDou=.FALSE.
         level=0 !no RWU, no root growth
         t=tQbcCh(1)!initial time =first time of the soil BC.
      ENDIF

! open balance file:
      OPEN (UNIT=10,FILE='out/balance.out',STATUS='UNKNOWN')
      WRITE(10,110)
110   FORMAT(' Absolute and relative mass balance error for water and solute transport ')
      CLOSE (10)
! open removal file
      OPEN (UNIT=10,FILE='out/remove.out',STATUS='UNKNOWN')
      WRITE(10,111)
111   FORMAT(' Total amount of water and solute removed by the root system, ',/,&
           ' by zero and first order reactions, and by drainage at the bottom of the domain. ')
      CLOSE(10)
! open Observation node file
      IF (ObsOK) THEN
         CALL ObsIni
      ENDIF
      IF (profOK) THEN
         OPEN (UNIT=121,FILE='out/ProfileTH.out',STATUS='UNKNOWN')
         OPEN (UNIT=122,FILE='out/ProfilePH.out',STATUS='UNKNOWN')
         OPEN (UNIT=123,FILE='out/ProfileS.out',STATUS='UNKNOWN')
         WRITE(form,'I3')nz
         WRITE (121,'(/''Averaged water content profile for each time step.'')')
         WRITE (121,'(/''Time   Z- water content'')')
         WRITE (121,'((12X)'//form//'(1X,f12.4))') (zmax-dzgrid*(i-1),i=1,nz)
         WRITE (122,'(/''Averaged water potential profile for each time step.'')')
         WRITE (122,'(/''Time   Z- water potential'')')
         WRITE (122,'((12X)'//form//'(1X,f12.4))') (zmax-dzgrid*(i-1),i=1,nz)
         WRITE (123,'(/''Total Sink profile for each time step.'')')
         WRITE (123,'(/''Time   Z- sink'')')
         WRITE (123,'((12X)'//form//'(1X,f12.4))') (zmax-dzgrid*(i-1),i=1,nz)
      ENDIF
      IF(lOrt) call IADMake(nPt,nElm,maxElm,maxbnd,IAD,IADN,IADD)
!--------------------------------------------------------------------------------

! first time step:
      dtOpt=dt
      tCallR=t+dtRoot
      dmroot=0.0_dp

! find next output time:
    1 kOuRoo=kOuRoo+1
      IF (kOuRoo.LE.nOuRoo) THEN
         IF (tOuRoo(kOuRoo).LE.t) GOTO 1
      ENDIF
    2 kOuFEM=kOuFEM+1
      IF (kOuFEM.LE.nOuFEM) THEN
         IF (tOuFEM(kOuFEM).LE.t) GOTO 2
      ENDIF
    3 kaxemg=kaxemg+1
      IF (kaxemg.LE.naxemg) THEN
         IF (tnewax(kaxemg).LT.t) GOTO 3
      ENDIF
! Z-profiles
      touProf=touProf+t !to start at initial time < rootin
      IF ((ProfOK).AND.(dtProf.NE.999)) THEN
    5    kouProf=kouProf+1
         IF (kouProf.LE.nouProf) THEN
            IF (touProf(kouProf).LE.t) GOTO 5
         ENDIF
      ENDIF
! probes
      touProbe=touProbe+t !to start at initial time < rootin
      IF ((ObsOK).AND.(dtProbe.NE.999)) THEN
    6    kouProbe=kouProbe+1
         IF (kouProbe.LE.nouProbe) THEN
            IF (touProbe(kouProbe).LE.t) GOTO 6
         ENDIF
      ENDIF
      t_begin=t

! if Doussan, calculate weighing function and matrices
      
      IF (lDou) THEN
         CALL SetupDou(t,naxes,level,dt)
         IF (level==4) THEN
         ! find next time change in root BC
    4       kBCr=kBCr+1
            IF (kBCr.LE.nBCr) THEN
               IF (tBCr(kBCr).LE.t) GOTO 4
            ENDIF
         ENDIF
         !initialize stress
         stressBC=.FALSE.
      ELSE
     ! get betaw and betac initial distribution:
         CALL BetDis(t)
     ! ...and normalize:
         CALL BetNrm
      ENDIF

! calculate parameters:
      CALL SetMat(0,hOld,hTemp,hNew,MatNum,theta)
      IF (lChem) THEN
         CALL ChInit(MatNum,hNew,dtMaxC,dt,theta,Dxy,Exy)
      ENDIF
      CALL CalcGeom
      CALL CalcWidthSWMS
      iCount=0
      CALL SubReg(MatNum,Conc,t,hNew,iCount,lChem,theta,level)

! initialize it1 (first iteration after a root growth)
      it1=.TRUE.
      itMaxRoot=itMax
      IF (oldT) then
         switchsolve=.false.
      ELSE
         switchsolve=.true.
      ENDIF

! --------------------------- Start of time loop ------------------------------
!call cpu_time(t0)
  200 CONTINUE
! time step adjustment:
      tOld=t
      dtOld=dt
      IF (kOuFEM.LE.nOuFEM) THEN
         tFEM=tOuFEM(kOuFEM)
      ELSE
         tFEM=tMax
      ENDIF

      IF (kOuRoo.LE.nOuRoo) THEN
         tRoo=tOuRoo(kOuRoo)
      ELSE
         tRoo=tMax
      ENDIF
      IF (level.EQ.4) THEN
         IF (kBCr.LE.nBCr) THEN
            tcBCr=tBCR(kBCr)
         ELSE
            tcBCr=tMax
         ENDIF
      ELSE
         tcBCr=tMax
      ENDIF
      IF ((ProfOK).AND.(dtProf.NE.999).AND.(kouProf.LE.nouProf)) THEN
         tProf=tOuProf(kouProf)
      ELSE
         tProf=tMax
      ENDIF
      IF ((ObsOK).AND.(dtProbe.NE.999).AND.(kouProbe.LE.nouProbe)) THEN
         tProbe=tOuProbe(kouProbe)
      ELSE
         tProbe=tMax
      ENDIF
      CALL TmCont(iter,t,dt,dtOpt,tCallR,tFEM,tRoo,dtMaxC,tcBCr,tProf,tProbe)

      t=t+dt
  222 CONTINUE
      WRITE (*,'(/'' Time ='',1pE12.5)') t
      iter=0
      iter_root=0
! calculate current relative stress from average soil strength:
      IF (level.EQ.3) CALL Stress(sAvg,cAvg,rs,concrs,lChem)

! calculate current values of Tpot and time-dependent B.C.:
      CALL SetTp(t,rs,concrs,LA,level,lChem)
      CALL SetBC(t,hNew,Kode)

! apply root extraction to FEprint *,'time2'M-sink term for each node:
      IF (.not.(switchSolve)) THEN
         IF (lDou) THEN! 
            WRITE (*,'(''+'',$)')
            CALL SolveRoot(hNew,hold,t,dt,it1,level,iter_root)
            CALL Setsnk(hNew,Conc,lDou)
            it1=.FALSE.
         ENDIF
      ENDIF

! initialize a few variable for solute transport:
      IF (lChem) THEN
         IF (iCount.EQ.0) CALL ChInit(MatNum,hNew,dtMaxC,dt,theta,Dxy,Exy)
      ENDIF
! solve soil water flow equation and calculate pressure heads:
      if (.not.(switchSolve)) then
         CALL Water(hNew,hOld,hTemp,t,dt,dtOpt,tOld,MatNum,Kode,iter,iter_tot,ReDo,theta,level,lDou, &
IAD,IADN,IADD,lOrt)
      elseif (switchSolve) then
         CALL Water2(hNew,hOld,hTemp,t,dt,dtOpt,tOld,MatNum,Kode,iter,iter_tot,ReDo,theta,theta_old,level,lDou, &
IAD,IADN,IADD,lOrt)
      endif
      IF (ReDo) THEN
         savelast=.false.
         GOTO 222
      ENDIF
! solve solute transport equation and calculate concentrations:
      IF (lChem) THEN
         CALL Solute(dt,t,Conc,tPulse,hOld,hNew,dtMaxC,MatNum,Kode,KodCB,theta,Dxy,Exy,theta_old)
      ENDIF

! calculate water and solute (if included) mass balance
      iCount=1
      CALL SubReg(MatNum,Conc,t,hNew,iCount,lChem,theta,level)

! calculate actual transpiration rate Tact (global variable in MODULE PlntData):
      IF ((level.GE.2).AND.(.not.(lDou))) CALL ActTrs
      IF (level.EQ.3) THEN
         DO ipl=1,nplant !the subroutines are not adjusted yet to work with multiuple plants
! translate transpired water into biomass:
            CALL Effncy(t,rs,concrs,W,lChem)
            delBM=W*Tact(ipl)*dt

! partition that goes to new roots is accumulated in dtroot
! until next root growth step:
            CALL Ratio(t,rs,concrs,RSR,lChem) 
            dr=RSR/(1.+RSR)*delBM
            dmroot=dmroot+dr

! remainder of delBM is partition that goes to shoot:
            dsh=delBM-dr
            mshoot=mshoot+dsh

! current leaf area:
            CALL Leaves(t,LAmsh)
            LA=LA+dsh*LAmsh
         END DO
      ENDIF

!save Doussan data before deleting of vector and root data
      IF ((ABS(t-tOuRoo(kOuRoo)).LE.0.001_dp*dt).AND.(level.eq.4)) THEN !MJ09
         ! Output for Doussan
         CALL OutDou(t,kOuRoo)
         write(*,*)'outRoo. correctly written'
      ENDIF
      root_growth=.false.
! if appropriate, grow root and update the distribution functions 'betaw'
! and 'betac' (extraction intensity per volume, [L^-3])
      IF (ABS(t-tCallR).LT.0.001_dp*dt) THEN
         root_growth=.true.

!update soil strength values:
         CALL Solstr(hNew,MatNum)

! if temperature input provided, update temperature values:
         IF (temp) CALL Temper(t)

!if nutrient deficiency/ion toxicity data are provided, update
! the corresponding impedance nodal values:
         IF (toxi) CALL ConTox(Conc)
!print *,'l1053'
! let root system grow:
         CALL Root(Conc,t,dmroot,mroot,sAvg,cAvg,grwfac,naxes,norder,kaxemg,level,temp,toxi)
         tCallR=tCallR+dtRoot

! update Dousssan weighing factors and matrices
	  IF (lDou) THEN
! delete the old vector of list
            IF (nrecold.NE.0) THEN
	        DO ipl=1, nplant
		    CALL SM_delete(plantmatrix(ip))!multipleroot
               END DO
            ENDIF
!            deallocate(listvec)
	     deallocate(nBC_irecn)
!            if ((ave) .or. (eqDis)) deallocate(voxel_node)
	     CALL SetupDou(t,naxes,level,dt)
	     it1=.TRUE.
	  ELSE
! get betaw and betac distribution:
            CALL BetDis(t)
! normalize above so that integral of betaw and betac over 3D-domain is equal to unity:
            CALL BetNrm
	  ENDIF
      ENDIF
!print *,'l1078'
!general output:
      CALL WriteLog(t,grwfac,sAvg,cAvg,mshoot,mroot,level)
!ouputs for observation probes
      IF (ObsOK) THEN
         IF (dtprobe.EQ.999) THEN
            CALL OutObsProbe(t,hNew,MatNum)
         ELSEIF (ABS(t-tOuProbe(kOuProbe)).LE.0.001_dp*dt) THEN
            CALL OutObsProbe(t,hNew,MatNum)
            kouProbe=kouProbe+1
         ENDIF
      ENDIF

!z-profiles
      IF (profOK) THEN
         IF (dtprof.EQ.999) THEN
            CALL Zprofiles(t,hnew,MatNum)
         ELSEIF (ABS(t-tOuProf(kOuProf)).LE.0.001_dp*dt) THEN
            CALL Zprofiles(t,hnew,MatNum)
            kouProf=kouProf+1
         ENDIF
      ENDIF

!FEM, root output at specified points in time: (independent of DoussanMat parameters
      IF (ABS(t-tOuFEM(kOuFEM)).LE.0.001_dp*dt) THEN
         CALL OutFEM(t,hNew,Conc,MatNum,kOuFEM)
  	  CALL FlxOut(kOuFEM,hNew) !SWMS3D
         kOuFEM=kOuFEM+1
      ENDIF

! check next BC time for root
      IF (level.eq.4) THEN
         IF (ABS(t-tBCr(kBCr)).LE.0.001_dp*dt) THEN
            kBCr=kBCr+1
         ENDIF
      ENDIF
      IF (ABS(t-tOuRoo(kOuRoo)).LE.0.001_dp*dt) THEN
         CALL OutRoo(t,naxes,mroot,mshoot,LA,sAvg,cAvg)
         kOuRoo=kOuRoo+1

! Output for Doussan
         IF (level.eq.4) THEN !if OutDou is called list may be empty, due to delete compleet
!	     CALL OutDou(t,kOuRoo)
      	  ENDIF
      ENDIF
      if (savelast) then
         CALL OutFEM(t,hNew,Conc,MatNum,0)
         CALL FlxOut(kOuFEM,hNew) !SWMS3D
         CALL OutDou(t,0)
         CALL OutRoo(t,naxes,mroot,mshoot,LA,sAvg,cAvg)
!         call cpu_time(t1)
         open (unit=9,File='SimulationTime.out',status='unknown')
         write(unit=9,FMT=*) t1-t0
         close(unit=9)
         write(*,*)'Simulation time =',t1-t0
         STOP
      endif
! end of simulation?:
      IF (ABS(t-tMax).LE.0.001_dp*dt) THEN
!         call cpu_time(t1)
         open (unit=9,File='SimulationTime.out',status='unknown')
         write(unit=9,FMT=*) t1-t0
         close(unit=9)
         write(*,*)'Simulation time =',t1-t0
         CALL Getout
      ENDIF

! pressure heads for new time level:
      hOld =hNew
      hTemp=hNew

201   CONTINUE
      GOTO 200
! ------------------------ End of time loop ----------------------------------
    END
! ==============================================================================
! Source file INPUT ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
SUBROUTINE Applic(hOld,hNew,hTemp,dt,tPulse,Kode,KodCB,MatNum,Conc,lOrt,ObsOK,profOK,no_root)
! nPt = total number of nodal points
! nBCpts = total number of nodal points with specified BC
! nElm = total number of elements
! ne* = number of elements (half-cuboids) in '*'-direction
! nel = nex*ney (number of elements per horizontal element-layer)
      USE Typedef
      USE RootData, only: rrt
      USE GridData
      USE tmctrl
      USE BoundData
      USE SolData
      USE CumData
      USE DomData
      USE ObsData
      USE doussanmat, only : old,oldT,counter2,PH_micro2,Phi_mat,h_mat2,ana_aan
      IMPLICIT NONE
      REAL(sp),intent(out) :: hOld(maxnod),hTemp(maxnod),hNew(maxnod)
      REAL(sp),intent(out) :: tpulse,dt
      REAL,intent(out) :: Conc(maxnod)
      REAL(sp) :: Qbcrec(mxBcCh)
      INTEGER(sp),intent(out)::KodCB(maxbdr+maxIrrig),Kode(maxnod),MatNum(maxnod)
      CHARACTER nodes*12,text*5
      INTEGER(sp) :: nfrdr,nh,nQ,nI,l,j,ise,ie,ip,i,idum,k,iL(4),err,ii,nLast,i_old,nl
      INTEGER(sp) :: RelEpsT,i_roottyp,i_reduc,integ_zprof,i2
      INTEGER(sp) :: node_temp(maxnod),i_noroot,i_continu
      INTEGER(sp) :: xIrrig(maxIrrig),yIrrig(maxIrrig),zIrrig(maxIrrig)
      REAL(sp) :: A11,A22,A33,A12,A13,A23,C11,C22,C33,ini,ixmin,ixymin,imin
      logical :: lOrt,ObsOK,profOK,firstOK,no_root

      nMat=0
      xmax=-1.E+30_dp
      xmin=+1.E+30_dp
      ymax=-1.E+30_dp
      ymin=+1.E+30_dp
      zmax=-1.E+30_dp
      zmin=+1.E+30_dp
      profOK=.FALSE.
      OPEN(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
      OPEN (Unit=10,FILE='in/control.in',STATUS='OLD',ERR=10)
      WRITE(15,'(//''++ General Information ++'')',advance='no')
      WRITE(15,'(/''-------------------------'')',advance='no')
      READ (10,*)
      READ (10,*)
      READ (10,*) Headln
      READ (10,*)
      READ (10,*)
      READ (10,*) LnUnit,TmUnit,MsUnit,CnUnit
      READ (10,*)
      READ (10,*)
      READ (10,*) itMax,itMaxRoot
      READ (10,*)
      READ (10,*)
      READ (10,*) RelEpsT,factorRelEps
      READ (10,*)
      READ (10,*)
      READ (10,*) epslonPH,epslonWC,epslonR,epslonS
      IF (RelEpsT .eq. 2) THEN
         RelEps=.FALSE.
      ELSEIF (RelEpsT .eq. 1) THEN
         WRITE(15,'(/''* Relative tolerance is used for WC, PH and SInk with value of'',3i5)',advance='no')1/factorRelEps
         RelEps=.TRUE.
      ELSE
         STOP 'Set RelEps properly in control.in'
      ENDIF
      IF (factorRelEps.eq.0) THEN
         STOP 'Control.in: factor for relative criterium may not be zero'
      ENDIF
      READ (10,*)
      READ (10,*)
      READ (10,*) dt,dtMin,dtMax,FacInc,FacDec,dtRoot !MJ09
      READ (10,*)
      READ (10,*)
      READ (10,*) nOuFEM
      IF (nOuFEM.GT.mxoFEM) STOP 'Input from  < control.in >  exceeds parameter "mxoFEM". Program terminated.'
      READ (10,*)
      READ (10,*)
      READ (10,*) (tOuFEM(i),i=1,nOuFEM)
      READ (10,*)
      READ (10,*)
      READ (10,*) nOuRoo
      IF (nOuRoo.GT.mxoRoo) STOP 'Input from  < control.in >  exceeds parameter "mxoRoo". Program terminated.'
      READ (10,*)
      READ (10,*)
      READ (10,*) (tOuRoo(i),i=1,nOuRoo)
      IF (nOuFEM.EQ.0) tMax=tOuRoo(nOuRoo)
      IF (nOuRoo.EQ.0) tMax=tOuFEM(nOuFEM)
      IF ((nOuFEM.NE.0).AND.(nOuRoo.NE.0)) tMax=MAX(tOuFEM(nOuFEM),tOuRoo(nOuRoo))
      READ(10,*)
      READ(10,*)
      READ(10,*) integ_zprof,dtprof
      IF (integ_zprof.EQ.1) THEN
          profOK=.TRUE.
          IF (dtprof<999) THEN
             nouProf=0
             ini=0
             DO WHILE (ini.LT.tmax)
                nouProf=nouProf+1
                tOuProf(nouProf)=ini+dtprof
                ini=tOuProf(nouProf)
                IF (nouProf.GT.mxProf) THEN
                   print *,'too small time step in z-profiles: only ',mxProf,' profiles will be kept'
                  GOTO 77
                ENDIF
            ENDDO
   77    ENDIF
      ENDIF
      READ(10,*)
      READ(10,*)
      READ(10,*) i_old
      READ(10,*)
      READ(10,*)
      READ(10,*) i_roottyp
      READ(10,*)
      READ(10,*)
      READ(10,*) i_reduc
      READ(10,*)
      READ(10,*)
      READ(10,*) i_noroot
      READ(10,*)
      READ(10,*)
      READ(10,*) i_continu
      IF (i_old .eq. 2) THEN
         old=.FALSE.
      ELSEIF (i_old .eq. 1) THEN
         old=.TRUE.
      ELSE
         STOP 'Set logical OLD way properly in Condroot.in'
      ENDIF
      oldT=old
      IF (i_roottyp .eq. 2) THEN
         rrt=.FALSE.
      ELSEIF (i_roottyp .eq. 1) THEN
         rrt=.TRUE.
      ELSE
         STOP 'Set parameter properly to switch Roottyp on (1) or off (2)'
      ENDIF
      IF (i_reduc .eq. 2) THEN
         ana_aan=.FALSE.
      ELSEIF (i_reduc .eq. 1) THEN
         ana_aan=.TRUE.
      ELSE
         STOP 'Set parameter properly for memory reducement (1) or off (2)'
      ENDIF
      IF (i_noroot .eq. 2) THEN
         no_root=.FALSE.
      ELSEIF (i_noroot .eq. 1) THEN
         no_root=.TRUE.
      ELSE
         STOP 'Set parameter for a no root simulation on (1) or off (2)'
      ENDIF
      IF ((ana_aan) .and. .not.(old)) STOP 'Cannot use microscopic model and memory reduction simultaneous'
      IF (i_continu .eq. 2) THEN
         continu=.FALSE.
      ELSEIF (i_continu .eq. 1) THEN
         continu=.TRUE.
         WRITE(15,'(/''Roots outside of the domain are automatically introduced on the other side of the domain'')',advance='no')
      ELSE
         STOP 'Set parameter for a continuous domain simulation on (1) or off (2)'
      ENDIF
      IF (continu.AND.(.NOT.old)) STOP 'Domain continuity currently only works with soil Averaging method'	!(Couvreur dec 2009)
      CLOSE (10)
      OPEN (Unit=10,FILE='in/grid.in',STATUS='OLD',ERR=20)
      WRITE(15,'(//''++ Soil Geometry ++'')',advance='no')
      WRITE(15,'(/''--------------------'')',advance='no')
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) nPt,nElm
      IF (nPt.GT.maxnod) STOP 'Input from  < grid.in >  exceeds parameter "maxnod". Program terminated.'
      IF (nElm.GT.maxelm) STOP 'Input from  < grid.in >  exceeds parameter "maxelm". Program terminated.'
!      lOrt=.false.
!      if(nBand.gt.10.or.nPt.gt.300) lOrt=.true.
      lOrt=.true.
      READ (10,*)
      READ (10,*)
      READ (10,*) nx,ny,nz,nex,ney,nez,dxGrid,dyGrid,dzGrid
      nel=nex*ney
      nl=nx*ny
      CALL IniGrid !initialize xgrid, ygrid,zgrid and scaling factors
      !CALL IniMatData! initialize A matrices
! elements:
      READ (10,*)
      READ (10,*)
      READ (10,*)
      DO 2 i=1,nElm
         READ(10,*) iDum,(elmnod(i,k),k=1,6),A11,A22,A33,A12,A13,A23,C11,C22,C33
         ConAxx(i)=C11*A11*A11+C22*A12*A12+C33*A13*A13
         ConAyy(i)=C11*A12*A12+C22*A22*A22+C33*A23*A23
         ConAzz(i)=C11*A13*A13+C22*A23*A23+C33*A33*A33
         ConAxy(i)=C11*A11*A12+C22*A12*A22+C33*A13*A23
         ConAxz(i)=C11*A11*A13+C22*A12*A23+C33*A13*A33
         ConAyz(i)=C11*A12*A13+C22*A22*A23+C33*A23*A33
    2 CONTINUE
      WRITE(15,'(/''* Number of elements in x, y and z directions:'',3i5)',advance='no')nex,ney,nez
      WRITE(15,'(/''* Number of nodes in x, y and z directions:'',3i5)',advance='no')nx,ny,nz
      CLOSE (10)
! nodes:
      OPEN (Unit=10,FILE='in/nodes.in',STATUS='OLD',ERR=30)
    3 READ (10,'(A5)') text
      IF (text.NE.'Node#') GOTO 3
      DO 1 i=1,nPt
         READ (Unit=10,iostat=err) iDum,MatNum(i),xGrid(i),yGrid(i), zGrid(i),hOld(i),Conc(i),Axy(i),Bxy(i),Dxy(i),Exy(i)
         IF (err.ne.0) THEN !Javaux, 2005 in case there are no scaling factors
	     READ (10,*) iDum,MatNum(i),xGrid(i),yGrid(i),zGrid(i), hOld(i),Conc(i)
            Axy(i)=1.
	     Bxy(i)=1.
	     Dxy(i)=1.
	     Exy(i)=1.
         ENDIF
         Kode(i)=0
         Q(i)=0.0_dp
         Qc(i)=0.0_dp
         hNew(i)=hOld(i)
         hTemp(i)=hOld(i)
         nMat=MAX(nMat,MatNum(i))
         xmin=MIN(xmin,xGrid(i))
         ymin=MIN(ymin,yGrid(i))
         zmin=MIN(zmin,zGrid(i))
         xmax=MAX(xmax,xGrid(i))
         ymax=MAX(ymax,yGrid(i))
         zmax=MAX(zmax,zGrid(i))
     1 CONTINUE
       DO i=1,nx
         xCol(i)=xmin+(i-1)*dxGrid
       ENDDO
       DO i=1,ny
         yCol(i)=ymin+(i-1)*dyGrid
       ENDDO
      CLOSE (10)
      IF (nmat.GT.maxmat) STOP 'Input from nodes file exceeds parameter "maxmat". Program terminated.'
      DO 301 iP=1,nPt
         ListNE(iP)=0
301   CONTINUE
      nBand=1
      DO 303 iE=1,nElm
         DO 302 iSE=1,3
          IF (mod(iE,2).eq.0) THEN
             iL(1)=3+iSE
             iL(2)=4+iSE-iSE/3*6
             iL(3)=5+iSE-iSE/2*6
             iL(4)=  iSE
          ELSEIF (mod(iE,2).eq.1) THEN
             iL(1)=3+iSE
             iL(2)=5+iSE-iSE/2*6
             iL(3)=4+iSE-iSE/3*6
             iL(4)=  iSE
          ENDIF
          i=elmnod(iE,iL(1))
          j=elmnod(iE,iL(2))
          k=elmnod(iE,iL(3))
          l=elmnod(iE,iL(4))
            ListNE(i)=ListNE(i)+1
            ListNE(j)=ListNE(j)+1
            ListNE(k)=ListNE(k)+1
            ListNE(l)=ListNE(l)+1
          if(abs(i-j).gt.nBand) nBand=abs(i-j)
          if(abs(i-k).gt.nBand) nBand=abs(i-k)
          if(abs(i-l).gt.nBand) nBand=abs(i-l)
          if(abs(j-k).gt.nBand) nBand=abs(j-k)
          if(abs(j-l).gt.nBand) nBand=abs(j-l)
          if(abs(k-l).gt.nBand) nBand=abs(k-l)
302      CONTINUE
303   CONTINUE
      nBand=nBand+1
      OPEN (Unit=10,FILE='in/bc.in',STATUS='OLD',ERR=40)
      WRITE(15,'(//''++ Soil Boundary Conditions ++'')',advance='no')
      WRITE(15,'(/''------------------------------'')',advance='no')
! (t)-BC [L^3/T]  (Kode(i)=-1):
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) nQ
      IF (nQ.GT.maxbdr) STOP 'Input from  < bc.in >  exceeds parameter "maxbdr". Program terminated.'
      READ (10,*)
      READ (10,*)
      DO i=1,nQ
         iBCPt(i)=i
      ENDDO
      DO 4 i=1,nQ
         Kode(iBCPt(i))=-1
    4 CONTINUE
      READ (10,*)
      READ (10,*)
      READ (10,*) nQbcCh
      IF (nQbcCh.GT.mxBcCh) STOP 'Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.'
      READ (10,*)
      READ (10,*) homogene				!(Couvreur jan 2010)
      READ (10,*)
      IF (homogene==1) THEN
         READ (10,*) (tQbcCh(i),Qbcrec(i),i=1,nQbcCh)
         Qbc(:,1)=Qbcrec
      ELSE
         DO j=1,nQbcCh
            READ (10,*) (Qbcrec(i),i=1,nQ+1)!on peut enlever une dimension à qbctot!!!
            tQbcCh(j)=Qbcrec(1)
            Qbc(j,1:nQ)=Qbcrec(2:nQ+1)
         ENDDO
      ENDIF
      IF (nQ==0) THEN
         WRITE(15,'(/''* No flux B.C.'')',advance='no')
      ELSE
         WRITE(15,'(/''* Nodes with flux B.C.:'',i5)',advance='no') nQ
         WRITE(15,'(/''* Number of time imposed flux B.C.:'',i5)',advance='no')nQbcCh
      ENDIF
! Irrigator [L^3/T]  (Kode(i)=-3):			(Couvreur jan 2010)
      READ (10,*)
      READ (10,*)
      READ (10,*) nI
      IF (nI.GT.maxIrrig) STOP 'Input from  < bc.in >  exceeds parameter "maxIrrig". Program terminated.'
      READ (10,*)
      READ (10,*) (xIrrig(i),yIrrig(i),zIrrig(i),i=1,nI)
      READ (10,*)
      READ (10,*)
      READ (10,*) nIbcCh
      IF (nIbcCh.GT.mxBcCh) STOP 'Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.'
      READ (10,*)
      READ (10,*) (tIbcCh(i),Ibc(i),i=1,nIbcCh)
      DO i=1,nI
         IF ((xIrrig(i).GT.xmax).OR.(xIrrig(i).LT.xmin)) STOP 'Irrigator out of soil domain (X). Program terminated.'
         IF ((yIrrig(i).GT.ymax).OR.(yIrrig(i).LT.ymin)) STOP 'Irrigator out of soil domain (Y). Program terminated.'
         IF ((zIrrig(i).GT.zmax).OR.(zIrrig(i).LT.zmin)) STOP 'Irrigator out of soil domain (Z). Program terminated.'
         ixmin=1+((xIrrig(i)-xmin)/dxgrid)
         ixymin=ixmin+((yIrrig(i)-ymin)/dygrid)*nx
         imin=ixymin+(abs(zIrrig(i)-zmax)/dzgrid)*nl
         IF (ixmin.NE.floor(ixmin)) STOP 'Irrigator position not equal to a node position (X). Program terminated.'
         IF (ixymin.NE.floor(ixymin)) STOP 'Irrigator position not equal to a node position (Y). Program terminated.'
         IF (imin.NE.floor(imin)) STOP 'Irrigator position not equal to a node position (Z). Program terminated.'
         iBCPt(nQ+i)=imin
      ENDDO
      DO i=1,nI
         Kode(iBCPt(nQ+i))=-3
      ENDDO
! h(t)-BC [L]  (Kode(i)=+1):
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) nh
      IF (nh.GT.maxbdr) STOP 'Input from  < bc.in >  exceeds parameter "maxbdr". Program terminated.'
      READ (10,*)
      READ (10,*)
      DO i=1,nh
         iBCPt(nQ+nI+i)=i
      ENDDO
      DO 5 i=1,nh
       Kode(iBCPt(nQ+nI+i))=+1
    5 CONTINUE
      READ (10,*)
      READ (10,*)
      READ (10,*) nhbcCh
      IF (nhbcCh.GT.mxBcCh) STOP 'Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.'
      READ (10,*)
      READ (10,*) (thbcCh(i),hbc(i),i=1,nhbcCh)
      IF (nh==0) THEN
         WRITE(15,'(/''* No head B.C.'')',advance='no')
      ELSE
         WRITE(15,'(/''* Nodes with head B.C.:'',i5)',advance='no')nh
         WRITE(15,'(/''* Number of time imposed head B.C.:'',i5)',advance='no')nhbcCh
      ENDIF
! free drainage-BC  (Kode(i)=-2):
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) nFrdr,nLast
      IF (nFrdr.GT.maxbdr) STOP 'Input from  < bc.in >  exceeds parameter "maxbdr". Program terminated.'
      READ (10,*)
      READ (10,*)
      ii=nLast-nFrdr+1
      DO i=1,nFrdr
         iBCPT(nQ+nI+nh+i)=ii
         ii=ii+1
      ENDDO
      IF (nFrdr==0) THEN
         WRITE(15,'(/''* No free drainage nodes.'')',advance='no')
      ELSE
         WRITE(15,'(/''* Bottom B.C. is free drainage'')',advance='no')
      ENDIF
      DO 6 i=1,nFrdr
         Kode(iBCPt(nQ+nI+nh+i))=-2
    6 CONTINUE
      nBCPts=nQ+nI+nh+nFrdr
      READ (10,*)
      READ (10,*)
      READ (10,*)
      do i=1,nBCPts
         KodCB(i)=-1
      enddo
      READ (10,*)
      READ (10,*)
      READ (10,*)      nCBnd1
      IF (nCBnd1.GT.mxBcCh) STOP 'Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.'
      READ (10,*)
      READ (10,*) (tCBnd1(i),CBnd1(i),i=1,nCBnd1)
      READ (10,*)
      READ (10,*)
      READ (10,*)      nCBnd2
      IF (nCBnd2.GT.mxBcCh) STOP 'Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.'
      READ (10,*)
      READ (10,*) (tCBnd2(i),CBnd2(i),i=1,nCBnd2)
      READ (10,*)
      READ (10,*)
      READ (10,*)      tPulse
      IF ((nCBnd1==0).AND.(nCBnd2==0)) THEN
         WRITE(15,'(/''* No solute transport BCs'')',advance='no')
      ELSEIF (nCBnd1.NE.0) THEN
         WRITE(15,'(/''* 1st time dependent solute transport B.C.: '',i5)',advance='no')nCBnd1
      ELSEIF (nCBnd2.NE.0) THEN
         WRITE(15,'(/''* 2nd time dependent solute transport B.C.: '',i5)',advance='no')nCBnd2
      ENDIF
      CLOSE (10)
!NEW JAVAUX
      WRITE (*,'('' ... looking for file  <Probes.in> ...'')')
      OPEN (Unit=10,FILE='in/Probes.in',STATUS='OLD',ERR=50)
      ObsOK=.TRUE.
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) npr,dtProbe
      WRITE(15,'(/''* Number of sampling probes:'',i5)')npr
      READ (10,*)
      READ (10,*)
      DO i=1,npr
         READ (10,*) Pr(i),Pt(i),CrP(i)
      ENDDO
      READ (10,*)
      READ (10,*)
      READ (10,*) (VarP(i),i=1,npr)
      READ (10,*)
      READ (10,*)
      READ (10,*) (distrP(i),i=1,npr)
      firstOK=.TRUE.
      DO i=1,npr
         IF (Pt(i)==4) THEN !user define node

            IF (firstOK) THEN
               READ (10,*)
               READ (10,*)
               firstOK=.FAlSE.
            ENDIF
            NodebyPr(i)=Crp(i)
            READ (10,*) (node_temp(i2),i2=1,CrP(i)+1)
	     DO i2=1,CrP(i)
	        NodePr(i,i2)=node_temp(i2+1)
	     ENDDO
         ENDIF
      ENDDO

!time
     IF (dtprobe.NE.999) THEN
          nouProbe=0
          ini=0
          DO WHILE (ini.LT.tmax)
             nouProbe=nouProbe+1
             tOuProbe(nouProbe)=ini+dtprobe
             ini=tOuProbe(nouProbe)
             IF (nouProbe.GT.mxProf) THEN
                print *,'too small time step in probe.in: only ',mxProf,' profiles will be kept'
               GOTO 78
             ENDIF
         ENDDO
   78 ENDIF
     CLOSE (10)

!run first calculations
      CALL NodebyProbe

      RETURN
   10 STOP 'File  < control.in >  not found -- program terminated.'
   20 STOP 'File  < grid.in >  not found -- program terminated.'
   30 STOP 'Input file for nodal values not found -- program terminated.'
   40 STOP 'File  < bc.in >  not found -- program terminated.'
   50  WRITE (*,'(/'' File  <Probes.in>  not found: no defined observation locations'')')
      END
!******************************************************************************
      SUBROUTINE DouIn(lDou)
      USE Typedef
      USE DoussanMat, only : ageLr,LrRoot,nLr,ageKh,Khroot,nKh,hx_min,ave,eqDis,stresfun,stresval1,stresval2
      IMPLICIT NONE
      LOGICAL,intent(out):: lDou
      INTEGER(sp) :: i,noRWU,i_ave,i_eqDis
      lDou=.FALSE.
!open simulation summary file
      OPEN(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
      WRITE(15,'(//''++ Hydraulic Data for roots ++'')',advance='no')
      WRITE(15,'(/''-----------------------------'')',advance='no')
      WRITE (*,'('' ... looking for file  <CondRoot.in> ...'')')
      OPEN (Unit=10,FILE='in/CondRoot.in',STATUS='OLD',ERR=3001)
!if CondRoot exists
      lDou = .TRUE.
      WRITE(15,'(/''* The Doussan algorithm will be used for the root water uptake (lDou=TRUE)'')',advance='no')
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) (nLr(i),i=1,3)
      READ (10,*)
      READ (10,*) (ageLr(1,i),LrRoot(1,i),i=1,nLr(1))  !main axes
      READ (10,*) (ageLr(2,i),LrRoot(2,i),i=1,nLr(2))  !secondary axes
      READ (10,*) (ageLr(3,i),LrRoot(3,i),i=1,nLr(3))  !tertiary axes
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) (nKh(i),i=1, 3)
      READ (10,*)
      READ (10,*) (ageKh(1,i),KhRoot(1,i),i=1,nKh(1))
      READ (10,*) (ageKh(2,i),KhRoot(2,i),i=1,nKh(2))
      READ (10,*) (ageKh(3,i),KhRoot(3,i),i=1,nKh(3))
!averaging method
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) i_ave
      READ (10,*)
      READ (10,*) i_eqDis
      IF (i_ave .eq. 2) THEN
         ave=.FALSE.
      ELSEIF (i_ave .eq. 1) THEN
         ave=.TRUE.
         write(15,'(/''* Averaging method: treated as one root (ave= TRUE)'')',advance='no')
      ELSE
         STOP 'Set AVERAGE properly in Condroot.in'
      ENDIF
      IF (i_eqDis .eq. 2) THEN
         eqDis=.FALSE.
          IF (i_ave .eq. 2) THEN
             write(15,'(/''* No averaging method will be used below voxel scale (ave= FALSE and eqDis= FALSE)'')',advance='no')
          ENDIF
      ELSEIF (i_eqDis .eq. 1) THEN
         eqDis=.TRUE.
         write(15,'(/''* Averaging method: use equal distance function (eqDis= TRUE)'')',advance='no')
      ELSE
         STOP 'Set eqDis parameter properly in Condroot.in'
      ENDIF
!stress functions
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)stresfun
      READ (10,*,ERR=20)
      IF (stresfun.EQ.1) THEN
          READ (10,*,ERR=20) hx_min
          hx_min=-ABS(hx_min)
          WRITE(15,'(/''* The limiting stress value is  '',1pE11.4)',advance='no') hx_min
      ELSE
          READ (10,*,ERR=20) stresval1,stresval2
          IF (stresfun.EQ.2) THEN
             WRITE(15,'(/''* A linear stress function will be used with values = '',1pE11.4,1pE11.4)',advance='no') stresval1,stresval2
          ELSE
             WRITE(15,'(/''* The Tuzet stress function will be used with values  = '',1pE11.4,1pE11.4)',advance='no') stresval1,stresval2
          ENDIF
      ENDIF
      if ((ave) .and. (eqDis)) STOP 'set one root and equidistant approach correctly; may not be both switched on'
      CLOSE (15)
      RETURN
      CLOSE (10)
3001  WRITE (15,'(/'' File  <CondRoot.in>  not found --'')')
      CLOSE (15)
20    STOP ' Root pressure head limitation not found in <CondRoot.in>. -- program terminated.'
      END
!*******************************************************************************
      SUBROUTINE ChemIn(lChem)
      USe Typedef
      Use SolData
      IMPLICIT NONE
      LOGICAL lChem
      INTEGER(sp) :: i,k
      lChem=.TRUE.
      OPEN(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
      WRITE(15,'(//''++ Solute Transport Information ++'')',advance='no')
      WRITE(15,'(/''----------------------------------'')')
      OPEN(Unit=10,File='in/chem.in',STATUS='OLD',ERR=10)
      READ (10,*)
      READ (10,*)
      READ (10,*)
      NLevel=1
      READ (10,*) epsi
      IF (epsi.LT.0.999_dp) NLevel=2
      READ (10,*)
      READ (10,*) PeCr
      READ (10,*)
      READ (10,*)
      DO 1 k=1,nMat
         READ (10,*,ERR=20) (ChPar(i,k),i=1,9)
    1 CONTINUE
      CLOSE (10)
      WRITE (15,'(/'' File < chem.in > found --'')',advance='no')
      WRITE (15,'('' -- simulation will include solute transport.'')',advance='no')
      RETURN
   10 WRITE (15,'('' File  < chem.in >  not found --'')',advance='no')
      WRITE (15,'('' -- simulation will not include solute transport. (lChem=False)'')',advance='no')
      lChem=.FALSE.
      RETURN
   20 STOP 'Insufficient data in  < chem.in >  -- program terminated.'
      END
!****************************************************************************
      SUBROUTINE SoilIn
      USE Typedef
      USE SolData
      USE WatFun
      IMPLICIT NONE
      REAL(sp) ::alh,h1,hn,htab1,htabn
      INTEGER(sp) ::i,k
      OPEN(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
      WRITE(15,'(//''++ Soil Information ++'')',advance='no')
      WRITE(15,'(/''----------------------'')',advance='no')
      OPEN(Unit=10,FILE='in/soil.in',STATUS='OLD',ERR=10)
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*,ERR=30)nMat,h1,hN
      IF (h1*hN.EQ..0_dp) THEN
         WRITE (*,'('' Input error -- please use non-zero values for hTab1 and hTab2 in  < soil.in >.'')')
         STOP 'Program terminated.'
      ENDIF
      WRITE(15,'(/''* Tabulated values between: '',1pE11.4,'' and '',1pE11.4)',advance='no') h1,hN
      WRITE(15,'(/''* Number of soil layers:'',1i7)',advance='no') nMat
      hTab1=-MIN(ABS(h1),ABS(hN))
      hTabN=-MAX(ABS(h1),ABS(hN))
      READ (10,*)
      READ (10,*)
      DO 1 k=1,nMat
         READ (10,*,ERR=20) (par(i,k),i=1,11)
    1 CONTINUE
      CLOSE (10)
! generate table for interpolation:
      IF (ABS(hTabN-hTab1).GT.0._dp) THEN
         alh1=LOG10(-hTab1)
         dlh=(LOG10(-hTabN)-alh1)/REAL(nTab-1)
         DO 2 i=1,nTab
            alh=alh1+REAL(i-1)*dlh
            hTab(i)=-(10._dp)**alh
    2    CONTINUE
         DO 3 k=1,nMat
            DO 3 i=1,nTab
               ConTab(i,k)=FKP(hTab(i),par(:,k))
               CapTab(i,k)=FCP(hTab(i),par(:,k))
	        TheTab(i,k)=FTh(hTab(i),par(:,k))
    3    CONTINUE
      ELSE
         WRITE (*,'('' Input error --  hTab1  and  hTab2  in  < soil.in>  cannot be equal.'')')
         STOP 'Program terminated.'
      ENDIF
      RETURN
   10 STOP 'File  < soil.in >  not found -- program terminated.'
   20 STOP 'Insufficient data in  < soil.in >  -- program terminated.'
   30 STOP 'Soil input file has been changed (23 Oct 08): now, add the number soil types at line 4 !'
      END
!***************************************************************************
      SUBROUTINE RootIn(naxes,norder,t,mroot,mshoot,LA,sAvg,cAvg,toxi)
      USE Typedef
      USE TempData
      USE RootData
      USE PlntData
      USE GridData
      USE GeoData
      USE ConData
      USE StrData
      USE DoussanMat, only: nplant
      IMPLICIT NONE

      REAL(sp) :: ranang(maxemg),mroot,mshoot,LA,maxZ
      INTEGER(sp) :: ifg(maxemg),naxes,naxtot,norder,iplant,ipl
      INTEGER(sp) :: iestbl,i,j,k,irec,igrow,ifive,timeRT
      REAL(sp) :: linlim,ranfac,cavg,savg,t
      LOGICAL toxi
      CHARACTER infile*12,ch*1

      OPEN(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
      WRITE(15,'(//''++ Root System ++'')',advance='no')
      WRITE(15,'(/''-----------------'')',advance='no')
      OPEN (UNIT=10,FILE='in/root.in',STATUS='OLD',ERR=10)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) naxemg
      IF (naxemg.GT.maxemg) STOP 'Input from  < root.in >  exceeds parameter "maxemg". Program terminated.'
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (tnewax(i),nnewax(i),i=1,naxemg)
      naxtot=0
      DO 1 i=1,naxemg
         ifg(i)=naxtot+1
         naxtot=naxtot+nnewax(i)
    1 CONTINUE
      IF (naxtot.GT.maxaxs) STOP 'Input from  < root.in >  exceeds parameter "maxaxs". Program terminated.'
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (geoaxs(ifg(i)),i=1,naxemg)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (ranang(i),i=1,naxemg)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (nangax(ifg(i)),i=1,naxemg)
      READ (10,*,ERR=20)
      DO 2 i=1,naxemg
       READ (10,*,ERR=20)  (tempax(ifg(i),j),angaxs(ifg(i),j),j=1,nangax(ifg(i)))
       DO 21 k=ifg(i)+1,ifg(i)+nnewax(i)-1
          geoaxs(k)=geoaxs(ifg(i))
          nangax(k)=nangax(ifg(i))
          DO 21 j=1,nangax(k)
             angaxs(k,j)=angaxs(ifg(i),j)
             tempax(k,j)=tempax(ifg(i),j)
   21    CONTINUE
       DO 2 k=ifg(i),ifg(i)+nnewax(i)-1
         ranfac=ran(irnd)
          DO 2 j=1,nangax(k)
          angaxs(k,j)=(angaxs(k,j)+(.5_dp-ranfac)*ranang(i))/180._dp*pi
    2 CONTINUE
  999 READ (10,'(1A1)',ERR=20) ch
      IF (ch.NE.'*') GOTO 999
      READ (10,*,ERR=20) geolat
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) nanglt
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (templt(i),anglat(i),i=1,nanglt)
      DO 3 j=1,nanglt
         anglat(i)=anglat(i)/180._dp*pi
    3 CONTINUE
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) norder
      IF (norder.GT.maxord) STOP 'Input from  < root.in >  exceeds parameter "maxord". Program terminated.'
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (nVch(i),i=1,norder)
      READ (10,*,ERR=20)
      DO 4 i=1,norder
         READ (10,*,ERR=20) (ageVch(i,j),Vch(i,j),j=1,nVch(i))
    4 CONTINUE
  998 READ (10,'(1A1)',ERR=20) ch
      IF (ch.NE.'*') GOTO 998
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (nMPLch(i),i=1,norder)
      READ (10,*,ERR=20)
      DO 5 i=1,norder
         READ (10,*,ERR=20) (sMPLch(i,j),MPLch(i,j),j=1,nMPLch(i))
    5 CONTINUE
  997 READ (10,'(1A1)',ERR=20) ch
      IF (ch.NE.'*') GOTO 997
      READ (10,*,ERR=20) tmin,topt,tmax
      trange=tmax-tmin
      tmid=(tmin+tmax)/2._dp
      IF (topt.LT.tmid) THEN
         expo=LOG(.5_dp)/LOG((topt-tmin)/trange)
      ELSE
         expo=LOG(.5_dp)/LOG((topt-tmax)/(-trange))
      ENDIF
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) toxi
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      IF (toxi) THEN
         READ (10,*,ERR=20) cmin,coptmi,coptma,cmax
      ELSE
         READ (10,*,ERR=20)
      ENDIF
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) simp
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) refgrd
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (strsen(i),i=1,norder)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (rdmang(i),i=1,norder)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (brlmax(i),i=1,norder)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (brspac(i),i=1,norder-1)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (brnang(i),i=1,norder-1)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (dtbrch(i),i=1,norder-1)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) p50,h50,p1,p2
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) CMm,VMax,xin,SpWgt,fk
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) h0,h1,h2,h3
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20)
      READ (10,*,ERR=20) (nUrf(i),i=1,norder)
      READ (10,*,ERR=20)
      DO 44 i=1,norder
         READ (10,*,ERR=20) (age(i,j),Urf(i,j),j=1,nUrf(i))
 44   CONTINUE
      CLOSE (10)
      h0=-ABS(h0)
      h1=-ABS(h1)
      h2=-ABS(h2)
      h3=-ABS(h3)
      DO 6 i=1,norder-1
         brnang(i)=brnang(i)/180._dp*pi
         rdmang(i)=rdmang(i)/180._dp*pi
    6 CONTINUE
      rdmang(norder)=rdmang(norder)/180._dp*pi
!RootTip
	  no_root_gwth=.FALSE.
	  IF (rrt) THEN
		no_root_gwth=.TRUE.
	    OPEN (UNIT=10,FILE='in/param.txt',STATUS='OLD',ERR=50) !param.txt
		READ (10,*)
              READ (10,*)
		READ (10,*)
		READ (10,*)
		READ (10,*) timeRT
! no multiple roots with RootTyp
              nplant=1
              xplant(1)=0
              yplant(1)=0

              CALL RunRootTip(naxes,timeRT)
		t=timeRT !implicit change of type INT->real

              WRITE(15,'(/''RootTyp was used for generating the root system'')')
              WRITE(15,'(/''Initial Time for simulation is '',1pE11.4)',advance='no') t
	  ELSE
              infile='in/RootSys'
              OPEN (UNIT=10,FILE=infile,STATUS='OLD',ERR=60)
              WRITE(15,'(/''* Root system imported from file '',1A12)',advance='no') infile
	!	WRITE (*, '(///'' Name of Root System Input File: '',$)')
	!	READ (*, '(1A12)') infile
		READ (10,*,ERR=40)
		READ (10,*,ERR=40) t
             WRITE(15,'(/''* Initial Time for simulation is '',1pE11.4)',advance='no') t
		READ (10,*,ERR=40)
		READ (10,*,ERR=40)
              READ (10,*,ERR=40) nplant
              WRITE(15,'(/''* Amount of different root systems : '',I5)',advance='no') nplant
		   READ (10,*,ERR=40)
		   READ (10,*,ERR=40)
                 DO i=1,nplant
                    READ (10,*,ERR=40) iplant,xplant(iplant),yplant(iplant)
                    WRITE(15,'(/''* seed locations for plant '',1I5, '' is '', 2F9.3)',advance='no')iplant,xplant(iplant),yplant(iplant)
                 ENDDO
		READ (10,*,ERR=40)
		READ (10,*,ERR=40)
		READ (10,*,ERR=40) mroot,mshoot,LA
		READ (10,*,ERR=40)
		READ (10,*,ERR=40)
		READ (10,*,ERR=40) sAvg,cAvg
		READ (10,*,ERR=40)
		READ (10,*,ERR=40)
		READ (10,*,ERR=40) naxes
		READ (10,*,ERR=40)
		READ (10,*,ERR=40)
		READ (10,*,ERR=40) nbr
		READ (10,*,ERR=40)
		READ (10,*,ERR=40)
		READ (10,*,ERR=40) nrec
		IF ((nbr.LT.1).OR.(nrec.LT.1))  GOTO 20
		READ (10,*,ERR=40)
		READ (10,*,ERR=40)
		READ (10,*,ERR=40)

    7	READ (10,*,ERR=40) irec,xs(irec),ys(irec),zs(irec),irecpr(irec),ordseg(irec),&
		     ibrseg(irec),seglen(irec),segsur(irec),segmas(irec)
		READ (10,*,ERR=40) timorg(irec)
	    IF (irec.LT.nrec) GOTO 7
          !  maxZ = maxval(zs(1:nrec))
         !  zs(1:nrec)=zs(1:nrec)-maxZ
!check first soil node locations
           !IF ((zgrid(npt).LT.0).AND.(zgrid(1).EQ.0)) THEN !soils surface=0 and bottom= negative value
           !    zgrid=zgrid-minval(zgrid)
           !    print *,'z direction has been automatically changed in the code to be at the soil surface equal to',(zgrid(1))
          ! ENDIF 
		IF ((zs(1).GT.zgrid(1)).OR.(zs(1).LT.zgrid(npt))) THEN
			WRITE (*,'(//'' Inconsistency -- Position of first root segment not within the spatial domain'')')
			WRITE (*,'('' as defined in  < nodes.in >.''/)')
                    	STOP 'Program terminated.'
		ENDIF
		IF (.not.(continu)) THEN			!(Couvreur dec 2009)
                 DO ipl=1,nplant
		      CALL CheckSize(ipl)
                 ENDDO
		ENDIF
		READ (10,*,ERR=40)
		READ (10,*,ERR=40)
		READ (10,*,ERR=40) ngrow
		IF (ngrow.LT.1) GOTO 20
		READ (10,*,ERR=40)
		READ (10,*,ERR=40)
		READ (10,*,ERR=40)
		READ (10,*,ERR=40)
    8	READ (10,*,ERR=40) igrow,xg(igrow),yg(igrow),zg(igrow),irecsg(igrow),ordgrw(igrow),ibrgrw(igrow),&
	&brlgth(igrow), iaxis(igrow)
		READ (10,*,ERR=40) ovrlen(igrow),nestbl(igrow)
		IF (nestbl(igrow).GT.0) THEN
			 ifive=0
   81		linlim=MIN(ifive+5,nestbl(igrow))
			READ (10,*,ERR=40) (timest(igrow,iestbl),iestbl=ifive+1,linlim)
			ifive=ifive+5
			IF (linlim.LT.nestbl(igrow)) GOTO 81
		ENDIF

		IF (igrow.LT.ngrow) GOTO 8
                zg(1:ngrow)=zg(1:ngrow)-maxZ
            !    write(*,*)'zg=',zg(1:ngrow)
		CLOSE (10)
		WRITE (*,'(///'' Simulation starts at time = '',F8.3,''.'')') t
		WRITE (15,'(/''The root system consists now of '',I5,'' branch(es)[including '',I3,'' axis(es)]'')')nbr,naxes
		WRITE (15,'('' or of a total of '',I5,'' segment(s) and '',I5,'' growing branch tip(s).'')')nrec,ngrow
		WRITE (*,'(/'' Total root mass is '',F9.3,'', total shoot mass '',F9.3,''.'')')mroot,mshoot
		WRITE (*,'('' Leaf area is '',F9.3,''.'')')LA
		WRITE (*,'(///70X,''<ENTER>''/)')
      ENDIF
      RETURN
   10 STOP 'File  < root.in >  not found -- program terminated.'
   20 STOP 'Data inconsistency in  < root.in >  -- program terminated.'
   30 STOP 'Root input file does not exist -- program terminated.'
   40 STOP 'Data error in root system file  -- program terminated.'
   50 STOP 'File param.txt not found --program terminated'
   60 STOP 'File RootSys not found --program terminated'
      END
!**************************************************************************************************
      SUBROUTINE PlntIn(level,lChem,toxi,lDou,t)
      USE PlntData
      USE tmctrl
      IMPLICIT NONE
      INTEGER (sp) :: level,i,err,typeBCrt,funBC
      REAL(sp) :: tBCrt1,tBCrt2,value,t
      LOGICAL lChem,toxi,lDou

!by default level=3 exept if lDou is on (then level=4)
      level=3
!read BCRoot.in
IF (lDou) THEN !CondRoot exists
! DOUSSAN with mix BC INPUTS-level 4 (Javaux 06)
! file with boundary conditions given for root in therms of PH or of flux
         WRITE (*,'('' ... looking for file  <BCroot.in> ...'')')
         OPEN (Unit=10,FILE='in/BCroot.in',STATUS='OLD',ERR=3001)
         READ (10,*)
         READ (10,*)
         READ (10,*)
         READ (10,*)
         READ (10,*) funBC
         READ (10,*)
         READ (10,*) nBCr
         IF (nBCr.GT.mxBcCh) THEN
              print *,' the number of boundary conditions for root exceeded the maximum. Simulation stopped.'
              STOP
         ENDIF
         READ (10,*)
         IF (funBC.LE.2) THEN !free format or constant value
             READ (10,*) (tBCr(i),typeBCr(i),BCroot(i),i=1,nBCr)
         ELSE
            READ (10,*) tBCrt1,typeBCrt,value
            tBCrt2=tMax
            DO i=1,nBCr
               tBCr(i)=(tBCrt2-tBCrt1)/(nBCr-1)*(i-1)+tBCrt1
               typeBCr(i)=typeBCrt
               IF (funBC.eq.3) THEN  !positive sin function
                  BCroot(i) = value*((sin(2*pi*(tBCr(i)-0.25)))/2+0.5)
               ELSEIF (funBC.eq.4) THEN !'jump' function
                  BCroot(i) = value*pi*sin((tBCr(i)-0.25)*2*pi)
	           IF (BCroot(i).LT.0) BCroot(i)=0.0
	        ELSE
                  WRITE (*,'('' ... ERROR in input file BCroot.in ...'')')
                  STOP
               ENDIF
            ENDDO
         ENDIF

         IF (tBCR(1).LT.t) THEN
            WRITE(*,'(/'' initial time for root BC lower than initial simulation time'')')
            STOP
         ENDIF
         WRITE (*,'(/'' BC-data for root found.'')')
         WRITE (*,'(/'' level 4 simulations will be performed'' )')
         level=4
         WRITE (*,'(///70X,''<ENTER>''/)')
 !        READ (*,*)
         RETURN
         CLOSE (10)
3001     WRITE (*,'(/'' File  <BCroot.in>  not found --'')')
      ENDIF
!---------------------------------------
      OPEN (Unit=10,FILE='in/plant.in',STATUS='OLD',ERR=10)
!    --- all plant input is expressed as piecewise linear functions ---
! no-stress potential transpiration rate per leaf area as a function of time:
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) ntTpLA
      IF (ntTpLA.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (tTpLA(i),TpLAc(i),i=1,ntTpLA)
! no-stress water use efficiency (dry mass gained per water mass transpired)
! as f(time):
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) ntW
      IF (ntW.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (tW(i),Wc(i),i=1,ntW)
! no-stress root/shoot ratio as a f(time):
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) ntRSR
      IF (ntRSR.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (tRSR(i),RSRc(i),i=1,ntRSR)
! transpiration reduction factor as a function of relative stress due soil strength:
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) nfTpLA
      IF (nfTpLA.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (sfTpLA(i),fTpLAc(i),i=1,nfTpLA)
! water use efficiency factor as a function of relative stress due soil strength:
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) nsfW
      IF (nsfW.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (sfW(i),fWc(i),i=1,nsfW)
! root/shoot-ratio factor as a function of relative stress due soil strength:
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) nsfRSR
      IF (nsfRSR.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (sfRSR(i),fRSRc(i),i=1,nsfRSR)
! relative stress as a function of soil strength:
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) ns
      IF (ns.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (sc(i),rsc(i),i=1,ns)
! transpiration reduction factor as a function of relative stress due solute concentration:
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) ncTpLA
      IF (ncTpLA.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (scTpLA(i),cTpLAc(i),i=1,ncTpLA)
! water use efficiency factor as a function of relative stress due solute concentration:
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) nscW
      IF (nscW.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (scW(i),cWc(i),i=1,nscW)
! root/shoot-ratio factor as a function of relative stress due solute concentration:
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) nscRSR
      IF (nscRSR.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (scRSR(i),cRSRc(i),i=1,nscRSR)
! relative stress as a function of solute concentration:
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) ncnc
      IF (ncnc.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (cncp(i),rscnc(i),i=1,ncnc)
! leaf area per dry shoot mass as a function of time:
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) ntLA
      IF (ntLA.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (tLA(i),LAc(i),i=1,ntLA)
      CLOSE (10)
      IF (h50.LT.9.E+20_dp) THEN
         IF (p50.LT.9.E+20_dp) THEN
            WRITE (*,'(/'' The van Genuchten (1978) expression will be used for the water'')')
            WRITE (*,'(''     extraction function. '')')
         ELSE
            WRITE (*,'(/'' The van Genuchten (1978) expression will be used for the water'')')
            WRITE (*,'(''     extraction function, but osmotic potential effects will be'')')
            WRITE (*,'(''     neglected. '')')
         ENDIF
      ELSE
         WRITE (*,'(//'' The Feddes et al. (1978) expression will be used for the water'')')
         WRITE (*,'(''     extraction function. '')')
      ENDIF
      IF (toxi) THEN
         WRITE (*,'(//'' Nutrient deficiency and/or ion toxicity effects will be included'')')
         WRITE (*,'(''     in the simulation. '')')
      ELSE
         WRITE (*,'(//'' Nutrient deficiency and/or ion toxicity effects will not be '')')
         WRITE (*,'(''     included in the simulation. '')')
      END IF
      WRITE (*,'(///70X,''<ENTER>''/)')
      RETURN
   10 WRITE (*,'(/'' File  < plant.in >  not found --'')')
      WRITE (*,'('' ... looking for file  < tpot.in > ...'')')
      OPEN (Unit=10,FILE='in/tpot.in',STATUS='OLD',ERR=20)
! potential (e.g., not limited by low soil water potential) transpiration rate
! as f(time):
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) nTpot
      IF (nTpot.GT.mxBcCh) GOTO 99
      READ (10,*)
      READ (10,*) (tTpot(i),Tpotc(i),i=1,nTpot)
      WRITE (*,'(/'' Tpot(time)-data found.'')')
      IF (lChem) THEN
         WRITE (*,'(/'' Interaction between root water - solute uptake / root growth'')')
         WRITE (*,'(''     and soil water content / soil strength will be simulated.'')')
      ELSE
         WRITE (*,'(/'' Interaction between root water uptake / root growth'')')
         WRITE (*,'(''     and soil water content / soil strength will be simulated.'')')
      END IF
      WRITE (*,'(/'' Shoot growth and assimilate distribution will not be considered. '')')
      WRITE (*,'(///70X,''<ENTER>''/)')
  !    READ (*,*)
      level=2
      IF (h50.LT.9.E+20_dp) THEN
         IF (p50.LT.9.E+20_dp) THEN
            WRITE (*,'(/'' The van Genuchten (1978) expression will be used for the water'')')
            WRITE      (*,'(''     extraction function. '')')
         ELSE
            WRITE (*,'(/'' The van Genuchten (1978) expression will be used for the water'')')
            WRITE (*,'(''     extraction function, but osmotic potential effects will be'')')
            WRITE (*,'(''     neglected. '')')
         ENDIF
      ELSE
         WRITE (*,'(//'' The Feddes et al. (1978) expression will be used for the water'')')
         WRITE (*,'(''     extraction function. '')')
      ENDIF
      IF (toxi) THEN
         WRITE (*,'(//'' Nutrient deficiency and/or ion toxicity effects will be included'')')
         WRITE (*,'(''     in the simulation. '')')
      ELSE
         WRITE (*,'(//'' Nutrient deficiency and/or ion toxicity effects will not be '')')
         WRITE (*,'(''     included in the simulation. '')')
      END IF
      WRITE (*,'(///70X,''<ENTER>''/)')
!      READ (*,*)
      RETURN
   20 WRITE (*,'(/'' File  < tpot.in >  not found --'')')
      WRITE (*,'('' -- will simulate root growth as affected by soil strength only.'')')
      WRITE (*,'(/'' Water and solute uptake, shoot growth and assimilate      '')')
      WRITE (*,'(''     distribution will not be considered.'')')
      IF (toxi) THEN
         IF (lChem) THEN
            WRITE (*,'(//'' Nutrient deficiency and/or ion toxicity effects will be included'')')
            WRITE (*,'(''     in the simulation. '')')
         ELSE
            WRITE (*,'(//'' Nutrient deficiency and/or ion toxicity effects will be included'')')
            WRITE (*,'(''     in the simulation using initial concentration.'')')
         END IF
      ELSE
         WRITE (*,'(//'' Nutrient deficiency and/or ion toxicity effects will not be '')')
         WRITE (*,'(''     included in the simulation. '')')
      END IF
      WRITE (*,'(///70X,''<ENTER>''/)')
 !     READ (*,*)
      level=1
      RETURN
   99 STOP 'Input from  < plant.in >  exceeds parameter "mxBcCh". Program terminated.'
      END
!***************************************************************************************
      SUBROUTINE TempIn(temp)
      USE TempData
      Implicit none
      LOGICAL temp
      INTEGER(sp) :: iz
      OPEN (Unit=10,FILE='in/temp.in',STATUS='OLD',ERR=10)
! soil temperature as a function of time and depth
! (piecewise linear function):
      READ (10,*)
      READ (10,*)
      READ (10,*)
      READ (10,*) nz_temp
      IF (nz_temp.GT.mxdpth) STOP 'Input from  < temp.in >  exceeds parameter"mxdpth". Program terminated.'
      READ (10,*)
      READ (10,*)
      READ (10,*) (depth(iz),iz=1,nz_temp)
      IF (ABS(depth(1)).GT.1.E-30_dp) THEN
       WRITE (*,'(//'' Inconsistency - First depth value in  < temp.in >  corresponds to soil surface'')')
       WRITE (*,'('' and must be zero.''/)')
       STOP 'Program terminated.'
      ENDIF
      READ (10,*)
      READ (10,*)
      nt=0
  100 nt=nt+1
      IF (nt.GT.mxtime) STOP 'Input from  < temp.in >  exceeds parameter"mxtime". Program terminated.'
      READ (10,*,END=20) time(nt),(temtim(nt,iz),iz=1,nz_temp)
      GOTO 100
   20 CLOSE (10)
      nt=nt-1
      WRITE (*,'(/'' Temperature data found.'')')
      WRITE (*,'('' Simulation will include temperature influence.'')')
      WRITE (*,'(///70X,''<ENTER>''/)')
  !    READ (*,*)
      temp=.TRUE.
      RETURN
   10 WRITE (*,'(/'' File  < temp.in >  not found --'')')
      WRITE (*,'('' -- will simulate root growth assuming optimal soil temperature conditions.'')')
      WRITE (*,'(///70X,''<ENTER>''/)')
  !    READ (*,*)
      temp=.FALSE.
      RETURN
      END
!********************************************************************
   SUBROUTINE NodebyProbe
!define the number of nodes located in one observation plane/probe window
   USE Typedef
   USE ObsData
   USE GridData
   USE DomData
   IMPLICIT NONE 
   INTEGER(sp) :: ip,delta,in,imin,i,iz,ix,iy

! observation planes
     DO ip=1,npr
        IF (Pt(ip)==3) THEN !horiz. plane perp. to Z
	     nodebyPr(ip)=nx*ny
            delta=int((abs(CrP(ip)-zmax))/dzgrid)

            imin=nodebyPr(ip)*delta+1
            in=imin
            i=0
            DO i=1,nodebyPr(ip)
               NodePr(ip,i)=in
               in=in+1
            ENDDO
        ELSEIF (Pt(ip)==2) THEN !vert. plane perp. to Y
	    nodebyPr(ip)=nx*nz
           delta=int((abs(CrP(ip)-ymin))/dygrid)
           imin=delta*nx+1
           i=0
           in=imin
           DO iz=1,nz
              DO ix=1,nx
			   i=i+1
		       NodePr(ip,i)=in
		       in=in+1
               ENDDO
			  in=in+(nx)*(ny-1)
           ENDDO
         ELSEIF (Pt(ip)==1) THEN ! vert plane perp to X
           nodebyPr(ip)=ny*nz
           delta=int((abs(CrP(ip)-xmin))/dxgrid)
           imin=delta+1
           i=0
           in=imin
           DO iz=1,nz
              DO iy=1,ny
                 i=i+1
                 NodePr(ip,i)=in
                 in=in+nx
              ENDDO
           ENDDO
        ELSEIF (Pt(ip).NE.4) THEN
          write(*,*)'Wrong plane direction in Probes.in'
       	STOP
        ENDIF
		 write(*,*)'plane',ip,'has ',nodebyPr(ip),'nodes.'
		 write(*,*) Pt(ip),ip,nodebyPr(ip),delta,imin,in
         IF(nodebyPr(ip)>1000) THEN
		    write(*,*) 'too many elements in one observation plane (>1000)'
		    STOP
		 ENDIF
	   ENDDO
	   RETURN
	   stop
	 END SUBROUTINE NodebyProbe
! ==============================================================================
! Source file ENVIRONMENTAL FUNCTIONS ||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
! current spatial soil strength distribution:
      SUBROUTINE Solstr(h,MatNum)
      USE GridData
      USE SolData
      USE WatFun
      USE StrData
      IMPLICIT NONE
      INTEGER(sp) :: MatNum(maxnod),i,m
      REAL(sp) ::  h(maxnod)
! s = soil strength
! par(11,M) = sMax = max. soil strength for each soil [=f(rhoB,sand content)]
      DO 1 i=1,nPt
         M=MatNum(i)
!!!!* assign current soil strength values to each node:
      s(i)=(1._dp-FThNrm(h(i),par(:,M)))*(1._dp-FThNrm(h(i),par(:,M)))*(1._dp-FThNrm(h(i),par(:,M)))*par(11,M)
!!!!* current local impedance factor due to soil strength:
         IF (s(i).GT.simp) THEN
            imps(i)=0.0_dp
         ELSE
            imps(i)=1._dp-s(i)/simp
         ENDIF
    1 CONTINUE
      RETURN
      END
!***********************************************************
! current spatial temperature distribution:
      SUBROUTINE Temper(t)
      USE TempData
      USE GridData
      implicit none
      REAL(sp) :: val1,val2,depfac,t,timfac
      INTEGER(sp) :: i,it,it1,it2,iz1,iz,iz2
!!!!* find the current time interval:
      IF (t.GE.time(nt)) THEN
         it1=nt
      ELSE
         it=nt
    1    it=it-1
         IF ((t.LT.time(it)).AND.(it.GT.1)) GOTO 1
         it1=it
         it2=it+1
         timfac=(t-time(it1))/(time(it2)-time(it1))
      ENDIF
      DO 10 i=1,nPt
!!!!* assign current temperature values to each node --
!!!!* find the appropriate depth interval:
         IF (zGrid(i).LE.zGrid(1)-depth(nz_temp)) THEN
            iz1=nz_temp
         ELSE
            iz=nz_temp
   11       iz=iz-1
            IF ((zGrid(i).GT.zGrid(1)-depth(iz)).AND.(iz.GT.1)) GOTO 11
            iz1=iz
            iz2=iz+1
            depfac=(zGrid(i)-(zGrid(1)-depth(iz1)))/((zGrid(1)-depth(iz2))-(zGrid(1)-depth(iz1)))
         ENDIF
!!!!* interpolate along depth- and time-coordinates:
         IF (iz1.EQ.nz_temp) THEN
            IF (it1.EQ.nt) THEN
               tem(i)=temtim(nt,nz_temp)
            ELSE
               tem(i)=timfac*(temtim(it2,nz_temp)-temtim(it1,nz_temp))+temtim(it1,nz_temp)
            ENDIF
         ELSE
            IF (it1.EQ.nt) THEN
               tem(i)=depfac*(temtim(nt,iz2)-temtim(nt,iz1))+temtim(nt,iz1)
            ELSE
               val1=depfac*(temtim(it1,iz2)-temtim(it1,iz1))+temtim(it1,iz1)
               val2=depfac*(temtim(it2,iz2)-temtim(it2,iz1))+temtim(it2,iz1)
               tem(i)=timfac*(val2-val1)+val1
            ENDIF
         ENDIF
!* current local impedance factor due to soil temperature:
         IF ((tem(i).GE.tmax).OR.(tem(i).LE.tmin)) THEN
            impt(i)=0.0_dp
         ELSE
            IF (topt.LT.tmid) THEN
               impt(i)=SIN(pi*((tem(i)-tmin)/trange)**expo)
            ELSE
               impt(i)=SIN(pi*((tem(i)-tmax)/(-trange))**expo)
            ENDIF
         ENDIF
   10 CONTINUE
      RETURN
      END
!**********************************************************************
! current nodal solute concentration impedance factor values
      SUBROUTINE ConTox(Conc)
      USE GridData
      USE ConData
      IMPLICIT NONE
      REAL,intent(in) :: Conc(maxnod)
      REAL(sp) :: c0,c1,c2,c3
      INTEGER(sp) :: i
!* current local impedance factor due to soil water solution concentration:
      c0=cmin
      c1=coptmi
      c2=coptma
      c3=cmax
      DO 10 i=1,nPt
         impc(i)=0.0_dp
         IF (Conc(i).GT.c0.AND.Conc(i).LT.c1) impc(i)=(Conc(i)-c0)/(c1-c0)
         IF (Conc(i).GE.c1.AND.Conc(i).LE.c2) impc(i)=1._dp
         IF (Conc(i).GT.c2.AND.Conc(i).LT.c3)impc(i)=(Conc(i)-c3)/(c2-c3)
10    CONTINUE
      RETURN
      END
!***********************************************************************************
! current local soil strength value:
      SUBROUTINE StrLoc(corner,sLoc)
      USE Typedef
      USE StrData
      IMPLICIT NONE
      REAL(sp):: sLoc
      INTEGER(sp):: corner(8)
      sLoc=(s(corner(1))+s(corner(2))+s(corner(3))+s(corner(4))+s(corner(5))+s(corner(6))+s(corner(7))+s(corner(8)))/8._dp
      RETURN
      END
!***********************************************************************************
! current local temperature value:
      SUBROUTINE TemLoc(temp,corner,tLoc)
      USE TempData
      IMPLICIT NONE
      INTEGER(sp):: corner(8)
      LOGICAL temp
      REAL(sp):: tLoc
      IF (temp) THEN
         tLoc=(tem(corner(1))+tem(corner(2))+tem(corner(3))+tem(corner(4))+tem(corner(5))&
              +tem(corner(6))+tem(corner(7))+tem(corner(8)))/8._dp
      ELSE
         tLoc=topt
      ENDIF
      RETURN
      END

! ==============================================================================
! Source file WATER    |||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
      SUBROUTINE WATER(hNew,hOld,hTemp,t,dt,dtOpt,tOld,MatNum,Kode,iter,iter_tot,ReDo,theta&
      ,level,lDou,IAD,IADN,IADD,lOrt)

      USE tmctrl
      USE GridData, only : epslonPH,epslonR,epslonWC,epslonS,RelEps,factorRelEps,&
nPt,Axy,Dxy,Exy,sinkOld,sink,itMax
      USE MatData
      USE DoussanMat, only : PHr,PHrOld,sinkRold,SinkR,PHrtemp,PHstemp,delta2, &
switchcriterion,savelast,stressBC,PHs,delta2old,hx_min,nplant
      USE CumData, only :wBalR
      USE RootData, only :nrec
      USE WatFun
      USE SolData
      IMPLICIT NONE

      REAL(sp),intent(in):: hOld(maxnod),told
      REAL(sp),intent(inout)::hTemp(maxnod),hNew(maxnod),dt,t,dtopt,theta(maxnod)
      INTEGER(sp),intent(in):: MatNum(maxnod),Kode(maxnod)
      INTEGER(sp),intent(inout) :: iter,iter_tot
      LOGICAL,intent(out) :: ReDo
      LOGICAL,intent(in) :: lDou

      REAL(sp)::eps_WC,tol_R,tol_S,tol_PH,theta_Old(maxnod)
      REAL(sp)::delta,delta3,delta2oldmax,maxDelta2,deltaPHr,deltaPHrold,deltaS
      REAL(dp)::err_mat(maxnod),deltaPHs,deltath, maxPhrtot, maxPhr(1:maxplant)
      REAl(sp) :: A1(maxbnd,maxnod),B1(maxnod),Ti,t0,t1,sl,him,dPHsold
      INTEGER(sp) :: IAD(maxbnd,maxnod),IADN(maxnod),IADD(maxnod),iter_root,M,NORTH,iter_soil
      INTEGER(sp) :: i,level,is,it,ipl
      LOGICAL :: ItCrit,stressBC_old,lOrt,Celia=.FALSE.
      REAl(dp) ::sol(maxnod)

!----------------------- initialisation ------------------------------------
      iter_root=0
      ReDo=.FALSE.
      theta_old=theta
      IF (lDou) then
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

!print *, hTemp(1:10)
!print *, hNew(1:10)
      
      IF (level==4) stressBC_old=stressBC
      CALL SetMat(1,hOld,hTemp,hNew,MatNum,theta)
      IF (Celia) THEN
print *,'celia'!
         CALL ResetCelia(Kode,hNew,hOld,dt,lOrt,IAD,IADN,IADD,B1,theta,theta_old,iter)
      ELSE
         CALL Reset(Kode,hNew,hOld,dt,theta,theta_old,lOrt,IAD,IADN,IADD,B1)
      ENDIF
      CALL Dirich(Kode,hNew,lOrt,IADD)
      CALL SolveIlut(sol)
! test for convergence:
      ItCrit=.FALSE.!true if no convergence for soil or root
      hTemp=hNew
      hNew=sol


      iter=iter+1
      iter_tot=iter_tot+1

      IF (RelEps) THEN
         tol_PH=epslonPH
         if (epslonPH .lt. abs(maxval(hNew(1:nPt)/factorRelEps))) tol_PH= abs(maxval(hNew(1:nPt)/factorRelEps))
      ENDIF
      deltaPHs=maxval(ABS(hNew(1:nPt)-hTemp(1:nPt)))
!print *, 'errPH=',deltaPHs
      IF (deltaPHs.GT.tol_PH) THEN !.OR.deltath.gt.eps_WC concelled 24.06.09
         WRITE (*,'(''x'',$)')
         ItCrit=.TRUE.
      ENDIF
      IF ((lDou).and.(.not.(itCrit))) THEN !for level 4 in case there is convergence for soil (itCrit==.false.)
         sinkOld=sink
         DO ipl=1,nplant
            PHrTemp(1:nrec+1,ipl)=PHr(1:nrec+1,ipl)
         ENDDO
         CALL SolveRoot(hNew,hold,t,dt,.false.,level,iter_root)
!relative or absolute error
         IF (RelEps) THEN 
            DO ipl=1,nplant										!Enlever le do/enddo?(Couvreur jan 2010)
               maxPhr(ipl)=maxval(PHr(1:nrec+1,1:nplant)/factorRelEps)
            ENDDO
            maxPhrtot=maxval(maxPhr(1:nplant))
            tol_R=epslonR
            if (epslonR .lt. abs(maxPhrtot)) tol_R= abs(maxPhrtot) ! we take the largest tolerance
         ENDIF

!check for all nodes if h_value is lower than relative or absolute critetrium
         maxPhr=0
         DO ipl=1,nplant										!Enlever le do/enddo?(Couvreur jan 2010)
            maxPhr(ipl)=maxval(abs(PHr(1:nrec+1,ipl)-PHrTemp(1:nrec+1,ipl)))
         ENDDO
         deltaPHr=maxval(maxPhr(1:nplant))							!Enlever maxval?(Couvreur jan 2010)
         deltaS=maxval(abs(SinkOld-sink))
         IF ((deltaPHr.GT.tol_R).OR.(deltaS.GT.tol_S)) THEN
            write(*,'(''+'',$)')
            ItCrit=.TRUE.! no convergence for root-> rerun soil and root.
         ENDIF
      ENDIF

      if ((savelast) .and. (.not.(iTCrit)) ) return

      IF (ItCrit) THEN
         IF (iter.LT.itMax) THEN
            GOTO 250
         ELSE IF (dt.LE.dtMin) THEN
            WRITE (*,'(''Cannot converge, have to stop here. Sorry.'')')
            savelast=.true.
            if (savelast) return
         ELSE
            WRITE (*,'(''!'',$)')
            WRITE (*,'('' Iter. tot.: '',I5,$)') iter_tot
            hNew =hOld
            hTemp=hOld
            DO ipl=1,nplant
               PHr(1:nrec+1,ipl)=PHrOld(1:nrec+1,ipl)
            ENDDO
            dt=MAX(dt/3._sp,dtMin) !changed javaux july 09
            dtOpt=dt
            t=tOld+dt
            IF (level==4) stressBC=stressBC_old
            ReDo=.TRUE.
            GOTO 200
         ENDIF
     ELSE
!iteration is OK just make a last update of soil WC& PH with the right sink term
     ENDIF
     CALL WatInf(Kode,dt)
! --------------------- End of FEM iteration loop -----------------------------
200  RETURN
     END SUBROUTINE WATER

!*******************************************************************************
     SUBROUTINE SOLVEILUT(sol)
      USE MatData
      USE GridData, only : nPt

      integer(sp) :: jwork(2*maxnod),lfil,maxits=1000,ierr,iw(maxnod),ju(numnz+1)
      real(dp) :: work(maxnod+1),droptol,vv(maxnod,20),eps,alu(numnz),jlu(numnz)
      integer(sp) ::ipar(16),nwk,ipl
      real(dp):: fpar(16),ww(maxnod*40),xran(maxnod)
      REAl(dp), intent(out) ::sol(maxnod)
      external cg,bcg,dbcg,bcgstab,tfqmr,gmres,fgmres,dqgmres
      external cgnr, fom, runrc, ilut,ilu0,ulid

      !set the parameters for the iterative solvers
      ipar(2) = 1      !1 is left preconditioning,2 right,3 both, 0 = no precon
      ipar(3) = 1      !1 -> convergence test based on residualnorms; 2->based on change in approximate solution
      ipar(4) = nPt*40      !# elements in array w
      ipar(5) = 16      !size of krylov subspace
      ipar(6) = maxits      !max number of matrix-vector operations
      fpar(1) = 1.0E-8!javaux 1507 chnaged from 1.0E-6      !relative tolerance
      fpar(2) = 1.0E-12 !javaux 1507 chnaged from 1.0E-10      !absolute tolerance
      lfil = 3
      droptol =1.0E-4!javaux 1507 chnaged from   1.0e-4
      nwk = 2*numNZ

      !preconditioner
      call ilut (nPt,A_sparse,JCOL,IROW,lfil,droptol,alu,jlu,ju,nwk,work,jwork,ierr)

      !solve
      xran=0.0 !initial guess
      call runrc(nPt,B,sol,ipar,fpar,ww,xran,A_sparse,JCOL,IROW,alu,jlu,ju,bcgstab)

     END SUBROUTINE SOLVEILUT

!**************************************************************************
   SUBROUTINE runrc(n,rhs,sol,ipar,fpar,wk,guess,a,ja,ia,&
    au,jau,ju,solver)
    !  USE MatData, only: B
      implicit none
      integer n,ipar(16),ia(n+1),ja(*),ju(*),jau(*)
      real*8 fpar(16),rhs(n),sol(n),guess(n),wk(*),a(*),au(*)
      !real(kind(1e0)) fpar(16)
      external solver
!c-----------------------------------------------------------------------
!c     the actual tester. It starts the iterative linear system solvers
!c     with a initial guess suppied by the user.
!c
!c     The structure {au, jau, ju} is assumed to have the output from
!c     the ILU* routines in ilut.f.
!c
!c-----------------------------------------------------------------------
!c     local variables
!c
      integer i, iou, its
      real*8 res, dnrm2
!c     real dtime, dt(2), time
!c     external dtime
      external dnrm2
      save its,res

!ipar(2) can be 0, 1, 2, please don't use 3
      if (ipar(2).gt.2) then
         print *, 'I can not do both left and right preconditioning.'
         return
      endif

!normal execution

      its = 0
      res = 0.0D0
      do i = 1, n
         sol(i) = guess(i)
      enddo
      ipar(1) = 0

  10  call solver(n,rhs,sol,ipar,fpar,wk)

!output the residuals

      if (ipar(7).ne.its) then
         its = ipar(7)
      endif
      res = fpar(5)
!c
      if (ipar(1).eq.1) then
         call amux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
         goto 10
      else if (ipar(1).eq.2) then
         call atmux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
         goto 10
      else if (ipar(1).eq.3 .or. ipar(1).eq.5) then
         call lusol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
         goto 10
      else if (ipar(1).eq.4 .or. ipar(1).eq.6) then
         call lutsol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
         goto 10
      else if (ipar(1).le.0) then
         if (ipar(1).eq.0) then
 !            WRITE (*,'(''*'',$)')
            !Iterative solver has satisfied convergence test.'
         else if (ipar(1).eq.-1) then
            print *, 'Iterative solver has iterated too many times.'
         else if (ipar(1).eq.-2) then
            print *, 'Iterative solver was not given enough work space.'
            print *, 'The work space should at least have ', ipar(4),&
           ' elements.'
         else if (ipar(1).eq.-3) then
            print *, 'Iterative sovler is facing a break-down.'
         else
            print *, 'Iterative solver terminated. code =', ipar(1)
         endif
      endif
      return
      end
!c-----end-of-runrc
  function distdot(n,x,ix,y,iy)
      integer n, ix, iy
      real*8 distdot, x(*), y(*), ddot
      external ddot
      distdot = ddot(n,x,ix,y,iy)
      return
      end
! ==============================================================================
! Source file WATER2    |||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
      SUBROUTINE WATER2(hNew,hOld,hTemp,t,dt,dtOpt,tOld,MatNum,Kode,iter,iter_tot,ReDo,theta,theta_old,level,lDou, &
IAD,IADN,IADD,lOrt)
      USE tmctrl
      USE GridData
      USE MatData
      USE DoussanMat, only : PHr,PHrOld,sinkRold,SinkR,PHrtemp,PHstemp,delta2, &
      switchcriterion,savelast,stressBC,PHs,delta2old,tempsinkR,itcrit_root,hx_min, &
      sinkRtemp,tcheck2,nplant
      USE CumData, only :wBalR
      USE RootData, only :nrec,ibrseg,irecpr
      IMPLICIT NONE

      REAL(sp),intent(in):: hOld(maxnod),told
      REAL(sp),intent(inout)::hTemp(maxnod),hNew(maxnod),dt,t,dtopt,theta(maxnod),theta_old(maxnod)
      INTEGER(sp),intent(in):: MatNum(maxnod),Kode(maxnod)
      INTEGER(sp),intent(inout) :: iter,iter_tot
      LOGICAL,intent(out) :: ReDo
      LOGICAL,intent(in) :: lDou
      LOGICAL :: ItCrit,stressBC_old,lOrt,ItCrit_soil
      INTEGER(sp):: i,level,is
      REAL(sp)::theta_isi(maxnod),theta_isn(maxnod)!,delta(1:nPt),delta3(1:nPt),mdelta,mdelta3
      REAL(sp)::delta,delta3,maxDelta2,minPhr(1:nplant),minPhrtot
      REAl(SP) :: A1(maxbnd,maxnod),B1(maxnod)!,sinkROld(0:nrec)
      integer(sp) :: NORTH,iter_soil !MNorth=4,MaxItO=200,North
      INTEGER(sp):: IAD(maxbnd,maxnod),IADN(maxnod),IADD(maxnod),ipl
      real(sp) :: oldTempSinkR!,diffSinkR,diffSinkROld

!initialization
      theta_isi = 0
      ReDo=.FALSE.
      theta_old=theta !check related to scaling factors)
      IF (lDou) then
         PHrOld=PHr
         sinkRold=sinkR
      ENDIF
      iter_soil = iter
      252 continue
         iter_soil =iter_soil+1
!solve soil first, with last known root status
         if (level==4) stressBC_old=stressBC
         CALL SetMat(1,hOld,hTemp,hNew,MatNum,theta)
         CALL Reset(Kode,hNew,hOld,dt,lOrt,IAD,IADN,IADD,B1)
         CALL Dirich(Kode,hNew,lOrt,IADD)
         if(lOrt) then
            call ILU (A,nPt,maxbnd,IAD,IADN,IADD,A1)
            North=4!north=0 (symmetric matrix -> is not the case)
            call OrthoMin(A,B1,B,nPt,maxbnd,maxnod,IAD,IADN,IADD,A1,North)
         else
            call Solve !(A,B,NumNP,MBand,MBandD)
         end if
! test for convergence:
         ItCrit=.FALSE.!true if no convergence for soil or root
         ItCrit_soil=.false.
         hTemp=hNew
         if(lOrt) B=B1
         hNew=B
         CALL SetMat(1,hOld,hTemp,hNew,MatNum,theta)
         theta_isn = theta
         IF (RelEps) THEN
            if (epslonPH .lt. abs(maxval(hOld(1:nPt)/factorRelEps))) epslonPH= abs(maxval(hOld(1:nPt)/factorRelEps))
         ENDIF
         delta3=0
         DO is=1,nPt
            delta=ABS(hNew(is)-hTemp(is))
            delta3=ABS(theta_isn(is)-theta_isi(is))
            IF (delta.GT.epslonPH.and.delta3.gt.epslonWC) THEN
               WRITE (*,'(''-*-'',$)')
               ItCrit=.TRUE.
               ItCrit_soil=.true.
               exit
            ENDIF
         ENDDO
         IF (ItCrit_soil) THEN
            IF (iter_soil .LT.itmax) THEN
               GOTO 252
            ELSE IF (dt.LE.dtMin) THEN
               WRITE (*,'(''Cannot converge, have to stop here. Sorry.'')')
               savelast=.true.
               if (savelast) return
            ELSE
               WRITE (*,'(''!'',$)')
               WRITE (*,'(/''Iter. tot.: '',I5,$)') iter_tot
               hNew(1:nPt) =hOld(1:nPt)
               hTemp(1:nPt)=hOld(1:nPt)
               PHr=PHrOld
               dt=MAX(dt/10._sp,dtMin)!3?
               dtOpt=dt
               t=tOld+dt
               if (level==4) stressBC=stressBC_old
!                  print *,'new',t,dt,dtOpt
                  ReDo=.TRUE.
                  GOTO 200
               ENDIF
         ENDIF
!print *,'na 1x soil'
!root in iterative process with soil
         IF (.not.(ItCrit_soil) .and. (lDou)) THEN
            !only to compare old vector set to zero (old = temp in this case)
            delta2old=0
            sinkRtemp=0
            251 CONTINUE
            iter =iter+1
            iter_tot=iter_tot+1
            !BACKUP OLD PHr, PHs, radialrootflow, axialrootflow, parameter Qi
            PHrTemp=PHr
            PHsTemp=PHs
            sinkRtemp=sinkR(0:nrec,1) !WRONG FOR MULTIPLE ROOTSlinbcg

            oldTempSinkR = maxval(tempSinkR)
            itCrit =.false.
            itCrit_root=.false.
            delta2old=delta2
            CALL SolveRoot(hNew,hold,t,dt,.false.,level,iter)
            IF (RelEps) THEN

               DO ipl=1,nplant
                  minPhr(ipl)=minval(PHrOld(1:nrec+1,ipl)/factorRelEps)
               END DO
                  minPhrtot=minval(minPhr(1:nplant))
               if (epslonR .lt. abs(minPhrtot)) epslonR= abs(minPhrtot)
            ENDIF

            DO i=1,nrec+1
               delta2(i)=abs(PHrTemp(i,1)-PHr(i,1))!WRONG for MULTOIPLE ROOTS
!write(*,*)'delta2=',delta2(i)
                  IF (delta2(i).gt.epslonR) THEN
                   if (switchcriterion .ne.1) then
                      tempSinkR(i) = abs(sinkRold(i-1,1) - sinkR(i-1,1)) !sinkR(0:nrec),temp(1:nrec+1)
                      if (delta2(i).lt.epslonR+0.1*epslonR) then
                         ItCrit=.false.
                         ItCrit_root=.false.
                         exit
                       else
                         ItCrit=.TRUE.! no convergence for root (sink based)
                         itCrit_ROOT=.TRUE.
                      endif
                   else
!                      write(*,*)'delta2=',delta2(i)
!                      write(*,*)'delta2old=',delta2old(i)
                      if ((delta2(i) - delta2old(i)).lt.0.001) then
                         ItCrit=.false.
                         itCrit_ROOT=.false.
                         exit
                      endif
                      ItCrit=.TRUE.! no convergence for root (sink based)
                      itCrit_ROOT=.TRUE.
                      exit
                   endif
                  ENDIF
            ENDDO
               maxDelta2=maxval(delta2)!maximum difference of convergence value
!               write(*,*)'switchCriterion',switchCriterion
!write(*,*)'min(sum....) sinkr',sum(abs(sinkR(0:nrec)))
!write(*,*)'min(sum....) sinkrTemp',sum(abs(sinkRtemp(0:nrec)))
               IF (switchcriterion .ne.1 .and. (ItCrit) .and. iter.gt.20 .and. &
 sum(abs(sinkR(0:nrec,1))).lt.sum(abs(sinkRtemp(0:nrec))) ) then !WRONG for MULTOIPLE ROOTS
!.and. maxval(tempSinkR).lt.epslonS) then !.and. sum(abs(sinkR(0:nrec,1))).lt.sum(abs(sinkRtemp(0:nrec))) ) then !maxval(tempSinkR).lt.epslonS) then
                     ItCrit=.false.
                     ItCrit_root=.false.
               ELSEIF (switchcriterion .ne.1 .and. (ItCrit)) then
!                IF (switchcriterion .ne.1 .and. (ItCrit)) then
                  do i=1,nrec+1
! !                  write(*,*)'delta2-delta2old',abs(delta2(i)-delta2old(i))
!                    ! write(*,*)'delta2=',delta2(i)
!                    ! write(*,*)'delta2old=',delta2old(i)
!                      write(*,*)'maxval(tempSinkR)',maxval(tempSinkR)
                     if (OldTempSinkR.gt. maxval(tempSinkR) .and. maxval(tempSinkR) .lt. epslonS) then
!!(abs(delta2(i)-delta2old(i)).lt.0.5) .and.
                        ItCrit=.false.
                        ItCrit_root=.false.
                        exit
                     endif
                  enddo
                ENDIF
             if ((savelast) .and. (.not.(iTCrit))) return
             IF (.not.(itCrit)) THEN
               theta_isi=theta
               CALL SetMat(1,hOld,hTemp,hNew,MatNum,theta)
               CALL Reset(Kode,hNew,hOld,dt,lOrt,IAD,IADN,IADD,B1)
               CALL Dirich(Kode,hNew,lOrt,IADD)
               if(lOrt) then
                  call ILU (A,nPt,maxbnd,IAD,IADN,IADD,A1)
                  North=4!north=0 (symmetric matrix -> is not the case)
                  call OrthoMin(A,B1,B,nPt,maxbnd,maxnod,IAD,IADN,IADD,A1,North)
               else
                  call Solve !(A,B,NumNP,MBand,MBandD)
               endif
               ! test for convergence:
       !        ItCrit=.FALSE.!true if no convergence for soil or root
       !        ItCrit_soil=.false.
               hTemp=hNew
               if(lOrt) B=B1
               hNew=B
               CALL SetMat(1,hOld,hTemp,hNew,MatNum,theta)
               theta_isn = theta
               IF (RelEps) THEN
                  if (epslonPH .lt. abs(maxval(hOld(1:nPt)/1000))) epslonPH= abs(maxval(hOld(1:nPt)/1000))
               ENDIF
               DO is=1,nPt
                  delta=ABS(hNew(is)-hTemp(is))
                  delta3=ABS(theta_isn(is)-theta_isi(is))
                  IF (delta.GT.epslonPH.and.delta3.gt.epslonWC) THEN
                     WRITE (*,'(''-2emalsoil-'',$)')
                     ItCrit=.TRUE.
                     exit
                  ENDIF
               ENDDO
            ENDIF
            IF (ItCrit) THEN
               IF (iter.LT.itMax) THEN
                  GOTO 251
               ELSEIF (dt.LE.dtMin) THEN
                  WRITE (*,'(''Cannot converge, have to stop here. Sorry.'')')
                  savelast=.true.
                  if (savelast) return
               ELSE
                  WRITE (*,'(''!'',$)')
                  WRITE (*,'(/''Iter. tot.: '',I5,$)') iter_tot
                     hNew =hOld
                     hTemp=hOld
                     PHr=PHrOld
!                  print *,'old',t,dt,'to',tOld
                  dt=MAX(dt/10._sp,dtMin)!3?
                  dtOpt=dt
                  t=tOld+dt
                  if (level==4) stressBC=stressBC_old
!                     print *,'new',t,dt,dtOpt
                     ReDo=.TRUE.
                     GOTO 200
               ENDIF
            ENDIF
         ENDIF
         IF (tcheck2) then
            dt=MAX(dt/1000._sp,dtMin)
            dtOpt=dt
            tcheck2=.false.
         ENDIF
!         do i=0,nrec
!               IF (PHs(i).lt.hx_min) then
!               !     PHsTop=PHs(i)
!                  savelast=.true.
!               ENDIF
!          enddo
      CALL WatInf(Kode,dt)
! --------------------- End of FEM iteration loop -----------------------------
200   RETURN
      END
!*********************************************************************
      SUBROUTINE SetBC(t,h,Kode)
      USE BoundData
      USE GridData
      USE CumData
      IMPLICIT NONE
      REAL(sp) :: h(maxnod),t,irrigflw,head,Qtot,wtot
      REAL(sp) , dimension (:) :: volflw(1+(homogene-1)*mxBcCh)
      INTEGER(sp) :: Kode(maxnod),i
! calculate current value for Q-BC's, if any:
      IF (nQbcCh.EQ.0) GOTO 25
      i=0
   22 i=i+1
      IF (i.GT.nQbcCh) THEN
         volflw=Qbc(nQbcCh,:)
      ELSE
         IF (t.GT.tQbcCh(i)) GOTO 22
         volflw=Qbc(i-1,:)+(t-tQbcCh(i-1))/(tQbcCh(i)-tQbcCh(i-1))*(Qbc(i,:)-Qbc(i-1,:))
      ENDIF
! calculate current value for I-BC's, if any:				!(Couvreur jan 2010)
   25 IF (nIbcCh.EQ.0) GOTO 30
      i=0
   28 i=i+1
      IF (i.GT.nIbcCh) THEN
         irrigflw=Ibc(nIbcCh)
      ELSE
         IF (t.GT.tIbcCh(i)) GOTO 28
         irrigflw=Ibc(i-1)+(t-tIbcCh(i-1))/(tIbcCh(i)-tIbcCh(i-1))*(Ibc(i)-Ibc(i-1))
      ENDIF
! calculate current value for h-BC's, if any:
   30 IF (nhbcCh.EQ.0) GOTO 40
      i=0
   33 i=i+1
      IF (i.GT.nhbcCh) THEN
         head=hbc(nhbcCh)
      ELSE
         IF (t.GT.thbcCh(i)) GOTO 33
         head=hbc(i-1)+(t-thbcCh(i-1))/(thbcCh(i)-thbcCh(i-1))*(hbc(i)-hbc(i-1))
      ENDIF
!calculate QBC*surface
   40 Qtot=0
      DO 44 i=1,nBCPts
         IF (Kode(iBCPt(i)).EQ.+1) THEN
            h(iBCPt(i))=head
         ELSEIF (Kode(iBCPt(i)).EQ.-1) THEN
            IF (homogene==1) THEN						!(Couvreur jan 2010)
	        Q(iBCPt(i))=volflw(1)*Width(i)
            ELSE
               Q(iBCPt(i))=volflw(i)*Width(i)
            ENDIF
	     Qtot=Qtot+Q(iBCPt(i))
            wtot=wtot+Width(i)
         ELSEIF (Kode(iBCPt(i)).EQ.-3) THEN
	     Q(iBCPt(i))=irrigflw*nx*ny*dxgrid*dygrid			!(Couvreur jan 2010)
	     Qtot=Qtot+Q(iBCPt(i))
            wtot=wtot+Width(i)
         ENDIF
   44 CONTINUE
! calculate current value for solute transport boundary conditions, if any:
      IF (nCBnd1.EQ.0) GOTO 50
      i=0
   55 i=i+1
      IF (i.GT.nCBnd1) THEN
         cBound(1)=CBnd1(nCBnd1)
      ELSE
         IF (t.GT.tCBnd1(i)) GOTO 55
         cBound(1)=CBnd1(i-1)+(t-tCBnd1(i-1))/(tCBnd1(i)-tCBnd1(i-1))*(CBnd1(i)-CBnd1(i-1))
      ENDIF
   50 IF (nCBnd2.EQ.0) GOTO 60
      i=0
   66 i=i+1
      IF (i.GT.nCBnd2) THEN
         cBound(2)=CBnd2(nCBnd2)
      ELSE
         IF (t.GT.tCBnd2(i)) GOTO 66
         cBound(2)=CBnd2(i-1)+(t-tCBnd2(i-1))/(tCBnd2(i)-tCBnd2(i-1))*(CBnd2(i)-CBnd2(i-1))
      ENDIF
   60 RETURN
      END
!*********************************************************************************
      SUBROUTINE CalcWidthSWMS
!Javaux May08 from SWMS_3D
      USE GridData
      USE DomData
      IMPLICIT NONE
      INTEGER(sp):: iB,ix,iy,iCord,iDiv,iRest
      REAL (sp):: WidthX,WidthY
      iB=0
      DO iy=1,ny
        DO ix=1,nx
          iB=iB+1
          if(ix.ne.1.and.ix.ne.nx) then
            WidthX=(xCol(ix+1)-xCol(ix-1))
          else if(ix.eq.1) then
            WidthX=(xCol(ix+1)-xCol(1))
	     if(continu) WidthX=(xCol(2)-xCol(nx)+nx*dxGrid)				!(Couvreur dec 2009)
          else if(ix.eq.nx) then
            WidthX=(xCol(nx)-xCol(nx-1))
	     if(continu) WidthX=(xCol(1)-xCol(nx-1)+nx*dxGrid)				!(Couvreur dec 2009)
          end if
          if(iy.ne.1.and.iy.ne.ny) then
            WidthY=(yCol(iy+1)-yCol(iy-1))
          else if(iy.eq.1) then
            WidthY=(yCol(iy+1)-yCol(1))
	     if(continu) WidthY=(yCol(2)-yCol(ny)+ny*dyGrid)				!(Couvreur dec 2009)
          else if(iy.eq.ny) then
            WidthY=(yCol(ny)-yCol(ny-1))
	     if(continu) WidthY=(yCol(1)-yCol(ny-1)+ny*dyGrid)				!(Couvreur dec 2009)
          end if
          iCord=ix+iy
          iRest=mod(iCord,2)
          iDiv=3
          if(iRest.eq.1) iDiv=6
          Width(iB)=WidthX*WidthY/iDiv
          iCord=ix+iy+nz
          iRest=mod(iCord,2)
          iDiv=3
          if(iRest.eq.0) iDiv=6
          Width(iB+ny*nx)=WidthX*WidthY/iDiv
	 ENDDO
      ENDDO
      RETURN
      END
!************************************************************************************
      SUBROUTINE SetMat(K,hOld,hTemp,hNew,MatNum,theta)

      USE GridData, only : nPt,Axy,Bxy,Dxy,Exy
      USE SolData
      USE WatFun
      IMPLICIT NONE

      REAL(sp), intent(in):: hNew(maxnod),hTemp(maxnod),hOld(maxnod)
      REAL(sp):: ci,cpi,Ti,t0,t1
      REAL(sp):: sl,him,hi2,hi1
      REAL(sp), intent(out) :: theta(maxnod)
      INTEGER(sp):: MatNum(maxnod),it,i,K,m

!$OMP PARALLEL DO private(M,hi1,hi2,hiM,iT,Sl,ci,Cpi,Ti)
      DO i=1,nPt
! calculate nodal conductivity values:
         IF(K.EQ.1) conO(i)=con(i)
         M=MatNum(i)
         hi1=hTemp(i)/Axy(i) !javaux
         hi2=hNew(i)/Axy(i)!javaux
         hiM=0.1_dp*hi1+0.9_dp*hi2
         IF ((hiM.GE.hTab(nTab)).AND.(hiM.LE.hTab(1))) THEN
            iT=INT((LOG10(-hiM)-alh1)/dlh)+1
            Sl=(ConTab(iT+1,M)-ConTab(iT,M))/(hTab(iT+1)-hTab(iT))
            ci=ConTab(iT,M)+Sl*(hiM-hTab(iT))
         ELSE
            ci=FKP(hiM,par(:,M))
         ENDIF
         con(i)=ci*Bxy(i)
         IF(K.EQ.0) conO(i)=con(i)

! calculate nodal capacity values:
         hi1=hOld(i)/Axy(i)
         hi2=hNew(i)/Axy(i)
         hiM=(hi2+hi1)/2
         IF ((hiM.GE.hTab(nTab)).AND.(hiM.LE.hTab(1))) THEN
            iT=INT((LOG10(-hiM)-alh1)/dlh)+1
            Sl =(CapTab(iT+1,M)-CapTab(iT,M))/(hTab(iT+1)-hTab(iT))
            Cpi=CapTab(iT,M)+Sl*(hiM-hTab(iT))
            Ti=TheTab(iT,M)+Sl*(hiM-hTab(iT))
         ELSE
            Cpi=FCP(hiM,par(:,M))
            Ti=Fth(hiM,par(:,M))
         ENDIF
         cap(i)=Cpi*Dxy(i)/Axy(i)
         theta(i)=par(2,M)*Exy(i)+(Ti-par(2,M))*Dxy(i)
     ENDDO
!$OMP END PARALLEL DO
!call cpu_time(t1)

      RETURN
      END SUBROUTINE SetMat !MJ09
!****************************************************************************
     SUBROUTINE CalcGeom
      USE GridData
      USE CumData
      USE DomData
      USE MatData
      IMPLICIT NONE

      INTEGER(sp) :: iE,iSE,i,j,k,l,ii,jj,iL(4)
      REAL(sp) :: bi(1:4),ci(1:4),di(1:4),ai(4),xGridi,xGridj,xGridk,xGridl,yGridi,yGridj,yGridk,yGridl
      REAL(sp) :: Cxx,Czz,Cyy,Cxy,Cxz,Cyz
      REAL(dp) :: Det
      INTEGER(sp) :: ex,ey

! Loop on elements:
      ex=0
      ey=1
      DO 200 iE=1,nElm  !half cubes
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
         ! Loop on subelements: 3 tetrahedrals per half soil cube
         DO 120 iSE=1,3
            IF (mod(iE,2).eq.0) THEN !(xGrid(elmnod(iE,1))-xGrid(elmnod(iE,2)).gt.0) then !4 1 2 8 5 6!
               iL(1)=3+iSE
               iL(2)=4+iSE-iSE/3*6
               iL(3)=5+iSE-iSE/2*6
               iL(4)=  iSE
            ELSEIF (mod(iE,2).eq.1) THEN !(xGrid(elmnod(iE,1))-xGrid(elmnod(iE,2)).lt.0) then !5 8 7 1 4 3 !
               iL(1)=3+iSE
               iL(2)=5+iSE-iSE/2*6
               iL(3)=4+iSE-iSE/3*6
               iL(4)=  iSE
            ENDIF
            i=elmnod(iE,iL(1)) !i,j,k,l are corner nodes of the tetraedal element, ise counts 3 which is equal to three tedrahedals which is equal to half a cubic.
            j=elmnod(iE,iL(2))
            k=elmnod(iE,iL(3))
            l=elmnod(iE,iL(4))
            !if master node then coefficient calculated as normal; secondly signs are needed for orientation
	     IF (((ex.GT.nex-2).OR.(ey.EQ.ney)).AND.continu) THEN					!"Bridge sub-element?" (Couvreur dec 2009)
		 xGridi=xGrid(i)
		 xGridj=xGrid(j)
		 xGridk=xGrid(k)
		 xGridl=xGrid(l)
		 yGridi=yGrid(i)
		 yGridj=yGrid(j)
		 yGridk=yGrid(k)
		 yGridl=yGrid(l)
	        IF (ex.GT.nex-2) THEN
	           IF (mod(iE,2).eq.0) THEN				!Second semi-cube  (modif. position of elmnod(iE,[1 3 4 6]))
		       IF (iSE.EQ.1) THEN
		          xGridi=xGridi+nex*dxgrid/2
		          xGridl=xGridl+nex*dxgrid/2
		          IF (iL(2).EQ.6) THEN
		             xGridj=xGridj+nex*dxgrid/2
		          ELSEIF (iL(3).EQ.6) THEN
		             xGridk=xGridk+nex*dxgrid/2
		          ENDIF
		       ENDIF
		       IF (iSE.EQ.2) THEN
		          xGridj=xGridj+nex*dxgrid/2
		          xGridk=xGridk+nex*dxgrid/2
		       ENDIF
		       IF (iSE.EQ.3) THEN
		          xGridi=xGridi+nex*dxgrid/2
		          xGridl=xGridl+nex*dxgrid/2
		          IF (iL(2).EQ.1) THEN
		             xGridj=xGridj+nex*dxgrid/2
		          ELSEIF (iL(3).EQ.1) THEN
		             xGridk=xGridk+nex*dxgrid/2
		          ENDIF
		       ENDIF
	           ELSEIF (mod(iE,2).eq.1) THEN				!First semi-cube   (modif. position of elmnod(iE,[2 5]))
		       IF (iSE.EQ.1) THEN
		          IF (iL(2).EQ.5) THEN
		             xGridj=xGridj+nex*dxgrid/2
		          ELSEIF (iL(3).EQ.5) THEN
		             xGridk=xGridk+nex*dxgrid/2
		          ENDIF
		       ENDIF
		       IF (iSE.EQ.2) THEN
		          xGridi=xGridi+nex*dxgrid/2
		          xGridl=xGridl+nex*dxgrid/2
		       ENDIF
		       IF (iSE.EQ.3) THEN
		          IF (iL(2).EQ.2) THEN
		             xGridj=xGridj+nex*dxgrid/2
		          ELSEIF (iL(3).EQ.2) THEN
		             xGridk=xGridk+nex*dxgrid/2
		          ENDIF
		       ENDIF
		    ENDIF
		 ENDIF
	        IF (ey.EQ.ney) THEN
	           IF (mod(iE,2).eq.0) THEN				!Second semi-cube  (modif. position of elmnod(iE,[1 4]))
		       IF (iSE.EQ.1) THEN
		          yGridi=yGridi+ney*dygrid
		          yGridl=yGridl+ney*dygrid
		       ENDIF
		       IF (iSE.EQ.2) THEN
		          IF (iL(2).EQ.1) THEN
		             yGridj=yGridj+ney*dygrid
		          ELSEIF (iL(3).EQ.1) THEN
		             yGridk=yGridk+ney*dygrid
		          ENDIF
		       ENDIF
		       IF (iSE.EQ.3) THEN
		          IF (iL(2).EQ.1) THEN
		             yGridj=yGridj+ney*dygrid
		          ELSEIF (iL(3).EQ.1) THEN
		             yGridk=yGridk+ney*dygrid
		          ENDIF
		       ENDIF
	           ELSEIF (mod(iE,2).eq.1) THEN				!First semi-cube   (modif. position of elmnod(iE,[2 3 5 6]))
		       IF (iSE.EQ.1) THEN
		          yGridj=yGridj+ney*dygrid
		          yGridk=yGridk+ney*dygrid
		       ENDIF
		       IF (iSE.EQ.2) THEN
		          yGridi=yGridi+ney*dygrid
		          yGridl=yGridl+ney*dygrid
		          IF (iL(2).EQ.6) THEN
		             yGridj=yGridj+ney*dygrid
		          ELSEIF (iL(3).EQ.6) THEN
		             yGridk=yGridk+ney*dygrid
		          ENDIF
		       ENDIF
		       IF (iSE.EQ.3) THEN
		          yGridi=yGridi+ney*dygrid
		          yGridl=yGridl+ney*dygrid
		          IF (iL(2).EQ.2) THEN
		             yGridj=yGridj+ney*dygrid
		          ELSEIF (iL(3).EQ.2) THEN
		             yGridk=yGridk+ney*dygrid
		          ENDIF
		       ENDIF
		    ENDIF
		 ENDIF
               bi(1)=-(yGridk-yGridj)*(zGrid(l)-zGrid(j))+(yGridl-yGridj)*(zGrid(k)-zGrid(j))
               bi(2)=+(yGridl-yGridk)*(zGrid(i)-zGrid(k))-(yGridi-yGridk)*(zGrid(l)-zGrid(k))
               bi(3)=-(yGridi-yGridl)*(zGrid(j)-zGrid(l))+(yGridj-yGridl)*(zGrid(i)-zGrid(l))
               bi(4)=+(yGridj-yGridi)*(zGrid(k)-zGrid(i))-(yGridk-yGridi)*(zGrid(j)-zGrid(i))
               ci(1)=+(xGridk-xGridj)*(zGrid(l)-zGrid(j))-(xGridl-xGridj)*(zGrid(k)-zGrid(j))
               ci(2)=-(xGridl-xGridk)*(zGrid(i)-zGrid(k))+(xGridi-xGridk)*(zGrid(l)-zGrid(k))
               ci(3)=+(xGridi-xGridl)*(zGrid(j)-zGrid(l))-(xGridj-xGridl)*(zGrid(i)-zGrid(l))
               ci(4)=-(xGridj-xGridi)*(zGrid(k)-zGrid(i))+(xGridk-xGridi)*(zGrid(j)-zGrid(i))
               di(1)=-(xGridk-xGridj)*(yGridl-yGridj)+(xGridl-xGridj)*(yGridk-yGridj)
               di(2)=+(xGridl-xGridk)*(yGridi-yGridk)-(xGridi-xGridk)*(yGridl-yGridk)
               di(3)=-(xGridi-xGridl)*(yGridj-yGridl)+(xGridj-xGridl)*(yGridi-yGridl)
               di(4)=+(xGridj-xGridi)*(yGridk-yGridi)-(xGridk-xGridi)*(yGridj-yGridi)
    ! coefficient a-> shape_function = a + bx +cy +dz
               ai(1)=xGridj*yGridk*zGrid(l) + xGridk*yGridl*zGrid(j) + xGridl*yGridj*zGrid(k) - &
                    xGridl*yGridk*zGrid(j) - xGridj*yGridl*zGrid(k) - xGridk*yGridj*zGrid(l)
               ai(2)=xGridi*yGridl*zGrid(k) + xGridk*yGridi*zGrid(l) + xGridl*yGridk*zGrid(i) - &
                    xGridi*yGridk*zGrid(l) - xGridk*yGridl*zGrid(i) - xGridl*yGridi*zGrid(k)
               ai(3)=xGridi*yGridj*zGrid(l) + xGridj*yGridl*zGrid(i) + xGridl*yGridi*zGrid(j) - &
                    xGridi*yGridl*zGrid(j) - xGridj*yGridi*zGrid(l) - xGridl*yGridj*zGrid(i)
               ai(4)=- xGridi*yGridj*zGrid(k) - xGridj*yGridk*zGrid(i) - xGridk*yGridi*zGrid(j) + &
                    xGridi*yGridk*zGrid(j) + xGridj*yGridi*zGrid(k) + xGridk*yGridj*zGrid(i)
               Det=(xGridl-xGridi)*bi(4)+(yGridl-yGridi)*ci(4)+(zGrid(l)-zGrid(i))*di(4)
            ELSE
               bi(1)=-(yGrid(k)-yGrid(j))*(zGrid(l)-zGrid(j))+(yGrid(l)-yGrid(j))*(zGrid(k)-zGrid(j))
               bi(2)=+(yGrid(l)-yGrid(k))*(zGrid(i)-zGrid(k))-(yGrid(i)-yGrid(k))*(zGrid(l)-zGrid(k))
               bi(3)=-(yGrid(i)-yGrid(l))*(zGrid(j)-zGrid(l))+(yGrid(j)-yGrid(l))*(zGrid(i)-zGrid(l))
               bi(4)=+(yGrid(j)-yGrid(i))*(zGrid(k)-zGrid(i))-(yGrid(k)-yGrid(i))*(zGrid(j)-zGrid(i))
               ci(1)=+(xGrid(k)-xGrid(j))*(zGrid(l)-zGrid(j))-(xGrid(l)-xGrid(j))*(zGrid(k)-zGrid(j))
               ci(2)=-(xGrid(l)-xGrid(k))*(zGrid(i)-zGrid(k))+(xGrid(i)-xGrid(k))*(zGrid(l)-zGrid(k))
               ci(3)=+(xGrid(i)-xGrid(l))*(zGrid(j)-zGrid(l))-(xGrid(j)-xGrid(l))*(zGrid(i)-zGrid(l))
               ci(4)=-(xGrid(j)-xGrid(i))*(zGrid(k)-zGrid(i))+(xGrid(k)-xGrid(i))*(zGrid(j)-zGrid(i))
               di(1)=-(xGrid(k)-xGrid(j))*(yGrid(l)-yGrid(j))+(xGrid(l)-xGrid(j))*(yGrid(k)-yGrid(j))
               di(2)=+(xGrid(l)-xGrid(k))*(yGrid(i)-yGrid(k))-(xGrid(i)-xGrid(k))*(yGrid(l)-yGrid(k))
               di(3)=-(xGrid(i)-xGrid(l))*(yGrid(j)-yGrid(l))+(xGrid(j)-xGrid(l))*(yGrid(i)-yGrid(l))
               di(4)=+(xGrid(j)-xGrid(i))*(yGrid(k)-yGrid(i))-(xGrid(k)-xGrid(i))*(yGrid(j)-yGrid(i))
    ! coefficient a-> shape_function = a + bx +cy +dz
               ai(1)=xGrid(j)*yGrid(k)*zGrid(l) + xGrid(k)*yGrid(l)*zGrid(j) + xGrid(l)*yGrid(j)*zGrid(k) - &
                    xGrid(l)*yGrid(k)*zGrid(j) - xGrid(j)*yGrid(l)*zGrid(k) - xGrid(k)*yGrid(j)*zGrid(l)
               ai(2)=xGrid(i)*yGrid(l)*zGrid(k) + xGrid(k)*yGrid(i)*zGrid(l) + xGrid(l)*yGrid(k)*zGrid(i) - &
                    xGrid(i)*yGrid(k)*zGrid(l) - xGrid(k)*yGrid(l)*zGrid(i) - xGrid(l)*yGrid(i)*zGrid(k)
               ai(3)=xGrid(i)*yGrid(j)*zGrid(l) + xGrid(j)*yGrid(l)*zGrid(i) + xGrid(l)*yGrid(i)*zGrid(j) - &
                    xGrid(i)*yGrid(l)*zGrid(j) - xGrid(j)*yGrid(i)*zGrid(l) - xGrid(l)*yGrid(j)*zGrid(i)
               ai(4)=- xGrid(i)*yGrid(j)*zGrid(k) - xGrid(j)*yGrid(k)*zGrid(i) - xGrid(k)*yGrid(i)*zGrid(j) + &
                    xGrid(i)*yGrid(k)*zGrid(j) + xGrid(j)*yGrid(i)*zGrid(k) + xGrid(k)*yGrid(j)*zGrid(i)
               Det=(xGrid(l)-xGrid(i))*bi(4)+(yGrid(l)-yGrid(i))*ci(4)+(zGrid(l)-zGrid(i))*di(4)
            ENDIF
            VE(iE,iSE)=abs(Det)/6.
            if (Det.eq.0 .or. VE(iE,iSE).eq.0) then
               write(*,*)'Ve=0',VE,'det=',Det,'element nr.',iE
               print*,'coo i',xgrid(i),ygrid(i),zgrid(i)
               print*,'coo j',xgrid(j),ygrid(j),zgrid(j)
               print*,'coo k',xgrid(k),ygrid(k),zgrid(k)
               print*,'coo l',xgrid(l),ygrid(l),zgrid(l)
               pause
            endif
            DO ii=1,4
               B1fact(iE,iSE,ii)=Cxz*bi(ii)+Cyz*ci(ii)+Czz*di(ii)
            END DO
            DO 110 ii=1,4
               DO 100 jj=1,4
                  E(iE,iSE,ii,jj)=Cxx*bi(ii)*bi(jj)+Cyy*ci(ii)*ci(jj)+Czz*di(ii)*di(jj) &
                  +Cxy*(bi(ii)*ci(jj)+bi(jj)*ci(ii))+Cxz*(bi(ii)*di(jj)+bi(jj)*di(ii)) &
                  +Cyz*(ci(ii)*di(jj)+ci(jj)*di(ii))
               100  CONTINUE
            110  CONTINUE
         120  CONTINUE
      200  CONTINUE
      RETURN
      END SUBROUTINE CalcGeom !Couvreur jan 2010
!****************************************************************************
     SUBROUTINE Reset(Kode,hNew,hOld,dt,theta,theta_old,lOrt,IAD,IADN,IADD)
      USE ParamData, only: maxbnd
      USE GridData
      USE PlntData, only: Tpot,Tact
      USE CumData
      USE SolData
      USE DomData
      USE MatData
      USE BoundData,only: iBCPt
      IMPLICIT NONE

      REAL(sp) :: hNew(maxnod),hold(maxnod),DS(maxnod)
      REAL(sp) :: theta(maxnod),theta_old(maxnod)
      REAL(sp) :: F(maxnod)
      INTEGER(sp) :: Kode(maxnod),iL(4),ind,ind1
      INTEGER(sp) :: iSE,i,j,k,l,iE,n,ii,jj,kk,iG,jG
      REAL(sp) :: cone,cape
      REAL(sp) :: fmul,bmul
      REAL(sp) :: amul,QN
      REAL(sp) :: dt,t1,t0
      REAL(dp) :: SinkE
      REAL(dp) :: B1(maxnod)
      REAL(sp) :: vol_elem_sink,VE_ori
      INTEGER(sp) :: IAD(maxbnd,maxnod),IADN(maxnod),IADD(maxnod)
      LOGICAL :: lOrt,yes

      do j=1,irow(npt+1)-1
         a_sparse(j) = 0.0d0
      enddo
      B(1:npt)=0.0
      B1=0.0
      DS=0.0
      F=0.0
      RootSkold=0.0
      vol_elem_sink=0
! Loop on elements:
      DO 200 iE=1,nElm  !half cubes
         ! Loop on subelements: 3 tetrahedrals per half soil cube
         DO 120 iSE=1,3
            IF (mod(iE,2).eq.0) THEN !(xGrid(elmnod(iE,1))-xGrid(elmnod(iE,2)).gt.0) then !4 1 2 8 5 6!
               iL(1)=3+iSE
               iL(2)=4+iSE-iSE/3*6
               iL(3)=5+iSE-iSE/2*6
               iL(4)=  iSE
            ELSEIF (mod(iE,2).eq.1) THEN !(xGrid(elmnod(iE,1))-xGrid(elmnod(iE,2)).lt.0) then !5 8 7 1 4 3 !
               iL(1)=3+iSE
               iL(2)=5+iSE-iSE/2*6
               iL(3)=4+iSE-iSE/3*6
               iL(4)=  iSE
            ENDIF
            i=elmnod(iE,iL(1)) !i,j,k,l are corner nodes of the tetraedal element, ise counts 3 which is equal to three tedrahedals which is equal to half a cubic.
            j=elmnod(iE,iL(2))
            k=elmnod(iE,iL(3))
            l=elmnod(iE,iL(4))
            CapE=(Cap(i)+Cap(j)+Cap(k)+Cap(l))/4
            ConE=(Con(i)+Con(j)+Con(k)+Con(l))/4
            AMul=ConE/VE(iE,iSE)/36
            BMul=ConE/6
            FMul=VE(iE,iSE)/20
            SinkE=(sink(i)+sink(j)+sink(k)+sink(l))/4
            DS(i)=DS(i)+FMul*(4*SinkE+sink(i))
            DS(j)=DS(j)+FMul*(4*SinkE+sink(j))
            DS(k)=DS(k)+FMul*(4*SinkE+sink(k))
            DS(l)=DS(l)+FMul*(4*SinkE+sink(l))
!flow leaving soil is average sink times their volume  over all the elements
            RootSkold=RootSkold+VE(iE,iSE)*SinkE
            if (tpot.ne.0) then!
                  !              write(*,*)'Rootsk=',RootSk,'SE',SinkE
                   !             write(*,*)'sink per node=',sink(k),sink(j),sink(l),sink(i)
               vol_elem_sink = vol_elem_sink + VE(iE,iSE)
            endif
            DO 110 ii=1,4
               iG=elmnod(iE,iL(ii))
               F(iG)=F(iG)+FMul*(4*CapE+Cap(iG)) !*5
               B1(iG)=B1(iG)+BMul*B1fact(iE,iSE,ii)
               DO 100 jj=1,4
                  jG=elmnod(iE,iL(jj))
                  call Find(iG,jG,kk,nPt,maxbnd,IAD,IADN)
                  ! normal value added to array A_sparse (CSR format)
                  A_sparse(IROW(iG)-1+kk)=A_sparse(IROW(iG)-1+kk)+AMul*E(iE,iSE,ii,jj)
               100  CONTINUE
            110  CONTINUE
         120  CONTINUE
      200  CONTINUE
       
      DO 240 N=1,nPt
         IF (Kode(N).gt.0) then !hbc
            QN=B1(N)+DS(N)!+F(N)*(theta(N)-theta_old(N))/dt !Tom had deleted the term +F(N) etc... same than SWMS
            ind=IROW(N)
            ind1=IROW(N+1)-1
            DO j=ind,ind1
               QN=QN+A_sparse(j)*hNew(JCOL(j))
               !if a row is part of a bound.cond. then get the values in A_sparse and multiply by the corresponding column values in hNew
            enddo
            Q(N)=QN
         elseif (Kode(N).EQ.0) then !no BC defined
            Q(N)=0
!         elseif (Kode(N).EQ.-1) then !rain (homogeneous or not)
         elseif (Kode(N).EQ.-2) then !free drainage
            Q(N)=-con(N)*Width(N)
         endif
      240 CONTINUE
! Complete construction of RHS vector and form effective matrix:
! use row index - 1 to get the number of previous entries, add up the entry number of the diagonal value to find the location of the diagonal value in A_sparse
      DO 310 i=1,nPt
         B(i)=Q(i)-B1(i)+hold(i)*F(i)/dt-DS(i)
         A_sparse(IROW(i)-1+IADD(i)) = A_sparse(IROW(i)-1+IADD(i)) +  F(i)/dt !CSR format
      310 CONTINUE
      RETURN
      END
!***********************************************************************
     SUBROUTINE ResetCelia(Kode,hNew,hOld,dt,lOrt,IAD,IADN,IADD,ThNew,ThOld,iter)

      USE ParamData, only: maxbnd
      USE GridData
      USE PlntData, only: Tpot,Tact
      USE CumData
      USE SolData
      USE DomData
      USE MatData
      USE BoundData,only: iBCPt
      IMPLICIT NONE

      REAL(sp), intent(in) :: hNew(maxnod),hold(maxnod),ThNew(maxnod),ThOld(maxnod)
      INTEGER, intent(in) :: iter
      REAL(sp) :: DS(maxnod)
      REAL(sp) :: F(maxnod),EE(4,4)
      INTEGER(sp) :: Kode(maxnod),iL(4),ind,ind1,teller,index2,index3
      INTEGER(sp) :: ise,k,l,nnz,nny,nnx,ib,i,j,ie,n,ii,jj,kk,nn,ig,jg
      integer(sp)::igg(12),jgg(12),sign_a,sign_b,sign_c,sign_d
      real(sp):: A_test(5,5)
      REAL(sp) :: cone,cape
      REAL(sp) :: bi(1:4),ci(1:4),di(1:4),ai(4),fmul,bmul
      REAL(sp) :: cxx,czz,cyy,cxy,cxz,cyz,amul,qn
      REAL(sp) :: dt,t1,t0,shape_func_1,shape_func_2,shape_func_3,shape_func_4
      REAL(dp) :: SinkE,det
      REAL(dp) :: B1(maxnod)
      REAL(sp) :: vol_elem_sink,VE_ori
      INTEGER(sp) :: IAD(maxbnd,maxnod),IADN(maxnod),IADD(maxnod),index,tel
      LOGICAL :: lOrt,yes

      a_sparse(1:irow(npt+1)-1)=0.00 
      B(1:npt)=0.0_dp
      B1(1:npt)=hNew
      DS=0.0 !is estimated at each iteration step
      F=0.0
      RootSkold=0.0
      vol_elem_sink=0

! Loop on elements:
      DO 200 iE=1,nElm  !half cubes
         Cxx=ConAxx(iE)
         Cyy=ConAyy(iE)
         Czz=ConAzz(iE)
         Cxy=ConAxy(iE)
         Cxz=ConAxz(iE)
         Cyz=ConAyz(iE)

! Loop on subelements: 3 tetrahedrals per half soil cube
         DO 120 iSE=1,3
            IF (xGrid(elmnod(iE,1))-xGrid(elmnod(iE,2)).gt.0) then !4 1 2 8 5 6!(mod(iE,2).eq.0) THEN
               iL(1)=3+iSE
               iL(2)=4+iSE-iSE/3*6
               iL(3)=5+iSE-iSE/2*6
               iL(4)=  iSE
            ELSEIF (xGrid(elmnod(iE,1))-xGrid(elmnod(iE,2)).lt.0) then !5 8 7 1 4 3 !(mod(iE,2).eq.1) THEN
               iL(1)=3+iSE
               iL(2)=5+iSE-iSE/2*6
               iL(3)=4+iSE-iSE/3*6
               iL(4)=  iSE
            ENDIF
            i=elmnod(iE,iL(1)) !i,j,k,l are corner nodes of the tetraedal element, ise counts 3 which is equal to three tedrahedals which is equal to half a cubic.
            j=elmnod(iE,iL(2))
            k=elmnod(iE,iL(3))
            l=elmnod(iE,iL(4))

!if master node then coefficient calculated as normal; secondly signs are needed for orientation
            bi(1)=-(yGrid(k)-yGrid(j))*(zGrid(l)-zGrid(j))+(yGrid(l)-yGrid(j))*(zGrid(k)-zGrid(j))
            bi(2)=+(yGrid(l)-yGrid(k))*(zGrid(i)-zGrid(k))-(yGrid(i)-yGrid(k))*(zGrid(l)-zGrid(k))
            bi(3)=-(yGrid(i)-yGrid(l))*(zGrid(j)-zGrid(l))+(yGrid(j)-yGrid(l))*(zGrid(i)-zGrid(l))
            bi(4)=+(yGrid(j)-yGrid(i))*(zGrid(k)-zGrid(i))-(yGrid(k)-yGrid(i))*(zGrid(j)-zGrid(i))
            ci(1)=+(xGrid(k)-xGrid(j))*(zGrid(l)-zGrid(j))-(xGrid(l)-xGrid(j))*(zGrid(k)-zGrid(j))
            ci(2)=-(xGrid(l)-xGrid(k))*(zGrid(i)-zGrid(k))+(xGrid(i)-xGrid(k))*(zGrid(l)-zGrid(k))
            ci(3)=+(xGrid(i)-xGrid(l))*(zGrid(j)-zGrid(l))-(xGrid(j)-xGrid(l))*(zGrid(i)-zGrid(l))
            ci(4)=-(xGrid(j)-xGrid(i))*(zGrid(k)-zGrid(i))+(xGrid(k)-xGrid(i))*(zGrid(j)-zGrid(i))
            di(1)=-(xGrid(k)-xGrid(j))*(yGrid(l)-yGrid(j))+(xGrid(l)-xGrid(j))*(yGrid(k)-yGrid(j))
            di(2)=+(xGrid(l)-xGrid(k))*(yGrid(i)-yGrid(k))-(xGrid(i)-xGrid(k))*(yGrid(l)-yGrid(k))
            di(3)=-(xGrid(i)-xGrid(l))*(yGrid(j)-yGrid(l))+(xGrid(j)-xGrid(l))*(yGrid(i)-yGrid(l))
            di(4)=+(xGrid(j)-xGrid(i))*(yGrid(k)-yGrid(i))-(xGrid(k)-xGrid(i))*(yGrid(j)-yGrid(i))
 ! coefficient a-> shape_function = a + bx +cy +dz
            ai(1)=xGrid(j)*yGrid(k)*zGrid(l) + xGrid(k)*yGrid(l)*zGrid(j) + xGrid(l)*yGrid(j)*zGrid(k) - &
                 xGrid(l)*yGrid(k)*zGrid(j) - xGrid(j)*yGrid(l)*zGrid(k) - xGrid(k)*yGrid(j)*zGrid(l)
            ai(2)=xGrid(i)*yGrid(l)*zGrid(k) + xGrid(k)*yGrid(i)*zGrid(l) + xGrid(l)*yGrid(k)*zGrid(i) - &
xGrid(i)*yGrid(k)*zGrid(l) - xGrid(k)*yGrid(l)*zGrid(i) - xGrid(l)*yGrid(i)*zGrid(k)
            ai(3)=xGrid(i)*yGrid(j)*zGrid(l) + xGrid(j)*yGrid(l)*zGrid(i) + xGrid(l)*yGrid(i)*zGrid(j) - &
                 xGrid(i)*yGrid(l)*zGrid(j) - xGrid(j)*yGrid(i)*zGrid(l) - xGrid(l)*yGrid(j)*zGrid(i)
            ai(4)=- xGrid(i)*yGrid(j)*zGrid(k) - xGrid(j)*yGrid(k)*zGrid(i) - xGrid(k)*yGrid(i)*zGrid(j) + &
                 xGrid(i)*yGrid(k)*zGrid(j) + xGrid(j)*yGrid(i)*zGrid(k) + xGrid(k)*yGrid(j)*zGrid(i)
            
            Det=(xGrid(l)-xGrid(i))*bi(4)+(yGrid(l)-yGrid(i))*ci(4)+(zGrid(l)-zGrid(i))*di(4)
            VE(iE,iSE)=abs(Det)/6. !absolute value??
            CapE=(Cap(i)+Cap(j)+Cap(k)+Cap(l))/4
            ConE=(Con(i)+Con(j)+Con(k)+Con(l))/4
            AMul=ConE/36./VE(iE,iSE)
            BMul=ConE/6
            FMul=VE(iE,iSE)/20
            SinkE=(sink(i)+sink(j)+sink(k)+sink(l))/4
            DS(i)=DS(i)+FMul*(4*SinkE+sink(i))
            DS(j)=DS(j)+FMul*(4*SinkE+sink(j))
            DS(k)=DS(k)+FMul*(4*SinkE+sink(k))
            DS(l)=DS(l)+FMul*(4*SinkE+sink(l))
            

!flow leaving soil is average sink times their volume  over all the elements
            RootSkold=RootSkold+VE(iE,iSE)*SinkE
             if (tpot.ne.0) then!
                vol_elem_sink = vol_elem_sink + VE(iE,iSE)
             endif

             DO 14 ii=1,4
                iG=elmnod(iE,iL(ii))
                F(iG)=F(iG)+FMul*5
                B(iG)=B(iG)+BMul*(Cxz*bi(ii)+Cyz*ci(ii)+Czz*di(ii))! B1 before!!
                DO 13 jj=1,4
                   jG=elmnod(iE,iL(jj))
                   EE(ii,jj)=Cxx*bi(ii)*bi(jj)+Cyy*ci(ii)*ci(jj)+Czz*di(ii)*di(jj) &
                   +Cxy*(bi(ii)*ci(jj)+bi(jj)*ci(ii))&
                   +Cxz*(bi(ii)*di(jj)+bi(jj)*di(ii)) &
                   +Cyz*(ci(ii)*di(jj)+ci(jj)*di(ii))
                   call Find(iG,jG,kk,nPt,maxbnd,IAD,IADN)
                   ! normal value added to array A_sparse (CSR format)
                   A_sparse(IROW(iG)-1+kk)=A_sparse(IROW(iG)-1+kk)+AMul*EE(ii,jj)
13             CONTINUE
14           CONTINUE
120      CONTINUE
200  CONTINUE

!Determine boundary fluxes

     DO 19 N=1,nPt
          IF (Kode(N).gt.0) then
             QN=B(N)+DS(N)+F(N)*(ThNew(N)-ThOld(N))/dt
             ind=IROW(N)
             ind1=IROW(N+1)-1
             DO j=ind,ind1
                QN=QN+A_sparse(j)*hNew(JCOL(j))
                !if a row is part of a bound.cond. then get the values in A_sparse and multiply by the corresponding column values in hNew
             ENDDO
             Q(N)=QN
          else IF (Kode(N).EQ.-2) then
             Q(N)=-con(N)*Width(N)
          else
             Q(N)=0
          endif
   19 CONTINUE

! Complete construction of RHS vector and form effective matrix:
! use row index - 1 to get the number of previous entries, add up the entry number of the diagonal value to find the location of the diagonal value in A_sparse
      DO 20 i=1,nPt
         A_sparse(IROW(i)-1+IADD(i)) = A_sparse(IROW(i)-1+IADD(i))+ F(i)*Cap(i)/dt !A_sparse(IROW(i)-1+IADD(i)) = A_sparse(IROW(i)-1+IADD(i)) +  F(i)/dt !CSR format!A_sparse(IROW(i)-1+IADD(i)) = A_sparse(IROW(i)-1+IADD(i)) +  F(i)/dt !CSR format
         B(i)=F(i)*Cap(i)*hNew(i)/dt-F(i)*(ThNew(i)-ThOld(i))/dt+ Q(i)-B(i)-DS(i)   !B(i)=Q(i)-B1(i)+hold(i)*F(i)/dt-DS(i)
20   CONTINUE

      RETURN
      END SUBROUTINE ResetCelia
!***********************************************************************

    SUBROUTINE Dirich(Kode,hNew,lOrt,IADD)
      USE GridData
      USE MatData
      Implicit None
      REAL(sp) :: hNew(maxnod)
      INTEGER(sp) :: Kode(maxnod),k,l,m,n,IADD(maxnod),ind,ind1,j
      logical lOrt
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
70    CONTINUE
      RETURN
      END
!*******************************************************************
      SUBROUTINE Solve
      USE GridData, only: nPt,nBand
      USE MatData

      IMPLICIT NONE
      REAL(dp) :: C
      REAL(sp) :: t
      INTEGER(sp) :: i,j,k,l,m,n

! Reduction:
      DO 30 N=1,nPt
         DO 20 M=2,nBand
            IF (ABS(A(M,N)).LT.1.E-30_dp) GOTO 20
            C=A(M,N)/A(1,N)
            I=N+M-1
            IF (I.GT.nPt) GOTO 20
            J=0
            DO 10 K=M,nBand
               J=J+1
10          IF (ABS(A(K,N)).GT.0._dp) A(J,I)=A(J,I)-C*A(K,N)
            A(M,N)=C
            B(I)=B(I)-A(M,N)*B(N)
20       CONTINUE
         B(N)=B(N)/A(1,N)
30    CONTINUE

! Back substitution:
      N=nPt
40    DO 50 K=2,nBand
         L=N+K-1
         IF (L.GT.nPt) GOTO 60
50    IF (ABS(A(K,N)).GT.0._dp) B(N)=B(N)-A(K,N)*B(L)
60    N=N-1
      IF (N.GT.0) GOTO 40
      RETURN
      END
!*****************************************************************************
      SUBROUTINE TmCont(iter,t,dt,dtOpt,tCallR,tFEM,tRoo,dtMaxC,tcBCr,tProf,tProbe)
      USE tmctrl
      IMPLICIT NONE
      REAL(sp), intent(out) :: dt
      REAL(sp), intent(in) :: t,dtMaxC,tCallR,tFEM,tRoo,tcBCr,tProf,tProbe
      REAL(sp), intent(inout) :: dtOpt
      REAL(sp)::tfix
      INTEGER(sp) :: iter

      dtMax=MIN(dtMax,dtMaxC)
      tFix=MIN(tcallr,tFEM,tRoo,tMax,tcBCr,tProf,tProbe)
!print *,tcallr,tFEM,tRoo,tMax,'tbcr',tcBCr,tProf,tProbe,'tfix',tfix
      IF (iter.LE.3.AND.(tFix-t).GE.FacInc*dtOpt)  dtOpt=MIN(dtMax,FacInc*dtOpt)!3
      IF (iter.GE.7)  dtOpt=MAX(dtMin,FacDec*dtOpt) !7
      dt=MIN(dtOpt,tFix-t)
      dt=MIN((tFix-t)/ANINT((tFix-t)/dt),dtMax)
      IF(tFix-t.NE.dt.and.dt.GT.(tFix-t)/2._dp) dt=(tFix-t)/2._dp
print *,iter,'dt=',dt

      RETURN
      END SUBROUTINE TmCont
!****************************************************************
      SUBROUTINE WatInf(Kode,dt)
      USE GridData
      USE CumData
      IMPLICIT NONE
      REAL(sp):: vMean(2),dt
      INTEGER(sp):: Kode(maxnod),i,j
      vMean(1)=0.0_dp
      vMean(2)=0.0_dp
      DO 13 i=1,nPt
         j=iabs(Kode(i))
         wCumA=wCumA+abs(Q(i))*dt
         IF (j.NE.0) THEN
            vMean(j)=vMean(j)-Q(i)
         ENDIF
13    CONTINUE
      wCumA=wCumA+abs(RootSk*dt)
      CumRt=CumRt+RootSk*dt
      VolSink = RootSk*dt
      wCumT=CumRt
	DO 14 j=1,2
         CumQ(j)=CumQ(j)+vMean(j)*dt
         VolQ = vMean(j)*dt
         wCumT=wCumT+CumQ(j)
14    CONTINUE
      RETURN
      END
! ==============================================================================
! Source file SINK |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
      SUBROUTINE BetDis(t)
      USE RootData
      USE PlntData
      USE GridData
      IMPLICIT NONE
      REAL(sp):: betcon(8)
      INTEGER(sp):: corner(8)
      LOGICAL intsec,split
      REAL(sp):: x1,x2,xA,xB,y1,y2,yA,yB,z1,z2,zA,zB
      REAL(sp):: xInt,yInt,zInt,xCent,yCent,zCent
      REAL (sp)::dis,sqrdis,sumB,weight,segage,srface
      REAL(sp):: radius,sgmpl,totlen,prvsur,t
      INTEGER (sp)::i,ibr,iprv,igrow,iseg,iorder,iurf,iface,ic,irec,imin,ipl,isub
      REAL(dp):: blengt
      PrvSur=0.0_dp
      TotSur=0.0_dp
      TotLen=0.0_dp
      split=.FALSE.
      ipl=0
      DO 1 iseg=1,nrec
         TotLen=TotLen+seglen(iseg)
         PrvSur=PrvSur+segsur(iseg)
    1      CONTINUE
      IF (TotLen.LT.1.E-20_dp) TotLen=1.E-20_dp
      IF (PrvSur.LT.1.E-20_dp) PrvSur=1.E-20_dp
          DO 3 i=1,nPt
         betaw(i)=0.0_dp
         betac(i)=0.0_dp
    3 CONTINUE
! now go through each root segment and update the node surface function:
      DO 2 ibr=1,nbr
! find the tip segment of the branch 'ibrnch'
         irec=nrec+1
   11    irec=irec-1
         IF (ibrseg(irec).NE.ibr) GOTO 11
         IF (seglen(irec).LT.1.E-20_dp) THEN
            iprv=irecpr(irec)
            IF (iprv.EQ.0) GOTO 2
            GOTO 101
         ELSE
            DO 34 igrow=1,ngrow
               IF (ibrgrw(igrow).EQ.ibr) THEN
                  iprv=irec
                  IF (iprv.EQ.0) GOTO 2
                  xA=xg(igrow)
                  yA=yg(igrow)
                  zA=zg(igrow)
                  GOTO 102
               ENDIF
   34       CONTINUE
         ENDIF
  100    iprv=irecpr(irec)
         IF (iprv.EQ.0) GOTO 2
  101    IF (ibrseg(iprv).EQ.ibrseg(irec)) THEN
            xA=xs(irec)+xplant(ipl)
            yA=ys(irec)+yplant(ipl)
            zA=zs(irec)
         ELSE
            GOTO 2
         ENDIF
! now find cuboid around segment:
  102    CALL Neighb(xA,yA,zA,corner,imin)
         x1=xgrid(corner(1))
         x2=xgrid(corner(4))
         y1=ygrid(corner(1))
         y2=ygrid(corner(4))
         z1=zgrid(corner(1))
         z2=zgrid(corner(5))
! now find the other end of the segment:
         xB=(xs(iprv))+xplant(ipl)
         yB=(ys(iprv))+yplant(ipl)
         zB=(zs(iprv))
! calculate segment surface:
         sgmpl=segmas(iprv)/seglen(iprv)
         radius=sqrt(sgmpl/SpWgt/pi)
         segsur(iprv)=2*pi*radius*seglen(iprv)
         blengt=seglen(iprv)
         srface=segsur(iprv)
         TotSur=TotSur+segsur(iprv)
! calculate segment weighing factor according to age:
         iorder=0
!  433    iorder=iorder+1
!         IF (ordseg(iprv).NE.iorder) GOTO 433
         iorder=ordseg(iprv)
         segage=t-timorg(iprv)
         IF (segage.GE.0.0_dp) THEN
            IF (segage.GE.age(iorder,nUrf(iorder))) THEN
               Weight=Urf(iorder,nUrf(iorder))
            ELSE
               iUrf=nUrf(iorder)
    4          iUrf=iUrf-1
               IF ((segage.LT.age(iorder,iUrf)).AND.(iUrf.GT.1)) GOTO 4
               Weight=Urf(iorder,iUrf)+(segage-age(iorder,iUrf))/&
                    (age(iorder,iUrf+1)-age(iorder,iUrf))*&
                    (Urf(iorder,iUrf+1)-Urf(iorder,iUrf))
            ENDIF
         ELSE
            Weight=1.0_dp
         ENDIF
         IF (intsec(xA,yA,zA,xB,yB,zB,irec,iprv,isub,ipl,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,iFace)) THEN
            split=.TRUE.
! calculate subsegment length and surface...
   41       blengt=SQRT((xInt-xA)*(xInt-xA)+(yInt-yA)*(yInt-yA)+(zInt-zA)*(zInt-zA))
            srface=blengt*segsur(iprv)/seglen(iprv)
! calculate subsegment center coordinates...
            xCent=xA+(xInt-xA)/2._dp
            yCent=yA+(yInt-yA)/2._dp
            zCent=zA+(zInt-zA)/2._dp
! and calculate the distribution for each node:
            sumB=0.0_dp
            DO 31 ic=1,8
             sqrdis=(xCent-xgrid(corner(ic)))*(xCent-xgrid(corner(ic)))+ &
                  (yCent-ygrid(corner(ic)))*(yCent-ygrid(corner(ic)))+ &
                  (zCent-zgrid(corner(ic)))*(zCent-zgrid(corner(ic)))
             dis=SQRT(sqrdis)
             IF (dis.LT.1.E-20_dp) dis=1.E-20_dp
             betcon(ic)=1._dp/dis
             sumB=sumB+betcon(ic)
   31       CONTINUE
            DO 32 ic=1,8
               betaw(corner(ic))=betaw(corner(ic))+betcon(ic)/sumB*Weight*(blengt/TotLen)
               betac(corner(ic))=betac(corner(ic))+betcon(ic)/sumB*Weight*(srface/PrvSur)
   32       CONTINUE
            xA=xInt
            yA=yInt
            zA=zInt
            IF(iFace.EQ.1) zInt=zInt-1.E-5_dp*dzGrid
            IF(iFace.EQ.2) zInt=zInt+1.E-5_dp*dzGrid
            IF(iFace.EQ.3) yInt=yInt-1.E-5_dp*dyGrid
            IF(iFace.EQ.4) yInt=yInt+1.E-5_dp*dyGrid
            IF(iFace.EQ.5) xInt=xInt-1.E-5_dp*dxGrid
            IF(iFace.EQ.6) xInt=xInt+1.E-5_dp*dxGrid
            CALL Neighb(xInt,yInt,zInt,corner,imin)
!calculate cuboid's corners coordinates:
              x1=xgrid(corner(1))
              x2=xgrid(corner(4))
              y1=ygrid(corner(1))
              y2=ygrid(corner(4))
              z1=zgrid(corner(1))
              z2=zgrid(corner(5))
!MJ1.1
              IF (intsec(xA,yA,zA,xB,yB,zB,irec,iprv,isub,ipl,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,iFace)) THEN
                  GOTO 41
              ELSE
! calculate subsegment length and surface...
                  split=.FALSE.
                  blengt=SQRT((xB-xA)*(xB-xA)+(yB-yA)*(yB-yA)+(zB-zA)*(zB-zA))
                  srface=blengt*segsur(iprv)/seglen(iprv)
              ENDIF
         ENDIF
!calculate segment center coordinates:
         xCent=xA+(xB-xA)/2._dp
         yCent=yA+(yB-yA)/2._dp
         zCent=zA+(zB-zA)/2._dp
!calculate the distribution for each node:
         sumB=0.0
         DO 21 ic=1,8
            sqrdis=(xCent-xgrid(corner(ic)))*(xCent-xgrid(corner(ic)))+&
                 (yCent-ygrid(corner(ic)))*(yCent-ygrid(corner(ic)))+&
                 (zCent-zgrid(corner(ic)))*(zCent-zgrid(corner(ic)))
            dis=SQRT(sqrdis)
            IF (dis.LT.1.E-20_dp) dis=1.E-20_dp
            betcon(ic)=1._dp/dis
            sumB=sumB+betcon(ic)
   21    CONTINUE
         DO 22 ic=1,8
            betaw(corner(ic))=betaw(corner(ic))+betcon(ic)/sumB*Weight*(blengt/TotLen)
            betac(corner(ic))=betac(corner(ic))+betcon(ic)/sumB*Weight*(srface/PrvSur)
   22    CONTINUE
         irec=iprv
         GOTO 100
    2 CONTINUE
      RETURN
      END SUBROUTINE BetDis
!********************************************************************************************
      SUBROUTINE BetNrm
      USE PlntDATA
      USE GridData
      IMPLICIT NONE
      INTEGER(sp):: i,j,k,l,iE,ise
      REAL(sp):: betwe,betce,VEl,Sbetac,Sbetaw
      Sbetaw=0.0_dp
      Sbetac=0.0_dp
      DO 1 iE=1,nElm
       DO 1 iSE=1,3
          i=elmnod(iE,3+iSE)
          j=elmnod(iE,4+iSE-iSE/3*6)
          k=elmnod(iE,5+iSE-iSE/2*6)
          l=elmnod(iE,  iSE)
          betwE=(betaw(i)+betaw(j)+betaw(k)+betaw(l))/4.
          Sbetaw=Sbetaw+betwE
          betcE=(betac(i)+betac(j)+betac(k)+betac(l))/4.
          Sbetac=Sbetac+betcE
    1 CONTINUE
      VEl=dxGrid*dyGrid*dzGrid/6._dp
      IF (Sbetaw.LT.1.E-20_dp) GOTO 100
      Sbetaw=Sbetaw*VEl
      DO 2 i=1,nPt
         betaw(i)=betaw(i)/Sbetaw
2     CONTINUE
100   IF (Sbetac.LT.1.E-20_dp) RETURN
      Sbetac=Sbetac*VEl
      DO 3 i=1,nPt
         betac(i)=betac(i)/Sbetac
3     CONTINUE
      RETURN
      END SUBROUTINE BetNrm
!**************************************************************************************
      LOGICAL FUNCTION intsec(xA,yA,zA,xB,yB,zB,ifoln,irecn,isub,ipl,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,iFace)
      USE Typedef
      USE GridData, only: xgrid,ygrid,dxgrid,dygrid,nex,ney,continu
      USE DoussanMat, only: transroot,nsub
      IMPLICIT NONE
      REAL (sp) ::xA,yA,zA,xB,xBp,yB,yBp,zB,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,xc
      REAL(sp) ::yc,zc,f
      INTEGER(sp) ::iFace,irecn,ifoln,isub,ipl
      intsec=.FALSE.
      IF (continu) THEN										!reconstruction of the continuous B-Ap segment (Couvreur dec 2009)
         IF (isub.EQ.1) THEN
            xBp=xB+(transroot(ifoln,1,nsub(ifoln,ipl),ipl)-transroot(irecn,1,1,ipl))*(nex*dxgrid/2)	!ifoln refers to the "Ap" node and irecn to the "B" node
	     yBp=yB+(transroot(ifoln,2,nsub(ifoln,ipl),ipl)-transroot(irecn,2,1,ipl))*(ney*dygrid)
         ELSE
            xBp=xB+(transroot(irecn,1,isub,ipl)-transroot(irecn,1,1,ipl))*(nex*dxgrid/2)
	     yBp=yB+(transroot(irecn,2,isub,ipl)-transroot(irecn,2,1,ipl))*(ney*dygrid)
         ENDIF
      ELSE
         xBp=xB
         yBp=yB
      ENDIF
      IF (zB.LT.z1) THEN
! may have intersection with x-y-plane of cube at z1
         f=(z1-zA)/(zB-zA)
         xc=xA+f*(xBp-xA)
         yc=yA+f*(yBp-yA)
         IF (((xc.GE.x1).AND.(xc.LE.x2).AND.yc.GE.y1).AND.(yc.LE.y2)) THEN
            xInt=xc
            yInt=yc
            zInt=z1
            iFace=1
            GOTO 1
         ENDIF
      ENDIF
      IF (zB.GT.z2) THEN
! may have intersection with x-y-plane of cube at z2
         f=(z2-zA)/(zB-zA)
         xc=xA+f*(xBp-xA)
         yc=yA+f*(yBp-yA)
         IF ((xc.GE.x1).AND.(xc.LE.x2).AND.(yc.GE.y1).AND.(yc.LE.y2)) THEN
            xInt=xc
            yInt=yc
            zInt=z2
            iFace=2
            GOTO 1
         ENDIF
      ENDIF
      IF (yBp.LT.y1) THEN
! may have intersection with x-z-plane of cube at y1
         f=(y1-yA)/(yBp-yA)
         xc=xA+f*(xBp-xA)
         zc=zA+f*(zB-zA)
         IF ((xc.GE.x1).AND.(xc.LE.x2).AND.(zc.GE.z1).AND.(zc.LE.z2)) THEN
            xInt=xc
            yInt=y1
            zInt=zc
            iFace=3
            GOTO 1
         ENDIF
      ENDIF
      IF (yBp.GT.y2) THEN
! may have intersection with x-z-plane of cube at y2
         f=(y2-yA)/(yBp-yA)
         xc=xA+f*(xBp-xA)
         zc=zA+f*(zB-zA)
         IF ((xc.GE.x1).AND.(xc.LE.x2).AND.(zc.GE.z1).AND.(zc.LE.z2)) THEN
            xInt=xc
            yInt=y2
            zInt=zc
            iFace=4
            GOTO 1
         ENDIF
      ENDIF
      IF (xBp.LT.x1) THEN
! may have intersection with y-z-plane of cube at x1
         f=(x1-xA)/(xBp-xA)
         yc=yA+f*(yBp-yA)
         zc=zA+f*(zB-zA)
         IF ((yc.GE.y1).AND.(yc.LE.y2).AND.(zc.GE.z1).AND.(zc.LE.z2)) THEN
            xInt=x1
            yInt=yc
            zInt=zc
            iFace=5
            GOTO 1
         ENDIF
      ENDIF
      IF (xBp.GT.x2) THEN
! may have intersection with y-z-plane of cube at x2
         f=(x2-xA)/(xBp-xA)
         yc=yA+f*(yBp-yA)
         zc=zA+f*(zB-zA)
         IF ((yc.GE.y1).AND.(yc.LE.y2).AND.(zc.GE.z1).AND.(zc.LE.z2)) THEN
            xInt=x2
            yInt=yc
            zInt=zc
            iFace=6
            GOTO 1
         ENDIF
      ENDIF
      RETURN
    1 intsec=.TRUE.
      RETURN
      END FUNCTION intsec
!===============================================================================
! Doussan model implementation through 3 Subroutines
!	- SetupDou
!       - segment
!	- SetBCroot
!       - conductroot
!	- SolveRoot
! and the adaptation of SetSink
!===============================================================================
      SUBROUTINE SetupDou(t,naxes,level,dt)
      USE TypeDef
      USE SparseMatrix
      USE DoussanMat, only: GH,Khr,Lr,Qi,Q_bc,curr_BCr,no_voxels,voxel_node,numNodes_voxel,voxel_no, &
	  cp_mean,ave,eqDis,old,h_mat2,nlibr,nbc_irecn,nbc_iprvn,loc_Q,curr_bctp,hx_min,Phi_mat,nBCn,inimat, &
         plantmatrix,nplant,transroot,nsub
      USE GridData, only:xgrid,ygrid,zgrid,dzgrid,dygrid,dxgrid,nex,ney,nElm,continu
      USE RootData, only: nrec, nbr, ibrseg, xs, ys, zs, irecpr, ngrow, ibrgrw, xg,yg,zg,segsur,seglen,xplant,yplant
      USE tmctrl, only: tOuRoo,t_begin
      IMPLICIT NONE
      INTEGER(sp) ::ibr,irecn,ifoln,igrow,iprvn,iBCn,naxes!as in DoussaMat
      INTEGER(sp) ::isub,ind,level,i,j,corner(1:8)
      INTEGER :: err,ipl !mutliple roots
      REAL(sp) ::t,xA,yA,zA,xB,yB,zB,x1,x2,y1,y2,z1,z2,dt
      REAL(sp) :: y(1:nLibr),h_mat(1:nLibr),arr(1:nLibr-1),cumsum(1:nLibr-1),PHtemp
      REAL(sp) :: theta(1:nLibr),K(1:nLibr),C(1:nLibr)
      LOGICAL :: n_apex,run

        
      allocate (nBC_irecn(1:naxes))
      nBC_irecn=0

! initialize the pointer of the first el of the list
!* initialize matrices
      CALL IniMat

      iBCn=0
!* Current axial and radial conductance matrix
      CALL ConductRoot(t)

!* go through different plants
      DO ipl=1,nplant
!* Current boundary conditions for root
         CALL SetBCroot(t,curr_BCr(ipl),curr_BCtp(ipl),level)
!multiple roots
         err=SM_allocate(plantmatrix(ipl), nrec+1, nrec+1)
         IF(err/=0) Stop 'Could not create plantmatrix'
!* go through each root segment and update the node surface function:
         DO ibr=1,nbr !all plants have same number of root segments!
            n_apex=.false.
!* find the tip segment of the branch 'ibr'
            irecn=nrec
            DO WHILE (ibrseg(irecn).NE.ibr)
	        irecn=irecn-1
            END DO
!            IF (irecn.eq.0) segsur(irecn)=segsur(irecn+1)
!* the first one we find is an apex
            IF (seglen(irecn)<1.E-20) THEN ! skip this segment too small to be taken
               xA=xs(irecn)+xplant(ipl)
               yA=ys(irecn)+yplant(ipl)
	        IF (continu) THEN
	           DO WHILE (xA.GE.(minval(xGrid)+nex*dxgrid/2))			! node translation if outside of the domain (Couvreur nov 2009)
	              xA=xA-nex*dxgrid/2
	              transroot(irecn,1,1,ipl)=transroot(irecn,1,1,ipl)-1
	           END DO
	           DO WHILE (xA.LT.(minval(xGrid)))
	              xA=xA+nex*dxgrid/2
	              transroot(irecn,1,1,ipl)=transroot(irecn,1,1,ipl)+1
	           END DO
	           DO WHILE (yA.GE.(minval(yGrid)+ney*dygrid))
	              yA=yA-ney*dygrid
	              transroot(irecn,2,1,ipl)=transroot(irecn,2,1,ipl)-1
	           END DO
	           DO WHILE (yA.LT.(minval(yGrid)))
	              yA=yA+ney*dygrid
	              transroot(irecn,2,1,ipl)=transroot(irecn,2,1,ipl)+1
	           END DO
                  DO igrow=1,ngrow
                     IF (ibrgrw(igrow)==ibr) transroot(nrec+igrow,:,1,ipl)=transroot(irecn,:,1,ipl)
		    END DO
                  nsub(irecn,ipl)=1
	        ENDIF
               zA=zs(irecn)
               ifoln=irecn ! following node ID
               irecn=irecpr(irecn) !current node ID
            ELSE
               n_apex=.true.!ifoln does not exist if not continu
               DO igrow=1,ngrow
                  IF (ibrgrw(igrow)==ibr) THEN !apex of branch ibr
                     xA=xg(igrow)+xplant(ipl)
                     yA=yg(igrow)+yplant(ipl)
		       IF (continu) THEN
	                 DO WHILE (xA.GE.(minval(xGrid)+nex*dxgrid/2))				! Translate the tip of n "domain length" if it is outside of the domain (Couvreur nov 2009)
	                    xA=xA-nex*dxgrid/2
	                    transroot(nrec+igrow,1,1,ipl)=transroot(nrec+igrow,1,1,ipl)-1
	                 END DO
	                 DO WHILE (xA.LT.(minval(xGrid)))
	                    xA=xA+nex*dxgrid/2
	                    transroot(nrec+igrow,1,1,ipl)=transroot(nrec+igrow,1,1,ipl)+1
	                 END DO
	                 DO WHILE (yA.GE.(minval(yGrid)+ney*dygrid))
	                    yA=yA-ney*dygrid
	                    transroot(nrec+igrow,2,1,ipl)=transroot(nrec+igrow,2,1,ipl)-1
	                 END DO
	                 DO WHILE (yA.LT.(minval(yGrid)))
	                    yA=yA+ney*dygrid
	                    transroot(nrec+igrow,2,1,ipl)=transroot(nrec+igrow,2,1,ipl)+1
	                 END DO
	                 ifoln=nrec+igrow                              ! Identify the row of transroot in which the eventual translations of the tip are recorded (Couvreur nov 2009)
                        nsub(ifoln,ipl)=1
		       ENDIF
                     zA=zg(igrow)
                  ENDIF
               END DO
            ENDIF

	     IF (irecn==0) THEN !there exists a branch ibr but not yet any segment!
	        run=.false.
	     ELSE
               run=.true.
	     ENDIF

!* then the rest of the branch up to the seed or the embranchment
            DO WHILE (run)
!* "upper" node
	        iprvn=irecpr(irecn)
!* location of the rear or "upper" end
               xB=xs(irecn)+xplant(ipl)
               yB=ys(irecn)+yplant(ipl)
	        IF (continu) THEN
	           DO WHILE (xB.GE.(minval(xGrid)+nex*dxgrid/2))	                ! Translate the segment of n "domain length" if it is outside of the domain (Couvreur nov 2009)
	              xB=xB-nex*dxgrid/2
     	              transroot(irecn,1,1,ipl)=transroot(irecn,1,1,ipl)-1
	           END DO
	           DO WHILE (xB.LT.(minval(xGrid)))
	              xB=xB+nex*dxgrid/2
	              transroot(irecn,1,1,ipl)=transroot(irecn,1,1,ipl)+1
	           END DO
	           DO WHILE (yB.GE.(minval(yGrid)+ney*dygrid))
	              yB=yB-ney*dygrid
	              transroot(irecn,2,1,ipl)=transroot(irecn,2,1,ipl)-1
	           END DO
	           DO WHILE (yB.LT.(minval(yGrid)))
	              yB=yB+ney*dygrid
	              transroot(irecn,2,1,ipl)=transroot(irecn,2,1,ipl)+1
	           END DO
	        ENDIF
               zB=zs(irecn)
!* calculate the gravity components z (always positive & maximum at soil surface)
               GH(irecn)=(zA+zB)/2.
!* calculate number of subsegment for node/segment irecn and corresponding weights
               CALL segment(xA,yA,zA,xB,yB,zB,ifoln,irecn,ipl)

! several changes -> multiple roots
	        err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,-Khr(irecn)/seglen(irecn),.false.)
               If(err/=0) Stop 'Could not insert element into plantmatrix'
!if apex (bottom part of root)
               IF (n_apex) THEN
                  Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn)/seglen(irecn)+Lr(irecn)*segsur(irecn),.false.)
                  If(err/=0) Stop 'Could not insert element into plantmatrix'
               ELSE
 	           Err=SM_set(Plantmatrix(ipl), irecn+1,ifoln+1,-Khr(ifoln)/seglen(ifoln),.false.) !row, col,value
                  Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn)/seglen(irecn)+Khr(ifoln)/seglen(ifoln)+Lr(irecn)*segsur(irecn),.false.)
                  If(err/=0) Stop 'Could not insert element into plantmatrix'
               ENDIF

! define 1st part of Q (Q=Qi.*PHsoil+Qbc) -> RHS
	        Qi(irecn,ipl)=Lr(irecn)*segsur(irecn)
! if reached the seed or the embranchement => change branch
               IF (iprvn==0) THEN!seed=first segment
	           GH(0)=zB
                  IF (curr_BCtp(ipl)==2) THEN!flux
	              Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn)/seglen(irecn),.false.)
		       Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn)/seglen(irecn),.false.)
                     If(err/=0) Stop 'Could not insert element into plantmatrix'
!own additive; position (2,2) changes -> what mathieu did is correct, only mathematically wrongly written!!!
!                     err=SM_set(plantmatrix(ipl), irecn+1,irecn+1,Khr(ifoln)/seglen(ifoln)+Lr(irecn)*segsur(irecn),.true.)
!own additive;position (2,1) = 0, changes into zero
!		        Err=SM_set(Plantmatrix(ipl), irecn+1,iprvn+1,0._dp,.true.)
!rhs has only 1 entry (first one)
                     Qi(iprvn,ipl)=0._dp
		       Q_bc(iprvn,ipl)=curr_BCr(ipl);!Q(0,0)
!                     Q_bc(irecn,ipl)=curr_BCr
                  ELSE IF (curr_BCtp(ipl)==1) THEN!PH
!iprvn+1,iprvn+1
	              Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,1._dp,.true.)!position (1,1)
                     Qi(iprvn,ipl)=0._dp
! true means overwrite (delete first), no addition
		       Err=SM_set(Plantmatrix(ipl), irecn+1,iprvn+1,0._dp,.true.)!irecn+1,iprvn+1
                     If(err/=0) Stop 'Could not insert element into plantmatrix'
!first entry in rhs is PH (inlc. gravity)
		       Q_bc(iprvn,ipl)=curr_BCr(ipl)+GH(iprvn)!Q(iprv)-BCr*Kh(irecn)/seglen(irecn);
!second entry is PH incl. gravity times these parameters
		       Q_bc(irecn,ipl)=(curr_BCr(ipl)+GH(iprvn))*Khr(irecn)/seglen(irecn)
                  ENDIF
                  iBCn=iBCn+1
	           nBC_irecn(iBCn)=irecn
                  run=.false.
      	           CALL segment(xB,yB,zB,xB,yB,zB,iprvn,iprvn,ipl) !MJ dec08 this line is there what for????
               ELSEIF (ibrseg(iprvn).NE.ibrseg(irecn)) THEN!start of the branch but not from the seed -> gaat van onder na boven, duz iprvn (nieuwe positie) is nu niet van die ene branch maar van een side branch
                  Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn)/seglen(irecn),.false.)!iprvn+1,iprvn+1
                  Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn)/seglen(irecn),.false.)!iprvn+1,irecn+1
                  If(err/=0) Stop 'Could not insert element into plantmatrix'
                  run=.false.
               ENDIF

! definition for the next run of the loop
               ifoln=irecn
               irecn=iprvn
!* location of the final node
               xA=xB				! A directly redefined with report to B in order not to have to calculate all the translations of A to the inside of the domain (Couvreur nov 2009)
               yA=yB 
               zA=zB
!* from here, not an apex
               n_apex=.false.
            END DO ! loop on branch nodes
         END DO ! loop on root branches
     

         if (.not.(old)) then
            IF ((ave) .or. (eqDis)) THEN
               DO i=1,nElm !no_voxels -> number of nodes in half a cubiod (imin)
                  if (no_voxels(i) .eq. 0) goto 40
                  DO j=1,no_voxels(i)
                     irecn = voxel_node(i,1,j)
                     isub = voxel_node(i,2,j)
                     corner=loc_Q(irecn,1:8,isub,ipl)
                     !coordinates of voxel
                     x1=xgrid(corner(1))
                     y1=ygrid(corner(1))
                     z1=zgrid(corner(1))
                     x2=x1+dxgrid						!Defining voxel position with report to corner(4),etc. could create problems in case of continuous domain (Couvreur dec 2009)
                     y2=y1+dygrid
                     z2=z1+dzgrid
                     if (ave) then
                        cp_mean(irecn,1,isub,ipl)=(x2+x1)/2
                        cp_mean(irecn,2,isub,ipl)=(y2+y1)/2
                        cp_mean(irecn,3,isub,ipl)=(z2+z1)/2
                     elseif (eqDis) then
                        numNodes_voxel(irecn,isub)=no_voxels(i) !total root nodes of cuboid in no_voxels; half cubiods have same imin
                     endif
                  ENDDO
       40         continue
               ENDDO
            ENDIF
            write(*,*)'t=',t
            write(*,*)'t_begin=',t_begin,t_begin+dt
            write(*,*)'old=',old
!if (.not.(old)) then
            if (t .le. t_begin+dt) then
!call cpu_time(t0)
               write(*,*)'Matric flux potential Library loaded'
!-----------------------------------------------------------------------
!Matric flux potential lookup table -> to obtain corresponding PH
!------------------------------------------------------------------------
!               if (w_sub(irecn,isub,ipl).eq.0) then
               PHtemp = 1e-5
!               else
!                  PHtemp=PH_root_sub(irecn,isub)+50
!               endif
               call logspace(hx_min,PHtemp,nLibr,y)
               h_mat = -y
               h_mat2=(h_mat(2:size(h_mat))+h_mat(1:size(h_mat)-1))/2
               call setmat_anaLibr(h_mat,theta,K,C)
! calculate matric flux potential integral K(h)dh -> numerical integration -> to linearize non-linear conductivity
!               cumsum = cumulitive sum
               arr=abs(h_mat(2:size(h_mat))-h_mat(1:size(h_mat)-1))*(K(2:size(K))+K(1:size(K)-1))/2
!$OMP PARALLEL DO
               do i=1,size(arr)
                  cumsum(i) = sum(arr(1:i))
               end do
!$OMP END PARALLEL DO
!               write(*,*)'cumsum',cumsum(1:100)
!               pause
               Phi_mat = cumsum
!               call cpu_time(t1)
!               write(*,*)'simutime',t1-t0
            endif
         endif
!* total number of embrach. to the seed
         nBCn=iBCn
      END DO ! loop on plants
! save transroot for the basic root nodes -> draw the root in matlab (Couvreur jan 2010)
      IF (continu) THEN
         OPEN (UNIT=444,FILE='out/Transroot.out',STATUS='UNKNOWN')
         WRITE (444,'(/''Basic root nodes translations for each plant.'')')
         DO ipl=1,nplant
         WRITE (444,'(/'' '')')
         WRITE (444,'(/''Plant '',i1)') ipl
         WRITE (444,'(/'' X  Y (times)'')')
            DO i=1,nrec+ngrow
               WRITE (444,'(i3, i3)') transroot(i,1,nsub(i,ipl),ipl), transroot(i,2,nsub(i,ipl),ipl)
            END DO
            WRITE (444,'(/'' '')')
            WRITE (444,'(/'' '')')
         END DO
         CLOSE(444)
      ENDIF
      END SUBROUTINE SetupDou
!*******************************************************************************
      SUBROUTINE segment(xA,yA,zA,xB,yB,zB,ifoln,irecn,ipl)
      USE Typedef
      USE GridData, only:xgrid,ygrid,zgrid,dzgrid,dygrid,dxgrid,nex,ney,continu
      USE DoussanMat, only: nsub,loc_Q,w_sub,sum_dis,w_dis,cent,Intc,no_voxels,voxel_no,voxel_node,indexValue,old,ave,eqdis,l_seg,l_sub,isubmax,transroot
      USE RootData, only: seglen,ibrseg,irecpr
      IMPLICIT NONE
      REAL(sp),intent(in) :: xA,yA,zA,xB,yB,zB
      INTEGER(sp),intent(in) :: ipl
      REAL(sp)::xIni,yIni,zIni,x1,x2,y1,y2,z1,z2,xAp,yAp,zAp
      REAL(sp) ::yint,xint,zint,xCent,yCent,zCent,sqrdis,dis
      REAL(dp)::blengt,blengtTot
      LOGICAL :: SPLIT,intsec
      INTEGER(sp) ::iface,iFaceOld,isub,corner(8),ic,ifoln,irecn,i,imin,reftransroot(2)

! number of subsegment for node irecn
      nsub(irecn,ipl)=0
      l_seg(irecn)=0
      SPLIT=.TRUE.
! initialization
      xAp=xA
      yAp=yA
      zAp=zA
      blengtTot=0
      reftransroot=transroot(ifoln,:,nsub(ifoln,ipl),ipl)			!in intc, transroot, etc., the node from rootsys (-> not created by instec) is written in the last columns
      DO WHILE(SPLIT)
         nsub(irecn,ipl)=nsub(irecn,ipl)+1
         isub=nsub(irecn,ipl)
         IF (isub.GT.isubmax) THEN
            print *,'Number of subsegments higher than maximum admitted, segment',irecn,'too long with report to the grid resolution'	!(Couvreur dec 2009)
            STOP
         ENDIF
!  find cuboid around the apical node of the segment:
         CALL Neighb(xAp,yAp,zAp,corner,imin)

         loc_Q(irecn,1:8,isub,ipl)=corner
         if (.not.(old)) then
            if ((ave) .or. (eqdis)) then
               if (voxel_no(irecn,isub) .eq. 0) then !remove double values from different branchings; if this root node is not processed process it
                  voxel_no(irecn,isub)=imin ! voxel number
                  no_voxels(imin)=no_voxels(imin)+1 !number of root nodes in a voxel
                  voxel_node(imin,1,no_voxels(imin))=irecn
                  voxel_node(imin,2,no_voxels(imin))=isub
                  if (no_voxels(imin) .gt. indexValue) then
                     print*,'Parameter indexValue in Segment, representing number of nodes in a voxel, has to be set larger'
                     STOP
                  endif
               endif
            endif
         endif

! calculate cuboid's corners coordinates:
         x1=xgrid(corner(1))
         y1=ygrid(corner(1))
         z1=zgrid(corner(1))
         x2=x1+dxgrid						!Valid in all cases (continu or not)  (Couvreur dec 2009)
         y2=y1+dygrid
         z2=z1+dzgrid

! check is segment in completely included in the cuboid
         IF (intsec(xAp,yAp,zAp,xB,yB,zB,ifoln,irecn,isub,ipl,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,iFace)) THEN !true when the root segment intersects the cuboïd
            SPLIT=.TRUE.
	  ELSE
            SPLIT=.FALSE.!just/still 1 loop and then no subsegmentation
            xInt=xB
            yInt=yB
            zInt=zB
         ENDIF

!correction of Ap position
         IF (isub.GT.1) THEN
	     IF(iFaceOld.EQ.1) zAp=zAp+5.E-5_dp*dzGrid
            IF(iFaceOld.EQ.2) zAp=zAp-5.E-5_dp*dzGrid
            IF(iFaceOld.EQ.3) yAp=yAp+5.E-5_dp*dyGrid
            IF(iFaceOld.EQ.4) yAp=yAp-5.E-5_dp*dyGrid
            IF(iFaceOld.EQ.5) xAp=xAp+5.E-5_dp*dxGrid
            IF(iFaceOld.EQ.6) xAp=xAp-5.E-5_dp*dxGrid
         ELSEIF (continu) THEN					!(Couvreur dec 2009)
            Intc(ifoln,1,nsub(ifoln,ipl),ipl)=xAp
            Intc(ifoln,2,nsub(ifoln,ipl),ipl)=yAp
            Intc(ifoln,3,nsub(ifoln,ipl),ipl)=zAp
	  ENDIF

! calculate (sub-)segment length
         blengt=SQRT((xInt-xAp)*(xInt-xAp)+(yInt-yAp)*(yInt-yAp)+ (zInt-zAp)*(zInt-zAp))
         blengtTot=blengtTot+blengt

! calculate relative length of this (sub)segment
         if (irecn.eq.0) then
            w_sub(irecn,isub,ipl)=blengt/seglen(irecn+1)
         else
            w_sub(irecn,isub,ipl)=blengt/seglen(irecn)!blengt is sp and seglen is dp...
         endif
         IF ((.not.(split)).and.(isub.eq.1.)) THEN
	     w_sub(irecn,isub,ipl)=1._dp
         ENDIF
         IF ((w_sub(irecn,isub,ipl).GT.(1.-1.E-7)).or.((.not.(split)).and.(isub.eq.1.))) THEN
 	     w_sub(irecn,isub,ipl)=1._dp
 	  ELSEIF (w_sub(irecn,isub,ipl).LT.1.E-7) THEN
 	     w_sub(irecn,isub,ipl)=0._dp
 	  ENDIF
         l_sub(irecn,isub,ipl)=blengt
         l_seg(irecn)=l_seg(irecn)+blengt

! calculate (sub)segment center coordinates...
         xCent=xAp+(xInt-xAp)/2.
         yCent=yAp+(yInt-yAp)/2.
         zCent=zAp+(zInt-zAp)/2.
         cent(irecn,1,isub,ipl)=xCent
         cent(irecn,2,isub,ipl)=yCent
         cent(irecn,3,isub,ipl)=zCent

! and calculate the distribution for each node:
         sum_dis(irecn,isub,ipl)=0.0_dp
         DO ic=0,1 							!Valid in all cases (continu or not)  (Couvreur dec 2009)
            w_dis(irecn,4*ic+1,isub,ipl)=1./SQRT((xCent-xgrid(corner(4*ic+1)))**2+(yCent-ygrid(corner(4*ic+1)))**2+(zCent-zgrid(corner(4*ic+1)))**2)
            w_dis(irecn,4*ic+2,isub,ipl)=1./SQRT((xCent-xgrid(corner(4*ic+1))-dxgrid)**2+(yCent-ygrid(corner(4*ic+1)))**2+(zCent-zgrid(corner(4*ic+1)))**2)
            w_dis(irecn,4*ic+3,isub,ipl)=1./SQRT((xCent-xgrid(corner(4*ic+1)))**2+(yCent-ygrid(corner(4*ic+1))-dygrid)**2+(zCent-zgrid(corner(4*ic+1)))**2)
            w_dis(irecn,4*ic+4,isub,ipl)=1./SQRT((xCent-xgrid(corner(4*ic+1))-dxgrid)**2+(yCent-ygrid(corner(4*ic+1))-dygrid)**2+(zCent-zgrid(corner(4*ic+1)))**2)
         END DO
         DO ic=1,8
            IF (w_dis(irecn,ic,isub,ipl).GT.1.E+20) w_dis(irecn,ic,isub,ipl)=1.E+20
            sum_dis(irecn,isub,ipl)=sum_dis(irecn,isub,ipl)+w_dis(irecn,ic,isub,ipl)
         END DO

!         DO ic=1,8
!             sqrdis=(xCent-xgrid(corner(ic)))*(xCent-xgrid(corner(ic)))+&
!                    (yCent-ygrid(corner(ic)))*(yCent-ygrid(corner(ic)))+&
!                    (zCent-zgrid(corner(ic)))*(zCent-zgrid(corner(ic)))
!             dis=SQRT(sqrdis)
!             IF (dis.LT.1.E-20) dis=1.E-20
!             w_dis(irecn,ic,isub,ipl)=1./dis
!             sum_dis(irecn,isub,ipl)=sum_dis(irecn,isub,ipl)+w_dis(irecn,ic,isub,ipl)
!!print *,'4085',irecn, isub,ic !problem here when maxZ is imposed at line 2064
!         END DO

! save Ap position (inside of the soil domain -> xmax and ymax of the soil boudaries not included) (Couvreur dec 2009)
         IF (continu.AND.(isub.GT.1)) THEN
            IF (xAp.EQ.minval(xGrid)+nex*dxgrid/2) THEN
               xAp=xAp-nex*dxgrid/2
               transroot(irecn,1,isub,ipl)=transroot(irecn,1,isub,ipl)-1
            ENDIF
            IF (yAp.EQ.minval(yGrid)+ney*dygrid) THEN
               yAp=yAp-ney*dygrid
               transroot(irecn,2,isub,ipl)=transroot(irecn,2,isub,ipl)-1
            ENDIF
            Intc(irecn,1,isub-1,ipl)=xAp
            Intc(irecn,2,isub-1,ipl)=yAp
            Intc(irecn,3,isub-1,ipl)=zAp
         ENDIF

         IF (SPLIT) THEN
! preparation for next loop
            xAp=xInt
            yAp=yInt
            zAp=zInt
	     iFaceOld=iFace
            IF(iFace.EQ.1) THEN
               zAp=zInt-5.E-5_dp*dzGrid
               IF(continu) transroot(irecn,:,isub+1,ipl)=reftransroot			!(Couvreur dec 2009)
            ELSEIF(iFace.EQ.2) THEN
               zAp=zInt+5.E-5_dp*dzGrid
               IF(continu) transroot(irecn,:,isub+1,ipl)=reftransroot			!(Couvreur dec 2009)
            ELSEIF(iFace.EQ.3) THEN
	        yAp=yInt-5.E-5_dp*dyGrid
		 IF (continu) THEN								!(Couvreur dec 2009)
                  IF (yAp.LT.minval(ygrid)) THEN
                     yAp=yAp+ney*dygrid
                     transroot(irecn,1,isub+1,ipl)=reftransroot(1)
                     transroot(irecn,2,isub+1,ipl)=reftransroot(2)+1
                     reftransroot=transroot(irecn,:,isub+1,ipl)
                  ELSE
                     transroot(irecn,:,isub+1,ipl)=reftransroot
                  ENDIF
               ENDIF
            ELSEIF(iFace.EQ.4) THEN
	        yAp=yInt+5.E-5_dp*dyGrid
		 IF (continu) THEN								!(Couvreur dec 2009)
                  IF (yAp.GT.minval(ygrid)+ney*dygrid) THEN
                     yAp=yAp-ney*dygrid
                     transroot(irecn,1,isub+1,ipl)=reftransroot(1)
                     transroot(irecn,2,isub+1,ipl)=reftransroot(2)-1
                     reftransroot=transroot(irecn,:,isub+1,ipl)
                  ELSE
                     transroot(irecn,:,isub+1,ipl)=reftransroot
                  ENDIF
               ENDIF
            ELSEIF(iFace.EQ.5) THEN
		 xAp=xInt-5.E-5_dp*dxGrid
		 IF (continu) THEN								!(Couvreur dec 2009)
                  IF (xAp.LT.minval(xgrid)) THEN
                     xAp=xAp+nex*dxgrid/2
                     transroot(irecn,1,isub+1,ipl)=reftransroot(1)+1
                     transroot(irecn,2,isub+1,ipl)=reftransroot(2)
                     reftransroot=transroot(irecn,:,isub+1,ipl)
                  ELSE
                     transroot(irecn,:,isub+1,ipl)=reftransroot
                  ENDIF
               ENDIF
            ELSEIF(iFace.EQ.6) THEN
		 xAp=xInt+5.E-5_dp*dxGrid
		 IF (continu) THEN								!(Couvreur dec 2009)
                  IF (xAp.GT.minval(xgrid)+nex*dxgrid/2) THEN
                     xAp=xAp-nex*dxgrid/2
                     transroot(irecn,1,isub+1,ipl)=reftransroot(1)-1
                     transroot(irecn,2,isub+1,ipl)=reftransroot(2)
                     reftransroot=transroot(irecn,:,isub+1,ipl)
                  ELSE
                     transroot(irecn,:,isub+1,ipl)=reftransroot
                  ENDIF
               ENDIF
	     ENDIF
         ELSEIF (continu.AND.(isub.GT.1)) THEN						!sort transroot in the same order as intc (Couvreur dec 2009)
            transroot(irecn,:,isub+1,ipl)=transroot(irecn,:,1,ipl)
            transroot(irecn,:,1:isub+1,ipl)=transroot(irecn,:,2:isub+2,ipl)
         ENDIF
      END DO
! Ap was on a face of the cube containig A. If Ap was on a boundary too, it is displaced out of the domain (because of the +/-5.E-5...).
! Then -> translation. Ap is now closer to B (on the same side and inside of the domain). (Couvreur dec 2009)

!save brench basis position									!(Couvreur dec 2009)
      IF (continu.AND.(ibrseg(irecpr(irecn)).NE.ibrseg(irecn))) THEN
         Intc(irecn,1,isub,ipl)=xB
         Intc(irecn,2,isub,ipl)=yB
         Intc(irecn,3,isub,ipl)=zB
      ENDIF
! correction of inacurrate seglen
      IF (nsub(irecn,ipl).EQ.1) THEN !nosplit
         w_sub(irecn,1,ipl)=1
      ELSE
         DO isub=1,nsub(irecn,ipl)
           w_sub(irecn,isub,ipl)=l_sub(irecn,isub,ipl)/l_seg(irecn)
         END DO
      ENDIF
      RETURN
      END SUBROUTINE segment
!****************************************************************
      SUBROUTINE SetBCroot(t,BCr,BCtp,level)
! estimate current BC for root system when using Doussan Sink term
      USE TypeDef
      USE PlntData, only : tBCr,typeBCr,BCroot,nBCr,Tpot
      IMPLICIT NONE
      REAL(sp), intent(in) :: t
      REAL(dp), intent(out) ::BCr
      INTEGER, intent(in)::level
      INTEGER (sp),intent(out) :: BCtp
      INTEGER (sp) :: ifc

      IF (level==4) THEN
! calculate current root UBC-value from input BC(time)-function:
         ifc=0
  201    ifc=ifc+1
         IF (ifc.GT.nBCr) THEN !beyond the last user-defined BC value it remains constant
            BCr=BCroot(nBCr)
	     BCtp=typeBCr(nBCr)
         ELSE
           IF (t.GE.tBCr(ifc)) GOTO 201 !try to find the first tBCR wwhich is beyond current t (>=1)
	    IF (ifc.EQ.1) THEN
	       BCr=BCroot(ifc)
              BCtp=typeBCr(ifc)
	    ELSEIF ((typeBCr(ifc-1)).EQ.typeBCr(ifc)) THEN
              BCr=BCroot(ifc-1)+(BCroot(ifc)-BCroot(ifc-1))*(t-tBCr(ifc-1))/(tBCr(ifc)-tBCr(ifc-1))
              BCtp=typeBCr(ifc-1)
	    ELSE
              BCr=BCroot(ifc-1)
              BCtp=typeBCr(ifc-1)
	    ENDIF
         ENDIF
      ELSE !lDou but with level 2 or 3
         BCr=Tpot
         BCtp=2
      ENDIF

      IF (BCtp.EQ.2) THEN
         BCr=-abs(BCr) !Javaux always negative flow at the collar
      ENDIF
      RETURN
      END SUBROUTINE SetBCroot
!**************************************************************************
      SUBROUTINE Rootstress(PHtop,BCtop,BCtptop,iterBC,BC_switch,Jcol,ipl)
      USE Typedef
      USE DoussanMat, only: curr_BCtp,stressBC,hx_min,stresfun,stresval1,stresval2
      USE GridData, only: RelEps,epslonR,factorRelEps

      IMPLICIT NONE
      REAL(dp), intent(inout):: BCtop
      INTEGER, intent(inout):: BCtptop
      INTEGER, intent(in):: ipl
      REAL(dp), intent(in):: PHtop,Jcol
      LOGICAL, intent(out) :: BC_switch,iterBC
      REAL(dp):: reduc,del_PHr

!check if the root collar abs(PH) is larger than abs(hx_min)+tolerance and adapt the collar BC
      IF (curr_BCtp(ipl)==2) THEN
         IF (RelEps) THEN 
            del_PHr=-max(abs(hx_min/factorRelEps),epslonR)
         ELSE
            del_PHr=-epslonR
         ENDIF
         IF ((stresfun.EQ.1).and.(PHtop<hx_min+del_PHr)) THEN
!top node at lower PH than allowed: start of stressed conditions
            print *,'stress in the collar xylem: change to PH BC, PHtop=',PHtop,' is lower than criterion ',hx_min,' +',del_PHr
	     stressBC=.true.
	     BCtptop=1
	     BCtop=hx_min
            BC_switch=.TRUE.
            iterBC=.TRUE.
         ELSEIF ((stresfun.EQ.2).and.(PHtop.LE.stresval1)) THEN !linear decrease
		 !to be checked.
            reduc=(PHtop-stresval2)/(stresval1-stresval2)
            IF (reduc.GT.1) THEN
               reduc=1
            ELSEIF (reduc.LT.0) THEN
               reduc=0
            ENDIF
            BCtop=Jcol*reduc
            BCtptop=2
            stressBC=.true.
            iterBC=.TRUE.
            print *,'flux reduction of factor',reduc
         ELSEIF (stresfun.EQ.3) THEN !Tuzet function
            reduc=(1+exp(stresval1*stresval2))/(1+exp(stresval1*(stresval2-PHtop)))
            BCtop=Jcol*reduc
            BCtptop=2
            print *,'flux reduction of factor',reduc
            IF (reduc.LT.1) THEN
               BC_switch=.TRUE.
               iterBC=.TRUE.
            ELSE
               BC_switch=.FALSE.
               iterBC=.FALSE.
            ENDIF
         ELSE
            BC_switch=.FALSE.
            iterBC=.FALSE.
         ENDIF
      ELSEIF ((curr_BCtp(ipl)==1).and.(abs(Jcol)>1.00001*abs(BCtop))) THEN !factor to avoid oscillation between stress and no stress
!check if end of stress (curr_BCtp is PH with PH is hxmin but BCtp_new is still flux)
         print *,'end of stress conditions, collar flux=',Jcol,' is larger than the prescribed flux', BCtop
!	  BCtptop=BCnew is already of type 2!
         stressBC=.FALSE.
         BC_switch=.FALSE.
         print *,'info',stressBC,BCtptop,BCtop, BC_switch,iterBC
      ELSE
         BC_switch=.FALSE.
         iterBC=.FALSE.
         IF (stressBC) print *,'no change: stress=',stressBC
      ENDIF
      END SUBROUTINE Rootstress
!*************************************************************
      SUBROUTINE ConductRoot(t)
      USE TypeDef
      USE DoussanMat, only : nLr,Lr,Lrroot,ageLr,nKh,Khr,Khroot,ageKh
      USE RootData, only: nrec,timorg,rrt,age,ordseg
      IMPLICIT NONE
      REAL(sp), intent (in) ::t
      REAL(sp)::segage
      INTEGER(sp)::irecn,iage,typ
! calculate segment lateral and long. conductivity according to age:
      DO irecn=1,nrec
         segage=t-timorg(irecn)
         typ=ordseg(irecn) !-1 !wrong type definition in the btable generated from tootyp
!         write(*,*)'timorg',timorg(irecn)
!         write(*,*)'typ=',typ
!         write(*,*)'irecn',irecn
!pause
	  if (typ.EQ.0) THEN  !type 0 is the seed+small segment
	     typ=1
         elseif (typ.GT.3) THEN!no more than 3 root types
	     typ=3
	  ENDIF
! radial conductivity matrix
         iage=1
	  DO WHILE (ageLr(typ,iage)<=segage)
            iage=iage+1
	  ENDDO
         IF (iage>nLr(typ)) THEN
	     Lr(irecn)=LrRoot(typ,nLr(typ))
         ELSEIF (iage==1) THEN
	     Lr(irecn)=LrRoot(typ,iage)
	  ELSE
            Lr(irecn)=LrRoot(typ,iage-1)+(LrRoot(typ,iage)-LrRoot(typ,iage-1))*(segage-ageLr(typ,iage-1))/&
	               (ageLr(typ,iage)-ageLr(typ,iage-1))
         ENDIF
!	   if (irecn>=185) then
!	      print *,'CONDUCTROOT',iage,segage,ageLr(iage),irecn,Lr(irecn)
!	   endif
! axial conductance matrix
         iage=1
	  DO WHILE (ageKh(typ,iage)<=segage)
            iage=iage+1
	  ENDDO
         IF (iage>nKh(typ)) THEN
	     Khr(irecn)=KhRoot(typ,nKh(typ))
	  ELSEIF (iage==1) THEN
	     Khr(irecn)=Khroot(typ,iage)
	  ELSE
            Khr(irecn)=KhRoot(typ,iage-1)+(KhRoot(typ,iage)-KhRoot(typ,iage-1))*(segage-ageKh(typ,iage-1))/&
	                (ageKh(typ,iage)-ageKh(typ,iage-1))
         ENDIF
!	   if (irecn>=185) then
!	      print *,'CONDUCTROOT2',iage,segage,irecn,Khr(irecn)
!	      stop
!	   endif
      END DO
      END SUBROUTINE ConductRoot
!******************************************************************************
      SUBROUTINE SolveRoot(hNew,hold,t,dt,it1,level,iter)
!calculate sink for lDou simulations
!solve DOUSSAN flow equations within the xylem
      USE TypeDef
      USE SparseMatrix
      USE NumericalRecipes
      USE Paramdata,only: MaxNod, pi
      USE GridData, only: xgrid,ygrid,zgrid, sink, betaw,dxgrid,dygrid,dzgrid,Wn,elmnod,nPt,nElm,RootSk,Axy,Bxy,Dxy,Exy,nex,ney,continu
      USE PlntData, only :Tpot,Tact
      USE DoussanMat, only: stressBC,Q_bc,Qi,PHs,Qd,PHr,GH,w_dis,nsub,sinkR,Jintot,nrecOld,loc_Q,nBCn,bcr_usr,counter2, &
      bctp_usr,curr_bcr,curr_bctp,Khr,w_sub,NBC_irecn,sum_dis,cent,Intc,Lr,ave,old,axialRootFlow,Joutr,PH_root_sub, &
      cp_mean,oldT,hx_min,plantmatrix,nplant !sa,ija,
      USE RootData, only: nrec,seglen,segsur,nbr,ibrseg,irecpr
      USE SolData
      USE tmctrl, only: t_begin
      IMPLICIT NONE

      REAL,intent(in) :: hNew(maxnod),hold(maxnod),t,dt
      LOGICAL,intent(in) :: it1
      INTEGER, intent(in):: level
      REAL (dp)::PHr_tot(1:nrec+1), err_loc,Qd2(1:nrec+1),dVolSk,Voltot,checksub
      REAL (dp)::maxerr=1.e-10_dp,BCr_new,BCr_old,Jcol,sumW,psub
      REAL (sp) :: K_bar,K_corner,K_interface,theta_single,K_single,C_single
      REAL (sp) :: partFlow(1:8),sumFraction,PHtemp,totalAxialFlow=0
      INTEGER(sp) ::old_BCtp,iter_loc,itol,itmaxl,iter,corner_ic,ibr
      INTEGER(sp) :: iprvn,irecn,ifoln,iBC,ic,BCtp_new,isub,err,ipl
      LOGICAl :: root_growth,repet,iterBC,BC_switch,n_apex,run
      !common /root_grow/ root_growth

      if (ave) then
         counter2=0
      endif

! initialisation
      sink=0._dp
      RootSk=0.0_dp

      DO ipl=1,nplant
! check current plant BC
         BCr_old=BCr_usr(ipl)
	  CALL SetBCroot(t,BCr_usr(ipl),BCtp_usr(ipl),level)
!calculation of Tpot if BCtp==2
         IF (BCtp_usr(ipl)==2) Tpot=BCr_usr(ipl)
!if stress at previous time step, then BCtp new=2 and currBctp=1
         IF ((stressBC).AND.((BCtp_usr(ipl)==1).OR.(BCr_old/=BCr_usr(ipl)))) THEN
!tocheck!!!with tuzet it will not work
!there was a stress with BC=2 but in the meantime user BC changed
            stressBC=.FALSE.
         ENDIF
         BCr_new=BCr_usr(ipl)
         BCtp_new=BCtp_usr(ipl)
         iterBC=.TRUE.
         BC_switch=.FALSE.
         DO WHILE (iterBC)
!adapt BC
            IF (((BCr_new/=curr_BCr(ipl)).AND.(.not.(stressBC))).OR.((BCtp_new/=curr_BCtp(ipl)).AND.(.not.(stressBC))).OR.(BC_switch)) THEN
!updating BC
	        curr_BCr(ipl)=BCr_new
	        old_BCtp=curr_BCtp(ipl)
	        curr_BCtp(ipl)=BCtp_new
!updating the BC part of the matrices Qbc and Cd
               iprvn=0!seed=first segment
               IF (curr_BCtp(ipl)==2) THEN  !Flux
                  IF (old_BCtp/=curr_BCtp(ipl)) THEN
                     DO iBC=1,nBCn  !all the brenches connected to the seed
                        irecn=nBC_irecn(iBC)!connected to the seed
                        IF (iBC==1) THEN
                           repet=.true.!delete previous value
                        ELSE
                           repet=.false.!add to previous value
                        ENDIF
                        err=SM_set(plantmatrix(ipl),iprvn+1,iprvn+1,Khr(irecn)/seglen(irecn),repet)!position (1,1)
                        err=SM_set(plantmatrix(ipl),iprvn+1,irecn+1,-Khr(irecn)/seglen(irecn),repet)!position (1,2)
                        err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,-Khr(irecn)/seglen(irecn),repet) !mathieu!position (2,1)
                        If(err/=0) Stop 'Could not insert element into plantmatrix'
                        Q_bc(irecn,ipl)=0._dp
                     ENDDO
                     Qi(iprvn,ipl)=0 !from seed has to be given, from irecn is already stated in DoussanMat
                  ENDIF
                  Q_bc(iprvn,ipl)=curr_BCr(ipl);!Q(0,0) a flow imposed
               ELSE IF (curr_BCtp(ipl)==1) THEN  !PH
!                  Cd(iprvn,iprvn)=1. !C(0,0)
                  IF (old_BCtp/=curr_BCtp(ipl)) THEN
                     err=SM_set(plantmatrix(ipl), iprvn+1,iprvn+1,1._dp,.true.)
                     If(err/=0) Stop 'Could not insert element into plantmatrix'
                  ENDIF
                  Qi(iprvn,ipl)=0._dp
                  Q_bc(iprvn,ipl)=(curr_BCr(ipl)+GH(iprvn)) !all the PH for root was given in total head!!
                  DO iBC=1,nBCn
                     irecn=nBC_irecn(iBC)!connected to the seed
                     IF (old_BCtp/=curr_BCtp(ipl)) THEN
                        err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,0._dp,.true.)
                        err=SM_set(plantmatrix(ipl),iprvn+1,irecn+1,0._dp,.true.)!is toch nul zowiezo..
!                        err=SM_set(plantmatrix(ipl),iprvn+1,iprvn+1,1._dp,.true.)
                        IF(err/=0) Stop 'Could not insert element into plantmatrix'
                     ENDIF
                     Q_bc(irecn,ipl)=(curr_BCr(ipl)+GH(iprvn))*Khr(irecn)/seglen(irecn)!all the PH for root was given in total head!!
                  END DO
               ENDIF
            ENDIF
            IF (root_growth) goto 201
            IF (oldT) goto 201
            IF (t.eq.t_begin+dt .and. iter .gt. 1 .or. t.gt.t_begin+dt) THEN
               CALL analytical_approach(hNew,ipl)
            ELSE
               old=.true.
            ENDIF
!for t=tstep. Initially no root flow is calculated, and not all the boundary conditions for the analytical problem are available, so solve
!first time step using the averaging way
!find current PH soil and create Q
201         IF ((old) .or. (root_growth)) THEN
               CALL average_method(hNew,ipl)
            ENDIF
            IF (nrec.eq.1) THEN!only seed
               IF (curr_BCtp(ipl)==2) THEN
                  Q_bc=Q_bc ! this must be checked MJ dec2008
               ENDIF
            ELSE
!put in format (1:nRec+1) for solving the system
               Qd2(1:nrec+1)=Qd(0:nrec,ipl)
               PHr_tot(1:nrec+1)=PHr(1:nrec+1,ipl)+GH(0:nrec)!total water potential
               IF (it1) THEN
                  PHr_tot(1:nrec+1)=PHs(0:nrec,ipl)   ! initial xylem PH -> could also be zero but then convergence should take longer (this should be the common case)
               ENDIF
               IF (curr_BCtp(ipl)==1) THEN
                  PHr_tot(1)=curr_BCr(ipl)+GH(0) !bound value should always be incorporated if a PH is imposed
               ENDIF
               itol=3
               itmaxl=(nrec+1)*3 !maximum iterations for biconjugate gradient-> needed in case of root growth
! solve with biconjugated gradient method (Numerical recipies)
               CALL linbcg(nrec+1,Qd2,PHr_tot,itol,maxerr,itmaxl,iter_loc,err_loc,ipl)
               IF (err_loc.LE.maxerr) THEN
                  nrecOld=nrec
               ELSE
                  print *,'No convergence in the root system '
                  stop
               ENDIF
            ENDIF

!====================== end solve root system ====================================
!----------------------calculate axial root flow-------------------------------
!for each branch cumulative flow is calculated
            DO ibr=1,nbr
               n_apex=.false.
!* find the tip segment of the branch 'ibr'
               irecn=nrec
               DO WHILE (ibrseg(irecn).NE.ibr)
	           irecn=irecn-1
	        END DO
!* the first one we find is an apex
               IF (seglen(irecn)<1.E-20) THEN ! skip this segment too small to be taken
	           ifoln=irecn ! following node ID
	           irecn=irecpr(irecn) !current node ID
               ELSE
	           n_apex=.true.!ifoln does not exist
               ENDIF
	        IF (irecn==0) THEN !there exists a branch ibr but not yet any segment!
	           run=.false.
	        ELSE
                  run=.true.
	        ENDIF
!* then the rest of the branch up to the seed or the embranchment
               DO WHILE (run)
!* "upper" node
	           iprvn=irecpr(irecn)
!	            IF (n_apex) THEN
		    axialRootFlow(irecn,ipl) = Khr(irecn)*(PHr_tot(irecn+1)-PHr_tot(iprvn+1))/seglen(irecn)
!                  ELSE
!		        axialRootFlow(irecn,ipl) = Khr(irecn)*(PHr_tot(irecn+1)-PHr_tot(iprvn+1))/seglen(irecn) !No difference between the two condition results???(Couvreur dec 2009)
!                  ENDIF
! if reached the seed or the embranchement => change branch
                  IF (iprvn==0.OR.(ibrseg(iprvn).NE.ibrseg(irecn))) run=.false.				!Code simplification (Couvreur dec 2009)
!                     axialRootFlow(irecn,ipl) = Khr(irecn)*(PHr_tot(irecn+1)-PHr_tot(iprvn+1))/seglen(irecn)
!                     run=.false.
!                  ELSEIF (ibrseg(iprvn).NE.ibrseg(irecn)) THEN!start of the branch but not from the seed -> gaat van onder na boven, duz iprvn (nieuwe positie) is nu niet van die ene branch maar van een side branch
! 	               axialRootFlow(irecn,ipl) = Khr(irecn)*(PHr_tot(irecn+1)-PHr_tot(iprvn+1))/seglen(irecn) !No difference between the two condition results???(Couvreur dec 2009)
!                     run=.false.
!                  ENDIF
                  if (irecn .eq. 0) then
                     if (curr_BCtp(ipl) == 2) then !flux
                        axialRootFlow(irecn,ipl) = Q_bc(irecn,ipl)
                     elseif (curr_BCtp(ipl) == 1) then
                        axialRootFlow(irecn,ipl) = Khr(irecn)*(PHr_tot(irecn+1)-PHr_tot(iprvn+1))/seglen(irecn+1)
                     endif
                     run=.false.
                  endif
!get total flow; which is addition of all flows at each large branch connecting to the seed		Quid des branches latérales ? (Couvreur dec 2009)
                  if (iprvn .eq. 0) then
                     totalAxialFlow = totalAxialFlow + axialRootFlow(irecn,ipl)
                  endif
! definition for the next run of the loop
                  ifoln=irecn
                  irecn=iprvn
!* from here, not an apex
                  n_apex=.false.
               END DO ! branch nodes loop
            END DO ! branches loop
!write(*,*)'total axialRootflow',totalAxialFlow
            totalAxialFlow = 0
!---------------------------------end axial root flow calculation-----------------
!***********************************************************************
! calculate SINK term
!**********************************************************************

! initialisation
            betaw=0.
            Jintot=0.
            Jcol=0.0_dp
            Voltot=0.0_dp
            DO irecn=0,nrec
               Joutr(irecn) = Qi(irecn,ipl)*(PHs(irecn,ipl)-PHr_tot(irecn+1)) ![L3/T] (soil is still from 0 to nrec while root is from 1 to nrec+1)
               sinkR(irecn,ipl)=Joutr(irecn) !root flow!!!!!!!!!!!! [L3/T]
               Jintot=Jintot+Joutr(irecn)/(dxgrid*dygrid*dzgrid*nsub(irecn,ipl))		!Why nsub? Why calculate Jintot? (Couvreur jan 2010)
               IF (t.eq.t_begin+dt .and. iter .gt. 1 .or. t.gt.t_begin+dt) THEN
                  checksub=0
                  DO isub=1,nsub(irecn,ipl)
                     dVolSk=0
                     DO ic=1,8
                        corner_ic=loc_Q(irecn,ic,isub,ipl)
                        !delta flux extracted
                        dVolSk=Joutr(irecn)*w_dis(irecn,ic,isub,ipl)/sum_dis(irecn,isub,ipl)*w_sub(irecn,isub,ipl)
                        !sink term
                        sink(corner_ic)=sink(corner_ic)+dVolSk/(dxgrid*dygrid*dzgrid)/Wn(corner_ic)
                        betaw(corner_ic)=1._dp
                        !sum of extracted fluxes
                        RootSk=RootSk+dVolSk
                        Voltot=Voltot+w_dis(irecn,ic,isub,ipl)/sum_dis(irecn,isub,ipl)*w_sub(irecn,isub,ipl)
                     END DO
                     checksub=checksub+w_sub(irecn,isub,ipl)
                  END DO !end sub-segments
!first time step
               ELSE 
                  DO isub=1,nsub(irecn,ipl)
                     DO ic=1,8
                        corner_ic=loc_Q(irecn,ic,isub,ipl)
                        sink(corner_ic)=sink(corner_ic)+Joutr(irecn)*w_sub(irecn,isub,ipl)*w_dis(irecn,ic,isub,ipl) &
/sum_dis(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/Wn(corner_ic)
                        betaw(corner_ic)=1._dp ! betaw is defined for the Reset subroutine
                     END DO
                  END DO
               ENDIF
               Jcol=Jcol+Joutr(irecn)
            END DO ! nodes loop
!actual transpiration
            Tact(ipl)=Jcol

!gravitational component subtracted again to represent data in outRoo.x without gravity
            PHr(1:nrec+1,ipl)=PHr_tot(1:nrec+1)-GH(0:nrec)
            PHs(0:nrec,ipl)=PHs(0:nrec,ipl)-GH(0:nrec)

            IF (ave) THEN
               PH_root_sub(0:nrec,1:20)=PH_root_sub(0:nrec,1:20)-cp_mean(0:nrec,3,1:20,ipl)
            ELSEIF (.not.(ave)) THEN
               PH_root_sub(0:nrec,1:20)=PH_root_sub(0:nrec,1:20)-cent(0:nrec,3,1:20,ipl)
            ENDIF
            IF (((curr_BCtp(ipl)==2).OR.(stressBC)).AND.(.not.(BC_switch))) THEN
!check stress a the collar in case of BC flux
!check flux at the collar in case of PH BC+stress
!if swicth at previous iterBC, no check
!print *,'checkstress', BC_switch
               CALL Rootstress(PHr(1,ipl),BCr_new,BCtp_new,iterBC,BC_switch,Jcol,ipl)! change BCtp_new=currBCtp to 1 if stress occurs
	     ELSE
	        iterBC=.FALSE.
	        BC_switch=.FALSE.
 	     ENDIF
         ENDDO ! do while iterBC
      ENDDO ! plant loop
      RETURN
      END SUBROUTINE SolveRoot
!***************************************************************************
SUBROUTINE analytical_approach(hNew,ipl)
!      USE OMP_LIB
      USE TypeDef
      USE Paramdata,only: MaxNod, pi
      USE GridData, only: xgrid,ygrid,zgrid,dxgrid,dygrid,dzgrid,elmnod,nPt,Axy,Bxy,Dxy,Exy
      USE DoussanMat, only :old,oldT,ave,eqDis,nsub,PHs,loc_q,cp_mean,voxel_no,no_voxels, &
 voxel_node,numnodes_voxel,h_mat2,hx_min,phr_sub,n,indexvalue,counter2,phi_mat, &
PH_micro2,ph_root_sub,Qd,w_dis,w_sub,cent,Intc,PHr,Lr,sum_dis,switchcriterion,Qi,Q_bc, &
tcheck,tcheck2,count_nodes,k_ave,B
      USE RootData, only: nrec,seglen,segsur,nbr,ibrseg,irecpr
      use SolData
      IMPLICIT NONE
      integer(sp) :: checker1=0,blaat,irecnTemp,isubTemp,ki
      REAL,intent(in) :: hNew(maxnod)
      REAL(sp) :: beta
      REAL(sp) :: y(1:n),h_mat(1:n),arr(1:n-1),Phi(1:n-1)
      REAL(sp),dimension(1:8) ::h_matK,temp2
      real(sp) :: K_k,C_k,theta_k,temp
      REAL(sp) :: cumsum(1:n-1),Phi_rend
      REAL(sp) :: z1,z2,frac1,frac2,xCent,yCent,Zcent,xInt,yInt,zInt
      REAl(sp) :: testx(2),testy(2),testz(2),xgrid_max,xgrid_min,ygrid_max,ygrid_min,mXY(2),point(1:4,1:2)
      REAL(sp) :: zgrid_max,zgrid_min,mZX(2),mZY(2)
      REAl(sp) :: PH_micro(1:n-1)
      REAL(sp) :: theta(1:n),K(1:n),C(1:n),rho,r0,r(1:n-1),rend
      REAL (dp) ::sumW,psub,q_root,q_out
      real (sp) :: pos_cent,testvar
      REAL(dp),dimension(1:8) ::pond_i,q_tot_c
      real(sp),dimension(1:4,1:4) :: dis,Flowpoint_rel,PHpoint_rel
      real(sp),dimension(1:4) :: mFlows,mPHs,PH_inter,Flow_inter
      real(sp) :: meanPHS,meanFlows
      INTEGER(sp) ::j,counter,l,w,Flow_corner(1:8,1:6)
      INTEGER(sp) ::corner(1:8),var1,kl,imin,ipl
      INTEGER(sp) :: irecn,i,isub,pos,ja=0,irecn_tmp,isub_tmp
      REAL(sp) :: z_angle,x_angle,y_angle,x1,x2,y1,y2,q_root_tmp(1:indexValue),pond_subTemp
      LOGICAl :: zalign,yalign,xalign,extra_checker=.false.,countja
!parameters for counting nodes that are under limiting conditions (ks < lr*r0)
      countja=.false.
      count_nodes=0
      DO irecn=0,nrec !loop over segments
         old=oldT
         x_angle=0
         y_angle=0
         z_angle=0
 	  !initialize soil PH matrix
         PHs(irecn,ipl)=0._dp
         DO isub=1,nsub(irecn,ipl) !loop over subsegments
            dis=0
            Flowpoint_rel=0
            PHpoint_rel=0
            mFlows=0
            mPHs=0
            PH_inter=0
            Flow_inter=0
!---------------------------calculate the flux in each soil node with Darcy's law--------
            if (.not.eqDis) then
               DO i=1,8 !for each corner of a cube
                  corner(i)=loc_Q(irecn,i,isub,ipl)
               END DO
               !determine nodes surrounding a soil node to determine the flow
               call FlowNeighb(corner,Flow_corner)
               call flux_node(hNew,Flow_corner,corner,irecn,isub,q_tot_c)
            elseif (eqDis) then !method C
               q_tot_c(1:8)=0
            endif
!---------------------------end flow calculation-----------------------------------
            !initialization and r0 definement
            checker1=0 !reset for each cp
            DO i=1,8
               pond_i(i)=w_dis(irecn,i,isub,ipl)
               corner(i)=loc_Q(irecn,i,isub,ipl)
            END DO
            zalign=.false.
            yalign=.false.
            xalign=.false.
!not always true: adjust accordingly
            if (irecn .eq. 0) then
               r0=5e-2 !set to minimal radius for root collar ID only-> not important
            else
               if (segsur(irecn) .eq. 0 .or. seglen(irecn) .eq. 0) then
                  r0 = 5e-2 !minimal root radius (in cm)
               else
                  r0 = 1/(2*pi)*segsur(irecn)/seglen(irecn)
               end if
            endif
!-----------------------------------------------------------------------
!
! REDUCTION from 3D -> 2D plane, currently averaging procedure is taken
! in axial direction (i.e. direction perpendicular to radial flow)  which is sufficient (checked)
!----------------------------------------------------------------------
            !for the one_root procedure; the center point of a segment is the middle of a soil voxel
            IF (ave) THEN
               xCent = cp_mean(irecn,1,isub,ipl)!cent(irecn,1,isub,ipl)
               yCent = cp_mean(irecn,2,isub,ipl)!cent(irecn,2,isub,ipl)
               zCent = cp_mean(irecn,3,isub,ipl)!cent(irecn,3,isub,ipl)
               imin = voxel_no(irecn,isub) !voxel numer
               if (counter2(imin) .ne. 0) then !check of conductivity drop in this voxel is already calculated
                  irecnTemp = voxel_node(imin,1,counter2(imin))
                  isubTemp = voxel_node(imin,2,counter2(imin))
                  pond_subTemp = w_sub(irecnTemp,isubTemp,ipl)
                  extra_checker=.true.
                  goto 210 !to calculate analytical solution only once for the large root node
               elseif (counter2(imin) .eq. 0 .and. no_voxels(imin) .ne. 0) then
                  do blaat=1,no_voxels(imin)
                     if (voxel_node(imin,1,blaat) .eq. irecn .and. voxel_node(imin,2,blaat) &
 .eq. isub .and. w_sub(irecn,isub,ipl).ne.0) counter2(imin)=blaat
                  enddo
               endif
            ELSEIF (.not.(ave)) THEN ! if not one root method each root node is a center point
               xCent = cent(irecn,1,isub,ipl)
               yCent = cent(irecn,2,isub,ipl)
               zCent = cent(irecn,3,isub,ipl)
            ENDIF
            !interception point x,y,zInt with soil cube
            xInt = Intc(irecn,1,isub,ipl)
            yInt = Intc(irecn,2,isub,ipl)
            zInt = Intc(irecn,3,isub,ipl)
!---------------------------derive x,y,z-allignment and segment angle------------------
!Get aligment of a root in the current voxel. If an intersection point of a root with the cube is identical to the center point of the root, pick a direction dependent on the surface it crosses
!otherwise calculate the segment angle. first from a z-directional point of view, if the branching angle is smaller then a given range then look if it is
!alligned in x-direction (using a projection of the z-coordinate on the 2D generated plain)
!if this angle doesnt comply to a criterion then the allignement is performed in y-direction
! weigh factor Wb is calculated but not needed in remainder calculations
!dimensions of the voxel, maxima and minima in each direction
            zgrid_max = maxval(zgrid(corner(1:8)))
            zgrid_min = minval(zgrid(corner(1:8)))
            xgrid_max = maxval(xgrid(corner(1:8)))
            xgrid_min = minval(xgrid(corner(1:8)))
            ygrid_max = maxval(ygrid(corner(1:8)))
            ygrid_min = minval(ygrid(corner(1:8)))
!get allignment direction
            if (xCent-xInt .eq. 0 .and. yCent-yInt .eq. 0 .and. zCent-zInt .eq. 0) then
               if (zCent .eq. zgrid_max .or. zCent .eq. zgrid_min) then
                  zalign=.true.
               elseif (yCent .eq. ygrid_max .or. yCent .eq. ygrid_min) then
                  yalign=.true.
               elseif (xCent .eq. xgrid_max .or. xCent .eq. xgrid_min) then
                  xalign=.true.
               endif
            else
!goneometrie; beta is radial angle betwen z difference and the length of the segment-> convert to degree
               beta = asin(abs(zInt-zCent)/sqrt((xInt-xCent)*(xInt-xCent)+ &
(yInt-yCent)*(yInt-yCent)+(zInt-zCent)*(zInt-zCent)))
               z_angle=abs(beta/(2*pi))*360
               if (z_angle .le. 90.1 .and. z_angle .ge. 45) then !certain error deviatin given to criterion
                  zalign=.true.
               else
                  zInt=zCent !projection in 2d plane to get angle in x-or y direction
                  beta = acos(abs(xInt-xCent)/sqrt((xInt-xCent)*(xInt-xCent)+ &
(yInt-yCent)*(yInt-yCent)+(zInt-zCent)*(zInt-zCent)))
                  x_angle=abs(beta/(2*pi))*360
                  x_angle=90-x_angle !becoz of cosine, angle interpreted from 0 to 45 and angle has to be adjusted for weight factor
                  if (x_angle .le. 90 .and. x_angle .ge. 45) then
                     xalign=.true.
                  else
                     yalign=.true.
                  endif
               endif
            endif
!-----------------------------------end allignment-----------------------------------
!for each allocated root direction in a voxel a checker1 is set which determines where the cp is located, in the voxel or at the edges. If it is located in the  ! voxel then the analytical solution is applied, otherwise the old method is sufficient
            !set checker1
            if (zalign) then
               if (abs(xCent-xgrid_max) .gt. r0 .and. abs(xCent-xgrid_min).gt.r0 &
 .and. abs(yCent-ygrid_max) .gt. r0 .and. abs(yCent-ygrid_min).gt.r0 ) then !in plane
                  checker1 = 1
               else !if it is at a boundary -> apply a simple method to obtain PH at the soil-root interface (for the moment)
                  checker1=0
               end if
            elseif (xalign) then
               if (abs(yCent-ygrid_max) .gt. r0 .and. abs(yCent-ygrid_min).gt. r0 &
.and. abs(zCent-zgrid_max) .gt. r0 .and. abs(zCent-zgrid_min).gt.r0 ) then !in plane
                  checker1 = 1
               else !if it is at a boundary -> apply a simple method to obtain PH at the soil-root interface (for the moment)
                  checker1=0
               end if
            elseif (yalign) then
               if (abs(xCent-xgrid_max) .gt. r0 .and. abs(xCent-xgrid_min).gt.r0 &
.and. abs(zCent-zgrid_max) .gt. r0 .and. abs(zCent-zgrid_min).gt.r0 ) then !in plane
                  checker1 = 1
               else !if it is at a boundary -> apply a simple method to obtain PH at the soil-root interface (for the moment)
                  checker1=0
               end if
            endif
!Average for the specific direction chosen from a 3D to a 2D mesh (for application of the 2D analytical solution). In direction plane is orientated
!To be clearer; using the allignment of the root node the PH and Flow boundary conditions
!are transferred to the 4 corner nodes spanning the 2D domain in which the
!analytical approach is applied
!
!               C1----------------C2
!		|		  |
!		|		  |
!		|		  |
!		|	o	  |
!		|		  |
!		|		  |
!		|		  |
!		C3----------------C4
!
! z-align, x-align, y-align are discussed; basicly the same routines;
            IF (checker1.eq.1 .and. (zalign)) THEN
               do i=1,4
                  z1 = zgrid(corner(i))
                  z2 = zgrid(corner(i+4))
                  frac1 = 1/abs(zCent-z1)
                  frac2 = 1/abs(zCent-z2)
                  if (abs(zCent-z1) .eq. 0) then
                     PH_inter(i) = hnew(corner(i))
 !ph bound. cond.
                     Flow_inter(i) = q_tot_c(i)
 !flow bound. cond.
                  elseif (abs(zCent-z2) .eq. 0) then
                     PH_inter(i) = hnew(corner(i+4))
                     Flow_inter(i) = q_tot_c(i+4)
                  else
                     PH_inter(i) = (frac1*hnew(corner(i)) + frac2*hnew(corner(i+4)))/(frac1+frac2)
                     Flow_inter(i) = (frac1*q_tot_c(i) + frac2*q_tot_c(i+4))/(frac1+frac2)
                     !corner -> inter (2D plane)
                  endif
               enddo
               testx(1) = abs(xgrid_max-xCent)
               testx(2) = abs(xgrid_min-xCent)
               mXY(1) = minval(testx) !minimum distance in x-dimension
               testy(1) = abs(ygrid_max-yCent)
               testy(2) = abs(ygrid_min-yCent)
               mXY(2) = minval(testy) !minimum distance in y-direction
! get minimum radius from root node to a 2D edge; to be used as outer radius
!for the analytical approach
               rend = minval(mXY)     !outer radius of the soil cylinder
               if (eqDis) then !get an equal radius dependent on the number of root nodes in this voxel
                  testvar=numNodes_voxel(irecn,isub)**(1d0/3d0)
                  rend = 1/2d0*sqrt((dzgrid)/floor(testvar))
               endif
               if (rend .lt. r0) then
                  ja=1
                  goto 203
               endif
!the dashed arrow line is set as the minimal distance in this example; i.e. rend (outer radius)
!             C1------0---------C2
!		|	|	  |
!		|	|	  |
!		|	|	  |
!		0<----->o---------0
!		|	|	  |
!		|	|	  |
!		|	|	  |
!		C3------0---------C4
!new coordinate system 2D based on the radius of the microscopic model (question is, how large is ratio rend/r0)
!so the outer radius has intersection points with the 2D plane as drawn above (the zeros's (0) are the
!intersection points with the 2D plane
!these intersection points are calculated here; for method B,C (see MS)
               if (eqDis) then !middle of voxel as point of departure for bnd.cond.
                  point(1,1) = (xgrid_max+xgrid_min)/2+rend
                  point(1,2) = (ygrid_max+ygrid_min)/2
                  point(2,1) = (xgrid_max+xgrid_min)/2
                  point(2,2) = (ygrid_max+ygrid_min)/2+rend
                  point(3,1) = (xgrid_max+xgrid_min)/2-rend
                  point(3,2) = (ygrid_max+ygrid_min)/2
                  point(4,1) = (xgrid_max+xgrid_min)/2
                  point(4,2) = (ygrid_max+ygrid_min)/2-rend
               else
                  point(1,1) = xCent+rend
                  point(1,2) = yCent
                  point(2,1) = xCent
                  point(2,2) = yCent+rend
                  point(3,1) = xCent-rend
                  point(3,2) = yCent
                  point(4,1) = xCent
                  point(4,2) = yCent-rend
               endif
!calculate the pressure head, based on a distance average, for each new coordinate surrounding the CP(center ponint, i.e. the root node)---- 2D!!!!!!!! z-alignment-> x and y-direction is evaluated
!so outer boundary conditions are mapped on the 4 points that represent the outer circle
!The new PH in each intersection point (zero) with 2D domain is mPHs, here only part of that solution is given
!due to distance saciling!! so "_rel" is intermediate value for PH or Flow (later these values are put to the correct valuesin mPHs,mFlow
               do j=1,size(point,1)
                  do i=1,4
                     dis(j,i) = SQRT((point(j,1)-xgrid(corner(i)))**2+(point(j,2)-ygrid(corner(i)))**2)
                     PHpoint_rel(j,i) = 1/dis(j,i)*PH_inter(i) !inter values being the values in the 4 corner nodes of the 2D domain
                     Flowpoint_rel(j,i) = 1/dis(j,i)*Flow_inter(i)
                     !inter -> point (distance)
                  end do
               end do
            ELSEIF (checker1.eq.1 .and. (xalign)) THEN
               kl=1
               do i=1,7,2
                  x1 = xgrid(corner(i))
                  x2 = xgrid(corner(i+1))
                  frac1 = 1/abs(xCent-x1)
                  frac2 = 1/abs(xCent-x2)
                  if (abs(xCent-x1) .eq. 0) then
                     PH_inter(kl) = hnew(corner(i))
                     Flow_inter(kl) = q_tot_c(i)
                  else if (abs(xCent-x2) .eq. 0) then
                     PH_inter(kl) = hnew(corner(i+1))
                     Flow_inter(kl) = q_tot_c(i+1)
                  else
                     PH_inter(kl) = (frac1*hnew(corner(i)) + frac2*hnew(corner(i+1)))/(frac1+frac2)
                     Flow_inter(kl) = (frac1*q_tot_c(i) + frac2*q_tot_c(i+1))/(frac1+frac2)
                     !corner -> inter (2D plane)
                  endif
                  kl=kl+1
               end do
               testz(1) = abs(zgrid_max-zCent)
               testz(2) = abs(zgrid_min-zCent)
               mZY(1) = minval(testz) !minimum distance in z-dimension
               testy(1) = abs(ygrid_max-yCent)
               testy(2) = abs(ygrid_min-yCent)
               mZY(2) = minval(testy) !minimum distance in y-direction
               rend = minval(mZY)     !outer radius of the soil cylinder
               if (eqDis) then
                  testvar=numNodes_voxel(irecn,isub)**(1d0/3d0)
                  rend = 1/2d0*sqrt((dxgrid)/floor(testvar))
               endif
               if (rend .lt. r0) then
                  ja=1
                  goto 203
               endif
               !new coordinate system 2D based on the radius of the microscopic model (question is, how large is ratio rend/r0)
               if (eqDis) then !middle of voxel as point of departure for bnd.cond.
                  point(1,1) = (zgrid_max+zgrid_min)/2+rend
                  point(1,2) = (ygrid_max+ygrid_min)/2
                  point(2,1) = (zgrid_max+zgrid_min)/2
                  point(2,2) = (ygrid_max+ygrid_min)/2+rend
                  point(3,1) = (zgrid_max+zgrid_min)/2-rend
                  point(3,2) = (ygrid_max+ygrid_min)/2
                  point(4,1) = (zgrid_max+zgrid_min)/2
                  point(4,2) = (ygrid_max+ygrid_min)/2-rend
               else
                  point(1,1) = zCent+rend
                  point(1,2) = yCent
                  point(2,1) = zCent
                  point(2,2) = yCent+rend
                  point(3,1) = zCent-rend
                  point(3,2) = yCent
                  point(4,1) = zCent
                  point(4,2) = yCent-rend
               endif
               do j=1,size(point,1)
                  do i=1,4
                     if (i.eq.1) kl=1 !y and z-coord of these nodes needed
                     if (i.eq.2) kl=3
                     if (i.eq.3) kl=5
                     if (i.eq.4) kl=7
                     dis(j,i) = SQRT((point(j,1)-zgrid(corner(kl)))**2+(point(j,2)-ygrid(corner(kl)))**2)
                     PHpoint_rel(j,i) = 1/dis(j,i)*PH_inter(i)
                     Flowpoint_rel(j,i) = 1/dis(j,i)*Flow_inter(i)
                     !inter -> point (distance)
                  end do
               end do
            ELSEIF (checker1.eq.1 .and. (yalign)) THEN
               kl=1
               do i=1,4
                  var1=i
                  IF (i.gt.2) var1=var1+2
                  y1 = xgrid(corner(var1))
                  y2 = xgrid(corner(var1+2))
                  frac1 = 1/abs(yCent-y1)
                  frac2 = 1/abs(yCent-y2)
                  if (abs(yCent-y1) .eq. 0) then
                     PH_inter(kl) = hnew(corner(var1))
                     Flow_inter(kl) = q_tot_c(var1)
                  elseif (abs(yCent-y2) .eq. 0) then
                     PH_inter(kl) = hnew(corner(var1+2))
                     Flow_inter(kl) = q_tot_c(var1+2)
                  else
                     PH_inter(kl) = (frac1*hnew(corner(var1)) + frac2*hnew(corner(var1+2)))/(frac1+frac2)
                     Flow_inter(kl) = (frac1*q_tot_c(var1) + &
frac2*q_tot_c(var1+2))/(frac1+frac2)
                     !corner -> inter (2D plane)
                  endif
                  kl=kl+1
               enddo
               testz(1) = abs(zgrid_max-zCent)
               testz(2) = abs(zgrid_min-zCent)
               mZX(1) = minval(testz) !minimum distance in z-dimension
               testx(1) = abs(xgrid_max-xCent)
               testx(2) = abs(xgrid_min-xCent)
               mZX(2) = minval(testx) !minimum distance in y-direction
               rend = minval(mZX)
               if (eqDis) then
                  testvar=numNodes_voxel(irecn,isub)**(1d0/3d0)
                  rend = 1/2d0*sqrt((dygrid)/floor(testvar))
               endif
               if (rend .lt. r0) then
                  ja=1
                  goto 203
               endif
               !new coordinate system 2D based on the radius of the microscopic model (question is, how large is ratio rend/r0)
               if (eqDis) then !middle of voxel as point of departure for bnd.cond.
                  point(1,1) = (zgrid_max+zgrid_min)/2+rend
                  point(1,2) = (xgrid_max+xgrid_min)/2
                  point(2,1) = (zgrid_max+zgrid_min)/2
                  point(2,2) = (xgrid_max+xgrid_min)/2+rend
                  point(3,1) = (zgrid_max+zgrid_min)/2-rend
                  point(3,2) = (xgrid_max+xgrid_min)/2
                  point(4,1) = (zgrid_max+zgrid_min)/2
                  point(4,2) = (xgrid_max+xgrid_min)/2-rend
               else
                  point(1,1) = zCent+rend
                  point(1,2) = xCent
                  point(2,1) = zCent
                  point(2,2) = xCent+rend
                  point(3,1) = zCent-rend
                  point(3,2) = xCent
                  point(4,1) = zCent
                  point(4,2) = xCent-rend
               endif
               do j=1,size(point,1)
                  do i=1,4
                     if (i.eq.1) kl=1
                     if (i.eq.2) kl=2
                     if (i.eq.3) kl=5
                     if (i.eq.4) kl=6
                     dis(j,i) = SQRT((point(j,1)-zgrid(corner(kl)))**2+(point(j,2)-xgrid(corner(kl)))**2)
                     PHpoint_rel(j,i) = 1/dis(j,i)*PH_inter(i)
                     Flowpoint_rel(j,i) = 1/dis(j,i)*Flow_inter(i)
                     !inter -> point (distance)
                  end do
               end do
            ENDIF
!IF CP is located in the plane (most of the cases)
            210 continue   !one root approach continues here if applicable
            if (checker1 .eq. 1 .or. (extra_checker)) then
               if (extra_checker) then
                  extra_checker=.false.
                  goto 200
               endif
!-----------------------------------------------------------------------------------------------------------------
! Determine the soil cylinder
!-----------------------------------------------------------------------------------------------------------------
!after outer boundary conditions at the new circle surrounding the cp are known for each allignment, calculate the average PH and flow for these boundary conditions
!these average PH and Flow values (meanPHs and meanFlows) at the outer radius are the boundary conditions
!for the analytical approach
               do j=1,4
                  mPHs(j) = sum(PHpoint_rel(j,:))/sum(1/dis(j,:))
                  mFlows(j) = sum(Flowpoint_rel(j,:))/sum(1/dis(j,:))
               end do
               meanPHs = (sum(mPHs(1:4))/4)
               meanFlows = (sum(mFlows(1:4))/4)
               !fractioned value should be used (subsegmentation)
               meanPHs = meanPHs*w_sub(irecn,isub,ipl)
               meanFlows = meanFlows*w_sub(irecn,isub,ipl)
!-----------------------------------------------------------------------------------------------------------------
!partitioning the radii used for the anlytical approach (in log scale); closer
!radii near the root, larger radii at outer radius
!-----------------------------------------------------------------------------------------------------------------
               !spatial grid
               call logspace(r0,rend,n-1,y)
               r=y(1:n-1)
               rho = rend/r0
!-----------------------------------------------------------------------
!
! Calculate BOUNDARY conditions for analytical solution: Phi_rend(based on meanPHs),q_end,q_root
!
!-----------------------------------------------------------------------------------------------------------------
! Calculate Phi_rend -> obtained from the PH boundary condition (the PH at outer radius of microscopic model)
!-----------------------------------------------------------------------------------------------------------------
               call logspace(hx_min,meanPHs,n,y)
               h_mat = -y
               call SetMat_ana(h_mat,theta,K,C)
! calculate matric flux potential integral K(h)dh -> numerical integration
               !cumsum = cumulitive sum
               arr=abs(h_mat(2:size(h_mat))-h_mat(1:size(h_mat)-1))*(K(2:size(K))+K(1:size(K)-1))/2
               do i=1,size(arr)
                  cumsum(i) = sum(arr(1:i))
               end do
               !last value from integral is at r=rend
               Phi_rend=cumsum(size(cumsum))
!--------------q_root------------------------
!get current voxel number; for a loop over all the nodes in current voxel, get their irecn and isub values; calculate q_root for all root nodes
               if (ave) then
                  imin = voxel_no(irecn,isub)
                  do ki=1,no_voxels(imin)
                     irecn_tmp = voxel_node(imin,1,ki)
                     isub_tmp = voxel_node(imin,2,ki)
                     PHr_sub(irecn_tmp,isub_tmp) = w_sub(irecn_tmp,isub_tmp,ipl)*PHr(irecn_tmp+1,ipl)
                     q_root_tmp(ki)= Lr(irecn_tmp) * (PH_root_sub(irecn_tmp,isub_tmp) &
- PHr_sub(irecn_tmp,isub_tmp))
                  enddo
!===============================================================================
! chose the way how q_root is defined:
! 1) sum of all root flows, consider as one sucking root
! PHr_sub is not used, then only to check wheter or not limiting conditions have been reached -> see PHi
                  PHr_sub(irecn,isub) = w_sub(irecn,isub,ipl)*PHr(irecn+1,ipl)
!1)
                  q_root = sum(q_root_tmp(1:ki)) !total uptake -> one big root
!!-------------------------------------------------------------------------------
               elseif (.not.(ave)) then
!PHxylem = PHr; PHr_sub = subsegment xylemPH
                  PHr_sub(irecn,isub) = w_sub(irecn,isub,ipl)*PHr(irecn+1,ipl)
!CALCULATE BOUNDARY CONDITION q_root
                  !should be fractioned
                  q_root= Lr(irecn) * (PH_root_sub(irecn,isub) - PHr_sub(irecn,isub))
               endif
!CALCULATE BOUNDARY CONDITION q_out
               !steady state -> condition for doussan matrix
               if (meanFlows .le. q_root*(r0/rend)) then
                  q_out = meanFlows
               else
                  q_out = q_root*(r0/rend)
               endif
               !de willigen, q_out=zero
               !if (eqDis) q_out=0
!**************************************************************************
!------------------------------ calculate PHI -----------------------------
!**************************************************************************
!-----------------------------------------------------
! actual PHi calculation; hx_min equals lower boundary integral from flux density -> Phi_r_root = 0
!-----------------------------------------------------
               if (PHr_sub(irecn,isub) .le. hx_min*w_sub(irecn,isub,ipl)) then !if PHsoil-root interface = PHxylem limiting -> flow is zero (ultimate case limiting root PH
                  Phi = ((Phi_rend + q_out*rend/2*log(1/rho) ) / &
( rho**2-1 + 2*rho**2*log(1/rho) ) ) * (r**2/r0**2-1 + 2*rho**2*log(r0/r)) &
+ rend*q_out*log(r/r0)
               else !always true, also when hx_min is reached for PHxylem. If finally PHsoil-root inter = PHxylem (=hx_min), see above
                  Phi = Phi_rend + (q_root*r0 - q_out*rend)* &
(r**2/r0**2/(2*(1-rho**2)) + rho**2/(1-rho**2)*(log(rend/r)-1/2.)) &
 + q_out*rend*log(r/rend)
               end if
               counter=1
               !determine PH from Phi in microscopic domain
               do l=1,size(Phi)
                  do w=counter,size(Phi_mat)
                     if (Phi_mat(w) .lt. Phi(l)) then
                        counter = w;
                     elseif (Phi_mat(w) .gt. Phi(l)) then
                        EXIT
                     end if
                  end do
                  PH_micro(l)=h_mat2(counter)
               end do
               PH_micro2(irecn,1:size(Phi),isub)=PH_micro
!water potential criterion; check B*K_ave (soil) with Lr*r0 (see MS)
               if (PH_micro(n-1)-PH_micro(1) .eq. 0) THEN
                  k_ave(irecn,isub) = 0
                  ja=1
                  goto 203
               else
                  k_ave(irecn,isub) = (Phi(n-1)-Phi(1))/(PH_micro(n-1)-PH_micro(1))
               endif
               B(irecn,isub)=2*(1-rho**2)/(-2*rho**2*(log(rho)-1/2.)-1)
               IF (k_ave(irecn,isub)*B(irecn,isub) .lt. Lr(irecn)*r0) then
                  countja=.true.
               endif
! set switchcriterion; if B*Ks > Lr*r0
! set switchcriterion; if B*Ks < Lr*r0
! set switchcriterion; if B*Ks +_ 5% < > Lr*r0
               IF ((ave) .and. switchcriterion.eq.1) PH_micro(1)=meanPHs+ &
r0/k_ave(irecn,isub)*q_out*rho*log(1/rho)+r0/ &
(k_ave(irecn,isub)*B(irecn,isub))*q_out*rho
               IF (k_ave(irecn,isub)*B(irecn,isub) .lt. Lr(irecn)*r0+Lr(irecn)*r0*5/100 .and. k_ave(irecn,isub)*B(irecn,isub) .gt. &
Lr(irecn)*r0-Lr(irecn)*r0*5/100) THEN
                  switchCriterion=2
                  PH_micro(1)=meanPHs/2+PHr(irecn,ipl)/2+1/(Lr(irecn)*2)*q_out*rho*log(1/rho)+ &
1/(Lr(irecn)*2)*q_out*rho
               ELSEIF (k_ave(irecn,isub)*B(irecn,isub) .lt. Lr(irecn)*r0-Lr(irecn)*r0*5/100) THEN
                  switchCriterion=3
                  PH_micro(1)=PHr(irecn,ipl)+B(irecn,isub)/Lr(irecn)*q_out*rho*log(1/rho)+1/Lr(irecn)*q_out*rho
               ENDIF
               IF (switchCriterion .ne. 1 .and. (.not.(tcheck)) ) then
                  tcheck2=.true.
                  tcheck=.true.
               ENDIF
               if (PH_micro(1).gt.-0.8) then
                  ja=1
                  goto 203
               endif
!if one big root approach is taken:
!
!calculate the radius r at which the real centerpoint is located. Do not take the z-height
!into account. Find at which position in the radius array this value matches (sort of)
!take that index value in PH_micro and assign that PH as interface PH of that center point of a root segment.
           200 IF (ave) THEN
                  pos_cent=sqrt((cent(irecn,1,isub,ipl)-cp_mean(irecn,1,isub,ipl))*(cent(irecn,1,isub,ipl)-cp_mean(irecn,1,isub,ipl))&
 + (cent(irecn,2,isub,ipl)-cp_mean(irecn,2,isub,ipl))*(cent(irecn,2,isub,ipl)-cp_mean(irecn,2,isub,ipl))) !calculate radius from original cp to center of a voxel in 2D plane
                  DO i=1,size(r)-1
                     IF (pos_cent .eq. 0_sp) THEN
                        pos=1
                        EXIT
                     ELSEIF (pos_cent .le. r(i+1) .and. pos_cent .ge. r(i)) THEN
                        pos=i+1
                        EXIT
                     ENDIF
                  ENDDO
                  if (checker1.eq.1) then
                     PH_root_sub(irecn,isub) = PH_micro2(irecn,pos,isub)
                  else
!denote PH_int to all root nodes within this voxel independent of distance
                     PH_root_sub(irecn,isub) = PH_micro2(irecnTemp,pos,isubTemp)/pond_subTemp*w_sub(irecn,isub,ipl) !take difference in w_sub between cp's into account
                     if (segsur(irecnTemp) .eq. 0 .or. seglen(irecnTemp) .eq. 0) then
                        r0 = 5e-2 !minimal root radius (in cm)
                     else
                        r0 = 1/(2*pi)*segsur(irecnTemp)/seglen(irecnTemp)
                     end if
                     IF (k_ave(irecnTemp,isubTemp)*B(irecnTemp,isubTemp) .lt. Lr(irecnTemp)*r0) then
                        countja=.true.  !to check how many root nodes are under "limiting conditions"
                     endif
                  endif
               ENDIF
               IF (.not.(ave)) THEN
                  PH_root_sub(irecn,isub) = PH_micro(1)
               ENDIF
!-----------------------------------------------------------------------
! Give every soil-root interface node its correct PH value
!-----------------------------------------------------------------------
               PH_root_sub(irecn,isub) = PH_root_sub(irecn,isub) + zCent*w_sub(irecn,isub,ipl)!calc PH at interface with gravity
            else!checker !some segments are divided in subsegments in which a cp is located at the edge of the plane (and/or cube),use average approach of original R-SWMS
               sumW=sum_dis(irecn,isub,ipl)
               psub=w_sub(irecn,isub,ipl)
               PH_root_sub(irecn,isub) = ((hnew(corner(1))+zgrid(corner(1)))*pond_i(1)+&
                          (hnew(corner(2))+zgrid(corner(2)))*pond_i(2)+(hnew(corner(3))+&
               zgrid(corner(3)))*pond_i(3)+(hnew(corner(4))+zgrid(corner(4)))*pond_i(4)+&
                          (hnew(corner(5))+zgrid(corner(5)))*pond_i(5)+(hnew(corner(6))+&
               zgrid(corner(6)))*pond_i(6)+(hnew(corner(7))+zgrid(corner(7)))*pond_i(7)+&
                          (hnew(corner(8))+zgrid(corner(8)))*pond_i(8))/sumW*psub
            end if
!---------------------------------------------------------------------------
! ph at interface
            PHs(irecn,ipl)= PHs(irecn,ipl) + PH_root_sub(irecn,isub)
203         IF (ja.eq.1) then !average method used in some cases
               if (irecn .eq. 0) then
                  r0=5e-2
               else
                  r0 = 1/(2*pi)*segsur(irecn)/seglen(irecn)
               endif
               temp=Lr(irecn)*r0
               do j=1,8
                  h_matK(j)=hnew(corner(j))
                  call SetMat_singlevalue(h_matK(j),theta_k,K_k,C_k)
                  temp2(j)=K_k
               enddo
               if (sum(temp2)/8 .lt. temp) then
                  countja=.true.
               endif
               ja=0
               sumW=sum_dis(irecn,isub,ipl)
               psub=w_sub(irecn,isub,ipl)
               PH_root_sub(irecn,isub)=((hnew(corner(1))+zgrid(corner(1)))*pond_i(1)+&
                          (hnew(corner(2))+zgrid(corner(2)))*pond_i(2)+(hnew(corner(3))+&
               zgrid(corner(3)))*pond_i(3)+(hnew(corner(4))+zgrid(corner(4)))*pond_i(4)+&
                          (hnew(corner(5))+zgrid(corner(5)))*pond_i(5)+(hnew(corner(6))+&
               zgrid(corner(6)))*pond_i(6)+(hnew(corner(7))+zgrid(corner(7)))*pond_i(7)+&
                          (hnew(corner(8))+zgrid(corner(8)))*pond_i(8))/sumW*psub
               PHs(irecn,ipl)=PHs(irecn,ipl)+PH_root_sub(irecn,isub)
            ENDIF
         END DO!end do-loop over subsegments
 ! calculate matrix Q (doussan matrix)
         Qd(irecn,ipl)=Qi(irecn,ipl)*PHs(irecn,ipl)+Q_bc(irecn,ipl)
         if (countja) count_nodes=count_nodes+1
         countja=.false.
      END DO
      END SUBROUTINE analytical_approach
!=====================================================================================
               SUBROUTINE flux_node(hNew,Flow_corner,corner,irecn,isub,q_tot_c)
               USE typedef
               use paramdata, only: maxnod
               USE GridData, only: xgrid,ygrid,zgrid,dxgrid,dygrid,dzgrid
          !     USE DoussanMAt, only:
               implicit none
               REAL,intent(in) :: hNew(maxnod)
               integer(sp) :: i,j,Flow_corner(1:8,1:6)!,temp(1:6),v
               REAL(dp) :: q_tot_c_t
               REAL(sp) :: K_rn
               REAL(sp) :: theta_single,K_single,C_single,Kn,Kcn,Kc,dz
               INTEGER(sp) ::corner(1:8),irecn,isub
               real(dp),intent(out) :: q_tot_c(1:8)
               q_tot_c(1:8) = 0
              ! loop over all soil nodes from cubic (sub cubic -> subsegment)
               do i=1,8
               !calculate soil hydraulic conductivity in the corner node and the surrounding node (in next loop)
                  call SetMat_singlevalue(hnew(corner(i)),theta_single,K_single,C_single)
                  Kc = K_single
                  !loop over all surrounding nodes
!                  do j=1,v-1
                  do j=1,6
!                  call SetMat_singlevalue(hnew(temp(j)),theta_single,K_single,C_single)
                     if (Flow_corner(i,j).eq.0) then
                        q_tot_c_t = 0
                        goto 40
                     endif
                     call SetMat_singlevalue(hnew(Flow_corner(i,j)),theta_single,K_single,C_single)
                     Kn = K_single
                     Kcn = (Kc+Kn)/2._sp !arithmetic mean soil hydraulic conductivity
                     !distance between the two nodes evaluated
                     if (xgrid(corner(i))-xgrid(Flow_corner(i,j)) .ne. 0) then
                        dz=dxgrid
                     else if (ygrid(corner(i))-ygrid(Flow_corner(i,j)) .ne. 0) then
                        dz=dygrid
                     else if (zgrid(corner(i))-zgrid(Flow_corner(i,j)) .ne. 0) then
                        dz=dzgrid
                        ! -> q_tot_c_t = 0 (no flow from z-direction taken into account
                     end if
                     !flow from corner node to surrounding node
                     q_tot_c_t = - Kcn * (hnew(corner(i)) - hnew(Flow_corner(i,j)))/dz
                     !take the sum over all the calculates flows from and to the corner node
                     !q_tot_c has size 8 (corner nodes of cube)
                  40 continue
!                    write(*,*)'q_tot_c_t=',q_tot_c_t
                     q_tot_c(i) = q_tot_c(i) + q_tot_c_t
                  end do
               end do
               return
               end subroutine flux_node
!=====================================================================================
   subroutine average_method(hNew,ipl)
     use typedef
     USE Paramdata,only: MaxNod,pi
     use doussanmat, only: PHs,nsub,w_sub,sum_dis,PH_root_sub,Qi,Qd,Q_bc,w_dis,loc_Q,Lr,count_nodes
     USE GridData, only: zgrid
     USE RootData, only: nrec,segsur,seglen
     implicit none
     integer(sp) :: isub,i,j
     REAL,intent(in) :: hNew(maxnod)
     REAL (dp) ::sumW,psub
     REAL(dp),dimension(1:8) ::pond_i
     REAL(sp),dimension(1:8) ::h_mat,temp2
     real(sp) :: K,C,theta,r0,temp
     INTEGER(sp) ::corner(1:8)
     INTEGER(sp) :: irecn,ipl
     logical :: countja

     countja=.false.
     count_nodes=0

     DO irecn=0,nrec
        if (irecn .eq. 0) then
           r0=5e-2
        else
           r0 = 1/(2*pi)*segsur(irecn)/seglen(irecn)
        endif
        temp=Lr(irecn)*r0
!initialize soil PH matrix
        PHs(irecn,ipl)=0._dp
        DO isub=1,nsub(irecn,ipl)
           DO i=1,8
              pond_i(i)=w_dis(irecn,i,isub,ipl)
              corner(i)=loc_Q(irecn,i,isub,ipl)
           END DO
           do j=1,8
              h_mat(j)=hnew(corner(j))
              call SetMat_singlevalue(h_mat(j),theta,K,C)
              temp2(j)=K
           enddo
           if (sum(temp2)/8 .lt. temp) then
              countja=.true.
           endif
           sumW=sum_dis(irecn,isub,ipl)!total inverse distance to node 
           psub=w_sub(irecn,isub,ipl)!subsegment/total segment ratio
           PH_root_sub(irecn,isub)=((hnew(corner(1))+zgrid(corner(1)))*pond_i(1)+&
                      (hnew(corner(2))+zgrid(corner(2)))*pond_i(2)+(hnew(corner(3))+&
           zgrid(corner(3)))*pond_i(3)+(hnew(corner(4))+zgrid(corner(4)))*pond_i(4)+&
                      (hnew(corner(5))+zgrid(corner(5)))*pond_i(5)+(hnew(corner(6))+&
           zgrid(corner(6)))*pond_i(6)+(hnew(corner(7))+zgrid(corner(7)))*pond_i(7)+&
                      (hnew(corner(8))+zgrid(corner(8)))*pond_i(8))/sumW*psub
	    PHs(irecn,ipl)=PHs(irecn,ipl)+PH_root_sub(irecn,isub)    
        END DO    !end do-loop over subsegments
! calculate matrix Qd
        Qd(irecn,ipl)=Qi(irecn,ipl)*PHs(irecn,ipl)+Q_bc(irecn,ipl)
        if (countja) count_nodes=count_nodes + 1
           countja=.false.
        END DO !end do-loop over irecn
        end subroutine average_method
!===================================================================================
!***********************************************************************************
      SUBROUTINE logspace(d1,d2,n,y)
      USE typedef
      IMPLICIT NONE
      integer(dp),intent(in) :: n
      integer(dp) :: i
      real(sp) :: qq(1:n),Y_(1:n)
      real(sp) :: d1, d2, d31, d32
      real(sp), intent(out) :: y(1:n)
      !logspace algorithm
      d31 = log10(abs(d1))
      d32 = log10(abs(d2))
!write(*,*)'d31',d31
!write(*,*)'d32',d32
!$OMP PARALLEL DO
      do i=1,n-1
         qq(i) =d31+(i-1)*(d32-d31)/(n-1)
      end do
!$OMP END PARALLEL DO
!pause
      !put in one array
      y_(1:n-1) = qq(1:n-1)
      y_(n) = d32
      !convert logscale to normal values
      y = (10)**y_(1:n)
      return
      end subroutine logspace
!***********************************************************************
      SUBROUTINE SetMat_ana(hTemp,theta,K,C)
      USE GridData
      USE SolData
      USE WatFun
      USE Doussanmat, only : n
      IMPLICIT NONE
      REAL(sp), intent(in) :: hTemp(1:n)
      REAL(sp):: Ti(1:n),t0,t1
      REAL(sp), intent(out) :: theta(1:n),C(1:n),K(1:n)
      INTEGER(sp):: i,M
!      call cpu_time(t0)
      M=nMat
!$OMP PARALLEL DO
      DO i=1,n
   !      write(*,*)'I am thread',omp_get_thread_num(),' of ', omp_get_num_threads(),'threads'
! nodal conductivity values:
         K(i)=FKP(hTemp(i),par(:,M))*Bxy(i)
         C(i)=FCP(hTemp(i),par(:,M))*Dxy(i)/Axy(i)
	  Ti(i)=Fth(hTemp(i),par(:,M))
         theta(i)=par(2,M)*Exy(i)+(Ti(i)-par(2,M))*Dxy(i) !!to be checked (MJ09)
      ENDDO
      RETURN
      END SUBROUTINE SetMat_ana
!***********************************************************************
      SUBROUTINE SetMat_anaLibr(hTemp,theta,K,C)
      USE GridData
      USE SolData
      USE WatFun
      USE Doussanmat, only: nLibr
      IMPLICIT NONE
      REAL(sp), intent(in) :: hTemp(1:nLibr)
      REAL(sp):: Ti(1:nLibr)
      REAL(sp), intent(out) :: theta(1:nLibr),C(1:nLibr),K(1:nLibr)
      INTEGER(sp):: i,M
      M=nMat
!$OMP PARALLEL DO
      DO i=1,nLibr
! nodal conductivity values:
         K(i)=FKP(hTemp(i),par(:,M))
         C(i)=FCP(hTemp(i),par(:,M))
	  Ti(i)=Fth(hTemp(i),par(:,M))
         theta(i)=par(2,M)+(Ti(i)-par(2,M))  !!to be checked (MJ09)
      enddo
!$OMP END PARALLEL DO
      RETURN
      END
!***********************************************************************
      SUBROUTINE SetMat_singlevalue(hTemp,theta,K,C)
      USE GridData
      USE SolData
      USE WatFun
      IMPLICIT NONE
      REAL(sp), intent(in) :: hTemp
      REAL(sp):: Ti
      REAL(sp), intent(out) :: theta,C,K
      INTEGER(sp):: i,M
      M=nMat
! nodal conductivity values:
      K=FKP(hTemp,par(:,M))
      C=FCP(hTemp,par(:,M))
      Ti=Fth(hTemp,par(:,M))
      theta=par(2,M)+(Ti-par(2,M))  !!to be checked (MJ09)
      RETURN
      END SUBROUTINE SetMat_singlevalue
!**************************************************************************************
      SUBROUTINE SetSnk(hNew,Conc,lDou)
      USE TypeDef
      USE Paramdata,only: MaxNod
      USE TempData, only: tem
      USE PlntData, only: Tpot,CMm,VMax,fk,xin,TotSur
      USE GridData, only: nPt,betaw,betac,sink,csink,xgrid,ygrid,zgrid
      IMPLICIT NONE
      REAL,intent(in) :: Conc(maxnod),hNew(maxnod)
      LOGICAL ::lDou
      REAL(sp) :: alpha,active
      INTEGER ::i

      DO i=1,nPt
         IF (.not.(lDou)) THEN
            IF (betaw(i).GT.1.E-20_dp) THEN
               sink(i)=alpha(hNew(i),Conc(i),tem(i))*betaw(i)*Tpot
            ELSE
               sink(i)=0.0_dp
            ENDIF
         ENDIF
         IF (betac(i).GT.1.E-20_dp) THEN
            active=TotSur*betac(i)*(VMax/(CMm+Conc(i))+fk)
         ELSE
            active=0.0_dp
         END IF
         csink(i)=xin*sink(i)+(1-xin)*active
      END DO

      RETURN
      END SUBROUTINE SetSnk
!**************************************************************************************
      FUNCTION alpha(h,Conc,tem)
      USE PlntData
      IMPLICIT NONE
      REAL(sp), intent(in) ::conc,tem,h
      REAL(sp):: alpha,alphap,alphah,p,rktemp
!!!!  Calculate osmotic potential [L] from soil solution concentration [M/L3].
!     To use the first or the second method simply render active the lines
!     where necessary quantities are calculated for the chosen expression and
!    comment out the lines related to the other expression; recompile the
!     program; beware of units and, if needed, make proper aggiustments.
!     Also, adjust for the valency of the chemical species studied if the
!     second method is adopted.
!**************************
!     For very dilute solutions the van't Hoff expression can be used (T.L.
!     Allen and R.M. Keefer, "Chemistry: Experiment and Theory", 1974,
!     Harper and Row Publishers, pp 283-285):
!                              p=CRT
!     where:
!     p is the osmotic potential [L], C is the concentration [M/L3], R is a
!     constant [(L3*pressure)/(M*temperature)], and T is the temperature
!     (K degrees)
!
!     first transform temperature degrees from Celsius to Kelvin...
!
        rKtemp=tem+273.15_dp
!
!     1 liter = 1000 cm3; 1 atm = 1033.6 cm H2O;
!     R = 0.082 (liter atm K-1 mole-1) = 84811.0144 (cm3 cm K-1 mole-1);
!
        p=(-Conc)*rKtemp*84811.0144_dp
!
!**********************************
!     The following constants are taken from "Diagnosis and Improvement of
!     Saline and Alkali Soils", USDA Agriculture Handbook no. 60, 1954:
!               1 atm = 0.36 millimhos/cm           (Fig 6, p 15);
!               1 millimhos/cm = -0.091 meq/liter   (Fig 4, p 12);
!     these values are valid for the range of EC that allows plant growth.
!    Therefore 1 atm = -0.3276 meq/liter; also:
!               1 atm = 1033.6 cm H2O;     and:
!               1 meq/liter = 1.0E+06/n M/cm3, being n the valency of the      ion
!                                                        in the soil solution.
!
!      p=Conc*(-0.3276)*1000000/1033.6
!
!********************************************************************************
      alpha=0.0_dp
      IF (h50.LT.9.E+20_dp) THEN
         alphah=1._dp/(1._dp+(h/h50)**p1)
         IF (p50.LT.9.E+20_dp) THEN
            alphap=1._dp/(1._dp+(p+p50)**p2)
            alpha=alphah*alphap
         ELSE
            alpha=alphah
         ENDIF
      ELSE
         IF ((h.GT.h3).AND.(h.LT.h2)) alpha=(h-h3)/(h2-h3)
         IF ((h.GE.h2).AND.(h.LE.h1)) alpha=1.0_dp
         IF ((h.GT.h1).AND.(h.LT.h0)) alpha=(h-h0)/(h1-h0)
      ENDIF
      RETURN
      END FUNCTION alpha
!********************************************************************************
! weigth for each node in fucntion of their neighbouors
      SUBROUTINE CalcWnodes
      USE GridData
      IMPLICIT NONE
      INTEGER(sp) :: iE,ic
! calculate actual overall transpiration rate from soil domain:
      Wn=0.0_dp
      DO 1 iE=1,nElm,2
! assign to Wn the number of times a cuboid corner node is called
       Wn(elmnod(iE,1))=Wn(elmnod(iE,1))+1
       Wn(elmnod(iE+1,6))=Wn(elmnod(iE+1,6))+1
       Wn(elmnod(iE,3))=Wn(elmnod(iE,3))+1
       Wn(elmnod(iE,2))=Wn(elmnod(iE,2))+1
       Wn(elmnod(iE,4))=Wn(elmnod(iE,4))+1
       Wn(elmnod(iE+1,3))=Wn(elmnod(iE+1,3))+1
       Wn(elmnod(iE,6))=Wn(elmnod(iE,6))+1
       Wn(elmnod(iE,5))=Wn(elmnod(iE,5))+1
    1 CONTINUE
      Wn=Wn/8._dp
      RETURN
      END SUBROUTINE CalcWnodes
! ==============================================================================
! Source file SOLUTE |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
      SUBROUTINE SOLUTE(dt,t,Conc,tPulse,hOld,hNew,dtMaxC,MatNum,Kode,KodCB,theta,Dxy,Exy,theta_old)
      USE BoundData
      USE GridData
      USE SolData
      USE WatFun
      USE CumData
      USE MatData, only: As,Bs
      IMPLICIT NONE

      REAL(sp),intent(in)::dt,t,tpulse,dtmaxc,hNew(maxnod),hOld(maxnod)
      REAL,intent(inout) :: Conc(maxnod)
      REAL(sp) :: theta(maxnod),theta_old(maxnod)
      INTEGER(sp),intent(in)::MatNum(maxnod),Kode(maxnod),KodCB(maxbdr)
      REAL(sp) :: VxE(4),VyE(4),VzE(4),Ac(maxnod),S(4,4),DS(maxnod)
      INTEGER(sp):: i,j,k,l,m,ie,ise,ic,j1,i1,j2,i2,ib,lev,iL(4),List(4)
      REAL(sp) :: Vzz,Vyy,Vxx,Al,Ak,Ai,Aj,aa,Exy(Maxnod),Dxy(MaxNod)
      REAL(sp) :: ec1,gce,fmul,rootch,vzee,vyee,vxee,cone,VEl,det,cayz
      REAL(sp) :: smul2,smul1,ace,fce,ec6,ec5,ec4,ec3,ec2,cazz,alf,cayy
      REAL(sp) :: dpom,caxz,caxy,caxx,bi(4),ci(4),di(4),F(maxnod)
!**********************************************************************************
! Initialization
      alf=1._dp-epsi
      IF (t.GT.tPulse) THEN
         DO 11 i=1,2
            cBound(i)=0.0_dp
11       CONTINUE
      ENDIF
      DO 13 i=1,Npt
         Bs(i) =0.0_dp
         Qc(i)=0.0_dp
         IF (epsi.LT.0.001_dp) THEN
            As(Maxbnd,i)=0.0_dp
         ELSE
            DO 12 j=1,2*Nband-1
               As(j,i)=0.0_dp
12          CONTINUE
         ENDIF
13    CONTINUE

      DO 21 Lev=1,NLevel
         IF (Lev.EQ.NLevel) THEN
            CALL Veloc(hNew)
            CALL Disper(MatNum,hNew,theta,Dxy,Exy)
            CALL PeCour(MatNum,dtMaxC,dt,hNew,theta)
         ELSE
            CALL Disper(MatNum,hNew,theta,Dxy,Exy)
         ENDIF
         DO 14 i=1,Npt
            M=MatNum(i)
            IF (Lev.NE.NLevel) THEN
             DPom=dt/6._dp/(theta_old(i)+ChPar(1,M)*ChPar(5,M))
               Dispxx(i)=Dispxx(i)+Vx(i)*Vx(i)*DPom
               Dispyy(i)=Dispyy(i)+Vy(i)*Vy(i)*DPom
               Dispzz(i)=Dispzz(i)+Vz(i)*Vz(i)*DPom
               Dispxy(i)=Dispxy(i)+Vx(i)*Vy(i)*DPom
               Dispxz(i)=Dispxz(i)+Vx(i)*Vz(i)*DPom
               Dispyz(i)=Dispyz(i)+Vy(i)*Vz(i)*DPom
            ELSE
               Ac(i)=-(theta_old(i)*alf+theta(i)*epsi)-ChPar(1,M)*ChPar(5,M)
             DPom=dt/6._dp/(theta(i)+ChPar(1,M)*ChPar(5,M))
               Dispxx(i)=Dispxx(i)-Vx(i)*Vx(i)*DPom
               Dispyy(i)=Dispyy(i)-Vy(i)*Vy(i)*DPom
               Dispzz(i)=Dispzz(i)-Vz(i)*Vz(i)*DPom
               Dispxy(i)=Dispxy(i)-Vx(i)*Vy(i)*DPom
               Dispxz(i)=Dispxz(i)-Vx(i)*Vz(i)*DPom
               Dispyz(i)=Dispyz(i)-Vy(i)*Vz(i)*DPom
               Gc(i)=ChPar(8,M)*theta(i)+ChPar(1,M)*ChPar(9,M)
               Fc(i)=ChPar(6,M)*theta(i)+ChPar(1,M)*ChPar(7,M)*ChPar(5,M)+sink(i)-csink(i)
            ENDIF
14       CONTINUE
         DO 15 i=1,Npt
            F(i)=0.0_sp
            IF (Lev.EQ.NLevel) DS(i)=0.0_sp
15       CONTINUE
! Loop on elements
         DO 19 iE=1,Nelm
          CAxx=ConAxx(iE)
          CAyy=ConAyy(iE)
          CAzz=ConAzz(iE)
          CAxy=ConAxy(iE)
          CAxz=ConAxz(iE)
          CAyz=ConAyz(iE)
! Loop on subelements
          DO 18 iSE=1,3
             iL(1)=3+iSE
             iL(2)=4+iSE-iSE/3*6
             iL(3)=5+iSE-iSE/2*6
             iL(4)=iSE
            i=elmnod(iE,iL(1))
            j=elmnod(iE,iL(2))
            k=elmnod(iE,iL(3))
            l=elmnod(iE,iL(4))
              List(1)=i
              List(2)=j
              List(3)=k
              List(4)=l
           bi(1)=-(yGrid(k)-yGrid(j))*(zGrid(l)-zGrid(j))+(yGrid(l)-yGrid(j))*(zGrid(k)-zGrid(j))
           bi(2)=+(yGrid(l)-yGrid(k))*(zGrid(i)-zGrid(k))-(yGrid(i)-yGrid(k))*(zGrid(l)-zGrid(k))
           bi(3)=-(yGrid(i)-yGrid(l))*(zGrid(j)-zGrid(l))+(yGrid(j)-yGrid(l))*(zGrid(i)-zGrid(l))
           bi(4)=+(yGrid(j)-yGrid(i))*(zGrid(k)-zGrid(i))- (yGrid(k)-yGrid(i))*(zGrid(j)-zGrid(i))
           ci(1)=+(xGrid(k)-xGrid(j))*(zGrid(l)-zGrid(j))-(xGrid(l)-xGrid(j))*(zGrid(k)-zGrid(j))
           ci(2)=-(xGrid(l)-xGrid(k))*(zGrid(i)-zGrid(k))+(xGrid(i)-xGrid(k))*(zGrid(l)-zGrid(k))
           ci(3)=+(xGrid(i)-xGrid(l))*(zGrid(j)-zGrid(l))-(xGrid(j)-xGrid(l))*(zGrid(i)-zGrid(l))
           ci(4)=-(xGrid(j)-xGrid(i))*(zGrid(k)-zGrid(i))+(xGrid(k)-xGrid(i)) *(zGrid(j)-zGrid(i))
           di(1)=-(xGrid(k)-xGrid(j))*(yGrid(l)-yGrid(j))+(xGrid(l)-xGrid(j))*(yGrid(k)-yGrid(j))
           di(2)=+(xGrid(l)-xGrid(k))*(yGrid(i)-yGrid(k))-(xGrid(i)-xGrid(k))*(yGrid(l)-yGrid(k))
           di(3)=-(xGrid(i)-xGrid(l))*(yGrid(j)-yGrid(l))+(xGrid(j)-xGrid(l))*(yGrid(i)-yGrid(l))
           di(4)=+(xGrid(j)-xGrid(i))*(yGrid(k)-yGrid(i))-(xGrid(k)-xGrid(i))*(yGrid(j)-yGrid(i))
           Det=(xGrid(l)-xGrid(i))*bi(4)+(yGrid(l)-yGrid(i))*ci(4)+(zGrid(l)-zGrid(i))*di(4)
           VEl=abs(Det)/6._sp
!           Calculate Velocities
            AA=1._sp/Det
            Ai=CAxx*Bi(1)+CAxy*Ci(1)+CAxz*Di(1)
            Aj=CAxx*Bi(2)+CAxy*Ci(2)+CAxz*Di(2)
            Ak=CAxx*Bi(3)+CAxy*Ci(3)+CAxz*Di(3)
            Al=CAxx*Bi(4)+CAxy*Ci(4)+CAxz*Di(4)
            if(Lev.eq.NLevel) then
              Vxx=AA*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k)+Al*hNew(l))+CAxz
            else
              Vxx=AA*(Ai*hOld(i)+Aj*hOld(j)+Ak*hOld(k)+Al*hOld(l))+CAxz
            end if
            Ai=CAxy*Bi(1)+CAyy*Ci(1)+CAyz*Di(1)
            Aj=CAxy*Bi(2)+CAyy*Ci(2)+CAyz*Di(2)
            Ak=CAxy*Bi(3)+CAyy*Ci(3)+CAyz*Di(3)
            Al=CAxy*Bi(4)+CAyy*Ci(4)+CAyz*Di(4)
            if(Lev.eq.NLevel) then
              Vyy=AA*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k)+Al*hNew(l))+CAyz
            else
              Vyy=AA*(Ai*hOld(i)+Aj*hOld(j)+Ak*hOld(k)+Al*hOld(l))+CAyz
            end if
            Ai=CAxz*Bi(1)+CAyz*Ci(1)+CAzz*Di(1)
            Aj=CAxz*Bi(2)+CAyz*Ci(2)+CAzz*Di(2)
            Ak=CAxz*Bi(3)+CAyz*Ci(3)+CAzz*Di(3)
            Al=CAxz*Bi(4)+CAyz*Ci(4)+CAzz*Di(4)
            if(Lev.eq.NLevel) then
              Vzz=AA*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k)+Al*hNew(l))+CAzz
            else
              Vzz=AA*(Ai*hOld(i)+Aj*hOld(j)+Ak*hOld(k)+Al*hOld(l))+CAzz
            end if
            if(Lev.ne.NLevel) then
              ConE=(ConO(i)+ConO(j)+ConO(k)+ConO(l))/4._dp
              VxE(1)=-ConO(i)*Vxx
              VxE(2)=-ConO(j)*Vxx
              VxE(3)=-ConO(k)*Vxx
              VxE(4)=-ConO(l)*Vxx
              VyE(1)=-ConO(i)*Vyy
              VyE(2)=-ConO(j)*Vyy
              VyE(3)=-ConO(k)*Vyy
              VyE(4)=-ConO(l)*Vyy
              VzE(1)=-ConO(i)*Vzz
              VzE(2)=-ConO(j)*Vzz
              VzE(3)=-ConO(k)*Vzz
              VzE(4)=-ConO(l)*Vzz
            else
              ConE=(Con(i)+Con(j)+Con(k)+Con(l))/4._dp
              VxE(1)=-Con(i)*Vxx
              VxE(2)=-Con(j)*Vxx
              VxE(3)=-Con(k)*Vxx
              VxE(4)=-Con(l)*Vxx
              VyE(1)=-Con(i)*Vyy
              VyE(2)=-Con(j)*Vyy
              VyE(3)=-Con(k)*Vyy
              VyE(4)=-Con(l)*Vyy
              VzE(1)=-Con(i)*Vzz
              VzE(2)=-Con(j)*Vzz
              VzE(3)=-Con(k)*Vzz
              VzE(4)=-Con(l)*Vzz
            end if
            VxEE=-ConE*Vxx
            VyEE=-ConE*Vyy
            VzEE=-ConE*Vzz

             IF (Lev.EQ.1) THEN
              RootCh=VEl*dt*(Conc(i)*csink(i)+Conc(j)*csink(j)+Conc(k)*csink(k)+Conc(l)*csink(l))/4._sp
              CumCh0=CumCh0-VEl*dt*(Gc(i)+Gc(j)+Gc(k)+Gc(l))/4._sp
              CumCh1=CumCh1-VEl*dt*((Fc(i)-sink(i))*Conc(i)+(Fc(j)-sink(j))*Conc(j)+&
                   (Fc(k)-sink(k))*Conc(k)+(Fc(l)-sink(l))*Conc(l))/4._sp-RootCh
              CumChR=CumChR+RootCh
             ENDIF
             FMul=VEl/5._dp
             GcE=(Gc(i)+Gc(j)+Gc(k)+Gc(l))/4._sp
             Ec1=(Dispxx(i)+Dispxx(j)+Dispxx(k)+Dispxx(l))/4._sp
             Ec2=(Dispyy(i)+Dispyy(j)+Dispyy(k)+Dispyy(l))/4._sp
             Ec3=(Dispzz(i)+Dispzz(j)+Dispzz(k)+Dispzz(l))/4._sp
             Ec4=(Dispxy(i)+Dispxy(j)+Dispxy(k)+Dispxy(l))/4._sp
             Ec5=(Dispxz(i)+Dispxz(j)+Dispxz(k)+Dispxz(l))/4._sp
             Ec6=(Dispyz(i)+Dispyz(j)+Dispyz(k)+Dispyz(l))/4._sp
             IF (Lev.EQ.NLevel) AcE=(Ac(i)+Ac(j)+Ac(k)+Ac(l))/4._sp
             FcE=(Fc(i)+Fc(j)+Fc(k)+Fc(l))/4._sp
             SMul1=-1._dp/VEl/36._sp
             SMul2=VEl/30._sp
             DO 17 j1=1,4
              i1=List(j1)
              F(i1)=F(i1)+FMul*(GcE+Gc(i1)/4._sp)
              IF(Lev.EQ.NLevel) DS(i1)=DS(i1)+FMul*(AcE+Ac(i1)/4._sp)
              DO 16 j2=1,4
               i2=List(j2)
               S(j1,j2)=SMul1*(Ec1*bi(j1)*bi(j2)+Ec2*ci(j1)*ci(j2)+ Ec3*di(j1)*di(j2)+&
                    Ec4*(bi(j1)*ci(j2)+ci(j1)*bi(j2))+Ec5*(bi(j1)*di(j2)+di(j1)*bi(j2))+&
                    Ec6*(ci(j1)*di(j2)+di(j1)*ci(j2)))
               S(j1,j2)=S(j1,j2)-(6*VEl/Det)*((bi(j2)/30._sp)*(VxEE+VxE(j1)/4._sp)+&
                    (ci(j2)/30._sp)*(VyEE+VyE(j1)/4._sp)+(di(j2)/30._sp)*(VzEE+VzE(j1)/4._sp))
               ic=1
               IF (i1.EQ.i2) ic=2
               S(j1,j2)=S(j1,j2)+SMul2*ic*(FcE+(Fc(i1)+Fc(i2))/4._sp)
               IF (Lev.NE.NLevel) THEN
                 Bs(i1)=Bs(i1)-alf*S(j1,j2)*Conc(i2)
               ELSE
                 iB=Nband+i2-i1
                 As(iB,i1)=As(iB,i1)+epsi*S(j1,j2)
               ENDIF
               IF (Lev.EQ.1.AND.Kode(i1).GT.0)Qc(i1)=Qc(i1)-S(j1,j2)*Conc(i2)
16                CONTINUE
17               CONTINUE
18          CONTINUE
19       CONTINUE
       DO 20 i=1,Npt
          M=MatNum(i)
          IF (Lev.EQ.1.AND.Kode(i).GT.0) Qc(i)=Qc(i)-F(i)
          IF (Lev.NE.NLevel) THEN
             Bs(i)=Bs(i)-alf*F(i)
          ELSE
             As(Nband,i)=As(Nband,i)+DS(i)/dt
             Bs(i)=Bs(i)+DS(i)/dt*Conc(i)-epsi*F(i)
          ENDIF
20       CONTINUE
21    CONTINUE
! Boundary condition
      CALL C_Bound(Kode,KodCB,Conc,dt,DS)
! Solve the global matrix equation for transport
      IF (epsi.LT.0.001_sp) THEN
         DO 22 i=1,Npt
            Bs(i)=Bs(i)/As(Nband,i)
22       CONTINUE
      ELSE
         CALL SolveT
      ENDIF
      DO 23 i=1,Npt
         Conc(i)=Bs(i)!sngl(B(i))
         IF (Conc(i).LT.1.0E-30_dp) Conc(i)=0.0_dp
23    CONTINUE
      CALL SolInf(dt,Kode)
      RETURN
      END SUBROUTINE SOLUTE
!*******************************************************************
      SUBROUTINE C_Bound(Kode,KodCB,Conc,dt,DS)
      USE BoundData
      USE GridData
      USE SolData
      USE CumData
      USE MatData
      IMPLICIT NONE
      REAL(sp),intent(in) :: dt,DS(maxnod)
      REAL,intent(inout) :: Conc(maxnod)
      INTEGER(sp),intent(inout) :: Kode(maxnod),KodCB(maxbdr)
      REAL(sp) :: alf,cBnd
      INTEGER(sp) :: cKod,i,j
      alf=1._dp-epsi
      DO 14 i=1,nPt
       IF (Kode(i).NE.0) THEN
          DO 11 j=1,nBCPts
             IF (iBCPt(j).EQ.i) THEN
                  IF (KodCB(j).GT.0) THEN
                       cKod=1
                     cBnd=cBound(KodCB(j))
                  ELSE
                       IF (Q(i).GT.0._dp) THEN
                                cKod=3
                        cBnd=cBound(-KodCB(j))
                         ELSE
                      cKod=2
                   ENDIF
                  ENDIF
                  GOTO 12
             ENDIF
11        CONTINUE
12        CONTINUE
! Dirichlet boundary condition
          IF (cKod.EQ.1) THEN
             Qc(i)=Qc(i)+Q(i)*(epsi*cBnd+alf*Conc(i))-DS(i)*(cBnd-Conc(i))/dt
             DO 13 j=1,2*NBand-1
                  A(j,i)=0.0_dp
13           CONTINUE
             A(NBand,i)=1._dp
             B(i)=cBnd
          ENDIF
! Neumann boundary condition
          IF (cKod.EQ.2) THEN
             Qc(i)=Q(i)*Conc(i)
          ENDIF
! Cauchy boundary condition
          IF (cKod.EQ.3) THEN
             B(i)=B(i)-Q(i)*(cBnd-alf*Conc(i))
             A(NBand,i)=A(NBand,i)-epsi*Q(i)
             Qc(i)=Q(i)*cBnd
          ENDIF
       ENDIF
14    CONTINUE
      RETURN
      END SUBROUTINE C_Bound
!**********************************************************************************
      SUBROUTINE ChInit(MatNum,hNew,dtMaxC,dt,theta,Dxy,Exy)
      USE GridData
      USE SolData
      USE WatFun
      IMPLICIT NONE
      INTEGER(sp),intent(in) ::  MatNum(maxnod)
      REAL(sp),intent(in) :: dt,hNew(maxnod),theta(maxnod),Dxy(maxnod),Exy(maxnod)
      REAL,intent(inout) :: dtMaxC
      INTEGER(sp) ::  i,M
      DO 11 i=1,nPt
        IF (NLevel.EQ.2) THEN
          M=MatNum(i)
          Gc(i)=ChPar(8,M)*theta(i)+ChPar(1,M)*ChPar(9,M)
          Fc(i)=ChPar(6,M)*theta(i)+ChPar(1,M)*ChPar(7,M)*ChPar(5,M)+sink(i)-csink(i)
        ENDIF
11    CONTINUE
      CALL Veloc(hNew)
      CALL Disper(MatNum,hNew,theta,Dxy,Exy)
      CALL PeCour(MatNum,dtMaxC,dt,hNew,theta)
      RETURN
      END SUBROUTINE ChInit
!*************************************************************************
      SUBROUTINE Veloc(hNew)
      USE GridData
      USE SolData
      USE CumData
      IMPLICIT NONE
      REAL(sp),intent(in) :: hNew(maxnod)
      INTEGER(sp) :: i,j,k,l,m,ise,ie,iL(4),List(4)
      REAL(sp) :: det,dl,dk,di,dj,ci,cj,ck,cl,bi,bj,bk,bl,ai,aj,ak,al
      REAL(sp) :: vxx,vyy,vzz,cayz,caxz,caxy,cazz,cayy,caxx
      DO 11 i=1,nPt
       Vx(i)=0.0_dp
       Vy(i)=0.0_dp
       Vz(i)=0.0_dp
11    CONTINUE
      DO 14 iE=1,nElm
       CAxx=ConAxx(iE)
       CAyy=ConAyy(iE)
       CAzz=ConAzz(iE)
       CAxy=ConAxy(iE)
       CAxz=ConAxz(iE)
       CAyz=ConAyz(iE)
       DO 13 iSE=1,3
          IF (mod(iE,2).eq.0) THEN
             iL(1)=3+iSE
             iL(2)=4+iSE-iSE/3*6
             iL(3)=5+iSE-iSE/2*6
             iL(4)=  iSE
          ELSEIF (mod(iE,2).eq.1) THEN
             iL(1)=3+iSE
             iL(2)=5+iSE-iSE/2*6
             iL(3)=4+iSE-iSE/3*6
             iL(4)=  iSE
          ENDIF
          i=elmnod(iE,iL(1))
          j=elmnod(iE,iL(2))
          k=elmnod(iE,iL(3))
          l=elmnod(iE,iL(4))
!          iL(1)=3+iSE
!          iL(2)=4+iSE-iSE/3*6
!          iL(3)=5+iSE-iSE/2*6
!          iL(4)=iSE
!          i=elmnod(iE,iL(1))
!          j=elmnod(iE,iL(2))
!          k=elmnod(iE,iL(3))
!          l=elmnod(iE,iL(4))
          List(1)=i
          List(2)=j
          List(3)=k
          List(4)=l
         bi=-(yGrid(k)-yGrid(j))*(zGrid(l)-zGrid(j))+(yGrid(l)-yGrid(j))*(zGrid(k)-zGrid(j))
         bj=+(yGrid(l)-yGrid(k))*(zGrid(i)-zGrid(k))-(yGrid(i)-yGrid(k))*(zGrid(l)-zGrid(k))
         bk=-(yGrid(i)-yGrid(l))*(zGrid(j)-zGrid(l))+(yGrid(j)-yGrid(l))*(zGrid(i)-zGrid(l))
         bl=+(yGrid(j)-yGrid(i))*(zGrid(k)-zGrid(i))-(yGrid(k)-yGrid(i))*(zGrid(j)-zGrid(i))
         ci=+(xGrid(k)-xGrid(j))*(zGrid(l)-zGrid(j))-(xGrid(l)-xGrid(j))*(zGrid(k)-zGrid(j))
         cj=-(xGrid(l)-xGrid(k))*(zGrid(i)-zGrid(k))+(xGrid(i)-xGrid(k))*(zGrid(l)-zGrid(k))
         ck=+(xGrid(i)-xGrid(l))*(zGrid(j)-zGrid(l))-(xGrid(j)-xGrid(l))*(zGrid(i)-zGrid(l))
         cl=-(xGrid(j)-xGrid(i))*(zGrid(k)-zGrid(i))+(xGrid(k)-xGrid(i))*(zGrid(j)-zGrid(i))
         di=-(xGrid(k)-xGrid(j))*(yGrid(l)-yGrid(j))+(xGrid(l)-xGrid(j))*(yGrid(k)-yGrid(j))
         dj=+(xGrid(l)-xGrid(k))*(yGrid(i)-yGrid(k))-(xGrid(i)-xGrid(k))*(yGrid(l)-yGrid(k))
         dk=-(xGrid(i)-xGrid(l))*(yGrid(j)-yGrid(l))+(xGrid(j)-xGrid(l))*(yGrid(i)-yGrid(l))
         dl=+(xGrid(j)-xGrid(i))*(yGrid(k)-yGrid(i))-(xGrid(k)-xGrid(i))*(yGrid(j)-yGrid(i))
         Det=(xGrid(l)-xGrid(i))*bl+(yGrid(l)-yGrid(i))*cl+(zGrid(l)-zGrid(i))*dl
          Ai=CAxx*bi+CAxy*ci+CAxz*di
          Aj=CAxx*bj+CAxy*cj+CAxz*dj
          Ak=CAxx*bk+CAxy*ck+CAxz*dk
          Al=CAxx*bl+CAxy*cl+CAxz*dl
            Vxx=Det*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k)+Al*hNew(l))
            Vxx=Vxx+CAxz
          Ai=CAxy*bi+CAyy*ci+CAyz*di
          Aj=CAxy*bj+CAyy*cj+CAyz*dj
          Ak=CAxy*bk+CAyy*ck+CAyz*dk
          Al=CAxy*bl+CAyy*cl+CAyz*dl
            Vyy=Det*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k)+Al*hNew(l))
            Vyy=Vyy+CAyz
          Ai=CAxz*bi+CAyz*ci+CAzz*di
          Aj=CAxz*bj+CAyz*cj+CAzz*dj
          Ak=CAxz*bk+CAyz*ck+CAzz*dk
          Al=CAxz*bl+CAyz*cl+CAzz*dl
            Vzz=Det*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k)+Al*hNew(l))
            Vzz=Vzz+CAzz
          DO 12 m=1,4
             l=List(m)
             Vx(l)=Vx(l)-Con(l)*Vxx
             Vy(l)=Vy(l)-Con(l)*Vyy
             Vz(l)=Vz(l)-Con(l)*Vzz
12          CONTINUE
13       CONTINUE
14    CONTINUE
      DO 15 i=1,nPt
       Vx(i)=Vx(i)/ListNE(i)
       Vy(i)=Vy(i)/ListNE(i)
       Vz(i)=Vz(i)/ListNE(i)
15    CONTINUE
      RETURN
      END SUBROUTINE Veloc
!****************************************************************
      SUBROUTINE Disper(MatNum,hNew,theta,Dxy,Exy)
      USE GridData
      USE SolData
      USE WatFun
      IMPLICIT NONE
      REAL(sp),intent(in) :: hNew(maxnod),theta(maxnod)
      INTEGER(sp),intent(in) :: MatNum(maxnod)
      REAL(sp) :: Tau,Vabs,Dxy(maxnod),Exy(maxnod),ThSati
      INTEGER(sp) :: i,M
      DO 11 i=1,nPt
       M=MatNum(i)
       ThSati=(par(3,M)-par(2,M))*Dxy(i)+par(2,M)*Exy(i)
       Tau=theta(i)**(7._sp/3._sp)/ThSati**2
       Vabs=sqrt(Vx(i)*Vx(i)+Vy(i)*Vy(i)+Vz(i)*Vz(i))
       IF (Vabs.GT.0._sp) THEN
          Dispxx(i)=ChPar(3,M)*Vx(i)*Vx(i)/Vabs+ChPar(4,M)*(Vz(i)*Vz(i)+Vy(i)*Vy(i))/Vabs+Fth(hNew(i),par(:,M))*ChPar(2,M)*Tau
          Dispyy(i)=ChPar(3,M)*Vy(i)*Vy(i)/Vabs+ChPar(4,M)*(Vx(i)*Vx(i)+Vz(i)*Vz(i))/Vabs+Fth(hNew(i),par(:,M))*ChPar(2,M)*Tau
          Dispzz(i)=ChPar(3,M)*Vz(i)*Vz(i)/Vabs+ChPar(4,M)*(Vx(i)*Vx(i)+Vy(i)*Vy(i))/Vabs+Fth(hNew(i),par(:,M))*ChPar(2,M)*Tau
          Dispxy(i)=(ChPar(3,M)-ChPar(4,M))*Vx(i)*Vy(i)/Vabs
          Dispxz(i)=(ChPar(3,M)-ChPar(4,M))*Vx(i)*Vz(i)/Vabs
          Dispyz(i)=(ChPar(3,M)-ChPar(4,M))*Vy(i)*Vz(i)/Vabs
       ELSE
          Dispxx(i)=Fth(hNew(i),par(:,M))*ChPar(2,M)*Tau
          Dispyy(i)=Fth(hNew(i),par(:,M))*ChPar(2,M)*Tau
          Dispzz(i)=Fth(hNew(i),par(:,M))*ChPar(2,M)*Tau
          Dispxy(i)=0.0_dp
          Dispxz(i)=0.0_dp
          Dispyz(i)=0.0_dp
       ENDIF
11    CONTINUE
      RETURN
      END SUBROUTINE Disper
!****************************************************************
      SUBROUTINE SolveT
      USE GridData
      USE MatData, only: As, Bs
      IMPLICIT NONE
      REAL(dp) :: P,C,Sum
      INTEGER(sp):: k,N1,i,kk,kc,ii,j,L,jj,M
      N1=Npt-1
      DO 12 k=1,N1
         P=1._dp/As(NBand,k)
         kk=k+1
         kc=NBand
         DO 11 i=kk,Npt
            kc=kc-1
            IF (kc.LT.0) goto 12
            C=-P*As(kc,i)
            As(kc,i)=C
            ii=kc+1
            L=kc+NBand-1
            DO 11 j=ii,L
             jj=j+NBand-kc
             As(j,i)=As(j,i)+C*As(jj,k)
11       CONTINUE
12    CONTINUE
      DO 14 i=2,Npt
       jj=NBand+1-i
       ii=1
       IF (jj.LT.0) THEN
          jj=1
          ii=i-NBand+1
       ENDIF
       Sum=0.0_dp
       DO 13 j=jj,NBand-1
          Sum=Sum+As(j,i)*Bs(ii)
          ii=ii+1
13     CONTINUE
       Bs(i)=Bs(i)+Sum
14    CONTINUE
      Bs(Npt)=Bs(Npt)/As(NBand,Npt)
      DO 16 k=1,N1
       i=Npt-k
       jj=i
       m=min(2*NBand-1,NBand+k)
       Sum=0.0_dp
       DO 15 j=NBand+1,m
          jj=jj+1
          Sum=Sum+As(j,i)*Bs(jj)
15     CONTINUE
       Bs(i)=(Bs(i)-Sum)/As(NBand,i)
16    CONTINUE
      RETURN
      END SUBROUTINE SolveT
!****************************************************************
      SUBROUTINE PeCour(MatNum,dtMaxC,dt,hNew,theta)
      USE GridData
      USE SolData
      USE WatFun
      IMPLICIT NONE
      INTEGER(sp), intent(in):: MatNum(maxnod)
      REAL(sp), intent(in) :: dt,hNew(maxnod),theta(maxnod)
      REAL(sp), intent(out):: dtMaxC
      INTEGER(sp)::iL(4),i,j,k,l,ie,ise
      REAL(sp):: xmax,xmin,zmin,zmax,ymin,ymax,pecx,pecy,delx
      REAL(sp):: r1,r2,r3,r4,rmin
      REAL(sp):: vzmax,vymax,vxmax,vze,vye,vxe,dze,dye,dxe,delz,dely
      REAL(sp):: courx,coury,courz,cour1,cour2,cour3,dt3,dt2,dt1,pecz
      Peclet=0.0_sp
      Courant=0.0_sp
      dtMaxC=1.e+30_sp
      DO 12 iE=1,Nelm
       DO 11 iSE=1,3
          PecX=99999._sp
          PecY=99999._sp
          PecZ=99999._sp
            dt1=1.e+30_sp
            dt2=1.e+30_sp
            dt3=1.e+30_sp
          iL(1)=3+iSE
          iL(2)=4+iSE-iSE/3*6
          iL(3)=5+iSE-iSE/2*6
          iL(4)=iSE
            i=elmnod(iE,iL(1))
            j=elmnod(iE,iL(2))
            k=elmnod(iE,iL(3))
            l=elmnod(iE,iL(4))
!            thetai=Fth(hNew(i),par(:,MatNum(i)))
!            thetaj=Fth(hNew(j),par(:,MatNum(j)))
 !           thetak=Fth(hNew(k),par(:,MatNum(k)))
 !           thetal=Fth(hNew(l),par(:,MatNum(l)))
!JAVAUX AMIN1 and AMAX! have been changed to min and max
          xmax=max(xGrid(i),xGrid(j),xGrid(k),xGrid(l))
          xmin=min(xGrid(i),xGrid(j),xGrid(k),xGrid(l))
          ymax=max(yGrid(i),yGrid(j),yGrid(k),yGrid(l))
          ymin=min(yGrid(i),yGrid(j),yGrid(k),yGrid(l))
          zmax=max(zGrid(i),zGrid(j),zGrid(k),zGrid(l))
          zmin=min(zGrid(i),zGrid(j),zGrid(k),zGrid(l))
            delX=xmax-xmin
            delY=ymax-ymin
            delZ=zmax-zmin
          DxE=(Dispxx(i)+Dispxx(j)+Dispxx(k)+Dispxx(l))/4
          DyE=(Dispyy(i)+Dispyy(j)+Dispyy(k)+Dispyy(l))/4
          DzE=(Dispzz(i)+Dispzz(j)+Dispzz(k)+Dispzz(l))/4
          VxE=abs(Vx(i)+Vx(j)+Vx(k)+Vx(l))/4
          VyE=abs(Vy(i)+Vy(j)+Vy(k)+Vy(l))/4
          VzE=abs(Vz(i)+Vz(j)+Vz(k)+Vz(l))/4
          IF (DxE.GT.0._sp) PecX=VxE*delX/DxE
          IF (DyE.GT.0._sp) PecY=VyE*delY/DyE
          IF (DzE.GT.0._sp) PecZ=VzE*delZ/DzE
!JAVAUX AMIN1 and AMAX! have been changed to min and max
          IF (PecX.NE.99999._sp) Peclet=max(Peclet,PecX)
          IF (PecY.NE.99999._sp) Peclet=max(Peclet,PecY)
          IF (PecZ.NE.99999._sp) Peclet=max(Peclet,PecZ)
!JAVAUX AMIN1 and AMAX! have been changed to min and max
          Peclet=min(Peclet,99999._sp)
          VxMax=max(abs(Vx(i))/theta(i),abs(Vx(j))/theta(j),abs(Vx(k))/theta(k),abs(Vx(l))/theta(l))
          VyMax=max(abs(Vy(i))/theta(i),abs(Vy(j))/theta(j),abs(Vy(k))/theta(k),abs(Vy(l))/theta(l))
          VzMax=max(abs(Vz(i))/theta(i),abs(Vz(j))/theta(j),abs(Vz(k))/theta(k),abs(Vz(l))/theta(l))
          R1=1._sp+ChPar(1,MatNum(i))*ChPar(5,MatNum(i))/theta(i)
          R2=1._sp+ChPar(1,MatNum(j))*ChPar(5,MatNum(j))/theta(j)
          R3=1._sp+ChPar(1,MatNum(k))*ChPar(5,MatNum(k))/theta(k)
          R4=1._sp+ChPar(1,MatNum(l))*ChPar(5,MatNum(l))/theta(l)
          RMin=min(R1,R2,R3,R4)
          CourX=VxMax*dt/(delX*RMin)
          CourY=VyMax*dt/(delY*RMin)
          CourZ=VzMax*dt/(delZ*RMin)
!JAVAUX AMIN1 and AMAX! have been changed to min and max
          Courant=max(Courant,CourX,CourY,CourZ)
          Cour1=1.0_sp
          Cour2=1.0_sp
          Cour3=1.0_sp
!JAVAUX AMIN1 and AMAX! have been changed to min and max
          IF(PecX.ne.99999._sp) Cour1=min(1._sp,PeCr/max(0.5_sp,PecX))
          IF(PecY.ne.99999._sp) Cour2=min(1._sp,PeCr/max(0.5_sp,PecY))
          IF(PecZ.ne.99999._sp) Cour3=min(1._sp,PeCr/max(0.5_sp,PecZ))
          IF(VxMax.gt.0._sp) dt1=Cour1*delX*RMin/VxMax
          IF(VyMax.gt.0._sp) dt2=Cour2*delY*RMin/VyMax
          IF(VzMax.gt.0._sp) dt3=Cour3*delZ*RMin/VzMax
!JAVAUX AMIN1 and AMAX! have been changed to min and max
          dtMaxC=min(dtMaxC,dt1,dt2,dt3)
11       CONTINUE
12    CONTINUE
      RETURN
      END SUBROUTINE PeCour
!*************************************************************************************
      SUBROUTINE SolInf(dt,Kode)
      USE GridData
      USE CumData
      IMPLICIT NONE
      REAL(sp):: sMean(2),dt
      INTEGER(sp):: Kode(maxnod),i,j
      sMean(1)=0.0_sp
      sMean(2)=0.0_sp
      DO 12 i=1,nPt
         j=iabs(Kode(i))
         IF (j.NE.0) THEN
            sMean(j)=sMean(j)-Qc(i) !Qc=nodal value of solute flux
         ENDIF
12    CONTINUE
      cCumA=abs(CumCh0)+abs(CumCh1)+abs(CumChR)
      cCumT=CumCh0+CumCh1+CumChR
	DO 13 j=1,2
         ChemS(j)=ChemS(j)+sMean(j)*dt
         cCumT=cCumT+ChemS(j)
         cCumA=cCumA+abs(ChemS(j))
13    CONTINUE
      RETURN
      END SUBROUTINE SolInf
! ==============================================================================
! Source file PLANT GROWTH |||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
! current relative stress due to soil strength and solute concentration:
      SUBROUTINE Stress(sAvg,cAvg,rs,concrs,lChem)
      USE PlntData
      USE StrData
      IMPLICIT NONE
      REAL(sp) :: sAvg,cAvg,rs,concrs
      INTEGER(sp) ::  k
      LOGICAL lChem
! interpolate along piecewise linear function to get relative stress, rs,
! corresponding to actual average soil strength sAvg;
! assume rs = 0.0 at zero soil strength;
! assume rs = 1.0 at sImp (soil strength at which growth ceases completely):
      IF (sAvg.GE.sImp) THEN
         rs=1.0_dp
         ELSE IF (sAvg.GE.sc(ns)) THEN
         rs=rsc(ns)+(sAvg-sc(ns))/(sImp-sc(ns))*(1.-rsc(ns))
      ELSE
         k=ns
    1    k=k-1
         IF (k.EQ.0) THEN
! soil strength smaller than smallest ss for which rs is specified:
            rs=sAvg/sc(1)*rsc(1)
         ELSE
            IF(sc(k).GT.sAvg) GOTO 1
            rs=rsc(k)+(sAvg-sc(k))/(sc(k+1)-sc(k))*(rsc(k+1)-rsc(k))
         ENDIF
      ENDIF
! now interpolate along piecewise linear function to get relative stress, concrs,
! corresponding to actual average solute concentration cAvg;
! assume concrs = 0.0 for average solute concentration within the optimal range;
! assume concrs = 1.0 for average solute concentration outside the range of
! minimum and maximum concentration allowed (at which growth ceases completely):
      IF (lChem) THEN
         k=0
   11    k=k+1
         IF (k.GT.ncnc) THEN
            concrs=rscnc(ncnc)
         ELSE
            IF (cAvg.GT.cncp(k)) GOTO 11
            concrs=rscnc(k-1)+(cAvg-cncp(k-1))/(cncp(k)-cncp(k-1))*(rscnc(k)-rscnc(k-1))
         ENDIF
      ELSE
         concrs=0.0_dp
      END IF
      RETURN
      END SUBROUTINE Stress
!***********************************************************************************************
! current potential transpiration:
      SUBROUTINE SetTp(t,rs,concrs,LA,level,lChem)
      USE PlntData
!      USE DoussanMat
      IMPLICIT NONE
      INTEGER(sp) :: level,ifc,k
      REAL(sp) :: LA,fctpla,ftpla,tpla,t,rs,concrs
      LOGICAL lChem
      IF (level.LE.1) GOTO 100!MJ09
      IF (level.EQ.2) GOTO 200 
      IF (level.EQ.3) GOTO 300
      IF (level.EQ.4) GOTO 400 !new javaux 05- doussan
! simulation level 1 !**
! Root uptake not considered:
  100 Tpot=0.0_dp
      RETURN
! simulation level 2 !**
! ! calculate current Tpot-value from input Tpot(time)-function:
  200 ifc=0
  201 ifc=ifc+1
      IF (ifc.GT.nTpot) THEN
         Tpot=Tpotc(nTpot)
      ELSE
         IF (t.GT.tTpot(ifc)) GOTO 201
         Tpot=Tpotc(ifc-1)+(Tpotc(ifc)-Tpotc(ifc-1))*(t-tTpot(ifc-1))/(tTpot(ifc)-tTpot(ifc-1))
      ENDIF
      RETURN
! simulation level 3 !**
! calculate current TpLA-value:
  300 ifc=0
  301 ifc=ifc+1
      IF (ifc.GT.ntTpLA) THEN
         TpLA=TpLAc(ntTpLA)
      ELSE
         IF (t.GT.tTpLA(ifc)) GOTO 301
         TpLA=TpLAc(ifc-1)+(TpLAc(ifc)-TpLAc(ifc-1))*(t-tTpLA(ifc-1))/(tTpLA(ifc)-tTpLA(ifc-1))
      ENDIF
! interpolate along piecewise linear function to get transp. reduction factor,
! fTpLA, corresponding to current relative stress due to soil strength, rs;
! by definition, fTpLA = 1.0 at zero stress (rs = 0.0)
      IF (rs.GE.sfTpLA(nfTpLA)) THEN
!relative stress greater than greatest rs for which fTpLA is specified:
         fTpLA=fTpLAc(nfTpLA)
      ELSE
         k=nfTpLA
  302    k=k-1
         IF (k.EQ.0) THEN
! relative stress smaller than smallest rs for which fTpLA is specified:
            fTpLA=1.0_sp+rs/sfTpLA(1)*(fTpLAc(1)-1._sp)
         ELSE
            IF (sfTpLA(k).GT.rs) GOTO 302
            fTpLA=fTpLAc(k)+(rs-sfTpLA(k))/(sfTpLA(k+1)-sfTpLA(k))*(fTpLAc(k+1)-fTpLAc(k))
         ENDIF
      ENDIF
! interpolate along piecewise linear function to get transp. reduction factor,
! fcTpLA, corresponding to current relative stress due to solute conc., concrs;
! by definition, fcTpLA = 1.0 at zero stress (concrs = 0.0)
      IF (lChem) THEN
         IF (concrs.GE.scTpLA(ncTpLA)) THEN
! relative stress greater than greatest concrs for which fTpLA is specified:
            fcTpLA=cTpLAc(ncTpLA)
         ELSE
            k=ncTpLA
  303       k=k-1
            IF (k.EQ.0) THEN
! relative stress smaller than smallest concrs for which fTpLA is specified:
               fcTpLA=1.0_sp+concrs/scTpLA(1)*(cTpLAc(1)-1._sp)
            ELSE
               IF (scTpLA(k).GT.concrs) GOTO 303
            fcTpLA=cTpLAc(k)+(concrs-scTpLA(k))/(scTpLA(k+1)-scTpLA(k))*(cTpLAc(k+1)-cTpLAc(k))
            ENDIF
         ENDIF
      ELSE
         fcTpLA=1.0_sp
      ENDIF
      Tpot=TpLA*MIN(fTpLA,fcTpLA)*LA
      RETURN
! simulation level 4 !**
400   Tpot=9999
      RETURN
      END SUBROUTINE SetTp
!*********************************************************************************************
! integrate uptake over spatial domain:
      SUBROUTINE ActTrs
      USE GridData
      USE PlntData, only : Tact
      IMPLICIT NONE
      REAL(dp) :: sum,sume
      INTEGER(sp) :: corner(8),ie,ic
! calculate actual overall transpiration rate from soil domain:
      sum=0.0_dp
      DO 1 iE=1,nElm,2
! assign cuboid corner nodes:
       corner(1)=elmnod(iE,1)
       corner(2)=elmnod(iE+1,6)
       corner(3)=elmnod(iE,3)
       corner(4)=elmnod(iE,2)
       corner(5)=elmnod(iE,4)
       corner(6)=elmnod(iE+1,3)
       corner(7)=elmnod(iE,6)
       corner(8)=elmnod(iE,5)
! add up average cuboid sink terms, integrate over volume:
       sumE=0.0_dp
       DO 11 ic=1,8
          sumE=sumE+sink(corner(ic))
   11    CONTINUE
       sum=sum+sumE
    1 CONTINUE
      Tact=sum*dxGrid*dyGrid*dzGrid/8._dp
      RETURN
      END SUBROUTINE ActTrs
!********************************************************************************
! current water use efficiency:
      SUBROUTINE Effncy(t,rs,concrs,W,lChem)
      USE PlntData
      IMPLICIT NONE
      REAL(sp) ::  t,rs,concrs,W,fcw,fw
      INTEGER(sp) ::  ifc,k
      LOGICAL lChem
! calculate current W:
      ifc=0
    1 ifc=ifc+1
      IF (ifc.GT.ntW) THEN
       W=Wc(ntW)
      ELSE
       IF (t.GT.tW(ifc)) GOTO 1
       W=Wc(ifc-1)+(Wc(ifc)-Wc(ifc-1))*(t-tW(ifc-1))/(tW(ifc)-tW(ifc-1))
      ENDIF
! interpolate along piecewise linear function to get water use efficiency reduction
! factor, fW, corresponding to current relative stress due to soil strength, rs;
! by definition, fW = 1.0 at zero stress (rs = 0.0)
      IF (rs.GE.sfW(nsfW)) THEN
! relative stress greater than greatest rs for which fW is specified:
         fW=fwc(nsfW)
      ELSE
         k=nsfW
    2    k=k-1
         IF (k.EQ.0) THEN
! relative stress smaller than smallest rs for which fW is specified:
            fW=1.0_sp+rs/sfW(1)*(fWc(1)-1._sp)
         ELSE
            IF (sfW(k).GT.rs) GOTO 2
            fW=fWc(k)+(rs-sfW(k))/(sfW(k+1)-sfW(k)) *(fWc(k+1)-fWc(k))
         ENDIF
      ENDIF
! interpolate along piecewise linear function to get water use efficiency
! reduction factor, fcW, corresponding to current relative stress, concrs;
! by definition, fcW = 1.0 at zero stress (concrs = 0.0)
      IF (lChem) THEN
         IF (concrs.GE.scW(nscW)) THEN
! relative stress greater than greatest concrs for which fcW is specified:
            fcW=cWc(nscW)
         ELSE
            k=nscW
    3       k=k-1
            IF (k.EQ.0) THEN
! relative stress smaller than smallest rs for which fW is specified:
               fcW=1.0_sp+concrs/scW(1)*(cWc(1)-1._sp)
            ELSE
               IF (scW(k).GT.concrs) GOTO 3
               fcW=cWc(k)+(concrs-scW(k))/(scW(k+1)-scW(k))*(cWc(k+1)-cWc(k))
            ENDIF
         ENDIF
      ELSE
         fcW=1.0_dp
      ENDIF
      W=W*MIN(fW,fcW)
      RETURN
      END SUBROUTINE Effncy
!********************************************************************************
! current root/shoot partitioning ratio:
      SUBROUTINE Ratio(t,rs,concrs,RSR,lChem)
      USE PlntData
      IMPLICIT NONE
      REAL (sp) :: t,rs,concrs,RSR,fcrsr,frsr
      INTEGER(sp) ::  ifc,k
      LOGICAL lChem
! calculate current RSR-values for each soil strength:
      ifc=0
    1 ifc=ifc+1
      IF (ifc.GT.ntRSR) THEN
         RSR=RSRc(ntRSR)
      ELSE
         IF (t.GT.tRSR(ifc)) GOTO 1
         RSR=RSRc(ifc-1)+(RSRc(ifc)-RSRc(ifc-1))*(t-tRSR(ifc-1))/(tRSR(ifc)-tRSR(ifc-1))
      ENDIF
! interpolate along piecewise linear function to get root/shoot ratio reduction
! factor, fRSR, corresponding to current relative stress due to soil strength, rs;
!by definition, fRSR = 1.0 at zero stress (rs = 0.0)
      IF (rs.GE.sfRSR(nsfRSR)) THEN
! relative stress greater than greatest rs for which fRSR is specified:
         fRSR=fRSRc(nsfRSR)
      ELSE
         k=nsfRSR
    2    k=k-1
         IF (k.EQ.0) THEN
! relative stress smaller than smallest rs for which fRSR is specified:
            fRSR=1._sp+rs/sfRSR(1)*(fRSRc(1)-1._sp)
         ELSE
            IF (sfRSR(k).GT.rs) GOTO 2
            fRSR=fRSRc(k)+(rs-sfRSR(k))/(sfRSR(k+1)-sfRSR(k))*(fRSRc(k+1)-fRSRc(k))
         ENDIF
      ENDIF
!interpolate along piecewise linear function to get root/shoot ratio reduction factor,
! fcRSR, corresponding to current relative stress due to solute concentration, concrs;
! by definition, fcRSR = 1.0 at zero stress (concrs = 0.0)
      IF (lChem) THEN
         IF (concrs.GE.scRSR(nscRSR)) THEN
! relative stress greater than greatest concrs for which fcRSR is specified:
            fcRSR=cRSRc(nscRSR)
         ELSE
            k=nscRSR
    3       k=k-1
            IF (k.EQ.0) THEN
! relative stress smaller than smallest concrs for which fcRSR is specified:
               fcRSR=1.0_sp+concrs/scRSR(1)*(cRSRc(1)-1._sp)
            ELSE
               IF (scRSR(k).GT.concrs) GOTO 3
               fcRSR=cRSRc(k)+(concrs-scRSR(k))/(scRSR(k+1)-scRSR(k))*(cRSRc(k+1)-cRSRc(k))
            ENDIF
         ENDIF
      ELSE
         fcRSR=1.0_dp
      ENDIF
      RSR=RSR*MAX(fRSR,fcRSR)
      RETURN
      END SUBROUTINE Ratio
!*****************************************************************
! current leaf area increase per new dry shoot mass:
      SUBROUTINE Leaves(t,LAmshv)
      USE PlntData
      IMPLICIT NONE
      REAL(sp) ::  LAmshv,LA,t
      INTEGER(sp) ::  nt,ifc
!calculate current LA/msh:
      ifc=0
    1 ifc=ifc+1
      IF (ifc.GT.ntLA) THEN
         LAmshv=LAc(ntLA)
      ELSE
         IF (t.GT.tLA(ifc)) GOTO 1
         LAmshv=LAc(ifc-1)+(LAc(ifc)-LAc(ifc-1))*(t-tLA(ifc-1))/(tLA(ifc)-tLA(ifc-1))
      ENDIF
      RETURN
      END SUBROUTINE Leaves
!===============================================================================
! Source file ROOT GROWTH with RootTyp |||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
   SUBROUTINE RunRootTip(naxes,t)
! run RootTip for a given time (in days) from 0
! return naxes:= primary roots in RooTyp
! call C-functions, which return arrays of variables
! most of these are then saved under RootData module
   USE Typedef
   USE RootData
   IMPLICIT NONE
   Real(dp) :: origin(3)=(/ 0.0,0.0,0.0 /)
   Integer(sp):: n,simtime, n_nodes, n_meris,n_br_crea,n_br_del, n_axes,t
   Integer(sp) ::maxlist,naxes,n_nodecorr,i,nn,ipl=1!runroottip currently doesn't work with more than one plant? It should be updated? (Couvreur dec 2009)
   Integer, allocatable :: node(:)
   Integer, allocatable :: prev_node(:)
   Integer, allocatable :: ordn(:)
   Integer, allocatable :: orda(:)
   Integer, allocatable :: agerec(:)
   Integer, allocatable :: axen(:)
   Integer, allocatable :: axeF(:)
   Integer, allocatable :: meris_seg(:)
   Integer, allocatable :: prev_nodea(:)
   Integer, allocatable :: axenf(:)
   Integer, allocatable :: axeC(:)
   Integer, allocatable :: axe_c2f(:)
   Real(dp), allocatable :: xrec(:)
   Real(dp), allocatable :: yrec(:)
   Real(dp), allocatable :: zrec(:)
   Real(dp), allocatable :: xa(:)
   Real(dp), allocatable :: ya(:)
   Real(dp), allocatable :: za(:)
   Real(dp), allocatable :: diamrec(:)
   Real(dp), allocatable :: diama(:)
   Real(dp), allocatable :: agea(:)
   Real(dp), allocatable :: wc(:)
   REAL(sp) :: treal

!initialize root function
   CALL init_RootType1(origin)
   print *,'RootTyp is running ...'
    Do simtime=1,t,1
      Call iterate_roottype1(simtime)
    END DO
   print *,'RootTyp converged'
   CALL number_of_nodes(n_nodes, n_br_crea,n_br_del) !C-function
   n_meris=n_br_crea-n_br_del
   !n_meris= all meristems
   !n_nodes= number maximum of nodes: some of them have been deleted...
   !print *,n_br_crea,n_br_del,n_meris,n_nodes
   !print *,'------------------------------'
   Allocate (node(n_nodes))
   Allocate (prev_node(n_nodes))
   Allocate (ordn(n_nodes))
   Allocate (xrec(n_nodes))
   Allocate (yrec(n_nodes))
   Allocate (zrec(n_nodes))
   Allocate (diamrec(n_nodes))
   Allocate (agerec(n_nodes)) !in days
   Allocate (axen(n_nodes))
   Allocate (axenf(n_nodes))
   Allocate (axeF(n_meris))
   Allocate (orda(n_meris))
   Allocate (xa(n_meris))
   Allocate (ya(n_meris))
   Allocate (za(n_meris))
   Allocate (diama(n_meris))
   Allocate (agea(n_meris))
   Allocate (meris_seg(n_meris))
   Allocate (prev_nodea(n_meris))
   Allocate (axe_c2f(n_meris))
   Allocate (axeC(n_meris))

! extract nodes C-function < interface
 CALL extract_nodes(node,prev_node,ordn,axen,xrec,yrec,zrec,diamrec,agerec,axeF,orda,xa,&
 & ya,za,diama,agea,meris_seg,prev_nodea,n_axes,axenf,axe_c2f,axeC,n_nodecorr) 
   nrec=n_nodecorr

!meaning of node is not clear?
   xs(1:nrec)=xrec
   ys(1:nrec)=yrec
   zs(1:nrec)=-zrec!roottip in mm
   !Z must decrease downward!!
   do n=1,nrec
     if (zs(n)>0.0_sp)	 zs(n)=-zs(n) !in case there is roots >soil surface
!	 print *,axenf(n),axenf(n-1)
	 if ((n.NE.1).AND.(axenf(n).NE.(axenf(n-1)))) then
	   if (axenf(n).NE.(axenf(n-1)+1)) then
	      print *,'axe',axenf(n-1)+1,'is missing'
	    endif
	 endif
   enddo

! create R-SWMS variables for root
   irecpr(1:nrec)=prev_node
   ordseg(1:nrec)=ordn+1. !in Roottip from 0 to 7, here from 1 to 8
   ibrseg(1:nrec)=axenf
   timorg(1:nrec)=agerec !day of creation+time real=age
   segdiam(1:nrec)=diamrec

!growing apices= "meristems with axes" in RootTyp
   ngrow=n_axes
   nbr=n_axes !number of apex=number of branhces
   xg(1:n_axes)=xa
   yg(1:n_axes)=ya
   zg(1:n_axes)=-za
   irecsg(1:n_axes)=prev_nodea
   ordgrw(1:n_axes)=orda+1
   ibrgrw(1:n_axes)=axeF

!check root nodes outside of the soil
   do n=1,ngrow
     if (zg(n)>0.0_sp) zg(n)=-zg(n) !in case there is roots >soil surface
	 if ((n.NE.1).AND.(axenf(n).NE.(axenf(n-1)))) then
	   if (axenf(n).NE.(axenf(n-1)+1)) then
	      print *,'axe',axenf(n-1)+1,'is missing'
	    endif
	 endif
   enddo

!estimate length and surface of segments and get naxes
    CALL Estimate_Seglen(naxes)

!simplifiy system
     DO n=1,10 !n iterations!
        !CALL OutRoo(real(nrec),naxes,0,0,0,0,0)
        CALL SimplifyRoot
        CALL Estimate_Seglen2
     ENDDO
     CALL OutRoo(real(999),naxes,0,0,0,0,0)
     CALL AdaptOutRoot
     CALL Estimate_Seglen2
     CALL close_c_interface() !C-function
     CALL finish_RootType1() !C-function
     treal=t
     print *,'number of nodes after',n,' iterations at time t= ', treal,' is ',nrec
     CALL OutRoo(treal,naxes,0,0,0,0,0)
     CALL CheckSize(ipl)
    END SUBROUTINE RunRootTip
!********************************************************************************
   SUBROUTINE Estimate_Seglen(naxes)
   USE typedef
   USE ParamData
   USE RootData, only : seglen,segsur,brlgth,segdiam,irecpr,nrec,nbr,ordgrw&
   &,ibrseg,zs,ys,xs,irecsg,num_seg,br_rec,connex,zg,yg,xg
   IMPLICIT NONE
   REAL(sp) :: xend,yend,zend,xtop,ytop,ztop
   INTEGER(sp) :: inode,n,brn,inode2,naxes,num
   connex(:)=.FALSE.
   inode=1
   naxes=0
   DO n=1,nbr !for all branches
      inode=irecsg(n) !inode =prec
      xend=xs(inode)!end note=apex
      yend=ys(inode)
      zend=zs(inode)
	  brn=ibrseg(inode)!branch number
      brlgth(brn)=0.0_dp
	  num=0
	  DO WHILE (brn==n)
		 inode2=irecpr(inode)
		 xtop=xs(inode2)
		 ytop=ys(inode2)
         ztop=zs(inode2)
		 seglen(inode)=SQRT((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
		 segsur(inode)=seglen(inode)*pi*segdiam(inode)!segment surface
		 brlgth(brn)=brlgth(brn)+seglen(inode)
		 num=num+1;
		 br_rec(brn)=inode2;!record at which this axe is connected
!next step:branching of previous node
		 brn=ibrseg(inode2)
		 xend=xtop
		 yend=ytop
		 zend=ztop
		 inode=inode2
      ENDDO
	  num_seg(n)=num
	 ! print *, 'numsegtest',num_seg(n)
	  IF (ordgrw(n)==1) naxes=naxes+1 !number of axes
	  connex(br_rec(n))=.TRUE. !logical which defines whtehr there is a connection to that node
   ENDDO
!   if ((inode-1).NE.(nrec)) print *,'problem 2',inode,nrec
   END SUBROUTINE Estimate_Seglen
!****************************************************************************
   SUBROUTINE SimplifyRoot
   USE typedef
   USE RootData, only : irecpr,nrec,ibrseg,irecsg,connex,xs,ys,zs,timorg,segdiam,ordseg,seglen,xg,yg,zg
   IMPLICIT NONE
   real(sp) :: xrec(1:nrec),yrec(1:nrec),zrec(1:nrec),agerec(1:nrec)
   real(dp) :: diamrec(1:nrec),ltot
   INTEGER(sp) :: ordn(1:nrec),axenf(1:nrec),n_nodecorr,nt
   INTEGER(sp) :: ibr,inew,iold,i_connex,prec(nrec),i2prec_old,iprec_old,n,oldi(nrec)
   INTEGER(sp) ::old2new(nrec),i,newtot,prec2(nrec),ibrprec,iprec(1:nrec)
   LOGICAL :: delnode

!initialisation
    xrec=xs(1:nrec)
    yrec=ys(1:nrec)
    zrec=zs(1:nrec)
    axenf=ibrseg(1:nrec)
    agerec=timorg(1:nrec)
    diamrec=segdiam(1:nrec)
    ordn=ordseg(1:nrec)
    iprec=irecpr(1:nrec)
    n_nodecorr=nrec

   inew=0
   iold=0!start with node 1==apex

!check each root
   DO WHILE (iold<n_nodecorr)
      inew=inew+1
      nt=inew !total number of nodes
      iold=iold+1
!if previous apical node than delete it
      delnode=.FALSE.
      if ((xrec(iold)==xg(axenf(iold))).and.(yrec(iold)==yg(axenf(iold))).and.&
	  &(zrec(iold)==zg(axenf(iold)))) then
         irecsg(axenf(iold))=inew
!if previous node is  on the same branch
         if (axenf(iold)==axenf(iold+1)) then
           old2new(iold)=inew
           iold=iold+1
          endif	 !else: just one apex connected to the next node node letion
         delnode=.TRUE.
      endif

	  ibrseg(inew)=axenf(iold) !br#
	  xs(inew)=xrec(iold)
	  ys(inew)=yrec(iold)
	  zs(inew)=zrec(iold)
	  ordseg(inew)=ordn(iold)
	  timorg(inew)=agerec(iold) !orig time
	  segdiam(inew)=diamrec(iold) !diam
	  oldi(inew)=iold !keep in mind the previous numerotation
	  old2new(iold)=inew
	  ltot=seglen(iold) !length
	  iprec_old=iprec(iold)!node before
	  i2prec_old=iprec(iprec_old) !2 nodes before
!	  print *,'i',inew,'iold',iold,'iprecold',iprec_old,connex(iprec_old)
     DO WHILE((axenf(iprec_old)==axenf(iold)).AND.(connex(iprec_old).EQV..FALSE.)&
	     &.AND.(iprec_old.NE.0).AND.(i2prec_old.NE.0).AND.(axenf(i2prec_old)==axenf(iold)).AND.(ltot<0.6)&
		 &.AND.((iprec_old).NE.irecsg(ibrseg(inew))).AND..NOT.(delnode)) !this is not the apical node
     !same branhces+ no connection to iprecold+ iprecold is not the seed and not 2 nodes before the seed
	     !we can delete iprecold
!		 print *,'deletednode',iprec_old,'iold',iold,'con',connex(iprec_old),i2prec_old,ltot
		 old2new(iprec_old)=9999
		 iold=iold+1
	     iprec_old=iprec(iold)!node before
	     i2prec_old=iprec(iprec_old)
         ltot=ltot+seglen(iold)
         delnode=.FALSE.
     ENDDO
     prec(inew)=iprec_old
     IF (irecsg(ibrseg(inew))==oldi(inew)) THEN
! if the previous node of the apex of the branch where node inew is is inew, then itmust be also updated
        irecsg(ibrseg(inew))=inew
      ENDIF
   ENDDO
   nrec=nt

!correct prec matrix
 Do i=1,nrec
    if (prec(i).NE.0_sp) THEN
       irecpr(i)=old2new(prec(i))
	   if (irecpr(i)==9999) then
  !         	   print *,'i',i,'prec old',prec(i),'prec new',old2new(prec(i))
!	   read(*,*)
	   endif
	else
       irecpr(i)=0
	endif
enddo
!write(*,'(5I5,3E12.2)') (n,irecpr(n),prec2(n),oldi(n),ibrseg(n),xs(n),ys(n),zs(n),n=1,nrec)
  END SUBROUTINE SimplifyRoot

!****************************************************************************
   SUBROUTINE AdaptOutRoot
!adapt root description to RSWMS
   USE TypeDef
   USE RootData, only : irecpr,nrec,ibrseg,irecsg,xs,ys,zs,timorg,segdiam,ordseg,num_seg,nbr
   IMPLICIT NONE
   real(sp):: xrec(1:nrec),yrec(1:nrec),zrec(1:nrec),agerec(1:nrec)
   real(dp) :: diamrec(1:nrec)
   INTEGER(sp):: prev(1:nrec),ordn(1:nrec),axenf(1:nrec),prec(1:nrec)
   INTEGER(sp) :: ibr2,ibr1,nn,i,n, old2new(nrec)

!initialisation
    xrec=xs(1:nrec)
    yrec=ys(1:nrec)
    zrec=zs(1:nrec)
    axenf=ibrseg(1:nrec)
    agerec=timorg(1:nrec)
    prev=irecpr(1:nrec)
    diamrec=segdiam(1:nrec)
    ibr1=1
    DO n=1,nbr !for all branches
	!get the number of nodes for that branch
      nn=num_seg(n)
	  ibr1=ibr1 !ibr1=ibr2+1 at the run  of teh next loop!!
	  ibr2=ibr1+nn-1
!adapt previous node to the meristem
	  irecsg(n)=ibr2
	  DO i=1,nn
	    xs(ibr2)=xrec(ibr1)
	    ys(ibr2)=yrec(ibr1)
	    zs(ibr2)=zrec(ibr1)
	    ibrseg(ibr2)=axenf(ibr1)
	    timorg(ibr2)=agerec(ibr1) !orig time
	    segdiam(ibr2)=diamrec(ibr1) !diam
           prec(ibr2)=prev(ibr1)
           old2new(ibr1)=ibr2
           ibr1=ibr1+1
           ibr2=ibr2-1
	  ENDDO
	ENDDO

!correct prec matrix
 Do i=1,nrec
    if (prec(i).NE.0_sp) THEN
       irecpr(i)=old2new(prec(i))
	else
       irecpr(i)=0
	endif
enddo
   END SUBROUTINE AdaptOutRoot
!*********************************************************************
   SUBROUTINE Estimate_Seglen2
   USE typedef
   USE ParamData
   USE RootData, only : seglen,segsur,brlgth,segdiam,irecpr,nrec,nbr,ordgrw&
   &,ibrseg,zs,ys,xs,irecsg,num_seg,br_rec,connex,zg,yg,xg
   IMPLICIT NONE
   REAL(sp) :: xend,yend,zend,xtop,ytop,ztop
   INTEGER(sp) :: inode,n,brn,inode2,naxes,num,connected2
   connex(:)=.FALSE.
   inode=1
   naxes=0
   DO n=1,nbr !for all branches
      inode=irecsg(n) !inode =apex<-meristem
      xend=xg(n)!end note=apex
      yend=yg(n)
      zend=zg(n)
	  brn=n !branch number
      brlgth(brn)=0.0_dp
	  num=0
!		 		 print *,'seglentest1',seglen(inode),inode,n,brn,xend
	  DO WHILE (brn==n)
!	  	 inode2=irecpr(inode)
		 xtop=xs(inode)
		 ytop=ys(inode)
         ztop=zs(inode)
		 seglen(inode)=SQRT((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2) !inode222
		 		 segsur(inode)=seglen(inode)*pi*segdiam(inode)!segment surface
	! length of node  is the length of the segment located on the top of it
	!length of node 1=0
	     if (seglen(inode)==0) then	!apical node directly connected to another banch
		      xtop=xs(irecpr(inode))
			  ytop=ys(irecpr(inode))
			  ztop=zs(irecpr(inode))
	         seglen(inode)=SQRT((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2) !inode222
			 segsur(inode)=seglen(inode)*pi*segdiam(inode)!segment surface
	      endif
		 brlgth(brn)=brlgth(brn)+seglen(inode)
         num=num+1
!		 print *,'seglentest1',brn,inode,irecpr(inode),xtop,xend,segsur(inode),seglen(inode)
		 xend=xtop
		 yend=ytop
		 zend=ztop
         inode=irecpr(inode)
		 		 connected2=inode
		 brn=ibrseg(inode)
      ENDDO
      num_seg(n)=num
	  connex(connected2)=.TRUE.
!	  print *,'final',num_seg(n)
    ENDDO
!   if ((inode-1).NE.(nrec)) print *,'problem 2',inode,nrec
   END SUBROUTINE Estimate_Seglen2
!****************************************************************************
   SUBROUTINE CheckSize(ipl)		!Checksize has to know wich plant is being checked in order to place it correctly with report to the grid (Couvreur dec 2009)
! check max/min position of roots as compared to the soil grid
   USE RootData
   USE GridData
   IMPLICIT NONE
   REAL(dp)::maxX,maxY,maxZ,maxXs,maxYs,maxZs,minX,minY,minZ,minXs,minYs,minZs
   INTEGER(sp)::ipl
   maxX=minval(xGrid)+nex*dxgrid/2		!Adapted to continuous and non continuous soil domain (Couvreur dec 2009)
   minX=minval(xGrid)
   maxY=minval(YGrid)+ney*dygrid
   minY=minval(YGrid)
   maxZ=maxval(ZGrid)
   minZ=minval(ZGrid)
   maxXs=maxval(xs(1:nrec))+xplant(ipl)
   minXs=minval(xs(1:nrec))+xplant(ipl)
   maxYs=maxval(Ys(1:nrec))+yplant(ipl)
   minYs=minval(Ys(1:nrec))+yplant(ipl)
   maxZs=maxval(Zs(1:nrec))
   minZs=minval(Zs(1:nrec))
   if (maxXs>=maxX) THEN
      print *,'X root too large'
      goto 20
   endif
   if (maxYs>=maxY) THEN
      print *,'Y root too large'
      goto 20
    endif
   if (maxZs>maxZ) THEN !
      print *,'Upper root node (=',maxZs,') is higher than soil max. z (=',maxZ,')'
      GOTO 20
   endif
   if (minXs<=minX) THEN
      print *,'X root too small'
      goto 20
   endif
   if (minYs<=minY) THEN
      print *,'Y root too small'
      goto 20
   endif
   if (minZs<minZ) THEN
      print *,'Lower root node (=',minZs,') is deeper than soil min. z (=',minZ,')'
      goto 20
   endif
    RETURN
    20 print *,maxXs,maxYs,maxZs,minXs,minYs,minZs
	print *,'Please re-run R-SWMS/RootTyp'
	STOP
   END SUBROUTINE CheckSize
!==============================================================================
! Source file ROOT GROWTH ||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
! root system growth:
      SUBROUTINE ROOT(Conc,t,dmroot,mroot,sAvg,cAvg,grwfac,naxes,norder,kaxemg,level,temp,toxi)
      USE RootData
      USE tmctrl
      IMPLICIT NONE
      REAL(sp),intent(in) ::t
      REAL,intent(in) :: Conc(maxnod)
      REAL(sp), intent(inout) ::dmroot,mroot
      REAL, intent(out) :: sAvg,cAvg,grwfac
      INTEGER, intent(in) ::norder,level
      INTEGER, intent(inout) ::kaxemg,naxes
      LOGICAL , intent(in) :: temp,toxi
      INTEGER (sp) :: corner(8),nrecol
      REAL(sp) ::  dx,dy,dz,cloc,space,sloc,newlen,newmas,MPL,sumgrw
      INTEGER(sp) :: ngrwnw,iseg,igrow,imin
      nrecol=nrec
      ngrwnw=ngrow
      sAvg=0.0_dp
      cAvg=0.0_dp
      sumgrw=0.0_dp
!check if it is time to originate new axis(es):
    1 IF (kaxemg.LE.naxemg) THEN
       IF (t.GT.tnewax(kaxemg)) THEN
          CALL Newaxs(t,sumgrw,sAvg,kaxemg,naxes,ngrwnw,temp,toxi)
          kaxemg=kaxemg+1
          GOTO 1
       ENDIF
      ENDIF
! apply potential growth to all growing branches:
      DO 10 igrow=1,ngrow
! develop new branches at all established points that are 'ripe':
       CALL Brdevl(t,sumgrw,sAvg,temp,toxi,igrow,ngrwnw)
! find surrounding grid points:
       CALL Neighb(xg(igrow),yg(igrow),zg(igrow),corner,imin)
! calculate representative soil strength value felt by root tip:
       CALL StrLoc(corner,sLoc)
       sAvg=sAvg+sLoc
! calculate branch age:
       CALL Brnage(t,igrow,age)
! calculate new heading angle and tentative segment length:
       CALL Nwcomp(igrow,corner,newlen,dx,dy,dz,dtRoot,age,temp,toxi)
! calculate tentative segment mass:
       CALL Maslen(sLoc,ordgrw(igrow),MPL)
       newmas=newlen*MPL
! make a new, tentative record:
       CALL Mkrecd(igrow,newlen,newmas,t)
! 'grow', tentatively:
       CALL Grow(igrow,dx,dy,dz,newlen)
! add up tentative need for dry mass:
       sumgrw=sumgrw+newmas
10    CONTINUE
! increase root mass by accumulated increment 'dmroot' or by potential growth:
      IF (level.EQ.3) THEN
       mroot=mroot+MIN(dmroot,sumgrw)
      ELSE
       mroot=mroot+sumgrw
      ENDIF
! average soil strength experienced by growing root system:
      sAvg=sAvg/REAL(ngrwnw)
! calculate growth factor from potential growth and available assimilates:
      IF ((level.EQ.3).AND.(sumgrw.GT.0.0_dp)) THEN
       grwfac=dmroot/sumgrw
      ELSE
       grwfac=0.0_dp
      ENDIF
! reset 'dmroot':
      dmroot=0.0_dp
      DO 20 igrow=1,ngrwnw
       IF (level.EQ.3) THEN
! adjust growth according to available assimilates:
          CALL Adjust(igrow,grwfac,newlen)
       ELSE
          newlen=seglen(irecsg(igrow))
          IF (newlen.GT.1.E-06_dp) THEN
             toosml(irecsg(igrow))=.FALSE.
             brlgth(igrow)=brlgth(igrow)+newlen
          ELSE
             toosml(irecsg(igrow))=.TRUE.
          ENDIF
       ENDIF
       IF ((.NOT.(stopgr(igrow))).AND.(.NOT.(toosml(irecsg(igrow)))).AND.(ordgrw(igrow).LT.norder)) THEN
! calculate branch spacing along new segment:
          CALL Spacng(igrow,space)
! ...and establish new branching points:
          CALL Establ(igrow,newlen,space,t,dtRoot)
       ENDIF
   20 CONTINUE
! remove all new segments smaller than 1.E-6:
      CALL Remove(nrecol)
      DO 30 igrow=1,ngrwnw
       CALL Boxlim(igrow,dtRoot)
       IF (stopgr(igrow)) THEN
! branch has reached its maximum length, make final record:
          CALL Mkrecd(igrow,0.0_sp,0.0_sp,t+dtRoot)
       ENDIF
30    CONTINUE
! remove branches that have stopped growth from list of growing branches
! and update the number of growing tips:
      CALL Update(ngrwnw)
! calculate average solute concentration felt by the root system:
        DO 35 iseg=1,nrec
           CALL      Neighb(xs(iseg),ys(iseg),zs(iseg),corner,imin)
           cLoc=(Conc(corner(1))+Conc(corner(2))+Conc(corner(3))+Conc(corner(4))+&
                 Conc(corner(5))+Conc(corner(6))+Conc(corner(7))+Conc(corner(8)))/8
           cAvg=cAvg+cLoc
   35        CONTINUE
      cAvg=cAvg/REAL(nrec)
      RETURN
      END SUBROUTINE ROOT
!********************************************************************************
! originate new axes:
      SUBROUTINE Newaxs(t,sumgrw,sAvg,kaxemg,naxes,ngrwnw,temp,toxi)
      Use RootData
      Use ParamData
      IMPLICIT NONE
      REAL(sp) :: newlen,newmas,MPL
      INTEGER (sp) :: corner(8),i
      REAL(sp) ::  sAvg,sumgrw,t,tloc,alpha,sloc,v
      REAL(sp) ::  dx,dy,dz,dxgeo,dygeo,dzgeo,dxstr,dystr,dzstr
      INTEGER (sp) :: kaxemg,naxes,ngrwnw,imin
      LOGICAL temp,toxi
! generate a total of 'nnewax' new axes:
      DO 10 i=1,nnewax(kaxemg)
! increment total axes number:
       naxes=naxes+1
! increment total branch number and total number of growing branches:
       nbr=nbr+1
       IF (ngrwnw.LT.maxgrw) THEN
          ngrwnw=ngrwnw+1
       ELSE
          WRITE(*,'(''Maximum number of growing tips -- PROGRAM TERMINATED.'')')
          STOP
       ENDIF
! assign values for axis number, branching order and branch number:
       iaxis(ngrwnw)=naxes
       ordgrw(ngrwnw)=1
       ibrgrw(ngrwnw)=nbr
! as of now, the new branch has no estbl. brnch. pts. itself:
       nestbl(ngrwnw)=0
! also, length is zero and no pending brnch. pts. exist:
       brlgth(ngrwnw)=0.0_sp
       ovrlen(ngrwnw)=-1._sp
! no segment behind the tip:
       irecsg(ngrwnw)=0
! emergence position:
       xg(ngrwnw)=xs(1)
       yg(ngrwnw)=ys(1)
       zg(ngrwnw)=zs(1)
! find surrounding grid points:
       CALL Neighb(xg(ngrwnw),yg(ngrwnw),zg(ngrwnw),corner,imin)
! initial heading components --
! get soil strength gradient components:
       CALL Ssgcom(corner,1._sp,strsen(1),dxstr,dystr,dzstr)
! current local temperature:
       CALL TemLoc(temp,corner,tLoc)
! geotropism component (betapr = preferrential heading angle with xy-plane):
       alpha=ran(irnd)*2._sp*pi
       CALL Geocom(tLoc,1,iaxis(ngrwnw),1._sp,alpha,dxgeo,dygeo,dzgeo)
! add up all components:
       dx=dxstr+dxgeo
       dy=dystr+dygeo
       dz=dzstr+dzgeo
! calculate representative soil strength value felt by root tip:
       CALL StrLoc(corner,sLoc)
       sAvg=sAvg+sLoc
! calculate length of the new segment from time left for growth:
       CALL Uniera(0.0_sp,1,v)
       CALL Length(v,t-tnewax(kaxemg),newlen,corner,temp,toxi)
! check if maximum length has been reached:
       IF (newlen.GE.brlmax(1)) THEN
          newlen=brlmax(1)
       ENDIF
! make sure new segment is not too small (prevent divison by zero):
      newlen=MAX(newlen,1.E-4_sp)
! calculate tentative segment mass:
       CALL Maslen(sLoc,ordgrw(ngrwnw),MPL)
       newmas=newlen*MPL
! make a new, tentative record:
       CALL Mkrecd(ngrwnw,newlen,newmas,tnewax(kaxemg))
! 'grow', tentatively:
       CALL Grow(ngrwnw,dx,dy,dz,newlen)
! add up tentative need for dry mass:
       sumgrw=sumgrw+newmas
10    CONTINUE
      RETURN
      END SUBROUTINE Newaxs
!****************************************************************************
! originate new sub-branches from branch 'igrow':
      SUBROUTINE Brdevl(t,sumgrw,sAvg,temp,toxi,igrow,ngrwnw)
      USE RootData
      USE ParamData
      IMPLICIT NONE
      INTEGER(sp) :: irec,iprev,igrow,iest,ngrwnw,kest
      REAL(sp) :: newlen,newmas,MPL
      REAL(sp) :: alpha,beta,gamma,delta,deltat,sloc,v
      REAL(sp) :: sAvg,sumgrw,t,test
      REAL(sp) :: dx,dy,dz,t1,t2,x1,x2,y1,y2,z1,z2
      INTEGER(sp) :: corner(8),imin
      LOGICAL temp,toxi
      kest=0
      DO 10 iest=1,nestbl(igrow)
       test=timest(igrow,iest)
       IF ((t-test).GT.dtbrch(ordgrw(igrow))) THEN
! we have an established point 'ripe' for branching:
          kest=kest+1
! increment total branch number and total number of growing branches:
          nbr=nbr+1
          IF (ngrwnw.LT.maxgrw) THEN
             ngrwnw=ngrwnw+1
          ELSE
             WRITE(*,'(/''Maximum number of growing tips -- PROGRAM TERMINATED.'')')
             STOP
          ENDIF
! assign values for axis number, branching order and branch number:
          iaxis(ngrwnw)=iaxis(igrow)
          ordgrw(ngrwnw)=ordgrw(igrow)+1
          ibrgrw(ngrwnw)=nbr
! as of now, the new branch has no estbl. brnch. pts. itself:
          nestbl(ngrwnw)=0
! also, length is zero and no pending brnch. pts. exist:
          brlgth(ngrwnw)=0.0_sp
          ovrlen(ngrwnw)=-1._sp
! find the segment where the branching occurs:
          irec=irecsg(igrow)
          IF (timorg(irec).LE.test) THEN
! have branching in the most recent segment:
             iprev=irec
             x2=xg(igrow)
             y2=yg(igrow)
             z2=zg(igrow)
             t2=t
          ELSE
! need to go back through the segment records:
    1        IF (timorg(irecpr(irec)).GT.test) THEN
                irec=irecpr(irec)
                GOTO 1
             ENDIF
             iprev=irecpr(irec)
             x2=xs(irec)
             y2=ys(irec)
             z2=zs(irec)
             t2=timorg(irec)
          ENDIF
! this is also the segment behind the tip:
          irecsg(ngrwnw)=iprev
          x1=xs(iprev)
          y1=ys(iprev)
          z1=zs(iprev)
          t1=timorg(iprev)
          dx=x2-x1
          dy=y2-y1
          dz=z2-z1
          deltat=t2-t1
! calculate exact branching position from time of establishment:
          xg(ngrwnw)=x1+(test-t1)/deltat*dx
          yg(ngrwnw)=y1+(test-t1)/deltat*dy
          zg(ngrwnw)=z1+(test-t1)/deltat*dz
! find surrounding grid points:
          CALL Neighb(xg(ngrwnw),yg(ngrwnw),zg(ngrwnw),corner,imin)
! calculate representative soil strength value felt by root tip:
          CALL StrLoc(corner,sLoc)
          sAvg=sAvg+sLoc
!calculate initial heading components:
          gamma=ran(irnd)*2._sp*pi
          delta=brnang(ordgrw(igrow))
          CALL Angchg(dx,dy,dz,alpha,beta,gamma,delta)
! calculate length of the new segment from time left for growth:
          CALL Uniera(0.0_sp,ordgrw(ngrwnw),v)
          CALL Length(v,t-dtbrch(ordgrw(igrow))-test,newlen,corner,temp,toxi)
! check if maximum length has been reached:
          IF (newlen.GE.brlmax(ordgrw(ngrwnw))) THEN
             newlen=brlmax(ordgrw(ngrwnw))
          ENDIF
! make sure new segment is not too small (prevent divison by zero):
          newlen=MAX(newlen,1.E-04_sp)
!calculate tentative segment mass:
          CALL Maslen(sLoc,ordgrw(ngrwnw),MPL)
          newmas=newlen*MPL
! make a new, tentative record:
          CALL Mkrecd(ngrwnw,newlen,newmas,test+dtbrch(ordgrw(igrow)))
! 'grow', tentatively:
          CALL Grow(ngrwnw,dx,dy,dz,newlen)
! add up tentative need for dry mass:
          sumgrw=sumgrw+newmas
       ENDIF
10    CONTINUE
! remove all developed points from the list of established branching points:
      nestbl(igrow)=nestbl(igrow)-kest
      DO 20 iest=1,nestbl(igrow)
       timest(igrow,iest)=timest(igrow,iest+kest)
20    CONTINUE
      RETURN
      END SUBROUTINE Brdevl
!***************************************************************************
! find 8 corner points of cube (double element) that contains root tip
! or segment originating or ending point:
      SUBROUTINE Neighb(x,y,z,corner,imin)
      USE GridData
      USE DomData
      IMPLICIT NONE
       REAL(sp):: delta1,x,y,z
      INTEGER(sp):: ixmin,ixymin,imin,corner(8),noutdomain
! nElm = total number of elements
!  ne* = number of elements (half-cuboids) in '*'-direction
!  nel = nex*ney (number of elements per horizontal element-layer)
! imin = element-ID# of odd-numbered half of the cuboid containing the root
!        tip or segment
! find imin by checking distance between tip coord's and cuboid center-points
! search along x-axis (start at first element):
      delta1=x-(xmin)!-1.E-6_dp  Why not to take the real position of the node ? Substracting 1.E-6 risks to take the node out of the domain/element (Couvreur dec 2009)
      IF (delta1.GE.nex*dxgrid/2) STOP 'Root node out of the soil domain (Xdirection)'		!Verification seems not to be necessary (Couvreur dec 2009)
      ixmin=1+floor(delta1/dxgrid)*2
! search along y-axis (start at ixmin):
      delta1=y-(ymin)!-1.E-6_dp
      IF (delta1.GE.ney*dygrid) STOP 'Root node out of the soil domain (Ydirection)'		!(Couvreur dec 2009)
      ixymin=ixmin+floor(delta1/dygrid)*nex
! search along z-axis (start at ixymin):
      delta1=abs(z-zmax)!+1.E-6_dp
      IF (delta1.GT.nez*dzgrid) STOP 'Root node out of the soil domain (Zdirection)'		!(Couvreur dec 2009)
      imin=ixymin+floor(delta1/dzgrid)*nel
!print *,'cornercheck',x,y,z,xmin,ymin,zmax,imin
!assign cuboid corner nodes:
      corner(1)=elmnod(imin,1)
      corner(2)=elmnod(imin+1,6)
      corner(3)=elmnod(imin,3)
      corner(4)=elmnod(imin,2)
      corner(5)=elmnod(imin,4)
      corner(6)=elmnod(imin+1,3)
      corner(7)=elmnod(imin,6)
      corner(8)=elmnod(imin,5)
      RETURN
      END SUBROUTINE Neighb
 !--------------------------------------------------------------
 SUBROUTINE FlowNeighb(corner,Flow_corner)
      USE GridData
      IMPLICIT NONE
      integer(sp) :: i,c,corner(8)
      integer(sp), intent(out) :: Flow_corner(1:8,1:6)
      Flow_corner=0
      do c=1,8
       Flow_corner(c,1) = corner(c) + 1
       Flow_corner(c,2) = corner(c) - 1
       Flow_corner(c,3) = corner(c) + 1 + ney
       Flow_corner(c,4) = corner(c) - 1 - ney
       Flow_corner(c,5) = corner(c) + (nPt/(nez+1))
       Flow_corner(c,6) = corner(c) - (nPt/(nez+1))
       do i=1,6
          if (Flow_corner(c,i) .lt. 0) Flow_corner(c,i)=0
       enddo
       end do
      RETURN
  END SUBROUTINE FlowNeighb
!**********************************************************************
! branch age:
      SUBROUTINE Brnage(t,igrow,tempage)
      USE RootData
      IMPLICIT NONE
      INTEGER*4 igrow,irec
      REAL(sp):: t,tempage
      irec=irecsg(igrow)
! go back through the segment records to first record of branch 'igrow':
    1 IF (irecpr(irec).NE.0) THEN
       IF (ordseg(irecpr(irec)).EQ.ordseg(irec)) THEN
          irec=irecpr(irec)
          GOTO 1
       ENDIF
      ENDIF
      tempage=t-timorg(irec)
      RETURN
      END SUBROUTINE Brnage
!******************************************************************************
! unimpeded elongation rate:
      SUBROUTINE Uniera(tempage,iorder,v)
      USE Typedef
      USE RootData
      IMPLICIT NONE
      REAL(sp)::  tempage,v
      INTEGER(sp):: ivch,iorder
! calculate unimpeded elongation rate as a function of order, age:
      IF (tempage.GE.agevch(iorder,nvch(iorder))) THEN
       v=vch(iorder,nvch(iorder))
      ELSE
       ivch=nvch(iorder)
    1    ivch=ivch-1
       IF ((tempage.LT.agevch(iorder,ivch)).AND.(ivch.GT.1)) GOTO 1
       v=vch(iorder,ivch)+(tempage-agevch(iorder,ivch))/(agevch(iorder,ivch+1)-&
            agevch(iorder,ivch))*(vch(iorder,ivch+1)-vch(iorder,ivch))
      ENDIF
      RETURN
      END SUBROUTINE Uniera
!******************************************************************************
! length of new segment:
      SUBROUTINE Length(v,dtRoot,newlen,corner,temp,toxi)
      USE TempData
      USE ConData
      USE StrData
      IMPLICIT NONE
      REAL(sp)::  v,newlen,cncimp,temimp,strimp,dtroot
      INTEGER(sp)::  corner(8)
      LOGICAL temp,toxi
! calculate length of new segment, taking location into account:
      strimp=(imps(corner(1))+imps(corner(2))+imps(corner(3))+imps(corner(4))+&
              imps(corner(5))+imps(corner(6))+imps(corner(7))+imps(corner(8)))/8._sp
      IF (temp) THEN
         temimp=(impt(corner(1))+impt(corner(2))+impt(corner(3))+impt(corner(4))+&
                 impt(corner(5))+impt(corner(6))+impt(corner(7))+impt(corner(8)))/8._sp
      ELSE
         temimp=1.
      ENDIF
      IF (toxi) THEN
         cncimp=(impc(corner(1))+impc(corner(2))+impc(corner(3))+impc(corner(4))+&
                 impc(corner(5))+impc(corner(6))+impc(corner(7))+impc(corner(8)))/8._sp
      ELSE
         cncimp=1._sp
      ENDIF
      newlen=(cncimp*strimp*temimp)*v*dtRoot
      RETURN
      END SUBROUTINE Length
!******************************************************************************
! soil strength gradient components of new heading vector:
      SUBROUTINE Ssgcom(corner,stndrd,strsen,dxstr,dystr,dzstr)
      USE GridData
      USE StrData
      IMPLICIT NONE
      REAL (sp):: stndrd,strsen,dxstr,dystr,dzstr,factor
      INTEGER(sp) :: corner(8)
! calculate, normalize, and weight soil strength gradient components:
      dxstr=(s(corner(1))+s(corner(3))+s(corner(5))+s(corner(7))-&
             s(corner(2))-s(corner(4))-s(corner(6))-s(corner(8)))/(4._sp*dxgrid)
      dystr=(s(corner(1))+s(corner(2))+s(corner(5))+s(corner(6))-&
             s(corner(3))-s(corner(4))-s(corner(7))-s(corner(8)))/(4._sp*dygrid)
      dzstr=(s(corner(1))+s(corner(2))+s(corner(3))+s(corner(4))-&
             s(corner(5))-s(corner(6))-s(corner(7))-s(corner(8)))/(4._sp*dzgrid)
      IF ((ABS(dxstr).GT.1.E-10_sp).OR.(ABS(dystr).GT.1.E-10_sp).OR.(ABS(dzstr).GT.1.E-10_sp)) THEN
         factor=stndrd/refgrd*strsen
         dxstr=factor*dxstr
         dystr=factor*dystr
         dzstr=factor*dzstr
      ENDIF
      RETURN
      END SUBROUTINE Ssgcom
!******************************************************************************
! geotropism components of new heading vector:
      SUBROUTINE Geocom(tLoc,iord,iax,stndrd,alpha,dxgeo,dygeo,dzgeo)
      USE typedef
      IMPLICIT NONE
      REAL(sp)::  tLoc,stndrd,alpha,dxgeo,dygeo,dzgeo,geotrp,betapr
      INTEGER (sp):: iord,iax
! calculate geotropism components
! (betapr = preferrential heading angle with xy-plane):
      IF (iord.LE.2) THEN
         CALL Prfang(tLoc,iord,iax,geotrp,betapr)
         dxgeo=stndrd*geotrp*COS(betapr)*COS(alpha)
         dygeo=stndrd*geotrp*COS(betapr)*SIN(alpha)
         dzgeo=stndrd*geotrp*SIN(betapr)
      ELSE
         dxgeo=0.0_dp
         dygeo=0.0_dp
         dzgeo=0.0_dp
      ENDIF
      RETURN
      END SUBROUTINE Geocom
!************************************************************!*
! preferrential growth angle with horizontal plane:
      SUBROUTINE Prfang(tLoc,iord,iax,geotrp,betapr)
      USE GeoData
      IMPLICIT NONE
      REAL(sp)::  geotrp,betapr,tLoc
      INTEGER(sp)::  iord,iax,igch
! interpolate along piecewise linear function to get
! preferred growth angle as a function of axis#, temperature:
      IF (iord.EQ.1) THEN
! use axis data:
       IF (tLoc.GE.tempax(iax,nangax(iax))) THEN
! temperature greater than greatest T for which betapr is specified --
! use values for greatest specified T:
          betapr=angaxs(iax,nangax(iax))
       ELSE
          iGch=nangax(iax)
    1       iGch=iGch-1
          IF (iGch.EQ.0) THEN
! temperature smaller than smallest T for which betapr is specified --
! use values for smallest specified T:
             betapr=angaxs(iax,1)
          ELSE
             IF (tLoc.LT.tempax(iax,iGch)) GOTO 1
             betapr=angaxs(iax,iGch)+(tLoc-tempax(iax,iGch))/&
                  (tempax(iax,iGch+1)-tempax(iax,iGch))*&
                  (angaxs(iax,iGch+1)-angaxs(iax,iGch))
          ENDIF
       ENDIF
       geotrp=geoaxs(iax)
      ELSE
! use main lateral data:
       IF (tLoc.GE.templt(nanglt)) THEN
! temperature greater than greatest T for which betapr is specified --
! use values for greatest specified T:
          betapr=anglat(nanglt)
       ELSE
          iGch=nanglt
    2       iGch=iGch-1
          IF (iGch.EQ.0) THEN
! temperature smaller than smallest T for which betapr is specified --
! use values for smallest specified T:
             betapr=anglat(1)
          ELSE
             IF (tLoc.LT.templt(iGch)) GOTO 2
             betapr=anglat(iGch)+(tLoc-templt(iGch))/&
                  (templt(iGch+1)-templt(iGch))*(anglat(iGch+1)-anglat(iGch))
          ENDIF
       ENDIF
       geotrp=geolat
      ENDIF
      RETURN
      END SUBROUTINE Prfang
!******************************************************************************
! length and heading vector components of new segmen:t
      SUBROUTINE Nwcomp(igrow,corner,newlen,dx,dy,dz,dtRoot,tempage,temp,toxi)
      USE RootData
      USE ParamData
      IMPLICIT NONE
      REAL(sp)::  newlen
      REAL(sp)::  alpha,beta,dxstr,dystr,dzstr,tloc,dxgeo,dygeo,dzgeo,v
      REAL(sp)::  dx,dy,dz,dtroot,tempage,stndrd,gamma,delta
      INTEGER(sp)::  corner(8),igrow
      LOGICAL temp,toxi
! use old segment length as normalizing standard for components:
      stndrd=seglen(irecsg(igrow))
! calculate the old heading components dx,dy,dz:
      dx=xg(igrow)-xs(irecsg(igrow))
      dy=yg(igrow)-ys(irecsg(igrow))
      dz=zg(igrow)-zs(irecsg(igrow))
! change old heading angle by some random amount without changing length:
      gamma=ran(irnd)*2._sp*pi
      delta=ran(irnd)*rdmang(ordgrw(igrow))
      CALL Angchg(dx,dy,dz,alpha,beta,gamma,delta)
! get soil strength gradient components:
      CALL Ssgcom(corner,stndrd,strsen(ordgrw(igrow)),dxstr,dystr,dzstr)
! current local temperature:
      CALL TemLoc(temp,corner,tLoc)
! get geotropism components:
      CALL Geocom(tLoc,ordgrw(igrow),iaxis(igrow),stndrd,alpha,dxgeo,dygeo,dzgeo)
! add up all components:
      dx=dx+dxstr+dxgeo
      dy=dy+dystr+dygeo
      dz=dz+dzstr+dzgeo
! calculate length of new segment, taking location into account:
      CALL Uniera(tempage,ordgrw(igrow),v)
      CALL Length(v,dtRoot,newlen,corner,temp,toxi)
! check if maximum length has been reached:
      IF (brlgth(igrow)+newlen.GE.brlmax(ordgrw(igrow))) THEN
       newlen=brlmax(ordgrw(igrow))-brlgth(igrow)
      ENDIF
      RETURN
      END SUBROUTINE Nwcomp
!******************************************************************************
      SUBROUTINE Angchg(dx,dy,dz,alpha,beta,gamma,delta)
      Use paramData
      IMPLICIT NONE
      REAL(sp) :: totlen,horlen,dx,dy,dz,alpha,beta,gamma,delta,a,r,c,d
      totlen=SQRT(dx*dx+dy*dy+dz*dz)
      horlen=SQRT(dx*dx+dy*dy)
      r=SIN(delta)*totlen
      a=COS(gamma)*r
      beta=ASIN(SIGN(MIN(1._sp,ABS(dz/totlen)),dz))
      IF (horlen.GT.1.E-20_sp) THEN
       alpha=ACOS(SIGN(MIN(1._sp,ABS(dx/horlen)),dx))
       IF (dy.LT.0.0_sp) alpha=2._sp*pi-alpha
       c=SIN(beta)*a
       d=SIN(gamma)*r
       dx=COS(delta)*dx-c/horlen*dx-SIN(alpha)*d
       dy=COS(delta)*dy-c/horlen*dy+COS(alpha)*d
      ELSE
       alpha=gamma
       dx=COS(alpha)*r
       dy=SIN(alpha)*r
      ENDIF
      dz=COS(delta)*dz+COS(beta)*a
      horlen=SQRT(dx*dx+dy*dy)
      beta=ASIN(SIGN(MIN(1._sp,ABS(dz/totlen)),dz))
      IF (horlen.GT.1.E-20_sp) THEN
       alpha=ACOS(SIGN(MIN(1._sp,ABS(dx/horlen)),dx))
       IF (dy.LT.0.0_sp) alpha=2._sp*pi-alpha
      ELSE
       alpha=gamma
      ENDIF
      RETURN
      END SUBROUTINE Angchg
!*****************************************************************
! mass per length:
      SUBROUTINE Maslen(sLoc,iorder,mpl)
      USE Typedef
      USE RootData
      IMPLICIT NONE
      REAL(sp) :: MPL,Sloc
      INTEGER(sp) :: iorder,implch
!calculate mass per length as a function of order, soil strength:
      IF (sLoc.GE.sMPLch(iorder,nMPLch(iorder))) THEN
       MPL=MPLch(iorder,nMPLch(iorder))
      ELSE
       iMPLch=nMPLch(iorder)
    1    iMPLch=iMPLch-1
       IF ((sLoc.LT.sMPLch(iorder,iMPLch)).AND.(iMPLch.GT.1)) GOTO 1
       MPL=MPLch(iorder,iMPLch)+(sLoc-sMPLch(iorder,iMPLch))/&
            (sMPLch(iorder,iMPLch+1)-sMPLch(iorder,iMPLch))*&
            (MPLch(iorder,iMPLch+1)-MPLch(iorder,iMPLch))
      ENDIF
      RETURN
      END SUBROUTINE Maslen
!*****************************************************************************
! new segment record:
      SUBROUTINE Mkrecd(igrow,length,mass,t)
      USE RootData
      USE PlntData
      IMPLICIT NONE
      INTEGER(sp) :: igrow
      REAL(sp) :: length,mass,t,radius,sgmpl
      IF (nrec.LT.maxrec) THEN
       nrec=nrec+1
      ELSE
       WRITE (*,'(/''Maximum number of root segments -- PROGRAM TERMINATED.'')')
       STOP
      ENDIF
      xs(nrec)=xg(igrow)
      ys(nrec)=yg(igrow)
      zs(nrec)=zg(igrow)
      irecpr(nrec)=irecsg(igrow)
      ordseg(nrec)=ordgrw(igrow)
      ibrseg(nrec)=ibrgrw(igrow)
      seglen(nrec)=length
      segmas(nrec)=mass
      timorg(nrec)=t
      if (length.eq.0) length=1e-7
      sgmpl=mass/length !new javaux
      radius=sqrt(sgmpl/SpWgt/pi)!new javaux
      segsur(nrec)=2*pi*radius*length!new javaux
      RETURN
      END  SUBROUTINE Mkrecd
!***************************************************************
! grow length 'newlen' along the heading vector (dx,dy,dz):
      SUBROUTINE Grow(igrow,dx,dy,dz,newlen)
      USE RootData
      IMPLICIT NONE
      REAL(sp) :: factor,dx,dy,dz,newlen
      INTEGER(sp) :: igrow
! calculate the actual length of components dx,dy,dz:
      factor=newlen/SQRT(dx*dx+dy*dy+dz*dz)
      dx=factor*dx
      dy=factor*dy
      dz=factor*dz
      xg(igrow)=xg(igrow)+dx
      yg(igrow)=yg(igrow)+dy
      zg(igrow)=zg(igrow)+dz
      irecsg(igrow)=nrec
      RETURN
      END SUBROUTINE Grow
!********************************************************************************
! adjust root growth to not exceed available assimilate:
      SUBROUTINE Adjust(igrow,grwfac,newlen)
      USE RootData
      IMPLICIT NONE
      INTEGER(sp) :: irec,igrow
      REAL(sp) :: newlen,grwfac,factor
! grwfac > 1.  indicates more assimilate being sent to root than can be used
! for growth under current conditions --
! assume that extra assimilate is exudated.
! find the segment behind the tip 'igrow':
      irec=irecsg(igrow)
! apply growth factor to tentative segment length:
      newlen=seglen(irec)*MIN(grwfac,1._sp)
! check if maximum length has been reached:
      IF (brlgth(igrow)+newlen.GE.brlmax(ordgrw(igrow))) THEN
       newlen=brlmax(ordgrw(igrow))-brlgth(igrow)
       stopgr(igrow)=.TRUE.
      ELSE
       stopgr(igrow)=.FALSE.
      ENDIF
! make sure first branch segments are not too small
! after being adjusted (prevent div. by zero):
      IF (irecpr(irec).EQ.0) THEN
       newlen=MAX(newlen,1.E-04_sp)
      ELSE
       IF (ordseg(irecpr(irec)).NE.ordseg(irec)) newlen=MAX(newlen,1.E-04_sp)
      ENDIF
      IF (newlen.GT.1.E-06_sp) THEN
       toosml(irec)=.FALSE.
! adjust mass of segment:
       segmas(irec)=segmas(irec)*MIN(grwfac,1._sp)
!calculate length correction factor:
       factor=newlen/seglen(irec)
! calculate exact position of tip:
       xg(igrow)=xs(irec)+factor*(xg(igrow)-xs(irec))
       yg(igrow)=ys(irec)+factor*(yg(igrow)-ys(irec))
       zg(igrow)=zs(irec)+factor*(zg(igrow)-zs(irec))
       brlgth(igrow)=brlgth(igrow)+newlen
! adjust length of segment:
       seglen(irec)=newlen
      ELSE
! we have a candidate for removal at the end of time step
       toosml(irec)=.TRUE.
      ENDIF
      RETURN
      END SUBROUTINE Adjust
!******************************************************************************
! distance between successive sub-branch origination points:
      SUBROUTINE Spacng(igrow,space)
      USE RootData
      IMPLICIT NONE
      INTEGER(sp) :: igrow
      REAL(sp) :: space
      space=brspac(ordgrw(igrow))
      RETURN
      END SUBROUTINE Spacng
!******************************************************************************
! establish new sub-branch origination points along branch 'igrow':
      SUBROUTINE Establ(igrow,newlen,space,t,dtRoot)
      USE RootData
      IMPLICIT NONE
      REAL(sp):: newlen,t,space,dtroot,over,smin,first,deltax
      INTEGER(sp):: igrow
      over=ovrlen(igrow)
      IF (over.GT.0._sp) THEN
       smin=MIN(over,space)
       IF (smin.GT.newlen) THEN
          first=-1._sp
          over=smin-newlen
       ELSE
          first=smin
       ENDIF
      ELSE
       IF (space.GT.newlen) THEN
          first=-1._sp
          over=space-newlen
       ELSE
          first=space
       ENDIF
      ENDIF
      IF ((first.GT.0._sp).AND.(nestbl(igrow).LT.maxest)) THEN
       deltax=first
    1    nestbl(igrow)=nestbl(igrow)+1
       timest(igrow,nestbl(igrow))=t+deltax/newlen*dtRoot
       deltax=deltax+space
       IF ((deltax.LE.newlen).AND.(nestbl(igrow).LT.maxest)) GOTO 1
       over=deltax-newlen
      ENDIF
      ovrlen(igrow)=over
      RETURN
      END SUBROUTINE Establ
!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!!!!*
! root has to stay within domain limits:
      SUBROUTINE Boxlim(igrow,dtRoot)
      USE RootData
      USE DomData
      IMPLICIT NONE
      INTEGER(sp) :: irec,igrow
      REAL(sp):: length,mass,newlen,newmas,dx,dy,dz
      REAL(sp):: factor,fractn,del,delmin,dtroot
      CHARACTER wall*4
! find the segment behind the tip 'igrow':
    1 irec=irecsg(igrow)
      IF (toosml(irec)) RETURN
! check if tip is outside the domain
! and if so, which side wall is intersected first:
      delmin=1.E+37_dp
      wall='none'
      IF (xg(igrow).LT.xmin) THEN
       del=(xmin-xs(irec))/(xg(igrow)-xs(irec))*seglen(irec)
       IF (del.LT.delmin) THEN
          delmin=del
          wall='xmin'
       ENDIF
      ENDIF
      IF (xg(igrow).GT.xmax) THEN
       del=(xmax-xs(irec))/(xg(igrow)-xs(irec))*seglen(irec)
       IF (del.LT.delmin) THEN
          delmin=del
          wall='xmax'
       ENDIF
      ENDIF
      IF (yg(igrow).LT.ymin) THEN
       del=(ymin-ys(irec))/(yg(igrow)-ys(irec))*seglen(irec)
       IF (del.LT.delmin) THEN
          delmin=del
          wall='ymin'
       ENDIF
      ENDIF
      IF (yg(igrow).GT.ymax) THEN
       del=(ymax-ys(irec))/(yg(igrow)-ys(irec))*seglen(irec)
       IF (del.LT.delmin) THEN
          delmin=del
          wall='ymax'
       ENDIF
      ENDIF
      IF (zg(igrow).LT.zmin) THEN
       del=(zmin-zs(irec))/(zg(igrow)-zs(irec))*seglen(irec)
       IF (del.LT.delmin) THEN
          delmin=del
          wall='zmin'
       ENDIF
      ENDIF
      IF (zg(igrow).GT.zmax) THEN
       del=(zmax-zs(irec))/(zg(igrow)-zs(irec))*seglen(irec)
       IF (del.LT.delmin) THEN
          delmin=del
          wall='zmax'
       ENDIF
      ENDIF
      IF (wall.EQ.'none') RETURN
! if we get to here, we know we have the tip outside the domain and need
! new heading components:
      dx=xg(igrow)-xs(irec)
      dy=yg(igrow)-ys(irec)
      dz=zg(igrow)-zs(irec)
      CALL Compon(wall(1:1),irec,dx,dy,dz)
      IF (delmin.GT.1.E-06_dp) THEN
! break up segment --
       length=seglen(irec)
       mass=segmas(irec)
       fractn=delmin/length
! part with unchanged heading:
       segmas(irec)=fractn*mass
       seglen(irec)=delmin
! pull back tip to start part with new heading:
       xg(igrow)=xs(irec)+fractn*(xg(igrow)-xs(irec))
       yg(igrow)=ys(irec)+fractn*(yg(igrow)-ys(irec))
       zg(igrow)=zs(irec)+fractn*(zg(igrow)-zs(irec))
! to prevent numerical problems, pull back exactly to wall
       IF (wall.EQ.'xmin') xg(igrow)=xmin
       IF (wall.EQ.'xmax') xg(igrow)=xmax
       IF (wall.EQ.'ymin') yg(igrow)=ymin
       IF (wall.EQ.'ymax') yg(igrow)=ymax
       IF (wall.EQ.'zmin') zg(igrow)=zmin
       IF (wall.EQ.'zmax') zg(igrow)=zmax
       newlen=length-seglen(irec)
       IF (newlen.GT.1.E-06_dp) THEN
          newmas=mass-segmas(irec)
          CALL Mkrecd(igrow,newlen,newmas,timorg(irec)+fractn*dtRoot)
          toosml(nrec)=.FALSE.
          CALL Grow(igrow,dx,dy,dz,newlen)
       ENDIF
      ELSE
! don't split, change heading only --
! calculate the actual length of components dx,dy,dz:
       factor=seglen(irec)/SQRT(dx*dx+dy*dy+dz*dz)
       dx=factor*dx
       dy=factor*dy
       dz=factor*dz
! new position of tip:
       xg(igrow)=xs(irec)+dx
       yg(igrow)=ys(irec)+dy
       zg(igrow)=zs(irec)+dz
      ENDIF
      GOTO 1
      END SUBROUTINE Boxlim
!******************************************************************************
      SUBROUTINE Compon(wall,irec,dx,dy,dz)
      USE RootData
      USE DomData
      IMPLICIT NONE
      INTEGER(sp):: irec
      REAL(sp):: dx,dy,dz
      CHARACTER wall*1
      IF (wall.EQ.'x') THEN
       dx=0.0_dp
       IF ((ABS(dy).LT.1.E-20_dp).AND.(ABS(dz).LT.1.E-20_dp)) THEN
          IF (ABS(ymax-ys(irec)).GT.ABS(ymin-ys(irec))) THEN
             dy=ymax-ys(irec)
          ELSE
             dy=ymin-ys(irec)
          ENDIF
          IF (ABS(zmax-zs(irec)).GT.ABS(zmin-zs(irec))) THEN
             dz=zmax-zs(irec)
          ELSE
             dz=zmin-zs(irec)
          ENDIF
       ENDIF
      ELSEIF(wall.EQ.'y') THEN
       dy=0.0_dp
       IF ((ABS(dx).LT.1.E-20_dp).AND.(ABS(dz).LT.1.E-20_dp)) THEN
          IF (ABS(xmax-xs(irec)).GT.ABS(xmin-xs(irec))) THEN
             dx=xmax-xs(irec)
          ELSE
             dx=xmin-xs(irec)
          ENDIF
          IF (ABS(zmax-zs(irec)).GT.ABS(zmin-zs(irec))) THEN
             dz=zmax-zs(irec)
          ELSE
             dz=zmin-zs(irec)
          ENDIF
       ENDIF
      ELSEIF(wall.EQ.'z') THEN
       dz=0.0_dp
       IF ((ABS(dx).LT.1.E-20_dp).AND.(ABS(dy).LT.1.E-20_dp)) THEN
          IF (ABS(xmax-xs(irec)).GT.ABS(xmin-xs(irec))) THEN
             dx=xmax-xs(irec)
          ELSE
             dx=xmin-xs(irec)
          ENDIF
          IF (ABS(ymax-ys(irec)).GT.ABS(ymin-ys(irec))) THEN
             dy=ymax-ys(irec)
          ELSE
             dy=ymin-ys(irec)
          ENDIF
       ENDIF
      ENDIF
      RETURN
      END SUBROUTINE Compon
!************************************************************************
! remove all segments from list that are too small:
      SUBROUTINE Remove(nrecol)
      USE RootData
      IMPLICIT NONE
      INTEGER(sp) :: nrecnw,nrecol,ifrom,ito,igrow
! set nrecnw equal to the current number of segment records
      nrecnw=nrec
      ito=0
      ifrom=nrecol
   10 ifrom=ifrom+1
       IF (toosml(ifrom)) THEN
! identify and pull back tip that belongs to the record being removed:
          igrow=0
   11     igrow=igrow+1
          IF (ibrgrw(igrow).NE.ibrseg(ifrom)) GOTO 11
          xg(igrow)=xs(ifrom)
          yg(igrow)=ys(ifrom)
          zg(igrow)=zs(ifrom)
          irecsg(igrow)=irecpr(ifrom)
! decrease total number of records by '1':
          nrec=nrec-1
! if this the first candidate for removal (to be overwritten),
! mark position:
          IF (ito.EQ.0) ito=ifrom
       ELSE IF (ito.GT.0) THEN
! can move from 'ifrom' because at least one removal candidate
! has been previously identified and marked as 'ito':
          toosml(ito)=toosml(ifrom)
          xs(ito)=xs(ifrom)
          ys(ito)=ys(ifrom)
          zs(ito)=zs(ifrom)
          irecpr(ito)=irecpr(ifrom)
          ordseg(ito)=ordseg(ifrom)
          ibrseg(ito)=ibrseg(ifrom)
          seglen(ito)=seglen(ifrom)
          segmas(ito)=segmas(ifrom)
          timorg(ito)=timorg(ifrom)
! tip that belongs to moved segment needs the new segment reference, ito:
          igrow=0
   12       igrow=igrow+1
          IF (ibrgrw(igrow).NE.ibrseg(ito)) GOTO 12
          irecsg(igrow)=ito
          ito=ito+1
       ENDIF
      IF (ifrom.LT.nrecnw) GOTO 10
! now nrec is the actual number of segment records
      RETURN
      END SUBROUTINE Remove
!*********************************************************************************
      SUBROUTINE Update(ngrwnw)
      USE RootData
      IMPLICIT NONE
      INTEGER(sp):: ito,ifrom,iest,ngrwnw
! may have more growing branch tips at the end of time step:
      ngrow=ngrwnw
      ito=0
      ifrom=0
   10 ifrom=ifrom+1
       IF (stopgr(ifrom)) THEN
          ngrow=ngrow-1
! if this the first candidate for removal (to be overwritten),
! mark position:
          IF (ito.EQ.0) ito=ifrom
       ELSE IF (ito.GT.0) THEN
! can move from 'ifrom' because at least one removal candidate
! has been previously identified and marked as 'ito':
          stopgr(ito)=stopgr(ifrom)
          xg(ito)=xg(ifrom)
          yg(ito)=yg(ifrom)
          zg(ito)=zg(ifrom)
          iaxis(ito)=iaxis(ifrom)
          irecsg(ito)=irecsg(ifrom)
          ordgrw(ito)=ordgrw(ifrom)
          ibrgrw(ito)=ibrgrw(ifrom)
          brlgth(ito)=brlgth(ifrom)
          ovrlen(ito)=ovrlen(ifrom)
          nestbl(ito)=nestbl(ifrom)
          DO 2 iest=1,nestbl(ito)
             timest(ito,iest)=timest(ifrom,iest)
2           CONTINUE
          ito=ito+1
       ENDIF
      IF (ifrom.LT.ngrwnw) GOTO 10
!now ngrow is the actual number of growing tips at beginning of next step
      RETURN
      END SUBROUTINE Update
! ==============================================================================
! Source file OUTPUT |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
      SUBROUTINE ObsIni

      USE Typedef
      USE ObsData
      USE GridData, only : xgrid,ygrid,zgrid
      IMPLICIT NONE

      INTEGER(sp) :: ip,minl,i
      CHARACTER filename*13,form*3,form2*2,form1*1

	DO ip=1,nPr
		 minl=MIN(nodebyPr(ip),500)
!define file names
            WRITE (filename,'(A13)')'out/Probe.  '
            IF (ip.LT.10) THEN
              WRITE (filename(11:11),'(I1)') ip
            ELSE
              WRITE (filename(11:12),'(I2)') ip
            ENDIF
	     OPEN (UNIT=10,FILE=filename,STATUS='UNKNOWN')
		!if more than 500 nodes
		IF (minl.NE.nodebyPr(ip)) THEN
                 IF (ip.LT.10) THEN
                     WRITE (filename(12:12),'(A1)')'b'
                 ELSE
                     WRITE (filename(13:13),'(A1)')'b'
                 ENDIF
				 OPEN(UNIT=11,FILE=filename,STATUS='UNKNOWN')
		 ENDIF
		! write titles
            WRITE (10,'(/''Observation nodes for each time step.'')')
	      IF (Pt(ip)==1) WRITE(10,'(/''Cross-section perpendicular to X axis at X = '',1pE11.4)')xgrid(NodePr(ip,1))
             IF (Pt(ip)==2) WRITE(10,'(/''Cross-section perpendicular to Y axis at Y = '',1pE11.4)')ygrid(NodePr(ip,1))
	      IF (Pt(ip)==3) WRITE(10,'(/''Cross-section  perpendicular to Z axis at Z = '',1pE11.4)')zgrid(NodePr(ip,1))
             IF (Pt(ip)==4) WRITE(10,'(/''Probe location defined by the user'')')
	      IF (DistrP(ip).EQ.1)  THEN
			WRITE(10,'(/''This probe contains'',I4,'' nodes'')') nodebyPr(ip)
                     IF (varP(ip)==2) WRITE (10,'(/''      Time   Averaged water content'')')
		       IF (varP(ip)==1) WRITE (10,'(/''      Time   Averaged water potential'')')
		       IF (varP(ip)==3) WRITE (10,'(/''      Time   Aver. water potential   Aver. water content'')')
		ELSE
		   IF (minl.NE.nodebyPr(ip)) THEN
			WRITE(10,'(/''This probe contains'',I4,'' nodes but only 500 are given in this file, the rest is in the file b.'')')nodebyPr(ip)
			WRITE(form,'I3')500
		   ELSE
		      WRITE(10,'(/''This plane contains'',I4,'' nodes'')') nodebyPr(ip)
                    IF (nodebyPr(ip).GT.99) THEN
                       WRITE(form,'I3')nodebyPr(ip)!only valid for larger than 100 numbers!!!!
                     ELSEIF (nodebyPr(ip).GT.9) THEN
                       WRITE(form2,'I2')nodebyPr(ip)
                    ELSE
                       WRITE(form1,'I1')nodebyPr(ip)
		  	ENDIF
                ENDIF
		  WRITE(10,'(/''Nodes ID'')')
                IF (nodebyPr(ip).GT.99) THEN
                     WRITE(10,'('//form//'(1X,I5))') (NodePr(ip,i),i=1,minl)
                ELSEIF (nodebyPr(ip).GT.9) THEN
                    WRITE(10,'('//form2//'(1X,I5))') (NodePr(ip,i),i=1,minl)
                ELSE
                    WRITE(10,'('//form1//'(1X,I5))') (NodePr(ip,i),i=1,minl)
		  ENDIF
		  IF (varP(ip)==2) WRITE(10,'(/''      Time   Node water content'')')
		  IF (varP(ip)==1) WRITE(10,'(/''      Time   Node water potential'')')
                IF (varP(ip)==3) WRITE(10,'(/''      Time   Node water potential (line 1) and node water content (line 2)'')')
!if more than 500 nodes
			IF (minl.NE.nodebyPr(ip)) THEN
                        WRITE (11,'(/''Observation nodes for each time step.'')')
		 		 IF (Pt(ip)==1) WRITE(11,'(/''plane perpendicular to X axis at X = '',1pE11.4)')Crp(ip)
				 IF (Pt(ip)==2) WRITE(11,'(/''plane perpendicular to Y axis at Y = '',1pE11.4)')Crp(ip)
		 		 IF (Pt(ip)==3) WRITE(11,'(/''plane perpendicular to Z axis at Z = '',1pE11.4)')Crp(ip)
                            IF (Pt(ip)==4) WRITE(11,'(/''probe  location defined by the user'')')
				 WRITE(11,'(/''This probe contains'',I4,'' nodes. Here are given the'',I4,'' nodes larger than 500'')') &
				 nodebyPr(ip),(NodebyPr(ip)-500)
				 WRITE (11,'(/''Nodes ID'')')
				 WRITE(form,'I3')(nodebyPr(ip)-500)
		 		 WRITE(11,'('//form//'(1X,I5))') (NodePr(ip,i),i=501,nodebyPr(ip))
				 IF (varP(ip)==2) WRITE (11,'(/''      Time   Node water content'')')
				 IF (varP(ip)==1) WRITE (11,'(/''      Time   Node water potential'')')
				 IF (varP(ip)==3) WRITE (11,'(/''      Time   Node water potential (line 1) and node water content (line 2)'')')
			    ENDIF
		     ENDIF
          ENDDO
      END SUBROUTINE ObsIni
!****************************************************************************************

      SUBROUTINE WriteLog(t,grwfac,sAvg,cAvg,mshoot,mroot,level)
      Use Typedef
      USE PlntData
      USE Doussanmat,only : count_nodes,PHr,nplant
      IMPLICIT NONE
      REAL(sp):: mshoot,mroot,t,grwfac,sAvg,cAvg
      INTEGER ::level,ipl
      CHARACTER :: na*10,file*8
      
      WRITE (file,'(A7)')'out/log'
      na='    n/a   '
      DO ipl=1,nplant
         WRITE (file(8:8),'(I1)') ipl
         OPEN (UNIT=10,FILE=file,STATUS='OLD',POSITION='APPEND')
         IF (level.EQ.1) WRITE (10,'(1pE10.3,1X,3(A10,1X),2(1pE10.3,1X),A10,1X,1pE10.3)') &
                  t,na,na,na,sAvg,cAvg,na,mroot
         IF (level.EQ.2) WRITE(10,'(3(1pE10.3,1X),A10,1X,2(1pE10.3,1X),A10,1X,1pE10.3)') &
                      t,Tpot,Tact(ipl),na,sAvg,cAvg,na,mroot
         IF (level.EQ.3) WRITE (10,'(8(1pE10.3,1X))')t,Tpot,Tact(ipl),grwfac,sAvg,cAvg,mshoot,mroot
         IF (level.EQ.4) WRITE (10,'(3(1pE11.4,1X),I6,1X,1pE11.4,1X,1pE11.4)')t,Tpot,Tact(ipl),count_nodes,PHr(1,ipl)
      END DO
      RETURN
      END SUBROUTINE WriteLog
!****************************************************************************************
      SUBROUTINE OutFEM(t,hNew,Conc,MatNum,kOuFEM)
      USE GridData
      USE SolData
      USE WatFun
      IMPLICIT NONE
      REAL,intent(in) :: Conc(maxnod)
      REAL(sp),intent(in):: hNew(maxnod),t
      INTEGER(sp), intent(in)::MatNum(maxnod),koufem
      INTEGER(sp):: corner(8),ie,i,m,ic
      CHARACTER file*13
      REAL(sp):: sum,sumc,sume,sumce,tht,solt,tottht,totsol

      WRITE (file,'(A11)')'out/outfem.'
      IF (kOuFEM.LT.10) THEN
       WRITE (file(12:12),'(I1)') kOuFEM
      ELSE
       WRITE (file(12:13),'(I2)') kOuFEM
      ENDIF
      OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
      WRITE (8,*) HeadLn
! calculate total water and solute volume within domain:
      sum=0.0_dp
      sumC=0.0_dp
      DO 1 iE=1,nElm,2
! assign cuboid corner nodes:
       corner(1)=elmnod(iE,1)
       corner(2)=elmnod(iE+1,6)
       corner(3)=elmnod(iE,3)
       corner(4)=elmnod(iE,2)
       corner(5)=elmnod(iE,4)
       corner(6)=elmnod(iE+1,3)
       corner(7)=elmnod(iE,6)
       corner(8)=elmnod(iE,5)
! average cuboid water and solute content:
       sumE=0.0_sp
       sumCE=0.0_sp
       DO 11 ic=1,8
          i=corner(ic)
          M=MatNum(i)
          tht=Fth(hNew(i),par(:,M))
          sumE=sumE+tht
          solt=(Fth(hNew(i),par(:,M))+ChPar(1,M)*ChPar(5,M))*Conc(i)
          sumCe=sumCE+solt
11       CONTINUE
       sum=sum+sumE
       sumC=sumC+sumCE
1     CONTINUE
      tottht= sum*dxGrid*dyGrid*dzGrid/8._sp
      totsol=sumC*dxGrid*dyGrid*dzGrid/8._sp
      WRITE(8,'(/''Total Water Volume at Time  '',F12.4,1X,A5,'' is'',1pE15.8,     ''.'')')&
           t,TmUnit,tottht
      WRITE(8,'(/''Total Solute Volume at Time '',F12.4,1X,A5,'' is'',1pE15.8,     ''.'')')&
           t,TmUnit,totsol
      WRITE (8,'(/''Length Unit is '',A5)')LnUnit
      WRITE (8,'(/''Node# Mater.#'',4X,''x'',9X,''y'',9X,''z'',11X,''h'',9X,''conc.'',7X,''theta'',7X,''wsink'',7X,''csink'')')
! write nodal values:
      DO 2 i=1,nPt
       M=MatNum(i)
       tht=Fth(hNew(i),par(:,M))
       WRITE (8,'(I5,6X,I2,3(1X,1pE9.2),5(1X,1pE11.4))')i,M,xGrid(i),yGrid(i),zGrid(i),hNew(i),Conc(i),tht,sink(i),csink(i)
2     CONTINUE
      CLOSE (8)
      RETURN
      END SUBROUTINE OutFEM
!*************************************************************************
      SUBROUTINE Zprofiles(t,hnew,MatNum)
      USE GridData
      USE SolData
      USE WatFun
      IMPLICIT NONE
      REAL(sp),intent(in) :: t,hNew(maxnod)
      INTEGER(sp), intent(in)::MatNum(maxnod)
      INTEGER(sp):: iz,M,prof,ixy,i
      CHARACTER file*12,form*5
      REAL(sp):: thz(1000),phz(1000),sz(1000),hh
      REAL(dp)::wci
       WRITE(form,'I3')(nx*ny)
!open files
        OPEN (UNIT=121,FILE='out/ProfileTH.out',STATUS='OLD',POSITION='APPEND')
        OPEN (UNIT=122,FILE='out/ProfilePH.out',STATUS='OLD',POSITION='APPEND')
        OPEN (UNIT=123,FILE='out/ProfileS.out',STATUS='OLD',POSITION='APPEND')
!calculate average/sum
      prof=0
      DO iz=1,nz
         thz(iz)=0.
         phz(iz)=0.
         sz(iz)=0.
         DO ixy=1,nx*ny
            i=ixy+prof
            M=MatNum(i)
            thz(iz)=thz(iz)+Fth(hNew(i),par(:,M))
            phz(iz)=phz(iz)+hNew(i)
            sz(iz)=sz(iz)+sink(i)
         ENDDO
         thz(iz)=thz(iz)/(nx*ny)
         phz(iz)=phz(iz)/(nx*ny)
         prof=prof+nx*ny
      ENDDO
         WRITE (121,'('//form//'(1X,F15.6))') t,(thz(i),i=1,nz)
         WRITE (122,'('//form//'(1X,F15.6))') t,(phz(i),i=1,nz)
         WRITE (123,'('//form//'(1X,F15.6))') t,(sz(i),i=1,nz)
close(121)
close(122)
close(123)
      END SUBROUTINE Zprofiles
!*************************************************************************
      SUBROUTINE OutObsProbe(t,hNew,MatNum)
      USE GridData
      USE SolData
      USE WatFun
      USE ObsData
      IMPLICIT NONE
	  INTEGER(sp)::ip,n,M,minl,i
	  INTEGER(sp), intent(in)::MatNum(maxnod)
	  REAL(sp),intent(in) :: t,hNew(maxnod)
	  CHARACTER file*13,form*5
	  REAL(sp):: tht(1000),pht(1000),AverPH,AverTH
      DO ip=1,nPr
!calculations
	 tht=0
	 pht=0
	 DO n=1,nodebyPr(ip)
	    M=MatNum(NodePr(ip,n))
           tht(n)=Fth(hNew(NodePr(ip,n)),par(:,M))
	    pht(n)=hNew(NodePr(ip,n))
!	    print *,n,tht(n),pht(n),ip,NodePr(ip,n)
	ENDDO

!define file names
        WRITE (file,'(A13)')'out/Probe.   '
        IF (ip.LT.10) THEN
          WRITE (file(11:11),'(I1)') ip
        ELSE
          WRITE (file(11:12),'(I2)') ip
        ENDIF
        OPEN (UNIT=10,FILE=file,STATUS='OLD',POSITION='APPEND')
		IF (DistrP(ip).EQ.1) THEN
		!average of each plane
           AverTH=sum(tht)/nodebyPr(ip)
		   AverPH=sum(pht)/nodebyPr(ip)
		   IF (VarP(ip)==1) WRITE (10,'(2(1X,F12.4))') t,AverPH
		   IF (VarP(ip)==2) WRITE (10,'(2(1X,F12.4))') t,AverTH
		   IF (VarP(ip)==3) WRITE (10,'(3(1X,F12.4))') t,AverPH,AverTH
		   close(10)
		ELSE
		!complete distribution asked
           minl=MIN(nodebyPr(ip),500)
		   WRITE(form,'I3')(minl+1)
		   IF (VarP(ip)==1) WRITE (10,'('//form//'(1X,F12.4))') t,(pht(i),i=1,minl)
		   IF (VarP(ip)==2) WRITE (10,'('//form//'(1X,F12.4))') t,(tht(i),i=1,minl)
		   IF (VarP(ip)==3) THEN
             WRITE (10,'('//form//'(1X,F12.4))') t,(pht(i),i=1,minl)
		     WRITE (10,'('//form//'(1X,F12.4))') t,(tht(i),i=1,minl)
		   ENDIF
           CLOSE(10)
!IF more than 500 nodes
           IF (minl.NE.nodebyPr(ip)) THEN
             IF (ip.LT.10) THEN
                WRITE (file(12:12),'(A1)')'b'
             ELSE
		        WRITE (file(13:13),'(A1)')'b'
             ENDIF
!			 print *,'la',pht(399),pht(501),pht(800),nodebyPr(ip)
		     OPEN(UNIT=11,FILE=file,STATUS='OLD',POSITION='APPEND')
		     WRITE(form,'I3')(nodebyPr(ip)-500)+1
		     IF (VarP(ip)==1) WRITE (11,'('//form//'(1X,F12.4))') t,(pht(i),i=501,nodebyPr(ip))
		     IF (VarP(ip)==2) WRITE (11,'('//form//'(1X,F12.4))') t,(tht(i),i=501,nodebyPr(ip))
		     IF (VarP(ip)==3) THEN
               WRITE (11,'('//form//'(1X,F12.4))') t,(pht(i),i=501,nodebyPr(ip))
		       WRITE (11,'('//form//'(1X,F12.4))') t,(tht(i),i=501,nodebyPr(ip))
		     ENDIF
		     CLOSE(11)
           ENDIF
         ENDIF
	  ENDDO
	  RETURN
      END SUBROUTINE OutObsProbe
!*************************************************************************
      SUBROUTINE OutRoo(t,naxes,mroot,mshoot,LA,sAvg,cAvg)
      USE RootData
      USE DoussanMat, only : nplant

      IMPLICIT NONE
      REAL(sp) :: mroot,mshoot,LA,t,sAvg,cAvg,linlim
      INTEGER(sp) :: irec,igrow,ifive,naxes,iestbl,i
      CHARACTER outfile*12

      WRITE (outfile,'(A4)')'out/'
      WRITE (outfile(5:12),'(F7.3)')t
      OPEN (UNIT=8,FILE=outfile,STATUS='UNKNOWN')
      WRITE (8,'(''Time:'')')
      WRITE (8,*) t
      WRITE (8,*)
      WRITE (8,'(''Number of seeds'')')
      WRITE (8,*) nplant
      WRITE (8,*)
      WRITE (8,'(''ID, X and Y coordinates of the seeds (one per line)'')')
      DO i=1,nplant
            WRITE (8,'(I5,3(1pE9.2))') i,xplant(i),yplant(i)
      ENDDO
      WRITE (8,*)
      WRITE (8,'(''Root DM, shoot DM, leaf area:'')')
      WRITE (8,*) mroot,mshoot,LA
      WRITE (8,*)
      WRITE (8,'(''Average soil strength and solute concentration experienced by root system:'')')
      WRITE (8,*) sAvg,cAvg
      WRITE (8,*)
      WRITE (8,'(''Total # of axes:'')')
      WRITE (8,*) naxes
      WRITE (8,*)
      WRITE (8,'(''Total # of branches, including axis(es):'')')
      WRITE (8,*) nbr
      WRITE (8,*)
      WRITE (8,'(''Total # of segment records:'')')
      WRITE (8,*) nrec
      WRITE (8,*)
      WRITE (8,'(''segID#'',4X,''x'',10X,''y'',10X,''z'',6X,''prev or '','' br#  length   surface  mass'')')
      WRITE (8,'(''origination time'')')
! write list of all segment records:
      DO 4 irec=1,nrec
       WRITE (8,'(I5,3(1X,1pE10.3),1X,I5,I2,I5,3(1pE9.2))')&
            irec,xs(irec),ys(irec),zs(irec),irecpr(irec),ordseg(irec),ibrseg(irec),&
            seglen(irec),segsur(irec),segmas(irec)
       WRITE (8,'(1pE11.4)') timorg(irec)
4     CONTINUE
      WRITE (8,*)
      WRITE (8,'(''Total # of growing branch tips:'')')
      WRITE (8,*) ngrow
      WRITE (8,*)
      WRITE (8,'(''tipID#'',4X,''xg'',10X,''yg'',10X,''zg'',6X,''sg.bhd.tp. '',''ord  br#  tot.br.lgth. axs#'')')
      WRITE (8,'(''overlength'',2X,''# of estblished points'')')
      WRITE (8,'(''time of establishing (-->)'')')
! write list of all growing tips:
      DO 5 igrow=1,ngrow
       WRITE(8,'(I5,3(1X,1pE11.4),1X,I5,6X,I2,1X,I5,1X,1pE11.4,3X,I3)')&
            igrow,xg(igrow),yg(igrow),zg(igrow),irecsg(igrow),ordgrw(igrow),&
            ibrgrw(igrow),brlgth(igrow),iaxis(igrow)
       WRITE (8,'(1pE11.4,1X,I5)')ovrlen(igrow),nestbl(igrow)
       IF (nestbl(igrow).GT.0) THEN
          ifive=0
51          linlim=MIN(ifive+5,nestbl(igrow))
          WRITE (8,'(5(1X,1pE13.6))')(timest(igrow,iestbl),iestbl=ifive+1,linlim)
          ifive=ifive+5
          IF (linlim.LT.nestbl(igrow)) GOTO 51
       ENDIF
5     CONTINUE
      CLOSE(8)
      RETURN
      END SUBROUTINE OutRoo
!*********************************************************************
      SUBROUTINE OutDou(t,kOuR)
      USE RootData
      USE DoussanMat, only: Lr,Khr,PHs,PHr,sinkR,axialRootFlow,Qi,nplant,Qd,Q_bc
      USE GridData, only: dxgrid,dygrid,dzgrid
      USE PlntData, only:Tact
      IMPLICIT NONE
      REAL(sp) :: t
      INTEGER(sp) :: irec,kouR,ipl
      CHARACTER file*14
      WRITE (file,'(A12)')'out/outRoot.'

      DO ipl=1,nplant
         WRITE (file(11:11),'(I1)') ipl
         IF (kOuR.LT.10) THEN
            WRITE (file(13:13),'(I1)') kOuR
         ELSE
            WRITE (file(13:14),'(I2)') kOuR
         ENDIF
         OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
         WRITE (8,'(''Time:'')')
         WRITE (8,*) t
         WRITE (8,*)
         WRITE (8,'(''Total # of segment records:'')')
         WRITE (8,*) nrec
         WRITE (8,*)
         WRITE (8,'(A6,3A11,A5,A6,A10,A10,A15,A13,A15,A15,3A5,A15)') &
            'segID#','x','y','z','br#','prev','Lr','Khr','PHinter','PHxylem','radialRootFlow','axialRootFlow','Qi','Qd','Q_bc','  segment radius'
         DO 4 irec=1,nrec
           WRITE (8,'(I6,3(1X,1pE10.3),I5,1X,I6,1X,E9.3,1X,E9.3,1X,E13.5,1X,E13.5,6(1X,E13.5))') irec,xs(irec)+xplant(ipl),ys(irec)+yplant(ipl),&
           zs(irec),ibrseg(irec),irecpr(irec),Lr(irec),Khr(irec),PHs(irec,ipl),PHr(irec+1,ipl),sinkR(irec,ipl),&
           axialRootFlow(irec,ipl),Qi(irec,ipl),Qd(irec,ipl),Q_bc(irec,ipl),segsur(irec)/seglen(irec)/2/pi
4        CONTINUE
         CLOSE(8)
      END DO
      RETURN
      END SUBROUTINE OutDou
!***************************************************************************
      SUBROUTINE SubReg(MatNum,Conc,t,hNew,iCount,lChem,theta,level)
      USE GridData
      USE SolData
      USE WatFun
      USE CumData
      USE Doussanmat, only : SinkR
      IMPLICIT NONE
      INTEGER(sp),intent(in) :: MatNum(maxnod),level
      LOGICAL,intent(in) :: lChem
      REAL(sp),intent(in) :: t,hNew(maxnod),theta(maxnod)
      REAL,intent(in) :: Conc(maxnod)
      INTEGER(sp) :: i,j,k,l,Mi,Mj,Mk,Ml,iSE,iE,icount,iL(4)
      REAL(sp) :: dl,cbalr,cc,cblat
      REAL(sp) ::cnewe,cl,bl,cel,wel
      REAL(dp) ::WatVol,ConVol,ww,wBalT,cBalT,DeltC,DeltW,WNewE
      OPEN (UNIT=10,FILE='out/balance.out',STATUS='OLD',POSITION='APPEND')

! initializing some variables
      WatVol=0.0
      DeltW=0.0
      ConVol=0.0
      DeltC=0.0
      DO 13 iE=1,nElm !loop over elements
         wEl=0.0
         cEl=0.0
         DO 12 iSE=1,3!loop over subelements
            IF (mod(iE,2).eq.0) THEN
               iL(1)=3+iSE
               iL(2)=4+iSE-iSE/3*6
               iL(3)=5+iSE-iSE/2*6
               iL(4)=  iSE
            ELSEIF (mod(iE,2).eq.1) THEN
               iL(1)=3+iSE
               iL(2)=5+iSE-iSE/2*6
               iL(3)=4+iSE-iSE/3*6
               iL(4)=  iSE
            ENDIF
            i=elmnod(iE,iL(1)) !i,j,k,l are corner nodes of the tetraedal element, ise counts 3 which is equal to three tedrahedals which is equal to half a cubic.
            j=elmnod(iE,iL(2))
            k=elmnod(iE,iL(3))
            l=elmnod(iE,iL(4))
            Mi=MatNum(i)
            Mj=MatNum(j)
            Mk=MatNum(k)
            Ml=MatNum(l)
!            Tact = VE(iE,iSE)*(sink(i)+sink(j)+sink(k)+sink(l))/4._dp
            WNewE=VE(iE,iSE)*(theta(i)+theta(j)+theta(k)+theta(l))/4
            WatVol=WatVol+WNewE
            wEl=wEl+WNewE
            IF (lChem) THEN
               CNewE=VE(iE,iSE)*((theta(i)+ChPar(1,Mi)*ChPar(5,Mi))*Conc(i)+(theta(j)+&
                  ChPar(1,Mj)*ChPar(5,Mj))*Conc(j)+(theta(k)+ChPar(1,Mk)*ChPar(5,Mk))*Conc(k)+&
                  (theta(l)+ChPar(1,Ml)*ChPar(5,Ml))*Conc(l))/4
               ConVol=ConVol+CNewE
               cEl=cEl+CNewE
            ENDIF
            IF (iSE.EQ.3) THEN
               IF (iCount.EQ.0) THEN
                  WatIn(iE)=wEl
                  IF (lChem) SolIn(iE)=cEl
               ELSE
                  DeltW=DeltW+abs(WatIn(iE)-wEl)
                  IF (lChem) DeltC=DeltC+abs(SolIn(iE)-cEl)
               ENDIF
            ENDIF
12       CONTINUE
13    CONTINUE
!     Mass balance calculation
      IF (iCount.EQ.0) THEN
         wVolI=WatVol
         IF (lChem) cVolI=ConVol
         IF (level.EQ.4) THEN !new javaux
            WRITE(10,135)
         ELSE
            WRITE(10,130)
         ENDIF
         IF (lChem) THEN
            WRITE(10,140) t,WatVol,ConVol
         ELSE
            WRITE(10,150) t,WatVol
         ENDIF
      ELSE
         wBalT=WatVol-wVolI+wCumT
         ww=max(DeltW,wCumA)
         IF (ww.GE.1.e-25_sp) wBalR=abs(wBalT)/ww*100    !wBalR=abs(wBalT)/WatVol*100._sp!wBalR=abs(wBalT)/ww*100._sp
         write(*,*)'wBalR=',wBalR
         IF (lChem) THEN
            cBalT=ConVol-cVolI+cCumT
            cc=max(DeltC,cCumA)
!javaux      cc=amax1(DeltC,cCumA)
            IF (cc.GE.1.e-25_sp) cBalR=abs(cBalT)/cc*100
         ENDIF
         IF (lChem) THEN
            WRITE(10,160) t,WatVol,wBalT,wBalR,ConVol,cBalT,cBalR,Peclet,Courant
         ELSE
            WRITE(10,170) t,WatVol,wBalT,wBalR,ww,wcumT,WatVol-wvolI,wcumA,deltW,RootSk,sum(SinkR)
         ENDIF
      ENDIF
      CLOSE (10)
      WatVolOld=WatVol
130   FORMAT(/,'    Time [T]   WatVol [V]  WatBalT [V]  WatBalR [%]   CncVol [M]  ',&
               'CncBalT [M]  CncBalR [%]    Peclet       Courant')
135   FORMAT(/,'    Time [T]   WatVol [V]  WatBalT [V]  WatBalR [%]  ww[V]    ',&
               'wCumT [V]  WatVol-wVolI[V]   wCumA[V]    Deltaw[V]   TotalRootFlow [V/T]',&
               '  TotRadialFlow [V/T]')
140   FORMAT(f12.4,1X,1pE12.5,27x,1pE12.5)
150   FORMAT(f12.4,1X,1pE12.5)
160   FORMAT(f12.4,1X,8(1pE12.5,1X))
170   FORMAT(f12.4,1X,10(1pE12.5,1X))
      OPEN (UNIT=10,FILE='out/remove.out',STATUS='OLD',POSITION='APPEND')
!3     READ (10,*,ERR=4)
!      GOTO 3
      IF (iCount.EQ.0) THEN
         WRITE(10,230)
      ELSE
         WRITE(10,240)  t,CumCh0,CumCh1,CumChR,ChemS(1),ChemS(2),CumRt ,CumQ(1) ,CumQ(2)
      ENDIF
230   FORMAT(/,'    Time [T]   CumCh0 [M]   CumCh1 [M]   CumChR [M] ChemS(1) [M]',&
               ' ChemS(2) [M]    CumRt [V]  CumQ(1) [V]  CumQ(2) [V]')
240   FORMAT(f12.4,1X,8(1pE12.5,1X))
      CLOSE(10)
      RETURN
      END SUBROUTINE SubReg
!************************************************************************
      SUBROUTINE FlxOut(kOuFEM,hNew)
      USE typedef
      USE Soldata
      USE GridData
      IMPLICIT NONE
      REAL(sp),intent(in):: hNew(maxnod)
      INTEGER(sp), intent(in)::koufem
      INTEGER(sp):: i
      CHARACTER file*13
      WRITE (file,'(A11)')'out/veloci.'

      IF (kOuFEM.LT.10) THEN
       WRITE (file(12:12),'(I1)') kOuFEM
      ELSE
       WRITE (file(12:13),'(I2)') kOuFEM
      ENDIF
      OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
      WRITE (8,*) HeadLn
      WRITE (8,'(/''Node#'',9X,''x'',9X,''y'',9X,''z'',11X,''Vx'',8X,''Vy'',8X,''Vz'')')
      CALL Veloc(hNew)
! write nodal values:
      DO 2 i=1,nPt
       WRITE (8,'(I5,8X,3(1X,1pE9.2),3(1X,1pE11.4))')i,xGrid(i),yGrid(i),zGrid(i),Vx(i),Vy(i),Vz(i)
2     CONTINUE
      CLOSE(8)
      RETURN
      END SUBROUTINE FlxOut
!****************************************************************************
      SUBROUTINE Getout
      WRITE(*,'(//'' End of simulation. Program terminated normally.''///)')
      STOP
      END SUBROUTINE Getout
! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! Source file ORTHOFEM.FOR |||||||||||||||||||||||||||||||||||||||||||||
!
!                            ORTHOFEM
!
!                           VERSION 1.02
!
!                      FORTRAN SUBROUTINES FOR
!                   ORTHOMIN OR CONJUGATE GRADIENT
!               MATRIX SOLUTION ON FINITE-ELEMENT GRIDS
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                           CARL A. MENDOZA
!
!                      WITH CONTRIBUTIONS FROM:
!                            RENE THERRIEN
!
!                     BASED ON AN ORIGINAL CODE BY:
!                         FRANK W. LETNIOWSKI
!
!              WATERLOO CENTRE FOR GROUNDWATER RESEARCH
!                       UNIVERSITY OF WATERLOO
!                         WATERLOO, ONTARIO
!                          CANADA, N2L 3G1
!
!                    LATEST UPDATE: JANUARY 1991
!
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                   COPYRIGHT (c) 1989, 1990, 1991
!                            E.A. SUDICKY
!                            C.A. MENDOZA
!              WATERLOO CENTRE FOR GROUNDWATER RESEARCH
!
!          DUPLICATION OF THIS PROGRAM, OR ANY PART THEREOF,
!          WITHOUT THE EXPRESS PERMISSION OF THE COPYRIGHT
!          HOLDERS IS STRICTLY FORBIDDEN
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!               Modified for SWMS_3D code by Jirka Simunek
!                            august 1994
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                             DISCLAIMER
!
!     ALTHOUGH GREAT CARE HAS BEEN TAKEN IN PREPARING THIS CODE
!     AND THE ACCOMPANYING DOCUMENTATION, THE AUTHORS CANNOT BE
!     HELD RESPONSIBLE FOR ANY ERRORS OR OMISSIONS.  THE USER IS
!     EXPECTED TO BE FAMILIAR WITH THE FINITE-ELEMENT METHOD,
!     PRECONDITIONED ITERATIVE TECHNIQUES AND FORTRAN PROGRAMING.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!             A USER'S GUIDE IS AVAILABLE -> CONSULT IT!
!
!     THESE SUBROUTINES SOLVE A BANDED (OR SPARSE) MATRIX USING:
!       - PRECONDITIONED ORTHOMIN FOR ASYMMETRIC MATRICES, OR
!       - PRECONDITIONED CONJUGATE GRADIENT FOR SYMMETRIC MATRICES
!         (FULL MATRIX STORAGE REQUIRED)
!
!     PRECONDITIONING IS BY INCOMPLETE LOWER-UPPER DECOMPOSITION
!       - ONLY ONE FACTORIZATION (GAUSSIAN ELIMINATION) IS PERFORMED
!       - EQUIVALENT TO DKR FACTORIZATION
!
!     THE SUBROUTINES ARE DESIGNED FOR FINITE-ELEMENT GRIDS
!       - ARBITRARY ELEMENT SHAPES AND NUMBERING MAY BE USED
!         - NUMBERING MAY, HOWEVER, AFFECT EFFICIENCY
!           - TRY TO MINIMIZE THE BANDWIDTH AS MUCH AS POSSIBLE
!         - ALL ELEMENTS MUST HAVE THE SAME NUMBER OF LOCAL NODES
!
!
!     THE FOLLOWING ROUTINES ARE CALLED FROM THE SOURCE PROGRAM:
!       IADMAKE (IN,NN,NE,NLN,MNLN,MAXNB,IAD,IADN,IADD)
!         -> ASSEMBLE ADJACENCY MATRIX
!       FIND (I,J,K,NN,MAXNB,IAD,IADN)
!         -> LOCATE MATRIX POSITION FOR A NODAL PAIR (ASSEMBLY)
!       ILU (R,NN,MAXNB,IAD,IADN,IADD,B)
!         -> DECOMPOSE GLOBAL MATRIX
!       ORTHOMIN (R,C,GT,NNR,MAXNB,MAXNN,IAD,IADN,IADD,B,VRV,
!                 RES,RQI,RQ,Q,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,
!                 NORTH,MNORTH,MAXIT)
!         -> SOLVE DECOMPOSED MATRIX
!
!     THESE ROUTINES CALL OTHER ROUTINES (LOCATED DIRECTLY BELOW THE
!     APPROPRIATE PRIMARY ROUTINE IN THE CODE)
!
!     THE FOLLOWING ARRAYS MUST BE DEFINED IN THE SOURCE PROGRAM
!     (THESE ARRAYS ARE PASSED TO THE SOLVER SUBROUTINES):
!
!     IN(MNLN,MAXNE) - INCIDENCE MATRIX (ELEMENTAL NODE DEFINITION)
!
!     GT(MAXNN)      - RIGHT-HAND-SIDE VECTOR
!     C(MAXNN)       - SOLUTION VECTOR
!     R(MAXNB,MAXNN) - GLOBAL MATRIX TO BE SOLVED
!
!     ARRAY DIMENSIONING PARAMETERS
!
!     MAXNN  - MAXIMUM NUMBER OF NODES
!     MAXNE  - MAXIMUM NUMBER OF ELEMENTS
!     MNLN   - MAXIMUM NUMBER OF LOCAL NODES (IN AN ELEMENT)
!     MAXNB  - MAXIMUM NUMBER OF NODES ADJACENT TO A PARTICULAR NODE
!              (INCLUDING ITSELF).
!            - IE. THE MAXIMUM NUMBER OF INDEPENDENT NODES THAT A
!              PARTICULAR NODE SHARES AN ELEMENT WITH.
!            - THIS WILL BE IDENTICALLY EQUIVALENT TO THE MAXIMUM
!              NUMBER OF NONZERO ENTRIES IN A ROW OF THE FULL MATRIX.
!     MNORTH - MAXIMUM NUMBER OF ORTHOGONALIZATIONS PERFORMED
!              (AT LEAST MNORTH = 1 REQUIRED FOR CONJUGATE GRADIENT)
!
!
!     ORTHOMIN ARRAY SPACE/VARIABLES
!
!     NORTH  - NUMBER OF ORTHOGONALIZATIONS TO PERFORM
!            - SET NORTH=0 FOR SYMMETRIC MATRICES (CONJUGATE GRADIENT)
!     ECNVRG - RESIDUAL CONVERGENCE TOLERANCE
!     ACNVRG - ABSOLUTE CONVERGENCE TOLERANCE
!     RCNVRG - RELATIVE CONVERGENCE TOLERANCE
!     MAXIT  - MAXIMUM NUMBER OF ITERATIONS TO PERFORM
!     ITERP  - NUMBER OF ITERATIONS PERFORMED
!
!     B(MAXNB,MAXNN) - ILU DECOMPOSED MATRIX
!     Q(MAXNN)   - SEARCH DIRECTION Q
!     RQ(MAXNN)  - PRODUCT OF R AND Q
!     VRV(MAXNN) - EITHER V OR PRODUCT OF R AND V
!     RES(MAXNN) - RESIDUAL
!
!     QI(MAXNN,MNORTH)  - STORAGE OF Q'S
!     RQI(MAXNN,MNORTH) - STORAGE OF PRODUCTS OF R AND Q
!     RQIDOT(MNORTH)    - STORAGE OF DOT PRODUCTS OF RQ AND RQ
!
!     RESV   - PREVIOUS VALUE OF RES V DOT PRODUCT (CONJUGATE GRADIENT)
!
!     IAD(MAXNB,MAXNN) - ADJACENCY MATRIX (NODAL CONNECTIONS)
!
!     IADN(MAXNN) - NUMBER OF ADJACENT NODES IN IAD (SELF-INCLUSIVE)
!
!     IADD(MAXNN) - POSITION OF DIAGONAL IN ADJACENCY MATRIX
!
!
!     OTHER PARAMETERS PASSED FROM SOURCE PROGRAM
!
!     NN  - NUMBER OF NODES
!     NE  - NUMBER OF ELEMENTS
!     NLN - NUMBER OF LOCAL NODES IN AN ELEMENT
!
!
!     APPROXIMATE REAL STORAGE SPACE FOR ORTHOMIN AND MATRIX EQUATION
!
!       ((6 + 2!MAXNB + 2!MNORTH)!MAXNN)!(8 BYTES)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IADMake(NumNP,NumEl,NumEld,MaxNB,IAD,IADN,IADD)
      use GridData
      use Matdata, only: numnz,irow,jcol,A_sparse
!     Generate the adjacency matrix for nodes from the element
!     indidence matrix
!     Requires subroutine Insert
      implicit none
      !dimension IAD(MaxNB,NumNP),IADN(NumNP),IADD(NumNP),iL(4)
      !integer e
      integer(sp) :: NumNP,NumEl,NumEld,MaxNB,i,j,k,l,kk,iSE
      integer(sp) :: IAD(MaxNB,NumNP),IADN(NumNP),IADD(NumNP),iL(4),iE
!     Determine independent adjacency within each element
!     version for SWMS_3D
      do 12 i=1,NumNP
        IADN(i)=0
        IADD(i)=0
        do 11 j=1,MaxNB
          IAD(j,i)=0
11      continue
12    continue
      do 14 iE=1,NumEl
!       Loop on subelements
        DO 13 iSE=1,3
          IF (mod(iE,2).eq.0) THEN
             iL(1)=3+iSE
             iL(2)=4+iSE-iSE/3*6
             iL(3)=5+iSE-iSE/2*6
             iL(4)=  iSE
          ELSEIF (mod(iE,2).eq.1) THEN
             iL(1)=3+iSE
             iL(2)=5+iSE-iSE/2*6
             iL(3)=4+iSE-iSE/3*6
             iL(4)=  iSE
          ENDIF
          i=elmnod(iE,iL(1))
          j=elmnod(iE,iL(2))
          k=elmnod(iE,iL(3))
          l=elmnod(iE,iL(4))
          call Insert(i,j,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(j,i,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(i,k,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(k,i,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(j,k,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(k,j,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(i,l,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(l,i,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(j,l,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(l,j,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(k,l,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(l,k,kk,NumNP,MaxNB,IAD,IADN)
13      continue
14    continue
!number of nonzero entries
      numnz=sum(IADN)+nPt !adjacencies + diagonal values
      JCOL=0
      IROW=0
      !first value is a zero in IROW
      IROW(1)=1
      do 15 i=1,NumNP
        call Insert(i,i,kk,NumNP,MaxNB,IAD,IADN)
        IADD(i)=kk
        IROW(1+i)=IROW(i)+IADN(i) !previous  irow value (at current i) + #entries in that row
        do j=1,IADN(i)
           JCOL(IROW(i)-1+j)=IAD(j,i) ! column value added to JCOL
        enddo
15    continue
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE INSERT (I,J,K,NN,MAXNB,IAD,IADN)
!     ADD J TO THE ADJACENCY LIST FOR I
!     RETURNS THE POSITION K WHERE IT HAS BEEN ADDED, OR WHERE IT
!     WAS ALREADY IN THE LIST.
      USE typedef
      implicit none
      INTEGER(sp) :: IAD(MAXNB,NN),IADN(NN)
      INTEGER(sp) :: I,J,K,NN,MAXNB,N,L,INODE
      !DIMENSION IAD(MAXNB,NN),IADN(NN)
      LOGICAL :: FOUND
      FOUND = .FALSE.
!     DETERMINE NUMBER OF NODES ALREADY IN ADJACENCY LIST
      N = IADN(I)
      K = N + 1
!     DETERMINE WHETHER ALREADY IN LIST
      DO 10 L=1,N
        INODE = IAD(L,I)
        IF (INODE.GE.J) THEN
          K = L
          IF (INODE.EQ.J) FOUND = .TRUE.
          GO TO 15
        ENDIF
   10 CONTINUE
   15 CONTINUE
!     PLACE IN LIST (NUMERICAL ORDER)
      IF (FOUND) THEN
       CONTINUE
      ELSE
        IF ((N+1).GT.MAXNB) THEN
          WRITE (* ,601) I,MAXNB
          WRITE (50,601) I,MAXNB
  601     FORMAT (//5X,'ERROR IN IADMAKE: NODE ',I5,' HAS > ' &
                 ,I5,' ADJACENCIES')
          STOP
        ENDIF
        IADN(I) = N + 1
        DO 20 L=(N+1),(K+1),(-1)
          IAD(L,I) = IAD(L-1,I)
   20   CONTINUE
       IAD(K,I) = J
      ENDIF
!write(*,*)'door insert'
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FIND (I,J,K,NN,MAXNB,IAD,IADN)
!     FOR NODE I, DETERMINE THE 'BAND' (K) RELATED TO ITS ADJACENCY TO
!     NODE J.
!     IF NODE NOT ADJACENT, RETURN 0 AS THE 'BAND'
      use typedef
!      DIMENSION IAD(MAXNB,NN),IADN(NN)
      integer(sp) :: IAD(MAXNB,NN),IADN(NN),NN,MAXNB,I,J,K,N,INODE,L
      K = 0
      N = IADN(I)
      DO 10 L=1,N
        INODE = IAD(L,I)
!     EXIT THE LOOP IF AT OR PAST THE REQUIRED POSITION
        IF (INODE.GE.J) THEN
          IF (INODE.EQ.J) K = L
          GO TO 20
        ENDIF
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ILU (R,NN,MAXNB,IAD,IADN,IADD,B)
!     INCOMPLETE LOWER-UPPER DECOMPOSITION OF MATRIX R INTO B
!     ONE STEP OF GAUSSIAN ELIMINATION PERFORMED
!     DIAGONAL DOMINANCE IS ASSUMED - NO PIVOTING PERFORMED
!     REQUIRES FUNCTION DU
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION R(MAXNB,NN),IAD(MAXNB,NN),IADN(NN),IADD(NN),B(MAXNB,NN)
      use typedef
      implicit none
      integer(sp) :: MAXNB,NN,I,J,N,K,ICUR,INODE,L,IAD(MAXNB,NN),IADN(NN),IADD(NN)
      real(sp) :: SUM,R(MAXNB,NN),D,DU
      real(sp), intent(out) :: B(MAXNB,NN)

!     INITIALIZE B
    B(1:MAXNB,1:NN)=0.
write(*,*)'nn=',NN
write(*,*)'maxnb=',maxnb
write(*,*)'sizeR=',size(R,1),size(R,2)



!     LOOP OVER NODES
      DO 20 I=1,NN
!     DETERMINE NUMBER OF BANDS/POSITION OF DIAGONAL IN THIS ROW
        N = IADN(I)
        K = IADD(I)
!write(*,*)'N=',N
!write(*,*)'K=',K
!     LOWER TRIANGULAR MATRIX
!write(*,*)'ok'
!write(*,*)'sizeR=',size(R,1),size(R,2)
        DO 30 J=1,(K-1)
!write(*,*)'j=',J
!write(*,*)'R --',R(1:11,i)
          SUM = R(J,I)
!write(*,*)'ok1'
          ICUR = IAD(J,I)
          DO 40 L=1,(J-1)
            INODE = IAD(L,I)
!write(*,*)'ok2'
            SUM = SUM - B(L,I)*DU(INODE,ICUR,NN,MAXNB,IAD,IADN,IADD,B)
!write(*,*)'ok3'
   40     CONTINUE
          B(J,I) = SUM
!write(*,*)'B(j,i)=',B(J,I)
   30   CONTINUE
!write(*,*)'lower'
!     DIAGONAL
        SUM = R(K,I)
        DO 50 L=1,(K-1)
          INODE = IAD(L,I)
          SUM = SUM - B(L,I)*DU(INODE,I,NN,MAXNB,IAD,IADN,IADD,B)
   50   CONTINUE
        D = 1.0D0/SUM
        B(K,I) = D
!write(*,*)'diagonal'
!     UPPER TRIANGULAR MATRIX
!       - ACTUALLY D*U TO OBTAIN UNIT DIAGONAL
        DO 60 J=(K+1),N
          SUM = R(J,I)
          ICUR = IAD(J,I)
          DO 70 L=1,(K-1)
            INODE = IAD(L,I)
            SUM = SUM - B(L,I)*DU(INODE,ICUR,NN,MAXNB,IAD,IADN,IADD,B)
   70     CONTINUE
          B(J,I) = D*SUM
   60   CONTINUE
!write(*,*)'upper'
!write(*,*)'I,NN',I,NN
   20 CONTINUE
write(*,*)'blaat'
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION DU (I,INODE,NN,MAXNB,IAD,IADN,IADD,B)
!     SEARCHES THE I'TH ROW OF THE UPPER DIAGONAL MATRIX
!     FOR AN ADJACENCY TO THE NODE 'INODE'
!     RETURNS CORRESPONDING VALUE OF B (OR ZERO)
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION IAD(MAXNB,NN),IADN(NN),IADD(NN),B(MAXNB,NN)
      use typedef
      implicit none
      integer(sp) :: MAXNB,NN,I,J,N,K,INODE,IAD(MAXNB,NN),IADN(NN),IADD(NN)
      real(sp) :: B(MAXNB,NN),TEMP,DU
      TEMP = 0.0D0
      N = IADN(I)
      K = IADD(I)
      IF (I.EQ.INODE) THEN
        TEMP = 1.0D0
        GO TO 20
      ENDIF
      DO 10 J=(K+1),N
        IF (INODE.EQ.IAD(J,I)) THEN
          TEMP = B(J,I)
          GO TO 20
        ENDIF
   10 CONTINUE
   20 CONTINUE
      DU = TEMP
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      SUBROUTINE ORTHOMIN (R,C,GT,NN,MAXNB,MAXNN,IAD,IADN,IADD,B,VRV,&
!                          RES,RQI,RQ,Q,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,&
!                          NORTH,MNORTH,MAXIT)
       SUBROUTINE OrthoMin(R,C,GT,NN,MAXNB,MAXNN,IAD,IADN,IADD,B,NORTH)
!     ORTHOMIN OR CONJUGATE GRADIENT ACCELERATION/SOLUTION
!     CONJUGATE GRADIENT (SYMMETRIC MATRIX) IF NORTH=0
!     (HOWEVER, NOTE THAT MNORTH MUST BE AT LEAST 1)
!     REQUIRES FUNCTIONS SDOT,SDOTK,SNRM
!     REQUIRES SUBROUTINES LUSOLV,MATM2,SAXPYK,SCOPY,SCOPYK
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION R(MAXNB,NN),C(NN),GT(NN),IAD(MAXNB,NN),IADN(NN)
!      DIMENSION IADD(NN),B(MAXNB,NN),VRV(NN),RES(NN),RQI(MAXNN,MNORTH)
!      DIMENSION RQ(NN),Q(NN),RQIDOT(MNORTH),QI(MAXNN,MNORTH)
       use typedef
       implicit none
       INTEGER(sp) ::  MNorth=4
       INTEGER(sp), intent(in):: IAD(MAXNB,NN),IADN(NN),IADD(NN),MAXNB,MAXNN,NN
       real(sp),intent(in) :: B(MAXNB,NN),GT(NN),R(MAXNB,NN)
       real(sp), intent(out) :: C(NN)
       REAl(sp) :: VRV(NN),RES(NN),RQI(MAXNN,4),RQ(NN)
       REAL(sp) :: Q(NN),QI(MAXNN,4),RQIDOT(4),ECNVRG=1e-6,RCNVRG=1e-6,ACNVRG=1e-6
       integer(sp) :: MaxIt=1000,North
       INTEGER(sp) :: iHomog,I,norcur,iter,K,ITERp
       REAL(sp) :: DOT,ALPHA,OMEGA,RQNORM,RESMAX,DXNORM,CNORM,XRATIO,SDOT,SDOTK,SNRM2,RESV
!     INITIALIZE RESIDUAL VECTOR
      CALL MATM2 (RES,R,C,NN,IAD,IADN,MAXNB)
!write(*,*)'res',res
!pause
!     Solution for homogeneous system of equations - Modified by Simunek
      iHomog=0
      DO 10 I=1,NN
        RES(I) = GT(I) - RES(I)
        if(abs(GT(i)).gt.1.d-300) iHomog=1
   10 CONTINUE
      if(iHomog.eq.0) then
        do 11 i=1,NN
	  C(i)=GT(i)
11      continue
        return
      end if
!write(*,*)'C',C
!pause
!     LOOP OVER A MAXIMUM OF MAXIT ITERATIONS
      NORCUR = 0
      DO 100 ITER=1,MAXIT
!     INVERT LOWER/UPPER MATRICES
        CALL SCOPY (NN,RES,VRV)
        CALL LUSOLV (NN,MAXNB,IAD,IADN,IADD,B,VRV)
!write(*,*)'VRV id not okay then ILU wrong',VRV
!pause
!     COPY V INTO Q
        CALL SCOPY (NN,VRV,Q)
!     CALCULATE PRODUCT OF R AND V
        CALL MATM2 (VRV,R,Q,NN,IAD,IADN,MAXNB)
!     COPY RV INTO RQ
        CALL SCOPY (NN,VRV,RQ)
!     RES V DOT PRODUCT (CONJUGATE GRADIENT)
        IF (NORTH.EQ.0) THEN
          DOT = SDOT(NN,RES,Q)
          IF (NORCUR.EQ.0) RESV = DOT
        ENDIF
!write(*,*)'DOT=',DOT
!pause
!     LOOP OVER PREVIOUS ORTHOGONALIZATIONS
        K = 1
   20   IF (K.GT.NORCUR) GO TO 30
!     DETERMINE WEIGHTING FACTOR (CONJUGATE GRADIENT)
          IF (NORTH.EQ.0) THEN
            ALPHA = DOT/RESV
            RESV = DOT
!write(*,*)'Alpha,RESV',ALPHA,RESV
!pause
!     DETERMINE WEIGHTING FACTOR (ORTHOMIN)
          ELSE
            DOT = SDOTK(NN,K,RQI,VRV,MAXNN,MNORTH)
            ALPHA = -DOT/RQIDOT(K)
          ENDIF
!     SUM TO OBTAIN NEW Q AND RQ
          CALL SAXPYK (NN,ALPHA,K,QI,Q,MAXNN,MNORTH)
          CALL SAXPYK (NN,ALPHA,K,RQI,RQ,MAXNN,MNORTH)
!write(*,*)'QI=',QI(1:NN,1)
!pause
!write(*,*)'RQI=',RQI(1:NN,1)
!pause
          K = K + 1
          GO TO 20
   30   CONTINUE
!     CALCULATE WEIGHTING FACTOR (CONJUGATE GRADIENT)
        IF (NORTH.EQ.0) THEN
          DOT = SDOT(NN,Q,RQ)
          OMEGA = RESV/DOT
!write(*,*)'omega=',OMEGA
!     CALCULATE WEIGHTING FACTOR (ORTHOMIN)
        ELSE
          DOT = SDOT(NN,RES,RQ)
          RQNORM = SDOT(NN,RQ,RQ)
          OMEGA = DOT/RQNORM
        ENDIF
!     SAVE VALUES FOR FUTURE ORTHOGONALIZATIONS
!write(*,*)'Q=',Q
!write(*,*)'RQ=',RQ
!pause
        NORCUR = NORCUR + 1
        IF (NORCUR.GT.NORTH) NORCUR = 1
        CALL SCOPYK (NN,NORCUR,Q,QI,MAXNN,MNORTH)
        CALL SCOPYK (NN,NORCUR,RQ,RQI,MAXNN,MNORTH)
        RQIDOT(NORCUR) = RQNORM
! write(*,*)'QI2=',QI(1:NN,1)
!pause
!write(*,*)'RQI2=',RQI(1:NN,1)
!pause
!     UPDATE SOLUTION/RESIDUAL VECTORS
        CALL SAXPYK (NN,OMEGA,NORCUR,QI,C,MAXNN,MNORTH)
        CALL SAXPYK (NN,-OMEGA,NORCUR,RQI,RES,MAXNN,MNORTH)
!     DETERMINE CONVERGENCE PARAMETERS
        RESMAX = SNRM2(NN,RES)
        DXNORM = ABS(OMEGA)*SNRM2(NN,Q)
        CNORM = SNRM2(NN,C)
        XRATIO = DXNORM/CNORM
!write(*,*)'resmax',resmax
!write(*,*)'DXNORM',dxnorm
!write(*,*)'cnorm',cnorm
!write(*,*)'xratio',xratio
!     ITERATION (DEBUG) OUTPUT
!     STOP ITERATING IF CONVERGED
        IF (RESMAX.LT.ECNVRG) GO TO 200
        IF (XRATIO.LT.RCNVRG) GO TO 200
        IF (DXNORM.LT.ACNVRG) GO TO 200
  100 CONTINUE
!     TERMINATE IF TOO MANY ITERATIONS
      WRITE (* ,602) MAXIT,RESMAX,XRATIO,DXNORM
      WRITE (70,602) MAXIT,RESMAX,XRATIO,DXNORM
 602  FORMAT (///5X,'ORTHOMIN TERMINATES -- TOO MANY ITERATIONS', &
               /8X,'MAXIT  = ',I5, &
               /8X,'RESMAX = ',E12.4, &
               /8X,'XRATIO = ',E12.4, &
               /8X,'DXNORM = ',E12.4)
      STOP
  200 CONTINUE
!     RETURN NUMBER OF ITERATIONS REQUIRED
      ITERP = ITER
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE LUSOLV (NN,MAXNB,IAD,IADN,IADD,B,VRV)
!     LOWER DIAGONAL MATRIX INVERSION BY FORWARD SUBSTITUTION
!     UPPER DIAGONAL MATRIX INVERSION BY BACKWARD SUBSTITUTION
!     LOWER/UPPER MATRICES ARE IN B
!     RIGHT-HAND-SIDE VECTOR IS IN VRV AT START
!     'SOLUTION' IS RETURNED IN VRV UPON EXIT
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION IAD(MAXNB,NN),IADN(NN),IADD(NN),B(MAXNB,NN),VRV(NN)
        use typedef
        implicit none
       INTEGER(sp):: IAD(MAXNB,NN),IADN(NN),IADD(NN),MAXNB,NN
       real(sp) :: B(MAXNB,NN),VRV(NN),SUM
       integer(sp) :: I,K,J,INODE,N
!     LOWER INVERSION
      DO 20 I=1,NN
        SUM = VRV(I)
        K = IADD(I)
        DO 30 J=1,(K-1)
          INODE = IAD(J,I)
          SUM = SUM - B(J,I)*VRV(INODE)
   30   CONTINUE
        VRV(I) = B(K,I)*SUM
   20 CONTINUE
!     UPPER INVERSION
      DO 40 I=NN,1,-1
        SUM = VRV(I)
        N = IADN(I)
        K = IADD(I)
        DO 50 J=(K+1),N
          INODE = IAD(J,I)
          SUM = SUM - B(J,I)*VRV(INODE)
   50   CONTINUE
        VRV(I) = SUM
   40 CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MATM2 (S1,R,P,NN,IAD,IADN,MAXNB)
!     MULTIPLY MATRIX R BY VECTOR P TO OBTAIN S1
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION S1(NN),P(NN),R(MAXNB,NN),IAD(MAXNB,NN),IADN(NN)
        use typedef
        implicit none
       INTEGER(sp):: IAD(MAXNB,NN),IADN(NN),MAXNB,NN
       real(sp) :: P(NN),SUM,R(MAXNB,NN),S1(NN)
       integer(sp) :: I,J,INODE,N
      DO 30 I=1,NN
        SUM = 0.0D0
        N = IADN(I)
        DO 40 J=1,N
          INODE = IAD(J,I)
          SUM = SUM + R(J,I)*P(INODE)
   40   CONTINUE
        S1(I) = SUM
   30 CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SDOT (NN,R,B)
!     OBTAIN DOT PRODUCT OF R AND B
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION R(NN),B(NN)
        use typedef
        implicit none
        real(sp) :: SDOT,R(NN),B(NN)
        integer(sp) :: NN,L
      SDOT = 0.0D0
      DO 100 L=1,NN
        SDOT = SDOT + R(L)*B(L)
100   CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SDOTK (NN,K,R,B,MAXNN,MNORTH)
!     OBTAIN DOT PRODUCT OF R AND B
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION R(MAXNN,MNORTH),B(NN)
        use typedef
        implicit none
        real(sp) :: SDOTK,B(NN),R(MAXNN,MNorth)
        integer(sp) :: NN,L,K,MNORTH,MAXNN
      SDOTK = 0.0D0
      DO 100 L=1,NN
        SDOTK = SDOTK + R(L,K)*B(L)
100   CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SNRM2 (NN,R)
!     COMPUTE MAXIMUM NORM OF R
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION R(NN)
      use typedef
      integer(sp) :: L,NN
      real(sp) :: R(NN),TEMP,SNRM2
      SNRM2 = 0.0D0
      DO 100 L=1,NN
        TEMP = ABS(R(L))
        IF (TEMP.GT.SNRM2) SNRM2 = TEMP
100   CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SAXPYK (NN,SA,K,FX,FY,MAXNN,MNORTH)
!     MULTIPLY VECTOR FX BY SCALAR SA AND ADD TO VECTOR FY
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION FX(MAXNN,MNORTH),FY(NN)
      use typedef
      integer(sp) :: I,NN,K,MAXNN,MNORTH
      REAL(sp):: FX(MAXNN,MNORTH),FY(NN),SA
      IF (NN.GT.0) THEN
        DO 100 I=1,NN
           FY(I) = SA*FX(I,K) + FY(I)
100     CONTINUE
      ENDIF
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SCOPY (NN,FX,FY)
!     COPY VECTOR FX INTO VECTOR FY
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION FX(NN),FY(NN)
        use typedef
        integer(sp) :: I,NN
        REAL(sp):: FX(NN),FY(NN)
      IF (NN.GT.0) THEN
        DO 100 I=1,NN
           FY(I) = FX(I)
100     CONTINUE
      ENDIF
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SCOPYK (NN,K,FX,FY,MAXNN,MNORTH)
!     COPY VECTOR FX INTO VECTOR FY
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION FX(N),FY(MAXNN,MNORTH)
        use typedef
      integer(sp) :: I,NN,K,MAXNN,MNORTH
      REAL(sp):: FY(MAXNN,MNORTH),FX(NN)
      IF (NN.GT.0) THEN
        DO 100 I=1,NN
           FY(I,K) = FX(I)
100     CONTINUE
      ENDIF
      RETURN
      END

