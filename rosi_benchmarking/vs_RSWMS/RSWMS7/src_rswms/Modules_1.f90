!************************************************************************
!> general information, definition of max. matrix sizes for pre-allocation
MODULE ParamData
  USE typedef
  IMPLICIT NONE
  INTEGER(ap),PARAMETER :: maxnod=662661,mxpnts=100,maxplant=2
  INTEGER(ap),PARAMETER :: maxbdr=13000,mxBcCh=13000,maxIrrig=5,mxBcCh1=50
  INTEGER(ap),PARAMETER :: mxtime=1000,mxdpth=100,maxelm= 640000
  INTEGER(ap),PARAMETER :: maxgrw=10000,maxrec=200000,maxest=5000
  INTEGER(ap),PARAMETER :: maxemg=50,maxord=3,maxmat=10,maxbnd=19
  INTEGER(ap),PARAMETER :: maxobs=30
  INTEGER(ap),PARAMETER :: maxParticle=10000 !solute particles inside root
  INTEGER(ap),PARAMETER:: maxAMG=5000000,maxamg2=5000000 !array size for sparse A matrix and sparse row index matrix
  INTEGER(ap) :: icount,iter_root,iter_tot=0,iter=5,last_out
  REAL(dp), PARAMETER :: pi=3.14159265358979323
  LOGICAL :: lvtk=.FALSE., lOutPartrace=.FALSE., ldirect=.FALSE.,lChem=.FALSE.,lretry=.FALSE.
  LOGICAL :: lOrt, lSalinity=.true.,lPartUp=.FALSE.
END MODULE ParamData
!************************************************************************
!> data needed for plant
!> root characteristics and Doussan input variables
MODULE PlntData
  USE typedef
  USE ParamData, ONLY: mxBcCh, mxBcCh1, maxplant
  IMPLICIT NONE
  INTEGER(ap) :: nTpot,ntTpLA,nfTpLA,ncTpLA,ntLA,ntRSR,nscRSR,t_ini,t_dev,t_mid,t_end
  INTEGER(ap) :: ntW,nsfW,nscW,ncnc,ns,nBCr,nsfRSR,typeBCr(mxBcCh)
  REAL(sp) :: tTpot(mxBcCh),Tpotc(mxBcCh)
  REAL(sp) :: tTpLA(mxBcCh),TpLAc(mxBcCh),LA,hlim,sf
  REAL(sp) :: sfTpLA(mxBcCh),fTpLAc(mxBcCh),scTpLA(mxBcCh)
  REAL(sp) :: tLA(mxBcCh),LAc(mxBcCh),cTpLAc(mxBcCh)
  REAL(sp) :: h50,p50,p1,p2,CMm,VMax,fk,xin,h0,h1,h2,h3
  REAL(sp) :: tRSR(mxBcCh),RSRc(mxBcCh)
  REAL(sp) :: sfRSR(mxBcCh1),fRSRc(mxBcCh1),scRSR(mxBcCh1)
  REAL(sp) :: sc(mxBcCh1),rsc(mxBcCh1),cncp(mxBcCh1),rscnc(mxBcCh1)
  REAL(sp) :: tW(mxBcCh1),Wc(mxBcCh1),cRSRc(mxBcCh1)
  REAL(sp) :: sfW(mxBcCh1),fWc(mxBcCh1),scW(mxBcCh1),cWc(mxBcCh1)
  REAL(sp) :: tBCr(mxBcCh),BCroot(mxBcCh)!Doussan in BC
  REAL(sp) :: TpLA, TpLA_pot
  REAL(sp) :: a_r,a_1,a_2
  REAL(dp) :: Tpot,Tact(maxplant),TotSur=0._dp,PHcollar,b_1,b_2,SpWgt
END MODULE PlntData
!************************************************************************
!> data needed for root
MODULE RootData
  USE typedef
  USE paramData, ONLY: maxord,maxplant,maxrec,maxgrw,maxord,maxest,mxpnts,maxemg
  IMPLICIT NONE
  INTEGER(ap) :: nbr,nrec,ngrow,naxemg,nvch(maxord),nMPLch(maxord)
  REAL(sp) :: xplant(maxplant),yplant(maxplant),condMP
  INTEGER(ap) :: ordseg(maxrec),ibrseg(0:maxrec),iaxis(maxgrw),num_seg(maxgrw),br_rec(maxgrw)
  INTEGER(ap) :: nUrf(maxord),nAQPc
  INTEGER(ap) :: naxes,norder
  real(8) :: rand !random number created in RSWMS34.f90
  !>INTEGER(ap), PARAMETER :: seed=7654321._sp
  INTEGER(ap) :: irecpr(maxrec),irecsg(maxgrw),nnewax(maxemg)
  INTEGER(ap) :: nestbl(maxgrw),ordgrw(maxgrw),ibrgrw(maxgrw),lPast
  REAL(dp) :: seglen(maxrec),segsur(0:maxrec),segmas(maxrec),segdiam(maxrec),brlgth(maxgrw),Krs,Kcomp,segrad(maxrec)
  REAL(dp) :: segSoluteMass(maxrec), segconc(maxrec), crossSectionSeg(maxrec), segvol(maxrec) !for SoluteRoot
  REAL(dp) :: sign_in,PH_crit,sAvg,cAvg,vol_root,vol_buff,size_buff,m_in,brlmax(maxord),f_rad(maxord)!,ovrlen(maxgrw)
  REAL(sp) :: xs(maxrec),ys(maxrec),zs(maxrec),timorg(maxrec),g1,g2,ovrtime(maxgrw)
  REAL(sp) :: xg(maxgrw),yg(maxgrw),zg(maxgrw)
  REAL(sp) :: agevch(maxord,mxpnts)
  REAL(sp) :: timest(maxgrw,maxest),vch(maxord,mxpnts)
  REAL(dp) :: sMPLch(maxord,mxpnts),MPLch(maxord,mxpnts),MPL(maxgrw)=0._dp
  REAL(sp) :: age(maxord,mxpnts),Urf(maxord,mxpnts)
  REAL(sp) :: tnewax(maxemg),inaxs(maxemg),strsen(maxord),rdmang(maxord)
  REAL(sp) :: brspac(maxord-1),brnang(maxord-1),dtbrch(maxord-1),tlim_dJvL,OmegaC
  REAL(sp) :: delbm,LAmsh,mroot,mshoot
  REAL(sp), ALLOCATABLE, DIMENSION(:) :: Rho,tPast,AQPh,AQPv
  REAL(sp), ALLOCATABLE, DIMENSION(:,:) :: sinkredDPast
  REAL(dp) :: mhorm(maxrec,maxplant),concol,mcol=0._dp,msign_notrans,csign_notrans,res_t,delta_h,mcol_i=0._dp
  LOGICAL :: maizeroottyp=.FALSE.,loliumroottyp=.FALSE.,wheatroottyp=.FALSE.,lrrt=.FALSE.,lrrs=.FALSE.,lSomma=.FALSE.
  LOGICAL :: lDou=.FALSE.,lCou=.FALSE.,lSUF=.FALSE.,lFed=.FALSE.,ldJvL=.FALSE.,lJarvis=.FALSE.,lUrf=.FALSE.
  LOGICAL :: lGap=.FALSE.,lAQPc=.FALSE.,l_secrad=.FALSE.,l_conduc=.FALSE.,lno_RWU=.FALSE.,stopgr(maxgrw)=.FALSE.
  LOGICAL :: lno_Archi=.FALSE.,lno_root_growth=.FALSE.,lRootTyp_growth=.FALSE.,lSomma_growth=.FALSE.,lUpdate_growth=.FALSE.
  LOGICAL :: lCalloc=.FALSE.,ltemp=.FALSE.,lSign=.FALSE.,lSign_inst=.FALSE.,ltoxi=.FALSE.,lSinkCube=.FALSE.
  LOGICAL(sp), ALLOCATABLE, DIMENSION(:) :: toosml,connex,l_SignOn!,stopgr
  LOGICAL :: it1
  REAL(sp) :: rsr,dr,dsh,grwfac,w,dmroot,rs,concrs,Hseq
  REAL(sp) :: sigma=1.0   
END MODULE RootData
!************************************************************************
!> Data needed for the soil grid 
MODULE GridData
  USE ParamData, ONLY: maxElm
  USE typedef
  IMPLICIT NONE
  REAL(sp), ALLOCATABLE, DIMENSION(:) :: xgrid,ygrid,zgrid,Axy,Bxy,Dxy,Exy,betaw,betac,Width
  REAL(sp), ALLOCATABLE, DIMENSION(:) :: Vn,VElm,HElm,HElmOld,RLD,RSD
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: SSF
  REAL(sp), ALLOCATABLE,DIMENSION(:) :: Wn,sink,csink,sinkOld,sink_cube,sink_cubeOld,csink_cube,betac_cube,betac_cube2
  REAL(sp), ALLOCATABLE,DIMENSION(:,:) :: Deter
  REAL(dp), ALLOCATABLE,DIMENSION(:,:,:) :: B1fact,Ax,Ay,Az
  REAL(sp), ALLOCATABLE,DIMENSION(:,:,:) :: bi,ci,di
  REAL(sp), ALLOCATABLE,DIMENSION(:,:,:,:) :: E
  INTEGER(ap), ALLOCATABLE,DIMENSION (:,:) :: elmnod
  INTEGER(ap), ALLOCATABLE,DIMENSION (:) :: subN,n_neigh
  INTEGER(ap), ALLOCATABLE, DIMENSION (:) :: IADN,IADD
  INTEGER(ap), ALLOCATABLE, DIMENSION (:,:) :: IAD,MacroList  
  REAL(sp) :: dxGrid,dyGrid,dzGrid,epslonPH,epslonWC,epslonR,factorRelEps,epslonS,dxSSF,dySSF,dzSSF,dxRLD,dyRLD,dzRLD,dxRho,dyRho,dzRho
  INTEGER(ap) ::nPt,nElm,nel,nBand,nBCPts,nex,ney,nez,itMax,itMaxRoot,nx,ny,nz,nexSSF,neySSF,nezSSF,nexRLD,neyRLD,nezRLD,nexRho,neyRho,nezRho,geom
  CHARACTER LnUnit*5,TmUnit*5,MsUnit*5,CnUnit*5
  REAL(sp) ::xCol(1000),yCol(1000),x_cent=0,y_cent=0
  REAL(sp) :: RootSk,checkSink,RootSkOld,rad_cyl=0
  LOGICAL :: RelEps,continu
  !subelement order of first cube->  subN(iE) = 1
  !subelement order of second  cube->  subN(iE) = 2
  !iL(corner_tetraheda,subelements,subN(iE))
  INTEGER(sp), PARAMETER :: iL(1:4,1:5,1:2) = RESHAPE([1, 4, 7, 3, &
       1, 6, 4, 2, &
       5, 6, 7, 1, &
       7, 6, 8, 4, &
       1, 6, 7, 4, &
       5, 2, 3, 1, &
       5, 6, 8, 2, &
       5, 8, 7, 3, &
       5, 2, 8, 3, &
       2, 8, 3, 4 ],SHAPE(iL))
  INTEGER(ap),SAVE::iadd_temp(maxElm*80)
  
CONTAINS
  SUBROUTINE IniGrid
    IMPLICIT NONE
    ALLOCATE (xgrid(1:nPt))
    ALLOCATE (ygrid(1:nPt))
    ALLOCATE (zgrid(1:nPt))
    ALLOCATE (Vn(1:nPt))
    ALLOCATE (VElm(1:nElm))
    ALLOCATE (HElm(1:nElm))
    ALLOCATE (HElmOld(1:nElm))
    ALLOCATE (subN(1:nPt))
    ALLOCATE (Axy(1:nPt))
    ALLOCATE (Bxy(1:nPt))
    ALLOCATE (Dxy(1:nPt))
    ALLOCATE (Exy(1:nPt))
    ALLOCATE (betaw(1:nPt))
    ALLOCATE (betac(1:nPt))
    ALLOCATE (sink(1:nPt))
    sink=0._dp
    ALLOCATE (sinkOld(1:nPt))
    sinkOld=0._dp
    ALLOCATE (sink_cube(1:nElm))
    sink_cube=0._dp
    ALLOCATE (sink_cubeOld(1:nElm))
    sink_cubeOld=0._dp
    ALLOCATE (csink_cube(1:nElm))
    csink_cube=0._dp
    ALLOCATE (betac_cube(1:nElm))
    betac_cube=0._dp
    ALLOCATE (betac_cube2(1:nElm))
    betac_cube2=0._dp
    ALLOCATE (csink(1:nPt))
    ALLOCATE (Wn(1:nPt))
    ALLOCATE (Width(1:nPt))
    ALLOCATE (elmnod(1:8,1:nElm))
    ALLOCATE (Deter(1:5,1:nElm))
    ALLOCATE (B1fact(1:4,1:5,1:nElm))
    ALLOCATE (Ax(1:4,1:5,1:nElm))
    ALLOCATE (Ay(1:4,1:5,1:nElm))
    ALLOCATE (Az(1:4,1:5,1:nElm))
    ALLOCATE (bi(1:4,1:5,1:nElm))
    ALLOCATE (ci(1:4,1:5,1:nElm))
    ALLOCATE (di(1:4,1:5,1:nElm))
    ALLOCATE (E(1:4,1:4,1:5,1:nElm))
    
  END SUBROUTINE IniGrid
END MODULE GridData
!*************************************************************************
!> Data needed for the Doussan root water uptake module
MODULE DoussanMat
  USE typedef
  USE paramData, ONLY: mxBcCh1
  USE GridData
  USE RootData, ONLY: ngrow,nrec,maxrec,maxgrw,lSomma_growth
  USE SparseMatrix !multiple roots
  IMPLICIT NONE
  
  ! Doussan matrices
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: PHr_sub,PH_root_sub,axialRootFlow,veloRoot
  REAL(sp), ALLOCATABLE, DIMENSION (:,:) :: k_ave,B,Qd,Qi,Q_bc1,Q_bc2,Q_bc
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: Lr,Khr,GH,qroot,l_seg,Khr_pot,Lr_pot,Inv_c1KxKr,Inv_ciiKrKr
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: Phi_mat,h_mat2,Joutr
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: delta2,delta2old,sinkRtemp
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: tempSinkR,curr_BCr,BCr_usr
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: PHs,PHr,PHrOld,PHrTemp,PHsTemp,PHo
  REAL(sp), ALLOCATABLE, DIMENSION (:,:) ::sinkR,SinkROld
  REAL(dp), ALLOCATABLE, DIMENSION (:) :: Phs_osmotic !osmotic head
  REAL(dp) :: KhRoot(1:3,mxBcCh1),LrRoot(1:3,mxBcCh1),Jintot,sinktot
  REAL(sp) :: ageKh(1:3,mxBcCh1),ageLr(1:3,mxBcCh1),hx_min,stresval1,stresval2,cavitb,cavitc
  REAL(sp), ALLOCATABLE, DIMENSION (:,:,:) :: PH_micro2
  REAL(sp), ALLOCATABLE, DIMENSION (:,:,:,:) :: w_dis,cent,cp_mean,Intc
  INTEGER(ap), ALLOCATABLE, DIMENSION (:,:,:) :: cube_i
  REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: w_sub,l_sub,sum_dis, beta_weight
  INTEGER(ap) :: nBCn,nKh(1:3),nLr(1:3),nmax,nrecOld,isubmax=25,nplant=1,solveroot_call=0
  INTEGER(ap) :: nLibr=10000,indexValue=1000,switchcriterion=1,n=50
  INTEGER(ap), ALLOCATABLE, DIMENSION(:) :: counter2,BCtp_usr,curr_BCtp
  INTEGER(ap), ALLOCATABLE, DIMENSION(:) :: nBC_irecn,nBC_iprvn,no_voxels
  INTEGER(ap), ALLOCATABLE, DIMENSION(:,:) :: numNodes_voxel,voxel_no,nsub
  INTEGER(ap), ALLOCATABLE, DIMENSION (:,:,:) :: voxel_node
  INTEGER(ap), ALLOCATABLE, DIMENSION (:,:,:,:) :: loc_Q,transroot,transtip
  INTEGER(ap), ALLOCATABLE, DIMENSION(:) ::iro,jco,jao,iao
  REAL(dp), ALLOCATABLE, DIMENSION (:)::aij,ao
  
  INTEGER(ap)::count_nodes,cavitfun, stresfun
  LOGICAL :: loop1=.TRUE.,stressBC=.FALSE.,ave,old,oldT,eqDis,ItCrit_root=.TRUE.,once=.TRUE.
  LOGICAL :: moment=.FALSE.,savelast=.FALSE.,switchSolve,tcheck=.FALSE.
  LOGICAL :: tcheck2=.FALSE.,ana_aan
  TYPE(SparseMatrixType), ALLOCATABLE, DIMENSION (:) :: plantmatrix !multiple roots
  
CONTAINS
  SUBROUTINE IniMat
    USE RootData, ONLY :nrec
    IMPLICIT NONE
    IF (loop1) THEN
       loop1=.FALSE.
    ELSE
       DEALLOCATE (voxel_no)
       DEALLOCATE (no_voxels)
       DEALLOCATE(numNodes_voxel)
       DEALLOCATE(voxel_node)
       DEALLOCATE(PH_micro2)
       DEALLOCATE(Phi_mat)
       DEALLOCATE(h_mat2)
       DEALLOCATE(delta2)
       DEALLOCATE(delta2old)
       DEALLOCATE(w_dis)
       DEALLOCATE(cube_i)
       DEALLOCATE(cent)
       DEALLOCATE(cp_mean)
       DEALLOCATE(Intc)
       DEALLOCATE(tempSinkR)
       DEALLOCATE(nsub)
       DEALLOCATE(w_sub)
       DEALLOCATE(beta_weight)
       DEALLOCATE(l_sub)
       DEALLOCATE(l_seg)
       DEALLOCATE(GH)
       DEALLOCATE(Joutr)
       DEALLOCATE(veloRoot)
       DEALLOCATE(axialRootFlow)
       DEALLOCATE(sinkR)
       DEALLOCATE(sinkRtemp)
       DEALLOCATE(Lr)
       DEALLOCATE(Lr_pot)
       DEALLOCATE(Khr)
       DEALLOCATE(Khr_pot)
       DEALLOCATE(PHs)
       DEALLOCATE(PHo)
       DEALLOCATE(PHs_osmotic)
       DEALLOCATE(PHr)
       DEALLOCATE(PHrOld)
       DEALLOCATE(SinkROld)
       DEALLOCATE(PHrTemp)
       DEALLOCATE(PHr_sub)
       DEALLOCATE(PHsTemp)
       DEALLOCATE(PH_root_sub)
       DEALLOCATE(loc_Q)
       DEALLOCATE(Qi)
       DEALLOCATE(Q_bc)
       DEALLOCATE(Q_bc1)
       DEALLOCATE(Q_bc2)
       DEALLOCATE(Qd)
       DEALLOCATE(sum_dis)
       DEALLOCATE(transroot)
       DEALLOCATE(transtip)
       DEALLOCATE(BCtp_usr)
       DEALLOCATE(curr_BCtp)
       DEALLOCATE(curr_BCr)
       DEALLOCATE(BCR_usr)
       DEALLOCATE(k_ave)
       DEALLOCATE(B)
       DEALLOCATE(plantmatrix)
       IF (ave) THEN
          DEALLOCATE(counter2)
       ENDIF
    ENDIF
    IF (ana_aan) THEN
       IF (ave) THEN

          ALLOCATE(counter2(1))
          counter2=0
       ENDIF
       ALLOCATE(PH_micro2(1,1,1)) !1:500 -> length r is lower than 500, so no. of compartments lower
       PH_micro2=0
       ALLOCATE (Phi_mat(1))
       Phi_mat = 0._sp
       ALLOCATE (h_mat2(1))
       h_mat2 = 0._sp
    ELSE
       IF (ave) THEN
          ALLOCATE(counter2(1:Nelm))

          counter2=0
       ENDIF
       ALLOCATE(PH_micro2(0:nrec,1:500,1:isubmax)) !1:500 -> length r is lower than 500, so no. of compartments lower
       PH_micro2=0
       ALLOCATE(numNodes_voxel(0:nrec,1:isubmax))
       numNodes_voxel=0
       ALLOCATE (Phi_mat(1:nLibr-1))
       Phi_mat = 0._sp
       ALLOCATE (h_mat2(1:nLibr-1))
       h_mat2 = 0._sp
    ENDIF
    IF (.NOT.(old)) THEN
       IF ((ave) .OR. (eqdis)) THEN
          ALLOCATE (voxel_no(0:nrec,1:isubmax))
          voxel_no=0
          ALLOCATE (no_voxels(1:nElm))
          no_voxels=0
          ALLOCATE (voxel_node(1:nElm,1:2,indexValue))
          voxel_node=0
          ALLOCATE(numNodes_voxel(0:nrec,1:isubmax))
          numNodes_voxel=0
       ELSE
          ALLOCATE(numNodes_voxel(1,1))
          numNodes_voxel=0
          ALLOCATE (voxel_no(1,1))
          voxel_no=0
          ALLOCATE (no_voxels(1))
          no_voxels=0
          ALLOCATE (voxel_node(1,1,1))
          voxel_node=0
       ENDIF
    ELSE
       ALLOCATE(numNodes_voxel(1,1))
       numNodes_voxel=0
       ALLOCATE (voxel_no(1,1))
       voxel_no=0
       ALLOCATE (no_voxels(1))
       no_voxels=0
       ALLOCATE (voxel_node(1,1,1))
       voxel_node=0
    ENDIF
    ALLOCATE(delta2(1:nrec+1))
    delta2=0
    ALLOCATE(delta2old(1:nrec+1))
    delta2old=0
    ALLOCATE(tempSinkR(1:nrec+1))
    tempSinkR=0
    ALLOCATE (l_seg(0:nrec))
    l_seg=0._dp
    ALLOCATE (w_sub(0:nrec,1:isubmax,1:nplant))
    w_sub=0._dp
    ALLOCATE (beta_weight(0:nrec,1:isubmax,1:nplant))
    beta_weight=0._dp
    ALLOCATE (l_sub(0:nrec,1:isubmax,1:nplant))
    l_sub=0._dp
    ALLOCATE (PHr_sub(0:nrec,1:isubmax))
    PHr_sub=0._dp
    ALLOCATE (PH_root_sub(0:nrec,1:isubmax))
    PH_root_sub=0._dp
    ALLOCATE (GH(0:nrec))
    GH=0._dp
    ALLOCATE (sinkR(0:nrec,1:nplant))
    sinkR=0._dp
    ALLOCATE (BCtp_usr(1:nplant))
    BCtp_usr=0
    ALLOCATE (curr_BCtp(1:nplant))
    curr_BCtp=0
    ALLOCATE (curr_BCr(1:nplant))
    curr_BCr=0
    ALLOCATE (BCr_usr(1:nplant))
    BCr_usr=0
    ALLOCATE (sinkRtemp(0:nrec))
    sinkRtemp=0._dp
    ALLOCATE (Joutr(0:nrec))
    Joutr=0._dp
    ALLOCATE (veloRoot(0:nrec,1:nplant))
    veloRoot=0._dp
    ALLOCATE (axialRootFlow(0:nrec,1:nplant))
    axialRootFlow=0._dp
    ALLOCATE (Lr(0:nrec))
    Lr=0._dp
    ALLOCATE (Lr_pot(0:nrec))
    Lr_pot=0._dp
    ALLOCATE (Khr(0:nrec))
    Khr=0._dp
    ALLOCATE (Khr_pot(0:nrec))
    Khr_pot=0._dp
    ALLOCATE (PHs(0:nrec,1:nplant))
    PHs=0._dp
    ALLOCATE (PHo(0:nrec,1:nplant))
    PHo=0._dp
    ALLOCATE (PHs_osmotic(nElm))
    PHs_osmotic=0._dp
    ALLOCATE (PHr(1:nrec+1,1:nplant))
    PHr=0._dp
    ALLOCATE (PHrOld(1:nrec+1,1:nplant))
    PHrOld=0._dp
    ALLOCATE (SinkROld(0:nrec,1:nplant))
    SinkROld=0._dp
    ALLOCATE (PHrTemp(1:nrec+1,1:nplant))
    PHrTemp=0._dp
    ALLOCATE (PHsTemp(0:nrec,1:nplant))
    PHsTemp=0._dp
    IF(lSomma_growth) THEN
       ALLOCATE (transroot(0:maxrec+maxgrw,1:2,1:isubmax,1:nplant))
       transroot=0
       ALLOCATE (transtip(0:maxgrw,1:2,1:isubmax,1:nplant))
       transtip=0
       ALLOCATE (nsub(0:maxrec+maxgrw,1:nplant))
       nsub=1
    ELSEIF(.NOT. lSomma_growth) THEN
       ALLOCATE (nsub(0:nrec+ngrow,1:nplant))
       nsub=1
       ALLOCATE (transroot(0:nrec+ngrow,1:2,1:isubmax,1:nplant))
       transroot=0
       ALLOCATE (transtip(0:ngrow,1:2,1:isubmax,1:nplant))
       transtip=0
    END IF
    ALLOCATE (loc_Q(0:nrec,1:8,1:isubmax,1:nplant))
    loc_Q=0
    ALLOCATE (w_dis(0:nrec,1:8,1:isubmax,1:nplant))
    w_dis=0._dp
    ALLOCATE (cube_i(0:nrec,1:isubmax,1:nplant))
    cube_i=0
    ALLOCATE (cent(0:nrec,1:3,1:isubmax,1:nplant))
    cent=0._dp
    ALLOCATE (cp_mean(0:nrec,1:3,1:isubmax,1:nplant))
    cp_mean=0._dp
    ALLOCATE (Intc(0:nrec+ngrow,1:3,1:isubmax,1:nplant))
    Intc=0._dp
    ALLOCATE (Qi(0:nrec,1:nplant))
    Qi=0._dp
    ALLOCATE (Q_bc(0:nrec,1:nplant))
    Q_bc=0._dp
    ALLOCATE (Q_bc1(0:nrec,1:nplant))
    Q_bc1=0._dp
    ALLOCATE (Q_bc2(0:nrec,1:nplant))
    Q_bc2=0._dp
    ALLOCATE (Qd(0:2*nrec+1,1:nplant))
    Qd=0._dp
    ALLOCATE (sum_dis(0:nrec,1:isubmax,1:nplant))
    sum_dis=0._dp
    ALLOCATE (k_ave(0:nrec,1:isubmax))
    k_ave=0._dp
    ALLOCATE (B(0:nrec,1:isubmax))
    B=0._dp
    ALLOCATE (plantmatrix(1:nplant))
  END SUBROUTINE IniMat
END MODULE DoussanMat
!*******************************************************************
!> Data needed for the solute transport inside the roots
MODULE SoluteRootMat
  USE typedef
  USE ParamData, only:  maxParticle,maxrec,mxBcCh1
  USE RootData, ONLY: nrec
  IMPLICIT NONE 

  Type Particle
     INTEGER(AP) :: ID
     INTEGER(ap) :: segNum
     REAL(dp) :: positionOld
     REAL(dp) :: position
     REAL(dp) :: mass
     REAL(dp) :: segLen
     REAL(dp) :: partOrig
     TYPE(Particle), POINTER :: prev
     TYPE(Particle), POINTER :: next
  END type Particle

  TYPE(Particle), POINTER :: firstP, pParticle
  LOGICAL :: loop2=.TRUE.,l_linSorb=.false.,l_freundSorb=.false.
  REAL(sp) :: particelMass(maxParticle),agePr(1:3,mxBcCh1),agePs(1:3,mxBcCh1),theta_R,rho_R, segsorb(maxrec),timestepfactor
  REAL(dp) :: PrRoot(1:3,mxBcCh1),PsRoot(1:3,mxBcCh1),sorp(2)=0, seg_upt(maxrec)
  INTEGER(ap)  :: totalParticleNum, irecfollow(maxrec,4),numfollow(maxrec),uptakeorder
  INTEGER(ap)  :: segNumParticle(maxParticle),nPerm(1:3),nPass(1:3)
 ! REAL(ap)  ::  particlePosition(maxParticle), particlePositionOld(maxParticle)
  REAL(dp), allocatable, DIMENSION (:) :: ccube,mupt_adv,mupt_diff, m_upt, fact,mass_elm
  REAL(dp), allocatable, DIMENSION (:) :: Perm,retard,frac_pass, segsolv, theta_elm,Tottransfer

  CONTAINS
    SUBROUTINE IniSolute(t)
      USE typedef
      USE RootData, ONLY: nrec,timorg,ordseg
      USE GridData, ONLY: nElm
      IMPLICIT NONE
      
      REAL(sp), INTENT(in) :: t
      REAL(sp) :: segage
      INTEGER(ap) :: iage,irecn,typ


      IF (.NOT. ALLOCATED(ccube)) THEN
         ALLOCATE(ccube(nElm))
         ccube = 0._dp
      END IF

      IF (.NOT. ALLOCATED(mass_elm)) THEN
         ALLOCATE(mass_elm(nElm))
         mass_elm = 0._dp
      END IF

      IF (.NOT. ALLOCATED(fact)) THEN
         ALLOCATE(fact(nElm))
        fact = 1._dp
      END IF

      IF (.NOT. ALLOCATED(Tottransfer)) THEN
         ALLOCATE(Tottransfer(1))
         Tottransfer=0._dp
      END IF

      IF (.NOT. ALLOCATED(m_upt)) THEN
         ALLOCATE(m_upt(nElm))
         m_upt = 0._dp
      END IF
      IF (.NOT. ALLOCATED(mupt_adv)) THEN
         ALLOCATE(mupt_adv(nElm))
         mupt_adv = 0._dp
      END IF
      IF (.NOT. ALLOCATED(mupt_diff)) THEN
         ALLOCATE(mupt_diff(nElm))
         mupt_diff = 0._dp
      END IF
      IF (.NOT. ALLOCATED(theta_elm)) THEN
         ALLOCATE(theta_elm(nElm))
         theta_elm = 0._dp
      END IF
      IF(.NOT. ALLOCATED(retard)) THEN
         ALLOCATE(retard(maxrec))
         retard = 1._dp
      END IF
     ! IF(.NOT. ALLOCATED(segsorb)) THEN
     !    ALLOCATE(segsorb(nrec))
     !    segsorb = 0._dp
     ! END IF
      IF(.NOT. ALLOCATED(segsolv)) THEN
         ALLOCATE(segsolv(nrec))
         segsolv = 0._dp
      END IF
      ! root permeabilities 
      IF(ALLOCATED(Perm)) DEALLOCATE(Perm)
      ALLOCATE(Perm(nrec))
      Perm = 0._dp
      IF(ALLOCATED(frac_pass)) DEALLOCATE(frac_pass)
      ALLOCATE(frac_pass(nrec))
      frac_pass = 0._dp

      !> root permeability matrices
      DO irecn=1,nrec
         segage = t-timorg(irecn)
         typ = ordseg(irecn)
         iage=1
         DO WHILE (agePr(typ,iage).le.segage)
            iage=iage+1
         ENDDO
         IF (iage>nPerm(typ)) THEN
            Perm(irecn)=PrRoot(typ,nPerm(typ))
         ELSEIF (iage==1) THEN
            Perm(irecn)=PrRoot(typ,iage)
         ELSE
            Perm(irecn)=PrRoot(typ,iage-1)+(PrRoot(typ,iage)-PrRoot(typ,iage-1))*(segage-agePr(typ,iage-1))/&
                 (agePr(typ,iage)-agePr(typ,iage-1))
         ENDIF
      END DO

      DO irecn=1,nrec
         segage = t-timorg(irecn)
         typ = ordseg(irecn)
         iage=1
         DO WHILE (agePs(typ,iage).le.segage)
            iage=iage+1
         ENDDO
         IF (iage>nPass(typ)) THEN
            frac_pass(irecn)=PsRoot(typ,nPass(typ))
         ELSEIF (iage==1) THEN
            frac_pass(irecn)=PsRoot(typ,iage)
         ELSE
            frac_pass(irecn)=PsRoot(typ,iage-1)+(PsRoot(typ,iage)-PsRoot(typ,iage-1))*(segage-agePs(typ,iage-1))/&
                 (agePs(typ,iage)-agePs(typ,iage-1))
         ENDIF
      END DO  

    END SUBROUTINE IniSolute
  
END MODULE SoluteRootMat
!*******************************************************************
MODULE NumericalRecipes
  USE Typedef
  
CONTAINS
  
  SUBROUTINE linbcg(n,b,x,itol,tol,itmax,iter,err,ipl)
    ! solve inverse matrix with preconditionned biconjugate gradient method      
    ! code taken from Numerical recipes in fortran 77, p. 79
    USE SparseMatrix
    USE DoussanMat, ONLY: plantmatrix
    IMPLICIT NONE
    
    INTEGER(ap), INTENT(in) :: itmax,itol,n,ipl
    INTEGER(ap), INTENT(out) :: iter
    REAL(dp), INTENT(inout):: tol,b(:),x(:)
    REAL (dp), INTENT(out) ::err
    INTEGER (sp) :: j
    REAL(dp) :: ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm,EPS
    REAL(dp) :: p(n),pp(n),r(n),rr(n),z(n),zz(n)
    PARAMETER (EPS=1.E-14_dp)
    
    iter=0     
    CALL SM_multiply(plantmatrix(ipl), x, r)
    DO 11 j=1,n
       r(j)=b(j)-r(j)
       rr(j)=r(j)
11     CONTINUE
       znrm=1._dp
       IF (itol.EQ.1) THEN
          bnrm=snrm(n,b,itol)
       ELSE IF (itol.EQ.2) THEN
          CALL SM_divide_by_diagonal(plantmatrix(ipl), b,z)
          bnrm=snrm(n,z,itol)
       ELSE IF (itol.EQ.3.OR.itol.EQ.4) THEN
          CALL SM_divide_by_diagonal(plantmatrix(ipl), b,z)
          bnrm=snrm(n,z,itol)
          CALL SM_divide_by_diagonal(plantmatrix(ipl), r,z)
          znrm=snrm(n,z,itol)
       ELSE
	  PRINT *,'illegal itol in linbcg'
ENDIF
CALL SM_divide_by_diagonal(plantmatrix(ipl), r,z)     

100 IF (iter.LE.itmax) THEN
   iter=iter+1
   zm1nrm=znrm
   CALL SM_divide_by_diagonal(plantmatrix(ipl), rr,zz)
   bknum=0._dp
   DO 12 j=1,n
      bknum=bknum+z(j)*rr(j)
12    CONTINUE
      IF(iter.EQ.1) THEN
         DO 13 j=1,n
            p(j)=z(j)
            pp(j)=zz(j)
13       CONTINUE
      ELSE
           bk=bknum/bkden
           DO 14 j=1,n
              p(j)=bk*p(j)+z(j)
              pp(j)=bk*pp(j)+zz(j)
14         CONTINUE
      ENDIF
      bkden=bknum
     ! call atimes(n,p,z,0)
      CALL SM_multiply(plantmatrix(ipl), p, z)
      akden=0._dp
      DO 15 j=1,n
          akden=akden+z(j)*pp(j)
15    CONTINUE
      ak=bknum/akden
      !         call atimes(n,pp,zz,1)
               
      CALL SM_multiply_transpose(plantmatrix(ipl), pp, zz)
      DO 16 j=1,n
         x(j)=x(j)+ak*p(j)
         r(j)=r(j)-ak*z(j)
         rr(j)=rr(j)-ak*zz(j)
16    CONTINUE  
      CALL SM_divide_by_diagonal(plantmatrix(ipl), r,z)
      IF(itol.EQ.1.OR.itol.EQ.2)THEN
         znrm=1._dp
         err=snrm(n,r,itol)/bnrm
      ELSE IF(itol.EQ.3.OR.itol.EQ.4)THEN
         znrm=snrm(n,z,itol)
      IF(ABS(zm1nrm-znrm).GT.EPS*znrm) THEN
          dxnrm=ABS(ak)*snrm(n,p,itol)
          err=znrm/ABS(zm1nrm-znrm)*dxnrm
      ELSE
         err=znrm/bnrm
         GOTO 100
      ENDIF
         xnrm=snrm(n,x,itol)
      IF(err.LE.0.5_dp*xnrm) THEN
         err=err/xnrm
      ELSE
         err=znrm/bnrm
         GOTO 100
      ENDIF
      ENDIF
      IF(err.GT.tol) THEN
        GOTO 100
      ENDIF
      ENDIF	
      RETURN      
    END SUBROUTINE linbcg
!  (C) Copr. 1986-92 Numerical Recipes Software '%12'%*Sim+).
!***************************************************************************
    FUNCTION snrm(n,sx,itol)
!> compute one or two norms of a vector sx(1:n), as signaled by itol. used by linbcg
!> from Numerical Recipes in F., p.81      
      INTEGER(ap):: n,itol,i,isamax
      REAL(dp):: sx(:),snrm
      
      IF (itol.LE.3)THEN
         snrm=0._dp
         DO 11 i=1,n
            snrm=snrm+sx(i)**2
11       CONTINUE
         snrm=SQRT(snrm)
      ELSE
         isamax=1
         DO 12 i=1,n
            IF(ABS(sx(i)).GT.ABS(sx(isamax))) isamax=i
12       CONTINUE
         snrm=ABS(sx(isamax))
      ENDIF
      RETURN
      END FUNCTION snrm

 END MODULE NumericalRecipes
!========================================================================
!> Data needed for time control
 MODULE tmctrl
   USE typedef
   IMPLICIT NONE
   INTEGER(ap), PARAMETER :: mxOut=450,mxProf=1000
   INTEGER(ap) :: nOut,nouProf,nouProbe
   INTEGER(ap) :: kout=0,kouprof=0,kouprobe=0,kbcr
   INTEGER(ap) ::kaxemg=0
   REAL(sp) :: tOut(mxOut),touProf(mxProf),touProbe(mxProf)
   REAL(sp) :: dtroot=1000.,dtMin,dtMax,FacInc,FacDec,tmax,t_begin,dtProf,dtProbe
   REAL(sp) :: dtMaxC=1.E+30_dp,tcallr,told,dtold,tpulse,dtopt
   REAL(sp) :: tcBCr,tProf,tProbe
   REAL(dp) :: time0=0., time1=0.,time2=0.,time3=0.
   LOGICAL :: tlevel,tlevel_soil
 END MODULE tmctrl
 !*******************************************************************
 !> Concentration data
 MODULE ConData
   USE typedef
   USE ParamData, ONLY: maxnod
   IMPLICIT NONE
   REAL(sp):: impc(maxnod),coptma,cmax,cmin,coptmi
 END MODULE ConData
 !********************************************************************
 !> data for boundary conditions
 MODULE BoundData
   USE typedef
   USE ParamData, ONLY: mxBcCh, maxbdr, maxIrrig
   IMPLICIT NONE
   REAL(sp) :: cBound(2),CBnd2(mxBcCh)
   REAL(sp) :: tQbcCh(mxBcCh),thbcCh(mxBcCh),hbc(mxBcCh)
   REAL(sp) :: tCBnd1(mxBcCh),CBnd1(mxBcCh),tCBnd2(mxBcCh)
   REAL(sp) :: xqmin1,xqmin2,xqmax1,xqmax2
   REAL(sp), DIMENSION (:,:) :: Qbc(mxBcCh,mxBcCh),tIbcCh(maxIrrig,mxBcCh),Ibc(maxIrrig,mxBcCh)
   INTEGER(ap) :: iBCPt(maxbdr+maxIrrig),nQbcCh,nIbcCh,nCBnd1,nhbcCh,nCBnd2,qfun,homogene=1
 END MODULE BoundData
 !********************************************************************
 !> data for temperature
 MODULE TempData
   USE ParamData, ONLY: mxtime, mxdpth, maxnod
   USE typedef
   IMPLICIT NONE
   REAL(sp):: time_S(mxtime),depth(mxdpth),time_TA(mxtime),time_PA(mxtime),time_PD(mxtime),tstart
   REAL(sp):: tem(maxnod),impt(maxnod),T_atm(mxtime),P_atm(mxtime),P_diff(mxtime)
   REAL(sp):: tempermin,topt,tempermax,trange,tmid,expo
   REAL(sp):: Tatm_usr,Patm_usr,Pdiff_usr
   REAL(sp), DIMENSION (:,:) :: temtim(mxtime,mxdpth)
   INTEGER(ap):: nz_tempS,nt_tempS,nt_tempA,nt_presA,nt_presD
 END MODULE TempData
 !********************************************************************
 MODULE CumData
   USE Typedef
   USE Paramdata, ONLY: maxelm, maxnod
   IMPLICIT NONE
   REAL(sp) :: cumCh0=0.,cumChr=0.,cumCh1=0.,CumRt=0.
   REAL(dp),DIMENSION(2):: CumQ=(/0.,0./),ChemS=(/0.,0./)
   REAL(dp) :: wCumA=0.,wCumT=0.,cCumA,cCumT,wVolI,cVolI
   REAL(sp) :: VolSink=0._dp, VolQ=0._dp, WatVolOld=0._dp,wBalR
   REAL(sp) :: ConAxx(maxelm),ConAyy(maxelm),ConAzz(maxelm)
   REAL(sp):: ConAxy(maxelm),ConAxz(maxelm),ConAyz(maxelm)
   INTEGER(ap) :: ListNE(maxnod)
   REAL(sp), ALLOCATABLE,DIMENSION (:) :: Qc, Q
   REAL(sp), ALLOCATABLE,DIMENSION (:) :: WatIn,SolIn
 END MODULE CumData
!********************************************************************
!> Data for soil domain
 MODULE DomData
   USE TypeDef
   IMPLICIT NONE
   REAL(sp):: xmin,xmax,ymin,ymax,zmin,zmax
 END MODULE DomData
 !********************************************************************

 MODULE GeoData
   USE TYPEDEF
   USE ParamData
   IMPLICIT NONE
   REAL(sp) :: geoaxs(maxemg),angaxs(maxemg,mxpnts)
   REAL(sp) :: tempax(maxemg,mxpnts)
   REAL(sp):: anglat(mxpnts),templt(mxpnts),geolat
   INTEGER(ap) ::nangax(maxemg),nanglt
 END MODULE GeoData
!********************************************************************
 MODULE MatData
   USE TypeDef
   USE GridData, ONLY :nband,nPt 
   IMPLICIT NONE
   INTEGER(ap) :: NumNZ, time_step=0
   INTEGER(ap), ALLOCATABLE, DIMENSION (:) ::jlu,IROW,JCOL
   REAL(dp), ALLOCATABLE, DIMENSION (:) :: alu
   REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: A
   REAL(dp), ALLOCATABLE, DIMENSION (:) :: A_sparse
   REAL(sp), ALLOCATABLE, DIMENSION (:,:) :: As
   REAL(sp), ALLOCATABLE, DIMENSION (:,:) :: A1
   REAL(dp), ALLOCATABLE, DIMENSION (:) :: B
   REAL(sp), ALLOCATABLE, DIMENSION (:) :: Bs
 END MODULE MatData
!********************************************************************
 !> output variables for observation probes
 MODULE ObsData
   USE TypeDef
   USE Paramdata, ONLY: maxobs
   ! npr: number of probes
   ! Pr: Probe ID
   ! Pt: plane direction (perp to X (1), to Y (2) or to Z (3)) or probe given by the user (4))
   ! CrP: crossing point or number of user nodes (if Pt=4)
   ! VarP: variable (WC=1, PH=2, both=3)
   ! distrP: 1=all, 2=average
   IMPLICIT NONE
   INTEGER(ap) :: npr, Pt(maxobs),nodebyPr(maxobs),nprof,Pr(maxobs)
   INTEGER(ap) :: NodePr(maxobs,1000) !1000=max number of node for a given plane
   INTEGER(ap) :: VarP(maxobs),distrP(maxobs),Pl(maxobs),varProf(maxobs)
   INTEGER(ap) ::CrP(maxobs)
   LOGICAL :: ObsOK,profOK
 END MODULE ObsData
 !********************************************************************
 !> Soil data
 MODULE SolData
   USE Typedef
   USE ParamData, ONLY: maxnod, maxIrrig, maxbdr, maxmat
   IMPLICIT NONE
   INTEGER(ap) :: NLevel,nMat,MatNum(maxnod)
   REAL(sp) :: ChPar(10,maxmat),par(11,maxmat),epsi,Peclet,Courant,ssMaxTab(maxmat)
   REAL(sp) :: PeCr
   REAL(sp), ALLOCATABLE,DIMENSION (:) :: Vx,Vy,Vz
   REAL(sp), ALLOCATABLE,DIMENSION (:) :: conO,con,cap
   REAL(sp), ALLOCATABLE,DIMENSION (:) :: Fc, Gc
   REAL(sp), ALLOCATABLE,DIMENSION (:) :: Dispxx,Dispyy,Dispzz
   REAL(sp), ALLOCATABLE,DIMENSION (:) :: Dispxy,Dispxz,Dispyz
   REAL(sp), ALLOCATABLE, DIMENSION(:) :: hOld,hTemp,hNew,Conc
   REAL(sp), ALLOCATABLE, DIMENSION(:) :: theta, theta_old
   REAL(dp), ALLOCATABLE, DIMENSION(:) :: concPar,ConcR
   REAL(dp), ALLOCATABLE, DIMENSION(:) :: massPar,massR
   REAL(dp):: sum_upt
   INTEGER(ap), ALLOCATABLE, DIMENSION(:)::Kode
   INTEGER(ap) ::KodCB(maxbdr+maxIrrig)
   LOGICAL :: soiltab=.false.,lTab=.false.,lMacro=.FALSE.
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: l_elmMacro
 END MODULE SolData
 !********************************************************************
 !> Rhizosphere Data
 MODULE RhizoData
   USE Typedef
   IMPLICIT NONE
   INTEGER(ap) :: rhizoModel
   REAL(dp) :: bulkPara(6), StaticRhizoPara(2) ! thtR, thtS, lambda, hcr, Ks, L;
! lambda_Static, hcr_Static 
   REAL(dp) :: RhizoPara(9) ! omega, beta, cTot, ni, d, tau0, gamma, thob, rhow
   REAL(dp) :: RhizoSpatialPara(2)
   REAL(sp), ALLOCATABLE, DIMENSION(:) :: thetaTot, thetaNonEq, thetaNonEqOld
   REAL(sp), ALLOCATABLE, DIMENSION(:) :: cTot_r, Rnorm, tauTht, hEqRhizo 
   LOGICAL :: lRhizo
 END MODULE RhizoData
 !********************************************************************
 !> Data for root growth components influenced by soil strength
 MODULE StrData
 USE TypeDef
 USE ParamData, ONLY: maxnod
 IMPLICIT NONE
 REAL (sp):: s(maxnod),condu(maxnod),imps(maxnod),ssmax,simp,refgrd
 END MODULE StrData
 !********************************************************************
 MODULE WatFun
   USE Typedef
   USE SolData, ONLY :nmat
   USE RhizoData
   IMPLICIT NONE
   
   INTEGER::nTab
   REAL(sp),ALLOCATABLE, DIMENSION (:) :: hTab,hTab_MFP
   REAL(sp), ALLOCATABLE, DIMENSION (:,:) :: TheTab,ConTab,CapTab,MFPTab
   REAL(sp) ::alh1,dlh
 CONTAINS

	! initialize matrix sizes for tabulated soil characteristics
	SUBROUTINE IniTab
	 IMPLICIT NONE
	 ALLOCATE (hTab(1:nTab))
	 ALLOCATE (TheTab(1:nTab,1:nmat))
	 ALLOCATE (ConTab(1:nTab,1:nmat))
	 ALLOCATE (CapTab(1:nTab,1:nmat))
	 ALLOCATE (MFPTab(1501,1:nmat))
	 ALLOCATE (hTab_MFP(1501))
	END SUBROUTINE IniTab
		  
	! current local theta (volumetric soil water content) values:
	REAL(sp) FUNCTION Fth(h,Parloc, SoilNode)
	  REAL(dp) :: thr,ths,THT, hcr
	  REAL(sp),INTENT(in) :: h,Parloc(:)
	  INTEGER :: SoilNode
      ! if rhizosphere then thr and ths are from rhizobulk
      IF (.NOT. lRhizo) THEN
    	  thr=Parloc(2)
	      ths=Parloc(3)
          IF (h.LT.0.0_dp) THEN
		    THT=FThNrm(h,Parloc)
            Fth=thr+THT*(ths-thr)
          ELSE
		    Fth=ths
          ENDIF
      ! In case Rhizosphere is considered   
      ELSE              
         hcr = bulkPara(4)
         thr=bulkPara(1)
         ths=bulkPara(2)
         IF (h .LT. hcr) THEN
            THT = FthNrm(h,Parloc)
            Fth = thr + THT*(ths-thr)
         ELSE
            Fth = ths
         ENDIF
      ENDIF
	  RETURN
	END FUNCTION Fth
	
    ! current local THETA (0...1 normalized theta) values:
	REAL(sp) FUNCTION FThNrm(h,Parloc)
	 REAL(dp) :: a,n,m,a2,n2, m2,w1,w2
	 REAL(dp) :: lambda, hcr
	 REAL(sp),INTENT(in) :: h,Parloc(:)
	 INTEGER :: MatMod
	 a=Parloc(4)
	 n=Parloc(5)
	 m=1._dp-1._dp/n
	 MatMod=INT(Parloc(1))
	 IF (lRhizo) THEN ! Rhizosphere is considered
		lambda = bulkPara(3)
		hcr    = bulkPara(4)
                MatMod = 3          
	ENDIF	
	 SELECT CASE(MatMod)
	 CASE (1) ! case1: Mualem Van Genuchten
		IF (h.LT.0.0_dp) THEN
		   FthNrm=(1._dp+(a*ABS(h))**n)**(-m)
		ELSE
		   FthNrm=1.0_dp
		ENDIF
	 CASE (2) ! case2: Dual porosity model
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
     CASE (3) ! Case 3, Brooks and Corey as define in Rhizo.in
        IF (h .LT. hcr) THEN
            FthNrm = (hcr/h)**lambda
        ELSE
            FthNrm = 1.0_dp
        ENDIF
	 CASE DEFAULT
		STOP 'Sorry, only three soil models possible at the moment. Please check <soil.in>.'   
	 END SELECT
	 RETURN
	END FUNCTION FThNrm

! current local THETA values (from soiltab):
REAL(sp) FUNCTION Fth_soiltab(h,M)
 REAL(sp),INTENT(in) :: h
 INTEGER(sp),INTENT(in) :: M
 INTEGER(sp) :: iT
 IF (h.LT.0.0_dp) THEN
	iT=1
	DO WHILE ((h.LT.hTab(iT+4)).AND.(iT.LE.nTab-8))
	   iT=iT+4
	END DO
	DO WHILE ((h.LT.hTab(iT+1)).AND.(iT.LE.nTab-2))
	   iT=iT+1
	END DO
	Fth_soiltab=TheTab(iT,M)+(TheTab(iT+1,M)-TheTab(iT,M))/(hTab(iT+1)-hTab(iT))*(h-hTab(iT))
 ELSE
	Fth_soiltab=TheTab(1,M)
 ENDIF
 RETURN
END FUNCTION Fth_soiltab
!*********************************************************************
REAL(sp) FUNCTION Fh_from_Th(Th,Parloc)
 REAL(dp) :: thr,ths,a,n,m
 REAL(sp),INTENT(in) :: Th,Parloc(:)
 INTEGER :: MatMod
 
 thr=Parloc(2)
 ths=Parloc(3)
 a=Parloc(4)
 n=Parloc(5)
 m=1._dp-1._dp/n
 MatMod=INT(Parloc(1))
 
 SELECT CASE(MatMod)
 CASE (1) !> + case1: Mualem Van Genuchten
	IF (Th.LT.ths) THEN
	   Fh_from_Th=-(((((Th-thr)/(ths-thr))**(-1.0_dp/m))-1.0_dp)**(1.0_dp/n))/a
	ELSE
	   Fh_from_Th=0.0_dp
	ENDIF
 CASE (2)  !> + case2: Dual porosity model
	STOP 'Can t currently use Fh_from_Th with dual porosity model'
 END SELECT
 RETURN
END FUNCTION Fh_from_Th
!*********************************************************************
!> current local matric flux potential value
REAL(sp) FUNCTION Fmfp_soiltab(h,M,Parloc)
 REAL(sp),INTENT(in) :: h,Parloc(:)
 INTEGER(sp),INTENT(in) :: M
 INTEGER(sp) :: iT
 IF (h.LT.0.0_dp) THEN
	iT=1
	DO WHILE ((h.LT.hTab_MFP(iT+100)).AND.(iT.LE.1500-200))
	   iT=iT+100
	END DO
	DO WHILE ((h.LT.hTab_MFP(iT+10)).AND.(iT.LE.1500-20))
	   iT=iT+10
	END DO
	DO WHILE ((h.LT.hTab_MFP(iT+1)).AND.(iT.LE.1500-2))
	   iT=iT+1
	END DO
	Fmfp_soiltab=MFPTab(iT,M)+(MFPTab(iT+1,M)-MFPTab(iT,M))/(hTab_MFP(iT+1)-hTab_MFP(iT))*(h-hTab_MFP(iT))
 ELSE
	Fmfp_soiltab=MFPTab(1,M)+h*Parloc(6)
 ENDIF
 RETURN
END FUNCTION Fmfp_soiltab
!*********************************************************************
REAL(sp) FUNCTION Fh_from_mfp_soiltab(mfp,M,Parloc)
 REAL(sp),INTENT(in) :: mfp,Parloc(:)
 INTEGER(sp),INTENT(in) :: M
 INTEGER(sp) :: iT
 IF (mfp.GT.0.0_dp) THEN
	iT=1
	DO WHILE ((mfp.LT.MFPTab(iT+100,M)).AND.(iT.LE.1500-200))
	   iT=iT+100
	END DO
	DO WHILE ((mfp.LT.MFPTab(iT+10,M)).AND.(iT.LE.1500-20))
	   iT=iT+10
	END DO
	DO WHILE ((mfp.LT.MFPTab(iT+1,M)).AND.(iT.LE.1500-2))
	   iT=iT+1
	END DO
	Fh_from_mfp_soiltab=hTab_MFP(iT)+(hTab_MFP(iT+1)-hTab_MFP(iT))/(MFPTab(iT+1,M)-MFPTab(iT,M))*(mfp-MFPTab(iT,M))
 ELSE
	Fh_from_mfp_soiltab=-15010.0_sp
 ENDIF
 RETURN
END FUNCTION Fh_from_mfp_soiltab
!********************************************************************
!> current local soil conductivity values   
REAL(sp) FUNCTION FKP(h,Parloc,soilNode)
 REAL(dp) ::  Ks,THT,lambda,a,a2,n,n2,m,m2,w1,w2,Sv1,Sv2,rNumer,rDenom
 REAL(sp),INTENT(in) :: h,Parloc(:)
 REAL(dp) :: L, ni,d, thtR,thtS, rhob, rhow,cw,  hcr
 REAL(sp) :: Kb,S, thtM, mu, thtBulk
 INTEGER :: MatMod, soilNode
 LOGICAL :: Ralloc =.TRUE.
 
 a=Parloc(4)
 n=Parloc(5)
 Ks=Parloc(6)
 lambda=Parloc(7)
 m=1._dp-1._dp/n
 MatMod=INT(Parloc(1))
 !> FKP in case of rhizosphere model
 IF (lRhizo) MatMod = 3
 SELECT CASE(MatMod)
 CASE (1) ! Mualem Van Genuchten
	IF (h.LT.0.0_dp) THEN
	   THT=(1._dp+(a*ABS(h))**n)**(-m)
	   FKP=Ks*(THT**lambda)*(1._dp-(1._dp-THT**(1._dp/m))**m)**2
	ELSE
	   FKP=Ks
	ENDIF
 CASE (2) !Dual porosity model
	w2=parloc(8)
	a2=parloc(9)
	n2=parloc(10)
	w1=1.0-w2
	m2=1._dp-1._dp/n2
	IF (h.LT.0.0_dp) THEN
	   Sv1=(1._dp+(a*ABS(h))**n)**(-m)
	   Sv2=(1._dp+(a2*ABS(h))**n2)**(-m2)
	   THT=w1*Sv1+w2*Sv2
	   rNumer=w1*a*(1._dp-(1._dp-Sv1**(1._dp/m))**m) +&
			w2*a2*(1._dp-(1._dp-Sv2**(1._dp/m2))**m2) 
	   rDenom=w1*a+w2*a2
	   FKP=Ks*(THT**lambda)*(rNumer/rDenom)**2
	ELSE
	   FKP=Ks
	ENDIF
  CASE (3) ! Condectivity for Rhizosphere model
        thtR =   bulkPara(1)
        thtS =   bulkPara(2)
        lambda = bulkPara(3)
        hcr =    bulkPara(4)
        Ks =     bulkPara(5)
        L  =     bulkPara(6)
        ni =     rhizoPara(4)
        d  =     rhizoPara(5)
        rhob =   rhizoPara(8)
        rhow =   rhizoPara(9)
        ! Total Saturation
        S = (thetaTot(soilNode)-thtR)/(thtS-thtR)
        ! Bulk conductivity
        IF (S .LT. 1.) THEN
            Kb = Ks*S**(2./lambda + 2. + L)
            ! mucilage water content thtM
            IF (.NOT. ALLOCATED(Rnorm)) THEN
                Ralloc = .false.
                ALLOCATE (Rnorm(soilNode))
                Rnorm = 1
            ENDIF 
            IF (Rnorm(soilNode) .LT. 1.) THEN
                IF (h .LT. hcr) THEN
                        thtBulk = (thtS-thtR)*(hcr/h)**lambda + thtR
                ELSE
                        thtBulk = thtS
                ENDIF
                thtM = ((thetaTot(soilNode)-thtBulk*Rnorm(soilNode)))/&
                        (1. - Rnorm(soilNode))
                cw = (cTot_r(soilNode)*rhob)/(rhow*thtM)
                mu = 1. + ni*cw**d
                FKP = 1./mu * Kb
            ELSE
                FKP = Kb
            ENDIF
        ELSE
           FKP = Ks
        ENDIF
 END SELECT
 IF (.NOT. Ralloc) THEN
        DEALLOCATE(Rnorm)
        Ralloc = .TRUE.
ENDIF
 RETURN
END FUNCTION FKP

! current local conductivity values (from soiltab):
REAL(sp) FUNCTION FKP_soiltab(h,M)
 REAL(sp),INTENT(in) :: h
 INTEGER(sp),INTENT(in) :: M
 INTEGER(sp) :: iT
 IF (h.LT.0.0_dp) THEN
	iT=1
	DO WHILE ((h.LT.hTab(iT+4)).AND.(iT.LE.nTab-8))
	   iT=iT+4
	END DO
	DO WHILE ((h.LT.hTab(iT+1)).AND.(iT.LE.nTab-2))
	   iT=iT+1
	END DO
	FKP_soiltab=(ConTab(iT,M)+(ConTab(iT+1,M)-ConTab(iT,M))/(hTab(iT+1)-hTab(iT))*(h-hTab(iT)))
 ELSE
	FKP_soiltab=ConTab(1,M)
 ENDIF
 RETURN
END FUNCTION FKP_soiltab

!********************************************************************
!> current local soil water capacity values
REAL(sp) FUNCTION FCP(h,Parloc)
 REAL(dp) :: thr,ths,a,n,m,C2a,C2b,m2,w1,W2,a2,n2
 INTEGER :: MatMod
 REAL(sp), INTENT(in) :: h,Parloc(:)
 REAL(dp) :: hcr, lambda 
 MatMod=INT(Parloc(1))
 thr=parloc(2)
 ths=parloc(3)
 a=parloc(4)
 n=parloc(5)
 m=1._dp-1._dp/n
 ! Rhizosphere
 IF (lRhizo) MatMod = 3
 SELECT CASE(MatMod)
 CASE (1) ! Mualem Van Genuchten
	IF (h.LT.0.0_dp) THEN
	   FCP=(ths-thr) *(a*n*m*((a*ABS(h))**(n-1._dp)))/((1._dp+(a*ABS(h))**n)**(m+1._dp))
	ELSE
	   FCP=0.0_dp
	ENDIF
 CASE (2) !Dual porosity model
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
 CASE (3) ! Rhizosphere model
        thr =    bulkPara(1)
        ths =    bulkPara(2)
        lambda = bulkPara(3)
        hcr =    bulkPara(4)
        IF (h .LT. hcr) THEN
                FCP = -(ths-thr)*lambda*(hcr/h)**lambda/h
        ELSE
                FCP = 0.0_dp
        ENDIF
         
 END SELECT
 RETURN
END FUNCTION FCP

! current local soil water capacity values (from soiltab):
REAL(sp) FUNCTION FCP_soiltab(h,M)
 REAL(sp),INTENT(in) :: h
 INTEGER(sp),INTENT(in) :: M
 INTEGER(sp) :: iT
 IF (h.LT.0.0_dp) THEN
	iT=1
	DO WHILE ((h.LT.hTab(iT+4)).AND.(iT.LE.nTab-8))
	   iT=iT+4
	END DO
	DO WHILE ((h.LT.hTab(iT+1)).AND.(iT.LE.nTab-2))
	   iT=iT+1
	END DO
	FCP_soiltab=(CapTab(iT,M)+(CapTab(iT+1,M)-CapTab(iT,M))/(hTab(iT+1)-hTab(iT))*(h-hTab(iT)))
 ELSE
	FCP_soiltab=0.0_dp
 ENDIF
 RETURN
END FUNCTION FCP_soiltab
END MODULE WatFun
 !********************************************************************
  MODULE disToRoot
  USE typedef
  USE RootData, ONLY: xs, ys, zs, nrec, timorg
  USE GridData, ONLY: xgrid, ygrid, zgrid, nPt, dxGrid, dyGrid, dzGrid
  USE RhizoData, ONLY: RhizoPara, RhizoSpatialPara, cTot_r, Rnorm
  USE DoussanMat, ONLY: loc_Q
  REAL(dp), ALLOCATABLE, DIMENSION (:,:) ::  rhizoMatrix

  CONTAINS
!> RStat is a Subroutine that calculate Rnorm and Ctot_r for a static root
!system
    SUBROUTINE RStat(t)
    IMPLICIT NONE
    REAL(sp) ::t
    REAL(sp) :: alpha,  r, dist, cTot,  RnormEXP
    INTEGER  ::  i, k, jj, count1
    INTEGER, allocatable, dimension(:) :: rhizoNodes
    ALLOCATE(cTot_r(nPt))
    ALLOCATE(Rnorm(nPt))
    ALLOCATE(rhizoNodes(nPt))
    !> Initalized parameters
    rhizoNodes = 0  
    Rnorm(1:nPt) = 1.
    !> Finding and storing soil nodes near the root
    DO i=0,size(loc_Q(:,1,1,1))-1
            DO jj=1,size(loc_Q(1,1,:,1))
                    DO k=1,8
                            IF(loc_Q(i,k,jj,1) .GT. 0) THEN
                                    rhizoNodes(loc_Q(i,k,jj,1)) = loc_Q(i,k,jj,1) ! Rhizo node true
                            ENDIF
                    ENDDO
            ENDDO
    ENDDO
    !> RhizoMatrix contain information of the soil node (:,1), min distance from the root (:,2), and the age of the closest root (:,3)
     count1=0
     DO i=1,size(rhizoNodes)
             IF(rhizoNodes(i) .GT. 0) count1 = count1+1
     ENDDO
     ALLOCATE(rhizoMatrix(count1,3))
     ! Initialized rhizoMatrix !
     rhizoMatrix(:,1) = 0
     rhizoMatrix(:,2) = 10000.
     rhizoMatrix(:,3) = 10000.
     k=1
     DO i=1,npt
       IF(rhizoNodes(i) .GT. 0) THEN
               rhizoMatrix(k,1) = rhizoNodes(i)
               DO jj=1,nrec
                       dist = SQRT((xs(jj)-xgrid(i))**2 +(ys(jj)-ygrid(i))**2 + (zs(jj)-zgrid(i))**2)
                       IF (dist .LT. rhizoMatrix(k,2)) THEN
                               rhizoMatrix(k,2) = dist
                               rhizoMatrix(k,3) = timorg(jj)
                       ENDIF
               ENDDO
               k = k +1
       ENDIF
     ENDDO

    cTot  = RhizoPara(3)
    alpha = RhizoSpatialPara(1)
    RnormEXP = RhizoSpatialPara(2)
    cTot_r = 0.
    !> CTot_r is a vector that contain the concentration of mucilage at each soil
    !node. It is calculated based on the Monod equation and has the Form:
    !Ctot_r = ctot * (timeOrigine/(maxTime/alpha + timeorigine)). 
    !alpha is an empirical parameter. The bigger alpha the concentration
    !approach cTot
    !faster. 

    DO i=1,count1
        cTot_r(INT(rhizoMatrix(i,1))) =cTot*rhizoMatrix(i,3)/(MAXVAL(rhizoMatrix(:,3))/alpha + rhizoMatrix(i,3))
        Rnorm(INT(rhizoMatrix(i,1))) = EXP(-RnormEXP*cTot_r(INT(rhizoMatrix(i,1))))
    ENDDO

    call GrowRhizo(t)
    OPEN (UNIT=15,FILE='out/Rnorm.out',STATUS='UNKNOWN')
    WRITE(15,*) 'soilNode   minDist    age cTot_r, Rnorm'
    DO i=1,size(rhizoMatrix(:,1))
            WRITE(15,*) rhizoMatrix(i,1), rhizoMatrix(i,2), rhizoMatrix(i,3),cTot_r(INT(rhizoMatrix(i,1))),Rnorm(INT(rhizoMatrix(i,1)))
    END DO
    CLOSE(15)

  END SUBROUTINE RStat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!** A subroutine to calculate the concentration and Rnorm for growing root. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  SUBROUTINE GrowRhizo(t)
    REAL(sp) :: t, relConc, timeOrigin
    INTEGER  :: i
    REAL(sp) :: alpha, cTot, RnormEXP

    cTot = RhizoPara(3)
    alpha = RhizoSpatialPara(1)
    RnormEXP = RhizoSpatialPara(2)

    DO i=1,size((rhizoMatrix(:,3)))
            timeOrigin = rhizoMatrix(i,3)
            relConc = relC(t, timeOrigin,alpha)
            cTot_r(INT(rhizoMatrix(i,1))) = relConc * cTot !*rhizoMatrix(i,3)/(MAXVAL(rhizoMatrix(:,3))/alpha + rhizoMatrix(i,3))
            Rnorm(INT(rhizoMatrix(i,1))) = EXP(-RnormEXP*cTot_r(INT(rhizoMatrix(i,1))))
   ENDDO

  END SUBROUTINE GrowRhizo

  FUNCTION relC(t,tOrig,alpha)
  ! A function that calculate the relative concentration of mucilage. 
  ! The function get the current simulation time and the time of origin
  ! and return the relative concentration fot that time. 

        REAL(sp) relC ! Relative concentration of mucilage
        REAL(sp) tDec, t, tOrig, alpha ! Time from which concentration start to decay
        tDec = tOrig + 5.0

        IF (t .LE. tOrig) THEN
                relC = 0.0
        ELSEIF (t .GT. tOrig .AND. t .LT. tDec) THEN
                relC = erf(t)
        ELSE
                relC = erf(t) -erf((t-tDec)/alpha)
        ENDIF

        RETURN
  END FUNCTION relC

END MODULE disToRoot
!********************************************************************
!> A Module to calculate the rhizopshrer static parameters. This is also used
!for the rhizosphere initial condition. 
MODULE RhizoStat
  USE typedef
  USE RhizoData, ONLY: bulkPara, RhizoPara

  CONTAINS
  !< Secant Methodi
  !*************************************************
  FUNCTION ERR(tht1, h1)
      REAL(dp) :: omega, rhob, rhow, cTot, beta, hcr, lambda
      REAL(sp) :: tht1, heq, h1, thtModel, ERR
      omega = RhizoPara(1); beta = RhizoPara(2); cTot = RhizoPara(3)
      rhob = RhizoPara(8); rhow = RhizoPara(9)
      lambda = bulkPara(3); hcr = bulkPara(4)
      heq = h1 + omega * (rhob*cTot/(rhow*tht1))**beta
      thtModel = RetCurve(heq,hcr,lambda)
      ERR = thtModel - tht1
      RETURN
  END FUNCTION ERR

  !> RetCurve return theta based on brooks and corey
  FUNCTION RetCurve(h, hcr, lambda)
    REAL(sp) :: h, S, RetCurve
    REAL(dp) :: hcr, lambda, thtR, thtS
    thtR = bulkPara(1); thtS = bulkPara(2)
    IF (h .LT. hcr) THEN
        S = (hcr/h)**lambda
    ELSE
        S = 1.
    ENDIF
    RetCurve = S*(thtS-thtR) + thtR
    RETURN
  END FUNCTION RetCurve
  !*************************************************
      
END MODULE RhizoStat
