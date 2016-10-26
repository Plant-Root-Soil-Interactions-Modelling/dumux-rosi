!===============================================================================
! Source file ROOT GROWTH with RootTyp |||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
   SUBROUTINE RunRootTip(t)
! run RootTip for a given time (in days) from 0
! return naxes:= primary roots in RooTyp
! call C-functions, which return arrays of variables
! most of these are then saved under RootData module
   USE Typedef
   USE RootData
   USE GridData, only: continu
   IMPLICIT NONE
   Real(dp) :: origin(3)=(/ 0.0,0.0,0.0 /)
   Integer(sp):: n,simtime, n_nodes, n_meris,n_br_crea,n_br_del, n_axes,t
   Integer(sp) ::maxlist,n_nodecorr,i,nn,ipl=1!runroottip currently doesn´t work with more than one plant. It should be updated... (Couvreur dec 2009)
   Integer, allocatable :: node(:)
   Integer, allocatable :: prev_node(:)
   Integer, allocatable :: ordn(:)
   Integer, allocatable :: orda(:)
   Integer, allocatable :: agerec(:)
   Integer, allocatable :: axeF(:)
   Integer, allocatable :: prev_nodea(:)
   Integer, allocatable :: axenf(:)
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
   Allocate (axenf(n_nodes))
   Allocate (axeF(n_meris))
   Allocate (orda(n_meris))
   Allocate (xa(n_meris))
   Allocate (ya(n_meris))
   Allocate (za(n_meris))
   Allocate (diama(n_meris))
   Allocate (agea(n_meris))
   Allocate (prev_nodea(n_meris))

! extract nodes C-function < interface
   CALL extract_nodes(node,prev_node,ordn,xrec,yrec,zrec,diamrec,agerec,axeF,orda,xa,&
 & ya,za,diama,agea,prev_nodea,n_axes,axenf,n_nodecorr) 
   nrec=n_nodecorr

!meaning of node is not clear?
   xs(1:nrec)=xrec(1:nrec)
   ys(1:nrec)=yrec(1:nrec)
   zs(1:nrec)=-zrec(1:nrec)!roottip in mm
   !Z must decrease downward!!
   do n=2,nrec
      if (zs(n)>0.0_sp)	 zs(n)=-zs(n) !in case there is roots >soil surface
      if ((n.NE.1).AND.(axenf(n).NE.(axenf(n-1)))) then
	  if (axenf(n).NE.(axenf(n-1)+1)) then
	     print *,'axe',axenf(n-1)+1,'is missing'
	  endif
      endif
   enddo

! create R-SWMS variables for root
   irecpr(1:nrec)=prev_node(1:nrec)
   ordseg(1:nrec)=ordn(1:nrec)+1. !in Roottip from 0 to 7, here from 1 to 8
   ibrseg(1:nrec)=axenf(1:nrec)
   timorg(1:nrec)=agerec(1:nrec) !day of creation+time real=age
   segdiam(1:nrec)=diamrec(1:nrec)

!growing apices= "meristems with axes" in RootTyp
   ngrow=n_axes
   nbr=n_axes !number of apex=number of branhces
   xg(1:n_axes)=xa(1:n_axes)
   yg(1:n_axes)=ya(1:n_axes)
   zg(1:n_axes)=-za(1:n_axes)
   irecsg(1:n_axes)=prev_nodea(1:n_axes)
   ordgrw(1:n_axes)=orda(1:n_axes)+1
   ibrgrw(1:n_axes)=axeF(1:n_axes)
   Allocate(connex(1:nrec))

!check root nodes outside of the soil
   if (zg(1)>0.0_sp) zg(1)=-zg(1)
   do n=2,ngrow
      if (zg(n)>0.0_sp) zg(n)=-zg(n) !in case there is roots >soil surface
      if (axenf(n).NE.(axenf(n-1))) then
	  if (axenf(n).NE.(axenf(n-1)+1)) then
	     print *,'axe',axenf(n-1)+1,'is missing'
	  endif
      endif
   enddo

!check root type (Couvreur nov 2010)
   IF (maxval(ordseg).EQ.13) THEN		!Recognizes old roots created by RootTyp
       maizeroottyp=.TRUE.
   ELSEIF (maxval(ordseg).EQ.5) THEN
       loliumroottyp=.TRUE.
   ELSEIF (maxval(ordseg).EQ.20) THEN
       wheatroottyp=.TRUE.
   ENDIF

!estimate length and surface of segments and get naxes
   CALL Estimate_Seglen
!   CALL OutRoo(real(0),naxes,0,0,0,0,0)
!simplifiy system
   DO n=1,2
      CALL SimplifyRoot
      CALL Estimate_Seglen2
   ENDDO
!   CALL OutRoo(real(999),naxes,0,0,0,0,0)
!   CALL AdaptOutRoot
!   CALL Estimate_Seglen2
   CALL close_c_interface() !C-function
   CALL finish_RootType1() !C-function
   treal=t
   print *,'number of nodes after',n,' iterations at time t= ', treal,' is ',nrec
   CALL OutRoo(treal,0,0,0,0,0)
   IF (.NOT.continu) CALL CheckSize(1)	!No need to check size if continuous domain. ipl is 1 because RootTyp can currently only simulate 1 root system (Couvreur feb 2010)
   END SUBROUTINE RunRootTip
!********************************************************************************
   SUBROUTINE Estimate_Seglen
   USE typedef
   USE ParamData
   USE RootData
   IMPLICIT NONE
   REAL(sp) :: xend,yend,zend,xtop,ytop,ztop,xtop2,ytop2,ztop2
   INTEGER(sp) :: inode,n,brn,inode2,num
   if(lrrt) connex(:)=.FALSE.
   inode=1
   naxes=0
   DO n=1,nbr !for all branches
      brlgth(n)=0.0_dp
      inode=irecsg(n) !inode =prec
      xtop=xs(inode)!end note=apex
      ytop=ys(inode)
      ztop=zs(inode)
      xend=xg(n)	!Length associated with a node is the length of the segment under the node (Couvreur feb 2010)
      yend=yg(n)
      zend=zg(n)
      seglen(inode)=SQRT((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
      if (seglen(inode)==0) then	!apical node "on" the meristem
         xs(inode)=xs(irecpr(inode))+(xtop-xs(irecpr(inode)))/2.0_dp
         ys(inode)=ys(irecpr(inode))+(ytop-ys(irecpr(inode)))/2.0_dp
         zs(inode)=zs(irecpr(inode))+(ztop-zs(irecpr(inode)))/2.0_dp
         xtop=xs(inode)
         ytop=ys(inode)
         ztop=zs(inode)
         seglen(inode)=SQRT((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
      endif
      brlgth(n)=brlgth(n)+seglen(inode)
      num=1
      IF (lrrt) then
         segsur(inode)=seglen(inode)*pi*segdiam(inode)!segment surface
      ELSE 
         segsur(inode)=seglen(inode)*pi*2._dp*segrad(inode) 
      END IF
      DO WHILE (ibrseg(irecpr(inode))==n.AND.irecpr(inode).GT.0)
         inode2=irecpr(inode)
         xend=xtop
         yend=ytop
         zend=ztop
         xtop=xs(inode2)
         ytop=ys(inode2)
         ztop=zs(inode2)
         seglen(inode2)=SQRT((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
         if (seglen(inode2)==0) then	!superposed nodes
            xs(inode2)=xs(irecpr(inode2))+(xtop-xs(irecpr(inode2)))/2.0_dp
            ys(inode2)=ys(irecpr(inode2))+(ytop-ys(irecpr(inode2)))/2.0_dp
            zs(inode2)=zs(irecpr(inode2))+(ztop-zs(irecpr(inode2)))/2.0_dp
            xtop=xs(inode2)
            ytop=ys(inode2)
            ztop=zs(inode2)
            seglen(inode2)=SQRT((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
         endif
         IF (lrrt) THEN
            segsur(inode2)=seglen(inode2)*pi*segdiam(inode2)!segment surface
         ELSE
            segsur(inode2)=seglen(inode2)*pi*2._dp*segrad(inode2)
         END IF
         brlgth(n)=brlgth(n)+seglen(inode2)
         num=num+1;
         inode=inode2
      ENDDO
      br_rec(n)=irecpr(inode)	!records at which node (on the father axis) this axis is connected
      num_seg(n)=num
      IF (maizeroottyp) THEN
         IF (ordgrw(n).LE.11) naxes=naxes+1 !number of axes ("principal roots")
      ELSEIF (loliumroottyp) THEN
         IF (ordgrw(n).LT.3) naxes=naxes+1
      ELSEIF (wheatroottyp) THEN
         IF (ordgrw(n).LE.18) naxes=naxes+1
      ELSEIF (ordgrw(n).LT.3) THEN					!Other roots (Couvreur nov 2010)
         naxes=naxes+1 !number of axes ("principal roots")
      ENDIF
      if (lrrt .and. (br_rec(n).NE.0)) connex(br_rec(n))=.TRUE. !logical which defines whether there is a connection to that node
   ENDDO
   RETURN
   END SUBROUTINE Estimate_Seglen
!****************************************************************************
   SUBROUTINE SimplifyRoot	!Redesigned for the new RootTyp with senescence (Couvreur may 2010)
   USE typedef
   USE RootData, only : irecpr,nrec,ibrseg,irecsg,connex,xs,ys,zs,timorg,segdiam,ordseg,seglen,xg,yg,zg
   USE GridData, only : dxGrid,dyGrid,dzGrid
   IMPLICIT NONE
   real(sp) :: xrec(1:nrec),yrec(1:nrec),zrec(1:nrec),agerec(1:nrec)
   real(dp) :: diamrec(1:nrec),ltot
   INTEGER(sp) :: ordn(1:nrec),axenf(1:nrec),resolution
   INTEGER(sp) :: ibr,inew,iold,i_connex,prec(nrec),iprec_old,n,oldi(nrec)
   INTEGER(sp) :: old2new(nrec),i,newtot,ibrprec,iprec(1:nrec)

!initialisation
   xrec=xs(1:nrec)
   yrec=ys(1:nrec)
   zrec=zs(1:nrec)
   axenf=ibrseg(1:nrec)
   agerec=timorg(1:nrec)
   diamrec=segdiam(1:nrec)
   ordn=ordseg(1:nrec)
   iprec=irecpr(1:nrec)
   resolution=min(2.,dxGrid,dyGrid,dzGrid)			!Adapt root simplification to grid resolution (threshold = 2 cm because root properties evolve with distance of that order) (Couvreur mar 2010)
   inew=1
   iold=1!start with node 2==proximal node (Couvreur may 2010)
   old2new(1)=1
!check each root
   DO WHILE (iold<nrec)
      inew=inew+1
      iold=iold+1
      iprec_old=iprec(iold)!node before
      irecpr(inew)=old2new(iprec(iold))
      ltot=seglen(iprec_old) !length
      DO WHILE ((axenf(iprec_old)==axenf(iold)).AND.(connex(iold).EQV..FALSE.)&
	        &.AND.(seglen(iold)<0.8*resolution).AND.(seglen(iold)+ltot.LE.resolution)&
		 &.AND.((iold).NE.irecsg(ibrseg(iold)))) !this is not the apical node
         !same branches + no connection to iold
         !(one more condition : the new segment is shorter than the grid resolution) (Couvreur feb 2010)
         old2new(iold)=999999
         ltot=ltot+seglen(iold)
         iold=iold+1
	  iprec_old=iprec(iold)!node before
      ENDDO
      ibrseg(inew)=axenf(iold) !br#
      xs(inew)=xrec(iold)
      ys(inew)=yrec(iold)
      zs(inew)=zrec(iold)
      ordseg(inew)=ordn(iold)
      timorg(inew)=agerec(iold) !orig time
      segdiam(inew)=diamrec(iold) !diam
      oldi(inew)=iold !keep in mind the previous numerotation
      old2new(iold)=inew
      IF (irecsg(ibrseg(inew))==oldi(inew)) THEN
! if the previous node of the apex of the branch where node inew is is inew, then itmust be also updated
         irecsg(ibrseg(inew))=inew
      ENDIF
   ENDDO
   nrec=inew
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
   USE RootData, only : seglen,segsur,brlgth,segdiam,irecpr,nbr,ordgrw&
   &,ibrseg,zs,ys,xs,irecsg,num_seg,br_rec,connex,zg,yg,xg,naxes
   IMPLICIT NONE
   REAL(sp) :: xend,yend,zend,xtop,ytop,ztop,xtop2,ytop2,ztop2
   INTEGER(sp) :: inode,n,brn,inode2,num,connected2
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
      DO WHILE (brn==n.AND.inode.GT.0)
	  xtop=xs(inode)
	  ytop=ys(inode)
         ztop=zs(inode)
         seglen(inode)=SQRT((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
	  !length of node is the length of the segment located on the top of it	!Wrong, at the basis (Couvreur feb 2010)
	  !length of node 1=0
         if (seglen(inode)==0) then	!superposed nodes
            xs(inode)=xs(irecpr(inode))+(xtop-xs(irecpr(inode)))/2.0_dp
            ys(inode)=ys(irecpr(inode))+(ytop-ys(irecpr(inode)))/2.0_dp
            zs(inode)=zs(irecpr(inode))+(ztop-zs(irecpr(inode)))/2.0_dp
            xtop=xs(inode)
            ytop=ys(inode)
            ztop=zs(inode)
	     seglen(inode)=SQRT((xtop-xend)**2+(ztop-zend)**2+(ytop-yend)**2)
         endif
	  segsur(inode)=seglen(inode)*pi*segdiam(inode)!segment surface
         brlgth(brn)=brlgth(brn)+seglen(inode)
         num=num+1
	  xend=xtop
	  yend=ytop
	  zend=ztop
         inode=irecpr(inode)
	  connected2=inode
	  brn=ibrseg(inode)
      ENDDO
      num_seg(n)=num
      IF (connected2.NE.0) connex(connected2)=.TRUE.
    ENDDO
    END SUBROUTINE Estimate_Seglen2
!****************************************************************************
   SUBROUTINE CheckSize(ipl)		!Checksize has to know wich plant is being checked in order to place it correctly with report to the grid (Couvreur dec 2009)
! check max/min position of roots as compared to the soil grid
   USE RootData
   USE GridData
   IMPLICIT NONE
   REAL(dp)::maxX,maxY,maxZ,maxXs,maxYs,maxZs,minX,minY,minZ,minXs,minYs,minZs
   INTEGER(sp)::ipl
   maxX=minval(xGrid)+nex*dxgrid		!Adapted to continuous and non continuous soil domain (Couvreur dec 2009)
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
   if (maxXs>maxX) THEN
      print *,'X root too large'
      goto 20
   endif
   if (maxYs>maxY) THEN
      print *,'Y root too large'
      goto 20
    endif
   if (maxZs>maxZ) THEN !
      print *,'Upper root node (=',maxZs,') is higher than soil max. z (=',maxZ,')'
      GOTO 20
   endif
   if (minXs<minX) THEN
      print *,'X root too small'
      goto 20
   endif
   if (minYs<minY) THEN
      print *,'Y root too small'
      goto 20
   endif
   if (minZs<minZ) THEN
      print *,'Lower root node (=',minZs,') is deeper than soil min. z (=',minZ,')'
      goto 20
   endif
    RETURN
    20 print *,'root max:',maxXs,maxYs,maxZs,'root min:',minXs,minYs,minZs,'soil:',maxX,maxY
	print *,'Please re-run R-SWMS/RootTyp'
	STOP
   END SUBROUTINE CheckSize
!****************************************************************************
