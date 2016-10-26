!==============================================================================
! Source file ROOT GROWTH ||||||||||||||||||||||||||||||||||||||||||||||||||||||
! contains the following subroutines
! - ROOT
! - Newaxs
! - Brdevl
! - Neighb
! - FlowNeighb
! - ElemNeighb
! - Brnage
! - Uniera
! - Length
! - Ssgcom1
! - Ssgcom2
! - SortMatrix
! - RotMatrix
! - CondDecision
! - Geocom
! - Prfang
! - Nwcomp
! - Angchg
! - Maslen
! - Mkrecd
! - Grow
! - Adjust
! - Spacng
! - Establ
! - Boxlim
! - BoxlimCylinder
! - Compon
! - Remove
! - Update
! - Secondary_Radial_Growth
! - SplitMacro
! - ShiftRoot
! ==============================================================================
!> main subroutine for root system growth
SUBROUTINE ROOT(t,dt)
  USE typedef
  USE ParamData, ONLY: pi
  USE RootData
  USE GridData, ONLY: nPt,geom,continu,nex,ney,dxgrid,dygrid
  USE tmctrl, ONLY: dtroot,kaxemg
  USE SolData, ONLY: conc,hNew,par,nmat
  USE DoussanMat, ONLY: nplant,Lr,transroot,transtip
  USE StrData, ONLY: simp

  IMPLICIT NONE
  REAL(sp),INTENT(in) ::t,dt
  REAL(dp) ::  cloc,sloc
  REAL(dp) :: newlen,newmas,sumgrw
  REAL(sp) :: dx,dy,dz,space,tempage
  INTEGER(ap) :: corner(8),nrecol
  INTEGER(ap) :: ngrwnw,iseg,igrow,imin,ipl
  !> \param t current simulation time
  !> \param dt current time step
  
  nrecol = nrec
  ngrwnw = ngrow
  sAvg   = 0.0_dp
  cAvg   = 0.0_dp
  sumgrw = 0.0_dp

  ! Secondary radial growth for all established segments
  IF(l_secrad) CALL Secondary_Radial_Growth(dt)

  ! check if it is time to originate new axis(es):
1 IF (kaxemg.LE.naxemg) THEN
     IF (t.GT.tnewax(kaxemg)) THEN
        CALL Newaxs(t,sumgrw,ngrwnw)
        kaxemg = kaxemg+1
        GOTO 1
     ENDIF
  ENDIF

  ! apply potential growth to all growing branches:
  tiploop:  DO igrow = 1,ngrow
     ! develop new branches at all established points that are 'ripe':
     CALL Brdevl(t,sumgrw,igrow,ngrwnw)
    ! find surrounding grid points:
     IF(continu) THEN
        CALL Neighb(xg(igrow)+transtip(igrow,1,1,1)*nex*dxgrid, &
             yg(igrow)+transtip(igrow,2,1,1)*ney*dygrid,zg(igrow),corner,imin)
     ELSE
        CALL Neighb(xg(igrow),yg(igrow),zg(igrow),corner,imin)
     END IF
     ! calculate representative soil strength value felt by root tip:
     CALL StrLoc(corner,sLoc)
     sAvg = sAvg+sLoc
     ! calculate branch age:
     CALL Brnage(t,igrow,tempage) 
     ! calculate new heading angle and tentative segment length:
     CALL Nwcomp(igrow,corner,newlen,dx,dy,dz,tempage)
     ! calculate tentative segment mass:
     CALL Maslen(igrow,sLoc,ordgrw(igrow))
     newmas = newlen*MPL(igrow)
     ! 'grow', tentatively:
     CALL Grow(igrow,dx,dy,dz,newlen,newmas,t)
     ! add up tentative need for dry mass:
     sumgrw = sumgrw+newmas
     ! make a new, tentative record:
  END DO tiploop

  ! if shoot is considered, increase root mass by accumulated increment 'dmroot' or by potential growth:
  IF (lCalloc) THEN
     mroot = mroot+MIN(dmroot,sumgrw)
  ELSE
     mroot = mroot+sumgrw
  ENDIF

  ! average soil strength experienced by growing root system:
  sAvg = sAvg/REAL(ngrwnw)
  ! calculate growth factor from potential growth and available assimilates:
  IF ((lCalloc).AND.(sumgrw.GT.0.0_dp)) THEN
     grwfac = dmroot/sumgrw
  ELSE
     grwfac = 0.0_dp
  ENDIF

  ! reset 'dmroot':
  dmroot = 0.0_dp 
  ALLOCATE(toosml(1:nrec+maxgrw))
  !ALLOCATE(stopgr(1:nbr+maxgrw))
  toosml = .FALSE.
  !stopgr=.FALSE.

  ! calculate total branch lengths
  DO igrow = 1,ngrwnw
     IF (lCalloc) THEN
        ! adjust growth according to available assimilates:
        CALL Adjust(igrow,newlen)
     ELSE
        newlen = seglen(irecsg(igrow))
        IF (newlen.GT.1.E-06_dp) THEN
           toosml(irecsg(igrow)) = .FALSE.
           brlgth(igrow) = brlgth(igrow)+newlen
        ELSE
           toosml(irecsg(igrow))=.TRUE.
        ENDIF
!!$        !DOES NOT WORK PROPERLY AT THE MOMENT!!!
!!$        IF (brlgth(igrow)+newlen.GE.brlmax(ordgrw(igrow))) THEN
!!$           newlen=brlmax(ordgrw(igrow))-brlgth(igrow)
!!$           stopgr(igrow)=.TRUE.
!!$        ELSE
!!$           stopgr(igrow)=.FALSE.
!!$        ENDIF
     ENDIF

     IF ((.NOT.(stopgr(igrow))).AND.(.NOT.(toosml(irecsg(igrow)))).AND.(ordgrw(igrow).LT.norder)) THEN
        ! calculate branch spacing along new segment:
        CALL Spacng(igrow,space)
        ! ...and establish new branching points:
        CALL Establ(igrow,newlen,space,t)
     ENDIF
  END DO 

  ! remove all new segments smaller than 1.E-6:
  CALL Remove(nrecol)
  DO igrow = 1,ngrwnw
     IF(geom.EQ.3 .AND. nmat.GT.1) then
        CALL Boxlim_Cylinder(igrow)
     ELSE
        CALL Boxlim(igrow)
     END IF
  ! remove all new segments smaller than 1.E-6:
!!$     IF (stopgr(igrow)) THEN
!!$        !< branch has reached its maximum length, make final record:
!!$        CALL Mkrecd(igrow,0._dp,0._dp,t+dtRoot)
!!$     ENDIF
  END DO

  ! remove branches that have stopped growth from list of growing branches
  ! and update the number of growing tips:
  CALL Update(ngrwnw)
  ! calculate average solute concentration felt by the root system:
  IF(ltoxi) THEN
     DO iseg = 1,nrec
        CALL Neighb(xs(iseg),ys(iseg),zs(iseg),corner,imin)
        cLoc = (Conc(corner(1))+Conc(corner(2))+Conc(corner(3))+Conc(corner(4))+&
             Conc(corner(5))+Conc(corner(6))+Conc(corner(7))+Conc(corner(8)))/8
        cAvg = cAvg+cLoc
     END DO
     cAvg = cAvg/REAL(nrec)
  ELSE
     cAvg = 0._dp
  END IF
  DEALLOCATE(toosml)
  !DEALLOCATE(stopgr)
  
  ! continuous domain; determine all the transroot for the new tips
  IF(continu .AND. .NOT.lDou) THEN
     DO igrow = 1,ngrow
        CALL tiptrans(xg(igrow),yg(igrow),igrow,1,1)
     END DO
  END IF
  
  RETURN
END SUBROUTINE ROOT
!********************************************************************************
!> originates new axes
SUBROUTINE Newaxs(t,sumgrw,ngrwnw)
  USE typedef
  USE RootData
  USE ParamData
  USE StrData, ONLY : ssmax
  USE tmctrl, ONLY: kaxemg
  USE DoussanMat, ONLY: transroot,transtip,nsub
  USE GridData, ONLY: continu,nex,ney,dxgrid,dygrid,dzgrid
  IMPLICIT NONE

  REAL(sp),INTENT(in) ::t
  REAL(dp), INTENT(inout) :: sumgrw
  REAL(dp) :: newmas,newlen,unimplen,sLoc  !,MPL
  REAL(sp) ::  tloc,alpha,v,tgrow
  REAL(sp) ::  dx,dy,dz,dxini,dyini,dzini,dxstr,dystr,dzstr,drivex,drivey,drivez
  REAL(sp), DIMENSION (3,3) :: condmatrix
  INTEGER(ap), INTENT(inout) :: ngrwnw
  INTEGER(ap) :: imin,i,corner(8)
  !> \param sumgrw total root dry mass
  !> \param ngrwnw number of root tips (new growth)
 
  ! generate a total of 'nnewax' new axes:
  DO i = 1,nnewax(kaxemg)
     ! increment total axes number:
     naxes = naxes+1
     ! increment total branch number and total number of growing branches:
     nbr = nbr+1
     IF (ngrwnw.LT.maxgrw) THEN
        ngrwnw = ngrwnw+1
     ELSE
        WRITE(*,'(''Maximum number of growing tips -- PROGRAM TERMINATED.'')')
        STOP
     ENDIF
     ! assign values for axis number, branching order and branch number:
     iaxis(ngrwnw)  = naxes
     ordgrw(ngrwnw) = 1
     ibrgrw(ngrwnw) = nbr
     ! as of now, the new branch has no estbl. brnch. pts. itself:
     nestbl(ngrwnw)=0
     ! also, length is zero and no pending brnch. pts. exist:
     brlgth(ngrwnw)  = 0.0_dp
     ovrtime(ngrwnw) = -1._sp
     ! no segment behind the tip:
     irecsg(ngrwnw) = 1
     ! emergence position:
     xg(ngrwnw) = xs(1)
     yg(ngrwnw) = ys(1)
     zg(ngrwnw) = zs(1)
     ! find surrounding grid points:
     IF(continu) THEN
        CALL tiptrans(xg(ngrwnw),yg(ngrwnw),ngrwnw,1,1)
        nsub(nrec+ngrwnw,1) = 1
        CALL Neighb(xg(ngrwnw)+transtip(ngrwnw,1,1,1)*nex*dxgrid, &
             yg(ngrwnw)+transtip(ngrwnw,2,1,1)*ney*dygrid,zg(ngrwnw),corner,imin)
     ELSE
        CALL Neighb(xg(ngrwnw),yg(ngrwnw),zg(ngrwnw),corner,imin)
     END IF
     ! choose approach and get soil strength gradient components: 
	 IF (l_conduc) THEN
	 CALL Ssgcom2(corner,condmatrix)
	 ELSE 
	 CALL Ssgcom1(corner,1._dp,strsen(1),dxstr,dystr,dzstr)
	 END IF 
	 
     ! current local temperature:
     if (ltemp) CALL TemLoc(corner,tLoc)


	 
     ! calculate the direction of the first sgement of the new axis:
     CALL RANDOM_NUMBER(rand)
	 alpha = rand*2._sp*pi
	 CALL Initial(1,alpha,inaxs(i),dxini,dyini,dzini)
     ! add up all components:
		dx = dxini
		dy = dyini
		dz = dzini

     ! calculate representative soil strength value felt by root tip:
     CALL StrLoc(corner,sLoc)
     sAvg = sAvg+sLoc
     ! calculate length of the new segment from time left for growth:
     CALL Uniera(0._sp,1,v)
     CALL Length(v,t-tnewax(kaxemg),drivex,drivey,drivez,condmatrix,newlen,unimplen,corner)
	 
     ! check if maximum length has been reached:
     IF (newlen.GE.brlmax(1))  newlen = brlmax(1)
     ! make sure new segment is not too small (prevent divison by zero):
     newlen = MAX(newlen,1.E-4_dp)
     ! calculate tentative segment mass:
     CALL Maslen(ngrwnw,sLoc,ordgrw(ngrwnw))
     newmas = newlen*MPL(ngrwnw)
     ! 'grow', tentatively:
     CALL Grow(ngrwnw,dx,dy,dz,newlen,newmas,tnewax(kaxemg))
     ! add up tentative need for dry mass:
     sumgrw = sumgrw+newmas
  END DO

  RETURN
END SUBROUTINE Newaxs
!****************************************************************************
!> originates new sub-branches from branch 'igrow'
!> branch needs a certain age to branch and there is a certain inter-branch time
SUBROUTINE Brdevl(t,sumgrw,igrow,ngrwnw)
  USE typedef
  USE RootData
  USE ParamData
  USE StrData, ONLY : ssmax
  USE GridData, ONLY: continu,nex,ney,dxgrid,dygrid
  USE DoussanMat, ONLY: transroot,transtip,nsub
  IMPLICIT NONE

  INTEGER(ap) :: irec,iprev,igrow,iest,ngrwnw,kest
  REAL(dp) :: newlen,unimplen,newmas,sumgrw,sLoc   !,MPL
  REAL(sp) :: alpha,beta,gamma,delta,deltat,v
  REAL(sp) :: t,test,tgrow,drivex,drivey,drivez
  REAL(sp) :: dx,dy,dz,t1,t2,x1,x2,y1,y2,z1,z2
  REAL(sp), DIMENSION (3,3) :: condmatrix
  INTEGER(ap) :: corner(8),imin
  !> \param t current time
  !> \param sumgrw total root dry mass
  !> \param igrow current tip or branch number
  !> \param ngrwnw number of root tips (new growth)
  
  kest = 0
  DO iest = 1,nestbl(igrow)
     test = timest(igrow,iest)
     IF ((t-test).GT.dtbrch(ordgrw(igrow))) THEN
        ! we have an established point 'ripe' for branching:
        kest = kest+1
        ! increment total branch number and total number of growing branches:
        nbr = nbr+1
        IF (ngrwnw.LT.maxgrw) THEN
           ngrwnw = ngrwnw+1
        ELSE
           WRITE(*,'(/''Maximum number of growing tips -- PROGRAM TERMINATED.'')')
           STOP
        ENDIF
        ! assign values for axis number, branching order and branch number:
        iaxis(ngrwnw)  = iaxis(igrow)
        ordgrw(ngrwnw) = ordgrw(igrow)+1
        ibrgrw(ngrwnw) = nbr
        ! as of now, the new branch has no estbl. brnch. pts. itself:
        nestbl(ngrwnw) = 0
        ! also, length is zero and no pending brnch. pts. exist:
        brlgth(ngrwnw)  = 0.0_dp
        ovrtime(ngrwnw) = -1._sp
        ! find the segment where the branching occurs:
        irec=irecsg(igrow)
        IF (timorg(irec).LE.test) THEN
           ! have branching in the most recent segment:
           iprev = irec
           x2    = xg(igrow)
           y2    = yg(igrow)
           z2    = zg(igrow)
           t2    = t
        ELSE
           ! need to go back through the segment records:
1          IF(timorg(irecpr(irec)).GT.test) THEN
           irec = irecpr(irec)
           GOTO 1
        ENDIF
        iprev = irecpr(irec)
        x2 = xs(irec)
        y2 = ys(irec)
        z2 = zs(irec)
        t2 = timorg(irec)
     ENDIF
     ! this is also the segment behind the tip:
     irecsg(ngrwnw) = iprev
     x1             = xs(iprev)
     y1             = ys(iprev)
     z1             = zs(iprev)
     t1             = timorg(iprev)
     dx             = x2-x1
     dy             = y2-y1
     dz             = z2-z1
     deltat         = t2-t1


     ! --> small error, but origination point of new branch should be at an existing root node for Doussan!
     ! branch originates at x2,y2,z2
     xg(ngrwnw) = x2
     yg(ngrwnw) = y2
     zg(ngrwnw) = z2

     ! find surrounding grid points:
     IF(continu) THEN
        CALL tiptrans(xg(ngrwnw),yg(ngrwnw),ngrwnw,1,1)
        nsub(nrec+ngrwnw,1) = 1
        CALL Neighb(xg(ngrwnw)+transtip(ngrwnw,1,1,1)*nex*dxgrid, &
             yg(ngrwnw)+transtip(ngrwnw,2,1,1)*ney*dygrid,zg(ngrwnw),corner,imin)
     ELSE
        CALL Neighb(xg(ngrwnw),yg(ngrwnw),zg(ngrwnw),corner,imin)
     END IF

     ! calculate representative soil strength value felt by root tip:
     CALL StrLoc(corner,sLoc)
     sAvg = sAvg+sLoc

     ! calculate initial heading components:
     CALL RANDOM_NUMBER(rand)
	 gamma = rand*2._sp*pi
     delta = brnang(ordgrw(igrow))
     CALL Angchg(dx,dy,dz,alpha,beta,gamma,delta)

     ! calculate length of the new segment from time left for growth:
     CALL Uniera(0.0_sp,ordgrw(ngrwnw),v)
     CALL Length(v,t-dtbrch(ordgrw(igrow))-test,drivex,drivey,drivez,condmatrix,newlen,unimplen,corner)

     ! check if maximum length has been reached:
     IF (newlen.GE.brlmax(ordgrw(ngrwnw))) THEN
        newlen = brlmax(ordgrw(ngrwnw))
     ENDIF

     ! make sure new segment is not too small (prevent divison by zero):
     newlen = MAX(newlen,1.E-04_dp)

     ! calculate tentative segment mass:
     CALL Maslen(ngrwnw,sLoc,ordgrw(ngrwnw))
     newmas = newlen*MPL(ngrwnw)

     ! 'grow', tentatively:
     CALL Grow(ngrwnw,dx,dy,dz,newlen,newmas,test+dtbrch(ordgrw(igrow)))
     ! add up tentative need for dry mass:
     sumgrw = sumgrw+newmas
  ENDIF
END DO
! remove all developed points from the list of established branching points:
nestbl(igrow) = nestbl(igrow)-kest
DO iest = 1,nestbl(igrow)
   timest(igrow,iest) = timest(igrow,iest+kest)
END DO

RETURN
END SUBROUTINE Brdevl
!***************************************************************************
!> finds 8 corner points of cube (double element) that contains root tip
!> or segment originating or ending point 
SUBROUTINE Neighb(x,y,z,corner,imin)
  USE typedef
  USE GridData
  USE DomData
  IMPLICIT NONE

  REAL(sp), INTENT(in) :: x,y,z
  REAL(sp):: delta1
  INTEGER(ap):: ixmin,ixymin
  INTEGER(ap), INTENT(out) :: imin,corner(8)
  !> \param x root x coordinate
  !> \param y root y coordinate
  !> \param z root z coordinate
  !> \param corner soil corner coordinates of voxel surrounding root coords  
  !> \param nElm = total number of elements
  !>  ne* = number of elements (half-cuboids) in '*'-direction
  !> \param nel = nex*ney (number of elements per horizontal element-layer)
  !> \param imin = element-ID# of odd-numbered half of the cuboid containing the root
  !>        tip or segment
  
  

  ! find imin by checking distance between tip coord´s and cuboid center-points
  delta1 = x-(xmin)
  IF (delta1.EQ.nex*dxgrid) THEN
     ixmin = nex
  ELSEIF (delta1.GT.nex*dxgrid) THEN
     print*, delta1, nex*dxgrid, nex, xmin,dxgrid, x
     STOP 'Root node out of the soil domain (Xdirection)'		
  ELSE
     ixmin = 1+FLOOR(delta1/dxgrid)
  END IF

  ! search along y-axis (start at ixmin):
  delta1 = y-(ymin)
  IF (delta1.EQ.ney*dygrid) THEN
     ixymin = ixmin + (FLOOR(delta1/dygrid)-1)*ney
  ELSEIF (delta1.GT.ney*dygrid) THEN
     print*, delta1, ney*dygrid, ney, ymin,dygrid, y 
     STOP 'Root node out of the soil domain (Ydirection)'	
  ELSE
     ixymin = ixmin+FLOOR(delta1/dygrid)*nex
  END IF
  
  ! search along z-axis (start at ixymin):
  delta1 = ABS(z-zmax)
  IF (delta1.EQ.nez*dzgrid) THEN
     imin = ixymin + (FLOOR(delta1/dzgrid)-1)*nel
  ELSEIF (delta1.GT.nez*dzgrid) THEN
     print*, delta1, nez*dzgrid, nez, dzgrid, z
     STOP 'Root node out of the soil domain (Zdirection)'
  ELSE
     imin = ixymin+FLOOR(delta1/dzgrid)*nel
  END IF
IF(imin.gt.nelm) print*, x,y,z

  ! assign cuboid corner nodes:
  corner(1) = elmnod(1,imin)
  corner(2) = elmnod(2,imin)
  corner(3) = elmnod(3,imin)
  corner(4) = elmnod(4,imin)
  corner(5) = elmnod(5,imin)
  corner(6) = elmnod(6,imin)
  corner(7) = elmnod(7,imin)
  corner(8) = elmnod(8,imin)

  RETURN
END SUBROUTINE Neighb
!--------------------------------------------------------------
!>  neighbouring nodes of each node of an element 
SUBROUTINE FlowNeighb(corner,Flow_corner)
  USE typedef
  USE GridData
  IMPLICIT NONE

  INTEGER(ap) :: i,c
  INTEGER(ap), INTENT(in) :: corner(8)
  INTEGER(ap), INTENT(out) :: Flow_corner(1:8,1:6)
  !> \param corner nodes of the centre element
  !> \param Flow_corner neighbouring nodes of the centre element
  
  Flow_corner=0
  DO c=1,8
     Flow_corner(c,1) = corner(c) + 1
     Flow_corner(c,2) = corner(c) - 1
     Flow_corner(c,3) = corner(c) + 1 + ney
     Flow_corner(c,4) = corner(c) - 1 - ney
     Flow_corner(c,5) = corner(c) + (nPt/(nez+1))
     Flow_corner(c,6) = corner(c) - (nPt/(nez+1))
     DO i=1,6
        IF (Flow_corner(c,i) .LT. 0) Flow_corner(c,i) = 0
     ENDDO
  END DO

  RETURN
END SUBROUTINE FlowNeighb
!--------------------------------------------------------------
!>  neighbouring elements of an element iElmx 
SUBROUTINE ElemNeighb(iElm,NeighElm)
  USE typedef
  USE GridData
  IMPLICIT NONE

  INTEGER(ap) :: ix,iy,iz,i
  INTEGER(ap), INTENT(in) :: iElm  
  INTEGER(ap), INTENT(out) :: NeighElm(1:6)
  !> \param iElm centre element
  !> \param NeighElm neighbouring elements (6)
  
  NeighElm = 0
  !> z-direction
  iz = int(ceiling(real(iElm)/real(nex)/real(ney)))
  IF (iz.LE.1) THEN
     NeighElm(5) = 0
  ELSE
     NeighElm(5) = iElm + (nELm/nez)
  END IF
  IF (iz.GE.nez) THEN
     NeighElm(6) = 0
  ELSE
     NeighElm(6) = iElm - (nElm/nez)
  END IF

  !> y-direction
  iy = int(ceiling(real(iElm-(iz-1)*nex*ney) /real(nex)))
  IF (iy.LE.1) THEN
     NeighElm(3) = 0
  ELSE
     NeighElm(3) = iElm - nex
  END IF
  IF (iy.GE.ney) THEN
     NeighElm(4) = 0
  ELSE
     NeighElm(4) = iElm + nex
  END IF

  !> x-direction
  ix = iElm - (iz-1)*nex*ney - (iy-1)*nex
  IF (ix.LE.1) THEN
     NeighElm(1) = 0
  ELSE
     NeighElm(1) = iElm - 1
  END IF
  IF (ix.GE.nex) THEN
     NeighElm(2) = 0
  ELSE
     NeighElm(2) = iElm + 1
  END IF  
  DO i = 1,6
     IF (NeighElm(i).LT.0 .OR. NeighElm(i).GT.nElm) NeighElm(i)=0
  ENDDO

  RETURN
END SUBROUTINE ElemNeighb
!**********************************************************************
!>  defines branch age 
SUBROUTINE Brnage(t,igrow,tempage)
  USE typedef
  USE RootData
  IMPLICIT NONE

  INTEGER*4 igrow,irec
  REAL(sp), INTENT(in) :: t
  REAL(sp), INTENT(out) :: tempage
  !> \param t current time
  !> \param igrow current branch or tip number
  !> \param tempage branch age (temporal parameter name)
  
  irec=irecsg(igrow)
  !> go back through the segment records to first record of branch 'igrow':
1 IF (irecpr(irec).NE.0) THEN
     IF (ordseg(irecpr(irec)).EQ.ordseg(irec)) THEN
        irec = irecpr(irec)
        GOTO 1
     ENDIF
  ENDIF
  tempage = t-timorg(irec)
  IF(tempage.lt.0._sp) THEN
     timorg(irec) = t
     tempage      = 0._sp
  END IF

  RETURN
END SUBROUTINE Brnage
!******************************************************************************
!> unimpeded (= maximum) elongation rate
SUBROUTINE Uniera(tempage,iorder,v)
  USE Typedef
  USE RootData
  IMPLICIT NONE

  REAL(sp), INTENT(in) ::  tempage
  REAL(sp), INTENT(out) :: v
  INTEGER(ap), INTENT(in) :: iorder
  INTEGER(ap):: ivch
  !> \param tempage branch age (temporal parameter name)
  !> \param iorder order of current brach
  !> \param v growth velocity
  
  ! calculate unimpeded elongation rate as a function of order, age:
  IF (tempage.GE.agevch(iorder,nvch(iorder))) THEN
     v = vch(iorder,nvch(iorder))
  ELSE
     ivch = nvch(iorder)
1    ivch = ivch-1
     IF ((tempage.LT.agevch(iorder,ivch)).AND.(ivch.GT.1)) GOTO 1
     v = vch(iorder,ivch)+(tempage-agevch(iorder,ivch))/(agevch(iorder,ivch+1)-&
          agevch(iorder,ivch))*(vch(iorder,ivch+1)-vch(iorder,ivch))
  ENDIF

  RETURN
END SUBROUTINE Uniera
!******************************************************************************
!> check whether the root is in a macropore or not 
SUBROUTINE CheckMacro(corner,lMacro)
  USE typedef
  USE StrData
  IMPLICIT NONE

  LOGICAL, INTENT(out) :: lMacro
  INTEGER(ap), INTENT(in) ::  corner(8)
  INTEGER(sp) ::  i
  !> \param lMacro
  !> \param corner soil nodes surrounding growing tip
 
  ! check
  lMacro=.false.
  DO i=1,8
	IF (s(corner(i)).LE.1.E-1) THEN
	lMacro=.true.
	EXIT
	ENDIF
  END DO

  RETURN
END SUBROUTINE CheckMacro
!******************************************************************************
!> length of new segment 
SUBROUTINE Length(v,tgrow,drivex,drivey,drivez,condmatrix,newlen,unimplen,corner)
  USE typedef
  USE TempData
  USE ConData
  USE StrData
  USE RootData, ONLY: ltemp,tnewax,ltoxi,l_conduc
  IMPLICIT NONE
  LOGICAL :: lMacro
  REAL(sp), INTENT(in) ::  v,tgrow, drivex, drivey, drivez
  REAL(sp), DIMENSION (3,3) , INTENT(in):: condmatrix
  REAL(sp) :: factord,drivexun, driveyun, drivezun
  REAL(sp) :: kmidx, kmidy, kmidz,seff 
  REAL(sp) :: cncimp,temimp,strimp
  REAL(dp), INTENT(out) :: newlen,unimplen
  INTEGER(ap), INTENT(in) ::  corner(8)

  !> \param v growth velocity
  !> \param tgrow length of current growth time step
  !> \param newlen length of new segment
  !> \param corner soil nodes surrounding growing tip
 
  ! calculate length of new segment, taking location into account, approach different for soil strength/conductance approach:
  
   IF (l_conduc) THEN 
	!> calculate a unit vector for the driving force F
	factord = (SQRT(drivex**2+drivey**2+drivez**2))	


    drivexun = drivex/factord
	driveyun = drivey/factord
	drivezun = drivez/factord


    kmidx = (condmatrix(1,1)*drivexun+condmatrix(2,1)*driveyun+condmatrix(3,1)*drivezun)
	kmidy = (condmatrix(1,2)*drivexun+condmatrix(2,2)*driveyun+condmatrix(3,2)*drivezun)
	kmidz = (condmatrix(1,3)*drivexun+condmatrix(2,3)*driveyun+condmatrix(3,3)*drivezun)


	
	seff = 1/(SQRT(kmidx**2+kmidy**2+kmidz**2))	
    strimp = 1-(seff/ssmax)
   ELSE
   
		CALL CheckMacro(corner,lMacro) !> root length shall be max length in macropore
		IF (lMacro) THEN
			strimp = max(imps(corner(1)),imps(corner(2)),imps(corner(3)),imps(corner(4)),&
			imps(corner(5)),imps(corner(6)),imps(corner(7)),imps(corner(8)))
		ELSE  
			strimp = (imps(corner(1))+imps(corner(2))+imps(corner(3))+imps(corner(4))+&
			imps(corner(5))+imps(corner(6))+imps(corner(7))+imps(corner(8)))/8._sp
		ENDIF
   ENDIF
  
  
  IF (ltemp) THEN
     temimp = (impt(corner(1))+impt(corner(2))+impt(corner(3))+impt(corner(4))+&
          impt(corner(5))+impt(corner(6))+impt(corner(7))+impt(corner(8)))/8._sp
  ELSE
     temimp = 1._sp
  ENDIF
  IF (ltoxi) THEN
     cncimp = (impc(corner(1))+impc(corner(2))+impc(corner(3))+impc(corner(4))+&
          impc(corner(5))+impc(corner(6))+impc(corner(7))+impc(corner(8)))/8._sp
  ELSE
     cncimp = 1._sp
  ENDIF
  
  newlen = DBLE((cncimp*strimp*temimp)*v*tgrow)
  unimplen = DBLE(v*tgrow)

  RETURN
END SUBROUTINE Length
!******************************************************************************
!> soil strength gradient components of new heading vector
SUBROUTINE Ssgcom1(corner,stndrd,senstr,dxstr,dystr,dzstr)
  USE typedef
  USE ParamData, ONLY: pi
  USE GridData
  USE StrData
  IMPLICIT NONE

  REAL(sp) :: factor,senstr,smax  
  REAL(dp), INTENT(in) :: stndrd
   REAL(sp), INTENT(out) :: dxstr,dystr,dzstr
  INTEGER(ap), INTENT(in) :: corner(8)

  !> \param corner soil nodes surrounding growing tip
  !> \param stndrd length of previous segment
  !> \param senstr sensitivity of growth direction to soil strength
  !> \param dxstr x component of new segment with respect to the soil strength gradient
  !> \param dystr y component of new segment with respect to the soil strength gradient
  !> \param dzstr z component of new segment with respect to the soil strength gradient

  smax = max(s(corner(1)),s(corner(2)),s(corner(3)),s(corner(4)),&
       s(corner(5)),s(corner(6)),s(corner(7)),s(corner(8)))


  !> calculate, normalize, and weight soil strength gradient components:
  dxstr = (s(corner(1))+s(corner(3))+s(corner(5))+s(corner(7))-&
       s(corner(2))-s(corner(4))-s(corner(6))-s(corner(8)))/(4._sp*smax*dxgrid)
  dystr = (s(corner(1))+s(corner(2))+s(corner(5))+s(corner(6))-&
       s(corner(3))-s(corner(4))-s(corner(7))-s(corner(8)))/(4._sp*smax*dygrid)
  dzstr = -(s(corner(1))+s(corner(2))+s(corner(3))+s(corner(4))-&
       s(corner(5))-s(corner(6))-s(corner(7))-s(corner(8)))/(4._sp*smax*dzgrid)


  IF ((ABS(dxstr).GT.1.E-10_sp).OR.(ABS(dystr).GT.1.E-10_sp).OR.(ABS(dzstr).GT.1.E-10_sp)) THEN
     dxstr = dxstr*senstr*stndrd
     dystr = dystr*senstr*stndrd
     dzstr = dzstr*senstr*stndrd
  ENDIF

  RETURN
END SUBROUTINE Ssgcom1
!******************************************************************************
!> soil strength gradient components of new heading vector (K-like approach)
SUBROUTINE Ssgcom2(corner,condmatrix)
  USE typedef
  USE GridData
  USE StrData
  USE RootData
  USE ParamData, ONLY: pi
  IMPLICIT NONE

  REAL(sp) :: smax,kxx_lin,kyy_lin,kzz_lin
  REAL(sp) :: kxx_skew_x,kyy_skew_x,kzz_skew_x
  REAL(sp) :: kxx_skew_y,kyy_skew_y,kzz_skew_y
  REAL(sp) :: kxx_skew_z,kyy_skew_z,kzz_skew_z  
  REAL(sp) :: ortho1,ortho2,skew_x1,skew_x2,skew_y1,skew_y2,skew_z1,skew_z2
  REAL(sp) :: kxx,kxy,kxz,kyx,kyy,kyz,kzx,kzy,kzz
  REAL(sp) :: kx1,ky1,kz1,kx2,ky2,kz2,phi
  INTEGER(sp) :: i
  INTEGER(ap), INTENT(in) :: corner(8)
  REAL(sp), DIMENSION (3,3) :: rotmatrix_x1,rotmatrix_y1,rotmatrix_z1,rotmatrix_x2,rotmatrix_y2,rotmatrix_z2
  REAL(sp), DIMENSION (3,3) :: condmatrix
  REAL(sp), DIMENSION (3,3) :: skewmatrix_in_x, skewmatrix_in_y, skewmatrix_in_z, linmatrix
  REAL(sp), DIMENSION (3,3) :: skewmatrix_mi
  REAL(sp), DIMENSION (3,3) :: skewmatrix_x, skewmatrix_y, skewmatrix_z,skewmatrix
  REAL(sp), DIMENSION (3,4) :: inmatrix
  !> \param corner soil nodes surrounding growing tip
  !> \param stndrd length of previous segment
  !> \param dxstr x component of new segment with respect to the soil Kstrength (conductivity like approach)
  !> \param dystr y component of new segment with respect to the soil Kstrength (conductivity like approach)
  !> \param dzstr z component of new segment with respect to the soil Kstrength (conductivity like approach)
  
   DO i=1,8
		condu(corner(i)) = 1/s(corner(i))		
	IF (condu(corner(i)).GE.1.E2_sp) THEN
		condu(corner(i)) = condMP 
		s(corner(i)) = 1/condMP
	ENDIF 	
  ENDDO
  
    smax = max(s(corner(1)),s(corner(2)),s(corner(3)),s(corner(4)),&
       s(corner(5)),s(corner(6)),s(corner(7)),s(corner(8)))
	
  
	   
  !> conductance1: axes of anisotropy aligned with the Cartesian Coordinate systm

  !> calculate averaged & normalized soil conductance per grid plane (harmonic mean of the resistance in each corner point)
  
	kx1=(condu(corner(1))+condu(corner(3))+condu(corner(5))+condu(corner(7)))/(4._sp)
	kx2=(condu(corner(2))+condu(corner(4))+condu(corner(6))+condu(corner(8)))/(4._sp)
	ky1=(condu(corner(1))+condu(corner(2))+condu(corner(5))+condu(corner(6)))/(4._sp)
	ky2=(condu(corner(3))+condu(corner(4))+condu(corner(7))+condu(corner(8)))/(4._sp)
	kz1=(condu(corner(1))+condu(corner(2))+condu(corner(3))+condu(corner(4)))/(4._sp)
	kz2=(condu(corner(5))+condu(corner(6))+condu(corner(7))+condu(corner(8)))/(4._sp)
	
	
	kxx_lin=((2._sp)/(1._sp/kx1+1._sp/kx2))
	kyy_lin=((2._sp)/(1._sp/ky1+1._sp/ky2))
	kzz_lin=((2._sp)/(1._sp/kz1+1._sp/kz2))

	
	IF ((ABS(kxx_lin).GT.1.E-10_sp).OR.(ABS(kyy_lin).GT.1.E-10_sp).OR.(ABS(kzz_lin).GT.1.E-10_sp)) THEN	
	kxx_lin= nint((kxx_lin) * 1.E4)/1.E4
	kyy_lin= nint((kyy_lin) * 1.E4)/1.E4
	kzz_lin= nint((kzz_lin) * 1.E4)/1.E4
	ENDIF
	
	
	linmatrix = (reshape((/kxx_lin,0._sp,0._sp,0._sp,kyy_lin,&
						0._sp,0._sp,0._sp,kzz_lin/),shape(linmatrix)))


  !> conductance2: x is aligned with the Cartesian x-axis, y&z in 45-degree angle

  !> calculate averaged & normalized conductance per cube section (harmonic mean of the conductance in each corner point)
 
	kx1=(condu(corner(1))+condu(corner(3))+condu(corner(5))+condu(corner(7)))/(4._sp)
	kx2=(condu(corner(2))+condu(corner(4))+condu(corner(6))+condu(corner(8)))/(4._sp)	
	ky1= (condu(corner(3))+condu(corner(4))+condu(corner(1))/2+condu(corner(2))/2+condu(corner(7))/2+condu(corner(8))/2)/(4._sp)
	ky2= (condu(corner(5))+condu(corner(6))+condu(corner(1))/2+condu(corner(2))/2+condu(corner(7))/2+condu(corner(8))/2)/(4._sp)	
	kz2= (condu(corner(1))+condu(corner(2))+condu(corner(3))/2+condu(corner(4))/2+condu(corner(5))/2+condu(corner(6))/2)/(4._sp)	
	kz2= (condu(corner(7))+condu(corner(8))+condu(corner(3))/2+condu(corner(4))/2+condu(corner(5))/2+condu(corner(6))/2)/(4._sp)
	
	!> calculate the conductivities along the diagonal principal axes
	kxx_skew_x=((2._sp)/(1._sp/kx1+1._sp/kx2))
	kyy_skew_x=((2._sp)/(1._sp/ky1+1._sp/ky2))
	kzz_skew_x=((2._sp)/(1._sp/kz1+1._sp/kz2))
	
	
	  !> normalize the 'soil conductivity' with the total conductance (from the axis aligned conductance tensor) 
	IF ((ABS(kxx_skew_x).GT.1.E-10_sp).OR.(ABS(kyy_skew_x).GT.1.E-10_sp).OR.(ABS(kzz_skew_x).GT.1.E-10_sp)) THEN	
	kxx_skew_x= nint((kxx_skew_x)  * 1.E4)/1.E4
	kyy_skew_x= nint((kyy_skew_x)  * 1.E4)/1.E4
	kzz_skew_x= nint((kzz_skew_x)  * 1.E4)/1.E4
	ENDIF
	
	
	!>create a conductivity matrix 
	
	skewmatrix_in_x = (reshape((/kxx_skew_x,0._sp,0._sp,0._sp,kyy_skew_x,&
						0._sp,0._sp,0._sp,kzz_skew_x/),shape(skewmatrix_in_x)))
						

	CALL RotMatrix(phi,rotmatrix_x1,rotmatrix_y1,rotmatrix_z1,rotmatrix_x2,rotmatrix_y2,rotmatrix_z2)

	 skewmatrix_mi = matmul(rotmatrix_x1,skewmatrix_in_x)
	 skewmatrix_x = matmul(skewmatrix_mi,rotmatrix_x2)

  !> conductance3: y is aligned with the Cartesian y-axis, x&z in 45-degree angle


  !> calculate averaged & normalized conductance per cube section (harmonic mean of the conductance in each corner point)
 	
	kx1= (condu(corner(1))+condu(corner(3))+condu(corner(4))/2+condu(corner(2))/2+condu(corner(5))/2+condu(corner(7))/2)/(4._sp)
	kx2= (condu(corner(6))+condu(corner(8))+condu(corner(4))/2+condu(corner(2))/2+condu(corner(5))/2+condu(corner(7))/2)/(4._sp)	
	ky1=(condu(corner(1))+condu(corner(2))+condu(corner(5))+condu(corner(6)))/(4._sp)
	ky2=(condu(corner(3))+condu(corner(4))+condu(corner(7))+condu(corner(8)))/(4._sp)	
	kz1= (condu(corner(2))+condu(corner(4))+condu(corner(1))/2+condu(corner(3))/2+condu(corner(6))/2+condu(corner(8))/2)/(4._sp)
	kz2= (condu(corner(5))+condu(corner(7))+condu(corner(1))/2+condu(corner(3))/2+condu(corner(6))/2+condu(corner(8))/2)/(4._sp)

	!> calculate the conductivities along the diagonal principal axes
	kxx_skew_y=((2._sp)/(1._sp/kx1+1._sp/kx2))
	kyy_skew_y=((2._sp)/(1._sp/ky1+1._sp/ky2))
	kzz_skew_y=((2._sp)/(1._sp/kz1+1._sp/kz2))
	
	
	  !> normalize the 'soil conductivity' with the total strength (from the orthogonal conductivity tensor) 
	IF ((ABS(kxx_skew_y).GT.1.E-10_sp).OR.(ABS(kyy_skew_y).GT.1.E-10_sp).OR.(ABS(kzz_skew_y).GT.1.E-10_sp)) THEN	
	kxx_skew_y= nint((kxx_skew_y)  * 1.E4)/1.E4
	kyy_skew_y= nint((kyy_skew_y)  * 1.E4)/1.E4
	kzz_skew_y= nint((kzz_skew_y)  * 1.E4)/1.E4
	ENDIF
	
	
	!>create a conductivity matrix 

	skewmatrix_in_y = (reshape((/kxx_skew_y,0._sp,0._sp,0._sp,kyy_skew_y,&
						0._sp,0._sp,0._sp,kzz_skew_y/),shape(skewmatrix_in_y)))		
					

	CALL RotMatrix(phi,rotmatrix_x1,rotmatrix_y1,rotmatrix_z1,rotmatrix_x2,rotmatrix_y2,rotmatrix_z2)
	
	 skewmatrix_mi = matmul(rotmatrix_y1,skewmatrix_in_y)
	 skewmatrix_y = matmul(skewmatrix_mi,rotmatrix_y2)

  !> conductance3: z is aligned with the Cartesian z-axis, x&y in 45-degree angle


  !> calculate averaged & normalized conductance per cube section (harmonic mean of the conductance in each corner point)
 
	kx1= (condu(corner(1))+condu(corner(5))+condu(corner(2))/2+condu(corner(6))/2+condu(corner(3))/2+condu(corner(7))/2)/(4._sp)	
	kx2= (condu(corner(4))+condu(corner(8))+condu(corner(2))/2+condu(corner(6))/2+condu(corner(3))/2+condu(corner(7))/2)/(4._sp)	
	ky1= (condu(corner(2))+condu(corner(6))+condu(corner(1))/2+condu(corner(5))/2+condu(corner(4))/2+condu(corner(8))/2)/(4._sp)		
	ky2= (condu(corner(3))+condu(corner(7))+condu(corner(1))/2+condu(corner(5))/2+condu(corner(4))/2+condu(corner(8))/2)/(4._sp)
	kz1=(condu(corner(1))+condu(corner(2))+condu(corner(3))+condu(corner(4)))/(4._sp)
	kz2=(condu(corner(5))+condu(corner(6))+condu(corner(7))+condu(corner(8)))/(4._sp)
	
	!> calculate the conductivities along the diagonal principal axes
	kxx_skew_z=((2._sp)/(1._sp/kx1+1._sp/kx2))
	kyy_skew_z=((2._sp)/(1._sp/ky1+1._sp/ky2))
	kzz_skew_z=((2._sp)/(1._sp/kz1+1._sp/kz2))
	

	  !> normalize the 'soil conductivity' with the total strength (from the orthogonal conductivity tensor) 
	IF ((ABS(kxx_skew_z).GT.1.E-10_sp).OR.(ABS(kyy_skew_z).GT.1.E-10_sp).OR.(ABS(kzz_skew_z).GT.1.E-10_sp)) THEN	
	kxx_skew_z= nint((kxx_skew_z)  * 1.E4)/1.E4
	kyy_skew_z= nint((kyy_skew_z)  * 1.E4)/1.E4
	kzz_skew_z= nint((kzz_skew_z)  * 1.E4)/1.E4
	ENDIF
	

	skewmatrix_in_z = (reshape((/kxx_skew_z,0._sp,0._sp,0._sp,kyy_skew_z,&
						0._sp,0._sp,0._sp,kzz_skew_z/),shape(skewmatrix_in_z)))	
		
		
	CALL RotMatrix(phi,rotmatrix_x1,rotmatrix_y1,rotmatrix_z1,rotmatrix_x2,rotmatrix_y2,rotmatrix_z2)
	
	 skewmatrix_mi = matmul(rotmatrix_z1,skewmatrix_in_z)
	 skewmatrix_z = matmul(skewmatrix_mi,rotmatrix_z2)

	
	!> Find the maximum values of all conductivity tensors
	inmatrix = (reshape((/kxx_lin,kyy_lin,kzz_lin,kxx_skew_x,kyy_skew_x,kzz_skew_x,&
						kxx_skew_y,kyy_skew_y,kzz_skew_y, kxx_skew_z,&
						kyy_skew_z, kzz_skew_z/),shape(inmatrix)))	

	CALL Sortmatrix(inmatrix,ortho1,skew_x1,skew_y1,skew_z1,ortho2,skew_x2,skew_y2,skew_z2)

	!>Decision Rule
	!> 1) take the tensor with the maximum cond in one of the three directions
	IF (skew_x1.GT.MAX(ortho1,skew_y1,skew_z1))  	THEN	
		condmatrix = skewmatrix_x	
	ELSEIF (skew_y1.GT.MAX(ortho1,skew_x1,skew_z1)) THEN
		condmatrix = skewmatrix_y
	ELSEIF (skew_z1.GT.MAX(ortho1,skew_x1,skew_y1))	THEN
		condmatrix = skewmatrix_z
	!> 2) if two or more of these values are the same, take the tensor with the highest second largest value	
	ELSEIF (skew_x2.GT.MAX(ortho2,skew_y2,skew_z2)) THEN	
		condmatrix = skewmatrix_x	
	ELSEIF (skew_y2.GT.MAX(ortho2,skew_x2,skew_z2)) THEN
		condmatrix = skewmatrix_y
	ELSEIF (skew_z2.GT.MAX(ortho2,skew_x2,skew_y2))	THEN
		condmatrix = skewmatrix_z	
	!> 3) if two of the skew conductivity tensors are the same AND if their maximum values are higher then the maximum value of the orthogonal tensor		
	ELSEIF ((ortho1.LT.MAX(skew_x1,skew_y1,skew_z1)).AND.&
               (MAX(skew_x1,skew_y1,skew_z1).NE.MIN(skew_x1,skew_y1,skew_z1))) THEN		
		CALL CondDecision1(skew_x1,skew_y1,skew_z1,skewmatrix_x,skewmatrix_y,skewmatrix_z,skewmatrix)
		condmatrix = skewmatrix
	!> 4) if the three skew conductivity tensors are the same and if their maximum value is hogher than that of the orthogonal tensor
	ELSEIF (ortho1.LT.MAX(skew_x1,skew_y1,skew_z1)) THEN
		CALL CondDecision2(skewmatrix_x,skewmatrix_y,skewmatrix_z,skewmatrix)
		condmatrix = skewmatrix
	ELSE
	!> in all other cases: take the orthogonal conductivity tensor! 
		condmatrix = linmatrix
	ENDIF

 !> take the k-values out of the conductivity matrix 
 
	kxx = condmatrix(1,1)
	kxy = condmatrix(2,1)
	kxz = condmatrix(3,1)
	kyx = condmatrix(1,2)
	kyy = condmatrix(2,2)
	kyz = condmatrix(3,2)
	kzx = condmatrix(1,3)
	kzy = condmatrix(2,3)
	kzz = condmatrix(3,3)
	
   
  RETURN
END SUBROUTINE Ssgcom2

!******************************************************************************
SUBROUTINE Sortmatrix(inmatrix,ortho1,skew_x1,skew_y1,skew_z1,ortho2,skew_x2,skew_y2,skew_z2)
  USE typedef
  IMPLICIT NONE

  REAL(sp), DIMENSION (3,4), INTENT(in) :: inmatrix
  REAL(sp), DIMENSION (3,4) :: smatrix
  REAL(sp), INTENT(out) :: ortho1, skew_x1, skew_y1, skew_z1
  REAL(sp), INTENT(out) :: ortho2, skew_x2, skew_y2, skew_z2
  REAL(sp) :: dummy
  INTEGER(sp) :: hh,ii,jj
  
  
  smatrix=inmatrix
  
  DO hh=1,4
	DO ii=1,3
		DO jj=ii,3
			IF (smatrix(ii,hh).LE.smatrix(jj,hh)) THEN 
				dummy = smatrix(jj,hh)
				smatrix(jj,hh) = smatrix(ii,hh)
				smatrix(ii,hh) = dummy
			ENDIF
		ENDDO
	ENDDO
ENDDO

	
	ortho1 =  smatrix(1,1)
	skew_x1 = smatrix(1,2)
	skew_y1 = smatrix(1,3)
	skew_z1 = smatrix(1,4)
	ortho2 =  smatrix(2,1)
	skew_x2 = smatrix(2,2)
	skew_y2 = smatrix(2,3)
	skew_z2 = smatrix(2,4)
	
END SUBROUTINE Sortmatrix
!******************************************************************************
!> Random decision which of two equivalent conductivities should be chosen
SUBROUTINE CondDecision1(skew_x1,skew_y1,skew_z1,skewmatrix_x,skewmatrix_y,skewmatrix_z,skewmatrix)
  USE typedef
  USE RootData, ONLY: rand
  IMPLICIT NONE

  REAL(sp) :: u,j
  REAL(sp), INTENT(in):: skew_x1,skew_y1,skew_z1
  REAL(sp), DIMENSION (3,3), INTENT(in) :: skewmatrix_x, skewmatrix_y, skewmatrix_z
  REAL(sp), DIMENSION (3,3), INTENT(out) :: skewmatrix
	CALL RANDOM_NUMBER(rand)
	u = rand
	j = FLOOR(2*u) 
	
	IF (skew_x1.EQ.skew_y1) THEN
			IF (j.EQ.0) THEN
			skewmatrix = skewmatrix_x
			ELSE
			skewmatrix = skewmatrix_y
			ENDIF
		ELSEIF 	(skew_x1.EQ.skew_z1) THEN
			IF (j.EQ.0) THEN
			skewmatrix = skewmatrix_x
			ELSE
			skewmatrix = skewmatrix_z
			ENDIF	
		ELSE		
			IF (j.EQ.0) THEN
			skewmatrix = skewmatrix_y
			ELSE
			skewmatrix = skewmatrix_z
			ENDIF	
	ENDIF 		


END SUBROUTINE CondDecision1
!******************************************************************************
!> Random decision which of three equivalent conductivity tensors shall be chosen
SUBROUTINE CondDecision2(skewmatrix_x,skewmatrix_y,skewmatrix_z,skewmatrix)
  USE typedef
  USE RootData, ONLY: rand
  IMPLICIT NONE

  REAL(sp) :: u,j
  REAL(sp), DIMENSION (3,3), INTENT(in) :: skewmatrix_x, skewmatrix_y, skewmatrix_z
  REAL(sp), DIMENSION (3,3), INTENT(out) :: skewmatrix

	CALL RANDOM_NUMBER(rand)
	u = rand
	j = FLOOR(3*u) 
	
		IF (j.EQ.0) THEN
			skewmatrix = skewmatrix_x
		ELSEIF (j.EQ.1) THEN
			skewmatrix = skewmatrix_y
		ELSE 
			skewmatrix = skewmatrix_z
		ENDIF


END SUBROUTINE CondDecision2
!******************************************************************************
!> Rotation Matrices for 45° rotation around the main axes
SUBROUTINE RotMatrix(phi,rotmatrix_x1,rotmatrix_y1,rotmatrix_z1,rotmatrix_x2,rotmatrix_y2,rotmatrix_z2)
  USE typedef
  USE ParamData, ONLY: pi
  IMPLICIT NONE

  REAL(sp) :: phi,theta,psi
  REAL(sp), DIMENSION (3,3), INTENT(out) :: rotmatrix_x1,rotmatrix_y1,rotmatrix_z1
  REAL(sp), DIMENSION (3,3), INTENT(out) :: rotmatrix_x2,rotmatrix_y2,rotmatrix_z2
  
  
  		phi = 45._sp/180._sp*pi
 
		rotmatrix_x1 = (reshape((/1._sp,0._sp,0._sp,0._sp,cos(phi),sin(phi),0._sp,(-sin(phi)),cos(phi)/),shape(rotmatrix_x1)))
		rotmatrix_x2 = (reshape((/1._sp,0._sp,0._sp,0._sp,cos(phi),-sin(phi),0._sp,(sin(phi)),cos(phi)/),shape(rotmatrix_x2)))
		rotmatrix_y1 = (reshape((/cos(phi),0._sp,(-sin(phi)),0._sp,1._sp,0._sp,sin(phi),0._sp,cos(phi)/),shape(rotmatrix_y1)))
		rotmatrix_y2 = (reshape((/cos(phi),0._sp,sin(phi),0._sp,1._sp,0._sp,(-sin(phi)),0._sp,cos(phi)/),shape(rotmatrix_y2)))
		rotmatrix_z1 = (reshape((/cos(phi),sin(phi),0._sp,(-sin(phi)),cos(phi),0._sp,0._sp,0._sp,1._sp/),shape(rotmatrix_z1)))	
		rotmatrix_z2 = (reshape((/cos(phi),(-sin(phi)),0._sp,sin(phi),cos(phi),0._sp,0._sp,0._sp,1._sp/),shape(rotmatrix_z2)))
		
		
END SUBROUTINE RotMatrix
!******************************************************************************
!> geotropism components of new heading vector
SUBROUTINE Geocom(tLoc,iord,iax,stndrd,alpha,dxgeo,dygeo,dzgeo)
  USE typedef
  IMPLICIT NONE

  REAL(sp), INTENT(in)::  alpha
  REAL(sp) :: tLoc,betapr,geotrp
  REAL(sp), INTENT(out) :: dxgeo,dygeo,dzgeo
  REAL(dp), INTENT(in) :: stndrd
  INTEGER (ap), INTENT(in) :: iord,iax
  !> \param tLoc average value of the current soil temperature values of the 8 corner nodes
  !> \param iord current branch order
  !> \param iax current axis number
  !> \param stndrd length of previous segment
  !> \param alpha azimuth angle
  !> \param dxstr x component of new segment with respect to the soil strength gradient
  !> \param dystr y component of new segment with respect to the soil strength gradient
  !> \param dzstr z component of new segment with respect to the soil strength gradient

  ! calculate geotropism components
  ! (betapr = preferrential heading angle with xy-plane):
  IF (iord.LE.4) THEN
     CALL Prfang(tLoc,iord,iax,geotrp,betapr)
     dxgeo = stndrd*geotrp*COS(betapr)*COS(alpha)
     dygeo = stndrd*geotrp*COS(betapr)*SIN(alpha)
     dzgeo = stndrd*geotrp*SIN(betapr)
  ELSE
     dxgeo = 0.0_sp
     dygeo = 0.0_sp
     dzgeo = 0.0_sp
  ENDIF

  RETURN
END SUBROUTINE Geocom
!************************************************************!*

!> initial direction of the new axis
SUBROUTINE Initial(iord,alpha,inaxs,dxini,dyini,dzini)
  USE typedef
  USE typedef
  USE GeoData
  USE RootData, ONLY: rand
  IMPLICIT NONE
  
  REAL(sp), INTENT(in)::  alpha,inaxs
  REAL(sp), INTENT(out) :: dxini,dyini,dzini
  INTEGER (ap), INTENT(in) :: iord
  !> \param tLoc average value of the current soil temperature values of the 8 corner nodes
  !> \param iord current branch order
  !> \param iax current axis number
  !> \param stndrd length of previous segment
  !> \param alpha azimuth angle
  !> \param dxini x component of the first segment of the new axis
  !> \param dyini y component of the first segment of the new axis
  !> \param dzini z component of the first segment of the new axis
  !> \param inaxs angle of the first segment of the new axis (from the horizontal plane)
	
	! calculate components of the initial angle of the new axis:  
  IF (iord.LE.4) THEN	 
     dxini = COS(inaxs)*SIN(inaxs)
     dyini = COS(inaxs)*SIN(alpha)
     dzini = SIN(inaxs)
  ELSE
     dxini = 0.0_sp
     dyini = 0.0_sp
     dzini = 0.0_sp
  ENDIF


  RETURN
END SUBROUTINE Initial
!************************************************************!*
!> preferential growth angle with horizontal plane
SUBROUTINE Prfang(tLoc,iord,iax,geotrp,betapr)
  USE typedef
  USE GeoData
  USE RootData, ONLY: rand
  IMPLICIT NONE

  REAL(sp), INTENT(in)::  tLoc
  REAL(sp), INTENT(out) :: geotrp,betapr
  INTEGER(ap), INTENT(in)::  iord,iax
  INTEGER(ap)::  igch
  !> \param tLoc average value of the current soil temperature values of the 8 corner nodes
  !> \param iord current branch order
  !> \param iax current axis number
  !> \param geotrp weighting factor for geotropism angle of growth direction vector
  !> \param betapr preferential vertical growth angle

  ! interpolate along piecewise linear function to get
  ! preferred growth angle as a function of axis#, temperature:
  IF (iord.EQ.1) THEN
     ! use axis data:
     IF (tLoc.GE.tempax(iax,nangax(iax))) THEN
        ! temperature greater than greatest T for which betapr is specified --
        ! use values for greatest specified T:
        betapr = angaxs(iax,nangax(iax))
     ELSE
        iGch = nangax(iax)
1       iGch = iGch-1
        IF (iGch.EQ.0) THEN
           ! temperature smaller than smallest T for which betapr is specified --
           ! use values for smallest specified T:
           betapr = angaxs(iax,1)
        ELSE
           IF (tLoc.LT.tempax(iax,iGch)) GOTO 1
           betapr = angaxs(iax,iGch)+(tLoc-tempax(iax,iGch))/&
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
        ! addition of random component - for iord=1 already done in Input.f90 
        betapr = anglat(nanglt)
     ELSE
        iGch = nanglt
2       iGch = iGch-1
        IF (iGch.EQ.0) THEN
           ! temperature smaller than smallest T for which betapr is specified --
           ! use values for smallest specified T:
           betapr = anglat(1)
        ELSE
           IF (tLoc.LT.templt(iGch)) GOTO 2
           betapr = anglat(iGch)+(tLoc-templt(iGch))/&
                (templt(iGch+1)-templt(iGch))*(anglat(iGch+1)-anglat(iGch))
           betapr = anglat(nanglt)
        ENDIF
     ENDIF
     geotrp = geolat
  ENDIF

  RETURN
END SUBROUTINE Prfang
!******************************************************************************
!> length and heading vector components of new segment 
SUBROUTINE Nwcomp(igrow,corner,newlen,dx,dy,dz,tempage)
  USE typedef
  USE RootData
  USE ParamData
  USE StrData, ONLY : ssmax
  USE tmctrl, ONLY: dtroot
  IMPLICIT NONE

  REAL(dp):: newlen,unimplen,stndrd
  REAL(sp), INTENT(out) :: dx,dy,dz
  REAL(sp)::  dxstr,dystr,dzstr,dxgeo,dygeo,dzgeo,drivex,drivey,drivez
  REAL(sp)::  alpha,beta,tLoc,v,tgrow
  REAL(sp), DIMENSION (3,3) :: condmatrix
  REAL(sp)::  tempage,gamma,stdev,u1,u2,r,theta,delta
  INTEGER(ap), INTENT(in)::  corner(8),igrow
  INTEGER(ap) ::  itip,iprev
  
  !> \param igrow current tip or branch number
  !> \param corner 8 surrounding soil nodes of growing tip
  !> \param newlen length of new segment 
  !> \param dx x component of new segment 
  !> \param dy y component of new segment 
  !> \param dz z component of new segment 
  !> \param tempage current branch age

 
  !> use old segment length as normalizing standard for components:
  stndrd = seglen(irecsg(igrow))

  !> calculate the old heading components dx,dy,dz:
  itip = irecsg(igrow)
  iprev = irecpr(itip)
  dx = xs(itip)-xs(iprev)
  dy = ys(itip)-ys(iprev)
  dz = zs(itip)-zs(iprev)
  ! dx = xg(igrow)-xs(irecsg(igrow))
  ! dy = yg(igrow)-ys(irecsg(igrow))
  ! dz = zg(igrow)-zs(irecsg(igrow))

  !> change old heading angle by some random amount without changing length:
  CALL RANDOM_NUMBER(rand)
  gamma = (2*rand-1)*2._sp*pi
  
 !> add a random angle with the expected value 0 and a calc stdev to the old heading angle
 !> calc random angle from a random sample from a normal(Gaussian) distribution scaled with the segment length and the range of rdmang
 CALL Uniera(tempage,ordgrw(igrow),v)
  CALL Length(v,dtroot,drivex,drivey,drivez,condmatrix,newlen,unimplen,corner)
	stdev = (unimplen)**0.5*(rdmang(ordgrw(igrow)))
    if (stdev.LT.0.0) THEN
       Print*, 'standard deviation must be positive'
    end IF
	
	CALL RANDOM_NUMBER(rand)
    u1 = rand
	CALL RANDOM_NUMBER(rand)
    u2 = rand
    r = sqrt( -2.0*log(u1) )
    theta = 2.0*pi*u2	
    delta = stdev*r*sin(theta)
  CALL Angchg(dx,dy,dz,alpha,beta,gamma,delta)

  !> get soil strength gradient components:
	 IF (l_conduc) THEN
	 CALL Ssgcom2(corner,condmatrix)
	 ELSE 
	 CALL Ssgcom1(corner,stndrd,strsen(ordgrw(igrow)),dxstr,dystr,dzstr)
	 END IF 

  !> current local temperature:
  if(ltemp) CALL TemLoc(corner,tLoc)

  !< get geotropism components:
  CALL Geocom(tLoc,ordgrw(igrow),iaxis(igrow),stndrd,alpha,dxgeo,dygeo,dzgeo)
  

  !> add up the different direction components according to the used approach (ssgcom1 or ssgcom2)
 IF (l_conduc) THEN 
 

	drivex = (dx+dxgeo)
	drivey = (dy+dygeo)
	drivez = (dz+dzgeo) 


 
  dx = (condmatrix(1,1)*drivex+condmatrix(2,1)*drivey+condmatrix(3,1)*drivez)
  dy = (condmatrix(1,2)*drivex+condmatrix(2,2)*drivey+condmatrix(3,2)*drivez)
  dz = (condmatrix(1,3)*drivex+condmatrix(2,3)*drivey+condmatrix(3,3)*drivez)
  

 ELSE 
  dx = (dxstr+dx+dxgeo)
  dy = (dystr+dy+dygeo)
  dz = (dzstr+dz+dzgeo) 
  
   
 END IF

  !> calculate length of new segment, taking location into account:
  CALL Uniera(tempage,ordgrw(igrow),v) 
  CALL Length(v,dtroot,drivex,drivey,drivez,condmatrix,newlen,unimplen,corner)

  !> check if maximum length has been reached:
  IF (brlgth(igrow)+newlen.GE.brlmax(ordgrw(igrow))) THEN
     newlen = brlmax(ordgrw(igrow))-brlgth(igrow)
     stopgr(igrow) = .TRUE.
  ELSE
     stopgr(igrow) = .FALSE.
  ENDIF

  RETURN
END SUBROUTINE Nwcomp
!******************************************************************************
!> changes the previous azimuth and polar angle by a random component each
SUBROUTINE Angchg(dx,dy,dz,alpha,beta,gamma,delta)
  USE typedef
  USE paramData
  IMPLICIT NONE

  REAL(sp), INTENT(inout) :: dx,dy,dz
  REAL(sp), INTENT(in) :: gamma,delta
  REAL(sp), INTENT(out) :: alpha,beta
  REAL(sp) :: totlen,horlen,a,r,c,d
  !> \param dx x component of new segment 
  !> \param dy y component of new segment 
  !> \param dz z component of new segment 
  !> \param alpha azimuth of the previous segment
  !> \param beta polar angle of the previous segment
  !> \param gamma perturbation of the azimuth
  !> \param delta perturbation of the polar angle

  totlen = SQRT(dx*dx+dy*dy+dz*dz)
  horlen = SQRT(dx*dx+dy*dy)
  r      = SIN(delta)*totlen
  a      = COS(gamma)*r
  beta   = ASIN(SIGN(MIN(1._sp,ABS(dz/totlen)),dz))
  IF (horlen.GT.1.E-20_sp) THEN
     alpha = ACOS(SIGN(MIN(1._sp,ABS(dx/horlen)),dx))
     IF (dy.LT.0.0_sp) alpha = 2._sp*pi-alpha
     c  = SIN(beta)*a
     d  = SIN(gamma)*r
     dx = COS(delta)*dx-c/horlen*dx-SIN(alpha)*d
     dy = COS(delta)*dy-c/horlen*dy+COS(alpha)*d
  ELSE
     alpha = gamma
     dx    = COS(alpha)*r
     dy    = SIN(alpha)*r
  ENDIF
  dz     = COS(delta)*dz+COS(beta)*a
  horlen = SQRT(dx*dx+dy*dy)
  beta   = ASIN(SIGN(MIN(1._sp,ABS(dz/totlen)),dz))
  IF (horlen.GT.1.E-20_sp) THEN
     alpha = ACOS(SIGN(MIN(1._sp,ABS(dx/horlen)),dx))
     IF (dy.LT.0.0_sp) alpha = 2._sp*pi-alpha
  ELSE
     alpha = gamma
  ENDIF

  RETURN
END SUBROUTINE Angchg
!*****************************************************************
!> mass per length
SUBROUTINE Maslen(igrow,sLoc,iorder)
  USE Typedef
  USE RootData, ONLY: MPL,MPLch,nMPLch,sMPLch
  IMPLICIT NONE

  REAL(dp), INTENT(in) :: sLoc
  INTEGER(ap), INTENT(in) :: igrow,iorder
  INTEGER(ap) :: iMPL
  !> \param igrow current growing number of tip or branch
  !> \param sLoc average soil strength of the 8 corner nodes around the growing tip
  !> \param iorder current branch order

  ! calculate mass per length as a function of order, soil strength:
  IF (sLoc.GE.sMPLch(iorder,nMPLch(iorder))) THEN
     MPL(igrow) = MPLch(iorder,nMPLch(iorder))
  ELSE
     iMPL = nMPLch(iorder)
1    iMPL = iMPL-1
     IF ((sLoc.LT.sMPLch(iorder,iMPL)).AND.(iMPL.GT.1)) GOTO 1
     MPL(igrow) = MPLch(iorder,iMPL)+(sLoc-sMPLch(iorder,iMPL))/&
          (sMPLch(iorder,iMPL+1)-sMPLch(iorder,iMPL))*&
          (MPLch(iorder,iMPL+1)-MPLch(iorder,iMPL))
  ENDIF

  RETURN
END SUBROUTINE Maslen
!*****************************************************************
!> grow last segment record; update tip and segment information 
!> merged old subroutines GROW and MKRECD
SUBROUTINE Grow(igrow,dx,dy,dz,length,mass,rectime)
  USE typedef
  USE RootData
  USE PlntData, ONLY: SpWgt
  USE ParamData, ONLY: pi
  USE GridData, ONLY: continu
  USE DoussanMat,ONLY: transroot,nsub
  IMPLICIT NONE

  REAL(sp) :: dx,dy,dz,factor,rectime
  REAL(dp), INTENT(in) :: length,mass
  INTEGER(ap) :: igrow
  !> \param igrow current growing number of tip or branch
  !> \param dx x component of new segment 
  !> \param dy y component of new segment 
  !> \param dz z component of new segment 
  !> \param length length of new segment
  !> \param mass (dry) mass of new segment 
  !> \param rectime origination time of new segment 

  ! increase number of node records 
  IF (nrec.LT.maxrec) THEN
     nrec = nrec+1
  ELSE
     WRITE (*,'(/''Maximum number of root segments -- PROGRAM TERMINATED.'')')
     STOP
  ENDIF 
  ! calculate the actual length of components dx,dy,dz:
  factor = REAL(length)/SQRT(dx*dx+dy*dy+dz*dz)
  dx     = factor*dx
  dy     = factor*dy
  dz     = factor*dz

  ! grow at last root node
  xs(nrec) = xg(igrow)+dx
  ys(nrec) = yg(igrow)+dy
  zs(nrec) = zg(igrow)+dz

  ! make record for new root node
  irecpr(nrec) = irecsg(igrow)
  ordseg(nrec) = ordgrw(igrow)
  ibrseg(nrec) = ibrgrw(igrow)
  seglen(nrec) = length
  segmas(nrec) = mass
  timorg(nrec) = rectime
  IF (length.EQ.0._dp) seglen(nrec) = 1.E-7_dp
  segrad(nrec) = SQRT(MPL(igrow)/SpWgt/pi)
  !print*,'nrec',nrec,'MPL',MPL(igrow),'rad',segrad(nrec)
  IF(segrad(nrec).LT.1.E-7_dp) segrad(nrec) = 1.E-7_dp
  segsur(nrec)          = 2._dp*pi*segrad(nrec)*seglen(nrec)
  crossSectionSeg(nrec) = segrad(nrec)*segrad(nrec)*pi

  IF(continu) THEN
     CALL roottrans(xs(nrec),ys(nrec),nrec,1,1)
     nsub(nrec,1) = 1
  END IF

  ! update tip information
  xg(igrow)     = xs(nrec)
  yg(igrow)     = ys(nrec)
  zg(igrow)     = zs(nrec)
  irecsg(igrow) = nrec
  RETURN
END SUBROUTINE Grow
!********************************************************************************
!> in case shoot growth is considered, adjust root growth to not exceed available assimilate
SUBROUTINE Adjust(igrow,newlen)
  USE typedef
  USE RootData, ONLY: stopgr,toosml,segmas,ordgrw,irecsg,seglen,xs,ys,zs,xg,yg,zg,brlgth,cAvg,brlmax,irecpr,ordseg,grwfac
  USE tmctrl, ONLY: dtroot
  IMPLICIT NONE

  INTEGER(ap) :: irec,igrow
  REAL(dp), INTENT(out) :: newlen
  REAL(sp) :: factor
  !> \param igrow number of growing tip or branch
  !> \param newlen length of new segment

  ! grwfac > 1.  indicates more assimilate being sent to root than can be used
  ! for growth under current conditions --
  ! assume that extra assimilate is exudated.
  ! find the segment behind the tip 'igrow':
  irec = irecsg(igrow)

  ! apply growth factor to tentative segment length:
  newlen = seglen(irec)*MIN(grwfac,1._dp)

  ! check if maximum length has been reached:
  IF (brlgth(igrow)+newlen.GE.brlmax(ordgrw(igrow))) THEN
     newlen = brlmax(ordgrw(igrow))-brlgth(igrow)
     stopgr(igrow) = .TRUE.
  ELSE
     stopgr(igrow) = .FALSE.
  ENDIF

  ! make sure first branch segments are not too small
  ! after being adjusted (prevent div. by zero):
  IF (irecpr(irec).EQ.0) THEN
     newlen = MAX(newlen,1.E-07_dp)
  ELSE
     IF (ordseg(irecpr(irec)).NE.ordseg(irec)) newlen = MAX(newlen,1.E-07_dp)
  ENDIF
  IF (newlen.GT.1.E-06_dp) THEN
     toosml(irec) = .FALSE.
     ! adjust mass of segment:
     segmas(irec) = segmas(irec)*MIN(grwfac,1._dp)
     ! calculate length correction factor:
     factor = REAL(newlen/seglen(irec))
     ! calculate exact position of tip:
     xg(igrow)     = xs(irec)+factor*(xg(igrow)-xs(irec))
     yg(igrow)     = ys(irec)+factor*(yg(igrow)-ys(irec))
     zg(igrow)     = zs(irec)+factor*(zg(igrow)-zs(irec))
     brlgth(igrow) = brlgth(igrow)+newlen
     ! adjust length of segment:
     seglen(irec) = newlen
  ELSE
     ! we have a candidate for removal at the end of time step
     toosml(irec) = .TRUE.
  ENDIF

  RETURN
END SUBROUTINE Adjust
!******************************************************************************
!> time difference between successive sub-branch origination points 
SUBROUTINE Spacng(igrow,space)
  USE typedef
  USE RootData
  IMPLICIT NONE

  INTEGER(ap), INTENT(in) :: igrow
  REAL(sp), INTENT(out) :: space
  !> \param igrow number of growing tip
  !> \param space age distance between two originating branches

  space = brspac(ordgrw(igrow))

  RETURN
END SUBROUTINE Spacng
!******************************************************************************
!> establish new sub-branch origination points along branch 'igrow' 
SUBROUTINE Establ(igrow,newlen,space,t)
  USE typedef
  USE RootData
  USE tmctrl, ONLY: dtroot
  IMPLICIT NONE

  REAL(sp), INTENT(in):: t,space
  REAL(sp) :: first,deltat,over,smin,v_old,v_new
  REAL(dp), INTENT(in):: newlen
  INTEGER(ap), INTENT(in):: igrow
  !> \param igrow number of growing branch or tip
  !> \param newlen length of new segment 
  !> \param space age distance between two originating branches
  !> \param t current simulation time

  over = ovrtime(igrow)  
  IF (over.GT.0._sp) THEN
     smin = MIN(over,space)
     IF (smin.GT.newlen) THEN
        first = -1._sp
        over = smin-dtroot  
     ELSE
        first = smin
     ENDIF
  ELSE
     IF (space.GT.dtroot) THEN 
        first = -1._sp
        over = space-dtroot 
     ELSE
        first = space
     ENDIF
  ENDIF
  IF ((first.GT.0._sp).AND.(nestbl(igrow).LT.maxest)) THEN
     deltat                      = first 
1    nestbl(igrow)               = nestbl(igrow)+1
     timest(igrow,nestbl(igrow)) = t+deltat 
     deltat                      = deltat+space
     IF ((deltat.LE.dtroot).AND.(nestbl(igrow).LT.maxest)) GOTO 1
     over = deltat-dtroot
  ENDIF
  ovrtime(igrow) = over

  RETURN
END SUBROUTINE Establ
!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!**!!!!*
!> root has to stay within domain limits
SUBROUTINE Boxlim(igrow)
  USE typedef
  USE RootData
  USE DomData
  USE GridData, ONLY: continu
  USE tmctrl, ONLY: dtroot
  IMPLICIT NONE

  INTEGER(ap), INTENT(in) :: igrow
  INTEGER(ap) :: irec,iret,imin
  REAL(sp):: dx,dy,dz,newtime
  REAL(dp):: factor,fractn,delmin,del,length,newlen,mass,newmas
  CHARACTER wall*4
  !> \param igrow number of growing tip

  ! find the segment behind the tip 'igrow':
  ! then check each of the boudaries (x,y,z)
1 iret = irecsg(igrow) ! tip node
  irec = irecpr(iret)  ! prev. node
  IF (toosml(iret)) RETURN

  ! check if tip is outside the domain
  ! and if so, which side wall is intersected first:
  delmin = 1.E+37_dp
  wall   = 'none'

  ! domain walls are the growth boundarys
  IF(.NOT.continu) THEN
     IF (xs(iret).LT.xmin) THEN
        del = (xmin-xs(irec))/(xs(iret)-xs(irec))*seglen(iret)
        IF (del.LT.delmin) THEN
           delmin = del
           wall   = 'xmin'
        ENDIF
     ENDIF
     IF (xs(iret).GT.xmax) THEN
        del = (xmax-xs(irec))/(xs(iret)-xs(irec))*seglen(iret)
        IF (del.LT.delmin) THEN
           delmin = del
           wall   = 'xmax'
        ENDIF
     ENDIF
     IF (ys(iret).LT.ymin) THEN
        del = (ymin-ys(irec))/(ys(iret)-ys(irec))*seglen(iret)
        IF (del.LT.delmin) THEN
           delmin = del
           wall   = 'ymin'
        ENDIF
     ENDIF
     IF (ys(iret).GT.ymax) THEN
        del = (ymax-ys(irec))/(ys(iret)-ys(irec))*seglen(iret)
        IF (del.LT.delmin) THEN
           delmin = del
           wall   = 'ymax'
        ENDIF
     ENDIF
  END IF
  IF (zs(iret).LT.zmin) THEN
     del = (zmin-zs(irec))/(zs(iret)-zs(irec))*seglen(iret)
     IF (del.LT.delmin) THEN
        delmin = del
        wall   = 'zmin'
     ENDIF
  ENDIF
  IF (zs(iret).GT.zmax) THEN
     del = (zmax-zs(irec))/(zs(iret)-zs(irec))*seglen(iret)
     IF (del.LT.delmin) THEN
        delmin = del
        wall   = 'zmax'
     ENDIF
  ENDIF
  
  IF (wall.EQ.'none') RETURN
  ! if we get to here (wall.NE. 'none'), we know we have the tip outside the domain and need
  ! new heading components:
  dx = xs(iret)-xs(irec)
  dy = ys(iret)-ys(irec)
  dz = zs(iret)-zs(irec)
  CALL Compon(wall(1:1),irec,dx,dy,dz)
  IF (delmin.GT.1.E-06_sp) THEN
     ! break up segment --
     length = seglen(iret)
     mass   = segmas(iret)
     fractn = delmin/length
     ! part with unchanged heading:
     segmas(iret) = fractn*mass
     seglen(iret) = delmin
     ! pull back tip to start part with new heading:
     xs(iret) = xs(irec)+fractn*(xs(iret)-xs(irec))
     ys(iret) = ys(irec)+fractn*(ys(iret)-ys(irec))
     zs(iret) = zs(irec)+fractn*(zs(iret)-zs(irec))
     
     ! to prevent numerical problems, pull back exactly to wall
     IF (wall.EQ.'xmin') xs(iret) = xmin 
     IF (wall.EQ.'xmax') xs(iret) = xmax
     IF (wall.EQ.'ymin') ys(iret) = ymin
     IF (wall.EQ.'ymax') ys(iret) = ymax
     IF (wall.EQ.'zmin') zs(iret) = zmin
     IF (wall.EQ.'zmax') zs(iret) = zmax

     ! update tip info
     xg(igrow) = xs(iret)
     yg(igrow) = ys(iret)
     zg(igrow) = zs(iret)

     newlen=length-seglen(iret)
     IF (newlen.GT.1.E-06_dp) THEN
        newtime = timorg(iret)+fractn*dtRoot !convert to sp!
        !toosml(nrec)=.FALSE.
        CALL Grow(igrow,dx,dy,dz,newlen,segmas(iret),newtime)
     ENDIF
  ELSE
     ! don't split, change heading only --
     ! calculate the actual length of components dx,dy,dz:
     factor = seglen(iret)/SQRT(dx*dx+dy*dy+dz*dz)
     dx     = factor*dx
     dy     = factor*dy
     dz     = factor*dz
     ! new position of tip:
     xs(iret) = xs(irec)+dx
     ys(iret) = ys(irec)+dy
     zs(iret) = zs(irec)+dz
     ! update tip info
     xg(igrow) = xs(iret)
     yg(igrow) = ys(iret)
     zg(igrow) = zs(iret)
  END IF
  GOTO 1

END SUBROUTINE Boxlim
 !******************************************************************************
!> for cylindrical soil domains --> growth along cylinder walls
!> if soil geometry is cylindrical, handle like domain edge 
SUBROUTINE Boxlim_Cylinder(igrow) 
  USE typedef
  USE RootData
  USE DomData
  USE tmctrl, ONLY: dtroot
  USE GridData, ONLY:rad_cyl,x_cent,y_cent,dxGrid,dyGrid
  IMPLICIT NONE

  INTEGER(ap), INTENT(in) :: igrow
  INTEGER(ap) :: irec,iret
  REAL(sp):: irad,newtime
  REAL(sp):: factor,fractn,dx,dy,dz
  REAL(sp):: azi,xg_cyl,yg_cyl
  REAL(dp):: length,mass,newmas,delmin,del,newlen
  CHARACTER wall*4
  !> \param igrow number of growing tip

  ! find the segment behind the tip 'igrow':
  ! loop is needed in case of new growth
  ! pull back tip to smallest possible radius
2 iret = irecsg(igrow) !> tip node
  irec = irecpr(iret)  !> prev. node
  IF (toosml(iret)) RETURN
  delmin = 1.E+37_dp
  wall   = 'none'


  ! define orientation of segment and tip node with respect to center, azimuth & quadrant
  irad = sqrt( (xs(iret)-x_cent)*(xs(iret)-x_cent) + (ys(iret)-y_cent)*(ys(iret)-y_cent) )
  azi  = ATAN2(ys(iret)-y_cent,xs(iret)-x_cent)

  ! old heading components
  dx     = xs(iret)-xs(irec)
  dy     = ys(iret)-ys(irec)
  dz     = zs(iret)-zs(irec)
  xg_cyl = xs(iret)
  yg_cyl = ys(iret)
  
  ! check if tip node is within cylinder (x and y)
  IF (irad .GT. (rad_cyl)/sqrt(2._sp)) THEN
     xg_cyl = cos(azi)*(rad_cyl-dxGrid)/sqrt(2._sp)
     yg_cyl = sin(azi)*(rad_cyl-dyGrid)/sqrt(2._sp)    
     wall   = 'cyli'
  END IF
  
  ! check in z
  IF (zs(iret) .LT. zmin) THEN
     del = (zmin-zs(irec))/(zs(iret)-zs(irec))*seglen(iret)
     IF (abs(del) .LT. delmin) THEN
        delmin = del
        wall   = 'zmin'
     END IF
  ENDIF
  IF (zs(iret) .GT. zmax) THEN
     del = (zmax-zs(irec))/(zs(iret)-zs(irec))*seglen(iret)
     IF (abs(del) .LT. delmin) THEN
        delmin = del
        wall   ='zmax'
     END IF
  ENDIF  
  

  IF (wall.EQ.'none') RETURN

  IF (wall.EQ.'cyli') THEN
     xs(iret) = xg_cyl
     ys(iret) = yg_cyl
     dx       = 0._dp
     dy       = 0._dp
     IF(abs(dz).LT.1.E-10) THEN
        IF (ABS(zmax-zs(irec)).GT.ABS(zmin-zs(irec))) THEN
           dz = zmax-zs(irec)
        ELSE
           dz = zmin-zs(irec)
        ENDIF
     END IF
  ELSEIF (wall.EQ.'zmin') THEN
     zs(iret) = zmin
     dz = 0._dp
     IF(abs(dx).LT.1.E-10 .AND. abs(dy).LT.1.E-10) THEN
        dx = x_cent-xs(irec)
        dy = y_cent-ys(irec)
     END IF
  ELSEIF (wall.EQ.'zmax') THEN
     zs(iret) = zmax
     dz = 0._dp
     IF(abs(dx).LT.1.E-10 .AND. abs(dy).LT.1.E-10) THEN
        dx = x_cent-xs(irec)
        dy = y_cent-ys(irec)
     END IF
  END IF

  ! length from segment node to wall
  delmin = sqrt( (xs(iret)-xs(irec))**2 + (ys(iret)-ys(irec))**2 + (zs(iret)-zs(irec))**2 )

  IF (delmin.GT.1.E-06_sp) THEN
     ! break up segment --

     newlen = seglen(iret)-delmin
     mass   = segmas(iret)
     fractn = delmin/seglen(iret) ! length until wall / total length

     ! part with unchanged heading:
     segmas(iret) = (1-fractn)*mass
     seglen(iret) = delmin

     IF (newlen.GT.1.E-06_dp) THEN
        newtime = timorg(irec)+(1-fractn)*dtRoot ! convert to sp
        CALL Grow(igrow,dx,dy,dz,newlen,segmas(iret),newtime)
     else
        ! don´t split, change heading only --
        ! calculate the actual length of components dx,dy,dz:
        factor = seglen(iret)/SQRT(dx*dx+dy*dy+dz*dz)
        dx     = factor*dx
        dy     = factor*dy
        dz     = factor*dz
        
        ! new position of tip:
        xs(iret) = xs(irec)+dx
        ys(iret) = ys(irec)+dy
        zs(iret) = zs(irec)+dz

        ! update tip info
        xg(igrow) = xs(iret)
        yg(igrow) = ys(iret)
        zg(igrow) = zs(iret)
       
     END IF
  ELSE   
     ! don´t split, change heading only --
     ! calculate the actual length of components dx,dy,dz:
     factor = seglen(iret)/SQRT(dx*dx+dy*dy+dz*dz)
     dx     = factor*dx
     dy     = factor*dy
     dz     = factor*dz
     
     ! new position of tip:
     xs(iret) = xs(irec)+dx
     ys(iret) = ys(irec)+dy
     zs(iret) = zs(irec)+dz
     
     ! update tip info
     xg(igrow) = xs(iret)
     yg(igrow) = ys(iret)
     zg(igrow) = zs(iret)

  END IF
  GOTO 2
END SUBROUTINE Boxlim_Cylinder
!******************************************************************************
!> reorients the part of the segment that was outside the domain back into the domain
SUBROUTINE Compon(wall,irec,dx,dy,dz)
  USE typedef
  USE RootData
  USE DomData
  IMPLICIT NONE

  INTEGER(ap),INTENT(in):: irec
  REAL(sp), INTENT(inout):: dx,dy,dz
  CHARACTER wall*1
  !> \param wall intersection side
  !> \param irec number of segment
  !> \param dx x component of growing tip
  !> \param dy y component of growing tip
  !> \param dz z component of growing tip

  IF (wall.EQ.'x') THEN
     dx = 0.0_sp
     IF ((ABS(dy).LT.1.E-10_dp).AND.(ABS(dz).LT.1.E-10_dp)) THEN
        IF (ABS(ymax-ys(irec)).GT.ABS(ymin-ys(irec))) THEN
           dy = ymax-ys(irec)
        ELSE
           dy = ymin-ys(irec)
        ENDIF
        IF (ABS(zmax-zs(irec)).GT.ABS(zmin-zs(irec))) THEN
           dz = zmax-zs(irec)
        ELSE
           dz = zmin-zs(irec)
        ENDIF
     ENDIF
  ELSEIF(wall.EQ.'y') THEN
     dy = 0.0_sp
     IF ((ABS(dx).LT.1.E-10_dp).AND.(ABS(dz).LT.1.E-10_dp)) THEN
        IF (ABS(xmax-xs(irec)).GT.ABS(xmin-xs(irec))) THEN
           dx = xmax-xs(irec)
        ELSE
           dx = xmin-xs(irec)
        ENDIF
        IF (ABS(zmax-zs(irec)).GT.ABS(zmin-zs(irec))) THEN
           dz = zmax-zs(irec)
        ELSE
           dz = zmin-zs(irec)
        ENDIF
     ENDIF
  ELSEIF(wall.EQ.'z') THEN
     dz = 0.0_sp
     IF ((ABS(dx).LT.1.E-10_dp).AND.(ABS(dy).LT.1.E-10_dp)) THEN
        IF (ABS(xmax-xs(irec)).GT.ABS(xmin-xs(irec))) THEN
           dx = xmax-xs(irec)
        ELSE
           dx = xmin-xs(irec)
        ENDIF
        IF (ABS(ymax-ys(irec)).GT.ABS(ymin-ys(irec))) THEN
           dy = ymax-ys(irec)
        ELSE
           dy = ymin-ys(irec)
        ENDIF
     ENDIF
  ENDIF

  RETURN
END SUBROUTINE Compon
!******************************************************************************
!> remove all segments from list that are too small 
SUBROUTINE Remove(nrecol)
  USE typedef
  USE RootData
  IMPLICIT NONE

  INTEGER(ap) :: nrecnw,nrecol,ifrom,ito,igrow
  !> \param nrecol total number of segment at the beginning of the call to ROOT

  ! set nrecnw equal to the current number of segment records
  nrecnw = nrec
  ito    = 0
  ifrom  = nrecol
10 ifrom = ifrom+1
  IF (toosml(ifrom)) THEN
     ! identify and pull back tip that belongs to the record being removed:
     igrow = 0
11   igrow = igrow+1
     IF (ibrgrw(igrow).NE.ibrseg(ifrom)) GOTO 11
     xg(igrow)     = xs(ifrom)
     yg(igrow)     = ys(ifrom)
     zg(igrow)     = zs(ifrom)
     irecsg(igrow) = irecpr(ifrom)
     ! decrease total number of records by '1':
     nrec = nrec-1
     ! if this the first candidate for removal (to be overwritten),
     ! mark position:
     IF (ito.EQ.0) ito = ifrom
  ELSE IF (ito.GT.0) THEN
     ! can move from 'ifrom' because at least one removal candidate
     ! has been previously identified and marked as 'ito':
     toosml(ito) = toosml(ifrom)
     xs(ito)     = xs(ifrom)
     ys(ito)     = ys(ifrom)
     zs(ito)     = zs(ifrom)
     irecpr(ito) = irecpr(ifrom)
     ordseg(ito) = ordseg(ifrom)
     ibrseg(ito) = ibrseg(ifrom)
     seglen(ito) = seglen(ifrom)
     segmas(ito) = segmas(ifrom)
     timorg(ito) = timorg(ifrom)
     ! tip that belongs to moved segment needs the new segment reference, ito:
     igrow = 0
12   igrow = igrow+1
     IF (ibrgrw(igrow).NE.ibrseg(ito)) GOTO 12
     irecsg(igrow) = ito
     ito = ito+1
  ENDIF
  IF (ifrom.LT.nrecnw) GOTO 10
  ! now nrec is the actual number of segment records
  RETURN
END SUBROUTINE Remove
!*********************************************************************************
!> Update root records after removal of segments 
SUBROUTINE Update(ngrwnw)
  USE typedef
  USE RootData
  USE GridData, ONLY: continu
  USE DoussanMat, ONLY: transroot
  IMPLICIT NONE

  INTEGER(ap):: ito,ifrom,iest,ngrwnw
  !> \param ngrwnw number of growing tips

  ! may have more growing branch tips at the end of time step:
  ngrow  = ngrwnw
  ito    = 0
  ifrom  = 0
10 ifrom = ifrom+1
  IF (stopgr(ifrom)) THEN
     ngrow = ngrow-1
     ! if this the first candidate for removal (to be overwritten),
     ! mark position:
     IF (ito.EQ.0) ito = ifrom
  ELSE IF (ito.GT.0) THEN
     ! can move from 'ifrom' because at least one removal candidate
     ! has been previously identified and marked as 'ito':
     stopgr(ito)  = stopgr(ifrom)
     xg(ito)      = xg(ifrom)
     yg(ito)      = yg(ifrom)
     zg(ito)      = zg(ifrom)
     iaxis(ito)   = iaxis(ifrom)
     irecsg(ito)  = irecsg(ifrom)
     ordgrw(ito)  = ordgrw(ifrom)
     ibrgrw(ito)  = ibrgrw(ifrom)
     brlgth(ito)  = brlgth(ifrom)
     ovrtime(ito) = ovrtime(ifrom)
     nestbl(ito)  = nestbl(ifrom)
     DO iest=1,nestbl(ito)
        timest(ito,iest) = timest(ifrom,iest)
     END DO
     IF(continu) transroot(ito,1:2,1,1) = transroot(ifrom,1:2,1,1)
     ito = ito+1
  ENDIF
  IF (ifrom.LT.ngrwnw) GOTO 10
  ! now ngrow is the actual number of growing tips at beginning of next step
  RETURN
END SUBROUTINE Update
!*********************************************************************************
!> additional secondary radial growth
SUBROUTINE Secondary_Radial_Growth(dt)
  USE typedef
  USE RootData, ONLY: nrec,segsur,f_rad,ordseg
  IMPLICIT NONE
  REAL(sp),INTENT(in) :: dt
  INTEGER(ap):: irec

  DO irec = 1,nrec
     segsur(irec) = segsur(irec) + segsur(irec)*(f_rad(ordseg(irec))-1)*dt
  END DO
  
END SUBROUTINE Secondary_Radial_Growth
!*********************************************************************************
