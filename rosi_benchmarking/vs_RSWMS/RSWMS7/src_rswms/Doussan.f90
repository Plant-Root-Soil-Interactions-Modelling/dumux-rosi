! ==============================================================================
! Source file DOUSSAN |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
!> calculates where roots intersect which plane of a soil voxel
LOGICAL FUNCTION intsec(xA,yA,zA,xB,yB,zB,ifoln,irecn,isub,ipl,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,iFace)
  USE Typedef
  USE GridData, ONLY: xgrid,ygrid,dxgrid,dygrid,nex,ney,continu
  USE DoussanMat, ONLY: transroot,nsub,isubmax
  IMPLICIT NONE
  REAL (sp) ::xA,yA,zA,xB,xBp,yB,yBp,zB,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,xc
  REAL(sp) ::yc,zc,f
  INTEGER(ap) ::iFace,irecn,ifoln,isub,ipl
  intsec=.FALSE.
  IF (continu) THEN										
     !> if continous domain, reconstruction of the continuous B-Ap segment
     IF (isub.EQ.1) THEN
        xBp=xB+(transroot(ifoln,1,nsub(ifoln,ipl),ipl)-transroot(irecn,1,1,ipl))*(nex*dxgrid)	!ifoln refers to the "Ap" node and irecn to the "B" node
        yBp=yB+(transroot(ifoln,2,nsub(ifoln,ipl),ipl)-transroot(irecn,2,1,ipl))*(ney*dygrid)
     ELSE
        xBp=xB+(transroot(irecn,1,isub,ipl)-transroot(irecn,1,1,ipl))*(nex*dxgrid)
        yBp=yB+(transroot(irecn,2,isub,ipl)-transroot(irecn,2,1,ipl))*(ney*dygrid)
     ENDIF
  ELSE
     xBp=xB
     yBp=yB
  ENDIF
  
  IF (zB.LT.z1) THEN
     !> may have intersection with x-y-plane of cube at z1
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
     !> may have intersection with x-y-plane of cube at z2
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
     !> may have intersection with x-z-plane of cube at y1
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
     !> may have intersection with x-z-plane of cube at y2
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
     !> may have intersection with y-z-plane of cube at x1
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
     !> may have intersection with y-z-plane of cube at x2
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
!> Doussan model implementation through 3 Subroutines
!>	- SetupDou
!>       - segment
!>	- SetBCroot
!>       - conductroot
!>	- SolveRoot
!> and the adaptation of SetSink
!===============================================================================
SUBROUTINE SetupDou(t,dt)
  USE TypeDef
  USE Paramdata,ONLY: ldirect,lChem,lretry,last_out,maxemg 
  USE SparseMatrix
  USE DoussanMat
  USE GridData, ONLY:xgrid,ygrid,zgrid,dzgrid,dygrid,dxgrid,nex,ney,nElm,continu,nPt,betac,elmnod,subN,iL
  USE RootData, ONLY: lCalloc,loliumroottyp,maizeroottyp,wheatroottyp,nbr,segsur,ibrseg,seglen,xplant,yplant,ibrgrw,zs,ys,xs,&
		zg,yg,xg,ordseg,timorg,irecpr,nurf,age,urf,norder,naxes,lDou,lFed,lCou,lno_Archi,crossSectionSeg,irecsg
  USE tmctrl, ONLY: t_begin,nOut,tout
  USE PlntData, ONLY: TotSur,Tpot
  IMPLICIT NONE
  INTEGER(ap) ::ibr,irecn,ifoln,igrow,iprvn,iBCn,typ
  INTEGER(ap) ::isub,i,j,corner(1:8),iUrf,iseg,iorder,k,l,iE,iSE
  INTEGER(ap) :: err,ipl,ibc 
  REAL(sp) ::t,xA,yA,zA,xB,yB,zB,x1,x2,y1,y2,z1,z2,dt,Weight,betce,VEl,Sbetac, segage
  REAL(dp) ::PrvSur
  REAL(sp) :: y(1:nLibr),h_mat(1:nLibr),arr(1:nLibr-1),cumsum(1:nLibr-1),PHtemp
  REAL(sp) :: theta(1:nLibr),Kh(1:nLibr),C(1:nLibr),temp_phr
  REAL(sp) :: rs,concrs
  LOGICAL :: n_apex,run
      
  WRITE (*,*) '... Analysing <RootSys> ...'    
  IF (lDou.OR.(lCou.AND.(.NOT.lno_Archi))) THEN  
     ALLOCATE (nBC_irecn(1:maxemg))
     nBC_irecn=0
  
     !> initialize the pointer of the first el of the list
     !> initialize matrices
     CALL IniMat
     iBCn=0

     !> shift rootsystem for continous domain - needed if macropores for ConductRoot
     IF(continu) THEN
        DO irecn=1,nrec
           CALL roottrans(xs(irecn),ys(irecn),irecn,1,1)
        END DO
        DO igrow=1,ngrow
           CALL tiptrans(xg(igrow),yg(igrow),igrow,1,1)
        END DO
        transroot(nrec+1:nrec+ngrow,:,1,1)=transtip(1:ngrow,:,1,1)
     END IF

     !> Current axial and radial conductance matrix
     CALL ConductRoot(t)
  
     !> prepare calculation of betaw and betac	
     ! IF (lChem) THEN
     PrvSur=0.0_dp
     DO iseg=1,nrec+1
        PrvSur=PrvSur+segsur(iseg)
     ENDDO
     IF (PrvSur.LT.1.E-20_dp) PrvSur=1.E-20_dp
     TotSur = PrvSur
     betac=0.0_dp
  
     !ENDIF

     !> go through different plants
  ELSEIF (lFed.AND.(.NOT.lno_Archi)) THEN	!For Doussan, that is done in IniMat
     ALLOCATE (transroot(0:nrec+ngrow,1:2,1:isubmax,1:nplant))
     transroot=0
     ALLOCATE (nsub(0:nrec+ngrow,1:nplant))
     nsub=0
     ALLOCATE (cube_i(0:nrec,1:isubmax,1:nplant))
     cube_i=0
     ALLOCATE (Intc(0:nrec+ngrow,1:3,1:isubmax,1:nplant))
     Intc=0._dp
     ALLOCATE (l_seg(0:nrec))					!actually not used
     l_seg=0._dp
     ALLOCATE (beta_weight(0:nrec,1:isubmax,1:nplant))		!actually not used
     beta_weight=0._dp
     ALLOCATE (w_sub(0:nrec,1:isubmax,1:nplant))			!actually not used
     w_sub=0._dp
     ALLOCATE (l_sub(0:nrec,1:isubmax,1:nplant))			!actually not used
     l_sub=0._dp
     ALLOCATE (cent(0:nrec,1:3,1:isubmax,1:nplant))		!actually not used
     cent=0._dp
     ALLOCATE (sum_dis(0:nrec,1:isubmax,1:nplant))
     sum_dis=0._dp
     ALLOCATE (w_dis(0:nrec,1:8,1:isubmax,1:nplant))
     w_dis=0._dp
     betaw=0._sp
  ENDIF

  IF (lretry) THEN
     t=tOut(last_out)
     WRITE (*,'(///'' Simulation continues at time = '',F8.3,''.'')') t
  ENDIF
  
  DO ipl=1,nplant
     !> Current boundary conditions for root
     IF (lDou) THEN
        IF (.NOT.lCalloc) THEN
           CALL SetBCroot(t,curr_BCr(ipl),curr_BCtp(ipl))
        ELSE
           IF (ipl.GE.2) STOP 'Assimilate allocation currently doesnt work with multiple plants'
           CALL Stress(rs,concrs)
           CALL SetTp(t,rs,concrs,ipl)
           curr_BCr(ipl)=BCr_usr(ipl)     ! -ABS(Tpot) 
           curr_BCtp(ipl)=2
        ENDIF
        !multiple roots
        err=SM_allocate(plantmatrix(ipl), nrec+1, nrec+1)
        IF(err/=0) STOP 'Could not create plantmatrix'
     ENDIF
     !> go through each root segment and update the node surface function:
     DO ibr=1,nbr ! all plants have same number of roots
        n_apex=.FALSE.
        !> find the tip segment of the branch 'ibr'
        irecn = irecsg(ibr)

        !> the first one we find is an apex
        IF (seglen(irecn)<1.E-20) THEN !> skip this segment too small to be taken
           xA=xs(irecn)+xplant(ipl)
           yA=ys(irecn)+yplant(ipl)
           zA=zs(irecn)
           IF (continu) THEN
              xA=xA+transroot(irecn,1,1,ipl)*dxgrid*nex
              yA=yA+transroot(irecn,2,1,ipl)*dygrid*ney
           END IF
           
           ifoln=irecn !> following node ID
           irecn=irecpr(irecn) !> current node ID

        ELSE !> segment long enough

           n_apex=.TRUE.!> ifoln does not exist if not continu
           DO igrow=1,ngrow
              IF (ibrgrw(igrow)==ibr) THEN !> apex of branch ibr
                 xA=xg(igrow)+xplant(ipl)
                 yA=yg(igrow)+yplant(ipl)
                 zA=zg(igrow)
                 IF (continu) THEN
                    xA=xA+transroot(nrec+igrow,1,1,ipl)*dxgrid*nex
                    yA=yA+transroot(nrec+igrow,2,1,ipl)*dygrid*ney
                    ifoln=nrec+igrow                             
                    nsub(nrec+igrow,ipl)=1
                 END IF
              ENDIF
           END DO
        ENDIF
        
        IF (irecn==0) THEN !> there exists a branch ibr but not yet any segment!
           run=.FALSE.
        ELSE
           run=.TRUE.
        ENDIF
        
        !> then the rest of the branch up to the seed of the embranchment
        DO WHILE (run)
           !> "upper" node
           iprvn=irecpr(irecn)
           !> location of the rear or "upper" end
           xB=xs(irecn)+xplant(ipl)
           yB=ys(irecn)+yplant(ipl)
           zB=zs(irecn)
           IF (continu) THEN
              xB=xB+transroot(irecn,1,1,ipl)*dxgrid*nex
              yB=yB+transroot(irecn,2,1,ipl)*dygrid*ney
           END IF
           IF (lDou.OR.(lCou.AND.(.NOT.lno_Archi))) THEN
              !> calculate the gravity components z (always positive & maximum at soil surface)
              GH(irecn)=(zA+zB)/2.
           ENDIF
           !> calculate segment weighing factor according to age:
           iorder=ordseg(irecn)
           IF (maizeroottyp) THEN
              IF (iorder.LT.12) THEN  		!In RootTyp, maize root types below 12 are principal roots (12 and 13 are lateral)
                 typ=1
              ELSE
                 typ=2
              ENDIF
           ELSEIF (loliumroottyp) THEN
              IF (iorder.LT.3) THEN  		!In RootTyp, lolium root types below 3 are principal roots
                 typ=1
              ELSEIF (iorder.LT.5) THEN
                 typ=2
              ELSE
                 typ=3
              ENDIF
           ELSEIF (wheatroottyp) THEN
              IF (iorder.LT.19) THEN  		!In RootTyp, wheat root types below 19 are principal roots
                 typ=1
              ELSEIF (iorder.LT.20) THEN
                 typ=2
              ELSE
                 typ=3
              ENDIF
           ELSE
              IF (iorder.EQ.0) THEN  !type 0 is the seed+small segment
                 typ=1
              ELSEIF (iorder.GT.3) THEN!no more than 3 root types
                 typ=3
              ELSE
                 typ=iorder
              ENDIF
           ENDIF
           
           segage=t-timorg(irecn)
           IF (lChem.and.(segage.GE.0.0_dp)) THEN
              IF (segage.GE.age(typ,nUrf(typ))) THEN
                 Weight=Urf(typ,nUrf(typ))
              ELSE
                 iUrf=nUrf(typ)
                 iUrf=iUrf-1
                 DO WHILE ((segage.LT.age(typ,iUrf)).AND.(iUrf.GT.1))
                    iUrf=iUrf-1
                 ENDDO
                 Weight=Urf(typ,iUrf)+(segage-age(typ,iUrf))/&
                      (age(typ,iUrf+1)-age(typ,iUrf))*&
                      (Urf(typ,iUrf+1)-Urf(typ,iUrf))
              ENDIF
           ELSE
              Weight=1.0_dp
           ENDIF

           !> calculate number of subsegment for node/segment irecn and corresponding weights
           CALL segment(xA,yA,zA,xB,yB,zB,ifoln,irecn,ipl,Weight,PrvSur)

           IF (lDou) THEN
              !> several changes -> multiple roots
              err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,-Khr(irecn)/seglen(irecn),.FALSE.)
              IF(err/=0) STOP 'Could not insert element into plantmatrix'
              !> if apex (bottom part of root)
              IF (n_apex) THEN
                 Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn)/seglen(irecn)+Lr(irecn)*segsur(irecn),.FALSE.)
                 IF(err/=0) STOP 'Could not insert element into plantmatrix'
              ELSE
                 Err=SM_set(Plantmatrix(ipl), irecn+1,ifoln+1,-Khr(ifoln)/seglen(ifoln),.FALSE.) !row, col,value
                 Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn)/seglen(irecn)+Khr(ifoln)/seglen(ifoln)+Lr(irecn)*segsur(irecn),.FALSE.)
                 IF(err/=0) STOP 'Could not insert element into plantmatrix'
              ENDIF

              !> define 1st part of Q (Q=Qi.*PHsoil+Qbc) -> RHS
              Qi(irecn,ipl)=Lr(irecn)*segsur(irecn)
              !if(irecn.lt.100) print*,'qi',irecn,qi(irecn,ipl)
           ENDIF
           !> if reached the seed or the embranchement => change branch
           IF (iprvn==0) THEN!seed=first segment
              IF (lDou) THEN
                 GH(0)=zB
                 iBCn=iBCn+1
                 nBC_irecn(iBCn)=irecn
                 IF (curr_BCtp(ipl)==2.and..not.(ldirect)) THEN!flux
                    Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn)/seglen(irecn),.FALSE.)
                    Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn)/seglen(irecn),.FALSE.)
                    IF(err/=0) STOP 'Could not insert element into plantmatrix'
                    !> own additive; position (2,2) changes
                    !> own additive;position (2,1) = 0, changes into zero
                    !> rhs has only 1 entry (first one)
                    Qi(iprvn,ipl)=0._dp
                    Q_bc1(iprvn,ipl)=curr_BCr(ipl) !Q(0,0)
                    !                     Q_bc(irecn,ipl)=curr_BCr
                 ELSE IF ((curr_BCtp(ipl)==2).and.(ldirect)) THEN  !flux + DIRECT
                    Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,1._dp,.TRUE.)
                    Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,0._dp,.TRUE.)
                    Err=SM_set(Plantmatrix(ipl), irecn+1,iprvn+1,0._dp,.TRUE.)
                    IF(err/=0) STOP 'Could not insert element into plantmatrix'
                    !own additive; position (2,2) changes 
                    !own additive;position (2,1) = 0, changes into zero
                    !rhs has only 1 entry (first one)
                    Qi(iprvn,ipl)=0._dp
                    !Q_bc1(iprvn,ipl)=curr_BCr(ipl)!Q(0,0)
                    !                     Q_bc(irecn,ipl)=curr_BCr
                    DO iBC=1,1  !> all the branches connected to the seed
                       irecn=nBC_irecn(iBC)
                       temp_phr=-(curr_BCr(ipl)/Khr(irecn)*seglen(irecn+1)+200.)

                       Qi(iprvn,ipl)=0._dp !> from seed has to be given, from irecn is already stated in DoussanMat
                       Q_bc1(iprvn,ipl)=temp_Phr*.1+GH(iprvn) !> use ph from last time step
                       Q_bc2(iprvn,ipl)=temp_Phr*2.+GH(iprvn)
                       Q_bc1(irecn,ipl)=(temp_Phr*.1+GH(iprvn))*Khr(irecn)/seglen(irecn) !> all the PH for root was given in total head!!
                       Q_bc2(irecn,ipl)=(temp_Phr*2+GH(iprvn))*Khr(irecn)/seglen(irecn)
                    ENDDO
                 ELSE IF (curr_BCtp(ipl)==1) THEN !PH
                    !iprvn+1,iprvn+1
                    Qi(iprvn,ipl)=0._dp
                    !> true means overwrite (delete first), no addition
                    !iprvn+1,iprvn+1
                    Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,1._dp,.TRUE.) !position (1,1)
                    Err=SM_set(Plantmatrix(ipl), irecn+1,iprvn+1,0._dp,.TRUE.) 
                    IF(err/=0) STOP 'Could not insert element into plantmatrix'
                    !> first entry in rhs is PH (inlc. gravity)
                    Q_bc1(iprvn,ipl)=curr_BCr(ipl)+GH(iprvn)!Q(iprv)-BCr*Kh(irecn)/seglen(irecn);
                    !> second entry is PH incl. gravity times these parameters
                    Q_bc1(irecn,ipl)=(curr_BCr(ipl)+GH(iprvn))*Khr(irecn)/seglen(irecn)
                    Q_bc2 = 0._dp
                 ENDIF
              ENDIF

              run=.FALSE.

           ELSEIF (ibrseg(iprvn).NE.ibrseg(irecn)) THEN !> start of the branch but not from the seed -> gaat van onder na boven, duz iprvn (nieuwe positie) is nu niet van die ene branch maar van een side branch
              IF (lDou) THEN
                 Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn)/seglen(irecn),.FALSE.)
                 Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn)/seglen(irecn),.FALSE.)
                 IF(err/=0) STOP 'Could not insert element into plantmatrix'
              ENDIF
              run=.FALSE.
           ENDIF

           !> definition for the next run of the loop
           ifoln=irecn
           irecn=iprvn
           !> location of the final node
           xA=xB	!> A directly redefined with report to B in order not to have to calculate all the translations of A to the inside of the domain
           yA=yB 
           zA=zB
           !> from here, not an apex
           n_apex=.FALSE.
        END DO !> loop on branch nodes
     END DO !> loop on root branches

     IF (lDou.OR.(lCou.AND.(.NOT.lno_Archi))) THEN
        IF (.NOT.(old)) THEN
           IF ((ave) .OR. (eqDis)) THEN
              DO i=1,nElm !> no_voxels -> number of nodes in  a cuboid (imin)
                 IF (no_voxels(i) .EQ. 0) GOTO 40
                 DO j=1,no_voxels(i)
                    irecn = voxel_node(i,1,j)
                    isub = voxel_node(i,2,j)
                    corner=loc_Q(irecn,1:8,isub,ipl)
                    !coordinates of voxel
                    x1=xgrid(corner(1))
                    y1=ygrid(corner(1))
                    z1=zgrid(corner(1))
                    x2=xgrid(corner(8))
                    y2=ygrid(corner(8))
                    z2=zgrid(corner(8))
                    IF (ave) THEN
                       cp_mean(irecn,1,isub,ipl)=(x2+x1)/2
                       cp_mean(irecn,2,isub,ipl)=(y2+y1)/2
                       cp_mean(irecn,3,isub,ipl)=(z2+z1)/2
                    ELSEIF (eqDis) THEN
                       numNodes_voxel(irecn,isub)=no_voxels(i) !> total root nodes of cuboid in no_voxels;  cubiods have same imin
                    ENDIF
                 ENDDO
 40              CONTINUE
              ENDDO
           ENDIF
           WRITE(*,*)'t=',t
           WRITE(*,*)'t_begin=',t_begin,t_begin+dt
           WRITE(*,*)'old=',old
           !if (.not.(old)) then
           IF (t .LE. t_begin+dt) THEN
              !call cpu_time(t0)
              WRITE(*,*)'Matric flux potential Library loaded'
              !-----------------------------------------------------------------------
              !> Matric flux potential lookup table -> to obtain corresponding PH
              !------------------------------------------------------------------------
              PHtemp = 1e-5
              CALL logspace(hx_min,PHtemp,nLibr,y)
              h_mat = -y
              h_mat2=(h_mat(2:SIZE(h_mat))+h_mat(1:SIZE(h_mat)-1))/2
              CALL setmat_anaLibr(h_mat,theta,Kh,C)
              !> calculate matric flux potential integral K(h)dh -> numerical integration -> to linearize non-linear conductivity
              arr=ABS(h_mat(2:SIZE(h_mat))-h_mat(1:SIZE(h_mat)-1))*(Kh(2:SIZE(Kh))+Kh(1:SIZE(Kh)-1))/2
              cumsum(SIZE(arr)) = SUM(arr(1:SIZE(arr)))
              Phi_mat = cumsum
           ENDIF
        ENDIF
        !> total number of embrach. to the seed
        nBCn=iBCn
     ENDIF
  END DO !> loop on plants

  IF (lChem) THEN		
     Sbetac=0.0_dp
     VEl=dxGrid*dyGrid*dzGrid/6._sp
     DO iE=1,nElm
        DO iSE=1,5
           i=elmnod(iL(1,iSE,subN(iE)),iE) !> i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
           j=elmnod(iL(2,iSE,subN(iE)),iE)
           k=elmnod(iL(3,iSE,subN(iE)),iE)
           l=elmnod(iL(4,iSE,subN(iE)),iE)
           betcE=(betac(i)+betac(j)+betac(k)+betac(l))/4.
           Sbetac=Sbetac+betcE
        ENDDO
     ENDDO
     IF (Sbetac.GT.1.E-20_dp) THEN
        Sbetac=Sbetac*VEl
        betac=betac/Sbetac
     ELSE
        WRITE(*,*)'Sbetac < 10^-20,  Sbetac = ',Sbetac
     ENDIF
  ENDIF

 END SUBROUTINE SetupDou
 !*******************************************************************************
 SUBROUTINE segment(xA,yA,zA,xB,yB,zB,ifoln,irecn,ipl,Weight,PrvSur)
   USE Typedef
   USE ParamData, ONLY: lChem
   USE GridData, ONLY:xgrid,ygrid,zgrid,dzgrid,dygrid,dxgrid,nex,ney,continu,betac,betaw
   USE DoussanMat, ONLY: nsub,loc_Q,w_sub,sum_dis,w_dis,cent,Intc,no_voxels,voxel_no,&
        voxel_node,indexValue,old,ave,eqdis,l_seg,l_sub,isubmax,transroot,cube_i, beta_weight,Lr
   USE RootData, ONLY: seglen,ibrseg,irecpr,segsur,lDou,lCou,lFed,lno_Archi
   USE SolData, ONLY: nMat,matNum,par,lMacro

   IMPLICIT NONE
   REAL(sp),INTENT(in) :: xA,yA,zA,xB,yB,zB
   INTEGER(ap),INTENT(in) :: ipl
   REAL(sp):: x1,x2,y1,y2,z1,z2,xAp,yAp,zAp
   REAL(sp) ::yint,xint,zint,xCent,yCent,zCent,weight,dummat
   REAL(dp)::blengt,blengtTot,srface,PrvSur
   LOGICAL :: SPLIT,intsec,lMix=.FALSE.
   INTEGER(ap) ::iface,iFaceOld,isub,corner(8),ic,ifoln,irecn,reftransroot(2),imin,iMat

   !> number of subsegment for node irecn
   !IF(continu) nsub(irecn,ipl)=0
   nsub(irecn,ipl)=0
   IF (lDou.OR.((lCou.OR.lFed).AND.(.NOT.lno_Archi))) THEN
      l_seg(irecn)=0
      blengtTot=0.00000001_dp
   ENDIF
   SPLIT=.TRUE.
   !> initialization
   xAp=xA
   yAp=yA
   zAp=zA
   IF(continu) THEN
      IF(ifoln.NE.0) THEN
         reftransroot=transroot(ifoln,:,nsub(ifoln,ipl),ipl)	!> in Intc, transroot, cube_i, etc., the node from rootsys (-> not created by instec) is written in the last columns
      ELSE
         reftransroot=transroot(ifoln,:,1,ipl)
      ENDIF
   ENDIF

   splitloop:   DO WHILE(SPLIT)
      nsub(irecn,ipl)=nsub(irecn,ipl)+1
      isub=nsub(irecn,ipl)
      IF (isub.GT.isubmax) THEN
         PRINT *,'Number of subsegments higher than maximum admitted, segment',irecn,'too long with report to the grid resolution'
         STOP
      ENDIF
      !  find cuboid around the apical node of the segment:
      CALL Neighb(xAp,yAp,zAp,corner,imin)
      cube_i(irecn,isub,ipl)=INT(imin)
      IF (lDou.OR.(lCou.AND.(.NOT.lno_Archi))) THEN
         Loc_q(irecn,1:8,isub,ipl)=corner
         IF (.NOT.(old)) THEN
            IF ((ave) .OR. (eqdis)) THEN
               IF (voxel_no(irecn,isub) .EQ. 0) THEN !> remove double values from different branchings; if this root node is not processed process it
                  voxel_no(irecn,isub)=imin !> voxel number
                  no_voxels(imin)=no_voxels(imin)+1 !> number of root nodes in a voxel
                  voxel_node(imin,1,no_voxels(imin))=irecn
                  voxel_node(imin,2,no_voxels(imin))=isub
                  IF (no_voxels(imin) .GT. indexValue) THEN
                     PRINT*,'Parameter indexValue in Segment, representing number of nodes in a voxel, has to be set larger'
                     STOP
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      !> calculate cuboid´s corners coordinates	
      !> Valid in all cases (continu or not) 
      x1=xgrid(corner(5))
      y1=ygrid(corner(5))
      z1=zgrid(corner(5))
      x2=x1+dxgrid					
      y2=y1+dygrid
      z2=z1+dzgrid



      !> check if segment is completely included in the cuboid
      IF (intsec(xAp,yAp,zAp,xB,yB,zB,ifoln,irecn,isub,ipl,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,iFace)) THEN !> true when the root segment intersects the cuboïd
         SPLIT=.TRUE.
      ELSE
         SPLIT=.FALSE. !> just/still 1 loop and then no subsegmentation
         xInt=xB
         yInt=yB
         zInt=zB
      ENDIF
      !> correction of Ap position 
      IF (isub.GT.1) THEN
         SELECT CASE(iFaceOld)
         CASE(1,3,5) 
            zAp=zAp+5.E-5_dp*dzGrid
         CASE(2,4,6) 
            zAp=zAp-5.E-5_dp*dzGrid
         END SELECT
      ELSEIF (continu) THEN	
         !> The "downward" extremity of a segment is the following node while the other (downward) parts belong to the current node 
         Intc(ifoln,1,nsub(ifoln,ipl),ipl)=xAp
         Intc(ifoln,2,nsub(ifoln,ipl),ipl)=yAp
         Intc(ifoln,3,nsub(ifoln,ipl),ipl)=zAp
      ENDIF

      IF (lDou.OR.((lCou.OR.lFed).AND.(.NOT.lno_Archi))) THEN
         !> calculate (sub-)segment length and surface
         blengt=SQRT((xInt-xAp)*(xInt-xAp)+(yInt-yAp)*(yInt-yAp)+ (zInt-zAp)*(zInt-zAp))
         blengtTot=blengtTot+blengt
         srface = 0.0_dp
         !TotSur = 0.0_sp
         IF (irecn.GT.0) THEN 
            srface=blengt*segsur(irecn)/seglen(irecn)
            beta_weight(irecn,isub,ipl)=srface/Prvsur
         ENDIF


         !> calculate relative length of this (sub)segment
         IF (irecn.EQ.0) THEN
            w_sub(irecn,isub,ipl)=blengt/seglen(irecn+1)
         ELSE
            w_sub(irecn,isub,ipl)=blengt/seglen(irecn)
         ENDIF
         IF ((.NOT.(split)).AND.(isub.EQ.1.)) THEN
            w_sub(irecn,isub,ipl)=1._dp
         ENDIF
         IF ((w_sub(irecn,isub,ipl).GT.(1.-1.E-7)).OR.((.NOT.(split)).AND.(isub.EQ.1.))) THEN
            w_sub(irecn,isub,ipl)=1._dp
         ELSEIF (w_sub(irecn,isub,ipl).LT.1.E-7) THEN
            w_sub(irecn,isub,ipl)=0._dp
         ENDIF
         l_sub(irecn,isub,ipl)=blengt
         l_seg(irecn)=l_seg(irecn)+blengt

         !> calculate (sub)segment center coordinates...
         xCent=xAp+(xInt-xAp)/2.
         yCent=yAp+(yInt-yAp)/2.
         zCent=zAp+(zInt-zAp)/2.
         cent(irecn,1,isub,ipl)=xCent
         cent(irecn,2,isub,ipl)=yCent
         cent(irecn,3,isub,ipl)=zCent

         !> and calculate the distribution for each node:
         sum_dis(irecn,isub,ipl)=0.0_dp
         DO ic=0,1 							!Valid in all cases (continu or not) 
            w_dis(irecn,4*ic+1,isub,ipl)=SQRT((xCent-xgrid(corner(4*ic+1)))**2+(yCent-ygrid(corner(4*ic+1)))**2+(zCent-zgrid(corner(4*ic+1)))**2)
            w_dis(irecn,4*ic+2,isub,ipl)=SQRT((xCent-xgrid(corner(4*ic+1))-dxgrid)**2+(yCent-ygrid(corner(4*ic+1)))**2+(zCent-zgrid(corner(4*ic+1)))**2)
            w_dis(irecn,4*ic+3,isub,ipl)=SQRT((xCent-xgrid(corner(4*ic+1)))**2+(yCent-ygrid(corner(4*ic+1))-dygrid)**2+(zCent-zgrid(corner(4*ic+1)))**2)
            w_dis(irecn,4*ic+4,isub,ipl)=SQRT((xCent-xgrid(corner(4*ic+1))-dxgrid)**2+(yCent-ygrid(corner(4*ic+1))-dygrid)**2+(zCent-zgrid(corner(4*ic+1)))**2)
         END DO

         DO ic=1,8
            IF (w_dis(irecn,ic,isub,ipl).LT.1.E-20_dp) w_dis(irecn,ic,isub,ipl)=1.E-20_dp
            w_dis(irecn,ic,isub,ipl)=1._dp/w_dis(irecn,ic,isub,ipl)
            sum_dis(irecn,isub,ipl)=sum_dis(irecn,isub,ipl)+w_dis(irecn,ic,isub,ipl)
            IF (lFed) betaw(corner(ic))=betaw(corner(ic))+w_dis(irecn,ic,isub,ipl)/sum_dis(irecn,isub,ipl)*Weight*(blengt/blengtTot)
            IF (lChem) betac(corner(ic))=betac(corner(ic))+w_dis(irecn,ic,isub,ipl)/sum_dis(irecn,isub,ipl)*Weight*(srface/PrvSur)
         END DO
      ENDIF

      !> save Ap position (inside of the soil domain -> xmax and ymax of the soil boudaries not included)
      IF (continu.AND.(isub.GT.1)) THEN
         IF (xAp.EQ.MINVAL(xGrid)+nex*dxgrid) THEN
            xAp=xAp-nex*dxgrid
            transroot(irecn,1,isub,ipl)=transroot(irecn,1,isub,ipl)-1
         ENDIF
         IF (yAp.EQ.MINVAL(yGrid)+ney*dygrid) THEN
            yAp=yAp-ney*dygrid
            transroot(irecn,2,isub,ipl)=transroot(irecn,2,isub,ipl)-1
         ENDIF
         Intc(irecn,1,isub-1,ipl)=xAp            
         Intc(irecn,2,isub-1,ipl)=yAp
         Intc(irecn,3,isub-1,ipl)=zAp
      ENDIF

      IF (SPLIT) THEN
         !> preparation for next loop
         xAp=xInt
         yAp=yInt
         zAp=zInt
         iFaceOld=iFace
         SELECT CASE (iFace)
         CASE(1) 
            zAp=zInt-5.E-5_sp*dzGrid
            IF(continu) transroot(irecn,:,isub+1,ipl)=reftransroot		
         CASE(2) 
            zAp=zInt+5.E-5_sp*dzGrid
            IF(continu) transroot(irecn,:,isub+1,ipl)=reftransroot	
         CASE(3) 
            yAp=yInt-5.E-5_sp*dyGrid
            IF (continu) THEN							
               IF (yAp.LT.MINVAL(ygrid)) THEN
                  yAp=yAp+ney*dygrid
                  transroot(irecn,1,isub+1,ipl)=reftransroot(1)
                  transroot(irecn,2,isub+1,ipl)=reftransroot(2)+1
                  reftransroot=transroot(irecn,:,isub+1,ipl)
               ELSE
                  transroot(irecn,:,isub+1,ipl)=reftransroot

               ENDIF
            ENDIF
         CASE(4) 
            yAp=yInt+5.E-5_sp*dyGrid
            IF (continu) THEN							
               IF (yAp.GT.MINVAL(ygrid)+ney*dygrid) THEN
                  yAp=yAp-ney*dygrid
                  transroot(irecn,1,isub+1,ipl)=reftransroot(1)
                  transroot(irecn,2,isub+1,ipl)=reftransroot(2)-1
                  reftransroot=transroot(irecn,:,isub+1,ipl)
               ELSE
                  transroot(irecn,:,isub+1,ipl)=reftransroot

               ENDIF
            ENDIF
         CASE(5) 
            xAp=xInt-5.E-5_sp*dxGrid
            IF (continu) THEN							
               IF (xAp.LT.MINVAL(xgrid)) THEN
                  xAp=xAp+nex*dxgrid
                  transroot(irecn,1,isub+1,ipl)=reftransroot(1)+1
                  transroot(irecn,2,isub+1,ipl)=reftransroot(2)
                  reftransroot=transroot(irecn,:,isub+1,ipl)
               ELSE
                  transroot(irecn,:,isub+1,ipl)=reftransroot

               ENDIF
            ENDIF
         CASE(6) 
            xAp=xInt+5.E-5_sp*dxGrid
            IF (continu) THEN							
               IF (xAp.GT.MINVAL(xgrid)+nex*dxgrid) THEN
                  xAp=xAp-nex*dxgrid
                  transroot(irecn,1,isub+1,ipl)=reftransroot(1)-1
                  transroot(irecn,2,isub+1,ipl)=reftransroot(2)
                  reftransroot=transroot(irecn,:,isub+1,ipl)
               ELSE
                  transroot(irecn,:,isub+1,ipl)=reftransroot
               ENDIF
            ENDIF
         END SELECT
      ELSEIF (continu.AND.(isub.GT.1)) THEN						
         !> sort transroot in the same order as intc
         transroot(irecn,:,isub+1,ipl)=transroot(irecn,:,1,ipl)
         transroot(irecn,:,1:isub+1,ipl)=transroot(irecn,:,2:isub+2,ipl)

      ENDIF
   END DO splitloop

   !> Ap was on a face of the cube containig A. If Ap was on a boundary too, it is displaced out of the domain (because of the +/-5.E-5...).
   !> Then -> translation. Ap is now closer to B (on the same side and inside of the domain).

   !> macropore
   do isub=1,nsub(irecn,ipl)
      corner=loc_q(irecn,:,isub,ipl)
      dummat=sum(MatNum(corner))/8._sp
      if(dummat.gt.1. .and. dummat.lt.2.) then
         do ic=1,8
            if(MatNum(corner(ic)).eq.2) w_dis(irecn,ic,isub,ipl)=0._dp
         end do
      end if
   end do

   !> save branch basis position 							
   IF (irecn.GT.0) THEN
      IF (continu.AND.(ibrseg(irecpr(irecn)).NE.ibrseg(irecn))) THEN
         Intc(irecn,1,isub,ipl)=xB
         Intc(irecn,2,isub,ipl)=yB
         Intc(irecn,3,isub,ipl)=zB
      ENDIF
   ENDIF
   IF (lDou) THEN
      !> correction of inacurrate seglen
      IF (nsub(irecn,ipl).EQ.1) THEN !nosplit
         w_sub(irecn,1,ipl)=1._dp
      ELSE
         DO isub=1,nsub(irecn,ipl)
            w_sub(irecn,isub,ipl)=l_sub(irecn,isub,ipl)/l_seg(irecn)
         END DO
      ENDIF
   ENDIF
   RETURN
 END SUBROUTINE segment
 !****************************************************************
 !> ### estimate current BC for root system when using Doussan Sink term ###
 SUBROUTINE SetBCroot(t,BCr,BCtp)
   USE TypeDef
   USE PlntData, ONLY : tBCr,typeBCr,BCroot,nBCr 
   IMPLICIT NONE
   REAL(sp), INTENT(in) :: t
   REAL(dp), INTENT(out) ::BCr
   INTEGER (sp),INTENT(out) :: BCtp
   INTEGER (sp) :: ifc

   !> calculate current root UBC-value from input BC(time)-function:
   ifc=0
 201 ifc=ifc+1
   IF (ifc.GT.nBCr) THEN !> beyond the last user-defined BC value it remains constant
      BCr=BCroot(nBCr)
      BCtp=typeBCr(nBCr)
   ELSE
      IF (t.GE.tBCr(ifc)) GOTO 201 !> try to find the first tBCR wwhich is beyond current t (>=1)
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
   IF (BCtp.EQ.2) THEN
      BCr=-ABS(BCr) !> always negative flow at the collar
   ENDIF
   RETURN
 END SUBROUTINE SetBCroot
 !**************************************************************************
!> ### Calculates stomatal closure ###
 SUBROUTINE Rootstress(PHtop,BCtop,BCtptop,iterBC,BC_switch,Jcol,ipl)
   USE Typedef
   USE DoussanMat, ONLY: curr_BCtp,stressBC,hx_min,stresfun,stresval1,stresval2
   USE GridData, ONLY: RelEps,epslonR,factorRelEps
   USE PlntData, ONLY: a_r,a_1,a_2
   USE TempData, ONLY: Tatm_usr,Patm_usr,Pdiff_usr
   USE RootData, ONLY: concol,csign_notrans,PH_crit,delta_h

   IMPLICIT NONE
   REAL(dp), INTENT(inout):: BCtop
   INTEGER, INTENT(inout):: BCtptop
   INTEGER, INTENT(in):: ipl
   REAL(dp), INTENT(in):: PHtop,Jcol
   LOGICAL, INTENT(out) :: BC_switch,iterBC
   REAL(dp):: reduc,del_PHr

   !> check if the root collar abs(PH) is larger than abs(hx_min)+tolerance and adapt the collar BC
   IF (curr_BCtp(ipl)==2) THEN
      IF (RelEps) THEN 
         del_PHr=-MAX(ABS(hx_min/factorRelEps),epslonR)
      ELSE
         del_PHr=-epslonR
      ENDIF
      IF ((stresfun.EQ.1).AND.(PHtop<hx_min+del_PHr)) THEN
         !> + stresfun=1: top node at lower PH than allowed: start of stressed conditions
         PRINT *,'stress in the collar xylem: change to PH BC, PHtop=',PHtop,' is lower than criterion ',hx_min,' +',del_PHr
         stressBC=.TRUE.
         BCtptop=1
         BCtop=hx_min
         BC_switch=.TRUE.
         iterBC=.TRUE.
      ELSEIF ((stresfun.EQ.2).AND.(PHtop.LE.stresval1)) THEN !> + stresfun=2: linear decrease
         reduc=(PHtop-stresval2)/(stresval1-stresval2)
         IF (reduc.GT.1) THEN
            reduc=1
         ELSEIF (reduc.LT.0) THEN
            reduc=0
         ENDIF
         BCtop=-ABS(Jcol)*reduc
         BCtptop=2
         stressBC=.TRUE.
         iterBC=.TRUE.
         PRINT *,'stresfun=2; flux reduction of factor',reduc
      ELSEIF (stresfun.EQ.3) THEN !> + stresfun=3: Tuzet function
         reduc=(1+EXP(stresval1*stresval2))/(1+EXP(stresval1*(stresval2-PHtop)))
         BCtop=-ABS(Jcol)*reduc
         BCtptop=2
         PRINT *,'stresfun=3; flux reduction of factor',reduc
         IF (reduc.LT.1) THEN
            BC_switch=.TRUE.
            iterBC=.TRUE.
         ELSE
            BC_switch=.FALSE.
            iterBC=.FALSE.
         END IF
      ELSEIF (stresfun.EQ.4) THEN !> + stresfun=4: additional signaling
            reduc= a_r + (1-a_r)*exp(-concol*a_1*exp(a_2*(abs(PHtop)-abs(PH_crit))))
            IF (reduc.LT.1) THEN
               PRINT *,'stresfun=4; flux reduction of factor',reduc
               BCtop=-ABS(Jcol)*reduc
               BCtptop=2
               BC_switch=.TRUE.
               iterBC=.TRUE.
            ELSE
               BC_switch=.FALSE.
               iterBC=.FALSE.
            ENDIF
      ELSEIF (stresfun.EQ.5) THEN !> + stresfun=5: additional signaling without particle tracker
         reduc= a_r + (1-a_r)*exp(-csign_notrans*a_1*exp(a_2*(abs(PHtop)-abs(PH_crit))))
         IF (reduc.LT.1) THEN
            BCtop=-ABS(Jcol)*reduc
            BCtptop=2
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

   ELSEIF ((curr_BCtp(ipl)==1).AND.(ABS(Jcol)>1.00001*ABS(BCtop))) THEN !factor to avoid oscillation between stress and no stress
      !> check if end of stress (curr_BCtp is PH with PH is hxmin but BCtp_new is still flux)
      PRINT *,'end of stress conditions, collar flux=',ABS(Jcol),' is larger than the prescribed flux', ABS(BCtop)
      !	  BCtptop=BCnew is already of type 2!
      stressBC=.FALSE.
      BC_switch=.FALSE.
      PRINT *,'info',stressBC,BCtptop,BCtop, BC_switch,iterBC

   ELSE
      BC_switch=.FALSE.
      iterBC=.FALSE.
      IF (stressBC) PRINT *,'no change: stress=',stressBC
   ENDIF

 END SUBROUTINE Rootstress
 !*************************************************************
!> ### allocate radial and axial root conductivities to root segments acc. to their age ###
 SUBROUTINE ConductRoot(t)
   USE TypeDef
   USE DoussanMat, ONLY : nLr,Lr,Lrroot,ageLr,nKh,Khr,Khroot,ageKh,Khr_pot,Lr_pot,cavitfun,transroot
   USE RootData, ONLY: nrec,timorg,age,ordseg,maizeroottyp,loliumroottyp,wheatroottyp,lGap,lAQPc,xs,ys,zs
   USE SolData, ONLY: MatNum,nmat,par,l_elmMacro
   USE GridData, ONLY: n_neigh,nex,ney,dxgrid,dygrid,continu
   IMPLICIT NONE

   REAL(sp), INTENT (in) :: t
   REAL(sp) :: segage,xC,yC
   INTEGER(ap) :: irecn,iage,typ,corner(8),imin,iMat,nM,c

   !> calculate segment lateral and long. conductivity according to age:
   DO irecn=1,nrec
      segage=t-timorg(irecn)
      typ=ordseg(irecn) 
      IF (maizeroottyp) THEN
         IF (typ.LT.12) THEN  		!In RootTyp, maize root types below 12 are principal roots (12 and 13 are lateral) 
            typ=1
         ELSE
            typ=2
         ENDIF
      ELSEIF (loliumroottyp) THEN
         IF (typ.LT.3) THEN  		!In RootTyp, lolium root types below 3 are principal roots
            typ=1
         ELSEIF (typ.LT.5) THEN
            typ=2
         ELSE
            typ=3
         ENDIF
      ELSEIF (wheatroottyp) THEN
         IF (typ.LT.19) THEN  		!In RootTyp, wheat root types below 19 are principal roots 
            typ=1
         ELSEIF (typ.LT.20) THEN
            typ=2
         ELSE
            typ=3
         ENDIF
      ELSE
         IF (typ.EQ.0) THEN  !type 0 is the seed+small segment
            typ=1
         ELSEIF (typ.GT.3) THEN!no more than 3 root types
            typ=3
         ENDIF
      ENDIF
      !> radial conductivity matrix
      iage=1
      DO WHILE (ageLr(typ,iage).le.segage)
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
      

      !> if root segments are within a non-conductive soil material --> set radial conductivity to zero, or scale Lr according to the amount of non-cond. material
      !> macropore
      IF(nMat.GT.1) THEN
         DO iMat=1,nMat 
            IF (par(6,iMat) .LE. 1.0E-05_sp) THEN !non-conductive
               IF(continu) THEN
                  xC=xs(irecn)+transroot(irecn,1,1,1)*nex*dxgrid
                  yC=ys(irecn)+transroot(irecn,2,1,1)*ney*dygrid
                  CALL Neighb(xC,yC,zs(irecn),corner,imin)
               ELSE
                  CALL Neighb(xs(irecn),ys(irecn),zs(irecn),corner,imin)
               END IF
               IF (ANY(MatNum(corner) .EQ. iMat))  Lr(irecn) = 0._dp 
            ELSEIF(par(11,iMat).LE.1E-5_sp)  THEN  !macropore
               IF(continu) THEN
                  xC=xs(irecn)+transroot(irecn,1,1,1)*nex*dxgrid
                  yC=ys(irecn)+transroot(irecn,2,1,1)*ney*dygrid
                  CALL Neighb(xC,yC,zs(irecn),corner,imin)
               ELSE
                  CALL Neighb(xs(irecn),ys(irecn),zs(irecn),corner,imin)
               END IF
               nM=8  !> number of mat.nodes within bulk material
               DO c=1,8
                  IF(MatNum(corner(c)) .EQ. iMat) nM=nM-1
               END DO
               Lr(irecn) = Lr(irecn)*nM/8._dp 
               IF (l_elmMacro(imin) .AND. n_neigh(imin).LT.1)  Lr(irecn) = 0.000001_dp
            END IF
         END DO
      END IF
        
      !> axial conductance matrix
      iage=1
      DO WHILE (ageKh(typ,iage).le.segage)
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
      IF (cavitfun.EQ.1) Khr_pot(irecn)=Khr(irecn)
      IF (lGap.OR.lAQPc) Lr_pot(irecn)=Lr(irecn)
   END DO

 END SUBROUTINE ConductRoot
!******************************************************************************
!> ### calculate sink for lDou simulations
!> solve DOUSSAN flow equations within the xylem ###
 SUBROUTINE SolveRoot(t,dt,it1,iter)
   USE TypeDef
   USE SparseMatrix
   USE NumericalRecipes
   USE Paramdata,ONLY: MaxNod, pi, ldirect,lChem
   USE GridData, ONLY: xgrid,ygrid,zgrid, sink, betaw,dxgrid,dygrid,dzgrid,Wn,elmnod,nPt,nElm,RootSk,Axy,Bxy,Dxy,Exy,nex,ney,continu, sink_cube, betac_cube, betac_cube2,n_neigh,MacroList
   USE PlntData, ONLY :Tpot,Tact,LA
   USE DoussanMat, ONLY: stressBC,Q_bc1,Q_bc2,Q_bc,Qi,PHs,Qd,PHr,GH,w_dis,nsub,sinkR,Jintot,nrecOld,loc_Q,nBCn,bcr_usr,counter2, &
        bctp_usr,curr_bcr,curr_bctp,Khr,w_sub,NBC_irecn,sum_dis,cent,Intc,Lr,ave,old,axialRootFlow,Joutr,PH_root_sub, &
        cp_mean,oldT,hx_min,plantmatrix,nplant, cube_i, beta_weight,jao,iao,ao,iro,jco,aij,phrold,cavitfun,cavitb,cavitc,Khr_pot,Lr_pot,solveroot_call,veloRoot
   USE RootData, ONLY: nrec,seglen,segsur,nbr,ibrseg,irecpr,lno_root_growth,lCalloc,crossSectionSeg,lGap,lAQPc
   USE SolData, only: hnew,MatNum,lMacro,l_elmMacro
   USE tmctrl, ONLY: t_begin, tlevel,dtroot
   IMPLICIT NONE

   REAL(sp) ,INTENT(in) :: t,dt
   LOGICAL,INTENT(in) :: it1
   REAL (dp)::PHr_tot(1:nrec+1),PH_direct(1:nrec+1), err_loc,Qd2(1:2*(nrec+1)),dVolSk,Voltot,checksub
   REAL (dp)::maxerr=1.e-10_dp,BCr_new,BCr_old,Jcol,rootsk2,sink_dummy,betac_dummy
   REAL (dp), SAVE::temp_phr
   REAL (sp) :: totalAxialFlow=0,Gap,AQPc
   REAL (sp) :: rs,concrs
   INTEGER(ap) ::old_BCtp,iter_loc,itol,itmaxl,iter,corner_ic,ibr,c_i,i,i_n,i_dummy
   INTEGER(ap) :: iprvn,irecn,ifoln,iBC,ic,BCtp_new,isub,err,ipl,num_elem
   LOGICAL :: repet,iterBC,BC_switch,n_apex,run
   LOGICAL ,SAVE::firstcall=.TRUE.

   IF(t.LE.dtRoot) phr_tot(1:nrec+1)=0._dp
   
   IF (ave)  counter2=0
   
   solveroot_call=solveroot_call+1     
   
   !> initialisation
   sink=0._dp
   RootSk=0.0_dp
   RootSk2=0.0_dp
   sink_cube=0.0_dp
   betac_cube=0.0_dp

   DO ipl=1,nplant    
      
      IF ((cavitfun.EQ.1.OR.lGap.OR.lAQPc).AND.solveroot_call.EQ.1) THEN							
         !> part added from setupdou allowing a change in the root axial conductance due to cavitiation during the simulation 
         err=SM_allocate(plantmatrix(ipl), nrec+1, nrec+1)
         IF(err/=0) Stop 'Could not create plantmatrix'
         DO ibr=1,nbr !> all plants have same number of root segments!
            n_apex=.false.
            !> find the tip segment of the branch 'ibr'
            irecn=nrec
            DO WHILE (ibrseg(irecn).NE.ibr)
               irecn=irecn-1
            END DO
            IF (seglen(irecn)<1.E-20) THEN ! skip this segment too small to be taken
               ifoln=irecn ! following node ID
               irecn=irecpr(irecn) !current node ID
            ELSE
               n_apex=.TRUE.!ifoln does not exist if not continu
            ENDIF
            IF (irecn==0) THEN !there exists a branch ibr but not yet any segment!
               run=.false.
            ELSE
               run=.true.
            ENDIF
            !> then the rest of the branch up to the seed or the embranchment
            DO WHILE (run)
               !> "upper" node
               iprvn=irecpr(irecn)
               !> several changes -> multiple roots
               err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,-Khr(irecn)/seglen(irecn),.false.)
               If(err/=0) Stop 'Could not insert element into plantmatrix'
               !> if apex (bottom part of root)
               IF (n_apex) THEN
                  Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn)/seglen(irecn)+Lr(irecn)*segsur(irecn),.false.)
                  If(err/=0) Stop 'Could not insert element into plantmatrix'
               ELSE
                  Err=SM_set(Plantmatrix(ipl), irecn+1,ifoln+1,-Khr(ifoln)/seglen(ifoln),.false.) !row, col,value
                  Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn)/seglen(irecn)+Khr(ifoln)/seglen(ifoln)+Lr(irecn)*segsur(irecn),.false.)
                  If(err/=0) Stop 'Could not insert element into plantmatrix'
               ENDIF
               !> define 1st part of Q (Q=Qi.*PHsoil+Qbc) -> RHS
               Qi(irecn,ipl)=Lr(irecn)*segsur(irecn)
               IF (iprvn==0) THEN!> seed=first segment
                  IF (curr_BCtp(ipl)==2) THEN !> flux BC
                     Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn)/seglen(irecn),.false.)
                     Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn)/seglen(irecn),.false.)
                     If(err/=0) Stop 'Could not insert element into plantmatrix'
                  ELSE IF (curr_BCtp(ipl)==1) THEN !> PH BC
                     Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,1._dp,.true.)!position (1,1)
                     ! true means overwrite (delete first), no addition
                     Err=SM_set(Plantmatrix(ipl), irecn+1,iprvn+1,0._dp,.true.)!irecn+1,iprvn+1
                     If(err/=0) Stop 'Could not insert element into plantmatrix'
                  ENDIF
                  run=.false.
               ELSEIF (ibrseg(iprvn).NE.ibrseg(irecn)) THEN !start of the branch but not from the seed -> gaat van onder na boven, duz iprvn (nieuwe positie) is nu niet van die ene branch maar van een side branch
                  Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn)/seglen(irecn),.false.)
                  Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn)/seglen(irecn),.false.)
                  If(err/=0) Stop 'Could not insert element into plantmatrix'
                  run=.false.
               ENDIF
               !> definition for the next run of the loop
               ifoln=irecn
               irecn=iprvn
               !> from here, not an apex
               n_apex=.false.
            END DO ! loop on branch nodes
         END DO ! loop on root branches
      ENDIF ! if cavitfun or gap or AQPc
      
      
      !> check current plant BC
      BCr_old=BCr_usr(ipl)
      IF (.NOT.lCalloc) THEN
         CALL SetBCroot(t,BCr_usr(ipl),BCtp_usr(ipl))
      ELSE
         IF (ipl.GE.2) STOP 'Assimilate allocation currently doesnt work with multiple plants'
         CALL Stress(rs,concrs)
         CALL SetTp(t,rs,concrs,ipl)
         BCr_usr(ipl)=ABS(Tpot)
         BCtp_usr(ipl)=2
      ENDIF
      !> calculation of Tpot if BCtp==2
      IF (BCtp_usr(ipl)==2) Tpot=+ABS(BCr_usr(ipl))
      !> if stress at previous time step, then BCtp new=2 and currBctp=1
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
         !> adapt BC
         IF (firstcall.OR.((BCr_new/=curr_BCr(ipl)).AND.(.NOT.(stressBC))).OR.((BCtp_new/=curr_BCtp(ipl)).AND.(.NOT.(stressBC))).OR.(BC_switch)) THEN
            !> updating BC
            curr_BCr(ipl)=BCr_new
            old_BCtp=curr_BCtp(ipl)
            curr_BCtp(ipl)=BCtp_new
            
            !> updating the BC part of the matrices Qbc and Cd
            iprvn=0 !> seed=first segment
            IF (curr_BCtp(ipl)==2.AND.ldirect) THEN  !> 1. Flux and direct solver
               IF (.NOT.tlevel) THEN

                  DO iBC=1,nBCn  !> all the branches connected to the seed
                     irecn=nBC_irecn(iBC) !connected to the seed
                     err=SM_set(plantmatrix(ipl), iprvn+1,iprvn+1,1._dp,.TRUE.)!position (1,1)
                     err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,0._dp,.TRUE.)!position (2,1)
                     err=SM_set(plantmatrix(ipl),iprvn+1,irecn+1,0._dp,.TRUE.)!position (1,2)
                     IF(err/=0) STOP 'Could not insert element into plantmatrix'
                  ENDDO
               ENDIF

               DO iBC=1,nBCn  !all the branches connected to the seed
                  irecn=nBC_irecn(iBC)
                  temp_phr=-(curr_BCr(ipl)/Khr(irecn)*seglen(irecn+1)+200.)

                  Qi(iprvn,ipl)=0._dp !> from seed has to be given, from irecn is already stated in DoussanMat
                  Q_bc1(iprvn,ipl)=temp_Phr*.1+GH(iprvn) !> use ph from last time step
                  Q_bc2(iprvn,ipl)=temp_Phr+GH(iprvn)
                  Q_bc1(irecn,ipl)=(temp_Phr*.1+GH(iprvn))*Khr(irecn)/seglen(irecn) !> all the PH for root was given in total head!!
                  Q_bc2(irecn,ipl)=(temp_Phr+GH(iprvn))*Khr(irecn)/seglen(irecn)
               ENDDO

            ELSE IF (curr_BCtp(ipl)==2.AND..NOT.ldirect) THEN  !> 2. Flux and iter solver
               IF (old_BCtp/=curr_BCtp(ipl)) THEN
                  DO iBC=1,nBCn  !all the branches connected to the seed
                     irecn=nBC_irecn(iBC)!connected to the seed
                     IF (iBC==1) THEN
                        repet=.TRUE.!delete previous value
                     ELSE
                        repet=.FALSE.!add to previous value
                     ENDIF
                     err=SM_set(plantmatrix(ipl),iprvn+1,iprvn+1,Khr(irecn)/seglen(irecn),repet)  !position (1,1)
                     err=SM_set(plantmatrix(ipl),iprvn+1,irecn+1,-Khr(irecn)/seglen(irecn),repet) !position (1,2)
                     err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,-Khr(irecn)/seglen(irecn),repet) !position (2,1)
                     IF(err/=0) STOP 'Could not insert element into plantmatrix'
                     Q_bc1(irecn,ipl)=0._dp
                  ENDDO
                  !> from seed has to be given, from irecn is already stated in DoussanMat
                  Qi(iprvn,ipl)=0._dp
               ENDIF
               Q_bc1(iprvn,ipl)=curr_BCr(ipl) !Q(0,0) a flow imposed

            ELSE IF (curr_BCtp(ipl)==1) THEN  !PH
               IF (old_BCtp/=curr_BCtp(ipl)) THEN
                  err=SM_set(plantmatrix(ipl), iprvn+1,iprvn+1,1._dp,.TRUE.)
                  IF(err/=0) STOP 'Could not insert element into plantmatrix'
               ENDIF
               Qi(iprvn,ipl)=0._dp
               !> all the PH for root was given in total head!!
               Q_bc1(iprvn,ipl)=(curr_BCr(ipl)+GH(iprvn)) 
               DO iBC=1,nBCn
                  irecn=nBC_irecn(iBC) !connected to the seed
                  IF (old_BCtp/=curr_BCtp(ipl)) THEN
                     err=SM_set(plantmatrix(ipl),iprvn+1,iprvn+1,1._dp,.TRUE.)
                     err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,0._dp,.TRUE.)
                     err=SM_set(plantmatrix(ipl),iprvn+1,irecn+1,0._dp,.TRUE.)
                     IF(err/=0) STOP 'Could not insert element into plantmatrix'
                  ENDIF
                  Q_bc1(irecn,ipl)=(curr_BCr(ipl)+GH(iprvn))*Khr(irecn)/seglen(irecn) !all the PH for root was given in total head!!
                  Q_bc2 = 0._dp
               END DO
            ENDIF
         ENDIF

         IF (firstcall) firstcall=.FALSE.

         IF (oldT) GOTO 201
         IF (t.EQ.t_begin+dt .AND. iter .GT. 1 .OR. t.GT.t_begin+dt) THEN
            CALL analytical_approach(ipl)
         ELSE
            old=.TRUE.
         ENDIF
         !> for t=tstep. Initially no root flow is calculated, and not all the boundary conditions for the analytical problem are available, so solve
         !> first time step using the averaging way
         !> find current PH soil and create Q
 201     IF (old) THEN
            CALL average_method(ipl)
         ENDIF
         IF (nrec.EQ.1) THEN!only seed
            IF (curr_BCtp(ipl)==2) THEN
               Q_bc1=Q_bc1 ! this must be checked MJ dec2008
            ENDIF
         ELSE
            !> put in format (1:nRec+1) for solving the system
            Qd2(1:nrec+1)=Qd(0:nrec,ipl)
            Qd2(nrec+2:) = Qd(nrec+1:,ipl)
            PHr_tot(1:nrec+1)=PHr(1:nrec+1,ipl)+GH(0:nrec) !total water potential
            !print*,'PHr',Phr(1:100,1)
            !print*,'qd2',Qd2(1:100)
            IF (it1) THEN
               PHr_tot(1:nrec+1)=PHs(0:nrec,ipl)   !> initial xylem PH -> could also be zero but then convergence should take longer (this should be the common case)
            ENDIF
            IF (curr_BCtp(ipl)==1) THEN
               PHr_tot(1)=curr_BCr(ipl)+GH(0) !> bound value should always be incorporated if a PH is imposed
            ENDIF
            itol=3
            itmaxl=(nrec+1)*3 !> maximum iterations for biconjugate gradient-> needed in case of root growth
            !> solve with biconjugated gradient method (Numerical recipies)

            IF (ldirect) THEN
               num_elem= SM_numberOfElements(plantmatrix(ipl))
               CALL SolveRootDirect(ipl,num_elem,Qd2,Ph_direct,temp_phr)
               phr_tot = ph_direct(1:nrec+1)
            ELSE
               CALL linbcg(nrec+1,Qd2(1:nrec+1),PHr_tot,itol,maxerr,itmaxl,iter_loc,err_loc,ipl)

!!$               print *, 'flux - linbcg: ', Khr(1)/seglen(1)*(Phr_tot(1)-phr_tot(2))
!!$               print *, 'flux - direct: ', Khr(1)/seglen(1)*(Ph_direct(1)-ph_direct(2))
!!$               print *, 'linbcg: ',phr_tot(1:3)
!!$               print *, 'mumps: ',ph_direct(1:3)
!!$               print *, '**************************'
!!$               print *, 'Norm: ',sqrt(sum(ph_direct**2))- sqrt(sum(phr_tot**2))
!!$               WRITE( *, * ) 'Diff is:  ',(abs(phr_tot(I))-abs(ph_direct(I)),I=1,size(phr_tot))
!!$               print *, 'Max: ', maxval(abs(abs(phr_tot)-abs(ph_direct)))
!!$               print *, 'Min: ', minval(abs(abs(phr_tot)-abs(ph_direct)))
         
               IF (err_loc.LE.maxerr) THEN
                  nrecOld=nrec
               ELSE
                  PRINT *,'No convergence in the root system'
                  STOP
               ENDIF
            ENDIF !endif direct solver or not

         ENDIF

             !> ====================== end solve root system ====================================
             !> ----------------------calculate axial root flow-------------------------------
 !> for each branch cumulative flow is calculated
         Branchloop: DO ibr=1,nbr
            n_apex=.FALSE.
            !> find the tip segment of the branch 'ibr'
            irecn=nrec
            DO WHILE (ibrseg(irecn).NE.ibr)
               irecn=irecn-1
            END DO
            !> the first one we find is an apex
            IF (seglen(irecn)<1.E-20) THEN !> skip this segment too small to be taken
               ifoln=irecn !> following node ID
               irecn=irecpr(irecn) !> current node ID
            ELSE
               n_apex=.TRUE.!> ifoln does not exist
            ENDIF
            IF (irecn==0) THEN !> there exists a branch ibr but not yet any segment!
               run=.FALSE.
            ELSE
               run=.TRUE.
            ENDIF
            !> then the rest of the branch up to the seed or the embranchment
            DO WHILE (run)
               !* "upper" node
               iprvn=irecpr(irecn)
               axialRootFlow(irecn,ipl) = Khr(irecn)*(PHr_tot(irecn+1)-PHr_tot(iprvn+1))/seglen(irecn)
               !> if reached the seed or the embranchement => change branch
               IF (iprvn==0.OR.(ibrseg(iprvn).NE.ibrseg(irecn))) run=.FALSE.			
               IF (irecn .EQ. 0) THEN
                  IF (curr_BCtp(ipl) == 2) THEN !> flux
                     axialRootFlow(irecn,ipl) = Q_bc1(irecn,ipl)
                  ENDIF
                  run=.FALSE.
               ENDIF
               !> get total flow; which is addition of all flows at each large branch connecting to the seed		Quid des branches latérales ? 
               IF (iprvn .EQ. 0) THEN
                  totalAxialFlow = totalAxialFlow + axialRootFlow(irecn,ipl)
               ENDIF

               !> get velocity for each segement -> needed in soluteRoot
               veloRoot(irecn,ipl) = axialRootFlow(irecn,ipl)/crossSectionSeg(irecn)

               !> definition for the next run of the loop
               ifoln=irecn
               irecn=iprvn
               !> from here, not an apex
               n_apex=.FALSE.
            END DO ! branch nodes loop
         END DO Branchloop 
         !write(*,*)'total axialRootflow',totalAxialFlow
         totalAxialFlow = 0
         !> ---------------------------------end axial root flow calculation-----------------
         !***********************************************************************
         !> calculate SINK term
         !**********************************************************************

         !> initialisation
         betaw=0.
         sink_cube=0.
         Jintot=0.
         Jcol=0.0_dp
         Voltot=0.0_dp
         DO irecn=0,nrec
            Joutr(irecn) = Qi(irecn,ipl)*(PHs(irecn,ipl)-PHr_tot(irecn+1)) !> \param Joutr [L3/T] (soil is still from 0 to nrec while root is from 1 to nrec+1
            sinkR(irecn,ipl)=Joutr(irecn) !> \param sinkR root flow! [L3/T]
            Jintot=Jintot+Joutr(irecn)/(dxgrid*dygrid*dzgrid*nsub(irecn,ipl))
            IF (t.EQ.t_begin+dt .AND. iter .GT. 1 .OR. t.GT.t_begin+dt) THEN
               checksub=0
               DO isub=1,nsub(irecn,ipl)
                  dVolSk=0.
                  c_i=cube_i(irecn,isub,ipl)
                  IF(irecn.EQ.0) c_i=cube_i(1,isub,ipl)

                  IF(lMacro .AND. l_elmMacro(c_i)) THEN
                     !> Uptake from voxels with mixed material is shifted to bulk soil
                     IF(n_neigh(c_i).gt.0) then
                        sink_dummy = Joutr(irecn)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/real(n_neigh(c_i))
                        betac_dummy = beta_weight(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/real(n_neigh(c_i))
                        DO i_n=1,6
                           IF(MacroList(c_i,i_n).GT.0) THEN
                              i_dummy = MacroList(c_i,i_n)
                              Sink_cube(i_dummy) = Sink_cube(i_dummy) + sink_dummy
                              betac_cube(i_dummy) = betac_cube(i_dummy) +  betac_dummy
                           END IF
                        END DO
                     ELSE
                        Sink_cube(c_i) = 0._dp
                        betac_cube(c_i) = 0._dp
                     END IF
                  ELSE
                     Sink_cube(c_i)=Sink_cube(c_i) + Joutr(irecn)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
                     betac_cube(c_i)=betac_cube(c_i) + beta_weight(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
                  END IF

                  DO ic=1,8
                     corner_ic=loc_Q(irecn,ic,isub,ipl)
                     IF(irecn.EQ.0) corner_ic=loc_Q(1,ic,isub,ipl)
                     !> delta flux extracted
                     dVolSk=Joutr(irecn)*w_dis(irecn,ic,isub,ipl)/sum_dis(irecn,isub,ipl)*w_sub(irecn,isub,ipl)
                     !sink term
                     sink(corner_ic)=sink(corner_ic)+dVolSk/(dxgrid*dygrid*dzgrid)/Wn(corner_ic)
                     betaw(corner_ic)=1._dp
                     !> sum of extracted fluxes

                     !RootSk=RootSk+dVolSk

                     Voltot=Voltot+w_dis(irecn,ic,isub,ipl)/sum_dis(irecn,isub,ipl)*w_sub(irecn,isub,ipl)
                  END DO
                  checksub=checksub+w_sub(irecn,isub,ipl)
               END DO !>end sub-segments

               !first time step
            ELSE 
               DO isub=1,nsub(irecn,ipl)
                  c_i=cube_i(irecn,isub,ipl)
                  IF(irecn.EQ.0) c_i=cube_i(1,isub,ipl)
                  IF(lMacro .AND. l_elmMacro(c_i)) THEN

                     !Uptake from voxels with mixed material is shifted to bulk soil
                     IF(n_neigh(c_i).gt.0) then
                        sink_dummy = Joutr(irecn)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/real(n_neigh(c_i))
                        betac_dummy = beta_weight(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/real(n_neigh(c_i))
                        DO i_n=1,6
                           IF(MacroList(c_i,i_n).GT.0) THEN
                              i_dummy = MacroList(c_i,i_n)
                              Sink_cube(i_dummy) = Sink_cube(i_dummy) + sink_dummy
                              betac_cube(i_dummy) = betac_cube(i_dummy) +  betac_dummy
                           END IF
                        END DO
                     ELSE
                        Sink_cube(c_i) = 0._dp
                        betac_cube(c_i) = 0._dp
                     END IF

                  ELSE
                     Sink_cube(c_i)=Sink_cube(c_i) + Joutr(irecn)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
                     betac_cube(c_i)=betac_cube(c_i) + beta_weight(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
                     !print*,'sink_cube',Sink_cube(c_i),'c_i',c_i,'irec',irecn
                  END IF
                  
                  DO ic=1,8
                     corner_ic=loc_Q(irecn,ic,isub,ipl)
                     IF(irecn.EQ.0) corner_ic=loc_Q(1,ic,isub,ipl)
                     sink(corner_ic)=sink(corner_ic)+Joutr(irecn)*w_sub(irecn,isub,ipl)*w_dis(irecn,ic,isub,ipl) &
                          /sum_dis(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/Wn(corner_ic)
                     betaw(corner_ic)=1._dp ! betaw is defined for the Reset subroutine
                  END DO
               END DO
            ENDIF
            Jcol=Jcol+Joutr(irecn)
         END DO ! nodes loop

         RootSk = SUM(sink_cube*dxgrid*dygrid*dzgrid)
         betac_cube2 = betac_cube/sum(betac_cube)

         !print *, 'sink: ', rootsk,jcol

         !> actual transpiration
         IF (Jcol.GE.0) THEN
            Tact(ipl)=Jcol
         ELSE
            Tact(ipl)=0
         END IF
         !> gravitational component subtracted again to represent data in outRoo.x without gravity
         PHr(1:nrec+1,ipl)=PHr_tot(1:nrec+1)-GH(0:nrec)
         PHs(0:nrec,ipl)=PHs(0:nrec,ipl)-GH(0:nrec)      

         IF (cavitfun.EQ.1) THEN								
            ! Khr still need to be expressed with an additional plant# dimension
            Khr(1:nrec)=Khr_pot(1:nrec)*EXP(-(-PHr(2:nrec+1,ipl)/cavitb)**cavitc)
         ENDIF
         IF (lGap.OR.lAQPc) THEN								
            DO irecn=1,nrec
                Lr(irecn)=Lr_pot(irecn)*Gap(PHs(irecn,ipl))*AQPc(PHs(irecn,ipl))	
                !> Gap is 1 in case lGap is .false., same for AQPc. If we wanted to be exact, the total factor should be 
                !> 2*Gap*AQPc/(Gap+AQPc), since the two conductances are in series
            ENDDO
         ENDIF
            
         IF (ave) THEN
            PH_root_sub(0:nrec,:)=PH_root_sub(0:nrec,:)-cp_mean(0:nrec,3,:,ipl)
         ELSEIF (.NOT.(ave)) THEN
            PH_root_sub(0:nrec,:)=PH_root_sub(0:nrec,:)-cent(0:nrec,3,:,ipl)
         ENDIF
         IF (.NOT.lCalloc) THEN
            IF (((curr_BCtp(ipl).EQ.2).OR.(stressBC)).AND.(.NOT.(BC_switch))) THEN
               !> check stress a the collar in case of BC flux
               !> check flux at the collar in case of PH BC+stress
               !> if swicth at previous iterBC, no check
               CALL Rootstress(PHr(1,ipl),BCr_new,BCtp_new,iterBC,BC_switch,Jcol,ipl)
               !> change BCtp_new=currBCtp to 1 if stress occurs
            ELSE
               iterBC=.FALSE.
               BC_switch=.FALSE.
            ENDIF
         ENDIF
      ENDDO ! do while iterBC !actual transpiration
   ENDDO ! plant loop
   RETURN
 END SUBROUTINE SolveRoot
!**************************************************************************************
 FUNCTION Gap(h)
   USE TypeDef
   USE RootData, ONLY: lGap,g1,g2
   IMPLICIT NONE
   REAL(dp), INTENT(in) ::h
   REAL(sp):: Gap
   Gap=1.0_sp
   IF (lGap) THEN
      IF ((h.LT.g1).AND.(h.GT.g2)) THEN
         Gap=(h-g2)/(g1-g2)
      ELSEIF (h.LE.g2) THEN
         Gap=1.E-9_sp
      ENDIF
   ENDIF
   RETURN
 END FUNCTION Gap
 !**************************************************************************************
 FUNCTION AQPc(h)
   USE TypeDef
   USE RootData, ONLY: lAQPc,AQPh,AQPv,nAQPc
   IMPLICIT NONE
   REAL(dp), INTENT(in) ::h
   REAL(sp):: AQPc
   INTEGER(sp):: ih
   AQPc=1.0_sp
   IF (lAQPc) THEN
      ih=1
      DO WHILE (AQPh(ih).GE.h)
         ih=ih+1
      ENDDO
      IF (ih.GT.nAQPc) THEN
         AQPc=AQPv(nAQPc)
      ELSEIF (ih.EQ.1) THEN
         AQPc=AQPv(1)
      ELSE
         AQPc=AQPv(ih-1)+(AQPv(ih)-AQPv(ih-1))*(h-AQPh(ih-1))/(AQPh(ih)-AQPh(ih-1))
      ENDIF
   ENDIF
   RETURN
 END FUNCTION AQPc
!******************************************************************************
 !> ### calculates Couvreur et al. (2012) RWU parameters
 !> solves DOUSSAN flow equations within the xylem ###
 SUBROUTINE SSFdis(kout)
   USE TypeDef
   USE SparseMatrix
   USE NumericalRecipes
   USE GridData, ONLY: dxgrid,dygrid,dzgrid,zGrid,nElm,nPt,nex,ney,nez,SSF,Wn,nx,ny,nz
   USE DoussanMat, ONLY: GH,w_dis,nsub,loc_Q,nBCn,Khr,w_sub,sum_dis,Lr,Lr_pot,ave,old,PH_root_sub,oldT,hx_min,plantmatrix,nplant,cube_i,Inv_ciiKrKr,Inv_c1KxKr,PHs,stresfun
   USE RootData, ONLY: nrec,seglen,segsur,nbr,ibrseg,irecpr,Krs,Kcomp,ldJvL,lSinkCube,lSUF
   USE SolData, ONLY: hNew
   IMPLICIT NONE
   
   REAL (dp)::PHr_tot1(1:nrec+1),PHr_tot2(1:nrec+1),err_loc,Q1(1:nrec+1),Q2(1:nrec+1),Joutr1(1:nrec),Joutr2(1:nrec)
   REAL (dp)::PHs1(0:nrec),PHs2(0:nrec),maxerr=1.e-10_dp,Jcol1,Jcol2,SSFt,SSFcum,sumPhiHsi,sumPhi,sumPhi2,sumHsi,sumHsi2,Beta(1:nx*ny*nz)
   REAL (sp)::Hsi,r,iz,Phi,PhiSUF(1:nrec),Hseq,Lr_fact(1:nrec),Gap,AQPc
   REAL(dp), ALLOCATABLE,DIMENSION(:,:) ::SSFunstockNew,SSFunstock
   REAL(sp), ALLOCATABLE,DIMENSION(:) ::sink_cube1,sink_cube2,sink1,sink2,Hs_cube2
   INTEGER(ap), ALLOCATABLE,DIMENSION(:) :: ind
   INTEGER(ap) ::old_BCtp,iter_loc,itol,itmaxl,iter,corner_ic,ibr,c_i,i,j,k,nk
   INTEGER(ap) :: iprvn,irecn,ifoln,ic,isub,err,ipl
   INTEGER(ap), INTENT(in)::kout
   LOGICAL :: n_apex,run
   CHARACTER :: file*17

   WRITE (*,*) '... Calculating SSF and Krs ...'
   IF (.NOT.lSUF.OR.kout.EQ.0) THEN
      WRITE (*,*) '... Solving Doussan system ...'
      ALLOCATE (sink_cube1(1:nElm))
      ALLOCATE (sink_cube2(1:nElm))
      ALLOCATE (Hs_cube2(1:nElm))
      sink_cube1=0._sp
      sink_cube2=0._sp
      Hs_cube2=0._sp
      ALLOCATE (sink1(1:nPt))
      ALLOCATE (sink2(1:nPt))
      sink1=0._sp
      sink2=0._sp
      !> PHs homogeneous when characterizing SSF and Krs
      PHs1(0:nrec)=-300.0_dp
      !> PHs heterogeneous when characterizing Kcomp
      PHs2(0:nrec)=-150.0_dp+2*GH(0:nrec)
      DO ipl=1,nplant      
         err=SM_allocate(plantmatrix(ipl), nrec+1, nrec+1)
         IF(err/=0) STOP 'Could not create plantmatrix'
         DO ibr=1,nbr !all plants have same number of roots!
            run=.TRUE.
            n_apex=.FALSE.
            !* find the tip segment of the branch 'ibr'
            irecn=nrec
            DO WHILE (ibrseg(irecn).NE.ibr)
               irecn=irecn-1
            END DO
            IF (seglen(irecn)<1.E-20) THEN ! skip this segment too small to be taken
               ifoln=irecn ! following node ID
               irecn=irecpr(irecn) !current node ID
            ELSE
               n_apex=.TRUE.!ifoln does not exist if not continu
            ENDIF

            !* then the rest of the branch up to the seed of the embranchment
            DO WHILE (run)
               !* "upper" node
               iprvn=irecpr(irecn)
               Q1(irecn+1)=Lr(irecn)*segsur(irecn)*PHs1(irecn)
               Q2(irecn+1)=Lr(irecn)*segsur(irecn)*PHs2(irecn)
               err=SM_set(plantmatrix(ipl),irecn+1,iprvn+1,-Khr(irecn)/seglen(irecn),.FALSE.)
               IF(err/=0) STOP 'Could not insert element into plantmatrix'
               !if apex (bottom part of root)
               IF (n_apex) THEN
                  Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn)/seglen(irecn)+Lr(irecn)*segsur(irecn),.FALSE.)		!Here Lr corresponds to Lr_pot if Solveroot was not called yet
                  IF(err/=0) STOP 'Could not insert element into plantmatrix'
               ELSE
                  Err=SM_set(Plantmatrix(ipl), irecn+1,ifoln+1,-Khr(ifoln)/seglen(ifoln),.FALSE.) !row, col,value
                  Err=SM_set(Plantmatrix(ipl), irecn+1,irecn+1,Khr(irecn)/seglen(irecn)+Khr(ifoln)/seglen(ifoln)+Lr(irecn)*segsur(irecn),.FALSE.)
                  IF(err/=0) STOP 'Could not insert element into plantmatrix'
               ENDIF
               IF (iprvn==0) THEN!seed=first segment
                  Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn)/seglen(irecn),.FALSE.)
                  Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn)/seglen(irecn),.FALSE.)
                  IF(err/=0) STOP 'Could not insert element into plantmatrix'
                  Q1(iprvn+1)=-100.0_dp
                  Q2(iprvn+1)=-100.0_dp
                  run=.FALSE.
               ELSEIF (ibrseg(iprvn).NE.ibrseg(irecn)) THEN!start of the branch but not from the seed
                  Err=SM_set(Plantmatrix(ipl), iprvn+1,iprvn+1,Khr(irecn)/seglen(irecn),.FALSE.)!iprvn+1,iprvn+1
                  Err=SM_set(Plantmatrix(ipl), iprvn+1,irecn+1,-Khr(irecn)/seglen(irecn),.FALSE.)!iprvn+1,irecn+1
                  IF(err/=0) STOP 'Could not insert element into plantmatrix'
                  run=.FALSE.
               ENDIF
               ! definition for the next run of the loop
               ifoln=irecn
               irecn=iprvn
               n_apex=.FALSE.
            END DO ! loop on branch nodes
         END DO ! loop on branches

         PHr_tot1(1:nrec+1)=PHs1(0:nrec)   ! initial xylem PH -> could also be zero but then convergence should take longer (this should be the common case)
         PHr_tot2(1:nrec+1)=PHs2(0:nrec)
         itol=3
         itmaxl=(nrec+1)*3 !maximum iterations for biconjugate gradient-> needed in case of root growth
         ! solve with biconjugated gradient method (Numerical recipies)
         CALL linbcg(nrec+1,Q1(1:nrec+1),PHr_tot1,itol,maxerr,itmaxl,iter_loc,err_loc,ipl)           
         IF (err_loc.GT.maxerr) THEN
            PRINT *,'No convergence in the root system, err_loc > maxerr',err_loc,maxerr
            STOP
         ENDIF
         CALL linbcg(nrec+1,Q2(1:nrec+1),PHr_tot2,itol,maxerr,itmaxl,iter_loc,err_loc,ipl)           
         IF (err_loc.GT.maxerr) THEN
            PRINT *,'No convergence in the root system, err_loc > maxerr',err_loc,maxerr
            STOP
         ENDIF
      
         !====================== end solve root system ====================================
         !***********************************************************************
         ! calculate SINK term
         !**********************************************************************
         ! initialisation
         Jcol1=0.0_dp
         Jcol2=0.0_dp
         DO irecn=1,nrec
            Joutr1(irecn) = Lr(irecn)*segsur(irecn)*(PHs1(irecn)-PHr_tot1(irecn+1)) ![L3/T] (soil is still from 0 to nrec while root is from 1 to nrec+1)
            Joutr2(irecn) = Lr(irecn)*segsur(irecn)*(PHs2(irecn)-PHr_tot2(irecn+1))
            IF (.NOT.lSUF) THEN
               DO isub=1,nsub(irecn,ipl)
                  c_i=cube_i(irecn,isub,ipl)
                  Sink_cube1(c_i)=Sink_cube1(c_i)+Joutr1(irecn)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
                  Sink_cube2(c_i)=Sink_cube2(c_i)+Joutr2(irecn)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
                  IF (.NOT.lSinkCube) THEN
                     DO ic=1,8
                        corner_ic=loc_Q(irecn,ic,isub,ipl)
                        sink1(corner_ic)=sink1(corner_ic)+Joutr1(irecn)*w_sub(irecn,isub,ipl)*w_dis(irecn,ic,isub,ipl) &
                             /sum_dis(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/Wn(corner_ic)
                        sink2(corner_ic)=sink2(corner_ic)+Joutr2(irecn)*w_sub(irecn,isub,ipl)*w_dis(irecn,ic,isub,ipl) &
                             /sum_dis(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)/Wn(corner_ic)
                     END DO
                  ENDIF
               END DO
            ENDIF
            Jcol1=Jcol1+Joutr1(irecn)
            Jcol2=Jcol2+Joutr2(irecn)
         END DO ! nodes loop

         Krs=Jcol1/(-300.0_dp-PHr_tot1(1))
         IF (.NOT.lSUF) THEN
            IF(.NOT.ALLOCATED(SSF)) ALLOCATE (SSF(1:nElm))
            SSF(1:nElm)=Sink_cube1(1:nElm)*dxgrid*dygrid*dzgrid/Jcol1
            IF (.NOT.lSinkCube) THEN
               Beta(1:nx*ny*nz)=sink1(1:nx*ny*nz)/Jcol1
            ENDIF
         ELSE		!> if lSUF (and kout=0)
            ALLOCATE (SSF(1:nrec))
            ALLOCATE (Inv_c1KxKr(1:nrec))
            SSF(1:nrec)=Joutr1(1:nrec)/Jcol1
            Inv_c1KxKr(1:nrec)=SSF(1:nrec)*Krs
         ENDIF
      END DO 		!plant loop
   ELSE		!> if lSUF and kout>0
      DO ipl=1,nplant 				!Actually currently not made to work with several plants
         WRITE (*,*) '... Using Inv_c1 ...'
         PHs(1:nrec,ipl)=0._dp
         DO irecn=1,nrec
            DO isub=1,nsub(irecn,ipl)
	        PHs(irecn,ipl)=PHs(irecn,ipl)+SUM(hnew(loc_Q(irecn,1:8,isub,ipl))*w_dis(irecn,1:8,isub,ipl)/SUM(w_dis(irecn,1:8,isub,ipl)))*w_sub(irecn,isub,ipl)/SUM(w_sub(irecn,1:nsub(irecn,ipl),ipl))
            END DO    !end do-loop over subsegments
            Lr_fact(irecn)=Gap(PHs(irecn,ipl))*AQPc(PHs(irecn,ipl))
            Lr(irecn)=Lr_pot(irecn)*Lr_fact(irecn)
         END DO
         Krs=DOT_PRODUCT(Inv_c1KxKr(1:nrec),Lr_fact(1:nrec))
         SSF(1:nrec)=Inv_c1KxKr(1:nrec)*Lr_fact(1:nrec)/Krs
         SSFcum=SUM(SSF(1:nrec))
      END DO
   ENDIF
      
   print *,'... Calculating Kcomp ...'
   DO ipl=1,nplant 
      IF (.NOT.lSUF.OR.kout.EQ.0) THEN
         IF (.NOT.lSUF) THEN
            ALLOCATE (SSFunstock(1:nElm,1:2))
            SSFunstock(1:nElm,2)=SSF				!> Will progressively unstock SSFi
            DO i=1,nElm
               SSFunstock(i,1)=i					!> Contains element numbers
               Hs_cube2(i)=-150.0_sp+(zGrid(i)+zGrid(i+nx*ny)) !> Works for both XY periodic and non-periodic domains
            END DO
            nk=nElm							!> \param nk=number of SSFi stocked in SSFunstock
            SSFt=0.0001_dp						
            !> \param SSFt Threshold value above which SSFi is considered in the estimation of Kcomp
            SSFcum=0.0_dp						!> \param SSFcum Cumulation of "accepted" SSFi
            j=0							!> j=number of SSFi already cumulated in SSFcum
            DO WHILE (SSFcum.LT.0.9_dp)
               SSFt=0.5_dp*SSFt+0.5_dp*SUM(SSFunstock(1:nk,2))/nk
               ALLOCATE (SSFunstockNew(1:nk,1:2))
               SSFunstockNew=0.0_dp
               k=0							!> k=number of SSFi stocked in SSFunstockNew
               DO i=1,nk
                  IF (SSFunstock(i,2).GE.SSFt) THEN
                     SSFcum=SSFcum+SSFunstock(i,2)
                     j=j+1
                     Phi=Sink_cube2(INT(SSFunstock(i,1)))*dxgrid*dygrid*dzgrid/SSFunstock(i,2)-Jcol2
                     Hsi=Hs_cube2(INT(SSFunstock(i,1)))+200.0_sp
                     sumPhiHsi=sumPhiHsi+Phi*Hsi
                     sumHsi=sumHsi+Hsi
                     sumHsi2=sumHsi2+Hsi**2
                     sumPhi=sumPhi+Phi
                     sumPhi2=sumPhi2+Phi**2
                  ELSEIF (SSF(i).GT.0.0_dp) THEN
                     SSFunstockNew(k+1,1)=SSFunstock(i,1)
                     SSFunstockNew(k+1,2)=SSFunstock(i,2)
                     k=k+1
                  ENDIF
               END DO
               nk=k
               DEALLOCATE (SSFunstock)
               ALLOCATE (SSFunstock(1:nk,1:2))
               SSFunstock=SSFunstockNew(1:nk,1:2)
               DEALLOCATE (SSFunstockNew)
            END DO
            Kcomp=(j*sumPhiHsi-sumHsi*sumPhi)/(j*sumHsi2-sumHsi**2)
            r=(sumPhiHsi-sumHsi*sumPhi/j)/SQRT((sumHsi2-sumHsi**2/j)*(sumPhi2-sumPhi**2/j))
         ELSE	!if lSUF
            ALLOCATE (ind(1:nrec))
            Hseq=DOT_PRODUCT(SSF(1:nrec),PHs2(1:nrec))
            k=0
            DO i=1,nrec
               IF (SSF(i).GT.0.5/nrec) THEN
                  PhiSUF(i)=Joutr2(i)/SSF(i)-Jcol2
                  k=k+1
                  ind(k)=i
               ELSE
                  PhiSUF(i)=0.0
               ENDIF
            END DO
            Kcomp=DOT_PRODUCT(PHs2(ind(1:k))-Hseq,PhiSUF(ind(1:k)))/DOT_PRODUCT(PHs2(ind(1:k))-Hseq,PHs2(ind(1:k))-Hseq)
            r=(DOT_PRODUCT(PhiSUF(ind(1:k)),PHs2(ind(1:k)))-SUM(PHs2(ind(1:k)))*SUM(PhiSUF(ind(1:k)))/k)/SQRT((DOT_PRODUCT(PhiSUF(ind(1:k)),PhiSUF(ind(1:k)))-SUM(PhiSUF(ind(1:k)))**2/k)*(DOT_PRODUCT(PHs2(ind(1:k)),PHs2(ind(1:k)))-SUM(PHs2(ind(1:k)))**2/k))
            SSFcum=SUM(SSF(ind(1:k)))
         ENDIF
         WRITE (*,*) 'Krs [cm³/d]',Krs
         WRITE (*,*) 'Kcomp [cm³/d]',Kcomp
         WRITE (*,*) 'Correlation coefficient [-]',r
         WRITE (*,*) 'Percentage of cumulated SSF considered in the calculation of Kcomp [%]',SSFcum*100
      ELSE		!if lSUF and kout>0
         PHs(1:nrec,ipl)=0._dp
         DO irecn=1,nrec
            DO isub=1,nsub(irecn,ipl)
	        PHs(irecn,ipl)=PHs(irecn,ipl)+SUM(hnew(loc_Q(irecn,1:8,isub,ipl))*w_dis(irecn,1:8,isub,ipl)/SUM(w_dis(irecn,1:8,isub,ipl)))*w_sub(irecn,isub,ipl)/SUM(w_sub(irecn,1:nsub(irecn,ipl),ipl))
            END DO    !end do-loop over subsegments
            Lr_fact(irecn)=Gap(PHs(irecn,ipl))*AQPc(PHs(irecn,ipl))
            Lr(irecn)=Lr_pot(irecn)*Lr_fact(irecn)
         END DO
         Krs=DOT_PRODUCT(Inv_c1KxKr(1:nrec),Lr_fact(1:nrec))
         SSF(1:nrec)=Inv_c1KxKr(1:nrec)*Lr_fact(1:nrec)/Krs
         Kcomp=Krs
      ENDIF

      !> writes Couvreur.out
      WRITE (file,'(A12)')'out/Couvreur'
      WRITE (file(13:13),'(I1)') ipl
      WRITE (file(14:14),'(A1)') '.'
      IF (kOut.LT.10) THEN
         WRITE (file(15:15),'(I1)') kOut
      ELSEIF (kOut.LT.100) THEN
         WRITE (file(15:16),'(I2)') kOut
      ELSE
         WRITE (file(15:17),'(I3)') kOut
      ENDIF
      OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
      WRITE (8,'(''********** COUVREUR ET AL. (2012) ROOT WATER UPTAKE MODEL PARAMETERS **********'')')
      WRITE (8,*)
      WRITE (8,'(''Stress function ? (1=yes,2=no)/ Threshold value of the stress function [hPa]'')')
      IF (stresfun.EQ.1) THEN
         WRITE (8,'(I1,3X,F10.2)') 1,hx_min
      ELSE
         WRITE (8,'(I1,3X,I1)') 2,0	!> non isohydric plants (Tuzet, etc. cannot currently be considered when using Couvreur RWU model)
      ENDIF
      WRITE (8,*)
      WRITE (8,'(''Dimension of the root water uptake model (1=soil domain,2=root system domain)'')')
      IF (lSUF) THEN
         WRITE (8,'(I1)') 2
      ELSE
         WRITE (8,'(I1)') 1
      ENDIF
      WRITE (8,*)
      WRITE (8,'(''Use of modified de Jong van Lier function for Hseq prediction ? (1=yes,2=no)'')')
      IF (ldJvL) THEN
         WRITE (8,'(I1)') 1
      ELSE
         WRITE (8,'(I1)') 2
      ENDIF
      WRITE (8,*)
      WRITE (8,'(''Equivalent conductance of the whole root system [cm³/hPa/day]'')')
      WRITE (8,'(F10.8)') Krs
      WRITE (8,*)
      WRITE (8,'(''Compensatory RWU conductance of the root system [cm³/hPa/day]'')')
      WRITE (8,'(F10.8)') Kcomp
      WRITE (8,*)
      IF (lSUF) THEN
         WRITE (8,'(''* Standard Uptake Fraction distribution'')')
         WRITE (8,'(''Number of root segments [-]'')')
         WRITE (8,'(I6)') nrec
         WRITE (8,'(''  SegID#    SUF'')')
         DO i=1,nrec
               WRITE (8,'(I6,3X,F12.10)') i,SSF(i)		!> Note that both SUF and SSF have the same name (SSF) in the script since they cannot be used in the same time
         ENDDO
         CLOSE (8)
      ELSEIF (lSinkCube) THEN
         WRITE (8,'(''* Standard Sink Fraction distribution'')')
         WRITE (8,'(''Number of elements in each dimension [-] and elements size [cm]'')')
         WRITE (8,'(I5,1X,I5,1X,I5,3(1X,F10.3))') nex,ney,nez,dxGrid,dyGrid,dzGrid
         WRITE (8,'(''  Element#    SSF'')')
         DO i=1,nElm
               WRITE (8,'(I6,3X,F12.10)') i,SSF(i)
         ENDDO
         CLOSE (8)
      ELSE
         WRITE (8,'(''* Standard Sink Fraction distribution'')')
         WRITE (8,'(''Number of nodes in each dimension [-] and elements size [cm]'')')
         WRITE (8,'(I5,1X,I5,1X,I5,3(1X,F10.3))') nx,ny,nz,dxGrid,dyGrid,dzGrid
         WRITE (8,'(''  Node#    BetaSSF'')')
         DO i=1,nx*ny*nz
            WRITE (8,'(I6,3X,F14.10)') i,Beta(nx*ny*nz+1-i)		!> Nodes # inverted for Hydrus
         ENDDO
         CLOSE (8)
         STOP 'SSF generated for Hydrus!'
      ENDIF
   ENDDO ! plant loop
   RETURN
 END SUBROUTINE SSFdis
!***************************************************************************
 SUBROUTINE analytical_approach(ipl)
   !      USE OMP_LIB
   USE TypeDef
   USE Paramdata,ONLY: MaxNod, pi,lSalinity
   USE GridData, ONLY: xgrid,ygrid,zgrid,dxgrid,dygrid,dzgrid,elmnod,nPt,Axy,Bxy,Dxy,Exy
   USE DoussanMat, ONLY :old,oldT,ave,eqDis,nsub,PHs,loc_q,cp_mean,voxel_no,no_voxels, &
        voxel_node,numnodes_voxel,h_mat2,hx_min,phr_sub,n,indexvalue,counter2,phi_mat, &
        PH_micro2,ph_root_sub,Qd,w_dis,w_sub,cent,Intc,PHr,Lr,sum_dis,switchcriterion,Qi,Q_bc, &
        tcheck,tcheck2,count_nodes,k_ave,B,Q_bc1,Q_bc2,cube_i,PHs_osmotic,PHo
   USE RootData, ONLY: nrec,seglen,segsur,nbr,ibrseg,irecpr,sigma
   USE SolData, Only: hnew, hold, matnum
   IMPLICIT NONE

   INTEGER(ap) :: checker1=0,blaat,irecnTemp,isubTemp,ki
   REAL(sp) :: beta
   REAL(sp) :: y(1:n),h_mat(1:n),arr(1:n-1),Phi(1:n-1)
   REAL(sp),DIMENSION(1:8) ::h_matK,temp2
   REAL(sp) :: K_k,C_k,theta_k,temp
   REAL(sp) :: cumsum(1:n-1),Phi_rend
   REAL(sp) :: z1,z2,frac1,frac2,xCent,yCent,Zcent,xInt,yInt,zInt
   REAL(sp) :: testx(2),testy(2),testz(2),xgrid_max,xgrid_min,ygrid_max,ygrid_min,mXY(2),point(1:4,1:2)
   REAL(sp) :: zgrid_max,zgrid_min,mZX(2),mZY(2)
   REAL(sp) :: PH_micro(1:n-1)
   REAL(sp) :: theta(1:n),K(1:n),C(1:n),rho,r0,r(1:n-1),rend
   REAL (dp) ::sumW,psub,q_root,q_out
   REAL (dp) :: pos_cent,testvar
   REAL(dp),DIMENSION(1:8) ::pond_i,q_tot_c
   REAL(sp),DIMENSION(1:4,1:4) :: dis,Flowpoint_rel,PHpoint_rel
   REAL(sp),DIMENSION(1:4) :: mFlows,mPHs,PH_inter,Flow_inter
   REAL(sp) :: meanPHS,meanFlows
   INTEGER(ap) ::j,counter,l,w,Flow_corner(1:8,1:6)
   INTEGER(ap) ::corner(1:8),var1,kl,imin,ipl
   INTEGER(ap) :: irecn,i,isub,pos,ja=0,irecn_tmp,isub_tmp
   REAL(sp) :: z_angle,x_angle,y_angle,x1,x2,y1,y2,q_root_tmp(1:indexValue),pond_subTemp
   LOGICAL :: zalign,yalign,xalign,extra_checker=.FALSE.,countja

   !> parameters for counting nodes that are under limiting conditions (ks < lr*r0)
   countja=.FALSE.
   count_nodes=0
   DO irecn=0,nrec !> loop over segments
      old=oldT
      x_angle=0
      y_angle=0
      z_angle=0
      !> initialize soil PH matrix
      PHs(irecn,ipl)=0._dp
      DO isub=1,nsub(irecn,ipl) !> loop over subsegments
         dis=0
         Flowpoint_rel=0
         PHpoint_rel=0
         mFlows=0
         mPHs=0
         PH_inter=0
         Flow_inter=0
         !> ---------------------------calculate the flux in each soil node with Darcy´s law--------
         IF (.NOT.eqDis) THEN
            DO i=1,8 !> for each corner of a cube
               corner(i)=loc_Q(irecn,i,isub,ipl)
            END DO
            !> determine nodes surrounding a soil node to determine the flow
            CALL FlowNeighb(corner,Flow_corner)
            CALL flux_node(Flow_corner,corner,q_tot_c)
         ELSEIF (eqDis) THEN !method C
            q_tot_c(1:8)=0
         ENDIF
         !> ---------------------------end flow calculation-----------------------------------
         !> initialization and r0 definement
         checker1=0 !> reset for each cp
         DO i=1,8
            pond_i(i)=w_dis(irecn,i,isub,ipl)
            corner(i)=loc_Q(irecn,i,isub,ipl)
         END DO
         zalign=.FALSE.
         yalign=.FALSE.
         xalign=.FALSE.
         !> not always true: adjust accordingly
         IF (irecn .EQ. 0) THEN
            r0=5e-2 !set to minimal radius for root collar ID only-> not important
         ELSE
            IF (segsur(irecn) .EQ. 0 .OR. seglen(irecn) .EQ. 0) THEN
               r0 = 5e-2 !> minimal root radius (in cm)
            ELSE
               r0 = 1/(2*pi)*segsur(irecn)/seglen(irecn)
            END IF
         ENDIF
         !-----------------------------------------------------------------------
         !
         !> REDUCTION from 3D -> 2D plane, currently averaging procedure is taken
         !> in axial direction (i.e. direction perpendicular to radial flow)  which is sufficient (checked)
         !----------------------------------------------------------------------
         !> for the one_root procedure; the center point of a segment is the middle of a soil voxel
         IF (ave) THEN
            xCent = cp_mean(irecn,1,isub,ipl)
            yCent = cp_mean(irecn,2,isub,ipl)
            zCent = cp_mean(irecn,3,isub,ipl)
            imin = voxel_no(irecn,isub) !voxel number
            IF (counter2(imin) .NE. 0) THEN !> check of conductivity drop in this voxel is already calculated
               irecnTemp = voxel_node(imin,1,counter2(imin))
               isubTemp = voxel_node(imin,2,counter2(imin))
               pond_subTemp = w_sub(irecnTemp,isubTemp,ipl)
               extra_checker=.TRUE.
               GOTO 210 !> to calculate analytical solution only once for the large root node
            ELSEIF (counter2(imin) .EQ. 0 .AND. no_voxels(imin) .NE. 0) THEN
               DO blaat=1,no_voxels(imin)
                  IF (voxel_node(imin,1,blaat) .EQ. irecn .AND. voxel_node(imin,2,blaat) &
                       .EQ. isub .AND. w_sub(irecn,isub,ipl).NE.0) counter2(imin)=blaat
               ENDDO
            ENDIF
         ELSEIF (.NOT.(ave)) THEN !> if not one root method each root node is a center point
            xCent = cent(irecn,1,isub,ipl)
            yCent = cent(irecn,2,isub,ipl)
            zCent = cent(irecn,3,isub,ipl)
         ENDIF
         !> interception point x,y,zInt with soil cube
         xInt = Intc(irecn,1,isub,ipl)
         yInt = Intc(irecn,2,isub,ipl)
         zInt = Intc(irecn,3,isub,ipl)
         !> ---------------------------derive x,y,z-allignment and segment angle------------------
         !> Get aligment of a root in the current voxel. If an intersection point of a root with the cube is 
         !> identical to the center point of the root, pick a direction dependent on the surface it crosses
         !> otherwise calculate the segment angle. first from a z-directional point of view, if the branching 
         !> angle is smaller then a given range then look if it is
         !> aligned in x-direction (using a projection of the z-coordinate on the 2D generated plain)
         !> if this angle doesnt comply to a criterion then the allignement is performed in y-direction
         !> weigh factor Wb is calculated but not needed in remainder calculations
         !> dimensions of the voxel, maxima and minima in each direction
         zgrid_max = MAXVAL(zgrid(corner(1:8)))
         zgrid_min = MINVAL(zgrid(corner(1:8)))
         xgrid_max = MAXVAL(xgrid(corner(1:8)))
         xgrid_min = MINVAL(xgrid(corner(1:8)))
         ygrid_max = MAXVAL(ygrid(corner(1:8)))
         ygrid_min = MINVAL(ygrid(corner(1:8)))
         !> get alignment direction
         IF (xCent-xInt .EQ. 0 .AND. yCent-yInt .EQ. 0 .AND. zCent-zInt .EQ. 0) THEN
            IF (zCent .EQ. zgrid_max .OR. zCent .EQ. zgrid_min) THEN
               zalign=.TRUE.
            ELSEIF (yCent .EQ. ygrid_max .OR. yCent .EQ. ygrid_min) THEN
               yalign=.TRUE.
            ELSEIF (xCent .EQ. xgrid_max .OR. xCent .EQ. xgrid_min) THEN
               xalign=.TRUE.
            ENDIF
         ELSE
            !> geometry; beta is radial angle betwen z difference and the length of the segment-> convert to degree
            beta = ASIN(ABS(zInt-zCent)/SQRT((xInt-xCent)*(xInt-xCent)+ &
                 (yInt-yCent)*(yInt-yCent)+(zInt-zCent)*(zInt-zCent)))
            z_angle=ABS(beta/(2*pi))*360
            IF (z_angle .LE. 90.1 .AND. z_angle .GE. 45) THEN !> certain error deviatin given to criterion
               zalign=.TRUE.
            ELSE
               zInt=zCent !> projection in 2d plane to get angle in x-or y direction
               beta = ACOS(ABS(xInt-xCent)/SQRT((xInt-xCent)*(xInt-xCent)+ &
                    (yInt-yCent)*(yInt-yCent)+(zInt-zCent)*(zInt-zCent)))
               x_angle=ABS(beta/(2*pi))*360
               x_angle=90-x_angle !> becoz of cosine, angle interpreted from 0 to 45 and angle has to be adjusted for weight factor
               IF (x_angle .LE. 90 .AND. x_angle .GE. 45) THEN
                  xalign=.TRUE.
               ELSE
                  yalign=.TRUE.
               ENDIF
            ENDIF
         ENDIF
         !> -----------------------------------end allignment-----------------------------------
         !> for each allocated root direction in a voxel a checker1 is set which determines where the cp is located, 
         !> in the voxel or at the edges. If it is located in the  ! voxel then the analytical solution is applied, 
         !> otherwise the old method is sufficient
         !> set checker1
         IF (zalign) THEN
            IF (ABS(xCent-xgrid_max) .GT. r0 .AND. ABS(xCent-xgrid_min).GT.r0 &
                 .AND. ABS(yCent-ygrid_max) .GT. r0 .AND. ABS(yCent-ygrid_min).GT.r0 ) THEN !in plane
               checker1 = 1
            ELSE !> if it is at a boundary -> apply a simple method to obtain PH at the soil-root interface (for the moment)
               checker1=0
            END IF
         ELSEIF (xalign) THEN
            IF (ABS(yCent-ygrid_max) .GT. r0 .AND. ABS(yCent-ygrid_min).GT. r0 &
                 .AND. ABS(zCent-zgrid_max) .GT. r0 .AND. ABS(zCent-zgrid_min).GT.r0 ) THEN !in plane
               checker1 = 1
            ELSE !> if it is at a boundary -> apply a simple method to obtain PH at the soil-root interface (for the moment)
               checker1=0
            END IF
         ELSEIF (yalign) THEN
            IF (ABS(xCent-xgrid_max) .GT. r0 .AND. ABS(xCent-xgrid_min).GT.r0 &
                 .AND. ABS(zCent-zgrid_max) .GT. r0 .AND. ABS(zCent-zgrid_min).GT.r0 ) THEN !in plane
               checker1 = 1
            ELSE !> if it is at a boundary -> apply a simple method to obtain PH at the soil-root interface (for the moment)
               checker1=0
            END IF
         ENDIF
!> Average for the specific direction chosen from a 3D to a 2D mesh (for application of the 2D analytical solution). 
!> In direction plane is orientated
!> To be clearer; using the allignment of the root node the PH and Flow boundary conditions
!> are transferred to the 4 corner nodes spanning the 2D domain in which the
!> analytical approach is applied
!>
!>              C1----------------C2
!>		|		  |
!>		|		  |
!>		|		  |
!>		|	o	  |
!>		|		  |
!>		|		  |
!>		|		  |
!>		C3----------------C4
!>
!> z-align, x-align, y-align are discussed; basicly the same routines;
         IF (checker1.EQ.1 .AND. (zalign)) THEN
            DO i=1,4
               z1 = zgrid(corner(i))
               z2 = zgrid(corner(i+4))
               frac1 = 1/ABS(zCent-z1)
               frac2 = 1/ABS(zCent-z2)
               IF (ABS(zCent-z1) .EQ. 0) THEN
                  PH_inter(i) = hnew(corner(i))
                  !> ph bound. cond.
                  Flow_inter(i) = q_tot_c(i)
                  !> flow bound. cond.
               ELSEIF (ABS(zCent-z2) .EQ. 0) THEN
                  PH_inter(i) = hnew(corner(i+4))
                  Flow_inter(i) = q_tot_c(i+4)
               ELSE
                  PH_inter(i) = (frac1*hnew(corner(i)) + frac2*hnew(corner(i+4)))/(frac1+frac2)
                  Flow_inter(i) = (frac1*q_tot_c(i) + frac2*q_tot_c(i+4))/(frac1+frac2)
                  !> corner -> inter (2D plane)
               ENDIF
            ENDDO
            testx(1) = ABS(xgrid_max-xCent)
            testx(2) = ABS(xgrid_min-xCent)
            mXY(1) = MINVAL(testx) !> minimum distance in x-dimension
            testy(1) = ABS(ygrid_max-yCent)
            testy(2) = ABS(ygrid_min-yCent)
            mXY(2) = MINVAL(testy) !> minimum distance in y-direction
!> get minimum radius from root node to a 2D edge; to be used as outer radius
!> for the analytical approach
            rend = MINVAL(mXY)     !> outer radius of the soil cylinder
            IF (eqDis) THEN !> get an equal radius dependent on the number of root nodes in this voxel
               testvar=numNodes_voxel(irecn,isub)**(1d0/3d0)
               rend = 1/2d0*SQRT((dzgrid)/FLOOR(testvar))
            ENDIF
            IF (rend .LT. r0) THEN
               ja=1
               GOTO 203
            ENDIF
!> the dashed arrow line is set as the minimal distance in this example; i.e. rend (outer radius)
!>             C1------0---------C2
!>		|	|	  |
!>		|	|	  |
!>		|	|	  |
!>		0<----->o---------0
!>		|	|	  |
!>		|	|	  |
!>		|	|	  |
!>		C3------0---------C4
!> new coordinate system 2D based on the radius of the microscopic model (question is, how large is ratio rend/r0)
!> so the outer radius has intersection points with the 2D plane as drawn above (the zeros´s (0) are the
!> intersection points with the 2D plane
!> these intersection points are calculated here; for method B,C (see MS)
            IF (eqDis) THEN !> middle of voxel as point of departure for bnd.cond.
               point(1,1) = (xgrid_max+xgrid_min)/2+rend
               point(1,2) = (ygrid_max+ygrid_min)/2
               point(2,1) = (xgrid_max+xgrid_min)/2
               point(2,2) = (ygrid_max+ygrid_min)/2+rend
               point(3,1) = (xgrid_max+xgrid_min)/2-rend
               point(3,2) = (ygrid_max+ygrid_min)/2
               point(4,1) = (xgrid_max+xgrid_min)/2
               point(4,2) = (ygrid_max+ygrid_min)/2-rend
            ELSE
               point(1,1) = xCent+rend
               point(1,2) = yCent
               point(2,1) = xCent
               point(2,2) = yCent+rend
               point(3,1) = xCent-rend
               point(3,2) = yCent
               point(4,1) = xCent
               point(4,2) = yCent-rend
            ENDIF
!> calculate the pressure head, based on a distance average, for each new coordinate surrounding the 
!> CP(center ponint, i.e. the root node)---- 2D!!!!!!!! z-alignment-> x and y-direction is evaluated
!> so outer boundary conditions are mapped on the 4 points that represent the outer circle
!> The new PH in each intersection point (zero) with 2D domain is mPHs, here only part of that solution is given
!> due to distance saciling!! so "_rel" is intermediate value for PH or Flow (later these values are put to 
!> the correct valuesin mPHs,mFlow
            DO j=1,SIZE(point,1)
               DO i=1,4
                  dis(j,i) = SQRT((point(j,1)-xgrid(corner(i)))**2+(point(j,2)-ygrid(corner(i)))**2)
                  PHpoint_rel(j,i) = 1/dis(j,i)*PH_inter(i) !> inter values being the values in the 4 corner nodes of the 2D domain
                  Flowpoint_rel(j,i) = 1/dis(j,i)*Flow_inter(i)
                  !> inter -> point (distance)
               END DO
            END DO
         ELSEIF (checker1.EQ.1 .AND. (xalign)) THEN
            kl=1
            DO i=1,7,2
               x1 = xgrid(corner(i))
               x2 = xgrid(corner(i+1))
               frac1 = 1/ABS(xCent-x1)
               frac2 = 1/ABS(xCent-x2)
               IF (ABS(xCent-x1) .EQ. 0) THEN
                  PH_inter(kl) = hnew(corner(i))
                  Flow_inter(kl) = q_tot_c(i)
               ELSE IF (ABS(xCent-x2) .EQ. 0) THEN
                  PH_inter(kl) = hnew(corner(i+1))
                  Flow_inter(kl) = q_tot_c(i+1)
               ELSE
                  PH_inter(kl) = (frac1*hnew(corner(i)) + frac2*hnew(corner(i+1)))/(frac1+frac2)
                  Flow_inter(kl) = (frac1*q_tot_c(i) + frac2*q_tot_c(i+1))/(frac1+frac2)
                  !> corner -> inter (2D plane)
               ENDIF
               kl=kl+1
            END DO
            testz(1) = ABS(zgrid_max-zCent)
            testz(2) = ABS(zgrid_min-zCent)
            mZY(1) = MINVAL(testz) !> minimum distance in z-dimension
            testy(1) = ABS(ygrid_max-yCent)
            testy(2) = ABS(ygrid_min-yCent)
            mZY(2) = MINVAL(testy) !> minimum distance in y-direction
            rend = MINVAL(mZY)     !> outer radius of the soil cylinder
            IF (eqDis) THEN
               testvar=numNodes_voxel(irecn,isub)**(1d0/3d0)
               rend = 1/2d0*SQRT((dxgrid)/FLOOR(testvar))
            ENDIF
            IF (rend .LT. r0) THEN
               ja=1
               GOTO 203
            ENDIF
            !> new coordinate system 2D based on the radius of the microscopic model (question is, how large is ratio rend/r0)
            IF (eqDis) THEN !> middle of voxel as point of departure for bnd.cond.
               point(1,1) = (zgrid_max+zgrid_min)/2+rend
               point(1,2) = (ygrid_max+ygrid_min)/2
               point(2,1) = (zgrid_max+zgrid_min)/2
               point(2,2) = (ygrid_max+ygrid_min)/2+rend
               point(3,1) = (zgrid_max+zgrid_min)/2-rend
               point(3,2) = (ygrid_max+ygrid_min)/2
               point(4,1) = (zgrid_max+zgrid_min)/2
               point(4,2) = (ygrid_max+ygrid_min)/2-rend
            ELSE
               point(1,1) = zCent+rend
               point(1,2) = yCent
               point(2,1) = zCent
               point(2,2) = yCent+rend
               point(3,1) = zCent-rend
               point(3,2) = yCent
               point(4,1) = zCent
               point(4,2) = yCent-rend
            ENDIF
            DO j=1,SIZE(point,1)
               DO i=1,4
                  IF (i.EQ.1) kl=1 !> y and z-coord of these nodes needed
                  IF (i.EQ.2) kl=3
                  IF (i.EQ.3) kl=5
                  IF (i.EQ.4) kl=7
                  dis(j,i) = SQRT((point(j,1)-zgrid(corner(kl)))**2+(point(j,2)-ygrid(corner(kl)))**2)
                  PHpoint_rel(j,i) = 1/dis(j,i)*PH_inter(i)
                  Flowpoint_rel(j,i) = 1/dis(j,i)*Flow_inter(i)
                  !> inter -> point (distance)
               END DO
            END DO
         ELSEIF (checker1.EQ.1 .AND. (yalign)) THEN
            kl=1
            DO i=1,4
               var1=i
               IF (i.GT.2) var1=var1+2
               y1 = xgrid(corner(var1))
               y2 = xgrid(corner(var1+2))
               frac1 = 1/ABS(yCent-y1)
               frac2 = 1/ABS(yCent-y2)
               IF (ABS(yCent-y1) .EQ. 0) THEN
                  PH_inter(kl) = hnew(corner(var1))
                  Flow_inter(kl) = q_tot_c(var1)
               ELSEIF (ABS(yCent-y2) .EQ. 0) THEN
                  PH_inter(kl) = hnew(corner(var1+2))
                  Flow_inter(kl) = q_tot_c(var1+2)
               ELSE
                  PH_inter(kl) = (frac1*hnew(corner(var1)) + frac2*hnew(corner(var1+2)))/(frac1+frac2)
                  Flow_inter(kl) = (frac1*q_tot_c(var1) + &
                       frac2*q_tot_c(var1+2))/(frac1+frac2)
                  !> corner -> inter (2D plane)
               ENDIF
               kl=kl+1
            ENDDO
            testz(1) = ABS(zgrid_max-zCent)
            testz(2) = ABS(zgrid_min-zCent)
            mZX(1) = MINVAL(testz) !> minimum distance in z-dimension
            testx(1) = ABS(xgrid_max-xCent)
            testx(2) = ABS(xgrid_min-xCent)
            mZX(2) = MINVAL(testx) !> minimum distance in y-direction
            rend = MINVAL(mZX)
            IF (eqDis) THEN
               testvar=numNodes_voxel(irecn,isub)**(1d0/3d0)
               rend = 1/2d0*SQRT((dygrid)/FLOOR(testvar))
            ENDIF
            IF (rend .LT. r0) THEN
               ja=1
               GOTO 203
            ENDIF
            !> new coordinate system 2D based on the radius of the microscopic model (question is, how large is ratio rend/r0)
            IF (eqDis) THEN !> middle of voxel as point of departure for bnd.cond.
               point(1,1) = (zgrid_max+zgrid_min)/2+rend
               point(1,2) = (xgrid_max+xgrid_min)/2
               point(2,1) = (zgrid_max+zgrid_min)/2
               point(2,2) = (xgrid_max+xgrid_min)/2+rend
               point(3,1) = (zgrid_max+zgrid_min)/2-rend
               point(3,2) = (xgrid_max+xgrid_min)/2
               point(4,1) = (zgrid_max+zgrid_min)/2
               point(4,2) = (xgrid_max+xgrid_min)/2-rend
            ELSE
               point(1,1) = zCent+rend
               point(1,2) = xCent
               point(2,1) = zCent
               point(2,2) = xCent+rend
               point(3,1) = zCent-rend
               point(3,2) = xCent
               point(4,1) = zCent
               point(4,2) = xCent-rend
            ENDIF
            DO j=1,SIZE(point,1)
               DO i=1,4
                  IF (i.EQ.1) kl=1
                  IF (i.EQ.2) kl=2
                  IF (i.EQ.3) kl=5
                  IF (i.EQ.4) kl=6
                  dis(j,i) = SQRT((point(j,1)-zgrid(corner(kl)))**2+(point(j,2)-xgrid(corner(kl)))**2)
                  PHpoint_rel(j,i) = 1/dis(j,i)*PH_inter(i)
                  Flowpoint_rel(j,i) = 1/dis(j,i)*Flow_inter(i)
                  !> inter -> point (distance)
               END DO
            END DO
         ENDIF
         !> IF CP is located in the plane (most of the cases)
210      CONTINUE   !> one root approach continues here if applicable
         IF (checker1 .EQ. 1 .OR. (extra_checker)) THEN
            IF (extra_checker) THEN
               extra_checker=.FALSE.
               GOTO 200
            ENDIF
!-----------------------------------------------------------------------------------------------------------------
!> Determine the soil cylinder
!> -----------------------------------------------------------------------------------------------------------------
!> after outer boundary conditions at the new circle surrounding the cp are known for each allignment, 
!> calculate the average PH and flow for these boundary conditions
!> these average PH and Flow values (meanPHs and meanFlows) at the outer radius are the boundary conditions
!> for the analytical approach
            DO j=1,4
               mPHs(j) = SUM(PHpoint_rel(j,:))/SUM(1/dis(j,:))
               mFlows(j) = SUM(Flowpoint_rel(j,:))/SUM(1/dis(j,:))
            END DO
            meanPHs = (SUM(mPHs(1:4))/4)
            meanFlows = (SUM(mFlows(1:4))/4)
            !> fractioned value should be used (subsegmentation)
            meanPHs = meanPHs*w_sub(irecn,isub,ipl)
            meanFlows = meanFlows*w_sub(irecn,isub,ipl)
!-----------------------------------------------------------------------------------------------------------------
!> partitioning the radii used for the anlytical approach (in log scale); closer
!> radii near the root, larger radii at outer radius
!> -----------------------------------------------------------------------------------------------------------------
            !> spatial grid
            CALL logspace(r0,rend,n-1,y)
            r=y(1:n-1)
            rho = rend/r0
!-----------------------------------------------------------------------
!
!> Calculate BOUNDARY conditions for analytical solution: Phi_rend(based on meanPHs),q_end,q_root
!>
!> -----------------------------------------------------------------------------------------------------------------
!>  Calculate Phi_rend -> obtained from the PH boundary condition (the PH at outer radius of microscopic model)
!> -----------------------------------------------------------------------------------------------------------------
            CALL logspace(hx_min,meanPHs,n,y)
            h_mat = -y
            CALL SetMat_ana(h_mat,theta,K,C)
            !> calculate matric flux potential integral K(h)dh -> numerical integration
            !> cumsum = cumulative sum
            arr=ABS(h_mat(2:SIZE(h_mat))-h_mat(1:SIZE(h_mat)-1))*(K(2:SIZE(K))+K(1:SIZE(K)-1))/2
            DO i=1,SIZE(arr)
               cumsum(i) = SUM(arr(1:i))
            END DO
            !> last value from integral is at r=rend
            Phi_rend=cumsum(SIZE(cumsum))
!--------------q_root------------------------
!> get current voxel number; for a loop over all the nodes in current voxel, get their irecn and isub values; calculate q_root for all root nodes
            IF (ave) THEN
               imin = voxel_no(irecn,isub)
               DO ki=1,no_voxels(imin)
                  irecn_tmp = voxel_node(imin,1,ki)
                  isub_tmp = voxel_node(imin,2,ki)
                  PHr_sub(irecn_tmp,isub_tmp) = w_sub(irecn_tmp,isub_tmp,ipl)*PHr(irecn_tmp+1,ipl)
                  q_root_tmp(ki)= Lr(irecn_tmp) * (PH_root_sub(irecn_tmp,isub_tmp) &
                       - PHr_sub(irecn_tmp,isub_tmp))
               ENDDO
!===============================================================================
!> chose the way how q_root is defined:
!> 1) sum of all root flows, consider as one sucking root
!> PHr_sub is not used, then only to check wheter or not limiting conditions have been reached -> see PHi
               PHr_sub(irecn,isub) = w_sub(irecn,isub,ipl)*PHr(irecn+1,ipl)
               !1)
               q_root = SUM(q_root_tmp(1:ki)) !> total uptake -> one big root
               !!-------------------------------------------------------------------------------
            ELSEIF (.NOT.(ave)) THEN
               !> PHxylem = PHr; PHr_sub = subsegment xylemPH
               PHr_sub(irecn,isub) = w_sub(irecn,isub,ipl)*PHr(irecn+1,ipl)
               !> ALCULATE BOUNDARY CONDITION q_root
               !> should be fractioned
               q_root= Lr(irecn) * (PH_root_sub(irecn,isub) - PHr_sub(irecn,isub))
            ENDIF
            !> CALCULATE BOUNDARY CONDITION q_out
            !> steady state -> condition for doussan matrix
            IF (meanFlows .LE. q_root*(r0/rend)) THEN
               q_out = meanFlows
            ELSE
               q_out = q_root*(r0/rend)
            ENDIF
            !de willigen, q_out=zero
            !if (eqDis) q_out=0
!**************************************************************************
!> ------------------------------ calculate PHI -----------------------------
!> **************************************************************************
!-----------------------------------------------------
!> actual PHi calculation; hx_min equals lower boundary integral from flux density -> Phi_r_root = 0
!> -----------------------------------------------------
            IF (PHr_sub(irecn,isub) .LE. hx_min*w_sub(irecn,isub,ipl)) THEN 
               !> if PHsoil-root interface = PHxylem limiting -> flow is zero (ultimate case limiting root PH
               Phi = ((Phi_rend + q_out*rend/2*LOG(1/rho) ) / &
                    ( rho**2-1 + 2*rho**2*LOG(1/rho) ) ) * (r**2/r0**2-1 + 2*rho**2*LOG(r0/r)) &
                    + rend*q_out*LOG(r/r0)
            ELSE !> always true, also when hx_min is reached for PHxylem. If finally PHsoil-root inter = PHxylem (=hx_min), see above
               Phi = Phi_rend + (q_root*r0 - q_out*rend)* &
                    (r**2/r0**2/(2*(1-rho**2)) + rho**2/(1-rho**2)*(LOG(rend/r)-1/2.)) &
                    + q_out*rend*LOG(r/rend)
            END IF
            counter=1
            !> determine PH from Phi in microscopic domain
            DO l=1,SIZE(Phi)
               DO w=counter,SIZE(Phi_mat)
                  IF (Phi_mat(w) .LT. Phi(l)) THEN
                     counter = w;
                  ELSEIF (Phi_mat(w) .GT. Phi(l)) THEN
                     EXIT
                  END IF
               END DO
               PH_micro(l)=h_mat2(counter)
            END DO
            PH_micro2(irecn,1:SIZE(Phi),isub)=PH_micro
            !> water potential criterion; check B*K_ave (soil) with Lr*r0 (see MS)
            IF (PH_micro(n-1)-PH_micro(1) .EQ. 0) THEN
               k_ave(irecn,isub) = 0
               ja=1
               GOTO 203
            ELSE
               k_ave(irecn,isub) = (Phi(n-1)-Phi(1))/(PH_micro(n-1)-PH_micro(1))
            ENDIF
            B(irecn,isub)=2*(1-rho**2)/(-2*rho**2*(LOG(rho)-1/2.)-1)
            IF (k_ave(irecn,isub)*B(irecn,isub) .LT. Lr(irecn)*r0) THEN
               countja=.TRUE.
            ENDIF
            !> set switchcriterion; if B*Ks > Lr*r0
            !> set switchcriterion; if B*Ks < Lr*r0
            !> set switchcriterion; if B*Ks +_ 5% < > Lr*r0
            IF ((ave) .AND. switchcriterion.EQ.1) PH_micro(1)=meanPHs+ &
                 r0/k_ave(irecn,isub)*q_out*rho*LOG(1/rho)+r0/ &
                 (k_ave(irecn,isub)*B(irecn,isub))*q_out*rho
            IF (k_ave(irecn,isub)*B(irecn,isub) .LT. Lr(irecn)*r0+Lr(irecn)*r0*5/100 .AND. k_ave(irecn,isub)*B(irecn,isub) .GT. &
                 Lr(irecn)*r0-Lr(irecn)*r0*5/100) THEN
               switchCriterion=2
               PH_micro(1)=meanPHs/2+PHr(irecn,ipl)/2+1/(Lr(irecn)*2)*q_out*rho*LOG(1/rho)+ &
                    1/(Lr(irecn)*2)*q_out*rho
            ELSEIF (k_ave(irecn,isub)*B(irecn,isub) .LT. Lr(irecn)*r0-Lr(irecn)*r0*5/100) THEN
               switchCriterion=3
               PH_micro(1)=PHr(irecn,ipl)+B(irecn,isub)/Lr(irecn)*q_out*rho*LOG(1/rho)+1/Lr(irecn)*q_out*rho
            ENDIF
            IF (switchCriterion .NE. 1 .AND. (.NOT.(tcheck)) ) THEN
               tcheck2=.TRUE.
               tcheck=.TRUE.
            ENDIF
            IF (PH_micro(1).GT.-0.8) THEN
               ja=1
               GOTO 203
            ENDIF
 !>if one big root approach is taken:
 !>
 !> calculate the radius r at which the real centerpoint is located. Do not take the z-height
 !> into account. Find at which position in the radius array this value matches (sort of)
 !> take that index value in PH_micro and assign that PH as interface PH of that center point of a root segment.
200         IF (ave) THEN
               pos_cent=SQRT((cent(irecn,1,isub,ipl)-cp_mean(irecn,1,isub,ipl))*(cent(irecn,1,isub,ipl)-cp_mean(irecn,1,isub,ipl))&
                    + (cent(irecn,2,isub,ipl)-cp_mean(irecn,2,isub,ipl))*(cent(irecn,2,isub,ipl)-cp_mean(irecn,2,isub,ipl))) 
               !> calculate radius from original cp to center of a voxel in 2D plane
               DO i=1,SIZE(r)-1
                  IF (pos_cent .EQ. 0_dp) THEN
                     pos=1
                     EXIT
                  ELSEIF (pos_cent .LE. r(i+1) .AND. pos_cent .GE. r(i)) THEN
                     pos=i+1
                     EXIT
                  ENDIF
               ENDDO
               IF (checker1.EQ.1) THEN
                  PH_root_sub(irecn,isub) = PH_micro2(irecn,pos,isub)
               ELSE
                  !> denote PH_int to all root nodes within this voxel independent of distance
                  PH_root_sub(irecn,isub) = PH_micro2(irecnTemp,pos,isubTemp)/pond_subTemp*w_sub(irecn,isub,ipl) 
                  !> take difference in w_sub between cps into account
                  IF (segsur(irecnTemp) .EQ. 0 .OR. seglen(irecnTemp) .EQ. 0) THEN
                     r0 = 5e-2 !> minimal root radius (in cm)
                  ELSE
                     r0 = 1/(2*pi)*segsur(irecnTemp)/seglen(irecnTemp)
                  END IF
                  IF (k_ave(irecnTemp,isubTemp)*B(irecnTemp,isubTemp) .LT. Lr(irecnTemp)*r0) THEN
                     countja=.TRUE.  !> to check how many root nodes are under "limiting conditions"
                  ENDIF
               ENDIF
            ENDIF
            IF (.NOT.(ave)) THEN
               PH_root_sub(irecn,isub) = PH_micro(1)
            ENDIF
!-----------------------------------------------------------------------
!> Give every soil-root interface node its correct PH value
!> -----------------------------------------------------------------------
            PH_root_sub(irecn,isub) = PH_root_sub(irecn,isub) + zCent*w_sub(irecn,isub,ipl)!calc PH at interface with gravity
         ELSE!checker 
            !> some segments are divided in subsegments in which a cp is located at the edge of the plane (and/or cube),use average approach of original R-SWMS
            sumW=sum_dis(irecn,isub,ipl)
            psub=w_sub(irecn,isub,ipl)
            PH_root_sub(irecn,isub) = ((hnew(corner(1))+zgrid(corner(1)))*pond_i(1)+&
                 (hnew(corner(2))+zgrid(corner(2)))*pond_i(2)+(hnew(corner(3))+&
                 zgrid(corner(3)))*pond_i(3)+(hnew(corner(4))+zgrid(corner(4)))*pond_i(4)+&
                 (hnew(corner(5))+zgrid(corner(5)))*pond_i(5)+(hnew(corner(6))+&
                 zgrid(corner(6)))*pond_i(6)+(hnew(corner(7))+zgrid(corner(7)))*pond_i(7)+&
                 (hnew(corner(8))+zgrid(corner(8)))*pond_i(8))/sumW*psub
         END IF
!---------------------------------------------------------------------------
!> ph at interface
         PHs(irecn,ipl)= PHs(irecn,ipl) + PH_root_sub(irecn,isub)

!> add salinity potential 
         if (lSalinity) then
            PHs(irecn,ipl)=PHs(irecn,ipl) + sigma*PHs_osmotic(cube_i(irecn,isub,ipl))*psub
            PHo(irecn,ipl)=PHo(irecn,ipl) + sigma*PHs_osmotic(cube_i(irecn,isub,ipl))*psub
         endif

203      IF (ja.EQ.1) THEN !> average method used in some cases
            IF (irecn .EQ. 0) THEN
               r0=5e-2
            ELSE
               r0 = 1/(2*pi)*segsur(irecn)/seglen(irecn)
            ENDIF
            temp=Lr(irecn)*r0
            DO j=1,8
               h_matK(j)=hnew(corner(j))
               CALL SetMat_singlevalue(h_matK(j),MatNum(corner(j)),theta_k,K_k,C_k,corner(j))
               temp2(j)=K_k
            ENDDO
            IF (SUM(temp2)/8 .LT. temp) THEN
               countja=.TRUE.
            ENDIF
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

            if(lSalinity) then
               PHs(irecn,ipl)=PHs(irecn,ipl) + sigma*PHs_osmotic(cube_i(irecn,isub,ipl))*psub
               PHo(irecn,ipl)=PHo(irecn,ipl) + sigma*PHs_osmotic(cube_i(irecn,isub,ipl))*psub
            endif

         ENDIF
      END DO !> end do-loop over subsegments
      !> calculate matrix Q (doussan matrix)
      IF (countja) count_nodes=count_nodes+1
      countja=.FALSE.
   END DO
   
   Qd(0:nrec,ipl) =Qi(:,ipl)*PHs(:,ipl)+Q_bc1(:,ipl)
   Qd(nrec+1:,ipl)=Qi(:,ipl)*PHs(:,ipl)+Q_bc2(:,ipl)
 END SUBROUTINE analytical_approach
!=====================================================================================
 SUBROUTINE flux_node(Flow_corner,corner,q_tot_c)
   USE typedef
   USE GridData, ONLY: xgrid,ygrid,zgrid,dxgrid,dygrid,dzgrid,nPt
   USE SolData, ONLY: MatNum, hnew
   IMPLICIT NONE
   
   INTEGER(ap) :: i,j,Flow_corner(1:8,1:6)
   REAL(dp) :: q_tot_c_t
   REAL(sp) :: theta_single,K_single,C_single,Kn,Kcn,Kc,dz
   INTEGER(ap) ::corner(1:8)
   REAL(dp),INTENT(out) :: q_tot_c(1:8)

   q_tot_c(1:8) = 0
   !> loop over all soil nodes from cubic (sub cubic -> subsegment)
   DO i=1,8
      !> calculate soil hydraulic conductivity in the corner node and the surrounding node (in next loop)
      CALL SetMat_singlevalue(hnew(corner(i)),MatNum(corner(i)),theta_single,K_single,C_single, corner(i))
      Kc = K_single
      !> loop over all surrounding nodes
      DO j=1,6
         IF (Flow_corner(i,j).EQ.0) THEN
            q_tot_c_t = 0
            GOTO 40
         ENDIF
         CALL SetMat_singlevalue(hnew(Flow_corner(i,j)),MatNum(corner(j)),theta_single,K_single,C_single,Flow_corner(i,j))
         Kn = K_single
         Kcn = (Kc+Kn)/2._dp !> arithmetic mean soil hydraulic conductivity
         !> distance between the two nodes evaluated
         IF (xgrid(corner(i))-xgrid(Flow_corner(i,j)) .NE. 0) THEN
            dz=dxgrid
         ELSE IF (ygrid(corner(i))-ygrid(Flow_corner(i,j)) .NE. 0) THEN
            dz=dygrid
         ELSE IF (zgrid(corner(i))-zgrid(Flow_corner(i,j)) .NE. 0) THEN
            dz=dzgrid
            !> -> q_tot_c_t = 0 (no flow from z-direction taken into account
         END IF
         !> flow from corner node to surrounding node
         q_tot_c_t = - Kcn * (hnew(corner(i)) - hnew(Flow_corner(i,j)))/dz
         !> take the sum over all the calculates flows from and to the corner node
         !> q_tot_c has size 8 (corner nodes of cube)
40       CONTINUE
         q_tot_c(i) = q_tot_c(i) + q_tot_c_t
      END DO
   END DO
   RETURN
 END SUBROUTINE flux_node
!=====================================================================================
 SUBROUTINE average_method(ipl)
   USE typedef
   USE Paramdata,ONLY: MaxNod,pi,lSalinity
   USE doussanmat, ONLY: PHs,nsub,w_sub,sum_dis,PH_root_sub,Qi,Qd,Q_bc1,Q_bc2,w_dis,&
        loc_Q,Lr,count_nodes,Q_bc,cube_i,PHs_osmotic,PHo
   USE GridData, ONLY: zgrid,nPt
   USE RootData, ONLY: nrec,segsur,seglen,sigma
   USE SolData, ONLY: MatNum, hnew
   IMPLICIT NONE

   INTEGER(ap) :: isub,i,j  
   REAL (dp) ::sumW,psub
   REAL(dp),DIMENSION(1:8) ::pond_i
   REAL(sp),DIMENSION(1:8) ::h_mat,temp2
   REAL(sp) :: K,C,theta,r0,temp
   INTEGER(ap) ::corner(1:8)
   INTEGER(ap) :: irecn,ipl,M
   LOGICAL :: countja
   
   countja=.FALSE.
   count_nodes=0
   DO irecn=1,nrec
      IF (irecn .EQ. 0) THEN
         r0=5e-2
      ELSE
         r0 = 1/(2*pi)*segsur(irecn)/seglen(irecn)
      ENDIF
      temp=Lr(irecn)*r0
      !> initialize soil PH matrix
      PHs(irecn,ipl)=0._dp
      PHo(irecn,ipl)=0._dp
      DO isub=1,nsub(irecn,ipl)
         pond_i(:)=w_dis(irecn,:,isub,ipl)
         corner(:)=loc_Q(irecn,:,isub,ipl)
         DO j=1,8
            h_mat(j)=hnew(corner(j))
            M=MatNum(corner(j))		
            CALL SetMat_singlevalue(h_mat(j),M,theta,K,C,corner(j))
            temp2(j)=K
         ENDDO
         IF (SUM(temp2)/8 .LT. temp) THEN
            countja=.TRUE.
         ENDIF
         sumW=sum_dis(irecn,isub,ipl) !> total inverse distance to node 
         psub=w_sub(irecn,isub,ipl)   !> subsegment length/total segment length ratio
         PH_root_sub(irecn,isub)=((hnew(corner(1))+zgrid(corner(1)))*pond_i(1)+&
              (hnew(corner(2))+zgrid(corner(2)))*pond_i(2)+(hnew(corner(3))+&
              zgrid(corner(3)))*pond_i(3)+(hnew(corner(4))+zgrid(corner(4)))*pond_i(4)+&
              (hnew(corner(5))+zgrid(corner(5)))*pond_i(5)+(hnew(corner(6))+&
              zgrid(corner(6)))*pond_i(6)+(hnew(corner(7))+zgrid(corner(7)))*pond_i(7)+&
              (hnew(corner(8))+zgrid(corner(8)))*pond_i(8))/sumW*psub
         PHs(irecn,ipl)=PHs(irecn,ipl)+PH_root_sub(irecn,isub)    
         
         if(lSalinity) then
            PHs(irecn,ipl)=PHs(irecn,ipl) + sigma*PHs_osmotic(cube_i(irecn,isub,ipl))*psub
            PHo(irecn,ipl)=PHo(irecn,ipl) + sigma*PHs_osmotic(cube_i(irecn,isub,ipl))*psub
         endif
         
         
      END DO    !> end do-loop over subsegments
      !> calculate matrix Qd

      IF (countja) count_nodes=count_nodes + 1
      countja=.FALSE.
   END DO !> end do-loop over irecn
   
   Qd(0:nrec,ipl) = Qi(:,ipl)*PHs(:,ipl)+Q_bc1(:,ipl)
   Qd(nrec+1:,ipl) = Qi(:,ipl)*PHs(:,ipl)+Q_bc2(:,ipl)
   
 END SUBROUTINE average_method
!===================================================================================
!***********************************************************************************
SUBROUTINE logspace(d1,d2,n,y)
  USE typedef
  IMPLICIT NONE

  INTEGER(sp),INTENT(in) :: n
  INTEGER(sp) :: i
  REAL(sp) :: qq(1:n),Y_(1:n)
  REAL(sp) :: d1, d2, d31, d32
  REAL(sp), INTENT(out) :: y(1:n)

  !>logspace algorithm
  d31 = LOG10(ABS(d1))
  d32 = LOG10(ABS(d2))

  DO i=1,n-1
     qq(i) =d31+(i-1)*(d32-d31)/(n-1)
  END DO

  !>put in one array
  y_(1:n-1) = qq(1:n-1)
  y_(n) = d32

  !> convert logscale to normal values
  y = (10)**y_(1:n)
  RETURN
END SUBROUTINE logspace
!***********************************************************************
SUBROUTINE SetMat_ana(hTemp,theta,K,C)
  USE GridData
  USE SolData, ONLY: par
  USE WatFun
  USE Doussanmat, ONLY : n
  USE RhizoData, ONLY : lRhizo
  IMPLICIT NONE

  REAL(sp), INTENT(in) :: hTemp(1:n)
  REAL(sp):: Ti(1:n)!,t0,t1
  REAL(sp), INTENT(out) :: theta(1:n),C(1:n),K(1:n)
  INTEGER(ap):: i,M
  IF (lRhizo) THEN
      STOP 'Someone callSetMat_ana while RhizoModel'
  ENDIF
  DO i=1,n
     M=nMat
     !> nodal conductivity values:
     K(i)=FKP(hTemp(i),par(:,M),1)*Bxy(i)
     C(i)=FCP(hTemp(i),par(:,M))*Dxy(i)/Axy(i)
     Ti(i)=Fth(hTemp(i),par(:,M),1)
     theta(i)=par(2,M)*Exy(i)+(Ti(i)-par(2,M))*Dxy(i) 
  END DO

  RETURN
END SUBROUTINE SetMat_ana
!***********************************************************************
SUBROUTINE SetMat_anaLibr(hTemp,theta,K,C)
  USE GridData
  USE SolData, ONLY: par
  USE WatFun
  USE Doussanmat, ONLY: nLibr
  USE RhizoData, ONLY : lRhizo
  IMPLICIT NONE

  REAL(sp), INTENT(in) :: hTemp(1:nLibr)
  REAL(sp), INTENT(out) :: theta(1:nLibr),C(1:nLibr),K(1:nLibr)
  INTEGER(ap):: i,M
  IF (lRhizo) THEN
      STOP 'Someone call SetMat_anaLibr while RhizoModel'
  ENDIF
  M=nMat
  DO i=1,nLibr
     !> nodal conductivity values:
     K(i)=FKP(hTemp(i),par(:,M),1)
     C(i)=FCP(hTemp(i),par(:,M))
     theta(i)=Fth(hTemp(i),par(:,M),i)
  END DO

  RETURN
END SUBROUTINE SetMat_anaLibr
!***********************************************************************
SUBROUTINE SetMat_singlevalue(hTemp,M,theta,K,C,soilNodeAddress)
  USE GridData
  USE SolData, ONLY: par,soiltab
  USE WatFun
  IMPLICIT NONE

  REAL(sp), INTENT(in) :: hTemp
  INTEGER(sp), INTENT(in) :: M
  INTEGER :: soilNodeAddress
  REAL(sp), INTENT(out) :: theta,C,K
  !> nodal conductivity values:
  IF (soiltab) THEN
     K=FKP_soiltab(hTemp,M)
     C=FCP_soiltab(hTemp,M)
     theta=Fth_soiltab(hTemp,M)
  ELSE
     K=FKP(hTemp,par(:,M), soilNodeAddress)
     C=FCP(hTemp,par(:,M))
     theta=Fth(hTemp,par(:,M), soilNodeAddress) 
  END IF
  RETURN
END SUBROUTINE SetMat_singlevalue
!**************************************************************************************
! weigth for each node in fucntion of their neighbouors
SUBROUTINE CalcWnodes
  USE typedef
  USE GridData
  IMPLICIT NONE

  INTEGER(ap) :: iE

  !> calculate actual overall transpiration rate from soil domain:
  Wn=0.0_dp
  DO iE=1,nElm
     !> assign to Wn the number of times a cuboid corner node is called
     Wn(elmnod(1,iE))=Wn(elmnod(1,iE))+1
     Wn(elmnod(2,iE))=Wn(elmnod(2,iE))+1
     Wn(elmnod(3,iE))=Wn(elmnod(3,iE))+1
     Wn(elmnod(4,iE))=Wn(elmnod(4,iE))+1
     Wn(elmnod(5,iE))=Wn(elmnod(5,iE))+1
     Wn(elmnod(6,iE))=Wn(elmnod(6,iE))+1
     Wn(elmnod(7,iE))=Wn(elmnod(7,iE))+1
     Wn(elmnod(8,iE))=Wn(elmnod(8,iE))+1
  END DO
  Wn=Wn/8._sp

  RETURN
END SUBROUTINE CalcWnodes
!********************************************************************************
!> ### for continuous domains determine length of the root that is growing back into the
!> other side of the soil domain ###
SUBROUTINE roottrans(xt,yt,inode,isub,ipl)
  USE typedef
  USE GridData, ONLY: nex,ney,dxgrid,dygrid
  USE DoussanMat, ONLY: transroot
  USE DomData
  IMPLICIT NONE

  REAL(sp), INTENT(in) :: xt,yt
  REAL(sp) :: xd,yd
  INTEGER(ap), INTENT(in) :: inode,ipl,isub
  
  xd = xt
  yd = yt
  transroot(inode,1:2,isub,ipl)=0
  DO WHILE (xd.GE.xmax) 
     xd=xd-nex*dxgrid
     transroot(inode,1,1,ipl)=transroot(inode,1,isub,ipl)-1
  END DO
  DO WHILE (xd.LT.xmin)
     xd=xd+nex*dxgrid
     transroot(inode,1,1,ipl)=transroot(inode,1,isub,ipl)+1
  END DO
  DO WHILE (yd.GE.ymax)
     yd=yd-ney*dygrid
     transroot(inode,2,1,ipl)=transroot(inode,2,isub,ipl)-1
  END DO
  DO WHILE (yd.LT.ymin)
     yd=yd+ney*dygrid
     transroot(inode,2,1,ipl)=transroot(inode,2,isub,ipl)+1
  END DO

END SUBROUTINE roottrans
!******************************************************************************** 
!> ### for continuous domains determine length of the root that is growing back into the
!> other side of the soil domain ###
SUBROUTINE tiptrans(xt,yt,inode,isub,ipl)
  USE typedef
  USE GridData, ONLY: nex,ney,dxgrid,dygrid
  USE DoussanMat, ONLY: transtip
  USE DomData
  IMPLICIT NONE

  REAL(sp), INTENT(in) :: xt,yt
  REAL(sp) :: xd,yd
  INTEGER(ap), INTENT(in) :: inode,ipl,isub

  xd = xt
  yd = yt
  transtip(inode,1:2,isub,ipl)=0
  DO WHILE (xd.GE.xmax)
     xd=xd-nex*dxgrid
     transtip(inode,1,1,ipl)=transtip(inode,1,isub,ipl)-1
  END DO
  DO WHILE (xd.LT.xmin)
     xd=xd+nex*dxgrid
     transtip(inode,1,1,ipl)=transtip(inode,1,isub,ipl)+1
  END DO
  DO WHILE (yd.GE.ymax)
     yd=yd-ney*dygrid
     transtip(inode,2,1,ipl)=transtip(inode,2,isub,ipl)-1
  END DO
  DO WHILE (yd.LT.ymin)
     yd=yd+ney*dygrid
     transtip(inode,2,1,ipl)=transtip(inode,2,isub,ipl)+1
  END DO

END SUBROUTINE tiptrans
!********************************************************************************
