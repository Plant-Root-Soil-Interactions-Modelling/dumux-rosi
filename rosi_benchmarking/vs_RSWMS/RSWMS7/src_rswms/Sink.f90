! ==============================================================================
! Source file SINK |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
SUBROUTINE BetDis(t)
  USE typedef
  USE RootData
  USE PlntData
  USE GridData
  IMPLICIT NONE
  REAL(sp):: betcon(8)
  INTEGER(ap):: corner(8)
  LOGICAL intsec,split
  REAL(sp):: x1,x2,xA,xB,y1,y2,yA,yB,z1,z2,zA,zB
  REAL(sp):: xInt,yInt,zInt,xCent,yCent,zCent
  REAL (sp)::dis,sqrdis,sumB,weight,segage,srface
  REAL(sp):: totlen,prvsur,t 
  INTEGER (sp)::i,ibr,iprv,igrow,iseg,iorder,iurf,iface,ic,irec,imin,ipl,isub,typ
  REAL(dp):: blengt

  PrvSur=0.0_dp
  TotSur=0.0_dp
  TotLen=0.0_dp
  split=.FALSE.
  ipl=1
  DO iseg=1,nrec
     TotLen=TotLen+seglen(iseg)
     PrvSur=PrvSur+segsur(iseg)
  END DO
  IF (TotLen.LT.1.E-20_dp) TotLen=1.E-20_dp
  IF (PrvSur.LT.1.E-20_dp) PrvSur=1.E-20_dp
  DO i=1,nPt
     betaw(i)=0.0_dp
     betac(i)=0.0_dp
  END DO
  !> go through each root segment and update the node surface function:
2 DO  ibr=1,nbr
     !> find the tip segment of the branch 'ibrnch'
     irec=nrec+1
11   irec=irec-1
     IF (ibrseg(irec).NE.ibr) GOTO 11
     IF (seglen(irec).LT.1.E-20_dp) THEN
        iprv=irecpr(irec)
        IF (iprv.EQ.0) GOTO 2
        GOTO 101
     ELSE
        DO  igrow=1,ngrow
           IF (ibrgrw(igrow).EQ.ibr) THEN
              iprv=irec
              IF (iprv.EQ.0) GOTO 2
              xA=xg(igrow)
              yA=yg(igrow)
              zA=zg(igrow)
              GOTO 102
           ENDIF
        END DO
     ENDIF
100  iprv=irecpr(irec)
     IF (iprv.EQ.0) GOTO 2
101  IF (ibrseg(iprv).EQ.ibrseg(irec)) THEN
        xA=xs(irec)+xplant(ipl)
        yA=ys(irec)+yplant(ipl)
        zA=zs(irec)
     ELSE
        GOTO 2
     ENDIF
     !> find cuboid around segment:
102  CALL Neighb(xA,yA,zA,corner,imin)
     x1=xgrid(corner(1))
     x2=xgrid(corner(4))
     y1=ygrid(corner(1))
     y2=ygrid(corner(4))
     z1=zgrid(corner(1))
     z2=zgrid(corner(5))
     !> find the other end of the segment:
     xB=(xs(iprv))+xplant(ipl)
     yB=(ys(iprv))+yplant(ipl)
     zB=(zs(iprv))
     !> calculate segment surface:
     blengt=seglen(iprv)
     srface=segsur(iprv)
     TotSur=TotSur+segsur(iprv)
     !> calculate segment weighing factor according to age:
     iorder=0
     iorder=ordseg(iprv)
     IF (maizeroottyp) THEN
        IF (iorder.LT.12) THEN  		!> In RootTyp, maize root types below 12 are principal roots (12 and 13 are lateral)
           typ=1
        ELSE
           typ=2
        ENDIF
     ELSEIF (loliumroottyp) THEN
        IF (iorder.LT.3) THEN  		!> In RootTyp, lolium root types below 3 are principal roots
           typ=1
        ELSEIF (iorder.LT.5) THEN
           typ=2
        ELSE
           typ=3
        ENDIF
     ELSEIF (wheatroottyp) THEN
        IF (iorder.LT.19) THEN  		!> In RootTyp, wheat root types below 19 are principal roots
           typ=1
        ELSEIF (iorder.LT.20) THEN
           typ=2
        ELSE
           typ=3
        ENDIF
     ELSE
        IF (iorder.EQ.0) THEN  !> type 0 is the seed+small segment
           typ=1
        ELSEIF (iorder.GT.3) THEN !> no more than 3 root types
           typ=3
        ELSE
           typ=iorder
        ENDIF
     ENDIF
     segage=t-timorg(iprv)
     IF (segage.GE.0.0_dp.AND.lUrf) THEN
        IF (segage.GE.age(typ,nUrf(typ))) THEN
           Weight=Urf(typ,nUrf(typ))
        ELSE
           iUrf=nUrf(typ)
4          iUrf=iUrf-1
           IF ((segage.LT.age(typ,iUrf)).AND.(iUrf.GT.1)) GOTO 4
           Weight=Urf(typ,iUrf)+(segage-age(typ,iUrf))/&
                (age(typ,iUrf+1)-age(typ,iUrf))*&
                (Urf(typ,iUrf+1)-Urf(typ,iUrf))
        ENDIF
     ELSE
        Weight=1.0_dp
     ENDIF
     IF (intsec(xA,yA,zA,xB,yB,zB,irec,iprv,isub,ipl,x1,y1,z1,x2,y2,z2,xInt,yInt,zInt,iFace)) THEN
        split=.TRUE.
        ! calculate subsegment length and surface...
41      blengt=SQRT((xInt-xA)*(xInt-xA)+(yInt-yA)*(yInt-yA)+(zInt-zA)*(zInt-zA))
        srface=blengt*segsur(iprv)/seglen(iprv)
        ! calculate subsegment center coordinates...
        xCent=xA+(xInt-xA)/2._dp
        yCent=yA+(yInt-yA)/2._dp
        zCent=zA+(zInt-zA)/2._dp
        ! and calculate the distribution for each node:
        sumB=0.0_sp
        DO ic=1,8
           sqrdis=(xCent-xgrid(corner(ic)))*(xCent-xgrid(corner(ic)))+ &
                (yCent-ygrid(corner(ic)))*(yCent-ygrid(corner(ic)))+ &
                (zCent-zgrid(corner(ic)))*(zCent-zgrid(corner(ic)))
           dis=SQRT(sqrdis)
           IF (dis.LT.1.E-20_dp) dis=1.E-20_sp
           betcon(ic)=1._sp/dis
           sumB=sumB+betcon(ic)
        END DO
        DO ic=1,8
           betaw(corner(ic))=betaw(corner(ic))+betcon(ic)/sumB*Weight*(blengt/TotLen)
           betac(corner(ic))=betac(corner(ic))+betcon(ic)/sumB*Weight*(srface/PrvSur)
        END DO
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
        !calculate cuboid´s corners coordinates:
        x1=xgrid(corner(1))
        x2=xgrid(corner(4))
        y1=ygrid(corner(1))
        y2=ygrid(corner(4))
        z1=zgrid(corner(1))
        z2=zgrid(corner(5))
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
     DO ic=1,8
        sqrdis=(xCent-xgrid(corner(ic)))*(xCent-xgrid(corner(ic)))+&
             (yCent-ygrid(corner(ic)))*(yCent-ygrid(corner(ic)))+&
             (zCent-zgrid(corner(ic)))*(zCent-zgrid(corner(ic)))
        dis=SQRT(sqrdis)
        IF (dis.LT.1.E-20_dp) dis=1.E-20_dp
        betcon(ic)=1._dp/dis
        sumB=sumB+betcon(ic)
     END DO
     DO ic=1,8
        betaw(corner(ic))=betaw(corner(ic))+betcon(ic)/sumB*Weight*(blengt/TotLen)
        betac(corner(ic))=betac(corner(ic))+betcon(ic)/sumB*Weight*(srface/PrvSur)
     END DO
     irec=iprv
     GOTO 100
  END DO
  RETURN
END SUBROUTINE BetDis
!********************************************************************************************
SUBROUTINE BetNrm
  USE typedef
  USE PlntDATA
  USE RootData, ONLY: OmegaC,lJarvis,lUrf
  USE DoussanMat, ONLY: stresfun
  USE GridData
  IMPLICIT NONE
  INTEGER(ap):: i,j,k,l,iE,ise,kOut,ipl
  REAL(dp):: betwe,betce,VEl,Sbetac,Sbetaw
  CHARACTER :: file*17
  kout=0
  ipl=1
  Sbetaw=0.0_dp
  Sbetac=0.0_dp
  
  DO iE=1,nElm
     DO iSE=1,5
        
        i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
        j=elmnod(iL(2,iSE,subN(iE)),iE)
        k=elmnod(iL(3,iSE,subN(iE)),iE)
        l=elmnod(iL(4,iSE,subN(iE)),iE)
        
        betwE=(betaw(i)+betaw(j)+betaw(k)+betaw(l))/4.
        Sbetaw=Sbetaw+betwE
        betcE=(betac(i)+betac(j)+betac(k)+betac(l))/4.
        Sbetac=Sbetac+betcE
     END DO
  END DO
  VEl=dxGrid*dyGrid*dzGrid/5._dp
  IF (Sbetaw.LT.1.E-20_dp) RETURN
  Sbetaw=Sbetaw*VEl
  DO  i=1,nPt
     betaw(i)=betaw(i)/Sbetaw
  END DO
  
  WRITE (file,'(A10)')'out/Feddes'
  WRITE (file(11:11),'(I1)') ipl
  WRITE (file(12:12),'(A1)') '.'
  IF (kout.LT.10) THEN
     WRITE (file(13:13),'(I1)') kout
  ELSEIF (kout.LT.100) THEN
     WRITE (file(13:14),'(I2)') kout
  ELSE
     WRITE (file(13:15),'(I3)') kout
  ENDIF
  OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
  WRITE (8,'(''********** FEDDES ROOT WATER UPTAKE MODEL PARMS **********'')')
  WRITE (8,*)
  WRITE (8,'(''Stress function ? (0=no,1=Feddes,2=van Genuchten)'')')
  WRITE (8,'(I1)') stresfun
  WRITE (8,*)
  WRITE (8,'(''Feddes stress function:'')')
  WRITE (8,'(''h0      h1      h2      h3	 [hPa]'')')
  WRITE (8,'(F5.1,2X,F5.1,2X,F8.1,2X,F9.1)') h0,h1,h2,h3
  WRITE (8,*)
  WRITE (8,'(''van Genuchten stress function:'')')
  WRITE (8,'(''p50[hPa]   h50[hPa]   p1[-]  p2[-]'')')
  WRITE (8,'(F8.1,2X,F8.1,2X,F5.1,2X,F5.1)') p50,h50,p1,p2
  WRITE (8,*)
  WRITE (8,'(''Use of Jarvis (1989) function for compensatory RWU prediction ? (1=yes,2=no)'')')
  IF (lJarvis) THEN
     WRITE (8,'(I1)') 1
  ELSE
     WRITE (8,'(I1)') 2
  ENDIF
  WRITE (8,*)
  WRITE (8,'(''Omega_c value [-] (1=uncompensated RWU,0=fully compensated RWU)'')')
  WRITE (8,'(F4.2)') OmegaC
  WRITE (8,*)
  WRITE (8,'(''Uptake Reduction Function (URF) [-] with segment age? Only if rRLD not given beside (1=yes,2=no)'')')
  IF (lUrf) THEN
     WRITE (8,'(I1)') 1
  ELSE
     WRITE (8,'(I1)') 2
  ENDIF
  WRITE (8,*)
  WRITE (8,'(''* Root Length distribution of the root system'')')
  WRITE (8,'(''Number of elements in each dimension [-] and elements size [cm]'')')
  WRITE (8,'(I5,1X,I5,1X,I5,3(1X,F10.3))') nex,ney,nez,dxGrid,dyGrid,dzGrid
  WRITE (8,'(''    Node#    betaw(RLD)'')')
  DO i=1,nPt
     WRITE (8,'(I6,3X,F12.10)') i,betaw(nPt+1-i)		!nodes order inversed for Hydrus
  ENDDO
  CLOSE (8)
  RETURN
END SUBROUTINE BetNrm
!********************************************************************************************
!> ### root length density distribution ###
SUBROUTINE RLDdis(t,kout)
  USE typedef
  USE ParamData, ONLY: pi
  USE RootData
  USE PlntData
  USE GridData, ONLY: RLD,RSD,continu,dxGrid,dyGrid,dzGrid,VElm,nElm,nex,ney,nez,xGrid,yGrid
  USE DoussanMat, ONLY: transroot,nsub,Intc,cube_i,stresfun
  IMPLICIT NONE
  INTEGER(ap):: corner(8)
  LOGICAL intsec
  REAL(sp):: x1,x2,xA,xB,y1,y2,yA,yB,z1,z2,zA,zB
  REAL(sp):: xInt,yInt,zInt
  REAL (sp):: Weight,segage
  REAL(sp):: t
  REAL(dp):: blengt,srface,TotLen,iz
  INTEGER (sp):: i,k,ifol,ibr,igrow,iseg,iorder,iurf,iface,ic,irec,imin,ipl,isub,typ,transrootA(2)
  INTEGER(ap), INTENT(in)::kout
  CHARACTER :: file*17
  
  WRITE (*,*) '... Calculating RLD ...'
  TotSur=0.0_dp
  TotLen=0.0_dp
  ipl=1
  ALLOCATE (RLD(1:nElm))
  RLD(1:nElm)=0.0_dp
  ALLOCATE (RSD(1:nElm))
  RSD(1:nElm)=0.0_dp
  !> go through each root segment and update the node surface function:
2 DO ibr=1,nbr						!> Loop on root branches
     !> find the tip segment of the branch 'ibrnch'
     irec=nrec+1
11   irec=irec-1
     IF (ibrseg(irec).NE.ibr) GOTO 11
     DO igrow=1,ngrow
        IF (ibrgrw(igrow).EQ.ibr) THEN
           ifol=nrec+igrow
           IF (irec.EQ.0) GOTO 2			!"if root sys top, next branch!"
           xA=xg(igrow)+xplant(ipl)
           yA=yg(igrow)+yplant(ipl)
           zA=zg(igrow)
           IF (continu) THEN
              transrootA(1:2)=transroot(nrec+igrow,1:2,nsub(nrec+igrow,ipl),ipl)
              xA=xA+transrootA(1)*(nex*dxgrid)
              yA=yA+transrootA(2)*(ney*dygrid)
           ENDIF
           GOTO 102
        ENDIF
     END DO
100  irec=irecpr(ifol)					!Loop on root segments
     IF (irec.EQ.0) GOTO 2				!"if root sys top, next branch!"
     IF (ibrseg(irec).EQ.ibrseg(ifol)) THEN		!"on the same brench?"
        xA=Intc(ifol,1,nsub(ifol,ipl),ipl)
        yA=Intc(ifol,2,nsub(ifol,ipl),ipl)
        zA=Intc(ifol,3,nsub(ifol,ipl),ipl)
     ELSE
        GOTO 2						!"if branch top, next branch!"	By doing that, we commit the error of ommitting the length of the first part of all branches
     ENDIF
102  IF (lUrf) THEN
        ! calculate segment weighing factor according to age:
        iorder=0
        iorder=ordseg(irec)
        IF (maizeroottyp) THEN
           IF (iorder.LT.12) THEN  			!In RootTyp, maize root types below 12 are principal roots (12 and 13 are lateral)
              typ=1
           ELSE
              typ=2
           ENDIF
        ELSEIF (loliumroottyp) THEN
           IF (iorder.LT.3) THEN  			!In RootTyp, lolium root types below 3 are principal roots
              typ=1
           ELSEIF (iorder.LT.5) THEN
              typ=2
           ELSE
              typ=3
           ENDIF
        ELSEIF (wheatroottyp) THEN
           IF (iorder.LT.19) THEN  			!In RootTyp, wheat root types below 19 are principal roots
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
        segage=t-timorg(irec)
        IF (segage.GE.0.0_dp) THEN
           IF (segage.GE.age(typ,nUrf(typ))) THEN
              Weight=Urf(typ,nUrf(typ))
           ELSE
              iUrf=nUrf(typ)
4             iUrf=iUrf-1
              IF ((segage.LT.age(typ,iUrf)).AND.(iUrf.GT.1)) GOTO 4
              Weight=Urf(typ,iUrf)+(segage-age(typ,iUrf))/&
                   (age(typ,iUrf+1)-age(typ,iUrf))*&
                   (Urf(typ,iUrf+1)-Urf(typ,iUrf))
           ENDIF
        ELSE
           Weight=1.0_dp
        ENDIF
     ELSE
        Weight=1.0_dp
     ENDIF
     DO isub=1,nsub(irec,ipl)				!Loop on sub-segments
        xB=Intc(irec,1,isub,ipl)
        yB=Intc(irec,2,isub,ipl)
        zB=Intc(irec,3,isub,ipl)
        ! calculate subsegment length and surface...
        IF (continu) THEN
           blengt=SQRT((xB-xA+(transrootA(1)-transroot(irec,1,isub,ipl))*nex*dxGrid)**2+(yB-yA+(transrootA(2)-transroot(irec,2,isub,ipl))*ney*dyGrid)**2+(zB-zA)**2)
        ELSE
           blengt=SQRT((xB-xA)**2+(yB-yA)**2+(zB-zA)**2)
        ENDIF
        i=0
42      IF (seglen(irec).GT.0.) THEN
           srface=blengt*segsur(irec+i)/seglen(irec+i)
        ELSE
           i=i-1						!Temporary solution for segment whose saved length is 0
           GOTO 42
        ENDIF
        RLD(cube_i(irec,isub,ipl))=RLD(cube_i(irec,isub,ipl))+blengt*Weight
        RSD(cube_i(irec,isub,ipl))=RSD(cube_i(irec,isub,ipl))+srface*Weight
        TotLen=TotLen+blengt*Weight
        TotSur=TotSur+srface*Weight
        xA=xB
        yA=yB
        zA=zB
        transrootA(1:2)=transroot(irec,1:2,isub,ipl)
     END DO						!Loop on sub-segments
     ifol=irec
     GOTO 100						!Loop on root segments
  END DO					!Loop on root branches
  RLD=RLD/(dxGrid*dyGrid*dzGrid)	!/TotLen
  RSD=RSD/(dxGrid*dyGrid*dzGrid)	!/TotSur
  
  IF (lFed) THEN
     WRITE (file,'(A10)')'out/Feddes'
     WRITE (file(11:11),'(I1)') ipl
     WRITE (file(12:12),'(A1)') '.'
     IF (kout.LT.10) THEN
        WRITE (file(13:13),'(I1)') kout
     ELSEIF (kout.LT.100) THEN
        WRITE (file(13:14),'(I2)') kout
     ELSE
        WRITE (file(13:15),'(I3)') kout
     ENDIF
     OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
     WRITE (8,'(''********** FEDDES ROOT WATER UPTAKE MODEL PARMS **********'')')
     WRITE (8,*)
     WRITE (8,'(''Stress function ? (0=no,1=Feddes,2=van Genuchten)'')')
     WRITE (8,'(I1)') stresfun
     WRITE (8,*)
     WRITE (8,'(''Feddes stress function:'')')
     WRITE (8,'(''h0      h1      h2      h3	 [hPa]'')')
     WRITE (8,'(F5.1,2X,F5.1,2X,F8.1,2X,F9.1)') h0,h1,h2,h3
     WRITE (8,*)
     WRITE (8,'(''van Genuchten stress function:'')')
     WRITE (8,'(''p50[hPa]   h50[hPa]   p1[-]  p2[-]'')')
     WRITE (8,'(F8.1,2X,F8.1,2X,F5.1,2X,F5.1)') p50,h50,p1,p2
     WRITE (8,*)
     WRITE (8,'(''Use of Jarvis (1989) function for compensatory RWU prediction ? (1=yes,2=no)'')')
     IF (lJarvis) THEN
        WRITE (8,'(I1)') 1
     ELSE
        WRITE (8,'(I1)') 2
     ENDIF
     WRITE (8,*)
     WRITE (8,'(''Omega_c value [-] (1=uncompensated RWU,0=fully compensated RWU)'')')
     WRITE (8,'(F4.2)') OmegaC
     WRITE (8,*)
     WRITE (8,'(''Uptake Reduction Function (URF) [-] with segment age? Only if rRLD not given beside (1=yes,2=no)'')')
     IF (lUrf) THEN
        WRITE (8,'(I1)') 1
     ELSE
        WRITE (8,'(I1)') 2
     ENDIF
     WRITE (8,*)
     WRITE (8,'(''* Root Length distribution of the root system'')')
     WRITE (8,'(''Number of elements in each dimension [-] and elements size [cm]'')')
     WRITE (8,'(I5,1X,I5,1X,I5,3(1X,F10.3))') nex,ney,nez,dxGrid,dyGrid,dzGrid
     WRITE (8,'(''    Element#    RLD'')')
     DO i=1,nElm
        WRITE (8,'(I6,3X,F12.10)') i,RLD(i)
     ENDDO
     CLOSE (8)

  ELSEIF (ldJvL) THEN
     WRITE (file,'(A8)')'out/dJvL'
     WRITE (file(9:9),'(I1)') ipl
     WRITE (file(10:10),'(A1)') '.'
     IF (kout.LT.10) THEN
        WRITE (file(11:11),'(I1)') kout
     ELSEIF (kout.LT.100) THEN
        WRITE (file(11:12),'(I2)') kout
     ELSE
        WRITE (file(11:13),'(I3)') kout
     ENDIF
     OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
     WRITE (8,'(''********** DE JONG VAN LIER APPROACH COUPLED TO COUVREUR ET AL. (2012) RWU MODEL **********'')')
     WRITE (8,*)
     WRITE (8,'(''Number of past days considered in the effective sink (if set to zero, only the last sink will be considered)'')')
     WRITE (8,'(F5.2)') tlim_dJvL
     WRITE (8,*)
     WRITE (8,'(''* Distribution of the geometrical parameter for depletion zone characterizing'')')
     WRITE (8,'(''Number of elements in each dimension [-] and elements size [cm]'')')
     WRITE (8,'(I5,1X,I5,1X,I5,3(1X,F10.3))') nex,ney,nez,dxGrid,dyGrid,dzGrid
     WRITE (8,'(''    Element#    Rho'')')
     WRITE (*,*) '... Calculating Rho ...'
     ALLOCATE (Rho(1:nElm))
     DO i=1,nElm
        IF (RLD(i).GT.0) THEN
           Rho(i)=4/((RSD(i)/RLD(i)/2/pi)**2-0.53**2/pi/RLD(i)+2*(1/(pi*RLD(i))+(RSD(i)/RLD(i)/2/pi)**2)*LOG(0.53*2*pi*RLD(i)/RSD(i)/SQRT(pi*RLD(i))))
        ELSE
           Rho(i)=999.9999
        ENDIF
        WRITE (8,'(I6,3X,F10.6)') i,Rho(i)
     ENDDO
     CLOSE (8)
  ENDIF
  RETURN
END SUBROUTINE RLDdis
!********************************************************************************************
!> ### Couvreur RWU model input; Couvreur.in ###
SUBROUTINE CouIn
  USE typedef
  USE GridData
  USE RootData, ONLY: Krs,Kcomp,ldJvL,Rho,tlim_dJvL,lno_Archi,lSUF,nrec
  USE DoussanMat, ONLY: hx_min,stresfun
  IMPLICIT NONE
  INTEGER(sp):: el,i_dJvL,i_UptakeDomain
  REAL(sp):: element
  
   OPEN (Unit=10,FILE='in/Couvreur.in',STATUS='OLD',ERR=10)
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) stresfun,i_UptakeDomain, i_dJvL
  READ (10,*)
  READ (10,*)
  IF (i_UptakeDomain.EQ.0) lSUF=.TRUE.
  IF (i_dJvL.EQ.1) ldJvL=.TRUE.
  IF (lno_Archi) THEN
     READ (10,*) Krs, Kcomp, hx_min
     READ (10,*)
     READ (10,*)	
     IF (lSUF) THEN
        READ (10,*) nrec
        READ (10,*)
	READ (10,*)
        ALLOCATE (SSF(1:nrec))
        DO el=1,nrec
           READ (10,*,ERR=20) element,SSF(el)
        ENDDO
     ELSE
	READ (10,*)nElm
        READ (10,*)
        READ (10,*)
        ALLOCATE (SSF(1:nElm))
        DO el=1,nElm
           READ (10,*,ERR=20) element,SSF(el)
        ENDDO
     ENDIF
  ENDIF
  CLOSE (10)
  IF (ldJvL) THEN
     CALL MFPC_Table
     OPEN (Unit=15,FILE='in/dJvL.in',STATUS='OLD',ERR=30)
     READ (15,*)
     READ (15,*)
     READ (15,*)
     READ (15,*,ERR=40) tlim_dJvL
     IF (lno_Archi) THEN
        READ (15,*)
        READ (15,*)
        READ (15,*)
        READ (15,*,ERR=40) nexRho,neyRho,nezRho,dxRho,dyRho,dzRho
        READ (15,*)
        ALLOCATE (Rho(1:nElm))
        DO el=1,nexRho*neyRho*nezRho
           READ (15,*,ERR=40) element,Rho(el)
        ENDDO
     ENDIF
     CLOSE (15)
  ENDIF
  RETURN
10 STOP 'File  < Couvreur.in >  not found -- program terminated.'
20 STOP 'Data inconsistency in  < Couvreur.in >  -- program terminated.'
30 STOP 'File  < dJvL.in >  not found -- program terminated.'
40 STOP 'Data inconsistency in  < dJvL.in >  -- program terminated.'
END SUBROUTINE CouIn
!********************************************************************************************
!> ### Input for Feddes RWU; Feddes.in ###
SUBROUTINE FedIn
  USE typedef
  USE GridData
  USE RootData, ONLY: lUrf,lJarvis,OmegaC,lno_Archi
  USE DoussanMat, ONLY: stresfun
  USE PlntData, ONLY: h0,h1,h2,h3,p50,h50,p1,p2
  IMPLICIT NONE
  INTEGER(sp):: el,element,i_Jarvis,i_Urf

  OPEN (Unit=10,FILE='in/Feddes.in',STATUS='OLD',ERR=10)
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) stresfun
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) h0,h1,h2,h3
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) p50,h50,p1,p2
  READ (10,*)
  READ (10,*)
  READ (10,*) i_Jarvis
  IF (i_Jarvis.EQ.1) lJarvis=.TRUE.
  READ (10,*)
  READ (10,*)
  READ (10,*) OmegaC
  READ (10,*)
  READ (10,*)
  READ (10,*) i_Urf
  IF (i_Urf.EQ.1) lUrf=.TRUE.
  IF (lno_Archi) THEN
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) nexRLD,neyRLD,nezRLD,dxRLD,dyRLD,dzRLD
     READ (10,*)
     ALLOCATE (RLD(1:nexRLD*neyRLD*nezRLD))
     DO el=1,nexRLD*neyRLD*nezRLD
        READ (10,*,ERR=20) element,RLD(el)
     ENDDO
  ENDIF
  CLOSE (10)
  RETURN
10 STOP 'File  < Feddes.in >  not found -- program terminated.'
20 STOP 'Data inconsistency in  < Feddes.in >  -- program terminated.'
END SUBROUTINE FedIn
!**************************************************************************************
!> ### Calculation of the sink term ###
SUBROUTINE SetSnk(t,BCr,BCtp)
  USE TypeDef
  USE Paramdata,ONLY: lChem
  USE TempData, ONLY: tem
  USE PlntData, ONLY: Tpot,Tact,PHcollar,CMm,VMax,fk,xin,TotSur
  USE GridData, ONLY: nPt,nElm,betaw,betac,sink,csink,xgrid,ygrid,zgrid,elmnod,subN,iL,SSF,Vn,VElm,HElm,dxGrid,dyGrid,dzGrid,RLD
  USE DoussanMat, ONLY : sink_cube,betac_cube,betac_cube2,csink_cube,stressBC,Lr,Lr_pot,Khr,Inv_c1KxKr,Inv_ciiKrKr,Joutr,w_sub,nsub,cube_i,nplant,PHs,loc_Q,w_dis,sum_dis
  USE RootData, ONLY: lCou,lSUF,lDou,lFed,lJarvis,OmegaC,Krs,Kcomp,lno_RWU,ldJvL,tlim_dJvL,Rho,lPast,tPast,sinkredDPast,nrec,seglen,segsur,Hseq
  USE SolData, ONLY: hnew,conc,par,theta,MatNum
  USE Watfun, ONLY: Fh_from_Th,Fh_from_mfp_soiltab,Fmfp_soiltab
  IMPLICIT NONE
  REAL(sp), ALLOCATABLE, DIMENSION(:) :: Phi,wPast,Mseq
  REAL(sp) :: alpha,active,active_cube,ConcE,t,HElm2(nElm),Lr_fact(1:nrec),Gap,AQPc
  REAL(dp), INTENT(inout) ::BCr
  INTEGER (sp),INTENT(inout) :: BCtp
  INTEGER(ap) ::i,j,k,l,iE,iSE,isub,irecn,ipl,c_i
  
  IF (lFed) THEN
     Tpot=Bcr
     DO i=1,nElm
        IF (RLD(i).GT.1.E-20_dp) THEN
           sink_cube(i)=RLD(i)/SUM(RLD)*ABS(Tpot)/VElm(i)*alpha(Fh_from_Th(SUM(theta(elmnod(1:8,i)))/8.0_sp,par(:,MatNum(elmnod(1,i)))),SUM(Conc(elmnod(1:8,i)))/8.0_sp,SUM(tem(elmnod(1:8,i)))/8.0_sp)
        ELSE
           sink_cube(i)=0.0_dp
        ENDIF
     END DO
     Tact(1)=DOT_PRODUCT(sink_cube(1:nElm),VElm)
     IF (lJarvis.AND.OmegaC.LT.1) THEN
        sink_cube=sink_cube/MAX(OmegaC,Tact(1)/ABS(Tpot))
        Tact(1)=DOT_PRODUCT(sink_cube(1:nElm),VElm)
     ENDIF

  ELSEIF (lCou.AND..NOT.lSUF) THEN
     ALLOCATE (Phi(1:nElm))
     IF (ldJvL) THEN
        ALLOCATE (wPast(1:lPast))
        ALLOCATE (Mseq(1:nElm))
        IF (lPast.GT.2) THEN
           wPast(1)=((tPast(2)-t)/tlim_dJvL)**2+2*((tPast(2)-t)/tlim_dJvL)+1
           wPast(2:lPast-1)=((tPast(3:lPast)-t)/tlim_dJvL)**2-((tPast(2:lPast-1)-t)/tlim_dJvL)**2+2*((tPast(3:lPast)-t)/tlim_dJvL)-2*((tPast(2:lPast-1)-t)/tlim_dJvL)
           wPast(lPast)=-((tPast(lPast)-t)/tlim_dJvL)**2-2*((tPast(lPast)-t)/tlim_dJvL)
        ELSEIF (lPast.EQ.2) THEN
           wPast(1)=0.0_sp
           wPast(2)=1.0_sp
        ELSE
           wPast(1)=1.0_sp
        ENDIF
        DO i=1,nElm
           !> Calculation of the layer's Hseq, here called HElm
           HElm(i)=Fh_from_mfp_soiltab(Fmfp_soiltab(Fh_from_Th(SUM(theta(elmnod(1:8,i)))/8.0_sp,par(:,MatNum(elmnod(1,i)))),MatNum(elmnod(1,i)),par(:,MatNum(elmnod(1,i))))-sink_cube(i)/Rho(i),MatNum(elmnod(1,i)),par(:,MatNum(elmnod(1,i))))+(zgrid(elmnod(1,i))+zgrid(elmnod(8,i)))/2.0_sp
        END DO
     ELSE
        DO i=1,nElm
           HElm(i)=Fh_from_Th(SUM(theta(elmnod(1:8,i)))/8.0_sp,par(:,MatNum(elmnod(1,i))))+(zgrid(elmnod(1,i))+zgrid(elmnod(8,i)))/2.0_sp	!> water potential of the mean node water content
        END DO
     ENDIF
     Hseq=DOT_PRODUCT(HElm(1:nElm),SSF)
     IF (BCtp.EQ.2) THEN
        PHcollar=-ABS(BCr)/Krs+Hseq
        Tpot=ABS(BCr)									![L³/T]
        stressBC=.FALSE.
        CALL RootStressCouvreur(PHcollar)
     ELSE
        PHCollar=BCr
     ENDIF
     IF (BCtp.EQ.1.OR.stressBC) THEN
        Tact(1)=Krs*(PHcollar-Hseq)		!> Convert the collar pressure head into the actual transpiration
     ELSE
        Tact(1)=Tpot									![L³/T]
     ENDIF
     Phi=Kcomp*(HElm-Hseq)
     sink_cube=ABS(Tact(1))*SSF/dxGrid/dyGrid/dzGrid+Phi*SSF/dxGrid/dyGrid/dzGrid

     IF (ldJvL) THEN
        !> Calculation of the layer's Hseq, here called HElm
        DO i=1,nElm
           HElm2(i)=Fh_from_mfp_soiltab(Fmfp_soiltab(Fh_from_Th(SUM(theta(elmnod(1:8,i)))/8.0_sp,par(:,MatNum(elmnod(1,i)))),MatNum(elmnod(1,i)),par(:,MatNum(elmnod(1,i))))-sink_cube(i)/Rho(i),MatNum(elmnod(1,i)),par(:,MatNum(elmnod(1,i))))+(zgrid(elmnod(1,i))+zgrid(elmnod(8,i)))/2.0_sp
        END DO
        Phi=Kcomp*(HElm-Hseq)
        sink_cube=0.5_sp*sink_cube+0.5_sp*(ABS(Tact(1))*SSF/dxGrid/dyGrid/dzGrid+Phi*SSF/dxGrid/dyGrid/dzGrid)
        HElm=(HElm+HElm2)/2.0_sp
     ENDIF
  ELSEIF (lSUF) THEN
     DO ipl=1,nplant 		! Code not ready for multiple plants yet
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
        print *,'SUM(SSF) Krs',SUM(SSF(1:nrec)),Krs
        SSF(1:nrec)=SSF(1:nrec)/SUM(SSF(1:nrec))
        Hseq=DOT_PRODUCT(PHs(1:nrec,ipl),SSF)
        IF (BCtp.EQ.2) THEN
           PHcollar=-ABS(BCr)/Krs+Hseq
           Tpot=ABS(BCr)				![L³/T]
           stressBC=.FALSE.
           CALL RootStressCouvreur(PHcollar)
        ELSE
           PHCollar=BCr
        ENDIF
        IF (BCtp.EQ.1.OR.stressBC) THEN
           Tact(1)=Krs*(PHcollar-Hseq)		![L³/T] Convert the collar pressure head into the actual transpiration 
        ELSE
           Tact(1)=Tpot					![L³/T]
        ENDIF
        Joutr(1:nrec)=(ABS(Tact(1))+Kcomp*(PHs(1:nrec,ipl)-Hseq))*SSF(1:nrec)
        sink_cube(1:nElm)=0.0
        DO irecn=1,nrec
           DO isub=1,nsub(irecn,ipl)
              c_i=cube_i(irecn,isub,ipl)
              sink_cube(c_i)=sink_cube(c_i)+Joutr(irecn)*w_sub(irecn,isub,ipl)/(dxgrid*dygrid*dzgrid)
           END DO
        END DO
     END DO
  ELSEIF (lno_RWU) THEN
     DO i=1,nPt
        sink(i)=0.0_dp
     END DO
  ENDIF
  IF (lChem) THEN
     DO i=1,nPt
        IF (betac(i).GT.1.E-20_dp) THEN
           active=TotSur*betac(i)*(VMax/(CMm+Conc(i))+fk)
        ELSE
           active=0.0_dp
        END IF
        csink(i) = xin*sink(i)+(1-xin)*active
     END DO
  ENDIF
  DO iE=1,nElm
     IF (lDou.and.lChem) THEN
        IF (betac_cube(iE).GT.1.E-20_dp) THEN
           DO iSE=1,5
              i=elmnod(iL(1,iSE,subN(iE)),iE)
              j=elmnod(iL(2,iSE,subN(iE)),iE)
              k=elmnod(iL(3,iSE,subN(iE)),iE)
              l=elmnod(iL(4,iSE,subN(iE)),iE)
              ConcE = (conc(i) + conc(j) + conc(k) + conc(l) )/4._sp
              active_cube=TotSur*betac_cube(iE)*(VMax/(CMm+ConcE)+fk)
           ENDDO
        ELSE
           active_cube=0.0_dp
        END IF
        csink_cube(iE) = xin*sink_cube(iE)+(1-xin)*active_cube
     ENDIF
  END DO
  RETURN
END SUBROUTINE SetSnk
!**************************************************************************************
FUNCTION alpha(h,Conc,tem)
  USE typedef
  USE PlntData
  USE RootData, ONLY: lFed,lDou
  USE DoussanMat, ONLY: stresfun
  IMPLICIT NONE
  REAL(sp), INTENT(in) ::conc,tem,h
  REAL(sp):: alpha,alphap,alphah,p,rktemp
  !>  Calculate osmotic potential [L] from soil solution concentration [M/L3].
  !>     To use the first or the second method simply render active the lines
  !>     where necessary quantities are calculated for the chosen expression and
  !>    comment out the lines related to the other expression; recompile the
  !>     program; beware of units and, if needed, make proper aggiustments.
  !>     Also, adjust for the valency of the chemical species studied if the
  !>     second method is adopted.
  !> **************************
  !>     For very dilute solutions the van´t Hoff expression can be used (T.L.
  !>     Allen and R.M. Keefer, "Chemistry: Experiment and Theory", 1974,
  !>     Harper and Row Publishers, pp 283-285):
  !>                              p=CRT
  !>     where:
  !>     p is the osmotic potential [L], C is the concentration [M/L3], R is a
  !>     constant [(L3*pressure)/(M*temperature)], and T is the temperature
  !>     (K degrees)
  !>
  !>     first transform temperature degrees from Celsius to Kelvin...
  !
  rKtemp=tem+273.15_dp
  !
  !>     1 liter = 1000 cm3; 1 atm = 1033.6 cm H2O;
  !>     R = 0.082 (liter atm K-1 mole-1) = 84811.0144 (cm3 cm K-1 mole-1);
  !
  p=(-Conc)*rKtemp*84811.0144_dp
  !
  !> **********************************
  !>     The following constants are taken from ´Diagnosis and Improvement of
  !>     Saline and Alkali Soils´, USDA Agriculture Handbook no. 60, 1954:
  !>               1 atm = 0.36 millimhos/cm           (Fig 6, p 15);
  !>               1 millimhos/cm = -0.091 meq/liter   (Fig 4, p 12);
  !>     these values are valid for the range of EC that allows plant growth.
  !>    Therefore 1 atm = -0.3276 meq/liter; also:
  !>               1 atm = 1033.6 cm H2O;     and:
  !>               1 meq/liter = 1.0E+06/n M/cm3, being n the valency of the      ion
  !>                                                        in the soil solution.
  !>
  !>      p=Conc*(-0.3276)*1000000/1033.6
  !>
  !> ********************************************************************************
  alpha=0.0_dp
  IF (stresfun.EQ.2) THEN
     alphah=1._dp/(1._dp+(h/h50)**p1)
     IF (p50.LT.9.E+20_dp) THEN
        alphap=1._dp/(1._dp+(p+p50)**p2)
        alpha=alphah*alphap
     ELSE
        alpha=alphah
     ENDIF
  ELSEIF (stresfun.EQ.1) THEN
     IF ((h.GT.h3).AND.(h.LT.h2)) THEN
        alpha=(h-h3)/(h2-h3)
     ELSEIF ((h.GE.h2).AND.(h.LE.h1)) THEN
        alpha=1.0_dp
     ELSEIF ((h.GT.h1).AND.(h.LT.h0)) THEN
        alpha=(h-h0)/(h1-h0)
     ELSE
        alpha=1.E-9_dp
     ENDIF
  ELSE
     alpha=1.0_dp
  ENDIF
  RETURN
END FUNCTION alpha
!**************************************************************************
SUBROUTINE RootstressCouvreur(PHtop)
  USE Typedef
  USE DoussanMat, ONLY: stressBC,hx_min,stresfun
  USE GridData, ONLY: RelEps,epslonR,factorRelEps
  
  IMPLICIT NONE
  REAL(dp), INTENT(inout):: PHtop
  REAL(dp):: del_PHr
  
  !> check if the root collar abs(PH) is larger than abs(hx_min)+tolerance and adapt the collar BC
  IF (RelEps) THEN 
     del_PHr=-MAX(ABS(hx_min/factorRelEps),epslonR)
  ELSE
     del_PHr=-epslonR
  ENDIF
  IF ((stresfun.EQ.1).AND.(PHtop<hx_min+del_PHr)) THEN
     !> top node at lower PH than allowed: start of stressed conditions
     PRINT *,'stress in the collar xylem: PHtop=',PHtop,' is lower than criterion ',hx_min
     stressBC=.TRUE.
     PHtop=hx_min
  ENDIF
END SUBROUTINE RootstressCouvreur
!********************************************************************************
