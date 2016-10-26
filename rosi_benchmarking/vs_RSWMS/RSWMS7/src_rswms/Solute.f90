! ==============================================================================
! Source file SOLUTE |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
      SUBROUTINE SOLUTE(dt,t,tPulse,dtMaxC,KodCB,IAD,IADN,IADD,iCount)
      USE BoundData
      !USE ParamData, ONLY: maxbnd
      USE GridData, ONLY: nband,nPt,nElm,elmnod,iL,subN,sink,deter,csink,ci,di,bi,ax,ay,az
      USE SolData, ONLY: Dispxx,Dispxy,Dispyy,Dispyz,Dispzz,Dispxz,Fc,Gc,conc,Nlevel,matnum,chpar,vx,vy,vz,epsi,theta_old,theta,cono,con,hnew,hold,kode
      USE WatFun
      USE CumData
      USE MatData, ONLY: As,Bs,A1
      USE DoussanMat, ONLY : sink_cube, betac_cube, betac_cube2, csink_cube
      IMPLICIT NONE

      REAL(sp),INTENT(in)::dt,t,tpulse
      REAL(sp),INTENT(inout) :: dtmaxc
      INTEGER(ap),INTENT(inout) :: KodCB(maxbdr),iCount
      REAL(sp) :: VxE(4),VyE(4),VzE(4),S(4,4)
      INTEGER(ap):: i,j,k,l,m,iE,iSE,ic,j1,i1,j2,i2,ib,lev,List(4),kk
      REAL(sp) :: Vzz,Vyy,Vxx
      REAL(sp) :: ec1,gce,fmul,rootch,vzee,vyee,vxee,cone,VEl,cayz,s_m
      REAL(sp) :: smul2,smul1,AcE,fce,ec6,ec5,ec4,ec3,ec2,CAzz,alf,CAyy
      REAL(sp) :: dpom,CAxz,CAxy,CAxx
      INTEGER(ap):: IAD(nband,nPt),IADN(nPt),IADD(nPt),north
      !REAL(sp)::Wx(4),Wy(4),Wz(4),W12,W13,W14,W23,W24,W34
      !REAL(sp)::A11,A12,A13,A14,A21,A22,A23,A24,A31,A32,A33,A34,A41,A42,A43,A44
      REAL(sp), ALLOCATABLE,DIMENSION (:) :: Ac
      REAL(sp), ALLOCATABLE,DIMENSION (:) :: B1, F ,DS
        
      ALLOCATE(As(nband,nPt))
      ALLOCATE(A1(nband,nPt))
      ALLOCATE(B1(nPt),F(nPt),DS(nPt))
      ALLOCATE(Ac(nPt),Bs(nPt))
      IF (.NOT.ALLOCATED(Dispxx)) THEN
         ALLOCATE(Dispxx(nPt),Dispyy(nPt),Dispzz(nPt))
         ALLOCATE(Dispxy(nPt),Dispxz(nPt),Dispyz(nPt))
      ENDIF
      IF(.NOT.ALLOCATED(Fc)) ALLOCATE(Fc(nPt),Gc(nPt))
      IF(.NOT.ALLOCATED(Qc)) ALLOCATE(Qc(nPt))
      Gc = 0._sp
      Fc = 0._sp
      Ac = 0._sp
      Qc = 0._sp
      F =  0._sp
      DS = 0._sp
      Bs = 0._sp
!**********************************************************************************
      IF (iCount.EQ.0) CALL ChInit(dtMaxC,dt)
      ! Initialuization
      alf=1._sp-epsi
     
      IF (t.GT.tPulse) THEN
         cBound=0.0_dp
      ENDIF
      Bs =0.0_dp
      B1 = Conc
      Qc=0.0_dp
      DO 13 i=1,Npt
         IF (epsi.LT.0.001_dp) THEN
            As(iadd(i),i)=0.0_dp
         ELSE
            DO 12 j=1,nband!2*Nband-1
               As(j,i)=0.0_sp
12          CONTINUE
         ENDIF
13    CONTINUE

      DO 21 Lev=1,NLevel
         IF (Lev.EQ.NLevel) THEN
            CALL Veloc
            CALL Disper
            CALL PeCour(dtMaxC,dt)
            !IF (lUpW) CALL WeFact
         ELSE
            CALL Disper
         ENDIF
    
         DO 14 i=1,Npt
            M=MatNum(i)
            IF (Lev.NE.NLevel) THEN
               !if(.not.lUpW) then
                  DPom=dt/6._sp/(theta_old(i)+ChPar(1,M)*ChPar(5,M))
                  Dispxx(i)=Dispxx(i)+Vx(i)*Vx(i)*DPom
                  Dispyy(i)=Dispyy(i)+Vy(i)*Vy(i)*DPom
                  Dispzz(i)=Dispzz(i)+Vz(i)*Vz(i)*DPom
                  Dispxy(i)=Dispxy(i)+Vx(i)*Vy(i)*DPom
                  Dispxz(i)=Dispxz(i)+Vx(i)*Vz(i)*DPom
                  Dispyz(i)=Dispyz(i)+Vy(i)*Vz(i)*DPom
               !end if
            ELSE
               Ac(i)=-(theta_old(i)*alf+theta(i)*epsi)-ChPar(1,M)*ChPar(5,M)
               !if(.not.lUpW) then
                  DPom=dt/6._sp/(theta(i)+ChPar(1,M)*ChPar(5,M))
                  Dispxx(i)=Dispxx(i)-Vx(i)*Vx(i)*DPom
                  Dispyy(i)=Dispyy(i)-Vy(i)*Vy(i)*DPom
                  Dispzz(i)=Dispzz(i)-Vz(i)*Vz(i)*DPom
                  Dispxy(i)=Dispxy(i)-Vx(i)*Vy(i)*DPom
                  Dispxz(i)=Dispxz(i)-Vx(i)*Vz(i)*DPom
                  Dispyz(i)=Dispyz(i)-Vy(i)*Vz(i)*DPom
               !end if
               Gc(i)=ChPar(8,M)*theta(i)+ChPar(1,M)*ChPar(9,M)
               Fc(i)=ChPar(6,M)*theta(i)+ChPar(1,M)*ChPar(7,M)*ChPar(5,M)-csink(i)+sink(i)
            ENDIF
14       CONTINUE
         
         DO 15 i=1,Npt    
            IF (Lev.EQ.NLevel) DS(i)=0.0_dp
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
          DO 18 iSE=1,5
             i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
             j=elmnod(iL(2,iSE,subN(iE)),iE)
             k=elmnod(iL(3,iSE,subN(iE)),iE)
             l=elmnod(iL(4,iSE,subN(iE)),iE)
             List(1)=i
             List(2)=j
             List(3)=k
             List(4)=l
             VEl=ABS(Deter(iSE,iE))/6._sp		!Deter, Ax, Ay, Az from calcgeom (Couvreur oct 2010)
!            Calculate Velocities
             IF (Lev.EQ.NLevel) THEN
               Vxx=(Ax(1,iSE,iE)*hNew(i)+Ax(2,iSE,iE)*hNew(j)+Ax(3,iSE,iE)*hNew(k)+Ax(4,iSE,iE)*hNew(l))/Deter(iSE,iE)+CAxz
               Vyy=(Ay(1,iSE,iE)*hNew(i)+Ay(2,iSE,iE)*hNew(j)+Ay(3,iSE,iE)*hNew(k)+Ay(4,iSE,iE)*hNew(l))/Deter(iSE,iE)+CAyz
               Vzz=(Az(1,iSE,iE)*hNew(i)+Az(2,iSE,iE)*hNew(j)+Az(3,iSE,iE)*hNew(k)+Az(4,iSE,iE)*hNew(l))/Deter(iSE,iE)+CAzz
             ELSE
               Vxx=(Ax(1,iSE,iE)*hOld(i)+Ax(2,iSE,iE)*hOld(j)+Ax(3,iSE,iE)*hOld(k)+Ax(4,iSE,iE)*hOld(l))/Deter(iSE,iE)+CAxz
               Vyy=(Ay(1,iSE,iE)*hOld(i)+Ay(2,iSE,iE)*hOld(j)+Ay(3,iSE,iE)*hOld(k)+Ay(4,iSE,iE)*hOld(l))/Deter(iSE,iE)+CAyz
               Vzz=(Az(1,iSE,iE)*hOld(i)+Az(2,iSE,iE)*hOld(j)+Az(3,iSE,iE)*hOld(k)+Az(4,iSE,iE)*hOld(l))/Deter(iSE,iE)+CAzz
             END IF
             IF (Lev.NE.NLevel) THEN
                ConE=(ConO(i)+ConO(j)+ConO(k)+ConO(l))/4._sp
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
             ELSE
                ConE=(Con(i)+Con(j)+Con(k)+Con(l))/4._sp
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
             END IF
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
             FMul=VEl/5._sp
             GcE=(Gc(i)+Gc(j)+Gc(k)+Gc(l))/4._sp
             Ec1=(Dispxx(i)+Dispxx(j)+Dispxx(k)+Dispxx(l))/4._sp
             Ec2=(Dispyy(i)+Dispyy(j)+Dispyy(k)+Dispyy(l))/4._sp
             Ec3=(Dispzz(i)+Dispzz(j)+Dispzz(k)+Dispzz(l))/4._sp
             Ec4=(Dispxy(i)+Dispxy(j)+Dispxy(k)+Dispxy(l))/4._sp
             Ec5=(Dispxz(i)+Dispxz(j)+Dispxz(k)+Dispxz(l))/4._sp
             Ec6=(Dispyz(i)+Dispyz(j)+Dispyz(k)+Dispyz(l))/4._sp
             IF (Lev.EQ.NLevel) AcE=(Ac(i)+Ac(j)+Ac(k)+Ac(l))/4._sp
             FcE=(Fc(i)+Fc(j)+Fc(k)+Fc(l))/4._sp
             SMul1=-1._sp/VEl/36._sp
             SMul2=VEl/30._sp
!!$                  if(lUpW) then
!!$                     W12=WeTab(1,iSE)
!!$                     W13=WeTab(2,iSE)
!!$                     W14=WeTab(3,iSE)
!!$                     W23=WeTab(4,iSE)
!!$                     W24=WeTab(5,iSE)
!!$                     W34=WeTab(6,iSE)
!!$                     A11=-2.*W12+2.*W14+2.*W13
!!$                     A12=-2.*W12+W14+W13
!!$                     A13=-W12+W14+2.*W13
!!$                     A14=-W12+2.*W14+W13
!!$                     A21=-W23+2.*W12+W24
!!$                     A22=-2.*W23+2.*W12+2.*W24
!!$                     A23=-2.*W23+W12+W24
!!$                     A24=-W23+W12+2.*W24
!!$                     A31=-W34+W23-2.*W13
!!$                     A32=-W34+2.*W23-W13
!!$                     A33=-2.*W34+2.*W23-2.*W13
!!$                     A34=-2.*W34+W23-W13
!!$                     A41=-2.*W14+W34-W24
!!$                     A42=-W14+W34-2.*W24
!!$                     A43=-W14+2.*W34-W24
!!$                     A44=-2.*W14+2.*W34-2.*W24
!!$                     Wx(1)=VxE(1)*A11+VxE(2)*A12+VxE(3)*A13+VxE(4)*A14
!!$                     Wx(2)=VxE(1)*A21+VxE(2)*A22+VxE(3)*A23+VxE(4)*A24
!!$                     Wx(3)=VxE(1)*A31+VxE(2)*A32+VxE(3)*A33+VxE(4)*A34
!!$                     Wx(4)=VxE(1)*A41+VxE(2)*A42+VxE(3)*A43+VxE(4)*A44
!!$                     Wy(1)=VyE(1)*A11+VyE(2)*A12+VyE(3)*A13+VyE(4)*A14
!!$                     Wy(2)=VyE(1)*A21+VyE(2)*A22+VyE(3)*A23+VyE(4)*A24
!!$                     Wy(3)=VyE(1)*A31+VyE(2)*A32+VyE(3)*A33+VyE(4)*A34
!!$                     Wy(4)=VyE(1)*A41+VyE(2)*A42+VyE(3)*A43+VyE(4)*A44
!!$                     Wz(1)=VzE(1)*A11+VzE(2)*A12+VzE(3)*A13+VzE(4)*A14
!!$                     Wz(2)=VzE(1)*A21+VzE(2)*A22+VzE(3)*A23+VzE(4)*A24
!!$                     Wz(3)=VzE(1)*A31+VzE(2)*A32+VzE(3)*A33+VzE(4)*A34
!!$                     Wz(4)=VzE(1)*A41+VzE(2)*A42+VzE(3)*A43+VzE(4)*A44
!!$           end if
             DO 17 j1=1,4
                i1=List(j1)
                !F(i1)=F(i1)+FMul*(GcE+Gc(i1)/4._sp) + Vel/4._sp * csink_cube(iE)!!f_n
                F(i1)=F(i1)+FMul*(GcE+Gc(i1)/4._sp) !!f_n
                IF(Lev.EQ.NLevel) DS(i1)=DS(i1)+FMul*(AcE+Ac(i1)/4._sp)!Q_nm
                S_m=0.0_dp
                DO 16 j2=1,4
                   i2=List(j2)
                   S(j1,j2)=SMul1*(Ec1*bi(j1,iSE,iE)*bi(j2,iSe,iE)+Ec2*ci(j1,iSE,iE)*ci(j2,iSE,iE)+ Ec3*di(j1,iSe,iE)*di(j2,iSe,iE)+&
                     Ec4*(bi(j1,iSE,iE)*ci(j2,iSE,iE)+ci(j1,iSE,iE)*bi(j2,iSE,iE))+Ec5*(bi(j1,iSE,iE)*di(j2,iSE,iE)+di(j1,iSE,iE)*bi(j2,iSE,iE))+&
                     Ec6*(ci(j1,iSE,iE)*di(j2,iSE,iE)+di(j1,iSE,iE)*ci(j2,iSE,iE)))
                   S(j1,j2)=S(j1,j2)-(6*VEl/Deter(iSE,iE))*((bi(j2,iSE,iE)/30._sp)*(VxEE+VxE(j1)/4._sp)+&  !(6*VEl/Deter(iE,iSE))*
                     (ci(j2,iSE,iE)/30._sp)*(VyEE+VyE(j1)/4._sp)+(di(j2,iSE,iE)/30._sp)*(VzEE+VzE(j1)/4._sp))

                   !if(lUpW) S(j1,j2)=S(j1,j2)-(bi(iE,iSE,j2)/240.*Wx(j1)+ci(iE,iSE,j2)/240.*Wy(j1)+di(iE,iSE,j2)/240.*Wz(j1))

                   ic=1
                   IF (i1.EQ.i2) ic=2
                   S(j1,j2)=S(j1,j2)+SMul2*ic*(FcE+(Fc(i1)+Fc(i2))/4._sp)
                   !S(j1,j2)=S(j1,j2)+ic*(SMul2*(FcE+(Fc(i1)+Fc(i2))/4._sp) + VEL/24._sp*sink_cube(iE))
                   S_m = S_m -(epsi*S(j1,j2))*Conc(i2)
                   IF (Lev.NE.NLevel) THEN
                      Bs(i1)=Bs(i1)-alf*S(j1,j2)*Conc(i2)
                   ELSE
                      CALL Find(i1,i2,kk,nPt,nband,IAD,IADN)
                      iB=kk
                      !iB=iadd_temp(k)
                      As(iB,i1)=As(iB,i1)+epsi*S(j1,j2)
                   ENDIF
                   IF (Lev.EQ.1.AND.Kode(i1).GT.0) Qc(i1)=Qc(i1)-S(j1,j2)*Conc(i2) 
16              CONTINUE
17                 CONTINUE
18        CONTINUE !sub-element loop
19       CONTINUE  !element loop

         DO 20 i=1,Npt
            M=MatNum(i)
            IF (Lev.EQ.1.AND.Kode(i).GT.0) Qc(i)=Qc(i)-F(i)
            IF (Lev.NE.NLevel) THEN
               Bs(i)=Bs(i)-alf*F(i)
            ELSE
               As(iadd(i),i)=As(iadd(i),i)+DS(i)/dt
               Bs(i)=Bs(i)+DS(i)/dt*Conc(i)-epsi*F(i)
            ENDIF
20       CONTINUE
21    CONTINUE
  
      ! Boundary condition 
      CALL C_Bound(KodCB,dt,DS,IADD)
      ! Solve the global matrix equation for transport
      IF (epsi.LT.0.001_dp) THEN
         DO 22 i=1,Npt
            Bs(i)=Bs(i)/As(iadd(i),i)
22       CONTINUE
      ELSE 
         CALL ILU (As,nPt,nband,IAD,IADN,IADD,A1)
         North=4!north=0 (symmetric matrix -> is not the case)
         CALL OrthoMin(As,B1,Bs,nPt,nband,nPt,IAD,IADN,IADD,A1,North)
         !CALL SOLVET
      ENDIF
      Conc=B1
      Bs = B1
      DO 23 i=1,Npt
        
         !sngl(B(i))
         IF (Conc(i).LT.1.0E-30_dp) Conc(i)=0.0_dp
23    CONTINUE
      CALL SolInf(dt)
! ====================================
      DEALLOCATE(As,A1,Bs)
      DEALLOCATE(B1,F,DS)
      DEALLOCATE(Ac,Fc,Gc,Qc)
      DEALLOCATE(Dispxx,Dispyy,Dispzz)
      DEALLOCATE(Dispxy,Dispxz,Dispyz)
      RETURN
      END SUBROUTINE SOLUTE
!*******************************************************************
      SUBROUTINE C_Bound(KodCB,dt,DS,IADD)
      USE BoundData
      USE GridData, ONLY: nBCPts
      USE SolData, ONLY: epsi,Kode,conc
      USE CumData
      USE MatData
      IMPLICIT NONE
      REAL(sp),INTENT(in) :: dt,DS(nPt)
      INTEGER(ap),INTENT(inout) ::KodCB(nPt)
      REAL(sp) :: alf,cBnd
      INTEGER(ap) :: cKod,i,j
      INTEGER(ap) :: IADD(nPt)
      alf=1._sp-epsi
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
             As(iadd(i),i)=1._dp
             Bs(i)=1._dp*cBnd
             
!!$             DO 13 j=1,2*NBand-1
!!$                  As(j,i)=0.0_dp
!!$13           CONTINUE
!!$             As(NBand,i)=1._dp
!!$             Bs(i)=cBnd
            
          ENDIF
! Neumann boundary condition
          
          IF (cKod.EQ.2) THEN
             Qc(i)=Q(i)*Conc(i)
            
          ENDIF
! Cauchy boundary condition

          IF (cKod.EQ.3) THEN
             Bs(i)=Bs(i)-Q(i)*(cBnd-alf*Conc(i))
             As(iadd(i),i)=As(iadd(i),i)-epsi*Q(i)
             Qc(i)=Q(i)*cBnd
          
          ENDIF
       ENDIF
14    CONTINUE
      RETURN
      END SUBROUTINE C_Bound
!**********************************************************************************
      SUBROUTINE ChInit(dtMaxC,dt)
      USE GridData
      USE SolData
      USE WatFun
      IMPLICIT NONE
      REAL(sp),INTENT(in) :: dt
      REAL,INTENT(inout) :: dtMaxC
      INTEGER(ap) ::  i,M
      IF (.NOT.ALLOCATED(Fc)) ALLOCATE(Fc(nPt),Gc(nPt))
     
         DO 11 i=1,nPt
            IF (NLevel.EQ.2) THEN
            M=MatNum(i)
            Gc(i)=ChPar(8,M)*theta(i)+ChPar(1,M)*ChPar(9,M)
            Fc(i)=ChPar(6,M)*theta(i)+ChPar(1,M)*ChPar(7,M)*ChPar(5,M)+(sink(i)-csink(i))
         ENDIF
11       CONTINUE
     
      CALL Veloc
      CALL Disper
      CALL PeCour(dtMaxC,dt)
      RETURN
      END SUBROUTINE ChInit
!*************************************************************************
      SUBROUTINE Veloc
      USE GridData
      USE SolData
      USE CumData
      IMPLICIT NONE
    
      INTEGER(ap) :: i,j,k,l,m,ise,ie,List(4)
      REAL(sp) :: vxx,vyy,vzz,cayz,caxz,caxy,cazz,cayy,caxx, VE, A
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
       DO 13  iSE=1,5
          i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
          j=elmnod(iL(2,iSE,subN(iE)),iE)
          k=elmnod(iL(3,iSE,subN(iE)),iE)
          l=elmnod(iL(4,iSE,subN(iE)),iE)
          List(1)=i
          List(2)=j
          List(3)=k
          List(4)=l
          VE=ABS(Deter(iSE,iE))/6._sp
          A=1./VE/6
            Vxx=A*(Ax(1,iSE,iE)*hNew(i)+Ax(2,iSE,iE)*hNew(j)+Ax(3,iSE,iE)*hNew(k)+Ax(4,iSE,iE)*hNew(l))	!Deter, Ax, Ay and Az are now globals (Couvreur mar 2010)
            Vxx=Vxx+CAxz
            Vyy=A*(Ay(1,iSE,iE)*hNew(i)+Ay(2,iSE,iE)*hNew(j)+Ay(3,iSE,iE)*hNew(k)+Ay(4,iSE,iE)*hNew(l))
            Vyy=Vyy+CAyz
            Vzz=A*(Az(1,iSE,iE)*hNew(i)+Az(2,iSE,iE)*hNew(j)+Az(3,iSE,iE)*hNew(k)+Az(4,iSE,iE)*hNew(l))
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
      SUBROUTINE Disper
      USE GridData
      USE SolData
      USE WatFun
      IMPLICIT NONE
     
      REAL(sp) :: Tau,Vabs,ThSati
      INTEGER(ap) :: i,M
      IF (.NOT.ALLOCATED(Dispxx)) THEN 
         ALLOCATE(Dispxx(nPt),Dispyy(nPt),Dispzz(nPt))
         ALLOCATE(Dispxy(nPt),Dispxz(nPt),Dispyz(nPt))
      ENDIF
      DO 11 i=1,nPt
       M=MatNum(i)
       IF (soiltab) THEN	!(Couvreur nov 2011)
          ThSati=(TheTab(1,M)-TheTab(nTab,M))*Dxy(i)+TheTab(nTab,M)*Exy(i)
       ELSE
          ThSati=(par(3,M)-par(2,M))*Dxy(i)+par(2,M)*Exy(i)
       ENDIF
       Tau=theta(i)**(7._dp/3._dp)/ThSati**2
       Vabs=SQRT(Vx(i)*Vx(i)+Vy(i)*Vy(i)+Vz(i)*Vz(i))
       IF (Vabs.GT.0._dp) THEN
          IF (soiltab) THEN	!(Couvreur nov 2011)
             Dispxx(i)=ChPar(3,M)*Vx(i)*Vx(i)/Vabs+ChPar(4,M)*(Vz(i)*Vz(i)+Vy(i)*Vy(i))/Vabs+Fth_soiltab(hNew(i),M)*ChPar(2,M)*Tau
             Dispyy(i)=ChPar(3,M)*Vy(i)*Vy(i)/Vabs+ChPar(4,M)*(Vx(i)*Vx(i)+Vz(i)*Vz(i))/Vabs+Fth_soiltab(hNew(i),M)*ChPar(2,M)*Tau
             Dispzz(i)=ChPar(3,M)*Vz(i)*Vz(i)/Vabs+ChPar(4,M)*(Vx(i)*Vx(i)+Vy(i)*Vy(i))/Vabs+Fth_soiltab(hNew(i),M)*ChPar(2,M)*Tau
          ELSE
             Dispxx(i)=ChPar(3,M)*Vx(i)*Vx(i)/Vabs+ChPar(4,M)*(Vz(i)*Vz(i)+Vy(i)*Vy(i))/Vabs+Fth(hNew(i),par(:,M),i)*ChPar(2,M)*Tau
             Dispyy(i)=ChPar(3,M)*Vy(i)*Vy(i)/Vabs+ChPar(4,M)*(Vx(i)*Vx(i)+Vz(i)*Vz(i))/Vabs+Fth(hNew(i),par(:,M),i)*ChPar(2,M)*Tau
             Dispzz(i)=ChPar(3,M)*Vz(i)*Vz(i)/Vabs+ChPar(4,M)*(Vx(i)*Vx(i)+Vy(i)*Vy(i))/Vabs+Fth(hNew(i),par(:,M),i)*ChPar(2,M)*Tau
          ENDIF
          Dispxy(i)=(ChPar(3,M)-ChPar(4,M))*Vx(i)*Vy(i)/Vabs
          Dispxz(i)=(ChPar(3,M)-ChPar(4,M))*Vx(i)*Vz(i)/Vabs
          Dispyz(i)=(ChPar(3,M)-ChPar(4,M))*Vy(i)*Vz(i)/Vabs
       ELSE
          IF (soiltab) THEN	!(Couvreur nov 2011)
             Dispxx(i)=Fth_soiltab(hNew(i),M)*ChPar(2,M)*Tau
             Dispyy(i)=Fth_soiltab(hNew(i),M)*ChPar(2,M)*Tau
             Dispzz(i)=Fth_soiltab(hNew(i),M)*ChPar(2,M)*Tau
          ELSE
             Dispxx(i)=Fth(hNew(i),par(:,M),i)*ChPar(2,M)*Tau
             Dispyy(i)=Fth(hNew(i),par(:,M),i)*ChPar(2,M)*Tau
             Dispzz(i)=Fth(hNew(i),par(:,M),i)*ChPar(2,M)*Tau
          ENDIF
          Dispxy(i)=0.0_dp
          Dispxz(i)=0.0_dp
          Dispyz(i)=0.0_dp
       ENDIF
11    CONTINUE
      RETURN
      END SUBROUTINE Disper
!****************************************************************
!!$      SUBROUTINE SolveT
!!$      USE GridData
!!$      USE MatData, ONLY: As, Bs
!!$      IMPLICIT NONE
!!$      REAL(sp) :: P,C,Sum
!!$      INTEGER(ap):: k,N1,i,kk,kc,ii,j,L,jj,M
!!$      N1=Npt-1
!!$      DO 12 k=1,N1
!!$         P=1._sp/As(NBand,k)
!!$         kk=k+1
!!$         kc=NBand
!!$         DO 11 i=kk,Npt
!!$            kc=kc-1
!!$            IF (kc.LE.0) GOTO 12
!!$            C=-P*As(kc,i)
!!$            As(kc,i)=C
!!$            ii=kc+1
!!$            L=kc+NBand-1
!!$            DO 11 j=ii,L
!!$             jj=j+NBand-kc
!!$             As(j,i)=As(j,i)+C*As(jj,k)
!!$11       CONTINUE
!!$12    CONTINUE
!!$      DO 14 i=2,Npt
!!$       jj=NBand+1-i
!!$       ii=1
!!$       IF (jj.LE.0) THEN
!!$          jj=1
!!$          ii=i-NBand+1
!!$       ENDIF
!!$       Sum=0.0_dp
!!$       DO 13 j=jj,NBand-1
!!$          Sum=Sum+As(j,i)*Bs(ii)
!!$          ii=ii+1
!!$13     CONTINUE
!!$       Bs(i)=Bs(i)+Sum
!!$14    CONTINUE
!!$      Bs(Npt)=Bs(Npt)/As(NBand,Npt)
!!$      DO 16 k=1,N1
!!$       i=Npt-k
!!$       jj=i
!!$       m=MIN(2*NBand-1,NBand+k)
!!$       Sum=0.0_dp
!!$       DO 15 j=NBand+1,m
!!$          jj=jj+1
!!$          Sum=Sum+As(j,i)*Bs(jj)
!!$15     CONTINUE
!!$       Bs(i)=(Bs(i)-Sum)/As(NBand,i)
!!$16    CONTINUE
!!$      RETURN
!!$      END SUBROUTINE SolveT
!****************************************************************
      SUBROUTINE PeCour(dtMaxC,dt)
      USE GridData
      USE SolData
      USE WatFun
      IMPLICIT NONE
      REAL(sp), INTENT(in) :: dt
      REAL(sp), INTENT(out):: dtMaxC
      INTEGER(ap)::i,j,k,l,ie,ise
      REAL(sp):: pecx,pecy,delx !,xmax,xmin,zmin,zmax,ymin,ymax
      REAL(sp):: r1,r2,r3,r4,rmin
      REAL(sp):: vzmax,vymax,vxmax,vze,vye,vxe,dze,dye,dxe,delz,dely
      REAL(sp):: courx,coury,courz,cour1,cour2,cour3,dt3,dt2,dt1,pecz
      Peclet=0.0_dp
      Courant=0.0_dp
      dtMaxC=1.e+30_dp
      DO 12 iE=1,Nelm
       DO 11 iSE=1,5
          PecX=99999._dp
          PecY=99999._dp
          PecZ=99999._dp
          dt1=1.e+30_dp
          dt2=1.e+30_dp
          dt3=1.e+30_dp

          i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
          j=elmnod(iL(2,iSE,subN(iE)),iE)
          k=elmnod(iL(3,iSE,subN(iE)),iE)
          l=elmnod(iL(4,iSE,subN(iE)),iE)
!            thetai=Fth(hNew(i),par(:,MatNum(i)))
!            thetaj=Fth(hNew(j),par(:,MatNum(j)))
!           thetak=Fth(hNew(k),par(:,MatNum(k)))
!           thetal=Fth(hNew(l),par(:,MatNum(l)))
!!JAVAUX AMIN1 and AMAX! have been changed to min and max
!          xmax=max(xGrid(i),xGrid(j),xGrid(k),xGrid(l))
!          xmin=min(xGrid(i),xGrid(j),xGrid(k),xGrid(l))
!          ymax=max(yGrid(i),yGrid(j),yGrid(k),yGrid(l))
!          ymin=min(yGrid(i),yGrid(j),yGrid(k),yGrid(l))
!          zmax=max(zGrid(i),zGrid(j),zGrid(k),zGrid(l))
!          zmin=min(zGrid(i),zGrid(j),zGrid(k),zGrid(l))
!            delX=xmax-xmin
!            delY=ymax-ymin
!            delZ=zmax-zmin
          delX=dxGrid				!Voxel size is constant and using xmax, xmin, etc causes problems in continuous (Couvreur oct 2010)
          delY=dyGrid
          delZ=dzGrid
          DxE=(Dispxx(i)+Dispxx(j)+Dispxx(k)+Dispxx(l))/4
          DyE=(Dispyy(i)+Dispyy(j)+Dispyy(k)+Dispyy(l))/4
          DzE=(Dispzz(i)+Dispzz(j)+Dispzz(k)+Dispzz(l))/4
          VxE=ABS(Vx(i)+Vx(j)+Vx(k)+Vx(l))/4
          VyE=ABS(Vy(i)+Vy(j)+Vy(k)+Vy(l))/4
          VzE=ABS(Vz(i)+Vz(j)+Vz(k)+Vz(l))/4
          IF (DxE.GT.0._dp) PecX=VxE*delX/DxE
          IF (DyE.GT.0._dp) PecY=VyE*delY/DyE
          IF (DzE.GT.0._dp) PecZ=VzE*delZ/DzE
!JAVAUX AMIN1 and AMAX! have been changed to min and max
          IF (PecX.NE.99999._dp) Peclet=MAX(Peclet,PecX)
          IF (PecY.NE.99999._dp) Peclet=MAX(Peclet,PecY)
          IF (PecZ.NE.99999._dp) Peclet=MAX(Peclet,PecZ)
!JAVAUX AMIN1 and AMAX! have been changed to min and max
          Peclet=MIN(Peclet,99999._sp)
          VxMax=MAX(ABS(Vx(i))/theta(i),ABS(Vx(j))/theta(j),ABS(Vx(k))/theta(k),ABS(Vx(l))/theta(l))
          VyMax=MAX(ABS(Vy(i))/theta(i),ABS(Vy(j))/theta(j),ABS(Vy(k))/theta(k),ABS(Vy(l))/theta(l))
          VzMax=MAX(ABS(Vz(i))/theta(i),ABS(Vz(j))/theta(j),ABS(Vz(k))/theta(k),ABS(Vz(l))/theta(l))
          R1=1._dp+ChPar(1,MatNum(i))*ChPar(5,MatNum(i))/theta(i)
          R2=1._dp+ChPar(1,MatNum(j))*ChPar(5,MatNum(j))/theta(j)
          R3=1._dp+ChPar(1,MatNum(k))*ChPar(5,MatNum(k))/theta(k)
          R4=1._dp+ChPar(1,MatNum(l))*ChPar(5,MatNum(l))/theta(l)
          RMin=MIN(R1,R2,R3,R4)
          CourX=VxMax*dt/(delX*RMin)
          CourY=VyMax*dt/(delY*RMin)
          CourZ=VzMax*dt/(delZ*RMin)
!JAVAUX AMIN1 and AMAX! have been changed to min and max
          Courant=MAX(Courant,CourX,CourY,CourZ)
          Cour1=1.0_dp
          Cour2=1.0_dp
          Cour3=1.0_dp
!JAVAUX AMIN1 and AMAX! have been changed to min and max
          IF(PecX.NE.99999._sp) Cour1=MIN(1._sp,PeCr/MAX(0.5_sp,PecX))
          IF(PecY.NE.99999._sp) Cour2=MIN(1._sp,PeCr/MAX(0.5_sp,PecY))
          IF(PecZ.NE.99999._sp) Cour3=MIN(1._sp,PeCr/MAX(0.5_sp,PecZ))
          IF(VxMax.GT.0._dp) dt1=Cour1*delX*RMin/VxMax
          IF(VyMax.GT.0._dp) dt2=Cour2*delY*RMin/VyMax
          IF(VzMax.GT.0._dp) dt3=Cour3*delZ*RMin/VzMax
!JAVAUX AMIN1 and AMAX! have been changed to min and max
          dtMaxC=MIN(dtMaxC,dt1,dt2,dt3)
11       CONTINUE
12    CONTINUE
      RETURN
      END SUBROUTINE PeCour
!***********************************************************************************
      SUBROUTINE SolInf(dt)
      USE GridData
      USE CumData
      USE SolData, ONLY: Kode
      IMPLICIT NONE
      REAL(sp):: sMean(2),dt
      INTEGER(ap):: i,j
      sMean(1)=0.0_dp
      sMean(2)=0.0_dp
      DO 12 i=1,nPt
         j=iabs(Kode(i))
         IF (j.NE.0) THEN
            sMean(j)=sMean(j)-Qc(i) !Qc=nodal value of solute flux, only boundary values!!
         ENDIF
12    CONTINUE
      cCumA=ABS(CumCh0)+ABS(CumCh1)+ABS(CumChR)
      cCumT=CumCh0+CumCh1+CumChR
	DO 13 j=1,2
         ChemS(j)=ChemS(j)+sMean(j)*dt
         cCumT=cCumT+ChemS(j)
         cCumA=cCumA+ABS(ChemS(j))
13    CONTINUE
      RETURN
      END SUBROUTINE SolInf
!***********************************************************************************
