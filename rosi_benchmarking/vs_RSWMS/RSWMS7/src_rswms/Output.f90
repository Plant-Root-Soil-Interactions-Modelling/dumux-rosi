! ==============================================================================
! Source file OUTPUT |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
!> generates output for observation nodes
SUBROUTINE ObsIni    
  USE Typedef
  USE ObsData
  USE GridData, ONLY : xgrid,ygrid,zgrid
  IMPLICIT NONE
  
  INTEGER(ap) :: ip,minl,i
  CHARACTER filename*13,form*3,form2*2,form1*1
  
  DO ip=1,nPr
     minl=MIN(nodebyPr(ip),500)
     !> define file names
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
           WRITE(form,'(I3)')500
        ELSE
           WRITE(10,'(/''This plane contains'',I4,'' nodes'')') nodebyPr(ip)
           IF (nodebyPr(ip).GT.99) THEN
              WRITE(form,'(I3)')nodebyPr(ip)!only valid for larger than 100 numbers!!!!
           ELSEIF (nodebyPr(ip).GT.9) THEN
              WRITE(form2,'(I2)')nodebyPr(ip)
           ELSE
              WRITE(form1,'(I1)')nodebyPr(ip)
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
           WRITE(form,'(I3)')(nodebyPr(ip)-500)
           WRITE(11,'('//form//'(1X,I5))') (NodePr(ip,i),i=501,nodebyPr(ip))
           IF (varP(ip)==2) WRITE (11,'(/''      Time   Node water content'')')
           IF (varP(ip)==1) WRITE (11,'(/''      Time   Node water potential'')')
           IF (varP(ip)==3) WRITE (11,'(/''      Time   Node water potential (line 1) and node water content (line 2)'')')
        ENDIF
     ENDIF
  ENDDO
END SUBROUTINE ObsIni
!****************************************************************************************
!> generates log file; output of log depends on modelled processes
SUBROUTINE WriteLog(t)
  USE Typedef
  USE PlntData
  USE RootData
  USE Doussanmat,ONLY : count_nodes,PHr,nplant,stresfun
  USE GridData,ONLY: rootsk,HElm,SSF,nElm
  USE ParamData, ONLY: lPartUp
  USE SoluteRootMat, Only: totalParticleNum
  USE SolData, ONLY: sum_upt
  IMPLICIT NONE
  REAL(sp):: t
  INTEGER ::ipl
  CHARACTER :: na*10,file*8
  !> \param t current time
  
  WRITE (file,'(A7)')'out/log'
  na='    n/a   '
  DO ipl=1,nplant
     WRITE (file(8:8),'(I1)') ipl
     OPEN (UNIT=10,FILE=file,STATUS='OLD',POSITION='APPEND')
     IF (lFed) THEN
        WRITE(10,'(3(1pE10.3,1X),A10,1X,2(1pE10.3,1X),A10,1X,1pE10.3)') &
             t,Tpot,Tact(ipl),na,sAvg,cAvg,na,mroot
     ELSEIF (lCalloc) THEN 
        WRITE (10,'(12(1pE11.4,1X),I6)')&
             t,Tpot,Tact(ipl),grwfac,sAvg,mshoot,mroot,LA,PHr(1,ipl),&
             concol,TpLA,TpLA_pot,count_nodes
     ELSEIF (lDou .AND. .NOT.lCalloc) THEN
        if (lSign .OR. lSign_inst) then
           WRITE (10,'(4(1pE11.4,1X),I6,1X,9(1pE11.4,1X))')&
                t,Tpot,Tact(ipl),sAvg,count_nodes,PHr(1,ipl),mcol,concol,&
                msign_notrans,csign_notrans,res_t
        elseif(lPartUp) then
           WRITE (10,'(4(1pE11.4,1X),I6,1X,3(1pE11.4,1X), 1I9,1X, 1pE11.4, 1X)')&
                t,Tpot,Tact(ipl),sAvg,count_nodes,PHr(1,ipl),mcol,concol, totalParticleNum, sum_upt
        else
           WRITE (10,'(4(1pE11.4,1X),I6,1X,1(1pE11.4,1X))')&
                t,Tpot,Tact(ipl),sAvg,count_nodes,PHr(1,ipl)
        endif
     ELSEIF (lCou) THEN
        WRITE (10,'(3(1pE11.4,1X),I2,1X,1pE11.4,1X,1pE11.4)')&
             t,Tpot,Tact(ipl),0,PHcollar,Hseq
     ENDIF
  END DO
  RETURN
END SUBROUTINE WriteLog
!****************************************************************************************
!> generates soil grid output, format is equal to nodes.in
SUBROUTINE OutFEM(t,kout)
  USE GridData
  USE SolData
  USE WatFun
  USE RhizoData
  IMPLICIT NONE
  
  REAL(sp),INTENT(in)::t
  INTEGER(ap), INTENT(in)::kout
  INTEGER(ap):: corner(8),ie,m,ic,i
  CHARACTER file*14
  REAL(sp):: sum,sumc,sume,sumce,tht_(nPt),solt,tottht,totsol,tht
  !> \param t current time
  !> \param kout current number for the FEM output
  
  WRITE (file,'(A11)')'out/outfem.'
  IF (kout.LT.10) THEN
     WRITE (file(12:12),'(I1)') kout
  ELSEIF (kout.LT.100) THEN
     WRITE (file(12:13),'(I2)') kout
  ELSE
     WRITE (file(12:14),'(I3)') kout	
  ENDIF
  OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
  ! calculate total water and solute volume within domain:
  sum=0.0_dp
  sumC=0.0_dp
  DO 1 iE=1,nElm!,2
     ! assign cuboid corner nodes:
     corner(1)=elmnod(1,iE)
     corner(2)=elmnod(2,iE)
     corner(3)=elmnod(3,iE)
     corner(4)=elmnod(4,iE)
     corner(5)=elmnod(5,iE)
     corner(6)=elmnod(6,iE)
     corner(7)=elmnod(7,iE)
     corner(8)=elmnod(8,iE)
     ! average cuboid water and solute content:
     sumE=0.0_sp
     sumCE=0.0_sp
     DO 11 ic=1,8
        i=corner(ic)
        M=MatNum(i)
        IF (soiltab) THEN
           tht=Fth_soiltab(hNew(i),M)
           solt=(Fth_soiltab(hNew(i),M)+ChPar(1,M)*ChPar(5,M))*Conc(i)
        ELSE
           tht=Fth(hNew(i),par(:,M),i)
           solt=(Fth(hNew(i),par(:,M),i)+ChPar(1,M)*ChPar(5,M))*Conc(i)
        ENDIF
        sumE=sumE+tht
        sumCe=sumCE+solt
11   ENDDO
     sum=sum+sumE
     sumC=sumC+sumCE
1 ENDDO
  tottht= sum*dxGrid*dyGrid*dzGrid/8._sp
  totsol=sumC*dxGrid*dyGrid*dzGrid/8._sp
  WRITE(8,'(/''Total Water Volume at Time  '',F12.4,1X,A5,'' is'',1pE15.8,     ''.'')')&
       t,TmUnit,tottht
  WRITE(8,'(/''Total Solute Volume at Time '',F12.4,1X,A5,'' is'',1pE15.8,     ''.'')')&
       t,TmUnit,totsol
  WRITE (8,'(/''Length Unit is '',A5)')LnUnit
  IF (.NOT. lRhizo) THEN
      WRITE (8,'(/''Node# Mater.#'',4X,''x'',9X,''y'',9X,''z'',11X,''h'',9X,''conc.'',7X,''theta'',7X,''wsink'',7X,''csink'')')
      IF (soiltab) THEN
         DO  i=1,nPt
            M=MatNum(i)
            tht_(i)=Fth_soiltab(hNew(i),M)
            WRITE (8,'(I7,6X,I2,3(1X,1pE11.4),5(1X,1pE11.4))')i,M,xGrid(i),yGrid(i),zGrid(i),hNew(i),Conc(i),tht_(i),sink(i),csink(i)
         ENDDO
      ELSE
         DO  i=1,nPt
            M=MatNum(i)
            tht_(i)=Fth(hNew(i),par(:,M),i)
            WRITE (8,'(I7,6X,I2,3(1X,1pE11.4),5(1X,1pE11.4))')i,M,xGrid(i),yGrid(i),zGrid(i),hNew(i),Conc(i),tht_(i),sink(i),csink(i)
         ENDDO
      ENDIF
  ELSE
      WRITE (8,'(/''Node# Mater.#'',4X,''x'',9X,''y'',9X,''z'',11X,''h'',9X,''conc.'',7X,''theta'',7X,''wsink'',7X,''csink'',7X,''thetaTot'',7X,''cond'')')
      DO i=1,nPt
         WRITE (8,'(I7,6X,I2,3(1X,1pE11.4),7(1X,1pE11.4))')i,M,xGrid(i),yGrid(i),zGrid(i),hNew(i),Conc(i),theta(i),sink(i),csink(i),&
                                                           thetaTot(i), con(i)
      ENDDO
  ENDIF
  CLOSE (8)
  RETURN
END SUBROUTINE OutFEM
!****************************************************************************************
!> in case of a continous domain, saves transroot for the basic root nodes -> draw the root in matlab
SUBROUTINE OutTrans(t,kout)
  USE typedef
  USE DoussanMat, ONLY: transroot,transtip,nplant,nsub
  USE RootData, ONLY: nrec,ngrow
  IMPLICIT NONE
  
  REAL(sp),INTENT(in)::t
  INTEGER(ap),INTENT(in)::kout
  INTEGER(ap) :: i,ipl
  CHARACTER file*17    
  !> \param t current time
  !> \param kout current number for the FEM output
  
  ipl = 1 !dummy
  WRITE (file,'(A14)')'out/transroot.'
  IF (kout.LT.10) THEN
     WRITE (file(15:15),'(I1)') kout
  ELSEIF (kout.LT.100) THEN
     WRITE (file(15:16),'(I2)') kout
  ELSE
     WRITE (file(15:17),'(I3)') kout	
  ENDIF
  OPEN (UNIT=444,FILE='out/Transroot.out',STATUS='UNKNOWN')
  WRITE (444,'(/''Basic root nodes translations for each plant.'')')
  DO ipl=1,nplant
     WRITE (444,'(/'' '')')
     WRITE (444,'(/''Plant '',i1)') ipl
     WRITE (444,'(/'' X  Y (times)'')')
     DO i=1,nrec
        WRITE (444,'(i3, i3)') transroot(i,1,nsub(i,ipl),ipl), transroot(i,2,nsub(i,ipl),ipl)
     END DO
     WRITE (444,'(/'' '')')
     WRITE (444,'(/'' '')')
     DO i=1,ngrow
        WRITE (444,'(i3, i3)') transtip(i,1,nsub(i,ipl),ipl), transtip(i,2,nsub(i,ipl),ipl)
     END DO
  END DO
  CLOSE(444)
  RETURN
END SUBROUTINE OutTrans
!*************************************************************************
!> if appropriate, writes depths profiles 
SUBROUTINE Zprofiles(t)
  USE GridData
  USE SolData
  USE WatFun
  IMPLICIT NONE
  REAL(sp),INTENT(in) :: t
  INTEGER(ap):: iz,M,prof,ixy,i
  CHARACTER form*5
  REAL(dp):: thz(1000),phz(1000),sz(1000)
  !> \param t current time
  
  WRITE(form,'(I3)')(nx*ny)
  ! open files
  OPEN (UNIT=121,FILE='out/ProfileTH.out',STATUS='OLD',POSITION='APPEND')
  OPEN (UNIT=122,FILE='out/ProfilePH.out',STATUS='OLD',POSITION='APPEND')
  OPEN (UNIT=123,FILE='out/ProfileS.out',STATUS='OLD',POSITION='APPEND')
  ! calculate average/sum
  prof=0
  DO iz=1,nz
     thz(iz)=0.
     phz(iz)=0.
     sz(iz)=0.
     DO ixy=1,nx*ny
        i=ixy+prof
        M=MatNum(i)
        IF (soiltab) THEN
           thz(iz)=thz(iz)+Fth_soiltab(hNew(i),M)
        ELSE
           thz(iz)=thz(iz)+Fth(hNew(i),par(:,M),i)
        ENDIF
        phz(iz)=phz(iz)+hNew(i)
        sz(iz)=sz(iz)+sink(i)
     ENDDO
     thz(iz)=thz(iz)/(nx*ny)
     phz(iz)=phz(iz)/(nx*ny)
     prof=prof+nx*ny
  ENDDO
  WRITE (121,'(1X,F15.6)') t,(thz(i),i=1,nz)
  WRITE (122,'(1X,F15.6)') t,(phz(i),i=1,nz)
  WRITE (123,'(1X,F15.6)') t,(sz(i),i=1,nz)
  CLOSE(121)
  CLOSE(122)
  CLOSE(123)
END SUBROUTINE Zprofiles
!*************************************************************************
!> output for observation probes
SUBROUTINE OutObsProbe(t)
  USE GridData
  USE SolData
  USE WatFun
  USE ObsData
  IMPLICIT NONE
  INTEGER(ap)::ip,n,M,minl,i
  REAL(sp),INTENT(in) :: t
  CHARACTER file*13,form*5
  REAL(sp):: tht(1000),pht(1000),AverPH,AverTH
  !> \param t current time
  
  DO ip=1,nPr
     ! calculations, theta and PH
     tht=0
     pht=0
     DO n=1,nodebyPr(ip)
        M=MatNum(NodePr(ip,n))
        IF (soiltab) THEN
           tht(n)=Fth_soiltab(hNew(NodePr(ip,n)),M)
        ELSE
           tht(n)=Fth(hNew(NodePr(ip,n)),par(:,M),NodePr(ip,n))
        ENDIF
        pht(n)=hNew(NodePr(ip,n))
     ENDDO
     
     ! define file names
     WRITE (file,'(A13)')'out/Probe.   '
     IF (ip.LT.10) THEN
        WRITE (file(11:11),'(I1)') ip
     ELSE
        WRITE (file(11:12),'(I2)') ip
     ENDIF
     OPEN (UNIT=10,FILE=file,STATUS='OLD',POSITION='APPEND')
     IF (DistrP(ip).EQ.1) THEN
        !> average of each plane
        AverTH=SUM(tht)/nodebyPr(ip)
        AverPH=SUM(pht)/nodebyPr(ip)
        IF (VarP(ip)==1) WRITE (10,'(2(1X,F12.4))') t,AverPH
        IF (VarP(ip)==2) WRITE (10,'(2(1X,F12.4))') t,AverTH
        IF (VarP(ip)==3) WRITE (10,'(3(1X,F12.4))') t,AverPH,AverTH
        CLOSE(10)
     ELSE
     ! complete distribution asked
     minl=MIN(nodebyPr(ip),500)
     WRITE(form,'(I3)')(minl+1)
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
        OPEN(UNIT=11,FILE=file,STATUS='OLD',POSITION='APPEND')
        WRITE(form,'(I3)')(nodebyPr(ip)-500)+1
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
!> in case of Somma root growth, writes root output; format is the same as RootSys 
SUBROUTINE OutRoo(t)
  USE typedef
  USE ParamData, ONLY : pi
  USE RootData
  USE PlntData, ONLY : LA
  USE DoussanMat, ONLY : nplant
  
  IMPLICIT NONE
  REAL(sp) :: t
  INTEGER(ap) :: irec,igrow,ifive,iestbl,i,linlim
  CHARACTER outfile*12
  !> \param t current time
  
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
  DO irec=1,nrec
     WRITE (8,'(I5,3(1X,1pE10.3),1X,I5,1X,I2,1X,I5,2(1X,1pE14.8),1X,1pE9.3)')&		
          irec,xs(irec),ys(irec),zs(irec),irecpr(irec),ordseg(irec),ibrseg(irec),&
          seglen(irec),segsur(irec),segmas(irec)
     WRITE (8,'(1pE11.4)') timorg(irec)
  END DO
  WRITE (8,*)
  WRITE (8,'(''Total # of growing branch tips:'')')
  WRITE (8,*) ngrow
  WRITE (8,*)
  WRITE (8,'(''tipID#'',4X,''xg'',10X,''yg'',10X,''zg'',6X,''sg.bhd.tp. '',''ord  br#  tot.br.lgth. axs#'')')
  WRITE (8,'(''overlength'',2X,''# of estblished points'')')
  
  WRITE (8,'(''time of establishing (-->)'')')
  ! write list of all growing tips:
  DO igrow=1,ngrow
     WRITE(8,'(I5,3(1X,1pE11.4),1X,I5,6X,I2,1X,I5,1X,1pE11.4,3X,I3)')&
          igrow,xg(igrow),yg(igrow),zg(igrow),irecsg(igrow),ordgrw(igrow),&
          ibrgrw(igrow),brlgth(igrow),iaxis(igrow)
     WRITE (8,'(1pE11.4,1X,I5)')ovrtime(igrow),nestbl(igrow)
     IF (nestbl(igrow).GT.0) THEN
        ifive=0
51      linlim=MIN(ifive+5,nestbl(igrow))
        WRITE (8,'(5(1X,1pE13.6))')(timest(igrow,iestbl),iestbl=ifive+1,linlim)
        ifive=ifive+5
        IF (linlim.LT.nestbl(igrow)) GOTO 51
     ENDIF
  END DO
  CLOSE(8)
  RETURN
END SUBROUTINE OutRoo
!*********************************************************************
!> writes a list of all root segments and their corresponding radial and axial
!> fluxes, conductivities, etc.: OutRoot.XX 
SUBROUTINE OutDou(t,kout)
  USE RootData
  USE DoussanMat, ONLY: Lr,Khr,PHs,PHr,sinkR,axialRootFlow,Qi,nplant,Qd,Q_bc1,veloRoot
  USE GridData, ONLY: dxgrid,dygrid,dzgrid
  USE PlntData, ONLY:Tact
  USE ParamData, ONLY: pi
  USE SoluteRootMat, ONLY : segsorb,l_linSorb,l_freundSorb, segsolv
  IMPLICIT NONE
  REAL(sp) :: t
  INTEGER(ap) :: irec,kout,ipl
  CHARACTER file*15
  WRITE (file,'(A12)')'out/outRoot.'
  !> \param t current time
  !> \param kout current number for the root output
  
  DO ipl=1,nplant
     WRITE (file(11:11),'(I1)') ipl
     IF (kout.LT.10) THEN
        WRITE (file(13:13),'(I1)') kout
     ELSEIF (kout.LT.100) THEN
        WRITE (file(13:14),'(I2)') kout
     ELSE
        WRITE (file(13:15),'(I3)') kout
     ENDIF
     OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
     WRITE (8,'(''Time:'')')
     WRITE (8,*) t
     WRITE (8,*)
     WRITE (8,'(''Total # of segment records:'')')
     WRITE (8,*) nrec
     WRITE (8,*)
     WRITE (8,90)

    If(l_linSorb.or.l_freundSorb)then
    DO irec=1,nrec
    WRITE (8,'(I6,3(1X,1pE10.3),I5,1X,I6,1X,E9.3,1X,E9.3,1X,E13.5,1X,E13.5,E13.5,8(1X,E13.5),E13.5)') &
             irec,xs(irec)+xplant(ipl),ys(irec)+yplant(ipl),&
             zs(irec),ibrseg(irec),irecpr(irec),Lr(irec),Khr(irec),PHs(irec,ipl),PHr(irec+1,ipl),sinkR(irec,ipl),&
             axialRootFlow(irec,ipl),Qi(irec,ipl),Qd(irec,ipl),Q_bc1(irec,ipl),segrad(irec),veloRoot(irec,ipl),&
             segconc(irec),segsorb(irec)
     END DO
     CLOSE(8)
90   FORMAT(' segID          x          y         z     br#    prev     Lr     Khr      Phinter     PHxylem     radialRootFlow     axialRootFlow      Qi        Qd       Q_bc       radius       veloRoot           segconc   segsorb ')
	Else
DO irec=1,nrec
	
        WRITE (8,'(I6,3(1X,1pE10.3),I5,1X,I6,1X,E9.3,1X,E9.3,1X,E13.5,1X,E13.5,E13.5,8(1X,E13.5))') &
             irec,xs(irec)+xplant(ipl),ys(irec)+yplant(ipl),&
             zs(irec),ibrseg(irec),irecpr(irec),Lr(irec),Khr(irec),PHs(irec,ipl),PHr(irec+1,ipl),sinkR(irec,ipl),&
             axialRootFlow(irec,ipl),Qi(irec,ipl),Qd(irec,ipl),Q_bc1(irec,ipl),segrad(irec),veloRoot(irec,ipl),&
             segconc(irec) 
     END DO
     CLOSE(8)
91    FORMAT(' segID          x          y         z     br#    prev     Lr     Khr      Phinter     PHxylem     radialRootFlow     axialRootFlow      Qi        Qd       Q_bc       radius       veloRoot       segconc   ')
End if 
  END DO
  RETURN
END SUBROUTINE OutDou
!***************************************************************************
!> writes balance.out and remove.out 
SUBROUTINE SubReg(t,iCount)
  USE ParamData, ONLY: lChem,lretry,last_out
  USE GridData
  USE SolData
  USE WatFun
  USE CumData
  USE Doussanmat, ONLY : SinkR
  USE Rootdata, ONLY : lDou,lCou,lFed,lno_RWU,lno_Archi,ldJvL,tlim_dJvL,lPast,tPast,sinkredDPast
  USE tmctrl, ONLY : tOut
  IMPLICIT NONE
  REAL(sp),INTENT(in) :: t
  INTEGER(ap) :: i,j,k,l,Mi,Mj,Mk,Ml,iSE,iE,iCount,icut
  REAL(sp) :: cbalr=0,cc,cbalt
  INTEGER(sp) :: xE,yE,zE
  REAL(sp) ::cnewe,cel,wel,DeltC,time
  REAL(dp) ::WatVol,ConVol,ww,wBalT,DeltW,WNewE,VE
  REAL(sp), ALLOCATABLE,DIMENSION (:) :: tPast_old
  REAL(sp), ALLOCATABLE,DIMENSION (:,:) :: sinkredD,sinkredDPast_old
  REAL(dp), ALLOCATABLE,DIMENSION (:) :: SSF_old,RLD_old
  !> \param t current time
  !> \param iCount at the first timestep iCount = 0 
  
  IF (.not.allocated(watin)) ALLOCATE(WatIn(nElm),SolIn(nElm))
  IF (ldJvL) ALLOCATE (sinkredD(1:nElm,1))
  IF (icount.EQ.0) THEN
     VElm(1:nElm)=dxGrid*dyGrid*dzGrid
     IF (lCou.AND.lno_Archi.AND.(nexSSF.NE.nex.OR.neySSF.NE.ney.OR.nezSSF.NE.nez)) THEN
        IF (iCount.EQ.0) THEN
           ALLOCATE (SSF_old(1:nexSSF*neySSF*nezSSF))
           SSF_old=SSF
           DEALLOCATE (SSF)
           ALLOCATE (SSF(1:nElm))
           SSF=0._dp
           DO zE=1,nezSSF
              DO yE=1,neySSF
                 DO xE=1,nexSSF
                    SSF(INT(FLOOR(((REAL(xE)-1)*nex)/nexSSF)+1+FLOOR(((REAL(yE)-1)*ney)/neySSF)*nex+FLOOR(((REAL(zE)-1)*nez)/nezSSF)*nex*ney))=SSF(INT(FLOOR(((REAL(xE)-1)*nex)/nexSSF)+1+FLOOR(((REAL(yE)-1)*ney)/neySSF)*nex+FLOOR(((REAL(zE)-1)*nez)/nezSSF)*nex*ney))+SSF_old(xE+(yE-1)*nexSSF+(zE-1)*nexSSF*neySSF)
                 END DO
              END DO
           END DO
        ENDIF
     ELSEIF (lFed.AND.lno_Archi.AND.(nexRLD.NE.nex.OR.neyRLD.NE.ney.OR.nezRLD.NE.nez)) THEN
        IF (iCount.EQ.0) THEN
           ALLOCATE (RLD_old(1:nexRLD*neyRLD*nezRLD))
           RLD_old=RLD
           DEALLOCATE (RLD)
           ALLOCATE (RLD(1:nElm))
           RLD=0._dp
           DO zE=1,nezRLD
              DO yE=1,neyRLD
                 DO xE=1,nexRLD
                    RLD(INT(FLOOR(((REAL(xE)-1)*nex)/nexRLD)+1+FLOOR(((REAL(yE)-1)*ney)/neyRLD)*nex+FLOOR(((REAL(zE)-1)*nez)/nezRLD)*nex*ney))=(RLD(INT(FLOOR(((REAL(xE)-1)*nex)/nexRLD)+1+FLOOR(((REAL(yE)-1)*ney)/neyRLD)*nex+FLOOR(((REAL(zE)-1)*nez)/nezRLD)*nex*ney))*dxGrid*dyGrid*dzGrid+RLD_old(xE+(yE-1)*nexRLD+(zE-1)*nexRLD*neyRLD)*dxRLD*dyRLD*dzRLD)/dxGrid*dyGrid*dzGrid
                 END DO
              END DO
           END DO
        ENDIF
     ENDIF
  ENDIF
  
  ! write balance.out
  OPEN (UNIT=10,FILE='out/balance.out',STATUS='OLD',POSITION='APPEND')
  
  ! initializing some variables
  WatVol=0.0_dp
  DeltW=0.0
  ConVol=0.0
  DeltC=0.0
  
  DO iE=1,nElm ! loop over elements
     wEl=0.0
     cEl=0.0
     DO iSE=1,5  
        i=elmnod(iL(1,iSE,subN(iE)),iE) ! i,j,k,l are corner nodes of the tetrahedal element, ise counts 5 which is equal to five tedrahedrals which is equal to a cubic.
        j=elmnod(iL(2,iSE,subN(iE)),iE)
        k=elmnod(iL(3,iSE,subN(iE)),iE)
        l=elmnod(iL(4,iSE,subN(iE)),iE)
        Mi=MatNum(i)
        Mj=MatNum(j)
        Mk=MatNum(k)
        Ml=MatNum(l)
        VE=ABS(Deter(iSE,iE))/6.
        WNewE=VE*(theta(i)+theta(j)+theta(k)+theta(l))/4 !V -> water mass balance (elements)
        WatVol=WatVol+WNewE !sum of V
        wEl=wEl+WNewE
        IF (lChem) THEN
           CNewE=VE*((theta(i)+ChPar(1,Mi)*ChPar(5,Mi))*Conc(i)+(theta(j)+&
                ChPar(1,Mj)*ChPar(5,Mj))*Conc(j)+(theta(k)+ChPar(1,Mk)*ChPar(5,Mk))*Conc(k)+&
                (theta(l)+ChPar(1,Ml)*ChPar(5,Ml))*Conc(l))/4  ! solute mass balance (elements)
           ConVol=ConVol+CNewE !sum of solute mass
           cEl=cEl+CNewE
        ENDIF
        IF (iSE.EQ.5) THEN
           IF (iCount.EQ.0) THEN
              WatIn(iE)=wEl
              IF (lChem) SolIn(iE)=cEl
           ELSE
              DeltW=DeltW+ABS(WatIn(iE)-wEl)
              IF (lChem) DeltC=DeltC+ABS(SolIn(iE)-cEl)
           ENDIF
        ENDIF
     END DO
     IF (ldJvL) sinkredD(iE,1)=sink_cube(iE)	
  END DO
  IF (ldJvL) THEN
     IF (iCount.EQ.0) THEN
        lPast=3							!Number of times saved in tPast, sinkRedD, etc.
        ALLOCATE (tPast(1:lPast))
        tPast=(/0._sp,t-0.0002_sp,t-0.0001_sp/)		!Times saved in sinkRedD
        ALLOCATE (sinkredDPast(SIZE(sinkredD,1),1:lPast))	!Dimensions are (number of soil elements ; number of times saved)
        sinkredDPast(:,1)=sinkredD(:,1)			!Contains the past values of sink_cube for the latest times, 1st column for the oldest time step.
        sinkredDPast(:,2)=sinkredD(:,1)
        sinkredDPast(:,3)=sinkredD(:,1)
     ELSE
        ALLOCATE (tPast_old(1:lPast))
        tPast_old=tPast
        ALLOCATE (sinkredDPast_old(SIZE(sinkredD,1),1:lPast))
        sinkredDPast_old=sinkredDPast				!Keep previous information in mind
        icut=0
14      icut=icut+1						!tPast and sinkredD are going to be pruned, beginning by the oldest time steps. "icut-1" is the number of pruned time steps.
        IF ((tPast(icut)-t).LT.(-tlim_dJvL)) THEN
           lPast=lPast-1
           IF (lPast.GT.0) THEN
              GOTO 14
           ELSE
              lPast=lPast+1
           ENDIF
        ELSE
           lPast=lPast+1
        ENDIF
        DEALLOCATE (tPast,sinkredDPast)
        ALLOCATE (tPast(1:lPast))
        ALLOCATE (sinkredDPast(SIZE(sinkredD,1),1:lPast))
        IF (lPast.GT.1) THEN
           tPast=(/tPast_old(icut:lPast+icut-2),t/)		!Here we already use the "new" value of lPast.
           sinkredDPast(:,1:lPast-1)=sinkredDPast_old(:,icut:lPast+icut-2)
           sinkredDPast(:,lPast)=sinkredD(:,1)
        ELSE
           tPast=t
           sinkredDPast(:,1)=sinkredD(:,1)
        ENDIF
     ENDIF
  ENDIF
  
  ! Mass balance calculation
  IF (iCount.EQ.0) THEN
     IF (.NOT.lretry) THEN
        wVolI=WatVol
        IF (lChem) cVolI=ConVol
        IF (lChem) THEN 
           WRITE(10,130)
        ELSEIF (lno_RWU) THEN
           WRITE(10,135)
        ELSEIF (lDou) THEN
           WRITE(10,136)
        ELSE
           WRITE(10,137)
        ENDIF
        IF (lChem) THEN
           WRITE(10,140) t,WatVol,ConVol
        ELSE
           WRITE(10,150) t,WatVol
        ENDIF
     ELSE
        CLOSE (10)
        OPEN (Unit=20,FILE='out/balance.out',STATUS='OLD',ERR=10)
        READ (20,*)
        READ (20,*)
        READ (20,*)
        IF (.NOT.lChem) THEN
           READ (20,*) time,wVolI
        ELSE
           READ (20,*) time,wVolI,cVolI
        ENDIF
3       READ (20,*) time,WatVol,wBalT,wBalR,ww,wCumT,WatVol,wCumA,deltW
        IF (time.NE.tOut(last_out)) GOTO 3
        CLOSE (20)
        ! write remove.out
        OPEN (Unit=20,FILE='out/remove.out',STATUS='OLD',ERR=20)
        READ (20,*)
        READ (20,*)
        READ (20,*)
        READ (20,*)
4       READ (20,*) time,CumCh0,CumCh1,CumChR,ChemS(1),ChemS(2),CumRt,CumQ(1),CumQ(2)
        IF (time.NE.tOut(last_out)) GOTO 4
        CLOSE (20)
     ENDIF
  ELSE
     wBalT=WatVol-wVolI+wCumT
     ww=MAX(DeltW,wCumA)
     IF (ww.GE.1.e-25_dp) wBalR=ABS(wBalT)/ww*100  
     WRITE(*,*)'wBalR=',wBalR
     IF (lChem) THEN
        cBalT=ConVol-cVolI+cCumT
        cc=MAX(DeltC,cCumA)
        IF (cc.GE.1.e-25_dp) cBalR=ABS(cBalT)/cc*100
        WRITE(*,*)'cBalR=',cBalR
        WRITE(10,160) t,WatVol,wBalT,wBalR,ConVol,cBalT,cBalR,Peclet,Courant,wcumT,wcumA,deltW,deltc 
     ELSEIF (lno_RWU) THEN
        WRITE(10,171) t,WatVol,wBalT,wBalR,peclet,courant,wcumT,WatVol-wvolI,wcumA,deltW
     ELSEIF (lDou) THEN
        WRITE(10,170) t,WatVol,wBalT,wBalR,ww,wcumT,WatVol-wvolI,wcumA,deltW,RootSk,SUM(SinkR)
     ELSE
        WRITE(10,172) t,WatVol,wBalT,wBalR,ww,wcumT,WatVol-wvolI,wcumA,deltW,RootSk
     ENDIF
  ENDIF
  CLOSE (10)
  WatVolOld=WatVol
130 FORMAT(/,'    Time [T]   WatVol [V]  WatBalT [V]  WatBalR [%]  CncVol [M]   ',&
         'CncBalT [M]  CncBalR [%]  Peclet       Courant      ',&
         'wCumT [V]    wCumA[V]   DeltaW[V]   DeltaC[V]')
135 FORMAT(/,'    Time [T]   WatVol [V]  WatBalT [V]  WatBalR [%] ' ,&
         'Peclet    Courant',&
         'wCumT [V]  WatVol-wVolI[V]   wCumA[V]    Deltaw[V] ')
136 FORMAT(/,'    Time [T]   WatVol [V]  WatBalT [V]  WatBalR [%]  ww [V]  ',&
         'wCumT [V]  WatVol-wVolI[V]   wCumA[V]    Deltaw[V]   TotalRootFlow [V/T]',&
         '  TotRadialFlow [V/T]')
137 FORMAT(/,'    Time [T]   WatVol [V]  WatBalT [V]  WatBalR [%]  ww [V]  ',&
         'wCumT [V]  WatVol-wVolI[V]   wCumA[V]    Deltaw[V]   TotalRootFlow [V/T]')
140 FORMAT(f12.4,1X,1pE12.5,27x,1pE12.5)
150 FORMAT(f12.4,1X,1pE12.5)
160 FORMAT(f12.4,1X,14(1pE12.5,1X))
170 FORMAT(f12.4,1X,10(1pE12.5,1X))
171 FORMAT(f12.4,1X,9(1pE12.5,1X))
172 FORMAT(f12.4,1X,9(1pE12.5,1X))
  OPEN (UNIT=10,FILE='out/remove.out',STATUS='OLD',POSITION='APPEND')
  IF (iCount.EQ.0) THEN
     IF (.NOT.lretry) THEN
        WRITE(10,230)
     ENDIF
  ELSE
     WRITE(10,240)  t,CumCh0,CumCh1,CumChR,ChemS(1),ChemS(2),CumRt ,CumQ(1) ,CumQ(2)
  ENDIF
230 FORMAT(/,'    Time [T]   CumCh0 [M]   CumCh1 [M]   CumChR [M] ChemS(1) [M]',&
         ' ChemS(2) [M]    CumRt [V]  CumQ(1) [V]  CumQ(2) [V]')
240 FORMAT(f12.4,1X,8(1pE12.5,1X))
  CLOSE(10)
  RETURN
10 STOP '< out/balance.out > not found  -- program terminated.'
20 STOP '< out/remove.out > not found  -- program terminated.'
END SUBROUTINE SubReg
!************************************************************************
!> writes veloci.XX, snkElm.XX, and betElm.XX 
SUBROUTINE FlxOut(kOut)
  USE typedef
  USE ParamData, ONLY: lOutPartrace,pi
  USE Soldata
  USE GridData
  USE DoussanMat, ONLY :betac_cube2
  USE PlntData, ONLY:TotSur
  USE RootData, ONLY:lDou,lCou,ldJvL,lFed
  USE Watfun, ONLY: Fh_from_Th
  IMPLICIT NONE
  
  INTEGER(ap), INTENT(in)::kout
  INTEGER(ap):: i
  CHARACTER file*14,  file2*14,  file3*14
  !> \param kOut current number for the FEM output  
  
  WRITE (file,'(A11)')'out/veloci.'
  WRITE (file2,'(A11)')'out/snkElm.'
  WRITE (file3,'(A11)')'out/betElm.'
  
  IF (kOut.LT.10) THEN
     WRITE (file(12:12),'(I1)') kOut
     WRITE (file2(12:12),'(I1)') kOut
     WRITE (file3(12:12),'(I1)') kOut
  ELSEIF (kout.LT.100) THEN
     WRITE (file(12:13),'(I2)') kOut
     WRITE (file2(12:13),'(I2)') kOut
     WRITE (file3(12:13),'(I2)') kOut
  ELSE
     WRITE (file(12:14),'(I3)') kOut
     WRITE (file2(12:14),'(I3)') kOut
     WRITE (file3(12:14),'(I3)') kOut	
  ENDIF

  OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
  WRITE (8,'(/''Node#'',9X,''x'',9X,''y'',9X,''z'',11X,''Vx'',8X,''Vy'',8X,''Vz'')')
  CALL Veloc
  
  ! write nodal values:
  DO i=1,nPt
     WRITE (8,'(I5,8X,3(1X,1pE11.4),3(1X,1pE11.4))')i,xGrid(i),yGrid(i),zGrid(i),Vx(i),Vy(i),Vz(i)
  END DO
  CLOSE(8)
  
  if (lOutPartrace) then
     OPEN (UNIT=50,FILE=file3,STATUS='UNKNOWN')
     WRITE (50,'(/''Element Nr.'',9X,''BetaElm'')' )
     DO i=1,nElm
        WRITE (50,*) i,  betac_cube2 (i)*TotSur
     END DO
     CLOSE(50)
  end if
  if (lOutPartrace.OR.lDou) then
     OPEN (UNIT=90,FILE=file2,STATUS='UNKNOWN')
     WRITE (90,'(/''Element Nr.'',9X,''SinkElm'')' )
     DO i=1,nElm
        WRITE (90,*) i,Sink_cube(i)
     END DO
     CLOSE(90)
  ELSEIF ((lCou.AND..NOT.ldJvL).OR.lFed) THEN
     OPEN (UNIT=90,FILE=file2,STATUS='UNKNOWN')
     WRITE (90,'(/''Element Nr.'',9X,''SinkElm'',9X,''HElm'')' )
     DO i=1,nElm
        WRITE (90,*) i,Sink_cube(i),HElm(i)
     END DO
     CLOSE(90)
  ELSEIF (lCou.AND.ldJvL) THEN
     OPEN (UNIT=90,FILE=file2,STATUS='UNKNOWN')
     WRITE (90,'(/''Element Nr.'',9X,''SinkElm'',9X,''Hint'',9X,''Hbulk'')' )
     DO i=1,nElm
        WRITE (90,*) i,Sink_cube(i),HElm(i),Fh_from_Th(SUM(theta(elmnod(1:8,i)))/8.0_sp,par(:,MatNum(elmnod(1,i))))+(zgrid(elmnod(1,i))+zgrid(elmnod(8,i)))/2.0_sp
     END DO
     CLOSE(90)
  end if
  RETURN
END SUBROUTINE FlxOut
!****************************************************************************
!> writes .vtk output for the soil grid; veloci.XX 
SUBROUTINE OutVTK(kOut)
  USE typedef
  USE Soldata
  USE GridData
  USE ParamData, ONLY: lchem
  USE RootData, ONLY: lSomma_growth
  USE StrData, ONLY: s
  USE RhizoData
  IMPLICIT NONE
  
  INTEGER(ap), INTENT(in)::kout
  INTEGER(ap):: i, k, j, ix, iy, iz
  CHARACTER file*22
  !> \param kOut current number for the FEM output
  
  WRITE (file,'(A15)')'out/vtk/veloci.'
  
  IF (kOut.LT.10) THEN
     WRITE (file(16:16),'(I1)') kOut
     WRITE (file(17:20),'(A4)') '.vtk'
  ELSEIF (kOut.LT.100) THEN
     WRITE (file(16:17),'(I2)') kOut
     WRITE (file(18:21),'(A4)') '.vtk'
  ELSE
     WRITE (file(16:18),'(I3)') kOut
     WRITE (file(19:22),'(A4)') '.vtk'
  ENDIF
  
  PRINT *, file
  
  OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
  WRITE (8,'(A26)')'# vtk DataFile Version 3.0' 
  WRITE (8,'(A12)')'model R-SWMS'
  WRITE (8,'(A5)')'ASCII'
  WRITE (8,'(A24)')'DATASET RECTILINEAR_GRID'
  
  IF (continu) THEN
     
     WRITE (8,'(A11,1X,3I6)')'DIMENSIONS ', nx+1, ny+1, nz
     
     WRITE (8,'(A14,1X,I6,1X,A5)')'X_COORDINATES ',nx+1, 'float' 
     DO i=1,nx
        WRITE (8,'(1X,1pE11.4)',advance="no") xgrid(i)
     END DO
     WRITE (8,'(1X,1pE11.4)',advance="no") xgrid(SIZE(xgrid))+dxgrid
     
     WRITE (8,'(/A14,1X,I6,1X,A5)')'Y_COORDINATES ',ny+1, 'float'
     DO i=1,ny*nx,nx
        WRITE (8,'(1pE11.4)',advance="no") ygrid(i)
     END DO
     WRITE (8,'(1X,1pE11.4)',advance="no") ygrid(SIZE(ygrid))+dygrid
     
     WRITE (8,'(/A14,1X,I6,1X,A5)')'Z_COORDINATES ',nz, 'float' 
     DO i=1,SIZE(zgrid),nx*ny
        WRITE (8,'(1pE11.4)',advance="no") zgrid(i)
     END DO
     
     WRITE (8,'(/A10,1X,I7)')'POINT_DATA', (nx+1)*(ny+1)*nz
     WRITE (8,'(A24)')'SCALARS velocity float 3'
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     
     CALL Veloc
     
     k=0
     j=0
     DO  iz=1,nz
        DO  iy=1,ny+1
           DO  ix=1,nx+1     
              k = ix + (iy-1)*nx + (iz-1)*nx*ny
              j = j +1
              IF (ix.EQ.nx+1) THEN
                 k = k - nx
              END IF
              IF  (iy.EQ.ny+1) THEN
                 k = k - (nx*ny)
              END IF
              WRITE (8,'(3(1X,1pE11.4))')Vx(k),Vy(k),Vz(k)
           END DO
        END DO
     END DO
     
     WRITE (8,'(A27)')'SCALARS pressurehead float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     k=0
     j=0
     DO  iz=1,nz
        DO  iy=1,ny+1
           DO  ix=1,nx+1     
              k = ix + (iy-1)*nx + (iz-1)*nx*ny
              j = j +1
              IF (ix.EQ.nx+1) THEN
                 k = k - nx
              END IF
              IF  (iy.EQ.ny+1) THEN
                 k = k - (nx*ny)
              END IF
              WRITE (8,'(1X,1pE11.4)',advance="no")hnew(k)
           END DO
        END DO
     END DO
     WRITE (8,*)''
     WRITE (8,'(A17)')'SCALARS wc float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     k=0
     j=0
     DO  iz=1,nz
        DO  iy=1,ny+1
           DO  ix=1,nx+1     
              k = ix + (iy-1)*nx + (iz-1)*nx*ny
              j = j +1
              IF (ix.EQ.nx+1) THEN
                 k = k - nx
              END IF
              IF  (iy.EQ.ny+1) THEN
                 k = k - (nx*ny)
              END IF
              WRITE (8,'(1X,1pE11.4)',advance="no")theta(k)
           END DO
        END DO
     END DO
     WRITE (8,*)''
     WRITE (8,'(A17)')'SCALARS mat integer '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     k=0
     j=0
     DO  iz=1,nz
        DO  iy=1,ny+1
           DO  ix=1,nx+1     
              k = ix + (iy-1)*nx + (iz-1)*nx*ny
              j = j +1
              IF (ix.EQ.nx+1) THEN
                 k = k - nx
              END IF
              IF  (iy.EQ.ny+1) THEN
                 k = k - (nx*ny)
              END IF
              WRITE (8,'(1X,I3)',advance="no") MatNum(k)
           END DO
        END DO
     END DO
     WRITE (8,*)''
     WRITE (8,'(A17)')'SCALARS Kode integer '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     k=0
     j=0
     DO  iz=1,nz
        DO  iy=1,ny+1
           DO  ix=1,nx+1     
              k = ix + (iy-1)*nx + (iz-1)*nx*ny
              j = j +1
              IF (ix.EQ.nx+1) THEN
                 k = k - nx
              END IF
              IF  (iy.EQ.ny+1) THEN
                 k = k - (nx*ny)
              END IF
              WRITE (8,'(1X,I3)',advance="no") Kode(k)
           END DO
        END DO
     END DO
     IF(lchem) THEN
        WRITE (8,*)''
        WRITE (8,'(A19)')'SCALARS conc float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        k=0
        j=0
        DO  iz=1,nz
           DO  iy=1,ny+1
              DO  ix=1,nx+1     
                 k = ix + (iy-1)*nx + (iz-1)*nx*ny
                 j = j +1
                 IF (ix.EQ.nx+1) THEN
                    k = k - nx
                 END IF
                 IF  (iy.EQ.ny+1) THEN
                    k = k - (nx*ny)
                 END IF
                 WRITE (8,'(1X,1pE11.4)',advance="no")conc(k)
              END DO
           END DO
        END DO
     END IF
     WRITE (8,*)''
     IF (lRhizo) THEN
        WRITE (8,'(A20)')'SCALARS thtTot float'
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        k=0
        j=0
        DO iz=1,nz
            DO iy=1,ny+1
                DO ix=1,nx+1
                    k = ix+ (iy-1)*nx +(iz-1)*nx*ny
                    j = j + 1
                    IF (ix .EQ. nx+1) THEN
                        k = k -nx
                    ENDIF
                    IF (iy .EQ. ny + 1) THEN
                        k = k - (nx*ny)
                    ENDIF
                    WRITE (8, '(1X,1pE11.4)', advance="no") thetaTot(k)
                ENDDO
            ENDDO
        ENDDO
    ENDIF
     IF (lRhizo) THEN
        WRITE (8,'(A19)')'SCALARS Rnorm float'
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        k=0
        j=0
        DO iz=1,nz
            DO iy=1,ny+1
                DO ix=1,nx+1
                    k = ix+ (iy-1)*nx +(iz-1)*nx*ny
                    j = j + 1
                    IF (ix .EQ. nx+1) THEN
                        k = k -nx
                    ENDIF
                    IF (iy .EQ. ny + 1) THEN
                        k = k - (nx*ny)
                    ENDIF
                    WRITE (8, '(1X,1pE11.4)', advance="no") Rnorm(k)
                ENDDO
            ENDDO
        ENDDO
    ENDIF
     IF (lRhizo) THEN
        WRITE (8,'(A18)')'SCALARS cTot float'
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        k=0
        j=0
        DO iz=1,nz
            DO iy=1,ny+1
                DO ix=1,nx+1
                    k = ix+ (iy-1)*nx +(iz-1)*nx*ny
                    j = j + 1
                    IF (ix .EQ. nx+1) THEN
                        k = k -nx
                    ENDIF
                    IF (iy .EQ. ny + 1) THEN
                        k = k - (nx*ny)
                    ENDIF
                    WRITE (8, '(1X,1pE11.4)', advance="no") cTot_r(k)
                ENDDO
            ENDDO
        ENDDO
    ENDIF
    
     IF (lRhizo) THEN
        WRITE (8,'(A18)')'SCALARS cond float'
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        k=0
        j=0
        DO iz=1,nz
            DO iy=1,ny+1
                DO ix=1,nx+1
                    k = ix+ (iy-1)*nx +(iz-1)*nx*ny
                    j = j + 1
                    IF (ix .EQ. nx+1) THEN
                        k = k -nx
                    ENDIF
                    IF (iy .EQ. ny + 1) THEN
                        k = k - (nx*ny)
                    ENDIF
                    WRITE (8, '(1X,1pE11.4)', advance="no") con(k)
                ENDDO
            ENDDO
        ENDDO
      ENDIF


!!$           WRITE (8,'(A19)')'SCALARS sink float '
!!$           WRITE (8,'(A20)')'LOOKUP_TABLE default'
!!$           k=0
!!$           j=0
!!$           DO  iz=1,nz
!!$              DO  iy=1,ny+1
!!$                 DO  ix=1,nx+1     
!!$                    k = ix + (iy-1)*nx + (iz-1)*nx*ny
!!$                    j = j +1
!!$                    IF (ix.EQ.nx+1) THEN
!!$                       k = k - nx
!!$                    END IF
!!$                    IF  (iy.EQ.ny+1) THEN
!!$                       k = k - (nx*ny)
!!$                    END IF
!!$                    WRITE (8,'(1X,1pE11.4)',advance="no")sink(k)
!!$                 END DO
!!$              END DO
!!$           END DO
     IF(lSomma_growth) THEN
        WRITE (8,'(A27)')'SCALARS SoilStrength float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        k=0
        j=0
        DO  iz=1,nz
           DO  iy=1,ny+1
              DO  ix=1,nx+1     
                 k = ix + (iy-1)*nx + (iz-1)*nx*ny
                 j = j +1
                 IF (ix.EQ.nx+1) THEN
                    k = k - nx
                 END IF
                 IF  (iy.EQ.ny+1) THEN
                    k = k - (nx*ny)
                 END IF
                 WRITE (8,'(1X,1pE11.4)',advance="no") s(k)
              END DO
           END DO
        END DO
     END IF
     IF(lchem) THEN
        WRITE (8,*)''
        WRITE (8,'(A20)')'SCALARS csink float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        k=0
        j=0
        DO  iz=1,nz
           DO  iy=1,ny+1
              DO  ix=1,nx+1     
                 k = ix + (iy-1)*nx + (iz-1)*nx*ny
                 j = j +1
                 IF (ix.EQ.nx+1) THEN
                    k = k - nx
                 END IF
                 IF  (iy.EQ.ny+1) THEN
                    k = k - (nx*ny)
                 END IF
                 WRITE (8,'(1X,1pE11.4)',advance="no")csink(k)
              END DO
           END DO
        END DO
     END IF         
     
  ELSE
     
     WRITE (8,'(A11,1X,3I6)')'DIMENSIONS ', nx, ny, nz
     
     WRITE (8,'(A14,1X,I6,1X,A5)')'X_COORDINATES ',nx, 'float' 
     DO i=1,nx
        WRITE (8,'(1X,1pE11.4)',advance="no") xgrid(i)
     END DO
     
     WRITE (8,'(/A14,1X,I6,1X,A5)')'Y_COORDINATES ',ny, 'float'
     DO i=1,ny*nx,nx
        WRITE (8,'(1pE11.4)',advance="no") ygrid(i)
     END DO
     
     WRITE (8,'(/A14,1X,I6,1X,A5)')'Z_COORDINATES ',nz, 'float' 
     DO i=1,SIZE(zgrid),nx*ny
        WRITE (8,'(1pE11.4)',advance="no") zgrid(i)
     END DO
     
     WRITE (8,'(/A15,1X,I7)')'POINT_DATA', nPt
     WRITE (8,'(A24)')'SCALARS velocity float 3'
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     
     CALL Veloc
     DO i=1,nPt
        WRITE (8,'(3(1X,1pE11.4))')Vx(i),Vy(i),Vz(i)
     END DO
     
     WRITE (8,'(A27)')'SCALARS pressurehead float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE11.4)',advance="no") hnew(1:nPt)
     
     WRITE (8,'(A17)')'SCALARS wc float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE11.4)',advance="no") theta(1:nPt)
     
     WRITE (8,'(A19)')'SCALARS conc  float'
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE11.4)',advance="no") conc(1:nPt)
     
     WRITE (8,'(A19)')'SCALARS mat integer'
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,I3)',advance="no") MatNum(1:nPt) 
     
!!$           WRITE (8,'(A19)')'SCALARS sink float '
!!$           WRITE (8,'(A20)')'LOOKUP_TABLE default'
!!$           WRITE (8,'(1X,1pE11.4)',advance="no") sink(1:nPt)
    IF (lRhizo) THEN
     WRITE (8, '(A20)')'SCALARS thtTot float'
     WRITE (8, '(A20)')'LOOKUP_TABLE default'
     WRITE (8, '(1X,1pE11.4)', advance="no") thetaTot(1:nPt)

     WRITE (8, '(A19)')'SCALARS Rnorm float'
     WRITE (8, '(A20)')'LOOKUP_TABLE default'
     WRITE (8, '(1X,1pE11.4)', advance="no") Rnorm(1:nPt)
     
     WRITE (8, '(A18)')'SCALARS cTot float'
     WRITE (8, '(A20)')'LOOKUP_TABLE default'
     WRITE (8, '(1X,1pE11.4)', advance="no") cTot_r(1:nPt)
   
     WRITE (8, '(A18)')'SCALARS cond float'
     WRITE (8, '(A20)')'LOOKUP_TABLE default'
     WRITE (8, '(1X,1pE11.4)', advance='no') con(1:nPt)
   ENDIF

     WRITE (8,'(A20)')'SCALARS csink float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE11.4)',advance="no") csink(1:nPt)
     
     WRITE (8,'(A20)')'SCALARS Kode integer '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,I3)',advance="no") Kode(1:nPt)
     
     IF(lSomma_growth) THEN
        WRITE (8,'(A27)')'SCALARS SoilStrength float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        WRITE (8,'(1X,1pE11.4)',advance="no") s(1:nPt)
     END IF
     
  ENDIF

  WRITE (8,*)''
  WRITE (8,'(/A9,1X,I7)')'CELL_DATA', nElm
  WRITE (8,'(A24)')'SCALARS sinkElm float'
  WRITE (8,'(A20)')'LOOKUP_TABLE default'
  WRITE (8,'(1X,1pE11.4)',advance="no") sink_cube(1:nElm)
  
  WRITE (8,'(A24)')'SCALARS csinkElm float'
  WRITE (8,'(A20)')'LOOKUP_TABLE default'
  WRITE (8,'(1X,1pE11.4)',advance="no") csink_cube(1:nElm)

  CLOSE (8)
  
  RETURN
END SUBROUTINE OutVTK
!************************************************************************
!> generates .vtk output for the root outRoot.XX.XXX 
SUBROUTINE OutDouVTK(kOut,t)
  USE typedef
  USE ParamData, ONLY: pi
  USE RootData
  USE DoussanMat, ONLY:  veloRoot,Lr,Khr,PHs,PHr,sinkR,axialRootFlow,Qi,nplant,Qd,Q_bc1,nsub,transroot,&
       transtip
  USE GridData, ONLY: dxgrid,dygrid,dzgrid,nx,ny,nz,continu
  USE PlntData, ONLY:Tact
  USE SoluteRootMat, ONLY: segsorb,l_linSorb,l_freundSorb
  IMPLICIT NONE
  
  INTEGER(ap):: j, sizen,  kk
  REAL(sp) ::  xr(nrec+ngrow), yr(nrec+ngrow),t
  INTEGER(ap) :: irec,kout,ipl,trans_i(2), trans_j(2)
  INTEGER(ap):: igrow,irecn
  CHARACTER file*29 
  !> \param t current simulation time
  !> \param kOut current number for the FEM output
  
  DO ipl=1,nplant
     
     WRITE (file,'(A18)')'out/vtk/outRootXX.'
     
     
     IF (ipl.LT.10) THEN
        WRITE (file(16:16),'(I1)') 0
        WRITE (file(17:17),'(I1)') ipl
     ELSE
        WRITE (file(16:17),'(I2)') ipl
     END IF
     IF (kOut.LT.10) THEN
        WRITE (file(19:19),'(I1)') kOut
        WRITE (file(20:24),'(A4)') '.vtk'
     ELSEIF (kOut.LT.100) THEN
        WRITE (file(19:20),'(I2)') kOut
        WRITE (file(21:25),'(A4)') '.vtk'
     ELSE
        WRITE (file(19:21),'(I3)') kOut
        WRITE (file(22:26),'(A4)') '.vtk'
     ENDIF
     
     PRINT *, file
     
     OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
     WRITE (8,'(A26)')'# vtk DataFile Version 3.0' 
     WRITE (8,'(A12)')'model R-SWMS'
     WRITE (8,'(A5)')'ASCII'
     WRITE (8,'(A16)')'DATASET POLYDATA'
     WRITE (8,'(A7,1X,I6,1X,A5)')'POINTS ', nrec+ngrow, 'float'
     
     
     IF (continu) THEN
        
        ! write out points root segments
        DO  irec=1,nrec
           xr(irec) = xs(irec)+xplant(ipl)+1.*transroot(irec,1,nsub(irec,ipl),ipl)*dxgrid*nx
           yr(irec) = ys(irec)+yplant(ipl)+1.*transroot(irec,2,nsub(irec,ipl),ipl)*dygrid*ny
           WRITE (8,'(1X,3(1X,E12.5))',advance="no") xr(irec),yr(irec),zs(irec)
           WRITE (8,*)''
        END DO
        
        ! write out points root tips
        DO  igrow=1,ngrow
           xr(nrec+igrow) = xg(igrow)+xplant(ipl)+1.*transtip(igrow,1,nsub(1,ipl),ipl)*dxgrid*nx
           yr(nrec+igrow) = yg(igrow)+yplant(ipl)+1.*transtip(igrow,2,nsub(1,ipl),ipl)*dygrid*ny
           WRITE (8,'(1X,3(1X,E12.5))',advance="no") xr(nrec+igrow),yr(nrec+igrow),zg(igrow)
           WRITE (8,*)''
        END DO 

        kk=0
        DO j=1,nrec !segments
           trans_i = transroot(irecpr(j),:,nsub(irecpr(j),ipl),ipl)
           trans_j = transroot(j,:,nsub(j,ipl),ipl)
           IF (ALL(trans_i.EQ.trans_j)) THEN
              if(irecpr(j).ne.0) then
                 kk = kk +1
              end if
              
           ENDIF
        end DO
        
        DO j=1,nbr! tips
           irecn=nrec
           DO WHILE (ibrseg(irecn).NE.j)
              irecn=irecn-1
           END DO
           DO igrow=1,ngrow
              IF (ibrgrw(igrow)==j) THEN
                 trans_i = transtip(igrow,:,nsub(nrec+igrow,ipl),ipl)
                 trans_j = transroot(irecn,:,nsub(irecn,ipl),ipl)
                 IF (ALL(trans_i.EQ.trans_j)) THEN
                    kk = kk +1
                 end if
              ENDIF
           END DO
        end DO
        
        
        WRITE (8,'(A6,1X,I8,1X,I9)')'LINES ',kk, kk*3
        
        DO j=1,nrec !tips
           trans_i = transroot(irecpr(j),:,nsub(irecpr(j),ipl),ipl)
           trans_j = transroot(j,:,nsub(j,ipl),ipl)
           IF (ALL(trans_i.EQ.trans_j)) THEN
              if(irecpr(j).ne.0) then
                 WRITE (8,'(3(1X,I8))',advance="no") 2, irecpr(j)-1, j-1
                 WRITE (8,*)''
              end if
              
           ENDIF
        ENDDO
        
        
        DO j=1,nbr !segments
           irecn=nrec
           DO WHILE (ibrseg(irecn).NE.j)
              irecn=irecn-1
           END DO
           DO igrow=1,ngrow
              IF (ibrgrw(igrow)==j) THEN
                 trans_i = transroot(nrec+igrow,:,nsub(nrec+igrow,ipl),ipl)
                 trans_j = transroot(irecn,:,nsub(irecn,ipl),ipl)
                 IF (ALL(trans_i.EQ.trans_j)) THEN
                    WRITE (8,'(3(1X,I8))',advance="no") 2, irecn-1, nrec+igrow-1
                    WRITE (8,*)''
                 end if
              ENDIF
           END DO
        end DO
        
        
        
     ELSE
        
        ! write points segments
        DO  irec=1,nrec
           WRITE (8,'(1X,3(1X,1pE10.3))',advance="no") xs(irec)+xplant(ipl),ys(irec)+yplant(ipl),zs(irec)
           WRITE (8,*)''
        END DO
        
        ! write points tips
        DO  igrow=1,ngrow
           WRITE (8,'(1X,3(1X,E12.5))',advance="no") xg(igrow)+xplant(ipl),yg(igrow)+yplant(ipl),zg(igrow)
           WRITE (8,*)''
        END DO

        sizen = 0
        DO j=1,nrec !segments
           if(irecpr(j).ne.0) then
              sizen = sizen +1
           ENDIF
        ENDDO
        
        DO j=1,nbr !tips
           irecn=nrec
           DO WHILE (ibrseg(irecn).NE.j)
              irecn=irecn-1
           END DO
           DO igrow=1,ngrow
              IF (ibrgrw(igrow)==j) THEN
                 sizen = sizen +1
              END IF
           end DO
        end do
        
        
        
        WRITE (8,'(A6,1X,I8,1X,I9)')'LINES ', sizen, sizen*3
        
        DO j=1,nrec !segments
           if(irecpr(j).ne.0) then
              WRITE (8,'(3(1X,I8))',advance="no") 2, irecpr(j)-1, j-1
              WRITE (8,*)''
           end if
        ENDDO
        
        DO j=1,nbr !tips
           irecn=nrec
           DO WHILE (ibrseg(irecn).NE.j)
              irecn=irecn-1
           END DO
           DO igrow=1,ngrow
              IF (ibrgrw(igrow)==j) THEN
                 WRITE (8,'(3(1X,I8))',advance="no") 2, irecn-1, nrec+igrow-1
                 WRITE (8,*)''
              ENDIF
           END DO
        end DO
        
     ENDIF !end if continu
         
     WRITE (8,'(A11,1X,I6)')'POINT_DATA ', nrec+ngrow
     IF(.NOT.lno_RWU) THEN
        WRITE (8,'(A17)')'SCALARS Lr float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        WRITE (8,'(1X,1pE11.4)',advance="no") Lr(1:nrec)
        WRITE (8,'(1X,1pE11.4)',advance="no") Lr(irecsg(1:ngrow))
        WRITE (8,*)''
        
        WRITE (8,'(A18)')'SCALARS Khr float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        WRITE (8,'(1X,1pE11.4)',advance="no") Khr(1:nrec)
        WRITE (8,'(1X,1pE11.4)',advance="no") Khr(irecsg(1:ngrow))
        WRITE (8,*)''
        
        WRITE (8,'(A21)')'SCALARS PHinter float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        WRITE (8,'(1X,1pE11.4)',advance="no") PHs(1:nrec,ipl)
        WRITE (8,'(1X,1pE11.4)',advance="no") PHs(irecsg(1:ngrow),ipl)
        WRITE (8,*)''
        
        WRITE (8,'(A21)')'SCALARS PHxylem float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        WRITE (8,'(1X,1pE11.4)',advance="no") PHr(2:nrec+1,ipl)
        WRITE (8,'(1X,1pE11.4)',advance="no") Phr(irecsg(1:ngrow)+1,ipl)
        WRITE (8,*)''
        
        WRITE (8,'(A28)')'SCALARS radialRootFlow float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        WRITE (8,'(1X,1pE11.4)',advance="no") sinkR(1:nrec,ipl)
        WRITE (8,'(1X,1pE11.4)',advance="no") sinkR(irecsg(1:ngrow),ipl)
        WRITE (8,*)''
        
        WRITE (8,'(A28)')'SCALARS axialRootFlow float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        WRITE (8,'(1X,1pE11.4)',advance="no") axialRootFlow(1:nrec,ipl)
        WRITE (8,'(1X,1pE11.4)',advance="no") axialRootFlow(irecsg(1:ngrow),ipl)
        WRITE (8,*)''
        
        WRITE (8,'(A17)')'SCALARS Qi float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        WRITE (8,'(1X,1pE11.4)',advance="no") Qi(1:nrec,ipl)
        WRITE (8,'(1X,1pE11.4)',advance="no") Qi(irecsg(1:ngrow),ipl)
        WRITE (8,*)''
        
        WRITE (8,'(A17)')'SCALARS Qd float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        WRITE (8,'(1X,1pE11.4)',advance="no") Qd(1:nrec,ipl)
        WRITE (8,'(1X,1pE11.4)',advance="no") Qd(irecsg(1:ngrow),ipl)
        WRITE (8,*)''
        
        WRITE (8,'(A17)')'SCALARS Qb float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        WRITE (8,'(1X,1pE11.4)',advance="no") Q_bc1(1:nrec,ipl)
        WRITE (8,'(1X,1pE11.4)',advance="no") Q_bc1(irecsg(1:ngrow),ipl)
        WRITE (8,*)''
        
        WRITE (8,'(A19)')'SCALARS rconc float '
        WRITE (8,'(A20)')'LOOKUP_TABLE default'
        WRITE (8,'(1X,1pE11.4)',advance="no") segconc(1:nrec)
        WRITE (8,'(1X,1pE11.4)',advance="no") segconc(irecsg(1:ngrow))
           WRITE (8,*)''

        IF(l_linSorb .or. l_freundSorb) THEN
           WRITE (8,'(A21)')'SCALARS segsorb float '
           WRITE (8,'(A20)')'LOOKUP_TABLE default'
           WRITE (8,'(1X,1pE11.4)',advance="no") segsorb(1:nrec)
           DO igrow=1,ngrow
              WRITE (8,'(1X,1pE11.4)',advance="no") 0E+00
           END DO
           WRITE (8,*)''
        END IF
     END IF
          
     WRITE (8,'(A19)')'SCALARS age float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE11.4)',advance="no") t-timorg(1:nrec)
     DO irec=1,ngrow
        WRITE (8,'(1X,1pE11.4)',advance="no") 0E+00
     END DO
     WRITE (8,*)''
     
     WRITE (8,'(A28)')'SCALARS segmentRadius float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     DO  irec=1,nrec
        WRITE (8,'(1X,1pE14.8)',advance="no") segrad(irec) 
     END DO
     DO  irec=1,ngrow
        WRITE (8,'(1X,1pE14.8)',advance="no") 0E+00 
     end DO
     WRITE (8,*)''
     
     CLOSE(8)
     
     RETURN
  END DO
END SUBROUTINE OutDouVTK
!************************************************************************
!> generates output in case partrace was running
SUBROUTINE PartraceOut(t)
  USE typedef
  USE Soldata
  USE GridData
  IMPLICIT NONE
  REAL(sp),INTENT(in)::t
  INTEGER(ap):: i,j,k,l,ix,iy,iz
  CHARACTER file*31,file2*19,file3*36
  LOGICAL,save::lfirst=.true.
  !> \param t current simulation time
  
  IF (lfirst) THEN
     OPEN (UNIT=8,FILE='out/Partrace/out.velocity.00001',STATUS='UNKNOWN')
     OPEN (UNIT=10,FILE='out/Partrace/out.water_content.00001',STATUS='UNKNOWN')
     
     WRITE (10,'(a14)') '#water content'
     WRITE (10,'(a15)') '#CPU: 1 from: 1'
     
     WRITE (8,'(a9)') '#velocity'
     WRITE (8,'(a15)') '#CPU: 1 from: 1'
     IF (continu) THEN
        WRITE (10,'(A26,I6,A9,I6)') '#Global no. of nodepoints: ', (nx+1)*(ny+1)*nz, '  local: ',(nx+1)*(ny+1)*nz
        WRITE (10,'(A22,I6,A8,I6)') '#Nodepoints x: 1  to: ',nx+1,'  from: ',nx+1
        WRITE (10,'(A22,I6,A8,I6)') '#Nodepoints y: 1  to: ',ny+1,'  from: ',ny+1
        WRITE (10,'(A22,I6,A8,I6)') '#Nodepoints z: 1  to: ',nz,'  from: ',nz
        
        WRITE (8,'(A26,I6,A9,I6)') '#Global no. of nodepoints: ', (nx+1)*(ny+1)*nz, '  local: ',(nx+1)*(ny+1)*nz
        WRITE (8,'(A22,I6,A8,I6)') '#Nodepoints x: 1  to: ',nx+1,'  from: ',nx+1
        WRITE (8,'(A22,I6,A8,I6)') '#Nodepoints y: 1  to: ',ny+1,'  from: ',ny+1
        WRITE (8,'(A22,I6,A8,I6)') '#Nodepoints z: 1  to: ',nz,'  from: ',nz
     ELSE
        WRITE (10,'(A26,I6,A9,I6)') '#Global no. of nodepoints: ', nx*ny*nz, '  local: ',nx*ny*nz
        WRITE (10,'(A22,I6,A8,I6)') '#Nodepoints x: 1  to: ',nx,'  from: ',nx
        WRITE (10,'(A22,I6,A8,I6)') '#Nodepoints y: 1  to: ',ny,'  from: ',ny
        WRITE (10,'(A22,I6,A8,I6)') '#Nodepoints z: 1  to: ',nz,'  from: ',nz
        
        WRITE (8,'(A26,I6,A9,I6)') '#Global no. of nodepoints: ', nx*ny*nz, '  local: ',nx*ny*nz
        WRITE (8,'(A22,I6,A8,I6)') '#Nodepoints x: 1  to: ',nx,'  from: ',nx
        WRITE (8,'(A22,I6,A8,I6)') '#Nodepoints y: 1  to: ',ny,'  from: ',ny
        WRITE (8,'(A22,I6,A8,I6)') '#Nodepoints z: 1  to: ',nz,'  from: ',nz
     ENDIF
     CLOSE(8)
     CLOSE(10)
     lfirst = .false.
  ENDIF
  
  WRITE (file,'(A31)')'out/Partrace/out.velocity.00001'
  WRITE (file2,'(A19)')'out/Partrace/geoPar'
  WRITE (file3,'(A36)')'out/Partrace/out.water_content.00001'
  
  OPEN (UNIT=8,FILE=file,STATUS='OLD',POSITION='APPEND')
  OPEN (UNIT=9,FILE=file2,STATUS='UNKNOWN')
  OPEN (UNIT=10,FILE=file3,STATUS='OLD',POSITION='APPEND')      
  
  WRITE (10,'(A19,F0.3)') '#simulation time:  ',t-100
  WRITE (8,'(A19,F0.3)') '#simulation time:  ',t-100
  
  CALL Veloc
  !> write nodal soil values:
  IF (continu) THEN
     k=0
     j=0
     DO  iz=nz,1,-1
        DO  iy=1,ny+1
           DO  ix=1,nx+1     
              k = ix + (iy-1)*nx + (iz-1)*nx*ny
              j = j +1
              IF (ix.EQ.nx+1) THEN
                 k = k - nx
              END IF
              IF  (iy.EQ.ny+1) THEN
                 k = k - (nx*ny)
              END IF
              WRITE (8,'(I7,8X,4(1X,1pE11.4))')j,Vx(k)/theta(k),Vy(k)/theta(k),Vz(k)/theta(k),theta(k)
              WRITE (9,'(I7,8X,3(1X,1pE11.4))')j,xGrid(k),yGrid(k),zGrid(k)
              WRITE (10,'(I7,8X,1(1X,1pE11.4))')j,theta(k)
           END DO
        END DO
     END DO
  ELSE
     l=0
     DO j=1,nz
        DO i=1,nx*ny
           k = nx*ny*nz - j*nx*ny +i
           l = l+1
           WRITE (8,'(I7,8X,4(1X,1pE11.4))')l,Vx(k)/theta(k),Vy(k)/theta(k),Vz(k)/theta(k),theta(k)
           WRITE (9,'(I7,8X,3(1X,1pE11.4))')l,xGrid(k),yGrid(k),zGrid(k)
           WRITE (10,'(I7,8X,1(1X,1pE11.4))')l,theta(k)
        END DO
     END DO
  ENDIF
  
  
  CLOSE(8)
  CLOSE(9)
  CLOSE(10)
  RETURN
END SUBROUTINE PartraceOut
!****************************************************************************
!> ### end of simulation ###
SUBROUTINE Getout
  WRITE(*,'(//'' End of simulation. Program terminated normally.''///)')
  STOP
END SUBROUTINE Getout
!****************************************************************************
!> UNUSED AT THE MOMENT 
SUBROUTINE OutFEM_VTK_Sub(kOut)
  USE GridData
  USE SolData
  USE WatFun
  IMPLICIT NONE
  
  INTEGER(sp), intent(in)::kout
  INTEGER(sp)::ie,i,ise,global
  CHARACTER file*20
  !> \param kout current number for the FEM output
  
  WRITE (file,'(A20)')'out/outfem_sub00.vtk'
  IF (kout.LT.10) THEN
     WRITE (file(16:16),'(I1)') kout
  ELSE
     WRITE (file(15:16),'(I2)') kout
  ENDIF
  
  PRINT *, file
  
  OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
  WRITE (8,'(A26)')'# vtk DataFile Version 3.0' 
  WRITE (8,'(A12)')'model R-SWMS'
  WRITE (8,'(A5)')'ASCII'
  WRITE (8,'(A25)')'DATASET UNSTRUCTURED_GRID'
  
  WRITE (8,'(/A6,1X,I7, 1X, A5)')'POINTS', nPt, 'float'
  do i=1,nPt
     WRITE (8,'(3(1X,1pE12.5))',advance="no") xgrid(i), ygrid(i) ,zGrid(i)
     WRITE (8,*)''
  end do
  
  
  WRITE (8,'(/A5,1X,I7,1X,I8)')'CELLS', nElm*5, 5*nElm*5
  do iE=1,nElm
     do iSE=1,5
        WRITE (8,'((1X,I6))',advance="no") 4
        do global=1,4
           WRITE (8,'((1X,I6))',advance="no") elmnod(iL(global,iSE,subN(iE)), iE)-1
        end do
        WRITE (8,*)''
     end do
  enddo
  
  WRITE (8,'(/A10,1X,I7)')'CELL_TYPES', nElm*5
  do iE=1,nElm
     do iSe=1,5
        WRITE (8,'(8(1X,I2))',advance="no") 10
     end do
  enddo
  
  WRITE (8,*)''
  WRITE (8,'(A11,1X,I7)')'POINT_DATA ', nPt
  WRITE (8,'(A30)')'SCALARS pressurehead float '
  WRITE (8,'(A20)')'LOOKUP_TABLE default'
  WRITE (8,'(1X,1pE15.6)',advance="no") hnew(1:nPt)
  
  WRITE (8,*)''
  WRITE (8,'(A30)')'SCALARS wc float '
  WRITE (8,'(A20)')'LOOKUP_TABLE default'
  WRITE (8,'(1X,1pE15.6)',advance="no") theta(1:nPt)
  
  WRITE (8,*)''
  WRITE (8,'(A30)')'SCALARS conc float '
  WRITE (8,'(A20)')'LOOKUP_TABLE default'
  WRITE (8,'(1X,1pE15.6)',advance="no") conc(1:nPt)
  
  WRITE (8,*)''
  WRITE (8,'(A30)')'SCALARS sink float '
  WRITE (8,'(A20)')'LOOKUP_TABLE default'
  WRITE (8,'(1X,1pE15.6)',advance="no") sink(1:nPt)
  
  WRITE (8,*)''
  WRITE (8,'(A30)')'SCALARS csink float '
  WRITE (8,'(A20)')'LOOKUP_TABLE default'
  WRITE (8,'(1X,1pE15.6)',advance="no") csink(1:nPt)
  
  WRITE (8,*)''
  WRITE (8,'(A11,1X,I7)')'CELL_DATA ', nElm*5
  WRITE (8,'(A30)')'SCALARS subN float '
  WRITE (8,'(A20)')'LOOKUP_TABLE default'
  do iE=1,nElm
     do iSe=1,5
        WRITE (8,'(1X,I7)',advance="no") subN(iE)
     end do
  enddo
  
  
  CLOSE (8)
  RETURN
END SUBROUTINE OutFEM_VTK_Sub
!****************************************************************************
!> generates .vtk outputs for particles within roots, in case of hormone transport; ParRoot.XX 
SUBROUTINE OutParticleVTK(kout,t)
  USE typedef
  USE RootData, ONLY: seglen,segsur,segconc,segSoluteMass,nrec,ngrow,nbr,xg,yg,zg,xs,ys,zs,xplant,yplant,irecpr,ibrseg,ibrgrw
  USE DoussanMat, ONLY:  veloRoot,Lr,Khr,PHs,PHr,sinkR,axialRootFlow,Qi,nplant,Qd,Q_bc1,nsub,transroot
  USE GridData, ONLY: dxgrid,dygrid,dzgrid,nx,ny,nz,continu
  USE SoluteRootMat, ONLY: Particle,firstP,totalParticleNum
  USE PlntData, ONLY:Tact
  
  IMPLICIT NONE
  
  REAL(sp) :: xr(nrec+ngrow), yr(nrec+ngrow), zr(nrec+ngrow)
  REAL(sp) :: seg_tipx(nrec),seg_tipy(nrec),seg_tipz(nrec)
  REAL(sp) :: t
  INTEGER(ap) :: irec,kout,ipl, rstype
  INTEGER(ap):: igrow,irecn,iseg,j
  CHARACTER file*29 
  TYPE(Particle) ,POINTER  :: P
  !> \param t current simulation time
  !> \param rOuR current number for the root output
 
  DO ipl=1,nplant
     
     WRITE (file,'(A18)')'out/vtk/ParRootXX.'
     
     
     IF (ipl.LT.10) THEN
        WRITE (file(16:16),'(I1)') 0
        WRITE (file(17:17),'(I1)') ipl
     ELSE
        WRITE (file(16:17),'(I2)') ipl
     END IF
     IF (kout.LT.10) THEN
        WRITE (file(19:19),'(I1)') kout
        WRITE (file(20:24),'(A4)') '.vtk'
     ELSEIF (kout.LT.100) THEN
        WRITE (file(19:20),'(I2)') kout
        WRITE (file(21:25),'(A4)') '.vtk'
     ELSE
        WRITE (file(19:21),'(I3)') kout
        WRITE (file(22:26),'(A4)') '.vtk'
     ENDIF
     
     PRINT *, file
     
     OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
     
     IF (continu) THEN
        
        ! write out points root segments
        DO  irec=1,nrec
           xr(irec) = xs(irec)+xplant(ipl)+1.*transroot(irec,1,nsub(irec,ipl),ipl)*dxgrid*nx
           yr(irec) = ys(irec)+yplant(ipl)+1.*transroot(irec,2,nsub(irec,ipl),ipl)*dygrid*ny
           zr(irec) = zs(irec)
           !WRITE (8,'(1X,3(1X,E12.5))',advance="no") xr(irec),yr(irec),zs(irec)
           !WRITE (8,*)''
        END DO
        
        ! write out points foot tips
        DO  igrow=1,ngrow
           xr(nrec+igrow) = xg(igrow)+xplant(ipl)+1.*transroot(nrec+igrow,1,nsub(1,ipl),ipl)*dxgrid*nx
           yr(nrec+igrow) = yg(igrow)+yplant(ipl)+1.*transroot(nrec+igrow,2,nsub(1,ipl),ipl)*dygrid*ny
           zr(nrec+igrow) = zg(igrow)
           !WRITE (8,'(1X,3(1X,E12.5))',advance="no") xr(nrec+igrow),yr(nrec+igrow),zg(igrow)
           !WRITE (8,*)''
        END DO
        
        
        
     ELSE !if not continue
        
        ! points segments
        DO  irec=1,nrec
           xr(irec) = xs(irec)+xplant(ipl)
           yr(irec) = ys(irec)+yplant(ipl)
           zr(irec) = zs(irec)
        END DO
        
        ! points tips
        DO  igrow=1,ngrow
           xr(nrec+igrow) = xg(igrow)+xplant(ipl)
           yr(nrec+igrow) = yg(igrow)+yplant(ipl)              
           zr(nrec+igrow) = zg(igrow)        
        END DO
        
     END IF !continu 
     
     DO j=1,nrec !segments
        IF(irecpr(j).NE.0) THEN
           seg_tipx(irecpr(j)) = (xr(j)-xr(irecpr(j)))/seglen(irecpr(j))
           seg_tipy(irecpr(j)) = (yr(j)-yr(irecpr(j)))/seglen(irecpr(j))
           seg_tipz(irecpr(j)) = (zr(j)-zr(irecpr(j)))/seglen(irecpr(j))
        END IF
     ENDDO
     
     DO j=1,nbr !tips
        irecn=nrec
        DO WHILE (ibrseg(irecn).NE.j)
           irecn=irecn-1
        END DO
        DO igrow=1,ngrow
           IF (ibrgrw(igrow)==j) THEN
              
              seg_tipx(irecn) = (xr(nrec+igrow)-xr(irecn))/seglen(irecn)
              seg_tipy(irecn) = (yr(nrec+igrow)-yr(irecn))/seglen(irecn)
              seg_tipz(irecn) = (zr(nrec+igrow)-zr(irecn))/seglen(irecn)
              
           ENDIF
        END DO
     END DO
     
     ALLOCATE(p)
     NULLIFY(p)
     
     WRITE (8,'(A26)')'# vtk DataFile Version 3.0' 
     WRITE (8,'(A12)')'model R-SWMS'
     WRITE (8,'(A5)')'ASCII'
     WRITE (8,'(A25)')'DATASET UNSTRUCTURED_GRID'
     WRITE (8,'(A7,1X,I8,1X,A5)')'POINTS ',totalParticleNum  , 'float'   
     p => firstP
     DO WHILE(ASSOCIATED(P))
        WRITE (8,'(3(1X,1pE11.4))',advance="no") xr(p%segnum)+seg_tipx(p%segnum)*(seglen(p%segnum)-p%position), &
             yr(p%segnum)+seg_tipy(p%segnum)*(seglen(p%segnum)-p%position), &
             zr(p%segnum)+seg_tipz(p%segnum)*(seglen(p%segnum)-p%position)
        WRITE (8,*)''
        P => P%NEXT
     END DO
     
     WRITE (8,'(A8,1X,I8,1X,I8)')'CELLS ',totalParticleNum ,totalParticleNum*2
     DO j=1,totalParticleNum 
        WRITE (8,'(2(1X,I7))',advance="no") 1,j-1
        WRITE (8,*)''
     END DO
     
     WRITE (8,'(A16,1X,I8)')'CELL_TYPES ',totalParticleNum 
     DO j=1,totalParticleNum 
        WRITE (8,'(1X,I7)',advance="no") 1
     END DO
     WRITE (8,*)''
     
     WRITE (8,'(A11,1X,I8)')'POINT_DATA ', totalParticleNum
     WRITE (8,'(A17)')'SCALARS ID float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     p=> firstP
     DO WHILE(ASSOCIATED(P))
        WRITE (8,'((1X,I8))',advance="no") p%ID
        P => P%NEXT
     END DO
     WRITE (8,*)''
     
     WRITE (8,'(A19)')'SCALARS Mass float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     p=> firstP
     DO WHILE(ASSOCIATED(P))
        WRITE (8,'((1X,1pE11.4))',advance="no") p%mass
        P => P%NEXT
     END DO
     WRITE (8,*)''
     
     WRITE (8,'(A21)')'SCALARS SegNum float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     p=> firstP
     DO WHILE(ASSOCIATED(P))
        WRITE (8,'((1X,I8))',advance="no") p%segnum
        P => P%NEXT
     END DO
     WRITE (8,*)''
     
     WRITE (8,'(A21)')'SCALARS PartAge float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     p=> firstP
     DO WHILE(ASSOCIATED(P))
        WRITE (8,'((1X,1pE11.4))',advance="no") t-p%partOrig
        P => P%NEXT
     END DO
     WRITE (8,*)''
     
     CLOSE(8)
     
     RETURN
  END DO
END SUBROUTINE OutParticleVTK
!****************************************************************************
