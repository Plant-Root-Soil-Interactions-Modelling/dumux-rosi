! ==============================================================================
! Source file INPUT ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
!> ### reads main inputs ###
SUBROUTINE Applic(dt)
!> \param nPt = total number of nodal points
!> \param nBCpts = total number of nodal points with specified BC
!> \param nElm = total number of elements
!> \param ne* = number of elements (half-cuboids) in '*'-direction
!> \param nel = nex*ney (number of elements per horizontal element-layer)
  USE iso_c_binding
  USE Typedef
  USE ParamData, ONLY: ldirect,lvtk,lOutPartrace,lChem,lretry,last_out,&
       lSalinity,maxmat,maxIrrig,maxbdr,mxBcCh,lPartUp
  USE RootData, ONLY: lno_root_growth,lno_RWU,ltemp,lCalloc,lrrs,lrrt,&
       lRootTyp_growth,lSomma_growth,lUpdate_growth,lCou,lDou,lno_Archi,&
       lFed,lSinkCube
  USE GridData
  USE tmctrl
  USE BoundData
  USE SolData, ONLY: hold, htemp, hnew, conc, Kode, matNum, nmat, KodCB, lTab
  USE CumData
  USE DomData
  USE ObsData 
  USE doussanmat, ONLY : old,oldT,counter2,PH_micro2,Phi_mat,h_mat2,ana_aan
  USE RhizoData, ONLY: lRhizo
  USE SoluteRootMat, ONLY: timestepfactor,uptakeorder
  IMPLICIT NONE
  
  REAL(sp) :: Qbcrec(mxBcCh)
  CHARACTER text*5
  CHARACTER infile*14
  INTEGER(ap) :: nh,nQ,nI,nFrDr,l,j,ise,ie,ip,i,idum,k,err,ii,nl,aa,dummy
  INTEGER(ap) :: integ_zprof,i2
  INTEGER(ap) :: assCode,hfun,nQ1
  INTEGER(ap) :: xIrrig(maxIrrig),yIrrig(maxIrrig),zIrrig(maxIrrig)
  REAL(sp) :: A11,A22,A33,A12,A13,A23,C11,C22,C33,ini,ixmin,ixymin,imin,dt,xqmin,xqmax,xhmin,xhmax
  LOGICAL :: lOrt=.FALSE.,firstOK=.FALSE.,ltop=.TRUE.,lFrDr=.FALSE.
  Integer(ap), ALLOCATABLE, DIMENSION(:) ::node_temp
  ALLOCATE(node_temp(nPt))
  nMat=0
  xmax=-1.E+30_dp
  xmin=+1.E+30_dp
  ymax=-1.E+30_dp
  ymin=+1.E+30_dp
  zmax=-1.E+30_dp
  zmin=+1.E+30_dp
  profOK=.FALSE.
  OPEN(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
  !> #### reads control.in ####
  OPEN (Unit=10,FILE='in/control.in',STATUS='OLD',ERR=10)
  WRITE(15,'(//''++ General Information ++'')',advance='no')
  WRITE(15,'(/''-------------------------'')',advance='no')
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) LnUnit,TmUnit,MsUnit,CnUnit
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) itMax,itMaxRoot
  READ (10,*)
  READ (10,*)
  READ (10,*) RelEps,factorRelEps
  READ (10,*)
  READ (10,*)
  READ (10,*) epslonPH,epslonWC,epslonR,epslonS
  IF (RelEps)   WRITE(15,'(/''* Relative tolerance is used for WC, PH and SInk with value of'',f5.3)',advance='no')1/factorRelEps
  IF (factorRelEps.EQ.0) THEN
     STOP 'Control.in: factor for relative criterium may not be zero'
  ENDIF
  READ (10,*)
  READ (10,*)
  READ (10,*) dt,dtMin,dtMax,FacInc,FacDec,dtRoot 
  READ (10,*)
  READ (10,*)
  READ (10,*) lretry,last_out
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) nOut
  IF (nOut.GT.mxOut) STOP 'Input from  < control.in >  exceeds parameter "mxOut". Program terminated.'
  READ (10,*)
  READ (10,*)
  READ (10,*) (tOut(i),i=1,nOut)
  tmax=tOut(nOut)
  READ (10,*)
  READ (10,*)
  READ (10,*) lvtk,lOutpartrace,profOK
  READ (10,*)
  READ (10,*)
  READ (10,*) dtprof
  IF(profOK) THEN
     IF (dtprof<999) THEN
        nouProf=0
        ini=0
        DO  WHILE (ini.LT.tmax)
           nouProf=nouProf+1
           tOuProf(nouProf)=ini+dtprof
           ini=tOuProf(nouProf)
           IF (nouProf.GT.mxProf) THEN
              PRINT *,'too small time step in z-profiles: only ',mxProf,' profiles will be kept'
              GOTO 77
           ENDIF
        END DO
     END IF
  END IF
77 READ(10,*)
  READ(10,*)
  READ(10,*)
  READ(10,*)
  READ(10,*) lno_RWU,lFed,lDou,lCou,lSinkCube
  READ(10,*)
  READ(10,*)
  READ(10,*)
  READ(10,*) lno_Archi,lrrs,lrrt
  READ(10,*)
  READ(10,*)
  READ(10,*)
  READ(10,*) lno_root_growth,lRootTyp_growth,lSomma_growth,lUpdate_growth
  READ(10,*)
  READ(10,*)
  READ(10,*)
  READ(10,*) lCalloc,lChem,ltemp,continu,lSalinity,lPartUp, lRhizo
  READ(10,*)
  READ(10,*)
  READ(10,*)
  READ(10,*) ldirect,old,ana_aan,ltab
  CLOSE (10)

  oldT=old
  
  !> check for compatibility of processes
  if((.not.lno_RWU) .and. lno_archi .and. .NOT.lCou)  STOP 'RWU models should always be associated with a root architecture'
  if(.not.lCou .and. .not.lSinkCube)  STOP 'Sink nodes currently only programmed to generate SSF for Hydrus -> uses Couvreur & archi from RootSys'
  if(.not.lno_RWU .and. .not.lFed .and. .not.lDou .and. .not.lCou) STOP 'Please set at least one RWU option in control.in to .TRUE.'
  if(lRootTyp_growth) STOP 'RootTyp root growth during the scenario run is not yet implemented'
  if(lCalloc .and. .not. (lSomma_growth.and.lDou)) STOP 'Assimilate allocation can currently only be associated with Somma root growth and Doussan root water uptake'
  if(continu) WRITE(15,'(/''* XY domain periodicity for soil water flow, root system architecture and root water flow'')',advance='no')
  if (continu.and.(.not.old)) STOP 'Domain continuity currently only works with soil Averaging method'
  if(ldirect .and. .not. lDou) STOP 'Direct Doussan can only be used with Doussan RWU model'
  if(ldirect .and. .not. lno_root_growth) STOP 'Direct Doussan can only be used without root growth'
  IF(ana_aan .AND. .NOT. old) STOP 'Cannot use microscopic model and memory reduction simultaneous'

  !> #### reads grid.in ####
  OPEN (Unit=10,FILE='in/grid.in',STATUS='OLD',ERR=20)
  WRITE(15,'(//''++ Soil Geometry ++'')',advance='no')
  WRITE(15,'(/''--------------------'')',advance='no')
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) nPt,nElm
  ALLOCATE(hOld(nPt))
  ALLOCATE(hTemp(nPt))
  ALLOCATE(hNew(nPt))
  ALLOCATE(conc(nPt))
  ALLOCATE(Kode(nPt))
  IF (nPt.GT.maxnod) STOP 'Input from  < grid.in >  exceeds parameter "maxnod". Program terminated.'
  IF (nElm.GT.maxelm) STOP 'Input from  < grid.in >  exceeds parameter "maxelm". Program terminated.'
  lOrt=.TRUE.
  READ (10,*)
  READ (10,*)
  READ (10,*) nx,ny,nz,nex,ney,nez,dxGrid,dyGrid,dzGrid
  nel=nex*ney
  nl=nx*ny
  CALL IniGrid !> initialize xgrid, ygrid,zgrid and scaling factors
  ! elements:
  READ (10,*)
  READ (10,*)
  READ (10,*)
  DO i=1,nElm
     READ(10,*) iDum,(elmnod(k,i),k=1,8),subN(i),A11,A22,A33,A12,A13,A23,C11,C22,C33
     ConAxx(i)=C11*A11*A11+C22*A12*A12+C33*A13*A13
     ConAyy(i)=C11*A12*A12+C22*A22*A22+C33*A23*A23
     ConAzz(i)=C11*A13*A13+C22*A23*A23+C33*A33*A33
     ConAxy(i)=C11*A11*A12+C22*A12*A22+C33*A13*A23
     ConAxz(i)=C11*A11*A13+C22*A12*A23+C33*A13*A33
         ConAyz(i)=C11*A12*A13+C22*A22*A23+C33*A23*A33
  END DO
  WRITE(15,'(/''* Number of elements in x, y and z directions:'',3i5)',advance='no')nex,ney,nez
  WRITE(15,'(/''* Number of nodes in x, y and z directions:'',3i5)',advance='no')nx,ny,nz
  WRITE(15,'(/''* Length of elements in x, y and z directions:'',3f7.2)',advance='no')dxGrid,dyGrid,dzGrid
  CLOSE (10)
  !> #### reads nodes.in ####
  IF (.NOT.lretry) THEN
     OPEN (Unit=10,FILE='in/nodes.in',STATUS='OLD',ERR=30)
  3  READ (10,'(A5)') text
     IF (text.NE.'Node#') GOTO 3
     DO i=1,nPt
        READ (Unit=10,iostat=err) iDum,MatNum(i),xGrid(i),yGrid(i), zGrid(i),hOld(i),Conc(i),Axy(i),Bxy(i),Dxy(i),Exy(i)
        IF (err.NE.0) THEN ! in case there are no scaling factors
           READ (10,*) iDum,MatNum(i),xGrid(i),yGrid(i),zGrid(i), hOld(i),Conc(i)
           Axy(i)=1.
           Bxy(i)=1.
           Dxy(i)=1.
           Exy(i)=1.
        ENDIF
        Kode(i)=0
        hNew(i)=hOld(i)
        hTemp(i)=hOld(i)
        nMat=MAX(nMat,MatNum(i))
        xmin=MIN(xmin,xGrid(i))
        ymin=MIN(ymin,yGrid(i))
        zmin=MIN(zmin,zGrid(i))
        xmax=MAX(xmax,xGrid(i))
        ymax=MAX(ymax,yGrid(i))
        zmax=MAX(zmax,zGrid(i))
     END DO
  ELSE
     IF (last_out.GT.0.AND.last_out.LT.10) THEN
        write(infile,'(A11,I1)') 'out/outfem.',last_out
     ELSEIF (last_out.LT.100) THEN
        write(infile,'(A11,I2)') 'out/outfem.',last_out
     ELSEIF (last_out.LT.1000) THEN
        write(infile,'(A11,I3)') 'out/outfem.',last_out
     ELSE
        STOP 'Last output # should be between 1 and 999, when retrying a simulation from an existing output'
     ENDIF
     OPEN (Unit=10,FILE=infile,STATUS='OLD',ERR=30)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     DO i=1,nPt ! in case there are no scaling factors
        READ (10,*) iDum,MatNum(i),xGrid(i),yGrid(i),zGrid(i), hOld(i),Conc(i)
        Axy(i)=1.
        Bxy(i)=1.
        Dxy(i)=1.
        Exy(i)=1.
        Kode(i)=0
        hNew(i)=hOld(i)
        hTemp(i)=hOld(i)
        nMat=MAX(nMat,MatNum(i))
        xmin=MIN(xmin,xGrid(i))
        ymin=MIN(ymin,yGrid(i))
        zmin=MIN(zmin,zGrid(i))
        xmax=MAX(xmax,xGrid(i))
        ymax=MAX(ymax,yGrid(i))
        zmax=MAX(zmax,zGrid(i))
     END DO
  ENDIF
  DO i=1,nx
     xCol(i)=xmin+(i-1)*dxGrid
  ENDDO
  DO i=1,ny
     yCol(i)=ymin+(i-1)*dyGrid
  ENDDO
  CLOSE (10)
  IF (nmat.GT.maxmat) STOP 'Input from nodes file exceeds parameter "maxmat". Program terminated.'
  DO iP=1,nPt
     ListNE(iP)=0
  END DO
  nBand=1	
  !> nBand calculated like in SWMS3D 
  DO  iE=1,nElm
     DO  iSE=1,5
        i=elmnod(iL(1,iSE,subN(iE)),iE)
        j=elmnod(iL(2,iSE,subN(iE)),iE)
        k=elmnod(iL(3,iSE,subN(iE)),iE)
        l=elmnod(iL(4,iSE,subN(iE)),iE)
        ListNE(i)=ListNE(i)+1
        ListNE(j)=ListNE(j)+1
        ListNE(k)=ListNE(k)+1
        ListNE(l)=ListNE(l)+1
        IF(ABS(i-j).GT.nBand) nBand=ABS(i-j)
        IF(ABS(i-k).GT.nBand) nBand=ABS(i-k)
        IF(ABS(i-l).GT.nBand) nBand=ABS(i-l)
        IF(ABS(j-k).GT.nBand) nBand=ABS(j-k)
        IF(ABS(j-l).GT.nBand) nBand=ABS(j-l)
        IF(ABS(k-l).GT.nBand) nBand=ABS(k-l)
     END DO
  END DO
  nBand=nBand+1
  PRINT *, 'band: ', nband

  !> #### reads mesh.in for additional geometrical information ####
  OPEN(Unit=10,FILE='in/mesh.in',STATUS='OLD',ERR=60)
  DO i=1,9
     READ (10,*)
  END DO
  READ (10,*) geom
  IF (geom.EQ.3) THEN
     READ (10,*)
     READ (10,*) rad_cyl
     ! grid center
     x_cent = (nex*dxGrid/2+xmin)
     y_cent = (ney*dyGrid/2+ymin)
  END IF
  CLOSE(10)

  !> #### reads bc.in; soil boundary conditions ####
  OPEN (Unit=10,FILE='in/bc.in',STATUS='OLD',ERR=40)
  WRITE(15,'(//''++ Soil Boundary Conditions ++'')',advance='no')
  WRITE(15,'(/''------------------------------'')',advance='no')
  ! flux BC: (t)-BC [L^3/T]  (Kode(i)=-1)
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) qfun,ltop
  READ (10,*)
  SELECT CASE(qfun)
     CASE(0) !> + case0: no flux
        READ (10,*) 
        nQ=0
     CASE(1) !> + case1: flux BC on total top surface
        READ(10,*)
        nQ=nx*ny
        IF (ltop) THEN
         DO i=1,nQ
           iBCPt(i)=i
         ENDDO
        ELSE
           ii=nPt-nQ+1
           DO i=1,nQ
              iBCPt(i)=ii
              ii=ii+1
           ENDDO
        END IF
        DO  i=1,nQ
           Kode(iBCPt(i))=-1
        END DO
     CASE(2) !> + case2: top domain is partially irrigated on a strip between xqmin and xqmax
        READ(10,*) xqmin, xqmax
        IF ((xqmin.LT.xmin).OR.(xqmax.GT.xmax)) STOP 'xmin or xmax for Q bc in <bc.in> not within soil domain. Program terminated.'
        xqmin=FLOOR(xqmin/dxgrid)*dxgrid
        xqmax=FLOOR((xqmax+dxgrid)/dxgrid)*dxgrid !incl. right nodes...
        nQ=(xqmax-xqmin)/dxgrid*ny
        aa=ABS(xqmax-xqmin)/dxgrid
        i=1
        DO j=1,ny
           DO k=1,aa
              iBCPt(i)=(xqmin-xmin)/dxgrid+k+(j-1)*nx
              i=i+1
           END DO
        END DO
        DO i=1,nQ
           Kode(iBCPt(i))=-1
        END DO
     CASE(3) !> + case3: top domain is partially irrigated on two strips between xqmin1 and xqmax1 and from xqmin2 and xqmax2
        homogene = 2 !non-homogeneous distribution of BC -- time and space dependent changes of flux
        READ(10,*) xqmin1, xqmax1, xqmin2, xqmax2
        IF ((xqmin1.LT.xmin).OR.(xqmin2.LT.xmin).OR.(xqmax1.GT.xmax).OR.(xqmax2.GT.xmax)) STOP 'boundaries for partial root zone irrigation for Q bc in <bc.in> not within soil domain. Program terminated.'
        IF((xqmin2.LT.xqmin1).OR.(xqmax1.GT.xqmax2)) STOP 'xmin1 & xmax1 should be smaller than xmin2 & xmax2, Program terminated from <bc.in> '
        IF(xqmin2.LE.xqmax1) STOP 'no double allocation of the same nodes possible for partial root zone irrigation; Program terminated from <bc.in>'
        xqmin1=FLOOR(xqmin1/dxgrid)*dxgrid
        xqmax1=FLOOR(xqmax1/dxgrid)*dxgrid 
        xqmin2=FLOOR((xqmin2+dxgrid)/dxgrid)*dxgrid 
        xqmax2=FLOOR((xqmax2+dxgrid)/dxgrid)*dxgrid 
        nQ1=(xqmax1-xqmin1)/dxgrid*ny
        aa=ABS(xqmax1-xqmin1)/dxgrid
        i=1
        DO j=1,ny
           DO k=1,aa
              iBCPt(i)=(xqmin1-xmin)/dxgrid+k+(j-1)*nx
              i=i+1
           END DO
        END DO         
        !=nQ1+1start allocating nodes from the end of the first part
        nQ=nQ1+(xqmax2-xqmin2)/dxgrid*ny
        aa=ABS(xqmax2-xqmin2)/dxgrid
        DO j=1,ny
           DO k=1,aa
              iBCPt(i)=(xqmin2-xmin)/dxgrid+k+(j-1)*nx
              i=i+1
           END DO
        END DO
        DO i=1,nQ
           Kode(iBCPt(i))=-1
        END DO
        CASE(4) !> + case4: irrigation only on material #1 (e.g. cylindrical domain where outer filling material = #2)
           READ(10,*)
           dummy = nx*ny
           nQ = 0
           DO i=1,dummy
              IF(MatNum(i).EQ.1) THEN
                 nQ = nQ+1
                 iBCPt(nQ)=i
              END IF
           ENDDO
           DO  i=1,nQ
              Kode(iBCPt(i))=-1
           END DO
     CASE DEFAULT 
        STOP 'Sorry, qfun in <bc.in> can only take a value between 0 and 4. Program terminated.'
     END SELECT
  READ (10,*)
  READ (10,*) nQbcCh
  IF (nQbcCh.GT.mxBcCh) STOP 'Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.'

  READ (10,*)
  IF(nQbcCh.EQ.0) THEN
     READ(10,*)
  ELSE
     IF (homogene.EQ.1) THEN
        DO i=1,nQbcCh
           READ (10,*) tQbcCh(i),Qbcrec(i)
        END DO
        Qbc(:,1)=Qbcrec
     ELSE 
        IF(qfun.EQ.4) THEN
           DO j=1,nQbcCh
              READ (10,*) (Qbcrec(i),i=1,3) !only 2 compartments possible
              tQbcCh(j)=Qbcrec(1)    !> \param tQbcCh time of a top flux BC
              Qbc(j,1:nQ1)=Qbcrec(2) !> \param Qbc top flux [L T-1]
              Qbc(j,nQ1+1:nQ)=Qbcrec(3)
           ENDDO
        ELSE
           DO j=1,nQbcCh
              READ (10,*) (Qbcrec(i),i=1,nQ+1) !on peut enlever une dimension à qbctot!
              tQbcCh(j)=Qbcrec(1)
              Qbc(j,1:nQ)=Qbcrec(2:nQ+1)
           ENDDO
        ENDIF
     ENDIF
  ENDIF
  IF (nQ.EQ.0) THEN
     WRITE(15,'(/''* No flux B.C.'')',advance='no')
  ELSE
     WRITE(15,'(/''* Nodes with flux B.C.:'',i5)',advance='no') nQ
     WRITE(15,'(/''* Number of time imposed flux B.C.:'',i5)',advance='no')nQbcCh
  ENDIF
  !> Irrigators [L/T]  (Kode(i)=-3):
  !> flux will be multiplied by surface of 1 voxel
  READ (10,*) 
  READ (10,*) 
  READ (10,*) 
  READ (10,*) nI
  IF (nI.GT.maxIrrig) STOP 'Input from  < bc.in >  exceeds parameter "maxIrrig". Program terminated.'
  READ (10,*)
  IF(nI.EQ.0)THEN
     READ(10,*)
  ELSE
     homogene=2
     DO i=1,nI
        READ (10,*) xIrrig(i),yIrrig(i),zIrrig(i)
     END DO
  END IF
  READ (10,*)
  READ (10,*) nIbcCh
  IF (nIbcCh.GT.mxBcCh) STOP 'Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.'
  READ (10,*)
  IF(nIBcCh.EQ.0)THEN
     READ(10,*)
  ELSE
     DO i=1,nI
        READ (10,*) (tIbcCh(i,j),Ibc(i,j),j=1,nIbcCh)
     END DO
  END IF
  DO i=1,nI
     IF ((xIrrig(i).GT.xmax).OR.(xIrrig(i).LT.xmin)) STOP 'Irrigator out of soil domain (X). Program terminated.'
     IF ((yIrrig(i).GT.ymax).OR.(yIrrig(i).LT.ymin)) STOP 'Irrigator out of soil domain (Y). Program terminated.'
     IF ((zIrrig(i).GT.zmax).OR.(zIrrig(i).LT.zmin)) STOP 'Irrigator out of soil domain (Z). Program terminated.'
     ixmin=1+((xIrrig(i)-xmin)/dxgrid)
     ixymin=ixmin+((yIrrig(i)-ymin)/dygrid)*nx
     imin=ixymin+(ABS(zIrrig(i)-zmax)/dzgrid)*nl
     IF (ixmin.NE.FLOOR(ixmin)) STOP 'Irrigator position not equal to a node position (X). Program terminated.'
     IF (ixymin.NE.FLOOR(ixymin)) STOP 'Irrigator position not equal to a node position (Y). Program terminated.'
     IF (imin.NE.FLOOR(imin)) STOP 'Irrigator position not equal to a node position (Z). Program terminated.'
     iBCPt(nQ+i)=INT(imin)
  ENDDO
  DO i=1,nI
     Kode(iBCPt(nQ+i))=-3
  ENDDO
  !> Pressure head BC; h(t)-BC [L]  (Kode(i)=+1):
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) hfun,ltop
  READ (10,*)
  IF ((ltop .AND. (qfun.GT.0)) .AND. (hfun.GT.0)) STOP 'Top boundary condition can either be flux or pressure head. Please redefine hfun and/or qfun in <bc.in>. Program terminated.'
  SELECT CASE(hfun)
     CASE(0) !> + Case0: no PH BC
        READ(10,*)
        nh = 0
     CASE (1) !> + Case1: total surface BC
        READ (10,*)
        IF (ltop) THEN
           nh=nx*ny
           DO i=1,nh
              iBCPt(nQ+nI+i)=i
           ENDDO
           DO i=1,nh
              Kode(iBCPt(nQ+nI+i))=+1
           END DO
        ELSE
           IF(qfun.EQ.4) THEN !irrigation on cylinder --> also cylindrical bottom bc 
              dummy = nx*ny
              nh = 0
              ii = nPt-dummy
              DO i=1,dummy
                 ii=ii+1
                 IF(MatNum(ii).EQ.1) THEN
                    nh = nh+1
                    iBCPt(nQ+nI+nh) = ii
                    Kode(iBCPt(nQ+nI+nh)) = +2
                 END IF
              ENDDO
           ELSE
              nh=nx*ny
              ii=nPt-nh+1
              DO i=1,nh
                 iBCPt(nQ+nI+i)=ii
                 ii=ii+1
              END DO
              DO i=1,nh
                 Kode(iBCPt(nQ+nI+i))=+2
              END DO
           END IF
        ENDIF
     CASE(2) !> + Case2: PH BC on partial surface area between xmin and xmax
       READ (10,*) xhmin,xhmax 
       IF ((xhmin.LT.xmin).OR.(xhmax.GT.xmax)) STOP 'xmin or xmax for PH bc in <bc.in> not within soil domain. Program terminated.'
       xhmin=FLOOR(xhmin/dxgrid)*dxgrid
       xhmax=FLOOR(xhmax+dxgrid/dxgrid)*dxgrid !incl. right nodes...
       nh=(xhmax-xhmin)/dxgrid*ny
       IF (nh.GT.maxbdr) STOP 'Input from  < bc.in >  exceeds parameter "maxbdr". Program terminated.'
       aa=ABS(xhmax-xhmin)/dxgrid
       IF (ltop) THEN
          i=1
          DO j=1,ny
             DO k=1,aa
                iBCPt(nQ+nI+i)=(xhmin-xmin)/dxgrid+k+(j-1)*nx
                i=i+1
             END DO
          END DO
          DO i=1,nh
             Kode(iBCPt(nQ+nI+i))=+1
          END DO
       ELSE
          i=1
          DO j=1,ny
             DO k=1,aa
                iBCPt(nQ+nI+i)=(nPt-nx*ny)+((xhmin-xmin)/dxgrid+k+(j-1)*nx)
                i=i+1
             END DO
          END DO
          DO i=1,nh
             Kode(iBCPt(nQ+nI+i))=+2
          END DO 
       END IF
    CASE DEFAULT 
       STOP 'Sorry, hfun in <bc.in> can only take a value of 0, 1 or 2. Program terminated.'
    END SELECT
  READ (10,*)
  READ (10,*) nhbcCh
  IF (nhbcCh.GT.mxBcCh) STOP 'Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.'
  READ (10,*)
  if (nhbcCh.EQ.0) then
     READ (10,*)
  else
     READ (10,*) (thbcCh(i),hbc(i),i=1,nhbcCh)
  end if
  IF (nh.EQ.0) THEN
     WRITE(15,'(/''* No head B.C.'')',advance='no')
  ELSE
     WRITE(15,'(/''* Nodes with head B.C.:'',i5)',advance='no')nh
     WRITE(15,'(/''* Number of time imposed head B.C.:'',i5)',advance='no')nhbcCh
  ENDIF
  !> free drainage-BC  (Kode(i)=-2):
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) lFrDr

  IF (.NOT.ltop .AND. lFrdr) STOP 'Bottom BC is already defined as pressure head -- no free drainage possible'
  if(lFrDr) then
     nFrDr = nx*ny
     ii = nx*ny*nz - nFrDr +1
     DO i=1,nFrDr
        iBCPt(nQ+nI+nh+i)=ii
        ii=ii+1
     END DO
     DO i=1,nFrDr
        Kode(iBCPt(nQ+nI+nh+i))=-2
     END DO
  else
     nfrdr = 0
  end if

  IF (.not. lFrDr) THEN
     WRITE(15,'(/''* No free drainage nodes.'')',advance='no')
  ELSE
     WRITE(15,'(/''* Bottom B.C. is free drainage'')',advance='no')
  ENDIF
!!$  DO i=1,nFrdr
!!$     Kode(iBCPt(nQ+nI+nh+i))=-2
!!$  END DO
  nBCPts=nQ+nI+nh+nFrdr
  READ (10,*)
  READ (10,*)
  READ (10,*)
  !> solute transport BC:
  !> the type of BC to be invoked for solute transport (KodCB>0 first type; KodCB<0 third type) 
  READ (10,*) assCode
  KodCB=0
  DO i=1,nBCPts/2
     KodCB(i)=assCode
  ENDDO
  READ (10,*)
  READ (10,*)
  READ (10,*) nCBnd1, nCBnd2
  IF (nCBnd1.GT.mxBcCh) STOP 'Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.'
  IF (nCBnd2.GT.mxBcCh) STOP 'Input from  < bc.in >  exceeds parameter "mxBcCh". Program terminated.'
  READ (10,*)
  if (nCBnd1.GT.0) then
     READ (10,*) (tCBnd1(i),CBnd1(i),i=1,nCBnd1)
  else
     READ (10,*)
  end if
  READ (10,*)
  if (nCBnd2.GT.0) then
     READ (10,*) (tCBnd2(i),CBnd2(i),i=1,nCBnd2)
  else
     READ (10,*)
  end if
  READ (10,*)
  READ (10,*)
  READ (10,*) tPulse
  IF ((nCBnd1==0).AND.(nCBnd2==0)) THEN
     WRITE(15,'(/''* No solute transport BCs'')',advance='no')
  ELSEIF (nCBnd1.NE.0) THEN
     WRITE(15,'(/''* 1st time dependent solute transport B.C.: '',i5)',advance='no')nCBnd1
  ELSEIF (nCBnd2.NE.0) THEN
     WRITE(15,'(/''* 2nd time dependent solute transport B.C.: '',i5)',advance='no')nCBnd2
  ENDIF
  CLOSE (10)


  ! reads timestepfactor from partrace input
  IF(lPartUp.OR.lSalinity) THEN     
     OPEN (Unit=10,FILE='Input_PARTRACE.pTraceInpV12',STATUS='OLD',ERR=70)
     DO i=1,31
        READ (10,*)
     END DO
     READ (10,*) timestepfactor
     DO i=1,53
        READ (10,*)
     END DO
     READ (10,*) uptakeorder
     CLOSE (10)
  END IF


  !> #### if input file is available, reads probes.in ####
  WRITE (*,'('' ... looking for file  <Probes.in> ...'')')
  OPEN (Unit=10,FILE='in/Probes.in',STATUS='OLD',ERR=111)
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
           firstOK=.FALSE.
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
           PRINT *,'too small time step in probe.in: only ',mxProf,' profiles will be kept'
           GOTO 78
        ENDIF
     ENDDO
  ENDIF
78 CLOSE (10)

  !run first calculations
  CALL NodebyProbe
  DEALLOCATE(node_temp)
  RETURN
10 STOP 'File  < control.in >  not found -- program terminated.'
20 STOP 'File  < grid.in >  not found -- program terminated.'
30 STOP 'Input file for nodal values not found -- program terminated.'
40 STOP 'File  < bc.in >  not found -- program terminated.'
60 STOP 'File  < mesh.in >  not found -- program terminated.'
70 STOP ' File  <Input_PARTRACE.pTraceInpV12>  not found --'
111 WRITE (*,'(/'' File  <Probes.in>  not found: no defined observation locations'')')
END SUBROUTINE Applic
!******************************************************************************
!> ### reads Doussan RWU related inputs ###
SUBROUTINE DouIn
  USE Typedef
  USE DoussanMat, ONLY : ageLr,LrRoot,nLr,ageKh,Khroot,nKh,hx_min,ave,eqDis,stresval1,&
       stresval2,cavitfun,cavitb,cavitc,stresfun
  USE PlntData, ONLY : a_r,a_1,a_2,b_1,b_2,h0,h1,h2,h3
  USE RootData, ONLY: sign_in,PH_crit,delta_h,lSign,lSign_inst,size_buff,vol_root,&
       vol_buff,lGap,g1,g2,lAQPc,AQPh,AQPv,nAQPc,seglen,segsur
  USE ParamData, ONLY : pi,lPartUp,lSalinity
  USE SoluteRootMat, ONLY : nPerm,agePr,PrRoot,nPass,agePs,PsRoot,sorp,l_linSorb,l_freundSorb, theta_R, rho_R
  IMPLICIT NONE
  INTEGER(ap) :: i,i_ave,i_eqDis,i_Gap,i_AQPc
  !open simulation summary file
  OPEN(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
  WRITE(15,'(//''++ Hydraulic Data for roots ++'')',advance='no')
  WRITE(15,'(/''-----------------------------'')',advance='no')
  WRITE (*,'('' ... looking for file  <CondRoot.in> ...'')')
  !> #### reads CondRoot.in ####
  OPEN (Unit=10,FILE='in/CondRoot.in',STATUS='OLD',ERR=3001)
  WRITE(15,'(/''* The Doussan algorithm will be used for the root water uptake (lDou=TRUE)'')',advance='no')
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) (nLr(i),i=1,3)
  READ (10,*)
  !> \param LrRoot radial root hydraulic conductivity [L/T/L]; per root order (must be defined for 3 orders)
  READ (10,*) (ageLr(1,i),LrRoot(1,i),i=1,nLr(1))  !main axes
  READ (10,*) (ageLr(2,i),LrRoot(2,i),i=1,nLr(2))  !secondary axes
  READ (10,*) (ageLr(3,i),LrRoot(3,i),i=1,nLr(3))  !tertiary axes
  READ (10,*)
  READ (10,*)
  READ (10,*) (nKh(i),i=1, 3)
  READ (10,*)
  !> \param KhRoot axial root hydraulic conductance [L4/T/L]; per root order (must be defined for 3 orders)
  READ (10,*) (ageKh(1,i),KhRoot(1,i),i=1,nKh(1))
  READ (10,*) (ageKh(2,i),KhRoot(2,i),i=1,nKh(2))
  READ (10,*) (ageKh(3,i),KhRoot(3,i),i=1,nKh(3))
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) (nPerm(i),i=1, 3)
  READ (10,*)
  !> \param PrRoot radial root solute permeability [L/T]; per root order (must be defined for 3 orders)
  IF(lPartUp.OR.lSalinity)THEN
     READ (10,*) (agePr(1,i),PrRoot(1,i),i=1,nPerm(1))
     READ (10,*) (agePr(2,i),PrRoot(2,i),i=1,nPerm(2))
     READ (10,*) (agePr(3,i),PrRoot(3,i),i=1,nPerm(3))
     READ (10,*)
     READ (10,*)
     READ (10,*) (nPass(i),i=1, 3)
     READ (10,*)
     !> \param PrRoot radial root solute permeability [L/T]; per root order (must be defined for 3 orders)
     READ (10,*) (agePs(1,i),PsRoot(1,i),i=1,nPass(1))
     READ (10,*) (agePs(2,i),PsRoot(2,i),i=1,nPass(2))
     READ (10,*) (agePs(3,i),PsRoot(3,i),i=1,nPass(3))
     READ (10,*)
     READ (10,*)
     READ (10,*) l_linSorb, l_freundSorb
     READ (10,*)
     IF(l_linSorb) THEN
        READ (10,*) sorp(1)
        READ (10,*)
        READ (10,*)
        READ (10,*) theta_R, rho_R
     ELSEIF (l_freundSorb) THEN
        READ (10,*) sorp(1), sorp(2)
        READ (10,*)
        READ (10,*)
        READ (10,*) theta_R, rho_R
     ELSE
        READ(10,*)
        READ(10,*)
        READ(10,*)
        READ(10,*)theta_R, rho_R
     END IF
  ELSE
     DO i=1,18
        READ(10,*)
     END DO
  END IF
  !averaging method
  READ(10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) i_ave, i_eqDis
  IF (i_ave .EQ. 0) THEN
     ave=.FALSE.
  ELSEIF (i_ave .EQ. 1) THEN
     ave=.TRUE.
     WRITE(15,'(/''* Averaging method: treated as one root (ave= TRUE)'')',advance='no')
  ELSE
     STOP 'Set AVERAGE properly in Condroot.in'
  ENDIF
  IF (i_eqDis .EQ. 0) THEN
     eqDis=.FALSE.
     IF (i_ave .EQ. 0) THEN
        WRITE(15,'(/''* No averaging method will be used below voxel scale (ave= FALSE and eqDis= FALSE)'')',advance='no')
     ENDIF
  ELSEIF (i_eqDis .EQ. 1) THEN
     eqDis=.TRUE.
     WRITE(15,'(/''* Averaging method: use equal distance function (eqDis= TRUE)'')',advance='no')
  ELSE
     STOP 'Set eqDis parameter properly in Condroot.in'
  ENDIF
  !> 1. stress functions
  READ (10,*,ERR=20)
  READ (10,*,ERR=20)
  READ (10,*,ERR=20) stresfun
  READ (10,*,ERR=20)
  SELECT CASE(stresfun)
     CASE(0) !> + case0: No plant water stress responses considered
        READ (10,*,ERR=20)
        WRITE(15,'(/''* No plant water stress is considered '')')
     CASE(1) !> + case1: Hydraulic regulation of stomata at a limiting stress value
        READ (10,*,ERR=20) hx_min !> \param hx_min limitin PH_collar
        hx_min=-ABS(hx_min)
        WRITE(15,'(/''* Hydraulic regulation of stomata with a limiting stress value = '',1pE11.4,1pE11.4)',advance='no') &
             hx_min
     CASE(2) !> + case2: Hydraulic regulation of stomata; linear stress function
        READ (10,*,ERR=20) stresval1,stresval2
        WRITE(15,'(/''* A linear stress function will be used with values = '',1pE11.4,1pE11.4)',advance='no') &
             stresval1,stresval2
     CASE(3) !> + case3: Hydraulic regulation of stomata; Tuzet et al. (2003) function
        READ (10,*,ERR=20) stresval1,stresval2
        WRITE(15,'(/''* The Tuzet stress function will be used with values  = '',1pE11.4,1pE11.4)',advance='no') &
             stresval1,stresval2
     CASE(4) !> + case4: Hydraulic and hormonal regulation of stomata; Tardieu et al. (1993) function; add. hormone transport
        lSign=.TRUE.
        READ (10,*,ERR=20) a_r,a_1,a_2,PH_crit,delta_h,sign_in,size_buff
        WRITE(15,'(/''* Tardieu&Davies stomatal conductance will be calculated with signal transport; '',2(1pE11.4,1X))')  &
             a_r,a_1,a_2
        IF(size_buff.GT.0) THEN 
           vol_buff = vol_root*size_buff
        ELSE
           vol_buff = segsur(1)*segsur(1)/seglen(1)/4/pi
        END IF
     CASE(5) !> + case5: Hydraulic and hormonal regulation of stomata; Tardieu et al. (1993) function; NO hormone transport
        lSign_inst=.TRUE.
        READ (10,*,ERR=20) a_r,a_1,a_2,PH_crit,delta_h,sign_in,size_buff
        WRITE(15,'(/''* Tardieu&Davies stomatal conductance will be calculated with instantaneous signaling; '',2(1pE11.4,1X))')  &
             a_r,a_1,a_2
        IF(size_buff.GT.0) THEN 
           vol_buff = vol_root*size_buff
        ELSE
           vol_buff = segsur(1)*segsur(1)/seglen(1)/4/pi
        END IF
     CASE DEFAULT
        STOP 'Sorry, stresfun in <bc.in> can only take values from 0-5. Program terminated.'     
  END SELECT
  READ (10,*)
  READ (10,*)
  READ (10,*,ERR=30) cavitfun !> 2. cavitational effect on root axial conductance
  IF (cavitfun.EQ.1) THEN
     WRITE (*,*) 'Cavitation''s effect on root axial conductance will be considered'
     READ (10,*)
     READ (10,*,ERR=30) cavitb,cavitc
  ELSE
     READ (10,*)
     READ (10,*)
  ENDIF
  READ (10,*)
  READ (10,*)
  READ (10,*,ERR=40) i_Gap !> 3. Air gap / Rhizosphere hydrophobicity effect on radial root conductivity
  IF (i_Gap.EQ.1) THEN
     WRITE (*,*) 'Air gap / rhizosphere hydrophobicity effect on radial conductivity will be considered'
     lGap=.TRUE.
     READ (10,*)
     READ (10,*,ERR=40) g1,g2
  ELSE
     READ (10,*)
     READ (10,*)
  ENDIF
  READ (10,*)
  READ (10,*)
  READ (10,*,ERR=50) i_AQPc, nAQPc !> 4. Aquaporin effect on radial root conductivity
  IF (i_AQPc.EQ.1) THEN
     WRITE (*,*) 'Aquaporin status''s effect on radial conductivity will be considered'
     lAQPc=.TRUE.
     READ (10,*)
     ALLOCATE (AQPh(1:nAQPc))
     ALLOCATE (AQPv(1:nAQPc))
     READ (10,*,ERR=50) (AQPh(i),AQPv(i),i=1,nAQPc)
  ELSE
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*)
  ENDIF
  IF ((ave) .AND. (eqDis)) STOP 'set one root and equidistant approach correctly; may not be both switched on'
  CLOSE (15)
  RETURN
  CLOSE (10)
3001 STOP ' File  <CondRoot.in>  not found --'
20 STOP ' Water stress information not found in <CondRoot.in>. -- program terminated.'
30 STOP ' Cavitation information not found in <CondRoot.in>. -- program terminated.'
40 STOP ' Gap / hydrophobicity information not found in <CondRoot.in>. -- program terminated.'
50 STOP ' AQP behavior information not found in <CondRoot.in>. -- program terminated.'
END SUBROUTINE DouIn
!*******************************************************************************
!> ### Reads solute related inputs ###
SUBROUTINE ChemIn
  USE Typedef
  USE SolData
  USE PlntData, ONLY: cMm,Vmax,xin,spWgt,fk
  USE RootData, ONLY: nUrf,age,Urf
  IMPLICIT NONE
  INTEGER(ap) :: i,k
  OPEN(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
  WRITE(15,'(//''++ Solute Transport Information ++'')',advance='no')
  WRITE(15,'(/''----------------------------------'')')
  OPEN(Unit=10,File='in/chem.in',STATUS='OLD',ERR=10)
  READ (10,*)
  READ (10,*)
  NLevel=1
  READ (10,*) epsi,PeCr
  IF (epsi.LT.0.999_dp) NLevel=2
  READ (10,*)
  READ (10,*)
  DO k=1,nMat
     READ (10,*,ERR=20) (ChPar(i,k),i=1,9)
  END DO
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) CMm,VMax,xin,fk
  READ (10,*)
  READ (10,*)
  READ (10,*,ERR=20) (nUrf(i),i=1,3)
  READ (10,*,ERR=20)
  DO  i=1,3
     READ (10,*,ERR=20) (age(i,k),Urf(i,k),k=1,nUrf(i))
  END DO
  CLOSE (10)
  WRITE (15,'('' -- simulation will include solute transport.'')',advance='no')
  CLOSE (15)
  RETURN
10 STOP 'File  < chem.in >  not found  -- program terminated'
20 STOP 'Insufficient data in  < chem.in >  -- program terminated.'
END SUBROUTINE ChemIn
!****************************************************************************
!> ### reads soil related inputs ###
SUBROUTINE SoilIn
  USE Typedef
  USE SolData
  USE WatFun
  USE RhizoData, ONLY: rhizoModel, lRhizo
  IMPLICIT NONE
  REAL(sp) ::alh,h1,hn,htab1,htabn
  INTEGER(ap) ::i,k
  
  OPEN(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
  WRITE(15,'(//''++ Soil Information ++'')',advance='no')
  WRITE(15,'(/''----------------------'')',advance='no')
  !> Check which  Rhizosphere model is considered.
  IF (lRhizo) THEN
          OPEN(Unit=11,FILE='in/Rhizo.in',STATUS='OLD',ERR=60)
          READ(11,*)
          READ(11,*)
         READ(11,*)
         READ(11,*,ERR=60) rhizoModel
         CLOSE(11)
  ENDIF
  !> #### reads soil.in ####
  OPEN(Unit=10,FILE='in/soil.in',STATUS='OLD',ERR=10)
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*,ERR=30)nMat,h1,hN,nTab
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
  DO k=1,nMat
      READ (10,*,ERR=20) (par(i,k),i=1,11)
  END DO
  !> Call for rhizosphere properties 
  IF (lRhizo)  CALL RhizoIn
 
  CLOSE (10)
          
  !> define if the soil types are contrasting (e.g. macropore)
  IF (nMat.GT.1) THEN
      DO k=1,nMat
          IF(par(11,k).LT.1E-5_sp) lMacro=.TRUE.
      END DO
  END IF
  CALL IniTab
  !> generate table for interpolation:
  IF (ABS(hTabN-hTab1).GT.0._dp) THEN
      alh1=LOG10(-hTab1)
      dlh=(LOG10(-hTabN)-alh1)/REAL(nTab-1)
      DO i=1,nTab
          alh=alh1+REAL(i-1)*dlh
          hTab(i)=-(10._sp)**alh
      END DO
      DO k=1,nMat
          WRITE(15,'(/''   PH    KH    C    TH '')') 
          DO i=1,nTab
              ConTab(i,k)=FKP(hTab(i),par(:,k),i)
              CapTab(i,k)=FCP(hTab(i),par(:,k))
              TheTab(i,k)=FTh(hTab(i),par(:,k),i)
              WRITE (15,'(4(1X,1pE10.3))') hTab(i),TheTab(i,k),CapTab(i,k),ConTab(i,k)
          END DO
      END DO
             
  ELSE
      WRITE (*,'('' Input error --  hTab1  and  hTab2  in  < >  cannot be equal.'')')
      STOP 'Program terminated.'
  ENDIF
  RETURN
          !> #### if available, reads soiltab.in ####
        10 OPEN(Unit=10,FILE='in/soiltab.in',STATUS='OLD',ERR=40)
           soiltab=.true.
          READ (10,*,ERR=50)
          READ (10,*,ERR=50)
          READ (10,*,ERR=50)
          READ (10,*,ERR=50) nMat, nTab
          READ (10,*,ERR=50) 
          READ (10,*,ERR=50) 
          CALL IniTab
          
          !> tabulated soil hydraulic properties
          !> col1=PH, col2=WC, col3=CH, col4=KH
          DO i=1,nTab
             READ (10,*,ERR=50) hTab(i),(TheTab(i,k),CapTab(i,k),ConTab(i,k),k=1,nMat)
          ENDDO
          READ (10,*,ERR=50)
          READ (10,*,ERR=50) (ssMaxTab(k),k=1,nMat)
          
          alh1=LOG10(ABS(hTab(1)))
          dlh=(LOG10(ABS(hTab(nTab)))-alh1)/REAL(nTab-1)
          WRITE(15,'(/''* Tabulated soil input file'')',advance='no') 
          WRITE(15,'(/''* Number of soil layers:'',1i7)',advance='no') nMat
          WRITE(15,'(/''* Number of tabulations:'',1i7)',advance='no') nTab
          CLOSE (10)
          RETURN
20 STOP 'Insufficient data in  < soil.in >  -- program terminated.'
30 STOP 'Soil input file has been changed (25Feb2011): now, add at the 4th line a 4th number giving the number of tabulations (by default, put 100) !'
40 STOP 'File  < soil.in > or < soiltab.in > not found -- program terminated.'
50 STOP 'Error in < soiltab.in > -- program terminated.'
60 STOP 'ERROR in <rhizo.in -- program terminated'
END SUBROUTINE SoilIn
!***************************************************************************
SUBROUTINE RhizoIn
!> A subroutine that read the rhizosphere parameters and calculate
!RhizoStaticPara
USE Typedef
USE RhizoData, ONLY: rhizoModel, bulkPara, RhizoPara, RhizoSpatialPara, StaticRhizoPara
USE RhizoStat, ONLY: RetCurve, ERR
IMPLICIT NONE
INTEGER :: i,j, itNum=0 , fitLen=15000, maxIter=10000, inx_hcr
REAL(dp) :: thtS, thtR, lambd, hcr, omega, beta, cTot, rhob, rhow
REAL(sp) ::  tht1, tht2, f1, f2, rt, dx, thtMid, sumERR, preERR
REAL(dp) :: lambd_stat, hcr_stat, epsi
REAL(sp), ALLOCATABLE, DIMENSION(:) :: tempTht, tempH, err1
LOGICAL :: converge = .FALSE.
ALLOCATE(tempTht(fitLen), tempH(fitLen))
tempH = -1.
epsi  = 0.0001
OPEN(Unit=11,FILE='in/Rhizo.in',STATUS='OLD',ERR=20)
READ(11,*)
READ(11,*)
READ(11,*)
READ(11,*,ERR=20) rhizoModel
READ(11,*)
READ(11,*)
READ(11,*,ERR=20) (bulkPara(i), i=1,6)
READ(11,*)
READ(11,*)
READ(11,*,ERR=20) (RhizoPara(i), i=1,9)
READ(11,*)
READ(11,*)
READ(11,*,ERR=20) (RhizoSpatialPara(i), i=1,2)
CLOSE(11)
thtR = bulkPara(1); thtS = bulkPara(2); lambd = bulkPara(3); hcr = bulkPara(4)
omega = RhizoPara(1); beta = RhizoPara(2); cTot = RhizoPara(3); rhob = RhizoPara(8); rhow = RhizoPara(9)
!> A procedure to calculate Rhizosphere static parameter. Fitst we use the
!Secant method to generate tempTht(tempH) relation and after we will fit that
!data into the brooks and cory model to find lambd_stat, hcr_stat. The static
!is also used for the initial condition. 
tht1 = thtR
tht2 = thtS
DO i=1,fitLen
    tempH(i) =  - i
    f1 = ERR(tht1, tempH(i))
    f2 = ERR(tht2, tempH(i))
    IF (f1 * f2 .GE. 0.) THEN
        tempTht(i) = thtS
    ELSE 
        IF (f1 < 0.) THEN
            rt = tht1
            dx = tht2 - tht1
        ELSE
            rt = tht2
            dx = tht1 - tht2
        ENDIF
        DO WHILE(.NOT. converge )
            dx = dx*0.5
            thtMid = rt + dx
            f2 = ERR(thtMid, tempH(i))
            IF (f2 .LE. 0.) rt = thtMid
            IF (abs(dx) .LT. epsi) THEN
                tempTht(i) = thtMid
                converge = .TRUE. 
            ENDIF
            itNum = itNum + 1
            IF (itNum .GE. maxIter) THEN
                STOP 'Max Iteration in the secant method reached. see RhizoIn subroutine'
            ENDIF
         ENDDO
    ENDIF
!    write(*,*) tempH(i), tempTht(i), RetCurve(tempH(i),hcr,lambd), itNum, abs(dx)
    converge = .FALSE.
    itNum = 0
    tht1 = thtR
    tht2 = thtS
            
ENDDO
!< Finding hcr of the static Rhizo
DO i=1,fitLen
    IF (tempTht(i) .NE. tempTht(i+1)) THEN
        inx_hcr = i+1
        StaticRhizoPara(2) = (tempH(i+1) + tempH(i))/2.
        GOTO 13
    ENDIF
    IF (i .EQ. fitLen)  STOP 'Did not found hcr Static - program terminated '

ENDDO
13 continue
!< A simple iterative procedure to find lambda rhizo. Note that lambda rhizo is
!always lower then lambda bulk and higher than zero. 
f1=0.
lambd_stat = lambd
converge = .FALSE.
itNum = 0
ALLOCATE(err1(fitLen-inx_hcr))
DO WHILE (.NOT. converge)
    DO j=1,(fitLen-inx_hcr)
        f1 = RetCurve(tempH(j+inx_hcr),hcr,lambd_stat)
        err1(j) = (f1-tempTht(j+inx_hcr))**2.0
    ENDDO
    preERR = sumERR
    sumERR = sum(err1)
    IF (itNum .EQ. 0)  preERR = sumERR + 1.
    IF (sumERR .LT. epsi .OR. preERR .LT. sumERR) THEN
        StaticRhizoPara(1) = lambd_stat
        converge = .TRUE.
    ENDIF
    lambd_stat = lambd_stat - 0.001
    itNum = itNum + 1
ENDDO
IF (sumERR .GE. 1.0) THEN
   STOP 'can not converge in RhizoIn.'
ENDIF
RETURN
20 STOP 'ERROR in <Rhizo.in> -- Program terminated'
END SUBROUTINE RhizoIn
!***************************************************************************
SUBROUTINE MFPC_table
  USE Typedef
  USE SolData
  USE WatFun
  USE RhizoData, ONLY: lRhizo 
  IMPLICIT NONE
  INTEGER(sp):: k,i
  IF(lRhizo) THEN
        WRITE(*,*) 'MFPC_table can not be use in rhizosphere modelling'
        STOP
  ENDIF
  hTab_MFP(1501)=-15010.0_sp
  DO  k=1,nMat
     MFPTab(1501,k)=0.0_sp
  END DO
  DO  i=1,1500
     hTab_MFP(1501-i)=-15010.0_sp+10.0_sp*i
     DO  k=1,nMat
        MFPTab(1501-i,k)=MFPTab(1502-i,k)+(FKP(hTab_MFP(1501-i),par(:,k),i)+FKP(hTab_MFP(1502-i),par(:,k),i))*10.0_sp/2.0_sp
     END DO
  END DO
  RETURN
END SUBROUTINE MFPC_Table
!***************************************************************************
!> ### reads root related inputs ###
SUBROUTINE RootIn(t)
  USE Typedef
  USE tmctrl
  USE ParamData
  USE TempData
  USE RootData
  USE PlntData
  USE GridData
  USE GeoData
  USE ConData
  USE StrData
  USE DoussanMat, ONLY: nplant
  IMPLICIT NONE
  
  REAL(sp) :: maxZ,dumR(maxemg),dumR1(maxemg,mxpnts)=0,dumR2(maxemg,mxpnts)=0
  INTEGER(ap) :: ifg(maxemg)=0,naxtot,iplant,ipl,dumI(maxemg)=0
  INTEGER(ap) :: iestbl,i,j,k,irec,igrow,ifive,timeRT,linlim
  REAL(sp) :: t
  CHARACTER infile*12,ch*1
    
  !> #### if RootTyp is used, reads param.txt ####
  IF (lrrt) THEN	
     OPEN (UNIT=10,FILE='in/param.txt',STATUS='OLD',ERR=50)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) timeRT
     !> no multiple roots with RootTyp
     nplant=1
     xplant(1)=0
     yplant(1)=0
     CALL RunRootTip(timeRT)
     t=timeRT !implicit change of type INT->real
     WRITE(15,'(/''RootTyp was used for generating the root system'')')
     WRITE(15,'(/''Initial Time for simulation is '',1pE11.4)',advance='no') t
     !> #### reads RootSys ####
  ELSEIF (lrrs) THEN	!Read RootSys and root.in --> root growth with Somma
     WRITE (*,*) '... Reading <in/RootSys> ...'
     infile='in/RootSys'
     OPEN (UNIT=10,FILE=infile,STATUS='OLD',ERR=60)
     WRITE(15,'(/''* Root system imported from file '',1A12)',advance='no') infile
     READ (10,*,ERR=40)
     !> 1. general parameters
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
        WRITE(15,'(/''* seed locations for plant '',1I5, '' is '', 2F9.3)',advance='no') iplant,xplant(iplant),yplant(iplant)
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
     !> 2. segment information
     DO i=1,nrec
        READ (10,*,ERR=40) irec,xs(irec),ys(irec),zs(irec),irecpr(irec),ordseg(irec),&
             ibrseg(irec),seglen(irec),segsur(irec),segmas(irec)
        READ (10,*,ERR=40) timorg(irec)
        !Cross Section of each segment - later needed for flow velocities & particle tracker
        segrad(irec) = segsur(irec)/2._dp/pi/seglen(irec)
        crossSectionSeg(irec) = segrad(irec)**2*pi
        segvol(irec) =  segrad(irec)**2*pi*seglen(irec) !root volume per root segment
        vol_root = vol_root +  segvol(irec) !total root volume
     END DO
     IF (MAXVAL(ordseg).EQ.13) THEN		!Recognizes old roots created by RootTyp (Couvreur nov 2010)
        maizeroottyp=.TRUE.
     ELSEIF (MAXVAL(ordseg).EQ.5) THEN
        loliumroottyp=.TRUE.
     ELSEIF (MAXVAL(ordseg).EQ.20) THEN
        wheatroottyp=.TRUE.
     ENDIF
     ibrseg(0)=0.
     segsur(0) =0.
     IF ((zs(1).GT.zgrid(1)).OR.(zs(1).LT.zgrid(npt))) THEN
        WRITE (*,'(//'' Inconsistency -- Position of first root segment not within the spatial domain'')')
        WRITE (*,'('' as defined in  < nodes.in >.''/)')
        STOP 'Program terminated.'
     ENDIF
     IF (.NOT.(continu)) THEN		
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
     !> 3. tip information
8    READ (10,*,ERR=40) igrow,xg(igrow),yg(igrow),zg(igrow),irecsg(igrow),ordgrw(igrow),ibrgrw(igrow),&
          brlgth(igrow), iaxis(igrow)
     READ (10,*,ERR=40) ovrtime(igrow),nestbl(igrow)
     IF (nestbl(igrow).GT.0) THEN
        ifive=0
81      linlim=MIN(ifive+5,nestbl(igrow))
        READ (10,*,ERR=40) (timest(igrow,iestbl),iestbl=ifive+1,linlim)
        ifive=ifive+5
        IF (linlim.LT.nestbl(igrow)) GOTO 81
     ENDIF
     maxZ=MAXVAL(ZGrid)
     IF (igrow.LT.ngrow) GOTO 8
     zg(1:ngrow)=zg(1:ngrow)-maxZ !??? maxZ?
     CLOSE (10)
  END IF

  !> #### if Somma growth, reads root.in  ####
  IF (lSomma_growth) THEN
     OPEN(UNIT=15,FILE='out/simul_sum.out',STATUS='OLD',POSITION='APPEND')
     WRITE(15,'(//''++ Root System ++'')',advance='no')
     WRITE(15,'(/''-----------------'')',advance='no')
     WRITE (*,'('' ... looking for file  <Root.in> ...'')')
     OPEN (UNIT=10,FILE='in/root.in',STATUS='OLD',ERR=10)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) naxemg
     IF (naxemg.GT.maxemg) STOP 'Input from  < root.in >  exceeds parameter "maxemg". Program terminated.'
     READ (10,*,ERR=20)
     DO i=1,naxemg
        READ (10,*,ERR=20) tnewax(i),nnewax(i), dumR(i)
     END DO
     naxtot=naxes
     ifg(1:naxtot)=1 !initial axes aways belong to axis group #1
     DO i=1,naxemg
        ifg(naxtot+1:naxtot+nnewax(i))=i
        naxtot=naxtot+nnewax(i)
     END DO
     DO i=1,naxtot
        inaxs(i)=dumR(ifg(i))/180._sp*pi
     END DO
     dumR = 0.
     IF (naxtot.GT.maxemg) STOP 'Input from  < root.in >  exceeds parameter "maxemg". Program terminated.'
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) (dumR(i),i=1,naxemg)
     DO i=1,naxtot
        geoaxs(i)=dumR(ifg(i))
     END DO
     dumR = 0.
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) (dumI(i),i=1,naxemg)
     DO i=1,naxtot
        nangax(i)=dumI(ifg(i))
     END DO
     READ (10,*,ERR=20)
     DO  i=1,naxemg
        READ (10,*,ERR=20) (dumR1(i,j),dumR2(i,j), j=1, dumI(i))
     END DO
     DO i=1,naxtot
        tempax(i,1:nangax(i)) = dumR1(ifg(i),1:nangax(i))
        angaxs(i,1:nangax(i)) = dumR2(ifg(i),1:nangax(i))
     END DO
     dumI = 0.
     dumR1 = 0
     dumR2 = 0
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) geolat
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) nanglt
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) (templt(i),anglat(i),i=1,nanglt)
     DO j=1,nanglt
        anglat(j)=anglat(j)/180._sp*pi
     END DO
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) norder
     IF (norder.GT.maxord) STOP 'Input from  < root.in >  exceeds parameter "maxord". Program terminated.'
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) (nVch(i),i=1,norder)
     READ (10,*,ERR=20)
     DO i=1,norder
        READ (10,*,ERR=20) (ageVch(i,j),Vch(i,j),j=1,nVch(i))
     END DO
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) SpWgt
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) (nMPLch(i),i=1,norder)
     READ (10,*,ERR=20)
     DO  i=1,norder
        READ (10,*,ERR=20) (sMPLch(i,j),MPLch(i,j),j=1,nMPLch(i))
     END DO
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) l_conduc
     READ (10,*,ERR=20)	 
     READ (10,*,ERR=20)	condMP
     READ (10,*,ERR=20)	 
     READ (10,*,ERR=20) (strsen(i),i=1,norder)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)     
     READ (10,*,ERR=20) (rdmang(i),i=1,norder)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) tempermin,topt,tempermax	
     trange=tempermax-tempermin
     tmid=(tempermin+tempermax)/2._sp
     IF (topt.LT.tmid) THEN
        expo=LOG(.5_sp)/LOG((topt-tempermin)/trange)
     ELSE
        expo=LOG(.5_sp)/LOG((topt-tempermax)/(-trange))
     ENDIF
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) ltoxi
     READ (10,*,ERR=20)
     IF (ltoxi) THEN
        READ (10,*,ERR=20) cmin,coptmi,coptma,cmax
     ELSE
        READ (10,*,ERR=20)
     ENDIF
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) (brlmax(i),i=1,norder)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) (brspac(i),i=1,norder-1)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) (brnang(i),i=1,norder-1)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) (dtbrch(i),i=1,norder-1)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20)
     READ (10,*,ERR=20) l_secrad
     READ (10,*,ERR=20)
     IF(l_secrad) READ (10,*,ERR=20) (f_rad(i),i=1,norder)
     CLOSE (10)
     DO i=1,norder-1
        brnang(i)=brnang(i)/180._sp*pi
        rdmang(i)=rdmang(i)/180._sp*pi
     END DO
     rdmang(norder)=rdmang(norder)/180._sp*pi
   ENDIF
   
   IF (.NOT.lretry) THEN
      WRITE (*,'(///'' Simulation starts at time = '',F8.3,''.'')') t
      WRITE (15,'(/''The root system consists now of '',I5,'' branch(es)[including '',I3,'' axis(es)]'')')nbr,naxes
      WRITE (15,'('' or of a total of '',I5,'' segment(s) and '',I5,'' growing branch tip(s).'')')nrec,ngrow
      WRITE (*,'(/'' Total root mass is '',F9.3,'', total shoot mass '',F9.3,''.'')')mroot,mshoot
      WRITE (*,'('' Leaf area is '',F9.3,''.'')')LA
      WRITE (*,'(///70X,''<ENTER>''/)') 
   ENDIF
   RETURN
10 STOP 'Data inconsistency in  < root.in > not found -- program terminated.'
20 STOP 'Data inconsistency in  < root.in >  -- program terminated.'
  !    30 STOP 'Root input file does not exist -- program terminated.'
40 STOP 'Data error in root system file  -- program terminated.'
50 STOP 'File param.txt not found --program terminated'
60 STOP 'File RootSys not found --program terminated'
END SUBROUTINE RootIn
!**************************************************************************************************
!> ### plant / shoot relevant inputs ###
SUBROUTINE PlntIn(t)
  USE typedef
  USE ParamData, ONLY: lChem,Pi,lretry,last_out
  USE PlntData
  USE RootData,ONLY : lDou,lCou,lFed,lSomma_growth,lno_RWU,lCalloc,lno_Archi,ltemp,ltoxi
  USE DoussanMat, ONLY: stresfun
  USE tmctrl
  IMPLICIT NONE
  INTEGER(ap) :: i,j,err,funBC,aa,bb,n,nt
  REAL(sp) :: t,BCd(mxbcch),typed(mxbcch)
  REAL(sp), ALLOCATABLE, DIMENSION(:) :: tBCrt1,dummy
  INTEGER(ap), ALLOCATABLE, DIMENSION(:) :: typeBCrt
  
  !read BCRoot.in
  IF  (.NOT.lno_RWU) THEN		!BCroot.in now read as long as the RWU process and no Carbon allocation are simulated, Tpot.in erased 
     ! file with boundary conditions given for root in therms of PH or of flux
     WRITE (*,'('' ... looking for file  <BCroot.in> ...'')')
     !> #### reads BCroot.in ####
     OPEN (Unit=10,FILE='in/BCroot.in',STATUS='OLD',ERR=3001)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) funBC 
     READ (10,*)
     READ (10,*) nBCr
	 
     IF (nBCr.GT.mxBcCh) THEN
        PRINT *,' the number of boundary conditions for root exceeded the maximum. Simulation stopped.'
        STOP
     ENDIF
     READ (10,*)
     IF (funBC.LE.2) THEN !> free format or constant value
        READ (10,*) (tBCr(i),typeBCr(i),BCroot(i),i=1,nBCr)
        IF (lFed) THEN
           DO i=1,nBCr
              IF (typeBCr(i).NE.2) STOP 'With Feddes RWU model, the upper BC at the root collar should be of flux type'
           ENDDO
        ENDIF
     ELSE
        ALLOCATE (tBCrt1(1:nBCr+1))
        ALLOCATE (dummy(1:nBCr+1))
        ALLOCATE (typeBCrt(1:nBCr+1))
        READ (10,*) (tBCrt1(i),typeBCrt(i),dummy(i),i=1,nBCr)
        IF (lFed) THEN
           DO i=1,nBCr
              IF (typeBCrt(i).NE.2) STOP 'With Feddes RWU model, the upper BC at the root collar should be of flux type'
           ENDDO
        ENDIF
        !> if runtime is longer than last tBCrin, the last BCr will be used
        IF (tBCrt1(nBCr).LT.tMax) THEN
           nBCr = nBCr+1
           tBCrt1(nBCr) = tMax
           dummy(nBCr) = dummy(nBCr-1)
        END IF
        !> if transpiration is not constant, interpolate daily values
        nt = 0
        DO i=2,nBCr
           aa = tBCrt1(i-1)
           bb = tBCrt1(i)
           DO j = 1,bb-aa
              nt = nt+1
              BCd(nt) = dummy(i-1) + j*(dummy(i)-dummy(i-1))/(tBCrt1(i)-tBCrt1(i-1))
              typed(nt) = typeBCrt(i-1)
           END DO
        END DO
        
        n = 1
        DO i=1,nt
           DO j=1,25 !> sub-daily interpolation steps (25/day)
              tBCr(n) = tBCrt1(1)+(i-1)+REAL(j)/25
              typeBCr(n) = typed(i)
              IF (funBC .EQ. 3) THEN
                 BCroot(n) = BCd(i)*(SIN(2*pi*(tBCr(n)-0.25))+1) 
                 n = n+1
              ELSE IF (funBC .EQ. 4) THEN
                 BCroot(n) = BCd(i)*pi*SIN((tBCr(n)-0.25)*2*pi)
                 IF (BCroot(n).LT.0) BCroot(n)=0.0
                 n = n+1
              ELSE
                 WRITE (*,'('' ... ERROR in input file BCroot.in ... (funBC id#)'')')
                 STOP
              END IF
           END DO
        END DO
        nBCr = n-1
     ENDIF
     WRITE (*,'(/'' BC-data for root found.'')')
     WRITE (*,'(///70X,''<ENTER>''/)')
10   CLOSE (10)
     
     IF (lno_Archi) THEN
        IF (lretry) THEN
           t=tOut(last_out)
           WRITE (*,'(///'' Simulation continues at time = '',F8.3,''.'')') t
        ELSE
           t=tBCr(1) !First time step is the first time of BCroot.in when RootSys or RootTyp architectures are not used 
        ENDIF
     ENDIF
     IF ((tBCR(1).LT.t).AND.(.NOT.lretry)) THEN
        PRINT *,'tBCR',tBCR(1),tBCR(2),tBCR(3),tBCR(4),tBCR(5),'t',t
        WRITE(*,'(/'' initial time for root BC lower than initial simulation time'')')
        STOP
     ENDIF
  ENDIF
  !> #### if shoot growth is considered, reads plant.in ####
  IF (lCalloc) THEN
     OPEN (Unit=10,FILE='in/plant.in',STATUS='OLD',ERR=11)
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
     !> no-stress water use efficiency (dry mass gained per water mass transpired)
     !> as f(time):
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) ntW
     IF (ntW.GT.mxBcCh) GOTO 99
     READ (10,*)
     READ (10,*) (tW(i),Wc(i),i=1,ntW)
     !> no-stress root/shoot ratio as a f(time):
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) ntRSR
     IF (ntRSR.GT.mxBcCh) GOTO 99
     READ (10,*)
     READ (10,*) (tRSR(i),RSRc(i),i=1,ntRSR)
     !> transpiration reduction factor as a function of relative stress due soil strength:
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) nfTpLA
     IF (nfTpLA.GT.mxBcCh) GOTO 99
     READ (10,*)
     READ (10,*) (sfTpLA(i),fTpLAc(i),i=1,nfTpLA)
     !> water use efficiency factor as a function of relative stress due soil strength:
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) nsfW
     IF (nsfW.GT.mxBcCh) GOTO 99
     READ (10,*)
     READ (10,*) (sfW(i),fWc(i),i=1,nsfW)
     !> root/shoot-ratio factor as a function of relative stress due soil strength:
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) nsfRSR
     IF (nsfRSR.GT.mxBcCh) GOTO 99
     READ (10,*)
     READ (10,*) (sfRSR(i),fRSRc(i),i=1,nsfRSR)
     !> relative stress as a function of soil strength:
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) ns
     IF (ns.GT.mxBcCh) GOTO 99
     READ (10,*)
     READ (10,*) (sc(i),rsc(i),i=1,ns)
     !> transpiration reduction factor as a function of relative stress due solute concentration:
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) ncTpLA
     IF (ncTpLA.GT.mxBcCh) GOTO 99
     READ (10,*)
     READ (10,*) (scTpLA(i),cTpLAc(i),i=1,ncTpLA)
     !> water use efficiency factor as a function of relative stress due solute concentration:
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) nscW
     IF (nscW.GT.mxBcCh) GOTO 99
     READ (10,*)
     READ (10,*) (scW(i),cWc(i),i=1,nscW)
     !> root/shoot-ratio factor as a function of relative stress due solute concentration:
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) nscRSR
     IF (nscRSR.GT.mxBcCh) GOTO 99
     READ (10,*)
     READ (10,*) (scRSR(i),cRSRc(i),i=1,nscRSR)
     !> relative stress as a function of solute concentration:
     READ (10,*)
     READ (10,*)
     READ (10,*)
     READ (10,*) ncnc
     IF (ncnc.GT.mxBcCh) GOTO 99
     READ (10,*)
     READ (10,*) (cncp(i),rscnc(i),i=1,ncnc)
     READ (10,*)
     READ (10,*)
     READ (10,*)
     !> leaf area per dry shoot mass as a function of time:
     READ (10,*) ntLA
     IF (ntLA.GT.mxBcCh) GOTO 99
     READ (10,*)
     READ (10,*) (tLA(i),LAc(i),i=1,ntLA)
     CLOSE (10)
END IF
  !> Options declaration; which processes are modeled
  IF (lSomma_growth) THEN
     IF (.NOT.lno_RWU) THEN
        IF (lChem) THEN
           WRITE (*,'(/'' Interaction between root water - solute uptake / root growth'')')
           WRITE (*,'(''     and soil water content / soil strength will be simulated.'')')
        ELSE
           WRITE (*,'(/'' Interaction between root water uptake / root growth'')')
           WRITE (*,'(''     and soil water content / soil strength will be simulated.'')')
        END IF
        IF (lCalloc) THEN
           WRITE (*,'(/'' Shoot growth and assimilate distribution will be considered. '')')
           WRITE (*,'(///70X,''<ENTER>''/)')
        ELSE
           WRITE (*,'(/'' Shoot growth and assimilate distribution will not be considered. '')')
           WRITE (*,'(///70X,''<ENTER>''/)')
        ENDIF
     ELSE
        WRITE (*,'('' -- will simulate root growth as affected by soil strength only.'')')
        WRITE (*,'(/'' Water and solute uptake, shoot growth and assimilate      '')')
        WRITE (*,'(''     distribution will not be considered.'')')
     ENDIF
     IF (ltoxi) THEN
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
     ENDIF
  ENDIF
  IF (lFed) THEN
     IF (stresfun.EQ.2) THEN
        IF (p50.LT.9.E+20_dp) THEN
           WRITE (*,'(/'' The van Genuchten (1978) expression will be used as water'')')
           WRITE      (*,'(''     extraction function. '')')
        ELSE
           WRITE (*,'(/'' The van Genuchten (1978) expression will be used for the water'')')
           WRITE (*,'(''     extraction function, but osmotic potential effects will be'')')
           WRITE (*,'(''     neglected. '')')
        ENDIF
     ELSEIF (stresfun.EQ.1) THEN
        WRITE (*,'(//'' The Feddes et al. (1978) expression will be used for the water'')')
        WRITE (*,'(''     extraction function. '')')
     ELSE
        WRITE (*,'(//'' The water extraction function will be directly proportional to the relative'')')
        WRITE (*,'(''     root length density. '')')
     ENDIF
  ELSEIF (lDou) THEN
     WRITE (*,'(//'' The Doussan et al. (1998) expression will be used as water extraction function'')')
  ELSEIF (lCou) THEN
     WRITE (*,'(//'' The Couvreur et al. (2012) expression will be used as water extraction function'')')
  ENDIF
  WRITE (*,'(///70X,''<ENTER>''/)')
  RETURN
99 STOP 'Input from  < plant.in >  exceeds parameter "mxBcCh". Program terminated.'
3001 STOP ' File  <BCroot.in>  not found --'
11 STOP ' File  < plant.in >  not found --'
END SUBROUTINE PlntIn
!***************************************************************************************
!> ### if temperature effect is considered, reads temp.in ###
SUBROUTINE TempIn
  USE TempData
  USE RootData, ONLY: ltemp,lSign
  IMPLICIT NONE
  INTEGER(ap) :: ii, jj

  OPEN (Unit=10,FILE='in/temp.in',STATUS='OLD',ERR=100)
  !> soil temperature as a function of time and depth
  !> (piecewise linear function):
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) nz_tempS,nt_tempS
  IF (nz_tempS.GT.mxdpth) STOP 'Input from  < temp.in >  exceeds parameter"mxdpth". Program terminated.'
  READ (10,*)
  READ (10,*)
  READ (10,*) (depth(ii),ii=1,nz_tempS)
  IF (ABS(depth(1)).GT.1.E-30_dp) THEN
     WRITE (*,'(//'' Inconsistency - Second depth value in  < temp.in >  corresponds to soil surface'')')
     WRITE (*,'('' and must be zero.''/)')
     STOP 'Program terminated.'
  ENDIF
  READ (10,*)
  READ (10,*)
  DO ii=1, nt_tempS
     READ (10,*) time_S(ii),(temtim(ii,jj),jj=1,nz_tempS)
  END DO
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) nt_tempA,nt_presA,nt_presD 
  READ (10,*)
  READ (10,*)
  READ (10,*) (time_TA(ii),T_atm(ii),ii=1,nt_TempA)
  READ (10,*)
  READ (10,*)
  READ (10,*) (time_PA(ii),P_atm(ii),ii=1,nt_presA)
  READ (10,*)
  READ (10,*)
  READ (10,*) (time_PD(ii),P_diff(ii),ii=1,nt_presD)
  CLOSE(10)
  IF(ltemp) THEN
     WRITE (*,'(/'' Temperature data found.'')')
     WRITE (*,'('' Simulation will include temperature influence for root growth.'')')
     WRITE (*,'(///70X,''<ENTER>''/)')
  ELSEIF(lSign) THEN
     WRITE (*,'(/'' Temperature data found.'')')
     WRITE (*,'('' Simulation will include signaling.'')')
     WRITE (*,'(///70X,''<ENTER>''/)')
  END IF
  RETURN
100 STOP ' File  < temp.in >  not found -- needed if temperature influence on root growth or signaling should be modeled!'
END SUBROUTINE TempIn
!********************************************************************
!> ### defines the number of nodes located in one observation plane / probe window ###
SUBROUTINE NodebyProbe
  USE Typedef
  USE ObsData
  USE GridData
  USE DomData
  IMPLICIT NONE 
  INTEGER(ap) :: ip,delta,in,imin,i,iz,ix,iy
  
  !> observation planes
  DO ip=1,npr
     IF (Pt(ip)==3) THEN !horiz. plane perp. to Z
        nodebyPr(ip)=nx*ny
        delta=INT((ABS(CrP(ip)-zmax))/dzgrid)
        imin=nodebyPr(ip)*delta+1
        in=imin
        i=0
        DO i=1,nodebyPr(ip)
           NodePr(ip,i)=in
           in=in+1
        ENDDO
     ELSEIF (Pt(ip)==2) THEN !vert. plane perp. to Y
        nodebyPr(ip)=nx*nz
        delta=INT((ABS(CrP(ip)-ymin))/dygrid)
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
        delta=INT((ABS(CrP(ip)-xmin))/dxgrid)
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
        WRITE(*,*)'Wrong plane direction in Probes.in'
        STOP
     ENDIF
     WRITE(*,*)'plane',ip,'has ',nodebyPr(ip),'nodes.'
     WRITE(*,*) Pt(ip),ip,nodebyPr(ip),delta,imin,in
     IF (nodebyPr(ip)>1000) THEN
        WRITE(*,*) 'too many elements in one observation plane (>1000)'
        STOP
     ENDIF
  ENDDO
  RETURN
END SUBROUTINE NodebyProbe
!*********************************************************************
!> ### Defines list with elements that are of mixed material ###
SUBROUTINE ListMacro
  USE typedef
  USE GridData
  USE SolData, ONLY: MatNum,l_elmMacro,par,nMat
  IMPLICIT NONE
  INTEGER(ap) :: iElm,matSum,i
  INTEGER(ap) :: NeighElm(1:6),cif,k,imat
  cif=0  !number of mixed material voxels
  ! identify macropore material
  DO k=1,nMat
     IF(par(11,k).LT.1E-5_sp) imat=k
  END DO
  ! identify mixed material elements
  DO iElm=1,nElm
     IF(ANY(MatNum(elmnod(:,iElm)).EQ.imat)) THEN
        l_elmMacro(iElm) = .TRUE.
        cif=cif+1
        CALL ElemNeighb(iElm,NeighElm)
        DO i = 1,6
           IF(NeighElm(i).GT.0)THEN
              IF(ANY(MatNum(elmnod(:,NeighElm(i))).EQ.imat)) THEN
                 MacroList(iElm,i)=0
              ELSE
                 MacroList(iElm,i)=NeighElm(i)
                 n_neigh(iElm) = n_neigh(iElm)+1
              END IF
           END IF
        END DO
     END IF
  END DO
  !print*, 'number of mixed material voxels', cif
END SUBROUTINE ListMacro
!*********************************************************************
