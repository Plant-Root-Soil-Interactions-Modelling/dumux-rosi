! ==============================================================================
! Source file PLANT GROWTH |||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
! current relative stress due to soil strength and solute concentration:
SUBROUTINE Stress(rs,concrs)
  USE ParamData, ONLY: lChem
  USE PlntData
  USE RootData, ONLY: sAvg,cAvg
  USE StrData
  IMPLICIT NONE
  REAL(sp), INTENT(out) :: rs,concrs
  INTEGER(ap) ::  k
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
11   k=k+1
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
SUBROUTINE SetTp(t,rs,concrs,ipl)
  USE ParamData, ONLY: lChem
  USE PlntData
  USE RootData, ONLY: concol
  USE DoussanMat, ONLY: PHr,BCr_usr
  USE TempData
  IMPLICIT NONE
  INTEGER(ap) :: ifc,k,ipl,aa
  REAL(sp), INTENT(in) :: t,rs,concrs
  REAL(sp) :: fcTpLA,fTpLA,fhTpLA
  
  CALL Atmosphere(t)

  ! calculate current TpLA-value:
  ifc=1
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
302  k=k-1
     IF (k.EQ.0) THEN
        ! relative stress smaller than smallest rs for which fTpLA is specified:
        fTpLA=1.0_dp+rs/sfTpLA(1)*(fTpLAc(1)-1._dp)
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
        ! fcTpLA=cTpLAc(ncTpLA)
     ELSE
        k=ncTpLA
303     k=k-1
        IF (k.EQ.0) THEN
           ! relative stress smaller than smallest concrs for which fTpLA is specified:
           fcTpLA=1.0_dp+concrs/scTpLA(1)*(cTpLAc(1)-1._dp)
        ELSE
           IF (scTpLA(k).GT.concrs) GOTO 303
           fcTpLA=cTpLAc(k)+(concrs-scTpLA(k))/(scTpLA(k+1)-scTpLA(k))*(cTpLAc(k+1)-cTpLAc(k))
        ENDIF
     ENDIF
  ELSE
     fcTpLA=1.0_dp
  ENDIF
  Tpot=TpLA*LA
  BCr_usr(ipl)=-ABS(Tpot*MIN(fTpLA,fcTpLA,fhTpLA)) !root BC is negative, while transpiration is positive... 
  RETURN
END SUBROUTINE SetTp
!*********************************************************************************************
! integrate uptake over spatial domain:
SUBROUTINE ActTrs
  USE GridData
  USE PlntData, ONLY : Tact
  IMPLICIT NONE
  REAL(sp) :: sum,sume
  INTEGER(ap) :: corner(8),ie,ic
  ! calculate actual overall trandpiration rate from soil domain:
  sum=0.0_dp
  DO  iE=1,nElm!,2
     ! assign cuboid corner nodes:
     corner(1)=elmnod(1,iE)
     corner(2)=elmnod(2,iE)
     corner(3)=elmnod(3,iE)
     corner(4)=elmnod(4,iE)
     corner(5)=elmnod(5,iE)
     corner(6)=elmnod(6,iE)
     corner(7)=elmnod(7,iE)
     corner(8)=elmnod(8,iE)
     ! add up average cuboid sink terms, integrate over volume:
     sumE=0.0_dp
     DO ic=1,8
        sumE=sumE+sink(corner(ic))
     END DO
     sum=sum+sumE
  END DO
  Tact=sum*dxGrid*dyGrid*dzGrid/8._dp
  RETURN
END SUBROUTINE ActTrs
!********************************************************************************
! current water use efficiency:
SUBROUTINE Effncy(t,rs,concrs,W)
  USE ParamData, ONLY: lChem
  USE PlntData
  IMPLICIT NONE
  REAL(sp) ::  t,rs,concrs,W,fcw,fw
  INTEGER(ap) ::  ifc,k
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
  ! reduction factor, fcW, corredponding to current relative stress, concrs;
  ! by definition, fcW = 1.0 at zero stress (concrs = 0.0)
  IF (lChem) THEN
     IF (concrs.GE.scW(nscW)) THEN
        ! relative stress greater than greatest concrs for which fcW is dpecified:
        fcW=cWc(nscW)
     ELSE
        k=nscW
3       k=k-1
        IF (k.EQ.0) THEN
           ! relative stress smaller than smallest rs for which fW is dpecified:
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
SUBROUTINE Ratio(t,rs,concrs,RSR)
  USE ParamData, ONLY: lChem
  USE PlntData
  IMPLICIT NONE
  REAL (sp) :: t,rs,concrs,RSR,fcrsr,frsr
  INTEGER(ap) ::  ifc,k
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
  REAL(sp) ::  LAmshv,t 
  INTEGER(ap) ::  ifc
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
!*****************************************************************





