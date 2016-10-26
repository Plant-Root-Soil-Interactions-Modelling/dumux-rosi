! ==============================================================================
! Source file ENVIRONMENTAL FUNCTIONS ||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
! current spatial soil strength distribution:
SUBROUTINE Solstr(h)
  USE typedef
  USE GridData
  USE SolData
  USE WatFun
  USE StrData
  USE RhizoData, ONLY : bulkPara,  lRhizo
  IMPLICIT NONE
  INTEGER(ap) :: i,m
  REAL(sp), INTENT(in) ::  h(nPt)
  REAL(sp) :: sat,dummy

  ! s = soil strength
  ! par(11,M) = sMax = max. soil strength for each soil [=f(rhoB,sand content)]
  DO i=1,nPt
     M=MatNum(i)
!!!!* assign current soil strength values to each node:
     IF (soiltab) THEN	!(Couvreur nov 2011)
        s(i)=((1._sp-(Fth_soiltab(h(i),M)-TheTab(nTab,M))/(TheTab(1,M)-TheTab(nTab,M)))**3)*ssMaxTab(M)
     ELSE
        IF (.NOT. lRhizo) THEN
			sat=(FTh(h(i),par(:,M),i)-par(2,M))/(par(3,M)-par(2,M)) 
		ELSE
		    sat=(FTh(h(i),par(:,M),i)-bulkPara(1))/(bulkPara(2)-bulkPara(1)) 
		ENDIF
        IF (sat .LE. 0._sp) sat=0.001_sp
        dummy=0.35*LOG10(abs(h(i))*0.1*sat)+0.93*par(11,M)+1.26 !kPa
        s(i)=10**dummy*0.001 !MPa
    ENDIF
!!!!* current local impedance factor due to soil strength:
    !local maximum penetrometer resistance, Bengough(2011)
    ssmax = 4._sp-2.33_sp*h(i)*0.0001_sp !MPa 
    IF(ssmax.LE.0._sp) ssmax = s(i)
    imps(i) = 1._sp-s(i)/ssmax

!     IF (s(i).GT.simp) THEN
!        imps(i)=0.0_sp
!     ELSE
!        imps(i)=1._sp-s(i)/simp
!     ENDIF

  END DO
  RETURN
END SUBROUTINE Solstr
!***********************************************************
! current spatial temperature distribution in the soil:
SUBROUTINE Temper(t)
  USE ParamData, ONLY: pi
  USE TempData
  USE GridData
  IMPLICIT NONE
  REAL(sp) :: val1,val2,depfac,t,timfac
  INTEGER(ap) :: i,it,it1,it2,iz1,iz,iz2
!!!!* find the current time interval:
  IF (t.GE.time_S(nt_tempS)) THEN
     it1=nt_tempS
  ELSE
     it=nt_tempS
1    it=it-1
     IF ((t.LT.time_S(it)).AND.(it.GT.1)) GOTO 1
     it1=it
     it2=it+1
     timfac=(t-time_S(it1))/(time_S(it2)-time_S(it1))
  ENDIF
  DO  i=1,nPt
!!!!* assign current temperature values to each node --
!!!!* find the appropriate depth interval:
     IF (zGrid(i).LE.zGrid(1)-depth(nz_tempS)) THEN
        iz1=nz_tempS
     ELSE
        iz=nz_tempS
11      iz=iz-1
        IF ((zGrid(i).GT.zGrid(1)-depth(iz)).AND.(iz.GT.1)) GOTO 11
        iz1=iz
        iz2=iz+1
        depfac=(zGrid(i)-(zGrid(1)-depth(iz1)))/((zGrid(1)-depth(iz2))-(zGrid(1)-depth(iz1)))
     ENDIF
!!!!* interpolate along depth- and time-coordinates:
     IF (iz1.EQ.nz_tempS) THEN
        IF (it1.EQ.nt_tempS) THEN
           tem(i)=temtim(nt_tempS,nz_tempS)
        ELSE
           tem(i)=timfac*(temtim(it2,nz_tempS)-temtim(it1,nz_tempS))+temtim(it1,nz_tempS)
        ENDIF
     ELSE
        IF (it1.EQ.nt_tempS) THEN
           tem(i)=depfac*(temtim(nt_tempS,iz2)-temtim(nt_tempS,iz1))+temtim(nt_tempS,iz1)
        ELSE
           val1=depfac*(temtim(it1,iz2)-temtim(it1,iz1))+temtim(it1,iz1)
           val2=depfac*(temtim(it2,iz2)-temtim(it2,iz1))+temtim(it2,iz1)
           tem(i)=timfac*(val2-val1)+val1
        ENDIF
     ENDIF
     !* current local impedance factor due to soil temperature:
     IF ((tem(i).GE.tempermax).OR.(tem(i).LE.tempermin)) THEN
        impt(i)=0.0_dp
     ELSE
        IF (topt.LT.tmid) THEN
           impt(i)=SIN(pi*((tem(i)-tempermin)/trange)**expo)
        ELSE
           impt(i)=SIN(pi*((tem(i)-tempermax)/(-trange))**expo)
        ENDIF
     ENDIF
  END DO
  RETURN
END SUBROUTINE Temper
!**********************************************************************
! current nodal solute concentration impedance factor values
SUBROUTINE ConTox
  USE GridData
  USE ConData
  USE SolData, ONLY: Conc
  IMPLICIT NONE
  REAL(sp) :: c0,c1,c2,c3
  INTEGER(ap) :: i
  !* current local impedance factor due to soil water solution concentration:
  c0=cmin
  c1=coptmi
  c2=coptma
  c3=cmax
  DO i=1,nPt
     impc(i)=0.0_dp
     IF (Conc(i).GT.c0.AND.Conc(i).LT.c1) impc(i)=(Conc(i)-c0)/(c1-c0)
     IF (Conc(i).GE.c1.AND.Conc(i).LE.c2) impc(i)=1._dp
     IF (Conc(i).GT.c2.AND.Conc(i).LT.c3)impc(i)=(Conc(i)-c3)/(c2-c3)
  END DO
  RETURN
END SUBROUTINE ConTox
!***********************************************************************************
! current local soil strength value:
SUBROUTINE StrLoc(corner,sLoc)
  USE Typedef
  USE StrData
  IMPLICIT NONE
  REAL(dp), INTENT(out):: sLoc
  INTEGER(ap), INTENT(in):: corner(8)
  sLoc=(s(corner(1))+s(corner(2))+s(corner(3))+s(corner(4))+s(corner(5))+&
       s(corner(6))+s(corner(7))+s(corner(8)))/8._dp
  RETURN
END SUBROUTINE StrLoc
!***********************************************************************************
! current local temperature value:
SUBROUTINE TemLoc(corner,tLoc)
  USE TempData
  USE RootData,ONLY: ltemp
  IMPLICIT NONE
  INTEGER(ap):: corner(8)
  REAL(sp):: tLoc
  IF (ltemp) THEN
     tLoc=(tem(corner(1))+tem(corner(2))+tem(corner(3))+tem(corner(4))+tem(corner(5))&
          +tem(corner(6))+tem(corner(7))+tem(corner(8)))/8._dp
  ELSE
     tLoc=topt
  ENDIF
  RETURN
END SUBROUTINE TemLoc
!***********************************************************************************
!Atmospheric conditions
SUBROUTINE Atmosphere(t)
USE TempData
IMPLICIT NONE
INTEGER(ap) :: aa
REAL(sp), INTENT(in) :: t

!interpolate atmospheric temperature, T_atm for time
IF (t.GE.time_TA(nt_tempS)) THEN
   Tatm_usr=T_atm(nt_tempS)
ELSE
   DO aa=1,nt_tempS
      IF (t.GE.time_TA(aa)) THEN
         Tatm_usr=(T_atm(aa+1)-T_atm(aa))/(time_TA(aa+1)-time_TA(aa))*t+T_atm(1)
      ENDIF
   ENDDO
ENDIF
Tatm_usr=Tatm_usr+273

!interpolate atmospheric pressure, P_atm for time
IF (t.GE.time_PA(nt_presA)) THEN
   Patm_usr=P_atm(nt_presA)
ELSE
   DO aa=1,nt_presA
      IF (t.GE.time_PA(aa)) THEN
         Patm_usr=(P_atm(aa+1)-P_atm(aa))/(time_PA(aa+1)-time_PA(aa))*t+P_atm(1)
      ENDIF
   ENDDO
ENDIF

!interpolate atmospheric pressure deficit, P_diff for time
IF (t.GE.time_PD(nt_presD)) THEN
   Pdiff_usr=P_diff(nt_presD)
ELSE
   DO aa=1,nt_presD
      IF (t.GE.time_PD(aa)) THEN
         Pdiff_usr=(P_diff(aa+1)-P_diff(aa))/(time_PD(aa+1)-time_PD(aa))*t+P_diff(1)
      ENDIF
   ENDDO
ENDIF

END SUBROUTINE Atmosphere
!***********************************************************************************

