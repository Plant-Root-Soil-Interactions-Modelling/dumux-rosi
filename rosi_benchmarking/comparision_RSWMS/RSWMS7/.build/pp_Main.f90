PROGRAM RSWMS_MAIN
  USE iso_c_binding
  USE typedef
  USE ParamData, ONLY: lvtk,lOutPartrace
  USE tmctrl
  USE SolData, ONLY: theta,conc,theta_old,ConO,Kode,con,hold,cap,hTemp,hnew
  USE CumData, ONLY: WatIn, SolIn,Q
  USE DoussanMat, ONLY: curr_BCtp,savelast,nBC_irecn,plantmatrix
  USE PlntData, ONLY: TotSur
  USE RootData, ONLY: lDou,lno_RWU
  USE sparsematrix
  USE RhizoData, ONLY: thetaNonEq, thetaNonEqOld,thetaTot, rhizoModel, lRhizo
  USE GridData, ONLY: nPt
  IMPLICIT NONE
  INCLUDE 'mpif.h' !if rswms runs on the cluster
  interface
     subroutine rearrange_r2p(a,b)
       Use Typedef
       real(dp), Intent(in)::a(:)
       real(dp),Intent(out)::b(:)
     end subroutine rearrange_r2p
  end interface
  INTEGER :: jj
  REAL(sp) ::t, dt
  LOGICAL :: tick=.TRUE.
  Integer :: error !if rswms runs on the cluster
  Call mpi_init(error) !if rswms runs on the cluster
  CALL RSWMS_INI(t,dt)
! Initiate Rhizo if rhizodynamics is considered
IF (lRhizo .AND. RhizoModel .EQ. 3) THEN
    DO jj=1,nPt
        thetaNonEqOld(jj) = thetaTot(jj)
    ENDDO
ENDIF
Timeloop: DO WHILE (.NOT.(ABS(t-tMax).LE.0.001_dp*dt))
   CALL CPU_TIME(time0)
   print*,'time for root growth, solute flow etc. -- start of new time step', time0-time3
   CALL timestep(t,dt)
   CALL RSWMS(t,dt)
   ! Updating thetaNonEqOld
   IF (lRhizo .AND. RhizoModel .EQ. 3) THEN
      DO jj=1,nPt
         thetaNonEqOld(jj) = thetaNonEq(jj)
      ENDDO
   ENDIF
      IF (savelast) THEN
         CALL OutFEM(t,0)
         IF(lvtk) CALL OutVTK(kOut)
         CALL FlxOut(kOut) !SWMS3D
         IF(lOutPartrace) CALL PartraceOut(t)
         IF (lDou) THEN
            CALL OutDou(t,0)
            IF(lvtk)CALL OutDouVTK(0)
         ENDIF
         OPEN (unit=9,File='SimulationTime.out',status='unknown')
         STOP
      ENDIF
      ! end of simulation?:
      IF (ABS(t-tMax).LE.0.001_dp*dt) THEN
         Print *, "end of sim"
         OPEN (unit=9,File='SimulationTime.out',status='unknown')
         DEALLOCATE(theta_old,hOld,hTemp,hNew,theta,conc)
         DeALLOCATE(Kode,Q)
         DEALLOCATE(conO,con,cap)
         DEALLOCATE(WatIn,SolIn)
         CALL Getout
      ENDIF
  END DO Timeloop
  Call MPI_FINALIZE(error) !if rswms runs on the cluster
END PROGRAM RSWMS_MAIN
!****************************************************************
