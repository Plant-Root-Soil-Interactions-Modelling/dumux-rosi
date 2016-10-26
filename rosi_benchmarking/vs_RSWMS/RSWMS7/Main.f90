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

#ifdef WITH_PARTRACE
  USE GridData, ONLY: continu,nex,ney,nez,nx,ny,nz,nPt,nElm,sink_cube,betac_cube2
  USE SolData, ONLY: Vx,Vy,Vz,concPar,concR, massPar, massR
  USE TempData, ONLY: tstart 
  USE ParamData, ONlY: lPartUp,lSalinity
  USE SoluteRootMat
#endif  

  IMPLICIT NONE
  INCLUDE 'mpif.h'          !if rswms runs on the cluster
  
  interface
     
     subroutine rearrange_r2p(a,b)
       Use Typedef
       real(dp), Intent(in)::a(:)
       real(dp),Intent(out)::b(:)
     end subroutine rearrange_r2p
  end interface
  
#ifdef WITH_PARTRACE
  INTEGER(ap):: i, iNodes
  REAL(dp) :: timeP,dummy_soil,dummy_root
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: theta_P, v1, v2, v3, sinkE, sinkE_temp, fact_temp, &
       sink_in
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: velocity_P
#endif  
  INTEGER :: jj
  REAL(sp) ::t, dt
  LOGICAL :: tick=.TRUE.

  Integer :: error       !if rswms runs on the cluster
  Call mpi_init(error)   !if rswms runs on the cluster
  CALL RSWMS_INI(t,dt)  

#ifdef WITH_PARTRACE
  CALL initPartrace()
  if(lPartUp) CALL IniSolute(t)

  IF (continu) THEN 
     iNodes = (nx+1)* (ny+1)* nz
  ELSE
     iNodes = nPt
  ENDIF
  
   IF (.NOT.(ALLOCATED(theta_P))) THEN
     ALLOCATE (theta_P(iNodes))
     ALLOCATE (velocity_P(iNodes,3))
     ALLOCATE (v1(iNodes))
     ALLOCATE (v2(iNodes))
     ALLOCATE (v3(iNodes))
     ALLOCATE (concPar(nElm))
     ALLOCATE (concR(nElm))
     ALLOCATE (sinkE(nElm))
     ALLOCATE (sinkE_temp(nElm))
     ALLOCATE (massR(nElm))
     ALLOCATE (massPar(nElm))
     if(lPartUp)  ALLOCATE (fact_temp(nElm))
  END IF

#endif

! Initiate Rhizo if rhizodynamics is considered
IF (lRhizo .AND. RhizoModel .EQ. 3) THEN
    DO jj=1,nPt
        thetaNonEqOld(jj) = thetaTot(jj)
    ENDDO
ENDIF

Timeloop:  DO WHILE (.NOT.(ABS(t-tMax).LE.0.001_dp*dt))

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
   
   
#ifdef WITH_PARTRACE

     CALL REARRANGE_RSWMS_TO_PARTRACE(velocity_P,theta_P,iNodes)
     
     v1 = velocity_P(:,1)
     v2 = velocity_P(:,2)
     v3 = velocity_P(:,3)

     if (.not. lPartUp) then
        sinkE = 0._dp !define sink for passive, active or 0 for exclusion
        CALL Rearrange_r2p(sinkE,sinkE_temp)
        CALL setpartracefields(inodes,v1,v2,v3,theta_P,sinkE_temp,nElm)
     else
        CALL Rearrange_r2p(fact, fact_temp)
        CALL setpartracefields(inodes,v1,v2,v3,theta_P,fact_temp,nElm)
     end if

     !print*, 'fact', fact

     dummy_soil =  sum(mass_elm)  ! solute mass in soil before calling partrace
     dummy_root = sum(m_upt)  ! solute mass taken up by roots
   
     timeP = t-tstart
     CALL runpartrace(timeP)
     CALL setConcentration(nElm,concPar)
     CALL Rearrange_p2r(concPar,concR)
    
     CALL setMass(nElm, MassPar)
     CALL Rearrange_p2r(massPar, massR)


     !print*, 'uptake from partrace', dummy_soil - sum(massR), dummy_soil, sum(massR)
     !print*, 'uptake by roots', dummy_root

     if(lSalinity) CALL Salinity

!!$     print*,'hier tout(kout)',t, tout(kout), 0.001_sp*dt, t-tOut(kOut)

#endif


      IF (savelast) THEN
         CALL OutFEM(t,0)
#ifdef WITH_PARTRACE
         IF(lvtk .AND. (lPartUp .OR. lSalinity)) CALL OutVTKCouple(kOut)
#else
         IF(lvtk) CALL OutVTK(kOut)
#endif
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
#ifdef WITH_PARTRACE
         DEALLOCATE(Vx,Vy,Vz)
         CALL closePartrace()
#endif
         CALL Getout

      ENDIF

        


  END DO Timeloop
  Call MPI_FINALIZE(error)    !if rswms runs on the cluster
 
     
END PROGRAM RSWMS_MAIN
!****************************************************************

