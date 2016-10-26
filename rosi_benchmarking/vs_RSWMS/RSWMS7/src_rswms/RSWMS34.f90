!Makefile .   To switch OpenMp on: pgf90 -mp -o executable file.f90
!                                   export OMP_NUM_THREADS= 2
!                                   ./executable
! To use openmp -> switch module OMP_LIB on!!! gfortran compiler can´t read OMP_LIB
! for gfortran only new compiler (v3) supports openmp
!******************************************************************************
!              Mathieu Javaux- Tom Schroeder- Valentin Couvreur-              *
!                     Natalie Schroeder- Katrin Huber-                        *
!                            Jan Vanderborght                                 *
!                         Fz-Juelich (Germany)                                *
!                                                                             *
!                                                                             *
!                                                                             *
!    R-SWMS_3D :An algorithm for three-dimensional, simultaneous modeling     *
!    of root growth, transient soil water flow, solute transport, and root    *
!    water and solute uptake.                                                 *
!    coupled with RootTyp (Pages et al.)                                      *
!                                                                             *
!                        Version 3.4.           August10                      *
!                                                                             *
!    Program Version for VAX - FORTRAN                                        *
!                                                                             *
!   based on the code of                                                      *
!   Francesca Somma & Volker Clausnitzer, University of California, Davis     *
!       (1998)                                                                *
!*****************************************************************************

!====================================================================
!> main program
!> ====================================================================

SUBROUTINE RSWMS(t,dt)
  
  USE typedef
  USE CumData
  USE PlntData
  USE tmctrl
  USE GridData
  USE RootData
  USE TempData
  USE MatData
  USE DoussanMat
  USE SparseMatrix
  USE soldata
  USE StrData, ONLY : ssmax
  USE ObsData
  USE DomData
  USE BoundData, ONLY : tQbcCh
  USE ParamData
  USE OMP_LIB
  USE SoluteRootMat
  USE ParticlesInRoot
  IMPLICIT NONE
  
  !variable declaration
  LOGICAL ReDo  
  REAL(sp) :: BCr,TR1
  !REAL(dp) ::t0,t1!delta2
  REAL(sp) :: dt,t
  !REAL(sp) :: ran
  INTEGER(ap) :: daytime(8),ipl,BCtp
  INTEGER(ap) :: imin,corner,igrow,irec,i
  CHARACTER form*3,file*8 
  
  
  !> --------------------------- Start of time loop ------------------------------
  
  solveroot_call=0
  !> time step adjustment:
  
222 CONTINUE
  !call outdouvtk(0,t)
  
  WRITE (*,'(/'' Time ='',1pE12.5)') t
  
  iter=0
  iter_root=0
  
  
  IF(.NOT.ALLOCATED(Q)) ALLOCATE(Q(nPt))

  !> Set soil boundary conditions
  CALL SetBC(t)					
    
  CALL CPU_TIME(time1)
  print*,'Time to set BC', time1-time0

  !> apply root extraction to FEM-sink term for each node:
  IF (.NOT.(switchSolve)) THEN
     WRITE (*,'(''+'',$)')
     IF (lDou) THEN
        CALL SolveRoot(t,dt,it1,iter_root)
        it1=.FALSE.
     ELSEIF ((.NOT.lno_RWU).AND.(.NOT.lCalloc)) THEN
        CALL SetBCroot(t,BCr,BCtp)
        CALL Setsnk(t,BCr,BCtp) !Also required when lDou because of solute uptake
     ENDIF
  END IF
   
  CALL CPU_TIME(time2)
  print*,'time to solve root', time2-time1

  !> solve soil water flow equation and calculate pressure heads:
  IF(lno_Archi) CALL IniMat
  CALL Water(t,dt,dtOpt,tOld,ReDo,IAD,IADN,IADD)

  CALL CPU_TIME(time3)
  print*,'time to solve water', time3-time2
  
  IF (ReDo) THEN
     savelast=.FALSE.
     GOTO 222
  ENDIF
  
  
  !> Calculate hormonal signaling
  IF(lSign .OR. lSign_inst) THEN
     CALL IniSolute(t)
     mhorm=0._dp
     msign_notrans=0._dp
     csign_notrans=0._dp
     TR1 = abs(PH_crit)-abs(delta_h)
     
     DO ipl=1, nplant
        DO irec=1,nrec
           IF(abs(Phr(irec+1,ipl)).GT. TR1) THEN
              mhorm(irec,ipl)=sign_in*(abs(Phr(irec+1,ipl))-TR1)*dt*segmas(irec)
           END IF
        END DO
     END DO
     msign_notrans=sum(mhorm) 
     IF (msign_notrans .GT. 0._dp) THEN
        mcol_i = mcol_i + msign_notrans
        csign_notrans = mcol_i/vol_buff
        mcol_i = mcol_i - csign_notrans*Tact(1)*dt
     END IF
     IF(lSign) CALL SoluteRoot(icount,dt,t)
  END IF

  !> Calculate solute uptake using partrace and SoluteRoot
  IF(lPartUp) THEN
     CALL IniSolute(t)
     mupt_diff = 0._dp
     mupt_adv = 0._dp
     CALL SoluteRoot(icount,dt,t)
  END IF

  
  !> solve solute transport equation and calculate concentrations:
  IF (lChem) THEN
     CALL Solute(dt,t,tPulse,dtMaxC,KodCB,IAD,IADN,IADD,icount)
  ENDIF
  !> calculate water and solute (if included) mass balance
  iCount=1
  CALL SubReg(t,iCount)
  !> calculate actual transpiration rate Tact (global variable in MODULE PlntData):
  IF (lCalloc.AND.(.NOT.(lDou))) CALL ActTrs
  IF (lCalloc) THEN
     DO ipl=1,nplant 
        IF (ipl.GE.2) STOP 'Assimilate allocation currently doesnt work with multiple plants'
        !> translate transpired water into biomass:
        CALL Effncy(t,rs,concrs,W)
        delBM=W*Tact(ipl)*dt
        
        !> partition that goes to new roots is accumulated in dtroot
        !> until next root growth step:
        CALL Ratio(t,rs,concrs,RSR) 
        dr=RSR/(1.+RSR)*delBM
        dmroot=dmroot+dr 
        
        !> remainder of delBM is partition that goes to shoot:
        dsh=delBM-dr
        mshoot=mshoot+dsh
        
        !> calculate current leaf area:
        CALL Leaves(t,LAmsh)
        LA=LA+dsh*LAmsh
     END DO
  ENDIF
  
  !> save Doussan data before deleting of vector and root data
  IF (ABS(t-tOut(kOut)).LE.0.001_sp*dt) THEN 
     !> Output for Doussan
     IF(lDou) THEN
        CALL OutDou(t,kout)
        WRITE(*,*)'outRoo. correctly written'
     END IF
     IF(.not.lno_root_growth) CALL OutRoo(t)
     IF(continu)  CALL OutTrans(t,kout)
     IF(lvtk .AND. (lDou .OR. .NOT.lno_root_growth)) THEN
        CALL OutDouVTK(kout,t)
        if(lSign .or. lPartUp)  CALL OutParticleVTK(kout,t)
     END IF
     IF(lPartUp .OR. lSalinity) CALL OutVTKCouple(kOut,concR)
     !kout = kout+1          
  ENDIF
  
  !> general output:
  if(.not.lno_RWU) CALL WriteLog(t)
  
  !> if appropriate, grow root and update the distribution functions 'betaw'
  !> and 'betac' (extraction intensity per volume, [L^-3])
  IF (ABS(t-tCallR).LT.0.001_dp*dt .AND. (.NOT.lno_root_growth)) THEN
     
     !update soil strength values:
     CALL Solstr(hNew)
     
     !> if temperature input provided, update temperature values:
     IF (ltemp) CALL Temper(t)
     
     !> if nutrient deficiency/ion toxicity data are provided, update
     !> the corresponding impedance nodal values:
     IF (ltoxi) CALL ConTox

     !> let root system grow:
     CALL Root(t,dt)
     tCallR=tCallR+dtRoot

     !> update Dousssan weighing factors and matrices
     IF (lDou) THEN
        !> delete the old vector of list
        IF (nrecold.NE.0) THEN
           DO ipl=1, nplant
              CALL SM_delete(plantmatrix(ipl))
           END DO
        ENDIF
        !            deallocate(listvec)
        ! reshape transroot for doussan
        IF(continu) transroot(nrec+1:nrec+ngrow,:,:,1)=transtip(1:ngrow,:,:,1)
        DEALLOCATE(nBC_irecn)
        !            if ((ave) .or. (eqDis)) deallocate(voxel_node)
        CALL SetupDou(t,dt)
        !it1=.TRUE.
        !ELSE
        !< get betaw and betac distribution:
        !   CALL BetDis(t)
        ! normalize above so that integral of betaw and betac over 3D-domain is equal to unity:
        !   CALL BetNrm
     ENDIF
  ENDIF

  !> in case of a static root system with ages smaller than runtime and variable conductivities over time
  IF(lUpdate_growth .AND. lDou) THEN
     !> delete the old vector of list
     IF (nrecold.NE.0) THEN
        DO ipl=1, nplant
           CALL SM_delete(plantmatrix(ipl))
        END DO
     ENDIF
     DEALLOCATE(nBC_irecn)
     CALL SetupDou(t,dt)
  END IF
  
  !> ouputs for observation probes
  IF (ObsOK) THEN
     IF (dtprobe.EQ.999) THEN
        CALL OutObsProbe(t)
     ELSEIF (ABS(t-tOuProbe(kOuProbe)).LE.0.001_dp*dt) THEN
        CALL OutObsProbe(t)
        kouProbe=kouProbe+1
     ENDIF
  ENDIF
  
  !> z-profiles
  IF (profOK) THEN
     IF (dtprof.EQ.999) THEN
        CALL Zprofiles(t)
     ELSEIF (ABS(t-tOuProf(kOuProf)).LE.0.001_dp*dt) THEN
        CALL Zprofiles(t)
        kouProf=kouProf+1
     ENDIF
  ENDIF
  
  !> FEM, root output at specified points in time: (independent of DoussanMat parameters
  IF (ABS(t-tOut(kOut)).LE.0.001_sp*dt) THEN
     !< in case macroscopic hydraulic properties need to be saved at each outfem time
     IF (lCou.AND.(lGap.OR.lAQPc)) CALL SSFdis(kOut) 
     CALL OutFEM(t,kOut)
     IF(lvtk .AND. .NOT. (lPartUp .OR. lSalinity)) CALL OutVTK(kOut)
     CALL FlxOut(kOut) 
     IF(lOutPartrace) CALL PartraceOut(t)
     kOut=kOut+1
  ENDIF
  
  !> check next BC time for root
  IF ((.NOT.lno_RWU).AND.(.NOT.lCalloc)) THEN
     IF (ABS(t-tBCr(kBCr)).LE.0.001_dp*dt) THEN
        kBCr=kBCr+1
     ENDIF
  ENDIF

  !> pressure heads for new time level:
  hOld =hNew
  hTemp=hNew

  !> ------------------------ End of time loop ----------------------------------
END SUBROUTINE rswms

!**************************************************************************
!> ### Initializing RSWMS ###
subroutine RSWMS_INI(t,dt)
  USE iso_fortran_env
  USE typedef
  USE CumData
  USE PlntData
  USE tmctrl
  USE GridData
  USE RootData
  USE TempData
  USE MatData
  USE DoussanMat
  USE SparseMatrix
  USE soldata
  USE RhizoData
  USE ObsData
  USE DomData,  ONLY :zmax
  USE BoundData, ONLY : tQbcCh
  USE ParamData
  USE disToRoot, ONLY: RStat
  !USE OMP_LIB
  
  IMPLICIT NONE
  REAL(sp), INTENT (out):: t,dt
  INTEGER(ap) :: TimeArray(3),i,ipl
  INTEGER(ap) :: daytime(8),irec,igrow
  INTEGER(ap) :: j,nn,un,istat,pid
  INTEGER(ap) :: ttt
  INTEGER, ALLOCATABLE,DIMENSION (:) :: seed
  REAL(sp) :: random
  CHARACTER form*3,file*8
  
  call RANDOM_SEED(size = nn)
  ALLOCATE(seed(nn))
		
  !> initiate random number generator:
  CALL itime(TimeArray)
  random = ran(TimeArray(1)+TimeArray(2)+TimeArray(3))
  CALL DATE_AND_TIME(values=daytime)
  ttt = INT(daytime(5)*3600+daytime(6)*60+daytime(7))  
  pid = getpid()
  ttt = ieor(ttt, int(pid, kind(ttt)))
  !seed = 0 !INT(daytime(5)*3600+daytime(6)*60+daytime(7))  
  !PRINT *,'need random number generator ',seed, ran(seed)
  
  do i = 1,nn
     seed(i) = lcg(ttt)
  end do
  Print*, 'seed', seed
  CALL RANDOM_SEED(put=seed)
  DEALLOCATE(seed)
  CALL RANDOM_NUMBER(rand)
	

  !> open simulation summary file
  OPEN (UNIT=15,FILE='out/simul_sum.out',STATUS='UNKNOWN')
  WRITE(15,112)
112 FORMAT('+++++ SIMULATION INPUT SUMMARY +++++ ')
  CLOSE(15)
  ObsOK=.FALSE.
  
!> -----------------------------------INPUT---------------------------------------
  !> get input for the specific problem (domain, BC, IC, control parameters):
  CALL Applic(dt)
  !allocate soil variables (module soldata)
  ALLOCATE(Vx(nPt),Vy(nPt),Vz(nPt))
  ALLOCATE(theta_old(nPt),theta(nPt))
  ALLOCATE(conO(nPt),con(nPt),cap(nPt))
  ALLOCATE(l_elmMacro(1:nElm))
  !> allocate thetaTot's if rhizosphere model. 
  IF (lRhizo) ALLOCATE(thetaTot(nPt), thetaNonEq(nPt), thetaNonEqOld(nPt), tauTht(nPt), hEqRhizo(nPt)) 
  l_elmMacro = .false.
  !> calculate Wn (new Javaux)
  CALL CalcWnodes
  IF (lChem)  CALL ChemIn  !> get solute transport information
  CALL SoilIn !> get soil material input and set up K,C-table for interpolation
  
  
  !> list for voxels of mixed material 
  IF(lDou .AND. lMacro) THEN
     ALLOCATE(MacroList(1:nElm,1:6))
     MacroList = 0
     ALLOCATE(n_neigh(1:nElm))
     n_neigh=0
     CALL ListMacro
  END IF

  IF (.NOT.(lno_Archi)) THEN

     !> initial root system and growth parameters:
     CALL RootIn(t)                  
     !> get Doussan model input information
     IF (lDou) THEN
        CALL DouIn
     ELSEIF (lFed) THEN
        CALL FedIn
        CALL SetupDou(t,dt)	!> Calculates the root system position as compared to the position of the plant collar and of the limits of the periodic domain if so.
        CALL RLDdis(t,kout)		!> Calculates the root length & surface density distribution
        !         CALL BetNrm
     ELSEIF (lCou) THEN
        WRITE(*,'(/''SSF, Krs and Kcomp generated from RootSys and CondRoot.in'')')
        CALL CouIn
        CALL DouIn
        CALL SetupDou(t,dt)
        CALL SSFdis(kout)
        IF (ldJvL) CALL RLDdis(t,kout)
     ELSEIF(.NOT. lDou .AND. lSomma_growth) THEN
        ALLOCATE(transroot(0:maxrec+maxgrw,1:2,1:isubmax,1:nplant))
        transroot=0
        ALLOCATE(transtip(0:maxgrw,1:2,1:isubmax,1:nplant))
        transtip=0 
        ALLOCATE(nsub(0:maxgrw+maxrec,1:nplant))
        nsub = 1
     END IF
  ELSEIF (lCou) THEN
     CALL CouIn			!Initial time determined by PlntIn
  ELSEIF (lFed) THEN
     CALL FedIn			!Initial time determined by PlntIn
  ENDIF
  
!> soil temperature over depth / in time and atmosperic conditions:
  IF (ltemp .OR. lCalloc)  CALL TempIn 
  
  IF ((.NOT.lno_RWU).AND.(.NOT.lretry)) THEN	!> log file created even for RWU models with no root architecture
     !> open log file:
     WRITE (file,'(A7)')'out/log'
     DO ipl=1,nplant
        WRITE (file(8:8),'(I1)') ipl
        OPEN (UNIT=10,FILE=file,STATUS='UNKNOWN')
        IF (lDou) THEN
           IF(lCalloc) THEN
              WRITE (10,90)
           ELSEIF (lSign .OR. lSign_inst) THEN
              WRITE (10,91)
           ELSEIF (lPartUp) THEN
              WRITE (10,94)
           ELSE
              WRITE (10,93)
           ENDIF
        ELSEIF (lCou) THEN
           WRITE (10,92)
        ENDIF
90      FORMAT('  Time       Tpot       Tact    grwfac     sAvg       mShoot     mRoot     LA     H_collar     sign_conc    TpLA    TpLA_pot     #stressed_Nodes ')
91      FORMAT('  Time       Tpot       Tact   sAvg     soil_lim_nodes    H_collar   mcol   signal_conc        msign_notrans    csign_notrans    res_time')
92      FORMAT('  Time       Tpot       Tact  soil_lim_nodes    H_collar      H_seq')
93      FORMAT('  Time       Tpot       Tact   sAvg    soil_lim_nodes    H_collar')
94      FORMAT('  Time       Tpot       Tact   sAvg     soil_lim_nodes    H_collar   mcol   segconc     #particle      m_upt')
        CLOSE (10)
     ENDDO
  ELSE
     t=tQbcCh(1)!initial time =first time of the soil BC.
  ENDIF
  
  !> get plant parameter input:
  CALL PlntIn(t)	!> Reads BCroot.in which is necessary even for RWU models without root architecture, and reads Plant.in which is necessary for assimilate allocation

  IF (.NOT.lretry) THEN
     !> open balance file:
     OPEN (UNIT=10,FILE='out/balance.out',STATUS='UNKNOWN')
     WRITE(10,110)
110  FORMAT(' Absolute and relative mass balance error for water and solute transport ')
     CLOSE (10)
     !> open removal file
     OPEN (UNIT=10,FILE='out/remove.out',STATUS='UNKNOWN')
     WRITE(10,111)
111  FORMAT(' Total amount of water and solute removed by the root system, ',/,&
          ' by zero and first order reactions, and by drainage at the bottom of the domain. ')
     CLOSE(10)
     !> open Observation node file
     IF (ObsOK) THEN
        CALL ObsIni
     ENDIF
     IF (profOK) THEN
        OPEN (UNIT=121,FILE='out/ProfileTH.out',STATUS='UNKNOWN')
        OPEN (UNIT=122,FILE='out/ProfilePH.out',STATUS='UNKNOWN')
        OPEN (UNIT=123,FILE='out/ProfileS.out',STATUS='UNKNOWN')
        WRITE(form,'(I3)')(nz)
        WRITE (121,'(/''Averaged water content profile for each time step.'')')
        WRITE (121,'(/''Time   Z- water content'')')
        WRITE (121,'((12X)'//form//'(1X,f12.4))') (zmax-dzgrid*(i-1),i=1,nz)
        WRITE (122,'(/''Averaged water potential profile for each time step.'')')
        WRITE (122,'(/''Time   Z- water potential'')')
        WRITE (122,'((12X)'//form//'(1X,f12.4))') (zmax-dzgrid*(i-1),i=1,nz)
        WRITE (123,'(/''Total Sink profile for each time step.'')')
        WRITE (123,'(/''Time   Z- sink'')')
        WRITE (123,'((12X)'//form//'(1X,f12.4))') (zmax-dzgrid*(i-1),i=1,nz)
     ENDIF
  ENDIF
  ALLOCATE(IADN(nPt),IADD(nPt))
  ALLOCATE(IAD(maxbnd,nPt))
  CALL IADMake(nPt,nElm,maxbnd,IAD,IADN,IADD)
  !--------------------------------------------------------------------------------
  !> first time step:
  dtOpt=dt
  if(.not.lno_root_growth) then
     tCallR=t+dtRoot
  else
     tCallR = 1e+30_sp
  end if
  if(lno_root_growth) tCallR = 1e+30
  dmroot=0.0_dp
  tstart = t
  
  !> find next output time:
1 kOut=kOut+1
  IF (kOut.LE.nOut) THEN
     IF (tOut(kOut).LE.t) GOTO 1
  ENDIF
2 kaxemg=kaxemg+1
  IF (kaxemg.LE.naxemg) THEN
     IF (tnewax(kaxemg).LT.t) GOTO 2
  ENDIF
  !> Z-profiles
  touProf=touProf+t !> to start at initial time < rootin
  IF ((ProfOK).AND.(dtProf.NE.999)) THEN
3    kouProf=kouProf+1
     IF (kouProf.LE.nouProf) THEN
        IF (touProf(kouProf).LE.t) GOTO 3
     ENDIF
  ENDIF
  !> probes
  touProbe=touProbe+t !> to start at initial time < rootin
  IF ((ObsOK).AND.(dtProbe.NE.999)) THEN
4    kouProbe=kouProbe+1
     IF (kouProbe.LE.nouProbe) THEN
        IF (touProbe(kouProbe).LE.t) GOTO 4
     ENDIF
  ENDIF
  t_begin=t
  
  !> if Doussan, calculate weighing function and matrices
  IF (lDou) THEN
     CALL SetupDou(t,dt)
     !> initialize stress
     stressBC=.FALSE.
  ELSEIF (lFed) THEN
     !> get betaw and betac initial distribution:
     CALL BetDis(t)
     !> ...and normalize:
     CALL BetNrm
  ENDIF
  IF ((.NOT.lno_RWU).AND.(.NOT.lCalloc)) THEN
     !> find next time change in root BC
5    kBCr=kBCr+1
     IF (kBCr.LE.nBCr) THEN
        IF (tBCr(kBCr).LE.t) GOTO 5
     ENDIF
  ENDIF
  
  !> calculate parameters:
  IF (lRhizo)  CALL RStat(t)
  CALL SetMat(0)
  
  CALL CalcGeom
  CALL CalcWidthSWMS_new
  iCount=0
  CALL SubReg(t,iCount)
  !> initialize it1 (first iteration after a root growth)
  it1=.TRUE.
  itMaxRoot=itMax
  IF (oldT) THEN
     switchsolve=.FALSE.
  ELSE
     switchsolve=.TRUE.
  ENDIF
  
  IF(lSign) THEN 
     ALLOCATE (l_SignOn(1:nrec+1))
     l_SignOn = .FALSE.
  END IF
  contains
    function lcg(s)
		INTEGER :: lcg
        INTEGER :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967)
              end if
              s = mod(s * 279470, 4294967)
              lcg = int(mod(s, int(huge(0), ap)), kind(0))
    end function lcg
  
end subroutine RSWMS_INI
!******************************************************************************
!> ### adjusts next time step ###
subroutine timestep(t,dt)
  
  USE tmctrl
  USE ParamData, ONLY: iter,mxbcch
  USE PlntData, ONLY: nBCr,tbcr
  USE RootData, ONLY: lCalloc,lno_RWU
  USE ObsData, ONLY: ObsOK, ProfOK
  IMPLICIT NONE

  REAL(sp) :: t, dt
  REAL(sp) :: tFemRoo

  tOld=t
  dtOld=dt
  IF (kOut.LE.nOut) THEN
     tFEMRoo=tOut(kOut)
  ELSE
     tFEMRoo=tMax
  ENDIF
  
  IF ((.NOT.lno_RWU).AND.(.NOT.lCalloc)) THEN
     IF (kBCr.LE.nBCr) THEN
        tcBCr=tBCR(kBCr)
     ELSE
        tcBCr=tMax
     ENDIF
  ELSE
     tcBCr=tMax
  ENDIF
  IF ((ProfOK).AND.(dtProf.NE.999).AND.(kouProf.LE.nouProf)) THEN
     tProf=tOuProf(kouProf)
  ELSE
     tProf=tMax
  ENDIF
  IF ((ObsOK).AND.(dtProbe.NE.999).AND.(kouProbe.LE.nouProbe)) THEN
     tProbe=tOuProbe(kouProbe)
  ELSE
     tProbe=tMax
  ENDIF
  CALL TmCont(iter,t,dt,dtOpt,tCallR,tFEMRoo,dtMaxC,tcBCr,tProf,tProbe)
  
  t=t+dt

END subroutine timestep
!**************************************************************************
