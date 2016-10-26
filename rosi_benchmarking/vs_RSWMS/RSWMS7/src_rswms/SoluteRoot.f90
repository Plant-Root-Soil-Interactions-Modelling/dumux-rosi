! ==============================================================================
! Source file SOLUTE ROOT ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
!>  \brief module for particle tracking inside the roots
!> Solute root model implementation through 9 Subroutines
!>	- MoveParticles
!>      - MoveSingleParticle
!>	- InsertParticle
!>      - RemoveParticle
!>	- ParticleToRoot
!>      - SoluteRoot (main)
!>      - OutParticleVTK (Output.f90)
!>      - UpdateRootSys
!>      - ReduceParticleMass

MODULE ParticlesInRoot

CONTAINS

    !> move a single particle by a given velocity
    !> calculates the new position, new segment number and
    !> removes particle form list if it leaves the root (collar)
  SUBROUTINE MoveSingleParticle(pFirst,P,ipl,dt,numVic,t)
    USE TypeDef
    USE DoussanMat, ONLY: veloRoot,nrec
    USE RootData, ONLY: irecpr,segSoluteMass,seglen
    USE SoluteRootMat, ONLY: Particle, totalParticleNum,irecfollow,numfollow,retard
    USE ParamData, ONLY: lPartUp 
    IMPLICIT NONE

    TYPE(Particle) ,POINTER ,INTENT(inout) :: P,pFirst
    TYPE(Particle) ,POINTER :: Pprev
    REAL(sp) ,INTENT(in) ::dt,t
    INTEGER(ap) ::ipl,numVic
    REAL(dp):: k, moveVec,tm,lm1,lmt
    LOGICAL :: lforw 
  	!> \param pFirst pointer to first particle in particlelist
	!> \param P pointer to the current particle
	!> \param ipl current root system (not used)
	!> \param dt time step size
	!> \param numVic number of particle (victims) which will be removed in the next time step
	!> \param t current simulation time


    tm  = 0. ! time a particle already was moved 
    lm1 = 0. ! moved length in segment
    lmt = 0. ! moved length within one time step

    DO WHILE((tm.LT.dt ))
       moveVec    =  veloRoot(P%segNum,ipl) *retard(P%segNum) * (dt-tm)
       P%position = P%positionOld + moveVec

       IF (moveVec .GE. 0) THEN ! forward movement
          lforw = .TRUE.    
          lm1 = seglen(P%segNum) - P%positionOld
          lmt = lmt + lm1
          tm  = tm + lm1/(veloRoot(P%segNum,ipl)*retard(P%segNum))  
          
          IF (P%position .GT. seglen(P%segNum)) THEN ! Change segment
             P%segNum = irecpr(P%segNum)
             
             IF (P%segNum.EQ.0)  THEN 
                CALL RemoveParticle(pFirst,P,Pprev,lforw,numVic,t)
                totalParticleNum =  totalParticleNum -1
                P => Pprev
                EXIT    
             END IF

             P%segLen = seglen(P%segNum)
             P%positionOld = 0. 
             

          END IF

 
       ELSE ! backward movement
          lforw = .FALSE.
          lm1   = P%positionOld
          lmt   = lmt + lm1
          tm    = tm + lm1/abs(veloRoot(P%segNum,ipl)*retard(P%segNum))
          
          IF(P%position .LT. 0) THEN !Change segment

             IF (numfollow(P%segNum).GT.0) THEN
                P%segNum = irecfollow(P%segNum,1) !CHANGE HERE FOR BIFURCATION!
                P%positionOld = seglen(P%segNum)
             ELSE !root tip reached
                IF(lPartUp) THEN   
                   !Particles get stuck in end of segment
                   p%position    = 0.0
                   P%positionOld = 0.0
                   P%segNum      = P%segNum
                   p%segLen      = segLen(P%segNum)
                ELSE    
                   ! hormonal signalling                
                   CALL RemoveParticle(pFirst,P,Pprev,lforw,numVic,t) 
                   totalParticleNum =  totalParticleNum -1
                   P => Pprev
!!$                   !if (associated(P)) then
!!$                   !   P%segLen = seglen(irecfollow(P%segNum,1))
!!$                   !   P%positionOld = P%segLen
!!$                   !end if
                   EXIT
                END IF
             END IF
          END IF

       END IF
 
    END DO

 
!!$    k = 0. 
!!$
!!$
!!$    DO WHILE((k.LT.1 ))
!!$       !print*,P%segNum
!!$       moveVec =  veloRoot(P%segNum,ipl) *retard(P%segNum) * dt
!!$       P%position = P%positionOld + moveVec /(1.-k)
!!$      
!!$       if (moveVec.ge.0) then
!!$          k = (seglen(P%segNum) -  P%positionOld)/moveVec*(1.-k)
!!$          lforw=.TRUE.
!!$       else
!!$          !backward
!!$          k = -P%positionOld/moveVec*(1.-k)
!!$          lforw=.FALSE.
!!$          IF (numfollow(P%segNum).GT.0) THEN
!!$             P%segNum = irecfollow(P%segNum,1)
!!$          ELSE
!!$             P%segNum =  P%segNum
!!$
!!$          END IF
!!$          !print*,'segnum - out',p%segNum
!!$       end if
!!$
!!$
!!$       IF (P%position.GT. seglen(P%segNum)) THEN !Forward
!!$            P%segNum = irecpr(P%segNum)
!!$
!!$          IF (P%segNum.EQ.0)  THEN  
!!$            !print*,'CALL ME FORWARD!', P%segNum, 'segNum', irecpr(P%segNum)
!!$             CALL RemoveParticle(pFirst,P,Pprev,lforw,numVic,t)
!!$             totalParticleNum =  totalParticleNum -1
!!$             P => Pprev
!!$             !print*,'segnum',P%segnum
!!$             EXIT    
!!$          END IF  
!!$          P%segLen = seglen(P%segNum)
!!$          P%positionOld = 0. 
!!$
!!$       ELSEIF (P%position.LT. 0) THEN !Backward
!!$
!!$          IF (irecfollow(P%segNum,1).EQ.0)  THEN 
!!$             ! print*,'CALL ME BAckWARd!'
!!$             IF(lPartUp) THEN   
!!$                !Particles get stuck in end of segment
!!$                p%position=0.0
!!$                P%positionOld=0.0
!!$                P%segNum=P%segNum
!!$                p%segLen = segLen(P%segNum)
!!$             ELSE
!!$
!!$                CALL RemoveParticle(pFirst,P,Pprev,lforw,numVic,t)  !for hormonal signalling
!!$                totalParticleNum =  totalParticleNum -1
!!$                !print*,'totalparticlenumMOVE',totalparticlenum
!!$                P => Pprev
!!$                !print*,'jetzt hier'
!!$                if (associated(p)) then
!!$                   P%segLen = seglen(irecfollow(P%segNum,1))
!!$                   P%positionOld = P%segLen
!!$                end if
!!$                EXIT
!!$             END IF
!!$          END IF
!!$       END IF
!!$
!!$
!!$    END DO

END SUBROUTINE MoveSingleParticle
!***************************************************************************
  !> go over list of particles and move each particle
  !> with a given velocity
SUBROUTINE MoveParticles(pFirst,dt,icount,numVic,t)
  USE TypeDef
  USE ParamData, ONLY:  maxParticle, pi
  USE DoussanMat, ONLY: veloRoot,nplant
  USE RootData, ONLY: seglen,segsur,segconc,segSoluteMass,nrec,lno_root_growth
  USE SoluteRootMat
  IMPLICIT NONE

  REAL(sp) ,INTENT(in) :: dt,t
  TYPE(Particle) ,POINTER ,INTENT(inout) :: pFirst
  INTEGER(ap) :: P,irec,ipl,numVic
  INTEGER(ap), INTENT(in) :: icount

  	!> \param pFirst pointer to first particle in particlelist
	!> \param icount icount=0 at first time step
	!> \param dt time step size
	!> \param numVic number of particle (victims) which will be removed in the next time step
	!> \param t current simulation time
  
  !call MakeFollowList(irecfollow)
    if(icount.EQ.0)  then
       call MakeFollowList
    elseif(.not.lno_root_growth) then
       call MakeFollowList
    endif
       

  DO ipl=1,nplant
     segSoluteMass = 0.
     pParticle => firstP !point to beginning of the list
     DO WHILE (ASSOCIATED(pParticle)) !end of list reached

        CALL MoveSingleParticle(pFirst,pParticle,ipl,dt,numVic,t)!move one particle
        IF (ASSOCIATED(pParticle)) THEN
           pParticle%positionOld = pParticle%position
           pParticle => pParticle%next !go to next link in the list
        ELSE
           pParticle => pFirst
        END IF
     END DO
  END DO

    !print*,'totalparticlenumout',totalparticlenum
END SUBROUTINE MoveParticles
!***************************************************************************
  !> this routine creates a new particle (new entry in linked list) 
  !> which is added at the beginning (left hand side) of the list
SUBROUTINE InsertParticle(pFirst,irec,pos,mass,t)
  USE TypeDef
  USE DoussanMat, ONLY: nplant
  USE RootData, ONLY: seglen,irecsg
  USE SoluteRootMat, ONLY: Particle
  IMPLICIT NONE

  INTEGER(ap),INTENT(in):: irec
  REAL(dp),INTENT(in):: mass
  REAL(dp),INTENT(in):: pos
  REAL(sp),INTENT(in):: t
  INTEGER(ap),SAVE :: particleID=0
  TYPE(Particle) ,POINTER,INTENT(inout) :: pFirst
  TYPE(Particle) ,POINTER :: pParticle
  	!> \param pFirst pointer to first particle in particlelist
	!> \param irec current root segment
	!> \param pos location of particle within root segment
	!> \param mass mass of current particle
	!> \param t current simulation time
  

  particleID =  particleID +1
  ALLOCATE(pParticle)      ! new particle

  NULLIFY(pParticle%prev)  ! prev of new particle points to Null
  pParticle%next => pFirst ! next of new particle is pFirst

  ! prev of pFirst points to new particle
  IF(ASSOCIATED(pFirst)) pFirst%prev => pParticle
  pFirst => pParticle      ! pFirst is now new particle

  pParticle%ID = particleID   ! ID
  pParticle%segNum = irec     ! segNum
  pParticle%positionOld = pos ! old position (same as new)
  pParticle%position = pos    ! new position 
  pParticle%mass = mass      ! particle mass
  pParticle%segLen = seglen(irec)  
  pParticle%partOrig = t      !origination time of particle


END SUBROUTINE InsertParticle
!***************************************************************************
  !> this routine removes a particle (entry in linked list) and 
  !> links the pointers of the prev and next particle (entry in list) 
  !> adds up particle masses that are removed at the collar
  !> calculates average residence time of particles
SUBROUTINE RemoveParticle(pFirst,pVictim,pVictim_prev,lforw,numVic,t)
  USE TypeDef
  USE SoluteRootMat, ONLY: Particle
  USE RootData, ONLY: m_in,res_t
  IMPLICIT NONE

  TYPE(Particle) ,POINTER ,INTENT(inout) :: pFirst
  TYPE(Particle) ,POINTER ,INTENT(inout) :: pVictim
  TYPE(Particle) ,POINTER ,INTENT(out) :: pVictim_prev
  REAL(sp),INTENT(in) :: t
  INTEGER(ap) :: numVic
  LOGICAL :: lforw
  !> \param pFirst pointer to first particle in particlelist
  !> \param pVictim particle to be removed from the system
  !> \param pVictim_prev previous particle from pVictim
  !> \param lforw true if upward transport/flow
  !> \param t current simulation time
  !> \param numVic number of particles to be removed from the system

  if (lforw)  then 
     m_in=m_in+pVictim%mass
     res_t=(res_t*numVic+t-pVictim%partOrig)/(numVic+1)
     numVic=numVic+1
!!$     if(pVictim%segNum.EQ.2) then
!!$        mdry=mdry+pVictim%mass
!!$     elseif(pVictim%segNum.EQ.3)then
!!$        mwet=mwet+pVictim%mass
!!$     endif
  endif

  !pVictim_prev => pVictim%prev

  IF (ASSOCIATED(pVictim%prev)) THEN
     ! pVictim has a prev entry
     pVictim%prev%next => pVictim%next
     pVictim_prev => pVictim%prev
  ELSE  
     !if pVictim = pFirst -> there is no prev: set pFirst to next entry
     !pFirst => pVictim%next 
     pVictim_prev => pFirst%prev
     pFirst => pVictim%next 
  END IF

  IF (ASSOCIATED(pVictim%next)) THEN
     ! pVictim has a next entry
     pVictim%next%prev => pVictim%prev
  ELSE
     ! if last entry of the list -> set next pointer from prev to Null
     ! but only if it is not the only (very last) particle 
     IF (ASSOCIATED(pVictim%prev)) NULLIFY(pVictim%prev%next)
  END IF

  ! remove pVictim
  DEALLOCATE(pVictim)

END SUBROUTINE RemoveParticle
!***************************************************************************
  !> this routine calculates the mass which is taken up by a root segemnt
  !> and insert paricles with this mass inside the root (domain)
SUBROUTINE ParticleToRoot(pFirst,ipl,dt,t)
  USE TypeDef
  USE DoussanMat, ONLY: nsub,loc_Q,w_sub,cube_i,Joutr,isubmax
  USE RootData, ONLY: nrec, nbr, ibrseg, xs, ys, zs, irecpr, ngrow, ibrgrw, xg,yg,zg,segsur,seglen,ordseg,timorg,&
       age,nUrf,Urf,segsur,lSign,mhorm,irecsg,segconc,mcol, lno_root_growth, SegSoluteMass
  USE GridData, ONLY: dzgrid,dygrid,dxgrid,csink_cube, nElm, elmnod
  USE Paramdata,ONLY:  pi,lPartUp
  USE SoluteRootMat
  USE SolData, ONLY: theta, massR, sum_upt
  IMPLICIT NONE
  TYPE(Particle) ,POINTER :: P
  TYPE(Particle) ,POINTER ,INTENT(inout) :: pFirst
  INTEGER(ap) ::ibr,irecn,ifoln,igrow,iprvn,ci,ip,isub,ipl,numPa, cor(8), iE, counter
  REAL(dp) :: xp, posi, m_red, mass_adv, mass_diff, VoxVol, massroot, balance, IniMass(nElm),massInsert, m_pot_elm(nElm), m_pot
  REAL(sp):: dt,t
  LOGICAL :: n_apex,run, corr_root, corr_soil  
  INTEGER(ap)::massex_pot(nrec*isubmax,2)     !tabel for potential mass exchange for each subsegment nrec*isubmax
  REAL(dp)::store_mpot(nrec*isubmax), corr_upt(nElm)
  !> \param pFirst pointer to first particle in particlelist 
  !> \param ipl current root system (not used)
  !> \param dt time step size
  !> \param t current simulation time
massex_pot=0._dp
store_mpot=0._dp
sum_upt=0._dp
seg_upt=0._dp
corr_upt=0._dp


!print*, massex_pot

  VoxVol=dzgrid*dygrid*dxgrid
  m_red  = 0._dp
  mass_adv  = 0._dp
  mass_diff = 0._dp
  mass_elm = massR
  m_upt = 0._dp
  corr_root = .FALSE.
  corr_soil = .FALSE.
  massInsert=0._dp
  m_pot_elm=0._dp
  counter=1




  IF(.NOT.lno_root_growth)THEN
     CALL UpdateRootSys
  ENDIF


  !Soil Voxel water content for solute concentration
  do  iE=1,nElm
     cor(1)=elmnod(1,iE)
     cor(2)=elmnod(2,iE)
     cor(3)=elmnod(3,iE)
     cor(4)=elmnod(4,iE)
     cor(5)=elmnod(5,iE)
     cor(6)=elmnod(6,iE)
     cor(7)=elmnod(7,iE)
     cor(8)=elmnod(8,iE)

     theta_elm(iE) = (theta(cor(1))+theta(cor(2))+theta(cor(3))+theta(cor(4))+theta(cor(5))+&
          theta(cor(6))+theta(cor(7))+theta(cor(8)))/8
     IniMass(iE) = mass_elm(iE)
     ccube(iE) = mass_elm(iE)/(VoxVol*theta_elm(iE))

  end do

  ! Get mass in root segments
  ALLOCATE(P)
  NULLIFY(P)

  P=> firstP     !pointer setzen aus solute root mat
  segSoluteMass=0._dp  
  DO WHILE(ASSOCIATED(P))
     segSoluteMass(P%segNum) = segSoluteMass(P%segNum) + P%Mass   !segsolutemass länge von maxseg
     P => P%NEXT
  END DO



  !numPa = number of particles which sum represets the input mass
  numPa = 50  !put this in input
  IF(.NOT.lSign) THEN
     DO ibr=1,nbr
        n_apex=.FALSE.
        !* find the tip segment of the branch 'ibr'
        irecn=nrec
        DO WHILE (ibrseg(irecn).NE.ibr)
           irecn=irecn-1
        END DO
        !* the first one we find is an apex
        IF (seglen(irecn)<1.E-20) THEN ! skip this segment too small to be taken
           ifoln=irecn ! following node ID
           irecn=irecpr(irecn) !current node ID
        ELSE
           n_apex=.TRUE.!ifoln does not exist
        ENDIF
        IF (irecn==0) THEN !there exists a branch ibr but not yet any segment!
           run=.FALSE.
        ELSE
           run=.TRUE.
        ENDIF
        !* then the rest of the branch up to the seed or the embranchment
        DO WHILE (run)
           !* "upper" node
           iprvn=irecpr(irecn)                 
           ! if reached the seed or the embranchement => change branch
           IF (iprvn.EQ.0 .OR. (ibrseg(iprvn).NE.ibrseg(irecn))) run=.FALSE.			

           xp = 0.
           DO isub=1,nsub(irecn,ipl)
              ci=cube_i(irecn,isub,ipl) !element number

              !IF ((Joutr(irecn).NE.0).AND.(ccube(ci).NE.0)) THEN

              IF (isub.EQ.1) THEN
                 xp = 0.                  
              ELSE
                 xp = xp + seglen(irecn)*w_sub(irecn,isub-1,ipl)
              END IF
            
              if(Joutr(irecn).ge.0) then
                 mass_adv = frac_pass(irecn)*Joutr(irecn)*ccube(ci)*w_sub(irecn,isub,ipl)*dt
              elseif (Joutr(irecn).lt.0) then     !mass is removed from the root with advective transport
                 mass_adv = frac_pass(irecn)*Joutr(irecn)*(segconc(irecn)-segsorb(irecn))*w_sub(irecn,isub,ipl)*dt
              end if

              mass_diff = (ccube(ci)-segconc(irecn)*retard(irecn))*segsur(irecn)* & 	
                   w_sub(irecn,isub,ipl)*Perm(irecn)*dt			!ccube is in [M/L³water]  
              !segconc and segsorb is in [M/L³tot] -> divide by theta_R to receive cw=[M/L³wat]

             m_pot=mass_adv+mass_diff             !potential uptake for subsegment
!rint*, size(massex_pot), 'size massex'
              massex_pot(counter,1)=irecn
              massex_pot(counter,2)=ci
              store_mpot(counter)=m_pot

              counter= counter+1

              m_red = mass_adv + mass_diff

              !print*,'nochmal mupt',m_upt(ci),'ci',ci,'m_red',m_red
              mupt_adv(ci) = mupt_adv(ci) + mass_adv
              mupt_diff(ci) = mupt_diff(ci) + mass_diff
              m_red = 0._dp
              mass_adv = 0._dp
              mass_diff = 0._dp

           END DO !sub segments



           ! definition for the next run of the loop
           ifoln=irecn
           irecn=iprvn
           !* from here, not an apex
           n_apex=.FALSE.
        END DO !nodes 
     END DO ! branches  



! potential uptake is calculated
! now check if m_pot<=m_elm
    
     
     do counter=1,size(store_mpot)
        iE=massex_pot(counter,2)
        if(iE.gt.0)then
           m_pot_elm(iE)=m_pot_elm(iE)+store_mpot(counter)
          ! print*,counter, 'counter',  iE, 'iE', m_pot_elm(iE), 'm_pot_elm', store_mpot(counter), 'mpot', mass_elm(iE)
        end if
     end do


! 2 cases for a positive potential uptake: uptake < elementmass, uptake > elementmass
     do iE=1, nElm

!print*, m_pot_Elm(iE), 'potential uptake', mass_elm(iE), 'mass elm'

        if (m_pot_elm(iE).gt.1e-22)then
 !PRINT*, 'potential uptake is feasible'
           if (m_pot_elm(iE).le.mass_elm(iE)) then
              !Insert particles!
              do counter=1,size(store_mpot)
                 if (massex_pot(counter,2).eq.iE)then
                    irecn=massex_pot(counter,1)
                    seg_upt(irecn)=seg_upt(irecn)+store_mpot(counter)
                 end if
              end do

            !  print*, 'insert particles',  seg_upt(1:nrec), 'seg_upt'
              m_upt(iE)=m_upt(iE)+m_pot_elm(iE)

           else if (m_pot_elm(iE).gt.mass_elm(iE)) then
              !scale uptake mass
              print*, 'SCALE UPTAKE'
              do counter=1,size(store_mpot)
                 if (massex_pot(counter,2).eq.iE) then
                    irecn=massex_pot(counter,1)
                    seg_upt(irecn)=seg_upt(irecn)+(mass_elm(iE)/m_pot_elm(iE)*store_mpot(counter))
                 end if
                 m_pot_elm(iE)=mass_elm(iE)
              end do
              m_upt(iE)=m_upt(iE)+mass_elm(iE) 
           end if
        end if

! if mass is exudated from the root, seg_upt decreases, m_upt increasees
! check that only as much mass is removed from the root as is present in the root
corr_upt=0._dp
        if (m_pot_elm(iE).lt.-1e-22)then
        !   PRINT*, 'REMOVE ROOT MASS'
           do counter=1,size(store_mpot)
              if (massex_pot(counter,2).eq.iE) then
                 irecn=massex_pot(counter,1)
                 seg_upt(irecn)=seg_upt(irecn)+store_mpot(counter)
                 
                 if (abs(seg_upt(irecn)).gt.segSoluteMass(irecn)) then
                    seg_upt(irecn)=-segSoluteMass(irecn)
                    corr_upt(iE)=corr_upt(iE)+abs(m_pot_elm(iE))+seg_upt(irecn)!(abs(seg_upt(irecn))-(abs(segSoluteMass(irecn))))
                    print*, 'CORRECT UPTAKE'   ! updated mass needs to be tranferred to soil
                 end if

              end if
              ! m_pot_elm(iE)=mass_elm(iE)
           end do
           m_upt(iE)=m_upt(iE)+m_pot_elm(iE)+corr_upt(iE)   !!!! check this
         
          
        end if

     end do


!print*, sum(seg_upt), 'sum_seg_upt', sum(store_mpot), 'm_pot', sum(m_pot_elm), 'sum_elm', sum(m_upt), 'm_upt'

!!! insert particles



!!$              ! insert particles in root

DO irecn=1,nrec
   IF (seg_upt(irecn).gt.1e-22) THEN

   !   print*, 'insertParticle'
      DO ip=1,numPa
         ! position of particle in root segement
         posi = (ip-1)*(seglen(irecn)*w_sub(irecn,isub,ipl)/numPa)+xp
         CALL InsertParticle(pFirst,irecn,posi,seg_upt(irecn)/numPa,t)
         massInsert = massInsert+m_red/numPa
         totalParticleNum =  totalParticleNum +1

      END DO

   ELSE IF  (seg_upt(irecn).le.-1e-22) THEN   ! mass leaves the root segment
      print*, 'ReduceMass'
      IF (segSoluteMass(irecn).LE.abs(seg_upt(irecn)))  corr_root=.TRUE.

      CALL ReduceParticleMass(seg_upt(irecn), irecn, corr_root)

   END IF
  
END DO 


! Sum mass uptake of root   -> written in log1 file

sum_upt=sum(m_upt)
!print*, 'sum_upt', sum_upt


     !Calculate fact to submit to ParTrace

     DO iE=1, nElm
        fact(iE)=1._dp
        if(abs(m_upt(iE)).GE.1e-22) THEN
           IF(IniMass(iE).LE.1e-22) THEN
              fact(iE)=1._dp 
           ELSE
              fact(iE) = (1-m_upt(iE)/IniMass(iE))**timestepfactor
           END IF
        ELSE
           fact(iE)=1._dp 
        END IF
     END DO

  ELSE						! hormornal signalling
     DO irecn=1,nrec
        if(mhorm(irecn,ipl).GT.1E-22) THEN
           DO isub=1,nsub(irecn,ipl)  
              IF (isub.EQ.1) THEN
                 xp = 0.                  
              ELSE
                 xp = xp + seglen(irecn)*w_sub(irecn,isub-1,ipl)
              END IF
              DO ip=1,numPa
                 ! position of particle in root segement
                 posi = (ip-1)*(seglen(irecn)*w_sub(irecn,isub,ipl)/numPa)+xp
                 CALL InsertParticle(pFirst,irecn,posi,mhorm(irecn,ipl)/numPa,t)
		 totalParticleNum =  totalParticleNum +1	
              END DO
           END DO
        end if
     END DO

  END IF
  ALLOCATE(P)
  NULLIFY(P)

  P=> firstP     !pointer setzen aus solute root mat
  segSoluteMass=0._dp  
  DO WHILE(ASSOCIATED(P))
     segSoluteMass(P%segNum) = segSoluteMass(P%segNum) + P%Mass   !segsolutemass länge von maxseg
     P => P%NEXT
  END DO



  massroot=0.
  do irecn=1, nrec
     massroot=massroot+segSolutemass(irecn)
  end do

END SUBROUTINE ParticleToRoot
!***************************************************************************

!> main subroutine of solute transport inside root
!> calculates for one time step the input and the moveing 
!> of the particles
SUBROUTINE SoluteRoot(iCount,dt,t)
  USE TypeDef
  USE SolData, ONLY: conc,concR
  USE GridData, ONLY: elmnod, nElm
  USE ParamData, ONLY:  maxParticle, pi, lPartUp
  USE RootData, ONLY: seglen,segsur,segconc,segSoluteMass,nrec,lSign,&
       m_in,mcol,vol_buff,concol,irecpr,res_t,segvol,vol_root
  USE SoluteRootMat
  USE DoussanMat, ONLY: axialRootFlow
  USE PlntData, ONLY: Tact
  IMPLICIT NONE

  INTEGER(ap),INTENT(in) :: icount
  REAL(sp),INTENT(in) :: dt,t
  TYPE(Particle) ,POINTER  :: P
  INTEGER(ap) :: irec, corner(8), iE, aa, numVic
  INTEGER(ap), SAVE :: kk=0
  !> \param icount icount=0 at first time step
  !> \param dt time step size
  !> \param t current simulation time

!!$  
!!$  !calculate concentration per cube if solute transport is calculated based on FEM solution
!!$  IF(.NOT.lSign .AND. lPartUp) THEN
!!$     ALLOCATE(ccube(nElm))
!!$     ccube = 0.
!!$     !DO iE=1,nElm
!!$        !corner(:)=elmnod(:,iE)
!!$        !ccube(iE) = SUM(conc(corner))/8. 
!!$     !END DO
!!$     ccube = concR
!!$  END IF

  ! calculate retardation factor for sorption
  if(l_linSorb .or. l_freundSorb) CALL CalculateRootRetard

  ! first move particles, then insert new particles
  m_in=0._dp
!!$  mdry=0._dp
!!$  mwet=0._dp
  numVic = 0
  res_t = 0._dp
 ! IF (lPartUp) vol_buff = vol_root	 !vol_root !buffer volume == total root volume
  IF ((icount.EQ.0) .OR. (totalParticleNum .EQ. 0)) NULLIFY(firstP)
  IF (ASSOCIATED(firstP)) CALL MoveParticles(firstP,dt,icount,numVic,t)
  CALL ParticleToRoot(firstP,1,dt,t)
     
  ALLOCATE(p)
  NULLIFY(p)

  ! calculate mass per root segment
  p=> firstP
  segSoluteMass = 0.
  DO WHILE(ASSOCIATED(P))
     segSoluteMass(P%segNum) = segSoluteMass(P%segNum) + P%Mass
     P => P%NEXT
  END DO

  !calculate concentration in root segment
  DO irec=1,nrec
     IF (seglen(irec) .GE. 1.E-20) THEN ! skip this segment too small to be taken
        segconc(irec) = segSoluteMass(irec)/segvol(irec)  !segconc in [M/L³tot]
     END IF
  END DO


  ! first root segment acts as buffer
  if(lPartUp) mcol = mcol + m_in !mass gain of buffer

  ! calculate cumulated conc. at root collar (for the last time step)

  if(lSign) then 
     concol = 0._dp ! calculate  conc. at root collar (for the last time step)
     concol = mcol/vol_buff
     mcol = max(mcol - concol*Tact(1)*dt,0._dp) !mass loss of buffer
  end if

END SUBROUTINE SoluteRoot
!***************************************************************************
  !> this routine stores a list of the connection of root segments
  SUBROUTINE MakeFollowList
  USE TypeDef
  USE RootData, ONLY: nrec, nbr, ibrseg, irecpr, ibrseg, seglen
  USE SoluteRootMat, ONLY: irecfollow,numfollow
  IMPLICIT NONE

  INTEGER(ap) ::ibr,irecn,ifoln,iprvn,a
  LOGICAL :: n_apex,run
  irecFollow = 0
  numFollow = 0

  DO ibr=1,nbr
     n_apex=.FALSE.
     !* find the tip segment of the branch 'ibr'
     irecn=nrec
     DO WHILE (ibrseg(irecn).NE.ibr)
        irecn=irecn-1
     END DO
     !* the first one we find is an apex
     IF (seglen(irecn)<1.E-20) THEN ! skip this segment too small to be taken
        ifoln=irecn ! following node ID
        irecn=irecpr(irecn) !current node ID
     ELSE
        n_apex=.TRUE.!ifoln does not exist
     ENDIF
     IF (irecn==0) THEN !there exists a branch ibr but not yet any segment! 
        run=.FALSE.
     ELSE
        run=.TRUE.
     ENDIF
     !* then the rest of the branch up to the seed or the embranchment
     DO WHILE (run)
        !* "upper" node
        iprvn=irecpr(irecn)                 
        ! if reached the seed or the embranchement => change branch
        IF (iprvn==0.OR.(ibrseg(iprvn).NE.ibrseg(irecn)))THEN !run=.FALSE.		
           EXIT
        ENDIF
        !number of following root segments (branches)
        a = numFollow(irecpr(irecn))+1
        !list which root segment follows
        irecfollow(irecpr(irecn),a) = irecn
        numFollow(irecpr(irecn)) = a
        ! definition for the next run of the loop
        ifoln=irecn
        irecn=iprvn
        !* from here, not an apex
        n_apex=.FALSE.
     END DO !nodes 
  END DO ! branches

END SUBROUTINE MakeFollowList
!***************************************************************************
!> calculates a retardation factor for the solute particles in case of sorption
SUBROUTINE CalculateRootRetard
  USE typedef
  USE RootData, ONLY: nrec,segconc
  USE SoluteRootMat, ONLY: retard, sorp, l_linSorb, l_freundSorb, segsorb , rho_R, theta_R, segsolv

  IMPLICIT NONE
  
  REAL(dp) :: epsilon = 1e-5, cw, r, c, rdata, k, n, f, ff,itsum, itold
  INTEGER(ap) :: ii,irecn
  
  ! calculates retardation within roots due to sorption
  ! assumes fully saturated roots; kd is not scaled with water content
  ! references in Partrace / particleclassSorp.cpp and dissertation of
  ! O. Neuendorf
!theta_R = 1
!rho_R = 1 

  if(l_linSorb) then
	do irecn=1,nrec 
     k = sorp(1)
     retard(irecn) = 1/(theta_R + rho_R*k)
     segsorb(irecn)=(1-retard(irecn))*segconc(irecn)

    end do 
	


 elseif(l_freundSorb) then

     do irecn=1,nrec
        c = segconc(irecn)			! in [M/L³tot]
        rdata = retard(irecn)
        k = sorp(1) !/rho_R*theta_R
        n = sorp(2)
		

!! change to water volume based concentrations

	c=c/theta_R				!in [M/L³wat]
!	print*, c, rdata, retard(irecn)

        if (rdata > 0._dp) then
           cw = rdata * c 
        else
           cw = 0.5 * c 
        end if

        do ii=1,50 !max. 50 iterations
           r = cw-(cw+k*cw**n-c)/(1._dp+n*k*cw**(n-1));
           
           if(r<0) then
              cw = cw/2
           else              
              if(abs((cw-r)/cw)<epsilon) then
                 cw = r
                 goto 10
              else
                 cw = r
              end if
           end if
        end do ! iterations
        
10      if(ii==50) then
           rdata = 1._dp
           
        else if  (c==0._dp) then
           rdata = 1._dp 

        else
           rdata = cw/c					!rdata in [-] (M/L³wat)
        end if

        retard(irecn) = rdata*theta_R			!retard in [-] [M/L³tot]
      
        segsorb(irecn) = (c-cw)*theta_R 		!segsorb in [M/L³tot]  
	segsolv(irecn)= cw				!segsolv in [M/L³tot]
		
     end do ! loop over nrec

end if ! sorption type


  
END SUBROUTINE CalculateRootRetard

!******************************
SUBROUTINE UpdateRootSys
  USE typedef
  USE RootData, ONLY: seglen, segsur, nrec, segvol, vol_root
  USE ParamData, ONLY: Pi


  IMPLICIT NONE

  REAL(dp) :: crossSectionSeg(nrec) , segrad(nrec)
  INTEGER(ap) :: ii,irec

!When root growth, update segment information:
!- segment volume
!- buffer volume



DO irec=1,nrec
segrad(irec) = segsur(irec)/2._dp/pi/seglen(irec)
crossSectionSeg(irec) = segrad(irec)**2*pi
segvol(irec) =  segrad(irec)**2*pi*seglen(irec) !root volume per root segment
vol_root = vol_root +  segvol(irec) !total root volume
END DO


END SUBROUTINE UpdateRootSys
!***************************************************************************
SUBROUTINE ReduceParticleMass(m_red, irecn, corr_root)
  USE typedef
  USE RootData, ONLY:segSoluteMass, nrec
  USE SoluteRootmat, ONLY:Particle, totalParticleNum, firstP

  IMPLICIT NONE
 
  !TYPE(Particle) ,POINTER ,INTENT(inout) ::  Pfirst
  TYPE(Particle) ,POINTER :: P
  REAL(dp)::red_fact(nrec) , m_red
  INTEGER(ap):: irecn,  num
  LOGICAL:: corr_root
  

!Calculate Reduction Factor red_fact 
  red_fact=1._dp


  ALLOCATE(P)
  NULLIFY(P)

  P => firstP   
  segSoluteMass=0._dp  
  DO WHILE(ASSOCIATED(P))
     segSoluteMass(P%segNum) = segSoluteMass(P%segNum) + P%Mass   !segsolutemass länge von maxseg
     P => P%NEXT

  END DO

  IF (corr_root) THEN
     red_fact(irecn)=0._dp
  ELSE IF (segSoluteMass(irecn).GE.1e-22) THEN
     red_fact(irecn)=1+(m_red/segSoluteMass(irecn))

  ELSEIF(segSoluteMass(irecn).LT.1e-22) THEN
     red_fact(irecn)=1.0_dp  ! produces mass!
 
  ENDIF

  ALLOCATE(P)
  NULLIFY(P)
  P=>firstP

  segSoluteMass=0._dp  
  Do While(associated(P))
     ! print*, P%Mass, 'particle mass', P%segNum, 'ParticleSegNum', red_fact(P%segNum), 'red_fact'
     P%Mass=P%Mass*red_fact(P%segNum)
     segSoluteMass(P%segNum) = segSoluteMass(P%segNum) + P%Mass
     P => P%Next
  end do
!print*,'reduce2',' segsolutemass',segsolutemass(irecn),'irecn',irecn
if(red_fact(irecn).eq.1) stop 

END SUBROUTINE ReduceParticleMass
!********************************************************************************
END MODULE ParticlesInRoot
