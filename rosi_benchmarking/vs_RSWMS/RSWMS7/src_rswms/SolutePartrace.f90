! ==============================================================================
! Source file SolutePartrace |||||||||||||||||||||||||||||||||||||||||||||||||||
! ==============================================================================
!> \brief rearrange the velocity and the water content (theta) 
!> 		  arrays from R-SWMS into the array ordering from PARTRACE
!> 		  R-SWMS starts with the upper left front node
!> 		  PARTRACE starts with the lower left front node
SUBROUTINE REARRANGE_RSWMS_TO_PARTRACE(velocity_P,theta_P,iNodes)
!> \param velocity_P velocity array needed by PARTRACE
!> \param theta_P water content array needed by PARTRACE
!> \param iNodes number of grid nodes 
  USE GridData, ONLY: continu,nx,ny,nz,nPt
  USE ParamData, ONLY: maxnod
  USE SolData, ONLY: Vx,Vy,Vz,theta, concR
  USE DoussanMat, ONLY: curr_BCtp
  USE typedef
  IMPLICIT NONE
  
  !REAL(sp), INTENT (inout):: sinkW(maxnod),sinkC(maxnod)
  REAL(dp),INTENT (out):: theta_P(iNodes)
  REAL(dp), INTENT (out)::velocity_P(iNodes,3)
  INTEGER(ap),  INTENT (in):: iNodes
  REAL(sp):: sinkW_r(nPt),sinkC_r(nPt),velo_r(nPt,3),theta_r(nPt)
  INTEGER(ap):: i,j,k,l,ix,iy,iz

  ! nessecary because Partrace orders elements and nodes different 
  ! from R-SWMS
  
  CALL Veloc
  velo_r(:,1) = Vx/theta
  velo_r(:,2) = Vy/theta
  velo_r(:,3) = Vz/theta
  theta_r = theta
   
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
              velocity_P(j,:) = velo_r(k,:)
              theta_P(j) = theta_r(k)
           END DO
        END DO
     END DO
     
  ELSE
     l=0
     DO j=1,nz
        DO 2 i=1,nx*ny
           k = nx*ny*nz - j*nx*ny +i
           l = l+1
           velocity_P(l,:) = velo_r(k,:)
           theta_P(l) = theta_r(k)
2       END DO
     END DO
  ENDIF

END SUBROUTINE REARRANGE_RSWMS_TO_PARTRACE
!****************************************************************
! soil element based variables from partrace to rswms order
SUBROUTINE Rearrange_p2r(var_P,var_R)

  USE GridData, ONLY: nex,ney,nez,nElm
  USE typedef
  USE Soldata
  IMPLICIT NONE

  REAL(dp), INTENT (in):: var_P(nElm)
  REAL(dp), INTENT (out):: var_R(nElm)
  INTEGER(ap):: i,j,k,l,ix,iy,iz

  ! nessecary because Partrace orders elements and nodes different 
  ! from R-SWMS
  ! rearrange concentration from partrace to rswms ordering

  l=0
  DO j=1,nez
     DO  i=1,nex*ney
        k = nex*ney*nez - j*nex*ney +i
        l = l+1
        var_R(l) = var_P(k)
     END DO
  END DO
  
END SUBROUTINE Rearrange_p2r
!****************************************************************
! soil element based variables from rswms to partrace order
SUBROUTINE Rearrange_r2p(var_R,var_P)

  USE GridData, ONLY: nex,ney,nez,nElm
  USE typedef
  USE Soldata
  IMPLICIT NONE
  
  REAL(dp), INTENT (in):: var_R(:)
  REAL(dp), INTENT (out):: var_P(:)
  INTEGER(ap):: i,j,k,l,ix,iy,iz

  ! nessecary because Partrace orders elements and nodes different 
  ! from R-SWMS

  l=0
  DO j=1,nez
     DO  i=1,nex*ney
        k = nex*ney*nez - j*nex*ney +i
        l = l+1
        var_P(k) = var_R(l)
     END DO
  END DO
  
END SUBROUTINE Rearrange_r2p
!****************************************************************
SUBROUTINE OutVTKCouple(kout)
  USE typedef
  USE ParamData, ONLY: lsalinity,lPartUp
  USE DoussanMat, ONLY: PHs_osmotic
  USE SoluteRootMat, ONLY: m_upt,ccube,mupt_adv,mupt_diff
  USE Soldata
  USE GridData
  IMPLICIT NONE
  INTEGER(ap), INTENT(in)::kout
  INTEGER(ap):: i, k, j, ix, iy, iz
  CHARACTER file*22
  WRITE (file,'(A16)')'out/vtk/velociC.'

  IF (kout.LT.10) THEN
     WRITE (file(16:16),'(I1)') kout
     WRITE (file(17:20),'(A4)') '.vtk'
  ELSEIF (kout.LT.100) THEN
     WRITE (file(16:17),'(I2)') kout
     WRITE (file(18:21),'(A4)') '.vtk'
  ELSE
     WRITE (file(16:18),'(I3)') kout
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

     WRITE (8,'(/A15,1X,I7)')'POINT_DATA', (nx+1)*(ny+1)*nz
     WRITE (8,'(A24)')'SCALARS velocity float 3'
     WRITE (8,'(A20)')'LOOKUP_TABLE default'

     CALL Veloc(hNew)

     k=0
     j=0
     DO  iz=1,nz
        DO  iy=1,ny+1
           DO  ix=1,nx+1     
              k = ix + (iy-1)*nx + (iz-1)*nx*nx
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

     WRITE (8,'(A30)')'SCALARS pressurehead float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     k=0
     j=0
     DO  iz=1,nz
        DO  iy=1,ny+1
           DO  ix=1,nx+1     
              k = ix + (iy-1)*nx + (iz-1)*nx*nx
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
     WRITE (8,'(A30)')'SCALARS wc float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     k=0
     j=0
     DO  iz=1,nz
        DO  iy=1,ny+1
           DO  ix=1,nx+1     
              k = ix + (iy-1)*nx + (iz-1)*nx*nx
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

     WRITE (8,'(A30)')'SCALARS conc float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     k=0
     j=0
     DO  iz=1,nz
        DO  iy=1,ny+1
           DO  ix=1,nx+1     
              k = ix + (iy-1)*nx + (iz-1)*nx*nx
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
!!$     WRITE (8,*)''
!!$     WRITE (8,'(A30)')'SCALARS sink float '
!!$     WRITE (8,'(A20)')'LOOKUP_TABLE default'
!!$     k=0
!!$     j=0
!!$     DO  iz=1,nz
!!$        DO  iy=1,ny+1
!!$           DO  ix=1,nx+1     
!!$              k = ix + (iy-1)*nx + (iz-1)*nx*nx
!!$              j = j +1
!!$              IF (ix.EQ.nx+1) THEN
!!$                 k = k - nx
!!$              END IF
!!$              IF  (iy.EQ.ny+1) THEN
!!$                 k = k - (nx*ny)
!!$              END IF
!!$              WRITE (8,'(1X,1pE11.4)',advance="no")sink(k)
!!$           END DO
!!$        END DO
!!$     END DO
     WRITE (8,*)''
     WRITE (8,'(A30)')'SCALARS csink float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     k=0
     j=0
     DO  iz=1,nz
        DO  iy=1,ny+1
           DO  ix=1,nx+1     
              k = ix + (iy-1)*nx + (iz-1)*nx*nx
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
     WRITE (8,*)''

     CALL Veloc(hNew)
     DO 2 i=1,nPt
        WRITE (8,'(3(1X,1pE11.4))')Vx(i),Vy(i),Vz(i)
2    END DO

     WRITE (8,'(A30)')'SCALARS pressurehead float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE11.4)',advance="no") hnew(1:nPt)
     WRITE (8,*)''

     WRITE (8,'(A30)')'SCALARS wc float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE11.4)',advance="no") theta(1:nPt)
     WRITE (8,*)''

     WRITE (8,'(A30)')'SCALARS conc float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE11.4)',advance="no") conc(1:nPt)
     WRITE (8,*)''

     WRITE (8,'(A30)')'SCALARS csink float '
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE11.4)',advance="no") csink(1:nPt)

  ENDIF
  WRITE (8,*)''
  WRITE (8,'(/A15,1X,I7)')'CELL_DATA', nElm
  WRITE (8,'(A24)')'SCALARS sinkElm float'
  WRITE (8,'(A20)')'LOOKUP_TABLE default'
  WRITE (8,'(1X,1pE11.4)',advance="no") sink_cube(1:nElm)
  WRITE (8,*)''

  WRITE (8,*)''
  !WRITE (8,'(/A15,1X,I7)')'CELL_DATA', nElm
  WRITE (8,'(A24)')'SCALARS concR float'
  WRITE (8,'(A20)')'LOOKUP_TABLE default'
  WRITE (8,'(1X,1pE13.4)',advance="no") concR(1:nElm)
  WRITE (8,*)''

    

  IF(lPartUp) THEN
     WRITE (8,'(A24)')'SCALARS ccube float'
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE13.4)',advance="no") ccube(1:nElm)

     WRITE (8,'(A24)')'SCALARS m_upt float'
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE13.4)',advance="no") m_upt(1:nElm)
     WRITE (8,*)''

     WRITE (8,'(A24)')'SCALARS mupt_adv float'
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE13.4)',advance="no") mupt_adv(1:nElm)
     WRITE (8,*)''

     WRITE (8,'(A24)')'SCALARS mupt_diff float'
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE13.4)',advance="no") mupt_diff(1:nElm)
     WRITE (8,*)''

  END IF

  
  IF (lsalinity) THEN
     WRITE (8,*)''
     WRITE (8,'(A27)')'SCALARS osmoticHead float'
     WRITE (8,'(A20)')'LOOKUP_TABLE default'
     WRITE (8,'(1X,1pE11.4)',advance="no") PHs_osmotic(1:nElm)
     WRITE (8,*)''
  END IF
  
  CLOSE (8)
  RETURN
END SUBROUTINE OutVTKCouple
!************************************************************************
! transformation from concentration [micromol / cm^3] to osmotic head [cm]
! Fac = factor for transformation
SUBROUTINE Salinity
  USE typedef
  USE DoussanMat, ONLY: PHs_osmotic
  USE SolData, ONLY: ConcR
  IMPLICIT NONE
  Real(dp) :: Fac

  Fac = -50.
  Phs_osmotic = Fac * concR 
 
  
  RETURN
END SUBROUTINE Salinity
!********************************************************************************
