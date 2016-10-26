!automatioc grid file for the forward model
!gfortran -O0 -o grid Grid_SWMS.f90
!**** generates 3D-element mesh ******

program grid
  IMPLICIT NONE
  Integer, Parameter :: sp=selected_real_kind(4)
  Integer, Parameter :: dp=selected_real_kind(8)
  !**** generates 3D-element mesh ******
  LOGICAL continu, bilayered
  INTEGER i,j,k,l,m,n,o,p,ii,jj,geom,width,ICtype,MatNum,iCord
  INTEGER nex,ney,nez,nx,ny,nz,ax,by,np,nl,ne,ie,ix,iy,iz,iRest,layers,Xgradient
  REAL(dp) radius,zlim,htop,conctr,r,hinit,v1,v2,v3,v4,v5
  REAL(dp) x,y,z,dx,dy,dz,xmin,ymin,zmax,wx
  LOGICAL ::lchem,lcalloc,ltemp
  jj=2

  !* check soil.in or soiltab.in for number of soil material				
  OPEN (Unit=10,FILE='../in/soil.in',STATUS='OLD',ERR=20)
  GOTO 21
20 OPEN (Unit=10,FILE='../in/soiltab.in',STATUS='OLD',ERR=11)
21 DO ii=1,3
     READ (10,*)
  ENDDO
  READ (10,*) layers
  IF (layers .eq. 2) THEN
     bilayered=.TRUE.
  ELSEIF (layers .eq. 1) THEN
     bilayered=.FALSE.
  ELSE
     STOP 'The number of soil layers can currently only be (1) or (2) for the building of nodes.in'
  ENDIF
  CLOSE(10)
    
!* check control.in for domain continuity and solute transport						
  OPEN (Unit=10,FILE='../in/control.in',STATUS='OLD',ERR=14)
  DO ii=1,50
     READ (10,*)
  ENDDO
  READ (10,*) lcalloc,lChem,ltemp,continu
  if(lchem)  WRITE (*,'(/'' Simulation will not include solute transport. (lChem=False)'')')
  IF (.not.continu) THEN
     WRITE (*,'(/''No continous soil domain, continu=FALSE'')')
  ELSEIF (continu) THEN
     WRITE (*,'(/'' Continous soil domain, continu=TRUE'')')
  ELSE
     STOP 'Set parameter for a continuous domain simulation in control.in to  ON (1) or OFF (2)'
  ENDIF
  CLOSE(10)
				
  OPEN (Unit=10,FILE='../in/mesh.in',STATUS='OLD',ERR=13)
  READ (10,*)
  READ (10,*)
  READ (10,*)
  READ (10,*) dx, dy, dz,nex, ney, nez,xmin, ymin, zmax
  READ (10,*)
  READ (10,*)
  READ (10,*) ICtype, htop
  READ (10,*)
  READ (10,*)
  READ (10,*) geom
  READ (10,*)
  SELECT CASE(geom)
  CASE(0)
     READ (10,*)
  CASE(1)
     READ (10,*) zlim 
  CASE(2)
     READ (10,*) width
  CASE(3)
     READ (10,*) radius  
  CASE DEFAULT
     READ (10,*)
  END SELECT
  READ (10,*)
  READ (10,*)
  IF (lchem) THEN
     READ (10,*) conctr   
  ELSE
     conctr=0
  END IF
  CLOSE(10)

    
  
  !* number of nodes in each direction:
  IF (continu) THEN
     nx=nex
     ny=ney
     nz=nez+1
  ELSE
     nx=nex+1
     ny=ney+1
     nz=nez+1
  ENDIF


  !* number of points per z-layer:
  nl=nx*ny

  !* total number of nodal points and of elements:
  np=nl*nz
  ne=nex*ney*nez
  
  OPEN (UNIT=10,FILE='../in/grid.in',STATUS='UNKNOWN')
  WRITE (10,'(''***** GRID ELEMENTS INFORMATION *****'')')
  WRITE (10,'(/''     nPt    nElm'')')
  WRITE (10,'(2(3X,I8))')np,ne
  WRITE (10,'(/''     nx     ny     nz'',  ''     nex     ney     nez'',  ''       dx          dy          dz'')')
  WRITE (10,'(6(3X,I5),3(1X,1pE11.4))') nx,ny,nz,nex,ney,nez,dx,dy,dz
  WRITE (10,'(/''* Element Information *'')')
  WRITE (10,'(''     iE       i     j     k     l     m     n    o    p   sub '',''    -----A-----  --C--'')')
  
  !* element nodes:
  ie=0
  do iz=1,nez
     do iy=1,ney
        do ix=1,nex
           ie=ie+1
           i=ix+nx*(iy-1)+nx*ny*(iz-1)
           j=i+1
           if ((continu).and.(ix.eq.nex)) j=j-nex
           k=i+nx
           if ((continu).and.(iy.eq.ney)) k=k-ney*nex
           l=k+1
           if ((continu).and.(ix.eq.nex)) l=l-nex
           
           m=i+nx*ny
           n=m+1
           if ((continu).and.(ix.eq.nex)) n=n-nex
           o=m+nx
           if ((continu).and.(iy.eq.ney)) o=o-ney*nex
           p=o+1
           if((continu).and. (ix.eq.nex)) p=p-nex
           
           iCord=ix+iy+iz
           
           iRest=mod(iCord,jj)
           if(iRest.eq.1) then
              WRITE (10,'(I7,3X,9(I7,1X),3X,6I2,1X,3I2)')ie,i,j,k,l,m,n,o,p,1,1,1,1,0,0,0,1,1,1
              
           else
              WRITE (10,'(I7,3X,9(I7,1X),3X,6I2,1X,3I2)')ie,i,j,k,l,m,n,o,p,2,1,1,1,0,0,0,1,1,1
           end if
        end do
     end do
  end do
  
  CLOSE (10)

  !* nodal coordinates:
  OPEN (UNIT=10,FILE='../in/nodes.in',STATUS='UNKNOWN')
  WRITE (10,'(//////''Node# Mater.#'',4X,''x'',11X,''y'',11X,''z'', 11X,''h'',11X,''C'')')

!************************************************************************
  i=0
  iz=0
  z=zmax+dz
100 iz=iz+1
  z=z-dz
  iy=0
  y=ymin-dy
110 iy=iy+1
  y=y+dy
  ix=0
  x=xmin-dx
120 ix=ix+1
  x=x+dx
  i=i+1
  
  IF (ICtype==1) THEN
     hinit=htop
  ELSEIF (ICtype==2) THEN
     Xgradient=0
     hinit=htop+zmax-z-Xgradient*cos((x*2*3.141592/(nex*dx))-(3.141592/2))*400 !linear initial PH distribution whatever the top z-value
  ELSE !(ICtype=3)
     OPEN (Unit=15,FILE='../SL_h1_T1/out_SL_h1_T1/outfem.134',STATUS='OLD',ERR=12)
     READ (15,*)
     READ (15,*)
     READ (15,*)
     READ (15,*)
     READ (15,*)
     READ (15,*)
     READ (15,*)
     READ (15,*)
     READ (15,*)
     READ (15,*) v1,v2,v3,v4,v5,hinit
  ENDIF
 
     IF (abs(x).LT.1.E-7) x = 0._sp
     IF (abs(y).LT.1.E-7) y = 0._sp
  !************************************************************************
  !Allocation of material properties to soil nodes
  MatNum=1
  IF (bilayered) THEN
     SELECT CASE(geom)
        CASE(0)
        CASE(1)
           IF(z.LT.zlim) MatNum=2
        CASE(2)
           wx=dx*width/2  
           IF((x.LT.wx) .AND. (x.GT.-wx)) MatNum=2 
        CASE(3)
           !mid point of x-y area
           ax = (nex*dx/2+xmin)
           by = (ney*dy/2+ymin)
           r=sqrt((x-ax)*(x-ax)+(y-by)*(y-by))
           if(r.ge.radius)  MatNum=2  
        END SELECT
     ENDIF

!************************************************************************
  WRITE (10,'(I7,6X,I2,3(1X,1pE11.4),1X,1pE11.4,1X,1pE11.4)') i,MatNum,x,y,z,hinit,Conctr
  IF (ix.LT.nx) GO TO 120
  IF (iy.LT.ny) GO TO 110
  IF (iz.LT.nz) GO TO 100

  

  CLOSE (10)
  CLOSE (15)

  
WRITE(*,'(//'' Grid was generated successfully.''///)')
RETURN

11 STOP 'Both of  < soil.in > and < soiltab.in >  files not found -- program terminated.'
12 STOP 'File  containing the initial conditions  not found -- program terminated.'
13 STOP 'File  <mesh.in>  not found -- program terminated.' 
14 STOP 'File  < control.in >  not found -- program terminated.'

END program grid
