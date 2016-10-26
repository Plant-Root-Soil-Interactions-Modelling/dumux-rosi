! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! Source file ORTHOFEM.FOR |||||||||||||||||||||||||||||||||||||||||||||
!
!                            ORTHOFEM
!
!                           VERSION 1.02
!
!                      FORTRAN SUBROUTINES FOR
!                   ORTHOMIN OR CONJUGATE GRADIENT
!               MATRIX SOLUTION ON FINITE-ELEMENT GRIDS
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                           CARL A. MENDOZA
!
!                      WITH CONTRIBUTIONS FROM:
!                            RENE THERRIEN
!
!                     BASED ON AN ORIGINAL CODE BY:
!                         FRANK W. LETNIOWSKI
!
!              WATERLOO CENTRE FOR GROUNDWATER RESEARCH
!                       UNIVERSITY OF WATERLOO
!                         WATERLOO, ONTARIO
!                          CANADA, N2L 3G1
!
!                    LATEST UPDATE: JANUARY 1991
!
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                   COPYRIGHT (c) 1989, 1990, 1991
!                            E.A. SUDICKY
!                            C.A. MENDOZA
!              WATERLOO CENTRE FOR GROUNDWATER RESEARCH
!
!          DUPLICATION OF THIS PROGRAM, OR ANY PART THEREOF,
!          WITHOUT THE EXPRESS PERMISSION OF THE COPYRIGHT
!          HOLDERS IS STRICTLY FORBIDDEN
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!               Modified for SWMS_3D code by Jirka Simunek
!                            august 1994
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                             DISCLAIMER
!
!     ALTHOUGH GREAT CARE HAS BEEN TAKEN IN PREPARING THIS CODE
!     AND THE ACCOMPANYING DOCUMENTATION, THE AUTHORS CANNOT BE
!     HELD RESPONSIBLE FOR ANY ERRORS OR OMISSIONS.  THE USER IS
!     EXPECTED TO BE FAMILIAR WITH THE FINITE-ELEMENT METHOD,
!     PRECONDITIONED ITERATIVE TECHNIQUES AND FORTRAN PROGRAMING.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!             A USER´S GUIDE IS AVAILABLE -> CONSULT IT!
!
!     THESE SUBROUTINES SOLVE A BANDED (OR SPARSE) MATRIX USING:
!       - PRECONDITIONED ORTHOMIN FOR ASYMMETRIC MATRICES, OR
!       - PRECONDITIONED CONJUGATE GRADIENT FOR SYMMETRIC MATRICES
!         (FULL MATRIX STORAGE REQUIRED)
!
!     PRECONDITIONING IS BY INCOMPLETE LOWER-UPPER DECOMPOSITION
!       - ONLY ONE FACTORIZATION (GAUSSIAN ELIMINATION) IS PERFORMED
!       - EQUIVALENT TO DKR FACTORIZATION
!
!     THE SUBROUTINES ARE DESIGNED FOR FINITE-ELEMENT GRIDS
!       - ARBITRARY ELEMENT SHAPES AND NUMBERING MAY BE USED
!         - NUMBERING MAY, HOWEVER, AFFECT EFFICIENCY
!           - TRY TO MINIMIZE THE BANDWIDTH AS MUCH AS POSSIBLE
!         - ALL ELEMENTS MUST HAVE THE SAME NUMBER OF LOCAL NODES
!
!
!     THE FOLLOWING ROUTINES ARE CALLED FROM THE SOURCE PROGRAM:
!       IADMAKE (IN,NN,NE,NLN,MNLN,MAXNB,IAD,IADN,IADD)
!         -> ASSEMBLE ADJACENCY MATRIX
!       FIND (I,J,K,NN,MAXNB,IAD,IADN)
!         -> LOCATE MATRIX POSITION FOR A NODAL PAIR (ASSEMBLY)
!       ILU (R,NN,MAXNB,IAD,IADN,IADD,B)
!         -> DECOMPOSE GLOBAL MATRIX
!       ORTHOMIN (R,C,GT,NNR,MAXNB,MAXNN,IAD,IADN,IADD,B,VRV,
!                 RES,RQI,RQ,Q,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,
!                 NORTH,MNORTH,MAXIT)
!         -> SOLVE DECOMPOSED MATRIX
!
!     THESE ROUTINES CALL OTHER ROUTINES (LOCATED DIRECTLY BELOW THE
!     APPROPRIATE PRIMARY ROUTINE IN THE CODE)
!
!     THE FOLLOWING ARRAYS MUST BE DEFINED IN THE SOURCE PROGRAM
!     (THESE ARRAYS ARE PASSED TO THE SOLVER SUBROUTINES):
!
!     IN(MNLN,MAXNE) - INCIDENCE MATRIX (ELEMENTAL NODE DEFINITION)
!
!     GT(MAXNN)      - RIGHT-HAND-SIDE VECTOR
!     C(MAXNN)       - SOLUTION VECTOR
!     R(MAXNB,MAXNN) - GLOBAL MATRIX TO BE SOLVED
!
!     ARRAY DIMENSIONING PARAMETERS
!
!     MAXNN  - MAXIMUM NUMBER OF NODES
!     MAXNE  - MAXIMUM NUMBER OF ELEMENTS
!     MNLN   - MAXIMUM NUMBER OF LOCAL NODES (IN AN ELEMENT)
!     MAXNB  - MAXIMUM NUMBER OF NODES ADJACENT TO A PARTICULAR NODE
!              (INCLUDING ITSELF).
!            - IE. THE MAXIMUM NUMBER OF INDEPENDENT NODES THAT A
!              PARTICULAR NODE SHARES AN ELEMENT WITH.
!            - THIS WILL BE IDENTICALLY EQUIVALENT TO THE MAXIMUM
!              NUMBER OF NONZERO ENTRIES IN A ROW OF THE FULL MATRIX.
!     MNORTH - MAXIMUM NUMBER OF ORTHOGONALIZATIONS PERFORMED
!              (AT LEAST MNORTH = 1 REQUIRED FOR CONJUGATE GRADIENT)
!
!
!     ORTHOMIN ARRAY SPACE/VARIABLES
!
!     NORTH  - NUMBER OF ORTHOGONALIZATIONS TO PERFORM
!            - SET NORTH=0 FOR SYMMETRIC MATRICES (CONJUGATE GRADIENT)
!     ECNVRG - RESIDUAL CONVERGENCE TOLERANCE
!     ACNVRG - ABSOLUTE CONVERGENCE TOLERANCE
!     RCNVRG - RELATIVE CONVERGENCE TOLERANCE
!     MAXIT  - MAXIMUM NUMBER OF ITERATIONS TO PERFORM
!     ITERP  - NUMBER OF ITERATIONS PERFORMED
!
!     B(MAXNB,MAXNN) - ILU DECOMPOSED MATRIX
!     Q(MAXNN)   - SEARCH DIRECTION Q
!     RQ(MAXNN)  - PRODUCT OF R AND Q
!     VRV(MAXNN) - EITHER V OR PRODUCT OF R AND V
!     RES(MAXNN) - RESIDUAL
!
!     QI(MAXNN,MNORTH)  - STORAGE OF Q´S
!     RQI(MAXNN,MNORTH) - STORAGE OF PRODUCTS OF R AND Q
!     RQIDOT(MNORTH)    - STORAGE OF DOT PRODUCTS OF RQ AND RQ
!
!     RESV   - PREVIOUS VALUE OF RES V DOT PRODUCT (CONJUGATE GRADIENT)
!
!     IAD(MAXNB,MAXNN) - ADJACENCY MATRIX (NODAL CONNECTIONS)
!
!     IADN(MAXNN) - NUMBER OF ADJACENT NODES IN IAD (SELF-INCLUSIVE)
!
!     IADD(MAXNN) - POSITION OF DIAGONAL IN ADJACENCY MATRIX
!
!
!     OTHER PARAMETERS PASSED FROM SOURCE PROGRAM
!
!     NN  - NUMBER OF NODES
!     NE  - NUMBER OF ELEMENTS
!     NLN - NUMBER OF LOCAL NODES IN AN ELEMENT
!
!
!     APPROXIMATE REAL STORAGE SPACE FOR ORTHOMIN AND MATRIX EQUATION
!
!       ((6 + 2!MAXNB + 2!MNORTH)!MAXNN)!(8 BYTES)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE IADMake(NumNP,NumEl,MaxNB,IAD,IADN,IADD)
      USE typedef
      USE GridData, ONLY: elmnod,iL,subN,nPt
      USE Matdata, ONLY: numnz,irow,jcol,A_sparse,alu,jlu
!     Generate the adjacency matrix for nodes from the element
!     indidence matrix
!     Requires subroutine Insert
      IMPLICIT NONE
      !dimension IAD(MaxNB,NumNP),IADN(NumNP),IADD(NumNP),iL(4)
      !integer e
      INTEGER(ap) :: NumNP,NumEl,MaxNB,i,j,k,l,kk,iSE
      INTEGER(ap) :: IAD(MaxNB,NumNP),IADN(NumNP),IADD(NumNP),iE
!     Determine independent adjacency within each element
!     version for SWMS_3D
        IADN=0
        IADD=0
        IAD=0
      DO 14 iE=1,NumEl
!       Loop on subelements
        DO 13 iSE=1,5
           i=elmnod(iL(1,iSE,subN(iE)),iE)!i,j,k,l are corner nodes of the tetraedal element, ise counts 5 which is equal to five tedrahedals which is equal to a cubic.
           j=elmnod(iL(2,iSE,subN(iE)),iE)
           k=elmnod(iL(3,iSE,subN(iE)),iE)
           l=elmnod(iL(4,iSE,subN(iE)),iE)
           CALL Insert(i,j,kk,NumNP,MaxNB,IAD,IADN)
           CALL Insert(j,i,kk,NumNP,MaxNB,IAD,IADN)
           CALL Insert(i,k,kk,NumNP,MaxNB,IAD,IADN)
           CALL Insert(k,i,kk,NumNP,MaxNB,IAD,IADN)
           CALL Insert(j,k,kk,NumNP,MaxNB,IAD,IADN)
           CALL Insert(k,j,kk,NumNP,MaxNB,IAD,IADN)
           CALL Insert(i,l,kk,NumNP,MaxNB,IAD,IADN)
           CALL Insert(l,i,kk,NumNP,MaxNB,IAD,IADN)
           CALL Insert(j,l,kk,NumNP,MaxNB,IAD,IADN)
           CALL Insert(l,j,kk,NumNP,MaxNB,IAD,IADN)
           CALL Insert(k,l,kk,NumNP,MaxNB,IAD,IADN)
           CALL Insert(l,k,kk,NumNP,MaxNB,IAD,IADN)
13        CONTINUE
14    CONTINUE
!number of nonzero entries
      numnz=SUM(IADN)+nPt !adjacencies + diagonal values
     
      ALLOCATE(jcol(numnz))
      ALLOCATE(a_sparse(numnz))
      ALLOCATE(irow(nPt+1))
    
      JCOL=0
      IROW=0
      !first value is a zero in IROW
      IROW(1)=1
      DO 15 i=1,NumNP
        CALL Insert(i,i,kk,NumNP,MaxNB,IAD,IADN)
        IADD(i)=kk
        IROW(1+i)=IROW(i)+IADN(i) !previous  irow value (at current i) + #entries in that row
        DO j=1,IADN(i)
           JCOL(IROW(i)-1+j)=IAD(j,i) ! column value added to JCOL
        ENDDO
15    CONTINUE       
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE INSERT (I,J,K,NN,MAXNB,IAD,IADN)
!     ADD J TO THE ADJACENCY LIST FOR I
!     RETURNS THE POSITION K WHERE IT HAS BEEN ADDED, OR WHERE IT
!     WAS ALREADY IN THE LIST.
      USE typedef
      IMPLICIT NONE
      INTEGER(ap) :: IAD(MAXNB,NN),IADN(NN)
      INTEGER(ap) :: I,J,K,NN,MAXNB,N,L,INODE
      !DIMENSION IAD(MAXNB,NN),IADN(NN)
      LOGICAL :: FOUND
      FOUND = .FALSE.
!     DETERMINE NUMBER OF NODES ALREADY IN ADJACENCY LIST
      N = IADN(I)
      K = N + 1
!     DETERMINE WHETHER ALREADY IN LIST
      DO 10 L=1,N
        INODE = IAD(L,I)
        IF (INODE.GE.J) THEN
          K = L
          IF (INODE.EQ.J) FOUND = .TRUE.
          GO TO 15
        ENDIF
   10 CONTINUE
   15 CONTINUE
!     PLACE IN LIST (NUMERICAL ORDER)
      IF (FOUND) THEN
       CONTINUE
      ELSE
        IF ((N+1).GT.MAXNB) THEN
          WRITE (* ,601) I,MAXNB
          WRITE (50,601) I,MAXNB
  601     FORMAT (//5X,'ERROR IN IADMAKE: NODE ',I5,' HAS > ' &
                 ,I5,' ADJACENCIES')
          STOP
        ENDIF
        IADN(I) = N + 1
        DO 20 L=(N+1),(K+1),(-1)
          IAD(L,I) = IAD(L-1,I)
   20   CONTINUE
       IAD(K,I) = J
      ENDIF
!write(*,*)'door insert'
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FIND (I,J,K,NN,MAXNB,IAD,IADN)
!     FOR NODE I, DETERMINE THE 'BAND' (K) RELATED TO ITS ADJACENCY TO
!     NODE J.
!     IF NODE NOT ADJACENT, RETURN 0 AS THE 'BAND'
      USE typedef
!      DIMENSION IAD(MAXNB,NN),IADN(NN)
     INTEGER(ap) :: NN,MAXNB,IAD(MAXNB,NN),IADN(NN),I,J,K,N,INODE,L
      K = 0
      N = IADN(I)
      DO 10 L=1,N
        INODE = IAD(L,I)
!     EXIT THE LOOP IF AT OR PAST THE REQUIRED POSITION
          IF (INODE.EQ.J) K = L
        IF (INODE.GE.J) THEN  
          GO TO 20
        ENDIF
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$      SUBROUTINE FIND2 (I,J,K,NN,MAXNB,IAD,N,DD)
!!$!     FOR NODE I, DETERMINE THE 'BAND' (K) RELATED TO ITS ADJACENCY TO
!!$!     NODE J.
!!$!     IF NODE NOT ADJACENT, RETURN 0 AS THE 'BAND'
!!$      USE typedef
!!$!      DIMENSION IAD(MAXNB,NN),IADN(NN)
!!$     INTEGER(ap) :: NN,MAXNB,IAD(MAXNB,NN),J,K,N,INODE,L,DD,start,I
!!$      K = 0
!!$      !N = IADN(I)
!!$      !DD = IADD(I)
!!$      IF (J.LT.I) THEN 
!!$         start = 1
!!$      ELSE
!!$         start = DD
!!$      END IF
!!$
!!$      DO 10 L=start,N
!!$        INODE = IAD(L,I)
!!$!     EXIT THE LOOP IF AT OR PAST THE REQUIRED POSITION
!!$          IF (INODE.EQ.J) K = L
!!$        IF (INODE.GE.J) THEN  
!!$          GO TO 20
!!$        ENDIF
!!$   10 CONTINUE
!!$   20 CONTINUE
!!$      RETURN
!!$      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ILU (R,NN,MAXNB,IAD,IADN,IADD,B)
!     INCOMPLETE LOWER-UPPER DECOMPOSITION OF MATRIX R INTO B
!     ONE STEP OF GAUSSIAN ELIMINATION PERFORMED
!     DIAGONAL DOMINANCE IS ASSUMED - NO PIVOTING PERFORMED
!     REQUIRES FUNCTION DU
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION R(MAXNB,NN),IAD(MAXNB,NN),IADN(NN),IADD(NN),B(MAXNB,NN)
      USE typedef
      IMPLICIT NONE
      INTEGER(ap) :: MAXNB,NN,I,J,N,K,ICUR,INODE,L,IAD(MAXNB,NN),IADN(NN),IADD(NN)
      REAL(sp) :: SUM,R(MAXNB,NN),D,DU
      REAL(sp), INTENT(out) :: B(MAXNB,NN)

!     INITIALIZE B
    B(1:MAXNB,1:NN)=0.
!write(*,*)'nn=',NN
!write(*,*)'maxnb=',maxnb
!write(*,*)'sizeR=',size(R,1),size(R,2)



!     LOOP OVER NODES
      DO 20 I=1,NN
!     DETERMINE NUMBER OF BANDS/POSITION OF DIAGONAL IN THIS ROW
        N = IADN(I)
        K = IADD(I)
!write(*,*)'N=',N
!write(*,*)'K=',K
!     LOWER TRIANGULAR MATRIX
!write(*,*)'ok'
!write(*,*)'sizeR=',size(R,1),size(R,2)
        DO 30 J=1,(K-1)
!write(*,*)'j=',J
!write(*,*)'R --',R(1:11,i)
          SUM = R(J,I)
!write(*,*)'ok1'
          ICUR = IAD(J,I)
          DO 40 L=1,(J-1)
            INODE = IAD(L,I)
!write(*,*)'ok2'
            SUM = SUM - B(L,I)*DU(INODE,ICUR,NN,MAXNB,IAD,IADN,IADD,B)
!write(*,*)'ok3'
   40     CONTINUE
          B(J,I) = SUM
!write(*,*)'B(j,i)=',B(J,I)
   30   CONTINUE
!write(*,*)'lower'
!     DIAGONAL
        SUM = R(K,I)
        DO 50 L=1,(K-1)
          INODE = IAD(L,I)
          SUM = SUM - B(L,I)*DU(INODE,I,NN,MAXNB,IAD,IADN,IADD,B)
   50   CONTINUE
        D = 1.0D0/SUM
        B(K,I) = D
!write(*,*)'diagonal'
!     UPPER TRIANGULAR MATRIX
!       - ACTUALLY D*U TO OBTAIN UNIT DIAGONAL
        DO 60 J=(K+1),N
          SUM = R(J,I)
          ICUR = IAD(J,I)
          DO 70 L=1,(K-1)
            INODE = IAD(L,I)
            SUM = SUM - B(L,I)*DU(INODE,ICUR,NN,MAXNB,IAD,IADN,IADD,B)
   70     CONTINUE
          B(J,I) = D*SUM
   60   CONTINUE
!write(*,*)'upper'
!write(*,*)'I,NN',I,NN
   20 CONTINUE
!write(*,*)'blaat'
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION DU (I,INODE,NN,MAXNB,IAD,IADN,IADD,B)
!     SEARCHES THE I´TH ROW OF THE UPPER DIAGONAL MATRIX
!     FOR AN ADJACENCY TO THE NODE 'INODE'
!     RETURNS CORRESPONDING VALUE OF B (OR ZERO)
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION IAD(MAXNB,NN),IADN(NN),IADD(NN),B(MAXNB,NN)
      USE typedef
      IMPLICIT NONE
      INTEGER(ap) :: MAXNB,NN,I,J,N,K,INODE,IAD(MAXNB,NN),IADN(NN),IADD(NN)
      REAL(sp) :: B(MAXNB,NN),TEMP,DU
      TEMP = 0.0D0
      N = IADN(I)
      K = IADD(I)
      IF (I.EQ.INODE) THEN
        TEMP = 1.0D0
        GO TO 20
      ENDIF
      DO 10 J=(K+1),N
        IF (INODE.EQ.IAD(J,I)) THEN
          TEMP = B(J,I)
          GO TO 20
        ENDIF
   10 CONTINUE
   20 CONTINUE
      DU = TEMP
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      SUBROUTINE ORTHOMIN (R,C,GT,NN,MAXNB,MAXNN,IAD,IADN,IADD,B,VRV,&
!                          RES,RQI,RQ,Q,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,&
!                          NORTH,MNORTH,MAXIT)
       SUBROUTINE OrthoMin(R,C,GT,NN,MAXNB,MAXNN,IAD,IADN,IADD,B,NORTH)
!      ORTHOMIN OR CONJUGATE GRADIENT ACCELERATION/SOLUTION
!      CONJUGATE GRADIENT (SYMMETRIC MATRIX) IF NORTH=0
!      (HOWEVER, NOTE THAT MNORTH MUST BE AT LEAST 1)
!      REQUIRES FUNCTIONS SDOT,SDOTK,SNRM
!      REQUIRES SUBROUTINES LUSOLV,MATM2,SAXPYK,SCOPY,SCOPYK
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION R(MAXNB,NN),C(NN),GT(NN),IAD(MAXNB,NN),IADN(NN)
!      DIMENSION IADD(NN),B(MAXNB,NN),VRV(NN),RES(NN),RQI(MAXNN,MNORTH)
!      DIMENSION RQ(NN),Q(NN),RQIDOT(MNORTH),QI(MAXNN,MNORTH)
       USE typedef
       IMPLICIT NONE
       INTEGER(ap) ::  MNorth=4
       INTEGER(ap), INTENT(in):: IAD(MAXNB,NN),IADN(NN),IADD(NN),MAXNB,MAXNN,NN
       REAL(sp),INTENT(in) :: B(MAXNB,NN),GT(NN),R(MAXNB,NN)
       REAL(sp), INTENT(out) :: C(NN)
       REAL(sp) :: VRV(NN),RES(NN),RQI(MAXNN,4),RQ(NN)
       REAL(sp) :: Q(NN),QI(MAXNN,4),RQIDOT(4)
       REAL(dp) :: ECNVRG=1e-6,RCNVRG=1e-6,ACNVRG=1e-6
       INTEGER(ap) :: MaxIt=1000,North
       INTEGER(ap) :: iHomog,I,norcur,iter,K,ITERp
       REAL(sp) :: DOT,ALPHA,OMEGA,RQNORM,RESMAX,DXNORM,CNORM,XRATIO,SDOT,SDOTK,SNRM2,RESV
!     INITIALIZE RESIDUAL VECTOR
      CALL MATM2 (RES,R,C,NN,IAD,IADN,MAXNB)
!write(*,*)'res',res
!pause
!     Solution for homogeneous system of equations - Modified by Simunek
      iHomog=0
      DO 10 I=1,NN
        RES(I) = GT(I) - RES(I)
        IF(ABS(GT(i)).GT.1.d-300) iHomog=1
   10 CONTINUE
      IF(iHomog.EQ.0) THEN
        DO 11 i=1,NN
	  C(i)=GT(i)
11      CONTINUE
        RETURN
      END IF
!write(*,*)'C',C
!pause
!     LOOP OVER A MAXIMUM OF MAXIT ITERATIONS
      NORCUR = 0
      DO 100 ITER=1,MAXIT
!     INVERT LOWER/UPPER MATRICES
        CALL SCOPY (NN,RES,VRV)
        CALL LUSOLV (NN,MAXNB,IAD,IADN,IADD,B,VRV)
!write(*,*)'VRV id not okay then ILU wrong',VRV
!pause
!     COPY V INTO Q
        CALL SCOPY (NN,VRV,Q)
!     CALCULATE PRODUCT OF R AND V
        CALL MATM2 (VRV,R,Q,NN,IAD,IADN,MAXNB)
!     COPY RV INTO RQ
        CALL SCOPY (NN,VRV,RQ)
!     RES V DOT PRODUCT (CONJUGATE GRADIENT)
        IF (NORTH.EQ.0) THEN
          DOT = SDOT(NN,RES,Q)
          IF (NORCUR.EQ.0) RESV = DOT
        ENDIF
!write(*,*)'DOT=',DOT
!pause
!     LOOP OVER PREVIOUS ORTHOGONALIZATIONS
        K = 1
   20   IF (K.GT.NORCUR) GO TO 30
!     DETERMINE WEIGHTING FACTOR (CONJUGATE GRADIENT)
          IF (NORTH.EQ.0) THEN
            ALPHA = DOT/RESV
            RESV = DOT
!write(*,*)'Alpha,RESV',ALPHA,RESV
!pause
!     DETERMINE WEIGHTING FACTOR (ORTHOMIN)
          ELSE
            DOT = SDOTK(NN,K,RQI,VRV,MAXNN,MNORTH)
            ALPHA = -DOT/RQIDOT(K)
          ENDIF
!     SUM TO OBTAIN NEW Q AND RQ
          CALL SAXPYK (NN,ALPHA,K,QI,Q,MAXNN,MNORTH)
          CALL SAXPYK (NN,ALPHA,K,RQI,RQ,MAXNN,MNORTH)
!write(*,*)'QI=',QI(1:NN,1)
!pause
!write(*,*)'RQI=',RQI(1:NN,1)
!pause
          K = K + 1
          GO TO 20
   30   CONTINUE
!     CALCULATE WEIGHTING FACTOR (CONJUGATE GRADIENT)
        IF (NORTH.EQ.0) THEN
          DOT = SDOT(NN,Q,RQ)
          OMEGA = RESV/DOT
!write(*,*)'omega=',OMEGA
!     CALCULATE WEIGHTING FACTOR (ORTHOMIN)
        ELSE
          DOT = SDOT(NN,RES,RQ)
          RQNORM = SDOT(NN,RQ,RQ)
          OMEGA = DOT/RQNORM
        ENDIF
!     SAVE VALUES FOR FUTURE ORTHOGONALIZATIONS
!write(*,*)'Q=',Q
!write(*,*)'RQ=',RQ
!pause
        NORCUR = NORCUR + 1
        IF (NORCUR.GT.NORTH) NORCUR = 1
        CALL SCOPYK (NN,NORCUR,Q,QI,MAXNN,MNORTH)
        CALL SCOPYK (NN,NORCUR,RQ,RQI,MAXNN,MNORTH)
        RQIDOT(NORCUR) = RQNORM
! write(*,*)'QI2=',QI(1:NN,1)
!pause
!write(*,*)'RQI2=',RQI(1:NN,1)
!pause
!     UPDATE SOLUTION/RESIDUAL VECTORS
        CALL SAXPYK (NN,OMEGA,NORCUR,QI,C,MAXNN,MNORTH)
        CALL SAXPYK (NN,-OMEGA,NORCUR,RQI,RES,MAXNN,MNORTH)
!     DETERMINE CONVERGENCE PARAMETERS
        RESMAX = SNRM2(NN,RES)
        DXNORM = ABS(OMEGA)*SNRM2(NN,Q)
        CNORM = SNRM2(NN,C)
        XRATIO = DXNORM/CNORM
!write(*,*)'resmax',resmax
!write(*,*)'DXNORM',dxnorm
!write(*,*)'cnorm',cnorm
!write(*,*)'xratio',xratio
!     ITERATION (DEBUG) OUTPUT
!     STOP ITERATING IF CONVERGED
        IF (RESMAX.LT.ECNVRG) GO TO 200
        IF (XRATIO.LT.RCNVRG) GO TO 200
        IF (DXNORM.LT.ACNVRG) GO TO 200
  100 CONTINUE
!     TERMINATE IF TOO MANY ITERATIONS
      WRITE (* ,602) MAXIT,RESMAX,XRATIO,DXNORM
      WRITE (70,602) MAXIT,RESMAX,XRATIO,DXNORM
 602  FORMAT (///5X,'ORTHOMIN TERMINATES -- TOO MANY ITERATIONS', &
               /8X,'MAXIT  = ',I5, &
               /8X,'RESMAX = ',E12.4, &
               /8X,'XRATIO = ',E12.4, &
               /8X,'DXNORM = ',E12.4)
      STOP
  200 CONTINUE
!     RETURN NUMBER OF ITERATIONS REQUIRED
      ITERP = ITER
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE LUSOLV (NN,MAXNB,IAD,IADN,IADD,B,VRV)
!     LOWER DIAGONAL MATRIX INVERSION BY FORWARD SUBSTITUTION
!     UPPER DIAGONAL MATRIX INVERSION BY BACKWARD SUBSTITUTION
!     LOWER/UPPER MATRICES ARE IN B
!     RIGHT-HAND-SIDE VECTOR IS IN VRV AT START
!     'SOLUTION' IS RETURNED IN VRV UPON EXIT
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION IAD(MAXNB,NN),IADN(NN),IADD(NN),B(MAXNB,NN),VRV(NN)
        USE typedef
        IMPLICIT NONE
       INTEGER(ap):: IAD(MAXNB,NN),IADN(NN),IADD(NN),MAXNB,NN
       REAL(sp) :: B(MAXNB,NN),VRV(NN),SUM
       INTEGER(ap) :: I,K,J,INODE,N
!     LOWER INVERSION
      DO 20 I=1,NN
        SUM = VRV(I)
        K = IADD(I)
        DO 30 J=1,(K-1)
          INODE = IAD(J,I)
          SUM = SUM - B(J,I)*VRV(INODE)
   30   CONTINUE
        VRV(I) = B(K,I)*SUM
   20 CONTINUE
!     UPPER INVERSION
      DO 40 I=NN,1,-1
        SUM = VRV(I)
        N = IADN(I)
        K = IADD(I)
        DO 50 J=(K+1),N
          INODE = IAD(J,I)
          SUM = SUM - B(J,I)*VRV(INODE)
   50   CONTINUE
        VRV(I) = SUM
   40 CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MATM2 (S1,R,P,NN,IAD,IADN,MAXNB)
!     MULTIPLY MATRIX R BY VECTOR P TO OBTAIN S1
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION S1(NN),P(NN),R(MAXNB,NN),IAD(MAXNB,NN),IADN(NN)
        USE typedef
        IMPLICIT NONE
       INTEGER(ap):: IAD(MAXNB,NN),IADN(NN),MAXNB,NN
       REAL(sp) :: P(NN),SUM,R(MAXNB,NN),S1(NN)
       INTEGER(ap) :: I,J,INODE,N
      DO 30 I=1,NN
        SUM = 0.0D0
        N = IADN(I)
        DO 40 J=1,N
          INODE = IAD(J,I)
          SUM = SUM + R(J,I)*P(INODE)
   40   CONTINUE
        S1(I) = SUM
   30 CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SDOT (NN,R,B)
!     OBTAIN DOT PRODUCT OF R AND B
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION R(NN),B(NN)
        USE typedef
        IMPLICIT NONE
        REAL(sp) :: SDOT,R(NN),B(NN)
        INTEGER(ap) :: NN,L
      SDOT = 0.0D0
      DO 100 L=1,NN
        SDOT = SDOT + R(L)*B(L)
100   CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SDOTK (NN,K,R,B,MAXNN,MNORTH)
!     OBTAIN DOT PRODUCT OF R AND B
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION R(MAXNN,MNORTH),B(NN)
        USE typedef
        IMPLICIT NONE
        REAL(sp) :: SDOTK,B(NN),R(MAXNN,MNorth)
        INTEGER(ap) :: NN,L,K,MNORTH,MAXNN
      SDOTK = 0.0D0
      DO 100 L=1,NN
        SDOTK = SDOTK + R(L,K)*B(L)
100   CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION SNRM2 (NN,R)
!     COMPUTE MAXIMUM NORM OF R
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION R(NN)
      USE typedef
      INTEGER(ap) :: L,NN
      REAL(sp) :: R(NN),TEMP,SNRM2
      SNRM2 = 0.0D0
      DO 100 L=1,NN
        TEMP = ABS(R(L))
        IF (TEMP.GT.SNRM2) SNRM2 = TEMP
100   CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SAXPYK (NN,SA,K,FX,FY,MAXNN,MNORTH)
!     MULTIPLY VECTOR FX BY SCALAR SA AND ADD TO VECTOR FY
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION FX(MAXNN,MNORTH),FY(NN)
      USE typedef
      INTEGER(ap) :: I,NN,K,MAXNN,MNORTH
      REAL(sp):: FX(MAXNN,MNORTH),FY(NN),SA
      IF (NN.GT.0) THEN
        DO 100 I=1,NN
           FY(I) = SA*FX(I,K) + FY(I)
100     CONTINUE
      ENDIF
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SCOPY (NN,FX,FY)
!     COPY VECTOR FX INTO VECTOR FY
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION FX(NN),FY(NN)
        USE typedef
        INTEGER(ap) :: I,NN
        REAL(sp):: FX(NN),FY(NN)
      IF (NN.GT.0) THEN
        DO 100 I=1,NN
           FY(I) = FX(I)
100     CONTINUE
      ENDIF
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SCOPYK (NN,K,FX,FY,MAXNN,MNORTH)
!     COPY VECTOR FX INTO VECTOR FY
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION FX(N),FY(MAXNN,MNORTH)
        USE typedef
      INTEGER(ap) :: I,NN,K,MAXNN,MNORTH
      REAL(sp):: FY(MAXNN,MNORTH),FX(NN)
      IF (NN.GT.0) THEN
        DO 100 I=1,NN
           FY(I,K) = FX(I)
100     CONTINUE
      ENDIF
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

