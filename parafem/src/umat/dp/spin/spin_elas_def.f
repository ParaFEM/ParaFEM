*****************************************************************************
**  UMAT, FOR ABAQUS/STANDARD FOR ELASTIC BEHAVIOUR USING THE DEFORMATION  **
**  GRADIENT. LARGE DEFORMATION FORMULATION FOR PLANE STRAIN,              **
**  AXI-SYMMETRIC AND 3D ELEMENTS                                          **
*****************************************************************************
*****************************************************************************
**
**
**
*USER SUBROUTINE
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
C
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C
      PARAMETER (M=3,N=3,ID=3,ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,
     +          SIX=6.D0, NINE=9.D0, TOLER=0.D-6)
C
      DIMENSION DFGI(M,N),DFGR(M,N),XIDEN(M,N),
     +          VEG(M,N),TVEG(M,N),DFRT(M,N),
     +          STR(M,N),STRR(M,N),DSTR(M,N),
     +          SPINW(M,N),WS(M,N),SW(M,N) 
C
C
C
C
C
C--------------------------------------------------------------------
C
      E = 210000.
      XNUE = 0.3
C
C    ZERO ARRAYS
C
      DO 20 I=1,M
      DO 20 J=1,N
        STR(I,J)=0.0
        STRR(I,J)=0.0
        DSTR(I,J)=0.0
        DFRT(I,J)=0.0
        DFGR(I,J)=0.0
        DFGI(I,J)=0.0
20    CONTINUE
C
      DO I=1,4
      DO J=1,4
      DDSDDE(I,J) = 0.0
      END DO
      END DO
      IF(NTENS.GT.4) THEN
      DO I=1,6
      DO J=1,6
      DDSDDE(I,J) = 0.0
      END DO
      END DO
      END IF  
C
C    WRITE STRESSES FROM PREVIOUS TIME STEP IN TO ARRAY STR
C   
      DO K=1,3
      STR(K,K) = STATEV(K)
      END DO
      STR(1,2) = STATEV(4)
      STR(2,1) = STATEV(4)
      IF(NTENS.GT.4)THEN
      STR(1,3) = STATEV(5)
      STR(3,1) = STATEV(5)
      STR(2,3) = STATEV(6)
      STR(3,2) = STATEV(6)
      END IF
C
C     DEFINE XIDENTITY MATRIX
C
      DO 50 I=1,M
      DO 50 J=1,N
        IF(I .EQ. J) THEN
          XIDEN(I,J)=1.0D0
        ELSE
          XIDEN(I,J)=0.0D0
        END IF
50    CONTINUE
C
C 
C    SET UP ELASTICITY MATRIX
C   
      EBULK3 = E/(ONE-TWO*XNUE)
      EG2 = E/(ONE+XNUE)
      EG = EG2/TWO
      ELAM = (EBULK3-EG2)/THREE 
C
C
      DO K1 = 1, 3
         DO K2 = 1, 3
         DDSDDE(K2,K1) = ELAM
         END DO
         DDSDDE(K1,K1) = EG2 + ELAM
      END DO
C
        DDSDDE(4,4) = EG
        IF(NTENS.GT.4) THEN
        DDSDDE(5,5) = EG
        DDSDDE(6,6) = EG
        END IF
C  
C    DETERMINE DEFORMATION GRADIENT RATE  
C
      IF(DTIME.GT.0.) THEN
      DO 110 I=1,M
      DO 110 J=1,N
        DFGR(I,J)=(DFGRD1(I,J)-DFGRD0(I,J))/DTIME
110   CONTINUE 
      END IF
C
C     DETERMINE VELOCITY GRADIENT     
C
      CALL KINVER(DFGRD0,DFGI)
C
      CALL KMLT(DFGR,DFGI,VEG)
C    
C    DETERMINE RATE OF DEFORMATION FROM VELOCITY GRADIENT  
C
      CALL KTRANS(VEG,TVEG)
C
      DO 120 I=1,M
      DO 120 J=1,N
      DFRT(I,J)=(ONE/TWO)*(VEG(I,J)+TVEG(I,J))
120   CONTINUE
C
C    DETERMINE JAUMANN STRESS RATE                         
C
      CALL KTRACE(DFRT,TRVAL)
C
      DO 160 I=1,M
      DO 160 J=1,N
        STRR(I,J)=((E*DFRT(I,J))/(ONE+XNUE))
     +  +((E*XNUE*TRVAL*XIDEN(I,J))/(ONE-XNUE-(TWO*XNUE*XNUE)))
160   CONTINUE
C
C    DETERMINE CONTINUUM SPIN
C
      DO 201 I=1,M
      DO 201 J=1,N
      SPINW(I,J)=(ONE/TWO)*(VEG(I,J)-TVEG(I,J))
  201 CONTINUE
C     
C    DETERMINE STRESS RATE WRT UNDEFORMED CONFIGURATION              
C
      CALL KMLT(SPINW,STR,WS)
      CALL KMLT(STR,SPINW,SW)
C
C
      DO 231 I=1,M
      DO 231 J=1,N
      STRR(I,J)=STRR(I,J)+WS(I,J)-SW(I,J)
C
C     .... AND  INTEGRATE         
C
      STR(I,J) = STR(I,J) + STRR(I,J)*DTIME
 231  CONTINUE
C
C    WRITE UPDATED STRESSES IN TO ABAQUS ARRAY STRESS
C
      DO I=1,3
      STRESS(I)=STR(I,I)
      END DO
      STRESS(4)=STR(1,2)
      IF(NTENS.GT.4) THEN
      STRESS(5)=STR(1,3)
      STRESS(6)=STR(2,3)
      END IF
C
C     STORE UPDATED STRESSES IN STATE VARIABLES
C
      DO K=1,3
      STATEV(K) = STR(K,K)
      END DO
      STATEV(4)=STR(1,2)
      IF(NTENS.GT.4) THEN
      STATEV(5)=STR(1,3)
      STATEV(6)=STR(2,3)
      END IF
C
C
C***************************************************
C
C
      RETURN
      END
**
***********************************************
**           UTILITY    SUBROUTINES           *
***********************************************
**************************************
**           INVERSE MATRIX          *
**************************************
*USER SUBROUTINE
      SUBROUTINE KINVER(DF,A)
C      
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (M=3,N=3)
      DIMENSION DF(M,N),A(M,N)
C
      DO 5 I=1,M
      DO 5 J=1,N
      A(I,J)=DF(I,J)
5     CONTINUE
      DO 10 K=1,M
      P=A(K,K)
      A(K,K)=1.
      DO 20 J=1,N
      A(K,J)=A(K,J)/P
20    CONTINUE
      DO 10 I=1,M
      IF(I .EQ. K) GO TO 10
      P=A(I,K)
      A(I,K)=0.
      DO 30 J=1,N
      A(I,J)=A(I,J)-A(K,J)*P
30    CONTINUE
10    CONTINUE
      RETURN
      END                                                  
**
**
**************************************
**         MULTIPLY MATRIX  1        *
**************************************
*USER SUBROUTINE
      SUBROUTINE KMLT(DM1,DM2,DM)
C      
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (M=3,N=3)
      DIMENSION DM1(M,N),DM2(M,N),DM(M,N)
C
      DO 10 I=1,M
      DO 10 J=1,N
      X=0.0
      DO 20 K=1,M 
      X=X+DM1(I,K)*DM2(K,J)
20    CONTINUE
      DM(I,J)=X
10    CONTINUE
      RETURN
      END
**
**
*************************************************************
**         MULTIPLY MATRIX   AXISYMM AND PLANE STRAIN       *
*************************************************************
*USER SUBROUTINE
      SUBROUTINE KMLT1(DM1,DM2,DM,NTENS)
C      
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (M=4)
C
      DIMENSION DM1(M,M),DM2(M),DM(M)
C
      DO 10 I=1,NTENS
      X=0.0
      DO 20 K=1,NTENS 
      Y=DM1(I,K)*DM2(K)
      IF (K.GE.4) Y = 2*Y
      X=X+Y
20    CONTINUE
      DM(I)=X
10    CONTINUE
      RETURN
      END
**
**
**************************************
**         MULTIPLY MATRIX   3D      *
**************************************
*USER SUBROUTINE
      SUBROUTINE KMLT2(DM1,DM2,DM,NTENS)
C      
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (M=6)
C
      DIMENSION DM1(M,M),DM2(M),DM(M)
C
      DO 10 I=1,NTENS
      X=0.0
      DO 20 K=1,NTENS 
      Y=DM1(I,K)*DM2(K)
      IF (K.GE.4) Y = 2*Y
      X=X+Y
20    CONTINUE
      DM(I)=X
10    CONTINUE
C      DO I=1,NTENS
C      DO J=1,NTENS
C      WRITE(*,*) DM1(I,J)
C      END DO
C      END DO
C      STOP
C      WRITE(*,*) DM2(1),DM2(2),DM2(3),DM2(4),DM2(5),DM2(6)
C      WRITE(*,*) DM(1),DM(2),DM(3),DM(4),DM(5),DM(6)
      RETURN
      END
**
**
***************************************
**          EFFECTIVE STRESS          *
**   (CONTRACTED MATRIX CALCULATION)  *
***************************************
*USER SUBROUTINE
      SUBROUTINE KEFFP(EFF1,VAL1)
C      
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (M=3,N=3)
      DIMENSION EFF1(M,N)
C
      X=0.0
      DO 10 I=1,M
      DO 10 J=1,N
       X=X+EFF1(I,J)*EFF1(I,J)
10    CONTINUE
      IF(X .LE. 0.0) GO TO 20 
      VAL1=DSQRT((3.0/2.0)*X)
20    RETURN
      END
**
**
***************************************
**   EFFECTIVE STRAIN (RATE)          *
**  (CONTRACTED MATRIX CALCULATION)   *
***************************************
*USER SUBROUTINE
      SUBROUTINE KEFFPS(EFF2,VAL2)
C      
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (M=3,N=3)
      DIMENSION EFF2(M,N)
C
      X=0.0
      DO 10 I=1,M
      DO 10 J=1,N
       X=X+EFF2(I,J)*EFF2(I,J)
10    CONTINUE
      IF(X .LE. 0.0) GO TO 20 
      VAL2=DSQRT((2.0/3.0)*X)
20    RETURN
      END
**
**
***************************************
**   TRACE OF MATRIX CALCULATION      *
***************************************
*USER SUBROUTINE
      SUBROUTINE KTRACE(DE,TVAL)
C      
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (M=3,N=3)
      DIMENSION DE(M,N)
C
      X=0.0
      DO 10 I=1,M
      DO 10 J=1,N
      IF(I .EQ. J) THEN
      X=X+DE(I,J)
      ELSE
      END IF
10    CONTINUE
      TVAL=X
      RETURN
      END
**
**
*****************************************************
**   DEVIATORIC STRESS CALCULATION    *
*****************************************************
*USER SUBROUTINE
      SUBROUTINE KDEVIA(STRSS,XIDENTY,DEVITO)
C      
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (M=3,N=3)
      DIMENSION STRSS(M,N),XIDENTY(M,N),DEVITO(M,N)
C
      X=0.0
      DO 10 I=1,M
      DO 10 J=1,N
      IF(I .EQ. J) THEN
      X=X+STRSS(I,J)
      ELSE
      END IF
10    CONTINUE
C
      DO 20 I=1,M
      DO 20 J=1,N
      IF(I .EQ. J) THEN
        DEVITO(I,J)=STRSS(I,J)-((1./3.)*X*XIDENTY(I,J))
      ELSE
        DEVITO(I,J)=STRSS(I,J)
      END IF
20    CONTINUE
      RETURN
      END
**
**
******************************************************
** TRANSPOSE OF MATRIX                               *
******************************************************
*USER SUBROUTINE
      SUBROUTINE KTRANS(ORIGN,TRAN)
C      
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (M=3,N=3)
      DIMENSION ORIGN(M,N),TRAN(M,N)
C
      DO 10 I=1,M
      DO 10 J=1,N
      TRAN(J,I)=ORIGN(I,J)
10    CONTINUE
      RETURN
      END
**
**
***********************************************
** CONTRACTED TENSOR PRODUCT CALCULATION      *
***********************************************
*USER SUBROUTINE
      SUBROUTINE KTENPROD(EFF1,EFF2,VAL1)
C      
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER (M=3,N=3)
      DIMENSION EFF1(M,N), EFF2(M,N)
C
      X=0.0
      DO 10 I=1,M
      DO 10 J=1,N
       X=X+EFF1(I,J)*EFF2(I,J)
10    CONTINUE
      IF(X .LE. 0.0) GO TO 20 
      VAL1=X
20    RETURN
      END
**
**
**
*USER SUBROUTINE
      SUBROUTINE DISP(U,KSTEP,KINC,TIME,NODE,JDOF)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION U(3), TIME(2)
C
      STRRAT = 0.2
C
C      U(1)=STRRAT*TIME(1)
      U(1) = DEXP(STRRAT*TIME(1)) - 1.
C
C
      RETURN
C
      END
**
