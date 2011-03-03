***********************************************************************************
**  UMAT, FOR ABAQUS/STANDARD INCORPORATING ELASTIC-PLASTIC LINEAR               **
**  ISOTROPIC HARDENING. LARGE DEFORMATION FORMULATION FOR PLANE STRAIN          **
**  AND AXI-SYMMETRIC ELEMENTS. IMPLICIT INTEGRATION WITH CONSISTENT JACOBIAN    **
***********************************************************************************
***********************************************************************************
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
     +          SIX=6.D0, NINE=9.D0, TOLER=1.D-5)
C
      DIMENSION XIDEN(M,N),XNV(4),DPSTRAN(4), STRESSOLD(4),
     +          DESTRAN(4), DSTRESS(4), XNDIR(M,N),
     +          STR(M,N),DSTR(M,N), DPSTRN(M,N)          
C
C
C
C
C
C--------------------------------------------------------------------
C
C     SPECIFY MATERIAL PROPERTIES
C
C
      E = 210000.0
      XNUE = 0.3
      SIGY0 = 240.
      h = 1206.0
C
C    ZERO MATRICES WITH UNSPECIFIED COMPONENTS
C
      DO I=1,M
      DO J=1,N
      STR(I,J)=0.0
      END DO
      END DO
C
C    RECOVER EFFECTIVE PLASTIC STRAIN, p, AND ISOTROPIC 
C    HARDENING VARIABLE, r,FROM PREVIOUS TIME STEP
C
      p = STATEV(1)
      r = STATEV(2)
C
C 
C    SET UP ELASTICITY MATRIX
C   
      EBULK3 = E/(ONE-TWO*XNUE)
      XK = EBULK3/THREE
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
C
C     DEFINE IDENTITY MATRIX
C
      DO 50 I=1,M
      DO 50 J=1,N
        IF(I .EQ. J) THEN
          XIDEN(I,J)=1.0D0
        ELSE
          XIDEN(I,J)=0.0D0
        END IF
  50  CONTINUE
C
C     SAVE STRESS AT BEGINNING OF TIME STEP IN STRESSOLD
C
      DO K=1,4
      STRESSOLD(K) = STRESS(K)
      END DO
C
C     OBTAIN TRIAL (ELASTIC) STRESS
C
      CALL KMLT1(DDSDDE,DSTRAN,DSTRESS,NTENS)
      DO K=1,NTENS
      STRESS(K) = STRESS(K) + DSTRESS(K)
      END DO
C
C     WRITE TRIAL STRESSES IN MATRIX FORM
C      
      DO K = 1,3
      STR(K,K) = STRESS(K)
      END DO
      STR(1,2) = STRESS(4)
      STR(2,1) = STRESS(4)
C
C     CALCULATE DEVIATORIC TRIAL STRESS
C 
      CALL KDEVIA(STR,XIDEN,DSTR) 
C
C     CALCULATE EFFECTIVE TRIAL STRESS
C
      CALL KEFFP(DSTR,PJ)
C
C    DETERMINE TRIAL FLOW DIRECTION FOR USE WITH THE JACOBIAN
C
      DO I=1,3
      DO J=1,3
      XNDIR(I,J) = DSTR(I,J)/PJ
      END DO
      END DO
C
C    ...AND WRITE IN VOIGT NOTATION
C
      DO K=1,3
      XNV(K) = XNDIR(K,K)
      END DO
      XNV(4) = XNDIR(1,2)
C
C    DETERMINE IF THE YIELD CONDITION IS SATISFIED
C
      ZY = PJ - r - SIGY0
C
      IF (ZY.GT.0.) THEN
C
C     USE NEWTON ITERATION TO DETERMINE EFFECTIVE PLASTIC STRAIN INCREMENT
C 
      r0 = r
      dp = 0.
      DO KNEWT=1,10
      RES = PJ-3.*EG*dp-r-SIGY0
      dp = dp + RES/(3.*EG+h)
      r = r0 + h*dp 
      IF(DABS(RES).LT.TOLER) GOTO 10
      END DO
   10 CONTINUE
C
C     DETERMINE THE INCREMENTS IN PLASTIC STRAIN
C
      DO I = 1,3
       DO J = 1,3
       DPSTRN(I,J) = (THREE/TWO)*dp*DSTR(I,J)/PJ
       END DO
      END DO 
C
C     WRITE THE STRAIN INCREMENTS IN VOIGT NOTATION (WITH ENGG SHEARS)
C
      DO K=1,3
      DPSTRAN(K) = DPSTRN(K,K)
      END DO
      DPSTRAN(4) = 2*DPSTRN(1,2)
C      
C     CALCULATE THE ELASTIC STRAIN INCREMENTS      
C      
      DO K=1,4
      DESTRAN(K)=DSTRAN(K)-DPSTRAN(K)
      END DO  
C
C     DETERMINE STRESS INCREMENT
C
       CALL KMLT1(DDSDDE,DESTRAN,DSTRESS,NTENS)
C
C      UPDATE THE STRESS, EFFECTIVE PLASTIC STRAIN
C      (NOTE: ISOTROPIC HARDENING VARIABLE ALREADY UPDATED)
C
       DO K = 1,NTENS
       STRESS(K) = STRESSOLD(K) + DSTRESS(K)
       END DO
C
       p = p + dp
C       
C     STORE UPDATED STATE VARIABLES
C
       STATEV(1) = p
       STATEV(2) = r
C
C
C    DETERMINE JACOBIAN
C
      XR = (PJ-3.*EG*dp)/PJ
      Q = (ONE/(ONE+3.*EG/h)-XR)*THREE/TWO
      DO I=1,3
      DO J=1,3
       DDSDDE(I,J) = 2*EG*Q*XNV(I)*XNV(J) 
     +  + (XK-EG*XR*TWO/THREE) + TWO*EG*XR*XIDEN(I,J)
      END DO
      END DO
      DO K=1,3
      DDSDDE(K,4) = TWO*EG*Q*XNV(K)*XNV(4)
      DDSDDE(4,K) = DDSDDE(K,4)
      END DO
      DDSDDE(4,4) = TWO*EG*Q*XNV(4)*XNV(4) + EG*XR
C
      END IF
C
      RETURN
      END
**
***********************************************
**           UTILITY    SUBROUTINES           *
***********************************************
**
**
***************************************************
**         MULTIPLY 4X4 MATRIX WITH 4X1 VECTOR    *
***************************************************
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
      X=X+Y
20    CONTINUE
      DM(I)=X
10    CONTINUE
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
**
********************************************
**         DOT PRODUCT OF TWO VECTORS      *
********************************************
*USER SUBROUTINE
      SUBROUTINE DOTPROD(DM1,DM2,DM,NTENS)
C      
      INCLUDE 'ABA_PARAM.INC'
C
C      PARAMETER (M=4)
C
      DIMENSION DM1(4),DM2(4)
C
      Y=0.0
      DO 20 K=1,NTENS 
      X=DM1(K)*DM2(K)
      Y=X+Y
   20 CONTINUE
      DM=Y
      RETURN
      END
**
      SUBROUTINE DYADICPROD(DM1,DM2,DM3,NTENS)
C      
      INCLUDE 'ABA_PARAM.INC'
C
C      PARAMETER (M=4)
C
      DIMENSION DM1(4),DM2(4),DM3(4,4)
C
      DO I=1,4
       DO J=1,4
        DM3(I,J) = DM1(I)*DM2(J)
       END DO
      END DO
C  
      RETURN
      END
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















