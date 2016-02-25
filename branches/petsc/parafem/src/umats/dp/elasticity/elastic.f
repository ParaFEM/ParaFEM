*****************************************************************************
**  UMAT FOR ABAQUS/STANDARD INCORPORATING ELASTIC BEHAVIOUR  FOR PLANE    **
**  STRAIN AND AXI-SYMMETRIC ELEMENTS.                                     **
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
      DIMENSION DSTRESS(4), DDS(4,4)
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
           DDS(K2,K1) = ELAM
         END DO
        DDS(K1,K1) = EG2 + ELAM
      END DO
C
        DDS(4,4) = EG
C
C
C     DETERMINE STRESS INCREMENT
C
C
       TRVAL = DSTRAN(1)+DSTRAN(2)+DSTRAN(3)
       DO K=1,3
       DSTRESS(K) = 2*EG*DSTRAN(K)+ELAM*TRVAL
       END DO
       DSTRESS(4) = EG*DSTRAN(4)
C
C      UPDATE STRESS
C
       DO K = 1,NTENS
       STRESS(K) = STRESS(K) + DSTRESS(K)
       END DO
C
C
C    DETERMINE JACOBIAN
C
C
      DO I=1,3
       DO J=1,3
        DDSDDE(I,J) = DDS(I,J)
       END DO
      END DO
      DDSDDE(4,4) = DDS(4,4)
C
C
      RETURN
      END
**


