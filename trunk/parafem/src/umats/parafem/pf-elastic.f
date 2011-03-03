*****************************************************************************
**  UMAT FOR ABAQUS/STANDARD INCORPORATING ELASTIC BEHAVIOUR FOR ISOTROPIC **
**  ISOTHERMAL ELASTICITY                                                  ** 
**                                                                         **
**  CANNOT BE USED FOR PLANE STRESS                                        **
**                                                                         **             
**  PROPS(1) - E                                                           **
**  PROPS(2) - NU                                                          **
**                                                                         **
**  Original version modified by                                           **
**                                                                         **
**  louise.lever@manchester.ac.uk                                          **
**  lee.margetts@manchester.ac.uk                                          **
**                                                                         **
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
C    INCLUDE 'ABA_PARAM.INC'
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
C------------------------------------------------------------------------------
C     
C    TEST FOR NUMBER OF DIMENSIONS
C
      IF(NDI.NE.3) THEN
         WRITE(7,*) 'THIS UMAT MAY ONLY BE USED FOR ELEMENTS
     1   WITH THREE DIRECT STRESS COMPONENTS'
         RETURN
      END IF
C
C
C------------------------------------------------------------------------------
C 
C    ELASTIC PROPERTIES
C   
      EMOD   = PROPS(1)
      ENU    = PROPS(2)
      EBULK3 = EMOD/(ONE-TWO*ENU)
      EG2    = EMOD/(ONE+ENU)
      EG     = EG2/TWO
      ELAM   = (EBULK3-EG2)/THREE 
C
C
C------------------------------------------------------------------------------
C
C    ELASTIC STIFFNESS
C
C    REAL ARRAY DDSDDE() IS THE "DEE" MATRIX IN S&G EDITION 4
C
      DO K1 = 1, NDI
        DO K2 = 1, NDI
          DDSDDE(K2,K1) = ELAM
        END DO
        DDSDDE(K1,K1)   = EG2 + ELAM
      END DO
C
      DO K1 = NDI+1, NTENS
        DDSDDE(K1,K1) = EG
      END DO
C        
C
C------------------------------------------------------------------------------
C
C     CALCULATE STRESS
C
C     REAL ARRAY DSTRAN() IS "MATMUL(BEE,ELD)" IN S&G EDITION 4
C     THIS LOOP IS THE EQUIVALENT OF THE LINE:
C
C       SIGMA=MATMUL(DEE,MATMUL(BEE,ELD))
C
       DO K1=1,NTENS
         DO K2=1,NTENS
           STRESS(K2)=STRESS(K2)+DDSDDE(K2,K1)*DSTRAN(K1)
         END DO
       END DO
C
       RETURN
       END
**


