MODULE PLASTICITY

  !------------------------------------------------------------------------------
  !/****h* ParaFEM-Shared/plasticity
  !*  NAME
  !*    plasticity module
  !*  FUNCTION
  !*    Contains the following subroutines and functions needed to solve 
  !*    plasticity problems:
  !*
  !*    Subroutine             Purpose
  !*
  !*    INVAR                  Returns stress invariants SIGMA, DSBAR and 
  !*                           Lode angle THETA from current stresses held in 
  !*                           STRESS.
  !*
  !*    VMFLOW                 Returns the von Mises flow vector VMFL. STRESS 
  !*                           holds the second deviatoric invariant (2D only).
  !*
  !*    FMACAT                 Intermediate step.
  !*
  !*    FMRMAT                 Returns matrix RMAT from the von Mises flow 
  !*                           vector VMFL, invariant DSBAR, plastic multiplier
  !*                           DLAM and elastic matrix DEE.
  !*
  !*    FORMAA                 Returns modified matrix DAATD from the von Mises
  !*                           flow vector VMFL and matrix RMAT.
  !*
  !*    VMPL                   Returns the plastic stress-strain matrix PL for
  !*                           a von Mises material. STRESS holds the stresses
  !*                           and DEE holds the elastic stress-strain matrix.
  !*
  !*    FORM_TEMP              Forms temp matrix for consistent return.
  !*
  !*    MOCOUF                 Returns the Mohr-Coulomb failure function, F, 
  !*                           from the strength parameters PHI and C and 
  !*                           stress invariants SIGM, DSBAR and THETA.
  !*
  !*    MOCOUQ                 Returns the plastic potentail terms DQ1, DQ2
  !*                           and DQ3 for a Mohr-Coulomb material from
  !*                           dilation angle PSI (in degrees) and invariants
  !*                           DSBAR and THETA.
  !*
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !****
  !*/

  USE precision
  
  CONTAINS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE INVAR(stress,sigm,dsbar,theta)
  
  !/****f* plasticity/invar
  !*  NAME
  !*    SUBROUTINE: invar
  !*  SYNOPSIS
  !*    Usage:      CALL invar(stress,sigm,dsbar,theta)
  !*  FUNCTION
  !*    Forms the stress invariants in 2-d or 3-d.
  !*  INPUTS
  !*    The following real array argument has the INTENT(IN) attribute:
  !*
  !*    stress(:)              : stress vector
  !*
  !*    The following scalar real arguments have the INTENT(OUT) attribute:
  !*
  !*    sigm                   : mean stress invariant
  !*    dsbar                  : shear stress invariant
  !*    theta                  : parameter in "theta" integrator
  !*
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/

  IMPLICIT NONE
  
  REAL(iwp),INTENT(in)  :: stress(:)
  REAL(iwp),INTENT(out) :: sigm,dsbar,theta

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  REAL(iwp) :: sx,sy,sz,txy,dx,dy,dz,xj3,sine
  REAL(iwp) :: s1,s2,s3,s4,s5,s6,ds1,ds2,ds3,d2,d3,sq3
  INTEGER   :: nst 
  
  nst = ubound(stress,1)
 
  SELECT CASE (nst)
  
    CASE(4)
  
    sx    = stress(1)
    sy    = stress(2)
    txy   = stress(3)
    sz    = stress(4)
    sigm  = (sx+sy+sz)/3._iwp
    dsbar = sqrt((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+6._iwp*txy**2)/sqrt(2._iwp)
    
    IF(dsbar<1.e-10_iwp) THEN
      theta = .0_iwp
    ELSE
      dx   = (2._iwp*sx-sy-sz)/3._iwp
      dy   = (2._iwp*sy-sz-sx)/3._iwp
      dz   = (2._iwp*sz-sx-sy)/3._iwp
      xj3  = dx*dy*dz-dz*txy**2
      sine = -13.5_iwp*xj3/dsbar**3
      
      IF(sine>1._iwp)  sine =  1._iwp
      IF(sine<-1._iwp) sine = -1._iwp
     
      theta = asin(sine)/3._iwp
    END IF

    CASE(6)
   
    sq3   = sqrt(3._iwp)
    s1    = stress(1)  
    s2    = stress(2)
    s3    = stress(3) 
    s4    = stress(4)
    s5    = stress(5)
    s6    = stress(6)
    sigm  = (s1+s2+s3)/3._iwp
    d2    = ((s1-s2)**2+(s2-s3)**2+(s3-s1)**2)/6._iwp+s4*s4+s5*s5+s6*s6
    ds1   = s1-sigm 
    ds2   = s2-sigm  
    ds3   = s3-sigm
    d3    = ds1*ds2*ds3-ds1*s5*s5-ds2*s6*s6-ds3*s4*s4+2._iwp*s4*s5*s6
    dsbar = sq3*sqrt(d2)

    IF(dsbar<1.e-10_iwp)THEN
      theta = 0._iwp
    ELSE
      sine = -3._iwp*sq3*d3/(2._iwp*sqrt(d2)**3)
      IF(sine>1._iwp)  sine =  1._iwp 
      IF(sine<-1._iwp) sine = -1._iwp
      theta = asin(sine)/3._iwp
    END IF
   
    CASE DEFAULT

    Print*,"Wrong size for nst in invar"

  END SELECT

  RETURN
  
  END SUBROUTINE INVAR

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE VMFLOW(stress,dsbar,vmfl)

  !/****f* plasticity/vmflow
  !*  NAME
  !*    SUBROUTINE: vmflow
  !*  SYNOPSIS
  !*    Usage:      CALL vmflow(stress,dsbar,vmfl)
  !*  FUNCTION
  !*    Forms the von Mises flow vector.
  !*  INPUTS
  !*    The following real array argument has the INTENT(IN) attribute:
  !*
  !*    stress                 : stress vector
  !*
  !*    The following scalar real argument has the INTENT(IN) attribute:
  !*
  !*    dsbar                  : shear stress invariant
  !*
  !*    The following real array argument has the INTENT(OUT) attribute:
  !*
  !*    vmfl                   : von Mises "flow" vector
  !*
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/

  IMPLICIT NONE

  REAL(iwp),INTENT(in)  :: stress(:),dsbar
  REAL(iwp),INTENT(out) :: vmfl(:)

  INTEGER   :: nst
  REAL(iwp) :: sigm
  
  nst = ubound(stress,1)
  
  SELECT CASE (nst)

    CASE(4)
    
    sigm    = (stress(1)+stress(2)+stress(4))/3._iwp
    vmfl(1) = stress(1)-sigm
    vmfl(2) = stress(2)-sigm
    vmfl(3) = stress(3)*2._iwp 
    vmfl(4) = stress(4)-sigm
    vmfl    = vmfl*1.5_iwp/dsbar

    CASE(6)
    
    sigm    = (stress(1)+stress(2)+stress(3))/3._iwp
    vmfl(1) = stress(1)-sigm
    vmfl(2) = stress(2)-sigm
    vmfl(3) = stress(3)-sigm
    vmfl(4) = stress(4)*2._iwp 
    vmfl(5) = stress(5)*2._iwp 
    vmfl(6) = stress(6)*2._iwp
    vmfl    = vmfl*1.5_iwp/dsbar

    CASE DEFAULT
    
    PRINT*,"Wrong size for nst in vmflow"

  END SELECT

  RETURN

  END SUBROUTINE VMFLOW

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE FORMAA(flow,nst,rmat,modee)
  
  !/****f* plasticity/formaa
  !*  NAME
  !*    SUBROUTINE: formaa
  !*  SYNOPSIS
  !*    Usage:      CALL formaa(flow,nst,rmat,modee)
  !*  FUNCTION
  !*    Modification to dee matrix.
  !*  INPUTS
  !*    The following real array arguments have the INTENT(IN) attribute:
  !*
  !*    flow(:)                : vector
  !*    rmat(:,:)              : array
  !*
  !*    The following scalar integer argument has the INTENT(IN) attribute:
  !*
  !*    nst                    : number of stress (strain) terms
  !*
  !*    The following real array argument has the INTENT(OUT) attribute:
  !*
  !*    modee(:,:)             : Modified DEE matrix
  !*
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/
  
  IMPLICIT NONE

  REAL(iwp),INTENT(in)  :: flow(:),rmat(:,:)
  REAL(iwp),INTENT(out) :: modee(:,:)
  INTEGER, INTENT(in)   :: nst
  REAL(iwp)             :: flowt(1,nst),flowa(nst,1)
  
  flowt(1,:) = flow 
  flowa(:,1) = flow
  modee      = matmul(matmul(matmul(rmat,flowa),flowt),rmat)

  RETURN
  
  END SUBROUTINE FORMAA

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE FMRMAT(vmfl,nst,dsbar,dlam,dee,temp,rmat)
  
  !/****f* plasticity/fmrmat
  !*  NAME
  !*    SUBROUTINE: fmrmat
  !*  SYNOPSIS
  !*    Usage:      CALL fmrmat(vmfl,nst,dsbar,dlam,dee,temp,rmat)
  !*  FUNCTION
  !*    Forms the r matrix.
  !*  INPUTS
  !*    The following real array arguments have the INTENT(IN) attribute:
  !*
  !*    vmfl(:)                : von Mises "flow" vector
  !*    dee(:,:)               : stress-strain matrix
  !*    temp(:,:)              : temporary array
  !*
  !*    The following scalar real arguments have the INTENT(IN) attribute:
  !*
  !*    dsbar                  : shear stress invariant
  !*    dlam                   : plastic multiplier
  !*
  !*    The following scalar integer argument has the INTENT(IN) attribute:
  !*
  !*    nst                    : number of stress (strain) terms
  !*
  !*    The following real array argument has the INTENT(OUT) attribute:
  !*
  !*    rmat(:,:)              : r matrix
  !*
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/

  IMPLICIT NONE
  
  REAL(iwp),INTENT(in)  :: vmfl(:),dsbar,dlam,dee(:,:),temp(:,:)
  REAL(iwp),INTENT(out) :: rmat(:,:)
  INTEGER,INTENT(in)    :: nst
  
!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER               :: i,j
  REAL(iwp)             :: con
  REAL(iwp)             :: acat(nst,nst),acatc(nst,nst),qmat(nst,nst)

  DO i=1,nst
    DO j=1,nst
      acat(i,j) = vmfl(i)*vmfl(j)
    END DO
  END DO
  
  acat  = (temp-acat)/dsbar
  acatc = matmul(dee,acat)
  qmat  = acatc*dlam

  DO i=1,nst
    qmat(i,i) = qmat(i,i)+1._iwp
  END DO
  
  DO i=1,nst
    
    con       = qmat(i,i)
    qmat(i,i) = 1._iwp
    qmat(i,:) = qmat(i,:)/con
    
    DO j=1,nst
      IF(j/=i) THEN
        con       = qmat(j,i)
        qmat(j,i) = 0.0_iwp
        qmat(j,:) = qmat(j,:) - qmat(i,:)*con
      END IF
    END DO
  END DO
  
  rmat = matmul(qmat,dee)

  RETURN

  END SUBROUTINE FMRMAT

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE FMACAT(vmfl,nst,temp,acat)

  !/****f* plasticity/fmacat
  !*  NAME
  !*    SUBROUTINE: fmacat
  !*  SYNOPSIS
  !*    Usage:      CALL fmacat(vmfl,nst,temp,acat)
  !*  FUNCTION
  !*    Intermediate step.
  !*  INPUTS
  !*    The following real array arguments have the INTENT(IN) attribute:
  !*
  !*    vmfl(:)                : von Mises "flow" vector
  !*    temp(:,:)              : temporary array
  !*
  !*    The following scalar integer argument has the INTENT(IN) attribute:
  !*
  !*    nst                    : number of stress (strain) terms
  !*
  !*    The following real array argument has the INTENT(OUT) attribute:
  !*
  !*    acat(:,:)              : intermediate array
  !*
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/
  
  IMPLICIT NONE

  REAL(iwp),INTENT(in)  :: vmfl(:),temp(:,:)
  INTEGER,INTENT(in)    :: nst
  REAL(iwp),INTENT(out) :: acat(:,:)
  INTEGER               :: i,j
  
  DO i=1,nst
    DO j=1,nst
      acat(i,j) = vmfl(i)*vmfl(j)
    END DO
  END DO
  
  acat = temp - acat

  RETURN
  
  END SUBROUTINE FMACAT

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE VMPL(e,v,stress,pl)

  !/****f* plasticity/vmpl
  !*  NAME
  !*    SUBROUTINE: vmpl
  !*  SYNOPSIS
  !*    Usage:      CALL vmpl(e,v,stress,pl)
  !*  FUNCTION
  !*    Forms plastic matrix for a von-Mises material
  !*  INPUTS
  !*    The following real array argument has the INTENT(IN) attribute:
  !*
  !*    stress(:)              : stress term increments
  !*
  !*    The following scalar real arguments have the INTENT(IN) attribute:
  !*
  !*    e                      : Young's modulus
  !*    v                      : Poisson's ratio
  !*
  !*    The following real array argument has the INTENT(OUT) attribute:
  !*
  !*    pl(:,:)                : plastic matrix
  !*
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/
  
  IMPLICIT NONE
  
  REAL(iwp),INTENT(in)  :: e,v,stress(:)
  REAL(iwp),INTENT(out) :: pl(:,:)
  
!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  REAL(iwp)             :: sx,sy,sz,txy,tyz,tzx,dsbar,ee,term(6)
  INTEGER               :: i,j,nst
  
  nst = ubound(stress,1)

  SELECT CASE (nst)
  
    CASE(4)
 
    sx      = stress(1)
    sy      = stress(2) 
    txy     = stress(3)
    sz      = stress(4)
    dsbar   = sqrt((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+6._iwp*txy**2)/sqrt(2._iwp)
    ee      = 1.5_iwp*e/((1._iwp+v)*dsbar*dsbar)
    term(1) = (2._iwp*sx-sy-sz)/3._iwp
    term(2) = (2._iwp*sy-sz-sx)/3._iwp
    term(3) = txy             
    term(4) = (2._iwp*sz-sx-sy)/3._iwp

    CASE(6)
    
    sx      = stress(1)
    sy      = stress(2)
    sz      = stress(3)
    txy     = stress(4)
    tyz     = stress(5)
    tzx     = stress(6)
    dsbar   = sqrt((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+                            &
              6._iwp*txy**2 + 6._iwp*tyz**2 + 6._iwp*tzx**2)/sqrt(2._iwp)
    ee      = 1.5_iwp*e/((1._iwp+v)*dsbar*dsbar)
    term(1) = (2._iwp*sx-sy-sz)/3._iwp
    term(2) = (2._iwp*sy-sz-sx)/3._iwp
    term(3) = (2._iwp*sz-sx-sy)/3._iwp
    term(4) = txy
    term(5) = tyz
    term(6) = tzx

  END SELECT
  
  DO i=1,nst
    DO j=1,nst
      pl(i,j) = term(i)*term(j)*ee
      pl(j,i) = pl(i,j)
    END DO
  END DO

  RETURN
  
  END SUBROUTINE VMPL 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE FORM_TEMP(temp)
  
  !/****f* plasticity/form_temp
  !*  NAME
  !*    SUBROUTINE: form_temp
  !*  SYNOPSIS
  !*    Usage:      CALL form_temp(temp)
  !*  FUNCTION
  !*    Forms temp matrix for consistent return.
  !*  INPUTS
  !*    The following real array argument has the INTENT(INOUT) attribute:
  !*
  !*    temp(:,:)              : temporary matrix
  !*
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/
  
  REAL(iwp),INTENT(inout) :: temp(:,:)
  INTEGER                 :: nst
  
  nst = ubound(temp,1)
  
  SELECT CASE (nst)

    CASE(4)
    
    temp      =  .0_iwp
    temp(1,1) =  1._iwp
    temp(2,2) =  1._iwp
    temp(4,4) =  1._iwp
    temp(3,3) =  3._iwp
    temp(1,2) =  -.5_iwp
    temp(2,1) =  -.5_iwp
    temp(1,4) =  -.5_iwp
    temp(4,1) =  -.5_iwp
    temp(2,4) =  -.5_iwp
    temp(4,2) =  -.5_iwp
    
    CASE(6)
  
    temp      =  .0_iwp
    temp(1,1) =  1._iwp
    temp(2,2) =  1._iwp
    temp(3,3) =  1._iwp
    temp(4,4) =  3._iwp
    temp(5,5) =  3._iwp
    temp(6,6) =  3._iwp
    temp(1,2) =  -.5_iwp
    temp(2,1) =  -.5_iwp
    temp(1,3) =  -.5_iwp
    temp(3,1) =  -.5_iwp
    temp(2,3) =  -.5_iwp
    temp(3,2) =  -.5_iwp
    
    CASE DEFAULT
    
    PRINT*,"Wrong size for nst in form_temp "
  
  END SELECT

  END SUBROUTINE FORM_TEMP
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE MOCOUF(phi,c,sigm,dsbar,theta,f)
  
  !/****f* plasticity/mocouf
  !*  NAME
  !*    SUBROUTINE: mocouf
  !*  SYNOPSIS
  !*    Usage:      CALL mocouf(phi,c,sigm,dsbar,theta,f)
  !*  FUNCTION
  !*    This subroutine calculates the value of the yield function for a 
  !*    Mohr-Coulomb material (phi in degrees).
  !*  INPUTS
  !*    The following scalar real arguments have the INTENT(IN) attribute:
  !*
  !*    phi                    : friction angle
  !*    c                      : cohesion
  !*    sigm                   : mean stress
  !*    dsbar                  : stress invariant
  !*    theta                  : parameter in "theta" integrator
  !*
  !*    The following scalar real argument has the INTENT(OUT) attribute:
  !*
  !*    f                      : current stress state
  !*
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/

  IMPLICIT NONE

  REAL(iwp),INTENT(in)  :: phi,c,sigm,dsbar,theta   
  REAL(iwp),INTENT(out) :: f

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  REAL(iwp)             :: phir,snph,csph,csth,snth
  REAL(iwp)             :: one=1.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d180=180.0_iwp

  phir = phi*d4*ATAN(one)/d180
  snph = SIN(phir) 
  csph = COS(phir) 
  csth = COS(theta)
  snth = SIN(theta)
  f    = snph*sigm+dsbar*(csth/SQRT(d3)-snth*snph/d3)-c*csph

  RETURN

  END SUBROUTINE MOCOUF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE MOCOUQ(psi,dsbar,theta,dq1,dq2,dq3)
  
  !/****f* plasticity/mocouq
  !*  NAME
  !*    SUBROUTINE: mocouq
  !*  SYNOPSIS
  !*    Usage:      CALL mocouq(psi,dsbar,theta,dq1,dq2,dq3)
  !*  FUNCTION
  !*    This subroutine forms the derivatives of a Mohr-Coulomb potential
  !*    function with respect to the three invariants (psi in degrees).
  !*  INPUTS
  !*    The following scalar real arguments have the INTENT(IN) attribute:
  !*
  !*    psi                    : angle of dilation
  !*    dsbar                  : shear stress invariant
  !*    theta                  : parameter in "theta" integrator
  !*
  !*    The following scalar real arguments have the INTENT(OUT) attribute:
  !*
  !*    dq1                    : Mohr-Coulomb plastic potential derivative
  !*    dq2                    : Mohr-Coulomb plastic potential derivative
  !*    dq3                    : Mohr-Coulomb plastic potential derivative
  !*  
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/

  IMPLICIT NONE
 
  REAL(iwp),INTENT(IN)  :: psi,dsbar,theta
  REAL(iwp),INTENT(OUT) :: dq1,dq2,dq3
  
!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  REAL(iwp)             :: psir,snth,snps,sq3,c1,csth,cs3th,tn3th,tnth
  REAL(iwp)             :: zero=0.0_iwp,pt49=0.49_iwp,pt5=0.5_iwp,one=1.0_iwp
  REAL(iwp)             :: d3=3.0_iwp,d4=4.0_iwp,d180=180.0_iwp
  
  psir  = psi*d4*ATAN(one)/d180 
  snth  = SIN(theta) 
  snps  = SIN(psir)
  sq3   = SQRT(d3)  
  dq1   = snps
  
  IF(ABS(snth).GT.pt49)THEN
    c1 = one
    IF(snth.LT.zero)c1 = -one
    dq2   = (sq3*pt5-c1*snps*pt5/sq3)*sq3*pt5/dsbar 
    dq3   = zero
  ELSE
    csth  = COS(theta)
    cs3th = COS(d3*theta)
    tn3th = TAN(d3*theta)
    tnth  = snth/csth
    dq2   = sq3*csth/dsbar*((one+tnth*tn3th)+snps*(tn3th-tnth)/sq3)*pt5
    dq3   = pt5*d3*(sq3*snth+snps*csth)/(cs3th*dsbar*dsbar)
  END IF
  
  RETURN
  
  END SUBROUTINE MOCOUQ

END MODULE PLASTICITY
