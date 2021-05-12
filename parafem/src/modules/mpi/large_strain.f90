MODULE LARGE_STRAIN

  !/****h* /large_strain
  !*  NAME
  !*    MODULE: large_strain
  !*  SYNOPSIS
  !*    Usage:      USE large_strain
  !*  FUNCTION
  !*    Contains subroutines required for large strains. These subroutines are 
  !*    part of an old version of ParaFEM. They need re-engineering to be made
  !*    compatible with the version of ParaFEM in the Google Code repository. 
  !*    
  !*    Subroutine             Purpose
  !*    KINE3D                 Compute Green-Lagrange deformation tensor (E)
  !*    VENANTKIRCHOFF         Compute second Piola-Kirchoff tensor (S)
  !*    PUSH2R                 Push forward a second order tensor
  !*    PUSH4R                 Push forward a 4 order tensor
  !*    VOIGT2TO1              Write components of matrix in vector form
  !*    VOIGT4TO2              Write components of tensor in matrix form
  !*    GET_GAUSS_POINTS       Compute derivatives of the shape functions 
  !*    CALC_NN_PP             Distribute nodes across processors
  !*    REST_NF                Find node freedom array from rest
  !*    NUM_TO_G2              Assign equation numbers to elements
  !*    CALC_NEQ               Find the number of equations in the problem
  !*    COMPUTE_NPES_PP        Find the value of npes_pp
  !*    PCG_VER1               Preconditioned conjugate gradient solver
  !*  AUTHOR
  !*    F. Calvo Plaza
  !*    V. Szeremi
  !*    L. Margetts
  !*  COPYRIGHT
  !*    2008-2011 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  This is a temporary module that supports program 
  !*  src/programs/5th_ed/xx7.f90
  !*/

  USE mpi_wrapper
  USE PRECISION
  USE MP_INTERFACE
  USE PARTITION
  USE GATHER_SCATTER
  USE LOADING
  USE ELEMENTS
  USE MATHS
  USE NEW_LIBRARY 
  
  CONTAINS

  SUBROUTINE KINE3D(igauss,auxm,coord,points,det,detF,beeF,defE,derivF,jacF)

    !/****f* large_strain/kine3D
    !*  NAME
    !*    SUBROUTINE: kine3D
    !*  SYNOPSIS
    !*    Usage:      CALL kine3D(igauss,auxm,coord,points,det,detF,beeF, &
    !*                            defE,derivF,jacF)
    !*  FUNCTION
    !*    Compute the Green-Lagrange deformation tensor (E) and other 
    !*    kinematics objects in a Gauss point of the deformed element
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    igauss            : Integer
    !*                      : Gauss point number
    !*
    !*    auxm(nod,ndim)    : Real
    !*                      : Nodal displacements
    !*
    !*    coord(nod,ndim)   : Real
    !*                      : Nodal coordinates
    !*
    !*    points(ndim,nip)  : Real
    !*                      : Gauss points coordinates in the reference
    !*                        element
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    det               : Real
    !*                      : Jacobian of the transformation from the reference
    !*                        element to the undeformed element
    !*
    !*    detF              : Real
    !*                      : Jacobian of the transformation from the 
    !*                        undeformed element to the deformed element
    !*
    !*    beeF(nst,ndof)    : Real
    !*                      : Shape function derivatives in the deformed
    !*                        element arranged in B array
    !*
    !*    defE(ndim,ndim)   : Real
    !*                      : Green-Lagrange deformation tensor
    !*
    !*    derivF(ndim,nod)  : Real
    !*                      : Shape function derivatives in the deformed
    !*                        element 
    !*
    !*    jacF(ndim,ndim)   : Real
    !*                      : Deformation gradient tensor
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    28.02.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    IMPLICIT NONE

    INTEGER,   INTENT(IN)  :: igauss
    REAL(iwp), INTENT(IN)  :: auxm(:,:), coord(:,:), points(:,:)
    REAL(iwp), INTENT(OUT) :: det, detF, beeF(:,:), defE(:,:), derivF(:,:), &
                              jacF(:,:)
    INTEGER :: ndim, nod
    REAL(iwp),ALLOCATABLE :: der(:,:), jac(:,:), deriv(:,:), jacFinv(:,:),  &
                             derivFtran(:,:), rightCG(:,:), jacinv(:,:)

    ndim = UBOUND(derivF,1) 
    nod  = UBOUND(derivF,2) 
    ALLOCATE(der(ndim,nod), jac(ndim,ndim), deriv(ndim,nod), jacinv(ndim,ndim))
    ALLOCATE(jacFinv(ndim,ndim), derivFtran(nod,ndim), rightCG(ndim,ndim))
      
!------------------------------------------------------------------------------
! 1. Transformation from the reference element to the undeformed element
!------------------------------------------------------------------------------

    CALL SHAPE_DERIVATIVES(igauss,points,der)
    jac = MATMUL(der,coord)
    CALL DETERMINANT_MATRIX_3x3(jac,det)
    CALL INVERT_MATRIX_3x3(jac,jacinv)
    deriv = MATMUL(jacinv,der)

!------------------------------------------------------------------------------
! 2. Transformation from the undeformed element to the deformed element
!------------------------------------------------------------------------------

    jacF = MATMUL(TRANSPOSE(auxm),TRANSPOSE(deriv))
    jacF(1,1) = jacF(1,1) + 1.0
    jacF(2,2) = jacF(2,2) + 1.0
    jacF(3,3) = jacF(3,3) + 1.0
    CALL DETERMINANT_MATRIX_3x3(jacF,detF)
    CALL INVERT_MATRIX_3x3(jacF,jacFinv)
    derivFtran = MATMUL(TRANSPOSE(deriv),jacFinv)
    derivF = TRANSPOSE(derivFtran)
!   CALL BEEMAT(derivF,beeF)
    CALL BEEMAT(beeF,derivF) ! correct argument order
    rightCG = MATMUL(TRANSPOSE(jacF),jacF)
    defE = 0.5*rightCG
    defE(1,1) = defE(1,1) - 0.5
    defE(2,2) = defE(2,2) - 0.5
    defE(3,3) = defE(3,3) - 0.5
      
  END SUBROUTINE KINE3D

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE VENANTKIRCHHOFF(defE,emod,nu,piolaS,cmat)

    !/****f* large_strain/venantkirchhoff
    !*  NAME
    !*    SUBROUTINE: venantkirchhoff
    !*  SYNOPSIS
    !*    Usage:      CALL venantkirchhoff(defE,emod,nu,piolaS,cmat)
    !*  FUNCTION
    !*    Compute the Second Piola-Kirchhoff tensor (S) and the material
    !*    elasticity tensor (cmat) from the Green-Lagrange deformation
    !*    tensor (E) for a Saint Venant-Kirchhoff material
    !*     
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    defE(ndim,ndim)           : Real
    !*                              : Green-Lagrange deformation tensor
    !*    
    !*    emod                      : Real
    !*                              : Young's modulus
    !*
    !*    nu                        : Real
    !*                              : Poisson's coefficient
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    piolaS(ndim,ndim)         : Real
    !*                              : Piola-Kirchhoff second stress tensor
    !*
    !*    cmat(ndim,ndim,ndim,ndim) : Real
    !*                              : Material elasticity tensor
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    10.04.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  For sake of programming convenience, inditial form is used for the
    !*  operations
    !*
    !*  defE, piolaS and cmat are symmetric, but in this version
    !*  saving memory is not the main objective, so the matrix
    !*  form is used instead the vector form (Voigt notation) 
    !*
    !*  WHAT FOLLOWS NOW IS JUST FOR SAKE OF CLARITY AND UNDERSTANDING:
    !*
    !*  The material elasticity tensor computed in the loop before is:      
    !*    cmat_ABCD = partial(S_AB)/partial(E_CD) (Holzapfel, pag 252)
    !*  with the following properties:
    !*    1) minor symmetries:
    !*       cmat_ABCD = cmat_BACD = cmat_ABDC = cmat_BADC
    !*    2) MAJOR symmetries:
    !*       cmat_ABCD = cmat_CDAB
    !*  Due to major + minor symmetries ----> 21 independent components
    !*  (Holzapfel, pags 252-253) 
    !*  In what follows, it can be seen that only 21 components are defined
    !*  independently:
    !*
    !*  cmat(1,1,1,1) = lambda + 2*mu 
    !*  cmat(1,1,1,2) = 0.0
    !*  cmat(1,1,1,3) = 0.0
    !*  cmat(1,1,2,1) = cmat(1,1,1,2)  !minor
    !*  cmat(1,1,2,2) = lambda 
    !*  cmat(1,1,2,3) = 0.0 
    !*  cmat(1,1,3,1) = cmat(1,1,1,3)  !minor
    !*  cmat(1,1,3,2) = cmat(1,1,2,3)  !minor
    !*  cmat(1,1,3,3) = lambda 
    !*    
    !*  cmat(1,2,1,1) = cmat(1,1,1,2)  !MAJOR
    !*  cmat(1,2,1,2) = mu
    !*  cmat(1,2,1,3) = 0.0
    !*  cmat(1,2,2,1) = cmat(1,2,1,2)  !minor
    !*  cmat(1,2,2,2) = 0.0 
    !*  cmat(1,2,2,3) = 0.0 
    !*  cmat(1,2,3,1) = cmat(1,2,1,3)  !minor
    !*  cmat(1,2,3,2) = cmat(1,2,2,3)  !minor
    !*  cmat(1,2,3,3) = 0.0
    !*
    !*  cmat(1,3,1,1) = cmat(1,1,1,3)  !MAJOR
    !*  cmat(1,3,1,2) = cmat(1,2,1,3)  !MAJOR 
    !*  cmat(1,3,1,3) = mu
    !*  cmat(1,3,2,1) = cmat(1,3,1,2)  !minor
    !*  cmat(1,3,2,2) = 0.0 
    !*  cmat(1,3,2,3) = 0.0 
    !*  cmat(1,3,3,1) = cmat(1,3,1,3)  !minor
    !*  cmat(1,3,3,2) = cmat(1,3,2,3)  !minor
    !*  cmat(1,3,3,3) = 0.0
    !*
    !*  cmat(2,1,1,1) = cmat(1,2,1,1)  !minor
    !*  cmat(2,1,1,2) = cmat(1,2,1,2)  !minor
    !*  cmat(2,1,1,3) = cmat(1,2,1,3)  !minor 
    !*  cmat(2,1,2,1) = cmat(1,2,2,1)  !minor
    !*  cmat(2,1,2,2) = cmat(1,2,2,2)  !minor 
    !*  cmat(2,1,2,3) = cmat(1,2,3,2)  !minor
    !*  cmat(2,1,3,1) = cmat(1,2,1,3)  !minor
    !*  cmat(2,1,3,2) = cmat(1,2,2,3)  !minor
    !*  cmat(2,1,3,3) = cmat(1,2,3,3)  !minor
    !*    
    !*  cmat(2,2,1,1) = cmat(1,1,2,2)  !MAJOR
    !*  cmat(2,2,1,2) = cmat(1,2,2,2)  !MAJOR
    !*  cmat(2,2,1,3) = cmat(1,3,2,2)  !MAJOR
    !*  cmat(2,2,2,1) = cmat(2,1,2,2)  !MAJOR
    !*  cmat(2,2,2,2) = lambda + 2*mu
    !*  cmat(2,2,2,3) = 0.0 
    !*  cmat(2,2,3,1) = cmat(2,2,1,3)  !minor
    !*  cmat(2,2,3,2) = cmat(2,2,2,3)  !minor
    !*  cmat(2,2,3,3) = lambda
    !*    
    !*  cmat(2,3,1,1) = cmat(1,1,2,3)  !MAJOR
    !*  cmat(2,3,1,2) = cmat(1,2,2,3)  !MAJOR 
    !*  cmat(2,3,1,3) = cmat(1,3,2,3)  !MAJOR 
    !*  cmat(2,3,2,1) = cmat(2,1,2,3)  !MAJOR
    !*  cmat(2,3,2,2) = cmat(2,2,2,3)  !MAJOR 
    !*  cmat(2,3,2,3) = mu
    !*  cmat(2,3,3,1) = cmat(1,3,2,3)  !MAJOR + minor
    !*  cmat(2,3,3,2) = cmat(2,3,2,3)  !minor
    !*  cmat(2,3,3,3) = 0.0
    !*    
    !*  cmat(3,1,1,1) = cmat(1,3,1,1)  !minor
    !*  cmat(3,1,1,2) = cmat(1,3,1,2)  !minor
    !*  cmat(3,1,1,3) = cmat(1,3,1,3)  !minor 
    !*  cmat(3,1,2,1) = cmat(1,3,2,1)  !minor
    !*  cmat(3,1,2,2) = cmat(1,3,2,2)  !minor 
    !*  cmat(3,1,2,3) = cmat(1,3,2,3)  !minor
    !*  cmat(3,1,3,1) = cmat(1,3,3,1)  !minor
    !*  cmat(3,1,3,2) = cmat(1,3,3,2)  !minor
    !*  cmat(3,1,3,3) = cmat(1,3,3,3)  !minor
    !*    
    !*  cmat(3,2,1,1) = cmat(2,3,1,1)  !minor
    !*  cmat(3,2,1,2) = cmat(2,3,1,2)  !minor
    !*  cmat(3,2,1,3) = cmat(2,3,1,3)  !minor
    !*  cmat(3,2,2,1) = cmat(2,3,2,1)  !minor
    !*  cmat(3,2,2,2) = cmat(2,3,2,2)  !minor 
    !*  cmat(3,2,2,3) = cmat(2,3,2,3)  !minor 
    !*  cmat(3,2,3,1) = cmat(2,3,3,1)  !minor
    !*  cmat(3,2,3,2) = cmat(2,3,3,2)  !minor
    !*  cmat(3,2,3,3) = cmat(2,3,3,3)  !minor 
    !*    
    !*  cmat(3,3,1,1) = cmat(1,1,3,3)  !MAJOR
    !*  cmat(3,3,1,2) = cmat(1,2,3,3)  !MAJOR 
    !*  cmat(3,3,1,3) = cmat(1,3,3,3)  !MAJOR 
    !*  cmat(3,3,2,1) = cmat(2,1,3,3)  !MAJOR
    !*  cmat(3,3,2,2) = cmat(2,2,3,3)  !MAJOR
    !*  cmat(3,3,2,3) = cmat(2,3,3,3)  !MAJOR
    !*  cmat(3,3,3,1) = cmat(3,1,3,3)  !MAJOR
    !*  cmat(3,3,3,2) = cmat(3,2,3,3)  !MAJOR
    !*  cmat(3,3,3,3) = lambda + 2*mu
    !*    
    !*  DO jA = 1,3
    !*    DO jB = 1,3
    !*      DO jC = 1,3
    !*        DO jD = 1,3
    !*          WRITE(22,*)jA,jB,jC,jD,"   ",cmat(jA,jB,jC,jD)
    !*        END DO
    !*      END DO
    !*    END DO
    !*  END DO
    !*    
    !*  Alternative forms to compute piolaS (just for information):
    !*  a) Tensorial/indicial form:  piolaS_AB = cmat_ABCD*defE_CD
    !*    (Belytschko, pag 228)
    !*    piolaS = 0.0
    !*    DO jA = 1,3
    !*      DO jB = 1,3
    !*        DO jC = 1,3
    !*          DO jD = 1,3
    !*            piolaS(jA,jB) = piolaS(jA,jB) + cmat(jA,jB,jC,jD)*defE(jC,jD)
    !*          END DO
    !*        END DO
    !*      END DO
    !*      WRITE(*,*)(piolaS(jA,jB),jB=1,3)
    !*    END DO
    !*    
    !*  b) Matricial form: D_11 = cmat_1111, D_12 = cmat_1122, etc... 
    !*    (Belytschko, pag 227): with the following Voigt notation
    !*    (slightly different from Belytschko, but the same that in FEAP):
    !*    11    ---> 1
    !*    22    ---> 2 
    !*    33    ---> 3
    !*    12,21 ---> 4 
    !*    23,32 ---> 5 
    !*    13,31 ---> 6 
    !*    
    !*    |S_11|   | D_11  D_12  D_13   0     0     0   | |E_11| 
    !*    |S_22|   |       D_22  D_23   0     0     0   | |E_22|
    !*    |S_33|   |             D_33   0     0     0   | |E_33|
    !*    |S_12| = |                   D_44   0     0   | |2*E_12|
    !*    |S_23|   |     Symmetry            D_55   0   | |2*E_23|
    !*    |S_13|   |                               D_66 | |2*E_13|
    !*    
    !*    piolaS(1,1) = cmat(1,1,1,1)*defE(1,1) + &   !D_11*defE(1,1)
    !*                  cmat(1,1,2,2)*defE(2,2) + &   !D_12*defE(2,2) 
    !*                  cmat(1,1,3,3)*defE(3,3)       !D_13*defE(3,3)
    !*    piolaS(2,2) = cmat(2,2,1,1)*defE(1,1) + &   !D_21*defE(1,1)
    !*                  cmat(2,2,2,2)*defE(2,2) + &   !etc...
    !*                  cmat(2,2,3,3)*defE(3,3)
    !*    piolaS(3,3) = cmat(3,3,1,1)*defE(1,1) + &
    !*                  cmat(3,3,2,2)*defE(2,2) + &
    !*                  cmat(3,3,3,3)*defE(3,3)
    !*    piolaS(1,2) = cmat(1,2,1,2)*defE(1,2)*2.0 
    !*    piolaS(2,3) = cmat(2,3,2,3)*defE(2,3)*2.0 
    !*    piolaS(1,3) = cmat(1,3,1,3)*defE(1,3)*2.0 
    !*    piolaS(2,1) = piolaS(2,1)
    !*    piolaS(3,2) = piolaS(2,3)
    !*    piolaS(3,1) = piolaS(1,3)
    !*    DO jA = 1,3
    !*      WRITE(*,*)(piolaS(jA,jB),jB=1,3)
    !*    END DO
    !*
    !*/

    IMPLICIT NONE

    REAL(iwp), INTENT(IN)  :: defE(:,:), emod, nu
    REAL(iwp), INTENT(OUT) :: piolaS(:,:), cmat(:,:,:,:)
    INTEGER   :: jA, jB, jC, jD
    REAL(iwp) :: traceE, lambda, mu, tensorI(3,3), delta(3,3)

 !-----------------------------------------------------------------------------
 ! 1.  Piola-Kirchhoff second stress tensor is computed using the following 
 !     indicial formula:
 !     piolaS_AB = lambda*tr(E)I_AB + 2*mu*E_AB 
 !     See (Belytschko, pag 228), (Bonet, pag 120)
 !-----------------------------------------------------------------------------

    lambda = nu*emod/((1.0+nu)*(1.0-2.0*nu))
    mu     = emod/(2.0*(1+nu))
    traceE = defE(1,1) + defE(2,2) + defE(3,3)
      
    tensorI      = 0.0    !2-order identity tensor
    tensorI(1,1) = 1.0
    tensorI(2,2) = 1.0
    tensorI(3,3) = 1.0
      
    DO jA = 1,3
      DO jB = 1,3
        piolaS(jA,jB) = lambda*traceE*tensorI(jA,jB) + 2.0*mu*defE(jA,jB)
      END DO
    END DO

!------------------------------------------------------------------------------
! 2.  Material elasticity tensor is computed using the following indicial
!     formula:
!     cmat_ABCD = lambda*delta_AB*delta_CD + 
!                 mu*(delta_AC*delta_BD + delta_AD*delta_BC)
!     See (Belytschko, pag 228)
!------------------------------------------------------------------------------

    delta      = 0.0      !Kronecker delta
    delta(1,1) = 1.0
    delta(2,2) = 1.0
    delta(3,3) = 1.0

    DO jA = 1,3
      DO jB = 1,3
        DO jC = 1,3
          DO jD = 1,3
            cmat(jA,jB,jC,jD) = lambda*tensorI(jA,jB)*tensorI(jC,jD) + &
                mu*(delta(jA,jC)*delta(jB,jD) + delta(jA,jD)*delta(jB,jC))
          END DO
        END DO
      END DO
    END DO
      
  END SUBROUTINE VENANTKIRCHHOFF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE PUSH2R(jac_X_x,defF,piolaS,sigma)

    !/****f* large_strain/push2r
    !*  NAME
    !*    SUBROUTINE: push2r
    !*  SYNOPSIS
    !*    Usage:      CALL push2r(jac_X_x,defF,piolaS,sigma)
    !*  FUNCTION
    !*    Push-forward a 2-order tensor
    !*    sigma = (1.0/jac_X_x)*defF*piolaS*defF^T
    !*    
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    jac_X_x              : Real
    !*                         : Jacobian from underformed to current element
    !*
    !*    defF(ndim,ndim)      : Real
    !*                         : Deformation gradient tensor
    !*
    !*    piolaS(ndim,ndim)    : Real
    !*                         : Piola-Kirchhoff second stress tensor
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    sigma(ndim,ndim)     : Real
    !*                         : Cauchy stress tensor
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    02.11.2006
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*    piolaS and sigma are symmetric, but in this first version
    !*    saving memory is not the main objective, so the matrix
    !*    form is used instead the vector form (Voigt notation) 
    !*
    !*/

    IMPLICIT NONE

    REAL(iwp), INTENT(IN)  :: jac_X_x, defF(:,:), piolaS(:,:)
    REAL(iwp), INTENT(OUT) :: sigma(:,:)
    INTEGER :: ia, ib, jA, jB

!------------------------------------------------------------------------------
! 1. Push-forward is computed using the indicial formula:
!    sigma_ab = (1.0/jac_X_x) * F_aA * F_bB * S_AB
!------------------------------------------------------------------------------

    sigma = 0.0
    DO ia = 1,3
      DO ib = 1,3
        DO jA = 1,3
          DO jB = 1,3
            sigma(ia,ib) = sigma(ia,ib) + &
                           defF(ia,jA)*defF(ib,jB)*piolaS(jA,jB)
          END DO
        END DO
        sigma(ia,ib) = sigma(ia,ib)/jac_X_x
      END DO
    END DO   
      
  END SUBROUTINE PUSH2R

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE PUSH4R(jac_X_x,defF,cmat,cspa)

    !/****f* large_strain/push4r
    !*  NAME
    !*    SUBROUTINE: push4r
    !*  SYNOPSIS
    !*    Usage:      CALL push4r(jac_X_x,defF,cmat,cspa)
    !*  FUNCTION
    !*    Push-forward a 4-order tensor
    !*    cspa = push-forward(cmat)
    !*    
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    jac_X_x                   : Real
    !*                              : Jacobian from underformed to current
    !*                                element
    !*
    !*    defF(ndim,ndim)           : Real
    !*                              : Deformation gradient tensor
    !*
    !*    cmat(ndim,ndim,ndim,ndim) : Real
    !*                              : Material elasticity tensor
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    cspa(ndim,ndim,ndim,ndim) : Real
    !*                              : Spatial elasticity tensor
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    02.11.2006
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  The elasticity tensor has minor symmetries and in this
    !*  case (hyperelasticity hypothesis) also major symmetries,
    !*  so many components are repeated. However, in this version
    !*  the objective is not to save memory, but to program it clear.

    IMPLICIT NONE

    REAL(iwp), INTENT(IN)  :: jac_X_x, defF(:,:), cmat(:,:,:,:)
    REAL(iwp), INTENT(OUT) :: cspa(:,:,:,:)
    INTEGER :: ia, ib, ic, id, jA, jB, jC, jD

!------------------------------------------------------------------------------
! 1. Push-forward is computed using the indicial formula:
!    cspa_abcd = (1.0/jac_X_x) * F_aA * F_bB * F_cC * F_dD * cmat_ABCD
!------------------------------------------------------------------------------
      
    cspa = 0.0
    DO ia = 1,3
      DO ib = 1,3
        DO ic = 1,3
          DO id = 1,3
            DO jA = 1,3
              DO jB = 1,3
                DO jC = 1,3
                  DO jD = 1,3
                    cspa(ia,ib,ic,id) = cspa(ia,ib,ic,id) + &
           defF(ia,jA)*defF(ib,jB)*defF(ic,jC)*defF(id,jD)*cmat(jA,jB,jC,jD)
                  END DO
                END DO
              END DO
            END DO
            cspa(ia,ib,ic,id) = cspa(ia,ib,ic,id)/jac_X_x
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE PUSH4R

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE VOIGT2TO1(sigma,sigma1C) 

    !/****f* large_strain/voigt2to1
    !*  NAME
    !*    SUBROUTINE: voigt2to1
    !*  SYNOPSIS
    !*    Usage:      CALL voigt2to1(sigma,sigma1C)
    !*  FUNCTION
    !*    Write in vector form the 6 independent components of a symmetric
    !*    matrix
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    sigma(ndim,ndim)       : Real
    !*                           : Symmetric matrix (i.e. Cauchy stress tensor)
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    sigma1C(nst)           : Real   
    !*                           : Vector containing the independent
    !*                             components of a symmetric matrix
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    28.02.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  Look Bonet's book (page 176), although the Voigt notation is changed
    !*  (5 and 6 components are swapped)
    !*
    !*/

    IMPLICIT NONE

    REAL(iwp), INTENT(IN)  :: sigma(:,:)
    REAL(iwp), INTENT(OUT) :: sigma1C(:)
    INTEGER :: ndim

    ndim = UBOUND(sigma,1) 
      
    IF (ndim==2) THEN
      sigma1C(1) = sigma(1,1)
      sigma1C(2) = sigma(2,2)
      sigma1C(3) = sigma(1,2)
    ELSE IF (ndim==3) THEN
      sigma1C(1) = sigma(1,1)
      sigma1C(2) = sigma(2,2)
      sigma1C(3) = sigma(3,3)
      sigma1C(4) = sigma(1,2)
      sigma1C(5) = sigma(2,3)
      sigma1C(6) = sigma(1,3)
     END IF
      
  END SUBROUTINE voigt2to1

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

  SUBROUTINE voigt4to2(cspa,deeF) 

    !/****f* large_strain/voigt4to2
    !*  NAME
    !*    SUBROUTINE: voigt4to2
    !*  SYNOPSIS
    !*    Usage:      CALL voigt4to2(cspa,deeF)
    !*  FUNCTION
    !*    Write in matrix form the 21 independent components of a 4-order
    !*    tensor with major and minor symmetries
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    cspa(ndim,ndim,ndim,ndim) : Real
    !*                              : 4-order tensor with major and minor
    !*                                symmetries (i.e. spatial elasticity
    !*                                tensor)
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    deeF(nst,nst)              : Real   
    !*                               : Matrix containing the independent
    !*                                 components of a 4-order tensor with
    !*                                 major and minor symmetries
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    28.02.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  Look Bonet's book (page 176)
    !*
    !*  The subroutine should be completed with the expression for 2D problems
    !*
    !*/

    IMPLICIT NONE

    REAL(iwp), INTENT(IN)  :: cspa(:,:,:,:)
    REAL(iwp), INTENT(OUT) :: deeF(:,:)
    INTEGER :: ndim

    ndim = UBOUND(cspa,1) 
      
    IF (ndim==3) THEN
      deeF(1,1) = cspa(1,1,1,1)
      deeF(1,2) = cspa(1,1,2,2)
      deeF(1,3) = cspa(1,1,3,3)
      deeF(1,4) = 0.5*(cspa(1,1,1,2) + cspa(1,1,2,1))
      deeF(1,5) = 0.5*(cspa(1,1,2,3) + cspa(1,1,3,2))
      deeF(1,6) = 0.5*(cspa(1,1,1,3) + cspa(1,1,3,1))

      deeF(2,1) = deeF(1,2)      !symmetry
      deeF(2,2) = cspa(2,2,2,2)
      deeF(2,3) = cspa(2,2,3,3)
      deeF(2,4) = 0.5*(cspa(2,2,1,2) + cspa(2,2,2,1))
      deeF(2,5) = 0.5*(cspa(2,2,2,3) + cspa(2,2,3,2))
      deeF(2,6) = 0.5*(cspa(2,2,1,3) + cspa(2,2,3,1))

      deeF(3,1) = deeF(1,3)      !symmetry
      deeF(3,2) = deeF(2,3)      !symmetry
      deeF(3,3) = cspa(3,3,3,3)
      deeF(3,4) = 0.5*(cspa(3,3,1,2) + cspa(3,3,2,1))
      deeF(3,5) = 0.5*(cspa(3,3,2,3) + cspa(3,3,3,2))
      deeF(3,6) = 0.5*(cspa(3,3,1,3) + cspa(3,3,3,1))

      deeF(4,1) = deeF(1,4)      !symmetry
      deeF(4,2) = deeF(2,4)      !symmetry
      deeF(4,3) = deeF(3,4)      !symmetry 
      deeF(4,4) = 0.5*(cspa(1,2,1,2) + cspa(1,2,2,1))
      deeF(4,5) = 0.5*(cspa(1,2,2,3) + cspa(1,2,3,2))
      deeF(4,6) = 0.5*(cspa(1,2,1,3) + cspa(1,2,3,1))

      deeF(5,1) = deeF(1,5)      !symmetry
      deeF(5,2) = deeF(2,5)      !symmetry
      deeF(5,3) = deeF(3,5)      !symmetry
      deeF(5,4) = deeF(4,5)      !symmetry
      deeF(5,5) = 0.5*(cspa(2,3,2,3) + cspa(2,3,3,2))
      deeF(5,6) = 0.5*(cspa(2,3,1,3) + cspa(2,3,3,1))

      deeF(6,1) = deeF(1,6)      !symmetry
      deeF(6,2) = deeF(2,6)      !symmetry
      deeF(6,3) = deeF(3,6)      !symmetry
      deeF(6,4) = deeF(4,6)      !symmetry
      deeF(6,5) = deeF(5,6)      !symmetry
      deeF(6,6) = 0.5*(cspa(1,3,1,3) + cspa(1,3,3,1))
    END IF

  END SUBROUTINE voigt4to2
  
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE GET_GAUSS_POINTS(element,gausspoints,weights)

    !/****f* large_strain/get_gauss_points
    !*  NAME
    !*    SUBROUTINE: get_gauss_points
    !*  SYNOPSIS
    !*    Usage:      CALL get_gauss_points(element,gausspoints,weights)
    !*  FUNCTION
    !*    Compute the derivatives of the shape functions at a Gauss point
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    element            : Character
    !*                       : Element type (e.g. "tetrahedron", "hexahedron")
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    gausspoints(ndim,nip)   : Real
    !*                            : Gauss points coordinates at the reference 
    !*                              element
    !*
    !*    weights(nip)            : Real
    !*                            : Weight of the Gauss points
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    24.01.2008
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  Although Gauss points don't need to be sorted, we want to follow the
    !*  criterium of giving the Gauss points in the same order as the element
    !*  connectivity when possible (e.g. 4 integration points for tetrahedra,
    !*  8 integration points for hexahedra)
    !*   
    !*  This routine has to be completed with elements for 1 and 2 dimensions  
    !*  and 3D elements with unusual number of nodes
    !*
    !*  NEED TO REPLACE THIS WITH THE SMITH AND GRIFFITHS VERSION
    !*/
     
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)  :: element
    REAL(iwp),    INTENT(OUT) :: gausspoints(:,:), weights(:) 
    INTEGER :: nip
    REAL(iwp) :: root3inv

    nip = UBOUND(gausspoints,2)

    SELECT CASE(element)
      CASE('tetrahedron')
        SELECT CASE(nip)
          CASE(1)
            weights(1)       = 1.0/6.0    !Gauss point 1, weight
            gausspoints(1,1) = 0.25       !Gauss point 1, coordinate x
            gausspoints(2,1) = 0.25       !Gauss point 1, coordinate y
            gausspoints(3,1) = 0.25       !Gauss point 1, coordinate z
          CASE(4)
            weights(:)       =  0.25/6.0
            gausspoints(1,1) = .58541020  !Gauss point 1, coordinate x
            gausspoints(2,1) = .13819660  !Gauss point 1, coordinate y
            gausspoints(3,1) = .13819660  !Gauss point 1, coordinate z
            gausspoints(1,2) = .13819660  !Gauss point 2, coordinate x
            gausspoints(2,2) = .58541020  !Gauss point 2, coordinate y
            gausspoints(3,2) = .13819660  !Gauss point 2, coordinate z
            gausspoints(1,3) = .13819660  !etc...
            gausspoints(2,3) = .13819660
            gausspoints(3,3) = .13819660
            gausspoints(1,4) = .13819660
            gausspoints(2,4) = .13819660
            gausspoints(3,4) = .58541020
          CASE DEFAULT
            WRITE(*,*)"Wrong number of integrating points!!"
        END SELECT
      CASE('hexahedron')
        SELECT CASE(nip)
          CASE(1)
            weights(1)       = 8.0        !Gauss point 1, weight
            gausspoints(1,1) = 0.0        !Gauss point 1, coordinate x
            gausspoints(2,1) = 0.0        !Gauss point 1, coordinate y
            gausspoints(3,1) = 0.0        !Gauss point 1, coordinate z
          CASE(8)
            root3inv = 1.0/sqrt(3.0)
            weights(:)       =  1.0
            gausspoints(1,1) =  root3inv  !Gauss point 1, coordinate x
            gausspoints(2,1) = -root3inv  !Gauss point 1, coordinate y
            gausspoints(3,1) = -root3inv  !Gauss point 1, coordinate z
            gausspoints(1,2) =  root3inv  !Gauss point 2, coordinate x
            gausspoints(2,2) =  root3inv  !Gauss point 2, coordinate y
            gausspoints(3,2) = -root3inv  !Gauss point 2, coordinate z
            gausspoints(1,3) = -root3inv  !etc...
            gausspoints(2,3) =  root3inv
            gausspoints(3,3) = -root3inv
            gausspoints(1,4) = -root3inv
            gausspoints(2,4) = -root3inv
            gausspoints(3,4) = -root3inv
            gausspoints(1,5) =  root3inv
            gausspoints(2,5) = -root3inv
            gausspoints(3,5) =  root3inv
            gausspoints(1,6) =  root3inv
            gausspoints(2,6) =  root3inv
            gausspoints(3,6) =  root3inv
            gausspoints(1,7) = -root3inv
            gausspoints(2,7) =  root3inv
            gausspoints(3,7) =  root3inv
            gausspoints(1,8) = -root3inv
            gausspoints(2,8) = -root3inv
            gausspoints(3,8) =  root3inv
          CASE DEFAULT
            WRITE(*,*)"Wrong number of integrating points!!"
        END SELECT
      CASE DEFAULT
        WRITE(*,*)"Wrong element type!!"
    END SELECT

    RETURN

  END SUBROUTINE GET_GAUSS_POINTS
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE CALC_NN_PP(g_num_pp,nn_pp,nn_start)
    
    !/****f* decomposition/calc_nn_pp
    !*  NAME
    !*    SUBROUTINE: calc_nn_pp
    !*  SYNOPSIS
    !*    Usage:      CALL calc_nn_pp(g_num_pp,nn_pp,nn_start)
    !*  FUNCTION
    !*    Returns number of nodes "nn_pp" and starting node number "nn_start"
    !*    on the calling process
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:    
    !*
    !*    g_num_pp(nod,nels_pp)  : Integer
    !*                           : Elements connectivity
    !*
    !*    The following arguments have the INTENT(OUT) attribute:    
    !*
    !*    nn_pp                  : Integer
    !*                           : Range of nodes on the calling process
    !*
    !*    nn_start               : Integer
    !*                           : First node number on the calling process
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    11.02.2008
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  What follows is an example to understand the routine:
    !*  
    !*  If the nodes of all the elements of a process are:
    !*  3, 4, 7, 8, 11, 12, 16, 17, 20, 25, 50, 51, 52, 53, 62, 65, 66, 69, 102
    !*
    !* then the result is:
    !*
    !* nn_start = 3 
    !* nn_end   = 102 
    !* nn       = 100 
    !*
    !* THIS SUBROUTINE IS ONLY USED IN PROGRAM XX7 AND NEEDS TO BE REMOVED 
    !* FROM THE LATEST VERSION OF PARAFEM
    !*/

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: g_num_pp(:,:)
    INTEGER, INTENT(OUT) :: nn_pp, nn_start
    INTEGER :: nels_pp, iel, nn_end, nodemin, nodemax

    nels_pp  = UBOUND(g_num_pp,2)
    nn_start = MINVAL(g_num_pp(:,1))
    nn_end   = MAXVAL(g_num_pp(:,1))

    DO iel = 2,nels_pp
      nodemin = MINVAL(g_num_pp(:,iel))
      nodemax = MAXVAL(g_num_pp(:,iel))
      IF(nodemin < nn_start) THEN
        nn_start = nodemin
      END IF
      IF(nodemax > nn_end) THEN
        nn_end   = nodemax
      END IF
    END DO
    nn_pp = (nn_end - nn_start) + 1

    RETURN

  END SUBROUTINE CALC_NN_PP
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE REST_NF(nn_start,rest,nf_pp)

    !/****f* large_strain/rest_nf
    !*  NAME
    !*    SUBROUTINE: rest_nf
    !*  SYNOPSIS
    !*    Usage:      CALL rest_nf(rest,nf_pp)
    !*  FUNCTION
    !*    Assign the global degree of freedom number (equation number) to the
    !*    degrees of freedom of the nodes in the process
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:    
    !*
    !*    nn_start            : Integer
    !*                        : First node number in a process
    !*
    !*    rest(nr,nodof+1)    : Integer
    !*                        : List of the nodes with some degree of freedom
    !*                          fixed. In the input, the degrees of freedom
    !*                          fixed have 1 and the degrees of freedom not
    !*                          fixed have 0. In the output, the degrees of
    !*                          freedom fixed have 0 and the degrees of freedom
    !*                          not fixed have the global equation number
    !*
    !*    The following arguments have the INTENT(OUT) attribute:    
    !*
    !*    nf_pp(nodof,nn_pp)  : Integer
    !*                        : Nodes where the equation number has been
    !*                          assigned to each degree of freedom of the
    !*                          nodes
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    07.02.2008
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  What follows is an example to understand the routine:
    !*  
    !*  Let's say that we enter with this array with: (nr=5), (0 free, 1 fixed)
    !*   Node     Disp_X  Disp_Y  Disp_Z
    !*    1         1       1       0
    !*    2         1       0       0
    !*    6         1       0       0
    !*   10         0       1       1
    !*   15         1       0       0
    !*  
    !*  This means that the displacement of node 1 has been fixed in the
    !*  direction X and Y but is free in Z. The displacement of node 2 has been
    !*  fixed in the direction X but is free in Y and Z, etc..
    !*  
    !*  If we use 3 processes, with this range of nodes on each process:
    !*    process 1:   1 -  9           (nn_ff =  9)
    !*    process 2:   4 - 13           (nn_ff = 10)
    !*    process 3:   8 - 16           (nn_ff =  9)
    !* 
    !*  the result on each process is:
    !* 
    !*                     PROCESS 1
    !*  Node1  Node2  Node3  Node4  Node5  Node6  Node7  Node8  Node9   
    !*    0      0      4      7     10      0     15      18    21 
    !*    0      2      5      8     11     13     16      19    22
    !*    1      3      6      9     12     14     17      20    23
    !* 
    !*                     PROCESS 2
    !*  Node4  Node5  Node6  Node7  Node8  Node9  Node10  Node11  Node12 Node13
    !*    7     10      0      15     18     21     24      25      28     31  
    !*    8     11     13      16     19     22      0      26      29     32
    !*    9     12     14      17     20     23      0      27      30     33
    !* 
    !*                     PROCESS 3
    !*  Node8  Node9  Node10  Node11  Node12  Node13  Node14  Node15  Node16
    !*    18     21     24      25      28      31      34       0      39
    !*    19     22      0      26      29      32      35      37      40
    !*    20     23      0      27      30      33      36      38      41
    !* 
    !*  THIS SUBROUTINE IS ONLY USED BY PROGRAM XX7 AND NEEDS TO BE REMOVED
    !*  FROM PARAFEM
    !*/    

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: nn_start, rest(:,:)
    INTEGER, INTENT(OUT) :: nf_pp(:,:)
    INTEGER :: i, j, k, nr, nodof, nn_pp, nn_end, numdof, contrest, noderest
    
    nr     = UBOUND(rest,1)
    nodof  = UBOUND(rest,2) - 1
    nn_pp  = UBOUND(nf_pp,2)
    nn_end = nn_start + nn_pp - 1 

    !-----------------------------------------------------------------------
    ! 1. Find global number of degrees of freedom before the first node of 
    !    the process
    !-----------------------------------------------------------------------

    numdof   = 0
    contrest = 1
    noderest = rest(contrest,1)

    DO i = 1,nn_start-1
      IF (i==noderest) THEN
        DO j = 2,nodof+1 
          IF(rest(contrest,j)/=1) THEN
            numdof = numdof + 1
          END IF
        END DO
        contrest = contrest + 1
        noderest = rest(contrest,1)
      ELSE
        numdof = numdof + nodof
      END IF
    END DO

    !-----------------------------------------------------------------------
    ! 2. Populate nf_pp with the global degrees of freedom 
    !-----------------------------------------------------------------------

    nf_pp = 0
    k = 0
    DO i = nn_start,nn_end
      k = k + 1
      IF (i==noderest) THEN
        DO j = 2,nodof+1 
          IF(rest(contrest,j)/=1) THEN
            numdof = numdof + 1
            nf_pp(j-1,k) = numdof
          END IF
        END DO
        contrest = contrest + 1
        noderest = rest(contrest,1)
      ELSE
        DO j = 2,nodof+1 
          numdof = numdof + 1
          nf_pp(j-1,k) = numdof
        END DO
      END IF
    END DO

  END SUBROUTINE REST_NF

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE NUM_TO_G2(g_num,nf,g_g,nn_start)

    !/****f* large_strain/num_to_g2
    !*  NAME
    !*    SUBROUTINE: num_to_g2
    !*  SYNOPSIS
    !*    Usage:      CALL num_to_g2(g_num,nf,g_g,nn_start)
    !*  FUNCTION
    !*    For one element, assigns the equation number to each degree of
    !*    freedom of that element
    !*  
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    g_num(nod)          : Integer
    !*                        : Element connectivity
    !*                          
    !*    nf(nodof,nn_pp)     : Integer
    !*                        : Nodes where the equation number has been
    !*                          assigned to each the degree of freedom of the
    !*                          nodes
    !*
    !*    nn_start            : Integer
    !*                        : Lowest node number in the process
    !*                          
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    g_g(ndof)           : Integer
    !*                        : Degrees of freedom sorted in the element
    !*  AUTHOR
    !*    *** ****
    !*  CREATION DATE
    !*    **.**.****
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2008
    !*        Available under commercial licence
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  What follows is to understand the routine using 3D 8-node hexahedra
    !*  with displacements as unknows:
    !*  
    !*  What the loop does: (inode is different for each i)
    !*   i=1 -----> k=3
    !*      g(1) = nf(1,inode)
    !*      g(2) = nf(2,inode)
    !*      g(3) = nf(3,inode)
    !*   i=2 -----> k=6
    !*      g(4) = nf(1,inode)
    !*      g(5) = nf(2,inode)
    !*      g(6) = nf(3,inode)
    !*   ...
    !*   etc...
    !*   ...
    !*   i=8 -----> k=24
    !*      g(22) = nf(1,inode)
    !*      g(23) = nf(2,inode)
    !*      g(24) = nf(3,inode)
    !*
    !*  THIS SUBROUTINE IS ONLY USED IN PROGRAM XX7 AND NEEDS TO BE REMOVED
    !*  FROM PARAFEM
    !*
    !*/
        
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: g_num(:), nf(:,:), nn_start 
    INTEGER, INTENT(OUT) :: g_g(:)
    INTEGER :: i, k, nod, nodof, inode

    nod   = UBOUND(g_num,1)
    nodof = UBOUND(nf,1)

    DO i = 1,nod
      inode = g_num(i) - nn_start + 1
      k = i*nodof
      g_g(k-nodof+1:k) = nf(:,inode)
    END DO

    RETURN

  END SUBROUTINE NUM_TO_G2   
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE CALC_NEQ(nn,rest,neq)

    !/****f* large_strain/calc_neq
    !*  NAME
    !*    SUBROUTINE: calc_neq
    !*  SYNOPSIS
    !*    Usage:      CALL calc_neq(nn,rest,neq)
    !*  FUNCTION
    !*    Compute the number of equations of the problem
    !*  
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    nn                  : Integer
    !*                        : Total number of nodes
    !*
    !*    rest(nr,nodof+1)    : Integer
    !*                        : List of the nodes with some degree of freedom
    !*                          fixed. In the input, the degrees of freedom
    !*                          fixed have 1 and the degrees of freedom not
    !*                          fixed have 0. In the output, the degrees of
    !*                          freedom fixed have 0 and the degrees of freedom
    !*                          not fixed have the global equation number
    !*
    !*    The following arguments have the INTENT(OUT) attribute:    
    !*
    !*    neq                 : Integer
    !*                        : Total number of equations
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    05.03.2008
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: nn, rest(:,:)
    INTEGER, INTENT(OUT) :: neq
    INTEGER :: i, j, nr, nodof
   
    nr    = UBOUND(rest,1)
    nodof = UBOUND(rest,2) - 1

    neq = 0
    IF (nr>0) THEN
      DO i = 1,nr
        DO j = 1,nodof
          IF (rest(i,j+1)==0) THEN
            neq = neq + 1
          END IF
        END DO
      END DO
    END IF
    neq = (nn-nr)*nodof + neq

  END SUBROUTINE CALC_NEQ

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE COMPUTE_NPES_PP(nels,neq,nn,npes,numpe,g_num_pp,rest,npes_pp)

    !/****f* large_strain/compute_npes_pp
    !*  NAME
    !*    SUBROUTINE: compute_npes_pp
    !*  SYNOPSIS
    !*    Usage:      CALL compute_npes_pp(nels,neq,nn,npes,numpe,g_num_pp, &
    !*                                     rest,npes_pp)
    !*  FUNCTION
    !*    Compute the integer npes_pp for processes comunication regarding
    !*    the distribution of equations across processes
    !*  
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    nels                 : Integer
    !*                           Total number of elements
    !*
    !*    neq                  : Integer
    !*                           Total number of equations
    !*
    !*    nn                   : Integer
    !*                           Total number of nodes
    !*
    !*    npes                 : Integer
    !*                           Number of processes
    !*
    !*    numpe                : Integer
    !*                           Process number
    !*
    !*    g_num_pp(nod,nels_pp): Integer
    !*                           Nodal connectivity
    !*
    !*    rest(nr,nodof+1)     : Integer
    !*                           Constrained degrees of freedom of nodes
    !*
    !*    The following array argument has the INTENT(OUT) attribute: 
    !*
    !*    npes_pp        : Integer
    !*                     Maximum number of processors for communications
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    05.03.2008
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  Based on master process collecting and doing the work
    !*
    !*  THIS SUBROUTINE IS ONLY USED BY PROGRAM XX7 AND NEEDS TO BE REMOVED
    !*  FROM PARAFEM
    !*/    

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: nels, neq, nn, npes, numpe, g_num_pp(:,:), rest(:,:)
    INTEGER, INTENT(OUT) :: npes_pp
    INTEGER :: i, j, k, inod, iel, nels_pp, nod, nr, nodof, bufsize, ier, &
               statu(MPI_STATUS_SIZE), bufsize1, bufsize2, countbound,    &
               countequa, actualproc, numproc, maxconec, jel
    INTEGER, ALLOCATABLE :: nodesnum_pp(:), nodesnum(:), counting(:), &
                            num_elem_pp(:), cone(:,:), equaproc(:),   &
                            equanodes(:), nodes01(:,:), vect12(:),    &
                            nodesconec(:,:), elemproc(:), elemacu(:)

    nod     = UBOUND(g_num_pp,1)
    nels_pp = UBOUND(g_num_pp,2)
    nr      = UBOUND(rest,1)
    nodof   = UBOUND(rest,2) - 1

    !----------------------------------------------------------------------
    ! 1. Compute the maximum number of elements shared by a node
    !----------------------------------------------------------------------

    ALLOCATE(nodesnum_pp(nn))
    IF (numpe==1) THEN
      ALLOCATE(nodesnum(nn))
      nodesnum    = 0
    END IF

    nodesnum_pp = 0
    DO iel = 1,nels_pp
      nodesnum_pp(g_num_pp(:,iel)) = nodesnum_pp(g_num_pp(:,iel)) + 1
    END DO

    bufsize = nn
    CALL MPI_REDUCE(nodesnum_pp,nodesnum,bufsize,MPI_INTEGER,MPI_SUM,0, &
                       MPI_COMM_WORLD,ier)

    IF(numpe==1) maxconec = MAXVAL(nodesnum(:))

    DEALLOCATE(nodesnum_pp)

    !----------------------------------------------------------------------
    ! 2. Master process populates array "nodesconec", containing for each node
    !    the list of elements that contain the node (element connectivity)
    !----------------------------------------------------------------------
    
    IF (numpe==1) THEN
      ALLOCATE(nodesconec(nn,maxconec), counting(nn), cone(nod,nels_pp))
      ALLOCATE(num_elem_pp(npes), elemacu(npes+1))
    END IF

    !----------------------------------------------------------------------
    ! 2.1 Master process receives the integer nels_pp from slave processes
    !----------------------------------------------------------------------

    IF (numpe==1) THEN
      num_elem_pp    = 0
      num_elem_pp(1) = nels_pp
      elemacu(1)     = 0
      elemacu(2)     = nels_pp
    END IF
    bufsize = 1

    DO i = 2,npes
      IF(numpe==i) THEN
        CALL MPI_SEND(nels_pp,bufsize,MPI_INTEGER,0,i,MPI_COMM_WORLD,ier)
      END IF
      IF(numpe==1) THEN
        CALL MPI_RECV(j,bufsize,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
        num_elem_pp(i) = j
        elemacu(i+1) = elemacu(i) + j
      END IF
    END DO
    
    !----------------------------------------------------------------------
    ! 2.2 Master process goes through its own elements
    !----------------------------------------------------------------------

    IF (numpe==1) THEN
      counting   = 0
      nodesconec = 0
      DO iel = 1,nels_pp
        counting(g_num_pp(:,iel)) = counting(g_num_pp(:,iel)) + 1
        DO inod = 1,nod
          nodesconec(g_num_pp(inod,iel),counting(g_num_pp(inod,iel))) = iel
        END DO
      END DO
    END IF

    !----------------------------------------------------------------------
    ! 2.3 Master process receives from slave processes the connectivity array
    !     and goes through slave processes' elements
    !----------------------------------------------------------------------
    DO i = 2,npes
      IF (numpe==i) THEN
        bufsize1 = nels_pp*nod
        CALL MPI_SEND(g_num_pp,bufsize1,MPI_INTEGER,0,i,MPI_COMM_WORLD,ier)
      END IF
      IF (numpe==1) THEN
        bufsize2 = num_elem_pp(i)*nod
        CALL MPI_RECV(cone,bufsize2,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier) 
        jel = 0
        DO iel = elemacu(i)+1,elemacu(i+1)
          jel = jel + 1
          counting(cone(:,jel)) = counting(cone(:,jel)) + 1
          DO inod = 1,nod
            nodesconec(cone(inod,jel),counting(cone(inod,jel))) = iel
          END DO
        END DO
      END IF
    END DO

    IF (numpe==1) THEN
      DEALLOCATE(counting, cone)
    END IF

    !----------------------------------------------------------------------
    ! 3. Master process populates "elemproc"/"equaproc", assigning to each
    !    element/equation  its process number 
    !----------------------------------------------------------------------

    IF (numpe==1) THEN
      ALLOCATE(elemproc(nels), equaproc(neq))

      CALL CALC_ELEMPROC(nels,npes,elemproc)
      CALL CALC_ELEMPROC(neq,npes,equaproc)
    END IF

    !----------------------------------------------------------------------
    ! 4. Master process populates "nodes01", assigning to each pair 
    !    (node,processor) 1 if there is some element in that processor
    !    containing the node and 0 if not.
    !----------------------------------------------------------------------
 
    IF (numpe==1) THEN
      ALLOCATE(nodes01(nn,npes))

      nodes01 = 0
      DO inod = 1,nn
        DO j = 1,nodesnum(inod)
          nodes01(inod,elemproc(nodesconec(inod,j))) = 1
        END DO
      END DO

      DEALLOCATE(nodesconec, nodesnum, elemproc)
    END IF

    !----------------------------------------------------------------------
    ! 5. Master process populates array "equanodes", assigning to each 
    !    equation number the node producing the equation
    !----------------------------------------------------------------------

    IF (numpe==1) THEN
      ALLOCATE(equanodes(neq))

      IF (nr>0) THEN
        countbound = 1
      END IF
      countequa = 0
      DO i = 1,nn
        IF (i==rest(countbound,1)) THEN
          DO j = 1,nodof
            IF (rest(countbound,j+1)==0) THEN
              countequa = countequa + 1
              equanodes(countequa) = i
            END IF
          END DO
          IF (countbound < nr) THEN
            countbound = countbound + 1
          END IF
        ELSE
          DO j = 1,nodof
            countequa = countequa + 1
            equanodes(countequa) = i
          END DO
        END IF
      END DO
    END IF

    !----------------------------------------------------------------------
    ! 6. Master process computes "npes_pp" by a loop through the number of
    !    equations
    !----------------------------------------------------------------------

    IF (numpe==1) THEN
      ALLOCATE(vect12(npes))

      npes_pp = 0
      actualproc = 1
      vect12 = 0

      DO i = 1,neq
        IF (equaproc(i) /= actualproc) THEN
          numproc = 0
          DO j = 1,npes
            IF (vect12(j)>0) THEN 
              numproc = numproc + 1
            END IF
          END DO
          npes_pp = MAX(npes_pp,numproc)
          actualproc = actualproc + 1
          vect12 = 0
        END IF
        DO inod = 1,nod
          vect12(:) = vect12(:) + nodes01(equanodes(i),:)
        END DO
      END DO

      numproc = 0
      DO j = 1,npes
        IF (vect12(j)>0) THEN 
          numproc = numproc + 1
        END IF
      END DO
      npes_pp = MAX(npes_pp,numproc)

      DEALLOCATE(equaproc, equanodes, nodes01, vect12)

    END IF

    !----------------------------------------------------------------------
    ! 7. Master process broadcasts "npes_pp" to slave processes
    !----------------------------------------------------------------------

    bufsize = 1
    CALL MPI_BCAST(npes_pp,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  END SUBROUTINE COMPUTE_NPES_PP

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE PCG_VER1(inewton,limit,tol,A_pp,b_pp,Minv_pp,rn0,x_pp,iters)

    !/****f* large_strain/pcg_ver1
    !*  NAME
    !*    SUBROUTINE: pcg_ver1
    !*  SYNOPSIS
    !*    Usage:      CALL pcg_ver1(inewton,limit,tol,A_pp,b_pp,Minv_pp, &
    !*                              rn0,x_pp,iters)
    !*  FUNCTION
    !*    Iterative solver (PCG) to solve a linear system of equations
    !*    Ver1: The residual r is always computed approximately as
    !*          r = r - alpha*q
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    inewton                  : Integer
    !*                             : Newton-Raphson iteration number
    !*
    !*    limit                    : Integer
    !*                             : Maximum number of PCG iterations allowed
    !*
    !*    tol                      : Real
    !*                             : Tolerance for PCG
    !*
    !*    A_pp(ntot,ntot,nels_pp)) : Real
    !*                             : Element stiffness matrices
    !*
    !*    b_pp(neq_pp)             : Real
    !*                             : Residual of the equations
    !*
    !*    Minv_pp(neq_pp)          : Real
    !*                             : Inverse of diagonal preconditioner
    !*
    !*    The following arguments have the INTENT(INOUT) attribute:
    !*
    !*    rn0                      : Real
    !*                             : Square norm of the residual in the first
    !*                               Newton-Raphson iteration
    !*
    !*    x_pp(neq_pp)             : Real
    !*                             : Solution of the system of equations
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    iters                    : Integer
    !*                             : Number of PCG iterations
    !*
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    01.06.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  The solution scheme is based on x = 0.0 as starting guess and 
    !*  is as follows:
    !*
    !*  r = b - A*x  (r = b  if x = 0.0 as starting guess)
    !*  d = Minv*r
    !*  delta_new = r*d
    !*  iters = 0
    !*  WHILE (iters<limit and sqrt(rn/rn0) < tol) DO
    !*    iters = iters + 1
    !*    q = A*d
    !*    alpha = delta_new/(d*q)
    !*    x = x + alpha*d
    !*    r = r - alpha*q
    !*    s = Minv*r
    !*    delta_old = delta_new
    !*    delta_new = r*s
    !*    beta = delta_new/delta_old
    !*    d = s + beta*d
    !*  END DO  
    !* 
    !*  THIS SUBROUTINE IS ONLY USED IN PROGRAM XX7 AND NEEDS TO BE REMOVED
    !*  FROM PARAFEM
    !*/

    IMPLICIT NONE

    INTEGER,   INTENT(IN)    :: inewton, limit
    REAL(iwp), INTENT(IN)    :: tol, b_pp(:), A_pp(:,:,:), Minv_pp(:)
    REAL(iwp), INTENT(INOUT) :: rn0, x_pp(:)
    INTEGER,   INTENT(OUT)   :: iters
    LOGICAL :: converged
    INTEGER :: nels_pp, neq_pp, ntot, iel
    REAL(iwp) :: alpha, beta, delta_old, delta_new
    REAL(iwp) :: tolr0, rn
    REAL(iwp), ALLOCATABLE :: d_pp(:), dmul_pp(:,:), qtemp_pp(:,:),  &
                              q_pp(:), s_pp(:), r_pp(:)

    !-----------------------------------------------------------------------
    ! 1. Allocate internal arrays
    !-----------------------------------------------------------------------

    neq_pp  = UBOUND(b_pp,1)
    ntot    = UBOUND(A_pp,1)
    nels_pp = UBOUND(A_pp,3)

    ALLOCATE (d_pp(neq_pp), dmul_pp(ntot,nels_pp), qtemp_pp(ntot,nels_pp),   &
              q_pp(neq_pp), s_pp(neq_pp), r_pp(neq_pp))

    !-----------------------------------------------------------------------
    ! 2. Compute the norm of the residual in the first Newto-Raphson iteration
    !    to compare with later iterations
    !-----------------------------------------------------------------------

    r_pp(:) = b_pp(:)

    IF (inewton==1) THEN
      rn0 = DOT_PRODUCT_P(r_pp,r_pp)
    END IF
    tolr0 = rn0*tol**2
	
    !-----------------------------------------------------------------------
    ! 3. Initialise scalars and arrays
    !-----------------------------------------------------------------------

    d_pp = Minv_pp * r_pp 
    delta_new = DOT_PRODUCT_P(r_pp,d_pp)

    iters     = 0
    converged = .FALSE.

    !-----------------------------------------------------------------------
    ! 4. PCG algorithm
    !-----------------------------------------------------------------------

    iterations  :  DO
      iters = iters + 1

      dmul_pp = .0_iwp
      CALL GATHER(d_pp,dmul_pp)
      DO iel = 1, nels_pp
        qtemp_pp(:,iel) = MATMUL(A_pp(:,:,iel),dmul_pp(:,iel))
      END DO
      q_pp = 0._iwp
      CALL SCATTER(q_pp,qtemp_pp)

      alpha = delta_new/DOT_PRODUCT_P(d_pp,q_pp)
      x_pp = x_pp + alpha*d_pp
      r_pp = r_pp - alpha*q_pp
      s_pp = Minv_pp*r_pp
      delta_old = delta_new
      delta_new = DOT_PRODUCT_P(r_pp,s_pp)
      beta = delta_new/delta_old
      d_pp = s_pp + beta*d_pp

      !---------------------------------------------------------------------
      ! 4.1 Check convergence. Master process reports results
      !---------------------------------------------------------------------

      rn = DOT_PRODUCT_P(r_pp,r_pp)
  
      IF (numpe==1) THEN
        WRITE(90,'(i7,3(1p,e17.7))')iters,SQRT(rn0),SQRT(rn),SQRT(rn/rn0)
!       CALL FLUSH(90)
      END IF

      IF (rn<tolr0) THEN
        converged = .TRUE.
      END IF
 
      IF (converged .OR. iters==limit) THEN
        EXIT
      END IF

    END DO iterations

    !-----------------------------------------------------------------------
    ! 5. Deallocate local arrays
    !-----------------------------------------------------------------------

    DEALLOCATE (d_pp,dmul_pp,qtemp_pp,q_pp,s_pp,r_pp)

    RETURN

  END SUBROUTINE PCG_VER1

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE SHAPE_FUNCTIONS(i,points,fun)

    !/****f* large_strain/shape_functions
    !*  NAME
    !*    SUBROUTINE: shape_functions
    !*  SYNOPSIS
    !*    Usage:      CALL shape_functions(i,points,fun)
    !*  FUNCTION
    !*    Compute the value of the shape functions at a Gauss point
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    i                  : Integer
    !*                       : Gauss point number
    !*
    !*    points(ndim,nip)   : Real
    !*                       : Gauss points coordinates at the reference 
    !*                         element
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    fun(nod)           : Real
    !*                       : Value of the shape functions at a 
    !*                         Gauss point
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    25.01.2008
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  This routine has to be completed with elements for 1 and 2 dimensions  
    !*  and 3D elements with unusual number of nodes
    !*   
    !*  THIS SUBROUTINE IS ONLY USED IN PROGRAM XX7 AND NEEDS TO BE REMOVED
    !*  FROM PARAFEM
    !*/
     
    INTEGER,   INTENT(IN)  :: i
    REAL(iwp), INTENT(IN)  :: points(:,:)
    REAL(iwp), INTENT(OUT) :: fun(:) 
    INTEGER :: ndim, nod
    REAL(iwp) :: xi, xip, xim, eta, etap, etam, zeta, zetap, zetam 

    ndim = UBOUND(points,1)
    nod  = UBOUND(fun,1)
    
    SELECT CASE(ndim)
      CASE(3) !three dimensional elements
        xi    = points(1,i)
        eta   = points(2,i)
        zeta  = points(3,i)
        xip   = 1._iwp + xi 
        etap  = 1._iwp + eta 
        zetap = 1._iwp + zeta
        xim   = 1._iwp - xi 
        etam  = 1._iwp - eta 
        zetam = 1._iwp - zeta
        SELECT CASE(nod)	 
          CASE(4) !4-node tetrahedra
            fun(1) = xi                       !N1
            fun(2) = eta                      !N2
            fun(3) = 1._iwp - xi - eta - zeta !N3
            fun(4) = zeta                     !N4
          CASE(8) !8-node hexahedra
            fun(1) = 0.125*xip*etam*zetam     !N1
            fun(2) = 0.125*xip*etap*zetam     !N2
            fun(3) = 0.125*xim*etap*zetam     !N3
            fun(4) = 0.125*xim*etam*zetam     !N4
            fun(5) = 0.125*xip*etam*zetap     !N5
            fun(6) = 0.125*xip*etap*zetap     !N6
            fun(7) = 0.125*xim*etap*zetap     !N7
            fun(8) = 0.125*xim*etam*zetap     !N8
          CASE DEFAULT
            WRITE(*,*)"Wrong number of nodes!!"
        END SELECT
      CASE DEFAULT
        WRITE(*,*)"Wrong number of dimensions!!"
    END SELECT

    RETURN

  END SUBROUTINE SHAPE_FUNCTIONS
  SUBROUTINE PLASTICITY(deeF,jacF,jacFinc,sigma1C,statev,lnstrainelas, &
   sigma,detF,statevar_num,iel,igauss)
    ! Updates the stresses and computes the spatial tangent operator from the 
    ! material contribution in the following form: Aijkl=(1/2J)*(D:L:B)ijkl
    ! D is the derivative of the Kirchhoff stress with respect to the 
    ! logarithmic strain (this format is the format of the tangent operator 
    ! in a UMAT, and this is the reason why it is adopted)
    ! L is the derivative of the logarithmic strain with respect to B, obtained
    ! in the SUBROUTINE LN_DERIV
    ! B is the derivative of B with respect to F, obtained in the SUBROUTINE 
    ! BDERIVF
    IMPLICIT NONE

    REAL(iwp), INTENT(IN) :: jacF(:,:), jacFinc(:,:), detF
    INTEGER, INTENT(IN) :: statevar_num
    REAL(iwp), INTENT(OUT) :: deeF(:,:), sigma(:,:), sigma1C(:)
    REAL(iwp), INTENT(INOUT) :: lnstrainelas(:), statev(:)
    REAL(iwp) :: d_tensor(9,9), lnderiv(9,9), strderb(9,9), bderiv(9,9),      &
     geom_comp(9,9), ddsdde(6,6), b_tensor_3x3(3,3), et1(6), et2(6), et3(6),  &
     b_tensor(6), dstran(6), btensor(6), lne1, lne2, lne3, e1, e2, e3
    REAL(iwp), PARAMETER :: half=0.5_iwp, one=1._iwp, two=2._iwp,             &
     tol=0.000001_iwp
    INTEGER :: ntens, iel, igauss
       
    ntens=6

    ! Calculate the previously converged elastic B strain tensor
    lnstrainelas=two*lnstrainelas
    
    CALL exp_tensor(lnstrainelas,b_tensor_3x3,e1,e2,e3,et1,et2,et3,tol)
    
    ! Calculate the trial elastic B strain tensor
    b_tensor_3x3=MATMUL(MATMUL(jacFinc,b_tensor_3x3),TRANSPOSE(jacFinc))
    
    b_tensor(1)=b_tensor_3x3(1,1)
    b_tensor(2)=b_tensor_3x3(2,2)
    b_tensor(3)=b_tensor_3x3(3,3)
    b_tensor(4)=b_tensor_3x3(1,2)
    b_tensor(5)=b_tensor_3x3(2,3)
    b_tensor(6)=b_tensor_3x3(1,3)
    
    CALL eigen(e1,e2,e3,b_tensor,et1,et2,et3,tol)

    ! Calculate logarithmic strain through the B tensor
    CALL ln_strain(lne1,lne2,lne3,e1,e2,e3,et1,et2,et3,lnstrainelas,tol)
    
    ! Assign the values so that the umat is as similar as possible to ABAQUS 
    ! umat
    lnstrainelas(4)=two*lnstrainelas(4)
    lnstrainelas(5)=two*lnstrainelas(5)
    lnstrainelas(6)=two*lnstrainelas(6)

    ! Calculate the Cauchy stress tensor and the consistent tangent operator 
    ! (material contribution)
    ! Calling UMAT, which consists in the same format as the one used in
    ! ABAQUS. The only difference is that the material properties should be
    ! defined inside of the subroutine
    sigma1C=0._iwp
    CALL umat(sigma1C,statev,ddsdde,lnstrainelas,dstran,ntens,statevar_num,iel,igauss)

    ! Convert Kirchhoff stress into Cauchy stress
    sigma1C=(one/detF)*sigma1C
    
    ! This is the updated elastic strain
    lnstrainelas(4)=half*lnstrainelas(4)
    lnstrainelas(5)=half*lnstrainelas(5)
    lnstrainelas(6)=half*lnstrainelas(6)

    ! This is the second order tensor of updated Cauchy stress
    sigma(1,1)=sigma1C(1)
    sigma(2,2)=sigma1C(2)
    sigma(3,3)=sigma1C(3)
    sigma(1,2)=sigma1C(4)
    sigma(2,3)=sigma1C(5)
    sigma(1,3)=sigma1C(6)
    sigma(2,1)=sigma(1,2)
    sigma(3,1)=sigma(1,3)
    sigma(3,2)=sigma(2,3)

    ! Transforms the D matrix from ABAQUS notation to a complete fourth-order
    ! tensor
    CALL abq_to_parafem(ddsdde,d_tensor)
    
    ! Computes the derivative of the logarithmic strain with respect to the 
    ! Left Cauchy-Green deformation tensor
    CALL ln_deriv(lnderiv,e1,e2,e3,lne1,lne2,lne3,et1,et2,et3,tol)
    
    ! Computes the derivative of the Left Cauchy-Green deformation tensor with
    ! respect to the deformation gradient
    CALL bderivf(b_tensor,bderiv)

    strderb=MATMUL(d_tensor,lnderiv)
    d_tensor=MATMUL(strderb,bderiv)
    
    ! Calculate the geometric component
    CALL geom_component(sigma1C,geom_comp)
    
    ! Assemble the tangent operator
    deeF=(d_tensor*(half/detF))+geom_comp
    
    ! Enforce symmetry (correcting minor errors)
    deeF=(deeF+TRANSPOSE(deeF))*half

  RETURN
  END SUBROUTINE PLASTICITY

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE EXP_TENSOR(tensor,tensor_exp_3x3,e1,e2,e3,et1,et2,et3,tol)
    ! Calculates the tensor exponential from a symmetric tensor argument
    IMPLICIT NONE

    REAL(iwp), INTENT(IN) :: tensor(:), tol
    REAL(iwp), INTENT(OUT) :: tensor_exp_3x3(:,:), et1(:), et2(:), et3(:),e1, &
     e2, e3
    REAL(iwp) :: tensor_exp(6)

    CALL eigen(e1,e2,e3,tensor,et1,et2,et3,tol)

    tensor_exp=DEXP(e1)*et1+DEXP(e2)*et2+DEXP(e3)*et3

    tensor_exp_3x3(1,1)=tensor_exp(1)
    tensor_exp_3x3(2,2)=tensor_exp(2)
    tensor_exp_3x3(3,3)=tensor_exp(3)
    tensor_exp_3x3(1,2)=tensor_exp(4)
    tensor_exp_3x3(2,3)=tensor_exp(5)
    tensor_exp_3x3(1,3)=tensor_exp(6)
    tensor_exp_3x3(2,1)=tensor_exp_3x3(1,2)
    tensor_exp_3x3(3,2)=tensor_exp_3x3(2,3)
    tensor_exp_3x3(3,1)=tensor_exp_3x3(1,3)

  RETURN
  END SUBROUTINE EXP_TENSOR

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE EIGEN(e1,e2,e3,tensor,et1,et2,et3,tol)
    ! This subroutine calculates the eigenvalues and eigenprojections of a    
    ! second order tensor
    IMPLICIT NONE

    REAL(iwp), INTENT(IN) :: tensor(:), tol
    REAL(iwp), INTENT(OUT) :: et1(:), et2(:), et3(:), e1, e2, e3
    REAL(iwp) :: tensor2(6), unit_tensor(6), scaled_tensor(6),                &
     scaled_tensor2(6), inv1, inv2, inv3, r, q, alpha, pi, tr2, sinv1, sinv2, &
     sinv3, str2, sca_tol
    REAL(iwp), PARAMETER :: zero=0._iwp, half=0.5_iwp, one=1._iwp, two=2._iwp,&
     three=3._iwp, nine=9._iwp, cons=5._iwp
    INTEGER :: i

    ! Scale the tolerance
    sca_tol=cons*tol
    
    ! Initialize the values of the eigenprojections
    et1=zero
    et2=zero
    et3=zero
    
    ! Create the squared tensor
    tensor2(1)=(tensor(1)**2)+(tensor(4)**2)+(tensor(6)**2)
    tensor2(2)=(tensor(4)**2)+(tensor(2)**2)+(tensor(5)**2)
    tensor2(3)=(tensor(6)**2)+(tensor(5)**2)+(tensor(3)**2)
    tensor2(4)=tensor(1)*tensor(4)+tensor(4)*tensor(2)+tensor(6)*tensor(5)
    tensor2(5)=tensor(4)*tensor(6)+tensor(2)*tensor(5)+tensor(5)*tensor(3)
    tensor2(6)=tensor(1)*tensor(6)+tensor(4)*tensor(5)+tensor(6)*tensor(3)
     
    ! Initialize the scaled tensor
    scaled_tensor=cons*tensor

    ! Create the scaled squared tensor
    scaled_tensor2(1)=(scaled_tensor(1)**2)+(scaled_tensor(4)**2)+                   &
     (scaled_tensor(6)**2)
    scaled_tensor2(2)=(scaled_tensor(4)**2)+(scaled_tensor(2)**2)+                   &
     (scaled_tensor(5)**2)
    scaled_tensor2(3)=(scaled_tensor(6)**2)+(scaled_tensor(5)**2)+                   &
     (scaled_tensor(3)**2)
    scaled_tensor2(4)=scaled_tensor(1)*scaled_tensor(4)+scaled_tensor(4)*            &
     scaled_tensor(2)+scaled_tensor(6)*scaled_tensor(5)
    scaled_tensor2(5)=scaled_tensor(4)*scaled_tensor(6)+scaled_tensor(2)*            &
     scaled_tensor(5)+scaled_tensor(5)*scaled_tensor(3)
    scaled_tensor2(6)=scaled_tensor(1)*scaled_tensor(6)+scaled_tensor(4)*            &
     scaled_tensor(5)+scaled_tensor(6)*scaled_tensor(3)
    
    ! Calculate the invariants of the tensor
    inv1=tensor(1)+tensor(2)+tensor(3)

    tr2=tensor2(1)+tensor2(2)+tensor2(3)
    
    inv2=half*((inv1**2)-tr2)
    
    inv3=(tensor(1)*tensor(2)*tensor(3)+two*tensor(4)*tensor(5)*tensor(6))-                     &
     ((tensor(6)**2)*tensor(2)+(tensor(5)**2)*tensor(1)+(tensor(4)**2)*tensor(3))
    
    ! Calculate the invariants of the scaled tensor
    sinv1=scaled_tensor(1)+scaled_tensor(2)+scaled_tensor(3)

    str2=scaled_tensor2(1)+scaled_tensor2(2)+scaled_tensor2(3)
    
    sinv2=half*((sinv1**2)-str2)
    
    sinv3=(scaled_tensor(1)*scaled_tensor(2)*scaled_tensor(3)+two*             &
     scaled_tensor(4)*scaled_tensor(5)*scaled_tensor(6))-                     &
     ((scaled_tensor(6)**2)*scaled_tensor(2)+(scaled_tensor(5)**2)*           &
     scaled_tensor(1)+(scaled_tensor(4)**2)*scaled_tensor(3))
    
    ! Calculate the eigenvalues
    r=(-two*(sinv1**3)+nine*sinv1*sinv2-27._iwp*sinv3)/54._iwp
    q=((sinv1**2)-three*sinv2)/nine
    
    IF (q<zero) THEN
      q=zero
    END IF

    alpha=r/DSQRT(q**3)

    ! Check for other extreme cases for q
    IF (alpha<-one) THEN
      alpha=-one
    ELSEIF (alpha>one) THEN
      alpha=one
!   ELSEIF (ISNAN(alpha)) THEN
!     alpha=one
!   needs replacing as ISNAN is not in all releases of Fortran
    END IF
    
    alpha=DACOS(alpha)
    pi=4._iwp*DATAN(one)

    e1=-two*DSQRT(q)*DCOS(alpha/three)+sinv1/three
    e2=-two*DSQRT(q)*DCOS((alpha+two*pi)/three)+sinv1/three
    e3=-two*DSQRT(q)*DCOS((alpha-two*pi)/three)+sinv1/three
    
    ! Calculate the eigenprojections
    ! Assign the value to the unit_tensor
    DO i=1,3
      unit_tensor(i)=one
      unit_tensor(i+3)=zero
    END DO

    ! If all eigenvalues are different
    IF ((DABS(e1-e2)>sca_tol).and.(DABS(e2-e3)>sca_tol).and.(DABS(e1-e3)>     &
     sca_tol)) THEN
      et1=(e1/(two*(e1**3)-sinv1*(e1**2)+sinv3))*(scaled_tensor2-(sinv1-e1)*        &
       scaled_tensor+(sinv3/e1)*unit_tensor)
      et2=(e2/(two*(e2**3)-sinv1*(e2**2)+sinv3))*(scaled_tensor2-(sinv1-e2)*        &
       scaled_tensor+(sinv3/e2)*unit_tensor)
      et3=(e3/(two*(e3**3)-sinv1*(e3**2)+sinv3))*(scaled_tensor2-(sinv1-e3)*        &
       scaled_tensor+(sinv3/e3)*unit_tensor)

    ! If two eigenvalues are equal
    ! If e1 and e2 are equal
    ELSEIF ((DABS(e1-e2)<sca_tol).and.(DABS(e2-e3)>sca_tol).and.(DABS(e1-e3)> &
     sca_tol)) THEN
      et3=(e3/(two*(e3**3)-sinv1*(e3**2)+sinv3))*(scaled_tensor2-(sinv1-e3)*        &
       scaled_tensor+(sinv3/e3)*unit_tensor)
      et1=unit_tensor-et3

    ! If e1 and e3 are equal
    ELSEIF ((DABS(e1-e3)<sca_tol).and.(DABS(e1-e2)>sca_tol).and.(DABS(e2-e3)>sca_tol)) THEN
      et2=(e2/(two*(e2**3)-sinv1*(e2**2)+sinv3))*(scaled_tensor2-(sinv1-e2)*        &
       scaled_tensor+(sinv3/e2)*unit_tensor)
      et1=unit_tensor-et2

    ! If e2 and e3 are equal
    ELSEIF ((DABS(e2-e3)<sca_tol).and.(DABS(e1-e2)>sca_tol).and.(DABS(e1-e3)>sca_tol)) THEN
      et1=(e1/(two*(e1**3)-sinv1*(e1**2)+sinv3))*(scaled_tensor2-(sinv1-e1)*        &
       scaled_tensor+(sinv3/e1)*unit_tensor)
      et2=unit_tensor-et1

    ! If all eigenvalues are equal
    ELSE
      et1=unit_tensor

    END IF

    e1=e1/cons
    e2=e2/cons
    e3=e3/cons
    
  RETURN
  END SUBROUTINE EIGEN

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE LN_STRAIN(lne1,lne2,lne3,e1,e2,e3,et1,et2,et3,lnstrainelas,tol)
    ! Computes the Left Cauchy-Green deformation tensor and then the          
    ! logarithmic strain tensor eigenvalues
    IMPLICIT NONE

    REAL(iwp), INTENT(IN) :: et1(:), et2(:), et3(:), e1, e2, e3, tol
    REAL(iwp), INTENT(OUT) :: lnstrainelas(:), lne1, lne2, lne3

    lne1=DLOG(e1)
    lne2=DLOG(e2)
    lne3=DLOG(e3)

    lnstrainelas=0.5_iwp*(lne1*et1+lne2*et2+lne3*et3)

  RETURN
  END SUBROUTINE LN_STRAIN

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  
  SUBROUTINE ABQ_TO_PARAFEM(tensor_6x6,tensor_9x9)
    ! Transforms the infinitesimal strain tangent operator from ABAQUS Voigt 
    ! notation to a full fourth order tangent operator so that the tensorial 
    ! operations can be carried out
    IMPLICIT NONE

    REAL(iwp), INTENT(IN) :: tensor_6x6(:,:)
    REAL(iwp), INTENT(OUT) :: tensor_9x9(:,:)

    tensor_9x9(1,1)=tensor_6x6(1,1)
    tensor_9x9(1,2)=tensor_6x6(1,4)
    tensor_9x9(1,3)=tensor_6x6(1,6)
    tensor_9x9(1,4)=tensor_6x6(1,4)
    tensor_9x9(1,5)=tensor_6x6(1,2)
    tensor_9x9(1,6)=tensor_6x6(1,5)
    tensor_9x9(1,7)=tensor_6x6(1,6)
    tensor_9x9(1,8)=tensor_6x6(1,5)
    tensor_9x9(1,9)=tensor_6x6(1,3)
    
    tensor_9x9(2,1)=tensor_6x6(4,1)
    tensor_9x9(2,2)=tensor_6x6(4,4)
    tensor_9x9(2,3)=tensor_6x6(4,6)
    tensor_9x9(2,4)=tensor_6x6(4,4)
    tensor_9x9(2,5)=tensor_6x6(4,2)
    tensor_9x9(2,6)=tensor_6x6(4,5)
    tensor_9x9(2,7)=tensor_6x6(4,6)
    tensor_9x9(2,8)=tensor_6x6(4,5)
    tensor_9x9(2,9)=tensor_6x6(4,3)
    
    tensor_9x9(3,1)=tensor_6x6(6,1)
    tensor_9x9(3,2)=tensor_6x6(6,4)
    tensor_9x9(3,3)=tensor_6x6(6,6)
    tensor_9x9(3,4)=tensor_6x6(6,4)
    tensor_9x9(3,5)=tensor_6x6(6,2)
    tensor_9x9(3,6)=tensor_6x6(6,5)
    tensor_9x9(3,7)=tensor_6x6(6,6)
    tensor_9x9(3,8)=tensor_6x6(6,5)
    tensor_9x9(3,9)=tensor_6x6(6,3)
    
    tensor_9x9(4,1)=tensor_6x6(4,1)
    tensor_9x9(4,2)=tensor_6x6(4,4)
    tensor_9x9(4,3)=tensor_6x6(4,6)
    tensor_9x9(4,4)=tensor_6x6(4,4)
    tensor_9x9(4,5)=tensor_6x6(4,2)
    tensor_9x9(4,6)=tensor_6x6(4,5)
    tensor_9x9(4,7)=tensor_6x6(4,6)
    tensor_9x9(4,8)=tensor_6x6(4,5)
    tensor_9x9(4,9)=tensor_6x6(4,3)
    
    tensor_9x9(5,1)=tensor_6x6(2,1)
    tensor_9x9(5,2)=tensor_6x6(2,4)
    tensor_9x9(5,3)=tensor_6x6(2,6)
    tensor_9x9(5,4)=tensor_6x6(2,4)
    tensor_9x9(5,5)=tensor_6x6(2,2)
    tensor_9x9(5,6)=tensor_6x6(2,5)
    tensor_9x9(5,7)=tensor_6x6(2,6)
    tensor_9x9(5,8)=tensor_6x6(2,5)
    tensor_9x9(5,9)=tensor_6x6(2,3)
    
    tensor_9x9(6,1)=tensor_6x6(5,1)
    tensor_9x9(6,2)=tensor_6x6(5,4)
    tensor_9x9(6,3)=tensor_6x6(5,6)
    tensor_9x9(6,4)=tensor_6x6(5,4)
    tensor_9x9(6,5)=tensor_6x6(5,2)
    tensor_9x9(6,6)=tensor_6x6(5,5)
    tensor_9x9(6,7)=tensor_6x6(5,6)
    tensor_9x9(6,8)=tensor_6x6(5,5)
    tensor_9x9(6,9)=tensor_6x6(5,3)
    
    tensor_9x9(7,1)=tensor_6x6(6,1)
    tensor_9x9(7,2)=tensor_6x6(6,4)
    tensor_9x9(7,3)=tensor_6x6(6,6)
    tensor_9x9(7,4)=tensor_6x6(6,4)
    tensor_9x9(7,5)=tensor_6x6(6,2)
    tensor_9x9(7,6)=tensor_6x6(6,5)
    tensor_9x9(7,7)=tensor_6x6(6,6)
    tensor_9x9(7,8)=tensor_6x6(6,5)
    tensor_9x9(7,9)=tensor_6x6(6,3)
    
    tensor_9x9(8,1)=tensor_6x6(5,1)
    tensor_9x9(8,2)=tensor_6x6(5,4)
    tensor_9x9(8,3)=tensor_6x6(5,6)
    tensor_9x9(8,4)=tensor_6x6(5,4)
    tensor_9x9(8,5)=tensor_6x6(5,2)
    tensor_9x9(8,6)=tensor_6x6(5,5)
    tensor_9x9(8,7)=tensor_6x6(5,6)
    tensor_9x9(8,8)=tensor_6x6(5,5)
    tensor_9x9(8,9)=tensor_6x6(5,3)
    
    tensor_9x9(9,1)=tensor_6x6(3,1)
    tensor_9x9(9,2)=tensor_6x6(3,4)
    tensor_9x9(9,3)=tensor_6x6(3,6)
    tensor_9x9(9,4)=tensor_6x6(3,4)
    tensor_9x9(9,5)=tensor_6x6(3,2)
    tensor_9x9(9,6)=tensor_6x6(3,5)
    tensor_9x9(9,7)=tensor_6x6(3,6)
    tensor_9x9(9,8)=tensor_6x6(3,5)
    tensor_9x9(9,9)=tensor_6x6(3,3)
    
  RETURN
  END SUBROUTINE ABQ_TO_PARAFEM

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE LN_DERIV(lnderiv,e1,e2,e3,lne1,lne2,lne3,et1,et2,et3,tol)
    ! Computes the derivative of the logarithmic strain with respect to the   
    ! Left Cauchy-Green deformation tensor (according to CS Jog)
    IMPLICIT NONE

    REAL(iwp), INTENT(IN) :: et1(:), et2(:), et3(:), e1, e2, e3, lne1, lne2,  &
     lne3, tol
    REAL(iwp), INTENT(OUT) :: lnderiv(:,:)
    REAL(iwp), PARAMETER :: one=1.0_iwp

    ! Calculate the derivative of the logarithmic strain with respect to the   
    ! Left Cauchy-Green tensor
    ! If all eigenvalues are different
    IF ((DABS(e1-e2)>tol).and.(DABS(e2-e3)>tol).and.(DABS(e1-e3)>tol)) THEN
      lnderiv=(one/e1)*tensor_p_iklj(et1,et1)+                                &
       (one/e2)*tensor_p_iklj(et2,et2)+                                       &
       (one/e3)*tensor_p_iklj(et3,et3)+                                       &
       ((lne1-lne2)/(e1-e2))*tensor_p_iklj(et1,et2)+                          &
       ((lne1-lne3)/(e1-e3))*tensor_p_iklj(et1,et3)+                          &
       ((lne2-lne1)/(e2-e1))*tensor_p_iklj(et2,et1)+                          &
       ((lne2-lne3)/(e2-e3))*tensor_p_iklj(et2,et3)+                          &
       ((lne3-lne1)/(e3-e1))*tensor_p_iklj(et3,et1)+                          &
       ((lne3-lne2)/(e3-e2))*tensor_p_iklj(et3,et2)

    ! If two eigenvalues are equal    
    ! If e1 and e2 are equal
    ELSEIF ((DABS(e1-e2)<tol).and.(DABS(e2-e3)>tol).and.(DABS(e1-e3)>tol)) THEN
      lnderiv=(one/e1)*tensor_p_iklj(et1,et1)+                                &
       (one/e3)*tensor_p_iklj(et3,et3)+                                       &
       ((lne1-lne3)/(e1-e3))*tensor_p_iklj(et1,et3)+                          &
       ((lne3-lne1)/(e3-e1))*tensor_p_iklj(et3,et1)
    
    ! If e1 and e3 are equal
    ELSEIF ((DABS(e1-e3)<tol).and.(DABS(e1-e2)>tol).and.(DABS(e2-e3)>tol)) THEN
      lnderiv=(one/e1)*tensor_p_iklj(et1,et1)+                                &
       (one/e2)*tensor_p_iklj(et2,et2)+                                       &
       ((lne1-lne2)/(e1-e2))*tensor_p_iklj(et1,et2)+                          &
       ((lne2-lne1)/(e2-e1))*tensor_p_iklj(et2,et1)

    ! If e2 and e3 are equal
    ELSEIF ((DABS(e2-e3)<tol).and.(DABS(e1-e2)>tol).and.(DABS(e1-e3)>tol)) THEN
      lnderiv=(one/e1)*tensor_p_iklj(et1,et1)+                                &
       (one/e2)*tensor_p_iklj(et2,et2)+                                       &
       ((lne1-lne2)/(e1-e2))*tensor_p_iklj(et1,et2)+                          &
       ((lne2-lne1)/(e2-e1))*tensor_p_iklj(et2,et1)

    ! If all eigenvalues are equal
    ELSE
      lnderiv=(one/e1)*tensor_p_iklj(et1,et1)
    END IF

  RETURN
  END SUBROUTINE LN_DERIV

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  FUNCTION TENSOR_P_IKLJ(tensora,tensorb)
    ! Produces the tensorial product of two second order tensors, in order to 
    ! form a fourth order tensor (three dimensions)
    IMPLICIT NONE

    REAL(iwp), INTENT(IN) :: tensora(:), tensorb(:)
    REAL(iwp) :: tensor_p_iklj(9,9)
    
    tensor_p_iklj(1,1)=tensora(1)*tensorb(1) !1111
    tensor_p_iklj(1,2)=tensora(4)*tensorb(1) !1121
    tensor_p_iklj(1,3)=tensora(6)*tensorb(1) !1131
    tensor_p_iklj(1,4)=tensora(1)*tensorb(4) !1112
    tensor_p_iklj(1,5)=tensora(4)*tensorb(4) !1122
    tensor_p_iklj(1,6)=tensora(6)*tensorb(4) !1132
    tensor_p_iklj(1,7)=tensora(1)*tensorb(6) !1113
    tensor_p_iklj(1,8)=tensora(4)*tensorb(6) !1123
    tensor_p_iklj(1,9)=tensora(6)*tensorb(6) !1133
    
    tensor_p_iklj(2,1)=tensor_p_iklj(1,2) !2111
    tensor_p_iklj(2,2)=tensora(2)*tensorb(1) !2121
    tensor_p_iklj(2,3)=tensora(5)*tensorb(1) !2131
    tensor_p_iklj(2,4)=tensora(4)*tensorb(4) !2112
    tensor_p_iklj(2,5)=tensora(2)*tensorb(4) !2122
    tensor_p_iklj(2,6)=tensora(5)*tensorb(4) !2132
    tensor_p_iklj(2,7)=tensora(4)*tensorb(6) !2113
    tensor_p_iklj(2,8)=tensora(2)*tensorb(6) !2123
    tensor_p_iklj(2,9)=tensora(5)*tensorb(6) !2133

    tensor_p_iklj(3,1)=tensor_p_iklj(1,3) !3111
    tensor_p_iklj(3,2)=tensor_p_iklj(2,3) !3121
    tensor_p_iklj(3,3)=tensora(3)*tensorb(1) !3131
    tensor_p_iklj(3,4)=tensora(6)*tensorb(4) !3112
    tensor_p_iklj(3,5)=tensora(5)*tensorb(4) !3122
    tensor_p_iklj(3,6)=tensora(3)*tensorb(4) !3132
    tensor_p_iklj(3,7)=tensora(6)*tensorb(6) !3113
    tensor_p_iklj(3,8)=tensora(5)*tensorb(6) !3123
    tensor_p_iklj(3,9)=tensora(3)*tensorb(6) !3133

    tensor_p_iklj(4,1)=tensor_p_iklj(1,4) !1211
    tensor_p_iklj(4,2)=tensor_p_iklj(2,4) !1221
    tensor_p_iklj(4,3)=tensor_p_iklj(3,4) !1231
    tensor_p_iklj(4,4)=tensora(1)*tensorb(2) !1212
    tensor_p_iklj(4,5)=tensora(4)*tensorb(2) !1222
    tensor_p_iklj(4,6)=tensora(6)*tensorb(2) !1232
    tensor_p_iklj(4,7)=tensora(1)*tensorb(5) !1213
    tensor_p_iklj(4,8)=tensora(4)*tensorb(5) !1223
    tensor_p_iklj(4,9)=tensora(6)*tensorb(5) !1233
    
    tensor_p_iklj(5,1)=tensor_p_iklj(1,5) !2211
    tensor_p_iklj(5,2)=tensor_p_iklj(2,5) !2221
    tensor_p_iklj(5,3)=tensor_p_iklj(3,5) !2231
    tensor_p_iklj(5,4)=tensor_p_iklj(4,5) !2212
    tensor_p_iklj(5,5)=tensora(2)*tensorb(2) !2222
    tensor_p_iklj(5,6)=tensora(5)*tensorb(2) !2232
    tensor_p_iklj(5,7)=tensora(4)*tensorb(5) !2213
    tensor_p_iklj(5,8)=tensora(2)*tensorb(5) !2223
    tensor_p_iklj(5,9)=tensora(5)*tensorb(5) !2233
    
    tensor_p_iklj(6,1)=tensor_p_iklj(1,6) !3211
    tensor_p_iklj(6,2)=tensor_p_iklj(2,6) !3221
    tensor_p_iklj(6,3)=tensor_p_iklj(3,6) !3231
    tensor_p_iklj(6,4)=tensor_p_iklj(4,6) !3212
    tensor_p_iklj(6,5)=tensor_p_iklj(5,6) !3222
    tensor_p_iklj(6,6)=tensora(3)*tensorb(2) !3232
    tensor_p_iklj(6,7)=tensora(6)*tensorb(5) !3213
    tensor_p_iklj(6,8)=tensora(5)*tensorb(5) !3223
    tensor_p_iklj(6,9)=tensora(3)*tensorb(5) !3233
    
    tensor_p_iklj(7,1)=tensor_p_iklj(1,7) !1311
    tensor_p_iklj(7,2)=tensor_p_iklj(2,7) !1321
    tensor_p_iklj(7,3)=tensor_p_iklj(3,7) !1331
    tensor_p_iklj(7,4)=tensor_p_iklj(4,7) !1312
    tensor_p_iklj(7,5)=tensor_p_iklj(5,7) !1322
    tensor_p_iklj(7,6)=tensor_p_iklj(6,7) !1332
    tensor_p_iklj(7,7)=tensora(1)*tensorb(3) !1313
    tensor_p_iklj(7,8)=tensora(4)*tensorb(3) !1323
    tensor_p_iklj(7,9)=tensora(6)*tensorb(3) !1333
    
    tensor_p_iklj(8,1)=tensor_p_iklj(1,8) !2311
    tensor_p_iklj(8,2)=tensor_p_iklj(2,8) !2321
    tensor_p_iklj(8,3)=tensor_p_iklj(3,8) !2331
    tensor_p_iklj(8,4)=tensor_p_iklj(4,8) !2312
    tensor_p_iklj(8,5)=tensor_p_iklj(5,8) !2322
    tensor_p_iklj(8,6)=tensor_p_iklj(6,8) !2332
    tensor_p_iklj(8,7)=tensor_p_iklj(7,8) !2313
    tensor_p_iklj(8,8)=tensora(2)*tensorb(3) !2323
    tensor_p_iklj(8,9)=tensora(5)*tensorb(3) !2333

    tensor_p_iklj(9,1)=tensor_p_iklj(1,9) !3311
    tensor_p_iklj(9,2)=tensor_p_iklj(2,9) !3321
    tensor_p_iklj(9,3)=tensor_p_iklj(3,9) !3331
    tensor_p_iklj(9,4)=tensor_p_iklj(4,9) !3312
    tensor_p_iklj(9,5)=tensor_p_iklj(5,9) !3322
    tensor_p_iklj(9,6)=tensor_p_iklj(6,9) !3332
    tensor_p_iklj(9,7)=tensor_p_iklj(7,9) !3313
    tensor_p_iklj(9,8)=tensor_p_iklj(8,9) !3323
    tensor_p_iklj(9,9)=tensora(3)*tensorb(3) !3333

  END FUNCTION TENSOR_P_IKLJ

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  
  SUBROUTINE BDERIVF(b_tensor,bderiv)
    ! Computes the derivative of the Left Cauchy-Green deformation tensor with
    ! respect to the deformation gradient
    IMPLICIT NONE

    REAL(iwp), INTENT(IN) :: b_tensor(:)
    REAL(iwp), INTENT(OUT) :: bderiv(:,:)
    REAL(iwp), PARAMETER :: zero=0.0_iwp, two=2.0_iwp
    
    bderiv(1,1)=two*b_tensor(1)
    bderiv(1,2)=b_tensor(4)
    bderiv(1,3)=b_tensor(6)
    bderiv(1,4)=b_tensor(4)
    bderiv(1,5)=zero
    bderiv(1,6)=zero
    bderiv(1,7)=b_tensor(6)
    bderiv(1,8)=zero
    bderiv(1,9)=zero
    
    bderiv(2,1)=zero
    bderiv(2,2)=b_tensor(1)
    bderiv(2,3)=zero
    bderiv(2,4)=b_tensor(1)
    bderiv(2,5)=two*b_tensor(4)
    bderiv(2,6)=b_tensor(6)
    bderiv(2,7)=zero
    bderiv(2,8)=b_tensor(6)
    bderiv(2,9)=zero
    
    bderiv(3,1)=zero
    bderiv(3,2)=zero
    bderiv(3,3)=b_tensor(1)
    bderiv(3,4)=zero
    bderiv(3,5)=zero
    bderiv(3,6)=b_tensor(4)
    bderiv(3,7)=b_tensor(1)
    bderiv(3,8)=b_tensor(4)
    bderiv(3,9)=two*b_tensor(6)
    
    bderiv(4,1)=two*b_tensor(4)
    bderiv(4,2)=b_tensor(2)
    bderiv(4,3)=b_tensor(5)
    bderiv(4,4)=b_tensor(2)
    bderiv(4,5)=zero
    bderiv(4,6)=zero
    bderiv(4,7)=b_tensor(5)
    bderiv(4,8)=zero
    bderiv(4,9)=zero
    
    bderiv(5,1)=zero
    bderiv(5,2)=b_tensor(4)
    bderiv(5,3)=zero
    bderiv(5,4)=b_tensor(4)
    bderiv(5,5)=two*b_tensor(2)
    bderiv(5,6)=b_tensor(5)
    bderiv(5,7)=zero
    bderiv(5,8)=b_tensor(5)
    bderiv(5,9)=zero
    
    bderiv(6,1)=zero
    bderiv(6,2)=zero
    bderiv(6,3)=b_tensor(4)
    bderiv(6,4)=zero
    bderiv(6,5)=zero
    bderiv(6,6)=b_tensor(2)
    bderiv(6,7)=b_tensor(4)
    bderiv(6,8)=b_tensor(2)
    bderiv(6,9)=two*b_tensor(5)
    
    bderiv(7,1)=two*b_tensor(6)
    bderiv(7,2)=b_tensor(5)
    bderiv(7,3)=b_tensor(3)
    bderiv(7,4)=b_tensor(5)
    bderiv(7,5)=zero
    bderiv(7,6)=zero
    bderiv(7,7)=b_tensor(3)
    bderiv(7,8)=zero
    bderiv(7,9)=zero
    
    bderiv(8,1)=zero
    bderiv(8,2)=b_tensor(6)
    bderiv(8,3)=zero
    bderiv(8,4)=b_tensor(6)
    bderiv(8,5)=two*b_tensor(5)
    bderiv(8,6)=b_tensor(3)
    bderiv(8,7)=zero
    bderiv(8,8)=b_tensor(3)
    bderiv(8,9)=zero
    
    bderiv(9,1)=zero
    bderiv(9,2)=zero
    bderiv(9,3)=b_tensor(6)
    bderiv(9,4)=zero
    bderiv(9,5)=zero
    bderiv(9,6)=b_tensor(5)
    bderiv(9,7)=b_tensor(6)
    bderiv(9,8)=b_tensor(5)
    bderiv(9,9)=two*b_tensor(3)

  RETURN
  END SUBROUTINE BDERIVF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE GEOM_COMPONENT(sigma1C,geom_comp)
    ! Calculates the geometric component of the corresponding integration point
    IMPLICIT NONE

    REAL(iwp), INTENT(IN) :: sigma1C(:)
    REAL(iwp), INTENT(OUT) :: geom_comp(:,:)
    REAL(iwp), PARAMETER :: zero=0.0_iwp
    
    geom_comp(1,1)=-sigma1C(1)
    geom_comp(1,2)=zero
    geom_comp(1,3)=zero
    geom_comp(1,4)=-sigma1C(4)
    geom_comp(1,5)=zero
    geom_comp(1,6)=zero
    geom_comp(1,7)=-sigma1C(6)
    geom_comp(1,8)=zero
    geom_comp(1,9)=zero
    
    geom_comp(2,1)=zero
    geom_comp(2,2)=-sigma1C(1)
    geom_comp(2,3)=zero
    geom_comp(2,4)=zero
    geom_comp(2,5)=-sigma1C(4)
    geom_comp(2,6)=zero
    geom_comp(2,7)=zero
    geom_comp(2,8)=-sigma1C(6)
    geom_comp(2,9)=zero
    
    geom_comp(3,1)=zero
    geom_comp(3,2)=zero
    geom_comp(3,3)=-sigma1C(1)
    geom_comp(3,4)=zero
    geom_comp(3,5)=zero
    geom_comp(3,6)=-sigma1C(4)
    geom_comp(3,7)=zero
    geom_comp(3,8)=zero
    geom_comp(3,9)=-sigma1C(6)
    
    geom_comp(4,1)=-sigma1C(4)
    geom_comp(4,2)=zero
    geom_comp(4,3)=zero
    geom_comp(4,4)=-sigma1C(2)
    geom_comp(4,5)=zero
    geom_comp(4,6)=zero
    geom_comp(4,7)=-sigma1C(5)
    geom_comp(4,8)=zero
    geom_comp(4,9)=zero
    
    geom_comp(5,1)=zero
    geom_comp(5,2)=-sigma1C(4)
    geom_comp(5,3)=zero
    geom_comp(5,4)=zero
    geom_comp(5,5)=-sigma1C(2)
    geom_comp(5,6)=zero
    geom_comp(5,7)=zero
    geom_comp(5,8)=-sigma1C(5)
    geom_comp(5,9)=zero
    
    geom_comp(6,1)=zero
    geom_comp(6,2)=zero
    geom_comp(6,3)=-sigma1C(4)
    geom_comp(6,4)=zero
    geom_comp(6,5)=zero
    geom_comp(6,6)=-sigma1C(2)
    geom_comp(6,7)=zero
    geom_comp(6,8)=zero
    geom_comp(6,9)=-sigma1C(5)
    
    geom_comp(7,1)=-sigma1C(6)
    geom_comp(7,2)=zero
    geom_comp(7,3)=zero
    geom_comp(7,4)=-sigma1C(5)
    geom_comp(7,5)=zero
    geom_comp(7,6)=zero
    geom_comp(7,7)=-sigma1C(3)
    geom_comp(7,8)=zero
    geom_comp(7,9)=zero
    
    geom_comp(8,1)=zero
    geom_comp(8,2)=-sigma1C(6)
    geom_comp(8,3)=zero
    geom_comp(8,4)=zero
    geom_comp(8,5)=-sigma1C(5)
    geom_comp(8,6)=zero
    geom_comp(8,7)=zero
    geom_comp(8,8)=-sigma1C(3)
    geom_comp(8,9)=zero
    
    geom_comp(9,1)=zero
    geom_comp(9,2)=zero
    geom_comp(9,3)=-sigma1C(6)
    geom_comp(9,4)=zero
    geom_comp(9,5)=zero
    geom_comp(9,6)=-sigma1C(5)
    geom_comp(9,7)=zero
    geom_comp(9,8)=zero
    geom_comp(9,9)=-sigma1C(3)
    
    RETURN
  END SUBROUTINE GEOM_COMPONENT

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE GEEF_CENTROID(auxm,coord,geeF0,jacF0,detF0,ndim,nod)
    ! Evaluates the deformation gradient, its determinant, and the discrete   
    ! spatial gradient operator at the element centroid
    IMPLICIT NONE

    INTEGER,   INTENT(IN)  :: ndim, nod
    REAL(iwp), INTENT(IN)  :: auxm(:,:), coord(:,:)
    REAL(iwp), INTENT(OUT) :: geeF0(:,:), jacF0(:,:), detF0

    REAL(iwp),ALLOCATABLE :: der(:,:), jac(:,:), deriv(:,:), jacinv(:,:),     &
                             jacFinv(:,:), derivF(:,:), derivFtran(:,:)

    ALLOCATE(der(ndim,nod), jac(ndim,ndim), deriv(ndim,nod),                  &
     jacinv(ndim,ndim), jacFinv(ndim,ndim), derivF(ndim,nod),                 &
     derivFtran(nod,ndim))
      
    der(1,1) =  0.125_iwp
    der(2,1) = -0.125_iwp
    der(3,1) = -0.125_iwp
    der(1,2) =  0.125_iwp
    der(2,2) =  0.125_iwp
    der(3,2) = -0.125_iwp
    der(1,3) = -0.125_iwp
    der(2,3) =  0.125_iwp
    der(3,3) = -0.125_iwp
    der(1,4) = -0.125_iwp
    der(2,4) = -0.125_iwp
    der(3,4) = -0.125_iwp
    der(1,5) =  0.125_iwp
    der(2,5) = -0.125_iwp
    der(3,5) =  0.125_iwp
    der(1,6) =  0.125_iwp
    der(2,6) =  0.125_iwp
    der(3,6) =  0.125_iwp
    der(1,7) = -0.125_iwp
    der(2,7) =  0.125_iwp
    der(3,7) =  0.125_iwp
    der(1,8) = -0.125_iwp
    der(2,8) = -0.125_iwp
    der(3,8) =  0.125_iwp

    jac = MATMUL(der,coord)
    CALL INVERT_MATRIX_3x3(jac,jacinv)
    deriv = MATMUL(jacinv,der)

    jacF0 = MATMUL(TRANSPOSE(auxm),TRANSPOSE(deriv))
    jacF0(1,1) = jacF0(1,1) + 1.0_iwp
    jacF0(2,2) = jacF0(2,2) + 1.0_iwp
    jacF0(3,3) = jacF0(3,3) + 1.0_iwp
    CALL DETERMINANT_MATRIX_3x3(jacF0,detF0)
    CALL INVERT_MATRIX_3x3(jacF0,jacFinv)
    
    derivFtran = MATMUL(TRANSPOSE(deriv),jacFinv)
    derivF = TRANSPOSE(derivFtran)
    CALL GEEMAT(derivF,geeF0)

  RETURN

  END SUBROUTINE GEEF_CENTROID

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE DEFGRA(igauss,auxm,coord,points,det,detF,beeF,geeF,jacF,ndim,nod)
    ! Calculates the deformation gradient and provides details about the      
    ! transformation between normal and natural coordinates as well as the 
    ! spatial discrete gradient operator
    IMPLICIT NONE

    INTEGER,   INTENT(IN)  :: igauss, ndim, nod
    REAL(iwp), INTENT(IN)  :: auxm(:,:), coord(:,:), points(:,:)
    REAL(iwp), INTENT(OUT) :: det, detF, beeF(:,:), jacF(:,:), geeF(:,:)

    REAL(iwp),ALLOCATABLE :: der(:,:), jac(:,:), deriv(:,:), jacFinv(:,:),    &
                             derivFtran(:,:), jacinv(:,:), derivF(:,:)

    ALLOCATE(der(ndim,nod), jac(ndim,ndim), deriv(ndim,nod),                  &
     jacinv(ndim,ndim), jacFinv(ndim,ndim), derivFtran(nod,ndim),             &
     derivF(ndim,nod))
      
    CALL SHAPE_DERIVATIVES(igauss,points,der)
    jac = MATMUL(der,coord)
    CALL DETERMINANT_MATRIX_3x3(jac,det)
    CALL INVERT_MATRIX_3x3(jac,jacinv)
    deriv = MATMUL(jacinv,der)

    jacF = MATMUL(TRANSPOSE(auxm),TRANSPOSE(deriv))
    jacF(1,1) = jacF(1,1) + 1.0_iwp
    jacF(2,2) = jacF(2,2) + 1.0_iwp
    jacF(3,3) = jacF(3,3) + 1.0_iwp
    CALL DETERMINANT_MATRIX_3x3(jacF,detF)
    CALL INVERT_MATRIX_3x3(jacF,jacFinv)

    derivFtran = MATMUL(TRANSPOSE(deriv),jacFinv)
    derivF = TRANSPOSE(derivFtran)
!   CALL BEEMAT(derivF,beeF) !- Francisco Calvo format
    CALL BEEMAT(beeF,derivF) !- Book format
    CALL GEEMAT(derivF,geeF)

  RETURN
  END SUBROUTINE DEFGRA

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE DEFGRAINC(igauss,auxm_inc,upd_coord,points,jacFinc,ndim,nod)
    ! Calculates the incremental deformation gradient
    IMPLICIT NONE

    INTEGER,   INTENT(IN)  :: igauss, nod, ndim
    REAL(iwp), INTENT(IN)  :: auxm_inc(:,:), upd_coord(:,:), points(:,:)
    REAL(iwp), INTENT(OUT) :: jacFinc(:,:)

    REAL(iwp), ALLOCATABLE :: der(:,:), jac(:,:), jacinv(:,:), deriv(:,:)     

    ALLOCATE(der(ndim,nod), jac(ndim,ndim), jacinv(ndim,ndim), deriv(ndim,nod))
      
    CALL SHAPE_DERIVATIVES(igauss,points,der)
    jac = MATMUL(der,upd_coord)
    CALL INVERT_MATRIX_3x3(jac,jacinv)
    deriv = MATMUL(jacinv,der)

    jacFinc = MATMUL(TRANSPOSE(auxm_inc),TRANSPOSE(deriv))
    jacFinc(1,1) = jacFinc(1,1) + 1.0_iwp
    jacFinc(2,2) = jacFinc(2,2) + 1.0_iwp
    jacFinc(3,3) = jacFinc(3,3) + 1.0_iwp

  RETURN
  END SUBROUTINE DEFGRAINC

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE GEEMAT(derivF,geeF)
    ! Calculates the discrete spatial gradient operator
    IMPLICIT NONE

    REAL(iwp),   INTENT(IN)  :: derivF(:,:)
    REAL(iwp), INTENT(OUT) :: geeF(:,:)
    INTEGER :: i, j, k, nod

    geeF=0.0_iwp
    nod=UBOUND(derivF,2)

    DO i=1,UBOUND(geeF,1)
      IF (i<4) THEN
        k=i
        DO j=1,nod
          geeF(i,(k+(3*(j-1))))=derivF(1,j)
        END DO
      ELSEIF ((i<7) .and. (i>3)) THEN
        k=i-3
        DO j=1,nod
          geeF(i,(k+(3*(j-1))))=derivF(2,j)
        END DO
      ELSE
        k=i-6
        DO j=1,nod
          geeF(i,(k+(3*(j-1))))=derivF(3,j)
        END DO
      END IF
    END DO

  RETURN
  END SUBROUTINE GEEMAT

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE DOUBLE_CONTRACTION_4(tensor_4_out,tensor_4_1,tensor_4_2)
    ! Performs the double contractions of two fourth order tensors (3x3x3x3)
    IMPLICIT NONE

    real(iwp), INTENT(IN) :: tensor_4_1(:,:,:,:), tensor_4_2(:,:,:,:)
    real(iwp), INTENT(OUT) :: tensor_4_out(:,:,:,:)
    integer :: i, j, k, l, m, n

    tensor_4_out=0.0_iwp

    DO i=1,3
      DO j=1,3
        DO k=1,3
          DO l=1,3
            DO m=1,3
              DO n=1,3
                tensor_4_out(i,j,k,l)=tensor_4_out(i,j,k,l)+                  &
                 tensor_4_1(i,j,m,n)*tensor_4_2(m,n,k,l)
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  RETURN
  END SUBROUTINE DOUBLE_CONTRACTION_4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE PARAFEM_TO_LARGE_STRAIN_VOIGT(large_strain_voigt,tensor)
    ! Transforms a fourth order tensor into a large strain Voigt notation (so 
    ! that the tensor product of a fourth order tensor with two non-symmetric 
    ! tensors give the same value as a traditional matrix multiplication)
    IMPLICIT NONE

    real(iwp), INTENT(OUT) :: large_strain_voigt(:,:)
    real(iwp), INTENT(IN) :: tensor(:,:,:,:)

    large_strain_voigt(1,1)=tensor(1,1,1,1)
    large_strain_voigt(1,2)=tensor(1,1,2,1)
    large_strain_voigt(1,3)=tensor(1,1,3,1)
    large_strain_voigt(1,4)=tensor(1,1,1,2)
    large_strain_voigt(1,5)=tensor(1,1,2,2)
    large_strain_voigt(1,6)=tensor(1,1,3,2)
    large_strain_voigt(1,7)=tensor(1,1,1,3)
    large_strain_voigt(1,8)=tensor(1,1,2,3)
    large_strain_voigt(1,9)=tensor(1,1,3,3)

    large_strain_voigt(2,1)=tensor(2,1,1,1)
    large_strain_voigt(2,2)=tensor(2,1,2,1)
    large_strain_voigt(2,3)=tensor(2,1,3,1)
    large_strain_voigt(2,4)=tensor(2,1,1,2)
    large_strain_voigt(2,5)=tensor(2,1,2,2)
    large_strain_voigt(2,6)=tensor(2,1,3,2)
    large_strain_voigt(2,7)=tensor(2,1,1,3)
    large_strain_voigt(2,8)=tensor(2,1,2,3)
    large_strain_voigt(2,9)=tensor(2,1,3,3)

    large_strain_voigt(3,1)=tensor(3,1,1,1)
    large_strain_voigt(3,2)=tensor(3,1,2,1)
    large_strain_voigt(3,3)=tensor(3,1,3,1)
    large_strain_voigt(3,4)=tensor(3,1,1,2)
    large_strain_voigt(3,5)=tensor(3,1,2,2)
    large_strain_voigt(3,6)=tensor(3,1,3,2)
    large_strain_voigt(3,7)=tensor(3,1,1,3)
    large_strain_voigt(3,8)=tensor(3,1,2,3)
    large_strain_voigt(3,9)=tensor(3,1,3,3)

    large_strain_voigt(4,1)=tensor(1,2,1,1)
    large_strain_voigt(4,2)=tensor(1,2,2,1)
    large_strain_voigt(4,3)=tensor(1,2,3,1)
    large_strain_voigt(4,4)=tensor(1,2,1,2)
    large_strain_voigt(4,5)=tensor(1,2,2,2)
    large_strain_voigt(4,6)=tensor(1,2,3,2)
    large_strain_voigt(4,7)=tensor(1,2,1,3)
    large_strain_voigt(4,8)=tensor(1,2,2,3)
    large_strain_voigt(4,9)=tensor(1,2,3,3)

    large_strain_voigt(5,1)=tensor(2,2,1,1)
    large_strain_voigt(5,2)=tensor(2,2,2,1)
    large_strain_voigt(5,3)=tensor(2,2,3,1)
    large_strain_voigt(5,4)=tensor(2,2,1,2)
    large_strain_voigt(5,5)=tensor(2,2,2,2)
    large_strain_voigt(5,6)=tensor(2,2,3,2)
    large_strain_voigt(5,7)=tensor(2,2,1,3)
    large_strain_voigt(5,8)=tensor(2,2,2,3)
    large_strain_voigt(5,9)=tensor(2,2,3,3)

    large_strain_voigt(6,1)=tensor(3,2,1,1)
    large_strain_voigt(6,2)=tensor(3,2,2,1)
    large_strain_voigt(6,3)=tensor(3,2,3,1)
    large_strain_voigt(6,4)=tensor(3,2,1,2)
    large_strain_voigt(6,5)=tensor(3,2,2,2)
    large_strain_voigt(6,6)=tensor(3,2,3,2)
    large_strain_voigt(6,7)=tensor(3,2,1,3)
    large_strain_voigt(6,8)=tensor(3,2,2,3)
    large_strain_voigt(6,9)=tensor(3,2,3,3)

    large_strain_voigt(7,1)=tensor(1,3,1,1)
    large_strain_voigt(7,2)=tensor(1,3,2,1)
    large_strain_voigt(7,3)=tensor(1,3,3,1)
    large_strain_voigt(7,4)=tensor(1,3,1,2)
    large_strain_voigt(7,5)=tensor(1,3,2,2)
    large_strain_voigt(7,6)=tensor(1,3,3,2)
    large_strain_voigt(7,7)=tensor(1,3,1,3)
    large_strain_voigt(7,8)=tensor(1,3,2,3)
    large_strain_voigt(7,9)=tensor(1,3,3,3)

    large_strain_voigt(8,1)=tensor(2,3,1,1)
    large_strain_voigt(8,2)=tensor(2,3,2,1)
    large_strain_voigt(8,3)=tensor(2,3,3,1)
    large_strain_voigt(8,4)=tensor(2,3,1,2)
    large_strain_voigt(8,5)=tensor(2,3,2,2)
    large_strain_voigt(8,6)=tensor(2,3,3,2)
    large_strain_voigt(8,7)=tensor(2,3,1,3)
    large_strain_voigt(8,8)=tensor(2,3,2,3)
    large_strain_voigt(8,9)=tensor(2,3,3,3)

    large_strain_voigt(9,1)=tensor(3,3,1,1)
    large_strain_voigt(9,2)=tensor(3,3,2,1)
    large_strain_voigt(9,3)=tensor(3,3,3,1)
    large_strain_voigt(9,4)=tensor(3,3,1,2)
    large_strain_voigt(9,5)=tensor(3,3,2,2)
    large_strain_voigt(9,6)=tensor(3,3,3,2)
    large_strain_voigt(9,7)=tensor(3,3,1,3)
    large_strain_voigt(9,8)=tensor(3,3,2,3)
    large_strain_voigt(9,9)=tensor(3,3,3,3)

  RETURN
  END SUBROUTINE PARAFEM_TO_LARGE_STRAIN_VOIGT

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE FBAR_COMPONENT(deeF,sigma,fbar_comp)
    ! Calculates the addition to the tangent operator of the F bar scheme
    IMPLICIT NONE

    real(iwp), INTENT(IN) :: deeF(:,:), sigma(:,:)
    real(iwp), INTENT(OUT) :: fbar_comp(:,:)
    real(iwp) ::  tensor(3,3,3,3), q_tensor(3,3,3,3), unit_tensor(3,3),       &
     a_ixi(3,3,3,3), ixi(3,3,3,3)
    real(iwp), PARAMETER :: zero=0._iwp, one=1._iwp, two=2._iwp, three=3._iwp
    integer :: i, j, k, l

    unit_tensor=zero
    DO i=1,3
      unit_tensor(i,i)=one
    END DO

    CALL large_strain_voigt_to_parafem(deeF,tensor)

    ixi=tensor_p_ijkl(unit_tensor,unit_tensor)

    CALL double_contraction_4(a_ixi,tensor,ixi)

    q_tensor=(one/three)*a_ixi-(two/three)*tensor_p_ijkl(sigma,unit_tensor)

    CALL parafem_to_large_strain_voigt(fbar_comp,q_tensor)

  RETURN
  END SUBROUTINE FBAR_COMPONENT

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE LARGE_STRAIN_VOIGT_TO_PARAFEM(large_strain_voigt,tensor)
    ! Transforms a large strain Voigt notation tensor to a full fourth order 
    ! tensor
    IMPLICIT NONE

    real(iwp), INTENT(IN) :: large_strain_voigt(:,:)
    real(iwp), INTENT(OUT) :: tensor(:,:,:,:)

    tensor(1,1,1,1)=large_strain_voigt(1,1)
    tensor(1,1,2,1)=large_strain_voigt(1,2)
    tensor(1,1,3,1)=large_strain_voigt(1,3)
    tensor(1,1,1,2)=large_strain_voigt(1,4)
    tensor(1,1,2,2)=large_strain_voigt(1,5)
    tensor(1,1,3,2)=large_strain_voigt(1,6)
    tensor(1,1,1,3)=large_strain_voigt(1,7)
    tensor(1,1,2,3)=large_strain_voigt(1,8)
    tensor(1,1,3,3)=large_strain_voigt(1,9)

    tensor(2,1,1,1)=large_strain_voigt(2,1)
    tensor(2,1,2,1)=large_strain_voigt(2,2)
    tensor(2,1,3,1)=large_strain_voigt(2,3)
    tensor(2,1,1,2)=large_strain_voigt(2,4)
    tensor(2,1,2,2)=large_strain_voigt(2,5)
    tensor(2,1,3,2)=large_strain_voigt(2,6)
    tensor(2,1,1,3)=large_strain_voigt(2,7)
    tensor(2,1,2,3)=large_strain_voigt(2,8)
    tensor(2,1,3,3)=large_strain_voigt(2,9)

    tensor(3,1,1,1)=large_strain_voigt(3,1)
    tensor(3,1,2,1)=large_strain_voigt(3,2)
    tensor(3,1,3,1)=large_strain_voigt(3,3)
    tensor(3,1,1,2)=large_strain_voigt(3,4)
    tensor(3,1,2,2)=large_strain_voigt(3,5)
    tensor(3,1,3,2)=large_strain_voigt(3,6)
    tensor(3,1,1,3)=large_strain_voigt(3,7)
    tensor(3,1,2,3)=large_strain_voigt(3,8)
    tensor(3,1,3,3)=large_strain_voigt(3,9)

    tensor(1,2,1,1)=large_strain_voigt(4,1)
    tensor(1,2,2,1)=large_strain_voigt(4,2)
    tensor(1,2,3,1)=large_strain_voigt(4,3)
    tensor(1,2,1,2)=large_strain_voigt(4,4)
    tensor(1,2,2,2)=large_strain_voigt(4,5)
    tensor(1,2,3,2)=large_strain_voigt(4,6)
    tensor(1,2,1,3)=large_strain_voigt(4,7)
    tensor(1,2,2,3)=large_strain_voigt(4,8)
    tensor(1,2,3,3)=large_strain_voigt(4,9)

    tensor(2,2,1,1)=large_strain_voigt(5,1)
    tensor(2,2,2,1)=large_strain_voigt(5,2)
    tensor(2,2,3,1)=large_strain_voigt(5,3)
    tensor(2,2,1,2)=large_strain_voigt(5,4)
    tensor(2,2,2,2)=large_strain_voigt(5,5)
    tensor(2,2,3,2)=large_strain_voigt(5,6)
    tensor(2,2,1,3)=large_strain_voigt(5,7)
    tensor(2,2,2,3)=large_strain_voigt(5,8)
    tensor(2,2,3,3)=large_strain_voigt(5,9)

    tensor(3,2,1,1)=large_strain_voigt(6,1)
    tensor(3,2,2,1)=large_strain_voigt(6,2)
    tensor(3,2,3,1)=large_strain_voigt(6,3)
    tensor(3,2,1,2)=large_strain_voigt(6,4)
    tensor(3,2,2,2)=large_strain_voigt(6,5)
    tensor(3,2,3,2)=large_strain_voigt(6,6)
    tensor(3,2,1,3)=large_strain_voigt(6,7)
    tensor(3,2,2,3)=large_strain_voigt(6,8)
    tensor(3,2,3,3)=large_strain_voigt(6,9)

    tensor(1,3,1,1)=large_strain_voigt(7,1)
    tensor(1,3,2,1)=large_strain_voigt(7,2)
    tensor(1,3,3,1)=large_strain_voigt(7,3)
    tensor(1,3,1,2)=large_strain_voigt(7,4)
    tensor(1,3,2,2)=large_strain_voigt(7,5)
    tensor(1,3,3,2)=large_strain_voigt(7,6)
    tensor(1,3,1,3)=large_strain_voigt(7,7)
    tensor(1,3,2,3)=large_strain_voigt(7,8)
    tensor(1,3,3,3)=large_strain_voigt(7,9)

    tensor(2,3,1,1)=large_strain_voigt(8,1)
    tensor(2,3,2,1)=large_strain_voigt(8,2)
    tensor(2,3,3,1)=large_strain_voigt(8,3)
    tensor(2,3,1,2)=large_strain_voigt(8,4)
    tensor(2,3,2,2)=large_strain_voigt(8,5)
    tensor(2,3,3,2)=large_strain_voigt(8,6)
    tensor(2,3,1,3)=large_strain_voigt(8,7)
    tensor(2,3,2,3)=large_strain_voigt(8,8)
    tensor(2,3,3,3)=large_strain_voigt(8,9)

    tensor(3,3,1,1)=large_strain_voigt(9,1)
    tensor(3,3,2,1)=large_strain_voigt(9,2)
    tensor(3,3,3,1)=large_strain_voigt(9,3)
    tensor(3,3,1,2)=large_strain_voigt(9,4)
    tensor(3,3,2,2)=large_strain_voigt(9,5)
    tensor(3,3,3,2)=large_strain_voigt(9,6)
    tensor(3,3,1,3)=large_strain_voigt(9,7)
    tensor(3,3,2,3)=large_strain_voigt(9,8)
    tensor(3,3,3,3)=large_strain_voigt(9,9)

  RETURN
  END SUBROUTINE LARGE_STRAIN_VOIGT_TO_PARAFEM

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  FUNCTION TENSOR_P_IJKL(tensora,tensorb)
    ! Produces the tensorial product of two second order tensors, in order to 
    ! form a fourth order tensor (three dimensions)
    IMPLICIT NONE

    REAL(iwp), INTENT(IN) :: tensora(:,:), tensorb(:,:)
    REAL(iwp) :: tensor_p_ijkl(3,3,3,3)
    integer :: i, j, k, l

    DO i=1,3
      DO j=1,3
        DO k=1,3
          DO l=1,3
            tensor_p_ijkl(i,j,k,l)=tensora(i,j)*tensorb(k,l)
          END DO
        END DO
      END DO
    END DO   

  END FUNCTION TENSOR_P_IJKL
  
  SUBROUTINE UMAT(stress,statev,ddsdde,stran,dstran,ntens,statevar_num,iel,igauss)

    ! This subroutine returns the updated stress, strain and the tangent 
    ! operator for the Eccentric-Ellipsoid model with (linear) isotropic 
    ! hardening and associative plastic flow rule
    IMPLICIT NONE
    
    REAL(iwp), INTENT(OUT) :: stress(:), ddsdde(:,:)
	REAL(iwp), INTENT(INOUT) :: stran(:), statev(:), dstran(:)
    REAL(iwp) :: scalar_term, eqplas, e, nu, yield_t, yield_c, f_lower_0,     &
     f_upper_0, syield, hard, sfs, zeta, stress_eq, plastic_mul, norm_flow,   &
     yield, ran_scalar, norm_solution, nen, eqplas_trial, res_piv, alpha,     &
     mprod, mprev, mder, alpha1, alpha2
    REAL(iwp) :: unit_tensor(6), ixi(6,6), ixi_sym(6,6), f_4(6,6), f_2(6),    &
     flow_dir(7), fs(6), res_strain(6), dnds(6,6), unit_7(7,7),               &
     strain_trial(6), residuals(8), results(8), jacobian(8,8), fd(7),         &
     inv_jacobian(8,8), e_comp(6,6), nxn(6,6), en(6), g_mat(7,7),             &
     flow_dir_stress(6), grad_flow(7,7), dir(8), results0(8)
    REAL(iwp), PARAMETER :: zero=0._iwp, one=1._iwp, two=2._iwp,              &
	 tol=0.000001_iwp, half=0.5_iwp, tol_nr=0.00000001_iwp, beta=0.0001_iwp,  &
     ls_const=0.1_iwp, four=4._iwp
    INTEGER, INTENT(IN) :: ntens, statevar_num
    INTEGER :: i, j, iter, iel, igauss, ls_iter
    INTEGER, PARAMETER :: max_nr_iter=50, max_ls_iter=25
     
	! Assign material properties (user defined)
    ! Yield properties obtained from Wolfram et al. 2012
    e=12700._iwp
    nu=0.3_iwp
    yield_t=52._iwp !(0.41%)
    yield_c=105._iwp !(0.83%)
    !zeta=0.2_iwp
    zeta=0.49_iwp
    !hard=0.001 ! Corresponds to 0% of the elastic slope
    hard=0.296_iwp ! Corresponds to 5% of the elastic slope in tension
    hard=0.038_iwp ! Corresponds to 5% of the elastic slope in compression
     
    ! Assign derived material properties
    f_upper_0=(yield_t+yield_c)/(two*yield_t*yield_c)
    f_lower_0=(one/two)*((one/yield_t)-(one/yield_c))
       
    ! Recover the previous equivalent plastic strain
    eqplas_trial=statev(1)

    ! Initializing variables
	ddsdde=zero
    stress=zero

    ! Compute the linear elastic isotropic stiffness matrix
    scalar_term=e/((one+nu)*(one-two*nu))

    ddsdde(1,1)=one-nu
    ddsdde(2,2)=one-nu
    ddsdde(3,3)=one-nu
    ddsdde(4,4)=(one-two*nu)/two
    ddsdde(5,5)=(one-two*nu)/two
    ddsdde(6,6)=(one-two*nu)/two
    ddsdde(1,2)=nu
    ddsdde(1,3)=nu
    ddsdde(2,1)=nu
    ddsdde(2,3)=nu
    ddsdde(3,1)=nu
    ddsdde(3,2)=nu

    ddsdde=scalar_term*ddsdde
    
    ! Calculate the predictor strain and stress
    DO i=1,ntens
      DO j=1,ntens
        stress(i)=stress(i)+ddsdde(i,j)*stran(j)
      END DO
    END DO
    
    ! Define the unit tensor, IxI and IxI_sym
    ixi=zero
    ixi_sym=zero
    
    DO i=1,3
      unit_tensor(i)=one
      unit_tensor(i+3)=zero
      DO j=1,3
        ixi(i,j)=one
      END DO
      
      ixi_sym(i,i)=one
      ixi_sym(i+3,i+3)=half
    END DO

    ! Define the fourth order tensor F_4 and the second order tensor F_2
    f_4=-zeta*(f_upper_0**2)*ixi+(zeta+one)*(f_upper_0**2)*ixi_sym
    f_2=f_lower_0*unit_tensor
        
    ! Calculate the equivalent stress
    fs=MATMUL(f_4,stress)
    fs(4:6)=four*fs(4:6)
    sfs=DOT_PROD(stress,fs,6)
    sfs=SQRT(sfs)
    stress_eq=sfs+DOT_PROD(f_2,stress,6)
        
    ! Calculate the equivalent yield stress
    syield=one+hard*eqplas_trial
    
    ! Determine if there is yielding
    ! =========================================================================
    IF ((stress_eq-syield)>=tol) THEN

      ! This material point is yielding, proceed with the return-mapping
      ! =======================================================================

      ! Initialise some matrices
      g_mat=zero
      g_mat(1:6,1:6)=ddsdde
      g_mat(7,7)=hard
      
      unit_7=zero
      DO i=1,7
        unit_7(i,i)=one
      END DO
      
      ! Initialize variables before the local Newton-Raphson loop
      plastic_mul=zero
      strain_trial=stran
      
      DO iter=1,max_nr_iter
      
        ! Warn if the maximum number of iterations has been reached
        IF (iter==max_nr_iter) THEN
          WRITE(*,*) 'Maximum local Newton-Raphson iterations have been reached'
        END IF
        
        ! Calculate the flow direction
        IF (iter.EQ.1) THEN
          flow_dir_stress=(fs/sfs)+f_2
      
          ! Compute the residuals 
          yield=sfs+DOT_PROD(f_2,stress,6)-(one+hard*(plastic_mul+eqplas_trial))
         
          ! Assemble the residual and the results vector
          residuals=zero
          results(1:6)=stran(1:6)  
          
          residuals(7)=zero
          results(7)=-eqplas_trial
          
          residuals(8)=yield
          results(8)=zero
        END IF
        
        ! Assemble the flow direction
        flow_dir(1:6)=flow_dir_stress(1:6)
        flow_dir(7)=one
         
        ! Calculate the Jacobian
        ! =====================================================================
         
        ! Calculate the derivative of the flow vector with respect to stress
        fs=MATMUL(f_4,stress)
        fs(4:6)=four*fs(4:6)
        
        dnds=(f_4/sfs)
        dnds(4:6,4:6)=four*dnds(4:6,4:6)
        
        dnds=dnds-TENSOR_PRODUCT_22(fs,(fs/(sfs**3)))
         
        grad_flow=zero
        grad_flow(1:6,1:6)=dnds
                 
        ! Assemble the jacobian
        jacobian(1:7,1:7)=unit_7+plastic_mul*MATMUL(grad_flow,g_mat)
        fd=MATMUL(flow_dir,g_mat)
        
        jacobian(8,1:7)=fd(1:7)
        jacobian(1:7,8)=flow_dir(1:7)        
        jacobian(8,8)=zero
        
        ! Invert the Jacobian
        CALL INVERSE(jacobian,inv_jacobian,8)
        ! =====================================================================
         
        ! Compute direction of advance
        dir=-MATMUL(inv_jacobian,residuals)
         
        ! Line search scheme
        ! =====================================================================
        
        ! Set up some initial results
        alpha=one
        mprod=half*DOT_PROD(residuals,residuals,8)
        mprev=mprod
        mder=-two*mprod
        results0=results
        
        DO ls_iter=1,max_ls_iter
        
          ! Update new results
          results=results0+alpha*dir
          plastic_mul=results(8)
          eqplas=-results(7)
          stran(1:6)=results(1:6)
          stress=MATMUL(ddsdde,stran)
          
          fs=MATMUL(f_4,stress)
          fs(4:6)=four*fs(4:6)
          sfs=DOT_PROD(stress,fs,6)
          sfs=SQRT(sfs)
          flow_dir_stress=(fs/sfs)+f_2
          
          res_strain=stran-strain_trial+plastic_mul*flow_dir_stress
          res_piv=-eqplas+eqplas_trial+plastic_mul
          yield=sfs+DOT_PROD(f_2,stress,6)-(one+hard*(plastic_mul+eqplas_trial))
         
          residuals(1:6)=res_strain(1:6)
          residuals(7)=res_piv
          residuals(8)=yield
          
          mprod=half*DOT_PROD(residuals,residuals,8)
          
          IF (mprod<=((one-two*ls_const*beta)*mprev)) THEN
            EXIT
          ELSE
            alpha1=ls_const*alpha
            alpha2=(-(alpha**2)*mder)/(two*(mprod-mprev-alpha*mder))
            
            IF (alpha1>=alpha2) THEN
              alpha=alpha1
            ELSE
              alpha=alpha2
            END IF
          END IF
        END DO
        ! =====================================================================
         
        ! Exit if convergence
        ! =====================================================================
        norm_solution=zero
        DO i=1,8
          norm_solution=norm_solution+(residuals(i)**2)
        END DO
        norm_solution=SQRT(norm_solution)
         
        IF (norm_solution<=tol_nr) THEN
          EXIT
        END IF
        ! =====================================================================
        
      END DO
      
      ! Update equivalent plastic strain     
      eqplas=eqplas_trial+plastic_mul

      ! Assemble tangent operator
      CALL INVERSE(ddsdde,e_comp,6)
      dnds=e_comp+plastic_mul*dnds
      CALL INVERSE(dnds,ddsdde,6)
      en=MATMUL(ddsdde,flow_dir_stress)
      nxn=tensor_product_22(en,en)
      
      DO i=4,6
        flow_dir(i)=half*flow_dir(i)
      END DO
      
      nen=double_contraction_22(flow_dir,en)+hard
      ddsdde=ddsdde-(one/nen)*nxn
      
      ! Update state variables
      statev(1)=eqplas

    END IF
    
    ! End of yielding
    ! =========================================================================
    
  RETURN
  END SUBROUTINE UMAT
  
  FUNCTION DOUBLE_CONTRACTION_22(vector_input_1,vector_input_2)
    ! This function performs a double contraction of two second order tensors, 
    ! in vector notation, as following, alpha=AijCij
    IMPLICIT NONE
    
	real(iwp), INTENT(IN) :: vector_input_1(:), vector_input_2(:)
    real(iwp) :: double_contraction_22
    integer :: i
    
    double_contraction_22=0._iwp
    
    DO i=1,3
      double_contraction_22=double_contraction_22+vector_input_1(i)*          &
       vector_input_2(i)
    END DO
    
    DO i=4,6
      double_contraction_22=double_contraction_22+2._iwp*vector_input_1(i)*   &
       vector_input_2(i)
    END DO
    
  RETURN
  END FUNCTION DOUBLE_CONTRACTION_22
  
  FUNCTION TENSOR_PRODUCT_22(vector_input_1,vector_input_2)
    ! This function calculates a fourth order tensor through the tensorial 
    ! product of two second order tensor, in matrix notation, as following
    ! Aijkl=BijCkl
    IMPLICIT NONE
    
	real(iwp), INTENT(IN) :: vector_input_1(:), vector_input_2(:)
    real(iwp) :: tensor_product_22(6,6)
    integer :: i, j
    
    DO i=1,6
      DO j=1,6
        tensor_product_22(i,j)=vector_input_1(i)*vector_input_2(j)
      END DO
    END DO
    
  RETURN
  END FUNCTION TENSOR_PRODUCT_22
  
  FUNCTION DOT_PROD(vector_input_1,vector_input_2,dimen)
    ! This function performs a double contraction of two second order tensors, 
    ! in vector notation, as following, dot_prod=uivi, with dimen being the 
    ! dimension of both vectors
    IMPLICIT NONE
    
	real(iwp), INTENT(IN) :: vector_input_1(:), vector_input_2(:)
    integer, INTENT(IN) :: dimen
    real(iwp) :: dot_prod
    integer :: i
    
    dot_prod=0._iwp
    
    DO i=1,dimen
      dot_prod=dot_prod+vector_input_1(i)*vector_input_2(i)
    END DO
       
  RETURN
  END FUNCTION DOT_PROD
  
  SUBROUTINE INVERSE(a,c,n)
    !This subroutine calculates the inverse of a nxn matrix
    ! Input is a(n,n)
    ! n is the dimension
    ! c is the inverse of a    
    IMPLICIT NONE

    INTEGER :: n, i, j, k  
    REAL(iwp) :: a(n,n), c(n,n), l(n,n), u(n,n), b(n), d(n), x(n)
    REAL(iwp) :: coeff
    
    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    l=0._iwp
    u=0._iwp
    b=0._iwp

    ! step 1: forward elimination
    DO k=1,(n-1)
      DO i=(k+1),n
        coeff=a(i,k)/a(k,k)
        l(i,k)=coeff
        DO j=(k+1),n
          a(i,j)=a(i,j)-coeff*a(k,j)
        END DO
      END DO
    END DO

    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    DO i=1,n
      l(i,i)=1._iwp
    END DO
    
    ! U matrix is the upper triangular part of A
    DO j=1,n
      DO i=1,j
        u(i,j)=a(i,j)
      END DO
    END DO

    ! Step 3: compute columns of the inverse matrix C
    DO k=1,n
      b(k)=1._iwp
      d(1)=b(1)
      
      ! Step 3a: Solve Ld=b using the forward substitution
      DO i=2,n
        d(i)=b(i)
        DO j=1,(i-1)
          d(i)=d(i)-l(i,j)*d(j)
        END DO
      END DO
      
      ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/u(n,n)
      DO i=(n-1),1,-1
        x(i)=d(i)
        DO j=n,(i+1),-1
        x(i)=x(i)-u(i,j)*x(j)
        END DO
        x(i)=x(i)/u(i,i)
      END DO
      
      ! Step 3c: fill the solutions x(n) into column k of C
      DO i=1,n
        c(i,k)=x(i)
      END DO
      b(k)=0._iwp
    END DO
    
  RETURN
  END SUBROUTINE INVERSE  

 SUBROUTINE UMAT1(stress,statev,ddsdde,stran,dstran,ntens,statevar_num,iel,igauss)

    REAL(iwp), INTENT(OUT) :: stress(:), ddsdde(:,:)
	REAL(iwp), INTENT(INOUT) :: stran(:), statev(:), dstran(:)
    REAL(iwp) :: scalar_term, e, nu
    REAL(iwp), PARAMETER :: zero=0._iwp, one=1._iwp, two=2._iwp
    INTEGER, INTENT(IN) :: ntens, statevar_num, iel, igauss
    INTEGER :: i, j
	
    ! Assign material properties
    e=12700._iwp
    nu=0.3_iwp
    
    ! Setting the fourth and second-order tensors to zero
	ddsdde=zero
    stress=zero

    ! Compute the linear elastic isotropic stiffness matrix
    scalar_term=e/((one+nu)*(one-two*nu))

    ddsdde(1,1)=one-nu
    ddsdde(2,2)=one-nu
    ddsdde(3,3)=one-nu
    ddsdde(4,4)=(one-two*nu)/two
    ddsdde(5,5)=(one-two*nu)/two
    ddsdde(6,6)=(one-two*nu)/two
    ddsdde(1,2)=nu
    ddsdde(1,3)=nu
    ddsdde(2,1)=nu
    ddsdde(2,3)=nu
    ddsdde(3,1)=nu
    ddsdde(3,2)=nu

    ddsdde=scalar_term*ddsdde

    ! Calculate the predictor strain and stress
    DO i=1,ntens
      DO j=1,ntens
        stress(i)=stress(i)+ddsdde(i,j)*stran(j)
      END DO
    END DO
    
  RETURN
  END SUBROUTINE UMAT1

END MODULE LARGE_STRAIN  
