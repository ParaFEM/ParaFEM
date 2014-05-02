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
  
END MODULE LARGE_STRAIN  
