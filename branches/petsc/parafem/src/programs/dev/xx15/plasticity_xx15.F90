MODULE plasticity_xx15
  
  USE PRECISION
  USE MATHS
  USE NEW_LIBRARY
  
CONTAINS
  
  SUBROUTINE PLASTICITY(deeF,jacF,jacFinc,sigma1C,statev,lnstrainelas,         &
                        sigma,detF,statevar_num,iel,igauss,noncon_flag,umat)
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
    INTEGER, INTENT(IN) :: statevar_num, iel, igauss
    REAL(iwp), INTENT(OUT) :: deeF(:,:), sigma(:,:), sigma1C(:)
    REAL(iwp), INTENT(INOUT) :: lnstrainelas(:), statev(:)
    LOGICAL, INTENT(INOUT)   :: noncon_flag
    INTERFACE
      SUBROUTINE umat(stress,statev,ddsdde,stran,dstran,ntens,                 &
                      statevar_num,iel,igauss,noncon_flag)
        USE PRECISION
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ntens
        INTEGER, INTENT(IN) :: statevar_num, iel, igauss
        REAL(iwp), INTENT(OUT) :: stress(:), ddsdde(:,:)
        REAL(iwp), INTENT(INOUT) :: stran(:), statev(:), dstran(:)
        LOGICAL, INTENT(INOUT)   :: noncon_flag
      END SUBROUTINE UMAT
    END INTERFACE
    
    REAL(iwp) :: d_tensor(9,9), lnderiv(9,9), strderb(9,9), bderiv(9,9),      &
     geom_comp(9,9), ddsdde(6,6), b_tensor_3x3(3,3), et1(6), et2(6), et3(6),  &
     b_tensor(6), dstran(6), btensor(6), lne1, lne2, lne3, e1, e2, e3
    REAL(iwp), PARAMETER :: half=0.5_iwp, one=1._iwp, two=2._iwp,             &
     tol=0.000001_iwp
    INTEGER :: ntens
       
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

    WRITE(*,'(A,3f10.6)') "K", lne1, lne2, lne3
    WRITE(*,'(A,6f10.6)') "K", lnstrainelas
    
    ! Calculate the Cauchy stress tensor and the consistent tangent operator 
    ! (material contribution)
    ! Calling UMAT, which consists in the same format as the one used in
    ! ABAQUS. The only difference is that the material properties should be
    ! defined inside of the subroutine
    sigma1C=0._iwp
    CALL umat(sigma1C,statev,ddsdde,lnstrainelas,dstran,ntens,statevar_num,    &
     iel,igauss,noncon_flag)

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
    
    WRITE(*,'(A,81f10.2)') "K", lnderiv

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

    REAL(iwp) :: d1, d2, d3
    REAL(iwp) :: et1_kopp(6), et2_kopp(6), et3_kopp(6), e1_kopp, e2_kopp, e3_kopp
    REAL(iwp) :: et1_lapack(6), et2_lapack(6), et3_lapack(6), e1_lapack, e2_lapack, e3_lapack

    CALL EIGEN_KOPP(e1_kopp,e2_kopp,e3_kopp,tensor,et1_kopp,et2_kopp,et3_kopp,tol)
    CALL EIGEN_LAPACK(e1_lapack,e2_lapack,e3_lapack,tensor,et1_lapack,et2_lapack,et3_lapack,tol)

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
    ELSEIF (ISNAN(alpha)) THEN
      alpha=one
    END IF
    
    alpha=DACOS(alpha)
    pi=4._iwp*DATAN(one)

    d1=-two*DSQRT(q)*DCOS(alpha/three)+sinv1/three
    d2=-two*DSQRT(q)*DCOS((alpha+two*pi)/three)+sinv1/three
    d3=-two*DSQRT(q)*DCOS((alpha-two*pi)/three)+sinv1/three

!!$    ! Unsorted
!!$    e1 = d1
!!$    e2 = d2
!!$    e3 = d3
    ! Sort the eigenvalues to match dsyev
    IF (d1 < d2) THEN
      IF (d1 < d3) THEN
        e1 = d1
        IF (d2 < d3) THEN
          e2 = d2; e3 = d3
        ELSE ! d3 <= d2
          e2 = d3; e3 = d2
        END IF
      ELSE ! d3 <= d1
        e1 = d3; e2 = d1; e3 = d2
      END IF
    ELSE ! d2 <= d1
      IF (d2 < d3) THEN
        e1 = d2
        IF (d1 < d3) THEN
          e2 = d1; e3 = d3
        ELSE ! d3 <= d1
          e2 = d3; e3 = d1
        END IF
      ELSE ! d3 <= d2
        e1 = d3; e2 = d2; e3 = d1
      END IF
    END IF
    
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

    WRITE(*,'(A,3f6.2)') "P", e1, e2, e3
    WRITE(*,'(A,6f6.2)') "P", et1
    WRITE(*,'(A,6f6.2)') "P", et2
    WRITE(*,'(A,6f6.2)') "P", et3
    WRITE(*,'(A,3f6.2)') "L", e1_lapack, e2_lapack, e3_lapack
    WRITE(*,'(A,6f6.2)') "L", et1_lapack
    WRITE(*,'(A,6f6.2)') "L", et2_lapack
    WRITE(*,'(A,6f6.2)') "L", et3_lapack
    WRITE(*,'(A,3f6.2)') "K", e1_kopp, e2_kopp, e3_kopp
    WRITE(*,'(A,6f6.2)') "K", et1_kopp
    WRITE(*,'(A,6f6.2)') "K", et2_kopp
    WRITE(*,'(A,6f6.2)') "K", et3_kopp

    CALL EIGEN_KOPP(e1,e2,e3,tensor,et1,et2,et3,tol)

  RETURN
  END SUBROUTINE EIGEN

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE EIGEN_LAPACK(e1,e2,e3,tensor,et1,et2,et3,tol)
    ! This subroutine calculates the eigenvalues and eigenprojections of a
    ! second order tensor using LAPACK.
    IMPLICIT NONE

    EXTERNAL dsyev

    REAL(iwp), INTENT(IN) :: tensor(:), tol
    REAL(iwp), INTENT(OUT) :: et1(:), et2(:), et3(:), e1, e2, e3

    INTEGER, PARAMETER :: lwork = 3*3-1
    DOUBLE PRECISION   :: a(3,3), w(3), work(lwork)
    INTEGER            :: info
    DOUBLE PRECISION   :: p1(6), p2(6), p3(6), t(3,3), v(3,1)

    a(1,1) = tensor(1)
    a(2,2) = tensor(2)
    a(3,3) = tensor(3)
    a(1,2) = tensor(4)
    a(2,3) = tensor(5)
    a(1,3) = tensor(6)

    CALL dsyev('V','U',3,A,3,w,work,lwork,info)

    e1 = w(1)
    e2 = w(2)
    e3 = w(3)

!!$    v(:,1) = a(:,1)
!!$    t = MATMUL(v,transpose(v))
!!$    et1 = (/t(1,1),t(2,2),t(3,3),t(1,2),t(2,3),t(1,3)/)
!!$    v(:,1) = a(:,2)
!!$    t = MATMUL(v,transpose(v))
!!$    et2 = (/t(1,1),t(2,2),t(3,3),t(1,2),t(2,3),t(1,3)/)
!!$    v(:,1) = a(:,3)
!!$    t = MATMUL(v,transpose(v))
!!$    et3 = (/t(1,1),t(2,2),t(3,3),t(1,2),t(2,3),t(1,3)/)
    v(:,1) = a(:,1)
    t = MATMUL(v,transpose(v))
    p1 = (/t(1,1),t(2,2),t(3,3),t(1,2),t(2,3),t(1,3)/)
    v(:,1) = a(:,2)
    t = MATMUL(v,transpose(v))
    p2 = (/t(1,1),t(2,2),t(3,3),t(1,2),t(2,3),t(1,3)/)
    v(:,1) = a(:,3)
    t = MATMUL(v,transpose(v))
    p3 = (/t(1,1),t(2,2),t(3,3),t(1,2),t(2,3),t(1,3)/)

    ! If all eigenvalues are different
    IF ((DABS(e1-e2)>tol).AND.(DABS(e2-e3)>tol).AND.(DABS(e1-e3)>tol)) THEN
      et1 = p1
      et2 = p2
      et3 = p3
    ! If two eigenvalues are equal
    ! If e1 and e2 are equal
    ELSEIF ((DABS(e1-e2)<tol).AND.(DABS(e2-e3)>tol).AND.(DABS(e1-e3)>tol)) THEN
      et1 = p1 + p2
      et2 = 0
      et3 = p3
    ! If e1 and e3 are equal
    ELSEIF ((DABS(e1-e3)<tol).AND.(DABS(e1-e2)>tol).AND.(DABS(e2-e3)>tol)) THEN
      et1 = p1 + p3
      et2 = p2
      et3 = 0
    ! If e2 and e3 are equal
    ELSEIF ((DABS(e2-e3)<tol).AND.(DABS(e1-e2)>tol).AND.(DABS(e1-e3)>tol)) THEN
      et1 = p1
      et2 = p2 + p3
      et3 = 0
    ! If all eigenvalues are equal
    ELSE
      et1 = p1 + p2 + p3
      et2 = 0
      et3 = 0
    END IF
  END SUBROUTINE EIGEN_LAPACK

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE EIGEN_KOPP(e1,e2,e3,tensor,et1,et2,et3,tol)
    ! This subroutine calculates the eigenvalues and eigenprojections of a
    ! second order tensor using Kopp's hybrid algorithm for real symmetric
    ! matrices.
    IMPLICIT NONE

    INTERFACE
      SUBROUTINE DSYEVH3(A, Q, W)
        DOUBLE PRECISION A(3,3)
        DOUBLE PRECISION Q(3,3)
        DOUBLE PRECISION W(3)
      END SUBROUTINE DSYEVH3
    END INTERFACE

    REAL(iwp), INTENT(IN) :: tensor(:), tol
    REAL(iwp), INTENT(OUT) :: et1(:), et2(:), et3(:), e1, e2, e3

    DOUBLE PRECISION :: a(3,3), q(3,3), w(3)
    DOUBLE PRECISION :: p1(6), p2(6), p3(6), r1(6), r2(6), r3(6), t(3,3), v(3,1)

    a(1,1) = tensor(1)
    a(2,2) = tensor(2)
    a(3,3) = tensor(3)
    a(1,2) = tensor(4)
    a(2,3) = tensor(5)
    a(1,3) = tensor(6)

    CALL dsyevh3(a,q,w)

    v(:,1) = q(:,1)
    t = MATMUL(v,transpose(v))
    r1 = (/t(1,1),t(2,2),t(3,3),t(1,2),t(2,3),t(1,3)/)
    v(:,1) = q(:,2)
    t = MATMUL(v,transpose(v))
    r2 = (/t(1,1),t(2,2),t(3,3),t(1,2),t(2,3),t(1,3)/)
    v(:,1) = q(:,3)
    t = MATMUL(v,transpose(v))
    r3 = (/t(1,1),t(2,2),t(3,3),t(1,2),t(2,3),t(1,3)/)

!!$!!$    ! Unsorted
!!$!!$    e1 = w(1)
!!$!!$    e2 = w(2)
!!$!!$    e3 = w(3)
!!$!!$    et1 = p1
!!$!!$    et2 = p2
!!$!!$    et3 = p3
!!$    ! Sort the eigenvalues to match dsyev
!!$    IF (w(1) < w(2)) THEN
!!$      IF (w(1) < w(3)) THEN
!!$        e1 = w(1); et1 = p1
!!$        IF (w(2) < w(3)) THEN
!!$          e2 = w(2); et2 = p2; e3 = w(3); et3 = p3
!!$        ELSE ! w(3) <= w(2)
!!$          e2 = w(3); et2 = p3; e3 = w(2); et3 = p2
!!$        END IF
!!$      ELSE ! w(3) <= w(1)
!!$        e1 = w(3); et1 = p3; e2 = w(1); et2 = p1; e3 = w(2); et3 = p2
!!$      END IF
!!$    ELSE ! w(2) <= w(1)
!!$      IF (w(2) < w(3)) THEN
!!$        e1 = w(2); et1 = p2
!!$        IF (w(1) < w(3)) THEN
!!$          e2 = w(1); et2 = p1; e3 = w(3); et3 = p3
!!$        ELSE ! w(3) <= w(1)
!!$          e2 = w(3); et2 = p3; e3 = w(1); et3 = p1
!!$        END IF
!!$      ELSE ! w(3) <= w(2)
!!$        e1 = w(3); et1 = p3; e2 = w(2); et2 = p2; e3 = w(1); et3 = p1
!!$      END IF
!!$    END IF

!!$    ! Unsorted
!!$    e1 = w(1)
!!$    e2 = w(2)
!!$    e3 = w(3)
!!$    et1 = p1
!!$    et2 = p2
!!$    et3 = p3
    ! Sort the eigenvalues to match dsyev
    IF (w(1) < w(2)) THEN
      IF (w(1) < w(3)) THEN
        e1 = w(1); p1 = r1
        IF (w(2) < w(3)) THEN
          e2 = w(2); p2 = r2; e3 = w(3); p3 = r3
        ELSE ! w(3) <= w(2)
          e2 = w(3); p2 = r3; e3 = w(2); p3 = r2
        END IF
      ELSE ! w(3) <= w(1)
        e1 = w(3); p1 = r3; e2 = w(1); p2 = r1; e3 = w(2); p3 = r2
      END IF
    ELSE ! w(2) <= w(1)
      IF (w(2) < w(3)) THEN
        e1 = w(2); p1 = r2
        IF (w(1) < w(3)) THEN
          e2 = w(1); p2 = r1; e3 = w(3); p3 = r3
        ELSE ! w(3) <= w(1)
          e2 = w(3); p2 = r3; e3 = w(1); p3 = r1
        END IF
      ELSE ! w(3) <= w(2)
        e1 = w(3); p1 = r3; e2 = w(2); p2 = r2; e3 = w(1); p3 = r1
      END IF
    END IF
    ! If all eigenvalues are different
    IF ((DABS(e1-e2)>tol).AND.(DABS(e2-e3)>tol).AND.(DABS(e1-e3)>tol)) THEN
      et1 = p1
      et2 = p2
      et3 = p3
    ! If two eigenvalues are equal
    ! If e1 and e2 are equal
    ELSEIF ((DABS(e1-e2)<tol).AND.(DABS(e2-e3)>tol).AND.(DABS(e1-e3)>tol)) THEN
      et1 = p1 + p2
      et2 = 0
      et3 = p3
    ! If e1 and e3 are equal
    ELSEIF ((DABS(e1-e3)<tol).AND.(DABS(e1-e2)>tol).AND.(DABS(e2-e3)>tol)) THEN
      et1 = p1 + p3
      et2 = p2
      et3 = 0
    ! If e2 and e3 are equal
    ELSEIF ((DABS(e2-e3)<tol).AND.(DABS(e1-e2)>tol).AND.(DABS(e1-e3)>tol)) THEN
      et1 = p1
      et2 = p2 + p3
      et3 = 0
    ! If all eigenvalues are equal
    ELSE
      et1 = p1 + p2 + p3
      et2 = 0
      et3 = 0
    END IF
  END SUBROUTINE EIGEN_KOPP

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

END MODULE plasticity_xx15
