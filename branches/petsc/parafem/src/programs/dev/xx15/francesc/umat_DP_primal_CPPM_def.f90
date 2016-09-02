MODULE MATERIAL

USE precision

CONTAINS

  SUBROUTINE UMAT(stress,statev,ddsdde,stran,dstran,ntens,statevar_num,iel,   &
   igauss,noncon_flag)

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
    LOGICAL, INTENT(INOUT) :: noncon_flag
     
	! Assign material properties (user defined)
    ! Yield properties obtained from Wolfram et al. 2012
    e=12700._iwp
    nu=0.3_iwp
    yield_t=52._iwp !(0.41%)
    yield_c=105._iwp !(0.83%)
    !zeta=0.2_iwp
    zeta=0.49_iwp
    !hard=0.001 ! Corresponds to 0% of the elastic slope
    !hard=0.296_iwp ! Corresponds to 5% of the elastic slope in tension
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
          noncon_flag=.TRUE.
          EXIT
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
  
END MODULE MATERIAL