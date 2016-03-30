MODULE MATERIAL

USE precision

CONTAINS

  SUBROUTINE UMAT(stress,statev,ddsdde,stran,dstran,ntens,nstatv,iel,igauss,noncon_flag)

    ! This subroutine returns the updated stress, strain and the tangent 
    ! operator for the integrated constitutive law
    real(iwp), INTENT(IN) :: dstran(:)
    real(iwp), INTENT(INOUT) :: stran(:), statev(:)
    real(iwp), INTENT(INOUT) :: stress(:), ddsdde(:,:)
    
    integer, INTENT(IN) :: ntens, nstatv, iel, igauss
    integer :: i, j, iter

    real(iwp) :: bulk_mod, shear_mod, scalar_term, smises, syield, eqplas, &
     tr, plastic_mul, e, nu, h, syield0, tol, sdev_norm, smises_trial, delta, tangent, &
     res_der, yield_fun, stress_inf
    real(iwp) :: stress_dev(6), proj_tensor(6,6), flow_flow(6,6), unit_tensor(6,6), &
     ixi_tensor(6,6), flow(6), term(6,6), unit_sym(6,6)
     
    real(iwp), PARAMETER :: zero=0._iwp, half=0.5_iwp, one=1._iwp, two=2._iwp, three=3._iwp, six=6._iwp
    integer, PARAMETER :: max_iter=25
    
    logical, INTENT(INOUT) :: noncon_flag
   
    ! Assign the tolerance to check converge of the newton-raphson
    tol=0.000001_iwp

    ! Assign material properties
    e=206900._iwp
    nu=0.29_iwp
    h=129.24_iwp
    syield0=450._iwp
    stress_inf=715._iwp
    delta=16.93_iwp

    ! Recover the previous equivalent plastic strain
    eqplas=statev(1)

    ! Setting the fourth and second-order tensors to zero
    proj_tensor=zero
    flow_flow=zero
    unit_tensor=zero
    ixi_tensor=zero
    ddsdde=zero
    stress=zero

    ! Assigning values to the shear and bulk modulus
    bulk_mod=e/(three*(one-two*nu))
    shear_mod=e/(two*(one+nu))

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

    stress_dev=stress
    tr=stress(1)+stress(2)+stress(3)
    stress_dev(1)=stress_dev(1)-(one/three)*tr
    stress_dev(2)=stress_dev(2)-(one/three)*tr
    stress_dev(3)=stress_dev(3)-(one/three)*tr

    ! Calculate the equivalent Von Mises stress
    smises_trial=((stress(1)-stress(2))**2)+((stress(2)-stress(3))**2)+ &
     ((stress(3)-stress(1))**2)
    DO i=4,6
      smises_trial=smises_trial+six*(stress(i)**2)
    END DO
    smises_trial=SQRT(smises_trial/two)

    ! Calculate the yield stress
    syield=(stress_inf-syield0)*(one-EXP(-delta*eqplas))+h*eqplas+syield0   

    ! Determine if there is yielding
    IF (smises_trial>(syield+tol)) THEN
 
      ! Calculate the plastic multiplicator through a newton-raphson scheme
      ! Initialize the plastic multiplicator and the residual function
      plastic_mul=zero
      yield_fun=smises_trial-syield

      DO iter=1,max_iter
        ! If convergence is not achieved
        IF (i==max_iter) THEN
          WRITE(*,*) 'Maximum local Newton-Raphson iterations have been reached'
          noncon_flag=.TRUE.
          EXIT
        END IF
        
        ! Calculate the hardening slope
        tangent=(stress_inf-syield0)*delta*EXP(-delta*(eqplas+plastic_mul))+h
      
        ! Calculate the residual derivative
        res_der=-three*shear_mod-tangent

        ! New guess for the plastic multiplicator
        plastic_mul=plastic_mul-(yield_fun/res_der)

        ! Update stresses
        DO i=1,3
          stress(i)=stress_dev(i)*(one-((plastic_mul*three*shear_mod)/smises_trial))+ &
           (one/three)*tr
        END DO

        DO i=4,6
          stress(i)=stress_dev(i)*(one-((plastic_mul*three*shear_mod)/smises_trial))
        END DO
        
        ! Calculate the equivalent Von Mises stress
        smises=((stress(1)-stress(2))**2)+((stress(2)-stress(3))**2)+ &
         ((stress(3)-stress(1))**2)
        DO i=4,6
          smises=smises+six*(stress(i)**2)
        END DO
        smises=SQRT(smises/two)

        ! Check for convergence
        syield=(stress_inf-syield0)*(one-EXP(-delta*(eqplas+plastic_mul)))+ &
         h*(eqplas+plastic_mul)+syield0
        yield_fun=smises-syield

        IF (yield_fun<tol) THEN
          EXIT
        END IF
      END DO

      ! Create the unit flow vector
      sdev_norm=((stress_dev(1)**2)+(stress_dev(2)**2)+(stress_dev(3)**2)+ &
       (two*stress_dev(4)**2)+(two*stress_dev(5)**2)+(two*stress_dev(6)**2))
      sdev_norm=SQRT(sdev_norm)
      flow=stress_dev/sdev_norm

      ! Update the equivalent plastic strain
      eqplas=eqplas+plastic_mul

      ! Update the elastic strains
      DO i=1,3
        stran(i)=stran(i)-plastic_mul*(SQRT(three/two))*flow(i)
      END DO

      DO i=4,6
        stran(i)=stran(i)-two*plastic_mul*(SQRT(three/two))*flow(i)
      END DO

      ! Assemble the tangent operator
      ! Create the IxI tensor
      DO i=1,3
        DO j=1,3
          ixi_tensor(i,j)=one
        END DO
      END DO

      ! Create the unit tensor
      DO i=1,3
        unit_tensor(i,i)=one
      END DO
      unit_tensor(4,4)=half
      unit_tensor(5,5)=half
      unit_tensor(6,6)=half

      DO i=1,6
        DO j=1,6
          ! Create the projection tensor
          proj_tensor(i,j)=unit_tensor(i,j)-(one/three)*ixi_tensor(i,j)

          ! Create the tensor product of the two flow vectors
          flow_flow(i,j)=flow(i)*flow(j)
          
          
          
          ! MAYBE FLOW_FLOW NEEDS A MODIFICATION HERE
          
          

          ! Assemble all the parts into the final tangent operator
          ddsdde(i,j)=ddsdde(i,j)-((plastic_mul*six*(shear_mod**2))/smises_trial)*proj_tensor(i,j)+    &
          six*(shear_mod**2)*((plastic_mul/smises_trial)-(one/(three*shear_mod+tangent)))*flow_flow(i,j)
        END DO
      END DO
    END IF

    ! Update the state variables
    statev(1)=eqplas

  RETURN
  END SUBROUTINE UMAT
  
END MODULE MATERIAL
