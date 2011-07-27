MODULE bem_lib
!    Library for EBE elasticity and potential problems
!    Ndof replaced by N_dof to avoid any conflict with global variables
USE precision;   IMPLICIT NONE
REAL(iwp):: Pi= 3.14159265359_iwp
CONTAINS

  SUBROUTINE form_s(GG,ell,kappa,omega,Gamma,s)
  ! forms the s vector in bicgstab(l)
  real(iwp),intent(in)::GG(:,:),kappa; integer,intent(in)::ell
  real(iwp),intent(out)::omega,Gamma(:),s(:)
  real(iwp)::HH(ell-1,ell-1),Gamma0(ell+1),p(ell-1),q(ell-1),  &
        Gamma1(ell+1),NGamma0,NGamma1,cosine
            HH= -GG(2:ell,2:ell); HH = inverse(HH)    
            p = MATMUL(HH,GG(2:ell,1))  ;   q = MATMUL(HH,GG(2:ell,ell+1))
            Gamma0(1) = 1.0_iwp; Gamma0(ell+1) = .0_iwp; Gamma0(2:ell) = p
            Gamma1(1) = .0_iwp; Gamma1(ell+1) = 1.0_iwp; Gamma1(2:ell) = q
            NGamma0   = DOT_PRODUCT(Gamma0,matmul(GG,Gamma0))
            NGamma1   = DOT_PRODUCT(Gamma1,matmul(GG,Gamma1))
            omega     = DOT_PRODUCT(Gamma0,matmul(GG,Gamma1))
            cosine = ABS(omega)/SQRT(ABS(NGamma0*NGamma1)); omega=omega/NGamma1
            if(cosine<kappa) omega = (kappa/cosine) * omega
            Gamma = Gamma0 - omega * Gamma1
            s(1:ell) = Gamma(2:ell+1)       ;  s(ell+1) = .0_iwp
  return
  END SUBROUTINE form_s

 
REAL(iwp) FUNCTION norm(x)
! L2 norm of vector x
IMPLICIT NONE
REAL(iwp),INTENT(IN)::x(:)
norm = SQRT(SUM(x**2))
END FUNCTION norm

FUNCTION inverse(matrix)
! Gauss Jordan for small matrices
IMPLICIT NONE
 REAL(iwp)::matrix(:,:)
 REAL(iwp)::inverse(UBOUND(matrix,1),UBOUND(matrix,2))
 INTEGER::i,k,n; REAL(iwp)::con; n=UBOUND(matrix,1)
 DO k=1,n
   con = matrix(k,k); matrix(k,k)=1._iwp
   matrix(k,:)=matrix(k,:)/con
    DO i=1,n
       IF(i/=k) THEN
          con = matrix(i,k);matrix(i,k)=.0_iwp
          matrix(i,:) = matrix(i,:)-matrix(k,:)*con
       END IF
    END DO
 END DO
 inverse = matrix
END FUNCTION inverse


FUNCTION UK(dxr,r,E,ny,Cdim)
!--------------------------------------------
!   FUNDAMENTAL SOLUTION FOR DISPLACEMENTS
!   isotropic material (Kelvin solution)
!--------------------------------------------
IMPLICIT NONE
REAL(iwp),INTENT(IN)                 :: dxr(:)        !   rx/r etc.
REAL(iwp),INTENT(IN)                 :: r             !   r
REAL(iwp),INTENT(IN)                 :: E             !   Young's modulus
REAL(iwp),INTENT(IN)                 :: ny            !   Poisson's ratio
INTEGER,INTENT(IN)	:: Cdim          !   Cartesian dimension
REAL(iwp)    :: UK(Cdim,Cdim)  !   Function returns array of same dim as dxr
REAL(iwp)    :: G,c,c1,onr,clog,conr  !   Temps
G= E/(2.0*(1+ny))
c1= 3.0 - 4.0*ny
SELECT CASE (Cdim)
        CASE (2)                                !     Two-dimensional solution
		c= 1.0/(8.0*Pi*G*(1.0 - ny))
		clog= -c1*LOG(r)
		UK(1,1)= c*(clog + dxr(1)*dxr(1))
		UK(1,2)= c*dxr(1)*dxr(2)
		UK(2,2)= c*(clog + dxr(2)*dxr(2))
		UK(2,1)= UK(1,2)
        CASE(3)                                 !      Three-dimensional solution
		c= 1.0/(16.0*Pi*G*(1.0 - ny))
		conr=c/r
		UK(1,1)= conr*(c1 + dxr(1)*dxr(1))
		UK(1,2)= conr*dxr(1)*dxr(2)
		UK(1,3)= conr*dxr(1)*dxr(3)
		UK(2,1)= UK(1,2)
		UK(2,2)= conr*(c1 + dxr(2)*dxr(2))
		UK(2,3)= conr*dxr(2)*dxr(3)
		UK(3,1)= UK(1,3)
		UK(3,2)= UK(2,3)
		UK(3,3)= conr*(c1 + dxr(3)*dxr(3))
	CASE DEFAULT
END SELECT
RETURN
END FUNCTION UK

FUNCTION TK(dxr,r,Vnor,ny,Cdim)
!--------------------------------------------
!   FUNDAMENTAL SOLUTION FOR TRACTIONS
!   isotropic material (Kelvin solution)
!--------------------------------------------
IMPLICIT NONE
REAL(iwp),INTENT(IN)                 :: dxr(:)         !   rx/r etc.
REAL(iwp),INTENT(IN)                 :: r              !   r
REAL(iwp),INTENT(IN)                 :: Vnor(:)        !   normal vector
REAL(iwp),INTENT(IN)                 :: ny             !   Poisson's ratio
INTEGER,INTENT(IN)      :: Cdim                   !   Cartesian dimension
REAL(iwp) :: TK(Cdim,Cdim)      !   Function returns array of same dim as dxr
REAL(iwp) :: c2,c3,costh,Conr   !   Temps
c3= 1.0 - 2.0*ny
Costh= DOT_PRODUCT (Vnor,dxr)
SELECT CASE (Cdim)
	CASE (2)        
		c2= 1.0/(4.0*Pi*(1.0 - ny))
		Conr= c2/r
                TK(1,1)= -(Conr*(c3 + 2.0*dxr(1)*dxr(1))*Costh)
                TK(1,2)= -(Conr*(2.0*dxr(1)*dxr(2)*Costh +   &
                c3*(Vnor(1)*dxr(2) - Vnor(2)*dxr(1))))
		TK(2,2)= -(Conr*(c3 + 2.0*dxr(2)*dxr(2))*Costh)
                TK(2,1)= -(Conr*(2.0*dxr(1)*dxr(2)*Costh +   &
                c3*(Vnor(2)*dxr(1) - Vnor(1)*dxr(2))))
	CASE(3)           !    Three-dimensional
		c2= 1.0/(8.0*Pi*(1.0 - ny))
		Conr= c2/r**2
TK(1,1)= -Conr*(c3 + 3.0*dxr(1)*dxr(1))*Costh
TK(1,2)= -Conr*(3.0*dxr(1)*dxr(2)*Costh - c3*(Vnor(2)*dxr(1) - Vnor(1)*dxr(2))) 
TK(1,3)= -Conr*(3.0*dxr(1)*dxr(3)*Costh - c3*(Vnor(3)*dxr(1) - Vnor(1)*dxr(3)))
TK(2,1)= -Conr*(3.0*dxr(1)*dxr(2)*Costh - c3*(Vnor(1)*dxr(2) - Vnor(2)*dxr(1)))
TK(2,2)= -Conr*(c3 + 3.0*dxr(2)*dxr(2))*Costh
TK(2,3)= -Conr*(3.0*dxr(2)*dxr(3)*Costh - c3*(Vnor(3)*dxr(2) - Vnor(2)*dxr(3)))
TK(3,1)= -Conr*(3.0*dxr(1)*dxr(3)*Costh - c3*(Vnor(1)*dxr(3) - Vnor(3)*dxr(1)))
TK(3,2)= -Conr*(3.0*dxr(2)*dxr(3)*Costh - c3*(Vnor(2)*dxr(3) - Vnor(3)*dxr(2)))
TK(3,3)= -Conr*(c3 + 3.0*dxr(3)*dxr(3))*Costh
	CASE DEFAULT
END SELECT
END FUNCTION TK

SUBROUTINE SK(TS,DXR,R,C2,C3)
!------------------------------------------------------------
!   KELVIN SOLUTION FOR STRESS  
!    TO BE MULTIPLIED WITH T
!------------------------------------------------------------
REAL(iwp), INTENT(OUT) :: TS(:,:)   ! Fundamental solution
REAL(iwp), INTENT(IN)  :: DXR(:)    ! rx , ry, rz
REAL(iwp), INTENT(IN)  :: R         ! r
REAL(iwp), INTENT(IN)  :: C2,C3     ! Elastic constants
INTEGER ::  Cdim      !   Cartesian dimension
INTEGER :: NSTRES  !  No. of stress components
INTEGER :: II(6), JJ(6) !  sequence of stresses in pseudo-vector
REAL(iwp)    :: A
INTEGER :: I,N,J,K
Cdim= UBOUND(DXR,1)
IF(CDIM == 2) THEN
   NSTRES= 3
   II(1:3)= (/1,2,1/)
   JJ(1:3)= (/1,2,2/)   
ELSE
   NSTRES= 6
   II= (/1,2,3,1,2,3/)
   JJ= (/1,2,3,2,3,1/)
END IF
Coor_directions:&
DO K=1,Cdim
	Stress_components:&
	DO N=1,NSTRES
				I= II(N)
				J= JJ(N)
				A= 0.
				IF(K .EQ. J) A= A + DXR(I)
				IF(I .EQ. J) A= A - DXR(K)
				IF(K .EQ. I) A= A + DXR(J)
				A= A*C3
				TS(N,K)= C2/R*(A + Cdim*DXR(I)*DXR(J)*DXR(K))
				IF(Cdim .EQ. 3) TS(N,K)= TS(N,K)/R
	END DO &
	Stress_components
END DO &
Coor_directions
RETURN
END SUBROUTINE SK

SUBROUTINE RK(US,DXR,R,VNORM,C3,C5,C6,C7,ny)
!------------------------------------------------------------
!    KELVIN SOLUTION FOR STRESS COMPUTATION
!    TO BE MULTIPLIED WITH U
!------------------------------------------------------------
REAL(iwp), INTENT(OUT) :: US(:,:)        ! Fundamental solution
REAL(iwp), INTENT(IN)  :: DXR(:)         ! rx , ry, rz
REAL(iwp), INTENT(IN)  :: R              ! r
REAL(iwp), INTENT(IN)  :: VNORM(:)       ! nx , ny , nz
REAL(iwp), INTENT(IN)  :: C3,C5,C6,C7,ny ! Elastic constants
INTEGER ::  Cdim   !   Cartesian dimension
INTEGER :: NSTRES  !   No. of stress components
INTEGER :: II(6), JJ(6) !  sequence of stresses in pseudo-vector
REAL(iwp)    :: costh, B,C,cny
INTEGER :: I,N,J,K
Cdim= UBOUND(DXR,1) 
IF(CDIM == 2) THEN
   NSTRES= 3
   II(1:3)= (/1,2,1/)
   JJ(1:3)= (/1,2,2/)   
ELSE
	NSTRES= 6                   
	II= (/1,2,3,1,2,3/)
	JJ= (/1,2,3,2,3,1/)
END IF  
COSTH= DOT_Product(dxr,vnorm) 
Cny= Cdim*ny
Coor_directions:&
DO K=1,Cdim
	Stress_components:&
	DO N=1,NSTRES 
				I= II(N)
				J= JJ(N)
				B= 0.
				C= 0.
				IF(I .EQ. J) B= Cdim*C3*DXR(K)
				IF(I .EQ. K) B= B + Cny*DXR(J)
				IF(J .EQ. K) B= B + Cny*DXR(I)
				B= COSTH *(B - C6*DXR(I)*DXR(J)*DXR(K))
				C= DXR(J)*DXR(K)*Cny
				IF(J .EQ.K) C= C + C3
				C= C*VNORM(I)
				B= B+C
				C= DXR(I)*DXR(K)*Cny
				IF(I .EQ. K) C=C + C3
				C= C*VNORM(J)
				B= B+C
				C= DXR(I)*DXR(J)*Cdim*C3
				IF(I .EQ. J) C= C - C7
				C= C*VNORM(K)
				US(N,K)= (B + C)*C5/R/R
				IF(Cdim .EQ. 3) US(N,K)= US(N,K)/R
	END DO &
	Stress_components
END DO &
Coor_directions 
RETURN
END SUBROUTINE RK

SUBROUTINE Trans_mat(v1,v2,v3, T)
!-----------------------------------------------
!  Computes Stress Transformation Matrix in 3-D
!----------------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(IN)     :: v1(3),v2(3),v3(3)       !  unit vectors in orthogonal directions
REAL(iwp), INTENT(OUT)    :: T(6,6)                  !  transformation matrix
REAL(iwp)                 :: v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z   ! temps
v1x= v1(1) ; v1y= v1(2) ; v1z= v1(3)
v2x= v2(1) ; v2y= v2(2) ; v2z= v2(3)
v3x= v3(1) ; v3y= v3(2) ; v3z= V3(3)
!   T?11
T(1,1)= v1x**2  ; T(1,2)= v2x**2 ; T(1,3)= v3x**2
T(2,1)= v1y**2 ;  T(2,2)= v2y**2 ; T(2,3)= v3y**2
T(3,1)= v1z**2 ;  T(3,2)= v2z**2 ; T(3,3)= v3z**2
!   T?21
T(1,4)= 2.0*v1y*v1x ; T(1,5)= 2.0*v2y*v2x ; T(1,6)= 2.0*v3y*v3x
T(2,4)= 2.0*v1y*v1z ; T(2,5)= 2.0*v2y*v2z ;  T(2,6)= 2.0*v3y*v3z
T(3,4)= 2.0*v1x*v1z ; T(3,5)= 2.0*v2x*v2z ;  T(3,6)= 2.0*v3x*v3z
!   T?12
T(4,1)= v1x*v2x ; T(4,2)= v2x*v3x ; T(4,3)= v1x*v3x
T(5,1)= v1y*v2y ; T(5,2)= v2y*v3y ; T(5,3)= v1y*v3y
T(6,1)= v1z*v2z ; T(6,2)= v2z*v3z ; T(6,3)= v1z*v3z
!   T?22
T(4,4)= v1x*v2y+v1y*v2x ; T(4,5)= v2x*v3y+v2y*v3x ; T(4,6)= v1x*v3y+v1y*v3x
T(5,4)= v1y*v2z+v1z*v2y ; T(5,5)= v2y*v3z+v2z*v3y ; T(5,6)= v1y*v3z+v1z*v3y
T(6,4)= v1x*v2z+v1z*v2x ; T(6,5)= v2x*v3z+v2z*v3x ; T(6,6)= v1x*v3z+v1z*v3x
RETURN
END  SUBROUTINE Trans_mat
 
SUBROUTINE D_mat(E,ny,D,Cdim)
!-----------------------------------
!   Computes isotropic D-matrix
!   Plane-strain (Cdim= 2)
!   or 3-D       (Cdim= 3)
!-----------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(IN)   :: E      !  Young's modulus
REAL(iwp), INTENT(IN)   :: ny     !  Poisson's ratio
INTEGER,INTENT(IN) :: Cdim   !  Cartesian Dimension
REAL(iwp), INTENT(OUT)  :: D(:,:) !  D-matrix
REAL(iwp)               :: c1,c2,G
c1= E*(1.0-ny)/( (1.0+ny)*(1.0-2.0*ny) )
c2= ny/(1.0-ny)
G = E/(2.0*(1.0+ny))
D = 0.0
SELECT CASE (Cdim)
	CASE (2)
		D(1,1)= 1.0  ; D(2,2)= 1.0
		D(2,1)= c2   ; D(1,2)= c2
		D(3,3)= G/c1
	CASE (3)         !    3-D
		D(1,1)= 1.0  ;  D(2,2)= 1.0   ;  D(3,3)= 1.0
		D(2,1)= c2   ;  D(1,3)= c2   ;  D(2,3)= c2
		D(1,2)= c2   ;  D(3,1)= c2   ;  D(3,2)= c2
		D(4,4)= G/c1 ;  D(5,5)= G/c1 ;  D(6,6)= G/c1
	CASE DEFAULT
END SELECT
D= c1*D
RETURN
END SUBROUTINE D_mat

SUBROUTINE D_mat_anis(D,E1,G1,E2,G2,ny2,Cdim)
!-----------------------------------
!   Computes an-isotropic D-matrix
!   Plane-strain (Cdim= 2)
!   or 3-D       (Cdim= 3)
!-----------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(OUT) :: D(:,:)!  D-matrix
REAL(iwp), INTENT(IN)  :: E1    ! Young's modulus, dir 1
REAL(iwp), INTENT(IN)  :: G1    ! Shear modulus , dir 1
REAL(iwp), INTENT(IN)  :: E2    ! Young's modulus, dir 2
REAL(iwp), INTENT(IN)  :: G2    ! Shear Modulus , dir 2
REAL(iwp), INTENT(IN)  :: ny2   ! Poisson's ratio, dir 2
INTEGER,INTENT(IN):: Cdim  !  Cartesian Dimension
REAL(iwp)              :: n     !  ratio E1/E2
REAL(iwp) :: cc,c1,c2,c3,c4,ny1   !   temps
ny1= 0.5*E1/G1 -1.0
n= E1/E2
cc= E2/(1.+ny1)/(1.-ny1-2.*n*ny2**2)
c1= n*(1.-n*ny2**2)*cc
c3= n*ny2*(1.0+ny1)*cc
c4= (1 - ny1**2)*cc
D= 0. !  only nonzero components of D are assigned
SELECT CASE (Cdim)
	CASE (2)          !    plane strain
		D(1,1)= c1 ; D(2,2)= c4
		D(1,2)= c3 ; D(2,1)= c3
		D(3,3)= G2
	CASE (3)          !    3-D
		c2= n*(ny1+n*ny2**2)*cc
		D(1,1)= C1 ; D(2,2)= c1 ; D(3,3)= c4
		D(1,2)= C2 ; D(1,3)= c3 ; D(2,3)= C3
		D(2,1)= C2 ; D(3,1)= c3 ; D(3,2)= C3
		D(4,4)= G1 ; D(5,5)= G2 ; D(6,6)= G2
	CASE DEFAULT
END SELECT
RETURN
END SUBROUTINE D_mat_anis

SUBROUTINE Normal_Jac(v3,Jac,xsi,eta,ldim,nodes,inci,coords)
!--------------------------------------------------------
!   Computes normal vector and Jacobian
!   at point with local coordinates xsi,eta
!--------------------------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(IN):: xsi,eta            ! intrinsic co-ordinates of point
INTEGER,INTENT(IN):: ldim             ! element dimension
INTEGER,INTENT(IN):: nodes            ! number of nodes
INTEGER,INTENT(IN):: inci(:)          ! element incidences
REAL(iwp), INTENT(IN)  :: coords(:,:)      ! node coordinates
REAL(iwp),INTENT(OUT):: v3(:)              ! Vector normal to point
REAL(iwp),INTENT(OUT):: Jac                ! Jacobian
REAL(iwp),ALLOCATABLE  :: DNi(:,:)         ! Derivatives of shape function
REAL(iwp),ALLOCATABLE  :: v1(:),v2(:)      ! Vectors in xsi,eta directions
INTEGER :: Cdim ,i                      ! Cartesian dimension
!    Cartesian dimension:
 Cdim= ldim+1
!    Allocate temporary arrays
ALLOCATE (DNi(nodes,ldim),V1(Cdim),V2(Cdim))
!    Compute derivatives of shape function
Call Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
!    Compute vectors in xsi (eta) direction(s)
DO I=1,Cdim
	V1(I)= DOT_PRODUCT(DNi(:,1),COORDS(I,:))
	IF(ldim == 2) THEN
		V2(I)= DOT_PRODUCT(DNi(:,2),COORDS(I,:))
	END IF
END DO
!    Compute normal vector
IF(ldim == 1) THEN
	v3(1)= V1(2)
	v3(2)= -V1(1)
	ELSE
  V3= Vector_ex(v1,v2)
END IF
!    Normalise
CAll Vector_norm(V3,Jac)
DEALLOCATE (DNi,V1,V2)
RETURN
END SUBROUTINE Normal_Jac

SUBROUTINE Tangent(v1,v2,xsi,eta,ldim,nodes,inci,coords)
!--------------------------------------------------------
!   Computes vectors tangent to BE
!   at point with local coordinates xsi,eta
!--------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN):: ldim             ! element dimension
REAL(iwp), INTENT(IN)  :: xsi,eta     ! intrinsic co-ordinates of point
INTEGER,INTENT(IN):: nodes            ! number of nodes
INTEGER,INTENT(IN):: inci(:)          ! element incidences
REAL(iwp), INTENT(IN)  :: coords(:,:)      ! node coordinates
REAL(iwp),INTENT(OUT)  :: v1(ldim+1),v2(ldim+1)    ! Vector normal to point
REAL(iwp),ALLOCATABLE  :: DNi(:,:)         ! Derivatives of shape function
!      REAL,ALLOCATABLE  :: v1(:),v2(:)      ! Vectors in xsi,eta directions
INTEGER :: Cdim ,i                      ! Cartesian dimension
!    Cartesian dimension:
Cdim= ldim+1
!    Allocate temporary arrays
ALLOCATE (DNi(nodes,ldim))
!    Compute derivatives of shape function
Call Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
!    Compute vectors in xsi (eta) direction(s)
DO i=1,Cdim
	v1(i)= DOT_PRODUCT(DNi(:,1),COORDS(i,:))
	IF(ldim == 2) THEN
		v2(i)= DOT_PRODUCT(DNi(:,2),COORDS(i,:))
	END IF
END DO
DEALLOCATE (DNi)
RETURN
END SUBROUTINE Tangent

Subroutine Serendip_func(Ni,xsi,eta,ldim,nodes,inci)
!---------------------------------
!   Computes Serendipity shape functions  Ni(xsi,eta)
!   for one and two-dimensional (linear/parabolic) finite boundary elements
!---------------------------------
REAL(iwp),INTENT(OUT):: Ni(:)      ! Array with shape function values at xsi,eta
REAL(iwp), INTENT(IN):: xsi,eta    ! intrinsic co-ordinates
INTEGER,INTENT(IN):: ldim     ! element dimension
INTEGER,INTENT(IN):: nodes    ! number of nodes
INTEGER,INTENT(IN):: inci(:)  ! element incidences
REAL(iwp)              :: mxs,pxs,met,pet   !  temporary variables
SELECT CASE (ldim)
	CASE (1)   !   one-dimensional element
		Ni(1)= 0.5*(1.0 - xsi) ;  Ni(2)= 0.5*(1.0 + xsi)
		IF(nodes == 2) RETURN  !  linear element finshed
		Ni(3)=  1.0 - xsi*xsi
		Ni(1)= Ni(1) - 0.5*Ni(3) ; Ni(2)= Ni(2) - 0.5*Ni(3)
	CASE(2)    !    two-dimensional element
		mxs= 1.0-xsi ; pxs= 1.0+xsi ; met= 1.0-eta ; pet= 1.0+eta
		Ni(1)= 0.25*mxs*met ; Ni(2)= 0.25*pxs*met
		Ni(3)= 0.25*pxs*pet ; Ni(4)= 0.25*mxs*pet
    IF(nodes == 4) RETURN   !  linear element finshed
    IF(Inci(5) > 0) THEN        !  zero node number means node is missing
			Ni(5)= 0.5*(1.0 -xsi*xsi)*met
			Ni(1)= Ni(1) - 0.5*Ni(5) ; Ni(2)= Ni(2) - 0.5*Ni(5)
    END IF
    IF(Inci(6) > 0) THEN
			Ni(6)= 0.5*(1.0 -eta*eta)*pxs
			Ni(2)= Ni(2) - 0.5*Ni(6) ;  Ni(3)= Ni(3) - 0.5*Ni(6)
    END IF
    IF(Inci(7) > 0) THEN
			Ni(7)= 0.5*(1.0 -xsi*xsi)*pet
			Ni(3)= Ni(3) - 0.5*Ni(7) ; Ni(4)= Ni(4) - 0.5*Ni(7)
    END IF
    IF(Inci(8) > 0) THEN
			Ni(8)= 0.5*(1.0 -eta*eta)*mxs
			Ni(4)= Ni(4) - 0.5*Ni(8) ; Ni(1)= Ni(1) - 0.5*Ni(8)
    END IF
CASE DEFAULT     !   error message
CALL Error_message('Element dimension not 1 or 2' )
END SELECT
RETURN
END SUBROUTINE Serendip_func

SUBROUTINE Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
!---------------------------------
!   Computes Derivatives ofSerendipity shape functions  Ni(xsi,eta)
!   for one and two-dimensional (linear/parabolic) finite boundary elements
!---------------------------------
IMPLICIT NONE
REAL(iwp),INTENT(OUT):: DNi(:,:) ! Array with shape function derivatives at xsi,eta
REAL(iwp),INTENT(IN):: xsi,eta     ! intrinsic co-ordinates
INTEGER,INTENT(IN):: ldim     ! element dimension
INTEGER,INTENT(IN):: nodes    ! number of nodes
INTEGER,INTENT(IN):: inci(:)  ! element incidences
REAL(iwp)            :: mxs,pxs,met,pet   !  temporary variables
SELECT CASE (ldim)
	CASE (1)   !   one-dimensional element
		DNi(1,1)= -0.5
		DNi(2,1)= 0.5
		IF(nodes == 2) RETURN  !  linear element finshed
		DNi(3,1)=  -2.0*xsi
		DNi(1,1)= DNi(1,1) - 0.5*DNi(3,1)
		DNi(2,1)= DNi(2,1) - 0.5*DNi(3,1)
	CASE (2)   !    two-dimensional element
		mxs= 1.0-xsi
		pxs= 1.0+xsi
		met= 1.0-eta
		pet= 1.0+eta
		DNi(1,1)= -0.25*met
		DNi(1,2)= -0.25*mxs
		DNi(2,1)=  0.25*met
		DNi(2,2)= -0.25*pxs
		DNi(3,1)=  0.25*pet
		DNi(3,2)=  0.25*pxs
		DNi(4,1)= -0.25*pet
		DNi(4,2)=  0.25*mxs
		IF(nodes == 4) RETURN  !  linear element finshed
		IF(Inci(5) > 0) THEN   !  zero node number means node is missing
			DNi(5,1)= -xsi*met
			DNi(5,2)= -0.5*(1.0 -xsi*xsi)
			DNi(1,1)= DNi(1,1) - 0.5*DNi(5,1)
			DNi(1,2)= DNi(1,2) - 0.5*DNi(5,2)
			DNi(2,1)= DNi(2,1) - 0.5*DNi(5,1)
			DNi(2,2)= DNi(2,2) - 0.5*DNi(5,2)
		END IF
		IF(Inci(6) > 0) THEN
			DNi(6,1)= 0.5*(1.0 -eta*eta)
			DNi(6,2)= -eta*pxs
			DNi(2,1)= DNi(2,1) - 0.5*DNi(6,1)
			DNi(2,2)= DNi(2,2) - 0.5*DNi(6,2)
			DNi(3,1)= DNi(3,1) - 0.5*DNi(6,1)
			DNi(3,2)= DNi(3,2) - 0.5*DNi(6,2)
		END IF
		IF(Inci(7) > 0) THEN
			DNi(7,1)= -xsi*pet
			DNi(7,2)= 0.5*(1.0 -xsi*xsi)
			DNi(3,1)= DNi(3,1) - 0.5*DNi(7,1)
			DNi(3,2)= DNi(3,2) - 0.5*DNi(7,2)
			DNi(4,1)= DNi(4,1) - 0.5*DNi(7,1)
			DNi(4,2)= DNi(4,2) - 0.5*DNi(7,2)
		END IF
		IF(Inci(8) > 0) THEN
			DNi(8,1)= -0.5*(1.0 -eta*eta)
			DNi(8,2)= -eta*mxs
			DNi(4,1)= DNi(4,1) - 0.5*DNi(8,1)
			DNi(4,2)= DNi(4,2) - 0.5*DNi(8,2)
			DNi(1,1)= DNi(1,1) - 0.5*DNi(8,1)
			DNi(1,2)= DNi(1,2) - 0.5*DNi(8,2)
		END IF
	CASE DEFAULT     !   error message
CALL Error_message('Element dimension not 1 or 2' )
END SELECT
RETURN
END SUBROUTINE Serendip_deriv
 
SUBROUTINE Cartesian(Ccor,Ni,ldim,elcor)
!--------------------------------------------------------
!   Computes Cartesian coordinates
!   at point with local coordinates xsi,eta
!--------------------------------------------------------
IMPLICIT NONE
REAL(iwp),INTENT(OUT)  :: Ccor(:)            ! Cart. coords of point xsi,eta
REAL(iwp),INTENT(IN)   :: Ni(:)              ! Shape functions at xsi,eta
REAL(iwp),INTENT(IN)  :: elcor(:,:)       ! element coordinates
INTEGER           :: ldim, Cdim,I             ! Cartesian dimension
!    Cartesian dimension:
Cdim= ldim+1
!    Compute vectors in xsi (eta) direction(s)
DO I=1,Cdim
	Ccor(I)= DOT_PRODUCT(Ni(:),Elcor(I,:))
END DO
RETURN
END SUBROUTINE Cartesian

SUBROUTINE Vector_norm(v3,Vlen)
!----------------------------------------
!   Normalise vector
!----------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(INOUT)  :: V3(:)               !     Vector to be normalised
REAL(iwp), INTENT(OUT)    :: Vlen                !     length of vector
Vlen= SQRT( SUM(v3*v3))
IF(Vlen == 0.) RETURN
V3= V3/Vlen
RETURN
END SUBROUTINE Vector_norm

FUNCTION Vector_ex(v1,v2)
!----------------------------------------
!   Returns vector x-product v1xv2
!   Where v1 and v2 are dimension 3
!----------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(IN) :: V1(3),V2(3)               !     Input
REAL(iwp)             :: Vector_ex(3)              !     Result
Vector_ex(1)=V1(2)*V2(3)-V2(2)*V1(3)
Vector_ex(2)=V1(3)*V2(1)-V1(1)*V2(3)
Vector_ex(3)=V1(1)*V2(2)-V1(2)*V2(1)
RETURN
END FUNCTION Vector_ex

REAL FUNCTION Dist(xa,xe,Cdim)
!----------------------------------------
!    Computes the distance between two points
!    with coordinates (xa,ya) and (xe,ye)
!------------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(IN)    :: xa(:)      !  coords of point 1
REAL(iwp), INTENT(IN)    :: xe(:)      !  coords of point 2
INTEGER, INTENT(IN) :: Cdim       !  Cartesian dimension
INTEGER             :: N
REAL(iwp)           :: SUMS
SUMS= 0.0_iwp
DO N=1,Cdim
	SUMS= SUMS + (xa(n)-xe(n))**2
END DO
Dist= SQRT(SUMS)
RETURN
END FUNCTION Dist

REAL FUNCTION Direc(xA,xE)
!--------------------------------------------------------
!  Computes the Direction-angle from point xA to point xE
!--------------------------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(IN)     :: xA(2)
REAL(iwp), INTENT(IN)     :: xE(2)
REAL(iwp)                 :: pi=3.1415926536_iwp
Direc=ATAN2((xE(2)-xA(2)),(xE(1)-xA(1)))
IF (Direc < 0.00000000) THEN
	Direc= Direc + 2*pi
END IF
RETURN
END FUNCTION Direc

SUBROUTINE Ortho(v3,v1,v2)
!----------------------------------------------------------
!     DETERMINES ORTHOGONAL VECTORS
!     V1, V2, V3
!     USING
!     V2 = V3 X VX
!     OR
!     V2 = V3 X VY     IF V3=VX
!     AND
!     V1 = V2 X V3
!--------------------------------------------------------
REAL(iwp), INTENT (IN) ::  v3(:)         !  Normal vector
REAL(iwp)              ::  V1(:), V2(:)  !  Orthogonal vectors
REAL(iwp)              ::  vx(3),vy(3)   !  Vectors in coordinate directions
vx= (/1.0,0.0,0.0/) ; vy= (/0.0,1.0,0.0/)
IF(ABS(V3(1)) + 0.005 .GE. 1.0) THEN
	v2= Vector_ex(v3,vy)
ELSE
	v2= Vector_ex(v3,vx)
END IF
V1= Vector_ex(v2,v3)
RETURN
END SUBROUTINE Ortho

REAL FUNCTION Min_dist(Elcor,xPi,Nodel,ldim,inci)
IMPLICIT NONE
REAL(iwp),INTENT(IN)   ::      Elcor(:,:)    !       Coordinates of Element
REAL(iwp),INTENT(IN)   ::      xPi(:)    !       Coordinates of Collocation point
REAL(iwp)                     ::      DET,A,B,C,D 
REAL(iwp)                     ::      xsi,eta,Dxsi,Deta
REAL(iwp)                     ::      DistPS,DistPS_N
REAL(iwp)                     ::      L,ELengx,ELenge
REAL(iwp),ALLOCATABLE         ::      Ni(:)   ! Shape function
REAL(iwp),ALLOCATABLE         ::      DNi(:,:)! Derivatives of shape function
REAL(iwp),ALLOCATABLE         ::      r(:),xS(:),Dxs(:,:)
INTEGER,INTENT(IN)            ::      Nodel,ldim
INTEGER,INTENT(IN)            ::      inci(:)
INTEGER                       ::      n,Cdim

Cdim= ldim + 1
ALLOCATE(Ni(Nodel),DNi(Nodel,ldim),r(Cdim),xS(Cdim),DxS(Cdim,ldim))
SELECT CASE(Cdim)
	CASE(2)
		xsi=0.0
		CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
		xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
		xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
		r= xS-xPi
		CALL Elength(L,Elcor,Nodel,ldim)
		DistPS=Dist(xPi(:),xS(:),Cdim)
		IF(((DistPS-L/2)/L) > 4.)THEN
			Min_dist=DistPS
			RETURN
		END IF
		DO n=1,20
			IF(n > 1)DistPS= DistPS_N
			CALL Serendip_deriv(DNi,xsi,eta,ldim,nodel,inci)
			DxS(1,1)= DOT_PRODUCT(DNi(:,1),Elcor(1,:))
			DxS(2,1)= DOT_PRODUCT(DNi(:,1),Elcor(2,:))
			DET= DxS(1,1)**2+DxS(2,1)**2
			Dxsi= -1/DET * DOT_PRODUCT(DxS(:,1),r(:))
			xsi= xsi+ Dxsi
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
			xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
			xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
			r= xS- xPi
			DistPS_N=Dist(xPi(:),xS(:),Cdim)
			IF(ABS((DistPS- DistPS_N)/DistPS_N) < 0.05)EXIT
		END DO
		IF(xsi >= -1.0 .and. xsi <= 1.0)THEN
			Min_dist=DistPS_N
			RETURN
		ELSE IF(xsi < -1.0)THEN
			Min_dist=Dist(Elcor(:,1),xPi(:),Cdim)
			RETURN
		ELSE IF(xsi > -1.0)THEN
			Min_dist=Dist(Elcor(:,2),xPi(:),Cdim)
			RETURN
		END IF
	CASE(3)
		xsi=0.0
		eta=0.0
		CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
		xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
		xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
		xS(3)= DOT_PRODUCT(Ni(:),Elcor(3,:))
		r= xS-xPi
  ELengx= Dist((Elcor(:,3)+Elcor(:,2))/2.,(Elcor(:,4)+Elcor(:,1))/2.,Cdim) ! Lxsi
  ELenge= Dist((Elcor(:,2)+Elcor(:,1))/2.,(Elcor(:,3)+Elcor(:,4))/2.,Cdim) ! Leta
		DistPS=Dist(xPi(:),xS(:),Cdim)
  IF(((DistPS-ELengx/2)/ELengx) > 4. .and.((DistPS-ELenge/2)/ELenge)> 4.)THEN
			Min_dist=DistPS
			RETURN
		END IF
		DO n=1,40
			IF(n > 1)DistPS= DistPS_N
			CALL Serendip_deriv(DNi,xsi,eta,ldim,nodel,inci)
			DxS(1,1)= DOT_PRODUCT(DNi(:,1),Elcor(1,:))
			DxS(2,1)= DOT_PRODUCT(DNi(:,1),Elcor(2,:))
			DxS(3,1)= DOT_PRODUCT(DNi(:,1),Elcor(3,:))
			DxS(1,2)= DOT_PRODUCT(DNi(:,2),Elcor(1,:))
			DxS(2,2)= DOT_PRODUCT(DNi(:,2),Elcor(2,:))
			DxS(3,2)= DOT_PRODUCT(DNi(:,2),Elcor(3,:))
			A= DxS(1,1)**2+DxS(2,1)**2+DxS(3,1)**2
                   B= DxS(1,1)*DxS(1,2)+DxS(2,1)*DxS(2,2)+DxS(3,1)*DxS(3,2)
			C=B
			D= DxS(1,2)**2+DxS(2,2)**2+DxS(3,2)**2
			DET= A*D - C*B
			Dxsi = -1/DET * DOT_PRODUCT(DxS(:,1),r(:))
			Deta = -1/DET * DOT_PRODUCT(DxS(:,2),r(:))
			xsi= xsi+ Dxsi
			eta= eta+ Deta
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
			xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
			xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
			xS(3)= DOT_PRODUCT(Ni(:),Elcor(3,:))
			r= xS- xPi
			DistPS_N=Dist(xPi(:),xS(:),Cdim)
			IF(ABS((DistPS- DistPS_N)/DistPS_N) < 0.01)EXIT
		END DO
  IF(xsi >= -1.0 .and. xsi <= 1.0 .and. eta >= -1.0 .and. eta <= 1.0)THEN
			Min_dist=DistPS_N
			RETURN
		ELSE IF(xsi < -1.0 .and. eta < -1.0)THEN
			Min_dist=Dist(Elcor(:,1),xPi(:),Cdim)
			RETURN
		ELSE IF(xsi > 1.0 .and. eta > 1.0)THEN
			Min_dist=Dist(Elcor(:,3),xPi(:),Cdim)
			RETURN
		ELSE IF(xsi > 1.0 .and. eta < -1.0)THEN
			Min_dist=Dist(Elcor(:,2),xPi(:),Cdim)
			RETURN
		ELSE IF(xsi < -1.0 .and. eta > 1.0)THEN
			Min_dist=Dist(Elcor(:,4),xPi(:),Cdim)
			RETURN
		ELSE IF(xsi >= -1.0 .and. xsi <= 1.0 .and. eta < -1.0)THEN
			eta=-1.0
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
			xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
			xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
			xS(3)= DOT_PRODUCT(Ni(:),Elcor(3,:))
			Min_dist=Dist(xS(:),xPi(:),Cdim)
			RETURN
		ELSE IF(xsi >= -1.0 .and. xsi <= 1.0 .and. eta > 1.0)THEN
			eta=1.0
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
			xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
			xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
			xS(3)= DOT_PRODUCT(Ni(:),Elcor(3,:))
			Min_dist=Dist(xS(:),xPi(:),Cdim)
			RETURN
		ELSE IF(eta >= -1.0 .and. eta <= 1.0 .and. xsi > 1.0)THEN
			xsi=1.0
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
			xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
			xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
			xS(3)= DOT_PRODUCT(Ni(:),Elcor(3,:))
			Min_dist=Dist(xS(:),xPi(:),Cdim)
			RETURN
		ELSE IF(eta >= -1.0 .and. eta <= 1.0 .and. xsi < -1.0)THEN
			xsi=-1.0
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
			xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
			xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
			xS(3)= DOT_PRODUCT(Ni(:),Elcor(3,:))
			Min_dist=Dist(xS(:),xPi(:),Cdim)
			RETURN
		END IF
END SELECT

END FUNCTION Min_dist

REAL FUNCTION Min_dist1(Elcor,xPi,Nodel,inci,ELengx,Elenge,ldim)
IMPLICIT NONE
REAL(iwp),INTENT(IN)   ::   Elcor(:,:)      !       Coordinates of Element
REAL(iwp),INTENT(IN)   ::   xPi(:)     !   Coordinates of Collocation point
REAL(iwp),INTENT(IN)   ::   ELengx,ELenge       !   Elementlength xsi and eta
REAL(iwp)              ::      DET,A,B,C,D,F1,F2 
REAL(iwp)              ::      xsi,eta,Dxsi,Deta
REAL(iwp)              ::      DistPS,DistPS_N
REAL(iwp)              ::      L
REAL(iwp)              ::      ERR
REAL(iwp),ALLOCATABLE  ::    Ni(:)                  ! Shape function
REAL(iwp),ALLOCATABLE        ::      DNi(:,:) ! Derivatives of shape function
REAL(iwp),ALLOCATABLE        ::      r(:),xS(:),Dxs(:,:)
INTEGER,INTENT(IN)::	Nodel
INTEGER,INTENT(IN)::	inci(:)
INTEGER,INTENT(IN)::	ldim
INTEGER           ::      n,cdim
cdim= ldim + 1
ALLOCATE(Ni(Nodel),DNi(Nodel,ldim),r(Cdim),xS(Cdim),DxS(Cdim,ldim))
SELECT CASE(Cdim)
  CASE(2)
	  xsi=0.0
		CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)  
                CALL Cartesian(xS,Ni,ldim,Elcor)
		r= xS-xPi
		DistPS=Dist(xPi(:),xS(:),Cdim)
		IF(((DistPS-ELengx/2)/Elengx) > 4.)THEN
			Min_dist1=DistPS-Elengx/2
			RETURN
		END IF
    DO n=1,40
      IF(n > 1)DistPS= DistPS_N
      CALL Serendip_deriv(DNi,xsi,eta,ldim,nodel,inci)
      DxS(1,1)= DOT_PRODUCT(DNi(:,1),Elcor(1,:))
      DxS(2,1)= DOT_PRODUCT(DNi(:,1),Elcor(2,:))
      DET= DxS(1,1)**2+DxS(2,1)**2
      Dxsi= -1/DET * DOT_PRODUCT(DxS(:,1),r(:))
      xsi= xsi+ Dxsi
      IF(ABS(xsi) > 1.) xsi= xsi/ABS(xsi)
      CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
      CALL Cartesian(xS,Ni,ldim,Elcor)
      r= xS- xPi
      DistPS_N= Dist(xPi(:),xS(:),Cdim)
      IF(DistPS_N > DistPS)THEN
!                xsi=xsi- Dxsi
        Min_dist1= DistPS
        RETURN
      END IF
      ERR= (DistPS- DistPS_N)/DistPS_N
      IF(ERR < 0.05)THEN 
        Min_dist1= DistPS_N
!             IF(numpe==1) WRITE(12,*)'n=',n
        RETURN
      END IF
		END DO
	CASE(3)
		xsi=0.0
		eta=0.0
		CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
    CALL Cartesian(xS,Ni,ldim,Elcor)
		r= xS-xPi
		DistPS=Dist(xPi(:),xS(:),Cdim)
 IF(((DistPS-ELengx/2)/ELengx) > 4. .and.((DistPS-ELenge/2)/ELenge)> 4.)THEN
			Min_dist1=DistPS
			RETURN
		END IF
    DO n=1,40
			IF(n > 1)DistPS= DistPS_N
			CALL Serendip_deriv(DNi,xsi,eta,ldim,nodel,inci)
			DxS(1,1)= DOT_PRODUCT(DNi(:,1),Elcor(1,:))
			DxS(2,1)= DOT_PRODUCT(DNi(:,1),Elcor(2,:))
			DxS(3,1)= DOT_PRODUCT(DNi(:,1),Elcor(3,:))
			DxS(1,2)= DOT_PRODUCT(DNi(:,2),Elcor(1,:))
			DxS(2,2)= DOT_PRODUCT(DNi(:,2),Elcor(2,:))
			DxS(3,2)= DOT_PRODUCT(DNi(:,2),Elcor(3,:))
			A= DxS(1,1)**2+DxS(2,1)**2+DxS(3,1)**2
                B= DxS(1,1)*DxS(1,2)+DxS(2,1)*DxS(2,2)+DxS(3,1)*DxS(3,2)
			C=B
			D= DxS(1,2)**2+DxS(2,2)**2+DxS(3,2)**2
			DET= A*D - C*B
			F1= DOT_PRODUCT(DxS(:,1),r(:))
      F2= DOT_PRODUCT(DxS(:,2),r(:))
      Dxsi = -1/DET * (F1*D - F2*B)
			Deta = -1/DET * (F2*A - F1*C)
			xsi= xsi+ Dxsi
			eta= eta+ Deta
      IF(ABS(xsi) > 1.) xsi= xsi/ABS(xsi)
      IF(ABS(eta) > 1.) eta= eta/ABS(eta)
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
      CALL Cartesian(xS,Ni,ldim,Elcor)
			r= xS- xPi
			DistPS_N=Dist(xPi(:),xS(:),Cdim)
      IF(DistPS_N > DistPS)THEN
      xsi=xsi- Dxsi
      eta=eta- Deta
      Min_dist1= DistPS
!                 IF(numpe==1)WRITE(12,*)'n=',n
!     IF(numpe==1)WRITE(12,*)'XSI=',xsi
!     IF(numpe==1)WRITE(12,*)'ETA=',eta
      RETURN
    END IF
    ERR= (DistPS- DistPS_N)/DistPS_N
    IF(ERR < 0.05)THEN
      Min_dist1= DistPS_N
  !   IF(numpe==1)WRITE(12,*)'n=',n
  !             IF(numpe==1)WRITE(12,*)'XSI=',xsi
  !             IF(numpe==1)WRITE(12,*)'ETA=',eta
      RETURN
    END IF
	END DO
END SELECT
DEALLOCATE(Ni,DNi,r,xS,DxS)
END FUNCTION Min_dist1


SUBROUTINE Elength(L,Elcor,nodes,ldim)
!------------------------------------------------
!   Computes the length of a boundary element
!----------------------------------------------
IMPLICIT NONE
REAL(iwp),INTENT (IN)   ::      Elcor(:,:)
INTEGER, INTENT (IN)	::	nodes, ldim
REAL(iwp),INTENT (OUT)  ::      L
REAL(iwp)               ::      B(ldim+1), distB3, distB2, p, a, c
INTEGER                 ::      Cdim
Cdim=ldim+1
SELECT CASE (Nodes)
	CASE (2)
		L=Dist(Elcor(:,1),Elcor(:,2),Cdim)
		RETURN
	CASE (3)
		B=(Elcor(:,1)+Elcor(:,2))/2.0
		distB3=Dist(B(:),Elcor(:,3),Cdim)
		distB2=Dist(B,Elcor(:,2),Cdim)
		IF (distB3/distB2 < 0.01) THEN
			L=Dist(Elcor(:,1),Elcor(:,2),Cdim)
			RETURN
		END IF
		IF (distB3/distB2 < 0.1) THEN
			c=distB3/distB2
   L=2.0*distB2*(1.0+2.0/3.0*c**2-2.0/5.0*c**4) !Length Parabel linearisiert
			RETURN
		END IF
                p=distB2**2/(2*DistB3)          !Parabel Parameter p=y**2/(2*x)
		a=SQRT(p**2+distB2**2)
   L=2.0*(distB2/(2.0*p)*a + p/2.0*LOG((distB2 + a)/p)) !Length Parabel exakt
		RETURN
	CASE DEFAULT
END SELECT
END SUBROUTINE Elength


SUBROUTINE Gauss_coor(Cor,Wi,Intord)
!------------------------------------
! 	 Returns Gauss coordinates and Weights for up to 8 Gauss points
!------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(OUT)  :: Cor(8)     !       Gauss point coordinates
REAL(iwp), INTENT(OUT)  :: Wi(8)      !       weigths
INTEGER,INTENT(IN) :: Intord	 !	 integration order
SELECT CASE (Intord)
	CASE (1)
		Cor(1)= 0.
		Wi(1) = 2.0
	CASE(2)
		Cor(1)= .577350269	; Cor(2)= -Cor(1)
		Wi(1) = 1.0 ;  Wi(2) = Wi(1)
	CASE(3)
		Cor(1)= .774596669	; Cor(2)= 0.0 ; Cor(3)= -Cor(1)
		Wi(1) = .555555555	; Wi(2) = .888888888 ; Wi(3) = Wi(1)
	CASE(4)
  Cor(1)= .861136311 ; Cor(2)= .339981043 ; Cor(3)= -Cor(2) ; Cor(4)= -Cor(1)
  Wi(1) = .347854845 ; Wi(2) = .652145154 ; Wi(3) = Wi(2) ; Wi(4) = Wi(1)
	CASE(5)
  Cor(1)= .9061798459 ; Cor(2)= .5384693101 ; Cor(3)= .0 ; Cor(4)= -Cor(2)
		Cor(5)= -Cor(1)
 Wi(1) = .236926885 ; Wi(2) = .478628670 ; Wi(3) = .568888888 ; Wi(4) = Wi(2)
		Wi(5) = Wi(1)
	CASE(6)
		Cor(1)= .932469514 ; Cor(2)= .661209386 ; Cor(3)= .238619186
		Cor(4)= -Cor(3) ;  Cor(5)= -Cor(2) ; Cor(6)= -Cor(1)
		Wi(1) = .171324492 ; Wi(2) = .360761573 ; Wi(3) = .467913934
		Wi(4) = Wi(3) ; Wi(5) = Wi(2) ; Wi(6) = Wi(1)
	CASE(7)
		Cor(1)= .949107912 ; Cor(2)= .741531185 ; Cor(3)= .405845151
		Cor(4)= 0.
		Cor(5)= -Cor(3) ;Cor(6)= -Cor(2) ;Cor(7)= -Cor(1)
		Wi(1) = .129484966 ; Wi(2) = .279705391 ; Wi(3) = .381830050
		Wi(4) = .417959183
		Wi(5) = Wi(3) ; Wi(6) = Wi(2) ; Wi(7) = Wi(1)
	CASE(8)
                Cor(1)= .960289856 ; Cor(2)= .796666477
                Cor(3)= .525532409 ; Cor(4)= .183434642
                Cor(5)= -Cor(4) ; Cor(6)= -Cor(3)
                Cor(7)= -Cor(2) ; Cor(8)= -Cor(1)
                Wi(1) = .101228536 ; Wi(2) = .222381034
                Wi(3) = .313706645 ;Wi(4) = .362683783
		Wi(5) = Wi(4) ; Wi(6) = Wi(3) ; Wi(7) = Wi(2) ; Wi(8) = Wi(1)
	CASE DEFAULT
CALL Error_Message('Gauss points not in range 1-8')
END SELECT
END SUBROUTINE Gauss_coor

SUBROUTINE Gauss_Laguerre_coor(Cor,Wi,Intord)
!------------------------------------
! 	Returns Gauss_Laguerre coordinates and Weights for up to 8 Gauss points
!------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(OUT)  :: Cor(8)     !       Gauss point coordinates
REAL(iwp), INTENT(OUT)  :: Wi(8)      !       weigths
INTEGER,INTENT(IN) :: Intord	 !	 integration order
SELECT CASE (Intord)
	CASE (1)
		Cor(1)= 0.5
		Wi(1) = 1.0
	CASE(2)
		Cor(1)= .112008806	; Cor(2)=.602276908
		Wi(1) = .718539319 ;	Wi(2) =.281460680
	CASE(3)
                Cor(1)= .063890793  ; Cor(2)= .368997063 ; Cor(3)= .766880303
                Wi(1) = .513404552  ; Wi(2) = .391980041 ; Wi(3) = .0946154065
	CASE(4)
                Cor(1)= .0414484801 ; Cor(2)=.245274914
                Cor(3)=.556165453 ; Cor(4)= .848982394
                Wi(1) = .383464068      ; Wi(2) =.386875317
                Wi(3) =.190435126 ; Wi(4) = .0392254871
	CASE(5)
                Cor(1)= .0291344721 ; Cor(2)= .173977213
                Cor(3)= .411702520; Cor(4)=.677314174
		Cor(5)= .894771361
                Wi(1) = .297893471      ; Wi(2) = .349776226
                Wi(3) =.234488290 ; Wi(4) = .0989304595
		Wi(5) = .0189115521
	CASE(6)
                Cor(1)= .0216340058 ; Cor(2)= .129583391
                Cor(3)= .314020449
                Cor(4)= .538657217      ; Cor(5)= .756915337
                Cor(6)=.922668851
		Wi(1) =  .238763662 ; Wi(2) =.308286573 ; Wi(3) =.245317426
		Wi(4) = .142008756 ; Wi(5) =.0554546223 ; Wi(6) =.0101689586
	CASE(7)
		Cor(1)= .0167193554 ; Cor(2)= .100185677 ; Cor(3)= .246294246
		Cor(4)= .433463493
                Cor(5)= .632350988  ; Cor(6)= .811118626 ; Cor(7)= .940848166
                Wi(1) = .196169389  ; Wi(2) = .270302644 ; Wi(3) = .239681873
		Wi(4) = .165775774
                Wi(5) = .0889432271 ; Wi(6) =.0331943043 ; Wi(7)=.00593278701
	CASE(8)
                Cor(1)= .0133202441 ; Cor(2)=.0797504290 ; Cor(3)= .197871029
                Cor(4)= .354153994
                Cor(5)= .529458575      ; Cor(6)= .701814529
                Cor(7)= .849379320 ; Cor(8)= .953326450
                Wi(1) = .164416604 ; Wi(2) = .237525610
                Wi(3) = .226841984 ;Wi(4) = .175754079
                Wi(5) = .112924030 ; Wi(6) =.0578722107
                Wi(7) =.0209790737 ;Wi(8) =.00368640710
	CASE DEFAULT
	CALL Error_Message('Gauss points not in range 1-8')
END SELECT
END SUBROUTINE Gauss_Laguerre_coor

INTEGER FUNCTION Ngaus(RonL,ne)
!-----------------------------------------------------
! 	 Function returns number of Gauss points needed 
! 	 to integrate a function 1/rn
!------------------------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(IN)::      RonL                    !       R/L
INTEGER , INTENT(IN) ::      ne            !       exponent (1,2,3)
REAL(iwp)            ::      Rlim(7)  !       array to store values of table
INTEGER              ::      n
SELECT CASE(ne)
	CASE(1)
         Rlim= (/1.6382, 0.6461, 0.3550, 0.2230, 0.1490, 0.1021, 0.0698/)
	CASE(2)
         Rlim= (/2.6230, 1.0276, 0.5779, 0.3783, 0.2679, 0.1986, 0.1512/)
	CASE(3)
         Rlim= (/3.5627, 1.3857, 0.7846, 0.5212, 0.3767, 0.2864, 0.2249/)
	CASE DEFAULT
END SELECT
Ngaus=0
DO	 N=1,7
	IF(RonL >= Rlim(N)) THEN
		Ngaus= N+1
		EXIT
	END IF
END DO 
IF (Ngaus == 0)THEN                   !       Point is too near the surface
	Ngaus=8
END IF
RETURN
END FUNCTION Ngaus

SUBROUTINE Integ2P (Elcor,Inci,Nodel,Ncol,xP,k,dUe,dTe,Ndest,Isym)
!--------------------------------------------------
! 	 Computes  [dT]e and [dU]e for 2-D potential problems
! 	 by numerical integration
!-------------------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(IN)    :: Elcor(:,:)     !       Element coordinates
INTEGER, INTENT(IN) :: Ndest(:,:)     !   Node destination vector
INTEGER, INTENT(IN) :: Inci(:)        !       Element Incidences
INTEGER, INTENT(IN) :: Nodel          !       No. of Element Nodes
INTEGER , INTENT(IN):: Ncol           !       Number of points
                                      ! Pi (coll. points)
INTEGER , INTENT(IN):: Isym
REAL(iwp) , INTENT(IN)       :: xP(:,:)    !  Array with coll. points coords.
REAL(iwp) , INTENT(IN)       :: k          ! Permeability
REAL(iwp) , INTENT(OUT)      :: dUe(:,:),dTe(:,:)    !       arrays
REAL(iwp) :: epsi= 1.0E-4_iwp        !        Small value for comparing coords
REAL(iwp) :: Eleng,Rmin,RonL,Glcor(8),Wi(8),Ni(Nodel),Vnorm(2),GCcor(2)
REAL(iwp) :: UP,Jac,dxr(2),TP,r,pi,c1,c2,xsi,eta,dxdxb
INTEGER :: i,m,n,Mi,nr,ldim,cdim,nreg,id
pi=3.14159265359_iwp
ldim= 1
cdim=ldim+1
CALL Elength(Eleng,Elcor,nodel,ldim)			! Element Length
dUe= 0.0 ; dTe= 0.0                   ! Clear arrays for summation
!-----------------------------------------------------------------
!         Integration off-diagonal coeff.  -> normal Gauss Quadrature
!-----------------------------------------------------------------
dUe= 0.0 ; dTe= 0.0         ! Clear arrays for summation
Colloc_points:	DO i=1,Ncol
!    (Elcor,xPi,Nodel,inci,ELengx,Elenge,ldim)
     Rmin= Min_dist1(Elcor,xP(:,i),Nodel,inci,ELeng,Eleng,ldim)
! Distance coll. point and element
     RonL= Rmin/Eleng                         !       R/L
     Mi= Ngaus(RonL,1)  !       Number of Gauss points for (1/r) singularity
     Call Gauss_coor(Glcor,Wi,Mi)            ! Assign coords/Weights
         Gauss_points:   DO m=1,Mi
             xsi= Glcor(m)
             CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)
                     !       Shape function value
             Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor)
                   ! Jacobian and normal
             CALL Cartesian(GCcor,Ni,ldim,elcor)  !  Cart. coords of Gauss pt
             r= Dist(GCcor,xP(:,i),cdim)          !     Dist. P,Q
             dxr= (GCcor-xP(:,i))/r               !       rx/r , ry/r
             UP= U(r,k,cdim) ; TP= T(r,dxr,Vnorm,cdim)    ! Kernels
               Node_points:    DO n=1,Nodel
                                  IF(Isym == 0)THEN
                                  iD= i
                                   ELSE
                                  iD= Ndest(i,1)    !  line number in array
                                  END IF
                                  IF (id == 0) CYCLE
                                  IF(Dist(Elcor(:,n),xP(:,i),cdim)>epsi) THEN
                                ! Only case where coords of n and Pi not same
                                  dUe(id,n)= dUe(id,n) + Ni(n)*UP*Jac*Wi(m)
                                  dTe(id,n)= dTe(id,n) + Ni(n)*TP*Jac*Wi(m)
                                  END IF
                END DO Node_points
        END DO Gauss_points
END DO Colloc_points
!------------------------------
! 		 Diagonal terms of dUe
!------------------------------
c1= 1/(2.0*pi*k)
Colloc_points1:	DO i=1,Ncol
		Node_points1: DO n=1,Nodel
                   IF(Isym == 0)THEN
                     iD= i
                   ELSE
                     iD= Ndest(i,1)              !  line number in array
                   END IF
                   IF (id == 0) CYCLE
                   IF(Dist(Elcor(:,n),xP(:,i),cdim) > Epsi) CYCLE
                    ! only do this when Pi is node n
                         Nreg=1
                         IF(n == 3) nreg= 2
                           Subregions: DO nr=1,Nreg
                                       Mi= 4
                                       Call Gauss_Laguerre_coor(Glcor,Wi,Mi)
			Gauss_points1:	DO m=1,Mi
                            SELECT CASE (n)
                             CASE (1)
                              xsi= 2.0*Glcor(m)-1.0
                              dxdxb= 2.0
                             CASE (2)
                              xsi= 1.0 -2.0*Glcor(m)
                              dxdxb= 2.0
                             CASE (3)
                              dxdxb= 1.0
                              IF(nr == 1) THEN
                               xsi= -Glcor(m)
                               ELSE
                               xsi= Glcor(m)
                              END IF
                             CASE DEFAULT
                             END SELECT
                         CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)
                              !       Shape function value
                    Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor)
                     ! Jacobian
                        dUe(id,n)= dUe(id,n) + Ni(n)*c1*Jac*dxdxb*Wi(m)
                        END DO Gauss_points1
                           END DO Subregions

                            Mi= 2  
                            Call Gauss_coor(Glcor,Wi,Mi)
                     ! Assign coords/Weights
		Gauss_points2:	DO m=1,Mi
                      SELECT CASE (n)
                       CASE (1)
                       c2=-LOG(Eleng)*c1
                       CASE (2)
                       c2=-LOG(Eleng)*c1
                       CASE (3)
                       c2=LOG(2/Eleng)*c1
                       CASE DEFAULT
                      END SELECT
                       xsi= Glcor(m)
                       CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)
                            !       Shape function value
                 Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor)
                        ! Jacobian and normal
                         dUe(id,n)= dUe(id,n) + Ni(n)*c2*Jac*Wi(m)
                 END DO Gauss_points2
            END DO Node_points1
         END DO Colloc_points1
    RETURN
END SUBROUTINE Integ2P

SUBROUTINE Integ2E(Elcor,Inci,Nodel,Ncol,xP,E,ny,dUe,dTe,Ndest,Isym)
!--------------------------------------------------
!    Computes  [dT]e and [dU]e for 2-D elasticity problems
!    by numerical integration
!-------------------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(IN)    :: Elcor(:,:)     !   Element coordinates
INTEGER, INTENT(IN) :: Ndest(:,:)     !   Node destination vector
INTEGER, INTENT(IN) :: Inci(:)        !   Element Incidences
INTEGER, INTENT(IN) :: Nodel          !   No. of Element Nodes
INTEGER , INTENT(IN):: Ncol           !   Number of points Pi (coll. points)
INTEGER , INTENT(IN):: Isym  
REAL(iwp), INTENT(IN)                :: E,ny           !   Elastic constants
REAL(iwp), INTENT(IN)    :: xP(:,:)        !   Array with coll. points coords.
REAL(iwp), INTENT(OUT)   :: dUe(:,:),dTe(:,:)
 !  arrays for storing element coefficients
REAL(iwp)    :: epsi= 1.0E-4_iwp  !    Small value for comparing coords
REAL(iwp)    :: Eleng,Rmin,RonL,Glcor(8),Wi(8),Ni(Nodel),Vnorm(2),GCcor(2)
REAL(iwp)    :: Jac,dxr(2),UP(2,2),TP(2,2), xsi, eta, r, dxdxb,Pi,C,C1
INTEGER   :: i,j,k,m,n,Mi,nr,ldim,cdim,iD,nD,Nreg
Pi=3.14159265359_iwp
C=(1.0+ny)/(4*Pi*E*(1.0-ny))			
ldim= 1                             ! Element dimension
cdim=ldim+1
CALL Elength(Eleng,Elcor,nodel,ldim)  ! Element Length
dUe= 0.0_iwp ; dTe= 0.0_iwp                 ! Clear arrays for summation
 Colloc_points: DO i=1,Ncol
         Rmin= Min_dist1(Elcor,xP(:,i),Nodel,inci,ELeng,Eleng,ldim)
          !  Distance coll. point and element
         RonL= Rmin/Eleng                   !  R/L
	!     Integration off-diagonal coeff.  -> normal Gauss Quadrature
           Mi= Ngaus(RonL,1)
              !  Number of Gauss points for (1/r) singularity
           Call Gauss_coor(Glcor,Wi,Mi)  ! Assign coords/Weights
		Gauss_points: DO m=1,Mi
                   xsi= Glcor(m)
                   CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)
                    !   Shape function value
                   Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor)
                    ! Jacobian and normal
                   CALL Cartesian(GCcor,Ni,ldim,elcor) ! Cart. coords of Gauss pt
                      r= Dist(GCcor,xP(:,i),cdim)               !  Dist. P,Q
                      dxr= (GCcor-xP(:,i))/r         !  rx/r , ry/r
                      UP= UK(dxr,r,E,ny,Cdim) ; TP= TK(dxr,r,Vnorm,ny,Cdim)
                       !  Kernels
			Node_points:	DO n=1,Nodel
				Direction_P:	DO j=1,2
                                  IF(Isym == 0)THEN
                                  iD= 2*(i-1) + j
                                  ELSE
                                  iD= Ndest(i,j)          !  line number in array
                                  END IF
                                  IF (id == 0) CYCLE
                            Direction_Q:    DO k= 1,2
                               nD= 2*(n-1) + k           !  column number in array
                               IF(Dist(Elcor(:,n),xP(:,i),cdim) > epsi) THEN
                                ! n and Pi not same
                              dUe(iD,nD)= dUe(iD,nD) + Ni(n)*UP(j,k)*Jac*Wi(m)
                              dTe(iD,nD)= dTe(iD,nD) + Ni(n)*TP(j,k)*Jac*Wi(m)
                               ELSE
                 dUe(iD,nD)= dUe(iD,nD) + Ni(n)*C*dxr(j)*dxr(k)*Jac*Wi(m)
                               !  non-log part of U
                               END IF
                            END DO Direction_Q
                              END DO Direction_P
                         END DO Node_points
                      END DO Gauss_points
                  END DO Colloc_points
!         IF(numpe==1)WRITE (12,*)'Ncol ',Ncol
!         IF(numpe==1)WRITE (12,*)DTe(:,1:6)
!      Diagonal terms of dUe

    	C= C*(3.0-4.0*ny)								
	Colloc_points1:	DO i=1,Ncol
             Node_points1: DO n=1,Nodel
                            IF(Dist(Elcor(:,n),xP(:,i),cdim) > Epsi) CYCLE
                             ! only do when Pi is node n
                             Nreg=1
                              IF (n == 3) nreg= 2
                                Subregions: DO nr=1,Nreg
                                      Mi= 4
                                  Call Gauss_Laguerre_coor(Glcor,Wi,Mi)
			Gauss_points1:	DO m=1,Mi
                               SELECT CASE (n)
                                CASE (1)
                                 xsi= 2.0*Glcor(m)-1.0
                                 dxdxb= 2.0
                                CASE (2)
                                 xsi= 1.0 -2.0*Glcor(m)
                                 dxdxb= 2.0
                                CASE (3)
                                 dxdxb= 1.0
                                 IF(nr == 1) THEN
                                   xsi= -Glcor(m)
                                 ELSE
                                   xsi= Glcor(m)
                                 END IF
                                CASE DEFAULT
                               END SELECT
                               CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)
                                 !   Shape function value
                    Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor)
                     ! Jacobian
                         Direction1:     DO j=1,2
                             IF(Isym == 0)THEN
                              iD= 2*(i-1) + j
                             ELSE
                              iD= Ndest(i,j)     !  line number in array
                             END IF
                             IF (id == 0) CYCLE           
                               nD= 2*(n-1) + j  !  column number in array
                          dUe(iD,nD)= dUe(iD,nD) + Ni(n)*C*Jac*dxdxb*Wi(m)
                         END DO Direction1
                      END DO Gauss_points1
                  END DO Subregions
                   Mi= 2
                   Call Gauss_coor(Glcor,Wi,Mi)
                     ! Assign coords/Weights
              Gauss_points2:  DO m=1,Mi
                    SELECT CASE (n)
                      CASE (1)
                       C1=-LOG(Eleng)*C
                      CASE (2)
                       C1=-LOG(Eleng)*C
                      CASE (3)
                       C1=LOG(2/Eleng)*C
                      CASE DEFAULT
                     END SELECT
                     xsi= Glcor(m)
                     CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)
                                            !       Shape function value
                   Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor)
                         ! Jacobian and normal
                         Direction2:     DO j=1,2
                           IF(Isym == 0)THEN
                             iD= 2*(i-1) + j
                           ELSE
                             iD= Ndest(i,j)       !  line number in array
                           END IF
                           IF (id == 0) CYCLE           
                           nD= 2*(n-1) + j        !  column number in array
                           dUe(iD,nD) = dUe(iD,nD) + Ni(n)*C1*Jac*Wi(m)
                         END DO Direction2
               END DO Gauss_points2
            END DO Node_points1
          END DO Colloc_points1
     RETURN
END SUBROUTINE Integ2E

SUBROUTINE Integ3(Elcor,Inci,Nodel,Ncol,xPi,N_dof,E,ny,ko,dUe,dTe,Ndest,Isym)
!--------------------------------------------------
!    Computes  [dT]e and [dU]e for 3-D problems
!    by numerical integration
!-------------------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(IN)    :: Elcor(:,:)     !   Element coordinates
INTEGER, INTENT(IN) :: Ndest(:,:)     !   Node destination vector
INTEGER, INTENT(IN) :: Inci(:)        !   Element Incidences
INTEGER, INTENT(IN) :: Nodel          !   No. of Element Nodes
INTEGER , INTENT(IN):: Ncol           !   Number of points Pi (coll. points)
REAL(iwp) , INTENT(IN)   :: xPi(:,:)       !   Array with coll. points coords.
INTEGER , INTENT(IN):: N_dof       !   Number dgrees of freedom /node (1 or 3)
INTEGER , INTENT(IN):: Isym          
REAL(iwp) , INTENT(IN) :: E,ny  !   Elastic constants (for elasticity problems)
REAL(iwp) , INTENT(IN) :: ko            
REAL(iwp) , INTENT(OUT)  :: dUe(:,:),dTe(:,:) !arrays for storing coefficients
REAL(iwp) :: Elengx,Elenge,Rmin,RLx,RLe,Glcorx(8),Wix(8),Glcore(8),Wie(8),Weit,r
REAL(iwp) :: Ni(Nodel),Vnorm(3),GCcor(3),dxr(3),Jac,Jacb,xsi,eta,xsib,etab
REAL(iwp) :: UP(N_dof,N_dof),TP(N_dof,N_dof)   !   Arrays for storing kernels
INTEGER :: i,m,n,k,ii,jj,ntr,Mi,Ki,id,nd,lnod,Ntri
INTEGER :: ldim= 2       !   Element dimension
INTEGER :: Cdim= 3       !   Cartesian dimension
ELengx= Dist((Elcor(:,3)+Elcor(:,2))/2.,(Elcor(:,4)+Elcor(:,1))/2.,Cdim)! Lxsi
ELenge= Dist((Elcor(:,2)+Elcor(:,1))/2.,(Elcor(:,3)+Elcor(:,4))/2.,Cdim)! Leta
dUe= 0.0 ; dTe= 0.0                 ! Clear arrays for summation
!---------------------------------------------------------------
!     Part 1 : Pi is not one of the element nodes
!---------------------------------------------------------------
Colloc_points:	DO i=1,Ncol
                  IF(.NOT. ALL(Inci /= i)) CYCLE
                      !  Check if incidence array contains i
               Rmin= Min_dist1(Elcor,xPi(:,i),Nodel,inci,ELengx,Elenge,ldim)
                !  Distance coll. point and element
               Mi= Ngaus(Rmin/Elengx,2)
                       !  Number of G.P. in xsi dir. for (1/r2) sing.
                Call Gauss_coor(Glcorx,Wix,Mi)
                    !  Assign coords/Weights xsi-direction
                 Ki= Ngaus(Rmin/Elenge,2)
                          !  Number of G.P. in eta dir. for (1/r2) sing.
                 Call Gauss_coor(Glcore,Wie,Ki)
                   !  Assign coords/Weights eta-direction
Gauss_points_xsi: DO m=1,Mi
                     xsi= Glcorx(m)
  Gauss_points_eta: DO k=1,Ki
                     eta= Glcore(k)
                     Weit= Wix(m)*Wie(k)
                     CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)
                      !   Shape function value
                     Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor)
                      ! Jacobian and normal
                     CALL Cartesian(GCcor,Ni,ldim,elcor)
                           ! Cart. coords of Gauss pt
                      r= Dist(GCcor,xPi(:,i),Cdim)             !  Dist. P,Q
                      dxr= (GCcor-xPi(:,i))/r             !  rx/r , ry/r
                      IF(N_dof .EQ. 1) THEN
                        UP= U(r,ko,Cdim) ; TP= T(r,dxr,Vnorm,Cdim)
                         !  Kernels Potential problem
                      ELSE
                      UP= UK(dxr,r,E,ny,Cdim) ; TP= TK(dxr,r,Vnorm,ny,Cdim)
                       !  Kernels elasticity
                      END IF
!                        Node_points:    DO n=1,Nodes
                           Direction_P:    DO ii=1,N_dof
                                            IF(Isym == 0)THEN
                                             iD= N_dof*(i-1) + ii
                                             !  line number in array
                                            ELSE
                                             iD= Ndest(i,ii)
                                              !  line number in array
                                            END IF
                                            IF (id == 0) CYCLE
                                  Direction_Q:    DO jj=1,N_dof
                                        Node_points:    DO n=1,Nodel
                                         nD= N_dof*(n-1) + jj
                                           !  column number in array
                          dUe(iD,nD)= dUe(iD,nD) + Ni(n)*UP(ii,jj)*Jac*Weit
                          dTe(iD,nD)= dTe(iD,nD) + Ni(n)*TP(ii,jj)*Jac*Weit
                                        END DO Node_points
                                 END DO Direction_Q
                              END DO Direction_P
!                         END DO Node_points
                       END DO Gauss_points_eta
              END DO Gauss_points_xsi
END DO Colloc_points
!---------------------------------------------------------
!     Part 1 : Pi is one of the element nodes
!---------------------------------------------------------
Colloc_points1: DO i=1,Ncol
                   lnod= 0
                   DO n= 1,Nodel      !   Determine which local node is Pi
                       IF(Inci(n) .EQ. i) THEN
                              lnod=n
                       END IF
                   END DO
                   IF(lnod .EQ. 0) CYCLE        !  None -> next Pi
                          Ntri= 2
                          IF(lnod > 4) Ntri=3   !  Number of triangles
                      Triangles:      DO ntr=1,Ntri
                                CALL Tri_RL(RLx,RLe,Elengx,Elenge,lnod,ntr)
                                Mi= Ngaus(RLx,2)
                             !  Number of G.P. in xsi dir. for (1/r2) sing.
!                               Mi= 8
                              Call Gauss_coor(Glcorx,Wix,Mi)
                                !  Assign coords/Weights xsi-direction
                                Ki= Ngaus(RLe,2)
                             !  Number of G.P. in eta dir. for (1/r2) sing.
!                               Ki= 8
                                Call Gauss_coor(Glcore,Wie,Ki)
                                !  Assign coords/Weights eta-direction
			Gauss_points_xsi1:	DO m=1,Mi
                                                  xsib= Glcorx(m)
				Gauss_points_eta1:	DO k=1,Ki
                                                          etab= Glcore(k)
                                                          Weit= Wix(m)*Wie(k)
                 CALL Trans_Tri(ntr,lnod,xsib,etab,xsi,eta,Jacb)
                 !  Coord transf from triang coords
                 CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)
                 !   Shape function value
                 Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor)
                 ! Jacobian and normal
                 Jac= Jac*Jacb
                 CALL Cartesian(GCcor,Ni,ldim,elcor)! Cart. coords of Gauss pt
                 r= Dist(GCcor,xPi(:,i),Cdim)   !  Dist. P,Q
                 dxr= (GCcor-xPi(:,i))/r        !  rx/r , ry/r
                 IF(N_dof .EQ. 1) THEN
                    UP= U(r,ko,Cdim) ; TP= T(r,dxr,Vnorm,Cdim)
                     !  Kernels Potential problem
                 ELSE
                    UP= UK(dxr,r,E,ny,Cdim) ; TP= TK(dxr,r,Vnorm,ny,Cdim)
                     !  Kernels elasticity
                 END IF
                  Direction_P1: DO ii=1,N_dof
                    IF(Isym == 0)THEN
                     iD= N_dof*(i-1) + ii        !  line number in array
                    ELSE
                    iD= Ndest(i,ii)             !  line number in array
                    END IF                                                                                                                           
                    IF (id == 0) CYCLE
                      Direction_Q1:  DO jj=1,N_dof
                            Node_points1: DO n=1,Nodel
                             nD= N_dof*(n-1) + jj    !  column number in array
                           dUe(iD,nD)= dUe(iD,nD) + Ni(n)*UP(ii,jj)*Jac*Weit
                           IF(Inci(n) /= i) THEN !   diagonal elements of dTe not computed
                            dTe(iD,nD)= dTe(iD,nD) + Ni(n)*TP(ii,jj)*Jac*Weit
                           END IF
                            END DO Node_points1
                      END DO Direction_Q1
                  END DO Direction_P1
                END DO Gauss_points_eta1
             END DO Gauss_points_xsi1
          END DO Triangles
     END DO Colloc_points1
   RETURN
END SUBROUTINE Integ3

SUBROUTINE Triangel_Coord(cor_tri,lnod,ntr)
!---------------------------------------------
!   Assigns local coordinates of triangular
!   subelements
!---------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: lnod,ntr      !  node, subelement no. 
REAL(iwp) , INTENT(OUT)  :: cor_tri(2,3)  !  xsi,eta of triangle nodes
REAL(iwp)                :: xsii(8),etai(8)
INTEGER                  :: nod_tri(3),n
SELECT CASE (ntr)
	CASE (1)
	SELECT CASE(lnod)
		CASE (1)
			nod_tri=(/2,3,1/)
		CASE (2)
			nod_tri=(/3,4,2/)
		CASE (3)
			nod_tri=(/1,2,3/)
		CASE (4)
			nod_tri=(/1,2,4/)
		CASE (5)
			nod_tri=(/4,1,5/)
		CASE (6)
			nod_tri=(/1,2,6/)
		CASE (7)
			nod_tri=(/4,1,7/)
		CASE (8)
			nod_tri=(/1,2,8/)
		CASE DEFAULT	
	END SELECT
	CASE (2)
	SELECT CASE(lnod)
		CASE (1)
			nod_tri=(/3,4,1/)
		CASE (2)
			nod_tri=(/4,1,2/)
		CASE (3)
			nod_tri=(/4,1,3/)
		CASE (4)
			nod_tri=(/2,3,4/)
		CASE (5)
			nod_tri=(/2,3,5/)
		CASE (6)
			nod_tri=(/3,4,6/)
		CASE (7)
			nod_tri=(/2,3,7/)
		CASE (8)
			nod_tri=(/3,4,8/)
		CASE DEFAULT
	END SELECT		
	CASE (3)
	SELECT CASE(lnod)
		CASE (5)
			nod_tri=(/3,4,5/)
		CASE (6)
			nod_tri=(/4,1,6/)
		CASE (7)
			nod_tri=(/1,2,7/)
		CASE (8)
			nod_tri=(/2,3,8/)
		CASE DEFAULT
	END SELECT		
	CASE DEFAULT
END SELECT
xsii=(/-1.0,1.0,1.0,-1.0,0.0,1.0,0.0,-1.0/)
etai=(/-1.0,-1.0,1.0,1.0,-1.0,0.0,1.0,0.0/)
cor_tri=0
DO n=1, 3
	cor_tri(1,n)= xsii(nod_tri(n))
	cor_tri(2,n)= etai(nod_tri(n))
END DO
END SUBROUTINE Triangel_Coord

SUBROUTINE Trans_Tri(ntr,lnod,xsib,etab,xsi,eta,Jacb)
!--------------------------------------------
!  Transforms from local triangle coordinates
!  to xsi,eta coordinates and computes 
!  the Jacobean of the transformation
!---------------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT (IN)       ::     xsib,etab     !  local coordinates
INTEGER, INTENT (IN)    ::    ntr,lnod      !   subelement-no, local node no.
REAL(iwp), INTENT (OUT)      ::     Jacb,xsi,eta  !   Jacobean and xsi,eta coords
REAL(iwp)                    ::      Nb(3),dNbdxb(3),dNbdeb(3),dxdxb,dxdeb,  &
                                dedxb,dedeb,cor_tri(2,3)
INTEGER                      ::      n,i
CALL Triangel_Coord(cor_tri,lnod,ntr)
!
!    Transform xsi-bar and eta-bar to xsi and eta
!
Nb(1)=0.25*(1.0+xsib)*(1.0-etab)
Nb(2)=0.25*(1.0+xsib)*(1.0+etab)
Nb(3)=0.5*(1.0-xsib)
xsi = 0.0
eta = 0.0
DO n=1,3
	xsi = xsi+Nb(n)*cor_tri(1,n) 
	eta = eta+Nb(n)*cor_tri(2,n) 
END DO
!
!    Jacobian of Transformation xsi-bar and eta-bar to xsi and eta
!
dNbdxb(1)=0.25*(1.0-etab)
dNbdxb(2)=0.25*(1.0+etab)
dNbdxb(3)=-0.5
dNbdeb(1)=-0.25*(1.0+xsib)
dNbdeb(2)=0.25*(1.0+xsib)
dNbdeb(3)=0.0
dxdxb=0.0
dedxb=0.0
dxdeb=0.0
dedeb=0.0
DO i=1,3
	dxdxb=dxdxb+dNbdxb(i)*cor_tri(1,i)
	dedxb=dedxb+dNbdxb(i)*cor_tri(2,i)
	dedeb=dedeb+dNbdeb(i)*cor_tri(2,i)
	dxdeb=dxdeb+dNbdeb(i)*cor_tri(1,i)
END DO
Jacb=dxdxb*dedeb-dedxb*dxdeb
END SUBROUTINE Trans_Tri

SUBROUTINE Tri_RL(RLx,RLe,Elengx,Elenge,lnod,ntr)
!---------------------------------------  
!  Computes ize of triangular sub-element
!  and Rmin
!---------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT (IN)   ::     Elengx,Elenge  !  length of element in xsi,eta dir
INTEGER, INTENT (IN)::      lnod,ntr       !  local node, subelement no.
REAL(iwp), INTENT (OUT)  ::     RLx,RLe        !  Lenghts of subelement
IF(lnod <= 4) THEN
	RLx=Elengx/Elenge
	RLe=Elenge/Elengx
END IF
IF (lnod == 5 .or. lnod == 7) THEN
	SELECT CASE (ntr)
		CASE (1)
			RLx=(Elengx/2.0)/Elenge
			RLe=Elenge/(Elengx/2.0)
		CASE (2)
			RLx=(Elengx/2.0)/Elenge
			RLe=Elenge/(Elengx/2.0)
		CASE (3)
			RLx=Elengx/Elenge
			RLe=Elenge/Elengx
	END SELECT
END IF
IF (lnod == 6 .or. lnod == 8) THEN
	SELECT CASE (ntr)
		CASE (1)
			RLx=Elengx/(Elenge/2)
			RLe=(Elenge/2)/Elengx
		CASE (2)
			RLx=Elengx/(Elenge/2)
			RLe=(Elenge/2)/Elengx
		CASE (3)
			RLx=Elengx/Elenge
			RLe=Elenge/Elengx
	END SELECT
END IF
END SUBROUTINE Tri_RL

  REAL FUNCTION U(r,k,Cdim)
  !-------------------------------
  !   Fundamental solution for Potential problems
  !   Temperature/Potential
  !------------------------------
  REAL(iwp),INTENT(IN)     ::  r     !   Distance between source and field point
  REAL(iwp),INTENT(IN)     ::  k     !   Conducivity
  INTEGER,INTENT(IN)  :: Cdim   !   Cartesian dimension (2-D,3-D)
  SELECT CASE (CDIM)
     CASE (2)          	!  Two-dimensional solution
        U= 1.0/(2.0*Pi*k)*LOG(1/r)
     CASE (3)          	!  Three-dimensional solution
        U= 1.0/(4.0*Pi*r*k)
     CASE DEFAULT
        U=0.0
        WRITE (12,*)'Cdim not equal 2 or 3 in Function U(...)'
  END SELECT
  END FUNCTION U

 REAL FUNCTION T(r,dxr,Vnorm,Cdim)
 !-------------------------------
 !   Fundamental solution for Potential problems
 !   Normal gradient
 !------------------------------
 REAL(iwp),INTENT(IN)::            r !   Distance between source and field point
 REAL(iwp),INTENT(IN)::       dxr(:) !   rx/r , ry/r , rz/r
 REAL(iwp),INTENT(IN)::     Vnorm(:) !   Normal vector
 INTEGER,INTENT(IN) ::    Cdim	!   Cartesian dimension
 SELECT CASE (Cdim)
    CASE (2)           	!  Two-dimensional solution
      T= -DOT_PRODUCT (Vnorm,dxr)/(2.0*Pi*r)
    CASE (3)           	!  Three-dimensional solution
      T= -DOT_PRODUCT (Vnorm,dxr)/(4.0*Pi*r*r)
    CASE DEFAULT
      T=0.0
      WRITE (12,*)'Cdim not equal 2 or 3 in Function U(...)'
 END SELECT
 END FUNCTION T
 FUNCTION dU(r,dxr,Cdim)
  !-------------------------------
  !   Derivatives of Fundamental solution for Potential problems
  !   Temperature/Potential
  !------------------------------
  REAL(iwp),INTENT(IN)::       r    !   Distance between source and field point
  REAL(iwp),INTENT(IN)::  dxr(:)    !   Distances in Cartesian directions divided by r
  REAL :: dU(UBOUND(dxr,1))    !   dU is array of same dim as dxr
  INTEGER ,INTENT(IN)            :: Cdim   !   Cartesian dimension (2-D,3-D)
  REAL(iwp) :: C
  SELECT CASE (CDIM)
     CASE (2)           !  Two-dimensional solution
      C=1/(2.0*Pi*r)
      dU(1)= C*dxr(1)
      dU(2)= C*dxr(2)
     CASE (3)           !  Three-dimensional solution
      C=1/(4.0*Pi*r**2)
      dU(1)= C*dxr(1)
      dU(2)= C*dxr(2)
      dU(3)= C*dxr(3)
     CASE DEFAULT
 END SELECT
 END FUNCTION dU
 FUNCTION dT(r,dxr,Vnorm,Cdim)
 !-------------------------------
 !   derivatives of the Fundamental solution for Potential problems
 !   Normal gradient
 !------------------------------
 INTEGER,INTENT(IN) :: Cdim     !   Cartesian dimension
 REAL(iwp),INTENT(IN)::        r     !   Distance between source and field point
 REAL(iwp),INTENT(IN)::    dxr(:)    !   Distances in Cartesian directions divided by R
 REAL(iwp),INTENT(IN)::  Vnorm(:)    !   Normal vector
 REAL(iwp) :: dT(UBOUND(dxr,1))      !   dT is array of same dim as dxr 
 REAL(iwp) :: C,COSTH
 COSTH= DOT_PRODUCT (Vnorm,dxr)
 SELECT CASE (Cdim)
    CASE (2)           !  Two-dimensional solution
     C= 1/(2.0*Pi*r**2)
     dT(1)= C*COSTH*dxr(1)
     dT(2)= C*COSTH*dxr(2)
    CASE (3)           !  Three-dimensional solution
     C= 3/(4.0*Pi*r**3)
     dT(1)= C*COSTH*dxr(1)
     dT(2)= C*COSTH*dxr(2)
     dT(3)= C*COSTH*dxr(3)
    CASE DEFAULT
 END SELECT
 END FUNCTION dT


SUBROUTINE BFLOW(Flow,xsi,eta,u,Inci,Elcor,k)
!----------------------------------------------
!    Computes flow vectors in direction tangential to the
!    Boundary 
!-----------------------------------------------
REAL(iwp) , INTENT(OUT)   :: Flow(:)  !  Flow vector
REAL(iwp) , INTENT(IN)    :: xsi,eta  !  intrinsic coordinates of point
REAL(iwp) , INTENT(IN)    :: u(:,:)     !  Nodal temperatures/potentials
INTEGER, INTENT (IN) :: Inci(:)  !  Element Incidences
REAL(iwp), INTENT (IN)    :: Elcor(:,:)  !  Element coordinates
REAL(iwp), INTENT (IN)    :: k        !  Conductivity   
REAL(iwp), ALLOCATABLE    :: Vxsi(:),Veta(:),DNi(:,:),V3(:) 
INTEGER :: Nodes,Cdim,Ldim
REAL(iwp) :: Jxsi,Jeta,Flows(2),v1(3),v2(3),CosA,CosB,CosG,CosT
Nodes= UBOUND(ELCOR,2)  !  Number of nodes
Cdim= UBOUND(ELCOR,1)   !  Cartesian Dimension
Ldim= Cdim-1            !  Local (element) dimension
ALLOCATE (Vxsi(cdim),Dni(Nodes,Ldim),v3(cdim))
IF(ldim > 1)	ALLOCATE (Veta(cdim))
!   Compute Vector(s) tangential to boundary surface
CALL Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
Vxsi(1)= Dot_Product(Dni(:,1),Elcor(1,:))
Vxsi(2)= Dot_Product(Dni(:,1),Elcor(2,:))
IF(Cdim == 2) THEN
	CALL Vector_norm(Vxsi,Jxsi)
	Flow(1)= -k*Dot_product(Dni(:,1),u(:,1))/Jxsi
ELSE
	Vxsi(3)= Dot_Product(Dni(:,1),Elcor(3,:))
	CALL Vector_norm(Vxsi,Jxsi)
	Veta(1)= Dot_Product(Dni(:,2),Elcor(1,:))
	Veta(2)= Dot_Product(Dni(:,2),Elcor(2,:))
	Veta(3)= Dot_Product(Dni(:,2),Elcor(3,:))
	CALL Vector_norm(Veta,Jeta)
	!   Flows in skew coordinate system
	Flows(1)= -k*Dot_product(Dni(:,1),u(:,1))/Jxsi
	Flows(2)= -k*Dot_product(Dni(:,2),u(:,1))/Jeta
	!   Orthoginal system
	v3= Vector_ex(Vxsi,Veta)
	Call Ortho(v3,v1,v2)
	CosA= DOT_Product(Vxsi,v1)
	CosB= DOT_Product(Veta,v2)
	CosG= DOT_Product(Veta,v1)
	CosT= DOT_Product(Vxsi,v2)
	Flow(1)= Flows(1)*CosA**2 + Flows(2)* CosG**2
	Flow(2)= Flows(1)*CosT**2 + Flows(2)* CosB**2
END IF
RETURN
END SUBROUTINE BFLOW

SUBROUTINE BStress(Stress,xsi,eta,u,t,Inci,Elcor,E,ny,IPS)
!----------------------------------------------
!    Computes stresses in a plane tangential to the
!    Boundary Element
!-----------------------------------------------
REAL(iwp) , INTENT(OUT)   :: Stress(:)!  Stress vector
REAL(iwp) , INTENT(IN)    :: xsi,eta  !  intrinsic coordinates of point
REAL(iwp) , INTENT(IN)    :: u(:,:)   !  Nodal displacements
REAL(iwp) , INTENT(IN)    :: t(:,:)   !  Nodal Tractions
INTEGER, INTENT (IN) :: Inci(:)  !  Element Incidences
REAL(iwp), INTENT (IN)    :: Elcor(:,:)  !  Element coordinates
REAL(iwp), INTENT (IN)    :: E,ny     !  Elastic constants   
INTEGER, INTENT (IN) ::  IPS ! IPS= 0 plane strain; =1 plane stress
REAL(iwp), ALLOCATABLE    :: Vxsi(:),Veta(:),DNi(:,:),Ni(:),trac(:)
INTEGER :: Nodes, Cdim, Ldim
REAL(iwp) :: Jxsi,Jeta,v1(3),v2(3),CosA, CosB, CosG, CosT,v3(3)
REAL(iwp) :: C1,C2,G,tn,ts,ts1,ts2
REAL(iwp) , ALLOCATABLE :: Dudxsi(:),Dudeta(:),Strain(:),Strains(:)
Nodes= UBOUND(ELCOR,2)  !  Number of nodes
Cdim = UBOUND(ELCOR,1)  !  Cartesian Dimension
Ldim= Cdim-1            !  Local (element) dimension
ALLOCATE (Vxsi(cdim),Veta(cdim),Dni(Nodes,Ldim),Ni(Nodes))
ALLOCATE (Dudxsi(Cdim),Dudeta(Cdim),trac(Cdim))
!   Compute Vector(s) tangential to boundary surface
CALL Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
CALL Serendip_func(Ni,xsi,eta,ldim,nodes,inci)
trac(1)= Dot_Product(Ni,t(:,1))
trac(2)= Dot_Product(Ni,t(:,2))
Vxsi(1)= Dot_Product(Dni(:,1),Elcor(1,:))
Vxsi(2)= Dot_Product(Dni(:,1),Elcor(2,:))
IF(Cdim == 2) THEN
		ALLOCATE (Strain(1))
		CALL Vector_norm(Vxsi,Jxsi)
		V3(1)=   Vxsi(2)
		V3(2)= - Vxsi(1)
		tn= Dot_Product(v3(1:2),trac)
		ts= Dot_Product(vxsi(1:2),trac)
		DuDxsi(1)= Dot_Product(Dni(:,1),u(:,1))
		DuDxsi(2)= Dot_Product(Dni(:,1),u(:,2))
		Strain(1)= Dot_Product(DuDxsi,Vxsi)/Jxsi
		IF(IPS == 2) THEN
		 Stress(1)= E*Strain(1) + ny*tn     !     plane stress
		ELSE
		 Stress(1)= 1/(1.0-ny)*(E/(1.0+ny)*Strain(1) - ny*tn)		! Plane strain
		END IF
		Stress(2)= tn
		Stress(3)= ts
ELSE
		ALLOCATE (Strain(3),Strains(3))
		trac(3)= Dot_Product(Ni,t(:,3))
		Vxsi(3)= Dot_Product(Dni(:,1),Elcor(3,:))
		CALL Vector_norm(Vxsi,Jxsi)
		Veta(1)= Dot_Product(Dni(:,2),Elcor(1,:))
		Veta(2)= Dot_Product(Dni(:,2),Elcor(2,:))
		Veta(3)= Dot_Product(Dni(:,2),Elcor(3,:))
		CALL Vector_norm(Veta,Jeta)
		V3= Vector_ex(Vxsi,veta)
		tn= Dot_Product(v3,trac)
		DuDxsi(1)= Dot_Product(Dni(:,1),u(:,1))
		DuDxsi(2)= Dot_Product(Dni(:,1),u(:,2))
		DuDxsi(3)= Dot_Product(Dni(:,1),u(:,3))
		DuDeta(1)= Dot_Product(Dni(:,2),u(:,1))
		DuDeta(2)= Dot_Product(Dni(:,2),u(:,2))
		DuDeta(3)= Dot_Product(Dni(:,2),u(:,3))
!   Strains in skew coordinate system
		Strains(1)= Dot_product(DuDxsi,Vxsi)/Jxsi
		Strains(2)= Dot_product(DuDeta,Veta)/Jeta
		Strains(3)= Dot_product(DuDeta,Vxsi)/Jeta + &
                            Dot_product(DuDxsi,Veta)/Jxsi
!   Orthogonal system
		v3= Vector_ex(Vxsi,Veta)
		CALL Ortho(v3,v1,v2)
		ts1= DOT_Product(v1,trac)
		ts2= DOT_Product(v2,trac)
		CosA= DOT_Product(Vxsi,v1)
		CosB= DOT_Product(Veta,v2)
		CosG= DOT_Product(Veta,v1)
		CosT= DOT_Product(Vxsi,v2)
!   Compute Strains
		Strain(1)= Strains(1)*CosA**2 + Strains(2)*CosG**2 + &
							 Strains(3)*CosA*CosG
		Strain(2)= Strains(1)*CosT**2 + Strains(2)*CosB**2 + &
                                                        Strains(3)*CosT*CosB
		Strain(3)= Strains(1)*CosG*CosT + Strains(2)*CosG*CosB + &
                                              Strain(3)*(CosA*CosB+CosG*CosT)
!   Compute stresses
		C1= E/(1.0-ny**2)  ;  C2= ny/(1.0-ny) ; G=E/(1.0-2*ny)
		Stress(1)= C1*(Strain(1)+ny*strain(2))+ C2*Tn
		Stress(2)= C1*(Strain(2)+ny*strain(1))+ C2*Tn
		Stress(3)= tn
		Stress(4)= G*Strain(3)
		Stress(5)= ts1
		Stress(6)= ts2
END IF
RETURN
End SUBROUTINE BStress

SUBROUTINE Stiffness_BEM(nr,xP,Nodel,N_dof,Ndofe,NodeR,   &
                         Ncode,NdofR,Ndofc,KBE,A,tc,Cdim,Elres_u,Elres_t,&
                  IncieR,LdesteR,Nbel,ListR,TypeR,Bcode,Con,E,ny,Ndest,Isym)
!---------------------------------------------
!    Computes the stiffness matrix
!    of a boundary element region
!    no symmetry
!--------------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT(INOUT):: xP(:,:)    !  Array of node coordinates
INTEGER, INTENT(IN):: nr  
INTEGER, INTENT(IN):: Ncode(:)   !  Global restraint code
INTEGER, INTENT(IN):: NdofR      
INTEGER, INTENT(IN):: Ndofc      !  No of interface degrees of freedom
INTEGER, INTENT(IN):: NodeR(:)    
INTEGER, INTENT(IN):: TypeR(:)    
INTEGER, INTENT(IN):: Cdim    
INTEGER, INTENT(IN):: IncieR(:,:) 
INTEGER, INTENT(IN):: LdesteR(:,:) 
INTEGER, INTENT(INOUT):: Bcode(:,:) 
INTEGER, INTENT(IN):: Nbel(:) 
INTEGER, INTENT(IN):: ListR(:,:) 
INTEGER, INTENT(IN):: Nodel 
INTEGER, INTENT(IN):: N_dof
INTEGER, INTENT(IN):: Ndofe 
INTEGER, INTENT(IN):: Isym 
INTEGER, INTENT(IN):: Ndest(:,:) 
REAL(iwp), INTENT(INOUT):: Elres_u(:,:),Elres_t(:,:)  
REAL(iwp), INTENT(INOUT)      :: E,ny,Con    
REAL(iwp), INTENT(OUT)  :: KBE(:,:) !  Stiffness matrix
REAL(iwp), INTENT(OUT)  :: A(:,:)   ! u due to únit values ui
REAL(iwp), INTENT(OUT)  :: tc(:)    ! interface tractions
!   temporal arrays :
REAL(iwp), ALLOCATABLE :: dUe(:,:),dTe(:,:),Diag(:,:)
REAL(iwp), ALLOCATABLE :: Lhs(:,:)
REAL(iwp), ALLOCATABLE :: Rhs(:),RhsM(:,:) ! right hand sides
REAL(iwp), ALLOCATABLE :: u1(:),u2(:,:)    ! results
REAL(iwp), ALLOCATABLE :: Elcor(:,:) 
REAL(iwp)              :: Scat,Scad
INTEGER					:: NdofF
INTEGER :: Dof,k,l,nel
INTEGER :: n,m,Pos,i,j,nd,ne
ALLOCATE(dTe(NdofR,Ndofe),dUe(NdofR,Ndofe))   
ALLOCATE(Diag(NdofR,N_dof))                   
ALLOCATE(Lhs(NdofR,NdofR),Rhs(NdofR),RhsM(NdofR,NdofR))
ALLOCATE(u1(NdofR),u2(NdofR,NdofR)) 
ALLOCATE(Elcor(Cdim,Nodel))        
!------------------------------------------
!     Scaling
!------------------------------------------
CALL Scal(E,xP,Elres_u,Elres_t,Cdim,Scad,Scat)
!----------------------------------------------------------------
!  Compute and assemble element coefficient matrices
!----------------------------------------------------------------
Lhs= 0.0
Diag= 0.0
Rhs= 0.0
RhsM= 0.0
Elements_1:&
DO Nel=1,Nbel(nr)
		ne= ListR(nr,Nel)
		Elcor(:,:)= xP(:,IncieR(ne,:))      !    gather element coords
		IF(Cdim == 2) THEN
                        IF(N_dof == 1) THEN
   CALL Integ2P(Elcor,IncieR(ne,:),Nodel,NodeR(nr),xP,Con,dUe,dTe,Ndest,Isym)
			ELSE
   CALL Integ2E(Elcor,IncieR(ne,:),Nodel,NodeR(nr),xP,E,ny,dUe,dTe,Ndest,Isym)  
			END IF
		ELSE
   CALL Integ3(Elcor,IncieR(ne,:),Nodel,NodeR(nr),xP,N_dof,E,ny,    &
                Con,dUe,dTe,Ndest,Isym)
		END IF
                CALL AssemblyMR(ne,N_dof,Ndofe,Nodel,Lhs,Rhs,RhsM,  &
                DTe,DUe,LdesteR(ne,:),Ncode,Bcode,Diag,Elres_u,Elres_t,Scad)
END DO &
Elements_1
!------------------------------------------------------------
!  Add azimuthal integral for infinite regions
!------------------------------------------------------------
IF(TypeR(nr) == 2) THEN
	DO m=1, NodeR(nr)
                DO n=1, N_dof
                        k=N_dof*(m-1)+n
			Diag(k,n) = Diag(k,n) + 1.0
		END DO
	END DO		 
END IF
!-------------------------------------------------------------
!  Add Diagonal coefficients
!-------------------------------------------------------------
Nodes_global: &
DO m=1,NodeR(nr)
	Degrees_of_Freedoms_node: &
        DO n=1, N_dof
                DoF = (m-1)*N_dof + n   !  global degree of freedom no.
                k = (m-1)*N_dof + 1     !  address in coeff. matrix (row)
                l = k + N_dof - 1       !  address in coeff. matrix (column)
                IF (NCode(DoF) == 1 .or. NCode(DoF) == 2) THEN
                        ! Dirichlet - Add Diagonal to Rhs
			Pos = 0
			Nel = 0
             !   get local degree of freedom no corresponding to global one
			Elements_all: &
			DO i=1,Nbel(nr)	
				ne= ListR(nr,i)
				Degrees_of_freedom_elem: &
				DO j=1,Ndofe
					IF (DoF == LdesteR(ne,j)) THEN
						Nel = ne
						Pos = j
						EXIT				 
					END IF
				END DO &
				Degrees_of_freedom_elem
			IF (Nel /= 0) EXIT 
			END DO &
			Elements_all
			Rhs(k:l) = Rhs(k:l) - Diag(k:l,n)*Elres_u (Nel,Pos)
			IF(NCode(DoF) == 2)THEN
                          RhsM(k:l,DoF) = RhsM(k:l,DoF) - Diag(k:l,n) / Scad       
			END IF
		ELSE
			Lhs(k:l,Dof)= Lhs(k:l,Dof) + Diag(k:l,n)
                        	! Neuman - Add to Lhs
		END IF	
	END DO &
	Degrees_of_Freedoms_node
END DO	&
Nodes_global	
!    Solve problem 
CALL Solve_Multi(Lhs,Rhs,RhsM,u1,u2)
!------------------------------------------
!		Back - Scaling
!------------------------------------------
DO N=1,NdofC
	u1(N)= u1(N) / Scat
	u2(N,:)= u2(N,:) / Scat
END DO
M=NdofC
NdofF= NdofR-NdofC
DO N=1,NdofF
	M=M+1
	IF(NCode(M) == 0) THEN
		u1(M)= u1(M) * Scad
		u2(M,:)= u2(M,:) * Scad
	ELSE
		u1(M)= u1(M) / Scat
		u2(M,:)= u2(M,:) / Scat
	END IF
END DO
	Elres_u(:,:)= Elres_u(:,:) * Scad
	Elres_t(:,:)= Elres_t(:,:) / Scat

!--------------------------------------
!  Gather element results due to 
!  zero Dirichlet conditions at the interface
!--------------------------------------
Elements2:	&
DO nel=1,Nbel(nr)
	ne= ListR(nr,nel)
	D_o_F1:		&
	DO nd=1,Ndofe
		IF(Ncode(LdesteR(ne,nd)) == 0) THEN
			Elres_u(ne,nd) =  u1(LdesteR(ne,nd))
		ELSE IF(Bcode(ne,nd) == 1 .or. Bcode(ne,nd) == 2) THEN
			Elres_t(ne,nd) =  u1(LdesteR(ne,nd))
		END IF
	END DO &
	D_o_F1
END DO &
Elements2

!------------------------------------
!   Gather stiffness matrix KBE and matrix A
!------------------------------------
Interface_DoFs: &
DO N=1,Ndofc
	KBE(N,:)= u2(N,:)
	tc(N)= u1(N)
END DO &
Interface_DoFs
A= 0.0
M=NdofC
Free_DoFs: &
DO N=1,NdofF
	M= M+1
	A(N,1:NdofC)= u2(M,:)
END DO &
Free_DoFs
DEALLOCATE (dUe,dTe,Diag,Lhs,Rhs,RhsM,u1,u2,Elcor)
RETURN
END SUBROUTINE Stiffness_BEM

SUBROUTINE Solve_Multi(Lhs,Rhs,RhsM,u,uM)
!---------------------------------------------
!    Solution of system of equations
!    by Gauss Elimination
!    for multple right hand sides
!---------------------------------------------
REAL(iwp) ::    Lhs(:,:)    !    Equation Left hand side
REAL(iwp) ::    Rhs(:)      !    Equation right hand side 1
REAL(iwp) ::    RhsM(:,:)   !    Equation right hand sides 2
REAL(iwp) ::    u(:)        !    Unknowns 1
REAL(iwp) ::    uM(:,:)     !    Unknowns 2
REAL(iwp) ::    FAC
INTEGER  M,Nrhs            !    Size of system
INTEGER  i,n,nr
M= UBOUND(RhsM,1) ; Nrhs= UBOUND(RhsM,2)
!  Reduction
Equation_n: &
DO n=1,M-1
   IF(ABS(Lhs(n,n)) < 1.0E-10) THEN
     CALL Error_Message('Singular Matrix')
   END IF
   Equation_i: &
	 DO i=n+1,M
     FAC= Lhs(i,n)/Lhs(n,n)
     Lhs(i,n+1:M)= Lhs(i,n+1:M) - Lhs(n,n+1:M)*FAC
     Rhs(i)= Rhs(i) - Rhs(n)*FAC
		 RhsM(i,:)= RhsM(i,:) - RhsM(n,:)*FAC
   END DO  & 
	 Equation_i
END DO &
Equation_n
!     Backsubstitution 
Unknown_1: &
DO n= M,1,-1	  
	 u(n)= -1.0/Lhs(n,n)*(SUM(Lhs(n , n+1:M)*u(n+1:M)) - Rhs(n))
END DO &
Unknown_1
Load_case: &
DO Nr=1,Nrhs
  Unknown_2: &
	DO n= M,1,-1	  
         uM(n,nr)= -1.0/Lhs(n,n)*(SUM(Lhs(n,n+1:M)*uM(n+1:M , nr)) - RhsM(n,nr))
  END DO &
	Unknown_2
END DO &
Load_case
RETURN
END SUBROUTINE Solve_Multi

SUBROUTINE AssemblyMR(Nel,N_dof,Ndofe,Nodel,Lhs,Rhs,RhsM,      &
                      DTe,DUe,Ldest,Ncode,Bcode,Diag,Elres_u,Elres_t,Scad)
!---------------------------------------------
!  Assembles Element contributions DTe , DUe
!  into global matrix Lhs, vector Rhs
!  and matrix RhsM 
!---------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN)      :: NEL
REAL(iwp)            :: Lhs(:,:) !  Eq.left hand side
REAL(iwp)            :: Rhs(:)   !  Right hand side
REAL(iwp)            :: RhsM(:,:) ! Matrix of right hand sides
REAL(iwp), INTENT(IN):: DTe(:,:),DUe(:,:)   !  Element arrays
REAL(iwp), INTENT(INOUT)                     :: Elres_u(:,:),Elres_t(:,:)  
INTEGER , INTENT(IN)    :: LDest(:) ! Element destination vector
INTEGER , INTENT(IN) :: NCode(:) ! Boundary code (global) 
INTEGER , INTENT(IN) :: BCode(:,:) ! Boundary code (global) 
INTEGER , INTENT(IN) :: N_dof
INTEGER , INTENT(IN) :: Ndofe
INTEGER , INTENT(IN) :: Nodel
REAL(iwp) :: Diag(:,:) ! Array containing diagonal coeff of DT
INTEGER :: n,Ncol,m,k,l
REAL(iwp) :: Scad
DoF_per_Element:&
DO m=1,Ndofe  
	Ncol=Ldest(m)      !   Column number 
	IF(BCode(nel,m) == 0) THEN    !   Neumann BC
		Rhs(:) = Rhs(:) + DUe(:,m)*Elres_t(nel,m)
!     The assembly of dTe depends on the global BC
		IF (NCode(Ldest(m)) == 0) THEN	
			Lhs(:,Ncol)=  Lhs(:,Ncol) + DTe(:,m)
		ELSE
			Rhs(:) = Rhs(:) - DTe(:,m) * Elres_u(nel,m)
		END IF
	ELSE IF(BCode(nel,m) == 1) THEN   !   Dirichlet BC
		Lhs(:,Ncol) = Lhs(:,Ncol) - DUe(:,m)
		Rhs(:)= Rhs(:) - DTe(:,m) * Elres_u(nel,m)
	END IF
	IF(BCode(nel,m) == 2) THEN   !   Interface
		Lhs(:,Ncol) = Lhs(:,Ncol) - DUe(:,m)
	END IF
	IF(NCode(Ldest(m)) == 2) THEN   !   Interface
		RhsM(:,Ncol)= RhsM(:,Ncol) - DTe(:,m) / Scad
	END IF
END DO &
DoF_per_Element
!				Sum of off-diagonal coefficients
DO n=1,Nodel
        DO k=1,N_dof
                l=(n-1)*N_dof+k
		Diag(:,k)= Diag(:,k) - DTe(:,l)
	END DO
END DO
	RETURN
END SUBROUTINE AssemblyMR

	SUBROUTINE Error_Message(TEXT)
  !---------------------------------------
  ! Writes an error message onto an error file
  ! and the console and terminates the program
  !--------------------------------------
   implicit none
   CHARACTER (LEN=*) TEXT
   LOGICAL :: EXST
   INQUIRE(FILE='ERR.DAT', EXIST= EXST)
   IF(EXST) THEN
    OPEN (UNIT=99,FILE='ERR.DAT',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND')
   ELSE
    OPEN (UNIT=99,FILE='ERR.DAT',STATUS='NEW',FORM='FORMATTED')
   END IF
    WRITE(99,'(A)') TEXT
!   CALL PERROR('Fatal Error, see file ERR.DAT')
   CALL EXIT(1)
  END SUBROUTINE Error_Message

  SUBROUTINE Solve(Lhs,Rhs,F)
  !---------------------------------------------
  !    Solution of system of equations
  !    by Gauss Elimination
  !---------------------------------------------
  IMPLICIT NONE
  REAL(iwp)  ::    Lhs(:,:)    !    Equation Left hand side
  REAL(iwp)  ::    Rhs(:)      !    Equation right hand side
  REAL(iwp)  ::    F(:)        !    Unknowns
  REAL(iwp)  ::    FAC
  INTEGER           ::              M           !    Size of system
  INTEGER           ::              i,n
  M= UBOUND(Lhs,1)
  !  Reduction
  Equation_n:DO n=1,M-1
                IF(Lhs(n,n) < 1.0E-14 .and. Lhs(n,n) > -1.0E-14) THEN
     CALL Error_Message('Singular Matrix')
                END IF
    Equation_i: DO i=n+1,M
                 FAC= Lhs(i,n)/Lhs(n,n)
                 Lhs(i,n+1:M)= Lhs(i,n+1:M) - Lhs(n,n+1:M)*FAC
                 Rhs(i)= Rhs(i) - Rhs(n)*FAC
                END DO   Equation_i
             END DO Equation_n
  !     Backsubstitution
  Unknown_n: DO n= M,1,-1
                F(n)= -1.0/Lhs(n,n)*(SUM(Lhs(n,n+1:M)*F(n+1:M)) - Rhs(n))
        END DO Unknown_n
  RETURN
  END SUBROUTINE Solve

 SUBROUTINE Assembly(Lhs,Rhs,DTe,DUe,Ldest,BCode,Ncode &
           ,Elres_u,Elres_te,Diag,Ndofe,N_dof,Nodel,Fac)
!---------------------------------------------
!  Assembles Element contributions DTe , DUe
!  into global matrix Lhs and vector Rhs
!  Also sums off-diagonal coefficients 
!  for the computation of diagonal coefficients
!---------------------------------------------
IMPLICIT NONE
REAL(iwp)            :: Lhs(:,:),Rhs(:)    ! Global arrays
REAL(iwp), INTENT(IN):: DTe(:,:),DUe(:,:)  ! Element arrays
INTEGER , INTENT(IN)    :: LDest(:)          ! Element destination vector
INTEGER , INTENT(IN)            :: BCode(:)  ! Boundary code(local)
INTEGER , INTENT(IN)            :: NCode(:)  ! Boundary code (global) 
INTEGER , INTENT(IN)            :: Ndofe     ! D.o.F´s / Elem
INTEGER , INTENT(IN)            :: N_dof      ! D.o.F´s / Node
INTEGER , INTENT(IN)            :: Nodel     ! Nodes/Element
REAL(iwp) , INTENT(IN)               :: Elres_u(:)       ! vector u for element
REAL(iwp) , INTENT(IN)               :: Elres_te(:)      ! vector t for element
REAL(iwp) , INTENT(IN)               :: Fac(:)           ! Mult. factors for symmetry  
REAL(iwp) :: Diag(:,:)          ! Array containing diagonal coeff of DT
INTEGER   :: n,Ncol,m,k,l
DoF_per_Element:&
DO m=1,Ndofe  
	Ncol=Ldest(m)      !   Column number 
	IF(BCode(m) == 0) THEN    !   Neumann BC
		Rhs(:) = Rhs(:) + DUe(:,m)*Elres_te(m)*Fac(m)
	!     The assembly of dTe depends on the global BC
		IF (NCode(Ldest(m)) == 0 .and. Ncol /= 0) THEN	
			Lhs(:,Ncol)=  Lhs(:,Ncol) + DTe(:,m)*Fac(m)
		ELSE
			Rhs(:) = Rhs(:) - DTe(:,m) * Elres_u(m)*Fac(m)
		END IF
	END IF
	IF(BCode(m) == 1) THEN   !   Dirichlet BC
		Lhs(:,Ncol) = Lhs(:,Ncol) - DUe(:,m)*Fac(m)
		Rhs(:)= Rhs(:) - DTe(:,m) * Elres_u(m)*Fac(m)
	END IF
END DO &
DoF_per_Element
!				Sum of off-diagonal coefficients
DO n=1,Nodel
        DO k=1,N_dof
                l=(n-1)*N_dof+k
		Diag(:,k)= Diag(:,k) - DTe(:,l)
	END DO
END DO
	RETURN
END SUBROUTINE Assembly

SUBROUTINE Mirror(Isym,nsy,Nodes,Elcor,Fac,Incie,Ldeste,Elres_te,Elres_ue &
                 ,Nodel,N_dof,Cdim)
!--------------------------------------------
!     Creates mirror image of element
!--------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN)			::	Isym       ! symmetry indicator 
INTEGER, INTENT(IN)			::  nsy        ! symmetry count
INTEGER, INTENT(IN)			::  nodes      ! highest node no
REAL(iwp), INTENT(IN OUT)            ::  Elcor(:,:) ! Coords (will be modified)
REAL(iwp), INTENT(IN OUT)            ::      Elres_te(:)! Tractions of Element 
REAL(iwp), INTENT(IN OUT)            ::      Elres_ue(:)! Displacements of Element 
REAL(iwp), INTENT(OUT)                               ::  Fac(:)     ! Multiplication factors 
INTEGER, INTENT(IN OUT)	::  Incie(:)	 ! Incidences     (will be
INTEGER, INTENT(IN OUT)	::	Ldeste(:)	 ! Destinations    modified)
INTEGER, INTENT(IN)			::  Nodel      ! Nodes per element
INTEGER, INTENT(IN)                     ::  N_dof       ! d.o.F. per Node 
INTEGER, INTENT(IN)			::  Cdim       ! Cartesian dimension
REAL(iwp)  :: TD(3) ! Transformation vector (diagonal elements of T)
INTEGER :: n,m,Ison1,Ison2,Ison3,i
IF(nsy == 1)RETURN
!     Assign coefficients of TD 
SELECT CASE (nsy-1)
	CASE(1)
		TD=(/-1.0,1.0,1.0/)
	CASE(2)
		TD=(/-1.0,-1.0,1.0/) 
	CASE(3)
		TD=(/1.0,-1.0,1.0/) 
	CASE(4)
		TD=(/1.0,1.0,-1.0/) 
	CASE(5)
		TD=(/-1.0,1.0,-1.0/) 
	CASE(6)
		TD=(/-1.0,-1.0,-1.0/) 
	CASE(7)
		TD=(/1.0,-1.0,-1.0/) 
END SELECT
!     generate coordinates and incidences
Nodes1: &
DO n=1,nodel
	Direction: &
	DO m=1,Cdim
		Elcor(m,n)= Elcor(m,n)*TD(m)
	END DO & 
	Direction 
	!	 Check if point is on any symmetry plane
	Ison1= 0 
	Ison2= 0 
	Ison3= 0
	IF(Elcor(1,n)==0.0) Ison1=1  
	IF(Elcor(2,n)==0.0) Ison2=1  
	IF(Cdim > 2 .AND. Elcor(3,n)==0.0) Ison3=1  
	!   only change incidences for unprimed nodes
	IF(ison1==1 .AND. nsy-1==1)CYCLE
	IF(ison2==1 .AND. nsy-1==3) CYCLE
	IF(ison1+ison2+ison3 > 1 .AND. nsy-1<4) CYCLE
	Incie(n)= Incie(n) + Nodes
END DO &
Nodes1
!     generate multiplication factors elast. Problems only
IF(N_dof > 1) THEN
	I=0
	Nodes2: &
	DO n=1,nodel
		Degrees_of_freedom1: &
                DO m=1,N_dof
			I=I+1
                        Fac(I)= TD(m)  ! Multiplication factor for symmetry
		END DO & 
		Degrees_of_freedom1
	END DO &
	Nodes2
END IF
!   Reverse destination vector for selected elem
SELECT CASE (nsy-1)
CASE (1,3,4,6)
        CALL Reverse(Incie,elcor,ldeste,Elres_te,Elres_ue,N_dof,Cdim,nodel)
CASE DEFAULT
END SELECT
RETURN
END SUBROUTINE Mirror

SUBROUTINE Reverse(Inci,elcor,ldest,Elres_te,Elres_ue,N_dof,Cdim,nodel)
!--------------------------------------
!  reverses incidences, destination vector 
!  and co-ordinates
!  so that outward normal is reversed
!--------------------------------------
IMPLICIT NONE
INTEGER, INTENT (INOUT) :: Inci(:)			!   Incidences
REAL(iwp), INTENT (INOUT)    :: Elcor(:,:)     !   Coordinates
REAL(iwp), INTENT (INOUT)    :: Elres_te(:)  !   Tractions of Element
REAL(iwp), INTENT (INOUT)    :: Elres_ue(:)  !   Displacements of Element 
INTEGER, INTENT (INOUT) :: Ldest(:)	    !   Destination vector
INTEGER, INTENT (IN)            :: N_dof  !   No of degrees of freedom per node
INTEGER, INTENT (IN)            :: Cdim  !   Cartesian dimension
INTEGER, INTENT (IN)            :: Nodel !   No of nodes per element
REAL(iwp), ALLOCATABLE               :: Elcort(:,:)             !   Temps
REAL(iwp), ALLOCATABLE               :: Elres_tet(:)            !   Temps
REAL(iwp), ALLOCATABLE               :: Elres_uet(:)            !   Temps
INTEGER, ALLOCATABLE		:: Incit(:),Ldestt(:)   !   Temps
INTEGER                         :: Node(8)          !    reversing sequence
INTEGER                         :: n,nc,Nchanges
ALLOCATE (Incit(nodel),Elcort(Cdim,nodel),     &
           Ldestt(Nodel*n_dof),Elres_tet(Nodel*n_dof),Elres_uet(Nodel*n_dof))
Incit= Inci
Elcort= Elcor
Ldestt= Ldest
Elres_tet= Elres_te
Elres_uet= Elres_ue
SELECT CASE (Cdim)
	CASE (2)    !   2-D problem
		Node(1:2)= (/2,1/) ;  Nchanges= 2
	CASE (3)    !    3-D  problem
		Node= (/1,4,3,2,8,7,6,5/) ; Nchanges= nodel !-1
END SELECT
Number_changes: &
DO n=1,Nchanges
	nc= Node(n)
	inci(n)= Incit(nc) ; Elcor(:,n)= Elcort(:,nc)
        Ldest(N_dof*(n-1)+1:N_dof*n)= Ldestt(N_dof*(nc-1)+1:N_dof*nc)
        Elres_te(N_dof*(n-1)+1:N_dof*n)=Elres_tet(N_dof*(nc-1)+1:N_dof*nc) 
        Elres_ue(N_dof*(n-1)+1:N_dof*n)=Elres_uet(N_dof*(nc-1)+1:N_dof*nc) 
END DO &
Number_changes
DEALLOCATE(Incit,Elcort,Ldestt,Elres_tet,Elres_uet)
RETURN
END SUBROUTINE Reverse

SUBROUTINE Jobin(Title,Cdim,N_dof,Toa,Nreg,Ltyp,Con,E,ny &
                ,Isym,nodel,nodes,Maxe)
!------------------------------------------------
!    Subroutine to read in basic job information
!------------------------------------------------
CHARACTER(LEN=80), INTENT(OUT):: Title
INTEGER, INTENT(OUT) :: Cdim,N_dof,Toa,Nreg,Ltyp,Isym,nodel
INTEGER, INTENT(OUT) :: Nodes,Maxe
REAL(iwp), INTENT(OUT)    :: Con,E,ny
READ(11,'(A80)') Title
WRITE(12,*)'Project:',Title
READ(11,*) Cdim
WRITE(12,*)'Cartesian_dimension:',Cdim
READ(11,*) N_dof        !    Degrees of freedom per node
IF(N_Dof == 1) THEN
        WRITE(12,*)'Potential Problem'
ELSE
        WRITE(12,*)'Elasticity Problem'
END IF
IF(N_dof == 2)THEN 
        READ(11,*) Toa   !Toa ....Type of analysis
                        ! (solid plane strain = 1,solid plane stress = 2)
	IF(Toa == 1)THEN
                WRITE(12,*)'Type of Analysis: Solid Plane Strain'
	ELSE
                WRITE(12,*)'Type of Analysis: Solid Plane Stress'
	END IF
END IF	
READ(11,*) Nreg       !   Type of region
IF(NReg == 1) THEN
        WRITE(12,*)'Finite Region'
ELSE
        WRITE(12,*)'Infinite Region'
END IF
READ(11,*) Isym       !   Symmetry code
SELECT CASE (isym)
CASE(0)
WRITE(12,*)'No symmetry'
CASE(1) 
WRITE(12,*)'Symmetry about y-z plane' 
CASE(2)
WRITE(12,*)'Symmetry about y-z and x-z planes'
CASE(3) 
WRITE(12,*)'Symmetry about all planes'
END SELECT
READ(11,*) Ltyp        !   Element type
IF(Ltyp == 1) THEN
WRITE(12,*)'Linear Elements'
ELSE
WRITE(12,*)'Quadratic Elements'
END IF
!     Determine number of nodes per element
IF(Cdim == 2) THEN    !    Line elements
 IF(Ltyp == 1) THEN
  Nodel= 2
 ELSE
  Nodel= 3
 END IF
ELSE                  !    Surface elements
 IF(Ltyp == 1) THEN
  Nodel= 4
 ELSE
  Nodel= 8
 END IF
END IF
!   Read properties
IF(N_dof == 1) THEN
        READ(11,*) Con
        WRITE(12,*)'Conductivity=',Con
ELSE	
        READ(11,*) E,ny
	IF(ToA == 2)THEN		!	Solid Plane Stress
		E=E*((1+2*ny)/(1+ny)**2)
		ny =ny/(1+ny)
	END IF			
        WRITE(12,*)'Modulus:',E
        WRITE(12,*)'Poissons ratio:',ny
END IF

READ(11,*) Nodes
WRITE(12,*)'Number of Nodes of System:',Nodes
READ(11,*) Maxe     
WRITE(12,*)'Number of Elements of System:', Maxe
RETURN
END SUBROUTINE Jobin

SUBROUTINE JobinMR(Title,Cdim,N_dof,Toa,Ltyp,Isym,nodel,nodes,maxe)
!------------------------------------------------
!    Subroutine to read in basic job information
!------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(OUT)   :: Cdim      !   Cartesian dimension
INTEGER                :: ldim      !   Dimension of Element
INTEGER, INTENT(OUT)   :: N_dof      !   No. of degeres of freedom per node
INTEGER, INTENT(OUT)   :: Toa       !   Type of analysis
                                    ! (plane strain = 1, plane stress = 2)
INTEGER                :: Ltyp      !   Element type(linear = 1, quadratic = 2)
INTEGER, INTENT(OUT)   :: Nodel     !   No. of nodes per element
INTEGER, INTENT(OUT)   :: Nodes     !   No. of nodes of System
INTEGER, INTENT(OUT)   :: Maxe      !   Number of Elements of System
INTEGER                :: Isym      !   Symmetry code
CHARACTER(LEN=80)      :: Title     !               Title of calculation
INTEGER                :: nr,nb

READ(11,'(A80)') Title
WRITE(12,*)'Project:',Title
READ(11,*) Cdim
WRITE(12,*)'Cartesian_dimension:',Cdim
ldim= Cdim - 1
READ(11,*) N_dof        !    Degrees of freedom per node
IF(N_Dof == 1) THEN
        WRITE(12,*)'Potential Problem'
ELSE
        WRITE(12,*)'Elasticity Problem'
END IF
IF(N_dof == 2)THEN 
        READ(11,*) Toa     !  Toa ....Type of analysis
                          !    (solid plane strain = 1,solid plane stress = 2)
	IF(Toa == 1)THEN
                WRITE(12,*)'Type of Analysis: Solid Plane Strain'
	ELSE
                WRITE(12,*)'Type of Analysis: Solid Plane Stress'
	END IF
END IF	
READ(11,*) Ltyp        !   Element type
IF(Ltyp == 1) THEN
WRITE(12,*)'Linear Elements'
ELSE
WRITE(12,*)'Quadratic Elements'
END IF
!     Determine number of nodes per element
IF(Cdim == 2) THEN    !    Line elements
 IF(Ltyp == 1) THEN
	Nodel= 2
 ELSE
	Nodel= 3
 END IF
ELSE                  !    Surface elements
 IF(Ltyp == 1) THEN
	Nodel= 4
 ELSE
	Nodel= 8
 END IF
END IF
READ(11,*) Nodes
WRITE(12,*)'Number of Nodes of System:',Nodes
READ(11,*) Maxe     
WRITE(12,*)'Number of Elements of System:', Maxe
END SUBROUTINE JobinMR

SUBROUTINE Reg_Info(Nregs,ToA,N_dof,TypeR,ConR,ER,nyR,Nbel,ListR)
!----------------------------------------------------------------
!    Subroutine to read in basic job information for each region
!----------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN) ::      Nregs        ! Number of Regions 
INTEGER,INTENT(IN) ::      ToA          ! Type of analysis
                       ! (solid plane strain = 1,solid plane stress = 2)
INTEGER,INTENT(INOUT)    ::      TypeR(:)
                     ! Type of BE-regions (1 == finite, 2 == Infinite)
INTEGER,INTENT(INOUT)  ::      Nbel(:)
                     ! Number of Boundary Elements each region
INTEGER,INTENT(INOUT)  ::      ListR(:,:)
                     ! List of Elementnumbers each region
INTEGER,INTENT(IN) ::      N_dof       ! No. of degeres of freedom per node
REAL(iwp),INTENT(INOUT) ::      ConR(:)    ! Conductivity of each region
REAL(iwp),INTENT(INOUT) ::      ER(:)      ! Youngsmodulus of regions
REAL(iwp),INTENT(INOUT) ::      nyR(:)     ! Poissons ratio of regions
INTEGER            ::      Isym       ! Symmetrycode
INTEGER            ::      nr,nb
ListR= 0
Region_loop: &
DO nr=1,Nregs
        WRITE(12,*)' Region ',nr
        READ(11,*) TypeR(nr)       !   Type of region
	IF(TypeR(nr) == 1) THEN
                WRITE(12,*)'Finite Region'
	ELSE
                WRITE(12,*)'Infinite Region'
	END IF
        READ(11,*) Isym       !   Symmetry code
	SELECT CASE (Isym)
	CASE(0)
        WRITE(12,*)'No symmetry'
	CASE(1) 
        WRITE(12,*)'Symmetry about y-z plane' 
	CASE(2)
        WRITE(12,*)'Symmetry about y-z and x-z planes'
	CASE(3) 
        WRITE(12,*)'Symmetry about all planes'
	END SELECT
	!   Read properties
        IF(N_dof == 1) THEN
                READ(11,*) ConR(nr)
                WRITE(12,*)'Conductivity=',ConR(nr)
	ELSE	
                READ(11,*) ER(nr),nyR(nr)
                IF(ToA == 2)THEN                        !   Solid Plane Stress
			ER(nr)=ER(nr)*((1+2*nyR(nr))/(1+nyR(nr))**2)
			nyR(nr) =nyR(nr)/(1+nyR(nr))
		END IF			
                WRITE(12,*)'Youngsmodulus:',ER(nr)
                WRITE(12,*)'Poissons ratio:',nyR(nr)
	END IF
        READ(11,*)Nbel(nr)
!	IF(Nbel(nr) > MaxeR)MaxeR= Nbel(nr)
        READ(11,*)(ListR(nr,nb),nb=1,Nbel(nr))
        WRITE(12,*) ' List of Boundary Elements: '
        WRITE(12,*)(ListR(nr,nb),nb=1,Nbel(nr))
END DO &
Region_loop
RETURN
END SUBROUTINE Reg_info

SUBROUTINE BCInput(Elres_u,Elres_t,Bcode,nodel,ndofe,n_dof) 
!-----------------------------------------------------------
!			Reads boundary conditions
!-----------------------------------------------------------
IMPLICIT NONE
REAL(iwp),INTENT(INOUT)    :: Elres_u(:,:)  !  Element results , u
REAL(iwp),INTENT(INOUT)    :: Elres_t(:,:)  !  Element results , t 
INTEGER,INTENT(INOUT) :: BCode(:,:)    !  Element BC´s
INTEGER,INTENT(IN)  :: nodel         !  Nodes per element
INTEGER,INTENT(IN)  :: ndofe         !  D.o.F. per Element
INTEGER,INTENT(IN)  :: n_dof          !  D.o.F per Node
INTEGER :: NE_u,NE_t,n,Nel,m,na  
WRITE(12,*)''
WRITE(12,*)'Elements with Dirichlet BC´s: '
WRITE(12,*)''
Elres_u = 0.0_iwp  ! Default prescribed values for u = 0.0
BCode = 0             ! Default BC= Neumann Condition			
READ(11,*)NE_u    
IF(NE_u > 0) THEN
DO n=1,NE_u
        READ(11,*) Nel,(Elres_u(Nel,m),m=1,Ndofe)                 
!       READ(11,*) Nel,(BCode(Nel,m),m=1,Ndofe)
END DO  
        BCode(Nel,:)=1
        WRITE(12,*)'Element ',Nel,'  Prescribed values: '
        CALL FLUSH(12)
        Na= 1
        DO M= 1,Nodel , Nodel - 1
                WRITE(12,*) Elres_u(Nel,na:na+n_dof-1)
                Na= na+N_dof
        END DO
END IF
WRITE(12,*)''
WRITE(12,*)'Elements with Neuman BC´s: '
WRITE(12,*)''
CALL FLUSH(12)
Elres_t(:,:) = 0.0_iwp   !   Default prescribed values = 0.0
READ(11,*)NE_t            
DO n=1,NE_t
        READ(11,*) Nel,(Elres_t(Nel,m),m=1,Ndofe)  
END DO   
        WRITE(12,*)'Element ',Nel,'  Prescribed values: '
	Na= 1
        DO M= 1,Nodel, Nodel - 1
                WRITE(12,*) Elres_t(Nel,na:na+n_dof-1)
                Na= na+N_dof
        END DO 
CALL FLUSH(12)
RETURN
END SUBROUTINE BCInput

SUBROUTINE Geomin(Nodes,Maxe,xp,Inci,Nodel,Cdim)
!------------------------------------
!   Inputs mesh geometry 
!-------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) ::	Nodes			!   Number of nodes
INTEGER, INTENT(IN) ::	Maxe			!   Number of elements
INTEGER, INTENT(IN) ::	Nodel			!   Number of Nodes of elements
INTEGER, INTENT(IN) ::	Cdim			!   Cartesian Dimension
REAL(iwp), INTENT(INOUT)   ::  xP(:,:)         !   Node co-ordinates
REAL(iwp)   ::      xmax(Cdim),xmin(Cdim),delta_x(Cdim)
INTEGER, INTENT(INOUT):: Inci(:,:) !   Element incidences
INTEGER               :: Node,Nel,M,n
!-------------------------------------------------------
!		Read Node Co-ordinates from Inputfile
!-------------------------------------------------------
DO Node=1,Nodes
 READ(11,*) (xP(M,Node),M=1,Cdim)
END DO
DO Node=1,Nodes,Nodes-1
 WRITE(12,'(A5,I5,A8,3F8.2)') 'Node ',Node,&
         '  Coor  ',(xP(M,Node),M=1,Cdim)
END DO

!-------------------------------------------------------
!		Read Incidences from Inputfile
!-------------------------------------------------------
 WRITE(12,*)''
 WRITE(12,*)'Incidences: '
 WRITE(12,*)''
        DO Nel=1,Maxe
           READ(11,*) (Inci(Nel,n),n=1,Nodel)
        END DO
DO Nel = 1 , Maxe, Maxe - 1
 WRITE(12,'(A3,I5,A8,4I5)')'EL ',Nel,'  Inci  ',Inci(Nel,:)
END DO 
CALL FLUSH(12)
RETURN
END SUBROUTINE Geomin

LOGICAL FUNCTION Match(Inci1,Inci2)
!-----------------------------------
!    Returns a value of TRUE if the incidences 
!    Inci1 and Inci2 match
!------------------------------------
IMPLICIT NONE
INTEGER, INTENT (IN) :: Inci1(:) !  1. incidence array
INTEGER, INTENT (IN) :: Inci2(:) !  2. incidence array
INTEGER :: Nodes,Node,N1,N2,Ncount
Nodes= UBOUND(Inci1,1)
Ncount= 0
Node_loop1: &
DO N1=1,Nodes
	Node= Inci1(n1)
	Node_loop2: &
	DO N2=1,Nodes
		IF(Node == Inci2(n2)) Ncount= Ncount+1
	END DO &
	Node_loop2
END DO &
Node_loop1
IF(Ncount == Nodes) THEN
 Match= .TRUE.
ELSE
 Match= .FALSE.
END IF
END FUNCTION Match

SUBROUTINE Destination(Isym,Ndest,Ldest,xP,Inci,Ndofs,nodes,N_dof,Nodel,Maxe)
!-------------------------------------------------------------------
! Determine Node destination vector and Element destination vector  
!-------------------------------------------------------------------
IMPLICIT NONE
REAL(iwp), INTENT (IN)         ::      xP(:,:)        !       Node co-ordinates
INTEGER, INTENT (IN OUT)  ::      Ndest(:,:)     ! Node destination vector
INTEGER, INTENT (IN OUT)  ::      Ldest(:,:)     ! Element destination vector
INTEGER, INTENT (IN OUT)  ::      Ndofs          ! DoF's of System
INTEGER, INTENT (IN)      ::      Inci(:,:)      !       Element Incidences
INTEGER, INTENT (IN)      ::      Isym,nodes,N_dof,Nodel,Maxe 
INTEGER                   ::      k,m,n,Nel,l
!---------------------------------------------------
!     Determine Node destination vector
!			Set Ndest == 0 if Point is on a symmetry plane
!---------------------------------------------------
!		no symmetry
IF(Isym == 0) THEN
	k=1
	Nodes0:	DO m=1, nodes
                                                DO n=1, N_dof                     
							Ndest(m,n)= k
							k=k+1
						END DO
					END DO Nodes0
!		y-z symmetry
ELSE IF(Isym == 1) THEN
	k=1
	Nodes1:	DO m=1, nodes
						IF(xP(1,m) == 0.0)THEN
							Ndest(m,1)= 0
						ELSE
							Ndest(m,1)= k
							k=k+1
						END IF
                                                DO n=2, N_dof
							Ndest(m,n)= k
							k=k+1
						END DO
					END DO Nodes1
!		x-z and y-z symmetry
ELSE IF(Isym == 2) THEN
	k=1
	Nodes2:	DO m=1, nodes
						IF(xP(1,m) == 0.0)THEN
							Ndest(m,1)= 0
						ELSE
							Ndest(m,1)= k
							k=k+1
						END IF
						IF(xP(2,m) == 0.0)THEN
							Ndest(m,2)= 0
						ELSE
							Ndest(m,2)= k
							k=k+1
						END IF
                                                IF (N_dof == 3) THEN
							Ndest(m,3)= k
							k=k+1
						END IF
					END DO Nodes2
!		x-y, x-z and y-z symmetry
ELSE
	k=1
	Nodes3: DO m=1, nodes
						IF(xP(1,m) == 0.0)THEN
							Ndest(m,1)= 0
						ELSE
							Ndest(m,1)= k
							k=k+1
						END IF
						IF(xP(2,m) == 0.0)THEN
							Ndest(m,2)= 0
						ELSE
							Ndest(m,2)= k
							k=k+1
						END IF
						IF(xP(3,m) == 0.0)THEN
							Ndest(m,3)= 0
						ELSE
							Ndest(m,3)= k
							k=k+1
						END IF
					END DO Nodes3
END IF						
Ndofs= k-1                                                     ! DoF's of System
!------------------------------------------
!     Determine Element destination vector
!---------------------------------------------
Elements:&
DO Nel=1,Maxe
	DO n=1,Nodel
                k= (n-1)*N_dof+1
                l= n*N_dof
		Ldest(Nel,k:l)=Ndest(Inci(Nel,n),:)	
	END DO		   
END DO &
Elements
END SUBROUTINE Destination

SUBROUTINE Scal(E,xP,Elres_u,Elres_t,Cdim,Scad,Scat)
IMPLICIT NONE
REAL(iwp), INTENT (INOUT)  ::      E                       !       Youngsmodulus
REAL(iwp), INTENT (INOUT)  ::      xP(:,:)                 !       Node co-ordinates
REAL(iwp), INTENT (INOUT)  ::      Elres_u(:,:)            ! Element results , u
REAL(iwp), INTENT (INOUT)  ::      Elres_t(:,:)            ! Element results , t
REAL(iwp), INTENT (OUT)    ::      Scad  
REAL(iwp), INTENT (OUT)    ::      Scat
INTEGER, INTENT(IN)   ::      Cdim                    !   Cartesian Dimension
REAL(iwp)                  ::      xmax(Cdim),xmin(Cdim),delta_x(Cdim)

!-------------------------------------------------------
!		Determine Scalefactor for Tractions
!		Scat ... 1/E
!-------------------------------------------------------
Scat= 1./E             ! Scalefactor for Tractions
E=1.0                  ! Scaled Youngsmodulus
Elres_t=Elres_t*Scat   ! Scaled prescribed Tractions by Scat
!----------------------------------------------------------
!		Determine Scalefactor for Displacements
!		Scad ... max. Distance in any co-ordinate direction
!----------------------------------------------------------
xmax(1)= MAXVAL(xp(1,:))
xmax(2)= MAXVAL(xp(2,:))
IF(Cdim == 3)xmax(3)= MAXVAL(xp(3,:))
xmin(1)= MINVAL(xp(1,:))
xmin(2)= MINVAL(xp(2,:))
IF(Cdim == 3)xmin(3)= MINVAL(xp(3,:))
delta_x= xmax - xmin
Scad= MAXVAL(delta_x)    !       Scad ... Scalefactor for Displacements 
xP=xP/Scad               ! Scaled Node co-ordinates
Elres_u=Elres_u/Scad     ! Scaled prescribed Displacements by Scad
END SUBROUTINE Scal




SUBROUTINE bicgstab(a,b,x,n,tol,its)
! bicgstab cast as a subroutine
! on input x is the starting guess , on output the solution
IMPLICIT NONE
REAL(iwp),INTENT(IN)::a(:,:),b(:),tol ; REAL(iwp),INTENT(IN OUT)::x(:)
INTEGER,INTENT(IN):: n,its    ;  INTEGER:: iters ; LOGICAL:: converged
REAL(iwp)::alpha,beta,rho0,rho1,w,zero=.0_iwp,one=1._iwp
REAL(iwp)::v(n),r(n),r0_hat(n),p(n),s(n),t(n)  ! local arrays
    iters = 0
     r = b - MATMUL(a,x)   ;   r0_hat = r
     rho0 = one ; alpha = one; w = one; v = zero; p = zero
     rho1 = DOT_PRODUCT(r0_hat,r)
!     WRITE(12,'(/,A)')"First few iterations" 
       DO
            iters = iters + 1 ; converged = norm(r) < tol * norm(b)
            IF(iters==its.OR. converged) EXIT
            beta = (rho1/rho0)*(alpha/w)
            p = r + beta*(p - w*v) ; v = MATMUL(a,p)
            alpha = rho1/DOT_PRODUCT(r0_hat,v)
            s = r - alpha*v ; t = MATMUL(A,s)
            w=DOT_PRODUCT(t,s)/DOT_PRODUCT(t,t) 
            rho0 = rho1; rho1 = -w*DOT_PRODUCT(r0_hat,t)
            x = x + alpha*p + w*s ; r = s - w*t
!            IF(iters<5) WRITE(12,'(3E12.4)') x
       END DO
WRITE(12,'(/,A,I5,A,/)')"It took BiCGSTAB ",iters,"  iterations to converge"
END SUBROUTINE bicgstab

 SUBROUTINE bicgstab_l(a,b,neq,x,x0,tol,limit,ell,kappa)
 IMPLICIT NONE
 REAL(iwp),INTENT(IN)::a(:,:),b(:),x0,tol,kappa
 REAL(iwp),INTENT(OUT)::x(:); INTEGER,INTENT(IN)::neq,limit,ell
 INTEGER:: iters,j,k
 REAL(iwp)::alpha,beta,rho,gama,omega,norm_r,r0_norm,error,          &
            cosine,NGamma0,NGamma1,one=1._iwp,zero=.0_iwp
 LOGICAL:: converged
 REAL(iwp)::s(ell+1),GG(ell+1,ell+1),HH(ell-1,ell-1),p(ell-1),q(ell-1),     &
            Gamma0(ell+1),Gamma1(ell+1),Gamma(ell+1),                 &
            rt(neq),y(neq),y1(neq),r(neq,ell+1),u(neq,ell+1)
!      initialisation phase
    x = x0   ;    y = x   ; y1 = MATMUL(a,y) ; y=y1; rt = b - y
    r=zero ; r(:,1) = rt   ;  u = zero  ; gama = one  ; omega=one ;  k = 0
    norm_r = norm(rt)   ;   r0_norm = norm_r   ;  error = one   ;  iters = 0
!     bicgstab(ell)  iterations
       iterations : DO
            iters = iters + 1    ;            converged = error  < tol
            IF(iters==limit.OR. converged) EXIT 
            gama = - omega*gama  ;  y = r(:,1)
            DO j = 1 , ell
               rho = DOT_PRODUCT(rt,y)  ;  beta = rho/gama
               u(:,1:j) = r(:,1:j) - beta * u(:,1:j)   ;      y = u(:,j)
               y1 = MATMUL(a,y); y=y1; u(:,j+1) = y 
               gama = DOT_PRODUCT(rt,y) ; alpha = rho/gama;x=x+ alpha * u(:,1)
               r(:,1:j) = r(:,1:j) - alpha * u(:,2:j+1)    
               y = r(:,j)  ;  y1 = MATMUL(a,y) ; y=y1 ; r(:,j+1) = y
            END DO
            GG = MATMUL(TRANSPOSE(r),r); HH=-GG(2:ell,2:ell); HH=inverse(HH)   
            p = MATMUL(HH,GG(2:ell,1))     ;   q = MATMUL(HH,GG(2:ell,ell+1))
            Gamma0(1) = one; Gamma0(ell+1) = zero; Gamma0(2:ell) = p
            Gamma1(1) = zero; Gamma1(ell+1) = one; Gamma1(2:ell) = q
            NGamma0   = DOT_PRODUCT(Gamma0,matmul(GG,Gamma0))
            NGamma1   = DOT_PRODUCT(Gamma1,matmul(GG,Gamma1))
            omega     = DOT_PRODUCT(Gamma0,matmul(GG,Gamma1))
            cosine = ABS(omega)/SQRT(ABS(NGamma0*NGamma1)); omega=omega/NGamma1
            IF(cosine<kappa) omega = (kappa/cosine) * omega
            Gamma = Gamma0 - omega * Gamma1
            s(1:ell) = Gamma(2:ell+1)       ;  s(ell+1) = zero 
            x = x - MATMUL(r,s);r(:,1)=MATMUL(r,Gamma);u(:,1)=MATMUL(u,Gamma)
            norm_r = norm(r(:,1))  ;  error = norm_r/r0_norm    ;  k = k + 1   
        END DO iterations
  WRITE(12,'(/,A,I5,A,/)')"It took BiCGSTAB_L ",iters," iterations to converge"
 END SUBROUTINE bicgstab_l


 SUBROUTINE rhs_and_diag(Rhs,DTe,DUe,Ldest,BCode,Ncode &
           ,Elres_u,Elres_te,Diag,Ndofe,N_dof,Nodel,Fac)
!---------------------------------------------
!  Assembles Element contributions DTe , DUe
!  into vector Rhs and sums off-diagonal coefficients 
!  for the computation of diagonal coefficients
!---------------------------------------------
IMPLICIT NONE
REAL(iwp)            :: Rhs(:)    ! Global vector
REAL(iwp), INTENT(IN):: DTe(:,:),DUe(:,:)  ! Element arrays
INTEGER , INTENT(IN)    :: LDest(:)          ! Element destination vector
INTEGER , INTENT(IN)            :: BCode(:)  ! Boundary code(local)
INTEGER , INTENT(IN)            :: NCode(:)  ! Boundary code (global) 
INTEGER , INTENT(IN)            :: Ndofe     ! D.o.F´s / Elem
INTEGER , INTENT(IN)            :: N_dof      ! D.o.F´s / Node
INTEGER , INTENT(IN)            :: Nodel     ! Nodes/Element
REAL(iwp) , INTENT(IN)               :: Elres_u(:)       ! vector u for element
REAL(iwp) , INTENT(IN)               :: Elres_te(:)      ! vector t for element
REAL(iwp) , INTENT(IN)               :: Fac(:)           ! Mult. factors for symmetry  
REAL(iwp) :: Diag(:,:)          ! Array containing diagonal coeff of DT
INTEGER   :: m,n,k,l,Ncol
DoF_per_Element:&
DO m=1,Ndofe  
	Ncol=Ldest(m)      !   Column number 
	IF(BCode(m) == 0) THEN    !   Neumann BC
		Rhs(:) = Rhs(:) + DUe(:,m)*Elres_te(m)*Fac(m)
	!     The assembly of dTe depends on the global BC
                IF (.not.(NCode(Ldest(m)) == 0 .and. Ncol /= 0)) THEN 
			Rhs(:) = Rhs(:) - DTe(:,m) * Elres_u(m)*Fac(m)
		END IF
	END IF
	IF(BCode(m) == 1) THEN   !   Dirichlet BC
		Rhs(:)= Rhs(:) - DTe(:,m) * Elres_u(m)*Fac(m)
	END IF
END DO &
DoF_per_Element
!				Sum of off-diagonal coefficients
DO n=1,Nodel
        DO k=1,N_dof
                l=(n-1)*N_dof+k
		Diag(:,k)= Diag(:,k) - DTe(:,l)
	END DO
END DO
	RETURN
END SUBROUTINE rhs_and_diag


 SUBROUTINE form_lhs(lhs,DTe,DUe,Ldeste,BCode,Ncode,Ndofe,Fac)
!---------------------------------------------
!  Finds Element contributions DTe , DUe to lhs in right combination
!---------------------------------------------
IMPLICIT NONE
REAL(iwp)            :: lhs(:,:)    ! lhs array
REAL(iwp), INTENT(IN):: DTe(:,:),DUe(:,:)  ! Element arrays
INTEGER , INTENT(IN)    :: LDeste(:)         ! Element destination vector
INTEGER , INTENT(IN)            :: BCode(:)  ! Boundary code(local)
INTEGER , INTENT(IN)            :: NCode(:)  ! Boundary code (global) 
INTEGER , INTENT(IN)            :: Ndofe     ! D.o.F´s / Elem
REAL(iwp) , INTENT(IN)               :: Fac(:)           ! Mult. factors for symmetry  
INTEGER   :: n,Ncol,m  ;  lhs = 0.0_iwp  
DoF_per_Element:&
DO m=1,Ndofe  
        Ncol=Ldeste(m)      !   Column number 
	IF(BCode(m) == 0) THEN    !   Neumann BC
                IF (NCode(Ldeste(m)) == 0 .and. Ncol /= 0) THEN 
                     lhs(:,m) =  DTe(:,m)*Fac(m)  
		END IF
	END IF
	IF(BCode(m) == 1) THEN   !   Dirichlet BC
                     lhs(:,m) =  - DUe(:,m)*Fac(m)    
	END IF
END DO &
DoF_per_Element
	RETURN
END SUBROUTINE form_lhs

SUBROUTINE get_km(Cdim,Nel,y,active_diag,g,qmul,km)
! Gets km and qmul via g
REAL(iwp),INTENT(IN)::y(:),active_diag(:,:)
REAL(iwp),INTENT(OUT)::qmul(:),km(:,:); INTEGER,INTENT(OUT)::g(:)
INTEGER,INTENT(IN)::Cdim,nel
        IF(Cdim==2) THEN
            g = (/2*Nel-1,2*Nel/);  qmul = y(g)
            km(1,:) = active_diag(2*Nel-1,:); km(2,:)=active_diag(2*Nel,:)
        ELSE IF(Cdim==3) THEN
            g = (/3*Nel-2,3*Nel-1,3*nel/);  qmul = y(g)
            km(1,:) = active_diag(3*Nel-2,:); km(2,:)=active_diag(3*Nel-1,:)
            km(3,:) = active_diag(3*Nel,:)
        END IF
END SUBROUTINE get_km

END MODULE bem_lib 
