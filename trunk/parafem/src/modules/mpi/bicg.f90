MODULE BICG

  !/****h* /bicg
  !*  NAME
  !*    MODULE: bicg
  !*  SYNOPSIS
  !*    Usage:      USE bicg
  !*  FUNCTION
  !*    Contains subroutines related to BiCGSTAB(l)
  !*    These subroutines are parallel and require MPI.
  !*    
  !*    Subroutine             Purpose
  !*
  !*    BICGSTABL_P            Parallel BiCGSTAB(l) solver
  !*  AUTHOR
  !*    L. Margetts
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2015 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/
  
  USE precision
  USE global_variables
  USE maths
  USE gather_scatter

  IMPLICIT NONE
  
  CONTAINS


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE BICGSTABL_P(b_pp,cjiters,cjits,cjtol,ell,ieq_start,no_pp,pmul_pp,  &
                       store_pp,storke_pp,x_pp)

USE precision
USE global_variables
USE maths

IMPLICIT NONE

INTEGER, INTENT(IN)      :: cjits,ell,ieq_start,no_pp(:)
INTEGER, INTENT(OUT)     :: cjiters
REAL(iwp), INTENT(IN)    :: b_pp(:),cjtol,store_pp(:),storke_pp(:,:,:)
REAL(iwp), INTENT(INOUT) :: pmul_pp(:,:),x_pp(:)

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

 INTEGER                :: i,j,k,l,iel
 INTEGER                :: fixed_freedoms_pp
 INTEGER                :: neq_pp,nels_pp,ntot ! any global?
 REAL(iwp),PARAMETER    :: zero = 0.0_iwp
 REAL(iwp),PARAMETER    :: one  = 1.0_iwp
 REAL(iwp)              :: gama,omega,norm_r,r0_norm,error,rho1,beta
 REAL(iwp),ALLOCATABLE  :: y_pp(:),y1_pp(:),rt_pp(:),r_pp(:,:),u_pp(:,:)
 REAL(iwp),ALLOCATABLE  :: utemp_pp(:,:),GG(:,:),s(:),Gamma(:)
 
 LOGICAL                :: cj_converged = .false.

!------------------------------------------------------------------------------
! 2. Allocate local arrays
!------------------------------------------------------------------------------

 neq_pp            = UBOUND(x_pp,1)
 nels_pp           = UBOUND(storke_pp,3)
 fixed_freedoms_pp = UBOUND(no_pp,1)
 ntot              = UBOUND(storke_pp,1)
 
 ALLOCATE(y_pp(neq_pp))           ; y_pp     = zero
 ALLOCATE(y1_pp(neq_pp))          ; y1_pp    = zero
 ALLOCATE(rt_pp(neq_pp))          ; rt_pp    = zero
 ALLOCATE(r_pp(neq_pp,ell+1))     ; r_pp     = zero
 ALLOCATE(u_pp(neq_pp,ell+1))     ; u_pp     = zero
 ALLOCATE(utemp_pp(ntot,nels_pp)) ; utemp_pp = zero
 ALLOCATE(GG(ell+1,ell+1))        ; GG       = zero
 ALLOCATE(s(ell+1))               ; s        = zero
 ALLOCATE(Gamma(ell+1))           ; Gamma    = zero
 
!------------------------------------------------------------------------------
! 3. Initialise BiCGSTAB(l)
!------------------------------------------------------------------------------

!IF(iters==1) x_pp=x0; pmul_pp=zero; y1_ppn=zero; y_pp=x_pp
!This works in context of p126, but not sure if it is general

!  x_pp=x0 - now an input value
   
   pmul_pp=zero; y1_pp=zero; y_pp=x_pp

   CALL gather(y_pp,pmul_pp)
   elements_3: DO iel=1,nels_pp  
     utemp_pp(:,iel)=MATMUL(storke_pp(:,:,iel),pmul_pp(:,iel))
   END DO elements_3
   CALL scatter(y1_pp,utemp_pp)
  
   DO i=1,fixed_freedoms_pp; k=no_pp(i)-ieq_start+1
     y1_pp(k)=y_pp(k)*store_pp(k)
   END DO; y_pp=y1_pp; rt_pp=b_pp-y_pp; r_pp=zero; r_pp(:,1)=rt_pp
   
   u_pp=zero; gama=one; omega=one; k=0; norm_r=norm_p(rt_pp)
   r0_norm=norm_r; error=one; cjiters=0
   
!------------------------------------------------------------------------------
! 4. BiCGSTAB(ell) iterations 
!------------------------------------------------------------------------------

   bicg_iterations: DO; cjiters=cjiters+1   
     cj_converged=error<cjtol; IF(cjiters==cjits.OR.cj_converged) EXIT
     gama=-omega*gama; y_pp=r_pp(:,1)
     DO j=1,ell
       rho1=DOT_PRODUCT_P(rt_pp,y_pp); beta=rho1/gama
       u_pp(:,1:j)=r_pp(:,1:j)-beta*u_pp(:,1:j)
       pmul_pp=zero; y_pp=u_pp(:,j); y1_pp=zero; CALL gather(y_pp,pmul_pp)
       elements_4: DO iel=1,nels_pp
         utemp_pp(:,iel)=MATMUL(storke_pp(:,:,iel),pmul_pp(:,iel))
       END DO elements_4; CALL scatter(y1_pp,utemp_pp)
       DO i=1,fixed_freedoms_pp; l=no_pp(i)-ieq_start+1
         y1_pp(l)=y_pp(l)*store_pp(l)
       END DO; y_pp=y1_pp; u_pp(:,j+1)=y_pp
       gama=DOT_PRODUCT_P(rt_pp,y_pp); alpha=rho1/gama
       x_pp=x_pp+alpha*u_pp(:,1)
       r_pp(:,1:j)=r_pp(:,1:j)-alpha*u_pp(:,2:j+1)
       pmul_pp=zero; y_pp=r_pp(:,j); y1_pp=zero; CALL gather(y_pp,pmul_pp)
       elements_5: DO iel=1,nels_pp
         utemp_pp(:,iel)=MATMUL(storke_pp(:,:,iel),pmul_pp(:,iel))
       END DO elements_5; CALL scatter(y1_pp,utemp_pp)
       DO i=1,fixed_freedoms_pp; l=no_pp(i)-ieq_start+1
         y1_pp(l)=y_pp(l)*store_pp(l)
       END DO; y_pp=y1_pp; r_pp(:,j+1)=y_pp
     END DO
     DO i=1,ell+1; DO j=1,ell+1
       GG(i,j)=DOT_PRODUCT_P(r_pp(:,i),r_pp(:,j))
     END DO; END DO
     CALL form_s(GG,ell,kappa,omega,Gamma,s) 
     x_pp=x_pp-MATMUL(r_pp,s); r_pp(:,1)=MATMUL(r_pp,Gamma)
     u_pp(:,1)=MATMUL(u_pp,Gamma); norm_r=norm_p(r_pp(:,1))
     error=norm_r/r0_norm; k=k+1
   END DO bicg_iterations            
  
!------------------------------------------------------------------------------
! 5. Clean up
!------------------------------------------------------------------------------

 DEALLOCATE(y_pp,y1_pp,rt_pp,r_pp,u_pp,utemp_pp,GG,s,Gamma)
 
 RETURN

END SUBROUTINE BICGSTABL_P

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END MODULE BICG
