MODULE BICG

  !/****h* /bicg
  !*  NAME
  !*    MODULE: bicg
  !*  SYNOPSIS
  !*    Usage:      USE bicg
  !*  FUNCTION
  !*    Contains subroutines required by BiCGSTAB(l)
  !*    These subroutines are parallel and require MPI.
  !*    
  !*    Subroutine             Purpose
  !*
  !*    FORM_S                 Forms the s vector in bicgstab(l)
  !*  AUTHOR
  !*    L. Margetts
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/
  
  USE precision
  USE maths
  USE gather_scatter

  IMPLICIT NONE
  
  CONTAINS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE FORM_S(GG,ell,kappa,omega,Gamma,s)

    !/****f* bicg/form_s
    !*  NAME
    !*    SUBROUTINE: form_s
    !*  SYNOPSIS
    !*    Usage:      CALL form_s(GG,ell,kappa,omega,gamma,s)
    !*  FUNCTION
    !*    Forms the s vector in bicgstab(l)
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    GG(:,:)                  : Real
    !*                             : GG Array
    !*
    !*    kappa                    : Real
    !*                             : Variable
    !*
    !*    ell                      : Integer
    !*                             : Variable
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    omega                    : Real
    !*                             : Variable
    !*    
    !*    Gamma(:)                 : Real   
    !*                             : Vector
    !*
    !*    s(:)                     : Real
    !*                             : Vector
    !*
    !*  AUTHOR
    !*    I.M. Smith
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    REAL(iwp),INTENT(IN)  :: GG(:,:),kappa
    INTEGER, INTENT(IN)   :: ell
    REAL(iwp),INTENT(OUT) :: omega,Gamma(:),s(:)
    REAL(iwp)             :: HH(ell-1,ell-1),Gamma0(ell+1),p(ell-1),q(ell-1),&
                             Gamma1(ell+1),NGamma0,NGamma1,cosine
 
    HH            = -GG(2:ell,2:ell)
 
    call invert(HH)  
 
    p             = matmul(HH,GG(2:ell,1))  
    q             = matmul(HH,GG(2:ell,ell+1))
    Gamma0(1)     = 1.0_iwp
    Gamma0(ell+1) = 0.0_iwp
    Gamma0(2:ell) = p
    Gamma1(1)     = 0.0_iwp
    Gamma1(ell+1) = 1.0_iwp
    Gamma1(2:ell) = q
    NGamma0       = dot_product(Gamma0,matmul(GG,Gamma0))
    NGamma1       = dot_product(Gamma1,matmul(GG,Gamma1))
    omega         = dot_product(Gamma0,matmul(GG,Gamma1))
    cosine        = abs(omega)/sqrt(abs(NGamma0*NGamma1))
    omega         = omega/NGamma1
 
    if(cosine<kappa) omega = (kappa/cosine) * omega
 
    Gamma         = Gamma0 - omega * Gamma1
    s(1:ell)      = Gamma(2:ell+1)       
    s(ell+1)      = 0.0_iwp
   
    RETURN
    
  END SUBROUTINE FORM_S
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

END MODULE BICG
