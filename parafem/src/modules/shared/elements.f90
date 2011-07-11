MODULE ELEMENTS

  !/****h* /elements
  !*  NAME
  !*    MODULE: elements
  !*  SYNOPSIS
  !*    Usage:      USE elements
  !*  FUNCTION
  !*    Contains subroutines used for computing element level matrices
  !*    
  !*    Subroutine             Purpose
  !*
  !*    SHAPE_FUN              Computes the shape functions
  !*    SHAPE_DER              Compute the derivatives of the shape functions
  !*    BEEMAT                 Compute the B matrix
  !*    SAMPLE                 Returns local coords of the integrating points
  !*    ECMAT                  Returns the consistent mass matrix
  !*    DEEMAT                 Compute the D matrix
  !*    SHAPE_DERIVATIVES      Compute the derivatives of the shape functions
  !* 
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*    L. Margetts
  !*  COPYRIGHT
  !*    2004-2011 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  
  !*  LM not sure about module name or module content
  !*
  !*/
  
  USE precision

  CONTAINS                                  

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE SHAPE_FUN(fun,points,i)

    !/****f* elements/shape_fun
    !*  NAME
    !*    SUBROUTINE: shape_fun
    !*  SYNOPSIS
    !*    Usage:      CALL shape_fun(fun,points,i)
    !*  FUNCTION
    !*    Compute the value of the shape functions at a Gauss point.
    !*    Source: "Programming the Finite Element Method"
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
    !*    The following argument has the INTENT(INOUT) attribute:
    !*
    !*    fun(nod)           : Real
    !*                       : Value of the shape functions at a 
    !*                         Gauss point
    !*  AUTHOR
    !*    I.M. Smith
    !*    D.V. Griffiths
    !*    L. Margetts
    !*  COPYRIGHT
    !*    (c) University of Manchester 2004-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)     :: i
    REAL(iwp),INTENT(IN)    :: points(:,:)
    REAL(iwp),INTENT(INOUT) :: fun(:)
    REAL(iwp)               :: eta,xi,etam,etap,xim,xip,zetam,zetap,c1,c2,c3    
    REAL(iwp),PARAMETER     :: one=1.0_iwp,two=2.0_iwp,d8=8.0_iwp
    REAL(iwp)               :: t1,t2,t3,t4,t5,t6,t7,t8,t9
    REAL(iwp)               :: zeta,xi0,eta0,zeta0
    INTEGER                 :: xii(20),etai(20),zetai(20),l,ndim,nod
    
    ndim = UBOUND(points,2)
    nod  = UBOUND(fun,1)  
    
    SELECT CASE (ndim)
      
      CASE(1) ! one dimensional cases
      
      xi = points(i,1)
      
        SELECT CASE(nod)
        
          CASE(2)
 
          t1 = -1.-xi 
          t2 =  1.-xi
          
          fun(1) =  t2/2. 
          fun(2) = -t1/2.
           
          CASE(3)
          
          t1 =  -1.-xi 
          t2 =  -xi 
          t3 =  1.-xi
          
          fun(1) = t2*t3/2. 
          fun(2) = -t1*t3 
          fun(3) = t1*t2/2.
          
          CASE(4)
          
          t1 = -1.-xi 
          t2 = -1./3.-xi 
          t3 = 1./3.-xi 
          t4 = 1.-xi
          
          fun(1) =  t2*t3*t4*9./16.  
          fun(2) = -t1*t3*t4*27./16.
          fun(3) =  t1*t2*t4*27./16. 
          fun(4) = -t1*t2*t3*9./16.
          
          CASE(5)
          
          t1 = -1.-xi 
          t2 = -0.5-xi 
          t3 = -xi 
          t4 = 0.5-xi 
          t5 = 1.-xi
          
          fun(1) =  t2*t3*t4*t5*2./3. 
          fun(2) = -t1*t3*t4*t5*8./3.
          fun(3) =  t1*t2*t4*t5*4. 
          fun(4) = -t1*t2*t3*t5*8./3.
          fun(5) =  t1*t2*t3*t4*2./3.
          
          CASE DEFAULT
          
          print*,"wrong number of nodes in shape_fun"
          
        END SELECT
        
        CASE(2) ! two dimensional cases
        
        c1   = points(i,1)
        c2   = points(i,2)
        c3   = 1.-c1-c2 
        xi   = points(i,1)
        eta  = points(i,2)
        etam =.25*(1.-eta)
        etap =.25*(1.+eta)
        xim  =.25*(1.-xi)
        xip  =.25*(1.+xi)
          
        SELECT CASE(nod)
        
          CASE(3)
          
          fun = (/c1,c3,c2/)  
          
          CASE(6)
          
          fun(1) = (2.*c1-1.)*c1 
          fun(6) = 4.*c1*c2 
          fun(5) = (2.*c2-1.)*c2
          fun(4) = 4.*c2*c3      
          fun(3) = (2.*c3-1.)*c3 
          fun(2) = 4.*c3*c1
            
          CASE(15)
          
          t1 = c1-.25  
          t2 = c1-.5 
          t3 = c1-.75   
          t4 = c2-.25
          t5 = c2-.5   
          t6 = c2-.75 
          t7 = c3-.25  
          t8 = c3-.5 
          t9 = c3-.75
           
          fun(1)  = 32./3.*c1*t1*t2*t3   
          fun(12) = 128./3.*c1*c2*t1*t2
          fun(11) = 64.*c1*c2*t1*t4     
          fun(10) = 128./3.*c1*c2*t4*t5
          fun(9)  = 32./3.*c2*t4*t5*t6   
          fun(8)  = 128./3.*c2*c3*t4*t5
          fun(7)  = 64.*c2*c3*t4*t7      
          fun(6)  = 128./3.*c2*c3*t7*t8
          fun(5)  = 32./3.*c3*t7*t8*t9   
          fun(4)  = 128./3.*c3*c1*t7*t8
          fun(3)  = 64.*c3*c1*t1*t7      
          fun(2)  = 128./3.*c3*c1*t1*t2
          fun(13) = 128.*c1*c2*t1*c3    
          fun(15) = 128.*c1*c2*c3*t4
          fun(14) = 128.*c1*c2*c3*t7      
          
          CASE(4)
          
          fun = (/4.*xim*etam,4.*xim*etap,4.*xip*etap,4.*xip*etam/)
          
          CASE(8)
          
          fun = (/4.*etam*xim*(-xi-eta-1.),32.*etam*xim*etap,&
                  4.*etap*xim*(-xi+eta-1.),32.*xim*xip*etap, &
                  4.*etap*xip*(xi+eta-1.), 32.*etap*xip*etam,&
                  4.*xip*etam*(xi-eta-1.), 32.*xim*xip*etam/)
            
          CASE(9)
          
          etam = eta - 1.
          etap = eta + 1.
          xim  = xi - 1.
          xip  = xi + 1.
          
          fun = (/.25*xi*xim*eta*etam,-.5*xi*xim*etap*etam,&
                 .25*xi*xim*eta*etap,-.5*xip*xim*eta*etap,&
                 .25*xi*xip*eta*etap,-.5*xi*xip*etap*etam,&
                 .25*xi*xip*eta*etam,-.5*xip*xim*eta*etam,xip*xim*etap*etam/)
          
          CASE DEFAULT
          
          print*,"wrong number of nodes in shape_fun"
          
        END SELECT
        
        CASE(3) ! three dimensional cases
        
        xi    = points(i,1)
        eta   = points(i,2)
        zeta  = points(i,3)
        etam  = one  - eta 
        xim   = one  - xi  
        zetam = one  - zeta
        etap  = eta  + one 
        xip   = xi   + one 
        zetap = zeta + one
         
        SELECT CASE(nod)
        
          CASE(4)
        
          fun(1)  = xi   
          fun(2)  = eta 
          fun(3)  = zeta 
          fun(4)  = one - fun(1) - fun(2) - fun(3)
          
          CASE(8)
        
          fun = (/0.125_iwp*xim*etam*zetam,0.125_iwp*xim*etam*zetap,          &
                  0.125_iwp*xip*etam*zetap,0.125_iwp*xip*etam*zetam,          &
                  0.125_iwp*xim*etap*zetam,0.125_iwp*xim*etap*zetap,          &
                  0.125_iwp*xip*etap*zetap,0.125_iwp*xip*etap*zetam/)
          
          CASE(14) !type 6 element
  
          fun(1) = (xi*eta+xi*zeta+two*xi+eta*zeta+two*eta+two*zeta+two)* &   
                   (xi-one)*(eta-one)*(zeta-one)/d8
          fun(2) =-(xi*eta-xi*zeta+two*xi-eta*zeta+two*eta-two*zeta+two)* &   
                   (xi-one)*(eta-one)*(zeta+one)/d8
          fun(3) =-(xi*eta-xi*zeta+two*xi+eta*zeta-two*eta+two*zeta-two)* &   
                   (xi+one)*(eta-one)*(zeta+one)/d8
          fun(4) = (xi*eta+xi*zeta+two*xi-eta*zeta-two*eta-two*zeta-two)* &   
                   (xi+one)*(eta-one)*(zeta-one)/d8
          fun(5) =-(xi+one)*(xi-one)*(eta-one)*(zeta+one)*(zeta-one)/two
          fun(6) =-(xi-one)*(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
          fun(7) = (xi+one)*(xi-one)*(eta+one)*(eta-one)*(zeta+one)/two
          fun(8) = (xi+one)*(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
          fun(9) =-(xi+one)*(xi-one)*(eta+one)*(eta-one)*(zeta-one)/two
          fun(10)= (xi*eta-xi*zeta-two*xi+eta*zeta+two*eta-two*zeta-two)*   &
                   (xi-one)*(eta+one)*(zeta-one)/d8
          fun(11)=-(xi*eta+xi*zeta-two*xi-eta*zeta+two*eta+two*zeta-two)*   &
                   (xi-one)*(eta+one)*(zeta+one)/d8
          fun(12)=-(xi*eta+xi*zeta-two*xi+eta*zeta-two*eta-two*zeta+two)*   &
                   (xi+one)*(eta+one)*(zeta+one)/d8
          fun(13)= (xi*eta-xi*zeta-two*xi-eta*zeta-two*eta+two*zeta+two)*   &
                   (xi+one)*(eta+one)*(zeta-one)/d8
          fun(14)= (xi+one)*(xi-one)*(eta+one)*(zeta+one)*(zeta-one)/two
        
          CASE(20)
          
          xii   = (/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
          etai  = (/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
          zetai = (/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
            
          DO l=1,20
            xi0   = xi*xii(l)
            eta0  = eta*etai(l)
            zeta0 = zeta*zetai(l)
            IF(l==4.or.l==8.or.l==16.or.l==20) THEN
              fun(l)=.25*(1.-xi*xi)*(1.+eta0)*(1.+zeta0)
            ELSE IF(l>=9.and.l<=12)THEN
              fun(l)=.25*(1.+xi0)*(1.-eta*eta)*(1.+zeta0)
            ELSE IF(l==2.or.l==6.or.l==14.or.l==18)THEN
              fun(l)=.25*(1.+xi0)*(1.+eta0)*(1.-zeta*zeta)
            ELSE
              fun(l)=.125*(1.+xi0)*(1.+eta0)*(1.+zeta0)*(xi0+eta0+zeta0-2)
            END IF
          END DO
          
          CASE DEFAULT
          
          print*,"wrong number of nodes in shape_fun"
          
        END SELECT
          
        CASE DEFAULT
        
        print*,"wrong number of dimensions in shape_fun"  
      
      END SELECT
      
    RETURN
    
  END SUBROUTINE SHAPE_FUN 
   
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
 
  SUBROUTINE SHAPE_DER(der,points,i)
  
    !/****f* elements/shape_der
    !*  NAME
    !*    SUBROUTINE: shape_der
    !*  SYNOPSIS
    !*    Usage:      CALL shape_der(der,points,i)
    !*  FUNCTION
    !*    Compute the derivatives of the shape functions at a Gauss point.
    !*    Source: "Programming the Finite Element Method"
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    i                  : Integer
    !*                       : Gauss point number
    !*
    !*    points(nip,ndim)   : Real
    !*                       : Gauss points coordinates at the reference 
    !*                         element
    !*
    !*  OUTPUTS
    !*    The following argument has the INTENT(INOUT) attribute:
    !*
    !*    der(ndim,nod)      : Real
    !*                       : Derivatives of the shape functions at a 
    !*                         Gauss point
    !*  AUTHOR
    !*    I.M. Smith
    !*    D.V. Griffiths
    !*    L. Margetts
    !*  COPYRIGHT
    !*    (c) University of Manchester 2004-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)     :: i
    REAL(iwp),INTENT(IN)    :: points(:,:)
    REAL(iwp),INTENT(INOUT) :: der(:,:)
    REAL(iwp)               :: eta,xi,zeta,xi0,eta0,zeta0 
    REAL(iwp)               :: etam,etap,xim,xip,c1,c2,c3 
    REAL(iwp)               :: x2p1,x2m1,e2p1,e2m1,zetam,zetap 
    REAL(iwp)               :: t1,t2,t3,t4,t5,t6,t7,t8,t9 
    REAL(iwp),PARAMETER     :: one=1.0_iwp,two=2.0_iwp,d4=4.0_iwp
    REAL(iwp),PARAMETER     :: d8=8.0_iwp, zero=0.0_iwp
    INTEGER                 :: xii(20),etai(20),zetai(20),l,ndim,nod
    
    ndim = UBOUND(der , 1)
    nod  = UBOUND(der , 2)
  
    SELECT CASE (ndim)
     
      CASE(1) ! one dimensional case
      
      xi = points(i,1)
      
      SELECT CASE (nod)
     
        CASE(2)
        
        der(1,1) = -0.5
        der(1,2) = 0.5
           
        CASE(3)
        
        t1 = -1.-xi 
        t2 = -xi  
        t3 = 1.-xi
             
        der(1,1) = -(t3+t2)/2.  
        der(1,2) = (t3+t1)    
        der(1,3) = -(t2+t1)/2.   
           
        CASE(4)
        
        t1 = -1.-xi 
        t2 = -1./3.-xi 
        t3 = 1./3.-xi 
        t4 = 1.-xi
        
        der(1,1) = -(t3*t4+t2*t4+t2*t3)*9./16.     
        der(1,2) = (t3*t4+t1*t4+t1*t3)*27./16. 
        der(1,3) = -(t2*t4+t1*t4+t1*t2)*27./16. 
        der(1,4) = (t2*t3+t1*t3+t1*t2)*9./16.   
           
        CASE(5)
        
        t1 = -1.-xi 
        t2 = -0.5-xi 
        t3 = -xi 
        t4 = 0.5-xi 
        t5 = 1.-xi
             
        der(1,1) = -(t3*t4*t5+t2*t4*t5+t2*t3*t5+t2*t3*t4)*2./3.   
        der(1,2) =  (t3*t4*t5+t1*t4*t5+t1*t3*t5+t1*t3*t4)*8./3.
        der(1,3) = -(t2*t4*t5+t1*t4*t5+t1*t2*t5+t1*t2*t4)*4. 
        der(1,4) =  (t2*t3*t5+t1*t3*t5+t1*t2*t5+t1*t2*t3)*8./3.
        der(1,5) = -(t2*t3*t4+t1*t3*t4+t1*t2*t4+t1*t2*t3)*2./3.
       
        CASE DEFAULT
        
        print*,"wrong number of nodes in shape_der"        
       
      END SELECT
     
      CASE(2)      ! two dimensional elements
      
      xi   = points(i,1)
      eta  = points(i,2) 
      c1   = xi 
      c2   = eta 
      c3   = 1.-c1-c2
      etam =.25*(1.-eta)
      etap =.25*(1.+eta)
      xim  =.25*(1.-xi)
      xip  =.25*(1.+xi)
      x2p1 = 2.*xi+1. 
      x2m1 = 2.*xi-1.
      e2p1 = 2.*eta+1. 
      e2m1 = 2.*eta-1.
       
      SELECT CASE (nod)
        
        CASE(3)
       
        der(1,1) =  1.
        der(1,3) =  0.
        der(1,2) = -1.
        der(2,1) =  0.
        der(2,3) =  1.
        der(2,2) = -1.
        
        CASE(6) 
        
        der(1,1) = 4.*c1-1. 
        der(1,6) = 4.*c2 
        der(1,5) = 0.  
        der(1,4) = -4.*c2
        der(1,3) = -(4.*c3-1.)
        der(1,2) = 4.*(c3-c1)
        der(2,1) = 0.
        der(2,6) = 4.*c1 
        der(2,5) = 4.*c2-1.
        der(2,4) = 4.*(c3-c2)
        der(2,3) = -(4.*c3-1.)  
        der(2,2) = -4.*c1
        
        CASE(15)  
        
        t1=c1-.25
        t2=c1-.5 
        t3=c1-.75  
        t4=c2-.25
        t5=c2-.5   
        t6=c2-.75 
        t7=c3-.25 
        t8=c3-.5 
        t9=c3-.75
        
        der(1,1)  = 32./3.*(t2*t3*(t1+c1)+c1*t1*(t3+t2))
        der(1,12) = 128./3.*c2*(t2*(t1+c1)+c1*t1) 
        der(1,11) = 64.*c2*t4*(t1+c1)
        der(1,10) = 128./3.*c2*t4*t5  
        der(1,9)  = 0. 
        der(1,8)  = -128./3.*c2*t4*t5
        der(1,7)  = -64.*c2*t4*(t7+c3) 
        der(1,6)  = -128./3.*c2*(t8*(t7+c3)+c3*t7)
        der(1,5)  = -32./3.*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
        der(1,4)  = 128./3.*(c3*t7*t8-c1*(t8*(t7+c3)+c3*t7))
        der(1,3)  = 64.*(c3*t7*(t1+c1)-c1*t1*(t7+c3))
        der(1,2)  = 128./3.*(c3*(t2*(t1+c1)+c1*t1)-c1*t1*t2)
        der(1,13) = 128.*c2*(c3*(t1+c1)-c1*t1) 
        der(1,15) = 128.*c2*t4*(c3-c1)
        der(1,14) = 128.*c2*(c3*t7-c1*(t7+c3))
        der(2,1)  = 0.0 
        der(2,12) = 128./3.*c1*t1*t2
        der(2,11) = 64.*c1*t1*(t4+c2)
        der(2,10) = 128./3.*c1*(t5*(t4+c2)+c2*t4)
        der(2,9)  = 32./3.*(t5*t6*(t4+c2)+c2*t4*(t6+t5))
        der(2,8)  = 128./3.*((c3*(t5*(t4+c2)+c2*t4))-c2*t4*t5)
        der(2,7)  = 64.*(c3*t7*(t4+c2)-c2*t4*(t7+c3))
        der(2,6)  = 128./3.*(c3*t7*t8-c2*(t8*(t7+c3)+c3*t7))
        der(2,5)  =-32./3.*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
        der(2,4)  =-128./3.*c1*(t8*(t7+c3)+c3*t7)
        der(2,3)  =-64.*c1*t1*(t7+c3)  
        der(2,2)  =-128./3.*c1*t1*t2
        der(2,13) =128.*c1*t1*(c3-c2)
        der(2,15) =128.*c1*(c3*(t4+c2)-c2*t4)
        der(2,14) =128.*c1*(c3*t7-c2*(c3+t7))        
        
        CASE(4)  
        
        der(1,1)  = -etam
        der(1,2)  = -etap
        der(1,3)  =  etap
        der(1,4)  =  etam
        der(2,1)  =  -xim
        der(2,2)  =   xim
        der(2,3)  =   xip
        der(2,4)  =  -xip
        
        CASE(8)
        
        der(1,1)  =  etam*(2.*xi+eta)
        der(1,2)  = -8.*etam*etap
        der(1,3)  =  etap*(2.*xi-eta)
        der(1,4)  = -4.*etap*xi
        der(1,5)  =  etap*(2.*xi+eta)
        der(1,6)  =  8.*etap*etam
        der(1,7)  =  etam*(2.*xi-eta)
        der(1,8)  = -4.*etam*xi
        der(2,1)  =  xim*(xi+2.*eta)
        der(2,2)  = -4.*xim*eta
        der(2,3)  =  xim*(2.*eta-xi)
        der(2,4)  =  8.*xim*xip
        der(2,5)  =  xip*(xi+2.*eta)
        der(2,6)  = -4.*xip*eta
        der(2,7)  =  xip*(2.*eta-xi)
        der(2,8)  = -8.*xim*xip   
       
        CASE(9)
        
        etam = eta - 1.
        etap = eta + 1.
        xim  = xi - 1.
        xip  = xi + 1.
         
        der(1,1)  =  .25*x2m1*eta*etam  
        der(1,2)  = -.5*x2m1*etap*etam
        der(1,3)  =  .25*x2m1*eta*etap  
        der(1,4)  =  -xi*eta*etap
        der(1,5)  =  .25*x2p1*eta*etap  
        der(1,6)  = -.5*x2p1*etap*etam
        der(1,7)  =  .25*x2p1*eta*etam  
        der(1,8)  =  -xi*eta*etam
        der(1,9)  =   2.*xi*etap*etam  
        der(2,1)  =  .25*xi*xim*e2m1
        der(2,2)  =  -xi*xim*eta        
        der(2,3)  =  .25*xi*xim*e2p1
        der(2,4)  =  -.5*xip*xim*e2p1   
        der(2,5)  =  .25*xi*xip*e2p1
        der(2,6)  =  -xi*xip*eta        
        der(2,7)  =  .25*xi*xip*e2m1
        der(2,8)  =  -.5*xip*xim*e2m1  
        der(2,9)  =   2.*xip*xim*eta
       
        CASE DEFAULT
        
        print*,"wrong number of nodes in shape_der"        
       
      END SELECT
      
      CASE(3)  ! three dimensional elements
      
      xi     = points(i,1)
      eta    = points(i,2)
      zeta   = points(i,3)
      etam   = one  - eta 
      xim    = one  - xi
      zetam  = one  - zeta
      etap   = eta  + one 
      xip    = xi   + one 
      zetap  = zeta + one
      
      SELECT CASE (nod)
       
        CASE(4)
        
        der(1:3,1:4) =  zero
        der(1,1)     =  one
        der(2,2)     =  one
        der(3,3)     =  one
        der(1,4)     = -1.0_iwp
        der(2,4)     = -1.0_iwp 
        der(3,4)     = -1.0_iwp
        
        CASE(8)
        
        der(1,1) = -0.125_iwp*etam*zetam    
        der(1,2) = -0.125_iwp*etam*zetap
        der(1,3) =  0.125_iwp*etam*zetap     
        der(1,4) =  0.125_iwp*etam*zetam
        der(1,5) = -0.125_iwp*etap*zetam    
        der(1,6) = -0.125_iwp*etap*zetap
        der(1,7) =  0.125_iwp*etap*zetap     
        der(1,8) =  0.125_iwp*etap*zetam
        der(2,1) = -0.125_iwp*xim*zetam     
        der(2,2) = -0.125_iwp*xim*zetap
        der(2,3) = -0.125_iwp*xip*zetap     
        der(2,4) = -0.125_iwp*xip*zetam
        der(2,5) =  0.125_iwp*xim*zetam      
        der(2,6) =  0.125_iwp*xim*zetap
        der(2,7) =  0.125_iwp*xip*zetap      
        der(2,8) =  0.125_iwp*xip*zetam
        der(3,1) = -0.125_iwp*xim*etam      
        der(3,2) =  0.125_iwp*xim*etam
        der(3,3) =  0.125_iwp*xip*etam       
        der(3,4) = -0.125_iwp*xip*etam
        der(3,5) = -0.125_iwp*xim*etap     
        der(3,6) =  0.125_iwp*xim*etap
        der(3,7) =  0.125_iwp*xip*etap       
        der(3,8) = -0.125_iwp*xip*etap 

        CASE(14) ! type 6 element
       
        der(1,1)  =  (two*xi*eta+two*xi*zeta+d4*xi+eta*zeta+eta+zeta)*          &
                     (eta-one)*(zeta-one)/d8
        der(1,2)  = -(two*xi*eta-two*xi*zeta+d4*xi-eta*zeta+eta-zeta)*          &
                     (eta-one)*(zeta+one)/d8
        der(1,3)  = -(two*xi*eta-two*xi*zeta+d4*xi+eta*zeta-eta+zeta)*          &
                     (eta-one)*(zeta+one)/d8
        der(1,4)  =  (two*xi*eta+two*xi*zeta+d4*xi-eta*zeta-eta-zeta)*          &
                     (eta-one)*(zeta-one)/d8
        der(1,5)  = -(eta-one)*(zeta+one)*(zeta-one)*xi 
        der(1,6)  = -(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
        der(1,7)  =  (eta+one)*(eta-one)*(zeta+one)*xi
        der(1,8)  =  (eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
        der(1,9)  = -(eta+one)*(eta-one)*(zeta-one)*xi  
        der(1,10) =  (two*xi*eta-two*xi*zeta-d4*xi+eta*zeta+eta-zeta)*          &
                     (eta+one)*(zeta-one)/d8
        der(1,11) = -(two*xi*eta+two*xi*zeta-d4*xi-eta*zeta+eta+zeta)*          &
                     (eta+one)*(zeta+one)/d8
        der(1,12) = -(two*xi*eta+two*xi*zeta-d4*xi+eta*zeta-eta-zeta)*          &
                     (eta+one)*(zeta+one)/d8
        der(1,13) =  (two*xi*eta-two*xi*zeta-d4*xi-eta*zeta-eta+zeta)*          &
                     (eta+one)*(zeta-one)/d8
        der(1,14) =  (eta+one)*(zeta+one)*(zeta-one)*xi
        der(2,1)  =  (two*xi*eta+xi*zeta+xi+two*eta*zeta+d4*eta+zeta)*          &
                     (xi-one)*(zeta-one)/d8                                  
        der(2,2)  = -(two*xi*eta-xi*zeta+xi-two*eta*zeta+d4*eta-zeta)*          &
                     (xi-one)*(zeta+one)/d8
        der(2,3)  = -(two*xi*eta-xi*zeta+xi+two*eta*zeta-d4*eta+zeta)*          &
                     (xi+one)*(zeta+one)/d8
        der(2,4)  =  (two*xi*eta+xi*zeta+xi-two*eta*zeta-d4*eta-zeta)*          &
                     (xi+one)*(zeta-one)/d8
        der(2,5)  = -(xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
        der(2,6)  = -(xi-one)*(zeta+one)*(zeta-one)*eta
        der(2,7)  =  (xi+one)*(xi-one)*(zeta+one)*eta
        der(2,8)  =  (xi+one)*(zeta+one)*(zeta-one)*eta
        der(2,9)  = -(xi+one)*(xi-one)*(zeta-one)*eta
        der(2,10) =  (two*xi*eta-xi*zeta-xi+two*eta*zeta+d4*eta-zeta)*          &
                     (xi-one)*(zeta-one)/d8
        der(2,11) = -(two*xi*eta+xi*zeta-xi-two*eta*zeta+d4*eta+zeta)*          &
                     (xi-one)*(zeta+one)/d8
        der(2,12) = -(two*xi*eta+xi*zeta-xi+two*eta*zeta-d4*eta-zeta)*          &
                     (xi+one)*(zeta+one)/d8
        der(2,13) =  (two*xi*eta-xi*zeta-xi-two*eta*zeta-d4*eta+zeta)*          &
                     (xi+one)*(zeta-one)/d8
        der(2,14) =  (xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
        der(3,1)  =  (xi*eta+two*xi*zeta+xi+two*eta*zeta+eta+d4*zeta)*          &
                     (xi-one)*(eta-one)/d8
        der(3,2)  = -(xi*eta-two*xi*zeta+xi-two*eta*zeta+eta-d4*zeta)*          &
                     (xi-one)*(eta-one)/d8
        der(3,3)  = -(xi*eta-two*xi*zeta+xi+two*eta*zeta-eta+d4*zeta)*          &
                     (xi+one)*(eta-one)/d8
        der(3,4)  =  (xi*eta+two*xi*zeta+xi-two*eta*zeta-eta-d4*zeta)*          &
                     (xi+one)*(eta-one)/d8
        der(3,5)  = -(xi+one)*(xi-one)*(eta-one)*zeta  
        der(3,6)  = -(xi-one)*(eta+one)*(eta-one)*zeta  
        der(3,7)  =  (xi+one)*(xi-one)*(eta+one)*(eta-one)/two
        der(3,8)  =  (xi+one)*(eta+one)*(eta-one)*zeta
        der(3,9)  = -(xi+one)*(xi-one)*(eta+one)*(eta-one)/two
        der(3,10) =  (xi*eta-two*xi*zeta-xi+two*eta*zeta+eta-d4*zeta)*          &
                     (xi-one)*(eta+one)/d8
        der(3,11) = -(xi*eta+two*xi*zeta-xi-two*eta*zeta+eta+d4*zeta)*          &
                     (xi-one)*(eta+one)/d8
        der(3,12) = -(xi*eta+two*xi*zeta-xi+two*eta*zeta-eta-d4*zeta)*          &
                     (xi+one)*(eta+one)/d8
        der(3,13) =  (xi*eta-two*xi*zeta-xi-two*eta*zeta-eta+d4*zeta)*          &
                     (xi+one)*(eta+one)/d8
        der(3,14) =  (xi+one)*(xi-one)*(eta+one)*zeta
      
        CASE(20)
        
        xii   = (/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
        etai  = (/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
        zetai = (/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
        
        DO l = 1,20
        
           xi0   = xi*xii(l)
           eta0  = eta*etai(l)
           zeta0 = zeta*zetai(l)
           
           IF(l==4.or.l==8.or.l==16.or.l==20) THEN
  
            der(1,l) = -.5*xi*(1.+eta0)*(1.+zeta0)
            der(2,l) =  .25*etai(l)*(1.-xi*xi)*(1.+zeta0)
            der(3,l) =  .25*zetai(l)*(1.-xi*xi)*(1.+eta0)
  
           ELSE IF(l>=9.and.l<=12)THEN
            der(1,l) =  .25*xii(l)*(1.-eta*eta)*(1.+zeta0)
            der(2,l) = -.5*eta*(1.+xi0)*(1.+zeta0)
            der(3,l) = .25*zetai(l)*(1.+xi0)*(1.-eta*eta)
           ELSE IF(l==2.or.l==6.or.l==14.or.l==18) THEN
            der(1,l) =  .25*xii(l)*(1.+eta0)*(1.-zeta*zeta)
            der(2,l) =  .25*etai(l)*(1.+xi0)*(1.-zeta*zeta)
            der(3,l) = -.5*zeta*(1.+xi0)*(1.+eta0)
           ELSE
            der(1,l) =  .125*xii(l)*(1.+eta0)*(1.+zeta0)*(2.*xi0+eta0+zeta0-1.)
            der(2,l) =  .125*etai(l)*(1.+xi0)*(1.+zeta0)*(xi0+2.*eta0+zeta0-1.)
            der(3,l) =  .125*zetai(l)*(1.+xi0)*(1.+eta0)*(xi0+eta0+2.*zeta0-1.)
           END IF
        END DO 
        
        CASE DEFAULT
        
        print*,"wrong number of nodes in shape_der"        
     
        END SELECT
    
      CASE DEFAULT
      
      print*,"wrong number of dimensions in shape_der"
    
    END SELECT
   
    RETURN
    
  END SUBROUTINE SHAPE_DER
 
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE BEEMAT(deriv,bee)

    !/****f* elements/beemat
    !*  NAME
    !*    SUBROUTINE: beemat
    !*  SYNOPSIS
    !*    Usage:      CALL beemat(deriv,bee)
    !*  FUNCTION
    !*    Build B-array for 3D elasticity or elastoplasticity (ih=3 or 4
    !*    respectively) or for 3D (ih=6)
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    deriv(ndim,nod)     : Real
    !*                        : Shape function derivatives in the undeformed
    !*                          mesh at a Gauss point
    !*
    !*    The following argument has the INTENT(INOUT) attribute:
    !*
    !*    bee(nst,ndof)       : Real
    !*                        : Shape function derivatives in the undeformed
    !*                          element arranged in B array
    !*
    !*  AUTHOR
    !*    I.M. Smith
    !*    D.V. Griffiths
    !*    L. Margetts
    !*  COPYRIGHT
    !*    (c) University of Manchester 2004-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  Need to replace with the subroutine in Smith and Griffiths Edition 4.
    !*/

    IMPLICIT NONE

    REAL(iwp), INTENT(IN)    :: deriv(:,:)
    REAL(iwp), INTENT(INOUT) :: bee(:,:)
    INTEGER                  :: k,l,m,n,ih,nod
    REAL(iwp)                :: x,y,z

    bee = 0.0_iwp
    nod = UBOUND(deriv,2)
    ih  = UBOUND(bee,1)

    SELECT CASE (ih)
      CASE(3,4)
        DO m = 1,nod
          k = 2*m
          l = k-1
          x = deriv(1,m)
          y = deriv(2,m)
          bee(1,l) = x
          bee(3,k) = x
          bee(2,k) = y
          bee(3,l) = y
        END DO
      CASE(6)
        DO m = 1,nod
          n = 3*m
          k = n-1
          l = k-1
          x = deriv(1,m)
          y = deriv(2,m)
          z = deriv(3,m)
          bee(1,l) = x
          bee(4,k) = x
          bee(6,n) = x
          bee(2,k) = y
          bee(4,l) = y
          bee(5,n) = y
          bee(3,n) = z
          bee(5,k) = z
          bee(6,l) = z
        END DO
      CASE DEFAULT
        PRINT*,'wrong dimension for nst in bee matrix'
    END SELECT

    RETURN

  END SUBROUTINE BEEMAT

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE SAMPLE(element,s,wt)

    !/****f* elements/sample
    !*  NAME
    !*    SUBROUTINE: sample
    !*  SYNOPSIS
    !*    Usage:      CALL sample(element,s,wt)
    !*  FUNCTION
    !*    Returns the local coordinates of the integrating points. 
    !*    Source: "Programming the Finite Element Method"
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    element             : Character
    !*                        : Element type
    !*                          
    !*	OUTPUTS
    !*    The following arguments have the INTENT(INOUT) attribute:
    !*
    !*    s(:,:)              : Real
    !*                        : Local coordinates of the integrating points
    !*
    !*    wt(:)               : Real
    !*	                      : Weighting function
    !*
    !*  AUTHOR
    !*    I.M. Smith
    !*    D.V. Griffiths
    !*    L. Margetts
    !*  COPYRIGHT
    !*    (c) University of Manchester 2004-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/
    
    IMPLICIT NONE
    
    REAL(iwp),INTENT(INOUT) :: s(:,:),wt(:)
    CHARACTER(*),INTENT(IN) :: element
    INTEGER                 :: nip 
    REAL(iwp)               :: root3, r15 , w(3),v(9),b,c
    
    root3 = 1.0_iwp/sqrt(3.0_iwp)
    r15   = 0.2_iwp*sqrt(15.0_iwp) 
    nip   = UBOUND(s,1) 

    w     = (/5.0_iwp/9.0_iwp,8.0_iwp/9.0_iwp,5.0_iwp/9.0_iwp/)
    v     = (/5.0_iwp/9.0_iwp*w,8.0_iwp/9.0_iwp*w,5.0_iwp/9.0_iwp*w/)
        
    SELECT CASE (element)
    
    CASE('line')
      
      SELECT CASE(nip)
      
        CASE(1)
        
        s(1,1) = 0.
        wt(1)  = 2.
        
        CASE(2)
             
        s(1,1) = root3
        s(2,1) = -s(1,1)
        wt(1)  = 1.
        wt(2)  = 1.
        
        CASE(3)
            
        s(1,1) = r15
        s(2,1) = .0     
        s(3,1) = -s(1,1)
        wt     = w
        
        CASE(4)
        
        s(1,1) = .861136311594053  
        s(2,1) = .339981043584856  
        s(3,1) = -s(2,1)  
        s(4,1) = -s(1,1)
        wt(1)  = .347854845137454 
        wt(2)  = .652145154862546 
        wt(3)  = wt(2)
        wt(4)  = wt(1)
                
        CASE(5)
            
        s(1,1) = .906179845938664
        s(2,1) = .538469310105683  
        s(3,1) = .0 
        s(4,1) = -s(2,1)
        s(5,1) = -s(1,1)
        wt(1)  = .236926885056189
        wt(2)  = .478628670499366
        wt(3)  = .568888888888889 
        wt(4)  = wt(2)
        wt(5)  = wt(1)
                
        CASE(6)
   
        s(1,1) = .932469514203152
        s(2,1) = .661209386466265 
        s(3,1) = .238619186083197
        s(4,1) = -s(3,1)
        s(5,1) = -s(2,1)
        s(6,1) = -s(1,1)
        wt(1)  = .171324492379170
        wt(2)  = .360761573048139 
        wt(3)  = .467913934572691
        wt(4)  = wt(3)
        wt(5)  = wt(2)
        wt(6)  = wt(1)
                 
        CASE DEFAULT
                    
        print*,"wrong number of integrating points for a line"
        
      END SELECT
      
      CASE('triangle') 
        
        SELECT CASE(nip)
                
        CASE(1)   ! for triangles weights multiplied by .5
        
        s(1,1) = 1./3.  
        s(1,2) = 1./3.  
        wt(1)  = .5
        
        CASE(3)
        
        s(1,1) = .5
        s(1,2) = .5
        s(2,1) = .5  
        s(2,2) = 0.
        s(3,1) = 0.  
        s(3,2) = .5
        
        wt(1)  = 1./3.  
        wt(2)  = wt(1) 
        wt(3)  = wt(1)   
        wt     = .5*wt
        
        CASE(6)
        
        s(1,1) = .816847572980459  
        s(1,2) = .091576213509771
        s(2,1) = s(1,2)
        s(2,2) = s(1,1) 
        s(3,1) = s(1,2)
        s(3,2) = s(1,2)
        s(4,1) = .108103018168070 
        s(4,2) = .445948490915965
        s(5,1) = s(4,2) 
        s(5,2) = s(4,1) 
        s(6,1) = s(4,2)  
        s(6,2) = s(4,2)
        
        wt(1)  = .109951743655322 
        wt(2)  = wt(1)  
        wt(3)  = wt(1)
        wt(4)  = .223381589678011 
        wt(5)  = wt(4)  
        wt(6)  = wt(4)    
        wt     = .5*wt
        
        CASE(7)
        
        s(1,1) = 1./3. 
        s(1,2) = 1./3.
        s(2,1) = .797426985353087 
        s(2,2) = .101286507323456
        s(3,1) = s(2,2) 
        s(3,2) = s(2,1) 
        s(4,1) = s(2,2) 
        s(4,2) = s(2,2)
        s(5,1) = .470142064105115 
        s(5,2) = .059715871789770
        s(6,1) = s(5,2) 
        s(6,2) = s(5,1)
        s(7,1) = s(5,1)
        s(7,2) = s(5,1)
   
        wt(1)  = .225 
        wt(2)  = .125939180544827 
        wt(3)  = wt(2)
        wt(4)  = wt(2)
        wt(5)  = .132394152788506
        wt(6)  = wt(5)      
        wt(7)  = wt(5)     
        wt     = .5*wt
        
        CASE(12)
  
        s(1,1)  = .873821971016996
        s(1,2)  = .063089014491502
        s(2,1)  = s(1,2) 
        s(2,2)  = s(1,1)
        s(3,1)  = s(1,2) 
        s(3,2)  = s(1,2)
        s(4,1)  = .501426509658179 
        s(4,2)  = .249286745170910
        s(5,1)  = s(4,2) 
        s(5,2)  = s(4,1)  
        s(6,1)  = s(4,2) 
        s(6,2)  = s(4,2)
        s(7,1)  = .636502499121399
        s(7,2)  = .310352451033785
        s(8,1)  = s(7,1) 
        s(8,2)  = .053145049844816 
        s(9,1)  = s(7,2) 
        s(9,2)  = s(7,1)
        s(10,1) = s(7,2)
        s(10,2) = s(8,2) 
        s(11,1) = s(8,2)
        s(11,2) = s(7,1)
        s(12,1) = s(8,2)
        s(12,2) = s(7,2)
        wt(1)   = .050844906370207
        wt(2)   = wt(1)
        wt(3)   = wt(1)
        wt(4)   =.116786275726379 
        wt(5)   = wt(4)
        wt(6)   = wt(4)
        wt(7)   =.082851075618374 
        wt(8:12)= wt(7)           
        wt      = .5*wt
        
        CASE(16)
        
        s(1,1)  = 1./3. 
        s(1,2)  = 1./3. 
        s(2,1)  =.658861384496478
        s(2,2)  =.170569307751761 
        s(3,1)  = s(2,2) 
        s(3,2)  = s(2,1)
        s(4,1)  = s(2,2)  
        s(4,2)  = s(2,2)
        s(5,1)  =.898905543365938 
        s(5,2)  =.050547228317031
        s(6,1)  = s(5,2)
        s(6,2)  = s(5,1) 
        s(7,1)  = s(5,2)  
        s(7,2)  = s(5,2)
        s(8,1)  =.081414823414554
        s(8,2)  =.459292588292723
        s(9,1)  = s(8,2)  
        s(9,2)  = s(8,1)
        s(10,1) = s(8,2)
        s(10,2) = s(8,2)
        s(11,1) =.008394777409958
        s(11,2) =.263112829634638
        s(12,1) = s(11,1)    
        s(12,2) =.728492392955404
        s(13,1) = s(11,2) 
        s(13,2) = s(11,1)  
        s(14,1) = s(11,2)
        s(14,2) = s(12,2)
        s(15,1) = s(12,2) 
        s(15,2) = s(11,1) 
        s(16,1) = s(12,2) 
        s(16,2) = s(11,2)
        
        wt(1)     = .144315607677787 
        wt(2)     = .103217370534718 
        wt(3)     = wt(2)
        wt(4)     = wt(2)
        wt(5)     = .032458497623198 
        wt(6)     = wt(5)  
        wt(7)     = wt(5)
        wt(8)     = .095091634267284 
        wt(9)     = wt(8)   
        wt(10)    = wt(8)
        wt(11)    = .027230314174435 
        wt(12:16) = wt(11)
        wt        = .5*wt
                
        CASE DEFAULT
                    
        print*,"wrong number of integrating points for a triangle"
        
      END SELECT
        
      CASE ('quadrilateral')
          
        SELECT CASE (nip)
        
        CASE(1)
            
        s(1,1) = .0 
        wt(1)  = 4.
              
        CASE(4)
            
        s(1,1) = -root3
        s(1,2) =  root3
        s(2,1) =  root3
        s(2,2) =  root3
        s(3,1) = -root3
        s(3,2) = -root3
        s(4,1) =  root3
        s(4,2) = -root3
                  
        wt = 1.0
           
        CASE(9)
        
        s(1:7:3,1) = -r15
        s(2:8:3,1) = .0
        s(3:9:3,1) =  r15
        s(1:3,2)   =  r15
        s(4:6,2)   =  .0 
        s(7:9,2)   = -r15
                     
        wt         = v
        
        CASE DEFAULT 
        
        print*,"wrong number of integrating points for a quadrilateral"
        
      END SELECT
      
      CASE('tetrahedron')    
        
        SELECT CASE(nip)
        
        CASE(1)          ! for tetrahedra weights multiplied by 1/6
        
        s(1,1)  = 0.25_iwp 
        s(1,2)  = 0.25_iwp 
        s(1,3)  = 0.25_iwp 
        
        wt(1)   = 1.0_iwp/6.0_iwp
                
        CASE(4)
        
        s(1,1)  = .58541020
        s(1,2)  = .13819660
        s(1,3)  = s(1,2)
        s(2,2)  = s(1,1)
        s(2,3)  = s(1,2)
        s(2,1)  = s(1,2)
        s(3,3)  = s(1,1)
        s(3,1)  = s(1,2)
        s(3,2)  = s(1,2)
        s(4,1)  = s(1,2)
        s(4,2)  = s(1,2)  
        s(4,3)  = s(1,2)

        wt(1:4) = .25/6.
                
        CASE(5)
            
        s(1,1)  = .25
        s(1,2)  = .25
        s(1,3)  = .25
        s(2,1)  = .5
        s(2,2)  = 1./6.
        s(2,3)  = s(2,2)
        s(3,2)  = .5
        s(3,3)  = 1./6.
        s(3,1)  = s(3,3)   
        s(4,3)  = .5
        s(4,1)  = 1./6. 
        s(4,2)  = s(4,1)
        s(5,1)  = 1./6.
        s(5,2)  = s(5,1)
        s(5,3)  = s(5,1) 

        wt(1)   = -.8  
        wt(2)   = 9./20.
        wt(3:5) = wt(2)  
        wt      = wt/6.  
                
        CASE(6)
        
        wt     =  4./3.        
        s(6,3) =  1.
        s(1,1) = -1.
        s(2,1) =  1.
        s(3,2) = -1. 
        s(4,2) =  1.
        s(5,3) = -1.
        
        CASE DEFAULT
        
        print*,"wrong number of integrating points for a tetrahedron"
        
      END SELECT
      
      CASE('hexahedron')
        
        SELECT CASE ( nip )
        
        CASE(1)
        
        s(1,1:3) = 0.0_iwp
        wt(1)    = 8.0_iwp
        
        CASE(8)   
            
        s(1,1) =  root3
        s(1,2) =  root3
        s(1,3) =  root3
        s(2,1) =  root3
        s(2,2) =  root3
        s(2,3) = -root3
        s(3,1) =  root3
        s(3,2) = -root3
        s(3,3) =  root3
        s(4,1) =  root3
        s(4,2) = -root3
        s(4,3) = -root3
        s(5,1) = -root3
        s(5,2) =  root3
        s(5,3) =  root3
        s(6,1) = -root3
        s(6,2) = -root3
        s(6,3) =  root3
        s(7,1) = -root3
        s(7,2) =  root3
        s(7,3) = -root3
        s(8,1) = -root3
        s(8,2) = -root3
        s(8,3) = -root3
       
        wt     = 1.0_iwp                                                
              
        CASE(14)
            
        b  = 0.795822426
        c  = 0.758786911
            
        wt(1:6) = 0.886426593
        wt(7:)  = 0.335180055
        s(1,1)  = -b 
        s(2,1)  =  b  
        s(3,2)  = -b 
        s(4,2)  =  b
        s(5,3)  = -b 
        s(6,3)  =  b
        s(7:,:) =  c
        s(7,1)  = -c 
        s(7,2)  = -c  
        s(7,3)  = -c 
        s(8,2)  = -c 
        s(8,3)  = -c
        s(9,1)  = -c 
        s(9,3)  = -c 
        s(10,3) = -c
        s(11,1) = -c
        s(11,2) = -c 
        s(12,2) = -c 
        s(13,1) = -c
        
        CASE(15)
        
        b  = 1.     
        c  = 0.674199862
            
        wt(1)    = 1.564444444 
        wt(2:7)  = 0.355555556  
        wt(8:15) = 0.537777778
            
        s(2,1)   = -b
        s(3,1)   =  b 
        s(4,2)   = -b
        s(5,2)   =  b
        s(6,3)   = -b  
        s(7,3)   =  b
        s(8:,:)  =  c
        s(8,1)   = -c
        s(8,2)   = -c 
        s(8,3)   = -c
        s(9,2)   = -c 
        s(9,3)   = -c
        s(10,1)  = -c
        s(10,3)  = -c 
        s(11,3)  = -c 
        s(12,1)  = -c
        s(12,2)  = -c
        s(13,2)  = -c  
        s(14,1)  = -c                          
                
        CASE(27)
            
        wt           = (/5./9.*v,8./9.*v,5./9.*v/)
        s(1:7:3,1)   = -r15
        s(2:8:3,1)   =   .0
        s(3:9:3,1)   =  r15
        s(1:3,3)     =  r15
        s(4:6,3)     =   .0 
        s(7:9,3)     = -r15
        s(1:9,2)     = -r15
        s(10:16:3,1) = -r15
        s(11:17:3,1) =   .0
        s(12:18:3,1) =  r15
        s(10:12,3)   =  r15
        s(13:15,3)   =   .0 
        s(16:18,3)   = -r15
        s(10:18,2)   =   .0   
        s(19:25:3,1) = -r15
        s(20:26:3,1) =   .0
        s(21:27:3,1) =  r15
        s(19:21,3)   =  r15
        s(22:24,3)   =   .0
        s(25:27,3)   = -r15
        s(19:27,2)   =  r15
         
        CASE DEFAULT
        
        print*,"wrong number of integrating points for a hexahedron" 
        
      END SELECT
      
      CASE DEFAULT
      
      print*,"not a valid element type" 
      
    END SELECT
    
    RETURN
    
  END SUBROUTINE SAMPLE

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE ECMAT(ecm,fun,ndof,nodof)

    !/****f* elements/ecmat
    !*  NAME
    !*    SUBROUTINE: ecmat
    !*  SYNOPSIS
    !*    Usage:      CALL ecmat(ecm,fun,ndof,nodof)
    !*  FUNCTION
    !*    Returns the consistent mass matrix ECM for an element with shape
    !*    functions FUN, NDOF freedoms and NODOF freedoms per node 
    !*    Source: "Programming the Finite Element Method"
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    fun(:)              : Real
    !*                        : Element shape functions
    !*
    !*    nodof               : Integer
    !*                        : Freedoms per node
    !*
    !*    ndof                : Integer
    !*                        : Number of degrees of freedom per element
    !*                          
    !*	OUTPUTS
    !*    The following arguments have the INTENT(INOUT) attribute:
    !*
    !*    ecm(:,:)            : Real
    !*                        : Element consistent mass matrix
    !*
    !*  AUTHOR
    !*    I.M. Smith
    !*    D.V. Griffiths
    !*  COPYRIGHT
    !*    (c) University of Manchester 2004-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    IMPLICIT NONE
    
    REAL(iwp),INTENT(IN)    :: fun(:)
    REAL(iwp),INTENT(INOUT) :: ecm(:,:)
    REAL(iwp)               :: nt(ndof,nodof),tn(nodof,ndof)
    INTEGER,INTENT(IN)      :: nodof,ndof
    INTEGER                 :: nod,i,j
    
    nod = ndof/nodof
    nt  = .0
    tn  = .0
    
    DO i = 1 , nod 
      DO j = 1 , nodof
        nt((i-1)*nodof+j,j) = fun(i)
        tn(j,(i-1)*nodof+j) = fun(i)
      END DO
    END DO
    
    ecm = matmul(nt,tn)
    
    RETURN
    
  END SUBROUTINE ECMAT
  
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE DEEMAT(e,v,dee)

    !/****f* elements/deemat
    !*  NAME
    !*    SUBROUTINE: deemat
    !*  SYNOPSIS
    !*    Usage:      CALL deemat(e,v,dee)
    !*  FUNCTION
    !*    Compute the material matrix (D-matrix) for a linear elastic material.
    !*    It's valid for plane strain (ih=3), axisymmetry or plane strain
    !*    elastoplasticity (ih=4) or 3D problems (ih=6)
    !*     
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    e                         : Real
    !*                              : Young's modulus
    !*
    !*    v                         : Real
    !*                              : Poisson's ratio
    !*
    !*    The following argument has the INTENT(INOUT) attribute:
    !*
    !*    dee(nst,nst)              : Real
    !*                              : Material matrix for linear elasticity
    !*  AUTHOR
    !*    Smith and Griffiths 4th Edition
    !*  COPYRIGHT
    !*    (c) University of Manchester 2004-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    IMPLICIT NONE
    INTEGER,PARAMETER       :: iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)    :: e,v
    REAL(iwp),INTENT(INOUT) :: dee(:,:)
    REAL(iwp)               :: v1,v2,c,vv
    REAL(iwp),PARAMETER     :: zero=0.0_iwp,pt5=0.5_iwp
    REAL(iwp),PARAMETER     :: one=1.0_iwp,two=2.0_iwp
    INTEGER                 :: i,ih
    
    dee = zero  
    ih  = UBOUND(dee,1)
    v1  = one-v
    c   = e/((one+v)*(one-two*v))

    SELECT CASE(ih)

    CASE(3)
      dee(1,1) = v1*c
      dee(2,2) = v1*c
      dee(1,2) = v*c
      dee(2,1) = v*c
      dee(3,3) = pt5*c*(one-two*v)
    CASE(4)
      dee(1,1) = v1*c
      dee(2,2) = v1*c
      dee(4,4) = v1*c
      dee(3,3) = pt5*c*(one-two*v) 
      dee(1,2) = v*c
      dee(2,1) = v*c
      dee(1,4) = v*c
      dee(4,1) = v*c
      dee(2,4) = v*c
      dee(4,2) = v*c
    CASE(6)
      v2       = v/(one-v)
      vv       = (one-two*v)/(one-v)*pt5
      DO i=1,3
        dee(i,i)=one
      END DO
      DO i=4,6
        dee(i,i)=vv
      END DO
      dee(1,2) = v2
      dee(2,1) = v2
      dee(1,3) = v2
      dee(3,1) = v2
      dee(2,3) = v2
      dee(3,2) = v2
      dee      = dee*e/(two*(one+v)*vv)
    CASE DEFAULT
      WRITE(*,*)'wrong size for dee matrix'
    END SELECT
 
  END SUBROUTINE DEEMAT

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  
  SUBROUTINE SHAPE_DERIVATIVES(i,points,der)

    !/****f* elements/shape_derivatives
    !*  NAME
    !*    SUBROUTINE: shape_derivatives
    !*  SYNOPSIS
    !*    Usage:      CALL shape_derivatives(i,points,der)
    !*  FUNCTION
    !*    Compute the derivatives of the shape functions at a Gauss point
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    i                  : Integer
    !*                       : Gauss point number
    !*
    !*    points(nip,ndim)   : Real
    !*                       : Gauss points coordinates at the reference 
    !*                         element
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    der(ndim,nod)      : Real
    !*                       : Derivatives of the shape functions at a 
    !*                         Gauss point
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    24.01.2008
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  In this routine, tetrahedra elements are supposed to be oriented
    !*  with the 4 node right-handed from the first 3 nodes.
    !*   
    !*  This routine has to be completed with elements for 1 and 2 dimensions  
    !*  and 3D elements with unusual number of nodes
    !*
    !*  A programming branch by Fran that needs to be removed from ParaFEM
    !*/
     
    INTEGER,   INTENT(IN)  :: i
    REAL(iwp), INTENT(IN)  :: points(:,:)
    REAL(iwp), INTENT(OUT) :: der(:,:) 
    INTEGER :: ndim, nod
    REAL(iwp) :: xi, xip, xim, eta, etap, etam, zeta, zetap, zetam 

    ndim = UBOUND(der,1)
    nod  = UBOUND(der,2)
    
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
          CASE(4) ! 4-node tetrahedra
            der(1,1) =  1._iwp  !N1_x
            der(2,1) =  0._iwp  !N1_y
            der(3,1) =  0._iwp  !N1_z
            der(1,2) =  0._iwp  !N2_x
            der(2,2) =  1._iwp  !N2_y
            der(3,2) =  0._iwp  !etc...
            der(1,3) = -1._iwp
            der(2,3) = -1._iwp
            der(3,3) = -1._iwp
            der(1,4) =  0._iwp
            der(2,4) =  0._iwp
            der(3,4) =  1._iwp
          CASE(8) !8-node hexahedra
            der(1,1) =  0.125_iwp*etam*zetam !N1_x
            der(2,1) = -0.125_iwp*xip*zetam  !N1_y
            der(3,1) = -0.125_iwp*xip*etam   !N1_z
            der(1,2) =  0.125_iwp*etap*zetam !N2_x
            der(2,2) =  0.125_iwp*xip*zetam  !N2_y
            der(3,2) = -0.125_iwp*xip*etap   !etc...
            der(1,3) = -0.125_iwp*etap*zetam
            der(2,3) =  0.125_iwp*xim*zetam
            der(3,3) = -0.125_iwp*xim*etap
            der(1,4) = -0.125_iwp*etam*zetam
            der(2,4) = -0.125_iwp*xim*zetam
            der(3,4) = -0.125_iwp*xim*etam
            der(1,5) =  0.125_iwp*etam*zetap
            der(2,5) = -0.125_iwp*xip*zetap
            der(3,5) =  0.125_iwp*xip*etam
            der(1,6) =  0.125_iwp*etap*zetap
            der(2,6) =  0.125_iwp*xip*zetap
            der(3,6) =  0.125_iwp*xip*etap
            der(1,7) = -0.125_iwp*etap*zetap
            der(2,7) =  0.125_iwp*xim*zetap
            der(3,7) =  0.125_iwp*xim*etap
            der(1,8) = -0.125_iwp*etam*zetap
            der(2,8) = -0.125_iwp*xim*zetap
            der(3,8) =  0.125_iwp*xim*etam
          CASE DEFAULT
            WRITE(*,*)"Wrong number of nodes!!"
        END SELECT
      CASE DEFAULT
        WRITE(*,*)"Wrong number of dimensions!!"
    END SELECT

    RETURN

  END SUBROUTINE SHAPE_DERIVATIVES
  
END MODULE ELEMENTS
