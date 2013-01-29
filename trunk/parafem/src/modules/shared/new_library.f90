MODULE NEW_LIBRARY

  !/****h* /new_library
  !*  NAME
  !*    MODULE: new_library
  !*  SYNOPSIS
  !*    Usage:      USE new_library
  !*  FUNCTION
  !*    Contains subroutines provided in the Smith and Griffiths text book
  !*    
  !*    Subroutine             Purpose
  !*    
  !*  AUTHORS
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*    L. Margetts
  !*  COPYRIGHT
  !*    2004-2013 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  Module required by the 5th Edition of "Programming the Finite Element
  !*  Method. Take care when modifying
  !*/
  
  USE precision
  USE geometry

  CONTAINS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

SUBROUTINE formm(stress,m1,m2,m3)
!
! This subroutine forms the derivatives of the invariants with respect to
! stress in 2- or 3-d.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::stress(:)
 REAL(iwp),INTENT(OUT)::m1(:,:),m2(:,:),m3(:,:)
 REAL(iwp)::sx,sy,txy,tyz,tzx,sz,dx,dy,dz,sigm,zero=0.0_iwp,one=1.0_iwp,  &
   two=2.0_iwp,three=3.0_iwp,six=6.0_iwp,nine=9.0_iwp
 INTEGER::nst,i,j
 nst=UBOUND(stress,1)
 SELECT CASE(nst)
 CASE(4)
   sx=stress(1)
   sy=stress(2)
   txy=stress(3)
   sz=stress(4)
   dx=(two*sx-sy-sz)/three
   dy=(two*sy-sz-sx)/three
   dz=(two*sz-sx-sy)/three
   sigm=(sx+sy+sz)/three
   m1=zero
   m2=zero
   m3=zero
   m1(1,1:2)=one
   m1(2,1:2)=one
   m1(4,1:2)=one
   m1(1,4)=one
   m1(4,4)=one
   m1(2,4)=one
   m1=m1/nine/sigm
   m2(1,1)=two/three
   m2(2,2)=two/three
   m2(4,4)= two/three
   m2(2,4)=-one/three
   m2(4,2)=-one/three
   m2(1,2)=-one/three
   m2(2,1)=-one/three
   m2(1,4)=-one/three
   m2(4,1)=-one/three
   m2(3,3)=two
   m3(3,3)=-dz
   m3(1:2,3)=txy/three
   m3(3,1:2)=txy/three
   m3(3,4)=-two*txy/three
   m3(4,3)=-two*txy/three
   m3(1,1)=dx/three
   m3(2,4)=dx/three
   m3(4,2)=dx/three
   m3(2,2)=dy/three
   m3(1,4)=dy/three
   m3(4,1)=dy/three
   m3(4,4)=dz/three
   m3(1,2)=dz/three
   m3(2,1)=dz/three
 CASE(6)
   sx=stress(1)
   sy=stress(2)    
   sz=stress(3)
   txy=stress(4)  
   tyz=stress(5) 
   tzx=stress(6)
   sigm=(sx+sy+sz)/three
   dx=sx-sigm  
   dy=sy-sigm 
   dz=sz-sigm
   m1=zero
   m2=zero
   m1(1:3,1:3)=one/(three*sigm)
   DO i=1,3 
     m2(i,i)=two 
     m2(i+3,i+3)=six 
   END DO
   m2(1,2)=-one
   m2(1,3)=-one 
   m2(2,3)=-one
   m3(1,1)=dx
   m3(1,2)=dz 
   m3(1,3)=dy 
   m3(1,4)=txy  
   m3(1,5)=-two*tyz
   m3(1,6)=tzx 
   m3(2,2)=dy 
   m3(2,3)=dx 
   m3(2,4)=txy
   m3(2,5)=tyz 
   m3(2,6)=-two*tzx 
   m3(3,3)=dz
   m3(3,4)=-two*txy
   m3(3,5)=tyz 
   m3(3,6)=tzx
   m3(4,4)=-three*dz 
   m3(4,5)=three*tzx
   m3(4,6)=three*tyz
   m3(5,5)=-three*dx
   m3(5,6)=three*txy 
   m3(6,6)=-three*dy
   DO i=1,6 
     DO j=i+1,6
       m1(j,i)=m1(i,j) 
       m2(j,i)=m2(i,j)   
       m3(j,i)=m3(i,j)
     END DO
   END DO
   m1=m1/three
   m2=m2/three
   m3=m3/three
 CASE DEFAULT
   WRITE(*,*)"nst size not recognised in formm"
 END SELECT
RETURN   
END SUBROUTINE formm

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

  SUBROUTINE BEEMAT(bee,deriv)

    !/****f* elements/beemat
    !*  NAME
    !*    SUBROUTINE: beemat
    !*  SYNOPSIS
    !*    Usage:      CALL beemat(bee,deriv)
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

  SUBROUTINE DEEMAT(dee,e,v)

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

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE FORMNF_CAVITY(nf,nftol,nodof,aa,bb,cc,nels,nxe,nze,problem)
  
    !/****f* structure_dof/formnf_cavity
    !*  NAME
    !*    SUBROUTINE: formnf_cavity
    !*  SYNOPSIS
    !*    Usage:      CALL formnf_cavity(nf,nftol,nodof,aa,bb,cc,nels,     &
    !*                                   nxe,nze,problem)
    !*  FUNCTION
    !*    Forms node freedom array NF for a rectangular lid driven cavity
    !*    of any size
    !*  
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    nels                : Integer
    !*                        : Total number of elements
    !* 
    !*    nxe                 : Integer
    !*                        : Number of elements in the x direction
    !* 
    !*    nye                 : Integer
    !*                        : Number of elements in the y direction
    !* 
    !*    problem             : Integer
    !*                        : Code for the problem type
    !*                          1 = full cavity
    !*                          2 = symmetric half cavity
    !*                          Any other value = undefined
    !*
    !*    rest(nr,nodof+1)    : Integer
    !*                        : List of the nodes with some degree of freedom
    !*                          fixed. In the input, the degrees of freedom
    !*                          fixed have 1 and the degrees of freedom not
    !*                          fixed have 0. In the output, the degrees of
    !*                          freedom fixed have 0 and the degrees of freedom
    !*                          not fixed have the global equation number
    !*
    !*    The following arguments have the INTENT(INOUT) attribute:    
    !*
    !*    nf(nn,nodof+1)      : Integer
    !*                        : Node freedom array
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    05.03.2008
    !*  COPYRIGHT
    !*    (c) University of Manchester
    !******
    !*  Need to modify with an error trap for problem.
    !*  Also need to deallocate allocated arrays.
    !*/


    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: nf(:,:)
    INTEGER, INTENT(IN)    :: nels,nxe,nze,problem
    INTEGER                :: i,j,n,nye,nn,nodof,iel
    INTEGER, ALLOCATABLE   :: num(:)
    REAL(iwp), INTENT(IN)  :: aa,bb,cc,nftol
    REAL(iwp)              :: minx,maxx,miny,maxy,minz,maxz
    REAL(iwp), ALLOCATABLE :: coord(:,:),g_coord(:,:)

    !-----------------------------------------------------------------------
    ! 1. Allocate arrays 
    !-----------------------------------------------------------------------
   
    nn = UBOUND(nf,2)
    
    ALLOCATE(coord(20,3),num(20),g_coord(3,nn))

    !-----------------------------------------------------------------------
    ! 2. Initiallize 
    !-----------------------------------------------------------------------
   
    minx = 0._iwp
    miny = 0._iwp
    minz = 0._iwp
    nye  = (nels/nxe)/nze
    maxx = aa * nxe
    maxy = bb * nye
    maxz = cc * nze
    
    DO i = 1,nn
      nf(1,i)   = 0
      nf(2:4,i) = 1
    END DO

    !-----------------------------------------------------------------------
    ! 3. Find coordinates of nodes and mark equations
    !-----------------------------------------------------------------------
     
    DO iel = 1,nels
        CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,num)
        DO i = 1,7,2
          nf(1,num(i)) = 1
        END DO
        DO i = 13,19,2
          nf(1,num(i)) = 1
        END DO
        g_coord(:,num) = transpose(coord)
    END DO
    
    
    node_freedoms: SELECT CASE(problem)
    
    !-----------------------------------------------------------------------
    ! 4. Find node freedoms for a full lid-driven cavity
    !-----------------------------------------------------------------------
    
    CASE (1)
    
    DO i = 1,nn
      IF(ABS(g_coord(1,i)-minx)<nftol.OR.abs(g_coord(1,i)-maxx)<nftol) THEN
        nf(2:4,i) = 0
      END IF
      IF(ABS(g_coord(2,i)-miny)<nftol.OR.abs(g_coord(2,i)-maxy)<nftol) THEN
        nf(2:4,i) = 0
      END IF
      IF(ABS(g_coord(3,i)-minz)<nftol) THEN
        nf(3:4,i) = 0
        nf(2,i)   = 1
      END IF
      IF(ABS(g_coord(3,i)+maxz)<nftol) THEN
        nf(2:4,i) = 0
        IF (ABS(g_coord(2,i)-miny)<nftol) THEN
          IF (ABS(g_coord(1,i)-minx)<nftol) nf(1,i) = 0
        END IF
      END IF
    END DO
    
    !-----------------------------------------------------------------------
    ! 5. Find node freedoms for a symmetric half lid-driven cavity
    !-----------------------------------------------------------------------
    
    CASE (2)
    
    DO i = 1,nn
      IF(ABS(g_coord(1,i)-minx)<nftol.OR.abs(g_coord(1,i)-maxx)<nftol) THEN
        nf(2:4,i) = 0
      END IF
      IF(ABS(g_coord(2,i)-miny)<nftol) THEN
        nf(2:4,i) = 0
      END IF
      IF(ABS(g_coord(2,i)-maxy)<nftol) nf(3,i) = 0
      IF(abs(g_coord(3,i)-minz)<nftol) THEN
        nf(2,i) = 1
        nf(3,i) = 0
        nf(4,i) = 0
      END IF
      IF(ABS(g_coord(3,i)+maxz)<nftol) THEN
        nf(2:4,i) = 0
        IF (ABS(g_coord(2,i)-miny)<nftol) THEN
          IF (ABS(g_coord(1,i)-minx)<nftol) nf(1,i) = 0
        END IF
      END IF
    END DO
    
    END SELECT node_freedoms

    !-----------------------------------------------------------------------
    ! 6. Number the equations
    !-----------------------------------------------------------------------
    
    n = 0
    DO j = 1,nn
      DO i = 1,nodof
        IF(nf(i,j)/=0) THEN
          n       = n+1
          nf(i,j) = n
        END IF
      END DO
    END DO  
   
  RETURN

  END SUBROUTINE FORMNF_CAVITY

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE FIXITY_DATA(g_coord,node,sense,val,lid_velocity,nftol,nn,problem)
  
    !/****f* structure_dof/fixity_data
    !*  NAME
    !*    SUBROUTINE: fixity_data
    !*  SYNOPSIS
    !*    Usage:      CALL fixity_data(g_coord,node,sense,val,lid_velocity,   &
    !*                                 nftol,nn,problem)
    !*  FUNCTION
    !*    Fixes the value of velocity on all nodes of the lid
    !*  
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    problem             : Integer
    !*                        : Code for the problem type.
    !*                          1 = full cavity
    !*                          2 = symmetric half cavity
    !*                          Any other value = undefined    
    !* 
    !*    g_coord(:,:)        : Real
    !*                        : Co-ordinates of the nodes
    !* 
    !*    lid_velocity        : Real
    !*                        : The prescribed velocity of the lid
    !* 
    !*    nftol               : Real
    !*                        : The tolerance used when searching for nodes
    !*                          on the basis of their position
    !*
    !*  OUTPUTS
    !*    The following arguments have the INTENT(OUT) attribute:
    !* 
    !*    node(:)             : Integer
    !*                        : A vector listing the nodes positioned on the
    !*                          lid that will have an imposed velocity
    !* 
    !*    sense(:)            : Integer
    !*                        : A vector identifying which of the x, y or z 
    !*                          freedoms will have an imposed velocity
    !*                          1 = x direction
    !*                          2 = y direction
    !*                          3 = z direction
    !*                          Any other value is currently undefined
    !* 
    !*    val(:)              : Real
    !*                        : A vector that holds the velocity values
    !* 
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    19.04.2002
    !*  COPYRIGHT
    !*    (c) University of Manchester
    !******
    !*  Need to modify with an error trap for problem.
    !*  Need to have error trap for sense too.
    !*/

    IMPLICIT NONE
    
    INTEGER,  INTENT(OUT) :: node(:),sense(:)
    INTEGER,  INTENT(IN)  :: problem
    INTEGER               :: i,l,nn
    REAL(iwp),INTENT(IN)  :: g_coord(:,:),lid_velocity,nftol
    REAL(iwp),INTENT(OUT) :: val(:)
    REAL(iwp)             :: minz,miny

    !-----------------------------------------------------------------------
    ! 1. Initiallize
    !-----------------------------------------------------------------------

    nn    = UBOUND(g_coord,2)
    l     = 0
    node  = 0
    sense = 0
    val   = 0._iwp
    miny  = 0._iwp
    minz  = 0._iwp

    lid_vel: SELECT CASE (problem)
    
    !-----------------------------------------------------------------------
    ! 2. Fix the fluid velocities for a full lid-driven cavity
    !-----------------------------------------------------------------------
 
    CASE (1)
  
    DO i = 1,nn
      IF(ABS(g_coord(3,i)-minz)<nftol) THEN
        l        = l+1
        node(l)  = i
        sense(l) = 2
        val(l)   = lid_velocity
      END IF
    END DO
    
    !-----------------------------------------------------------------------
    ! 3. Fix the fluid velocities for a half lid-driven cavity
    !-----------------------------------------------------------------------

    CASE (2)
  
    DO i=1,nn
      IF(ABS(g_coord(3,i)-minz)<nftol) THEN
        l=l+1
        node(l)=i
        sense(l)=2
        val(l)=lid_velocity
      END IF
    END DO

    END SELECT lid_vel

    RETURN

  END SUBROUTINE fixity_data
    
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE FMKEPARP(keparp,c12,c32,c42)

    !/****f* fluid_flow/fmkeparp
    !*  NAME
    !*    SUBROUTINE: fmkeparp
    !*  SYNOPSIS
    !*    Usage:      CALL fmkeparp(keparp,c12,c32,c42)
    !*  FUNCTION
    !*    Forms sub-parent 'stiffness' matrix related to the pressure component
    !*    in 3-d for u-p-v-w version of Navier Stokes. KE incomplete without 
    !*    KEPARV and C11. c.f. FMKEPARV and ST_C11 (stores C11 by element)
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:    
    !*
    !*    c12(:,:)            : Real
    !*                        : The c12 submatrix from the 3D Navier-Stokes
    !*                          element stiffness matrix
    !*
    !*    c32(:,:)            : Real
    !*                        : The c32 submatrix from the 3D Navier-Stokes
    !*                          element stiffness matrix
    !*
    !*    c42(:,:)            : Real
    !*                        : The c42 submatrix from the 3D Navier-Stokes
    !*                          element stiffness matrix
    !*    
    !*  OUTPUTS
    !*    The following arguments have the INTENT(OUT) attribute:    
    !*
    !*    keparp(:,:)         : Real
    !*                        : The sub-parent 'stiffness' matrix
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    24.04.2002
    !*  COPYRIGHT
    !*    (c) University of Manchester
    !******
    !*/    

    IMPLICIT NONE
    
    REAL(iwp),INTENT(IN)  :: c12(:,:),c32(:,:),c42(:,:)
    REAL(iwp),INTENT(OUT) :: keparp(:,:)
    INTEGER               :: nod,nodf
    
    nodf = 8
    nod  = 20
    
    keparp(1:nod,1:nodf)                 = c12 
    keparp(nod+1:nod+nod,1:nodf)         = c32
    keparp(nod+nod+1:nod+nod+nod,1:nodf) = c42
    
    RETURN
    
  END SUBROUTINE FMKEPARP

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
 

  SUBROUTINE FMKEPARV(keparv,c21,c23,c24)

    !/****f* fluid_flow/fmkeparv
    !*  NAME
    !*    SUBROUTINE: fmkeparv
    !*  SYNOPSIS
    !*    Usage:      CALL fmkeparv(keparv,c21,c23,c24)
    !*  FUNCTION
    !*    Forms sub-parent 'stiffness' matrix related to the velocity component
    !*    in 3-d for u-p-v-w version of Navier Stokes. KE incomplete without 
    !*    KEPARP and C11. c.f. FMKEPARP and ST_C11 (stores C11 by element)
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:    
    !*
    !*    c21(:,:)            : Real
    !*                        : The c21 submatrix from the 3D Navier-Stokes
    !*                          element stiffness matrix
    !*
    !*    c23(:,:)            : Real
    !*                        : The c23 submatrix from the 3D Navier-Stokes
    !*                          element stiffness matrix
    !*
    !*    c24(:,:)            : Real
    !*                        : The c24 submatrix from the 3D Navier-Stokes
    !*                          element stiffness matrix
    !*    
    !*  OUTPUTS
    !*    The following arguments have the INTENT(OUT) attribute:    
    !*
    !*    keparv(:,:)         : Real
    !*                        : The sub-parent 'stiffness' matrix
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    24.04.2002
    !*  COPYRIGHT
    !*    (c) University of Manchester
    !******
    !*/    

     IMPLICIT NONE
     
     REAL(iwp),INTENT(IN)  :: c21(:,:),c23(:,:),c24(:,:)
     REAL(iwp),INTENT(OUT) :: keparv(:,:)
     INTEGER               :: nod,nodf
     
     nodf =  8
     nod  = 20
     
     keparv(1:nodf,1:nod)                 = c21 
     keparv(1:nodf,nod+1:nod+nod)         = c23
     keparv(1:nodf,nod+nod+1:nod+nod+nod) = c24
  
     RETURN
     
   END SUBROUTINE FMKEPARV
 
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine formnf_upvw(nf,nftol,nodof,aa,bb,cc,nels,nxe,nze,problem)
!   forms nf array for lid-driven cavity problem of any size
!   rectangular cavity only upvw order
 implicit none
 integer, intent(inout)::nf(:,:)
 integer, intent(in)::nels,nxe,nze,problem
 integer::i,j,n,nye,nn,nodof,iel
 integer, allocatable::num(:)
 real(iwp), intent(in)::aa,bb,cc,nftol
 real(iwp)::minx,maxx,miny,maxy,minz,maxz
 real(iwp), allocatable :: coord(:,:),g_coord(:,:)
 nn=ubound(nf,2)
 allocate (coord(20,3),num(20),g_coord(3,nn))
 minx=0.; miny=0.; minz=0.; nye=(nels/nxe)/nze
 maxx=aa*nxe; maxy=bb*nye; maxz=cc*nze
 do i=1,nn
     nf(1:4,i)=1;nf(2,i)=0
 end do
 do iel=1,nels
     call geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,num)
     do i=1,7,2; nf(2,num(i))=1; end do
     do i=13,19,2; nf(2,num(i))=1; end do
     g_coord(:,num)=transpose(coord)
 end do
 node_freedoms: select case(problem)
 !---- lid-driven cavity ------------------------------------------------------
 case (1)
  do i=1,nn
   if(abs(g_coord(1,i)-minx)<nftol.or.abs(g_coord(1,i)-maxx)<nftol) then
        nf(1,i)=0; nf(3,i)=0; nf(4,i)=0
   end if
   if(abs(g_coord(2,i)-miny)<nftol.or.abs(g_coord(2,i)-maxy)<nftol) then
     nf(1,i)=0; nf(3,i)=0; nf(4,i)=0
   end if
   if(abs(g_coord(3,i)-minz)<nftol) then
     nf(1,i)=1; nf(3,i)=0; nf(4,i)=0
   end if
   if(abs(g_coord(3,i)+maxz)<nftol) then
     nf(1,i)=0; nf(3,i)=0; nf(4,i)=0
     if (abs(g_coord(2,i)-miny)<nftol) then
     if (abs(g_coord(1,i)-minx)<nftol) nf(2,i)=0
     end if
   end if
 end do
 !---- lid-driven 'half' cavity -----------------------------------------------
 case (2)
  do i=1,nn
   if(abs(g_coord(1,i)-minx)<nftol.or.abs(g_coord(1,i)-maxx)<nftol) then
        nf(1,i)=0; nf(3,i)=0; nf(4,i)=0
   end if
   if(abs(g_coord(2,i)-miny)<nftol) then
     nf(1,i)=0; nf(3,i)=0; nf(4,i)=0
   end if
   if(abs(g_coord(2,i)-maxy)<nftol) nf(3,i)=0
   if(abs(g_coord(3,i)-minz)<nftol) then
     nf(1,i)=1; nf(3,i)=0; nf(4,i)=0
   end if
   if(abs(g_coord(3,i)+maxz)<nftol) then
     nf(1,i)=0; nf(3,i)=0; nf(4,i)=0
     if (abs(g_coord(2,i)-miny)<nftol) then
     if (abs(g_coord(1,i)-minx)<nftol) nf(2,i)=0
     end if
   end if
 end do
  end select node_freedoms
 n=0
   do j=1,nn; do i=1,nodof
     if(nf(i,j)/=0) then
     n=n+1;nf(i,j)=n
   end if
   end do; end do  
return
end subroutine formnf_upvw

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine fixity_upvw(g_coord,node,sense,val,lid_velocity,nftol,nn,problem)
!fixes value of lid velocity on all nodes of lid - upvw order
implicit none
integer,intent(out)::node(:),sense(:)
integer::i,l,nn,problem
real(iwp),intent(in)::g_coord(:,:),lid_velocity,nftol
real(iwp),intent(out)::val(:)
real(iwp)::minz,miny
nn=ubound(g_coord,2)
l=0;node=0;sense=0;val=0.;miny=0.;minz=.0
lid_vel: select case(problem)
 !---- lid-driven cavity ------------------------------------------------------
 case (1)
  do i=1,nn
    if(abs(g_coord(3,i)-minz)<nftol) then
       l=l+1; node(l)=i; sense(l)=1; val(l)=lid_velocity
    end if
  end do
!---- lid-driven 'half' cavity -----------------------------------------------
 case (2)
  do i=1,nn
    if(abs(g_coord(2,i)-miny)<nftol) then
       l=l+1; node(l)=i; sense(l)=1; val(l)=lid_velocity
    end if
  end do
! do l=1,nn; write(11,'(2i5,e12.4)') node(l),sense(l),val(l); end do
end select lid_vel
return
end subroutine fixity_upvw

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

! Modified From Smith and Griffiths:

subroutine formupvw(storke,iel,c11,c12,c21,c23,c32,c24,c42)
!   forms unsymmetrical 'stiffness' matrix in 3-d for
!   u-p-v-w version of navier stokes
 implicit none
 real(iwp),intent(in)::c11(:,:),c12(:,:),c21(:,:),c23(:,:),c32(:,:),     &
                  c24(:,:),c42(:,:)
 real(iwp),intent(out):: storke(:,:,:); integer:: nod,nodf,ntot,iel
 nod=ubound(c11,1); nodf=ubound(c21,1); ntot=nod+nodf+nod+nod
 storke(1:nod,1:nod,iel) = c11; storke(1:nod,nod+1:nod+nodf,iel) = c12 
 storke(nod+1:nod+nodf,1:nod,iel) = c21
 storke(nod+1:nod+nodf,nod+nodf+1:nod+nodf+nod,iel) = c23
 storke(nod+1:nod+nodf,nod+nodf+nod+1:ntot,iel) = c24 
 storke(nod+nodf+1:nod+nodf+nod,nod+1:nod+nodf,iel) = c32
 storke(nod+nodf+1:nod+nodf+nod,nod+nodf+1:nod+nodf+nod,iel) = c11 
 storke(nod+nodf+nod+1:ntot,nod+1:nod+nodf,iel) = c42
 storke(nod+nodf+nod+1:ntot,nod+nodf+nod+1:ntot,iel) = c11     
end subroutine formupvw
 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE REARRANGE(rest)

    !/****f* steering/rearrange
    !*  NAME
    !*    SUBROUTINE: rearrange
    !*  SYNOPSIS
    !*    Usage:      CALL rearrange(rest)
    !*  FUNCTION
    !*    Modifies REST array, converting restrained flags into restrained
    !*    equation numbers.   
    !*  INPUTS
    !*    The following argument has the INTENT(INOUT) attribute:
    !*
    !*    rest                : Integer
    !*                        : Array of restrained nodes
    !*                        : Possible input values are 1 or 0
    !*                        : 1 - unrestrained
    !*                        : 0 - restrained
    !*  AUTHOR
    !*    I.M Smith
    !*  COPYRIGHT
    !*    (c) University of Manchester 2004-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/
    
    IMPLICIT NONE
    
    INTEGER, INTENT(INOUT) :: rest(:,:)
    INTEGER                :: nr, nodof, i, k, m

    nr    = ubound(rest,1)
    nodof = ubound(rest,2) - 1

    m = 0

    DO i = 1, nr
      DO k = 2, nodof+1
        IF (rest(i,k)/=0) THEN
          m = m + 1
          rest(i,k) = m
        END IF
      END DO
    END DO

    DO i = 1, nr
      k = rest(i,1)
      DO m = 2, nodof+1
        IF (rest(i,m)/=0) rest(i,m) = rest(i,m) + nodof*(k-i)
      END DO
    END DO

  END SUBROUTINE rearrange

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  subroutine rearrange_2(rest)
  implicit none; integer,intent(inout)::rest(:,:)
  integer::nr,i,m; nr=ubound(rest,1)
    DO i=1,nr ; rest(i,2) = rest(i,1) - i; END DO
    DO i=1,nr;  IF(rest(i,2)/=0) THEN ; m=i-1; EXIT ; END IF; END DO
    DO i=m,nr;  rest(i,2) = rest(i,2) + 1; END DO
  end subroutine rearrange_2

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE find_g3(num,g,rest)

    !/****f* steering/find_g3
    !*  NAME
    !*    SUBROUTINE: find_g3
    !*  SYNOPSIS
    !*    Usage:      CALL find_g3(num,g,rest)
    !*  FUNCTION
    !*    Finds g from node numbers and restraints "rest"
    !*  INPUTS
    !*    The following argument has the INTENT(INOUT) attribute:
    !*
    !*    rest                : Integer
    !*                        : Array of restrained nodes
    !*                        : Possible input values are 1 or 0
    !*                        : 1 - unrestrained
    !*                        : 0 - restrained
    !*  AUTHOR
    !*    I.M Smith
    !*  COPYRIGHT
    !*    (c) University of Manchester 2004-2010
    !******
    !*  Seems to have a bug 
    !*  Smith and Griffiths "Programming the Finite Element Method", Edition 4
    !*
    !*/

    INTEGER, INTENT(IN)  :: rest(:,:),num(:)
    INTEGER, INTENT(OUT) :: g(:)
    INTEGER              :: i,j,l,nod,nodof,s1,s2,s3,only,first,last,half
    LOGICAL              :: found
    
    nod   = ubound(num,1) 
    nodof = ubound(rest,2) - 1

    DO i=1,nod 

      l     = num(i)
      first = 1 
      last  = ubound(rest,1)

      DO
        IF(first==last) EXIT ! silly array or converged
        half = (first + last)/2
        IF(l<=rest(half,1)) THEN
          last  = half      ! discard second half
        ELSE
          first = half + 1  ! discard first half
        END IF
      END DO

      only  = first
      found = (l==rest(only,1))

      IF(found) THEN
        DO j = 1 , nodof
          g(nodof*i-nodof+j) = rest(only,j+1)  
        END DO
      ELSE
        IF(only==ubound(rest,1)) THEN
          k = 0
        ELSE
          k = 1
        END IF 
        DO
          s1 = only - k
          IF(sum(rest(s1,2:))/=0.OR.s1==1) EXIT
          k = k + 1
        END DO
        s2 = maxval(rest(s1,2:))    
        IF(only==ubound(rest,1)) THEN
          s3 = l - rest(s1,1) - (k)
        ELSE
          s3 = l - rest(s1,1) - (k - 1)
        END IF
        DO j=1,nodof
          g(nodof*i-nodof+j) = s2+(nodof*s3)-(nodof-j)
        END DO
      END IF

    END DO
    
  END SUBROUTINE find_g3

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
  subroutine find_g(num,g,rest)
  ! Finds g from node numbers and restraints "rest"
  ! Only for uncoupled problems
  implicit none
  integer,intent(in)::num(:),rest(:,:); integer,intent(out)::g(:)
  integer::inc,i,j,k,l,count,nr,nod,nodof
  nr=ubound(rest,1); nod=ubound(num,1);nodof=ubound(rest,2)-1
  g=1;inc=0
  do i=1,nod
   count=0;l=num(i)
   do j=1,nr
      if(l>rest(j,1)) then
         do k=2,nodof+1; if(rest(j,k)==0)count=count+1;end do
      end if
   end do
   do k=2,nodof+1
    inc=inc+1
    do j=1,nr
     if(l==rest(j,1).and.rest(j,k)==0) then
       g(inc)=0; count=count+1
     end if
    end do
    if(g(inc)/=0) g(inc) = l*nodof-count-(nodof+1-k)
   end do
  end do
  end subroutine find_g
 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  subroutine find_g4(num,g,rest)
 ! 4th try - binary search  for one degree of freedom per node
  integer,intent(in)::rest(:,:),num(:); integer,intent(out)::g(:)
  integer:: i,l,nod,only,first,last,half ,s1,s2,s3
  logical::found       ; nod = ubound(num,1) ; nodof = ubound(rest,2) - 1
  do i=1,nod ; l = num(i)
   first = 1 ;  last = ubound(rest,1)
    do
     if(first==last) exit ! silly array
     half = (first + last)/2
     if(l<=rest(half,1)) then
       last = half  ! discard second half
     else
       first = half + 1  ! discard first half
     end if
    end do
   only = first; found = (l==rest(only,1))
   IF(found) THEN ;  g(i) = 0
   ELSE ; s1 = only -1; s2 = rest(s1,2); s3 = l - rest(s1,1)
    g(i) = s2 + s3 - 1
   END IF
  end do
  end subroutine find_g4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE g_t_g_ns(nod,g_t,g)
  
    !/****f* structure_dof/g_t_g_ns
    !*  NAME
    !*    SUBROUTINE: g_t_g_ns
    !*  SYNOPSIS
    !*    Usage:      CALL g_t_g_ns(nod,g_t,g)
    !*  FUNCTION
    !*    Finds g from g_t (Navier - Stokes)
    !*  INPUTS
    !*  AUTHOR
    !*    Ian M. Smith
    !*  COPYRIGHT
    !*    (c) University of Manchester
    !******
    !*  
    !*/

    IMPLICIT NONE
    
    INTEGER,INTENT(IN)  :: nod,g_t(:)
    INTEGER,INTENT(OUT) :: g(:)
    INTEGER             :: i
    
    DO i=1,nod
      g(i)    = g_t(4*i-3)
      g(i+28) = g_t(4*i-1)
      g(i+48) = g_t(4*i) 
    END DO
    
    DO i=1,4 
      g(i+20) = g_t(8*i-6)
      g(i+24) = g_t(8*i+42)
    END DO
    
  END SUBROUTINE g_t_g_ns

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------  

  SUBROUTINE FIND_NO_PP_TEMP(G_NUM_PP,G_G_PP,IEQ_START,NDIM,NEQ_PP,NODE,      &
                             SENSE,NREQ_PP,NREQ_PP_START,NO_PP_TEMP)

    !/****f* steering/find_no_pp_temp
    !*  NAME
    !*    SUBROUTINE: find_no_pp_temp
    !*  SYNOPSIS
    !*    Usage:      CALL find_no_pp_temp(g_num_pp,g_g_pp,node,sense,        &
    !*                                     no_pp_temp,nreq_pp)
    !*  FUNCTION
    !*    Forms vector of loaded equations NO_PP_TEMP from vector of loaded 
    !*    nodes NODE using equation numbers stored in G_G_PP. The subroutine
    !*    also returns the number of freedoms that have applied loads or fixed
    !*    displacements on the calling processor. Will work for serial and 
    !*    parallel programs.
    !*  INPUTS
    !*    Node numbers in NODE can be in any order.
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    11.11.2011
    !*  COPYRIGHT
    !*    (c) University of Manchester 2011-2012
    !******
    !*  A much better strategy would be to know which elements the fixed nodes
    !*  belong to. This could be provided in the input deck. The equation 
    !*  numbers could then be found more quickly. At the moment this subroutine
    !*  searches blindly to determine which element holds the node and then 
    !*  picks the equation number from g_g_pp.
    !*/

    IMPLICIT NONE
    INTEGER,INTENT(IN)    :: g_num_pp(:,:)
    INTEGER,INTENT(IN)    :: g_g_pp(:,:)
    INTEGER,INTENT(IN)    :: sense(:)
    INTEGER,INTENT(IN)    :: node(:)
    INTEGER,INTENT(IN)    :: ndim           
    INTEGER,INTENT(IN)    :: ieq_start
    INTEGER,INTENT(IN)    :: neq_pp
    INTEGER,INTENT(INOUT) :: no_pp_temp(:)
    INTEGER,INTENT(OUT)   :: nreq_pp        ! restrained freedoms per processor
    INTEGER,INTENT(OUT)   :: nreq_pp_start  
    INTEGER               :: i,iel,j        ! loop counters
    INTEGER               :: nels_pp        ! number of elements per processor
    INTEGER               :: nod            ! number of nodes per element
    INTEGER               :: nr             ! number of restrained nodes
    INTEGER               :: pos            ! position in g_g_pp array
    INTEGER,ALLOCATABLE   :: found(:)
 
!------------------------------------------------------------------------------
! 1. Initiallize variables
!------------------------------------------------------------------------------
 
    nels_pp       = UBOUND(g_num_pp,2)
    nod           = UBOUND(g_num_pp,1)
    nr            = UBOUND(node,1)
    nreq_pp       = 0
    nreq_pp_start = 0
    no_pp_temp    = 0
    pos           = 0

    ALLOCATE(found(nr))
    found         = 0

!------------------------------------------------------------------------------
! 2. Loop over restrained nodes to find restrained equations
!------------------------------------------------------------------------------    

    outer:   DO i = 1,nr
    middle:    DO iel = 1,nels_pp
    inner:       DO j = 1,nod
                   IF(node(i) == g_num_pp(j,iel)) THEN
                     found(i)            = 1
                     nreq_pp             = nreq_pp + 1
                     pos                 = ((j-1)*ndim) + sense(i)
                     no_pp_temp(nreq_pp) = g_g_pp(pos,iel)
                     EXIT middle ! found so look for next equation
                   END IF
                 END DO inner
               END DO middle
             END DO outer

!------------------------------------------------------------------------------
! 3. Delete equations that are not present on the local processor
!------------------------------------------------------------------------------

   DO i = 1,nr
     IF(no_pp_temp(i) < ieq_start) THEN
       found(i) = 0
     ELSE IF (no_pp_temp(i) >= ieq_start + neq_pp) THEN
       found(i) = 0
     END IF
   END DO

!------------------------------------------------------------------------------
! 4. Find the starting fixed equation on the processor
!------------------------------------------------------------------------------

    find: DO i = 1, nr
            IF(found(i) == 0) THEN
              nreq_pp_start = nreq_pp_start + 1
            ELSE IF(found(i) == 1) THEN
              nreq_pp_start = nreq_pp_start + 1
              EXIT find
            END IF
          END DO find
      
    RETURN
  
  END SUBROUTINE FIND_NO_PP_TEMP
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------ 

  SUBROUTINE REINDEX (ieq_start,no,no_local_temp,num_no,          &
                                    no_index_start,neq_pp)
  
    !/****f* steering/reindex
    !*  NAME
    !*    SUBROUTINE: reindex
    !*  SYNOPSIS
    !*    Usage:      CALL reindex(ieq_start,no,no_local_temp,             &
    !*                                         num_no,no_index_start,neq_pp)
    !*  FUNCTION
    !*    Creates a local index of loaded equations. Will function in both
    !*    serial and MPI based programs.
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    ieq_start           : Integer
    !*                        : The starting equation number on the processor
    !*
    !*    neq_pp              : Integer
    !*                          Number of equations per processor
    !*
    !*    no(:)               : Integer
    !*                        : List of the nodes with some degree of freedom
    !*                          loaded.
    !*
    !*    The following arguments have the INTENT(INOUT) attribute:    
    !*
    !*    no_local_temp(:)    : Integer
    !*                        : A temporary array
    !*
    !*  OUTPUTS
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    num_no              : Integer
    !*                        : The number of local fixed nodes
    !*
    !*    no_index_start      : Integer
    !*                        : The starting position in the array no_local_temp
    !*    
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    24.04.2002
    !*  COPYRIGHT
    !*    (c) University of Manchester 2002-2010
    !******
    !*  Note: This looks more like REINDEX_FIXED_EQUATIONS. Rename?
    !*
    !*/
    
    IMPLICIT NONE
  
    INTEGER, INTENT(IN)     :: ieq_start, neq_pp, no(:)
    INTEGER, INTENT(INOUT)  :: no_local_temp(:)
    INTEGER, INTENT(OUT)    :: num_no, no_index_start
    INTEGER                 :: fixed_nodes, i
  
    fixed_nodes    = UBOUND(no,1)
    no_index_start = 0
    num_no         = 0
  
    DO i = 1,fixed_nodes
      IF (no(i) < ieq_start) THEN
        CYCLE
      ELSE IF (no(i) >= ieq_start + neq_pp) THEN
        EXIT
      ELSE IF (no_index_start==0) THEN
        no_index_start   = i
        no_local_temp(1) = no(i)
        num_no           = 1
      ELSE
        num_no                = num_no + 1
        no_local_temp(num_no) = no(i)
      END IF
    END DO

  END SUBROUTINE REINDEX
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------  
  
  SUBROUTINE ABAQUS2SG(element,g_num)

    !/****f* steering/abaqus2sg
    !*  NAME
    !*    SUBROUTINE: abaqus2sg
    !*  SYNOPSIS
    !*    Usage:      CALL abaqus2sg(element,g_num)
    !*  FUNCTION
    !*    Converts the node steering array from the ABAQUS node ordering to 
    !*    the ordering in Smith and Griffiths "Programming the Finite Element
    !*    Method", 4th Edition. Works with both serial and parallel programs.
    !*
    !*    For 20-node hexahedra, the conversion is:
    !*
    !*    S&G   :   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
    !*    Abaqus:   1  7 19 13  3  5 17 15  8 12 20  9  4 11 16 10  2  6 18 14 
    !*  INPUTS
    !*
    !*  AUTHOR
    !*    L. Margetts
    !*  COPYRIGHT
    !*    (c) University of Manchester 2009-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/
  
    IMPLICIT NONE
    INTEGER, INTENT(INOUT)        :: g_num(:,:)
    INTEGER, ALLOCATABLE          :: temp(:)
    INTEGER                       :: iel
    INTEGER                       :: nod
    INTEGER                       :: nels
    CHARACTER(LEN=15), INTENT(IN) :: element

!------------------------------------------------------------------------------
! 1. Find dimensions of problem and allocate temp
!------------------------------------------------------------------------------

    nels = UBOUND(g_num,2)
    nod  = UBOUND(g_num,1)
    
!------------------------------------------------------------------------------
! 2. Identify which type of element is being dealt with and reorder the nodes
!------------------------------------------------------------------------------

    SELECT CASE (element)
    
    CASE('hexahedron')
    
      SELECT CASE (nod)
      
      CASE(20)
   
      ALLOCATE(temp(nod))
      temp = 0

      DO iel = 1,nels
        temp          = g_num(:,iel)
        g_num(1,iel)  = temp(4)
        g_num(2,iel)  = temp(12)
        g_num(3,iel)  = temp(1)
        g_num(4,iel)  = temp(9)
        g_num(5,iel)  = temp(2)
        g_num(6,iel)  = temp(10)
        g_num(7,iel)  = temp(3)
        g_num(8,iel)  = temp(11)
        g_num(9,iel)  = temp(20)
        g_num(10,iel) = temp(17)
        g_num(11,iel) = temp(18)
        g_num(12,iel) = temp(19)
        g_num(13,iel) = temp(8)
        g_num(14,iel) = temp(16)
        g_num(15,iel) = temp(5)
        g_num(16,iel) = temp(13)
        g_num(17,iel) = temp(6)
        g_num(18,iel) = temp(14)
        g_num(19,iel) = temp(7)
        g_num(20,iel) = temp(15)
      END DO
        
      DEALLOCATE(temp)
     
      CASE(8)
      
      ALLOCATE(temp(nod))
      temp = 0

      DO iel = 1,nels
        temp          = g_num(:,iel)
        g_num(1,iel)  = temp(1)
        g_num(2,iel)  = temp(5)
        g_num(3,iel)  = temp(6)
        g_num(4,iel)  = temp(2)
        g_num(5,iel)  = temp(4)
        g_num(6,iel)  = temp(8)
        g_num(7,iel)  = temp(7)
        g_num(8,iel)  = temp(3)
      END DO
     
      DEALLOCATE(temp)
      
      CASE default
      
        PRINT*
        PRINT*, "This element type not supported in ABAQUS2SG"
        PRINT*, "Program aborting"
        PRINT*

        STOP
 
      END SELECT
    
    CASE('tetrahedron')
      
      SELECT CASE (nod)
      
      CASE(4)

      ALLOCATE(temp(nod))
      
      PRINT*
      PRINT*, "This is untested for stress elements"
      
      DO iel = 1,nels
        temp          = g_num(:,iel)
!        g_num(1,iel)  = temp(1)
!        g_num(2,iel)  = temp(3)
!        g_num(3,iel)  = temp(4)
!        g_num(4,iel)  = temp(2)
        g_num(1,iel)  = temp(1)
        g_num(2,iel)  = temp(3)
        g_num(3,iel)  = temp(2)
        g_num(4,iel)  = temp(4)
      END DO

      DEALLOCATE(temp)

      CASE default
        
        PRINT*
        PRINT*, "Wrong number of nodes for a tetrahedron"
        PRINT*, "Program aborting"
        PRINT*
        
        STOP

      END SELECT
    
    CASE default
    
      PRINT*
      PRINT*, "Wrong type of element in subroutine ABAQUS2SG"
      PRINT*, "Program aborting"
      PRINT*
     
      STOP
 
    END SELECT
    
    RETURN
    
    END SUBROUTINE ABAQUS2SG
    
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

END MODULE NEW_LIBRARY
