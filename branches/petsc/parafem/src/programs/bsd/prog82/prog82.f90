PROGRAM EBE_BEM
!------------------------------------------------------------------------------
!  General purpose BEM program for solving elasticity problems 
!  This version uses EBE with BiCGSTAB(l) in serial
!------------------------------------------------------------------------------
USE bem_lib        ;    IMPLICIT NONE  ! N_dof replaces Ndof
INTEGER, ALLOCATABLE :: Inci(:,:)  !  Element Incidences
INTEGER, ALLOCATABLE :: BCode(:,:), NCode(:) !  Element BC´s
INTEGER, ALLOCATABLE :: Ldest(:,:) !  Element destination vector
INTEGER, ALLOCATABLE :: Ndest(:,:) !  Node destination vector
REAL(iwp), ALLOCATABLE :: Elres_u(:,:)  !  Element results , u
REAL(iwp), ALLOCATABLE :: Elres_t(:,:)  !  Element results , t
REAL(iwp), ALLOCATABLE :: Elcor(:,:)    !  Element coordinates
REAL(iwp), ALLOCATABLE :: xP(:,:)       !  Node co-ordinates
REAL(iwp), ALLOCATABLE :: dUe(:,:),dTe(:,:),lhs(:,:),Diag(:,:),pmul(:)
REAL(iwp), ALLOCATABLE :: km(:,:),qmul(:)
REAL(iwp), ALLOCATABLE :: store_dUe(:,:,:),store_dTe(:,:,:)
REAL(iwp), ALLOCATABLE :: F(:)     !  global RHS
REAL(iwp), ALLOCATABLE :: u1(:)    !  global vector of unknowns 
CHARACTER (LEN=80) :: Title
INTEGER :: Cdim,Node,m,n,Istat,Nodel,Nel,N_dof,Toa
INTEGER :: Nreg,Ltyp,Nodes,Maxe,Ndofe,Ndofs,ndg,NE_u,NE_t                
INTEGER :: nod,nd,i,j,k,l,DoF,Pos,Isym,nsym,nsy
INTEGER :: its,iters,ell
REAL(iwp),ALLOCATABLE    :: Fac(:)     !  Factors for symmetry
REAL(iwp),ALLOCATABLE    :: Elres_te(:),Elres_ue(:)   
INTEGER,ALLOCATABLE :: Incie(:)   !  Incidences for one element
INTEGER,ALLOCATABLE :: Ldeste(:),g(:)  !  Destination vector 1 elem
REAL(iwp) :: Con,E,ny,Scat,Scad,tol,kappa,alpha,beta,rho,gama,       &
             omega,norm_r,r0_norm,error,one=1._iwp,zero=.0_iwp
LOGICAL:: converged
REAL(iwp),ALLOCATABLE::s(:),GG(:,:),Gamma(:),                        &
                       rt(:),y(:),y1(:),r(:,:),uu(:,:)
!-----------------------------------------------------
!   Read job information
!-----------------------------------------------------
OPEN (UNIT=11,FILE='prog82.dat',FORM='FORMATTED') !  Input
OPEN (UNIT=12,FILE='prog82.res',FORM='FORMATTED')!  Output
Call Jobin(Title,Cdim,N_dof,Toa,Nreg,Ltyp,Con,E,ny,&
           Isym,nodel,nodes,maxe)
Nsym= 2**Isym   !   number of symmetry loops
ALLOCATE(xP(Cdim,Nodes))   !  Array for node coordinates
ALLOCATE(Inci(Maxe,Nodel)) !  Array for incidences
CALL Geomin(Nodes,Maxe,xp,Inci,Nodel,Cdim)
Ndofe= Nodel*N_dof   !    Total degrees of freedom of element
ALLOCATE(BCode(Maxe,Ndofe))      !    Element Boundary codes
ALLOCATE(Elres_u(Maxe,Ndofe),Elres_t(Maxe,Ndofe))       
CALL BCinput(Elres_u,Elres_t,Bcode,nodel,ndofe,n_dof)
READ(11,*) tol,its,ell,kappa ! BiCGStab  data
ALLOCATE(Ldest(maxe,Ndofe))  ! Elem. destination vector
ALLOCATE(Ndest(Nodes,N_dof))
!---------------------------------------------------------------------
!     Determine Node destination vector and Element destination vector 
!---------------------------------------------------------------------
CALL Destination(Isym,Ndest,Ldest,xP,Inci,Ndofs,nodes,N_dof,Nodel,Maxe)
!---------------------------------------------------------------------
!     Determine global Boundary code vector
!---------------------------------------------
ALLOCATE(NCode(Ndofs))            
NCode=0
DoF_o_System: &
DO  nd=1,Ndofs
    DO Nel=1,Maxe
        DO m=1,Ndofe
           IF (nd == Ldest(Nel,m) .and. NCode(nd) == 0) THEN 
               NCode(nd)= NCode(nd)+BCode(Nel,m)
           END IF
        END DO
    END DO
END DO &
DoF_o_System
IF(N_dof ==1)E= Con
CALL Scal(E,xP(:,:),Elres_u(:,:),Elres_t(:,:),Cdim,Scad,Scat)
ALLOCATE(dTe(Ndofs,Ndofe),dUe(Ndofs,Ndofe),lhs(Ndofs,Ndofe))! Elem.coef.matrices
ALLOCATE(store_dTe(Maxe,Ndofs,Ndofe),store_dUe(Maxe,Ndofs,Ndofe)) ! store els
ALLOCATE(Diag(Ndofs,N_dof))          ! Diagonal coefficients
ALLOCATE(F(Ndofs),u1(Ndofs))         ! global arrays
ALLOCATE(Elcor(Cdim,Nodel))             !  Elem. Coordinates
ALLOCATE(Fac(Ndofe))                 !  Factor for symmetric elements
ALLOCATE(Incie(Nodel))               !  Element incidences
ALLOCATE(Ldeste(Ndofe),pmul(Ndofe),km(N_dof,N_dof),g(N_dof),qmul(N_dof))!dest.
ALLOCATE(Elres_te(Ndofe),Elres_ue(Ndofe))    ! Tractions of Element
!----------------------------------------------------------------
!  Compute element coefficient matrices
!----------------------------------------------------------------
Lhs=zero; F = zero; u1 = zero; Diag = zero;store_dUe = zero; store_dTe = zero
Elements_1:&
DO Nel=1,Maxe
        Symmetry_loop:&
        DO nsy= 1,Nsym
                Elcor(:,:)= xP(:,Inci(Nel,:))  !   gather element coordinates
                Incie= Inci(nel,:)             !   incidences
                Ldeste= Ldest(nel,:)           !   and destinations
                Fac(1:nodel*n_dof)= 1.0_iwp
                Elres_te(:)=Elres_t(Nel,:)
                IF(Isym > 0) THEN
                   CALL Mirror(Isym,nsy,Nodes,Elcor,Fac,    &
                         Incie,Ldeste,Elres_te,Elres_ue,nodel,n_dof,Cdim) 
                END IF
                IF(Cdim == 2) THEN
                        IF(N_dof == 1) THEN
                          CALL Integ2P(Elcor,Incie,Nodel,Nodes,       &
                                xP,Con,dUe,dTe,Ndest,Isym)
                        ELSE
                          CALL Integ2E(Elcor,Incie,Nodel,Nodes,       &
                                xP,E,ny,dUe,dTe,Ndest,Isym)
                        END IF
                ELSE
                   CALL Integ3(Elcor,Incie,Nodel,Nodes,xP,N_dof &
                                    ,E,ny,Con,dUe,dTe,Ndest,Isym)    
                END IF
!  Now build global F and diag but not LHS
                CALL rhs_and_diag(F,DTe,DUe,Ldeste,BCode(Nel,:),Ncode &
                        ,Elres_u(Nel,:),Elres_te,Diag,Ndofe,N_dof,Nodel,Fac)
        END DO &
        Symmetry_loop
        store_dUe(Nel,:,:) = dUe;    store_dTe(Nel,:,:) = dTe
END DO &
Elements_1
!------------------------------------------------------------
!  Add azimuthal integral for infinite regions  
!------------------------------------------------------------
IF(Nreg == 2) THEN
        DO m=1, Nodes
                DO n=1, N_dof
                        IF(Ndest(m,n) == 0)CYCLE
                        k=Ndest(m,n)
                        Diag(k,n) = Diag(k,n) + 1.0_iwp
                END DO
        END DO           
END IF   
!-------------------------------------------------------------
!  Store active Diagonal coefficients     
!-------------------------------------------------------------
DO m=1,Ndofs            ! Loop over collocation points
   Nod=0
   DO n=1, Nodes
      DO l=1,N_dof
         IF (m == Ndest(n,l))THEN
              Nod=n   ;       EXIT
         END IF  
      END DO
      IF (Nod /= 0)EXIT
   END DO
   DO k=1,N_dof
      DoF=Ndest(Nod,k)
      IF(DoF /= 0) THEN
         IF(NCode(DoF) == 1) THEN
               Nel=0    ;     Pos=0
           DO i=1,Maxe
              DO j=1,Ndofe
                 IF(DoF == Ldest(i,j))THEN
                    Nel=i   ;    Pos=j   ;    EXIT
                 END IF
              END DO
              IF(Nel /= 0)EXIT
           END DO
           F(m) = F(m) - Diag(m,k) * Elres_u(Nel,Pos)
         END IF
       END IF
   END DO
END DO   
!---------------------------------------------------------
!   Solve system of equations  element by element
!---------------------------------------------------------
ALLOCATE(s(ell+1),GG(ell+1,ell+1),Gamma(ell+1),                 &
            rt(Ndofs),y(Ndofs),y1(Ndofs),r(Ndofs,ell+1),uu(Ndofs,ell+1))
!      initialisation phase
     u1 = zero   ;    y = u1    ;    y1 = zero
Elements_2 : DO Nel = 1 , Maxe
     Dte = store_DTe(Nel,:,:);  DUe = store_DUe(Nel,:,:)
     Ldeste = Ldest(Nel,:); pmul = y(Ldeste)   
     CALL form_lhs(lhs,DTe,DUe,Ldeste,BCode(Nel,:),Ncode,Ndofe,Fac)
     y1 = y1 + MATMUL(lhs,pmul) 
END DO Elements_2 
DO i = 1 , nodes  
     CALL get_km(Cdim,i,y,Diag,g,qmul,km)
     y1(g) = y1(g) + MATMUL(km,qmul)
END DO    
    y=y1; rt = F - y
    r=zero ; r(:,1) = rt   ;  uu = zero  ; gama = one  ; omega=one 
    norm_r = norm(rt)   ;   r0_norm = norm_r   ;  error = one   ;  iters = 0
!     bicgstab(ell)  iterations
iterations : DO
   iters = iters + 1    ;            converged = error  < tol  
   IF(iters==its.OR. converged) EXIT 
   gama = - omega*gama  ;  y = r(:,1)
    DO j = 1 , ell
               rho = DOT_PRODUCT(rt,y)  ;  beta = rho/gama
               uu(:,1:j) = r(:,1:j) - beta * uu(:,1:j)   ;  y = uu(:,j)
               y1 = zero
         Elements_3: DO Nel = 1 , Maxe
               Dte = store_DTe(Nel,:,:);  DUe = store_DUe(Nel,:,:)
               Ldeste = Ldest(Nel,:); pmul = y(Ldeste)
               CALL form_lhs(lhs,DTe,DUe,Ldeste,BCode(Nel,:),Ncode,Ndofe,Fac)
               y1 = y1 + MATMUL(lhs,pmul)   
         END DO Elements_3  
         DO i = 1 , nodes   
               CALL get_km(Cdim,i,y,Diag,g,qmul,km)
               y1(g) = y1(g) + MATMUL(km,qmul)
         END DO       
               y=y1; uu(:,j+1) = y
               gama = DOT_PRODUCT(rt,y); alpha = rho/gama
               u1=u1+ alpha * uu(:,1)
               r(:,1:j) = r(:,1:j) - alpha * uu(:,2:j+1)    
               y = r(:,j)
               y1 = zero
         Elements_4: DO Nel = 1 , Maxe
                Dte = store_DTe(Nel,:,:);  DUe = store_DUe(Nel,:,:)
                Ldeste = Ldest(Nel,:);  pmul = y(Ldeste)
                CALL form_lhs(lhs,DTe,DUe,Ldeste,BCode(Nel,:),Ncode,Ndofe,Fac)
                y1 = y1 + MATMUL(lhs,pmul)  
         END DO Elements_4   
         DO i = 1 , nodes  
                CALL get_km(Cdim,i,y,Diag,g,qmul,km)
                y1(g) = y1(g) + MATMUL(km,qmul)
         END DO     
               y=y1 ; r(:,j+1) = y
    END DO
    GG = MATMUL(TRANSPOSE(r),r)
    CALL form_s(gg,ell,kappa,omega,gamma,s)
    u1 = u1 - MATMUL(r,s);r(:,1)=MATMUL(r,Gamma)
    uu(:,1)=MATMUL(uu,Gamma)
    norm_r = norm(r(:,1))  ;  error = norm_r/r0_norm       
END DO iterations
  WRITE(12,'(/,A,I5,A,/)')"It took BiCGSTAB_L ",iters," iterations to converge"
!   Gather Element results from global result vector u1
Elements_5:     &
DO nel=1,maxe , maxe - 1 
   D_o_F1:         &
   DO nd=1,Ndofe
      IF(Ldest(nel,nd) /= 0)THEN
         IF(NCode(Ldest(nel,nd)) == 0) THEN
            Elres_u(nel,nd) = Elres_u(nel,nd) + u1(Ldest(nel,nd))
         ELSE
            Elres_t(nel,nd) = Elres_t(nel,nd) + u1(Ldest(nel,nd))
         END IF
      END IF
   END DO &
   D_o_F1  
   Elres_u(nel,:)= Elres_u(nel,:) * Scad
   Elres_t(nel,:)= Elres_t(nel,:) / Scat
   WRITE(12,'(24F12.5)') (Elres_u(nel,m), m=1,Ndofe)
   WRITE(12,'(24F12.5)') (Elres_t(nel,m), m=1,Ndofe)
END DO &
Elements_5
END PROGRAM EBE_BEM
