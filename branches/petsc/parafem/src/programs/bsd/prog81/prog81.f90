PROGRAM General_purpose_BEM
!------------------------------------------------------------------------------
!  General purpose BEM program for solving elasticity and potential problems
!  This version uses iterative equation solution by BiCGStab(l) in serial
!------------------------------------------------------------------------------
USE bem_lib      !  contains precision
USE precision 
IMPLICIT NONE      !  Ndof changed to N_dof
INTEGER, ALLOCATABLE :: Inci(:,:)  !  Element Incidences
INTEGER, ALLOCATABLE :: BCode(:,:), NCode(:) !  Element BC´s
INTEGER, ALLOCATABLE :: Ldest(:,:) !  Element destination vector
INTEGER, ALLOCATABLE :: Ndest(:,:) !  Node destination vector
REAL(iwp), ALLOCATABLE :: Elres_u(:,:)  !  Element results , u
REAL(iwp), ALLOCATABLE :: Elres_t(:,:)  !  Element results , t
REAL(iwp), ALLOCATABLE :: Elcor(:,:)    !  Element coordinates
REAL(iwp), ALLOCATABLE :: xP(:,:)       !  Node co-ordinates
REAL(iwp), ALLOCATABLE :: dUe(:,:),dTe(:,:),Diag(:,:)
REAL(iwp), ALLOCATABLE :: Lhs(:,:),F(:)
REAL(iwp), ALLOCATABLE :: u1(:)    !  global vector of unknowns 
CHARACTER (LEN=80) :: Title
INTEGER :: Cdim,Node,m,n,Istat,Nodel,Nel,N_dof,Toa
INTEGER :: Nreg,Ltyp,Nodes,Maxe,Ndofe,Ndofs,ndg,NE_u,NE_t               
INTEGER :: nod,nd,i,j,k,l,DoF,Pos,Isym,nsym,nsy,its,ell
REAL(iwp),ALLOCATABLE    :: Fac(:)     !  Factors for symmetry
REAL(iwp),ALLOCATABLE    :: Elres_te(:),Elres_ue(:)   
INTEGER,ALLOCATABLE :: Incie(:)   !  Incidences for one element
INTEGER,ALLOCATABLE :: Ldeste(:)  !  Destination vector 1 elem
REAL(iwp) :: Con,E,ny,Scat,Scad,tol,kappa
!-----------------------------------------------------
!   Read job information
!-----------------------------------------------------
OPEN (UNIT=11,FILE='prog81.dat',FORM='FORMATTED') !  Input
OPEN (UNIT=12,FILE='prog81.res',FORM='FORMATTED') !  Output
Call Jobin(Title,Cdim,N_dof,Toa,Nreg,Ltyp,Con,E,ny,&
           Isym,nodel,nodes,maxe)
Nsym= 2**Isym   !   number of symmetry loops
ALLOCATE(xP(Cdim,Nodes))   !  Array for node coordinates
xP = 0.0_iwp
ALLOCATE(Inci(Maxe,Nodel)) !  Array for incidences
Inci = 0
CALL Geomin(Nodes,Maxe,xp,Inci,Nodel,Cdim)
CALL flush(12)
Ndofe= Nodel*N_dof   !    Total degrees of freedom of element
ALLOCATE(BCode(Maxe,Ndofe))      !    Element Boundary codes
ALLOCATE(Elres_u(Maxe,Ndofe),Elres_t(Maxe,Ndofe))       
BCode = 0 ; Elres_u = 0.0_iwp ; Elres_t = 0.0_iwp
CALL BCinput(Elres_u,Elres_t,Bcode,nodel,ndofe,n_dof)
READ(11,*) tol,its,ell,kappa    ! data for bicgstab(l)   
ALLOCATE(Ldest(maxe,Ndofe))  ! Elem. destination vector
ALLOCATE(Ndest(Nodes,N_dof))
Ldest = 0; Ndest = 0
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
ALLOCATE(dTe(Ndofs,Ndofe),dUe(Ndofs,Ndofe))        ! Elem. coef. matrices
ALLOCATE(Diag(Ndofs,N_dof))                        ! Diagonal coefficients
ALLOCATE(Lhs(Ndofs,Ndofs),F(Ndofs),u1(Ndofs))      ! global arrays
ALLOCATE(Elcor(Cdim,Nodel))             !  Elem. Coordinates
ALLOCATE(Fac(Ndofe))                 !  Factor for symmetric elements
ALLOCATE(Incie(Nodel))               !  Element incidences
ALLOCATE(Ldeste(Ndofe))                      !  Element destination
ALLOCATE(Elres_te(Ndofe),Elres_ue(Ndofe))    ! Tractions of Element
!----------------------------------------------------------------
!  Compute element coefficient matrices
!----------------------------------------------------------------
Lhs(:,:) = 0.0_iwp; F(:) = 0.0_iwp; u1(:) = 0.0_iwp; Diag(:,:) = 0.0_iwp
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
                  CALL Mirror(Isym,nsy,Nodes,Elcor,Fac,   &
                  Incie,Ldeste,Elres_te,Elres_ue,nodel,n_dof,Cdim) 
                END IF
                IF(Cdim == 2) THEN
                        IF(N_dof == 1) THEN
                          CALL Integ2P(Elcor,Incie,Nodel,Nodes, &
                          xP,Con,dUe,dTe,Ndest,Isym)
                        ELSE
                          CALL Integ2E(Elcor,Incie,Nodel,Nodes, &
                          xP,E,ny,dUe,dTe,Ndest,Isym)
                        END IF
                ELSE
                   CALL Integ3(Elcor,Incie,Nodel,Nodes,xP,N_dof &
                                    ,E,ny,Con,dUe,dTe,Ndest,Isym)    
                END IF
                CALL Assembly(Lhs,F,DTe,DUe,Ldeste,BCode(Nel,:),Ncode &
                        ,Elres_u(Nel,:),Elres_te,Diag,Ndofe,N_dof,Nodel,Fac)     
        END DO &
        Symmetry_loop
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
!  Add Diagonal coefficients     !  penalties?    
!-------------------------------------------------------------
DO m=1,Ndofs            ! Loop over collocation points
   Nod=0
   DO n=1, Nodes
      DO l=1,N_dof
         IF (m == Ndest(n,l))THEN
             Nod=n     ;                EXIT
         END IF  
      END DO
      IF (Nod /= 0)EXIT
   END DO
   DO k=1,N_dof
      DoF=Ndest(Nod,k)
      IF(DoF /= 0) THEN
         IF(NCode(DoF) == 1) THEN
            Nel=0    ;    Pos=0
            DO i=1,Maxe
               DO j=1,Ndofe
                  IF(DoF == Ldest(i,j))THEN
                    Nel=i   ;    Pos=j  ;    EXIT
                  END IF
               END DO
            IF(Nel /= 0)EXIT
            END DO
            F(m) = F(m) - Diag(m,k) * Elres_u(Nel,Pos)
         ELSE
           Lhs(m,DoF)= Lhs(m,DoF) + Diag(m,k)
         END IF
       END IF
    END DO
END DO
!---------------------------------------------------------
!   Solve system of equations iteratively
!---------------------------------------------------------
CALL bicgstab_l(Lhs,F,Ndofs,u1,0.0_iwp,tol,its,ell,kappa)
!   Gather Element results from global result vector u1
Elements_2:     &
DO nel=1,maxe,maxe - 1  
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
Elements_2
END PROGRAM General_purpose_BEM

