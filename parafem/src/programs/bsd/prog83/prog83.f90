PROGRAM PARALLEL_BEM
!------------------------------------------------------------------------------
!  General purpose BEM program for solving elasticity problems 
!  This version is parallel and uses BiCGSTAB(l). It must be linked to an
!  MPI library.
!------------------------------------------------------------------------------
USE bem_lib_p; USE precision; USE timing; USE maths; USE mp_interface
USE global_variables; USE gather_scatter; USE bicg
IMPLICIT NONE  !  Ndof changed to N_dof,Maxe is nels
INTEGER, ALLOCATABLE :: Inci(:,:)  !  Element Incidences
INTEGER, ALLOCATABLE :: BCode(:,:), NCode(:) !  Element BC´s
INTEGER, ALLOCATABLE :: Ldest(:,:) !  Element destination vector
INTEGER, ALLOCATABLE :: Ndest(:,:) !  Node destination vector
REAL(iwp), ALLOCATABLE :: Elres_u(:,:)  !  Element results , u
REAL(iwp), ALLOCATABLE :: Elres_t(:,:)  !  Element results , t
REAL(iwp), ALLOCATABLE :: Elcor(:,:)    !  Element coordinates
REAL(iwp), ALLOCATABLE :: xP(:,:)       !  Node co-ordinates
REAL(iwp), ALLOCATABLE :: dUe(:,:),dTe(:,:),lhs(:,:),Diag(:,:),pmul(:)
REAL(iwp), ALLOCATABLE :: km(:,:),qmul(:),Diag1(:,:) 
REAL(iwp), ALLOCATABLE :: store_dUe_pp(:,:,:),store_dTe_pp(:,:,:)
REAL(iwp), ALLOCATABLE :: F(:),F1(:)     !  global RHS
REAL(iwp), ALLOCATABLE :: u1(:),y_cop(:) ! global vector of unknowns 
CHARACTER (LEN=80) :: Title
CHARACTER (LEN=50) :: job_name='prog83'
INTEGER :: Cdim,m,n,Nodel,Nel,N_dof,Toa,N_tot
INTEGER :: Nreg,Ltyp,Nodes,Maxe,Ndofe,Ndofs               
INTEGER :: nod,nd,i,j,k,l,DoF,Pos,Isym,nsym,nsy
INTEGER :: its,iters,ell,nels,ielpe
INTEGER :: partitioner=1
REAL(iwp),ALLOCATABLE    :: Fac(:)     !  Factors for symmetry
REAL(iwp),ALLOCATABLE    :: Elres_te(:),Elres_ue(:)   
INTEGER,ALLOCATABLE :: Incie(:)   !  Incidences for one element
INTEGER,ALLOCATABLE :: Ldeste(:),g(:)  !  Destination vector 1 elem
REAL(iwp) :: Con,E,ny,Scat,Scad,tol,kappa,alpha,beta,rho,gama,       &
             omega,norm_r,r0_norm,error,one=1._iwp,zero=.0_iwp
LOGICAL:: converged
REAL(iwp),ALLOCATABLE::s(:),GG(:,:),Gamma(:),rt(:),y(:),y1(:),r(:,:),&
                       uu(:,:),timest(:)
ALLOCATE(timest(20))
timest=0.0_iwp
timest(1) = elap_time(); CALL find_pe_procs(numpe,npes)
!-----------------------------------------------------
!   Read job information
!-----------------------------------------------------
OPEN (UNIT=11,FILE='job.dat',FORM='FORMATTED',ACTION='READ') !  Input
OPEN (UNIT=21,FILE='coord.dat',FORM='FORMATTED',ACTION='READ') !  Input
OPEN (UNIT=31,FILE='plane.dat',FORM='FORMATTED',ACTION='READ') !  Input
OPEN (UNIT=41,FILE='load.dat',FORM='FORMATTED',ACTION='READ') !  Input
IF(numpe==1)OPEN(UNIT=12,FILE='prog83.res',FORM='FORMATTED',ACTION='WRITE')!O/P
IF(numpe==1) WRITE(12,*) "This job ran on ",npes," processors"
Call Jobin(Title,Cdim,N_dof,Toa,Nreg,Ltyp,Con,E,ny,&
           Isym,nodel,nodes,nels)
Nsym= 2**Isym   !   number of symmetry loops
ALLOCATE(xP(Cdim,Nodes))   !  Array for node coordinates
ALLOCATE(Inci(nels,Nodel)) !  Array for incidences
CALL Geomin(Nodes,nels,xp,Inci,Nodel,Cdim)
Ndofe= Nodel*N_dof   !    Total degrees of freedom of element
ALLOCATE(BCode(nels,Ndofe))      !    Element Boundary codes
ALLOCATE(Elres_u(nels,Ndofe),Elres_t(nels,Ndofe))       
CALL BCinput(Elres_u,Elres_t,Bcode,nodel,ndofe,n_dof)
READ(11,*) tol,its,ell,kappa      ! BiCGStab data
ALLOCATE(Ldest(nels,Ndofe))  ! Elem. destination vector
ALLOCATE(Ndest(Nodes,N_dof))
!---------------------------------------------------------------------
!     Determine Node destination vector and Element destination vector 
!---------------------------------------------------------------------
CALL Destination(Isym,Ndest,Ldest,xP,Inci,Ndofs,nodes,N_dof,Nodel,nels)
!---------------------------------------------------------------------
!     Determine global Boundary code vector
!---------------------------------------------------------------------
ALLOCATE(NCode(Ndofs))   
CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)
IF(numpe==1) WRITE(12,*) "Elements on first processor ",nels_pp
NCode=0
DoF_o_System: &
DO  nd=1,Ndofs
    DO Nel=1,nels 
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
ALLOCATE(store_dTe_pp(Ndofs,Ndofe,nels_pp),        &
         store_dUe_pp(Ndofs,Ndofe,nels_pp)) ! store el matrices on procs
ALLOCATE(Diag(Ndofs,N_dof),Diag1(Ndofs,N_dof))!Diag cos 
ALLOCATE(F(Ndofs),u1(Ndofs),F1(Ndofs))         ! global arrays
ALLOCATE(Elcor(Cdim,Nodel))             !  Elem. Coordinates
ALLOCATE(Fac(Ndofe))                 !  Factor for symmetric elements
ALLOCATE(Incie(Nodel))               !  Element incidences
ALLOCATE(Ldeste(Ndofe),pmul(Ndofe),km(N_dof,N_dof),g(N_dof),qmul(N_dof))!dest.
ALLOCATE(Elres_te(Ndofe),Elres_ue(Ndofe))    ! Tractions of Element
!----------------------------------------------------------------
!  Compute element coefficient matrices
!----------------------------------------------------------------
Lhs=zero; F1 = zero; u1 = zero; Diag1 = zero;    N_tot = Ndofs*N_dof
store_dUe_pp = zero; store_dTe_pp = zero    ; ielpe = iel_start
Elements_1:&
DO Nel=1,nels_pp  
        Symmetry_loop:&
        DO nsy= 1,Nsym
                Elcor(:,:)= xP(:,Inci(ielpe,:)) !  gather element coordinates
                Incie= Inci(ielpe,:)             !   incidences
                Ldeste= Ldest(ielpe,:)           !   and destinations
                Fac(1:nodel*n_dof)= 1.0_iwp
                Elres_te(:)=Elres_t(ielpe,:)
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
                CALL rhs_and_diag(F1,DTe,DUe,Ldeste,BCode(ielpe,:),Ncode &
                        ,Elres_u(ielpe,:),Elres_te,Diag1,Ndofe,N_dof,Nodel,Fac)
        END DO &
        Symmetry_loop           ;    ielpe = ielpe + 1
        store_dUe_pp(:,:,Nel) = dUe;    store_dTe_pp(:,:,Nel) = dTe
END DO &
Elements_1
CALL MPI_ALLREDUCE(F1,F,Ndofs,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
CALL MPI_ALLREDUCE(Diag1,Diag,N_tot,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
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
           DO i=1,nels
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
IF(numpe==1) WRITE(12,*) "Time before eq solution is ", elap_time()-timest(1)
!---------------------------------------------------------
!   Solve system of equations  element by element
!---------------------------------------------------------
ALLOCATE(s(ell+1),GG(ell+1,ell+1),Gamma(ell+1),y_cop(Ndofs),              &
            rt(Ndofs),y(Ndofs),y1(Ndofs),r(Ndofs,ell+1),uu(Ndofs,ell+1))
!      initialisation phase
   u1 = zero ; y = u1 ; y_cop = y; y1 = zero ; neq = Ndofs ; ielpe = iel_start
Elements_2 : DO Nel = 1 , nels_pp 
     Dte = store_DTe_pp(:,:,Nel);  DUe = store_DUe_pp(:,:,Nel)
     Ldeste = Ldest(ielpe,:); pmul = y(Ldeste)
     CALL form_lhs(lhs,DTe,DUe,Ldeste,BCode(Nel,:),Ncode,Ndofe,Fac)
     y1 = y1 + MATMUL(lhs,pmul) ;  ielpe = ielpe + 1   
END DO Elements_2 
CALL MPI_ALLREDUCE(y1,y,neq,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)  
DO i = 1 , nodes   
     CALL get_km(Cdim,i,y_cop,Diag,g,qmul,km)
     y(g) = y(g) + MATMUL(km,qmul)   
END DO    
    rt = F - y
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
               y1 = zero     ; y_cop = y  ;  ielpe = iel_start   
         Elements_3: DO Nel = 1 , nels_pp 
               Dte = store_DTe_pp(:,:,Nel);  DUe = store_DUe_pp(:,:,Nel)
               Ldeste = Ldest(ielpe,:); pmul = y(Ldeste)
               CALL form_lhs(lhs,DTe,DUe,Ldeste,BCode(Nel,:),Ncode,Ndofe,Fac)
               y1 = y1 + MATMUL(lhs,pmul)  ;  ielpe = ielpe + 1 
         END DO Elements_3 
         CALL MPI_ALLREDUCE(y1,y,neq,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)   
         DO i = 1 , nodes    
               CALL get_km(Cdim,i,y_cop,Diag,g,qmul,km)
               y(g) = y(g) + MATMUL(km,qmul)  
         END DO     
               uu(:,j+1) = y
               gama = DOT_PRODUCT(rt,y); alpha = rho/gama
               u1=u1+ alpha * uu(:,1)
               r(:,1:j) = r(:,1:j) - alpha * uu(:,2:j+1)    
               y = r(:,j)
               y1 = zero    ; y_cop = y ;   ielpe = iel_start
         Elements_4: DO Nel = 1 , nels_pp 
                Dte = store_DTe_pp(:,:,Nel);  DUe = store_DUe_pp(:,:,Nel)
                Ldeste = Ldest(ielpe,:);  pmul = y(Ldeste)
                CALL form_lhs(lhs,DTe,DUe,Ldeste,BCode(Nel,:),Ncode,Ndofe,Fac)
                y1 = y1 + MATMUL(lhs,pmul) ;   ielpe = ielpe + 1   
         END DO Elements_4 
         CALL MPI_ALLREDUCE(y1,y,neq,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier) 
         DO i = 1 , nodes   
                CALL get_km(Cdim,i,y_cop,Diag,g,qmul,km)     
                y(g) = y(g) + MATMUL(km,qmul)               
         END DO   
                r(:,j+1) = y
    END DO
    GG = MATMUL(TRANSPOSE(r),r)
    CALL form_s(gg,ell,kappa,omega,gamma,s)
    u1 = u1 - MATMUL(r,s);r(:,1)=MATMUL(r,Gamma)
    uu(:,1)=MATMUL(uu,Gamma)
    norm_r = norm(r(:,1))  ;  error = norm_r/r0_norm       
END DO iterations
  IF(numpe==1) WRITE(12,'(/,A,I5,A,/)')      &
                    "It took BiCGSTAB_L ",iters," iterations to converge"
!   Gather Element results from global result vector u1
Elements_5:     &
DO nel=1,nels , nels - 1  
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
   IF(numpe==1) WRITE(12,'(24F12.5)') (Elres_u(nel,m), m=1,Ndofe)
   IF(numpe==1) WRITE(12,'(24F12.5)') (Elres_t(nel,m), m=1,Ndofe)
END DO &
Elements_5
IF(numpe==1) WRITE(12,*) "This analysis took ", elap_time() - timest(1)
CALL shutdown()
END PROGRAM PARALLEL_BEM
