PROGRAM p123         
!------------------------------------------------------------------------------
!      program 7.5 three dimensional analysis of Laplace's equation
!      using 8-node brick elements, preconditioned conjugate gradient solver
!      only integrate one element , diagonal preconditioner diag_precon
!      parallel version  ;   central loaded or fixed freedom   ;  box_bc  
!------------------------------------------------------------------------------
  USE new_library; USE geometry_lib; USE precision; USE utility;USE mp_module
  USE  timing    ; USE global_variables1;  USE gather_scatter6;  IMPLICIT NONE
 !  ndof, nels, neq , ntot  are now global variables - not declared
 INTEGER::nxe,nye,nze,nn,nr,nip,nodof=1,nod=8, nres, is , it ,       &
          i,j,k,ndim=3,iters,limit,iel,num_no,no_index_start,        &
          neq_temp,nn_temp , loaded_freedoms, fixed_freedoms
 REAL(iwp)::aa,bb,cc,kx,ky,kz,det,tol,up,alpha,beta,q,penalty=1.e20_iwp 
 CHARACTER(LEN=15)::element= 'hexahedron';    LOGICAL :: converged          
!-------------------------- dynamic arrays-------------------------------------
 REAL(iwp),ALLOCATABLE :: points(:,:),kc(:,:),coord(:,:), weights(:),        &
                         p_g_co_pp(:,:,:), jac(:,:), der(:,:), deriv(:,:),   &
                         col(:,:),row(:,:),kcx(:,:),kcy(:,:),kcz(:,:),       &
                         diag_precon_pp(:),p_pp(:),r_pp(:),x_pp(:),          &
                         xnew_pp(:),u_pp(:),pmul_pp(:,:),utemp_pp(:,:),      &
                         d_pp(:),diag_precon_tmp(:,:),val(:),val_f(:),       &
                         store_pp(:),eld(:)
 INTEGER, ALLOCATABLE :: rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp(:,:),no(:),&
                         no_f(:),no_local_temp(:),no_local_temp_f(:),no_local(:)
!-------------------------input and initialisation---------------------------
 timest(1) = elap_time( )   ;    CALL find_pe_procs(numpe,npes)
 IF (numpe==npes) THEN
  OPEN (10,FILE='p123.dat',STATUS=    'OLD',ACTION='READ')
  READ (10,*) nels,nxe,nze,nip,aa,bb,cc,kx,ky,kz, tol,limit ,             &
              loaded_freedoms,fixed_freedoms
 END IF
 CALL bcast_inputdata_p123(numpe,npes,nels,nxe,nze,nip,aa,bb,cc,kx,ky,kz, &
                            tol,limit,loaded_freedoms,fixed_freedoms)
 CALL calc_nels_pp   ;   ndof=nod*nodof ; ntot=ndof ; nye = nels/nxe/nze  
      neq_temp = 0; nn_temp = 0 ; nr=(nxe+1)*(nye+1)+(nxe+1)*nze+nye*nze  
 ALLOCATE (points(nip,ndim),coord(nod,ndim),jac(ndim,ndim),kc(ntot,ntot),    &
           der(ndim,nod),deriv(ndim,nod),rest(nr,nodof+1),kcx(ntot,ntot),    &
           g(ntot),pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp),col(ntot,1), &
           p_g_co_pp(nod,ndim,nels_pp),g_num_pp(nod,nels_pp),weights(nip),   &
           num(nod),g_g_pp(ntot,nels_pp),no(1),kcy(ntot,ntot),val(1),       &
           no_local_temp(1),row(1,ntot),diag_precon_tmp(ntot,nels_pp),      &
           val_f(1),eld(ntot),no_f(1),no_local_temp_f(1),kcz(ntot,ntot))
      CALL box_bc8(nxe,nye,nze,rest); ielpe=iel_start; nres=nxe*(nze-1)+1
      IF(loaded_freedoms>0) THEN; no = nres;   val = 10.; END IF
      IF(fixed_freedoms>0)THEN; no_f=nres; val_f = 100. ; END IF 
      CALL sample(element,points,weights); CALL rearrange_2(rest)  
!---------- loop the elements for global cordinates etc ---------------------
     elements_0: DO iel = 1 , nels_pp
             CALL geometry_8bxz(ielpe,nxe,nze,aa,bb,cc,coord,num)
             CALL find_g4(num,g,rest) ; g_num_pp(:,iel) = num
             p_g_co_pp(:,:,iel) = coord; g_g_pp(:,iel) = g; ielpe = ielpe+1 
             i = MAXVAL(g); j = MAXVAL(num)
             IF(i>neq_temp)neq_temp = i; IF(j>nn_temp)nn_temp = j
     END DO elements_0
   neq = reduce(neq_temp)  ;  nn = reduce(nn_temp)
   CALL calc_neq_pp        ;  CALL make_ggl(g_g_pp); diag_precon_tmp = .0_iwp
   DO i=1,neq_pp; IF(nres==ieq_start+i-1) THEN; it = numpe;is = i; END IF
   END DO
   ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),         &
            u_pp(neq_pp),diag_precon_pp(neq_pp),d_pp(neq_pp),store_pp(neq_pp))
            r_pp = .0_iwp; p_pp = .0_iwp; x_pp = .0_iwp; xnew_pp = .0_iwp
            diag_precon_pp = .0_iwp ;   store_pp = .0_iwp
!---------- element stiffness integration and build the preconditioner---------
      iel=1;CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,num)
      kcx = .0_iwp; kcy = .0_iwp; kcz = .0_iwp
      gauss_pts_1:  DO i=1,nip
                CALL shape_der (der,points,i) ; jac = MATMUL(der,coord)
                det=determinant(jac);CALL invert(jac);deriv = MATMUL(jac,der)
                row(1,:) = deriv(1,:); eld=deriv(1,:); col(:,1) = eld
                kcx = kcx + MATMUL(col,row)*det*weights(i)
                row(1,:) = deriv(2,:); eld=deriv(2,:); col(:,1) = eld
                kcy = kcy + MATMUL(col,row)*det*weights(i)
                row(1,:) = deriv(3,:); eld=deriv(3,:); col(:,1) = eld
                kcz = kcz + MATMUL(col,row)*det*weights(i)
      END DO gauss_pts_1        ; kc = kcx*kx + kcy*ky + kcz*kz 
     elements_1: DO iel = 1,nels_pp  
       DO k=1,ntot;diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel)+kc(k,k);END DO    
     END DO elements_1
     CALL scatter(diag_precon_pp,diag_precon_tmp);DEALLOCATE(diag_precon_tmp)
     IF(numpe==it)THEN
      OPEN (11,FILE='p123.res',STATUS='REPLACE',ACTION='WRITE')              
      WRITE(11,'(A,I5,A)') "This job ran on ", npes , "  processors" 
      WRITE(11,'(A)') "Global coordinates and node numbers "
      DO i = 1 , nels_pp,nels_pp-1                                         
                 WRITE(11,'(A,I8)')"Element ",i ; num = g_num_pp(:,i)
      DO k = 1,nod;WRITE(11,'(A,I8,3E12.4)')                               &
                 "  Node",num(k),p_g_co_pp(k,:,i); END DO
      END DO
      WRITE(11,'(A,3(I8,A))') "There are ",nn," nodes",nr," restrained and",&
                           neq," equations"
      WRITE(11,*) "Time after setup is   :", elap_time( ) - timest(1)
     END IF
!-------------------- get starting r-----------------------------------------
    IF(loaded_freedoms>0) THEN
     CALL reindex_fixed_nodes(ieq_start,no,no_local_temp,num_no,no_index_start)
     ALLOCATE(no_local(1:num_no)) ; no_local = no_local_temp(1:num_no)
     DEALLOCATE(no_local_temp)
        DO i = 1 , num_no
         r_pp(no_local(i)-ieq_start+1) = val(no_index_start + i - 1)
        END DO
    END IF            ;     q = SUM_P(r_pp)
    IF(numpe==it)THEN  ;  WRITE(11,'(A,E12.4)') "The total load is ", q
    END IF
    IF(fixed_freedoms>0) THEN
     CALL reindex_fixed_nodes(ieq_start,no_f,no_local_temp_f,        &
                              num_no,no_index_start)
     ALLOCATE(no_local(1:num_no)) ; no_local = no_local_temp_f(1:num_no)
     DEALLOCATE(no_local_temp_f)
        DO i = 1 , num_no        ; j=no_local(i) - ieq_start + 1
         diag_precon_pp(j)=diag_precon_pp(j) + penalty 
         r_pp(j) = diag_precon_pp(j) *  val_f(no_index_start + i - 1)
         store_pp(j) =  diag_precon_pp(j)
        END DO
    END IF
    diag_precon_pp=1._iwp/diag_precon_pp;d_pp=diag_precon_pp*r_pp; p_pp=d_pp
!--------------------preconditioned c. g. iterations---------------------------
       iters = 0
     iterations  :      DO 
             iters = iters + 1     ;    u_pp = 0._iwp  ; pmul_pp = .0_iwp
       CALL gather(p_pp,pmul_pp) 
       elements_2 : DO iel = 1, nels_pp
                    utemp_pp(:,iel) = MATMUL(kc,pmul_pp(:,iel)) 
       END DO elements_2  ;       CALL scatter(u_pp,utemp_pp)
    IF(fixed_freedoms>0) THEN
        DO i = 1 , num_no;  j = no_local(i)-ieq_start+1
         u_pp(j)=p_pp(j) * store_pp(j)
        END DO
    END IF
!--------------------------pcg equation solution-------------------------------
           up=DOT_PRODUCT_P(r_pp,d_pp); alpha= up/ DOT_PRODUCT_P(p_pp,u_pp)
           xnew_pp = x_pp + p_pp* alpha ; r_pp=r_pp - u_pp*alpha
           d_pp = diag_precon_pp*r_pp ;  beta=DOT_PRODUCT_P(r_pp,d_pp)/up
           p_pp=d_pp+p_pp*beta    
           CALL checon_par(xnew_pp,x_pp,tol,converged,neq_pp)    
           IF(converged .OR. iters==limit) EXIT
     END DO iterations
     IF(numpe==it)THEN
       WRITE(11,'(A,I5)')"The number of iterations to convergence was  ",iters 
       WRITE(11,'(A)')   "The  potentials are   :"
       WRITE(11,'(A)') "   Freedom       Potential"
       DO i = 1 , 4
        WRITE(11,'(I5,A,E12.4)') nres+i-1, "     ", xnew_pp(is+i-1)
       END DO
     END IF
  IF(numpe==it) WRITE(11,*) "This analysis took   ", elap_time( ) - timest(1)
  CALL MPI_FINALIZE(ier)
 END PROGRAM p123
