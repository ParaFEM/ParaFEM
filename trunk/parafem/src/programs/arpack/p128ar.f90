PROGRAM p128d      
!------------------------------------------------------------------------------
!      program 10.4 eigenvalues and eigenvectors of a cuboidal elastic
!      solid in 3d using  uniform 8-node hexahedral elements ; Arnoldi "dp" 
!      for lumped mass this is done element by element : parallel version
!------------------------------------------------------------------------------
 USE new_library ; USE geometry_lib ; USE precision!; USE timing    
 USE utility; USE mp_module; USE global_variables1; USE gather_scatter6
 IMPLICIT NONE  ! ndof,nels,neq,ntot are global - not declared
 INTEGER::nxe,nye,nze,nn,nr,nip,nodof=3,nod=8,nst=6,i,j,k,iel,ndim=3,iters, &
          model,nconv,ncv,nev,                                              &
          ido,ierr,info,iparam(11),ipntr(11),ishfts,lworkl,maxitr  
 REAL(iwp)::aa,bb,cc,rho,e,nu,det  , sigma,tol
 REAL(iwp),PARAMETER :: zero = 0.0_iwp , one = 1.0_iwp  
 CHARACTER (LEN=15) :: element = 'hexahedron'
 CHARACTER :: bmat*1, which*2 ;  LOGICAL :: rvec
!--------------------------- dynamic arrays------------------------------------
 REAL(iwp),ALLOCATABLE:: points(:,:),dee(:,:),coord(:,:),vdiag(:),resid(:), &  
                         fun(:),jac(:,:),der(:,:),deriv(:,:),weights(:),  & 
                         bee(:,:),km(:,:),emm(:,:),ecm(:,:),utemp(:), &
                         diag(:),pmul(:),d(:,:),v(:,:),workl(:),workd(:),   &
                         p_g_co_pp(:,:,:),diag1(:),udiag(:)   
 INTEGER, ALLOCATABLE  :: rest(:,:), g(:), num(:) , g_num_pp(:,:) ,         &
                          g_g_pp(:,:)
 LOGICAL, ALLOCATABLE  :: select(:)   ; INCLUDE  'debug.h'
!----------------------input and initialisation-----------------------------
 ! timest(1) = elap_time( ) 
 CALL find_pe_procs(numpe,npes)
 IF(numpe==npes)  OPEN (11,FILE='p128d.res',STATUS='REPLACE',ACTION='WRITE')
 OPEN (10,FILE='p128d.dat',STATUS=    'OLD',ACTION='READ')
 READ (10,*) nels,nxe,nze,nip,aa,bb,cc,rho,e,nu,                          &
                    nev,ncv,bmat,which,tol,maxitr
 CALL calc_nels_pp ;  ndof=nod*nodof  ; ntot = ndof ; nn = 0; neq = 0
 nr = (nxe + 1) * (nze + 1)  ; nye = nels/nxe/nze ; lworkl = ncv*(ncv+8)
 ndigit = -3; logfil = 6; msgets = 0; mseupd = 0 ; msaitr = 0
 msaupd = 1;  msaup2 = 0; mseigt = 0; msapps = 0   
 info = 0;    ishfts = 1; model = 1 ; ido = 0    
 iparam(1) = ishfts;  iparam(3) = maxitr; iparam(7) = model
 ALLOCATE (rest(nr,nodof+1),points(nip,ndim),pmul(ntot),          &
            coord(nod,ndim),fun(nod),jac(ndim,ndim), weights(nip),           &
            g_num_pp(nod,nels_pp),der(ndim,nod),deriv(ndim,nod),dee(nst,nst),&
            num(nod),km(ntot,ntot),g(ntot),g_g_pp(ntot,nels_pp),             &
            ecm(ntot,ntot),utemp(ntot),d(ncv,2),workl(lworkl),select(ncv),   &
            p_g_co_pp(nod,ndim,nels_pp),bee(nst,ntot),emm(ntot,ntot))    
   rest = 0 ;  DO i=1,nr; rest(i,1) = i; END DO  ;ielpe = iel_start
!------------------ loop the elements to set up global arrays  ----------------
 elements_1 : DO iel =1,nels_pp
               CALL geometry_8bxz(ielpe,nxe,nze,aa,bb,cc,coord,num)
               CALL find_g(num,g,rest) ; g_num_pp(:,iel)=num
               p_g_co_pp(:,:,iel)=coord; g_g_pp(:,iel)=g ; ielpe = ielpe + 1
 END DO elements_1
 nn = (nxe+1)*(nze+1)*(nye+1)  ;  neq = (nn - nr)*nodof 
 IF(numpe==npes) THEN
    WRITE(11,'(A,I5,A)') "This job ran on ",npes, "   processors"
    WRITE(11,'(A)') "Global coordinates and node numbers"
    DO i=1,nels_pp,nels_pp-1; WRITE(11,'(A,I8)')"Element ",i
      num = g_num_pp(:,i)
      DO k=1,nod;WRITE(11,'(A,I8,3E12.4)')                                &
          "  Node",num(k),p_g_co_pp(k,:,i); END DO
    END DO
    WRITE(11,'(A,3(I8,A))') "There are ",nn," nodes",nr," restrained  and",  &
                            neq," equations"
!   WRITE(11,*) "Time after setup is  :", elap_time( ) - timest(1)
 END IF
   ALLOCATE  ( v(neq,ncv),diag(neq),udiag(neq),vdiag(neq),workd(3*neq),    &
               resid(neq),diag1(neq))   
   diag1 = zero    
   v     = zero
   diag  = zero
   udiag = zero
   vdiag = zero
   workd = zero
   resid = zero
   diag1 = zero

   CALL sample( element, points, weights); CALL deemat(dee,e,nu)
  
!--------------- element stiffness integration and assembly--------------------
  elements_2: DO iel=1,nels_pp     ! allows for different km,emm  
                coord = p_g_co_pp(:,:,iel) ;  g = g_g_pp(:,iel)
                km = zero ; emm = zero
                integrating_pts_1:  DO i=1,nip
                  CALL shape_fun(fun,points,i)
                  CALL shape_der(der,points,i); jac=MATMUL(der,coord)
                  det= determinant(jac) ; CALL invert(jac) 
                  deriv = MATMUL(jac,der);CALL beemat(bee,deriv)
                  km= km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
                  CALL ecmat(ecm,fun,ntot,nodof);emm=emm+ecm*det*weights(i)*rho
                END DO integrating_pts_1
                DO k=1,ntot
                   IF(g(k)/=0) diag1(g(k))=diag1(g(k))+SUM(emm(k,:))
                END DO
  END DO elements_2
  CALL MPI_ALLREDUCE(diag1,diag,neq,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)   
 PRINT*, "Element stiffness"
!------------------------------find eigenvalues--------------------------------
  iters = 0; diag = one / sqrt(diag) ! diag holds l**(-1/2) 
    DO
       iters = iters + 1
       PRINT *, "Iters = ", iters, " calling DSAUPD"
       CALL dsaupd(ido,bmat,neq,which,nev,tol,resid,                 &
                   ncv,v,neq,iparam,ipntr,workd,workl,lworkl,info) 
       IF( ido /=-1 .AND. ido /= 1) EXIT
       diag1 = zero;  vdiag = workd(ipntr(1):ipntr(1) + neq - 1) * diag
      elements_3 : DO iel = 1 , nels_pp
                      g = g_g_pp(:,iel)
                      DO i=1,ntot
                       IF(g(i)==0)pmul(i)=zero;IF(g(i)/=0)pmul(i)=vdiag(g(i))
                      END DO
                      utemp = MATMUL(km,pmul)   
!                      CALL dgemv('n',ntot,ntot,one,km,ntot,pmul,1,zero,utemp,1) 
                      DO i=1,ntot
                       IF(g(i)/=0) diag1(g(i)) = diag1(g(i)) + utemp(i)
                      END DO
      END DO elements_3
      CALL MPI_ALLREDUCE(diag1,udiag,neq,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      udiag = udiag * diag ; workd(ipntr(2):ipntr(2) + neq -1) = udiag
    END DO
    IF(numpe==npes) WRITE(11,'(A,I8,A)') "It took ",iters,"  iterations"
! Either we have convergence or there is an error
  IF(info < 0) THEN
    IF(numpe==npes) WRITE(11,*) "Fatal error in dsaupd "
  ELSE
   rvec = .TRUE.
   PRINT *, "Calling DSEUPD"
   CALL dseupd(rvec, 'All',select, d , v , neq, sigma,               &
               bmat, neq, which, nev, tol, resid,                    &
               ncv, v, neq,iparam,ipntr,workd, workl, lworkl, ierr)
   IF(ierr/=0) THEN
    IF(numpe==npes) THEN
     WRITE(11,*) "Fatal error in dseupd "
     WRITE(11,*) "Should compute residuals "
    END IF
   END IF
   IF(numpe==npes) THEN
     WRITE(11,'(A)') "The eigenvalues are  :"
     WRITE(11,'(5E12.4)') d(1:nev,1)
     WRITE(11,'(A)') "The eigenvectors are  :"
     DO i = 1 , nev
        udiag(:) = v(1:neq,i)  ; udiag = udiag * diag  
        WRITE(11,'("Eigenvector number  ",I4," is: ")')  i
        WRITE(11,'(6E12.4)') udiag(1:6)
     END DO
   END IF
  END IF
! IF(numpe==npes) WRITE(11,*) "This analysis took  :", elap_time( ) - timest(1)
 CALL shutdown() 
END PROGRAM p128d 
