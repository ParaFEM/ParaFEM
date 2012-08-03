  PROGRAM p124      
!------------------------------------------------------------------------------
!      Program 12.4 conduction equation using 8-node hexahedral elements;
!      pcg version implicit; integration in time using 'theta' method
!      parallel version
!------------------------------------------------------------------------------

 USE precision  ; USE global_variables ; USE mp_interface
 USE input      ; USE output           ; USE loading
 USE timing     ; USE maths            ; USE gather_scatter
 USE partition  ; USE elements         ; USE steering
 USE geometry   ; USE pcg
 
 IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

! neq,ntot are now global variables - not declared

 INTEGER, PARAMETER  :: nodof=1,ndim=3
 INTEGER             :: nod
 INTEGER             :: nxe,nye,nze,nn,nr,nip,neq_temp,nn_temp,i,j
 INTEGER             :: k,iel,nstep,npri,nres,iters,limit,it,is,nlen
 INTEGER             :: loaded_nodes,fixed_freedoms
 INTEGER             :: argc,iargc,meshgen,partitioner
 INTEGER             :: nels,ndof,ielpe,npes_pp
 REAL(iwp)           :: aa,bb,cc,kx,ky,kz,det,theta,dtim,real_time
 REAL(iwp)           :: val0 = 100.0_iwp
 REAL(iwp)           :: tol,alpha,beta,up,big
 REAL(iwp),PARAMETER :: zero = 0.0_iwp
 CHARACTER(LEN=15)   :: element
 CHARACTER(LEN=50)   :: fname,job_name,label
 CHARACTER(LEN=50)   :: program_name='p124'
 LOGICAL             :: converged = .false.

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------ 

 REAL(iwp),ALLOCATABLE :: loads_pp(:),u_pp(:),p_pp(:),points(:,:),kay(:,:)
 REAL(iwp),ALLOCATABLE :: coord(:,:),fun(:),jac(:,:),der(:,:),deriv(:,:)
 REAL(iwp),ALLOCATABLE :: weights(:),d_pp(:),kc(:,:),pm(:,:),funny(:,:)
 REAL(iwp),ALLOCATABLE :: p_g_co_pp(:,:,:),storka_pp(:,:,:)
 REAL(iwp),ALLOCATABLE :: storkb_pp(:,:,:),x_pp(:),xnew_pp(:),pmul_pp(:,:)
 REAL(iwp),ALLOCATABLE :: utemp_pp(:,:),diag_precon_pp(:),diag_precon_tmp(:,:)
 REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:),timest(:)
 INTEGER,ALLOCATABLE   :: rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp(:,:)
 
!------------------------------------------------------------------------------
! 3. Read job_name from the command line. 
!    Read control data, mesh data, boundary and loading conditions
!------------------------------------------------------------------------------ 

  ALLOCATE(timest(25))
  timest    = zero
  timest(1) = elap_time()

  CALL find_pe_procs(numpe,npes)
  
  PRINT *, "FIND_PE_PROCS on processor ", numpe, " of ", npes

  argc = iargc()
  IF(argc /= 1) CALL job_name_error(numpe,program_name)
  CALL GETARG(1,job_name)

  CALL read_p124(job_name,numpe,dtim,element,fixed_freedoms,kx,ky,kz,limit,   &
                 loaded_nodes,meshgen,nels,nip,nn,nod,npri,nr,nstep,          &
                 partitioner,theta,tol)

  CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)

  ndof = nod*nodof
  ntot = ndof

  ALLOCATE(g_num_pp(nod,nels_pp))
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp))
  ALLOCATE(rest(nr,nodof+1))

  g_num_pp   = 0
  g_coord_pp = zero
  rest       = 0

  timest(2) = elap_time()

  CALL read_g_num_pp2(job_name,iel_start,nn,npes,numpe,g_num_pp)
  timest(3) = elap_time()

  IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
  timest(4) = elap_time()

  CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
  timest(5) = elap_time()

  CALL read_rest(job_name,numpe,rest)
  timest(6) = elap_time()

  IF(numpe==1) PRINT *, " *** Read input data in: ", timest(6)-timest(1)," s"
  
! nn_temp=0
! neq_temp=0
! nye=nels/nxe/nze
! nr=(nxe+1)*(nye+1)+(nxe+1)*nze+nye*nze
! ielpe=iel_start

!------------------------------------------------------------------------------
! 4. Allocate dynamic arrays used in main program
!------------------------------------------------------------------------------

 ALLOCATE (points(nip,ndim),weights(nip),kay(ndim,ndim),                      &
           coord(nod,ndim),fun(nod),jac(ndim,ndim),der(ndim,nod),g(ntot),     &
           deriv(ndim,nod),pm(ntot,ntot),                                     &
           kc(ntot,ntot),funny(1,nod),num(nod),                               &
           g_g_pp(ntot,nels_pp),storka_pp(ntot,ntot,nels_pp),                 &
           utemp_pp(ntot,nels_pp),storkb_pp(ntot,ntot,nels_pp),               &
           pmul_pp(ntot,nels_pp))

  IF(numpe==1) PRINT *, " *** Allocated dynamic arrays in: ",                 &
                          elap_time()-timest(6)," s"

!------------------------------------------------------------------------------
! 5. Loop the elements to find the steering array and the number of equations
!    to solve.
!------------------------------------------------------------------------------

  CALL rearrange_2(rest)  

  g_g_pp = 0

  elements_1: DO iel = 1, nels_pp
    CALL find_g4(g_num_pp(:,iel),g_g_pp(:,iel),rest)
  END DO elements_1
   
  neq = 0

  elements_2: DO iel = 1, nels_pp
    i = MAXVAL(g_g_pp(:,iel))
    IF(i > neq) neq = i
  END DO elements_2

  neq  = MAX_INTEGER_P(neq)

  timest(7) = elap_time()

!------------------------------------------------------------------------------
! 6. Create interprocessor communication tables
!------------------------------------------------------------------------------

  CALL calc_neq_pp          
  CALL calc_npes_pp(npes,npes_pp)
  CALL make_ggl2(npes_pp,npes,g_g_pp)

  nres = nxe*(nze-1) + 1

  DO i = 1,neq_pp
    IF(nres==ieq_start+i-1) THEN
      it = numpe; is = i
    END IF
  END DO

  timest(8) = elap_time()

  IF(numpe==1) PRINT *, " *** Created ggl in: ", timest(8)-timest(7), " s"

!------------------------------------------------------------------------------
! 7. Allocate arrays dimensioned by neq_pp 
!------------------------------------------------------------------------------

  ALLOCATE(loads_pp(neq_pp),diag_precon_pp(neq_pp),u_pp(neq_pp),              &
           d_pp(neq_pp),p_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp))

  loads_pp  = zero ; diag_precon_pp = zero ; u_pp = zero
  d_pp      = zero ; p_pp           = zero ; x_pp = zero ; xnew_pp = zero

  timest(9) = elap_time()

  IF(numpe==1) PRINT *, " *** Allocated arrays dimensioned by neq_pp in: ",   &
                          timest(9)-timest(8), " s"  

  
! kay=0.0_iwp
! kay(1,1)=kx
! kay(2,2)=ky
! kay(3,3)=kz
! CALL box_bc8(nxe,nye,nze,rest)

!------------------------------------------------------------------------------
! 8. Element stiffness integration and storage
!------------------------------------------------------------------------------

 CALL sample(element,points,weights)

 storka_pp = zero 
 storkb_pp = zero
 kay       = zero
 kay(1,1)  = kx
 kay(2,2)  = ky
 kay(3,3)  = kz
 
 elements_3: DO iel=1,nels_pp

   kc = zero ; pm = zero

   gauss_pts: DO i=1,nip
     CALL shape_der(der,points,i)
     CALL shape_fun(fun,points,i)
     funny(1,:) = fun(:)
     jac        = MATMUL(der,g_coord_pp(:,:,iel))
     det        = determinant(jac)
     CALL invert(jac)
     deriv      = MATMUL(jac,der)
     kc         = kc + MATMUL(MATMUL(TRANSPOSE(deriv),kay),deriv)*det*weights(i)
     pm         = pm + MATMUL(TRANSPOSE(funny),funny)*det*weights(i) 
   END DO gauss_pts
   
   storka_pp(:,:,iel)=pm+kc*theta*dtim
   storkb_pp(:,:,iel)=pm-kc*(1._iwp-theta)*dtim
   
 END DO elements_3

 timest(10) = elap_time()
 
!------------------------------------------------------------------------------
! 9. Build the diagonal preconditioner
!------------------------------------------------------------------------------ 
  
  ALLOCATE(diag_precon_tmp(ntot,nels_pp))
  diag_precon_tmp = zero

  elements_4: DO iel = 1,nels_pp 
    DO k=1,ntot
      diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel)+storka_pp(k,k,iel)
    END DO
  END DO elements_4
 
 CALL scatter(diag_precon_pp,diag_precon_tmp)

 DEALLOCATE(diag_precon_tmp)

 diag_precon_pp=1._iwp/diag_precon_pp ! needs moving

!------------------------------------------------------------------------------
! 10. Read in the initial conditions and assign to equations
!------------------------------------------------------------------------------

 loads_pp = val0    ! needs to be read in from file
 pmul_pp  = .0_iwp

 IF(numpe==it)THEN
   fname = job_name(1:INDEX(job_name, " ")-1) // ".res"
   OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
   WRITE(11,'(A)')"    Time     Pressure  Iterations"
 END IF

!------------------------------------------------------------------------------
! 11. Time-stepping loop
!------------------------------------------------------------------------------

 timesteps: DO j=1,nstep
 
   real_time = j*dtim
   u_pp      = zero

   CALL gather(loads_pp,pmul_pp)
   elements_5: DO iel=1,nels_pp   
     utemp_pp(:,iel)=MATMUL(storkb_pp(:,:,iel),pmul_pp(:,iel))
   END DO elements_5
   CALL scatter(u_pp,utemp_pp)

   loads_pp=u_pp

!------------------------------------------------------------------------------
! 12. Solve simultaneous equations by pcg
!------------------------------------------------------------------------------

   d_pp = diag_precon_pp*loads_pp
   p_pp = d_pp
   x_pp = zero

   iters = 0

   iterations: DO 

     iters   = iters+1
     u_pp    = zero
     pmul_pp = zero

     CALL gather(p_pp,pmul_pp)  
     elements_6: DO iel=1,nels_pp
       utemp_pp(:,iel)=MATMUL(storka_pp(:,:,iel),pmul_pp(:,iel)) 
     END DO elements_6
     CALL scatter(u_pp,utemp_pp)

!------------------------------------------------------------------------------
! 13. PCG equation solution
!------------------------------------------------------------------------------

     up       = DOT_PRODUCT_P(loads_pp,d_pp)
     alpha    = up/DOT_PRODUCT_P(p_pp,u_pp)
     xnew_pp  = x_pp+p_pp*alpha
     loads_pp = loads_pp-u_pp*alpha
     d_pp     = diag_precon_pp*loads_pp
     beta     = DOT_PRODUCT_P(loads_pp,d_pp)/up
     p_pp     = d_pp+p_pp*beta
     u_pp     = xnew_pp

     CALL checon_par(xnew_pp,tol,converged,x_pp)
     IF(converged.OR.iters==limit)EXIT

   END DO iterations

   loads_pp=xnew_pp

   IF(j/npri*npri==j.AND.numpe==1)WRITE(11,'(2E12.4,I5)')real_time,      &
     loads_pp(is),iters

 END DO timesteps

! Smith and Griffiths output

! DO i=1,neq_pp
!   IF(nres==ieq_start+i-1)THEN
!     it=numpe
!     is=i
!   END IF
! END DO

 IF(numpe==it)THEN
   WRITE(11,'(A,I5,A)')"This job ran on  ",npes,"  processors"  
   WRITE(11,'(A)')"Global coordinates and node numbers"
   DO i=1,nels_pp,nels_pp - 1                                        
     WRITE(11,'(A,I8)')"Element ",i
     num=g_num_pp(:,i)
     DO k=1,nod
       WRITE(11,'(A,I8,3E12.4)')"  Node",num(k),p_g_co_pp(k,:,i)
     END DO
   END DO
   WRITE(11,'(A,3(I8,A))')"There are ",nn," nodes",nr," restrained and",  &
     neq," equations"
   WRITE(11,*)"Time after setup  is   :",elap_time()-timest(1)
 END IF

 IF(numpe==it) WRITE(11,*)"This analysis took  :",elap_time()-timest(1)

 CALL shutdown() 

END PROGRAM p124
