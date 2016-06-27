PROGRAM xx16       
!------------------------------------------------------------------------- 
!      Program xx16 three dimensional analysis of an elastic solid
!      using 20-node brick elements, preconditioned conjugate gradient
!      solver; diagonal preconditioner diag_precon; parallel version
!      loaded_nodes only; tests native mode for MPI on Xeon Phi
!------------------------------------------------------------------------- 

!USE mpi_wrapper  !remove comment for serial compilation
 USE precision; USE global_variables; USE mp_interface; USE input
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE steering; USE new_library

!-------------------------------------------------------------------------
! 1. Modules not USE'd in the 5th edition of the book
!-------------------------------------------------------------------------

 USE omp_lib
 USE, INTRINSIC :: ISO_C_BINDING ! to output to C binary file

 IMPLICIT NONE

!-------------------------------------------------------------------------
! 2. Declare variables used in the main program
!-------------------------------------------------------------------------

! neq,ntot are now global variables - must not be declared

 INTEGER,PARAMETER::testcase=4
 INTEGER,PARAMETER::nodof=3,ndim=3,nst=6
 INTEGER::loaded_nodes,iel,i,j,k,iters,limit,nn,nr,nip,nod,nels,ndof,    &
   npes_pp,node_end,node_start,nodes_pp,meshgen,partitioner,nlen,        &
   section,fnum
 REAL(iwp),PARAMETER::zero=0.0_iwp
 REAL(iwp),PARAMETER::one=1.0_iwp
 REAL(iwp)::dot_temp
 REAL(iwp)::e,v,det,tol,up,alpha,beta,q; LOGICAL::converged=.false.
 REAL(iwp)::gigaflops,gigaflops_1,gigaflops_2
 REAL(iwp)::flop,flop_1,flop_2
 REAL(iwp):: DDOT ! BLAS MKL
 CHARACTER(LEN=50) :: argv
 CHARACTER(LEN=15) :: element
 CHARACTER(LEN=6)  :: ch 
 CHARACTER(LEN=80) :: cbuffer
!LOGICAL           :: io_binary=.false.
 LOGICAL           :: io_binary=.true. 
 LOGICAL           :: sym_storkm=.true.
!LOGICAL           :: sym_storkm=.false.

!-------------------------------------------------------------------------
! 3. Declare dynamic arrays
!-------------------------------------------------------------------------

 REAL(iwp),ALLOCATABLE::points(:,:),dee(:,:),weights(:),val(:,:),        &
   disp_pp(:),g_coord_pp(:,:,:),jac(:,:),der(:,:),deriv(:,:),bee(:,:),   &
   storkm_pp(:,:,:),eps(:),sigma(:),diag_precon_pp(:),p_pp(:),r_pp(:),   &
   x_pp(:),xnew_pp(:),u_pp(:),pmul_pp(:,:),utemp_pp(:,:),d_pp(:),        &
   timest(:),timest_min(:),timest_max(:),diag_precon_tmp(:,:),           &
   eld_pp(:,:),temp(:),utemp(:),pad1(:),pad2(:),km(:,:)
 REAL(iwp),ALLOCATABLE :: vstorkm_pp(:,:)
 INTEGER,ALLOCATABLE   :: rest(:,:),g_num_pp(:,:),g_g_pp(:,:),node(:)

!-------------------------------------------------------------------------
! 4. Read job name from the command line.
!    Read control data, mesh data, boundary and loading conditions
!-------------------------------------------------------------------------
 
 ALLOCATE(timest(24),timest_min(24),timest_max(24))
 timest=zero; timest_min=zero; timest_max=zero
 CALL find_pe_procs(numpe,npes)
 timest(1)=MPI_Wtime()
 CALL getname(argv,nlen) 
 CALL read_p121(argv,numpe,e,element,limit,loaded_nodes,meshgen,nels,    &
   nip,nn,nod,nr,partitioner,tol,v)
 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ndof=nod*nodof; ntot=ndof

 IF(sym_storkm) THEN
   j=ntot
   DO i=1,ntot-1
     j=j+ntot-i
   END DO
   ALLOCATE(vstorkm_pp(j,nels_pp),km(ntot,ntot))
 ELSE
   ALLOCATE(storkm_pp(ntot,ntot,nels_pp),km(ntot,ntot))
 END IF

 ALLOCATE(g_num_pp(nod, nels_pp),g_coord_pp(nod,ndim,nels_pp),           &
   rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
 CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)
 IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
 CALL read_rest(argv,numpe,rest)
 timest(2)=MPI_Wtime()-timest(1)

!-------------------------------------------------------------------------
! 5. Allocate dynamic arrays used in the main program
!-------------------------------------------------------------------------

 ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim),der(ndim,nod),    &
   deriv(ndim,nod),bee(nst,ntot),weights(nip),eps(nst),sigma(nst),       &
   g_g_pp(ntot,nels_pp),utemp(ntot))

!-------------------------------------------------------------------------
! 6. Loop the elements to find the steering array and the number of 
!    equations to solve
!-------------------------------------------------------------------------

 CALL rearrange(rest)
 g_g_pp = 0; neq = 0

 timest(15)=MPI_Wtime()

 elements_0: DO iel=1,nels_pp
   CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
 END DO elements_0

 timest(16)=MPI_Wtime()-timest(15)

 neq=MAXVAL(g_g_pp); neq=max_p(neq)

 DEALLOCATE(rest)

!-------------------------------------------------------------------------
! 7. Create interprocessor communication tables
!-------------------------------------------------------------------------

 CALL calc_neq_pp
 CALL calc_npes_pp(npes,npes_pp)
 CALL make_ggl(npes_pp,npes,g_g_pp)

!-------------------------------------------------------------------------
! 8. Allocate arrays dimensioned by neq_pp
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
! 9. Element stiffness integration and storage
!-------------------------------------------------------------------------

 timest(13)=MPI_Wtime()
 dee=zero; CALL deemat(dee,e,v); CALL sample(element,points,weights)

 IF(sym_storkm) THEN
   vstorkm_pp       = zero
   DO iel=1,nels_pp
     km = zero
     DO i=1,nip
       CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
       det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
       CALL beemat(bee,deriv)
       km(1:ntot,1:ntot)=km(1:ntot,1:ntot) +                            &
                   MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)   
     END DO
     j = 1
     DO i = 1,ntot
       k = j+ntot-i
       vstorkm_pp(j:k,iel)=km(i:,i)  
       j = k+1
     END DO
   END DO
 ELSE
   storkm_pp=zero
   elements_1: DO iel=1,nels_pp
     gauss_pts_1: DO i=1,nip
       CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
       det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
       CALL beemat(bee,deriv)
       storkm_pp(1:ntot,1:ntot,iel)=storkm_pp(1:ntot,1:ntot,iel) +             &
                      MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)   
     END DO gauss_pts_1
   END DO elements_1
 END IF
 
 timest(14)=MPI_Wtime()-timest(13)

!-------------------------------------------------------------------------
! 10. Build the diagonal preconditioner
!-------------------------------------------------------------------------

 timest(17)=MPI_Wtime()

 ALLOCATE(diag_precon_tmp(ntot,nels_pp)); diag_precon_tmp = zero
 ALLOCATE(diag_precon_pp(neq_pp));        diag_precon_pp  = zero
 
 IF(sym_storkm) THEN
   DO iel = 1,nels_pp 
     j = 1
     DO i = 1,ndof
       diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel) + vstorkm_pp(j,iel)
       j = j + ndof - i + 1
     END DO
   END DO
 ELSE   
   elements_2: DO iel=1,nels_pp ; DO i=1,ndof
     diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel)+storkm_pp(i,i,iel)
   END DO;  END DO elements_2
 END IF
 
 CALL scatter(diag_precon_pp,diag_precon_tmp); DEALLOCATE(diag_precon_tmp)
 timest(18)=MPI_Wtime()-timest(17)

 IF(numpe==1)THEN
   OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
   fnum = 11; section = 1
   CALL write_summary(fnum,neq,nn,npes,nr,numpe,section,testcase,timest)
 END IF 

 IF(numpe==1) PRINT *, "Write summary"
!-------------------------------------------------------------------------
! 11. Read in loaded nodes and get starting r_pp
!-------------------------------------------------------------------------

 ALLOCATE(r_pp(neq_pp)) ; r_pp = zero

 timest(19)=MPI_Wtime()
 IF(loaded_nodes>0) THEN
   ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes)); node=0; val=zero
   CALL read_loads(argv,numpe,node,val)
   CALL load(g_g_pp,g_num_pp,node,val,r_pp(1:)); q=SUM_P(r_pp(1:))
   IF(numpe==1) WRITE(11,'(A,E12.4)') " The total load is:               ",q
   DEALLOCATE(node,val)
 END IF
 timest(20)=MPI_Wtime()-timest(19)
 DEALLOCATE(g_g_pp)

 IF(numpe==1) PRINT *, "Read in loaded nodes"
!-------------------------------------------------------------------------
! 12. Invert the preconditioner
!-------------------------------------------------------------------------

 diag_precon_pp=1._iwp/diag_precon_pp

 IF(numpe==1) PRINT *, "Inverted preconditioner"
!-------------------------------------------------------------------------
! 13. Initialise preconditioned conjugate gradient solver
!-------------------------------------------------------------------------

 ALLOCATE(utemp_pp(ntot,nels_pp)) ; utemp_pp = zero
 ALLOCATE(pmul_pp(ntot,nels_pp))  ; pmul_pp  = zero
 ALLOCATE(p_pp(neq_pp))           ; p_pp     = zero
 ALLOCATE(x_pp(neq_pp))           ; x_pp     = zero
 ALLOCATE(xnew_pp(neq_pp))        ; xnew_pp  = zero
 ALLOCATE(u_pp(neq_pp))           ; u_pp     = zero
 ALLOCATE(d_pp(neq_pp))           ; d_pp     = zero

 IF(numpe==1) PRINT *, "Allocated neq_pp arrays"

 d_pp = diag_precon_pp*r_pp
 p_pp = d_pp
 x_pp = zero

!-------------------------------------------------------------------------
! 14. Preconditioned conjugate gradient iterations
!-------------------------------------------------------------------------

 iters=0

 timest(3)=MPI_Wtime()

 iterations: DO 
   iters=iters+1; u_pp=zero; pmul_pp=zero; utemp_pp=zero
   IF(numpe==1) PRINT *,"Iters = ",iters

   timest(4)=MPI_Wtime()
   CALL gather(p_pp,pmul_pp(1:ntot,:))
   timest(5)=timest(5)+MPI_Wtime()-timest(4)

!-------------------------------------------------------------------------
! 15. Elements loop
!    
!    Compare (i) MATMUL, (ii) DGEMV and a nested loop with (iii) DDOT
!
!    See sites.duke.edu/vamvanij/2013/08/29/hw2-matrix-multiplication
!-------------------------------------------------------------------------

   timest(6)=MPI_Wtime()

  SELECT CASE(testcase)

     CASE(1) ! Matmul

     DO iel=1,nels_pp
       utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
     END DO


     CASE(2) ! DGEMV

     DO iel=1,nels_pp
       CALL dgemv('n',ntot,ntot,1.0_iwp,storkm_pp(1:ntot,1:ntot,iel),ntot,    &
                  pmul_pp(1:ntot,iel),1,0.0_iwp,utemp_pp(1:ntot,iel),1)
     END DO


     CASE(3) ! DGEMV + OMP

     !$OMP PARALLEL
     !$OMP DO PRIVATE(iel)     
     DO iel=1,nels_pp
       CALL dgemv('n',ntot,ntot,1.0_iwp,storkm_pp(1:ntot,1:ntot,iel),ntot,    &
                  pmul_pp(1:ntot,iel),1,0.0_iwp,utemp_pp(1:ntot,iel),1)
     END DO
     !$OMP END DO
     !$OMP END PARALLEL

     CASE(4) ! DSYMV with all element stiffness matrices symmetric

     !$OMP PARALLEL
     !$OMP DO PRIVATE(iel)     
     DO iel=1,nels_pp
       CALL dspmv('L',ntot,1.0,vstorkm_pp(:,iel),pmul_pp(:,iel),1,0.0, &
                   utemp_pp(:,iel),1)    
     END DO
     !$OMP END DO
     !$OMP END PARALLEL
      
     ! PRINT *, "Test case ",testcase, " not implemented. Program aborted"
     ! CALL SHUTDOWN()        

     CASE(5) ! Minimal storage of KM in elements loop

        PRINT *, "Test case ",testcase, " not implemented. Program aborted"
        CALL SHUTDOWN()        

     CASE(6) ! Matrix-matrix multiply

     km = storkm_pp(1:ntot,1:ntot,1)
     CALL dgemm('N','N',ntot,nels_pp,ntot,one,km,ntot,pmul_pp,ntot,zero,      &
                 utemp_pp,ntot)

     CASE(7) ! Symmetric matrix-matrix multiply
   
     km = storkm_pp(1:ntot,1:ntot,1) 
     CALL dsymm('L','l',ntot,nels_pp,one,km,ntot,pmul_pp,ntot,one,utemp_pp,ntot)


     CASE DEFAULT

        PRINT *, "Test case ",testcase, " not recognised. Program aborted"
        CALL SHUTDOWN()        

   END SELECT

   timest(7)=timest(7)+MPI_Wtime()-timest(6)
 
   timest(8)=MPI_Wtime()

   CALL scatter(u_pp,utemp_pp(1:ntot,:))

   timest(9)=timest(9)+MPI_Wtime()-timest(8)

!-------------------------------------------------------------------------
!-------------------------- pcg equation solution ------------------------
!-------------------------------------------------------------------------

   timest(10)= MPI_Wtime()

   IF(testcase==3) THEN 
    
     ! add in some more openmp

     up       = DOT_PRODUCT_P(r_pp,d_pp)
     dot_temp = DOT_PRODUCT_P(p_pp,u_pp)
     alpha    = up/dot_temp

     !$OMP PARALLEL
     !$OMP WORKSHARE

     xnew_pp  = x_pp+p_pp*alpha
     r_pp     = r_pp-u_pp*alpha
     d_pp     = diag_precon_pp*r_pp

     !$OMP END WORKSHARE
     !$OMP END PARALLEL

     dot_temp = DOT_PRODUCT_P(r_pp,d_pp)
     beta     = dot_temp/up


     !$OMP PARALLEL
     !$OMP WORKSHARE

     p_pp     = d_pp+p_pp*beta

     !$OMP END WORKSHARE
     !$OMP END PARALLEL

   ELSE

     up       = DOT_PRODUCT_P(r_pp,d_pp)
     dot_temp = DOT_PRODUCT_P(p_pp,u_pp)
     alpha    = up/dot_temp
     xnew_pp  = x_pp+p_pp*alpha
     r_pp     = r_pp-u_pp*alpha
     d_pp     = diag_precon_pp*r_pp
     dot_temp = DOT_PRODUCT_P(r_pp,d_pp)
     beta     = dot_temp/up
     p_pp     = d_pp+p_pp*beta

   END IF

   timest(11)= timest(11)+MPI_Wtime()-timest(10)

   CALL checon_par(xnew_pp,tol,converged,x_pp)    
   IF(converged.OR.iters==limit)EXIT
 END DO iterations

  timest(12)=MPI_Wtime()-timest(3)  

 IF(sym_storkm) THEN
   DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,vstorkm_pp,pmul_pp,     &
              utemp_pp) 

 ELSE
   DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,storkm_pp,pmul_pp,      &
              utemp_pp) 
 END IF

!--------------- recover stresses at centroidal gauss point --------------
 IF(numpe==1) THEN
   WRITE(11,'(A,I6)')    " Iterations to convergence was:         ",iters
   WRITE(11,'(A,E12.4)') " The central nodal displacement:  ",                &
                          xnew_pp(1)
   WRITE(11,'(/A)')" The Centroid point stresses for element 1 are"
 END IF 
 ALLOCATE(eld_pp(ntot,nels_pp)); eld_pp=zero; points=zero; nip=1; iel=1
 CALL gather(xnew_pp(1:),eld_pp); DEALLOCATE(xnew_pp)
 gauss_pts_2: DO i=1,nip
   CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
   CALL invert(jac); deriv=MATMUL(jac,der); CALL beemat(bee,deriv)
   eps=MATMUL(bee,eld_pp(:,iel)); sigma=MATMUL(dee,eps)
   IF(numpe==1.AND.i==1) THEN
     WRITE(11,'(A,I5)')" Point ",i ; WRITE(11,'(6E12.4)') sigma
   END IF
 END DO gauss_pts_2; DEALLOCATE(g_coord_pp)
!------------------------ write out displacements ------------------------

 CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)

 IF(numpe==1) THEN

   IF(io_binary) THEN
     OPEN(12,file=argv(1:nlen)//".bin.ensi.DISPL-000001",status='replace',&
          form='unformatted',access='stream')
     cbuffer="Alya Ensight Gold --- Vector per-node variable file"
     WRITE(12) cbuffer
     cbuffer="part"
     WRITE(12) cbuffer
     WRITE(12) int(1,kind=c_int)
     cbuffer="coordinates"
     WRITE(12) cbuffer
   ELSE
     WRITE(ch,'(I6.6)') numpe
     OPEN(12,file=argv(1:nlen)//".ensi.DISPL-"//ch,status='replace',      &
          action='write')
     WRITE(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
     WRITE(12,'(A/A/A)') "part", "     1","coordinates"
   END IF

 END IF

 ALLOCATE(disp_pp(nodes_pp*ndim)); disp_pp=zero
 timest(21)=MPI_Wtime()
 CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,          &
                    node_start,node_end,eld_pp,disp_pp,1)
 timest(22)=MPI_Wtime()-timest(21)
 timest(23)=MPI_Wtime()

 DEALLOCATE(g_num_pp)
 ALLOCATE(temp(nodes_pp)) ; temp = zero

 IF(io_binary) THEN
   DO i=1,ndim ; temp=zero
     DO j=1,nodes_pp; k=i+(ndim*(j-1)); temp(j)=disp_pp(k); END DO
     CALL dismsh_ensi_pb2(12,1,nodes_pp,npes,numpe,1,temp)
   END DO ; IF(numpe==1) CLOSE(12)
 ELSE
   DO i=1,ndim ; temp=zero
     DO j=1,nodes_pp; k=i+(ndim*(j-1)); temp(j)=disp_pp(k); END DO
     CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp)
   END DO ; IF(numpe==1) CLOSE(12)
 END IF

 timest(24)=MPI_Wtime()-timest(23)
 
!--------------------------- gigaflops in pcg -----------------------------
 flop   = zero ; gigaflops   = zero
 flop_1 = zero ; gigaflops_1 = zero
 flop_2 = zero ; gigaflops_2 = zero

! elements_3 loop
  flop_1      = real(iters)*(real(nels)*(2._iwp*real(ntot)*real(ntot)))  
  gigaflops_1 = flop_1/1.E+09_iwp

! pcg process
  flop_2      = real(iters)*((2._iwp*real(neq))+(2._iwp*real(neq))+       &
                (2._iwp*real(neq))+(2._iwp*real(neq))+real(neq)+          &
                (2._iwp*real(neq))+(2._iwp*real(neq)))
  gigaflops_2 = flop_2/1.E+09_iwp     ! solution

! total flops

  gigaflops   = gigaflops_1 + gigaflops_2

! maximum and minimum timings

  CALL MPI_ALLREDUCE(timest,timest_min,24,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ier)

  CALL MPI_ALLREDUCE(timest,timest_max,24,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ier)

!-------------------------------------------------------------------------------
 IF(numpe==1)THEN
 
   WRITE(11,'(/2A/)')    " WALL CLOCK TIME                               MAX", &
                         "       MIN"
   WRITE(11,'(A,2F10.4)')" Time to read input:                    ",      &
                          timest_max(2),timest_min(2) 
   WRITE(11,'(A,2F10.4)')" Time to find g:                        ",      &
                          timest_max(16),timest_min(16)  
   WRITE(11,'(A,2F10.4)')" Time to compute stiffness:             ",      &
                          timest_max(14),timest_min(14)  
   WRITE(11,'(A,2F10.4)')" Time to compute preconditioner:        ",      &
                          timest_max(18),timest_min(18) 
   WRITE(11,'(A,2F10.4)')" Time to get starting r:                ",      &
                          timest_max(20),timest_min(20) 
   WRITE(11,'(A,2F10.4)')" Time in gather:                        ",      &
                          timest_max(5),timest_min(5)  
   WRITE(11,'(A,2F10.4)')" Time in elements loop:                 ",      &
                          timest_max(7),timest_min(7)  
   WRITE(11,'(A,2F10.4)')" Time in scatter:                       ",      &
                          timest_max(9),timest_min(9)  
   WRITE(11,'(A,2F10.4)')" Time in PCG process:                   ",      &
                          timest_max(11),timest_min(11)  
   WRITE(11,'(A,2F10.4)')" Time in scatter nodes:                 ",      &
                          timest_max(22),timest_min(22)  
   WRITE(11,'(A,2F10.4/)')" Time to write displacements:           ",     &
                          timest_max(24),timest_min(24)  

   WRITE(11,'(A,2F10.4)')" Total time to solve equations:         ",      &
                          timest_max(12),timest_min(12)  
   WRITE(11,'(A,2F10.4)')" Total time for the analysis:           ",      &
                          MPI_WTIME()-timest(1)  

   WRITE(11,'(/A/)')     " PERFORMANCE IN SOLVER "

   WRITE(11,'(A,F10.4)') " GFLOPS in elements loop:                         ", &
                          gigaflops_1  
   WRITE(11,'(A,F10.4)') " GFLOPS in pcg process:                           ", &
                          gigaflops_2  
   WRITE(11,'(A,F10.4/)')" Total GFLOPS in solver:                          ", &
                          gigaflops  

   WRITE(11,'(A,2F10.4)')" GFLOPS per second in elements loop:    ",           &
                          gigaflops_1/timest_max(7),gigaflops_1/timest_min(7)  
   WRITE(11,'(A,2F10.4)')" GFLOPS per second in PCG process:      ",           &
                          gigaflops_2/timest_max(11),gigaflops_2/timest_min(11) 
   WRITE(11,'(A,2F10.4)')" GFLOPS per second in solver:           ",           &
                          gigaflops/timest_max(12),gigaflops/timest_min(12)  

 END IF

 CALL SHUTDOWN() 

CONTAINS

 SUBROUTINE write_summary(fnum,neq,nn,npes,nr,numpe,section,testcase,timest)

   IMPLICIT none

   INTEGER,INTENT(IN) :: fnum,neq,nn,npes,nr,numpe,section,testcase
   REAL,INTENT(IN)    :: timest(:)
   CHARACTER(LEN=8)   :: nthreads
   CHARACTER(LEN=60)  :: date
   CHARACTER(LEN=60)  :: hostname

!--------------------------------------------------------------------------
! 1. Retrieve information about the environment
!
!    Run the following command before running the job to get the date and
!    time the job was executed or submitted
!
!    $export MYDATE=$(date)
!--------------------------------------------------------------------------

   CALL get_environment_variable("OMP_NUM_THREADS",nthreads)
   CALL get_environment_variable("MYDATE",date)
   CALL get_environment_variable("HOSTNAME",hostname)

   IF(section==1) THEN 
     
     WRITE(fnum,'(/2A)') " ",TRIM(date)

     WRITE(fnum,'(/A)') " BASIC JOB INFORMATION"

     IF(testcase==1) WRITE(fnum,'(/A/)') " Test case 1: MATMUL"
     IF(testcase==2) WRITE(fnum,'(/A/)') " Test case 2: DGEMV"
     IF(testcase==3) WRITE(fnum,'(/A/)') " Test case 3: DGEMV + OPENMP"
     IF(testcase==4) WRITE(fnum,'(/A/)') " Test case 4: DSYMV + OPENMP"
     IF(testcase==5) WRITE(fnum,'(/A/)') " Test case 5: MINIMAL KM + OPENMP"
     IF(testcase==6) WRITE(fnum,'(/A/)') " Test case 6: DGEMM"
     IF(testcase==7) WRITE(fnum,'(/A/)') " Test case 7: DSYMM"

     WRITE(fnum,'(2A)')   " This job run on host ", TRIM(hostname)

     WRITE(fnum,'(/A,I7)')" Number of MPI processes:              ",npes
     WRITE(fnum,'(2A)')   " Number of OMP threads:                      ",    &
                            TRIM(nthreads)
     WRITE(fnum,'(A,I7)') " Number of elements:                   ",nels
     WRITE(fnum,'(A,I7)') " Number of nodes:                      ",nn
     WRITE(fnum,'(A,I7)') " Number of equations:                  ",neq

   ELSE

     PRINT *, " Error in subroutine write_summary"

   END IF

 RETURN
 END SUBROUTINE

END PROGRAM xx16
