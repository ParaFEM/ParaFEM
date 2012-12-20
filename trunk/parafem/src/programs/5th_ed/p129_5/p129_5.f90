PROGRAM p129       
!------------------------------------------------------------------------------
!      Program 12.9 forced vibration of a 3d elastic solid
!      Lumped or consistent mass
!      Implicit integration by theta method : parallel version
!------------------------------------------------------------------------------

  USE precision     ; USE global_variables ; USE mp_interface
  USE input         ; USE output           ; USE loading
  USE timing        ; USE maths            ; USE gather_scatter
  USE partition     ; USE elements         ; USE steering
  USE pcg

  IMPLICIT NONE
 
!------------------------------------------------------------------------------ 
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

 ! neq,ntot are now global variables - not declared 

 INTEGER               :: nn,nr,nip,nodof=3,nod,nst=6
 INTEGER               :: i,j,k,iel,ndim=3,nstep,npri,iters,limit
 INTEGER               :: ndof,nels,npes_pp
 INTEGER               :: node_end,node_start,nodes_pp,loaded_nodes
 INTEGER               :: argc,iargc ! command line arguments 
 INTEGER               :: meshgen    ! mesh type 1 = S&G ; 2 = Abaqus
 INTEGER               :: partitioner! type 1 = S&G ; type 2 = external
 REAL(iwp)             :: e,v,det,rho,alpha1,beta1,omega,theta,period
 REAL(iwp)             :: pi,dtim,volume,c1,c2,c3,c4,real_time,tol,big,up
 REAL(iwp)             :: alpha,beta,tload     
 REAL(iwp),PARAMETER   :: zero = 0.0_iwp    
 CHARACTER(LEN=15)     :: element
 CHARACTER(LEN=50)     :: program_name='p129g' 
 CHARACTER(LEN=50)     :: fname
 CHARACTER(LEN=50)     :: job_name
 CHARACTER(LEN=50)     :: label
 CHARACTER(LEN=6)      :: step
 LOGICAL               :: consistent=.TRUE.,converged
 
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

 REAL(iwp),ALLOCATABLE :: loads_pp(:),points(:,:),dee(:,:),fun(:),jac(:,:)
 REAL(iwp),ALLOCATABLE :: der(:,:),deriv(:,:),weights(:),bee(:,:),km(:,:)
 REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:),x1_pp(:),d1x1_pp(:),d2x1_pp(:)
 REAL(iwp),ALLOCATABLE :: emm(:,:),ecm(:,:),x0_pp(:),d1x0_pp(:),d2x0_pp(:)
 REAL(iwp),ALLOCATABLE :: store_km_pp(:,:,:),vu_pp(:),store_mm_pp(:,:,:)
 REAL(iwp),ALLOCATABLE :: u_pp(:),p_pp(:),d_pp(:),x_pp(:),xnew_pp(:)
 REAL(iwp),ALLOCATABLE :: pmul_pp(:,:),utemp_pp(:,:),diag_precon_pp(:)
 REAL(iwp),ALLOCATABLE :: diag_precon_tmp(:,:),temp_pp(:,:,:),disp_pp(:)
 REAL(iwp),ALLOCATABLE :: loadNodeValue(:,:),fext_pp(:),ans(:,:),val(:,:)
 REAL(iwp),ALLOCATABLE :: timest(:)
 INTEGER,ALLOCATABLE   :: rest(:,:),g(:),g_num_pp(:,:),g_g_pp(:,:)
 INTEGER,ALLOCATABLE   :: loadNodeNum(:),ttliters(:),node(:)       
 
!------------------------------------------------------------------------------
! 3. Read job_name from the command line. 
!    Read control data, mesh data, boundary and loading conditions. 
!------------------------------------------------------------------------------

  ALLOCATE(timest(20))
  timest    = zero 
  timest(1) = elap_time()

  CALL find_pe_procs(numpe,npes)
  argc = iargc()
  IF (argc /= 1) CALL job_name_error(numpe,program_name)
  CALL GETARG(1, job_name)

  CALL read_p129(job_name,numpe,alpha1,beta1,e,element,limit,loaded_nodes,    &
                 meshgen,nels,nip,nn,nod,npri,nr,nstep,omega,partitioner,     &
                 rho,theta,tol,v)
   
  CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)

  ndof = nod*nodof
  ntot = ndof
 
  ALLOCATE(g_num_pp(nod, nels_pp)) 
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp)) 
  ALLOCATE(rest(nr,nodof+1)) 
 
  g_num_pp  = 0
  g_coord_pp= zero
  rest      = 0

  CALL read_g_num_pp(job_name,iel_start,nels,nn,numpe,g_num_pp)
  IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
  CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
  CALL read_rest(job_name,numpe,rest)
   
!------------------------------------------------------------------------------
! 4. Allocate dynamic arrays used in main program
!------------------------------------------------------------------------------
 
 ALLOCATE(points(nip,ndim),g(ntot),fun(nod),dee(nst,nst),jac(ndim,ndim),      &
          weights(nip),der(ndim,nod),deriv(ndim,nod),bee(nst,ntot),           &
          g_g_pp(ntot,nels_pp),emm(ntot,ntot),ecm(ntot,ntot),                 &
          store_km_pp(ntot,ntot,nels_pp),utemp_pp(ntot,nels_pp),              &
          pmul_pp(ntot,nels_pp),store_mm_pp(ntot,ntot,nels_pp),               &
          temp_pp(ntot,ntot,nels_pp),diag_precon_tmp(ntot,nels_pp),           &
          ans(3,nstep),ttliters(nstep))
 
 timest(2) = elap_time()

!------------------------------------------------------------------------------
! 5. Loop the elements to find the steering array and the number of equations
!    to solve.
!------------------------------------------------------------------------------

  CALL rearrange(rest)
  
  g_g_pp = 0

  elements_1: DO iel = 1, nels_pp
    CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest)
  END DO elements_1

  neq = 0
  
  elements_2: DO iel = 1, nels_pp  
    i = MAXVAL(g_g_pp(:,iel))
    IF(i > neq) neq = i
  END DO elements_2  

  neq = MAX_INTEGER_P(neq)
 
  timest(3) = elap_time()
  
!------------------------------------------------------------------------------
! 6. Create interprocessor communication tables
!------------------------------------------------------------------------------

  CALL calc_neq_pp
  CALL calc_npes_pp(npes,npes_pp)
  CALL make_ggl(npes_pp,npes,g_g_pp)
 
  timest(4) = elap_time()
  
!------------------------------------------------------------------------------
! 7. Allocate arrays dimensioned by neq_pp 
!------------------------------------------------------------------------------

 ALLOCATE(x0_pp(neq_pp),d1x0_pp(neq_pp),x1_pp(neq_pp),vu_pp(neq_pp),          &
          diag_precon_pp(neq_pp),u_pp(neq_pp),d2x0_pp(neq_pp),                &
          loads_pp(neq_pp),d1x1_pp(neq_pp),d2x1_pp(neq_pp),d_pp(neq_pp),      &
          p_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),fext_pp(neq_pp))

  x0_pp          = zero  ;  d1x0_pp        = zero  ;  x1_pp    = zero 
  vu_pp          = zero  ;  diag_precon_pp = zero  ;  u_pp     = zero
  d2x0_pp        = zero  ;  loads_pp       = zero  ;  d1x1_pp  = zero  
  d2x1_pp        = zero  ;  d_pp           = zero  ;  p_pp     = zero
  x_pp           = zero  ;  xnew_pp        = zero  ;  fext_pp  = zero 

  timest(5) = elap_time()

!------------------------------------------------------------------------------
! 8. Element stiffness integration and storage
!------------------------------------------------------------------------------ 
 
 CALL deemat(e,v,dee)
 CALL sample(element,points,weights)
 
 store_km_pp     = zero
 store_mm_pp     = zero
 diag_precon_tmp = zero
 pi              = ACOS(-1._iwp)
 period          = 2._iwp * pi/omega
 dtim            = period/20._iwp
 c1              = (1._iwp - theta) * dtim
 c2              = beta1 - c1
 c3              = alpha1 + 1._iwp/(theta * dtim)
 c4              = beta1 + theta * dtim

 elements_3: DO iel=1,nels_pp

   volume = zero
   emm    = zero
   ecm    = zero

   gauss_points_1: DO i=1,nip     
     CALL shape_der(der,points,i)
     jac   = MATMUL(der,g_coord_pp(:,:,iel))
     det   = determinant(jac)
     CALL invert(jac)
     deriv = matmul(jac,der)
     CALL beemat(deriv,bee)
     store_km_pp(:,:,iel) = store_km_pp(:,:,iel) +                            & 
                            MATMUL(MATMUL(TRANSPOSE(bee),dee),bee) *          &
                            det*weights(i)
     volume               = volume + det*weights(i)
     CALL shape_fun(fun,points,i)
     IF(consistent)THEN
       CALL ecmat(ecm,fun,ntot,nodof)
       ecm = ecm*det*weights(i)*rho
       emm = emm+ecm
     END IF
   END DO gauss_points_1   
   
   IF(.NOT.consistent)THEN
     DO i=1,ntot
       emm(i,i)=volume*rho/13._iwp
     END DO
     DO i=1,19,6
       emm(i,i)=emm(4,4)*.125_iwp
     END DO
     DO i=2,20,6
       emm(i,i)=emm(4,4)*.125_iwp
     END DO
     DO i=3,21,6
       emm(i,i)=emm(4,4)*.125_iwp
     END DO
     DO i=37,55,6
       emm(i,i)=emm(4,4)*.125_iwp
     END DO
     DO i=38,56,6
       emm(i,i)=emm(4,4)*.125_iwp
     END DO
     DO i=39,57,6
       emm(i,i)=emm(4,4)*.125_iwp
     END DO
   END IF
 
   store_mm_pp(:,:,iel)=emm
 
 END DO elements_3

 timest(6) = elap_time()
 
!------------------------------------------------------------------------------
! 9. Build the diagonal preconditioner
!------------------------------------------------------------------------------ 

 diag_precon_tmp = zero

 elements_4: DO iel=1,nels_pp
   DO k=1,ntot
     diag_precon_tmp(k,iel) = diag_precon_tmp(k,iel)    +                     &
                              store_mm_pp(k,k,iel) * c3 +                     &
                              store_km_pp(k,k,iel) * c4
   END DO
 END DO elements_4
 
 CALL scatter(diag_precon_pp,diag_precon_tmp)
 
 diag_precon_pp=1._iwp/diag_precon_pp
 
 DEALLOCATE(diag_precon_tmp)

 timest(7) = elap_time()

!-------------------------------------------------------------------------------
! 10. Read in applied forces and assign to equations
!-------------------------------------------------------------------------------
 
  IF(loaded_nodes > 0) THEN

    ALLOCATE(node(loaded_nodes))
    ALLOCATE(val(ndim,loaded_nodes))
    
    val    = zero
    node   = 0

    CALL read_loads(job_name,numpe,node,val)
    CALL load(g_g_pp,g_num_pp,node,val,fext_pp(1:))

    tload = SUM_P(fext_pp(1:))

    DEALLOCATE(node)
    DEALLOCATE(val)

  END IF

  timest(8) = elap_time()
  
!------------------------------------------------------------------------------
! 8. Allocate disp_pp array and open file to write displacement output
!------------------------------------------------------------------------------

  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
  ALLOCATE(disp_pp(nodes_pp*ndim))
    
  disp_pp = zero 
  
  IF(numpe==1) THEN
    fname   = job_name(1:INDEX(job_name, " ")-1)//".dis"
    OPEN(24, file=fname, status='replace', action='write')
    label   = "*DISPLACEMENT"  
  END IF
  
!------------------------------------------------------------------------------
! 11. Initial conditions
!------------------------------------------------------------------------------

 x0_pp     = zero
 d1x0_pp   = zero
 d2x0_pp   = zero
 real_time = zero
 ans       = zero
 ttliters  = 0

!------------------------------------------------------------------------------
! 12. Time stepping loop
!------------------------------------------------------------------------------ 

 timesteps: DO j=1,nstep
   
   timest(9) = elap_time()
   real_time = real_time + dtim
   loads_pp  = zero
   u_pp      = zero
   vu_pp     = zero
   
   elements_5: DO iel=1,nels_pp    ! gather for rhs multiply
     temp_pp(:,:,iel)=store_km_pp(:,:,iel)*c2+store_mm_pp(:,:,iel)*c3
   END DO elements_5
   
   CALL gather(x0_pp,pmul_pp)
   
   DO iel=1,nels_pp
     utemp_pp(:,iel)=MATMUL(temp_pp(:,:,iel),pmul_pp(:,iel))
   END DO
   
   CALL scatter(u_pp,utemp_pp)

!------------------------------------------------------------------------------
! 13. Velocity
!------------------------------------------------------------------------------

   temp_pp = store_mm_pp/theta
   CALL gather(d1x0_pp,pmul_pp)
   
   DO iel=1,nels_pp
     utemp_pp(:,iel)=MATMUL(temp_pp(:,:,iel),pmul_pp(:,iel))
   END DO
   
   CALL scatter(vu_pp,utemp_pp)   ! doesn't add to last u_pp
  
   loads_pp = fext_pp * (theta*dtim*cos(omega*real_time)+c1*                  &
                         cos(omega*(real_time-dtim))) 

   loads_pp = u_pp + vu_pp + loads_pp

!------------------------------------------------------------------------------
! 14. Solve simultaneous equations by PCG
!------------------------------------------------------------------------------

   d_pp  = diag_precon_pp*loads_pp
   p_pp  = d_pp
   x_pp  = zero
   iters = 0
   
   iterations: DO 

     iters   = iters + 1
     u_pp    = zero
     vu_pp   = zero
     temp_pp = store_mm_pp*c3+store_km_pp*c4
     
     CALL gather(p_pp,pmul_pp)
     
     elements_6: DO iel=1,nels_pp
       utemp_pp(:,iel)=MATMUL(temp_pp(:,:,iel),pmul_pp(:,iel))
     END DO elements_6
     
     CALL scatter(u_pp,utemp_pp)
     
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
   
   x1_pp     = xnew_pp 
   d1x1_pp   = (x1_pp-x0_pp)/(theta*dtim)-d1x0_pp*(1._iwp-theta)/theta
   d2x1_pp   = (d1x1_pp-d1x0_pp)/(theta*dtim)-d2x0_pp*(1._iwp-theta)/theta
   x0_pp     = x1_pp
   d1x0_pp   = d1x1_pp
   d2x0_pp   = d2x1_pp
   timest(10)= timest(10) + (elap_time() - timest(9))
   timest(9) = elap_time()
   
   utemp_pp  = zero
   CALL gather(x1_pp(1:),utemp_pp)
   CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,              &
                      node_start,node_end,utemp_pp,disp_pp,1)
   CALL write_nodal_variable(label,24,j,nodes_pp,npes,numpe,ndim,disp_pp)      
   
   timest(11)  = timest(11) + (elap_time() - timest(9))
   ans(1,j)    = real_time
   ans(2,j)    = cos(omega*real_time)
   ans(3,j)    = x1_pp(neq_pp)
   ttliters(j) = iters
 
 END DO timesteps
 
 timest(12) = elap_time()
 
!------------------------------------------------------------------------------
! 15. Output performance data and end program
!------------------------------------------------------------------------------ 

 CALL WRITE_P129(ans,job_name,neq,nn,npes,nr,numpe,timest,tload,ttliters)
 
 CLOSE(24)

 CALL shutdown() 

END PROGRAM p129
