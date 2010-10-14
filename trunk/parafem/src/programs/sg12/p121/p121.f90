PROGRAM p121         
!------------------------------------------------------------------------------ 
!      Program 12.1 three dimensional analysis of an elastic solid
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

  INTEGER,PARAMETER     :: nodof=3,ndim=3,nst=6
  INTEGER               :: loaded_nodes,iel,i,iters,limit
  INTEGER               :: nn,nr,nip,nod,nels,ndof,npes_pp
  INTEGER               :: node_end,node_start,nodes_pp
  INTEGER               :: argc,iargc,meshgen
  REAL(iwp)             :: e,v,det,tol,up,alpha,beta,tload
  REAL(iwp),PARAMETER   :: zero = 0.0_iwp
  CHARACTER(LEN=15)     :: element
  CHARACTER(LEN=50)     :: program_name='p121'
  CHARACTER(LEN=50)     :: fname,job_name,label
  LOGICAL               :: converged = .false.

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE :: points(:,:),dee(:,:),weights(:),val(:,:),disp_pp(:)
  REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:),jac(:,:),der(:,:),deriv(:,:)
  REAL(iwp),ALLOCATABLE :: bee(:,:),storkm_pp(:,:,:),eld(:),eps(:),sigma(:)
  REAL(iwp),ALLOCATABLE :: diag_precon_pp(:),p_pp(:),r_pp(:),x_pp(:),xnew_pp(:)
  REAL(iwp),ALLOCATABLE :: u_pp(:),pmul_pp(:,:),utemp_pp(:,:),d_pp(:),timest(:)
  REAL(iwp),ALLOCATABLE :: diag_precon_tmp(:,:),eld_pp(:,:),tensor_pp(:,:,:)
  INTEGER,  ALLOCATABLE :: rest(:,:),g_num_pp(:,:),g_g_pp(:,:),node(:)

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

  CALL read_p121(job_name,numpe,e,element,limit,loaded_nodes,meshgen,nels,nip,&
                 nn,nod,nr,tol,v)

  CALL calc_nels_pp(nels)

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

  ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim),                      &
           der(ndim,nod),deriv(ndim,nod),eld_pp(ntot,nels_pp),bee(nst,ntot),  &
           storkm_pp(ntot,ntot,nels_pp),eld(ntot),eps(nst),sigma(nst),        &
           pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp),                      &
           weights(nip),diag_precon_tmp(ntot,nels_pp),g_g_pp(ntot,nels_pp))

  timest(2) = elap_time()

!------------------------------------------------------------------------------
! 5. Loop the elements to find the steering array and the number of equations
!    to solve.
!------------------------------------------------------------------------------

  CALL rearrange(rest)
  
  g_g_pp = 0

  elements_1: DO iel = 1, nels_pp
    CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
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

  ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),            &
           u_pp(neq_pp),d_pp(neq_pp),diag_precon_pp(neq_pp))

  p_pp    = zero  ;  r_pp = zero  ;  x_pp = zero
  xnew_pp = zero  ;  u_pp = zero  ;  d_pp = zero  ; diag_precon_pp = zero

  timest(5) = elap_time()

!------------------------------------------------------------------------------
! 8. Element stiffness integration and storage
!------------------------------------------------------------------------------

  CALL deemat(e,v,dee)
  CALL sample(element,points,weights)
 
  storkm_pp       = zero
  diag_precon_tmp = zero
 
  elements_3: DO iel=1,nels_pp
    gauss_pts_1: DO i=1,nip
      CALL shape_der(der,points,i)
      jac   = MATMUL(der,g_coord_pp(:,:,iel))
      det   = determinant(jac)
      CALL invert(jac)
      deriv = MATMUL(jac,der)
      CALL beemat(deriv,bee)
      storkm_pp(:,:,iel)   = storkm_pp(:,:,iel) +                             &
                             MATMUL(MATMUL(TRANSPOSE(bee),dee),bee) *         &
                             det*weights(i)   
    END DO gauss_pts_1
  END DO elements_3
  
  timest(6) = elap_time()

!------------------------------------------------------------------------------
! 9. Build the diagonal preconditioner
!------------------------------------------------------------------------------

  diag_precon_tmp = zero
 
  elements_4: DO iel = 1,nels_pp 
    DO i = 1,ndof
      diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel) + storkm_pp(i,i,iel)
    END DO
  END DO elements_4

  CALL scatter(diag_precon_pp,diag_precon_tmp)

  DEALLOCATE(diag_precon_tmp)

  timest(7) = elap_time()

!------------------------------------------------------------------------------
! 10. Get starting r_pp
!------------------------------------------------------------------------------
 
  IF(loaded_nodes > 0) THEN


    ALLOCATE(node(loaded_nodes))
    ALLOCATE(val(ndim,loaded_nodes))
    
    val    = zero
    node   = 0

    CALL read_loads(job_name,numpe,node,val)
    CALL load(g_g_pp,g_num_pp,node,val,r_pp(1:))

    tload = SUM_P(r_pp(1:))

    DEALLOCATE(node)
    DEALLOCATE(val)

  END IF
  
  timest(8) = elap_time()

  diag_precon_pp = 1._iwp/diag_precon_pp
  d_pp           = diag_precon_pp*r_pp
  p_pp           = d_pp

!------------------------------------------------------------------------------
! 11. Preconditioned conjugate gradient iterations
!------------------------------------------------------------------------------

  iters = 0

  iterations: DO 
    iters   = iters + 1
    u_pp    = zero
    pmul_pp = zero

    CALL gather(p_pp,pmul_pp)
    elements_5: DO iel=1,nels_pp
      utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
    END DO elements_5
    CALL scatter(u_pp,utemp_pp)

    up      = DOT_PRODUCT_P(r_pp,d_pp)
    alpha   = up/DOT_PRODUCT_P(p_pp,u_pp)
    xnew_pp = x_pp + p_pp*alpha
    r_pp    = r_pp - u_pp*alpha
    d_pp    = diag_precon_pp*r_pp
    beta    = DOT_PRODUCT_P(r_pp,d_pp)/up
    p_pp    = d_pp + p_pp*beta  

    CALL checon_par(xnew_pp,tol,converged,x_pp)    
    IF(converged.OR.iters==limit)EXIT
  END DO iterations

  timest(9) = elap_time()

!------------------------------------------------------------------------------
! 12. Recover stresses at centroidal gauss point
!------------------------------------------------------------------------------

  ALLOCATE(tensor_pp(nst,nip,nels_pp))
  tensor_pp = zero
  nip       = 1
  eld_pp    = zero
  
  CALL gather(xnew_pp,eld_pp)

  elements_6: DO iel=1,nels_pp
    gauss_pts_2: DO i=1,nip 
      CALL shape_der(der,points,i)
      jac                = MATMUL(der,g_coord_pp(:,:,iel))
      CALL invert(jac)
      deriv              = MATMUL(jac,der)
      CALL beemat(deriv,bee)
      eps                = MATMUL(bee,eld_pp(:,1))
      sigma              = MATMUL(dee,eps)
      tensor_pp(:,i,iel) = sigma
    END DO gauss_pts_2
  END DO elements_6
  
  timest(10) = elap_time()

!------------------------------------------------------------------------------
! 13. Output results
!------------------------------------------------------------------------------

  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
  ALLOCATE(disp_pp(nodes_pp*ndim))
    
  disp_pp = zero

  IF(numpe==1) THEN
    fname   = job_name(1:INDEX(job_name, " ")-1)//".dis"
    OPEN(24, file=fname, status='replace', action='write')
  END IF

  label     = "*DISPLACEMENT"
  utemp_pp  = zero
  CALL gather(xnew_pp(1:),utemp_pp)
  CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,              &
                     node_start,node_end,utemp_pp,disp_pp,1)
  CALL write_nodal_variable(label,24,1,nodes_pp,npes,numpe,ndim,disp_pp)

  IF(numpe==1) CLOSE(24)

  timest(11) = elap_time()
   
!------------------------------------------------------------------------------
! 14. Output performance data
!------------------------------------------------------------------------------

  CALL WRITE_P121(iters,job_name,neq,nn,npes,nr,numpe,timest,tload)
 
  CALL shutdown() 
 
END PROGRAM p121
