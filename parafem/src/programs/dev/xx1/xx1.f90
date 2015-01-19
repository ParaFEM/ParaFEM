PROGRAM xx1         
!------------------------------------------------------------------------------ 
! Program XX.1 Three dimensional analysis of an elastic solid using load 
!              control or displacement control. Compare with P121.
!
!              Support for Abaqus UMAT
!              Some binary I/O. 
!              Support for multiple material types.
!              Vector storage for symmetric element stiffness matrix.
!              Global strain energy computation.
!------------------------------------------------------------------------------ 
                                 
  USE precision     ; USE global_variables ; USE mp_interface
  USE input         ; USE output           ; USE loading
  USE timing        ; USE maths            ; USE gather_scatter
  USE partition     ; USE elements         ; USE steering        ; USE pcg
  USE new_library

  USE, INTRINSIC :: ISO_C_BINDING ! to output C binary file
  
  IMPLICIT NONE

!------------------------------------------------------------------------------ 
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

 ! neq,ntot are now global variables - not declared

  INTEGER,PARAMETER     :: nodof=3,ndim=3,nst=6
  INTEGER               :: loaded_nodes,fixed_freedoms,iel,i,j,k,l,idx1,idx2
  INTEGER               :: iters,limit,nn,nr,nip,nod,nels,ndof,npes_pp
  INTEGER               :: node_end,node_start,nodes_pp
  INTEGER               :: argc,iargc,meshgen,partitioner=1
  INTEGER               :: fixed_freedoms_pp,fixed_freedoms_start
  INTEGER               :: nprops,np_types
  REAL(iwp)             :: e,v,det,tol,up,alpha,beta,tload
  REAL(iwp)             :: ewt_pp ! total strain energy on local core
  REAL(iwp)             :: ewt    ! total strain energy for whole domain
  REAL(iwp)             :: ewe    ! strain energy for element
  REAL(iwp),PARAMETER   :: zero=0.0_iwp
  REAL(iwp),PARAMETER   :: penalty=1.0e20_iwp
  CHARACTER(LEN=15)     :: element
  CHARACTER(LEN=50)     :: program_name='xx1'
  CHARACTER(LEN=50)     :: job_name,label
  LOGICAL               :: converged=.false.
  LOGICAL               :: sym_storkm=.true.
! LOGICAL               :: sym_storkm=.false.
! LOGICAL               :: io_binary=.false.
  LOGICAL               :: io_binary=.true.

!------------------------------------------------------------------------------ 
! 1a. Declare variables used in the UMAT
!
!     ntens == nst
!     noel  == iel
!
!------------------------------------------------------------------------------

! INTEGER,PARAMETER     :: ndi=3,nprops=2
  INTEGER               :: nshr,nstatv,npt,layer,kspt,kstep,kinc
  REAL(iwp)             :: sse,spd,scd,rpl,dtime,temp,dtemp,pnewdt,celent
  CHARACTER(LEN=80)     :: cmname
  CHARACTER(LEN=80)     :: cbuffer
  CHARACTER(LEN=50)     :: fname   ! should be obsolete

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE :: points(:,:),dee(:,:),weights(:),val(:,:),disp_pp(:)
  REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:),jac(:,:),der(:,:),deriv(:,:)
  REAL(iwp),ALLOCATABLE :: bee(:,:),storkm_pp(:,:,:),eld(:,:),eps(:),sigma(:)
  REAL(iwp),ALLOCATABLE :: vstorkm_pp(:,:),km(:,:) ! dsymmv form
  REAL(iwp),ALLOCATABLE :: diag_precon_pp(:),p_pp(:),r_pp(:),x_pp(:),xnew_pp(:)
  REAL(iwp),ALLOCATABLE :: u_pp(:),pmul_pp(:,:),utemp_pp(:,:),d_pp(:),timest(:)
  REAL(iwp),ALLOCATABLE :: diag_precon_tmp(:,:),eld_pp(:,:),tensor_pp(:,:,:)
  REAL(iwp),ALLOCATABLE :: valf(:),store_pp(:)
  REAL(iwp),ALLOCATABLE :: fun(:),shape_integral_pp(:,:)
  REAL(iwp),ALLOCATABLE :: stress_integral_pp(:,:),stressnodes_pp(:)
  REAL(iwp),ALLOCATABLE :: principal_integral_pp(:,:),princinodes_pp(:)
  REAL(iwp),ALLOCATABLE :: principal(:),reacnodes_pp(:)
  REAL(iwp),ALLOCATABLE :: tempres(:),prop(:,:),ewea(:,:)
  INTEGER,  ALLOCATABLE :: rest(:,:),g_num_pp(:,:),g_g_pp(:,:),node(:)
  INTEGER,  ALLOCATABLE :: no(:),no_pp(:),no_pp_temp(:),sense(:)
  INTEGER,  ALLOCATABLE :: etype_pp(:)

!------------------------------------------------------------------------------
! 2a. Declare dynamic arrays required by UMAT
!
!     stress(ntens)       == sigma(nst)
!     ddsdde(ntens,ntens) == dee(nst,nst)
!     dstran(ntens)       == eps(nst)
!     coords()            == points()
!------------------------------------------------------------------------------
  
  REAL(iwp),ALLOCATABLE :: statev(:),ddsddt(:),drplde(:),drpldt(:)
  REAL(iwp),ALLOCATABLE :: stran(:),time(:),predef(:),dpred(:),props(:)
  REAL(iwp),ALLOCATABLE :: drot(:,:),dfgrd0(:,:),dfgrd1(:,:)

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

  CALL read_xx1(job_name,numpe,e,element,fixed_freedoms,limit,loaded_nodes, &
                 meshgen,nels,nip,nn,nod,nr,partitioner,tol,v)
 
  nprops   = 2   ! needs to go in .dat file and be read using READ_XX1
  np_types = 2   ! needs to go in .dat file and be read using READ_XX1
 
  CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)

  ndof = nod*nodof
  ntot = ndof
 
  ALLOCATE(g_num_pp(nod, nels_pp)) 
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp)) 
  ALLOCATE(rest(nr,nodof+1)) 
  ALLOCATE(etype_pp(nels_pp))
  ALLOCATE(prop(nprops,np_types))
 
  g_num_pp  = 0
  g_coord_pp= zero
  rest      = 0
  etype_pp  = 0
  prop      = 0

  timest(2) = elap_time()

  IF(io_binary) THEN
    CALL read_g_num_pp_be(job_name,iel_start,nn,npes,numpe,g_num_pp)
  ELSE 
    CALL read_g_num_pp(job_name,iel_start,nn,npes,numpe,g_num_pp)
  END IF

  timest(3) = elap_time()

  IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
  timest(4) = elap_time()

  IF(io_binary) THEN
    CALL read_g_coord_pp_be(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
  ELSE
    CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
  END IF

  timest(5) = elap_time()

  CALL read_rest(job_name,numpe,rest)

  timest(6) = elap_time()

  fname = job_name(1:LEN_TRIM(job_name)) // ".mat"
  CALL read_materialValue(prop,fname,numpe,npes)

  IF(io_binary) THEN
    CALL read_etype_pp_be(job_name,npes,numpe,etype_pp)
  ELSE
    PRINT *, "Different material types not supported for ASCII files"
    CALL shutdown()
  END IF

  timest(7) = elap_time()
   
  PRINT *, "Read job data on process ", numpe
 
!------------------------------------------------------------------------------
! 4. Allocate dynamic arrays used in main program
!------------------------------------------------------------------------------

  ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim),principal(ndim),      &
           der(ndim,nod),deriv(ndim,nod),bee(nst,ntot),                       &
           eld(ntot,1),eps(nst),sigma(nst),                                   &
           pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp),                      &
           weights(nip),g_g_pp(ntot,nels_pp),fun(nod))

  IF(sym_storkm) THEN
    j=ntot
    DO i=1,ntot-1
      j=j+ntot-i
    END DO
    ALLOCATE(vstorkm_pp(j,nels_pp),km(ntot,ntot))
    vstorkm_pp=zero ; km=zero
  ELSE
    ALLOCATE(storkm_pp(ntot,ntot,nels_pp)); storkm_pp=zero
  END IF

  PRINT *, "Allocated dynamic arrays on process ", numpe

!------------------------------------------------------------------------------
! 4a. Allocate dynamic arrays used in the UMAT
!------------------------------------------------------------------------------

  ALLOCATE(props(nprops))

  props     = zero
  props(1)  = e
  props(2)  = v

!------------------------------------------------------------------------------
! 5. Loop the elements to find the steering array and the number of equations
!    to solve.
!------------------------------------------------------------------------------

  CALL rearrange(rest)
  
  g_g_pp = 0 ; neq = 0

  elements_1: DO iel = 1, nels_pp
    CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
!   CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest)
  END DO elements_1

  neq = MAXVAL(g_g_pp); neq= MAX_P(neq)
 
  timest(8) = elap_time()

  PRINT *, "Looped the elements to find the steering array on process ", numpe

!------------------------------------------------------------------------------
! 6. Create interprocessor communication tables
!------------------------------------------------------------------------------

  CALL calc_neq_pp
  CALL calc_npes_pp(npes,npes_pp)
  CALL make_ggl(npes_pp,npes,g_g_pp)
 
  timest(9) = elap_time()

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  PRINT *, "Created interprocessor communication tables on process ", numpe

!------------------------------------------------------------------------------
! 7. Allocate arrays dimensioned by neq_pp 
!------------------------------------------------------------------------------
  
  ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),            &
           u_pp(neq_pp),d_pp(neq_pp),diag_precon_pp(neq_pp))

  p_pp    = zero  ;  r_pp = zero  ;  x_pp = zero
  xnew_pp = zero  ;  u_pp = zero  ;  d_pp = zero  ; diag_precon_pp = zero

  timest(10) = elap_time()

  PRINT *, "Allocated arrays dimensioned by neq_pp on process ", numpe

!------------------------------------------------------------------------------
! 8. Element stiffness integration and storage
!
!    Note that the UMAT subroutine is currently outside the loop
!    elements_3. This is only correct when all the elements are of the same
!    type and have the same material properties.
!------------------------------------------------------------------------------

  CALL sample(element,points,weights)

! dee = zero              ! comment out if different material types
! CALL deemat(dee,e,v)

! dee = zero
! CALL umat(sigma,statev,dee,sse,spd,scd,rpl,ddsddt,drplde,drpldt,stran,      &
!           eps,time,dtime,temp,dtemp,predef,dpred,cmname,ndi,nshr,           &
!           nst,nstatv,props,nprops,points,drot,pnewdt,celent,dfgrd0,dfgrd1,  &
!           iel,npt,layer,kspt,kstep,kinc)

  IF(sym_storkm) THEN
    vstorkm_pp       = zero
    elements_3: DO iel=1,nels_pp
      km = zero
      e  = prop(1,etype_pp(iel))  ! if each element has own properties
      v  = prop(2,etype_pp(iel))
      dee = zero
      CALL deemat(dee,e,v)
      gauss_pts_1: DO i=1,nip
        CALL shape_der(der,points,i)
        jac   = MATMUL(der,g_coord_pp(:,:,iel))
        det   = determinant(jac)
        CALL invert(jac)
        deriv = MATMUL(jac,der)
        CALL beemat(bee,deriv)
        km    = km + MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
      END DO gauss_pts_1
      j = 1
      DO i = 1,ntot
        k = j+ntot-i
        vstorkm_pp(j:k,iel)=km(i:,i)  
        j = k+1
      END DO
    END DO elements_3
  ELSE
    storkm_pp        = zero
    elements_3a: DO iel=1,nels_pp
      PRINT *, "etype_pp for iel ", iel+iel_start-1, " = ",etype_pp(iel)
      e  = prop(1,etype_pp(iel))  ! if each element has own properties
      v  = prop(2,etype_pp(iel))
      dee = zero
      CALL deemat(dee,e,v)
      gauss_pts_1a: DO i=1,nip
        CALL shape_der(der,points,i)
        jac   = MATMUL(der,g_coord_pp(:,:,iel))
        det   = determinant(jac)
        CALL invert(jac)
        deriv = MATMUL(jac,der)
        CALL beemat(bee,deriv)
        storkm_pp(:,:,iel)   = storkm_pp(:,:,iel) +                           &
                             MATMUL(MATMUL(TRANSPOSE(bee),dee),bee) *         &
                             det*weights(i)   
      END DO gauss_pts_1a
    END DO elements_3a
  END IF
 
  timest(11) = elap_time()

!------------------------------------------------------------------------------
! 9. Build the diagonal preconditioner
!------------------------------------------------------------------------------
  
  ALLOCATE(diag_precon_tmp(ntot,nels_pp))
  diag_precon_tmp = zero

  IF(sym_storkm) THEN
    elements_4: DO iel = 1,nels_pp 
      j = 1
      DO i = 1,ndof
        diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel) + vstorkm_pp(j,iel)
        j = j + ndof - i + 1
      END DO
    END DO elements_4
  ELSE 
    elements_4a: DO iel = 1,nels_pp 
      DO i = 1,ndof
        diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel) + storkm_pp(i,i,iel)
      END DO
    END DO elements_4a
  END IF

  CALL scatter(diag_precon_pp,diag_precon_tmp)

  DEALLOCATE(diag_precon_tmp)

  timest(12) = elap_time()

!------------------------------------------------------------------------------
! 10. Read in fixed nodal displacements and assign to equations
!------------------------------------------------------------------------------

  fixed_freedoms_pp = 0

  IF(fixed_freedoms > 0) THEN

    ALLOCATE(node(fixed_freedoms),no(fixed_freedoms),valf(fixed_freedoms),    &
             no_pp_temp(fixed_freedoms),sense(fixed_freedoms))

    node = 0 ; no = 0 ; no_pp_temp = 0 ; sense = 0 ; valf = zero

    CALL read_fixed(job_name,numpe,node,sense,valf)
    CALL find_no(node,rest,sense,no)
!   CALL reindex_fixed_nodes(ieq_start,no,no_pp_temp,fixed_freedoms_pp,       &
    CALL reindex(ieq_start,no,no_pp_temp,fixed_freedoms_pp,                   &
                             fixed_freedoms_start,neq_pp)

    ALLOCATE(no_pp(fixed_freedoms_pp),store_pp(fixed_freedoms_pp))

    no_pp    = 0
    store_pp = zero
    no_pp    = no_pp_temp(1:fixed_freedoms_pp)

    DEALLOCATE(node,no,sense,no_pp_temp)

  END IF

  DEALLOCATE(rest)

!------------------------------------------------------------------------------
! 11. Read in loaded nodes and get starting r_pp
!------------------------------------------------------------------------------
 
  IF(loaded_nodes > 0) THEN

    ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes))
    
    val  = zero ; node = 0

    CALL read_loads(job_name,numpe,node,val)
    CALL load(g_g_pp,g_num_pp,node,val,r_pp(1:))

    tload = SUM_P(r_pp(1:))

    DEALLOCATE(node,val)

  ELSE

    tload = zero

  END IF
  
  DEALLOCATE(g_g_pp)
  
  timest(13) = elap_time()

!------------------------------------------------------------------------------
! 12. Invert the preconditioner.
!     If there are fixed freedoms, first apply a penalty
!------------------------------------------------------------------------------

  IF(fixed_freedoms_pp > 0) THEN
    DO i = 1, fixed_freedoms_pp
       j                 = no_pp(i) - ieq_start + 1
       diag_precon_pp(j) = diag_precon_pp(j) + penalty
       store_pp(i)       = diag_precon_pp(j)
    END DO
  END IF

  diag_precon_pp = 1._iwp/diag_precon_pp

!------------------------------------------------------------------------------
! 13. Initiallize preconditioned conjugate gradient
!------------------------------------------------------------------------------

  IF(fixed_freedoms_pp > 0) THEN
    DO i = 1,fixed_freedoms_pp
       j       = no_pp(i) - ieq_start + 1
       k       = fixed_freedoms_start + i - 1
       r_pp(j) = store_pp(i) * valf(k)
    END DO
  END IF

  d_pp  = diag_precon_pp*r_pp
  p_pp  = d_pp
  x_pp  = zero

!------------------------------------------------------------------------------
! 14. Preconditioned conjugate gradient iterations
!------------------------------------------------------------------------------

  iters = 0

  iterations: DO 
    iters    = iters + 1
    u_pp     = zero
    pmul_pp  = zero
    utemp_pp = zero

    CALL gather(p_pp,pmul_pp)

    IF(sym_storkm) THEN
      elements_5: DO iel=1,nels_pp
        CALL dspmv('L',ntot,1.0,vstorkm_pp(:,iel),pmul_pp(:,iel),1,0.0, &
                   utemp_pp(:,iel),1)    
      END DO elements_5
    ELSE
      elements_5a: DO iel=1,nels_pp
        utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
      END DO elements_5a
    END IF

    CALL scatter(u_pp,utemp_pp)

    IF(fixed_freedoms_pp > 0) THEN
      DO i = 1, fixed_freedoms_pp
        j       = no_pp(i) - ieq_start + 1
        u_pp(j) = p_pp(j) * store_pp(i)
      END DO
    END IF

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

  DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,pmul_pp) 

!
! reduce memory usage, but needed for strain energy calculation
! 
! IF(sym_storkm) THEN
!   DEALLOCATE(vstorkm_pp)
! ELSE
!   DEALLOCATE(storkm_pp)
! END IF
!
 
  timest(14) = elap_time()

!------------------------------------------------------------------------------
! 15. Print out displacements, stress, principal stress and reactions
!------------------------------------------------------------------------------

  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
  
  timest(15) = elap_time()

  IF(numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ")-1)//".dis"
    OPEN(24, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".str"
    OPEN(25, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".pri"
    OPEN(26, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".vms"
    OPEN(27, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".rea"
    OPEN(28, file=fname, status='replace', action='write')
  END IF

!------------------------------------------------------------------------------
! 16a. Displacements
!------------------------------------------------------------------------------

  ALLOCATE(eld_pp(ntot,nels_pp)) 
  eld_pp = zero
  CALL gather(xnew_pp(1:),eld_pp)
  DEALLOCATE(xnew_pp)

  ALLOCATE(disp_pp(nodes_pp*ndim))
  disp_pp = zero

  label   = "*DISPLACEMENT"

  CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,              &
                     node_start,node_end,eld_pp,disp_pp,1)

  timest(16) = elap_time()

  IF(io_binary) THEN

    IF(numpe==1) THEN
     
      fname = job_name(1:INDEX(job_name, " ")-1)//".bin.ensi.DISPL-000001"
      OPEN(40,file=fname,status='replace',action='write',                     & 
           form='unformatted',access='stream')

      cbuffer="Alya Ensight Gold --- Vector per-node variable file"
      WRITE(40) cbuffer
      cbuffer="part"        ; WRITE(40) cbuffer
      WRITE(40) int(1,kind=c_int)
      cbuffer="coordinates" ; WRITE(40) cbuffer 
   
    END IF
      
    ALLOCATE(tempres(nodes_pp))
    tempres = zero

    DO i=1,ndim; tempres=zero        
      DO j=1,nodes_pp
        k=i+(ndim*(j-1))
        tempres(j)=disp_pp(k)
      END DO
      CALL dismsh_ensi_pb2(40,1,nodes_pp,npes,numpe,1,tempres)
    END DO

    DEALLOCATE(tempres)

    CLOSE(40)

  ELSE 

    CALL write_nodal_variable(label,24,1,nodes_pp,npes,numpe,ndim,disp_pp)

  END IF

  DEALLOCATE(disp_pp)

  timest(17) = elap_time()

  IF(numpe==1) CLOSE(24)
  
!------------------------------------------------------------------------------
! 16b. Stresses
!------------------------------------------------------------------------------

!  ALLOCATE(shape_integral_pp(nod,nels_pp))
!  ALLOCATE(stress_integral_pp(nod*nst,nels_pp))
!  ALLOCATE(stressnodes_pp(nodes_pp*nst))
!  ALLOCATE(principal_integral_pp(nod*nodof,nels_pp))
!  ALLOCATE(princinodes_pp(nodes_pp*nodof))
!  ALLOCATE(reacnodes_pp(nodes_pp*nodof))
  
!  shape_integral_pp     = zero
!  stress_integral_pp    = zero
!  stressnodes_pp        = zero
!  principal_integral_pp = zero  
!  princinodes_pp        = zero
!  reacnodes_pp          = zero
!  utemp_pp              = zero
!  dee                   = zero
!  sigma                 = zero

!  DO iel = 1,nels_pp
!    DO i = 1,nip
!      CALL shape_der(der,points,i)
!      jac   = MATMUL(der,g_coord_pp(:,:,iel))
!      det   = DETERMINANT(jac) 
!      CALL invert(jac)
!      deriv = MATMUL(jac,der)
!      CALL beemat(bee,deriv)
!      eps   = MATMUL(bee,eld_pp(:,iel))
!      sigma = MATMUL(dee,eps) ! umat replaces this line
!      dee   = zero
!      sigma = zero
!      CALL umat(sigma,statev,dee,sse,spd,scd,rpl,ddsddt,drplde,drpldt,stran,  &
!               eps,time,dtime,temp,dtemp,predef,dpred,cmname,ndi,nshr,       &
!               nst,nstatv,props,nprops,points,drot,pnewdt,celent,dfgrd0,     &
!               dfgrd1,iel,npt,layer,kspt,kstep,kinc)
!      CALL PRINCIPALSTRESS3D(sigma,principal)
!      utemp_pp(:,iel) = utemp_pp(:,iel) +                                     &
!                        MATMUL(TRANSPOSE(bee),sigma)*det*weights(i)

!      CALL shape_fun(fun,points,i)

!      DO j = 1,nod
!        idx1 = (j-1)*nst
!        idx2 = (j-1)*nodof
!        shape_integral_pp(j,iel) = shape_integral_pp(j,iel) +                 &
!                                   fun(j)*det*weights(i)
!        DO k = 1,nst
!          stress_integral_pp(idx1+k,iel) = stress_integral_pp(idx1+k,iel) +   &
!                                           fun(j)*sigma(k)*det*weights(i)
!        END DO
!        DO k = 1,nodof
!          principal_integral_pp(idx2+k,iel) = principal_integral_pp(idx2+k,iel) + &
!                                              fun(j)*principal(k)*det*weights(i)
!        END DO
!      END DO
!    END DO !gauss
!  END DO !elements
    
!------------------------------------------------------------------------------
! 16c. Stress
!------------------------------------------------------------------------------

!  label = "*STRESS"

!  CALL nodal_projection(npes,nn,nels_pp,g_num_pp,nod,nst,nodes_pp,            &
!                        node_start,node_end,shape_integral_pp,                &
!                        stress_integral_pp,stressnodes_pp)
!  CALL write_nodal_variable(label,25,1,nodes_pp,npes,numpe,nst,               &
!                            stressnodes_pp)
                            
!  DEALLOCATE(stress_integral_pp,stressnodes_pp)

!  IF(numpe==1) CLOSE(25)

!  PRINT *, "Write Stress"

!------------------------------------------------------------------------------
! 16d. Principal stress
!------------------------------------------------------------------------------
  
!  label = "*PRINCIPAL STRESS"
  
!  CALL nodal_projection(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,          &
!                        node_start,node_end,shape_integral_pp,                &
!                        principal_integral_pp,princinodes_pp)
!  CALL write_nodal_variable(label,26,1,nodes_pp,npes,numpe,nodof,             &
!                            princinodes_pp)
                            
!  DEALLOCATE(principal_integral_pp)

!  IF(numpe==1) CLOSE(26)

!------------------------------------------------------------------------------
! 16e. Von Mises stress (rho_v)
!      rho_v = sqrt( ( (rho1-rho2)^2 + (rho2-rho3)^2 + (rho1-rho3)^2 ) / 2 )
!------------------------------------------------------------------------------
  
!  label = "*MISES STRESS"
  
!  DO i = 1,nodes_pp
!    j = ((i-1)*nodof)+1
!    k = j + 1
!    l = j + 2
!    princinodes_pp(j) = SQRT(((princinodes_pp(j)-princinodes_pp(k)) **2 +     &
!                              (princinodes_pp(k)-princinodes_pp(l)) **2 +     &
!                              (princinodes_pp(j)-princinodes_pp(l)) **2)      &
!                              * 0.5_iwp)
!    princinodes_pp(k:l) = zero
!  END DO

!  CALL write_nodal_variable(label,27,1,nodes_pp,npes,numpe,nodof,             &
!                            princinodes_pp)
                            
!  DEALLOCATE(princinodes_pp)

!  IF(numpe==1) CLOSE(27)

!------------------------------------------------------------------------------
! 16f. Reactions
!------------------------------------------------------------------------------

!  label = "*NODAL REACTIONS"
!  CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,             &
!                     node_start,node_end,utemp_pp,reacnodes_pp,0)
!  CALL write_nodal_variable(label,28,1,nodes_pp,npes,numpe,nodof,             &
!                            reacnodes_pp)
!  DEALLOCATE(reacnodes_pp,shape_integral_pp)

!  IF(numpe==1) CLOSE(28)

  timest(18) = elap_time()

!------------------------------------------------------------------------------
! 16g. Strain energy
!
!      http://solidmechanics.org/text/Chapter7_2/Chapter7_2.htm
!
!      Element strain energy W = 1/2 * uT * Ke * u
!
!      Where Ke = element stiffness, uT = Transpose of element displacements 
!      and u = element displacements
!
!      If this section turns out to be slow, it can be optimised
!------------------------------------------------------------------------------

  ewt_pp     = zero
  ewt        = zero
  utemp_pp   = zero
  eld        = zero

  ALLOCATE(ewea(1,1))

  DO iel = 1,nels_pp
    IF(sym_storkm) THEN ! rebuild km
      j  = 1
      km = zero
      DO i = 1,ntot
         k = j+ntot-i
         km(i:,i) = vstorkm_pp(j:k,iel)
         km(i,i:) = vstorkm_pp(j:k,iel)
         j = k+1
      END DO
    ELSE
      km = storkm_pp(:,:,iel)
    END IF 
    ewea     = zero
    eld(:,1) = eld_pp(:,iel)
    ewea     = MATMUL(MATMUL(TRANSPOSE(eld),km),eld)*0.5_iwp 
    ewt_pp = ewt_pp + ewea(1,1)
  END DO

  CALL MPI_ALLREDUCE(ewt_pp,ewt,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)

  DEALLOCATE(ewea)

  timest(19) = elap_time()

!------------------------------------------------------------------------------
! 17. Output performance data
!------------------------------------------------------------------------------

  IF(numpe==1) THEN
    CALL WRITE_XX1(ewt,fixed_freedoms,iters,job_name,loaded_nodes,neq,nn,     &
                   npes,nr,numpe,timest,tload)
  END IF
 
  CALL MPI_BARRIER(mpi_comm_world,ier)
  
  CALL shutdown() 
 
END PROGRAM XX1
