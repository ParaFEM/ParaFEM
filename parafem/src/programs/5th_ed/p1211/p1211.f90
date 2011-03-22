 PROGRAM p1211
!------------------------------------------------------------------------------
!      program 12.11 three-d strain of an elastic-plastic(von Mises) solid
!      using 20-node brick elements; consistent return method, pcg parallel
!      displacement control, load control, adaptive increment sizes
!------------------------------------------------------------------------------

 USE precision ; USE global_variables ; USE mp_interface
 USE input     ; USE output           ; USE loading
 USE timing    ; USE maths            ; USE gather_scatter
 USE partition ; USE elements         ; USE steering
 USE pcg       ; USE plasticity

 IMPLICIT NONE

!-------------------------------------------------------------------------------
! 1. Declare program variables.
!-------------------------------------------------------------------------------
! fixed_freedoms_pp    = num_no in book
! fixed_freedoms_start = no_index_start in book
! neq,ntot are use associated

 INTEGER, PARAMETER   :: nodof=3,nst=6,ndim=3
 INTEGER, PARAMETER   :: nprops=3              ! no of material properties 
 INTEGER, PARAMETER   :: nod=20                ! no of nodes per element 

 INTEGER              :: argc,iargc,meshgen
 INTEGER              :: cjiters
 INTEGER              :: cjits
 INTEGER              :: cjtot
 INTEGER              :: iel,i,j,k,ii,jj,iy    ! counters
 INTEGER              :: fileIndexLength
 INTEGER              :: fixed_freedoms        ! number of fixed displacements
 INTEGER              :: fixed_freedoms_pp     ! number on local processor
 INTEGER              :: fixed_freedoms_start  ! starting freedom number
 INTEGER              :: loaded_nodes          ! number of loaded nodes
 INTEGER              :: loadIncrement         ! increment counter
 INTEGER              :: loadIncrementMax      ! maximum permitted load incs
 INTEGER              :: nels,nn,nr,nip
 INTEGER              :: node_end
 INTEGER              :: nodes_pp
 INTEGER              :: node_start
 INTEGER              :: np_types              ! no of property types
 INTEGER              :: npes_pp
 INTEGER              :: numSteps              ! number of load steps
 INTEGER              :: outputIncrement       ! increment counter for results
 INTEGER              :: plasiters
 INTEGER              :: plasitersMin          ! for adaptive load stepping
 INTEGER              :: plasitersMax          ! for adaptive load stepping

 REAL(iwp),PARAMETER  :: zero    = 0.0_iwp
 REAL(iwp),PARAMETER  :: one     = 1.0_iwp
 REAL(iwp),PARAMETER  :: penalty = 1.e20_iwp

 REAL(iwp)            :: e
 REAL(iwp)            :: v
 REAL(iwp)            :: sbary
 REAL(iwp)            :: dlam,dslam,dsbar
 REAL(iwp)            :: lode_theta
 REAL(iwp)            :: det,fnew,ff,fstiff
 REAL(iwp)            :: top,bot,sigm,tload
 REAL(iwp)            :: tloads
 REAL(iwp)            :: plastol,cjtol
 REAL(iwp)            :: ltol
 REAL(iwp)            :: fftol
 REAL(iwp)            :: up,alpha,beta
 REAL(iwp)            :: tinc,tplas,maxDisp
 REAL(iwp)            :: factor                ! load increment size
 REAL(iwp)            :: fractionApplied       ! fraction of total load
 REAL(iwp)            :: tolerance             ! convergence test
 REAL(iwp)            :: bdyld_l2n             ! l2 norm of bdylds_pp
 REAL(iwp)            :: ddyld_l2n             ! l2 norm of ddylds_pp

 CHARACTER(LEN=50)    :: program_name='p1211'
 CHARACTER(LEN=15)    :: element
 CHARACTER(LEN=50)    :: fname
 CHARACTER(LEN=50)    :: fname_dat
 CHARACTER(LEN=50)    :: fname_loads
 CHARACTER(LEN=50)    :: fname_results
 CHARACTER(LEN=50)    :: fname_dis
 CHARACTER(LEN=50)    :: job_name
 CHARACTER(LEN=50)    :: text
 CHARACTER(LEN=6)     :: step
 CHARACTER(LEN=6)     :: increment
 
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

 REAL(iwp),ALLOCATABLE :: loads_pp(:),points(:,:),bdylds_pp(:)
 REAL(iwp),ALLOCATABLE :: pmul_pp(:,:),dee(:,:),jac(:,:),weights(:)
 REAL(iwp),ALLOCATABLE :: oldis_pp(:),der(:,:),deriv(:,:),bee(:,:),km(:,:)
 REAL(iwp),ALLOCATABLE :: eld(:),eps(:),sigma(:),bload(:),eload(:),elso(:)
 REAL(iwp),ALLOCATABLE :: ddylds_pp(:),dl(:,:),dl_old(:,:),dload(:)
 REAL(iwp),ALLOCATABLE :: stressv(:),qinc(:),store_pp(:),dtemp_pp(:,:)
 REAL(iwp),ALLOCATABLE :: p_pp(:),x_pp(:),xnew_pp(:),u_pp(:),utemp_pp(:,:)
 REAL(iwp),ALLOCATABLE :: d_pp(:),diag_precon_tmp(:,:)
 REAL(iwp),ALLOCATABLE :: vmfl(:),caflow(:),dsigma(:),ress(:),rmat(:,:)
 REAL(iwp),ALLOCATABLE :: acat(:,:),acatc(:,:),qmat(:,:),qinva(:),daatd(:,:)
 REAL(iwp),ALLOCATABLE :: vmflq(:),vmfla(:),qinvr(:),vmtemp(:,:),temp(:,:)
 REAL(iwp),ALLOCATABLE :: disp_pp(:,:),xnewel_pp(:,:)
 REAL(iwp),ALLOCATABLE :: xnewnodes_pp(:),fext_pp(:),timest(:)
 REAL(iwp),ALLOCATABLE :: storkm_pp(:,:,:), storkm_pp_old(:,:,:)
 REAL(iwp),ALLOCATABLE :: tensor_pp(:,:,:), tensor_pp_old(:,:,:)
 REAL(iwp),ALLOCATABLE :: totd_pp(:), totd_pp_old(:)
 REAL(iwp),ALLOCATABLE :: diag_precon_pp(:), diag_precon_pp_old(:)
 REAL(iwp),ALLOCATABLE :: prop(:,:)          ! element properties matrix
 REAL(iwp),ALLOCATABLE :: loadNodeValue(:,:) ! nodal values of applied forces
 REAL(iwp),ALLOCATABLE :: loadEqnValue_pp(:) ! eqn values of applied forces
 REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:)  ! nodal coordinates
 REAL(iwp),ALLOCATABLE :: val(:)             ! prescribed load/disp values

 INTEGER, ALLOCATABLE  :: g_g_pp(:,:)
 INTEGER, ALLOCATABLE  :: g_num_pp(:,:)
 INTEGER, ALLOCATABLE  :: no(:)             ! freedoms to be loaded/fixed 
 INTEGER, ALLOCATABLE  :: node(:)           ! loaded nodes vector
 INTEGER, ALLOCATABLE  :: no_pp(:)          ! freedoms on local processor
 INTEGER, ALLOCATABLE  :: no_pp_temp(:)     ! temporary array
 INTEGER, ALLOCATABLE  :: rest(:,:)
 INTEGER, ALLOCATABLE  :: sense(:)          ! sense of freedoms to be fixed
 INTEGER, ALLOCATABLE  :: etype_pp(:)       ! element property type vector
 INTEGER, ALLOCATABLE  :: loadNodeNum(:)    ! node numbers where forces applied
 INTEGER, ALLOCATABLE  :: loadEqnNum_pp(:)  ! eqn numbers where forces applied
 INTEGER, ALLOCATABLE  :: loadStepInfo(:,:) ! (1,j) = fixed_nodes at step j
                                            ! (2,j) = loaded_nodes at step j

!------------------------------------------------------------------------------
! 3. Declare logical variables. If displacement_control == .false. then the 
!    program will use load_control 
!------------------------------------------------------------------------------

 LOGICAL               :: plastic_converged     = .false.
 LOGICAL               :: cj_converged          = .false.
 LOGICAL               :: printLast             = .false.
 LOGICAL               :: displacement_control  = .false.

!------------------------------------------------------------------------------
! 4. Start timer and initialize MPI
!------------------------------------------------------------------------------

  ALLOCATE(timest(20))
  timest(1) = elap_time()
  CALL find_pe_procs(numpe,npes)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier) !- Barrier inserted for debugging

!------------------------------------------------------------------------------
! 5. Read job_name from command line and set filenames for I/O
!------------------------------------------------------------------------------

  argc = iargc()

! argc=1 !hack for terra

  IF (argc /= 1) CALL job_name_error(numpe,program_name)
  
  call GETARG(1, job_name)

  fname_dat      = job_name(1:INDEX(job_name, " ")-1) // ".dat"
  fname_results  = job_name(1:INDEX(job_name, " ")-1) // ".res"
  fname_dis      = job_name(1:INDEX(job_name, " ")-1) // ".dis"

!------------------------------------------------------------------------------
! 6. Master processor reads data from control file job_name.dat and distributes
!    data to all slave processors
!------------------------------------------------------------------------------

  IF(numpe==1) THEN
    OPEN (10, file=fname_dat, status='old', action='read')
    READ (10,*) element,meshgen,                                               &
                nels,nn,nr,nip,                                                &
                plasitersMax,plasitersMin,loadIncrementMax,cjits,plastol,      &
                cjtol,fftol,ltol,np_types,numSteps
    OPEN (11, file=fname_results, status='replace', action='write')
  END IF

  CALL bcast_inputdata_p1211(numpe,npes,element,meshgen,nels,nn,nr,nip,        &
                             plasitersMax,plasitersMin,loadIncrementMax,       &
                             cjits,plastol,cjtol,fftol,ltol,np_types,numSteps)
    
  CALL check_inputdata_p1211(numpe,npes,element,meshgen,nels,nn,nr,nip,        &
                             plasitersMax,plasitersMin,loadIncrementMax,       &
                             cjits,plastol,cjtol,fftol,ltol,np_types,numSteps)

  ntot     = nod * nodof

  ALLOCATE(loadStepInfo(2,numSteps))
  loadStepInfo = 0

  IF(numpe==1) THEN
    DO i = 1, numsteps
      READ(10,*) loadStepInfo(1,i)
      READ(10,*) loadStepInfo(2,i)
    END DO
  END IF

  CALL MPI_BCAST(loadStepInfo,2*numSteps,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  IF(numpe==1) CLOSE(10)

!------------------------------------------------------------------------------
! 7. Create array to store the sizes of successful increments
!------------------------------------------------------------------------------

  ALLOCATE(qinc(loadIncrementMax))
  qinc = zero

!-------------------------------------------------------------------------------
! 8. The master processor imports the mesh, material properties and restraints. 
!    Each processor stores the information that is needed locally.
!-------------------------------------------------------------------------------

  ! Elements
  
  CALL calc_nels_pp(nels)
! fname = fname_base(1:INDEX(fname_base, " ")-1) // ".psize"
! CALL readall_nels_pp_fname(fname)

  ALLOCATE(g_num_pp(nod, nels_pp)) 
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp)) 
  ALLOCATE(prop(nprops,np_types))
  ALLOCATE(etype_pp(nels_pp))
  ALLOCATE(rest(nr,nodof+1))

  g_num_pp   = 0
  g_coord_pp = zero
  prop       = zero
  etype_pp   = 0
  rest       = 0

! CALL read_g_num_pp(job_name,iel_start,nels,nn,numpe,g_num_pp)
  CALL read_g_num_pp2(job_name,iel_start,nn,npes,numpe,g_num_pp)

  IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)

  CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
  CALL read_rest(job_name,numpe,rest)

  fname     = job_name(1:INDEX(job_name, " ")-1) // ".d"
  CALL read_materialID_pp(etype_pp,fname,nn,ndim,nels,nod,iel_start,          &
                          numpe,npes)

  fname     = job_name(1:INDEX(job_name, " ")-1) // ".mat"
  CALL read_materialValue(prop,fname,numpe,npes)

  CALL read_rest(job_name,numpe,rest)

  ! I/O error capture
  ! Will only work on master processor. Needs implementing properly

  DO iel = 1,nels_pp
    IF(etype_pp(iel) > np_types) THEN
      IF(numpe == 1) THEN
        WRITE(11,'(/2A)') "------------------------------------------",        &
                          "----------"
        WRITE(11,'(A)')   "Program has detected a fatal error"
        WRITE(11,'(/A,I6,A/)') "  Material ",etype_pp(iel),                    &
                               " has no properties" 
        WRITE(11,'(A)')   "Analysis aborted"
        WRITE(11,'(2A/)') "------------------------------------------",        &
                          "----------"
!       CALL FLUSH_(11)
      END IF
      CALL shutdown()
    END IF
  END DO

!-------------------------------------------------------------------------------
! 9. Allocate arrays
!-------------------------------------------------------------------------------

  ALLOCATE (points(nip,ndim),weights(nip),g_g_pp(ntot,nels_pp),dee(nst,nst),   &
            stressv(nst),jac(ndim,ndim),dtemp_pp(ntot,nels_pp),elso(nst),      &
            der(ndim,nod),deriv(ndim,nod),utemp_pp(ntot,nels_pp),              &
            bee(nst,ntot),km(ntot,ntot),eld(ntot),eps(nst),sigma(nst),         &
            bload(ntot),eload(ntot),vmfl(nst),qinvr(nst),dload(ntot),ress(nst),&
            pmul_pp(ntot,nels_pp),caflow(nst),dsigma(nst),dl(nip,nels_pp),     &
            rmat(nst,nst),acat(nst,nst),acatc(nst,nst),qmat(nst,nst),          &
            qinva(nst),daatd(nst,nst),vmflq(nst),vmfla(nst),vmtemp(1,nst),     &
            storkm_pp(ntot,ntot,nels_pp), storkm_pp_old(ntot,ntot,nels_pp),    &
            tensor_pp(nst,nip,nels_pp), tensor_pp_old(nst,nip,nels_pp),        &
            dl_old(nip,nels_pp))

!-------------------------------------------------------------------------------
! 9. Create node freedom array, find element steering and neq
!-------------------------------------------------------------------------------

  CALL rearrange(rest) ! does REST need to be in ascending node order?

  g_g_pp = 0

  elements_0a: DO iel = 1, nels_pp
    CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest)
  END DO elements_0a

  neq = 0
  elements_0b: DO iel = 1, nels_pp  
    i = MAXVAL(g_g_pp(:,iel))
    IF(i > neq) neq = i
  END DO elements_0b  

  neq = max_integer_p(neq)
  
  CALL calc_neq_pp

  IF(numpe==1) THEN
    WRITE(11,'(A,I12)')    "Number of elements on processor one:    ", nels_pp
    WRITE(11,'(A,I12)')    "Number of equations in problem:         ", neq
    WRITE(11,'(A,I12,/)')  "Number of equations per processor:      ", neq_pp
!   CALL FLUSH_(11)
  END IF 

  CALL calc_npes_pp(npes,npes_pp)
  CALL make_ggl(npes_pp,npes,g_g_pp)

!------------------------------------------------------------------------------ 
! 10. Allocate and initialize arrays dimensioned by NEQ_PP 
!------------------------------------------------------------------------------
 
  ALLOCATE(loads_pp(neq_pp),bdylds_pp(neq_pp),oldis_pp(neq_pp),p_pp(neq_pp),   &
           totd_pp(neq_pp),diag_precon_pp_old(neq_pp),                         &
           xnew_pp(neq_pp),u_pp(neq_pp),diag_precon_pp(neq_pp),d_pp(neq_pp),   &
           ddylds_pp(neq_pp),totd_pp_old(neq_pp),fext_pp(neq_pp),x_pp(neq_pp))

  loads_pp           = zero
  bdylds_pp          = zero           
  oldis_pp           = zero
  totd_pp            = zero 
  totd_pp_old        = zero
  tensor_pp          = zero
  tensor_pp_old      = zero
  p_pp               = zero
  xnew_pp            = zero
  diag_precon_pp     = zero
  diag_precon_pp_old = zero
  fext_pp            = zero
  dl                 = zero
  dl_old             = zero

!------------------------------------------------------------------------------
! 18. Create temp array needed for FMRMAT and FMACAT
!------------------------------------------------------------------------------

  ALLOCATE(temp(nst,nst))
  temp = zero
  CALL form_temp(temp)
  
!------------------------------------------------------------------------------ 
! 11. Record time taken to import data and setup analysis
!------------------------------------------------------------------------------
 
  IF(numpe==1) THEN
    WRITE(11,'(A,F10.4/)')  "Time after setup:                         ",     &
                             elap_time( ) - timest(1)
    WRITE(11,'(2A)')    "------------------------------------------",         &
                        "----------"
    WRITE(11,'(2A)')    "                   ANALYSIS DATA          ",         &
                        "          "
    WRITE(11,'(2A)')    "------------------------------------------",         &
                        "----------"
!   CALL FLUSH_(11)
  END IF 

!-----------------------------------------------------------------------------
! 12. Initial element stiffness integration before entering load stepping loop
!-----------------------------------------------------------------------------

  CALL sample(element,points,weights)
  
  storkm_pp     = zero
  storkm_pp_old = zero

  elements_2: DO iel = 1 , nels_pp
              
              e          = prop(2,etype_pp(iel))
              v          = prop(3,etype_pp(iel))

              CALL deemat(e,v,dee) ! Fran's version needs checking

              gauss_pts_1: DO i =1 , nip    
                CALL shape_der (der,points,i)
                jac   = MATMUL(der,g_coord_pp(:,:,iel))
                det   = determinant(jac)
                CALL invert(jac)
                deriv = MATMUL(jac,der)
                CALL beemat (deriv,bee)           
                storkm_pp(:,:,iel) =  storkm_pp(:,:,iel)  +                   &
                                      MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)  &
                                     *det*weights(i)
             END DO gauss_pts_1
             
  END DO elements_2

!-------------------------------------------------------------------------------
! 13. Compute initial diagonal preconditioner, needs inverting later
!-------------------------------------------------------------------------------
 
  ALLOCATE(diag_precon_tmp(ntot,nels_pp))
  diag_precon_tmp    = zero
  diag_precon_pp     = zero
  diag_precon_pp_old = zero

  elements_1: DO iel = 1,nels_pp 
    DO k = 1,ntot
      diag_precon_tmp(k,iel) = diag_precon_tmp(k,iel) + storkm_pp(k,k,iel)
    END DO
  END DO elements_1

  CALL scatter(diag_precon_pp,diag_precon_tmp)

  DEALLOCATE(diag_precon_tmp)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 14. Start load step loop
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
  load_steps : DO ii = 1 , numSteps 
  
  IF(numpe==1) THEN
    WRITE(11,'(/,A,I12)') "LOAD STEP                               ", ii
!   CALL FLUSH_(11)
  END IF

  qinc            = zero ! collect successful increment sizes for each load step
  fext_pp         = zero ! set external loads vector to zero
  fixed_freedoms  = loadStepInfo(1,ii)
  loaded_nodes    = loadStepInfo(2,ii)

  IF(fixed_freedoms > 0 .AND. loaded_nodes > 0) THEN
    IF(numpe == 1) THEN
      WRITE(11,'(/2A)') "------------------------------------------",      &
                        "----------"
      WRITE(11,'(A)')   "Program has detected a fatal error"
      WRITE(11,'(/A)') "  Applying fixed_freedoms and loaded_nodes at the"
      WRITE(11,'(A/)') "  same time is not currently supported."
      WRITE(11,'(A)')   "Analysis aborted"
      WRITE(11,'(2A/)') "------------------------------------------",       &
                        "----------"
!     CALL FLUSH_(11)
    END IF
    CALL shutdown()
  END IF

!-------------------------------------------------------------------------------
! 15. Read in fixed nodal displacements and assign to equations 
!-------------------------------------------------------------------------------
  
  IF(fixed_freedoms > 0) THEN
    
    displacement_control = .true. 

    IF(numpe==1) THEN
      WRITE(11,'(A,I12)')   "Number of fixed freedoms:               ",        &
                             fixed_freedoms
      WRITE(11,'(A,I12)')   "Program using displacement control      "
!     CALL FLUSH_(11)
    END IF

    ALLOCATE(node(fixed_freedoms))
    ALLOCATE(no(fixed_freedoms))
    ALLOCATE(no_pp_temp(fixed_freedoms))
    ALLOCATE(sense(fixed_freedoms))
    ALLOCATE(val(fixed_freedoms))

    node        = 0
    no          = 0
    no_pp_temp  = 0
    sense       = 0
    val         = zero

    CALL getFileNumber(step,ii)
    fname = job_name(1:INDEX(job_name, " ")-1) // "_" //                       &
            step(1:INDEX(step, " ")-1)  // ".fix"
    CALL read_fixed(fname,numpe,sense,node,val)
    CALL find_no(node,rest,sense,no)
    CALL reindex_fixed_nodes(ieq_start,no,no_pp_temp,fixed_freedoms_pp,        &
                             fixed_freedoms_start,neq_pp)

    ALLOCATE(no_pp(fixed_freedoms_pp),store_pp(fixed_freedoms_pp))
    no_pp    = 0 
    store_pp = 0
    no_pp    = no_pp_temp(1:fixed_freedoms_pp)   
 
    DEALLOCATE(node)
    DEALLOCATE(no)
    DEALLOCATE(sense)
    DEALLOCATE(no_pp_temp)

  END IF

!-------------------------------------------------------------------------------
! 16. Read in applied forces and assign to equations
!-------------------------------------------------------------------------------
 
  IF(loaded_nodes > 0) THEN

    displacement_control = .false.

    IF(numpe==1) THEN
      WRITE(11,'(A,I12)') "Number of loaded nodes:                 ",          &
                           loaded_nodes
      WRITE(11,'(A,I12)') "Program using load control              "
!     CALL FLUSH_(11)
    END IF

    CALL getFileNumber(step,ii)
 
    fname_loads = job_name(1:INDEX(job_name, " ")-1) // "_" //                 &
!                 step(1:INDEX(step, " ")-1)  // ".lds"
                  step(1:INDEX(step, " ")-1)

    ALLOCATE(loadNodeNum(loaded_nodes))
    ALLOCATE(loadNodeValue(ndim,loaded_nodes))
    
    loadNodeValue = zero
    loadNodeNum   = 0

    CALL read_loads(fname_loads,numpe,loadNodeNum,loadNodeValue)

    CALL load(g_g_pp,g_num_pp,loadNodeNum,loadNodeValue,fext_pp(1:))

    tload = SUM_P(fext_pp(1:))
    IF(numpe==1) THEN
       WRITE(11,'(A,E14.6)') "Total load to be applied:             ",         &
                              tload
    END IF
    tload = 0.0_iwp

    DEALLOCATE(loadNodeNum)
    DEALLOCATE(loadNodeValue)

  END IF

!------------------------------------------------------------------------------
! 17. Invert the preconditioner. 
!     If there are fixed nodes, first apply a penalty.
!------------------------------------------------------------------------------
 
  IF(fixed_freedoms_pp>0) THEN
    DO i=1,fixed_freedoms_pp
      j                 = no_pp(i) - ieq_start + 1
      diag_precon_pp(j) = diag_precon_pp(j) + penalty
      store_pp(i)       = diag_precon_pp(j)
    END DO
  END IF

  IF(ii==1) THEN
    diag_precon_pp     = one/diag_precon_pp
    diag_precon_pp_old = diag_precon_pp
  END IF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 18. Start load increment loop
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
  fractionApplied     = 0.0_iwp   ! fraction of the load applied
  factor              = 1.0_iwp   ! fraction of load to be applied
  iy                  = 1
  loadIncrement       = 0
  outputIncrement     = 0

  load_increments : do

  loadIncrement      = loadIncrement + 1

!------------------------------------------------------------------------------
! 19. Save old data in case current load increment is abandoned and zero arrays
!     used in plasticity loop
!------------------------------------------------------------------------------

  tensor_pp_old      = tensor_pp
  storkm_pp_old      = storkm_pp
  diag_precon_pp_old = diag_precon_pp
  totd_pp_old        = totd_pp
  dl_old             = dl

  ! Arrays to zero

  dee                = zero
  elso               = zero
  stressv            = zero
  jac                = zero
  der                = zero
  bee                = zero
  eld                = zero
  eps                = zero
  sigma              = zero
  bload              = zero
  eload              = zero
  vmfl               = zero
  qinvr              = zero
  dload              = zero
  caflow             = zero
  dsigma             = zero
  ress               = zero
  rmat               = zero
  acat               = zero
  acatc              = zero
  qmat               = zero
  qinva              = zero
  daatd              = zero
  vmflq              = zero
  vmfla              = zero
  vmtemp             = zero
   
!------------------------------------------------------------------------------
! 20. Terminate analysis if load increments reaches maximum
!------------------------------------------------------------------------------

  IF(outputIncrement > loadIncrementMax - 1) THEN
    IF(numpe==1) THEN
      WRITE(11,'(/2A)')   "----------------------------------------",         &
                          "------------"
      WRITE(11,'(A,I12)') "Number of load increments greater than: ",         &
                           loadIncrementMax 
      WRITE(11,'(A)')     "Analysis terminated by program"
      WRITE(11,'(A)')     "Try reducing the applied load "
      WRITE(11,'(A,E12.4)') "This analysis took:                     ",       &
                             elap_time( ) - timest(1)
      WRITE(11,'(2A)')    "----------------------------------------",         &
                          "------------"
      WRITE(11,'(/A)') "If running the same analysis again, reduce the "
      WRITE(11,'(A/)') "solution time by using these step sizes:       "
      DO i=1,outputIncrement
        WRITE(11,'(A,I5,A,F12.8)') "Load step ", i, " ", qinc(i)
      END DO
!     CALL FLUSH_(11)
    END IF
    CALL shutdown()
  END IF

  IF(numpe==1) THEN
    WRITE(11,'(/A,I12)')  "LOAD INCREMENT NUMBER                   ",         &
                           loadIncrement
    WRITE(11,'(A,F13.5)') "Current increment size:                ", factor
    IF(displacement_control) THEN 
      WRITE(11,'(/A,I12)')"    Iteration        tload       tloads"
    ELSE
      WRITE(11,'(/A,I12)')"    Iteration    bdyld_l2n    ddyld_l2n    Tolerance"
    END IF
!   CALL FLUSH_(11)
  END IF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 21. Start plastic iteration loop
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
   plasiters = 0
   loads_pp  = zero 
   cjtot     = 0
   tinc      = elap_time()

   plastic_iterations: DO   
    
     plasiters = plasiters + 1

     IF(plasiters==1) THEN

       IF(fixed_freedoms_pp>0) THEN
         DO i = 1,fixed_freedoms_pp
                     j = no_pp(i) - ieq_start + 1
                     k = fixed_freedoms_start + i - 1
           loads_pp(j) = store_pp(i) * val(k) * factor
         END DO
       END IF

       IF(loaded_nodes > 0) loads_pp = fext_pp * factor

     END IF

     IF(plasiters>1) THEN

       loads_pp = zero 
       loads_pp = loads_pp + bdylds_pp

       IF(fixed_freedoms_pp > 0) THEN

         DO i = 1,fixed_freedoms_pp
                     j = no_pp(i) - ieq_start + 1
           loads_pp(j) = zero
         END DO

       END IF

     END IF

     utemp_pp  = zero
     bdylds_pp = zero
     ddylds_pp = zero
     d_pp      = diag_precon_pp * loads_pp
     p_pp      = d_pp
     x_pp      = zero 

!------------------------------------------------------------------------------
! 22. Solve the simultaneous equations by pcg
!------------------------------------------------------------------------------
 
     cjiters = 0
     
     conjugate_gradients:  DO
       
       cjiters = cjiters + 1 
       u_pp    = zero   
       pmul_pp = zero
       
       CALL gather(p_pp,pmul_pp)
       elements_3 : DO iel = 1 , nels_pp
         utemp_pp(:,iel)=MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
!        CALL dgemv('n',ntot,ntot,1.0,km,ntot,                                 &
!                    pmul_pp(:,iel),1,0.0,utemp_pp(:,iel),1)
       END DO elements_3  
       CALL scatter(u_pp,utemp_pp)

       IF(fixed_freedoms_pp>0) THEN
         DO i = 1,fixed_freedoms_pp
           j = no_pp(i) - ieq_start + 1
           IF(plasiters==1) THEN 
             u_pp(j) = p_pp(j) * store_pp(i)
           ELSE 
             u_pp(j) = zero
           END IF
         END DO
       END IF        

       up           = DOT_PRODUCT_P(loads_pp,d_pp)
       alpha        = up/DOT_PRODUCT_P(p_pp,u_pp)
       xnew_pp      = x_pp + p_pp* alpha
       loads_pp     = loads_pp - u_pp*alpha
       d_pp         = diag_precon_pp*loads_pp 
       beta         = DOT_PRODUCT_P(loads_pp,d_pp)/up
       p_pp         = d_pp + p_pp * beta   

       cj_converged = .TRUE.
       
       CALL checon_par(xnew_pp,cjtol,cj_converged,x_pp)
       
       IF(cj_converged) EXIT
       IF(cjiters==cjits) THEN
         IF(numpe==1) THEN
           WRITE(11,'(/2A)')   "----------------------------------------",     &
                               "------------"
           WRITE(11,'(A)')     "Program has detected a fatal error"
           WRITE(11,'(A)')     "  Number of iterations used in PCG has reached "
           WRITE(11,'(A,I8)')  "  the limit defined by the user: ", cjits
           WRITE(11,'(A)')     "  Try increasing the limit CJITS or check the "
           WRITE(11,'(A)')     "  input deck for errors "
           WRITE(11,'(A)')     "Analysis aborted"
           WRITE(11,'(2A)')    "----------------------------------------",     &
                               "------------"
         END IF
         CALL SHUTDOWN()
       END IF
     END DO conjugate_gradients

     cjtot = cjtot + cjiters  

!------------------------------------------------------------------------------
! 23. Compute element stresses
!------------------------------------------------------------------------------

     tplas          = elap_time()
     loads_pp       = xnew_pp     
     pmul_pp        = zero  
     CALL gather(xnew_pp,pmul_pp)
     utemp_pp       = zero
     dtemp_pp       = zero

!------------------------------------------------------------------------------
! 24. Go round the elements
!------------------------------------------------------------------------------

     storkm_pp     = zero

     elements_4: DO iel = 1 , nels_pp
       bload     = zero
       dload     = zero   

       eld = pmul_pp(:,iel)
  
       sbary      = prop(1,etype_pp(iel))
       e          = prop(2,etype_pp(iel))
       v          = prop(3,etype_pp(iel))

!------------------------------------------------------------------------------
! 25. Go round the Gauss Points
!------------------------------------------------------------------------------
 
       gauss_points_2 : DO i = 1 , nip
          elso    = zero
          CALL shape_der(der,points,i)
          jac     = MATMUL(der,g_coord_pp(:,:,iel))
          det     = determinant(jac)
          CALL invert(jac)
          deriv   = MATMUL(jac,der)
          CALL beemat (deriv,bee)
          eps     = MATMUL(bee,eld)
          CALL deemat(e,v,dee)
          stressv = tensor_pp(:,i,iel)
          CALL invar(stressv,sigm,dsbar,lode_theta) 

          ff = dsbar - sbary                           
          IF(ff>fftol) THEN
            dlam   = dl(i,iel) 
            CALL vmflow(stressv,dsbar,vmfl)
            CALL fmrmat(vmfl,nst,dsbar,dlam,dee,temp,rmat)
            caflow = MATMUL(rmat,vmfl)
            bot    = DOT_PRODUCT(vmfl,caflow)
            CALL formaa(vmfl,nst,rmat,daatd)
            dee    = rmat - daatd/bot
          END IF
        
          sigma    = MATMUL(dee,eps)
          stressv  = sigma + tensor_pp( : , i , iel)

          CALL invar(stressv,sigm,dsbar,lode_theta)                            

!------------------------------------------------------------------------------
! 26. Check whether yield is violated
!------------------------------------------------------------------------------
 
          fnew   = dsbar - sbary  
          fstiff = fnew

          IF (fnew>=zero) THEN
            CALL deemat(e,v,dee) ! Fran's version needs checking 
            CALL vmflow(stressv,dsbar,vmfl)
            caflow    = MATMUL(dee,vmfl)
            bot       = DOT_PRODUCT(vmfl,caflow)
            dlam      = fnew/bot
            elso      = caflow*dlam 
            stressv   = tensor_pp( : , i , iel) + sigma - elso
            CALL invar(stressv,sigm,dsbar,lode_theta)
            fnew      = dsbar-sbary
            
            iterate_on_fnew : DO
            
              CALL vmflow(stressv,dsbar,vmfl)
              caflow = MATMUL(dee,vmfl)*dlam
              ress   = stressv - (tensor_pp(: , i , iel) +sigma - caflow)
              CALL fmacat(vmfl,nst,temp,acat)
              acat   = acat / dsbar
              acatc  = MATMUL(dee,acat)
              qmat   = acatc*dlam
            
              DO k = 1, nst
                qmat(k,k) = qmat(k,k) + one
              END DO
              
              CALL invert(qmat)
              vmtemp(1,:) = vmfl
              vmtemp      = MATMUL(vmtemp,qmat)
              vmflq       = vmtemp(1,:)
              top         = DOT_PRODUCT(vmflq,ress)
              vmtemp      = MATMUL(vmtemp,dee)
              vmfla       = vmtemp(1,:) 
              bot         = DOT_PRODUCT(vmfla,vmfl) 
              dslam       = (fnew - top)/bot
              qinvr       = MATMUL(qmat,ress)
              qinva       = MATMUL(matmul(qmat,dee),vmfl)
              dsigma      = -qinvr - qinva*dslam
              stressv     = stressv + dsigma
              CALL invar(stressv,sigm,dsbar,lode_theta)
              fnew        = dsbar - sbary
              dlam        = dlam + dslam

              IF (fnew<ltol) EXIT

            END DO iterate_on_fnew

            dl(i,iel) = dlam
            elso      = tensor_pp( : , i , iel) + sigma - stressv
            eload     = MATMUL(elso,bee)
            bload     = bload+eload*det*weights(i)
             
            CALL vmflow(stressv,dsbar,vmfl)
            CALL fmrmat(vmfl,nst,dsbar,dlam,dee,temp,rmat)
             
            caflow    = MATMUL(rmat,vmfl)
            bot       = DOT_PRODUCT(vmfl,caflow)
             
            CALL formaa(vmfl,nst,rmat,daatd)
          
            dee       = rmat - daatd/bot
          
          END IF
          
          IF(fstiff<zero) CALL deemat(e,v,dee)
          
          storkm_pp(:,:,iel) = storkm_pp(:,:,iel) +                           &
                               MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)         &
                               *det* weights(i)

!------------------------------------------------------------------------------
! 27. Update the Gauss Point stresses
!------------------------------------------------------------------------------

          tensor_pp( : , i , iel) = tensor_pp( : , i , iel) + sigma - elso
          stressv                 = tensor_pp( : , i , iel)
          eload                   = matmul(stressv,bee)
          dload                   = dload+eload*det*weights(i)
          
    END DO gauss_points_2

!------------------------------------------------------------------------------
! 28. Compute the total bodyloads vector
!------------------------------------------------------------------------------
 
    utemp_pp(:,iel) = utemp_pp(:,iel) + bload
    dtemp_pp(:,iel) = dtemp_pp(:,iel) + dload
 
  END DO elements_4 

!------------------------------------------------------------------------------
! 29. Update the preconditioner
!------------------------------------------------------------------------------
 
  ALLOCATE(diag_precon_tmp(ntot,nels_pp))
  
  diag_precon_pp     = zero
  diag_precon_tmp    = zero
  
  elements_5: DO iel = 1,nels_pp 
    DO k=1,ntot
      diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel) + storkm_pp(k,k,iel)
    END DO
  END DO elements_5
  
  CALL scatter(diag_precon_pp,diag_precon_tmp)
  
  DEALLOCATE(diag_precon_tmp)
  
  IF(fixed_freedoms_pp > 0) THEN
    DO i = 1,fixed_freedoms_pp
       j = no_pp(i) - ieq_start + 1
       diag_precon_pp(j) = diag_precon_pp(j) + penalty
       store_pp(i)       = diag_precon_pp(j)
    END DO
  END IF

  diag_precon_pp     = one/diag_precon_pp

!------------------------------------------------------------------------------
! 30. Test for convergence in the plastic iteration loop. Load control is
!     based on M.A. Crisfield, "Non-linear Finite Element Analysis of Solids 
!     and Structures", Volume 1, Wiley, p39. Displacement control is based on
!     Smith and Griffiths Edition 4.
!------------------------------------------------------------------------------

  bdylds_pp = zero
  ddylds_pp = zero

  CALL scatter(bdylds_pp,utemp_pp) 
  CALL scatter(ddylds_pp,dtemp_pp)
 
  tload      = SUM_P(ddylds_pp)  

  IF(plasiters==1) plastic_converged = .FALSE.

  IF(displacement_control) THEN
    tloads   = ABS(SUM_P(bdylds_pp)) 
    IF(tloads < plastol) plastic_converged = .TRUE.
    IF(numpe==1) THEN
      WRITE(11,'(A,I12,A,E12.4,A,E12.4)')" ",plasiters," ",tload,              &
                                         " ", tloads
    END IF 
  ELSE ! load control
    ddyld_l2n  = NORM_P(ddylds_pp)  
    bdyld_l2n  = NORM_P(bdylds_pp)
    tolerance  = bdyld_l2n/ddyld_l2n
    IF(tolerance < plastol) plastic_converged = .TRUE.
    IF(numpe==1) THEN
      WRITE(11,'(A,I12,A,E12.4,A,E12.4,A,E12.4)')" ",plasiters," ",bdyld_l2n,  &
                                                 " ", ddyld_l2n, " ", tolerance
    END IF 
  END IF

  totd_pp     = totd_pp + loads_pp           

  IF(plastic_converged .OR. plasiters==plasitersMax) EXIT
                                      
 END DO plastic_iterations
 
!------------------------------------------------------------------------------
! 31. Automatic refinement of load increments
!------------------------------------------------------------------------------

 iy              = iy - 1
 fractionApplied = fractionApplied + factor ! fraction of load applied
 
! Test whether the previous increment has converged. If it did not, we need to 
! reduce the size of the load increment and reapply. If it did, but too quickly
! i.e. plasiters < plasitersMin, we need to increase the size of the next load
! increment. These two scenarios are covered in SUBROUTINE adaptiveStepSize

 IF(plasiters == plasitersMax) THEN
   ! restore vectors to their old values
   storkm_pp       = storkm_pp_old
   diag_precon_pp  = diag_precon_pp_old
   tensor_pp       = tensor_pp_old
   totd_pp         = totd_pp_old
   dl              = dl_old
   fractionApplied = fractionApplied - factor
   iy              = iy + 1
 END IF
 
!------------------------------------------------------------------------------
! 32. Output results
!------------------------------------------------------------------------------
 
    IF(plasiters < plasitersMax) THEN

      outputIncrement       = outputIncrement + 1
      qinc(outputIncrement) = factor

      IF(numpe==1) THEN
       
        CALL getFileNumber(increment, outputIncrement)

        fname_dis = job_name(1:INDEX(job_name, " ")-1)  //                    &
                  "_"  // step(1:INDEX(step, " ")-1)  //                      &
                  "_i" // increment(1:INDEX(increment, " ")-1) // ".dis" 

        OPEN(24, file=fname_dis, status='replace', action='write')

      END IF

      CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
      ALLOCATE(xnewel_pp(ntot,nels_pp))
      ALLOCATE(xnewnodes_pp(nodes_pp*ndim))
 
      xnewel_pp    = zero 
      xnewnodes_pp = zero
      text         = "*DISPLACEMENT"

      CALL gather(totd_pp(1:),xnewel_pp)
      CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,          &
                         node_start,node_end,xnewel_pp,xnewnodes_pp,1)
      CALL write_nodal_variable(text,24,outputIncrement,nodes_pp,npes,numpe,  &
                                ndim,xnewnodes_pp)

      DEALLOCATE(xnewel_pp)
      DEALLOCATE(xnewnodes_pp)

      IF(numpe==1) THEN
!       CALL FLUSH_(24)
        CLOSE(24)
      END IF
 
    END IF

!------------------------------------------------------------------------------
! 33. Check whether step size needs adjusting
!------------------------------------------------------------------------------
 
    IF(iy /= 0) THEN
      CALL adaptiveStepSize(iy, factor, plasiters, plasitersMin,              &
                         plasitersMax)
    END IF

!------------------------------------------------------------------------------
! 34. Print data useful for checking how the analysis progressed.
!------------------------------------------------------------------------------

    maxDisp = 0
    maxDisp = MAXABSVAL_P(totd_pp)
 
    IF(numpe==1) THEN
      WRITE(11,'(/A,E12.4)')"Time for load increment:                ",        &
                             elap_time() - tinc
      WRITE(11,'(A,E12.4)') "Maximum displacement:                   ",        &
                             maxDisp
      WRITE(11,'(A,E12.4)') "Load applied:                           ",        &
                             tload
      WRITE(11,'(A,E12.4)') "sigma x:                                ",        &
                             tensor_pp(1,1,1)
      WRITE(11,'(A,E12.4)') "sigma y:                                ",        &
                             tensor_pp(2,1,1)
      WRITE(11,'(A,E12.4)') "sigma z:                                ",        &
                             tensor_pp(3,1,1)
      WRITE(11,'(A,I12)')   "Number of cj iterations:                ", cjtot
      WRITE(11,'(A,I12)')   "Number of plastic iterations:           ",        &
                             plasiters
      WRITE(11,'(A,F12.4/)')"cj iters / plastic iteration:           ",        &
                          &  real(cjtot)/real(plasiters)
    END IF

    IF(iy == 0 .OR. factor < 0.00001) THEN
      IF(numpe==1) THEN
        WRITE(11,'(/2A)')     "----------------------------------------",     &
                              "------------"
        WRITE(11,'(A,F12.4)') "Fraction of load applied:               ",     &
                               fractionApplied

        IF(factor < 0.00001) THEN
          WRITE(11,'(/A)')                                                    &
          "NOTE: Program ended as last automatic increment was "
          WRITE(11,'(A/)')                                                    & 
          "      less than 0.001% of the total applied load."
        ELSE
          WRITE(11,'(A)')       "All load increments completed"
        END IF

        WRITE(11,'(A,E12.4)') "This analysis took:                     ",   &
                               elap_time( ) - timest(1)
        WRITE(11,'(2A)')      "----------------------------------------",     &
                              "------------"
        WRITE(11,'(2A/)')     "----------------------------------------",      &
                              "------------"
      END IF
      EXIT
    END IF

    IF(plasiters == plasitersMax) THEN
      IF(numpe==1) THEN
        WRITE(11,'(/2A)')   "----------------------------------------",        &
                            "------------"
        WRITE(11,'(A)')     "Load increment aborted"
        WRITE(11,'(A)')     "Restarting with smaller increments"
      END IF
    ELSE
      IF(numpe==1) THEN
        WRITE(11,'(/2A)')   "----------------------------------------",        &
                            "------------"
        WRITE(11,'(A)')     "Load increment completed"
        WRITE(11,'(A)')     "Starting next increment"
      END IF
    END IF

    IF(numpe==1) THEN
      WRITE(11,'(A,I12)')   "Load increments remaining:              ", iy
      WRITE(11,'(A,F12.4)') "Size of remaining load increments:      ", factor
      WRITE(11,'(A,F12.4)') "Fraction of load applied:               ",        &
                             fractionApplied
      WRITE(11,'(2A)')      "----------------------------------------",        &
                            "------------"
      WRITE(11,'(2A/)')     "----------------------------------------",        &
                            "------------"
    END IF

!   CALL FLUSH_(11)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 35. End load increment loop
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
  END DO load_increments

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 36. End load step loop
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
   IF(numpe==1) THEN
     WRITE(11,'(/A)') "If running the same analysis again, reduce the "
     WRITE(11,'(A/)') "solution time by using these step sizes:       "
     DO i=1,outputIncrement
       WRITE(11,'(A,I5,A,F12.8)') "Load step ", i, " ", qinc(i)
     END DO
     WRITE(11,'(/2A/)')     "----------------------------------------",        &
                            "------------"
   END IF
    
 END DO load_steps

!------------------------------------------------------------------------------
! 37. End program
!------------------------------------------------------------------------------

 CLOSE(11)
 CALL shutdown()
 

END PROGRAM p1211
