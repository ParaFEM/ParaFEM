PROGRAM xx12      
!------------------------------------------------------------------------------
!      Program XX.12 Three dimensional analysis of conduction equation using 
!                    8-node hexahedral elements; pcg version implicit;  
!                    integration in time using 'theta' method parallel version
!
!------------------------------------------------------------------------------
  
  USE precision  ; USE global_variables ; USE mp_interface
  USE input      ; USE output           ; USE loading
  USE timing     ; USE maths            ; USE gather_scatter
  USE partition  ; USE elements         ; USE steering
  USE geometry   ; USE pcg              ; USE new_library
  
  USE, INTRINSIC :: ISO_C_BINDING ! to output C binary file
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------
  
! neq,ntot are now global variables - not declared
  
  INTEGER, PARAMETER  :: ndim=3,nodof=1,nprops=5
  INTEGER             :: nod,nn,nr,nip
  INTEGER             :: j_glob,j_loc,j_step,j_step2,j_temp,j_chk,j_chk2
  INTEGER             :: j_npri,j_npri_chk
  INTEGER             :: i,j,k,l,limit,iel,red_blk
  INTEGER             :: nxe,nye,nze,neq_temp,nn_temp,npp
  INTEGER             :: nstep,nstep_tot,npri,npri_chk,nres,it,is,nlen
  INTEGER             :: node_end,node_start,nodes_pp
  INTEGER             :: loaded_freedoms,fixed_freedoms,loaded_nodes
  INTEGER             :: fixed_freedoms_pp,fixed_freedoms_start
  INTEGER             :: loaded_freedoms_pp,loaded_freedoms_start
  INTEGER             :: nels,ndof,ielpe,npes_pp
  INTEGER             :: argc,iargc,meshgen,partitioner
  INTEGER             :: np_types,ntime,el_print,i_o
  INTEGER             :: prog,tz
  REAL(iwp)           :: aa,bb,cc,kx,ky,kz,det,theta,dtim,real_time
  REAL(iwp)           :: tol,alpha,beta,up,big,q
  REAL(iwp)           :: rho,cp,val0
  REAL(iwp),PARAMETER :: zero = 0.0_iwp,penalty=1.e20_iwp
  REAL(iwp),PARAMETER :: t0 = 0.0_iwp
  CHARACTER(LEN=15)   :: element,chk
  CHARACTER(LEN=50)   :: fname,job_name,label,stepnum
  CHARACTER(LEN=50)   :: program_name='xx12'
!  character (len=512) command
  LOGICAL             :: converged = .false.
  LOGICAL             :: solid=.true.
  CHARACTER(LEN=80)   :: cbuffer
  
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------ 
  
  REAL(iwp),ALLOCATABLE :: loads_pp(:),u_pp(:),p_pp(:),points(:,:),kay(:,:)
  REAL(iwp),ALLOCATABLE :: coord(:,:),fun(:),jac(:,:),der(:,:),deriv(:,:)
  REAL(iwp),ALLOCATABLE :: weights(:),d_pp(:),kc(:,:),pm(:,:),funny(:,:)
  REAL(iwp),ALLOCATABLE :: p_g_co_pp(:,:,:),storka_pp(:,:,:),storkb_pp(:,:,:)
  REAL(iwp),ALLOCATABLE :: x_pp(:),xnew_pp(:),pmul_pp(:,:)
  REAL(iwp),ALLOCATABLE :: utemp_pp(:,:),diag_precon_pp(:),diag_precon_tmp(:,:)
  REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:),timest(:)
  REAL(iwp),ALLOCATABLE :: disp_pp(:),eld_pp(:,:)
  REAL(iwp),ALLOCATABLE :: val(:,:),val_f(:),store_pp(:),r_pp(:)
  REAL(iwp),ALLOCATABLE :: kcx(:,:),kcy(:,:),kcz(:,:)
  REAL(iwp),ALLOCATABLE :: eld(:),col(:,:),row(:,:),storkc_pp(:,:,:)
  REAL(iwp),ALLOCATABLE :: prop(:,:),amp(:),tempres(:),timesteps_real(:,:)
  INTEGER,ALLOCATABLE   :: rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp(:,:)
  INTEGER,ALLOCATABLE   :: no(:),no_pp(:),no_f_pp(:),no_pp_temp(:),sense(:)
  INTEGER,ALLOCATABLE   :: node(:),iters(:),iters_tot(:),timesteps_int(:,:)
  INTEGER,ALLOCATABLE   :: etype_pp(:),nf(:,:),oldlds(:),g_coord(:,:)

!------------------------------------------------------------------------------
! *. Details of some xx12 specific variables
!------------------------------------------------------------------------------ 
  
  !/*
  !*  INTEGERS
  !*    j_glob      : Global time step counter
  !*    j_loc       : Time step counter within time section
  !*    j_step      : Time section counter
  !*    j_step2     : Time section counter used for restarting from checkpoint
  !*    j_temp      : Temporary counter used to calculate current j_step
  !*    j_chk       : Global time step counter retrieved from checkpoint
  !*    j_chk2      : Logical flag to initialise if restarting from checkpoint
  !*    j_npri      : File numbering output for ENSIGHT GOLD format
  !*    j_npri_chk  : Counter to trigger checkpointing
  !*    nstep       : Number of timesteps for given time section
  !*    nstep_tot   : Total number of timesteps
  !*    npri        : Frequency of results output
  !*    npri_chk    : Frequency of checkpointing
  !*    ntime       : Number of time sections
  !*    el_print    : Element number output for quick verification of results
  !*    i_o         : Output format (ParaFEM binary, ascii or ENSIGHT binary)
  !*  
  !*  REALS
  !*    dtim        : Time step size
  !*    real_time   : Current time in solution
  !*  
  !*/
  
!------------------------------------------------------------------------------
! 3. Read job_name from the command line. 
!    Read control data, mesh data, boundary and loading conditions
!------------------------------------------------------------------------------ 
  
  ALLOCATE(timest(25))
  timest    = zero
  timest(1) = elap_time()
  
  CALL find_pe_procs(numpe,npes)
  
!  PRINT *, "FIND_PE_PROCS on processor ", numpe, " of ", npes
  
  argc = iargc()
  IF(argc /= 1) CALL job_name_error(numpe,program_name)
  CALL GETARG(1,job_name)
  
  CALL read_xx12(job_name,numpe,element,fixed_freedoms,limit,                 &
                 loaded_nodes,meshgen,nels,nip,nn,nod,ntime,nr,               &
                 partitioner,theta,tol,np_types,chk,val0,el_print,i_o)
  
  CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)
  
  ndof = nod*nodof
  ntot = ndof
  
  ALLOCATE(g_num_pp(nod,nels_pp))
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp))
  IF (nr>0) ALLOCATE(rest(nr,nodof+1))
  ALLOCATE(etype_pp(nels_pp))
  ALLOCATE(prop(nprops,np_types))
  ALLOCATE(timesteps_real(ntime,1),timesteps_int(ntime,3))
  
  g_num_pp       = 0
  g_coord_pp     = zero
  IF (nr>0) rest = 0
  etype_pp       = 0
  prop           = zero
  q              = zero
  
  timest(2) = elap_time()
  
! CALL read_g_num_pp2(job_name,iel_start,nn,npes,numpe,g_num_pp)
!  CALL read_elements(job_name,iel_start,nn,npes,numpe,etype_pp,g_num_pp)
  IF(numpe==1) PRINT *, "CALL read_g_num_pp_be"
  CALL read_g_num_pp_be(job_name,iel_start,nn,npes,numpe,g_num_pp)
  IF(numpe==1) PRINT *, "CALL read_etype_pp"
!  CALL read_etype_pp(job_name,npes,numpe,etype_pp)
  CALL read_etype_pp_be(job_name,npes,numpe,etype_pp)
  timest(3) = elap_time()
  
  IF(meshgen == 2)THEN
    IF(numpe==1) PRINT *,"Calling abaqus2sg, meshgen = ",meshgen
    CALL abaqus2sg(element,g_num_pp)
  END IF
  timest(4) = elap_time()
  
!  CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
  CALL read_g_coord_pp_be(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
  timest(5) = elap_time()
  
  IF (nr>0) CALL read_rest(job_name,numpe,rest)
  timest(6) = elap_time()
  
  IF(numpe==1) PRINT *, "np_types = ", np_types
  
  fname = job_name(1:INDEX(job_name, " ")-1) // ".mat"  
  CALL read_materialValue(prop,fname,numpe,npes)

  fname = job_name(1:INDEX(job_name, " ")-1) // ".time"  
  CALL read_timesteps(timesteps_real,timesteps_int,fname,numpe,npes)
  
  nstep_tot = 0
  npp = 0
  DO i=1,ntime
    nstep_tot = nstep_tot + timesteps_int(i,1)
    npp = npp + timesteps_int(i,1)/timesteps_int(i,2)
  END DO
  IF(numpe==1) PRINT *,"nstep_tot = ",nstep_tot
  
  IF(numpe==1) PRINT *, " *** Read input data in: ", timest(6)-timest(1)," s"
  
!------------------------------------------------------------------------------
! 4. Allocate dynamic arrays used in main program
!------------------------------------------------------------------------------
  
  ALLOCATE (points(nip,ndim),weights(nip),kay(ndim,ndim),                     &
            coord(nod,ndim),fun(nod),jac(ndim,ndim),der(ndim,nod),g(ntot),    &
            deriv(ndim,nod),pm(ntot,ntot),                                    &
            kc(ntot,ntot),funny(1,nod),num(nod),                              &
            g_g_pp(ntot,nels_pp),storka_pp(ntot,ntot,nels_pp),                &
            utemp_pp(ntot,nels_pp),storkb_pp(ntot,ntot,nels_pp),              &
            pmul_pp(ntot,nels_pp),iters(ntime),iters_tot(ntime))
  ALLOCATE (kcx(ntot,ntot),kcy(ntot,ntot),kcz(ntot,ntot),                     &
            eld(ntot),col(ntot,1),row(1,ntot),storkc_pp(ntot,ntot,nels_pp))
  ALLOCATE (amp(nstep_tot))
  
  IF(numpe==1) PRINT *, " *** Allocated dynamic arrays in: ",                 &
                          elap_time()-timest(6)," s"
  
!------------------------------------------------------------------------------
! 5. Loop the elements to find the steering array and the number of equations
!    to solve.
!------------------------------------------------------------------------------
  
  IF (nr>0) CALL rearrange_2(rest)  
  
  g_g_pp = 0
  
  ! When nr = 0, g_num_pp and g_g_pp are identical
  IF(nr>0) THEN
    elements_1: DO iel = 1, nels_pp
      CALL find_g4(g_num_pp(:,iel),g_g_pp(:,iel),rest)
    END DO elements_1
    DEALLOCATE(rest)
  ELSE
    g_g_pp = g_num_pp
  END IF
  
  neq = 0
  
  elements_2: DO iel = 1, nels_pp
    i = MAXVAL(g_g_pp(:,iel))
    IF(i > neq) neq = i
  END DO elements_2
  
  neq = MAX_INTEGER_P(neq)
  
  timest(7) = elap_time()
  IF(numpe==1) PRINT *, "End of 5"
  
!------------------------------------------------------------------------------
! 6. Create interprocessor communication tables
!------------------------------------------------------------------------------
  
  CALL calc_neq_pp          
  CALL calc_npes_pp(npes,npes_pp)
  CALL make_ggl(npes_pp,npes,g_g_pp)
  
  nres = el_print
  
  DO i = 1,neq_pp
    IF(nres==ieq_start+i-1) THEN
      it = numpe; is = i
      IF(numpe==it) PRINT *, " *** it = ", it, " is = ", i
    END IF
  END DO
  
  timest(8) = elap_time()
  IF(numpe==1) PRINT *, "End of 6"
  
!------------------------------------------------------------------------------
! 7. Allocate arrays dimensioned by neq_pp 
!------------------------------------------------------------------------------
  
  ALLOCATE(loads_pp(neq_pp),diag_precon_pp(neq_pp),u_pp(neq_pp),d_pp(neq_pp), &
           p_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),r_pp(neq_pp))
  
  loads_pp  = zero ; diag_precon_pp = zero ; u_pp = zero ; r_pp    = zero
  d_pp      = zero ; p_pp           = zero ; x_pp = zero ; xnew_pp = zero
  
  timest(9) = elap_time()
  IF(numpe==1) PRINT *, "End of 7"
  
!------------------------------------------------------------------------------
! 8. Element stiffness integration and storage
!------------------------------------------------------------------------------
  
  CALL sample(element,points,weights)
  
  red_blk=1
  SELECT CASE (chk)
    CASE('start')
      IF(numpe==1) PRINT*, "Initialising new run"
      j_chk   = 0
      j_glob  = 1
      j_step2 = 1
      j_chk2  = 0
    CASE('restart')
      IF(numpe==1) PRINT*, "Restarting from checkpoint"
      j_chk2  = 1
      CALL read_x_pp(job_name,npes,numpe,j_chk,x_pp)
      j_glob = j_chk+1
      j_temp = 0
      DO i=1,ntime
        j_temp  = j_temp + timesteps_int(i,1)
        j_step2 = i
        IF(j_temp>j_glob-1)EXIT 
      END DO
      IF(j_step2>1) j_chk = (j_glob-1) - (j_temp-timesteps_int(j_step2,1))
    CASE default
      IF(numpe==1) PRINT*, "Invalid checkpoint flag in .dat"
      IF(numpe==1) PRINT*, "Program aborting"
  END SELECT
  
!------------------------------------------------------------------------------
! *. Start timesection loop
!------------------------------------------------------------------------------ 

  iters_tot = 0
  iters = 0
  
  timesections: DO j_step=j_step2,ntime
    j_npri_chk = 0
    dtim       = timesteps_real(j_step,1)
    nstep      = timesteps_int(j_step,1)
    npri       = timesteps_int(j_step,2)
    npri_chk   = timesteps_int(j_step,3)
    
    storka_pp = zero 
    storkb_pp = zero
    
    elements_3: DO iel=1,nels_pp
      
      kay       = zero
      kay(1,1)  = prop(1,etype_pp(iel))  ! kx
      kay(2,2)  = prop(2,etype_pp(iel))  ! ky
      kay(3,3)  = prop(3,etype_pp(iel))  ! kz
      rho       = prop(4,etype_pp(iel))  ! rho
      cp        = prop(5,etype_pp(iel))  ! cp
      
      kc = zero ; pm = zero
      
      gauss_pts: DO i=1,nip
        CALL shape_der(der,points,i)
        CALL shape_fun(fun,points,i)
        funny(1,:) = fun(:)
        jac        = MATMUL(der,g_coord_pp(:,:,iel))
        det        = determinant(jac)
        CALL invert(jac)
        deriv      = MATMUL(jac,der)
        kc         = kc + MATMUL(MATMUL(TRANSPOSE(deriv),kay),deriv)*         &
                                                                 det*weights(i)
        pm         = pm + MATMUL(TRANSPOSE(funny),funny)*det*weights(i)*rho*cp
!        !----------------------------------------------------------------------
!        !- Added to check if storekc_pp would be the same as xx11
!        row(1,:) = deriv(1,:)
!        eld      = deriv(1,:)
!        col(:,1) = eld
!        kcx      = kcx + MATMUL(col,row)*det*weights(i)
!        row(1,:) = deriv(2,:)
!        eld      = deriv(2,:)
!        col(:,1) = eld
!        kcy      = kcy + MATMUL(col,row)*det*weights(i)
!        row(1,:) = deriv(3,:)
!        eld      = deriv(3,:)
!        col(:,1) = eld
!        kcz      = kcz + MATMUL(col,row)*det*weights(i)
!        !----------------------------------------------------------------------
      END DO gauss_pts
      
      storka_pp(:,:,iel)=pm+kc*theta*dtim
      storkb_pp(:,:,iel)=pm-kc*(1._iwp-theta)*dtim
!      storkc_pp(:,:,iel) = kcx*kx + kcy*ky + kcz*kz
!      
!      !-- To run in xx11 mode set prog=11, in xx12 prog=12   
!      prog=12
!      IF(prog==11) storka_pp = storkc_pp
      
    END DO elements_3
    
    timest(10) = elap_time()
    IF(numpe==1) PRINT *, "End of 8"
    
!------------------------------------------------------------------------------
! 9. Build the diagonal preconditioner
!------------------------------------------------------------------------------
    
    ALLOCATE(diag_precon_tmp(ntot,nels_pp))
    diag_precon_tmp = zero
    diag_precon_pp  = zero
    
    elements_4: DO iel = 1,nels_pp 
      DO k=1,ntot
        diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel)+storka_pp(k,k,iel)
      END DO
    END DO elements_4
    
    CALL scatter(diag_precon_pp,diag_precon_tmp)
    
    DEALLOCATE(diag_precon_tmp)
    
    timest(11) = elap_time()
    IF(numpe==1) PRINT *, "End of 9" 
    
!------------------------------------------------------------------------------
! 10. Allocate disp_pp array and open file to write temperature output
!------------------------------------------------------------------------------
    
    IF(j_step==1 .OR. j_chk2==1)THEN
      IF(numpe==it)THEN !---Open file-temp outputs-specified node  
        fname = job_name(1:INDEX(job_name, " ")-1) // ".ttr2"
        OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
      END IF
      
      CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
      ALLOCATE(disp_pp(nodes_pp))
      ALLOCATE(eld_pp(ntot,nels_pp))
      
      IF(numpe==1) THEN
        IF(i_o==2)THEN !---Open file-temp outputs-ascii ParaFEM
          fname   = job_name(1:INDEX(job_name, " ")-1)//".ttr"
          OPEN(24, file=fname, status='replace', action='write')
        END IF
        IF(i_o==1)THEN !---Open file-tempoutputs-binary ParaFEM
          fname   = job_name(1:INDEX(job_name, " ")-1)//".ttrb"
          OPEN(25, file=fname, status='replace', action='write',              &
               access='sequential', form='unformatted')
        END IF
        fname   = job_name(1:INDEX(job_name, " ")-1)//".npp"
        OPEN(26, file=fname, status='replace', action='write')
        label   = "*TEMPERATURE"  
        WRITE(26,*) nn
!        WRITE(26,*) nstep/npri
        !-Unchecked
        WRITE(26,*) npp
        WRITE(26,*) npes
        IF(i_o==3.AND.j_glob==1)THEN!---Open file-temp out-binary ENSIGHT GOLD
          fname = job_name(1:INDEX(job_name, " ")-1)//".bin.ensi.NDTTR-000001"
          OPEN(27,file=fname,status='replace',action='write',                 &
          form='unformatted',access='stream')
           
          cbuffer="Alya Ensight Gold --- Scalar per-node variable file"
          WRITE(27) cbuffer
          cbuffer="part"        ; WRITE(27) cbuffer
          WRITE(27) int(1,kind=c_int)
          cbuffer="coordinates" ; WRITE(27) cbuffer 
        END IF
      END IF
    END IF
    
    IF(numpe==1) PRINT *, "End of 10"
    
!------------------------------------------------------------------------------
! 11. Read in fixed nodal temperatures and assign to equations
!------------------------------------------------------------------------------
  
    IF((fixed_freedoms > 0 .AND. j_glob == 1) .OR.                            &
       (fixed_freedoms > 0 .AND. j_chk2 == 1)) THEN
      
      ALLOCATE(node(fixed_freedoms),no(fixed_freedoms),                       &
               no_pp_temp(fixed_freedoms),sense(fixed_freedoms))
      ALLOCATE(val_f(fixed_freedoms))
      
      node  = 0 ; no = 0 ; no_pp_temp = 0 ; sense = 0
      val_f = zero
      
      CALL read_fixed(job_name,numpe,node,sense,val_f)
      CALL find_no2(g_g_pp,g_num_pp,node,sense,no)
      
!      PRINT *, "After find_no2 no = ", no, " on PE ", numpe
      
      CALL reindex(ieq_start,no,no_pp_temp,                                   &
                   fixed_freedoms_pp,fixed_freedoms_start,neq_pp)
      
!      PRINT *, "After reindex no = ", no, " on PE ", numpe
      
      ALLOCATE(no_f_pp(fixed_freedoms_pp),store_pp(fixed_freedoms_pp))
      
      no_f_pp  = 0 
      store_pp = zero
      no_f_pp  = no_pp_temp(1:fixed_freedoms_pp)
      
      DEALLOCATE(node,no,sense,no_pp_temp)
      
    END IF
    
    IF(fixed_freedoms == 0) fixed_freedoms_pp = 0
    
    timest(12) = elap_time()
    IF(numpe==1) PRINT *, "End of 11"
    
!------------------------------------------------------------------------------
! 12. Invert the preconditioner. 
!     If there are fixed freedoms, first apply a penalty
!------------------------------------------------------------------------------
    
    IF(fixed_freedoms_pp > 0) THEN
      DO i = 1,fixed_freedoms_pp
        l =  no_f_pp(i) - ieq_start + 1
        diag_precon_pp(l) = diag_precon_pp(l) + penalty
        store_pp(i)       = diag_precon_pp(l)
      END DO
    END IF
    
    diag_precon_pp = 1._iwp/diag_precon_pp
    
    IF(numpe==1) PRINT *, "End of 12"
    
!------------------------------------------------------------------------------
! 13. Read in loaded nodes and get starting r_pp
!------------------------------------------------------------------------------
      
      loaded_freedoms = loaded_nodes ! hack
      IF(loaded_freedoms==0) loaded_freedoms_pp=0
      IF((loaded_freedoms > 0 .AND. j_glob == 1) .OR.                         &
         (loaded_freedoms > 0 .AND. j_chk2 == 1)) THEN
        
        ALLOCATE(node(loaded_freedoms),val(nodof,loaded_freedoms))
        ALLOCATE(no_pp_temp(loaded_freedoms))
        
        val = zero ; node = 0
        
        CALL read_amplitude(job_name,numpe,nstep_tot,amp)
        CALL read_loads(job_name,numpe,node,val)
        CALL reindex(ieq_start,node,no_pp_temp,loaded_freedoms_pp,            &
                     loaded_freedoms_start,neq_pp)
        
        ALLOCATE(no_pp(loaded_freedoms_pp))
        
        no_pp    = no_pp_temp(1:loaded_freedoms_pp)
        
        DEALLOCATE(no_pp_temp)
        DEALLOCATE(node)
      
      END IF
    
    timest(12) = elap_time()
    IF(numpe==1) PRINT *, "End of 13"
    
!------------------------------------------------------------------------------
! 14. Start time stepping loop
!------------------------------------------------------------------------------
    
!    iters_tot = 0
    IF(j_chk2==0)j_chk=0
    timesteps: DO j_loc=j_chk+1,nstep
      
      real_time = 0
      j_npri = 0
      IF(j_step>1)THEN
        DO i = 1,j_step-1
          real_time = real_time + timesteps_real(i,1)*timesteps_int(i,1)
          j_npri = j_npri + timesteps_int(i,1)/timesteps_int(i,2)
        END DO
      END IF
      real_time = real_time + j_loc*dtim
      j_npri = j_npri + j_loc/npri + 1
      
      timest(15) = elap_time()
      IF(numpe==1) PRINT *, "End of 14"
      
!------------------------------------------------------------------------------
! 15. Apply loads (sources and/or sinks) supplied as a boundary value
!------------------------------------------------------------------------------
      
      loads_pp  = zero
      
      IF(loaded_freedoms_pp > 0) THEN
        DO i = 1, loaded_freedoms_pp
          IF(amp(j_loc)==0.0)THEN
            loads_pp(no_pp(i)-ieq_start+1) = val(loaded_freedoms_start+i-1,1) &
                                             *dtim*(1.0E-34)
          ELSE
            loads_pp(no_pp(i)-ieq_start+1) = val(loaded_freedoms_start+i-1,1) &
                                             *dtim*amp(j_loc)
          END IF
        END DO
      END IF
      
      q = q + SUM_P(loads_pp)

      IF(numpe==1) PRINT *, "End of 15"
      
!------------------------------------------------------------------------------
! 16. Compute RHS of time stepping equation, using storkb_pp, then add 
!     result to loads
!------------------------------------------------------------------------------
      
      u_pp              = zero
      pmul_pp           = zero
      utemp_pp          = zero
      
      IF(j_glob/=1) THEN
        
        IF(fixed_freedoms_pp > 0) THEN
          DO i = 1, fixed_freedoms_pp
            l       = no_f_pp(i) - ieq_start + 1
            k       = fixed_freedoms_start + i - 1
            x_pp(l) = val_f(k)
          END DO
        END IF
        
        CALL gather(x_pp,pmul_pp)
        elements_2a: DO iel=1,nels_pp
          utemp_pp(:,iel)=MATMUL(storkb_pp(:,:,iel),pmul_pp(:,iel))
        END DO elements_2a
        CALL scatter(u_pp,utemp_pp)
        
        IF(fixed_freedoms_pp > 0) THEN
          DO i = 1, fixed_freedoms_pp
            l       = no_f_pp(i) - ieq_start + 1
            k       = fixed_freedoms_start + i - 1
            u_pp(l) = store_pp(i)*val_f(k)
          END DO
        END IF
        
        loads_pp = loads_pp + u_pp
        
        IF(numpe==1) PRINT *, "End of 16"
        
      ELSE
        
!------------------------------------------------------------------------------
! 17. Set initial temperature
!------------------------------------------------------------------------------
        
        x_pp = val0
        
        IF(fixed_freedoms_pp > 0) THEN
          DO i = 1, fixed_freedoms_pp
            l       = no_f_pp(i) - ieq_start + 1
            k       = fixed_freedoms_start + i - 1
            x_pp(l) = val_f(k)
          END DO
        END IF
        
        CALL gather(x_pp,pmul_pp)
        elements_2c: DO iel=1,nels_pp
          utemp_pp(:,iel)=MATMUL(storka_pp(:,:,iel),pmul_pp(:,iel))
        END DO elements_2c
        CALL scatter(u_pp,utemp_pp)

! Added lines - seems to fix error at fixed nodes when t=dt
! However, only works when simulation also contains loaded nodes
! If fixed nodes only, PCG iterates until iterlimit
!        IF(fixed_freedoms_pp > 0) THEN
!          DO i = 1, fixed_freedoms_pp
!            l       = no_f_pp(i) - ieq_start + 1
!            k       = fixed_freedoms_start + i - 1
!            u_pp(l) = store_pp(i)*val_f(k)
!          END DO
!        END IF
        
        loads_pp = loads_pp + u_pp
        
        IF(numpe==1) PRINT *, "End of 17"
        
!------------------------------------------------------------------------------
! 18. Output "results" at t=0
!------------------------------------------------------------------------------
        
        eld_pp   = zero
        disp_pp  = zero
        tz       = 0
        CALL gather(x_pp(1:),eld_pp)
        
        CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,       &
                           node_start,node_end,eld_pp,disp_pp,1)
        
        IF(numpe==it)THEN
          !---Write file-temp outputs-specified node
          WRITE(11,'(E12.4,8E19.8)')t0,disp_pp(is)
        END IF
        
        IF(i_o==1)THEN !---Write file-temp outputs-binary ParaFEM
          CALL write_nodal_variable_binary(label,25,tz,nodes_pp,npes,         &
                                           numpe,nodof,disp_pp)
        END IF
        IF(i_o==2)THEN !---Write file-temp outputs-ascii ParaFEM
          CALL write_nodal_variable2(label,24,tz,nodes_pp,npes,               &
                                     numpe,nodof,disp_pp)
        END IF
        
        IF(i_o==3)THEN !---Write file-temp outputs-binary ENSIGHT GOLD
          ALLOCATE(tempres(nodes_pp))
          tempres = zero
          DO l=1,nodes_pp
            tempres(l)=disp_pp(l)
          END DO
          CALL dismsh_ensi_pb2(27,j_glob,nodes_pp,npes,numpe,nodof,tempres)
          
          DEALLOCATE(tempres)
          IF(numpe==1) CLOSE (27)  
        END IF
        
      END IF ! From section 16
      
      IF(numpe==1) PRINT *, "End of 18"
      
!------------------------------------------------------------------------------
! 19. Initialize PCG process
! 
!     When x = 0._iwp p and r are just loads but in general p=r=loads-A*x,
!     so form r = A*x. Here, use LHS part of the transient equation storka_pp
!------------------------------------------------------------------------------
      r_pp              = zero
      pmul_pp           = zero
      utemp_pp          = zero
!     x_pp              = zero
      
      CALL gather(x_pp,pmul_pp)
      elements_2b: DO iel=1,nels_pp
        utemp_pp(:,iel)=MATMUL(storka_pp(:,:,iel),pmul_pp(:,iel))
      END DO elements_2b
      CALL scatter(r_pp,utemp_pp)
      
      IF(fixed_freedoms_pp > 0) THEN
        DO i = 1, fixed_freedoms_pp
          l       = no_f_pp(i) - ieq_start + 1
          k       = fixed_freedoms_start + i - 1
          r_pp(l) = store_pp(i)*val_f(k)
        END DO
      END IF
      
      r_pp = loads_pp - r_pp
      d_pp = diag_precon_pp*r_pp
      p_pp = d_pp
      
      IF(numpe==1) PRINT *, "End of 19"
!------------------------------------------------------------------------------
! 20. Solve simultaneous equations by pcg
!------------------------------------------------------------------------------
      
      iters(j_step) = 0
      
      iterations: DO
        
        iters(j_step) = iters(j_step)+1
        
        u_pp     = zero
        pmul_pp  = zero
        utemp_pp = zero
        
        CALL gather(p_pp,pmul_pp)
        elements_6: DO iel=1,nels_pp
          utemp_pp(:,iel)=MATMUL(storka_pp(:,:,iel),pmul_pp(:,iel))
        END DO elements_6
        CALL scatter(u_pp,utemp_pp)
        
!        IF(numpe==1) PRINT *, "End of 20"
        
!------------------------------------------------------------------------------
! 21. PCG equation solution
!------------------------------------------------------------------------------
        
        IF(fixed_freedoms_pp > 0) THEN
          DO i = 1, fixed_freedoms_pp
            l       = no_f_pp(i) - ieq_start + 1
            u_pp(l) = p_pp(l) * store_pp(i)
          END DO
        END IF
        
        up       = DOT_PRODUCT_P(r_pp,d_pp)
        alpha    = up/DOT_PRODUCT_P(p_pp,u_pp)
        xnew_pp  = x_pp+p_pp*alpha
        r_pp     = r_pp-u_pp*alpha
        d_pp     = diag_precon_pp*r_pp
        beta     = DOT_PRODUCT_P(r_pp,d_pp)/up
        p_pp     = d_pp+p_pp*beta
        
        CALL checon_par(xnew_pp,tol,converged,x_pp)
        IF(converged.OR.iters(j_step)==limit)EXIT
        
      END DO iterations
      
      timest(13) = timest(13) + (elap_time() - timest(15))
      timest(16) = elap_time()
      
      IF(numpe==1) PRINT *, "End of 21"
      
!------------------------------------------------------------------------------
! 22. Output results and checkpoint
!------------------------------------------------------------------------------
      
      IF(j_loc/npri*npri==j_loc)THEN
      j_npri_chk=j_npri_chk+1
        
        eld_pp   = zero
        disp_pp  = zero
        CALL gather(xnew_pp(1:),eld_pp)
        
        CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,       &
                           node_start,node_end,eld_pp,disp_pp,1)
        
        IF(numpe==it)THEN
          !---Write file-temp outputs-specified node
          WRITE(11,'(E12.4,8E19.8)')real_time,disp_pp(is)
        END IF      
        
        IF(i_o==1)THEN !---Write file-temp outputs-binary ParaFEM
          CALL write_nodal_variable_binary(label,25,j_glob,nodes_pp,npes,     &
                                           numpe,nodof,disp_pp)
        END IF
        IF(i_o==2)THEN !---Write file-temp outputs-ascii ParaFEM
          CALL write_nodal_variable2(label,24,j_glob,nodes_pp,npes,           &
                                     numpe,nodof,disp_pp)
        END IF
        
        IF(i_o==3)THEN !---Write file-temp outputs-binary ENSIGHT GOLD
          IF(numpe==1)THEN
            WRITE(stepnum,'(I0.6)') j_npri
            fname=job_name(1:INDEX(job_name," ")-1)//".bin.ensi.NDTTR-"       &
                  //stepnum
            OPEN(27,file=fname,status='replace',action='write',               &
            form='unformatted',access='stream')
          
            cbuffer="Alya Ensight Gold --- Scalar per-node variable file"
            WRITE(27) cbuffer
            cbuffer="part"        ; WRITE(27) cbuffer
            WRITE(27) int(1,kind=c_int)
            cbuffer="coordinates" ; WRITE(27) cbuffer 
          END IF
          
          ALLOCATE(tempres(nodes_pp))
          tempres = zero
          DO l=1,nodes_pp
            tempres(l)=disp_pp(l)
          END DO
          !-Is tempres needed? Would passing disp_pp directly work the same?
          CALL dismsh_ensi_pb2(27,j_glob,nodes_pp,npes,numpe,nodof,tempres)
          DEALLOCATE(tempres)
          IF(numpe==1) CLOSE (27)
        END IF
        
        IF(numpe==1) PRINT *, "Time ", real_time, "Iters ", iters(j_step)
        
        !--Write checkpoint file
        IF(j_npri_chk==npri_chk)THEN
          j_npri_chk=0
          IF(numpe==1 .AND. red_blk==1)                                       &
             PRINT *, "Checkpoint: j =",j_glob,", red_blk = blk"
          IF(numpe==1 .AND. red_blk==-1)                                      &
             PRINT *, "Checkpoint: j =",j_glob,", red_blk = red"
          IF(red_blk==1)THEN
            IF(numpe==1)THEN
              fname = job_name(1:INDEX(job_name, " ")-1)//".chk_blk"
              OPEN(28,file=fname,status='replace',action='write',             &
                   form='unformatted',access='stream')
              WRITE(28) int(j_glob,kind=c_int)
            END IF
          END IF
          IF(red_blk==-1)THEN
            IF(numpe==1)THEN
              fname = job_name(1:INDEX(job_name, " ")-1)//".chk_red"
              OPEN(28,file=fname,status='replace',action='write',             &
                   form='unformatted',access='stream')
              WRITE(28) int(j_glob,kind=c_int)
            END IF
          END IF
          CALL write_x_pp(label,28,j_glob,nodes_pp,npes,numpe,nodof,x_pp)
          IF(numpe==1) CLOSE (28)
          red_blk=red_blk*(-1)
        END IF
      END IF
      
      timest(14) = timest(14) + (elap_time() - timest(16))
      iters_tot(j_step) = iters_tot(j_step) + iters(j_step)
      
      j_chk2  = 0
      j_glob=j_glob+1
    END DO timesteps
  END DO timesections

  timest(13) = timest(12) + timest(13)
  timest(14) = timest(13) + timest(14)
  
  IF(numpe==1) PRINT *, "End of 22"
  
  IF(numpe==1)THEN
    CLOSE(11)
    CLOSE(24)
  END IF
  
  IF(numpe==1) PRINT *, "Timest ", timest
  
  CALL WRITE_XX12(fixed_freedoms,job_name,loaded_freedoms,neq,nn,npes,nr,     &
                  numpe,timest,q,iters,iters_tot,tol,val0,ntime,              &
                  timesteps_int,timesteps_real)
  
  CALL shutdown()
  
END PROGRAM xx12
