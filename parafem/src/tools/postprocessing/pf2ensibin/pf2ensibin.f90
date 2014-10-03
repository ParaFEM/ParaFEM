PROGRAM pf2ensibin      
!------------------------------------------------------------------------------
!      Program pf2ensibin Tool to produce binary ensight gold formatted  
!                         input files.  
!                         Under development, currently for 4-node tetrahedra
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
  
  INTEGER, PARAMETER  :: ndim=3,nodof=1,nprops=5
  INTEGER             :: nod,nn,nr,nip
  INTEGER             :: i,j,k,l,m,n,iters,iters_tot,limit,iel
  INTEGER             :: nxe,nye,nze,neq_temp,nn_temp
  INTEGER             :: nstep,npri,nres,it,is,nlen
  INTEGER             :: node_end,node_start,nodes_pp
  INTEGER             :: loaded_freedoms,fixed_freedoms,loaded_nodes
  INTEGER             :: fixed_freedoms_pp,fixed_freedoms_start
  INTEGER             :: loaded_freedoms_pp,loaded_freedoms_start
  INTEGER             :: nels,ndof,ielpe,npes_pp
  INTEGER             :: argc,iargc,meshgen,partitioner
  INTEGER             :: np_types,el_print,i_o
  INTEGER             :: prog,tz,prnwidth,remainder
  REAL(iwp)           :: aa,bb,cc,kx,ky,kz,det,theta,dtim,real_time
  !REAL(iwp)           :: val0 = 100.0_iwp
  REAL(iwp)           :: tol,alpha,beta,up,big,q
  REAL(iwp)           :: rho,cp,val0,etype
  REAL(iwp),PARAMETER :: zero = 0.0_iwp,penalty=1.e20_iwp
  REAL(iwp),PARAMETER :: t0 = 0.0_iwp
!  CHARACTER(LEN=15)   :: element
  CHARACTER(LEN=50)   :: fname,job_name,label
  CHARACTER(LEN=50)   :: program_name='xx12'
  LOGICAL             :: converged = .false.
  LOGICAL             :: solid=.true.
  CHARACTER(LEN=80)   :: cbuffer
  CHARACTER(LEN=15), INTENT(IN) :: argv,element 
  
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------ 
  
  REAL(iwp),ALLOCATABLE :: loads_pp(:),u_pp(:),p_pp(:),points(:,:),kay(:,:)
  REAL(iwp),ALLOCATABLE :: coord(:,:),fun(:),jac(:,:),der(:,:),deriv(:,:)
  REAL(iwp),ALLOCATABLE :: weights(:),d_pp(:),kc(:,:),pm(:,:),funny(:,:)
  REAL(iwp),ALLOCATABLE :: p_g_co_pp(:,:,:),storka_pp(:,:,:)
  REAL(iwp),ALLOCATABLE :: storkb_pp(:,:,:),x_pp(:),xnew_pp(:),pmul_pp(:,:)
  REAL(iwp),ALLOCATABLE :: utemp_pp(:,:),diag_precon_pp(:),diag_precon_tmp(:,:)
  REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:),timest(:),g_coord(:,:)
  REAL(iwp),ALLOCATABLE :: disp_pp(:),eld_pp(:,:)
  REAL(iwp),ALLOCATABLE :: val(:,:),val_f(:),store_pp(:),r_pp(:)
  REAL(iwp),ALLOCATABLE :: kcx(:,:),kcy(:,:),kcz(:,:)
  REAL(iwp),ALLOCATABLE :: eld(:),col(:,:),row(:,:),storkc_pp(:,:,:)
  REAL(iwp),ALLOCATABLE :: prop(:,:),amp(:),tempres(:)
  INTEGER,ALLOCATABLE   :: rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp(:,:),no(:)
  INTEGER,ALLOCATABLE   :: no_pp(:),no_f_pp(:),no_pp_temp(:)
  INTEGER,ALLOCATABLE   :: sense(:),node(:)
  INTEGER,ALLOCATABLE   :: etype_pp(:),nf(:,:),oldlds(:)
  
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
  
  CALL read_xx12(job_name,numpe,dtim,element,fixed_freedoms,limit,            &
                 loaded_nodes,meshgen,nels,nip,nn,nod,npri,nr,nstep,          &
                 partitioner,theta,tol,np_types,val0,el_print,i_o)
  
  CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)
  
  ndof = nod*nodof
  ntot = ndof
  
  ALLOCATE(g_num_pp(nod,nels_pp))
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp))
  IF (nr>0) ALLOCATE(rest(nr,nodof+1))
  ALLOCATE(etype_pp(nels_pp))
  ALLOCATE(prop(nprops,np_types))
  
  g_num_pp       = 0
  g_coord_pp     = zero
  IF (nr>0) rest = 0
  etype_pp       = 0
  prop           = zero
  q              = zero
  
  timest(2) = elap_time()
  
! CALL read_g_num_pp2(job_name,iel_start,nn,npes,numpe,g_num_pp)
  CALL read_elements(job_name,iel_start,nn,npes,numpe,etype_pp,g_num_pp)
  timest(3) = elap_time()
  
  IF(meshgen == 2)THEN
    IF(numpe==1) PRINT *,"Calling abaqus2sg, meshgen = ",meshgen
    CALL abaqus2sg(element,g_num_pp)
  END IF
  timest(4) = elap_time()
  
  CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
  timest(5) = elap_time()
  
  IF (nr>0) CALL read_rest(job_name,numpe,rest)
  timest(6) = elap_time()
  
  IF(numpe==1) PRINT *, "np_types = ", np_types
  
  fname = job_name(1:INDEX(job_name, " ")-1) // ".mat"  
  CALL read_materialValue(prop,fname,numpe,npes)
  
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
            pmul_pp(ntot,nels_pp))
  ALLOCATE (kcx(ntot,ntot),kcy(ntot,ntot),kcz(ntot,ntot),                     &
            eld(ntot),col(ntot,1),row(1,ntot),storkc_pp(ntot,ntot,nels_pp))
  ALLOCATE (amp(nstep))
  
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
  
  !nres = 11488731 ! 11488731 25% model or 118564 5% model
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
  
  ALLOCATE(loads_pp(neq_pp),diag_precon_pp(neq_pp),u_pp(neq_pp),              &
           d_pp(neq_pp),p_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),r_pp(neq_pp))
  
  loads_pp  = zero ; diag_precon_pp = zero ; u_pp = zero ; r_pp    = zero
  d_pp      = zero ; p_pp           = zero ; x_pp = zero ; xnew_pp = zero
  
  timest(9) = elap_time()
  IF(numpe==1) PRINT *, "End of 7"
  
!------------------------------------------------------------------------------
! 8. Element stiffness integration and storage
!------------------------------------------------------------------------------
  
  CALL sample(element,points,weights)
  
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
      kc         = kc + MATMUL(MATMUL(TRANSPOSE(deriv),kay),deriv)*det*weights(i)
      pm         = pm + MATMUL(TRANSPOSE(funny),funny)*det*weights(i)*rho*cp

    END DO gauss_pts
    
    storka_pp(:,:,iel)=pm+kc*theta*dtim
    storkb_pp(:,:,iel)=pm-kc*(1._iwp-theta)*dtim
    
  END DO elements_3
  
  timest(10) = elap_time()
  IF(numpe==1) PRINT *, "End of 8"
  
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
  
  timest(11) = elap_time()
  IF(numpe==1) PRINT *, "End of 9" 

!------------------------------------------------------------------------------
! 10. Write case file
!------------------------------------------------------------------------------
  
  fname   = job_name(1:INDEX(job_name, " ")-1)//".bin.ensi.case"
  OPEN(12,FILE=fname)

  WRITE(12,'(A/A)')    "#", "# Post-processing file generated by subroutine &
                             &WRITE_ENSI in "
  WRITE(12,'(A,A,/A)') "#"," Smith, Griffiths and Margetts, 'Programming the &
                             &Finite Element Method',","# Wiley, 2013."        
  WRITE(12,'(A/A/A)')  "#","# Ensight Gold Format","#"
  WRITE(12,'(2A/A)')   "# Problem name: ",job_name(1:INDEX(job_name, " ")-1),"#"
  WRITE(12,'(A/A/A)')  "FORMAT","type:  ensight gold","GEOMETRY"
  WRITE(12,'(2A/A)')   "model: 1  ",job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.geo',"VARIABLE"
  WRITE(12,'(2A)')     "scalar per element:  material      ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.MATID'
  WRITE(12,'(2A)')     "scalar per element:  material_kx   ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_1_kx'
  WRITE(12,'(2A)')     "scalar per element:  material_ky   ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_2_ky'
  WRITE(12,'(2A)')     "scalar per element:  material_kz   ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_3_kz'
  WRITE(12,'(2A)')     "scalar per element:  material_rho  ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_4_rho'
  WRITE(12,'(2A)')     "scalar per element:  material_cp   ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_5_cp'
  WRITE(12,'(2A)')     "scalar per node: 1   temperature   ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.NDTTR-******'
!  IF(solid) THEN
!    WRITE(12,'(2A)')   "scalar per node:     restraint     ",                &
!                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.NDBND'
!    WRITE(12,'(2A)')   "vector per node:     displacement  ",                &
!                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.DISPL-******'
!  ELSE
!    WRITE(12,'(2A)')   "scalar per node:     pressure      ",                &
!                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.PRESSURE-******'
!  END IF
  WRITE(12,'(2A)')     "scalar per node:     fixed_nodes   ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.NDFIX'
  WRITE(12,'(2A)')     "scalar per node:     loaded_nodes  ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.NDLDS'
  WRITE(12,'(A/A)')     "TIME","time set:     1"
  WRITE(12,'(A,I5)')    "number of steps:",nstep/npri+1
  WRITE(12,'(A,I5)')    "filename start number:",1
  WRITE(12,'(A,I5)')    "filename increment:",1
  WRITE(12,'(A)')       "time values:"
  
  prnwidth = 5
  j=1
  WRITE(12,'(E12.5)',ADVANCE='no') zero
  timesteps: DO i=1,nstep
    real_time = i*dtim
    IF(i/npri*npri==i)THEN
      WRITE(12,'(E12.5)',ADVANCE='no') real_time
      j=j+1
    END IF
    IF(j==prnwidth)THEN
        WRITE(12,*)''
        j=0
    END IF
  END DO timesteps
  WRITE(12,*)''

!  remainder = mod(nstep/npri,prnwidth)
!  n         = ((nstep/npri) - remainder)/prnwidth
!  IF(nstep/npri<=prnwidth) THEN
!    DO i=0,nstep-1,npri
!      IF(i==nstep-1) THEN
!        WRITE(12,'(E12.5)') i*dtim
!      ELSE
!        WRITE(12,'(E12.5)',ADVANCE='no') i*dtim
!      END IF
!    END DO
!  ELSE
!    IF(remainder==0) THEN
!      DO j=0,n-1
!        m = j*prnwidth
!        l = j*prnwidth+prnwidth
!        WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
!      END DO
!    ELSE
!!     DO j=1,n-1
!      DO j=0,n-1
!        m = j*prnwidth
!        l = j*prnwidth+prnwidth
!        WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
!      END DO
!      m = n*prnwidth
!      l = n*prnwidth+remainder
!      DO i=m,l-1
!        IF(i==l-1) THEN
!          WRITE(12,'(E12.5)') dtim*i
!        ELSE
!          WRITE(12,'(E12.5)',ADVANCE='no') dtim*i
!        END IF
!      END DO
!    END IF
!  END IF
 
  CLOSE(12)
  
  IF(numpe==1) PRINT *, "End of 10" 
  
!----------------------------------------------------------------------------
! 11. Write geometry file
!----------------------------------------------------------------------------

  fname   = job_name(1:INDEX(job_name, " ")-1)//".bin.ensi.geo"
  OPEN(13,FILE=fname,STATUS="REPLACE",FORM="UNFORMATTED",                   &
                ACTION="WRITE", ACCESS="STREAM")

  cbuffer = "C Binary"                     ; WRITE(13) cbuffer
  cbuffer = "Problem name: "//job_name(1:INDEX(job_name, " ")-1) ; WRITE(13) cbuffer
  cbuffer = "Geometry files"               ; WRITE(13) cbuffer
  cbuffer = "node id off"                  ; WRITE(13) cbuffer
  cbuffer = "element id off"               ; WRITE(13) cbuffer
!  cbuffer = "node id given"                  ; WRITE(13) cbuffer
!  cbuffer = "element id given"               ; WRITE(13) cbuffer
  cbuffer = "part"                         ; WRITE(13) cbuffer
  WRITE(13) int(1,kind=c_int)
!  IF(ndim==2) THEN 
!     cbuffer = "2d-mesh"                   ; WRITE(13) cbuffer
!  END IF
  IF(ndim==3) THEN
     cbuffer = "Volume"                    ; WRITE(13) cbuffer
!     cbuffer = "Volume Mesh"                    ; WRITE(13) cbuffer
  END IF
  cbuffer = "coordinates"                  ; WRITE(13) cbuffer
  
  ALLOCATE(g_coord(ndim,nn))
  g_coord = 0.0_iwp
!  CALL READ_NODES(fname,nn,nn_start,numpe,g_coord)
  fname   = job_name(1:INDEX(job_name, " ")-1)//".d"
  CALL READ_NODES(fname,nn,1,numpe,g_coord)

!! Print to stdout to check values
!  IF(numpe==1) PRINT *, "g_coord"  
  WRITE(13) int(nn,kind=c_int)
  DO j=1,ndim
    DO i=1,nn
!      IF(numpe==1) PRINT *, g_coord(j,i)
      WRITE(13) real(g_coord(j,i),kind=c_float)
    END DO
  END DO
  
!  IF(ndim==2) THEN ! ensight requires zeros for the z-ordinate
!    DO i=1,nn
!      WRITE(13,'(A)') " 0.00000E+00" ! needs fixing for binary
!    END DO
!  END IF
  
!  SELECT CASE(element)
!    CASE('triangle')
!      SELECT CASE(nod)
!!         CASE(3)
!!           WRITE(13,'(A/I10)') "tria3", nels
!!           DO i = 1,nels
!!             WRITE(13,'(3I10)')g_num_pp(3,i),g_num_pp(2,i),g_num_pp(1,i)
!!           END DO
!        CASE DEFAULT
!          WRITE(13,'(A)')   "# Element type not recognised"
!      END SELECT
!    CASE('quadrilateral')
!      SELECT CASE(nod)
!        CASE(4)
!          WRITE(13,'(A/I10)') "quad4", nels
!          DO i = 1,nels
!            WRITE(13,'(4I10)')g_num_pp(1,i),g_num_pp(4,i),g_num_pp(3,i),g_num_pp(2,i)
!          END DO
!        CASE(8)
!          WRITE(13,'(A/I10)') "quad8", nels
!          DO i = 1,nels
!            WRITE(13,'(8I10)')g_num_pp(1,i),g_num_pp(7,i),g_num_pp(5,i),g_num_pp(3,i),   &
!                              g_num_pp(8,i),g_num_pp(6,i),g_num_pp(4,i),g_num_pp(2,i)
!          END DO
!        CASE DEFAULT
!          WRITE(13,'(A)')   "# Element type not recognised"
!      END SELECT
!    CASE('hexahedron')
!      SELECT CASE(nod)
!        CASE(8)
!          cbuffer = "hexa8"       ; WRITE(13) cbuffer
!          WRITE(13) int(nels,kind=c_int)
!          DO i = 1,nels
!            WRITE(13) int(g_num_pp(1,i),kind=c_int),int(g_num_pp(4,i),kind=c_int),&
!                      int(g_num_pp(8,i),kind=c_int),int(g_num_pp(5,i),kind=c_int),&
!                      int(g_num_pp(2,i),kind=c_int),int(g_num_pp(3,i),kind=c_int),&
!                      int(g_num_pp(7,i),kind=c_int),int(g_num_pp(6,i),kind=c_int)
!          END DO
!        CASE(20)
!          cbuffer = "hexa20"       ; WRITE(13) cbuffer
!          WRITE(13) int(nels,kind=c_int)
!          DO i = 1,nels
!            WRITE(13)                                                       &
!              int(g_num_pp(1,i),kind=c_int), int(g_num_pp(7,i),kind=c_int),       &
!              int(g_num_pp(19,i),kind=c_int),int(g_num_pp(13,i),kind=c_int),      &
!              int(g_num_pp(3,i),kind=c_int),int(g_num_pp(5,i),kind=c_int),        &
!              int(g_num_pp(17,i),kind=c_int),int(g_num_pp(15,i),kind=c_int),      &
!              int(g_num_pp(8,i),kind=c_int),int(g_num_pp(12,i),kind=c_int),       &
!              int(g_num_pp(20,i),kind=c_int),int(g_num_pp(9,i),kind=c_int),       &
!              int(g_num_pp(4,i),kind=c_int),int(g_num_pp(11,i),kind=c_int),       &
!              int(g_num_pp(16,i),kind=c_int),int(g_num_pp(10,i),kind=c_int),      &
!              int(g_num_pp(2,i),kind=c_int),int(g_num_pp(6,i),kind=c_int),        &
!              int(g_num_pp(18,i),kind=c_int),int(g_num_pp(14,i),kind=c_int) 
!          END DO
!        CASE DEFAULT
!          cbuffer = "# Element type not recognised" ; WRITE(13) cbuffer
!      END SELECT
!    CASE('tetrahedron')
!      SELECT CASE(nod)
!        CASE(4)
          cbuffer = "tetra4" ; WRITE(13) cbuffer
          WRITE(13) int(nels,kind=c_int)
          DO i = 1,nels
            WRITE(13) int(g_num_pp(1,i),kind=c_int),int(g_num_pp(3,i),kind=c_int), &
                      int(g_num_pp(2,i),kind=c_int),int(g_num_pp(4,i),kind=c_int)
!            WRITE(13) int(g_num_pp(1,i),kind=c_int),int(g_num_pp(2,i),kind=c_int), &
!                      int(g_num_pp(3,i),kind=c_int),int(g_num_pp(4,i),kind=c_int)
!!-----------Print to stdout to check values
!            IF(numpe==1) PRINT *, g_num_pp(1,i),g_num_pp(3,i),g_num_pp(2,i),       &
!                                  g_num_pp(4,i)
          END DO
!        CASE DEFAULT
!          cbuffer = "# Element type not recognised" ; WRITE(13)
!      END SELECT
!    CASE DEFAULT
!      cbuffer = "# Element type not recognised" ; WRITE(13)
!  END SELECT
  
  CLOSE(13)
  
  IF(numpe==1) PRINT *, "End of 11"
  
!-----------------------------------------------------------------------------
! 12. Write file containing material IDs
!-----------------------------------------------------------------------------
  
  OPEN(14,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.MATID',         &
          STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")
  
  cbuffer = "Alya Ensight Gold --- Scalar per-element variable file"
  WRITE(14) cbuffer
  cbuffer = "part"  ; WRITE(14) cbuffer
  WRITE(14) int(1,kind=c_int)
  
!  SELECT CASE(element)
!    CASE('triangle')
!      SELECT CASE(nod) 
!        CASE(3)
!          cbuffer = "tria3" ; WRITE(14) cbuffer
!        CASE DEFAULT
!          WRITE(14,'(A)') "# Element type not recognised"
!      END SELECT
!    CASE('quadrilateral')
!      SELECT CASE(nod) 
!        CASE(4)
!          cbuffer = "quad4" ; WRITE(14) cbuffer
!        CASE(8)
!          cbuffer = "quad8" ; WRITE(14) cbuffer
!        CASE DEFAULT
!          WRITE(14,'(A)') "# Element type not recognised"
!      END SELECT
!    CASE('hexahedron')
!      SELECT CASE(nod) 
!        CASE(8)
!          cbuffer = "hexa8"   ; WRITE(14) cbuffer
!        CASE(20)
!          cbuffer = "hexa20"  ; WRITE(14) cbuffer
!        CASE DEFAULT
!          WRITE(14,'(A)') "# Element type not recognised"
!      END SELECT
!    CASE('tetrahedron')
!      SELECT CASE(nod)
!        CASE(4)
          cbuffer = "tetra4"  ; WRITE(14) cbuffer
!        CASE DEFAULT
!        WRITE(14,'(A)') "# Element type not recognised"
!      END SELECT
!    CASE DEFAULT
!      WRITE(14,'(A)')   "# Element type not recognised"
!  END SELECT
  DO i = 1,nels
    etype=etype_pp(i)*1.0
    WRITE(14) real(etype,kind=c_float)
!    WRITE(14) int(etype_pp(i),kind=c_int) 
!    IF(numpe==1) PRINT *, etype_pp(i)
  END DO
  
  CLOSE(14)
  
  IF(numpe==1) PRINT *, "End of 12"
  
!-----------------------------------------------------------------------------
! 13. Write files containing material properties
!-----------------------------------------------------------------------------
  
  OPEN(15,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_1_kx',    &
          STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")
  
  cbuffer = "Alya Ensight Gold --- Scalar per-element variable file"
  WRITE(15) cbuffer
  cbuffer = "part"  ; WRITE(15) cbuffer
  WRITE(15) int(1,kind=c_int)
  cbuffer = "tetra4"  ; WRITE(15) cbuffer
  
  DO iel = 1,nels
    WRITE(15) real(prop(1,etype_pp(iel)),kind=c_float)
  END DO
  
  CLOSE(15)

  OPEN(16,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_2_ky',    &
          STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")
  
  cbuffer = "Alya Ensight Gold --- Scalar per-element variable file"
  WRITE(16) cbuffer
  cbuffer = "part"  ; WRITE(16) cbuffer
  WRITE(16) int(1,kind=c_int)
  cbuffer = "tetra4"  ; WRITE(16) cbuffer
  
  DO iel = 1,nels
    WRITE(16) real(prop(2,etype_pp(iel)),kind=c_float)
  END DO
  
  CLOSE(16)

  OPEN(17,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_3_kz',    &
          STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")
  
  cbuffer = "Alya Ensight Gold --- Scalar per-element variable file"
  WRITE(17) cbuffer
  cbuffer = "part"  ; WRITE(17) cbuffer
  WRITE(17) int(1,kind=c_int)
  cbuffer = "tetra4"  ; WRITE(17) cbuffer
  
  DO iel = 1,nels
    WRITE(17) real(prop(3,etype_pp(iel)),kind=c_float)
  END DO
  
  CLOSE(17)

  OPEN(18,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_4_rho',    &
          STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")
  
  cbuffer = "Alya Ensight Gold --- Scalar per-element variable file"
  WRITE(18) cbuffer
  cbuffer = "part"  ; WRITE(18) cbuffer
  WRITE(18) int(1,kind=c_int)
  cbuffer = "tetra4"  ; WRITE(18) cbuffer
  
  DO iel = 1,nels
    WRITE(18) real(prop(4,etype_pp(iel)),kind=c_float)
  END DO
  
  CLOSE(18)

  OPEN(19,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_5_cp',    &
          STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")
  
  cbuffer = "Alya Ensight Gold --- Scalar per-element variable file"
  WRITE(19) cbuffer
  cbuffer = "part"  ; WRITE(19) cbuffer
  WRITE(19) int(1,kind=c_int)
  cbuffer = "tetra4"  ; WRITE(19) cbuffer
  
  DO iel = 1,nels
    WRITE(19) real(prop(5,etype_pp(iel)),kind=c_float)
  END DO
  
  CLOSE(19)
  
  IF(numpe==1) PRINT *, "End of 13"
  
!-----------------------------------------------------------------------------
! 14. Write file containing fixed nodes
!-----------------------------------------------------------------------------  
  
  IF(fixed_freedoms > 0) THEN
    ALLOCATE(node(fixed_freedoms),sense(fixed_freedoms),val_f(fixed_freedoms))
    IF(numpe==1)THEN
      fname = job_name(1:INDEX(job_name, " ")-1) // ".fix"
      OPEN(20,FILE=fname,STATUS='OLD',ACTION='READ')
      DO i = 1,fixed_freedoms
        READ(20,*)node(i),sense(i),val_f(i)
      END DO
      CLOSE(20)
    END IF
    
    OPEN(21,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.NDFIX',         &
            STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")
    
    cbuffer = "Alya Ensight Gold --- Scalar per-partial-node variable file"
    WRITE(21) cbuffer
    cbuffer = "part"  ; WRITE(21) cbuffer
    WRITE(21) int(1,kind=c_int)
    cbuffer = "coordinates"  ; WRITE(21) cbuffer
    
    j=1
    DO i=1,nn
      IF(i==node(j)) THEN
        WRITE(21) real(val_f(j),kind=c_float)
        j=j+1
      ELSE
        WRITE(21) real(zero,kind=c_float)
      END IF
    END DO
    
    DEALLOCATE(node)
    CLOSE(21)
  END IF
  IF(numpe==1) PRINT *, "End of 14"
  
!-----------------------------------------------------------------------------
! 15. Write file containing loaded nodes
!-----------------------------------------------------------------------------
  
  loaded_freedoms = loaded_nodes ! hack
  IF(loaded_freedoms > 0) THEN
    ALLOCATE(node(loaded_freedoms),val(nodof,loaded_freedoms))
    IF(numpe==1)THEN
      fname = job_name(1:INDEX(job_name, " ")-1) // ".lds"
      OPEN(22, FILE=fname, STATUS='OLD', ACTION='READ')
      DO i = 1,loaded_nodes 
        READ(22,*) node(i),val(:,i)
      END DO
      CLOSE(22)
    END IF
  
    OPEN(23,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.NDLDS',         &
            STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")
    
    cbuffer = "Alya Ensight Gold --- Scalar per-partial-node variable file"
    WRITE(23) cbuffer
    cbuffer = "part"  ; WRITE(23) cbuffer
    WRITE(23) int(1,kind=c_int)
    cbuffer = "coordinates"  ; WRITE(23) cbuffer
    
    j=1
    DO i=1,nn
      IF(i==node(j)) THEN
        WRITE(23) real(val(1,j),kind=c_float)
!        IF(numpe==1) PRINT *,val(1,j)
        j=j+1
      ELSE
        WRITE(23) real(zero,kind=c_float)
!        IF(numpe==1) PRINT *,0.0
      END IF
    END DO
    
    DEALLOCATE(node)
    CLOSE(23)
  END IF
    
  IF(numpe==1) PRINT *, "End of 15"
  
  CALL shutdown()
  
END PROGRAM pf2ensibin
