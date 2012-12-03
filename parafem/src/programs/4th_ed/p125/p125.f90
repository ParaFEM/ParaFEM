PROGRAM p125      
!-------------------------------------------------------------------------
!   Program 8.4 conduction equation on a 3-d box volume using 8-node
!   hexahedral elements and a simple explicit algorithm : parallel version
!   box_bc : write on processor it at freedom nres
!-------------------------------------------------------------------------
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

 INTEGER,PARAMETER    :: nodof=1,ndim=3
 INTEGER              :: nels,ndof,npes_pp,nn,nr,nip,nod,i,j,k,iel
 INTEGER              :: neq_temp,nn_temp,nstep,npri,nres,it,is,nlen
 INTEGER              :: argc,iargc,meshgen,partitioner
 INTEGER              :: loaded_nodes,fixed_freedoms
 REAL(iwp)            :: kx,ky,kz,det,dtim,val0,real_time
 REAL(iwp),PARAMETER  :: zero=0.0_iwp
 CHARACTER(LEN=15)    :: element
 CHARACTER(LEN=50)    :: fname,job_name,label
 CHARACTER(LEN=50)    :: program_name='p125'

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------ 

 REAL(iwp),ALLOCATABLE :: loads_pp(:),points(:,:),kay(:,:)
 REAL(iwp),ALLOCATABLE :: jac(:,:),der(:,:),deriv(:,:),weights(:),kc(:,:)
 REAL(iwp),ALLOCATABLE :: pm(:,:),funny(:,:),globma_pp(:)
 REAL(iwp),ALLOCATABLE :: fun(:),store_pm_pp(:,:,:),newlo_pp(:),mass(:)
 REAL(iwp),ALLOCATABLE :: globma_tmp(:,:),pmul_pp(:,:),utemp_pp(:,:)
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
 
  argc = iargc()
  IF(argc /= 1) CALL job_name_error(numpe,program_name)
  CALL GETARG(1,job_name)

  CALL read_p125(job_name,numpe,dtim,element,fixed_freedoms,kx,ky,kz,         &
                 loaded_nodes,meshgen,nels,nip,nn,nod,npri,nr,nres,nstep,     &
                 partitioner,val0)

! READ(10,*)nels,nxe,nze,nip,aa,bb,cc,kx,ky,kz,dtim,nstep,npri,val0

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

!------------------------------------------------------------------------------
! 4. Allocate dynamic arrays used in main program
!------------------------------------------------------------------------------

  ALLOCATE (points(nip,ndim),weights(nip),kay(ndim,ndim),                     &
            jac(ndim,ndim),                                                   &
            der(ndim,nod),deriv(ndim,nod),kc(ntot,ntot),                      &
            g(ntot),funny(1,nod),num(nod),g_g_pp(ntot,nels_pp),               &
            store_pm_pp(ntot,ntot,nels_pp),mass(ntot),fun(nod),pm(ntot,ntot), &
            globma_tmp(ntot,nels_pp),pmul_pp(ntot,nels_pp),                   &
            utemp_pp(ntot,nels_pp))

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

  DO i = 1,neq_pp
    IF(nres==ieq_start+i-1) THEN
      it = numpe; is = i
    END IF
  END DO

  timest(8) = elap_time()

  IF(numpe==it)THEN
    fname=job_name(1:INDEX(job_name," ")-1) // ".res"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')              
    WRITE(11,'(A,I5,A)')"This job ran on ",npes,"  processors" 
    WRITE(11,'(A)')"Global coordinates and node numbers" 
    DO i=1,nels_pp,nels_pp-1
      WRITE(11,'(A,I8)')"Element ",i
      num=g_num_pp(:,i)
      DO k=1,nod
        WRITE(11,'(A,I8,3E12.4)')"  Node",num(k),g_coord_pp(k,:,i)
      END DO
    END DO
    WRITE(11,'(A,3(I8,A))')"There are ",nn," nodes",nr," restrained and",  &
      neq," equations"
    WRITE(11,*)"Time after setup  is  :",elap_time()-timest(1)
  END IF

!------------------------------------------------------------------------------
! 7. Allocate arrays dimensioned by neq_pp 
!------------------------------------------------------------------------------

  ALLOCATE(loads_pp(neq_pp),newlo_pp(neq_pp),globma_pp(neq_pp))

  loads_pp = zero; newlo_pp = zero ; globma_pp = zero

!------------------------------------------------------------------------------
! 8. Element stiffness integration and invert mass
!------------------------------------------------------------------------------
  
  CALL sample(element,points,weights)

  globma_tmp = .0_iwp
  kay        = 0.0_iwp
  kay(1,1)   = kx
  kay(2,2)   = ky
  kay(3,3)   = kz

  elements_3: DO iel=1,nels_pp
    kc=0.0_iwp
    pm=0.0_iwp
    gauss_pts: DO i=1,nip
      CALL shape_der(der,points,i)
      CALL shape_fun(fun,points,i)
      funny(1,:)=fun(:)
      jac=MATMUL(der,g_coord_pp(:,:,iel))
      det=determinant(jac)
      CALL invert(jac)
      deriv=MATMUL(jac,der)
      kc=kc+MATMUL(MATMUL(TRANSPOSE(deriv),kay),deriv)*det*weights(i)
      pm=pm+MATMUL(TRANSPOSE(funny),funny)*det*weights(i) 
    END DO gauss_pts
    DO i=1,ntot
      mass(i)=sum(pm(i,:))
    END DO
    pm=.0_iwp
    DO i=1,ntot
      pm(i,i)=mass(i)
    END DO
    store_pm_pp(:,:,iel)=pm-kc*dtim
    DO i=1,ntot
      globma_tmp(i,iel)=globma_tmp(i,iel)+mass(i)
    END DO
  END DO elements_3
 
 IF(numpe==it)                                                            &
   WRITE(11,*)"Time after element integration is :",elap_time()-timest(1)
 
 CALL scatter(globma_pp,globma_tmp)
 globma_pp = 1._iwp/globma_pp

 loads_pp  = val0
 
 DEALLOCATE(globma_tmp)

!------------------------------------------------------------------------------
! 9. Time stepping recursion
!------------------------------------------------------------------------------

 IF(numpe==it)THEN
   WRITE(11,'(A)')"    Time     Pressure"
 END IF
 
 timesteps: DO j=1,nstep
   
   real_time = j*dtim              

!---------------  go round the elements  ---------------------------------

   utemp_pp = .0_iwp
   pmul_pp  = .0_iwp
   CALL gather(loads_pp,pmul_pp)
   elements_4: DO iel=1,nels_pp  
     pm = store_pm_pp(:,:,iel) 
     utemp_pp(:,iel) = utemp_pp(:,iel)+MATMUL(pm,pmul_pp(:,iel)) 
   END DO elements_4
   CALL scatter(newlo_pp,utemp_pp)

   loads_pp = newlo_pp*globma_pp
   newlo_pp = .0_iwp 

   IF(j/npri*npri==j.AND.numpe==it)                                       &
     WRITE(11,'(2E12.4)')real_time,loads_pp(is) 

 END DO timesteps
 
 IF(numpe==it)WRITE(11,*)"This analysis took  :",elap_time()-timest(1)
 
 CALL shutdown() 

END PROGRAM p125
