PROGRAM p1210      
!------------------------------------------------------------------------------
!  Program 12.10 forced vibration of an elastic-plastic(Von Mises) solid
!  Viscoplastic strain method, lumped mass, explicit integration
!
!  Original program from Smith I.M. and Griffiths D.V. "Programming the Finite 
!  Element Method", Edition 4, Wiley, 2004.
!
!  Modified by Lee Margetts 29-03-2010
!------------------------------------------------------------------------------

  USE precision     ; USE global_variables ; USE mp_interface
  USE input         ; USE output           ; USE loading
  USE timing        ; USE maths            ; USE gather_scatter
  USE partition     ; USE elements         ; USE steering
  USE pcg           ; USE plasticity 
 
  IMPLICIT NONE

!------------------------------------------------------------------------------ 
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  ! neq,ntot are now global variables - not declared 

  INTEGER,PARAMETER     :: nodof=3,ndim=3,nst=6
  INTEGER               :: nn,nr,nip,loaded_nodes,nres
  INTEGER               :: i,j,k,jj,iel,nstep,npri,num_no
  INTEGER               :: no_index_start,is,it,nlen,ndof,nels,npes_pp
  INTEGER               :: node_end,node_start,nodes_pp
  INTEGER               :: argc,iargc,meshgen
  REAL(iwp)             :: rho,dtim,e,v,det,sbary,pload,sigm,f,fnew,fac
  REAL(iwp)             :: volume,sbar,dsbar,lode_theta,real_time,tload
  REAL(iwp),PARAMETER   :: zero = 0.0_iwp    
  CHARACTER(LEN=15)     :: element
  CHARACTER(LEN=50)     :: program_name='p1210' 
  CHARACTER(LEN=50)     :: fname,job_name,label
  CHARACTER(LEN=6)      :: step
 
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE :: points(:,:),bdylds_pp(:),x1_pp(:),d1x1_pp(:)
  REAL(iwp),ALLOCATABLE :: stressv(:),pl(:,:),emm(:),d2x1_pp(:)
  REAL(iwp),ALLOCATABLE :: tensor_pp(:,:,:),etensor_pp(:,:,:),val(:),dee(:,:)
  REAL(iwp),ALLOCATABLE :: mm_pp(:),jac(:,:),weights(:),der(:,:),deriv(:,:)
  REAL(iwp),ALLOCATABLE :: bee(:,:),eld(:),eps(:),sigma(:),bload(:),eload(:)
  REAL(iwp),ALLOCATABLE :: mm_tmp(:,:),p_g_co_pp(:,:,:),pmul_pp(:,:)
  REAL(iwp),ALLOCATABLE :: utemp_pp(:,:),disp_pp(:),fext_pp(:),timest(:)
  INTEGER,ALLOCATABLE   :: rest(:,:),no(:),no_local(:),g_num_pp(:,:)
  INTEGER,ALLOCATABLE   :: g_g_pp(:,:),no_local_temp(:),node(:)
   
!------------------------------------------------------------------------------
! 3. Read input data and initialise problem. The master processor imports the 
!    mesh, material properties and restraints. Each processor stores the 
!    information that is needed locally.
!------------------------------------------------------------------------------

  ALLOCATE(timest(20))
  timest    = zero 
  timest(1) = elap_time()
  
  CALL find_pe_procs(numpe,npes)
  argc = iargc()
  IF (argc /= 1) CALL job_name_error(numpe,program_name)
  CALL GETARG(1, job_name) 

  CALL read_p1210(job_name,numpe,dtim,e,element,loaded_nodes,meshgen,nels,    &
                  nip,nn,nod,npri,nr,nstep,pload,rho,sbary,v)
 
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

  ALLOCATE(points(nip,ndim),weights(nip),dee(nst,nst),pmul_pp(ntot,nels_pp),  &
           tensor_pp(nst,nip,nels_pp),no(loaded_nodes),pl(nst,nst),           &
           etensor_pp(nst,nip,nels_pp),jac(ndim,ndim),der(ndim,nod),          &
           deriv(ndim,nod),bee(nst,ntot),eld(ntot),eps(nst),                  &
           sigma(nst),emm(ntot),bload(ntot),eload(ntot),stressv(nst),         &
           g_g_pp(ntot,nels_pp),mm_tmp(ntot,nels_pp),utemp_pp(ntot,nels_pp))

 timest(2)  = elap_time()

!-------------------------------------------------------------------------------
! 5. Create node freedom array, find element steering and neq
!-------------------------------------------------------------------------------

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

 ALLOCATE(bdylds_pp(neq_pp),x1_pp(neq_pp),d1x1_pp(neq_pp),d2x1_pp(neq_pp),    &
          mm_pp(neq_pp),fext_pp(neq_pp))
          
 bdylds_pp = zero  ;  x1_pp     = zero  ;   d1x1_pp   = zero
 d2x1_pp   = zero  ;  mm_pp     = zero  ;   fext_pp   = zero
 
 timest(5) = elap_time()

!------------------------------------------------------------------------------
! 8. Calculate diagonal mass matrix 
!------------------------------------------------------------------------------
 
 mm_tmp     = zero
 
 CALL sample(element,points,weights)
 
 elements_1: DO iel=1,nels_pp
   
   volume = zero
   
   gauss_pts_1: DO i=1,nip
     CALL shape_der (der,points,i)
     jac    = MATMUL(der,g_coord_pp(:,:,iel))
     det    = determinant(jac)
     volume = volume+det*weights(i)*rho
   END DO gauss_pts_1
   
   emm           = volume/13._iwp
   emm(1:19:6)   = emm(4)*.125_iwp
   emm(2:20:6)   = emm(4)*.125_iwp
   emm(3:21:6)   = emm(4)*.125_iwp
   emm(37:55:6)  = emm(4)*.125_iwp
   emm(38:56:6)  = emm(4)*.125_iwp
   emm(39:57:6)  = emm(4)*.125_iwp
   mm_tmp(:,iel) = mm_tmp(:,iel)+emm 

 END DO elements_1

 CALL scatter(mm_pp,mm_tmp)
 DEALLOCATE(mm_tmp)

 timest(6) = elap_time()
 
!-------------------------------------------------------------------------------
! 9. Read in applied forces and assign to equations
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
  
  timest(7) = elap_time()

!------------------------------------------------------------------------------
! 10. Allocate disp_pp array and open results file
!------------------------------------------------------------------------------

  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
  ALLOCATE(disp_pp(nodes_pp*ndim))
    
  disp_pp = zero
  
  IF(numpe==1) THEN
    fname   = job_name(1:INDEX(job_name, " ")-1)//".dis"
    OPEN(12, file=fname, status='replace', action='write')
  END IF

!------------------------------------------------------------------------------
! 11. Explicit integration loop
!------------------------------------------------------------------------------

  tensor_pp  = zero
  etensor_pp = zero
  real_time  = zero

  time_steps: DO jj=1,nstep

    real_time = real_time+dtim
    x1_pp     = x1_pp+(d1x1_pp+d2x1_pp*dtim*.5_iwp)*dtim
    bdylds_pp = zero

!------------------------------------------------------------------------------
! 12. Element stress-strain relationship
!------------------------------------------------------------------------------
   
    timest(8) = elap_time()
    pmul_pp  = zero
    utemp_pp = zero
   
    CALL gather(x1_pp,pmul_pp)

    elements_2: DO iel=1,nels_pp          

      bload = zero
      eld   = pmul_pp(:,iel)

      gauss_pts_2: DO i=1,nip
        dee     = zero
        CALL deemat(e,v,dee)
        CALL shape_der(der,points,i)
        jac     = MATMUL(der,p_g_co_pp(:,:,iel))
        det     = determinant(jac)
        CALL invert(jac)
        deriv   = MATMUL(jac,der)
        CALL beemat(deriv,bee)
        eps     = MATMUL(bee,pmul_pp(:,iel))
        eps     = eps-etensor_pp(:,i,iel)
        sigma   = MATMUL(dee,eps)
        stressv = sigma+tensor_pp(:,i,iel)
        CALL invar(stressv,sigm,dsbar,lode_theta)
        fnew    = dsbar-sbary

!------------------------------------------------------------------------------
! 13. Check whether yield is violated
!------------------------------------------------------------------------------

        IF(fnew>=.0_iwp)THEN
          stressv = tensor_pp(:,i,iel)
          CALL invar(stressv,sigm,sbar,lode_theta)
          f       = sbar-sbary
          fac     = fnew/(fnew-f)
          stressv = tensor_pp(:,i,iel)+(1._iwp-fac)*sigma
          CALL vmpl(e,v,stressv,pl)
          dee     = dee-fac*pl
        END IF
       
        sigma = MATMUL(dee,eps)
        sigma = sigma+tensor_pp(:,i,iel)
        eload = MATMUL(sigma,bee)
        bload = bload+eload*det*weights(i)
 
!------------------------------------------------------------------------------
! 14. Update the gauss points
!------------------------------------------------------------------------------

        tensor_pp(:,i,iel)  = sigma 
        etensor_pp(:,i,iel) = etensor_pp(:,i,iel)+eps
 
      END DO gauss_pts_2      

      utemp_pp(:,iel) = utemp_pp(:,iel)-bload    

    END DO elements_2

    CALL scatter(bdylds_pp,utemp_pp)

    bdylds_pp = bdylds_pp+fext_pp*2._iwp
    bdylds_pp = bdylds_pp/mm_pp
    d1x1_pp   = d1x1_pp+(d2x1_pp+bdylds_pp)*.5_iwp*dtim
    d2x1_pp   = bdylds_pp

    timest(9) = timest(9) + (elap_time() - timest(8))
    timest(8) = elap_time()
    
!------------------------------------------------------------------------------
! 15. Output displacements
!------------------------------------------------------------------------------

    IF(jj==jj/npri*npri) THEN 

      disp_pp   = zero 
      utemp_pp  = zero
      label     = "*DISPLACEMENT"  
      CALL gather(x1_pp(1:),utemp_pp)
      CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,            &
                         node_start,node_end,utemp_pp,disp_pp,1)
      CALL write_nodal_variable(label,12,jj,nodes_pp,npes,numpe,ndim,disp_pp) 
                           
    END IF

    timest(10) = timest(10) + (elap_time() - timest(8))
     
  END DO time_steps

  timest(11) = elap_time()
 
!------------------------------------------------------------------------------
! 16. Output debugging information and performance data
!------------------------------------------------------------------------------

  CALL WRITE_P1210(job_name,neq,nn,npes,nr,numpe,timest,tload)
  
  IF(numpe==1) CLOSE(12)      

  CALL shutdown() 

END PROGRAM p1210
