program XX7
!------------------------------------------------------------------------------
!      program xx7:  elastic analysis, large deformations
!------------------------------------------------------------------------------

  USE precision
  USE global_variables
  USE mp_interface
  USE input
  USE output
  USE gather_scatter
  USE partition
  USE maths
  USE timing
  USE large_strain
  

  IMPLICIT NONE

  ! neq, ntot declared in global_variables
  INTEGER :: nels, nn, nr, nip, nodof=3, nod, nst=6, loaded_nodes, nn_pp
  INTEGER :: nf_start, fmt=1, i, j, k, ndim=3, iters, limit, iel, nn_start
  INTEGER :: num_load_steps, iload, igauss, dimH, inewton, jump, npes_pp
  INTEGER :: partitioner=1
  REAL(iwp) :: e, v, det, tol, maxdiff, tol2, detF
  REAL(iwp) :: energy, energy1, rn0
  CHARACTER(len=15) :: element
  CHARACTER(len=30) :: text

  CHARACTER(len=50) :: fname_base, fname
  INTEGER :: argc, iargc
  INTEGER :: fixed_nodes, numfix_pp, fixdim, writetimes=0
  REAL(iwp) :: timest(20)
 
  LOGICAL :: converged, timewrite=.TRUE.

  !-------------------------- dynamic arrays--------------------------------
  REAL(iwp),ALLOCATABLE:: points(:,:),coord(:,:),weights(:),diag_precon_pp(:),&
                          r_pp(:),xnew_pp(:),bee(:,:),load_value(:,:),        &
                          diag_precon_tmp(:,:), storekm_pp(:,:,:),            &
                          g_coord_pp(:,:),disp(:,:),g_coord(:,:),val_pp(:),   &
                          disp_pp(:,:), res_pp(:)
  REAL(iwp),ALLOCATABLE:: fext_pp(:), fextpiece_pp(:), deltax_pp(:),          &
                          fint_pp(:), kmat_elem(:,:), kgeo_elem(:,:),         &
                          xnewel_pp(:,:), jacF(:,:), auxm(:,:),               &
                          derivFtran(:,:), derivF(:,:), jacFinv(:,:),         &
                          beeF(:,:), rightCG(:,:), defE(:,:), piolaS(:,:),    &
                          cmat(:,:,:,:), sigma(:,:), cspa(:,:,:,:),           &
                          sigma1C(:), storefint_pp(:,:), deeF(:,:),           &
                          geomH(:,:), fixed_value(:), fixval_pp(:),           &
                          fixvalpiece_pp(:), elemdisp(:)
  INTEGER,ALLOCATABLE  :: num(:),g_num(:,:),g_num_pp(:,:),g_g_pp(:,:),        &
                          load_node(:),rest(:,:),nf_pp(:,:),no_pp(:)
  INTEGER,ALLOCATABLE  :: comp(:,:), fixed_node(:), fixed_dof(:),             &
                          fixelem_pp(:), fixdof_pp(:)
 
  INTEGER :: nodes_pp, node_start, node_end, idx1, idx2
  REAL(iwp),ALLOCATABLE :: value_shape(:),xnewnodes_pp(:),shape_integral_pp(:,:)
  REAL(iwp),ALLOCATABLE :: stress_integral_pp(:,:), stressnodes_pp(:)
  REAL(iwp),ALLOCATABLE :: principal_integral_pp(:,:), princinodes_pp(:)
  REAL(iwp),ALLOCATABLE :: principal(:), reacnodes_pp(:)

  !-------------------------------------------------------------------------
  ! 1. Input and initialisation
  !-------------------------------------------------------------------------

  timest(1) = ELAP_TIME( )

  CALL FIND_PE_PROCS(numpe,npes)

  argc = IARGC()
!     Output:  argc i: Number of arguments in the command line,
!                      ./par131 arg1 arg2 arg3 ... argn

!  If I have forgotten to write the name of the file in the command
!  line, the program is stopped
  IF (argc /= 1) THEN
    IF (numpe==npes) THEN
      WRITE(*,*) "Need name of filename_base!!"
    END IF
    CALL SHUTDOWN
    STOP
  END IF

!     Input:  1: The first argument in the command line (arg1) 
  CALL GETARG(1, fname_base)
!     Output: fname_base ch: Name of the file (maybe par131) 
 
  fname = fname_base(1:INDEX(fname_base," ")-1) // ".dat"
  CALL READ_DATA_XX7(fname,numpe,nels,nn,nr,loaded_nodes,fixed_nodes,nip,     &
                     limit,tol,e,v,nod,num_load_steps,jump,tol2)

  IF (nels < npes) THEN
    CALL SHUTDOWN
    WRITE(*,*)"Error: less elements than processors"
    STOP
  END IF

  IF(numpe==npes) THEN
    fname = fname_base(1:INDEX(fname_base," ")-1) // ".res"
    OPEN (11, file=fname, status='replace', action='write')
  END IF

  IF (nod==8) THEN
    element='hexahedron'
    dimH=8
  ELSE IF (nod==4) THEN
    element='tetrahedron'
    dimH=4
  END IF

!------------------------------------------------------------------------------
! 1. Get integration Gauss points and weights in the element
!------------------------------------------------------------------------------
  
  ALLOCATE(points(ndim,nip), weights(nip))
  CALL GET_GAUSS_POINTS(element,points,weights)

!------------------------------------------------------------------------------
! 2. Import and distribute mesh
!------------------------------------------------------------------------------

  CALL CALC_NELS_PP(fname_base,nels,npes,numpe,partitioner,nels_pp)

  ALLOCATE(g_num_pp(nod, nels_pp)) 

  fname = fname_base(1:INDEX(fname_base, " ")-1) // ".d"
! CALL READ_ELEMENTS(fname,iel_start,nels,nn,numpe,g_num_pp)
  CALL READ_ELEMENTS_2(fname,npes,nn,numpe,g_num_pp)

  CALL CALC_NN_PP(g_num_pp,nn_pp,nn_start)

  ALLOCATE(g_coord_pp(ndim, nn_pp))
  CALL READ_NODES(fname,nn,nn_start,numpe,g_coord_pp)

!-----------------------------------------------------------------------------

  ntot = nod * nodof

  CALL CALC_NODES_PP(nn,npes,numpe,node_end,node_start,nodes_pp)


  ALLOCATE(coord(nod,ndim))
  ALLOCATE(bee(nst,ntot))
  ALLOCATE(num(nod))
  ALLOCATE(g_g_pp(ntot,nels_pp))
  ALLOCATE(rest(nr,nodof+1))
  ALLOCATE(load_value(ndim,loaded_nodes))
  ALLOCATE(load_node(loaded_nodes))
  ALLOCATE(nf_pp(nodof,nn_pp))
  ALLOCATE(storekm_pp(ntot,ntot,nels_pp))
  ALLOCATE(kmat_elem(ntot,ntot))
  ALLOCATE(kgeo_elem(ntot,ntot))
  ALLOCATE(xnewel_pp(ntot,nels_pp))
  ALLOCATE(comp(nod,ndim))
  ALLOCATE(jacF(ndim,ndim))
  ALLOCATE(auxm(nod,ndim))
  ALLOCATE(jacFinv(ndim,ndim))
  ALLOCATE(derivFtran(nod,ndim))
  ALLOCATE(derivF(ndim,nod))
  ALLOCATE(beeF(nst,ntot))
  ALLOCATE(rightCG(ndim,ndim))
  ALLOCATE(defE(ndim,ndim))
  ALLOCATE(piolaS(ndim,ndim))
  ALLOCATE(cmat(ndim,ndim,ndim,ndim))
  ALLOCATE(sigma(ndim,ndim))
  ALLOCATE(cspa(ndim,ndim,ndim,ndim))
  ALLOCATE(sigma1C(nst))
  ALLOCATE(storefint_pp(ntot,nels_pp))
  ALLOCATE(deeF(nst,nst))
  ALLOCATE(geomH(dimH,dimH))
  ALLOCATE(principal(ndim))
  ALLOCATE(elemdisp(ntot))
  ALLOCATE(value_shape(nod))

  !-------------------------------------------------------------------------
  ! 3. Read and organise the structure of degrees of freedom
  !-------------------------------------------------------------------------

  fname = fname_base(1:INDEX(fname_base, " ")-1) // ".bnd"
  CALL READ_RESTRAINTS(fname,numpe,rest)

  CALL REST_NF(nn_start,rest,nf_pp)

  DO iel = 1,nels_pp
    CALL NUM_TO_G2(g_num_pp(:,iel), nf_pp, g_g_pp(:,iel), nn_start)
  END DO

  CALL CALC_NEQ(nn,rest,neq)

  CALL COMPUTE_NPES_PP(nels,neq,nn,npes,numpe,g_num_pp,rest,npes_pp)

  CALL CALC_NEQ_PP(nels)

  CALL MAKE_GGL(npes_pp,npes,g_g_pp)

  timest(2) = ELAP_TIME()

!------------------------------------------------------------------------------
! 4. Read and distribute essential boundary conditions
!------------------------------------------------------------------------------

  numfix_pp = 0

  IF (fixed_nodes>0) THEN

    ALLOCATE(fixed_node(fixed_nodes))
    ALLOCATE(fixed_dof(fixed_nodes))
    ALLOCATE(fixed_value(fixed_nodes))

    CALL READ_FIXED(fname_base,numpe,fixed_node,fixed_dof,fixed_value) ! new order

    IF (element=='hexahedron') THEN
      fixdim=4
    ELSE IF (element=='tetrahedron') THEN
      fixdim=20
    END IF

    ALLOCATE(fixelem_pp(fixdim*fixed_nodes))
    ALLOCATE(fixdof_pp(fixdim*fixed_nodes))
    ALLOCATE(fixval_pp(fixdim*fixed_nodes))
    ALLOCATE(fixvalpiece_pp(fixdim*fixed_nodes))
    
	fixelem_pp = 0
	fixdof_pp  = 0
	fixval_pp  = 0.0_iwp
    DO i = 1,fixed_nodes
      DO iel = 1,nels_pp
        DO j = 1,nod
          IF (fixed_node(i)==g_num_pp(j,iel)) THEN
            numfix_pp = numfix_pp + 1
            fixelem_pp(numfix_pp) = iel
            fixdof_pp(numfix_pp) = (j-1)*ndim + fixed_dof(i)
            fixval_pp(numfix_pp) = fixed_value(i)
          END IF
        END DO
      END DO
    END DO

    DEALLOCATE(fixed_node, fixed_dof, fixed_value)

    fixvalpiece_pp(:) = fixval_pp(:)/FLOAT(num_load_steps)

  END IF

  !-------------------------------------------------------------------------
  ! 5. Read and distribute natural boundary conditions
  !-------------------------------------------------------------------------

  ALLOCATE(fextpiece_pp(0:neq_pp))
  ALLOCATE(fext_pp(0:neq_pp))
  fext_pp = 0._iwp

  IF (loaded_nodes>0) THEN

    CALL READ_LOADS(fname_base,numpe,load_node,load_value)

    CALL LOAD_2(nn_start,g_num_pp,load_node,load_value,nf_pp,fext_pp(1:))
!   Originally CALL LOAD - but this differs from subroutine in ParaFEM
   
    DEALLOCATE(load_node,load_value)

  END IF

  fextpiece_pp(1:) = fext_pp(1:)/FLOAT(num_load_steps)

  !-------------------------------------------------------------------------
  ! 6. Allocate arrays
  !-------------------------------------------------------------------------

  ALLOCATE(r_pp(0:neq_pp), xnew_pp(0:neq_pp), diag_precon_pp(0:neq_pp))
  ALLOCATE(res_pp(0:neq_pp), deltax_pp(0:neq_pp), fint_pp(0:neq_pp))

  ALLOCATE(diag_precon_tmp(ntot,nels_pp))

  !-------------------------------------------------------------------------
  ! 7. Initialise the solution vector to 0.0
  !-------------------------------------------------------------------------

  xnew_pp = 0._iwp

  timest(3) = ELAP_TIME()

  ! Vector comp to compute F (gradient of deformation)
  DO i = 1,nod
    DO j = 1,ndim
      comp(i,j) = (i-1)*ndim + j
    END DO
  END DO

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!----------------------- Start Load LOOP --------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  DO iload = 1,num_load_steps
    converged = .FALSE.
    fext_pp(1:) = FLOAT(iload)*fextpiece_pp(1:)
    IF (numfix_pp>0) THEN
      fixval_pp(:) = FLOAT(iload)*fixvalpiece_pp(:)
    END IF
    
!------------------------------------------------------------------------------
!----------------------- Start Newton-Raphson iterations ----------------------
!------------------------------------------------------------------------------
    inewton = 0
    iterations: DO
      inewton = inewton + 1

      storefint_pp = 0._iwp

      CALL GATHER(xnew_pp(1:),xnewel_pp)
      IF (inewton>1 .and. numfix_pp>0) THEN
        DO i = 1,numfix_pp
          xnewel_pp(fixdof_pp(i),fixelem_pp(i)) = fixval_pp(i)
        END DO
      END IF

      IF (iload>1 .and. inewton==1 .and. numfix_pp>0) THEN
        DO i = 1,numfix_pp
          xnewel_pp(fixdof_pp(i),fixelem_pp(i))=FLOAT(iload-1)*fixvalpiece_pp(i)
        END DO
      END IF

      DO iel = 1,nels_pp
        kmat_elem = 0._iwp
        kgeo_elem = 0._iwp
        DO i = 1,nod
          num(i) = g_num_pp(i,iel) - nn_start + 1
        END DO
        coord = TRANSPOSE(g_coord_pp(:,num))
        auxm(:,1) = xnewel_pp(comp(:,1),iel)
        auxm(:,2) = xnewel_pp(comp(:,2),iel)
        auxm(:,3) = xnewel_pp(comp(:,3),iel)
        DO igauss = 1,nip
          CALL kine3D(igauss,auxm,coord,points,det,detF,beeF,defE,derivF,jacF)
          CALL venantkirchhoff(defE,e,v,piolaS,cmat)
          CALL push2r(detF,jacF,piolaS,sigma)
          CALL push4r(detF,jacF,cmat,cspa)
          CALL voigt2to1(sigma,sigma1C)
          CALL voigt4to2(cspa,deeF)
           
          storefint_pp(:,iel) = storefint_pp(:,iel) +         &
             MATMUL(TRANSPOSE(beeF),sigma1C)*det*detF*weights(igauss)
          
          kmat_elem(:,:) = kmat_elem(:,:) +     &
             MATMUL(MATMUL(TRANSPOSE(beeF),deeF),beeF)*det*detF*weights(igauss)

          geomH = MATMUL(MATMUL(TRANSPOSE(derivF),sigma),derivF)
          DO i = 1,dimH
            DO j = 1,dimH
              kgeo_elem(3*i-2,3*j-2) = kgeo_elem(3*i-2,3*j-2) + & 
                  geomH(i,j)*det*detF*weights(igauss)
              kgeo_elem(3*i-1,3*j-1) = kgeo_elem(3*i-1,3*j-1) + & 
                  geomH(i,j)*det*detF*weights(igauss)
              kgeo_elem(3*i,3*j) = kgeo_elem(3*i,3*j)         + & 
                  geomH(i,j)*det*detF*weights(igauss)
            END DO
          END DO
        END DO
        storekm_pp(:,:,iel) = kmat_elem(:,:) + kgeo_elem(:,:)
      END DO


      IF (inewton==1 .and. numfix_pp>0) THEN
        DO i = 1,numfix_pp
	      DO j = 1,ntot
	        IF (g_g_pp(j,fixelem_pp(i))>0) THEN
              storefint_pp(j,fixelem_pp(i)) =            &
              storefint_pp(j,fixelem_pp(i)) +            &
              fixvalpiece_pp(i)*storekm_pp(j,fixdof_pp(i),fixelem_pp(i))
	 	    END IF
		  END DO
	    END DO
	  END IF
		  
      fint_pp(:) = .0_iwp
      CALL SCATTER(fint_pp(1:),storefint_pp)

      r_pp(1:) = fext_pp(1:) - fint_pp(1:)       !Residual
      r_pp(0) = .0_iwp
	  
! Compute maxdiff of residual 
      maxdiff =  MAXABSVAL_P(r_pp(1:))

! Normalise residual vector and stiffness matrix for pcg
      IF (maxdiff == 0.0) THEN
        EXIT
      END IF

      diag_precon_tmp = .0_iwp
      DO iel = 1,nels_pp
        DO k = 1,ntot 
          diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel) + storekm_pp(k,k,iel)
        END DO
      END DO
  
!     Input: diag_precon_tmp r(ntot,nels_pp): Diagonal preconditioner at 
!                                             element level
      diag_precon_pp(:) = .0_iwp
      CALL SCATTER(diag_precon_pp(1:),diag_precon_tmp)
!     Output: diag_precon_pp r(1:neq_pp) Diagonal preconditioner assembled


!--------------------invert the preconditioner --------------------------------
      diag_precon_pp(1:) = 1._iwp/diag_precon_pp(1:)
      diag_precon_pp(0)  = .0_iwp

!---------------------------------------------------------------------------
!------------------------------- Solve using PCG ---------------------------
!---------------------------------------------------------------------------
      deltax_pp = .0_iwp
      res_pp    = r_pp

      CALL PCG_VER1(inewton,limit,tol,storekm_pp,r_pp(1:), &
                    diag_precon_pp(1:),rn0,deltax_pp(1:),iters)

      WRITE(91,*)iload,inewton,iters
      CALL FLUSH(91)

      xnew_pp(1:) = xnew_pp(1:) + deltax_pp(1:)
      xnew_pp(0) = .0_iwp

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

! Check convergence for Newton-Raphson iterations 
      energy = ABS(DOT_PRODUCT_P(res_pp(1:),deltax_pp(1:)))
      IF (inewton==1) THEN
       energy1 = energy
      END IF
      WRITE(*,2000)iload,inewton,energy,energy/energy1
      IF (inewton>1) THEN
        IF ((energy/energy1)<=tol2) THEN
          converged = .TRUE.
        END IF 
      END IF 
      IF(converged .OR. inewton==100) THEN
        EXIT
      END IF

    END DO iterations
!------------------------------------------------------------------------------
!------------------------- End Newton-Raphson iterations ----------------------
!------------------------------------------------------------------------------
    IF (numpe==npes) THEN
      WRITE(11,'(a,i3,a,f12.4,a,i4,a)') "Time after load step ",iload,": ", &
      ELAP_TIME() - timest(1),"     (",inewton," iterations )"
    END IF

    IF (iload==1) THEN
      ALLOCATE(disp_pp(ndim,nn_pp))
    END IF

    IF (iload==num_load_steps) THEN
      DEALLOCATE(diag_precon_tmp)
    END IF

!----------------------------------------------------------------------------
!-----------------------------print out results -----------------------------
!----------------------------------------------------------------------------
    IF (numpe==1) THEN
      IF (iload==1) THEN
        fname = fname_base(1:INDEX(fname_base, " ")-1) // ".dis"
        OPEN(24, file=fname, status='replace', action='write')
        fname = fname_base(1:INDEX(fname_base, " ")-1) // ".str"
        OPEN(25, file=fname, status='replace', action='write')
        fname = fname_base(1:INDEX(fname_base, " ")-1) // ".pri"
        OPEN(26, file=fname, status='replace', action='write')
        fname = fname_base(1:INDEX(fname_base, " ")-1) // ".rea"
        OPEN(27, file=fname, status='replace', action='write')
      END IF
    END IF

!-----print out displacements, stress, principal stress and reactions -------
    IF (MOD(iload,jump)==0) THEN
      
	  writetimes = writetimes + 1
      IF(timewrite) THEN
	    timest(4) = ELAP_TIME( )
	  END IF

      ALLOCATE(xnewnodes_pp(nodes_pp*nodof))
	  ALLOCATE(shape_integral_pp(nod,nels_pp))
	  ALLOCATE(stress_integral_pp(nod*nst,nels_pp))
	  ALLOCATE(stressnodes_pp(nodes_pp*nst))
	  ALLOCATE(principal_integral_pp(nod*nodof,nels_pp))
	  ALLOCATE(princinodes_pp(nodes_pp*nodof))
	  ALLOCATE(reacnodes_pp(nodes_pp*nodof))

      CALL GATHER(xnew_pp(1:),xnewel_pp)
      IF (numfix_pp > 0) THEN
        DO i = 1,numfix_pp
          xnewel_pp(fixdof_pp(i),fixelem_pp(i)) = fixval_pp(i)
        END DO
      END IF

      storefint_pp=0._iwp

      shape_integral_pp     = 0._iwp
      stress_integral_pp    = 0._iwp
      principal_integral_pp = 0._iwp
      DO iel = 1, nels_pp
        DO i = 1, nod
          num(i) = g_num_pp(i,iel) - nn_start + 1
        END DO
        coord = TRANSPOSE(g_coord_pp(:,num))
        auxm(:,1) = xnewel_pp(comp(:,1),iel)
        auxm(:,2) = xnewel_pp(comp(:,2),iel)
        auxm(:,3) = xnewel_pp(comp(:,3),iel)
        DO igauss = 1,nip
          CALL kine3D(igauss,auxm,coord,points,det,detF,beeF,defE,derivF,jacF)
          CALL venantkirchhoff(defE,e,v,piolaS,cmat)
          CALL push2r(detF,jacF,piolaS,sigma)
          CALL voigt2to1(sigma,sigma1C)
          CALL principalstress3D(sigma1C,principal)
          storefint_pp(:,iel) = storefint_pp(:,iel) +         &
                    MATMUL(TRANSPOSE(beeF),sigma1C)*det*detF*weights(igauss)
        
          CALL SHAPE_FUNCTIONS(igauss,points,value_shape)
        
          DO i = 1,nod
            idx1 = (i-1)*nst 
            idx2 = (i-1)*nodof
            shape_integral_pp(i,iel) = shape_integral_pp(i,iel) + &
                                       value_shape(i)*det*weights(igauss)
            DO j = 1,nst
              stress_integral_pp(idx1+j,iel) = stress_integral_pp(idx1+j,iel) +&
                        value_shape(i)*sigma1C(j)*det*weights(igauss)
            END DO
            DO j = 1,nodof
              principal_integral_pp(idx2+j,iel) =               &
                        principal_integral_pp(idx2+j,iel) +     &
                        value_shape(i)*principal(j)*det*weights(igauss)
            END DO
          END DO
        END DO
      END DO

      text = "*DISPLACEMENT"
      CALL SCATTER_NODES(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp, &
              node_start,node_end,xnewel_pp,xnewnodes_pp,1)
      CALL WRITE_NODAL_VARIABLE(text,24,iload,nodes_pp,npes,numpe,nodof, &
                                xnewnodes_pp)
	  DEALLOCATE(xnewnodes_pp)

      text = "*STRESS"
      CALL NODAL_PROJECTION(npes,nn,nels_pp,g_num_pp,nod,nst,nodes_pp,  &
       node_start,node_end,shape_integral_pp,stress_integral_pp,stressnodes_pp)
      CALL WRITE_NODAL_VARIABLE(text,25,iload,nodes_pp,npes,numpe,nst,   &
                                stressnodes_pp)
      DEALLOCATE(stress_integral_pp,stressnodes_pp)

      text = "*PRINCIPAL STRESS"
      CALL NODAL_PROJECTION(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp, &
     node_start,node_end,shape_integral_pp,principal_integral_pp,princinodes_pp)
      CALL WRITE_NODAL_VARIABLE(text,26,iload,nodes_pp,npes,numpe,nodof, &
                                princinodes_pp)
      DEALLOCATE(principal_integral_pp,princinodes_pp)
	 
      text = "*NODAL REACTIONS"
      CALL SCATTER_NODES(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp, &
              node_start,node_end,storefint_pp,reacnodes_pp,0)
      CALL WRITE_NODAL_VARIABLE(text,27,iload,nodes_pp,npes,numpe,nodof, &
                                reacnodes_pp)
      DEALLOCATE(reacnodes_pp)
	  
	  DEALLOCATE(shape_integral_pp)

      IF(timewrite) THEN
        timest(5) = elap_time( )
        timewrite = .FALSE.
      END IF

    END IF  !printing
  END DO !iload

  IF(numpe==1) THEN
    CLOSE(24)
    CLOSE(25)
    CLOSE(26)
    CLOSE(27)
  END IF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------- End Load LOOP --------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  IF (numpe==npes) THEN
    WRITE(11,'(a,i5,a)') "This job ran on ",npes," processors"
    WRITE(11,'(A,3(I8,A))')"There are ",nn," nodes",nels," elements and ",&
                           neq," equations"
    WRITE(11,*) "Time after the mesh   :", timest(2) - timest(1)
    WRITE(11,*) "Time after boundary conditions  :", timest(3) - timest(1)
    WRITE(11,*) "This analysis took  :", elap_time( ) - timest(1)
    WRITE(11,*) "Time to write results (each time) :", timest(5) - timest(4)
    WRITE(11,*) "Time inside the load loop  :", elap_time( ) - timest(3) - &
                                        writetimes*(timest(5)-timest(4))
    CALL FLUSH(11)
    CLOSE(11)
  END IF

!   Formats
  2000 format(' Energy  ',i3,1p,i3,1p,e25.15,1p,e25.15) 


!---------------------------------- shutdown ----------------------------------
  CALL SHUTDOWN()

 END PROGRAM XX7
