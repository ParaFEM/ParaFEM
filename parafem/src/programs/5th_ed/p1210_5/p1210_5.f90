PROGRAM p1210      
!-------------------------------------------------------------------------
!  Program 12.10 forced vibration of an elastic-plastic(Von Mises) solid
!  Viscoplastic strain method, lumped mass, explicit integration
!-------------------------------------------------------------------------
 USE precision; USE global_variables; USE mp_interface; USE input
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE elements; USE steering; USE plasticity 
 IMPLICIT NONE
! neq,ntot are now global variables - not declared 
 INTEGER,PARAMETER::nodof=3,ndim=3,nst=6
 INTEGER::nn,nr,nip,loaded_nodes,nres=1,nod,i,j,k,ii,jj,iel,nstep,npri,  &
   num_no,no_index_start,is,it,nlen,ndof,nels,npes_pp,node_end,          &
   node_start,nodes_pp,argc,iargc,meshgen,partitioner=1
 REAL(iwp)::rho,dtim,e,v,det,sbary,sigm,f,fnew,fac,volume,sbar,          &
   dsbar,lode_theta,real_time,tload,pload
 REAL(iwp),PARAMETER::zero=0.0_iwp    
 CHARACTER(LEN=15)::element,io_type; CHARACTER(LEN=50)::argv
 CHARACTER(LEN=6)::ch
!------------------------- dynamic arrays --------------------------------
 REAL(iwp),ALLOCATABLE::points(:,:),bdylds_pp(:),x1_pp(:),d1x1_pp(:),    &
   stressv(:),pl(:,:),emm(:),d2x1_pp(:),tensor_pp(:,:,:),                &
   etensor_pp(:,:,:),val(:,:),dee(:,:),mm_pp(:),jac(:,:),weights(:),     &
   der(:,:),deriv(:,:),bee(:,:),eld(:),eps(:),sigma(:),bload(:),         &
   eload(:),mm_tmp(:,:),g_coord_pp(:,:,:),pmul_pp(:,:),utemp_pp(:,:),    &
   disp_pp(:),fext_pp(:),temp(:)
 INTEGER,ALLOCATABLE::timest(:),rest(:,:),no(:),no_local(:),node(:),     &
   g_num_pp(:,:),g_g_pp(:,:),no_local_temp(:)
!----------------------- input and initialisation ------------------------
 ALLOCATE(timest(20)); timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen)
 CALL read_p1210(argv,numpe,dtim,e,element,loaded_nodes,meshgen,nels,    &
   nip,nn,nod,npri,nr,nres,nstep,partitioner,pload,rho,sbary,v)
 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ndof=nod*nodof; ntot=ndof
 ALLOCATE(g_num_pp(nod, nels_pp),g_coord_pp(nod,ndim,nels_pp),           & 
   rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
 CALL read_g_num_pp(argv,iel_start,nels,nn,numpe,g_num_pp)
 IF(meshgen == 2) THEN
   CALL abaqus2sg(element,g_num_pp); PRINT *, "Renumbered nodes"
 END IF
 CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
 CALL read_rest(argv,numpe,rest)
 ALLOCATE(points(nip,ndim),weights(nip),dee(nst,nst),                    &
   pmul_pp(ntot,nels_pp),tensor_pp(nst,nip,nels_pp),no(loaded_nodes),    &
   pl(nst,nst),etensor_pp(nst,nip,nels_pp),jac(ndim,ndim),der(ndim,nod), &
   deriv(ndim,nod),bee(nst,ntot),eld(ntot),eps(nst),sigma(nst),emm(ntot),&
   bload(ntot),eload(ntot),stressv(nst),g_g_pp(ntot,nels_pp),            &
   mm_tmp(ntot,nels_pp),utemp_pp(ntot,nels_pp))
!----------  find the steering array and equations per process -----------
 CALL rearrange(rest); g_g_pp=0; neq=0
 elements_0: DO iel=1,nels_pp
   CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest)
 END DO elements_0
 neq=MAXVAL(g_g_pp); neq=MAX_INTEGER_P(neq); CALL calc_neq_pp
 CALL calc_npes_pp(npes,npes_pp); CALL make_ggl(npes_pp,npes,g_g_pp)
 DO i=1,neq_pp; IF(nres==ieq_start+i-1)THEN;it=numpe;is=i;END IF;END DO
 IF(numpe==it) THEN
   OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
   WRITE(11,'(A,I6,A)') "This job ran on ",npes," processes"
   WRITE(11,'(A,3(I12,A))') "There are ",nn," nodes ",nr,                &
     " restrained and ", neq," equations"
   WRITE(11,'(A,F10.4)') "Time after setup was:",elap_time()-timest(1)
 END IF 
 CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
 ALLOCATE(disp_pp(nodes_pp*ndim),temp(nodes_pp)); disp_pp=zero; temp=0
 ALLOCATE(bdylds_pp(neq_pp),x1_pp(neq_pp),d1x1_pp(neq_pp),mm_pp(neq_pp), &     
   d2x1_pp(neq_pp),fext_pp(neq_pp)); bdylds_pp=zero; x1_pp=zero
   d1x1_pp=zero; d2x1_pp=zero; mm_pp=zero; fext_pp=zero
!-------------------- calculate diagonal mass matrix ---------------------
 mm_tmp=zero; CALL sample(element,points,weights)
 elements_1: DO iel=1,nels_pp
   volume=zero
   gauss_pts_1: DO i=1,nip
     CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
     det=determinant(jac); volume=volume+det*weights(i)*rho
   END DO gauss_pts_1
   emm=volume/13._iwp; emm(1:19:6)=emm(4)*.125_iwp
   emm(2:20:6)=emm(4)*.125_iwp; emm(3:21:6)=emm(4)*.125_iwp
   emm(37:55:6)=emm(4)*.125_iwp; emm(38:56:6)=emm(4)*.125_iwp
   emm(39:57:6)=emm(4)*.125_iwp; mm_tmp(:,iel)=mm_tmp(:,iel)+emm 
 END DO elements_1
 CALL scatter(mm_pp,mm_tmp); DEALLOCATE(mm_tmp)
!----------------------------------- loads -------------------------------
 IF(loaded_nodes>0) THEN
   ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes)); val=zero; node=0
   CALL read_loads(argv,numpe,node,val)
   CALL load(g_g_pp,g_num_pp,node,val,fext_pp(1:))
   tload = SUM_P(fext_pp(1:)); DEALLOCATE(node,val)
 END IF
!---------------------- explicit integration loop ------------------------
 tensor_pp=zero; etensor_pp=zero; real_time=zero
 time_steps: DO jj=1,nstep
    real_time=real_time+dtim; bdylds_pp=zero
    x1_pp=x1_pp+(d1x1_pp+d2x1_pp*dtim*.5_iwp)*dtim
!------------------ element stress-strain relationship -------------------
    pmul_pp=zero; utemp_pp=zero; CALL gather(x1_pp,pmul_pp)
    elements_2: DO iel=1,nels_pp          
      bload=zero; eld=pmul_pp(:,iel)
      gauss_pts_2: DO i=1,nip
        dee=zero; CALL deemat(e,v,dee); CALL shape_der(der,points,i)
        jac=MATMUL(der,g_coord_pp(:,:,iel)); det=determinant(jac)
        CALL invert(jac); deriv=MATMUL(jac,der); CALL beemat(deriv,bee)
        eps=MATMUL(bee,pmul_pp(:,iel)); eps=eps-etensor_pp(:,i,iel)
        sigma=MATMUL(dee,eps); stressv=sigma+tensor_pp(:,i,iel)
        CALL invar(stressv,sigm,dsbar,lode_theta); fnew=dsbar-sbary
!------------------ check whether yield is violated -----------------------
        IF(fnew>=.0_iwp)THEN
          stressv=tensor_pp(:,i,iel)
          CALL invar(stressv,sigm,sbar,lode_theta)
          f=sbar-sbary; fac=fnew/(fnew-f)
          stressv=tensor_pp(:,i,iel)+(1._iwp-fac)*sigma
          CALL vmpl(e,v,stressv,pl); dee=dee-fac*pl
        END IF
        sigma=MATMUL(dee,eps); sigma=sigma+tensor_pp(:,i,iel)
        eload=MATMUL(sigma,bee); bload=bload+eload*det*weights(i)
!--------------------- update the gauss points ---------------------------
        tensor_pp(:,i,iel)=sigma
        etensor_pp(:,i,iel)=etensor_pp(:,i,iel)+eps
      END DO gauss_pts_2; utemp_pp(:,iel)=utemp_pp(:,iel)-bload    
    END DO elements_2
    CALL scatter(bdylds_pp,utemp_pp); bdylds_pp=bdylds_pp+fext_pp*pload
    bdylds_pp=bdylds_pp/mm_pp
    d1x1_pp=d1x1_pp+(d2x1_pp+bdylds_pp)*.5_iwp*dtim; d2x1_pp=bdylds_pp
!---------------------- output displacements -----------------------------
    IF(jj==jj/npri*npri) THEN 
      IF(numpe==it) THEN
        WRITE(11,'(4E12.4)')real_time,x1_pp(is),d1x1_pp(is),d2x1_pp(is)
      END IF
      IF(numpe==1) THEN;  WRITE(ch,'(I6.6)') jj
        OPEN(12,file=argv(1:nlen)//".ensi.DISPL-"//ch,status='replace',  &
             action='write'); WRITE(12,'(A)')                            &
        "Alya Ensight Gold --- Vector per-node variable file"
        WRITE(12,'(A/A/A)') "part", "     1","coordinates"
      END IF; disp_pp=zero; utemp_pp=zero
      CALL gather(x1_pp(1:),utemp_pp)
      CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,     &
                         node_start,node_end,utemp_pp,disp_pp,1)
      DO i=1,ndim ; temp=zero
        DO ii=1,nodes_pp; k=i+(ndim*(ii-1)); temp(ii)=disp_pp(k); END DO
        CALL dismsh_ensi_p(12,jj,nodes_pp,npes,numpe,1,temp)
      END DO ; IF(numpe==1) CLOSE(12)
    END IF
  END DO time_steps
  IF(numpe==it) THEN
    WRITE(11,'(A,F10.4)')"This analysis took  :",elap_time()-timest(1)
  END IF; CALL shutdown() 
END PROGRAM p1210
