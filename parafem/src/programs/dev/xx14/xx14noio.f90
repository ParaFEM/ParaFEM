!$Id$
PROGRAM xx14noio
!-----------------------------------------------------------------------
! Anton Shterenlikht (University of Bristol)
! Luis Cebamanos (EPCC, University of Edinburgh)
!
! Program xx14 - linking ParaFEM with CGPACK, specifically
! modifying p121 from 5th edition to link with the cgca module.
!
! 12.1 is a three dimensional analysis of an elastic solid
! using 20-node brick elements, preconditioned conjugate gradient
! solver; diagonal preconditioner diag_precon; parallel version
! loaded_nodes only
!
! This version uses Cray extensions, which are all in TS 18508
! "Additional Parallel Features in Fortran", WG5/N2074.
! I think it's just CO_SUM.
!
! The CA CS is aligned with the FE CS, but not the same origin.
! Anyway, all cells are within the FE model, so mapping is easy.
!-----------------------------------------------------------------------
!USE mpi_wrapper  !remove comment for serial compilation

!*** CGPACK part *****************************************************72
! The CGPACK module must be used
use cgca
!*** end CGPACK part *************************************************72

USE precision
USE global_variables
USE mp_interface
USE input
USE output
USE loading
USE timing
USE maths
USE gather_scatter
USE steering
USE new_library

IMPLICIT NONE

! neq, ntot are now global variables - must not be declared

INTEGER, PARAMETER :: nodof = 3, ndim = 3, nst = 6
INTEGER :: loaded_nodes, iel, i, j, k, iters, limit, nn, nr, nip, nod, &
   nels, ndof, npes_pp, node_end, node_start, nodes_pp, meshgen,       &
   partitioner, nlen

REAL( iwp ), PARAMETER :: zero=0.0_iwp
REAL( iwp ) :: e, v, det, tol, up, alpha, beta, q

LOGICAL :: converged=.false.

CHARACTER( LEN=50 ) :: argv
CHARACTER( LEN=15 ) :: element
CHARACTER( LEN=6 ) :: ch

!---------------------------- dynamic arrays -------------------------72
REAL( iwp ), ALLOCATABLE :: points(:,:), dee(:,:), weights(:),         &
   val(:,:),                                                           &
   disp_pp(:), g_coord_pp(:,:,:), jac(:,:), der(:,:), deriv(:,:),      &
   bee(:,:), storkm_pp(:,:,:), eps(:), sigma(:), diag_precon_pp(:),    &
   p_pp(:), r_pp(:), x_pp(:), xnew_pp(:), u_pp(:), pmul_pp(:,:),       &
   utemp_pp(:,:), d_pp(:), timest(:), diag_precon_tmp(:,:),            &
   eld_pp(:,:), temp(:), tot_r_pp(:)
INTEGER, ALLOCATABLE :: rest(:,:), g_num_pp(:,:), g_g_pp(:,:), node(:)

!*** CGPACK part *****************************************************72
! CGPACK parameters
integer, parameter :: cgca_linum=5 ! number of loading iterations
logical( kind=ldef ), parameter :: cgca_yesdebug = .true.,             &
 cgca_nodebug = .false.
real( kind=rdef ), parameter :: cgca_zero = 0.0_rdef,                  &
 cgca_one = 1.0_rdef,                                                  &
 ! cleavage stress on 100, 110, 111 planes for BCC,
 ! see the manual for derivation, GPa.
 cgca_scrit(3) = (/ 1.05e1_rdef, 1.25e1_rdef, 4.90e1_rdef /)

! CGPACK variables
integer( kind=idef ) ::                                                &
   cgca_ir(3),            & ! coarray codimensions
   cgca_img,              &
   cgca_nimgs,            &
   cgca_ng,               & ! number of grains in the whole model
   cgca_clvg_iter,        & ! number of cleavage iterations
!   cgca_lc(3),            & ! local coordinates of a cell with its image
!   cgca_lowr(3),          & ! local coordinates of the lower box corner
!   cgca_uppr(3),          & ! local coordinates of the upper box corner
!   cgca_iflag,            & ! 1 box in, 2 box out, 3 neither 
   cgca_liter               ! load iteration number
integer( kind=iarr ) :: cgca_c(3) ! coarray dimensions
integer( kind=iarr ), allocatable :: cgca_space(:,:,:,:) [:,:,:]

real( kind=rdef ) ::                                                   &
  cgca_qual,             & ! quality
  cgca_bsz(3),           & ! the given and the updated "box" size
  cgca_origin(3),        & ! origin of the "box" cs, in FE cloads
  cgca_rot(3,3),         & ! rotation tensor *from* FE cs *to* CA cs
  cgca_dm,               & ! mean grain size, linear dim, phys units
  cgca_res,              & ! resolutions, cells per grain
  cgca_bcol(3),          & ! lower phys. coords of the coarray on image
  cgca_bcou(3),          & ! upper phys. coords of the coarray on image
  cgca_stress(3,3),      & ! stress tensor
  cgca_length,           & ! fracture length scale
  cgca_time_inc,         & ! time increment
  cgca_lres,             & ! linear resolution, cells per unit of length
  cgca_charlen,          & ! characteristic element length
  cgca_fracvol             ! volume (number) of fractured cells per img
real( kind=rdef ), allocatable :: cgca_grt(:,:,:)[:,:,:]
logical( kind=ldef ) :: cgca_solid
! logical( kind=ldef ) :: cgca_lflag
character( len=6 ) :: cgca_citer
!*** end CGPACK part *************************************************72


!------------------------ input and initialisation ---------------------
ALLOCATE( timest( 20 ) )
timest = zero
timest( 1 ) = elap_time()

!*    Get the rank of the processes and the total number of processes
!* intent( out ):
!*          numpe - integer, process number (rank)
!*           npes - integer, total number of processes (size)
CALL find_pe_procs( numpe, npes )
!
! dummy
! lines
! to make
! the line numbers
! match with
! xx14std
!
!
!
!
!

!*    Returns the base name of the data file.
!* intent( out ):
!*           argv - character(*), data file base name
!*           nlen - integer, number of characters in data file base name
CALL getname( argv, nlen )

!*    Master processor reads the general data for the problem
!     and broadcasts it to the slave processors.
!  in:
!        argv - character, file name to read from
!       numpe - MPI rank
!           e - Young's modulus
!     element - character, element type
!       limit - max number of iterations
!loaded_nodes - number of nodes with applied forces
!     meshgen - mesh numbering scheme
!        nels - total number of elements
!         nip - number of Gauss points
!          nn - total number of nodes in the mesh
!         nod - number of nodes per element
!         tol - tolerance
!           v - Poisson's ratio
!https://code.google.com/p/parafem/source/browse/trunk/parafem/src/modules/mpi/input.f90
CALL read_p121( argv, numpe, e, element, limit, loaded_nodes, meshgen, &
                nels, nip, nn, nod, nr, partitioner, tol, v )

! Calculates the number of elements, nels_pp, assigned to each
! processor.
! It is effectively a very naive method of mesh partitioning.
! The subroutine also computes, iel_start, the first element number
! on each processor. iel_start and nels_pp are the external variables
! modified by this subroutine. Note that they are global variables
! that are not passed through the list of arguments.
!
! https://code.google.com/p/parafem/source/browse/trunk/parafem/src/modules/mpi/gather_scatter.f90
!
! nels_pp is indeed a global var, defined in
! https://code.google.com/p/parafem/source/browse/trunk/parafem/src/modules/shared/global_variables.f90
CALL calc_nels_pp( argv, nels, npes, numpe, partitioner, nels_pp )

!   nod - number of nodes per element
! nodof = 3 (see the beginning of this file), probably
!           degrees of freedom per node.
ndof = nod * nodof

! ntot - a global variable,
! the total number of degrees of freedom per element
ntot = ndof

! g_num_pp(nod,nels_pp) - integer, elements connectivity
! g_coord_pp - global coord?
!https://sourceforge.net/p/parafem/code/HEAD/tree/trunk/parafem/src/modules/mpi/input.f90
ALLOCATE( g_num_pp( nod, nels_pp ) )
allocate( g_coord_pp( nod, ndim, nels_pp ) )
allocate( rest( nr,nodof+1) )

g_num_pp = 0
g_coord_pp = zero
rest = 0

!https://sourceforge.net/p/parafem/code/HEAD/tree/trunk/parafem/src/modules/mpi/input.f90
CALL read_g_num_pp( argv, iel_start, nn, npes, numpe, g_num_pp )

IF ( meshgen == 2 ) CALL abaqus2sg( element, g_num_pp )
CALL read_g_coord_pp( argv, g_num_pp, nn, npes, numpe, g_coord_pp )
CALL read_rest( argv, numpe, rest )

timest(2) = elap_time()

ALLOCATE( points(nip,ndim) )
allocate( dee(nst,nst) )
allocate( jac(ndim,ndim) )
allocate( der(ndim,nod) )
allocate( deriv( ndim, nod ) )
allocate( bee( nst, ntot ) )
allocate( weights( nip ) )
allocate( eps( nst ) )
allocate( sigma(nst) )
allocate( storkm_pp( ntot, ntot, nels_pp ) )
allocate( pmul_pp( ntot, nels_pp ) )
allocate( utemp_pp(ntot,nels_pp) )
allocate( g_g_pp(ntot,nels_pp) )

IF ( numpe==1 ) THEN
   open(  11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='write' )
   write( 11,'(A,I7,A)') "This job ran on ",npes," processes"
   write( 11,'(A,3(I12,A))') "There are ",nn," nodes", nr, &
        " restrained and ",neq," equations"
   write( 11,'(A,F10.4)') "Time to read input is:", timest(2)-timest(1)
   write( 11,'(A,F10.4)') "Time after setup is:", elap_time()-timest(1)
END IF
!*** end of ParaFEM input and initialisation *************************72


!*** CGPACK part *****************************************************72
! *** CGPACK first executable statement ***
! In this test set the number of images via the env var,
! or simply as an argument to aprun.
! the code must be able to cope with any value >= 1.
  cgca_img = this_image()
cgca_nimgs = num_images()

! Initialise random number seed.
! Need to do this before cgca_pfem_cenc, where there order of comms
! is chosen at *random* to even the remote access pattern.
!
! Argument:
! .false. - no debug output
!  .true. - with debug output
call cgca_irs( .false. )

! dump CGPACK parameters and some ParaFEM settings
if ( cgca_img .eq. 1 ) then
  call cgca_pdmp
  write (*,*) "Young's mod:", e, "Poisson's ratio", v
end if

! physical dimensions of the box, must be the same
! units as in the ParaFEM.
! Must be fully within the FE model, which for xx14
! is a cube with lower bound at (0,0,-10), and the
! upper bound at (10,10,0)
cgca_bsz = (/ 10.0, 10.0, 10.0 /)

! Origin of the box cs, in the same units.
! This gives the upper limits of the box at 0+10=10, 0+10=10, -10+10=0
! all within the FE model.
cgca_origin = (/ 0.0, 0.0, -10.0 /)

! Rotation tensor *from* FE cs *to* CA cs.
! The box cs is aligned with the box.
cgca_rot         = cgca_zero
cgca_rot( 1, 1 ) = 1.0
cgca_rot( 2, 2 ) = 1.0
cgca_rot( 3, 3 ) = 1.0

! mean grain size, also mm
cgca_dm = 8.0e-1_rdef

! resolution, cells per grain
cgca_res = 1.0e5_rdef

! cgpack length scale, also in mm
! Equivalent to crack propagation distance per unit of time,
! i.e. per second. Let's say 1 km/s = 1.0e3 m/s = 1.0e6 mm/s. 
cgca_length = 1.0e6_rdef

! In p121_medium, each element is 0.25 x 0.25 x 0.25 mm, so
! the charlen must be bigger than that.
cgca_charlen = 0.4

! each image calculates the coarray grid dimensions
call cgca_gdim( cgca_nimgs, cgca_ir, cgca_qual )

! calculate the resolution and the actual phys dimensions of the box
! subroutine cgca_cadim( bsz, res, dm, ir, c, lres, ng )
call cgca_cadim( cgca_bsz, cgca_res, cgca_dm, cgca_ir, cgca_c,         &
                 cgca_lres, cgca_ng )

! dump some stats from img 1
if (cgca_img .eq. 1 ) then
  write ( *, "(9(a,i0),tr1,g0,tr1,g0,3(a,g0),a)" )                     &
    "img: ", cgca_img  , " nimgs: ", cgca_nimgs,                       &
     " ("  , cgca_c (1), ","       , cgca_c (2), ",", cgca_c (3),      &
     ")["  , cgca_ir(1), ","       , cgca_ir(2), ",", cgca_ir(3),      &
     "] "  , cgca_ng   ,                                               &
    cgca_qual, cgca_lres,                                              &
         " (", cgca_bsz(1), ",", cgca_bsz(2), ",", cgca_bsz(3), ")"
  write (*,*) "dataset sizes for ParaView", cgca_c*cgca_ir
  write (*,"(a,es10.2)") "Total cells in the model",                   &
              product( real(cgca_c) * real(cgca_ir) )
end if

! allocate space coarray with 2 layers
!subroutine cgca_as( l1, u1, l2, u2, l3, u3,                           &
!             col1, cou1, col2, cou2, col3, props, coarray )
call cgca_as( 1, cgca_c(1),  1, cgca_c(2),  1, cgca_c(3),              &
              1, cgca_ir(1), 1, cgca_ir(2), 1, 2, cgca_space )

! calculate the phys. dim. of the coarray on each image
!subroutine cgca_imco( space, lres, bcol, bcou )
call cgca_imco( cgca_space, cgca_lres, cgca_bcol, cgca_bcou )

! dump box lower and upper corners from every image
write ( *,"(a,i0,2(a,3(es9.2,tr1)))" ) "img ", cgca_img,               &
       " bcol: ", cgca_bcol, "bcou: ", cgca_bcou

! and now in FE cs:
!write ( *,"(a,i0,2(a,3(es9.2,tr1)),a)" ) "img: ", cgca_img,            &
!   " FE bcol: (",                                                      &
!    matmul( transpose( cgca_rot ),cgca_bcol ) + cgca_origin,           &
!  ") FE bcou: (",                                                      &
!    matmul( transpose( cgca_rot ),cgca_bcou ) + cgca_origin, ")"

! confirm that image number .eq. MPI process number
write (*,*) "img",cgca_img," <-> MPI proc", numpe

! allocate the tmp centroids array: cgca_pfem_centroid_tmp%r ,
! an allocatable array component of a coarray variable of derived type
call cgca_pfem_ctalloc( ndim, nels_pp )

! set the centroids array component on this image, no remote calls.
! first dim - coord, 1,2,3
! second dim - element number, always starting from 1
! g_coord_pp is allocated as g_coord_pp( nod, ndim, nels_pp )
cgca_pfem_centroid_tmp%r = sum( g_coord_pp(:,:,:), dim=1 ) / nod

         ! set cgca_pfem_centroid_tmp[*]%r
sync all ! must add execution segment
         ! use cgca_pfem_centroid_tmp[*]%r

! Set lcentr private arrays on every image. Choose one of the two
! routines that do this:
! - cgca_pfem_cenc - uses all-to-all algorithm.
! - cgca_pfem_map  - uses CO_SUM, CO_MAX and *large* tmp arrays
! Both routines have identical sets of input arguments.
call cgca_pfem_cenc( cgca_origin, cgca_rot, cgca_bcol, cgca_bcou )
!call cgca_pfem_map( cgca_origin, cgca_rot, cgca_bcol, cgca_bcou )

! Dump lcentr for debug
! *** a lot *** of data
!call cgca_pfem_lcentr_dump

! Allocate cgca_pfem_integrity%i(:), array component of a coarray of
! derived type. Allocating *local* array.
! i is set to 1.0 on allocation.
call cgca_pfem_integalloc( nels_pp )
 
! Allocate the Young's modulus 2D array
call cgca_pfem_ealloc( nip, nels_pp )
  
! initially set the Young's modulus to "e" everywhere
cgca_pfem_enew = e

! Dump the FE centroids in CA cs to stdout.
! Obviously, this is an optional step, just for debug.
! Can be called from any or all images.
! call cgca_pfem_cendmp ! NO DUMP IN THIS VERSION - TAKES TOO LONG!!!

! Generate microstructure

! allocate rotation tensors
call cgca_art( 1, cgca_ng, 1, cgca_ir(1), 1, cgca_ir(2), 1, cgca_grt )

! initialise space
cgca_space( :, :, :, cgca_state_type_grain ) = cgca_liquid_state
cgca_space( :, :, :, cgca_state_type_frac  ) = cgca_intact_state

! nuclei, sync all inside
! last argument:
! .false. - no debug output
!  .true. - with debug output
call cgca_nr( cgca_space, cgca_ng, .false. )

! assign rotation tensors, sync all inside
call cgca_rt( cgca_grt )

! solidify
!subroutine cgca_sld( coarray, periodicbc, iter, heartbeat, solid )
! second argument:
!  .true. - periodic BC
! .false. - no periodic BC
! ===>>> implicit sync all inside <<<===
call cgca_sld( cgca_space, .false., 0, 10, cgca_solid )

! initiate grain boundaries
call cgca_igb( cgca_space )

! smoothen the GB, several iterations,
! halo exchange,
! sync needed following halo exchange
call cgca_gbs( cgca_space )
call cgca_hxi( cgca_space )

sync all

call cgca_gbs( cgca_space )
call cgca_hxi( cgca_space )

sync all

! update grain connectivity, local routine, no sync needed
call cgca_gcu( cgca_space )

! set a single crack nucleus in the centre of the x1=x2=0 edge
if ( cgca_img .eq. 1 ) then
  cgca_space( 1, 1, cgca_c(3)/2, cgca_state_type_frac )                &
              [ 1, 1, cgca_ir(3)/2 ] = cgca_clvg_state_100_edge
end if

         ! cgca_space changed locally on every image
sync all !
         ! cgca_space used

! Now can deallocate the temp array cgca_pfem_centroid_tmp%r.
! Could've done this earlier, but best to wait until sync all is
! required, to avoid extra sync.
call cgca_pfem_ctdalloc

! img 1 dumps space arrays to files
!if ( cgca_img .eq. 1 ) write (*,*) "dumping model to files"
!call cgca_pswci( cgca_space, cgca_state_type_grain, "zg0.raw" )
!call cgca_pswci( cgca_space, cgca_state_type_frac,  "zf0.raw" )
!if ( cgca_img .eq. 1 ) write (*,*) "finished dumping model to files"

! debug
!sync all
!stop

! allocate the stress array component of cgca_pfem_stress coarray
! subroutine cgca_pfem_salloc( nels_pp, intp, comp )
call cgca_pfem_salloc( nels_pp, nip, nst )

! just in case, not sure what is coming next...
sync all
!*** end CGPACK part *************************************************72


  !----------  find the steering array and equations per process ---------
  CALL rearrange(rest)
  g_g_pp = 0
  neq = 0

  ! map equations to elements
  ! http://parafem.googlecode.com/svn/trunk/parafem/src/modules/shared/new_library.f90
  elements_0: DO iel=1,nels_pp
     CALL find_g3( g_num_pp(:,iel), g_g_pp(:,iel), rest )
  END DO elements_0

  neq = MAXVAL(g_g_pp)
  neq = max_p(neq)
  CALL calc_neq_pp
  CALL calc_npes_pp( npes, npes_pp )
  CALL make_ggl( npes_pp, npes, g_g_pp )

  ALLOCATE( p_pp(neq_pp), source=zero )
  ! loads in the current increment, a fraction of tot_r_pp
  allocate( r_pp(neq_pp), source=zero )
  ! total loads
  allocate( tot_r_pp(neq_pp), source=zero )
  allocate( x_pp(neq_pp), source=zero )
  allocate( xnew_pp(neq_pp), source=zero )
  allocate( u_pp(neq_pp), source=zero )
  allocate( d_pp(neq_pp), source=zero )
  allocate( diag_precon_pp(neq_pp), source=zero )

  ! do this only once per program
  IF ( loaded_nodes > 0 ) THEN
     ALLOCATE( node(loaded_nodes), source=0 )
     allocate( val(ndim, loaded_nodes), source=zero )
     CALL read_loads( argv, numpe, node, val )
     CALL load( g_g_pp, g_num_pp, node, val, tot_r_pp(1:) )
     DEALLOCATE( node )
     deallocate( val )
  end if
  DEALLOCATE(g_g_pp)

  ! allocate arrays outside of the loop
  ALLOCATE( eld_pp(ntot,nels_pp), source=zero )

  ! Since the Young's modulus can be updated,
  ! all below must be inside the loading iterations!
  ! Make sure not to allocate/deallocate arrays within the loop
  load_iter: do cgca_liter=1, cgca_linum

     !------ element stiffness integration and build the preconditioner ---72

     ! deemat needs to go into a loop elements_1
     ! need an element array of E

     dee = zero

     CALL sample( element, points, weights )
     storkm_pp = zero

     elements_1: DO iel=1,nels_pp
        gauss_pts_1: DO i=1,nip
           ! from cgca_m3pfem.f90:
           ! allocate( cgca_pfem_enew( nip, nels_pp ) )
           CALL deemat( dee, cgca_pfem_enew(i,iel), v )
           CALL shape_der( der, points, i )
           jac = MATMUL( der, g_coord_pp(:,:,iel) )
           det = determinant(jac)
           CALL invert(jac)
           deriv = MATMUL( jac, der )
           CALL beemat( bee, deriv )
           storkm_pp(:,:,iel) = storkm_pp(:,:,iel) +                             &
                MATMUL( MATMUL( TRANSPOSE(bee), dee ), bee ) * det * weights(i)
        END DO gauss_pts_1
     END DO elements_1

     ALLOCATE( diag_precon_tmp( ntot, nels_pp ), source=zero )

     elements_2: DO iel=1,nels_pp
        DO i=1,ndof
           diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel) + storkm_pp(i,i,iel)
        END DO
     END DO elements_2

     CALL scatter( diag_precon_pp, diag_precon_tmp )
     DEALLOCATE( diag_precon_tmp )

     !----------------------------- get starting r --------------------------

     ! loading increases from 1/cgca_linum of tot_r_pp
     ! to tot_r_pp
     r_pp(1:) = tot_r_pp(1:) * cgca_liter/cgca_linum

     IF ( loaded_nodes > 0 ) THEN
        q = SUM_P( r_pp(1:) )
        IF ( numpe==1 ) then
           write (11, *) "Load iter:", cgca_liter
           write (11,'(A,E12.4)') "The total load is:", q
        end if
     END IF

     diag_precon_pp = 1._iwp/diag_precon_pp
     d_pp = diag_precon_pp*r_pp
     p_pp = d_pp
     x_pp = zero

     !--------------------- preconditioned cg iterations --------------------
     iters=0
     timest(3) = elap_time()

     iterations: DO
        iters = iters + 1
        u_pp = zero
        pmul_pp = zero
        utemp_pp = zero
        CALL gather(p_pp,pmul_pp)

        elements_3: DO iel=1,nels_pp
           utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
           !    CALL dgemv('n',ntot,ntot,1.0_iwp,storkm_pp(:,:,iel),ntot,         &
           !               pmul_pp(:,iel),1,0.0_iwp,utemp_pp(:,iel),1)
        END DO elements_3

        CALL scatter(u_pp,utemp_pp)

        !-------------------------- pcg equation solution ----------------------
        up = DOT_PRODUCT_P(r_pp,d_pp)
        alpha = up/DOT_PRODUCT_P(p_pp,u_pp)
        xnew_pp = x_pp+p_pp*alpha
        r_pp = r_pp-u_pp*alpha
        d_pp = diag_precon_pp*r_pp
        beta = DOT_PRODUCT_P(r_pp,d_pp)/up
        p_pp = d_pp+p_pp*beta
        CALL checon_par(xnew_pp,tol,converged,x_pp)
        IF ( converged .OR. iters == limit ) EXIT
     END DO iterations

     IF ( numpe==1 ) THEN
        write(11,'(A,I6)')"The number of iterations to convergence was ",iters
        write(11,'(A,F10.4)')"Time to solve equations was  :",                &
             elap_time()-timest(3)
        write(11,'(A,E12.4)')"The central nodal displacement is :",xnew_pp(1)
     END IF

     !--------------- recover stresses at centroidal gauss point ------------

     CALL gather(xnew_pp(1:),eld_pp)

     elmnts: DO iel = 1, nels_pp
        intpts: DO i = 1, nip
           call deemat(dee, cgca_pfem_enew( i, iel ), v)
! Compute the derivatives of the shape functions at a Gauss point.
! http://parafem.googlecode.com/svn/trunk/parafem/src/modules/shared/new_library.f90
           CALL shape_der(der,points,i)
           jac = MATMUL(der,g_coord_pp(:,:,iel))
           CALL invert(jac)
           deriv = MATMUL(jac,der)
           CALL beemat(bee,deriv)
           eps = MATMUL(bee,eld_pp(:,iel))
           sigma = MATMUL(dee,eps)
           !    write (*,*) "MPI rank", numpe, "el", iel,  &
           !                "int. point", i, "stress", sigma

           ! set the stress array on this image
           ! from cgca_m3pfem.f90:
           ! allocate( cgca_pfem_stress%stress( nels_pp, intp, comp ),
           cgca_pfem_stress%stress( iel, i, : ) = sigma

        END DO intpts
     end do elmnts
!*** end ParaFEM part ************************************************72

!*** CGPACK part *****************************************************72
! debug: dump stresses to stdout
!call cgca_pfem_sdmp
! dump stresses from last image for element 1
if ( cgca_img .eq. cgca_nimgs ) then
  do i = 1, nip
    write (*,"(2(a,i0),a,6es10.2)") "img: ", cgca_nimgs,               &
          " FE 1 int. p. ", i, " stress: ",                            &
          cgca_pfem_stress % stress( 1 , i , : )
  end do
end if

! all images sync here
sync all

! Calculate the mean stress tensor per image
!   subroutine cgca_pfem_simg( simg )
!   real( kind=rdef ), intent(out) :: simg(3,3)
call cgca_pfem_simg( cgca_stress )
write (*,*) "img:", cgca_img, " mean s tensor:", cgca_stress

! all images wait for each other, to make sure the stress arrays
! are not modified until all images calculate their mean values
sync all

! no real time increments in this problem
! I use the inverse of the length scale,
! which gives 1mm of crack propagation per increment maximum.
! I then can multiply it by a factor, e.g. a factor of 3 will mean
! that I get 3mm max ( 1/3 of the model ) per load increment.
cgca_time_inc = 2 * 1.0_rdef / cgca_length

! run cleavage for a correct number of iterations, which is a function
! of the characteristic length and the time increment
cgca_clvg_iter = nint( cgca_length * cgca_lres * cgca_time_inc )
if ( cgca_img .eq. 1 ) write (*,*) "load inc:", cgca_liter,            &
                                   "clvg iter:", cgca_clvg_iter

! Cray can use a version with CO_SUM.
! subroutine cgca_clvgp( coarray, rt, t, scrit, sub, gcus, &
!   periodicbc, iter, heartbeat, debug )
! sync all inside
! lower the crit stresses by a factor of 100.
call cgca_clvgp( cgca_space, cgca_grt, cgca_stress,                    &
     0.01_rdef * cgca_scrit, cgca_clvgsd, cgca_gcupdn, .false.,        &
     cgca_clvg_iter, 10, cgca_yesdebug )

!if ( cgca_img .eq. 1 ) write (*,*) "dumping model to file"
write ( cgca_citer, "(i0)" ) cgca_liter
!call cgca_pswci( cgca_space, cgca_state_type_frac,                     &
!                 "zf"//trim( cgca_citer )//".raw" )
!if ( cgca_img .eq. 1 ) write (*,*) "finished dumping model to file"
     
! Calculate number (volume) of fractured cells on each image.
! cgca_fracvol is a local, non-coarray, array, so no sync needed.
call cgca_fv( cgca_space, cgca_fracvol )

write (*,*) "img:", cgca_img, "fracvol:", cgca_fracvol

! calculate integrity, remote write
call cgca_pfem_intcalc1( cgca_c, cgca_fracvol )

! dump min integrity for all FE stored on this image
write (*,*) "img:", cgca_img, "min. integrity:",                            &
  minval( cgca_pfem_integrity % i )

! wait for integrity on all images to be calculated
sync all

! Young's modulus need to be updated on each image, local arrays only.
! The original Young's modulus value must be given as input.
! For each FE stored on this image.
call cgca_pfem_uym( e, nels_pp )

! dump updated Young's modulus
write (*,*) "img:", cgca_img, "*min* Young's mod:",                    &
            minval( cgca_pfem_enew )

sync all

!*** end CGPACK part *************************************************72


!*** ParaFEM part ****************************************************72

! end loading iterations
end do load_iter

! deallocate all arrays, moved from inside the loop
DEALLOCATE( p_pp )
deallocate( r_pp )
deallocate( x_pp )
deallocate( u_pp )
deallocate( d_pp )
deallocate( diag_precon_pp )
deallocate( storkm_pp )
deallocate( pmul_pp )
DEALLOCATE( xnew_pp )
DEALLOCATE( g_coord_pp )

!------------------------ write out displacements ----------------------

  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)

  IF ( numpe==1 ) THEN
     write(ch,'(I6.6)') numpe
     open( 12, file=argv(1:nlen)//".ensi.DISPL-"//ch,status='replace',    &
          action='write')
     write(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
     write(12,'(A/A/A)') "part", "     1","coordinates"
  END IF

  ALLOCATE( disp_pp(nodes_pp*ndim), source=zero )
  allocate( temp(nodes_pp), source=zero )
  CALL scatter_nodes( npes, nn, nels_pp, g_num_pp, nod, ndim, nodes_pp,  &
       node_start, node_end, eld_pp, disp_pp, 1 )
  DO i=1,ndim
     temp=zero
     DO j=1,nodes_pp
        k = i+(ndim*(j-1))
        temp(j) = disp_pp(k)
     END DO
     CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp)
  END DO

  IF ( numpe==1 ) then
     write( 11, '(A,F10.4)') "This analysis took: ", elap_time()-timest(1)
     close( 11 )
     close( 12 )
  end if
!*** end ParaFEM part ************************************************72


!*** CGPACK part *****************************************************72
! deallocate all CGPACK arrays
call cgca_ds(  cgca_space )
call cgca_drt( cgca_grt )
call cgca_pfem_sdalloc
call cgca_pfem_edalloc
call cgca_pfem_integdalloc
!*** end CGPACK part *************************************************72


!*** ParaFEM part ****************************************************72
CALL SHUTDOWN()
END PROGRAM xx14noio
