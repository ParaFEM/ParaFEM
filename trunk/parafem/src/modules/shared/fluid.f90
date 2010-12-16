MODULE FLUID

  !------------------------------------------------------------------------
  ! Contains subroutines for fluid flow problems
  !
  !      Authors  - I.M Smith, D.V. Griffiths, L. Margetts
  !      
  !      History
  !      20/02/09 - Created
  !-------------------------------------------------------------------------

  ! Routines in this module:
  !   FORMNF_CAVITY
  !   FIXITY_DATA
  !   FMKEPARP
  !   FMKEPARV
  !   FORMNF_UPVW
  !   FIXITY_UPVW
  !   FORMUPVW

  !---
  !--- Default modules
  !---
  
  USE precision
  USE geometry

  CONTAINS
  
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE FORMNF_CAVITY(nf,nftol,nodof,aa,bb,cc,nels,nxe,nze,problem)
  
    !/****f* structure_dof/formnf_cavity
    !*  NAME
    !*    SUBROUTINE: formnf_cavity
    !*  SYNOPSIS
    !*    Usage:      CALL formnf_cavity(nf,nftol,nodof,aa,bb,cc,nels,     &
    !*                                   nxe,nze,problem)
    !*  FUNCTION
    !*    Forms node freedom array NF for a rectangular lid driven cavity
    !*    of any size
    !*  
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    nels                : Integer
    !*                        : Total number of elements
    !* 
    !*    nxe                 : Integer
    !*                        : Number of elements in the x direction
    !* 
    !*    nye                 : Integer
    !*                        : Number of elements in the y direction
    !* 
    !*    problem             : Integer
    !*                        : Code for the problem type
    !*                          1 = full cavity
    !*                          2 = symmetric half cavity
    !*                          Any other value = undefined
    !*
    !*    rest(nr,nodof+1)    : Integer
    !*                        : List of the nodes with some degree of freedom
    !*                          fixed. In the input, the degrees of freedom
    !*                          fixed have 1 and the degrees of freedom not
    !*                          fixed have 0. In the output, the degrees of
    !*                          freedom fixed have 0 and the degrees of freedom
    !*                          not fixed have the global equation number
    !*
    !*    The following arguments have the INTENT(INOUT) attribute:    
    !*
    !*    nf(nn,nodof+1)      : Integer
    !*                        : Node freedom array
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    05.03.2008
    !*  COPYRIGHT
    !*    (c) University of Manchester
    !******
    !*  Need to modify with an error trap for problem.
    !*  Also need to deallocate allocated arrays.
    !*/


    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: nf(:,:)
    INTEGER, INTENT(IN)    :: nels,nxe,nze,problem
    INTEGER                :: i,j,n,nye,nn,nodof,iel
    INTEGER, ALLOCATABLE   :: num(:)
    REAL(iwp), INTENT(IN)  :: aa,bb,cc,nftol
    REAL(iwp)              :: minx,maxx,miny,maxy,minz,maxz
    REAL(iwp), ALLOCATABLE :: coord(:,:),g_coord(:,:)

    !-----------------------------------------------------------------------
    ! 1. Allocate arrays 
    !-----------------------------------------------------------------------
   
    nn = UBOUND(nf,2)
    
    ALLOCATE(coord(20,3),num(20),g_coord(3,nn))

    !-----------------------------------------------------------------------
    ! 2. Initiallize 
    !-----------------------------------------------------------------------
   
    minx = 0._iwp
    miny = 0._iwp
    minz = 0._iwp
    nye  = (nels/nxe)/nze
    maxx = aa * nxe
    maxy = bb * nye
    maxz = cc * nze
    
    DO i = 1,nn
      nf(1,i)   = 0
      nf(2:4,i) = 1
    END DO

    !-----------------------------------------------------------------------
    ! 3. Find coordinates of nodes and mark equations
    !-----------------------------------------------------------------------
     
    DO iel = 1,nels
        CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,num)
        DO i = 1,7,2
          nf(1,num(i)) = 1
        END DO
        DO i = 13,19,2
          nf(1,num(i)) = 1
        END DO
        g_coord(:,num) = transpose(coord)
    END DO
    
    
    node_freedoms: SELECT CASE(problem)
    
    !-----------------------------------------------------------------------
    ! 4. Find node freedoms for a full lid-driven cavity
    !-----------------------------------------------------------------------
    
    CASE (1)
    
    DO i = 1,nn
      IF(ABS(g_coord(1,i)-minx)<nftol.OR.abs(g_coord(1,i)-maxx)<nftol) THEN
        nf(2:4,i) = 0
      END IF
      IF(ABS(g_coord(2,i)-miny)<nftol.OR.abs(g_coord(2,i)-maxy)<nftol) THEN
        nf(2:4,i) = 0
      END IF
      IF(ABS(g_coord(3,i)-minz)<nftol) THEN
        nf(3:4,i) = 0
        nf(2,i)   = 1
      END IF
      IF(ABS(g_coord(3,i)+maxz)<nftol) THEN
        nf(2:4,i) = 0
        IF (ABS(g_coord(2,i)-miny)<nftol) THEN
          IF (ABS(g_coord(1,i)-minx)<nftol) nf(1,i) = 0
        END IF
      END IF
    END DO
    
    !-----------------------------------------------------------------------
    ! 5. Find node freedoms for a symmetric half lid-driven cavity
    !-----------------------------------------------------------------------
    
    CASE (2)
    
    DO i = 1,nn
      IF(ABS(g_coord(1,i)-minx)<nftol.OR.abs(g_coord(1,i)-maxx)<nftol) THEN
        nf(2:4,i) = 0
      END IF
      IF(ABS(g_coord(2,i)-miny)<nftol) THEN
        nf(2:4,i) = 0
      END IF
      IF(ABS(g_coord(2,i)-maxy)<nftol) nf(3,i) = 0
      IF(abs(g_coord(3,i)-minz)<nftol) THEN
        nf(2,i) = 1
        nf(3,i) = 0
        nf(4,i) = 0
      END IF
      IF(ABS(g_coord(3,i)+maxz)<nftol) THEN
        nf(2:4,i) = 0
        IF (ABS(g_coord(2,i)-miny)<nftol) THEN
          IF (ABS(g_coord(1,i)-minx)<nftol) nf(1,i) = 0
        END IF
      END IF
    END DO
    
    END SELECT node_freedoms

    !-----------------------------------------------------------------------
    ! 6. Number the equations
    !-----------------------------------------------------------------------
    
    n = 0
    DO j = 1,nn
      DO i = 1,nodof
        IF(nf(i,j)/=0) THEN
          n       = n+1
          nf(i,j) = n
        END IF
      END DO
    END DO  
   
  RETURN

  END SUBROUTINE FORMNF_CAVITY

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE FIXITY_DATA(g_coord,node,sense,val,lid_velocity,nftol,nn,problem)
  
    !/****f* structure_dof/fixity_data
    !*  NAME
    !*    SUBROUTINE: fixity_data
    !*  SYNOPSIS
    !*    Usage:      CALL fixity_data(g_coord,node,sense,val,lid_velocity,   &
    !*                                 nftol,nn,problem)
    !*  FUNCTION
    !*    Fixes the value of velocity on all nodes of the lid
    !*  
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    problem             : Integer
    !*                        : Code for the problem type.
    !*                          1 = full cavity
    !*                          2 = symmetric half cavity
    !*                          Any other value = undefined    
    !* 
    !*    g_coord(:,:)        : Real
    !*                        : Co-ordinates of the nodes
    !* 
    !*    lid_velocity        : Real
    !*                        : The prescribed velocity of the lid
    !* 
    !*    nftol               : Real
    !*                        : The tolerance used when searching for nodes
    !*                          on the basis of their position
    !*
    !*  OUTPUTS
    !*    The following arguments have the INTENT(OUT) attribute:
    !* 
    !*    node(:)             : Integer
    !*                        : A vector listing the nodes positioned on the
    !*                          lid that will have an imposed velocity
    !* 
    !*    sense(:)            : Integer
    !*                        : A vector identifying which of the x, y or z 
    !*                          freedoms will have an imposed velocity
    !*                          1 = x direction
    !*                          2 = y direction
    !*                          3 = z direction
    !*                          Any other value is currently undefined
    !* 
    !*    val(:)              : Real
    !*                        : A vector that holds the velocity values
    !* 
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    19.04.2002
    !*  COPYRIGHT
    !*    (c) University of Manchester
    !******
    !*  Need to modify with an error trap for problem.
    !*  Need to have error trap for sense too.
    !*/

    IMPLICIT NONE
    
    INTEGER,  INTENT(OUT) :: node(:),sense(:)
    INTEGER,  INTENT(IN)  :: problem
    INTEGER               :: i,l,nn
    REAL(iwp),INTENT(IN)  :: g_coord(:,:),lid_velocity,nftol
    REAL(iwp),INTENT(OUT) :: val(:)
    REAL(iwp)             :: minz,miny

    !-----------------------------------------------------------------------
    ! 1. Initiallize
    !-----------------------------------------------------------------------

    nn    = UBOUND(g_coord,2)
    l     = 0
    node  = 0
    sense = 0
    val   = 0._iwp
    miny  = 0._iwp
    minz  = 0._iwp

    lid_vel: SELECT CASE (problem)
    
    !-----------------------------------------------------------------------
    ! 2. Fix the fluid velocities for a full lid-driven cavity
    !-----------------------------------------------------------------------
 
    CASE (1)
  
    DO i = 1,nn
      IF(ABS(g_coord(3,i)-minz)<nftol) THEN
        l        = l+1
        node(l)  = i
        sense(l) = 2
        val(l)   = lid_velocity
      END IF
    END DO
    
    !-----------------------------------------------------------------------
    ! 3. Fix the fluid velocities for a half lid-driven cavity
    !-----------------------------------------------------------------------

    CASE (2)
  
    DO i=1,nn
      IF(ABS(g_coord(3,i)-minz)<nftol) THEN
        l=l+1
        node(l)=i
        sense(l)=2
        val(l)=lid_velocity
      END IF
    END DO

    END SELECT lid_vel

    RETURN

  END SUBROUTINE fixity_data
    
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE FMKEPARP(keparp,c12,c32,c42)

    !/****f* fluid_flow/fmkeparp
    !*  NAME
    !*    SUBROUTINE: fmkeparp
    !*  SYNOPSIS
    !*    Usage:      CALL fmkeparp(keparp,c12,c32,c42)
    !*  FUNCTION
    !*    Forms sub-parent 'stiffness' matrix related to the pressure component
    !*    in 3-d for u-p-v-w version of Navier Stokes. KE incomplete without 
    !*    KEPARV and C11. c.f. FMKEPARV and ST_C11 (stores C11 by element)
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:    
    !*
    !*    c12(:,:)            : Real
    !*                        : The c12 submatrix from the 3D Navier-Stokes
    !*                          element stiffness matrix
    !*
    !*    c32(:,:)            : Real
    !*                        : The c32 submatrix from the 3D Navier-Stokes
    !*                          element stiffness matrix
    !*
    !*    c42(:,:)            : Real
    !*                        : The c42 submatrix from the 3D Navier-Stokes
    !*                          element stiffness matrix
    !*    
    !*  OUTPUTS
    !*    The following arguments have the INTENT(OUT) attribute:    
    !*
    !*    keparp(:,:)         : Real
    !*                        : The sub-parent 'stiffness' matrix
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    24.04.2002
    !*  COPYRIGHT
    !*    (c) University of Manchester
    !******
    !*/    

    IMPLICIT NONE
    
    REAL(iwp),INTENT(IN)  :: c12(:,:),c32(:,:),c42(:,:)
    REAL(iwp),INTENT(OUT) :: keparp(:,:)
    INTEGER               :: nod,nodf
    
    nodf = 8
    nod  = 20
    
    keparp(1:nod,1:nodf)                 = c12 
    keparp(nod+1:nod+nod,1:nodf)         = c32
    keparp(nod+nod+1:nod+nod+nod,1:nodf) = c42
    
    RETURN
    
  END SUBROUTINE FMKEPARP

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
 

  SUBROUTINE FMKEPARV(keparv,c21,c23,c24)

    !/****f* fluid_flow/fmkeparv
    !*  NAME
    !*    SUBROUTINE: fmkeparv
    !*  SYNOPSIS
    !*    Usage:      CALL fmkeparv(keparv,c21,c23,c24)
    !*  FUNCTION
    !*    Forms sub-parent 'stiffness' matrix related to the velocity component
    !*    in 3-d for u-p-v-w version of Navier Stokes. KE incomplete without 
    !*    KEPARP and C11. c.f. FMKEPARP and ST_C11 (stores C11 by element)
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:    
    !*
    !*    c21(:,:)            : Real
    !*                        : The c21 submatrix from the 3D Navier-Stokes
    !*                          element stiffness matrix
    !*
    !*    c23(:,:)            : Real
    !*                        : The c23 submatrix from the 3D Navier-Stokes
    !*                          element stiffness matrix
    !*
    !*    c24(:,:)            : Real
    !*                        : The c24 submatrix from the 3D Navier-Stokes
    !*                          element stiffness matrix
    !*    
    !*  OUTPUTS
    !*    The following arguments have the INTENT(OUT) attribute:    
    !*
    !*    keparv(:,:)         : Real
    !*                        : The sub-parent 'stiffness' matrix
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    24.04.2002
    !*  COPYRIGHT
    !*    (c) University of Manchester
    !******
    !*/    

     IMPLICIT NONE
     
     REAL(iwp),INTENT(IN)  :: c21(:,:),c23(:,:),c24(:,:)
     REAL(iwp),INTENT(OUT) :: keparv(:,:)
     INTEGER               :: nod,nodf
     
     nodf =  8
     nod  = 20
     
     keparv(1:nodf,1:nod)                 = c21 
     keparv(1:nodf,nod+1:nod+nod)         = c23
     keparv(1:nodf,nod+nod+1:nod+nod+nod) = c24
  
     RETURN
     
   END SUBROUTINE FMKEPARV
 
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine formnf_upvw(nf,nftol,nodof,aa,bb,cc,nels,nxe,nze,problem)
!   forms nf array for lid-driven cavity problem of any size
!   rectangular cavity only upvw order
 implicit none
 integer, intent(inout)::nf(:,:)
 integer, intent(in)::nels,nxe,nze,problem
 integer::i,j,n,nye,nn,nodof,iel
 integer, allocatable::num(:)
 real(iwp), intent(in)::aa,bb,cc,nftol
 real(iwp)::minx,maxx,miny,maxy,minz,maxz
 real(iwp), allocatable :: coord(:,:),g_coord(:,:)
 nn=ubound(nf,2)
 allocate (coord(20,3),num(20),g_coord(3,nn))
 minx=0.; miny=0.; minz=0.; nye=(nels/nxe)/nze
 maxx=aa*nxe; maxy=bb*nye; maxz=cc*nze
 do i=1,nn
     nf(1:4,i)=1;nf(2,i)=0
 end do
 do iel=1,nels
     call geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,num)
     do i=1,7,2; nf(2,num(i))=1; end do
     do i=13,19,2; nf(2,num(i))=1; end do
     g_coord(:,num)=transpose(coord)
 end do
 node_freedoms: select case(problem)
 !---- lid-driven cavity ------------------------------------------------------
 case (1)
  do i=1,nn
   if(abs(g_coord(1,i)-minx)<nftol.or.abs(g_coord(1,i)-maxx)<nftol) then
        nf(1,i)=0; nf(3,i)=0; nf(4,i)=0
   end if
   if(abs(g_coord(2,i)-miny)<nftol.or.abs(g_coord(2,i)-maxy)<nftol) then
     nf(1,i)=0; nf(3,i)=0; nf(4,i)=0
   end if
   if(abs(g_coord(3,i)-minz)<nftol) then
     nf(1,i)=1; nf(3,i)=0; nf(4,i)=0
   end if
   if(abs(g_coord(3,i)+maxz)<nftol) then
     nf(1,i)=0; nf(3,i)=0; nf(4,i)=0
     if (abs(g_coord(2,i)-miny)<nftol) then
     if (abs(g_coord(1,i)-minx)<nftol) nf(2,i)=0
     end if
   end if
 end do
 !---- lid-driven 'half' cavity -----------------------------------------------
 case (2)
  do i=1,nn
   if(abs(g_coord(1,i)-minx)<nftol.or.abs(g_coord(1,i)-maxx)<nftol) then
        nf(1,i)=0; nf(3,i)=0; nf(4,i)=0
   end if
   if(abs(g_coord(2,i)-miny)<nftol) then
     nf(1,i)=0; nf(3,i)=0; nf(4,i)=0
   end if
   if(abs(g_coord(2,i)-maxy)<nftol) nf(3,i)=0
   if(abs(g_coord(3,i)-minz)<nftol) then
     nf(1,i)=1; nf(3,i)=0; nf(4,i)=0
   end if
   if(abs(g_coord(3,i)+maxz)<nftol) then
     nf(1,i)=0; nf(3,i)=0; nf(4,i)=0
     if (abs(g_coord(2,i)-miny)<nftol) then
     if (abs(g_coord(1,i)-minx)<nftol) nf(2,i)=0
     end if
   end if
 end do
  end select node_freedoms
 n=0
   do j=1,nn; do i=1,nodof
     if(nf(i,j)/=0) then
     n=n+1;nf(i,j)=n
   end if
   end do; end do  
return
end subroutine formnf_upvw

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine fixity_upvw(g_coord,node,sense,val,lid_velocity,nftol,nn,problem)
!fixes value of lid velocity on all nodes of lid - upvw order
implicit none
integer,intent(out)::node(:),sense(:)
integer::i,l,nn,problem
real(iwp),intent(in)::g_coord(:,:),lid_velocity,nftol
real(iwp),intent(out)::val(:)
real(iwp)::minz,miny
nn=ubound(g_coord,2)
l=0;node=0;sense=0;val=0.;miny=0.;minz=.0
lid_vel: select case(problem)
 !---- lid-driven cavity ------------------------------------------------------
 case (1)
  do i=1,nn
    if(abs(g_coord(3,i)-minz)<nftol) then
       l=l+1; node(l)=i; sense(l)=1; val(l)=lid_velocity
    end if
  end do
!---- lid-driven 'half' cavity -----------------------------------------------
 case (2)
  do i=1,nn
    if(abs(g_coord(2,i)-miny)<nftol) then
       l=l+1; node(l)=i; sense(l)=1; val(l)=lid_velocity
    end if
  end do
! do l=1,nn; write(11,'(2i5,e12.4)') node(l),sense(l),val(l); end do
end select lid_vel
return
end subroutine fixity_upvw

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

! Modified From Smith and Griffiths:

subroutine formupvw(storke,iel,c11,c12,c21,c23,c32,c24,c42)
!   forms unsymmetrical 'stiffness' matrix in 3-d for
!   u-p-v-w version of navier stokes
 implicit none
 real(iwp),intent(in)::c11(:,:),c12(:,:),c21(:,:),c23(:,:),c32(:,:),     &
                  c24(:,:),c42(:,:)
 real(iwp),intent(out):: storke(:,:,:); integer:: nod,nodf,ntot,iel
 nod=ubound(c11,1); nodf=ubound(c21,1); ntot=nod+nodf+nod+nod
 storke(1:nod,1:nod,iel) = c11; storke(1:nod,nod+1:nod+nodf,iel) = c12 
 storke(nod+1:nod+nodf,1:nod,iel) = c21
 storke(nod+1:nod+nodf,nod+nodf+1:nod+nodf+nod,iel) = c23
 storke(nod+1:nod+nodf,nod+nodf+nod+1:ntot,iel) = c24 
 storke(nod+nodf+1:nod+nodf+nod,nod+1:nod+nodf,iel) = c32
 storke(nod+nodf+1:nod+nodf+nod,nod+nodf+1:nod+nodf+nod,iel) = c11 
 storke(nod+nodf+nod+1:ntot,nod+1:nod+nodf,iel) = c42
 storke(nod+nodf+nod+1:ntot,nod+nodf+nod+1:ntot,iel) = c11     
end subroutine formupvw
 
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

END MODULE FLUID
