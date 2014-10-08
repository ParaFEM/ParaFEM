PROGRAM pf2ensibin_geo
!------------------------------------------------------------------------------
!      Program pf2ensibin_geo  Tool to produce binary ensight gold formatted  
!                              finite element mesh input file. Under 
!                              development, currently for 4-node tetrahedra
!------------------------------------------------------------------------------
  
!  USE precision  ; USE global_variables ; USE mp_interface
!  USE input      ; USE output           ; USE loading
!  USE timing     ; USE maths            ; USE gather_scatter
!  USE partition  ; USE elements         ; USE steering
!  USE geometry   ; USE pcg              ; USE new_library
  
  USE, INTRINSIC :: ISO_C_BINDING ! to output C binary file
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------
  
  INTEGER,PARAMETER   :: iwp=SELECTED_REAL_KIND(15)
  INTEGER, PARAMETER  :: ndim=3,nodof=1,nprops=5
  INTEGER             :: meshgen,partitioner,np_types,nels,nn,nr,nip,nod
  INTEGER             :: loaded_nodes,fixed_freedoms,nstep,npri
  INTEGER             :: limit,el_print,i_o
  INTEGER             :: i,j,k,prnwidth
  REAL(iwp),PARAMETER :: zero = 0.0_iwp
  REAL(iwp)           :: val0,dtim,theta,tol,etype,real_time
  CHARACTER(LEN=50)   :: fname,job_name
  CHARACTER(LEN=15)   :: element
  CHARACTER(LEN=80)   :: cbuffer
  
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------
  
  REAL(iwp),ALLOCATABLE :: temp_real(:),prop(:,:),val_f(:),val(:,:)
  INTEGER,ALLOCATABLE   :: temp_int(:),node(:),sense(:)
  
!------------------------------------------------------------------------------
! 3. Read job_name from the command line.
!    Read control data
!------------------------------------------------------------------------------
  
  CALL GETARG(1,job_name)
  
  PRINT *, "Starting conversion of ParaFEM input files"

  fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
  OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
  READ(10,*) element,meshgen,partitioner,np_types,                            &
             nels,nn,nr,nip,nod,loaded_nodes,fixed_freedoms,                  &
             val0,                                                            &
             dtim,nstep,npri,theta,                                           &
             tol,limit,el_print,i_o
  CLOSE(10)
  
!  fname   = job_name(1:INDEX(job_name, " ")-1)//".test"
!  OPEN(11,FILE=fname)
!  WRITE(11,*) element,meshgen,partitioner,np_types,                            &
!             nels,nn,nr,nip,nod,loaded_nodes,fixed_freedoms,                  &
!             val0,                                                            &
!             dtim,nstep,npri,theta,                                           &
!             tol,limit,el_print,i_o
!  CLOSE(11)
  
!----------------------------------------------------------------------------
! 4. Read material properties
!----------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ")-1) // ".mat"
  OPEN(11,FILE=fname,STATUS='OLD',ACTION='READ')
  ALLOCATE(prop(nprops+1,np_types))
  prop = zero
  !-Skip header  
  READ(11,*)
  READ(11,*)
  DO i=1,np_types
    READ(11,*) prop(:,i)
  END DO
  
!----------------------------------------------------------------------------
! 5. Open files
!----------------------------------------------------------------------------

  fname   = job_name(1:INDEX(job_name, " ")-1)//".d"
  OPEN(12,FILE=fname,STATUS='OLD',ACTION='READ')
  
  fname   = job_name(1:INDEX(job_name, " ")-1)//".bin.ensi.geo"
  OPEN(13,FILE=fname,STATUS="REPLACE",FORM="UNFORMATTED",                   &
                ACTION="WRITE", ACCESS="STREAM")
  
  OPEN(14,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.MATID',         &
          STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")
  
  OPEN(15,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_1_kx',    &
          STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

  OPEN(16,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_2_ky',    &
          STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

  OPEN(17,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_3_kz',    &
          STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

  OPEN(18,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_4_rho',    &
          STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

  OPEN(19,FILE=job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_5_cp',    &
          STATUS="REPLACE",FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")
  
!----------------------------------------------------------------------------
! 6. Read/Write geometry, MATID and material property files
!----------------------------------------------------------------------------
  
  !-Write .geo header
  cbuffer = "C Binary"                     ; WRITE(13) cbuffer
  cbuffer = "Problem name: "//job_name(1:INDEX(job_name, " ")-1) ; WRITE(13) cbuffer
  cbuffer = "Geometry files"               ; WRITE(13) cbuffer
  cbuffer = "node id off"                  ; WRITE(13) cbuffer
  cbuffer = "element id off"               ; WRITE(13) cbuffer
  cbuffer = "part"                         ; WRITE(13) cbuffer
  WRITE(13) int(1,kind=c_int)
  cbuffer = "Volume"                       ; WRITE(13) cbuffer
  cbuffer = "coordinates"                  ; WRITE(13) cbuffer
  WRITE(13) int(nn,kind=c_int)
  
  !- Write file headers
  DO i=14,19
    cbuffer = "Alya Ensight Gold --- Scalar per-element variable file"
    WRITE(i) cbuffer
    cbuffer = "part"  ; WRITE(i) cbuffer
    WRITE(i) int(1,kind=c_int)
    cbuffer = "tetra4"  ; WRITE(i) cbuffer
  END DO

!  !-Write .MATID header
!  cbuffer = "Alya Ensight Gold --- Scalar per-element variable file"
!  WRITE(14) cbuffer
!  cbuffer = "part"  ; WRITE(14) cbuffer
!  WRITE(14) int(1,kind=c_int)
!  cbuffer = "tetra4"  ; WRITE(14) cbuffer
!  
!  !-Write .ELMAT headers
!  cbuffer = "Alya Ensight Gold --- Scalar per-element variable file"
!  WRITE(15) cbuffer
!  cbuffer = "part"  ; WRITE(15) cbuffer
!  WRITE(15) int(1,kind=c_int)
!  cbuffer = "tetra4"  ; WRITE(15) cbuffer

  PRINT *, "Writing nodal coordinates"
  !-Read/Write nodal coordinates  
  ALLOCATE(temp_real(ndim+1))
  temp_real = zero
  DO j=1,ndim
  !-Skip section header  
    READ(12,*)
    READ(12,*)
    DO i=1,nn
      READ(12,*) temp_real
      WRITE(13) real(temp_real(j+1),kind=c_float)
    END DO
    IF(j/=ndim) REWIND(12)
  END DO
  DEALLOCATE(temp_real)
  !-Write element subheader in .geo
  cbuffer = "tetra4" ; WRITE(13) cbuffer
  WRITE(13) int(nels,kind=c_int)

  PRINT *, "Writing elements and material properties"
  ALLOCATE(temp_int(ndim+6))
  temp_int = 0
  !-Skip section header  
  READ(12,*)
  DO i = 1,nels
    READ(12,*) temp_int
    WRITE(13) int(temp_int(5),kind=c_int)
    WRITE(13) int(temp_int(6),kind=c_int)
    WRITE(13) int(temp_int(7),kind=c_int)
    WRITE(13) int(temp_int(8),kind=c_int)
    etype=temp_int(9)*1.0
    WRITE(14) real(etype,kind=c_float)
    WRITE(15) real(prop(2,temp_int(9)),kind=c_float)
    WRITE(16) real(prop(3,temp_int(9)),kind=c_float)
    WRITE(17) real(prop(4,temp_int(9)),kind=c_float)
    WRITE(18) real(prop(5,temp_int(9)),kind=c_float)
    WRITE(19) real(prop(6,temp_int(9)),kind=c_float)
  END DO
  DEALLOCATE(temp_int)
  
  CLOSE(12);CLOSE(13);CLOSE(14);CLOSE(15);CLOSE(16);CLOSE(17);CLOSE(18);CLOSE(19)
  
!----------------------------------------------------------------------------
! 7. Read/Write file containing fixed nodes
!----------------------------------------------------------------------------
  
  IF(fixed_freedoms > 0) THEN
    PRINT *, "Writing fixed nodes"
    ALLOCATE(node(fixed_freedoms),sense(fixed_freedoms),val_f(fixed_freedoms))
    fname = job_name(1:INDEX(job_name, " ")-1) // ".fix"
    OPEN(20,FILE=fname,STATUS='OLD',ACTION='READ')
    DO i = 1,fixed_freedoms
      READ(20,*)node(i),sense(i),val_f(i)
    END DO
    CLOSE(20)
  
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
  
!----------------------------------------------------------------------------
! 8. Read/Write file containing loaded nodes
!----------------------------------------------------------------------------
  
  IF(loaded_nodes > 0) THEN
    PRINT *, "Writing loaded nodes"
    ALLOCATE(node(loaded_nodes),val(nodof,loaded_nodes))
    fname = job_name(1:INDEX(job_name, " ")-1) // ".lds"
    OPEN(22, FILE=fname, STATUS='OLD', ACTION='READ')
    DO i = 1,loaded_nodes 
      READ(22,*) node(i),val(:,i)
    END DO
    CLOSE(22)

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
        j=j+1
      ELSE
        WRITE(23) real(zero,kind=c_float)
      END IF
    END DO
    
    DEALLOCATE(node)
    CLOSE(23)
  END IF
  
!------------------------------------------------------------------------------
! 9. Write case file
!------------------------------------------------------------------------------
  
  PRINT *, "Writing case file"
  fname   = job_name(1:INDEX(job_name, " ")-1)//".bin.ensi.case"
  OPEN(24,FILE=fname)

  WRITE(24,'(A/A)')    "#", "# Post-processing file generated by subroutine &
                             &WRITE_ENSI in "
  WRITE(24,'(A,A,/A)') "#"," Smith, Griffiths and Margetts, 'Programming the &
                             &Finite Element Method',","# Wiley, 2013."        
  WRITE(24,'(A/A/A)')  "#","# Ensight Gold Format","#"
  WRITE(24,'(2A/A)')   "# Problem name: ",job_name(1:INDEX(job_name, " ")-1),"#"
  WRITE(24,'(A/A/A)')  "FORMAT","type:  ensight gold","GEOMETRY"
  WRITE(24,'(2A/A)')   "model: 1  ",job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.geo',"VARIABLE"
  WRITE(24,'(2A)')     "scalar per element:  material      ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.MATID'
  WRITE(24,'(2A)')     "scalar per element:  material_kx   ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_1_kx'
  WRITE(24,'(2A)')     "scalar per element:  material_ky   ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_2_ky'
  WRITE(24,'(2A)')     "scalar per element:  material_kz   ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_3_kz'
  WRITE(24,'(2A)')     "scalar per element:  material_rho  ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_4_rho'
  WRITE(24,'(2A)')     "scalar per element:  material_cp   ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.ELMAT_5_cp'
  WRITE(24,'(2A)')     "scalar per node: 1   temperature   ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.NDTTR-******'
  WRITE(24,'(2A)')     "scalar per node:     fixed_nodes   ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.NDFIX'
  WRITE(24,'(2A)')     "scalar per node:     loaded_nodes  ",                &
                        job_name(1:INDEX(job_name, " ")-1)//'.bin.ensi.NDLDS'
  WRITE(24,'(A/A)')     "TIME","time set:     1"
  WRITE(24,'(A,I5)')    "number of steps:",nstep/npri+1
  WRITE(24,'(A,I5)')    "filename start number:",1
  WRITE(24,'(A,I5)')    "filename increment:",1
  WRITE(24,'(A)')       "time values:"
  
  prnwidth = 5
  j=1
  WRITE(24,'(E12.5)',ADVANCE='no') zero
  timesteps: DO i=1,nstep
    real_time = i*dtim
    IF(i/npri*npri==i)THEN
      WRITE(24,'(E12.5)',ADVANCE='no') real_time
      j=j+1
    END IF
    IF(j==prnwidth)THEN
        WRITE(24,*)''
        j=0
    END IF
  END DO timesteps
  WRITE(24,*)''
  CLOSE(24)
  
  PRINT *, "Conversion complete"
  
  CALL shutdown()
  
END PROGRAM pf2ensibin_geo
