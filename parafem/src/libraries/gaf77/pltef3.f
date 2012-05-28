c  ********************************************************************
c  *                                                                  *
c  *                       subroutine pltfld                          *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, Dalhousie University, Feb 19, 2003
c  Latest Update: Feb 19, 2003
c
c  PURPOSE  to output a plane of the 3-D log-elastic field to a file in
c           DISPLAY format (for display on a PostScript printer).
c
c  This routine takes an array (3-D) of data and dumps a 2-D plane through
c  the data to a file having a format readable by DISPLAY. Arguments are
c  as follows;
c
c   istat	unit number connected to a file to which error messages and
c		warnings are to be written. (input)
c
c     job	character string containing the title of the run. (input)
c
c    sub1	character string containing the subtitle of the run. (input)
c
c    sub2	character string containing the sub-subtitle of the run.
c		(input)
c
c    efld	real array of size at least nxe x nye x nze containing the
c		data to be displayed graphically. It is assumed that the
c		leading dimensions of efld are (nxe,nye). (input)
c
c    nrfx	integer giving the leading dimension of the efld array.
c		(input)
c
c    nrfy	integer giving the second dimension of the efld array.
c		(input)
c
c     nxe	integer containing the x-dimension size of the field to
c		display. (input)
c
c     nye	integer containing the y-dimension size of the field to
c		display. (input)
c
c     nze	integer containing the z-dimension size of the field to
c		display. (input)
c
c      xd	real value containing the physical dimension of the data in
c		the X (1st index) direction. (input)
c
c      yd	real value containing the physical dimension of the data in
c		the Y (2nd index) direction. (input)
c
c      zd	real value containing the physical dimension of the data in
c		the Z (3rd index) direction. (input)
c
c   ieplt	integer vector of length at least 3 containing the element
c		indices corresponding to the plane over which the random field
c		plot is to be drawn. The elements of ieplt correspond
c		to the (x,y,z) element indices and zero values in this vector
c		correspond to all elements in that direction. Only one of the
c		components is expected to be non-zero. So, for example, the
c		vector (0,26,0) means that the x-z plane at y-element = 26 is
c		to be drawn. If a 50 x 50 x 30 element domain has two
c		footings, one centered at x-node = 26 and y-node = 15 and
c		the other centered at x-node = 26 and y-node = 35, then
c		normally the y-z plane through x-node = 26 would be drawn
c		to see the random field under the two footings. Thus, in this
c		case, idplt might be (25,0,0) or (26,0,0) (depending on which
c		side of the nodal plane you want to be on). The z-indices are
c		measured from the top surface of the soil mass, thus z = 1 is
c		the top layer of the soil mass and z = nze is the bottom
c		layer, etc. (input)
c
c   title	character string containing the title of the plot. (input)
c
c   iefld	unit number to which output is sent. (input)
c
c  NOTES:
c	1) in RSETL3D, the Z-direction is negative when pointing downwards
c	   into the soil mass while the z-index starts counting at the top
c	   of the soil mass and increases downwards. This has been reflected
c	   in the mapping of the random field to the finite element mesh, so
c	   planes here can be dumped in the usual order.
c-----------------------------------------------------------------------------
      subroutine pltef3(istat,job,sub1,sub2,efld,nrfx,nrfy,
     >                  nxe,nye,nze,xd,yd,zd,ieplt,
     >                  title,iefld)
      real*4 efld(nrfx,nrfy,*)
      integer ieplt(*)
      character*(*) job, sub1, sub2, title
      character*128 sub0

   1  format(3a)
   2  format(2e13.5)
   4  format(2i8)
   5  format(6e12.4)
c					string lengths
      lj = lnblnk(job)
      l1 = lnblnk(sub1)
      l2 = lnblnk(sub2)
      lt = lnblnk(title)
c					construct sub0
      if( lt .ge. 128 ) then
         sub0 = title(1:128)
      else
         sub0 = title(1:lt)
         if( ieplt(1) .gt. 0 ) then
            call prins1(sub0(lt+1:),', x-plane = %i',ieplt(1))
         elseif( ieplt(2) .gt. 0 ) then
            call prins1(sub0(lt+1:),', y-plane = %i',ieplt(2))
         elseif( ieplt(3) .gt. 0 ) then
            call prins1(sub0(lt+1:),', z-plane = %i',ieplt(3))
         endif
      endif
      lt = lnblnk(sub0)
c					write out common header
      write(iefld,1) job(1:lj)
      write(iefld,1) sub0(1:lt)
      write(iefld,1) sub1(1:l1)
      write(iefld,1) sub2(1:l2)
      write(iefld,1) '2'
c					plot y-z plane?
      if( ieplt(1) .gt. 0 ) then
         i = ieplt(1)
         if( i .gt. nxe ) then
            write(istat,1)'Warning: illegal plane index provided for ran
     >dom field plot.'
            write(istat,1)'         See plot output file for details.'
            call print1(iefld,
     >      '*** ERROR: illegal y-z plane x-index: ieplt(1) = %i%n',i)
            call print1(iefld,
     >      '           must be between 0 and %i.%n',nxe)
            close(iefld)
            return
         endif
         write(iefld,1) 'y'
         write(iefld,1) 'z'
         write(iefld,1) 'log-E'
         write(iefld,2) 0., yd
         write(iefld,2) -zd, 0.
         write(iefld,4) nye, nze
         write(iefld,5) ((efld(i,j,k), j = 1, nye), k = 1, nze)

c					plot x-z plane?
      elseif( ieplt(2) .gt. 0 ) then
         j = ieplt(2)
         if( j .gt. nye ) then
            write(istat,1)'Warning: illegal plane index provided for ran
     >dom field plot.'
            write(istat,1)'         See plot output file for details.'
            call print1(iefld,
     >      '*** ERROR: illegal x-z plane y-index: ieplt(2) = %i%n',j)
            call print1(iefld,
     >      '           must be between 0 and %i.%n',nye)
            close(iefld)
            return
         endif
         write(iefld,1) 'x'
         write(iefld,1) 'z'
         write(iefld,1) 'log-E'
         write(iefld,2) 0., xd
         write(iefld,2) -zd, 0.
         write(iefld,4) nxe, nze
         write(iefld,5) ((efld(i,j,k), i = 1, nxe), k = 1, nze)

c					plot x-y plane?
      elseif( ieplt(3) .gt. 0 ) then
         k = ieplt(3)
         if( k .gt. nze ) then
            write(istat,1)'Warning: illegal plane index provided for ran
     >dom field plot.'
            write(istat,1)'         See plot output file for details.'
            call print1(iefld,
     >      '*** ERROR: illegal x-y plane z-index: ieplt(3) = %i%n',k)
            call print1(iefld,
     >      '           must be between 0 and %i.%n',nze)
            close(iefld)
            return
         endif
         write(iefld,1) 'x'
         write(iefld,1) 'y'
         write(iefld,1) 'log-E'
         write(iefld,2) 0., xd
         write(iefld,2) 0., yd
         write(iefld,4) nxe, nye
         write(iefld,5) ((efld(i,j,k), i = 1, nxe), j = 1, nye)
      else
         write(istat,1)'Warning: zero plot plane specified for random fi
     >eld plot.'
         write(istat,1)'         Check data file. Plotting abandoned.'
      endif
      close(iefld)

      return
      end
