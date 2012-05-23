c  *********************************************************************
c  *                                                                   *
c  *                         subroutine chksd                          *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, DalTech, Mar  3, 2000
c  Latest Update: Mar  3, 2000
c
c  PURPOSE  checks that the seeds used across a variety of runs are all
c           different
c
c  DESCRIPTION
c  This routine looks through a sequence of nb generator seeds and checks
c  to make sure that they are all different. If one seed matches another,
c  the two base filenames are written to standard output along with a
c  warning.
c
c  ARGUMENTS
c
c   kseed	integer vector of length at least nb containing the sequence
c		of seeds to check. This sequence is ordered from smallest
c		to largest in this routine. (input/output)
c
c     ifl	character string vector of length at least nb containing the
c		names of the input files corresponding to the above seeds.
c		The filenames have equal seeds are printed to standard output
c		along with a warning. (input)
c
c     ind	integer vector of length at least nb used as an indexing
c		array. On output, ind(1) is the index of the smallest
c		element of kseed, ind(2) is the index of the next smallest
c		element of kseed, and so on. Thus
c
c		kseed(ind(1)) < kseed(ind(2)) < ... < kseed(ind(nb))
c
c		(output)
c
c      nb	integer containing the number of seeds provided. (input)
c
c  REVISION HISTORY:
c
c-------------------------------------------------------------------------
      subroutine chksd(kseed,ifl,ind,nb)
      integer kseed(*), ind(*)
      character*(*) ifl(*)

   1  format(5a)
c						first sort kseed
      call isorti(kseed,ind,nb)
c						now check to see if any equal
      do 10 i = 2, nb
         i1 = ind(i-1)
         i2 = ind(i)
         if( kseed(i1) .eq. kseed(i2) ) then
            l1 = lnblnk(ifl(i1))
            l2 = lnblnk(ifl(i2))  
            write(6,1)'Warning: generator seeds in the following two fil
     >es are equal!'
            write(6,1)'     1) ',ifl(i1)(1:l1)
            write(6,1)'     2) ',ifl(i2)(1:l2)
         endif
  10  continue

      return
      end
