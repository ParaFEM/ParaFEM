c  *********************************************************************
c  *                                                                   *
c  *                         subroutine ofile                          *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, Dalhousie University, Jul  3, 2006
c  Latest Update: Jul  3, 2006
c
c  PURPOSE  creates an output file name having the same name as the provided
c           input file name, but with a new extension
c
c  DESCRIPTION
c  This routine takes an input file name, infile, finds its last `.' character
c  and replaces the extension with the new extension provided. If the new and
c  old extensions are identical, the output file name will be the same as the
c  input file name.
c
c  If the input file name does not have an extension, then the new extension
c  is added to the end of the input name to create the output name.
c
c  ARGUMENTS
c
c  infile	character string containing the input file name. (input)
c
c  outfil	character string which, on output, will contain the filename
c		with the new extension. (output)
c
c     ext	character string containing the extension which is to replace
c		the input file's extension to create the outfil name. If the
c		input file does not have an extension, then the extension is
c		added to the end of the infile name. NOTE: the extension
c		string is assumed to include the leading '.'. (input)
c
c  REVISION HISTORY:
c
c-------------------------------------------------------------------------
      subroutine ofile(infile,outfil,ext)
      character*(*) infile, outfil, ext

   1  format(3a)

      il   = lnblnk( infile )
      do 10 i = il, 2, -1
         if( infile(i:i) .eq. '.' ) then
            ib = i - 1
            go to 20
         endif
  10  continue
      ib = il

  20  lo = len(outfil)
      le = lnblnk(ext)
      if( ib+le .gt. lo ) then
         write(6,1)'Error: output filename length greater than output fi
     >lename string.'
         outfil = 'out'
         outfil(4:) = ext
         lo = lnblnk(outfil)
         write(6,1)'       output filename being set to ', outfil(1:lo)
         return
      endif

      outfil = infile
      outfil(ib+1:) = ext

      return
      end
