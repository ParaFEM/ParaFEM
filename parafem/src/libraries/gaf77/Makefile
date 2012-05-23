# Creates a library of double and single precision F77 math, simulation,
# and I/O routines callable from other f77 programs. See README, Summary.doc,
# and the modules themselves for details. To include, ensure that $LIBDIR
# is in your loader's search path and add the option "-lGAF77" to the f77
# command line.
# Written by Gordon A. Fenton, Princeton University, 1988.
#
# Latest Update: Jul 30, 2006
#
# REVISION HISTORY:
# 3.58	inserted contents of block data `dgblck.f' directly into dphi, dphinv,
#	and dgmcdf.f; dgblck.f is no longer needed. (May 18/97)
# 3.59  made fft routines consistent with real*4 versions. (Jul 2/97)
# 3.60	added Makefile for Lahey Fortan for PC's (Aug 18/97)
# 3.61	corrected arguments in calls to dfft1d in dpspct and dxspct (Dec 1/98)
# 3.62	improved accuracy of dbetnv, added dfdist, dfdinv, dchdst (Dec 4/98)
# 3.63	corrected error flag in dinvsm (Mar 30/99)
# 3.64	added dadfit, dstatm, and dsort, now use *'s for dummy
#	dimensions (Jun 9/99)
# 3.98	added etime_.c, fdate_.c, and getpid_.c for Linux emulation of BSD
#	(Apr 12/97)
# 3.99	inserted contents of block data `gblck.f' directly into phi and phinv;
#	gblck.f is no longer needed. (May 18/97)
# 4.00	revised printv module, adding prftm?.f and getfsp.f modules,
#	moved Linux BSD modules (see above) to separate library. (Jun 9/97)
# 4.01	corrected error in fft3d.f (Jul 2/97)
# 4.02	corrected error rendering powers greater than 9 in prfmte.f (Jul 28/97)
# 4.03	added Makefile for Lahey Fortan for PC's (Aug 18/97)
# 4.04	added fdist, fdinv, and chdst (Dec 4/98)
# 4.05	added dadfit and dsort (Jun 9/99)
# 4.33	moved Linux port BSD modules to separate library, concatenated to GAF77
#	(Jun 10/97)
# 4.34	added analysis routines, dcor1d, dsdf1d, dvfn1d, vregrs, and xregrs
#	(Jul 9/97)
# 4.35	provided additional histogram routines hstbkt and hstprd (Aug 12/97)
# 4.36	added Makefile for Lahey Fortan for PC's (Aug 18/97)
# 4.37	allowed multi-digit subscript in hstprd.f (Aug 24/97)
# 4.38	added 1-D Fractal-Scale variance function (dlfsc1.f) (Nov 19/97)
# 4.39	added analysis routine dcov1d (Dec 11/97)
# 4.40	corrected fractal-scale model dlfsc1, deleted Lahey Makefile (Jan 8/98)
# 4.41	corrected endpoints in dvfn1d (Mar 6/98)
# 4.42	now use *'s for dummy dimensions (Jun 9/99)
# 5.00	combined GAF77, DGAF77, and GAFsim into one library (Jun 9/99)
# 5.01	getarg is now an intrinsic function (Jul 27/99)
# 5.02	added etime, fdate, getlog, getpid, for those compilers which lack
#	these (ie, f77 under RedHat 5.2 and lower) (Jul 30/99)
# 5.03	replaced \\'s with char(92) for portability (Dec 1/99)
# 5.04	revised prints to add blanks to end of string (Jan 25/00)
# 5.05	added lognormal dist to adfit and dadfit (Feb 28/00)
# 5.06	added isort, dsorti, sorti, isorti; various sorters (Mar 3/00)
# 5.07	corrected bug in prins? routines when f format turns into e (Mar 4/00)
# 5.08	improved accuracy of adfit and statms (Mar 7/00)
# 5.10	modified LAS routines, now using Gauss quadrature (Jun 15/00)
# 5.11	corrected errors in adfit (Aug 15/00)
# 5.12	minor corrections to chifit and hstprd to reflect `pvalue' (Aug 27/00)
# 5.13	modified hstbkt.f to allow non-normalization (Dec 8/00)
# 5.14	modified LAS routines, use Gauss quadrature for all but overlapping
#	domains (Apr 9/01)
# 5.15	modified LAS routines dcvaa?.f and dcvab?.f - now use Gauss quadrature
#	to compute all local averages (May 9/01)
# 5.20	modified print? routines - eliminated '\' escapes and moved w.d
#	formats to between % and 'f', 'e', or 'i' (as in C). This is to
#	make libGAF77 compatible with Win??'s Visual Fortran. (Oct 4/01)
#	NOTE: under Winxx, LDOS in prints.f must be changed to .true.
# 5.21	added print0 and prins0 routines (Sep 8/02)
# 5.22	use 16 pt quadrature in dcvaa3 and dcvab3 (Mar 21/03)
# 5.23	use 16 pt quadrature in dcvaa1, dcvab1, dcvaa2, dcvab2, dlavx2,
#	and dlavx3 for variance/covariance calcs (Apr 8/03)
# 5.24	corrected smoothing algorithm in dsdf1d (Apr 28/03)
# 5.25	corrected trailing %space implementation in prints/v; use 20pt
#	Gauss quadrature and correct simplified cases in dlavx3, dcvab3,
#	and dcvaa3. (Jul 2/03)
# 5.26	forward FFT used in dpspct.f and dsdf1d.f, which no longer divides
#	by N, so SDF estimates were incorrect by a factor of N*N. (Apr 5/04)
# 5.27	renamed scans?.f files to nscan?.f and made them functions returning
#	number of items read. (Jul 13/04)
# 5.28	corrected string pointer js in nscanv.f. (Aug 20/04)
# 5.29	corrected limited sig digits in nrdfp and nrddp (Nov 26/04)
# 5.30	made X-Y the isotropic plane in dlafs3.f (Jan 11/05)
# 5.31	corrected computation of Beta dist in hstprd.f (Feb 16/05)
# 5.32	added srtgrp.f and chksd.f from old sim/rutils (Mar 17/05)
# 5.33	increased MXK from 16 to 256 in las1g. (Aug 10/05)
# 5.34	corrected parameter n for idist 6 in histp and hstprd (Nov 3/05)
# 5.35	added julian date function, corrected fft documentation (Apr 26/06)
# 5.36	added gdate, the inverse of julian (May 5/06)
# 5.37	implemented init = -1 in las1g (Jun 3/06)
# 5.38	added dbest1 and dbest2, corrected dpspct and dsdf1d (Jun 19/06)
# 5.39	added ofile, to create output file names, added test for zeroes to
#	statms and adfit/dadfit. (Jul 30/06)
#----------------------------------------------------------------------------

Version = 5.39

FOR     = df
FFLAGS  = /optimize:5
INSTALL = move
DELETE  = del
#				LIBNAME is the name of the library
LIBNAME	= GAF77
#				LIBDIR is where the library is to be installed
LIBDIR	= d:\gordon\f77\lib

FILES =	AtxB.f    AxB.f     AxBt.f    addot.f   adfit.f   adot.f    argfn.f   \
	argfp.f   argin.f   argop.f   aswap.f   averg.f   axpy.f    best1.f   \
	best2.f   beta.f    betnv.f   bsect.f   cabs.f    cgelu.f   chdst.f   \
	chifit.f  chinv.f   chksd.f   chol2.f   chol2h.f  choll.f   cholu.f   \
	clear.f   cmplx.f   condn.f   corn2d.f  corn3d.f  cov1d.f   cpplu.f   \
	dAtxB.f   dAxB.f    dAxBt.f   dadfit.f  daverg.f  daxpy.f   dbeta.f   \
	dbetnv.f  dbsect.f  dcabs.f   dcgelu.f  dchdst.f  dchinv.f  dchol2.f  \
	dcholl.f  dcholu.f  dclear.f  dcmplx.f  dcondn.f  dcor1d.f  dcov1d.f  \
	dcpplu.f  dcvaa1.f  dcvaa2.f  dcvaa3.f  dcvab1.f  dcvab2.f  dcvab3.f  \
	dcvit1.f  dcvit2.f  dcvit3.f  dcvmt2.f  dcvmt3.f  ddot.f    deigen.f  \
	deignv.f  dexpnd.f  dfacld.f  dfactr.f  dfdinv.f  dfdist.f  dfft1d.f  \
	dfft2d.f  dfft3d.f  dflpos.f  dfxpnt.f  dgamln.f  dgamma.f  dgausv.f  \
	dgjinv.f  dgmcdf.f  dgmdst.f  dgminv.f  dgselu.f  dgsint.f  dimprv.f  \
	dinslv.f  dinvsm.f  dlace1.f  dlafr1.f  dlafr2.f  dlafs3.f  dlavx1.f  \
	dlavx2.f  dlavx3.f  dlfsc1.f  dllinv.f  dlnvmp.f  dlsep2.f  dlsep3.f  \
	dlsfr2.f  dlsfr3.f  dlsmp1.f  dlspx1.f  dlspx2.f  dlspx3.f  dluinv.f  \
	dmaxnt.f  dminmx.f  dmodnr.f  dmodsc.f  dnewtr.f  dnmprv.f  dnsolv.f  \
	dntddm.f  dnvgss.f  dodmax.f  dorthH.f  dot3d.f   dot3db.f  dpder1.f  \
	dpder2.f  dphi.f    dphinv.f  dpnorm.f  dpolyd.f  dpolyn.f  dprtme.f  \
	dprtmt.f  dpsolv.f  dpspct.f  dqrgde.f  dqrgrs.f  dqrgsq.f  dqrhde.f  \
	dqrhim.f  dqrhsq.f  dqtran.f  drdfpt.f  dround.f  dsdf1d.f  dsecnt.f  \
	dsifa.f   dsisl.f   dslvch.f  dslvld.f  dsmall.f  dsmslv.f  dsort.f   \
	dsorti.f  dspgss.f  dspint.f  dspln3.f  dspplu.f  dsqnml.f  dstatm.f  \
	dstudi.f  dstudt.f  dsvdec.f  dswap.f   dtoesv.f  dtpint.f  dtridg.f  \
	dtrilu.f  dtrisv.f  dvdots.f  dvfn1d.f  dvgath.f  dvgaus.f  dvnorm.f  \
	dxspct.f  dzmatx.f  edge3d.f  eigen.f   eignv.f   expnd.f   facld.f   \
	factr.f   fdinv.f   fdist.f   fexps1.f  fexps2.f  fexps3.f  fft1d.f   \
	fft1g.f   fft2d.f   fft2g.f   fft3d.f   fft3g.f   flpos.f   formc1.f  \
	fosps1.f  fxpnt.f   gamln.f   gamma.f   gausv.f   gblck.f   gencmb.f  \
	getfsp.f  getop.f   gjinv.f   gmcdf.f   gmdst.f   gminv.f   grcdf.f   \
	gselu.f   gsint.f   histp.f   hstbkt.f  hstprd.f  iamax.f   ichrd.f   \
	idamax.f  idchrd.f  imprv.f   inslv.f   intr2d.f  intr3d.f  invsm.f   \
	iseed.f   isort.f   isorti.f  izero.f   las1g.f   las1i.f   las2g.f   \
	las2i.f   las3g.f   las3i.f   lcomp.f   lczer.f   ldczer.f  ldigit.f  \
	llinv.f   lnblnk.f  lnvmp.f   luinv.f   maxnt.f   minmx.f   modnr.f   \
	modsc.f   nchoos.f  nctarg.f  ndigit.f  newtr.f   nmprv.f   nrdarg.f  \
	nrddp.f   nrdfp.f   nrdint.f  nrdlog.f  nscan1.f  nscan2.f  nscan3.f  \
	nscan4.f  nscan5.f  nscan6.f  nscanv.f  nsolv.f   ntddm.f   ntwom.f   \
	nvgss.f   odmax.f   openfl.f  orthH.f   pder1.f   pder2.f   phi.f     \
	phinv.f   plan3d.f  pnorm.f   polydd.f  polyn.f   prfmte.f  prfmtf.f  \
	prfmti.f  prins0.f  prins1.f  prins2.f  prins3.f  prins4.f  prins5.f  \
	prins6.f  prins7.f  prins8.f  print0.f  print1.f  print2.f  print3.f  \
	print4.f  print5.f  print6.f  print7.f  print8.f  prints.f  printv.f  \
	prtme.f   prtmt.f   psfmte.f  psfmtf.f  psfmti.f  psolv.f   pspct.f   \
	qrgde.f   qrgrs.f   qrgsq.f   qrhde.f   qrhim.f   qrhsq.f   qtran.f   \
	randf.f   randu.f   round.f   secnt.f   second.f  side2d.f  side3d.f  \
	sifa.f    sisl.f    slvch.f   slvld.f   small.f   smslv.f   sort.f    \
	sorti.f   spgss.f   spint.f   spln3.f   spplu.f   sqnml.f   srtgrp.f  \
	statms.f  studi.f   studt.f   svdec.f   tbm2g.f   tbm3g.f   texp2d.f  \
	texp3d.f  thin1d.f  toesv.f   tpint.f   tridg.f   trilu.f   trisv.f   \
	uuinv.f   vcorel.f  vgath.f   vgaus.f   vnorm.f   vregrs.f  xregrs.f  \
	xspct.f   zmatx.f   derf.f    julian.f  dbest1.f  dbest2.f  ofile.f


OBJS  =	AtxB.obj    AxB.obj     AxBt.obj    addot.obj   adfit.obj   adot.obj    argfn.obj   \
	argfp.obj   argin.obj   argop.obj   aswap.obj   averg.obj   axpy.obj    best1.obj   \
	best2.obj   beta.obj    betnv.obj   bsect.obj   cabs.obj    cgelu.obj   chdst.obj   \
	chifit.obj  chinv.obj   chksd.obj   chol2.obj   chol2h.obj  choll.obj   cholu.obj   \
	clear.obj   cmplx.obj   condn.obj   corn2d.obj  corn3d.obj  cov1d.obj   cpplu.obj   \
	dAtxB.obj   dAxB.obj    dAxBt.obj   dadfit.obj  daverg.obj  daxpy.obj   dbeta.obj   \
	dbetnv.obj  dbsect.obj  dcabs.obj   dcgelu.obj  dchdst.obj  dchinv.obj  dchol2.obj  \
	dcholl.obj  dcholu.obj  dclear.obj  dcmplx.obj  dcondn.obj  dcor1d.obj  dcov1d.obj  \
	dcpplu.obj  dcvaa1.obj  dcvaa2.obj  dcvaa3.obj  dcvab1.obj  dcvab2.obj  dcvab3.obj  \
	dcvit1.obj  dcvit2.obj  dcvit3.obj  dcvmt2.obj  dcvmt3.obj  ddot.obj    deigen.obj  \
	deignv.obj  dexpnd.obj  dfacld.obj  dfactr.obj  dfdinv.obj  dfdist.obj  dfft1d.obj  \
	dfft2d.obj  dfft3d.obj  dflpos.obj  dfxpnt.obj  dgamln.obj  dgamma.obj  dgausv.obj  \
	dgjinv.obj  dgmcdf.obj  dgmdst.obj  dgminv.obj  dgselu.obj  dgsint.obj  dimprv.obj  \
	dinslv.obj  dinvsm.obj  dlace1.obj  dlafr1.obj  dlafr2.obj  dlafs3.obj  dlavx1.obj  \
	dlavx2.obj  dlavx3.obj  dlfsc1.obj  dllinv.obj  dlnvmp.obj  dlsep2.obj  dlsep3.obj  \
	dlsfr2.obj  dlsfr3.obj  dlsmp1.obj  dlspx1.obj  dlspx2.obj  dlspx3.obj  dluinv.obj  \
	dmaxnt.obj  dminmx.obj  dmodnr.obj  dmodsc.obj  dnewtr.obj  dnmprv.obj  dnsolv.obj  \
	dntddm.obj  dnvgss.obj  dodmax.obj  dorthH.obj  dot3d.obj   dot3db.obj  dpder1.obj  \
	dpder2.obj  dphi.obj    dphinv.obj  dpnorm.obj  dpolyd.obj  dpolyn.obj  dprtme.obj  \
	dprtmt.obj  dpsolv.obj  dpspct.obj  dqrgde.obj  dqrgrs.obj  dqrgsq.obj  dqrhde.obj  \
	dqrhim.obj  dqrhsq.obj  dqtran.obj  drdfpt.obj  dround.obj  dsdf1d.obj  dsecnt.obj  \
	dsifa.obj   dsisl.obj   dslvch.obj  dslvld.obj  dsmall.obj  dsmslv.obj  dsort.obj   \
	dsorti.obj  dspgss.obj  dspint.obj  dspln3.obj  dspplu.obj  dsqnml.obj  dstatm.obj  \
	dstudi.obj  dstudt.obj  dsvdec.obj  dswap.obj   dtoesv.obj  dtpint.obj  dtridg.obj  \
	dtrilu.obj  dtrisv.obj  dvdots.obj  dvfn1d.obj  dvgath.obj  dvgaus.obj  dvnorm.obj  \
	dxspct.obj  dzmatx.obj  edge3d.obj  eigen.obj   eignv.obj   expnd.obj   facld.obj   \
	factr.obj   fdinv.obj   fdist.obj   fexps1.obj  fexps2.obj  fexps3.obj  fft1d.obj   \
	fft1g.obj   fft2d.obj   fft2g.obj   fft3d.obj   fft3g.obj   flpos.obj   formc1.obj  \
	fosps1.obj  fxpnt.obj   gamln.obj   gamma.obj   gausv.obj   gblck.obj   gencmb.obj  \
	getfsp.obj  getop.obj   gjinv.obj   gmcdf.obj   gmdst.obj   gminv.obj   grcdf.obj   \
	gselu.obj   gsint.obj   histp.obj   hstbkt.obj  hstprd.obj  iamax.obj   ichrd.obj   \
	idamax.obj  idchrd.obj  imprv.obj   inslv.obj   intr2d.obj  intr3d.obj  invsm.obj   \
	iseed.obj   isort.obj   isorti.obj  izero.obj   las1g.obj   las1i.obj   las2g.obj   \
	las2i.obj   las3g.obj   las3i.obj   lcomp.obj   lczer.obj   ldczer.obj  ldigit.obj  \
	llinv.obj   lnblnk.obj  lnvmp.obj   luinv.obj   maxnt.obj   minmx.obj   modnr.obj   \
	modsc.obj   nchoos.obj  nctarg.obj  ndigit.obj  newtr.obj   nmprv.obj   nrdarg.obj  \
	nrddp.obj   nrdfp.obj   nrdint.obj  nrdlog.obj  nscan1.obj  nscan2.obj  nscan3.obj  \
	nscan4.obj  nscan5.obj  nscan6.obj  nscanv.obj  nsolv.obj   ntddm.obj   ntwom.obj   \
	nvgss.obj   odmax.obj   openfl.obj  orthH.obj   pder1.obj   pder2.obj   phi.obj     \
	phinv.obj   plan3d.obj  pnorm.obj   polydd.obj  polyn.obj   prfmte.obj  prfmtf.obj  \
	prfmti.obj  prins0.obj  prins1.obj  prins2.obj  prins3.obj  prins4.obj  prins5.obj  \
	prins6.obj  prins7.obj  prins8.obj  print0.obj  print1.obj  print2.obj  print3.obj  \
	print4.obj  print5.obj  print6.obj  print7.obj  print8.obj  prints.obj  printv.obj  \
	prtme.obj   prtmt.obj   psfmte.obj  psfmtf.obj  psfmti.obj  psolv.obj   pspct.obj   \
	qrgde.obj   qrgrs.obj   qrgsq.obj   qrhde.obj   qrhim.obj   qrhsq.obj   qtran.obj   \
	randf.obj   randu.obj   round.obj   secnt.obj   second.obj  side2d.obj  side3d.obj  \
	sifa.obj    sisl.obj    slvch.obj   slvld.obj   small.obj   smslv.obj   sort.obj    \
	sorti.obj   spgss.obj   spint.obj   spln3.obj   spplu.obj   sqnml.obj   srtgrp.obj  \
	statms.obj  studi.obj   studt.obj   svdec.obj   tbm2g.obj   tbm3g.obj   texp2d.obj  \
	texp3d.obj  thin1d.obj  toesv.obj   tpint.obj   tridg.obj   trilu.obj   trisv.obj   \
	uuinv.obj   vcorel.obj  vgath.obj   vgaus.obj   vnorm.obj   vregrs.obj  xregrs.obj  \
	xspct.obj   zmatx.obj   derf.obj    julian.obj  dbest1.obj  dbest2.obj  ofile.obj


all: lib install

lib: $(OBJS)
	lib /out:$(LIBNAME).lib $(OBJS)

install: lib
	$(DELETE) $(LIBDIR)\$(LIBNAME).lib
	$(INSTALL) $(LIBNAME).lib $(LIBDIR)

clean:
	$(DELETE) *.obj
