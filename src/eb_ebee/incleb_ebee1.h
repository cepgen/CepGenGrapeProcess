*     File eb_ebee/incleb_ebee1.h : Sat Mar 18 19:45:05 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*

      include 'inclc.inc'
      include 'inclg.inc'

      parameter (loutgo =  2, lincom =  1)
      parameter (lantip = -1, lprtcl =  1)
      parameter (lscalr =  1)
      parameter (lepexa =  2, lepexv =  3)
      parameter (lepina =  4, lepinv =  4)
      parameter (lextrn =  2, lintrn =  4)

* table of amplitudes
      parameter (neb_ebeegraph =50)
      parameter (neb_ebeeextn  =6)


* number of color base
      parameter (ncbase =1)

* number of all helicity states
      parameter (leb_ebeeag    =64)

      common /amprck/kmngr,  kmnext, kmlag
      common /amprcs/smpref, smproc
      character*80   smpref, smproc

      parameter (neb_ebeegrpsq = neb_ebeegraph*neb_ebeegraph)
      common /aeb_ebeemslct/jselg(neb_ebeegraph),jgraph,jgluon,jhiggs

* Color string information
      common /aeb_ebeemcsti/ kmcbas, kmcbmx, kmcstr(2), icinfo(6), icolst
      common /aeb_ebeemgrph/agcwrk(0:leb_ebeeag-1),agc(0:leb_ebeeag-1),
     &              aprop,ancp(neb_ebeegraph),ansp(0:neb_ebeegraph)
     &             ,cfmtx
      common /aeb_ebeemgrpi/igraph(neb_ebeegraph)
      complex*16 agc, agcwrk, aprop

* Momenta of external particles
      common /aeb_ebeemextr/peeb_ebee1(4),peeb_ebee2(4),peeb_ebee3(4),
     &               peeb_ebee4(4),peeb_ebee5(4),peeb_ebee6(4),
     &               prod(neb_ebeeextn, neb_ebeeextn)

      common /aeb_ebeemextp/pteb_ebee2u,exeb_ebee2u,cfeb_ebee2u,
     &               pteb_ebee4m,exeb_ebee4m,cfeb_ebee4m,pteb_ebee9u,
     &               exeb_ebee9u,cfeb_ebee9u,pteb_ebee17m,exeb_ebee17m,
     &               cfeb_ebee17m,pteb_ebee33m,exeb_ebee33m,
     &               cfeb_ebee33m,pteb_ebee62m,exeb_ebee62m,
     &               cfeb_ebee62m

      real*8     pteb_ebee2u(4,2), exeb_ebee2u(1)
      complex*16 cfeb_ebee2u(2,2)
      real*8     pteb_ebee4m(4,2), exeb_ebee4m(1)
      complex*16 cfeb_ebee4m(2,2)
      real*8     pteb_ebee9u(4,2), exeb_ebee9u(1)
      complex*16 cfeb_ebee9u(2,2)
      real*8     pteb_ebee17m(4,2), exeb_ebee17m(1)
      complex*16 cfeb_ebee17m(2,2)
      real*8     pteb_ebee33m(4,2), exeb_ebee33m(1)
      complex*16 cfeb_ebee33m(2,2)
      real*8     pteb_ebee62m(4,2), exeb_ebee62m(1)
      complex*16 cfeb_ebee62m(2,2)

* Normalization
      common /aeb_ebeemdbgg/fknorm,fkcall,nkcall

* Calculated table of amplitudes
      common /aeb_ebeematbl/av, lt, indexg
      complex*16 av(0:leb_ebeeag-1)
      integer    lt(0:neb_ebeeextn), indexg(neb_ebeeextn)

* Spin average
      common/aeb_ebeemspin/aspin,aident,jhs(neb_ebeeextn),
     &               jhe(neb_ebeeextn),jcpol(neb_ebeeextn)