*     File ebb_ebbee/inclebb_ebbee1.h : Sat Mar 18 19:45:05 2000
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
      parameter (nebb_ebbeegraph =50)
      parameter (nebb_ebbeeextn  =6)


* number of color base
      parameter (ncbase =1)

* number of all helicity states
      parameter (lebb_ebbeeag    =64)

      common /amprck/kmngr,  kmnext, kmlag

      parameter (nebb_ebbeegrpsq = nebb_ebbeegraph*nebb_ebbeegraph)
      common /aebb_ebbeemslct/jselg(nebb_ebbeegraph),jgraph,jgluon,
     &               jhiggs

* Color string information
      common /aebb_ebbeemcsti/ kmcbas, kmcbmx, kmcstr(2), icinfo(6), icolst
      common /aebb_ebbeemgrph/agcwrk(0:lebb_ebbeeag-1),agc(0:lebb_ebbeeag-1),
     &              aprop,ancp(nebb_ebbeegraph),ansp(0:nebb_ebbeegraph)
     &             ,cfmtx
      common /aebb_ebbeemgrpi/igraph(nebb_ebbeegraph)
      complex*16 agc, agcwrk, aprop

* Momenta of external particles
      common /aebb_ebbeemextr/peebb_ebbee1(4),peebb_ebbee2(4),
     &               peebb_ebbee3(4),peebb_ebbee4(4),peebb_ebbee5(4),
     &               peebb_ebbee6(4),
     &               prod(nebb_ebbeeextn, nebb_ebbeeextn)

      common /aebb_ebbeemextp/ptebb_ebbee2u,exebb_ebbee2u,cfebb_ebbee2u,
     &               ptebb_ebbee4m,exebb_ebbee4m,cfebb_ebbee4m,
     &               ptebb_ebbee9u,exebb_ebbee9u,cfebb_ebbee9u,
     &               ptebb_ebbee17m,exebb_ebbee17m,cfebb_ebbee17m,
     &               ptebb_ebbee33m,exebb_ebbee33m,cfebb_ebbee33m,
     &               ptebb_ebbee62m,exebb_ebbee62m,cfebb_ebbee62m

      real*8     ptebb_ebbee2u(4,2), exebb_ebbee2u(1)
      complex*16 cfebb_ebbee2u(2,2)
      real*8     ptebb_ebbee4m(4,2), exebb_ebbee4m(1)
      complex*16 cfebb_ebbee4m(2,2)
      real*8     ptebb_ebbee9u(4,2), exebb_ebbee9u(1)
      complex*16 cfebb_ebbee9u(2,2)
      real*8     ptebb_ebbee17m(4,2), exebb_ebbee17m(1)
      complex*16 cfebb_ebbee17m(2,2)
      real*8     ptebb_ebbee33m(4,2), exebb_ebbee33m(1)
      complex*16 cfebb_ebbee33m(2,2)
      real*8     ptebb_ebbee62m(4,2), exebb_ebbee62m(1)
      complex*16 cfebb_ebbee62m(2,2)

* Normalization
      common /aebb_ebbeemdbgg/fknorm,fkcall,nkcall

* Calculated table of amplitudes
      common /aebb_ebbeematbl/av, lt, indexg
      complex*16 av(0:lebb_ebbeeag-1)
      integer    lt(0:nebb_ebbeeextn), indexg(nebb_ebbeeextn)

* Spin average
      common/aebb_ebbeemspin/aspin,aident,jhs(nebb_ebbeeextn),
     &               jhe(nebb_ebbeeextn),jcpol(nebb_ebbeeextn)
