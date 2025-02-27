*     File ep_eXee/inclep_eXee1.h : Sat Mar 18 19:45:45 2000
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
      parameter (nep_eXeegraph =16)
      parameter (nep_eXeeextn  =6)


* number of all helicity states
      parameter (lep_eXeeag    =64)

      common /amprck/kmngr,  kmnext, kmlag

      parameter (nep_eXeegrpsq = nep_eXeegraph*nep_eXeegraph)
      common /aep_eXeemslct/jselg(nep_eXeegraph),jgraph,jgluon,jhiggs

* Color string information
      common /aep_eXeemcsti/ kmcbas, kmcbmx, icinfo(6), icolst
      common /aep_eXeemgrph/agcwrk(0:lep_eXeeag-1),agc(0:lep_eXeeag-1),
     &              aprop,ancp(nep_eXeegraph),ansp(0:nep_eXeegraph)
     &             ,cfmtx
      common /aep_eXeemgrpi/igraph(nep_eXeegraph)
      complex*16 agc, agcwrk, aprop

* Momenta of external particles
      common /aep_eXeemextr/peep_eXee1(4),peep_eXee2(4),peep_eXee3(4),
     &               peep_eXee4(4),peep_eXee5(4),peep_eXee6(4),
     &               prod(nep_eXeeextn, nep_eXeeextn)

      common /aep_eXeemextp/ptep_eXee2v,exep_eXee2v,cfep_eXee2v,
     &               ptep_eXee4m,exep_eXee4m,cfep_eXee4m,ptep_eXee9x,
     &               exep_eXee9x,cfep_eXee9x,ptep_eXee17m,exep_eXee17m,
     &               cfep_eXee17m,ptep_eXee33m,exep_eXee33m,
     &               cfep_eXee33m,ptep_eXee62m,exep_eXee62m,
     &               cfep_eXee62m

      real*8     ptep_eXee2v(4,2), exep_eXee2v(1)
      complex*16 cfep_eXee2v(2,2)
      real*8     ptep_eXee4m(4,2), exep_eXee4m(1)
      complex*16 cfep_eXee4m(2,2)
      real*8     ptep_eXee9x(4,2), exep_eXee9x(1)
      complex*16 cfep_eXee9x(2,2)
      real*8     ptep_eXee17m(4,2), exep_eXee17m(1)
      complex*16 cfep_eXee17m(2,2)
      real*8     ptep_eXee33m(4,2), exep_eXee33m(1)
      complex*16 cfep_eXee33m(2,2)
      real*8     ptep_eXee62m(4,2), exep_eXee62m(1)
      complex*16 cfep_eXee62m(2,2)

* Normalization
      common /aep_eXeemdbgg/fknorm,fkcall,nkcall

* Calculated table of amplitudes
      common /aep_eXeematbl/av, lt, indexg
      complex*16 av(0:lep_eXeeag-1)
      integer    lt(0:nep_eXeeextn), indexg(nep_eXeeextn)

* Spin average
      common/aep_eXeemspin/aspin,aident,jhs(nep_eXeeextn),
     &               jhe(nep_eXeeextn),jcpol(nep_eXeeextn)
