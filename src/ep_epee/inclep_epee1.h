*     File ep_epee/inclep_epee1.h : Sat Mar 18 19:44:57 2000
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
      parameter (nep_epeegraph =16)
      parameter (nep_epeeextn  =6)


* number of all helicity states
      parameter (lep_epeeag    =64)

      common /amprck/kmngr,  kmnext, kmlag

      parameter (nep_epeegrpsq = nep_epeegraph*nep_epeegraph)
      common /aep_epeemslct/jselg(nep_epeegraph),jgraph,jgluon,jhiggs

* Color string information
      common /aep_epeemcsti/ kmcbas, kmcbmx, icinfo(6), icolst
      common /aep_epeemgrph/agcwrk(0:lep_epeeag-1),agc(0:lep_epeeag-1),
     &              aprop,ancp(nep_epeegraph),ansp(0:nep_epeegraph)
     &             ,cfmtx
      common /aep_epeemgrpi/igraph(nep_epeegraph)
      complex*16 agc, agcwrk, aprop

* Momenta of external particles
      common /aep_epeemextr/peep_epee1(4),peep_epee2(4),peep_epee3(4),
     &               peep_epee4(4),peep_epee5(4),peep_epee6(4),
     &               prod(nep_epeeextn, nep_epeeextn)

      common /aep_epeemextp/ptep_epee2v,exep_epee2v,cfep_epee2v,
     &               ptep_epee4m,exep_epee4m,cfep_epee4m,ptep_epee9v,
     &               exep_epee9v,cfep_epee9v,ptep_epee17m,exep_epee17m,
     &               cfep_epee17m,ptep_epee33m,exep_epee33m,
     &               cfep_epee33m,ptep_epee62m,exep_epee62m,
     &               cfep_epee62m

      real*8     ptep_epee2v(4,2), exep_epee2v(1)
      complex*16 cfep_epee2v(2,2)
      real*8     ptep_epee4m(4,2), exep_epee4m(1)
      complex*16 cfep_epee4m(2,2)
      real*8     ptep_epee9v(4,2), exep_epee9v(1)
      complex*16 cfep_epee9v(2,2)
      real*8     ptep_epee17m(4,2), exep_epee17m(1)
      complex*16 cfep_epee17m(2,2)
      real*8     ptep_epee33m(4,2), exep_epee33m(1)
      complex*16 cfep_epee33m(2,2)
      real*8     ptep_epee62m(4,2), exep_epee62m(1)
      complex*16 cfep_epee62m(2,2)

* Normalization
      common /aep_epeemdbgg/fknorm,fkcall,nkcall

* Calculated table of amplitudes
      common /aep_epeematbl/av, lt, indexg
      complex*16 av(0:lep_epeeag-1)
      integer    lt(0:nep_epeeextn), indexg(nep_epeeextn)

* Spin average
      common/aep_epeemspin/aspin,aident,jhs(nep_epeeextn),
     &               jhe(nep_epeeextn),jcpol(nep_epeeextn)
