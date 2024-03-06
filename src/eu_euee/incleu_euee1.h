*     File eu_euee/incleu_euee1.h : Sat Mar 18 19:44:59 2000
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
      parameter (neu_eueegraph =48)
      parameter (neu_eueeextn  =6)


* number of color base
      parameter (ncbase =1)

* number of all helicity states
      parameter (leu_eueeag    =64)

      common /amprck/kmngr,  kmnext, kmlag
      common /amprcs/smpref, smproc
      character*80   smpref, smproc

      parameter (neu_eueegrpsq = neu_eueegraph*neu_eueegraph)
      common /aeu_eueemslct/jselg(neu_eueegraph),jgraph,jgluon,jhiggs

* Color string information
      common /aeu_eueemcsti/ kmcbas, kmcbmx, kmcstr(2), icinfo(6), icolst
      common /aeu_eueemgrph/agcwrk(0:leu_eueeag-1),agc(0:leu_eueeag-1),
     &              aprop,ancp(neu_eueegraph),ansp(0:neu_eueegraph)
     &             ,cfmtx
      common /aeu_eueemgrpi/igraph(neu_eueegraph)
      complex*16 agc, agcwrk, aprop

* Momenta of external particles
      common /aeu_eueemextr/peeu_euee1(4),peeu_euee2(4),peeu_euee3(4),
     &               peeu_euee4(4),peeu_euee5(4),peeu_euee6(4),
     &               prod(neu_eueeextn, neu_eueeextn)

      common /aeu_eueemextp/pteu_euee2p,exeu_euee2p,cfeu_euee2p,
     &               pteu_euee4m,exeu_euee4m,cfeu_euee4m,pteu_euee9p,
     &               exeu_euee9p,cfeu_euee9p,pteu_euee17m,exeu_euee17m,
     &               cfeu_euee17m,pteu_euee33m,exeu_euee33m,
     &               cfeu_euee33m,pteu_euee62m,exeu_euee62m,
     &               cfeu_euee62m

      real*8     pteu_euee2p(4,2), exeu_euee2p(1)
      complex*16 cfeu_euee2p(2,2)
      real*8     pteu_euee4m(4,2), exeu_euee4m(1)
      complex*16 cfeu_euee4m(2,2)
      real*8     pteu_euee9p(4,2), exeu_euee9p(1)
      complex*16 cfeu_euee9p(2,2)
      real*8     pteu_euee17m(4,2), exeu_euee17m(1)
      complex*16 cfeu_euee17m(2,2)
      real*8     pteu_euee33m(4,2), exeu_euee33m(1)
      complex*16 cfeu_euee33m(2,2)
      real*8     pteu_euee62m(4,2), exeu_euee62m(1)
      complex*16 cfeu_euee62m(2,2)

* Normalization
      common /aeu_eueemdbgg/fknorm,fkcall,nkcall

* Calculated table of amplitudes
      common /aeu_eueematbl/av, lt, indexg
      complex*16 av(0:leu_eueeag-1)
      integer    lt(0:neu_eueeextn), indexg(neu_eueeextn)

* Spin average
      common/aeu_eueemspin/aspin,aident,jhs(neu_eueeextn),
     &               jhe(neu_eueeextn),jcpol(neu_eueeextn)
