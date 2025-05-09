*     File eub_eubee/incleub_eubee1.h : Sat Mar 18 19:44:59 2000
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
      parameter (neub_eubeegraph =48)
      parameter (neub_eubeeextn  =6)


* number of color base
      parameter (ncbase =1)

* number of all helicity states
      parameter (leub_eubeeag    =64)

      common /amprck/kmngr,  kmnext, kmlag

      parameter (neub_eubeegrpsq = neub_eubeegraph*neub_eubeegraph)
      common /aeub_eubeemslct/jselg(neub_eubeegraph),jgraph,jgluon,
     &               jhiggs

* Color string information
      common /aeub_eubeemcsti/ kmcbas, kmcbmx, kmcstr(2), icinfo(6), icolst
      common /aeub_eubeemgrph/agcwrk(0:leub_eubeeag-1),agc(0:leub_eubeeag-1),
     &              aprop,ancp(neub_eubeegraph),ansp(0:neub_eubeegraph)
     &             ,cfmtx
      common /aeub_eubeemgrpi/igraph(neub_eubeegraph)
      complex*16 agc, agcwrk, aprop

* Momenta of external particles
      common /aeub_eubeemextr/peeub_eubee1(4),peeub_eubee2(4),
     &               peeub_eubee3(4),peeub_eubee4(4),peeub_eubee5(4),
     &               peeub_eubee6(4),
     &               prod(neub_eubeeextn, neub_eubeeextn)

      common /aeub_eubeemextp/pteub_eubee2p,exeub_eubee2p,cfeub_eubee2p,
     &               pteub_eubee4m,exeub_eubee4m,cfeub_eubee4m,
     &               pteub_eubee9p,exeub_eubee9p,cfeub_eubee9p,
     &               pteub_eubee17m,exeub_eubee17m,cfeub_eubee17m,
     &               pteub_eubee33m,exeub_eubee33m,cfeub_eubee33m,
     &               pteub_eubee62m,exeub_eubee62m,cfeub_eubee62m

      real*8     pteub_eubee2p(4,2), exeub_eubee2p(1)
      complex*16 cfeub_eubee2p(2,2)
      real*8     pteub_eubee4m(4,2), exeub_eubee4m(1)
      complex*16 cfeub_eubee4m(2,2)
      real*8     pteub_eubee9p(4,2), exeub_eubee9p(1)
      complex*16 cfeub_eubee9p(2,2)
      real*8     pteub_eubee17m(4,2), exeub_eubee17m(1)
      complex*16 cfeub_eubee17m(2,2)
      real*8     pteub_eubee33m(4,2), exeub_eubee33m(1)
      complex*16 cfeub_eubee33m(2,2)
      real*8     pteub_eubee62m(4,2), exeub_eubee62m(1)
      complex*16 cfeub_eubee62m(2,2)

* Normalization
      common /aeub_eubeemdbgg/fknorm,fkcall,nkcall

* Calculated table of amplitudes
      common /aeub_eubeematbl/av, lt, indexg
      complex*16 av(0:leub_eubeeag-1)
      integer    lt(0:neub_eubeeextn), indexg(neub_eubeeextn)

* Spin average
      common/aeub_eubeemspin/aspin,aident,jhs(neub_eubeeextn),
     &               jhe(neub_eubeeextn),jcpol(neub_eubeeextn)
