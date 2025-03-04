*     File eub_eubtt/incleub_eubtt1.h : Sat Mar 18 19:45:00 2000
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
      parameter (neub_eubttgraph =25)
      parameter (neub_eubttextn  =6)


* number of color base
      parameter (ncbase =1)

* number of all helicity states
      parameter (leub_eubttag    =64)

      common /amprck/kmngr,  kmnext, kmlag

      parameter (neub_eubttgrpsq = neub_eubttgraph*neub_eubttgraph)
      common /aeub_eubttmslct/jselg(neub_eubttgraph),jgraph,jgluon,
     &               jhiggs

* Color string information
      common /aeub_eubttmcsti/ kmcbas, kmcbmx, kmcstr(2), icinfo(6), icolst
      common /aeub_eubttmgrph/agcwrk(0:leub_eubttag-1),agc(0:leub_eubttag-1),
     &              aprop,ancp(neub_eubttgraph),ansp(0:neub_eubttgraph)
     &             ,cfmtx
      common /aeub_eubttmgrpi/igraph(neub_eubttgraph)
      complex*16 agc, agcwrk, aprop

* Momenta of external particles
      common /aeub_eubttmextr/peeub_eubtt1(4),peeub_eubtt2(4),
     &               peeub_eubtt3(4),peeub_eubtt4(4),peeub_eubtt5(4),
     &               peeub_eubtt6(4),
     &               prod(neub_eubttextn, neub_eubttextn)

      common /aeub_eubttmextp/pteub_eubtt2p,exeub_eubtt2p,cfeub_eubtt2p,
     &               pteub_eubtt4m,exeub_eubtt4m,cfeub_eubtt4m,
     &               pteub_eubtt9p,exeub_eubtt9p,cfeub_eubtt9p,
     &               pteub_eubtt17m,exeub_eubtt17m,cfeub_eubtt17m,
     &               pteub_eubtt33o,exeub_eubtt33o,cfeub_eubtt33o,
     &               pteub_eubtt62o,exeub_eubtt62o,cfeub_eubtt62o

      real*8     pteub_eubtt2p(4,2), exeub_eubtt2p(1)
      complex*16 cfeub_eubtt2p(2,2)
      real*8     pteub_eubtt4m(4,2), exeub_eubtt4m(1)
      complex*16 cfeub_eubtt4m(2,2)
      real*8     pteub_eubtt9p(4,2), exeub_eubtt9p(1)
      complex*16 cfeub_eubtt9p(2,2)
      real*8     pteub_eubtt17m(4,2), exeub_eubtt17m(1)
      complex*16 cfeub_eubtt17m(2,2)
      real*8     pteub_eubtt33o(4,2), exeub_eubtt33o(1)
      complex*16 cfeub_eubtt33o(2,2)
      real*8     pteub_eubtt62o(4,2), exeub_eubtt62o(1)
      complex*16 cfeub_eubtt62o(2,2)

* Normalization
      common /aeub_eubttmdbgg/fknorm,fkcall,nkcall

* Calculated table of amplitudes
      common /aeub_eubttmatbl/av, lt, indexg
      complex*16 av(0:leub_eubttag-1)
      integer    lt(0:neub_eubttextn), indexg(neub_eubttextn)

* Spin average
      common/aeub_eubttmspin/aspin,aident,jhs(neub_eubttextn),
     &               jhe(neub_eubttextn),jcpol(neub_eubttextn)
