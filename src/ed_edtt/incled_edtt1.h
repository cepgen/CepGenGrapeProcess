*     File ed_edtt/incled_edtt1.h : Sat Mar 18 19:45:00 2000
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
      parameter (ned_edttgraph =25)
      parameter (ned_edttextn  =6)


* number of color base
      parameter (ncbase =1)

* number of all helicity states
      parameter (led_edttag    =64)

      common /amprck/kmngr,  kmnext, kmlag

      parameter (ned_edttgrpsq = ned_edttgraph*ned_edttgraph)
      common /aed_edttmslct/jselg(ned_edttgraph),jgraph,jgluon,jhiggs

* Color string information
      common /aed_edttmcsti/ kmcbas, kmcbmx, kmcstr(2), icinfo(6), icolst
      common /aed_edttmgrph/agcwrk(0:led_edttag-1),agc(0:led_edttag-1),
     &              aprop,ancp(ned_edttgraph),ansp(0:ned_edttgraph)
     &             ,cfmtx
      common /aed_edttmgrpi/igraph(ned_edttgraph)
      complex*16 agc, agcwrk, aprop

* Momenta of external particles
      common /aed_edttmextr/peed_edtt1(4),peed_edtt2(4),peed_edtt3(4),
     &               peed_edtt4(4),peed_edtt5(4),peed_edtt6(4),
     &               prod(ned_edttextn, ned_edttextn)

      common /aed_edttmextp/pted_edtt2s,exed_edtt2s,cfed_edtt2s,
     &               pted_edtt4m,exed_edtt4m,cfed_edtt4m,pted_edtt9s,
     &               exed_edtt9s,cfed_edtt9s,pted_edtt17m,exed_edtt17m,
     &               cfed_edtt17m,pted_edtt33o,exed_edtt33o,
     &               cfed_edtt33o,pted_edtt62o,exed_edtt62o,
     &               cfed_edtt62o

      real*8     pted_edtt2s(4,2), exed_edtt2s(1)
      complex*16 cfed_edtt2s(2,2)
      real*8     pted_edtt4m(4,2), exed_edtt4m(1)
      complex*16 cfed_edtt4m(2,2)
      real*8     pted_edtt9s(4,2), exed_edtt9s(1)
      complex*16 cfed_edtt9s(2,2)
      real*8     pted_edtt17m(4,2), exed_edtt17m(1)
      complex*16 cfed_edtt17m(2,2)
      real*8     pted_edtt33o(4,2), exed_edtt33o(1)
      complex*16 cfed_edtt33o(2,2)
      real*8     pted_edtt62o(4,2), exed_edtt62o(1)
      complex*16 cfed_edtt62o(2,2)

* Normalization
      common /aed_edttmdbgg/fknorm,fkcall,nkcall

* Calculated table of amplitudes
      common /aed_edttmatbl/av, lt, indexg
      complex*16 av(0:led_edttag-1)
      integer    lt(0:ned_edttextn), indexg(ned_edttextn)

* Spin average
      common/aed_edttmspin/aspin,aident,jhs(ned_edttextn),
     &               jhe(ned_edttextn),jcpol(ned_edttextn)
