*     File ed_edmm/incled_edmm1.h : Sat Mar 18 19:45:00 2000
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
      parameter (ned_edmmgraph =25)
      parameter (ned_edmmextn  =6)


* number of color base
      parameter (ncbase =1)

* number of all helicity states
      parameter (led_edmmag    =64)

      common /amprck/kmngr,  kmnext, kmlag

      parameter (ned_edmmgrpsq = ned_edmmgraph*ned_edmmgraph)
      common /aed_edmmmslct/jselg(ned_edmmgraph),jgraph,jgluon,jhiggs

* Color string information
      common /aed_edmmmcsti/ kmcbas, kmcbmx, kmcstr(2), icinfo(6), icolst
      common /aed_edmmmgrph/agcwrk(0:led_edmmag-1),agc(0:led_edmmag-1),
     &              aprop,ancp(ned_edmmgraph),ansp(0:ned_edmmgraph)
     &             ,cfmtx
      common /aed_edmmmgrpi/igraph(ned_edmmgraph)
      complex*16 agc, agcwrk, aprop

* Momenta of external particles
      common /aed_edmmmextr/peed_edmm1(4),peed_edmm2(4),peed_edmm3(4),
     &               peed_edmm4(4),peed_edmm5(4),peed_edmm6(4),
     &               prod(ned_edmmextn, ned_edmmextn)

      common /aed_edmmmextp/pted_edmm2s,exed_edmm2s,cfed_edmm2s,
     &               pted_edmm4m,exed_edmm4m,cfed_edmm4m,pted_edmm9s,
     &               exed_edmm9s,cfed_edmm9s,pted_edmm17m,exed_edmm17m,
     &               cfed_edmm17m,pted_edmm33n,exed_edmm33n,
     &               cfed_edmm33n,pted_edmm62n,exed_edmm62n,
     &               cfed_edmm62n

      real*8     pted_edmm2s(4,2), exed_edmm2s(1)
      complex*16 cfed_edmm2s(2,2)
      real*8     pted_edmm4m(4,2), exed_edmm4m(1)
      complex*16 cfed_edmm4m(2,2)
      real*8     pted_edmm9s(4,2), exed_edmm9s(1)
      complex*16 cfed_edmm9s(2,2)
      real*8     pted_edmm17m(4,2), exed_edmm17m(1)
      complex*16 cfed_edmm17m(2,2)
      real*8     pted_edmm33n(4,2), exed_edmm33n(1)
      complex*16 cfed_edmm33n(2,2)
      real*8     pted_edmm62n(4,2), exed_edmm62n(1)
      complex*16 cfed_edmm62n(2,2)

* Normalization
      common /aed_edmmmdbgg/fknorm,fkcall,nkcall

* Calculated table of amplitudes
      common /aed_edmmmatbl/av, lt, indexg
      complex*16 av(0:led_edmmag-1)
      integer    lt(0:ned_edmmextn), indexg(ned_edmmextn)

* Spin average
      common/aed_edmmmspin/aspin,aident,jhs(ned_edmmextn),
     &               jhe(ned_edmmextn),jcpol(ned_edmmextn)
