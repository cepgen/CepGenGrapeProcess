*     File edb_edbtt/incledb_edbtt1.h : Sat Mar 18 19:45:01 2000
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
      parameter (nedb_edbttgraph =25)
      parameter (nedb_edbttextn  =6)


* number of color base
      parameter (ncbase =1)

* number of all helicity states
      parameter (ledb_edbttag    =64)

      common /amprck/kmngr,  kmnext, kmlag
      common /amprcs/smpref, smproc
      character*80   smpref, smproc

      parameter (nedb_edbttgrpsq = nedb_edbttgraph*nedb_edbttgraph)
      common /aedb_edbttmslct/jselg(nedb_edbttgraph),jgraph,jgluon,
     &               jhiggs

* Color string information
      common /aedb_edbttmcsti/ kmcbas, kmcbmx, kmcstr(2), icinfo(6), icolst
      common /aedb_edbttmgrph/agcwrk(0:ledb_edbttag-1),agc(0:ledb_edbttag-1),
     &              aprop,ancp(nedb_edbttgraph),ansp(0:nedb_edbttgraph)
     &             ,cfmtx
      common /aedb_edbttmgrpi/igraph(nedb_edbttgraph)
      complex*16 agc, agcwrk, aprop

* Momenta of external particles
      common /aedb_edbttmextr/peedb_edbtt1(4),peedb_edbtt2(4),
     &               peedb_edbtt3(4),peedb_edbtt4(4),peedb_edbtt5(4),
     &               peedb_edbtt6(4),
     &               prod(nedb_edbttextn, nedb_edbttextn)

      common /aedb_edbttmextp/ptedb_edbtt2s,exedb_edbtt2s,cfedb_edbtt2s,
     &               ptedb_edbtt4m,exedb_edbtt4m,cfedb_edbtt4m,
     &               ptedb_edbtt9s,exedb_edbtt9s,cfedb_edbtt9s,
     &               ptedb_edbtt17m,exedb_edbtt17m,cfedb_edbtt17m,
     &               ptedb_edbtt33o,exedb_edbtt33o,cfedb_edbtt33o,
     &               ptedb_edbtt62o,exedb_edbtt62o,cfedb_edbtt62o

      real*8     ptedb_edbtt2s(4,2), exedb_edbtt2s(1)
      complex*16 cfedb_edbtt2s(2,2)
      real*8     ptedb_edbtt4m(4,2), exedb_edbtt4m(1)
      complex*16 cfedb_edbtt4m(2,2)
      real*8     ptedb_edbtt9s(4,2), exedb_edbtt9s(1)
      complex*16 cfedb_edbtt9s(2,2)
      real*8     ptedb_edbtt17m(4,2), exedb_edbtt17m(1)
      complex*16 cfedb_edbtt17m(2,2)
      real*8     ptedb_edbtt33o(4,2), exedb_edbtt33o(1)
      complex*16 cfedb_edbtt33o(2,2)
      real*8     ptedb_edbtt62o(4,2), exedb_edbtt62o(1)
      complex*16 cfedb_edbtt62o(2,2)

* Normalization
      common /aedb_edbttmdbgg/fknorm,fkcall,nkcall

* Calculated table of amplitudes
      common /aedb_edbttmatbl/av, lt, indexg
      complex*16 av(0:ledb_edbttag-1)
      integer    lt(0:nedb_edbttextn), indexg(nedb_edbttextn)

* Spin average
      common/aedb_edbttmspin/aspin,aident,jhs(nedb_edbttextn),
     &               jhe(nedb_edbttextn),jcpol(nedb_edbttextn)
