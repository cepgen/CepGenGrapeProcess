*     File ep_eptt/aep_epttmptbl.f : Sat Mar 18 19:44:58 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aep_epttmptbl
      implicit real*8(a-h,o-z)
      include 'inclep_eptt1.h'
      include 'incl2.inc'
      include 'inclk.inc'
      include 'inclep_epttp.h'
*-----------------------------------------------------------------------
      jgraph = 0

* External lines
      call smextf(lincom,amp,pfep_eptt2,ptep_eptt2v,cfep_eptt2v)
      exep_eptt2v(1) = lprtcl
      call smextf(loutgo,amlp(1),pfep_eptt4,ptep_eptt4m,cfep_eptt4m)
      exep_eptt4m(1) = lantip
      call smextf(loutgo,amp,pfep_eptt9,ptep_eptt9v,cfep_eptt9v)
      exep_eptt9v(1) = lprtcl
      call smextf(lincom,amlp(1),pfep_eptt17,ptep_eptt17m,cfep_eptt17m)
      exep_eptt17m(1) = lantip
      call smextf(loutgo,amlp(3),pfep_eptt33,ptep_eptt33o,cfep_eptt33o)
      exep_eptt33o(1) = lprtcl
      call smextf(lincom,amlp(3),pfep_eptt62,ptep_eptt62o,cfep_eptt62o)
      exep_eptt62o(1) = lantip


* Buffer clear for amplitude
      do 200 ih = 0 , lep_epttag-1
         agc(ih) = (0.0d0,0.0d0)
  200 continue

      call aep_epttmpt0

      return
      end
************************************************************************
      subroutine aep_epttmpt0
      implicit real*8(a-h,o-z)

      include 'inclep_eptt1.h'
      include 'incl2.inc'
      include 'inclk.inc'
*-----------------------------------------------------------------------
      if(jselg(1) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 1
        call aep_epttg1
      endif

      if(jselg(2) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 2
        call aep_epttg2
      endif

      if(jselg(3) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 3
        call aep_epttg3
      endif

      if(jselg(4) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 4
        call aep_epttg4
      endif

      if(jselg(5) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 5
        call aep_epttg5
      endif

      if(jselg(6) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 6
        call aep_epttg6
      endif

      if(jselg(7) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 7
        call aep_epttg7
      endif

      if(jselg(8) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 8
        call aep_epttg8
      endif

      return
      end
