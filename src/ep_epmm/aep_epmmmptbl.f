*     File ep_epmm/aep_epmmmptbl.f : Sat Mar 18 19:44:58 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aep_epmmmptbl
      implicit real*8(a-h,o-z)
      include 'inclep_epmm1.h'
      include 'inclk.inc'
      include 'inclep_epmmp.h'
*-----------------------------------------------------------------------
      jgraph = 0

* External lines
      call smextf(lincom,amp,pfep_epmm2,ptep_epmm2v,cfep_epmm2v)
      exep_epmm2v(1) = lprtcl
      call smextf(loutgo,amlp(1),pfep_epmm4,ptep_epmm4m,cfep_epmm4m)
      exep_epmm4m(1) = lantip
      call smextf(loutgo,amp,pfep_epmm9,ptep_epmm9v,cfep_epmm9v)
      exep_epmm9v(1) = lprtcl
      call smextf(lincom,amlp(1),pfep_epmm17,ptep_epmm17m,cfep_epmm17m)
      exep_epmm17m(1) = lantip
      call smextf(loutgo,amlp(2),pfep_epmm33,ptep_epmm33n,cfep_epmm33n)
      exep_epmm33n(1) = lprtcl
      call smextf(lincom,amlp(2),pfep_epmm62,ptep_epmm62n,cfep_epmm62n)
      exep_epmm62n(1) = lantip


* Buffer clear for amplitude
      do 200 ih = 0 , lep_epmmag-1
         agc(ih) = (0.0d0,0.0d0)
  200 continue

      call aep_epmmmpt0

      return
      end
************************************************************************
      subroutine aep_epmmmpt0
      implicit real*8(a-h,o-z)

      include 'inclep_epmm1.h'
      include 'inclk.inc'
*-----------------------------------------------------------------------
      if(jselg(1) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 1
        call aep_epmmg1
      endif

      if(jselg(2) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 2
        call aep_epmmg2
      endif

      if(jselg(3) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 3
        call aep_epmmg3
      endif

      if(jselg(4) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 4
        call aep_epmmg4
      endif

      if(jselg(5) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 5
        call aep_epmmg5
      endif

      if(jselg(6) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 6
        call aep_epmmg6
      endif

      if(jselg(7) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 7
        call aep_epmmg7
      endif

      if(jselg(8) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 8
        call aep_epmmg8
      endif

      return
      end
