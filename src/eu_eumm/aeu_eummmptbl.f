*     File eu_eumm/aeu_eummmptbl.f : Sat Mar 18 19:44:59 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aeu_eummmptbl
      implicit real*8(a-h,o-z)
      include 'incleu_eumm1.h'
      include 'inclk.inc'
      include 'incleu_eummp.h'
*-----------------------------------------------------------------------
      jgraph = 0

* External lines
      call smextf(lincom,amuq(1),pfeu_eumm2,pteu_eumm2p,cfeu_eumm2p)
      exeu_eumm2p(1) = lprtcl
      call smextf(loutgo,amlp(1),pfeu_eumm4,pteu_eumm4m,cfeu_eumm4m)
      exeu_eumm4m(1) = lantip
      call smextf(loutgo,amuq(1),pfeu_eumm9,pteu_eumm9p,cfeu_eumm9p)
      exeu_eumm9p(1) = lprtcl
      call smextf(lincom,amlp(1),pfeu_eumm17,pteu_eumm17m,cfeu_eumm17m)
      exeu_eumm17m(1) = lantip
      call smextf(loutgo,amlp(2),pfeu_eumm33,pteu_eumm33n,cfeu_eumm33n)
      exeu_eumm33n(1) = lprtcl
      call smextf(lincom,amlp(2),pfeu_eumm62,pteu_eumm62n,cfeu_eumm62n)
      exeu_eumm62n(1) = lantip


* Buffer clear for amplitude
      do 200 ih = 0 , leu_eummag-1
         agc(ih) = (0.0d0,0.0d0)
  200 continue

      call aeu_eummmpt0

      return
      end
************************************************************************
      subroutine aeu_eummmpt0
      implicit real*8(a-h,o-z)

      include 'incleu_eumm1.h'
      include 'inclk.inc'
*-----------------------------------------------------------------------
      if(jselg(1) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 1
        call aeu_eummg1
      endif

      if(jselg(2) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 2
        call aeu_eummg2
      endif

      if(jselg(3) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 3
        call aeu_eummg3
      endif

      if(jselg(4) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 4
        call aeu_eummg4
      endif

      if(jselg(5) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 5
        call aeu_eummg5
      endif

      if(jselg(6) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 6
        call aeu_eummg6
      endif

      if(jselg(7) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 7
        call aeu_eummg7
      endif

      if(jselg(8) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 8
        call aeu_eummg8
      endif

      if(jselg(9) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 9
        call aeu_eummg9
      endif

      if(jselg(10) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 10
        call aeu_eummg10
      endif

      if(jselg(11) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 11
        call aeu_eummg11
      endif

      if(jselg(12) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 12
        call aeu_eummg12
      endif

      if(jselg(13) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 13
        call aeu_eummg13
      endif

      if(jselg(14) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 14
        call aeu_eummg14
      endif

      if(jselg(15) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 15
        call aeu_eummg15
      endif

      if(jselg(16) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 16
        call aeu_eummg16
      endif

      if(jhiggs .ne. 0) then
      if(jselg(17) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 17
        call aeu_eummg17
      endif
      endif

      if(jselg(18) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 18
        call aeu_eummg18
      endif

      if(jselg(19) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 19
        call aeu_eummg19
      endif

      if(jselg(20) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 20
        call aeu_eummg20
      endif

      if(jselg(21) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 21
        call aeu_eummg21
      endif

      if(jselg(22) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 22
        call aeu_eummg22
      endif

      if(jselg(23) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 23
        call aeu_eummg23
      endif

      if(jselg(24) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 24
        call aeu_eummg24
      endif

      if(jselg(25) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 25
        call aeu_eummg25
      endif

      return
      end
