*     File ecb_ecbtt/aecb_ecbttmptbl.f : Sat Mar 18 19:45:04 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aecb_ecbttmptbl
      implicit real*8(a-h,o-z)
      include 'inclecb_ecbtt1.h'
      include 'inclk.inc'
      include 'inclecb_ecbttp.h'
*-----------------------------------------------------------------------
      jgraph = 0

* External lines
      call smextf(loutgo,amuq(2),pfecb_ecbtt2,ptecb_ecbtt2q,cfecb_ecbtt2q)
      execb_ecbtt2q(1) = lantip
      call smextf(loutgo,amlp(1),pfecb_ecbtt4,ptecb_ecbtt4m,cfecb_ecbtt4m)
      execb_ecbtt4m(1) = lantip
      call smextf(lincom,amuq(2),pfecb_ecbtt9,ptecb_ecbtt9q,cfecb_ecbtt9q)
      execb_ecbtt9q(1) = lantip
      call smextf(lincom,amlp(1),pfecb_ecbtt17,ptecb_ecbtt17m,cfecb_ecbtt17m)
      execb_ecbtt17m(1) = lantip
      call smextf(loutgo,amlp(3),pfecb_ecbtt33,ptecb_ecbtt33o,cfecb_ecbtt33o)
      execb_ecbtt33o(1) = lprtcl
      call smextf(lincom,amlp(3),pfecb_ecbtt62,ptecb_ecbtt62o,cfecb_ecbtt62o)
      execb_ecbtt62o(1) = lantip


* Buffer clear for amplitude
      do 200 ih = 0 , lecb_ecbttag-1
         agc(ih) = (0.0d0,0.0d0)
  200 continue

      call aecb_ecbttmpt0

      return
      end
************************************************************************
      subroutine aecb_ecbttmpt0
      implicit real*8(a-h,o-z)

      include 'inclecb_ecbtt1.h'
      include 'inclk.inc'
*-----------------------------------------------------------------------
      if(jselg(1) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 1
        call aecb_ecbttg1
      endif

      if(jselg(2) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 2
        call aecb_ecbttg2
      endif

      if(jselg(3) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 3
        call aecb_ecbttg3
      endif

      if(jselg(4) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 4
        call aecb_ecbttg4
      endif

      if(jselg(5) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 5
        call aecb_ecbttg5
      endif

      if(jselg(6) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 6
        call aecb_ecbttg6
      endif

      if(jselg(7) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 7
        call aecb_ecbttg7
      endif

      if(jselg(8) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 8
        call aecb_ecbttg8
      endif

      if(jselg(9) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 9
        call aecb_ecbttg9
      endif

      if(jselg(10) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 10
        call aecb_ecbttg10
      endif

      if(jselg(11) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 11
        call aecb_ecbttg11
      endif

      if(jselg(12) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 12
        call aecb_ecbttg12
      endif

      if(jhiggs .ne. 0) then
      if(jselg(13) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 13
        call aecb_ecbttg13
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(14) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 14
        call aecb_ecbttg14
      endif
      endif

      if(jselg(15) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 15
        call aecb_ecbttg15
      endif

      if(jselg(16) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 16
        call aecb_ecbttg16
      endif

      if(jselg(17) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 17
        call aecb_ecbttg17
      endif

      if(jselg(18) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 18
        call aecb_ecbttg18
      endif

      if(jhiggs .ne. 0) then
      if(jselg(19) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 19
        call aecb_ecbttg19
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(20) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 20
        call aecb_ecbttg20
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(21) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 21
        call aecb_ecbttg21
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(22) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 22
        call aecb_ecbttg22
      endif
      endif

      if(jhiggs .ne. 0) then
      if(igauzb .ne. 0) then
      if(jselg(23) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 23
        call aecb_ecbttg23
      endif
      endif
      endif

      if(jselg(24) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 24
        call aecb_ecbttg24
      endif

      if(jselg(25) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 25
        call aecb_ecbttg25
      endif

      if(jhiggs .ne. 0) then
      if(jselg(26) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 26
        call aecb_ecbttg26
      endif
      endif

      if(jselg(27) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 27
        call aecb_ecbttg27
      endif

      if(jselg(28) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 28
        call aecb_ecbttg28
      endif

      if(jhiggs .ne. 0) then
      if(jselg(29) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 29
        call aecb_ecbttg29
      endif
      endif

      if(jselg(30) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 30
        call aecb_ecbttg30
      endif

      if(jselg(31) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 31
        call aecb_ecbttg31
      endif

      if(jselg(32) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 32
        call aecb_ecbttg32
      endif

      if(jselg(33) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 33
        call aecb_ecbttg33
      endif

      if(jhiggs .ne. 0) then
      if(jselg(34) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 34
        call aecb_ecbttg34
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(35) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 35
        call aecb_ecbttg35
      endif
      endif

      return
      end
