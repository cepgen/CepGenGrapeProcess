*     File etb_etbtt/aetb_etbttmptbl.f : Sat Mar 18 19:45:07 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aetb_etbttmptbl
      implicit real*8(a-h,o-z)
      include 'incletb_etbtt1.h'
      include 'inclk.inc'
      include 'incletb_etbttp.h'
*-----------------------------------------------------------------------
      jgraph = 0

* External lines
      call smextf(loutgo,amuq(3),pfetb_etbtt2,ptetb_etbtt2r,cfetb_etbtt2r)
      exetb_etbtt2r(1) = lantip
      call smextf(loutgo,amlp(1),pfetb_etbtt4,ptetb_etbtt4m,cfetb_etbtt4m)
      exetb_etbtt4m(1) = lantip
      call smextf(lincom,amuq(3),pfetb_etbtt9,ptetb_etbtt9r,cfetb_etbtt9r)
      exetb_etbtt9r(1) = lantip
      call smextf(lincom,amlp(1),pfetb_etbtt17,ptetb_etbtt17m,cfetb_etbtt17m)
      exetb_etbtt17m(1) = lantip
      call smextf(loutgo,amlp(3),pfetb_etbtt33,ptetb_etbtt33o,cfetb_etbtt33o)
      exetb_etbtt33o(1) = lprtcl
      call smextf(lincom,amlp(3),pfetb_etbtt62,ptetb_etbtt62o,cfetb_etbtt62o)
      exetb_etbtt62o(1) = lantip


* Buffer clear for amplitude
      do 200 ih = 0 , letb_etbttag-1
         agc(ih) = (0.0d0,0.0d0)
  200 continue

      call aetb_etbttmpt0

      return
      end
************************************************************************
      subroutine aetb_etbttmpt0
      implicit real*8(a-h,o-z)

      include 'incletb_etbtt1.h'
      include 'inclk.inc'
*-----------------------------------------------------------------------
      if(jselg(1) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 1
        call aetb_etbttg1
      endif

      if(jselg(2) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 2
        call aetb_etbttg2
      endif

      if(jselg(3) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 3
        call aetb_etbttg3
      endif

      if(jselg(4) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 4
        call aetb_etbttg4
      endif

      if(jselg(5) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 5
        call aetb_etbttg5
      endif

      if(jselg(6) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 6
        call aetb_etbttg6
      endif

      if(jselg(7) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 7
        call aetb_etbttg7
      endif

      if(jselg(8) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 8
        call aetb_etbttg8
      endif

      if(jselg(9) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 9
        call aetb_etbttg9
      endif

      if(jselg(10) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 10
        call aetb_etbttg10
      endif

      if(jselg(11) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 11
        call aetb_etbttg11
      endif

      if(jselg(12) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 12
        call aetb_etbttg12
      endif

      if(jhiggs .ne. 0) then
      if(jselg(13) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 13
        call aetb_etbttg13
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(14) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 14
        call aetb_etbttg14
      endif
      endif

      if(jselg(15) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 15
        call aetb_etbttg15
      endif

      if(jselg(16) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 16
        call aetb_etbttg16
      endif

      if(jselg(17) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 17
        call aetb_etbttg17
      endif

      if(jselg(18) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 18
        call aetb_etbttg18
      endif

      if(jhiggs .ne. 0) then
      if(jselg(19) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 19
        call aetb_etbttg19
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(20) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 20
        call aetb_etbttg20
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(21) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 21
        call aetb_etbttg21
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(22) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 22
        call aetb_etbttg22
      endif
      endif

      if(jhiggs .ne. 0) then
      if(igauzb .ne. 0) then
      if(jselg(23) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 23
        call aetb_etbttg23
      endif
      endif
      endif

      if(jselg(24) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 24
        call aetb_etbttg24
      endif

      if(jselg(25) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 25
        call aetb_etbttg25
      endif

      if(jhiggs .ne. 0) then
      if(jselg(26) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 26
        call aetb_etbttg26
      endif
      endif

      if(jselg(27) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 27
        call aetb_etbttg27
      endif

      if(jselg(28) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 28
        call aetb_etbttg28
      endif

      if(jhiggs .ne. 0) then
      if(jselg(29) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 29
        call aetb_etbttg29
      endif
      endif

      if(jselg(30) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 30
        call aetb_etbttg30
      endif

      if(jselg(31) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 31
        call aetb_etbttg31
      endif

      if(jselg(32) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 32
        call aetb_etbttg32
      endif

      if(jselg(33) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 33
        call aetb_etbttg33
      endif

      if(jhiggs .ne. 0) then
      if(jselg(34) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 34
        call aetb_etbttg34
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(35) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 35
        call aetb_etbttg35
      endif
      endif

      return
      end
