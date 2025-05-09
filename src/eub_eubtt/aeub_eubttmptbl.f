*     File eub_eubtt/aeub_eubttmptbl.f : Sat Mar 18 19:45:00 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aeub_eubttmptbl
      implicit real*8(a-h,o-z)
      include 'incleub_eubtt1.h'
      include 'inclk.inc'
      include 'incleub_eubttp.h'
*-----------------------------------------------------------------------
      jgraph = 0

* External lines
      call smextf(loutgo,amuq(1),pfeub_eubtt2,pteub_eubtt2p,cfeub_eubtt2p)
      exeub_eubtt2p(1) = lantip
      call smextf(loutgo,amlp(1),pfeub_eubtt4,pteub_eubtt4m,cfeub_eubtt4m)
      exeub_eubtt4m(1) = lantip
      call smextf(lincom,amuq(1),pfeub_eubtt9,pteub_eubtt9p,cfeub_eubtt9p)
      exeub_eubtt9p(1) = lantip
      call smextf(lincom,amlp(1),pfeub_eubtt17,pteub_eubtt17m,cfeub_eubtt17m)
      exeub_eubtt17m(1) = lantip
      call smextf(loutgo,amlp(3),pfeub_eubtt33,pteub_eubtt33o,cfeub_eubtt33o)
      exeub_eubtt33o(1) = lprtcl
      call smextf(lincom,amlp(3),pfeub_eubtt62,pteub_eubtt62o,cfeub_eubtt62o)
      exeub_eubtt62o(1) = lantip


* Buffer clear for amplitude
      do 200 ih = 0 , leub_eubttag-1
         agc(ih) = (0.0d0,0.0d0)
  200 continue

      call aeub_eubttmpt0

      return
      end
************************************************************************
      subroutine aeub_eubttmpt0
      implicit real*8(a-h,o-z)

      include 'incleub_eubtt1.h'
      include 'inclk.inc'
*-----------------------------------------------------------------------
      if(jselg(1) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 1
        call aeub_eubttg1
      endif

      if(jselg(2) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 2
        call aeub_eubttg2
      endif

      if(jselg(3) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 3
        call aeub_eubttg3
      endif

      if(jselg(4) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 4
        call aeub_eubttg4
      endif

      if(jselg(5) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 5
        call aeub_eubttg5
      endif

      if(jselg(6) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 6
        call aeub_eubttg6
      endif

      if(jselg(7) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 7
        call aeub_eubttg7
      endif

      if(jselg(8) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 8
        call aeub_eubttg8
      endif

      if(jselg(9) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 9
        call aeub_eubttg9
      endif

      if(jselg(10) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 10
        call aeub_eubttg10
      endif

      if(jselg(11) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 11
        call aeub_eubttg11
      endif

      if(jselg(12) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 12
        call aeub_eubttg12
      endif

      if(jselg(13) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 13
        call aeub_eubttg13
      endif

      if(jselg(14) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 14
        call aeub_eubttg14
      endif

      if(jselg(15) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 15
        call aeub_eubttg15
      endif

      if(jselg(16) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 16
        call aeub_eubttg16
      endif

      if(jhiggs .ne. 0) then
      if(jselg(17) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 17
        call aeub_eubttg17
      endif
      endif

      if(jselg(18) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 18
        call aeub_eubttg18
      endif

      if(jselg(19) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 19
        call aeub_eubttg19
      endif

      if(jselg(20) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 20
        call aeub_eubttg20
      endif

      if(jselg(21) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 21
        call aeub_eubttg21
      endif

      if(jselg(22) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 22
        call aeub_eubttg22
      endif

      if(jselg(23) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 23
        call aeub_eubttg23
      endif

      if(jselg(24) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 24
        call aeub_eubttg24
      endif

      if(jselg(25) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 25
        call aeub_eubttg25
      endif

      return
      end
