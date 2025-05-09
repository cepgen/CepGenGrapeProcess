*     File eub_eubee/aeub_eubeemptbl.f : Sat Mar 18 19:44:59 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aeub_eubeemptbl
      implicit real*8(a-h,o-z)
      include 'incleub_eubee1.h'
      include 'inclk.inc'
      include 'incleub_eubeep.h'
*-----------------------------------------------------------------------
      jgraph = 0

* External lines
      call smextf(loutgo,amuq(1),pfeub_eubee2,pteub_eubee2p,cfeub_eubee2p)
      exeub_eubee2p(1) = lantip
      call smextf(loutgo,amlp(1),pfeub_eubee4,pteub_eubee4m,cfeub_eubee4m)
      exeub_eubee4m(1) = lantip
      call smextf(lincom,amuq(1),pfeub_eubee9,pteub_eubee9p,cfeub_eubee9p)
      exeub_eubee9p(1) = lantip
      call smextf(lincom,amlp(1),pfeub_eubee17,pteub_eubee17m,cfeub_eubee17m)
      exeub_eubee17m(1) = lantip
      call smextf(loutgo,amlp(1),pfeub_eubee33,pteub_eubee33m,cfeub_eubee33m)
      exeub_eubee33m(1) = lprtcl
      call smextf(lincom,amlp(1),pfeub_eubee62,pteub_eubee62m,cfeub_eubee62m)
      exeub_eubee62m(1) = lantip


* Buffer clear for amplitude
      do 200 ih = 0 , leub_eubeeag-1
         agc(ih) = (0.0d0,0.0d0)
  200 continue

      call aeub_eubeempt0

      return
      end
************************************************************************
      subroutine aeub_eubeempt0
      implicit real*8(a-h,o-z)

      include 'incleub_eubee1.h'
      include 'inclk.inc'
*-----------------------------------------------------------------------
      if(jselg(1) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 1
        call aeub_eubeeg1
      endif

      if(jselg(2) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 2
        call aeub_eubeeg2
      endif

      if(jselg(3) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 3
        call aeub_eubeeg3
      endif

      if(jselg(4) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 4
        call aeub_eubeeg4
      endif

      if(jselg(5) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 5
        call aeub_eubeeg5
      endif

      if(jselg(6) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 6
        call aeub_eubeeg6
      endif

      if(jselg(7) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 7
        call aeub_eubeeg7
      endif

      if(jselg(8) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 8
        call aeub_eubeeg8
      endif

      if(jselg(9) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 9
        call aeub_eubeeg9
      endif

      if(jselg(10) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 10
        call aeub_eubeeg10
      endif

      if(jselg(11) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 11
        call aeub_eubeeg11
      endif

      if(jselg(12) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 12
        call aeub_eubeeg12
      endif

      if(jselg(13) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 13
        call aeub_eubeeg13
      endif

      if(jselg(14) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 14
        call aeub_eubeeg14
      endif

      if(jselg(15) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 15
        call aeub_eubeeg15
      endif

      if(jselg(16) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 16
        call aeub_eubeeg16
      endif

      if(jselg(17) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 17
        call aeub_eubeeg17
      endif

      if(jselg(18) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 18
        call aeub_eubeeg18
      endif

      if(jselg(19) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 19
        call aeub_eubeeg19
      endif

      if(jselg(20) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 20
        call aeub_eubeeg20
      endif

      if(jselg(21) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 21
        call aeub_eubeeg21
      endif

      if(jselg(22) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 22
        call aeub_eubeeg22
      endif

      if(jselg(23) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 23
        call aeub_eubeeg23
      endif

      if(jselg(24) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 24
        call aeub_eubeeg24
      endif

      if(jselg(25) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 25
        call aeub_eubeeg25
      endif

      if(jselg(26) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 26
        call aeub_eubeeg26
      endif

      if(jselg(27) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 27
        call aeub_eubeeg27
      endif

      if(jselg(28) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 28
        call aeub_eubeeg28
      endif

      if(jselg(29) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 29
        call aeub_eubeeg29
      endif

      if(jselg(30) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 30
        call aeub_eubeeg30
      endif

      if(jselg(31) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 31
        call aeub_eubeeg31
      endif

      if(jselg(32) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 32
        call aeub_eubeeg32
      endif

      if(jselg(33) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 33
        call aeub_eubeeg33
      endif

      if(jselg(34) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 34
        call aeub_eubeeg34
      endif

      if(jselg(35) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 35
        call aeub_eubeeg35
      endif

      if(jselg(36) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 36
        call aeub_eubeeg36
      endif

      if(jselg(37) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 37
        call aeub_eubeeg37
      endif

      if(jselg(38) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 38
        call aeub_eubeeg38
      endif

      if(jselg(39) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 39
        call aeub_eubeeg39
      endif

      if(jselg(40) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 40
        call aeub_eubeeg40
      endif

      if(jselg(41) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 41
        call aeub_eubeeg41
      endif

      if(jselg(42) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 42
        call aeub_eubeeg42
      endif

      if(jselg(43) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 43
        call aeub_eubeeg43
      endif

      if(jselg(44) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 44
        call aeub_eubeeg44
      endif

      if(jselg(45) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 45
        call aeub_eubeeg45
      endif

      if(jselg(46) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 46
        call aeub_eubeeg46
      endif

      if(jselg(47) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 47
        call aeub_eubeeg47
      endif

      if(jselg(48) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 48
        call aeub_eubeeg48
      endif

      return
      end
