*     File ecb_ecbee/aecb_ecbeemptbl.f : Sat Mar 18 19:45:04 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aecb_ecbeemptbl
      implicit real*8(a-h,o-z)
      include 'inclecb_ecbee1.h'
      include 'incl2.inc'
      include 'inclk.inc'
      include 'inclecb_ecbeep.h'
*-----------------------------------------------------------------------
      jgraph = 0

* External lines
      call smextf(loutgo,amuq(2),pfecb_ecbee2,ptecb_ecbee2q,cfecb_ecbee2q)
      execb_ecbee2q(1) = lantip
      call smextf(loutgo,amlp(1),pfecb_ecbee4,ptecb_ecbee4m,cfecb_ecbee4m)
      execb_ecbee4m(1) = lantip
      call smextf(lincom,amuq(2),pfecb_ecbee9,ptecb_ecbee9q,cfecb_ecbee9q)
      execb_ecbee9q(1) = lantip
      call smextf(lincom,amlp(1),pfecb_ecbee17,ptecb_ecbee17m,cfecb_ecbee17m)
      execb_ecbee17m(1) = lantip
      call smextf(loutgo,amlp(1),pfecb_ecbee33,ptecb_ecbee33m,cfecb_ecbee33m)
      execb_ecbee33m(1) = lprtcl
      call smextf(lincom,amlp(1),pfecb_ecbee62,ptecb_ecbee62m,cfecb_ecbee62m)
      execb_ecbee62m(1) = lantip


* Buffer clear for amplitude
      do 200 ih = 0 , lecb_ecbeeag-1
         agc(ih) = (0.0d0,0.0d0)
  200 continue

      call aecb_ecbeempt0

      return
      end
************************************************************************
      subroutine aecb_ecbeempt0
      implicit real*8(a-h,o-z)

      include 'inclecb_ecbee1.h'
      include 'incl2.inc'
      include 'inclk.inc'
*-----------------------------------------------------------------------
      if(jselg(1) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 1
        call aecb_ecbeeg1
      endif

      if(jselg(2) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 2
        call aecb_ecbeeg2
      endif

      if(jselg(3) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 3
        call aecb_ecbeeg3
      endif

      if(jselg(4) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 4
        call aecb_ecbeeg4
      endif

      if(jselg(5) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 5
        call aecb_ecbeeg5
      endif

      if(jselg(6) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 6
        call aecb_ecbeeg6
      endif

      if(jselg(7) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 7
        call aecb_ecbeeg7
      endif

      if(jselg(8) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 8
        call aecb_ecbeeg8
      endif

      if(jselg(9) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 9
        call aecb_ecbeeg9
      endif

      if(jselg(10) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 10
        call aecb_ecbeeg10
      endif

      if(jselg(11) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 11
        call aecb_ecbeeg11
      endif

      if(jselg(12) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 12
        call aecb_ecbeeg12
      endif

      if(jselg(13) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 13
        call aecb_ecbeeg13
      endif

      if(jselg(14) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 14
        call aecb_ecbeeg14
      endif

      if(jselg(15) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 15
        call aecb_ecbeeg15
      endif

      if(jselg(16) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 16
        call aecb_ecbeeg16
      endif

      if(jselg(17) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 17
        call aecb_ecbeeg17
      endif

      if(jselg(18) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 18
        call aecb_ecbeeg18
      endif

      if(jselg(19) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 19
        call aecb_ecbeeg19
      endif

      if(jselg(20) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 20
        call aecb_ecbeeg20
      endif

      if(jselg(21) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 21
        call aecb_ecbeeg21
      endif

      if(jselg(22) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 22
        call aecb_ecbeeg22
      endif

      if(jselg(23) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 23
        call aecb_ecbeeg23
      endif

      if(jselg(24) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 24
        call aecb_ecbeeg24
      endif

      if(jselg(25) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 25
        call aecb_ecbeeg25
      endif

      if(jselg(26) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 26
        call aecb_ecbeeg26
      endif

      if(jselg(27) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 27
        call aecb_ecbeeg27
      endif

      if(jselg(28) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 28
        call aecb_ecbeeg28
      endif

      if(jselg(29) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 29
        call aecb_ecbeeg29
      endif

      if(jselg(30) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 30
        call aecb_ecbeeg30
      endif

      if(jselg(31) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 31
        call aecb_ecbeeg31
      endif

      if(jselg(32) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 32
        call aecb_ecbeeg32
      endif

      if(jhiggs .ne. 0) then
      if(jselg(33) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 33
        call aecb_ecbeeg33
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(34) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 34
        call aecb_ecbeeg34
      endif
      endif

      if(jselg(35) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 35
        call aecb_ecbeeg35
      endif

      if(jselg(36) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 36
        call aecb_ecbeeg36
      endif

      if(jselg(37) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 37
        call aecb_ecbeeg37
      endif

      if(jselg(38) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 38
        call aecb_ecbeeg38
      endif

      if(jselg(39) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 39
        call aecb_ecbeeg39
      endif

      if(jselg(40) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 40
        call aecb_ecbeeg40
      endif

      if(jselg(41) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 41
        call aecb_ecbeeg41
      endif

      if(jselg(42) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 42
        call aecb_ecbeeg42
      endif

      if(jselg(43) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 43
        call aecb_ecbeeg43
      endif

      if(jselg(44) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 44
        call aecb_ecbeeg44
      endif

      if(jselg(45) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 45
        call aecb_ecbeeg45
      endif

      if(jselg(46) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 46
        call aecb_ecbeeg46
      endif

      if(jselg(47) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 47
        call aecb_ecbeeg47
      endif

      if(jselg(48) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 48
        call aecb_ecbeeg48
      endif

      if(jselg(49) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 49
        call aecb_ecbeeg49
      endif

      if(jselg(50) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 50
        call aecb_ecbeeg50
      endif

      return
      end