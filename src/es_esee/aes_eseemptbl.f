*     File es_esee/aes_eseemptbl.f : Sat Mar 18 19:45:01 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aes_eseemptbl
      implicit real*8(a-h,o-z)
      include 'incles_esee1.h'
      include 'incl2.inc'
      include 'inclk.inc'
      include 'incles_eseep.h'
*-----------------------------------------------------------------------
      jgraph = 0

* External lines
      call smextf(lincom,amdq(2),pfes_esee2,ptes_esee2t,cfes_esee2t)
      exes_esee2t(1) = lprtcl
      call smextf(loutgo,amlp(1),pfes_esee4,ptes_esee4m,cfes_esee4m)
      exes_esee4m(1) = lantip
      call smextf(loutgo,amdq(2),pfes_esee9,ptes_esee9t,cfes_esee9t)
      exes_esee9t(1) = lprtcl
      call smextf(lincom,amlp(1),pfes_esee17,ptes_esee17m,cfes_esee17m)
      exes_esee17m(1) = lantip
      call smextf(loutgo,amlp(1),pfes_esee33,ptes_esee33m,cfes_esee33m)
      exes_esee33m(1) = lprtcl
      call smextf(lincom,amlp(1),pfes_esee62,ptes_esee62m,cfes_esee62m)
      exes_esee62m(1) = lantip


* Buffer clear for amplitude
      do 200 ih = 0 , les_eseeag-1
         agc(ih) = (0.0d0,0.0d0)
  200 continue

      call aes_eseempt0

      return
      end
************************************************************************
      subroutine aes_eseempt0
      implicit real*8(a-h,o-z)

      include 'incles_esee1.h'
      include 'incl2.inc'
      include 'inclk.inc'
*-----------------------------------------------------------------------
      if(jselg(1) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 1
        call aes_eseeg1
      endif

      if(jselg(2) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 2
        call aes_eseeg2
      endif

      if(jselg(3) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 3
        call aes_eseeg3
      endif

      if(jselg(4) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 4
        call aes_eseeg4
      endif

      if(jselg(5) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 5
        call aes_eseeg5
      endif

      if(jselg(6) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 6
        call aes_eseeg6
      endif

      if(jselg(7) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 7
        call aes_eseeg7
      endif

      if(jselg(8) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 8
        call aes_eseeg8
      endif

      if(jselg(9) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 9
        call aes_eseeg9
      endif

      if(jselg(10) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 10
        call aes_eseeg10
      endif

      if(jselg(11) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 11
        call aes_eseeg11
      endif

      if(jselg(12) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 12
        call aes_eseeg12
      endif

      if(jselg(13) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 13
        call aes_eseeg13
      endif

      if(jselg(14) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 14
        call aes_eseeg14
      endif

      if(jselg(15) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 15
        call aes_eseeg15
      endif

      if(jselg(16) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 16
        call aes_eseeg16
      endif

      if(jselg(17) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 17
        call aes_eseeg17
      endif

      if(jselg(18) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 18
        call aes_eseeg18
      endif

      if(jselg(19) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 19
        call aes_eseeg19
      endif

      if(jselg(20) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 20
        call aes_eseeg20
      endif

      if(jselg(21) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 21
        call aes_eseeg21
      endif

      if(jselg(22) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 22
        call aes_eseeg22
      endif

      if(jselg(23) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 23
        call aes_eseeg23
      endif

      if(jselg(24) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 24
        call aes_eseeg24
      endif

      if(jselg(25) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 25
        call aes_eseeg25
      endif

      if(jselg(26) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 26
        call aes_eseeg26
      endif

      if(jselg(27) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 27
        call aes_eseeg27
      endif

      if(jselg(28) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 28
        call aes_eseeg28
      endif

      if(jselg(29) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 29
        call aes_eseeg29
      endif

      if(jselg(30) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 30
        call aes_eseeg30
      endif

      if(jselg(31) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 31
        call aes_eseeg31
      endif

      if(jselg(32) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 32
        call aes_eseeg32
      endif

      if(jselg(33) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 33
        call aes_eseeg33
      endif

      if(jselg(34) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 34
        call aes_eseeg34
      endif

      if(jselg(35) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 35
        call aes_eseeg35
      endif

      if(jselg(36) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 36
        call aes_eseeg36
      endif

      if(jselg(37) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 37
        call aes_eseeg37
      endif

      if(jselg(38) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 38
        call aes_eseeg38
      endif

      if(jselg(39) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 39
        call aes_eseeg39
      endif

      if(jselg(40) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 40
        call aes_eseeg40
      endif

      if(jselg(41) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 41
        call aes_eseeg41
      endif

      if(jselg(42) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 42
        call aes_eseeg42
      endif

      if(jselg(43) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 43
        call aes_eseeg43
      endif

      if(jselg(44) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 44
        call aes_eseeg44
      endif

      if(jselg(45) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 45
        call aes_eseeg45
      endif

      if(jselg(46) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 46
        call aes_eseeg46
      endif

      if(jselg(47) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 47
        call aes_eseeg47
      endif

      if(jselg(48) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 48
        call aes_eseeg48
      endif

      return
      end
