*     File ep_epmm/aep_epmmmps.f : Sat Mar 18 19:44:58 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
*             Graph No. 1 - 1
*         Generated No. 1
************************************************************************
      subroutine aep_epmmg1
      implicit real*8(a-h,o-z)

      include 'inclep_epmm1.h'
      include 'inclk.inc'
      include 'inclep_epmmp.h'
*-----------------------------------------------------------------------
      common /amwork/cfep_epmm15m,av6,av7,av8,av9,av10,av11,
     &               evep_epmm11d,eqep_epmm11d,exep_epmm15m,
     &               ptep_epmm15m,evep_epmm31d,eqep_epmm31d
      common /amwori/lt6,lt7,lt8,lt9,lt10,lt11
*     3632 (3632) + 108 (108) bytes used

      integer    lt6(0:3),lt7(0:3),lt8(0:3),lt9(0:3),lt10(0:4),
     &           lt11(0:5)
      real*8     evep_epmm11d(lepina),eqep_epmm11d(4,lepina),
     &           exep_epmm15m(2),ptep_epmm15m(4,3),
     &           evep_epmm31d(lepina),eqep_epmm31d(4,lepina)
      complex*16 cfep_epmm15m(2,4)
      complex*16 av6(lextrn*lextrn*lepina)
      complex*16 av7(lintrn*lextrn*lepina)
      complex*16 av8(lextrn*lintrn*lepina)
      complex*16 av9(lextrn*lextrn*lepina)
      complex*16 av10(lextrn*lextrn*lintrn*lextrn)
      complex*16 av11(lextrn*lepina*lextrn*lextrn*lextrn)
      complex*16 atmp
*-----------------------------------------------------------------------

* Denominators of propagators
      aprop = 1.0d0
      call snprpd(pphase,aprop,vnep_epmm11,
     &      ama**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm15,
     &      amlp(1)**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm31,
     &      ama**2,0.0d0)

* Internal momenta
      call smintv(lepina,ama,pfep_epmm11,eqep_epmm11d,evep_epmm11d,
     &            vnep_epmm11,igauab)
      call smintf(amlp(1),pfep_epmm15,vnep_epmm15,exep_epmm15m,
     &            ptep_epmm15m,cfep_epmm15m)
      call smintv(lepina,ama,pfep_epmm31,eqep_epmm31d,evep_epmm31d,
     &            vnep_epmm31,igauab)

* Vertices (10)

*     6(0): pfep_epmm9 p
*     6(1): pfep_epmm2 p
*     6(2): pfep_epmm11 photon
      call stppv(2212,lextrn,amp,pfep_epmm2,+1,vnep_epmm2,exep_epmm2v,
     &            cfep_epmm2v,ptep_epmm2v,2212,lextrn,amp,pfep_epmm9,-1,
     &            vnep_epmm9,exep_epmm9v,cfep_epmm9v,ptep_epmm9v,22,
     &            lepina,ama,pfep_epmm11,+1,vnep_epmm11,exep_epmm11d,
     &            eqep_epmm11d,capp,lt6,av6)

*     7(0): pfep_epmm4 electron
*     7(1): pfep_epmm15 electron
*     7(2): pfep_epmm11 photon
      call smffv(lintrn,lextrn,lepina,exep_epmm15m,exep_epmm4m,amlp(1),
     &           amlp(1),call(1,1),cfep_epmm15m,cfep_epmm4m,
     &           ptep_epmm15m,ptep_epmm4m,eqep_epmm11d,lt7,av7)

*     8(0): pfep_epmm15 electron
*     8(1): pfep_epmm17 electron
*     8(2): pfep_epmm31 photon
      call smffv(lextrn,lintrn,lepina,exep_epmm17m,exep_epmm15m,amlp(1),
     &           amlp(1),call(1,1),cfep_epmm17m,cfep_epmm15m,
     &           ptep_epmm17m,ptep_epmm15m,eqep_epmm31d,lt8,av8)

*     9(0): pfep_epmm33 muon
*     9(1): pfep_epmm62 muon
*     9(2): pfep_epmm31 photon
      call smffv(lextrn,lextrn,lepina,exep_epmm62n,exep_epmm33n,amlp(2),
     &           amlp(2),call(1,2),cfep_epmm62n,cfep_epmm33n,
     &           ptep_epmm62n,ptep_epmm33n,eqep_epmm31d,lt9,av9)

      call smconv(lt6,lt7,3,3,evep_epmm11d,av6,av7,lt10,av10)
      call smconf(lt8,lt10,2,3,exep_epmm15m,av8,av10,lt11,av11)
      call smconv(lt9,lt11,3,2,evep_epmm31d,av9,av11,lt,av)

      sym = + 1.0d0
      cf  = + 1.0d0
      aprop         = cf*sym/aprop

      indexg(1) = 6
      indexg(2) = 5
      indexg(3) = 4
      indexg(4) = 1
      indexg(5) = 3
      indexg(6) = 2

      call aep_epmmmpord(lt, av, indexg, agcwrk)

      ancp(jgraph) = 0.0d0
*     nbase = 1
      do 500 ih = 0 , lep_epmmag-1
         atmp    = agcwrk(ih)*aprop
         agc(ih) = agc(ih) + atmp
         ancp(jgraph) = ancp(jgraph) + atmp*conjg(atmp)
  500 continue

      return
      end
************************************************************************
*             Graph No. 2 - 1
*         Generated No. 2
************************************************************************
      subroutine aep_epmmg2
      implicit real*8(a-h,o-z)

      include 'inclep_epmm1.h'
      include 'inclk.inc'
      include 'inclep_epmmp.h'
*-----------------------------------------------------------------------
      common /amwork/cfep_epmm15m,av6,av7,av8,av9,av10,av11,
     &               evep_epmm11d,eqep_epmm11d,exep_epmm15m,
     &               ptep_epmm15m,evep_epmm31c,eqep_epmm31c
      common /amwori/lt6,lt7,lt8,lt9,lt10,lt11
*     3632 (3632) + 108 (108) bytes used

      integer    lt6(0:3),lt7(0:3),lt8(0:3),lt9(0:3),lt10(0:4),
     &           lt11(0:5)
      real*8     evep_epmm11d(lepina),eqep_epmm11d(4,lepina),
     &           exep_epmm15m(2),ptep_epmm15m(4,3),
     &           evep_epmm31c(lepinv),eqep_epmm31c(4,lepinv)
      complex*16 cfep_epmm15m(2,4)
      complex*16 av6(lextrn*lextrn*lepina)
      complex*16 av7(lintrn*lextrn*lepina)
      complex*16 av8(lextrn*lintrn*lepinv)
      complex*16 av9(lextrn*lextrn*lepinv)
      complex*16 av10(lextrn*lextrn*lintrn*lextrn)
      complex*16 av11(lextrn*lepinv*lextrn*lextrn*lextrn)
      complex*16 atmp
*-----------------------------------------------------------------------

* Denominators of propagators
      aprop = 1.0d0
      call snprpd(pphase,aprop,vnep_epmm11,
     &      ama**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm15,
     &      amlp(1)**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm31,
     &      amz**2,amz*agz)

* Internal momenta
      call smintv(lepina,ama,pfep_epmm11,eqep_epmm11d,evep_epmm11d,
     &            vnep_epmm11,igauab)
      call smintf(amlp(1),pfep_epmm15,vnep_epmm15,exep_epmm15m,
     &            ptep_epmm15m,cfep_epmm15m)
      call smintv(lepinv,amz,pfep_epmm31,eqep_epmm31c,evep_epmm31c,
     &            vnep_epmm31,igauzb)

* Vertices (10)

*     6(0): pfep_epmm9 p
*     6(1): pfep_epmm2 p
*     6(2): pfep_epmm11 photon
      call stppv(2212,lextrn,amp,pfep_epmm2,+1,vnep_epmm2,exep_epmm2v,
     &            cfep_epmm2v,ptep_epmm2v,2212,lextrn,amp,pfep_epmm9,-1,
     &            vnep_epmm9,exep_epmm9v,cfep_epmm9v,ptep_epmm9v,22,
     &            lepina,ama,pfep_epmm11,+1,vnep_epmm11,exep_epmm11d,
     &            eqep_epmm11d,capp,lt6,av6)

*     7(0): pfep_epmm4 electron
*     7(1): pfep_epmm15 electron
*     7(2): pfep_epmm11 photon
      call smffv(lintrn,lextrn,lepina,exep_epmm15m,exep_epmm4m,amlp(1),
     &           amlp(1),call(1,1),cfep_epmm15m,cfep_epmm4m,
     &           ptep_epmm15m,ptep_epmm4m,eqep_epmm11d,lt7,av7)

*     8(0): pfep_epmm15 electron
*     8(1): pfep_epmm17 electron
*     8(2): pfep_epmm31 z
      call smffv(lextrn,lintrn,lepinv,exep_epmm17m,exep_epmm15m,amlp(1),
     &           amlp(1),czll(1,1),cfep_epmm17m,cfep_epmm15m,
     &           ptep_epmm17m,ptep_epmm15m,eqep_epmm31c,lt8,av8)

*     9(0): pfep_epmm33 muon
*     9(1): pfep_epmm62 muon
*     9(2): pfep_epmm31 z
      call smffv(lextrn,lextrn,lepinv,exep_epmm62n,exep_epmm33n,amlp(2),
     &           amlp(2),czll(1,2),cfep_epmm62n,cfep_epmm33n,
     &           ptep_epmm62n,ptep_epmm33n,eqep_epmm31c,lt9,av9)

      call smconv(lt6,lt7,3,3,evep_epmm11d,av6,av7,lt10,av10)
      call smconf(lt8,lt10,2,3,exep_epmm15m,av8,av10,lt11,av11)
      call smconv(lt9,lt11,3,2,evep_epmm31c,av9,av11,lt,av)

      sym = + 1.0d0
      cf  = + 1.0d0
      aprop         = cf*sym/aprop

      indexg(1) = 6
      indexg(2) = 5
      indexg(3) = 4
      indexg(4) = 1
      indexg(5) = 3
      indexg(6) = 2

      call aep_epmmmpord(lt, av, indexg, agcwrk)

      ancp(jgraph) = 0.0d0
*     nbase = 1
      do 500 ih = 0 , lep_epmmag-1
         atmp    = agcwrk(ih)*aprop
         agc(ih) = agc(ih) + atmp
         ancp(jgraph) = ancp(jgraph) + atmp*conjg(atmp)
  500 continue

      return
      end
************************************************************************
*             Graph No. 3 - 1
*         Generated No. 3
************************************************************************
      subroutine aep_epmmg3
      implicit real*8(a-h,o-z)

      include 'inclep_epmm1.h'
      include 'inclk.inc'
      include 'inclep_epmmp.h'
*-----------------------------------------------------------------------
      common /amwork/cfep_epmm26m,av6,av7,av8,av9,av10,av11,
     &               evep_epmm11d,eqep_epmm11d,exep_epmm26m,
     &               ptep_epmm26m,evep_epmm31d,eqep_epmm31d
      common /amwori/lt6,lt7,lt8,lt9,lt10,lt11
*     3632 (3632) + 108 (108) bytes used

      integer    lt6(0:3),lt7(0:3),lt8(0:3),lt9(0:3),lt10(0:4),
     &           lt11(0:5)
      real*8     evep_epmm11d(lepina),eqep_epmm11d(4,lepina),
     &           exep_epmm26m(2),ptep_epmm26m(4,3),
     &           evep_epmm31d(lepina),eqep_epmm31d(4,lepina)
      complex*16 cfep_epmm26m(2,4)
      complex*16 av6(lextrn*lextrn*lepina)
      complex*16 av7(lextrn*lintrn*lepina)
      complex*16 av8(lintrn*lextrn*lepina)
      complex*16 av9(lextrn*lextrn*lepina)
      complex*16 av10(lextrn*lextrn*lextrn*lintrn)
      complex*16 av11(lextrn*lepina*lextrn*lextrn*lextrn)
      complex*16 atmp
*-----------------------------------------------------------------------

* Denominators of propagators
      aprop = 1.0d0
      call snprpd(pphase,aprop,vnep_epmm11,
     &      ama**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm26,
     &      amlp(1)**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm31,
     &      ama**2,0.0d0)

* Internal momenta
      call smintv(lepina,ama,pfep_epmm11,eqep_epmm11d,evep_epmm11d,
     &            vnep_epmm11,igauab)
      call smintf(amlp(1),pfep_epmm26,vnep_epmm26,exep_epmm26m,
     &            ptep_epmm26m,cfep_epmm26m)
      call smintv(lepina,ama,pfep_epmm31,eqep_epmm31d,evep_epmm31d,
     &            vnep_epmm31,igauab)

* Vertices (10)

*     6(0): pfep_epmm9 p
*     6(1): pfep_epmm2 p
*     6(2): pfep_epmm11 photon
      call stppv(2212,lextrn,amp,pfep_epmm2,+1,vnep_epmm2,exep_epmm2v,
     &            cfep_epmm2v,ptep_epmm2v,2212,lextrn,amp,pfep_epmm9,-1,
     &            vnep_epmm9,exep_epmm9v,cfep_epmm9v,ptep_epmm9v,22,
     &            lepina,ama,pfep_epmm11,+1,vnep_epmm11,exep_epmm11d,
     &            eqep_epmm11d,capp,lt6,av6)

*     7(0): pfep_epmm26 electron
*     7(1): pfep_epmm17 electron
*     7(2): pfep_epmm11 photon
      call smffv(lextrn,lintrn,lepina,exep_epmm17m,exep_epmm26m,amlp(1),
     &           amlp(1),call(1,1),cfep_epmm17m,cfep_epmm26m,
     &           ptep_epmm17m,ptep_epmm26m,eqep_epmm11d,lt7,av7)

*     8(0): pfep_epmm4 electron
*     8(1): pfep_epmm26 electron
*     8(2): pfep_epmm31 photon
      call smffv(lintrn,lextrn,lepina,exep_epmm26m,exep_epmm4m,amlp(1),
     &           amlp(1),call(1,1),cfep_epmm26m,cfep_epmm4m,
     &           ptep_epmm26m,ptep_epmm4m,eqep_epmm31d,lt8,av8)

*     9(0): pfep_epmm33 muon
*     9(1): pfep_epmm62 muon
*     9(2): pfep_epmm31 photon
      call smffv(lextrn,lextrn,lepina,exep_epmm62n,exep_epmm33n,amlp(2),
     &           amlp(2),call(1,2),cfep_epmm62n,cfep_epmm33n,
     &           ptep_epmm62n,ptep_epmm33n,eqep_epmm31d,lt9,av9)

      call smconv(lt6,lt7,3,3,evep_epmm11d,av6,av7,lt10,av10)
      call smconf(lt8,lt10,1,4,exep_epmm26m,av8,av10,lt11,av11)
      call smconv(lt9,lt11,3,2,evep_epmm31d,av9,av11,lt,av)

      sym = + 1.0d0
      cf  = + 1.0d0
      aprop         = cf*sym/aprop

      indexg(1) = 6
      indexg(2) = 5
      indexg(3) = 2
      indexg(4) = 1
      indexg(5) = 3
      indexg(6) = 4

      call aep_epmmmpord(lt, av, indexg, agcwrk)

      ancp(jgraph) = 0.0d0
*     nbase = 1
      do 500 ih = 0 , lep_epmmag-1
         atmp    = agcwrk(ih)*aprop
         agc(ih) = agc(ih) + atmp
         ancp(jgraph) = ancp(jgraph) + atmp*conjg(atmp)
  500 continue

      return
      end
************************************************************************
*             Graph No. 4 - 1
*         Generated No. 4
************************************************************************
      subroutine aep_epmmg4
      implicit real*8(a-h,o-z)

      include 'inclep_epmm1.h'
      include 'inclk.inc'
      include 'inclep_epmmp.h'
*-----------------------------------------------------------------------
      common /amwork/cfep_epmm26m,av6,av7,av8,av9,av10,av11,
     &               evep_epmm11d,eqep_epmm11d,exep_epmm26m,
     &               ptep_epmm26m,evep_epmm31c,eqep_epmm31c
      common /amwori/lt6,lt7,lt8,lt9,lt10,lt11
*     3632 (3632) + 108 (108) bytes used

      integer    lt6(0:3),lt7(0:3),lt8(0:3),lt9(0:3),lt10(0:4),
     &           lt11(0:5)
      real*8     evep_epmm11d(lepina),eqep_epmm11d(4,lepina),
     &           exep_epmm26m(2),ptep_epmm26m(4,3),
     &           evep_epmm31c(lepinv),eqep_epmm31c(4,lepinv)
      complex*16 cfep_epmm26m(2,4)
      complex*16 av6(lextrn*lextrn*lepina)
      complex*16 av7(lextrn*lintrn*lepina)
      complex*16 av8(lintrn*lextrn*lepinv)
      complex*16 av9(lextrn*lextrn*lepinv)
      complex*16 av10(lextrn*lextrn*lextrn*lintrn)
      complex*16 av11(lextrn*lepinv*lextrn*lextrn*lextrn)
      complex*16 atmp
*-----------------------------------------------------------------------

* Denominators of propagators
      aprop = 1.0d0
      call snprpd(pphase,aprop,vnep_epmm11,
     &      ama**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm26,
     &      amlp(1)**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm31,
     &      amz**2,amz*agz)

* Internal momenta
      call smintv(lepina,ama,pfep_epmm11,eqep_epmm11d,evep_epmm11d,
     &            vnep_epmm11,igauab)
      call smintf(amlp(1),pfep_epmm26,vnep_epmm26,exep_epmm26m,
     &            ptep_epmm26m,cfep_epmm26m)
      call smintv(lepinv,amz,pfep_epmm31,eqep_epmm31c,evep_epmm31c,
     &            vnep_epmm31,igauzb)

* Vertices (10)

*     6(0): pfep_epmm9 p
*     6(1): pfep_epmm2 p
*     6(2): pfep_epmm11 photon
      call stppv(2212,lextrn,amp,pfep_epmm2,+1,vnep_epmm2,exep_epmm2v,
     &            cfep_epmm2v,ptep_epmm2v,2212,lextrn,amp,pfep_epmm9,-1,
     &            vnep_epmm9,exep_epmm9v,cfep_epmm9v,ptep_epmm9v,22,
     &            lepina,ama,pfep_epmm11,+1,vnep_epmm11,exep_epmm11d,
     &            eqep_epmm11d,capp,lt6,av6)

*     7(0): pfep_epmm26 electron
*     7(1): pfep_epmm17 electron
*     7(2): pfep_epmm11 photon
      call smffv(lextrn,lintrn,lepina,exep_epmm17m,exep_epmm26m,amlp(1),
     &           amlp(1),call(1,1),cfep_epmm17m,cfep_epmm26m,
     &           ptep_epmm17m,ptep_epmm26m,eqep_epmm11d,lt7,av7)

*     8(0): pfep_epmm4 electron
*     8(1): pfep_epmm26 electron
*     8(2): pfep_epmm31 z
      call smffv(lintrn,lextrn,lepinv,exep_epmm26m,exep_epmm4m,amlp(1),
     &           amlp(1),czll(1,1),cfep_epmm26m,cfep_epmm4m,
     &           ptep_epmm26m,ptep_epmm4m,eqep_epmm31c,lt8,av8)

*     9(0): pfep_epmm33 muon
*     9(1): pfep_epmm62 muon
*     9(2): pfep_epmm31 z
      call smffv(lextrn,lextrn,lepinv,exep_epmm62n,exep_epmm33n,amlp(2),
     &           amlp(2),czll(1,2),cfep_epmm62n,cfep_epmm33n,
     &           ptep_epmm62n,ptep_epmm33n,eqep_epmm31c,lt9,av9)

      call smconv(lt6,lt7,3,3,evep_epmm11d,av6,av7,lt10,av10)
      call smconf(lt8,lt10,1,4,exep_epmm26m,av8,av10,lt11,av11)
      call smconv(lt9,lt11,3,2,evep_epmm31c,av9,av11,lt,av)

      sym = + 1.0d0
      cf  = + 1.0d0
      aprop         = cf*sym/aprop

      indexg(1) = 6
      indexg(2) = 5
      indexg(3) = 2
      indexg(4) = 1
      indexg(5) = 3
      indexg(6) = 4

      call aep_epmmmpord(lt, av, indexg, agcwrk)

      ancp(jgraph) = 0.0d0
*     nbase = 1
      do 500 ih = 0 , lep_epmmag-1
         atmp    = agcwrk(ih)*aprop
         agc(ih) = agc(ih) + atmp
         ancp(jgraph) = ancp(jgraph) + atmp*conjg(atmp)
  500 continue

      return
      end
************************************************************************
*             Graph No. 5 - 1
*         Generated No. 5
************************************************************************
      subroutine aep_epmmg5
      implicit real*8(a-h,o-z)

      include 'inclep_epmm1.h'
      include 'inclk.inc'
      include 'inclep_epmmp.h'
*-----------------------------------------------------------------------
      common /amwork/cfep_epmm43n,av6,av7,av8,av9,av10,av11,
     &               evep_epmm11d,eqep_epmm11d,exep_epmm43n,
     &               ptep_epmm43n,evep_epmm20d,eqep_epmm20d
      common /amwori/lt6,lt7,lt8,lt9,lt10,lt11
*     3632 (3632) + 108 (108) bytes used

      integer    lt6(0:3),lt7(0:3),lt8(0:3),lt9(0:3),lt10(0:4),
     &           lt11(0:5)
      real*8     evep_epmm11d(lepina),eqep_epmm11d(4,lepina),
     &           exep_epmm43n(2),ptep_epmm43n(4,3),
     &           evep_epmm20d(lepina),eqep_epmm20d(4,lepina)
      complex*16 cfep_epmm43n(2,4)
      complex*16 av6(lextrn*lextrn*lepina)
      complex*16 av7(lintrn*lextrn*lepina)
      complex*16 av8(lextrn*lintrn*lepina)
      complex*16 av9(lextrn*lextrn*lepina)
      complex*16 av10(lextrn*lextrn*lintrn*lextrn)
      complex*16 av11(lextrn*lepina*lextrn*lextrn*lextrn)
      complex*16 atmp
*-----------------------------------------------------------------------

* Denominators of propagators
      aprop = 1.0d0
      call snprpd(pphase,aprop,vnep_epmm11,
     &      ama**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm43,
     &      amlp(2)**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm20,
     &      ama**2,0.0d0)

* Internal momenta
      call smintv(lepina,ama,pfep_epmm11,eqep_epmm11d,evep_epmm11d,
     &            vnep_epmm11,igauab)
      call smintf(amlp(2),pfep_epmm43,vnep_epmm43,exep_epmm43n,
     &            ptep_epmm43n,cfep_epmm43n)
      call smintv(lepina,ama,pfep_epmm20,eqep_epmm20d,evep_epmm20d,
     &            vnep_epmm20,igauab)

* Vertices (10)

*     6(0): pfep_epmm9 p
*     6(1): pfep_epmm2 p
*     6(2): pfep_epmm11 photon
      call stppv(2212,lextrn,amp,pfep_epmm2,+1,vnep_epmm2,exep_epmm2v,
     &            cfep_epmm2v,ptep_epmm2v,2212,lextrn,amp,pfep_epmm9,-1,
     &            vnep_epmm9,exep_epmm9v,cfep_epmm9v,ptep_epmm9v,22,
     &            lepina,ama,pfep_epmm11,+1,vnep_epmm11,exep_epmm11d,
     &            eqep_epmm11d,capp,lt6,av6)

*     7(0): pfep_epmm33 muon
*     7(1): pfep_epmm43 muon
*     7(2): pfep_epmm11 photon
      call smffv(lintrn,lextrn,lepina,exep_epmm43n,exep_epmm33n,amlp(2),
     &           amlp(2),call(1,2),cfep_epmm43n,cfep_epmm33n,
     &           ptep_epmm43n,ptep_epmm33n,eqep_epmm11d,lt7,av7)

*     8(0): pfep_epmm43 muon
*     8(1): pfep_epmm62 muon
*     8(2): pfep_epmm20 photon
      call smffv(lextrn,lintrn,lepina,exep_epmm62n,exep_epmm43n,amlp(2),
     &           amlp(2),call(1,2),cfep_epmm62n,cfep_epmm43n,
     &           ptep_epmm62n,ptep_epmm43n,eqep_epmm20d,lt8,av8)

*     9(0): pfep_epmm4 electron
*     9(1): pfep_epmm17 electron
*     9(2): pfep_epmm20 photon
      call smffv(lextrn,lextrn,lepina,exep_epmm17m,exep_epmm4m,amlp(1),
     &           amlp(1),call(1,1),cfep_epmm17m,cfep_epmm4m,
     &           ptep_epmm17m,ptep_epmm4m,eqep_epmm20d,lt9,av9)

      call smconv(lt6,lt7,3,3,evep_epmm11d,av6,av7,lt10,av10)
      call smconf(lt8,lt10,2,3,exep_epmm43n,av8,av10,lt11,av11)
      call smconv(lt9,lt11,3,2,evep_epmm20d,av9,av11,lt,av)

      sym = + 1.0d0
      cf  = + 1.0d0
      aprop         = cf*sym/aprop

      indexg(1) = 4
      indexg(2) = 2
      indexg(3) = 6
      indexg(4) = 1
      indexg(5) = 3
      indexg(6) = 5

      call aep_epmmmpord(lt, av, indexg, agcwrk)

      ancp(jgraph) = 0.0d0
*     nbase = 1
      do 500 ih = 0 , lep_epmmag-1
         atmp    = agcwrk(ih)*aprop
         agc(ih) = agc(ih) + atmp
         ancp(jgraph) = ancp(jgraph) + atmp*conjg(atmp)
  500 continue

      return
      end
************************************************************************
*             Graph No. 6 - 1
*         Generated No. 6
************************************************************************
      subroutine aep_epmmg6
      implicit real*8(a-h,o-z)

      include 'inclep_epmm1.h'
      include 'inclk.inc'
      include 'inclep_epmmp.h'
*-----------------------------------------------------------------------
      common /amwork/cfep_epmm43n,av6,av7,av8,av9,av10,av11,
     &               evep_epmm11d,eqep_epmm11d,exep_epmm43n,
     &               ptep_epmm43n,evep_epmm20c,eqep_epmm20c
      common /amwori/lt6,lt7,lt8,lt9,lt10,lt11
*     3632 (3632) + 108 (108) bytes used

      integer    lt6(0:3),lt7(0:3),lt8(0:3),lt9(0:3),lt10(0:4),
     &           lt11(0:5)
      real*8     evep_epmm11d(lepina),eqep_epmm11d(4,lepina),
     &           exep_epmm43n(2),ptep_epmm43n(4,3),
     &           evep_epmm20c(lepinv),eqep_epmm20c(4,lepinv)
      complex*16 cfep_epmm43n(2,4)
      complex*16 av6(lextrn*lextrn*lepina)
      complex*16 av7(lintrn*lextrn*lepina)
      complex*16 av8(lextrn*lintrn*lepinv)
      complex*16 av9(lextrn*lextrn*lepinv)
      complex*16 av10(lextrn*lextrn*lintrn*lextrn)
      complex*16 av11(lextrn*lepinv*lextrn*lextrn*lextrn)
      complex*16 atmp
*-----------------------------------------------------------------------

* Denominators of propagators
      aprop = 1.0d0
      call snprpd(pphase,aprop,vnep_epmm11,
     &      ama**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm43,
     &      amlp(2)**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm20,
     &      amz**2,amz*agz)

* Internal momenta
      call smintv(lepina,ama,pfep_epmm11,eqep_epmm11d,evep_epmm11d,
     &            vnep_epmm11,igauab)
      call smintf(amlp(2),pfep_epmm43,vnep_epmm43,exep_epmm43n,
     &            ptep_epmm43n,cfep_epmm43n)
      call smintv(lepinv,amz,pfep_epmm20,eqep_epmm20c,evep_epmm20c,
     &            vnep_epmm20,igauzb)

* Vertices (10)

*     6(0): pfep_epmm9 p
*     6(1): pfep_epmm2 p
*     6(2): pfep_epmm11 photon
      call stppv(2212,lextrn,amp,pfep_epmm2,+1,vnep_epmm2,exep_epmm2v,
     &            cfep_epmm2v,ptep_epmm2v,2212,lextrn,amp,pfep_epmm9,-1,
     &            vnep_epmm9,exep_epmm9v,cfep_epmm9v,ptep_epmm9v,22,
     &            lepina,ama,pfep_epmm11,+1,vnep_epmm11,exep_epmm11d,
     &            eqep_epmm11d,capp,lt6,av6)

*     7(0): pfep_epmm33 muon
*     7(1): pfep_epmm43 muon
*     7(2): pfep_epmm11 photon
      call smffv(lintrn,lextrn,lepina,exep_epmm43n,exep_epmm33n,amlp(2),
     &           amlp(2),call(1,2),cfep_epmm43n,cfep_epmm33n,
     &           ptep_epmm43n,ptep_epmm33n,eqep_epmm11d,lt7,av7)

*     8(0): pfep_epmm43 muon
*     8(1): pfep_epmm62 muon
*     8(2): pfep_epmm20 z
      call smffv(lextrn,lintrn,lepinv,exep_epmm62n,exep_epmm43n,amlp(2),
     &           amlp(2),czll(1,2),cfep_epmm62n,cfep_epmm43n,
     &           ptep_epmm62n,ptep_epmm43n,eqep_epmm20c,lt8,av8)

*     9(0): pfep_epmm4 electron
*     9(1): pfep_epmm17 electron
*     9(2): pfep_epmm20 z
      call smffv(lextrn,lextrn,lepinv,exep_epmm17m,exep_epmm4m,amlp(1),
     &           amlp(1),czll(1,1),cfep_epmm17m,cfep_epmm4m,
     &           ptep_epmm17m,ptep_epmm4m,eqep_epmm20c,lt9,av9)

      call smconv(lt6,lt7,3,3,evep_epmm11d,av6,av7,lt10,av10)
      call smconf(lt8,lt10,2,3,exep_epmm43n,av8,av10,lt11,av11)
      call smconv(lt9,lt11,3,2,evep_epmm20c,av9,av11,lt,av)

      sym = + 1.0d0
      cf  = + 1.0d0
      aprop         = cf*sym/aprop

      indexg(1) = 4
      indexg(2) = 2
      indexg(3) = 6
      indexg(4) = 1
      indexg(5) = 3
      indexg(6) = 5

      call aep_epmmmpord(lt, av, indexg, agcwrk)

      ancp(jgraph) = 0.0d0
*     nbase = 1
      do 500 ih = 0 , lep_epmmag-1
         atmp    = agcwrk(ih)*aprop
         agc(ih) = agc(ih) + atmp
         ancp(jgraph) = ancp(jgraph) + atmp*conjg(atmp)
  500 continue

      return
      end
************************************************************************
*             Graph No. 7 - 1
*         Generated No. 7
************************************************************************
      subroutine aep_epmmg7
      implicit real*8(a-h,o-z)

      include 'inclep_epmm1.h'
      include 'inclk.inc'
      include 'inclep_epmmp.h'
*-----------------------------------------------------------------------
      common /amwork/cfep_epmm53n,av6,av7,av8,av9,av10,av11,
     &               evep_epmm11d,eqep_epmm11d,exep_epmm53n,
     &               ptep_epmm53n,evep_epmm20d,eqep_epmm20d
      common /amwori/lt6,lt7,lt8,lt9,lt10,lt11
*     3632 (3632) + 108 (108) bytes used

      integer    lt6(0:3),lt7(0:3),lt8(0:3),lt9(0:3),lt10(0:4),
     &           lt11(0:5)
      real*8     evep_epmm11d(lepina),eqep_epmm11d(4,lepina),
     &           exep_epmm53n(2),ptep_epmm53n(4,3),
     &           evep_epmm20d(lepina),eqep_epmm20d(4,lepina)
      complex*16 cfep_epmm53n(2,4)
      complex*16 av6(lextrn*lextrn*lepina)
      complex*16 av7(lextrn*lintrn*lepina)
      complex*16 av8(lintrn*lextrn*lepina)
      complex*16 av9(lextrn*lextrn*lepina)
      complex*16 av10(lextrn*lextrn*lextrn*lintrn)
      complex*16 av11(lextrn*lepina*lextrn*lextrn*lextrn)
      complex*16 atmp
*-----------------------------------------------------------------------

* Denominators of propagators
      aprop = 1.0d0
      call snprpd(pphase,aprop,vnep_epmm11,
     &      ama**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm53,
     &      amlp(2)**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm20,
     &      ama**2,0.0d0)

* Internal momenta
      call smintv(lepina,ama,pfep_epmm11,eqep_epmm11d,evep_epmm11d,
     &            vnep_epmm11,igauab)
      call smintf(amlp(2),pfep_epmm53,vnep_epmm53,exep_epmm53n,
     &            ptep_epmm53n,cfep_epmm53n)
      call smintv(lepina,ama,pfep_epmm20,eqep_epmm20d,evep_epmm20d,
     &            vnep_epmm20,igauab)

* Vertices (10)

*     6(0): pfep_epmm9 p
*     6(1): pfep_epmm2 p
*     6(2): pfep_epmm11 photon
      call stppv(2212,lextrn,amp,pfep_epmm2,+1,vnep_epmm2,exep_epmm2v,
     &            cfep_epmm2v,ptep_epmm2v,2212,lextrn,amp,pfep_epmm9,-1,
     &            vnep_epmm9,exep_epmm9v,cfep_epmm9v,ptep_epmm9v,22,
     &            lepina,ama,pfep_epmm11,+1,vnep_epmm11,exep_epmm11d,
     &            eqep_epmm11d,capp,lt6,av6)

*     7(0): pfep_epmm53 muon
*     7(1): pfep_epmm62 muon
*     7(2): pfep_epmm11 photon
      call smffv(lextrn,lintrn,lepina,exep_epmm62n,exep_epmm53n,amlp(2),
     &           amlp(2),call(1,2),cfep_epmm62n,cfep_epmm53n,
     &           ptep_epmm62n,ptep_epmm53n,eqep_epmm11d,lt7,av7)

*     8(0): pfep_epmm33 muon
*     8(1): pfep_epmm53 muon
*     8(2): pfep_epmm20 photon
      call smffv(lintrn,lextrn,lepina,exep_epmm53n,exep_epmm33n,amlp(2),
     &           amlp(2),call(1,2),cfep_epmm53n,cfep_epmm33n,
     &           ptep_epmm53n,ptep_epmm33n,eqep_epmm20d,lt8,av8)

*     9(0): pfep_epmm4 electron
*     9(1): pfep_epmm17 electron
*     9(2): pfep_epmm20 photon
      call smffv(lextrn,lextrn,lepina,exep_epmm17m,exep_epmm4m,amlp(1),
     &           amlp(1),call(1,1),cfep_epmm17m,cfep_epmm4m,
     &           ptep_epmm17m,ptep_epmm4m,eqep_epmm20d,lt9,av9)

      call smconv(lt6,lt7,3,3,evep_epmm11d,av6,av7,lt10,av10)
      call smconf(lt8,lt10,1,4,exep_epmm53n,av8,av10,lt11,av11)
      call smconv(lt9,lt11,3,2,evep_epmm20d,av9,av11,lt,av)

      sym = + 1.0d0
      cf  = + 1.0d0
      aprop         = cf*sym/aprop

      indexg(1) = 4
      indexg(2) = 2
      indexg(3) = 5
      indexg(4) = 1
      indexg(5) = 3
      indexg(6) = 6

      call aep_epmmmpord(lt, av, indexg, agcwrk)

      ancp(jgraph) = 0.0d0
*     nbase = 1
      do 500 ih = 0 , lep_epmmag-1
         atmp    = agcwrk(ih)*aprop
         agc(ih) = agc(ih) + atmp
         ancp(jgraph) = ancp(jgraph) + atmp*conjg(atmp)
  500 continue

      return
      end
************************************************************************
*             Graph No. 8 - 1
*         Generated No. 8
************************************************************************
      subroutine aep_epmmg8
      implicit real*8(a-h,o-z)

      include 'inclep_epmm1.h'
      include 'inclk.inc'
      include 'inclep_epmmp.h'
*-----------------------------------------------------------------------
      common /amwork/cfep_epmm53n,av6,av7,av8,av9,av10,av11,
     &               evep_epmm11d,eqep_epmm11d,exep_epmm53n,
     &               ptep_epmm53n,evep_epmm20c,eqep_epmm20c
      common /amwori/lt6,lt7,lt8,lt9,lt10,lt11
*     3632 (3632) + 108 (108) bytes used

      integer    lt6(0:3),lt7(0:3),lt8(0:3),lt9(0:3),lt10(0:4),
     &           lt11(0:5)
      real*8     evep_epmm11d(lepina),eqep_epmm11d(4,lepina),
     &           exep_epmm53n(2),ptep_epmm53n(4,3),
     &           evep_epmm20c(lepinv),eqep_epmm20c(4,lepinv)
      complex*16 cfep_epmm53n(2,4)
      complex*16 av6(lextrn*lextrn*lepina)
      complex*16 av7(lextrn*lintrn*lepina)
      complex*16 av8(lintrn*lextrn*lepinv)
      complex*16 av9(lextrn*lextrn*lepinv)
      complex*16 av10(lextrn*lextrn*lextrn*lintrn)
      complex*16 av11(lextrn*lepinv*lextrn*lextrn*lextrn)
      complex*16 atmp
*-----------------------------------------------------------------------

* Denominators of propagators
      aprop = 1.0d0
      call snprpd(pphase,aprop,vnep_epmm11,
     &      ama**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm53,
     &      amlp(2)**2,0.0d0)
      call snprpd(pphase,aprop,vnep_epmm20,
     &      amz**2,amz*agz)

* Internal momenta
      call smintv(lepina,ama,pfep_epmm11,eqep_epmm11d,evep_epmm11d,
     &            vnep_epmm11,igauab)
      call smintf(amlp(2),pfep_epmm53,vnep_epmm53,exep_epmm53n,
     &            ptep_epmm53n,cfep_epmm53n)
      call smintv(lepinv,amz,pfep_epmm20,eqep_epmm20c,evep_epmm20c,
     &            vnep_epmm20,igauzb)

* Vertices (10)

*     6(0): pfep_epmm9 p
*     6(1): pfep_epmm2 p
*     6(2): pfep_epmm11 photon
      call stppv(2212,lextrn,amp,pfep_epmm2,+1,vnep_epmm2,exep_epmm2v,
     &            cfep_epmm2v,ptep_epmm2v,2212,lextrn,amp,pfep_epmm9,-1,
     &            vnep_epmm9,exep_epmm9v,cfep_epmm9v,ptep_epmm9v,22,
     &            lepina,ama,pfep_epmm11,+1,vnep_epmm11,exep_epmm11d,
     &            eqep_epmm11d,capp,lt6,av6)

*     7(0): pfep_epmm53 muon
*     7(1): pfep_epmm62 muon
*     7(2): pfep_epmm11 photon
      call smffv(lextrn,lintrn,lepina,exep_epmm62n,exep_epmm53n,amlp(2),
     &           amlp(2),call(1,2),cfep_epmm62n,cfep_epmm53n,
     &           ptep_epmm62n,ptep_epmm53n,eqep_epmm11d,lt7,av7)

*     8(0): pfep_epmm33 muon
*     8(1): pfep_epmm53 muon
*     8(2): pfep_epmm20 z
      call smffv(lintrn,lextrn,lepinv,exep_epmm53n,exep_epmm33n,amlp(2),
     &           amlp(2),czll(1,2),cfep_epmm53n,cfep_epmm33n,
     &           ptep_epmm53n,ptep_epmm33n,eqep_epmm20c,lt8,av8)

*     9(0): pfep_epmm4 electron
*     9(1): pfep_epmm17 electron
*     9(2): pfep_epmm20 z
      call smffv(lextrn,lextrn,lepinv,exep_epmm17m,exep_epmm4m,amlp(1),
     &           amlp(1),czll(1,1),cfep_epmm17m,cfep_epmm4m,
     &           ptep_epmm17m,ptep_epmm4m,eqep_epmm20c,lt9,av9)

      call smconv(lt6,lt7,3,3,evep_epmm11d,av6,av7,lt10,av10)
      call smconf(lt8,lt10,1,4,exep_epmm53n,av8,av10,lt11,av11)
      call smconv(lt9,lt11,3,2,evep_epmm20c,av9,av11,lt,av)

      sym = + 1.0d0
      cf  = + 1.0d0
      aprop         = cf*sym/aprop

      indexg(1) = 4
      indexg(2) = 2
      indexg(3) = 5
      indexg(4) = 1
      indexg(5) = 3
      indexg(6) = 6

      call aep_epmmmpord(lt, av, indexg, agcwrk)

      ancp(jgraph) = 0.0d0
*     nbase = 1
      do 500 ih = 0 , lep_epmmag-1
         atmp    = agcwrk(ih)*aprop
         agc(ih) = agc(ih) + atmp
         ancp(jgraph) = ancp(jgraph) + atmp*conjg(atmp)
  500 continue

      return
      end