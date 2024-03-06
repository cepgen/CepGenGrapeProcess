*     File ec_ecmm/aec_ecmmkptbl.f : Sat Mar 18 19:45:03 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aec_ecmmkptbl
      implicit real*8(a-h,o-z)

      include 'inclec_ecmm1.h'
      include 'inclk.inc'
      include 'inclec_ecmmp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfec_ecmm2(i) =  + peec_ecmm1(i)
        pfec_ecmm4(i) =  + peec_ecmm2(i)
        pfec_ecmm6(i) =  + peec_ecmm1(i) + peec_ecmm2(i)
        pfec_ecmm8(i) =  - peec_ecmm3(i)
        pfec_ecmm10(i) =  + peec_ecmm1(i) - peec_ecmm3(i)
        pfec_ecmm12(i) =  + peec_ecmm2(i) - peec_ecmm3(i)
        pfec_ecmm14(i) =  + peec_ecmm1(i) + peec_ecmm2(i) - peec_ecmm3(i)
        pfec_ecmm16(i) =  - peec_ecmm4(i)
        pfec_ecmm18(i) =  + peec_ecmm1(i) - peec_ecmm4(i)
        pfec_ecmm20(i) =  + peec_ecmm2(i) - peec_ecmm4(i)
        pfec_ecmm22(i) =  + peec_ecmm1(i) + peec_ecmm2(i) - peec_ecmm4(i)
        pfec_ecmm24(i) =  - peec_ecmm3(i) - peec_ecmm4(i)
        pfec_ecmm26(i) =  + peec_ecmm1(i) - peec_ecmm3(i) - peec_ecmm4(i)
        pfec_ecmm28(i) =  + peec_ecmm2(i) - peec_ecmm3(i) - peec_ecmm4(i)
        pfec_ecmm30(i) =  + peec_ecmm5(i) + peec_ecmm6(i)
        pfec_ecmm32(i) =  - peec_ecmm5(i)
        pfec_ecmm34(i) =  + peec_ecmm1(i) - peec_ecmm5(i)
        pfec_ecmm36(i) =  + peec_ecmm2(i) - peec_ecmm5(i)
        pfec_ecmm38(i) =  + peec_ecmm1(i) + peec_ecmm2(i) - peec_ecmm5(i)
        pfec_ecmm40(i) =  - peec_ecmm3(i) - peec_ecmm5(i)
        pfec_ecmm42(i) =  + peec_ecmm1(i) - peec_ecmm3(i) - peec_ecmm5(i)
        pfec_ecmm44(i) =  + peec_ecmm2(i) - peec_ecmm3(i) - peec_ecmm5(i)
        pfec_ecmm46(i) =  + peec_ecmm4(i) + peec_ecmm6(i)
        pfec_ecmm48(i) =  - peec_ecmm4(i) - peec_ecmm5(i)
        pfec_ecmm50(i) =  + peec_ecmm1(i) - peec_ecmm4(i) - peec_ecmm5(i)
        pfec_ecmm52(i) =  + peec_ecmm2(i) - peec_ecmm4(i) - peec_ecmm5(i)
        pfec_ecmm54(i) =  + peec_ecmm3(i) + peec_ecmm6(i)
        pfec_ecmm56(i) =  - peec_ecmm3(i) - peec_ecmm4(i) - peec_ecmm5(i)
        pfec_ecmm58(i) =  - peec_ecmm2(i) + peec_ecmm6(i)
        pfec_ecmm60(i) =  - peec_ecmm1(i) + peec_ecmm6(i)
        pfec_ecmm62(i) =  + peec_ecmm6(i)

        pfec_ecmm3(i) = - pfec_ecmm2(i)
        pfec_ecmm5(i) = - pfec_ecmm4(i)
        pfec_ecmm7(i) = - pfec_ecmm6(i)
        pfec_ecmm9(i) = - pfec_ecmm8(i)
        pfec_ecmm11(i) = - pfec_ecmm10(i)
        pfec_ecmm13(i) = - pfec_ecmm12(i)
        pfec_ecmm15(i) = - pfec_ecmm14(i)
        pfec_ecmm17(i) = - pfec_ecmm16(i)
        pfec_ecmm19(i) = - pfec_ecmm18(i)
        pfec_ecmm21(i) = - pfec_ecmm20(i)
        pfec_ecmm23(i) = - pfec_ecmm22(i)
        pfec_ecmm25(i) = - pfec_ecmm24(i)
        pfec_ecmm27(i) = - pfec_ecmm26(i)
        pfec_ecmm29(i) = - pfec_ecmm28(i)
        pfec_ecmm31(i) = - pfec_ecmm30(i)
        pfec_ecmm33(i) = - pfec_ecmm32(i)
        pfec_ecmm35(i) = - pfec_ecmm34(i)
        pfec_ecmm37(i) = - pfec_ecmm36(i)
        pfec_ecmm39(i) = - pfec_ecmm38(i)
        pfec_ecmm41(i) = - pfec_ecmm40(i)
        pfec_ecmm43(i) = - pfec_ecmm42(i)
        pfec_ecmm45(i) = - pfec_ecmm44(i)
        pfec_ecmm47(i) = - pfec_ecmm46(i)
        pfec_ecmm49(i) = - pfec_ecmm48(i)
        pfec_ecmm51(i) = - pfec_ecmm50(i)
        pfec_ecmm53(i) = - pfec_ecmm52(i)
        pfec_ecmm55(i) = - pfec_ecmm54(i)
        pfec_ecmm57(i) = - pfec_ecmm56(i)
        pfec_ecmm59(i) = - pfec_ecmm58(i)
        pfec_ecmm61(i) = - pfec_ecmm60(i)
        pfec_ecmm63(i) = - pfec_ecmm62(i)
  100 continue

      vnec_ecmm2 =  + amass2(1)
      vnec_ecmm4 =  + amass2(2)
      vnec_ecmm6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vnec_ecmm8 =  + amass2(3)
      vnec_ecmm10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vnec_ecmm12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vnec_ecmm14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vnec_ecmm16 =  + amass2(4)
      vnec_ecmm18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vnec_ecmm20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vnec_ecmm22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vnec_ecmm24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vnec_ecmm26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vnec_ecmm28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vnec_ecmm30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vnec_ecmm32 =  + amass2(5)
      vnec_ecmm34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vnec_ecmm36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vnec_ecmm38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vnec_ecmm40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vnec_ecmm42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vnec_ecmm44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vnec_ecmm46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vnec_ecmm48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vnec_ecmm50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vnec_ecmm52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vnec_ecmm54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vnec_ecmm56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vnec_ecmm58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vnec_ecmm60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vnec_ecmm62 =  + amass2(6)
      vnec_ecmm3 = vnec_ecmm2
      vnec_ecmm5 = vnec_ecmm4
      vnec_ecmm7 = vnec_ecmm6
      vnec_ecmm9 = vnec_ecmm8
      vnec_ecmm11 = vnec_ecmm10
      vnec_ecmm13 = vnec_ecmm12
      vnec_ecmm15 = vnec_ecmm14
      vnec_ecmm17 = vnec_ecmm16
      vnec_ecmm19 = vnec_ecmm18
      vnec_ecmm21 = vnec_ecmm20
      vnec_ecmm23 = vnec_ecmm22
      vnec_ecmm25 = vnec_ecmm24
      vnec_ecmm27 = vnec_ecmm26
      vnec_ecmm29 = vnec_ecmm28
      vnec_ecmm31 = vnec_ecmm30
      vnec_ecmm33 = vnec_ecmm32
      vnec_ecmm35 = vnec_ecmm34
      vnec_ecmm37 = vnec_ecmm36
      vnec_ecmm39 = vnec_ecmm38
      vnec_ecmm41 = vnec_ecmm40
      vnec_ecmm43 = vnec_ecmm42
      vnec_ecmm45 = vnec_ecmm44
      vnec_ecmm47 = vnec_ecmm46
      vnec_ecmm49 = vnec_ecmm48
      vnec_ecmm51 = vnec_ecmm50
      vnec_ecmm53 = vnec_ecmm52
      vnec_ecmm55 = vnec_ecmm54
      vnec_ecmm57 = vnec_ecmm56
      vnec_ecmm59 = vnec_ecmm58
      vnec_ecmm61 = vnec_ecmm60
      vnec_ecmm63 = vnec_ecmm62
      return
      end
