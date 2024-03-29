*     File ebb_ebbmm/aebb_ebbmmkptbl.f : Sat Mar 18 19:45:06 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aebb_ebbmmkptbl
      implicit real*8(a-h,o-z)

      include 'inclebb_ebbmm1.h'
      include 'inclk.inc'
      include 'inclebb_ebbmmp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfebb_ebbmm2(i) =  + peebb_ebbmm1(i)
        pfebb_ebbmm4(i) =  + peebb_ebbmm2(i)
        pfebb_ebbmm6(i) =  + peebb_ebbmm1(i) + peebb_ebbmm2(i)
        pfebb_ebbmm8(i) =  - peebb_ebbmm3(i)
        pfebb_ebbmm10(i) =  + peebb_ebbmm1(i) - peebb_ebbmm3(i)
        pfebb_ebbmm12(i) =  + peebb_ebbmm2(i) - peebb_ebbmm3(i)
        pfebb_ebbmm14(i) =  + peebb_ebbmm1(i) + peebb_ebbmm2(i) - peebb_ebbmm3(i)
        pfebb_ebbmm16(i) =  - peebb_ebbmm4(i)
        pfebb_ebbmm18(i) =  + peebb_ebbmm1(i) - peebb_ebbmm4(i)
        pfebb_ebbmm20(i) =  + peebb_ebbmm2(i) - peebb_ebbmm4(i)
        pfebb_ebbmm22(i) =  + peebb_ebbmm1(i) + peebb_ebbmm2(i) - peebb_ebbmm4(i)
        pfebb_ebbmm24(i) =  - peebb_ebbmm3(i) - peebb_ebbmm4(i)
        pfebb_ebbmm26(i) =  + peebb_ebbmm1(i) - peebb_ebbmm3(i) - peebb_ebbmm4(i)
        pfebb_ebbmm28(i) =  + peebb_ebbmm2(i) - peebb_ebbmm3(i) - peebb_ebbmm4(i)
        pfebb_ebbmm30(i) =  + peebb_ebbmm5(i) + peebb_ebbmm6(i)
        pfebb_ebbmm32(i) =  - peebb_ebbmm5(i)
        pfebb_ebbmm34(i) =  + peebb_ebbmm1(i) - peebb_ebbmm5(i)
        pfebb_ebbmm36(i) =  + peebb_ebbmm2(i) - peebb_ebbmm5(i)
        pfebb_ebbmm38(i) =  + peebb_ebbmm1(i) + peebb_ebbmm2(i) - peebb_ebbmm5(i)
        pfebb_ebbmm40(i) =  - peebb_ebbmm3(i) - peebb_ebbmm5(i)
        pfebb_ebbmm42(i) =  + peebb_ebbmm1(i) - peebb_ebbmm3(i) - peebb_ebbmm5(i)
        pfebb_ebbmm44(i) =  + peebb_ebbmm2(i) - peebb_ebbmm3(i) - peebb_ebbmm5(i)
        pfebb_ebbmm46(i) =  + peebb_ebbmm4(i) + peebb_ebbmm6(i)
        pfebb_ebbmm48(i) =  - peebb_ebbmm4(i) - peebb_ebbmm5(i)
        pfebb_ebbmm50(i) =  + peebb_ebbmm1(i) - peebb_ebbmm4(i) - peebb_ebbmm5(i)
        pfebb_ebbmm52(i) =  + peebb_ebbmm2(i) - peebb_ebbmm4(i) - peebb_ebbmm5(i)
        pfebb_ebbmm54(i) =  + peebb_ebbmm3(i) + peebb_ebbmm6(i)
        pfebb_ebbmm56(i) =  - peebb_ebbmm3(i) - peebb_ebbmm4(i) - peebb_ebbmm5(i)
        pfebb_ebbmm58(i) =  - peebb_ebbmm2(i) + peebb_ebbmm6(i)
        pfebb_ebbmm60(i) =  - peebb_ebbmm1(i) + peebb_ebbmm6(i)
        pfebb_ebbmm62(i) =  + peebb_ebbmm6(i)

        pfebb_ebbmm3(i) = - pfebb_ebbmm2(i)
        pfebb_ebbmm5(i) = - pfebb_ebbmm4(i)
        pfebb_ebbmm7(i) = - pfebb_ebbmm6(i)
        pfebb_ebbmm9(i) = - pfebb_ebbmm8(i)
        pfebb_ebbmm11(i) = - pfebb_ebbmm10(i)
        pfebb_ebbmm13(i) = - pfebb_ebbmm12(i)
        pfebb_ebbmm15(i) = - pfebb_ebbmm14(i)
        pfebb_ebbmm17(i) = - pfebb_ebbmm16(i)
        pfebb_ebbmm19(i) = - pfebb_ebbmm18(i)
        pfebb_ebbmm21(i) = - pfebb_ebbmm20(i)
        pfebb_ebbmm23(i) = - pfebb_ebbmm22(i)
        pfebb_ebbmm25(i) = - pfebb_ebbmm24(i)
        pfebb_ebbmm27(i) = - pfebb_ebbmm26(i)
        pfebb_ebbmm29(i) = - pfebb_ebbmm28(i)
        pfebb_ebbmm31(i) = - pfebb_ebbmm30(i)
        pfebb_ebbmm33(i) = - pfebb_ebbmm32(i)
        pfebb_ebbmm35(i) = - pfebb_ebbmm34(i)
        pfebb_ebbmm37(i) = - pfebb_ebbmm36(i)
        pfebb_ebbmm39(i) = - pfebb_ebbmm38(i)
        pfebb_ebbmm41(i) = - pfebb_ebbmm40(i)
        pfebb_ebbmm43(i) = - pfebb_ebbmm42(i)
        pfebb_ebbmm45(i) = - pfebb_ebbmm44(i)
        pfebb_ebbmm47(i) = - pfebb_ebbmm46(i)
        pfebb_ebbmm49(i) = - pfebb_ebbmm48(i)
        pfebb_ebbmm51(i) = - pfebb_ebbmm50(i)
        pfebb_ebbmm53(i) = - pfebb_ebbmm52(i)
        pfebb_ebbmm55(i) = - pfebb_ebbmm54(i)
        pfebb_ebbmm57(i) = - pfebb_ebbmm56(i)
        pfebb_ebbmm59(i) = - pfebb_ebbmm58(i)
        pfebb_ebbmm61(i) = - pfebb_ebbmm60(i)
        pfebb_ebbmm63(i) = - pfebb_ebbmm62(i)
  100 continue

      vnebb_ebbmm2 =  + amass2(1)
      vnebb_ebbmm4 =  + amass2(2)
      vnebb_ebbmm6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vnebb_ebbmm8 =  + amass2(3)
      vnebb_ebbmm10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vnebb_ebbmm12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vnebb_ebbmm14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vnebb_ebbmm16 =  + amass2(4)
      vnebb_ebbmm18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vnebb_ebbmm20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vnebb_ebbmm22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vnebb_ebbmm24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vnebb_ebbmm26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vnebb_ebbmm28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vnebb_ebbmm30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vnebb_ebbmm32 =  + amass2(5)
      vnebb_ebbmm34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vnebb_ebbmm36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vnebb_ebbmm38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vnebb_ebbmm40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vnebb_ebbmm42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vnebb_ebbmm44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vnebb_ebbmm46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vnebb_ebbmm48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vnebb_ebbmm50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vnebb_ebbmm52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vnebb_ebbmm54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vnebb_ebbmm56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vnebb_ebbmm58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vnebb_ebbmm60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vnebb_ebbmm62 =  + amass2(6)
      vnebb_ebbmm3 = vnebb_ebbmm2
      vnebb_ebbmm5 = vnebb_ebbmm4
      vnebb_ebbmm7 = vnebb_ebbmm6
      vnebb_ebbmm9 = vnebb_ebbmm8
      vnebb_ebbmm11 = vnebb_ebbmm10
      vnebb_ebbmm13 = vnebb_ebbmm12
      vnebb_ebbmm15 = vnebb_ebbmm14
      vnebb_ebbmm17 = vnebb_ebbmm16
      vnebb_ebbmm19 = vnebb_ebbmm18
      vnebb_ebbmm21 = vnebb_ebbmm20
      vnebb_ebbmm23 = vnebb_ebbmm22
      vnebb_ebbmm25 = vnebb_ebbmm24
      vnebb_ebbmm27 = vnebb_ebbmm26
      vnebb_ebbmm29 = vnebb_ebbmm28
      vnebb_ebbmm31 = vnebb_ebbmm30
      vnebb_ebbmm33 = vnebb_ebbmm32
      vnebb_ebbmm35 = vnebb_ebbmm34
      vnebb_ebbmm37 = vnebb_ebbmm36
      vnebb_ebbmm39 = vnebb_ebbmm38
      vnebb_ebbmm41 = vnebb_ebbmm40
      vnebb_ebbmm43 = vnebb_ebbmm42
      vnebb_ebbmm45 = vnebb_ebbmm44
      vnebb_ebbmm47 = vnebb_ebbmm46
      vnebb_ebbmm49 = vnebb_ebbmm48
      vnebb_ebbmm51 = vnebb_ebbmm50
      vnebb_ebbmm53 = vnebb_ebbmm52
      vnebb_ebbmm55 = vnebb_ebbmm54
      vnebb_ebbmm57 = vnebb_ebbmm56
      vnebb_ebbmm59 = vnebb_ebbmm58
      vnebb_ebbmm61 = vnebb_ebbmm60
      vnebb_ebbmm63 = vnebb_ebbmm62
      return
      end
