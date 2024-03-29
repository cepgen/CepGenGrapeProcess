*     File ed_edmm/aed_edmmkptbl.f : Sat Mar 18 19:45:00 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aed_edmmkptbl
      implicit real*8(a-h,o-z)

      include 'incled_edmm1.h'
      include 'inclk.inc'
      include 'incled_edmmp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfed_edmm2(i) =  + peed_edmm1(i)
        pfed_edmm4(i) =  + peed_edmm2(i)
        pfed_edmm6(i) =  + peed_edmm1(i) + peed_edmm2(i)
        pfed_edmm8(i) =  - peed_edmm3(i)
        pfed_edmm10(i) =  + peed_edmm1(i) - peed_edmm3(i)
        pfed_edmm12(i) =  + peed_edmm2(i) - peed_edmm3(i)
        pfed_edmm14(i) =  + peed_edmm1(i) + peed_edmm2(i) - peed_edmm3(i)
        pfed_edmm16(i) =  - peed_edmm4(i)
        pfed_edmm18(i) =  + peed_edmm1(i) - peed_edmm4(i)
        pfed_edmm20(i) =  + peed_edmm2(i) - peed_edmm4(i)
        pfed_edmm22(i) =  + peed_edmm1(i) + peed_edmm2(i) - peed_edmm4(i)
        pfed_edmm24(i) =  - peed_edmm3(i) - peed_edmm4(i)
        pfed_edmm26(i) =  + peed_edmm1(i) - peed_edmm3(i) - peed_edmm4(i)
        pfed_edmm28(i) =  + peed_edmm2(i) - peed_edmm3(i) - peed_edmm4(i)
        pfed_edmm30(i) =  + peed_edmm5(i) + peed_edmm6(i)
        pfed_edmm32(i) =  - peed_edmm5(i)
        pfed_edmm34(i) =  + peed_edmm1(i) - peed_edmm5(i)
        pfed_edmm36(i) =  + peed_edmm2(i) - peed_edmm5(i)
        pfed_edmm38(i) =  + peed_edmm1(i) + peed_edmm2(i) - peed_edmm5(i)
        pfed_edmm40(i) =  - peed_edmm3(i) - peed_edmm5(i)
        pfed_edmm42(i) =  + peed_edmm1(i) - peed_edmm3(i) - peed_edmm5(i)
        pfed_edmm44(i) =  + peed_edmm2(i) - peed_edmm3(i) - peed_edmm5(i)
        pfed_edmm46(i) =  + peed_edmm4(i) + peed_edmm6(i)
        pfed_edmm48(i) =  - peed_edmm4(i) - peed_edmm5(i)
        pfed_edmm50(i) =  + peed_edmm1(i) - peed_edmm4(i) - peed_edmm5(i)
        pfed_edmm52(i) =  + peed_edmm2(i) - peed_edmm4(i) - peed_edmm5(i)
        pfed_edmm54(i) =  + peed_edmm3(i) + peed_edmm6(i)
        pfed_edmm56(i) =  - peed_edmm3(i) - peed_edmm4(i) - peed_edmm5(i)
        pfed_edmm58(i) =  - peed_edmm2(i) + peed_edmm6(i)
        pfed_edmm60(i) =  - peed_edmm1(i) + peed_edmm6(i)
        pfed_edmm62(i) =  + peed_edmm6(i)

        pfed_edmm3(i) = - pfed_edmm2(i)
        pfed_edmm5(i) = - pfed_edmm4(i)
        pfed_edmm7(i) = - pfed_edmm6(i)
        pfed_edmm9(i) = - pfed_edmm8(i)
        pfed_edmm11(i) = - pfed_edmm10(i)
        pfed_edmm13(i) = - pfed_edmm12(i)
        pfed_edmm15(i) = - pfed_edmm14(i)
        pfed_edmm17(i) = - pfed_edmm16(i)
        pfed_edmm19(i) = - pfed_edmm18(i)
        pfed_edmm21(i) = - pfed_edmm20(i)
        pfed_edmm23(i) = - pfed_edmm22(i)
        pfed_edmm25(i) = - pfed_edmm24(i)
        pfed_edmm27(i) = - pfed_edmm26(i)
        pfed_edmm29(i) = - pfed_edmm28(i)
        pfed_edmm31(i) = - pfed_edmm30(i)
        pfed_edmm33(i) = - pfed_edmm32(i)
        pfed_edmm35(i) = - pfed_edmm34(i)
        pfed_edmm37(i) = - pfed_edmm36(i)
        pfed_edmm39(i) = - pfed_edmm38(i)
        pfed_edmm41(i) = - pfed_edmm40(i)
        pfed_edmm43(i) = - pfed_edmm42(i)
        pfed_edmm45(i) = - pfed_edmm44(i)
        pfed_edmm47(i) = - pfed_edmm46(i)
        pfed_edmm49(i) = - pfed_edmm48(i)
        pfed_edmm51(i) = - pfed_edmm50(i)
        pfed_edmm53(i) = - pfed_edmm52(i)
        pfed_edmm55(i) = - pfed_edmm54(i)
        pfed_edmm57(i) = - pfed_edmm56(i)
        pfed_edmm59(i) = - pfed_edmm58(i)
        pfed_edmm61(i) = - pfed_edmm60(i)
        pfed_edmm63(i) = - pfed_edmm62(i)
  100 continue

      vned_edmm2 =  + amass2(1)
      vned_edmm4 =  + amass2(2)
      vned_edmm6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vned_edmm8 =  + amass2(3)
      vned_edmm10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vned_edmm12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vned_edmm14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vned_edmm16 =  + amass2(4)
      vned_edmm18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vned_edmm20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vned_edmm22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vned_edmm24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vned_edmm26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vned_edmm28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vned_edmm30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vned_edmm32 =  + amass2(5)
      vned_edmm34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vned_edmm36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vned_edmm38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vned_edmm40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vned_edmm42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vned_edmm44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vned_edmm46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vned_edmm48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vned_edmm50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vned_edmm52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vned_edmm54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vned_edmm56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vned_edmm58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vned_edmm60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vned_edmm62 =  + amass2(6)
      vned_edmm3 = vned_edmm2
      vned_edmm5 = vned_edmm4
      vned_edmm7 = vned_edmm6
      vned_edmm9 = vned_edmm8
      vned_edmm11 = vned_edmm10
      vned_edmm13 = vned_edmm12
      vned_edmm15 = vned_edmm14
      vned_edmm17 = vned_edmm16
      vned_edmm19 = vned_edmm18
      vned_edmm21 = vned_edmm20
      vned_edmm23 = vned_edmm22
      vned_edmm25 = vned_edmm24
      vned_edmm27 = vned_edmm26
      vned_edmm29 = vned_edmm28
      vned_edmm31 = vned_edmm30
      vned_edmm33 = vned_edmm32
      vned_edmm35 = vned_edmm34
      vned_edmm37 = vned_edmm36
      vned_edmm39 = vned_edmm38
      vned_edmm41 = vned_edmm40
      vned_edmm43 = vned_edmm42
      vned_edmm45 = vned_edmm44
      vned_edmm47 = vned_edmm46
      vned_edmm49 = vned_edmm48
      vned_edmm51 = vned_edmm50
      vned_edmm53 = vned_edmm52
      vned_edmm55 = vned_edmm54
      vned_edmm57 = vned_edmm56
      vned_edmm59 = vned_edmm58
      vned_edmm61 = vned_edmm60
      vned_edmm63 = vned_edmm62
      return
      end
