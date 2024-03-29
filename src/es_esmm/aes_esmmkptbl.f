*     File es_esmm/aes_esmmkptbl.f : Sat Mar 18 19:45:02 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aes_esmmkptbl
      implicit real*8(a-h,o-z)

      include 'incles_esmm1.h'
      include 'inclk.inc'
      include 'incles_esmmp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfes_esmm2(i) =  + pees_esmm1(i)
        pfes_esmm4(i) =  + pees_esmm2(i)
        pfes_esmm6(i) =  + pees_esmm1(i) + pees_esmm2(i)
        pfes_esmm8(i) =  - pees_esmm3(i)
        pfes_esmm10(i) =  + pees_esmm1(i) - pees_esmm3(i)
        pfes_esmm12(i) =  + pees_esmm2(i) - pees_esmm3(i)
        pfes_esmm14(i) =  + pees_esmm1(i) + pees_esmm2(i) - pees_esmm3(i)
        pfes_esmm16(i) =  - pees_esmm4(i)
        pfes_esmm18(i) =  + pees_esmm1(i) - pees_esmm4(i)
        pfes_esmm20(i) =  + pees_esmm2(i) - pees_esmm4(i)
        pfes_esmm22(i) =  + pees_esmm1(i) + pees_esmm2(i) - pees_esmm4(i)
        pfes_esmm24(i) =  - pees_esmm3(i) - pees_esmm4(i)
        pfes_esmm26(i) =  + pees_esmm1(i) - pees_esmm3(i) - pees_esmm4(i)
        pfes_esmm28(i) =  + pees_esmm2(i) - pees_esmm3(i) - pees_esmm4(i)
        pfes_esmm30(i) =  + pees_esmm5(i) + pees_esmm6(i)
        pfes_esmm32(i) =  - pees_esmm5(i)
        pfes_esmm34(i) =  + pees_esmm1(i) - pees_esmm5(i)
        pfes_esmm36(i) =  + pees_esmm2(i) - pees_esmm5(i)
        pfes_esmm38(i) =  + pees_esmm1(i) + pees_esmm2(i) - pees_esmm5(i)
        pfes_esmm40(i) =  - pees_esmm3(i) - pees_esmm5(i)
        pfes_esmm42(i) =  + pees_esmm1(i) - pees_esmm3(i) - pees_esmm5(i)
        pfes_esmm44(i) =  + pees_esmm2(i) - pees_esmm3(i) - pees_esmm5(i)
        pfes_esmm46(i) =  + pees_esmm4(i) + pees_esmm6(i)
        pfes_esmm48(i) =  - pees_esmm4(i) - pees_esmm5(i)
        pfes_esmm50(i) =  + pees_esmm1(i) - pees_esmm4(i) - pees_esmm5(i)
        pfes_esmm52(i) =  + pees_esmm2(i) - pees_esmm4(i) - pees_esmm5(i)
        pfes_esmm54(i) =  + pees_esmm3(i) + pees_esmm6(i)
        pfes_esmm56(i) =  - pees_esmm3(i) - pees_esmm4(i) - pees_esmm5(i)
        pfes_esmm58(i) =  - pees_esmm2(i) + pees_esmm6(i)
        pfes_esmm60(i) =  - pees_esmm1(i) + pees_esmm6(i)
        pfes_esmm62(i) =  + pees_esmm6(i)

        pfes_esmm3(i) = - pfes_esmm2(i)
        pfes_esmm5(i) = - pfes_esmm4(i)
        pfes_esmm7(i) = - pfes_esmm6(i)
        pfes_esmm9(i) = - pfes_esmm8(i)
        pfes_esmm11(i) = - pfes_esmm10(i)
        pfes_esmm13(i) = - pfes_esmm12(i)
        pfes_esmm15(i) = - pfes_esmm14(i)
        pfes_esmm17(i) = - pfes_esmm16(i)
        pfes_esmm19(i) = - pfes_esmm18(i)
        pfes_esmm21(i) = - pfes_esmm20(i)
        pfes_esmm23(i) = - pfes_esmm22(i)
        pfes_esmm25(i) = - pfes_esmm24(i)
        pfes_esmm27(i) = - pfes_esmm26(i)
        pfes_esmm29(i) = - pfes_esmm28(i)
        pfes_esmm31(i) = - pfes_esmm30(i)
        pfes_esmm33(i) = - pfes_esmm32(i)
        pfes_esmm35(i) = - pfes_esmm34(i)
        pfes_esmm37(i) = - pfes_esmm36(i)
        pfes_esmm39(i) = - pfes_esmm38(i)
        pfes_esmm41(i) = - pfes_esmm40(i)
        pfes_esmm43(i) = - pfes_esmm42(i)
        pfes_esmm45(i) = - pfes_esmm44(i)
        pfes_esmm47(i) = - pfes_esmm46(i)
        pfes_esmm49(i) = - pfes_esmm48(i)
        pfes_esmm51(i) = - pfes_esmm50(i)
        pfes_esmm53(i) = - pfes_esmm52(i)
        pfes_esmm55(i) = - pfes_esmm54(i)
        pfes_esmm57(i) = - pfes_esmm56(i)
        pfes_esmm59(i) = - pfes_esmm58(i)
        pfes_esmm61(i) = - pfes_esmm60(i)
        pfes_esmm63(i) = - pfes_esmm62(i)
  100 continue

      vnes_esmm2 =  + amass2(1)
      vnes_esmm4 =  + amass2(2)
      vnes_esmm6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vnes_esmm8 =  + amass2(3)
      vnes_esmm10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vnes_esmm12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vnes_esmm14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vnes_esmm16 =  + amass2(4)
      vnes_esmm18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vnes_esmm20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vnes_esmm22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vnes_esmm24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vnes_esmm26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vnes_esmm28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vnes_esmm30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vnes_esmm32 =  + amass2(5)
      vnes_esmm34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vnes_esmm36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vnes_esmm38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vnes_esmm40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vnes_esmm42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vnes_esmm44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vnes_esmm46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vnes_esmm48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vnes_esmm50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vnes_esmm52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vnes_esmm54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vnes_esmm56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vnes_esmm58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vnes_esmm60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vnes_esmm62 =  + amass2(6)
      vnes_esmm3 = vnes_esmm2
      vnes_esmm5 = vnes_esmm4
      vnes_esmm7 = vnes_esmm6
      vnes_esmm9 = vnes_esmm8
      vnes_esmm11 = vnes_esmm10
      vnes_esmm13 = vnes_esmm12
      vnes_esmm15 = vnes_esmm14
      vnes_esmm17 = vnes_esmm16
      vnes_esmm19 = vnes_esmm18
      vnes_esmm21 = vnes_esmm20
      vnes_esmm23 = vnes_esmm22
      vnes_esmm25 = vnes_esmm24
      vnes_esmm27 = vnes_esmm26
      vnes_esmm29 = vnes_esmm28
      vnes_esmm31 = vnes_esmm30
      vnes_esmm33 = vnes_esmm32
      vnes_esmm35 = vnes_esmm34
      vnes_esmm37 = vnes_esmm36
      vnes_esmm39 = vnes_esmm38
      vnes_esmm41 = vnes_esmm40
      vnes_esmm43 = vnes_esmm42
      vnes_esmm45 = vnes_esmm44
      vnes_esmm47 = vnes_esmm46
      vnes_esmm49 = vnes_esmm48
      vnes_esmm51 = vnes_esmm50
      vnes_esmm53 = vnes_esmm52
      vnes_esmm55 = vnes_esmm54
      vnes_esmm57 = vnes_esmm56
      vnes_esmm59 = vnes_esmm58
      vnes_esmm61 = vnes_esmm60
      vnes_esmm63 = vnes_esmm62
      return
      end
