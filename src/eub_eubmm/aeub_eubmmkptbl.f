*     File eub_eubmm/aeub_eubmmkptbl.f : Sat Mar 18 19:45:00 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aeub_eubmmkptbl
      implicit real*8(a-h,o-z)

      include 'incleub_eubmm1.h'
      include 'inclk.inc'
      include 'incleub_eubmmp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfeub_eubmm2(i) =  + peeub_eubmm1(i)
        pfeub_eubmm4(i) =  + peeub_eubmm2(i)
        pfeub_eubmm6(i) =  + peeub_eubmm1(i) + peeub_eubmm2(i)
        pfeub_eubmm8(i) =  - peeub_eubmm3(i)
        pfeub_eubmm10(i) =  + peeub_eubmm1(i) - peeub_eubmm3(i)
        pfeub_eubmm12(i) =  + peeub_eubmm2(i) - peeub_eubmm3(i)
        pfeub_eubmm14(i) =  + peeub_eubmm1(i) + peeub_eubmm2(i) - peeub_eubmm3(i)
        pfeub_eubmm16(i) =  - peeub_eubmm4(i)
        pfeub_eubmm18(i) =  + peeub_eubmm1(i) - peeub_eubmm4(i)
        pfeub_eubmm20(i) =  + peeub_eubmm2(i) - peeub_eubmm4(i)
        pfeub_eubmm22(i) =  + peeub_eubmm1(i) + peeub_eubmm2(i) - peeub_eubmm4(i)
        pfeub_eubmm24(i) =  - peeub_eubmm3(i) - peeub_eubmm4(i)
        pfeub_eubmm26(i) =  + peeub_eubmm1(i) - peeub_eubmm3(i) - peeub_eubmm4(i)
        pfeub_eubmm28(i) =  + peeub_eubmm2(i) - peeub_eubmm3(i) - peeub_eubmm4(i)
        pfeub_eubmm30(i) =  + peeub_eubmm5(i) + peeub_eubmm6(i)
        pfeub_eubmm32(i) =  - peeub_eubmm5(i)
        pfeub_eubmm34(i) =  + peeub_eubmm1(i) - peeub_eubmm5(i)
        pfeub_eubmm36(i) =  + peeub_eubmm2(i) - peeub_eubmm5(i)
        pfeub_eubmm38(i) =  + peeub_eubmm1(i) + peeub_eubmm2(i) - peeub_eubmm5(i)
        pfeub_eubmm40(i) =  - peeub_eubmm3(i) - peeub_eubmm5(i)
        pfeub_eubmm42(i) =  + peeub_eubmm1(i) - peeub_eubmm3(i) - peeub_eubmm5(i)
        pfeub_eubmm44(i) =  + peeub_eubmm2(i) - peeub_eubmm3(i) - peeub_eubmm5(i)
        pfeub_eubmm46(i) =  + peeub_eubmm4(i) + peeub_eubmm6(i)
        pfeub_eubmm48(i) =  - peeub_eubmm4(i) - peeub_eubmm5(i)
        pfeub_eubmm50(i) =  + peeub_eubmm1(i) - peeub_eubmm4(i) - peeub_eubmm5(i)
        pfeub_eubmm52(i) =  + peeub_eubmm2(i) - peeub_eubmm4(i) - peeub_eubmm5(i)
        pfeub_eubmm54(i) =  + peeub_eubmm3(i) + peeub_eubmm6(i)
        pfeub_eubmm56(i) =  - peeub_eubmm3(i) - peeub_eubmm4(i) - peeub_eubmm5(i)
        pfeub_eubmm58(i) =  - peeub_eubmm2(i) + peeub_eubmm6(i)
        pfeub_eubmm60(i) =  - peeub_eubmm1(i) + peeub_eubmm6(i)
        pfeub_eubmm62(i) =  + peeub_eubmm6(i)

        pfeub_eubmm3(i) = - pfeub_eubmm2(i)
        pfeub_eubmm5(i) = - pfeub_eubmm4(i)
        pfeub_eubmm7(i) = - pfeub_eubmm6(i)
        pfeub_eubmm9(i) = - pfeub_eubmm8(i)
        pfeub_eubmm11(i) = - pfeub_eubmm10(i)
        pfeub_eubmm13(i) = - pfeub_eubmm12(i)
        pfeub_eubmm15(i) = - pfeub_eubmm14(i)
        pfeub_eubmm17(i) = - pfeub_eubmm16(i)
        pfeub_eubmm19(i) = - pfeub_eubmm18(i)
        pfeub_eubmm21(i) = - pfeub_eubmm20(i)
        pfeub_eubmm23(i) = - pfeub_eubmm22(i)
        pfeub_eubmm25(i) = - pfeub_eubmm24(i)
        pfeub_eubmm27(i) = - pfeub_eubmm26(i)
        pfeub_eubmm29(i) = - pfeub_eubmm28(i)
        pfeub_eubmm31(i) = - pfeub_eubmm30(i)
        pfeub_eubmm33(i) = - pfeub_eubmm32(i)
        pfeub_eubmm35(i) = - pfeub_eubmm34(i)
        pfeub_eubmm37(i) = - pfeub_eubmm36(i)
        pfeub_eubmm39(i) = - pfeub_eubmm38(i)
        pfeub_eubmm41(i) = - pfeub_eubmm40(i)
        pfeub_eubmm43(i) = - pfeub_eubmm42(i)
        pfeub_eubmm45(i) = - pfeub_eubmm44(i)
        pfeub_eubmm47(i) = - pfeub_eubmm46(i)
        pfeub_eubmm49(i) = - pfeub_eubmm48(i)
        pfeub_eubmm51(i) = - pfeub_eubmm50(i)
        pfeub_eubmm53(i) = - pfeub_eubmm52(i)
        pfeub_eubmm55(i) = - pfeub_eubmm54(i)
        pfeub_eubmm57(i) = - pfeub_eubmm56(i)
        pfeub_eubmm59(i) = - pfeub_eubmm58(i)
        pfeub_eubmm61(i) = - pfeub_eubmm60(i)
        pfeub_eubmm63(i) = - pfeub_eubmm62(i)
  100 continue

      vneub_eubmm2 =  + amass2(1)
      vneub_eubmm4 =  + amass2(2)
      vneub_eubmm6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vneub_eubmm8 =  + amass2(3)
      vneub_eubmm10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vneub_eubmm12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vneub_eubmm14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vneub_eubmm16 =  + amass2(4)
      vneub_eubmm18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vneub_eubmm20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vneub_eubmm22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vneub_eubmm24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vneub_eubmm26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vneub_eubmm28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vneub_eubmm30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vneub_eubmm32 =  + amass2(5)
      vneub_eubmm34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vneub_eubmm36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vneub_eubmm38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vneub_eubmm40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vneub_eubmm42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vneub_eubmm44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vneub_eubmm46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vneub_eubmm48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vneub_eubmm50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vneub_eubmm52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vneub_eubmm54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vneub_eubmm56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vneub_eubmm58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vneub_eubmm60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vneub_eubmm62 =  + amass2(6)
      vneub_eubmm3 = vneub_eubmm2
      vneub_eubmm5 = vneub_eubmm4
      vneub_eubmm7 = vneub_eubmm6
      vneub_eubmm9 = vneub_eubmm8
      vneub_eubmm11 = vneub_eubmm10
      vneub_eubmm13 = vneub_eubmm12
      vneub_eubmm15 = vneub_eubmm14
      vneub_eubmm17 = vneub_eubmm16
      vneub_eubmm19 = vneub_eubmm18
      vneub_eubmm21 = vneub_eubmm20
      vneub_eubmm23 = vneub_eubmm22
      vneub_eubmm25 = vneub_eubmm24
      vneub_eubmm27 = vneub_eubmm26
      vneub_eubmm29 = vneub_eubmm28
      vneub_eubmm31 = vneub_eubmm30
      vneub_eubmm33 = vneub_eubmm32
      vneub_eubmm35 = vneub_eubmm34
      vneub_eubmm37 = vneub_eubmm36
      vneub_eubmm39 = vneub_eubmm38
      vneub_eubmm41 = vneub_eubmm40
      vneub_eubmm43 = vneub_eubmm42
      vneub_eubmm45 = vneub_eubmm44
      vneub_eubmm47 = vneub_eubmm46
      vneub_eubmm49 = vneub_eubmm48
      vneub_eubmm51 = vneub_eubmm50
      vneub_eubmm53 = vneub_eubmm52
      vneub_eubmm55 = vneub_eubmm54
      vneub_eubmm57 = vneub_eubmm56
      vneub_eubmm59 = vneub_eubmm58
      vneub_eubmm61 = vneub_eubmm60
      vneub_eubmm63 = vneub_eubmm62
      return
      end
