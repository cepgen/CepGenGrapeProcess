*     File edb_edbmm/aedb_edbmmkptbl.f : Sat Mar 18 19:45:01 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aedb_edbmmkptbl
      implicit real*8(a-h,o-z)

      include 'incledb_edbmm1.h'
      include 'inclk.inc'
      include 'incledb_edbmmp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfedb_edbmm2(i) =  + peedb_edbmm1(i)
        pfedb_edbmm4(i) =  + peedb_edbmm2(i)
        pfedb_edbmm6(i) =  + peedb_edbmm1(i) + peedb_edbmm2(i)
        pfedb_edbmm8(i) =  - peedb_edbmm3(i)
        pfedb_edbmm10(i) =  + peedb_edbmm1(i) - peedb_edbmm3(i)
        pfedb_edbmm12(i) =  + peedb_edbmm2(i) - peedb_edbmm3(i)
        pfedb_edbmm14(i) =  + peedb_edbmm1(i) + peedb_edbmm2(i) - peedb_edbmm3(i)
        pfedb_edbmm16(i) =  - peedb_edbmm4(i)
        pfedb_edbmm18(i) =  + peedb_edbmm1(i) - peedb_edbmm4(i)
        pfedb_edbmm20(i) =  + peedb_edbmm2(i) - peedb_edbmm4(i)
        pfedb_edbmm22(i) =  + peedb_edbmm1(i) + peedb_edbmm2(i) - peedb_edbmm4(i)
        pfedb_edbmm24(i) =  - peedb_edbmm3(i) - peedb_edbmm4(i)
        pfedb_edbmm26(i) =  + peedb_edbmm1(i) - peedb_edbmm3(i) - peedb_edbmm4(i)
        pfedb_edbmm28(i) =  + peedb_edbmm2(i) - peedb_edbmm3(i) - peedb_edbmm4(i)
        pfedb_edbmm30(i) =  + peedb_edbmm5(i) + peedb_edbmm6(i)
        pfedb_edbmm32(i) =  - peedb_edbmm5(i)
        pfedb_edbmm34(i) =  + peedb_edbmm1(i) - peedb_edbmm5(i)
        pfedb_edbmm36(i) =  + peedb_edbmm2(i) - peedb_edbmm5(i)
        pfedb_edbmm38(i) =  + peedb_edbmm1(i) + peedb_edbmm2(i) - peedb_edbmm5(i)
        pfedb_edbmm40(i) =  - peedb_edbmm3(i) - peedb_edbmm5(i)
        pfedb_edbmm42(i) =  + peedb_edbmm1(i) - peedb_edbmm3(i) - peedb_edbmm5(i)
        pfedb_edbmm44(i) =  + peedb_edbmm2(i) - peedb_edbmm3(i) - peedb_edbmm5(i)
        pfedb_edbmm46(i) =  + peedb_edbmm4(i) + peedb_edbmm6(i)
        pfedb_edbmm48(i) =  - peedb_edbmm4(i) - peedb_edbmm5(i)
        pfedb_edbmm50(i) =  + peedb_edbmm1(i) - peedb_edbmm4(i) - peedb_edbmm5(i)
        pfedb_edbmm52(i) =  + peedb_edbmm2(i) - peedb_edbmm4(i) - peedb_edbmm5(i)
        pfedb_edbmm54(i) =  + peedb_edbmm3(i) + peedb_edbmm6(i)
        pfedb_edbmm56(i) =  - peedb_edbmm3(i) - peedb_edbmm4(i) - peedb_edbmm5(i)
        pfedb_edbmm58(i) =  - peedb_edbmm2(i) + peedb_edbmm6(i)
        pfedb_edbmm60(i) =  - peedb_edbmm1(i) + peedb_edbmm6(i)
        pfedb_edbmm62(i) =  + peedb_edbmm6(i)

        pfedb_edbmm3(i) = - pfedb_edbmm2(i)
        pfedb_edbmm5(i) = - pfedb_edbmm4(i)
        pfedb_edbmm7(i) = - pfedb_edbmm6(i)
        pfedb_edbmm9(i) = - pfedb_edbmm8(i)
        pfedb_edbmm11(i) = - pfedb_edbmm10(i)
        pfedb_edbmm13(i) = - pfedb_edbmm12(i)
        pfedb_edbmm15(i) = - pfedb_edbmm14(i)
        pfedb_edbmm17(i) = - pfedb_edbmm16(i)
        pfedb_edbmm19(i) = - pfedb_edbmm18(i)
        pfedb_edbmm21(i) = - pfedb_edbmm20(i)
        pfedb_edbmm23(i) = - pfedb_edbmm22(i)
        pfedb_edbmm25(i) = - pfedb_edbmm24(i)
        pfedb_edbmm27(i) = - pfedb_edbmm26(i)
        pfedb_edbmm29(i) = - pfedb_edbmm28(i)
        pfedb_edbmm31(i) = - pfedb_edbmm30(i)
        pfedb_edbmm33(i) = - pfedb_edbmm32(i)
        pfedb_edbmm35(i) = - pfedb_edbmm34(i)
        pfedb_edbmm37(i) = - pfedb_edbmm36(i)
        pfedb_edbmm39(i) = - pfedb_edbmm38(i)
        pfedb_edbmm41(i) = - pfedb_edbmm40(i)
        pfedb_edbmm43(i) = - pfedb_edbmm42(i)
        pfedb_edbmm45(i) = - pfedb_edbmm44(i)
        pfedb_edbmm47(i) = - pfedb_edbmm46(i)
        pfedb_edbmm49(i) = - pfedb_edbmm48(i)
        pfedb_edbmm51(i) = - pfedb_edbmm50(i)
        pfedb_edbmm53(i) = - pfedb_edbmm52(i)
        pfedb_edbmm55(i) = - pfedb_edbmm54(i)
        pfedb_edbmm57(i) = - pfedb_edbmm56(i)
        pfedb_edbmm59(i) = - pfedb_edbmm58(i)
        pfedb_edbmm61(i) = - pfedb_edbmm60(i)
        pfedb_edbmm63(i) = - pfedb_edbmm62(i)
  100 continue

      vnedb_edbmm2 =  + amass2(1)
      vnedb_edbmm4 =  + amass2(2)
      vnedb_edbmm6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vnedb_edbmm8 =  + amass2(3)
      vnedb_edbmm10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vnedb_edbmm12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vnedb_edbmm14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vnedb_edbmm16 =  + amass2(4)
      vnedb_edbmm18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vnedb_edbmm20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vnedb_edbmm22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vnedb_edbmm24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vnedb_edbmm26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vnedb_edbmm28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vnedb_edbmm30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vnedb_edbmm32 =  + amass2(5)
      vnedb_edbmm34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vnedb_edbmm36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vnedb_edbmm38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vnedb_edbmm40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vnedb_edbmm42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vnedb_edbmm44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vnedb_edbmm46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vnedb_edbmm48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vnedb_edbmm50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vnedb_edbmm52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vnedb_edbmm54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vnedb_edbmm56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vnedb_edbmm58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vnedb_edbmm60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vnedb_edbmm62 =  + amass2(6)
      vnedb_edbmm3 = vnedb_edbmm2
      vnedb_edbmm5 = vnedb_edbmm4
      vnedb_edbmm7 = vnedb_edbmm6
      vnedb_edbmm9 = vnedb_edbmm8
      vnedb_edbmm11 = vnedb_edbmm10
      vnedb_edbmm13 = vnedb_edbmm12
      vnedb_edbmm15 = vnedb_edbmm14
      vnedb_edbmm17 = vnedb_edbmm16
      vnedb_edbmm19 = vnedb_edbmm18
      vnedb_edbmm21 = vnedb_edbmm20
      vnedb_edbmm23 = vnedb_edbmm22
      vnedb_edbmm25 = vnedb_edbmm24
      vnedb_edbmm27 = vnedb_edbmm26
      vnedb_edbmm29 = vnedb_edbmm28
      vnedb_edbmm31 = vnedb_edbmm30
      vnedb_edbmm33 = vnedb_edbmm32
      vnedb_edbmm35 = vnedb_edbmm34
      vnedb_edbmm37 = vnedb_edbmm36
      vnedb_edbmm39 = vnedb_edbmm38
      vnedb_edbmm41 = vnedb_edbmm40
      vnedb_edbmm43 = vnedb_edbmm42
      vnedb_edbmm45 = vnedb_edbmm44
      vnedb_edbmm47 = vnedb_edbmm46
      vnedb_edbmm49 = vnedb_edbmm48
      vnedb_edbmm51 = vnedb_edbmm50
      vnedb_edbmm53 = vnedb_edbmm52
      vnedb_edbmm55 = vnedb_edbmm54
      vnedb_edbmm57 = vnedb_edbmm56
      vnedb_edbmm59 = vnedb_edbmm58
      vnedb_edbmm61 = vnedb_edbmm60
      vnedb_edbmm63 = vnedb_edbmm62
      return
      end
