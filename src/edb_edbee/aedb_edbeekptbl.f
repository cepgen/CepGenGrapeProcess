*     File edb_edbee/aedb_edbeekptbl.f : Sat Mar 18 19:45:01 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aedb_edbeekptbl
      implicit real*8(a-h,o-z)

      include 'incledb_edbee1.h'
      include 'inclk.inc'
      include 'incledb_edbeep.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfedb_edbee2(i) =  + peedb_edbee1(i)
        pfedb_edbee4(i) =  + peedb_edbee2(i)
        pfedb_edbee6(i) =  + peedb_edbee1(i) + peedb_edbee2(i)
        pfedb_edbee8(i) =  - peedb_edbee3(i)
        pfedb_edbee10(i) =  + peedb_edbee1(i) - peedb_edbee3(i)
        pfedb_edbee12(i) =  + peedb_edbee2(i) - peedb_edbee3(i)
        pfedb_edbee14(i) =  + peedb_edbee1(i) + peedb_edbee2(i) - peedb_edbee3(i)
        pfedb_edbee16(i) =  - peedb_edbee4(i)
        pfedb_edbee18(i) =  + peedb_edbee1(i) - peedb_edbee4(i)
        pfedb_edbee20(i) =  + peedb_edbee2(i) - peedb_edbee4(i)
        pfedb_edbee22(i) =  + peedb_edbee1(i) + peedb_edbee2(i) - peedb_edbee4(i)
        pfedb_edbee24(i) =  - peedb_edbee3(i) - peedb_edbee4(i)
        pfedb_edbee26(i) =  + peedb_edbee1(i) - peedb_edbee3(i) - peedb_edbee4(i)
        pfedb_edbee28(i) =  + peedb_edbee2(i) - peedb_edbee3(i) - peedb_edbee4(i)
        pfedb_edbee30(i) =  + peedb_edbee5(i) + peedb_edbee6(i)
        pfedb_edbee32(i) =  - peedb_edbee5(i)
        pfedb_edbee34(i) =  + peedb_edbee1(i) - peedb_edbee5(i)
        pfedb_edbee36(i) =  + peedb_edbee2(i) - peedb_edbee5(i)
        pfedb_edbee38(i) =  + peedb_edbee1(i) + peedb_edbee2(i) - peedb_edbee5(i)
        pfedb_edbee40(i) =  - peedb_edbee3(i) - peedb_edbee5(i)
        pfedb_edbee42(i) =  + peedb_edbee1(i) - peedb_edbee3(i) - peedb_edbee5(i)
        pfedb_edbee44(i) =  + peedb_edbee2(i) - peedb_edbee3(i) - peedb_edbee5(i)
        pfedb_edbee46(i) =  + peedb_edbee4(i) + peedb_edbee6(i)
        pfedb_edbee48(i) =  - peedb_edbee4(i) - peedb_edbee5(i)
        pfedb_edbee50(i) =  + peedb_edbee1(i) - peedb_edbee4(i) - peedb_edbee5(i)
        pfedb_edbee52(i) =  + peedb_edbee2(i) - peedb_edbee4(i) - peedb_edbee5(i)
        pfedb_edbee54(i) =  + peedb_edbee3(i) + peedb_edbee6(i)
        pfedb_edbee56(i) =  - peedb_edbee3(i) - peedb_edbee4(i) - peedb_edbee5(i)
        pfedb_edbee58(i) =  - peedb_edbee2(i) + peedb_edbee6(i)
        pfedb_edbee60(i) =  - peedb_edbee1(i) + peedb_edbee6(i)
        pfedb_edbee62(i) =  + peedb_edbee6(i)

        pfedb_edbee3(i) = - pfedb_edbee2(i)
        pfedb_edbee5(i) = - pfedb_edbee4(i)
        pfedb_edbee7(i) = - pfedb_edbee6(i)
        pfedb_edbee9(i) = - pfedb_edbee8(i)
        pfedb_edbee11(i) = - pfedb_edbee10(i)
        pfedb_edbee13(i) = - pfedb_edbee12(i)
        pfedb_edbee15(i) = - pfedb_edbee14(i)
        pfedb_edbee17(i) = - pfedb_edbee16(i)
        pfedb_edbee19(i) = - pfedb_edbee18(i)
        pfedb_edbee21(i) = - pfedb_edbee20(i)
        pfedb_edbee23(i) = - pfedb_edbee22(i)
        pfedb_edbee25(i) = - pfedb_edbee24(i)
        pfedb_edbee27(i) = - pfedb_edbee26(i)
        pfedb_edbee29(i) = - pfedb_edbee28(i)
        pfedb_edbee31(i) = - pfedb_edbee30(i)
        pfedb_edbee33(i) = - pfedb_edbee32(i)
        pfedb_edbee35(i) = - pfedb_edbee34(i)
        pfedb_edbee37(i) = - pfedb_edbee36(i)
        pfedb_edbee39(i) = - pfedb_edbee38(i)
        pfedb_edbee41(i) = - pfedb_edbee40(i)
        pfedb_edbee43(i) = - pfedb_edbee42(i)
        pfedb_edbee45(i) = - pfedb_edbee44(i)
        pfedb_edbee47(i) = - pfedb_edbee46(i)
        pfedb_edbee49(i) = - pfedb_edbee48(i)
        pfedb_edbee51(i) = - pfedb_edbee50(i)
        pfedb_edbee53(i) = - pfedb_edbee52(i)
        pfedb_edbee55(i) = - pfedb_edbee54(i)
        pfedb_edbee57(i) = - pfedb_edbee56(i)
        pfedb_edbee59(i) = - pfedb_edbee58(i)
        pfedb_edbee61(i) = - pfedb_edbee60(i)
        pfedb_edbee63(i) = - pfedb_edbee62(i)
  100 continue

      vnedb_edbee2 =  + amass2(1)
      vnedb_edbee4 =  + amass2(2)
      vnedb_edbee6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vnedb_edbee8 =  + amass2(3)
      vnedb_edbee10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vnedb_edbee12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vnedb_edbee14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vnedb_edbee16 =  + amass2(4)
      vnedb_edbee18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vnedb_edbee20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vnedb_edbee22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vnedb_edbee24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vnedb_edbee26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vnedb_edbee28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vnedb_edbee30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vnedb_edbee32 =  + amass2(5)
      vnedb_edbee34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vnedb_edbee36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vnedb_edbee38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vnedb_edbee40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vnedb_edbee42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vnedb_edbee44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vnedb_edbee46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vnedb_edbee48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vnedb_edbee50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vnedb_edbee52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vnedb_edbee54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vnedb_edbee56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vnedb_edbee58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vnedb_edbee60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vnedb_edbee62 =  + amass2(6)
      vnedb_edbee3 = vnedb_edbee2
      vnedb_edbee5 = vnedb_edbee4
      vnedb_edbee7 = vnedb_edbee6
      vnedb_edbee9 = vnedb_edbee8
      vnedb_edbee11 = vnedb_edbee10
      vnedb_edbee13 = vnedb_edbee12
      vnedb_edbee15 = vnedb_edbee14
      vnedb_edbee17 = vnedb_edbee16
      vnedb_edbee19 = vnedb_edbee18
      vnedb_edbee21 = vnedb_edbee20
      vnedb_edbee23 = vnedb_edbee22
      vnedb_edbee25 = vnedb_edbee24
      vnedb_edbee27 = vnedb_edbee26
      vnedb_edbee29 = vnedb_edbee28
      vnedb_edbee31 = vnedb_edbee30
      vnedb_edbee33 = vnedb_edbee32
      vnedb_edbee35 = vnedb_edbee34
      vnedb_edbee37 = vnedb_edbee36
      vnedb_edbee39 = vnedb_edbee38
      vnedb_edbee41 = vnedb_edbee40
      vnedb_edbee43 = vnedb_edbee42
      vnedb_edbee45 = vnedb_edbee44
      vnedb_edbee47 = vnedb_edbee46
      vnedb_edbee49 = vnedb_edbee48
      vnedb_edbee51 = vnedb_edbee50
      vnedb_edbee53 = vnedb_edbee52
      vnedb_edbee55 = vnedb_edbee54
      vnedb_edbee57 = vnedb_edbee56
      vnedb_edbee59 = vnedb_edbee58
      vnedb_edbee61 = vnedb_edbee60
      vnedb_edbee63 = vnedb_edbee62
      return
      end
