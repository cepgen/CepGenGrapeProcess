*     File eu_euee/aeu_eueekptbl.f : Sat Mar 18 19:44:59 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aeu_eueekptbl
      implicit real*8(a-h,o-z)

      include 'incleu_euee1.h'
      include 'inclk.inc'
      include 'incleu_eueep.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfeu_euee2(i) =  + peeu_euee1(i)
        pfeu_euee4(i) =  + peeu_euee2(i)
        pfeu_euee6(i) =  + peeu_euee1(i) + peeu_euee2(i)
        pfeu_euee8(i) =  - peeu_euee3(i)
        pfeu_euee10(i) =  + peeu_euee1(i) - peeu_euee3(i)
        pfeu_euee12(i) =  + peeu_euee2(i) - peeu_euee3(i)
        pfeu_euee14(i) =  + peeu_euee1(i) + peeu_euee2(i) - peeu_euee3(i)
        pfeu_euee16(i) =  - peeu_euee4(i)
        pfeu_euee18(i) =  + peeu_euee1(i) - peeu_euee4(i)
        pfeu_euee20(i) =  + peeu_euee2(i) - peeu_euee4(i)
        pfeu_euee22(i) =  + peeu_euee1(i) + peeu_euee2(i) - peeu_euee4(i)
        pfeu_euee24(i) =  - peeu_euee3(i) - peeu_euee4(i)
        pfeu_euee26(i) =  + peeu_euee1(i) - peeu_euee3(i) - peeu_euee4(i)
        pfeu_euee28(i) =  + peeu_euee2(i) - peeu_euee3(i) - peeu_euee4(i)
        pfeu_euee30(i) =  + peeu_euee5(i) + peeu_euee6(i)
        pfeu_euee32(i) =  - peeu_euee5(i)
        pfeu_euee34(i) =  + peeu_euee1(i) - peeu_euee5(i)
        pfeu_euee36(i) =  + peeu_euee2(i) - peeu_euee5(i)
        pfeu_euee38(i) =  + peeu_euee1(i) + peeu_euee2(i) - peeu_euee5(i)
        pfeu_euee40(i) =  - peeu_euee3(i) - peeu_euee5(i)
        pfeu_euee42(i) =  + peeu_euee1(i) - peeu_euee3(i) - peeu_euee5(i)
        pfeu_euee44(i) =  + peeu_euee2(i) - peeu_euee3(i) - peeu_euee5(i)
        pfeu_euee46(i) =  + peeu_euee4(i) + peeu_euee6(i)
        pfeu_euee48(i) =  - peeu_euee4(i) - peeu_euee5(i)
        pfeu_euee50(i) =  + peeu_euee1(i) - peeu_euee4(i) - peeu_euee5(i)
        pfeu_euee52(i) =  + peeu_euee2(i) - peeu_euee4(i) - peeu_euee5(i)
        pfeu_euee54(i) =  + peeu_euee3(i) + peeu_euee6(i)
        pfeu_euee56(i) =  - peeu_euee3(i) - peeu_euee4(i) - peeu_euee5(i)
        pfeu_euee58(i) =  - peeu_euee2(i) + peeu_euee6(i)
        pfeu_euee60(i) =  - peeu_euee1(i) + peeu_euee6(i)
        pfeu_euee62(i) =  + peeu_euee6(i)

        pfeu_euee3(i) = - pfeu_euee2(i)
        pfeu_euee5(i) = - pfeu_euee4(i)
        pfeu_euee7(i) = - pfeu_euee6(i)
        pfeu_euee9(i) = - pfeu_euee8(i)
        pfeu_euee11(i) = - pfeu_euee10(i)
        pfeu_euee13(i) = - pfeu_euee12(i)
        pfeu_euee15(i) = - pfeu_euee14(i)
        pfeu_euee17(i) = - pfeu_euee16(i)
        pfeu_euee19(i) = - pfeu_euee18(i)
        pfeu_euee21(i) = - pfeu_euee20(i)
        pfeu_euee23(i) = - pfeu_euee22(i)
        pfeu_euee25(i) = - pfeu_euee24(i)
        pfeu_euee27(i) = - pfeu_euee26(i)
        pfeu_euee29(i) = - pfeu_euee28(i)
        pfeu_euee31(i) = - pfeu_euee30(i)
        pfeu_euee33(i) = - pfeu_euee32(i)
        pfeu_euee35(i) = - pfeu_euee34(i)
        pfeu_euee37(i) = - pfeu_euee36(i)
        pfeu_euee39(i) = - pfeu_euee38(i)
        pfeu_euee41(i) = - pfeu_euee40(i)
        pfeu_euee43(i) = - pfeu_euee42(i)
        pfeu_euee45(i) = - pfeu_euee44(i)
        pfeu_euee47(i) = - pfeu_euee46(i)
        pfeu_euee49(i) = - pfeu_euee48(i)
        pfeu_euee51(i) = - pfeu_euee50(i)
        pfeu_euee53(i) = - pfeu_euee52(i)
        pfeu_euee55(i) = - pfeu_euee54(i)
        pfeu_euee57(i) = - pfeu_euee56(i)
        pfeu_euee59(i) = - pfeu_euee58(i)
        pfeu_euee61(i) = - pfeu_euee60(i)
        pfeu_euee63(i) = - pfeu_euee62(i)
  100 continue

      vneu_euee2 =  + amass2(1)
      vneu_euee4 =  + amass2(2)
      vneu_euee6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vneu_euee8 =  + amass2(3)
      vneu_euee10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vneu_euee12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vneu_euee14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vneu_euee16 =  + amass2(4)
      vneu_euee18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vneu_euee20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vneu_euee22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vneu_euee24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vneu_euee26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vneu_euee28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vneu_euee30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vneu_euee32 =  + amass2(5)
      vneu_euee34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vneu_euee36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vneu_euee38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vneu_euee40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vneu_euee42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vneu_euee44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vneu_euee46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vneu_euee48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vneu_euee50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vneu_euee52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vneu_euee54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vneu_euee56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vneu_euee58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vneu_euee60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vneu_euee62 =  + amass2(6)
      vneu_euee3 = vneu_euee2
      vneu_euee5 = vneu_euee4
      vneu_euee7 = vneu_euee6
      vneu_euee9 = vneu_euee8
      vneu_euee11 = vneu_euee10
      vneu_euee13 = vneu_euee12
      vneu_euee15 = vneu_euee14
      vneu_euee17 = vneu_euee16
      vneu_euee19 = vneu_euee18
      vneu_euee21 = vneu_euee20
      vneu_euee23 = vneu_euee22
      vneu_euee25 = vneu_euee24
      vneu_euee27 = vneu_euee26
      vneu_euee29 = vneu_euee28
      vneu_euee31 = vneu_euee30
      vneu_euee33 = vneu_euee32
      vneu_euee35 = vneu_euee34
      vneu_euee37 = vneu_euee36
      vneu_euee39 = vneu_euee38
      vneu_euee41 = vneu_euee40
      vneu_euee43 = vneu_euee42
      vneu_euee45 = vneu_euee44
      vneu_euee47 = vneu_euee46
      vneu_euee49 = vneu_euee48
      vneu_euee51 = vneu_euee50
      vneu_euee53 = vneu_euee52
      vneu_euee55 = vneu_euee54
      vneu_euee57 = vneu_euee56
      vneu_euee59 = vneu_euee58
      vneu_euee61 = vneu_euee60
      vneu_euee63 = vneu_euee62
      return
      end