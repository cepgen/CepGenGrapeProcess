*     File eb_ebtt/aeb_ebttkptbl.f : Sat Mar 18 19:45:05 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aeb_ebttkptbl
      implicit real*8(a-h,o-z)

      include 'incleb_ebtt1.h'
      include 'inclk.inc'
      include 'incleb_ebttp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfeb_ebtt2(i) =  + peeb_ebtt1(i)
        pfeb_ebtt4(i) =  + peeb_ebtt2(i)
        pfeb_ebtt6(i) =  + peeb_ebtt1(i) + peeb_ebtt2(i)
        pfeb_ebtt8(i) =  - peeb_ebtt3(i)
        pfeb_ebtt10(i) =  + peeb_ebtt1(i) - peeb_ebtt3(i)
        pfeb_ebtt12(i) =  + peeb_ebtt2(i) - peeb_ebtt3(i)
        pfeb_ebtt14(i) =  + peeb_ebtt1(i) + peeb_ebtt2(i) - peeb_ebtt3(i)
        pfeb_ebtt16(i) =  - peeb_ebtt4(i)
        pfeb_ebtt18(i) =  + peeb_ebtt1(i) - peeb_ebtt4(i)
        pfeb_ebtt20(i) =  + peeb_ebtt2(i) - peeb_ebtt4(i)
        pfeb_ebtt22(i) =  + peeb_ebtt1(i) + peeb_ebtt2(i) - peeb_ebtt4(i)
        pfeb_ebtt24(i) =  - peeb_ebtt3(i) - peeb_ebtt4(i)
        pfeb_ebtt26(i) =  + peeb_ebtt1(i) - peeb_ebtt3(i) - peeb_ebtt4(i)
        pfeb_ebtt28(i) =  + peeb_ebtt2(i) - peeb_ebtt3(i) - peeb_ebtt4(i)
        pfeb_ebtt30(i) =  + peeb_ebtt5(i) + peeb_ebtt6(i)
        pfeb_ebtt32(i) =  - peeb_ebtt5(i)
        pfeb_ebtt34(i) =  + peeb_ebtt1(i) - peeb_ebtt5(i)
        pfeb_ebtt36(i) =  + peeb_ebtt2(i) - peeb_ebtt5(i)
        pfeb_ebtt38(i) =  + peeb_ebtt1(i) + peeb_ebtt2(i) - peeb_ebtt5(i)
        pfeb_ebtt40(i) =  - peeb_ebtt3(i) - peeb_ebtt5(i)
        pfeb_ebtt42(i) =  + peeb_ebtt1(i) - peeb_ebtt3(i) - peeb_ebtt5(i)
        pfeb_ebtt44(i) =  + peeb_ebtt2(i) - peeb_ebtt3(i) - peeb_ebtt5(i)
        pfeb_ebtt46(i) =  + peeb_ebtt4(i) + peeb_ebtt6(i)
        pfeb_ebtt48(i) =  - peeb_ebtt4(i) - peeb_ebtt5(i)
        pfeb_ebtt50(i) =  + peeb_ebtt1(i) - peeb_ebtt4(i) - peeb_ebtt5(i)
        pfeb_ebtt52(i) =  + peeb_ebtt2(i) - peeb_ebtt4(i) - peeb_ebtt5(i)
        pfeb_ebtt54(i) =  + peeb_ebtt3(i) + peeb_ebtt6(i)
        pfeb_ebtt56(i) =  - peeb_ebtt3(i) - peeb_ebtt4(i) - peeb_ebtt5(i)
        pfeb_ebtt58(i) =  - peeb_ebtt2(i) + peeb_ebtt6(i)
        pfeb_ebtt60(i) =  - peeb_ebtt1(i) + peeb_ebtt6(i)
        pfeb_ebtt62(i) =  + peeb_ebtt6(i)

        pfeb_ebtt3(i) = - pfeb_ebtt2(i)
        pfeb_ebtt5(i) = - pfeb_ebtt4(i)
        pfeb_ebtt7(i) = - pfeb_ebtt6(i)
        pfeb_ebtt9(i) = - pfeb_ebtt8(i)
        pfeb_ebtt11(i) = - pfeb_ebtt10(i)
        pfeb_ebtt13(i) = - pfeb_ebtt12(i)
        pfeb_ebtt15(i) = - pfeb_ebtt14(i)
        pfeb_ebtt17(i) = - pfeb_ebtt16(i)
        pfeb_ebtt19(i) = - pfeb_ebtt18(i)
        pfeb_ebtt21(i) = - pfeb_ebtt20(i)
        pfeb_ebtt23(i) = - pfeb_ebtt22(i)
        pfeb_ebtt25(i) = - pfeb_ebtt24(i)
        pfeb_ebtt27(i) = - pfeb_ebtt26(i)
        pfeb_ebtt29(i) = - pfeb_ebtt28(i)
        pfeb_ebtt31(i) = - pfeb_ebtt30(i)
        pfeb_ebtt33(i) = - pfeb_ebtt32(i)
        pfeb_ebtt35(i) = - pfeb_ebtt34(i)
        pfeb_ebtt37(i) = - pfeb_ebtt36(i)
        pfeb_ebtt39(i) = - pfeb_ebtt38(i)
        pfeb_ebtt41(i) = - pfeb_ebtt40(i)
        pfeb_ebtt43(i) = - pfeb_ebtt42(i)
        pfeb_ebtt45(i) = - pfeb_ebtt44(i)
        pfeb_ebtt47(i) = - pfeb_ebtt46(i)
        pfeb_ebtt49(i) = - pfeb_ebtt48(i)
        pfeb_ebtt51(i) = - pfeb_ebtt50(i)
        pfeb_ebtt53(i) = - pfeb_ebtt52(i)
        pfeb_ebtt55(i) = - pfeb_ebtt54(i)
        pfeb_ebtt57(i) = - pfeb_ebtt56(i)
        pfeb_ebtt59(i) = - pfeb_ebtt58(i)
        pfeb_ebtt61(i) = - pfeb_ebtt60(i)
        pfeb_ebtt63(i) = - pfeb_ebtt62(i)
  100 continue

      vneb_ebtt2 =  + amass2(1)
      vneb_ebtt4 =  + amass2(2)
      vneb_ebtt6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vneb_ebtt8 =  + amass2(3)
      vneb_ebtt10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vneb_ebtt12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vneb_ebtt14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vneb_ebtt16 =  + amass2(4)
      vneb_ebtt18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vneb_ebtt20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vneb_ebtt22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vneb_ebtt24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vneb_ebtt26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vneb_ebtt28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vneb_ebtt30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vneb_ebtt32 =  + amass2(5)
      vneb_ebtt34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vneb_ebtt36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vneb_ebtt38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vneb_ebtt40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vneb_ebtt42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vneb_ebtt44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vneb_ebtt46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vneb_ebtt48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vneb_ebtt50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vneb_ebtt52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vneb_ebtt54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vneb_ebtt56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vneb_ebtt58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vneb_ebtt60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vneb_ebtt62 =  + amass2(6)
      vneb_ebtt3 = vneb_ebtt2
      vneb_ebtt5 = vneb_ebtt4
      vneb_ebtt7 = vneb_ebtt6
      vneb_ebtt9 = vneb_ebtt8
      vneb_ebtt11 = vneb_ebtt10
      vneb_ebtt13 = vneb_ebtt12
      vneb_ebtt15 = vneb_ebtt14
      vneb_ebtt17 = vneb_ebtt16
      vneb_ebtt19 = vneb_ebtt18
      vneb_ebtt21 = vneb_ebtt20
      vneb_ebtt23 = vneb_ebtt22
      vneb_ebtt25 = vneb_ebtt24
      vneb_ebtt27 = vneb_ebtt26
      vneb_ebtt29 = vneb_ebtt28
      vneb_ebtt31 = vneb_ebtt30
      vneb_ebtt33 = vneb_ebtt32
      vneb_ebtt35 = vneb_ebtt34
      vneb_ebtt37 = vneb_ebtt36
      vneb_ebtt39 = vneb_ebtt38
      vneb_ebtt41 = vneb_ebtt40
      vneb_ebtt43 = vneb_ebtt42
      vneb_ebtt45 = vneb_ebtt44
      vneb_ebtt47 = vneb_ebtt46
      vneb_ebtt49 = vneb_ebtt48
      vneb_ebtt51 = vneb_ebtt50
      vneb_ebtt53 = vneb_ebtt52
      vneb_ebtt55 = vneb_ebtt54
      vneb_ebtt57 = vneb_ebtt56
      vneb_ebtt59 = vneb_ebtt58
      vneb_ebtt61 = vneb_ebtt60
      vneb_ebtt63 = vneb_ebtt62
      return
      end
