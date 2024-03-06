*     File ed_edtt/aed_edttkptbl.f : Sat Mar 18 19:45:00 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aed_edttkptbl
      implicit real*8(a-h,o-z)

      include 'incled_edtt1.h'
      include 'inclk.inc'
      include 'incled_edttp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfed_edtt2(i) =  + peed_edtt1(i)
        pfed_edtt4(i) =  + peed_edtt2(i)
        pfed_edtt6(i) =  + peed_edtt1(i) + peed_edtt2(i)
        pfed_edtt8(i) =  - peed_edtt3(i)
        pfed_edtt10(i) =  + peed_edtt1(i) - peed_edtt3(i)
        pfed_edtt12(i) =  + peed_edtt2(i) - peed_edtt3(i)
        pfed_edtt14(i) =  + peed_edtt1(i) + peed_edtt2(i) - peed_edtt3(i)
        pfed_edtt16(i) =  - peed_edtt4(i)
        pfed_edtt18(i) =  + peed_edtt1(i) - peed_edtt4(i)
        pfed_edtt20(i) =  + peed_edtt2(i) - peed_edtt4(i)
        pfed_edtt22(i) =  + peed_edtt1(i) + peed_edtt2(i) - peed_edtt4(i)
        pfed_edtt24(i) =  - peed_edtt3(i) - peed_edtt4(i)
        pfed_edtt26(i) =  + peed_edtt1(i) - peed_edtt3(i) - peed_edtt4(i)
        pfed_edtt28(i) =  + peed_edtt2(i) - peed_edtt3(i) - peed_edtt4(i)
        pfed_edtt30(i) =  + peed_edtt5(i) + peed_edtt6(i)
        pfed_edtt32(i) =  - peed_edtt5(i)
        pfed_edtt34(i) =  + peed_edtt1(i) - peed_edtt5(i)
        pfed_edtt36(i) =  + peed_edtt2(i) - peed_edtt5(i)
        pfed_edtt38(i) =  + peed_edtt1(i) + peed_edtt2(i) - peed_edtt5(i)
        pfed_edtt40(i) =  - peed_edtt3(i) - peed_edtt5(i)
        pfed_edtt42(i) =  + peed_edtt1(i) - peed_edtt3(i) - peed_edtt5(i)
        pfed_edtt44(i) =  + peed_edtt2(i) - peed_edtt3(i) - peed_edtt5(i)
        pfed_edtt46(i) =  + peed_edtt4(i) + peed_edtt6(i)
        pfed_edtt48(i) =  - peed_edtt4(i) - peed_edtt5(i)
        pfed_edtt50(i) =  + peed_edtt1(i) - peed_edtt4(i) - peed_edtt5(i)
        pfed_edtt52(i) =  + peed_edtt2(i) - peed_edtt4(i) - peed_edtt5(i)
        pfed_edtt54(i) =  + peed_edtt3(i) + peed_edtt6(i)
        pfed_edtt56(i) =  - peed_edtt3(i) - peed_edtt4(i) - peed_edtt5(i)
        pfed_edtt58(i) =  - peed_edtt2(i) + peed_edtt6(i)
        pfed_edtt60(i) =  - peed_edtt1(i) + peed_edtt6(i)
        pfed_edtt62(i) =  + peed_edtt6(i)

        pfed_edtt3(i) = - pfed_edtt2(i)
        pfed_edtt5(i) = - pfed_edtt4(i)
        pfed_edtt7(i) = - pfed_edtt6(i)
        pfed_edtt9(i) = - pfed_edtt8(i)
        pfed_edtt11(i) = - pfed_edtt10(i)
        pfed_edtt13(i) = - pfed_edtt12(i)
        pfed_edtt15(i) = - pfed_edtt14(i)
        pfed_edtt17(i) = - pfed_edtt16(i)
        pfed_edtt19(i) = - pfed_edtt18(i)
        pfed_edtt21(i) = - pfed_edtt20(i)
        pfed_edtt23(i) = - pfed_edtt22(i)
        pfed_edtt25(i) = - pfed_edtt24(i)
        pfed_edtt27(i) = - pfed_edtt26(i)
        pfed_edtt29(i) = - pfed_edtt28(i)
        pfed_edtt31(i) = - pfed_edtt30(i)
        pfed_edtt33(i) = - pfed_edtt32(i)
        pfed_edtt35(i) = - pfed_edtt34(i)
        pfed_edtt37(i) = - pfed_edtt36(i)
        pfed_edtt39(i) = - pfed_edtt38(i)
        pfed_edtt41(i) = - pfed_edtt40(i)
        pfed_edtt43(i) = - pfed_edtt42(i)
        pfed_edtt45(i) = - pfed_edtt44(i)
        pfed_edtt47(i) = - pfed_edtt46(i)
        pfed_edtt49(i) = - pfed_edtt48(i)
        pfed_edtt51(i) = - pfed_edtt50(i)
        pfed_edtt53(i) = - pfed_edtt52(i)
        pfed_edtt55(i) = - pfed_edtt54(i)
        pfed_edtt57(i) = - pfed_edtt56(i)
        pfed_edtt59(i) = - pfed_edtt58(i)
        pfed_edtt61(i) = - pfed_edtt60(i)
        pfed_edtt63(i) = - pfed_edtt62(i)
  100 continue

      vned_edtt2 =  + amass2(1)
      vned_edtt4 =  + amass2(2)
      vned_edtt6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vned_edtt8 =  + amass2(3)
      vned_edtt10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vned_edtt12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vned_edtt14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vned_edtt16 =  + amass2(4)
      vned_edtt18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vned_edtt20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vned_edtt22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vned_edtt24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vned_edtt26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vned_edtt28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vned_edtt30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vned_edtt32 =  + amass2(5)
      vned_edtt34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vned_edtt36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vned_edtt38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vned_edtt40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vned_edtt42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vned_edtt44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vned_edtt46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vned_edtt48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vned_edtt50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vned_edtt52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vned_edtt54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vned_edtt56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vned_edtt58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vned_edtt60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vned_edtt62 =  + amass2(6)
      vned_edtt3 = vned_edtt2
      vned_edtt5 = vned_edtt4
      vned_edtt7 = vned_edtt6
      vned_edtt9 = vned_edtt8
      vned_edtt11 = vned_edtt10
      vned_edtt13 = vned_edtt12
      vned_edtt15 = vned_edtt14
      vned_edtt17 = vned_edtt16
      vned_edtt19 = vned_edtt18
      vned_edtt21 = vned_edtt20
      vned_edtt23 = vned_edtt22
      vned_edtt25 = vned_edtt24
      vned_edtt27 = vned_edtt26
      vned_edtt29 = vned_edtt28
      vned_edtt31 = vned_edtt30
      vned_edtt33 = vned_edtt32
      vned_edtt35 = vned_edtt34
      vned_edtt37 = vned_edtt36
      vned_edtt39 = vned_edtt38
      vned_edtt41 = vned_edtt40
      vned_edtt43 = vned_edtt42
      vned_edtt45 = vned_edtt44
      vned_edtt47 = vned_edtt46
      vned_edtt49 = vned_edtt48
      vned_edtt51 = vned_edtt50
      vned_edtt53 = vned_edtt52
      vned_edtt55 = vned_edtt54
      vned_edtt57 = vned_edtt56
      vned_edtt59 = vned_edtt58
      vned_edtt61 = vned_edtt60
      vned_edtt63 = vned_edtt62
      return
      end
