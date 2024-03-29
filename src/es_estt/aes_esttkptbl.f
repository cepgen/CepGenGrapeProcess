*     File es_estt/aes_esttkptbl.f : Sat Mar 18 19:45:02 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aes_esttkptbl
      implicit real*8(a-h,o-z)

      include 'incles_estt1.h'
      include 'inclk.inc'
      include 'incles_esttp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfes_estt2(i) =  + pees_estt1(i)
        pfes_estt4(i) =  + pees_estt2(i)
        pfes_estt6(i) =  + pees_estt1(i) + pees_estt2(i)
        pfes_estt8(i) =  - pees_estt3(i)
        pfes_estt10(i) =  + pees_estt1(i) - pees_estt3(i)
        pfes_estt12(i) =  + pees_estt2(i) - pees_estt3(i)
        pfes_estt14(i) =  + pees_estt1(i) + pees_estt2(i) - pees_estt3(i)
        pfes_estt16(i) =  - pees_estt4(i)
        pfes_estt18(i) =  + pees_estt1(i) - pees_estt4(i)
        pfes_estt20(i) =  + pees_estt2(i) - pees_estt4(i)
        pfes_estt22(i) =  + pees_estt1(i) + pees_estt2(i) - pees_estt4(i)
        pfes_estt24(i) =  - pees_estt3(i) - pees_estt4(i)
        pfes_estt26(i) =  + pees_estt1(i) - pees_estt3(i) - pees_estt4(i)
        pfes_estt28(i) =  + pees_estt2(i) - pees_estt3(i) - pees_estt4(i)
        pfes_estt30(i) =  + pees_estt5(i) + pees_estt6(i)
        pfes_estt32(i) =  - pees_estt5(i)
        pfes_estt34(i) =  + pees_estt1(i) - pees_estt5(i)
        pfes_estt36(i) =  + pees_estt2(i) - pees_estt5(i)
        pfes_estt38(i) =  + pees_estt1(i) + pees_estt2(i) - pees_estt5(i)
        pfes_estt40(i) =  - pees_estt3(i) - pees_estt5(i)
        pfes_estt42(i) =  + pees_estt1(i) - pees_estt3(i) - pees_estt5(i)
        pfes_estt44(i) =  + pees_estt2(i) - pees_estt3(i) - pees_estt5(i)
        pfes_estt46(i) =  + pees_estt4(i) + pees_estt6(i)
        pfes_estt48(i) =  - pees_estt4(i) - pees_estt5(i)
        pfes_estt50(i) =  + pees_estt1(i) - pees_estt4(i) - pees_estt5(i)
        pfes_estt52(i) =  + pees_estt2(i) - pees_estt4(i) - pees_estt5(i)
        pfes_estt54(i) =  + pees_estt3(i) + pees_estt6(i)
        pfes_estt56(i) =  - pees_estt3(i) - pees_estt4(i) - pees_estt5(i)
        pfes_estt58(i) =  - pees_estt2(i) + pees_estt6(i)
        pfes_estt60(i) =  - pees_estt1(i) + pees_estt6(i)
        pfes_estt62(i) =  + pees_estt6(i)

        pfes_estt3(i) = - pfes_estt2(i)
        pfes_estt5(i) = - pfes_estt4(i)
        pfes_estt7(i) = - pfes_estt6(i)
        pfes_estt9(i) = - pfes_estt8(i)
        pfes_estt11(i) = - pfes_estt10(i)
        pfes_estt13(i) = - pfes_estt12(i)
        pfes_estt15(i) = - pfes_estt14(i)
        pfes_estt17(i) = - pfes_estt16(i)
        pfes_estt19(i) = - pfes_estt18(i)
        pfes_estt21(i) = - pfes_estt20(i)
        pfes_estt23(i) = - pfes_estt22(i)
        pfes_estt25(i) = - pfes_estt24(i)
        pfes_estt27(i) = - pfes_estt26(i)
        pfes_estt29(i) = - pfes_estt28(i)
        pfes_estt31(i) = - pfes_estt30(i)
        pfes_estt33(i) = - pfes_estt32(i)
        pfes_estt35(i) = - pfes_estt34(i)
        pfes_estt37(i) = - pfes_estt36(i)
        pfes_estt39(i) = - pfes_estt38(i)
        pfes_estt41(i) = - pfes_estt40(i)
        pfes_estt43(i) = - pfes_estt42(i)
        pfes_estt45(i) = - pfes_estt44(i)
        pfes_estt47(i) = - pfes_estt46(i)
        pfes_estt49(i) = - pfes_estt48(i)
        pfes_estt51(i) = - pfes_estt50(i)
        pfes_estt53(i) = - pfes_estt52(i)
        pfes_estt55(i) = - pfes_estt54(i)
        pfes_estt57(i) = - pfes_estt56(i)
        pfes_estt59(i) = - pfes_estt58(i)
        pfes_estt61(i) = - pfes_estt60(i)
        pfes_estt63(i) = - pfes_estt62(i)
  100 continue

      vnes_estt2 =  + amass2(1)
      vnes_estt4 =  + amass2(2)
      vnes_estt6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vnes_estt8 =  + amass2(3)
      vnes_estt10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vnes_estt12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vnes_estt14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vnes_estt16 =  + amass2(4)
      vnes_estt18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vnes_estt20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vnes_estt22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vnes_estt24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vnes_estt26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vnes_estt28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vnes_estt30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vnes_estt32 =  + amass2(5)
      vnes_estt34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vnes_estt36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vnes_estt38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vnes_estt40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vnes_estt42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vnes_estt44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vnes_estt46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vnes_estt48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vnes_estt50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vnes_estt52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vnes_estt54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vnes_estt56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vnes_estt58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vnes_estt60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vnes_estt62 =  + amass2(6)
      vnes_estt3 = vnes_estt2
      vnes_estt5 = vnes_estt4
      vnes_estt7 = vnes_estt6
      vnes_estt9 = vnes_estt8
      vnes_estt11 = vnes_estt10
      vnes_estt13 = vnes_estt12
      vnes_estt15 = vnes_estt14
      vnes_estt17 = vnes_estt16
      vnes_estt19 = vnes_estt18
      vnes_estt21 = vnes_estt20
      vnes_estt23 = vnes_estt22
      vnes_estt25 = vnes_estt24
      vnes_estt27 = vnes_estt26
      vnes_estt29 = vnes_estt28
      vnes_estt31 = vnes_estt30
      vnes_estt33 = vnes_estt32
      vnes_estt35 = vnes_estt34
      vnes_estt37 = vnes_estt36
      vnes_estt39 = vnes_estt38
      vnes_estt41 = vnes_estt40
      vnes_estt43 = vnes_estt42
      vnes_estt45 = vnes_estt44
      vnes_estt47 = vnes_estt46
      vnes_estt49 = vnes_estt48
      vnes_estt51 = vnes_estt50
      vnes_estt53 = vnes_estt52
      vnes_estt55 = vnes_estt54
      vnes_estt57 = vnes_estt56
      vnes_estt59 = vnes_estt58
      vnes_estt61 = vnes_estt60
      vnes_estt63 = vnes_estt62
      return
      end
