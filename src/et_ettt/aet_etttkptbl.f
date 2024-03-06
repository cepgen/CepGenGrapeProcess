*     File et_ettt/aet_etttkptbl.f : Sat Mar 18 19:45:07 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aet_etttkptbl
      implicit real*8(a-h,o-z)

      include 'inclet_ettt1.h'
      include 'inclk.inc'
      include 'inclet_etttp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfet_ettt2(i) =  + peet_ettt1(i)
        pfet_ettt4(i) =  + peet_ettt2(i)
        pfet_ettt6(i) =  + peet_ettt1(i) + peet_ettt2(i)
        pfet_ettt8(i) =  - peet_ettt3(i)
        pfet_ettt10(i) =  + peet_ettt1(i) - peet_ettt3(i)
        pfet_ettt12(i) =  + peet_ettt2(i) - peet_ettt3(i)
        pfet_ettt14(i) =  + peet_ettt1(i) + peet_ettt2(i) - peet_ettt3(i)
        pfet_ettt16(i) =  - peet_ettt4(i)
        pfet_ettt18(i) =  + peet_ettt1(i) - peet_ettt4(i)
        pfet_ettt20(i) =  + peet_ettt2(i) - peet_ettt4(i)
        pfet_ettt22(i) =  + peet_ettt1(i) + peet_ettt2(i) - peet_ettt4(i)
        pfet_ettt24(i) =  - peet_ettt3(i) - peet_ettt4(i)
        pfet_ettt26(i) =  + peet_ettt1(i) - peet_ettt3(i) - peet_ettt4(i)
        pfet_ettt28(i) =  + peet_ettt2(i) - peet_ettt3(i) - peet_ettt4(i)
        pfet_ettt30(i) =  + peet_ettt5(i) + peet_ettt6(i)
        pfet_ettt32(i) =  - peet_ettt5(i)
        pfet_ettt34(i) =  + peet_ettt1(i) - peet_ettt5(i)
        pfet_ettt36(i) =  + peet_ettt2(i) - peet_ettt5(i)
        pfet_ettt38(i) =  + peet_ettt1(i) + peet_ettt2(i) - peet_ettt5(i)
        pfet_ettt40(i) =  - peet_ettt3(i) - peet_ettt5(i)
        pfet_ettt42(i) =  + peet_ettt1(i) - peet_ettt3(i) - peet_ettt5(i)
        pfet_ettt44(i) =  + peet_ettt2(i) - peet_ettt3(i) - peet_ettt5(i)
        pfet_ettt46(i) =  + peet_ettt4(i) + peet_ettt6(i)
        pfet_ettt48(i) =  - peet_ettt4(i) - peet_ettt5(i)
        pfet_ettt50(i) =  + peet_ettt1(i) - peet_ettt4(i) - peet_ettt5(i)
        pfet_ettt52(i) =  + peet_ettt2(i) - peet_ettt4(i) - peet_ettt5(i)
        pfet_ettt54(i) =  + peet_ettt3(i) + peet_ettt6(i)
        pfet_ettt56(i) =  - peet_ettt3(i) - peet_ettt4(i) - peet_ettt5(i)
        pfet_ettt58(i) =  - peet_ettt2(i) + peet_ettt6(i)
        pfet_ettt60(i) =  - peet_ettt1(i) + peet_ettt6(i)
        pfet_ettt62(i) =  + peet_ettt6(i)

        pfet_ettt3(i) = - pfet_ettt2(i)
        pfet_ettt5(i) = - pfet_ettt4(i)
        pfet_ettt7(i) = - pfet_ettt6(i)
        pfet_ettt9(i) = - pfet_ettt8(i)
        pfet_ettt11(i) = - pfet_ettt10(i)
        pfet_ettt13(i) = - pfet_ettt12(i)
        pfet_ettt15(i) = - pfet_ettt14(i)
        pfet_ettt17(i) = - pfet_ettt16(i)
        pfet_ettt19(i) = - pfet_ettt18(i)
        pfet_ettt21(i) = - pfet_ettt20(i)
        pfet_ettt23(i) = - pfet_ettt22(i)
        pfet_ettt25(i) = - pfet_ettt24(i)
        pfet_ettt27(i) = - pfet_ettt26(i)
        pfet_ettt29(i) = - pfet_ettt28(i)
        pfet_ettt31(i) = - pfet_ettt30(i)
        pfet_ettt33(i) = - pfet_ettt32(i)
        pfet_ettt35(i) = - pfet_ettt34(i)
        pfet_ettt37(i) = - pfet_ettt36(i)
        pfet_ettt39(i) = - pfet_ettt38(i)
        pfet_ettt41(i) = - pfet_ettt40(i)
        pfet_ettt43(i) = - pfet_ettt42(i)
        pfet_ettt45(i) = - pfet_ettt44(i)
        pfet_ettt47(i) = - pfet_ettt46(i)
        pfet_ettt49(i) = - pfet_ettt48(i)
        pfet_ettt51(i) = - pfet_ettt50(i)
        pfet_ettt53(i) = - pfet_ettt52(i)
        pfet_ettt55(i) = - pfet_ettt54(i)
        pfet_ettt57(i) = - pfet_ettt56(i)
        pfet_ettt59(i) = - pfet_ettt58(i)
        pfet_ettt61(i) = - pfet_ettt60(i)
        pfet_ettt63(i) = - pfet_ettt62(i)
  100 continue

      vnet_ettt2 =  + amass2(1)
      vnet_ettt4 =  + amass2(2)
      vnet_ettt6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vnet_ettt8 =  + amass2(3)
      vnet_ettt10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vnet_ettt12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vnet_ettt14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vnet_ettt16 =  + amass2(4)
      vnet_ettt18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vnet_ettt20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vnet_ettt22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vnet_ettt24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vnet_ettt26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vnet_ettt28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vnet_ettt30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vnet_ettt32 =  + amass2(5)
      vnet_ettt34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vnet_ettt36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vnet_ettt38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vnet_ettt40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vnet_ettt42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vnet_ettt44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vnet_ettt46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vnet_ettt48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vnet_ettt50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vnet_ettt52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vnet_ettt54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vnet_ettt56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vnet_ettt58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vnet_ettt60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vnet_ettt62 =  + amass2(6)
      vnet_ettt3 = vnet_ettt2
      vnet_ettt5 = vnet_ettt4
      vnet_ettt7 = vnet_ettt6
      vnet_ettt9 = vnet_ettt8
      vnet_ettt11 = vnet_ettt10
      vnet_ettt13 = vnet_ettt12
      vnet_ettt15 = vnet_ettt14
      vnet_ettt17 = vnet_ettt16
      vnet_ettt19 = vnet_ettt18
      vnet_ettt21 = vnet_ettt20
      vnet_ettt23 = vnet_ettt22
      vnet_ettt25 = vnet_ettt24
      vnet_ettt27 = vnet_ettt26
      vnet_ettt29 = vnet_ettt28
      vnet_ettt31 = vnet_ettt30
      vnet_ettt33 = vnet_ettt32
      vnet_ettt35 = vnet_ettt34
      vnet_ettt37 = vnet_ettt36
      vnet_ettt39 = vnet_ettt38
      vnet_ettt41 = vnet_ettt40
      vnet_ettt43 = vnet_ettt42
      vnet_ettt45 = vnet_ettt44
      vnet_ettt47 = vnet_ettt46
      vnet_ettt49 = vnet_ettt48
      vnet_ettt51 = vnet_ettt50
      vnet_ettt53 = vnet_ettt52
      vnet_ettt55 = vnet_ettt54
      vnet_ettt57 = vnet_ettt56
      vnet_ettt59 = vnet_ettt58
      vnet_ettt61 = vnet_ettt60
      vnet_ettt63 = vnet_ettt62
      return
      end
