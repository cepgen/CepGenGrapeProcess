*     File etb_etbtt/aetb_etbttkptbl.f : Sat Mar 18 19:45:07 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aetb_etbttkptbl
      implicit real*8(a-h,o-z)

      include 'incletb_etbtt1.h'
      include 'inclk.inc'
      include 'incletb_etbttp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfetb_etbtt2(i) =  + peetb_etbtt1(i)
        pfetb_etbtt4(i) =  + peetb_etbtt2(i)
        pfetb_etbtt6(i) =  + peetb_etbtt1(i) + peetb_etbtt2(i)
        pfetb_etbtt8(i) =  - peetb_etbtt3(i)
        pfetb_etbtt10(i) =  + peetb_etbtt1(i) - peetb_etbtt3(i)
        pfetb_etbtt12(i) =  + peetb_etbtt2(i) - peetb_etbtt3(i)
        pfetb_etbtt14(i) =  + peetb_etbtt1(i) + peetb_etbtt2(i) - peetb_etbtt3(i)
        pfetb_etbtt16(i) =  - peetb_etbtt4(i)
        pfetb_etbtt18(i) =  + peetb_etbtt1(i) - peetb_etbtt4(i)
        pfetb_etbtt20(i) =  + peetb_etbtt2(i) - peetb_etbtt4(i)
        pfetb_etbtt22(i) =  + peetb_etbtt1(i) + peetb_etbtt2(i) - peetb_etbtt4(i)
        pfetb_etbtt24(i) =  - peetb_etbtt3(i) - peetb_etbtt4(i)
        pfetb_etbtt26(i) =  + peetb_etbtt1(i) - peetb_etbtt3(i) - peetb_etbtt4(i)
        pfetb_etbtt28(i) =  + peetb_etbtt2(i) - peetb_etbtt3(i) - peetb_etbtt4(i)
        pfetb_etbtt30(i) =  + peetb_etbtt5(i) + peetb_etbtt6(i)
        pfetb_etbtt32(i) =  - peetb_etbtt5(i)
        pfetb_etbtt34(i) =  + peetb_etbtt1(i) - peetb_etbtt5(i)
        pfetb_etbtt36(i) =  + peetb_etbtt2(i) - peetb_etbtt5(i)
        pfetb_etbtt38(i) =  + peetb_etbtt1(i) + peetb_etbtt2(i) - peetb_etbtt5(i)
        pfetb_etbtt40(i) =  - peetb_etbtt3(i) - peetb_etbtt5(i)
        pfetb_etbtt42(i) =  + peetb_etbtt1(i) - peetb_etbtt3(i) - peetb_etbtt5(i)
        pfetb_etbtt44(i) =  + peetb_etbtt2(i) - peetb_etbtt3(i) - peetb_etbtt5(i)
        pfetb_etbtt46(i) =  + peetb_etbtt4(i) + peetb_etbtt6(i)
        pfetb_etbtt48(i) =  - peetb_etbtt4(i) - peetb_etbtt5(i)
        pfetb_etbtt50(i) =  + peetb_etbtt1(i) - peetb_etbtt4(i) - peetb_etbtt5(i)
        pfetb_etbtt52(i) =  + peetb_etbtt2(i) - peetb_etbtt4(i) - peetb_etbtt5(i)
        pfetb_etbtt54(i) =  + peetb_etbtt3(i) + peetb_etbtt6(i)
        pfetb_etbtt56(i) =  - peetb_etbtt3(i) - peetb_etbtt4(i) - peetb_etbtt5(i)
        pfetb_etbtt58(i) =  - peetb_etbtt2(i) + peetb_etbtt6(i)
        pfetb_etbtt60(i) =  - peetb_etbtt1(i) + peetb_etbtt6(i)
        pfetb_etbtt62(i) =  + peetb_etbtt6(i)

        pfetb_etbtt3(i) = - pfetb_etbtt2(i)
        pfetb_etbtt5(i) = - pfetb_etbtt4(i)
        pfetb_etbtt7(i) = - pfetb_etbtt6(i)
        pfetb_etbtt9(i) = - pfetb_etbtt8(i)
        pfetb_etbtt11(i) = - pfetb_etbtt10(i)
        pfetb_etbtt13(i) = - pfetb_etbtt12(i)
        pfetb_etbtt15(i) = - pfetb_etbtt14(i)
        pfetb_etbtt17(i) = - pfetb_etbtt16(i)
        pfetb_etbtt19(i) = - pfetb_etbtt18(i)
        pfetb_etbtt21(i) = - pfetb_etbtt20(i)
        pfetb_etbtt23(i) = - pfetb_etbtt22(i)
        pfetb_etbtt25(i) = - pfetb_etbtt24(i)
        pfetb_etbtt27(i) = - pfetb_etbtt26(i)
        pfetb_etbtt29(i) = - pfetb_etbtt28(i)
        pfetb_etbtt31(i) = - pfetb_etbtt30(i)
        pfetb_etbtt33(i) = - pfetb_etbtt32(i)
        pfetb_etbtt35(i) = - pfetb_etbtt34(i)
        pfetb_etbtt37(i) = - pfetb_etbtt36(i)
        pfetb_etbtt39(i) = - pfetb_etbtt38(i)
        pfetb_etbtt41(i) = - pfetb_etbtt40(i)
        pfetb_etbtt43(i) = - pfetb_etbtt42(i)
        pfetb_etbtt45(i) = - pfetb_etbtt44(i)
        pfetb_etbtt47(i) = - pfetb_etbtt46(i)
        pfetb_etbtt49(i) = - pfetb_etbtt48(i)
        pfetb_etbtt51(i) = - pfetb_etbtt50(i)
        pfetb_etbtt53(i) = - pfetb_etbtt52(i)
        pfetb_etbtt55(i) = - pfetb_etbtt54(i)
        pfetb_etbtt57(i) = - pfetb_etbtt56(i)
        pfetb_etbtt59(i) = - pfetb_etbtt58(i)
        pfetb_etbtt61(i) = - pfetb_etbtt60(i)
        pfetb_etbtt63(i) = - pfetb_etbtt62(i)
  100 continue

      vnetb_etbtt2 =  + amass2(1)
      vnetb_etbtt4 =  + amass2(2)
      vnetb_etbtt6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vnetb_etbtt8 =  + amass2(3)
      vnetb_etbtt10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vnetb_etbtt12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vnetb_etbtt14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vnetb_etbtt16 =  + amass2(4)
      vnetb_etbtt18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vnetb_etbtt20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vnetb_etbtt22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vnetb_etbtt24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vnetb_etbtt26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vnetb_etbtt28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vnetb_etbtt30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vnetb_etbtt32 =  + amass2(5)
      vnetb_etbtt34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vnetb_etbtt36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vnetb_etbtt38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vnetb_etbtt40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vnetb_etbtt42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vnetb_etbtt44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vnetb_etbtt46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vnetb_etbtt48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vnetb_etbtt50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vnetb_etbtt52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vnetb_etbtt54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vnetb_etbtt56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vnetb_etbtt58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vnetb_etbtt60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vnetb_etbtt62 =  + amass2(6)
      vnetb_etbtt3 = vnetb_etbtt2
      vnetb_etbtt5 = vnetb_etbtt4
      vnetb_etbtt7 = vnetb_etbtt6
      vnetb_etbtt9 = vnetb_etbtt8
      vnetb_etbtt11 = vnetb_etbtt10
      vnetb_etbtt13 = vnetb_etbtt12
      vnetb_etbtt15 = vnetb_etbtt14
      vnetb_etbtt17 = vnetb_etbtt16
      vnetb_etbtt19 = vnetb_etbtt18
      vnetb_etbtt21 = vnetb_etbtt20
      vnetb_etbtt23 = vnetb_etbtt22
      vnetb_etbtt25 = vnetb_etbtt24
      vnetb_etbtt27 = vnetb_etbtt26
      vnetb_etbtt29 = vnetb_etbtt28
      vnetb_etbtt31 = vnetb_etbtt30
      vnetb_etbtt33 = vnetb_etbtt32
      vnetb_etbtt35 = vnetb_etbtt34
      vnetb_etbtt37 = vnetb_etbtt36
      vnetb_etbtt39 = vnetb_etbtt38
      vnetb_etbtt41 = vnetb_etbtt40
      vnetb_etbtt43 = vnetb_etbtt42
      vnetb_etbtt45 = vnetb_etbtt44
      vnetb_etbtt47 = vnetb_etbtt46
      vnetb_etbtt49 = vnetb_etbtt48
      vnetb_etbtt51 = vnetb_etbtt50
      vnetb_etbtt53 = vnetb_etbtt52
      vnetb_etbtt55 = vnetb_etbtt54
      vnetb_etbtt57 = vnetb_etbtt56
      vnetb_etbtt59 = vnetb_etbtt58
      vnetb_etbtt61 = vnetb_etbtt60
      vnetb_etbtt63 = vnetb_etbtt62
      return
      end
