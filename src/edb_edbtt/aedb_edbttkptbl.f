*     File edb_edbtt/aedb_edbttkptbl.f : Sat Mar 18 19:45:01 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aedb_edbttkptbl
      implicit real*8(a-h,o-z)

      include 'incledb_edbtt1.h'
      include 'inclk.inc'
      include 'incledb_edbttp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfedb_edbtt2(i) =  + peedb_edbtt1(i)
        pfedb_edbtt4(i) =  + peedb_edbtt2(i)
        pfedb_edbtt6(i) =  + peedb_edbtt1(i) + peedb_edbtt2(i)
        pfedb_edbtt8(i) =  - peedb_edbtt3(i)
        pfedb_edbtt10(i) =  + peedb_edbtt1(i) - peedb_edbtt3(i)
        pfedb_edbtt12(i) =  + peedb_edbtt2(i) - peedb_edbtt3(i)
        pfedb_edbtt14(i) =  + peedb_edbtt1(i) + peedb_edbtt2(i) - peedb_edbtt3(i)
        pfedb_edbtt16(i) =  - peedb_edbtt4(i)
        pfedb_edbtt18(i) =  + peedb_edbtt1(i) - peedb_edbtt4(i)
        pfedb_edbtt20(i) =  + peedb_edbtt2(i) - peedb_edbtt4(i)
        pfedb_edbtt22(i) =  + peedb_edbtt1(i) + peedb_edbtt2(i) - peedb_edbtt4(i)
        pfedb_edbtt24(i) =  - peedb_edbtt3(i) - peedb_edbtt4(i)
        pfedb_edbtt26(i) =  + peedb_edbtt1(i) - peedb_edbtt3(i) - peedb_edbtt4(i)
        pfedb_edbtt28(i) =  + peedb_edbtt2(i) - peedb_edbtt3(i) - peedb_edbtt4(i)
        pfedb_edbtt30(i) =  + peedb_edbtt5(i) + peedb_edbtt6(i)
        pfedb_edbtt32(i) =  - peedb_edbtt5(i)
        pfedb_edbtt34(i) =  + peedb_edbtt1(i) - peedb_edbtt5(i)
        pfedb_edbtt36(i) =  + peedb_edbtt2(i) - peedb_edbtt5(i)
        pfedb_edbtt38(i) =  + peedb_edbtt1(i) + peedb_edbtt2(i) - peedb_edbtt5(i)
        pfedb_edbtt40(i) =  - peedb_edbtt3(i) - peedb_edbtt5(i)
        pfedb_edbtt42(i) =  + peedb_edbtt1(i) - peedb_edbtt3(i) - peedb_edbtt5(i)
        pfedb_edbtt44(i) =  + peedb_edbtt2(i) - peedb_edbtt3(i) - peedb_edbtt5(i)
        pfedb_edbtt46(i) =  + peedb_edbtt4(i) + peedb_edbtt6(i)
        pfedb_edbtt48(i) =  - peedb_edbtt4(i) - peedb_edbtt5(i)
        pfedb_edbtt50(i) =  + peedb_edbtt1(i) - peedb_edbtt4(i) - peedb_edbtt5(i)
        pfedb_edbtt52(i) =  + peedb_edbtt2(i) - peedb_edbtt4(i) - peedb_edbtt5(i)
        pfedb_edbtt54(i) =  + peedb_edbtt3(i) + peedb_edbtt6(i)
        pfedb_edbtt56(i) =  - peedb_edbtt3(i) - peedb_edbtt4(i) - peedb_edbtt5(i)
        pfedb_edbtt58(i) =  - peedb_edbtt2(i) + peedb_edbtt6(i)
        pfedb_edbtt60(i) =  - peedb_edbtt1(i) + peedb_edbtt6(i)
        pfedb_edbtt62(i) =  + peedb_edbtt6(i)

        pfedb_edbtt3(i) = - pfedb_edbtt2(i)
        pfedb_edbtt5(i) = - pfedb_edbtt4(i)
        pfedb_edbtt7(i) = - pfedb_edbtt6(i)
        pfedb_edbtt9(i) = - pfedb_edbtt8(i)
        pfedb_edbtt11(i) = - pfedb_edbtt10(i)
        pfedb_edbtt13(i) = - pfedb_edbtt12(i)
        pfedb_edbtt15(i) = - pfedb_edbtt14(i)
        pfedb_edbtt17(i) = - pfedb_edbtt16(i)
        pfedb_edbtt19(i) = - pfedb_edbtt18(i)
        pfedb_edbtt21(i) = - pfedb_edbtt20(i)
        pfedb_edbtt23(i) = - pfedb_edbtt22(i)
        pfedb_edbtt25(i) = - pfedb_edbtt24(i)
        pfedb_edbtt27(i) = - pfedb_edbtt26(i)
        pfedb_edbtt29(i) = - pfedb_edbtt28(i)
        pfedb_edbtt31(i) = - pfedb_edbtt30(i)
        pfedb_edbtt33(i) = - pfedb_edbtt32(i)
        pfedb_edbtt35(i) = - pfedb_edbtt34(i)
        pfedb_edbtt37(i) = - pfedb_edbtt36(i)
        pfedb_edbtt39(i) = - pfedb_edbtt38(i)
        pfedb_edbtt41(i) = - pfedb_edbtt40(i)
        pfedb_edbtt43(i) = - pfedb_edbtt42(i)
        pfedb_edbtt45(i) = - pfedb_edbtt44(i)
        pfedb_edbtt47(i) = - pfedb_edbtt46(i)
        pfedb_edbtt49(i) = - pfedb_edbtt48(i)
        pfedb_edbtt51(i) = - pfedb_edbtt50(i)
        pfedb_edbtt53(i) = - pfedb_edbtt52(i)
        pfedb_edbtt55(i) = - pfedb_edbtt54(i)
        pfedb_edbtt57(i) = - pfedb_edbtt56(i)
        pfedb_edbtt59(i) = - pfedb_edbtt58(i)
        pfedb_edbtt61(i) = - pfedb_edbtt60(i)
        pfedb_edbtt63(i) = - pfedb_edbtt62(i)
  100 continue

      vnedb_edbtt2 =  + amass2(1)
      vnedb_edbtt4 =  + amass2(2)
      vnedb_edbtt6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vnedb_edbtt8 =  + amass2(3)
      vnedb_edbtt10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vnedb_edbtt12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vnedb_edbtt14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vnedb_edbtt16 =  + amass2(4)
      vnedb_edbtt18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vnedb_edbtt20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vnedb_edbtt22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vnedb_edbtt24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vnedb_edbtt26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vnedb_edbtt28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vnedb_edbtt30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vnedb_edbtt32 =  + amass2(5)
      vnedb_edbtt34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vnedb_edbtt36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vnedb_edbtt38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vnedb_edbtt40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vnedb_edbtt42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vnedb_edbtt44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vnedb_edbtt46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vnedb_edbtt48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vnedb_edbtt50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vnedb_edbtt52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vnedb_edbtt54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vnedb_edbtt56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vnedb_edbtt58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vnedb_edbtt60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vnedb_edbtt62 =  + amass2(6)
      vnedb_edbtt3 = vnedb_edbtt2
      vnedb_edbtt5 = vnedb_edbtt4
      vnedb_edbtt7 = vnedb_edbtt6
      vnedb_edbtt9 = vnedb_edbtt8
      vnedb_edbtt11 = vnedb_edbtt10
      vnedb_edbtt13 = vnedb_edbtt12
      vnedb_edbtt15 = vnedb_edbtt14
      vnedb_edbtt17 = vnedb_edbtt16
      vnedb_edbtt19 = vnedb_edbtt18
      vnedb_edbtt21 = vnedb_edbtt20
      vnedb_edbtt23 = vnedb_edbtt22
      vnedb_edbtt25 = vnedb_edbtt24
      vnedb_edbtt27 = vnedb_edbtt26
      vnedb_edbtt29 = vnedb_edbtt28
      vnedb_edbtt31 = vnedb_edbtt30
      vnedb_edbtt33 = vnedb_edbtt32
      vnedb_edbtt35 = vnedb_edbtt34
      vnedb_edbtt37 = vnedb_edbtt36
      vnedb_edbtt39 = vnedb_edbtt38
      vnedb_edbtt41 = vnedb_edbtt40
      vnedb_edbtt43 = vnedb_edbtt42
      vnedb_edbtt45 = vnedb_edbtt44
      vnedb_edbtt47 = vnedb_edbtt46
      vnedb_edbtt49 = vnedb_edbtt48
      vnedb_edbtt51 = vnedb_edbtt50
      vnedb_edbtt53 = vnedb_edbtt52
      vnedb_edbtt55 = vnedb_edbtt54
      vnedb_edbtt57 = vnedb_edbtt56
      vnedb_edbtt59 = vnedb_edbtt58
      vnedb_edbtt61 = vnedb_edbtt60
      vnedb_edbtt63 = vnedb_edbtt62
      return
      end
