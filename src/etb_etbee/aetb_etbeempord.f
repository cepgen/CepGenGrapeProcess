*     File etb_etbee/aetb_etbeempord.f : Sat Mar 18 19:45:07 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aetb_etbeempord(lvt, avt, iord, at)
      implicit real*8(a-h,o-z)
      include 'incletb_etbee1.h'
      include 'inclk.inc'
      integer    lvt(0:netb_etbeeextn), iord(netb_etbeeextn)
      complex*16 avt(0:letb_etbeeag-1), at(0:letb_etbeeag-1)
*
      integer j(netb_etbeeextn), jv(netb_etbeeextn)
*-----------------------------------------------------------------------
      ibas  = 1
      do 10 i = 1, netb_etbeeextn
        jv(iord(i)) = ibas
        ibas    = ibas*lvt(i)
        j(i)    = jhs(i)
   10 continue

      it = 0
  100 continue
        iv = 0
        do 110 i = 1, netb_etbeeextn
          iv = iv +   jv(i)*j(i)
  110   continue

        at(it) = avt(iv)
        it = it + 1

        ii = 1
  120   continue
          j(ii) = j(ii) + 1
          if(j(ii).gt.jhe(ii)) then
            j(ii) = jhs(ii)
            ii = ii + 1
            if(ii.le.netb_etbeeextn) then
              goto 120
            else
              goto 190
            endif
          endif
        goto 100
  190 continue
      return
      end

