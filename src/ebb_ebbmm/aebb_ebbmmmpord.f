*     File ebb_ebbmm/aebb_ebbmmmpord.f : Sat Mar 18 19:45:06 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aebb_ebbmmmpord(lvt, avt, iord, at)
      implicit real*8(a-h,o-z)
      include 'inclebb_ebbmm1.h'
      include 'inclk.inc'
      integer    lvt(0:nebb_ebbmmextn), iord(nebb_ebbmmextn)
      complex*16 avt(0:lebb_ebbmmag-1), at(0:lebb_ebbmmag-1)
*
      integer j(nebb_ebbmmextn), jv(nebb_ebbmmextn)
*-----------------------------------------------------------------------
      ibas  = 1
      do 10 i = 1, nebb_ebbmmextn
        jv(iord(i)) = ibas
        ibas    = ibas*lvt(i)
        j(i)    = jhs(i)
   10 continue

      it = 0
  100 continue
        iv = 0
        do 110 i = 1, nebb_ebbmmextn
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
            if(ii.le.nebb_ebbmmextn) then
              goto 120
            else
              goto 190
            endif
          endif
        goto 100
  190 continue
      return
      end

