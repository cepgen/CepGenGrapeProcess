*     File eu_euee/aeu_eueempsum.f : Sat Mar 18 19:44:59 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
*     CF version 
************************************************************************
      subroutine aeu_eueempsum(ans)
      implicit real*8(a-h,o-z)
      include 'incleu_euee1.h'
      include 'inclk.inc'
      real*8     ans
      complex*16 anc
*-----------------------------------------------------------------------
      jsum = 1
      do 10 ie = 1, neu_eueeextn
        jsum = jsum*(jhe(ie) - jhs(ie) + 1)
   10 continue

      anc = 0

*     sum over helicities
      call AgcPol(agc)
      do 100 ind = 0, jsum - 1

*       sum over color bases
        anc = anc  + agc(ind)*conjg(agc(ind))
  100 continue
      ans = cfmtx*anc

      return
      end

