*     File ep_eXmm/aep_eXmmmclr.f : Sat Mar 18 19:45:45 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aep_eXmmmclr
      implicit real*8(a-h,o-z)

      include 'inclep_eXmm1.h'
      include 'inclk.inc'

           do 300 igr = 0, nep_eXmmgraph
             ansp(igr) = 0.0d0
  300      continue
           do 310 igr = 1, nep_eXmmgraph
             ancp(igr) = 0.0d0
  310      continue
           fkcall = 0
           nkcall = 0

      return
      end
