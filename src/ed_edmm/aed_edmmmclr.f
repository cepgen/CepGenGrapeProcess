*     File ed_edmm/aed_edmmmclr.f : Sat Mar 18 19:45:00 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aed_edmmmclr
      implicit real*8(a-h,o-z)

      include 'incled_edmm1.h'
      include 'inclk.inc'

           do 300 igr = 0, ned_edmmgraph
             ansp(igr) = 0.0d0
  300      continue
           do 310 igr = 1, ned_edmmgraph
             ancp(igr) = 0.0d0
  310      continue
           fkcall = 0
           nkcall = 0

      return
      end
