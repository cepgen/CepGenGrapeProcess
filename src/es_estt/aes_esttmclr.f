*     File es_estt/aes_esttmclr.f : Sat Mar 18 19:45:02 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aes_esttmclr
      implicit real*8(a-h,o-z)

      include 'incles_estt1.h'
      include 'inclk.inc'

           do 300 igr = 0, nes_esttgraph
             ansp(igr) = 0.0d0
  300      continue
           do 310 igr = 1, nes_esttgraph
             ancp(igr) = 0.0d0
  310      continue
           fkcall = 0
           nkcall = 0

      return
      end
