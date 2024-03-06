      subroutine TAswap(A,B)
      implicit NONE
*------------ Argument ------------
      double precision  A,B
*----------------------------------
*------- Local variables -------
      double precision  Dtmp
*-------------------------------
      Dtmp = A
      A = B
      B = Dtmp
      return
      end
