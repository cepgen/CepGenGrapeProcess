      subroutine spstr4(kreg, nextn, lag, ans, p, agc)
      implicit real*8 (a-h,o-z)
      integer  kreg, nextn, lag
      real*8   ans, p(4,nextn)
      complex*16  agc(0:lag-1)
      parameter (MAXREG = 6)
      parameter (mxextn = 10)
      parameter (mxlag  = 6561)
      complex*16  grcevt
      common /grc4sp/ answrk(0:MAXREG-1)
     &               ,wrkp(4, mxextn,   0:MAXREG-1)
     &               ,grcevt(0:mxlag-1, 0:MAXREG-1)
*-----------------------------------------------------------------------
      if (kreg .GT. MAXREG-1) then
        write(6,*) '!!!Error in spstr4!!!'
        write(6,*) '  ---> kreg(=',kreg,') > MAXREG-1(=',MAXREG-1,')'
        STOP
      endif
      if (nextn .GT. mxextn) then
        write(6,*) '!!!Error in spstr4!!!'
        write(6,*) '  ---> nextn(=',nextn,') > mxextn(=',mxextn,')'
        STOP
      endif
      if (lag .GT. mxlag-1) then
        write(6,*) '!!!Error in spstr4!!!'
        write(6,*) '  ---> lag(=',lag,') > mxlag-1(=',mxlag-1,')'
        STOP
      endif
      answrk(kreg) = ans
      do i = 1 , nextn
       do j = 1 , 4
        wrkp(j,i, kreg) = p(j,i)
       enddo
      enddo
      do ih = 0, lag-1
       grcevt(ih, kreg) = agc( ih )
      enddo
      return
      end
