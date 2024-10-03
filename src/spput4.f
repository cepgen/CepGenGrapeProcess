      subroutine spput4(nreg, nextn, lag, p, agc)
      implicit real*8 (a-h,o-z)
      integer  nreg, nextn, lag
      real*8   p(4,nextn)
      complex*16  agc(0:lag-1)
      parameter (MAXREG = 6)
      parameter (mxextn = 10)
      parameter (mxlag  = 6561)
      complex*16  grcevt
      common /grc4sp/ answrk(0:MAXREG-1)
     &               ,wrkp(4, mxextn,   0:MAXREG-1)
     &               ,grcevt(0:mxlag-1, 0:MAXREG-1)
      dimension  cratio(0:MAXREG-1)
      common /sp4vec/ vec(4,mxextn)
      real*8    CepGen_random_number
      external  CepGen_random_number
*-----------------------------------------------------------------------
* ------- selection ------------------
      ireg = 0
      if( nreg .gt. 1 ) then
        allsum = 0.0d0
        do i = 0, nreg-1
           allsum = allsum + answrk(i)
        enddo
        tmpsum = 0.0d0
        do i = 0, nreg-1
           tmpsum = tmpsum + answrk(i)
           cratio(i) = tmpsum/allsum
        enddo
        cran = CepGen_random_number()
        do 3000 i = 0, nreg-1
           if( cratio(i) .gt. cran ) then
               ireg = i
               goto 3000
           endif
 3000   continue
      endif
      do i = 1 , nextn
       do j = 1 , 4
        p(j,i) = wrkp(j,i, ireg)
       enddo
      enddo
      do i = 1 , nextn
       do j = 1 , 4
        vec(j,i) = p(j,i)
       enddo
      enddo
      do ih = 0, lag-1
         agc(ih) = grcevt(ih, ireg)
      enddo
      return
      end
