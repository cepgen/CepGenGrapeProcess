*     File ebb_ebbee/fncebb_ebbee.f : Sat Mar 18 19:45:05 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      function fncebb_ebbee(x)
      implicit real*8(a-h,o-z)
      real*8   fncebb_ebbee

      include 'inclebb_ebbee1.h'
      include 'inclk.inc'

      parameter ( mxdim = 50 )
      common / loop0 / loop
      common / bparm1 / xl(mxdim),xu(mxdim),ndim,nwild,
     &                 ig(mxdim),ncall
      common / bparm2 / acc1,acc2,itmx1,itmx2

      real*8   x(mxdim)

      real*8     ansm, ans

      real*8     xx(mxdim),p(4,nebb_ebbeeextn),pp(nebb_ebbeeextn,nebb_ebbeeextn)
      common /sp4vec/ vec(4,nebb_ebbeeextn)
*=======================================================================
*          initialization
*=======================================================================

      ansum = 0.0d0

      do  i = 1, ndim
         xx(i) = x(i)
      enddo

      nreg  = 1
      dft   = 0.d0
*=======================================================================
*          kinematics
*=======================================================================

      kreg = 0
      do 1000 ireg = 1 , mxreg

         if(ireg .gt. nreg) go to 1000

        call kinem_4f(nebb_ebbeeextn, xx, p, pp, yacob,nreg,ireg,jump)

*-----------------------------------------------------------------------
*          reset the temporal buffer for the region 1
*-----------------------------------------------------------------------
      if(ireg .eq. 1) then
          dft = 0.d0
          do k = 1, nebb_ebbeeextn
           do j = 1, 4
            vec(j,k) = 0.d0
           enddo
          enddo
      endif

         if(jump .ne. 0) go to 1000

*-----------------------------------------------------------------------
*          for user cut
*-----------------------------------------------------------------------

*--> call usrcut(jump)
*--> if(jump .ne. 0) goto 1000

*-----------------------------------------------------------------------
*          four momenta of external particles
*-----------------------------------------------------------------------

      do i = 1, 4
*     1: initial b-bar mass=amdq(3)
      peebb_ebbee1(i) = p(i, 1)

*     2: initial positron mass=amlp(1)
      peebb_ebbee2(i) = p(i, 2)

*     3: final b-bar mass=amdq(3)
      peebb_ebbee3(i) = p(i, 3)

*     4: final positron mass=amlp(1)
      peebb_ebbee4(i) = p(i, 4)

*     5: final electron mass=amlp(1)
      peebb_ebbee5(i) = p(i, 5)

*     6: final positron mass=amlp(1)
      peebb_ebbee6(i) = p(i, 6)

      enddo

*-----------------------------------------------------------------------
*          inner products of momenta of external particles
*-----------------------------------------------------------------------

      do j = 1, nebb_ebbeeextn
       do i = 1, nebb_ebbeeextn
        prod(i, j) = pp(i, j)
       enddo
      enddo

*=======================================================================
*          amplitude calculation
*=======================================================================
        call aebb_ebbeekptbl
*       ===================
        call kinemp_4f( ireg )
*       ===================
        call aebb_ebbeemptbl
*       ===================
        call aebb_ebbeempsum(ansm)
*       ===================

        fknorm = yacob*aspin
        ans    = ansm*fknorm/aident
        ansum  = ansum + ans

        call spstr4(kreg, nebb_ebbeeextn, lebb_ebbeeag, ans, p, agc)

*-----------------------------------------------------------------------
*          save four momenta and probabilities of the region 1
*-----------------------------------------------------------------------
      if(ireg .eq. 1) then
          dft = ans
          do k = 1, nebb_ebbeeextn
           do j = 1, 4
            vec(j,k) = p(j,k)
           enddo
          enddo
      endif

 1000 continue
      call spput4(kreg, nebb_ebbeeextn, lebb_ebbeeag, p, agc)

      fncebb_ebbee = ansum

      return
      end
