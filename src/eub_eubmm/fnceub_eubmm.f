*     File eub_eubmm/fnceub_eubmm.f : Sat Mar 18 19:45:00 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      function fnceub_eubmm(x)
      implicit real*8(a-h,o-z)
      real*8   fnceub_eubmm

      include 'incleub_eubmm1.h'
      include 'inclk.inc'

      parameter ( mxdim = 50 )
      common / loop0 / loop
      common / bparm1 / xl(mxdim),xu(mxdim),ndim,nwild,
     &                 ig(mxdim),ncall
      common / bparm2 / acc1,acc2,itmx1,itmx2

      real*8   x(mxdim)

      real*8     ansm, ans

      real*8     xx(mxdim),p(4,neub_eubmmextn),pp(neub_eubmmextn,neub_eubmmextn)
      common /sp4vec/ vec(4,neub_eubmmextn)
*=======================================================================
*          initialization
*=======================================================================

      ansum = 0.0d0

      do i = 1, ndim
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

        call kinem_4f(neub_eubmmextn, xx, p, pp, yacob,nreg,ireg,jump)

*-----------------------------------------------------------------------
*          reset the temporal buffer for the region 1
*-----------------------------------------------------------------------
      if(ireg .eq. 1) then
          dft = 0.d0
          do k = 1, neub_eubmmextn
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
*     1: initial u-bar mass=amuq(1)
      peeub_eubmm1(i) = p(i, 1)

*     2: initial positron mass=amlp(1)
      peeub_eubmm2(i) = p(i, 2)

*     3: final u-bar mass=amuq(1)
      peeub_eubmm3(i) = p(i, 3)

*     4: final positron mass=amlp(1)
      peeub_eubmm4(i) = p(i, 4)

*     5: final muon mass=amlp(2)
      peeub_eubmm5(i) = p(i, 5)

*     6: final anti-muon mass=amlp(2)
      peeub_eubmm6(i) = p(i, 6)

      enddo

*-----------------------------------------------------------------------
*          inner products of momenta of external particles
*-----------------------------------------------------------------------

      do j = 1, neub_eubmmextn
       do i = 1, neub_eubmmextn
        prod(i, j) = pp(i, j)
       enddo
      enddo

*=======================================================================
*          amplitude calculation
*=======================================================================
        call aeub_eubmmkptbl
*       ===================
        call kinemp_4f( ireg )
*       ===================
        call aeub_eubmmmptbl
*       ===================
        call aeub_eubmmmpsum(ansm)
*       ===================

        fknorm = yacob*aspin
        ans    = ansm*fknorm/aident
        ansum  = ansum + ans

        call spstr4(kreg, neub_eubmmextn, leub_eubmmag, ans, p, agc)

*-----------------------------------------------------------------------
*          save four momenta and probabilities of the region 1
*-----------------------------------------------------------------------
      if(ireg .eq. 1) then
          dft = ans
          do k = 1, neub_eubmmextn
           do j = 1, 4
            vec(j,k) = p(j,k)
           enddo
          enddo
      endif

 1000 continue
      call spput4(kreg, neub_eubmmextn, leub_eubmmag, p, agc)

      fnceub_eubmm = ansum

      return
      end
