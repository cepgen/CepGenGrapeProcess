*     File edb_edbmm/fncedb_edbmm.f : Sat Mar 18 19:45:01 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      function fncedb_edbmm(x)
      implicit real*8(a-h,o-z)
      real*8   fncedb_edbmm

      include 'incledb_edbmm1.h'
      include 'inclk.inc'

      parameter ( mxdim = 50 )
      common / loop0 / loop
      common / bparm1 / xl(mxdim),xu(mxdim),ndim,nwild,
     &                 ig(mxdim),ncall
      common / bparm2 / acc1,acc2,itmx1,itmx2

      real*8   x(mxdim)

      real*8     ansm, ans

      real*8     xx(mxdim),p(4,nedb_edbmmextn),pp(nedb_edbmmextn,nedb_edbmmextn)
      common /sp4vec/ vec(4,nedb_edbmmextn)
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

        call kinem_4f(nedb_edbmmextn, xx, p, pp, yacob,nreg,ireg,jump)

*-----------------------------------------------------------------------
*          reset the temporal buffer for the region 1
*-----------------------------------------------------------------------
      if(ireg .eq. 1) then
          dft = 0.d0
          do k = 1, nedb_edbmmextn
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
*     1: initial d-bar mass=amdq(1)
      peedb_edbmm1(i) = p(i, 1)

*     2: initial positron mass=amlp(1)
      peedb_edbmm2(i) = p(i, 2)

*     3: final d-bar mass=amdq(1)
      peedb_edbmm3(i) = p(i, 3)

*     4: final positron mass=amlp(1)
      peedb_edbmm4(i) = p(i, 4)

*     5: final muon mass=amlp(2)
      peedb_edbmm5(i) = p(i, 5)

*     6: final anti-muon mass=amlp(2)
      peedb_edbmm6(i) = p(i, 6)

      enddo

*-----------------------------------------------------------------------
*          inner products of momenta of external particles
*-----------------------------------------------------------------------

      do j = 1, nedb_edbmmextn
       do i = 1, nedb_edbmmextn
        prod(i, j) = pp(i, j)
       enddo
      enddo

*=======================================================================
*          amplitude calculation
*=======================================================================
        call aedb_edbmmkptbl
*       ===================
        call kinemp_4f( ireg )
*       ===================
        call aedb_edbmmmptbl
*       ===================
        call aedb_edbmmmpsum(ansm)
*       ===================

        fknorm = yacob*aspin
        ans    = ansm*fknorm/aident
        ansum  = ansum + ans

        call spstr4(kreg, nedb_edbmmextn, ledb_edbmmag, ans, p, agc)

*-----------------------------------------------------------------------
*          save four momenta and probabilities of the region 1
*-----------------------------------------------------------------------
      if(ireg .eq. 1) then
          dft = ans
          do k = 1, nedb_edbmmextn
           do j = 1, 4
            vec(j,k) = p(j,k)
           enddo
          enddo
      endif

 1000 continue
      call spput4(kreg, nedb_edbmmextn, ledb_edbmmeag, p, agc)

      fncedb_edbmm = ansum

      return
      end
