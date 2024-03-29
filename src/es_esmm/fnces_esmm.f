*     File es_esmm/fnces_esmm.f : Sat Mar 18 19:45:02 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      function fnces_esmm(x)
      implicit real*8(a-h,o-z)
      real*8   fnces_esmm

      include 'incles_esmm1.h'
      include 'inclk.inc'

      parameter ( mxdim = 50 )
      common / loop0 / loop
      common / bparm1 / xl(mxdim),xu(mxdim),ndim,nwild,
     &                 ig(mxdim),ncall
      common / bparm2 / acc1,acc2,itmx1,itmx2

      real*8   x(mxdim)

      real*8     ansm, ans

      real*8     xx(mxdim),p(4,nes_esmmextn),pp(nes_esmmextn,nes_esmmextn)
      common /sp4vec/ vec(4,nes_esmmextn)
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

        call kinem_4f(nes_esmmextn, xx, p, pp, yacob,nreg,ireg,jump)

*-----------------------------------------------------------------------
*          reset the temporal buffer for the region 1
*-----------------------------------------------------------------------
      if(ireg .eq. 1) then
          dft = 0.d0
          do k = 1, nes_esmmextn
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
*     1: initial s mass=amdq(2)
      pees_esmm1(i) = p(i, 1)

*     2: initial positron mass=amlp(1)
      pees_esmm2(i) = p(i, 2)

*     3: final s mass=amdq(2)
      pees_esmm3(i) = p(i, 3)

*     4: final positron mass=amlp(1)
      pees_esmm4(i) = p(i, 4)

*     5: final muon mass=amlp(2)
      pees_esmm5(i) = p(i, 5)

*     6: final anti-muon mass=amlp(2)
      pees_esmm6(i) = p(i, 6)

      enddo

*-----------------------------------------------------------------------
*          inner products of momenta of external particles
*-----------------------------------------------------------------------

      do j = 1, nes_esmmextn
       do i = 1, nes_esmmextn
        prod(i, j) = pp(i, j)
       enddo
      enddo

*=======================================================================
*          amplitude calculation
*=======================================================================
        call aes_esmmkptbl
*       ===================
        call kinemp_4f( ireg )
*       ===================
        call aes_esmmmptbl
*       ===================
        call aes_esmmmpsum(ansm)
*       ===================

        fknorm = yacob*aspin
        ans    = ansm*fknorm/aident
        ansum  = ansum + ans

        call spstr4(kreg, nes_esmmextn, les_esmmag, ans, p, agc)

*-----------------------------------------------------------------------
*          save four momenta and probabilities of the region 1
*-----------------------------------------------------------------------
      if(ireg .eq. 1) then
          dft = ans
          do k = 1, nes_esmmextn
           do j = 1, 4
            vec(j,k) = p(j,k)
           enddo
          enddo
      endif

 1000 continue
      call spput4(kreg, nes_esmmextn, les_esmmag, p, agc)

      fnces_esmm = ansum

      return
      end
