*     File eub_eubee/geub_eubeefinit.f : Sat Mar 18 19:44:59 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine geub_eubeefinit
      implicit real*8(a-h,o-z)

      include 'incleub_eubee1.h'
      include 'graepia.inc'
      include 'inclk.inc'
      common /chcntl/ jwidth, jtgamm
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
* Propagator = (i)**(-pphase)
      pphase = -1

*-----------------------------------------------------------------------
* Parameters for the process
      kmngr  = neub_eubeegraph
      kmnext = neub_eubeeextn
      kmlag  = leub_eubeeag

*-----------------------------------------------------------------------
* Gauge parametes (default is unitary gauge)
      igau00 = 0
      igauab = 0
      igauwb = 0
      igauzb = 0
      igaugb = 0
      agauge(igau00) = 1.0d0
      agauge(igauab) = 1.0d0
      agauge(igauwb) = 1.0d0
      agauge(igauzb) = 1.0d0
      agauge(igaugb) = 1.0d0

* Spin average
      aspin = 1.0d0

*     1: initial u-bar mass=amuq(1)
      jhs(1) = 0
      jhe(1) = lextrn - 1
      aspin = aspin/dble(jhe(1) - jhs(1)+1)

*     2: initial positron mass=amlp(1)
      jhs(2) = 0
      jhe(2) = lextrn - 1
      aspin = aspin/dble(jhe(2) - jhs(2)+1)

*     3: final u-bar mass=amuq(1)
      jhs(3) = 0
      jhe(3) = lextrn - 1

*     4: final positron mass=amlp(1)
      jhs(4) = 0
      jhe(4) = lextrn - 1

*     5: final electron mass=amlp(1)
      jhs(5) = 0
      jhe(5) = lextrn - 1

*     6: final positron mass=amlp(1)
      jhs(6) = 0
      jhe(6) = lextrn - 1

* Flag of cyclic polarization

* Initial color average
      acolav = 1.0d0/3.0d0
      aspin = aspin*acolav

* Symmetry factor for final identical particles
      aident = 1.0d0
*     aident = 1.0d0
      aspin = aspin/aident

* graph selection
      do 10 n1 = 1, neub_eubeegraph
        i = int( (n1-1)/30 ) + 1
        ibit = n1 - 30*(i-1) - 1
        jselg(n1) = Ibtest(jgra_flag(i),ibit)
        if (print_flag .EQ. 1) then
          write(6,*)  n1, ' :', jselg(n1)
        endif
   10 continue

      jgraph = 0


*     Graph selection flags
      jgluon = 1

      jhiggs = 1

*     Running width (0) or fixed width(1) in CHANEL
      jwidth = 1

      jtgamm = 1

* weight for the color base
      cfmtx = 3.d0


* Color string information
      kmcbas    = 1
      kmcbmx    = 2
      kmcstr(1) = 1
      kmcstr(2) = 3

*-----------------------------------------------------------------------
* Particle name
      kmprtc(1) = 'u-bar'
      kmprtl(1) = 5
      kmprtc(2) = 'positron'
      kmprtl(2) = 8
      kmprtc(3) = 'u-bar'
      kmprtl(3) = 5
      kmprtc(4) = 'positron'
      kmprtl(4) = 8
      kmprtc(5) = 'electron'
      kmprtl(5) = 8
      kmprtc(6) = 'positron'
      kmprtl(6) = 8

* Masses of external particles
      amass1(1) = amuq(1)
      amass1(2) = amlp(1)
      amass1(3) = amuq(1)
      amass1(4) = amlp(1)
      amass1(5) = amlp(1)
      amass1(6) = amlp(1)
      amass2(1) = amass1(1)**2
      amass2(2) = amass1(2)**2
      amass2(3) = amass1(3)**2
      amass2(4) = amass1(4)**2
      amass2(5) = amass1(5)**2
      amass2(6) = amass1(6)**2

* Charge*3
      kcharg(1) =  -2
      kcharg(2) =   3
      kcharg(3) =  -2
      kcharg(4) =   3
      kcharg(5) =  -3
      kcharg(6) =   3

* KFcode
      kfcode(1) =  -2
      kfcode(2) = -11
      kfcode(3) =  -2
      kfcode(4) = -11
      kfcode(5) =  11
      kfcode(6) = -11

      return
      end
