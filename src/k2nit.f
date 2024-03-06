C##########################################################################
      SUBROUTINE K2NIT
C---------------------------------------------------------------------
C   GRACE System Library File
C   KINEM No. : 4009
C   Date      : 1996.04.14  ; Modified from 4007 to fit for two-photon
C                            processes.
C             : 1996.07.23  ; To fit SF.
C   Author    : Y.Kurihara
C---------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'inclk.inc'
      common/kmcntl_4f/iresns(4),icos3,icosq3,icos5,isr,iswap,ident
     &                ,iphi6,ieeee,i34,itag,isym
      data iprint/0/
*-----------------------------------------------------------------------
      COMMON/KINEM1/S(3),W(3),FACT
      COMMON/CUT001_4f/COSCUT(2,4),ENGYCT(2,4),AMASCT(2,6),ARESNS(2,4)
     .,opncut,swapm2
*-----------------------------------------------------------------------
      include 'graepia.inc'
      integer         jproc
       common /amjprc/jproc
*-----------------------------------------------------------------------
      cstag1=cos(0.5d0*rad)
      cstag2=cos( 4.d0*rad)
      cstag2=cos( 6.d0*rad)
      if     (coscut(2,1).ge. cstag1        .and.
     .        coscut(1,2).gt.-cstag1 ) then
       itag=1
       isym=0
       i34=3
      else if(coscut(2,1).lt. cstag1        .and.
     .        coscut(1,2).le.-cstag1 ) then
       itag=1
       isym=0
       i34=4
      else if(coscut(2,1).gt. cstag1        .and.
     .        coscut(1,2).lt.-cstag1 ) then
       itag=0
       isym=1
       i34=-999
      else
       if(coscut(2,1).lt. cstag2        .and.
     .    coscut(1,2).gt.-cstag2 ) then
        itag=2
       else
        itag=1
       end if
        isym=0
       if     (coscut(2,1).eq.-coscut(1,2)) then
        isym=1
        i34=-999
       else if(coscut(2,1).gt.-coscut(1,2)) then
        i34=3
       else
        i34=4
       end if
      end if
      isym = isym_34      !!!!!!
      i34  = ii34         !!!!!!
      write(6,*)'   itag,isym,i34 =',itag,isym,i34
      if (Neps_p.LT.1 .OR. Neps_p.GT.100) then
         write(6,*) '!!!Error in KINIT_4f!!!'
         write(6,*) '  ---> Invalid NEPSP =', Neps_p
         write(6,*) '  ---> It should be in [1,100].'
         write(6,*) '  ---> Good-bye!'
         STOP
      endif
* << Particle 5-6 >>
      IRESNS(2)  = Ireso56
      if (IRESNS(2).GE.1) then
        ARESNS(1,2)=amz
        ARESNS(2,2)=agz
      endif
* << Particle 4-5-6 >>
      IRESNS(3)  = Ireso456
      if (IRESNS(3).GE.1) then
        ARESNS(1,3)=Mas_reso456     ! If you want to treat narrow resonance,
        ARESNS(2,3)=Wid_reso456     ! set resonance mass and width.
      endif
*     -------
      icos3 = IcosP3
*     -------
      icos5 = IcosP5
*     -------
      iphi6=0
*     -------
      icosq3 = 2
*     ----------
      ident = 0
*     ---------
      write(6,*) '   Ireso56,Ireso456 =', IRESNS(2),IRESNS(3)
      if (Rnn_MJ56 .LE. 0) then
        write(6,*) '!!!Error in KINIT_4f!!!'
        write(6,*) '  ---> NNMJ56 =', Rnn_MJ56
        write(6,*) '  ---> It should be >0.'
        write(6,*) '  ---> Good-bye!'
        STOP
      endif
      Inn =  int(Rnn_MJ56) +1
      Dnn = dble(Rnn_MJ56) +1D0
      Diff_nn = Dnn - dble(Inn)
      if (Rnn_456 .LE. 1.) then
        write(6,*) '!!!Error in KINIT_4f!!!'
        write(6,*) '  ---> RNN456 =', Rnn_456
        write(6,*) '  ---> It should be >1.0.'
        write(6,*) '  ---> Good-bye!'
        STOP
      endif
*################# For Multi-Jacobian method #################
      if (IRESNS(2) .NE. 40)  GOTO 1999
*>>> QED
      weight_mj(0) = weimj_QED
      num_reso_mj = 0
*>>> Z^0
      if ( .false.
     +     .OR. jgra_sel.EQ. 4
     +     .OR. jgra_sel.EQ. 5
     +     .OR. jgra_sel.EQ.14
     +   ) then
         KF_mj = 23
         num_reso_mj = num_reso_mj + 1
         if (num_reso_mj .GT. max_reso_mj)  GOTO 1100
         call AMG_RESO( KF_mj
     &                 ,mass_mj(num_reso_mj), width_mj(num_reso_mj) )
         weight_mj(num_reso_mj)  = weimj_Z0   !!!
            if ( .false.
     +         .OR. mass_mj(num_reso_mj)  .LE.0
     +         .OR. width_mj(num_reso_mj) .LE.0
     +         .OR. weight_mj(num_reso_mj).LE.0
     +         )  GOTO 1200
         mass2_mj(num_reso_mj) = mass_mj(num_reso_mj)**2
      endif
*>>> Normalization of weights
      if (num_reso_mj .GT. 0) then
        sumsum = 0D0
        do i = 0, num_reso_mj
          sumsum = sumsum + weight_mj(i)
        enddo
        sum = 0D0
        do i = 0, num_reso_mj
          sum = sum + weight_mj(i)
          weisel_mj(i) = sum/sumsum
          weight_mj(i) = weight_mj(i)/sumsum
        enddo
        write(6,*) ' '
        write(6,*) 'MJ>> # of resonance particles =', num_reso_mj
        write(6,*) 'MJ>> Masses  ='
     &                ,(real(mass_mj(i)),i=1,num_reso_mj)
        write(6,*) 'MJ>> Widths  ='
     &                ,(real(width_mj(i)),i=1,num_reso_mj)
        write(6,*) 'MJ>> Weights ='
     &                ,(real(weight_mj(i)),i=0,num_reso_mj)
        write(6,*) 'MJ>> WeiSELs ='
     &                ,(real(weisel_mj(i)),i=0,num_reso_mj)
        write(6,*) ' '
      endif
      GOTO 1999
 1100 continue
      write(6,*) '!!!Error in KINIT_4f!!!'
      write(6,*) '  ---> num_reso_mj(=', num_reso_mj, ') exceeds'
     &            //' max_reso_mj(=', max_reso_mj, ').'
      write(6,*) '  ---> Please inform the author!'
      write(6,*) '  ---> Good-bye!'
      STOP
 1200 continue
      write(6,*) '!!!Error in KINIT_4f!!!'
      write(6,*) '  ---> Check parameters for resonance particles.'
      write(6,*) '  ---> num_reso_mj =', num_reso_mj
      write(6,*) '  ---> Mass   =', mass_mj(num_reso_mj)
      write(6,*) '  ---> Width  =', width_mj(num_reso_mj)
      write(6,*) '  ---> Weight =', weight_mj(num_reso_mj)
      write(6,*) '  ---> Good-bye!'
      STOP
 1999 continue
*#############################################################
      RETURN
      END
