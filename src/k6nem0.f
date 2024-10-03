      SUBROUTINE K6NEM0(NEXTRN,X,PE,PP,YACOB,NREG,IREG,JUMP)
      IMPLICIT REAL* 8(A-H,O-Z)
************************************************************************
      INTEGER NEXTRN
      PARAMETER ( MXDIM = 50 )
      COMMON / BPARM1 / XL(MXDIM),XU(MXDIM),NDIM,NWILD,
     &                 IG(MXDIM),NCALL
      REAL*8  X(MXDIM)
      REAL*8  PE(4,NEXTRN), PP(NEXTRN,NEXTRN)
      REAL*8  YACOB
      INTEGER NREG, IREG
      INTEGER JUMP
      INCLUDE 'inclk.inc'
      common/ktest/qx56
      common/atest/ttt1,uuu1,ttt2,uuu2,tttt,uuuu,sss2,sss1,q56,q12,q22
     .,d1,g1,t1,u1,t2,u2
      common/kmcntl_4f/iresns(4),icos3,icosq3,icos5,isr,iswap,ident
     &                ,iphi6,ieeee,i34,itag,isym
*--------------------------------------------------------------------
      COMMON/KINEM1/S(3),W(3),FACT
      COMMON/CUT001_4f/COSCUT(2,4),ENGYCT(2,4),AMASCT(2,6),ARESNS(2,4)
     .,opncut,swapm2
      REAL*8 PBOST1(4),PBOST2(4),PK1(4),TT1(4),TT2(4),PDUM(4)
      DATA PDUM/4*1.D0/
      common/TAISR/yy2,y2max
*-----------------------------------------------------------------------
      data nev/0/
*-----------------------------------------------------------------------
      double precision  P1_lab,P2_lab, E1_lab,E2_lab
     &                 ,Pcms_lab(3),Ecms_lab(3)
     &                 ,GAMMAcms_lab(3),BETGAMcms_lab(3)
     &                 ,vec_isr(4)
       common /GEP_LAB/ P1_lab,P2_lab, E1_lab,E2_lab
     &                 ,Pcms_lab,Ecms_lab
     &                 ,GAMMAcms_lab,   BETGAMcms_lab
     &                 ,vec_isr
*-----------------------------------------------------------------------
      include 'graepia.inc'
      double precision  CepGen_random_number, Flux_factor
       external         CepGen_random_number, Flux_factor
*-----------------------------------------------------------------------
      JUMP   = 0
      YACOB  = 0.D0
      AJACOB = 1
      ipoint = 1
      csut11=coscut(1,1)
      csut21=coscut(2,1)
      csut12=coscut(1,2)
      csut22=coscut(2,2)
      csut13=coscut(1,3)
      csut23=coscut(2,3)
      csut14=coscut(1,4)
      csut24=coscut(2,4)
      if(isr.eq.1) then
       nn=NN_ISR
       xx2min = ((amass1(3)+amass1(4)+amass1(5)+amass1(6))**2-amass2(1))
     &          /(s(2)-amass2(1))
       y2max  = 1.d0-xx2min
       yy2    = y2max * X(NDIM-ipol_Ebeam)**nn
       if (yy2.lt.1D-300) yy2=1D-300
       S(3)    = (S(2)-amass2(1))*(1.d0-yy2) +amass2(1)
       W(3)    = sqrt(  max( S(3), 0.D0 )  )
       if ( W(3) .LE. (amass1(1)+amass1(2)) )  GOTO 999
       vec_isr(1) = 0.D0
       vec_isr(2) = 0.D0
       vec_isr(3) = -(S(2)-amass2(1))/2.D0/W(2) *yy2
       vec_isr(4) = abs(vec_isr(3))
       Pz_cms = -vec_isr(3)
       E_cms  = sqrt(  vec_isr(4)*vec_isr(4)  +  S(3)  )
       Pcms_lab(3) = GAMMAcms_lab(2)*Pz_cms + BETGAMcms_lab(2)*E_cms
       Ecms_lab(3) = GAMMAcms_lab(2)*E_cms + BETGAMcms_lab(2)*Pz_cms
       GAMMAcms_lab(3)  = Ecms_lab(3)/W(3)
       BETGAMcms_lab(3) = Pcms_lab(3)/W(3)
       Pz_cms = vec_isr(3)
       E_cms  = vec_isr(4)
       vec_isr(3) = GAMMAcms_lab(2)*Pz_cms + BETGAMcms_lab(2)*E_cms
       vec_isr(4) = abs(vec_isr(3))
cccccccc
      else
       S(3) = S(2)
       W(3) = W(2)
       radi=1
       do i=1,4
        vec_isr(i) = 0.D0
       enddo
       Pcms_lab(3) = Pcms_lab(2)
       Ecms_lab(3) = Ecms_lab(2)
       GAMMAcms_lab(3)  = GAMMAcms_lab(2)
       BETGAMcms_lab(3) = BETGAMcms_lab(2)
      endif
cccccccc
      ipoint = 2
      E1     = (S(3)+AMASS2(1)-AMASS2(2))/2.D0/W(3)
      P1     = SQRT( max( E1-AMASS1(1), 0.D0 ) * (E1+AMASS1(1)) )
      E2     = (S(3)+AMASS2(2)-AMASS2(1))/2.D0/W(3)
      P2     = P1
*Particle 1
      PE(1,1) = 0
      PE(2,1) = 0
      PE(3,1) = P1
      PE(4,1) = E1
*Particle 2
      PE(1,2) = 0
      PE(2,2) = 0
      PE(3,2) =-P2
      PE(4,2) = E2
      IF(ISYM.EQ.1) THEN
       IF(X(4).LT.0.5D0) THEN
        X4=X(4)*2.D0
        I34=3
       ELsE
        X4=2.D0-X(4)*2.D0
        I34=4
       END IF
       AJACOB=AJACOB*2.D0
      ELSE
       X4=X(4)
      END IF
***************
*    QX56     *
*    X=3 or 4 *
***************
       QX56MN=(     AMASS1(7-I34)+AMASS1(5)+AMASS1(6))**2
       QX56MX=(W(3)-AMASS1(  I34)                    )**2
       IF (QX56MX.LE.QX56MN) GOTO 999
       IF (QX56MN.LE.0)      GOTO 999
       IF     (IRESNS(3).EQ. 0) THEN
        QX56=QX56MN+(QX56MX-QX56MN)*X(1)
        AJACOB=AJACOB*(QX56MX-QX56MN)
       ELSE IF(IRESNS(3).EQ. 1) THEN
        ZM    =ARESNS(1,3)
        ZM2   =ZM*ZM
        ZMG   =ARESNS(1,3)*ARESNS(2,3)
        THEMIN=ATAN((QX56MN-ZM2)/ZMG)
        THEMAX=ATAN((QX56MX-ZM2)/ZMG)
        THE   =THEMIN+(THEMAX-THEMIN)*X(1)
        QX56  =ZMG*TAN(THE)+ZM2
        AJACOB=AJACOB*(THEMAX-THEMIN)*((QX56-ZM2)**2+ZMG**2)/ZMG
       ELSE IF(IRESNS(3).EQ.-1) THEN
        IF(QX56MN.LE.0) GOTO 999
        QX56=QX56MN*(QX56MX/QX56MN)**X(1)
        AJACOB=AJACOB*QX56*LOG(QX56MX/QX56MN)
       ELSE IF(IRESNS(3).EQ.-2) THEN
        IF(QX56MN.LE.0) GOTO 999
        t_min = sqrt(QX56MN)
        t_max = sqrt(QX56MX)
        t_diff = t_max - t_min
        t_val = t_min + t_diff * X(1)
        QX56 = t_val**2
        AJACOB=AJACOB* 2.D0*t_val*t_diff
       ELSE IF(IRESNS(3).EQ. 40) THEN
         DD = QX56MX - QX56MN
         QX56 = QX56MN + DD*X(1)**Rnn_456
             AJACOB = AJACOB *DD *Rnn_456 *X(1)**(Rnn_456-1D0)
       ELSE IF(IRESNS(3).EQ. 2) THEN
        IF(X(1).LT.0.5D0) THEN
         X1=X(1)*2.D0
         AJACOB=AJACOB*2.D0
         IF(QX56MN.LE.0) GOTO 999
          QX56MX=(QX56MX+QX56MN)/2.D0
          QX56=QX56MN*(QX56MX/QX56MN)**X1
          AJACOB=AJACOB*QX56*LOG(QX56MX/QX56MN)
         ELSE
          X1=(1.D0-X(1))*2.D0
          AJACOB=AJACOB*2.D0
          QX56MN=(QX56MX+QX56MN)/2.D0
          QX56=QX56MN+(QX56MX-QX56MN)*X1
          AJACOB=AJACOB*(QX56MX-QX56MN)
         END IF
        ELSE
         WRITE(6,*)' IRESNS(3) =',IRESNS(3),' is NOT supported. '
         stop
        END IF
        AX56=QX56
        EX   = (S(3)+AMASS2(I34)-QX56)/2.D0/W(3)
        PX   = SQRT( max( EX-AMASS1(I34), 0.D0 ) * (EX+AMASS1(I34)) )
        EX56 = (S(3)+QX56-AMASS2(I34))/2.D0/W(3)
        DE=(AMASS2(I34-2)-AMASS2(5-I34)+QX56-AMASS2(I34))/2.D0/W(3)
        DP=(DE*(EX+PE(4,I34-2))-AMASS2(I34-2)+AMASS2(I34))/(PX+ABS(PE(3,I34-2)))
***************
*    COSX     *
*    X=3 or 4 *
***************
      ipoint = 3
       IF     (ICOS3.EQ.0) THEN
        COSX = COSCUT(1,I34-2)+(COSCUT(2,I34-2)-COSCUT(1,I34-2))*X(3)
        AJACOB=AJACOB*(COSCUT(2,I34-2)-COSCUT(1,I34-2))
       ELSE IF(ICOS3.EQ.2) THEN
        IF(I34.EQ.3) THEN
         ETA=1.D0-COSCUT(2,I34-2)
         DMIN=2.D0*(AMASS2(I34  )*(DE/(EX+PX)+AMASS2(I34  )/2.D0/(EX+PX)**2)
     .             -AMASS2(I34-2)*(DP/(E1+P1)+AMASS2(I34-2)/2.D0/(E1+P1)**2))
     .       +2.D0*ETA*P1*PX
         ETA=1.D0-COSCUT(1,I34-2)
         DMAX=2.D0*(AMASS2(I34  )*(DE/(EX+PX)+AMASS2(I34  )/2.D0/(EX+PX)**2)
     .             -AMASS2(I34-2)*(DP/(E1+P1)+AMASS2(I34-2)/2.D0/(E1+P1)**2))
     .       +2.D0*ETA*P1*PX
         IF(DMIN.GT.DMAX) GOTO 999
         nev=nev+1
         if(dmin.le.0) goto 999
         D=DMIN*(DMAX/DMIN)**X(3)
         AJACOB=AJACOB*D*LOG(DMAX/DMIN)/P1/PX/2.D0
         PP13=(D+AMASS2(1)+AMASS2(I34))/2.D0
         XXXX=AMASS2(I34  )*(DE/(EX+PX)+AMASS2(I34  )/2.D0/(EX+PX)**2)
     .-       AMASS2(I34-2)*(DP/(E1+P1)+AMASS2(I34-2)/2.D0/(E1+P1)**2)-D/2.D0
         COSX=1.D0+XXXX/P1/PX
         SINX=SQRT(-XXXX/P1/PX*(2.D0+XXXX/P1/PX))
         TTT1=-D
        ELSE
         ETA=1.D0+COSCUT(1,I34-2)
         DMIN=2.D0*(AMASS2(I34  )*(DE/(EX+PX)+AMASS2(I34  )/2.D0/(EX+PX)**2)
     .             -AMASS2(I34-2)*(DP/(E2+P2)+AMASS2(I34-2)/2.D0/(E2+P2)**2))
     .       +2.D0*ETA*P2*PX
         ETA=1.D0+COSCUT(2,I34-2)
         DMAX=2.D0*(AMASS2(I34  )*(DE/(EX+PX)+AMASS2(I34  )/2.D0/(EX+PX)**2)
     .             -AMASS2(I34-2)*(DP/(E2+P2)+AMASS2(I34-2)/2.D0/(E2+P2)**2))
     .       +2.D0*ETA*P2*PX
         IF(DMIN.GT.DMAX) GOTO 999
         D=DMIN*(DMAX/DMIN)**X(3)
         AJACOB=AJACOB*D*LOG(DMAX/DMIN)/P2/PX/2.D0
         PP24=(D+AMASS2(2)+AMASS2(I34))/2.D0
         XXXX=AMASS2(I34  )*(DE/(EX+PX)+AMASS2(I34  )/2.D0/(EX+PX)**2)
     .       -AMASS2(I34-2)*(DP/(E2+P2)+AMASS2(I34-2)/2.D0/(E2+P2)**2)-D/2.D0
         COSX=1.D0+XXXX/P2/PX
         SINX=SQRT(-XXXX/P1/PX*(2.D0+XXXX/P1/PX))
         COSX=-COSX
         TTT2=-D
        END IF
       ELSE
        WRITE(6,*)' ICOS3 =',ICOS3,' is NOT supported. '
        stop
       END IF
       IF(ABS(COSX).GT.1) GOTO 999
***************
*    PHIX     *
*    X=3 or 4 *
***************
      if (Lpol_Ebeam .AND. ipol_Ebeam.GT.0) then
       PHIX = X(NDIM)*2D0*PI
      else
       if (ROT_kinem_flag) then
        PHIX = CepGen_random_number() *2D0*PI
       else
        PHIX = 0D0
       endif
      endif
      AJACOB=AJACOB*2.D0*PI
      PHI1 = PHIX
*Particle X
      PE(1,I34) = PX*SINX*COS(PHIX)
      PE(2,I34) = PX*SINX*SIN(PHIX)
      PE(3,I34) = PX*COSX
      PE(4,I34) = EX
*Particle 456
      PBOST1(1) =-PX*SINX*COS(PHIX)
      PBOST1(2) =-PX*SINX*SIN(PHIX)
      PBOST1(3) =-PX*COSX
      PBOST1(4) = EX56
***************
*    Q56      *
***************
      ipoint = 4
       Q56MN=MAX((           AMASS1(5)+AMASS1(6))**2,AMASCT(1,2)**2)
       Q56MX=MIN((SQRT(aX56)-AMASS1(7-I34)      )**2,AMASCT(2,2)**2)
       IF(Q56MX.LE.Q56MN) GOTO 999
       IF     (IRESNS(2).EQ. 0) THEN
        Q56=Q56MN+(Q56MX-Q56MN)*X(2)
        AJACOB=AJACOB*(Q56MX-Q56MN)
       ELSE IF(IRESNS(2).EQ. 40) THEN
        DD = Q56MX - Q56MN
        if (num_reso_mj .EQ. 0) then
         if (Diff_nn .EQ. 0) then
          Q56 = Q56MN + DD*X(2)**Inn
          AJACOB = AJACOB *DD*Dnn*X(2)**(Inn-1)
         else
          Q56 = Q56MN + DD*X(2)**Dnn
          AJACOB = AJACOB *DD*Dnn*X(2)**(Dnn-1D0)
         endif
        else
         RN_mj = CepGen_random_number()
         do j = 0, num_reso_mj
          if (RN_mj .LT. weisel_mj(j)) GOTO 4100
         enddo
         write(6,*) '!!!Warning in KINEM_4f!!!'
         write(6,*) '  ---> Inconsistent Weight_mj members!'
         write(6,*) '  ---> Weights_mj =',(weight_mj(i),i=0,num_reso_mj)
         write(6,*) '  ---> Good-bye!'
         STOP
 4100    continue
         if (j .EQ. 0) then
          if (Diff_nn .EQ. 0) then
           Q56 = Q56MN + DD*X(2)**Inn
          else
           Q56 = Q56MN + DD*X(2)**Dnn
          endif
         else
          AMG = mass_mj(j)*width_mj(j)
          AM2 = mass2_mj(j)
          FFmin = atan((Q56MN-AM2)/AMG)
          FFmax = atan((Q56MX-AM2)/AMG)
C          dFF = FFmax - FFmin
          Q56 = AM2 + AMG*tan(FFmin+(FFmax-FFmin)*X(2))
         endif
         Density = 0D0
         if (weight_mj(0) .GT. 0) then
          XX_inv = ((Q56-Q56MN)/DD)**(1D0/Dnn)
          if (Diff_nn .EQ. 0) then
           Y_inv = DD*Dnn*XX_inv**(Inn-1D0)
          else
           Y_inv = DD*Dnn*XX_inv**(Dnn-1D0)
          endif
          Density = Density + weight_mj(0)/Y_inv
         endif
         do j = 1, num_reso_mj
          AMG = mass_mj(j)*width_mj(j)
          AM2 = mass2_mj(j)
          FFmin = atan((Q56MN-AM2)/AMG)
          FFmax = atan((Q56MX-AM2)/AMG)
C          dFF = FFmax - FFmin
          Y_inv = ((Q56-AM2)**2 + AMG*AMG) *(FFmax-FFmin)/AMG
          Density = Density + weight_mj(j)/Y_inv
         enddo
         AJACOB = AJACOB *(1D0/Density)
        endif
       ELSE IF(IRESNS(2).EQ. 1) THEN
        ZM    =ARESNS(1,2)
        ZM2   =ZM*ZM
        ZMG   =ARESNS(1,2)*ARESNS(2,2)
        THEMIN=ATAN((Q56MN-ZM2)/ZMG)
        THEMAX=ATAN((Q56MX-ZM2)/ZMG)
        THE   =THEMIN+(THEMAX-THEMIN)*X(2)
        Q56   =ZMG*TAN(THE)+ZM2
        AJACOB=AJACOB*(THEMAX-THEMIN)*((Q56-ZM2)**2+ZMG**2)/ZMG
       ELSE IF(IRESNS(2).EQ.-1) THEN
        IF(Q56MN.LE.0) goto 999
        Q56=Q56MN*(Q56MX/Q56MN)**X(2)
        AJACOB=AJACOB*Q56*LOG(Q56MX/Q56MN)
       ELSE IF(IRESNS(2).EQ.-2) THEN
        IF(Q56MN.LE.0) GOTO 999
        t_min = sqrt(Q56MN)
        t_max = sqrt(Q56MX)
        t_diff = t_max - t_min
        t_val = t_min + t_diff * X(2)
        Q56 = t_val**2
        AJACOB=AJACOB* 2.D0*t_val*t_diff
       ELSE IF(IRESNS(2).EQ. 2) THEN
        IF(X(2).LT.0.5D0) THEN
         X2=X(2)*2.D0
         q56MX=(Q56MX+Q56MN)/2.D0
         IF(Q56MN.LE.0) GOTO 999
         Q56=Q56MN*(Q56MX/Q56MN)**X2
         AJACOB=AJACOB*Q56*LOG(Q56MX/Q56MN)*2.d0
        ELSE
         X2=2.D0-X(2)*2.D0
         Q56MN=(Q56MX+Q56MN)/2.D0
         Q56=Q56MN+(Q56MX-Q56MN)*X2
         AJACOB=AJACOB*(Q56MX-Q56MN)*2.d0
        END IF
       ELSE IF(IRESNS(2).EQ. 3) THEN
        thre = thresh56
        IF(X(2).LT.thre) THEN
         X2=X(2)/thre
         Q56MX=MIN(Q56MX,(ARESNS(1,2)-3*ARESNS(2,2))**2)
         IF(Q56MX.LT.Q56MN) GOTO 999
         IF(Q56MN.LE.0)     GOTO 999
         Q56=Q56MN*(Q56MX/Q56MN)**X2
         AJACOB=AJACOB*Q56*LOG(Q56MX/Q56MN)/thre
        ELSE
         X2=(x(2)-thre)/(1.d0-thre)
         Q56MN=MAX(Q56MN,(ARESNS(1,2)-3*ARESNS(2,2))**2)
         IF(Q56MX.LT.Q56MN) GOTO 999
         ZM    =ARESNS(1,2)
         ZM2   =ZM*ZM
         ZMG   =ARESNS(1,2)*ARESNS(2,2)
         THEMIN=ATAN((Q56MN-ZM2)/ZMG)
         THEMAX=ATAN((Q56MX-ZM2)/ZMG)
         THE   =THEMIN+(THEMAX-THEMIN)*X2
         Q56   =ZMG*TAN(THE)+ZM2
         AJACOB=AJACOB*(THEMAX-THEMIN)*((Q56-ZM2)**2+ZMG**2)/ZMG/(1.d0-thre)
        END IF
       ELSE IF(IRESNS(2).EQ. 4) THEN
        thre = thresh56
        IF(X(2).LT.thre) THEN
          X2=X(2)/thre
          Q56MX=MIN(Q56MX,(ARESNS(1,2)-3*ARESNS(2,2))**2)
          IF(Q56MX.LT.Q56MN) GOTO 999
          IF(Q56MN.LE.0)     GOTO 999
          t_min = sqrt(Q56MN)
          t_max = sqrt(Q56MX)
          t_diff = t_max - t_min
          t_val = t_min + t_diff * X2
          Q56 = t_val**2
          AJACOB=AJACOB* 2.D0*t_val*t_diff/thre
        ELSE
          X2=(x(2)-thre)/(1.d0-thre)
          Q56MN=MAX(Q56MN,(ARESNS(1,2)-3*ARESNS(2,2))**2)
          IF(Q56MX.LT.Q56MN) GOTO 999
          ZM    =ARESNS(1,2)
          ZM2   =ZM*ZM
          ZMG   =ARESNS(1,2)*ARESNS(2,2)
          THEMIN=ATAN((Q56MN-ZM2)/ZMG)
          THEMAX=ATAN((Q56MX-ZM2)/ZMG)
          THE   =THEMIN+(THEMAX-THEMIN)*X2
          Q56   =ZMG*TAN(THE)+ZM2
          AJACOB=AJACOB*(THEMAX-THEMIN)*((Q56-ZM2)**2+ZMG**2)/ZMG/(1.d0-thre)
        END IF
       ELSE IF(IRESNS(2).EQ. 5) THEN
        thre = thresh56
        IF(X(2).LT.thre) THEN
         X2=X(2)/thre
         Q56MX=MIN(Q56MX,(ARESNS(1,2)-3*ARESNS(2,2))**2)
         IF(Q56MX.LT.Q56MN) GOTO 999
         IF(Q56MN.LE.0)     GOTO 999
         Q56=Q56MN+(Q56MX-Q56MN)*X2
         AJACOB=AJACOB*(Q56MX-Q56MN)/thre
        ELSE
         X2=(x(2)-thre)/(1.d0-thre)
         Q56MN=MAX(Q56MN,(ARESNS(1,2)-3*ARESNS(2,2))**2)
         IF(Q56MX.LT.Q56MN) GOTO 999
         ZM    =ARESNS(1,2)
         ZM2   =ZM*ZM
         ZMG   =ARESNS(1,2)*ARESNS(2,2)
         THEMIN=ATAN((Q56MN-ZM2)/ZMG)
         THEMAX=ATAN((Q56MX-ZM2)/ZMG)
         THE   =THEMIN+(THEMAX-THEMIN)*X2
         Q56   =ZMG*TAN(THE)+ZM2
         AJACOB=AJACOB*(THEMAX-THEMIN)*((Q56-ZM2)**2+ZMG**2)/ZMG/(1.d0-thre)
        END IF
       ELSE
        WRITE(6,*)' IRESNS(2) =',IRESNS(2),' is NOT supported. '
        stop
       END IF
       EX  = (QX56+AMASS2(7-I34)-Q56)/2.D0/SQRT(QX56)
       PX  = SQRT( max( EX-AMASS1(7-I34), 0.D0 ) * (EX+AMASS1(7-I34)) )
       E56 = (QX56+Q56-AMASS2(7-I34))/2.D0/SQRT(QX56)
       if (EX-AMASS1(7-I34) .LE. 0) then
C        write(6,*) '!!!Warning in K6NEM0!!!'
C        write(6,*) '  ---> EX            =', EX
C        write(6,*) '  ---> I34           =', I34
C        write(6,*) '  ---> AMASS1(7-I34) =', AMASS1(7-I34)
C        write(6,*) '  ---> Calculated PX =', PX
C        write(6,*) '  ---> This event is rejected.'
C        write(6,*) '  (User need not to worry about this.)'
        GOTO 999
       endif
***************
*    COSX     *
***************
      ipoint = 5
      X7=X(7)
      CALL LABTOK(PE(1,5-I34),PBOST1,pe(1,5-I34),PDUM,PK1,PBOST2)
      EE1=PK1(4)
      PP1=SQRT( (EE1-AMASS1(5-I34))*(EE1+AMASS1(5-I34)) )
      DE=EE1-EX
      DP=(DE*(EX+EE1)-AMASS2(5-I34)+AMASS2(7-I34))/(PP1+PX)
      TTMIN=2.D0*(AMASS2(7-I34)*(DE/(EX +PX )+AMASS2(7-I34)/2.D0/(EX +PX )**2)
     .           -AMASS2(5-I34)*(DP/(EE1+PP1)+AMASS2(5-I34)/2.D0/(EE1+PP1)**2))
      TTMAX=2.D0*(AMASS2(7-I34)*(DE/(EX +PX )+AMASS2(7-I34)/2.D0/(EX +PX )**2)
     .           -AMASS2(5-I34)*(DP/(EE1+PP1)+AMASS2(5-I34)/2.D0/(EE1+PP1)**2))
     .     +2.D0*2.D0*PX*PP1
      TMAX=TTMAX
      TMIN=TTMIN
      IF(TMAX.LT.TMIN) GOTO 999
      DTT=TMAX-TMIN
      IF(abs(coscut(i34-2,5-i34)).lt.cos(4.d0*rad)) THEN
       TT = TMIN+DTT*X7
      ELSE
       EPS=10D0**(-Neps_p)
       YMIN=-LOG( 1.D0+2.D0/EPS )/2
       YMAX= LOG( 1.D0+2.D0/EPS )/2
       Y=YMIN+(YMAX-YMIN)*X7
       Y7 = ((2.0D0+EPS)*EXP(2.0D0*Y) - EPS)/(EXP(2.0D0*Y) + 1.0D0 )/2.0D0
       AJACOB=AJACOB*(YMAX-YMIN)/COSH(Y)**2*(1+EPS)/2.D0
       TT = TMIN+DTT*Y7
      END IF
      AJACOB=AJACOB*DTT/PX/PP1/2.D0
      PPPP=(TT+AMASS2(5-I34)+AMASS2(7-I34))/2.D0
      XXXX=AMASS2(7-I34)*(DE/(EX +PX )+AMASS2(7-I34)/2.D0/(EX +PX )**2)
     .    -AMASS2(5-I34)*(DP/(EE1+PP1)+AMASS2(5-I34)/2.D0/(EE1+PP1)**2)-TT/2.D0
      IF(-XXXX/PP1/PX*(2.D0+XXXX/PP1/PX).GT.0.D0) THEN
       COSX=1.D0+XXXX/PX/PP1
       SINX=SQRT(-XXXX/PP1/PX*(2.D0+XXXX/PP1/PX))
      ELSE
       COSX=1.D0
       SINX=0.D0
      END IF
      IF(I34.EQ.3) THEN
       PP24=PPPP
       TTT2=-TT
      ELSE
       PP13=PPPP
       TTT1=-TT
      END IF
***************
*    PHIX     *
***************
      PHIX = X4*2*PI+PHI1
      AJACOB=AJACOB*2.D0*PI
*Particle X in px-p5-p6 CM frame
      PE(1,7-I34) = PX*SINX*COS(PHIX)
      PE(2,7-I34) = PX*SINX*SIN(PHIX)
      PE(3,7-I34) = PX*COSX
      PE(4,7-I34) = EX
*Particle 56 in px-p5-p6 CM frame
      PBOST2(1) =-PX*SINX*COS(PHIX)
      PBOST2(2) =-PX*SINX*SIN(PHIX)
      PBOST2(3) =-PX*COSX
      PBOST2(4) = E56
* boost from x-5-6 CM frame to 1-2 CM frame (Lab-frame)
      CALL KTOLAB(PE(1,5-I34),PBOST1,pe(1,7-I34),PBOST2,
     .            PE(1,7-I34),PBOST2)
      PP34=PE(4,4)*PE(4,3)-PE(3,4)*PE(3,3)-PE(1,4)*PE(1,3)-PE(2,4)*PE(2,3)
      Q34=AMASS2(3)+AMASS2(4)+2.D0*PP34
**********************************************************************
      IF(ISYM.EQ.1) THEN
       if(I34.eq.3 .and. pbost2(3).gt.0) goto 999
       if(I34.eq.4 .and. pbost2(3).lt.0) goto 999
      END IF
**********************************************************************
      IF(ICOS5.EQ.2) THEN
       TT1(1)=PE(1,I34-2)-PE(1,I34)
       TT1(2)=PE(2,I34-2)-PE(2,I34)
       TT1(3)=PE(3,I34-2)-PE(3,I34)
       TT1(4)=PE(4,I34-2)-PE(4,I34)
       TT2(1)=PE(1,I34-2)-PE(1,I34)
       TT2(2)=PE(2,I34-2)-PE(2,I34)
       TT2(3)=PE(3,I34-2)-PE(3,I34)
       TT2(4)=PE(4,I34-2)-PE(4,I34)
       PK1(1)=PE(1,5-I34)-PE(1,7-I34)
       PK1(2)=PE(2,5-I34)-PE(2,7-I34)
       PK1(3)=PE(3,5-I34)-PE(3,7-I34)
       PK1(4)=PE(4,5-I34)-PE(4,7-I34)
       CALL LABTOK(TT2,PBOST2,TT1,PK1, TT1,PK1)
      ELSE IF(ICOS5.EQ.1) THEN
       TT1(1)=PE(1,5-I34)
       TT1(2)=PE(2,5-I34)
       TT1(3)=PE(3,5-I34)
       TT1(4)=PE(4,5-I34)
       TT2(1)=PE(1,5-I34)
       TT2(2)=PE(2,5-I34)
       TT2(3)=PE(3,5-I34)
       TT2(4)=PE(4,5-I34)
       CALL LABTOK(TT2,PBOST2,TT1,PDUM, TT1,PK1)
      END IF
***************
*    COS5     *
***************
      ipoint = 6
      E5 = (Q56+AMASS2(5)-AMASS2(6))/2.D0/SQRT(Q56)
      P5 = SQRT( max( E5-AMASS1(5), 0.D0 ) * (E5+AMASS1(5)) )
      E6 = (Q56+AMASS2(6)-AMASS2(5))/2.D0/SQRT(Q56)
      P6 = P5
      ITU= 0
      IF(ICOS5.EQ.0) THEN
       COS5 = -1+2*X(5)
       AJACOB=AJACOB*2
       sin5=sqrt( max( 1.d0-cos5, 0.D0 ) * (1.d0+cos5) )
      ELSE IF(ICOS5.EQ.1) THEN
       EPS=1.D-2
       YMIN=-LOG( 1.D0+2.D0/EPS )/2
       YMAX= LOG( 1.D0+2.D0/EPS )/2
       Y=YMIN+(YMAX-YMIN)*X(5)
       COS5=(1+EPS)*TANH(Y)
       AJACOB=AJACOB*(YMAX-YMIN)/COSH(Y)**2*(1+EPS)
       sin5=sqrt( max( 1.d0-cos5, 0.D0 ) * (1.d0+cos5) )
      ELSE IF(ICOS5.EQ.2) THEN
       IF(X(5).LT.0.5D0) THEN
        X5=X(5)*2.D0
        ET1=TT1(4)
        PT1=SQRT( ET1**2+D)
        T2MIN=D-AMASS2(5)+2.D0*(ET1*AMASS2(5)/(E5+P5)-P5*D/(ET1+PT1))
        T2MAX=D-AMASS2(5)+2.D0*(ET1*E5       )
        IF(T2MAX.Le.T2MIN) GOTO 999
        IF(T2MIN.LE.0) GOTO 999
        T2=T2MIN*(T2MAX/T2MIN)**X5
        AJACOB=AJACOB*T2*LOG(T2MAX/T2MIN)/PT1/P5
        XXXX = T2+AMASS2(5)-D-2.D0*(ET1*AMASS2(5)/(E5+P5)-P5*D/(ET1+PT1))
        COS5=1.D0-XXXX/P5/PT1/2.D0
        SIN5=SQRT(XXXX/P5/PT1/2.D0*(2.D0-XXXX/P5/PT1/2.D0))
        U2=-(AMASS2(5)+AMASS2(6)-D-TT-Q56+T2)
        ITU=1
       ELSE
        X5=2.D0-X(5)*2.D0
        ET1=TT1(4)
        PT1=SQRT( ET1**2+D)
        U2MIN=D-AMASS2(6)+2.D0*(ET1*AMASS2(6)/(E6+P6)-P6*D/(ET1+PT1))
        U2MAX=D-AMASS2(6)+2.D0*(ET1*E6       )
        IF(U2MAX.LT.U2MIN) GOTO 999
        IF(U2MIN.LE.0) GOTO 999
        U2=U2MIN*(U2MAX/U2MIN)**X5
        AJACOB=AJACOB*U2*LOG(U2MAX/U2MIN)/PT1/P6
        XXXX = U2+AMASS2(6)-D-2.D0*(ET1*AMASS2(6)/(E6+P6)-P6*D/(ET1+PT1))
        COS6=1.D0-XXXX/P6/PT1/2.D0
        SIN5=SQRT(XXXX/P6/PT1/2.D0*(2.D0-XXXX/P6/PT1/2.D0))
        COS5=-COS6
        T2=-(AMASS2(5)+AMASS2(6)-D-TT-Q56+U2)
        ITU=2
       END IF
      ELSE
       WRITE(6,*)' ICOS5 =',ICOS5,' is NOT supported. '
       stop
      END IF
      IF(ABS(COS5).GT.1) GOTO 999
***************
*   PHI_P5    *
***************
      PHI5 = X(6)*2*PI+PHI1
      AJACOB=AJACOB*2*PI
*Particle 5 in 5-6 CM frame
      PE(1,5) = P5*SIN5*COS(PHI5)
      PE(2,5) = P5*SIN5*SIN(PHI5)
      PE(3,5) = P5*COS5
      PE(4,5) = E5
*Particle 6 in 5-6 CM frame
      PE(1,6) =-PE(1,5)
      PE(2,6) =-PE(2,5)
      PE(3,6) =-PE(3,5)
      PE(4,6) = E6
* boost from 5-6 CM frame to 1-2 CM frame (Lab-frame)
      IF(ICOS5.le.1) THEN
       CALL WTOLAB(PE(1,5),PE(1,6),PBOST2, PE(1,5),PE(1,6))
      ELSE
       CALL KTOLAB(TT2,PBOST2,PE(1,5),PE(1,6), PE(1,5),PE(1,6))
      END IF
      ipoint = 7
      ipoint = 71
13    continue
12    continue
*Set invariants
      DO I1=1,6
       DO J1=1,6
        PP(I1,J1)=PE(4,I1)*PE(4,J1)-PE(3,I1)*PE(3,J1)
     .-           PE(1,I1)*PE(1,J1)-PE(2,I1)*PE(2,J1)
       ENDDO
      ENDDO
      PP(1,2) = (S(3) -AMASS2(1)-AMASS2(2))/2.D0
      PP(2,1) = (S(3) -AMASS2(1)-AMASS2(2))/2.D0
      PP(1,3)=PP13
      PP(3,1)=PP13
      PP(2,4)=PP24
      PP(4,2)=PP24
      PP(3,4)=(Q34-AMASS2(3)-AMASS2(4))/2.D0
      PP(4,3)=(Q34-AMASS2(3)-AMASS2(4))/2.D0
      PP(5,6)=(Q56-AMASS2(5)-AMASS2(6))/2.D0
      PP(6,5)=(Q56-AMASS2(5)-AMASS2(6))/2.D0
      PP(7-I34,5)=(-Q56-AMASS2(7-I34)+QX56-2*PP(7-I34,6))/2.D0
      PP(5,7-I34)=(-Q56-AMASS2(7-I34)+QX56-2*PP(7-I34,6))/2.D0
      IF(I34.EQ.3) THEN
       SSS1=AMASS2(3)+Q56+2.D0*(PP(3,5)+PP(3,6))
       SSS2=QX56
      ELSE
       SSS1=QX56
       SSS2=AMASS2(4)+Q56+2.D0*(PP(4,5)+PP(4,6))
      END IF
      UUU1=AMASS2(1)+AMASS2(3)+TTT2+Q56-SSS1-TTT1
      UUU2=AMASS2(2)+AMASS2(4)+TTT1+Q56-SSS2-TTT2
      IF(ICOS5.EQ.2) THEN
       IF(I34.EQ.3) THEN
        TTTT=-T2
        UUUU=-U2
       ELSE
        TTTT=-U2
        UUUU=-T2
       END IF
      END IF
      colmbf=1.d0
      FACT = Flux_factor( S(3), amass2(1), amass2(2), jump )
      if (jump.EQ.1) then
       YACOB = 0
       RETURN
      endif
*Set jacobian
      YACOB = FACT*AJACOB/(2*PI)**2/(32*PI2)**3
     .      * aBETA(AMASS2(I34  )/S(3), aX56     /S(3))
     .      * aBETA(AMASS2(7-I34)/aX56, Q56      /aX56)
     .      * aBETA(AMASS2(5    )/Q56,  AMASS2(6)/Q56 )
      RETURN
999   CONTINUE
      JUMP=1
      yacob=0
      coscut(1,1)=csut11
      coscut(2,1)=csut21
      coscut(1,2)=csut12
      coscut(2,2)=csut22
      coscut(1,3)=csut13
      coscut(2,3)=csut23
      coscut(1,4)=csut14
      coscut(2,4)=csut24
      RETURN
      END
*-----------------------------------------------------------------------
      function abeta(z1,z2)
      implicit real*8(a-h,o-z)
      abeta=sqrt( max(1d0-2d0*(z1+z2)+(z1-z2)**2,0d0) )
      return
      end
