*-----------------------------------------------------------------------
*     SMFFV revised 94/04/08 by T.Ishikawa
*-----------------------------------------------------------------------
      SUBROUTINE SMFFV(L2,L1,LV,EW2,EW1,AM2,AM1,CPL,CE2,CE1,
     &                 PS2,PS1,EP,LT,AV)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
*   * dummy array size.
      PARAMETER (LTSIZE=20,LASIZE = 1024)
      INTEGER    L1, L2, LV
      DIMENSION  EW1(1), EW2(1)

      COMPLEX*16 CPL(2)
      COMPLEX*16 CE1(2,L1), CE2(2,L2)
      DIMENSION  PS1(4,3), PS2(4,3)

      DIMENSION  EP(4,LV)
      COMPLEX*16 AV(0:LASIZE)

      INTEGER    LT(0:LTSIZE)

      COMPLEX*16 AVT(4,2,2,2,2), CPN
      DIMENSION     CPR(2)

*--------------------------- Entry point -------------------------------

      IF( CPL(1) .NE. ZERO) THEN
          CPN      = CPL(1)
          CPR(1) = ONE
          CPR(2) = CPL(2)/CPL(1)
      ELSEIF( CPL(2).NE.ZERO) THEN
          CPN      = CPL(2)
          CPR(2) = ONE
          CPR(1) = CPL(1)/CPL(2)
      ELSE
          CPN      = ZERO
          CPR(1) = ONE
          CPR(2) = ONE
      ENDIF

      IF(EW1(1).GE.ZERO) THEN
         K1 = 3
      ELSE
         K1 = 1
      ENDIF
      IF(EW2(1).GE.ZERO) THEN
         K2 = 3
      ELSE
         K2 = 1
      ENDIF

      LT(0) =  3
      LT(1) =  L2
      LT(2) =  L1
      LT(3) =  LV

      IF( AM1.GT.0.0D0 )THEN
          IF( AM2.GT.0.0D0 )THEN

              DO LP = 1, LV
                 CALL FFVMM1(LP,K1,K2,CPR(1),CPR(2),
     .                      CE1(1,2),CE2(1,1),PS1(1,1),PS1(1,2),
     .                      PS2(1,1),PS2(1,2),EP(1,LP),AVT(1,1,1,1,1))
              ENDDO

              IF(L1.EQ.4) THEN
                 DO LP = 1, LV
                    CALL FFVMM2(LP,CPR(1),CPR(2),CE2(1,1),PS1(1,3),
     &                        EP(1,LP),AVT(1,1,1,2,1))
                 ENDDO
              ENDIF

              IF(L2.EQ.4) THEN
                 DO LP = 1, LV
                    CALL FFVMM3(LP,CPR(1),CPR(2),CE1(1,2),PS2(1,3),
     &                          EP(1,LP),AVT(1,1,1,1,2))
                 ENDDO
              ENDIF

              IF(L1.EQ.4 .AND. L2.EQ.4) THEN
                 DO LP = 1, LV
                    CALL FFVMM4(LP,CPR(1),CPR(2),
     .                          EP(1,LP),AVT(1,1,1,2,2))
                 ENDDO
              ENDIF

          ELSE

              DO LP = 1, LV
                 CALL FFVMM5(LP,K1,K2,CPR(1),CPR(2),
     .                      CE1(1,2),CE2(1,1),PS1(1,1),PS1(1,2),
     .                      PS2(1,1),PS2(1,2),EP(1,LP),AVT(1,1,1,1,1))
              ENDDO

              IF(L1.EQ.4) THEN
                 DO LP = 1, LV
                    CALL FFVMM6(LP,CPR(1),CPR(2),CE2(1,1),PS1(1,3),
     &                        EP(1,LP),AVT(1,1,1,2,1))
                 ENDDO
              ENDIF

              IF(L2.EQ.4) THEN
                 DO LP = 1, LV
                    CALL FFVMM3(LP,CPR(1),CPR(2),CE1(1,2),PS2(1,3),
     &                          EP(1,LP),AVT(1,1,1,1,2))
                 ENDDO
              ENDIF

              IF(L1.EQ.4 .AND. L2.EQ.4) THEN
                 DO LP = 1, LV
                    CALL FFVMM4(LP,CPR(1),CPR(2),
     .                          EP(1,LP),AVT(1,1,1,2,2))
                 ENDDO
              ENDIF
          ENDIF
      ELSE
          IF( AM2.GT.0.0D0 )THEN
              DO LP = 1, LV
                 CALL FFVMM7(LP,K1,K2,CPR(1),CPR(2),
     .                      CE1(1,2),CE2(1,1),PS1(1,1),PS1(1,2),
     .                      PS2(1,1),PS2(1,2),EP(1,LP),AVT(1,1,1,1,1))
              ENDDO

              IF(L1.EQ.4) THEN
                 DO LP = 1, LV
                    CALL FFVMM2(LP,CPR(1),CPR(2),CE2(1,1),PS1(1,3),
     &                        EP(1,LP),AVT(1,1,1,2,1))
                 ENDDO
              ENDIF

              IF(L2.EQ.4) THEN
                 DO LP = 1, LV
                    CALL FFVMM8(LP,CPR(1),CPR(2),CE1(1,2),PS2(1,3),
     &                          EP(1,LP),AVT(1,1,1,1,2))
                 ENDDO
              ENDIF

              IF(L1.EQ.4 .AND. L2.EQ.4) THEN
                 DO LP = 1, LV
                    CALL FFVMM4(LP,CPR(1),CPR(2),
     .                          EP(1,LP),AVT(1,1,1,2,2))
                 ENDDO
              ENDIF

          ELSE

              DO LP = 1, LV
                 CALL FFVMM9(LP,K1,K2,CPR(1),CPR(2),
     .                      CE1(1,2),CE2(1,1),PS1(1,1),PS1(1,2),
     .                      PS2(1,1),PS2(1,2),EP(1,LP),AVT(1,1,1,1,1))
              ENDDO

              IF(L1.EQ.4) THEN
                 DO LP = 1, LV
                    CALL FFVMM6(LP,CPR(1),CPR(2),CE2(1,1),PS1(1,3),
     &                        EP(1,LP),AVT(1,1,1,2,1))
                 ENDDO
              ENDIF

              IF(L2.EQ.4) THEN
                 DO LP = 1, LV
                    CALL FFVMM8(LP,CPR(1),CPR(2),CE1(1,2),PS2(1,3),
     &                          EP(1,LP),AVT(1,1,1,1,2))
                 ENDDO
              ENDIF

              IF(L1.EQ.4 .AND. L2.EQ.4) THEN
                 DO LP = 1, LV
                    CALL FFVMM4(LP,CPR(1),CPR(2),
     .                          EP(1,LP),AVT(1,1,1,2,2))
                 ENDDO
              ENDIF

          ENDIF
      ENDIF

      IA = 0
      DO IL  = 1, LV
       DO IP1 = 1, L1/2
        DO IL1 = 1, 2
         DO IP2 = 1, L2/2
          DO IL2 = 1, 2
           AV(IA) = CPN*AVT(IL, IL1, IL2, IP1, IP2)
           IA = IA + 1
*           print *,'smffv i,l,p1,l1,p2,l2',ia,il,ip1,il1,ip2,il2,
*           AVT(IL, IL1, IL2, IP1, IP2)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END
