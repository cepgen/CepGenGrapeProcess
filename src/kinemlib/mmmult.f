       SUBROUTINE MMMULT(R1,R2,R3)
       REAL*8 R1(3,3),R2(3,3),R3(3,3),RDUMMY(3,3)
       CALL MEMCLR(RDUMMY,9,0)
       DO I=1,3
        DO J=1,3
         DO K=1,3
          RDUMMY(I,J)=RDUMMY(I,J)+R1(I,K)*R2(K,J)
         ENDDO
        ENDDO
       ENDDO
       DO I=1,3
        DO J=1,3
         R3(I,J)=RDUMMY(I,J)
        ENDDO
       ENDDO
C
       RETURN
       END
