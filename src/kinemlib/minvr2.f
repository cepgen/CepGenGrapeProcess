       SUBROUTINE MINVR2(M1,M2)
       REAL*8 M1(3,3),M2(3,3),M3(3,3)
       DO I=1,3
        DO J=1,3
         M3(I,J)=M1(J,I)
        ENDDO
       ENDDO
       DO I=1,3
        DO J=1,3
         M2(I,J)=M3(I,J)
        ENDDO
       ENDDO
       RETURN
       END
