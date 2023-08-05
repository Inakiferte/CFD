**==rantest.spg  processed by SPAG 4.52O  at 15:46 on 28 Mar 1996
 
      SUBROUTINE RANTEST(Iseed)
c
c     test and initialize the random number generator
c
      IMPLICIT NONE
      INTEGER Iseed, i
      DOUBLE PRECISION RANF
 
      CALL RANSET(Iseed)
C----------------------------------------------------------
      WRITE(99,*)' ******** Begin test random numbers *********'
      DO i = 1, 5
CC         PRINT *, ' i,ranf() ', i, RANF(Iseed)
         WRITE(99,'(a4,i4,a10,f20.10)')' i=',i,' ranf()=',RANF(Iseed)
      END DO
      WRITE(99,*)' ******** End test random numbers ***********'
C----------------------------------------------------------
      RETURN
      END
