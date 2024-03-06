      integer function  Ibtest(Ivar,bit)
      implicit NONE
*----------- Argument -----------
      integer   Ivar,bit
*--------------------------------
*-------- Local variables --------
      integer*4  Ivar4,bit4
*---------------------------------
      Ivar4 = Ivar
      bit4  = bit
      if (btest(Ivar4,bit4)) then
         Ibtest = 1
       else
         Ibtest = 0
      endif
      return
      end
