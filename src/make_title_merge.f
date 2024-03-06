      subroutine  make_title_merge(
     &                  process, lpair, merge, Ce           ! Inputs
     &                 ,Ctitle                              ! Output
     &                            )
      implicit NONE
* --------------- Argument ---------------
      integer        process, lpair, merge        ! Input
      character*(*)  Ce                           ! Input
      character*(*)  Ctitle                       ! Output
* ----------------------------------------
* -------- Local variables --------
      character  Clpair(3)*5
* ---------------------------------
      Clpair(1) = 'e+ e-'
      Clpair(2) = 'm+ m-'
      Clpair(3) = 't+ t-'
      if (process .EQ. 3) then
        if     (merge .EQ. 12)       then        ! u+u^
          Ctitle = Ce // 'uu^'         // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 34)       then        ! d+d^
          Ctitle = Ce // 'dd^'         // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 56)       then        ! s+s^
          Ctitle = Ce // 'ss^'         // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 78)       then        ! c+c^
          Ctitle = Ce // 'cc^'         // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 90)       then        ! b+b^
          Ctitle = Ce // 'bb^'         // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 1234)     then        ! u+u^+d+d^
          Ctitle = Ce // 'uu^dd^'      // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 123456)   then        ! u+u^+d+d^+s+s^
          Ctitle = Ce // 'uu^dd^ss^'   // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 12345678)   then      ! u+u^+d+d^+s+s^+c+c^
          Ctitle = Ce // 'uu^-cc^'     // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 1234567890) then      ! u+u^+d+d^+s+s^+c+c^+b+b^
          Ctitle = Ce // 'uu^-bb^'     // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 17) then              ! u+c
          Ctitle = Ce // 'uc'          // ' -> '
     &          // Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 28) then              ! u^+c^
          Ctitle = Ce // 'u^c^'        // ' -> '
     &          // Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 35) then              ! d+s
          Ctitle = Ce // 'ds'          // ' -> '
     &          // Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 46) then              ! d^+s^
          Ctitle = Ce // 'd^s^'        // ' -> '
     &          // Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 359) then             ! d+s+b
          Ctitle = Ce // 'dsb'         // ' -> '
     &          // Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 460) then             ! d^+s^+b^
          Ctitle = Ce // 'd^s^b^'      // ' -> '
     &          // Ce // 'q'           // ' ' // Clpair(lpair)
        else
           write(6,*) '!!!Error in Make_Title_merge!!!'
           write(6,*) ' ---> Not-supported merge(merge=',merge,')'
           write(6,*) ' ---> Good-bye!'
           STOP
        endif
      else
        write(6,*) '!!!Error in Make_Title_merge!!!'
        write(6,*) ' ---> Not-supported process(process=',process,')'
        write(6,*) ' ---> Good-bye!'
        STOP
      endif
      return
      end
