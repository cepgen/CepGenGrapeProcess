      subroutine  Add_Graph(flag, num_flag, Igraph)
      implicit NONE
* --------------- Argument ---------------
      integer   num_flag              ! Input
      integer   flag(num_flag)        ! Input/Output
      integer   Igraph                ! Input
* ----------------------------------------
* ---------- Local variables ----------
      integer   Iflag, Igraph_local
* -------------------------------------
      Iflag = int( (Igraph-1)/30 ) + 1
      if (Iflag .GT. num_flag) then
        write(6,*) '!!!Error in Add_Graph!!!'
        write(6,*) ' ---> Iflag    =', Iflag
        write(6,*) ' ---> num_flag =', num_flag
        write(6,*) ' ---> Iflag should be <= num_flag.'
        write(6,*) ' ---> Good-bye!'
        STOP
      endif
      Igraph_local = mod(Igraph-1,30) + 1      ! from 1 to 30
      if (btest(flag(Iflag),Igraph_local-1)) then
        write(6,*) ' '
        write(6,*) '!!!Error in Add_Graph!!!'
        write(6,*) ' ---> Graph(',Igraph,') is already fired.'
        write(6,*) ' ---> This call is ignored.'
        write(6,*) ' '
      endif
      flag(Iflag) = flag(Iflag) + 2**(Igraph_local-1)
      return
      end
*=====================================================================
      subroutine  Sub_Graph(flag, num_flag, Igraph)
      implicit NONE
* --------------- Argument ---------------
      integer   num_flag              ! Input
      integer   flag(num_flag)        ! Input/Output
      integer   Igraph                ! Input
* ----------------------------------------
* ---------- Local variables ----------
      integer   Iflag, Igraph_local
* -------------------------------------
      Iflag = int( (Igraph-1)/30 ) + 1
      if (Iflag .GT. num_flag) then
        write(6,*) '!!!Error in Sub_Graph!!!'
        write(6,*) ' ---> Iflag    =', Iflag
        write(6,*) ' ---> num_flag =', num_flag
        write(6,*) ' ---> Iflag should be <= num_flag.'
        write(6,*) ' ---> Good-bye!'
        STOP
      endif
      Igraph_local = mod(Igraph-1,30) + 1      ! from 1 to 30
      if (.NOT.btest(flag(Iflag),Igraph_local-1)) then
        write(6,*) ' '
        write(6,*) '!!!Error in Sub_Graph!!!'
        write(6,*) ' ---> Graph(',Igraph,') is not yet fired.'
        write(6,*) ' ---> This call is ignored.'
        write(6,*) ' '
      endif
      flag(Iflag) = flag(Iflag) - 2**(Igraph_local-1)
      return
      end
