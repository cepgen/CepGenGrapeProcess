      subroutine  START_GRAPE(istep)
      implicit NONE
* --------- Argument ---------
      integer  istep
* ----------------------------
* ------ Local variables ------
      character  c(30)*64, cside(8)*8,cleft*8,cright*8, cstep*64
      integer    i,j,iclast,icend(2)
* -----------------------------
      c( 1)='########################################################'
      c( 2)='********************************************************'
      c( 3)='++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      c( 4)='   GGGGGG     RRRRRR        A      PPPPPPP     EEEEEE   '
      c( 5)='  G      G   R      R      A A     P      P   E         '
      c( 6)='  G          R       R    A   A    P       P  E         '
      c( 7)='  G  GGGGGG  R  RRRR     A     A   PPPPPPPP   EEEEEEE   '
      c( 8)='  G     G G  R   R       AAAAAAA   P          E         '
      c( 9)='  G     G G  R    R     A       A  P          E         '
      c(10)='   GGGGG  G  R      RR  A       A  P          EEEEEEEE  '
      c(11)='                                                        '
      c(12)='  GRAce-based generator for Proton-Electron collisions  '
      c(13)='                                                        '
      c(14)='               GRAPE-Dilepton_version1.1k               '
      c(15)='                     ^^^^^^^^                           '
      c(16)='                      Mar.27, 2003                      '
      c(17)='  Comments/bug-report to Tetsuo ABE (tabe@post.kek.jp)  '
      c(18)='++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      c(19)='********************************************************'
      c(20)='########################################################'
      iclast = 20
      cside(1) = '######'
       cside(2) = '##**++'
       cside(3) = '++**##'
        cside(4) = '##****'
        cside(5) = '****##'
         cside(6) = '##**++'
         cside(7) = '++**##'
      do 110 j = 1, 64
          if (c(1)(j:j) .EQ. ' ')  then
            icend(1) = j-1
            GOTO 121
          endif
 110  continue
 121  do 120 j = 1, 8
          if (cside(1)(j:j) .EQ. ' ') then
            icend(2) = j-1
            GOTO 210
          endif
 120  continue
 210  write(6,*) ' '
      write(6,*) ' '
      do 200 i = 1, iclast
        if      ((i.eq.1).or.(i.eq.iclast)) then
           cleft  = cside(1)
           cright = cside(1)
         elseif ((i.eq.2).or.(i.eq.(iclast-1))) then
           cleft  = cside(4)
           cright = cside(5)
         elseif ((i.eq.3).or.(i.eq.(iclast-2))) then
           cleft  = cside(6)
           cright = cside(7)
         else
           cleft  = cside(2)
           cright = cside(3)
        endif
        write(6,*) '     ' // cleft(1:icend(2)) // c(i)(1:icend(1))
     &                     // cright(1:icend(2))
 200  continue
      if (istep.eq.1) then
         cstep = '<<<<<<<<<< This is an INTEGRATION step. >>>>>>>>>>'
       elseif (istep.eq.2) then
         cstep = '<<<<<<< This is an EVENT-GENERATION step. >>>>>>>'
       else
         write(6,*) '!!!Error in START_GRAPE!!!'
         write(6,*) ' ---> Invalid istep value(=',istep,' )'
         write(6,*) ' ---> Good-bye!'
         STOP
      endif
      write(6,*) ' '
      write(6,*)  '               ' // cstep
      write(6,*) ' '
      return
      end
