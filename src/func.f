*     File func.f : Sat Mar 18 19:44:57 2000
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
************************************************************************
      function func(x)
      implicit real*8(a-h,o-z)
      common /amjprc/jproc
      parameter ( mxdim = 50 )
      real*8   x(mxdim)
      include 'funcs.inc'
      if(jproc .eq. 1) then
          func = fncep_epee(x)
      endif
      if(jproc .eq. 2) then
          func = fncep_epmm(x)
      endif
      if(jproc .eq. 3) then
          func = fncep_eptt(x)
      endif
      if(jproc .eq. 4) then
          func = fncep_eXee(x)
      endif
      if(jproc .eq. 5) then
          func = fncep_eXmm(x)
      endif
      if(jproc .eq. 6) then
          func = fncep_eXtt(x)
      endif
      if(jproc .eq. 7) then
          func = fnceu_euee(x)
      endif
      if(jproc .eq. 8) then
          func = fnceu_eumm(x)
      endif
      if(jproc .eq. 9) then
          func = fnceu_eutt(x)
      endif
      if(jproc .eq. 10) then
          func = fnceub_eubee(x)
      endif
      if(jproc .eq. 11) then
          func = fnceub_eubmm(x)
      endif
      if(jproc .eq. 12) then
          func = fnceub_eubtt(x)
      endif
      if(jproc .eq. 13) then
          func = fnced_edee(x)
      endif
      if(jproc .eq. 14) then
          func = fnced_edmm(x)
      endif
      if(jproc .eq. 15) then
          func = fnced_edtt(x)
      endif
      if(jproc .eq. 16) then
          func = fncedb_edbee(x)
      endif
      if(jproc .eq. 17) then
          func = fncedb_edbmm(x)
      endif
      if(jproc .eq. 18) then
          func = fncedb_edbtt(x)
      endif
      if(jproc .eq. 19) then
          func = fnces_esee(x)
      endif
      if(jproc .eq. 20) then
          func = fnces_esmm(x)
      endif
      if(jproc .eq. 21) then
          func = fnces_estt(x)
      endif
      if(jproc .eq. 22) then
          func = fncesb_esbee(x)
      endif
      if(jproc .eq. 23) then
          func = fncesb_esbmm(x)
      endif
      if(jproc .eq. 24) then
          func = fncesb_esbtt(x)
      endif
      if(jproc .eq. 25) then
          func = fncec_ecee(x)
      endif
      if(jproc .eq. 26) then
          func = fncec_ecmm(x)
      endif
      if(jproc .eq. 27) then
          func = fncec_ectt(x)
      endif
      if(jproc .eq. 28) then
          func = fncecb_ecbee(x)
      endif
      if(jproc .eq. 29) then
          func = fncecb_ecbmm(x)
      endif
      if(jproc .eq. 30) then
          func = fncecb_ecbtt(x)
      endif
      if(jproc .eq. 31) then
          func = fnceb_ebee(x)
      endif
      if(jproc .eq. 32) then
          func = fnceb_ebmm(x)
      endif
      if(jproc .eq. 33) then
          func = fnceb_ebtt(x)
      endif
      if(jproc .eq. 34) then
          func = fncebb_ebbee(x)
      endif
      if(jproc .eq. 35) then
          func = fncebb_ebbmm(x)
      endif
      if(jproc .eq. 36) then
          func = fncebb_ebbtt(x)
      endif
      if(jproc .eq. 37) then
          func = fncet_etee(x)
      endif
      if(jproc .eq. 38) then
          func = fncet_etmm(x)
      endif
      if(jproc .eq. 39) then
          func = fncet_ettt(x)
      endif
      if(jproc .eq. 40) then
          func = fncetb_etbee(x)
      endif
      if(jproc .eq. 41) then
          func = fncetb_etbmm(x)
      endif
      if(jproc .eq. 42) then
          func = fncetb_etbtt(x)
      endif
      return
      end
