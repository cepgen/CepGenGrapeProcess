*     File amclr.f : Sat Mar 18 19:44:57 2000
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
************************************************************************
      subroutine amclr
      implicit real*8(a-h,o-z)
      common /amjprc/jproc
      if(jproc .eq. 1) then
          call aep_epeemclr
      endif
      if(jproc .eq. 2) then
          call aep_epmmmclr
      endif
      if(jproc .eq. 3) then
          call aep_epttmclr
      endif
      if(jproc .eq. 4) then
          call aep_eXeemclr
      endif
      if(jproc .eq. 5) then
          call aep_eXmmmclr
      endif
      if(jproc .eq. 6) then
          call aep_eXttmclr
      endif
      if(jproc .eq. 7) then
          call aeu_eueemclr
      endif
      if(jproc .eq. 8) then
          call aeu_eummmclr
      endif
      if(jproc .eq. 9) then
          call aeu_euttmclr
      endif
      if(jproc .eq. 10) then
          call aeub_eubeemclr
      endif
      if(jproc .eq. 11) then
          call aeub_eubmmmclr
      endif
      if(jproc .eq. 12) then
          call aeub_eubttmclr
      endif
      if(jproc .eq. 13) then
          call aed_edeemclr
      endif
      if(jproc .eq. 14) then
          call aed_edmmmclr
      endif
      if(jproc .eq. 15) then
          call aed_edttmclr
      endif
      if(jproc .eq. 16) then
          call aedb_edbeemclr
      endif
      if(jproc .eq. 17) then
          call aedb_edbmmmclr
      endif
      if(jproc .eq. 18) then
          call aedb_edbttmclr
      endif
      if(jproc .eq. 19) then
          call aes_eseemclr
      endif
      if(jproc .eq. 20) then
          call aes_esmmmclr
      endif
      if(jproc .eq. 21) then
          call aes_esttmclr
      endif
      if(jproc .eq. 22) then
          call aesb_esbeemclr
      endif
      if(jproc .eq. 23) then
          call aesb_esbmmmclr
      endif
      if(jproc .eq. 24) then
          call aesb_esbttmclr
      endif
      if(jproc .eq. 25) then
          call aec_eceemclr
      endif
      if(jproc .eq. 26) then
          call aec_ecmmmclr
      endif
      if(jproc .eq. 27) then
          call aec_ecttmclr
      endif
      if(jproc .eq. 28) then
          call aecb_ecbeemclr
      endif
      if(jproc .eq. 29) then
          call aecb_ecbmmmclr
      endif
      if(jproc .eq. 30) then
          call aecb_ecbttmclr
      endif
      if(jproc .eq. 31) then
          call aeb_ebeemclr
      endif
      if(jproc .eq. 32) then
          call aeb_ebmmmclr
      endif
      if(jproc .eq. 33) then
          call aeb_ebttmclr
      endif
      if(jproc .eq. 34) then
          call aebb_ebbeemclr
      endif
      if(jproc .eq. 35) then
          call aebb_ebbmmmclr
      endif
      if(jproc .eq. 36) then
          call aebb_ebbttmclr
      endif
      if(jproc .eq. 37) then
          call aet_eteemclr
      endif
      if(jproc .eq. 38) then
          call aet_etmmmclr
      endif
      if(jproc .eq. 39) then
          call aet_etttmclr
      endif
      if(jproc .eq. 40) then
          call aetb_etbeemclr
      endif
      if(jproc .eq. 41) then
          call aetb_etbmmmclr
      endif
      if(jproc .eq. 42) then
          call aetb_etbttmclr
      endif
      return
      end
