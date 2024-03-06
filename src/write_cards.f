***********************************************************
*     Getting Input_parameters from Control_card File     *
*                    written by T.Abe                     *
*                   on Aug. 20 in 1998                    *
***********************************************************
* Modified by T.Abe on May 23, 2002
***********************************************
      subroutine Write_cards
      implicit NONE
* -------- Argument --------
* --------------------------
* -------------------- FFREAD stuff --------------------
* ------------------------------------------------------
      include 'graepia.inc'
* ------ Local variables ------
      integer   i, N
* ============ Printing user-defined control_cards ============
      write(6,*) ' '
      write(6,1100)     '   KFLBEAM     ', KF_Lbeam
      write(6,2100)     '   EBEAM       ', P_e_beam
      write(6,2100)     '   PBEAM       ', P_p_beam
      write(6,2300)     '   EPOL        ', (Ebeam_pol(i),i=1,3)
      write(6,1100)     '   PROCESS     ', process
      write(6,1100)     '   LPAIR       ', lpair
      write(6,1100)     '   ISR         ', Isr_flag
      write(6,1100)     '   NNISR       ', NN_ISR
      write(6,1100)     '   ISRSCALE    ', ISR_scale
      write(6,2100)     '   FACTISR     ', Factor_ISR
      write(6,1100)     '   QFLV        ', qflv
      write(6,1100)     '   MERGE       ', merge
      write(6,1100)     '   QCDSCALE    ', IQCD_scale
      write(6,1100)     '   NGROUP      ', Ngroup
      write(6,1100)     '   NSET        ', Nset
      write(6,1100)     '   STRF        ', Istrf
      N = min(num_gra_flg, 5)
      write(6,1500)     '   GRAFLG      ', (jgra_flag(i),i=1,N)
      if (num_gra_flg .GT. 5) then
         N = min(num_gra_flg, 10)
         write(6,1500)  '               ', (jgra_flag(i),i=6,N)
      endif
      if (num_gra_flg .GT.10) then
         N = min(num_gra_flg, 15)
         write(6,1500)  '               ', (jgra_flag(i),i=11,N)
      endif
      if (num_gra_flg .GT.15) then
         N = min(num_gra_flg, 20)
         write(6,1500)  '               ', (jgra_flag(i),i=16,N)
      endif
      if (num_gra_flg .GT.20) then
         N = num_gra_flg   ! min(num_gra_flg, 25)
         write(6,1500)  '               ', (jgra_flag(i),i=21,N)
      endif
      write(6,*)   '======= Electroweak Dilepton Production ======='
      write(6,1100)     '  GRASEL       ', jgra_sel
      write(6,2100)     '  QEDMJWEI     ', weimj_QED
      write(6,2100)     '  Z0MJWEI      ', weimj_Z0
      write(6,2100)     '  HMJWEI       ', weimj_H
      write(6,2100)     '  HMASS        ', am_higgs
      write(6,2100)     '  HWIDTH       ', ag_higgs
      write(6,*) ' '
      write(6,1100)     '  VMMODEL      ', Model_VM
      write(6,2100)     '  VMTSLOPE     ', t_slope_VM
      write(6,100)      '  VMEEINT      ', lee_int_VM
      write(6,130)      '  VMHELI       ', (helicity_VM(i),i=-1,+1)
      write(6,*) ' '
      write(6,100)      '  JPPROD       ', lJPprod
      write(6,*)        '  JPTYPE       ', (lJPtype(i),i=1,num_jptype)
      write(6,*)        '  JPBCORR      ', (bcorr_jp(i),i=1,num_jptype)
      write(6,*)        '  JPMJWEI      ', (weimj_jp(i),i=1,num_jptype)
      write(6,*)        '  JPPHASEP     ', (phaseP_jp(i),i=1,num_jptype)
      write(6,*)        '  JPPHASEG     ', (phaseG_jp(i),i=1,num_jptype)
      write(6,*)        '  JPGRAPH      ', (lJPgra(i),i=1,num_jpgra)
      write(6,*) ' '
      write(6,100)      '  YYPROD       ', lYYprod
      write(6,*)        '  YYTYPE       ', (lYYtype(i),i=1,num_yytype)
      write(6,*)        '  YYBCORR      ', (bcorr_yy(i),i=1,num_yytype)
      write(6,*)        '  YYMJWEI      ', (weimj_yy(i),i=1,num_yytype)
      write(6,*)        '  YYPHASEP     ', (phaseP_yy(i),i=1,num_yytype)
      write(6,*)        '  YYPHASEG     ', (phaseG_yy(i),i=1,num_yytype)
      write(6,*)        '  YYGRAPH      ', (lYYgra(i),i=1,num_yygra)
      write(6,*)   '========== BASES parameters =========='
      write(6,1100)     '  ITMX1        ', num_it_grid
      write(6,1100)     '  ITMX2        ', num_it_integ
      write(6,2100)     '  ACC1         ', acc1_card
      write(6,2100)     '  ACC2         ', acc2_card
      write(6,1100)     '  NCALL        ', num_call
      write(6,1100)     '  SCALEXQ      ', Iscale_xq
      write(6,1100)     '  ISYM34       ', isym_34
      write(6,1100)     '  I34          ', ii34
      write(6,1100)     '  RESNS56      ', Ireso56
      write(6,2100)     '  NNMJ56       ', Rnn_MJ56
      write(6,2100)     '  THRES56      ', thresh56
      write(6,1100)     '  RESNS456     ', Ireso456
      write(6,2100)     '  RNN456       ', Rnn_456
      write(6,2100)     '  MAS456       ', Mas_reso456
      write(6,2100)     '  WID456       ', Wid_reso456
      write(6,1100)     '  ICOS3        ', IcosP3
      write(6,1100)     '  ICOS5        ', IcosP5
      write(6,1100)     '  NREG         ', Nregion
      write(6,1100)     '  NEPSP        ', Neps_p
      write(6,2100)     '  RSCALEMN     ', Rscale_Mn
      write(6,*)   '======================================'
      write(6,*)   '========== SPRING parameters =========='
      write(6,1100)     '  NGEN         ', Ngen
      write(6,1100)     '  NMOD         ', Nmod
      write(6,1100)     '  MXTRY        ', mxtry
      write(6,1100)     '  PSISR        ', PS_isr
      write(6,1100)     '  PSFSR        ', PS_fsr
      write(6,1100)     '  PSBRA        ', PS_bra
      write(6,1100)     '  PSSUP        ', PS_sup
      write(6,1100)     '  PYDECAY      ', PY_decay
      write(6,1100)     '  PRIPT        ', Ipri_Pt
      write(6,*)   '======================================'
      write(6,*)   '========== Output of Generated Events =========='
      write(6,100)      '  PYLIST       ', LIST_flag
      write(6,1100)     '  NLIST        ', Nlist
      write(6,100)      '  NTPYT        ', NTPYT_flag
      write(6,100)      '  NTVEC        ', NTVEC_flag
      write(6,100)      '  ASCII        ', ASC_flag
      write(6,*)   '================================================'
      write(6,*)   '================== Misc. =================='
      write(6,2100)     '  EMINISR      ', E_min_isr
      write(6,100)      '  ROTKIN       ', ROT_kinem_flag
      write(6,100)      '  ROTPY        ', ROT_pyrand_flag
      write(6,1100)     '  QELAX        ', Qela_decay
      write(6,100)      '  WRTCARDS     ', Wrt_cards
      write(6,100)      '  RNDGEN       ', LRND_flag
      write(6,1100)     '  FRAMEAMP     ', Frame_amp
      write(6,*)   '==========================================='
      write(6,*)   '=============== << Cuts >> ==============='
      write(6,2200)     '  XXRNGME      ',  (x_Range_ME(i),i=1,2)
      write(6,2200)     '  XXRNGOB      ',  (x_Range_OB(i),i=1,2)
      write(6,2200)     '  YYRNGME      ',  (y_Range_ME(i),i=1,2)
      write(6,2200)     '  YYRNGOB      ',  (y_Range_OB(i),i=1,2)
      write(6,2200)     '  Q2RNGME      ', (Q2_Range_ME(i),i=1,2)
      write(6,2200)     '  Q2RNGOB      ', (Q2_Range_OB(i),i=1,2)
      write(6,2200)     '  WWRNGME      ',  (W_Range_ME(i),i=1,2)
      write(6,2200)     '  WWRNGOB      ',  (W_Range_OB(i),i=1,2)
      write(6,*) ' ---------- '
      write(6,2200)     '  MHADCUT      ', (W_cut(i),i=1,2)
      write(6,*) ' ---------- '
      write(6,2200)     '  Q2P          ', (Q2p_cut(i),i=1,2)
      write(6,*) ' ---------- '
      write(6,2200)     '  UCUT         ', (u_cut(i),i=1,2)
      write(6,*) ' ---------- '
      write(6,2400)     '  THMIN        ', (theta_min(i),i=1,4)
      write(6,2400)     '  THMAX        ', (theta_max(i),i=1,4)
      write(6,2400)     '  EMIN         ', (E_min(i),    i=1,4)
      write(6,2400)     '  EMAX         ', (E_max(i),    i=1,4)
      write(6,2400)     '  PMIN         ', (P_min(i),    i=1,4)
      write(6,2400)     '  PMAX         ', (P_max(i),    i=1,4)
      write(6,2400)     '  PTMIN        ', (Pt_min(i),   i=1,4)
      write(6,2400)     '  PTMAX        ', (Pt_max(i),   i=1,4)
      write(6,*) ' ---------- '
      write(6,100)      '  PTMCTFLG     ', lPtMAX_cut
      write(6,2200)     '  THPTMCT      ', (the_PtMAX_cut(i),i=1,2)
      write(6,2200)     '  PTMXCT       ', (PtMAX_cut(i),i=1,2)
      write(6,*) ' ---------- '
      write(6,2200)     '  MASSLL       ', (Mass56_cut(i),i=1,2)
      write(6,2200)     '               ', (Mass56_cut(i),i=3,4)
      write(6,2200)     '  MASSELL      ', (MassELL_cut(i),i=1,2)
      write(6,2200)     '  MASSQLL      ', (MassQLL_cut(i),i=1,2)
      write(6,2200)     '               ', (MassQLL_cut(i),i=3,4)
      write(6,*)   ' ------ Hot selection ------ '
      write(6,100)      '  HOTFLG       ', lHOT_flag
      write(6,100)      '   LVETO       ', lveto
      write(6,1100)     '   IVETO       ', Iveto
      write(6,2400)     '    THEVTMIN   ', (the_veto_min(i), i=1,4)
      write(6,2400)     '    THEVTMAX   ', (the_veto_max(i), i=1,4)
      write(6,2400)     '    EVTMIN     ', (ene_veto_min(i), i=1,4)
      write(6,2400)     '    EVTMAX     ', (ene_veto_max(i), i=1,4)
      write(6,2400)     '    PTVTMIN    ', (Pt_veto_min(i),  i=1,4)
      write(6,2400)     '    PTVTMAX    ', (Pt_veto_max(i),  i=1,4)
      write(6,1100)     '  IVISI        ', Ivisi
      write(6,2400)     '    THEVMIN    ', (the_visi_min(i), i=1,4)
      write(6,2400)     '    THEVMAX    ', (the_visi_max(i), i=1,4)
      write(6,2400)     '    EVMIN      ', (ene_visi_min(i), i=1,4)
      write(6,2400)     '    EVMAX      ', (ene_visi_max(i), i=1,4)
      write(6,2400)     '    PTVMIN     ', (Pt_visi_min(i),  i=1,4)
      write(6,2400)     '    PTVMAX     ', (Pt_visi_max(i),  i=1,4)

      write(6,100)      '   L3DCUT      ', l3D_cut
      write(6,2300)     '    TH3DMIN    ', (the_3D_min(i), i=1,3)
      write(6,2300)     '    TH3DMAX    ', (the_3D_max(i), i=1,3)
      write(6,2300)     '    E3DMIN     ', (E_3D_min(i),   i=1,3)
      write(6,2300)     '    E3DMAX     ', (E_3D_max(i),   i=1,3)
      write(6,2300)     '    P3DMIN     ', (P_3D_min(i),   i=1,3)
      write(6,2300)     '    P3DMAX     ', (P_3D_max(i),   i=1,3)
      write(6,2300)     '    PT3DMIN    ', (Pt_3D_min(i),  i=1,3)
      write(6,2300)     '    PT3DMAX    ', (Pt_3D_max(i),  i=1,3)
      write(6,2200)     '    M12CUT3D   ', (Mass12_3D(i),  i=1,2)
      write(6,2200)     '    Q2CT3DME   ', (Q2_3D_ME(i),   i=1,2)
      write(6,2200)     '    Q2CT3DOB   ', (Q2_3D_OB(i),   i=1,2)
      write(6,2200)     '    WCT3DME    ', (W_3D_ME(i),   i=1,2)
      write(6,2200)     '    WCT3DOB    ', (W_3D_OB(i),   i=1,2)
      write(6,2200)     '    XCT3DME    ', (X_3D_ME(i),   i=1,2)
      write(6,2200)     '    XCT3DOB    ', (X_3D_OB(i),   i=1,2)
      write(6,2200)     '    YCT3DME    ', (Y_3D_ME(i),   i=1,2)
      write(6,2200)     '    YCT3DOB    ', (Y_3D_OB(i),   i=1,2)

      write(6,100)      '   L3ECUT      ', l3E_cut
      write(6,2300)     '    TH3EMIN    ', (the_3E_min(i), i=1,3)
      write(6,2300)     '    TH3EMAX    ', (the_3E_max(i), i=1,3)
      write(6,2300)     '    E3EMIN     ', (E_3E_min(i),   i=1,3)
      write(6,2300)     '    E3EMAX     ', (E_3E_max(i),   i=1,3)
      write(6,2300)     '    P3EMIN     ', (P_3E_min(i),   i=1,3)
      write(6,2300)     '    P3EMAX     ', (P_3E_max(i),   i=1,3)
      write(6,2300)     '    PT3EMIN    ', (Pt_3E_min(i),  i=1,3)
      write(6,2300)     '    PT3EMAX    ', (Pt_3E_max(i),  i=1,3)
      write(6,2200)     '    M12CUT3E   ', (Mass12_3E(i),  i=1,2)
      write(6,2200)     '    Q2CT3EME   ', (Q2_3E_ME(i),   i=1,2)
      write(6,2200)     '    Q2CT3EOB   ', (Q2_3E_OB(i),   i=1,2)
      write(6,2200)     '    WCT3EME    ', (W_3E_ME(i),   i=1,2)
      write(6,2200)     '    WCT3EOB    ', (W_3E_OB(i),   i=1,2)
      write(6,2200)     '    XCT3EME    ', (X_3E_ME(i),   i=1,2)
      write(6,2200)     '    XCT3EOB    ', (X_3E_OB(i),   i=1,2)
      write(6,2200)     '    YCT3EME    ', (Y_3E_ME(i),   i=1,2)
      write(6,2200)     '    YCT3EOB    ', (Y_3E_OB(i),   i=1,2)

      write(6,100)      '   L3FCUT      ', l3F_cut
      write(6,2300)     '    TH3FMIN    ', (the_3F_min(i), i=1,3)
      write(6,2300)     '    TH3FMAX    ', (the_3F_max(i), i=1,3)
      write(6,2300)     '    E3FMIN     ', (E_3F_min(i),   i=1,3)
      write(6,2300)     '    E3FMAX     ', (E_3F_max(i),   i=1,3)
      write(6,2300)     '    P3FMIN     ', (P_3F_min(i),   i=1,3)
      write(6,2300)     '    P3FMAX     ', (P_3F_max(i),   i=1,3)
      write(6,2300)     '    PT3FMIN    ', (Pt_3F_min(i),  i=1,3)
      write(6,2300)     '    PT3FMAX    ', (Pt_3F_max(i),  i=1,3)
      write(6,2200)     '    M12CUT3F   ', (Mass12_3F(i),  i=1,2)
      write(6,2200)     '    Q2CT3FME   ', (Q2_3F_ME(i),   i=1,2)
      write(6,2200)     '    Q2CT3FOB   ', (Q2_3F_OB(i),   i=1,2)
      write(6,2200)     '    WCT3FME    ', (W_3F_ME(i),   i=1,2)
      write(6,2200)     '    WCT3FOB    ', (W_3F_OB(i),   i=1,2)
      write(6,2200)     '    XCT3FME    ', (X_3F_ME(i),   i=1,2)
      write(6,2200)     '    XCT3FOB    ', (X_3F_OB(i),   i=1,2)
      write(6,2200)     '    YCT3FME    ', (Y_3F_ME(i),   i=1,2)
      write(6,2200)     '    YCT3FOB    ', (Y_3F_OB(i),   i=1,2)


      write(6,100)      '   L2EVISIA    ', l2e_visiA_flag
      write(6,2300)     '     THE2E     ', (the_2e(i),i=1,3)
      write(6,2200)     '     EMIN2E    ', (E_min_2e(i),i=1,2)
      write(6,100)      '   L2EVISIB    ', l2e_visiB_flag
      write(6,2500)     '     THE2EB    ', (the_2eB(i),i=1,5)
      write(6,2400)     '     EMIN2EB   ', (E_min_2eB(i),i=1,4)
      write(6,100)      '   L3EVISI     ', l3e_visi_flag
      write(6,2400)     '     THE3E     ', (the_3e(i),i=1,4)
      write(6,2100)     '     EMIN3E    ', E_min_3e
      write(6,100)      '   LPTVISI     ', lPt_visi_flag
      write(6,2300)     '     THEMINPT  ', (the_min_pt(i),i=1,3)
      write(6,2300)     '     THEMAXPT  ', (the_max_pt(i),i=1,3)
      write(6,2300)     '     PTMINPT   ', (Pt_min_pt(i),i=1,3)
      write(6,100)      '   LSCATL      ', lscattL_flag
      write(6,2200)     '     THESCATL  ', (theta_scattL(i),i=1,2)
      write(6,2200)     '     ESCATL    ', (E_scattL(i),    i=1,2)
      write(6,2200)     '     THEPRODL  ', ( theta_prodL(i),i=1,2)
      write(6,2200)     '     PPRODL    ', ( P_prodL(i),    i=1,2)
      write(6,*)   ' --------------------------- '
      write(6,*)   '=========================================='
      write(6,*)   '========== Parameters in Quasi-elastic =========='
      write(6,2100)     '   RK          ', rk
      write(6,2100)     '   ANBD        ', A_NBD
      write(6,2100)     '   BNBD        ', B_NBD
      write(6,2100)     '   CNBD        ', C_NBD
      write(6,2100)     '   PUQ         ', P_u
      write(6,2100)     '   PDQ         ', P_d
      write(6,2100)     '   PSQ         ', P_s
      write(6,2100)     '   PUD         ', P_ud_bar
      write(6,2100)     '   SUPPI0      ', Supp_pi0
      write(6,2100)     '   ASLOPE      ', A_slope
      write(6,2100)     '   BSLOPE      ', B_slope
      write(6,2100)     '   CSLOPE      ', C_slope
      write(6,2100)     '   REDS        ', red_slope_s
      write(6,*)   '================================================='
 100  format(A15, L10)
 130  format(A15, 3L5)
 1100 format(A15,   I10)
 1500 format(A15, 5(I10,1X))
 2100 format(A15,   G10.4)
 2200 format(A15, 2(G10.4,2X))
 2300 format(A15, 3(G10.4,2X))
 2400 format(A15, 4(G10.4,2X))
 2500 format(A15, 5(G10.4,2X))
* =============================================================
      return
      end
