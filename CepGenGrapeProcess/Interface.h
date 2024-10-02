/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGenGrapeProcess_Interface_h
#define CepGenGrapeProcess_Interface_h

#include <array>
#include <complex>
#include <cstddef>

struct InitialStateCuts {
  std::array<float, 2> x_range, dum1, y_range, dum2, q2_range, dum3, w_range, dum4;
};

struct FinalStateCuts {
  int lcut;
  std::array<float, 3> the_min, the_max, e_min, e_max, p_min, p_max, pt_min, pt_max;
  std::array<float, 2> mass12, q2_me, q2_ob, w_me, w_ob, x_me, x_ob, y_me, y_ob;
};

template <size_t Ntype, size_t Ngra>
struct VectorMesonInfo {
  int lprod;
  std::array<int, Ntype> ltype;
  std::array<int, Ngra> lgra;
  std::array<float, Ntype> bcorr, weimj, phasep, phaseg;
};

#ifdef __cplusplus
extern "C" {
#endif
void amparm_();
double flux_factor_(double& s, double& mass2_1, double& mass2_2, int& jump);
double func_(double x[50]);
void get_graph_flag_();
void get_proc_();
void gfinit_();
void k2nit_();
void prmass_();
void start_grape_(int&);
void write_cards_();

extern struct {
  int jproc;
} amjprc_;

static constexpr size_t mxdim = 50;
extern struct {
  std::array<double, mxdim> xl, xu;
  int ndim, nwild;
  std::array<int, mxdim> ig;
  int ncall;
} bparm1_;

extern struct {
  int isrpol;
} cmn_isrpol_;

// masses and width of particles
extern struct {
  std::array<double, 38> mass;
} smmass_;

extern struct {
  std::array<double, 26> width;
} smgmma_;

extern struct {
  double pi, pi2, rad, gevpb, alpha, alphas, alpha0;
} smcnst_;

extern struct {
  std::array<double, 2> coscut[4], engyct[4], amasct[6], aresns[4];
  double opncut, swapm2;
} cut001_4f_;

extern struct {
  int kf_lbeam, ipol_ebeam;
  float p_e_beam, p_p_beam;
  std::array<float, 3> ebeam_pol;
  int lpol_ebeam;
} gep_beam_;

extern struct {
  std::array<float, 2> q2p_cut;
  std::array<float, 4> theta_min, theta_max, e_min, e_max, p_min, p_max, pt_min, pt_max;
  int lptmax_cut;
  std::array<float, 2> ptmax_cut, the_ptmax_cut, e_scattl, theta_scattl, p_prodl, theta_prodl;
  std::array<float, 2> mass56_cut[2];
  std::array<float, 2> massell_cut;
  std::array<float, 2> massqll_cut[2];
  std::array<float, 4> the_visi_min, the_visi_max, ene_visi_min, ene_visi_max, pt_visi_min, pt_visi_max;
  int ivisi;
  std::array<float, 3> the_2e;
  std::array<float, 2> e_min_2e;
  std::array<float, 4> the_3e;
  float e_min_3e;
  std::array<float, 3> the_min_pt, the_max_pt, pt_min_pt;
  std::array<float, 5> the_2eb;
  std::array<float, 4> e_min_2eb;
  int lhot_flag, l2e_visia_flag, l2e_visib_flag, l3e_visi_flag, lpt_visi_flag, lscattl_flag;
  std::array<float, 2> u_cut, w_cut;
  int lveto, iveto;
  std::array<float, 4> the_veto_min, the_veto_max, ene_veto_min, ene_veto_max, pt_veto_min, pt_veto_max;
  InitialStateCuts cut_me, cut_ob;
  FinalStateCuts cut_3d, cut_3e, cut_3f;
} gep_cut_;

extern struct {
  int jgra_sel;
  float weimj_qed, weimj_z0, weimj_h, am_higgs, ag_higgs;
} gep_ew_;

extern struct {
  std::array<double, 4> pe_str[12];
} gep_frm_;

// J/psi(psi) production
static constexpr size_t num_jptype = 6, num_jpgra = 2;
extern VectorMesonInfo<num_jptype, num_jpgra> gep_jp_;

extern struct {
  int isym_34, ii34, iscale_xq, ireso56, ireso456, nregion, icosp3, icosp5, neps_p;
  float thresh56, rnn_mj56, mas_reso456, wid_reso456, rscale_mn, rnn_456;
} gep_kine_;

extern struct {
  double p1_lab, p2_lab, e1_lab, e2_lab;
  std::array<double, 3> pcms_lab, ecms_lab, gammacms_lab, betgamcms_lab;
  std::array<double, 4> vec_isr;
} gep_lab_;

extern struct {
  int icount, print_flag, nmod, nlist, frame_amp, qela_decay;
  int list_flag, ntpyt_flag, asc_flag, ntvec_flag, wrt_cards, lrnd_flag;
} gep_misc_;

extern struct {
  int process, lpair, qflv, isr_flag, nn_isr, isr_scale;
  float factor_isr;
  std::array<int, 10> jgra_flag;
  int ngroup, nset, istrf;
  float wmin, wmax;
  int merge, iqcd_scale;
} gep_proc_;

// quasi-elastic process
extern struct {
  float rk, a_nbd, b_nbd, c_nbd, p_u, p_d, p_s, p_ud_bar, supp_pi0, a_slope, b_slope, c_slope, red_slope_s;
} gep_qela_;

// vector meson production
extern struct {
  int model_vm;
  float t_slope_vm;
  int lee_int_vm;
  std::array<int, 3> helicity_vm;
} gep_vm_;

// Upsilon production
static constexpr size_t num_yytype = 6, num_yygra = 2;
extern VectorMesonInfo<num_yytype, num_yygra> gep_yy_;

extern struct {
  std::array<double, 3> s, w;
  double fact;
} kinem1_;

extern struct {
  std::array<int, 10> kcharg, kfcode;
} kminfo_;

// regions multiplicity
extern struct {
  int mxreg;
} kmreg_;

// masses of external particles
extern struct {
  std::array<double, 10> amass1, amass2;
} kmmass_;

// phase of propagators
extern struct {
  std::complex<double> pphase;
} smpphs_;

extern struct {
  std::complex<double> cqed, cqcd, czww, caww, cggg, cwwaa, cwwza, cwwzz, cwwww, cgggg, cwhm, cwhp, cwym, cwyp, czpm,
      capm, czhy, chww, chzz, cwzm, cwam, cwzp, cwap, cwwhh, czzhh, cwzhm, cwzhp, cwahm, cwahp, cwwyy, czzyy, cwzym,
      cwzyp, cwaym, cwayp, cwwpm, czzpm, caapm, czapm, chpm, chyy, chhh, chhhh, cyyyy, cpmyy, chhpm, chhyy, cpmpm;
  std::array<std::complex<double>, 2> cwln[3], cwnl[3], cwdu[3], cwud[3], call[3], cauu[3], cadd[3], capp, capx,
      cznn[3], czll[3], czuu[3], czdd[3], cguu[3], cgdd[3], cmdu[3], cpud[3], chll[3], chuu[3], chdd[3], cyuu[3],
      cydd[3];
  std::complex<double> cwczp, cwcmz, cwcap, cwcma, cwczm, cwcpz, cwcam, cwcpa, czcmm, czcpp, cacmm, cacpp, cgcgg, cpczp,
      cpcap, cpcmz, cmczm, cmcam, cmcpz, cycmm, cycpp, chcmm, chcpp, chczz, c1zww, c1aww, c1ggg, c1wwaa, c1wwza, c1wwzz,
      c1wwww, c1gggg, c1whm, c1whp, c1wym, c1wyp, c1zpm, c1apm, c1zhy, c1hww, c1hzz, c1wzm, c1wam, c1wzp, c1wap, c1wwhh,
      c1zzhh, c1wzhm, c1wzhp, c1wahm, c1wahp, c1wwyy, c1zzyy, c1wzym, c1wzyp, c1waym, c1wayp, c1wwpm, c1zzpm, c1aapm,
      c1zapm, c1hpm, c1hyy, c1hhh, c1hhhh, c1yyyy, c1pmyy, c1hhpm, c1hhyy, c1pmpm;
  std::array<std::complex<double>, 2> c1wln[3], c1wnl[3], c1wdu[3], c1wud[3], c1all[3], c1auu[3], c1add[3], c1app,
      c1apx, c1znn[3], c1zll[3], c1zuu[3], c1zdd[3], c1guu[3], c1gdd[3], c1mdu[3], c1pud[3], c1hll[3], c1huu[3],
      c1hdd[3], c1yuu[3], c1ydd[3];
  std::complex<double> c1wczp, c1wcmz, c1wcap, c1wcma, c1wczm, c1wcpz, c1wcam, c1wcpa, c1zcmm, c1zcpp, c1acmm, c1acpp,
      c1gcgg, c1pczp, c1pcap, c1pcmz, c1mczm, c1mcam, c1mcpz, c1ycmm, c1ycpp, c1hcmm, c1hcpp, c1hczz;
  std::array<std::complex<double>, 4> c1aa, c1ww, c1za, c1zz, c1gg;
  std::array<std::complex<double>, 2> c1pp, c1yy, c1xx;
  std::complex<double> c1wxpm, c1wxmp, c1zy, c1ay;
  std::array<std::complex<double>, 8> c1nene, c1nmnm, c1ntnt, c1elel, c1mumu, c1tata, c1uquq, c1cqcq, c1tqtq, c1dqdq,
      c1sqsq, c1bqbq, c1cuquq, c1ccqcq, c1ctqtq, c1cdqdq, c1csqsq, c1cbqbq;
  std::complex<double> c1czcz, c1caca, c1cmcm, c1cpcp, c1cgcg, c1ahy, c1hza, c1haa, c1zahh, c1aahh, c1zayy, c1aayy;
  std::array<std::complex<double>, 2> c1ane, c1anm, c1ant;
} smcplc_;

extern struct {
  std::array<double, 4> vec[10];
} sp4vec_;

// for PDFlib
extern struct {
  double xmin, xmax, q2min, q2max;
} w50513_;
#ifdef __cplusplus
}
#endif
#endif
