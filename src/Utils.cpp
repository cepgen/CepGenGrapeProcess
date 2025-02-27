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

#include <CepGen/Core/Exception.h>
#include <CepGen/Modules/RandomGeneratorFactory.h>
#include <CepGen/Modules/StructureFunctionsFactory.h>
#include <CepGen/Physics/Constants.h>
#include <CepGen/Physics/Utils.h>
#include <CepGen/StructureFunctions/Parameterisation.h>
#include <CepGen/Utils/RandomGenerator.h>

#include <cmath>

#include "CepGenGrapeProcess/Interface.h"
#include "CepGenGrapeProcess/Types.h"
#include "CepGenGrapeProcess/Utils.h"

using namespace cepgen;
using namespace std::string_literals;

namespace grape {
  int convertLeptonId(pdgid_t pdgid) {
    switch (pdgid) {
      case PDG::electron:
        return 1;
      case PDG::muon:
        return 2;
      case PDG::tau:
        return 3;
      default:
        throw CG_FATAL("Grape:convertLeptonId") << "Invalid lepton pair: " << pdgid << " is not (yet) supported.";
    }
  }

  std::string processName(int jproc) { return "invalid process"s; }

  MiscBlockFiller::MiscBlockFiller(const ParametersList& params) : BlockFiller(params) {
    gep_misc_.icount = 0;
    gep_misc_.print_flag = steer<bool>("debug");
    gep_misc_.nmod = steer<int>("NMOD");
    gep_misc_.nlist = steer<int>("NLIST");
    gep_misc_.frame_amp = steer<int>("FRAMEAMP");
    gep_misc_.qela_decay = steer<int>("QELAX");
    gep_misc_.list_flag = steer<bool>("PYLIST");
    gep_misc_.ntpyt_flag = steer<bool>("NTPYT");
    gep_misc_.asc_flag = steer<bool>("ASCII");
    gep_misc_.ntvec_flag = steer<bool>("NTVEC");
    gep_misc_.wrt_cards = steer<bool>("WRTCARDS");
    gep_misc_.lrnd_flag = steer<bool>("RNDGEN");
  }

  ParametersDescription MiscBlockFiller::description() {
    auto desc = ParametersDescription();
    desc.add<int>("NMOD", 1000).setDescription("event number printout period");
    desc.add<int>("NLIST", 10).setDescription("PYJETS event printout period");
    desc.add<int>("FRAMEAMP", 1).setDescription("(1=ep, 2=gp)");
    desc.add<int>("QELAX", 1);
    desc.add<bool>("PYLIST", false);
    desc.add<bool>("NTPYT", false).setDescription("ntuple event generation");
    desc.add<bool>("ASCII", false);
    desc.add<bool>("NTVEC", false);
    desc.add<bool>("WRTCARDS", false);
    desc.add<bool>("RNDGEN", false);
    return desc;
  }

  ProcBlockFiller::ProcBlockFiller(const ParametersList& params) : BlockFiller(params) {
    gep_proc_.process = (int)steerAs<int, ProtonMode>("PROCESS");
    gep_proc_.lpair = convertLeptonId(steer<ParticleProperties>("LPAIR").pdgid);
    gep_proc_.qflv = (int)steer<ParticleProperties>("QFLV").pdgid;
    gep_proc_.isr_flag = (int)steerAs<int, LeptonISRMode>("ISR");
    gep_proc_.nn_isr = steer<int>("NNISR");
    gep_proc_.isr_scale = steer<double>("ISRSCALE");
    gep_proc_.factor_isr = steer<double>("FACTISR");
    fillArray<int>("GRAFLG", gep_proc_.jgra_flag);
    gep_proc_.ngroup = steer<int>("NGROUP");
    gep_proc_.nset = steer<int>("NSET");
    gep_proc_.istrf = steer<int>("STRF");
    gep_proc_.merge = (int)steerAs<int, DISMergeMode>("MERGE");
    gep_proc_.iqcd_scale = steer<int>("QCDSCALE");
  }

  ParametersDescription ProcBlockFiller::description() {
    auto desc = ParametersDescription();
    desc.addAs<int, ProtonMode>("PROCESS", ProtonMode::elastic).setDescription("proton vertex mode");
    desc.addAs<int, pdgid_t>("LPAIR", PDG::muon).setDescription("lepton pair considered");
    desc.addAs<int, pdgid_t>("QFLV", PDG::up).setDescription("scattered quark in DIS");
    desc.addAs<int, LeptonISRMode>("ISR", LeptonISRMode::sfMethod).setDescription("ISR mode for the lepton beam");
    desc.addAs<int, DISMergeMode>("MERGE", DISMergeMode::off)
        .setDescription("merging mode in the DIS process")
        .allow((int)DISMergeMode::off)
        .allow((int)DISMergeMode::u_ubar_d_dbar)
        .allow((int)DISMergeMode::u_ubar_d_dbar_s_sbar)
        .allow((int)DISMergeMode::u_ubar_d_dbar_s_sbar_c_cbar)
        .allow((int)DISMergeMode::u_ubar_d_dbar_s_sbar_c_cbar_b_bbar)
        .allow((int)DISMergeMode::u_c)
        .allow((int)DISMergeMode::ubar_cbar)
        .allow((int)DISMergeMode::d_s)
        .allow((int)DISMergeMode::dbar_sbar)
        .allow((int)DISMergeMode::d_s_b)
        .allow((int)DISMergeMode::dbar_sbar_bbar);
    ;
    desc.add<int>("NNISR", -1);
    desc.add<double>("ISRSCALE", 1.);
    desc.add<double>("FACTISR", 1.);
    desc.add<std::vector<double> >("GRAFLG", std::vector<double>(gep_proc_.jgra_flag.size(), 0.));
    desc.add<int>("NGROUP", 4).setDescription("PDF set author group from PDFlib");
    desc.add<int>("NSET", 32).setDescription("PDF set from PDFlib");
    desc.add<int>("STRF", 1);
    desc.add<int>("QCDSCALE", 1);
    return desc;
  }

  EWBlockFiller::EWBlockFiller(const ParametersList& params) : BlockFiller(params) {
    gep_ew_.jgra_sel = (int)steerAs<int, EwkProductionMode>("GRASEL");
    gep_ew_.weimj_qed = steer<double>("CNTMJWEI");
    gep_ew_.weimj_z0 = steer<double>("Z0MJWEI");
    gep_ew_.weimj_h = steer<double>("H0MJWEI");
    if (const auto mh = steer<double>("H0MASS"); mh >= 0.)
      gep_ew_.am_higgs = mh;
    else
      gep_ew_.am_higgs = PDG::get().mass(25);
    if (const auto wh = steer<double>("H0WIDTH"); wh >= 0.)
      gep_ew_.ag_higgs = wh;
    else
      gep_ew_.ag_higgs = PDG::get().width(25);
  }

  ParametersDescription EWBlockFiller::description() {
    auto desc = ParametersDescription();
    desc.setDescription("Electroweak dilepton production parameters");
    desc.addAs<int, EwkProductionMode>("GRASEL", EwkProductionMode::qed)
        .setDescription(
            "Feynman diagram selection (1: BH, 2: BH+ee interf., 3: QED w/o Z0, 3: EW, 13: QED Compton, 14: Z0)");
    desc.add<double>("CNTMJWEI", 1.);
    desc.add<double>("Z0MJWEI", 0.);
    desc.add<double>("H0MJWEI", 0.);
    desc.add<double>("H0MASS", -1.).setDescription("Higgs boson mass in GeV");
    desc.add<double>("H0WIDTH", -1.).setDescription("Higgs boson width in GeV");
    return desc;
  }

  VMBlocksFiller::VMBlocksFiller(const ParametersList& params) : BlockFiller(params) {
    gep_vm_.model_vm = steer<int>("VMMODEL");
    gep_vm_.t_slope_vm = steer<double>("VMTSLOPE");
    gep_vm_.lee_int_vm = steer<bool>("VMEEINT");
    fillArray<int>("VMHELI", gep_vm_.helicity_vm);
    gep_jp_.lprod = steer<bool>("JPPROD");
    fillArray<int>("JPGRAPH", gep_jp_.lgra);
    fillArray<int>("JPTYPE", gep_jp_.ltype);
    fillArray<double>("JPBCORR", gep_jp_.bcorr);
    fillArray<double>("JPMJWEI", gep_jp_.weimj);
    fillArray<double>("JPPHASEP", gep_jp_.phasep);
    fillArray<double>("JPPHASEG", gep_jp_.phaseg);
    gep_yy_.lprod = steer<bool>("YYPROD");
    fillArray<int>("YYGRAPH", gep_yy_.lgra);
    fillArray<int>("YYTYPE", gep_yy_.ltype);
    fillArray<double>("YYBCORR", gep_yy_.bcorr);
    fillArray<double>("YYMJWEI", gep_yy_.weimj);
    fillArray<double>("YYPHASEP", gep_yy_.phasep);
    fillArray<double>("YYPHASEG", gep_yy_.phaseg);
  }

  ParametersDescription VMBlocksFiller::description() {
    auto desc = ParametersDescription();
    desc.setDescription("Vector meson production parameters");
    desc.add<int>("VMMODEL", 0);
    desc.add<double>("VMTSLOPE", 0.);
    desc.add<bool>("VMEEINT", true);
    desc.add<std::vector<int> >("VMHELI", {1, 1, 1});
    desc.add<bool>("JPPROD", false);
    desc.add<std::vector<int> >("JPGRAPH", std::vector<int>(num_jpgra, 0));
    desc.add<std::vector<int> >("JPTYPE", std::vector<int>(num_jptype, 0));
    desc.add<std::vector<double> >("JPBCORR", std::vector<double>(num_jptype, 0.));
    desc.add<std::vector<double> >("JPMJWEI", std::vector<double>(num_jptype, 0.));
    desc.add<std::vector<double> >("JPPHASEP", std::vector<double>(num_jptype, 0.));
    desc.add<std::vector<double> >("JPPHASEG", std::vector<double>(num_jptype, 0.));
    desc.add<bool>("YYPROD", false);
    desc.add<std::vector<int> >("YYGRAPH", std::vector<int>(num_yygra, 0));
    desc.add<std::vector<int> >("YYTYPE", std::vector<int>(num_yytype, 0));
    desc.add<std::vector<double> >("YYBCORR", std::vector<double>(num_yytype, 0.));
    desc.add<std::vector<double> >("YYMJWEI", std::vector<double>(num_yytype, 0.));
    desc.add<std::vector<double> >("YYPHASEP", std::vector<double>(num_yytype, 0.));
    desc.add<std::vector<double> >("YYPHASEG", std::vector<double>(num_yytype, 0.));
    return desc;
  }

  CutsBlockFiller::CutsBlockFiller(const ParametersList& params) : BlockFiller(params) {
    saveLimits(steer<Limits>("XXRNGME"), gep_cut_.cut_me.x_range);
    saveLimits(steer<Limits>("YYRNGME"), gep_cut_.cut_me.y_range);
    saveLimits(steer<Limits>("Q2RNGME"), gep_cut_.cut_me.q2_range);
    saveLimits(steer<Limits>("WWRNGME"), gep_cut_.cut_me.w_range);
    saveLimits(steer<Limits>("XXRNGOB"), gep_cut_.cut_ob.x_range);
    saveLimits(steer<Limits>("YYRNGOB"), gep_cut_.cut_ob.y_range);
    saveLimits(steer<Limits>("Q2RNGOB"), gep_cut_.cut_ob.q2_range);
    saveLimits(steer<Limits>("WWRNGOB"), gep_cut_.cut_ob.w_range);
    saveLimits(steer<Limits>("Q2P"), gep_cut_.q2p_cut);
    fillLimitsFromPair<double>(gep_cut_.theta_min, gep_cut_.theta_max, "THMIN", "THMAX");
    fillLimitsFromPair<double>(gep_cut_.e_min, gep_cut_.e_max, "EMIN", "EMAX");
    fillLimitsFromPair<double>(gep_cut_.p_min, gep_cut_.p_max, "PMIN", "PMAX");
    fillLimitsFromPair<double>(gep_cut_.pt_min, gep_cut_.pt_max, "PTMIN", "PTMAX");
    gep_cut_.lptmax_cut = steer<bool>("PTMCTFLG");
    saveLimits(steer<Limits>("PTMXCT"), gep_cut_.ptmax_cut);
    saveLimits(steer<Limits>("THPTMCT"), gep_cut_.the_ptmax_cut);
    saveLimits(steer<Limits>("MASSELL"), gep_cut_.massell_cut);
    {
      size_t i = 0;
      const auto& massll = steer<std::vector<Limits> >("MASSLL");
      for (const auto& lim : massll) {
        if (i >= 2)
          throw CG_FATAL("grape:CutsBlockFiller")
              << "Invalid number of cuts defined for MASSLL parameter ; expecting at most 2, given: " << massll << ".";
        saveLimits(lim, gep_cut_.mass56_cut[i++]);
      }
    }
    {
      size_t i = 0;
      const auto& massqll = steer<std::vector<Limits> >("MASSQLL");
      for (const auto& lim : massqll) {
        if (i >= 2)
          throw CG_FATAL("grape:CutsBlockFiller")
              << "Invalid number of cuts defined for MASSQLL parameter ; expecting at most 2, given: " << massqll
              << ".";
        saveLimits(lim, gep_cut_.massqll_cut[i++]);
      }
    }
    gep_cut_.ivisi = steer<int>("IVISI");
    fillLimitsFromPair<double>(gep_cut_.the_visi_min, gep_cut_.the_visi_max, "THEVMIN", "THEVMAX");
    fillLimitsFromPair<double>(gep_cut_.ene_visi_min, gep_cut_.ene_visi_max, "EVMIN", "EVMAX");
    fillLimitsFromPair<double>(gep_cut_.pt_visi_min, gep_cut_.pt_visi_max, "PTVMIN", "PTVMAX");
    gep_cut_.lveto = steer<bool>("LVETO");
    gep_cut_.iveto = steer<int>("IVETO");
    fillLimitsFromPair<double>(gep_cut_.the_veto_min, gep_cut_.the_veto_max, "THEVTMIN", "THEVTMAX");
    fillLimitsFromPair<double>(gep_cut_.ene_veto_min, gep_cut_.ene_veto_max, "EVTMIN", "EVTMAX");
    fillLimitsFromPair<double>(gep_cut_.pt_veto_min, gep_cut_.pt_veto_max, "PTVTMIN", "PTVTMAX");
    for (auto& cuts_vs_name : std::map<std::string, FinalStateCuts*>{
             {"D", &gep_cut_.cut_3d}, {"E", &gep_cut_.cut_3e}, {"F", &gep_cut_.cut_3f}}) {
      const auto& name = cuts_vs_name.first;
      auto* cut = cuts_vs_name.second;
      cut->lcut = steer<bool>("L3"s + name + "CUT");
      fillLimitsFromPair<double>(cut->the_min, cut->the_max, "TH3"s + name + "MIN", "TH3"s + name + "MIN");
      fillLimitsFromPair<double>(cut->e_min, cut->e_max, "E3"s + name + "MIN", "E3"s + name + "MIN");
      fillLimitsFromPair<double>(cut->p_min, cut->p_max, "P3"s + name + "MIN", "P3"s + name + "MIN");
      fillLimitsFromPair<double>(cut->pt_min, cut->pt_max, "PT3"s + name + "MIN", "PT3"s + name + "MIN");
      fillLimitsFromValue<double>(cut->mass12, "M12CUT3"s + name, LimitSide::min);
      fillLimitsFromValue<double>(cut->q2_me, "Q2CT3"s + name + "ME", LimitSide::min);
      fillLimitsFromValue<double>(cut->q2_ob, "Q2CT3"s + name + "OB", LimitSide::min);
      fillLimitsFromValue<double>(cut->w_me, "WCT3"s + name + "ME", LimitSide::min);
      fillLimitsFromValue<double>(cut->w_ob, "WCT3"s + name + "OB", LimitSide::min);
      saveLimits(steer<Limits>("XCT3"s + name + "ME"), cut->x_me);
      saveLimits(steer<Limits>("XCT3"s + name + "OB"), cut->x_ob);
      saveLimits(steer<Limits>("YCT3"s + name + "ME"), cut->y_me);
      saveLimits(steer<Limits>("YCT3"s + name + "OB"), cut->y_ob);
    }
    gep_cut_.lhot_flag = steer<bool>("HOTFLG");
    gep_cut_.l2e_visia_flag = steer<bool>("L2EVISIA");
    fillArray<double>("THE2E", gep_cut_.the_2e);
    fillArray<double>("EMIN2E", gep_cut_.e_min_2e);
    gep_cut_.l2e_visib_flag = steer<bool>("L2EVISIB");
    fillArray<double>("THE2EB", gep_cut_.the_2eb);
    fillArray<double>("EMIN2EB", gep_cut_.e_min_2eb);
    gep_cut_.l3e_visi_flag = steer<bool>("L3EVISI");
    fillArray<double>("THE3E", gep_cut_.the_3e);
    gep_cut_.e_min_3e = steer<double>("EMIN3E");
    gep_cut_.lpt_visi_flag = steer<bool>("LPTVISI");
    fillLimitsFromPair<double>(gep_cut_.the_min_pt, gep_cut_.the_max_pt, "THEMINPT", "THEMAXPT");
    fillLimitsFromValue<double>(gep_cut_.pt_min_pt, "PTMINPT", LimitSide::min);
    gep_cut_.lscattl_flag = steer<bool>("LSCATL");
    saveLimits(steer<Limits>("ESCATL"), gep_cut_.e_scattl);
    saveLimits(steer<Limits>("PPRODL"), gep_cut_.p_prodl);
    saveLimits(steer<Limits>("THESCATL"), gep_cut_.theta_scattl);
    saveLimits(steer<Limits>("PPRODL"), gep_cut_.theta_prodl);
    saveLimits(steer<Limits>("UCUT"), gep_cut_.u_cut);
    saveLimits(steer<Limits>("MHAD"), gep_cut_.w_cut);
  }

  ParametersDescription CutsBlockFiller::description() {
    auto desc = ParametersDescription();
    for (const auto& name : {"ME"s, "OB"s}) {
      desc.add<Limits>("XXRNG" + name, {-1., 1.});
      desc.add<Limits>("YYRNG" + name, {-1., 1.});
      desc.add<Limits>("Q2RNG" + name, {-1.})
          .setDescription("negative squared mom.transfer range at electron vertex ("s +
                          (name == "ME" ? "without" : "with") + " ISR)");
      desc.add<Limits>("WWRNG" + name, {-1.});
    }
    desc.add<Limits>("Q2P", {0.}).setDescription("negative squared mom.transfer range at proton vertex");
    desc.add("THMIN", std::vector<double>(4, 0.))
        .setDescription(
            "minimum theta range for (scattered proton/quark, scattered electron, produced leptons 1, and 2)");
    desc.add("THMAX", std::vector<double>(4, 180.))
        .setDescription(
            "maximum theta range for (scattered proton/quark, scattered electron, produced leptons 1, and 2)");
    desc.add("EMIN", std::vector<double>(4, 0.))
        .setDescription(
            "minimum energy range for (scattered proton/quark, scattered electron, produced leptons 1, and 2)");
    desc.add("EMAX", std::vector<double>(4, kMaxValue))
        .setDescription(
            "maximum energy range for (scattered proton/quark, scattered electron, produced leptons 1, and 2)");
    desc.add("PMIN", std::vector<double>(4, 0.))
        .setDescription(
            "minimum momentum range for (scattered proton/quark, scattered electron, produced leptons 1, and 2)");
    desc.add("PMAX", std::vector<double>(4, kMaxValue))
        .setDescription(
            "maximum momentum range for (scattered proton/quark, scattered electron, produced leptons 1, and 2)");
    desc.add("PTMIN", std::vector<double>(4, 0.))
        .setDescription(
            "minimum transverse momentum range for (scattered proton/quark, scattered electron, produced leptons 1, "
            "and 2)");
    desc.add("PTMAX", std::vector<double>(4, kMaxValue))
        .setDescription(
            "maximum transverse momentum range for (scattered proton/quark, scattered electron, produced leptons 1, "
            "and 2)");
    desc.add<bool>("PTMCTFLG", true);
    desc.add<Limits>("PTMXCT", {0.}).setDescription("two-lepton transverse momentum range");
    desc.add<Limits>("THPTMCT", {0., 180.}).setDescription("two-lepton theta range");
    desc.add("MASSLL", std::vector<Limits>(2, Limits{})).setDescription("two-lepton minimal mass");
    desc.add<Limits>("MASSELL", {0.}).setDescription("two-lepton+scattered electron mass range");
    desc.add("MASSQLL", std::vector<Limits>(2, Limits{5.}))
        .setDescription("two-lepton+scattered quark minimal mass (for DIS process)");
    desc.add<int>("IVISI", -1).setDescription("number of particles required to satisfy (X)V(MIN/MAX) cuts");
    desc.add("THEVMIN", std::vector<double>(4, 0.))
        .setDescription(
            "minimum theta range in lab frame for IVISI (scattered proton/quark, scattered electron, produced leptons "
            "1, and 2)");
    desc.add("THEVMAX", std::vector<double>(4, 180.))
        .setDescription(
            "maximum theta range in lab frame for IVISI (scattered proton/quark, scattered electron, produced leptons "
            "1, and 2)");
    desc.add("EVMIN", std::vector<double>(4, 0.))
        .setDescription(
            "minimum energy range in lab frame for IVISI (scattered proton/quark, scattered electron, produced leptons "
            "1, and 2)");
    desc.add("EVMAX", std::vector<double>(4, kMaxValue))
        .setDescription(
            "maximum energy range in lab frame for IVISI (scattered proton/quark, scattered electron, produced leptons "
            "1, and 2)");
    desc.add("PTVMIN", std::vector<double>(4, 0.))
        .setDescription(
            "minimum transverse momentum range in lab frame for IVISI (scattered proton/quark, scattered electron, "
            "produced leptons 1, and 2)");
    desc.add("PTVMAX", std::vector<double>(4, kMaxValue))
        .setDescription(
            "maximum transverse momentum range in lab frame for IVISI (scattered proton/quark, scattered electron, "
            "produced leptons 1, and 2)");
    desc.add<bool>("LVETO", false);
    desc.add<int>("IVETO", 2);
    desc.add("THEVTMIN", std::vector<double>(4, 0.));
    desc.add("THEVTMAX", std::vector<double>(4, 180.));
    desc.add("EVTMIN", std::vector<double>(4, 0.));
    desc.add("EVTMAX", std::vector<double>(4, kMaxValue));
    desc.add("PTVTMIN", std::vector<double>(4, 0.));
    desc.add("PTVTMAX", std::vector<double>(4, kMaxValue));
    for (const auto& name : {"D"s, "E"s, "F"s}) {
      desc.add<bool>("L3"s + name + "CUT", false);
      desc.add("TH3"s + name + "MIN", std::vector<double>(3, 0.));
      desc.add("TH3"s + name + "MAX", std::vector<double>(3, 180.));
      desc.add("E3"s + name + "MIN", std::vector<double>(3, 0.));
      desc.add("E3"s + name + "MAX", std::vector<double>(3, kMaxValue));
      desc.add("P3"s + name + "MIN", std::vector<double>(3, 0.));
      desc.add("P3"s + name + "MAX", std::vector<double>(3, kMaxValue));
      desc.add("PT3"s + name + "MIN", std::vector<double>(3, 0.));
      desc.add("PT3"s + name + "MAX", std::vector<double>(3, kMaxValue));
      desc.add<Limits>("M12CUT3"s + name, {0.});
      for (const auto& name2 : {"ME"s, "OB"s}) {
        desc.add<Limits>("Q2CT3"s + name + name2, {0.});
        desc.add<Limits>("WCT3"s + name + name2, {0.});
        desc.add<Limits>("XCT3"s + name + name2, {0.});
        desc.add<Limits>("YCT3"s + name + name2, {0.});
      }
    }
    desc.add<bool>("HOTFLG", true);
    desc.add<bool>("L2EVISIA", false);
    desc.add<std::vector<double> >("THE2E", {0., 120., 180.});
    desc.add<std::vector<double> >("EMIN2E", {0., 0.});
    desc.add<bool>("L2EVISIB", false);
    desc.add<std::vector<double> >("THE2EB", {0., 0., 117., 180., 180.});
    desc.add<std::vector<double> >("EMIN2EB", {0., 0., 0., 0.});
    desc.add<bool>("L3EVISI", false);
    desc.add<std::vector<double> >("THE3E", {0., 0., 180., 180.});
    desc.add<double>("EMIN3E", 0.);
    desc.add<bool>("LPTVISI", false);
    desc.add("THEMINPT", std::vector<double>(3, 0.));
    desc.add("THEMAXPT", std::vector<double>(3, 180.));
    desc.add("PTMINPT", std::vector<double>(3, 0.));
    desc.add<bool>("LSCATL", false);
    desc.add<Limits>("ESCATL", {0.});
    desc.add<Limits>("PPRODL", {0.});
    desc.add<Limits>("THESCATL", {0., 180.});
    desc.add<Limits>("PPRODL", {0., 180.});
    desc.add<Limits>("UCUT", {-999.});
    desc.add<Limits>("MHAD", {1.08}).setDescription("hadronic system mass range, in GeV");
    return desc;
  }

  QelaBlockFiller::QelaBlockFiller(const ParametersList& params) : BlockFiller(params) {
    gep_qela_.rk = steer<double>("RK");
    gep_qela_.a_nbd = steer<double>("ANBD");
    gep_qela_.b_nbd = steer<double>("BNBD");
    gep_qela_.c_nbd = steer<double>("CNBD");
    gep_qela_.p_u = steer<double>("PUQ");
    gep_qela_.p_d = steer<double>("PDQ");
    gep_qela_.p_s = steer<double>("PSQ");
    gep_qela_.p_ud_bar = steer<double>("PUD");
    gep_qela_.supp_pi0 = steer<double>("SUPPI0");
    gep_qela_.a_slope = steer<double>("ASLOPE");
    gep_qela_.b_slope = steer<double>("BSLOPE");
    gep_qela_.c_slope = steer<double>("CSLOPE");
    gep_qela_.red_slope_s = steer<double>("REDS");
  }

  ParametersDescription QelaBlockFiller::description() {
    auto desc = ParametersDescription();
    desc.setDescription("Quasi-elastic production parameters");
    desc.add<double>("RK", 0.05);
    desc.add<double>("ANBD", 0.176);
    desc.add<double>("BNBD", 0.43132);
    desc.add<double>("CNBD", 0.86224);
    desc.add<double>("PUQ", 0.5);
    desc.add<double>("PDQ", 0.5);
    desc.add<double>("PSQ", 0.07);
    desc.add<double>("PUD", 0.02);
    desc.add<double>("SUPPI0", 0.4);
    desc.add<double>("ASLOPE", 0.00175);
    desc.add<double>("BSLOPE", 0.00167);
    desc.add<double>("CSLOPE", 0.353);
    desc.add<double>("REDS", 0.6);
    return desc;
  }

  KineBlockFiller::KineBlockFiller(const ParametersList& params) : BlockFiller(params) {
    gep_kine_.isym_34 = steer<int>("ISYM34");
    gep_kine_.ii34 = steer<int>("I34");
    gep_kine_.iscale_xq = steer<int>("SCALEXQ");
    gep_kine_.ireso56 = steer<int>("RESNS56");
    gep_kine_.ireso456 = steer<int>("RESNS456");
    gep_kine_.nregion = steer<int>("NREG");
    gep_kine_.icosp3 = steer<int>("ICOS3");
    gep_kine_.icosp5 = steer<int>("ICOS5");
    gep_kine_.neps_p = steer<int>("NEPSP");
    gep_kine_.thresh56 = steer<double>("THRES56");
    gep_kine_.rnn_mj56 = steer<double>("NCNT");
    gep_kine_.mas_reso456 = steer<double>("MAS456");
    gep_kine_.wid_reso456 = steer<double>("WID456");
    gep_kine_.rscale_mn = steer<double>("RSCALEMN");
    gep_kine_.rnn_456 = steer<double>("RNN456");
  }

  ParametersDescription KineBlockFiller::description() {
    auto desc = ParametersDescription();
    desc.setDescription("Process kinematics parameters");
    desc.add<int>("ISYM34", 1122);
    desc.add<int>("I34", 3);
    desc.add<int>("SCALEXQ", 2);
    desc.add<int>("RESNS56", 1122);
    desc.add<int>("RESNS456", -1);
    desc.add<int>("NREG", 2);
    desc.add<int>("ICOS3", 2);
    desc.add<int>("ICOS5", 2);
    desc.add<int>("NEPSP", -1);
    desc.add<double>("THRES56", 0.9);
    desc.add<double>("NCNT", 2.);
    desc.add<double>("MAS456", 0.);
    desc.add<double>("WID456", 0.);
    desc.add<double>("RSCALEMN", 1.);
    desc.add<double>("RNN456", 2.);
    return desc;
  }

  void initialiseConstants() {
    smcnst_.pi = M_PI;
    smcnst_.pi2 = M_PI * M_PI;
    smcnst_.rad = M_PI / 180.;
    smcnst_.gevpb = constants::GEVM2_TO_PB;
    smcnst_.alpha = constants::ALPHA_EM;
    smcnst_.alphas = constants::ALPHA_QCD;
    smcnst_.alpha0 = constants::ALPHA_EM;
    // initialise the masses
    size_t i = 0;
    for (const auto& part :
         {24,  23,     22,    21,      25,      24,      23,   12,     14,     16,     11,      13,     15,
          2,   4,      6,     1,       2,       3,       2212, 2112,   24,     24,     23,      22,     21,
          443, 100443, 30443, 9000443, 9010443, 9020443, 553,  100553, 200553, 300553, 9000553, 9010553})
      smmass_.mass[i++] = PDG::get().mass(part);
    i = 0;
    for (const auto& part :
         {24, 23,  25,     24,    23,      2,       4,       6,   1,      2,      3,      24,      24,
          23, 443, 100443, 30443, 9000443, 9010443, 9020443, 553, 100553, 200553, 300553, 9000553, 9010553})
      smgmma_.width[i++] = PDG::get().width(part);
  }

  void tuneKinematics() {
    if (gep_kine_.isym_34 == 1122) {
      if (gep_proc_.process != (int)ProtonMode::dis || (gep_ew_.jgra_sel > 0 && gep_ew_.jgra_sel <= 2)) {
        gep_kine_.isym_34 = 0;
        gep_kine_.ii34 = 3;
      } else
        gep_kine_.isym_34 = 1;
    }
    if (gep_kine_.neps_p < 0)
      gep_kine_.neps_p = gep_proc_.process == (int)ProtonMode::dis ? 11 : 12;
    if (gep_proc_.nn_isr < 0)
      gep_proc_.nn_isr = gep_proc_.process == (int)ProtonMode::dis ? 55 : 56;
    if (gep_kine_.ireso56 == 1122) {
      if (gep_ew_.jgra_sel <= 3)
        gep_kine_.ireso56 = -1;
      else if (gep_ew_.jgra_sel == 4) {
        gep_kine_.ireso56 = 3;
        gep_kine_.thresh56 = 0.4;
      } else if (gep_ew_.jgra_sel == 14) {
        gep_kine_.ireso56 = 3;
        gep_kine_.thresh56 = 0.;
      } else
        gep_kine_.ireso56 = -1;
    }
  }

  std::unique_ptr<cepgen::utils::RandomGenerator> random_number_generator;
  std::unique_ptr<cepgen::strfun::Parameterisation> structure_functions;
}  // namespace grape

extern "C" {
double cepgen_random_number_(void) {
  if (!grape::random_number_generator)
    grape::random_number_generator = cepgen::RandomGeneratorFactory::get().build("stl");
  return grape::random_number_generator->uniform();
}

void cepgen_hybrid_structure_functions_(double& q2, double& w, double& w1, double& w2) {
  if (!grape::structure_functions)
    grape::structure_functions = cepgen::StructureFunctionsFactory::get().build("grapeHybrid");
  static const auto mp2 = cepgen::PDG::get().mass(cepgen::PDG::proton);
  const auto xbj = cepgen::utils::xBj(q2, mp2, w * w);
  w1 = grape::structure_functions->W1(xbj, q2);
  w2 = grape::structure_functions->W2(xbj, q2);
}

void taswap_(double& lhs, double& rhs) { std::swap(lhs, rhs); }
bool ibtest_(int& variable, int& bit) { return (variable >> bit) & 0x1; }
}
