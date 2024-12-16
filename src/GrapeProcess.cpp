/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
 *                2003  Tetsuo Abe <tabe@post.kek.jp>
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
#include <CepGen/Event/Event.h>
#include <CepGen/Modules/ProcessFactory.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Process/Process.h>
#include <CepGen/Utils/StreamCollector.h>
#include <CepGen/Utils/String.h>

#include "CepGenGrapeProcess/Interface.h"
#include "CepGenGrapeProcess/Types.h"
#include "CepGenGrapeProcess/Utils.h"

using namespace cepgen;

class GrapeProcess final : public cepgen::proc::Process {
public:
  explicit GrapeProcess(const ParametersList& params)
      : proc::Process(params),
        debug_(steer<bool>("debug")),
        pair_(steer<ParticleProperties>("LPAIR")),
        electron_pol_(steerAs<double, float>("electronPolarisation")),
        electron_pol_polar_(steerAs<double, float>("electronPolarisationPolarAngle")),
        electron_pol_azimuthal_(steerAs<double, float>("electronPolarisationAzimuthalAngle")),
        vec_electron_pol_(Momentum::fromPThetaPhiE(electron_pol_, electron_pol_polar_, electron_pol_azimuthal_)),
        me_(PDG::get().mass(PDG::electron)),
        me2_(me_ * me_) {
    if (!Limits{0., 1.}.contains(electron_pol_))
      throw CG_FATAL("GrapeProcess") << "Electron polarisation should be between 0 and 1. User-provided value: "
                                     << electron_pol_ << ".";
    std::fill(m_x_.begin(), m_x_.end(), 0.);
    {                // call the Fortran initialisation subroutine
      int mode = 2;  // event generation mode used to retrieve the functional form
      start_grape_(mode);
    }

    gep_beam_.ebeam_pol = {electron_pol_, electron_pol_polar_, electron_pol_azimuthal_};

    grape::MiscBlockFiller{params_};
    grape::initialiseConstants();  // setmas in Grape
    amparm_();
    {
      std::string logging;
      auto sc = utils::StreamCollector(logging);
      prmass_();
      CG_DEBUG("GrapeProcess") << "From prmass: " << logging;
    }

    grape::ProcBlockFiller{params_};
    grape::EWBlockFiller{params_};
    grape::VMBlocksFiller{params_};
    grape::CutsBlockFiller{params_};
    grape::QelaBlockFiller{params_};
    grape::KineBlockFiller{params_};
  }

  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new GrapeProcess(*this)); }

  void addEventContent() override {
    proc::Process::setEventContent({{Particle::Role::IncomingBeam1, {PDG::electron}},
                                    {Particle::Role::IncomingBeam2, {PDG::proton}},
                                    {Particle::Role::Parton1, {PDG::photon}},
                                    {Particle::Role::Parton2, {PDG::photon}},
                                    {Particle::Role::OutgoingBeam1, {PDG::electron}},
                                    {Particle::Role::OutgoingBeam2, {PDG::proton}},
                                    {Particle::Role::CentralSystem, {+(spdgid_t)pair_.pdgid, -(spdgid_t)pair_.pdgid}}});
  }
  void prepareKinematics() override {
    grape::tuneKinematics();

    get_proc_();
    get_graph_flag_();

    gfinit_();
    initialiseKinematics();

    for (size_t i = 0; i < 7; ++i)
      defineVariable(m_x_[i], Mapping::linear, {0., 1.}, utils::format("x_%zu", i));
    if (gep_proc_.process != (int)grape::ProtonMode::elastic)
      defineVariable(m_x_[ndim()], Mapping::linear, {0., 1.}, "x_inel");
    if (gep_proc_.isr_flag == (int)grape::LeptonISRMode::sfMethod)
      defineVariable(m_x_[ndim()], Mapping::linear, {0., 1.}, "x_isr");
    bparm1_.ndim = ndim();

    kmreg_.mxreg = 2;  // two regions defined for cuts

    CG_INFO("GrapeProcess:prepareKinematics") << event().oneWithRole(Particle::Role::IncomingBeam1).pdgId() << " beam.";
    if (gep_beam_.kf_lbeam == 11) {  // do a few swaps for electron beams
      for (const auto& id : {1, 3, 4, 5})
        kminfo_.kfcode[id] *= -1;
      for (size_t i = 0; i < 3; ++i) {
        std::swap(smcplc_.czll[i][0], smcplc_.czll[i][1]);
        for (size_t j = 0; j < 2; ++j) {
          smcplc_.call[i][j] *= -1;
          smcplc_.czll[i][j] *= -1;
        }
      }
    }
    CG_INFO("GrapeProcess") << "Grape process initialised for dimension-" << ndim() << " phase space volume.\n\t"
                            << "Internal process number ('jproc'): " << amjprc_.jproc << ".";

    if (debug_)
      write_cards_();  // dump the steering card-type information into the output stream
  }
  double computeWeight() override { return func_(m_x_.data()); }
  void fillKinematics() override {
    pA() = Momentum(sp4vec_.vec[1].data());
    pB() = Momentum(sp4vec_.vec[0].data());
    pX() = Momentum(sp4vec_.vec[3].data());
    pY() = Momentum(sp4vec_.vec[2].data());
    pc(0) = Momentum(sp4vec_.vec[4].data());
    pc(1) = Momentum(sp4vec_.vec[5].data());
  }

  static ParametersDescription description() {
    auto desc = proc::Process::description();
    desc.setDescription("γγ → l⁺l¯ (Grape)");
    desc.add<bool>("debug", true).setDescription("debugging mode (Fortran blocks printout)");

    desc.add<double>("electronPolarisation", 0.)
        .setDescription("degree of polarisation of the electron beam (between 0 and 1)");
    desc.add<double>("electronPolarisationPolarAngle", 0.)
        .setDescription("polar angle of the electron beam polarization vector");
    desc.add<double>("electronPolarisationAzimuthalAngle", 0.)
        .setDescription("azimuthal angle of the electron beam polarization vector");

    desc += grape::MiscBlockFiller::description();
    desc += grape::ProcBlockFiller::description();
    desc += grape::EWBlockFiller::description();
    desc += grape::VMBlocksFiller::description();
    desc += grape::CutsBlockFiller::description();
    desc += grape::QelaBlockFiller::description();
    desc += grape::KineBlockFiller::description();
    return desc;
  }

private:
  /// Initialise interfacing common blocks' kinematics definitions
  /// \note kinit_4f in Grape
  inline void initialiseKinematics() {
    gep_beam_.p_e_beam = pA().p() * 1.e3;  // EBEAM expressed in MeV/c
    gep_beam_.p_p_beam = pB().p() * 1.e3;  // PBEAM expressed in MeV/c
    gep_beam_.kf_lbeam =
        event().oneWithRole(Particle::Role::IncomingBeam1).integerPdgId();  // KFLBEAM, KF code of the lepton beam

    // 1=proton, 2=electron
    gep_lab_.p1_lab = pB().p();
    gep_lab_.p2_lab = pA().p();
    gep_lab_.e1_lab = pB().energy();
    gep_lab_.e2_lab = pA().energy();
    for (size_t i = 0; i < 3; ++i) {
      kinem1_.s[i] = mp2_ + me2_ + 2. * (gep_lab_.e1_lab * gep_lab_.e2_lab + gep_lab_.p1_lab * gep_lab_.p2_lab);
      kinem1_.w[i] = std::sqrt(std::max(kinem1_.s[i], 0.));
      gep_lab_.pcms_lab[i] = gep_lab_.p1_lab - gep_lab_.p2_lab;
      gep_lab_.ecms_lab[i] = gep_lab_.e1_lab + gep_lab_.e2_lab;
      gep_lab_.gammacms_lab[i] = gep_lab_.ecms_lab[i] / kinem1_.w[i];
      gep_lab_.betgamcms_lab[i] = gep_lab_.pcms_lab[i] / kinem1_.w[i];
    }
    CG_INFO("GrapeProcess:prepareKinematics") << "\n"
                                              << "********** Information (in Lab. frame) **********\n"
                                              << "                              (in unit of GeV)\n"
                                              << "  P of electrons    = " << gep_lab_.p2_lab << "\n"
                                              << "  P of protons      = " << gep_lab_.p1_lab << "\n"
                                              << "  Mass of electron  = " << me_ << "\n"
                                              << "  Mass of proton    = " << mp_ << "\n"
                                              << "  sqrt(S)           = " << kinem1_.w[0] << "\n"
                                              << "  P of CMS          = " << gep_lab_.pcms_lab[0] << "\n"
                                              << "  E of CMS          = " << gep_lab_.ecms_lab[0] << "\n"
                                              << "  gamma of CMS      = " << gep_lab_.gammacms_lab[0] << "\n"
                                              << "  beta*gamma of CMS = " << gep_lab_.betgamcms_lab[0] << "\n"
                                              << "  masses            = " << kmmass_.amass1 << "\n"
                                              << "***********************************************";
    {  // check on centre-of-mass energy (should be sufficient to create the final state particles)
      float totmas = 0;
      for (size_t i = 2; i < 6; ++i)
        totmas += kmmass_.amass1[i];
      if (kinem1_.w[0] < totmas)
        throw CG_ERROR("GrapeProcess:.initialiseKinematics")
            << "CM energy = " << kinem1_.w[0] << " GeV < sum of final particles mass = " << totmas << " GeV.";
      grape::saveLimits(Limits{1.08}, gep_proc_.wmin, gep_proc_.wmax);
      gep_proc_.wmax = std::min(gep_proc_.wmax, float(kinem1_.w[0]) - totmas);
    }

    grape::saveLimits(Limits{0., 1.}, w50513_.xmin, w50513_.xmax);
    grape::saveLimits(Limits{0., 1.e4}, w50513_.q2min, w50513_.q2max);

    // angular cuts in CMS
    cut001_4f_.opncut = 0.;
    double angcut = 0.;
    grape::saveLimits(Limits{-std::cos(angcut), std::cos(angcut)}, cut001_4f_.coscut[0][0], cut001_4f_.coscut[0][1]);
    grape::saveLimits(Limits{-std::cos(angcut), std::cos(angcut)}, cut001_4f_.coscut[1][0], cut001_4f_.coscut[1][1]);
    //angcut=10.*smcnst_.rad;
    //angcut=0.;
    grape::saveLimits(Limits{-std::cos(angcut), std::cos(angcut)}, cut001_4f_.coscut[2][0], cut001_4f_.coscut[2][1]);
    grape::saveLimits(Limits{-std::cos(angcut), std::cos(angcut)}, cut001_4f_.coscut[3][0], cut001_4f_.coscut[3][1]);

    // energy  cuts in CMS
    const auto max_beam_energy = std::max(gep_lab_.e1_lab, gep_lab_.e2_lab);
    for (size_t i = 0; i < 4; ++i)
      grape::saveLimits(
          Limits{kmmass_.amass1[i + 2], max_beam_energy}, cut001_4f_.engyct[i][0], cut001_4f_.engyct[i][1]);

    grape::saveLimits(  // cut on invariant mass (particles 3-4)
        Limits{kmmass_.amass1[2] + kmmass_.amass1[3], kinem1_.w[0] - kmmass_.amass1[4] - kmmass_.amass1[5]},
        cut001_4f_.amasct[0][0],
        cut001_4f_.amasct[0][1]);

    const auto lim_56 = Limits{
        std::max(gep_cut_.mass56_cut[0][0], float(kmmass_.amass1[4] + kmmass_.amass1[5])),
        std::min(gep_proc_.lpair == 1 ? gep_cut_.mass56_cut[1][1] : gep_cut_.mass56_cut[0][1], float(kinem1_.w[0]))};
    if (!lim_56.valid())  // checking upper and lower values
      throw CG_FATAL("GrapeProcess") << "Invalid mass cut for particle 5-6: " << lim_56 << ".";
    grape::saveLimits(  // particle 5-6
        lim_56,
        cut001_4f_.amasct[1][0],
        cut001_4f_.amasct[1][1]);
    grape::saveLimits(  // particle 3-5
        Limits{kmmass_.amass1[2] + kmmass_.amass1[4], kinem1_.w[0] - kmmass_.amass1[3] - kmmass_.amass1[5]},
        cut001_4f_.amasct[2][0],
        cut001_4f_.amasct[2][1]);
    grape::saveLimits(  // particle 4-6
        Limits{kmmass_.amass1[3] + kmmass_.amass1[5], kinem1_.w[0] - kmmass_.amass1[2] - kmmass_.amass1[4]},
        cut001_4f_.amasct[3][0],
        cut001_4f_.amasct[3][1]);
    grape::saveLimits(  // particle 3-6
        Limits{kmmass_.amass1[2] + kmmass_.amass1[5], kinem1_.w[0] - kmmass_.amass1[3] - kmmass_.amass1[4]},
        cut001_4f_.amasct[4][0],
        cut001_4f_.amasct[4][1]);
    grape::saveLimits(  // particle 4-5
        Limits{kmmass_.amass1[3] + kmmass_.amass1[4], kinem1_.w[0] - kmmass_.amass1[2] - kmmass_.amass1[5]},
        cut001_4f_.amasct[5][0],
        cut001_4f_.amasct[5][1]);

    cmn_isrpol_.isrpol = static_cast<grape::LeptonISRMode>(gep_proc_.isr_flag) == grape::LeptonISRMode::sfMethod;
    k2nit_();
    {
      int jump = 0;  //FIXME
      double mp2 = mp2_, me2 = me2_;
      kinem1_.fact = constants::GEVM2_TO_PB * flux_factor_(kinem1_.s[0], mp2, me2, jump);
    }
    gep_beam_.lpol_ebeam = electron_pol_ != 0.;
    gep_beam_.ipol_ebeam = gep_beam_.lpol_ebeam && electron_pol_polar_ != 0.;
  }
  const bool debug_;
  const ParticleProperties pair_;
  const float electron_pol_, electron_pol_polar_, electron_pol_azimuthal_;

  const Momentum vec_electron_pol_;  ///< electron polarisation vector
  const double me_;                  ///< electron mass
  const double me2_;                 ///< electron squared mass;
  std::array<double, 10> m_x_;       ///< mapped variables
};
// register process
REGISTER_PROCESS("grape", GrapeProcess);
