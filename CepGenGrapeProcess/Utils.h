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

#ifndef CepGenGrapeProcess_Utils_h
#define CepGenGrapeProcess_Utils_h

#include <CepGen/Core/SteeredObject.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Utils/Limits.h>

namespace grape {
  static constexpr double kMaxValue = 1.e20;
  int convertLeptonId(cepgen::pdgid_t pdgid);
  template <typename T>
  class BlockFiller : public cepgen::SteeredObject<T> {
  public:
    using cepgen::SteeredObject<T>::SteeredObject;

  protected:
    template <typename U, typename V = U, size_t N = 0>
    void fillArray(const std::string& key, std::array<V, N>& var) const {
      auto vec = cepgen::SteeredObject<T>::template steer<std::vector<U> >(key);
      if (vec.empty())
        vec = std::vector<U>(N, cepgen::SteeredObject<T>::template steer<U>(key));
      else if (vec.size() < N)
        for (size_t i = vec.size(); i < N; ++i)
          vec.emplace_back(*vec.rbegin());
      std::copy(vec.begin(), vec.end(), var.begin());
    }
    enum struct LimitSide { min = 0, max = 1 };
    template <typename U, typename V = U, size_t N = 0>
    void fillLimitsFromValue(std::array<V, N>& var,
                             const std::string& key = "",
                             const LimitSide& side = LimitSide::min) const {
      if (!key.empty())
        fillArray<U, V, N>(key, var);
      else
        std::fill(var.begin(), var.end(), (side == LimitSide::min ? -1 : +1) * kMaxValue);
    }
    template <typename U, typename V = U, size_t N = 0>
    void fillLimitsFromPair(std::array<V, N>& min_var,
                            std::array<V, N>& max_var,
                            const std::string& min_key,
                            const std::string& max_key = "") const {
      fillLimitsFromValue<U, V, N>(min_var, min_key, LimitSide::min);
      fillLimitsFromValue<U, V, N>(max_var, max_key, LimitSide::max);
    }
  };
  struct MiscBlockFiller : BlockFiller<MiscBlockFiller> {
    explicit MiscBlockFiller(const cepgen::ParametersList&);
    static cepgen::ParametersDescription description();
  };
  struct ProcBlockFiller : BlockFiller<ProcBlockFiller> {
    explicit ProcBlockFiller(const cepgen::ParametersList&);
    static cepgen::ParametersDescription description();
  };
  /// Initialisation tool for EWK dilepton production parameters
  struct EWBlockFiller : BlockFiller<EWBlockFiller> {
    explicit EWBlockFiller(const cepgen::ParametersList&);
    static cepgen::ParametersDescription description();
  };
  /// Initialisation tool for vector meson production parameters
  struct VMBlocksFiller : BlockFiller<VMBlocksFiller> {
    explicit VMBlocksFiller(const cepgen::ParametersList&);
    static cepgen::ParametersDescription description();
  };
  struct CutsBlockFiller : BlockFiller<CutsBlockFiller> {
    explicit CutsBlockFiller(const cepgen::ParametersList&);
    static cepgen::ParametersDescription description();
  };
  struct QelaBlockFiller : BlockFiller<QelaBlockFiller> {
    explicit QelaBlockFiller(const cepgen::ParametersList&);
    static cepgen::ParametersDescription description();
  };
  struct KineBlockFiller : BlockFiller<KineBlockFiller> {
    explicit KineBlockFiller(const cepgen::ParametersList&);
    static cepgen::ParametersDescription description();
  };
  /// Export limits into common blocks' variables
  template <typename T>
  inline void saveLimits(const cepgen::Limits& lim, T& min, T& max) {
    min = lim.hasMin() ? lim.min() : -kMaxValue;
    max = lim.hasMax() ? lim.max() : +kMaxValue;
  }
  /// Export limits into common blocks' variables
  template <typename T>
  inline void saveLimits(const cepgen::Limits& lim, std::array<T, 2>& minmax) {
    saveLimits(lim, minmax[0], minmax[1]);
  }
  /// Initialise interfacing common blocks' particles properties and constants
  /// \note setmas in Grape
  void initialiseConstants();
  /// Tune the kinematics parameters
  /// \note kinem_auto_tune in Grape
  void tuneKinematics();
}  // namespace grape

#endif
