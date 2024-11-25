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

#include <CepGen/Modules/StructureFunctionsFactory.h>
#include <CepGen/Physics/Utils.h>
#include <CepGen/StructureFunctions/Parameterisation.h>

#include <cmath>

namespace grape {
  class HybridStructureFunctions : public cepgen::strfun::Parameterisation {
  public:
    explicit HybridStructureFunctions(const cepgen::ParametersList& params)
        : cepgen::strfun::Parameterisation(params),
          resonances_strfun_(
              cepgen::StructureFunctionsFactory::get().build(steer<cepgen::ParametersList>("resonances"))),
          continuum_strfun_(cepgen::StructureFunctionsFactory::get().build(steer<cepgen::ParametersList>("continuum"))),
          w2_cut_(std::pow(steer<double>("wCut"), 2)) {}

    static cepgen::ParametersDescription description() {
      auto desc = cepgen::strfun::Parameterisation::description();
      desc.setDescription("Grape hybrid structure functions");
      desc.add("resonances", cepgen::StructureFunctionsFactory::get().describeParameters("FioreBrasse"));
      desc.add("continuum", cepgen::StructureFunctionsFactory::get().describeParameters("ALLM97"));
      desc.add<double>("wCut", 2.).setDescription("mass value to transition from resonances to continuum regions");
      return desc;
    }

    void eval() override {
      const double w2 = cepgen::utils::mX2(args_.xbj, args_.q2, mp2_);
      if (w2 < w2_cut_) {
        setF2(resonances_strfun_->F2(args_.xbj, args_.q2));
        setFL(resonances_strfun_->FL(args_.xbj, args_.q2));
      } else {
        setF2(continuum_strfun_->F2(args_.xbj, args_.q2));
        setFL(continuum_strfun_->FL(args_.xbj, args_.q2));
      }
    }

  private:
    const std::unique_ptr<cepgen::strfun::Parameterisation> resonances_strfun_, continuum_strfun_;
    const double w2_cut_;
  };
}  // namespace grape
using grape::HybridStructureFunctions;
REGISTER_STRFUN("grapeHybrid", 1001, HybridStructureFunctions);
