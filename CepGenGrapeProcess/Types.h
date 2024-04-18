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

#ifndef CepGenGrapeProcess_Types_h
#define CepGenGrapeProcess_Types_h

namespace grape {
  enum struct EwkProductionMode { directBetheHeitler = 1, betheHeitler = 2, qed = 3, ew = 4, onlyCO = 13, onlyZ0 = 14 };
  enum struct LeptonISRMode { none = 0, sfMethod = 1, qedPartonShowerMethod = 2 };
  enum struct ProtonMode { elastic = 1, quasiElastic = 2, dis = 3 };
}  // namespace grape

#endif
