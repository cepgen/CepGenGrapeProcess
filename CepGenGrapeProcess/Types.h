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
  /// Feynman diagram selection
  enum struct EwkProductionMode {
    betheHeitlerWithoutEEInteference = 1,  ///< 2-photon B-H (w/o e+e- interference in dielectron case)
    betheHeitler = 2,                      ///< 2-photon B-H (with e+e- interference in dielectron case)
    qed = 3,                               ///< QED diagrams (all the diagrams w/o the Z0 contribution)
    ew = 4,                                ///< All diagrams
    qedCompton = 13,                       ///< QED-Compton type diagrams
    onlyZ0 = 14                            ///< Z0 production diagrams
  };
  enum struct LeptonISRMode { none = 0, sfMethod = 1, qedPartonShowerMethod = 2 };
  /// Process type of the proton vertex
  enum struct ProtonMode { elastic = 1, quasiElastic = 2, dis = 3 };
  /// Merging mode for the DIS process
  enum struct DISMergeMode {
    off = 0,
    u_ubar_d_dbar = 1234,
    u_ubar_d_dbar_s_sbar = 123456,
    u_ubar_d_dbar_s_sbar_c_cbar = 12345678,
    u_ubar_d_dbar_s_sbar_c_cbar_b_bbar = 1234567890,
    u_c = 17,
    ubar_cbar = 28,
    d_s = 35,
    dbar_sbar = 46,
    d_s_b = 359,
    dbar_sbar_bbar = 460
  };
}  // namespace grape

#endif
