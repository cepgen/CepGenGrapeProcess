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

#ifndef CepGenGrapeProcess_Process_h
#define CepGenGrapeProcess_Process_h

#include <functional>

namespace grape {
  class Process {
  public:
    static Process fromProcessId(int);

    std::function<void()> clear;
    std::function<void()> initialise;
    std::function<double(double)> run;
    std::function<void()> finalise;
  };
}  // namespace grape

#endif

extern "C" {
void gep_epeefinit_();
void gep_epmmfinit_();
void gep_epttfinit_();
void gep_eXeefinit_();
void gep_eXmmfinit_();
void gep_eXttfinit_();
void geu_eueefinit_();
void geu_eummfinit_();
void geu_euttfinit_();
void geub_eubeefinit_();
void geub_eubmmfinit_();
void geub_eubttfinit_();
void ged_edeefinit_();
void ged_edmmfinit_();
void ged_edttfinit_();
void gedb_edbeefinit_();
void gedb_edbmmfinit_();
void gedb_edbttfinit_();
void ges_eseefinit_();
void ges_esmmfinit_();
void ges_esttfinit_();
void gesb_esbeefinit_();
void gesb_esbmmfinit_();
void gesb_esbttfinit_();
void gec_eceefinit_();
void gec_ecmmfinit_();
void gec_ecttfinit_();
void gecb_ecbeefinit_();
void gecb_ecbmmfinit_();
void gecb_ecbttfinit_();
void geb_ebeefinit_();
void geb_ebmmfinit_();
void geb_ebttfinit_();
void gebb_ebbeefinit_();
void gebb_ebbmmfinit_();
void gebb_ebbttfinit_();
void get_eteefinit_();
void get_etmmfinit_();
void get_etttfinit_();
void getb_etbeefinit_();
void getb_etbmmfinit_();
void getb_etbttfinit_();
}
