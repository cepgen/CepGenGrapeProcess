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

#include "CepGenGrapeProcess/Process.h"

extern "C" {
extern void aep_epeemclr_();
extern void gep_epeefinit_();
extern double fncep_epee_(double*);

extern void aep_epmmmclr_();
extern void gep_epmmfinit_();
extern double fncep_epmm_(double*);

extern void aep_epttmclr_();
extern void gep_epttfinit_();
extern double fncep_eptt_(double*);

extern void aep_exeemclr_();
extern void gep_exeefinit_();
extern double fncep_exee_(double*);

extern void aep_exmmmclr_();
extern void gep_exmmfinit_();
extern double fncep_exmm_(double*);

extern void aep_exttmclr_();
extern void gep_exttfinit_();
extern double fncep_extt_(double*);

extern void aeu_eueemclr_();
extern void geu_eueefinit_();
extern double fnceu_euee_(double*);

extern void aeu_eummmclr_();
extern void geu_eummfinit_();
extern double fnceu_eumm_(double*);

extern void aeu_euttmclr_();
extern void geu_euttfinit_();
extern double fnceu_eutt_(double*);

extern void aeub_eubeemclr_();
extern void geub_eubeefinit_();
extern double fnceub_eubee_(double*);

extern void aeub_eubmmmclr_();
extern void geub_eubmmfinit_();
extern double fnceub_eubmm_(double*);

extern void aeub_eubttmclr_();
extern void geub_eubttfinit_();
extern double fnceub_eubtt_(double*);

extern void aed_edeemclr_();
extern void ged_edeefinit_();
extern double fnced_edee_(double*);

extern void aed_edmmmclr_();
extern void ged_edmmfinit_();
extern double fnced_edmm_(double*);

extern void aed_edttmclr_();
extern void ged_edttfinit_();
extern double fnced_edtt_(double*);

extern void aedb_edbeemclr_();
extern void gedb_edbeefinit_();
extern double fncedb_edbee_(double*);

extern void aedb_edbmmmclr_();
extern void gedb_edbmmfinit_();
extern double fncedb_edbmm_(double*);

extern void aedb_edbttmclr_();
extern void gedb_edbttfinit_();
extern double fncedb_edbtt_(double*);

extern void aes_eseemclr_();
extern void ges_eseefinit_();
extern double fnces_esee_(double*);

extern void aes_esmmmclr_();
extern void ges_esmmfinit_();
extern double fnces_esmm_(double*);

extern void aes_esttmclr_();
extern void ges_esttfinit_();
extern double fnces_estt_(double*);

extern void aesb_esbeemclr_();
extern void gesb_esbeefinit_();
extern double fncesb_esbee_(double*);

extern void aesb_esbmmmclr_();
extern void gesb_esbmmfinit_();
extern double fncesb_esbmm_(double*);

extern void aesb_esbttmclr_();
extern void gesb_esbttfinit_();
extern double fncesb_esbtt_(double*);

extern void aec_eceemclr_();
extern void gec_eceefinit_();
extern double fncec_ecee_(double*);

extern void aec_ecmmmclr_();
extern void gec_ecmmfinit_();
extern double fncec_ecmm_(double*);

extern void aec_ecttmclr_();
extern void gec_ecttfinit_();
extern double fncec_ectt_(double*);

extern void aecb_ecbeemclr_();
extern void gecb_ecbeefinit_();
extern double fncecb_ecbee_(double*);

extern void aecb_ecbmmmclr_();
extern void gecb_ecbmmfinit_();
extern double fncecb_ecbmm_(double*);

extern void aecb_ecbttmclr_();
extern void gecb_ecbttfinit_();
extern double fncecb_ecbtt_(double*);

extern void aeb_ebeemclr_();
extern void geb_ebeefinit_();
extern double fnceb_ebee_(double*);

extern void aeb_ebmmmclr_();
extern void geb_ebmmfinit_();
extern double fnceb_ebmm_(double*);

extern void aeb_ebttmclr_();
extern void geb_ebttfinit_();
extern double fnceb_ebtt_(double*);

extern void aebb_ebbeemclr_();
extern void gebb_ebbeefinit_();
extern double fncebb_ebbee_(double*);

extern void aebb_ebbmmmclr_();
extern void gebb_ebbmmfinit_();
extern double fncebb_ebbmm_(double*);

extern void aebb_ebbttmclr_();
extern void gebb_ebbttfinit_();
extern double fncebb_ebbtt_(double*);

extern void aet_eteemclr_();
extern void get_eteefinit_();
extern double fncet_etee_(double*);

extern void aet_etmmmclr_();
extern void get_etmmfinit_();
extern double fncet_etmm_(double*);

extern void aet_etttmclr_();
extern void get_etttfinit_();
extern double fncet_ettt_(double*);

extern void aetb_etbeemclr_();
extern void getb_etbeefinit_();
extern double fncetb_etbee_(double*);

extern void aetb_etbmmmclr_();
extern void getb_etbmmfinit_();
extern double fncetb_etbmm_(double*);

extern void aetb_etbttmclr_();
extern void getb_etbttfinit_();
extern double fncetb_etbtt_(double*);
}

using namespace grape;

std::function<void()> initialise;
std::function<double(double)> run;
std::function<void()> finalise;

Process Process::fromProcessId(int jproc) {
  Process process;
  switch (jproc) {
    case 1:
      process.name = "ep_epee";
      process.description = "p positron  --> p positron electron positron";
      process.clear = aep_epeemclr_;
      process.initialise = gep_epeefinit_;
      process.run = fncep_epee_;
      return process;
    case 2:
      process.name = "ep_epmm";
      process.description = "p positron  --> p positron muon anti-muon";
      process.clear = aep_epmmmclr_;
      process.initialise = gep_epmmfinit_;
      process.run = fncep_epmm_;
      return process;
    case 3:
      process.name = "ep_eptt";
      process.description = "p positron  --> p positron tau anti-tau";
      process.clear = aep_epttmclr_;
      process.initialise = gep_epttfinit_;
      process.run = fncep_eptt_;
      return process;
    case 4:
      process.name = "ep_exee";
      process.description = "p positron  --> xxx positron electron positron";
      process.clear = aep_exeemclr_;
      process.initialise = gep_exeefinit_;
      process.run = fncep_exee_;
      return process;
    case 5:
      process.name = "ep_exmm";
      process.description = "p positron  --> xxx positron muon anti-muon";
      process.clear = aep_exmmmclr_;
      process.initialise = gep_exmmfinit_;
      process.run = fncep_exmm_;
      return process;
    case 6:
      process.name = "ep_extt";
      process.description = "p positron  --> xxx positron tau anti-tau";
      process.clear = aep_exttmclr_;
      process.initialise = gep_exttfinit_;
      process.run = fncep_extt_;
      return process;
    case 7:
      process.name = "eu_euee";
      process.description = "u positron  --> u positron electron positron";
      process.clear = aeu_eueemclr_;
      process.initialise = geu_eueefinit_;
      process.run = fnceu_euee_;
      return process;
    case 8:
      process.name = "eu_eumm";
      process.description = "u positron  --> u positron muon anti-muon";
      process.clear = aeu_eummmclr_;
      process.initialise = geu_eummfinit_;
      process.run = fnceu_eumm_;
      return process;
    case 9:
      process.name = "eu_eutt";
      process.description = "u positron  --> u positron tau anti-tau";
      process.clear = aeu_euttmclr_;
      process.initialise = geu_euttfinit_;
      process.run = fnceu_eutt_;
      return process;
    case 10:
      process.name = "eub_eubee";
      process.description = "u-bar positron  --> u-bar positron electron positron";
      process.clear = aeub_eubeemclr_;
      process.initialise = geub_eubeefinit_;
      process.run = fnceub_eubee_;
      return process;
    case 11:
      process.name = "eub_eubmm";
      process.description = "u-bar positron  --> u-bar positron muon anti-muon";
      process.clear = aeub_eubmmmclr_;
      process.initialise = geub_eubmmfinit_;
      process.run = fnceub_eubmm_;
      return process;
    case 12:
      process.name = "eub_eubtt";
      process.description = "u-bar positron  --> u-bar positron tau anti-tau";
      process.clear = aeub_eubttmclr_;
      process.initialise = geub_eubttfinit_;
      process.run = fnceub_eubtt_;
      return process;
    case 13:
      process.name = "ed_edee";
      process.description = "d positron  --> d positron electron positron";
      process.clear = aed_edeemclr_;
      process.initialise = ged_edeefinit_;
      process.run = fnced_edee_;
      return process;
    case 14:
      process.name = "ed_edmm";
      process.description = "d positron  --> d positron muon anti-muon";
      process.clear = aed_edmmmclr_;
      process.initialise = ged_edmmfinit_;
      process.run = fnced_edmm_;
      return process;
    case 15:
      process.name = "ed_edtt";
      process.description = "d positron  --> d positron tau anti-tau";
      process.clear = aed_edttmclr_;
      process.initialise = ged_edttfinit_;
      process.run = fnced_edtt_;
      return process;
    case 16:
      process.name = "edb_edbee";
      process.description = "d-bar positron  --> d-bar positron electron positron";
      process.clear = aedb_edbeemclr_;
      process.initialise = gedb_edbeefinit_;
      process.run = fncedb_edbee_;
      return process;
    case 17:
      process.name = "edb_edbmm";
      process.description = "d-bar positron  --> d-bar positron muon anti-muon";
      process.clear = aedb_edbmmmclr_;
      process.initialise = gedb_edbmmfinit_;
      process.run = fncedb_edbmm_;
      return process;
    case 18:
      process.name = "edb_edbtt";
      process.description = "d-bar positron  --> d-bar positron tau anti-tau";
      process.clear = aedb_edbttmclr_;
      process.initialise = gedb_edbttfinit_;
      process.run = fncedb_edbtt_;
      return process;
    case 19:
      process.name = "es_esee";
      process.description = "s positron  --> s positron electron positron";
      process.clear = aes_eseemclr_;
      process.initialise = ges_eseefinit_;
      process.run = fnces_esee_;
      return process;
    case 20:
      process.name = "es_esmm";
      process.description = "s positron  --> s positron muon anti-muon";
      process.clear = aes_esmmmclr_;
      process.initialise = ges_esmmfinit_;
      process.run = fnces_esmm_;
      return process;
    case 21:
      process.name = "es_estt";
      process.description = "s positron  --> s positron tau anti-tau";
      process.clear = aes_esttmclr_;
      process.initialise = ges_esttfinit_;
      process.run = fnces_estt_;
      return process;
    case 22:
      process.name = "esb_esbee";
      process.description = "s-bar positron  --> s-bar positron electron positron";
      process.clear = aesb_esbeemclr_;
      process.initialise = gesb_esbeefinit_;
      process.run = fncesb_esbee_;
      return process;
    case 23:
      process.name = "esb_esbmm";
      process.description = "s-bar positron  --> s-bar positron muon anti-muon";
      process.clear = aesb_esbmmmclr_;
      process.initialise = gesb_esbmmfinit_;
      process.run = fncesb_esbmm_;
      return process;
    case 24:
      process.name = "esb_esbtt";
      process.description = "s-bar positron  --> s-bar positron tau anti-tau";
      process.clear = aesb_esbttmclr_;
      process.initialise = gesb_esbttfinit_;
      process.run = fncesb_esbtt_;
      return process;
    case 25:
      process.name = "ec_ecee";
      process.description = "c positron  --> c positron electron positron";
      process.clear = aec_eceemclr_;
      process.initialise = gec_eceefinit_;
      process.run = fncec_ecee_;
      return process;
    case 26:
      process.name = "ec_ecmm";
      process.description = "c positron  --> c positron muon anti-muon";
      process.clear = aec_ecmmmclr_;
      process.initialise = gec_ecmmfinit_;
      process.run = fncec_ecmm_;
      return process;
    case 27:
      process.name = "ec_ectt";
      process.description = "c positron  --> c positron tau anti-tau";
      process.clear = aec_ecttmclr_;
      process.initialise = gec_ecttfinit_;
      process.run = fncec_ectt_;
      return process;
    case 28:
      process.name = "ecb_ecbee";
      process.description = "c-bar positron  --> c-bar positron electron positron";
      process.clear = aecb_ecbeemclr_;
      process.initialise = gecb_ecbeefinit_;
      process.run = fncecb_ecbee_;
      return process;
    case 29:
      process.name = "ecb_ecbmm";
      process.description = "c-bar positron  --> c-bar positron muon anti-muon";
      process.clear = aecb_ecbmmmclr_;
      process.initialise = gecb_ecbmmfinit_;
      process.run = fncecb_ecbmm_;
      return process;
    case 30:
      process.name = "ecb_ecbtt";
      process.description = "c-bar positron  --> c-bar positron tau anti-tau";
      process.clear = aecb_ecbttmclr_;
      process.initialise = gecb_ecbttfinit_;
      process.run = fncecb_ecbtt_;
      return process;
    case 31:
      process.name = "eb_ebee";
      process.description = "b positron  --> b positron electron positron";
      process.clear = aeb_ebeemclr_;
      process.initialise = geb_ebeefinit_;
      process.run = fnceb_ebee_;
      return process;
    case 32:
      process.name = "eb_ebmm";
      process.description = "b positron  --> b positron muon anti-muon";
      process.clear = aeb_ebmmmclr_;
      process.initialise = geb_ebmmfinit_;
      process.run = fnceb_ebmm_;
      return process;
    case 33:
      process.name = "eb_ebtt";
      process.description = "b positron  --> b positron tau anti-tau";
      process.clear = aeb_ebttmclr_;
      process.initialise = geb_ebttfinit_;
      process.run = fnceb_ebtt_;
      return process;
    case 34:
      process.name = "ebb_ebbee";
      process.description = "b-bar positron  --> b-bar positron electron positron";
      process.clear = aebb_ebbeemclr_;
      process.initialise = gebb_ebbeefinit_;
      process.run = fncebb_ebbee_;
      return process;
    case 35:
      process.name = "ebb_ebbmm";
      process.description = "b-bar positron  --> b-bar positron muon anti-muon";
      process.clear = aebb_ebbmmmclr_;
      process.initialise = gebb_ebbmmfinit_;
      process.run = fncebb_ebbmm_;
      return process;
    case 36:
      process.name = "ebb_ebbtt";
      process.description = "b-bar positron  --> b-bar positron tau anti-tau";
      process.clear = aebb_ebbttmclr_;
      process.initialise = gebb_ebbttfinit_;
      process.run = fncebb_ebbtt_;
      return process;
    case 37:
      process.name = "et_etee";
      process.description = "t positron  --> t positron electron positron";
      process.clear = aet_eteemclr_;
      process.initialise = get_eteefinit_;
      process.run = fncet_etee_;
      return process;
    case 38:
      process.name = "et_etmm";
      process.description = "t positron  --> t positron muon anti-muon";
      process.clear = aet_etmmmclr_;
      process.initialise = get_etmmfinit_;
      process.run = fncet_etmm_;
      return process;
    case 39:
      process.name = "et_ettt";
      process.description = "t positron  --> t positron tau anti-tau";
      process.clear = aet_etttmclr_;
      process.initialise = get_etttfinit_;
      process.run = fncet_ettt_;
      return process;
    case 40:
      process.name = "etb_etbee";
      process.description = "t-bar positron  --> t-bar positron electron positron";
      process.clear = aetb_etbeemclr_;
      process.initialise = getb_etbeefinit_;
      process.run = fncetb_etbee_;
      return process;
    case 41:
      process.name = "etb_etbmm";
      process.description = "t-bar positron  --> t-bar positron muon anti-muon";
      process.clear = aetb_etbmmmclr_;
      process.initialise = getb_etbmmfinit_;
      process.run = fncetb_etbmm_;
      return process;
    case 42:
      process.name = "etb_etbtt";
      process.description = "t-bar positron  --> t-bar positron tau anti-tau";
      process.clear = aetb_etbttmclr_;
      process.initialise = getb_etbttfinit_;
      process.run = fncetb_etbtt_;
      return process;
  }
  throw CG_FATAL("grape:Process") << "Failed to retrieve methods for process #" << jproc << ".";
}
