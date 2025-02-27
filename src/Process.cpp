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

extern void aep_epmmmclr_();
extern void gep_epmmfinit_();

extern void aep_epttmclr_();
extern void gep_epttfinit_();

extern void aep_exeemclr_();
extern void gep_exeefinit_();

extern void aep_exmmmclr_();
extern void gep_exmmfinit_();

extern void aep_exttmclr_();
extern void gep_exttfinit_();

extern void aeu_eueemclr_();
extern void geu_eueefinit_();

extern void aeu_eummmclr_();
extern void geu_eummfinit_();

extern void aeu_euttmclr_();
extern void geu_euttfinit_();

extern void aeub_eubeemclr_();
extern void geub_eubeefinit_();

extern void aeub_eubmmmclr_();
extern void geub_eubmmfinit_();

extern void aeub_eubttmclr_();
extern void geub_eubttfinit_();

extern void aed_edeemclr_();
extern void ged_edeefinit_();

extern void aed_edmmmclr_();
extern void ged_edmmfinit_();

extern void aed_edttmclr_();
extern void ged_edttfinit_();

extern void aedb_edbeemclr_();
extern void gedb_edbeefinit_();

extern void aedb_edbmmmclr_();
extern void gedb_edbmmfinit_();

extern void aedb_edbttmclr_();
extern void gedb_edbttfinit_();

extern void aes_eseemclr_();
extern void ges_eseefinit_();

extern void aes_esmmmclr_();
extern void ges_esmmfinit_();

extern void aes_esttmclr_();
extern void ges_esttfinit_();

extern void aesb_esbeemclr_();
extern void gesb_esbeefinit_();

extern void aesb_esbmmmclr_();
extern void gesb_esbmmfinit_();

extern void aesb_esbttmclr_();
extern void gesb_esbttfinit_();

extern void aec_eceemclr_();
extern void gec_eceefinit_();

extern void aec_ecmmmclr_();
extern void gec_ecmmfinit_();

extern void aec_ecttmclr_();
extern void gec_ecttfinit_();

extern void aecb_ecbeemclr_();
extern void gecb_ecbeefinit_();

extern void aecb_ecbmmmclr_();
extern void gecb_ecbmmfinit_();

extern void aecb_ecbttmclr_();
extern void gecb_ecbttfinit_();

extern void aeb_ebeemclr_();
extern void geb_ebeefinit_();

extern void aeb_ebmmmclr_();
extern void geb_ebmmfinit_();

extern void aeb_ebttmclr_();
extern void geb_ebttfinit_();

extern void aebb_ebbeemclr_();
extern void gebb_ebbeefinit_();

extern void aebb_ebbmmmclr_();
extern void gebb_ebbmmfinit_();

extern void aebb_ebbttmclr_();
extern void gebb_ebbttfinit_();

extern void aet_eteemclr_();
extern void get_eteefinit_();

extern void aet_etmmmclr_();
extern void get_etmmfinit_();

extern void aet_etttmclr_();
extern void get_etttfinit_();

extern void aetb_etbeemclr_();
extern void getb_etbeefinit_();

extern void aetb_etbmmmclr_();
extern void getb_etbmmfinit_();

extern void aetb_etbttmclr_();
extern void getb_etbttfinit_();
}

using namespace grape;

std::function<void()> initialise;
std::function<double(double)> run;
std::function<void()> finalise;

Process Process::fromProcessId(int jproc) {
  Process process;
  switch (jproc) {
    case 1:
      process.clear = aep_epeemclr_;
      process.initialise = gep_epeefinit_;
      return process;
    case 2:
      process.clear = aep_epmmmclr_;
      process.initialise = gep_epmmfinit_;
      return process;
    case 3:
      process.clear = aep_epttmclr_;
      process.initialise = gep_epttfinit_;
      return process;
    case 4:
      process.clear = aep_exeemclr_;
      process.initialise = gep_exeefinit_;
      return process;
    case 5:
      process.clear = aep_exmmmclr_;
      process.initialise = gep_exmmfinit_;
      return process;
    case 6:
      process.clear = aep_exttmclr_;
      process.initialise = gep_exttfinit_;
      return process;
    case 7:
      process.clear = aeu_eueemclr_;
      process.initialise = geu_eueefinit_;
      return process;
    case 8:
      process.clear = aeu_eummmclr_;
      process.initialise = geu_eummfinit_;
      return process;
    case 9:
      process.clear = aeu_euttmclr_;
      process.initialise = geu_euttfinit_;
      return process;
    case 10:
      process.clear = aeub_eubeemclr_;
      process.initialise = geub_eubeefinit_;
      return process;
    case 11:
      process.clear = aeub_eubmmmclr_;
      process.initialise = geub_eubmmfinit_;
      return process;
    case 12:
      process.clear = aeub_eubttmclr_;
      process.initialise = geub_eubttfinit_;
      return process;
    case 13:
      process.clear = aed_edeemclr_;
      process.initialise = ged_edeefinit_;
      return process;
    case 14:
      process.clear = aed_edmmmclr_;
      process.initialise = ged_edmmfinit_;
      return process;
    case 15:
      process.clear = aed_edttmclr_;
      process.initialise = ged_edttfinit_;
      return process;
    case 16:
      process.clear = aedb_edbeemclr_;
      process.initialise = gedb_edbeefinit_;
      return process;
    case 17:
      process.clear = aedb_edbmmmclr_;
      process.initialise = gedb_edbmmfinit_;
      return process;
    case 18:
      process.clear = aedb_edbttmclr_;
      process.initialise = gedb_edbttfinit_;
      return process;
    case 19:
      process.clear = aes_eseemclr_;
      process.initialise = ges_eseefinit_;
      return process;
    case 20:
      process.clear = aes_esmmmclr_;
      process.initialise = ges_esmmfinit_;
      return process;
    case 21:
      process.clear = aes_esttmclr_;
      process.initialise = ges_esttfinit_;
      return process;
    case 22:
      process.clear = aesb_esbeemclr_;
      process.initialise = gesb_esbeefinit_;
      return process;
    case 23:
      process.clear = aesb_esbmmmclr_;
      process.initialise = gesb_esbmmfinit_;
      return process;
    case 24:
      process.clear = aesb_esbttmclr_;
      process.initialise = gesb_esbttfinit_;
      return process;
    case 25:
      process.clear = aec_eceemclr_;
      process.initialise = gec_eceefinit_;
      return process;
    case 26:
      process.clear = aec_ecmmmclr_;
      process.initialise = gec_ecmmfinit_;
      return process;
    case 27:
      process.clear = aec_ecttmclr_;
      process.initialise = gec_ecttfinit_;
      return process;
    case 28:
      process.clear = aecb_ecbeemclr_;
      process.initialise = gecb_ecbeefinit_;
      return process;
    case 29:
      process.clear = aecb_ecbmmmclr_;
      process.initialise = gecb_ecbmmfinit_;
      return process;
    case 30:
      process.clear = aecb_ecbttmclr_;
      process.initialise = gecb_ecbttfinit_;
      return process;
    case 31:
      process.clear = aeb_ebeemclr_;
      process.initialise = geb_ebeefinit_;
      return process;
    case 32:
      process.clear = aeb_ebmmmclr_;
      process.initialise = geb_ebmmfinit_;
      return process;
    case 33:
      process.clear = aeb_ebttmclr_;
      process.initialise = geb_ebttfinit_;
      return process;
    case 34:
      process.clear = aebb_ebbeemclr_;
      process.initialise = gebb_ebbeefinit_;
      return process;
    case 35:
      process.clear = aebb_ebbmmmclr_;
      process.initialise = gebb_ebbmmfinit_;
      return process;
    case 36:
      process.clear = aebb_ebbttmclr_;
      process.initialise = gebb_ebbttfinit_;
      return process;
    case 37:
      process.clear = aet_eteemclr_;
      process.initialise = get_eteefinit_;
      return process;
    case 38:
      process.clear = aet_etmmmclr_;
      process.initialise = get_etmmfinit_;
      return process;
    case 39:
      process.clear = aet_etttmclr_;
      process.initialise = get_etttfinit_;
      return process;
    case 40:
      process.clear = aetb_etbeemclr_;
      process.initialise = getb_etbeefinit_;
      return process;
    case 41:
      process.clear = aetb_etbmmmclr_;
      process.initialise = getb_etbmmfinit_;
      return process;
    case 42:
      process.clear = aetb_etbttmclr_;
      process.initialise = getb_etbttfinit_;
      return process;
  }
  throw CG_FATAL("grape:Process") << "Failed to retrieve methods for process #" << jproc << ".";
}
