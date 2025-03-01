/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include <stdio.h>
#include <sys/stat.h>

#include <cstring>

#include "JetScapeLogger.h"
#include "FreestreamMilneWrapper.h"

using namespace std;

// Register the module with the base class
RegisterJetScapeModule<FreestreamMilneWrapper>
    FreestreamMilneWrapper::reg("FreestreamMilne");

FreestreamMilneWrapper::FreestreamMilneWrapper() {
  preequilibrium_status_ = NOT_STARTED;
  SetId("Freestream-Milne");
}

FreestreamMilneWrapper::~FreestreamMilneWrapper() {
  if (preequilibrium_status_ != NOT_STARTED)
    delete fsmilne_ptr;
}

void FreestreamMilneWrapper::InitializePreequilibrium(
    PreEquilibriumParameterFile parameter_list) {
  JSINFO << "Initialize freestream-milne ...";
  VERBOSE(8);

  std::string input_file = GetXMLElementText(
      {"Preequilibrium", "FreestreamMilne", "freestream_input_file"});
  //is this necessary? if we just force the user to have the 'freestream_input' file in the correct directory

  fsmilne_ptr = new FREESTREAMMILNE();
  struct parameters *params = fsmilne_ptr->configure(input_file.c_str());

  double tau0 = GetXMLElementDouble(
      {"Preequilibrium", "tau0"});
  double taus = GetXMLElementDouble(
      {"Preequilibrium", "taus"});

  params->TAU0 = tau0;
  params->DTAU = taus - tau0;
  tau_hydro_ = taus;
}

void FreestreamMilneWrapper::EvolvePreequilibrium() {
  VERBOSE(8);
  JSINFO << "Initialize energy density profile in freestream-milne ...";
  // grab initial energy density from vector from initial state module
  std::vector<double> entropy_density =
      ini->GetEntropyDensityDistribution(); //note that this is the energy density when read by freestream-milne, not actually the entropy density!
  std::vector<float> entropy_density_float(entropy_density.begin(),
                                           entropy_density.end());
  fsmilne_ptr->initialize_from_vector(entropy_density_float);
  preequilibrium_status_ = INIT;
  if (preequilibrium_status_ == INIT) {
    JSINFO << "running freestream-milne ...";
    // evolve the medium via freestreaming
    fsmilne_ptr->run_freestream_milne();
    preequilibrium_status_ = DONE;
  }
  // now prepare to send the resulting hydro variables to the hydro module by coping hydro vectors to Preequilibrium base class members
  fsmilne_ptr->output_to_vectors(e_, P_, utau_, ux_, uy_, ueta_, pi00_, pi01_,
                                 pi02_, pi03_, pi11_, pi12_, pi13_, pi22_,
                                 pi23_, pi33_, bulk_Pi_, tau_hydro_,rho_b_,q0_,q1_,q2_,q3_);
  cout << "Hydro matching time: " << tau_hydro_ << endl;
}
