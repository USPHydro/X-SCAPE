#include <stdio.h>
#include <sys/stat.h>

#include <cstring>

#include "JetScapeLogger.h"
#include "AMPTGenesisWrapper.h"

using namespace std;

// Register the module with the base class
RegisterJetScapeModule<AMPTGenesisWrapper> AMPTGenesisWrapper::reg("CustomModuleAMPTGenesis");

AMPTGenesisWrapper::AMPTGenesisWrapper() {
  preequilibrium_status_ = NOT_STARTED;
  SetId("CustomModuleAMPTGenesis");
}

AMPTGenesisWrapper::~AMPTGenesisWrapper() {
  if (preequilibrium_status_ != NOT_STARTED)
    delete genesis_ptr;
}

void AMPTGenesisWrapper::InitializePreequilibrium(
    PreEquilibriumParameterFile parameter_list) {
  JSINFO << "Initialize AMPTGenesis ...";
  VERBOSE(8);
  using_ampt = true;
  genesis_ptr = new AMPTGenesis();
  std::string input_file = GetXMLElementText(
      {"Preequilibrium", "CustomModuleAMPTGenesis", "freestream_input_file"});
  //is this necessary? if we just force the user to have the 'freestream_input' file in the correct directory
  genesis_ptr->input_folder_path = input_file;
  genesis_ptr->output_file_path = "output_smearer.dat";
  int grid_points_x  = GetXMLElementInt(
      {"Preequilibrium", "CustomModuleAMPTGenesis", "grid_points_x"});
  int grid_points_y  = GetXMLElementInt(
      {"Preequilibrium", "CustomModuleAMPTGenesis", "grid_points_y"});
  int grid_points_eta = GetXMLElementInt(
      {"Preequilibrium", "CustomModuleAMPTGenesis", "grid_points_eta"});
  double grid_size_x = GetXMLElementDouble(
      {"Preequilibrium", "CustomModuleAMPTGenesis", "grid_size_x"});
  double grid_size_y = GetXMLElementDouble(
      {"Preequilibrium", "CustomModuleAMPTGenesis", "grid_size_y"});
  double grid_size_eta = GetXMLElementDouble(
      {"Preequilibrium", "CustomModuleAMPTGenesis", "grid_size_eta"});

  double smearing_constant = GetXMLElementDouble(
      {"Preequilibrium", "CustomModuleAMPTGenesis", "smearing_constant"});
  double smearing_sigma_r =GetXMLElementDouble(
      {"Preequilibrium", "CustomModuleAMPTGenesis", "smearing_sigma_r"});
  double smearing_sigma_eta =GetXMLElementDouble(
      {"Preequilibrium", "CustomModuleAMPTGenesis", "smearing_sigma_eta"});

  double sample_radius_xy = GetXMLElementDouble(
      {"Preequilibrium", "CustomModuleAMPTGenesis", "sample_radius_xy"});
  double sample_radius_eta = GetXMLElementDouble(
      {"Preequilibrium", "CustomModuleAMPTGenesis", "sample_radius_eta"});


  double tau_switch = GetXMLElementDouble(
      {"Preequilibrium", "taus"});
  genesis_ptr->tau0 = tau_switch;

  genesis_ptr->smearing_k = smearing_constant;
  genesis_ptr->nx = grid_points_x;
  genesis_ptr->ny = grid_points_y;
  genesis_ptr->neta = grid_points_eta;
  genesis_ptr->Lx = grid_size_x;
  genesis_ptr->Ly = grid_size_y;
  genesis_ptr->Leta = grid_size_eta;
  genesis_ptr->sigma_r = smearing_sigma_r;
  genesis_ptr->sigma_eta = smearing_sigma_eta;
  genesis_ptr->rxy = sample_radius_xy;
  genesis_ptr->reta = sample_radius_eta;
  xmax_fs  =grid_size_x;
  ymax_fs  =grid_size_y;
  etamax_fs=grid_size_eta;
  dx_fs    =grid_size_x/(grid_points_x - 1);
  dy_fs    =grid_size_y/(grid_points_y - 1);
  deta_fs  =grid_size_eta/(grid_points_eta - 1);
  neta_fs  =grid_points_eta;

}

void AMPTGenesisWrapper::EvolvePreequilibrium() {
  VERBOSE(8);
  JSINFO << "Initialize energy density profile in AMPT-Genesis ...";
  // grab initial energy density from vector from initial state module
  preequilibrium_status_ = INIT;
  if (preequilibrium_status_ == INIT) {
    JSINFO << "running AMPTGenesis ...";
    // evolve the medium via freestreaming
    genesis_ptr->run_genesis();
    preequilibrium_status_ = DONE;
    JSINFO << "finished running AMPTGenesis ...";
  }
  // now prepare to send the resulting hydro variables to the hydro module by coping hydro vectors to Preequilibrium base class members
  genesis_ptr->output_to_vectors(e_, P_, utau_, ux_, uy_, ueta_, pi00_, pi01_,
                                 pi02_, pi03_, pi11_, pi12_, pi13_, pi22_,
                                 pi23_, pi33_, bulk_Pi_,tau_hydro_,rho_b_,q0_,q1_,q2_,q3_);
    std::cout << " Match time ... " << genesis_ptr->tau0<< std::endl;
    std::cout << "nx from wrapper : " << e_.size()/genesis_ptr->neta << std::endl;
}

