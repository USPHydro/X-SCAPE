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
// -----------------------------------------
// This is a wrapper for iSpectraSampler (iSS) with the JETSCAPE framework
// -----------------------------------------

#ifndef ISPECTRASAMPLERWRAPPER_H
#define ISPECTRASAMPLERWRAPPER_H

#include <memory>

#include "SoftParticlization.h"
#include "iSS.h"

using namespace Jetscape;

class iSpectraSamplerWrapper : public SoftParticlization {
private:
  tinyxml2::XMLElement *iSS_xml_;

  int statusCode_;
  std::unique_ptr<iSS> iSpectraSampler_ptr_;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<iSpectraSamplerWrapper> reg;

  //Cache the XML input
  std::string input_file;
  std::string table_path;
  std::string particle_table_path;
  std::string working_path;
  int hydro_mode;
  int number_of_repeated_sampling;
  int afterburner_type;
  int echoLevel;
  int output_samples_into_files;
  int use_OSCAR_format;
  int use_gzip_format;
  int use_binary_format;
  int store_samples_in_memory;
  int turn_on_shear;
  int turn_on_bulk;
  int turn_on_rhob;
  int turn_on_diff;
  int include_deltaf_shear;
  int include_deltaf_bulk;
  int bulk_deltaf_kind;
  int include_deltaf_diffusion;
  int restrict_deltaf;
  int MC_sampling;
  double deltaf_max_ratio;
  int f0_is_not_small;
  int calculate_vn;
  int RegVisYield;
  int include_spectators;
  int sample_upto_desired_particle_number;
  void set_iSpectraSampler_ptr_parameters();

public:
  iSpectraSamplerWrapper();
  ~iSpectraSamplerWrapper();

  void CalculateTime();
  void ExecTime();

  void InitTask();
  void ExecuteTask();
  void ClearTask();
  void ClearHadronList();
  void WriteTask(weak_ptr<JetScapeWriter> w);

  int getSurfCellVector();
  void PassHadronListToJetscape();
  void PassHadronListToJetscapeSameEvent();
};

#endif // ISPECTRASAMPLERWRAPPER_H
