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

#ifndef TRENTOINITIAL_H
#define TRENTOINITIAL_H

#include <tuple>
#include <memory>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

#include "H5Cpp.h"

#include "fwd_decl.h"
#include "JetScapeModuleBase.h"
#include "collider.h"
#include "InitialState.h"
#include "JetScapeLogger.h"

using OptDesc = po::options_description;
using VarMap = po::variables_map;
using namespace H5;
using namespace trento;

namespace Jetscape {

typedef struct {
  double impact_parameter;
  double num_participant;
  double num_binary_collisions;
  double total_entropy;
  std::map<int, double> ecc; // order, eccentricity
  std::map<int, double> psi; // order, participant_plane
  double xmid, ymid;
} EventInfo;


////////////////////////// Trento Initial Condition Results ////////////////////
#define NMAX 6
#define MMAX 6
typedef struct{
    double b;                       // impact parameter [fm]
    int Npart;                      // number of participants
    int Ncoll;                      // number of collisions
    double E;                       // Total energy [GeV]
    double S;                       // Total entropy (as compued via ideal gas EOS)
    double S_hotQCD;                // Total entropy (as compued via hotQCD EOS)
    double E_entropy;               // Total energy (assumes trento output is in entropy density and use ideal EOS to get energy density)
    double re_ecc_p;                // Real part of mom. anisotropy computed with ideal gas EOS
    double im_ecc_p;                // Imag part of mom. anisotropy computed with ideal gas EOS
    double re_ecc_p_hotQCD;         // Real part of mom. anisotropy computed with hotQCD EOS
    double im_ecc_p_hotQCD;         // Imag part of mom. anisotropy computed with hotQCD EOS
    double re_ecc_p_prime;          // Real part of mom. anisotropy computed with ideal gas EOS - ideal hydro contribution only
    double im_ecc_p_prime;          // Imag part of mom. anisotropy computed with ideal gas EOS - ideal hydro contribution only
    double re_ecc_p_prime_hotQCD;   // Real part of mom. anisotropy computed with hotQCD EOS - ideal hydro contribution only
    double im_ecc_p_prime_hotQCD;   // Imag part of mom. anisotropy computed with hotQCD EOS - ideal hydro contribution only
    double R[NMAX];                 // <R^n>^{1/n} (n=1,..,NMAX)
    double eps[NMAX][MMAX+1];       // <eps_{n,m}> (n=1,..,NMAX, m=0,..,MMAX)
    double psi[NMAX][MMAX+1];       // <psi_{n,m}> (n=1,..,NMAX, m=0,..,MMAX)
} trento_event_info;


/**The output data format (from http://qcd.phy.duke.edu/trento/usage.html#output-options):
 * The grid will always be a square N × N array, with N = ceil(2*max/step).
 * So e.g. the default settings (max = 10 fm, step = 0.2 fm) imply a 100 × 100 grid.
 * The ceiling function ensures that the number of steps is always rounded up,
 * so e.g. given max = 10 fm and step 0.3 fm, the grid will be 67 × 67.
 * In this case, the actual grid max will be marginally increased (max = nsteps*step/2).
**/

////////////////////////// Trento Initial Condition Wrapper //////////////////////
class TrentoInitial : public InitialState {
public:
  // Initialize from XML configuration
  TrentoInitial();
  ~TrentoInitial();

  void InitTask();
  void ExecuteTask();
  void ClearTask();
  void WriteToHDF5();

  //Auxiliary functions to output to HDF5
  trento_event_info get_ecc(const std::vector<double> &eps,
                                    const  std::vector<double> &u0,
                                    const  std::vector<double> &ux,
                                    const  std::vector<double> &uy,
                                    const  std::vector<double> &Pi,
                                    const  std::vector<double> &pixx,
                                    const  std::vector<double> &piyy,
                                    const  std::vector<double> &pitautau,
                                    const  std::vector<double> &pixy,
                                    const  double &tau,
                                    const  std::array<int,2> &size, const std::array<double,2> &step);
  void output_hdf5(std::string filename,
                   std::vector<trento_event_info>& evt_vec,
                   std::vector<double>& entropy_density_distribution);

  struct RangeFailure : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };
  EventInfo info_;

private:
  std::shared_ptr<trento::Collider> TrentoGen_;
  std::pair<double, double> GenCenTab(std::string proj, std::string targ,
                                      VarMap var_map, int cL, int cH);
  /// The output instance.
  // Output output_;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<TrentoInitial> reg;
  int event_counter;
};

} // end namespace Jetscape

#endif // TRENTOINITIAL_H
