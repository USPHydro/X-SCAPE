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

#ifndef PREEQUILDYNAMICS_H
#define PREEQUILDYNAMICS_H

#include <vector>
#include "InitialState.h"
#include "JetScapeModuleBase.h"
#include "RealType.h"

#include "H5Cpp.h"

namespace Jetscape {
// Flags for preequilibrium dynamics status.
enum PreequilibriumStatus { NOT_STARTED, INIT, DONE, ERR };

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
} preequilibrium_event_info;

class PreEquilibriumParameterFile {
public:
  // preequilibrium dynamics parameters file name.
  char *preequilibrium_input_filename;
};

// Interface for the Preequilibrium Dynamics of the medium
class PreequilibriumDynamics : public JetScapeModuleBase {
private:
  PreEquilibriumParameterFile parameter_list_;
  // record preequilibrium start and end proper time [fm/c]
  real preequilibrium_tau_0_, preequilibrium_tau_max_;
  int event_counter;
  std::string output_filename_;

  void WriteToHDF5();

  //Auxiliary functions to output to HDF5
  preequilibrium_event_info get_ecc(const std::vector<double> &eps,
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
                   std::vector<preequilibrium_event_info>& evt_vec,
                   std::vector<double>& entropy_density_distribution);

public:
  PreequilibriumDynamics();

  virtual ~PreequilibriumDynamics();

  /** Reads the input parameters from the XML file under the tag <Preequilibrium>. Uses JetScapeSingnalManager Instance to retrive the Initial State Physics information. Calls InitializeHydro(parameter_list) and InitTask(); This explicit call can be used for actual initialization of modules such as @a Brick, @a MpiMusic, or @a OSU-HYDRO if attached as a @a polymorphic class. It also initializes the tasks within the current module.
    @sa Read about @a polymorphism in C++. Override Init (not InitTask) here as sub-tasks are called as well.
    */
  void Init() override;

  /** Calls EvolvePreequilibrium(); This explicit call can be used for actual execution of Preequilibrium evolution defined in the modules such as @a Brick, @a MpiMusic, or @a OSU-HYDRO if attached as a @a polymorphic class. It also execute the tasks within the current module.
    @sa Read about @a polymorphism in C++.
    */
  void ExecuteTask();

  void ClearTask();

  virtual void
  InitializePreequilibrium(PreEquilibriumParameterFile parameter_list) {}
  virtual void EvolvePreequilibrium() {}

  // add initial state shared pointer
  /** A pointer of type InitialState class.
    */
  std::shared_ptr<InitialState> ini;

  PreEquilibriumParameterFile &GetParameterList() { return parameter_list_; }

  int GetPreequilibriumStatus() { return (preequilibrium_status_); }

  // @return Start time (or tau) for hydrodynamic evolution
  real GetPreequilibriumStartTime() { return (preequilibrium_tau_0_); }

  // @return End time (or tau) for hydrodynamic evolution.
  real GetPreequilibriumEndTime() { return (preequilibrium_tau_max_); }

  // record preequilibrium running status
  PreequilibriumStatus preequilibrium_status_;

  std::vector<double> e_;
  std::vector<double> P_;
  std::vector<double> utau_;
  std::vector<double> ux_;
  std::vector<double> uy_;
  std::vector<double> ueta_;
  std::vector<double> pi00_;
  std::vector<double> pi01_;
  std::vector<double> pi02_;
  std::vector<double> pi03_;
  std::vector<double> pi11_;
  std::vector<double> pi12_;
  std::vector<double> pi13_;
  std::vector<double> pi22_;
  std::vector<double> pi23_;
  std::vector<double> pi33_;
  std::vector<double> bulk_Pi_;
  double tau_hydro_;
  std::vector<double> rho_b_;
  std::vector<double> q0_;
  std::vector<double> q1_;
  std::vector<double> q2_;
  std::vector<double> q3_;

  double xmax_fs;
  double ymax_fs;
  double etamax_fs;
  double dx_fs;
  double dy_fs;
  double deta_fs;
  double neta_fs;
  bool using_ampt = false;
};

} // end namespace Jetscape

#endif
