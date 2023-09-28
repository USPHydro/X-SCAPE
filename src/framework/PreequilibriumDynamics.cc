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

#include "PreequilibriumDynamics.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include <string>

#include <iostream>
#include "music/eos.h"

using namespace std;

#define MAGENTA "\033[35m"

namespace Jetscape {

PreequilibriumDynamics::PreequilibriumDynamics() {
  VERBOSE(8);
  SetId("PreequilibriumDynamics");
}

PreequilibriumDynamics::~PreequilibriumDynamics() {
  VERBOSE(8);
  disconnect_all();
}

void PreequilibriumDynamics::Init() {
  JetScapeModuleBase::InitTask();
  event_counter = 0;

  JSINFO << "Initialize PreequilibriumDynamics : " << GetId() << " ...";

  VERBOSE(8);

  // this is grabbing the initial entropy density ?
  ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  if (!ini) {
    JSWARN << "No initialization module, try: "
           << "auto trento = make_shared<TrentoInitial>(); "
           << "jetscape->Add(trento);";
  }

  InitializePreequilibrium(parameter_list_);

  InitTask();
  InitTasks();
}

void PreequilibriumDynamics::ExecuteTask() {
  VERBOSE(2) << "Run Preequilibrium : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();

  if (ini) {
    VERBOSE(3) << "length of entropy density vector="
               << ini->GetEntropyDensityDistribution().size();
  }

  EvolvePreequilibrium();
  WriteToHDF5();
}

void PreequilibriumDynamics::ClearTask() {
  e_.clear();
  P_.clear();
  utau_.clear();
  ux_.clear();
  uy_.clear();
  ueta_.clear();
  pi00_.clear();
  pi01_.clear();
  pi02_.clear();
  pi03_.clear();
  pi11_.clear();
  pi12_.clear();
  pi13_.clear();
  pi22_.clear();
  pi23_.clear();
  pi33_.clear();
  bulk_Pi_.clear();
  rho_b_.clear();
  q0_.clear();
  q1_.clear();
  q2_.clear();
  q3_.clear();
}


void PreequilibriumDynamics::WriteToHDF5() {

  //Assumes square grid
  int nx = ini->GetXSize();
  int ny = ini->GetYSize();
  double dx = ini->GetXStep();
  double dy = ini->GetYStep();

  std::vector<double> u0;
  for (int ix=0; ix<nx; ++ix){
    for (int iy=0; iy<ny; ++iy){
      int index = iy+ny*ix;
      u0.push_back(sqrt(1+pow(ux_[index],2)+pow(uy_[index],2)+pow(tau_hydro_*ueta_[index],2)));
    }
  }

  std::array<int,2> size = {nx,ny};
  std::array<double,2> step = {dx,dy};

  preequilibrium_event_info results_trento = get_ecc(e_,utau_, ux_, uy_, bulk_Pi_,
                                  pi11_, pi22_, pi00_, pi12_, tau_hydro_, size, step);

  std::string filename = GetXMLElementText({"outputFilename"})+"_"+std::to_string(event_counter)+".fs.h5";
  std::vector<preequilibrium_event_info> evt_vec = {results_trento};

  JSINFO << "Writing IC stats to HDF5 file: " << filename;

  output_hdf5(filename, evt_vec, e_);
  JSINFO << "Finished writing IC stats to HDF5 file: " << filename;
  event_counter++;

} // end namespace Jetscape

//Output the event to HDF5
void PreequilibriumDynamics::output_hdf5(std::string filename,
    std::vector<preequilibrium_event_info>& evt_vec,
    std::vector<double>& entropy_density_distribution){
    //Prep the HDF5 file structure
    hsize_t dim[1];
    dim[0] = evt_vec.size();
    // the length of dim
    int rank = sizeof(dim) / sizeof(hsize_t);

    hsize_t entropy_size[1] = {entropy_density_distribution.size()};

    //Array sizes
    hsize_t dim_R[1] = {NMAX};
    hsize_t dim_ecc[2] = {NMAX, MMAX+1};

    // defining the datatype to pass HDF55
    H5::CompType mtype(sizeof(preequilibrium_event_info));
    mtype.insertMember("b", HOFFSET(preequilibrium_event_info,b), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("npart", HOFFSET(preequilibrium_event_info, Npart), H5::PredType::NATIVE_INT);
    mtype.insertMember("ncoll", HOFFSET(preequilibrium_event_info, Ncoll), H5::PredType::NATIVE_INT);
    mtype.insertMember("E", HOFFSET(preequilibrium_event_info, E), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("S", HOFFSET(preequilibrium_event_info, S), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("S_hotQCD", HOFFSET(preequilibrium_event_info, S_hotQCD), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("E_entropy", HOFFSET(preequilibrium_event_info, E_entropy), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("re_ecc_p", HOFFSET(preequilibrium_event_info, re_ecc_p ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("im_ecc_p", HOFFSET(preequilibrium_event_info, im_ecc_p ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("re_ecc_p_hotQCD", HOFFSET(preequilibrium_event_info, re_ecc_p_hotQCD ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("im_ecc_p_hotQCD", HOFFSET(preequilibrium_event_info, im_ecc_p_hotQCD ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("re_ecc_p_prime", HOFFSET(preequilibrium_event_info, re_ecc_p_prime ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("im_ecc_p_prime", HOFFSET(preequilibrium_event_info, im_ecc_p_prime ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("re_ecc_p_prime_hotQCD", HOFFSET(preequilibrium_event_info, re_ecc_p_prime_hotQCD ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("im_ecc_p_prime_hotQCD", HOFFSET(preequilibrium_event_info, im_ecc_p_prime_hotQCD ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("R", HOFFSET(preequilibrium_event_info, R ), H5::ArrayType(H5::PredType::NATIVE_DOUBLE, 1, dim_R));
    mtype.insertMember("eps", HOFFSET(preequilibrium_event_info, eps ), H5::ArrayType(H5::PredType::NATIVE_DOUBLE, 2, dim_ecc));
    mtype.insertMember("psi", HOFFSET(preequilibrium_event_info, psi ), H5::ArrayType(H5::PredType::NATIVE_DOUBLE, 2, dim_ecc));

    // preparation of a dataset and a file.
    H5::DataSpace space(rank, dim);
    H5::DataSpace space_entropy(rank, entropy_size);
    H5::H5File *file = new H5::H5File(filename, H5F_ACC_TRUNC);
    H5::DataSet *dataset = new H5::DataSet(file->createDataSet("events", mtype, space));
    H5::DataSet *entropy_profile = new H5::DataSet(file->createDataSet("entropy_profile",
                                                   H5::ArrayType(H5::PredType::NATIVE_DOUBLE, 1, entropy_size ), space));
    // Write
    dataset->write(evt_vec.data(), mtype);
    entropy_profile->write(entropy_density_distribution.data(), H5::ArrayType(H5::PredType::NATIVE_DOUBLE, 1, entropy_size ));

    // Close
    file->close();
    delete dataset;
    delete file;
}


preequilibrium_event_info PreequilibriumDynamics::get_ecc(const std::vector<double> &eps,
                                    const  std::vector<double> &u0,
                                    const  std::vector<double> &ux,
                                    const  std::vector<double> &uy,
                                    const  std::vector<double> &Pi,
                                    const  std::vector<double> &pixx,
                                    const  std::vector<double> &piyy,
                                    const  std::vector<double> &pitautau,
                                    const  std::vector<double> &pixy,
                                    const  double &tau,
                                    const  std::array<int,2> &size, const std::array<double,2> &step){


    //---------------------
    // Initializations
    //---------------------

    ///EOS codes
    // 0: Ideal gas
    // 1: eosQ
    // 2: s95p-v1 (Huovinen/Petreczky)
    // 3: s95p with partial chemical equilibrium and chem. f.o. at 150 MeV
    // 4: s95p with partial chemical equilibrium and chem. f.o. at 155 MeV
    // 5: s95p with partial chemical equilibrium and chem. f.o. at 160 MeV
    // 6: s95p with partial chemical equilibrium and chem. f.o. at 165 MeV
    // 7: s95p-v1.2 (for UrQMD)
    // 8: WB parametrizarion
    // 9: hot_qcd
    // 91: hot_qcd with SMASH hadron ressonance gas
    // 10: neos
    // 11: neos3
    // 12: neos_b
    // 13: neos_bs
    // 14: neos_bqs
    // 15: neos_bqs (muB = 0.9)
    // 17: BEST
    // 19: EOS at finite muB from UH Collaboration

    const int eos_code = 91;
    EOS eos(eos_code);
    EOS eos_ideal(0);
    /// IMPORTANT: MUSIC eos class uses fm^-4 units, while we use GeV/fm^3

    //Will store the results for numerator and denominator of <ecc_{n,m}>
    std::array<std::complex<double>,NMAX*(MMAX+1)> num;
    std::array<double,NMAX> den;
    for (int n = 1; n<=NMAX; ++n){
        den[n-1] = 0;
        for (int m = 0; m<=MMAX; ++m){
            int ind = m + (MMAX+1)*(n-1);
            num[ind] = std::complex<double>(0,0);
        }
    }

    double total_energy = 0; //Total energy
    double total_ideal_entropy = 0; //Total entropy (using ideal EOS)
    double total_eos_entropy = 0; //Total entropy (using hotQCD EOS)
    double total_energy_p = 0; //Total energy (computed assuming input is entropy density and using ideal EOS)

    double ecc_p_conformal = 0; //Momentum anisotropy (uses only ideal Tmunu and conformal EOS)
    double ecc_p = 0; //Momentum anisotropy (uses only ideal Tmunu and HotQCD EOS)
    double ecc_p_prime_conformal = 0; //Momentum anisotropy (uses full Tmunu and conformal EOS)
    double ecc_p_prime = 0; //Momentum anisotropy (uses full Tmunu and HotQCD EOS)

    std::complex<double> num_ecc_p_conformal = 0;
    std::complex<double> num_ecc_p_eos = 0;
    std::complex<double> num_ecc_p_prime_eos = 0;
    std::complex<double> num_ecc_p_prime_conformal = 0;

    double den_ecc_p_conformal = 0;
    double den_ecc_p_eos = 0;
    double den_ecc_p_prime_eos = 0;
    double den_ecc_p_prime_conformal = 0;

    //Auxiliary variables
    std::complex<double> mean_pos(0,0);
    double const vol = step[0]*step[1]; //cell volume
    //---------------------
    // A) Assume we have energy density profile
    //---------------------

    //Get total energy and mean center

    for(int ix=0; ix<size[0]; ++ix){
        double x = (ix-size[0]/2.)*step[0];
        for(int iy=0; iy<size[1]; ++iy){
            double y = (iy-size[1]/2.)*step[1];
            int index = iy+size[1]*ix;
            double e = eps[index];
            double gamma = u0[index];
            total_energy += e*gamma;
            total_ideal_entropy += eos_ideal.get_entropy(e/Pythia8::HBARC,0)*gamma;
            total_eos_entropy += eos.get_entropy(e/Pythia8::HBARC,0.)*gamma;
            total_energy_p += Pythia8::HBARC*eos_ideal.get_s2e(e,0)*gamma;
            mean_pos += std::complex<double>(x,y)*e*gamma;
        }
    }
    mean_pos *= vol;
    total_energy *= vol;
    total_ideal_entropy *= vol;
    total_eos_entropy *= vol;
    total_energy_p *= vol;
    mean_pos /= total_energy;


    double E_tot = 0;
    //Compute <R^n> and <R^n e^{i m phi}> for n=1..NMAX and m=0..MMAX
    for(int ix=0; ix<size[0]; ++ix){
        double x = (ix-size[0]/2.)*step[0];
        for(int iy=0; iy<size[1]; ++iy){
            double y = (iy-size[1]/2.)*step[1];
            int index = iy+size[1]*ix;
            std::complex<double> z(x,y);
            z = z - mean_pos; //Recenter event
            double e = eps[index];
            double P_confomal = Pythia8::HBARC*eos_ideal.get_pressure(e/Pythia8::HBARC,0.); // e/3
            double P = Pythia8::HBARC*eos.get_pressure(e/Pythia8::HBARC,0.);
            double r = std::abs(z);
            double phi = std::arg(z);
            for (int n = 1; n<=NMAX; ++n){
                double aux = pow(r,n)*eps[index]*u0[index];
                den[n-1] += aux;
                for (int m = 0; m<=MMAX; ++m){
                    int ind =  m + (MMAX+1)*(n-1);
                    num[ind] += exp(std::complex<double>(0,m*phi))*aux;
               }
            }

            E_tot += (eps[index]+P_confomal + Pi[index] )*u0[index]*u0[index]-P_confomal - Pi[index] + pitautau[index];

            //Ideal Txx and Tyy using conformal eos
            double Txx_ideal_conformal = 4*(eps[index])*pow(ux[index],2)/3 + eps[index]/3;
            double Tyy_ideal_conformal = 4*(eps[index])*pow(uy[index],2)/3 + eps[index]/3;
            double Txy_ideal_conformal = 4*(eps[index])*ux[index]*uy[index]/3;

            //Ideal Txx and Tyy using MUSIC eos
            double Txx_ideal_eos = (eps[index] + P)*pow(ux[index],2) - eps[index]/3;
            double Tyy_ideal_eos = (eps[index] + P)*pow(uy[index],2) - eps[index]/3;
            double Txy_ideal_eos = 4*(eps[index] + P)*ux[index]*uy[index];

            //Txx and Tyy using conformal eos
            double Txx_conformal = Txx_ideal_conformal + Pi[index]*(1+pow(ux[index],2)) + pixx[index];
            double Tyy_conformal = Tyy_ideal_conformal + Pi[index]*(1+pow(uy[index],2)) + piyy[index];
            double Txy_conformal = Tyy_ideal_conformal + Pi[index]*ux[index]*uy[index] + pixy[index];

            //Txx and Tyy using MUSIC eos
            double Txx_eos = Txx_ideal_eos + Pi[index]*(1+pow(ux[index],2)) + pixx[index];
            double Tyy_eos = Txx_ideal_eos + Pi[index]*(1+pow(uy[index],2)) + piyy[index];
            double Txy_eos = Txy_ideal_eos + Pi[index]*ux[index]*uy[index] + pixy[index];

            num_ecc_p_conformal += std::complex<double>(Txx_ideal_conformal - Tyy_ideal_conformal,
                                                2*Txy_ideal_conformal);
            num_ecc_p_eos += std::complex<double>(Txx_ideal_eos - Tyy_ideal_eos,
                                                2*Txy_ideal_eos);
            num_ecc_p_prime_conformal += std::complex<double>(Txx_conformal - Tyy_conformal,
                                                2*Txy_conformal);
            num_ecc_p_prime_eos += std::complex<double>(Txx_eos - Tyy_eos,
                                                2*Txy_eos);

            den_ecc_p_conformal += Txx_ideal_conformal + Tyy_ideal_conformal;
            den_ecc_p_eos += Txx_ideal_eos + Tyy_ideal_eos;
            den_ecc_p_prime_conformal += Txx_conformal + Tyy_conformal;
            den_ecc_p_prime_eos += Txx_eos + Tyy_eos;
        }
    }
    E_tot *= vol;
    total_energy_p *= vol;

    preequilibrium_event_info results;
    results.E = E_tot*tau;
    results.S = total_ideal_entropy*tau;
    results.E_entropy = total_energy_p*tau;
    results.S_hotQCD = total_eos_entropy*tau;
    results.re_ecc_p = std::real(num_ecc_p_conformal/den_ecc_p_conformal);
    results.im_ecc_p = std::imag(num_ecc_p_conformal/den_ecc_p_conformal);
    results.re_ecc_p_hotQCD = std::real(num_ecc_p_eos/den_ecc_p_eos);
    results.re_ecc_p_hotQCD = std::imag(num_ecc_p_eos/den_ecc_p_eos);
    results.re_ecc_p_prime = std::real(num_ecc_p_prime_conformal/den_ecc_p_prime_conformal);
    results.im_ecc_p_prime = std::imag(num_ecc_p_prime_conformal/den_ecc_p_prime_conformal);
    results.re_ecc_p_prime_hotQCD = std::real(num_ecc_p_prime_eos/den_ecc_p_prime_eos);
    results.im_ecc_p_prime_hotQCD = std::imag(num_ecc_p_prime_eos/den_ecc_p_prime_eos);
    for (int n=1; n <= NMAX; ++n){
        results.R[n-1] = pow(den[(n-1)]*vol/total_energy,1./n);
        for (int m=0; m <= MMAX; ++m){
            int ind =  m + (MMAX+1)*(n-1);
            std::complex<double> ecc = num[ind]/den[n-1];
            results.eps[n-1][m] = abs(ecc);
            results.psi[n-1][m] = arg(ecc)+M_PI/m;
        }
    }
    return results;
  }  

} // end namespace Jetscape
