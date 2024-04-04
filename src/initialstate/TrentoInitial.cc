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

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <boost/bind.hpp>
#include <boost/tokenizer.hpp>
#include <algorithm>
#include <functional>
#include <string>
#include "JetScapeLogger.h"
#include "music/eos.h"
#include "TrentoInitial.h"

namespace Jetscape {

// Register the module with the base class
RegisterJetScapeModule<TrentoInitial> TrentoInitial::reg("TrentoInitial");

namespace {
/// @brief Tokenize a string.  The tokens will be separated by each non-quoted
///        space or equal character.  Empty tokens are removed.
///
/// @param input The string to tokenize.
///
/// @return Vector of tokens.
std::vector<std::string> tokenize(const std::string &input) {
  typedef boost::escaped_list_separator<char> separator_type;
  separator_type separator("\\",    // The escape characters.
                           "= ",    // The separator characters.
                           "\"\'"); // The quote characters.

  // Tokenize the intput.
  boost::tokenizer<separator_type> tokens(input, separator);

  // Copy non-empty tokens from the tokenizer into the result.
  std::vector<std::string> result;
  copy_if(tokens.begin(), tokens.end(), std::back_inserter(result),
          !boost::bind(&std::string::empty, _1));
  return result;
}
} // end namespace

// See header for explanation.
TrentoInitial::TrentoInitial() : InitialState() { SetId("Trento"); }

TrentoInitial::~TrentoInitial() = default;

void TrentoInitial::InitTask() {
  JSINFO << " Initialzie TRENTo initial condition ";

  // TRENTO OPTION DESK
  event_counter = 0;
  using namespace trento;
  using OptDesc = po::options_description;
  using VecStr = std::vector<std::string>;
  OptDesc main_opts{};
  main_opts.add_options()(
      "projectile",
      po::value<VecStr>()
          ->required()
          ->notifier( // use a lambda to verify there are exactly two projectiles
              [](const VecStr &projectiles) {
                if (projectiles.size() != 2)
                  throw po::required_option{"projectile"};
              }),
      "projectile symbols")("number-events", po::value<int>()->default_value(1),
                            "number of events");

  // Make all main arguments positional.
  po::positional_options_description positional_opts{};
  positional_opts.add("projectile", 2).add("number-events", 1);

  using VecPath = std::vector<fs::path>;
  OptDesc general_opts{"general options"};
  general_opts.add_options()("help,h", "show this help message and exit")(
      "version", "print version information and exit")(
      "bibtex", "print bibtex entry and exit")
      // ("default-config", "print a config file with default settings and exit")
      ("config-file,c", po::value<VecPath>()->value_name("FILE"),
       "configuration file\n(can be passed multiple times)");

  OptDesc output_opts{"output options"};
  output_opts.add_options()("quiet,q", po::bool_switch(),
                            "do not print event properties to stdout")(
      "output,o", po::value<fs::path>()->value_name("PATH"),
      "HDF5 file or directory for text files")(
      "no-header", po::bool_switch(), "do not write headers to text files");

  OptDesc phys_opts{"physical options"};
  phys_opts.add_options()(
      "reduced-thickness,p",
      po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
      "reduced thickness parameter")(
      "fluctuation,k",
      po::value<double>()->value_name("FLOAT")->default_value(1., "1"),
      "gamma fluctuation shape parameter")(
      "nucleon-width,w",
      po::value<double>()->value_name("FLOAT")->default_value(.5, "0.5"),
      "Gaussian nucleon width [fm]")(
      "nucleon-min-dist,d",
      po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
      "minimum nucleon-nucleon distance [fm]")(
      "mean-coeff,m",
      po::value<double>()->value_name("FLOAT")->default_value(1., "1."),
      "rapidity mean coefficient")(
      "std-coeff,s",
      po::value<double>()->value_name("FLOAT")->default_value(3., "3."),
      "rapidity std coefficient")(
      "skew-coeff,t",
      po::value<double>()->value_name("FLOAT")->default_value(0., "0."),
      "rapidity skew coefficient")(
      "skew-type,r", po::value<int>()->value_name("INT")->default_value(1, "1"),
      "rapidity skew type: 1: relative, 2: absolute, other: no skew")(
      "jacobian,j",
      po::value<double>()->value_name("FLOAT")->default_value(0.8, "0.8"),
      "<pt>/<mt> used in Jacobian")(
      "normalization,n",
      po::value<double>()->value_name("FLOAT")->default_value(1., "1"),
      "normalization factor");

  OptDesc coll_opts{"collision options"};
  coll_opts.add_options()(
      "beam-energy,e",
      po::value<double>()->value_name("FLOAT")->default_value(2760, "2760"),
      "collision beam energy sqrt(s) [GeV], initializes cross section")(
      "cross-section,x",
      po::value<double>()->value_name("FLOAT")->default_value(-1, "off"),
      "manual inelastic nucleon-nucleon cross section sigma_NN [fm^2]")(
      "b-min", po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
      "minimum impact parameter [fm]")(
      "b-max",
      po::value<double>()->value_name("FLOAT")->default_value(-1., "auto"),
      "maximum impact parameter [fm]")(
      "npart-min", po::value<int>()->value_name("INT")->default_value(0, "0"),
      "minimum Npart cut")("npart-max",
                           po::value<int>()->value_name("INT")->default_value(
                               std::numeric_limits<int>::max(), "INT_MAX"),
                           "maximum Npart cut")(
      "s-min", po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
      "minimum entropy cut")(
      "s-max",
      po::value<double>()->value_name("FLOAT")->default_value(
          std::numeric_limits<double>::max(), "DOUBLE_MAX"),
      "maxmimum entropy cut")(
      "random-seed",
      po::value<int64_t>()->value_name("INT")->default_value(-1, "auto"),
      "random seed")(
      "ncoll,b", po::bool_switch(),
      "calculate # of binary collision and binary collision density");

  OptDesc grid_opts{"grid options"};
  grid_opts.add_options()(
      "xy-max",
      po::value<double>()->value_name("FLOAT")->default_value(10., "10.0"),
      "xy max [fm]\n(transverse grid from -max to +max)")(
      "xy-step",
      po::value<double>()->value_name("FLOAT")->default_value(0.2, "0.2"),
      "transverse step size [fm]")(
      "eta-max",
      po::value<double>()->value_name("FLOAT")->default_value(0.0, "0.0"),
      "pseudorapidity max \n(eta grid from -max to +max)")(
      "eta-step",
      po::value<double>()->value_name("FLOAT")->default_value(0.5, "0.5"),
      "pseudorapidity step size");

  // Make a meta-group containing all the option groups except the main
  // positional options (don't want the auto-generated usage info for those).
  OptDesc usage_opts{};
  usage_opts.add(general_opts)
      .add(output_opts)
      .add(phys_opts)
      .add(coll_opts)
      .add(grid_opts);

  // Now a meta-group containing _all_ options.
  OptDesc all_opts{};
  all_opts.add(usage_opts).add(main_opts);

  // Will be used several times.
  const std::string usage_str{
      "usage: trento [options] projectile projectile [number-events = 1]\n"};
  const std::string usage_str3d{
      "To operate in 3D mode, make sure --eta-max is nonzero.\n"};

  // NOW LETS FILL IN THE OPTION DESK
  auto phy_opts = GetXMLElement({"IS", "Trento", "PhysicsInputs"});
  auto cut_opts = GetXMLElement({"IS", "Trento", "CutInputs"});
  auto trans_opts = GetXMLElement({"IS", "Trento", "TransInputs"});
  auto longi_opts = GetXMLElement({"IS", "Trento", "LongiInputs"});

  double xymax = GetXMax(), dxy = GetXStep();
  double etamax = GetZMax(), deta = GetZStep();

  auto random_seed = (*GetMt19937Generator())();
  //TEMPORARY FOR TESTING
  //auto random_seed = 1;
  //TEMPORARY
  JSINFO << "Random seed used for Trento " << random_seed;

  std::string proj(phy_opts->Attribute("projectile"));
  std::string targ(phy_opts->Attribute("target"));
  double sqrts = std::atof(phy_opts->Attribute("sqrts"));
  double cross_section = std::atof(phy_opts->Attribute("cross-section"));
  double normalization = std::atof(phy_opts->Attribute("normalization"));

  int cen_low = std::atoi(cut_opts->Attribute("centrality-low"));
  int cen_high = std::atoi(cut_opts->Attribute("centrality-high"));

  double p = std::atof(trans_opts->Attribute("reduced-thickness"));
  double k = std::atof(trans_opts->Attribute("fluctuation"));
  double w = std::atof(trans_opts->Attribute("nucleon-width"));
  double d = std::atof(trans_opts->Attribute("nucleon-min-dist"));

  double mean = std::atof(longi_opts->Attribute("mean-coeff"));
  double var = std::atof(longi_opts->Attribute("std-coeff"));
  double skew = std::atof(longi_opts->Attribute("skew-coeff"));
  int skew_type = std::atof(longi_opts->Attribute("skew-type"));
  double J = std::atof(longi_opts->Attribute("jacobian"));

  std::string options1 =
      +" --random-seed " + std::to_string(random_seed) + " --cross-section " +
      std::to_string(cross_section) + " --beam-energy " +
      std::to_string(sqrts) + " --reduced-thickness " + std::to_string(p) +
      " --fluctuation " + std::to_string(k) + " --nucleon-width " +
      std::to_string(w) + " --nucleon-min-dist " + std::to_string(d) +
      " --mean-coeff " + std::to_string(mean) + " --std-coeff " +
      std::to_string(var) + " --skew-coeff " + std::to_string(skew) +
      " --skew-type " + std::to_string(skew_type) + " --jacobian " +
      std::to_string(J) + " --quiet ";
  std::string options2 = " --normalization " + std::to_string(normalization) +
                         " --ncoll " // calcualte # of binary collision
                         + " --xy-max " + std::to_string(xymax) +
                         " --xy-step " + std::to_string(dxy) + " --eta-max " +
                         std::to_string(etamax) + " --eta-step " +
                         std::to_string(deta);
  // Handle centrality table, not normzlized, default grid, 2D (fast) !!!
  std::string cmd_basic = proj + " " + targ + " 10000 " + options1;
  VarMap var_map_basic{};
  po::store(po::command_line_parser(tokenize(cmd_basic))
                .options(all_opts)
                .positional(positional_opts)
                .run(),
            var_map_basic);

  std::string options_cut = "";
  if (cen_low == 0 && cen_high == 100) {
    JSINFO << "TRENTo Minimum Biased Mode Generates 0-100(%) of nuclear "
              "inelastic cross-section";
  } else {
    auto Ecut = GenCenTab(proj, targ, var_map_basic, cen_low, cen_high);
    double Ehigh = Ecut.first * normalization; // rescale the cut
    double Elow = Ecut.second * normalization; // rescale the cut

    JSINFO << "The total energy density cut for centrality = [" << cen_low
           << ", " << cen_high << "] (%) is:";
    JSINFO << Elow << "<dE/deta(eta=0)<" << Ehigh;
    options_cut = " --s-max " + std::to_string(Ehigh) + " --s-min " +
                  std::to_string(Elow);
    // Set trento configuration
  }
  std::string cmd =
      proj + " " + targ + " 1 " + options1 + options2 + options_cut;
  JSINFO << cmd;
  VarMap var_map{};
  po::store(po::command_line_parser(tokenize(cmd))
                .options(all_opts)
                .positional(positional_opts)
                .run(),
            var_map);
  TrentoGen_ = std::make_shared<trento::Collider>(var_map);
  SetRanges(xymax, xymax, etamax);
  SetSteps(dxy, dxy, deta);
  JSINFO << "TRENTo set";
}

bool compare_E(trento::records r1, trento::records r2) {
  return r1.mult > r2.mult;
}

std::pair<double, double> TrentoInitial::GenCenTab(std::string proj,
                                                   std::string targ,
                                                   VarMap var_map, int cL,
                                                   int cH) {
  // Terminate for nonsense
  if (cL < 0 || cL > 100 || cH < 0 || cH > 100 || cH < cL) {
    JSWARN << "Wrong centrality cuts! To be terminated.";
    exit(-1);
  }
  // These are all the parameters that could change the shape of centrality tables
  // Normalization prefactor parameter is factorized
  // They form a table header
  trento::Collider another_collider(var_map);
  double beamE = var_map["beam-energy"].as<double>();
  double xsection = var_map["cross-section"].as<double>();
  double pvalue = var_map["reduced-thickness"].as<double>();
  double fluct = var_map["fluctuation"].as<double>();
  double nuclw = var_map["nucleon-width"].as<double>();
  double dmin = var_map["nucleon-min-dist"].as<double>();
  char buffer[512];
  std::sprintf(buffer, "%s-%s-E-%1.0f-X-%1.2f-p-%1.2f-k-%1.2f-w-%1.2f-d-%1.2f",
               proj.c_str(), targ.c_str(), beamE, xsection, pvalue, fluct,
               nuclw, dmin);
  std::string header(buffer);
  JSINFO << "TRENTO centrality table header: " << header;
  // Create headering string hash tage for these parameter combination
  // Use this tag as a unique table filename for this specific parameter set
  std::hash<std::string> hash_function;
  size_t header_hash = hash_function(header);
  JSINFO << "Hash tag for this header: " << header_hash;
  // create dir incase it does not exist
  std::system("mkdir -p ./trento_data");
  char filename[512];
  std::sprintf(filename, "./trento_data/%zu", header_hash);
  // Step1: check it a table exist
  std::ifstream infile(filename);
  double Etab[101];
  double buff1, buff2;
  std::string line;
  if (infile.good()) {
    JSINFO << "The required centrality table exists. Load the table.";
    int i = 0;
    while (std::getline(infile, line)) {
      if (line[0] != '#') {
        std::istringstream iss(line);
        iss >> buff1 >> buff2 >> Etab[i];
        i++;
      }
    }
    infile.close();
  } else {
    JSINFO << "TRENTo is generating new centrality table for this new "
              "parameter set";
    JSINFO << "It may take 10(s) to 1(min).";

    another_collider.run_events();
    // Get all records and sort according to totoal energy
    auto event_records = another_collider.all_records();
    std::sort(event_records.begin(), event_records.end(), compare_E);
    // write centrality table
    int nstep = std::ceil(event_records.size() / 100);
    std::ofstream fout(filename);
    fout << "#\tproj\ttarj\tsqrts\tx\tp\tk\tw\td\n"
         << "#\t" << proj << "\t" << targ << "\t" << beamE << "\t" << xsection
         << "\t" << pvalue << "\t" << fluct << "\t" << nuclw << "\t" << dmin
         << "\n"
         << "#\tcen_L\tcen_H\tun-normalized total density\n";
    Etab[0] = 1e10;
    for (int i = 1; i < 100; i += 1) {
      auto ee = event_records[i * nstep];
      fout << i - 1 << "\t" << i << "\t" << ee.mult << std::endl;
      Etab[i] = ee.mult;
    }
    auto ee = event_records.back();
    fout << 99 << "\t" << 100 << "\t" << ee.mult << std::endl;
    Etab[100] = ee.mult;
    fout.close();
  }
  JSINFO << "#########" << Etab[cL] << " " << Etab[cH];
  return std::make_pair(Etab[cL], Etab[cH]);
}

void TrentoInitial::ExecuteTask() {
  JSINFO << " Exec TRENTo initial condition ";
  TrentoGen_->run_events();

  JSINFO << " TRENTo event info: ";
  auto tmp_event = TrentoGen_->expose_event();
  info_.impact_parameter = TrentoGen_->all_records().back().b;
  info_.num_participant = tmp_event.npart();
  info_.num_binary_collisions = tmp_event.ncoll();
  info_.total_entropy = tmp_event.multiplicity();
  info_.ecc = tmp_event.eccentricity();
  info_.psi = tmp_event.participant_plane();
  info_.xmid =
      -GetXMax() + tmp_event.mass_center_index().first * tmp_event.dxy();
  info_.ymid =
      -GetYMax() + tmp_event.mass_center_index().second * tmp_event.dxy();
  JSINFO << "b\tnpart\tncoll\tET\t(x-com, y-com) (fm)";
  JSINFO << info_.impact_parameter << "\t" << info_.num_participant << "\t"
         << info_.num_binary_collisions << "\t" << info_.total_entropy << "\t"
         << "(" << info_.xmid << ", " << info_.ymid << ")";

  JSINFO << " Load TRENTo density and ncoll density to JETSCAPE memory ";
  auto density_field = tmp_event.density_grid();
  auto ncoll_field = tmp_event.TAB_grid();
  JSINFO << density_field.num_elements() << " density elements";
  for (int i = 0; i < density_field.num_elements(); i++) {
    entropy_density_distribution_.push_back(density_field.data()[i]);
  }
  JSINFO << ncoll_field.num_elements() << " ncoll elements";
  for (int i = 0; i < ncoll_field.num_elements(); i++) {
    num_of_binary_collisions_.push_back(ncoll_field.data()[i]);
  }
  JSINFO << " TRENTO event generated and loaded ...";

  WriteToHDF5();
}

void TrentoInitial::ClearTask() {
  VERBOSE(2) << " : Finish creating initial condition ";
  entropy_density_distribution_.clear();
  num_of_binary_collisions_.clear();
}

void TrentoInitial::WriteToHDF5() {
  // Write out the initial condition
  std::vector<double> trento_profile = GetEntropyDensityDistribution();
  std::vector<double> ic_u0 = trento_profile;
  std::vector<double> ic_ux = trento_profile;
  std::vector<double> ic_uy = trento_profile;
  std::vector<double> ic_Pi = trento_profile;
  std::vector<double> ic_pixx = trento_profile;
  std::vector<double> ic_piyy = trento_profile;
  std::vector<double> ic_pixy = trento_profile;

  int nx = GetXSize();
  int ny = GetYSize();
  double dx = GetXStep();
  double dy = GetYStep();

  //Fluid initially at rest and without viscous components
  for (int ix=0; ix<nx;++ix)
  for (int iy=0; iy<ny;++iy) {
    ic_u0[iy+ix*ny] = 1.;
    ic_ux[iy+ix*ny] = 0.;
    ic_uy[iy+ix*ny] = 0.;
    ic_Pi[iy+ix*ny] = 0.;
    ic_pixx[iy+ix*ny] = 0.;
    ic_piyy[iy+ix*ny] = 0.;
  };


  std::array<int,2> size = {nx,ny};
  std::array<double,2> step = {dx,dy};
  std::array<double,2> step_fs = {dx,dy};
  trento_event_info results_trento = get_ecc(trento_profile,ic_u0, ic_ux, ic_uy, ic_Pi,
                                  ic_pixx, ic_piyy, ic_piyy, ic_pixy, 1.0, size, step);

  std::string filename = GetXMLElementText({"outputFilename"})+"_"+std::to_string(event_counter)+".ic.h5";
  std::vector<trento_event_info> evt_vec = {results_trento};

  JSINFO << "Writing IC stats to HDF5 file: " << filename;

  output_hdf5(filename, evt_vec, trento_profile);
  JSINFO << "Finished writing IC stats to HDF5 file: " << filename;
  event_counter++;

} // end namespace Jetscape

//Output the event to HDF5
void TrentoInitial::output_hdf5(std::string filename,
    std::vector<trento_event_info>& evt_vec,
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
    H5::CompType mtype(sizeof(trento_event_info));
    mtype.insertMember("b", HOFFSET(trento_event_info,b), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("npart", HOFFSET(trento_event_info, Npart), H5::PredType::NATIVE_INT);
    mtype.insertMember("ncoll", HOFFSET(trento_event_info, Ncoll), H5::PredType::NATIVE_INT);
    mtype.insertMember("E", HOFFSET(trento_event_info, E), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("S", HOFFSET(trento_event_info, S), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("S_hotQCD", HOFFSET(trento_event_info, S_hotQCD), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("E_entropy", HOFFSET(trento_event_info, E_entropy), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("re_ecc_p", HOFFSET(trento_event_info, re_ecc_p ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("im_ecc_p", HOFFSET(trento_event_info, im_ecc_p ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("re_ecc_p_hotQCD", HOFFSET(trento_event_info, re_ecc_p_hotQCD ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("im_ecc_p_hotQCD", HOFFSET(trento_event_info, im_ecc_p_hotQCD ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("re_ecc_p_prime", HOFFSET(trento_event_info, re_ecc_p_prime ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("im_ecc_p_prime", HOFFSET(trento_event_info, im_ecc_p_prime ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("re_ecc_p_prime_hotQCD", HOFFSET(trento_event_info, re_ecc_p_prime_hotQCD ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("im_ecc_p_prime_hotQCD", HOFFSET(trento_event_info, im_ecc_p_prime_hotQCD ), H5::PredType::NATIVE_DOUBLE);
    mtype.insertMember("R", HOFFSET(trento_event_info, R ), H5::ArrayType(H5::PredType::NATIVE_DOUBLE, 1, dim_R));
    mtype.insertMember("eps", HOFFSET(trento_event_info, eps ), H5::ArrayType(H5::PredType::NATIVE_DOUBLE, 2, dim_ecc));
    mtype.insertMember("psi", HOFFSET(trento_event_info, psi ), H5::ArrayType(H5::PredType::NATIVE_DOUBLE, 2, dim_ecc));

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


trento_event_info TrentoInitial::get_ecc(const std::vector<double> &eps,
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

    trento_event_info results;
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
}
