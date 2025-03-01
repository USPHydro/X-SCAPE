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

// Create a pythia collision at a specified point and return the two inital hard partons
#include "PythiaIsrGun.h"
#include <sstream>
#include <iostream>
#include <fstream>
#define MAGENTA "\033[35m"

using namespace std;

// Register the module with the base class
RegisterJetScapeModule<PythiaIsrGun> PythiaIsrGun::reg("PythiaIsrGun");

PythiaIsrGun::~PythiaIsrGun() { VERBOSE(8); }

void PythiaIsrGun::InitTask() {


  JSDEBUG << "Initialize PythiaIsrGun";
  VERBOSE(8);

  // Show initialization at INFO level
  readString("Init:showProcesses = off");
  readString("Init:showChangedSettings = off");
  readString("Init:showMultipartonInteractions = on");
  readString("Init:showChangedParticleData = off");
  if (JetScapeLogger::Instance()->GetInfo()) {
    readString("Init:showProcesses = on");
    readString("Init:showChangedSettings = on");
    readString("Init:showMultipartonInteractions = on");
    readString("Init:showChangedParticleData = on");
  }

  // No event record printout.
  readString("Next:numberShowInfo = 0");
  readString("Next:numberShowProcess = 0");
  readString("Next:numberShowEvent = 0");

  // Standard settings 
  readString("HardQCD:all = on"); // will repeat this line in the xml for demonstration
  // readString("HardQCD:gg2gg = on");
  // readString("HardQCD:gg2qqbar = on");
  // readString("HardQCD:qg2qg = on");
  // readString("HardQCD:qq2qq = on");
  // readString("HardQCD:qqbar2gg = on");
  // readString("HardQCD:qqbar2qqbarNew = on");
  readString("HardQCD:nQuarkNew = 3"); // Number Of Quark flavours
  readString("MultipartonInteractions:processLevel = 0"); 
  readString("MultipartonInteractions:nQuarkIn = 3"); // Number Of Quark flavours
 
  // readString("HardQCD:gg2ccbar = off");
  // readString("HardQCD:qqbar2ccbar = off");
  // readString("HardQCD:hardccbar = off");
  // readString("HardQCD:gg2bbbar = off");
  // readString("HardQCD:qqbar2bbbar = off");

  //  readString("HardQCD:gg2ccbar = on"); // switch on heavy quark channel
  //readString("HardQCD:qqbar2ccbar = on");
  readString("HadronLevel:Decay = off");
  readString("HadronLevel:all = on");
  readString("PartonLevel:ISR = off");
  readString("PartonLevel:MPI = on");
  //readString("PartonLevel:FSR = on");
  readString("PromptPhoton:all=on");
  readString("WeakSingleBoson:all=off");
  readString("WeakDoubleBoson:all=off");

  // For parsing text
  stringstream numbf(stringstream::app | stringstream::in | stringstream::out);
  numbf.setf(ios::fixed, ios::floatfield);
  numbf.setf(ios::showpoint);
  numbf.precision(1);
  stringstream numbi(stringstream::app | stringstream::in | stringstream::out);

  std::string s = GetXMLElementText({"Hard", "PythiaGun", "name"});
  SetId(s);
  // cout << s << endl;

  // SC: read flag for FSR
  FSR_on = GetXMLElementInt({"Hard", "PythiaGun", "FSR_on"});
  if (FSR_on)
    readString("PartonLevel:FSR = on");
  else
    readString("PartonLevel:FSR = off");

  pTHatMin = GetXMLElementDouble({"Hard", "PythiaGun", "pTHatMin"});
  pTHatMax = GetXMLElementDouble({"Hard", "PythiaGun", "pTHatMax"});

  flag_useHybridHad = GetXMLElementInt({"Hard", "PythiaGun", "useHybridHad"});

  JSINFO << MAGENTA << "Pythia Gun with FSR_on: " << FSR_on;
  JSINFO << MAGENTA << "Pythia Gun with " << pTHatMin << " < pTHat < "
         << pTHatMax;
  JSINFO << MAGENTA << "Use hybrid hadronization? " << flag_useHybridHad;

  numbf.str("PhaseSpace:pTHatMin = ");
  numbf << pTHatMin;
  readString(numbf.str());
  numbf.str("PhaseSpace:pTHatMax = ");
  numbf << pTHatMax;
  readString(numbf.str());

  // random seed
  // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
  tinyxml2::XMLElement *RandomXmlDescription = GetXMLElement({"Random"});
  readString("Random:setSeed = on");
  numbi.str("Random:seed = ");
  unsigned int seed = 0;
  if (RandomXmlDescription) {
    tinyxml2::XMLElement *xmle =
        RandomXmlDescription->FirstChildElement("seed");
    if (!xmle)
      throw std::runtime_error("Cannot parse xml");
    xmle->QueryUnsignedText(&seed);
  } else {
    JSWARN << "No <Random> element found in xml, seeding to 0";
  }
  VERBOSE(7) << "Seeding pythia to " << seed;
  numbi << seed;
  readString(numbi.str());

  // Species
  readString("Beams:idA = 2212");
  readString("Beams:idB = 2212");

  // Energy
  eCM = GetXMLElementDouble({"Hard", "PythiaGun", "eCM"});
  numbf.str("Beams:eCM = ");
  numbf << eCM;
  readString(numbf.str());

  std::stringstream lines;
  lines << GetXMLElementText({"Hard", "PythiaGun", "LinesToRead"}, false);
  int i = 0;
  while (std::getline(lines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue; // skip empty lines
    VERBOSE(7) << "Also reading in: " << s;
    readString(s);
  }

  // And initialize
  if (!init()) { // Pythia>8.1
    throw std::runtime_error("Pythia init() failed.");
  }

    std::ofstream sigma_printer;
    sigma_printer.open(printer, std::ios::trunc);


}

void PythiaIsrGun::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  JetScapeTask::WriteTasks(w);
}

void PythiaIsrGun::ExecuteTask() {
  VERBOSE(1) << "Run Hard Process : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();
  //Reading vir_factor from xml for MATTER
  double vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});

  bool flag62 = false;
  vector<Pythia8::Particle> p62;

  // sort by pt
  struct greater_than_pt {
    inline bool operator()(const Pythia8::Particle &p1,
                           const Pythia8::Particle &p2) {
      return (p1.pT() > p2.pT());
    }
  };


    FourVector p_p;

  int NSamplings = 0;
  ReDoSampling:
  do {
    NSamplings++;
    next();
    p62.clear();
 
    std::vector<int> IndexToSkip;

    for (int parid = 0; parid < event.size(); parid++) {
      if (parid < 3)
        continue; // 0, 1, 2: total event and beams
      Pythia8::Particle &particle = event[parid];

      if (!(particle.isGluon() ||
            particle.isQuark())) { // Getting rid of diquark
        if (particle.status() == -31) {
          IndexToSkip.push_back(particle.daughter1());
          IndexToSkip.push_back(particle.daughter2());
          IndexToSkip.push_back(event[particle.daughter1()].mother1());
          IndexToSkip.push_back(event[particle.daughter1()].mother2());
        } else if (particle.status() == -33) {
          IndexToSkip.push_back(particle.mother1());
          IndexToSkip.push_back(particle.mother2());
          IndexToSkip.push_back(event[particle.mother1()].daughter1());
          IndexToSkip.push_back(event[particle.mother1()].daughter2());
        }
      }
    }

    for (auto &ToSkip : IndexToSkip) {
      JSWARN << " Skipping non-parton index " << ToSkip;
    }
      if (!printer.empty()){
            std::ofstream sigma_printer;
            sigma_printer.open(printer, std::ios::out | std::ios::app);

            sigma_printer << "sigma = " << GetSigmaGen() << " Err =  " << GetSigmaErr() << endl ;
            //sigma_printer.close();

//      JSINFO << BOLDYELLOW << " sigma = " << GetSigmaGen() << " sigma err = " << GetSigmaErr() << " printer = " << printer << " is " << sigma_printer.is_open() ;
    };

    // pTarr[0]=0.0; pTarr[1]=0.0;
    // pindexarr[0]=0; pindexarr[1]=0;

    for (int parid = 0; parid < event.size(); parid++) {
      if (parid < 3)
        continue; // 0, 1, 2: total event and beams
        
      Pythia8::Particle &particle = event[parid];
      if (!FSR_on) {
          
          if ( !( (particle.status() == -21) || (particle.status() == -23) || (particle.status() == -31) || (particle.status() == -33) )) continue ;
          for(auto &ToSkip: IndexToSkip) {
            if(ToSkip == parid){
              goto SkipParton;
            }
          }
          // if ( (particle.status() > -21)||(particle.status()<-23) ) continue ;

        // only accept particles after MPI
        //if (particle.status() != 62)
          //continue;
        // only accept gluons and quarks
        // Also accept Gammas to put into the hadron's list
        //if (fabs(particle.id()) > 5 &&
          //  (particle.id() != 21 && particle.id() != 22))
          //continue;

        // reject rare cases of very soft particles that don't have enough e to get
        // reasonable virtuality
        //if (particle.pT() < 1.0 / sqrt(vir_factor))
          //continue;

        //if(particle.id()==22) cout<<"########this is a photon!######" <<endl;
        // accept
      } else { // FSR_on true: use Pythia vacuum shower instead of MATTER
        if (!particle.isFinal())
          continue;
        // only accept gluons and quarks
        // Also accept Gammas to put into the hadron's list
        if (fabs(particle.id()) > 5 &&
            (particle.id() != 21 && particle.id() != 22))
          continue;
      }
      
        VERBOSE(1) << MAGENTA << " particle from pythiagun id = " << particle.id() << " pz = " << particle.pz() << " px = " << particle.px() << " py = " << particle.py() << " E =  " << particle.e() <<  " status = " << particle.status() << " idex "<< particle.index() << 
        " Color " << particle.col() << " " << particle.acol() << " Mothers " << particle.mother1() << " " << particle.mother2() << " daughter " << particle.daughter1() << " " << particle.daughter2();

        p62.push_back(particle);

        SkipParton:;
    }

    // if you want at least 2
    if (p62.size() < 2)
      continue;
    //if ( p62.size() < 1 ) continue;

    // Now have all candidates, sort them
    // sort by pt
    // std::sort(p62.begin(), p62.end(), greater_than_pt());
    // // check...
    // for (auto& p : p62 ) cout << p.pT() << endl;

    flag62 = true;

  } while (!flag62);


  // // Roll for a starting point
  // // See: https://stackoverflow.com/questions/15039688/random-generator-from-vector-with-probability-distribution-in-c
  // std::random_device device;
  // std::mt19937 engine(device()); // Seed the random number engine

    FourVector x_p;
    
  if (!ini)
  {
    JSINFO << BOLDYELLOW << "No initial state module, setting the starting location to "
                  "0. Make sure to add e.g. trento before PythiaIsrGun.";
  }
  else
  {
    double t, x, y, z;
    bool pass = false;
//  while (!pass)
  //  {
        ini->SampleABinaryCollisionPoint(t, x, y, z);
        // JSINFO << MAGENTA << " pass at time " << t << " with x = " << x << " y = " << y << " z = " << z << "  ? " ;
        //cin >> pass;
    //}
    x_p.Set(x,y,z,t);
      ini->OutputHardCollisionPosition(t, x, y, z);
  }
    
  // Loop through particles
  // Accept them all

  double initial_state_label = -1 ;
  double final_state_label = 1 ;
  ini->pTHat.resize((p62.size())/4);
  int hCounter = 0;
  SetTotalMomentumPositive(0.0);
  SetTotalMomentumNegative(0.0);
  SetTotalMomentumFractionPositive(0.0);
  SetTotalMomentumFractionNegative(0.0);
  double TotalEnergyOfInitialStatePartons = 0.0;

  for (int np = 0; np < p62.size(); ++np)
  // for (int np = p62.size()-1; np >= 0 ; --np)
  {
    Pythia8::Particle &particle = p62.at(np);

      if (particle.status()==-21 || particle.status()==-31 )
      {
          TotalEnergyOfInitialStatePartons += particle.e();
          if(particle.pz() >= 0.0) {
            SetTotalMomentumPositive(GetTotalMomentumPositive() + particle.e());
            SetTotalMomentumFractionPositive(GetTotalMomentumFractionPositive() + (particle.e() + particle.pz() ) / ( 0.94 * eCM));
            }
          else {
            SetTotalMomentumNegative(GetTotalMomentumNegative() +particle.e());
            SetTotalMomentumFractionNegative(GetTotalMomentumFractionNegative() + (particle.e() - particle.pz() ) / ( 0.94 * eCM));
            }
      }
  }

  VERBOSE(2) << "Negative Partons Momentum "<< GetTotalMomentumNegative()
          << " Positive Partons Momentum "<< GetTotalMomentumPositive()
          << " eCM " << eCM
          << " TotalEnergyOfInitialStatePartons = " << TotalEnergyOfInitialStatePartons;
  if(GetTotalMomentumFractionNegative() >= 1.0 || GetTotalMomentumFractionPositive() >= 1.){
    JSINFO << "Redoing Pythia Sampling Since MPI energy is larger than eCM/2.1 ";
    if(NSamplings < 1000){
      goto ReDoSampling;
    }
    JSWARN << "Negative Partons Momentum Fraction "<< GetTotalMomentumFractionNegative()
           << " Positive Partons Momentum Fraction "<< GetTotalMomentumFractionPositive()
           << " eCM " << eCM
           << " TotalEnergyOfInitialStatePartons = " << TotalEnergyOfInitialStatePartons;
    throw std::runtime_error("Pythia Isr Gun outputs more energy in the MPI partons than eCM");
  }

  for (int np = 0; np < p62.size(); ++np)
  // for (int np = p62.size()-1; np >= 0 ; --np)
  {
    Pythia8::Particle &particle = p62.at(np);

    // VERBOSE(7) << "Adding particle with pid = " << particle.id()
    //            << ", pT = " << particle.pT() << ", y = " << particle.y()
    //            << ", phi = " << particle.phi() << ", e = " << particle.e();

    // JSINFO<< MAGENTA << " at x=" << x_p.x() << ", y=" << x_p.y() << ", z=" << x_p.z() << ", t = " << x_p.t();

      int label = 0;
      int stat = 0;
      
      if (particle.status()==-21 || particle.status()==-31 )
      {
          label = initial_state_label;
          initial_state_label--;
          stat = -1000; // raw initial state status, must go to an initial state module

      }
      if (particle.status()==-23 || particle.status()==-33)
      {
          label = final_state_label;
          final_state_label++;
          stat = 1000; // raw final state status, must go to a final state module with virtuality generation. 
          if( (label-1) % 2 == 0){
            ini->pTHat[(label-1)/2] = particle.pT();
          }
      }

      FourVector p_p(particle.px(),particle.py(),particle.pz(),particle.e());
      
    if (flag_useHybridHad != 1) {
        AddParton(make_shared<Parton>(label, particle.id(), stat, p_p,x_p));
        VERBOSE(1) << BOLDYELLOW << " Pythia particle eta = " << particle.eta() << " pz = " << particle.pz() << " pT = " << particle.pT() << " phi = "  << particle.phi();
    } else {
      auto ptn = make_shared<Parton>(label, particle.id(), stat, p_p, x_p);
      ptn->set_color(particle.col());
      ptn->set_anti_color(particle.acol()); 
      ptn->set_max_color(GetMax_ColorPerShower() * (np + 1));
      AddParton(ptn);

    }
  }

  // Getting Number of hard partons
  int NPP = p62.size();
  SetMax_Color(GetMax_ColorPerShower() * NPP);
  FourVector Zeros(0,0,0,0);

  ini->CollisionNegativeMomentum = std::vector<FourVector>(NPP/2,Zeros);
  ini->CollisionPositiveMomentum = std::vector<FourVector>(NPP/2,Zeros);
  ini->CollisionNegativeRotatedMomentum = std::vector<FourVector>(NPP/2,Zeros);
  ini->CollisionPositiveRotatedMomentum = std::vector<FourVector>(NPP/2,Zeros);
  ini->ClearHardPartonMomentum();

  VERBOSE(8) << GetNHardPartons();

  //REMARK: Check why this has to be called explictly, something wrong with generic recursive execution!!????

  // JetScapeTask::ExecuteTasks();
}
