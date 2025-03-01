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

#include "JetScape.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include "JetEnergyLossManager.h"
#include "FluidDynamics.h"
#include "XscapeBanner.h"
#include "InitialState.h"
#include "PreequilibriumDynamics.h"
#include "JetEnergyLoss.h"
#include "CausalLiquefier.h"
#include "SoftParticlization.h"
#include "QueryHistory.h"

#ifdef USE_HEPMC
  #include "JetScapeWriterHepMC.h"
  #ifdef USE_ROOT
    #include "JetScapeWriterRootHepMC.h"
  #endif
  #ifdef USE_ROOT
    #include "JetScapeWriterRootHepMC.h"
  #endif
#endif

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

namespace Jetscape {

/** Default constructor to create the main task of the JetScape framework. It sets the total number of events to 1.
   * By default, hydro events are used only once
   */
JetScape::JetScape()
    : JetScapeModuleBase(), n_events(1), n_events_printout(100), reuse_hydro_(false), n_reuse_hydro_(1),
      liquefier(nullptr), fEnableAutomaticTaskListDetermination(true) {
  VERBOSE(8);
  SetId("primary");
}

JetScape::~JetScape() {
  VERBOSE(8);
  //JetScapeSignalManager::Instance()->Clear();
  //not needed, use weak_ptr in JetScapeSignalManager class (=not owning)
}

void JetScape::Show() { ShowXscapeBanner(); }

//________________________________________________________________
void JetScape::Init() {
  Show();

  JSINFO << BOLDRED << "Initialize JetScape ...";

  // Set the names of the XML files in the JetScapeXML instance
  JetScapeXML::Instance()->OpenXMLMainFile(GetXMLMainFileName());
  JetScapeXML::Instance()->OpenXMLUserFile(GetXMLUserFileName());
  JSINFO << "================================================================";

  // Check whether XML elements in the User file are not included in the Main file
  CompareElementsFromXML();

  // Read some general parameters from the XML configuration file
  ReadGeneralParametersFromXML();

  // Loop through the XML User file elements to determine the task list, if enabled
  if (fEnableAutomaticTaskListDetermination) {
    DetermineTaskListFromXML();
    DetermineWritersFromXML();
    JSINFO
        << "================================================================";
  }

  QueryHistory::Instance()->AddMainTask(shared_from_this());

  // Has to be called explicitly since not really fully recursively (if ever needed)
  // So --> JetScape is "Task Manager" of all modules ...
  JSINFO << "Found " << GetNumberOfTasks() << " Modules Initialize them ... ";
  SetPointers();
  JSINFO << "Calling JetScape InitTasks()...";
  JetScapeTask::InitTasks();
}

//________________________________________________________________
void JetScape::recurseToBuild(std::vector<std::string> &elems, tinyxml2::XMLElement *mElement)
{
  tinyxml2::XMLElement *nextElement = mElement->FirstChildElement();

  if (nextElement != nullptr)
    recurseToBuild(elems, nextElement);

  nextElement = mElement->NextSiblingElement();

  if (nextElement != nullptr)
    recurseToBuild(elems, nextElement);

  elems.push_back(mElement->Name());
}

//________________________________________________________________
void JetScape::recurseToSearch(std::vector<std::string> &elems, tinyxml2::XMLElement *uElement)
{
  tinyxml2::XMLElement *nextElement = uElement->FirstChildElement();

  if (nextElement != nullptr)
    recurseToSearch(elems, nextElement);

  nextElement = uElement->NextSiblingElement();

  if (nextElement != nullptr)
    recurseToSearch(elems, nextElement);

  if (!std::binary_search(elems.begin(), elems.end(), uElement->Name())) {
    JSWARN << "User XML tag <" << uElement->Name() << "> is unrecognized. A default entry in the main XML file is required.";
    exit(-1);
  }
}

//________________________________________________________________
void JetScape::CompareElementsFromXML() {

  tinyxml2::XMLElement *uElement = JetScapeXML::Instance()->GetXMLRootUser()->FirstChildElement();
  tinyxml2::XMLElement *mElement = JetScapeXML::Instance()->GetXMLRootMain()->FirstChildElement();

  std::vector<std::string> elems;

  recurseToBuild(elems, mElement);
  sort(elems.begin(), elems.end());
  recurseToSearch(elems, uElement);
}

//________________________________________________________________
void JetScape::ReadGeneralParametersFromXML() {

  // Debug level
  std::string log_debug = GetXMLElementText({"debug"});
  if ((int)log_debug.find("off") >= 0)
    JetScapeLogger::Instance()->SetDebug(false);
  VERBOSE(1) << "JetScape Debug = " << log_debug;

  // Remark
  std::string log_remark = GetXMLElementText({"remark"});
  if ((int)log_remark.find("on") >= 0)
    JetScapeLogger::Instance()->SetRemark(true);
  VERBOSE(1) << "JetScape Remark = " << log_remark;

  // Verbose level
  int m_vlevel = GetXMLElementInt({"vlevel"});
  if (m_vlevel > 0) {
    JetScapeLogger::Instance()->SetVerboseLevel(m_vlevel);
    VERBOSE(1) << "JetScape Verbose Level = " << m_vlevel;
  }

  // Flag for automatic task list determination from User XML
  std::string enableAutomaticTaskListDetermination =
      GetXMLElementText({"enableAutomaticTaskListDetermination"});
  if ((int)enableAutomaticTaskListDetermination.find("true") >= 0) {
    fEnableAutomaticTaskListDetermination = true;
    VERBOSE(1)
        << "Enable automatic task list determination from User XML: True.";
  } else if ((int)enableAutomaticTaskListDetermination.find("false") >= 0) {
    fEnableAutomaticTaskListDetermination = false;
    VERBOSE(1)
        << "Enable automatic task list determination from User XML: False.";
  }

  // nEvents
  int nEvents = GetXMLElementInt({"nEvents"});
  if (nEvents) {
    SetNumberOfEvents(nEvents);
    JSINFO << "nEvents = " << nEvents;
  }
  n_events_printout = GetXMLElementInt({"nEvents_printout"});

  // Set whether to reuse hydro
  std::string reuseHydro = GetXMLElementText({"setReuseHydro"});
  if ((int)reuseHydro.find("true") >= 0) {
    SetReuseHydro(true);
    JSINFO << "Reuse Hydro: " << reuseHydro;
  }
  int nReuseHydro = GetXMLElementInt({"nReuseHydro"});
  if (nReuseHydro) {
    SetNReuseHydro(nReuseHydro);
    JSINFO << "nReuseHydro: " << nReuseHydro;
  }

  // Set up helper. Mostly used for random numbers
  // Needs the XML reader singleton set up
  JetScapeTaskSupport::ReadSeedFromXML();

  JSDEBUG << "JetScape Debug from XML = " << log_debug;
  JSDEBUG << "JetScape Remark from XML = " << log_remark;
}

//________________________________________________________________
void JetScape::DetermineTaskListFromXML() {

  // First, check for Liquefier and create it if so (since it needs to be passed to other modules)
  VERBOSE(2) << "Checking if Liquifier should be created...";
  tinyxml2::XMLElement *elementXML =
      (tinyxml2::XMLElement *)JetScapeXML::Instance()
          ->GetXMLRootUser()
          ->FirstChildElement();
  while (elementXML) {
    std::string elementName = elementXML->Name();
    if (elementName == "Liquefier") {
      liquefier = make_shared<CausalLiquefier>();
      JSINFO << "Created liquefier.";
    }
    elementXML = elementXML->NextSiblingElement();
  }

  // Loop through and create all modules
  tinyxml2::XMLElement *element =
      (tinyxml2::XMLElement *)JetScapeXML::Instance()
          ->GetXMLRootUser()
          ->FirstChildElement();
  while (element) {
    std::string elementName = element->Name();
    VERBOSE(2) << "Parsing element: " << elementName;

    // Initial state
    if (elementName == "IS") {

      tinyxml2::XMLElement *childElement =
          (tinyxml2::XMLElement *)element->FirstChildElement();
      while (childElement) {
        std::string childElementName = childElement->Name();
        VERBOSE(2) << "Parsing childElement: " << childElementName;

        //   - Trento
        if (childElementName == "Trento") {
          auto trento = JetScapeModuleFactory::createInstance("TrentoInitial");
          if (trento) {
            Add(trento);
            JSINFO << "JetScape::DetermineTaskList() -- Initial State: Added "
                      "Trento module to task list.";
          }
        }
        //   - Read initial conditions from file
        else if (childElementName == "initial_profile_path") {
#ifdef USE_HDF5
          auto initial =
              JetScapeModuleFactory::createInstance("InitialFromFile");
          if (initial) {
            Add(initial);
            JSINFO << "JetScape::DetermineTaskList() -- Initial state: Added "
                      "InitialFromFile to task list.";
          }
#else
          JSWARN << "InitialFromFile is attempted to be added, but HDF5 is not "
                    "installed!";
#endif
        } else if (childElementName == "IPGlasma") {
          auto ipglasma = JetScapeModuleFactory::createInstance("IPGlasma");
          if (ipglasma) {
            Add(ipglasma);
            JSINFO << "JetScape::DetermineTaskList() -- Initial State: Added "
                      "IPGlasma module to task list.";
          }
        } else if (childElementName == "initial_Ncoll_list") {
          auto initial =
              JetScapeModuleFactory::createInstance("NcollListFromFile");
          if (initial) {
            Add(initial);
            JSINFO << "JetScape::DetermineTaskList() -- Initial state: Added "
                      "NcollListFromFile to task list.";
          }
        }
        //   - Custom module
        else if (((int)childElementName.find("CustomModule") >= 0)) {
          auto customModule =
              JetScapeModuleFactory::createInstance(childElementName);
          if (customModule) {
            Add(customModule);
            JSINFO << "JetScape::DetermineTaskList() -- Initial state: Added "
                   << childElementName << " to task list.";
          }
        }

        childElement = childElement->NextSiblingElement();
      }
    }

    // Hard process
    else if (elementName == "Hard") {

      tinyxml2::XMLElement *childElement =
          (tinyxml2::XMLElement *)element->FirstChildElement();
      while (childElement) {
        std::string childElementName = childElement->Name();
        VERBOSE(2) << "Parsing childElement: " << childElementName;

        //   - PGun
        if (childElementName == "PGun") {
          auto pGun = JetScapeModuleFactory::createInstance(childElementName);
          if (pGun) {
            Add(pGun);
            JSINFO << "JetScape::DetermineTaskList() -- Hard Process: Added "
                      "PGun to task list.";
          }
        }
        //   - PythiaGun
        else if (childElementName == "PythiaGun") {
          auto pythiaGun =
              JetScapeModuleFactory::createInstance(childElementName);
          if (pythiaGun) {
            Add(pythiaGun);
            JSINFO << "JetScape::DetermineTaskList() -- Hard Process: Added "
                      "PythiaGun to task list.";
          }
        } else if (((int)childElementName.find("CustomModule") >= 0)) {
          auto customModule =
              JetScapeModuleFactory::createInstance(childElementName);
          if (customModule) {
            Add(customModule);
            JSINFO << "JetScape::DetermineTaskList() -- Hard Process: Added "
                   << childElementName << " to task list.";
          }
        }

        childElement = childElement->NextSiblingElement();
      }

    }

    // Pre-equilibrium
    else if (elementName == "Preequilibrium") {

      tinyxml2::XMLElement *childElement =
          (tinyxml2::XMLElement *)element->FirstChildElement();
      while (childElement) {
        std::string childElementName = childElement->Name();
        VERBOSE(2) << "Parsing childElement: " << childElementName;

        //    - NullPreDynamics
        if (childElementName == "NullPreDynamics") {
          auto predynamics =
              JetScapeModuleFactory::createInstance(childElementName);
          if (predynamics) {
            Add(predynamics);
            JSINFO << "JetScape::DetermineTaskList() -- PreDynamics: Added "
                      "NullPreDynamics to task list.";
          }
        } else if (childElementName == "Glasma") {
          auto predynamics =
              JetScapeModuleFactory::createInstance(childElementName);
          if (predynamics) {
            Add(predynamics);
            JSINFO << "JetScape::DetermineTaskList() -- PreDynamics: Added "
                      "Glasma to task list.";
          }
        } else if (childElementName == "FreestreamMilne") {
        //    - FreestreamMilne
#ifdef USE_FREESTREAM
          auto predynamics =
              JetScapeModuleFactory::createInstance(childElementName);
          if (predynamics) {
            Add(predynamics);
            JSINFO << "JetScape::DetermineTaskList() -- PreDynamics: Added "
                      "FreestreamMilne to task list.";
          }
#else
          JSWARN << "FreestreamMilne is attempted to be added, but freestream "
                    "is not installed!";
#endif
        } else if (((int)childElementName.find("CustomModule") >= 0)) {
        //   - Custom module
          auto customModule =
              JetScapeModuleFactory::createInstance(childElementName);
          if (customModule) {
            Add(customModule);
            JSINFO << "JetScape::DetermineTaskList() -- PreDynamics: Added "
                   << childElementName << " to task list.";
          }
        }

        childElement = childElement->NextSiblingElement();
      }

    }

    // Hydro
    else if (elementName == "Hydro") {

      // First, check if liquefier should be added (Note: Can't use GetXMLElementText(), since that only works for unique tags)
      VERBOSE(2) << "Checking if liquefer should be added: Hydro";
      bool bAddLiquefier = false;
      tinyxml2::XMLElement *childElementLiquefier =
          (tinyxml2::XMLElement *)element->FirstChildElement();
      while (childElementLiquefier) {
        std::string childElementName = childElementLiquefier->Name();
        VERBOSE(2) << "Parsing childElementLiq: " << childElementName;
        if (childElementName == "AddLiquefier") {
          std::string strAddLiquefier = childElementLiquefier->GetText();
          if ((int)strAddLiquefier.find("true") >= 0) {
            bAddLiquefier = true;
            VERBOSE(1) << "Add liquefier to Hydro: True.";
          } else {
            VERBOSE(1) << "Add liquefier to Hydro: False.";
          }
        }
        childElementLiquefier = childElementLiquefier->NextSiblingElement();
      }

      // Loop through elements to look for specific hydro module
      tinyxml2::XMLElement *childElement =
          (tinyxml2::XMLElement *)element->FirstChildElement();
      while (childElement) {
        std::string childElementName = childElement->Name();
        VERBOSE(2) << "Parsing childElement: " << childElementName;

        //   - Brick
        if (childElementName == "Brick") {
          auto hydro = JetScapeModuleFactory::createInstance(childElementName);
          if (hydro) {
            Add(hydro);
            JSINFO << "JetScape::DetermineTaskList() -- Hydro: Added Brick to "
                      "task list.";
            SetModuleId(childElement, hydro);
            if (bAddLiquefier) {
              dynamic_pointer_cast<FluidDynamics>(hydro)->add_a_liquefier(
                  liquefier);
              JSINFO << "JetScape::DetermineTaskList() -- Hydro: Added "
                        "liquefier to Brick.";
            }
          }
        }
        //   - Gubser
        else if (childElementName == "Gubser") {
          auto hydro = JetScapeModuleFactory::createInstance(childElementName);
          if (hydro) {
            Add(hydro);
            JSINFO << "JetScape::DetermineTaskList() -- Hydro: Added Gubser to "
                      "task list.";
            SetModuleId(childElement, hydro);
            if (bAddLiquefier) {
              dynamic_pointer_cast<FluidDynamics>(hydro)->add_a_liquefier(
                  liquefier);
              JSINFO << "JetScape::DetermineTaskList() -- Hydro: Added "
                        "liquefier to Gubser.";
            }
          }
        }
        //   - hydro_from_file
        else if (childElementName == "hydro_from_file") {
          auto hydro = JetScapeModuleFactory::createInstance("HydroFromFile");
          if (hydro) {
            Add(hydro);
            JSINFO << "JetScape::DetermineTaskList() -- Hydro: Added "
                      "hydro_from_file to task list.";
            SetModuleId(childElement, hydro);
            if (bAddLiquefier) {
              dynamic_pointer_cast<FluidDynamics>(hydro)->add_a_liquefier(
                  liquefier);
              JSINFO << "JetScape::DetermineTaskList() -- Hydro: Added "
                        "liquefier to hydro_from_file.";
            }
          }
        }
        //   - MUSIC
        else if (childElementName == "MUSIC") {
#ifdef USE_MUSIC
          auto hydro = JetScapeModuleFactory::createInstance("MUSIC");
          if (hydro) {
            Add(hydro);
            JSINFO << "JetScape::DetermineTaskList() -- Hydro: Added MUSIC to "
                      "task list.";
            SetModuleId(childElement, hydro);
            if (bAddLiquefier) {
              dynamic_pointer_cast<FluidDynamics>(hydro)->add_a_liquefier(
                  liquefier);
              JSINFO << "JetScape::DetermineTaskList() -- Hydro: Added "
                        "liquefier to MUSIC.";
            }
          }
#else
          JSWARN << "MUSIC is attempted to be added, but it is not installed!";
#endif
        }
        //   - CLVisc
        else if (childElementName == "CLVisc") {
          auto hydro = JetScapeModuleFactory::createInstance("CLVisc");
          if (hydro) {
            Add(hydro);
            JSINFO << "JetScape::DetermineTaskList() -- Hydro: Added CLVisc to "
                      "task list.";
          } else {
            JSWARN
                << "CLVisc is attempted to be added, but it is not installed!";
          }
        }

        //   - Custom module
        else if (((int)childElementName.find("CustomModule") >= 0)) {
          auto customModule =
              JetScapeModuleFactory::createInstance(childElementName);
          if (customModule) {
            Add(customModule);
            JSINFO << "JetScape::DetermineTaskList() -- Hydro: Added "
                   << childElementName << " to task list.";
            SetModuleId(childElement, customModule);
            if (bAddLiquefier) {
              dynamic_pointer_cast<FluidDynamics>(customModule)
                  ->add_a_liquefier(liquefier);
              JSINFO << "JetScape::DetermineTaskList() -- Hydro: Added "
                        "liquefier to CustomModule.";
            }
          }
        }

        childElement = childElement->NextSiblingElement();
      }
    }

    // Eloss
    else if (elementName == "Eloss") {

      auto jlossmanager = make_shared<JetEnergyLossManager>();
      auto jloss = make_shared<JetEnergyLoss>();

      // Check if liquefier should be added, and add it if so
      std::string strAddLiquefier =
          GetXMLElementText({"Eloss", "AddLiquefier"});
      if ((int)strAddLiquefier.find("true") >= 0) {
        jloss->add_a_liquefier(liquefier);
        JSINFO << "JetScape::DetermineTaskList() -- Added liquefier to Eloss.";
      } else {
        VERBOSE(1) << "Add liquefier to Eloss: False.";
      }

      // Loop through and add Eloss modules
      tinyxml2::XMLElement *childElement =
          (tinyxml2::XMLElement *)element->FirstChildElement();
      while (childElement) {
        std::string childElementName = childElement->Name();
        VERBOSE(2) << "Parsing childElement: " << childElementName;

        //   - Matter
        if (childElementName == "Matter") {
          auto matter = JetScapeModuleFactory::createInstance(childElementName);
          if (matter) {
            jloss->Add(
                matter); // Note: if you use Matter, it MUST come first (to set virtuality)
            JSINFO << "JetScape::DetermineTaskList() -- Eloss: Added Matter to "
                      "Eloss list.";
          }
        }
        //   - LBT
        else if (childElementName == "Lbt") {
          auto lbt = JetScapeModuleFactory::createInstance(childElementName);
          if (lbt) {
            jloss->Add(
                lbt); // go to 3rd party and ./get_lbtTab before adding this module
            JSINFO << "JetScape::DetermineTaskList() -- Eloss: Added LBT to "
                      "Eloss list.";
          }
        }
        //   - Martini
        else if (childElementName == "Martini") {
          auto martini =
              JetScapeModuleFactory::createInstance(childElementName);
          if (martini) {
            jloss->Add(martini);
            JSINFO << "JetScape::DetermineTaskList() -- Eloss: Added Martini "
                      "to Eloss list.";
          }
        }
        //   - AdS/CFT
        else if (childElementName == "AdSCFT") {
          auto adscft = JetScapeModuleFactory::createInstance(childElementName);
          if (adscft) {
            jloss->Add(adscft);
            JSINFO << "JetScape::DetermineTaskList() -- Eloss: Added AdS/CFT "
                      "to Eloss list.";
          }
        }
        //   - Custom module
        else if (((int)childElementName.find("CustomModule") >= 0)) {
          auto customModule =
              JetScapeModuleFactory::createInstance(childElementName);
          if (customModule) {
            jloss->Add(customModule);
            JSINFO << "JetScape::DetermineTaskList() -- Eloss: Added "
                   << childElementName << " to Eloss list.";
          }
        }

        childElement = childElement->NextSiblingElement();
      }

      jlossmanager->Add(jloss);
      Add(jlossmanager);
    }

    // Jet Hadronization
    else if (elementName == "JetHadronization") {

      // Create hadronization manager and module
      auto hadroMgr = make_shared<HadronizationManager>();
      auto hadro = make_shared<Hadronization>();

      // Determine type of hadronization module, and add it
      std::string hadronizationName =
          element->FirstChildElement("name")->GetText();
      if (hadronizationName == "colored") {
        auto hadroModule =
            JetScapeModuleFactory::createInstance("ColoredHadronization");
        if (hadroModule) {
          hadro->Add(hadroModule);
          JSINFO << "JetScape::DetermineTaskList() -- JetHadronization: Added "
                    "ColoredHadronization to task list.";
        }
      } else if (hadronizationName == "colorless") {
        auto hadroModule =
            JetScapeModuleFactory::createInstance("ColorlessHadronization");
        if (hadroModule) {
          hadro->Add(hadroModule);
          JSINFO << "JetScape::DetermineTaskList() -- JetHadronization: Added "
                    "ColorlessHadronization to task list.";
        }
      } else if (hadronizationName == "hybrid") {
        auto hadroModule =
            JetScapeModuleFactory::createInstance("HybridHadronization");
        if (hadroModule) {
          hadro->Add(hadroModule);
          JSINFO << "JetScape::DetermineTaskList() -- JetHadronization: Added "
                    "HybridHadronization to task list.";
        }
      }
      //   - Custom module
      else if (((int)hadronizationName.find("CustomModule") >= 0)) {
        auto customModule =
            JetScapeModuleFactory::createInstance(hadronizationName);
        if (customModule) {
          hadro->Add(customModule);
          JSINFO << "JetScape::DetermineTaskList() -- JetHadronization: Added "
                 << hadronizationName << " to task list.";
        }
      }

      hadroMgr->Add(hadro);
      Add(hadroMgr);
    }

    // Soft Particlization
    else if (elementName == "SoftParticlization") {

      tinyxml2::XMLElement *childElement =
          (tinyxml2::XMLElement *)element->FirstChildElement();
      while (childElement) {
        std::string childElementName = childElement->Name();
        VERBOSE(2) << "Parsing childElement: " << childElementName;

        //    - iSS
        if (childElementName == "iSS") {
#ifdef iSpectraSampler
          auto iSSmodule =
              JetScapeModuleFactory::createInstance(childElementName);
          if (iSSmodule) {
            tinyxml2::XMLElement *issElement =
                (tinyxml2::XMLElement *)childElement->FirstChildElement();
            while (issElement) {
              std::string issElementName = issElement->Name();
              if( issElementName ==  "Perform_resonance_decays"){
                bool perform_decays_flag = std::stoi(issElement->GetText()) == 1;
                std::dynamic_pointer_cast<SoftParticlization>(iSSmodule)->SetPerformDecays(perform_decays_flag);
              }
              issElement = issElement->NextSiblingElement();
            }
            Add(iSSmodule);
            JSINFO << "JetScape::DetermineTaskList() -- SoftParticlization: "
                      "Added iSS to task list.";
          }
#else
          JSWARN << "iSS is attempted to be added, but iSS is not installed!";
#endif
        }
        //   - Custom module
        else if (((int)childElementName.find("CustomModule") >= 0)) {
          auto customModule =
              JetScapeModuleFactory::createInstance(childElementName);
          if (customModule) {
            Add(customModule);
            JSINFO
                << "JetScape::DetermineTaskList() -- SoftParticlization: Added "
                << childElementName << " to task list.";
          }
        }

        childElement = childElement->NextSiblingElement();
      }
    }

    // Afterburner
    else if (elementName == "Afterburner") {

      tinyxml2::XMLElement *childElement =
          (tinyxml2::XMLElement *)element->FirstChildElement();
      while (childElement) {
        std::string childElementName = childElement->Name();
        VERBOSE(2) << "Parsing childElement: " << childElementName;

        //    - SMASH
        if (childElementName == "SMASH") {
#ifdef USE_SMASH
          auto smashModule =
              JetScapeModuleFactory::createInstance(childElementName);
          if (smashModule) {
            Add(smashModule);
            JSINFO << "JetScape::DetermineTaskList() -- Afterburner: Added "
                      "SMASH to task list.";
          }
#else
          JSWARN
              << "SMASH is attempted to be added, but SMASH is not installed!";
#endif
        }
        //   - Custom module
        else if (((int)childElementName.find("CustomModule") >= 0)) {
          auto customModule =
              JetScapeModuleFactory::createInstance(childElementName);
          if (customModule) {
            Add(customModule);
            JSINFO << "JetScape::DetermineTaskList() -- Afterburner: Added "
                   << childElementName << " to task list.";
          }
        }

        childElement = childElement->NextSiblingElement();
      }
    }

    // Parton printer
    else if (elementName == "PartonPrinter") {

      auto partonPrinter = JetScapeModuleFactory::createInstance(elementName);
      if (partonPrinter) {
        Add(partonPrinter);
        JSINFO << "JetScape::DetermineTaskList() -- Added PartonPrinter to "
                  "task list.";
      }
    }
    else if (elementName == "HadronPrinter") {
      auto hadronPrinter = JetScapeModuleFactory::createInstance(elementName);
      if (hadronPrinter) {
        Add(hadronPrinter);
        JSINFO << "JetScape::DetermineTaskList() -- Added HadronPrinter to "
                  "task list.";
      }
    }

    else {
      VERBOSE(2) << "Nothing to do.";
    }

    element = element->NextSiblingElement();
  }
}

//________________________________________________________________
void JetScape::SetModuleId(tinyxml2::XMLElement *moduleElement,
                           shared_ptr<JetScapeModuleBase> module) {

  tinyxml2::XMLElement *childElement =
      (tinyxml2::XMLElement *)moduleElement->FirstChildElement();
  while (childElement) {
    std::string childElementName = childElement->Name();
    if (childElementName == "name") {
      std::string name = childElement->GetText();
      module->SetId(name);
      JSINFO << "Set ID to: " << name;
    }
    childElement = childElement->NextSiblingElement();
  }
}

//________________________________________________________________
void JetScape::DetermineWritersFromXML() {

  // Get file output name to write to (without file extension, except if custom writer)
  std::string outputFilename = GetXMLElementText({"outputFilename"});

  // Copy string in order to set file extensions for each type
  std::string outputFilenameAscii = outputFilename;
  std::string outputFilenameAsciiGZ = outputFilename;
  std::string outputFilenameHepMC = outputFilename;
  std::string outputFilenameRootHepMC = outputFilename;
  std::string outputFilenameFinalStatePartonsAscii = outputFilename;
  std::string outputFilenameFinalStateHadronsAscii = outputFilename;

  // Check if each writer is enabled, and if so add it to the task list
  CheckForWriterFromXML("JetScapeWriterAscii",
                        outputFilenameAscii.append(".dat"));
  CheckForWriterFromXML("JetScapeWriterAsciiGZ",
                        outputFilenameAsciiGZ.append(".dat.gz"));
  CheckForWriterFromXML("JetScapeWriterHepMC",
                        outputFilenameHepMC.append(".hepmc"));
  CheckForWriterFromXML("JetScapeWriterRootHepMC",
                        outputFilenameRootHepMC.append(".root"));
  CheckForWriterFromXML("JetScapeWriterFinalStatePartonsAscii",
                        outputFilenameFinalStatePartonsAscii.append("_final_state_partons.dat"));
  CheckForWriterFromXML("JetScapeWriterFinalStateHadronsAscii",
                        outputFilenameFinalStateHadronsAscii.append("_final_state_hadrons.dat"));

  // Check for custom writers
  tinyxml2::XMLElement *element =
      (tinyxml2::XMLElement *)JetScapeXML::Instance()
          ->GetXMLRootUser()
          ->FirstChildElement();
  while (element) {
    std::string elementName = element->Name();
    VERBOSE(2) << "Parsing element: " << elementName;

    if (((int)elementName.find("CustomWriter") >= 0)) {
      CheckForWriterFromXML(elementName.c_str(), outputFilename);
    }
    element = element->NextSiblingElement();
  }
}

//________________________________________________________________
void JetScape::CheckForWriterFromXML(const char *writerName,
                                     std::string outputFilename) {

  std::string enableWriter = GetXMLElementText({writerName});
  VERBOSE(2) << "Parsing writer: " << writerName;
  if ((int)enableWriter.find("on") >= 0) {
    VERBOSE(2) << "Writer is on.";
    auto writer = JetScapeModuleFactory::createInstance(writerName);
    if (writer) {
      dynamic_pointer_cast<JetScapeWriter>(writer)->SetOutputFileName(
          outputFilename);
      Add(writer);
      JSINFO << "JetScape::DetermineTaskList() -- " << writerName << " ("
             << outputFilename.c_str() << ") added to task list.";
    }
    // Manually create HepMC writer if it is enabled, since JetScapeModuleFactor::map_type assumes single inheritance
    // from JetScapeModuleBase -- but JetScapeWriterHepMC has multiple inheritance
    else if (strcmp(writerName, "JetScapeWriterHepMC") == 0) {
#ifdef USE_HEPMC
      VERBOSE(2) << "Manually creating JetScapeWriterHepMC (due to multiple "
                    "inheritance)";
      auto writer = std::make_shared<JetScapeWriterHepMC>(outputFilename);
      Add(writer);
      JSINFO << "JetScape::DetermineTaskList() -- " << writerName << " ("
             << outputFilename.c_str() << ") added to task list.";
#else
      JSINFO << "JetScapeWriterHepMC is not enabled. Please, compile with HepMC enabled.";
#endif
    } else if (strcmp(writerName, "JetScapeWriterRootHepMC") == 0) {
#ifdef USE_HEPMC
#ifdef USE_ROOT
      VERBOSE(2) << "Manually creating JetScapeWriterRootHepMC (due to multiple "
                    "inheritance)";
      auto writer = std::make_shared<JetScapeWriterRootHepMC>(outputFilename);
      Add(writer);
      JSINFO << "JetScape::DetermineTaskList() -- " << writerName << " ("
             << outputFilename.c_str() << ") added to task list.";
#else
      JSINFO << "JetScapeWriterRootHepMC is not enabled. Please, compile with ROOT enabled.";
#endif
#else
      JSINFO << "JetScapeWriterRootHepMC is not enabled. Please, compile with HepMC and ROOT enabled.";
#endif
    } else {
      VERBOSE(2) << "Writer is NOT created...";
    }
  } else {
    VERBOSE(2) << "Writer is off.";
  }
}

// JP: Maybe better to move to the SignalManager ...
// kind of cluncky, maybe a better way ... ?
// Handle signal/slots in JetScape hence avoid passing pointers to sub tasks ...
void JetScape::SetPointers() {
  // to get hydro pointer for signals, use signal?
  JSINFO << "Set Hydro,JetEnergylossManager and IS Pointers for "
         << "SignalManager to create Signal/Slots";

  bool hydro_pointer_is_set = false;
  bool bulk_pointer_is_set = false;
  bool iss_pointer_is_set = false;

  for (auto it : GetTaskList()) {
    if (dynamic_pointer_cast<InitialState>(it)) {
      JetScapeSignalManager::Instance()->SetInitialStatePointer(
          dynamic_pointer_cast<InitialState>(it));
    } else if (dynamic_pointer_cast<PreequilibriumDynamics>(it)) {
      JetScapeSignalManager::Instance()->SetPreEquilibriumPointer(
          dynamic_pointer_cast<PreequilibriumDynamics>(it));
    } else if (dynamic_pointer_cast<FluidDynamics>(it) &&
               !hydro_pointer_is_set) {
      JetScapeSignalManager::Instance()->SetHydroPointer(
          dynamic_pointer_cast<FluidDynamics>(it));
      hydro_pointer_is_set = true;
    } else if (dynamic_pointer_cast<BulkDynamicsManager>(it) &&
               !bulk_pointer_is_set) {
      JetScapeSignalManager::Instance()->SetBulkDynamicsManagerPointer(
          dynamic_pointer_cast<BulkDynamicsManager>(it));
      bulk_pointer_is_set = true;
    } else if (dynamic_pointer_cast<JetEnergyLossManager>(it) && !dynamic_pointer_cast<IsrManager>(it)) {
      JetScapeSignalManager::Instance()->SetJetEnergyLossManagerPointer(
          dynamic_pointer_cast<JetEnergyLossManager>(it));
    } else if (dynamic_pointer_cast<HardProcess>(it)) {
      JetScapeSignalManager::Instance()->SetHardProcessPointer(
          dynamic_pointer_cast<HardProcess>(it));
    } else if (dynamic_pointer_cast<JetScapeWriter>(it) && it->GetActive()) {
      JetScapeSignalManager::Instance()->SetWriterPointer(
          dynamic_pointer_cast<JetScapeWriter>(it));
    } else if (dynamic_pointer_cast<PartonPrinter>(it)) {
      JetScapeSignalManager::Instance()->SetPartonPrinterPointer(
          dynamic_pointer_cast<PartonPrinter>(it));
    } else if (dynamic_pointer_cast<SoftParticlization>(it) &&
               !iss_pointer_is_set) {
      JetScapeSignalManager::Instance()->SetSoftParticlizationPointer(
          dynamic_pointer_cast<SoftParticlization>(it));
      iss_pointer_is_set = true;
    } else if (dynamic_pointer_cast<HadronizationManager>(it)) {
      JetScapeSignalManager::Instance()->SetHadronizationManagerPointer(
										dynamic_pointer_cast<HadronizationManager>(it));
    } else if (dynamic_pointer_cast<HadronPrinter>(it)) {
      JetScapeSignalManager::Instance()->SetHadronPrinterPointer(
										dynamic_pointer_cast<HadronPrinter>(it));
				}
  }
}

void JetScape::SetPerEventExecFlags(bool start_of_event)
{
  for (auto it : GetTaskList()) {

    auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);

    if (module && !module->IsTimeStepped() && !std::dynamic_pointer_cast<JetScapeWriter>(it)) {

      if (std::dynamic_pointer_cast<HadronizationManager>(module) || std::dynamic_pointer_cast<Afterburner>(module)) {
        if (start_of_event)
          module->SetActive(false);
        else
          module->SetActive(true);
      }
      else {
        if (start_of_event)
          module->SetActive(taskOrgActiveMap.find(module->GetTaskNumber())->second);
        else
          module->SetActive(false);
      }
   }
  }
}

void JetScape::ResetPerEventExecFlags()
{
  for (auto it : GetTaskList()) {
    auto module = std::dynamic_pointer_cast<JetScapeModuleBase>(it);

    if (module && !module->IsTimeStepped() && !std::dynamic_pointer_cast<JetScapeWriter>(it))
      module->SetActive(taskOrgActiveMap.find(module->GetTaskNumber())->second);
  }
}

void JetScape::Exec() {
  JSINFO << BOLDRED << "Run JetScape ...";
  JSINFO << BOLDRED << "Number of Events = " << GetNumberOfEvents();

  // JetScapeTask::ExecuteTasks(); Has to be called explicitly since not really fully recursively (if ever needed)
  // --> JetScape is "Task Manager" of all modules ...

  // Simple way of passing the writer module pointer
  vector<weak_ptr<JetScapeWriter>> vWriter;

  for (auto it : GetTaskList()) {
    if (dynamic_pointer_cast<JetScapeWriter>(it)) {
      if (it->GetActive()) {
        vWriter.push_back(dynamic_pointer_cast<JetScapeWriter>(it));
      }
    }
    if (dynamic_pointer_cast<JetScapeModuleBase>(it) && it->GetActive()) {
      dynamic_pointer_cast<JetScapeModuleBase>(it)->CheckExec();
    }

    //JP: Maybe a map (unordered) might not be truly necessary, since there is an
    // order, maybe a vector<bool> would suffice, certainly this is more generic
    // using the task number ... (TBD)
    taskOrgActiveMap.emplace(it->GetTaskNumber(), it->GetActive());
  }

  for (int i = 0; i < GetNumberOfEvents(); i++) {
    if (i % n_events_printout == 0) {
      JSINFO << BOLDRED << "Run Event # = " << i;
    }
    VERBOSE(1) << BOLDRED << "Run Event # = " << i;
    JSDEBUG << "Found " << GetNumberOfTasks() << " Modules Execute them ... ";


    // Execute and run per time step for modules if implemented ...
    if (ClockUsed())
    {
      VERBOSE(3)<<"Main Clock Reset ...";

      // Do per event execution except for Hadronization and Afterburner (if not timestepped) ...
      // Set the proper pre per event active etc flags first ...
      SetPerEventExecFlags(true);
      ExecuteTasks();

      GetMainClock()->Reset();

      JetScapeModuleBase::InitPerEventTasks();

      //JP: Quick and dirty to see all tasks ... make recursive if needed
      QueryHistory::Instance()->UpdateTaskMap();

      //QueryHistory::Instance()->PrintTaskMap();
      //JP: Quick and dirty to see all tasks ... make recursive ...
      //QueryHistory::Instance()->PrintTasks();

      do {

        VERBOSE(3)<< BOLDRED << "Current Main Clock Time = "<<GetMainClock()->GetCurrentTime()<<" dT = "<<GetMainClock()->GetDeltaT();

        // quick and dirty here ... (mainly for curiosity ...) and only for CalculateTime assuming that in any case this should take the longest and/or avoiding issues
        // via data exchanges in executive part ... as said, just a first quick look at how things could work out .,,.
        bool multiTask = false;

        //JP: silly to do everything per time step, everything with task map and task vectors could be done in InitPerTimeStep ...
        //   --> to be changed!!!

        if (multiTask) {

          int nTasks = GetNumberOfTasks();
          int nCPUs = thread::hardware_concurrency();

          std::vector<thread> threads;

          int nMaxThreads = nCPUs * 2;
          int n = 0;

          std::vector<std::weak_ptr<JetScapeTask>> vTaskMulti;
          std::vector<std::weak_ptr<JetScapeTask>> vTask;

          VERBOSE(2) << " Use multi-threading: (max) # of threads = # of CPU's "
               << nCPUs << " (found) * 2";

           // JP: also quick and dirty via task map from QueryHistory instance ...
          auto tMap = QueryHistory::Instance()->GetTaskMap();

          for(const auto &x: tMap){
            //if (std::dynamic_pointer_cast<JetScapeModuleBase>(x.second.lock()))
            if (x.second.lock()->GetMultiThread()) {
              //DEBUG
              //cout<<x.second.lock()->GetId()<<endl;
              vTaskMulti.push_back(x.second);
            }
            else
              vTask.push_back(x.second);
          }

          //DEBUG:
          //cout<<vTaskMulti.size()<<" "<<vTask.size()<<endl;

          for (auto t : vTaskMulti)
          {
            if (std::dynamic_pointer_cast<JetScapeModuleBase>(t.lock()))
              threads.push_back(thread(&JetScapeModuleBase::CalculateTime, dynamic_pointer_cast<JetScapeModuleBase>(t.lock())));
          }

          for (auto t : vTask)
          {
            if (std::dynamic_pointer_cast<JetScapeModuleBase>(t.lock()))
              //threads.push_back(thread(&JetScapeModuleBase::CalculateTime, dynamic_pointer_cast<JetScapeModuleBase>(t.lock())));
              std::dynamic_pointer_cast<JetScapeModuleBase>(t.lock())->CalculateTime();
          }

          for (auto &th : threads)
            th.join();

          threads.clear(); vTaskMulti.clear(); vTask.clear();
        }
        else
          JetScapeModuleBase::CalculateTimeTasks();

        JetScapeModuleBase::ExecTimeTasks();

      } while (GetMainClock()->Tick());

      // Follow upo with per event execution of Hadronization and Afterburner (if not timestepped)  ...
      // Set the proper pre per event active etc flags first ...
      SetPerEventExecFlags(false);
      ExecuteTasks();

      // Reset per event flags to orginal state to allow ClearTasks etc to
      // be executed properly and as expected ...
      ResetPerEventExecFlags();
    }
    else
    {
      ExecuteTasks();
      //JP: Quick and dirty to see all tasks ... make recursive if needed
      QueryHistory::Instance()->UpdateTaskMap();
    }

    // Then hand around the collection of writers and ask
    // modules to write what they like
    // Sequence of events:
    // -- writer->Exec is called and redirects to WriteEvent, which starts a new event line
    // -- any remaining exec's finish
    // -- all modules write their headers
    // -- Now all header info is known to the writers, so write out the header
    // -- all other Write()'s are being called
    // the result still confuses me. It's in the best possible order but it shouldn't be.

    // collect module header data
    for (auto w : vWriter) {
      auto f = w.lock();
      if (f) {
        JetScapeTask::CollectHeaders(w);
      }
    }
    // official header
    for (auto w : vWriter) {
      auto f = w.lock();
      if (f) {
        f->WriteHeaderToFile();
      }
    }

    // event data
    for (auto w : vWriter) {
      auto f = w.lock();
      if (f) {
        JetScapeTask::WriteTasks(w);
      }
    }

    // Finalize
    for (auto w : vWriter) {
      auto f = w.lock();
      if (f) {
        f->WriteEvent();
      }
    }

    // JP: If task not active then per time step is active (see above), which could lead to issues with hydro resuse. Follow up!
    // JS: New is timestepped flag should resolve this issue. Anything to undo below?
    // For reusal, deactivate task after it has finished
    // but before it gets cleaned up.
    if (reuse_hydro_) {
      if (n_reuse_hydro_ <= 0) {
        JSWARN << " reuse_hydro is set, but n_reuse_hydro = " << n_reuse_hydro_;
        throw std::runtime_error("Incompatible reusal settings.");
      }
      // Check if iMatter/ISR is used
      bool imatter_is_used = false;
      for (auto it : GetTaskList()) {
        if (it->GetId() == "PythiaGun"){
          for(auto itt : it->GetTaskList()){
            if (itt->GetId() == "IsrManager") {
              VERBOSE(1) << " iMatter is used with reuse_hydro,"
                      << " so initial state is rerun for each event.";
              imatter_is_used = true;
              break;
            }
          }
        }
      }
      bool hydro_pointer_is_set = false;
      bool iss_pointer_is_set = false;
      for (auto it : GetTaskList()) {
        if (!dynamic_pointer_cast<FluidDynamics>(it) &&
            !dynamic_pointer_cast<PreequilibriumDynamics>(it) &&
            !dynamic_pointer_cast<InitialState>(it) &&
            !dynamic_pointer_cast<SoftParticlization>(it)) {
          continue;
        }

        // IS: For ISR+3DGlauber, the initial state 3D Glauber is not used 
        // if imatter is not used, then initial state is rerun
        // This behavior must be rethaught for Au+Au 
        // where we would expect the InitialState to be run per hydro event
        if (imatter_is_used && dynamic_pointer_cast<InitialState>(it)) {
          continue;
        }

        if (dynamic_pointer_cast<FluidDynamics>(it))
          if(dynamic_pointer_cast<FluidDynamics>(it)->IsTimeStepped()) {
            JSWARN << " Reusing hydro with per time stepped = true not allowed!";
            throw std::runtime_error("Reusing hydro with per time stepped = true not allowed.");
          }

        // only deactivate the first hydro
        if (dynamic_pointer_cast<FluidDynamics>(it) && hydro_pointer_is_set) {
          continue;
        }

        if (i % n_reuse_hydro_ == n_reuse_hydro_ - 1) {
          JSDEBUG << " i was " << i
                  << " i%n_reuse_hydro_ = " << i % n_reuse_hydro_
                  << " --> ACTIVATING";
          it->SetActive(true);
          if (dynamic_pointer_cast<FluidDynamics>(it)) {
            hydro_pointer_is_set = true;
          }
        } else {
          JSDEBUG << " i was " << i
                  << " i%n_reuse_hydro_ = " << i % n_reuse_hydro_
                  << " --> DE-ACTIVATING";
          it->SetActive(false);
          if (dynamic_pointer_cast<FluidDynamics>(it)) {
            hydro_pointer_is_set = true;
          }
        }
        // Do the soft hadronization only at once 
        if (dynamic_pointer_cast<SoftParticlization>(it))
          if(dynamic_pointer_cast<SoftParticlization>(it)->IsTimeStepped()) {
            JSWARN << " Reusing hydro with per time stepped = true not allowed!";
            throw std::runtime_error("Reusing hydro with per time stepped = true not allowed.");
          }

        // only deactivate the first iSS
        if (dynamic_pointer_cast<SoftParticlization>(it) && iss_pointer_is_set) {
          continue;
        }

        if (i % n_reuse_hydro_ == n_reuse_hydro_ - 1) {
          JSDEBUG << " i was " << i
                  << " i%n_reuse_hydro_ = " << i % n_reuse_hydro_
                  << " --> ACTIVATING";
          it->SetActive(true);
          if (dynamic_pointer_cast<SoftParticlization>(it)) {
            iss_pointer_is_set = true;
          }
        } else {
          JSDEBUG << " i was " << i
                  << " i%n_reuse_hydro_ = " << i % n_reuse_hydro_
                  << " --> DE-ACTIVATING";
          it->SetActive(false);
          if (dynamic_pointer_cast<SoftParticlization>(it)) {
            iss_pointer_is_set = true;
          }
        }
      }
    }

    // Now clean up, only affects active tasks
    JetScapeModuleBase::ClearTasks();

    //have to call this after writer and call explciitly the clear functions
    //in finish per event, because like writer, clear only for active tasks ...
    //have to think a bit more how to make this workflow more consistent ...
    if (ClockUsed())
      JetScapeModuleBase::FinishPerEventTasks();

    IncrementCurrentEvent();
  }
}

void JetScape::FinishTask() {
  JSINFO << BOLDBLACK << "JetScape finished after " << GetNumberOfEvents()
         << " events!";
  JSDEBUG << "More infos wrap up/saving to file/closing file ...";
}

} // end namespace Jetscape
