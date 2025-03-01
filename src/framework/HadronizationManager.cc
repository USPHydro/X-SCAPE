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

#include "HadronizationManager.h"
#include "JetScapeLogger.h"
#include "JetScapeSignalManager.h"
#include <string>
#include "Hadronization.h"

#include <iostream>
#include <vector>

#include <thread>
#include <future>

using namespace std;

namespace Jetscape {

HadronizationManager::HadronizationManager() {
  SetId("HadronizationManager");
  GetFinalPartonListConnected = false;
  GetHadronListConnected = false;
  VERBOSE(8);
}

HadronizationManager::~HadronizationManager() {
  ClearTask();

  // The tasks are hadronization modules
  if (GetNumberOfTasks() > 0)
    EraseTaskLast();
}

void HadronizationManager::ClearTask() {
  JSDEBUG << "Hadronization task List ...";

  hd.clear();

  int n = GetNumberOfTasks();
  for (int i = 1; i < n; i++)
    EraseTaskLast();

  JetScapeSignalManager::Instance()->CleanUp();
  VERBOSE(8) << hd.size();
}

void HadronizationManager::Init() {
  JSINFO << "Initialize Hadronization Manager ...";

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid Hadronization Manager modules found ...";
    exit(-1);
  }

  JSINFO << "Found " << GetNumberOfTasks()
         << " Hadronization Manager Tasks/Modules Initialize them ... ";
  JetScapeTask::InitTasks();

  JSINFO << "Connect HadronizationManager Signal to Energy Loss Manager ...";
  JetScapeSignalManager::Instance()->ConnectGetFinalPartonListSignal(shared_from_this());
}

void HadronizationManager::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  JetScapeTask::WriteTasks(w);
}

void HadronizationManager::ExecuteTask() {

  VERBOSE(2) << "Run Hadronization Manager ...";
  JSDEBUG << "Task Id = " << this_thread::get_id();

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid Hadronization modules found ...";
    exit(-1);
  }

  CreateSignalSlots();

  if (GetGetFinalPartonListConnected()) {
    GetFinalPartonList(hd);
    hadrons.clear();
    GetHadronList(hadrons);
    VERBOSE(2) << " There are " << hd.size()
               << " partons ready for hadronization";
    VERBOSE(2) << " There are already " << hadrons.size() << " hadrons";
   
    
    for (auto it : GetTaskList()) {
      dynamic_pointer_cast<Hadronization>(it)->AddInPartons(hd);
      dynamic_pointer_cast<Hadronization>(it)->AddInHadrons(hadrons);
    }
  } else {
    VERBOSE(2) << " There are no partons ready for recombination";
  }
}

void HadronizationManager::CreateSignalSlots() {
  for (auto it : GetTaskList()) {
    for (auto it2 : it->GetTaskList()) {
      if (!dynamic_pointer_cast<Hadronization>(it2)
               ->GetTransformPartonsConnected()) {
        JetScapeSignalManager::Instance()->ConnectTransformPartonsSignal(
            dynamic_pointer_cast<Hadronization>(it),
            dynamic_pointer_cast<Hadronization>(it2));
      }
      if (!dynamic_pointer_cast<Hadronization>(it2)
               ->GetGetHydroHyperSurfaceConnected()) {
        JetScapeSignalManager::Instance()->ConnectGetHydroHyperSurfaceSignal(
            dynamic_pointer_cast<Hadronization>(it2));
      }
      if (!dynamic_pointer_cast<Hadronization>(it2)
               ->GetGetHydroCellSignalConnected()) {
        JetScapeSignalManager::Instance()->ConnectGetHydroCellSignal(
            dynamic_pointer_cast<Hadronization>(it2));
      }
    }
  }

  JetScapeSignalManager::Instance()->PrintTransformPartonsSignalMap();
}

void HadronizationManager::GetHadrons(vector<shared_ptr<Hadron>>& signal){
	//signal = outHadrons;
	signal.clear();
	// foreach hadronizon object tasks
	for(shared_ptr<JetScapeTask> it : GetTaskList()){
		vector<shared_ptr<Hadron>> tempHadronList;
		JetScapeTask *jet = it.get();
		Hadronization *hit = (Hadronization *) jet;
		tempHadronList = hit->GetHadrons();
		for(auto hadron : tempHadronList){
			signal.push_back(hadron);
		}
	}
}

void HadronizationManager::DeleteHadrons() {
  // foreach hadronizon object tasks
	for(shared_ptr<JetScapeTask> it : GetTaskList()){
		JetScapeTask *jet = it.get();
		Hadronization *hit = (Hadronization *) jet;
		hit->DeleteHadrons();
	}
}

} // namespace Jetscape
