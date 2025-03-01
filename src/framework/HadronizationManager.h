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

#ifndef HADRONIZATIONMANAGER_H
#define HADRONIZATIONMANAGER_H

#include "JetScapeTask.h"
#include "JetClass.h"
#include "sigslot.h"

#include "Hadronization.h"

#include <vector>

namespace Jetscape {

//    : public JetScapeTask,
class HadronizationManager
    : public JetScapeModuleBase,
      public std::enable_shared_from_this<HadronizationManager> {

public:
  HadronizationManager();
  virtual ~HadronizationManager();

  // Override Init (not InitTask) here as sub-tasks are called manually
  virtual void Init() override;
  virtual void ClearTask();
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);

  /** Override Exec() (and not ExecuteTask) here as function takes care of calling
   * the subtasks itself with some checks beforehand.
   */
  virtual void ExecuteTask();

  int GetNumSignals();

  void CreateSignalSlots();

  //get Hadrons from Hadronization submodules
  void GetHadrons(vector<shared_ptr<Hadron>>& signal);

  // deletes the hadrons from the different hadronization modules
  // this is used in the case of hadronization hadrons in the afterburner
  // otherwise these hadrons are printed to file and the same hadrons will be 
  // modified in the transport and printed again
  void DeleteHadrons();

  sigslot::signal1<vector<shared_ptr<Hadron>> &> GetHadronList; //get Hadrons from HardProcess NOT Hadronization submodules

  sigslot::signal1<vector<vector<shared_ptr<Parton>>> &> GetFinalPartonList;

  void SetGetFinalPartonListConnected(bool m_GetFinalPartonListConnected) {
    GetFinalPartonListConnected = m_GetFinalPartonListConnected;
  }
  const bool GetGetFinalPartonListConnected() {
    return GetFinalPartonListConnected;
  }

  void SetGetHadronListConnected(bool m_GetHadronListConnected) {
    GetHadronListConnected = m_GetHadronListConnected;
  }
  const bool GetGetHadronListConnected() { return GetHadronListConnected; }


private:
  bool GetFinalPartonListConnected;
  bool GetHadronListConnected;
  vector<vector<shared_ptr<Parton>>> hd;
  vector<shared_ptr<Hadron>> hadrons;
};

} // end namespace Jetscape

#endif
