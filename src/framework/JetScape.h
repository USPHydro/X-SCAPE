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

#ifndef JETSCAPE_H
#define JETSCAPE_H

#include "JetScapeLogger.h"
#include "JetScapeTaskSupport.h"
#include "JetScapeModuleBase.h"
#include "CausalLiquefier.h"
#include <unordered_map>

namespace Jetscape {

class JetScape : public JetScapeModuleBase,
                 public std::enable_shared_from_this<JetScape>
{

public:
  /** Default constructor to create the main task of the JetScape framework.
  */
  JetScape();

  /** This is a destructor for a JetScape.
   */
  virtual ~JetScape();

  /** This function initializes the main task of the JetScape framework.
   * As it calls JetScapeTask::InitTasks() function specifcally to initialize
   * the attached modules/tasks Init (not InitTask) is used here.
  */
  void Init() override;

  /** This function executes the modules/tasks of the main JetScapeTask for all
   * the events. It also calls "GetPartons()" function to print parton shower,
   * and "WriteTasks()" function to store the data in the XML file.
   * It overrides the usual Exec() (and not ExecuteTask()) as it also takes
   * care of calling its subtasks.
   */
  void Exec() override;

  void FinishTask() override;

  /** This function sets the total number of events to "m_n_events".
   */
  void SetNumberOfEvents(int m_n_events) { n_events = m_n_events; }

  /** This function returns the total number of events.
   */
  int GetNumberOfEvents() { return n_events; }

  /** Controls whether to reuse a hydro event (for speedup).
      The number of times is controled by SetNReuseHydro
   */
  inline void SetReuseHydro(const bool reuse_hydro) {
    reuse_hydro_ = reuse_hydro;
  }
  /** Returns whether hydro events are reused.
   */
  inline bool GetReuseHydro() const { return reuse_hydro_; }

  /** Controls number of times a hydro event gets reused.
      Reusal has to be explicitly turned on by SetReuseHydro.
      Turn it on first to avoid a warning.
   */
  inline void SetNReuseHydro(const unsigned int n_reuse_hydro) {
    if (!GetReuseHydro()) {
      JSWARN << "Number of hydro reusals set, but reusal not turned on.";
      JSWARN << "Try jetscape->SetReuseHydro (true);";
    }
    n_reuse_hydro_ = n_reuse_hydro;
  }
  inline unsigned int GetNReuseHydro() const { return n_reuse_hydro_; }

protected:
  void CompareElementsFromXML();
  void recurseToBuild(std::vector<std::string> &elems, tinyxml2::XMLElement *mElement);
  void recurseToSearch(std::vector<std::string> &elems, tinyxml2::XMLElement *uElement);
  void ReadGeneralParametersFromXML();
  void DetermineTaskListFromXML();
  void DetermineWritersFromXML();
  void CheckForWriterFromXML(const char *writerName,
                             std::string outputFilename);
  void SetModuleId(tinyxml2::XMLElement *moduleElement,
                   shared_ptr<JetScapeModuleBase> module);

  void SetPointers();

  /** Function to set the per event execution active flag so that
  if hadronization and Afterburner are attached and not per time step executed,
  that they will be automatically executed after the per time step modules are finished
  So currently possible workflow automatically executed correctly is:
  per event -> per timestep -> per event
   */
  void SetPerEventExecFlags(bool start_of_event);
  /** Function to reset the per event execution active flags to its orginal state
   */
  void ResetPerEventExecFlags();


  void Show();
  int n_events;
  int n_events_printout;

  bool reuse_hydro_;
  unsigned int n_reuse_hydro_;

  std::shared_ptr<CausalLiquefier> liquefier;

 // Option to automatically determine the task list from the XML file,
 // rather than manually calling JetScapeTask::Add() in the run macro.
  bool fEnableAutomaticTaskListDetermination;

  // list to store original SetActive flag settings
  std::unordered_multimap<int , bool > taskOrgActiveMap;

};

} // end namespace Jetscape

#endif
