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

#include "JetScapeTask.h"
#include "JetScapeTaskSupport.h"
#include "JetScapeLogger.h"

#include "JetEnergyLoss.h"

#include <iostream>

using namespace std;

namespace Jetscape {

JetScapeTask::JetScapeTask() {
  active_exec = true;
  id = "";
  my_task_number_ = JetScapeTaskSupport::Instance()->RegisterTask();
  VERBOSE(9);
}

JetScapeTask::~JetScapeTask() {
  VERBOSE(9);
  JSDEBUG << "Deleting task with id=" << GetId()
          << " and TaskNumber= " << GetMyTaskNumber();
}

void JetScapeTask::Init() {
  InitTask();
  InitTasks();
}

void JetScapeTask::InitTasks() {
  VERBOSE(7) << " : # Subtasks = " << tasks.size();

  for (auto it : tasks) {
    JSDEBUG << "Initalizing " << it->GetId();
    it->Init();
  }
}

void JetScapeTask::Exec() {
  ExecuteTask();
  ExecuteTasks();
}

void JetScapeTask::ExecuteTask() { VERBOSE(7); }


void JetScapeTask::ExecuteTasks() {
  VERBOSE(7) << " : # Subtasks = " << tasks.size();
  for (auto it : tasks) {
    if (it->active_exec) {
      JSDEBUG << "Executing " << it->GetId();
      it->Exec();
	  }
  }
}

void JetScapeTask::Clear() {
  ClearTask();
  ClearTasks();
}

void JetScapeTask::ClearTasks() {
  VERBOSE(7) << " : # Subtasks = " << tasks.size();
  for (auto it : tasks) {
    if (it->active_exec) {
      JSDEBUG << "Clearing " << it->GetId();
      it->Clear();
    }
  }
}

void JetScapeTask::Finish() {
  FinishTask();
  FinishTasks();
}

void JetScapeTask::FinishTasks() {
  VERBOSE(7) << " : # Subtasks = " << tasks.size();
  for (auto it : tasks) {
    if (it->active_exec) {
      JSDEBUG << "Finishing " << it->GetId();
      it->Finish();
    }
  }
}

void JetScapeTask::WriteTasks(weak_ptr<JetScapeWriter> w) {
  //VERBOSE(10);
  if (active_exec) {
    for (auto it : tasks)
      it->WriteTask(w);
  }
}

void JetScapeTask::CollectHeaders(weak_ptr<JetScapeWriter> w) {
  //VERBOSE(10);
  if (active_exec) {
    for (auto it : tasks)
      it->CollectHeader(w);
  }
}

void JetScapeTask::Add(shared_ptr<JetScapeTask> m_tasks) {
  tasks.push_back(m_tasks);
}

} // end namespace Jetscape
