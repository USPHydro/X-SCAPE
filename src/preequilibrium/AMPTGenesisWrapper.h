#ifndef AMPTGenesis_H
#define AMPTGenesis_H

#include "PreequilibriumDynamics.h"
#include "AMPTGenesis.cpp"

using namespace Jetscape;

class AMPTGenesisWrapper : public PreequilibriumDynamics
{

  public:
    AMPTGenesisWrapper();
    ~AMPTGenesisWrapper();

    void InitializePreequilibrium(PreEquilibriumParameterFile parameter_list);
    void EvolvePreequilibrium();
  private:
  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
    static RegisterJetScapeModule<AMPTGenesisWrapper> reg;
  //int mode; //!< records running mode
    AMPTGenesis *genesis_ptr;




};





#endif // AMPTGenesis_H
