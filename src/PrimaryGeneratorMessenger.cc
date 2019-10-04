

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIdirectory.hh"


PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* prim)
  :G4UImessenger(),
  fPrimaryGenerator(prim) 
{
  fPrimDir = new G4UIdirectory("/energy/");
  fPrimDir->SetGuidance("Select folding energy and spatial distribution type.");

  fcmd = new G4UIcmdWithAnInteger("/energy/setDistributionType",this);
  fcmd->SetParameterName("0-Near point source, 1-Far away point source, 2-Circular source. ",true);
  fcmd->SetDefaultValue(0);
  fcmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDcmd = new G4UIcmdWithADouble("/energy/setFoldingEnergy",this);
  fDcmd->SetParameterName("Set folding energy in keV",true);
  fDcmd->SetDefaultValue(100.);
  fDcmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fD2cmd = new G4UIcmdWithADouble("/energy/setEventAngle",this);
  fD2cmd->SetParameterName("Set event angular size within (0, 45) degrees",true);
  fD2cmd->SetDefaultValue(45);
  fD2cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fPrimDir;
  delete fcmd;
  delete fDcmd;
  delete fD2cmd;
}


void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, 
					    G4String newValue)
{

  if(command == fcmd){
    fPrimaryGenerator->SetDistType(std::stoi(newValue));
  }    	  

  if(command == fDcmd){
    fPrimaryGenerator->SetFoldingEnergy(std::stod(newValue));
  }
  
  if(command == fD2cmd){
    fPrimaryGenerator->SetEventAngle(std::stod(newValue));
  }

}
