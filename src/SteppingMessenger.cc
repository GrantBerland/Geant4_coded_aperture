

#include "SteppingMessenger.hh"

#include "SteppingAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"


SteppingMessenger::SteppingMessenger(SteppingAction* step)
	:G4UImessenger(),fSteppingAction(step) 
{
  fPrimDir = new G4UIdirectory("/dataCollection/");
  fPrimDir->SetGuidance("Set file names for results files.");

  fcmd1 = new G4UIcmdWithAString("/dataCollection/setHitFileName",this);
  fcmd1->SetParameterName("Enter file name, file will appear in ../data/ directory.",true);
  fcmd1->SetDefaultValue("hits.csv");
  fcmd1->AvailableForStates(G4State_PreInit, G4State_Idle);

  fcmd2 = new G4UIcmdWithAString("/dataCollection/setSignalFileName",this);
  fcmd2->SetParameterName("Enter file name.",true);
  fcmd2->SetDefaultValue("signal.csv");
  fcmd2->AvailableForStates(G4State_PreInit, G4State_Idle);


}



SteppingMessenger::~SteppingMessenger()
{
  delete fPrimDir;
  delete fcmd1;
  delete fcmd2;
}


void SteppingMessenger::SetNewValue(G4UIcommand* command, 
					    G4String newValue)
{

  G4String fullFilePath = "../data/";
  if(command == fcmd1){
    fullFilePath += newValue;
    fSteppingAction->SetHitFileName(fullFilePath);
  }    	  

  if(command == fcmd2){
    fullFilePath += newValue;
    fSteppingAction->SetSignalFileName(fullFilePath);
  }

}
