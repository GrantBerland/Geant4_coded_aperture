

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
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

  // Sets fPhotonFilename
  fScmd=new G4UIcmdWithAString("/energy/setPhotonFileName",this);
  fScmd->SetParameterName("Enter file name to draw photons from.",true);
  fScmd->SetDefaultValue("test1_photons.csv");
  fScmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  // Sets fRadioSourceType
  fScmd2=new G4UIcmdWithAnInteger("/energy/setRadioSource",this);
  fScmd2->SetParameterName("Enter isotope name.",true);
  fScmd2->SetDefaultValue(0);
  fScmd2->AvailableForStates(G4State_PreInit, G4State_Idle);

}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fPrimDir;
  delete fcmd;
  delete fDcmd;
  delete fD2cmd;
  delete fScmd;
  delete fScmd2;
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

  if(command == fScmd){
    G4String path = "../analysis/photonFiles/";
    path += newValue;
    fPrimaryGenerator->SetPhotonFilename(path);
  }

  if(command == fScmd2){
    fPrimaryGenerator->SetRadioSourceType(std::stoi(newValue));
  }

}
