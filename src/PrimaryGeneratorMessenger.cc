

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

  fcmd2 = new G4UIcmdWithAnInteger("/energy/setEnergyDistributionType",this);
  fcmd2->SetParameterName("0-Monoenergetic, 1-Exponential",true);
  fcmd2->SetDefaultValue(0);
  fcmd2->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fDcmd = new G4UIcmdWithADouble("/energy/setFoldingEnergy",this);
  fDcmd->SetParameterName("Set folding energy in keV",true);
  fDcmd->SetDefaultValue(100.);
  fDcmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDcmd2 = new G4UIcmdWithADouble("/energy/setSourceDistance",this);
  fDcmd2->SetParameterName("Set source distance in cm",true);
  fDcmd2->SetDefaultValue(20.);
  fDcmd2->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fPrimDir;
  delete fcmd;
  delete fcmd2;
  delete fDcmd;
  delete fDcmd2;
}


void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, 
					    G4String newValue)
{
  if(command == fcmd){
    fPrimaryGenerator->SetDistType(std::stoi(newValue));
  }    	  

  if(command == fcmd2){
    fPrimaryGenerator->SetEnergyDistribution(std::stoi(newValue));
  }    	  
  
  if(command == fDcmd){
    fPrimaryGenerator->SetFoldingEnergy(std::stod(newValue));
  }

  if(command == fDcmd2){
    fPrimaryGenerator->SetSourceDistance(std::stod(newValue));
  }
}
