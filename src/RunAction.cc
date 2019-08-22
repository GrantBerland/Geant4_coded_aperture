//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: RunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
// #include "Run.hh"
// #include "DetectorAnalysis.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
// #include "HistoManager.hh"

#include "Randomize.hh"
#include <ctime>

//#include <fstream>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


RunAction::RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{

  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
/*
  long seeds[2];
  time_t systime = time(NULL);
  
  seeds[0] = (long) systime;
  seeds[1] = (long) (systime*G4UniformRand());
  G4Random::setTheSeeds(seeds);
  std::cout << "Seeds set now: " << seeds[0] << ", " << seeds[1] << std::endl;
*/
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
/*
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

// Compute dose = total energy deposit in a run and its variance
G4double edep  = fEdep.GetValue();
G4double edep2 = fEdep2.GetValue();

G4double rms = edep2 - edep*edep/nofEvents;
if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;
*/
  // Print (or write to file)
  //

     //G4AnalysisManager* man = G4AnalysisManager::Instance();
     //man->Write();
     //man->CloseFile();
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEdep(G4double edep)
{
  fEdep += edep;
  fEdep2 += edep*edep;
}

void RunAction::LogEntry(G4double)
{
/*
  (*asciiFile) << std::setiosflags(std::ios::fixed)
   << std::setprecision(3)
   << std::setiosflags(std::ios::right)
   << std::setw(10);
  (*asciiFile) << edep << G4endl;
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
