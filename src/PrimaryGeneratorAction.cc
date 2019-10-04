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
// $Id: PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "PrimaryGeneratorMessenger.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include <fstream>
#include <stdexcept>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fPI(3.14159265358979323846),
  fDeg2Rad(3.14159265358979323846/180.),
  lossConeAngleDeg(64.),
  photonPhiLimitDeg(25.),  
  fDistType(0),
  fE0(100.),
  electronParticle(0),
  photonParticle(0),
  fPrimaryMessenger()
{

  fParticleGun  = new G4ParticleGun();

  fPrimaryMessenger = new PrimaryGeneratorMessenger(this);

  electronParticle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  
  photonParticle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");

  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fPrimaryMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Method called to populate member variables with number of 
  // particles to generate
  //CalculateParticlesToGenerate();

  // Selects electron for particle type
  fParticleGun->SetParticleDefinition(photonParticle);

  G4double x, y, z;
  G4double xDir, yDir, zDir;
  G4double theta, phi, R;

  G4double energy;
  G4double narrowingOffset;
  G4double detectorSize;

  detectorSize  = 40*2;  // mm , units assigned in switch-case

  do{
    energy = -(fE0-50) * std::log(1 - G4UniformRand()) * keV;
  } while(energy < 50.*keV);
  
  switch(fDistType)
  {
	case 0: // point source, near
  		x = 10.*mm;
		y = 5.*mm;
  		z = -20.*cm;

  		narrowingOffset = 0.4;
  		xDir = G4UniformRand()*narrowingOffset-narrowingOffset/2.;
  		yDir = G4UniformRand()*narrowingOffset-narrowingOffset/2.;
  		zDir = 1;
		break;

	case 1: // point source, infinitely far
		x = G4UniformRand()*detectorSize - detectorSize/2.; 
		x *= mm;
		y = G4UniformRand()*detectorSize - detectorSize/2.;
		y *= mm;
		z = -20.*cm;

		xDir = yDir = 0;
		zDir = 1;
		break;
	case 2: // structured circle
		R = std::sqrt(G4UniformRand() * 30.) * mm;
		theta = G4UniformRand() * 2. * fPI;
		x = 10.*mm + R * std::cos(theta);
		y = 5.*mm - R * std::sin(theta);
		z = -20.*cm;
		
		xDir = yDir = 0;
		zDir = 1;
		break;

	case 3: // Rotated plane, infinitely far away
		theta = -30. * fDeg2Rad; // pi/18 rad = 10 deg
		
		x = G4UniformRand()*detectorSize - detectorSize/2.; 
		x *= mm;

		x -= 5*cm;

		y = G4UniformRand()*detectorSize - detectorSize/2.;
		y *= mm;
		z = -20.*cm;

		xDir = std::cos(theta) + std::sin(theta); 
		yDir = 0; 
		zDir = -std::sin(theta) + std::cos(theta);
		break;

	case 4: // Uniformly distributed phi angle
		
		R = 20.*cm;

  		// Angle of particle about field line
		theta = G4UniformRand()*2.*fPI;
		
		// Zenith angle of particles wrt to detector
  		phi = std::acos(1 - 2 * G4UniformRand()*
		  (1.-std::cos(photonPhiLimitDeg * fDeg2Rad))/2.);
  		
		// Position on surface of sphere
		x = R * std::sin(phi) * std::cos(theta);
  		y = R * std::sin(phi) * std::sin(theta);
		z = -R * std::cos(phi);  		
	
		// Particle direction 
		xDir = G4UniformRand()*2. - 1.; // |x| ~ U[-1, 1]
  		yDir = G4UniformRand()*2. - 1.; // |y| ~ U[-1, 1]

		// |z| ~ U[0, phi_limit]
		/*
		zDir = std::cos(G4UniformRand()* 
				(photonPhiLimitDeg * fDeg2Rad));
		*/
		zDir = std::sqrt(xDir * xDir + yDir * yDir)/
		std::tan(G4UniformRand()*photonPhiLimitDeg * fDeg2Rad);
		break;

	default:
		throw std::invalid_argument("Choose distribution type!");
  }

  fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
  fParticleGun->SetParticleMomentumDirection(
		  G4ThreeVector(xDir, yDir, zDir));
  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->GeneratePrimaryVertex(anEvent);

}

