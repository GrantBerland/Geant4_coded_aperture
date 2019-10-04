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
// $Id: DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4GenericPolycone.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4IntersectionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 30*cm, env_sizeZ = 30*cm;

    // Material: Vacuum
    //TODO: check pressures, environment for Van Allen belt altitudes
  G4Material* vacuum_material = new G4Material("Vacuum",
              1.0 , 1.01*g/mole, 1.0E-25*g/cm3,
              kStateGas, 2.73*kelvin, 3.0E-18*pascal );
  
  // CZT for detector
  G4Element* Cd = new G4Element("Cadmium","Cd",48., 112.41*g/mole);
  G4Element* Zn = new G4Element("Zinc","Zn", 30., 65.38*g/mole);
  G4Element* Te = new G4Element("Tellurium","Te", 52., 127.60*g/mole);
  G4Material* CZT = new G4Material("CZT", 5.8*g/cm3, 3);
  CZT->AddElement(Cd, 48*perCent);
  CZT->AddElement(Zn, 2*perCent);
  CZT->AddElement(Te, 50*perCent);

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  // G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        vacuum_material,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //
  // Envelope
  //
  G4Box* solidEnv =
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

  G4LogicalVolume* logicEnv =
    new G4LogicalVolume(solidEnv,            //its solid
                        vacuum_material,             //its material
                        "Envelope");         //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  G4AssemblyVolume* detectorAssembly = new G4AssemblyVolume();
  
  G4RotationMatrix Rm;
  G4ThreeVector    Tm;
  G4Transform3D    Tr;

  
  
  G4double boxXY 	   = 4.*cm;
  G4double boxZ  	   = 1.5*mm;
  G4double aperatureSquare = 0.2*cm;
  G4double ap_det_spacing  = 20.*mm;
  G4double detectorXY      = 40.*mm;
  G4double detectorZ       = 5.*mm;
  G4double windowThickness = 0.5*mm;  // 2 windows, each 0.5 mm
  G4double shieldingBoxThickness = 1.*cm;

  // added dimension to "fill the gap" between detectors
  G4Box* aperature_base = new G4Box("Aperature-base",
		   		    (boxXY+2.*mm)/2.,
				    (boxXY+2.*mm)/2.,
				    boxZ/2.);
  
  G4Box* window = new G4Box("Window",
		   	    (boxXY+2.*mm)/2.,
			    (boxXY+2.*mm)/2.,
			    windowThickness/2.);

  
  G4Box* outerShieldingBox = new G4Box("Outer-shielding",
		   	    (boxXY*2.+shieldingBoxThickness)/2.,
			    (boxXY*2.+shieldingBoxThickness)/2.,
			    5.*cm/2.);
  
  
  G4Box* innerSubtractionBox = new G4Box("Inner-sub",
		   	    (boxXY*2.)/2.,
			    (boxXY*2.)/2.,
			    (5.*cm + 1.*cm)/2.);
  
  
  G4RotationMatrix* rotm = new G4RotationMatrix();   

  G4Box* coded_box = new G4Box("Coded-box",
		  		aperatureSquare/2.,
				aperatureSquare/2.,
				boxZ+5.*mm);
  
  
  G4UnionSolid* swapSolid;
  G4String placementXY_str; 
  G4double placementX, placementY; 
  G4String token;
  std::ifstream placementFile("coded_aperture_array.txt", std::ios_base::in);
  
  // Get number of lines in file
  int numberOfBoxes = 0;
  while(getline(placementFile, placementXY_str, '\n'))
    { numberOfBoxes++; }
  
  placementFile.close();


  // Reopen file to start from first line
  placementFile.open("coded_aperture_array.txt", std::ios_base::in);
  getline(placementFile, placementXY_str, '\n');
  
  token = placementXY_str.substr(
		  0, 
  		  placementXY_str.find(',')); 
  
  placementX = std::stod(token);
  
  token = placementXY_str.substr(
		  placementXY_str.find(',')+1, 
		  placementXY_str.find('\n'));
  
  placementY = std::stod(token);
  
  
  G4UnionSolid* coded_boxes = new G4UnionSolid("Combined-boxes",
		  				coded_box,
						coded_box,
						rotm,
						G4ThreeVector(
							placementX*cm,
							placementY*cm,
							0.)); 
  
  // starts at 1 since logicAp1 uses first line of file 
  for(int i=1; i<numberOfBoxes; i++)
  {

    getline(placementFile, placementXY_str, '\n');

    token = placementXY_str.substr(
		  0, 
  		  placementXY_str.find(',')); 
    
    placementX = std::stod(token); 
 
    token = placementXY_str.substr(
		  placementXY_str.find(',')+1, 
		  placementXY_str.find('\n'));
    
    placementY = std::stod(token); 

    swapSolid = new G4UnionSolid("Aperature-base",
	  			   coded_boxes,
	  			   coded_box,
	  			   rotm,
	  			   G4ThreeVector(placementX*cm,
					         placementY*cm,
						 0.));
 
    coded_boxes = swapSolid;
  }
  
  placementFile.close();

  // Subtraction solids
  G4SubtractionSolid* logicAp1 = 
	    new G4SubtractionSolid("Aperature-base",
	  			   aperature_base,
	  			   coded_boxes,
	  			   rotm,
	  			   G4ThreeVector(0.,0.,0.));
  
  G4SubtractionSolid* shieldingBox = 
	    new G4SubtractionSolid("Shielding-Box",
	  			   outerShieldingBox,
	  			   innerSubtractionBox,
	  			   rotm,
	  			   G4ThreeVector(0.,0.,0.));

  // Logical volumes
  G4LogicalVolume* logic_aperature_base =
    new G4LogicalVolume(logicAp1,            //its solid
                        nist->FindOrBuildMaterial("G4_W"), // material
                        "Aperature-base");         //its name

  
  G4LogicalVolume* logic_window =
    new G4LogicalVolume(window,            //its solid
                        nist->FindOrBuildMaterial("G4_Be"), // material
                        "Window");         //its name

  
  G4LogicalVolume* logic_shieldingBox =
    new G4LogicalVolume(shieldingBox,            //its solid
                        nist->FindOrBuildMaterial("G4_W"), // material
                        "Shielding-Box");         //its name

  // Assembly method
  
  // Window 1 (in front of aperture)
  Tm.setX(0.); Tm.setY(0.); Tm.setZ(-(windowThickness+boxZ)/2.);
  Tr = G4Transform3D(Rm, Tm); 
  
  detectorAssembly->AddPlacedVolume(logic_window, Tr);
  
  // Window 2 (between aperture and detector)
  Tm.setX(0.); Tm.setY(0.); Tm.setZ((windowThickness+boxZ)/2.);
  Tr = G4Transform3D(Rm, Tm); 
  
  detectorAssembly->AddPlacedVolume(logic_window, Tr);
 
  // Coded aperture unioned subtraction solid
  Tm.setX(0.); Tm.setY(0.); Tm.setZ(0.);
  Tr = G4Transform3D(Rm, Tm); 
  
  detectorAssembly->AddPlacedVolume(logic_aperature_base, Tr);
  
  
  G4Box* detectorBox = new G4Box("Detector",
		  		    0.5*detectorXY,
				    0.5*detectorXY,
				    0.5*detectorZ);

  G4LogicalVolume* logicDetector = new G4LogicalVolume(detectorBox,
							CZT,
							"Detector");

  Tm.setX(0.); Tm.setY(0.); Tm.setZ(ap_det_spacing);
  Tr = G4Transform3D(Rm, Tm); 
  
  detectorAssembly->AddPlacedVolume(logicDetector, Tr);


  G4int pm1[4] = {1, -1, 1, -1};
  G4int pm2[4] = {1, 1, -1, -1};
  G4double dimX = -2.1*cm;
  G4double dimZ = -2.1*cm;
  
  unsigned int numberDetectors = 4;
  for(unsigned int i=0; i<numberDetectors; i++)
  {
    Tm.setX(pm1[i]*dimX); Tm.setY(pm2[i]*dimZ); Tm.setZ(0.);
    Tr = G4Transform3D(Rm, Tm); 

    // Place assembly in world (or envelope)
    detectorAssembly->MakeImprint(logicEnv, Tr);
  }
  
 
  // Outer shielding box placement
  new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0.,0.,2*cm), 
		      logic_shieldingBox,            //its logical volume
                      "Shielding-Box",               //its name
                      logicEnv,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  // always return the physical World
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
