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
  G4double env_sizeXY = 20.*cm, env_sizeZ = 40.*cm;

    // Material: Vacuum
  G4Material* vacuum_material = new G4Material("Vacuum",
              1.0 , 1.01*g/mole, 1.0E-25*g/cm3,
              kStateGas, 2.73*kelvin, 3.0E-18*pascal );

  // GSFC REF high vacuum at 1e-6 Torr
  G4Material* GSFC_vacuum = new G4Material("Vacuum",
              1.0 , 1.01*g/mole, 1.0E-9*g/cm3,
              kStateGas, 2.73*kelvin, 1.33E-4*pascal );
  
  G4Material* air_material = new G4Material("Air",
              7.0 , 28.97*g/mole, 0.001004*g/cm3,
              kStateGas, 290*kelvin, 82000.*pascal );


  // CZT for detector
  G4Element* Cd = new G4Element("Cadmium","Cd",48., 112.41*g/mole);
  G4Element* Zn = new G4Element("Zinc","Zn", 30., 65.38*g/mole);
  G4Element* Te = new G4Element("Tellurium","Te", 52., 127.60*g/mole);
  G4Material* CZT = new G4Material("CZT", 5.8*g/cm3, 3);
  CZT->AddElement(Cd, 48*perCent);
  CZT->AddElement(Zn, 2*perCent);
  CZT->AddElement(Te, 50*perCent);


  // Nylon 12 for W-nylon layer
  G4Element* C = new G4Element("Carbon", "C", 6., 12.011*g/mole);
  G4Element* H = new G4Element("Hydrogen", "H", 1., 1.008*g/mole);
  G4Element* N = new G4Element("Nitrogen", "N", 7., 14.007*g/mole);
  G4Element* O = new G4Element("Oxygen", "O", 8., 15.999*g/mole);

  G4Material* nylon12   = new G4Material("Nylon12", 1.01*g/cm3, 4);
  G4int natoms;
  nylon12->AddElement(C, natoms=12);
  nylon12->AddElement(H, natoms=23);
  nylon12->AddElement(N, natoms=1);
  nylon12->AddElement(O, natoms=1);

  G4Material* W       = nist->FindOrBuildMaterial("G4_W");
  G4Material* W_nylon = new G4Material("W_nylon", 11.*g/cm3, 2);
  W_nylon->AddMaterial(W, 54.6*perCent);
  W_nylon->AddMaterial(nylon12, 45.4*perCent);


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
                        GSFC_vacuum,           //its material
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
  // FIX ME
  G4double aperatureSquare = 0.2*cm/2.;
  //G4double aperatureSquare = 0.22*cm;
  G4double ap_det_spacing  = 20.*mm;
  G4double detectorXY      = 40.*mm;
  G4double detectorZ       = 5.*mm;
  G4double windowThickness = 2.*mm;

  G4double pixelSize      = 2.5*mm;

  // added dimension to "fill the gap" between detectors
  G4Box* aperature_base = new G4Box("Aperature-base",
		   		    (boxXY+2.*mm)/2.,
				    (boxXY+2.*mm)/2.,
				    boxZ/2.);
  
  G4Box* window = new G4Box("Window",
		   	    (boxXY+2.*mm)/2.,
			    (boxXY+2.*mm)/2.,
			    windowThickness/2.);

  
  G4double collimatorHeight = 17.*mm;
  G4VSolid* collimatorBlock = new G4Box("Collimator",
                                   0.5*collimatorHeight,
                                   0.5*1.*mm,
                                   0.5*2*detectorXY);

  G4RotationMatrix* rotm = new G4RotationMatrix();   

  G4Box* coded_box = new G4Box("Coded-box",
		  		aperatureSquare/2.,
				aperatureSquare/2.,
				boxZ+5.*mm);
  
  
  G4UnionSolid* swapSolid;
  G4String placementXY_str; 
  G4double placementX, placementY; 
  G4String token, filename;

  filename = "NTHT_MURA_array.txt";
  std::ifstream placementFile(filename, std::ios_base::in);
  
  // Get number of lines in file
  int numberOfBoxes = 0;
  while(getline(placementFile, placementXY_str, '\n'))
    { numberOfBoxes++; }

  placementFile.close();


  // Reopen file to start from first line
  placementFile.open(filename, std::ios_base::in);
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
  numberOfBoxes = 0; 
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
  

   G4Box* S1 = new G4Box("S1", 2.5*mm/2, 10*cm/2, 6*cm/2);
   G4Box* S2 = new G4Box("S2", 6.*mm/2,  10*cm/2, 6*cm/2);
   G4Box* S3 = new G4Box("S3", 13.5*mm/2,10*cm/2, 6*cm/2);
   G4Box* S4 = new G4Box("S4", 1.*mm/2,  10*cm/2, 6*cm/2);


   G4LogicalVolume* logicS1 = new G4LogicalVolume(S1,
		   		nist->FindOrBuildMaterial("G4_Sn"),
				"S1");
   G4LogicalVolume* logicS2 = new G4LogicalVolume(S2,
		   		W_nylon,
				"S2");
   G4LogicalVolume* logicS3 = new G4LogicalVolume(S3,
	  		nist->FindOrBuildMaterial("G4_POLYETHYLENE"),
				"S3");
   
   G4LogicalVolume* logicS4 = new G4LogicalVolume(S4,
	  		nist->FindOrBuildMaterial("G4_Al"),
				"S4");

   G4double layerZshift = 1.*cm;
   new G4PVPlacement(0,                     	  //no rotation
                      G4ThreeVector(-4.5*cm,0.,layerZshift), 
		      logicS1,              //its logical volume
                      "S1",               //its name
                      logicEnv,                   //its mother  volume
                      false,                      //no boolean operation
                      0,                          //copy number
                      checkOverlaps);        	  //overlaps checking
   
   new G4PVPlacement(0,                     	  //no rotation
                      G4ThreeVector(-4.5*cm-2.5*mm/2-6.*mm/2,
			      0.,
			      layerZshift), 
		      logicS2,              //its logical volume
                      "S2",               //its name
                      logicEnv,                   //its mother  volume
                      false,                      //no boolean operation
                      0,                          //copy number
                      checkOverlaps);        	  //overlaps checking
   
   new G4PVPlacement(0,                     	  //no rotation
                      G4ThreeVector(-4.5*cm-2.5*mm/2-6.*mm-13.5*mm/2,
			      0.,
			      layerZshift), 
		      logicS3,              //its logical volume
                      "S3",               //its name
                      logicEnv,                   //its mother  volume
                      false,                      //no boolean operation
                      0,                          //copy number
                      checkOverlaps);        	  //overlaps checking

   new G4PVPlacement(0,                     	  //no rotation
                G4ThreeVector(-4.5*cm-2.5*mm/2-6.*mm-13.5*mm-2.*mm/2,
			      0.,
			      layerZshift), 
		      logicS4,              //its logical volume
                      "S4",               //its name
                      logicEnv,                   //its mother  volume
                      false,                      //no boolean operation
                      0,                          //copy number
                      checkOverlaps);        	  //overlaps checking
   // Logical volumes
  G4LogicalVolume* logic_aperature_base =
    new G4LogicalVolume(logicAp1,            //its solid
                        nist->FindOrBuildMaterial("G4_W"), // material
                        "Aperature-base");         //its name

  
  G4LogicalVolume* logic_window =
    new G4LogicalVolume(window,            //its solid
                        nist->FindOrBuildMaterial("G4_Be"), // material
                        "Window");         //its name

  
    
  G4RotationMatrix* rotmCol = new G4RotationMatrix();
  rotmCol->rotateX(90.*deg);

  G4VSolid* collimatorUnion = new G4UnionSolid("Collimator",
                                   collimatorBlock,
                                   collimatorBlock,
                                   rotmCol,
                                   G4ThreeVector());
  rotmCol->rotateX(-90.*deg);
  
  G4LogicalVolume* logicCollimator = new G4LogicalVolume(collimatorUnion,
                                      nist->FindOrBuildMaterial("G4_W"),
                                      "Collimator");
  
  
  
  
  // Assembly method
  
  // Window 1 (in front of aperture)
  Tm.setX(0.); Tm.setY(0.); Tm.setZ(-(windowThickness+boxZ)/2.);
  Tr = G4Transform3D(Rm, Tm); 
  
  detectorAssembly->AddPlacedVolume(logic_window, Tr);
  
  // Coded aperture unioned subtraction solid
  Tm.setX(0.); Tm.setY(0.); Tm.setZ(0.);
  Tr = G4Transform3D(Rm, Tm); 
  
  detectorAssembly->AddPlacedVolume(logic_aperature_base, Tr);
  

  G4VSolid* pixelBlock = new G4Box("Pixel",
                                0.5*pixelSize,
                                0.5*detectorZ,
                                0.5*pixelSize);

  G4AssemblyVolume* pixelAssembly = new G4AssemblyVolume();

  G4RotationMatrix RmT;
  G4ThreeVector    TmT;
  G4Transform3D    TrT;

  G4String pixelName;

  G4int pm1[4] = {1, -1, 1, -1};
  G4int pm2[4] = {1, 1, -1, -1};
  G4double dimX = -2.1*cm;
  G4double dimZ = -2.1*cm;

  G4LogicalVolume* logicPixel;

  // Create pixelated detector
  G4int nameCounter = 0;
  G4int numPixels = 16;
  for(G4int i=0; i<numPixels; i++){
    for(G4int j=0; j<numPixels; j++){

        pixelName = "P";
        pixelName += std::to_string(nameCounter);
        pixelName += "_i";
        pixelName += std::to_string(i);
        pixelName += "_j";
        pixelName += std::to_string(j);
        nameCounter++;

        logicPixel = new G4LogicalVolume(pixelBlock,
                                   CZT,
                                   pixelName);

        TmT.setX(pixelSize*(i-numPixels));
        TmT.setZ(pixelSize*(j-numPixels));
        //TmT.setY(-1.25*cm+0.1*mm);
	TmT.setY(0.);

        TrT = G4Transform3D(RmT, TmT);

        pixelAssembly->AddPlacedVolume(logicPixel, TrT);
    }
  }


  Tm.setX(2.1*cm); Tm.setY(-2.1*cm); Tm.setZ(2.25*cm);
  Rm.rotateX(90.*deg);
  Tr = G4Transform3D(Rm, Tm); 
  Rm.rotateX(-90.*deg);
	  
  // Insert pixel assembly on top of detector
  detectorAssembly->AddPlacedAssembly(pixelAssembly, Tr);	  
  
  unsigned int numberDetectors = 4;
  for(unsigned int i=0; i<numberDetectors; i++)
  {
  	  
    Tm.setX(pm1[i]*dimX); Tm.setY(pm2[i]*dimZ); Tm.setZ(0.);
    Tr = G4Transform3D(Rm, Tm); 

    // Place assembly in world (or envelope)
    detectorAssembly->MakeImprint(logicEnv, Tr);
  }
  
 
  // Collimator
  rotmCol->rotateY(90.*deg);
  new G4PVPlacement(rotmCol,                     //no rotation
                      G4ThreeVector(0.,0.,1.*cm), 
		      logicCollimator,            //its logical volume
                      "Collimator",               //its name
                      logicEnv,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking



  // Linepair test apparatus
  /*
  G4Box* collimatorSourceBlock = new G4Box("Source-box",
		  			3.*cm/2.,
					3.*cm/2.,
					3.*cm/2.);

  G4Box* collimatorSubBlock = new G4Box("Source-sub-box",
		  			2.*cm/2.,
					2.*cm/2.,
					2.*cm/2.);

  G4SubtractionSolid* sourceCollimatorBlock = new G4SubtractionSolid(
		  			"Source-box",
					collimatorSourceBlock,
					collimatorSubBlock,
					new G4RotationMatrix(),
					G4ThreeVector(0,0,1.*cm));


  G4LogicalVolume* logic_collimatorBlock = new G4LogicalVolume(
		  			sourceCollimatorBlock,
					nist->FindOrBuildMaterial("G4_W"),
					"Source-box");


  new G4PVPlacement(0,                     	  //no rotation
                      G4ThreeVector(0.,0.,-50.5*cm), 
		      logic_collimatorBlock,      //its logical volume
                      "Source-box",               //its name
                      logicEnv,                   //its mother  volume
                      false,                      //no boolean operation
                      0,                          //copy number
                      checkOverlaps);        	  //overlaps checking

  */



  G4Box* FaradayCup_box = new G4Box("FC",
		  		    15.*cm/2.,
				    15.*cm/2.,
				    9.65*mm/2.);




  G4LogicalVolume* logicFC = new G4LogicalVolume(FaradayCup_box,
		  			nist->FindOrBuildMaterial("G4_Al"),
					"FC"); 

  
  G4double xShift = 0. * cm;
  G4double yShift = 0. * cm;
  
  new G4PVPlacement(0,                     	  // No rotation
                    G4ThreeVector(xShift,yShift,-20.*cm), 
	      	    logicFC,                 // Logical volume
                    "FC",                     // Name
                    logicEnv,                     // Mother volume
                    false,                        // Boolean operation
                    0,                            // Copy number
                    checkOverlaps);        	  // Overlaps checking
  /*
  // Linepair test has equal line thickness with interline distance
  G4double line_thickness = 4.*mm;
  G4double LP_spacing     = 2.*line_thickness;

  G4Box* LP_box = new G4Box("LP-box",
		  	    5.*cm/2.,
			    5.*cm/2.,
			    1.*cm/2.);

  
  G4Box* LP_sub_line = new G4Box("LP-line",
		  		3.*cm/2.,
				line_thickness/2.,
				3.*cm/2.);
  


  G4UnionSolid* LP_lines = new G4UnionSolid("LP-lines",
		  			LP_sub_line,
					LP_sub_line,
					new G4RotationMatrix(),
					G4ThreeVector(0,LP_spacing,0));  
  
  LP_lines = new G4UnionSolid("LP-lines",
		  		LP_lines,
				LP_sub_line,
				new G4RotationMatrix(),
				G4ThreeVector(0,-LP_spacing,0));

  G4int numLPlines = 2;
  for(G4int i = 2; i < numLPlines; i++)
  {
  
    LP_lines = new G4UnionSolid("LP-lines",
		  		LP_lines,
				LP_sub_line,
				new G4RotationMatrix(),
				G4ThreeVector(0,-i*LP_spacing,0));

    LP_lines = new G4UnionSolid("LP-lines",
		  		LP_lines,
				LP_sub_line,
				new G4RotationMatrix(),
				G4ThreeVector(0,i*LP_spacing,0));


  }
  
  
  G4SubtractionSolid* LP_box_line = new G4SubtractionSolid(
		  			"LP-box",
					LP_box,
					LP_lines,
					new G4RotationMatrix(),
					G4ThreeVector(0,0,0));


  
  G4LogicalVolume* logic_LP_box = new G4LogicalVolume(LP_box_line,
		  			nist->FindOrBuildMaterial("G4_W"),
					"LP-box"); 
  
  G4double xShift = -2. * cm;
  G4double yShift = -2. * cm;
  
  new G4PVPlacement(0,                     	  //no rotation
                      G4ThreeVector(xShift,yShift,-20.*cm), 
		      logic_LP_box,              //its logical volume
                      "LP-box",               //its name
                      logicEnv,                   //its mother  volume
                      false,                      //no boolean operation
                      0,                          //copy number
                      checkOverlaps);        	  //overlaps checking
  
  */ 
  // always return the physical World
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
