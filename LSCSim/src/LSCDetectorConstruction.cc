#include "LSCSim/LSCDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RegionStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4String.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"
#include "G4UnionSolid.hh"
#include "GLG4Sim/GLG4param.hh"
#include "LSCSim/LSCPMTSD.hh"
#include "LSCSim/LSC_PMT_LogicalVolume.hh"

using namespace std;

LSCPMTSD * LSCDetectorConstruction::fPmt_SD = NULL;

LSCDetectorConstruction::LSCDetectorConstruction()
  : G4VUserDetectorConstruction(),
    G4UImessenger()
{
  fDetectorDir = new G4UIdirectory("/LSC/det/");

  fGeomCheck = 0;

  fGeomCheckOptCmd = new G4UIcmdWithAnInteger("/LSC/det/geomcheck", this);
  fMaterialDataFileCmd = new G4UIcmdWithAString("/LSC/det/materialdata", this);
  fGeometryDataFileCmd = new G4UIcmdWithAString("/LSC/det/geometrydata", this);
  fPMTPositionDataFileCmd = new G4UIcmdWithAString("/LSC/det/pmtposdata", this);
  fWhichDetectorCmd = new G4UIcmdWithAString("/LSC/det/detector", this);

  fGeometryDataFile = "";
  fPMTPositionDataFile = "";
  fMaterialDataFile = "";
  fWhichDetector = "LS";
}

LSCDetectorConstruction::~LSCDetectorConstruction()
{
  delete fDetectorDir;
  delete fGeomCheckOptCmd;
  delete fMaterialDataFileCmd;
  delete fGeometryDataFileCmd;
  delete fPMTPositionDataFileCmd;
}

void LSCDetectorConstruction::SetNewValue(G4UIcommand * command,
                                          G4String newValues)
{
  if (command == fGeomCheckOptCmd) {
    istringstream is(newValues);
    is >> fGeomCheck;
  }
  if (command == fGeometryDataFileCmd) {
    if (fGeometryDataFile.empty()) fGeometryDataFile = newValues;
  }
  if (command == fMaterialDataFileCmd) {
    if (fMaterialDataFile.empty()) fMaterialDataFile = newValues;
    G4cout << fMaterialDataFile << G4endl;
  }
  if (command == fPMTPositionDataFileCmd) {
    if (fPMTPositionDataFile.empty()) fPMTPositionDataFile = newValues;
  }
  if (command == fWhichDetectorCmd) { fWhichDetector = newValues; }
}

G4VPhysicalVolume * LSCDetectorConstruction::Construct()
{
  ConstructMaterials();
  return ConstructDetector();
}

G4VPhysicalVolume * LSCDetectorConstruction::ConstructDetector()
{
  GLG4param & geom_db(GLG4param::GetDB());

  if (fGeometryDataFile.empty()) {
    G4String msg = "Error, geometery data file could not be opened.\n";
    G4cerr << msg << G4endl;
    G4Exception("LSCDetectorConstruction::LSCDetectorConstruction", "",
                FatalException, msg);
  }
  else {
    geom_db.ReadFile(fGeometryDataFile.c_str());
  }

  // Visualization attributes
  auto visair = new G4VisAttributes();
  visair->SetForceSolid(true);
  visair->SetColour(G4Colour(0.3, 0.3, 0.3, 0.3));

  auto visssteel = new G4VisAttributes();
  visssteel->SetForceSolid(true);
  // visssteel->SetColour(G4Colour::Black());  
  visssteel->SetColour(G4Colour(0.64, 0.7, 0.73, 1));

  auto viswater = new G4VisAttributes();
  viswater->SetForceSolid(true);
  // viswater->SetColour(G4Colour(0,0,1,0.7)); // blue
  viswater->SetColour(G4Colour(0.03, 0.49, 0.69, 0.7));


  // World (Rock)
  G4double worldX = cm * geom_db["worldx"];
  G4double worldY = cm * geom_db["worldy"];
  G4double worldZ = cm * geom_db["worldz"];
  auto WorldBox = new G4Box("WorldBox", worldX / 2, worldY / 2, worldZ / 2);
  auto WorldLog = new G4LogicalVolume(WorldBox, G4Material::GetMaterial("Rock"),
                                      "WorldLog", 0, 0, 0);
  WorldLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto WorldPhys = new G4PVPlacement(0, G4ThreeVector(), WorldLog, "WorldPhys",
                                     0, false, 0, fGeomCheck);

  // Cavern (Top)
  G4double cavernR = cm * geom_db["cavern_arch_radius"];
  G4double cavernH = cm * geom_db["cavern_arch_height"];
  G4double cavernT = cm * geom_db["cavern_arch_thickness"];
  G4double cavernA = deg * geom_db["cavern_arch_angle"];
  auto CavernTubs = new G4Tubs("CavernTubs", cavernR - cavernH / 8, cavernR,
                                cavernH / 2, 0 * deg, cavernA * 2);

  // Cavern (Bottom)
  G4double cavernX = cm * geom_db["cavern_arch_height"];
  G4double cavernY = cm * geom_db["cavern_arch_height"];
  G4double cavernZ = cm * geom_db["cavern_height"];
  auto CavernBox = new G4Box("CavernBox", cavernX / 2, cavernY / 2,
                              (cavernZ - cavernT) / 2.);
  auto CavernRot = new G4RotationMatrix();
  CavernRot->rotateZ(-cavernA);
  CavernRot->rotateY(-90 * deg);

  // Pit
  G4double pitR = cm * geom_db["pit_radius"];
  G4double pitH = cm * geom_db["pit_height"];
  auto PitTubs = new G4Tubs("PitTubs", 0, pitR, pitH / 2, 0, 360 * deg);

  // Cavern (Union of pit, tunnel, top and bottom caverns)
  auto CavernUni = new G4UnionSolid(
      "CavernUni", CavernBox, CavernTubs,
      G4Transform3D(*CavernRot, G4ThreeVector(0, 0,
                                              CavernBox->GetZHalfLength() -
                                                  (cavernR - cavernH / 10))));
  CavernUni = new G4UnionSolid(
      "CavernUni", PitTubs, CavernUni, 0,
      G4ThreeVector(0, 0, pitH / 2. + CavernBox->GetZHalfLength()));

  auto CavernLog = new G4LogicalVolume(
      CavernUni, G4Material::GetMaterial("Air"), "CavernLog");
  CavernLog->SetVisAttributes(visair);
  auto CavernPhys =
      new G4PVPlacement(0, G4ThreeVector(), CavernLog, "CavernPhys", WorldLog,
                        false, 0, fGeomCheck);

  // Veto
  G4VPhysicalVolume * VetoLiquidPhys = nullptr;
  if (fWhichDetector !=  "PROTO") {
    G4double vetoR = cm * geom_db["veto_radius"];
    G4double vetoH = cm * geom_db["veto_height"];
    G4double vetoT = cm * geom_db["veto_thickness"];
    auto VetoTankTubs =
      new G4Tubs("VetoTankTubs", 0, vetoR, vetoH / 2, 0, 360 * deg);
    auto VetoTankLog = new G4LogicalVolume(
					   VetoTankTubs, G4Material::GetMaterial("Steel"), "VetoTankLog", 0, 0, 0);
    // VetoTankLog->SetVisAttributes(G4Colour::White());
    VetoTankLog->SetVisAttributes(visssteel);
    auto VetoTankPhys =
      new G4PVPlacement(0, G4ThreeVector(), VetoTankLog, "VetoTankPhys",
                        CavernLog, false, 0, fGeomCheck);
                        //WorldLog, false, 0, fGeomCheck);

    auto VetoLiquidTubs = new G4Tubs("VetoLiquidTubs", 0, vetoR - vetoT,
				     vetoH / 2 - vetoT, 0, 360 * deg);
    auto VetoLiquidLog =
      new G4LogicalVolume(VetoLiquidTubs, G4Material::GetMaterial("Water"),
                          "VetoLiquidLog", 0, 0, 0);
    // VetoLiquidLog->SetVisAttributes(G4Colour(0, 0, 1, 0.1));
    VetoLiquidLog->SetVisAttributes(viswater);
    VetoLiquidPhys =
      new G4PVPlacement(0, G4ThreeVector(), VetoLiquidLog, "VetoLiquidPhys",
                        VetoTankLog, false, 0, fGeomCheck);

    new G4LogicalBorderSurface("veto_logsurf1", VetoTankPhys, VetoLiquidPhys,
			       Stainless_opsurf);
    new G4LogicalBorderSurface("veto_logsurf2", VetoLiquidPhys, VetoTankPhys,
			       Stainless_opsurf);
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // --- PMT sensitive detector
  ///////////////////////////////////////////////////////////////////////////
  G4SDManager * fSDman = G4SDManager::GetSDMpointer();
  LSCPMTSD * pmtSDInner = new LSCPMTSD("/LSC/PMT/inner");
  fSDman->AddNewDetector(pmtSDInner);

  if (fWhichDetector == "LSCC")
    ConstructDetector_LSC_Cylinder(VetoLiquidPhys, pmtSDInner, geom_db);
  else if (fWhichDetector == "LSCS")
    ConstructDetector_LSC_Sphere(VetoLiquidPhys, pmtSDInner, geom_db);    
  else if (fWhichDetector == "PROTO")
    ConstructDetector_Prototype(WorldPhys, pmtSDInner, geom_db);

  return WorldPhys;
}
