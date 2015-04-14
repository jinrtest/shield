
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "ExN01DetectorConstruction.h"
#include "ExN01PhysicsList.h"
#include "ExN01PrimaryGeneratorAction.h"

#include "convertGeant4ToShield.h"

int main()
{
  // Construct the default run manager
  //
  printf("Test1\n");
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  //
  G4VUserDetectorConstruction* detector = new ExN01DetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new ExN01PhysicsList;
  runManager->SetUserInitialization(physics);

  // set mandatory user action class
  //
  G4VUserPrimaryGeneratorAction* gen_action = new ExN01PrimaryGeneratorAction;
  runManager->SetUserAction(gen_action);
  printf("Test2\n");

//   // Initialize G4 kernel
//   //
  runManager->Initialize();
  geant4toshield::printGeant4ShieldData(convertGeant4ToShield(runManager));
// 
//   // Get the pointer to the UI manager and set verbosities
//   //
//   G4UImanager* UI = G4UImanager::GetUIpointer();
//   UI->ApplyCommand("/run/verbose 1");
//   UI->ApplyCommand("/event/verbose 1");
//   UI->ApplyCommand("/tracking/verbose 1");
// 
//   // Start a run
//   //
//   G4int numberOfEvent = 3;
//   runManager->BeamOn(numberOfEvent);
// 
//   // Job termination
//   //
//   // Free the store: user actions, physics_list and detector_description are
//   //                 owned and deleted by the run manager, so they should not
//   //                 be deleted in the main() program !
//   //
  delete runManager;

  return 0;
}


