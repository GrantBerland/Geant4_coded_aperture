

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class SteppingAction;

class G4UIdirectory;
class G4UIcmdWithAString;

class SteppingMessenger : public G4UImessenger
{
public:
  SteppingMessenger(SteppingAction* );
  virtual ~SteppingMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);


private:
  SteppingAction*            fSteppingAction;
  G4UIdirectory*             fPrimDir;
  G4UIcmdWithAString*        fcmd1;
  G4UIcmdWithAString*        fcmd2;

};

#endif
