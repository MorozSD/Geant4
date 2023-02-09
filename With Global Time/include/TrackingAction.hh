#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"
//#include "G4ThreeVector.hh"
#include "globals.hh"

//class PrimaryGeneratorAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackingAction : public G4UserTrackingAction {

  public:  
    TrackingAction();
   ~TrackingAction() override;
   
    virtual void  PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);
     void Abs(G4double dl);
    //void  PreUserTrackingAction(const G4Track*) override;
   // void  PostUserTrackingAction(const G4Track*) override;
 private:
  //EventAction* fEventAction = nullptr;
  G4double fE = 0.;
};

  //G4double fE = 0.;
//  private:
    ///PrimaryGeneratorAction* fPrimary;

   // parameters for generator action #3
  //  G4ThreeVector fNewUz;
//};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

