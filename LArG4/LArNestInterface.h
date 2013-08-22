#ifndef LARG4_LARNESTINT
#define LARG4_LARNESTINT

//#include "LArVoxelReadout.h"
#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/globals.hh"
#include "TLorentzVector.h"

namespace larg4
{

  // inheritance here is strictly to get DriftIonizationElectrons() method.
  class LArNestInterface //: public LArVoxelReadout
  {
    //    friend void LArVoxelReadout::DriftIonizationElectrons(flni.CurrentStepEnergyDeposit()/MeV, totstep.mag()/cm, midPoint, g4time, trackID );

  public:
    LArNestInterface();

    double CurrentStepEnergyDeposit(); // takes the value from the current G4 step and internally sets the fraction of that energy that goes into ionization and photon production
    int NumberIonizationElectrons() {return fNumIonElectrons;}; // returns the number of ionization electrons created by the energy deposit
    int NumberScintillationPhotons() {return fNumScintPhotons;}; // returns the number of scintillation photons created
    std::string G4Stepper() {return fLastG4Stepper;};
    TLorentzVector    Last4Location() {return fLastStep4Vec;};

    void SetNumIonEl(int &dum) {fNumIonElectrons = dum;};
    void SetNumSciPhi(int &dum) {fNumScintPhotons = dum;};
    void SetEnergyDep(double &dum) { fEnergyDeposit = dum;};
    void SetG4Stepper(std::string &dum) { fLastG4Stepper = dum;};
    void SetLast4Location(TLorentzVector &dum) { fLastStep4Vec = dum;};

    void Reset() {fEnergyDeposit=0.0; fNumIonElectrons=0; fNumScintPhotons=0; fLastG4Stepper="NA"; fLastStep4Vec.SetXYZT(0.,0.,0.,0.);} ;
    // clears the stored value of the energy deposition to be ready for the next step


  private:
    double fEnergyDeposit;
    int fNumIonElectrons;
    int fNumScintPhotons;
    std::string fLastG4Stepper;
    TLorentzVector fLastStep4Vec;

  };

  // fwd-declare
  LArNestInterface& getLArNest();
  //  LArNestInterface& getLArNest(larg4::LArVoxelReadout& l);
  
}

#endif // end LARG4_LARNESTINT
