#include <LArNestInterface.h>


namespace larg4
{
  
  LArNestInterface::LArNestInterface() : fNumIonElectrons(0), fNumScintPhotons(0) {}
  double LArNestInterface::CurrentStepEnergyDeposit () 
  {
    return fEnergyDeposit;
  }

 

  // gets and returns singleton. There is only one instance of this class ever.
  LArNestInterface& getLArNest() {
    static LArNestInterface thisLNI;
    return thisLNI;
  }  

  //  LArNestInterface& getLArNest(larg4::LArVoxelReadout& l) {
  //  static LArNestInterface thisLNI(l); // copy constructor
  //  return thisLNI;
  //}  

}
