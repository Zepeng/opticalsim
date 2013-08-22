////////////////////////////////////////////////////////////////////////
/// \file OpDetPhotonTable.h
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// 
// This class holds a collection of PMT hits to be stored
// into an event as produced by the fast scintillation process.
//
// When a scintillating particle is stepped, it may generate
// one or more detected photons in each optical detector, 
// depending upon its location.
//
// For the fast scintillation process, the likelihood of generating
// a detected photo given a position in the detector is looked up.
// If one is detected, a hit in the relevant optical detector at
// an appropriate time (with ~20ns error bars for flight time) is 
// stored in this table, which is eventually read out at the end
// of the event by LArG4_module and stored in the event. 
//
// For slow scintillation / cerenkov processes, photons are
// generated and stepped about the detector, and if one steps
// into a volume designated at sensitive by Geant4, a simphoton
// is added to this table.
//
// The two sources can be distinguished by looking at the
// SetInSD flag of the OnePhoton object.
//
// Ben Jones, MIT, 11/10/12
//

#include "Geant4/G4PhysicalVolumeStore.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include <map>

#ifndef LBNE10KTOPDETPHOTONTABLE_h
#define LBNE10KTOPDETPHOTONTABLE_h 1


namespace larg4 {
  class LBNE10ktOpDetPhotonTable
    {
    public:
      ~LBNE10ktOpDetPhotonTable(){}
      static LBNE10ktOpDetPhotonTable * Instance();
   
      void AddPhoton( std::map<int, std::map<int, int>>* StepPhoton);
      
      std::map<int, std::map<int, int> >   GetPhotons();
      std::map<int, int>&                  GetPhotonsForOpChannel(int opcannel);
      
      void ClearTable();
      
    protected:
      LBNE10ktOpDetPhotonTable();

    private:
      std::map<int, std::map<int,int> > fDetectedPhotons;
      
    };

}


#endif
