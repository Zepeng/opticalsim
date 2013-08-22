// Christie Chiu and Ben Jones, MIT, 2012
//
// This is an analyzer module which writes the raw optical
// detector pulses for each PMT to an output file
//


#ifndef OpHitAna_H
#define OpHitAna_H 1

// LArSoft includes
#include "Geometry/Geometry.h"
#include "RawData/OpDetPulse.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"

// ROOT classes.
class TH1D; 
class TTree;

// C++ includes
#include <cstring>
#include <vector>

namespace opdet {
 
  class OpHitAna : public art::EDAnalyzer{
  public:
 
    // Standard constructor and destructor for an ART module.
    OpHitAna(const fhicl::ParameterSet&);
    virtual ~OpHitAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histogram we'll write.
    void beginJob();

    // The analyzer routine, called once per event. 
    void analyze (const art::Event&); 

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fInputModule;              // Input tag for OpDet collection
    float fSampleFreq;                     // in MHz
    float fTimeBegin;                      // in us
    float fTimeEnd;                        // in us
    short fPEheight;                       // in ADC counts

    // Flags to enable or disable output of debugging TH1 / TH2s
    bool fMakeHistPerChannel;
    bool fMakeHistAllChannels;
    bool fMakeHitTree;
 
    // Output TTree and its branch variables
    TTree * fHitTree;
    Int_t   fEventID;
    Int_t   fOpChannel;
    Float_t fPeakTime;
    Float_t fNPe;

    
  };

} 

#endif // OpHitAna_H




// OpHitAna_module.cc

// This is a short program required by the ART framework.  It enables
// a program (OpHitAna, in this case) to be called as a module
// from a .fcl file. It is unlikely that you'll have to make any
// changes to this file.

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet {
  DEFINE_ART_MODULE(OpHitAna);
}


// OpHitAna.cxx

// LArSoft includes
#include "RawData/OpDetPulse.h"
#include "RecoBase/OpHit.h"
#include "OpticalDetector/OpDigiProperties.h"


// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TH1.h"
#include "TTree.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// C++ Includes
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  OpHitAna::OpHitAna(fhicl::ParameterSet const& pset)
  {
    
    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get<std::string>("InputModule");

    art::ServiceHandle<OpDigiProperties> odp;
    fTimeBegin  = odp->TimeBegin();
    fTimeEnd    = odp->TimeEnd();
    fSampleFreq = odp->SampleFreq();


    fMakeHistPerChannel     = pset.get<bool>("MakeHistPerChannel");
    fMakeHistAllChannels    = pset.get<bool>("MakeHistAllChannels");
    fMakeHitTree     = pset.get<bool>("MakeHitTree");

    // If required, make TTree for output
        
    if(fMakeHitTree)
      {
	art::ServiceHandle<art::TFileService> tfs;
	fHitTree = tfs->make<TTree>("HitTree","HitTree");
	fHitTree->Branch("EventID",   &fEventID,   "EventID/I");
	fHitTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
	fHitTree->Branch("PeakTime",  &fPeakTime,  "PeakTime/F"); 
	fHitTree->Branch("NPe",       &fNPe,       "NPe/F");
      }


  }

  //-----------------------------------------------------------------------
  // Destructor
  OpHitAna::~OpHitAna() 
  {}
   
  //-----------------------------------------------------------------------
  void OpHitAna::beginJob()
  {
  }
   

  //-----------------------------------------------------------------------
  void OpHitAna::analyze(const art::Event& evt) 
  {
    
    // Create a handle for our vector of pulses
    art::Handle< std::vector< recob::OpHit > > HitHandle;

    // Create string for histogram name
    char HistName[50];
   

    // Read in HitHandle
    evt.getByLabel(fInputModule, HitHandle);

    // Access ART's TFileService, which will handle creating and writing
    // histograms for us.
    art::ServiceHandle<art::TFileService> tfs;
    
    art::ServiceHandle<geo::Geometry> geom;
    int NOpDet = geom->NOpDet();

    std::vector<TH1D*> HitHist;
    
    sprintf(HistName, "Event %d AllOpDets", evt.id().event());
    
    TH1D * AllHits;
    if(fMakeHistAllChannels)
      {
	AllHits = tfs->make<TH1D>(HistName, ";t (ns);", 
				  int((fTimeEnd - fTimeBegin) * fSampleFreq), 
				  fTimeBegin * 1000., 
				  fTimeEnd * 1000.);
      }
    
    for(int i=0; i!=NOpDet; ++i)
      {
	
	sprintf(HistName, "Event %d OpDet %i", evt.id().event(), i);
	if(fMakeHistPerChannel) HitHist.push_back ( tfs->make<TH1D>(HistName, ";t (ns);", 
							      int((fTimeEnd - fTimeBegin) * fSampleFreq), 
							      fTimeBegin * 1000., 
							      fTimeEnd * 1000.));
	
      }


    fEventID = evt.id().event();


    // For every OpHit in the vector
    for(unsigned int i = 0; i < HitHandle->size(); ++i)
      {
	// Get OpHit
	art::Ptr< recob::OpHit > TheHitPtr(HitHandle, i);
	recob::OpHit TheHit = *TheHitPtr;
	
	fOpChannel = TheHit.OpChannel();
	fNPe       = TheHit.PE();
        fPeakTime  = TheHit.PeakTime();
	
	if(fMakeHitTree)
	  fHitTree->Fill();

	if(fMakeHistPerChannel)
	  {
	    if(fOpChannel>int(HitHist.size()))
	      {
		mf::LogError("OpHitAna")<<"Error : Trying to fill channel out of range, " << fOpChannel;
	      }
	    HitHist[fOpChannel]->Fill(fPeakTime,fNPe);
	  }
	
	if(fMakeHistAllChannels) AllHits->Fill(fPeakTime, fNPe);
	
      }

  }

} // namespace opdet


