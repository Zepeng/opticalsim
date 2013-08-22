// Christie Chiu and Ben Jones, MIT, 2012
//
// This is an analyzer module which writes the raw optical
// detector pulses for each PMT to an output file
//


#ifndef OpFlashAna_H
#define OpFlashAna_H 1

// LArSoft includes
#include "Geometry/Geometry.h"
#include "RawData/OpDetPulse.h"

#include "TTree.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"

// ROOT classes.
class TH1D;  // 1-dimensional histogram class

// C++ includes
#include <cstring>
#include <vector>

namespace opdet {
 
  class OpFlashAna : public art::EDAnalyzer{
  public:
 
    // Standard constructor and destructor for an ART module.
    OpFlashAna(const fhicl::ParameterSet&);
    virtual ~OpFlashAna();

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
    
    float fYMin, fYMax, fZMin, fZMax;

    int PosHistYRes, PosHistZRes;

    bool fMakeFlashTimeHist;
    bool fMakeFlashPosHist;
    bool fMakePerFlashHists;

    bool fMakePerFlashTree;
    bool fMakePerOpDetTree;
 
    TTree * fPerFlashTree;
    TTree * fPerOpDetTree;
   
    Int_t fEventID;
    Int_t fFlashID;
    Int_t fOpChannel;
    Float_t fFlashTime; 
    Float_t fTotalPe;
    Float_t fNPe;
    Float_t fYCenter;
    Float_t fYWidth;
    Float_t fZCenter;
    Float_t fZWidth;
    
  };

} 

#endif // OpFlashAna_H




// OpFlashAna_module.cc

// This is a short program required by the ART framework.  It enables
// a program (OpFlashAna, in this case) to be called as a module
// from a .fcl file. It is unlikely that you'll have to make any
// changes to this file.

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet {
  DEFINE_ART_MODULE(OpFlashAna);
}


// OpFlashAna.cxx

// LArSoft includes
#include "RawData/OpDetPulse.h"
#include "RecoBase/OpFlash.h"
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
#include "TH2.h"
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
  OpFlashAna::OpFlashAna(fhicl::ParameterSet const& pset)
  {
    
    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get<std::string>("InputModule");


    art::ServiceHandle<OpDigiProperties> odp;
    fTimeBegin  = odp->TimeBegin();
    fTimeEnd    = odp->TimeEnd();
    fSampleFreq = odp->SampleFreq();

    
   
    fYMin =  pset.get<float>("YMin");
    fYMax =  pset.get<float>("YMax");
    fZMin =  pset.get<float>("ZMin");
    fZMax =  pset.get<float>("ZMax");
    
    fMakeFlashTimeHist = pset.get<bool>("MakeFlashTimeHist");
    fMakeFlashPosHist  = pset.get<bool>("MakeFlashPosHist");
    fMakePerFlashHists = pset.get<bool>("MakePerFlashHists");

    fMakePerFlashTree =  pset.get<bool>("MakePerFlashTree");
    fMakePerOpDetTree =  pset.get<bool>("MakePerOpDetTree");

    PosHistYRes=100;
    PosHistZRes=100;

    if(fMakePerOpDetTree)
      {
	fPerOpDetTree = new TTree("PerOpDetTree","PerOpDetTree");
	fPerOpDetTree->Branch("EventID",   &fEventID,    "EventID/I");
	fPerOpDetTree->Branch("FlashID",   &fFlashID,    "FlashID/I");
	fPerOpDetTree->Branch("OpChannel", &fOpChannel,  "OpChannel/I");
	fPerOpDetTree->Branch("FlashTime", &fFlashTime,  "FlashTime/F");
	fPerOpDetTree->Branch("NPe",       &fNPe,        "NPe/F");
      }
    if(fMakePerFlashTree)
      {
	fPerFlashTree = new TTree("PerFlashTree","PerFlashTree");
	fPerFlashTree->Branch("EventID",   &fEventID,    "EventID/I");
	fPerFlashTree->Branch("FlashID",   &fFlashID,    "FlashID/I");
	fPerFlashTree->Branch("YCenter",   &fYCenter,    "YCenter/F");
	fPerFlashTree->Branch("ZCenter",   &fZCenter,    "ZCenter/F");
	fPerFlashTree->Branch("YWidth",    &fYWidth,     "YWidth/F");
	fPerFlashTree->Branch("ZWidth",    &fZWidth,     "ZWidth/F");
	fPerFlashTree->Branch("FlashTime", &fFlashTime,  "FlashTime/F");
	fPerFlashTree->Branch("TotalPe",   &fTotalPe,    "TotalPe/F");
      }
    
    fFlashID=0;
  }

  //-----------------------------------------------------------------------
  // Destructor
  OpFlashAna::~OpFlashAna() 
  {}
   
  //-----------------------------------------------------------------------
  void OpFlashAna::beginJob()
  {
  }
   

  //-----------------------------------------------------------------------
  void OpFlashAna::analyze(const art::Event& evt) 
  {
    
    // Create a handle for our vector of pulses
    art::Handle< std::vector< recob::OpFlash > > FlashHandle;

    // Create string for histogram name
    char HistName[50];
    
    fFlashID=0;

    // Read in HitHandle
    evt.getByLabel(fInputModule, FlashHandle);

    // Access ART's TFileService, which will handle creating and writing
    // histograms for us.
    art::ServiceHandle<art::TFileService> tfs;
    
    std::vector<TH1D*> FlashHist;
    
    fEventID=evt.id().event();

    sprintf(HistName, "Event %d Flash Times", evt.id().event());
    TH1D * FlashTimes;
    if(fMakeFlashTimeHist)
      {
	FlashTimes = tfs->make<TH1D>(HistName, ";t (ns);", 
					    int((fTimeEnd - fTimeBegin) * fSampleFreq), 
					    fTimeBegin * 1000., 
					    fTimeEnd * 1000.);
      }

    TH2D * FlashPositions;
    if(fMakeFlashPosHist)
      {
	sprintf(HistName, "Event %d All Flashes YZ", evt.id().event());
	
        FlashPositions = tfs->make<TH2D>(HistName, ";y ;z ", 
						PosHistYRes, fYMin, fYMax,
						PosHistZRes, fZMin, fZMax);
      }
    
    art::ServiceHandle<geo::Geometry> geom;
    unsigned int NOpDet = geom->NOpDet();
    
    // For every OpFlash in the vector
    for(unsigned int i = 0; i < FlashHandle->size(); ++i)
      {

	// Get OpFlash
	art::Ptr< recob::OpFlash > TheFlashPtr(FlashHandle, i);
	recob::OpFlash TheFlash = *TheFlashPtr;
	


	fFlashTime = TheFlash.Time();
	fTotalPe=0;
	fFlashID++;
	
	TH2D * ThisFlashPosition;
	if(fMakePerFlashHists)
	  {
	    sprintf(HistName, "Event %d t = %f", evt.id().event(), fFlashTime);
	    FlashHist.push_back ( tfs->make<TH1D>(HistName, ";OpChannel;PE", 
						  NOpDet, 0, NOpDet));
	
	    sprintf(HistName, "Event %d Flash %f YZ", evt.id().event(), fFlashTime);
	    
	    ThisFlashPosition = tfs->make<TH2D>(HistName, ";y ;z ", 
						PosHistYRes, fYMin, fYMax,
						PosHistZRes, fZMin, fZMax);
	  }
	fYCenter = TheFlash.YCenter();
	fZCenter = TheFlash.ZCenter();
	fYWidth = TheFlash.YWidth();
	fZWidth = TheFlash.ZWidth();
	
	for(unsigned int j=0; j!=NOpDet; ++j)
	  {
	    if(fMakePerFlashHists) FlashHist.at(FlashHist.size()-1)->Fill(j, TheFlash.PE(j));
	    fNPe = TheFlash.PE(j);
	    fOpChannel=j;
	    
	    if(fMakePerOpDetTree) fPerOpDetTree->Fill();
	    
	    fTotalPe+=fNPe; 
	  }

	for(int y=0; y!=PosHistYRes; ++y)
	  for(int z=0; z!=PosHistZRes; ++z)
	    {
	      float ThisY = fYMin + (fYMax-fYMin)*float(y)/PosHistYRes + 0.0001;
	      float ThisZ = fZMin + (fZMax-fZMin)*float(z)/PosHistZRes + 0.0001;
	      if (fMakePerFlashHists) ThisFlashPosition->Fill(ThisY, ThisZ, fTotalPe * exp(-pow((ThisY-fYCenter)/fYWidth,2)/2.-pow((ThisZ-fZCenter)/fZWidth,2)/2.));
	      if (fMakeFlashPosHist) FlashPositions->Fill(ThisY, ThisZ, fTotalPe * exp(-pow((ThisY-fYCenter)/fYWidth,2)-pow((ThisZ-fZCenter)/fZWidth,2)));
	      
	    }
	
	if(fMakeFlashTimeHist) FlashTimes->Fill(fFlashTime, fTotalPe);
	
	if(fMakePerFlashTree)  fPerFlashTree->Fill();
	
      }

  }

} // namespace opdet


