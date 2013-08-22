// AnalysisExample_module.cc
// A basic "skeleton" to read in art::Event records from a file,
// access their information, and do something with them. 

// See
// <https://cdcvs.fnal.gov/redmine/projects/larsoftsvn/wiki/Using_the_Framework>
// for a description of the ART classes used here.

// Almost everything you see in the code below may have to be changed
// by you to suit your task. The example task is to make histograms
// and n-tuples related to dE/dx of particle tracks in the detector.

// As you try to understand why things are done a certain way in this
// example ("What's all this stuff about 'auto const&'?"), it will help
// to read ADDITIONAL_NOTES.txt in the same directory as this file.

#ifndef AnalysisExample_Module
#define AnalysisExample_Module

// LArSoft includes
#include "Simulation/SimChannel.h"
#include "Simulation/SimPhotons.h"
#include "Simulation/LArG4Parameters.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "Geometry/Geometry.h"
#include "Geometry/OpDetGeo.h"
#include "RawData/OpDetPulse.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "Utilities/GeometryUtilities.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// C++ Includes
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstring>

namespace AnalysisExample {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class AnalysisExample : public art::EDAnalyzer 
  {
  public:
 
    // Standard constructor and destructor for an ART module.
    explicit AnalysisExample(fhicl::ParameterSet const& pset);
    virtual ~AnalysisExample();

    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun(const art::Run& run);

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event. 
    void analyze (const art::Event& evt); 

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fSimulationProducerLabel; // The name of the producer that tracked simulated particles through the detector
    std::string fHitProducerLabel;        // The name of the producer that created hits
    int fSelectedPDG;                     // PDG code of particle we'll focus on
    double fBinSize;                      // For dE/dx work: the value of dx. 

    // Pointers to the histograms we'll create. 
    //TH1D* fPDGCodeHist;
    TH1D* fEventHist;
    TH1F* fPhotonHist;
    //TH1D* fMomentumHist;
    TH1D* fTrackLengthHist;
    //TH1D* fTrackXHist;
    //TH1D* fTrackYHist;
    //TH1D* fTrackZHist;
    TH1D* fDeviateHist;
    TH3D* fTrackVertexHist;
    TH1D* fAngleHist;
    // The n-tuples we'll create.
    TTree* fSimulationNtuple;
    TTree* fReconstructionNtuple;

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;
    int fPDG;
    int fTrackID;
    // Arrays for 4-vectors: (x,y,z,t) and (Px,Py,Pz,E).
    double fStartXYZT[4];
    double fEndXYZT[4];
    double fStartPE[4];
    double fEndPE[4];
    // Number of dE/dx bins.
    int fNdEdxBins;
    // Why "double*"? Because this is going to be an array whose length
    // (fNdEdxBins) we don't yet know.
    double* fdEdxBins;

    // Other variables that will be shared between different methods.
    art::ServiceHandle<geo::Geometry> fGeometry;       // pointer to Geometry service

  }; // class AnalysisExample


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // Constructor
  AnalysisExample::AnalysisExample(fhicl::ParameterSet const& parameterSet)
  {
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  // Destructor
  AnalysisExample::~AnalysisExample() 
  {}
   
  //-----------------------------------------------------------------------
  void AnalysisExample::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    fPhotonHist = tfs->make<TH1F>("photons",";photon times;", 1500000, 0, 1500000);
  }
   
  //-----------------------------------------------------------------------
  void AnalysisExample::beginRun(const art::Run& run)
  {
    // How to convert from number of electrons to GeV. If we're
    // getting this from a data file, The ultimate source of this
    // conversion factor is
    // ${SRT_PUBLIC_CONTEXT}/SimpleTypesAndConstants/PhysicalConstants.h.
    art::ServiceHandle<sim::LArG4Parameters> larParameters;
  }

  //-----------------------------------------------------------------------
  void AnalysisExample::reconfigure(fhicl::ParameterSet const& p)
  {
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fSimulationProducerLabel = p.get< std::string >("SimulationLabel");
    fHitProducerLabel        = p.get< std::string >("HitLabel");
    fSelectedPDG             = p.get< int         >("PDGcode");
    fBinSize                 = p.get< double      >("BinSize");
    return;
  }

  //-----------------------------------------------------------------------
  void AnalysisExample::analyze(const art::Event& event) 
  {

    /*
    TVector3 fTPCCenter = fGeometry->GetTPCFrontFaceCenter();
    std::cout << fTPCCenter.X() << std::endl;
    std::cout << fTPCCenter.Y() << std::endl;
    std::cout << fTPCCenter.Z() << std::endl;
    TVector3 fTPCDimensions = TVector3(fGeometry->DetHalfWidth(), fGeometry->DetHalfHeight(), fGeometry->DetLength());
    std::cout << fTPCDimensions.X() << std::endl;
    std::cout << fTPCDimensions.Y() << std::endl;
    std::cout << fTPCDimensions.Z() << std::endl;*/
    art::ServiceHandle<util::LArProperties> LArProp;
    art::ServiceHandle<util::DetectorProperties> detprop;

    
    double OpDetCenter[3];
    unsigned int o=0; unsigned int c=0;

    std::cout << fGeometry->Ncryostats() << std::endl;
    for(int OpChan=0; OpChan < fGeometry->NOpChannels(); OpChan++)
    {
     fGeometry->OpChannelToCryoOpDet(OpChan,o,c);
     fGeometry->Cryostat(c).OpDet(o).GetCenter(OpDetCenter);
     std::cout << OpDetCenter[0] << "   " << OpDetCenter[1] <<  "   " << OpDetCenter[2] << std::endl;
    }


    art::Handle< std::vector<sim::LBNE10ktPhotons> > photonHandle;
    event.getByLabel("largeant", photonHandle);

    std::cout << (*photonHandle).size() << std::endl;

    for(auto const& photon : (*photonHandle) )
    {
      std::cout << "OpChannel:" << photon.OpChannel << std::endl;
      std::map<int, int> DetMap = photon.DetectedPhotons;
      for(auto it = DetMap.begin(); it!= DetMap.end(); it++)
      {
        fPhotonHist->Fill(it->first, it->second);
        //std::cout << it->first << "  " << it->second << std::endl;
      }

    }
    
    art::Handle< std::vector<raw::OpDetPulse> > WaveformHandle;
    event.getByLabel("opdigi",WaveformHandle);
    for(unsigned int i = 0; i < WaveformHandle->size(); ++i)
    {
      art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle, i);
      raw::OpDetPulse ThePulse(*ThePulsePtr);

      auto onepulse = ThePulse.Waveform();
      //std::cout << onepulse.size() << std::endl;
      //std::cout << (pulse.Waveform()).size() << std::endl;
      for(unsigned int binNum = 0; binNum < onepulse.size(); ++binNum)
      {
        if(onepulse.at(binNum) > 0) std::cout << onepulse.at(binNum) << std::endl;
      }    
      std::cout << ThePulse.OpChannel() << std::endl;
    }
   
    return;
  }

  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see AnalysisExample.fcl for more information.
  DEFINE_ART_MODULE(AnalysisExample);

} // namespace AnalysisExample

#endif // AnalysisExample_Module

