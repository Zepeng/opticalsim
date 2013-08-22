// Christie Chiu and Ben Jones, MIT, 2012
//
// This module finds large pulses in an optical detector
// output waveform and produces OpHit objects representing
// discrete N.PE pulses. 
//  
// These OpHits are stored into the event for further
// analysis.
//

// OpHitFinder.h

#ifndef OpHitFinder_H
#define OpHitFinder_H 1

// LArSoft includes
#include "Geometry/Geometry.h"
#include "RawData/OpDetPulse.h"

// ART includes.
#include "art/Framework/Core/EDProducer.h"

// ROOT classes.
class TH1D;  // 1-dimensional histogram class

// C++ includes
#include <cstring>
#include <vector>

namespace opdet {
 
  class OpHitFinder : public art::EDProducer{
  public:
 
    // Standard constructor and destructor for an ART module.
    OpHitFinder(const fhicl::ParameterSet&);
    virtual ~OpHitFinder();

    // This method is called once, at the start of the job. In this
    // example, it will define the histogram we'll write.
    void beginJob();

    // The producer routine, called once per event. 
    void produce (art::Event&); 

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fInputModule;              // Input tag for OpDet collection
    float fSampleFreq;                     // in MHz
    float fPEArea;                         // in ADC counts
    float fHitThreshold;
    float fIsolationFrac;
    float fTimeBegin;
    float fTimeEnd;
    int   fSamplesPerBin;
    float fZeroSupThresh;
  };

} 

#endif // OpHitFinder_H




// OpHitFinder_module.cc

// This is a short program required by the ART framework.  It enables
// a program (OpHitFinder, in this case) to be called as a module
// from a .fcl file. It is unlikely that you'll have to make any
// changes to this file.

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet {
  DEFINE_ART_MODULE(OpHitFinder);
}


// OpHitFinder.cxx

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
#include "TF1.h"
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
  OpHitFinder::OpHitFinder(fhicl::ParameterSet const& pset)
  {
    // Infrastructure piece
    produces<std::vector< recob::OpHit> >();

    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get<std::string>("InputModule");

    // Get binning parameters from .fcl
    fSamplesPerBin = pset.get<int>("SamplesPerBin");
    fIsolationFrac = pset.get<float>("IsolationFrac");
    fHitThreshold  = pset.get<float>("HitThreshold");
    fZeroSupThresh = pset.get<float>("ZeroSupThresh");
    fPEArea        = pset.get<float>("PEArea");

    art::ServiceHandle<OpDigiProperties> odp;
    fTimeBegin  = odp->TimeBegin();
    fTimeEnd    = odp->TimeEnd();
    fSampleFreq = odp->SampleFreq();

  }

  //-----------------------------------------------------------------------
  // Destructor
  OpHitFinder::~OpHitFinder() 
  {}
   
  //-----------------------------------------------------------------------
  void OpHitFinder::beginJob()
  {
  }
   

  //-----------------------------------------------------------------------
  void OpHitFinder::produce(art::Event& evt) 
  {
    // Infrastructure piece
    std::unique_ptr<std::vector< recob::OpHit > > StoragePtr (new std::vector<recob::OpHit>);
    
    // Create a handle for our vector of pulses
    art::Handle< std::vector< raw::OpDetPulse > > WaveformHandle;
    
    
    // Create string for histogram name
    //  char HistName[50];
    //  char FitName[50];
    //  char FitEq[200];

    // Read in WaveformHandle
    evt.getByLabel(fInputModule, WaveformHandle);


    std::vector<int> Binned1;
    std::vector<int> Binned2;
    std::vector<int> Unipolar;
    
    // For every OpDet waveform in the vector given by WaveformHandle
    for(unsigned int i = 0; i < WaveformHandle->size(); ++i)
      {
	art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle, i);
	raw::OpDetPulse ThePulse(*ThePulsePtr);
	
	unsigned int Samples = ThePulse.Waveform().size();
	
	Binned1.clear();
	Binned2.clear();
	Unipolar.clear();

	Binned1.resize(int( Samples / fSamplesPerBin + 1)); 
	Binned2.resize(int( Samples / fSamplesPerBin + 1)); 
	Unipolar.resize(Samples); 
		       
	

	for(unsigned int binNum = 0; binNum < ThePulse.Waveform().size(); ++binNum)
	  {
	    // Fill histogram with waveform
	    // Raw signal is positive, so to prevent wraparound we set very negative 
	    //   values to unsigned short (this is probably not the best fix)

	    Unipolar[binNum]          = (ThePulse.Waveform()[binNum]);
	    if(binNum>1) Unipolar[binNum] += Unipolar[binNum-1];
	    if(Unipolar[binNum] < fZeroSupThresh) Unipolar[binNum]=0;


	    Binned1[int(binNum / fSamplesPerBin)          ] += Unipolar[binNum];
	    Binned2[int((binNum / fSamplesPerBin) + 0.5)  ] += Unipolar[binNum];
	    
	  }

	std::vector<int> FoundHits1;
	std::vector<int> FoundHits2;

	for(size_t binNum=0; binNum!=Binned1.size(); ++binNum)
	  {
	    if((binNum>0)&&((binNum+1)<Samples))
	      {

		if( Binned1[binNum]!=0)
		  {
		    if( (Binned1[binNum-1]==0) || (Binned1[binNum-1]/Binned1[binNum] < fIsolationFrac))
		      
		      if( (Binned1[binNum+1]==0)|| (Binned1[binNum+1]/Binned1[binNum] < fIsolationFrac))
			if(Binned2[binNum] > (fHitThreshold*fPEArea)) FoundHits1.push_back(binNum);
		  }
		if(Binned2[binNum]!=0)
		  {
		    if( (Binned2[binNum-1]==0) || (Binned2[binNum-1]/Binned2[binNum] < fIsolationFrac))
		      if( (Binned2[binNum+1]==0) || (Binned2[binNum+1]/Binned2[binNum] < fIsolationFrac))
			if(Binned2[binNum] > (fHitThreshold*fPEArea)) FoundHits2.push_back(binNum);		
		  }
	      }
	  }
	
	// we store hits in a map like this, indexed by
	//  peak time.  Therefore the same hit found by
	//  both binning schemes only counts one.

	std::map<int, recob::OpHit> HitsByTime;
	

	// Loop through making the hit objects
	for(size_t hit=0; hit!=FoundHits1.size(); ++hit)
	  {
	    //  mf::LogInfo("OpHitFinder")<<"Found hit in hist 1: " << FoundHits1[hit]<<" " <<Binned1[FoundHits1[hit]]<<std::endl;
	    unsigned int LowBin  = FoundHits1.at(hit ) * fSamplesPerBin;
	    double Area          =  0;
	    double Amplitude     =  0;
	    double PeakTime      =  0;
	    double Width         =  0;
	    double PE            =  0;

	    double PeakTimeError =  0;
	    double AreaError     =  0;
	    double AmplitudeError=  0;
	    double WidthError    =  0;
	    double PEError       =  0;

	    double HitVar        =  0;
	    
	    double sumQ2=0, sumQ=0;
	    int PeakBin = LowBin;
	    
	    for(size_t binNum = LowBin; binNum!=Unipolar.size(); ++binNum)
	      {
		double ThisVal = Unipolar[binNum];
		sumQ  += ThisVal*binNum;
		sumQ2 += pow(ThisVal,2)*binNum;
		Area  += ThisVal;
		if(ThisVal>Amplitude)
		  {
		    Amplitude = ThisVal;
		    PeakBin   = binNum; 
		  }	     
		if( int(binNum) - int(PeakBin) > fSamplesPerBin) break;
		
	      }
	    
	    //	    HitMean = sumQ / float(Area);
	    HitVar  = pow(sumQ2 - pow(sumQ,2),0.5)/float(Area);
	   
	    PeakTime      = fTimeBegin*1000. + PeakBin *                1000. / fSampleFreq;
	    PeakTimeError = HitVar / pow(Area,0.5) * 1000. / fSampleFreq;
	    Width         = HitVar *                 1000. / fSampleFreq;
	    PE            = Area / fPEArea;
	    
	    if(Area>0)
	      {
		HitsByTime[PeakBin] = recob::OpHit(ThePulse.OpChannel(),
						   PeakTime,      Width,      Area,      Amplitude,      PE,
						   PeakTimeError, WidthError, AreaError, AmplitudeError, PEError);
	      }
	    //  mf::LogInfo("OpHitFinder")<<"Hit found in hist 1 : " <<PeakTime<<" " << Area<<std::endl; 
	  }


	for(size_t hit=0; hit!=FoundHits2.size(); ++hit)
	  {
	    //  mf::LogInfo("OpHitFinder")<<"Found hit in hist 2: " << FoundHits2[hit]<<" " <<Binned2[FoundHits2[hit]]<<std::endl;
	    unsigned int LowBin  = std::max(FoundHits2.at(hit ) * fSamplesPerBin - fSamplesPerBin/2, 0);
	    double Area          =  0;
	    double Amplitude     =  0;
	    double PeakTime      =  0;
	    double Width         =  0;
	    double PE            =  0;

	    double PeakTimeError =  0;
	    double AreaError     =  0;
	    double AmplitudeError=  0;
	    double WidthError    =  0;
	    double PEError       =  0;

	    double HitVar        =  0;
	    
	    double sumQ2=0, sumQ=0;
	    int PeakBin = LowBin;
	    
	    for(size_t binNum = LowBin; binNum!=Unipolar.size(); ++binNum)
	      {
		double ThisVal = Unipolar[binNum];
		sumQ  += ThisVal*binNum;
		sumQ2 += pow(ThisVal,2)*binNum;
		Area  += ThisVal;
		if(ThisVal>Amplitude)
		  {
		    Amplitude = ThisVal;
		    PeakBin   = binNum; 
		  }	     
		if( (int(binNum) - int(PeakBin)) > fSamplesPerBin) break;
		
	      }
	    
	    //   HitMean = sumQ / float(Area);
	    HitVar  = pow(sumQ2 - pow(sumQ,2),0.5)/float(Area);
	   
	    PeakTime      = fTimeBegin*1000. + PeakBin *                1000. / (fSampleFreq);
	    PeakTimeError = HitVar / pow(Area,0.5) * 1000. / (fSampleFreq);
	    Width         = HitVar *                 1000. / (fSampleFreq);

	    PE            = Area / fPEArea;
	    
	    if(Area>0)
	      {
		HitsByTime[PeakBin] = recob::OpHit(ThePulse.OpChannel(),
						   PeakTime,      Width,      Area,      Amplitude,      PE,
						   PeakTimeError, WidthError, AreaError, AmplitudeError, PEError);
	      }
	  }

	for(std::map<int, recob::OpHit>::const_iterator it = HitsByTime.begin();
	    it!=HitsByTime.end();  it++)
	  {
	    StoragePtr->push_back(it->second);
	    
	    //   mf::LogInfo("OpHitFinder") << "Storing hit in opdet " << it->second.OpChannel() <<
	    //  " at time " << it->second.PeakTime() << " with area " << it->second.Area();
	    
	  }

      }


    evt.put(std::move(StoragePtr));

    return;
  }



} // namespace opdet


