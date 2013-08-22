// Christie Chiu and Ben Jones, MIT, 2012
//
// This module finds large pulses in an optical detector
// output waveform and produces OpHit objects representing
// discrete N.PE pulses. 
//  
// These OpHits are stored into the event for further
// analysis.
//

// OpLowIntensityHitFinder.h

#ifndef OpLowIntensityHitFinder_H
#define OpLowIntensityHitFinder_H 1

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
 
  class OpLowIntensityHitFinder : public art::EDProducer{
  public:
 
    // Standard constructor and destructor for an ART module.
    OpLowIntensityHitFinder(const fhicl::ParameterSet&);
    virtual ~OpLowIntensityHitFinder();

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
    float fTimeBegin;                      // in us
    float fTimeEnd;                        // in us
    short fPEheight;                       // in ADC counts

  };

} 

#endif // OpLowIntensityHitFinder_H




// OpLowIntensityHitFinder_module.cc

// This is a short program required by the ART framework.  It enables
// a program (OpLowIntensityHitFinder, in this case) to be called as a module
// from a .fcl file. It is unlikely that you'll have to make any
// changes to this file.

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet {
  DEFINE_ART_MODULE(OpLowIntensityHitFinder);
}


// OpLowIntensityHitFinder.cxx

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
  OpLowIntensityHitFinder::OpLowIntensityHitFinder(fhicl::ParameterSet const& pset)
  {
    // Infrastructure piece
    produces<std::vector< recob::OpHit> >();

    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get<std::string>("InputModule");

    art::ServiceHandle<OpDigiProperties> odp;
    fTimeBegin  = odp->TimeBegin();
    fTimeEnd    = odp->TimeEnd();
    fSampleFreq = odp->SampleFreq();

    // Get 1PE pulse parameters from .fcl
    fPEheight = pset.get<short>("PEheight");
  }

  //-----------------------------------------------------------------------
  // Destructor
  OpLowIntensityHitFinder::~OpLowIntensityHitFinder() 
  {}
   
  //-----------------------------------------------------------------------
  void OpLowIntensityHitFinder::beginJob()
  {
  }
   

  //-----------------------------------------------------------------------
  void OpLowIntensityHitFinder::produce(art::Event& evt) 
  {
    // Infrastructure piece
    std::unique_ptr<std::vector< recob::OpHit > > StoragePtr (new std::vector<recob::OpHit>);
    
    // Create a handle for our vector of pulses
    art::Handle< std::vector< raw::OpDetPulse > > WaveformHandle;

    // Create string for histogram name
    char HistName[50];
    char FitName[50];
    char FitEq[200];

    // Read in WaveformHandle
    evt.getByLabel(fInputModule, WaveformHandle);

    // Access ART's TFileService, which will handle creating and writing
    // histograms for us.
    art::ServiceHandle<art::TFileService> tfs;


    // For every OpDet waveform in the vector given by WaveformHandle
    for(unsigned int i = 0; i < WaveformHandle->size(); ++i)
      {
	// Get OpDetPulse
	art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle, i);
	raw::OpDetPulse ThePulse = *ThePulsePtr;

	// Make an instance of histogram and its pointer, changing all units to ns
	// Notice that histogram axis is in ns but binned by 1/64MHz;
	//   appropriate conversions are made from beginning and end time 
	//   in us, and frequency in MHz.
	sprintf(HistName, "Event %d OpDet %i", evt.id().event(), ThePulse.OpChannel());
	TH1D* myHistPointer= tfs->make<TH1D>(HistName, ";t (ns);", 
					     int((fTimeEnd - fTimeBegin) * fSampleFreq), 
					     fTimeBegin * 1000., 
					     fTimeEnd * 1000.);


	// Make vectors to keep track of gaussian fit ranges
	std::vector<unsigned int> fitBeginBin;
	std::vector<unsigned int> fitEndBin;
	std::vector<unsigned int> pulseMaxBin;


	// Initiate state for finding gaussian fit ranges
	// state 0: looking for first minimum (fit beginning)
	// state 1: looking for peak maximum
	// state 2: looking for second minimum (fit end)
	short state = 0;

	for(unsigned int binNum = 0; binNum < ThePulse.Waveform().size(); ++binNum)
	  {
	    // Fill histogram with waveform
	    // Raw signal is positive, so to prevent wraparound we set very negative 
	    //   values to unsigned short (this is probably not the best fix)

	    if (binNum < 3)
	      {
		myHistPointer->SetBinContent( binNum, 0);
		continue;
	      }

	    else
	      myHistPointer->SetBinContent( binNum,
					    (double) ThePulse.Waveform()[binNum] + 
					    (double) myHistPointer->GetBinContent(binNum-3));

	    // Looking for pulse beginning
	    if (state == 0)
	      {
		// If signal begins increasing
		if (myHistPointer->GetBinContent(binNum-3) < myHistPointer->GetBinContent(binNum))
		  { 
		    fitBeginBin.push_back(binNum-3);
		    state = 1; 
		  }


		// Else if signal begins decreasing
		else if (myHistPointer->GetBinContent(binNum-3) > myHistPointer->GetBinContent(binNum) && 
			 myHistPointer->GetBinContent(binNum-3) > fPEheight)
		  {
		    fitBeginBin.push_back(binNum-3);
		    pulseMaxBin.push_back(binNum-2);
		    state = 2; 
		  }

	      }

	    
	    // Looking for pulse peak
	    else if (state == 1)
	      {
		// If signal begins decreasing and peak is high enough
		if (myHistPointer->GetBinContent(binNum-3) > myHistPointer->GetBinContent(binNum) && 
		    myHistPointer->GetBinContent(binNum-3) > fPEheight)
		  { 
		    pulseMaxBin.push_back(binNum-3);
		    state = 2;
		  }

		// If we have an inflection point far from peak, take note
		if (binNum > 9)
		  {
		    if (myHistPointer->GetBinContent(binNum-6) - 2*myHistPointer->GetBinContent(binNum-3) + myHistPointer->GetBinContent(binNum) > 
			0 &&
			myHistPointer->GetBinContent(binNum-9) - 2*myHistPointer->GetBinContent(binNum-6) + myHistPointer->GetBinContent(binNum-3) < 
			0 && 
			myHistPointer->GetBinContent(binNum-3) > fPEheight &&
			myHistPointer->GetBinContent(binNum-6) - myHistPointer->GetBinContent(binNum-3) >= 
			(myHistPointer->GetBinContent(binNum-9) - myHistPointer->GetBinContent(binNum-6))/1.1)
		      { 
			if (pulseMaxBin.empty())
			  {
			    fitEndBin.push_back(binNum-6);
			    fitBeginBin.push_back(binNum-6);
			    pulseMaxBin.push_back(binNum-3);
			    continue;
			  }
			else if (binNum - pulseMaxBin.back() > 15)
			  {
			    fitEndBin.push_back(binNum-6);
			    fitBeginBin.push_back(binNum-6);
			    pulseMaxBin.push_back(binNum-3);
			  }
		      }		
		  }


	      }
	    

	    // Looking for pulse end
	    else if (state == 2)
	      {
		// If signal begins increasing
		if (myHistPointer->GetBinContent(binNum-3) < myHistPointer->GetBinContent(binNum))
		  { 
		    fitEndBin.push_back(binNum-3);
		    state = 0; 
		  }

		// If we have an inflection point far from peak, take note
		if (myHistPointer->GetBinContent(binNum-6) - 2*myHistPointer->GetBinContent(binNum-3) + myHistPointer->GetBinContent(binNum) < 
		    0 &&
		    myHistPointer->GetBinContent(binNum-9) - 2*myHistPointer->GetBinContent(binNum-6) + myHistPointer->GetBinContent(binNum-3) > 
		    0 && 
		    myHistPointer->GetBinContent(binNum-3) > fPEheight && 
		    binNum - pulseMaxBin.back() > 5 && 
		    myHistPointer->GetBinContent(binNum-6) - myHistPointer->GetBinContent(binNum-3) <= 
		    (myHistPointer->GetBinContent(binNum-9) - myHistPointer->GetBinContent(binNum-6))/1.1)
		  { 
		    fitEndBin.push_back(binNum-6);
		    fitBeginBin.push_back(binNum-6);
		    pulseMaxBin.push_back(binNum-3);
		  }		
	      }

	  } // for each bin

	// ignore mostly cut off pulses at the end; include mostly complete pulses
	if (state == 1) 
	  fitBeginBin.pop_back();
	if (state == 2)
	  fitEndBin.push_back(ThePulse.Waveform().size());


	// Now that we have the regions to fit gaussians to, fit gaussians!

	// Array to store fit parameters
	Double_t par[3*fitBeginBin.size()];

	for (unsigned int peakNum = 0; peakNum < fitBeginBin.size(); ++peakNum)
	  {
	    // Make fit title
	    sprintf(FitName, "Event %d OpDet %i Peak %d", evt.id().event(), ThePulse.OpChannel(), peakNum);

	    if (peakNum == 0)
	      sprintf(FitEq, "gaus(0)");
	    else
	      sprintf(FitEq, "%s+gaus(%d)", FitEq, peakNum*3);

	    // Fit to individual gaussians to get seeds
	    TF1* myFitPointer = new TF1(FitName, "gaus", 
					((double) fitBeginBin[peakNum] + (fTimeBegin * fSampleFreq)) / fSampleFreq * 1000., 
					((double) fitEndBin[peakNum] + (fTimeBegin * fSampleFreq)) / fSampleFreq * 1000.);
	    myFitPointer->SetParLimits(0, 0, 1000000);
	    myFitPointer->SetParLimits(2, 1000 / fSampleFreq, 1000000);
	    myHistPointer->TH1::Fit(FitName, "QR+");
	    myFitPointer->GetParameters(&par[3*peakNum]);

	  } // for each peak


	// Now do a total, sum of peaks fit! Only if we have peaks
	if (fitBeginBin.size() == 0)
	  continue;

	sprintf(FitName, "Event %d OpDet %i fit", evt.id().event(), ThePulse.OpChannel());
	TF1* myFitPointer = new TF1(FitName, FitEq, 
				    ((double) (fTimeBegin * fSampleFreq)) / fSampleFreq * 1000., 
				    //((double) fitBeginBin[0] + (fTimeBegin * fSampleFreq)) / fSampleFreq * 1000., 
				    ((double) fitEndBin[ThePulse.Waveform().size()] + (fTimeBegin * fSampleFreq)) / fSampleFreq * 1000.);

	// Set parameter seeds as fit parameters
	myFitPointer->SetParameters(par);


	// Set parameter limits based on what we know 
	for (unsigned int i = 0; i < fitBeginBin.size(); ++i)
	  {
	    myFitPointer->SetParLimits(3*i, 0., 1000000);
	    myFitPointer->SetParLimits(3*i+1, 
	    			       (fitBeginBin[i] + (fTimeBegin * fSampleFreq)) / fSampleFreq * 1000.,
	    			       (fitEndBin[i] + (fTimeBegin * fSampleFreq)) / fSampleFreq * 1000.);
				       //(pulseMaxBin[i] - 3 + (fTimeBegin * fSampleFreq)) / fSampleFreq * 1000.,
				       //(pulseMaxBin[i] + 3 + (fTimeBegin * fSampleFreq)) / fSampleFreq * 1000.);
	    myFitPointer->SetParLimits(3*i+2, 1000 / fSampleFreq, 100000/ fSampleFreq);
	    //myFitPointer->SetParLimits(3*i+2, 20, 1000000);
	  }

	myFitPointer->SetLineColor(kBlue);
	myHistPointer->TH1::Fit(FitName, "QR+");

	// Add up total number of photons for OpDet signal
	int photonSum = 0;
	mf::LogInfo ("Peak Info") << "Event: " << evt.id().event() << 
	  "   OpDet: " << ThePulse.OpChannel();


	for (unsigned int i = 0; i < fitBeginBin.size(); ++i)
	  {
	    // Save peak information to event
	    StoragePtr->push_back( recob::OpHit( ThePulse.OpChannel(),
						 (double) myFitPointer->GetParameter(3*i+1),
						 (double) myFitPointer->GetParameter(3*i+2),
						 (double) myFitPointer->GetParameter(3*i) * myFitPointer->GetParameter(3*i+2),
						 (double) myFitPointer->GetParameter(3*i),
						 (double) myFitPointer->GetParameter(3*i) * myFitPointer->GetParameter(3*i+2) / 11572.5764,
						 (double) myFitPointer->GetParError(3*i+1),
						 (double) myFitPointer->GetParError(3*i+2),
						 (double) pow( pow( myFitPointer->GetParError(3*i) * myFitPointer->GetParameter(3*i+2), 2) +
							       pow( myFitPointer->GetParError(3*i+2) * myFitPointer->GetParameter(3*i), 2), 0.5),
						 (double) myFitPointer->GetParError(3*i),
						 (double) pow( pow( myFitPointer->GetParError(3*i) * myFitPointer->GetParameter(3*i+2), 2) +
							       pow( myFitPointer->GetParError(3*i+2) * myFitPointer->GetParameter(3*i), 2), 0.5) / 11572.5764));

	    // Print peak information to screen
	    photonSum += (int) myFitPointer->GetParameter(3*i) * myFitPointer->GetParameter(3*i+2) / 11572.5764 + 0.5;
	    mf::LogInfo ("Peak Info") << 
	      "   Time: " << myFitPointer->GetParameter(3*i+1) <<
	      "   No. PE: " << (int) (myFitPointer->GetParameter(3*i) * myFitPointer->GetParameter(3*i+2) / 11572.5764 + 0.5);
	  }

	mf::LogInfo ("Peak Info") << "     TOTAL NO. PE: " << photonSum;

	
      } // for each OpDet in SimPhotonsCollection

    evt.put(std::move(StoragePtr));

    return;
  }


} // namespace opdet


