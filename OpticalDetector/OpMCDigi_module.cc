// \file OpMCDigi.h 
// \author Ben Jones and Christie Chiu, MIT, Sept 2012
//   bjpjones@mit.edu, cschiu@mit.edu
//
// This module starts from MC truth sim::OnePhoton objects
// and produces a digitized waveform.
//
// It is assumed that the electronics response is linear,
// and the 1PE waveform can be described by a discreate 
// response shape.  The many PE response is then the linear
// superposition of the relevant function at the appropriate
// arrival times.
//


#include "art/Framework/Core/EDProducer.h"
#include "RawData/OpDetPulse.h"
#include "Simulation/sim.h"
#include "Geometry/Geometry.h"
#include "OpticalDetector/OpDigiProperties.h"


#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoisson.h"


// ROOT includes.
#include <Rtypes.h>
#ifndef OpMCDigi_h
#define OpMCDigi_h 1



namespace opdet {

  class OpMCDigi : public art::EDProducer{
    public:
      
      OpMCDigi(const fhicl::ParameterSet&);
      virtual ~OpMCDigi();
      
      void produce(art::Event&);
      
      void beginJob();
      
     
    private:
      
      // The parameters we'll read from the .fcl file.
      std::string fInputModule;              // Input tag for OpDet collection
      float fSampleFreq;                     // in MHz
      float fTimeBegin;                      // in us
      float fTimeEnd;                        // in us
      float fQE;                             // quantum efficiency of opdet
      float fSaturationScale;                // adc count w/ saturation occurs
    
      float fDarkRate;                      // Noise rate in Hz
    
      std::vector<double> fSinglePEWaveform;
    
      CLHEP::RandFlat    * fFlatRandom;
      CLHEP::RandPoisson * fPoissonRandom;
    

   
    void AddTimedWaveform (int time, std::vector<double>& OldPulse, std::vector<double>& NewPulse);

  };
}

#endif



////////////////////////////////////////////////////////////////////////
/// \file  OpMCDigi_module.cc
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  bjpjones
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet{

  DEFINE_ART_MODULE(OpMCDigi);

}//end namespace opdet
////////////////////////////////////////////////////////////////////////

// \file OpMCDigi.cxx  
// \author Ben Jones and Christie Chiu, MIT 2010
//
//

// FMWK includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Simulation/SimListUtils.h"
#include "Simulation/SimPhotons.h"
#include "RawData/OpDetPulse.h"

// ROOT includes
#include <TH1D.h>
#include <TF1.h>
#include <TTree.h>
#include <TRandom.h>

// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>


// Debug flag; only used during code development.
const bool debug = true;

namespace opdet {
  

  OpMCDigi::OpMCDigi(fhicl::ParameterSet const& pset)
  {
    // Infrastructure piece
    produces<std::vector< raw::OpDetPulse> >();


    // Input Module and histogram parameters come from .fcl
    fInputModule = pset.get<std::string>("InputModule");
 
    art::ServiceHandle<OpDigiProperties> odp;
    fTimeBegin  = odp->TimeBegin();
    fTimeEnd    = odp->TimeEnd();
    fSampleFreq = odp->SampleFreq();

    fQE              = pset.get<double>("QE");
    fDarkRate        = pset.get<double>("DarkRate");
    fSaturationScale = pset.get<double>("SaturationScale");

    // Initialize toy waveform vector fSinglePEWaveform
    fSinglePEWaveform = odp->SinglePEWaveform();


    
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());
    createEngine(seed);

    // Sample a random fraction of detected photons 
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    fFlatRandom       = new CLHEP::RandFlat(engine);
    fPoissonRandom    = new CLHEP::RandPoisson(rng->getEngine());




  }
  
  //-------------------------------------------------


  void OpMCDigi::beginJob()
  {
  }


  //-------------------------------------------------
  
  OpMCDigi::~OpMCDigi() 
  {
  }
  

  //-------------------------------------------------


  void OpMCDigi::AddTimedWaveform (int binTime, std::vector<double>& OldPulse, std::vector<double>& NewPulse)
  {
 
    if( (binTime + NewPulse.size() ) > OldPulse.size())
      OldPulse.resize(binTime + NewPulse.size());
    
    // Add shifted NewWaveform to Waveform at pointer
    for(size_t i = 0; i!=NewPulse.size(); ++i)
      {
        OldPulse.at(binTime+i) += NewPulse.at(i);
      }
  }	  




  //-------------------------------------------------

  void OpMCDigi::produce(art::Event& evt)
  {
    // Infrastructure piece
    std::unique_ptr<std::vector< raw::OpDetPulse > >  StoragePtr (new std::vector<raw::OpDetPulse>);

    
    // Read in the Sim Photons
    sim::SimPhotonsCollection ThePhotCollection = sim::SimListUtils::GetSimPhotonsCollection(evt,fInputModule);

    
    double TimeBegin_ns  = fTimeBegin  *  1000;
    double TimeEnd_ns    = fTimeEnd    *  1000;
    double SampleFreq_ns = fSampleFreq /  1000;

    int nSamples = ( TimeEnd_ns-TimeBegin_ns)*SampleFreq_ns;
    
    art::ServiceHandle<geo::Geometry> geom;
    int NOpChannels = geom->NOpChannels();


    // This vector will store all the waveforms we will make
    std::vector<std::vector<double> > PulsesFromDetPhotons(NOpChannels,std::vector<double>(nSamples,0.0)); 
    
    // For every OpDet:
    for(sim::SimPhotonsCollection::const_iterator itOpDet=ThePhotCollection.begin(); itOpDet!=ThePhotCollection.end(); itOpDet++)
      {
	const sim::SimPhotons* ThePhot=itOpDet->second;
	
	int Ch = ThePhot->OpChannel();
	
	// For every photon in the hit:
	for(sim::SimPhotons::const_iterator itPhot = ThePhot->begin(); itPhot!=ThePhot->end(); itPhot++)
	  {
	    // Sample a random subset according to QE
	    if(fFlatRandom->fire(1.0)<=fQE)
	      {
		
		// Convert photon arrival time to the appropriate bin, dictated by fSampleFreq. Photon arrival time is in ns, beginning time in us, and sample frequency in MHz. Notice that we have to accommodate for the beginning time
		if((itPhot->Time > TimeBegin_ns) && (itPhot->Time < TimeEnd_ns))
		  {
		    int binTime = int((itPhot->Time - TimeBegin_ns) * SampleFreq_ns);		
	
		    // Call function to add waveforms to our pulse
		    AddTimedWaveform( binTime, PulsesFromDetPhotons[Ch], fSinglePEWaveform );
		    
		  }
	      } // random QE cut
	  } // for each Photon in SimPhotons
	
      }

    art::ServiceHandle<art::RandomNumberGenerator> rng;

    
    // Create vector of output objects, add dark noise and apply
    //  saturation

    std::vector<raw::OpDetPulse*> ThePulses(NOpChannels);
    for(int iCh=0; iCh!=NOpChannels; ++iCh)
      {
	PulsesFromDetPhotons[iCh].resize((TimeEnd_ns - TimeBegin_ns) * SampleFreq_ns);

	// Add dark noise
	double MeanDarkPulses = fDarkRate * (fTimeEnd-fTimeBegin) / 1000000;

	unsigned int NumberOfPulses = fPoissonRandom->fire(MeanDarkPulses);
	
	for(size_t i=0; i!=NumberOfPulses; ++i)
	  {
	    double PulseTime = (fTimeEnd-fTimeBegin)*fFlatRandom->fire(1.0);
	    int binTime = int(PulseTime * fSampleFreq);
	    
	    AddTimedWaveform( binTime, PulsesFromDetPhotons[iCh], fSinglePEWaveform );
	  }

	// Apply saturation for large signals
	for(size_t i=0; i!=PulsesFromDetPhotons[iCh].size(); ++i)
	  {
	    if(PulsesFromDetPhotons[iCh].at(i)>fSaturationScale) PulsesFromDetPhotons[iCh].at(i) = fSaturationScale;
	  }

	// Produce ADC pulse of integers rather than doubles
	
	std::vector<short> shortvec;
	
	for(size_t i=0; i!=PulsesFromDetPhotons[iCh].size(); ++i)
	  {
	    // Throw randoms to fairly sample +ve and -ve side of doubles
	    int ThisSample = PulsesFromDetPhotons[iCh].at(i);
	    if(ThisSample>0)
	      {
		if(fFlatRandom->fire(1.0) > (ThisSample - int(ThisSample)))
		  shortvec.push_back(int(ThisSample));
		else
		  shortvec.push_back(int(ThisSample)+1);
	      }
	    else
	      {
		if(fFlatRandom->fire(1.0) >  (int(ThisSample)-ThisSample))
		  shortvec.push_back(int(ThisSample));
		else
		  shortvec.push_back(int(ThisSample)-1);
	      }
	  }

	
	StoragePtr->push_back(	raw::OpDetPulse( iCh, shortvec ,0, fTimeBegin) );
	
      } // for each OpDet in SimPhotonsCollection
    

    evt.put(std::move(StoragePtr));
  }  
}

