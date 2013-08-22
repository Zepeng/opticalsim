// \file StitchPMTFrames.h 
// \author Ben Jones, MIT, Jan 2013
//   bjpjones@mit.edu
//


#include "art/Framework/Core/EDProducer.h"
#include "RawData/OpDetPulse.h"
#include "Simulation/sim.h"
#include "OpticalDetector/OpDigiProperties.h"


#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoisson.h"


// ROOT includes.
#include <Rtypes.h>
#ifndef StitchPMTFrames_h
#define StitchPMTFrames_h 1



namespace opdet {

  class StitchPMTFrames : public art::EDProducer{
    public:
      
      StitchPMTFrames(const fhicl::ParameterSet&);
      virtual ~StitchPMTFrames();
      
      void produce(art::Event&);
      
      void beginJob();
      
     
    private:
      
      // The parameters we'll read from the .fcl file.
      std::string   fInputModule;
      int           fFrameSize;
      int           fADCOffset;

  };
}

#endif



////////////////////////////////////////////////////////////////////////
/// \file  StitchPMTFrames_module.cc
///
/// \author  bjpjones
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet{

  DEFINE_ART_MODULE(StitchPMTFrames);

}//end namespace opdet
////////////////////////////////////////////////////////////////////////

// \file StitchPMTFrames.cxx  
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
  

  StitchPMTFrames::StitchPMTFrames(fhicl::ParameterSet const& pset)
  {
    // Infrastructure piece
    produces<std::vector< raw::OpDetPulse> >();

    // Get frame size specification
    fFrameSize = pset.get<int>("FrameSize");    

    // Get offset value
    fADCOffset = pset.get<int>("ADCOffset");    

    // Input filename read from fcl
    fInputModule = pset.get<std::string>("InputModule");    
  
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());
    createEngine(seed);

  }
  
  //-------------------------------------------------


  void StitchPMTFrames::beginJob()
  {
  }


  //-------------------------------------------------
  
  StitchPMTFrames::~StitchPMTFrames() 
  {
  }
  

  //-------------------------------------------------

  void StitchPMTFrames::produce(art::Event& evt)
  {
    // Infrastructure piece
    std::unique_ptr<std::vector< raw::OpDetPulse > >  StoragePtr (new std::vector<raw::OpDetPulse>);
 
    // vector to store pulses to go into event
    std::vector<raw::OpDetPulse*> PulsesToStore;

    // Create a handle for our vector of input pulses
    art::Handle< std::vector< raw::OpDetPulse > > WaveformHandle;

    // Read in WaveformHandle
    evt.getByLabel(fInputModule, WaveformHandle);


    std::map<unsigned short, std::vector<raw::OpDetPulse* > > OrganizedByChannel;
    std::map<unsigned short, unsigned int>  EarliestFramePerChannel;

    for(unsigned int i = 0; i < WaveformHandle->size(); ++i)
      {
	// Sort pulses by channel
	raw::OpDetPulse * ThePulsePtr = new raw::OpDetPulse(WaveformHandle->at(i));
	unsigned short ThisChannel = ThePulsePtr->OpChannel(); 
	OrganizedByChannel[ThisChannel].push_back(ThePulsePtr);	

	// And find the one that comes first in time in each channel
	
	if( (OrganizedByChannel.size()==1) ||
	    (EarliestFramePerChannel[ThisChannel] > ThePulsePtr->PMTFrame() ) )
	  EarliestFramePerChannel[ThisChannel] = ThePulsePtr->PMTFrame();
      }
    
    for(std::map<unsigned short, std::vector<raw::OpDetPulse* > >::iterator it = OrganizedByChannel.begin(); it!=OrganizedByChannel.end(); ++it)
      {
	// For each channel, make a new pulse
	raw::OpDetPulse * ThisPulse = new raw::OpDetPulse(it->first);
	PulsesToStore.push_back(ThisPulse );
	
	// Its frame number will be the earliest of those found for this channel
	ThisPulse->SetPMTFrame(EarliestFramePerChannel[it->first]);
	ThisPulse->SetFirstSample(0);
	
	for(size_t ipulsethischan = 0; ipulsethischan!=it->second.size(); ++ipulsethischan)
	  {
	    unsigned int Samples     = it->second.at(ipulsethischan)->Samples();
	    unsigned int FrameOffset = fFrameSize * (it->second.at(ipulsethischan)->PMTFrame() - ThisPulse->PMTFrame());
	    unsigned int FirstSample = it->second.at(ipulsethischan)->FirstSample();
	    unsigned int ThisPulseStart = FrameOffset + FirstSample;
	    

	    if( (ThisPulseStart+Samples) > ThisPulse->Waveform().size() )
	      {
		ThisPulse->Waveform().resize(ThisPulseStart + Samples);

	      }
	    
	    for(size_t isamp=0; isamp!=Samples; ++isamp)
	      {
		ThisPulse->Waveform().at(ThisPulseStart+isamp)+=
		  it->second.at(ipulsethischan)->Waveform().at(isamp) - fADCOffset;
	      }
   
	  }
      }
    
        

    // Delete the copies we made of the pulses
    for(std::map<unsigned short, std::vector<raw::OpDetPulse* > >::iterator it=OrganizedByChannel.begin(); it!=OrganizedByChannel.end(); ++it)
      {
	for(size_t i=0; i!=it->second.size(); ++i)
	  { 
	    if(it->second.at(i ) )
	      delete it->second.at(i);
	  }
      }
    OrganizedByChannel.clear();



    // Store the stitched pulses into the event
    for(size_t ichan=0; ichan!=PulsesToStore.size(); ++ichan)
      {
      	StoragePtr->push_back(*PulsesToStore.at(ichan));
      }

    evt.put(std::move(StoragePtr));
  }  
}

