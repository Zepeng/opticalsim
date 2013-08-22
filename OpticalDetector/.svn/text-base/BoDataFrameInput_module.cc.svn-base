// \file BoDataFrameInput.h 
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
#ifndef BoDataFrameInput_h
#define BoDataFrameInput_h 1



namespace opdet {

  class BoDataFrameInput : public art::EDProducer{
    public:
      
      BoDataFrameInput(const fhicl::ParameterSet&);
      virtual ~BoDataFrameInput();
      
      void produce(art::Event&);
      
      void beginJob();
      
     
    private:
      
      // The parameters we'll read from the .fcl file.
      std::string   fInputFile;
      std::ifstream fTextFile; 

  };
}

#endif



////////////////////////////////////////////////////////////////////////
/// \file  BoDataFrameInput_module.cc
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  bjpjones
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet{

  DEFINE_ART_MODULE(BoDataFrameInput);

}//end namespace opdet
////////////////////////////////////////////////////////////////////////

// \file BoDataFrameInput.cxx  
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
  

  BoDataFrameInput::BoDataFrameInput(fhicl::ParameterSet const& pset)
  {
    // Infrastructure piece
    produces<std::vector< raw::OpDetPulse> >();

    // Input filename read from fcl
    fInputFile = pset.get<std::string>("InputFile");    
    std::string FullFilePath(""); 
    cet::search_path sp("FW_SEARCH_PATH");
    if( !sp.find_file(fInputFile, FullFilePath) )
      throw cet::exception("BoDataFrameInput") << "Unable to find optical data file in " << sp.to_string();


    // We want to throw away the first command, as it will be "evt"
    std::string sub;
    getline(fTextFile, sub, ',');
    if(sub.substr(0,3)!=std::string("evt"))
      {
	mf::LogInfo("BoDataFrameInput")<<"Warning: first command in text file is not an evt block as expected. Trying to persevere anyway : " << sub; 
      }




    // Open file for reading
    fTextFile.open(FullFilePath);

    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());
    createEngine(seed);

  }
  
  //-------------------------------------------------


  void BoDataFrameInput::beginJob()
  {
  }


  //-------------------------------------------------
  
  BoDataFrameInput::~BoDataFrameInput() 
  {
  }
  

  //-------------------------------------------------

  void BoDataFrameInput::produce(art::Event& evt)
  {
    // Infrastructure piece
    std::unique_ptr<std::vector< raw::OpDetPulse > >  StoragePtr (new std::vector<raw::OpDetPulse>);
 
    // Create vector of pointers to our waveform pulses for each OpDet
    std::vector<raw::OpDetPulse*> ThePulses;

    if(!fTextFile)
      {
	mf::LogError("BoDataFrameInput")<<"Error in reading input file : " << fInputFile;
      }

    // Build an array of comma separated commands

    mf::LogInfo("BoDataFrameInput")<<"Reading optical pulses for event " <<evt.id().event() ;
    std::vector<std::string>  Commands;
    std::vector<unsigned int> Values;

    bool ContinueRead = true;
    std::string sub("");
    std::string line("");

    unsigned int CurrentChannel=0;


    unsigned int ChannelsThisEvent=0;
    unsigned int FEMModule=0;
    unsigned int TriggerID=0;
    unsigned int EventFrameNumber=0;
    unsigned int TriggerTypeThisChannel=0;
    unsigned int PMTFrameNumber=0;

    unsigned int FirstSample=0;

    std::map<int, std::vector<raw::OpDetPulse*> > PulsesThisEvent;
    
    while(fTextFile.good() && ContinueRead)
      {    
	getline(fTextFile, line);
	std::stringstream linestream(line);
	while(linestream.good() && ContinueRead)
	  {
	    getline(linestream, sub, ',');
	    std::stringstream ss(sub);
	    std::string Command("");
	    int Value=0;
	    ss>>Command;
	    ss>>Value;
	    mf::LogInfo("BoDataFrameInput")<<"Parsed a single line as  C: " << Command<< "  V: "<<Value;
	    
	    // Structural commands
	    if(Command=="evt")
	      {
		// This command specifies the beginning of the next readout event.
		//  if we find it, finish processing for now.
		mf::LogInfo("BoDataFrameInput")<<"Found start of next event, ID = " <<  Value<<", finish frame readout";
		ContinueRead=false;
	      }
	    else if(Command=="nch")
	      {
		// Number of channels in this event readout
		if(ChannelsThisEvent==0)
		  ChannelsThisEvent=Value;
		else
		  mf::LogInfo("BoDataFrameInput")<<"Confused by data input: nch specified twice for the same event. Persevering anyway...";	
	      }
	    else if(Command=="mod")
	      {
		// ID of the FEM module
		FEMModule=Value;
		Value=FEMModule; // This line just suppresses a compiler warning
	      }
	    else if(Command=="tid")
	      {
		// Trigger type (cosmic, michel, etc)
		TriggerID=Value;
		Value = TriggerID; // This line just suppresses a compiler warning
	      }
	    else if(Command=="efr")
	      {
		// The event frame number - see longer note under pfr
		EventFrameNumber=Value;
	      }


	    // Per readout window commands
	    
	    else if(Command=="chn")
	      {
		// chn tells us the channel number on the shaper board
		mf::LogInfo("BoDataFrameInput")<<"Beginning channel " <<Value;
	
		CurrentChannel = Value;
		PulsesThisEvent[CurrentChannel].push_back( new raw::OpDetPulse(CurrentChannel) );	
		
		FirstSample    = 0;
		PMTFrameNumber = 0;
	      }
	    else if(Command=="cid")
	      {
		// This gives the trigger type. Encoding not totally clear,
		// but we know 16 = external trigger
		TriggerTypeThisChannel=Value;
		Value =TriggerTypeThisChannel; // This line just suppresses a compiler warning
	      }
	    else if(Command=="smp")
	      {
		// smp tells us which sample to start at within the window
		FirstSample=Value;
		PulsesThisEvent[CurrentChannel].at(PulsesThisEvent[CurrentChannel].size()-1)->SetFirstSample(FirstSample);
	      }
	    else if(Command=="pfr")
	      {
		// This part is a bit complicated.
		// The PMT frame number tells us which readout frame 
		//  this signal occupies.
		// The PMT frame always occurs after the event frame.
		// The pfr value is the last 3 bits of the PMT frame number
		// So the goal is to find the first number after efr with
		//  these last 3 bits.

		// In practice : 
		//  if the last 3 bits of efr are > pfr
		//    the PMT frame is (all but 3 lsb of)efr + 8 + pfr  
		//  if the last 3 bits of efr are < pfr
		//    the PMT frame is (all but 3 lsb of)efr + pfr

		// This mask gets the 3LSB 
		int mask = 0x007;

		int EvLowBits =  EventFrameNumber & mask;
		int EvHighBits = EventFrameNumber & ~mask;

		if(Value < EvLowBits)
		  {
		    PMTFrameNumber = EvHighBits + 8 + Value;
		  }
		else if(Value > EvLowBits)
		  {
		    PMTFrameNumber = EvHighBits + Value;
		  }
		else if(Value == EvLowBits)
		  {
		    PMTFrameNumber = EvHighBits + EvLowBits;
		  }
		else
		  {
		    mf::LogError("BoDataFrameInput") << "pfr routine confused - this should be impossible";
		  }
		
		PulsesThisEvent[CurrentChannel].at(PulsesThisEvent[CurrentChannel].size()-1)->SetPMTFrame(PMTFrameNumber);
		
	      }
	    else if(Command=="adc")
	      {
		// Add one sample to the pulse object
		PulsesThisEvent[CurrentChannel].at(PulsesThisEvent[CurrentChannel].size()-1)->Waveform().push_back(Value);
	      }
	  }  
      }    
    

    // Now pack away the pulses we found into the event
    for(std::map<int, std::vector<raw::OpDetPulse*> >::const_iterator it = PulsesThisEvent.begin();
	it!=PulsesThisEvent.end(); ++it)
      {
	for(size_t pulse=0; pulse!=it->second.size(); ++pulse)
	  {
	    StoragePtr->push_back( *(it->second.at(pulse)) );
	  }
      }
    
    evt.put(std::move(StoragePtr));
  }  
}

