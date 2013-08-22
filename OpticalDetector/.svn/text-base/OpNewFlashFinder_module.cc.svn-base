// Ben Jones, MIT, 2012
//
// This module finds large pulses in an optical detector
// output waveform and produces OpHit objects representing
// discrete N.PE pulses. 
//  
// These OpHits are stored into the event for further
// analysis.
//

// OpNewFlashFinder.h

#ifndef OpNewFlashFinder_H
#define OpNewFlashFinder_H 1

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/OpDetGeo.h"
#include "RawData/OpDetPulse.h"

// ART includes.
#include "art/Framework/Core/EDProducer.h"

// C++ includes
#include <cstring>
#include <vector>

namespace opdet {
 
  class OpNewFlashFinder : public art::EDProducer{
  public:
 
    // Standard constructor and destructor for an ART module.
    OpNewFlashFinder(const fhicl::ParameterSet&);
    virtual ~OpNewFlashFinder();

    // This method is called once, at the start of the job. In this
    // example, it will define the histogram we'll write.
    void beginJob();

    // The producer routine, called once per event. 
    void produce (art::Event&); 

  private:

    // The parameters we'll read from the .fcl file.
    std::string fInputModule;              // Input tag for OpDet collection
    float fSampleFreq;                     // in MHz
    float fTimeBegin;
    float fTimeEnd;
    int   fSamplesPerBin;
    int   fMaxHitLength;
    int   fFastCompSamples;
    int   fPreSamples;
    float fZeroSupThresh;                  // As fraction of PE
    float fPEArea;                         // in ADC counts
    float fPEAmp;                         // in ADC counts
  
    float fFlashThreshold; 
    float fHitThreshold;
    
    float fIsolationFrac;
      
    float fBeamStartTime;                  // in us
    float fBeamSpillLength;                // in us



  };

} 

#endif // OpNewFlashFinder_H




// OpNewFlashFinder_module.cc

// This is a short program required by the ART framework.  It enables
// a program (OpNewFlashFinder, in this case) to be called as a module
// from a .fcl file. It is unlikely that you'll have to make any
// changes to this file.

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet {
  DEFINE_ART_MODULE(OpNewFlashFinder);
}


// OpNewFlashFinder.cxx

// LArSoft includes
#include "RawData/OpDetPulse.h"
#include "RecoBase/OpHit.h"
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
  OpNewFlashFinder::OpNewFlashFinder(fhicl::ParameterSet const& pset)
  {
    // Infrastructure piece
    produces<std::vector< recob::OpFlash> >();
    produces<std::vector< recob::OpHit> >();

    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get<std::string>("InputModule");

    art::ServiceHandle<OpDigiProperties> odp;
    fTimeBegin  = odp->TimeBegin();
    fTimeEnd    = odp->TimeEnd();
    fSampleFreq = odp->SampleFreq();
    fPEArea     = odp->GetSPECumulativeArea();   // Uncomment if unipolar to bipolar inversion offline
    // fPEArea     = odp->GetSPEArea();          // Uncomment if unipolar to bipolar inversion online

    fPEAmp     = odp->GetSPECumulativeAmplitude();   // Uncomment if unipolar to bipolar inversion offline
    // fPEAmp     = odp->GetSPEAmplitude();          // Uncomment if unipolar to bipolar inversion online


    // Get binning parameters from .fcl
    fSamplesPerBin  = pset.get<int>("SamplesPerBin");                // For broadly binned sum pulse, how many samples go in each bin
    fMaxHitLength   = pset.get<int>("MaxHitLength");                 // How long to look after bin start for pulse activity 
    fFastCompSamples= pset.get<int>("FastCompSamples");              // Number of samples counting as fast for fast/total ratio
    fPreSamples     = pset.get<int>("PreSamples");                   // Number of samples before the peak to count
    fIsolationFrac  = pset.get<float>("IsolationFrac");              // Fraction of central bin allowed in side bins
    fHitThreshold   = pset.get<float>("HitThreshold")  * fPEArea;    // Minimum activity (in PE) in one channel to contribute
    fFlashThreshold = pset.get<float>("FlashThreshold")* fPEArea;    // Minimum activity (in PE) in whole flash
    fZeroSupThresh  = pset.get<float>("ZeroSupThresh") * fPEAmp;     // Threshold to suppress low activity bins in inverted pulse
  
    fBeamStartTime   = pset.get<float>("BeamStartTime");             // Beam window start time (us)
    fBeamSpillLength = pset.get<float>("BeamSpillLength");           // Beam window length (us)


  }

  //-----------------------------------------------------------------------
  // Destructor
  OpNewFlashFinder::~OpNewFlashFinder() 
  {}
   
  //-----------------------------------------------------------------------
  void OpNewFlashFinder::beginJob()
  {
  }
   

  //-----------------------------------------------------------------------
  void OpNewFlashFinder::produce(art::Event& evt) 
  {

    // Infrastructure pieces
    std::unique_ptr<std::vector< recob::OpHit > >   HitPtr (new std::vector<recob::OpHit >);
    std::unique_ptr<std::vector< recob::OpFlash > > FlashPtr (new std::vector<recob::OpFlash >);
								
    // Create a handle for our vector of pulses
    art::Handle< std::vector< raw::OpDetPulse > > WaveformHandle;
    
    // Read in WaveformHandle
    evt.getByLabel(fInputModule, WaveformHandle);
    
    std::vector<std::vector<int> > UnipolarWaveforms;
    std::vector<int>               OpChannelIDs;
    
    std::vector<int> Binned1;
    std::vector<int> Binned2;
    
    
    
    // This first thing we do is:
    // - Invert all the waveforms to unipolar
    // - Add them all up into a single collective pulse, with broad binning

    //    mf::LogInfo("OpNewFlashFinder")<<"Producing unipolar and binned pulses";

     // For every OpDet waveform in the vector given by WaveformHandle
     for(unsigned int i = 0; i < WaveformHandle->size(); ++i)
       {

 	art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle, i);
 	raw::OpDetPulse ThePulse(*ThePulsePtr);
	
 	unsigned int Samples = ThePulse.Waveform().size();

 	if(Samples > Binned1.size())
 	  {
 	    Binned1.resize(int( Samples / fSamplesPerBin + 1)); 
 	    Binned2.resize(int( Samples / fSamplesPerBin + 1)); 
	  }
	
	OpChannelIDs.push_back(ThePulse.OpChannel());

 	UnipolarWaveforms.push_back(std::vector<int>(Samples));
 	size_t ThisIndex = UnipolarWaveforms.size()-1;
	
 	for(unsigned int binNum = 0; binNum < ThePulse.Waveform().size(); ++binNum)
 	  {

 	    // We make the unipolar histograms per channel as the cumulative 
 	    //  pulse to that point

 	    UnipolarWaveforms[ThisIndex][binNum]  
 	      = (ThePulse.Waveform()[binNum]);

	    // Form unipolar pulse, and add a slow decay to prevent baseline drift.
	    //  Its ugly to do this, but we'll know better what to do once
	    //  we have a real shaper simulation.

 	    if(binNum>1) 
	      {
		UnipolarWaveforms[ThisIndex][binNum] 
		  += UnipolarWaveforms[ThisIndex][binNum-1];
		if( abs(UnipolarWaveforms[ThisIndex][binNum]-UnipolarWaveforms[ThisIndex][binNum-1])<0.01)
		  UnipolarWaveforms[ThisIndex][binNum]*=0.95;
	      }

 	    // To prevent floating point or noise affecting the baseline, we zero
 	    //  suppress samples much below 1pe.
	    	    
 	    if(UnipolarWaveforms[ThisIndex][binNum] < (fZeroSupThresh)) 
 	      UnipolarWaveforms[ThisIndex][binNum]=0;

	   	   	    
 	    Binned1[int((float)binNum / (float)fSamplesPerBin)          ] += UnipolarWaveforms[ThisIndex][binNum];
 	    Binned2[int(((float)binNum / (float)fSamplesPerBin) + 0.5)  ] += UnipolarWaveforms[ThisIndex][binNum];
	    
 	  }
       }


     //
     // Next we'll run through that broadly binned pulse looking for
     //  isolated hits.
     //


     //    mf::LogInfo("OpNewFlashFinder")<<"Searching for isolated flashes";

     std::vector<int> FoundHits1;
     std::vector<int> FoundHits2;

     for(size_t binNum=0; binNum!=Binned1.size(); ++binNum)
       {
	 //first, check Binned1
	 float val_this = (float)Binned1[binNum];
	 float val_prev = 0;
	 if(binNum>0) val_prev = (float)Binned1[binNum-1];
	 float val_next = 0;
	 if((binNum+1)<Binned1.size()) val_next = (float)Binned1[binNum+1];
	 
	 //std::cout << "Binned1 (binNum=" << binNum << "):  " << val_prev << " " << val_this << " " << val_next << std::endl;

	 if( val_this!=0)
	   {	     
	     if( ( val_prev==0 || (val_prev/val_this < fIsolationFrac) ) &&
		 ( val_next==0 || (val_next/val_this < fIsolationFrac) ) &&
		 val_this > (fFlashThreshold) )
	       FoundHits1.push_back(binNum);
	   }
	 
	 //now, check Binned2
	 val_this = (float)Binned2[binNum];
	 val_prev = 0;
	 if(binNum>0) val_prev = (float)Binned2[binNum-1];
	 val_next = 0;
	 if((binNum+1)<Binned2.size()) val_next = (float)Binned2[binNum+1];

	 //std::cout << "Binned2 (binNum=" << binNum << "):  " << val_prev << " " << val_this << " " << val_next << std::endl;

	 if( val_this!=0)
	   {	     
	     if( ( val_prev==0 || (val_prev/val_this < fIsolationFrac) ) &&
		 ( val_next==0 || (val_next/val_this < fIsolationFrac) ) &&
		 val_this > (fFlashThreshold) )
	       FoundHits2.push_back(binNum);
	   }
       }

     //    mf::LogInfo("OpNewFlashFinder")<<"Organizing ophits into maps";

     // we store hit times for each flash in a map, indexed by
     //  flash and opchannel.  We will go through removing 
     //  duplicates later
     //
     // Flashes from the second binning scheme will be given negative ID's.
     //  these ID's are transient and do not get saved.
     //
     std::map<int, std::map<int, int> > HitTimes;
     
     
     // We store the actual hits we find in a map indexed by their channel and time
     std::map<int, std::map<int, recob::OpHit> > AllOpHits;
     
     
     // Loop over all waveforms, and look to see if they have
     //  activity in the vicinity of the flash we found.
     for(size_t wf=0; wf!=UnipolarWaveforms.size(); ++wf)
       {	    
	 // Check whether this wf has activity for each flash 
	 for(size_t flash=0; flash!=FoundHits1.size(); ++flash)
	   {
	     unsigned int LowBin   = FoundHits1.at(flash ) * fSamplesPerBin;
	     double Area           =  0;
	     double Amplitude      =  0;
	     double PeakTime       =  0;
	     double Width          =  0;
	     double PE             =  0;
	     
	     double HitVar         =  0;
	     double FastComponent  =  0;
	     double FastToTotal    =  0;

	     double sumQ2=0, sumQ=0;
	     int PeakBin = LowBin;
	     
	     for(size_t binNum = LowBin; binNum!=(LowBin+fSamplesPerBin); ++binNum)
	       {
		 if(UnipolarWaveforms[wf][binNum]>Amplitude)
		   {
		     Amplitude = UnipolarWaveforms[wf][binNum];
		     PeakBin   = binNum;
		   }
	       }
	     if(Amplitude>0)
	       {
		 for(int binNum = PeakBin-fPreSamples; binNum!=PeakBin+fMaxHitLength; ++binNum)
		   {
		     if( (binNum>0) && (binNum < int(UnipolarWaveforms[wf].size()) ) )
		       {
			 double ThisVal = UnipolarWaveforms[wf][binNum];
			 sumQ  += ThisVal*binNum;
			 sumQ2 += pow(ThisVal,2)*binNum;
			 Area  += ThisVal;
			 
			 if( (binNum - PeakBin) < fFastCompSamples ) FastComponent += UnipolarWaveforms[wf][binNum];
			 
		       } 
		 
		   }

		 
		 HitVar      = pow(sumQ2 - pow(sumQ,2),0.5)/float(Area);  
		 PeakTime    = fTimeBegin*1000. + PeakBin *    1000. / fSampleFreq;
		 Width       = HitVar *                        1000. / fSampleFreq;
		 PE          = Area / fPEArea;
		 FastToTotal = FastComponent / Area;
		
		 if(Area>(fHitThreshold))
		   {
		     int ThisChannel = OpChannelIDs[wf];
		     HitTimes[flash][ThisChannel]    = PeakBin;
		     AllOpHits[ThisChannel][PeakBin] =  recob::OpHit( ThisChannel,
								      PeakTime,  
								      Width,
								      Area,      
								      Amplitude,      
								      PE,
								      FastToTotal);
		 
		 
		   }
	       }
	   }
       
	 
	 // Check for activity in the second binning scheme
	 for(size_t flash=0; flash!=FoundHits2.size(); ++flash)
	   {
	     unsigned int LowBin  = std::max(FoundHits2.at(flash ) * fSamplesPerBin - fSamplesPerBin/2, 0);
	     double Area           =  0;
             double Amplitude      =  0;
             double PeakTime       =  0;
             double Width          =  0;
             double PE             =  0;

             double HitVar         =  0;
             double FastComponent  =  0;
             double FastToTotal    =  0;

             double sumQ2=0, sumQ=0;
             int PeakBin = LowBin;

             for(size_t binNum = LowBin; binNum!=(LowBin+fSamplesPerBin); ++binNum)
	       {
                 if(UnipolarWaveforms[wf][binNum]>Amplitude)
                   {
                     Amplitude = UnipolarWaveforms[wf][binNum];
                     PeakBin   = binNum;
                   }
               }
             if(Amplitude>0)
	       {
                 for(int binNum = PeakBin-fPreSamples; binNum!=PeakBin+fMaxHitLength; ++binNum)
                   {
                     if( (binNum>0) && (binNum < int(UnipolarWaveforms[wf].size()) ) )
                       {
                         double ThisVal = UnipolarWaveforms[wf][binNum];
                         sumQ  += ThisVal*binNum;
                         sumQ2 += pow(ThisVal,2)*binNum;
                         Area  += ThisVal;

                         if( (binNum - PeakBin) < fFastCompSamples ) FastComponent += UnipolarWaveforms[wf][binNum];

                       }

                   }


                 HitVar      = pow(sumQ2 - pow(sumQ,2),0.5)/float(Area);
                 PeakTime    = fTimeBegin*1000. + PeakBin *    1000. / fSampleFreq;
                 Width       = HitVar *                        1000. / fSampleFreq;
                 PE          = Area / fPEArea;
                 FastToTotal = FastComponent / Area;

                 if(Area>(fHitThreshold))
                   {
                     int ThisChannel = OpChannelIDs[wf];
                     HitTimes[flash][ThisChannel]    = PeakBin;
		     AllOpHits[ThisChannel][PeakBin] =  recob::OpHit( ThisChannel,
                                                                      PeakTime,
                                                                      Width,
                                                                      Area,
                                                                      Amplitude,
                                                                      PE,
                                                                      FastToTotal);

		     
		   }
	       }	     
	   }
       }
	   
     


    // Now we have to remove duplicates.  
    //  The strategy will be : the first flash making a claim to a hit gets dibbs.
    //  Then any subsequent flashes claiming that hit get merged with the first one.
    
    // Thess maps keeps track of which hits have been claimed by flashes
    std::map<int, std::map<int, bool  > > HitClaimed;
    std::map<int, std::map<int, int   > > HitClaimedByFlash;
    std::map<int, std::vector<int> >      FlashesToCombine;
    
    bool StillCombining=true;
    //   mf::LogInfo("OpNewFlashFinder")<<"Entering loop"<<std::endl;
    
    while(StillCombining)
      { 
	StillCombining=false;

	//	mf::LogInfo("OpNewFlashFinder")<<"Inside while"<<std::endl;
	for(std::map<int, std::map<int, int> >::iterator itFl=HitTimes.begin(); itFl!=HitTimes.end(); ++itFl)
	  {
	     
	    //	     mf::LogInfo("OpNewFlashFinder")<<"Inside for"<<std::endl;
	    // If we don't find anything to combine, we can escape the loop
	   
	    
	    int FlashID = itFl->first;
	    //    mf::LogInfo("OpNewFlashFinder")<<"Starting with " << FlashID<<std::endl;

	    for(std::map<int, int>::iterator itCh = itFl->second.begin();
		itCh!=itFl->second.end(); ++itCh)
	      {
		int OpChannelID = itCh->first;
		int HitTime     = itCh->second;
		
		// Check if this hit has already been claimed

		//		mf::LogInfo("OpNewFlashFinder")<<"Checking hit " <<OpChannelID<<" " <<HitTime<<std::endl;
		if((HitClaimed[OpChannelID][HitTime]==true))
		  {
		    // Lookup who already wants this hit
		    int PreviousClaimant = HitClaimedByFlash[OpChannelID][HitTime];
		    //   mf::LogInfo("OpNewFlashFinder")<<"Is claimed by " << PreviousClaimant<<std::endl;
		    // If it is already this flash, no need to do anything, otherwise
		    //  absorb new flash into old one
		    if(FlashID!=PreviousClaimant)
		      {
			//	 mf::LogInfo("OpNewFlashFinder")<<"Combining " << PreviousClaimant << " with " << FlashID<<std::endl;
			 // Take all of the new flashes hits into the previous flash
			 for(std::map<int, int>::iterator itEat = HitTimes[FlashID].begin();
			     itEat!=HitTimes[FlashID].end();
			     ++itEat)
			   HitTimes[PreviousClaimant][itEat->first] = HitTimes[FlashID][itEat->first];
			 
			 // Replace all references to the second flash in the claims table
			 for(std::map<int,std::map<int, int> >::iterator itClaims=HitClaimedByFlash.begin(); itClaims!=HitClaimedByFlash.end(); ++itClaims)
			   {
			     for(std::map<int, int>::iterator itClaims2 = itClaims->second.begin();
				 itClaims2!=itClaims->second.end(); ++itClaims2)
			       if(itClaims2->second==FlashID)
				 itClaims2->second = PreviousClaimant;
			   }
			 
			 // Erase the new flash ID, get out of the iterator loop
			 // but run through again to do more combining
			 
			 //		 mf::LogInfo("OpNewFlashFinder")<<"Combination made, escaping loop " <<std::endl;
			 HitTimes.erase(FlashID);
			 StillCombining=true;
			 break;
		       }
		       
		   }
		 else
		   {
		     // If this hit is unclaimed, claim it and carry on
		     HitClaimed[OpChannelID][HitTime]=true;
		     HitClaimedByFlash[OpChannelID][HitTime]=FlashID;
		   }
	       }
	     // If we just made a combination, escape the loop
	     //  and start again (since we just screwed up the
	     //  iterator)
	     if(StillCombining==true) break;     
	   }
       }
     

    //  mf::LogInfo("OpNewFlashFinder")<<"Producing final flashes";

     // Alright, we are almost there.  Now we just have to make the flashes.

     // Need a couple of things from the geometry
     art::ServiceHandle<geo::Geometry> geom;
     int NOpChannels = geom->NOpChannels();
     size_t Nplanes = geom->Nplanes();

 
     // These determine the beam window
     double LowerTime = fBeamStartTime * 1000;
     double UpperTime = (fBeamStartTime + fBeamSpillLength)*1000;

     double FastToTotal = 0;

     
     for(std::map<int, std::map<int, int> >::const_iterator itFl=HitTimes.begin(); itFl!=HitTimes.end(); ++itFl)     
       {
	 std::vector<double> PEs(NOpChannels,0);
	 double AveTime=0;
	 double TotalPE=0;

	 // These variables keep track of spatial extent of flash
	 //  in different directions	 
	 double sumy  = 0,  sumz  = 0;
	 double sumy2 = 0,  sumz2 = 0;

	 std::vector<double> sumw(Nplanes,0) ;
	 std::vector<double> sumw2(Nplanes,0);



	 for(std::map<int,int>::const_iterator itHit = itFl->second.begin();
	     itHit!=itFl->second.end(); ++itHit)
	   {
	     int OpChannelID = itHit->first;
	     int HitTime     = itHit->second;
	     
	     double ThisPE = AllOpHits[OpChannelID][HitTime].PE();
	     TotalPE          += ThisPE;
	     PEs[OpChannelID] = ThisPE;
	     
	     AveTime += AllOpHits[OpChannelID][HitTime].PeakTime() * ThisPE;
	     FastToTotal +=AllOpHits[OpChannelID][HitTime].FastToTotal() * ThisPE;
	     
	     // Find out where this hit is in space, and keep the 
	     //  averages adding up
	     unsigned int o=0; unsigned int c=0;
	     geom->OpChannelToCryoOpDet(OpChannelID,o,c);
	     
	     double xyz[3];
	     geom->Cryostat(c).OpDet(o).GetCenter(xyz);

	     for(size_t p=0; p!=Nplanes; p++)
	       {
		 unsigned int w = geom->NearestWire(xyz,p);
		 sumw[p]  += w * ThisPE;
		 sumw2[p] += pow(w,2) * ThisPE;
	       }
	     
	     sumy+=xyz[1] * ThisPE; sumy2+=pow(xyz[1],2) * ThisPE;
	     sumz+=xyz[2] * ThisPE; sumz2+=pow(xyz[2],2) * ThisPE;

	     
	   }

	 // Take averages and SDs where necessary

	 AveTime /= TotalPE;
	 FastToTotal /=TotalPE;
	 
	 double meany = sumy / TotalPE;
	 double meanz = sumz / TotalPE;
	 
	 double widthy = pow(sumy2 -  pow(sumy,2)/TotalPE, 0.5) / pow(TotalPE,0.5);
	 double widthz = pow(sumz2 -  pow(sumz,2)/TotalPE, 0.5) / pow(TotalPE,0.5);

	 std::vector<double> WireCenters(Nplanes,0);
	 std::vector<double> WireWidths(Nplanes,0);

	 for(size_t p=0; p!=Nplanes; ++p)
	   {
	     WireCenters[p] = sumw[p]/TotalPE;
	     WireWidths[p]  = pow(sumw2[p] - pow(sumw[p],2)/TotalPE, 0.5) / pow(TotalPE,0.5);
	   }

	 // Next thing to do is determine if we're on the beam time
	 int IsOnBeamTime = 0;
	 if((AveTime>LowerTime) && (AveTime<UpperTime))
	   IsOnBeamTime = 1;
	 
	 // Finally, store the flash.  Huzzah!
	 FlashPtr->push_back( recob::OpFlash( AveTime,
					      PEs, 
					      IsOnBeamTime,	
					      FastToTotal,
					      meany, 
					      widthy, 
					      meanz, 
					      widthz, 
					      WireCenters, 
					      WireWidths ));

	 
       }


     // Empty all the maps we accumulated
     HitTimes.clear();
     AllOpHits.clear();
     HitClaimed.clear();
     HitClaimedByFlash.clear();
     FlashesToCombine.clear();

     //mf::LogInfo("OpNewFlashFinder")<<"Storing " << FlashPtr->size() << " flash objects into event";

     // And stick both collections into the event
     evt.put(std::move(FlashPtr));
     evt.put(std::move(HitPtr));

     
     return;
  }



} // namespace opdet


