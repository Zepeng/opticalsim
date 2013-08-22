// OpFlashFinder, Christie Chiu, MIT 2012
//
// This module associates flashes of light in different PMTs
// into OpFlash objects, which represent flashes in different
// OpDets which occur at similar times, and therefore are
// likely from the same interaction in the detector.
//


// OpFlashFinder.h

#ifndef OpFlashFinder_H
#define OpFlashFinder_H 1

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/OpDetGeo.h"
#include "RawData/OpDetPulse.h"

// ART includes.
#include "art/Framework/Core/EDProducer.h"

// ROOT classes.
class TH1D;  // 1-dimensional histogram class

// C++ includes
#include <cstring>
#include <vector>

namespace opdet {
 
  class OpFlashFinder : public art::EDProducer{
  public:
 
    // Standard constructor and destructor for an ART module.
    OpFlashFinder(const fhicl::ParameterSet&);
    virtual ~OpFlashFinder();

    // This method is called once, at the start of the job. In this
    // example, it will define the histogram we'll write.
    void beginJob();

    // The producer routine, called once per event. 
    void produce (art::Event&); 

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fInputModule;              // Input tag for OpHit collection
    float fSampleFreq;                     // in MHz
    float fTimeBegin;                      // in us
    float fTimeEnd;                        // in us
    float fTimeWidth;

    float fBeamStartTime;                  // in us
    float fBeamSpillLength;                // in us
 
    unsigned int fMultiplicityCondition;
   
  };

} 

#endif // OpFlashFinder_H



// OpFlashFinder_module.cc

// This is a short program required by the ART framework.  It enables
// a program (OpFlashFinder, in this case) to be called as a module
// from a .fcl file. It is unlikely that you'll have to make any
// changes to this file.

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet {
  DEFINE_ART_MODULE(OpFlashFinder);
}


// OpFlashFinder.cxx

#ifndef __CINT__

// LArSoft includes
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
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH1.h"


// C++ Includes
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  OpFlashFinder::OpFlashFinder(fhicl::ParameterSet const& pset)
  {
    // Infrastructure piece
    produces<std::vector< recob::OpFlash> >();

    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get<std::string>("InputModule");

    
    art::ServiceHandle<OpDigiProperties> odp;
    fTimeBegin  = odp->TimeBegin();
    fTimeEnd    = odp->TimeEnd();
    fSampleFreq = odp->SampleFreq();

    fMultiplicityCondition = pset.get<unsigned int>("MultiplicityCondition");

    fBeamStartTime   = pset.get<float>("BeamStartTime");
    fBeamSpillLength = pset.get<float>("BeamSpillLength");
    
    fTimeWidth = pset.get<float>("TimeWidth") * 1000;
  }

  //-----------------------------------------------------------------------
  // Destructor
  OpFlashFinder::~OpFlashFinder() 
  {}
   
  //-----------------------------------------------------------------------
  void OpFlashFinder::beginJob()
  {
  }
   

  //-----------------------------------------------------------------------
  void OpFlashFinder::produce(art::Event& evt) 
  {
    // Infrastructure piece
    std::unique_ptr<std::vector< recob::OpFlash > > StoragePtr (new std::vector<recob::OpFlash>);

    // Get Geometry
    art::ServiceHandle< geo::Geometry > geom;

    // Create a handle for our vector of pulses
    art::Handle< std::vector< recob::OpHit > > PeakHandle;

    // Read in PeakHandle
    evt.getByLabel(fInputModule, PeakHandle);
    std::vector<art::Ptr< recob::OpHit > > PeakPtrs;
    art::fill_ptr_vector(PeakPtrs, PeakHandle);

    // Create vector to store sets of matching hits
    std::vector< std::vector<recob::OpHit> > SubeventSets;




    // First scan: compute PE averages for each OpDet; fill histograms
    std::vector< double > AvgPE(geom->NOpDet(), 0);
    std::vector< int > nAvgPE(geom->NOpDet(), 0);

    for(unsigned int i = 0; i < PeakPtrs.size(); ++i)
      {
	AvgPE.at(PeakPtrs[i]->OpChannel()) += PeakPtrs[i]->PE();
	nAvgPE.at(PeakPtrs[i]->OpChannel())++;

      }
    
    for(unsigned int i = 0; i < AvgPE.size(); ++i)
	AvgPE[i] = AvgPE[i] / nAvgPE[i];


    // Second scan: compute standard deviations for each OpDet
    std::vector< double > StdevPE(geom->NOpDet(), 0);

    for(unsigned int i = 0; i < PeakPtrs.size(); ++i)
      {
	StdevPE.at(PeakPtrs[i]->OpChannel()) += 
	  std::pow( PeakPtrs[i]->PE() - AvgPE.at(PeakPtrs[i]->OpChannel()), 2);
      }

    for(unsigned int i = 0; i < StdevPE.size(); ++i)
	StdevPE[i] = std::pow( StdevPE[i] / nAvgPE[i], 0.5);





    // Third scan: filter & sort peaks into subevents
    // While we still have peaks to sort into subevents
    while (PeakPtrs.size() > 0)
      {
	// Begin by getting OpHit with greatest #PE
	unsigned int iBigPeak = 0;

	// Scan through looking for index of biggest peak
	for(unsigned int i = 0; i < PeakPtrs.size(); ++i)
	  {
	    if( PeakPtrs[i]->PE() > PeakPtrs[iBigPeak]->PE())
	      iBigPeak = i;
	  }

	if ( PeakPtrs.size() == 0)
	  continue;

	// Get OpHit with greatest #PE, and remove from list
	recob::OpHit ThePeak = *PeakPtrs[iBigPeak];
	PeakPtrs.erase( PeakPtrs.begin() + iBigPeak );

	// Create vector for this set of matching hits
	std::vector<recob::OpHit> SubeventPtr;
	SubeventPtr.push_back(ThePeak);

	// Scan through rest of vector looking for other hits that 
	// match in time
	for(unsigned int j = 0; j < PeakPtrs.size(); ++j)
	  {

	    if(abs(ThePeak.PeakTime() - PeakPtrs[j]->PeakTime()) <= fTimeWidth)
	      {
		SubeventPtr.push_back( recob::OpHit(*PeakPtrs[j]) );
		PeakPtrs.erase( PeakPtrs.begin() + j );
	      }
	  }
	
	// Save vector of matching hits in SubeventSets, if it is large enough
	if (SubeventPtr.size() >= fMultiplicityCondition)
	  SubeventSets.push_back(SubeventPtr);
	SubeventPtr.clear();
      }
    





    // Extract #PE and average time for each group of peaks (subevent)
    std::vector< std::vector< double > > SubeventSetsPE;
    std::vector< double > SubeventSetsTotalPE;
    std::vector< double > SubeventSetsTime;

    // Scan through our subevents, saving only the number of PE and avg time
    for( unsigned int i = 0; i < SubeventSets.size(); ++i )
      {
	// Create a vector, nth entry is nth OpDet
	std::vector< double >* SubeventPE = 
	  new std::vector< double>(geom->NOpDet(), 0);

	// Set variable to compute mean arrival time
	double AvgPeakTime = 0;
	double TotalPE = 0;

	for( unsigned int j = 0; j < SubeventSets[i].size(); ++j)
	  {
	    // Fill vector with PE for each OpDet
	    SubeventPE -> at( SubeventSets[i][j].OpChannel() ) = 
	      SubeventSets[i][j].PE();
	  
	    // Add arrival time
	    AvgPeakTime += SubeventSets[i][j].PeakTime() * SubeventSets[i][j].PE();
	    TotalPE+= SubeventSets[i][j].PE();
	  }

	// Fill time vector with arrival time
	SubeventSetsTime.push_back(AvgPeakTime / TotalPE);

	SubeventSetsTotalPE.push_back(TotalPE);

	// Fill SubeventsPE with OpDet PE vector
	SubeventSetsPE.push_back(*SubeventPE);
      }





    // Now scan through subevents and merge if necessary
    int mergeAgain = 1;

    while (mergeAgain)
      {
	// Turn flag off to iterate another time through
	mergeAgain = 0;

	// For subevent A
	for (unsigned int i = 0; i < SubeventSetsPE.size(); ++i)
	  {
	    // For subevent B
	    for (unsigned int j = i+1; j < SubeventSetsPE.size(); ++j)
	      {
		// Times have to match up decently for merging
		if (abs(SubeventSetsTime[i] - SubeventSetsTime[j]) > fTimeWidth)
		  continue;
		
		bool CanMerge=true;

		// For each PMT in the subevents
		for (unsigned int k = 0; k < SubeventSetsPE[i].size(); ++k)
		  {
		    // If a PMT saw events in both, they can't be merged
		    if (SubeventSetsPE[i][k] != 0 && SubeventSetsPE[j][k] != 0)
		      CanMerge=false;
		  }
		
		if(CanMerge)
		  {
		    for (unsigned int l = 0; l < SubeventSetsPE[i].size(); ++l)
		      {
			if (SubeventSetsPE[i][l] == 0)
			  {
			    SubeventSetsPE[i][l] = SubeventSetsPE[j][l];
			    
			    // This part needs to be weighted average
			    SubeventSetsTime[i] = (SubeventSetsTime[i]*SubeventSetsTotalPE[i] + SubeventSetsTime[j]*SubeventSetsTotalPE[j])/(SubeventSetsTotalPE[i]+SubeventSetsTotalPE[j]);			
			    SubeventSetsTotalPE[i]+=SubeventSetsTotalPE[j];
			  }
		      }
		  }
					
		// Delete subevent B
		SubeventSetsPE.erase(SubeventSetsPE.begin() + j);
		SubeventSetsTime.erase(SubeventSetsTime.begin() + j);
		SubeventSetsTotalPE.erase(SubeventSetsTotalPE.begin() + j);
		
		// Turn flag on because we may need to merge again
		mergeAgain = 1;
	      } // For subevent B
	  } // For subevent A
      } // While there are still events to be merged


    double LowerTime = fBeamStartTime * 1000;
    double UpperTime = (fBeamStartTime + fBeamSpillLength)*1000;

    // Print everything to screen so that we can see what's going on
    // Also save subevent information to event
 
    size_t Nplanes = geom->Nplanes();

    for (unsigned int i = 0; i < SubeventSetsPE.size(); ++i)
      {

	double sumy  = 0,  sumz  = 0;
	double sumy2 = 0,  sumz2 = 0;
	
	std::vector<double> sumw(Nplanes,0) ;
	std::vector<double> sumw2(Nplanes,0);
	
	double totalPE = 0;


	for( unsigned int j = 0; j < SubeventSetsPE[i].size(); ++j )
	  {
	    
	    if(SubeventSetsPE[i][j] > 0)
	      {
		//	mf::LogInfo ("Peak Info") << "     OpDet: " << j <<
		//	  "   PE: " << SubeventSetsPE[i][j];

		double PE = SubeventSetsPE[i][j];

		unsigned int o=0; unsigned int c=0;		

		geom->OpChannelToCryoOpDet(j,o,c);
		
		double xyz[3];
		geom->Cryostat(c).OpDet(o).GetCenter(xyz);

		for(size_t p=0; p!=Nplanes; p++)
		  {
		    ///\todo: really need to specify TPC and cryostat as well here
		    unsigned int w = geom->NearestWire(xyz,p);
		    sumw[p]  += w * PE;
		    sumw2[p] += pow(w,2) * PE;
		  }

		sumy+=xyz[1] * PE; sumy2+=pow(xyz[1],2) * PE;
		sumz+=xyz[2] * PE; sumz2+=pow(xyz[2],2) * PE;
		
		totalPE+=PE;
		
	      }

	  }
	
	double meany = sumy / totalPE;
	double meanz = sumz / totalPE;
	
	double widthy = pow(sumy2 -  pow(sumy,2)/totalPE, 0.5) / pow(totalPE,0.5);
	double widthz = pow(sumz2 -  pow(sumz,2)/totalPE, 0.5) / pow(totalPE,0.5);

	std::vector<double> WireCenters(Nplanes,0);
	std::vector<double> WireWidths(Nplanes,0);

	for(size_t p=0; p!=Nplanes; ++p)
	  {
	    WireCenters[p] = sumw[p]/totalPE;
	    WireWidths[p]  = pow(sumw2[p] - pow(sumw[p],2)/totalPE, 0.5) / pow(totalPE,0.5);
	  }


	int IsOnBeamTime = 0;
	if((SubeventSetsTime[i]>LowerTime) && (SubeventSetsTime[i]<UpperTime))
	  IsOnBeamTime = 1;
	
	//	mf::LogInfo("OpFlashFinder")<<"FlashFinderDebug : " << IsOnBeamTime<< " " << LowerTime << " " << SubeventSetsTime[i] << " " << UpperTime<<std::endl;
	
	StoragePtr->push_back( recob::OpFlash( SubeventSetsTime[i],
					       SubeventSetsPE[i], IsOnBeamTime, meany, widthy, meanz, widthz, WireCenters, WireWidths ) );


      }



    evt.put(std::move(StoragePtr));

    return;
  }


} // namespace oprecon

#endif
