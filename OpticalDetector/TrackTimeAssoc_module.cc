//  This module makes a set of hypotheses for track light profiles
//   based on beizer tracks in the event.
//
//  This module will likely be subsumed into the OpFlashFinder module
//   eventually. But at present they are being developped in parallel.
//
//
// \file TrackTimeAssoc.cxx
// \author Ben Jones and Christie Chiu, MIT 2010
//
//



#include "art/Framework/Core/EDProducer.h"
#include "RawData/OpDetPulse.h"

// ROOT includes.
#include <Rtypes.h>
#ifndef TrackTimeAssoc_h
#define TrackTimeAssoc_h 1



namespace trkf{
  class BezierTrack;

}


namespace opdet {

  class TrackTimeAssoc : public art::EDProducer{
  public:
    
    TrackTimeAssoc(const fhicl::ParameterSet&);
    virtual ~TrackTimeAssoc();
    
    void produce(art::Event&);
    void reconfigure(fhicl::ParameterSet const& p);
      
    std::vector<double> GetHypotheses(trkf::BezierTrack* BTrack);
    void PrintHypotheses(std::vector<std::vector<double> > TrackHypotheses);
    double GetLikelihood(std::vector<double> signal, std::vector<double> hypothesis);

    
    void beginJob();
    
    
  private:
    std::string fTrackModuleLabel;
    std::string fFlashModuleLabel;
    int         fBezierResolution;
    int         fdQdxView;
    int         fVerbosity;
  };
}

#endif




////////////////////////////////////////////////////////////////////////
/// \file  TrackTimeAssoc_module.cc
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  bjpjones
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace opdet{

  DEFINE_ART_MODULE(TrackTimeAssoc);

}//end namespace opdet
////////////////////////////////////////////////////////////////////////



//  This module makes a set of hypotheses for track light profiles
//   based on beizer tracks in the event.
//
//  This module will likely be subsumed into the OpFlashFinder module
//   eventually. But at present they are being developped in parallel.
//
//
// \file TrackTimeAssoc.cxx
// \author Ben Jones and Christie Chiu, MIT 2010
//
//

// LArSoft includes
#include "Geometry/Geometry.h"
#include "PhotonPropagation/PhotonVisibilityService.h"
#include "RecoObjects/BezierTrack.h"
#include "RecoBase/OpFlash.h"

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
#include "TH2D.h"

// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

// Debug flag; only used during code development.
const bool debug = true;

namespace opdet {

  //-------------------------------------------------

  TrackTimeAssoc::TrackTimeAssoc(fhicl::ParameterSet const& pset)
  {
    // Infrastructure piece
    //    produces<std::vector< raw::OpDetPulse> >();

    this->reconfigure(pset);
   }


  //-------------------------------------------------

  void TrackTimeAssoc::reconfigure(fhicl::ParameterSet const& pset)
  {
    fTrackModuleLabel = pset.get<std::string>("TrackModuleLabel");
    fBezierResolution = pset.get<int>("BezierResolution");
    fdQdxView         = pset.get<double>("dQdxView");
  }


  //-------------------------------------------------

  void TrackTimeAssoc::beginJob()
  {
  }



  //-------------------------------------------------

  TrackTimeAssoc::~TrackTimeAssoc()
  {
  }

  // Get a hypothesis for the light collected for a bezier track
  std::vector<double> TrackTimeAssoc::GetHypotheses(trkf::BezierTrack* Btrack)
  {
    std::vector<double> ReturnVector;

    art::ServiceHandle<geo::Geometry> geom;
    ReturnVector.resize(geom->NOpDet());
    
    art::ServiceHandle<phot::PhotonVisibilityService> pvs;

    float TrackLength = Btrack->GetLength();

    double xyz[3];
    for (int b=0; b!=fBezierResolution; b++)
      {
	float s               = float(b) / float(fBezierResolution);
	float dQdx            = Btrack->GetdQdx(s, fdQdxView);
	
	Btrack->GetTrackPoint(s,xyz);
	std::vector<float>* PointVisibility = pvs->GetAllVisibilities(xyz);
	float LightAmount = dQdx*TrackLength/float(fBezierResolution);
	
	for(size_t OpDet =0; OpDet!=PointVisibility->size();  OpDet++)
	  {
	    ReturnVector.at(OpDet)+= PointVisibility->at(OpDet) * LightAmount;
	  }
      }
    return ReturnVector;
  }

  //-------------------------------------------------


  void TrackTimeAssoc::produce(art::Event& evt)
  {

    static int EventNo=0;
    


    // Read in flashes from the event
    art::Handle< std::vector<recob::OpFlash> > flashh;
    evt.getByLabel(fFlashModuleLabel, flashh);
    std::vector<art::Ptr<recob::OpFlash> > Flashes;
    for(unsigned int i=0; i < flashh->size(); ++i)
      {
	art::Ptr<recob::OpFlash> flash(flashh,i);
        Flashes.push_back(flash);
      }

    // Read in tracks from the event
    art::Handle< std::vector<recob::Track> > trackh;
    evt.getByLabel(fTrackModuleLabel, trackh);
    std::vector<art::Ptr<recob::Track> >  Tracks;
    for(unsigned int i=0; i < trackh->size(); ++i)
      {
	art::Ptr<recob::Track> track(trackh,i);
        Tracks.push_back(track);
      }


    // Use these to produce Bezier tracks
    std::vector<trkf::BezierTrack*> BTracks;
    BTracks.clear();
    for(size_t i=0; i!=Tracks.size(); i++)
      BTracks.push_back(new trkf::BezierTrack(*Tracks.at(i)));


    std::stringstream ss("");
    ss.flush();
    ss<<"Likelihood matrix"<<EventNo;
    TH2D * LL = new TH2D(ss.str().c_str(),"Chi2 for match; TrackID, FlashID", Tracks.size(), -0.00001, Tracks.size(), Flashes.size(),-0.00001,Flashes.size());
      
    ss.flush();
    ss<<"OnBeamChi2Ratio"<<EventNo;

    TH1D * OnBeamChi2Ratio = new TH1D(ss.str().c_str(), "On Beam Chi2/Chi2 min; TrackID; Relative Chi2", Tracks.size(),-0.0001, Tracks.size()); 
 
    ss.flush();
    ss<<"OnBeamChi2"<<EventNo;

    TH1D * OnBeamChi2 = new TH1D(ss.str().c_str(), "On Beam Chi2; TrackID; Relative Chi2", Tracks.size(), -0.0001, Tracks.size()); 
    


    art::ServiceHandle<geo::Geometry> geom;
    size_t NOpDets = geom->NOpDet();
    
    std::map<int, bool> OnBeamFlashes;
    
    std::vector<std::vector<double> > TrackHypotheses;
    std::vector<std::vector<double> > FlashShapes;
       

    // For each track
    for (size_t i=0; i!=BTracks.size(); ++i)
      {
	TrackHypotheses.push_back(GetHypotheses(BTracks.at(i)));
      }
    
    for(size_t f=0; f!=Flashes.size(); ++f)
      {
	std::vector<double> ThisFlashShape(NOpDets);
	for(size_t i=0; i!=NOpDets; ++i)
	  ThisFlashShape[i]=Flashes.at(f)->PE(i);
	FlashShapes.push_back(ThisFlashShape);
	if(Flashes.at(f)->OnBeamTime()) OnBeamFlashes[f]=true;
      }
	

    art::ServiceHandle<art::TFileService> tfs;
    
    
    std::map<int, std::map<int, double> > Chi2Map;
    for(size_t i=0; i!=TrackHypotheses.size(); ++i)
      {
	double BestChi2ThisTrack       = 0;
	double BestBeamChi2ThisTrack   = 0;
	for(size_t j=0; j!=FlashShapes.size(); ++j)	    
	  {
	    
	    double Chi2 = GetLikelihood(FlashShapes.at(i),TrackHypotheses.at(j));
	    Chi2Map[i][j]=Chi2;
	    
	    LL->Fill(i,j,Chi2);
	    if(Chi2<BestChi2ThisTrack) 
	      BestChi2ThisTrack = Chi2; 
	    if(OnBeamFlashes[j]) 
	      if(Chi2<BestBeamChi2ThisTrack)
		BestBeamChi2ThisTrack=Chi2;
	  }
	OnBeamChi2->Fill(i,BestBeamChi2ThisTrack);
	OnBeamChi2Ratio->Fill(i, BestBeamChi2ThisTrack / BestChi2ThisTrack);
      }
    
 
    
    EventNo++;
  
  }


//--------------------------------------------------
  
  double TrackTimeAssoc::GetLikelihood(std::vector<double> signal, std::vector<double> hypothesis)
  {
       
    double SignalIntegral;
    double HypoIntegral;

    for(size_t i=0; i!=signal.size(); ++i)
      {
	SignalIntegral+=signal.at(i);
	HypoIntegral+=hypothesis.at(i);
      }

    double NormFactor = SignalIntegral/HypoIntegral;
    
    double Chi2=0;
    
    for(size_t i=0; i!=signal.size(); ++i)
      {    
	if(hypothesis.at(i)!=0)
	  Chi2 += (NormFactor * hypothesis.at(i) - signal.at(i))/pow(NormFactor*hypothesis.at(i),0.5);	
      }
    return Chi2;
  }
  

  //--------------------------------------------------

  void TrackTimeAssoc::PrintHypotheses(std::vector<std::vector<double> > TrackHypotheses)
  {
    // List the light hypotheses per track, per PMT
    for (size_t i=0; i!=TrackHypotheses.size(); ++i)
      {
	mf::LogVerbatim("TrackTimeAssoc")<< "Visbility for track " << i <<std::endl;

	for(size_t j=0; j!=TrackHypotheses.at(i).size(); ++j)
	  {
	    mf::LogVerbatim("TrackTimeAssoc") << "Signal at PMT " << j << ", "  << TrackHypotheses.at(i).at(j)<<std::endl;
	    
	  }
      }
    
  }

}


