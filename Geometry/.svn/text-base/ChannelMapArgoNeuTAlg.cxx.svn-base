////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapArgoNeuTAlg.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "Geometry/ChannelMapArgoNeuTAlg.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

namespace geo{


  //----------------------------------------------------------------------------
  // Define sort order for cryostats in argoneut configuration
  static bool sortCryoArgoNeuT(const CryostatGeo* c1, const CryostatGeo* c2)
  {
    double xyz1[3] = {0.}, xyz2[3] = {0.};
    double local[3] = {0.}; 
    c1->LocalToWorld(local, xyz1);
    c2->LocalToWorld(local, xyz2);

    return xyz1[0] < xyz2[0];   
  }


  //----------------------------------------------------------------------------
  // Define sort order for tpcs in argoneut configuration.
  static bool sortTPCArgoNeuT(const TPCGeo* t1, const TPCGeo* t2) 
  {
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double local[3] = {0.};
    t1->LocalToWorld(local, xyz1);
    t2->LocalToWorld(local, xyz2);

    // sort TPCs according to x
    if(xyz1[0] < xyz2[0]) return true;

    return false;
  }


  //----------------------------------------------------------------------------
  // Define sort order for planes in argoneut configuration
  static bool sortPlaneArgoNeuT(const PlaneGeo* p1, const PlaneGeo* p2) 
  {
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double local[3] = {0.};
    p1->LocalToWorld(local, xyz1);
    p2->LocalToWorld(local, xyz2);

    // drift direction is negative, plane number increases in drift direction
    return xyz1[0] > xyz2[0];
  }


  //----------------------------------------------------------------------------
  bool sortWireArgoNeuT(WireGeo* w1, WireGeo* w2){
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};

    w1->GetCenter(xyz1); w2->GetCenter(xyz2);

    if( xyz1[2] < xyz2[2] ) return true; 

    return false;
  }


  //----------------------------------------------------------------------------
  ChannelMapArgoNeuTAlg::ChannelMapArgoNeuTAlg()
  {
  }

  //----------------------------------------------------------------------------
  ChannelMapArgoNeuTAlg::~ChannelMapArgoNeuTAlg()
  {
  }

  //----------------------------------------------------------------------------
  void ChannelMapArgoNeuTAlg::Initialize(std::vector<geo::CryostatGeo*> & cgeo)
  {
    fNcryostat = cgeo.size();
    
    mf::LogInfo("ChannelMapArgoNeuTAlg") << "Initializing...";
  

    std::sort(cgeo.begin(), cgeo.end(), sortCryoArgoNeuT);
    for(size_t c = 0; c < cgeo.size(); ++c) 
      cgeo[c]->SortSubVolumes(sortTPCArgoNeuT, sortPlaneArgoNeuT, sortWireArgoNeuT);
    

    fNTPC.resize(fNcryostat);
    fWireCounts.resize(fNcryostat);
    fFirstWireProj.resize(fNcryostat);
    fOrthVectorsY.resize(fNcryostat);
    fOrthVectorsZ.resize(fNcryostat);
    fPlaneBaselines.resize(fNcryostat);
    fWiresPerPlane.resize(fNcryostat);
    fFirstChannelInNextPlane.resize(fNcryostat);
    fFirstChannelInThisPlane.resize(fNcryostat);

    fTopChannel = 0;

    int RunningTotal = 0;

    for(unsigned int cs = 0; cs != fNcryostat; ++cs){
      
      fNTPC[cs] = cgeo[cs]->NTPC();
      
      // Size up all the vectors 
      fWireCounts[cs]             .resize(fNTPC[cs]);
      fFirstWireProj[cs] 	  .resize(fNTPC[cs]);
      fOrthVectorsY[cs]  	  .resize(fNTPC[cs]);
      fOrthVectorsZ[cs]  	  .resize(fNTPC[cs]);
      fPlaneBaselines[cs]	  .resize(fNTPC[cs]);
      fWiresPerPlane[cs] 	  .resize(fNTPC[cs]);
      fFirstChannelInThisPlane[cs].resize(fNTPC[cs]);
      fFirstChannelInNextPlane[cs].resize(fNTPC[cs]);

      for(unsigned int TPCCount = 0; TPCCount != fNTPC[cs]; ++TPCCount){
	unsigned int PlanesThisTPC = cgeo[cs]->TPC(TPCCount).Nplanes();
	fWireCounts[cs][TPCCount]   .resize(PlanesThisTPC);
	fFirstWireProj[cs][TPCCount].resize(PlanesThisTPC);
	fOrthVectorsY[cs][TPCCount] .resize(PlanesThisTPC);
	fOrthVectorsZ[cs][TPCCount] .resize(PlanesThisTPC);
	
	for(unsigned int PlaneCount = 0; PlaneCount != PlanesThisTPC; ++PlaneCount){
	  double ThisWirePitch = cgeo[cs]->TPC(TPCCount).WirePitch(0, 1, PlaneCount);
	  fWireCounts[cs][TPCCount][PlaneCount] = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Nwires();
	  
	  double  WireCentre1[3] = {0.,0.,0.};
	  double  WireCentre2[3] = {0.,0.,0.};
	  
	  double th = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Wire(0).ThetaZ();
	  double sth = sin(th), cth = cos(th);
	  
	  cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Wire(0).GetCenter(WireCentre1,0);
	  cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Wire(1).GetCenter(WireCentre2,0);
	  
	  // figure out if we need to flip the orthogonal vector 
	  // (should point from wire n -> n+1)
	  double OrthY = cth, OrthZ = -sth;
	  if(((WireCentre2[1] - WireCentre1[1])*OrthY 
	      + (WireCentre2[2] - WireCentre1[2])*OrthZ) < 0){
	    OrthZ *= -1;
	    OrthY *= -1;
	  }
	  
	  // Overall we are trying to build an expression that looks like
	  //  int NearestWireNumber = round((worldPos.OrthVector - FirstWire.OrthVector)/WirePitch);     
	  // That runs as fast as humanly possible.
	  //
	  // Casting to an int is much faster than c rounding commands like floor().  We have to add 0.5
	  // to account for rounding up as well as down.  floor(A) ~ (int)(A+0.5).  We account for the
	  // 0.5 in the first wire constant to avoid adding it every time.
	  //
	  // We can also predivide everything by the wire pitch so we don't do this in the loop
	  //
	  // Putting this together into the useful constants we will use later per plane and tpc:
	  fOrthVectorsY[cs][TPCCount][PlaneCount] = OrthY / ThisWirePitch;
	  fOrthVectorsZ[cs][TPCCount][PlaneCount] = OrthZ / ThisWirePitch;
	  
	  fFirstWireProj[cs][TPCCount][PlaneCount]  = WireCentre1[1]*OrthY + WireCentre1[2]*OrthZ;
	  fFirstWireProj[cs][TPCCount][PlaneCount] /= ThisWirePitch;
	  fFirstWireProj[cs][TPCCount][PlaneCount] -= 0.5;
	  
	  // now to count up wires in each plane and get first channel in each plane
	  int WiresThisPlane = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Nwires();	
	  fWiresPerPlane[cs] .at(TPCCount).push_back(WiresThisPlane);
	  fPlaneBaselines[cs].at(TPCCount).push_back(RunningTotal);
	  
	  RunningTotal += WiresThisPlane;

	  fFirstChannelInThisPlane[cs].at(TPCCount).push_back(fTopChannel);
	  fTopChannel += WiresThisPlane;
	  fFirstChannelInNextPlane[cs].at(TPCCount).push_back(fTopChannel);

	}// end loop over planes
      }// end loop over TPCs
    }// end loop over cryostats

    // calculate the total number of channels in the detector
    fNchannels = fTopChannel;

    LOG_DEBUG("ChannelMapArgoNeuT") << "# of channels is " << fNchannels;


    return;
  }
   
  //----------------------------------------------------------------------------
  void ChannelMapArgoNeuTAlg::Uninitialize()
  {
  }

  //----------------------------------------------------------------------------
  std::vector<geo::WireID> ChannelMapArgoNeuTAlg::ChannelToWire(uint32_t channel)  const
  {
    std::vector< geo::WireID > AllSegments; 
    unsigned int cstat = 0;
    unsigned int tpc   = 0;
    unsigned int plane = 0;
    unsigned int wire  = 0;
    
    // first check if this channel ID is legal
    if(channel > fTopChannel)
      throw cet::exception("Geometry") << "ILLEGAL CHANNEL ID for channel " << channel;
    
    // then go find which plane, tpc and cryostat it is in from the information we stored earlier
    for(unsigned int csloop = 0; csloop != fNcryostat; ++csloop){
      for(unsigned int tpcloop = 0; tpcloop != fNTPC[csloop]; ++tpcloop){
	for(unsigned int planeloop = 0; planeloop != fFirstChannelInNextPlane[csloop][tpcloop].size(); ++planeloop){
	  if(channel < fFirstChannelInNextPlane[csloop][tpcloop][planeloop]){
	    cstat = csloop;
	    tpc   = tpcloop;
	    plane = planeloop;
	    wire  = channel - fFirstChannelInThisPlane[cstat][tpcloop][planeloop];
	    break;
	  }	    
	}// end plane loop
      }// end tpc loop
    }// end cryostat loop

    geo::WireID CodeWire(cstat, tpc, plane, wire);
    
    AllSegments.push_back(CodeWire);

    return AllSegments;
  }

  //----------------------------------------------------------------------------
  uint32_t ChannelMapArgoNeuTAlg::Nchannels() const
  {
    return fNchannels;
  }
  
  //----------------------------------------------------------------------------
  WireID ChannelMapArgoNeuTAlg::NearestWireID(const TVector3& worldPos,
					      unsigned int    PlaneNo,
					      unsigned int    TPCNo,
					      unsigned int    cstat) const
  {

    // This part is the actual calculation of the nearest wire number, where we assume
    //  uniform wire pitch and angle within a wireplane
    
    unsigned int NearestWireNumber = (unsigned int)(worldPos[1]*fOrthVectorsY[cstat][TPCNo][PlaneNo] 
						    + worldPos[2]*fOrthVectorsZ[cstat][TPCNo][PlaneNo]      
						    - fFirstWireProj[cstat][TPCNo][PlaneNo]);
    
    unsigned int wireNumber = NearestWireNumber;
    
    // If we are outside of the wireplane range, throw an exception
    // (this response maintains consistency with the previous
    // implementation based on geometry lookup)
    if(NearestWireNumber < 0 ||
       NearestWireNumber >= fWireCounts[cstat][TPCNo][PlaneNo]){
      if(NearestWireNumber < 0 ) wireNumber = 0;
      else                       wireNumber = fWireCounts[cstat][TPCNo][PlaneNo] - 1;
      
      throw cet::exception("Can't Find Nearest Wire") << "for position (" 
						      << worldPos[0] << ","
						      << worldPos[1] << "," 
						      << worldPos[2] << ")"
						      << " approx wire number # " 
						      << wireNumber;
    }
    
    WireID wid(cstat, TPCNo, PlaneNo, wireNumber);
    return wid;

  }
  
  //----------------------------------------------------------------------------
  // This method returns the channel number, assuming the numbering scheme
  // is heirachical - that is, channel numbers run in order, for example:
  //                                             (Ben J Oct 2011)                   
  //                    Wire1     | 0
  //           Plane1 { Wire2     | 1
  //    TPC1 {          Wire3     | 2
  //           Plane2 { Wire1     | 3   increasing channel number
  //                    Wire2     | 4     (with no gaps)
  //    TPC2 { Plane1 { Wire1     | 5
  //           Plane2 { Wire1     | 6
  //                    Wire2     v 7
  //
  uint32_t ChannelMapArgoNeuTAlg::PlaneWireToChannel(unsigned int plane,
						     unsigned int wire,
						     unsigned int tpc,
						     unsigned int cstat) const
  {
    // This is the actual lookup part - first make sure coordinates are legal
    if(tpc   < fNTPC[cstat] &&
       plane < fWiresPerPlane[cstat][tpc].size() &&
       wire  < fWiresPerPlane[cstat][tpc][plane] ){
      // if the channel has legal coordinates, its ID is given by the wire
      // number above the number of wires in lower planes, tpcs and cryostats
      return fPlaneBaselines[cstat][tpc][plane] + wire;
    }
    else{  
      // if the coordinates were bad, throw an exception
      throw cet::exception("ChannelMapArgoNeuTAlg") << "NO CHANNEL FOUND for cryostat,tpc,plane,wire: "
						    << cstat << "," << tpc << "," << plane << "," << wire;
    }
    
    // made it here, that shouldn't happen, return UINT_MAX
    mf::LogWarning("ChannelMapArgoNeuTAlg") << "should not be at the point in the function, returning "
					    << "UINT_MAX";
    return UINT_MAX;

  }


  //----------------------------------------------------------------------------
  const SigType_t ChannelMapArgoNeuTAlg::SignalType(uint32_t const channel) const
  {

    SigType_t sigt;

    if(       channel <  fFirstChannelInThisPlane[0][0][1]     ){ sigt = kInduction;  }
    else if( (channel >= fFirstChannelInThisPlane[0][0][1]) &&
             (channel <  fFirstChannelInNextPlane[0][0][1])    ){ sigt = kCollection; }
    else{    mf::LogWarning("BadChannelSignalType") << "Channel " << channel
						    << " not given signal type." << std::endl;         }
            

    return sigt;
  }


  //----------------------------------------------------------------------------
  const View_t ChannelMapArgoNeuTAlg::View(uint32_t const channel) const
  {

    View_t view;
    
    if(       channel <  fFirstChannelInNextPlane[0][0][0]     ) view = geo::kU; 
    else if( (channel >= fFirstChannelInThisPlane[0][0][1]) &&
	     (channel <  fFirstChannelInNextPlane[0][0][1])    ) view = geo::kV; 
    else if( (channel >= fFirstChannelInThisPlane[0][0][2]) &&
	     (channel <  fFirstChannelInNextPlane[0][0][2])    ) view = geo::kZ; 
    else
      mf::LogWarning("BadChannelSignalType") << "Channel " << channel
					     << " not given view type.";

    return view;
  }  


} // namespace
