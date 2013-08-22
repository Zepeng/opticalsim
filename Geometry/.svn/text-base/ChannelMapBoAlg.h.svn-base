////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapBoAlg.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  tjyang@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_CHANNELBOMAPALG_H
#define GEO_CHANNELBOMAPALG_H

#include <vector>
#include <stdint.h>

#include "Geometry/ChannelMapAlg.h"

namespace geo{

  class ChannelMapBoAlg : public ChannelMapAlg{

  public:

    ChannelMapBoAlg();
    ~ChannelMapBoAlg();
    
    void                     Initialize(std::vector<geo::CryostatGeo*> & cgeo);
    void                     Uninitialize();
    std::vector<WireID>      ChannelToWire(uint32_t channel)           const;
    uint32_t                 Nchannels()                               const;
    WireID                   NearestWireID(const TVector3& worldPos,
					   unsigned int    PlaneNo,
					   unsigned int    TPCNo,
					   unsigned int    cstat)      const;
    uint32_t                 PlaneWireToChannel(unsigned int plane,
						unsigned int wire,
						unsigned int tpc,
						unsigned int cstat)    const;
   const View_t              View( uint32_t const channel )            const;
   const SigType_t           SignalType( uint32_t const channel )      const;
    
  private:
    
    unsigned int                                         fNcryostat;      ///< number of cryostats in the detector
    uint32_t                                             fNchannels;      ///< number of channels in the detector
    uint32_t                                             fTopChannel;     ///< book keeping highest channel #
    std::vector<unsigned int>                            fNTPC;           ///< number of TPCs in each cryostat
    std::vector<std::vector<std::vector<float>>>         fFirstWireProj;  ///< Distance (0,0,0) to first wire 	 
                                                                          ///< along orth vector per plane per TPC
    std::vector<std::vector<std::vector<float>>>         fOrthVectorsY;   ///< Unit vectors orthogonal to wires in
    std::vector<std::vector<std::vector<float>>>         fOrthVectorsZ;   ///< each plane - stored as 2 components
                                                                          ///< to avoid having to invoke any bulky
                                                                          ///< TObjects / CLHEP vectors etc	 
    std::vector<std::vector<std::vector<float>>>         fWireCounts;     ///< Number of wires in each plane - for
                                                                          ///< range checking after calculation   
    std::vector<std::vector<std::vector<unsigned int>>>  fPlaneBaselines; ///< The number of wires in all the 
                                                                          ///< tpcs and planes up to this one 
                                                                          ///< in the heirachy
    std::vector<std::vector<std::vector<unsigned int>>>  fWiresPerPlane;  ///< The number of wires in this plane 
                                                                          ///< in the heirachy
  };

}
#endif // GEO_CHANNELMAPBOALG_H

