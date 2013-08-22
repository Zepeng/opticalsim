////////////////////////////////////////////////////////////////////////
/// \file OpDetPhotonTable.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Implementation of the OpDetPhotonTable class.
//
// See comments in the OpDetPhotonTable.h file.
//
// Ben Jones, MIT, 11/12/12
//


#include "LBNE10ktOpDetPhotonTable.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/OpDetGeo.h"

namespace larg4 {
  LBNE10ktOpDetPhotonTable * TheLBNE10ktOpDetPhotonTable;
  
  //--------------------------------------------------
  LBNE10ktOpDetPhotonTable::LBNE10ktOpDetPhotonTable()
  {
    fDetectedPhotons.clear();
  }

  //--------------------------------------------------
  LBNE10ktOpDetPhotonTable * LBNE10ktOpDetPhotonTable::Instance()
  {
    if(!TheLBNE10ktOpDetPhotonTable){
      TheLBNE10ktOpDetPhotonTable = new LBNE10ktOpDetPhotonTable;
    }
    return TheLBNE10ktOpDetPhotonTable;  
  }


  //--------------------------------------------------
  void LBNE10ktOpDetPhotonTable::AddPhoton(std::map<int, std::map<int, int>>* StepPhotonTable)
  {
    //for(std::map<int,std::map<int, int>>::iterator it=StepPhotonTable->begin(); it!=StepPhotonTable->end(); ++it)
    for(auto it = StepPhotonTable->begin(); it!=StepPhotonTable->end(); it++)
    {
      //int OpChannel = it->first;
      for(auto in_it = it->second.begin(); in_it!=it->second.end(); in_it++)
      {
        fDetectedPhotons[it->first][in_it->first]+= in_it->second;
      }
    }
    //if(!fDetectedPhotons[opchannel]) 
    //  {
    //    for(size_t i = 0; i < photon->Time.size(); i++)
    //    {
    //      fDetectedPhotons[opchannel].push_back(photon->Time.at(i));
    //    }
	//fDetectedPhotons[opchannel]->SetChannel(opchannel);
    //  }
    //fDetectedPhotons[opchannel]->push_back(*photon);
    //    mf::LogInfo("OpDetPhotonTable") << "Registering detection of a photon in opchannel " <<opchannel<<std::endl;
    
  }


  //--------------------------------------------------
  void LBNE10ktOpDetPhotonTable::ClearTable()
  {
    for(std::map<int,std::map<int, int>>::iterator it=fDetectedPhotons.begin(); it!=fDetectedPhotons.end(); ++it)
      delete &(it->second);
    fDetectedPhotons.clear();
  }

  //--------------------------------------------------
  std::map<int, std::map<int, int>> LBNE10ktOpDetPhotonTable::GetPhotons()
  {
    return fDetectedPhotons;
  }

  //--------------------------------------------------
  std::map<int,int>& LBNE10ktOpDetPhotonTable::GetPhotonsForOpChannel(int opchannel)
  {
    return fDetectedPhotons[opchannel];
  }
  

}
