////////////////////////////////////////////////////////////////////////
/// \file  LArG4.cxx
/// \brief Use Geant4 to run the LArSoft detector simulation
///
/// \version $Id: LArG4.cxx,v 1.22 2010/07/20 06:08:30 bjpjones Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

/// This a module.  It has the following functions:
///
/// - Initialize Geant4 physics, detector geometry, and other
///   processing.
///
/// - Accept sim::MCTruth objects from the MC branch of the FMWK
///   Event structure.
///
/// - Pass the primary particles to the Geant4 simulation to calculate
///   "truth" information for the detector response.
///
/// - Pass the truth information to the DetSim branch of the FMWK event.

#ifndef LARG4_LARG4_H
#define LARG4_LARG4_H 1

#include "G4Base/G4Helper.h"
#include "G4Base/ConvertMCTruthToG4.h"

#include <cstring>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/Assns.h"
#include "cetlib/exception.h"
#include "cetlib/search_path.h"

// LArSoft Includes
#include "LArG4/LArVoxelReadoutGeometry.h"
#include "LArG4/PhysicsList.h"
#include "LArG4/ParticleListAction.h"
#include "LArG4/G4BadIdeaAction.h"
#include "LArG4/OpDetSensitiveDetector.h"
#include "LArG4/OpDetReadoutGeometry.h"
#include "LArG4/LArStackingAction.h"
#include "LArG4/LArVoxelReadout.h"
#include "LArG4/MaterialPropertyLoader.h"
#include "LArG4/OpDetPhotonTable.h"
#include "LArG4/LBNE10ktOpDetPhotonTable.h"
#include "Simulation/LArG4Parameters.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationBase/MCTruth.h"
#include "Simulation/ParticleList.h"
#include "Simulation/SimPhotons.h"
#include "Simulation/SimChannel.h"
#include "Geometry/Geometry.h"
#include "G4Base/DetectorConstruction.h"
#include "G4Base/UserActionManager.h"

// G4 Includes
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4UImanager.hh"
#include "Geant4/G4VUserDetectorConstruction.hh"
#include "Geant4/G4VUserPrimaryGeneratorAction.hh"
#include "Geant4/G4VUserPhysicsList.hh"
#include "Geant4/G4UserRunAction.hh"
#include "Geant4/G4UserEventAction.hh"
#include "Geant4/G4UserTrackingAction.hh"
#include "Geant4/G4UserSteppingAction.hh"
#include "Geant4/G4UserStackingAction.hh"
#include "Geant4/G4VisExecutive.hh"
#include "Geant4/G4VUserPhysicsList.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4LogicalVolumeStore.hh"
#include "Geant4/Randomize.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/globals.hh"

// ROOT Includes

// C++ Includes
#include <iostream>
#include <cstring>
#include <sys/stat.h>

// Forward declarations
class G4RunManager;
class G4UImanager;
class G4VisExecutive;

///Geant4 interface
namespace larg4 {  
 
  // Forward declarations within namespace.
  class LArVoxelListAction;
  class ParticleListAction;
  
  class LArG4 : public art::EDProducer{
  public:
 
    /// Standard constructor and destructor for an FMWK module.
    explicit LArG4(fhicl::ParameterSet const& pset);
    virtual ~LArG4();

    /// The main routine of this module: Fetch the primary particles
    /// from the event, simulate their evolution in the detctor, and
    /// produce the detector response.
    void produce (art::Event& evt); 
    void beginJob();

  private:

    g4b::G4Helper*       fG4Help;                   ///< G4 interface object					   
    LArVoxelListAction*  flarVoxelListAction;  	    ///< Geant4 user action to accumulate LAr voxel information.
    ParticleListAction*  fparticleListAction;  	    ///< Geant4 user action to particle information.		   

    std::string          fG4PhysListName;           ///< predefined physics list to use if not making a custom one
    std::string          fG4MacroPath;              ///< directory path for Geant4 macro file to be 
                                                    ///< executed before main MC processing.
    bool                 fdumpParticleList;         ///< Whether each event's sim::ParticleList will be displayed.
    bool                 fLBNE10ktOpFast;
    std::string          fOpDetSensitiveVolumeName; ///< Name of volume with which to associate OpDet 
                                                    ///< sensitive detectors	
    int                  fSmartStacking;            ///< Whether to instantiate and use class to 
                                                    ///< dictate how tracks are put on stack.	

  };

} // namespace LArG4

namespace larg4 {

  //----------------------------------------------------------------------
  // Constructor
  LArG4::LArG4(fhicl::ParameterSet const& pset)
    : fG4Help                (0)
    , flarVoxelListAction    (0)
    , fparticleListAction    (0)
    , fG4PhysListName        (pset.get< std::string >("G4PhysListName","larg4::PhysicsList"))
    , fdumpParticleList      (pset.get< bool        >("DumpParticleList")                   )
    , fSmartStacking         (pset.get< int         >("SmartStacking",0)                    )
    , fLBNE10ktOpFast        (pset.get<bool         >("LBNE10ktOpFast",false)                     )
  {
    LOG_DEBUG("LArG4") << "Debug: LArG4()";

    // get the random number seed, use a random default if not specified
    // in the configuration file.
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());
    // setup the random number service for Geant4, the "G4Engine" label is a
    // special tag setting up a global engine for use by Geant4/CLHEP
    createEngine(seed, "G4Engine");


    produces< std::vector<simb::MCParticle>>();
    if(!fLBNE10ktOpFast) produces< std::vector<sim::SimPhotons> >();
    //produces< std::vector<int> >();
    else produces< std::vector<sim::LBNE10ktPhotons> >();
    produces< std::vector<sim::SimChannel> >();
    produces< art::Assns<simb::MCTruth, simb::MCParticle> >();

    // constructor decides if initialized value is a path or an environment variable
    cet::search_path sp("FW_SEARCH_PATH");

    sp.find_file(pset.get< std::string >("GeantCommandFile"), fG4MacroPath);
    struct stat sb;
    if (fG4MacroPath.empty() || stat(fG4MacroPath.c_str(), &sb)!=0)
      // failed to resolve the file name
      throw cet::exception("NoG4Macro") << "G4 macro file "
					<< fG4MacroPath
					<< " not found!";

  }

  //----------------------------------------------------------------------
  // Destructor
  LArG4::~LArG4()
  {
    if(fG4Help) delete fG4Help;
  }

  //----------------------------------------------------------------------
  void LArG4::beginJob()
  {
    art::ServiceHandle<geo::Geometry> geom;

    fG4Help = new g4b::G4Helper(fG4MacroPath, fG4PhysListName);
    fG4Help->ConstructDetector(geom->GDMLFile());

    // Get the logical volume store and assign material properties
    MaterialPropertyLoader* MPL = new MaterialPropertyLoader();
    MPL->GetPropertiesFromServices();
    MPL->UpdateGeometry(G4LogicalVolumeStore::GetInstance());

    // Tell the detector about the parallel LAr voxel geometry.
    std::vector<G4VUserParallelWorld*> pworlds;

    // make a parallel world for each TPC in the detector
    pworlds.push_back( new LArVoxelReadoutGeometry("LArVoxelReadoutGeometry") );
    pworlds.push_back( new OpDetReadoutGeometry( geom->OpDetGeoName() ));

    fG4Help->SetParallelWorlds(pworlds);

    // Intialize G4 physics and primary generator action
    fG4Help->InitPhysics();

    // Use the UserActionManager to handle all the Geant4 user hooks.
    g4b::UserActionManager* uaManager = g4b::UserActionManager::Instance();

    // User-action class for accumulating LAr voxels.
    art::ServiceHandle<sim::LArG4Parameters> lgp;

    // UserAction for getting past a bug in v4.9.4.p02 of Geant4.
    // This action will not be used once the bug has been fixed
    // The techniques used in this UserAction are not to be repeated
    // as in general they are a very bad idea, ie they take a const
    // pointer and jump through hoops to change it
    G4BadIdeaAction *bia = new G4BadIdeaAction(fSmartStacking);
    uaManager->AddAndAdoptAction(bia);

    // User-action class for accumulating particles and trajectories
    // produced in the detector.
    fparticleListAction = new ParticleListAction(lgp->ParticleKineticEnergyCut(),
						 lgp->StoreTrajectories(),
						 lgp->KeepEMShowerDaughters());
    uaManager->AddAndAdoptAction(fparticleListAction);

    // UserActionManager is now configured so continue G4 initialization
    fG4Help->SetUserAction();

    // With an enormous detector with lots of rock ala LAr34 (nee LAr20)
    // we need to be smarter about stacking.
    if (fSmartStacking>0){
      G4UserStackingAction* stacking_action = new LArStackingAction(fSmartStacking);
      fG4Help->GetRunManager()->SetUserAction(stacking_action);
    }

  }

  //----------------------------------------------------------------------
  void LArG4::produce(art::Event& evt)
  {
    LOG_DEBUG("LArG4") << "produce()";

    // loop over the lists and put the particles and voxels into the event as collections
    std::unique_ptr< std::vector<simb::MCParticle> > partCol  (new std::vector<simb::MCParticle  >);
    std::unique_ptr< std::vector<sim::SimChannel>  > scCol    (new std::vector<sim::SimChannel>);
    std::unique_ptr< std::vector<sim::SimPhotons>  > PhotonCol(new std::vector<sim::SimPhotons>);
    //std::unique_ptr< std::vector<int>  > LBNE10ktPhotonCol(new std::vector< int>);
    std::unique_ptr< std::vector<sim::LBNE10ktPhotons>  > LBNE10ktPhotonCol(new std::vector<sim::LBNE10ktPhotons>);
    std::unique_ptr< art::Assns<simb::MCTruth, simb::MCParticle> > tpassn(new art::Assns<simb::MCTruth, simb::MCParticle>);

    // Fetch the lists of LAr voxels and particles.
    art::ServiceHandle<sim::LArG4Parameters> lgp;

    // Clear the detected photon table
    OpDetPhotonTable::Instance()->ClearTable();

    // reset the track ID offset as we have a new collection of interactions
    fparticleListAction->ResetTrackIDOffset();

    //look to see if there is any MCTruth information for this
    //event
    std::vector< art::Handle< std::vector<simb::MCTruth> > > mclists;
    evt.getManyByType(mclists);

    // Need to process Geant4 simulation for each interaction separately.
    for(size_t mcl = 0; mcl < mclists.size(); ++mcl){

      art::Handle< std::vector<simb::MCTruth> > mclistHandle = mclists[mcl];

      for(size_t i = 0; i < mclistHandle->size(); ++i){
	art::Ptr<simb::MCTruth> mct(mclistHandle,i);

	LOG_DEBUG("LArG4") << *(mct.get());

	// The following tells Geant4 to track the particles in this interaction.
	fG4Help->G4Run(mct);

	const sim::ParticleList& particleList = *( fparticleListAction->GetList() );
	for(auto i = particleList.begin(); i != particleList.end(); ++i){
	  // copy the particle so that it isnt const
	  simb::MCParticle p(*(*i).second);
	  partCol->push_back(p);

	  util::CreateAssn(*this, evt, *(partCol.get()), mct, *(tpassn.get()));
	}

	// Has the user request a detailed dump of the output objects?
	if (fdumpParticleList){
	  mf::LogInfo("LArG4") << "Dump sim::ParticleList; size()="
			       << particleList.size() << "\n"
			       << particleList;
	}

      }

    }// end loop over interactions
    
    // get the electrons from the LArVoxelReadout sensitive detector
    // Get the sensitive-detector manager.
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    
    // Find the sensitive detector with the name "LArVoxelSD".
    OpDetSensitiveDetector *theOpDetDet = dynamic_cast<OpDetSensitiveDetector*>(sdManager->FindSensitiveDetector("OpDetSensitiveDetector"));
  
    // Store the contents of the detected photon table
    //
    if(theOpDetDet)
    {
      if(!fLBNE10ktOpFast)
      {      
        LOG_DEBUG("Optical") << "Storing OpDet Hit Collection in Event";
      
        std::map<int, sim::SimPhotons*> ThePhotons = OpDetPhotonTable::Instance()->GetPhotons();
      
        if(ThePhotons.size() > 0){
            for(auto const& it : ThePhotons){
                sim::SimPhotons ph(*(it.second));
                PhotonCol->push_back(ph);
            }
        }
      }
      else
      {
        LOG_DEBUG("Optical") << "Storing OpDet Hit Collection in Event";

        std::map<int, std::map<int, int> > ThePhotons = LBNE10ktOpDetPhotonTable::Instance()->GetPhotons();

        if(ThePhotons.size() > 0)
        {
          for(std::map<int, std::map<int, int>>::const_iterator it = ThePhotons.begin(); it!= ThePhotons.end(); it++)
          {
            sim::LBNE10ktPhotons ph;
            ph.OpChannel = it->first;
            ph.DetectedPhotons = it->second;
            LBNE10ktPhotonCol->push_back(ph);
            //LBNE10ktPhotonCol->push_back(-1);
            //LBNE10ktPhotonCol->push_back(it->first);
            //for(std::map<int, int>::const_iterator in_it = it->second.begin(); in_it!= it->second.end(); in_it++)
            //{
            //  LBNE10ktPhotonCol->push_back(in_it->first);
            //  LBNE10ktPhotonCol->push_back(in_it->second);
            //}
          }
        }
      }
    }
      
    // only put the sim::SimChannels into the event once, not once for every
    // MCTruth in the event

    art::ServiceHandle<geo::Geometry> geom;
    for(unsigned int c = 0; c < geom->Ncryostats(); ++c){

      // map to keep track of which channels we already have SimChannels for in scCol
      // remake this map on each cryostat as channels ought not to be shared between cryostats, just
      // between TPC's
 
      std::map<unsigned int, unsigned int>  channelToscCol;

      unsigned int ntpcs =  geom->Cryostat(c).NTPC();
      for(unsigned int t = 0; t < ntpcs; ++t){
	std::string name("LArVoxelSD");
	char nums[32];
	sprintf(nums, "_Cryostat%u_TPC%d", c, t);
	name += nums;

	G4VSensitiveDetector* sd = sdManager->FindSensitiveDetector(name);
	// If this didn't work, then a sensitive detector with
	// the name "LArVoxelSD" does not exist.
	if ( !sd ){
	  throw cet::exception("LArG4") << "Sensitive detector '"
					<< name 
					<< "' does not exist";
	}

	// Convert the G4VSensitiveDetector* to a LArVoxelReadout*.
	LArVoxelReadout *larVoxelReadout = dynamic_cast<LArVoxelReadout*>(sd);

	// If this didn't work, there is a "LArVoxelSD" detector, but
	// it's not a LArVoxelReadout object.
	if ( !larVoxelReadout ){
	  throw cet::exception("LArG4") << "Sensitive detector '"
					<< name
					<< "' is not a LArVoxelReadout object";
	}

	LOG_DEBUG("LArG4") << "now put the SimChannels in the event";

	larVoxelReadout->MakeSimChannelVector();
	std::vector<sim::SimChannel> channels = larVoxelReadout->GetSimChannels();
	for(auto const& sc : channels){

	  // push sc onto scCol but only if we haven't already put something in scCol for this channel.
	  // if we have, then merge the ionization deposits.  Skip the check if we only have one TPC

	  if (ntpcs > 1) {
	    unsigned int ichan = sc.Channel();
	    std::map<unsigned int, unsigned int>::iterator itertest = channelToscCol.find(ichan);
	    if (itertest == channelToscCol.end()) {
	      channelToscCol[ichan] = scCol->size();
              scCol->push_back(sc);
	    }
	    else {
	      unsigned int idtest = itertest->second;
	      auto const tdcideMap = sc.TDCIDEMap();
   	      for(auto const& tdcide : tdcideMap){
                 for(auto const& ide : tdcide.second){
 	           double xyz[3] = {ide.x, ide.y, ide.z};
 	           scCol->at(idtest).AddIonizationElectrons(ide.trackID,
 				       tdcide.first,
 				       ide.numElectrons,
 				       xyz,
 				       ide.energy);
 	         } // end loop to add ionization electrons to  scCol->at(idtest)
 	      }// end loop over tdc to vector<sim::IDE> map
	    } // end if check to see if we've put SimChannels in for ichan yet or not
	  }
	  else {
	    scCol->push_back(sc);
	  } // end of check if we only have one TPC (skips check for multiple simchannels if we have just one TPC)
	} // end loop over simchannels for this TPC
	larVoxelReadout->ClearSimChannelVector();
      } // end loop over tpcs
    }// end loop over cryostats

    evt.put(std::move(scCol));
    evt.put(std::move(partCol));
    if(!fLBNE10ktOpFast) evt.put(std::move(PhotonCol));
    else evt.put(std::move(LBNE10ktPhotonCol));
    evt.put(std::move(tpassn));

    return;
  }

} // namespace LArG4

namespace larg4 {

  DEFINE_ART_MODULE(LArG4);

} // namespace LArG4

#endif // LARG4_LARG4_H
