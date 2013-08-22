////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadout.cxx
/// \brief A Geant4 sensitive detector that accumulates voxel information.
///
/// \version $Id: LArVoxelReadout.cxx,v 1.4 2009/08/12 20:52:14 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

#include "LArG4/LArVoxelReadout.h"
#include "LArG4/ParticleListAction.h"
#include "Utilities/DetectorProperties.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geant4/G4HCofThisEvent.hh"
#include "Geant4/G4TouchableHistory.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/Randomize.hh"
#include "Geant4/G4Poisson.hh"

#include <iostream>
#include <ctime>
#include <vector>

#include <stdio.h>

namespace larg4 {

  // Constructor.  Note that we force the name of this sensitive
  // detector to be the value expected by LArVoxelListAction.
  LArVoxelReadout::LArVoxelReadout(std::string const& name)
    : G4VSensitiveDetector(name), flniG4 (getLArNest())
      // grab/instantiate the LArNest singleton.
  {
    // Initialize values for the electron-cluster calculation.
    fChannelMap.clear();
    fChannels.clear();

    art::ServiceHandle<util::DetectorProperties> detprop;
    fSampleRate       = detprop->SamplingRate();
    fTriggerOffset    = detprop->TriggerOffset();

    sscanf(name.c_str(),"%*19s%1u%*4s%u",&fCstat,&fTPC);
  }

  //---------------------------------------------------------------------------------------
  LArVoxelReadout::~LArVoxelReadout() {}

  //---------------------------------------------------------------------------------------
  // Called at the start of each event.
  void LArVoxelReadout::Initialize(G4HCofThisEvent*)
  {
    fElectronLifetime      = fLarpHandle->ElectronLifetime();
    for (int i = 0; i<3; ++i)
      fDriftVelocity[i]    = fLarpHandle->DriftVelocity(fLarpHandle->Efield(i),
							fLarpHandle->Temperature())/1000.;
    double density         = fLarpHandle->Density(fLarpHandle->Temperature());
    fEfield                = fLarpHandle->Efield();
    fElectronClusterSize   = fLgpHandle->ElectronClusterSize();
    fLongitudinalDiffusion = fLgpHandle->LongitudinalDiffusion();
    fTransverseDiffusion   = fLgpHandle->TransverseDiffusion();
    fGeVToElectrons        = fLgpHandle->GeVToElectrons();
    fRecombA               = fLgpHandle->RecombA();
    // the recombination coefficient is in g/(MeVcm^2), but
    // we report energy depositions in MeV/cm, need to divide
    // Recombk from the LArG4Parameters service by the density
    // of the argon we got above.
    fRecombk               = fLgpHandle->Recombk()/density;
    fModBoxA               = fLgpHandle->ModBoxA();
    fModBoxB               = fLgpHandle->ModBoxB()/density;
    fDontDriftThem         = fLgpHandle->DisableWireplanes();
    fUseModBoxRecomb       = fLgpHandle->UseModBoxRecomb();

    fArgon39DecayRate      = fLarpHandle->Argon39DecayRate();
    fNESTOn                = fLgpHandle->NESTOn();
    fnStepsNest            = 0;
    fnStepsTot             = 0;


    mf::LogInfo("LArVoxelReadout")  << " e lifetime: "        << fElectronLifetime
				    << "\n Efield: "          << fLarpHandle->Efield()
				    << "\n Temperature: "     << fLarpHandle->Temperature()
				    << "\n Drift velocity: "  << fDriftVelocity[0]
				    <<" "<<fDriftVelocity[1]<<" "<<fDriftVelocity[2]
				    << "\n Drift electrons: " << fDontDriftThem
                                    << "\n Argon 39 Decay Rate: " << fArgon39DecayRate
                                    << "\n Use Modified Box Recombination: "<< fUseModBoxRecomb;
  }
  //---------------------------------------------------------------------------------------
  // Called at the end of each event.
  void LArVoxelReadout::EndOfEvent(G4HCofThisEvent*)
  {

    // report fraction of good NEST steps, if enabled
      if (fNESTOn) 
	mf::LogInfo("LArVoxelReadou::EndOfEvent")  << " Fraction of good NEST steps: "        << (double) fnStepsNest/fnStepsTot;
    // put in Argon 39 radioactive decays.  
    /// \todo -- make optical flashes to go along with these decays

    // get readout window from DetectorProperties service
    // get Argon 39 rate from LArG4 parameter
    // call Drift Ionilzation Electrons
    // negative track number (-39? -3939?)

    // break each radioactive decay into two steps to allow hits on two wires.
    // distribute the hits throughout the drift volume 

    // consider calling util::LarProperties::Eloss and ElossVar
    // to get more precise range and fluctuate it randomly.  Probably doesn't matter much
 
    if (fArgon39DecayRate > 0)
      {
	const geo::TPCGeo &tpcg = fGeoHandle->TPC(fTPC,fCstat);


	art::ServiceHandle<util::DetectorProperties> detprop;
	double samplerate  = detprop->SamplingRate(); // in ns/tick
	int readwindow = detprop->ReadOutWindowSize(); // in clock ticks
	double timewindow = samplerate*readwindow;
	double width = tpcg.HalfWidth()*2.0;
	double height = tpcg.HalfHeight()*2.0;
	double length = tpcg.Length();
	double tpcvolume = length*width*height;

	double e39 = 0.565;  // Beta energy in MeV
	double dx = 0.565/2.12; // range assuming MIP (use Eloss to get a better range)
	double d2 = dx/2.0;
	double expected_decays = tpcvolume*fArgon39DecayRate*timewindow*1E-9;  // the 1E-9 is to convert the ns to seconds

	// TPC global position

	double xyz1[3] = {0.,0.,0.};
	double xyz2[3] = {0.,0.,0.};
	tpcg.LocalToWorld(xyz1,xyz2);  // gives a point in the middle
	xyz2[0] -= width/2.0;
	xyz2[1] -= height/2.0;
	xyz2[2] -= length/2.0;

	int ndecays = G4Poisson(expected_decays);
	for (int i=0;i<ndecays;++i)
	  {
	    // randomize the simulation time and position

	    double xyz3[3]={0.,0.,0.};

	    bool scc = false;
	    xyz3[0] = G4UniformRand();
	    xyz3[1] = G4UniformRand();
	    xyz3[2] = G4UniformRand();

	    xyz3[0] = (xyz3[0]*width  + xyz2[0])*cm;
	    xyz3[1] = (xyz3[1]*height + xyz2[1])*cm;
	    xyz3[2] = (xyz3[2]*length + xyz2[2])*cm;

	    G4ThreeVector midpoint1(xyz3[0],xyz3[1],xyz3[2]);
	    G4ThreeVector midpoint2;

	    scc = false;
	    for (int itr=0;itr<20;itr++)  // fixed number of tries to make sure we stay in the TPC for the second endpoint.
	      {
	        G4ThreeVector decpath(1.0,0.0,0.0);  
	        decpath.setMag(d2*cm);  // displacement between midpoints is half the total length
	        double angle = G4UniformRand()*TMath::Pi()*2.0;
	        decpath.setPhi(angle);
		double ct=2.0*(G4UniformRand()-0.5);
		if (ct<-1.0) ct=-1.0;
		if (ct>1.0) ct=1.0;
	        angle = TMath::ACos(ct);
	        decpath.setTheta(angle);
	        midpoint2 = midpoint1 + decpath;
		if (  (midpoint2.x() > xyz2[0]*cm && midpoint2.x() < (xyz2[0]+width)*cm) &&
		      (midpoint2.y() > xyz2[1]*cm && midpoint2.y() < (xyz2[1]+height)*cm) &&
		      (midpoint2.z() > xyz2[2]*cm && midpoint2.z() < (xyz2[2]+length)*cm) )
		  {
		    scc = true;
		    break;
		  }
	      }

	    if (scc)
	      {
		DriftIonizationElectrons(e39/2.0,d2,midpoint1,0.0,-(UINT_MAX - 10),firstrad,readwindow);
		DriftIonizationElectrons(e39/2.0,d2,midpoint2,0.0,-(UINT_MAX - 10),subsequentrad,readwindow);
	      }
	  }
      }
  }

  //---------------------------------------------------------------------------------------
  void LArVoxelReadout::clear()
  {
  }

  //---------------------------------------------------------------------------------------
  void LArVoxelReadout::MakeSimChannelVector()
  {
    // iterate over the map and fill the vector of channels
    std::map<uint32_t, sim::SimChannel>::const_iterator itr = fChannelMap.begin();
    while( itr != fChannelMap.end() ){

      fChannels.push_back((itr->second));
      itr++;
    }

    return;
  }

  //--------------------------------------------------------------------------------------
  void LArVoxelReadout::ClearSimChannelVector()
  {
    fChannelMap.clear();
    fChannels.clear();
  }

  //---------------------------------------------------------------------------------------
  // Called for each step.
  G4bool LArVoxelReadout::ProcessHits( G4Step* step, G4TouchableHistory* )
  {
    // All work done for the "parallel world" "box of voxels" in
    // LArVoxelReadoutGeometry makes this a fairly simple routine.
    // First, the usual check for non-zero energy:

    // Only process the hit if the step is inside the active volume and
    // it has deposited energy.  The hit being inside the active volume
    // is virtually sure to happen because the LArVoxelReadoutGeometry
    // that this class makes use of only has voxels for inside the TPC.

    // The step can be no bigger than the size of the voxel,
    // because of the geometry set up in LArVoxelGeometry and the
    // transportation set up in PhysicsList.  Find the mid-point
    // of the step.


    //-----------------------------------------------------------------------------------
    G4double energy = step->GetTotalEnergyDeposit();
    if ( energy > 0 && !fDontDriftThem ){

      G4ThreeVector midPoint = 0.5*( step->GetPreStepPoint()->GetPosition()
				     + step->GetPostStepPoint()->GetPosition() );
      double g4time = step->GetPreStepPoint()->GetGlobalTime();

      // Find the Geant4 track ID for the particle responsible for depositing the
      // energy.  if we are only storing primary EM shower particles, and this energy
      // is from a secondary etc EM shower particle, the ID returned is the primary
      const int trackID = ParticleListAction::GetCurrentTrackID();

      // Now calculate the number of ionization electrons produced in this step
      // and drift them to the readout wires
      G4ThreeVector totstep = step->GetPostStepPoint()->GetPosition();
      totstep -= step->GetPreStepPoint()->GetPosition();

      // Note that if there is no particle ID for this energy deposit, the
      // trackID will be sim::NoParticleId.
      fnStepsTot++;
      if (!fNESTOn) 
	DriftIonizationElectrons(energy/MeV, totstep.mag()/cm, midPoint, g4time, trackID);
      else if (fNESTOn && flniG4.NumberIonizationElectrons()>0 && 
	       (!flniG4.G4Stepper().compare("NEST"))
	       )
	{
	  fnStepsNest++;
	  // This is predicated upon NEST being called first, and _then_ this
	  // stepper for each step. 
	  // Same call here as !fNESTOn, but someday it may not be, ...
	  DriftIonizationElectrons(energy/MeV, totstep.mag()/cm, midPoint, g4time, trackID);

	  /*
	  	  std::cout << "LArVoxelReadout: NumIonElectrons calculated:  " << flniG4.NumberIonizationElectrons() << std::endl;
	  	  std::cout << "LArVoxelReadout: NumScintPhotons calculated:  " << flniG4.NumberScintillationPhotons() << std::endl;
	  	  std::cout << "\t in " << flniG4.G4Stepper() << " at following 4-Location: " << std::endl;
	  	  std::cout << "LArVoxelReadout: " << flniG4.Last4Location()[0] <<  "," << flniG4.Last4Location()[1] <<   "," << flniG4.Last4Location()[2] << "," << flniG4.Last4Location()[3]  << "." << std::endl;
	  	  std::cout << "LArVoxelReadout: " << " energy deposited (MeV), step length (cm) " << energy/MeV <<  "," << totstep.mag()/cm  << "." << std::endl;
	  */
	}
      
    }

    std::string tmp("LArVoxelReadout");
    flniG4.SetG4Stepper(tmp);
	    
    //    flniG4.Reset();
    return true;
  }

  //----------------------------------------------------------------------------
  // energy is passed in with units of MeV, dx has units of cm
  void LArVoxelReadout::DriftIonizationElectrons(double energy,
						 double dx,
						 G4ThreeVector stepMidPoint,
						 const double simTime,
						 int trackID,
						 Radio_t radiological,
						 unsigned int tickmax
						 )
  {

    // Get the recombination factor for this voxel - Nucl.Instrum.Meth.A523:275-286,2004
    // R = A/(1 + (dE/dx)*k)
    // dE/dx is given by the voxel energy deposition, but have to convert it to MeV/cm
    // from GeV/voxel width
    // A = 0.800 +/- 0.003
    // k = (0.097+/-0.001) g/(MeVcm^2)
    // the dx depends on the trajectory of the step
    // k should be divided by the density as our dE/dx is in MeV/cm,
    // the division is handled in the constructor when we set fRecombk
    // B.Baller: Add Modified Box recombination - ArgoNeuT result submitted to JINST

    // This routine gets called frequently, once per every particle
    // traveling through every voxel. Use whatever tricks we can to
    // increase its execution speed.

    static double LifetimeCorr_const = -1000. * fElectronLifetime;
    static double nElectrons_const   = fGeVToElectrons  *1.e-3;
    static double LDiff_const        = std::sqrt(2.*fLongitudinalDiffusion);
    static double TDiff_const        = std::sqrt(2.*fTransverseDiffusion);
    static double RecipDriftVel[3]   = {1./fDriftVelocity[0], 
					1./fDriftVelocity[1], 
					1./fDriftVelocity[2]};

    static int radiologicaltdcoffset=0;  // static so remembers for subsequent calls for radiological sim

    // Map of electrons to store - catalogued by map[channel][tdc]
    static std::map<uint32_t, std::map<unsigned int,double> >  ElectronsToStore;
    static std::map<uint32_t, std::map<unsigned int,double> >  EnergyToStore;

    double recomb;
    double dEdx = energy/dx;
    // Guard against spurious values of dE/dx. Note: assumes density of LAr
    if(dEdx < 1.) dEdx = 1.;
    if(fUseModBoxRecomb) {
      if (dx){
	double Xi = fModBoxB * dEdx / fEfield;
	recomb = log(fModBoxA + Xi) / Xi;
      }
      else recomb = 0;
    } else {
      recomb = fRecombA/(1. + dEdx * fRecombk / fEfield);
    }

    static double xyz1[3] = {0.};

    double xyz[3] = {stepMidPoint.x() / cm,
		     stepMidPoint.y() / cm,
		     stepMidPoint.z() / cm};

    // Already know which TPC we're in because each LArVoxelReadout instance is
    // associated with just one

    try{
      const geo::TPCGeo &tpcg = fGeoHandle->TPC(fTPC,fCstat);

      // X drift distance - the drift direction can be either in
      // the positive or negative direction, so use std::abs
      double XDrift = std::abs(stepMidPoint.x()/cm - tpcg.PlaneLocation(0)[0]);

      if(XDrift < 0.) return;

      // Drift time (nano-sec)
      double TDrift             = XDrift * RecipDriftVel[0];
      if (tpcg.Nplanes() == 2){// special case for ArgoNeuT (plane 0 is the second wire plane)
	TDrift = ((XDrift - tpcg.PlanePitch(0,1)) * RecipDriftVel[0] 
		  + tpcg.PlanePitch(0,1) * RecipDriftVel[1]);
      }
	  
      double lifetimecorrection = TMath::Exp(TDrift / LifetimeCorr_const);
      double nElectrons         = lifetimecorrection * energy * recomb * nElectrons_const;
      if (fNESTOn) 
	{
	  // Just grab the number of electrons as produced by NEST and
	  // attenuate by lifetime, not recomb, as that's already been done in NEST.
	  nElectrons = flniG4.NumberIonizationElectrons() * lifetimecorrection;
	  energy = flniG4.CurrentStepEnergyDeposit();
	}

      // Longitudinal & transverse diffusion sigma (cm)
      double SqrtT    = std::sqrt(TDrift);
      double LDiffSig = SqrtT * LDiff_const;
      double TDiffSig = SqrtT * TDiff_const;
      int nClus       = 1 + (int)(nElectrons / fElectronClusterSize);

      // Compute arrays of values as quickly as possible.
      std::vector< double > XDiff(nClus);
      std::vector< double > YDiff(nClus);
      std::vector< double > ZDiff(nClus);
      std::vector< double > nElDiff(nClus, fElectronClusterSize);
      std::vector< double > nEnDiff(nClus, fElectronClusterSize);

      // Drift the last cluster with smaller size
      nElDiff[nClus-1] = nElectrons - (nClus-1)*fElectronClusterSize;

      for(size_t xx = 0; xx < nElDiff.size(); ++xx){
	if(nElectrons > 0) nEnDiff[xx] = energy/nElectrons*nElDiff[xx];
	else               nEnDiff[xx] = 0.;
      }

      // Smear drift times by x position and drift time
      G4RandGauss::shootArray( nClus, &XDiff[0], 0., LDiffSig);

      // Smear the Y,Z position by the transverse diffusion
      G4RandGauss::shootArray( nClus, &YDiff[0], stepMidPoint.y()/cm,TDiffSig);
      G4RandGauss::shootArray( nClus, &ZDiff[0], stepMidPoint.z()/cm,TDiffSig);

      // make a collection of electrons for each plane
      for(size_t p = 0; p < tpcg.Nplanes(); ++p){

	double Plane0Pitch = tpcg.Plane0Pitch(p);
	
	// "-" sign is because Plane0Pitch output is positive. Andrzej
	xyz1[0] = tpcg.PlaneLocation(0)[0] - Plane0Pitch;
	
	// Drift nClus electron clusters to the induction plane
	for(int k = 0; k < nClus; ++k){
	  // Correct drift time for longitudinal diffusion and plane
	  double TDiff = TDrift + XDiff[k] * RecipDriftVel[0];
	  // Take into account different Efields between planes
	  // Also take into account special case for ArgoNeuT where Nplanes = 2.
	  for (size_t ip = 0; ip<p; ++ip){
	    TDiff += tpcg.PlanePitch(ip,ip+1) * RecipDriftVel[tpcg.Nplanes()==3?ip+1:ip+2];
	  }
	  xyz1[1] = YDiff[k];
	  xyz1[2] = ZDiff[k];
	  
	  /// \todo think about effects of drift between planes
	  
	  // grab the nearest channel to the xyz position
	  try{
	    uint32_t channel = fGeoHandle->NearestChannel(xyz1, p, fTPC, fCstat);
	    
	    /// \todo check on what happens if we allow the tdc value to be
	    /// \todo beyond the end of the expected number of ticks
	    // Add potential decay/capture/etc delay effect, simTime.
	    unsigned int tdc = (unsigned int)((TDiff + simTime) / fSampleRate);
	    tdc += fTriggerOffset;

	    // on the first time we decide a tdc for a radiological, throw a random time for it.
	    // save it for other clusters, planes, and second (or subsequent) calls to this function

            if (radiological != notradiological)
	      {
	        if (p == 0 && k == 0)
		  {
                    double tmptdc = ((double) tickmax) * G4UniformRand();
	            if (radiological == firstrad) radiologicaltdcoffset = tmptdc - tdc;
		  }
		if ( (int) tdc >= -radiologicaltdcoffset)  {
		    tdc += radiologicaltdcoffset;
		  }
		else {
		  continue;  // skip the radiological deposition if it results in a negative tdc.
		  // should only happen due to diffusion or wire plane spacing.
		}
	      }
	    
	    // Add electrons produced by each cluster to the map
	    EnergyToStore[channel][tdc]    += nEnDiff[k];
	    ElectronsToStore[channel][tdc] += nElDiff[k];
	  }
	  catch(cet::exception &e){
	    mf::LogWarning("LArVoxelReadout") << "unable to drift electrons from point ("
					      << xyz[0] << "," << xyz[1] << "," << xyz[2]
					      << ") with exception " << e;
	  }
	} // end loop over clusters
      } // end loop over planes
      
      // Now store them in SimChannels
      
      // check if the current channel is already in the map, otherwise add it
      for(auto const& itread : ElectronsToStore){
	
	uint32_t channel = itread.first;
	std::map<uint32_t, sim::SimChannel>::iterator itchannelmap = fChannelMap.find(channel);
	
	if( itchannelmap != fChannelMap.end() ){
	  for(auto const& itreadinner : itread.second)
	    itchannelmap->second.AddIonizationElectrons(trackID,
							itreadinner.first,
							itreadinner.second,
							xyz,
							EnergyToStore[channel][itreadinner.first]);
	}
	else{
	  sim::SimChannel sc(channel);
	  for(auto const& itreadinner : itread.second)
	    sc.AddIonizationElectrons(trackID,
				      itreadinner.first,
				      itreadinner.second,
				      xyz,
				      EnergyToStore[channel][itreadinner.first]);
	  
	  fChannelMap[channel] = sc;
	}
      }
      ElectronsToStore.clear();
      EnergyToStore.clear();
      
    } // end try intended to catch points where TPC can't be found
    catch(cet::exception &e){
      mf::LogWarning("LArVoxelReadout") << "step cannot be found in a TPC\n"
					<< e;
    }
    
    return;
  }

  //---------------------------------------------------------------------------------------
  // Never used but still have to be defined for G4
  void LArVoxelReadout::DrawAll()  {}
  void LArVoxelReadout::PrintAll() {}

} // namespace larg4
