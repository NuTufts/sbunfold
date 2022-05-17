#include <iostream>
#include <string>

// ROOT
#include "TFile.h"

// larlite
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/mctruth.h"

int main( int nargs, char** argv )
{
  std::cout << "collect test data" << std::endl;
  std::string inputfile = argv[1];
  std::string outfile   = argv[2];

  larlite::storage_manager ioll( larlite::storage_manager::kREAD );
  ioll.add_in_filename( inputfile );
  ioll.open();

  int nentries = ioll.get_entries();
  std::cout << "number of entries: " << nentries << std::endl;

  TFile* out = new TFile( outfile.c_str(), "NEW" );
  std::vector< float > feature_v(6,0); // (mu px, mu py, mu pz, p px, p py, p pz)
  int mode = 0;
  float Enu = 0;
  float Evis = 0;
  TTree datatree("data","Test data tree");
  datatree.Branch("mode",&mode,"mode/I");
  datatree.Branch("Enu", &Enu, "Enu/F");
  datatree.Branch("Evis", &Evis, "Evis/F");
  datatree.Branch("feature_v", &feature_v);
  
  for (int ientry=0; ientry<nentries; ientry++) {
    ioll.go_to(ientry);

    // clear feature vector
    for (int i=0; i<6; i++) {
      feature_v[i] = 0.;
    }
    
    larlite::event_mctruth* ev_mctruth
      = (larlite::event_mctruth*)ioll.get_data(larlite::data::kMCTruth,"generator");
    if ( ev_mctruth->size()>0 ) {
      auto const& mctruth = ev_mctruth->at(0);
      auto const& mcnu = mctruth.GetNeutrino();
      std::cout << "Neutrino Event: CCNC=" << mcnu.CCNC() << " mode=" << mcnu.Mode() << std::endl;
      mode = mcnu.Mode();
      int nmu = 0;
      int nproton = 0;
      Enu = mcnu.Nu().Trajectory().front().E();
      Evis = 0;
      for ( auto const& part : mctruth.GetParticles() ) {
	int status = part.StatusCode();
	int pdgcode = part.PdgCode();
	auto const& traj = part.Trajectory();
	if ( traj.size()==0 )
	  continue;
	float ke = traj.front().E()-part.Mass();
	std::cout << "[" << status << "] pdgcode=" << pdgcode << " E=" << traj.front().E() << " KE=" << ke << std::endl;
	if ( status==1 ) {
	  // final state particle
	  if ( abs(pdgcode)==13 ) {
	    feature_v[0] = traj.front().Px();
	    feature_v[1] = traj.front().Py();
	    feature_v[2] = traj.front().Pz();
	    nmu++;
	    Evis += ke;
	  }
	  else if (abs(pdgcode)==2212) {
	    feature_v[3] = traj.front().Px();
	    feature_v[4] = traj.front().Py();
	    feature_v[5] = traj.front().Pz();
	    nproton++;
	    Evis += ke;
	  }
	}
      }

      if ( nmu==1 && nproton==1 && mcnu.CCNC()==0 )
	datatree.Fill();
      
    }//end of if has truth
  }//end of entry loop

  ioll.close();
  
  datatree.Write();
  out->Close();
  
  return 0;
}
