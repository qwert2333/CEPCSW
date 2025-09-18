#ifndef TRACK_CREATOR_H
#define TRACK_CREATOR_H

#include "CyberDataCol.h"
#include "Tools/Algorithm.h"
#include "Algorithm/TrackExtrapolatingAlg.h"
#include "TVector3.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include <fstream>

namespace Cyber{
  class TrackCreator{

  public: 

    //initialize a CaloHitCreator
    TrackCreator( const Settings& m_settings );
    ~TrackCreator() { delete m_TrkExtraAlg; };
   
    StatusCode CreateTracks( CyberDataCol& m_DataCol, 
                             std::vector<DataHandle<edm4hep::TrackCollection>*>& r_TrackCols, 
                             DataHandle<edm4hep::MCRecoTrackParticleAssociationCollection>* r_MCParticleTrkCol ); 
   
    StatusCode CreateTracksFromMCParticle(CyberDataCol& m_DataCol, 
                                          DataHandle<edm4hep::MCParticleCollection>& r_MCParticleCol);

    StatusCode SelectGoodTrack(std::vector<std::shared_ptr<Cyber::Track>>& trkCol);
    StatusCode Reset(){};

  private: 
    const std::map<int, int> PDGIDs = {
      {0, -11},
      {1, -13},
      {2, 211},
      {3, 321},
      {4, 2212},
    };
    const Cyber::Settings  settings; 
    Cyber::Algorithm*      m_TrkExtraAlg; 
    Cyber::Settings        m_TrkExtraSettings;  

    static bool compTrkIP( std::shared_ptr<Cyber::Track> trk1, std::shared_ptr<Cyber::Track> trk2 ){ 
      float IP1 = sqrt(trk1->getD0()*trk1->getD0() + trk1->getZ0()*trk1->getZ0() );
      float IP2 = sqrt(trk2->getD0()*trk2->getD0() + trk2->getZ0()*trk2->getZ0() );
      return IP1<IP2;
    }

    static bool compTrkP( std::shared_ptr<Cyber::Track> trk1, std::shared_ptr<Cyber::Track> trk2 ){
      return trk1->getMomentum() > trk2->getMomentum();
    }

    TMVA::Reader *mva_rdr = new TMVA::Reader("Silent");  

  };
};
#endif
