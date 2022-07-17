#ifndef TRACK_CREATOR_H
#define TRACK_CREATOR_H

#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"
#include "Algorithm/TrackExtrapolatingAlg.h"
#include "TVector3.h"

namespace PandoraPlus{
  class TrackCreator{

  public: 

    //initialize a CaloHitCreator
    TrackCreator( const Settings& m_settings );
    ~TrackCreator() { delete m_TrkExtraAlg; };
   
    StatusCode CreateTracks( PandoraPlusDataCol& m_DataCol, 
                             std::vector<DataHandle<edm4hep::TrackCollection>*>& r_TrackCols ); 
   
   
    StatusCode Reset(){};

  private: 
    const PandoraPlus::Settings  settings; 
    PandoraPlus::Algorithm*      m_TrkExtraAlg; 
    PandoraPlus::Settings        m_TrkExtraSettings;  


  };
};
#endif
