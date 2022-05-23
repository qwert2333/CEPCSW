#ifndef TRACK_CREATOR_H
#define TRACK_CREATOR_H

#include "PandoraPlusDataCol.h"
#include "TVector3.h"

class TrackCreator{

public: 

  class Settings{
  public:
    Settings(){};

    std::vector<std::string>  m_trackCollections; 

  };
  
  //initialize a CaloHitCreator
  TrackCreator( const Settings& m_settings );
  ~TrackCreator() {};

  StatusCode CreateTracks( PandoraPlusDataCol& m_DataCol, std::vector<DataHandle<edm4hep::TrackCollection>*> m_TrkCol); 

  StatusCode Reset() { return StatusCode::SUCCESS; };

private: 

  const Settings settings; 
  

};
#endif
