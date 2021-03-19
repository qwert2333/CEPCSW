#ifndef TRACK_CREATOR_H
#define TRACK_CREATOR_H

#include "PandoraPlusDataCol.h"


class TrackCreator{

public: 

  class Settings{
  public:
    Settings(){};

  };
  
  //initialize a CaloHitCreator
  TrackCreator( Settings& settings ){};
  ~TrackCreator();

  StatusCode ReadSettings( Settings& settings ){ return StatusCode::SUCCESS; };

  StatusCode GetTracks(PandoraPlusDataCol& dataCol ){ return StatusCode::SUCCESS; };



  void Reset(){};

private: 

  



};
#endif
