#ifndef VERTEX_CREATOR_H
#define VERTEX_CREATOR_H

#include "PandoraPlusDataCol.h"


class VertexCreator{

public: 

  class Settings{
  public:
    Settings(){};

  };
  
  //initialize a CaloHitCreator
  VertexCreator( Settings& settings ){};
  ~VertexCreator();

  StatusCode ReadSettings( Settings& settings ){ return StatusCode::SUCCESS; };

  StatusCode GetVertex(PandoraPlusDataCol& dataCol ){ return StatusCode::SUCCESS; };



  void Reset(){};

private: 

  



};
#endif
