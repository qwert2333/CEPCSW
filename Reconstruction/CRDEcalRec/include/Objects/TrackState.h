#ifndef _TRACKSTATE_H
#define _TRACKSTATE_H
#include "TVector3.h"

namespace PandoraPlus{

  class TrackState{
  public: 
    TrackState() {}
    ~TrackState() {};
    void Clear() {};

    int location; 
    float D0;
    float Z0; 
    float phi0; 
    float Kappa;  
    float tanLambda; 
    float Omega;   //Kappa = omega*1000/(0.3*B[T])
    TVector3 referencePoint; 

    static const int AtOther = 0 ; // any location other than the ones defined below
    static const int AtIP = 1 ;
    static const int AtFirstHit = 2 ;
    static const int AtLastHit = 3 ;
    static const int AtCalorimeter = 4 ;
    static const int AtVertex = 5 ;
    static const int LastLocation = AtVertex  ;
    

  };
};
#endif
