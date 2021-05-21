#ifndef TRACK_CREATOR_H
#define TRACK_CREATOR_H

#include "PandoraPlusDataCol.h"
#include "TVector3.h"

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

  StatusCode GetTracks(PandoraPlusDataCol& dataCol ){ 
    if(dataCol.MCParticleCol.size()==0) return StatusCode::SUCCESS;
    dataCol.ClearTrack();

    std::vector<CRDEcalEDM::Track> m_trkCol; m_trkCol.clear();
    for(int imc=0; imc<dataCol.MCParticleCol.size(); imc++){
      edm4hep::MCParticle m_MCp = dataCol.MCParticleCol[imc];
      if( fabs(m_MCp.getPDG())==211 || fabs(m_MCp.getPDG())==321 || fabs(m_MCp.getPDG())==2212 ){
        TVector3 mc_p(m_MCp.getMomentum()[0], m_MCp.getMomentum()[1], m_MCp.getMomentum()[2]);
        double phi0 = mc_p.Phi() + PI/2.; 
        if( phi0>2*PI ) phi0 = phi0-2*PI;
        if( phi0<0 )    phi0 = phi0+2*PI;
        double tanL = tan(PI/2. - mc_p.Theta());

        CRDEcalEDM::Track m_trk;
        m_trk.setHelix(0., 230., phi0, 0., tanL);   //d0, z0, phi0, kappa, tanLambda. 
        m_trk.setPID( m_MCp.getPDG() );
        m_trkCol.push_back(m_trk);
      }
    }
    dataCol.TrackCol = m_trkCol; 
    return StatusCode::SUCCESS; 
  };


  void Reset(){};

private: 

  

};
#endif
