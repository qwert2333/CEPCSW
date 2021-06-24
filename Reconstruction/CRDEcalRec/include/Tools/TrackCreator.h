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
  ~TrackCreator() {};

  StatusCode ReadSettings( Settings& settings ){ return StatusCode::SUCCESS; };

  StatusCode GetTracks(PandoraPlusDataCol& dataCol ){ 
    if(dataCol.MCParticleCol.size()==0) return StatusCode::SUCCESS;
    dataCol.ClearTrack();

    std::vector<CRDEcalEDM::Track> m_trkCol; m_trkCol.clear();
    for(int imc=0; imc<dataCol.MCParticleCol.size(); imc++){
      edm4hep::MCParticle m_MCp = dataCol.MCParticleCol[imc];
      if( fabs(m_MCp.getCharge())!=0 ){ //e+-, mu+-, pi+-, k+-, p/pbar
        TVector3 mc_p(m_MCp.getMomentum()[0], m_MCp.getMomentum()[1], m_MCp.getMomentum()[2]);
        double phi0 = mc_p.Phi() + PI/2.; 
        if( phi0>2*PI ) phi0 = phi0-2*PI;
        if( phi0<0 )    phi0 = phi0+2*PI;
        double tanL = tan(PI/2. - mc_p.Theta());
        double d0 = sqrt( m_MCp.getVertex()[0]*m_MCp.getVertex()[0] + m_MCp.getVertex()[1]*m_MCp.getVertex()[1] ); 
        double z0 = m_MCp.getVertex()[2]; 

        CRDEcalEDM::Track m_trk;
        m_trk.setVertex( m_MCp.getVertex()[0], m_MCp.getVertex()[1], m_MCp.getVertex()[2] );
        m_trk.setHelix(d0, z0, phi0, 0., tanL);   //d0, z0, phi0, kappa, tanLambda.
        m_trk.setPID( m_MCp.getPDG() );
        m_trk.setMomentum( mc_p );
        m_trk.setCharge( (int)m_MCp.getCharge() );
        m_trkCol.push_back(m_trk);
      }
    }
    dataCol.TrackCol = m_trkCol; 
    return StatusCode::SUCCESS; 
  };


  StatusCode MatchTrkEcalRelation( PandoraPlusDataCol& dataCol ){
    if(dataCol.TrackCol.size()==0) return StatusCode::SUCCESS;
    for(int itrk=0; itrk<dataCol.TrackCol.size(); itrk++){
      CRDEcalEDM::Track m_trk = dataCol.TrackCol[itrk];
      for(int ib=0; ib<dataCol.BlockVec.size(); ib++){
        if( dataCol.BlockVec[ib].MatchTrk( m_trk ) ) dataCol.BlockVec[ib].addTrk(m_trk);
      }
    }
    return StatusCode::SUCCESS;
  }

  void Reset(){};

private: 

  

};
#endif
