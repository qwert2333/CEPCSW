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

    //Create a truth track from MC particle. 
    std::vector<CRDEcalEDM::Track> m_trkCol; m_trkCol.clear();
    for(int imc=0; imc<dataCol.MCParticleCol.size(); imc++){
      edm4hep::MCParticle m_MCp = dataCol.MCParticleCol[imc];
      if( fabs(m_MCp.getCharge())!=0 ){ //e+-, mu+-, pi+-, k+-, p/pbar
        TVector3 mc_p(m_MCp.getMomentum()[0], m_MCp.getMomentum()[1], m_MCp.getMomentum()[2]);
        
        ///////
        double momentum_initial = sqrt(m_MCp.getMomentum()[0]*m_MCp.getMomentum()[0]+m_MCp.getMomentum()[1]*m_MCp.getMomentum()[1]+m_MCp.getMomentum()[2]*m_MCp.getMomentum()[2]);
        double theta_initial = mc_p.Theta();
        double phi_initial = mc_p.Phi(); 
            
        CRDEcalEDM::Track m_trk;
        m_trk.setHelix(momentum_initial,theta_initial,phi_initial);
        m_trk.setExtrapolation_points();
        m_trk.setExtrapolation_front_face();
        m_trk.setVertex( m_MCp.getVertex()[0], m_MCp.getVertex()[1], m_MCp.getVertex()[2] );
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
