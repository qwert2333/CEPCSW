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

  StatusCode GetTracks(PandoraPlusDataCol& dataCol, const edm4hep::TrackCollection& const_TrkCol){
    std::vector<CRDEcalEDM::Track> m_trkCol; m_trkCol.clear(); 
    for(int itrk=0; itrk<const_TrkCol.size(); itrk++){
      CRDEcalEDM::Track m_trk;
      m_trk.setHelix( const_TrkCol[itrk].getTrackStates(0).phi, 
                      const_TrkCol[itrk].getTrackStates(0).D0, 
                      const_TrkCol[itrk].getTrackStates(0).Z0, 
                      const_TrkCol[itrk].getTrackStates(0).omega *1000 / (0.3*3.) ,  //NOTE: hard-coding 3T B-field here!  
                      const_TrkCol[itrk].getTrackStates(0).tanLambda ); 
      m_trk.setExtrapolation_points();
      m_trk.setExtrapolation_front_face();
      m_trk.setVertex( const_TrkCol[itrk].getTrackStates(0).referencePoint[0], const_TrkCol[itrk].getTrackStates(0).referencePoint[1], const_TrkCol[itrk].getTrackStates(0).referencePoint[2] );
      m_trk.setCharge( const_TrkCol[itrk].getTrackStates(0).omega>0 ? 1 : -1 );
      //m_trk.setMomentum( mc_p ); //TODO: Need to set the momentum by hand, or write a function to calculate. 
      m_trkCol.push_back(m_trk);
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
