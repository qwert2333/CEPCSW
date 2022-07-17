#ifndef TRACK_CREATOR_C
#define TRACK_CREATOR_C

#include "Tools/TrackCreator.h"

namespace PandoraPlus{

  TrackCreator::TrackCreator(const Settings& m_settings) : settings( m_settings ){

  };

  StatusCode TrackCreator::CreateTracks( PandoraPlusDataCol& m_DataCol, std::vector<DataHandle<edm4hep::TrackCollection>*>& r_TrackCols ){
    if(r_TrackCols.size()==0 || settings.map_stringVecPars.at("trackCollections").size()==0) StatusCode::SUCCESS;

    //Save readin collections
    m_DataCol.collectionMap_Track.clear(); 
    for(int icol=0; icol<r_TrackCols.size(); icol++){
      const edm4hep::TrackCollection* const_TrkCol = r_TrackCols[icol]->get(); 

      std::vector<edm4hep::Track> m_TrkCol; m_TrkCol.clear();
      for(unsigned int icol=0; icol<const_TrkCol->size(); icol++){
        edm4hep::Track m_trk = const_TrkCol->at(icol);
        m_TrkCol.push_back(m_trk);
      }

      m_DataCol.collectionMap_Track[ settings.map_stringVecPars.at("trackCollections")[icol] ] = m_TrkCol; 
    }


    //Convert to local objects
    std::vector<PandoraPlus::Track*> m_trkCol; m_trkCol.clear();

    for(auto iter : m_DataCol.collectionMap_Track){
    auto const_TrkCol = iter.second; 
    for(int itrk=0; itrk<const_TrkCol.size(); itrk++){
      PandoraPlus::Track* m_trk = new PandoraPlus::Track();

      for(int its=0; its<const_TrkCol[itrk].trackStates_size(); its++){
        PandoraPlus::TrackState m_trkst;
        m_trkst.D0 = const_TrkCol[itrk].getTrackStates(its).D0;
        m_trkst.Z0 = const_TrkCol[itrk].getTrackStates(its).Z0;
        m_trkst.phi0 = const_TrkCol[itrk].getTrackStates(its).phi;
        m_trkst.tanLambda = const_TrkCol[itrk].getTrackStates(its).tanLambda;
        m_trkst.Omega = const_TrkCol[itrk].getTrackStates(its).omega;
        m_trkst.Kappa = m_trkst.Omega*1000./(0.3*settings.map_floatPars.at("BField"));   
        m_trkst.location = const_TrkCol[itrk].getTrackStates(its).location;
        m_trkst.referencePoint.SetXYZ( const_TrkCol[itrk].getTrackStates(its).referencePoint[0],
                                       const_TrkCol[itrk].getTrackStates(its).referencePoint[1],
                                       const_TrkCol[itrk].getTrackStates(its).referencePoint[2] );

        m_trk->AddTrackState( m_trkst );
      }
      m_trkCol.push_back(m_trk);
      m_DataCol.bk_TrackCol.push_back( m_trk );
    }}
    m_DataCol.TrackCol = m_trkCol;

    //Track extrapolation
    //  Write settings: geometry description
    //  m_TrkExtraSettings.map_floatPar["Nlayers"] = 28;

    m_TrkExtraAlg = new TrackExtrapolatingAlg();
    m_TrkExtraAlg->ReadSettings(m_TrkExtraSettings);
    m_TrkExtraAlg->Initialize();
    m_TrkExtraAlg->RunAlgorithm( m_DataCol );
    m_TrkExtraAlg->ClearAlgorithm();

    return StatusCode::SUCCESS;
  }

};
#endif
