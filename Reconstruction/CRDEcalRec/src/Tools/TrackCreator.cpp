#ifndef TRACK_CREATOR_C
#define TRACK_CREATOR_C

#include "Tools/TrackCreator.h"

TrackCreator::TrackCreator(const Settings& m_settings) : settings( m_settings ){

} 

StatusCode CreateTracks( PandoraPlusDataCol& m_DataCol ){

  if(settings.m_trackCollections.size()==0) StatusCode::SUCCESS;
  m_DataCol.collectionMap_Track.clear(); 

  for(int icol=0; icol<settings.m_trackCollections.size(); icol++){
    std::string col_name = settings.m_trackCollections[icol];

    DataHandle<edm4hep::TrackCollection>    r_TrkCol{col_name, Gaudi::DataHandle::Reader, this};
    const edm4hep::TrackCollection*         const_TrkCol = r_TrkCol.get(); 

    std::vector<edm4hep::Track> m_TrkCol; m_TrkCol.clear();
    for(unsigned int itrk=0; itrk<const_TrkCol->size(); itrk++){
      edm4hep::Track m_trk = const_TrkCol->at(itrk);
      m_TrkCol.push_back(m_trk);
    }

    m_DataCol.collectionMap_Track[col_name] = m_TrkCol; 
  }

  return StatusCode::SUCCESS;
}


#endif
