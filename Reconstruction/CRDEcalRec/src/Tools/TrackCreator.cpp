#ifndef TRACK_CREATOR_C
#define TRACK_CREATOR_C

#include "Tools/TrackCreator.h"

TrackCreator::TrackCreator(const Settings& m_settings) : settings( m_settings ){

} 

StatusCode TrackCreator::CreateTracks( PandoraPlusDataCol& m_DataCol, std::vector<DataHandle<edm4hep::TrackCollection>*> m_TrkCol ){

  if(settings.m_trackCollections.size()==0 || m_TrkCol.size()==0) StatusCode::SUCCESS;
  m_DataCol.collectionMap_Track.clear(); 

  //Readin from collection
  for(int icol=0; icol<m_TrkCol.size(); icol++){
    std::string col_name = settings.m_trackCollections[icol];
    const edm4hep::TrackCollection* const_TrkCol = m_TrkCol[icol]->get(); 

    for(unsigned int itrk=0; itrk<const_TrkCol->size(); itrk++) m_DataCol.collectionMap_Track[col_name].push_back( const_TrkCol->at(itrk) ) ; 
  }

  //Convert to local object


  return StatusCode::SUCCESS;
}


#endif
