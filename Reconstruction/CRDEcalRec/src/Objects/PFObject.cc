#ifndef PFOBJECT_C
#define PFOBJECT_C

#include "Objects/PFObject.h"
namespace PandoraPlus{

  void PFObject::Clear()
  {
    m_tracks.clear();
    m_ecal_clusters.clear();
    m_hcal_clusters.clear();
  }

  void PFObject::addTrack(const Track* _track){
    if( find( m_tracks.begin(), m_tracks.end(), _track)!=m_tracks.end() ){
      std::cout<<"ERROR: attempt to add an existing track into PFO! Skip it "<<std::endl;
    }
    else{
      m_tracks.push_back(_track);
    }
  }

  void PFObject::addECALCluster(const Calo3DCluster* _ecal_cluster){
    if( find( m_ecal_clusters.begin(), m_ecal_clusters.end(), _ecal_cluster)!=m_ecal_clusters.end() ){
      std::cout<<"ERROR: attempt to add an existing ECAL cluster into PFO! Skip it "<<std::endl;
    }
    else{
      m_ecal_clusters.push_back(_ecal_cluster);
    }
  }

  void PFObject::addHCALCluster(const Calo3DCluster* _hcal_cluster){
    if( find( m_hcal_clusters.begin(), m_hcal_clusters.end(), _hcal_cluster)!=m_hcal_clusters.end() ){
      std::cout<<"ERROR: attempt to add an existing HCAL cluster into PFO! Skip it "<<std::endl;
    }
    else{
      m_hcal_clusters.push_back(_hcal_cluster);
    }
  }

};
#endif