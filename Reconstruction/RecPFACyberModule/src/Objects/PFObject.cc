#ifndef PFOBJECT_C
#define PFOBJECT_C

#include "Objects/PFObject.h"
namespace Cyber{

  void PFObject::Clear()
  {
    m_pid = 0;
    m_tracks.clear();
    m_ecal_clusters.clear();
    m_hcal_clusters.clear();
  }

  std::shared_ptr<Cyber::PFObject> PFObject::Clone() const{
    std::shared_ptr<Cyber::PFObject> m_newpfo = std::make_shared<Cyber::PFObject>();
    m_newpfo->setTrack(m_tracks);
    m_newpfo->setECALCluster(m_ecal_clusters);
    m_newpfo->setHCALCluster(m_hcal_clusters);
    m_newpfo->setPID(m_pid);

    return m_newpfo;
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

  double PFObject::getECALClusterEnergy() const{
    if(m_ecal_clusters.size()==0) return 0.;

    double sumEn = 0;
    for(int i=0; i<m_ecal_clusters.size(); i++) sumEn += m_ecal_clusters[i]->getLongiE();
    return sumEn;

  }

  double PFObject::getHCALClusterEnergy() const{
    if(m_hcal_clusters.size()==0) return 0.;

    double sumEn = 0;
    for(int i=0; i<m_hcal_clusters.size(); i++) sumEn += m_hcal_clusters[i]->getHitsE();
    return sumEn;
  }

  double PFObject::getTrackMomentum() const{
    if(m_tracks.size()==0) return -99;

    double maxP = -1;
    for(int i=0; i<m_tracks.size(); i++){
      if(m_tracks[i]->getMomentum()>maxP) maxP = m_tracks[i]->getMomentum();
    }
    return maxP;
  }


};
#endif
